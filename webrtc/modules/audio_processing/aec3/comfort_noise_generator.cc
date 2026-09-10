/*
 *  Copyright (c) 2017 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "modules/audio_processing/aec3/comfort_noise_generator.h"

// Defines WEBRTC_ARCH_X86_FAMILY, used below.
#include "rtc_base/system/arch.h"

#if defined(WEBRTC_ARCH_X86_FAMILY)
#include <emmintrin.h>
#endif
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numeric>

#include "common_audio/signal_processing/include/signal_processing_library.h"
#include "modules/audio_processing/aec3/vector_math.h"
#include "rtc_base/checks.h"

namespace webrtc {

namespace {

// Computes the noise floor value that matches a WGN input of noise_floor_dbfs.
float GetNoiseFloorFactor(float noise_floor_dbfs) {
  // kdBfsNormalization = 20.f*log10(32768.f).
  constexpr float kdBfsNormalization = 90.30899869919436f;
  return 64.f * powf(10.f, (kdBfsNormalization + noise_floor_dbfs) * 0.1f);
}

// Table of sqrt(2) * sin(2*pi*i/32).
constexpr float kSqrt2Sin[32] = {
    +0.0000000f, +0.2758994f, +0.5411961f, +0.7856950f, +1.0000000f,
    +1.1758756f, +1.3065630f, +1.3870398f, +1.4142136f, +1.3870398f,
    +1.3065630f, +1.1758756f, +1.0000000f, +0.7856950f, +0.5411961f,
    +0.2758994f, +0.0000000f, -0.2758994f, -0.5411961f, -0.7856950f,
    -1.0000000f, -1.1758756f, -1.3065630f, -1.3870398f, -1.4142136f,
    -1.3870398f, -1.3065630f, -1.1758756f, -1.0000000f, -0.7856950f,
    -0.5411961f, -0.2758994f};

}  // namespace

ComfortNoiseGenerator::ComfortNoiseGenerator(const EchoCanceller3Config& config,
                                             Aec3Optimization optimization,
                                             size_t num_capture_channels)
    : optimization_(optimization),
      seed_(42),
      num_capture_channels_(num_capture_channels),
      noise_floor_(GetNoiseFloorFactor(config.comfort_noise.noise_floor_dbfs)),
      N2_initial_(
          std::make_unique<std::vector<std::array<float, kFftLengthBy2Plus1>>>(
              num_capture_channels_)),
      Y2_smoothed_(num_capture_channels_),
      N2_(num_capture_channels_) {
  for (size_t ch = 0; ch < num_capture_channels_; ++ch) {
    (*N2_initial_)[ch].fill(0.f);
    Y2_smoothed_[ch].fill(0.f);
    N2_[ch].fill(1.0e6f);
  }
}

ComfortNoiseGenerator::~ComfortNoiseGenerator() = default;

void ComfortNoiseGenerator::Compute(
    bool saturated_capture,
    rtc::ArrayView<const std::array<float, kFftLengthBy2Plus1>>
        capture_spectrum,
    rtc::ArrayView<FftData> lower_band_noise,
    rtc::ArrayView<FftData> upper_band_noise) {
  const auto& Y2 = capture_spectrum;

  if (!saturated_capture) {
    bool update_n2 = (N2_counter_ > 50);
    bool has_initial = (N2_initial_ != nullptr);
    if (has_initial && ++N2_counter_ == 1000) {
      N2_initial_.reset();
      has_initial = false;
    }

    for (size_t ch = 0; ch < num_capture_channels_; ++ch) {
      float* y2_sm = Y2_smoothed_[ch].data();
      const float* y2 = Y2[ch].data();
      float* n2 = N2_[ch].data();
      float* n2_init = has_initial ? (*N2_initial_)[ch].data() : nullptr;

      constexpr float kSmFactor = 0.90018f;
      constexpr float kN2Factor = 0.10002f;
      for (size_t k = 0; k < kFftLengthBy2Plus1; ++k) {
        // 1. Smooth Y2: 0.9f * sm + 0.1f * y_val
        float y_val = y2[k];
        float sm = y2_sm[k];
        sm = FastFloatAddPos(FastFloatMul(0.9f, sm), FastFloatMul(0.1f, y_val));
        y2_sm[k] = sm;

        // 2. Update N2 from Y2_smoothed
        float cur_n2 = n2[k];
        if (update_n2) {
          cur_n2 = (sm < cur_n2) ? FastFloatAddPos(FastFloatMul(kSmFactor, sm), FastFloatMul(kN2Factor, cur_n2))
                                 : FastFloatMul(cur_n2, 1.0002f);
        }
        if (cur_n2 < noise_floor_) {
          cur_n2 = noise_floor_;
        }
        n2[k] = cur_n2;

        // 3. Update N2_initial from N2 if active
        if (n2_init) {
          float cur_init = n2_init[k];
          if (cur_n2 > cur_init) {
            cur_init = FastFloatAddPos(FastFloatMul(0.999f, cur_init), FastFloatMul(0.001f, cur_n2));
          } else {
            cur_init = cur_n2;
          }
          if (cur_init < noise_floor_) {
            cur_init = noise_floor_;
          }
          n2_init[k] = cur_init;
        }
      }
    }
  }
}

void ComfortNoiseGenerator::GenerateComfortNoise(
    rtc::ArrayView<FftData> lower_band_noise,
    rtc::ArrayView<FftData> upper_band_noise) {
  const auto& N2 = N2_initial_ ? (*N2_initial_) : N2_;

  for (size_t ch = 0; ch < num_capture_channels_; ++ch) {
    FftData* N_low = &lower_band_noise[ch];
    const float* n2 = N2[ch].data();

    // Compute square root spectrum using FastFloatSqrt.
    std::array<float, kFftLengthBy2Plus1> N;
    for (size_t k = 0; k < kFftLengthBy2Plus1; ++k) {
      N[k] = FastFloatSqrt(n2[k]);
    }

#if defined(WEBRTC_ARCH_X86_FAMILY) || defined(WEBRTC_HAS_NEON)
    FftData* N_high = &upper_band_noise[ch];
    constexpr float kOneByNumBands = 1.f / (kFftLengthBy2Plus1 / 2 + 1);
    constexpr int kFftLengthBy2Plus1By2 = kFftLengthBy2Plus1 / 2;
    const float high_band_noise_level =
        std::accumulate(N.begin() + kFftLengthBy2Plus1By2, N.end(), 0.f) *
        kOneByNumBands;
    N_high->re[0] = N_high->re[kFftLengthBy2] = 0.f;
#endif

    N_low->re[0] = N_low->re[kFftLengthBy2] = 0.f;
    for (size_t k = 1; k < kFftLengthBy2; k++) {
      constexpr int kIndexMask = 32 - 1;
      seed_ = (seed_ * 69069 + 1) & (0x80000000 - 1);
      int i = seed_ >> 26;

      const float x = kSqrt2Sin[i];
      const float y = kSqrt2Sin[(i + 8) & kIndexMask];

      N_low->re[k] = FastFloatMul(N[k], x);
      N_low->im[k] = FastFloatMul(N[k], y);

#if defined(WEBRTC_ARCH_X86_FAMILY) || defined(WEBRTC_HAS_NEON)
      N_high->re[k] = high_band_noise_level * x;
      N_high->im[k] = high_band_noise_level * y;
#endif
    }
  }
}

}  // namespace webrtc
