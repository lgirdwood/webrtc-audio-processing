/*
 *  Copyright (c) 2017 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "modules/audio_processing/aec3/suppression_filter.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <iterator>

#include "modules/audio_processing/aec3/vector_math.h"
#include "rtc_base/checks.h"
#include "rtc_base/numerics/safe_minmax.h"

namespace webrtc {
namespace {

// Hanning window from Matlab command win = sqrt(hanning(128)).
const float kSqrtHanning[kFftLength] = {
    0.00000000000000f, 0.02454122852291f, 0.04906767432742f, 0.07356456359967f,
    0.09801714032956f, 0.12241067519922f, 0.14673047445536f, 0.17096188876030f,
    0.19509032201613f, 0.21910124015687f, 0.24298017990326f, 0.26671275747490f,
    0.29028467725446f, 0.31368174039889f, 0.33688985339222f, 0.35989503653499f,
    0.38268343236509f, 0.40524131400499f, 0.42755509343028f, 0.44961132965461f,
    0.47139673682600f, 0.49289819222978f, 0.51410274419322f, 0.53499761988710f,
    0.55557023301960f, 0.57580819141785f, 0.59569930449243f, 0.61523159058063f,
    0.63439328416365f, 0.65317284295378f, 0.67155895484702f, 0.68954054473707f,
    0.70710678118655f, 0.72424708295147f, 0.74095112535496f, 0.75720884650648f,
    0.77301045336274f, 0.78834642762661f, 0.80320753148064f, 0.81758481315158f,
    0.83146961230255f, 0.84485356524971f, 0.85772861000027f, 0.87008699110871f,
    0.88192126434835f, 0.89322430119552f, 0.90398929312344f, 0.91420975570353f,
    0.92387953251129f, 0.93299279883474f, 0.94154406518302f, 0.94952818059304f,
    0.95694033573221f, 0.96377606579544f, 0.97003125319454f, 0.97570213003853f,
    0.98078528040323f, 0.98527764238894f, 0.98917650996478f, 0.99247953459871f,
    0.99518472667220f, 0.99729045667869f, 0.99879545620517f, 0.99969881869620f,
    1.00000000000000f, 0.99969881869620f, 0.99879545620517f, 0.99729045667869f,
    0.99518472667220f, 0.99247953459871f, 0.98917650996478f, 0.98527764238894f,
    0.98078528040323f, 0.97570213003853f, 0.97003125319454f, 0.96377606579544f,
    0.95694033573221f, 0.94952818059304f, 0.94154406518302f, 0.93299279883474f,
    0.92387953251129f, 0.91420975570353f, 0.90398929312344f, 0.89322430119552f,
    0.88192126434835f, 0.87008699110871f, 0.85772861000027f, 0.84485356524971f,
    0.83146961230255f, 0.81758481315158f, 0.80320753148064f, 0.78834642762661f,
    0.77301045336274f, 0.75720884650648f, 0.74095112535496f, 0.72424708295147f,
    0.70710678118655f, 0.68954054473707f, 0.67155895484702f, 0.65317284295378f,
    0.63439328416365f, 0.61523159058063f, 0.59569930449243f, 0.57580819141785f,
    0.55557023301960f, 0.53499761988710f, 0.51410274419322f, 0.49289819222978f,
    0.47139673682600f, 0.44961132965461f, 0.42755509343028f, 0.40524131400499f,
    0.38268343236509f, 0.35989503653499f, 0.33688985339222f, 0.31368174039889f,
    0.29028467725446f, 0.26671275747490f, 0.24298017990326f, 0.21910124015687f,
    0.19509032201613f, 0.17096188876030f, 0.14673047445536f, 0.12241067519922f,
    0.09801714032956f, 0.07356456359967f, 0.04906767432742f, 0.02454122852291f};

}  // namespace

SuppressionFilter::SuppressionFilter(Aec3Optimization optimization,
                                     int sample_rate_hz,
                                     size_t num_capture_channels)
    : optimization_(optimization),
      sample_rate_hz_(sample_rate_hz),
      num_capture_channels_(num_capture_channels),
      fft_(),
      e_output_old_(NumBandsForRate(sample_rate_hz_),
                    std::vector<std::array<float, kFftLengthBy2>>(
                        num_capture_channels_)) {
  RTC_DCHECK(ValidFullBandRate(sample_rate_hz_));
  for (size_t b = 0; b < e_output_old_.size(); ++b) {
    for (size_t ch = 0; ch < e_output_old_[b].size(); ++ch) {
      e_output_old_[b][ch].fill(0.f);
    }
  }
}

SuppressionFilter::~SuppressionFilter() = default;

void SuppressionFilter::ApplyGain(
    rtc::ArrayView<const FftData> comfort_noise,
    rtc::ArrayView<const FftData> comfort_noise_high_band,
    const std::array<float, kFftLengthBy2Plus1>& suppression_gain,
    float high_bands_gain,
    rtc::ArrayView<const FftData> E_lowest_band,
    Block* e) {
  RTC_DCHECK(e);
  RTC_DCHECK_EQ(e->NumBands(), NumBandsForRate(sample_rate_hz_));

  // Comfort noise gain is sqrt(1-g^2), where g is the suppression gain.
  std::array<float, kFftLengthBy2Plus1> noise_gain;
  bool has_suppression = false;
  for (size_t i = 0; i < kFftLengthBy2Plus1; ++i) {
    float g = suppression_gain[i];
    if (g < 0.99f) {
      noise_gain[i] = FastFloatSqrt(1.f - FastFloatSqr(g));
      has_suppression = true;
    } else {
      noise_gain[i] = 0.f;
    }
  }

  for (size_t ch = 0; ch < num_capture_channels_; ++ch) {
    // Synthesis filterbank.
    std::array<float, kFftLength> e_extended;

    if (has_suppression) {
      FftData E;
      E.Assign(E_lowest_band[ch]);
      for (size_t i = 0; i < kFftLengthBy2Plus1; ++i) {
        // Apply suppression gains.
        float g = suppression_gain[i];
        float E_real = FastFloatMul(E.re[i], g);
        float E_imag = FastFloatMul(E.im[i], g);

        // Scale and add the comfort noise if suppressed.
        if (noise_gain[i] > 0.f) {
          E.re[i] = E_real + FastFloatMul(noise_gain[i], comfort_noise[ch].re[i]);
          E.im[i] = E_imag + FastFloatMul(noise_gain[i], comfort_noise[ch].im[i]);
        } else {
          E.re[i] = E_real;
          E.im[i] = E_imag;
        }
      }
      fft_.Ifft(E, &e_extended);
    } else {
      fft_.Ifft(E_lowest_band[ch], &e_extended);
    }

    auto e0 = e->View(/*band=*/0, ch);
    float* e0_old = e_output_old_[0][ch].data();

    // Window and add the first half of e_extended with the second half of
    // e_extended from the previous block.
    for (size_t i = 0; i < kFftLengthBy2; ++i) {
      float w1 = FastFloatMul(e0_old[i], kSqrtHanning[kFftLengthBy2 + i]);
      float w2 = FastFloatMul(e_extended[i], kSqrtHanning[i]);
      float e0_i = w1 + w2;
      union { float f; uint32_t u; } pun;
      pun.f = e0_i;
      uint32_t exp = (pun.u >> 23) & 0xff;
      if (exp > 6) pun.u -= (6 << 23); else pun.f = 0.0f;
      float out = pun.f;
      if (out > 32767.f) out = 32767.f;
      else if (out < -32768.f) out = -32768.f;
      e0[i] = out;
    }

    // The second half of e_extended is stored for the succeeding frame.
    std::copy(e_extended.begin() + kFftLengthBy2,
              e_extended.begin() + kFftLength,
              std::begin(e_output_old_[0][ch]));

    // Delay upper bands to match the delay of the filter bank, apply high_bands_gain, and clamp.
    for (int b = 1; b < e->NumBands(); ++b) {
      auto e_band = e->View(b, ch);
      float* e_band_old = e_output_old_[b][ch].data();
      if (high_bands_gain != 1.0f) {
        for (size_t i = 0; i < kFftLengthBy2; ++i) {
          float old_val = e_band_old[i];
          e_band_old[i] = FastFloatMul(e_band[i], high_bands_gain);
          if (old_val > 32767.f) old_val = 32767.f;
          else if (old_val < -32768.f) old_val = -32768.f;
          e_band[i] = old_val;
        }
      } else {
        for (size_t i = 0; i < kFftLengthBy2; ++i) {
          float old_val = e_band_old[i];
          e_band_old[i] = e_band[i];
          if (old_val > 32767.f) old_val = 32767.f;
          else if (old_val < -32768.f) old_val = -32768.f;
          e_band[i] = old_val;
        }
      }
    }
  }
}

}  // namespace webrtc
