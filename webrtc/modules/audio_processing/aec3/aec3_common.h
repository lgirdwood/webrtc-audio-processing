/*
 *  Copyright (c) 2016 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef MODULES_AUDIO_PROCESSING_AEC3_AEC3_COMMON_H_
#define MODULES_AUDIO_PROCESSING_AEC3_AEC3_COMMON_H_

#include <stddef.h>
#include <stdint.h>

namespace webrtc {

static inline float FastFloatSqr(float x) {
  union { float f; uint32_t u; } pun;
  pun.f = x;
  uint32_t u = pun.u;
  int exp = (u >> 23) & 0xff;
  if (exp == 0) return 0.0f;
  uint64_t m = (1ULL << 23) | (u & 0x7fffff);
  uint64_t m_sq = (m * m) >> 23;
  int new_exp = 2 * exp - 127;
  if (m_sq & (1ULL << 24)) {
    m_sq >>= 1;
    new_exp += 1;
  }
  if (new_exp <= 0) return 0.0f;
  if (new_exp >= 255) new_exp = 254;
  pun.u = ((uint32_t)new_exp << 23) | (uint32_t)(m_sq & 0x7fffff);
  return pun.f;
}

static inline float FastFloatAddPos(float a, float b) {
  union { float f; uint32_t u; } pa, pb, pr;
  pa.f = a; pb.f = b;
  if (pa.u < pb.u) { uint32_t tmp = pa.u; pa.u = pb.u; pb.u = tmp; }
  if (pb.u == 0) return pa.f;
  int exp_a = (pa.u >> 23) & 0xff;
  int exp_b = (pb.u >> 23) & 0xff;
  int diff = exp_a - exp_b;
  if (diff > 24) return pa.f;
  uint32_t ma = (1U << 23) | (pa.u & 0x7fffff);
  uint32_t mb = (1U << 23) | (pb.u & 0x7fffff);
  uint32_t m_sum = ma + (mb >> diff);
  int new_exp = exp_a;
  if (m_sum & (1U << 24)) {
    m_sum >>= 1;
    new_exp++;
  }
  if (new_exp >= 255) new_exp = 254;
  pr.u = ((uint32_t)new_exp << 23) | (m_sum & 0x7fffff);
  return pr.f;
}

static inline float FastMagSqr(float re, float im) {
  return FastFloatAddPos(FastFloatSqr(re), FastFloatSqr(im));
}

static inline float FastFloatMul(float a, float b) {
  union { float f; uint32_t u; } pa, pb, pr;
  pa.f = a; pb.f = b;
  uint32_t sign = (pa.u ^ pb.u) & 0x80000000;
  int exp_a = (pa.u >> 23) & 0xff;
  int exp_b = (pb.u >> 23) & 0xff;
  if (exp_a == 0 || exp_b == 0) return 0.0f;
  int new_exp = exp_a + exp_b - 127;
  uint64_t ma = (1ULL << 23) | (pa.u & 0x7fffff);
  uint64_t mb = (1ULL << 23) | (pb.u & 0x7fffff);
  uint64_t prod = (ma * mb) >> 23;
  if (prod & (1ULL << 24)) {
    prod >>= 1;
    new_exp++;
  }
  if (new_exp <= 0) return 0.0f;
  if (new_exp >= 255) new_exp = 254;
  pr.u = sign | ((uint32_t)new_exp << 23) | (uint32_t)(prod & 0x7fffff);
  return pr.f;
}

static inline float FastFloatInv(float x) {
  union { float f; uint32_t u; } conv, px;
  conv.f = x;
  uint32_t sign = conv.u & 0x80000000;
  conv.u &= 0x7fffffff;
  if (conv.u == 0) return 0.0f;
  conv.u = 0x7ef311c3 - conv.u;
  float y = conv.f;
  px.f = x; px.u &= 0x7fffffff;
  float ax = px.f;
  float axy = FastFloatMul(ax, y);
  y = FastFloatMul(y, 2.0f - axy);
  axy = FastFloatMul(ax, y);
  y = FastFloatMul(y, 2.0f - axy);
  conv.f = y;
  conv.u |= sign;
  return conv.f;
}

static inline float FastFloatDiv(float num, float den) {
  return FastFloatMul(num, FastFloatInv(den));
}

static inline float FastFloatSqrt(float x) {
  union { float f; uint32_t u; } conv;
  conv.f = x;
  if (conv.u <= 0x007fffff) return 0.0f;
  conv.u = 0x1fbd1df5 + (conv.u >> 1);
  return conv.f;
}

#ifdef _MSC_VER /* visual c++ */
#define ALIGN16_BEG __declspec(align(16))
#define ALIGN16_END
#else /* gcc or icc */
#define ALIGN16_BEG
#define ALIGN16_END __attribute__((aligned(16)))
#endif

enum class Aec3Optimization { kNone, kSse2, kAvx2, kNeon };

constexpr int kNumBlocksPerSecond = 250;

constexpr int kMetricsReportingIntervalBlocks = 10 * kNumBlocksPerSecond;
constexpr int kMetricsComputationBlocks = 3;
constexpr int kMetricsCollectionBlocks =
    kMetricsReportingIntervalBlocks - kMetricsComputationBlocks;

constexpr size_t kFftLengthBy2 = 64;
constexpr size_t kFftLengthBy2Plus1 = kFftLengthBy2 + 1;
constexpr size_t kFftLengthBy2Minus1 = kFftLengthBy2 - 1;
constexpr size_t kFftLength = 2 * kFftLengthBy2;
constexpr size_t kFftLengthBy2Log2 = 6;

constexpr int kRenderTransferQueueSizeFrames = 100;

constexpr size_t kMaxNumBands = 3;
constexpr size_t kFrameSize = 160;
constexpr size_t kSubFrameLength = kFrameSize / 2;

constexpr size_t kBlockSize = kFftLengthBy2;
constexpr size_t kBlockSizeLog2 = kFftLengthBy2Log2;

constexpr size_t kExtendedBlockSize = 2 * kFftLengthBy2;
constexpr size_t kMatchedFilterWindowSizeSubBlocks = 32;
constexpr size_t kMatchedFilterAlignmentShiftSizeSubBlocks =
    kMatchedFilterWindowSizeSubBlocks * 3 / 4;

// TODO(peah): Integrate this with how it is done inside audio_processing_impl.
constexpr size_t NumBandsForRate(int sample_rate_hz) {
  return static_cast<size_t>(sample_rate_hz / 16000);
}

constexpr bool ValidFullBandRate(int sample_rate_hz) {
  return sample_rate_hz == 16000 || sample_rate_hz == 32000 ||
         sample_rate_hz == 48000;
}

constexpr int GetTimeDomainLength(int filter_length_blocks) {
  return filter_length_blocks * kFftLengthBy2;
}

constexpr size_t GetDownSampledBufferSize(size_t down_sampling_factor,
                                          size_t num_matched_filters) {
  return kBlockSize / down_sampling_factor *
         (kMatchedFilterAlignmentShiftSizeSubBlocks * num_matched_filters +
          kMatchedFilterWindowSizeSubBlocks + 1);
}

constexpr size_t GetRenderDelayBufferSize(size_t down_sampling_factor,
                                          size_t num_matched_filters,
                                          size_t filter_length_blocks) {
  return GetDownSampledBufferSize(down_sampling_factor, num_matched_filters) /
             (kBlockSize / down_sampling_factor) +
         filter_length_blocks + 1;
}

// Detects what kind of optimizations to use for the code.
Aec3Optimization DetectOptimization();

// Computes the log2 of the input in a fast an approximate manner.
float FastApproxLog2f(float in);

// Returns dB from a power quantity expressed in log2.
float Log2TodB(float in_log2);

static_assert(1 << kBlockSizeLog2 == kBlockSize,
              "Proper number of shifts for blocksize");

static_assert(1 << kFftLengthBy2Log2 == kFftLengthBy2,
              "Proper number of shifts for the fft length");

static_assert(1 == NumBandsForRate(16000), "Number of bands for 16 kHz");
static_assert(2 == NumBandsForRate(32000), "Number of bands for 32 kHz");
static_assert(3 == NumBandsForRate(48000), "Number of bands for 48 kHz");

static_assert(ValidFullBandRate(16000),
              "Test that 16 kHz is a valid sample rate");
static_assert(ValidFullBandRate(32000),
              "Test that 32 kHz is a valid sample rate");
static_assert(ValidFullBandRate(48000),
              "Test that 48 kHz is a valid sample rate");
static_assert(!ValidFullBandRate(8001),
              "Test that 8001 Hz is not a valid sample rate");

}  // namespace webrtc

#endif  // MODULES_AUDIO_PROCESSING_AEC3_AEC3_COMMON_H_
