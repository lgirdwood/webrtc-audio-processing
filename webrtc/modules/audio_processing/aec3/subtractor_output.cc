/*
 *  Copyright (c) 2018 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "modules/audio_processing/aec3/subtractor_output.h"

#include <numeric>

namespace webrtc {

SubtractorOutput::SubtractorOutput() = default;
SubtractorOutput::~SubtractorOutput() = default;

void SubtractorOutput::Reset() {
  s_refined.fill(0.f);
  s_coarse.fill(0.f);
  e_refined.fill(0.f);
  e_coarse.fill(0.f);
  E_refined.re.fill(0.f);
  E_refined.im.fill(0.f);
  E2_refined.fill(0.f);
  E2_coarse.fill(0.f);
  e2_refined = 0.f;
  e2_coarse = 0.f;
  s2_refined = 0.f;
  s2_coarse = 0.f;
  y2 = 0.f;
}

void SubtractorOutput::ComputeMetrics(rtc::ArrayView<const float> y) {
  float y2_sum = 0.f;
  float e2_sum = 0.f;
  float s2_sum = 0.f;
  float s_max = 0.f;
  for (size_t i = 0; i < kBlockSize; ++i) {
    y2_sum = FastFloatAddPos(y2_sum, FastFloatSqr(y[i]));
    e2_sum = FastFloatAddPos(e2_sum, FastFloatSqr(e_refined[i]));
    float s_val = s_refined[i];
    s2_sum = FastFloatAddPos(s2_sum, FastFloatSqr(s_val));
    float abs_s = std::abs(s_val);
    if (abs_s > s_max) s_max = abs_s;
  }
  y2 = y2_sum;
  e2_refined = e2_sum;
  s2_refined = s2_sum;
  e2_coarse = e2_sum;
  s2_coarse = s2_sum;
  s_refined_max_abs = s_max;
  s_coarse_max_abs = s_max;
}

}  // namespace webrtc
