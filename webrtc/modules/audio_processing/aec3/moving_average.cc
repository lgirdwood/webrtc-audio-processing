
/*
 *  Copyright (c) 2018 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "modules/audio_processing/aec3/moving_average.h"

#include <algorithm>
#include <functional>

#include "modules/audio_processing/aec3/aec3_common.h"
#include "rtc_base/checks.h"

namespace webrtc {
namespace aec3 {

MovingAverage::MovingAverage(size_t num_elem, size_t mem_len)
    : num_elem_(num_elem),
      mem_len_(mem_len - 1),
      scaling_(1.0f / static_cast<float>(mem_len)),
      memory_(num_elem * mem_len_, 0.f),
      mem_index_(0) {
  RTC_DCHECK(num_elem_ > 0);
  RTC_DCHECK(mem_len > 0);
}

MovingAverage::~MovingAverage() = default;

void MovingAverage::Average(rtc::ArrayView<const float> input,
                            rtc::ArrayView<float> output) {
  RTC_DCHECK_EQ(input.size(), num_elem_);
  RTC_DCHECK_EQ(output.size(), num_elem_);

  if (mem_len_ == 0) {
    std::copy(input.begin(), input.end(), output.begin());
    return;
  }

  const float* in = input.data();
  float* out = output.data();
  float* cur_mem = &memory_[mem_index_ * num_elem_];

  if (mem_len_ == 3) {
    const float* m0 = &memory_[0];
    const float* m1 = &memory_[num_elem_];
    const float* m2 = &memory_[num_elem_ * 2];
    for (size_t k = 0; k < num_elem_; ++k) {
      float v = in[k];
      float s = FastFloatAddPos(FastFloatAddPos(FastFloatAddPos(v, m0[k]), m1[k]), m2[k]);
      cur_mem[k] = v;
      out[k] = FastFloatMul(s, scaling_);
    }
  } else {
    for (size_t k = 0; k < num_elem_; ++k) {
      float v = in[k];
      float s = v;
      for (size_t m = 0; m < mem_len_; ++m) {
        s = FastFloatAddPos(s, memory_[m * num_elem_ + k]);
      }
      cur_mem[k] = v;
      out[k] = FastFloatMul(s, scaling_);
    }
  }

  mem_index_ = (mem_index_ + 1);
  if (mem_index_ >= mem_len_) {
    mem_index_ = 0;
  }
}

}  // namespace aec3
}  // namespace webrtc
