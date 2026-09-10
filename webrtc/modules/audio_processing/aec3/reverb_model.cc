/*
 *  Copyright (c) 2018 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "modules/audio_processing/aec3/reverb_model.h"

#include <stddef.h>

#include <algorithm>
#include <functional>

#include "api/array_view.h"
#include "modules/audio_processing/aec3/aec3_common.h"

namespace webrtc {

ReverbModel::ReverbModel() {
  Reset();
}

ReverbModel::~ReverbModel() = default;

void ReverbModel::Reset() {
  reverb_.fill(0.);
}

void ReverbModel::UpdateReverbNoFreqShaping(
    rtc::ArrayView<const float> power_spectrum,
    float power_spectrum_scaling,
    float reverb_decay) {
  if (reverb_decay > 0) {
    float combined_scale = FastFloatMul(power_spectrum_scaling, reverb_decay);
    for (size_t k = 0; k < power_spectrum.size(); ++k) {
      float a = FastFloatMul(reverb_[k], reverb_decay);
      float b = FastFloatMul(power_spectrum[k], combined_scale);
      reverb_[k] = FastFloatAddPos(a, b);
    }
  }
}

void ReverbModel::UpdateReverb(
    rtc::ArrayView<const float> power_spectrum,
    rtc::ArrayView<const float> power_spectrum_scaling,
    float reverb_decay) {
  if (reverb_decay > 0) {
    for (size_t k = 0; k < power_spectrum.size(); ++k) {
      float a = FastFloatMul(reverb_[k], reverb_decay);
      float b = FastFloatMul(FastFloatMul(power_spectrum[k], power_spectrum_scaling[k]), reverb_decay);
      reverb_[k] = FastFloatAddPos(a, b);
    }
  }
}

}  // namespace webrtc
