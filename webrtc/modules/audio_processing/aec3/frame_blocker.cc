/*
 *  Copyright (c) 2016 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "modules/audio_processing/aec3/frame_blocker.h"

#include <cstring>
#include "modules/audio_processing/aec3/aec3_common.h"
#include "rtc_base/checks.h"

namespace webrtc {

FrameBlocker::FrameBlocker(size_t num_bands, size_t num_channels)
    : num_bands_(num_bands),
      num_channels_(num_channels),
      buffer_data_(num_bands * num_channels * kBlockSize, 0.f),
      buffer_size_(num_bands * num_channels, 0) {
  RTC_DCHECK_LT(0, num_bands);
  RTC_DCHECK_LT(0, num_channels);
}

FrameBlocker::~FrameBlocker() = default;

void FrameBlocker::InsertSubFrameAndExtractBlock(
    const std::vector<std::vector<rtc::ArrayView<float>>>& sub_frame,
    Block* block) {
  RTC_DCHECK(block);
  RTC_DCHECK_EQ(num_bands_, block->NumBands());
  RTC_DCHECK_EQ(num_bands_, sub_frame.size());
  for (size_t band = 0; band < num_bands_; ++band) {
    RTC_DCHECK_EQ(num_channels_, block->NumChannels());
    RTC_DCHECK_EQ(num_channels_, sub_frame[band].size());
    for (size_t channel = 0; channel < num_channels_; ++channel) {
      size_t idx = band * num_channels_ + channel;
      RTC_DCHECK_GE(kBlockSize - 16, buffer_size_[idx]);
      RTC_DCHECK_EQ(kSubFrameLength, sub_frame[band][channel].size());

      float* buf = &buffer_data_[idx * kBlockSize];
      const size_t cur_sz = buffer_size_[idx];
      const size_t samples_to_block = kBlockSize - cur_sz;
      const float* sub_frame_src = sub_frame[band][channel].data();
      float* block_dst = block->View(band, channel).data();

      std::memcpy(block_dst, buf, cur_sz * sizeof(float));
      std::memcpy(block_dst + cur_sz, sub_frame_src, samples_to_block * sizeof(float));
      const size_t rem = kSubFrameLength - samples_to_block;
      std::memcpy(buf, sub_frame_src + samples_to_block, rem * sizeof(float));
      buffer_size_[idx] = rem;
    }
  }
}

bool FrameBlocker::IsBlockAvailable() const {
  return kBlockSize == buffer_size_[0];
}

void FrameBlocker::ExtractBlock(Block* block) {
  RTC_DCHECK(block);
  RTC_DCHECK_EQ(num_bands_, block->NumBands());
  RTC_DCHECK_EQ(num_channels_, block->NumChannels());
  RTC_DCHECK(IsBlockAvailable());
  for (size_t band = 0; band < num_bands_; ++band) {
    for (size_t channel = 0; channel < num_channels_; ++channel) {
      size_t idx = band * num_channels_ + channel;
      RTC_DCHECK_EQ(kBlockSize, buffer_size_[idx]);
      float* buf = &buffer_data_[idx * kBlockSize];
      std::memcpy(block->View(band, channel).data(), buf, kBlockSize * sizeof(float));
      buffer_size_[idx] = 0;
    }
  }
}

}  // namespace webrtc
