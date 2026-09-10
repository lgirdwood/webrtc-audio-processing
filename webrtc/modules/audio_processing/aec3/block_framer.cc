/*
 *  Copyright (c) 2016 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "modules/audio_processing/aec3/block_framer.h"

#include <cstring>
#include <algorithm>

#include "modules/audio_processing/aec3/aec3_common.h"
#include "rtc_base/checks.h"

namespace webrtc {

BlockFramer::BlockFramer(size_t num_bands, size_t num_channels)
    : num_bands_(num_bands),
      num_channels_(num_channels),
      buffer_data_(num_bands * num_channels * kBlockSize, 0.f),
      buffer_size_(num_bands * num_channels, kBlockSize) {
  RTC_DCHECK_LT(0, num_bands);
  RTC_DCHECK_LT(0, num_channels);
}

BlockFramer::~BlockFramer() = default;

// All the constants are chosen so that the buffer is either empty or has enough
// samples for InsertBlockAndExtractSubFrame to produce a frame. In order to
// achieve this, the InsertBlockAndExtractSubFrame and InsertBlock methods need
// to be called in the correct order.
void BlockFramer::InsertBlock(const Block& block) {
  RTC_DCHECK_EQ(num_bands_, block.NumBands());
  RTC_DCHECK_EQ(num_channels_, block.NumChannels());
  for (size_t band = 0; band < num_bands_; ++band) {
    for (size_t channel = 0; channel < num_channels_; ++channel) {
      size_t idx = band * num_channels_ + channel;
      RTC_DCHECK_EQ(0, buffer_size_[idx]);
      float* buf = &buffer_data_[idx * kBlockSize];
      const float* block_src = block.View(band, channel).data();
      std::memcpy(buf, block_src, kBlockSize * sizeof(float));
      buffer_size_[idx] = kBlockSize;
    }
  }
}

void BlockFramer::InsertBlockAndExtractSubFrame(
    const Block& block,
    std::vector<std::vector<rtc::ArrayView<float>>>* sub_frame) {
  RTC_DCHECK(sub_frame);
  RTC_DCHECK_EQ(num_bands_, block.NumBands());
  RTC_DCHECK_EQ(num_channels_, block.NumChannels());
  RTC_DCHECK_EQ(num_bands_, sub_frame->size());
  for (size_t band = 0; band < num_bands_; ++band) {
    RTC_DCHECK_EQ(num_channels_, (*sub_frame)[0].size());
    for (size_t channel = 0; channel < num_channels_; ++channel) {
      size_t idx = band * num_channels_ + channel;
      RTC_DCHECK_LE(kSubFrameLength, buffer_size_[idx] + kBlockSize);
      RTC_DCHECK_GE(kBlockSize, buffer_size_[idx]);
      RTC_DCHECK_EQ(kSubFrameLength, (*sub_frame)[band][channel].size());

      float* buf = &buffer_data_[idx * kBlockSize];
      const size_t cur_sz = buffer_size_[idx];
      const size_t samples_to_frame = kSubFrameLength - cur_sz;
      float* sub_frame_dst = (*sub_frame)[band][channel].data();
      const float* block_src = block.View(band, channel).data();

      std::memcpy(sub_frame_dst, buf, cur_sz * sizeof(float));
      std::memcpy(sub_frame_dst + cur_sz, block_src,
                  samples_to_frame * sizeof(float));
      const size_t rem = kBlockSize - samples_to_frame;
      std::memcpy(buf, block_src + samples_to_frame,
                  rem * sizeof(float));
      buffer_size_[idx] = rem;
    }
  }
}

}  // namespace webrtc
