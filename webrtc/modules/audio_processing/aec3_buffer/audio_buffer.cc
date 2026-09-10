#include "modules/audio_processing/aec3_buffer/audio_buffer.h"

namespace webrtc {

AudioBuffer::AudioBuffer(size_t input_rate,
                         size_t input_num_channels,
                         size_t buffer_rate,
                         size_t buffer_num_channels,
                         size_t output_rate,
                         size_t output_num_channels)
    : num_channels_(buffer_num_channels),
      num_bands_(buffer_rate >= 48000 ? 3 : (buffer_rate >= 32000 ? 2 : 1)) {
  data_.resize(num_channels_);
  band_ptrs_.resize(num_channels_);
  channel_ptrs_.resize(num_channels_);

  for (size_t ch = 0; ch < num_channels_; ++ch) {
    data_[ch].resize(num_bands_);
    band_ptrs_[ch].resize(num_bands_);
    for (size_t b = 0; b < num_bands_; ++b) {
      data_[ch][b].resize(kSplitBandSize, 0.0f);
      band_ptrs_[ch][b] = data_[ch][b].data();
    }
    channel_ptrs_[ch] = data_[ch][0].data();
  }
}

AudioBuffer::~AudioBuffer() = default;

}  // namespace webrtc
