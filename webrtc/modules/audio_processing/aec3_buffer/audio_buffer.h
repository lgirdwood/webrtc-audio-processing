#ifndef MODULES_AUDIO_PROCESSING_AUDIO_BUFFER_H_
#define MODULES_AUDIO_PROCESSING_AUDIO_BUFFER_H_

#include <stddef.h>
#include <stdint.h>
#include <vector>

namespace webrtc {

class AudioBuffer {
 public:
  static constexpr int kSplitBandSize = 160;

  AudioBuffer(size_t input_rate,
              size_t input_num_channels,
              size_t buffer_rate,
              size_t buffer_num_channels,
              size_t output_rate,
              size_t output_num_channels);

  virtual ~AudioBuffer();

  size_t num_channels() const { return num_channels_; }
  size_t num_bands() const { return num_bands_; }
  size_t num_frames() const { return num_bands_ * kSplitBandSize; }
  size_t num_frames_per_band() const { return kSplitBandSize; }

  float* const* split_bands(size_t channel) {
    return band_ptrs_[channel].data();
  }
  const float* const* split_bands_const(size_t channel) const {
    return band_ptrs_[channel].data();
  }

  float* const* channels() {
    return channel_ptrs_.data();
  }
  const float* const* channels_const() const {
    return channel_ptrs_.data();
  }

 private:
  size_t num_channels_;
  size_t num_bands_;
  std::vector<std::vector<std::vector<float>>> data_;
  std::vector<std::vector<float*>> band_ptrs_;
  std::vector<float*> channel_ptrs_;
};

}  // namespace webrtc

#endif  // MODULES_AUDIO_PROCESSING_AUDIO_BUFFER_H_
