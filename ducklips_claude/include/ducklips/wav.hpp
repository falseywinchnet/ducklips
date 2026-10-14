#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace ducklips {

// Interleaved samples in [-1, 1].
struct AudioBuffer final {
    std::vector<float> samples{};
    std::size_t channels{2};
    double sample_rate{48000.0};
};

// Writes little-endian RIFF/WAVE IEEE float32, the format GUI.Forms
// AudioClip::load_wav accepts for stereo 48 kHz. Returns false on I/O failure.
[[nodiscard]] bool write_wav_float(const std::string& path, const AudioBuffer& audio);
// Reads PCM16, PCM24 or float32 WAVE files. Returns false and leaves output
// unchanged when the file is unreadable or in another format.
[[nodiscard]] bool read_wav(const std::string& path, AudioBuffer& output);

}
