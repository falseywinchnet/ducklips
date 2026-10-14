#pragma once

#include <string>
#include <vector>

namespace voice {

// Writes mono 16-bit PCM at the given rate. Returns false on failure.
[[nodiscard]] bool write_wave(const std::string& path, const std::vector<float>& samples, int rate);

}  // namespace voice
