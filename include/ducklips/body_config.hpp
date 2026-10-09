#pragma once
#include "ducklips/voice.hpp"
#include <istream>
#include <string>

namespace ducklips {
// Strict, transactional body configuration. SI units, key=value, # comments.
// Optional preset must be the first setting. Duplicate/unknown/nonfinite or
// out-of-range fields fail with a line number; output stays unchanged.
[[nodiscard]] bool read_body_config(std::istream& input, VoiceProfile& output, std::string& error);
[[nodiscard]] bool load_body_config(const std::string& path, VoiceProfile& output, std::string& error);
}
