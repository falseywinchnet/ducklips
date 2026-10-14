#include "voice/wave_file.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>

namespace voice {

namespace {

void put_u32(std::FILE* file, std::uint32_t value) {
  unsigned char bytes[4] = {static_cast<unsigned char>(value), static_cast<unsigned char>(value >> 8U),
                            static_cast<unsigned char>(value >> 16U), static_cast<unsigned char>(value >> 24U)};
  std::fwrite(bytes, 1, 4, file);
}

void put_u16(std::FILE* file, std::uint16_t value) {
  unsigned char bytes[2] = {static_cast<unsigned char>(value), static_cast<unsigned char>(value >> 8U)};
  std::fwrite(bytes, 1, 2, file);
}

}  // namespace

bool write_wave(const std::string& path, const std::vector<float>& samples, int rate) {
  std::FILE* file = std::fopen(path.c_str(), "wb");
  if (file == nullptr) {
    return false;
  }
  std::uint32_t data_bytes = static_cast<std::uint32_t>(samples.size() * 2U);
  std::fwrite("RIFF", 1, 4, file);
  put_u32(file, 36U + data_bytes);
  std::fwrite("WAVEfmt ", 1, 8, file);
  put_u32(file, 16U);
  put_u16(file, 1U);
  put_u16(file, 1U);
  put_u32(file, static_cast<std::uint32_t>(rate));
  put_u32(file, static_cast<std::uint32_t>(rate) * 2U);
  put_u16(file, 2U);
  put_u16(file, 16U);
  std::fwrite("data", 1, 4, file);
  put_u32(file, data_bytes);
  for (float sample : samples) {
    float clamped = std::clamp(sample, -1.0F, 1.0F);
    std::int16_t value = static_cast<std::int16_t>(std::lround(clamped * 32767.0F));
    put_u16(file, static_cast<std::uint16_t>(value));
  }
  bool good = std::ferror(file) == 0;
  return std::fclose(file) == 0 && good;
}

}  // namespace voice
