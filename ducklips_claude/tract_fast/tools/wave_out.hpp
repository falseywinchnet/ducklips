#pragma once

// Output helpers shared by the tools: the 96 -> 48 kHz decimation filter and
// a 16-bit mono WAV writer.

#include "tract/units.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <vector>

namespace {

// Two Butterworth biquads at 18 kHz for the 96 -> 48 kHz decimation.
struct Biquad {
  double b0 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
  double z1 = 0.0;
  double z2 = 0.0;
};

inline Biquad lowpass(double corner, double rate, double q) {
  double w = std::tan(tract::kPi * corner / rate);
  double norm = 1.0 / (1.0 + w / q + w * w);
  Biquad filter;
  filter.b0 = w * w * norm;
  filter.b1 = 2.0 * filter.b0;
  filter.b2 = filter.b0;
  filter.a1 = 2.0 * (w * w - 1.0) * norm;
  filter.a2 = (1.0 - w / q + w * w) * norm;
  return filter;
}

inline double run(Biquad& filter, double input) {
  double output = filter.b0 * input + filter.z1;
  filter.z1 = filter.b1 * input - filter.a1 * output + filter.z2;
  filter.z2 = filter.b2 * input - filter.a2 * output;
  return output;
}

inline void write_wave(const char* path, const std::vector<float>& samples, int rate) {
  std::ofstream out(path, std::ios::binary);
  std::uint32_t data_bytes = static_cast<std::uint32_t>(samples.size() * 2);
  out.write("RIFF", 4);
  std::uint32_t chunk = 36 + data_bytes;
  out.write(reinterpret_cast<const char*>(&chunk), 4);
  out.write("WAVEfmt ", 8);
  std::uint32_t format_size = 16;
  std::uint16_t format = 1;
  std::uint16_t channels = 1;
  std::uint32_t sample_rate = static_cast<std::uint32_t>(rate);
  std::uint32_t byte_rate = sample_rate * 2;
  std::uint16_t align = 2;
  std::uint16_t bits = 16;
  out.write(reinterpret_cast<const char*>(&format_size), 4);
  out.write(reinterpret_cast<const char*>(&format), 2);
  out.write(reinterpret_cast<const char*>(&channels), 2);
  out.write(reinterpret_cast<const char*>(&sample_rate), 4);
  out.write(reinterpret_cast<const char*>(&byte_rate), 4);
  out.write(reinterpret_cast<const char*>(&align), 2);
  out.write(reinterpret_cast<const char*>(&bits), 2);
  out.write("data", 4);
  out.write(reinterpret_cast<const char*>(&data_bytes), 4);
  for (float sample : samples) {
    double clipped = std::fmax(-1.0, std::fmin(1.0, static_cast<double>(sample)));
    std::int16_t value = static_cast<std::int16_t>(std::lround(clipped * 32767.0));
    out.write(reinterpret_cast<const char*>(&value), 2);
  }
}

}  // namespace
