// Render a sheet tube (Astra's engine, ported): a single design with its
// inlet reservoir and open outlet, or a connected .system.
// usage: sheet_render DESIGN.txt|SYSTEM.system SECONDS OUT.wav [rate]
// The physics runs at `rate` (default 24000, Astra's); the observed outlet
// pressure difference is written as 32-bit float at 48 kHz (raw values, as
// her render: (p - p_prev) * rate / 48000 / 10, linear interpolation).

#include "tract/sheet_tube.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace {

void write_float_wave(const char* path, const std::vector<float>& samples, int rate) {
  std::ofstream out(path, std::ios::binary);
  std::uint32_t data_bytes = static_cast<std::uint32_t>(samples.size() * 4);
  std::uint32_t chunk = 36 + data_bytes;
  std::uint32_t format_size = 16;
  std::uint16_t format = 3;
  std::uint16_t channels = 1;
  std::uint32_t sample_rate = static_cast<std::uint32_t>(rate);
  std::uint32_t byte_rate = sample_rate * 4;
  std::uint16_t align = 4;
  std::uint16_t bits = 32;
  out.write("RIFF", 4);
  out.write(reinterpret_cast<const char*>(&chunk), 4);
  out.write("WAVEfmt ", 8);
  out.write(reinterpret_cast<const char*>(&format_size), 4);
  out.write(reinterpret_cast<const char*>(&format), 2);
  out.write(reinterpret_cast<const char*>(&channels), 2);
  out.write(reinterpret_cast<const char*>(&sample_rate), 4);
  out.write(reinterpret_cast<const char*>(&byte_rate), 4);
  out.write(reinterpret_cast<const char*>(&align), 2);
  out.write(reinterpret_cast<const char*>(&bits), 2);
  out.write("data", 4);
  out.write(reinterpret_cast<const char*>(&data_bytes), 4);
  out.write(reinterpret_cast<const char*>(samples.data()), static_cast<std::streamsize>(data_bytes));
}

}  // namespace

int main(int count, char** values) {
  if (count < 4) {
    std::fprintf(stderr, "usage: sheet_render DESIGN.txt|SYSTEM.system SECONDS OUT.wav [rate]\n");
    return 2;
  }
  std::string path = values[1];
  double seconds = std::strtod(values[2], nullptr);
  double rate = count > 4 ? std::strtod(values[4], nullptr) : 24000.0;
  bool is_system = path.size() > 7 && path.substr(path.size() - 7) == ".system";
  tract::SheetSystem system(rate);
  std::vector<tract::SheetTube> single;
  bool radiating = false;
  if (is_system) {
    if (!system.load(path)) {
      std::fprintf(stderr, "cannot load %s\n", path.c_str());
      return 1;
    }
  } else {
    tract::SheetDesign design;
    if (!tract::read_sheet_design(path, design)) {
      std::fprintf(stderr, "cannot read %s\n", path.c_str());
      return 1;
    }
    single.emplace_back(design, rate);
    radiating = design.radiating;
  }
  long steps = static_cast<long>(seconds * rate);
  int per_output = static_cast<int>(std::lround(48000.0 / rate));
  std::vector<float> samples;
  samples.reserve(static_cast<std::size_t>(steps * per_output + 2));
  double previous_outlet = 0.0;
  double last = 0.0;
  std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
  for (long step = 0; step < steps; ++step) {
    double outlet = 0.0;
    if (is_system) {
      outlet = system.step();
      if (system.fault()) {
        std::fprintf(stderr, "fault at %.4f s\n", step / rate);
        return 1;
      }
    } else {
      tract::SheetTube& tube = single.front();
      if (!tube.step()) {
        std::fprintf(stderr, "fault at %.4f s\n", step / rate);
        return 1;
      }
      outlet = radiating ? tube.radiated_flow() * 1.0e3 : tube.pressure(tube.count() - 1);
    }
    // Radiating: the rate of change of the radiated flow (the sound); else
    // Astra's observed outlet pressure difference.
    double value = (outlet - previous_outlet) * rate / 48000.0 / 10.0;
    previous_outlet = outlet;
    for (int k = 0; k < per_output; ++k) {
      double share = static_cast<double>(k) / per_output;
      samples.push_back(static_cast<float>(last + share * (value - last)));
    }
    last = value;
  }
  double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
  write_float_wave(values[3], samples, 48000);
  long contacts = is_system ? 0 : single.front().contacts();
  std::fprintf(stderr, "%s: %.2f s in %.2f s (%.1fx real time), contacts %ld\n", values[3], seconds, elapsed,
               seconds / elapsed, contacts);
  return 0;
}
