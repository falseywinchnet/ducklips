// Render a sheet-tube design whose posture is HELD, not fixed: the breath
// (reservoir pressure), the cord's adduction (its ring support's offset) and
// stretch (its elasticity) and every support section's ring are coordinates
// of a SomaticHold (Astra's motor pools, sensing delay and relaxing strain
// memory), ticked at 1 kHz; the tissue's own strain memory relaxes inside the
// sheet. A human cannot hold a posture perfectly; nothing here is periodic.
//
// usage: sheet_hold DESIGN.txt SECONDS OUT.wav [options]
//   --rate R          physics rate (default 48000)
//   --no-hold         the design's static posture (the comparison)
//   --noise N         retained motor-pool fluctuation (default 0.08)
//   --common C        common-drive share of effort (default 0.02)
//   --relaxation G    share of posture in the relaxing memory (default 0.15)
//   --tissue F        tissue strain-relaxation share (default 0.1; 0 = off)
//   --cord-rest M     how far the relaxed cord's ring stands open, m (default 1e-3)
//   --cord-scale M    adduction span per unit effort, m (default 1e-3)
//   --stretch-rest E  the relaxed cord's elasticity factor (default 0.6)
//   --stretch-scale S elasticity span per unit effort (default 0.4)
//   --neutral A       the relaxed tract's area, m^2 (default 3e-4)
//   --support-scale M tract ring span per unit effort, m (default 2e-3)
//   --breath-scale P  breath span per unit effort, Pa (default 600)
// Effort is measured from the relaxed posture (an open cord, slack stretch,
// a neutral tract, no breath): a coordinate held far from it works harder
// and wanders more. Every coordinate starts relaxed and is reached.
//   --keep-wander     keep the design's reservoir wander and variation
//   --trace FILE      per-millisecond held values (breath, adduction, stretch, first ring, cord gap)
//   --seed N
// Writes 32-bit float at 48 kHz: the rate of change of the radiated flow.

#include "tract/sheet_tube.hpp"
#include "tract/somatic_hold.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

struct Options {
  double rate = 48000.0;
  bool hold = true;
  double noise = 0.08;
  double relaxation = 0.15;
  double common = 0.02;
  double tissue = 0.1;
  double cord_rest = 1e-3;
  double cord_scale = 1e-3;
  double stretch_rest = 0.6;
  double stretch_scale = 0.4;
  double neutral = 3e-4;
  double support_scale = 2e-3;
  double breath_scale = 600.0;
  bool keep_wander = false;
  bool churn = false;  // with --no-hold: re-set the cord's unchanged offset every ms (a test)
  const char* trace = nullptr;
  std::uint32_t seed = 9001U;
};

bool read_options(int count, char** values, Options& options) {
  for (int i = 4; i < count; ++i) {
    const char* name = values[i];
    bool has_value = i + 1 < count;
    if (std::strcmp(name, "--no-hold") == 0) {
      options.hold = false;
    } else if (std::strcmp(name, "--churn") == 0) {
      options.churn = true;
    } else if (std::strcmp(name, "--keep-wander") == 0) {
      options.keep_wander = true;
    } else if (has_value && std::strcmp(name, "--rate") == 0) {
      options.rate = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--noise") == 0) {
      options.noise = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--common") == 0) {
      options.common = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--relaxation") == 0) {
      options.relaxation = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--tissue") == 0) {
      options.tissue = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--cord-rest") == 0) {
      options.cord_rest = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--stretch-rest") == 0) {
      options.stretch_rest = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--neutral") == 0) {
      options.neutral = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--cord-scale") == 0) {
      options.cord_scale = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--stretch-scale") == 0) {
      options.stretch_scale = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--support-scale") == 0) {
      options.support_scale = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--breath-scale") == 0) {
      options.breath_scale = std::strtod(values[++i], nullptr);
    } else if (has_value && std::strcmp(name, "--trace") == 0) {
      options.trace = values[++i];
    } else if (has_value && std::strcmp(name, "--seed") == 0) {
      options.seed = static_cast<std::uint32_t>(std::strtoul(values[++i], nullptr, 10));
    } else {
      std::fprintf(stderr, "unknown option %s\n", name);
      return false;
    }
  }
  return true;
}

// Running mean and spread of a held value, for the report.
struct Spread {
  double sum = 0.0;
  double square = 0.0;
  long count = 0;
};

void add(Spread& spread, double value) {
  spread.sum += value;
  spread.square += value * value;
  ++spread.count;
}

double mean(const Spread& spread) {
  return spread.count > 0 ? spread.sum / spread.count : 0.0;
}

double deviation(const Spread& spread) {
  if (spread.count < 2) {
    return 0.0;
  }
  double mean = spread.sum / spread.count;
  return std::sqrt(std::fmax(0.0, spread.square / spread.count - mean * mean));
}

// The voice source from the glottal flow (cm^3/s): mean and peak flow, the
// share of time closed (below 5% of the peak), and the fastest closing
// (maximum flow declination rate, cm^3/s per ms).
void report_source(const std::vector<float>& flow, double rate) {
  if (flow.size() < 2) {
    return;
  }
  std::vector<float> sorted(flow);
  std::sort(sorted.begin(), sorted.end());
  double peak = sorted[sorted.size() * 99 / 100];
  double sum = 0.0;
  long closed = 0;
  double steepest = 0.0;
  for (std::size_t i = 0; i < flow.size(); ++i) {
    sum += flow[i];
    if (flow[i] < 0.05 * peak) {
      ++closed;
    }
    if (i > 0) {
      steepest = std::max(steepest, static_cast<double>(flow[i - 1] - flow[i]) * rate / 1000.0);
    }
  }
  std::fprintf(stderr, "source: mean flow %.0f, peak %.0f cm3/s, closed %.0f%%, MFDR %.0f cm3/s/ms\n",
               sum / flow.size(), peak, 100.0 * closed / flow.size(), steepest);
}

}  // namespace

int main(int count, char** values) {
  if (count < 4) {
    std::fprintf(stderr, "usage: sheet_hold DESIGN.txt SECONDS OUT.wav [options]\n");
    return 2;
  }
  Options options;
  if (!read_options(count, values, options)) {
    return 2;
  }
  tract::SheetDesign design;
  if (!tract::read_sheet_design(values[1], design)) {
    std::fprintf(stderr, "cannot read %s\n", values[1]);
    return 1;
  }
  double breath_target = design.pressure;
  if (options.hold && !options.keep_wander) {
    design.wander = 0.0;
    design.variation = 0.0;
  }
  if (options.hold) {
    design.pressure = 0.0;  // the held breath sets it
  }
  double seconds = std::strtod(values[2], nullptr);
  double rate = options.rate;
  tract::SheetTube tube(design, rate);
  tube.set_outlet(design.radiating ? tract::SheetTube::End::radiating : tract::SheetTube::End::open);
  if (options.tissue > 0.0) {
    tube.set_strain_relaxation(options.tissue, 0.1, 0.6, options.seed + 17U);
  }

  tract::HoldSettings settings;
  settings.noise = options.noise;
  settings.relaxation = options.relaxation;
  settings.common = options.common;
  tract::SomaticHold hold;
  hold.reset(options.seed, settings);
  int breath = hold.add(0.0, options.breath_scale, 0.0, 2500.0, 4);
  hold.set_target(breath, breath_target, 20.0);
  int cord_band = -1;
  for (int b = 0; b < static_cast<int>(design.bands.size()); ++b) {
    if (design.bands[static_cast<std::size_t>(b)].role == "cord") {
      cord_band = b;
    }
  }
  int adduction = -1;
  int stretch = -1;
  if (cord_band >= 0) {
    double offset = design.bands[static_cast<std::size_t>(cord_band)].offset;
    adduction = hold.add(offset + options.cord_rest, options.cord_scale, offset - 0.002,
                         offset + options.cord_rest + 0.002, 3);
    hold.set_target(adduction, offset, 40.0);
    stretch = hold.add(options.stretch_rest, options.stretch_scale, 0.3, 3.0, 3);
    hold.set_target(stretch, 1.0, 40.0);
  }
  // Every support section, and the lips, hold their rings.
  std::vector<int> ring_bands;
  std::vector<int> ring_coordinates;
  for (int b = 0; b < static_cast<int>(design.bands.size()); ++b) {
    const tract::SheetBand& band = design.bands[static_cast<std::size_t>(b)];
    if (band.role != "support" && band.role != "lips") {
      continue;
    }
    int group = band.role == "lips" ? 1 : 0;
    // Relaxed, the ring stands at the neutral tract's radius.
    double relaxed = std::sqrt(options.neutral / 3.14159265358979323846) - band.radius;
    int index = hold.add(relaxed, options.support_scale, band.offset - 0.6 * band.radius,
                         std::fmax(relaxed, band.offset) + 0.01, group);
    if (index < 0) {
      break;
    }
    hold.set_target(index, band.offset, 40.0);
    ring_bands.push_back(b);
    ring_coordinates.push_back(index);
  }

  long steps = static_cast<long>(seconds * rate);
  int per_output = static_cast<int>(std::lround(48000.0 / rate));
  int per_tick = static_cast<int>(std::lround(rate * tract::SomaticHold::kTick));
  std::vector<float> samples;
  samples.reserve(static_cast<std::size_t>(steps * per_output + 2));
  double previous_outlet = 0.0;
  double last = 0.0;
  Spread breath_spread;
  Spread adduction_spread;
  Spread stretch_spread;
  Spread ring_spread;
  // The glottal flow after 0.5 s, for the source report.
  std::vector<float> glottal;
  glottal.reserve(static_cast<std::size_t>(seconds * rate) + 1);
  double smallest_gap = 1.0;  // narrowest cord radius in this millisecond, m
  std::FILE* trace = options.trace != nullptr ? std::fopen(options.trace, "w") : nullptr;
  std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
  if (!options.hold) {
    tube.set_pressure(breath_target);
  }
  for (long step = 0; step < steps; ++step) {
    if (!options.hold && options.churn && cord_band >= 0 && step % per_tick == 0) {
      tube.set_band_offset(cord_band, design.bands[static_cast<std::size_t>(cord_band)].offset);
    }
    if (options.hold && step % per_tick == 0) {
      hold.tick();
      tube.set_pressure(hold.value(breath));
      if (cord_band >= 0) {
        tube.set_band_offset(cord_band, hold.value(adduction));
        tube.set_band_elasticity(cord_band, hold.value(stretch));
      }
      for (std::size_t r = 0; r < ring_bands.size(); ++r) {
        tube.set_band_offset(ring_bands[r], hold.value(ring_coordinates[r]));
      }
      if (trace != nullptr) {
        std::fprintf(trace, "%.3f %.1f %.2f %.4f %.2f %.2f\n", static_cast<double>(step) / rate, hold.value(breath),
                     cord_band >= 0 ? 1e6 * hold.value(adduction) : 0.0, cord_band >= 0 ? hold.value(stretch) : 0.0,
                     ring_coordinates.empty() ? 0.0 : 1e6 * hold.value(ring_coordinates.front()),
                     2e6 * smallest_gap);
        smallest_gap = 1.0;
      }
      if (static_cast<double>(step) / rate > 0.3) {
        add(breath_spread, hold.value(breath));
        if (cord_band >= 0) {
          add(adduction_spread, hold.value(adduction));
          add(stretch_spread, hold.value(stretch));
        }
        if (!ring_coordinates.empty()) {
          add(ring_spread, hold.value(ring_coordinates.front()));
        }
      }
    }
    if (trace != nullptr) {
      smallest_gap = std::fmin(smallest_gap, tube.cord_radius());
    }
    if (!tube.step()) {
      std::fprintf(stderr, "fault at %.4f s\n", step / rate);
      return 1;
    }
    if (static_cast<double>(step) / rate >= 0.5) {
      glottal.push_back(static_cast<float>(tube.glottal_flow() * 1e6));
    }
    double outlet = design.radiating ? tube.radiated_flow() * 1.0e3 : tube.pressure(tube.count() - 1);
    double value = (outlet - previous_outlet) * rate / 48000.0 / 10.0;
    previous_outlet = outlet;
    for (int k = 0; k < per_output; ++k) {
      double share = static_cast<double>(k) / per_output;
      samples.push_back(static_cast<float>(last + share * (value - last)));
    }
    last = value;
  }
  if (trace != nullptr) {
    std::fclose(trace);
  }
  double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
  write_float_wave(values[3], samples, 48000);
  std::fprintf(stderr, "%s: %.2f s in %.2f s (%.1fx real time), contacts %ld\n", values[3], seconds, elapsed,
               seconds / elapsed, tube.contacts());
  report_source(glottal, rate);
  if (options.hold) {
    std::fprintf(stderr,
                 "held mean (sd): breath %.0f (%.1f) Pa, adduction %.1f (%.2f) um, stretch %.3f (%.4f), "
                 "first ring %.1f (%.1f) um\n",
                 mean(breath_spread), deviation(breath_spread), 1e6 * mean(adduction_spread),
                 1e6 * deviation(adduction_spread), mean(stretch_spread), deviation(stretch_spread),
                 1e6 * mean(ring_spread), 1e6 * deviation(ring_spread));
  }
  return 0;
}
