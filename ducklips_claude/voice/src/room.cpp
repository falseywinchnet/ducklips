#include "voice/room.hpp"

#include <array>
#include <cmath>

namespace voice {

namespace {

constexpr std::int32_t kLines = 4;
// Mutually prime delays (ms) of a small room's early paths.
constexpr std::array<double, kLines> kDelayMs = {29.7, 37.1, 41.1, 43.7};

double next_uniform(std::uint64_t& state) {
  state ^= state << 13U;
  state ^= state >> 7U;
  state ^= state << 17U;
  return static_cast<double>(state >> 11U) * (1.0 / 9007199254740992.0) - 0.5;
}

// One biquad section (direct form I), RBJ cookbook coefficients.
struct Biquad {
  double b0 = 1.0;
  double b1 = 0.0;
  double b2 = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
};

Biquad low_shelf(double rate, double hz, double gain_db) {
  double a = std::pow(10.0, gain_db / 40.0);
  double omega = 2.0 * 3.14159265358979323846 * hz / rate;
  double cosine = std::cos(omega);
  double alpha = std::sin(omega) / 2.0 * std::sqrt(2.0);  // shelf slope 1
  double root = 2.0 * std::sqrt(a) * alpha;
  double a0 = (a + 1.0) + (a - 1.0) * cosine + root;
  Biquad filter;
  filter.b0 = a * ((a + 1.0) - (a - 1.0) * cosine + root) / a0;
  filter.b1 = 2.0 * a * ((a - 1.0) - (a + 1.0) * cosine) / a0;
  filter.b2 = a * ((a + 1.0) - (a - 1.0) * cosine - root) / a0;
  filter.a1 = -2.0 * ((a - 1.0) + (a + 1.0) * cosine) / a0;
  filter.a2 = ((a + 1.0) + (a - 1.0) * cosine - root) / a0;
  return filter;
}

Biquad peak(double rate, double hz, double gain_db) {
  double a = std::pow(10.0, gain_db / 40.0);
  double omega = 2.0 * 3.14159265358979323846 * hz / rate;
  double alpha = std::sin(omega) / (2.0 * 0.7);  // Q 0.7: about two octaves between half-gain points
  double a0 = 1.0 + alpha / a;
  Biquad filter;
  filter.b0 = (1.0 + alpha * a) / a0;
  filter.b1 = -2.0 * std::cos(omega) / a0;
  filter.b2 = (1.0 - alpha * a) / a0;
  filter.a1 = -2.0 * std::cos(omega) / a0;
  filter.a2 = (1.0 - alpha / a) / a0;
  return filter;
}

void apply_biquad(std::vector<float>& samples, const Biquad& filter) {
  double x1 = 0.0;
  double x2 = 0.0;
  double y1 = 0.0;
  double y2 = 0.0;
  for (float& sample : samples) {
    double x0 = sample;
    double y0 = filter.b0 * x0 + filter.b1 * x1 + filter.b2 * x2 - filter.a1 * y1 - filter.a2 * y2;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
    sample = static_cast<float>(y0);
  }
}

}  // namespace

void apply_room(std::vector<float>& samples, double rate, const RoomOptions& options) {
  if (samples.empty()) {
    return;
  }
  if (options.proximity_db != 0.0) {
    apply_biquad(samples, low_shelf(rate, options.proximity_hz, options.proximity_db));
  }
  if (options.presence_db != 0.0) {
    apply_biquad(samples, peak(rate, options.presence_hz, options.presence_db));
  }
  std::array<std::vector<double>, kLines> lines;
  std::array<std::size_t, kLines> positions{};
  std::array<double, kLines> gains{};
  std::array<double, kLines> damping{};
  for (std::int32_t line = 0; line < kLines; ++line) {
    std::size_t length = static_cast<std::size_t>(kDelayMs[static_cast<std::size_t>(line)] * 0.001 * rate);
    lines[static_cast<std::size_t>(line)].assign(length, 0.0);
    // Gain giving a 60 dB decay over the reverberation time.
    gains[static_cast<std::size_t>(line)] =
        std::pow(10.0, -3.0 * (static_cast<double>(length) / rate) / options.reverb_seconds);
  }
  double send = std::pow(10.0, options.reverb_level_db / 20.0);
  double peak = 0.0;
  for (float value : samples) {
    peak = std::max(peak, static_cast<double>(std::abs(value)));
  }
  double noise = peak * std::pow(10.0, options.noise_level_db / 20.0);
  std::uint64_t state = 0x2545F4914F6CDD1DULL ^ options.seed;
  double pink = 0.0;
  double absorb = 1.0 - std::exp(-2.0 * 3.14159265358979323846 * 4000.0 / rate);  // walls absorb highs
  for (float& sample : samples) {
    std::array<double, kLines> outs{};
    for (std::int32_t line = 0; line < kLines; ++line) {
      std::size_t index = static_cast<std::size_t>(line);
      outs[index] = lines[index][positions[index]];
    }
    // Hadamard mixing keeps the network lossless before the decay gains.
    std::array<double, kLines> mixed = {
        0.5 * (outs[0] + outs[1] + outs[2] + outs[3]), 0.5 * (outs[0] - outs[1] + outs[2] - outs[3]),
        0.5 * (outs[0] + outs[1] - outs[2] - outs[3]), 0.5 * (outs[0] - outs[1] - outs[2] + outs[3])};
    double dry = sample;
    double tail = 0.25 * (outs[0] + outs[1] + outs[2] + outs[3]);
    for (std::int32_t line = 0; line < kLines; ++line) {
      std::size_t index = static_cast<std::size_t>(line);
      double input = gains[index] * mixed[index] + send * dry;
      damping[index] += absorb * (input - damping[index]);
      lines[index][positions[index]] = damping[index];
      positions[index] = (positions[index] + 1) % lines[index].size();
    }
    pink += 0.05 * (next_uniform(state) - pink);
    sample = static_cast<float>(dry + tail + noise * 6.0 * pink);
  }
}

}  // namespace voice
