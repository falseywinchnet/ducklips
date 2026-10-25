#include "tract/variation.hpp"

#include <cmath>

namespace tract {

Variation::Variation(int channels, std::uint32_t seed) : random_(seed != 0U ? seed : 1U) {
  states_.resize(static_cast<std::size_t>(channels));
  for (State& state : states_) {
    state.phase = uniform();
    state.frequency = 0.35 + 0.9 * uniform();
    state.wait = 0.2 + uniform();
  }
}

double Variation::uniform() {
  std::uint32_t x = random_;
  x ^= x << 13U;
  x ^= x >> 17U;
  x ^= x << 5U;
  random_ = x;
  return static_cast<double>(x) / 4294967296.0;
}

void Variation::advance(double dt) {
  double noise_gain = 1.0 - std::exp(-dt / 0.08);
  double saw_gain = 1.0 - std::exp(-dt / 0.045);
  double value_gain = 1.0 - std::exp(-dt / 0.025);
  for (State& state : states_) {
    // Drift: a new random target now and then, approached smoothly.
    state.hold -= dt;
    if (state.hold <= 0.0) {
      state.target = 2.0 * uniform() - 1.0;
      state.hold = 0.08 + 0.32 * uniform();
    }
    state.noise += (state.target - state.noise) * noise_gain;
    // Slow build-up and recovery.
    state.phase = std::fmod(state.phase + dt * state.frequency, 1.0);
    state.saw += (2.0 * state.phase - 1.0 - state.saw) * saw_gain;
    // An occasional excursion that rises and subsides (peak 1 at x = 1/4).
    state.wait -= dt;
    state.age += dt;
    if (state.wait <= 0.0) {
      state.age = 0.0;
      state.duration = 0.16 + 0.4 * uniform();
      state.sign = uniform() < 0.5 ? -1.0 : 1.0;
      state.wait = state.duration + 0.15 + 1.2 * uniform();
    }
    double x = state.age / state.duration;
    state.bubble = x < 1.0 ? state.sign * (256.0 / 27.0) * x * std::pow(1.0 - x, 3.0) : 0.0;
    double mixed = 0.4 * state.noise + 0.25 * state.saw + 0.35 * state.bubble;
    state.value += (mixed - state.value) * value_gain;
  }
}

}  // namespace tract
