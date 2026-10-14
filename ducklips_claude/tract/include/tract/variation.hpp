#pragma once

// Slow, history-bearing variation of physical causes (muscle effort, tissue
// elasticity, breath), never of the output waveform. Ported from Astra's sheet
// tube (variation.mjs at 42c7337): each channel mixes 40% smoothed random
// drift (new target every 80-400 ms, smoothed over 80 ms), 25% smoothed
// sawtooth (0.35-1.25 Hz) and 35% occasional "bubbles" (one finite excursion
// that rises and subsides over 160-560 ms), then smooths the mix over 25 ms.
// Values stay within [-1, 1]; callers scale them by a bounded depth.

#include <cstdint>
#include <vector>

namespace tract {

class Variation {
 public:
  Variation(int channels, std::uint32_t seed);

  // Advance every channel by dt seconds (control rate).
  void advance(double dt);
  double value(int channel) const { return states_[static_cast<std::size_t>(channel)].value; }
  int channels() const { return static_cast<int>(states_.size()); }

 private:
  struct State {
    double phase = 0.0;
    double frequency = 0.35;
    double noise = 0.0;
    double target = 0.0;
    double hold = 0.0;
    double saw = 0.0;
    double bubble = 0.0;
    double age = 2.0;
    double duration = 0.3;
    double sign = 1.0;
    double wait = 0.2;
    double value = 0.0;
  };

  double uniform();

  std::uint32_t random_;
  std::vector<State> states_;
};

}  // namespace tract
