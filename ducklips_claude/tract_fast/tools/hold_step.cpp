// The step response of one held coordinate: from rest (0) to 1 at a rate.
// usage: hold_step RATE [noise] — prints ms and value until 200 ms.

#include "tract/somatic_hold.hpp"

#include <cstdio>
#include <cstdlib>

int main(int count, char** values) {
  double rate = count > 1 ? std::strtod(values[1], nullptr) : 160.0;
  tract::HoldSettings settings;
  settings.noise = count > 2 ? std::strtod(values[2], nullptr) : 0.0;
  settings.common = 0.0;
  settings.fastest = 1000.0;
  tract::SomaticHold hold;
  hold.reset(7U, settings);
  int index = hold.add(0.0, 1.0, -3.0, 3.0, 0);
  hold.set_target(index, 0.0, rate);
  hold.settle();
  hold.set_target(index, 1.0, rate);
  for (int ms = 1; ms <= 200; ++ms) {
    hold.tick();
    if (ms % 5 == 0) {
      std::printf("%d %.3f\n", ms, hold.value(index));
    }
  }
  return 0;
}
