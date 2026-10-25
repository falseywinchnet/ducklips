// Print the cord's air and walls sample by sample over a window.
// usage: sheet_probe DESIGN.txt FROM_SECONDS TO_SECONDS [rate] [every]
// Columns: t, glottal flow (cm^3/s), then for the junction below, each cord
// node and the junction above: pressure (Pa) and radius (um).

#include "tract/sheet_tube.hpp"

#include <cstdio>
#include <cstdlib>

int main(int count, char** values) {
  if (count < 4) {
    std::fprintf(stderr, "usage: sheet_probe DESIGN.txt FROM TO [rate] [every]\n");
    return 2;
  }
  tract::SheetDesign design;
  if (!tract::read_sheet_design(values[1], design)) {
    return 1;
  }
  double from = std::strtod(values[2], nullptr);
  double to = std::strtod(values[3], nullptr);
  double rate = count > 4 ? std::strtod(values[4], nullptr) : 48000.0;
  long every = count > 5 ? std::strtol(values[5], nullptr, 10) : 1;
  tract::SheetTube tube(design, rate);
  int first = -1;
  int last = -1;
  for (int node = 0; node < tube.count(); ++node) {
    if (design.bands[static_cast<std::size_t>(tube.band_of(node))].role == "cord") {
      if (first < 0) {
        first = node;
      }
      last = node;
    }
  }
  long steps = static_cast<long>(to * rate);
  for (long step = 0; step < steps; ++step) {
    if (!tube.step()) {
      std::fprintf(stderr, "fault at %.5f\n", step / rate);
      return 1;
    }
    double t = static_cast<double>(step) / rate;
    if (t < from || step % every != 0) {
      continue;
    }
    std::printf("%.5f %8.2f", t, 1e6 * tube.glottal_flow());
    for (int node = first - 1; node <= last + 1; ++node) {
      std::printf("  %8.1f %7.1f", tube.pressure(node), 1e6 * tube.radius(node));
    }
    std::printf("\n");
  }
  return 0;
}
