// Transfer function of a sheet tube's tract: with no breath and the cord at
// rest, a flow impulse enters the first node above the cord; the flow
// radiated at the lips is recorded and its spectral peaks printed.
// usage: sheet_transfer DESIGN.txt [rate]

#include "tract/sheet_tube.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

double level_db(const std::vector<double>& response, double hz, double rate) {
  std::complex<double> sum(0.0, 0.0);
  double omega = 2.0 * 3.141592653589793 * hz / rate;
  for (std::size_t n = 0; n < response.size(); ++n) {
    sum += response[n] * std::polar(1.0, -omega * static_cast<double>(n));
  }
  return 20.0 * std::log10(std::abs(sum) + 1e-30);
}

}  // namespace

int main(int count, char** values) {
  if (count < 2) {
    std::fprintf(stderr, "usage: sheet_transfer DESIGN.txt [rate]\n");
    return 2;
  }
  double rate = count > 2 ? std::strtod(values[2], nullptr) : 48000.0;
  tract::SheetDesign design;
  if (!tract::read_sheet_design(values[1], design)) {
    std::fprintf(stderr, "cannot read %s\n", values[1]);
    return 1;
  }
  design.pressure = 0.0;
  design.wander = 0.0;
  design.variation = 0.0;
  tract::SheetTube tube(design, rate);
  int cord_band = -1;
  for (int band = 0; band < tube.band_count(); ++band) {
    if (design.bands[static_cast<std::size_t>(band)].role == "cord") {
      cord_band = band;
    }
  }
  int node = cord_band >= 0 ? tube.band_first(cord_band) + tube.band_points(cord_band) : 0;
  std::vector<double> response;
  double previous = 0.0;
  long steps = static_cast<long>(0.25 * rate);
  for (long step = 0; step < steps; ++step) {
    if (step == 0) {
      tube.inject_flow(node, 1.0e-6 * rate);
    }
    tube.step();
    double flow = tube.radiated_flow();
    response.push_back((flow - previous) * rate);
    previous = flow;
  }
  std::vector<double> grid;
  for (double hz = 100.0; hz <= 6000.0; hz += 10.0) {
    grid.push_back(level_db(response, hz, rate));
  }
  double top = -1e9;
  for (double level : grid) {
    top = std::fmax(top, level);
  }
  std::printf("source node %d of %d\n", node, tube.count());
  for (std::size_t i = 1; i + 1 < grid.size(); ++i) {
    if (grid[i] > grid[i - 1] && grid[i] >= grid[i + 1] && grid[i] > top - 40.0) {
      std::size_t low = i;
      while (low > 0 && grid[low] > grid[i] - 3.0) {
        --low;
      }
      std::size_t high = i;
      while (high + 1 < grid.size() && grid[high] > grid[i] - 3.0) {
        ++high;
      }
      std::printf("peak %5.0f Hz  %+6.1f dB  bandwidth %4.0f Hz\n", 100.0 + 10.0 * i, grid[i] - top,
                  10.0 * static_cast<double>(high - low));
    }
  }
  return 0;
}
