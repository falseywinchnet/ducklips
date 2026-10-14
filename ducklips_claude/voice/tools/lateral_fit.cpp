// Fits the tongue body and lips of /l/ around a fixed tip contact so the
// tract's resonances match a dark American /l/ (bdl: F1 about 380, F2 about
// 950, F3 about 2850 Hz, with a deep valley between F2 and F3). The tract is
// the body's own geometry with the jaw raised as the planner raises it for
// oral consonants; resonances come from a lossy plane-wave chain matrix
// (closed form per frequency, no iteration).
//
// usage: lateral_fit [F1 F2 F3]

#include "voice/geometry.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

using Complex = std::complex<double>;

constexpr double kPi = 3.14159265358979323846;

// Volume-velocity transfer in dB at one frequency (as tools/tract_codebook).
double response_db(const std::array<double, voice::kOralSections>& areas, double hz) {
  double omega = 2.0 * kPi * hz;
  double attenuation = 1.2e-5 * std::sqrt(hz) + 0.25 / hz;
  Complex wavenumber(omega / voice::kSoundSpeed, -attenuation);
  double lip = std::max(areas[voice::kOralSections - 1], 0.02);
  double radius = std::sqrt(lip / kPi);
  double ka = omega / voice::kSoundSpeed * radius;
  Complex pressure = voice::kRhoC / lip * Complex(0.5 * ka * ka, 8.0 * ka / (3.0 * kPi));
  Complex flow(1.0, 0.0);
  for (int index = voice::kOralSections - 1; index >= 0; --index) {
    double impedance = voice::kRhoC / std::max(areas[static_cast<std::size_t>(index)], 0.005);
    Complex angle = wavenumber * voice::kSectionLength;
    Complex c = std::cos(angle);
    Complex s = std::sin(angle);
    Complex next_pressure = c * pressure + Complex(0.0, 1.0) * impedance * s * flow;
    Complex next_flow = Complex(0.0, 1.0) * s / impedance * pressure + c * flow;
    pressure = next_pressure;
    flow = next_flow;
  }
  return -20.0 * std::log10(std::abs(flow) + 1e-30);
}

int resonances(const std::array<double, voice::kOralSections>& areas, std::array<double, 3>& found) {
  constexpr double kStep = 10.0;
  double previous = response_db(areas, 150.0);
  double here = response_db(areas, 150.0 + kStep);
  int count = 0;
  for (double hz = 150.0 + 2.0 * kStep; hz < 4500.0 && count < 3; hz += kStep) {
    double next = response_db(areas, hz);
    if (here > previous && here >= next) {
      found[static_cast<std::size_t>(count)] = hz - kStep;
      ++count;
    }
    previous = here;
    here = next;
  }
  return count;
}

}  // namespace

int main(int count, char** values) {
  double target[3] = {380.0, 950.0, 2850.0};
  if (count >= 4) {
    for (int index = 0; index < 3; ++index) {
      target[index] = std::strtod(values[index + 1], nullptr);
    }
  }
  double best = 1.0e9;
  voice::TractShape best_shape;
  std::array<double, 3> best_found{};
  std::array<double, voice::kOralSections> areas{};
  const double tip_areas[4] = {0.12, 0.18, 0.25, 0.35};
  for (double tip_area : tip_areas) {
    for (double place = 3.0; place <= 13.01; place += 0.25) {
      for (int a = 0; a < 14; ++a) {
        double tongue_area = 0.25 * std::pow(3.5 / 0.25, a / 13.0);
        for (int l = 0; l < 8; ++l) {
          double lip_area = 1.0 * std::pow(6.0 / 1.0, l / 7.0);
          voice::TractShape shape;
          shape.tongue_place = place;
          shape.tongue_area = tongue_area;
          shape.tip_place = 16.7;
          shape.tip_root = std::sqrt(tip_area);
          shape.lip_root = std::sqrt(lip_area);
          shape.jaw_closure = 1.0;
          voice::tract_areas(shape, areas);
          std::array<double, 3> found{};
          if (resonances(areas, found) < 3) {
            continue;
          }
          double cost = 0.0;
          for (int index = 0; index < 3; ++index) {
            double error = std::log(found[static_cast<std::size_t>(index)] / target[index]);
            cost += error * error;
          }
          if (cost < best) {
            best = cost;
            best_shape = shape;
            best_found = found;
          }
        }
      }
    }
  }
  std::printf("L tongue %.2f area %.3f lips %.2f tip area %.2f -> %.0f %.0f %.0f (target %.0f %.0f %.0f) cost %.4f\n",
              best_shape.tongue_place, best_shape.tongue_area, best_shape.lip_root * best_shape.lip_root,
              best_shape.tip_root * best_shape.tip_root, best_found[0], best_found[1], best_found[2], target[0],
              target[1], target[2], best);
  return 0;
}
