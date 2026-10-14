// A codebook of mouth shapes and their sound: for a grid of articulator
// settings, the tract areas (the body's own geometry) and the transfer from
// glottal flow to radiated pressure on a 1/24-octave grid (100 Hz - 7.6 kHz).
// The transfer is a lossy plane-wave chain matrix with piston radiation:
// closed form per frequency, no iteration.
//
// usage: tract_codebook OUT_PREFIX
// Writes OUT_PREFIX_params.npy (configs x 5: tongue place, tongue area,
// lip area, protrusion, jaw closure) and OUT_PREFIX_response.npy
// (configs x frequencies, dB).

#include "time_frequency.hpp"
#include "voice/geometry.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

namespace {

using Complex = std::complex<double>;

constexpr double kPi = 3.14159265358979323846;

std::vector<double> grid_frequencies() {
  std::vector<double> hz;
  for (int index = 0;; ++index) {
    double value = 100.0 * std::pow(2.0, index / 24.0);
    if (value > 7600.0) {
      break;
    }
    hz.push_back(value);
  }
  return hz;
}

// Radiated pressure per unit glottal flow, in dB, for one area function.
double response_db(const std::array<double, voice::kOralSections>& areas, double hz) {
  double omega = 2.0 * kPi * hz;
  double rho_c = voice::kRhoC;
  // Losses: viscous/thermal rising with sqrt(f), walls falling with f.
  double attenuation = 1.2e-5 * std::sqrt(hz) + 0.25 / hz;
  Complex wavenumber(omega / voice::kSoundSpeed, -attenuation);
  double lip = std::max(areas[voice::kOralSections - 1], 0.02);
  double radius = std::sqrt(lip / kPi);
  double ka = omega / voice::kSoundSpeed * radius;
  Complex radiation = rho_c / lip * Complex(0.5 * ka * ka, 8.0 * ka / (3.0 * kPi));
  Complex pressure = radiation;  // per unit lip flow
  Complex flow(1.0, 0.0);
  for (int index = voice::kOralSections - 1; index >= 0; --index) {
    double area = std::max(areas[static_cast<std::size_t>(index)], 0.005);
    double impedance = rho_c / area;
    Complex angle = wavenumber * voice::kSectionLength;
    Complex c = std::cos(angle);
    Complex s = std::sin(angle);
    Complex next_pressure = c * pressure + Complex(0.0, 1.0) * impedance * s * flow;
    Complex next_flow = Complex(0.0, 1.0) * s / impedance * pressure + c * flow;
    pressure = next_pressure;
    flow = next_flow;
  }
  // Radiated far-field pressure ~ omega * lip flow; glottal flow = |flow|.
  return 20.0 * std::log10(omega / std::abs(flow) + 1e-30);
}

}  // namespace

int main(int count, char** values) {
  if (count < 2) {
    std::fprintf(stderr, "usage: tract_codebook OUT_PREFIX\n");
    return 2;
  }
  std::vector<double> hz = grid_frequencies();
  std::vector<float> parameters;
  std::vector<float> responses;
  std::array<double, voice::kOralSections> areas{};
  int configs = 0;
  for (double place = 3.0; place <= 15.01; place += 0.5) {
    for (int a = 0; a < 12; ++a) {
      double tongue_area = 0.15 * std::pow(4.5 / 0.15, a / 11.0);
      for (int l = 0; l < 12; ++l) {
        double lip_area = 0.1 * std::pow(12.0 / 0.1, l / 11.0);
        for (double protrusion : {0.0, 0.5, 1.0}) {
          for (double jaw : {0.0, 0.5, 1.0}) {
            voice::TractShape shape;
            shape.tongue_place = place;
            shape.tongue_area = tongue_area;
            shape.lip_root = std::sqrt(lip_area);
            shape.protrusion = protrusion;
            shape.jaw_closure = jaw;
            voice::tract_areas(shape, areas);
            parameters.push_back(static_cast<float>(place));
            parameters.push_back(static_cast<float>(tongue_area));
            parameters.push_back(static_cast<float>(lip_area));
            parameters.push_back(static_cast<float>(protrusion));
            parameters.push_back(static_cast<float>(jaw));
            for (double frequency : hz) {
              responses.push_back(static_cast<float>(response_db(areas, frequency)));
            }
            ++configs;
          }
        }
      }
    }
  }
  std::string prefix = values[1];
  bool good = tf::write_npy(prefix + "_params.npy", parameters, configs, 5) &&
              tf::write_npy(prefix + "_response.npy", responses, configs, static_cast<std::int32_t>(hz.size()));
  std::fprintf(stderr, "%d configurations x %zu frequencies\n", configs, hz.size());
  return good ? 0 : 1;
}
