// Fits /l/ to the human band profile instead of three formants: the target
// is bdl's mean /l/ minus his mean /ah/ (1/6-octave bands, middle slice), so
// the voice source and radiation cancel; the model is the same difference
// between two closed-form transfers through the body's own geometry
// (lossy plane-wave chain matrices, as tools/lateral_fit). The /l/ may carry
// a blind side branch opening just behind the tip contact (the pocket over
// the blade that the lateral channels bypass), whose quarter-wave zero
// deepens the valley above F2. Prints the best fit without and with it.
//
// usage: lateral_profile PROFILE.txt
//   PROFILE.txt: lines "l ref v1 v2 ...", "ah ref ...", "hz c1 c2 ..."
//   (tools/phone_gap-style band means, as written by out/profile.py)

#include "voice/geometry.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Complex = std::complex<double>;
using Areas = std::array<double, voice::kOralSections>;

constexpr double kPi = 3.14159265358979323846;
constexpr double kLowestHz = 250.0;
constexpr double kHighestHz = 5000.0;

struct Branch {
  double length = 0.0;  // cm; 0 means no branch
  double area = 0.0;    // cm^2
  int section = 0;      // oral section the branch opens into
};

// Volume-velocity transfer (dB) at one frequency, with an optional closed
// side branch as a shunt admittance at one section.
double response_db(const Areas& areas, const Branch& branch, double hz) {
  double omega = 2.0 * kPi * hz;
  double attenuation = 1.2e-5 * std::sqrt(hz) + 0.25 / hz;
  Complex wavenumber(omega / voice::kSoundSpeed, -attenuation);
  double lip = std::max(areas[voice::kOralSections - 1], 0.02);
  double radius = std::sqrt(lip / kPi);
  double ka = omega / voice::kSoundSpeed * radius;
  Complex pressure = voice::kRhoC / lip * Complex(0.5 * ka * ka, 8.0 * ka / (3.0 * kPi));
  Complex flow(1.0, 0.0);
  Complex branch_admittance(0.0, 0.0);
  if (branch.length > 0.0) {
    // Closed tube: input admittance j (A / rho c) tan(k l).
    branch_admittance = Complex(0.0, 1.0) * (branch.area / voice::kRhoC) * std::tan(wavenumber * branch.length);
  }
  for (int index = voice::kOralSections - 1; index >= 0; --index) {
    if (branch.length > 0.0 && index == branch.section) {
      flow += branch_admittance * pressure;
    }
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

// Band means of the transfer: three frequencies per band, power-averaged.
void band_levels(const Areas& areas, const Branch& branch, const std::vector<double>& centres,
                 std::vector<double>& levels) {
  for (std::size_t band = 0; band < centres.size(); ++band) {
    double power = 0.0;
    for (int offset = -1; offset <= 1; ++offset) {
      double hz = centres[band] * std::pow(2.0, offset / 18.0);
      power += std::pow(10.0, response_db(areas, branch, hz) / 10.0);
    }
    levels[band] = 10.0 * std::log10(power / 3.0);
  }
}

double profile_cost(const std::vector<double>& model, const std::vector<double>& target,
                    const std::vector<bool>& used) {
  double mean = 0.0;
  int count = 0;
  for (std::size_t band = 0; band < model.size(); ++band) {
    if (used[band]) {
      mean += model[band] - target[band];
      ++count;
    }
  }
  mean /= count;
  double cost = 0.0;
  for (std::size_t band = 0; band < model.size(); ++band) {
    if (used[band]) {
      double error = model[band] - target[band] - mean;
      cost += error * error;
    }
  }
  return std::sqrt(cost / count);
}

bool read_profile(const char* path, std::vector<double>& lateral, std::vector<double>& vowel,
                  std::vector<double>& centres) {
  std::ifstream input(path);
  std::string line;
  while (std::getline(input, line)) {
    std::istringstream fields(line);
    std::string phone;
    std::string side;
    fields >> phone;
    std::vector<double>* target = nullptr;
    if (phone == "hz") {
      target = &centres;
    } else {
      fields >> side;
      if (side != "ref") {
        continue;
      }
      target = phone == "l" ? &lateral : phone == "ah" ? &vowel : nullptr;
    }
    if (target == nullptr) {
      continue;
    }
    double value = 0.0;
    while (fields >> value) {
      target->push_back(value);
    }
  }
  return !lateral.empty() && lateral.size() == vowel.size() && centres.size() >= lateral.size();
}

struct Fit {
  double cost = 1.0e9;
  voice::TractShape shape;
  Branch branch;
};

}  // namespace

int main(int count, char** values) {
  if (count < 2) {
    std::fprintf(stderr, "usage: lateral_profile PROFILE.txt\n");
    return 2;
  }
  std::vector<double> lateral;
  std::vector<double> vowel;
  std::vector<double> centres;
  if (!read_profile(values[1], lateral, vowel, centres)) {
    std::fprintf(stderr, "cannot read %s\n", values[1]);
    return 1;
  }
  std::size_t bands = lateral.size();
  centres.resize(bands);
  std::vector<double> target(bands, 0.0);
  std::vector<bool> used(bands, false);
  for (std::size_t band = 0; band < bands; ++band) {
    target[band] = lateral[band] - vowel[band];
    used[band] = centres[band] >= kLowestHz && centres[band] <= kHighestHz;
  }
  // The /ah/ reference: the planner's AH posture (place 8.5, area 0.722,
  // lips 5.13), jaw open.
  voice::TractShape reference_shape;
  reference_shape.tongue_place = 8.5;
  reference_shape.tongue_area = 0.722;
  reference_shape.lip_root = std::sqrt(5.13);
  Areas reference_areas{};
  voice::tract_areas(reference_shape, reference_areas);
  std::vector<double> reference(bands, 0.0);
  Branch none;
  band_levels(reference_areas, none, centres, reference);

  const double tip_places[5] = {14.6, 15.2, 15.6, 16.2, 16.7};
  const double tip_areas[4] = {0.15, 0.3, 0.45, 0.7};
  const double branch_lengths[7] = {0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0};
  const double branch_areas[3] = {0.2, 0.5, 1.0};
  std::array<Fit, 2> best{};
  std::vector<double> model(bands, 0.0);
  Areas areas{};
  for (double tip_place : tip_places) {
    for (double tip_area : tip_areas) {
      for (double place = 4.0; place <= 12.01; place += 1.0) {
        for (int a = 0; a < 7; ++a) {
          double tongue_area = 0.18 * std::pow(3.0 / 0.18, a / 6.0);
          for (int l = 0; l < 4; ++l) {
            double lip_area = 1.0 * std::pow(5.0, l / 3.0);
            voice::TractShape shape;
            shape.tongue_place = place;
            shape.tongue_area = tongue_area;
            shape.tip_place = tip_place;
            shape.tip_root = std::sqrt(tip_area);
            shape.lip_root = std::sqrt(lip_area);
            shape.jaw_closure = 1.0;
            voice::tract_areas(shape, areas);
            int contact = static_cast<int>(tip_place / voice::kSectionLength);
            for (double length : branch_lengths) {
              for (double branch_area : branch_areas) {
                if (length == 0.0 && branch_area != branch_areas[0]) {
                  continue;
                }
                Branch branch;
                branch.length = length;
                branch.area = branch_area;
                branch.section = std::max(contact - 1, 0);
                band_levels(areas, branch, centres, model);
                for (std::size_t band = 0; band < bands; ++band) {
                  model[band] -= reference[band];
                }
                double cost = profile_cost(model, target, used);
                Fit& slot = best[length > 0.0 ? 1 : 0];
                if (cost < slot.cost) {
                  slot.cost = cost;
                  slot.shape = shape;
                  slot.branch = branch;
                }
              }
            }
          }
        }
      }
    }
  }
  const char* labels[2] = {"without branch", "with branch"};
  for (int which = 0; which < 2; ++which) {
    const Fit& fit = best[static_cast<std::size_t>(which)];
    std::printf("%-15s rms %.2f dB: tongue %.2f area %.3f tip %.2f tip area %.2f lips %.2f branch %.1f cm x %.1f cm2\n",
                labels[which], fit.cost, fit.shape.tongue_place, fit.shape.tongue_area, fit.shape.tip_place,
                fit.shape.tip_root * fit.shape.tip_root, fit.shape.lip_root * fit.shape.lip_root, fit.branch.length,
                fit.branch.area);
    voice::tract_areas(fit.shape, areas);
    band_levels(areas, fit.branch, centres, model);
    std::printf("  band  target  model\n");
    double mean = 0.0;
    int used_count = 0;
    for (std::size_t band = 0; band < bands; ++band) {
      model[band] -= reference[band];
      if (used[band]) {
        mean += model[band] - target[band];
        ++used_count;
      }
    }
    mean /= used_count;
    for (std::size_t band = 0; band < bands; band += 2) {
      if (used[band]) {
        std::printf("  %5.0f %6.1f %6.1f\n", centres[band], target[band], model[band] - mean);
      }
    }
  }
  return 0;
}
