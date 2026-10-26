// The reversible jet's source (kernels_jet.jet_source, VoiceBody's
// reversible_jet) for the native kernels (peel_core.cpp, forward_core.cpp),
// to the bit: the numba expressions operand for operand, compiled with
// -ffp-contract=off, the same libm exp and the square root instruction.
#ifndef JET_CORE_HPP
#define JET_CORE_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace jet_core {

using Index = std::int64_t;

constexpr double kPi = 3.141592653589793;      // numpy.pi
constexpr double kVelocityFinite = 1e100;      // kernels_jet.VELOCITY_FINITE

// The jet's fixed inputs (borrowed).
struct JetParameters {
  double bound;          // JETBOUND's speed (0: none)
  double closing;        // REVJETCLOSE's fade (0: none)
  const double* noise;   // noise_links x noise_steps (kernels_jet.noise_table), C order
  Index noise_links;
  Index noise_steps;
  Index seed;
};

// A link's filter and coefficient cache: four values each, wherever they
// are stored (a state row's columns, or a layer field each).
struct JetState {
  double* filter[4];
  double* cache[4];     // corner, its two gains, the root of the passed power
};

inline std::uint64_t split_mix(std::uint64_t key) {
  key = key + 0x9E3779B97F4A7C15ULL;
  key = (key ^ (key >> 30)) * 0xBF58476D1CE4E5B9ULL;
  key = (key ^ (key >> 27)) * 0x94D049BB133111EBULL;
  std::uint64_t mixed = key ^ (key >> 31);
  return mixed;
}

// kernels_jet.hashed_gauss.
inline double hashed_gauss(Index seed, Index link, Index step) {
  double total = 0.0;
  for (std::uint64_t draw = 0; draw < 4; ++draw) {
    std::uint64_t key = static_cast<std::uint64_t>(seed) ^
                        ((static_cast<std::uint64_t>(step) * 4096ULL + static_cast<std::uint64_t>(link)) * 4ULL + draw);
    double uniform = static_cast<double>(split_mix(key) >> 11) * (1.0 / 9007199254740992.0);
    total += uniform;
  }
  double gauss = (total - 2.0) * std::sqrt(3.0);
  return gauss;
}

// The source from the shaped noise (filter3) once the filter has run:
// kernels_jet.jet_source's last part (0 before the narrowest node; the
// silent value with the gate shut).
inline double jet_strength(const JetParameters& jet, Index link, Index narrowest, double area, double velocity,
                           double gain, double filter3, double root_passed) {
  if (narrowest < 0 || link < narrowest) {
    return 0.0;
  }
  if (gain == 0.0 && velocity <= kVelocityFinite) {
    double silent = (gain * 0.0) * filter3 / root_passed;
    return silent;
  }
  double diameter = std::sqrt(4.0 * area / kPi);
  double reynolds = velocity * diameter * 1.14e-3 / 1.86e-4;
  double onset = 1800.0 / std::max(reynolds, 1.0);
  double turbulent = std::max(1.0 - onset * onset, 0.0);
  double fade = std::min(std::max((0.6 - area) / 0.3, 0.0), 1.0);
  if (jet.closing > 0.0) {
    double closed = 0.25 * jet.closing;
    fade *= std::min(std::max((area - closed) / (jet.closing - closed), 0.0), 1.0);
  }
  double dynamic = 0.5 * 1.14e-3 * velocity * velocity * turbulent * fade * fade * (velocity / 2000.0);
  double distance = static_cast<double>(link - narrowest) * (35000.0 / 96000.0);
  double strength = gain * dynamic * std::sqrt(area) * std::exp(-distance / 1.0);
  double source = strength * filter3 / root_passed;
  return source;
}


// kernels_jet.jet_source: the link's source for `step` (its filter advanced).
inline double jet_source(const JetParameters& jet, Index link, double previous_flow, Index narrowest, double area,
                         const JetState& held, double gain, Index step) {
  double velocity = std::abs(previous_flow) / area;
  if (jet.bound > 0.0) {
    velocity = std::min(velocity, jet.bound);
  }
  double slit = std::max(area / 1.2, 0.01);
  double corner = std::min(std::max(0.2 * velocity / slit, 1500.0), 12000.0);
  double corner_gain = 0.0;
  double floor_gain = 0.0;
  double root_passed = 0.0;
  if (corner == *held.cache[0]) {
    corner_gain = *held.cache[1];
    floor_gain = *held.cache[2];
    root_passed = *held.cache[3];
  } else {
    corner_gain = 1.0 - std::exp(-2.0 * kPi * corner / 96000.0);
    floor_gain = 1.0 - std::exp(-2.0 * kPi * corner / (3.0 * 96000.0));
    double pole = 1.0 - corner_gain;
    double gain_squared = corner_gain * corner_gain;
    double pole_term = 1.0 - pole * pole;
    double passed = gain_squared * gain_squared * (1.0 + pole * pole) / (pole_term * pole_term * pole_term);
    root_passed = std::sqrt(passed);
    *held.cache[0] = corner;
    *held.cache[1] = corner_gain;
    *held.cache[2] = floor_gain;
    *held.cache[3] = root_passed;
  }
  double white = 0.0;
  if (link < jet.noise_links && step < jet.noise_steps) {
    white = jet.noise[link * jet.noise_steps + step];
  } else {
    white = hashed_gauss(jet.seed, link, step);
  }
  double* filter0 = held.filter[0];
  double* filter1 = held.filter[1];
  double* filter2 = held.filter[2];
  double* filter3 = held.filter[3];
  *filter0 += floor_gain * (white - *filter0);
  double once = white - *filter0;
  *filter1 += floor_gain * (once - *filter1);
  *filter2 += corner_gain * ((once - *filter1) - *filter2);
  *filter3 += corner_gain * (*filter2 - *filter3);
  double source = jet_strength(jet, link, narrowest, area, velocity, gain, *filter3, root_passed);
  return source;
}

// reverse_core._narrowest of one row's areas (stride `stride` doubles):
// the narrowest node from node 10 under 0.6 cm^2 (-1: none) and its area
// floored at 1e-4.
inline void narrowest_of(const double* areas, Index stride, Index nodes, Index& node_out, double& area_out) {
  Index narrowest = -1;
  double least = 0.6;
  for (Index node = 10; node < nodes; ++node) {
    double area = areas[node * stride];
    if (area < least) {
      least = area;
      narrowest = node;
    }
  }
  node_out = narrowest;
  area_out = std::max(least, 1e-4);
}

}  // namespace jet_core

#endif
