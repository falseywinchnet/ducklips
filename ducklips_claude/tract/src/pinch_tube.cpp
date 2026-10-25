#include "tract/pinch_tube.hpp"

#include "tract/units.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

constexpr int kSubsteps = 8;
constexpr double kJetRecovery = 0.2;   // share of the kinetic pressure recovered past the pinch
constexpr double kInertiaDecay = 0.95;  // the air mass's trapezoidal step, as the cord's
constexpr double kInertanceFloor = 0.002;  // cm^2: the air's mass is taken on openings at least this

double half_rho_over(double area) {
  double a = std::max(area, kMinimumArea);
  return 0.5 * kAirDensity / (a * a);
}

}  // namespace

void PinchTube::configure(const PinchSettings& settings, std::uint32_t seed) {
  settings_ = settings;
  random_ = seed != 0U ? seed : 1U;
  rest_gap_[0] = settings_.outer_gap;
  rest_gap_[1] = settings_.inner_gap;
  rest_gap_[2] = 0.01;
  rest_gap_[3] = 0.001;
  rest_gap_[4] = settings_.inner_gap;
  rest_gap_[5] = settings_.outer_gap;
  y_.fill(0.0);
  v_.fill(0.0);
  drift_.fill(0.0);
  eddy_.fill(0.0);
  modes_ = std::clamp(static_cast<int>(std::lround(settings_.string_modes)), 0, kMaxStringModes);
  if (modes_ > 0) {
    // The standing waves at the points' middles, the mode force a uniform
    // pressure makes (the shape summed over the wall), and each mode's loss.
    double point_width = settings_.width / static_cast<double>(kStringPoints);
    for (int mode = 0; mode < modes_; ++mode) {
      std::size_t n = static_cast<std::size_t>(mode);
      double weight = 0.0;
      for (int point = 0; point < kStringPoints; ++point) {
        double x = (static_cast<double>(point) + 0.5) / static_cast<double>(kStringPoints);
        double shape = std::sin(static_cast<double>(mode + 1) * kPi * x);
        shape_[n][static_cast<std::size_t>(point)] = shape;
        weight += shape * point_width * settings_.section_length;
      }
      pressure_weight_[n] = weight;
      mode_loss_[n] = settings_.loss * std::pow(static_cast<double>(mode + 1), settings_.string_loss_growth);
    }
    for (int section = 0; section < kSections; ++section) {
      mode_y_[static_cast<std::size_t>(section)].fill(0.0);
      mode_v_[static_cast<std::size_t>(section)].fill(0.0);
    }
    string_openings();
    return;
  }
  areas();
}

void PinchTube::set_rest(double lower_gap, double upper_gap, double leak) {
  rest_gap_[2] = lower_gap;
  rest_gap_[3] = upper_gap;
  leak_ = std::max(leak, 0.0);
}

void PinchTube::section_factors(int section, double& mass, double& elasticity, double& support) const {
  const double masses[kSections] = {settings_.s1_mass, settings_.s2_mass, settings_.s3_mass,
                                    settings_.s4_mass, settings_.s5_mass, settings_.s6_mass};
  const double elastics[kSections] = {settings_.s1_elasticity, settings_.s2_elasticity, settings_.s3_elasticity,
                                      settings_.s4_elasticity, settings_.s5_elasticity, settings_.s6_elasticity};
  const double supports[kSections] = {settings_.s1_support, settings_.s2_support, settings_.s3_support,
                                      settings_.s4_support, settings_.s5_support, settings_.s6_support};
  std::size_t s = static_cast<std::size_t>(section);
  mass = masses[s];
  elasticity = elastics[s];
  support = supports[s];
}

// Each section's opening (rest gap plus both walls' displacement) and its
// area; the narrowest.
void PinchTube::areas() {
  narrowest_ = 1.0e9;
  narrowest_index_ = 0;
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    opening_[s] = rest_gap_[s] + 2.0 * y_[s];
    area_[s] = settings_.width * std::max(opening_[s], 0.0);
    if (area_[s] < narrowest_) {
      narrowest_ = area_[s];
      narrowest_index_ = section;
    }
  }
}

// The pressure on each section's walls: from the inlet, less the kinetic
// pressure of the speed there and the share of the air mass's pressure and
// the friction met on the way, up to the narrowest section; past it the
// jet has separated and the walls feel the outlet's pressure.
void PinchTube::pressures(double p_in, double p_out, double area_below) {
  double total_inertance = 0.0;
  for (int section = 0; section < kSections; ++section) {
    total_inertance += settings_.section_length / std::max(area_[static_cast<std::size_t>(section)] + leak_, kInertanceFloor);
  }
  double speed_term = flow_ * std::abs(flow_);
  double inertance_so_far = 0.0;
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double a = area_[s] + leak_;
    inertance_so_far += 0.5 * settings_.section_length / std::max(a, kInertanceFloor);
    if (section <= narrowest_index_) {
      double kinetic = speed_term * (half_rho_over(a) - half_rho_over(area_below));
      double inertial = total_inertance > 0.0 ? inertial_pressure_ * inertance_so_far / total_inertance : 0.0;
      pressure_[s] = p_in - kinetic - inertial;
    } else {
      pressure_[s] = p_out;
    }
    inertance_so_far += 0.5 * settings_.section_length / std::max(a, kInertanceFloor);
  }
}

double PinchTube::gauss() {
  double sum = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    random_ ^= random_ << 13U;
    random_ ^= random_ >> 7U;
    random_ ^= random_ << 17U;
    sum += static_cast<double>(random_ >> 11U) / 9007199254740992.0;
  }
  return (sum - 2.0) * std::sqrt(3.0);
}

// The walls over one substep: the air's pressure on each wall, its support
// to rest, its elasticity to the neighbouring walls (stretched by the
// tension, wandering), its loss, and contact where opposite walls meet.
void PinchTube::advance_walls(double h) {
  double wall_area = settings_.width * settings_.section_length;
  std::array<double, kSections> mass{};
  std::array<double, kSections> elastic{};
  std::array<double, kSections> support{};
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double m = 0.0;
    double e = 0.0;
    double k = 0.0;
    section_factors(section, m, e, k);
    mass[s] = settings_.surface_density * wall_area * m;
    elastic[s] = settings_.elasticity * e * tension_ * std::max(1.0 + drift_[s], 0.05);
    support[s] = settings_.support * k * tension_;
  }
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double force = pressure_[s] * wall_area - support[s] * y_[s];
    double stiffness = support[s];
    if (section > 0) {
      double between = 0.5 * (elastic[s] + elastic[s - 1]);
      force += between * (y_[s - 1] - y_[s]);
      stiffness += between;
    }
    if (section + 1 < kSections) {
      double between = 0.5 * (elastic[s] + elastic[s + 1]);
      force += between * (y_[s + 1] - y_[s]);
      stiffness += between;
    }
    double loss = settings_.loss;
    double overlap = std::max(-opening_[s], 0.0) * 0.5;
    if (overlap > 0.0) {
      force += settings_.contact_stiffness * overlap;
      stiffness += settings_.contact_stiffness;
      loss = std::max(loss, settings_.contact_loss);
    }
    force -= 2.0 * loss * std::sqrt(mass[s] * stiffness) * v_[s];
    v_[s] += h * force / mass[s];
  }
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    y_[s] += h * v_[s];
  }
  areas();
}

// The strings' shapes across the slit: each point's displacement and speed
// from the modes, its opening (rest gap, widened toward one end by the
// taper, plus both walls), and each section's area and narrowest.
void PinchTube::string_openings() {
  double point_width = settings_.width / static_cast<double>(kStringPoints);
  narrowest_ = 1.0e9;
  narrowest_index_ = 0;
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double area = 0.0;
    double least = 1.0e9;
    for (int point = 0; point < kStringPoints; ++point) {
      std::size_t j = static_cast<std::size_t>(point);
      double y = 0.0;
      double v = 0.0;
      for (int mode = 0; mode < modes_; ++mode) {
        std::size_t n = static_cast<std::size_t>(mode);
        y += mode_y_[s][n] * shape_[n][j];
        v += mode_v_[s][n] * shape_[n][j];
      }
      point_y_[s][j] = y;
      point_v_[s][j] = v;
      double x = (static_cast<double>(point) + 0.5) / static_cast<double>(kStringPoints);
      double opening = rest_gap_[s] + settings_.string_taper * x + 2.0 * y;
      point_opening_[s][j] = opening;
      area += point_width * std::max(opening, 0.0);
      least = std::min(least, opening);
    }
    opening_[s] = least;
    area_[s] = area;
    if (area < narrowest_) {
      narrowest_ = area;
      narrowest_index_ = section;
    }
  }
}

// The strings over one substep. Each mode, per unit of its antinode
// displacement, has half the wall's mass, support and neighbour elasticity
// (the shape's square averaged over the slit), plus the string's tension
// stiffness growing as n^2. The air's pressure drives the odd modes; where
// opposite walls meet, the contact at each point drives every mode by the
// shape there: the pluck.
void PinchTube::advance_strings(double h) {
  double wall_area = settings_.width * settings_.section_length;
  double point_share = 1.0 / static_cast<double>(kStringPoints);
  std::array<double, kSections> mass{};
  std::array<double, kSections> elastic{};
  std::array<double, kSections> support{};
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double m = 0.0;
    double e = 0.0;
    double k = 0.0;
    section_factors(section, m, e, k);
    mass[s] = settings_.surface_density * wall_area * m;
    elastic[s] = settings_.elasticity * e * tension_ * std::max(1.0 + drift_[s], 0.05);
    support[s] = settings_.support * k * tension_;
  }
  std::array<double, kStringPoints> contact_force{};
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double point_mass = mass[s] * point_share;
    double point_contact = settings_.contact_stiffness * point_share;
    double contact_damping = 2.0 * settings_.contact_loss * std::sqrt(point_mass * point_contact);
    bool touching = false;
    for (int point = 0; point < kStringPoints; ++point) {
      std::size_t j = static_cast<std::size_t>(point);
      double overlap = std::max(-point_opening_[s][j], 0.0) * 0.5;
      contact_force[j] = overlap > 0.0 ? point_contact * overlap - contact_damping * point_v_[s][j] : 0.0;
      touching = touching || overlap > 0.0;
    }
    double modal_mass = 0.5 * mass[s];
    for (int mode = 0; mode < modes_; ++mode) {
      std::size_t n = static_cast<std::size_t>(mode);
      double wavenumber = static_cast<double>(mode + 1) * kPi / settings_.width;
      double string = settings_.string_tension * tension_ * wavenumber * wavenumber * 0.5 * settings_.width;
      double stiffness = 0.5 * support[s] + string;
      double force = pressure_[s] * pressure_weight_[n] - stiffness * mode_y_[s][n];
      if (section > 0) {
        double between = 0.25 * (elastic[s] + elastic[s - 1]);
        force += between * (mode_y_[s - 1][n] - mode_y_[s][n]);
        stiffness += between;
      }
      if (section + 1 < kSections) {
        double between = 0.25 * (elastic[s] + elastic[s + 1]);
        force += between * (mode_y_[s + 1][n] - mode_y_[s][n]);
        stiffness += between;
      }
      if (touching) {
        for (int point = 0; point < kStringPoints; ++point) {
          std::size_t j = static_cast<std::size_t>(point);
          force += contact_force[j] * shape_[n][j];
        }
      }
      force -= 2.0 * mode_loss_[n] * std::sqrt(modal_mass * stiffness) * mode_v_[s][n];
      mode_v_[s][n] += h * force / modal_mass;
    }
  }
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    for (int mode = 0; mode < modes_; ++mode) {
      std::size_t n = static_cast<std::size_t>(mode);
      mode_y_[s][n] += h * mode_v_[s][n];
    }
  }
  string_openings();
}

void PinchTube::step(double wave_below, double area_below, double wave_above, double area_above, double& out_below,
                     double& out_above) {
  // Tension from the intended pitch (glides over about 20 ms), and each
  // section's elasticity wander (one pole at wander_corner Hz), advanced
  // once per step.
  double wanted = (pitch_ / settings_.natural_pitch) * (pitch_ / settings_.natural_pitch);
  tension_ += (1.0 - std::exp(-kStepSeconds / 0.02)) * (wanted - tension_);
  if (settings_.wander > 0.0) {
    double pole = std::exp(-kTwoPi * settings_.wander_corner * kStepSeconds);
    double innovation = settings_.wander * std::sqrt(1.0 - pole * pole);
    for (int section = 0; section < kSections; ++section) {
      double m = 0.0;
      double e = 0.0;
      double k = 0.0;
      section_factors(section, m, e, k);
      double tamp = std::exp(-settings_.wander_tamping * (k - 1.0));
      std::size_t s = static_cast<std::size_t>(section);
      drift_[s] = pole * drift_[s] + innovation * tamp * gauss();
    }
  }
  // The flow through the chain from the arriving waves:
  //   2 (w_below - w_above) = (Z_below + Z_above + R + 2 L / h) U + K U |U| + (air mass memory),
  // L the air's mass in the sections, R the slits' friction, K the kinetic
  // pressure into the narrowest section, less the jet's recovery past it.
  double z_below = kRhoC / std::max(area_below, kMinimumArea);
  double z_above = kRhoC / std::max(area_above, kMinimumArea);
  double inertance = 0.0;
  double viscous = 0.0;
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double a = std::max(area_[s] + leak_, kMinimumArea);
    inertance += kAirDensity * settings_.section_length / std::max(a, kInertanceFloor);
    viscous += 12.0 * kAirViscosity * settings_.section_length * settings_.width * settings_.width / (a * a * a);
  }
  double pinch = narrowest_ + leak_;
  double kinetic = half_rho_over(pinch) - half_rho_over(area_below) -
                   kJetRecovery * (half_rho_over(pinch) - half_rho_over(area_above));
  kinetic = std::max(kinetic, 0.1 * half_rho_over(pinch));
  double mass_term = 2.0 * inertance / kStepSeconds;
  double linear = z_below + z_above + viscous + mass_term;
  double drive = 2.0 * (wave_below - wave_above) + mass_term * flow_ + kInertiaDecay * inertial_pressure_;
  double magnitude = 2.0 * std::abs(drive) / (linear + std::sqrt(linear * linear + 4.0 * kinetic * std::abs(drive)));
  double flow = drive >= 0.0 ? magnitude : -magnitude;
  inertial_pressure_ = mass_term * (flow - flow_) - kInertiaDecay * inertial_pressure_;
  flow_ = flow;
  double p_in = 2.0 * wave_below - z_below * flow_;
  double p_out = 2.0 * wave_above + z_above * flow_;
  out_below = p_in - wave_below;
  out_above = p_out - wave_above;
  pressure_below_ = p_in;
  pressure_above_ = p_out;
  // The walls under these pressures, and past the pinch the jet's eddies:
  // their pressure scales with the jet's kinetic pressure, so the walls are
  // shaken hardest while the jet runs fastest and not at all while shut.
  pressures(p_in, p_out, area_below);
  if (settings_.jet_turbulence > 0.0) {
    double pole = std::exp(-kTwoPi * settings_.turbulence_corner * kStepSeconds);
    double innovation = std::sqrt(1.0 - pole * pole);
    double jet = flow_ * flow_ * half_rho_over(narrowest_ + leak_);
    for (int section = narrowest_index_ + 1; section < kSections; ++section) {
      std::size_t s = static_cast<std::size_t>(section);
      eddy_[s] = pole * eddy_[s] + innovation * gauss();
      pressure_[s] += settings_.jet_turbulence * jet * eddy_[s];
    }
  }
  double h = kStepSeconds / static_cast<double>(kSubsteps);
  for (int substep = 0; substep < kSubsteps; ++substep) {
    if (modes_ > 0) {
      advance_strings(h);
    } else {
      advance_walls(h);
    }
  }
  contact_ = 0.0;
  for (int section = 0; section < kSections; ++section) {
    contact_ = std::max(contact_, std::max(-opening_[static_cast<std::size_t>(section)], 0.0) * 0.5);
  }
}

}  // namespace tract
