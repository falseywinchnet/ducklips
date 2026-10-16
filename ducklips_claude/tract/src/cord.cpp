#include "tract/cord.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

// Static fold properties (VocalTractLab 2.4 TriangularGlottis defaults).
constexpr double kRestThickness[2] = {0.24, 0.06};      // cm
constexpr double kMass[2] = {0.12, 0.03};               // g
constexpr double kSpring[2] = {80000.0, 8000.0};        // dyn/cm
constexpr double kContactSpring[2] = {240000.0, 24000.0};
constexpr double kCoupling = 25000.0;                   // dyn/cm, constant in tension
constexpr double kRestLength = 1.3;                     // cm
constexpr double kInletLength = 0.05;                   // cm
constexpr double kOutletLength = 0.01;                  // cm
constexpr double kNaturalPitch = 130.98;                // Hz at Q = 1 (M01 speaker file)
constexpr double kPitchPerTension = 132.07;             // Hz per unit Q
constexpr double kExitRecovery = 0.2;                   // share of kinetic pressure recovered leaving the glottis

double half_rho_over(double area) {
  double a = std::max(area, kMinimumArea);
  return 0.5 * kAirDensity / (a * a);
}

}  // namespace

CordPosture modal_posture() { return CordPosture{0.00495, 0.0004, 0.0002}; }
CordPosture breathy_posture() { return CordPosture{0.012, 0.011, 0.004}; }
CordPosture voiced_fricative_posture() { return CordPosture{0.02035, 0.01965, 0.0002}; }
CordPosture voiceless_plosive_posture() { return CordPosture{0.0501, 0.0501, 0.0002}; }
CordPosture aspirate_posture() { return CordPosture{0.0802, 0.07985, 0.0104}; }
CordPosture voiceless_fricative_posture() { return CordPosture{0.10225, 0.10225, 0.0002}; }

Cord::Cord() : posture_(modal_posture()), guide_(modal_posture()), target_(modal_posture()) {}

void Cord::set_target(const CordPosture& posture, double move_seconds) {
  target_ = posture;
  // First-order approach; move_seconds is about three time constants.
  // Two stages of half the time constant each (see advance_posture).
  double tau = std::max(move_seconds / 3.0, kStepSeconds);
  glide_gain_ = 1.0 - std::exp(-kStepSeconds / (0.5 * tau));
}

void Cord::advance_posture() {
  // Every command arrives once per 1 ms tick. Followed by one first-order
  // stage, its velocity jumped at every tick: with the folds held shut while
  // the breath rose, the leak through them carried tones at 1, 2 and 3 kHz
  // out through the nose. Two stages (each half the time constant) keep the
  // velocity continuous.
  guide_.lower_rest += glide_gain_ * (target_.lower_rest - guide_.lower_rest);
  guide_.upper_rest += glide_gain_ * (target_.upper_rest - guide_.upper_rest);
  guide_.chink += glide_gain_ * (target_.chink - guide_.chink);
  posture_.lower_rest += glide_gain_ * (guide_.lower_rest - posture_.lower_rest);
  posture_.upper_rest += glide_gain_ * (guide_.upper_rest - posture_.upper_rest);
  posture_.chink += glide_gain_ * (guide_.chink - posture_.chink);
  double offset_gain = 1.0 - std::exp(-kStepSeconds / 0.001);
  rest_offset_guide_ += offset_gain * (rest_offset_target_ - rest_offset_guide_);
  rest_offset_ += offset_gain * (rest_offset_guide_ - rest_offset_);
  // Tension from the intended pitch (open loop, as calibrated). The muscles'
  // unsteadiness arrives through tension_factor_ (slow variation and tremor
  // from the body); Klatt's three-sine flutter was a repeating pattern and is
  // gone.
  // The muscle's force cannot jump: tension follows its command with a 3 ms
  // time constant (commands arrive at the 1 kHz control rate, and stepping
  // the tension there modulated the voice at 1 kHz).
  // The swing stretches the fold along its length: the peak of each swing,
  // held and released over elongation_release_ (seconds, about a cycle).
  deflection_ = std::max(x_[0] * x_[0], deflection_ * std::exp(-kStepSeconds / elongation_release_));
  // Below the reference swing the stretch fades only to elongation_floor_:
  // as the breath dies the folds stop before their pitch falls further.
  double stretch = std::max((1.0 + elongation_ * deflection_) / (1.0 + elongation_ * elongation_reference_),
                            elongation_floor_);
  double natural = fold_tube_ ? fold_.natural_pitch : kNaturalPitch;
  double per_tension = fold_tube_ ? fold_.pitch_per_tension : kPitchPerTension;
  double anticipated = 1.0;
  if (fold_tube_ && elongation_ > 0.0) {
    double swing = std::max(fold_.swing_per_pascal * (expected_breath_ - fold_.swing_onset), 0.0);
    anticipated = std::max((1.0 + elongation_ * swing) / (1.0 + elongation_ * elongation_reference_),
                           elongation_floor_);
  }
  double command =
      std::max(1.0 + (pitch_ - natural) / per_tension, 0.05) * tension_factor_ * stretch / anticipated;
  double tension_gain = 1.0 - std::exp(-kStepSeconds / 0.0015);
  tension_guide_ += tension_gain * (command - tension_guide_);
  tension_ += tension_gain * (tension_guide_ - tension_);
  length_ = kRestLength * std::sqrt(tension_);
  thickness_[0] = kRestThickness[0] / std::sqrt(tension_);
  thickness_[1] = kRestThickness[1] / std::sqrt(tension_);
}

// Open and closed parts of each edge on the wedge-shaped rest gap.
void Cord::compute_geometry() {
  double rest[2] = {posture_.lower_rest + rest_offset_, posture_.upper_rest + rest_offset_};
  for (int edge = 0; edge < 2; ++edge) {
    double back = rest[edge] + x_[edge];  // at the vocal processes
    double front = rest[edge] >= 0.0 ? x_[edge] : back;  // at the anterior commissure
    double contact_length = 0.0;
    double mean_contact_z = 0.5 * length_;
    if (back > 0.0 && front > 0.0) {
      open_length_[edge] = length_;
      mean_width_[edge] = back + front;
    } else if (back <= 0.0 && front <= 0.0) {
      open_length_[edge] = 0.0;
      mean_width_[edge] = 0.0;
      contact_length = length_;
    } else {
      // Open toward the back, closed toward the front.
      double apex = rest[edge] > 0.0 ? length_ * (1.0 + x_[edge] / rest[edge]) : 0.0;
      apex = std::clamp(apex, 0.0, length_);
      open_length_[edge] = apex;
      mean_width_[edge] = std::max(back, 0.0);
      contact_length = length_ - apex;
      mean_contact_z = 0.5 * (apex + length_);
    }
    alpha_[edge] = contact_length / length_;
    rest_star_[edge] = rest[edge] >= 0.0 ? rest[edge] * (1.0 - mean_contact_z / length_) : rest[edge];
    section_area_[edge] = open_length_[edge] * mean_width_[edge] + std::max(posture_.chink, 0.0);
  }
  area_ = std::max(std::min(section_area_[0], section_area_[1]), kMinimumArea);
}

// Backward Euler for both edges together: springs, coupling, contact and
// damping implicit, pressure forces from the previous step's flow.
void Cord::advance_mechanics() {
  constexpr double h = kStepSeconds;
  double mass[2];
  double spring[2];
  double contact[2];
  double damping[2];
  double force[2];
  for (int edge = 0; edge < 2; ++edge) {
    mass[edge] = kMass[edge] / tension_;
    spring[edge] = kSpring[edge] * tension_;
    contact[edge] = kContactSpring[edge] * tension_;
    // Absolute damping does not change with tension (m k is constant);
    // contact adds up to one more critical damping as the edges press.
    damping[edge] = 2.0 * (damping_ratio_[edge] + alpha_[edge]) * std::sqrt(kMass[edge] * kSpring[edge]);
  }
  force[0] = pressure_lower_ * open_length_[0] * thickness_[0] +
             0.25 * (pressure_sub_ + pressure_lower_) * kInletLength * length_;
  force[1] = pressure_upper_ * open_length_[1] * thickness_[1] +
             0.25 * (pressure_supra_ + pressure_upper_) * kOutletLength * length_ +
             surface_pressure_ * (open_length_[1] * thickness_[1] + 0.5 * kOutletLength * length_);
  double a = mass[0] + damping[0] * h + h * h * (spring[0] + contact[0] * alpha_[0]) + kCoupling * h * h;
  double b = -kCoupling * h * h;
  double d = mass[1] + damping[1] * h + h * h * (spring[1] + contact[1] * alpha_[1]) + kCoupling * h * h;
  double e = force[0] * h * h + 2.0 * mass[0] * x_[0] - mass[0] * x_previous_[0] + damping[0] * h * x_[0] -
             h * h * contact[0] * alpha_[0] * rest_star_[0];
  double f = force[1] * h * h + 2.0 * mass[1] * x_[1] - mass[1] * x_previous_[1] + damping[1] * h * x_[1] -
             h * h * contact[1] * alpha_[1] * rest_star_[1];
  double determinant = a * d - b * b;
  double next0 = (e * d - b * f) / determinant;
  double next1 = (a * f - b * e) / determinant;
  x_previous_[0] = x_[0];
  x_previous_[1] = x_[1];
  x_[0] = next0;
  x_[1] = next1;
}

// Unit-variance noise (sum of four uniforms) for the kinematic jitter.
double Cord::jitter_gauss() {
  double sum = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    jitter_random_ ^= jitter_random_ << 13U;
    jitter_random_ ^= jitter_random_ >> 17U;
    jitter_random_ ^= jitter_random_ << 5U;
    sum += static_cast<double>(jitter_random_) / 4294967296.0;
  }
  return (sum - 2.0) * std::sqrt(3.0);
}

// Kinematic folds: prescribed edge swings (see set_kinematic), the opening
// integrated along the fold with a soft contact.
void Cord::advance_kinematics() {
  constexpr double kKinematicRestLength = 1.6;   // cm at kKinematicRestPitch
  constexpr double kKinematicRestPitch = 120.0;  // Hz
  constexpr double kKinematicThickness = 0.45;   // cm at rest length
  constexpr int kSamples = 16;
  double frequency = std::max(pitch_ * tension_factor_, 40.0);
  length_ = kKinematicRestLength * std::sqrt(frequency / kKinematicRestPitch);
  double thickness = kKinematicThickness * kKinematicRestLength / length_;
  thickness_[0] = 0.5 * thickness;
  thickness_[1] = 0.5 * thickness;
  // The swing follows the transglottal pressure, smoothed (two poles, 25 Hz).
  double pole = 1.0 - std::exp(-kTwoPi * 25.0 * kStepSeconds);
  transglottal_low_ += pole * ((pressure_sub_ - pressure_supra_) - transglottal_low_);
  transglottal_lower_ += pole * (transglottal_low_ - transglottal_lower_);
  double drive = std::max(transglottal_lower_, 0.0);
  double amplitude = amplitude_gain_ * std::sqrt(drive);
  amplitude *= std::clamp((drive - 1000.0) / 1000.0, 0.0, 1.0);  // 100-200 Pa
  phase_ += kTwoPi * frequency * kStepSeconds + jitter_ * std::sqrt(kStepSeconds) * jitter_gauss();
  if (phase_ > kTwoPi) {
    phase_ -= kTwoPi;
  }
  double rest[2] = {posture_.lower_rest + rest_offset_, posture_.upper_rest + rest_offset_};
  double edge_phase[2] = {phase_, phase_ - phase_lag_};
  double softness = std::max(contact_softness_, 1.0e-5);
  double step = length_ / kSamples;
  for (int edge = 0; edge < 2; ++edge) {
    double swing = amplitude * std::cos(edge_phase[edge]);
    x_[edge] = swing;
    double open_sum = 0.0;
    double closed = 0.0;
    for (int sample = 0; sample < kSamples; ++sample) {
      double z = (sample + 0.5) / kSamples;  // anterior commissure 0 .. vocal process 1
      double half_gap = rest[edge] * z + swing * std::sin(kPi * z);
      // Soft contact: the half gap passes smoothly through zero to zero
      // (tissue compresses; no corner at closure).
      double soft = softness * std::log1p(std::exp(std::min(half_gap / softness, 50.0)));
      if (half_gap / softness > 50.0) {
        soft = half_gap;
      }
      open_sum += 2.0 * soft * step;
      if (half_gap <= 0.0) {
        closed += 1.0;
      }
    }
    alpha_[edge] = closed / kSamples;
    open_length_[edge] = length_ * (1.0 - alpha_[edge]);
    section_area_[edge] = open_sum + std::max(posture_.chink, 0.0);
  }
  area_ = std::max(std::min(section_area_[0], section_area_[1]), kMinimumArea);
}

// Each slice's opening along the fold on the wedge-shaped rest gap (zero at
// the anterior commissure), interpolated from the lower edge's rest at the
// bottom slice to the upper edge's at the top: as compute_geometry, per slice.
void Cord::fold_geometry() {
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    double share = static_cast<double>(slice) / static_cast<double>(kFoldSlices - 1);
    double rest = posture_.lower_rest + share * (posture_.upper_rest - posture_.lower_rest) + rest_offset_;
    double x = cover_q_[slice];
    double back = rest + x;
    double front = rest >= 0.0 ? x : back;
    double contact_length = 0.0;
    double mean_contact_z = 0.5 * length_;
    double open_length = 0.0;
    double mean_width = 0.0;
    if (back > 0.0 && front > 0.0) {
      open_length = length_;
      mean_width = back + front;
    } else if (back <= 0.0 && front <= 0.0) {
      contact_length = length_;
    } else {
      double apex = rest > 0.0 ? length_ * (1.0 + x / rest) : 0.0;
      apex = std::clamp(apex, 0.0, length_);
      open_length = apex;
      mean_width = std::max(back, 0.0);
      contact_length = length_ - apex;
      mean_contact_z = 0.5 * (apex + length_);
    }
    slice_open_[slice] = open_length;
    slice_alpha_[slice] = contact_length / length_;
    slice_star_[slice] = rest >= 0.0 ? rest * (1.0 - mean_contact_z / length_) : rest;
    slice_area_[slice] = open_length * mean_width + std::max(posture_.chink, 0.0);
  }
}

double Cord::contact() const {
  if (!fold_tube_) {
    return alpha_[0];
  }
  double most = 0.0;
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    most = std::max(most, slice_alpha_[slice]);
  }
  return most;
}

// The slices' openings as the air sees them: the half gap along the fold
// (the wedge rest plus the slice's displacement) through a soft contact,
// summed over eight points along the length; once per step, after the
// pieces have moved. (The pieces' contact forces use the plain geometry.)
void Cord::fold_air_areas() {
  constexpr int kPoints = 8;
  double softness = std::max(fold_.contact_softness, 1.0e-6);
  double step = length_ / static_cast<double>(kPoints);
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    double share = static_cast<double>(slice) / static_cast<double>(kFoldSlices - 1);
    double rest = posture_.lower_rest + share * (posture_.upper_rest - posture_.lower_rest) + rest_offset_;
    double x = cover_q_[slice];
    double open_sum = 0.0;
    for (int point = 0; point < kPoints; ++point) {
      double z = (static_cast<double>(point) + 0.5) / static_cast<double>(kPoints);  // 0 front .. 1 vocal process
      double half_gap = (rest >= 0.0 ? rest * z : rest) + x;
      double ratio = half_gap / softness;
      double soft = ratio > 30.0 ? half_gap : softness * std::log1p(std::exp(std::max(ratio, -30.0)));
      open_sum += 2.0 * soft * step;
    }
    slice_area_[slice] = open_sum + std::max(posture_.chink, 0.0);
  }
}

// The fold tube's pieces, eight steps per step (1.3 us), symplectic Euler
// (every mode under 2 kHz: omega h < 0.02). Forces from the air's pressures
// of the last step; contact and its damping where the slice is closed.
void Cord::advance_fold_tube() {
  double q = std::max(tension_, 0.05);
  double slices = static_cast<double>(kFoldSlices);
  double slice_depth = fold_.depth / std::sqrt(q) / slices;
  double cover_mass = fold_.cover_mass / (slices * q);
  double body_mass = fold_.body_mass / (slices * q);
  double cover_k = fold_.cover_stiffness * q / slices;
  double coupling_k = fold_.coupling_stiffness * q / slices;
  double support_k = fold_.support_stiffness * q / slices;
  double vertical_k = fold_.vertical_stiffness * q;
  double contact_k = fold_.contact_stiffness * q / slices;
  double cover_c = 2.0 * fold_.cover_loss * std::sqrt(cover_mass * (cover_k + coupling_k));
  double support_c = 2.0 * fold_.support_loss * std::sqrt(body_mass * support_k);
  double contact_c = 2.0 * std::sqrt(cover_mass * contact_k);
  constexpr double h = kStepSeconds / static_cast<double>(kFoldSubsteps);
  double pole = fold_.muscle_corner > 0.0 ? std::exp(-kTwoPi * fold_.muscle_corner * h) : 0.0;
  double innovation = fold_.muscle_fluctuation * std::sqrt(1.0 - pole * pole);
  double cover_force[kFoldSlices];
  double body_force[kFoldSlices];
  for (int substep = 0; substep < kFoldSubsteps; ++substep) {
    if (fold_.muscle_fluctuation > 0.0) {
      muscle_ = pole * muscle_ + innovation * jitter_gauss();
    }
    double muscle_k = support_k * (1.0 + muscle_);
    fold_geometry();
    for (int slice = 0; slice < kFoldSlices; ++slice) {
      double x = cover_q_[slice];
      double force = slice_pressure_[slice] * slice_open_[slice] * slice_depth;
      if (slice == 0) {
        force += 0.25 * (pressure_sub_ + slice_pressure_[0]) * kInletLength * length_;
      }
      if (slice == kFoldSlices - 1) {
        force += 0.25 * (pressure_supra_ + slice_pressure_[slice]) * kOutletLength * length_ +
                 surface_pressure_ * (slice_open_[slice] * slice_depth + 0.5 * kOutletLength * length_);
      }
      force -= cover_k * x + coupling_k * (x - body_q_[slice]);
      if (slice > 0) {
        force -= vertical_k * (x - cover_q_[slice - 1]);
      }
      if (slice + 1 < kFoldSlices) {
        force -= vertical_k * (x - cover_q_[slice + 1]);
      }
      double alpha = slice_alpha_[slice];
      if (alpha > 0.0) {
        force -= contact_k * alpha * (x + slice_star_[slice]) + contact_c * alpha * cover_v_[slice];
      }
      double relative = cover_v_[slice] - body_v_[slice];
      cover_force[slice] = force - cover_c * relative;
      body_force[slice] = -muscle_k * body_q_[slice] - coupling_k * (body_q_[slice] - x) - support_c * body_v_[slice] +
                          cover_c * relative;
    }
    for (int slice = 0; slice < kFoldSlices; ++slice) {
      cover_v_[slice] += h * cover_force[slice] / cover_mass;
      body_v_[slice] += h * body_force[slice] / body_mass;
      cover_q_[slice] += h * cover_v_[slice];
      body_q_[slice] += h * body_v_[slice];
    }
  }
  fold_geometry();
  fold_air_areas();
  double narrowest = slice_area_[0];
  double closed = 0.0;
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    narrowest = std::min(narrowest, slice_area_[slice]);
    closed = std::max(closed, slice_alpha_[slice]);
  }
  x_[0] = cover_q_[0];
  x_[1] = cover_q_[kFoldSlices - 1];
  section_area_[0] = slice_area_[0];
  section_area_[1] = slice_area_[kFoldSlices - 1];
  alpha_[0] = closed;
  alpha_[1] = slice_alpha_[kFoldSlices - 1];
  open_length_[0] = slice_open_[0];
  open_length_[1] = slice_open_[kFoldSlices - 1];
  thickness_[0] = 0.5 * fold_.depth / std::sqrt(q);
  thickness_[1] = thickness_[0];
  area_ = std::max(narrowest, kMinimumArea);
}

// The air's pressure on each slice: Bernoulli's from the inlet while the
// jet holds to the walls, the outlet's past the point where it separates
// (the first slice beyond the narrowest that is 20% wider; none: the exit).
// Suction is held to one transglottal drop below the outlet (in a gap of a
// few microns the flow is viscous and Bernoulli's law no longer holds).
void Cord::fold_pressures(double flow, double p_below, double p_above, double area_below) {
  int narrowest = 0;
  for (int slice = 1; slice < kFoldSlices; ++slice) {
    if (slice_area_[slice] < slice_area_[narrowest]) {
      narrowest = slice;
    }
  }
  int separation = kFoldSlices;
  for (int slice = narrowest + 1; slice < kFoldSlices; ++slice) {
    if (slice_area_[slice] > 1.2 * slice_area_[narrowest]) {
      separation = slice;
      break;
    }
  }
  double speed_term = flow * std::abs(flow);
  double floor = p_above - std::abs(p_below - p_above);
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    double pressure = slice < separation
                          ? p_below - speed_term * (half_rho_over(slice_area_[slice]) - half_rho_over(area_below))
                          : p_above;
    slice_pressure_[slice] = std::max(pressure, floor);
  }
}

void Cord::step(double wave_below, double area_below, double wave_above, double area_above, double& out_below,
                double& out_above) {
  advance_posture();
  if (fold_tube_) {
    advance_fold_tube();
  } else if (kinematic_) {
    advance_kinematics();
  } else {
    advance_mechanics();
    compute_geometry();
  }

  // Flow through the glottis from the incident waves:
  //   2 (w_below - w_above) = (Z_below + Z_above + R_viscous) U + kinetic U|U|.
  // Kinetic: full Bernoulli loss into each contraction, and only 20%
  // recovered at each expansion leaving the glottis (VocalTractLab).
  double z_below = kRhoC / std::max(area_below, kMinimumArea);
  double z_above = kRhoC / std::max(area_above, kMinimumArea);
  double a0 = std::max(section_area_[0], kMinimumArea);
  double a1 = std::max(section_area_[1], kMinimumArea);
  double kinetic = 0.0;
  double viscous = 0.0;
  if (fold_tube_) {
    // Into the narrowest slice without loss, out of it as a jet (20%
    // recovered); Poiseuille along every slice.
    kinetic = half_rho_over(area_) - half_rho_over(area_below) -
              kExitRecovery * (half_rho_over(area_) - half_rho_over(area_above));
    double slice_depth = 2.0 * thickness_[0] / static_cast<double>(kFoldSlices);
    for (int slice = 0; slice < kFoldSlices; ++slice) {
      double a = std::max(slice_area_[slice], kMinimumArea);
      viscous += 12.0 * kAirViscosity * length_ * length_ * slice_depth / (a * a * a);
    }
  } else {
    kinetic = half_rho_over(a0) - half_rho_over(area_below);
    if (a1 < a0) {
      kinetic += half_rho_over(a1) - half_rho_over(a0);
    } else {
      kinetic -= kExitRecovery * (half_rho_over(a0) - half_rho_over(a1));
    }
    kinetic -= kExitRecovery * (half_rho_over(a1) - half_rho_over(area_above));
    // Viscous loss of the slit (Poiseuille between parallel plates).
    viscous = 12.0 * kAirViscosity * length_ * length_ *
              (thickness_[0] / (a0 * a0 * a0) + thickness_[1] / (a1 * a1 * a1));
  }
  kinetic = std::max(kinetic, 0.1 * half_rho_over(area_));
  // The air inside the glottis has mass: rho l / A for each edge's section
  // and the short inlet and outlet. Implicit in the flow, so the flow lags
  // the opening (the air column must be accelerated).
  // (The air mass is taken on areas of at least kInertanceFloor: as the
  // folds shut, rho l / A grew without bound within a step and the
  // trapezoidal step answered with a one-sample pressure spike at every
  // closure, heard as a crackle through voiced closures and nasals.)
  double inertance = 0.0;
  if (fold_tube_) {
    double slice_depth = 2.0 * thickness_[0] / static_cast<double>(kFoldSlices);
    for (int slice = 0; slice < kFoldSlices; ++slice) {
      inertance += kAirDensity * slice_depth / std::max(slice_area_[slice], inertance_floor_);
    }
    inertance += kAirDensity * (kInletLength + kOutletLength) / std::max(area_, inertance_floor_);
  } else {
    inertance = kAirDensity * (thickness_[0] / std::max(a0, inertance_floor_) +
                               thickness_[1] / std::max(a1, inertance_floor_) +
                               (kInletLength + kOutletLength) / std::max(area_, inertance_floor_));
  }
  // Trapezoidal (lossless) step for the air mass: its pressure
  //   p_n = (2 L / h)(U_n - U_{n-1}) - g p_{n-1},
  // with g slightly below 1 so a shut glottis (huge L, U held at 0) cannot
  // keep a Nyquist-rate pressure alive. (Backward Euler, as first written,
  // turns the inertance into a resistance growing as f^2: about 42 dyn s/cm^5
  // at 5 kHz, a third of the tract's own impedance, and the voiced spectrum
  // lost 20-30 dB at 4-7.5 kHz while phonating against the same body at rest.)
  constexpr double kInertiaDecay = 0.95;
  double mass_term = 2.0 * inertance / kStepSeconds;
  double linear = z_below + z_above + viscous + mass_term;
  double drive = 2.0 * (wave_below - wave_above) + mass_term * flow_ + kInertiaDecay * inertial_pressure_;
  double magnitude =
      2.0 * std::abs(drive) / (linear + std::sqrt(linear * linear + 4.0 * kinetic * std::abs(drive)));
  double flow = drive >= 0.0 ? magnitude : -magnitude;
  inertial_pressure_ = mass_term * (flow - flow_) - kInertiaDecay * inertial_pressure_;
  flow_ = flow;
  double heard = flow + flux_noise_;  // the flux the airways receive
  double p_below = 2.0 * wave_below - z_below * heard;
  // (Pressure at the glottal inlet and outlet; the inertial pressure is
  // taken across the glottis itself.)
  double p_above = 2.0 * wave_above + z_above * heard;
  out_below = p_below - wave_below;
  out_above = p_above - wave_above;

  // Pressures on the edges for the next step (the open parts feel them).
  double speed_term = flow * std::abs(flow);
  pressure_sub_ = p_below;
  pressure_supra_ = p_above;
  if (fold_tube_) {
    fold_pressures(flow, p_below, p_above, area_below);
    return;
  }
  pressure_lower_ = p_below - speed_term * (half_rho_over(a0) - half_rho_over(area_below));
  if (a1 < a0) {
    pressure_upper_ = pressure_lower_ - speed_term * (half_rho_over(a1) - half_rho_over(a0));
  } else {
    pressure_upper_ = pressure_lower_ + kExitRecovery * speed_term * (half_rho_over(a0) - half_rho_over(a1));
  }
}

}  // namespace tract
