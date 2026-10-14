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
  double command = std::max(1.0 + (pitch_ - kNaturalPitch) / kPitchPerTension, 0.05) * tension_factor_ * stretch;
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

void Cord::step(double wave_below, double area_below, double wave_above, double area_above, double& out_below,
                double& out_above) {
  advance_posture();
  advance_mechanics();
  compute_geometry();

  // Flow through the glottis from the incident waves:
  //   2 (w_below - w_above) = (Z_below + Z_above + R_viscous) U + kinetic U|U|.
  // Kinetic: full Bernoulli loss into each contraction, and only 20%
  // recovered at each expansion leaving the glottis (VocalTractLab).
  double z_below = kRhoC / std::max(area_below, kMinimumArea);
  double z_above = kRhoC / std::max(area_above, kMinimumArea);
  double a0 = std::max(section_area_[0], kMinimumArea);
  double a1 = std::max(section_area_[1], kMinimumArea);
  double kinetic = half_rho_over(a0) - half_rho_over(area_below);
  if (a1 < a0) {
    kinetic += half_rho_over(a1) - half_rho_over(a0);
  } else {
    kinetic -= kExitRecovery * (half_rho_over(a0) - half_rho_over(a1));
  }
  kinetic -= kExitRecovery * (half_rho_over(a1) - half_rho_over(area_above));
  kinetic = std::max(kinetic, 0.1 * half_rho_over(area_));
  // Viscous loss of the slit (Poiseuille between parallel plates).
  double viscous = 12.0 * kAirViscosity * length_ * length_ *
                   (thickness_[0] / (a0 * a0 * a0) + thickness_[1] / (a1 * a1 * a1));
  // The air inside the glottis has mass: rho l / A for each edge's section
  // and the short inlet and outlet. Implicit in the flow, so the flow lags
  // the opening (the air column must be accelerated).
  double inertance = kAirDensity * (thickness_[0] / a0 + thickness_[1] / a1 +
                                    (kInletLength + kOutletLength) / area_);
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
  pressure_lower_ = p_below - speed_term * (half_rho_over(a0) - half_rho_over(area_below));
  if (a1 < a0) {
    pressure_upper_ = pressure_lower_ - speed_term * (half_rho_over(a1) - half_rho_over(a0));
  } else {
    pressure_upper_ = pressure_lower_ + kExitRecovery * speed_term * (half_rho_over(a0) - half_rho_over(a1));
  }
}

}  // namespace tract
