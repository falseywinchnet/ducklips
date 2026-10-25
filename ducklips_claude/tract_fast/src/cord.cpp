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

// One edge's swing (share of the amplitude) under VocalTractLab 2019's skew
// and double pulsing, from its phase in the cosine convention (peak at 0)
// over the phase's 4 pi: their phase alpha = cosine phase + pi / 2, taken
// from the trough at 1.5 pi; each cycle warped at its knee; the first cycle
// (1.5 pi to 3.5 pi) the double-pulsed one.
double Cord::shaped_swing(double cosine_phase) const {
  constexpr double kFourPi = 2.0 * kTwoPi;
  double from_trough = std::fmod(cosine_phase + 0.5 * kPi - 1.5 * kPi, kFourPi);
  if (from_trough < 0.0) {
    from_trough += kFourPi;
  }
  int cycle = from_trough < kTwoPi ? 0 : 1;
  double within = from_trough - static_cast<double>(cycle) * kTwoPi;
  double knee = kPi * (1.0 + skew_);
  double warped = within < knee ? within * kPi / knee : kPi + (within - knee) * kPi / (kTwoPi - knee);
  double phase = 1.5 * kPi + static_cast<double>(cycle) * kTwoPi + warped;
  double swing = std::sin(phase);
  if (cycle == 0) {
    swing = (1.0 - 0.5 * double_pulsing_) * (swing + 1.0) - 1.0;
  }
  return swing;
}

// Geometric glottis 2025: one edge's displacement (cm) at its phase (cosine
// convention, peak at 0), with the exponential collision and pulse shape
// (see set_geometric_2025).
double Cord::collided_swing(double cosine_phase, double amplitude, double rest) const {
  constexpr double kEpsilon = 1.0e-6;
  double t = std::fmod(cosine_phase, kTwoPi);
  if (t < 0.0) {
    t += kTwoPi;
  }
  t /= kTwoPi;
  double value = amplitude * std::cos(kTwoPi * t);
  if (amplitude <= 0.0) {
    return value;
  }
  double touch = -std::max(rest, 0.0) / kPi;  // A2: where the folds begin to touch
  double settle = touch - 0.2 * amplitude;    // A3
  double pulse_shape = std::clamp(pulse_shape_ + cycle_now_[0], -1.0, 1.0);
  double t2 = 0.5 - kEpsilon;
  double t3 = 0.5 + kEpsilon;
  if (touch > -amplitude) {
    t2 = std::acos(touch / amplitude) / kTwoPi;
    t3 = 1.0 - t2;
  }
  double t1 = t2;
  double t4 = t3;
  if (pulse_shape < 0.0) {
    t1 = t2 + 0.1 * pulse_shape;
  } else {
    t4 = t3 + 0.1 * pulse_shape;
  }
  if (t < t1 || t >= t4) {
    return value;
  }
  double value_t1 = amplitude * std::cos(kTwoPi * t1);
  double slope_t1 = -amplitude * kTwoPi * std::sin(kTwoPi * t1);
  double tau_t1 = -(value_t1 - settle) / slope_t1;
  double value_t4 = amplitude * std::cos(kTwoPi * t4);
  double slope_t4 = -amplitude * kTwoPi * std::sin(kTwoPi * t4);
  double tau_t4 = (value_t4 - settle) / slope_t4;
  double closing = settle + (value_t1 - settle) * std::exp(-(t - t1) / tau_t1);
  double opening = settle + (value_t4 - settle) * std::exp(-(t4 - t) / tau_t4);
  return std::max(value, std::max(closing, opening));
}

// The pressure above the glottis through two second-order Butterworth
// low-pass sections at 25 Hz (bilinear, prewarped), for the swing.
double Cord::supraglottal_low(double pressure) {
  double k = std::tan(kPi * 25.0 * kStepSeconds);
  double norm = 1.0 / (1.0 + std::sqrt(2.0) * k + k * k);
  double b0 = k * k * norm;
  double a1 = 2.0 * (k * k - 1.0) * norm;
  double a2 = (1.0 - std::sqrt(2.0) * k + k * k) * norm;
  double x = pressure;
  for (int section = 0; section < 2; ++section) {
    // Transposed direct form II: y = b0 x + s0; s0 = 2 b0 x - a1 y + s1; s1 = b0 x - a2 y.
    double y = b0 * x + supra_filter_[section][0];
    supra_filter_[section][0] = 2.0 * b0 * x - a1 * y + supra_filter_[section][1];
    supra_filter_[section][1] = b0 * x - a2 * y;
    x = y;
  }
  return x;
}

// Geometric glottis 2025: length, thickness and swing; each edge's collided
// displacement; the areas along the fold at 32 points cut at the midline.
void Cord::advance_geometric_2025(double frequency) {
  const double kRestLength2025 = rest_length_2025_;        // cm (JD2, M01: 1.6; W02: 1.2)
  const double kRestThickness2025 = rest_thickness_2025_;  // cm (0.45; 0.40)
  const double kRestPitch2025 = rest_pitch_2025_;          // Hz (120; 220)
  constexpr int kSamples = 32;
  length_ = kRestLength2025 * std::sqrt(frequency / kRestPitch2025);
  double thickness = kRestThickness2025 * kRestLength2025 / length_;
  thickness_[0] = 0.5 * thickness;
  thickness_[1] = 0.5 * thickness;
  double supra = supraglottal_low(pressure_supra_);
  double transglottal = std::max(expected_breath_ * kPascal - supra, 0.0);
  double amplitude = amplitude_gain_ * (kRestLength2025 / 1.6) * std::sqrt(transglottal);
  if (transglottal < 1000.0) {
    amplitude = 0.0;
  } else if (transglottal < 2000.0) {
    amplitude *= (transglottal - 1000.0) / 1000.0;
  }
  amplitude *= relative_amplitude_ * abduction_scale_;
  // This cycle's variation, gliding across the cycle (cycle_now_[0], the
  // pulse shape's, is read by collided_swing).
  double progress = std::clamp(phase_ / kTwoPi, 0.0, 1.0);
  for (int k = 0; k < 4; ++k) {
    cycle_now_[k] = cycle_from_[k] + progress * (cycle_to_[k] - cycle_from_[k]);
  }
  amplitude *= std::max(1.0 + cycle_now_[2], 0.0);
  double rest[2] = {posture_.lower_rest + rest_offset_ + cycle_now_[3], posture_.upper_rest + rest_offset_ + cycle_now_[3]};
  double edge_phase[2] = {phase_, phase_ + kTwoPi - (phase_lag_ + cycle_now_[1])};
  double step = length_ / static_cast<double>(kSamples - 1);
  for (int edge = 0; edge < 2; ++edge) {
    double displacement = collided_swing(edge_phase[edge], amplitude, rest[edge]);
    x_[edge] = displacement;
    double area = std::max(posture_.chink, 0.0);
    double previous = 0.0;
    int closed = 0;
    for (int sample = 1; sample <= kSamples; ++sample) {
      double z = static_cast<double>(sample) / static_cast<double>(kSamples);
      double x = rest[edge] * z + displacement * std::sin(kPi * z);
      if (x < 0.0) {
        ++closed;
      }
      if (softness_2025_ > 0.0) {
        x = x / softness_2025_ > 50.0 ? x : softness_2025_ * std::log1p(std::exp(x / softness_2025_));
      } else if (x < 0.0) {
        x = 0.0;
      }
      area += step * (previous + x);
      previous = x;
    }
    alpha_[edge] = static_cast<double>(closed) / static_cast<double>(kSamples);
    open_length_[edge] = length_ * (1.0 - alpha_[edge]);
    section_area_[edge] = area;
  }
  area_ = std::max(std::min(section_area_[0], section_area_[1]), kMinimumArea);
}

// Kinematic folds: prescribed edge swings (see set_kinematic), the opening
// integrated along the fold with a soft contact.
void Cord::advance_kinematics() {
  constexpr double kKinematicRestLength = 1.6;   // cm at kKinematicRestPitch
  constexpr double kKinematicRestPitch = 120.0;  // Hz
  constexpr double kKinematicThickness = 0.45;   // cm at rest length
  constexpr int kSamples = 16;
  double frequency = std::max(pitch_ * tension_factor_, 40.0);
  if (flutter_ > 0.0) {
    flutter_time_ += kStepSeconds;
    double waves = std::sin(kTwoPi * 12.7 * flutter_time_) + std::sin(kTwoPi * 7.1 * flutter_time_) +
                   std::sin(kTwoPi * 4.7 * flutter_time_);
    frequency = std::max(frequency + (flutter_ / 50.0) * (frequency / 100.0) * waves, 40.0);
  }
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
  if (length_scaling_ != 0.0) {
    amplitude *= std::pow(kKinematicRestLength / length_, length_scaling_);
  }
  phase_ += kTwoPi * frequency * kStepSeconds + jitter_ * std::sqrt(kStepSeconds) * jitter_gauss();
  if (phase_ > kTwoPi) {
    phase_ -= kTwoPi;
    second_cycle_ = !second_cycle_;
    if (cycle_sd_[0] > 0.0 || cycle_sd_[1] > 0.0 || cycle_sd_[2] > 0.0 || cycle_sd_[3] > 0.0) {
      for (int k = 0; k < 4; ++k) {
        cycle_from_[k] = cycle_to_[k];
        cycle_to_[k] = cycle_sd_[k] > 0.0 ? cycle_sd_[k] * jitter_gauss() : 0.0;
      }
    }
  }
  if (geometric_2025_) {
    advance_geometric_2025(frequency);
    return;
  }
  double rest[2] = {posture_.lower_rest + rest_offset_, posture_.upper_rest + rest_offset_};
  double edge_phase[2] = {phase_, phase_ - phase_lag_};
  double softness = std::max(contact_softness_, 1.0e-5);
  double step = length_ / kSamples;
  for (int edge = 0; edge < 2; ++edge) {
    double swing = amplitude * std::cos(edge_phase[edge]);
    if (skew_ != 0.0 || double_pulsing_ > 0.0) {
      swing = amplitude * shaped_swing(edge_phase[edge] + (second_cycle_ ? kTwoPi : 0.0));
    }
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

// The fold's parallel strips along its length (fold_.strips: 1 to
// kMaxStrips, front first), each with its share of the length and its own
// stiffness and mass factors (FoldTubeSettings::strip*_stiffness, _mass).
int Cord::strip_count() const {
  return std::clamp(static_cast<int>(std::lround(fold_.strips)), 1, kMaxStrips);
}

double Cord::strip_stiffness(int strip) const {
  if (strip_count() == 1) {
    return 1.0;
  }
  return strip == 0 ? fold_.strip1_stiffness : strip == 1 ? fold_.strip2_stiffness : fold_.strip3_stiffness;
}

void Cord::slice_factors(int slice, double& mass, double& elasticity, double& support) const {
  const double masses[kFoldSlices] = {fold_.slice1_mass, fold_.slice2_mass, fold_.slice3_mass,
                                      fold_.slice4_mass, fold_.slice5_mass, fold_.slice6_mass};
  const double elasticities[kFoldSlices] = {fold_.slice1_elasticity, fold_.slice2_elasticity, fold_.slice3_elasticity,
                                            fold_.slice4_elasticity, fold_.slice5_elasticity, fold_.slice6_elasticity};
  const double supports[kFoldSlices] = {fold_.slice1_support, fold_.slice2_support, fold_.slice3_support,
                                        fold_.slice4_support, fold_.slice5_support, fold_.slice6_support};
  mass = masses[slice];
  elasticity = elasticities[slice];
  support = supports[slice];
}

double Cord::strip_mass(int strip) const {
  if (strip_count() == 1) {
    return 1.0;
  }
  return strip == 0 ? fold_.strip1_mass : strip == 1 ? fold_.strip2_mass : fold_.strip3_mass;
}

// Each strip's opening at each slice on the wedge-shaped rest gap (zero at
// the anterior commissure, front coordinate z 0..1), interpolated from the
// lower edge's rest at the bottom slice to the upper edge's at the top: the
// half gap rest z + x along the strip's part of the length [z0, z1]. A
// slice's open length, closed share and area are its strips' sums.
void Cord::fold_geometry() {
  int strips = strip_count();
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    double share = static_cast<double>(slice) / static_cast<double>(kFoldSlices - 1);
    double rest = posture_.lower_rest + share * (posture_.upper_rest - posture_.lower_rest) + rest_offset_;
    double slice_open = 0.0;
    double slice_closed = 0.0;
    double slice_area = 0.0;
    for (int strip = 0; strip < strips; ++strip) {
      double z0 = static_cast<double>(strip) / static_cast<double>(strips);
      double z1 = static_cast<double>(strip + 1) / static_cast<double>(strips);
      double strip_length = (z1 - z0) * length_;
      double x = cover_q_[strip][slice];
      double back = rest >= 0.0 ? rest * z1 + x : rest + x;
      double front = rest >= 0.0 ? rest * z0 + x : rest + x;
      double open_length = 0.0;
      double mean_width = 0.0;
      // The contact's centre, measured from the back (vocal process).
      double contact_from_back = 0.5 * ((1.0 - z0) + (1.0 - z1)) * length_;
      if (back > 0.0 && front > 0.0) {
        open_length = strip_length;
        mean_width = back + front;
      } else if (back <= 0.0 && front <= 0.0) {
        open_length = 0.0;
      } else {
        // Open from where the gap closes to the strip's back end.
        double apex = rest > 0.0 ? length_ * (1.0 + x / rest) : 0.0;
        apex = std::clamp(apex, (1.0 - z1) * length_, (1.0 - z0) * length_);
        open_length = apex - (1.0 - z1) * length_;
        mean_width = std::max(back, 0.0);
        contact_from_back = 0.5 * (apex + (1.0 - z0) * length_);
      }
      strip_open_[strip][slice] = open_length;
      strip_alpha_[strip][slice] = (strip_length - open_length) / strip_length;
      strip_star_[strip][slice] = rest >= 0.0 ? rest * (1.0 - contact_from_back / length_) : rest;
      slice_open += open_length;
      slice_closed += strip_length - open_length;
      slice_area += open_length * mean_width;
    }
    slice_open_[slice] = slice_open;
    slice_alpha_[slice] = slice_closed / length_;
    slice_area_[slice] = slice_area + std::max(posture_.chink, 0.0);
  }
}

double Cord::contact() const {
  if (pinch_on_) {
    return pinch_.contact();
  }
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
// (the wedge rest plus the displacement of the strip there) through a soft
// contact, summed over eight points along the length per strip; once per
// step, after the pieces have moved. (The pieces' contact forces use the
// plain geometry.)
void Cord::fold_air_areas() {
  constexpr int kPointsPerStrip = 8;
  int strips = strip_count();
  int points = kPointsPerStrip * strips;
  double softness = std::max(fold_.contact_softness, 1.0e-6);
  double step = length_ / static_cast<double>(points);
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    double share = static_cast<double>(slice) / static_cast<double>(kFoldSlices - 1);
    double rest = posture_.lower_rest + share * (posture_.upper_rest - posture_.lower_rest) + rest_offset_;
    double open_sum = 0.0;
    for (int point = 0; point < points; ++point) {
      double z = (static_cast<double>(point) + 0.5) / static_cast<double>(points);  // 0 front .. 1 vocal process
      double x = cover_q_[point / kPointsPerStrip][slice];
      // With fold_.string_mode on, the displacement follows a string's first
      // mode along the fold, (pi/2) sin(pi z) (fixed at the commissure and
      // the vocal process, the same mean): it closes progressively, not
      // all along its length at once.
      double shape = fold_.string_mode > 0.5 ? kHalfPi * std::sin(kPi * z) : 1.0;
      double half_gap = (rest >= 0.0 ? rest * z : rest) + x * shape;
      double ratio = half_gap / softness;
      double soft = ratio > 30.0 ? half_gap : softness * std::log1p(std::exp(std::max(ratio, -30.0)));
      open_sum += 2.0 * soft * step;
    }
    slice_area_[slice] = open_sum + std::max(posture_.chink, 0.0);
  }
}

// The penetration at which a nonlinear contact equals the linear spring (cm).
constexpr double kContactDepth = 0.01;

// The fold tube's pieces, eight steps per step (1.3 us), symplectic Euler
// (every mode under 2 kHz: omega h < 0.02). Forces from the air's pressures
// of the last step; contact and its damping where the slice is closed. Each
// strip is a column of slices of cover, body and support with its share of
// the length (so of the mass, the stiffness and the air's force) times its
// own factors; neighbouring strips pull on each other's cover and body
// along the fold (fold_.strip_coupling).
void Cord::advance_fold_tube() {
  double q = std::max(tension_, 0.05);
  double slices = static_cast<double>(kFoldSlices);
  int strips = strip_count();
  double slice_depth = fold_.depth / std::sqrt(q) / slices;
  double cover_mass = fold_.cover_mass / (slices * q);
  double body_mass = fold_.body_mass / (slices * q);
  double cover_k = fold_.cover_stiffness * q / slices;
  double coupling_k = fold_.coupling_stiffness * q / slices;
  double support_k = fold_.support_stiffness * q / slices;
  double vertical_k = fold_.vertical_stiffness * q;
  double contact_k = fold_.contact_stiffness * q / slices;
  double along_k = fold_.strip_coupling * q / slices;
  constexpr double h = kStepSeconds / static_cast<double>(kFoldSubsteps);
  double pole = fold_.muscle_corner > 0.0 ? std::exp(-kTwoPi * fold_.muscle_corner * h) : 0.0;
  double innovation = fold_.muscle_fluctuation * std::sqrt(1.0 - pole * pole);
  double eta = fold_.strain_stiffening;
  double eta_vertical = fold_.vertical_stiffening;
  // Per strip: its share of the length, masses and stiffnesses.
  double length_share[kMaxStrips];
  double strip_cover_mass[kMaxStrips];
  double strip_body_mass[kMaxStrips];
  double strip_k[kMaxStrips];
  double strip_cover_c[kMaxStrips];
  double strip_support_c[kMaxStrips];
  double strip_contact_c[kMaxStrips];
  for (int strip = 0; strip < strips; ++strip) {
    length_share[strip] = 1.0 / static_cast<double>(strips);
    double mass = length_share[strip] * strip_mass(strip);
    strip_cover_mass[strip] = cover_mass * mass;
    strip_body_mass[strip] = body_mass * mass;
    strip_k[strip] = length_share[strip] * strip_stiffness(strip);
    strip_cover_c[strip] =
        2.0 * fold_.cover_loss * std::sqrt(strip_cover_mass[strip] * (cover_k + coupling_k) * strip_k[strip]);
    strip_support_c[strip] = 2.0 * fold_.support_loss * std::sqrt(strip_body_mass[strip] * support_k * strip_k[strip]);
    strip_contact_c[strip] = 2.0 * std::sqrt(strip_cover_mass[strip] * contact_k * length_share[strip]);
  }
  // Each slice's own pieces (fold_.slice*_mass, _elasticity, _support) and
  // its wandering elasticity, advanced once per step.
  double slice_mass[kFoldSlices];
  double slice_elastic[kFoldSlices];
  double slice_support[kFoldSlices];
  double wander_pole = std::exp(-kTwoPi * fold_.elasticity_corner * kStepSeconds);
  double wander_innovation = fold_.elasticity_wander * std::sqrt(1.0 - wander_pole * wander_pole);
  for (int slice = 0; slice < kFoldSlices; ++slice) {
    slice_factors(slice, slice_mass[slice], slice_elastic[slice], slice_support[slice]);
    if (fold_.elasticity_wander > 0.0) {
      double tamp = std::exp(-fold_.elasticity_tamping * (slice_support[slice] - 1.0));
      elasticity_drift_[slice] = wander_pole * elasticity_drift_[slice] + wander_innovation * tamp * jitter_gauss();
      slice_elastic[slice] *= std::max(1.0 + elasticity_drift_[slice], 0.05);
    }
  }
  double cover_force[kMaxStrips][kFoldSlices];
  double body_force[kMaxStrips][kFoldSlices];
  for (int substep = 0; substep < kFoldSubsteps; ++substep) {
    if (fold_.muscle_fluctuation > 0.0) {
      muscle_ = pole * muscle_ + innovation * jitter_gauss();
    }
    double muscle_k = support_k * (1.0 + muscle_);
    fold_geometry();
    for (int strip = 0; strip < strips; ++strip) {
      double k_share = strip_k[strip];
      double air_share = length_share[strip];
      for (int slice = 0; slice < kFoldSlices; ++slice) {
        double x = cover_q_[strip][slice];
        double force = slice_pressure_[slice] * strip_open_[strip][slice] * slice_depth;
        if (slice == 0) {
          force += air_share * 0.25 * (pressure_sub_ + slice_pressure_[0]) * kInletLength * length_;
        }
        if (slice == kFoldSlices - 1) {
          force += air_share * 0.25 * (pressure_supra_ + slice_pressure_[slice]) * kOutletLength * length_ +
                   surface_pressure_ * (strip_open_[strip][slice] * slice_depth + air_share * 0.5 * kOutletLength * length_);
        }
        // Tissue stiffens with strain (fold_.strain_stiffening, eta: k (1 +
        // eta x^2), as in Story & Titze's body-cover springs).
        double stretch = x - body_q_[strip][slice];
        double elastic = k_share * slice_elastic[slice];
        force -= elastic * (cover_k * x * (1.0 + eta * x * x) + coupling_k * stretch * (1.0 + eta * stretch * stretch));
        if (slice > 0) {
          double shear = x - cover_q_[strip][slice - 1];
          force -= k_share * vertical_k * shear * (1.0 + eta_vertical * shear * shear);
        }
        if (slice + 1 < kFoldSlices) {
          double shear = x - cover_q_[strip][slice + 1];
          force -= k_share * vertical_k * shear * (1.0 + eta_vertical * shear * shear);
        }
        double body_x = body_q_[strip][slice];
        double body_along = 0.0;
        if (strip > 0) {
          force -= along_k * (x - cover_q_[strip - 1][slice]);
          body_along -= along_k * (body_x - body_q_[strip - 1][slice]);
        }
        if (strip + 1 < strips) {
          force -= along_k * (x - cover_q_[strip + 1][slice]);
          body_along -= along_k * (body_x - body_q_[strip + 1][slice]);
        }
        double alpha = strip_alpha_[strip][slice];
        if (alpha > 0.0) {
          // Tissue pressed on tissue (fold_.contact_exponent: 1 a linear
          // spring; 1.5 Hertzian, soft bodies in contact), equal to the linear
          // spring at kContactDepth of penetration.
          double penetration = x + strip_star_[strip][slice];
          double pressed = fold_.contact_exponent == 1.0
                               ? penetration
                               : (penetration < 0.0 ? -1.0 : 1.0) * kContactDepth *
                                     std::pow(std::abs(penetration) / kContactDepth, fold_.contact_exponent);
          force -= air_share * contact_k * alpha * pressed + strip_contact_c[strip] * alpha * cover_v_[strip][slice];
        }
        double relative = cover_v_[strip][slice] - body_v_[strip][slice];
        // Damping follows each slice's own pieces (critical share fixed).
        double cover_c = strip_cover_c[strip] * std::sqrt(slice_elastic[slice]);
        double support_c = strip_support_c[strip] * std::sqrt(slice_mass[slice] * slice_support[slice]);
        cover_force[strip][slice] = force - cover_c * relative;
        body_force[strip][slice] = -k_share * slice_support[slice] * muscle_k * body_x * (1.0 + eta * body_x * body_x) +
                                   elastic * coupling_k * stretch * (1.0 + eta * stretch * stretch) + body_along -
                                   support_c * body_v_[strip][slice] + cover_c * relative;
      }
    }
    for (int strip = 0; strip < strips; ++strip) {
      for (int slice = 0; slice < kFoldSlices; ++slice) {
        cover_v_[strip][slice] += h * cover_force[strip][slice] / strip_cover_mass[strip];
        body_v_[strip][slice] += h * body_force[strip][slice] / (strip_body_mass[strip] * slice_mass[slice]);
        cover_q_[strip][slice] += h * cover_v_[strip][slice];
        body_q_[strip][slice] += h * body_v_[strip][slice];
      }
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
  // The edges' displacement: the strips' mean, weighted by their length.
  x_[0] = 0.0;
  x_[1] = 0.0;
  for (int strip = 0; strip < strips; ++strip) {
    x_[0] += length_share[strip] * cover_q_[strip][0];
    x_[1] += length_share[strip] * cover_q_[strip][kFoldSlices - 1];
  }
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
  if (pinch_on_) {
    // The posture's rest positions (each side's) are the pinch sections'
    // rest gaps; its chink, the leak beside them.
    pinch_.set_rest(2.0 * posture_.lower_rest, 2.0 * posture_.upper_rest, posture_.chink);
    pinch_.set_pitch(pitch_);
    pinch_.step(wave_below, area_below, wave_above, area_above, out_below, out_above);
    flow_ = pinch_.flow();
    area_ = pinch_.area();
    pressure_sub_ = pinch_.pressure_below();
    pressure_supra_ = pinch_.pressure_above();
    tension_ = pinch_.tension();
    return;
  }
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
    if (vtl_flow_) {
      // VocalTractLab's empirical loss: the larger of a steady-flow law and
      // an acoustic one at 2250 Hz, per cm of each section.
      double sections[2] = {a0, a1};
      viscous = 0.0;
      for (int edge = 0; edge < 2; ++edge) {
        double steady = 0.03 * std::pow(1.0 / sections[edge], 2.6);
        double acoustic = 0.0065 * std::sqrt(2250.0) / sections[edge];
        viscous += thickness_[edge] * std::max(steady, acoustic);
      }
    }
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
  } else if (vtl_flow_) {
    inertance = kAirDensity * (thickness_[0] / std::max(a0, vtl_inertance_floor_) +
                               thickness_[1] / std::max(a1, vtl_inertance_floor_));
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
