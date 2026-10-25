#pragma once

// The folds forced by a given returning wave, as sustain FORCEWAVE runs them,
// as a templated copy for forward-mode derivatives.
//
// OBSERVED (src/voice_body.cpp VoiceBody::step, 2026-10-20): with forcing on,
// the cord reads the forced wave in place of the vocal tube's outgoing wave,
// and the only other vocal-tube quantity it reads is vocal_.area(0) (as
// area_above and in forced_flow_). The flow it passes therefore depends on
// the hypothesis only through two per-step inputs, the forced wave and the
// first node's area, plus the lungs, the trachea, the cord and the body's
// hypothesis-independent causes (variation, tremor, lung ramp). The vocal
// tube, nose and branches still step in sustain but nothing they do reaches
// the folds; they are not copied here.
//
// Copied (no edits to src/) for the configuration posture_recover.render
// uses: the triangular glottis (no fold tube, kinematic or pinch mode), no
// aspiration, frication, flux or fold turbulence (their gains are 0, so
// their contributions are exact zeros), pitch feedback 0.08 per closure, no
// ear latency, the default trachea (33 ring-wall nodes, no bronchial tree),
// lung port scale 1. Expressions are kept as the originals write them so the
// double instantiation can be compared bit for bit (the reference is built
// with -O3 -mcpu=native -ffp-contract=fast; so is this).

#include "dual.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

namespace tract_ad {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kStepRate = 96000.0;
constexpr double kStepSeconds = 1.0 / kStepRate;
constexpr long kControlSteps = 12;
constexpr double kControlSeconds = static_cast<double>(kControlSteps) * kStepSeconds;
constexpr double kAirDensity = 1.14e-3;
constexpr double kSoundSpeed = 35000.0;
constexpr double kRhoC = kAirDensity * kSoundSpeed;
constexpr double kAirViscosity = 1.86e-4;
constexpr double kNodeLength = kSoundSpeed / kStepRate;
constexpr double kPascal = 10.0;
constexpr double kMinimumArea = 1.0e-4;

// A checkpoint's double state promoted to S (zero tangents).
template <class S>
S promote_value(double value) {
  return Constant<S>::make(value);
}

template <class S>
void promote_vector(const std::vector<double>& source, std::vector<S>& target) {
  target.resize(source.size());
  for (std::size_t index = 0; index < source.size(); ++index) {
    target[index] = Constant<S>::make(source[index]);
  }
}

// ---- tract::Variation (src/variation.cpp), hypothesis-independent ----

struct VariationState {
  double phase = 0.0;
  double frequency = 0.35;
  double noise = 0.0;
  double target = 0.0;
  double hold = 0.0;
  double saw = 0.0;
  double bubble = 0.0;
  double age = 2.0;
  double duration = 0.3;
  double sign = 1.0;
  double wait = 0.2;
  double value = 0.0;
};

class VariationCopy {
 public:
  VariationCopy(int channels, std::uint32_t seed) : random_(seed != 0U ? seed : 1U) {
    states_.resize(static_cast<std::size_t>(channels));
    for (VariationState& state : states_) {
      state.phase = uniform();
      state.frequency = 0.35 + 0.9 * uniform();
      state.wait = 0.2 + uniform();
    }
  }

  double value(int channel) const { return states_[static_cast<std::size_t>(channel)].value; }

  void advance(double dt) {
    double noise_gain = 1.0 - std::exp(-dt / 0.08);
    double saw_gain = 1.0 - std::exp(-dt / 0.045);
    double value_gain = 1.0 - std::exp(-dt / 0.025);
    for (VariationState& state : states_) {
      state.hold -= dt;
      if (state.hold <= 0.0) {
        state.target = 2.0 * uniform() - 1.0;
        state.hold = 0.08 + 0.32 * uniform();
      }
      state.noise += (state.target - state.noise) * noise_gain;
      state.phase = std::fmod(state.phase + dt * state.frequency, 1.0);
      state.saw += (2.0 * state.phase - 1.0 - state.saw) * saw_gain;
      state.wait -= dt;
      state.age += dt;
      if (state.wait <= 0.0) {
        state.age = 0.0;
        state.duration = 0.16 + 0.4 * uniform();
        state.sign = uniform() < 0.5 ? -1.0 : 1.0;
        state.wait = state.duration + 0.15 + 1.2 * uniform();
      }
      double x = state.age / state.duration;
      state.bubble = x < 1.0 ? state.sign * (256.0 / 27.0) * x * std::pow(1.0 - x, 3.0) : 0.0;
      double mixed = 0.4 * state.noise + 0.25 * state.saw + 0.35 * state.bubble;
      state.value += (mixed - state.value) * value_gain;
    }
  }

 private:
  double uniform() {
    std::uint32_t x = random_;
    x ^= x << 13U;
    x ^= x >> 17U;
    x ^= x << 5U;
    random_ = x;
    return static_cast<double>(x) / 4294967296.0;
  }

  std::uint32_t random_;
  std::vector<VariationState> states_;
};

// ---- tract::series_element (src/tube.cpp) ----

inline double viscous_resistance_per_cm_value(double area) {
  double safe = std::max(area, kMinimumArea);
  return std::max(0.03 * std::pow(safe, -2.6), 0.0065 / safe * std::sqrt(2250.0));
}

template <class S>
S viscous_resistance_per_cm(const S& area) {
  S safe = ad_max(area, kMinimumArea);
  S steep = 0.03 * ad_pow(safe, -2.6);
  S layer = 0.0065 / safe * std::sqrt(2250.0);
  return ad_max(steep, layer);
}

template <class S>
S series_element(const S& forward, const S& backward, const S& z, const S& resistance, double source,
                 const S& first_area, const S& second_area, S& to_second, S& to_first) {
  S drive = 2.0 * (forward - backward) + source;
  S narrow = value_of(drive) >= 0.0 ? first_area : second_area;
  S wide = value_of(drive) >= 0.0 ? second_area : first_area;
  S kinetic = Constant<S>::make(0.0);
  if (value_of(wide) > value_of(narrow)) {
    S difference = 1.0 / ad_max(narrow, kMinimumArea) - 1.0 / wide;
    kinetic = 0.5 * kAirDensity * difference * difference;
  }
  S linear = 2.0 * z + resistance;
  S magnitude = 2.0 * ad_abs(drive) / (linear + ad_sqrt(linear * linear + 4.0 * kinetic * ad_abs(drive)));
  S flow = value_of(drive) >= 0.0 ? magnitude : -magnitude;
  to_second = backward + z * flow;
  to_first = forward - z * flow;
  return flow;
}

// ---- The trachea: tract::Tube (src/tube.cpp) for a static ring-wall tube
// with end ports only (no side ports, radiation, opening, seal, sources or
// hop loss; its own areas follow its linings: no static admittance) ----

template <class S>
class TracheaCopy {
  template <class T>
  friend class TracheaCopy;

 public:
  explicit TracheaCopy(const std::vector<double>& rest_areas) {
    std::size_t n = rest_areas.size();
    std::size_t links = n > 0 ? n - 1 : 0;
    // ring_wall() (src/tube.cpp).
    const double cover_density = 1.04;
    const double cover_thickness = 0.08;
    const double axial_tension = 6000.0;
    const double hoop_tension = 6000.0;
    const double cover_loss = 0.15;
    const double coupling_stiffness = 2.0e9;
    const double body_density = 1.04;
    const double body_depth = 0.5;
    const double support_stiffness = 5.0e6;
    const double support_loss = 1.0;
    rest_radius_.assign(n, 0.0);
    surface_.assign(n, 0.0);
    cover_mass_.assign(n, 0.0);
    body_mass_.assign(n, 0.0);
    hoop_k_.assign(n, 0.0);
    coupling_k_.assign(n, 0.0);
    support_k_.assign(n, 0.0);
    cover_c_.assign(n, 0.0);
    body_c_.assign(n, 0.0);
    support_offset_.assign(n, 0.0);
    frame_flow_.assign(n, 0.0);
    axial_k_.assign(links, 0.0);
    diag_.assign(n, 0.0);
    upper_.assign(n, 0.0);
    scratch_.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      // Tube::compile_materials and compile_node (support factor 1).
      rest_radius_[i] = std::sqrt(std::max(rest_areas[i], kMinimumArea) / kPi);
      double radius = rest_radius_[i];
      double perimeter = kTwoPi * radius;
      double surface = perimeter * kNodeLength;
      surface_[i] = surface;
      cover_mass_[i] = cover_density * cover_thickness * surface;
      body_mass_[i] = body_density * body_depth * surface;
      hoop_k_[i] = hoop_tension * surface / (radius * radius);
      coupling_k_[i] = coupling_stiffness * surface;
      support_k_[i] = support_stiffness * surface * 1.0;
      double cover_local = hoop_k_[i] + coupling_k_[i];
      double reduced = cover_mass_[i] * body_mass_[i] / (cover_mass_[i] + body_mass_[i]);
      cover_c_[i] = 2.0 * cover_loss * std::sqrt(reduced * cover_local);
      body_c_[i] = 2.0 * support_loss * std::sqrt(body_mass_[i] * support_k_[i]);
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
      double left_k = axial_tension * kTwoPi * rest_radius_[i] / kNodeLength;
      double right_k = axial_tension * kTwoPi * rest_radius_[i + 1] / kNodeLength;
      axial_k_[i] = left_k + right_k > 0.0 ? 2.0 * left_k * right_k / (left_k + right_k) : 0.0;
    }
    // The lining's tridiagonal matrix and its factorisation are constant for
    // a static tube (Tube::step_lining recomputes the same values each step).
    constexpr double h = kStepSeconds;
    constexpr double quarter = 0.25 * h * h;
    for (std::size_t i = 0; i < n; ++i) {
      double k_diagonal = hoop_k_[i] + coupling_k_[i];
      if (i > 0) {
        k_diagonal += axial_k_[i - 1];
      }
      if (i + 1 < n) {
        k_diagonal += axial_k_[i];
      }
      diag_[i] = cover_mass_[i] + quarter * k_diagonal;
      upper_[i] = i + 1 < n ? -quarter * axial_k_[i] : 0.0;
    }
    S zero = Constant<S>::make(0.0);
    cover_q_.assign(n, zero);
    cover_v_.assign(n, zero);
    body_q_.assign(n, zero);
    body_v_.assign(n, zero);
    rhs_.assign(n, zero);
    area_.assign(n, zero);
    pressure_.assign(n, zero);
    link_area_.assign(links, zero);
    forward_.assign(links, zero);
    backward_.assign(links, zero);
    forward_next_.assign(links, zero);
    backward_next_.assign(links, zero);
    forward_smooth_.assign(links, zero);
    backward_smooth_.assign(links, zero);
    link_resistance_.assign(links, zero);
    resistance_step_.assign(links, zero);
    link_flow_.assign(links, zero);
    left_incoming_ = zero;
    left_port_area_ = zero;
    right_incoming_ = zero;
    right_port_area_ = zero;
    outgoing_left_ = zero;
    outgoing_right_ = zero;
    refresh_geometry();
  }

  long count() const { return static_cast<long>(rest_radius_.size()); }

  // This tube's state from a double copy's (a checkpoint, zero tangents).
  void promote_from(const TracheaCopy<double>& other) {
    rest_radius_ = other.rest_radius_;
    surface_ = other.surface_;
    cover_mass_ = other.cover_mass_;
    body_mass_ = other.body_mass_;
    hoop_k_ = other.hoop_k_;
    coupling_k_ = other.coupling_k_;
    support_k_ = other.support_k_;
    cover_c_ = other.cover_c_;
    body_c_ = other.body_c_;
    support_offset_ = other.support_offset_;
    frame_flow_ = other.frame_flow_;
    axial_k_ = other.axial_k_;
    diag_ = other.diag_;
    upper_ = other.upper_;
    scratch_ = other.scratch_;
    promote_vector(other.cover_q_, cover_q_);
    promote_vector(other.cover_v_, cover_v_);
    promote_vector(other.body_q_, body_q_);
    promote_vector(other.body_v_, body_v_);
    promote_vector(other.rhs_, rhs_);
    promote_vector(other.area_, area_);
    promote_vector(other.pressure_, pressure_);
    promote_vector(other.link_area_, link_area_);
    promote_vector(other.forward_, forward_);
    promote_vector(other.backward_, backward_);
    promote_vector(other.forward_next_, forward_next_);
    promote_vector(other.backward_next_, backward_next_);
    promote_vector(other.forward_smooth_, forward_smooth_);
    promote_vector(other.backward_smooth_, backward_smooth_);
    promote_vector(other.link_resistance_, link_resistance_);
    promote_vector(other.resistance_step_, resistance_step_);
    promote_vector(other.link_flow_, link_flow_);
    resistance_countdown_ = other.resistance_countdown_;
    left_incoming_ = promote_value<S>(other.left_incoming_);
    left_port_area_ = promote_value<S>(other.left_port_area_);
    right_incoming_ = promote_value<S>(other.right_incoming_);
    right_port_area_ = promote_value<S>(other.right_port_area_);
    outgoing_left_ = promote_value<S>(other.outgoing_left_);
    outgoing_right_ = promote_value<S>(other.outgoing_right_);
  }
  const S& area(long node) const { return area_[static_cast<std::size_t>(node)]; }
  const S& outgoing_left() const { return outgoing_left_; }
  const S& outgoing_right() const { return outgoing_right_; }
  const S& pressure(long node) const { return pressure_[static_cast<std::size_t>(node)]; }

  void set_left_port(const S& incoming, const S& port_area) {
    left_incoming_ = incoming;
    left_port_area_ = port_area;
  }

  void set_right_port(const S& incoming, const S& port_area) {
    right_incoming_ = incoming;
    right_port_area_ = port_area;
  }

  void step() {
    step_lining();
    refresh_geometry();
    step_air();
  }

 private:
  void step_lining() {
    std::size_t n = rest_radius_.size();
    constexpr double h = kStepSeconds;
    for (std::size_t i = 0; i < n; ++i) {
      S force = coupling_k_[i] * body_q_[i];
      rhs_[i] = cover_mass_[i] * (2.0 * cover_q_[i] + h * cover_v_[i]) + 0.5 * h * h * force;
    }
    scratch_[0] = upper_[0] / diag_[0];
    rhs_[0] = rhs_[0] / diag_[0];
    for (std::size_t i = 1; i < n; ++i) {
      double denominator = diag_[i] - upper_[i - 1] * scratch_[i - 1];
      scratch_[i] = upper_[i] / denominator;
      rhs_[i] = (rhs_[i] - upper_[i - 1] * rhs_[i - 1]) / denominator;
    }
    for (std::size_t i = n - 1; i > 0; --i) {
      rhs_[i - 1] = rhs_[i - 1] - scratch_[i - 1] * rhs_[i];
    }
    for (std::size_t i = 0; i < n; ++i) {
      S old_q = cover_q_[i];
      S new_q = rhs_[i] - old_q;
      cover_v_[i] = 2.0 * (new_q - old_q) / h - cover_v_[i];
      cover_q_[i] = new_q;
    }
    for (std::size_t i = 0; i < n; ++i) {
      double stiffness = coupling_k_[i] + support_k_[i];
      S pull = coupling_k_[i] * cover_q_[i] + support_k_[i] * support_offset_[i] - stiffness * body_q_[i] +
               cover_c_[i] * cover_v_[i];
      S velocity = (body_mass_[i] * body_v_[i] + h * pull) /
                   (body_mass_[i] + h * (body_c_[i] + cover_c_[i]) + h * h * stiffness);
      body_v_[i] = velocity;
      body_q_[i] = body_q_[i] + h * velocity;
    }
  }

  // Tube::refresh_geometry for a frame at rest: no ramp, no support glide,
  // no node moving (its rate stays 0, so the frame's pushed air is 0); the
  // air's areas from the lining's displacement; the viscous law glided.
  void refresh_geometry() {
    std::size_t n = rest_radius_.size();
    for (std::size_t i = 0; i < n; ++i) {
      S radius = ad_max(rest_radius_[i] + cover_q_[i], 0.0);
      area_[i] = ad_max(kPi * radius * radius, kMinimumArea);
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
      link_area_[i] = 0.5 * (area_[i] + area_[i + 1]);
      if (resistance_countdown_ == 0) {
        S target = 1.0 * viscous_resistance_per_cm(link_area_[i]) * kNodeLength;
        resistance_step_[i] = (target - link_resistance_[i]) / 16.0;
      }
      link_resistance_[i] = link_resistance_[i] + resistance_step_[i];
    }
    resistance_countdown_ = resistance_countdown_ == 0 ? 15 : resistance_countdown_ - 1;
  }

  void step_air() {
    std::size_t n = rest_radius_.size();
    constexpr double h = kStepSeconds;
    constexpr double kContactRadius = 0.0025;
    constexpr double kContactFrequency = 1500.0;
    for (std::size_t i = 0; i < n; ++i) {
      S admittance = Constant<S>::make(0.0);
      S twice_sum = Constant<S>::make(0.0);
      S from_left = Constant<S>::make(0.0);
      S from_right = Constant<S>::make(0.0);
      if (i > 0) {
        S y = link_area_[i - 1] / kRhoC;
        from_left = forward_[i - 1];
        admittance = admittance + y;
        twice_sum = twice_sum + 2.0 * y * from_left;
      } else if (value_of(left_port_area_) > 0.0) {
        S y = left_port_area_ / kRhoC;
        from_left = left_incoming_;
        admittance = admittance + y;
        twice_sum = twice_sum + 2.0 * y * from_left;
      }
      if (i + 1 < n) {
        S y = link_area_[i] / kRhoC;
        from_right = backward_[i];
        admittance = admittance + y;
        twice_sum = twice_sum + 2.0 * y * from_right;
      } else if (value_of(right_port_area_) > 0.0) {
        S y = right_port_area_ / kRhoC;
        from_right = right_incoming_;
        admittance = admittance + y;
        twice_sum = twice_sum + 2.0 * y * from_right;
      }
      if (value_of(admittance) <= 0.0) {
        admittance = Constant<S>::make(kMinimumArea / kRhoC);
      }
      double surface = surface_[i];
      twice_sum = twice_sum - frame_flow_[i];
      S resistance = surface * surface / admittance + cover_c_[i];
      S drive = surface * twice_sum / admittance + cover_c_[i] * body_v_[i];
      S penetration = kContactRadius - (rest_radius_[i] + cover_q_[i]);
      if (value_of(penetration) > 0.0) {
        double omega = kTwoPi * kContactFrequency;
        drive = drive + cover_mass_[i] * omega * omega * penetration;
        resistance = resistance + 2.0 * cover_mass_[i] * omega;
      }
      S velocity = cover_v_[i];
      S next = ((cover_mass_[i] - 0.5 * h * resistance) * velocity + h * drive) / (cover_mass_[i] + 0.5 * h * resistance);
      S middle = 0.5 * (velocity + next);
      if (value_of(rest_radius_[i] + cover_q_[i]) < -kContactRadius && value_of(next) < 0.0) {
        cover_q_[i] = Constant<S>::make(-kContactRadius - rest_radius_[i]);
        next = Constant<S>::make(0.0);
        middle = 0.5 * velocity;
      }
      cover_v_[i] = next;
      S p = (twice_sum - surface * middle) / admittance;
      pressure_[i] = p;
      if (i > 0) {
        // Hop low-pass with k = wave_smoothing_ = 0 (the trachea's default).
        S wave = p - from_left;
        double k = 0.0;
        backward_smooth_[i - 1] = (1.0 - k) * wave + k * backward_smooth_[i - 1];
        backward_next_[i - 1] = backward_smooth_[i - 1];
      } else {
        outgoing_left_ = p - from_left;
      }
      if (i + 1 < n) {
        S wave = p - from_right;
        double k = 0.0;
        forward_smooth_[i] = (1.0 - k) * wave + k * forward_smooth_[i];
        forward_next_[i] = forward_smooth_[i];
      } else {
        outgoing_right_ = p - from_right;
      }
    }
    forward_.swap(forward_next_);
    backward_.swap(backward_next_);
    for (std::size_t i = 0; i + 1 < n; ++i) {
      S to_right = Constant<S>::make(0.0);
      S to_left = Constant<S>::make(0.0);
      link_flow_[i] = series_element(forward_[i], backward_[i], kRhoC / link_area_[i], link_resistance_[i], 0.0,
                                     area_[i], area_[i + 1], to_right, to_left);
      forward_[i] = to_right;
      backward_[i] = to_left;
    }
  }

  std::vector<double> rest_radius_;
  std::vector<double> surface_;
  std::vector<double> cover_mass_;
  std::vector<double> body_mass_;
  std::vector<double> hoop_k_;
  std::vector<double> coupling_k_;
  std::vector<double> support_k_;
  std::vector<double> cover_c_;
  std::vector<double> body_c_;
  std::vector<double> support_offset_;
  std::vector<double> frame_flow_;
  std::vector<double> axial_k_;
  std::vector<double> diag_;
  std::vector<double> upper_;
  std::vector<double> scratch_;
  std::vector<S> cover_q_;
  std::vector<S> cover_v_;
  std::vector<S> body_q_;
  std::vector<S> body_v_;
  std::vector<S> rhs_;
  std::vector<S> area_;
  std::vector<S> pressure_;
  std::vector<S> link_area_;
  std::vector<S> forward_;
  std::vector<S> backward_;
  std::vector<S> forward_next_;
  std::vector<S> backward_next_;
  std::vector<S> forward_smooth_;
  std::vector<S> backward_smooth_;
  std::vector<S> link_resistance_;
  std::vector<S> resistance_step_;
  std::vector<S> link_flow_;
  int resistance_countdown_ = 0;
  S left_incoming_;
  S left_port_area_;
  S right_incoming_;
  S right_port_area_;
  S outgoing_left_;
  S outgoing_right_;
};

// ---- tract::Cord (src/cord.cpp), the triangular glottis ----

template <class S>
struct PostureCopy {
  S lower_rest;
  S upper_rest;
  S chink;
};

template <class S>
class CordCopy {
  template <class T>
  friend class CordCopy;

 public:
  CordCopy() {
    posture_ = modal();
    guide_ = modal();
    target_ = modal();
    S zero = Constant<S>::make(0.0);
    for (int edge = 0; edge < 2; ++edge) {
      x_[edge] = zero;
      x_previous_[edge] = zero;
      open_length_[edge] = zero;
      mean_width_[edge] = zero;
      alpha_[edge] = zero;
      rest_star_[edge] = zero;
      section_area_[edge] = zero;
    }
    area_ = zero;
    flow_ = zero;
    inertial_pressure_ = zero;
    pressure_sub_ = zero;
    pressure_supra_ = zero;
    pressure_lower_ = zero;
    pressure_upper_ = zero;
  }

  static PostureCopy<S> modal() {
    PostureCopy<S> posture;
    posture.lower_rest = Constant<S>::make(0.00495);
    posture.upper_rest = Constant<S>::make(0.0004);
    posture.chink = Constant<S>::make(0.0002);
    return posture;
  }

  void set_target(const PostureCopy<S>& posture, double move_seconds) {
    target_ = posture;
    double tau = std::max(move_seconds / 3.0, kStepSeconds);
    glide_gain_ = 1.0 - std::exp(-kStepSeconds / (0.5 * tau));
  }

  void set_pitch(double hz) { pitch_ = hz; }

  // This cord's state from a double copy's (a checkpoint, zero tangents).
  void promote_from(const CordCopy<double>& other) {
    PostureCopy<double> postures[3] = {other.posture_, other.guide_, other.target_};
    PostureCopy<S>* targets[3] = {&posture_, &guide_, &target_};
    for (int index = 0; index < 3; ++index) {
      targets[index]->lower_rest = promote_value<S>(postures[index].lower_rest);
      targets[index]->upper_rest = promote_value<S>(postures[index].upper_rest);
      targets[index]->chink = promote_value<S>(postures[index].chink);
    }
    glide_gain_ = other.glide_gain_;
    pitch_ = other.pitch_;
    tension_factor_ = other.tension_factor_;
    rest_offset_ = other.rest_offset_;
    rest_offset_target_ = other.rest_offset_target_;
    rest_offset_guide_ = other.rest_offset_guide_;
    surface_pressure_ = other.surface_pressure_;
    flux_noise_ = other.flux_noise_;
    elongation_ = other.elongation_;
    elongation_reference_ = other.elongation_reference_;
    elongation_floor_ = other.elongation_floor_;
    tension_ = other.tension_;
    tension_guide_ = other.tension_guide_;
    length_ = other.length_;
    inertance_floor_ = other.inertance_floor_;
    for (int edge = 0; edge < 2; ++edge) {
      damping_ratio_[edge] = other.damping_ratio_[edge];
      thickness_[edge] = other.thickness_[edge];
      x_[edge] = promote_value<S>(other.x_[edge]);
      x_previous_[edge] = promote_value<S>(other.x_previous_[edge]);
      open_length_[edge] = promote_value<S>(other.open_length_[edge]);
      mean_width_[edge] = promote_value<S>(other.mean_width_[edge]);
      alpha_[edge] = promote_value<S>(other.alpha_[edge]);
      rest_star_[edge] = promote_value<S>(other.rest_star_[edge]);
      section_area_[edge] = promote_value<S>(other.section_area_[edge]);
    }
    area_ = promote_value<S>(other.area_);
    flow_ = promote_value<S>(other.flow_);
    inertial_pressure_ = promote_value<S>(other.inertial_pressure_);
    pressure_sub_ = promote_value<S>(other.pressure_sub_);
    pressure_supra_ = promote_value<S>(other.pressure_supra_);
    pressure_lower_ = promote_value<S>(other.pressure_lower_);
    pressure_upper_ = promote_value<S>(other.pressure_upper_);
  }

  void set_holding(double tension_factor, double rest_offset) {
    tension_factor_ = tension_factor;
    rest_offset_target_ = rest_offset;
  }

  const S& flow() const { return flow_; }
  double contact() const { return value_of(alpha_[0]); }

  // Diagnosis: the cord's state values (edges, contact, areas, flow, rest
  // gaps, tension) into row[0..11].
  void trace(double* row) const {
    row[0] = value_of(x_[0]);
    row[1] = value_of(x_[1]);
    row[2] = value_of(alpha_[0]);
    row[3] = value_of(alpha_[1]);
    row[4] = value_of(section_area_[0]);
    row[5] = value_of(section_area_[1]);
    row[6] = value_of(flow_);
    row[7] = value_of(posture_.lower_rest) + rest_offset_;
    row[8] = value_of(posture_.upper_rest) + rest_offset_;
    row[9] = tension_;
    row[10] = value_of(pressure_lower_);
    row[11] = value_of(pressure_upper_);
  }

  void step(const S& wave_below, const S& area_below, const S& wave_above, const S& area_above, S& out_below,
            S& out_above) {
    advance_posture();
    advance_mechanics();
    compute_geometry();
    S z_below = kRhoC / ad_max(area_below, kMinimumArea);
    S z_above = kRhoC / ad_max(area_above, kMinimumArea);
    S a0 = ad_max(section_area_[0], kMinimumArea);
    S a1 = ad_max(section_area_[1], kMinimumArea);
    S kinetic = half_rho_over(a0) - half_rho_over(area_below);
    if (value_of(a1) < value_of(a0)) {
      kinetic = kinetic + (half_rho_over(a1) - half_rho_over(a0));
    } else {
      kinetic = kinetic - kExitRecovery * (half_rho_over(a0) - half_rho_over(a1));
    }
    kinetic = kinetic - kExitRecovery * (half_rho_over(a1) - half_rho_over(area_above));
    S viscous = 12.0 * kAirViscosity * length_ * length_ * (thickness_[0] / (a0 * a0 * a0) + thickness_[1] / (a1 * a1 * a1));
    kinetic = ad_max(kinetic, 0.1 * half_rho_over(area_));
    S inertance = kAirDensity * (thickness_[0] / ad_max(a0, inertance_floor_) + thickness_[1] / ad_max(a1, inertance_floor_) +
                                 (kInletLength + kOutletLength) / ad_max(area_, inertance_floor_));
    constexpr double kInertiaDecay = 0.95;
    S mass_term = 2.0 * inertance / kStepSeconds;
    S linear = z_below + z_above + viscous + mass_term;
    S drive = 2.0 * (wave_below - wave_above) + mass_term * flow_ + kInertiaDecay * inertial_pressure_;
    S magnitude = 2.0 * ad_abs(drive) / (linear + ad_sqrt(linear * linear + 4.0 * kinetic * ad_abs(drive)));
    S flow = value_of(drive) >= 0.0 ? magnitude : -magnitude;
    inertial_pressure_ = mass_term * (flow - flow_) - kInertiaDecay * inertial_pressure_;
    flow_ = flow;
    S heard = flow + flux_noise_;
    S p_below = 2.0 * wave_below - z_below * heard;
    S p_above = 2.0 * wave_above + z_above * heard;
    out_below = p_below - wave_below;
    out_above = p_above - wave_above;
    S speed_term = flow * ad_abs(flow);
    pressure_sub_ = p_below;
    pressure_supra_ = p_above;
    pressure_lower_ = p_below - speed_term * (half_rho_over(a0) - half_rho_over(area_below));
    if (value_of(a1) < value_of(a0)) {
      pressure_upper_ = pressure_lower_ - speed_term * (half_rho_over(a1) - half_rho_over(a0));
    } else {
      pressure_upper_ = pressure_lower_ + kExitRecovery * speed_term * (half_rho_over(a0) - half_rho_over(a1));
    }
  }

 private:
  static constexpr double kRestThickness[2] = {0.24, 0.06};
  static constexpr double kMass[2] = {0.12, 0.03};
  static constexpr double kSpring[2] = {80000.0, 8000.0};
  static constexpr double kContactSpring[2] = {240000.0, 24000.0};
  static constexpr double kCoupling = 25000.0;
  static constexpr double kRestLength = 1.3;
  static constexpr double kInletLength = 0.05;
  static constexpr double kOutletLength = 0.01;
  static constexpr double kNaturalPitch = 130.98;
  static constexpr double kPitchPerTension = 132.07;
  static constexpr double kExitRecovery = 0.2;

  static S half_rho_over(const S& area) {
    S a = ad_max(area, kMinimumArea);
    return 0.5 * kAirDensity / (a * a);
  }

  void advance_posture() {
    guide_.lower_rest = guide_.lower_rest + glide_gain_ * (target_.lower_rest - guide_.lower_rest);
    guide_.upper_rest = guide_.upper_rest + glide_gain_ * (target_.upper_rest - guide_.upper_rest);
    guide_.chink = guide_.chink + glide_gain_ * (target_.chink - guide_.chink);
    posture_.lower_rest = posture_.lower_rest + glide_gain_ * (guide_.lower_rest - posture_.lower_rest);
    posture_.upper_rest = posture_.upper_rest + glide_gain_ * (guide_.upper_rest - posture_.upper_rest);
    posture_.chink = posture_.chink + glide_gain_ * (guide_.chink - posture_.chink);
    double offset_gain = 1.0 - std::exp(-kStepSeconds / 0.001);
    rest_offset_guide_ += offset_gain * (rest_offset_target_ - rest_offset_guide_);
    rest_offset_ += offset_gain * (rest_offset_guide_ - rest_offset_);
    // Elongation is off for the triangular glottis (elongation_ 0): the
    // stretch is exactly 1 and the deflection it would read is unused.
    double stretch = std::max((1.0 + elongation_ * 0.0) / (1.0 + elongation_ * elongation_reference_), elongation_floor_);
    double natural = kNaturalPitch;
    double per_tension = kPitchPerTension;
    double anticipated = 1.0;
    double command = std::max(1.0 + (pitch_ - natural) / per_tension, 0.05) * tension_factor_ * stretch / anticipated;
    double tension_gain = 1.0 - std::exp(-kStepSeconds / 0.0015);
    tension_guide_ += tension_gain * (command - tension_guide_);
    tension_ += tension_gain * (tension_guide_ - tension_);
    length_ = kRestLength * std::sqrt(tension_);
    thickness_[0] = kRestThickness[0] / std::sqrt(tension_);
    thickness_[1] = kRestThickness[1] / std::sqrt(tension_);
  }

  void compute_geometry() {
    S rest[2] = {posture_.lower_rest + rest_offset_, posture_.upper_rest + rest_offset_};
    for (int edge = 0; edge < 2; ++edge) {
      S back = rest[edge] + x_[edge];
      S front = value_of(rest[edge]) >= 0.0 ? x_[edge] : back;
      S contact_length = Constant<S>::make(0.0);
      S mean_contact_z = Constant<S>::make(0.5 * length_);
      if (value_of(back) > 0.0 && value_of(front) > 0.0) {
        open_length_[edge] = Constant<S>::make(length_);
        mean_width_[edge] = back + front;
      } else if (value_of(back) <= 0.0 && value_of(front) <= 0.0) {
        open_length_[edge] = Constant<S>::make(0.0);
        mean_width_[edge] = Constant<S>::make(0.0);
        contact_length = Constant<S>::make(length_);
      } else {
        S apex = value_of(rest[edge]) > 0.0 ? length_ * (1.0 + x_[edge] / rest[edge]) : Constant<S>::make(0.0);
        apex = ad_clamp(apex, 0.0, length_);
        open_length_[edge] = apex;
        mean_width_[edge] = ad_max(back, 0.0);
        contact_length = length_ - apex;
        mean_contact_z = 0.5 * (apex + length_);
      }
      alpha_[edge] = contact_length / length_;
      rest_star_[edge] = value_of(rest[edge]) >= 0.0 ? rest[edge] * (1.0 - mean_contact_z / length_) : rest[edge];
      section_area_[edge] = open_length_[edge] * mean_width_[edge] + ad_max(posture_.chink, 0.0);
    }
    area_ = ad_max(ad_min(section_area_[0], section_area_[1]), kMinimumArea);
  }

  void advance_mechanics() {
    constexpr double h = kStepSeconds;
    double mass[2];
    double spring[2];
    double contact[2];
    S damping[2];
    S force[2];
    for (int edge = 0; edge < 2; ++edge) {
      mass[edge] = kMass[edge] / tension_;
      spring[edge] = kSpring[edge] * tension_;
      contact[edge] = kContactSpring[edge] * tension_;
      damping[edge] = 2.0 * (damping_ratio_[edge] + alpha_[edge]) * std::sqrt(kMass[edge] * kSpring[edge]);
    }
    force[0] = pressure_lower_ * open_length_[0] * thickness_[0] +
               0.25 * (pressure_sub_ + pressure_lower_) * kInletLength * length_;
    force[1] = pressure_upper_ * open_length_[1] * thickness_[1] +
               0.25 * (pressure_supra_ + pressure_upper_) * kOutletLength * length_ +
               surface_pressure_ * (open_length_[1] * thickness_[1] + 0.5 * kOutletLength * length_);
    S a = mass[0] + damping[0] * h + h * h * (spring[0] + contact[0] * alpha_[0]) + kCoupling * h * h;
    double b = -kCoupling * h * h;
    S d = mass[1] + damping[1] * h + h * h * (spring[1] + contact[1] * alpha_[1]) + kCoupling * h * h;
    S e = force[0] * h * h + 2.0 * mass[0] * x_[0] - mass[0] * x_previous_[0] + damping[0] * h * x_[0] -
          h * h * contact[0] * alpha_[0] * rest_star_[0];
    S f = force[1] * h * h + 2.0 * mass[1] * x_[1] - mass[1] * x_previous_[1] + damping[1] * h * x_[1] -
          h * h * contact[1] * alpha_[1] * rest_star_[1];
    S determinant = a * d - b * b;
    S next0 = (e * d - b * f) / determinant;
    S next1 = (a * f - b * e) / determinant;
    x_previous_[0] = x_[0];
    x_previous_[1] = x_[1];
    x_[0] = next0;
    x_[1] = next1;
  }

  PostureCopy<S> posture_;
  PostureCopy<S> guide_;
  PostureCopy<S> target_;
  double glide_gain_ = 1.0;
  double pitch_ = 120.0;
  double tension_factor_ = 1.0;
  double rest_offset_ = 0.0;
  double rest_offset_target_ = 0.0;
  double rest_offset_guide_ = 0.0;
  double surface_pressure_ = 0.0;
  double flux_noise_ = 0.0;
  double damping_ratio_[2] = {0.1, 0.6};
  double elongation_ = 0.0;
  double elongation_reference_ = 0.0;
  double elongation_floor_ = 0.0;
  double tension_ = 1.0;
  double tension_guide_ = 1.0;
  double thickness_[2] = {0.24, 0.06};
  double length_ = 1.3;
  double inertance_floor_ = 0.01;
  S x_[2];
  S x_previous_[2];
  S open_length_[2];
  S mean_width_[2];
  S alpha_[2];
  S rest_star_[2];
  S section_area_[2];
  S area_;
  S flow_;
  S inertial_pressure_;
  S pressure_sub_;
  S pressure_supra_;
  S pressure_lower_;
  S pressure_upper_;
};

// ---- tract::VoiceBody::step in forcing mode (src/voice_body.cpp), with
// sustain's settings (tools/sustain.cpp via posture_recover.render) ----

template <class S>
class ForcedFolds {
  template <class T>
  friend class ForcedFolds;

 public:
  // `vocal_count`: the vocal tube's node count (it sets the variation's
  // channel count, whose draws interleave). `pitch`, `lung` (Pa): sustain's
  // arguments.
  ForcedFolds(int vocal_count, double pitch, double lung)
      : trachea_(trachea_areas()), variation_(kFirstWallChannel + vocal_count, kSeed) {
    tremor_scale_ = tremor_normalisation();
    intended_pitch_ = pitch;
    cord_.set_pitch(pitch);
    cord_.set_target(CordCopy<S>::modal(), 0.001);
    set_lung_pressure(lung, 0.08);
    zero_ = Constant<S>::make(0.0);
    from_cord_below_ = zero_;
    from_cord_above_ = zero_;
    forced_flow_ = zero_;
  }

  // The ear's gain per closure (BodyLayout::pitch_feedback, 0.08 as sustain
  // runs it). For diagnosis only: 0 removes the closures' whole-step timing
  // from the folds (see hear_pitch).
  void set_pitch_feedback(double gain) { pitch_feedback_ = gain; }

  // sustain FOLDKEYS: the folds' posture target set at a 1 ms tick (move 1 ms).
  void set_fold_target(const PostureCopy<S>& posture) { cord_.set_target(posture, 0.001); }

  // One step: the folds read `wave` as the vocal tube's returning wave and
  // `area0` as its first node's area (this step's, after the tube's step);
  // returns the flow they pass into the tube (VoiceBody::forced_flow).
  // This body's state from a double copy's (a checkpoint taken at a step:
  // the tangents start there at zero).
  void promote_from(const ForcedFolds<double>& other) {
    trachea_.promote_from(other.trachea_);
    cord_.promote_from(other.cord_);
    variation_ = other.variation_;
    zero_ = promote_value<S>(0.0);
    from_cord_below_ = promote_value<S>(other.from_cord_below_);
    from_cord_above_ = promote_value<S>(other.from_cord_above_);
    forced_flow_ = promote_value<S>(other.forced_flow_);
    lung_flow_slow_ = promote_value<S>(other.lung_flow_slow_);
    tick_ = other.tick_;
    intended_pitch_ = other.intended_pitch_;
    pitch_feedback_ = other.pitch_feedback_;
    pitch_correction_ = other.pitch_correction_;
    was_closed_ = other.was_closed_;
    last_closure_ = other.last_closure_;
    for (int channel = 0; channel < kTremorChannels; ++channel) {
      tremor_slow_[channel] = other.tremor_slow_[channel];
      tremor_slower_[channel] = other.tremor_slower_[channel];
      tremor_fast_[channel] = other.tremor_fast_[channel];
      tremor_faster_[channel] = other.tremor_faster_[channel];
      tremor_value_[channel] = other.tremor_value_[channel];
    }
    tremor_scale_ = other.tremor_scale_;
    tremor_random_ = other.tremor_random_;
    for (int cause = 0; cause < 3; ++cause) {
      cause_from_[cause] = other.cause_from_[cause];
      cause_to_[cause] = other.cause_to_[cause];
    }
    breath_factor_ = other.breath_factor_;
    lung_ = other.lung_;
    lung_target_ = other.lung_target_;
    lung_ramp_ = other.lung_ramp_;
    lung_ramp_from_ = other.lung_ramp_from_;
    lung_ramp_left_ = other.lung_ramp_left_;
    lung_gain_ = other.lung_gain_;
  }

  // Diagnosis: the cord's trace (CordCopy::trace) and the ear's pitch
  // correction into row[0..12].
  void trace(double* row) const {
    cord_.trace(row);
    row[12] = pitch_correction_;
  }

  S step(const S& wave, const S& area0) {
    if (tick_ % kControlSteps == 0) {
      control_tick();
    }
    double share = static_cast<double>(tick_ % kControlSteps) / static_cast<double>(kControlSteps);
    ++tick_;
    if (lung_ramp_left_ > 0) {
      lung_ramp_ += (lung_target_ - lung_ramp_from_) / static_cast<double>(kControlSteps);
      --lung_ramp_left_;
    } else {
      lung_ramp_ = lung_target_;
    }
    lung_ += lung_gain_ * (lung_ramp_ - lung_);
    cord_.set_holding(cause_from_[0] + share * (cause_to_[0] - cause_from_[0]),
                      cause_from_[1] + share * (cause_to_[1] - cause_from_[1]));
    double breath_target = cause_from_[2] + share * (cause_to_[2] - cause_from_[2]);
    breath_factor_ += (1.0 - std::exp(-kStepSeconds / 0.002)) * (breath_target - breath_factor_);
    double breath = lung_ * breath_factor_;
    // (aspiration_sample: aspiration, flux and fold turbulence are 0, so it
    // sets the cord's flux noise and surface pressure to exact zeros.)
    S port_area = trachea_.area(0) * kLungPortScale;
    S impedance = kRhoC / port_area;
    S lung_flow = (0.5 * breath - trachea_.outgoing_left()) / impedance;
    lung_flow_slow_ = lung_flow_slow_ + (1.0 - std::exp(-kTwoPi * 20.0 * kStepSeconds)) * (lung_flow - lung_flow_slow_);
    S compensation = ad_max(impedance - kAirwayResistance, 0.0) * lung_flow_slow_;
    trachea_.set_left_port(0.5 * (breath + compensation), port_area);
    trachea_.set_right_port(from_cord_below_, trachea_.area(trachea_.count() - 1));
    trachea_.step();
    S returning = wave;
    forced_flow_ = (from_cord_above_ - returning) * area0 / kRhoC;
    cord_.step(trachea_.outgoing_right(), trachea_.area(trachea_.count() - 1), returning, area0, from_cord_below_,
               from_cord_above_);
    hear_pitch();
    return forced_flow_;
  }

 private:
  static constexpr int kTremorChannels = 3;
  static constexpr int kFirstWallChannel = 4;
  static constexpr std::uint32_t kSeed = 73021U;
  static constexpr double kAirwayResistance = 1.5;
  static constexpr double kLungPortScale = 1.0;
  // BodyLayout defaults (include/tract/voice_body.hpp).
  static constexpr double kBreathVariation = 0.03;
  static constexpr double kTensionVariation = 0.015;
  static constexpr double kPostureVariation = 0.0015;
  static constexpr double kTremor = 0.005;
  static constexpr double kAdductorTremor = 0.0007;
  static constexpr double kBreathTremor = 0.015;
  static constexpr double kPitchFeedback = 0.08;
  static constexpr double kPitchCorrectionLimit = 0.3;

  static std::vector<double> trachea_areas() {
    // vowel_layout: 12 cm of rings, wider just below the larynx.
    std::vector<double> areas;
    long trachea_count = std::lround(12.0 / kNodeLength);
    for (long node = 0; node < trachea_count; ++node) {
      double from_top = (trachea_count - node - 0.5) * kNodeLength;
      double area = from_top < 1.0 ? 4.0 - from_top : 2.5;
      areas.push_back(std::max(area, 2.5));
    }
    return areas;
  }

  static double tremor_normalisation() {
    double slow = 0.0;
    double slower = 0.0;
    double fast = 0.0;
    double faster = 0.0;
    double slow_gain = 1.0 - std::exp(-kTwoPi * 2.0 * kControlSeconds);
    double fast_gain = 1.0 - std::exp(-kTwoPi * 15.0 * kControlSeconds);
    std::uint32_t random = 0x9E3779B9U;
    double sum = 0.0;
    constexpr int kWarm = 5000;
    constexpr int kTicks = 100000;
    for (int tick = 0; tick < kWarm + kTicks; ++tick) {
      double uniform_sum = 0.0;
      for (int draw = 0; draw < 4; ++draw) {
        random ^= random << 13U;
        random ^= random >> 17U;
        random ^= random << 5U;
        uniform_sum += static_cast<double>(random) / 4294967296.0;
      }
      double white = (uniform_sum - 2.0) * std::sqrt(3.0);
      slow += slow_gain * (white - slow);
      slower += slow_gain * (slow - slower);
      fast += fast_gain * (white - fast);
      faster += fast_gain * (fast - faster);
      if (tick >= kWarm) {
        sum += (faster - slower) * (faster - slower);
      }
    }
    return 1.0 / std::sqrt(sum / kTicks);
  }

  void set_lung_pressure(double pascal, double move_seconds) {
    lung_ramp_from_ = lung_ramp_;
    lung_target_ = pascal * kPascal;
    lung_ramp_left_ = static_cast<int>(kControlSteps);
    double tau = std::max(move_seconds / 3.0, kStepSeconds);
    lung_gain_ = 1.0 - std::exp(-kStepSeconds / tau);
  }

  double tremor_white() {
    double sum = 0.0;
    for (int draw = 0; draw < 4; ++draw) {
      tremor_random_ ^= tremor_random_ << 13U;
      tremor_random_ ^= tremor_random_ >> 17U;
      tremor_random_ ^= tremor_random_ << 5U;
      sum += static_cast<double>(tremor_random_) / 4294967296.0;
    }
    return (sum - 2.0) * std::sqrt(3.0);
  }

  double tremor_band(int channel) {
    double white = tremor_white();
    double slow_gain = 1.0 - std::exp(-kTwoPi * 2.0 * kControlSeconds);
    double fast_gain = 1.0 - std::exp(-kTwoPi * 15.0 * kControlSeconds);
    tremor_slow_[channel] += slow_gain * (white - tremor_slow_[channel]);
    tremor_slower_[channel] += slow_gain * (tremor_slow_[channel] - tremor_slower_[channel]);
    tremor_fast_[channel] += fast_gain * (white - tremor_fast_[channel]);
    tremor_faster_[channel] += fast_gain * (tremor_fast_[channel] - tremor_faster_[channel]);
    return tremor_scale_ * (tremor_faster_[channel] - tremor_slower_[channel]);
  }

  // VoiceBody::control_tick, the parts that reach the folds and lungs.
  void control_tick() {
    variation_.advance(kControlSeconds);
    for (int channel = 0; channel < kTremorChannels; ++channel) {
      tremor_value_[channel] = tremor_band(channel);
    }
    for (int cause = 0; cause < 3; ++cause) {
      cause_from_[cause] = cause_to_[cause];
    }
    cause_to_[0] = std::exp(pitch_correction_) *
                   (1.0 + kTensionVariation * variation_.value(1) + kTremor * tremor_value_[0]);
    cause_to_[1] = kPostureVariation * variation_.value(2) + kAdductorTremor * tremor_value_[1];
    cause_to_[2] = 1.0 + kBreathVariation * variation_.value(0) + kBreathTremor * tremor_value_[2];
  }

  // VoiceBody::hear_pitch (no ear latency): the closure instants trim the
  // tension. Piecewise constant in the hypothesis (whole steps): carried as
  // a plain double, no derivative.
  void hear_pitch() {
    bool closed = cord_.contact() > 0.0;
    if (closed && !was_closed_) {
      if (last_closure_ >= 0) {
        double heard = kStepRate / static_cast<double>(tick_ - last_closure_);
        double ratio = heard / intended_pitch_;
        bool voiced = lung_ > 2500.0 && ratio > 0.707 && ratio < 1.414;
        if (voiced) {
          pitch_correction_ += pitch_feedback_ * std::log(intended_pitch_ / heard);
          pitch_correction_ = std::clamp(pitch_correction_, -kPitchCorrectionLimit, kPitchCorrectionLimit);
        }
      }
      last_closure_ = tick_;
    }
    was_closed_ = closed;
  }

  TracheaCopy<S> trachea_;
  CordCopy<S> cord_;
  VariationCopy variation_;
  S zero_;
  S from_cord_below_;
  S from_cord_above_;
  S forced_flow_;
  S lung_flow_slow_ = Constant<S>::make(0.0);
  long tick_ = 0;
  double intended_pitch_ = 120.0;
  double pitch_feedback_ = kPitchFeedback;
  double pitch_correction_ = 0.0;
  bool was_closed_ = false;
  long last_closure_ = -1;
  double tremor_slow_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_slower_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_fast_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_faster_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_value_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_scale_ = 1.0;
  std::uint32_t tremor_random_ = 0x2545F491U;
  double cause_from_[3] = {1.0, 0.0, 1.0};
  double cause_to_[3] = {1.0, 0.0, 1.0};
  double breath_factor_ = 1.0;
  double lung_ = 0.0;
  double lung_target_ = 0.0;
  double lung_ramp_ = 0.0;
  double lung_ramp_from_ = 0.0;
  int lung_ramp_left_ = 0;
  double lung_gain_ = 1.0;
};

}  // namespace tract_ad
