#include "tract/tube.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace tract {

namespace {

double harmonic_mean(double a, double b) {
  return a + b > 0.0 ? 2.0 * a * b / (a + b) : 0.0;
}

// Frame speed (cm/s) above which a node's pieces are recompiled every step.
constexpr double kFastFrame = 0.5;

constexpr double kSmoothingArea = 3.0;        // cm^2: where the area-scaled hop loss equals the tube's k
constexpr double kMostSmoothing = 0.6;

}  // namespace

// Viscous loss per cm of tube (VocalTractLab 2.4, Constants.h): the larger of
// a steep narrow-duct law and a fixed-frequency boundary-layer law.
// The hop loss of a link when it follows the cross-section
// (Tube::set_smoothing_by_area): boundary-layer friction and heat
// conduction per unit length go as perimeter over area, 2 sqrt(pi / A) for a
// round duct, so the hop's k scales as sqrt(kSmoothingArea / A), equal to the
// tube's k at kSmoothingArea and at most kMostSmoothing.
double Tube::hop_smoothing(double area) const {
  return std::min(wave_smoothing_ * std::sqrt(kSmoothingArea / std::max(area, kMinimumArea)), kMostSmoothing);
}

double viscous_resistance_per_cm(double area) {
  double safe = std::max(area, kMinimumArea);
  return std::max(0.03 * std::pow(safe, -2.6), 0.0065 / safe * std::sqrt(2250.0));
}

//   2 (a - b) + S = (2 Z + R) U + K U |U|, solved in closed form; the waves
// then leave as b + Z U (toward the second side) and a - Z U (toward the first).
double series_element(double forward, double backward, double z, double resistance, double source,
                      double first_area, double second_area, double& to_second, double& to_first) {
  double drive = 2.0 * (forward - backward) + source;
  double narrow = drive >= 0.0 ? first_area : second_area;
  double wide = drive >= 0.0 ? second_area : first_area;
  double kinetic = 0.0;
  if (wide > narrow) {
    double difference = 1.0 / std::max(narrow, kMinimumArea) - 1.0 / wide;
    kinetic = 0.5 * kAirDensity * difference * difference;
  }
  double linear = 2.0 * z + resistance;
  double magnitude = 2.0 * std::abs(drive) / (linear + std::sqrt(linear * linear + 4.0 * kinetic * std::abs(drive)));
  double flow = drive >= 0.0 ? magnitude : -magnitude;
  to_second = backward + z * flow;
  to_first = forward - z * flow;
  return flow;
}

SegmentMaterial soft_wall() {
  SegmentMaterial material;
  material.name = "soft";
  return material;
}

SegmentMaterial firm_wall() {
  SegmentMaterial material;
  material.name = "firm";
  material.axial_tension = 6000.0;
  material.hoop_tension = 6000.0;
  material.coupling_stiffness = 1.0e9;
  material.support_stiffness = 4.0e5;
  return material;
}

SegmentMaterial hard_wall() {
  SegmentMaterial material;
  material.name = "hard";
  material.cover_thickness = 0.08;
  material.axial_tension = 10000.0;
  material.hoop_tension = 10000.0;
  material.coupling_stiffness = 3.0e9;
  material.body_depth = 0.5;
  material.support_stiffness = 5.0e7;
  return material;
}

SegmentMaterial ring_wall() {
  SegmentMaterial material;
  material.name = "ring";
  material.cover_thickness = 0.08;
  material.axial_tension = 6000.0;
  material.hoop_tension = 6000.0;
  material.coupling_stiffness = 2.0e9;
  material.body_depth = 0.5;
  material.support_stiffness = 5.0e6;
  return material;
}

SegmentMaterial lip_wall() {
  SegmentMaterial material;
  material.name = "lip";
  material.cover_thickness = 0.15;
  material.axial_tension = 2000.0;
  material.hoop_tension = 2000.0;
  material.coupling_stiffness = 3.0e8;
  material.body_depth = 1.2;
  material.support_stiffness = 8.0e4;
  return material;
}

Tube::Tube(std::vector<SegmentMaterial> materials, std::vector<TubeNode> nodes)
    : materials_(std::move(materials)), nodes_(std::move(nodes)) {
  std::size_t n = nodes_.size();
  std::size_t links = n > 0 ? n - 1 : 0;
  for (std::vector<double>* field : {&rest_radius_, &radius_target_, &radius_goal_, &radius_ramp_, &radius_rate_, &frame_flow_, &frame_low_, &frame_lower_, &surface_, &cover_mass_, &body_mass_, &hoop_k_, &coupling_k_,
                                     &support_k_, &cover_c_, &body_c_, &support_offset_, &cover_q_, &cover_v_,
                                     &body_q_, &body_v_, &diag_, &upper_, &rhs_, &scratch_, &area_, &pressure_,
                                     &wall_flow_}) {
    field->assign(n, 0.0);
  }
  support_factor_.assign(n, 1.0);
  support_target_.assign(n, 1.0);
  fast_node_.assign(n, 0);
  for (std::vector<double>* field : {&axial_k_, &link_area_, &forward_, &backward_, &forward_next_, &backward_next_, &forward_smooth_, &backward_smooth_, &forward_raw_, &backward_raw_,
                                     &series_, &link_resistance_, &link_flow_, &resistance_step_, &hop_k_}) {
    field->assign(links, 0.0);
  }
  viscous_area_.assign(links, std::numeric_limits<double>::quiet_NaN());
  viscous_value_.assign(links, 0.0);
  denominator_.assign(n, 0.0);
  body_denominator_.assign(n, 0.0);
  for (std::vector<double>* field : {&link_y_, &link_z_, &kinetic_right_, &kinetic_left_}) {
    field->assign(links, 0.0);
  }
  kinetic_first_area_.assign(links, std::numeric_limits<double>::quiet_NaN());
  kinetic_second_area_.assign(links, std::numeric_limits<double>::quiet_NaN());
  junction_admittance_.assign(n, std::numeric_limits<double>::quiet_NaN());
  for (std::vector<double>* field : {&junction_resistance_, &junction_plus_, &junction_minus_}) {
    field->assign(n, 0.0);
  }
  junction_dirty_.assign(n, 1);
  lining_node_dirty_.assign(n, 0);
  moved_node_.assign(n, 0);
  // (link_area_ starts at 0: every link's quantities are taken at the first
  // geometry refresh.)
  active_count_ = static_cast<Index>(n);
  compile_materials();
  refresh_geometry();
}

void Tube::set_opening(double position, double area) {
  if (position < 0.0) {
    opening_node_ = -1;
    opening_position_ = -1.0;
    opening_position_target_ = -1.0;
    return;
  }
  if (opening_position_ < 0.0) {
    opening_position_ = position;
    opening_area_ = std::max(area, kMinimumArea);
  }
  opening_position_target_ = position;
  opening_area_target_ = std::max(area, kMinimumArea);
}

// The load's node pair from its glided place; its state moves along with
// it when the pair moves by one node.
void Tube::follow_opening() {
  const double glide = 1.0 - std::exp(-kStepSeconds / 0.0005);
  opening_position_ += glide * (opening_position_target_ - opening_position_);
  opening_area_ += glide * (opening_area_target_ - opening_area_);
  Index last = count() - 2;
  Index node = std::clamp(static_cast<Index>(std::floor(opening_position_)), static_cast<Index>(0), last);
  if (node != opening_node_ && opening_node_ >= 0) {
    if (node == opening_node_ + 1) {
      opening_flow_[0] = opening_flow_[1];
      opening_pressure_[0] = opening_pressure_[1];
      opening_flow_[1] = 0.0;
      opening_pressure_[1] = pressure_[static_cast<std::size_t>(std::min(node + 1, count() - 1))];
    } else if (node == opening_node_ - 1) {
      opening_flow_[1] = opening_flow_[0];
      opening_pressure_[1] = opening_pressure_[0];
      opening_flow_[0] = 0.0;
      opening_pressure_[0] = pressure_[static_cast<std::size_t>(node)];
    } else {
      opening_flow_[0] = 0.0;
      opening_flow_[1] = 0.0;
    }
  }
  opening_node_ = node;
  opening_share_ = std::clamp(opening_position_ - static_cast<double>(node), 0.0, 1.0);
}

void Tube::set_active_length(double length_cm) {
  Index total = static_cast<Index>(nodes_.size());
  double nodes = std::max(length_cm, 2.0 * kNodeLength) / kNodeLength;
  Index whole = std::clamp(static_cast<Index>(std::floor(nodes)), static_cast<Index>(2), total);
  // Hysteresis: lengthen only once a quarter node past the next whole node,
  // so a length dithering at a boundary does not toggle a node.
  // Nor while the open end is shut (lips sealed): a node joined behind a
  // closure carried its huge impedance into the new waves (a click in the
  // /p/ of "jumps"); the length waits on the end's air mass meanwhile.
  bool end_open = area_[static_cast<std::size_t>(active_count_ - 1)] > 0.05;
  bool lengthen = nodes >= static_cast<double>(active_count_) + 1.25 && end_open;
  if (whole < active_count_ || lengthen) {
    Index previous = active_count_;
    active_count_ = whole;
    // Links newly joined carry the old open end's state on: waves that give
    // the old end's pressure and the flow it was radiating, and the new nodes
    // start at that pressure. (Joined silent, every lengthening dropped the
    // end's pressure for a step: a click 10-21 dB over the speech, which the
    // output normalisation then took as the peak.)
    if (previous >= 1 && active_count_ > previous) {
      std::size_t end = static_cast<std::size_t>(previous - 1);
      double pressure = pressure_[end];
      double flow = radiated_flow_;
      for (Index link = previous - 1; link + 1 < active_count_; ++link) {
        std::size_t index = static_cast<std::size_t>(link);
        double z = kRhoC / std::max(link_area_[index], kMinimumArea);
        forward_[index] = 0.5 * (pressure + z * flow);
        backward_[index] = 0.5 * (pressure - z * flow);
        forward_smooth_[index] = forward_[index];
        backward_smooth_[index] = backward_[index];
        pressure_[index + 1] = pressure;
      }
    }
  }
  end_extra_ = std::clamp(nodes - static_cast<double>(active_count_), 0.0, 2.5) * kNodeLength;
}

void Tube::update_rest_geometry(int ramp_steps) {
  int steps = std::max(ramp_steps, 1);
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    radius_goal_[i] = std::sqrt(std::max(nodes_[i].rest_area, kMinimumArea) / kPi);
    radius_ramp_[i] = (radius_goal_[i] - radius_target_[i]) / steps;
    // (Pieces are recompiled in refresh_geometry, every step while the node
    // or its support moves. Recompiled here, once per control tick, and
    // held between, every moving wall's masses and stiffnesses changed in a
    // staircase: lines at 8 and 16 kHz, 12 dB over the spectrum around them
    // in 40% of the frames, the owner's "faint ringing".)
  }
  ramp_left_ = steps;
}

void Tube::settle_geometry() {
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    radius_target_[i] = std::sqrt(std::max(nodes_[i].rest_area, kMinimumArea) / kPi);
    radius_goal_[i] = radius_target_[i];
    rest_radius_[i] = radius_target_[i];
    radius_rate_[i] = 0.0;
  }
  ramp_left_ = 0;
  compile_materials();
  refresh_geometry();
}

void Tube::set_rest_area(Index node, double area) {
  nodes_[static_cast<std::size_t>(node)].rest_area = std::max(area, 0.0);
}

void Tube::set_support_offset(Index node, double offset_cm) {
  support_offset_[static_cast<std::size_t>(node)] = offset_cm;
}

void Tube::set_left_port(double incoming, double port_area) {
  left_incoming_ = incoming;
  left_port_area_ = port_area;
}

void Tube::set_right_port(double incoming, double port_area) {
  right_incoming_ = incoming;
  right_port_area_ = port_area;
}

int Tube::add_side_port(Index node) {
  side_node_.push_back(node);
  side_incoming_.push_back(0.0);
  side_area_.push_back(0.0);
  side_outgoing_.push_back(0.0);
  return static_cast<int>(side_node_.size()) - 1;
}

void Tube::set_side_port(int port, double incoming, double port_area) {
  side_incoming_[static_cast<std::size_t>(port)] = incoming;
  side_area_[static_cast<std::size_t>(port)] = std::max(port_area, 0.0);
}

void Tube::add_series_source(Index node, double source) {
  series_[static_cast<std::size_t>(node)] += source;
}

// Masses, stiffnesses and damping of the three pieces at each node. They
// scale with the lining area of the node, which follows the rest geometry.
void Tube::compile_materials() {
  std::size_t n = nodes_.size();
  for (std::size_t i = 0; i < n; ++i) {
    radius_target_[i] = std::sqrt(std::max(nodes_[i].rest_area, kMinimumArea) / kPi);
    radius_goal_[i] = radius_target_[i];
    rest_radius_[i] = radius_target_[i];
    compile_node(i);
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    compile_link(i);
  }
}

// One node's pieces at its current rest radius.
void Tube::compile_node(std::size_t i) {
  {
    const SegmentMaterial& material = materials_[static_cast<std::size_t>(nodes_[i].material)];
    double radius = rest_radius_[i];
    double perimeter = kTwoPi * radius;
    double surface = perimeter * kNodeLength;
    surface_[i] = surface;
    cover_mass_[i] = material.cover_density * material.cover_thickness * surface;
    body_mass_[i] = material.body_density * material.body_depth * surface;
    // A tensioned ring of radius r resists radial displacement with a
    // pressure T / r^2 per unit displacement.
    hoop_k_[i] = material.hoop_tension * surface / (radius * radius);
    coupling_k_[i] = material.coupling_stiffness * surface;
    support_k_[i] = material.support_stiffness * surface * support_factor_[i];
    // Internal friction between the lining and the body acts on their
    // relative motion only; the support alone damps the body against the
    // frame. (Applied to absolute velocities, the friction braked the whole
    // wall at about 3x VocalTractLab's measured wall resistance and widened
    // the closed-tract F1 of /a/ to 360 Hz; rigid walls give 60.)
    double cover_local = hoop_k_[i] + coupling_k_[i];
    double reduced = cover_mass_[i] * body_mass_[i] / (cover_mass_[i] + body_mass_[i]);
    cover_c_[i] = 2.0 * material.cover_loss * std::sqrt(reduced * cover_local);
    body_c_[i] = 2.0 * material.support_loss * std::sqrt(body_mass_[i] * support_k_[i]);
  }
  lining_node_dirty_[i] = 1;
  lining_partial_ = true;
  junction_dirty_[i] = 1;
}

void Tube::scale_axial(double scale) {
  for (SegmentMaterial& material : materials_) {
    material.axial_tension *= scale;
  }
  for (std::size_t i = 0; i + 1 < nodes_.size(); ++i) {
    compile_link(i);
  }
}

void Tube::compile_link(std::size_t i) {
  {
    const SegmentMaterial& left = materials_[static_cast<std::size_t>(nodes_[i].material)];
    const SegmentMaterial& right = materials_[static_cast<std::size_t>(nodes_[i + 1].material)];
    // Edge coupling: the sheet's axial tension acts across the shared ring
    // edge, per unit circumference, over the node spacing.
    double left_k = left.axial_tension * kTwoPi * rest_radius_[i] / kNodeLength;
    double right_k = right.axial_tension * kTwoPi * rest_radius_[i + 1] / kNodeLength;
    double axial = harmonic_mean(left_k, right_k);
    if (!(axial == axial_k_[i])) {
      axial_k_[i] = axial;
      lining_dirty_ = true;
    }
  }
}

// The viscous law at link i's area (viscous_resistance_per_cm), taken
// again only when the area changed.
double Tube::link_viscous(std::size_t i) {
  if (!(link_area_[i] == viscous_area_[i])) {
    viscous_area_[i] = link_area_[i];
    viscous_value_[i] = viscous_resistance_per_cm(link_area_[i]);
  }
  return viscous_value_[i];
}

// Current air geometry from the lining's displacement.
void Tube::refresh_geometry() {
  if (air_frozen_ && static_admittance_) {
    refresh_first_area();   // (tube_step.cpp)
    return;
  }
  if (frame_settled()) {
    // (The frame, its supports and its pushed air at rest: only the areas
    // and the links follow the linings; tube_step.cpp.)
    refresh_settled_geometry();
    return;
  }
  std::size_t n = nodes_.size();
  // The frame follows its target critically damped, corner 160 Hz (about
  // 2 ms lag): quick enough for a stop release to pop. (At 60 Hz, tried
  // while hunting a 1 kHz buzz that came from the velum port instead, the
  // lips parted over 5 ms and the bursts were 10 dB weak.)
  constexpr double omega = kTwoPi * 160.0;
  bool moved = false;
  if (ramp_left_ > 0) {
    for (std::size_t i = 0; i < n; ++i) {
      radius_target_[i] += radius_ramp_[i];
    }
    if (--ramp_left_ == 0) {
      radius_target_ = radius_goal_;
    }
  }
  const double support_glide = 1.0 - std::exp(-kStepSeconds / 0.002);
  for (std::size_t i = 0; i < n; ++i) {
    double factor_gap = support_target_[i] - support_factor_[i];
    bool support_moving = std::abs(factor_gap) > 1.0e-7;
    if (support_moving) {
      if (support_log_glide_ > 0.0) {
        double log_glide = 1.0 - std::exp(-kStepSeconds / support_log_glide_);
        support_factor_[i] = std::exp(std::log(std::max(support_factor_[i], 1.0e-6)) +
                                      log_glide * (std::log(std::max(support_target_[i], 1.0e-6)) -
                                                   std::log(std::max(support_factor_[i], 1.0e-6))));
      } else {
        support_factor_[i] += support_glide * factor_gap;
      }
    }
    double offset = radius_target_[i] - rest_radius_[i];
    bool radius_moving = std::abs(offset) >= 1.0e-9 || std::abs(radius_rate_[i]) >= 1.0e-6;
    if (radius_moving) {
      radius_rate_[i] += kStepSeconds * (omega * omega * offset - 2.0 * omega * radius_rate_[i]);
      rest_radius_[i] = std::max(rest_radius_[i] + kStepSeconds * radius_rate_[i], 1.0e-4);
      fast_node_[i] = std::abs(radius_rate_[i]) > kFastFrame;
      moved = true;
      moved_node_[i] = 1;
    } else {
      radius_rate_[i] = 0.0;
      fast_node_[i] = 0;
      moved_node_[i] = 0;
    }
    // A moving node's pieces follow it every step. (Every 8 steps they
    // jumped at 12 kHz while the lips or tongue moved: a 12 kHz ring at
    // every closure and release, "a high pitched faint ringing about twice a
    // second"; once per control tick, lines at 8 and 16 kHz.)
    if (radius_moving || support_moving) {
      compile_node(i);
    }
  }
  if (moved) {
    // (A link's coupling depends on its two nodes' radii alone: only the
    // links beside a node that moved are compiled again.)
    for (std::size_t i = 0; i + 1 < n; ++i) {
      if (moved_node_[i] || moved_node_[i + 1]) {
        compile_link(i);
      }
    }
  }
  // Air the frame displaces; over a conserved span its net is removed in
  // proportion to each node's lining (the incompressible tongue moves air
  // along the tract; it does not make or destroy volume).
  double net = 0.0;
  double lining = 0.0;
  // The frame's motion is slow; it reaches the air through two poles at
  // 250 Hz. (Its rate follows a target that ramps anew at every 1 kHz
  // command, and the kinks pumped tones at 1, 2 and 3 kHz, 31-37 dB over
  // their neighbourhood, out of a sealed mouth through the nose while the
  // tongue moved during "my"'s /m/.)
  for (std::size_t i = 0; i < n; ++i) {
    // (Nodes past quiet_from_ are outside air: their changing area moves
    // no wall and pushes no air. Masked before the two poles, so a node
    // turning quiet or not fades in or out; switched after them, the step
    // in the pushed air came out as a two-sample spike 20 times the /f/
    // around it.)
    double pushed = static_cast<Index>(i) >= quiet_from_ ? 0.0 : surface_[i] * radius_rate_[i];
    frame_low_[i] += frame_gain_ * (pushed - frame_low_[i]);
    frame_lower_[i] += frame_gain_ * (frame_low_[i] - frame_lower_[i]);
    frame_flow_[i] = frame_lower_[i];
    if (static_cast<Index>(i) >= conserved_first_ && static_cast<Index>(i) <= conserved_last_) {
      net += frame_flow_[i];
      lining += surface_[i];
    }
  }
  if (lining > 0.0) {
    for (std::size_t i = 0; i < n; ++i) {
      if (static_cast<Index>(i) >= conserved_first_ && static_cast<Index>(i) <= conserved_last_) {
        frame_flow_[i] -= net * surface_[i] / lining;
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    double radius = std::max(rest_radius_[i] + (static_admittance_ ? 0.0 : cover_q_[i]), 0.0);
    area_[i] = std::max(kPi * radius * radius, kMinimumArea);
  }
  if (area_limit_) {
    constexpr double kOpeningCurve = 480.0e6;  // cm^2 / s^4
    constexpr double kClosingCurve = 40.0e6;
    constexpr double kStaticLimit = 400.0 * kStepSeconds;  // cm^2 per step
    if (limited_area_.size() != n) {
      limited_area_ = area_;
    }
    for (std::size_t i = 0; i < n; ++i) {
      double old = limited_area_[i];
      double delta = std::clamp(area_[i] - old, -kStaticLimit, kStaticLimit);
      double curve = area_[i] > old ? kOpeningCurve : kClosingCurve;
      double t = std::sqrt(std::sqrt(old / curve)) + kStepSeconds;
      double reach = curve * t * t * t * t - old;
      delta = std::clamp(delta, -reach, reach);
      area_[i] = std::max(old + delta, kMinimumArea);
      limited_area_[i] = area_[i];
    }
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double link_area = 0.5 * (area_[i] + area_[i + 1]);
    if (!(link_area == link_area_[i])) {
      link_area_[i] = link_area;
      link_y_[i] = link_area / kRhoC;
      link_z_[i] = kRhoC / link_area;
    }
    if (smoothing_by_area_ && (resistance_countdown_ == 0 || fast_node_[i] || fast_node_[i + 1])) {
      hop_k_[i] = hop_smoothing(link_area_[i]);
    }
    if (fast_node_[i] || fast_node_[i + 1]) {
      // Next to a fast-moving node (closure, release): exact every step.
      link_resistance_[i] = viscous_scale_ * link_viscous(i) * kNodeLength;
      resistance_step_[i] = 0.0;
      continue;
    }
    if (resistance_countdown_ == 0) {
      double target = viscous_scale_ * link_viscous(i) * kNodeLength;
      resistance_step_[i] = (target - link_resistance_[i]) / 16.0;
    }
    link_resistance_[i] += resistance_step_[i];
  }
  // The viscous law (a power of the area) is evaluated every 16 steps (it
  // was a fifth of the running time) and glides between evaluations: held
  // and then jumped, it pulsed at 6 kHz whenever a constriction moved.
  resistance_countdown_ = resistance_countdown_ == 0 ? 15 : resistance_countdown_ - 1;
}

void Tube::step() {
  if (air_frozen_) {
    refresh_geometry();
    return;
  }
  step_lining();
  refresh_geometry();
  step_air();
}

}  // namespace tract
