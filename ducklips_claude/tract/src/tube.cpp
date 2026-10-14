#include "tract/tube.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace tract {

namespace {

double harmonic_mean(double a, double b) {
  return a + b > 0.0 ? 2.0 * a * b / (a + b) : 0.0;
}

// Frame speed (cm/s) above which a node's pieces are recompiled every step.
constexpr double kFastFrame = 0.5;

// Radius below which the lining is in contact with the opposite wall.
constexpr double kContactRadius = 0.0025;  // cm

}  // namespace

// Viscous loss per cm of tube (VocalTractLab 2.4, Constants.h): the larger of
// a steep narrow-duct law and a fixed-frequency boundary-layer law.
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
  fast_node_.assign(n, 0);
  for (std::vector<double>* field : {&axial_k_, &link_area_, &forward_, &backward_, &forward_next_, &backward_next_, &forward_smooth_, &backward_smooth_,
                                     &series_, &link_resistance_, &link_flow_, &resistance_step_}) {
    field->assign(links, 0.0);
  }
  compile_materials();
  refresh_geometry();
}

void Tube::update_rest_geometry(int ramp_steps) {
  int steps = std::max(ramp_steps, 1);
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    radius_goal_[i] = std::sqrt(std::max(nodes_[i].rest_area, kMinimumArea) / kPi);
    radius_ramp_[i] = (radius_goal_[i] - radius_target_[i]) / steps;
    // Support factors take effect here; masses follow the radius each step.
    compile_node(i);
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
}

void Tube::compile_link(std::size_t i) {
  {
    const SegmentMaterial& left = materials_[static_cast<std::size_t>(nodes_[i].material)];
    const SegmentMaterial& right = materials_[static_cast<std::size_t>(nodes_[i + 1].material)];
    // Edge coupling: the sheet's axial tension acts across the shared ring
    // edge, per unit circumference, over the node spacing.
    double left_k = left.axial_tension * kTwoPi * rest_radius_[i] / kNodeLength;
    double right_k = right.axial_tension * kTwoPi * rest_radius_[i + 1] / kNodeLength;
    axial_k_[i] = harmonic_mean(left_k, right_k);
  }
}

// Current air geometry from the lining's displacement.
void Tube::refresh_geometry() {
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
  for (std::size_t i = 0; i < n; ++i) {
    double offset = radius_target_[i] - rest_radius_[i];
    if (std::abs(offset) < 1.0e-9 && std::abs(radius_rate_[i]) < 1.0e-6) {
      radius_rate_[i] = 0.0;
      fast_node_[i] = 0;
      continue;
    }
    radius_rate_[i] += kStepSeconds * (omega * omega * offset - 2.0 * omega * radius_rate_[i]);
    rest_radius_[i] = std::max(rest_radius_[i] + kStepSeconds * radius_rate_[i], 1.0e-4);
    // A slowly moving node's pieces are recompiled every 8 steps; a fast one
    // (a closure or release) every step. (Every 8 steps for all, its pieces
    // jumped at 12 kHz while the lips or tongue moved fast: a 12 kHz ring,
    // decaying over 15 ms, at every closure and release — the owner heard
    // "a high pitched faint ringing about twice a second".)
    bool fast = std::abs(radius_rate_[i]) > kFastFrame;
    fast_node_[i] = fast;
    if (fast || compile_countdown_ == 0) {
      compile_node(i);
      fast_link_ = fast_link_ || fast;
    }
    moved = true;
  }
  bool recompiled = compile_countdown_ == 0 || fast_link_;
  fast_link_ = false;
  compile_countdown_ = compile_countdown_ == 0 ? 7 : compile_countdown_ - 1;
  if (moved && recompiled) {
    for (std::size_t i = 0; i + 1 < n; ++i) {
      compile_link(i);
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
    frame_low_[i] += frame_gain_ * (surface_[i] * radius_rate_[i] - frame_low_[i]);
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
    double radius = std::max(rest_radius_[i] + cover_q_[i], 0.0);
    area_[i] = std::max(kPi * radius * radius, kMinimumArea);
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    link_area_[i] = 0.5 * (area_[i] + area_[i + 1]);
    if (fast_node_[i] || fast_node_[i + 1]) {
      // Next to a fast-moving node (closure, release): exact every step.
      link_resistance_[i] = viscous_resistance_per_cm(link_area_[i]) * kNodeLength;
      resistance_step_[i] = 0.0;
      continue;
    }
    if (resistance_countdown_ == 0) {
      double target = viscous_resistance_per_cm(link_area_[i]) * kNodeLength;
      resistance_step_[i] = (target - link_resistance_[i]) / 16.0;
    }
    link_resistance_[i] += resistance_step_[i];
  }
  // The viscous law (a power of the area) is evaluated every 16 steps (it
  // was a fifth of the running time) and glides between evaluations: held
  // and then jumped, it pulsed at 6 kHz whenever a constriction moved.
  resistance_countdown_ = resistance_countdown_ == 0 ? 15 : resistance_countdown_ - 1;
}

// The elastic piece: implicit midpoint over the whole lining at once,
//   (M + h^2 K / 4)(q1 + q0) = M (2 q0 + h v0) + h^2 f / 2,
// with K tridiagonal (hoop + coupling on the diagonal, axial tension between
// neighbours) and f the pull of the mass piece. This is the exact sum of all
// of the lining's modes for the step, solved by one tridiagonal sweep.
// The mass piece then follows by backward Euler, held by the support.
void Tube::step_lining() {
  std::size_t n = nodes_.size();
  if (n == 0) {
    return;
  }
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
    double force = coupling_k_[i] * body_q_[i];
    rhs_[i] = cover_mass_[i] * (2.0 * cover_q_[i] + h * cover_v_[i]) + 0.5 * h * h * force;
  }
  // Thomas algorithm (symmetric tridiagonal).
  scratch_[0] = upper_[0] / diag_[0];
  rhs_[0] /= diag_[0];
  for (std::size_t i = 1; i < n; ++i) {
    double denominator = diag_[i] - upper_[i - 1] * scratch_[i - 1];
    scratch_[i] = upper_[i] / denominator;
    rhs_[i] = (rhs_[i] - upper_[i - 1] * rhs_[i - 1]) / denominator;
  }
  for (std::size_t i = n - 1; i > 0; --i) {
    rhs_[i - 1] -= scratch_[i - 1] * rhs_[i];
  }
  for (std::size_t i = 0; i < n; ++i) {
    double old_q = cover_q_[i];
    double new_q = rhs_[i] - old_q;
    cover_v_[i] = 2.0 * (new_q - old_q) / h - cover_v_[i];
    cover_q_[i] = new_q;
  }
  for (std::size_t i = 0; i < n; ++i) {
    double stiffness = coupling_k_[i] + support_k_[i];
    double pull = coupling_k_[i] * cover_q_[i] + support_k_[i] * support_offset_[i] - stiffness * body_q_[i] +
                  cover_c_[i] * cover_v_[i];
    double velocity = (body_mass_[i] * body_v_[i] + h * pull) /
                      (body_mass_[i] + h * (body_c_[i] + cover_c_[i]) + h * h * stiffness);
    body_v_[i] = velocity;
    body_q_[i] += h * velocity;
  }
}

// Air: every node is a junction of its two air paths (and any port). Waves
// arriving this step set the junction pressure together with the lining's
// motion, solved as one scalar implicit exchange per node (no global solve):
// the lining's velocity and the air pressure it sees are found together.
void Tube::step_air() {
  std::size_t n = nodes_.size();
  constexpr double h = kStepSeconds;
  total_wall_flow_ = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double admittance = 0.0;
    double twice_sum = 0.0;
    double from_left = 0.0;
    double from_right = 0.0;
    if (i > 0) {
      double y = link_area_[i - 1] / kRhoC;
      from_left = forward_[i - 1];
      admittance += y;
      twice_sum += 2.0 * y * from_left;
    } else if (left_port_area_ > 0.0) {
      double y = left_port_area_ / kRhoC;
      from_left = left_incoming_;
      admittance += y;
      twice_sum += 2.0 * y * from_left;
    }
    if (i + 1 < n) {
      double y = link_area_[i] / kRhoC;
      from_right = backward_[i];
      admittance += y;
      twice_sum += 2.0 * y * from_right;
    } else if (radiating_right_) {
      // Lip radiation: a resistance R in parallel with an inertance L,
      // discretized by the trapezoid rule for the inertance's flow.
      double a = std::max(area_[i], kMinimumArea);
      double resistance = 128.0 * kRhoC / (9.0 * kPi * kPi * a);
      double inertance = 8.0 * kAirDensity / (3.0 * kPi * std::sqrt(kPi * a));
      double gain = 0.5 * h / inertance;
      admittance += 1.0 / resistance + gain;
      twice_sum -= radiation_inertance_flow_ + gain * radiation_previous_pressure_;
    } else if (right_port_area_ > 0.0) {
      double y = right_port_area_ / kRhoC;
      from_right = right_incoming_;
      admittance += y;
      twice_sum += 2.0 * y * from_right;
    }
    // Side branches joined at this node.
    for (std::size_t port = 0; port < side_node_.size(); ++port) {
      if (side_node_[port] == static_cast<Index>(i) && side_area_[port] > 0.0) {
        double y = side_area_[port] / kRhoC;
        admittance += y;
        twice_sum += 2.0 * y * side_incoming_[port];
      }
    }
    if (admittance <= 0.0) {
      admittance = kMinimumArea / kRhoC;
    }
    // The lining: an outward velocity draws air out of the junction; so
    // does the frame as articulation widens the node.
    double surface = surface_[i];
    twice_sum -= frame_flow_[i];
    double resistance = surface * surface / admittance + cover_c_[i];
    double drive = surface * twice_sum / admittance + cover_c_[i] * body_v_[i];
    double velocity = cover_v_[i];
    double next = ((cover_mass_[i] - 0.5 * h * resistance) * velocity + h * drive) /
                  (cover_mass_[i] + 0.5 * h * resistance);
    double middle = 0.5 * (velocity + next);
    // Contact with the opposite wall: the lining cannot pass its axis.
    if (rest_radius_[i] + cover_q_[i] < kContactRadius && next < 0.0) {
      cover_q_[i] = kContactRadius - rest_radius_[i];
      next = 0.0;
      middle = 0.5 * velocity;
    }
    cover_v_[i] = next;
    double p = (twice_sum - surface * middle) / admittance;
    pressure_[i] = p;
    for (std::size_t port = 0; port < side_node_.size(); ++port) {
      if (side_node_[port] == static_cast<Index>(i)) {
        side_outgoing_[port] = side_area_[port] > 0.0 ? p - side_incoming_[port] : 0.0;
      }
    }
    wall_flow_[i] = -surface * middle;
    total_wall_flow_ += surface * middle;
    if (i > 0) {
      // Each hop passes a gentle one-pole low-pass (wave_smoothing_ k:
      // y = (1 - k) x + k y'): no loss at DC, a loss growing with frequency
      // as the boundary layers' and heat conduction's do, which the
      // frequency-flat viscous term leaves out.
      double wave = p - from_left;
      backward_smooth_[i - 1] = (1.0 - wave_smoothing_) * wave + wave_smoothing_ * backward_smooth_[i - 1];
      backward_next_[i - 1] = backward_smooth_[i - 1];
    } else {
      outgoing_left_ = p - from_left;
    }
    if (i + 1 < n) {
      double wave = p - from_right;
      forward_smooth_[i] = (1.0 - wave_smoothing_) * wave + wave_smoothing_ * forward_smooth_[i];
      forward_next_[i] = forward_smooth_[i];
    } else if (radiating_right_) {
      double a = std::max(area_[i], kMinimumArea);
      double resistance_r = 128.0 * kRhoC / (9.0 * kPi * kPi * a);
      double inertance = 8.0 * kAirDensity / (3.0 * kPi * std::sqrt(kPi * a));
      radiation_inertance_flow_ += 0.5 * h / inertance * (p + radiation_previous_pressure_);
      radiation_previous_pressure_ = p;
      radiated_flow_ = p / resistance_r + radiation_inertance_flow_;
    } else {
      outgoing_right_ = p - from_right;
    }
  }
  std::swap(forward_, forward_next_);
  std::swap(backward_, backward_next_);
  scatter_links();
  std::fill(series_.begin(), series_.end(), 0.0);
}

// Each link carries a series element between its two travelling waves:
// viscous resistance, the kinetic loss of a jet expanding out of a
// narrowing (Borda-Carnot, (rho/2) U^2 (1/A_narrow - 1/A_wide)^2, only in
// the direction the air expands) and any series turbulence source S. With a
// wave a heading right and b heading left on a line of impedance Z,
//   2 (a - b) + S = (2 Z + R) U + K U |U|,
// solved in closed form; the waves then leave as b + Z U (right) and
// a - Z U (left). Pressure therefore builds behind a constriction and the
// flow through it is limited as Bernoulli says. (Losses applied as plain
// wave attenuation, as before, absorb energy but never reflect it: no
// pressure builds behind a closure.)
void Tube::scatter_links() {
  std::size_t n = nodes_.size();
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double to_right = 0.0;
    double to_left = 0.0;
    link_flow_[i] = series_element(forward_[i], backward_[i], kRhoC / link_area_[i], link_resistance_[i], series_[i],
                                   area_[i], area_[i + 1], to_right, to_left);
    forward_[i] = to_right;
    backward_[i] = to_left;
  }
}

void Tube::step() {
  step_lining();
  refresh_geometry();
  step_air();
}

}  // namespace tract
