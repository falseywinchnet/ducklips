// The tube's step: the lining's implicit midpoint step, the air's
// junctions and the links' series elements (Tube::step_lining,
// factor_lining, step_air, scatter_links).
//
// Pinned arithmetic. This file is compiled with -ffp-contract=off
// (CMakeLists.txt): the compiler fuses no multiply-add here, and each one
// the reference build fused (-ffp-contract=fast, read from its object code,
// 2026-10-20) is written as std::fma. So the roundings are the reference's
// by construction, not by the compiler's choice: under -ffp-contract=fast
// the choice of what to fuse moved with unrelated edits (keeping one
// quantity of step_air in a cache left the reference's paired vector divide
// unpaired, the compiler fused the drive's add, and every render moved by
// about 1e-9 relative; MEASURED, NOTES_engine.md). A "pragma clang fp
// contract(off)" does not stop the backend under -ffp-contract=fast.
// Quantities of the geometry alone are kept between steps while their
// inputs hold (equal by ==): the same inputs give the same values.

#include "tract/tube.hpp"

#include "tract/lanes.hpp"   // <arm_neon.h> on ARM; the same names on SSE2 elsewhere

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstring>
#include <utility>

namespace tract {

namespace {

// Radius below which the lining is in contact with the opposite wall.
constexpr double kContactRadius = 0.0025;  // cm
constexpr double kContactFrequency = 1500.0;  // Hz, the contact spring on the lining's mass

}  // namespace

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
  constexpr double half_h_squared = 0.5 * h * h;
  if (lining_dirty_) {
    factor_lining();
  } else if (lining_partial_) {
    refactor_lining();
  }
  // (The arrays as local pointers: each is its own vector, so none aliases
  // another, and the loops keep them in registers.)
  const double* __restrict coupling_k = coupling_k_.data();
  const double* __restrict support_k = support_k_.data();
  const double* __restrict support_offset = support_offset_.data();
  const double* __restrict cover_mass = cover_mass_.data();
  const double* __restrict cover_c = cover_c_.data();
  const double* __restrict body_mass = body_mass_.data();
  const double* __restrict body_denominator = body_denominator_.data();
  double* __restrict rhs = rhs_.data();
  double* __restrict cover_q = cover_q_.data();
  double* __restrict cover_v = cover_v_.data();
  double* __restrict body_q = body_q_.data();
  double* __restrict body_v = body_v_.data();
  for (std::size_t i = 0; i < n; ++i) {
    // cover_mass (2 q + h v) + (h^2 / 2) force, force = coupling body_q.
    double force = coupling_k[i] * body_q[i];
    double momentum = std::fma(cover_v[i], h, 2.0 * cover_q[i]);
    rhs[i] = std::fma(force, half_h_squared, cover_mass[i] * momentum);
  }
  if (diagonal_ && solve_diagonal(n)) {
    // (Solved without the elimination's chain: see solve_diagonal.)
  } else {
    solve_tridiagonal(n);
  }
  for (std::size_t i = 0; i < n; ++i) {
    double old_q = cover_q[i];
    double new_q = rhs[i] - old_q;
    cover_v[i] = 2.0 * (new_q - old_q) / h - cover_v[i];
    cover_q[i] = new_q;
  }
  for (std::size_t i = 0; i < n; ++i) {
    // pull = coupling cover_q + support offset - stiffness body_q + cover_c cover_v.
    double stiffness = coupling_k[i] + support_k[i];
    double pull = std::fma(coupling_k[i], cover_q[i], support_k[i] * support_offset[i]);
    pull = std::fma(-stiffness, body_q[i], pull);
    pull = std::fma(cover_c[i], cover_v[i], pull);
    double velocity = std::fma(pull, h, body_mass[i] * body_v[i]) / body_denominator[i];
    body_v[i] = velocity;
    body_q[i] = std::fma(velocity, h, body_q[i]);
  }
}

// Whether refresh_geometry's frame part would change nothing this step: no
// ramp and no area limit, and at every node the frame at rest (its offset
// under the 1e-9 cm the follower sleeps below and its rate exactly 0, so
// the follower writes back the zeros it holds), the support's glide asleep
// (gap at most 1e-7), the pushed air's two poles and flow exactly 0, the
// fast and moved marks clear. Then the support, radius and frame loops leave
// every value as it is (0 + g (0 - 0) = 0), no node or link is compiled, and
// refresh_settled_geometry does the rest. (The trachea's frame never moves:
// MEASURED in NOTES_residual_speed.md.)
bool Tube::frame_settled() const {
  if (ramp_left_ != 0 || area_limit_) {
    return false;
  }
  std::size_t n = nodes_.size();
  for (std::size_t i = 0; i < n; ++i) {
    bool at_rest = std::abs(radius_target_[i] - rest_radius_[i]) < 1.0e-9 && radius_rate_[i] == 0.0;
    bool supported = std::abs(support_target_[i] - support_factor_[i]) <= 1.0e-7;
    bool still = frame_low_[i] == 0.0 && frame_lower_[i] == 0.0 && frame_flow_[i] == 0.0;
    bool unmarked = fast_node_[i] == 0 && moved_node_[i] == 0;
    if (!(at_rest && supported && still && unmarked)) {
      return false;
    }
  }
  return true;
}

// refresh_geometry when frame_settled: its area loop and link loop as they
// are there (no fast node: the viscous glide and the area-scaled hop loss
// follow the countdown alone). The fused build computes the glide's step as
// one multiply-subtract, (scale law) kNodeLength - R (fnmsub, read from its
// object code 2026-10-25): written so here; no other product is added
// (bindings/engine_identity.py: the renders identical).
void Tube::refresh_settled_geometry() {
  std::size_t n = nodes_.size();
  for (std::size_t i = 0; i < n; ++i) {
    double radius = std::max(rest_radius_[i] + (static_admittance_ ? 0.0 : cover_q_[i]), 0.0);
    area_[i] = std::max(kPi * radius * radius, kMinimumArea);
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double link_area = 0.5 * (area_[i] + area_[i + 1]);
    if (!(link_area == link_area_[i])) {
      link_area_[i] = link_area;
      link_y_[i] = link_area / kRhoC;
      link_z_[i] = kRhoC / link_area;
    }
    if (smoothing_by_area_ && resistance_countdown_ == 0) {
      hop_k_[i] = hop_smoothing(link_area_[i]);
    }
    if (resistance_countdown_ == 0) {
      double scaled_law = viscous_scale_ * link_viscous(i);
      resistance_step_[i] = std::fma(scaled_law, kNodeLength, -link_resistance_[i]) / 16.0;
    }
    link_resistance_[i] += resistance_step_[i];
  }
  resistance_countdown_ = resistance_countdown_ == 0 ? 15 : resistance_countdown_ - 1;
}

// refresh_geometry for a frozen tube with static admittance
// (set_air_frozen): the ramp of every node's target (its count is shared)
// and the first node's frame, whose area is all a caller then reads (the
// first area depends on no other node's frame). The fused build computes
// the follower as fma(fma(offset, omega^2, rate (-2 omega)), h, rate) and
// the radius as fma(rate, h, radius) (read from its object code
// 2026-10-25); written so here.
void Tube::refresh_first_area() {
  constexpr double omega = kTwoPi * 160.0;
  if (ramp_left_ > 0) {
    std::size_t n = nodes_.size();
    for (std::size_t i = 0; i < n; ++i) {
      radius_target_[i] += radius_ramp_[i];
    }
    if (--ramp_left_ == 0) {
      radius_target_ = radius_goal_;
    }
  }
  double offset = radius_target_[0] - rest_radius_[0];
  bool radius_moving = std::abs(offset) >= 1.0e-9 || std::abs(radius_rate_[0]) >= 1.0e-6;
  if (radius_moving) {
    double damping = radius_rate_[0] * (-2.0 * omega);
    double pull = std::fma(offset, omega * omega, damping);
    radius_rate_[0] = std::fma(pull, kStepSeconds, radius_rate_[0]);
    rest_radius_[0] = std::max(std::fma(radius_rate_[0], kStepSeconds, rest_radius_[0]), 1.0e-4);
    fast_node_[0] = std::abs(radius_rate_[0]) > 0.5;
    moved_node_[0] = 1;
  } else {
    radius_rate_[0] = 0.0;
    fast_node_[0] = 0;
    moved_node_[0] = 0;
  }
  double radius = std::max(rest_radius_[0] + 0.0, 0.0);
  area_[0] = std::max(kPi * radius * radius, kMinimumArea);
}

// Thomas algorithm (symmetric tridiagonal), its elimination kept: the
// forward sweep's divisions each wait for the one before (the lining's
// step ran at 2.1 instructions a cycle, MEASURED on the vocal tube).
void Tube::solve_tridiagonal(std::size_t n) {
  const double* __restrict upper = upper_.data();
  const double* __restrict scratch = scratch_.data();
  const double* __restrict denominator = denominator_.data();
  double* __restrict rhs = rhs_.data();
  double carried = rhs[0] / diag_[0];
  rhs[0] = carried;
  for (std::size_t i = 1; i < n; ++i) {
    carried = std::fma(-upper[i - 1], carried, rhs[i]) / denominator[i];
    rhs[i] = carried;
  }
  for (std::size_t i = n - 1; i > 0; --i) {
    carried = std::fma(-scratch[i - 1], carried, rhs[i - 1]);
    rhs[i - 1] = carried;
  }
}

// The same solve when the lining has no axial coupling (every upper entry
// and every elimination factor a zero: AXIAL=0's vocal tube, REVNOSE's
// nose). Then fma(-0, carried, r) is r itself whenever r is not zero and
// carried is finite, and each node's division stands alone: all are taken
// at once, and the conditions checked after. Where one fails (a zero's sign
// or a non-finite value could differ) the sweep is not used and the
// elimination runs as it was (false returned, rhs untouched). Exactly the
// reference's values either way.
bool Tube::solve_diagonal(std::size_t n) {
  const double* __restrict denominator = denominator_.data();
  double* __restrict rhs = rhs_.data();
  double* __restrict quotient = quotient_.data();
  quotient[0] = rhs[0] / diag_[0];
  for (std::size_t i = 1; i < n; ++i) {
    quotient[i] = rhs[i] / denominator[i];
  }
  bool exact = std::isfinite(quotient[n - 1]);
  for (std::size_t i = 1; i < n; ++i) {
    // forward: rhs[i] kept by its fma; backward: quotient[i - 1] kept by its.
    exact = exact && rhs[i] != 0.0 && std::isfinite(quotient[i - 1]) && quotient[i - 1] != 0.0;
  }
  if (!exact) {
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    rhs[i] = quotient[i];
  }
  return true;
}

// The lining's matrix, its elimination and the bodies' backward-Euler
// denominators, from the compiled pieces (step_lining's arithmetic before
// they were kept): diag = cover_mass + (h^2 / 4) K_ii, upper = -(h^2 / 4)
// K_i,i+1, the elimination's denominators diag_i - upper_i-1 scratch_i-1,
// and body_mass + h (body_c + cover_c) + h^2 (coupling + support).
void Tube::factor_lining() {
  std::size_t n = nodes_.size();
  constexpr double h = kStepSeconds;
  constexpr double quarter = 0.25 * h * h;
  constexpr double negative_quarter = -quarter;
  for (std::size_t i = 0; i < n; ++i) {
    double k_diagonal = hoop_k_[i] + coupling_k_[i];
    if (i > 0) {
      k_diagonal += axial_k_[i - 1];
    }
    if (i + 1 < n) {
      k_diagonal += axial_k_[i];
    }
    diag_[i] = std::fma(k_diagonal, quarter, cover_mass_[i]);
    upper_[i] = i + 1 < n ? negative_quarter * axial_k_[i] : 0.0;
  }
  scratch_[0] = upper_[0] / diag_[0];
  for (std::size_t i = 1; i < n; ++i) {
    double denominator = std::fma(-upper_[i - 1], scratch_[i - 1], diag_[i]);
    denominator_[i] = denominator;
    scratch_[i] = upper_[i] / denominator;
  }
  for (std::size_t i = 0; i < n; ++i) {
    double stiffness = coupling_k_[i] + support_k_[i];
    double damping = body_c_[i] + cover_c_[i];
    body_denominator_[i] = std::fma(stiffness, h * h, std::fma(damping, h, body_mass_[i]));
  }
  diagonal_ = true;
  for (std::size_t i = 0; i < n; ++i) {
    diagonal_ = diagonal_ && upper_[i] == 0.0 && scratch_[i] == 0.0;
  }
  quotient_.resize(n);
  std::fill(lining_node_dirty_.begin(), lining_node_dirty_.end(), 0);
  lining_dirty_ = false;
  lining_partial_ = false;
}

// factor_lining for the nodes whose pieces were recompiled (a moving
// frame's), its arithmetic node by node: their diagonal and bodies'
// denominators, and the elimination from the first of them on while its
// factor changes (a node whose factor comes out the same leaves the next
// node's as it was). The coupling (upper entries) is not touched: a change
// of it factors the whole lining again.
void Tube::refactor_lining() {
  std::size_t n = nodes_.size();
  constexpr double h = kStepSeconds;
  constexpr double quarter = 0.25 * h * h;
  bool carry = false;
  for (std::size_t i = 0; i < n; ++i) {
    bool own = lining_node_dirty_[i] != 0;
    if (!own && !carry) {
      continue;
    }
    if (own) {
      double k_diagonal = hoop_k_[i] + coupling_k_[i];
      if (i > 0) {
        k_diagonal += axial_k_[i - 1];
      }
      if (i + 1 < n) {
        k_diagonal += axial_k_[i];
      }
      diag_[i] = std::fma(k_diagonal, quarter, cover_mass_[i]);
      double stiffness = coupling_k_[i] + support_k_[i];
      double damping = body_c_[i] + cover_c_[i];
      body_denominator_[i] = std::fma(stiffness, h * h, std::fma(damping, h, body_mass_[i]));
      lining_node_dirty_[i] = 0;
    }
    double previous = scratch_[i];
    if (i == 0) {
      scratch_[0] = upper_[0] / diag_[0];
    } else {
      double denominator = std::fma(-upper_[i - 1], scratch_[i - 1], diag_[i]);
      denominator_[i] = denominator;
      scratch_[i] = upper_[i] / denominator;
    }
    carry = std::memcmp(&previous, &scratch_[i], sizeof(double)) != 0;
  }
  lining_partial_ = false;
}

// Air: every node is a junction of its two air paths (and any port). Waves
// arriving this step set the junction pressure together with the lining's
// motion, solved as one scalar implicit exchange per node (no global solve):
// the lining's velocity and the air pressure it sees are found together.
void Tube::step_air() {
  if (plain_air()) {
    step_air_plain();
    return;
  }
  std::size_t n = static_cast<std::size_t>(active_count_);
  constexpr double h = kStepSeconds;
  constexpr double half_h = 0.5 * h;
  if (opening_node_ >= 0) {
    radiated_flow_ = 0.0;
  }
  if (opening_position_ >= 0.0) {
    follow_opening();
  }
  double share_now = opening_share_;
  // The opening's load (R parallel L for its area), at full weight (only
  // read with an opening): R, the trapezoid's gain on L, and the load's
  // admittance 1 / R + gain.
  double opening_resistance = 1.0;
  double slot_gain = 0.0;
  double slot_admittance = 0.0;
  if (opening_node_ >= 0) {
    opening_resistance = radiation_resistance_scale_ * 128.0 * kRhoC / (9.0 * kPi * kPi * opening_area_);
    double opening_inertance =
        radiation_mass_scale_ * 8.0 * kAirDensity / (3.0 * kPi * std::sqrt(kPi * opening_area_));
    slot_gain = half_h / opening_inertance;
    slot_admittance = 1.0 / opening_resistance + slot_gain;
  }
  // Lip radiation: a resistance R in parallel with an inertance L (the end
  // correction's air mass plus the tube past the last node), discretized by
  // the trapezoid rule for the inertance's flow; kept while the last active
  // node's area, the end's length and the scales hold.
  if (radiating_right_ && n > 0) {
    double a = std::max(area_[n - 1], kMinimumArea);
    if (!(a == radiation_area_ && end_extra_ == radiation_extra_ && radiation_mass_scale_ == radiation_mass_used_ &&
          radiation_resistance_scale_ == radiation_resistance_used_)) {
      radiation_area_ = a;
      radiation_extra_ = end_extra_;
      radiation_mass_used_ = radiation_mass_scale_;
      radiation_resistance_used_ = radiation_resistance_scale_;
      radiation_r_ = radiation_resistance_scale_ * 128.0 * kRhoC / (9.0 * kPi * kPi * a);
      radiation_inertance_ =
          radiation_mass_scale_ * 8.0 * kAirDensity / (3.0 * kPi * std::sqrt(kPi * a)) + kAirDensity * end_extra_ / a;
    }
  }
  double radiation_gain = radiating_right_ && n > 0 ? half_h / radiation_inertance_ : 0.0;
  double radiation_admittance = radiating_right_ && n > 0 ? 1.0 / radiation_r_ + radiation_gain : 0.0;
  // The members the loop reads as locals, and the arrays as local pointers
  // (each its own vector): a store to any double array may otherwise alias
  // a double member, and every member was read again after every store.
  const bool radiating = radiating_right_;
  const bool by_area = smoothing_by_area_;
  const double smoothing = wave_smoothing_;
  const double shelf = wave_shelf_;
  const double left_area = left_port_area_;
  const double left_incoming = left_incoming_;
  const double right_area = right_port_area_;
  const double right_incoming = right_incoming_;
  const double radiation_r = radiation_r_;
  const Index opening_node = opening_node_;
  // The side ports at each node, in port order: a chain from the node's
  // first (-1: none), built once a step (a scan of every port at every node,
  // twice, was a fifth of this loop's instructions).
  const std::size_t ports = side_node_.size();
  const double* __restrict side_area = side_area_.data();
  if (port_first_.size() != nodes_.size()) {
    port_first_.assign(nodes_.size(), -1);
  }
  port_next_.resize(ports);
  for (std::size_t port = ports; port-- > 0;) {
    Index node = side_node_[port];
    if (node >= 0 && static_cast<std::size_t>(node) < n) {
      port_next_[port] = port_first_[static_cast<std::size_t>(node)];
      port_first_[static_cast<std::size_t>(node)] = static_cast<int>(port);
    }
  }
  const int* __restrict port_first = port_first_.data();
  const int* __restrict port_next = port_next_.data();
  const double* __restrict side_incoming = side_incoming_.data();
  double* __restrict side_outgoing = side_outgoing_.data();
  const double* __restrict link_y = link_y_.data();
  const double* __restrict forward = forward_.data();
  const double* __restrict backward = backward_.data();
  const double* __restrict surface_of = surface_.data();
  const double* __restrict frame_flow = frame_flow_.data();
  const double* __restrict cover_mass = cover_mass_.data();
  const double* __restrict cover_c = cover_c_.data();
  const double* __restrict body_v = body_v_.data();
  const double* __restrict rest_radius = rest_radius_.data();
  const double* __restrict hop_k = hop_k_.data();
  double* __restrict cover_q = cover_q_.data();
  double* __restrict cover_v = cover_v_.data();
  double* __restrict pressure = pressure_.data();
  double* __restrict wall_flow = wall_flow_.data();
  double* __restrict junction_admittance = junction_admittance_.data();
  double* __restrict junction_resistance = junction_resistance_.data();
  double* __restrict junction_plus = junction_plus_.data();
  double* __restrict junction_minus = junction_minus_.data();
  char* __restrict junction_dirty = junction_dirty_.data();
  double* __restrict forward_smooth = forward_smooth_.data();
  double* __restrict backward_smooth = backward_smooth_.data();
  double* __restrict forward_raw = forward_raw_.data();
  double* __restrict backward_raw = backward_raw_.data();
  double* __restrict forward_next = forward_next_.data();
  double* __restrict backward_next = backward_next_.data();
  double total_wall = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    double admittance = 0.0;
    double twice_sum = 0.0;
    double from_left = 0.0;
    double from_right = 0.0;
    if (i > 0) {
      double y = link_y[i - 1];
      from_left = forward[i - 1];
      admittance += y;
      twice_sum = std::fma(2.0 * y, from_left, twice_sum);
    } else if (left_area > 0.0) {
      double y = left_area / kRhoC;
      from_left = left_incoming;
      admittance += y;
      twice_sum = std::fma(2.0 * y, from_left, twice_sum);
    }
    if (i + 1 < n) {
      double y = link_y[i];
      from_right = backward[i];
      admittance += y;
      twice_sum = std::fma(2.0 * y, from_right, twice_sum);
    } else if (radiating) {
      admittance += radiation_admittance;
      twice_sum -= std::fma(radiation_previous_pressure_, radiation_gain, radiation_inertance_flow_);
    } else if (right_area > 0.0) {
      double y = right_area / kRhoC;
      from_right = right_incoming;
      admittance += y;
      twice_sum = std::fma(2.0 * y, from_right, twice_sum);
    }
    // The opening's share at this node: R / w in parallel with L / w.
    int slot = -1;
    double weight = 0.0;
    if (opening_node >= 0) {
      if (static_cast<Index>(i) == opening_node) {
        slot = 0;
        weight = 1.0 - share_now;
      } else if (static_cast<Index>(i) == opening_node + 1) {
        slot = 1;
        weight = share_now;
      }
    }
    // (Each slot keeps the full load's inertance flow; the node draws its
    // share of it, so a share shrinking to nothing fades out, no step.)
    if (slot >= 0) {
      admittance = std::fma(slot_admittance, weight, admittance);
      double drawn = std::fma(slot_gain, opening_pressure_[slot], opening_flow_[slot]);
      twice_sum = std::fma(-weight, drawn, twice_sum);
    }
    // Side branches joined at this node.
    for (int port = port_first[i]; port >= 0; port = port_next[port]) {
      if (side_area[port] > 0.0) {
        double y = side_area[port] / kRhoC;
        admittance += y;
        twice_sum = std::fma(2.0 * y, side_incoming[port], twice_sum);
      }
    }
    if (admittance <= 0.0) {
      admittance = kMinimumArea / kRhoC;
    }
    // The lining: an outward velocity draws air out of the junction; so
    // does the frame as articulation widens the node.
    double surface = surface_of[i];
    twice_sum -= frame_flow[i];
    // The lining's resistance at this junction, surface^2 / admittance +
    // cover_c, and its implicit step's factors m + (h / 2) R and m - (h / 2)
    // R: kept while the admittance and the node's pieces hold.
    if (junction_dirty[i] || !(admittance == junction_admittance[i])) {
      double fresh = surface * surface / admittance + cover_c[i];
      double half_resistance = half_h * fresh;
      junction_resistance[i] = fresh;
      junction_plus[i] = half_resistance + cover_mass[i];
      junction_minus[i] = cover_mass[i] - half_resistance;
      junction_admittance[i] = admittance;
      junction_dirty[i] = 0;
    }
    double plus = junction_plus[i];
    double minus = junction_minus[i];
    double drive = surface * twice_sum / admittance + cover_c[i] * body_v[i];
    // Contact with the opposite wall: tissue meeting tissue, a stiff,
    // critically damped spring (about 1.5 kHz on the lining's mass) past the
    // contact radius. (A hard stop zeroed the lining's velocity in one step:
    // at every closure each glottal pulse pushed the walls apart and they
    // slammed back, a click every cycle, 84 a second through /dh/.)
    double radius = rest_radius[i] + cover_q[i];
    double penetration = kContactRadius - radius;
    if (penetration > 0.0) {
      constexpr double omega = kTwoPi * kContactFrequency;
      drive = std::fma(penetration, cover_mass[i] * omega * omega, drive);
      double resistance = std::fma(2.0 * cover_mass[i], omega, junction_resistance[i]);
      double half_resistance = half_h * resistance;
      plus = half_resistance + cover_mass[i];
      minus = cover_mass[i] - half_resistance;
    }
    double velocity = cover_v[i];
    double next = std::fma(h, drive, velocity * minus) / plus;
    double middle = 0.5 * (velocity + next);
    // (Never through the axis.)
    if (radius < -kContactRadius && next < 0.0) {
      cover_q[i] = -kContactRadius - rest_radius[i];
      next = 0.0;
      middle = 0.5 * velocity;
    }
    cover_v[i] = next;
    double drawn = surface * middle;
    double p = (twice_sum - drawn) / admittance;
    pressure[i] = p;
    if (slot >= 0) {
      opening_flow_[slot] = std::fma(slot_gain, p + opening_pressure_[slot], opening_flow_[slot]);
      radiated_flow_ = std::fma(weight, p / opening_resistance + opening_flow_[slot], radiated_flow_);
      opening_pressure_[slot] = p;
    }
    for (int port = port_first[i]; port >= 0; port = port_next[port]) {
      side_outgoing[port] = side_area[port] > 0.0 ? p - side_incoming[port] : 0.0;
    }
    wall_flow[i] = -drawn;
    total_wall += drawn;
    if (i > 0) {
      // Each hop passes a gentle one-pole low-pass (wave_smoothing_ k:
      // y = (1 - k) x + k y'): no loss at DC, a loss growing with frequency
      // as the boundary layers' and heat conduction's do, which the
      // frequency-flat viscous term leaves out.
      double wave = p - from_left;
      double k = by_area ? hop_k[i - 1] : smoothing;
      if (shelf > 0.0) {
        double gain = (1.0 - k) / (1.0 - shelf);
        backward_smooth[i - 1] =
            std::fma(gain, std::fma(-shelf, backward_raw[i - 1], wave), k * backward_smooth[i - 1]);
        backward_raw[i - 1] = wave;
      } else {
        backward_smooth[i - 1] = std::fma(wave, 1.0 - k, k * backward_smooth[i - 1]);
      }
      backward_next[i - 1] = backward_smooth[i - 1];
    } else {
      outgoing_left_ = p - from_left;
    }
    if (i + 1 < n) {
      double wave = p - from_right;
      double k = by_area ? hop_k[i] : smoothing;
      if (shelf > 0.0) {
        double gain = (1.0 - k) / (1.0 - shelf);
        forward_smooth[i] = std::fma(gain, std::fma(-shelf, forward_raw[i], wave), k * forward_smooth[i]);
        forward_raw[i] = wave;
      } else {
        forward_smooth[i] = std::fma(wave, 1.0 - k, k * forward_smooth[i]);
      }
      forward_next[i] = forward_smooth[i];
    } else if (radiating) {
      radiation_inertance_flow_ =
          std::fma(p + radiation_previous_pressure_, radiation_gain, radiation_inertance_flow_);
      radiation_previous_pressure_ = p;
      radiated_flow_ = p / radiation_r + radiation_inertance_flow_;
    } else {
      outgoing_right_ = p - from_right;
    }
  }
  total_wall_flow_ = total_wall;
  for (std::size_t port = 0; port < ports; ++port) {
    Index node = side_node_[port];
    if (node >= 0 && static_cast<std::size_t>(node) < n) {
      port_first_[static_cast<std::size_t>(node)] = -1;
    }
  }
  std::swap(forward_, forward_next_);
  std::swap(backward_, backward_next_);
  scatter_links();
  std::fill(series_.begin(), series_.end(), 0.0);
}

// ---- The plain tube's air (the trachea): step_air and scatter_links for a
// tube with its two end ports open and nothing else, nodes two to a
// vector. Each value is step_air's expression on the same inputs: the
// same operations in the same order, each lane rounding as the scalar
// (a vector fma is the scalar std::fma, per lane); the junction's lining
// factors and the links' kinetic coefficients are taken every step (they
// are pure functions of inputs that, when unchanged, give the values
// kept); the drawn flows are summed in node order after. MEASURED in
// NOTES_residual_speed.md. ----

// Whether step_air's general loop reduces to the plain one: no opening, no
// side ports, no radiation, no seal, the hop loss constant and a plain
// one-pole (no shelf), every node active, both end ports open.
bool Tube::plain_air() const {
  std::size_t n = nodes_.size();
  bool ends = left_port_area_ > 0.0 && right_port_area_ > 0.0;
  bool bare = side_node_.empty() && opening_node_ < 0 && opening_position_ < 0.0 && !radiating_right_;
  bool open = seal_area_ <= 0.0 && !smoothing_by_area_ && !(wave_shelf_ > 0.0);
  return ends && bare && open && n >= 4 && static_cast<std::size_t>(active_count_) == n;
}

// step_air's loop body for node i of a plain tube (its ends included).
void Tube::air_node_plain(std::size_t i, std::size_t n, double* drawn_out) {
  constexpr double h = kStepSeconds;
  constexpr double half_h = 0.5 * h;
  double admittance = 0.0;
  double twice_sum = 0.0;
  double from_left = 0.0;
  double from_right = 0.0;
  if (i > 0) {
    double y = link_y_[i - 1];
    from_left = forward_[i - 1];
    admittance += y;
    twice_sum = std::fma(2.0 * y, from_left, twice_sum);
  } else {
    double y = left_port_area_ / kRhoC;
    from_left = left_incoming_;
    admittance += y;
    twice_sum = std::fma(2.0 * y, from_left, twice_sum);
  }
  if (i + 1 < n) {
    double y = link_y_[i];
    from_right = backward_[i];
    admittance += y;
    twice_sum = std::fma(2.0 * y, from_right, twice_sum);
  } else {
    double y = right_port_area_ / kRhoC;
    from_right = right_incoming_;
    admittance += y;
    twice_sum = std::fma(2.0 * y, from_right, twice_sum);
  }
  if (admittance <= 0.0) {
    admittance = kMinimumArea / kRhoC;
  }
  double surface = surface_[i];
  twice_sum -= frame_flow_[i];
  if (junction_dirty_[i] || !(admittance == junction_admittance_[i])) {
    double fresh = surface * surface / admittance + cover_c_[i];
    double half_resistance = half_h * fresh;
    junction_resistance_[i] = fresh;
    junction_plus_[i] = half_resistance + cover_mass_[i];
    junction_minus_[i] = cover_mass_[i] - half_resistance;
    junction_admittance_[i] = admittance;
    junction_dirty_[i] = 0;
  }
  double plus = junction_plus_[i];
  double minus = junction_minus_[i];
  double drive = surface * twice_sum / admittance + cover_c_[i] * body_v_[i];
  double radius = rest_radius_[i] + cover_q_[i];
  double penetration = kContactRadius - radius;
  if (penetration > 0.0) {
    constexpr double omega = kTwoPi * kContactFrequency;
    drive = std::fma(penetration, cover_mass_[i] * omega * omega, drive);
    double resistance = std::fma(2.0 * cover_mass_[i], omega, junction_resistance_[i]);
    double half_resistance = half_h * resistance;
    plus = half_resistance + cover_mass_[i];
    minus = cover_mass_[i] - half_resistance;
  }
  double velocity = cover_v_[i];
  double next = std::fma(h, drive, velocity * minus) / plus;
  double middle = 0.5 * (velocity + next);
  if (radius < -kContactRadius && next < 0.0) {
    cover_q_[i] = -kContactRadius - rest_radius_[i];
    next = 0.0;
    middle = 0.5 * velocity;
  }
  cover_v_[i] = next;
  double drawn = surface * middle;
  double p = (twice_sum - drawn) / admittance;
  pressure_[i] = p;
  wall_flow_[i] = -drawn;
  drawn_out[i] = drawn;
  double k = wave_smoothing_;
  if (i > 0) {
    double wave = p - from_left;
    backward_smooth_[i - 1] = std::fma(wave, 1.0 - k, k * backward_smooth_[i - 1]);
    backward_next_[i - 1] = backward_smooth_[i - 1];
  } else {
    outgoing_left_ = p - from_left;
  }
  if (i + 1 < n) {
    double wave = p - from_right;
    forward_smooth_[i] = std::fma(wave, 1.0 - k, k * forward_smooth_[i]);
    forward_next_[i] = forward_smooth_[i];
  } else {
    outgoing_right_ = p - from_right;
  }
}

// air_node_plain for the interior nodes i and i + 1 (lanes 0 and 1), when
// neither touches the opposite wall and both admittances are positive (the
// caller has checked; otherwise it runs air_node_plain on each).
void Tube::air_pair_plain(std::size_t i, double* drawn_out) {
  constexpr double h = kStepSeconds;
  constexpr double half_h = 0.5 * h;
  float64x2_t zero = vdupq_n_f64(0.0);
  float64x2_t two = vdupq_n_f64(2.0);
  float64x2_t y_left = vld1q_f64(&link_y_[i - 1]);
  float64x2_t from_left = vld1q_f64(&forward_[i - 1]);
  float64x2_t y_right = vld1q_f64(&link_y_[i]);
  float64x2_t from_right = vld1q_f64(&backward_[i]);
  float64x2_t admittance = vaddq_f64(vaddq_f64(zero, y_left), y_right);
  float64x2_t twice_sum = vfmaq_f64(zero, vmulq_f64(two, y_left), from_left);
  twice_sum = vfmaq_f64(twice_sum, vmulq_f64(two, y_right), from_right);
  float64x2_t surface = vld1q_f64(&surface_[i]);
  twice_sum = vsubq_f64(twice_sum, vld1q_f64(&frame_flow_[i]));
  float64x2_t cover_c = vld1q_f64(&cover_c_[i]);
  float64x2_t cover_mass = vld1q_f64(&cover_mass_[i]);
  float64x2_t fresh = vaddq_f64(vdivq_f64(vmulq_f64(surface, surface), admittance), cover_c);
  float64x2_t half_resistance = vmulq_f64(vdupq_n_f64(half_h), fresh);
  float64x2_t plus = vaddq_f64(half_resistance, cover_mass);
  float64x2_t minus = vsubq_f64(cover_mass, half_resistance);
  vst1q_f64(&junction_resistance_[i], fresh);
  vst1q_f64(&junction_plus_[i], plus);
  vst1q_f64(&junction_minus_[i], minus);
  vst1q_f64(&junction_admittance_[i], admittance);
  junction_dirty_[i] = 0;
  junction_dirty_[i + 1] = 0;
  float64x2_t drive = vaddq_f64(vdivq_f64(vmulq_f64(surface, twice_sum), admittance),
                                vmulq_f64(cover_c, vld1q_f64(&body_v_[i])));
  float64x2_t velocity = vld1q_f64(&cover_v_[i]);
  float64x2_t next = vdivq_f64(vfmaq_f64(vmulq_f64(velocity, minus), vdupq_n_f64(h), drive), plus);
  float64x2_t middle = vmulq_f64(vdupq_n_f64(0.5), vaddq_f64(velocity, next));
  vst1q_f64(&cover_v_[i], next);
  float64x2_t drawn = vmulq_f64(surface, middle);
  float64x2_t p = vdivq_f64(vsubq_f64(twice_sum, drawn), admittance);
  vst1q_f64(&pressure_[i], p);
  vst1q_f64(&wall_flow_[i], vnegq_f64(drawn));
  vst1q_f64(&drawn_out[i], drawn);
  float64x2_t k = vdupq_n_f64(wave_smoothing_);
  float64x2_t pass = vdupq_n_f64(1.0 - wave_smoothing_);
  float64x2_t backward_wave = vsubq_f64(p, from_left);
  float64x2_t backward_smooth = vfmaq_f64(vmulq_f64(k, vld1q_f64(&backward_smooth_[i - 1])), backward_wave, pass);
  vst1q_f64(&backward_smooth_[i - 1], backward_smooth);
  vst1q_f64(&backward_next_[i - 1], backward_smooth);
  float64x2_t forward_wave = vsubq_f64(p, from_right);
  float64x2_t forward_smooth = vfmaq_f64(vmulq_f64(k, vld1q_f64(&forward_smooth_[i])), forward_wave, pass);
  vst1q_f64(&forward_smooth_[i], forward_smooth);
  vst1q_f64(&forward_next_[i], forward_smooth);
}

void Tube::step_air_plain() {
  std::size_t n = nodes_.size();
  if (drawn_.size() != n) {
    drawn_.assign(n, 0.0);
  }
  double* drawn = drawn_.data();
  air_node_plain(0, n, drawn);
  std::size_t i = 1;
  while (i + 2 < n) {
    // (Nodes i and i + 1, both interior.)
    bool positive = link_y_[i - 1] + link_y_[i] > 0.0 && link_y_[i] + link_y_[i + 1] > 0.0;
    bool apart = kContactRadius - (rest_radius_[i] + cover_q_[i]) <= 0.0 &&
                 kContactRadius - (rest_radius_[i + 1] + cover_q_[i + 1]) <= 0.0;
    if (positive && apart) {
      air_pair_plain(i, drawn);
    } else {
      air_node_plain(i, n, drawn);
      air_node_plain(i + 1, n, drawn);
    }
    i += 2;
  }
  while (i < n) {
    air_node_plain(i, n, drawn);
    i += 1;
  }
  double total_wall = 0.0;
  for (std::size_t node = 0; node < n; ++node) {
    total_wall += drawn[node];
  }
  total_wall_flow_ = total_wall;
  std::swap(forward_, forward_next_);
  std::swap(backward_, backward_next_);
  scatter_links_plain();
  std::fill(series_.begin(), series_.end(), 0.0);
}

// scatter_links without a seal, links two to a vector (each lane
// scatter_links' expressions; the kinetic coefficients taken every step).
void Tube::scatter_links_plain() {
  std::size_t n = nodes_.size();
  float64x2_t zero = vdupq_n_f64(0.0);
  float64x2_t one = vdupq_n_f64(1.0);
  float64x2_t two = vdupq_n_f64(2.0);
  float64x2_t least = vdupq_n_f64(kMinimumArea);
  float64x2_t coefficient = vdupq_n_f64(0.5 * kAirDensity);
  std::size_t i = 0;
  while (i + 2 < n) {
    float64x2_t first = vld1q_f64(&area_[i]);
    float64x2_t second = vld1q_f64(&area_[i + 1]);
    // std::max(area, kMinimumArea): (area < kMinimumArea) ? kMinimumArea : area.
    float64x2_t first_floored = vbslq_f64(vcltq_f64(first, least), least, first);
    float64x2_t second_floored = vbslq_f64(vcltq_f64(second, least), least, second);
    float64x2_t right_difference = vsubq_f64(vdivq_f64(one, first_floored), vdivq_f64(one, second));
    float64x2_t left_difference = vsubq_f64(vdivq_f64(one, second_floored), vdivq_f64(one, first));
    float64x2_t right_k = vmulq_f64(vmulq_f64(coefficient, right_difference), right_difference);
    float64x2_t left_k = vmulq_f64(vmulq_f64(coefficient, left_difference), left_difference);
    right_k = vbslq_f64(vcgtq_f64(second, first), right_k, zero);
    left_k = vbslq_f64(vcgtq_f64(first, second), left_k, zero);
    vst1q_f64(&kinetic_first_area_[i], first);
    vst1q_f64(&kinetic_second_area_[i], second);
    vst1q_f64(&kinetic_right_[i], right_k);
    vst1q_f64(&kinetic_left_[i], left_k);
    float64x2_t arriving_forward = vld1q_f64(&forward_[i]);
    float64x2_t arriving_backward = vld1q_f64(&backward_[i]);
    float64x2_t z = vld1q_f64(&link_z_[i]);
    float64x2_t drive = vaddq_f64(vmulq_f64(two, vsubq_f64(arriving_forward, arriving_backward)),
                                  vld1q_f64(&series_[i]));
    uint64x2_t onward = vcgeq_f64(drive, zero);
    float64x2_t kinetic = vbslq_f64(onward, right_k, left_k);
    float64x2_t linear = vaddq_f64(vmulq_f64(two, z), vld1q_f64(&link_resistance_[i]));
    float64x2_t size = vabsq_f64(drive);
    float64x2_t root = vsqrtq_f64(vfmaq_f64(vmulq_f64(vmulq_f64(vdupq_n_f64(4.0), kinetic), size), linear, linear));
    float64x2_t magnitude = vdivq_f64(vmulq_f64(two, size), vaddq_f64(linear, root));
    float64x2_t flow = vbslq_f64(onward, magnitude, vnegq_f64(magnitude));
    vst1q_f64(&link_flow_[i], flow);
    float64x2_t carried = vmulq_f64(z, flow);
    vst1q_f64(&forward_[i], vaddq_f64(arriving_backward, carried));
    vst1q_f64(&backward_[i], vsubq_f64(arriving_forward, carried));
    i += 2;
  }
  while (i + 1 < n) {
    double first_area = area_[i];
    double second_area = area_[i + 1];
    if (!(first_area == kinetic_first_area_[i] && second_area == kinetic_second_area_[i])) {
      kinetic_first_area_[i] = first_area;
      kinetic_second_area_[i] = second_area;
      kinetic_right_[i] = 0.0;
      kinetic_left_[i] = 0.0;
      if (second_area > first_area) {
        double difference = 1.0 / std::max(first_area, kMinimumArea) - 1.0 / second_area;
        kinetic_right_[i] = 0.5 * kAirDensity * difference * difference;
      }
      if (first_area > second_area) {
        double difference = 1.0 / std::max(second_area, kMinimumArea) - 1.0 / first_area;
        kinetic_left_[i] = 0.5 * kAirDensity * difference * difference;
      }
    }
    double arriving_forward = forward_[i];
    double arriving_backward = backward_[i];
    double z = link_z_[i];
    double drive = 2.0 * (arriving_forward - arriving_backward) + series_[i];
    double kinetic = drive >= 0.0 ? kinetic_right_[i] : kinetic_left_[i];
    double linear = 2.0 * z + link_resistance_[i];
    double magnitude = 2.0 * std::abs(drive) /
                       (linear + std::sqrt(std::fma(linear, linear, 4.0 * kinetic * std::abs(drive))));
    double flow = drive >= 0.0 ? magnitude : -magnitude;
    link_flow_[i] = flow;
    double carried = z * flow;
    forward_[i] = arriving_backward + carried;
    backward_[i] = arriving_forward - carried;
    i += 1;
  }
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
  // series_element inlined: the reference build fused only the square
  // under the root.
  std::size_t n = static_cast<std::size_t>(active_count_);
  const double seal = seal_area_;
  const double* __restrict area = area_.data();
  const double* __restrict link_z = link_z_.data();
  const double* __restrict link_resistance = link_resistance_.data();
  const double* __restrict series = series_.data();
  double* __restrict forward = forward_.data();
  double* __restrict backward = backward_.data();
  double* __restrict link_flow = link_flow_.data();
  double* __restrict kinetic_right = kinetic_right_.data();
  double* __restrict kinetic_left = kinetic_left_.data();
  double* __restrict kinetic_first_area = kinetic_first_area_.data();
  double* __restrict kinetic_second_area = kinetic_second_area_.data();
  if (seal_blend_ && seal_previous_least_.size() + 1 < n) {
    seal_previous_least_.resize(n - 1, std::numeric_limits<double>::quiet_NaN());
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double first_area = area[i];
    double second_area = area[i + 1];
    double open_share = 1.0;
    if (seal > 0.0 && seal_blend_) {
      double least = std::min(first_area, second_area);
      double before = seal_previous_least_[i];
      seal_previous_least_[i] = least;
      if (std::isnan(before)) {
        before = least;
      }
      if (least > seal && before <= seal) {
        open_share = (least - seal) / (least - before);
      } else if (least <= seal && before > seal) {
        open_share = (before - seal) / (before - least);
      } else if (least <= seal) {
        open_share = 0.0;
      }
      if (open_share <= 0.0) {
        link_flow[i] = 0.0;
        double arriving_right = forward[i];
        forward[i] = backward[i];
        backward[i] = arriving_right;
        continue;
      }
    } else if (seal > 0.0 && (first_area <= seal || second_area <= seal)) {
      // Sealed: no flow, each wave turned back.
      link_flow[i] = 0.0;
      double arriving_right = forward[i];
      forward[i] = backward[i];
      backward[i] = arriving_right;
      continue;
    }
    // The kinetic coefficients for flow toward each side (series_element's
    // jet loss out of a narrowing), kept while the link's two node areas
    // hold.
    if (!(first_area == kinetic_first_area[i] && second_area == kinetic_second_area[i])) {
      kinetic_first_area[i] = first_area;
      kinetic_second_area[i] = second_area;
      kinetic_right[i] = 0.0;
      kinetic_left[i] = 0.0;
      if (second_area > first_area) {
        double difference = 1.0 / std::max(first_area, kMinimumArea) - 1.0 / second_area;
        kinetic_right[i] = 0.5 * kAirDensity * difference * difference;
      }
      if (first_area > second_area) {
        double difference = 1.0 / std::max(second_area, kMinimumArea) - 1.0 / first_area;
        kinetic_left[i] = 0.5 * kAirDensity * difference * difference;
      }
    }
    double arriving_forward = forward[i];
    double arriving_backward = backward[i];
    double z = link_z[i];
    double drive = 2.0 * (arriving_forward - arriving_backward) + series[i];
    double kinetic = drive >= 0.0 ? kinetic_right[i] : kinetic_left[i];
    double linear = 2.0 * z + link_resistance[i];
    double magnitude = 2.0 * std::abs(drive) /
                       (linear + std::sqrt(std::fma(linear, linear, 4.0 * kinetic * std::abs(drive))));
    double flow = drive >= 0.0 ? magnitude : -magnitude;
    if (open_share < 1.0) {
      flow *= open_share;
    }
    link_flow[i] = flow;
    double carried = z * flow;
    forward[i] = arriving_backward + carried;
    backward[i] = arriving_forward - carried;
  }
}

}  // namespace tract
