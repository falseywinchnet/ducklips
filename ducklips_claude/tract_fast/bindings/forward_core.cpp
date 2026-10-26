// A hidden window's forward run (observer._forward_steps) to the bit: the
// numba kernel's step, expression for expression (-ffp-contract=off, the
// same libm and square root), with each node's and link's loop as there.
// Its cost was 25 ns a node-step in numba (MEASURED in
// NOTES_residual_speed.md); here the per-step views and the array
// reference counts are gone, and the jet's coefficient cache is per link
// (a pure function of the corner: the values do not depend on it).

#include "peel_core.h"
#include "jet_core.hpp"

#include "../include/tract/lanes.hpp"   // <arm_neon.h> on ARM; the same names on SSE2 elsewhere

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

using Index = std::int64_t;

constexpr double kRhoC = 1.14e-3 * 35000.0;    // observer.RHO_C, marching.RHO_C
constexpr double kStep = 1.0 / 96000.0;
constexpr double kPi = 3.141592653589793;
constexpr Index kColumns = 12;
constexpr Index kState = 19;                   // observer.STATE_COLUMNS

// marching._radiation.
void radiation(double area, double end_extra, double& resistance, double& gain) {
  resistance = 1.0 * 128.0 * kRhoC / (9.0 * (kPi * kPi) * area);
  double inertance = 1.0 * 8.0 * 1.14e-3 / (3.0 * kPi * std::sqrt(kPi * area)) + 1.14e-3 * end_extra / area;
  gain = 0.5 / 96000.0 / inertance;
}

// marching._lining on the node entry `walls` (geometry columns).
void lining(const double* walls, double& q, double& v, double& bq, double& bv) {
  constexpr double h = kStep;
  constexpr double quarter = 0.25 * h * h;
  double mc = walls[4];
  if (mc <= 0.0) {
    return;
  }
  double mb = walls[5];
  double hoop = walls[6];
  double coupling = walls[7];
  double support = walls[8];
  double cc = walls[9];
  double bc = walls[10];
  double diag = mc + quarter * (hoop + coupling);
  double rhs = mc * (2.0 * q + h * v) + 0.5 * h * h * coupling * bq;
  double new_q = rhs / diag - q;
  v = 2.0 * (new_q - q) / h - v;
  q = new_q;
  double stiffness = coupling + support;
  double pull = coupling * q - stiffness * bq + cc * v;
  bv = (mb * bv + h * pull) / (mb + h * (bc + cc) + h * h * stiffness);
  bq += h * bv;
}

// observer._junction_with: the node's pressure, its cover velocity stepped.
double junction(const double* walls, double known, double admittance, double& v, double bv) {
  constexpr double h = kStep;
  double mc = walls[4];
  if (mc <= 0.0) {
    double plain = known / admittance;
    return plain;
  }
  double surface = walls[3];
  double cc = walls[9];
  double resistance = surface * surface / admittance + cc;
  double drive = surface * known / admittance + cc * bv;
  double following = ((mc - 0.5 * h * resistance) * v + h * drive) / (mc + 0.5 * h * resistance);
  double middle = 0.5 * (v + following);
  v = following;
  double pressure = (known - surface * middle) / admittance;
  return pressure;
}

// observer.seal_share for link i at `step`.
double seal_share(const double* geometry, Index n, Index step, Index i, double seal, bool blend) {
  if (seal <= 0.0) {
    return 1.0;
  }
  const double* row = geometry + step * n * kColumns;
  double least = std::min(row[i * kColumns], row[(i + 1) * kColumns]);
  if (!blend) {
    return least <= seal ? 0.0 : 1.0;
  }
  double before = least;
  if (step > 0) {
    const double* earlier = row - n * kColumns;
    before = std::min(earlier[i * kColumns], earlier[(i + 1) * kColumns]);
  }
  if (least > seal && before <= seal) {
    double opening = (least - seal) / (least - before);
    return opening;
  }
  if (least <= seal && before > seal) {
    double closing = (before - seal) / (before - least);
    return closing;
  }
  return least <= seal ? 0.0 : 1.0;
}

// Two values from two places as a vector (lane 0 the first).
float64x2_t lanes(double low, double high) {
  float64x2_t pair = vsetq_lane_f64(high, vdupq_n_f64(low), 1);
  return pair;
}

// lining() on nodes i and i + 1 (their entries walls, walls + 12; their
// state rows node, node + 19), lane for lane (both linings present).
void lining_pair(const double* walls, double* node) {
  constexpr double h = kStep;
  constexpr double quarter = 0.25 * h * h;
  constexpr double half_h_squared = 0.5 * h * h;
  constexpr double h_squared = h * h;
  const double* other = walls + kColumns;
  double* next = node + kState;
  float64x2_t mc = lanes(walls[4], other[4]);
  float64x2_t mb = lanes(walls[5], other[5]);
  float64x2_t hoop = lanes(walls[6], other[6]);
  float64x2_t coupling = lanes(walls[7], other[7]);
  float64x2_t support = lanes(walls[8], other[8]);
  float64x2_t cc = lanes(walls[9], other[9]);
  float64x2_t bc = lanes(walls[10], other[10]);
  float64x2_t q = lanes(node[6], next[6]);
  float64x2_t v = lanes(node[7], next[7]);
  float64x2_t bq = lanes(node[8], next[8]);
  float64x2_t bv = lanes(node[9], next[9]);
  float64x2_t diag = vaddq_f64(mc, vmulq_f64(vdupq_n_f64(quarter), vaddq_f64(hoop, coupling)));
  float64x2_t momentum = vaddq_f64(vmulq_f64(vdupq_n_f64(2.0), q), vmulq_f64(vdupq_n_f64(h), v));
  float64x2_t pulled = vmulq_f64(vmulq_f64(vdupq_n_f64(half_h_squared), coupling), bq);
  float64x2_t rhs = vaddq_f64(vmulq_f64(mc, momentum), pulled);
  float64x2_t new_q = vsubq_f64(vdivq_f64(rhs, diag), q);
  v = vsubq_f64(vdivq_f64(vmulq_f64(vdupq_n_f64(2.0), vsubq_f64(new_q, q)), vdupq_n_f64(h)), v);
  q = new_q;
  float64x2_t stiffness = vaddq_f64(coupling, support);
  float64x2_t pull = vaddq_f64(vsubq_f64(vmulq_f64(coupling, q), vmulq_f64(stiffness, bq)), vmulq_f64(cc, v));
  float64x2_t body_top = vaddq_f64(vmulq_f64(mb, bv), vmulq_f64(vdupq_n_f64(h), pull));
  float64x2_t body_bottom = vaddq_f64(vaddq_f64(mb, vmulq_f64(vdupq_n_f64(h), vaddq_f64(bc, cc))),
                                      vmulq_f64(vdupq_n_f64(h_squared), stiffness));
  bv = vdivq_f64(body_top, body_bottom);
  bq = vaddq_f64(bq, vmulq_f64(vdupq_n_f64(h), bv));
  node[6] = vgetq_lane_f64(q, 0);
  next[6] = vgetq_lane_f64(q, 1);
  node[7] = vgetq_lane_f64(v, 0);
  next[7] = vgetq_lane_f64(v, 1);
  node[8] = vgetq_lane_f64(bq, 0);
  next[8] = vgetq_lane_f64(bq, 1);
  node[9] = vgetq_lane_f64(bv, 0);
  next[9] = vgetq_lane_f64(bv, 1);
}

// The junctions of interior nodes i and i + 1 (both with linings): their
// pressures into state rows [10], the covers' velocities stepped.
void junction_pair(const double* row, const double* y, double* state, Index i) {
  constexpr double h = kStep;
  const double* here = row + i * kColumns;
  const double* there = here + kColumns;
  double* node = state + i * kState;
  double* next = node + kState;
  float64x2_t y_left = lanes(y[i - 1], y[i]);
  float64x2_t y_right = lanes(y[i], y[i + 1]);
  float64x2_t arriving = lanes(state[(i - 1) * kState], node[0]);
  float64x2_t returning = lanes(node[1], next[1]);
  float64x2_t frame = lanes(here[1], there[1]);
  float64x2_t two = vdupq_n_f64(2.0);
  float64x2_t known = vsubq_f64(vaddq_f64(vmulq_f64(vmulq_f64(two, y_left), arriving),
                                          vmulq_f64(vmulq_f64(two, y_right), returning)), frame);
  float64x2_t admittance = vaddq_f64(y_left, y_right);
  float64x2_t mc = lanes(here[4], there[4]);
  float64x2_t surface = lanes(here[3], there[3]);
  float64x2_t cc = lanes(here[9], there[9]);
  float64x2_t v = lanes(node[7], next[7]);
  float64x2_t bv = lanes(node[9], next[9]);
  float64x2_t resistance = vaddq_f64(vdivq_f64(vmulq_f64(surface, surface), admittance), cc);
  float64x2_t drive = vaddq_f64(vdivq_f64(vmulq_f64(surface, known), admittance), vmulq_f64(cc, bv));
  float64x2_t half_resistance = vmulq_f64(vdupq_n_f64(0.5 * h), resistance);
  float64x2_t top = vaddq_f64(vmulq_f64(vsubq_f64(mc, half_resistance), v), vmulq_f64(vdupq_n_f64(h), drive));
  float64x2_t following = vdivq_f64(top, vaddq_f64(mc, half_resistance));
  float64x2_t middle = vmulq_f64(vdupq_n_f64(0.5), vaddq_f64(v, following));
  float64x2_t pressure = vdivq_f64(vsubq_f64(known, vmulq_f64(surface, middle)), admittance);
  node[7] = vgetq_lane_f64(following, 0);
  next[7] = vgetq_lane_f64(following, 1);
  node[10] = vgetq_lane_f64(pressure, 0);
  next[10] = vgetq_lane_f64(pressure, 1);
}

// The hop's one-poles on links i and i + 1.
void hop_pair(double* state, Index i, double shelf, double hop, double hop_zero) {
  double* first = state + i * kState;
  double* second = first + kState;
  double* third = second + kState;
  float64x2_t gain = vdupq_n_f64(shelf);
  float64x2_t pole = vdupq_n_f64(hop);
  float64x2_t zero = vdupq_n_f64(hop_zero);
  float64x2_t backward = vsubq_f64(lanes(first[10], second[10]), lanes(first[1], second[1]));
  float64x2_t backward_smooth = vaddq_f64(vmulq_f64(gain, vsubq_f64(backward, vmulq_f64(zero, lanes(first[4], second[4])))),
                                          vmulq_f64(pole, lanes(first[2], second[2])));
  float64x2_t forward = vsubq_f64(lanes(second[10], third[10]), lanes(first[0], second[0]));
  float64x2_t forward_smooth = vaddq_f64(vmulq_f64(gain, vsubq_f64(forward, vmulq_f64(zero, lanes(first[5], second[5])))),
                                         vmulq_f64(pole, lanes(first[3], second[3])));
  first[2] = vgetq_lane_f64(backward_smooth, 0);
  second[2] = vgetq_lane_f64(backward_smooth, 1);
  first[4] = vgetq_lane_f64(backward, 0);
  second[4] = vgetq_lane_f64(backward, 1);
  first[3] = vgetq_lane_f64(forward_smooth, 0);
  second[3] = vgetq_lane_f64(forward_smooth, 1);
  first[5] = vgetq_lane_f64(forward, 0);
  second[5] = vgetq_lane_f64(forward, 1);
}

// The series elements of links i and i + 1, both fully open (seal share
// exactly 1: no share multiplies the flow). The narrow and wide sides are
// selected by the drive's sign per lane and the kinetic coefficient set to
// 0 where the wide side is not wider (computed for both lanes: the areas
// are finite and floored, nothing exceptional arises in a discarded one).
void scatter_pair(const double* row, const double* y, double* state, Index i, bool jet, bool kinetic) {
  double* first = state + i * kState;
  double* second = first + kState;
  float64x2_t one = vdupq_n_f64(1.0);
  float64x2_t two = vdupq_n_f64(2.0);
  float64x2_t zero = vdupq_n_f64(0.0);
  float64x2_t z = vdivq_f64(one, lanes(y[i], y[i + 1]));
  float64x2_t smooth_forward = lanes(first[2], second[2]);
  float64x2_t smooth_backward = lanes(first[3], second[3]);
  float64x2_t drive = vmulq_f64(two, vsubq_f64(smooth_forward, smooth_backward));
  if (jet) {
    drive = vaddq_f64(drive, lanes(first[18], second[18]));
  }
  float64x2_t area_first = lanes(row[i * kColumns], row[(i + 1) * kColumns]);
  float64x2_t area_second = lanes(row[(i + 1) * kColumns], row[(i + 2) * kColumns]);
  uint64x2_t onward = vcgeq_f64(drive, zero);
  float64x2_t narrow = vbslq_f64(onward, area_first, area_second);
  float64x2_t wide = vbslq_f64(onward, area_second, area_first);
  float64x2_t k = zero;
  if (kinetic) {
    float64x2_t floor = vdupq_n_f64(1e-4);
    float64x2_t floored = vbslq_f64(vcltq_f64(narrow, floor), floor, narrow);
    float64x2_t difference = vsubq_f64(vdivq_f64(one, floored), vdivq_f64(one, wide));
    float64x2_t coefficient = vmulq_f64(vmulq_f64(vdupq_n_f64(0.5 * 1.14e-3), difference), difference);
    k = vbslq_f64(vcgtq_f64(wide, narrow), coefficient, zero);
  }
  float64x2_t linear = vaddq_f64(vmulq_f64(two, z), lanes(row[i * kColumns + 2], row[(i + 1) * kColumns + 2]));
  float64x2_t size = vabsq_f64(drive);
  float64x2_t root = vsqrtq_f64(vaddq_f64(vmulq_f64(linear, linear), vmulq_f64(vmulq_f64(vdupq_n_f64(4.0), k), size)));
  float64x2_t magnitude = vdivq_f64(vmulq_f64(two, size), vaddq_f64(linear, root));
  float64x2_t flow = vbslq_f64(onward, magnitude, vnegq_f64(magnitude));
  float64x2_t carried = vmulq_f64(z, flow);
  float64x2_t arriving_next = vaddq_f64(smooth_backward, carried);
  float64x2_t returning = vsubq_f64(smooth_forward, carried);
  first[1] = vgetq_lane_f64(returning, 0);
  second[1] = vgetq_lane_f64(returning, 1);
  first[0] = vgetq_lane_f64(arriving_next, 0);
  second[0] = vgetq_lane_f64(arriving_next, 1);
  first[13] = vgetq_lane_f64(flow, 0);
  second[13] = vgetq_lane_f64(flow, 1);
}

}  // namespace

extern "C" void forward_steps(const forward_inputs* inputs) {
  const forward_inputs in = *inputs;
  Index n = in.nodes;
  double* state = in.state;
  std::vector<double> y(static_cast<std::size_t>(n - 1));
  std::vector<double> caches(static_cast<std::size_t>(4 * (n - 1)), std::nan(""));   // per link (see the head)
  double shelf = (1.0 - in.hop) / (1.0 - in.hop_zero);
  jet_core::JetParameters jet;
  jet.bound = in.jet_bound;
  jet.closing = in.jet_closing;
  jet.noise = in.noise;
  jet.noise_links = in.noise_links;
  jet.noise_steps = in.noise_steps;
  jet.seed = in.jet_seed;
  for (Index step = in.t0; step < in.t1; ++step) {
    const double* row = in.geometry + step * n * kColumns;
    const double* before_row = in.geometry + std::max<Index>(step - 1, 0) * n * kColumns;
    if (in.jet_gain > 0.0) {
      // Each link's source from its own flow the step before.
      Index narrowest = -1;
      double least = 0.0;
      jet_core::narrowest_of(before_row, kColumns, n, narrowest, least);
      double gain = in.jet_gain;
      if (in.has_gate) {
        gain = in.jet_gain * in.gate[std::min(step, in.gate_length - 1)];
      }
      for (Index i = 0; i + 1 < n; ++i) {
        double* link_state = state + i * kState;
        jet_core::JetState held;
        for (Index c = 0; c < 4; ++c) {
          held.filter[c] = link_state + 14 + c;
          held.cache[c] = caches.data() + i * 4 + c;
        }
        link_state[18] = jet_core::jet_source(jet, i, link_state[13], narrowest, least, held, gain, step);
      }
    }
    Index lining_index = 0;
    while (lining_index + 1 < n) {
      const double* walls = before_row + lining_index * kColumns;
      if (walls[4] > 0.0 && walls[kColumns + 4] > 0.0) {
        lining_pair(walls, state + lining_index * kState);
      } else {
        for (Index j = lining_index; j < lining_index + 2; ++j) {
          double* node = state + j * kState;
          lining(before_row + j * kColumns, node[6], node[7], node[8], node[9]);
        }
      }
      lining_index += 2;
    }
    if (lining_index < n) {
      double* node = state + lining_index * kState;
      lining(before_row + lining_index * kColumns, node[6], node[7], node[8], node[9]);
    }
    for (Index i = 0; i + 1 < n; ++i) {
      y[static_cast<std::size_t>(i)] = 0.5 * (row[i * kColumns] + row[(i + 1) * kColumns]) / kRhoC;
    }
    Index node_index = 0;
    while (node_index < n) {
      Index i = node_index;
      if (i >= 1 && i + 2 < n && row[i * kColumns + 4] > 0.0 && row[(i + 1) * kColumns + 4] > 0.0) {
        junction_pair(row, y.data(), state, i);   // (interior nodes i, i + 1)
        node_index += 2;
        continue;
      }
      node_index += 1;
      double* node = state + i * kState;
      double known = 0.0;
      double admittance = 0.0;
      if (i == 0) {
        known = in.flow_in[step] + 2.0 * y[0] * state[1] - row[1];
        admittance = y[0];
      } else if (i < n - 1) {
        known = 2.0 * y[static_cast<std::size_t>(i - 1)] * state[(i - 1) * kState] +
                2.0 * y[static_cast<std::size_t>(i)] * node[1] - row[i * kColumns + 1];
        admittance = y[static_cast<std::size_t>(i - 1)] + y[static_cast<std::size_t>(i)];
      } else {
        double resistance = 0.0;
        double gain = 0.0;
        radiation(row[(n - 1) * kColumns], in.end_extra, resistance, gain);
        known = 2.0 * y[static_cast<std::size_t>(n - 2)] * state[(n - 2) * kState] - state[11] - gain * state[12] -
                row[(n - 1) * kColumns + 1];
        admittance = y[static_cast<std::size_t>(n - 2)] + 1.0 / resistance + gain;
      }
      node[10] = junction(row + i * kColumns, known, admittance, node[7], node[9]);
    }
    double resistance = 0.0;
    double gain = 0.0;
    radiation(row[(n - 1) * kColumns], in.end_extra, resistance, gain);
    double* last = state + (n - 1) * kState;
    state[11] += gain * (last[10] + state[12]);
    state[12] = last[10];
    in.lips_out[step] = last[10] / resistance + state[11];
    in.port_pressure[step] = state[10];
    Index link_index = 0;
    while (link_index + 2 < n) {
      hop_pair(state, link_index, shelf, in.hop, in.hop_zero);
      link_index += 2;
    }
    for (Index i = link_index; i + 1 < n; ++i) {
      double* link_state = state + i * kState;
      double* next = link_state + kState;
      double wave = link_state[10] - link_state[1];
      link_state[2] = shelf * (wave - in.hop_zero * link_state[4]) + in.hop * link_state[2];
      link_state[4] = wave;
      wave = next[10] - link_state[0];
      link_state[3] = shelf * (wave - in.hop_zero * link_state[5]) + in.hop * link_state[3];
      link_state[5] = wave;
    }
    Index scatter_index = 0;
    while (scatter_index + 1 < n) {
      Index i = scatter_index;
      if (i + 2 < n) {
        double share_first = seal_share(in.geometry, n, step, i, in.seal, in.seal_blend != 0);
        double share_second = seal_share(in.geometry, n, step, i + 1, in.seal, in.seal_blend != 0);
        if (share_first == 1.0 && share_second == 1.0) {
          scatter_pair(row, y.data(), state, i, in.jet_gain > 0.0, in.kinetic != 0);
          scatter_index += 2;
          continue;
        }
      }
      scatter_index += 1;
      double* link_state = state + i * kState;
      double first = row[i * kColumns];
      double second = row[(i + 1) * kColumns];
      double share = seal_share(in.geometry, n, step, i, in.seal, in.seal_blend != 0);
      if (share <= 0.0) {
        link_state[0] = link_state[3];
        link_state[1] = link_state[2];
        link_state[13] = 0.0;
        continue;
      }
      double z = 1.0 / y[static_cast<std::size_t>(i)];
      double drive = 2.0 * (link_state[2] - link_state[3]);
      if (in.jet_gain > 0.0) {
        drive += link_state[18];
      }
      double narrow = first;
      double wide = second;
      if (!(drive >= 0.0)) {
        narrow = second;
        wide = first;
      }
      double k = 0.0;
      if (in.kinetic && wide > narrow) {
        double difference = 1.0 / std::max(narrow, 1e-4) - 1.0 / wide;
        k = 0.5 * 1.14e-3 * difference * difference;
      }
      double linear = 2.0 * z + row[i * kColumns + 2];
      double magnitude = 2.0 * std::abs(drive) / (linear + std::sqrt(linear * linear + 4.0 * k * std::abs(drive)));
      double flow = drive >= 0.0 ? magnitude : -magnitude;
      if (share < 1.0) {
        flow *= share;
      }
      double arriving_next = link_state[3] + z * flow;
      link_state[1] = link_state[2] - z * flow;
      link_state[0] = arriving_next;
      link_state[13] = flow;
    }
  }
}
