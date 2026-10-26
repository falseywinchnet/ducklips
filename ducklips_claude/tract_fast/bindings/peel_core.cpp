// The tube peeled from the sound, time-major: observer.StatesPeel's and
// reverse_core.SidePeel's numbers, to the bit (peel_core.h has the
// schedule). Measurements: NOTES_residual_speed.md.
//
// Pinned arithmetic: compiled with -ffp-contract=off (CMakeLists.txt), and
// every expression below is the numba kernel's, operand for operand and in
// its left-to-right grouping (Python and C++ both associate * / + -
// leftwards). numba fuses no multiply-add and calls the same libm (exp) and
// the same square root instruction; min and max select as Python's do
// (std::min(a, b) is (b < a) ? b : a, Python's min(a, b) keeps a unless b <
// a; likewise max). The constants numba folds (h * h, 0.5 * h, 9 pi^2) are
// the products of the same doubles, folded or not.
//
// Two layers a lane each (layer_pair): a layer-step is a chain of about 100
// cycles (four divisions in series) in about 275 scalar instructions, and
// the processor's window held about two layers' chains at once (MEASURED
// 92 cycles a layer-step, 3.0 instructions a cycle, before the pairing).
// A pair carries two layers' chains in one instruction stream: the same
// IEEE operations per lane (vector divide, multiply, add and subtract round
// as their scalar forms). Where a lane's treatment differs (the start, the
// record, the open end's neighbour, node 0, a node without lining, the
// sound's ends) the pair runs as two scalar layers (layer_scalar).

#include "peel_core.h"
#include "jet_core.hpp"

#include "../include/tract/lanes.hpp"   // <arm_neon.h> on ARM; the same names on SSE2 elsewhere

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

using Index = std::int64_t;

constexpr double kRhoC = 1.14e-3 * 35000.0;    // reverse_core.RHO_C, observer.RHO_C, marching.RHO_C (the same double)
constexpr double kStep = 1.0 / 96000.0;        // the kernels' h
constexpr double kPi = 3.141592653589793;      // numpy.pi
constexpr Index kGeometryColumns = 12;
constexpr Index kStateColumns = 19;            // observer.STATE_COLUMNS
constexpr Index kRing = 4;                     // steps a wave ring holds (a reader needs two behind the writer)
constexpr Index kRingMask = kRing - 1;
// Prefetch distance and stride (MEASURED on the M4, Delta's 39,360 x 43
// geometry: 8 to 64 wavefronts ahead within 6.0-6.6 ns a layer-step, the
// 128-byte stride 0.2-0.5 ns better than 64; NOTES_residual_speed.md).
constexpr Index kPrefetchAhead = 32;           // wavefronts (about 8 us at 0.26 us a wavefront: past memory's latency)
constexpr Index kCacheLine = 128;              // bytes: this processor's cache line

// The layers' fields, each an array over the layers (structure of arrays:
// two neighbouring layers are one vector load).
enum LayerField : Index {
  kSmooth = 0,          // the hop's smoothed wave on link j
  kRawBefore,           // its last unsmoothed wave (STATES' shelf)
  kPreviousFlow,        // link j's flow at its last step
  kFilter0,             // the jet source's shaping filter, four values
  kFilter1,
  kFilter2,
  kFilter3,
  kSource,              // the source link j carried at its last step (STATES)
  kUnsmoothed,          // node j's pressure filter's unsmoothed wave (STATES)
  kCache0,              // the jet's coefficient cache: corner, two gains, root of passed power
  kCache1,
  kCache2,
  kCache3,
  kRightPrevious,       // the wave that entered link j at its last step
  kLeftPrevious,        // the wave it sent back then
  kEnteredNow,          // within a wavefront: the wave that entered link j (pair_link to pair_node)
  kLeftNow,             // and the wave it sent back
  kLayerFields
};

// The packed state of a peel of n nodes (doubles from the start). Layer
// arrays have `stride` entries (n - 1 used); node arrays `stride + 1` (n).
struct StateLayout {
  Index stride;     // n - 1 rounded up to even
  Index layers;     // kLayerFields x stride
  Index lining;     // 4 x (stride + 2): q, v, bq, bv over the nodes
  Index node_stride;
  Index ends;       // 2: the radiation's inertance flow, the last end pressure
  Index arriving;   // kRing x stride: A_j(t), the wave node j + 1 sends into link j, slot t mod kRing
  Index pressure;   // kRing x node_stride: P_i(t), node i's pressure, slot t mod kRing
  Index size;
};

StateLayout layout_for(Index n) {
  StateLayout layout;
  layout.stride = ((n - 1) + 1) / 2 * 2;
  layout.node_stride = layout.stride + 2;
  layout.layers = 0;
  layout.lining = layout.layers + kLayerFields * layout.stride;
  layout.ends = layout.lining + 4 * layout.node_stride;
  layout.arriving = layout.ends + 2;
  layout.pressure = layout.arriving + kRing * layout.stride;
  layout.size = layout.pressure + kRing * layout.node_stride;
  return layout;
}

// The wavefront delay of layer j's link behind the open end.
Index link_delay(Index n, Index j) {
  Index delay = 1 + 2 * (n - 2 - j);
  return delay;
}

// The peel's fixed inputs, copied once a call into this local record, and
// its state's arrays as typed borrows (valid for the call).
struct Peel {
  peel_inputs in;
  Index n;
  Index count;
  StateLayout layout;
  double shelf;     // STATES: (1 - hop) / (1 - hop_zero)
  double* field[kLayerFields];
  double* q;
  double* v;
  double* bq;
  double* bv;
  double* ends;
  double* arriving;
  double* pressure;
  // The jet's narrowest node and its area by geometry row (rows mod
  // jet_mask + 1; reverse_core._narrowest of the row's areas), filled as
  // the open end reaches each row: a link at step t reads row t - 2's
  // (observer.StatesPeel's narrowest[t - 1]), at most 2 n rows behind.
  Index* jet_node;
  double* jet_area;
  Index jet_mask;

  Peel(const peel_inputs& inputs, double* state)
      : in(inputs), n(inputs.nodes), count(inputs.steps), layout(layout_for(inputs.nodes)),
        shelf((1.0 - inputs.hop) / (1.0 - inputs.hop_zero)) {
    for (Index f = 0; f < kLayerFields; ++f) {
      field[f] = state + layout.layers + f * layout.stride;
    }
    q = state + layout.lining;
    v = q + layout.node_stride;
    bq = v + layout.node_stride;
    bv = bq + layout.node_stride;
    ends = state + layout.ends;
    arriving = state + layout.arriving;
    pressure = state + layout.pressure;
    jet_node = nullptr;
    jet_area = nullptr;
    jet_mask = 0;
  }

  const double* row(Index t, Index node) const {
    const double* entry = in.geometry + (t * n + node) * kGeometryColumns;
    return entry;
  }

  // The rings are slotted by the wavefront that wrote the value: A_j(t) by
  // node j + 1 at wavefront t + d_j - 1, P_i(t) by node i at t + d_i + 1
  // (the open end, i = n - 1, d = -1: wavefront t). Its readers come one
  // and two wavefronts later. Within a wavefront every layer's values share
  // a slot, so two neighbouring layers' are adjacent: one vector access.
  double& arriving_at(Index j, Index t) const {
    Index slot = (t + link_delay(n, j) - 1) & kRingMask;
    return arriving[slot * layout.stride + j];
  }
  double& pressure_at(Index i, Index t) const {
    Index slot = (t + link_delay(n, i) + 1) & kRingMask;
    return pressure[slot * layout.node_stride + i];
  }
};

// The helpers below are expanded into the wavefront loop (always_inline):
// a call taking the local Peel's address would let every later store
// through a double pointer alias its fields, which then are read again.

// kernels_jet.jet_source for link j (jet_core.hpp; its filter and cache in
// the layer fields).
[[gnu::always_inline]] inline double jet_source(const Peel& peel, Index link, double previous_flow, Index narrowest,
                                                double area, double gain, Index step) {
  const peel_inputs& in = peel.in;
  jet_core::JetParameters jet;
  jet.bound = in.jet_bound;
  jet.closing = in.jet_closing;
  jet.noise = in.noise;
  jet.noise_links = in.noise_links;
  jet.noise_steps = in.noise_steps;
  jet.seed = in.jet_seed;
  jet_core::JetState held;
  for (Index c = 0; c < 4; ++c) {
    held.filter[c] = peel.field[kFilter0 + c] + link;
    held.cache[c] = peel.field[kCache0 + c] + link;
  }
  double source = jet_core::jet_source(jet, link, previous_flow, narrowest, area, held, gain, step);
  return source;
}

// The jet source link j carries at step t (STATES at t_start: the
// forward state's), or the layer's last one when there is no jet.
[[gnu::always_inline]] inline double link_source(const Peel& peel, int variant, Index j, Index t, Index made) {
  const peel_inputs& in = peel.in;
  double source = peel.field[kSource][j];
  if (in.jet_gain > 0.0 && t > 0) {
    if (variant == PEEL_VARIANT_STATES && t == in.t_start) {
      source = in.init[j * kStateColumns + 18];
    } else {
      double gain = in.jet_gain;
      if (in.has_gate) {
        gain = in.jet_gain * in.gate[std::min(made, in.gate_length - 1)];
      }
      // (observer.StatesPeel's narrowest[made] is the narrowest of row made - 1.)
      Index slot = std::max<Index>(made - 1, 0) & peel.jet_mask;
      source = jet_source(peel, j, peel.field[kPreviousFlow][j], peel.jet_node[slot], peel.jet_area[slot], gain, made);
    }
  }
  return source;
}

// marching._kinetic / reverse_core._kinetic_term (the same expression).
[[gnu::always_inline]] inline double kinetic_term(double flow, double first, double second) {
  double narrow = first;
  double wide = second;
  if (!(flow >= 0.0)) {
    narrow = second;
    wide = first;
  }
  if (wide <= narrow) {
    return 0.0;
  }
  double difference = 1.0 / std::max(narrow, 1e-4) - 1.0 / wide;
  double term = 0.5 * (0.5 * 1.14e-3 * difference * difference) * flow * std::abs(flow);
  return term;
}

// marching._radiation (= reverse_core._radiation_gain): the open end's
// resistance and the trapezoid's gain on its inertance.
[[gnu::always_inline]] inline void radiation(double area, double end_extra, double& resistance, double& gain) {
  resistance = 1.0 * 128.0 * kRhoC / (9.0 * (kPi * kPi) * area);
  double inertance = 1.0 * 8.0 * 1.14e-3 / (3.0 * kPi * std::sqrt(kPi * area)) + 1.14e-3 * end_extra / area;
  gain = 0.5 / 96000.0 / inertance;
}

// reverse_core._lining_row on node i with the geometry row `walls`.
[[gnu::always_inline]] inline void lining_row(const double* walls, Index i, double* q, double* v, double* bq,
                                              double* bv) {
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
  double rhs = mc * (2.0 * q[i] + h * v[i]) + 0.5 * h * h * coupling * bq[i];
  double new_q = rhs / diag - q[i];
  v[i] = 2.0 * (new_q - q[i]) / h - v[i];
  q[i] = new_q;
  double stiffness = coupling + support;
  double pull = coupling * q[i] - stiffness * bq[i] + cc * v[i];
  bv[i] = (mb * bv[i] + h * pull) / (mb + h * (bc + cc) + h * h * stiffness);
  bq[i] += h * bv[i];
}

// reverse_core._middle_row: the cover's velocity stepped on the node's
// pressure (geometry row `walls`); returns its midpoint.
[[gnu::always_inline]] inline double middle_row(const double* walls, Index i, double pressure, double* v,
                                                const double* bv) {
  constexpr double h = kStep;
  double mc = walls[4];
  if (mc <= 0.0) {
    return 0.0;
  }
  double surface = walls[3];
  double cc = walls[9];
  double following = (mc * v[i] + h * surface * pressure - h * cc * (0.5 * v[i] - bv[i])) / (mc + 0.5 * h * cc);
  double middle = 0.5 * (v[i] + following);
  v[i] = following;
  return middle;
}

// reverse_core._narrowest of geometry row t into the jet ring.
[[gnu::always_inline]] inline void jet_row(const Peel& peel, Index t) {
  Index slot = t & peel.jet_mask;
  jet_core::narrowest_of(peel.row(t, 0), kGeometryColumns, peel.n, peel.jet_node[slot], peel.jet_area[slot]);
}

// The open end at step t (observer._states_open_end / reverse_core._side_open_end).
[[gnu::always_inline]] inline void open_end(const Peel& peel, Index t) {
  const peel_inputs& in = peel.in;
  Index n = peel.n;
  double* q = peel.q;
  double* v = peel.v;
  double* bq = peel.bq;
  double* bv = peel.bv;
  double flow_l = peel.ends[0];
  double previous = peel.ends[1];
  double resistance = 0.0;
  double gain = 0.0;
  const double* last = peel.row(t, n - 1);
  radiation(last[0], in.end_extra, resistance, gain);
  double y = 0.5 * (peel.row(t, n - 2)[0] + last[0]) / kRhoC;
  double pressure = (in.sound[t] - flow_l - gain * previous) / (1.0 / resistance + gain);
  lining_row(peel.row(std::max<Index>(t - 1, 0), n - 1), n - 1, q, v, bq, bv);
  double middle = middle_row(last, n - 1, pressure, v, bv);
  double known = pressure * (y + 1.0 / resistance + gain) + last[3] * middle;
  double arriving = (known + flow_l + gain * previous + last[1]) / (2.0 * y);
  flow_l += gain * (pressure + previous);
  previous = pressure;
  peel.pressure_at(n - 1, t) = pressure;
  peel.arriving_at(n - 2, t) = arriving;
  if (in.variant == PEEL_VARIANT_STATES && t == in.record) {
    double* recorded = in.recorded;
    recorded[(n - 1) * kStateColumns + 6] = q[n - 1];
    recorded[(n - 1) * kStateColumns + 7] = v[n - 1];
    recorded[(n - 1) * kStateColumns + 8] = bq[n - 1];
    recorded[(n - 1) * kStateColumns + 9] = bv[n - 1];
    recorded[(n - 1) * kStateColumns + 10] = pressure;
    recorded[11] = flow_l;
    recorded[12] = previous;
  }
  peel.ends[0] = flow_l;
  peel.ends[1] = previous;
}

// Link j at step t: its smoothed wave, flow and jet, the waves that entered
// it and left it (returned in entered_out, left_out).
template <int Variant, bool Jet>
[[gnu::always_inline]] inline void link_step(const Peel& peel, Index j, Index t, double& entered_out,
                                             double& left_out) {
  const peel_inputs& in = peel.in;
  Index made = std::max<Index>(t - 1, 0);
  double smooth = peel.field[kSmooth][j];
  if (Variant == PEEL_VARIANT_STATES) {
    if (!(t == in.t_start && in.t_start > 0)) {
      double raw = 0.0;
      if (t > 0) {
        raw = peel.pressure_at(j + 1, t - 1) - peel.arriving_at(j, t - 1);
      }
      smooth = peel.shelf * (raw - in.hop_zero * peel.field[kRawBefore][j]) + in.hop * smooth;
      peel.field[kRawBefore][j] = raw;
    }
  } else {
    double raw = 0.0;
    if (t > 0) {
      raw = peel.pressure_at(j + 1, t - 1) - peel.arriving_at(j, t - 1);
    }
    smooth = (1.0 - in.hop) * raw + in.hop * smooth;
  }
  const double* here = peel.row(made, j);
  double first = here[0];
  double second = here[kGeometryColumns];
  double z = kRhoC / (0.5 * (first + second));
  double arriving = peel.arriving_at(j, t);
  double flow = (arriving - smooth) / z;
  double entered = smooth + (z + 0.5 * here[2]) * flow;
  if (in.kinetic) {
    entered += kinetic_term(flow, first, second);
  }
  double source = peel.field[kSource][j];
  if (Jet) {
    source = link_source(peel, Variant, j, t, made);
    if (t > 0) {
      entered -= 0.5 * source;
    }
  }
  peel.field[kPreviousFlow][j] = flow;
  peel.field[kSmooth][j] = smooth;
  peel.field[kSource][j] = source;
  double left = entered - z * flow;
  if (Variant == PEEL_VARIANT_STATES && t == in.record + 1) {
    double* recorded = in.recorded + j * kStateColumns;
    recorded[3] = smooth;
    recorded[5] = peel.field[kRawBefore][j];
    recorded[13] = flow;
    recorded[14] = peel.field[kFilter0][j];
    recorded[15] = peel.field[kFilter1][j];
    recorded[16] = peel.field[kFilter2][j];
    recorded[17] = peel.field[kFilter3][j];
    recorded[18] = source;
    if (in.record >= in.t_start && in.record + 1 < peel.count) {
      recorded[0] = arriving;
      recorded[1] = left;
      recorded[2] = entered;
    }
  }
  entered_out = entered;
  left_out = left;
}

// Node j at step t: its pressure from link j's waves at t and t + 1
// (`following`: 0 past the sound), its lining and junction, the wave it
// sends into link j - 1 (node 0: the flow entering and its pressure).
template <int Variant>
[[gnu::always_inline]] inline void node_step(const Peel& peel, Index j, Index t, double following, double entered,
                                             double left) {
  const peel_inputs& in = peel.in;
  double pressure = 0.0;
  if (Variant == PEEL_VARIANT_STATES) {
    double unsmoothed = (following - in.hop * entered) / peel.shelf + in.hop_zero * peel.field[kUnsmoothed][j];
    peel.field[kUnsmoothed][j] = unsmoothed;
    pressure = unsmoothed + left;
    if (t == in.record) {
      in.recorded[j * kStateColumns + 4] = unsmoothed;
    }
  } else {
    pressure = (following - in.hop * entered) / (1.0 - in.hop) + left;
  }
  peel.pressure_at(j, t) = pressure;
  lining_row(peel.row(std::max<Index>(t - 1, 0), j), j, peel.q, peel.v, peel.bq, peel.bv);
  const double* here = peel.row(t, j);
  double middle = middle_row(here, j, pressure, peel.v, peel.bv);
  double right = 0.5 * (here[0] + here[kGeometryColumns]) / kRhoC;
  if (j > 0) {
    double left_y = 0.5 * (here[-kGeometryColumns] + here[0]) / kRhoC;
    double known = 0.0;
    if (Variant == PEEL_VARIANT_STATES) {
      known = pressure * (left_y + right) + here[3] * middle;
    } else {
      double side = 0.0;   // (no side port: SidePeel's port_node -1)
      known = pressure * (left_y + right + side) + here[3] * middle;
    }
    double arriving = (known + here[1] - 2.0 * right * left) / (2.0 * left_y);
    peel.arriving_at(j - 1, t) = arriving;
  } else {
    double flow = pressure * right + here[3] * middle + here[1] - 2.0 * right * left;
    in.flow_out[t] = flow;
    in.pressure_out[t] = pressure;
  }
  if (Variant == PEEL_VARIANT_STATES && t == in.record) {
    double* recorded = in.recorded + j * kStateColumns;
    recorded[6] = peel.q[j];
    recorded[7] = peel.v[j];
    recorded[8] = peel.bq[j];
    recorded[9] = peel.bv[j];
    recorded[10] = pressure;
  }
}

// Layer j at the wavefront where its link is at step t (node j at t - 1),
// each part where it runs.
template <int Variant, bool Jet>
[[gnu::always_inline]] inline void layer_scalar(const Peel& peel, Index j, Index t) {
  Index count = peel.count;
  Index t_start = peel.in.t_start;
  double entered = 0.0;
  double left = 0.0;
  double following = 0.0;
  if (t < count) {
    link_step<Variant, Jet>(peel, j, t, entered, left);
    following = entered;
  }
  Index node_time = t - 1;
  if (node_time >= t_start) {
    node_step<Variant>(peel, j, node_time, following, peel.field[kRightPrevious][j], peel.field[kLeftPrevious][j]);
  }
  if (t < count) {
    peel.field[kRightPrevious][j] = entered;
    peel.field[kLeftPrevious][j] = left;
  }
}

// Two doubles as a vector (lane 0 the first).
[[gnu::always_inline]] inline float64x2_t lanes(double low, double high) {
  float64x2_t pair = vsetq_lane_f64(high, vdupq_n_f64(low), 1);
  return pair;
}

// Whether layers j - 1 (lane 0, link at step t - 2) and j (lane 1, link at
// step t) can run as a pair: both links and nodes run, past the start's and
// the record's steps, both nodes above node 0, both linings present in the
// rows read. (A false here runs them as scalar layers: the same values.)
template <int Variant>
[[gnu::always_inline]] inline bool pair_runs(const Peel& peel, Index j, Index t) {
  const peel_inputs& in = peel.in;
  Index low_t = t - 2;
  if (j < 2 || t >= peel.count || low_t < in.t_start + 2) {
    return false;
  }
  if (Variant == PEEL_VARIANT_STATES) {
    Index record = in.record;
    if (record >= low_t - 2 && record <= t + 1) {
      return false;   // (a link at record + 1 or a node at record in either lane)
    }
  }
  const double* low_now = peel.row(low_t - 1, j - 1);
  const double* low_before = peel.row(low_t - 2, j - 1);
  const double* high_now = peel.row(t - 1, j);
  const double* high_before = peel.row(t - 2, j);
  bool lined = low_now[4] > 0.0 && low_before[4] > 0.0 && high_now[4] > 0.0 && high_before[4] > 0.0;
  return lined;
}

// A lane mask's complement.
[[gnu::always_inline]] inline uint64x2_t mask_not(uint64x2_t mask) {
  uint64x2_t complement = vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(mask)));
  return complement;
}

// The kinetic half-term of two links (kinetic_term per lane): both
// directions' coefficients from the geometry, then the one the flow's sign
// selects (the engine's scatter_links keeps both directions' coefficients
// likewise). Each lane's selected coefficient is kinetic_term's
// expression on the same areas; a lane whose wide side is not wider adds
// +0.0, as kinetic_term's literal. (Both coefficients are computed for
// every lane: the areas are finite and floored, so no exceptional value
// arises in the one discarded; the flow's own chain then waits on no
// division: MEASURED in NOTES_residual_speed.md.)
[[gnu::always_inline]] inline float64x2_t kinetic_pair(float64x2_t flow, float64x2_t first, float64x2_t second) {
  float64x2_t floor = vdupq_n_f64(1e-4);
  float64x2_t one = vdupq_n_f64(1.0);
  float64x2_t coefficient = vdupq_n_f64(0.5 * 1.14e-3);
  // std::max(narrow, 1e-4): (narrow < 1e-4) ? 1e-4 : narrow.
  float64x2_t first_floored = vbslq_f64(vcltq_f64(first, floor), floor, first);
  float64x2_t second_floored = vbslq_f64(vcltq_f64(second, floor), floor, second);
  float64x2_t forward_difference = vsubq_f64(vdivq_f64(one, first_floored), vdivq_f64(one, second));
  float64x2_t backward_difference = vsubq_f64(vdivq_f64(one, second_floored), vdivq_f64(one, first));
  float64x2_t forward_k = vmulq_f64(vmulq_f64(coefficient, forward_difference), forward_difference);
  float64x2_t backward_k = vmulq_f64(vmulq_f64(coefficient, backward_difference), backward_difference);
  // Forward (flow >= 0): narrow first, wide second; it adds where !(wide <= narrow).
  uint64x2_t forward_adds = mask_not(vcleq_f64(second, first));
  uint64x2_t backward_adds = mask_not(vcleq_f64(first, second));
  uint64x2_t forward = vcgeq_f64(flow, vdupq_n_f64(0.0));
  float64x2_t k = vbslq_f64(forward, forward_k, backward_k);
  uint64x2_t adds = vbslq_u64(forward, forward_adds, backward_adds);
  float64x2_t term = vmulq_f64(vmulq_f64(vmulq_f64(vdupq_n_f64(0.5), k), flow), vabsq_f64(flow));
  float64x2_t added = vbslq_f64(adds, term, vdupq_n_f64(0.0));
  return added;
}

// The jet sources of links j - 1 (step t - 2) and j (step t), lanes 0 and
// 1 (pair_runs holds: neither at the start): link_source's values. Where
// both links' coefficient caches hold their corner (the corner sits on its
// floor nearly always), the two filters advance in one vector, each lane
// jet_core::jet_source's expressions; the source's last part per lane in
// scalar (jet_core::jet_strength, its branches kept). Otherwise each link
// runs link_source.
template <int Variant>
[[gnu::always_inline]] inline float64x2_t jet_pair(const Peel& peel, Index j, Index t) {
  const peel_inputs& in = peel.in;
  Index jl = j - 1;
  Index made = t - 1;
  Index made_l = t - 3;
  Index slot = std::max<Index>(made - 1, 0) & peel.jet_mask;
  Index slot_l = std::max<Index>(made_l - 1, 0) & peel.jet_mask;
  float64x2_t area = lanes(peel.jet_area[slot_l], peel.jet_area[slot]);
  float64x2_t velocity = vdivq_f64(vabsq_f64(vld1q_f64(peel.field[kPreviousFlow] + jl)), area);
  if (in.jet_bound > 0.0) {
    float64x2_t bound = vdupq_n_f64(in.jet_bound);
    velocity = vbslq_f64(vcltq_f64(bound, velocity), bound, velocity);   // std::min(velocity, bound)
  }
  float64x2_t slit = vdivq_f64(area, vdupq_n_f64(1.2));
  float64x2_t least_slit = vdupq_n_f64(0.01);
  slit = vbslq_f64(vcltq_f64(slit, least_slit), least_slit, slit);     // std::max(slit, 0.01)
  float64x2_t corner = vdivq_f64(vmulq_f64(vdupq_n_f64(0.2), velocity), slit);
  float64x2_t lowest = vdupq_n_f64(1500.0);
  float64x2_t highest = vdupq_n_f64(12000.0);
  corner = vbslq_f64(vcltq_f64(corner, lowest), lowest, corner);
  corner = vbslq_f64(vcltq_f64(highest, corner), highest, corner);
  uint64x2_t held = vceqq_f64(corner, vld1q_f64(peel.field[kCache0] + jl));
  if (vgetq_lane_u64(held, 0) == 0 || vgetq_lane_u64(held, 1) == 0) {
    double low_source = link_source(peel, Variant, jl, t - 2, made_l);
    double high_source = link_source(peel, Variant, j, t, made);
    return lanes(low_source, high_source);
  }
  float64x2_t corner_gain = vld1q_f64(peel.field[kCache1] + jl);
  float64x2_t floor_gain = vld1q_f64(peel.field[kCache2] + jl);
  float64x2_t root_passed = vld1q_f64(peel.field[kCache3] + jl);
  double white_low = 0.0;
  double white_high = 0.0;
  if (jl < in.noise_links && made_l < in.noise_steps) {
    white_low = in.noise[jl * in.noise_steps + made_l];
  } else {
    white_low = jet_core::hashed_gauss(in.jet_seed, jl, made_l);
  }
  if (j < in.noise_links && made < in.noise_steps) {
    white_high = in.noise[j * in.noise_steps + made];
  } else {
    white_high = jet_core::hashed_gauss(in.jet_seed, j, made);
  }
  float64x2_t white = lanes(white_low, white_high);
  float64x2_t filter0 = vld1q_f64(peel.field[kFilter0] + jl);
  float64x2_t filter1 = vld1q_f64(peel.field[kFilter1] + jl);
  float64x2_t filter2 = vld1q_f64(peel.field[kFilter2] + jl);
  float64x2_t filter3 = vld1q_f64(peel.field[kFilter3] + jl);
  filter0 = vaddq_f64(filter0, vmulq_f64(floor_gain, vsubq_f64(white, filter0)));
  float64x2_t once = vsubq_f64(white, filter0);
  filter1 = vaddq_f64(filter1, vmulq_f64(floor_gain, vsubq_f64(once, filter1)));
  filter2 = vaddq_f64(filter2, vmulq_f64(corner_gain, vsubq_f64(vsubq_f64(once, filter1), filter2)));
  filter3 = vaddq_f64(filter3, vmulq_f64(corner_gain, vsubq_f64(filter2, filter3)));
  vst1q_f64(peel.field[kFilter0] + jl, filter0);
  vst1q_f64(peel.field[kFilter1] + jl, filter1);
  vst1q_f64(peel.field[kFilter2] + jl, filter2);
  vst1q_f64(peel.field[kFilter3] + jl, filter3);
  jet_core::JetParameters jet;
  jet.bound = in.jet_bound;
  jet.closing = in.jet_closing;
  jet.noise = in.noise;
  jet.noise_links = in.noise_links;
  jet.noise_steps = in.noise_steps;
  jet.seed = in.jet_seed;
  double gain_low = in.jet_gain;
  double gain_high = in.jet_gain;
  if (in.has_gate) {
    gain_low = in.jet_gain * in.gate[std::min(made_l, in.gate_length - 1)];
    gain_high = in.jet_gain * in.gate[std::min(made, in.gate_length - 1)];
  }
  double low_source = jet_core::jet_strength(jet, jl, peel.jet_node[slot_l], vgetq_lane_f64(area, 0),
                                             vgetq_lane_f64(velocity, 0), gain_low, vgetq_lane_f64(filter3, 0),
                                             vgetq_lane_f64(root_passed, 0));
  double high_source = jet_core::jet_strength(jet, j, peel.jet_node[slot], vgetq_lane_f64(area, 1),
                                              vgetq_lane_f64(velocity, 1), gain_high, vgetq_lane_f64(filter3, 1),
                                              vgetq_lane_f64(root_passed, 1));
  return lanes(low_source, high_source);
}

// Layers j - 1 and j as lanes 0 and 1 (pair_runs holds), the link half of
// layer_scalar's expressions (link j at step t, j - 1 at t - 2), operation
// for operation on both lanes; the jet per lane in scalar (its branches
// kept); and the nodes' linings (independent of this wavefront's links).
// The waves that entered and left go to the kEnteredNow and kLeftNow
// fields for pair_node.
template <int Variant, bool Jet>
[[gnu::always_inline]] inline void pair_link(const Peel& peel, Index j, Index t) {
  constexpr double h = kStep;
  const peel_inputs& in = peel.in;
  Index jl = j - 1;
  Index tl = t - 2;            // lane 0's link step
  Index made = t - 1;          // lane 1's geometry row (= its node's step)
  Index made_l = tl - 1;
  // The smoothed wave.
  // (P_j(t - 3) and P_{j+1}(t - 1) share a slot; so do A_{j-1}(t - 3) and A_j(t - 1).)
  float64x2_t pressure_above = vld1q_f64(&peel.pressure_at(jl + 1, tl - 1));
  float64x2_t arriving_before = vld1q_f64(&peel.arriving_at(jl, tl - 1));
  float64x2_t raw = vsubq_f64(pressure_above, arriving_before);
  float64x2_t smooth = vld1q_f64(peel.field[kSmooth] + jl);
  float64x2_t hop = vdupq_n_f64(in.hop);
  if (Variant == PEEL_VARIANT_STATES) {
    float64x2_t raw_before = vld1q_f64(peel.field[kRawBefore] + jl);
    float64x2_t shelved = vsubq_f64(raw, vmulq_f64(vdupq_n_f64(in.hop_zero), raw_before));
    smooth = vaddq_f64(vmulq_f64(vdupq_n_f64(peel.shelf), shelved), vmulq_f64(hop, smooth));
    vst1q_f64(peel.field[kRawBefore] + jl, raw);
  } else {
    smooth = vaddq_f64(vmulq_f64(vdupq_n_f64(1.0 - in.hop), raw), vmulq_f64(hop, smooth));
  }
  // The flow and the waves.
  const double* low_row = peel.row(made_l, jl);
  const double* high_row = peel.row(made, j);
  float64x2_t first = lanes(low_row[0], high_row[0]);
  float64x2_t second = lanes(low_row[kGeometryColumns], high_row[kGeometryColumns]);
  float64x2_t z = vdivq_f64(vdupq_n_f64(kRhoC), vmulq_f64(vdupq_n_f64(0.5), vaddq_f64(first, second)));
  float64x2_t arriving = vld1q_f64(&peel.arriving_at(jl, tl));
  float64x2_t flow = vdivq_f64(vsubq_f64(arriving, smooth), z);
  float64x2_t resistance = lanes(low_row[2], high_row[2]);
  float64x2_t series = vaddq_f64(z, vmulq_f64(vdupq_n_f64(0.5), resistance));
  float64x2_t entered = vaddq_f64(smooth, vmulq_f64(series, flow));
  if (in.kinetic) {
    entered = vaddq_f64(entered, kinetic_pair(flow, first, second));
  }
  float64x2_t source = vld1q_f64(peel.field[kSource] + jl);
  if (Jet) {
    source = jet_pair<Variant>(peel, j, t);
    entered = vsubq_f64(entered, vmulq_f64(vdupq_n_f64(0.5), source));
  }
  vst1q_f64(peel.field[kPreviousFlow] + jl, flow);
  vst1q_f64(peel.field[kSmooth] + jl, smooth);
  vst1q_f64(peel.field[kSource] + jl, source);
  float64x2_t left = vsubq_f64(entered, vmulq_f64(z, flow));
  vst1q_f64(peel.field[kEnteredNow] + jl, entered);
  vst1q_f64(peel.field[kLeftNow] + jl, left);
  // The nodes' linings (node j at step t - 1, j - 1 at t - 3), on the
  // geometry of the step before each node's.
  Index node_t = t - 1;
  Index node_tl = tl - 1;
  const double* low_walls = peel.row(node_tl - 1, jl);
  const double* high_walls = peel.row(node_t - 1, j);
  float64x2_t mc = lanes(low_walls[4], high_walls[4]);
  float64x2_t mb = lanes(low_walls[5], high_walls[5]);
  float64x2_t hoop = lanes(low_walls[6], high_walls[6]);
  float64x2_t coupling = lanes(low_walls[7], high_walls[7]);
  float64x2_t support = lanes(low_walls[8], high_walls[8]);
  float64x2_t cc = lanes(low_walls[9], high_walls[9]);
  float64x2_t bc = lanes(low_walls[10], high_walls[10]);
  float64x2_t q = vld1q_f64(peel.q + jl);
  float64x2_t v = vld1q_f64(peel.v + jl);
  float64x2_t bq = vld1q_f64(peel.bq + jl);
  float64x2_t bv = vld1q_f64(peel.bv + jl);
  constexpr double quarter = 0.25 * h * h;
  constexpr double half_h_squared = 0.5 * h * h;
  constexpr double h_squared = h * h;
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
  vst1q_f64(peel.q + jl, q);
  vst1q_f64(peel.v + jl, v);
  vst1q_f64(peel.bq + jl, bq);
  vst1q_f64(peel.bv + jl, bv);
}

// The node half of the pair (after every pair_link of the wavefront):
// node j at step t - 1 and j - 1 at t - 3, their pressures from the waves
// pair_link left, the cover's velocity, the junction and the waves into the
// links below.
template <int Variant>
[[gnu::always_inline]] inline void pair_node(const Peel& peel, Index j, Index t) {
  constexpr double h = kStep;
  const peel_inputs& in = peel.in;
  Index jl = j - 1;
  Index node_t = t - 1;
  Index node_tl = t - 3;
  float64x2_t entered = vld1q_f64(peel.field[kEnteredNow] + jl);
  float64x2_t left = vld1q_f64(peel.field[kLeftNow] + jl);
  float64x2_t entered_before = vld1q_f64(peel.field[kRightPrevious] + jl);
  float64x2_t left_before = vld1q_f64(peel.field[kLeftPrevious] + jl);
  float64x2_t hop = vdupq_n_f64(in.hop);
  float64x2_t lifted = vsubq_f64(entered, vmulq_f64(hop, entered_before));
  float64x2_t pressure;
  if (Variant == PEEL_VARIANT_STATES) {
    float64x2_t unsmoothed_before = vld1q_f64(peel.field[kUnsmoothed] + jl);
    float64x2_t unsmoothed = vaddq_f64(vdivq_f64(lifted, vdupq_n_f64(peel.shelf)),
                                       vmulq_f64(vdupq_n_f64(in.hop_zero), unsmoothed_before));
    vst1q_f64(peel.field[kUnsmoothed] + jl, unsmoothed);
    pressure = vaddq_f64(unsmoothed, left_before);
  } else {
    pressure = vaddq_f64(vdivq_f64(lifted, vdupq_n_f64(1.0 - in.hop)), left_before);
  }
  vst1q_f64(&peel.pressure_at(jl, node_tl), pressure);
  // The cover's velocity on the pressure (the node's own step's row).
  const double* low_here = peel.row(node_tl, jl);
  const double* high_here = peel.row(node_t, j);
  float64x2_t v = vld1q_f64(peel.v + jl);
  float64x2_t bv = vld1q_f64(peel.bv + jl);
  float64x2_t mc_now = lanes(low_here[4], high_here[4]);
  float64x2_t surface = lanes(low_here[3], high_here[3]);
  float64x2_t cc_now = lanes(low_here[9], high_here[9]);
  float64x2_t h_vector = vdupq_n_f64(h);
  float64x2_t pushed = vmulq_f64(vmulq_f64(h_vector, surface), pressure);
  float64x2_t damped = vmulq_f64(vmulq_f64(h_vector, cc_now), vsubq_f64(vmulq_f64(vdupq_n_f64(0.5), v), bv));
  float64x2_t cover_top = vsubq_f64(vaddq_f64(vmulq_f64(mc_now, v), pushed), damped);
  float64x2_t cover_bottom = vaddq_f64(mc_now, vmulq_f64(vdupq_n_f64(0.5 * h), cc_now));
  float64x2_t following = vdivq_f64(cover_top, cover_bottom);
  float64x2_t middle = vmulq_f64(vdupq_n_f64(0.5), vaddq_f64(v, following));
  vst1q_f64(peel.v + jl, following);
  // The junction, the wave into the link below.
  float64x2_t here_area = lanes(low_here[0], high_here[0]);
  float64x2_t above_area = lanes(low_here[kGeometryColumns], high_here[kGeometryColumns]);
  float64x2_t below_area = lanes(low_here[-kGeometryColumns], high_here[-kGeometryColumns]);
  float64x2_t half = vdupq_n_f64(0.5);
  float64x2_t rho_c = vdupq_n_f64(kRhoC);
  float64x2_t right = vdivq_f64(vmulq_f64(half, vaddq_f64(here_area, above_area)), rho_c);
  float64x2_t left_y = vdivq_f64(vmulq_f64(half, vaddq_f64(below_area, here_area)), rho_c);
  float64x2_t admittance = vaddq_f64(left_y, right);
  if (Variant != PEEL_VARIANT_STATES) {
    admittance = vaddq_f64(admittance, vdupq_n_f64(0.0));   // (left + right + side, side 0)
  }
  float64x2_t known = vaddq_f64(vmulq_f64(pressure, admittance), vmulq_f64(surface, middle));
  float64x2_t frame = lanes(low_here[1], high_here[1]);
  float64x2_t two = vdupq_n_f64(2.0);
  float64x2_t sent = vsubq_f64(vaddq_f64(known, frame), vmulq_f64(vmulq_f64(two, right), left_before));
  float64x2_t below_wave = vdivq_f64(sent, vmulq_f64(two, left_y));
  vst1q_f64(&peel.arriving_at(jl - 1, node_tl), below_wave);
  vst1q_f64(peel.field[kRightPrevious] + jl, entered);
  vst1q_f64(peel.field[kLeftPrevious] + jl, left);
}

// Wavefronts begin..target - 1.
template <int Variant, bool Jet>
void advance_variant(const peel_inputs& inputs, double* state, Index begin, Index target, double* snapshots,
                     Index snap_every, Index snap_slots) {
  Peel peel(inputs, state);
  Index n = peel.n;
  Index count = peel.count;
  Index t_start = inputs.t_start;
  Index complete = peel_complete(count, n);
  Index size = peel.layout.size;
  std::vector<Index> pair_heads(static_cast<std::size_t>(n));   // the wavefront's pairs, by their upper layer
  // The jet ring (a power of two above 2 n + 8 rows), its rows before
  // `begin` read from the geometry (the rows the deepest link reads first).
  Index ring_size = 16;
  while (ring_size < 2 * n + 8) {
    ring_size *= 2;
  }
  std::vector<Index> jet_nodes(static_cast<std::size_t>(Jet ? ring_size : 1), -1);
  std::vector<double> jet_areas(static_cast<std::size_t>(Jet ? ring_size : 1), 1.0);
  peel.jet_node = jet_nodes.data();
  peel.jet_area = jet_areas.data();
  peel.jet_mask = ring_size - 1;
  if (Jet) {
    Index first_row = std::max<Index>(begin - ring_size + 1, 0);
    for (Index row = first_row; row < std::min(begin, count); ++row) {
      jet_row(peel, row);
    }
  }
  for (Index wave = begin; wave < target; ++wave) {
    if (snap_every > 0 && wave % snap_every == 0 && wave > begin && wave / snap_every < snap_slots) {
      std::memcpy(snapshots + (wave / snap_every) * size, state, static_cast<std::size_t>(size) * sizeof(double));
    }
    // The geometry row the open end reaches kPrefetchAhead wavefronts on,
    // into the second-level cache: the layers then sweep each row from the
    // last node to the first over 2 n wavefronts (2 n rows, 370 kB for 44
    // nodes, stay there), but as some 90 interleaved strided streams, which
    // the hardware prefetchers do not follow (MEASURED: 5.6 ns a
    // layer-step with the geometry in cache, 10.9 ns from memory, before
    // this prefetch, 6.0 ns after; NOTES_residual_speed.md).
    Index ahead = wave + kPrefetchAhead;
    if (ahead < count) {
      const char* row_bytes = reinterpret_cast<const char*>(peel.row(ahead, 0));
      Index bytes = n * kGeometryColumns * static_cast<Index>(sizeof(double));
      for (Index offset = 0; offset < bytes; offset += kCacheLine) {
        __builtin_prefetch(row_bytes + offset, 0, 2);
      }
    }
    if (Jet && wave < count) {
      jet_row(peel, wave);
    }
    if (wave >= t_start && wave < count) {
      open_end(peel, wave);
    }
    // The links (and linings) of every layer, then the nodes: within a
    // wavefront no layer reads another's output, so the two loops'
    // iterations are independent (each loop's chain half a layer-step's).
    Index heads = 0;
    Index j = n - 2;
    while (j >= 0) {
      Index t = wave - link_delay(n, j);
      if (t < t_start) {
        break;   // (deeper layers are further behind)
      }
      if (t > count) {
        j -= 1;
        continue;
      }
      if (pair_runs<Variant>(peel, j, t)) {
        pair_link<Variant, Jet>(peel, j, t);
        pair_heads[heads] = j;
        heads += 1;
        j -= 2;
      } else {
        layer_scalar<Variant, Jet>(peel, j, t);
        j -= 1;
      }
    }
    for (Index k = 0; k < heads; ++k) {
      Index head = pair_heads[k];
      pair_node<Variant>(peel, head, wave - link_delay(n, head));
    }
  }
  if (snap_every > 0 && target % snap_every == 0 && target > begin && target < complete &&
      target / snap_every < snap_slots) {
    std::memcpy(snapshots + (target / snap_every) * size, state, static_cast<std::size_t>(size) * sizeof(double));
  }
}

}  // namespace

extern "C" {

int64_t peel_state_size(int64_t nodes) {
  StateLayout layout = layout_for(nodes);
  return layout.size;
}

int64_t peel_complete(int64_t steps, int64_t nodes) {
  // Node 0's last step (steps - 1) runs at wavefront steps - 1 + d_0 + 1.
  int64_t complete = steps + link_delay(nodes, 0) + 1;
  return complete;
}

void peel_begin(const peel_inputs* inputs, double* state) {
  const peel_inputs& in = *inputs;
  Index n = in.nodes;
  StateLayout layout = layout_for(n);
  std::fill(state, state + layout.size, 0.0);
  Peel peel(in, state);
  for (Index j = 0; j < layout.stride; ++j) {
    peel.field[kCache0][j] = std::nan("");
    peel.field[kCache1][j] = std::nan("");
    peel.field[kCache2][j] = std::nan("");
    peel.field[kCache3][j] = std::nan("");
  }
  if (in.variant == PEEL_VARIANT_STATES && in.t_start > 0) {
    // observer._states_begin.
    const double* init = in.init;
    for (Index i = 0; i < n; ++i) {
      peel.q[i] = init[i * kStateColumns + 6];
      peel.v[i] = init[i * kStateColumns + 7];
      peel.bq[i] = init[i * kStateColumns + 8];
      peel.bv[i] = init[i * kStateColumns + 9];
    }
    peel.ends[0] = init[11];
    peel.ends[1] = init[12];
    for (Index j = 0; j + 1 < n; ++j) {
      peel.field[kSmooth][j] = init[j * kStateColumns + 3];
      peel.field[kRawBefore][j] = init[j * kStateColumns + 5];
      peel.field[kPreviousFlow][j] = init[j * kStateColumns + 13];
      peel.field[kFilter0][j] = init[j * kStateColumns + 14];
      peel.field[kFilter1][j] = init[j * kStateColumns + 15];
      peel.field[kFilter2][j] = init[j * kStateColumns + 16];
      peel.field[kFilter3][j] = init[j * kStateColumns + 17];
      peel.field[kUnsmoothed][j] = init[j * kStateColumns + 4];
    }
  }
}

void peel_advance(const peel_inputs* inputs, double* state, int64_t begin, int64_t target, double* snapshots,
                  int64_t snap_every, int64_t snap_slots) {
  bool jet = inputs->jet_gain > 0.0;
  if (inputs->variant == PEEL_VARIANT_STATES && jet) {
    advance_variant<PEEL_VARIANT_STATES, true>(*inputs, state, begin, target, snapshots, snap_every, snap_slots);
  } else if (inputs->variant == PEEL_VARIANT_STATES) {
    advance_variant<PEEL_VARIANT_STATES, false>(*inputs, state, begin, target, snapshots, snap_every, snap_slots);
  } else if (jet) {
    advance_variant<PEEL_VARIANT_SIDE, true>(*inputs, state, begin, target, snapshots, snap_every, snap_slots);
  } else {
    advance_variant<PEEL_VARIANT_SIDE, false>(*inputs, state, begin, target, snapshots, snap_every, snap_slots);
  }
}

}  // extern "C"
