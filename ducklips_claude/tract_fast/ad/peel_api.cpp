// C entry points (ctypes) for the moving tube's geometry and peel with
// forward-mode tangents (analysis/gradients.py). Built with
// -ffp-contract=off, as numba compiles reverse_core.py.
//
// Arrays are row-major float64; tangent arrays [index][direction].

#include "geometry_peel.hpp"

#include <algorithm>
#include <cstring>
#include <thread>
#include <vector>

namespace {

using tract_ad::compile_lining;
using tract_ad::Constant;
using tract_ad::Dual;
using tract_ad::kinetic_term;
using tract_ad::LayerCheckpoints;
using tract_ad::lining_one;
using tract_ad::LiningRow;
using tract_ad::link_track;
using tract_ad::node_track;
using tract_ad::NodeTrack;
using tract_ad::value_of;
using tract_ad::wall_middle;

struct PeelInputs {
  long count = 0;        // steps of sound peeled
  long steps_total = 0;  // steps of geometry (the folds' render may run past the sound)
  int n = 0;             // nodes
  const double* rest = nullptr;       // [steps_total][n], the AREASEQ schedule (keyframe_schedule)
  const double* materials = nullptr;  // [n][9] (material_table)
  double viscous = 1.0;
  double hop = 0.0;
  double end_extra = 0.0;
  const double* flow_out = nullptr;   // [count]
  bool engine_phase = false;          // link resistances glided as the engine's (link_track), not reverse_core's
};

// Every node's primal frame and every link's primal resistance.
struct PrimalGeometry {
  std::vector<NodeTrack<double>> nodes;
  std::vector<std::vector<double>> links;
};

// The node tracks a layer reads: this block's own for a node its directions
// move, the primal one for any other. Four slots, node mod 4: a layer reads
// nodes j - 1 .. j + 1 and the layers descend.
template <class S>
class TrackCache {
 public:
  TrackCache(const PeelInputs& inputs, const PrimalGeometry& primal, const std::vector<char>& moved,
             const double* rest_dot, int directions, int first)
      : inputs_(inputs), primal_(primal), moved_(moved), rest_dot_(rest_dot), directions_(directions), first_(first) {
    for (int slot = 0; slot < 4; ++slot) {
      slot_node_[slot] = -1;
    }
  }

  bool moved(int node) const { return moved_[static_cast<std::size_t>(node)] != 0; }

  void ensure(int node) {
    if (!moved(node)) {
      return;
    }
    int slot = node % 4;
    if (slot_node_[slot] == node) {
      return;
    }
    node_track<S>(inputs_.steps_total, inputs_.n, node, inputs_.rest, inputs_.materials + 9 * node, rest_dot_,
                  directions_, first_, slots_[slot]);
    slot_node_[slot] = node;
  }

  // Link j's resistance track (call after ensure(j), ensure(j + 1)).
  void ensure_link(int link) {
    link_ = link;
    if (!moved(link) && !moved(link + 1)) {
      link_moved_ = false;
      return;
    }
    link_moved_ = true;
    NodeTrack<S> promoted_left;
    NodeTrack<S> promoted_right;
    const NodeTrack<S>& left = track_for(link, promoted_left);
    const NodeTrack<S>& right = track_for(link + 1, promoted_right);
    link_track<S>(inputs_.steps_total, left, right, inputs_.viscous, link_resistance_, inputs_.engine_phase);
  }

  S area(int node, long t) const {
    if (moved(node)) {
      return slots_[node % 4].area[static_cast<std::size_t>(t)];
    }
    return Constant<S>::make(primal_area(node, t));
  }

  S radius(int node, long t) const {
    if (moved(node)) {
      return slots_[node % 4].radius[static_cast<std::size_t>(t)];
    }
    return Constant<S>::make(primal_radius(node, t));
  }

  S frame(int node, long t) const {
    if (moved(node)) {
      return slots_[node % 4].frame[static_cast<std::size_t>(t)];
    }
    return Constant<S>::make(primal_frame(node, t));
  }

  S resistance(long t) const {
    if (link_moved_) {
      return link_resistance_[static_cast<std::size_t>(t)];
    }
    return Constant<S>::make(primal_resistance(t));
  }

  double primal_area(int node, long t) const {
    return primal_.nodes[static_cast<std::size_t>(node)].area[static_cast<std::size_t>(t)];
  }
  double primal_radius(int node, long t) const {
    return primal_.nodes[static_cast<std::size_t>(node)].radius[static_cast<std::size_t>(t)];
  }
  double primal_frame(int node, long t) const {
    return primal_.nodes[static_cast<std::size_t>(node)].frame[static_cast<std::size_t>(t)];
  }
  double primal_resistance(long t) const {
    return primal_.links[static_cast<std::size_t>(link_)][static_cast<std::size_t>(t)];
  }

 private:
  // A node's track as S for link_track: the cached one, or the primal one
  // promoted into `promoted`.
  const NodeTrack<S>& track_for(int node, NodeTrack<S>& promoted) const {
    if (moved(node)) {
      return slots_[node % 4];
    }
    const NodeTrack<double>& source = primal_.nodes[static_cast<std::size_t>(node)];
    std::size_t steps = source.area.size();
    promoted.area.resize(steps);
    promoted.fast = source.fast;
    for (std::size_t step = 0; step < steps; ++step) {
      promoted.area[step] = Constant<S>::make(source.area[step]);
    }
    return promoted;
  }

  const PeelInputs& inputs_;
  const PrimalGeometry& primal_;
  const std::vector<char>& moved_;
  const double* rest_dot_;
  int directions_;
  int first_;
  NodeTrack<S> slots_[4];
  int slot_node_[4];
  int link_ = 0;
  bool link_moved_ = false;
  std::vector<S> link_resistance_;
};

// A layer's geometry reads: G = double where none of the layer's nodes is
// moved by this block (no tangent arithmetic on the tube's constants), G = S
// where one is.
template <class S, bool kMoved>
struct GeometryRead;

template <class S>
struct GeometryRead<S, false> {
  using Type = double;
  static double area(const TrackCache<S>& cache, int node, long t) { return cache.primal_area(node, t); }
  static double radius(const TrackCache<S>& cache, int node, long t) { return cache.primal_radius(node, t); }
  static double frame(const TrackCache<S>& cache, int node, long t) { return cache.primal_frame(node, t); }
  static double resistance(const TrackCache<S>& cache, long t) { return cache.primal_resistance(t); }
};

template <class S>
struct GeometryRead<S, true> {
  using Type = S;
  static S area(const TrackCache<S>& cache, int node, long t) { return cache.area(node, t); }
  static S radius(const TrackCache<S>& cache, int node, long t) { return cache.radius(node, t); }
  static S frame(const TrackCache<S>& cache, int node, long t) { return cache.frame(node, t); }
  static S resistance(const TrackCache<S>& cache, long t) { return cache.resistance(t); }
};

// _radiation_gain (scales 1, 1): the open end's resistance and the
// trapezoidal gain of its inertance.
template <class G>
void radiation_gain(const G& area, double end_extra, G& resistance, G& gain) {
  const double pi = 3.14159265358979323846;
  const double rho_c = 1.14e-3 * 35000.0;
  resistance = 1.0 * 128.0 * rho_c / (9.0 * (pi * pi) * area);
  G inertance = 1.0 * 8.0 * 1.14e-3 / (3.0 * pi * tract_ad::ad_sqrt(pi * area)) + 1.14e-3 * end_extra / area;
  gain = 0.5 / 96000.0 / inertance;
}

// The peel's arrays for one block.
template <class S>
struct PeelArrays {
  std::vector<S> arriving;
  std::vector<S> pressure;
  std::vector<S> rightward;
  std::vector<S> leftward;
};

// The lips: the pressure from the load, the lining with it, the wave that
// arrived (the first loop of _peel_moving_side). From step t0 on (t0 > 0:
// a block whose tangents are zero before t0): the state at t0 replayed in
// double on the primal geometry and the stored primal outputs.
template <class S, bool kMoved>
void lips_layer(const PeelInputs& inputs, const TrackCache<S>& geometry, const LayerCheckpoints* outputs, long t0,
                PeelArrays<S>& arrays) {
  using G = typename GeometryRead<S, kMoved>::Type;
  const long count = inputs.count;
  const int n = inputs.n;
  const double rho_c = 1.14e-3 * 35000.0;
  const double* material = inputs.materials + 9 * (n - 1);
  S q = Constant<S>::make(0.0);
  S v = Constant<S>::make(0.0);
  S bq = Constant<S>::make(0.0);
  S bv = Constant<S>::make(0.0);
  S flow_l = Constant<S>::make(0.0);
  S previous = Constant<S>::make(0.0);
  if (t0 > 0) {
    const double* stored_arriving = outputs->arriving.data() + static_cast<std::size_t>(n - 1) * count;
    const double* stored_pressure = outputs->pressure.data() + static_cast<std::size_t>(n - 1) * count;
    double replay_q = 0.0;
    double replay_v = 0.0;
    double replay_bq = 0.0;
    double replay_bv = 0.0;
    double replay_flow = 0.0;
    double replay_previous = 0.0;
    LiningRow<double> replay_before = compile_lining(geometry.primal_radius(n - 1, 0), material);
    for (long t = 0; t < t0; ++t) {
      double resistance = 0.0;
      double gain = 0.0;
      radiation_gain(geometry.primal_area(n - 1, t), inputs.end_extra, resistance, gain);
      double pressure = stored_pressure[t];
      LiningRow<double> now = compile_lining(geometry.primal_radius(n - 1, t), material);
      lining_one(t == 0 ? now : replay_before, replay_q, replay_v, replay_bq, replay_bv);
      wall_middle(now, pressure, replay_v, replay_bv);
      replay_flow = replay_flow + gain * (pressure + replay_previous);
      replay_previous = pressure;
      replay_before = now;
      arrays.pressure[static_cast<std::size_t>(t)] = Constant<S>::make(pressure);
      arrays.arriving[static_cast<std::size_t>(t)] = Constant<S>::make(stored_arriving[t]);
    }
    q = Constant<S>::make(replay_q);
    v = Constant<S>::make(replay_v);
    bq = Constant<S>::make(replay_bq);
    bv = Constant<S>::make(replay_bv);
    flow_l = Constant<S>::make(replay_flow);
    previous = Constant<S>::make(replay_previous);
  }
  LiningRow<G> before = compile_lining(GeometryRead<S, kMoved>::radius(geometry, n - 1, t0 > 0 ? t0 - 1 : 0), material);
  for (long t = t0; t < count; ++t) {
    std::size_t index = static_cast<std::size_t>(t);
    G resistance;
    G gain;
    radiation_gain(GeometryRead<S, kMoved>::area(geometry, n - 1, t), inputs.end_extra, resistance, gain);
    G y = 0.5 * (GeometryRead<S, kMoved>::area(geometry, n - 2, t) + GeometryRead<S, kMoved>::area(geometry, n - 1, t)) /
          rho_c;
    arrays.pressure[index] = (inputs.flow_out[t] - flow_l - gain * previous) / (1.0 / resistance + gain);
    LiningRow<G> now = compile_lining(GeometryRead<S, kMoved>::radius(geometry, n - 1, t), material);
    lining_one(t == 0 ? now : before, q, v, bq, bv);
    S middle = wall_middle(now, arrays.pressure[index], v, bv);
    S known = arrays.pressure[index] * (y + 1.0 / resistance + gain) + now.surface * middle;
    arrays.arriving[index] = (known + flow_l + gain * previous + GeometryRead<S, kMoved>::frame(geometry, n - 1, t)) /
                             (2.0 * y);
    flow_l = flow_l + gain * (arrays.pressure[index] + previous);
    previous = arrays.pressure[index];
    before = now;
  }
}

// Layer j (link j and node j) of _peel_moving_side, from step t0 on (as
// lips_layer: the arrays hold the layer above's output; for t < t0 this
// layer's stored primal output is copied, and its running state replayed).
template <class S, bool kMoved>
void inner_layer(const PeelInputs& inputs, const TrackCache<S>& geometry, int j, const LayerCheckpoints* outputs,
                 long t0, PeelArrays<S>& arrays) {
  using G = typename GeometryRead<S, kMoved>::Type;
  const long count = inputs.count;
  const double hop = inputs.hop;
  const double rho_c = 1.14e-3 * 35000.0;
  std::vector<S>& arriving = arrays.arriving;
  std::vector<S>& pressure = arrays.pressure;
  std::vector<S>& rightward = arrays.rightward;
  std::vector<S>& leftward = arrays.leftward;
  S smooth = Constant<S>::make(0.0);
  if (t0 > 0) {
    double replay = 0.0;
    for (long t = 0; t < t0; ++t) {
      double raw = t > 0 ? value_of(pressure[static_cast<std::size_t>(t - 1)]) -
                               value_of(arriving[static_cast<std::size_t>(t - 1)])
                         : 0.0;
      replay = (1.0 - hop) * raw + hop * replay;
    }
    smooth = Constant<S>::make(replay);
  }
  for (long t = t0; t < count; ++t) {
    std::size_t index = static_cast<std::size_t>(t);
    S raw = t > 0 ? pressure[index - 1] - arriving[index - 1] : Constant<S>::make(0.0);
    smooth = (1.0 - hop) * raw + hop * smooth;
    long made = std::max(t - 1, 0L);
    G first = GeometryRead<S, kMoved>::area(geometry, j, made);
    G second = GeometryRead<S, kMoved>::area(geometry, j + 1, made);
    G z = rho_c / (0.5 * (first + second));
    S flow = (arriving[index] - smooth) / z;
    S entered = smooth + (z + 0.5 * GeometryRead<S, kMoved>::resistance(geometry, made)) * flow;
    entered = entered + kinetic_term(flow, first, second);
    rightward[index] = entered;
    leftward[index] = entered - z * flow;
  }
  for (long t = t0; t < count; ++t) {
    std::size_t index = static_cast<std::size_t>(t);
    S following = t + 1 < count ? rightward[index + 1] : Constant<S>::make(0.0);
    pressure[index] = (following - hop * rightward[index]) / (1.0 - hop) + leftward[index];
  }
  const double* material = inputs.materials + 9 * j;
  S q = Constant<S>::make(0.0);
  S v = Constant<S>::make(0.0);
  S bq = Constant<S>::make(0.0);
  S bv = Constant<S>::make(0.0);
  if (t0 > 0) {
    const double* stored_arriving = outputs->arriving.data() + static_cast<std::size_t>(j) * count;
    const double* stored_pressure = outputs->pressure.data() + static_cast<std::size_t>(j) * count;
    double replay_q = 0.0;
    double replay_v = 0.0;
    double replay_bq = 0.0;
    double replay_bv = 0.0;
    LiningRow<double> replay_before = compile_lining(geometry.primal_radius(j, 0), material);
    for (long t = 0; t < t0; ++t) {
      LiningRow<double> now = compile_lining(geometry.primal_radius(j, t), material);
      lining_one(t == 0 ? now : replay_before, replay_q, replay_v, replay_bq, replay_bv);
      wall_middle(now, stored_pressure[t], replay_v, replay_bv);
      replay_before = now;
      pressure[static_cast<std::size_t>(t)] = Constant<S>::make(stored_pressure[t]);
      arriving[static_cast<std::size_t>(t)] = Constant<S>::make(stored_arriving[t]);
    }
    q = Constant<S>::make(replay_q);
    v = Constant<S>::make(replay_v);
    bq = Constant<S>::make(replay_bq);
    bv = Constant<S>::make(replay_bv);
  }
  LiningRow<G> before = compile_lining(GeometryRead<S, kMoved>::radius(geometry, j, t0 > 0 ? t0 - 1 : 0), material);
  for (long t = t0; t < count; ++t) {
    std::size_t index = static_cast<std::size_t>(t);
    LiningRow<G> now = compile_lining(GeometryRead<S, kMoved>::radius(geometry, j, t), material);
    lining_one(t == 0 ? now : before, q, v, bq, bv);
    S middle = wall_middle(now, pressure[index], v, bv);
    G right = 0.5 * (GeometryRead<S, kMoved>::area(geometry, j, t) + GeometryRead<S, kMoved>::area(geometry, j + 1, t)) /
              rho_c;
    if (j > 0) {
      G left = 0.5 * (GeometryRead<S, kMoved>::area(geometry, j - 1, t) + GeometryRead<S, kMoved>::area(geometry, j, t)) /
               rho_c;
      S known = pressure[index] * (left + right + 0.0) + now.surface * middle;
      arriving[index] = (known + GeometryRead<S, kMoved>::frame(geometry, j, t) - 2.0 * right * leftward[index]) /
                        (2.0 * left);
    } else {
      arriving[index] = pressure[index] * right + now.surface * middle + GeometryRead<S, kMoved>::frame(geometry, 0, t) -
                        2.0 * right * leftward[index];
    }
    before = now;
  }
}

// The peel (_peel_moving, no side, no jet) from `start_layer` (n - 1: the
// lips; j < n - 1: from the stored output of layer j + 1) to the glottis,
// each layer from its own first step (`first_step` less one step a layer:
// a layer reads one step ahead). With `save`, each layer's output is stored
// (the primal pass: start_layer n - 1, first_step 0).
template <class S>
void peel_layers(const PeelInputs& inputs, TrackCache<S>& geometry, int start_layer, long first_step,
                 LayerCheckpoints* save, const LayerCheckpoints* load, PeelArrays<S>& arrays) {
  const long count = inputs.count;
  const int n = inputs.n;
  arrays.arriving.resize(static_cast<std::size_t>(count));
  arrays.pressure.resize(static_cast<std::size_t>(count));
  arrays.rightward.resize(static_cast<std::size_t>(count));
  arrays.leftward.resize(static_cast<std::size_t>(count));
  int first_layer = start_layer;
  long layers_done = 0;
  if (start_layer == n - 1) {
    geometry.ensure(n - 1);
    geometry.ensure(n - 2);
    long t0 = std::min(std::max(first_step - 4, 0L), count);
    if (geometry.moved(n - 1) || geometry.moved(n - 2)) {
      lips_layer<S, true>(inputs, geometry, load, t0, arrays);
    } else {
      lips_layer<S, false>(inputs, geometry, load, t0, arrays);
    }
    if (save != nullptr) {
      for (long t = 0; t < count; ++t) {
        save->arriving[static_cast<std::size_t>(n - 1) * count + t] = value_of(arrays.arriving[static_cast<std::size_t>(t)]);
        save->pressure[static_cast<std::size_t>(n - 1) * count + t] = value_of(arrays.pressure[static_cast<std::size_t>(t)]);
      }
    }
    first_layer = n - 2;
    layers_done = 1;
  } else {
    // From the stored primal output of the layer above (zero tangents).
    const double* stored_arriving = load->arriving.data() + static_cast<std::size_t>(start_layer + 1) * count;
    const double* stored_pressure = load->pressure.data() + static_cast<std::size_t>(start_layer + 1) * count;
    for (long t = 0; t < count; ++t) {
      arrays.arriving[static_cast<std::size_t>(t)] = Constant<S>::make(stored_arriving[t]);
      arrays.pressure[static_cast<std::size_t>(t)] = Constant<S>::make(stored_pressure[t]);
    }
  }
  for (int j = first_layer; j >= 0; --j) {
    geometry.ensure(j + 1);
    geometry.ensure(j);
    if (j > 0) {
      geometry.ensure(j - 1);
    }
    geometry.ensure_link(j);
    long t0 = std::min(std::max(first_step - 4 - layers_done, 0L), count);
    bool moved = geometry.moved(j) || geometry.moved(j + 1) || (j > 0 && geometry.moved(j - 1));
    if (moved) {
      inner_layer<S, true>(inputs, geometry, j, load, t0, arrays);
    } else {
      inner_layer<S, false>(inputs, geometry, j, load, t0, arrays);
    }
    ++layers_done;
    if (save != nullptr) {
      for (long t = 0; t < count; ++t) {
        save->arriving[static_cast<std::size_t>(j) * count + t] = value_of(arrays.arriving[static_cast<std::size_t>(t)]);
        save->pressure[static_cast<std::size_t>(j) * count + t] = value_of(arrays.pressure[static_cast<std::size_t>(t)]);
      }
    }
  }
}

PrimalGeometry primal_geometry(const PeelInputs& inputs) {
  PrimalGeometry primal;
  primal.nodes.resize(static_cast<std::size_t>(inputs.n));
  for (int node = 0; node < inputs.n; ++node) {
    node_track<double>(inputs.steps_total, inputs.n, node, inputs.rest, inputs.materials + 9 * node, nullptr, 0, 0,
                       primal.nodes[static_cast<std::size_t>(node)]);
  }
  primal.links.resize(static_cast<std::size_t>(inputs.n - 1));
  for (int link = 0; link + 1 < inputs.n; ++link) {
    link_track<double>(inputs.steps_total, primal.nodes[static_cast<std::size_t>(link)],
                       primal.nodes[static_cast<std::size_t>(link + 1)], inputs.viscous,
                       primal.links[static_cast<std::size_t>(link)], inputs.engine_phase);
  }
  return primal;
}

// A block's directions first .. first + width - 1: the nodes they move (any
// nonzero rest tangent), the outermost layer that reads one of them, and
// the first step their geometry moves (its first nonzero tick's).
void block_reach(const PeelInputs& inputs, const double* rest_dot, long ticks, int directions, int first, int width,
                 std::vector<char>& moved, int& start_layer, long& first_step) {
  int n = inputs.n;
  moved.assign(static_cast<std::size_t>(n), 0);
  start_layer = -1;
  first_step = inputs.steps_total;
  for (long tick = 0; tick < ticks; ++tick) {
    for (int node = 0; node < n; ++node) {
      const double* row = rest_dot + (tick * n + node) * static_cast<long>(directions) + first;
      for (int k = 0; k < width; ++k) {
        if (row[k] != 0.0) {
          moved[static_cast<std::size_t>(node)] = 1;
          first_step = std::min(first_step, tick * 96);
        }
      }
    }
  }
  for (int node = 0; node < n; ++node) {
    if (moved[static_cast<std::size_t>(node)]) {
      int layer = node >= n - 2 ? n - 1 : node + 1;
      start_layer = std::max(start_layer, layer);
    }
  }
}

// Everything the tangent blocks read; a thread takes every `stride`-th block
// from `offset`.
struct BlockWork {
  const PeelInputs* inputs = nullptr;
  const PrimalGeometry* primal = nullptr;
  const LayerCheckpoints* outputs = nullptr;
  const double* rest_dot = nullptr;
  long ticks = 0;
  int directions = 0;
  const int* block_starts = nullptr;
  int blocks = 0;
  int offset = 0;
  int stride = 1;
  double* implied_dot = nullptr;
  double* pressure_dot = nullptr;
  double* area0_dot = nullptr;
};

template <int N>
void tangent_block(const BlockWork& work, int first, int width) {
  const PeelInputs& inputs = *work.inputs;
  int directions = work.directions;
  std::vector<char> moved;
  int start_layer = -1;
  long first_step = 0;
  block_reach(inputs, work.rest_dot, work.ticks, directions, first, width, moved, start_layer, first_step);
  if (start_layer < 0) {
    return;  // these directions move no node: their tangents stay 0
  }
  PeelArrays<Dual<N>> arrays;
  TrackCache<Dual<N>> geometry(inputs, *work.primal, moved, work.rest_dot, directions, first);
  peel_layers<Dual<N>>(inputs, geometry, start_layer, first_step, nullptr, work.outputs, arrays);
  for (long t = 0; t < inputs.count; ++t) {
    double* implied_row = work.implied_dot + t * static_cast<long>(directions) + first;
    double* pressure_row = work.pressure_dot + t * static_cast<long>(directions) + first;
    for (int k = 0; k < width; ++k) {
      implied_row[k] = arrays.arriving[static_cast<std::size_t>(t)].d[k];
      pressure_row[k] = arrays.pressure[static_cast<std::size_t>(t)].d[k];
    }
  }
  if (moved[0]) {
    NodeTrack<Dual<N>> glottis;
    node_track<Dual<N>>(inputs.steps_total, inputs.n, 0, inputs.rest, inputs.materials, work.rest_dot, directions, first,
                        glottis);
    for (long step = 0; step < inputs.steps_total; ++step) {
      double* area_row = work.area0_dot + step * static_cast<long>(directions) + first;
      for (int k = 0; k < width; ++k) {
        area_row[k] = glottis.area[static_cast<std::size_t>(step)].d[k];
      }
    }
  }
}

void tangent_blocks(const BlockWork* work) {
  for (int block = work->offset; block < work->blocks; block += work->stride) {
    int first = work->block_starts[block];
    int width = work->block_starts[block + 1] - first;
    if (width <= 0) {
      continue;
    }
    if (width <= 8) {
      tangent_block<8>(*work, first, width);
    } else if (width <= 16) {
      tangent_block<16>(*work, first, width);
    } else {
      tangent_block<32>(*work, first, width);
    }
  }
}

void tangent_threads(const BlockWork& shared, int threads) {
  if (threads <= 1) {
    tangent_blocks(&shared);
    return;
  }
  std::vector<BlockWork> shares(static_cast<std::size_t>(threads), shared);
  std::vector<std::thread> workers;
  for (int thread = 0; thread < threads; ++thread) {
    shares[static_cast<std::size_t>(thread)].offset = thread;
    shares[static_cast<std::size_t>(thread)].stride = threads;
    workers.emplace_back(tangent_blocks, &shares[static_cast<std::size_t>(thread)]);
  }
  for (std::thread& worker : workers) {
    worker.join();
  }
}

}  // namespace

extern "C" {

// The geometry and peel: implied flow entering and the first node's pressure
// (count), the first node's area (steps_total); with directions > 0 and
// rest_dot ([tick][node][direction], ticks = ceil(steps_total / 96)), their
// tangents (zeroed here first), the directions in `blocks` blocks
// (block_starts: blocks + 1 ascending starts, each block at most 32 wide)
// shared among `threads` threads. engine_phase 1: the link resistances
// glided as the engine's (link_track), not reverse_core's.
int sj_peel(long count, long steps_total, int n, const double* rest, const double* materials, double viscous,
            double hop, double end_extra, const double* flow_out, int directions, const double* rest_dot, long ticks,
            const int* block_starts, int blocks, int threads, int engine_phase, double* implied, double* pressure_out,
            double* area0, double* implied_dot, double* pressure_dot, double* area0_dot) {
  for (int block = 0; block < blocks; ++block) {
    if (block_starts[block + 1] - block_starts[block] > 32 || block_starts[block + 1] < block_starts[block]) {
      return 1;
    }
  }
  PeelInputs inputs;
  inputs.count = count;
  inputs.steps_total = steps_total;
  inputs.n = n;
  inputs.rest = rest;
  inputs.materials = materials;
  inputs.viscous = viscous;
  inputs.hop = hop;
  inputs.end_extra = end_extra;
  inputs.flow_out = flow_out;
  inputs.engine_phase = engine_phase != 0;
  PrimalGeometry primal = primal_geometry(inputs);
  LayerCheckpoints outputs;
  outputs.arriving.assign(static_cast<std::size_t>(n) * count, 0.0);
  outputs.pressure.assign(static_cast<std::size_t>(n) * count, 0.0);
  std::vector<char> nothing(static_cast<std::size_t>(n), 0);
  TrackCache<double> geometry(inputs, primal, nothing, nullptr, 0, 0);
  PeelArrays<double> arrays;
  peel_layers<double>(inputs, geometry, n - 1, 0, &outputs, nullptr, arrays);
  std::memcpy(implied, arrays.arriving.data(), sizeof(double) * static_cast<std::size_t>(count));
  std::memcpy(pressure_out, arrays.pressure.data(), sizeof(double) * static_cast<std::size_t>(count));
  for (long step = 0; step < steps_total; ++step) {
    area0[step] = primal.nodes[0].area[static_cast<std::size_t>(step)];
  }
  if (directions <= 0) {
    return 0;
  }
  std::memset(implied_dot, 0, sizeof(double) * static_cast<std::size_t>(count) * directions);
  std::memset(pressure_dot, 0, sizeof(double) * static_cast<std::size_t>(count) * directions);
  std::memset(area0_dot, 0, sizeof(double) * static_cast<std::size_t>(steps_total) * directions);
  BlockWork work;
  work.inputs = &inputs;
  work.primal = &primal;
  work.outputs = &outputs;
  work.rest_dot = rest_dot;
  work.ticks = ticks;
  work.directions = directions;
  work.block_starts = block_starts;
  work.blocks = blocks;
  work.implied_dot = implied_dot;
  work.pressure_dot = pressure_dot;
  work.area0_dot = area0_dot;
  tangent_threads(work, threads);
  return 0;
}

// The primal geometry as GEODUMP's first three columns (area, pushed air,
// link resistance) per [step][node], for comparison with reverse_core.
int sj_geometry(long steps_total, int n, const double* rest, const double* materials, double viscous, int engine_phase,
                double* out) {
  PeelInputs inputs;
  inputs.count = 0;
  inputs.steps_total = steps_total;
  inputs.n = n;
  inputs.rest = rest;
  inputs.materials = materials;
  inputs.viscous = viscous;
  inputs.engine_phase = engine_phase != 0;
  PrimalGeometry primal = primal_geometry(inputs);
  for (long step = 0; step < steps_total; ++step) {
    for (int node = 0; node < n; ++node) {
      double* row = out + (step * n + node) * 3;
      row[0] = primal.nodes[static_cast<std::size_t>(node)].area[static_cast<std::size_t>(step)];
      row[1] = primal.nodes[static_cast<std::size_t>(node)].frame[static_cast<std::size_t>(step)];
      row[2] = node + 1 < n ? primal.links[static_cast<std::size_t>(node)][static_cast<std::size_t>(step)] : 0.0;
    }
  }
  return 0;
}

}  // extern "C"
