#pragma once

// The moving tube's geometry and its peel, as templated copies of
// analysis/reverse_core.py's _moving_geometry and _peel_moving (no side
// port, no jet, the kinetic term on), for forward-mode derivatives.
//
// Structure used (OBSERVED in reverse_core.py): each node's frame (rest
// radius follower, its pushed air, its lining's constants) depends only on
// that node's own rest-area schedule; a link's viscous resistance on its two
// nodes. A parameter that moves only some nodes' schedules (one node's log
// area in one phone) leaves every other node's geometry exactly unchanged,
// so its tangent is zero there. The peel works from the lips inward one node
// per layer, so such a parameter's tangent is zero in every layer outside
// the node's reach: a block of directions starts at the outermost layer any
// of its nodes reaches, from the primal layer state stored there.
//
// Expressions follow the numba code's order; the translation unit is built
// with -ffp-contract=off (numba does not fuse multiply-adds) so the double
// instantiation can be compared with it bit for bit.

#include "dual.hpp"

#include <cmath>
#include <vector>

namespace tract_ad {

constexpr double kPeelRhoC = 1.14e-3 * 35000.0;

template <class S>
struct LiningRow {
  S surface;
  S cover_mass;
  S body_mass;
  S hoop_k;
  S coupling_k;
  S support_k;
  S cover_c;
  S body_c;
};

// _compile_lining (support factor 1): a node's lining from its radius and
// material row (MATERIALS order: cover density, thickness, hoop tension,
// cover loss, coupling, body density, depth, support, support loss).
template <class S>
LiningRow<S> compile_lining(const S& radius, const double* material) {
  LiningRow<S> row;
  row.surface = 2.0 * 3.14159265358979323846 * radius * (35000.0 / 96000.0);
  row.cover_mass = material[0] * material[1] * row.surface;
  row.body_mass = material[5] * material[6] * row.surface;
  row.hoop_k = material[2] * row.surface / (radius * radius);
  row.coupling_k = material[4] * row.surface;
  row.support_k = material[7] * row.surface;
  S reduced = row.cover_mass * row.body_mass / (row.cover_mass + row.body_mass);
  row.cover_c = 2.0 * material[3] * ad_sqrt(reduced * (row.hoop_k + row.coupling_k));
  row.body_c = 2.0 * material[8] * ad_sqrt(row.body_mass * row.support_k);
  return row;
}

// _viscous_link: a link's viscous resistance for its area.
template <class S>
S viscous_link(const S& area, double scale) {
  S safe = ad_max(area, 1e-4);
  S steep = 0.03 * ad_pow(safe, -2.6);
  S layer = 0.0065 / safe * std::sqrt(2250.0);
  return scale * ad_max(steep, layer) * (35000.0 / 96000.0);
}

// One node's frame over the steps: radius, air area, pushed air (GEODUMP
// columns 11, 0 and 1) and whether it moves fast (the viscous links next to
// it then exact every step).
template <class S>
struct NodeTrack {
  std::vector<S> radius;
  std::vector<S> area;
  std::vector<S> frame;
  std::vector<char> fast;
};

// The rest area a step reads, with its tangent for this block's directions
// (rest_dot: [tick][node][direction], one tick per 96 steps).
template <class S>
struct RestReader;

template <>
struct RestReader<double> {
  static double read(const double* rest, long step, int n, int node, const double*, int, int) {
    return rest[step * n + node];
  }
};

template <int N>
struct RestReader<Dual<N>> {
  static Dual<N> read(const double* rest, long step, int n, int node, const double* rest_dot, int directions,
                      int first) {
    Dual<N> out;
    out.v = rest[step * n + node];
    const double* row = rest_dot + ((step / 96) * n + node) * static_cast<long>(directions);
    for (int k = 0; k < N; ++k) {
      int direction = first + k;
      out.d[k] = direction < directions ? row[direction] : 0.0;
    }
    return out;
  }
};

// _moving_geometry's per-node part for node `node`.
template <class S>
void node_track(long steps, int n, int node, const double* rest, const double* material, const double* rest_dot,
                int directions, int first, NodeTrack<S>& track) {
  const double h = 1.0 / 96000.0;
  const double omega = 2.0 * 3.14159265358979323846 * 160.0;
  const double gain = 1.0 - std::exp(-2.0 * 3.14159265358979323846 * 250.0 * h);
  track.radius.resize(static_cast<std::size_t>(steps));
  track.area.resize(static_cast<std::size_t>(steps));
  track.frame.resize(static_cast<std::size_t>(steps));
  track.fast.resize(static_cast<std::size_t>(steps));
  S first_rest = RestReader<S>::read(rest, 0, n, node, rest_dot, directions, first);
  S target = ad_sqrt(ad_max(first_rest, 1e-4) / 3.14159265358979323846);
  S radius = target;
  S goal = target;
  S ramp = Constant<S>::make(0.0);
  S rate = Constant<S>::make(0.0);
  S low = Constant<S>::make(0.0);
  S lower = Constant<S>::make(0.0);
  bool fast = false;
  int ramp_left = 0;
  for (long step = 0; step < steps; ++step) {
    if (step % 12 == 0) {
      S now = RestReader<S>::read(rest, step, n, node, rest_dot, directions, first);
      goal = ad_sqrt(ad_max(now, 1e-4) / 3.14159265358979323846);
      ramp = (goal - target) / 12.0;
      ramp_left = 12;
    }
    if (ramp_left > 0) {
      target = target + ramp;
      ramp_left -= 1;
      if (ramp_left == 0) {
        target = goal;
      }
    }
    S offset = target - radius;
    if (std::abs(value_of(offset)) >= 1e-9 || std::abs(value_of(rate)) >= 1e-6) {
      rate = rate + h * (omega * omega * offset - 2.0 * omega * rate);
      radius = ad_max(radius + h * rate, 1e-4);
      fast = std::abs(value_of(rate)) > 0.5;
    } else {
      rate = Constant<S>::make(0.0);
      fast = false;
    }
    S surface = 2.0 * 3.14159265358979323846 * radius * (35000.0 / 96000.0);
    S pushed = surface * rate;
    low = low + gain * (pushed - low);
    lower = lower + gain * (low - lower);
    std::size_t index = static_cast<std::size_t>(step);
    track.radius[index] = radius;
    track.frame[index] = lower;
    track.area[index] = ad_max(3.14159265358979323846 * radius * radius, 1e-4);
    track.fast[index] = fast ? 1 : 0;
  }
  (void)material;
}

// _moving_geometry's link part: link between `left` and `right` tracks.
// With `engine_phase`, as the engine's Tube instead (OBSERVED in
// src/tube.cpp: its constructor's refresh_geometry already took one glide
// step and set the countdown to 15, so the engine glides one step ahead of
// reverse_core and refreshes the viscous law at steps 15, 31, ..., where
// reverse_core refreshes at 0, 16, ...).
template <class S>
void link_track(long steps, const NodeTrack<S>& left, const NodeTrack<S>& right, double viscous, std::vector<S>& out,
                bool engine_phase) {
  out.resize(static_cast<std::size_t>(steps));
  S resistance = Constant<S>::make(0.0);
  S resistance_step = Constant<S>::make(0.0);
  int countdown = 0;
  if (engine_phase && steps > 0) {
    S link = 0.5 * (left.area[0] + right.area[0]);
    resistance_step = (viscous_link(link, viscous) - resistance) / 16.0;
    resistance = resistance + resistance_step;
    countdown = 15;
  }
  for (long step = 0; step < steps; ++step) {
    std::size_t index = static_cast<std::size_t>(step);
    S link = 0.5 * (left.area[index] + right.area[index]);
    if (left.fast[index] || right.fast[index]) {
      resistance = viscous_link(link, viscous);
      resistance_step = Constant<S>::make(0.0);
    } else {
      if (countdown == 0) {
        resistance_step = (viscous_link(link, viscous) - resistance) / 16.0;
      }
      resistance = resistance + resistance_step;
    }
    out[index] = resistance;
    countdown = countdown == 0 ? 15 : countdown - 1;
  }
}

// _lining_one for a node with lining `walls` (the step before's). State S,
// lining constants G (double for a node this block does not move: its
// constants then carry no tangent and cost no tangent arithmetic).
template <class S, class G>
void lining_one(const LiningRow<G>& walls, S& q, S& v, S& bq, S& bv) {
  const double h = 1.0 / 96000.0;
  const double quarter = 0.25 * h * h;
  G diag = walls.cover_mass + quarter * (walls.hoop_k + walls.coupling_k);
  S rhs = walls.cover_mass * (2.0 * q + h * v) + 0.5 * h * h * walls.coupling_k * bq;
  S new_q = rhs / diag - q;
  v = 2.0 * (new_q - q) / h - v;
  q = new_q;
  G stiffness = walls.coupling_k + walls.support_k;
  S pull = walls.coupling_k * q - stiffness * bq + walls.cover_c * v;
  bv = (walls.body_mass * bv + h * pull) / (walls.body_mass + h * (walls.body_c + walls.cover_c) + h * h * stiffness);
  bq = bq + h * bv;
}

// _wall_middle: with the node's pressure known, the cover's step and its
// mean velocity over the step (v updated).
template <class S, class G>
S wall_middle(const LiningRow<G>& walls, const S& pressure, S& v, const S& bv) {
  const double h = 1.0 / 96000.0;
  S following = (walls.cover_mass * v + h * walls.surface * pressure - h * walls.cover_c * (0.5 * v - bv)) /
                (walls.cover_mass + 0.5 * h * walls.cover_c);
  S middle = 0.5 * (v + following);
  v = following;
  return middle;
}

// _kinetic_term: the Borda-Carnot half-term of a link (flow S, areas G).
template <class S, class G>
S kinetic_term(const S& flow, const G& first_area, const G& second_area) {
  const G& narrow = value_of(flow) >= 0.0 ? first_area : second_area;
  const G& wide = value_of(flow) >= 0.0 ? second_area : first_area;
  if (value_of(wide) <= value_of(narrow)) {
    return Constant<S>::make(0.0);
  }
  G difference = 1.0 / ad_max(narrow, 1e-4) - 1.0 / wide;
  return 0.5 * (0.5 * 1.14e-3 * difference * difference) * flow * ad_abs(flow);
}

// The peel's arrays at the start of a layer (the layer above's output).
struct LayerCheckpoints {
  std::vector<double> arriving;  // [layer][t]
  std::vector<double> pressure;
};

// Source of a node's track for the peel: computed for this block, or the
// primal one promoted (a node this block's directions never move).
template <class S>
struct TrackSlot {
  int node = -1;
  NodeTrack<S> track;
};

}  // namespace tract_ad
