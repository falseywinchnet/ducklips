// The moving geometry (reverse_core._geometry_steps, Tube::refresh_geometry's
// replica) to the bit, and per step the least area, the first node's area
// and the reversible jet's narrowest node, gathered while each row is
// written (three later passes over the 96-byte rows read, between them, the
// whole geometry again: MEASURED in NOTES_residual_speed.md). The viscous
// law (a pow) is kept by link area, as the engine keeps it.
//
// Pinned arithmetic as peel_core.cpp: -ffp-contract=off, the numba
// kernel's expressions operand for operand, the same libm (exp, pow, sqrt).

#include "peel_core.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

using Index = std::int64_t;

constexpr double kPi = 3.141592653589793;         // numpy.pi
constexpr Index kColumns = 12;                     // a GEODUMP row's values per node
constexpr Index kStateColumns = 10;                // reverse_core.GEOMETRY_STATE
constexpr Index kMaterialColumns = 9;
constexpr Index kFirstJetNode = 10;                // reverse_core._narrowest's first node

// reverse_core._viscous_link.
double viscous_link(double area, double scale) {
  double safe = std::max(area, 1e-4);
  double law = std::max(0.03 * std::pow(safe, -2.6), 0.0065 / safe * std::sqrt(2250.0));
  double resistance = scale * law * (35000.0 / 96000.0);
  return resistance;
}

// reverse_core._compile_lining_at: node i's lining constants for `radius`.
void compile_lining(double radius, const double* material, double* row) {
  double surface = 2.0 * kPi * radius * (35000.0 / 96000.0);
  double cover_mass = material[0] * material[1] * surface;
  double body_mass = material[5] * material[6] * surface;
  double hoop_k = material[2] * surface / (radius * radius);
  double coupling_k = material[4] * surface;
  double support_k = material[7] * surface;
  double reduced = cover_mass * body_mass / (cover_mass + body_mass);
  row[3] = surface;
  row[4] = cover_mass;
  row[5] = body_mass;
  row[6] = hoop_k;
  row[7] = coupling_k;
  row[8] = support_k;
  row[9] = 2.0 * material[3] * std::sqrt(reduced * (hoop_k + coupling_k));
  row[10] = 2.0 * material[8] * std::sqrt(body_mass * support_k);
  row[11] = radius;
}

}  // namespace

extern "C" void geometry_steps(const geometry_inputs* inputs) {
  const geometry_inputs in = *inputs;
  Index n = in.nodes;
  constexpr double h = 1.0 / 96000.0;
  constexpr double omega = 2.0 * kPi * 160.0;
  const double gain = 1.0 - std::exp(-2.0 * kPi * 250.0 * h);
  double* state = in.state;
  // The viscous law by link area (a pure function of it: taken again only
  // when the area changed, as the engine's Tube::link_viscous keeps it).
  std::vector<double> law_area(static_cast<std::size_t>(n), std::nan(""));
  std::vector<double> law_value(static_cast<std::size_t>(n), 0.0);
  Index ramp_left = static_cast<Index>(in.counters[0]);
  Index countdown = static_cast<Index>(in.counters[1]);
  for (Index step = in.begin; step < in.end; ++step) {
    if (in.snap_every > 0 && step % in.snap_every == 0) {
      Index index = step / in.snap_every;
      std::copy(state, state + n * kStateColumns, in.snaps + index * n * kStateColumns);
      in.snap_counters[index * 2] = static_cast<double>(ramp_left);
      in.snap_counters[index * 2 + 1] = static_cast<double>(countdown);
    }
    const double* rest = in.rest + (step / in.rest_every) * n;
    double* row = in.out + step * n * kColumns;
    if (step % 12 == 0) {
      for (Index i = 0; i < n; ++i) {
        double* node = state + i * kStateColumns;
        node[2] = std::sqrt(std::max(rest[i], 1e-4) / kPi);
        node[3] = (node[2] - node[0]) / 12.0;
      }
      ramp_left = 12;
    }
    if (ramp_left > 0) {
      for (Index i = 0; i < n; ++i) {
        double* node = state + i * kStateColumns;
        node[0] += node[3];
      }
      ramp_left -= 1;
      if (ramp_left == 0) {
        for (Index i = 0; i < n; ++i) {
          double* node = state + i * kStateColumns;
          node[0] = node[2];
        }
      }
    }
    for (Index i = 0; i < n; ++i) {
      double* node = state + i * kStateColumns;
      double* entry = row + i * kColumns;
      double offset = node[0] - node[1];
      if (std::abs(offset) >= 1e-9 || std::abs(node[4]) >= 1e-6) {
        node[4] += h * (omega * omega * offset - 2.0 * omega * node[4]);
        node[1] = std::max(node[1] + h * node[4], 1e-4);
        node[7] = std::abs(node[4]) > 0.5 ? 1.0 : 0.0;
      } else {
        node[4] = 0.0;
        node[7] = 0.0;
      }
      if (step > in.begin && node[1] == entry[11 - n * kColumns]) {
        const double* before = entry - n * kColumns;   // (the step before's entry for node i)
        for (Index c = 3; c < kColumns; ++c) {
          entry[c] = before[c];
        }
      } else {
        compile_lining(node[1], in.materials + i * kMaterialColumns, entry);
      }
      double pushed = entry[3] * node[4];
      node[5] += gain * (pushed - node[5]);
      node[6] += gain * (node[5] - node[6]);
      entry[1] = node[6];
      entry[0] = std::max(kPi * node[1] * node[1], 1e-4);
    }
    for (Index i = 0; i + 1 < n; ++i) {
      double* node = state + i * kStateColumns;
      double link = 0.5 * (row[i * kColumns] + row[(i + 1) * kColumns]);
      bool exact = node[7] != 0.0 || state[(i + 1) * kStateColumns + 7] != 0.0;
      if (exact || countdown == 0) {
        if (!(link == law_area[static_cast<std::size_t>(i)])) {
          law_area[static_cast<std::size_t>(i)] = link;
          law_value[static_cast<std::size_t>(i)] = viscous_link(link, in.viscous);
        }
      }
      if (exact) {
        node[8] = law_value[static_cast<std::size_t>(i)];
        node[9] = 0.0;
      } else {
        if (countdown == 0) {
          node[9] = (law_value[static_cast<std::size_t>(i)] - node[8]) / 16.0;
        }
        node[8] += node[9];
      }
      row[i * kColumns + 2] = node[8];
    }
    row[(n - 1) * kColumns + 2] = 0.0;
    countdown = countdown == 0 ? 15 : countdown - 1;
    // The row's least area (observer._least_area: numpy's minimum, NaN once
    // any is NaN), its first area, its narrowest jet node (reverse_core._narrowest).
    double least = row[0];
    for (Index i = 1; i < n; ++i) {
      double value = row[i * kColumns];
      if (value < least || value != value) {
        if (least == least) {
          least = value;
        }
      }
    }
    in.least[step] = least;
    in.first_area[step] = row[0];
    Index narrowest = -1;
    double narrowest_area = 0.6;
    for (Index i = kFirstJetNode; i < n; ++i) {
      double value = row[i * kColumns];
      if (value < narrowest_area) {
        narrowest_area = value;
        narrowest = i;
      }
    }
    in.jet_node[step] = narrowest;
    in.jet_area[step] = std::max(narrowest_area, 1e-4);
  }
  in.counters[0] = static_cast<double>(ramp_left);
  in.counters[1] = static_cast<double>(countdown);
}
