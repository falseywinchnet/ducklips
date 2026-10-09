#pragma once
#include "ducklips/acoustic_network.hpp"
#include <array>
#include <span>
#include <vector>

namespace ducklips {
// Fourth-order, two-stage Gauss reference for the explicitly assumed canonical
// lossless network q'=-B G(t) lambda, lambda'=B^T C(t)^-1 q.
// Keeps physical charge and momentum continuous. No modal reprojection.
// Branches may ground a node; other branches must form a tree (parallel edges
// allowed). Fixed topology; positive C, nonnegative G; optional constant
// nodal flow and series pressure over a step.
// Dissipation, when needed, is supplied separately by the owning integrator.
// CanonicalSpeech owns source feedback and losses around this acoustic step.
class MovingTreeReference final {
public:
    void prepare(std::size_t nodes, std::span<const TreeEdge> branches);
    // Stage metrics at t+(1/2 +/- sqrt(3)/6)*dt. No allocation after prepare.
    // All spans must have their prepared sizes; dt finite and positive.
    // The caller supplies finite positive compliances and nonnegative G.
    void advance(std::span<const double> c1, std::span<const double> c2,
                 std::span<const double> g1, std::span<const double> g2, double dt,
                 std::span<const double> nodal_flow = {},
                 std::span<const double> series_pressure = {}) noexcept;
    // Gauss-weighted pressure over the last step; conjugate to held nodal flow.
    std::vector<double> average_pressure{};
    std::vector<double> charge{};
    std::vector<double> momentum{};
private:
    using Matrix = std::array<double,4>;
    using Pair = std::array<double,2>;
    std::vector<TreeEdge> branches_{};
    std::vector<std::size_t> order_{}, parent_{};
    std::vector<Matrix> diagonal_{}, coupling_{}, inverse_{};
    std::vector<Pair> stage_{};
};
}
