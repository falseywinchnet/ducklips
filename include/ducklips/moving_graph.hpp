#pragma once
#include "ducklips/acoustic_network.hpp"
#include <span>
#include <vector>

namespace ducklips {
// General-graph reference for the same canonical equations as MovingTreeReference:
// q' = -B G lambda + nodal_flow, lambda' = B^T C^-1 q + series_pressure.
// Cycles, parallel branches, ground ports and disconnected components retain
// every physical state. Positive C and nonnegative G, fixed prepared topology.
// As in MovingTreeReference, physical contact-state removal and dissipation
// belong to the owner; zero G stops transport but does not erase momentum.
// Dense pivoted solve of the two Gauss-stage pressure equations: intended as
// a correctness/integration reference for small networks, not a fast renderer.
class MovingGraphReference final {
public:
    void prepare(std::size_t nodes,std::span<const TreeEdge> branches);
    // No allocation after preparation. Invalid input or failed linear solve
    // returns false without changing accepted state or average port outputs.
    [[nodiscard]] bool advance(std::span<const double> c1,std::span<const double> c2,
        std::span<const double> g1,std::span<const double> g2,double dt,
        std::span<const double> nodal_flow={},std::span<const double> series_pressure={}) noexcept;
    std::vector<double> charge{},momentum{},average_pressure{},average_flow{};
private:
    bool solve() noexcept;
    std::size_t nodes_{},dimension_{};
    std::vector<TreeEdge> branches_{};
    std::vector<double> matrix_{},original_{},right_{},original_right_{},stage_{};
    std::vector<double> candidate_charge_{},candidate_momentum_{},candidate_pressure_{},candidate_flow_{};
};
}
