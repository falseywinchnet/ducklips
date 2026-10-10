#pragma once
#include "ducklips/acoustic_network.hpp"
#include <array>
#include <span>
#include <vector>
namespace ducklips {
struct RelaxationMetric final {double storage{},rate{};};
struct RelaxationGraphStep final {
    double energy_change{},physical_dissipation{},numerical_dissipation{},external_work{},energy_residual{},air_residual{};
};
// Coupled pressure, transport and passive relaxation. Fixed positive node
// compliance and nonnegative branch inverse inertia during one step.
// Node memory: z'=rate*(storage*P-z), q' includes -z'.
// Branch memory: z'=rate*(storage*Q-z), p' includes -z'.
// Energy of either memory is z*z/(2*storage). (0,0) denotes absent zero state.
// Two-stage Radau IIA; local memory and branch stages are eliminated before
// the 2*nodes pressure solve. No large matrix over all memory coordinates.
// Numerical damping is separately accounted, never called physical loss.
class RelaxationGraphReference final {
public:
    void prepare(std::size_t nodes,std::span<const TreeEdge> branches,
        std::span<const std::size_t> node_owners={},std::span<const std::size_t> branch_owners={});
    [[nodiscard]] bool advance(std::span<const double> compliance,std::span<const double> conductance,
        std::span<const double> resistance,std::span<const RelaxationMetric> node_metrics,
        std::span<const RelaxationMetric> branch_metrics,double dt,
        std::span<const double> nodal_flow={},std::span<const double> series_pressure={}) noexcept;
    std::vector<double> charge,momentum,node_memory,branch_memory,average_pressure,average_flow;
    std::vector<double> node_dissipation,branch_dissipation,resistance_dissipation;
    RelaxationGraphStep last{};
private:
    using Pair=std::array<double,2>;using Block=std::array<double,4>;
    [[nodiscard]] bool solve() noexcept;
    std::size_t nodes_{};
    std::vector<TreeEdge> edges_;
    std::vector<std::size_t> node_owners_,branch_owners_;
    std::vector<double> matrix_,original_,right_,original_right_,pressure_,next_q_,next_p_,next_node_,next_branch_;
    std::vector<double> trial_node_heat_,trial_branch_heat_,trial_resistance_heat_;
    std::vector<Block> node_kernel_,branch_kernel_,admittance_,node_decay_,branch_decay_;
    std::vector<Pair> free_flow_,flow_,dq_,dp_,node_drive_,branch_drive_;
};
}
