#pragma once
#include "ducklips/glottal_storage.hpp"
#include "ducklips/glottal_kinetics.hpp"
#include "ducklips/moving_graph.hpp"

namespace ducklips {
struct GlottalPassages final {
    // Fixed passage coefficients. These do not yet follow the fold geometry.
    std::array<double,3> inertance{};       // Pa s^2/m^3, strictly positive
    std::array<double,3> resistance{};      // Pa s/m^3, nonnegative
    std::array<double,3> quadratic_loss{};  // Pa s^2/m^6, nonnegative
};
struct GraphGlottalStep final {
    std::array<double,3> flow{},law_residual{};
    std::array<double,2> kinetic_mass_effort{};
    double load_effort{},reservoir_effort{},reservoir_work{},network_work{},passage_loss{};
    double energy_change{},energy_residual{}; // chambers+tissue+passages, excluding network storage
    std::size_t iterations{};
};
// Pressure-driven glottal storage connected to a total-flow acoustic graph.
// Three inertive/resistive passages solve their own flows, jointly
// with gas/tissue and the graph's average load. Graph charge is equivalent
// ambient-density mass volume; its pressure is the linearized mass effort.
// Supports fixed controls and mass/geometry-derived open passages.
// Neither mode is a jet-separation or complete closure model.
class GraphGlottalStorage final {
public:
    void prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
        const FoldParameters& tissue,const AirProperties& air,double tension,
        const FoldState& initial,const std::array<double,2>& pressure,const GlottalPassages& passages);
    // Half-cell passage inertia follows stored mass; viscous resistance
    // follows the same rectangular volumes. Positive open chambers only.
    void prepare_articulated(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
        const FoldParameters& tissue,const AirProperties& air,double tension,
        const FoldState& initial,const std::array<double,2>& pressure,double outlet_loss=1.);
    [[nodiscard]] bool advance(std::span<const double> c1,std::span<const double> c2,
        std::span<const double> g1,std::span<const double> g2,double dt,double reservoir_pressure,
        double rest_area,const std::array<double,4>& external_force={},
        std::span<const double> nodal_flow={},std::span<const double> series_pressure={}) noexcept;
    [[nodiscard]] bool source_energy(double& result) const noexcept;
    // Physical states may be restored by the owner. Failure preserves all
    // accepted state and diagnostics; no allocations during advance.
    MovingGraphReference network{};
    GlottalStorage chambers{};
    std::array<double,3> momentum{};
    GraphGlottalStep last{};
private:
    MovingGraphReference free_{},response_{};
    GasStorage gas_{};
    GlottalPassages passages_{};
    GlottalKinetics kinetics_{};
    bool articulated_{};
    std::vector<double> unit_flow_{};
    std::size_t source_node_{};
};
}
