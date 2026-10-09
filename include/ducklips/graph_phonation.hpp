#pragma once
#include "ducklips/graph_orifice.hpp"
#include "ducklips/sources.hpp"

namespace ducklips {
// Implicit pressure feedback between the existing two-mass folds and the
// total-flow graph port. A scalar solve makes the midpoint aperture agree
// with the pressure drop across the glottal orifice, excluding inlet inertia.
// Fold integration retains the existing symplectic Euler mechanics. This is
// not a new calibrated fold model or a joint mechanical energy theorem.
class GraphPhonation final {
public:
    void prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,const FoldParameters& parameters);
    [[nodiscard]] bool advance(std::span<const double> c1,std::span<const double> c2,
        std::span<const double> g1,std::span<const double> g2,double dt,
        double rest_area,double tension,double reservoir_pressure,double inlet_inertance,
        double surface_pressure,const AirProperties& air,
        std::span<const double> nodal_flow={},std::span<const double> series_pressure={}) noexcept;
    // All accepted state and diagnostics remain unchanged on rejection.
    GraphOrifice source{};
    TwoMassGlottis folds{};
    double midpoint_area{},transglottal_pressure{},coupling_residual{};
    std::size_t iterations{};
private:
    GraphOrifice trial_{};
};
}
