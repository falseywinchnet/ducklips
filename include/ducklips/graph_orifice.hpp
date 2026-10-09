#pragma once
#include "ducklips/moving_graph.hpp"
#include "ducklips/physics.hpp"

namespace ducklips {
// Total-flow source port for a moving, fixed-topology pressure network.
// The aperture and reservoir pressure are held for one integration step.
// Solves U against the work-conjugate average load P_free + Z U. Geometry
// displacement belongs in nodal_flow; do not add a separate mean-flow cavity.
// This is a quasi-steady orifice, not a dynamic vocal-fold model.
class GraphOrifice final {
public:
    void prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node);
    [[nodiscard]] bool advance(std::span<const double> c1,std::span<const double> c2,
        std::span<const double> g1,std::span<const double> g2,double dt,
        double area,double reservoir_pressure,double inlet_inertance,const AirProperties& air,
        std::span<const double> nodal_flow={},std::span<const double> series_pressure={}) noexcept;
    // Owner supplies/continues canonical state and handles physical contacts.
    // On rejection these accepted states, averages and source outputs survive.
    MovingGraphReference state{};
    // inlet_inertance is the boundary-to-centre air path, integrated by
    // implicit midpoint. Its momentum is retained across steps. Zero means
    // an ideal nodal port. A sealed aperture removes this passage state.
    double momentum{},flow{},load_pressure{},load_impedance{};
private:
    MovingGraphReference free_{},response_{};
    std::vector<double> unit_flow_{};
    std::size_t source_node_{};
};
}
