#pragma once
#include "ducklips/continuous_glottis.hpp"
#include "ducklips/moving_graph.hpp"

namespace ducklips {
struct GraphContinuousGlottalStep final {
    double free_pressure{},load_slope{},pressure{},flow{};
    double source_port_work{},graph_port_work{},port_residual{};
    double graph_external_work{};
};
// Implicit work-conjugate coupling of the continuous source to a linear
// acoustic graph. Graph pressure is its linearized mass-effort variable;
// the common flow is reference-volume mass flux, not finite-density volume
// flux. No independent waveform or lagged pressure drives the tract.
class GraphContinuousGlottis final {
public:
    void prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
        const FoldParameters&,const AirProperties&,double tension,const FoldState&,
        const std::array<double,2>& effort,const std::array<double,3>& flow={},
        double outlet_loss=1.,ContourPressureProfile profile=ContourPressureProfile::flat_caps);
    // Prescribed stage metrics, constant external drives, one source step.
    // No allocation after prepare. Failure preserves network, source and last.
    // A caller subdividing contact steps must recompute BOTH graph stages and
    // source feedback for each substep. Moving metrics also do geometry work:
    // graph_external_work + graph_port_work alone is not their energy change.
    [[nodiscard]] bool advance(std::span<const double> c1,std::span<const double> c2,
        std::span<const double> g1,std::span<const double> g2,double dt,double reservoir_pressure,
        const std::array<double,4>& external_force={},std::span<const double> nodal_flow={},
        std::span<const double> series_pressure={},
        ContinuousGlottalUpdate update=ContinuousGlottalUpdate::symmetric) noexcept;
    MovingGraphReference network{};
    ContinuousGlottis source{};
    GraphContinuousGlottalStep last{};
private:
    MovingGraphReference free_{},response_{};
    std::vector<double> unit_flow_{};
    std::size_t source_node_{};
};
}
