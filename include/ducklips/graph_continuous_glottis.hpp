#pragma once
#include "ducklips/continuous_glottis.hpp"
#include "ducklips/moving_graph.hpp"
#include "ducklips/radiation_boundary.hpp"
#include "ducklips/slit_thermal.hpp"
#include "ducklips/slit_viscous.hpp"
#include "ducklips/relaxation_graph.hpp"

namespace ducklips {
struct RadiationGeometry final {
    std::size_t branch{};
    double first_conductance{},last_conductance{},midpoint_area{},last_area{};
};
struct ThermalSite final {std::size_t node{};double compliance{},half_gap{};};
struct ThermalGeometry final {double compliance{},half_gap{};};
struct ViscousGeometry final {double segment_inertance{},branch_conductance{},half_gap{};};
struct ViscousSite final {std::size_t branch{};ViscousGeometry geometry{};};
struct JointGraphGeometry final {
    std::span<const double> first_compliance,midpoint_compliance,last_compliance;
    std::span<const double> first_conductance,midpoint_conductance,last_conductance;
    std::span<const ThermalGeometry> midpoint_thermal,last_thermal;
    std::span<const ViscousGeometry> midpoint_viscous,last_viscous;
    double midpoint_radiation_area{},last_radiation_area{};
};
struct GraphContinuousGlottalStep final {
    std::size_t source_closures{},source_reopenings{};
    double free_pressure{},load_slope{},pressure{},flow{};
    double source_port_work{},graph_port_work{},port_residual{};
    double graph_external_work{},source_geometry_work{};
    double radiated_energy{},radiation_geometry_work{},radiation_energy_residual{};
    double thermal_dissipation{},thermal_geometry_work{},thermal_energy_residual{},thermal_air_residual{};
    double viscous_dissipation{},viscous_geometry_work{},viscous_energy_residual{};
    double graph_geometry_work{},passive_numerical_dissipation{},passive_energy_change{},passive_energy_residual{},passive_air_residual{};
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
        double outlet_loss=1.,ContourPressureProfile profile=ContourPressureProfile::flat_caps,double outlet_inertance=0.);
    // Selected nodes with an explicit parallel-wall thermal cross section.
    // Initializes zero thermal charge; equilibrate explicitly for an initially
    // isothermal pressure state. Thermal stored_air() belongs in total mass.
    void prepare_thermal(const AirProperties&,std::span<const ThermalSite>);
    // Distinct physical segments may share a branch. Their inertances must
    // sum to no more than that branch's existing plug-flow inertia.
    // Zero-conductance branches may retain finite local profile sites; the
    // tuple (0,0,0) explicitly retires an absent site's memory and energy metric.
    void prepare_viscous(const AirProperties&,std::span<const ViscousSite>);
    // Call after configuring loss sites. All passive memories, pressure and
    // flow then share one solve. Midpoint freezing of prescribed geometry;
    // endpoint metric work and exact absence are separately accounted.
    void prepare_joint(std::size_t radiation_branch=ground_node,RelaxationOrder order=RelaxationOrder::third);
    [[nodiscard]] bool advance_joint(const JointGraphGeometry&,double dt,double reservoir_pressure,
        const std::array<double,4>& external_force={},std::span<const double> nodal_flow={},
        std::span<const double> series_pressure={},double outlet_inertance=-1.,
        ContinuousGlottalUpdate source_update=ContinuousGlottalUpdate::relaxation) noexcept;
    // Prescribed stage metrics, constant external drives, one source step.
    // No allocation after prepare. Failure preserves network, source and last.
    // A caller subdividing contact steps must recompute BOTH graph stages and
    // source feedback for each substep. Moving metrics also do geometry work:
    // graph_external_work + graph_port_work alone is not their energy change.
    [[nodiscard]] bool advance(std::span<const double> c1,std::span<const double> c2,
        std::span<const double> g1,std::span<const double> g2,double dt,double reservoir_pressure,
        const std::array<double,4>& external_force={},std::span<const double> nodal_flow={},
        std::span<const double> series_pressure={},
        ContinuousGlottalUpdate update=ContinuousGlottalUpdate::symmetric,double outlet_inertance=-1.,const RadiationGeometry* radiation_geometry=nullptr,
        std::span<const ThermalGeometry> thermal_geometry={},std::span<const ViscousGeometry> viscous_geometry={}) noexcept;
    MovingGraphReference network{};
    ContinuousGlottis source{};
    RadiationBoundary radiation{}; // Explicitly prepare before enabling a radiation step.
    std::vector<SlitThermal> thermal{};
    std::vector<SlitViscous> viscous{};
    GraphContinuousGlottalStep last{};
private:
    friend class JointAccuracy;
    bool joint_ready_{};
    std::size_t joint_radiation_branch_{ground_node};
    std::vector<TreeEdge> joint_edges_;
    RelaxationGraphReference joint_free_,joint_response_,joint_final_;
    std::vector<RelaxationMetric> joint_node_metrics_,joint_branch_metrics_;
    std::vector<double> joint_resistance_,joint_drive_;
    MovingGraphReference free_{},response_{};
    std::vector<double> unit_flow_{};
    std::vector<SlitThermal> thermal_trial_{};
    std::vector<std::size_t> thermal_nodes_{};
    std::vector<SlitViscous> viscous_trial_{};
    std::vector<std::size_t> viscous_branches_{};
    std::size_t source_node_{};
};
}
