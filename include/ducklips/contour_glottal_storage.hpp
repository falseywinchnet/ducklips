#pragma once
#include "ducklips/contour_kinetics.hpp"
#include "ducklips/glottal_storage.hpp"
#include "ducklips/fold_contact.hpp"
#include "ducklips/moving_graph.hpp"
namespace ducklips {
struct ContourGlottalStep final {
    std::array<double,3> flow{},law_residual{},constraint_effort{};
    std::array<double,4> force_residual{},kinetic_effort{};
    std::array<double,2> pressure{},area_effort{},contact_area_effort{};
    double reservoir_work{},load_work{},actuator_work{},external_work{},damping_loss{},passage_loss{};
    double numerical_loss{},energy_change{},energy_residual{};
    std::size_t iterations{};
};
// Joint gas, tissue and passage-momentum solve for a medial contour with
// positive chamber storage. Elastic contact may close either cap.
// Ports carry ambient-reference mass volume. Outlet effort is free_load +
// load_slope * average_outflow, allowing exact elimination of a linear tract.
// Symmetric steps preserve topology. Dissipative and exponentially fitted
// relaxation updates cross cap contact with explicit numerical-loss ledgers.
// Both modes require positive chamber storage and zero closed-cap momentum.
class ContourGlottalStorage final {
public:
    [[nodiscard]] bool prepare(const FoldParameters&,const AirProperties&,double tension,
        const FoldState&,const std::array<double,2>& pressure,double outlet_loss=1.) noexcept;
    [[nodiscard]] bool advance(double dt,double reservoir_pressure,double free_load,double load_slope,
        double rest_area,const std::array<double,4>& external_force={},
        ContourUpdate update=ContourUpdate::symmetric) noexcept;
    [[nodiscard]] bool energy(double&) const noexcept;
    [[nodiscard]] bool configuration(std::array<double,2>& area,ContourKineticState&) const noexcept;
    GlottalStorageState state{};
    std::array<double,3> momentum{};
    ContourGlottalStep last{};
private:
    double mechanical_energy(const std::array<double,4>&,const std::array<double,4>&) const noexcept;
    GasStorage gas_{};
    ContourKinetics kinetics_{};
    FoldContact contact_{};
    std::array<double,4> tissue_mass_{},stiffness_{},damping_{};
    std::array<double,2> coupling_{};
    double width_{};
};
class GraphContourGlottal final {
public:
    void prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
        const FoldParameters&,const AirProperties&,double tension,const FoldState&,
        const std::array<double,2>& pressure,double outlet_loss=1.);
    [[nodiscard]] bool advance(std::span<const double> c1,std::span<const double> c2,
        std::span<const double> g1,std::span<const double> g2,double dt,double reservoir_pressure,double rest_area,
        const std::array<double,4>& external_force={},std::span<const double> nodal_flow={},std::span<const double> series_pressure={},
        ContourUpdate update=ContourUpdate::symmetric) noexcept;
    MovingGraphReference network{};
    ContourGlottalStorage source{};
private:
    MovingGraphReference free_{},response_{};
    std::vector<double> unit_flow_{};
    std::size_t source_node_{};
};
}
