#pragma once
#include "ducklips/gas_storage.hpp"
#include "ducklips/sources.hpp"

namespace ducklips {
struct GlottalStorageState final {
    FoldState folds{};
    std::array<double,2> mass{}; // m3 of air at reference density
};
struct GlottalStorageStep final {
    std::array<double,2> average_pressure{},average_mass_effort{};
    std::array<double,3> flow_effort{}; // conjugate to lower-in, intercell, upper-out flow
    double port_work{},wall_work{},actuator_work{},external_work{},damping_loss{},energy_change{},energy_residual{};
    std::size_t iterations{};
};
// Two positive-volume rectangular glottal chambers attached to the four
// existing tissue masses (one strip). Gas pressure comes from stored mass
// and actual volume. Its force and displacement do matching work. Port
// flows are held inputs, in ambient-reference m3/s; inlet and outlet need
// not agree. This component has no jet/inertance/closure transport model.
class GlottalStorage final {
public:
    [[nodiscard]] bool prepare(const FoldParameters& folds,const AirProperties& air,double tension,
        const FoldState& initial,const std::array<double,2>& pressure) noexcept;
    // Flow signs: into lower cell, lower to upper cell, out of upper cell.
    // Force ordering: lower-left, lower-right, upper-left, upper-right.
    // Positive masses/volumes only. Failure preserves state and receipt.
    [[nodiscard]] bool advance(const std::array<double,3>& flow,double rest_area,
        const std::array<double,4>& external_force,double dt) noexcept;
    [[nodiscard]] bool energy(double& result) const noexcept;
    [[nodiscard]] bool volumes(std::array<double,2>& result) const noexcept;
    GlottalStorageState state{};
    GlottalStorageStep last{};
private:
    double mechanical_energy(const std::array<double,4>& x,const std::array<double,4>& velocity) const noexcept;
    GasStorage gas_{};
    std::array<double,4> tissue_mass_{},stiffness_{},damping_{};
    std::array<double,2> coupling_{},thickness_{},wall_area_{};
    double length_{};
};
}
