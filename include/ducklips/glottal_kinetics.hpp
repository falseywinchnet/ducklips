#pragma once
#include "ducklips/sources.hpp"

namespace ducklips {
struct GlottalKineticState final {
    std::array<double,3> inverse_inertance{};
    std::array<double,2> mass_effort{};
    double energy{};
};
struct GlottalKineticTransfer final {
    std::array<double,3> flow{};
    std::array<double,2> mass_effort{};
    double energy_change{},chain_residual{};
};
// Three serial dual passages: lower half, lower/upper halves, upper half.
// For fixed longitudinal thickness h, half-cell inertance is rho0*h^2/(2M).
// M is ambient-reference mass volume, not instantaneous geometric volume.
// K=lambda^2/(2L(M)); its mass derivative contributes kinetic enthalpy.
class GlottalKinetics final {
public:
    [[nodiscard]] bool configure(const FoldParameters& folds,const AirProperties& air,double outlet_loss) noexcept;
    [[nodiscard]] bool evaluate(const std::array<double,2>& mass,const std::array<double,3>& momentum,GlottalKineticState& result) const noexcept;
    [[nodiscard]] bool transfer(const std::array<double,2>& mass0,const std::array<double,2>& mass1,
        const std::array<double,3>& momentum0,const std::array<double,3>& momentum1,GlottalKineticTransfer& result) const noexcept;
    // Parallel-plate laminar resistance, converted to reference-mass flow.
    // Optional symmetric terminal kinetic-head loss, not a separation model.
    [[nodiscard]] bool losses(const std::array<double,2>& mass,const std::array<double,2>& volume,
        std::array<double,3>& resistance,std::array<double,3>& quadratic) const noexcept;
private:
    double density_{},viscosity_{},width_{},outlet_loss_{};
    std::array<double,2> thickness_{};
};
}
