#pragma once
#include "ducklips/physics.hpp"

namespace ducklips {
struct GasEnergy final {
    double energy{};       // J, available compression energy relative to ambient
    double pressure{};     // Pa, relative thermodynamic pressure (wall effort)
    double mass_effort{};  // Pa, enthalpy per ambient-reference volume
};
struct GasTransfer final {
    double mass_effort{},pressure{};
    double pressure_volume_slope{}; // d discrete pressure / d final volume, fixed masses
};
// Isentropic ideal-gas storage, with the configured small-signal bulk modulus
// rho0*c0^2. Mass is expressed as equivalent volume at ambient density.
// The mass-flow and wall-volume efforts are distinct derivatives of one
// energy. Positive finite mass and volume only; no disappearing-cell floor.
class GasStorage final {
public:
    [[nodiscard]] bool configure(const AirProperties& air) noexcept;
    [[nodiscard]] bool evaluate(double mass,double volume,GasEnergy& result) const noexcept;
    // Unit-volume constitutive state directly from enthalpy effort. Avoids
    // rounding a tiny density departure through 1 before computing energy.
    [[nodiscard]] bool evaluate_effort(double effort,GasEnergy&,double& relative_density) const noexcept;
    [[nodiscard]] bool transfer(double mass0,double volume0,double mass1,double volume1,GasTransfer& result) const noexcept;
    [[nodiscard]] bool mass_at_pressure(double volume,double pressure,double& mass) const noexcept;
private:
    double density_energy(double relative_density) const noexcept;
    double density_departure_energy(double departure) const noexcept;
    double bulk_{},gamma_{};
};
}
