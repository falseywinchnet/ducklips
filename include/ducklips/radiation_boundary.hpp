#pragma once
#include "ducklips/physics.hpp"
namespace ducklips {
struct RadiationBoundaryStep final {
    double radiated_energy{},energy_change{},energy_residual{};
};
// Passive one-pole aperture load: Z(s)=M*s/(1+s*M/R), where
// M=rho*end_correction*sqrt(A/pi)/A and R=M*M/radiation_coefficient.
// Its low-frequency impedance is M*s + radiation_coefficient*omega^2.
// This is the existing monopole/end-correction approximation, not a full
// high-frequency radiation solution for an arbitrary lip shape.
class RadiationBoundary final {
public:
    [[nodiscard]] bool prepare(const AirProperties&,double area,double end_correction=.7,double radiation_factor=1.4) noexcept;
    // Preserve canonical momentum through positive-area geometry changes.
    // At an exactly closed aperture retire its absent radiation state.
    [[nodiscard]] bool change_area(double area,double& geometry_work) noexcept;
    // Exact exchange of duct and radiation momentum at fixed geometry.
    // duct_conductance is reciprocal duct inertia; no state allocation.
    [[nodiscard]] bool advance(double& duct_momentum,double duct_conductance,double dt) noexcept;
    [[nodiscard]] double area() const noexcept{return area_;}
    [[nodiscard]] double inverse_inertance() const noexcept{return inverse_;}
    [[nodiscard]] double resistance() const noexcept;
    [[nodiscard]] double energy() const noexcept{return .5*inverse_*momentum*momentum;}
    double momentum{};
    RadiationBoundaryStep last{};
private:
    double area_{},inverse_{},inverse_scale_{},coefficient_{};
};
}
