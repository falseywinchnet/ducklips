#pragma once
#include "ducklips/contour_pressure_flow.hpp"
namespace ducklips {
// Dimensionless plane-channel conductance, linearized isothermal BGK gas,
// fully diffuse stationary walls. Valid positive finite delta only.
[[nodiscard]] bool slit_bgk_conductance(double delta,double& result) noexcept;
struct ContourPressureLoss final {
    std::array<double,25> resistance{}; // dissipated power = rate^T R rate
};
// Quasisteady local parallel-slit closure on the continuous contour field.
// This is an explicit constitutive approximation, not full moving-wall BGK.
// Configure identically to the kinetics that prepared the supplied flow state.
class ContourPressureDissipation final {
public:
    [[nodiscard]] bool configure(const FoldParameters&,const AirProperties&) noexcept;
    [[nodiscard]] bool prepare(const ContourPressureFlow&,ContourPressureLoss&) const noexcept;
    [[nodiscard]] bool power(const ContourPressureFlow&,const ContourPressureLoss&,const std::array<double,5>& rate,double&) const noexcept;
private:
    [[nodiscard]] bool integrate(const ContourPressureFlow&,double start,double end,int mapping,double lo,double hi,std::array<double,25>&,const std::array<double,25>* budget=nullptr,int depth=0) const noexcept;
    ContourPressureKinetics kinetics_{};
    std::array<double,2> thickness_{};
    double density_{},thermal_speed_{},gamma_{},viscosity_{},width_{};
};
}
