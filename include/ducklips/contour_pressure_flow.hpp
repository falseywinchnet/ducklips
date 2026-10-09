#pragma once
#include "ducklips/contour_pressure_storage.hpp"
namespace ducklips {
struct ContourPressureFlow final {
    std::array<double,2> area{},effort{},support{};
    // Columns: Qin, integrated Qmiddle, Qout, lower area rate, upper area rate.
    std::array<double,10> effort_rate{};
    std::array<double,25> inertia{}; // K = rate^T inertia rate / 2
    // For a moving contact endpoint, the no-slip slit loss diverges as
    // (area_rate^T contact_loss area_rate) log(1/axial_cutoff).
    // This is a divergence coefficient, not a finite replacement loss.
    std::array<double,4> contact_loss{};
};
struct ContourFlowSample final {
    double area{},relative_density{};
    std::array<double,5> flow{}; // Q(z) = flow dot rate
};
// Strong local continuity for the continuous contour gas field. Includes
// axial fluid inertia induced by moving walls and density redistribution.
// This component does not advance generalized momenta or regularize contact.
class ContourPressureKinetics final {
public:
    [[nodiscard]] bool configure(const FoldParameters&,const AirProperties&) noexcept;
    [[nodiscard]] bool prepare(const std::array<double,2>& area,const std::array<double,2>& effort,ContourPressureFlow&) const noexcept;
    [[nodiscard]] bool sample(const ContourPressureFlow&,double z,ContourFlowSample&) const noexcept;
    [[nodiscard]] bool energy(const ContourPressureFlow&,const std::array<double,5>& rate,double&) const noexcept;
private:
    [[nodiscard]] bool integrate(const ContourPressureFlow&,double lo,double hi,std::array<double,25>&,int depth=0) const noexcept;
    ContourPressureStorage storage_{};
    std::array<double,2> thickness_{};
    double density_{},bulk_{},gamma_{},viscosity_{},width_{};
};
}
