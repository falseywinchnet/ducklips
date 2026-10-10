#pragma once
#include "ducklips/physics.hpp"
#include <array>
#include <cstddef>
namespace ducklips {
struct SlitViscousStep final {double dissipated_energy{},energy_change{},energy_residual{};};
// No-slip, stationary parallel-wall Stokes impedance of one axial segment.
// The graph owns the complete branch's plug-flow inertia (1/conductance).
// This object adds that segment's wall drag and transverse profile memory,
// not another plug-flow inertia. Several segments may share one branch.
class SlitViscous final {
public:
    static constexpr std::size_t modes=24;
    [[nodiscard]] bool prepare(const AirProperties&,double segment_inertance,double branch_conductance,double half_gap) noexcept;
    [[nodiscard]] bool change_geometry(double segment_inertance,double branch_conductance,double half_gap,double& work) noexcept;
    [[nodiscard]] bool equilibrate(double flow) noexcept;
    [[nodiscard]] bool advance(double& branch_momentum,double dt) noexcept;
    [[nodiscard]] double energy() const noexcept;
    [[nodiscard]] double segment_inertance() const noexcept{return inertance_;}
    [[nodiscard]] double branch_conductance() const noexcept{return conductance_;}
    [[nodiscard]] double half_gap() const noexcept{return half_gap_;}
    std::array<double,modes> momentum{};
    SlitViscousStep last{};
private:
    static constexpr std::size_t dimension=modes+1;
    [[nodiscard]] bool diagonalize(double ratio) noexcept;
    std::array<double,dimension> rate_{};
    std::array<double,dimension*dimension> vectors_{};
    double inertance_{},conductance_{},half_gap_{},viscosity_{};
};
}
