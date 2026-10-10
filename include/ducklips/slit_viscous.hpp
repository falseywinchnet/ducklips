#pragma once
#include "ducklips/physics.hpp"
#include <array>
#include <cstddef>
namespace ducklips {
struct SlitViscousStep final {double dissipated_energy{},energy_change{},energy_residual{},blocked_impulse{};};
// No-slip, stationary parallel-wall Stokes impedance of one axial segment.
// The graph owns the complete branch's plug-flow inertia (1/conductance).
// This object adds that segment's wall drag and transverse profile memory,
// not another plug-flow inertia. Several segments may share one branch.
// A zero branch conductance constrains through-flow while any surviving local
// profile relaxes. The constraint absorbs blocked_impulse with zero flow work.
// The caller owns the branch's canonical momentum: at zero conductance it is
// an unobservable gauge, left unchanged here and reset by the graph's closure
// handling before reopening. This object owns only local profile memory.
// The exact tuple (0,0,0) denotes an absent segment, not a finite inertance.
// Changing to absence retires its memory and returns its metric-energy change;
// reopening starts that segment with no resurrected velocity profile.
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
