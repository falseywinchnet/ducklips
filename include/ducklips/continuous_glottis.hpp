#pragma once
#include "ducklips/contour_pressure_loss.hpp"
#include "ducklips/fold_contact.hpp"
namespace ducklips {
struct ContinuousGlottalState final {
    FoldState folds{};
    std::array<double,2> mass{},effort{},area{};
    std::array<double,3> flow{};
    // Canonical momenta conjugate to inlet/middle/outlet transferred volume
    // and the four tissue positions. Closed-port coordinates are eliminated;
    // their reported momenta and velocities are zero.
    std::array<double,7> momentum{};
    double energy{};
};
struct ContinuousGlottalStep final {
    std::array<double,3> flow{};
    double reservoir_work{},load_work{},external_work{},fluid_loss{},tissue_loss{};
    double energy_change{},energy_residual{},maximum_scaled_residual{};
    std::size_t iterations{};
};
// Continuous-pressure, shared fluid/tissue Hamiltonian reference. The current
// advance holds rest area and port topology fixed; a crossing is rejected
// transactionally. No disappearing chamber, independent passage inertia,
// source waveform, or prescribed tissue oscillation is used.
class ContinuousGlottis final {
public:
    [[nodiscard]] bool prepare(const FoldParameters&,const AirProperties&,double tension,
        const FoldState&,const std::array<double,2>& effort,
        const std::array<double,3>& flow={},double outlet_loss=1.) noexcept;
    [[nodiscard]] bool advance(double dt,double reservoir_pressure,double free_load,double load_slope,
        const std::array<double,4>& external_force={}) noexcept;
    [[nodiscard]] bool state(ContinuousGlottalState&) const noexcept;
    ContinuousGlottalStep last{};
private:
    struct Configuration final {
        std::array<double,4> z{},potential_gradient{};
        std::array<double,49> metric{},inverse{};
        ContourPressureState gas{};
        ContourPressureFlow flow{};
        double potential{};
    };
    [[nodiscard]] std::array<double,4> coordinates(const std::array<double,7>&) const noexcept;
    [[nodiscard]] bool configuration(const std::array<double,4>&,Configuration&) const noexcept;
    [[nodiscard]] double mechanical(const std::array<double,7>&) const noexcept;
    [[nodiscard]] bool resistance(const Configuration&,std::array<double,49>&) const noexcept;
    ContourPressureStorage storage_{};
    ContourPressureKinetics kinetics_{};
    ContourPressureDissipation loss_{};
    FoldContact contact_{};
    GasStorage gas_{};
    std::array<double,7> position_{},momentum_{};
    std::array<double,2> initial_mass_{},effort_guess_{};
    std::array<double,4> mass_{},stiffness_{},damping_{};
    std::array<double,2> coupling_{};
    std::array<bool,7> active_{};
    double rest_area_{},width_{},density_{},bulk_{},gamma_{},outlet_loss_{};
};
}
