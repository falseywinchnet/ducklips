#pragma once
#include "ducklips/physics.hpp"
#include <array>
#include <cstddef>
namespace ducklips {
struct SlitThermalStep final {
    double dissipated_energy{},energy_change{},energy_residual{},stored_air_residual{};
};
// Linear heat diffusion to isothermal parallel walls. The graph's primary
// charge is C_ad P; charge[] holds the additional reference-volume mass
// from thermal relaxation. Total stored air includes BOTH. No through-flow.
// Geometry remains positive; disappearing storage requires explicit transport.
class SlitThermal final {
public:
    static constexpr std::size_t modes=24;
    [[nodiscard]] bool prepare(const AirProperties&,double compliance,double half_gap) noexcept;
    [[nodiscard]] bool change_geometry(double compliance,double half_gap,double& work) noexcept;
    [[nodiscard]] bool equilibrate(double pressure) noexcept;
    [[nodiscard]] bool advance(double& primary_charge,double dt) noexcept;
    [[nodiscard]] double stored_air() const noexcept;
    [[nodiscard]] double energy() const noexcept;
    [[nodiscard]] double compliance() const noexcept{return compliance_;}
    [[nodiscard]] double half_gap() const noexcept{return half_gap_;}
    [[nodiscard]] double memory_storage(std::size_t mode) const noexcept;
    [[nodiscard]] double memory_rate(std::size_t mode) const noexcept;
    std::array<double,modes> charge{};
    SlitThermalStep last{};
private:
    std::array<double,modes> capacity_{},scale_{},rate_{};
    std::array<double,modes*modes> vectors_{};
    double compliance_{},half_gap_{},diffusivity_{},capacity_sum_{};
};
}
