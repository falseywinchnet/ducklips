#pragma once
#include "ducklips/glottal_contour.hpp"
#include "ducklips/gas_storage.hpp"
#include "ducklips/sources.hpp"
namespace ducklips {
// Alternative explicit axial shape assumptions; neither is measured anatomy.
enum class ContourPressureProfile { flat_caps, linear_span };
struct ContourPressureState final {
    // Coefficients of continuous specific-enthalpy effort, in Pa.
    std::array<double,2> effort{},mass{},area_effort{};
    // Symmetric matrices, row-major. mass_area = d mass / d control area.
    std::array<double,4> compliance{},mass_area{};
    double volume{},energy{},coenergy{};
};
// Continuous gas pressure over the positive part of the medial contour.
// The two shape functions are the same linear hat functions as the contour.
// mass[i] = integral A*r*phi[i] dz; their sum is total reference-mass volume.
// These are overlapping moments, NOT the old disjoint chamber masses.
// This storage component does not yet supply momentum or a speech source.
class ContourPressureStorage final {
public:
    // Supremum of the affine effort-field difference over the UNION of both
    // occupied supports. Includes contact limits, not extrapolated values in
    // tissue. The hat-coordinate interval is the same for both axial profiles.
    // Returns infinity for invalid or empty supports. Units: Pa of effort.
    [[nodiscard]] static double effort_difference(const std::array<double,2>& first_area,
        const std::array<double,2>& first_effort,const std::array<double,2>& last_area,
        const std::array<double,2>& last_effort) noexcept;
    [[nodiscard]] bool configure(const FoldParameters&,const AirProperties&,ContourPressureProfile=ContourPressureProfile::flat_caps) noexcept;
    [[nodiscard]] bool evaluate(const std::array<double,2>& area,
        const std::array<double,2>& effort,ContourPressureState&) const noexcept;
    // Restricted integrals for continuity reconstruction. Empty subintervals
    // return zero moments; full storage still requires positive air volume.
    [[nodiscard]] bool interval(const std::array<double,2>& area,
        const std::array<double,2>& effort,double first,double last,ContourPressureState&) const noexcept;
    [[nodiscard]] bool recover(const std::array<double,2>& area,
        const std::array<double,2>& mass,const std::array<double,2>& guess,
        ContourPressureState&) const noexcept;
private:
    GasStorage gas_{};
    std::array<double,2> thickness_{},basis_{};
    double bulk_{},gamma_{};
};
}
