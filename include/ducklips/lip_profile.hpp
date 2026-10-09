#pragma once
namespace ducklips {
// The existing axial lip-slot constriction, expressed as an area field rather
// than one midpoint sample per tube section. SI units; x runs inward from lips.
// No new anatomy or minimum aperture is introduced by this representation.
struct LipProfile final {
    double length{},base_area{},activation{},place{},core{},width{},target_area{};
    [[nodiscard]] bool valid() const noexcept;
    [[nodiscard]] double minimum_area() const noexcept;
    [[nodiscard]] double area(double x) const noexcept;
    // Integrates positive air, including a partly occupied section. Empty
    // sections return volume=0 and the interval midpoint as unused centroid.
    [[nodiscard]] bool measure(double first,double last,double& volume,double& centroid) const noexcept;
    // 1 / [rho integral(dx/A)]. Any closed point on the path seals it exactly.
    // A zero conductance does not imply that adjoining storage is empty.
    [[nodiscard]] bool conductance(double first,double last,double density,double& result) const noexcept;
};
}
