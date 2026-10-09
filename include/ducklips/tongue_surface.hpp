#pragma once
#include "ducklips/tongue_flow.hpp"

namespace ducklips {
// Owns the material map and a contour used ONLY to bracket exact surface
// intersections. Returned gaps are refined on the mapped circle itself.
class ConservativeTongueSurface final {
public:
    bool configure(Point2 centre,const TongueGesture& gesture) noexcept;
    [[nodiscard]] TongueMaterialPoint material(Point2 point) const noexcept;
    [[nodiscard]] double ray_gap(Point2 wall,Point2 inward,Point2 translation) const noexcept;
    [[nodiscard]] double wall_penetration(Point2 translation) const noexcept;
    [[nodiscard]] double admissible_pull(Point2 translation,Point2 pull) const noexcept;
    [[nodiscard]] const char* error() const noexcept { return error_; }
private:
    static constexpr std::size_t contour_count=128;
    ConservativeTongueMap map_{};
    Point2 centre_{};
    double command_{};
    std::array<Point2,contour_count> contour_{};
    const char* error_{nullptr};
};
}
