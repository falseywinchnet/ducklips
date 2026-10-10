#pragma once
#include "ducklips/tongue_flow.hpp"

namespace ducklips {
// Owns the material map and a contour used ONLY to bracket exact surface
// intersections. Returned gaps are refined on the mapped circle itself.
class ConservativeTongueSurface final {
public:
    bool configure(Point2 centre,const TongueGesture& gesture) noexcept;
    // Transport the blade attachment by deforming nearby body material.
    // Its moving support must fit between the oral wall and lower anchor.
    bool move_attachment(Point2 start,Point2 displacement) noexcept;
    [[nodiscard]] TongueMaterialPoint material(Point2 point) const noexcept;
    [[nodiscard]] double ray_gap(Point2 wall,Point2 inward) const noexcept;
    // Engine representation: intersect the carried polygon directly. The
    // mapped-circle refinement remains the offline geometric reference.
    void set_polygon_queries(bool enabled) noexcept { polygon_queries_=enabled; }
    [[nodiscard]] double wall_penetration() const noexcept;
    [[nodiscard]] const char* error() const noexcept { return error_; }
private:
    bool update_contour() noexcept;
    static constexpr std::size_t contour_count=128;
    ConservativeTongueMap map_{};
    TongueMaterialTranslation attachment_map_{};
    bool attachment_active_{false};
    bool polygon_queries_{false};
    Point2 centre_{};
    double command_{};
    std::array<Point2,contour_count> contour_{};
    const char* error_{nullptr};
};
}
