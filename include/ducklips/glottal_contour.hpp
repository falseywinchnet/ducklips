#pragma once
#include <array>
namespace ducklips {
struct GlottalContourMeasure final {
    double volume{},centroid{};
    std::array<double,2> area_derivative{}; // dV/d signed lower/upper control area
};
struct GlottalContourPassage final {
    bool open{};
    // Integrals split at the chamber boundary. Combine with each chamber's
    // own density when constructing inertia/resistance; no ambient assumption.
    std::array<double,2> inverse_area{},inverse_area_cubed{};
};
struct GlottalContourGeometry final {
    std::array<GlottalContourMeasure,2> chamber{};
    std::array<double,3> aperture{}; // inlet, chamber seam, outlet
    std::array<GlottalContourPassage,3> passage{}; // caps to air centroids, and between centroids
};
// Explicit geometry reference: signed medial opening is constant from each
// cap to its tissue centre, and affine between the two centres. Air is the
// positive part. This is a stated profile assumption, not measured anatomy.
// Sealing a passage does not erase adjoining wedge-shaped air storage.
struct GlottalContour final {
    double lower_thickness{},upper_thickness{},lower_area{},upper_area{};
    [[nodiscard]] bool valid() const noexcept;
    [[nodiscard]] double signed_area(double z) const noexcept;
    [[nodiscard]] bool measure(double first,double last,GlottalContourMeasure& result) const noexcept;
    [[nodiscard]] bool passage(double first,double last,GlottalContourPassage& result) const noexcept;
    [[nodiscard]] bool geometry(GlottalContourGeometry& result) const noexcept;
};
}
