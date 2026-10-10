#pragma once
#include "ducklips/lip_profile.hpp"
#include <span>
namespace ducklips {
struct LipViscousSegment final {double position{},inertance{},half_gap{};};
struct LipViscousMeasure final {
    double conductance{},quadrature_inertance{},relative_inertance_correction{};
};
// Positive quadrature of local slit impedances along an open axial lip path.
// Normalize the dx/A quadrature measure to the profile's independently
// integrated inertia, preserving the graph's existing plug-flow term exactly.
// Orders 4, 8, 16, 32. Output identities follow fixed reference abscissae.
// Closed paths and invalid geometry fail without changing accepted outputs.
// Spatial accuracy must be checked under order refinement; this is not a
// uniform-gap approximation or a complete moving-wall Stokes PDE.
[[nodiscard]] bool lip_viscous_segments(const LipProfile&,double first,double last,double slit_width,double density,
    std::span<LipViscousSegment>,LipViscousMeasure&) noexcept;
}
