#pragma once
#include "ducklips/tongue.hpp"
#include <span>
#include <vector>

namespace ducklips {
using PlanarContour = std::span<const Point2>;
// Finite-volume geometry reference for an explicitly extruded planar body.
// Coordinates and extrusion depth are metres. Domain and solid contours must
// be simple closed polygons (last vertex is not repeated); solids may overlap.
// Their UNION is excluded, so overlapping body/blade/pad tissue is not counted
// twice. This is not the empirical upper-gap-to-area speech geometry.
// A grid cell can contain disconnected pieces of air. These geometric totals
// must NOT be assigned one common pressure without resolving those pieces.
class PlanarAirspace final {
public:
    void prepare(std::span<const double> x_edges,std::span<const double> y_edges,
                 std::size_t maximum_edges);
    // No allocation after prepare. Geometry is integrated between all polygon
    // vertices, edge crossings and grid crossings; each remaining strip is
    // linear, so its midpoint integral is exact for the supplied polygons.
    // False means invalid input/capacity; output must then be discarded.
    bool build(PlanarContour domain,std::span<const PlanarContour> solids,double depth) noexcept;
    // Row-major cell volumes and free portions of internal grid faces.
    // Outer grid faces are left zero: a caller must prescribe wall/port BCs.
    std::vector<double> volume{},vertical_face{},horizontal_face{};
    [[nodiscard]] std::size_t columns() const noexcept {return x_.size()-1;}
    [[nodiscard]] std::size_t rows() const noexcept {return y_.size()-1;}
private:
    struct Segment {Point2 a{},b{};};
    struct Interval {double low{},high{};};
    void polygon_intervals(PlanarContour polygon,double x,std::vector<Interval>& out) noexcept;
    void free_intervals(PlanarContour domain,std::span<const PlanarContour> solids,double x) noexcept;
    std::vector<double> x_{},y_{},cuts_{},hits_{};
    std::vector<Segment> edges_{};
    std::vector<Interval> domain_{},solid_{},free_{},left_{};
    std::size_t maximum_edges_{};
};
}
