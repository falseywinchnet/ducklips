#pragma once
#include "ducklips/tongue.hpp"
#include <span>
#include <vector>

namespace ducklips {
using PlanarContour = std::span<const Point2>;
struct PlanarRegion final {
    std::size_t cell{};
    double volume{};
    Point2 interior{}; // representative point in one constituent open strip
};
struct PlanarRegionLink final {
    std::size_t first{},second{},face{};
    bool vertical{};
    double area{};
};
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
    // Optional component decomposition. Preparation allocates bounded scratch.
    // build_regions also builds the ordinary totals; on failure discard both.
    // Region indices describe this pose only, not moving-state identities.
    void prepare_regions(std::size_t piece_capacity,std::size_t link_capacity);
    bool build_regions(PlanarContour domain,std::span<const PlanarContour> solids,double depth) noexcept;
    std::vector<PlanarRegion> regions{};
    std::vector<PlanarRegionLink> region_links{};
    // Row-major cell volumes and free portions of internal grid faces.
    // Outer grid faces are left zero: a caller must prescribe wall/port BCs.
    std::vector<double> volume{},vertical_face{},horizontal_face{};
    [[nodiscard]] std::size_t columns() const noexcept {return x_.size()-1;}
    [[nodiscard]] std::size_t rows() const noexcept {return y_.size()-1;}
private:
    struct Segment {Point2 a{},b{};};
    struct Interval {double low{},high{};std::size_t low_edge{},high_edge{};};
    struct Hit {double y{};std::size_t edge{};};
    struct Piece {std::size_t cell{},parent{},region{},slice{};double low_left{},high_left{},low_right{},high_right{},volume{};Point2 interior{};};
    struct Range {std::size_t begin{},end{};};
    std::size_t piece_root(std::size_t piece) noexcept;
    double boundary_y(std::size_t edge,double x) const noexcept;
    void polygon_intervals(PlanarContour polygon,std::size_t first_edge,double x,std::vector<Interval>& out) noexcept;
    void free_intervals(PlanarContour domain,std::span<const PlanarContour> solids,double x) noexcept;
    std::vector<double> x_{},y_{},cuts_{};
    std::vector<Hit> hits_{};
    std::vector<Segment> edges_{};
    std::vector<Interval> domain_{},solid_{},free_{},left_{};
    std::vector<Piece> pieces_{};
    std::vector<Range> previous_rows_{},current_rows_{};
    std::size_t piece_capacity_{},link_capacity_{};
    std::size_t maximum_edges_{};
};
}
