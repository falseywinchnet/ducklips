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
// Open part of a selected DOMAIN edge, attached to one interior air region.
// This is a geometric attachment, not a glottal/lip flow or radiation law.
struct PlanarBoundaryPort final {
    std::size_t edge{},region{};
    double area{};
    // Area centroid of all openings for this edge/region pair. It may lie in
    // a covered gap between openings; it is not an interior witness point.
    Point2 centre{},outward_normal{};
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
    // Selected domain-edge indices are unique. Portions covered by solid
    // tissue (including coincident contact from outside) are sealed. Each
    // remaining positive-area portion retains its own adjacent region.
    bool build_regions(PlanarContour domain,std::span<const PlanarContour> solids,double depth,
                       std::span<const std::size_t> boundary_edges = {}) noexcept;
    // Linear motion of corresponding contour vertices, followed by fresh
    // region/volume/aperture extraction. This does NOT blend those outputs.
    // The caller supplies matching vertex identities and simple intermediate
    // polygons. It does not establish a muscle law or moving-state identity.
    bool build_between(PlanarContour first_domain,std::span<const PlanarContour> first_solids,
        PlanarContour last_domain,std::span<const PlanarContour> last_solids,double fraction,double depth,
        std::span<const std::size_t> boundary_edges = {}) noexcept;
    // Actual sampled contours for diagnostics: domain first, then solids.
    // Valid after successful build_between, until next preparation/motion call.
    [[nodiscard]] std::span<const PlanarContour> interpolated_contours() const noexcept {return motion_contours_;}
    std::vector<PlanarRegion> regions{};
    std::vector<PlanarRegionLink> region_links{};
    std::vector<PlanarBoundaryPort> boundary_ports{};
    // Row-major cell volumes and free portions of internal grid faces.
    // Outer grid faces are left zero: a caller must prescribe wall/port BCs.
    std::vector<double> volume{},vertical_face{},horizontal_face{};
    [[nodiscard]] std::size_t columns() const noexcept {return x_.size()-1;}
    [[nodiscard]] std::size_t rows() const noexcept {return y_.size()-1;}
private:
    struct Segment {Point2 a{},b{};};
    struct Interval {double low{},high{};std::size_t low_edge{},high_edge{};};
    struct Hit {double y{};std::size_t edge{};};
    struct Piece {std::size_t cell{},parent{},region{},slice{};double low_left{},high_left{},low_right{},high_right{},volume{};Point2 interior{};
        double x_left{},x_right{};std::size_t low_edge{},high_edge{};};
    struct Range {std::size_t begin{},end{};};
    std::size_t piece_root(std::size_t piece) noexcept;
    double boundary_y(std::size_t edge,double x) const noexcept;
    void polygon_intervals(PlanarContour polygon,std::size_t first_edge,double x,std::vector<Interval>& out) noexcept;
    void free_intervals(PlanarContour domain,std::span<const PlanarContour> solids,double x) noexcept;
    bool attach_boundaries(PlanarContour domain,std::span<const PlanarContour> solids,double depth,
                           std::span<const std::size_t> boundary_edges) noexcept;
    std::vector<double> x_{},y_{},cuts_{};
    std::vector<Hit> hits_{};
    std::vector<Segment> edges_{};
    std::vector<Interval> domain_{},solid_{},free_{},left_{};
    std::vector<Interval> boundary_blocked_{};
    std::vector<Point2> motion_points_{};
    std::vector<PlanarContour> motion_contours_{};
    std::vector<Piece> pieces_{};
    std::vector<Range> previous_rows_{},current_rows_{};
    std::size_t piece_capacity_{},link_capacity_{};
    std::size_t maximum_edges_{};
};
}
