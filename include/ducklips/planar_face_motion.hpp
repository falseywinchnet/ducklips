#pragma once
#include "ducklips/planar_airspace.hpp"

namespace ducklips {
// A fixed grid face in an extruded planar airway. The normal coordinate is x
// for a vertical face, y otherwise; low/high are tangential coordinates.
struct PlanarFace final {bool vertical{};double coordinate{},low{},high{};};
struct PlanarFaceWindow final {
    double first{},last{}; // normalized time, not seconds
    double average_area{},error_estimate{}; // m^2, quadrature estimate, not a rigorous bound
};
// Partition a single face's motion at changes in its contour intersections.
// This is a geometric opening history, NOT region identity or a flow law.
// Matching vertex identities and simple intermediate polygons are required.
class PlanarFaceMotion final {
public:
    void prepare(std::size_t maximum_edges,std::size_t maximum_events);
    // No allocations after preparation. False invalidates all outputs.
    // area_tolerance is an absolute error target for the normalized integral.
    bool build(PlanarContour first_domain,std::span<const PlanarContour> first_solids,
        PlanarContour last_domain,std::span<const PlanarContour> last_solids,
        PlanarFace face,double depth,double area_tolerance=1e-13) noexcept;
    std::vector<PlanarFaceWindow> windows{};
    double average_area{},error_estimate{};
private:
    struct Edge {Point2 a{},b{},da{},db{};std::size_t contour{};};
    struct Interval {double low{},high{};};
    struct Crossing {long double n[3]{},d[2]{};};
    bool add_event(long double time) noexcept;
    bool roots(const long double* polynomial,int degree) noexcept;
    double area(double fraction) noexcept;
    void intervals(std::size_t contour,double fraction,bool positive,std::vector<Interval>& output) noexcept;
    bool integrate(double first,double last,double tolerance,unsigned level,double& integral,double& error) noexcept;
    std::vector<Edge> edges_{};
    std::vector<Crossing> crossings_{};
    std::vector<double> events_{},hits_{};
    std::vector<Interval> domain_left_{},domain_right_{},solids_{};
    std::size_t maximum_edges_{},maximum_events_{},contours_{};
    PlanarFace face_{};double depth_{};
};
}
