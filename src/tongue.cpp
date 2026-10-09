#include "ducklips/tongue.hpp"
#include "ducklips/physics.hpp"
#include "ducklips/tongue_surface.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace ducklips {
namespace {

constexpr double no_hit = std::numeric_limits<double>::infinity();
constexpr double bend_start = MidsagittalTongue::pharynx_height;
constexpr double bend_end = MidsagittalTongue::pharynx_height + 0.5 * pi * MidsagittalTongue::bend_radius;
constexpr double blade_angle = 25.0 * pi / 180.0; // blade root on the body, above horizontal

struct Knot final {
    double x{};
    double value{};
};

// Outer wall distance above the oral midline against x (cm): soft palate,
// palatal dome, alveolar ridge, upper incisors, lip region.
constexpr std::array<Knot, 6> palate_knots{{
    {2.0, 0.9}, {4.5, 1.1}, {6.0, 0.9}, {7.2, 0.55}, {8.4, 0.5}, {9.4, 0.6},
}};

double interpolate(const std::array<Knot, 6>& knots, double x) noexcept {
    if (x <= knots[0].x) {
        return knots[0].value;
    }
    for (std::size_t index = 1; index < knots.size(); ++index) {
        if (x <= knots[index].x) {
            double t = (x - knots[index - 1].x) / (knots[index].x - knots[index - 1].x);
            double value = knots[index - 1].value + t * (knots[index].value - knots[index - 1].value);
            return value;
        }
    }
    return knots[knots.size() - 1].value;
}

Point2 add(const Point2& a, const Point2& b) noexcept {
    Point2 sum{.x = a.x + b.x, .y = a.y + b.y};
    return sum;
}

Point2 scaled(const Point2& a, double factor) noexcept {
    Point2 product{.x = a.x * factor, .y = a.y * factor};
    return product;
}

double dot(const Point2& a, const Point2& b) noexcept {
    double value = a.x * b.x + a.y * b.y;
    return value;
}

double norm(const Point2& a) noexcept {
    double value = std::sqrt(dot(a, a));
    return value;
}

// Distance along a unit ray from origin to a circle; 0 when the origin is
// inside it, infinity when the ray misses.
double ray_circle(const Point2& origin, const Point2& direction, const Point2& centre, double radius) noexcept {
    Point2 offset{.x = origin.x - centre.x, .y = origin.y - centre.y};
    double c = dot(offset, offset) - radius * radius;
    if (c <= 0.0) {
        return 0.0;
    }
    double b = dot(offset, direction);
    double discriminant = b * b - c;
    if (discriminant < 0.0) {
        return no_hit;
    }
    double t = -b - std::sqrt(discriminant);
    if (t < 0.0) {
        return no_hit;
    }
    return t;
}

// Distance along a unit ray to a segment; 0 when the segment lies on the
// wall's outer side of the ray origin (pressed into the wall).
double ray_segment(const Point2& origin, const Point2& direction, const Point2& start, const Point2& end) noexcept {
    Point2 edge{.x = end.x - start.x, .y = end.y - start.y};
    double denominator = direction.x * edge.y - direction.y * edge.x;
    if (std::fabs(denominator) < 1.0e-12) {
        return no_hit;
    }
    Point2 offset{.x = start.x - origin.x, .y = start.y - origin.y};
    double t = (offset.x * edge.y - offset.y * edge.x) / denominator;
    double u = (offset.x * direction.y - offset.y * direction.x) / denominator;
    if (u < 0.0 || u > 1.0) {
        return no_hit;
    }
    double distance = std::max(t, 0.0);
    return distance;
}

// Smooth minimum of two distances (cm): soft tissue surfaces meet without
// corners, which keeps the area function differentiable for the posture fit.
double smooth_minimum(double first, double second) noexcept {
    constexpr double softness = 0.04; // cm
    if (first == no_hit) {
        return second;
    }
    if (second == no_hit) {
        return first;
    }
    double low = std::min(first, second);
    double high = std::max(first, second);
    double blended = low - softness * std::log1p(std::exp(-(high - low) / softness));
    return blended;
}

Point2 blade_root(const Point2& centre) noexcept {
    Point2 root{.x = centre.x + MidsagittalTongue::body_radius * std::cos(blade_angle),
                .y = centre.y + MidsagittalTongue::body_radius * std::sin(blade_angle)};
    return root;
}

}

Point2 MidsagittalTongue::midline(double arc) noexcept {
    if (arc <= bend_start) {
        Point2 point{.x = 0.0, .y = arc};
        return point;
    }
    if (arc <= bend_end) {
        double angle = pi - (arc - bend_start) / bend_radius;
        Point2 point{.x = bend_radius + bend_radius * std::cos(angle), .y = bend_start + bend_radius * std::sin(angle)};
        return point;
    }
    Point2 point{.x = bend_radius + (arc - bend_end), .y = bend_start + bend_radius};
    return point;
}

Point2 MidsagittalTongue::normal(double arc) noexcept {
    if (arc <= bend_start) {
        Point2 direction{.x = -1.0, .y = 0.0};
        return direction;
    }
    if (arc <= bend_end) {
        double angle = pi - (arc - bend_start) / bend_radius;
        Point2 direction{.x = std::cos(angle), .y = std::sin(angle)};
        return direction;
    }
    Point2 direction{.x = 0.0, .y = 1.0};
    return direction;
}

double MidsagittalTongue::outer_offset(double arc) noexcept {
    if (arc <= bend_end) {
        return 0.9;
    }
    double offset = interpolate(palate_knots, midline(arc).x);
    return offset;
}

double MidsagittalTongue::floor_offset(double arc, const TonguePosture& posture) noexcept {
    // Deepest inner boundary when the tongue body is out of the way: the
    // tongue root and vallecula below the bend (advanced root widens the
    // pharynx), the floor of the mouth, and the lower incisors at the front;
    // the jaw opens the oral part.
    double jaw = posture.jaw;
    if (arc <= bend_end) {
        // The epiglottis and lower tongue root bound the hypopharynx closely
        // (the piriform recesses beside it are a separate side branch); the
        // boundary recedes toward the oropharynx.
        double lower = 0.55;
        double upper = 1.4;
        double blend = std::clamp((arc - 4.0) / 2.5, 0.0, 1.0);
        double depth = lower + blend * (upper - lower);
        double root = depth * (1.0 + 0.45 * std::clamp(posture.root, -1.0, 1.0));
        return root;
    }
    double x = midline(arc).x;
    if (x >= 8.0) {
        double teeth = 0.35 + 1.0 * std::clamp(jaw, 0.0, 1.0);
        return teeth;
    }
    double floor = 1.0 + 1.4 * std::clamp(jaw, 0.0, 1.0);
    return floor;
}

double MidsagittalTongue::area_from_gap(double gap, double arc) noexcept {
    double width = std::max(gap, 0.0);
    // Pharyngeal and oral cross-sections relate to the midsagittal gap by
    // different power laws; blended across the bend.
    double pharyngeal = 1.6 * std::pow(width, 1.35);
    double oral = 1.9 * std::pow(width, 1.45);
    double blend = std::clamp((arc - bend_start) / (bend_end - bend_start), 0.0, 1.0);
    double area = pharyngeal + blend * (oral - pharyngeal);
    return area;
}

double MidsagittalTongue::gap_from_area(double area, double arc) noexcept {
    double blend = std::clamp((arc - bend_start) / (bend_end - bend_start), 0.0, 1.0);
    double factor = 1.6 + blend * (1.9 - 1.6);
    double exponent = 1.35 + blend * (1.45 - 1.35);
    double gap = std::pow(std::max(area, 0.0) / factor, 1.0 / exponent);
    return gap;
}

TongueShape MidsagittalTongue::pose(const TonguePosture& posture, const TongueGesture& body,
                                    const TongueGesture& tip, const ConservativeTongueSurface* surface) noexcept {
    TongueShape shape{};
    double jaw = std::clamp(posture.jaw, 0.0, 1.0);
    // The jaw carries the tongue down with it.
    Point2 centre{.x = posture.body_x + posture.offset_x, .y = posture.body_y - 0.9 * jaw + posture.offset_y};

    // Body gesture: the dorsum moves along the wall normal at its place until
    // the gap reaches the target; activation above 1 presses further.
    if (body.activation > 0.0 && surface==nullptr) {
        double arc = length - body.place_from_lips;
        Point2 wall = add(midline(arc), scaled(normal(arc), outer_offset(arc)));
        Point2 inward = scaled(normal(arc), -1.0);
        double gap_now = ray_circle(wall, inward, centre, body_radius);
        if (gap_now == no_hit) {
            gap_now = norm(Point2{.x = wall.x - centre.x, .y = wall.y - centre.y}) - body_radius;
        }
        double reach = std::clamp(body.activation, 0.0, 1.0);
        double amount = reach * reach * (3.0 - 2.0 * reach);
        double press = 0.4 * std::max(body.activation - 1.0, 0.0);
        double shift = amount * (gap_now - body.gap) + press;
        centre = add(centre, scaled(normal(arc), shift));
    }

    Point2 root = blade_root(centre);
    if (surface!=nullptr) {
        const TongueMaterialPoint mapped_root=surface->material(root),mapped_centre=surface->material(centre);
        if (!mapped_root.valid || !mapped_centre.valid) { shape.valid=false;return shape; }
        root=mapped_root.position;centre=mapped_centre.position;shape.surface=surface;
    }
    double rest_angle = (-20.0 + 50.0 * std::clamp(posture.blade, -1.0, 1.0)) * pi / 180.0;
    Point2 rest_tip{.x = root.x + blade_length * std::cos(rest_angle), .y = root.y + blade_length * std::sin(rest_angle)};
    Point2 tip_point = rest_tip;
    if (tip.activation > 0.0) {
        double arc = length - tip.place_from_lips;
        Point2 wall_normal = normal(arc);
        Point2 wall = add(midline(arc), scaled(wall_normal, outer_offset(arc)));
        // Both rest_tip and target are pad centres. Offset the target
        // surface by its radius BEFORE interpolation; subtracting a whole
        // radius afterwards teleports the pad as activation leaves zero.
        Point2 target = add(wall, scaled(wall_normal, -(tip.gap + tip_radius)));
        double reach = std::clamp(tip.activation, 0.0, 1.0);
        double amount = reach * reach * (3.0 - 2.0 * reach);
        tip_point = add(rest_tip, scaled(Point2{.x = target.x - rest_tip.x, .y = target.y - rest_tip.y}, amount));
        // Pressing above full activation pushes the pad into the ridge.
        double press = 0.4 * std::max(tip.activation - 1.0, 0.0);
        tip_point = add(tip_point, scaled(wall_normal, press));
        // A laminal gesture leaves the tip down, its pad against the lingual
        // face of the lower incisors, and makes the constriction with the
        // blade (the area-domain narrowing does that).
        double down = std::clamp(tip.tip_down, 0.0, 1.0);
        if (down > 0.0) {
            double teeth_top = bend_start + bend_radius - (0.35 + jaw);
            Point2 behind_teeth{.x = incisors_x - tip_radius - 0.05, .y = teeth_top - tip_radius};
            Point2 laminal = add(rest_tip, scaled(Point2{.x = behind_teeth.x - rest_tip.x,
                                                         .y = behind_teeth.y - rest_tip.y}, amount));
            tip_point = add(scaled(tip_point, 1.0 - down), scaled(laminal, down));
        }
        // One tongue: a blade too short to reach drags the body after it.
        for (std::size_t iteration = 0; iteration < 4; ++iteration) {
            Point2 reach_vector{.x = tip_point.x - root.x, .y = tip_point.y - root.y};
            double span = norm(reach_vector);
            if (span <= blade_length) {
                break;
            }
            double excess = span - blade_length;
            Point2 pull = scaled(reach_vector, excess / span);
            double fraction=1;
            if (surface) { fraction=surface->admissible_pull(shape.surface_translation,pull);pull=scaled(pull,fraction); }
            centre = add(centre, pull);
            shape.drag += excess*fraction;
            shape.drag_vector = add(shape.drag_vector, pull);
            if (surface!=nullptr) { root=add(root,pull);shape.surface_translation=add(shape.surface_translation,pull); }
            else { root = blade_root(centre); }
            if (fraction<1) {
                const Point2 remaining{tip_point.x-root.x,tip_point.y-root.y};
                const double distance=norm(remaining);
                if (distance>blade_length) {
                    shape.tip_target_error=distance-blade_length;
                    tip_point=add(root,scaled(remaining,blade_length/distance));
                }
                break;
            }
        }
    }
    shape.body_centre = centre;
    if (surface) { shape.body_wall_penetration=surface->wall_penetration(shape.surface_translation); }
    shape.blade_root = root;
    shape.tip = tip_point;
    return shape;
}

void MidsagittalTongue::gaps(const TongueShape& shape, const TonguePosture& posture, std::span<const double> arcs,
                             std::span<double> output) noexcept {
    for (std::size_t index = 0; index < arcs.size(); ++index) {
        double arc = arcs[index];
        Point2 wall_normal = normal(arc);
        double offset = outer_offset(arc);
        Point2 wall = add(midline(arc), scaled(wall_normal, offset));
        Point2 inward = scaled(wall_normal, -1.0);
        double gap = offset + floor_offset(arc, posture);
        const double body_gap=shape.surface?shape.surface->ray_gap(wall,inward,shape.surface_translation):ray_circle(wall,inward,shape.body_centre,body_radius);
        if (std::isnan(body_gap)) { output[index]=body_gap;continue; }
        gap = smooth_minimum(gap, body_gap);
        gap = smooth_minimum(gap, ray_segment(wall, inward, shape.blade_root, shape.tip));
        gap = smooth_minimum(gap, ray_circle(wall, inward, shape.tip, tip_radius));
        output[index] = std::max(gap, 0.0);
    }
}

void MidsagittalTongue::floor_space(const TongueShape& shape, const TonguePosture& posture,
                                    std::span<const double> arcs, std::span<double> upper_gap,
                                    std::span<double> lower_gap, std::span<double> open_above) noexcept {
    (void)upper_gap;
    const double oral_height = bend_start + bend_radius; // y of the oral midline
    const double pad_bottom = shape.tip.y - tip_radius;
    const double pad_front = shape.tip.x + tip_radius;
    const double underside_back = std::min(frenulum_x, shape.tip.x - underside_reach);
    const double spacing = arcs.size() > 1 ? arcs[1] - arcs[0] : 1.0;
    // Pass 1: the height of the space under the tongue at each section
    // centre, and how much of the section the tongue covers (the tip moving
    // over a section uncovers or covers it gradually; a space that appeared
    // whole as the pad passed a section centre was a click, MEASURED).
    for (std::size_t index = 0; index < arcs.size(); ++index) {
        lower_gap[index] = 0.0;
        open_above[index] = 0.0;
        double arc = arcs[index];
        if (arc <= bend_end) {
            continue;
        }
        double x = midline(arc).x;
        double cover = std::clamp((pad_front - x) / spacing + 0.5, 0.0, 1.0);
        if (x < floor_back_x || x >= incisors_x || cover <= 0.0 || x <= underside_back) {
            continue;
        }
        double floor_y = oral_height - floor_offset(arc, posture);
        double place = std::min(x, shape.tip.x);
        // The blade is a slab: its underside lies one thickness below its
        // upper surface (1.2 cm at the root, the pad's diameter at the
        // tip), and the frenulum holds it no higher than the line from the
        // pad down to its contact with the floor. A blade drawn as that line
        // alone was paper thin and left up to 4 cm^2 below it (MEASURED:
        // /sh/ noise peak fell from 2.7 to 1.5 kHz).
        double t = (place - underside_back) / std::max(shape.tip.x - underside_back, 1.0e-6);
        double tether_y = floor_y + std::clamp(t, 0.0, 1.0) * (pad_bottom - floor_y);
        double slab_y = pad_bottom;
        double run = shape.tip.x - shape.blade_root.x;
        if (std::fabs(run) > 1.0e-6) {
            double along = std::clamp((place - shape.blade_root.x) / run, 0.0, 1.0);
            double top_y = shape.blade_root.y + along * (shape.tip.y - shape.blade_root.y);
            double thickness = blade_root_thickness + along * (2.0 * tip_radius - blade_root_thickness);
            slab_y = top_y - thickness;
            if (place < shape.blade_root.x) {
                slab_y = floor_y; // the body rests on the floor
            }
        }
        double height = std::max(std::min(tether_y, slab_y) - floor_y, 0.0);
        lower_gap[index] = height * cover;
        open_above[index] = cover; // temporary: the coverage
    }
    // Pass 2: the space opens forward into the mouth through its front face,
    // as far as the next section forward is not covered.
    double next_cover = 0.0;
    for (std::size_t index = arcs.size(); index-- > 0;) {
        double cover = open_above[index];
        double face = 0.0;
        if (cover > 0.0) {
            face = lower_gap[index] / cover * (1.0 - next_cover);
        }
        open_above[index] = face;
        next_cover = cover;
    }
}

}
