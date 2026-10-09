#include "geometry.hpp"
#include <algorithm>
#include <cmath>

namespace boundary {

Vec3 operator+(const Vec3& a, const Vec3& b) noexcept {
    Vec3 result{a.x + b.x, a.y + b.y, a.z + b.z};
    return result;
}

Vec3 operator-(const Vec3& a, const Vec3& b) noexcept {
    Vec3 result{a.x - b.x, a.y - b.y, a.z - b.z};
    return result;
}

Vec3 operator*(double scale, const Vec3& a) noexcept {
    Vec3 result{scale * a.x, scale * a.y, scale * a.z};
    return result;
}

double dot(const Vec3& a, const Vec3& b) noexcept {
    double result = a.x * b.x + a.y * b.y + a.z * b.z;
    return result;
}

Vec3 cross(const Vec3& a, const Vec3& b) noexcept {
    Vec3 result{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    return result;
}

double norm(const Vec3& a) noexcept {
    double result = std::sqrt(dot(a, a));
    return result;
}

TissueProperties tissue_properties(Tissue tissue) noexcept {
    TissueProperties properties{};
    switch (tissue) {
    case Tissue::tooth:
        properties.elastic = false;
        properties.reflectance = 0.999;
        break;
    case Tissue::palate:
        // Mucosa on bone: heavy foundation stiffness.
        properties.elastic = true;
        properties.surface_density = 15.0;
        properties.foundation_stiffness = 1.0e8;
        properties.tension = 200.0;
        properties.damping_ratio = 0.3;
        properties.reflectance = 0.995;
        break;
    case Tissue::tongue:
        properties.elastic = true;
        properties.surface_density = 21.0;
        properties.foundation_stiffness = 8.45e6;
        properties.tension = 50.0;
        properties.damping_ratio = 0.4;
        properties.reflectance = 0.995;
        properties.density_gradient = -0.5; // thinner toward the tip
        break;
    case Tissue::lip:
        properties.elastic = true;
        properties.surface_density = 12.0;
        properties.foundation_stiffness = 4.0e6;
        properties.tension = 30.0;
        properties.damping_ratio = 0.3;
        properties.reflectance = 0.995;
        properties.density_gradient = -0.4; // thinner toward the vermilion
        break;
    case Tissue::mouth_floor:
        properties.elastic = true;
        properties.surface_density = 21.0;
        properties.foundation_stiffness = 8.45e6;
        properties.tension = 40.0;
        properties.damping_ratio = 0.4;
        properties.reflectance = 0.995;
        break;
    }
    return properties;
}

namespace {

constexpr double cm = 1.0e-2;

Vec3 point_cm(double x, double y, double z) {
    Vec3 result{x * cm, y * cm, z * cm};
    return result;
}

Vec3 bilinear(const std::array<Vec3, 4>& c, double s, double t) {
    Vec3 bottom = (1.0 - s) * c[0] + s * c[1];
    Vec3 top = (1.0 - s) * c[3] + s * c[2];
    Vec3 result = (1.0 - t) * bottom + t * top;
    return result;
}

std::size_t cell_count(double length, double patch_size) {
    double cells = std::round(length / patch_size);
    std::size_t result = 1;
    if (cells > 1.0) {
        result = static_cast<std::size_t>(cells);
    }
    return result;
}

// Adds a planar quadrilateral subdivided into patches. air_direction is any
// vector pointing from the surface into the air; it orients the normal.
// Degenerate surfaces (an edge shorter than 1 micrometre) are skipped and
// return the current surface count without adding.
std::size_t add_surface(Scene& scene, const std::string& name, const Vec3& c0, const Vec3& c1,
                        const Vec3& c2, const Vec3& c3, const Vec3& air_direction, Tissue tissue,
                        double patch_size) {
    std::size_t index = scene.surfaces.size();
    double length_u = std::max(norm(c1 - c0), norm(c2 - c3));
    double length_v = std::max(norm(c3 - c0), norm(c2 - c1));
    if (length_u < 1.0e-6 || length_v < 1.0e-6) {
        return index;
    }
    Surface surface{};
    surface.name = name;
    surface.corners = {c0, c1, c2, c3};
    surface.tissue = tissue;
    surface.cells_u = cell_count(length_u, patch_size);
    surface.cells_v = cell_count(length_v, patch_size);
    surface.first_patch = scene.patches.size();
    Vec3 normal = cross(c1 - c0, c3 - c0);
    normal = (1.0 / norm(normal)) * normal;
    if (dot(normal, air_direction) > 0.0) {
        normal = -1.0 * normal;
    }
    surface.normal_into_wall = normal;
    for (std::size_t v = 0; v < surface.cells_v; ++v) {
        for (std::size_t u = 0; u < surface.cells_u; ++u) {
            double s0 = static_cast<double>(u) / static_cast<double>(surface.cells_u);
            double s1 = static_cast<double>(u + 1) / static_cast<double>(surface.cells_u);
            double t0 = static_cast<double>(v) / static_cast<double>(surface.cells_v);
            double t1 = static_cast<double>(v + 1) / static_cast<double>(surface.cells_v);
            Patch patch{};
            patch.surface = index;
            patch.cell_u = u;
            patch.cell_v = v;
            patch.corners[0] = bilinear(surface.corners, s0, t0);
            patch.corners[1] = bilinear(surface.corners, s1, t0);
            patch.corners[2] = bilinear(surface.corners, s1, t1);
            patch.corners[3] = bilinear(surface.corners, s0, t1);
            patch.centroid = bilinear(surface.corners, 0.5 * (s0 + s1), 0.5 * (t0 + t1));
            patch.normal_into_wall = normal;
            Vec3 diagonal_a = patch.corners[2] - patch.corners[0];
            Vec3 diagonal_b = patch.corners[3] - patch.corners[1];
            patch.area = 0.5 * norm(cross(diagonal_a, diagonal_b));
            scene.patches.push_back(patch);
        }
    }
    scene.surfaces.push_back(surface);
    return index;
}

// Appends the top row (largest cell_v) of a surface to a source list.
void add_top_row(const Scene& scene, std::size_t surface_index, std::vector<std::size_t>& sources) {
    if (surface_index >= scene.surfaces.size()) {
        return;
    }
    const Surface& surface = scene.surfaces[surface_index];
    std::size_t row = surface.cells_v - 1;
    for (std::size_t u = 0; u < surface.cells_u; ++u) {
        std::size_t patch = surface.first_patch + row * surface.cells_u + u;
        sources.push_back(patch);
    }
}

// A forward-facing frame (air toward +x) in the plane x around the opening
// y in [-half_width, half_width], z in [bottom, top], extent wide on every side.
void add_baffle(Scene& scene, double x, double half_width, double bottom, double top, double extent, Tissue tissue,
                double patch_size) {
    const Vec3 forward{1.0, 0.0, 0.0};
    const double outer = half_width + extent;
    add_surface(scene, "baffle_top", Vec3{x, -outer, top}, Vec3{x, outer, top}, Vec3{x, outer, top + extent},
                Vec3{x, -outer, top + extent}, forward, tissue, patch_size);
    add_surface(scene, "baffle_bottom", Vec3{x, -outer, bottom - extent}, Vec3{x, outer, bottom - extent},
                Vec3{x, outer, bottom}, Vec3{x, -outer, bottom}, forward, tissue, patch_size);
    add_surface(scene, "baffle_left", Vec3{x, -outer, bottom}, Vec3{x, -half_width, bottom}, Vec3{x, -half_width, top},
                Vec3{x, -outer, top}, forward, tissue, patch_size);
    add_surface(scene, "baffle_right", Vec3{x, half_width, bottom}, Vec3{x, outer, bottom}, Vec3{x, outer, top},
                Vec3{x, half_width, top}, forward, tissue, patch_size);
}

double linear(double a, double b, double fraction) {
    double result = a + (b - a) * fraction;
    return result;
}

}

SibilantPosture s_posture() {
    SibilantPosture posture{};
    posture.name = "s";
    posture.teeth_x = 1.0;
    posture.tooth_thickness = 0.2;
    posture.lip_length = 0.8;
    posture.cavity_width = 2.0;
    posture.lip_width = 2.0;
    posture.palate_z_back = 0.45;
    posture.palate_z_front = 0.7;
    posture.tongue_z_back = 0.25;
    posture.tongue_z_tip = 0.05;
    posture.tip_x = 1.0;
    posture.floor_z = 0.05;
    posture.upper_edge_z = 0.3;
    posture.lower_edge_z = 0.15;
    posture.lip_top_z = 0.45;
    posture.lip_bottom_z = -0.1;
    return posture;
}

SibilantPosture sh_posture() {
    SibilantPosture posture{};
    posture.name = "sh";
    posture.teeth_x = 2.0;
    posture.tooth_thickness = 0.2;
    posture.lip_length = 1.2;   // rounded and protruded
    posture.cavity_width = 2.0;
    posture.lip_width = 1.4;
    posture.palate_z_back = 0.7;
    posture.palate_z_front = 0.9;
    posture.tongue_z_back = 0.35;
    posture.tongue_z_tip = 0.25;
    posture.tip_x = 1.4;        // raised tip; sublingual pocket from 1.4 to 2.0
    posture.floor_z = -0.6;
    posture.upper_edge_z = 0.4;
    posture.lower_edge_z = 0.15;
    posture.lip_top_z = 0.55;
    posture.lip_bottom_z = -0.2;
    return posture;
}

void build_sibilant_scene(const SibilantPosture& p, double patch_size, double face_extent, double face_patch_scale,
                          Scene& scene) {
    scene = Scene{};
    scene.name = p.name;
    scene.source_kind = SourceKind::pressure_dipole;
    const double hw = 0.5 * p.cavity_width;
    const double hl = 0.5 * p.lip_width;
    const double xt = p.teeth_x;
    const double xf = p.teeth_x + p.tooth_thickness;
    const double xe = xf + p.lip_length;
    const double xtip = p.tip_x;
    const bool pocket = xtip < xt - 1.0e-9;
    const double palate_at_tip = linear(p.palate_z_back, p.palate_z_front, xtip / xt);
    double floor_at_teeth = p.tongue_z_tip;
    if (pocket) {
        floor_at_teeth = p.floor_z;
    }
    const Vec3 up{0.0, 0.0, 1.0};
    const Vec3 down{0.0, 0.0, -1.0};
    const Vec3 forward{1.0, 0.0, 0.0};
    const Vec3 backward{-1.0, 0.0, 0.0};

    add_surface(scene, "palate", point_cm(0.0, -hw, p.palate_z_back), point_cm(xt, -hw, p.palate_z_front),
                point_cm(xt, hw, p.palate_z_front), point_cm(0.0, hw, p.palate_z_back), down, Tissue::palate,
                patch_size);
    add_surface(scene, "tongue_blade", point_cm(0.0, -hw, p.tongue_z_back), point_cm(xtip, -hw, p.tongue_z_tip),
                point_cm(xtip, hw, p.tongue_z_tip), point_cm(0.0, hw, p.tongue_z_back), up, Tissue::tongue,
                patch_size);
    add_surface(scene, "constriction_face", point_cm(0.0, -hw, p.tongue_z_back), point_cm(0.0, hw, p.tongue_z_back),
                point_cm(0.0, hw, p.palate_z_back), point_cm(0.0, -hw, p.palate_z_back), forward, Tissue::tongue,
                patch_size);
    if (pocket) {
        add_surface(scene, "tongue_tip_underside", point_cm(xtip, -hw, p.floor_z), point_cm(xtip, hw, p.floor_z),
                    point_cm(xtip, hw, p.tongue_z_tip), point_cm(xtip, -hw, p.tongue_z_tip), forward,
                    Tissue::tongue, patch_size);
        add_surface(scene, "sublingual_floor", point_cm(xtip, -hw, p.floor_z), point_cm(xt, -hw, p.floor_z),
                    point_cm(xt, hw, p.floor_z), point_cm(xtip, hw, p.floor_z), up, Tissue::mouth_floor,
                    patch_size);
    }
    for (std::size_t side = 0; side < 2; ++side) {
        double sign = 1.0;
        if (side == 1) {
            sign = -1.0;
        }
        double y = sign * hw;
        Vec3 inward{0.0, -sign, 0.0};
        add_surface(scene, "side_wall_front", point_cm(0.0, y, p.tongue_z_back), point_cm(xtip, y, p.tongue_z_tip),
                    point_cm(xtip, y, palate_at_tip), point_cm(0.0, y, p.palate_z_back), inward, Tissue::tooth,
                    patch_size);
        if (pocket) {
            add_surface(scene, "side_wall_pocket", point_cm(xtip, y, p.floor_z), point_cm(xt, y, p.floor_z),
                        point_cm(xt, y, p.palate_z_front), point_cm(xtip, y, palate_at_tip), inward,
                        Tissue::tooth, patch_size);
        }
        // Closed dental arch beside a lip opening narrower than the cavity.
        add_surface(scene, "teeth_side_panel", point_cm(xt, sign * hl, floor_at_teeth), point_cm(xt, y, floor_at_teeth),
                    point_cm(xt, y, p.palate_z_front), point_cm(xt, sign * hl, p.palate_z_front), backward,
                    Tissue::tooth, patch_size);
        Vec3 lip_inward{0.0, -sign, 0.0};
        add_surface(scene, "incisal_gap_side", point_cm(xt, sign * hl, p.lower_edge_z), point_cm(xf, sign * hl, p.lower_edge_z),
                    point_cm(xf, sign * hl, p.upper_edge_z), point_cm(xt, sign * hl, p.upper_edge_z), lip_inward,
                    Tissue::tooth, patch_size);
        add_surface(scene, "lip_corner", point_cm(xf, sign * hl, p.lip_bottom_z), point_cm(xe, sign * hl, p.lip_bottom_z),
                    point_cm(xe, sign * hl, p.lip_top_z), point_cm(xf, sign * hl, p.lip_top_z), lip_inward,
                    Tissue::lip, patch_size);
    }
    add_surface(scene, "upper_incisor_back", point_cm(xt, -hl, p.upper_edge_z), point_cm(xt, hl, p.upper_edge_z),
                point_cm(xt, hl, p.palate_z_front), point_cm(xt, -hl, p.palate_z_front), backward, Tissue::tooth,
                patch_size);
    std::size_t lower_back = add_surface(scene, "lower_incisor_back", point_cm(xt, -hl, floor_at_teeth),
                                         point_cm(xt, hl, floor_at_teeth), point_cm(xt, hl, p.lower_edge_z),
                                         point_cm(xt, -hl, p.lower_edge_z), backward, Tissue::tooth, patch_size);
    add_surface(scene, "upper_incisal_edge", point_cm(xt, -hl, p.upper_edge_z), point_cm(xf, -hl, p.upper_edge_z),
                point_cm(xf, hl, p.upper_edge_z), point_cm(xt, hl, p.upper_edge_z), down, Tissue::tooth, patch_size);
    add_surface(scene, "lower_incisal_edge", point_cm(xt, -hl, p.lower_edge_z), point_cm(xf, -hl, p.lower_edge_z),
                point_cm(xf, hl, p.lower_edge_z), point_cm(xt, hl, p.lower_edge_z), up, Tissue::tooth, patch_size);
    add_surface(scene, "upper_incisor_front", point_cm(xf, -hl, p.upper_edge_z), point_cm(xf, hl, p.upper_edge_z),
                point_cm(xf, hl, p.lip_top_z), point_cm(xf, -hl, p.lip_top_z), forward, Tissue::tooth, patch_size);
    std::size_t lower_front = add_surface(scene, "lower_incisor_front", point_cm(xf, -hl, p.lip_bottom_z),
                                          point_cm(xf, hl, p.lip_bottom_z), point_cm(xf, hl, p.lower_edge_z),
                                          point_cm(xf, -hl, p.lower_edge_z), forward, Tissue::tooth, patch_size);
    add_surface(scene, "upper_lip", point_cm(xf, -hl, p.lip_top_z), point_cm(xe, -hl, p.lip_top_z),
                point_cm(xe, hl, p.lip_top_z), point_cm(xf, hl, p.lip_top_z), down, Tissue::lip, patch_size);
    add_surface(scene, "lower_lip", point_cm(xf, -hl, p.lip_bottom_z), point_cm(xe, -hl, p.lip_bottom_z),
                point_cm(xe, hl, p.lip_bottom_z), point_cm(xf, hl, p.lip_bottom_z), up, Tissue::lip, patch_size);

    // Turbulent jet striking the lower incisors: a fluctuating force on the
    // tooth, i.e. equal and opposite surface pressures on its two faces.
    add_top_row(scene, lower_back, scene.source_positive);
    add_top_row(scene, lower_front, scene.source_negative);

    if (face_extent > 0.0) {
        add_baffle(scene, xe * cm, hl * cm, p.lip_bottom_z * cm, p.lip_top_z * cm, face_extent, Tissue::lip,
                   face_patch_scale * patch_size);
    }

    const double lip_mid_z = 0.5 * (p.lip_top_z + p.lip_bottom_z);
    scene.microphone_distance = 0.15;
    scene.microphone = point_cm(xe + 15.0, 0.0, lip_mid_z);

    // 1-D control: area function from the same numbers. The sublingual
    // pocket is simply added cross-section, which is all a 1-D tube can do.
    const double dx = 0.2;
    std::size_t cells = static_cast<std::size_t>(std::round(xe / dx));
    scene.tube.cell_length = dx * cm;
    scene.tube.area.assign(cells, 0.0);
    for (std::size_t i = 0; i < cells; ++i) {
        double x = (static_cast<double>(i) + 0.5) * dx;
        double area_cm2 = 0.0;
        if (x < xt) {
            double palate = linear(p.palate_z_back, p.palate_z_front, x / xt);
            double floor = p.floor_z;
            if (x < xtip) {
                floor = linear(p.tongue_z_back, p.tongue_z_tip, x / xtip);
            }
            area_cm2 = p.cavity_width * (palate - floor);
        } else if (x < xf) {
            area_cm2 = p.lip_width * (p.upper_edge_z - p.lower_edge_z);
        } else {
            area_cm2 = p.lip_width * (p.lip_top_z - p.lip_bottom_z);
        }
        scene.tube.area[i] = area_cm2 * cm * cm;
    }
    scene.tube.source_face = static_cast<std::size_t>(std::round(xf / dx));
}

void build_duct_scene(double length, double side, double patch_size, bool closed_end, double flange_extent,
                      Scene& scene) {
    scene = Scene{};
    scene.name = "duct";
    if (closed_end) {
        scene.name = "box";
    }
    scene.source_kind = SourceKind::wall_acceleration;
    const double h = 0.5 * side;
    const Vec3 up{0.0, 0.0, 1.0};
    const Vec3 down{0.0, 0.0, -1.0};
    const Vec3 forward{1.0, 0.0, 0.0};
    Vec3 a{0.0, -h, 0.0};
    Vec3 b{0.0, h, 0.0};
    Vec3 c{0.0, h, side};
    Vec3 d{0.0, -h, side};
    Vec3 along{length, 0.0, 0.0};
    std::size_t piston = add_surface(scene, "piston", a, b, c, d, forward, Tissue::tooth, patch_size);
    add_surface(scene, "floor", a, a + along, b + along, b, up, Tissue::tooth, patch_size);
    add_surface(scene, "ceiling", d, d + along, c + along, c, down, Tissue::tooth, patch_size);
    Vec3 inward_left{0.0, 1.0, 0.0};
    Vec3 inward_right{0.0, -1.0, 0.0};
    add_surface(scene, "wall_left", a, a + along, d + along, d, inward_left, Tissue::tooth, patch_size);
    add_surface(scene, "wall_right", b, b + along, c + along, c, inward_right, Tissue::tooth, patch_size);
    if (closed_end) {
        Vec3 backward{-1.0, 0.0, 0.0};
        add_surface(scene, "end_cap", a + along, b + along, c + along, d + along, backward, Tissue::tooth, patch_size);
    } else if (flange_extent > 0.0) {
        add_baffle(scene, length, h, 0.0, side, flange_extent, Tissue::tooth, 2.0 * patch_size);
    }
    const Surface& source_surface = scene.surfaces[piston];
    std::size_t count = source_surface.cells_u * source_surface.cells_v;
    for (std::size_t k = 0; k < count; ++k) {
        scene.source_positive.push_back(source_surface.first_patch + k);
    }
    scene.microphone_distance = 0.15;
    scene.microphone = Vec3{length + 0.15, 0.0, h};
    const double dx = 0.002;
    std::size_t cells = static_cast<std::size_t>(std::round(length / dx));
    scene.tube.cell_length = dx;
    scene.tube.area.assign(cells, side * side);
    scene.tube.source_face = 0;
}

}
