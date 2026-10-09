#pragma once
// Boundary-surface geometry for the sibilant prototype: the front of the
// mouth as planar quadrilateral surfaces, each unfolded into a grid of
// patches (strips along u, joined at their edges along v). SI units.
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace boundary {

struct Vec3 final {
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

[[nodiscard]] Vec3 operator+(const Vec3& a, const Vec3& b) noexcept;
[[nodiscard]] Vec3 operator-(const Vec3& a, const Vec3& b) noexcept;
[[nodiscard]] Vec3 operator*(double scale, const Vec3& a) noexcept;
[[nodiscard]] double dot(const Vec3& a, const Vec3& b) noexcept;
[[nodiscard]] Vec3 cross(const Vec3& a, const Vec3& b) noexcept;
[[nodiscard]] double norm(const Vec3& a) noexcept;

enum class Tissue { tooth, palate, tongue, lip, mouth_floor };

// Wall material. Soft-tissue values follow Ishizaka-French-Flanagan cheek
// measurements (2.1 g/cm^2, 845000 dyn/cm^3); tension and the reflectance
// that stands in for viscothermal boundary-layer loss are calibrated guesses.
struct TissueProperties final {
    bool elastic{false};
    double surface_density{0.0};      // kg/m^2
    double foundation_stiffness{0.0}; // N/m^3
    double tension{0.0};              // N/m, membrane tension within and between strips
    double damping_ratio{0.0};
    double reflectance{1.0};          // pressure reflection per collocation, (0, 1]
    double density_gradient{0.0};     // relative change of surface density from u=0 to u=1
};

[[nodiscard]] TissueProperties tissue_properties(Tissue tissue) noexcept;

// A planar quadrilateral. c0->c1 spans u, c0->c3 spans v. The normal points
// from the air into the wall.
struct Surface final {
    std::string name{};
    std::array<Vec3, 4> corners{};
    Tissue tissue{Tissue::tooth};
    std::size_t cells_u{1};
    std::size_t cells_v{1};
    std::size_t first_patch{0};
    Vec3 normal_into_wall{};
};

// One collocation element. Patch index = first_patch + cell_v * cells_u + cell_u.
struct Patch final {
    std::size_t surface{0};
    std::size_t cell_u{0};
    std::size_t cell_v{0};
    std::array<Vec3, 4> corners{};
    Vec3 centroid{};
    Vec3 normal_into_wall{};
    double area{0.0};
};

enum class SourceKind { pressure_dipole, wall_acceleration };

// The same geometry as a 1-D area function for the control model: cells of
// equal length from the closed (constriction) end to the open (lip) end.
struct TubeProfile final {
    double cell_length{0.0};          // m
    std::vector<double> area{};       // m^2, one per cell
    std::size_t source_face{0};       // face index (0 = closed end) of the series source
};

struct Scene final {
    std::string name{};
    std::vector<Surface> surfaces{};
    std::vector<Patch> patches{};
    std::vector<std::size_t> source_positive{};
    std::vector<std::size_t> source_negative{};
    SourceKind source_kind{SourceKind::pressure_dipole};
    Vec3 microphone{};
    double microphone_distance{0.0};  // m, from the lip exit
    TubeProfile tube{};
};

// Midsagittal-plus-width description in centimetres. x runs forward from the
// constriction exit, y is lateral, z is vertical. The incisors are a slab of
// tooth_thickness at x = teeth_x; the jet leaves through the slit between the
// incisal edges and strikes the lower incisors.
struct SibilantPosture final {
    std::string name{};
    double teeth_x{1.0};
    double tooth_thickness{0.2};
    double lip_length{0.8};
    double cavity_width{2.0};
    double lip_width{2.0};
    double palate_z_back{0.45};
    double palate_z_front{0.7};
    double tongue_z_back{0.25};
    double tongue_z_tip{0.05};
    double tip_x{1.0};                // tongue tip; a sublingual pocket lies between tip_x and teeth_x
    double floor_z{0.05};             // floor of the sublingual pocket
    double upper_edge_z{0.3};
    double lower_edge_z{0.15};
    double lip_top_z{0.45};
    double lip_bottom_z{-0.1};
};

[[nodiscard]] SibilantPosture s_posture();
[[nodiscard]] SibilantPosture sh_posture();

// Builds surfaces, patches, jet source patches, microphone and the 1-D
// profile. patch_size is the target patch edge in metres. Replaces scene.
// face_extent > 0 adds the face as a baffle: a forward-facing frame of that
// width (m) around the lip opening, meshed at face_patch_scale * patch_size.
void build_sibilant_scene(const SibilantPosture& posture, double patch_size, double face_extent,
                          double face_patch_scale, Scene& scene);
// A rigid rectangular duct driven by a piston at x = 0, open at x = length
// unless closed_end caps it (a closed box: no radiation, no exit).
void build_duct_scene(double length, double side, double patch_size, bool closed_end, double flange_extent,
                      Scene& scene);

}
