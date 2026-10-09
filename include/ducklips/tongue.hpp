#pragma once
#include <cstddef>
#include <span>

namespace ducklips {
class ConservativeTongueSurface;

// Midsagittal geometry of the reference (adult male, 17.5 cm) vocal tract, in
// centimetres. x points forward, y up; the glottis is at the origin.
//
// The airway midline rises 7 cm through the pharynx, turns forward through a
// quarter circle of radius 2 cm under the soft palate, and runs forward under
// the hard palate to the lips. Each midline point has a normal pointing to
// the fixed outer wall (posterior pharyngeal wall, velum, palate, alveolar
// ridge). The inner side is the tongue: a body (a circle of tissue) with a
// blade of fixed length ending in the tip, all lowered by the jaw.
//
// A section's cross-sectional area comes from the gap between the outer wall
// and the first inner surface along the normal (Heinz-Stevens style power
// law, a modelling choice). Consonants place the tip or press the body
// against the wall; if the blade cannot reach, the body is dragged along.
struct Point2 final {
    double x{};
    double y{};
};

struct TonguePosture final {
    double body_x{2.8}; // cm, body centre before the jaw lowers it
    double body_y{6.0}; // cm
    double jaw{0.3};    // 0 closed .. 1 wide open
    double blade{0.0};  // -1 low .. +1 raised toward the palate (tip angle -70..+30 degrees)
    double root{0.0};   // -1 retracted .. +1 advanced tongue root (pharynx width)
    double offset_x{};  // cm, dynamic displacement of the body (the lingering drag)
    double offset_y{};  // cm
};

// A tip or body gesture in reference units.
struct TongueGesture final {
    double activation{};      // 0 rest, 1 at the target, >1 pressing into the wall
    double place_from_lips{}; // cm along the airway from the lip opening
    double gap{};             // cm target gap to the wall; 0 seals
    double tip_down{};        // tip gestures: 0 the tip goes to the place, 1 it stays behind the lower incisors
    double width{0.9};        // cm, shoulder extent for conservative body gestures
};

struct TongueShape final {
    const ConservativeTongueSurface* surface{}; // borrowed until the owner rebuilds its next pose
    Point2 surface_translation{};
    double body_wall_penetration{}; // cm, sampled conservative body contour
    double tip_target_error{}; // cm, unreachable target after body contact and blade reach constraints
    bool valid{true};
    Point2 body_centre{};
    Point2 blade_root{};
    Point2 tip{};
    double drag{};        // cm the body was pulled by the blade constraint
    Point2 drag_vector{}; // cm, the same pull as a displacement
};

class MidsagittalTongue final {
public:
    static constexpr double length = 17.5;          // cm, glottis to lip opening
    static constexpr double pharynx_height = 7.0;   // cm
    static constexpr double bend_radius = 2.0;      // cm
    static constexpr double body_radius = 2.7;      // cm
    static constexpr double blade_length = 2.2;     // cm
    static constexpr double tip_radius = 0.35;      // cm, the rounded pad of tissue at the tip

    [[nodiscard]] static Point2 midline(double arc) noexcept;
    [[nodiscard]] static Point2 normal(double arc) noexcept;
    [[nodiscard]] static double outer_offset(double arc) noexcept;
    [[nodiscard]] static double floor_offset(double arc, const TonguePosture& posture) noexcept;
    // Area (cm^2) of a gap (cm) at a place on the midline, and the inverse.
    [[nodiscard]] static double area_from_gap(double gap, double arc) noexcept;
    [[nodiscard]] static double gap_from_area(double area, double arc) noexcept;

    // Positions the tongue for a posture and its gestures.
    [[nodiscard]] static TongueShape pose(const TonguePosture& posture, const TongueGesture& body,
                                          const TongueGesture& tip, const ConservativeTongueSurface* surface=nullptr) noexcept;
    // Writes the gap (cm) at each arc position (cm from the glottis).
    static void gaps(const TongueShape& shape, const TonguePosture& posture, std::span<const double> arcs,
                     std::span<double> output) noexcept;

    // The floor of the mouth under the tongue. The tongue is a body in the
    // mouth, not its bottom: the floor (lowered by the jaw) lies beneath it
    // and the lower incisors close the space in front. The underside of the
    // blade runs from the bottom of the tip pad back and down to the
    // frenulum, which ties it to the floor; behind that the tongue rests on
    // the floor. Where the blade is lifted there is air below it.
    static constexpr double frenulum_x = 5.8;   // cm, where the underside meets the floor at most
    static constexpr double underside_reach = 0.8; // cm, shortest sloping underside behind the tip
    static constexpr double incisors_x = 8.0;   // cm, inner face of the lower incisors
    static constexpr double floor_back_x = 4.0; // cm, no air under the tongue behind this
    static constexpr double blade_root_thickness = 1.2; // cm, the blade slab at its root
    // For each arc position: the gap below the tongue down to the floor (cm,
    // 0 where the tongue rests on the floor or there is no floor space), and
    // the height (cm) of the face through which that space opens forward
    // into the mouth (at the last position under the tongue, 0 elsewhere).
    // upper_gap is left as it is.
    static void floor_space(const TongueShape& shape, const TonguePosture& posture, std::span<const double> arcs,
                            std::span<double> upper_gap, std::span<double> lower_gap,
                            std::span<double> open_above) noexcept;
};

}
