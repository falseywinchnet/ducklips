#pragma once
#include "ducklips/acoustic_network.hpp"
#include "ducklips/articulation_speeds.hpp"
#include "ducklips/tongue.hpp"
#include "ducklips/tongue_surface.hpp"
#include "ducklips/voice.hpp"
#include "ducklips/lip_profile.hpp"
#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {

// A vowel posture is a tongue position: body centre x, body centre y (cm,
// reference midsagittal frame), jaw opening (0..1), blade raising (-1..1)
// and tongue-root advancement (-1..1), plus the lips. Lip
// values are reference (17.5 cm tract) values; the model scales them.
inline constexpr std::size_t tongue_parameter_count = 5;
inline constexpr std::size_t anatomy_parameter_count = 4;
inline constexpr std::size_t shape_parameter_count = tongue_parameter_count + anatomy_parameter_count;
inline constexpr std::size_t body_x_parameter = 0;
inline constexpr std::size_t body_y_parameter = 1;
inline constexpr std::size_t jaw_parameter = 2;
inline constexpr std::size_t blade_parameter = 3;
inline constexpr std::size_t root_parameter = 4;
// Parameters 5..8: smooth log-area residuals (tapered cosine modes 1..4 over
// the tongue region) standing in for anatomy the coarse tongue geometry
// does not capture. Penalized so the tongue explains what it can.

struct VowelShape final {
    std::array<double, shape_parameter_count> coefficients{2.8, 6.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double lip_area{3.0e-4};   // m^2
    double lip_protrusion{0.0}; // m
};

enum class ConstrictionSlot : std::size_t { lips = 0, tongue_tip = 1, tongue_body = 2 };
inline constexpr std::size_t constriction_slot_count = 3;

// A local narrowing laid over the vowel posture. Places and widths are
// distances along the reference tract measured from the lip opening.
struct ConstrictionGesture final {
    double place_from_lips{};       // m
    double width{0.005};            // m, Gaussian half-width beyond the flat core
    double target_area{};           // m^2 at full activation; 0 closes
    double obstacle_from_lips{-1.0}; // m; >= 0 when the jet strikes an obstacle (teeth)
    double tip_down{0.0};            // tongue tip: 0 apical, 1 laminal (ConstrictionSpec::tip_down)
};

// Third-order critically damped approach to a target. A first-order
// "command" stage feeds a critically damped second-order body, all three
// poles at the same rate. From rest the motion starts with zero velocity and
// zero acceleration (a small ramp-up), never overshoots, and eases into the
// target (a small ramp-down). angular_rate is what a cue asks for;
// speed_limit is what the tissue allows (0 = unlimited); the slower wins.
struct CriticallyDampedTracker final {
    double value{};
    double velocity{};
    double command{};
    double target{};
    double angular_rate{40.0}; // rad/s
    double speed_limit{};      // rad/s

    void advance(double step) noexcept;
    void place(double position) noexcept;
};

// A constriction slot is one articulator (lips, tongue tip, tongue body). Its
// place, target area and width move continuously when a new gesture arrives
// while the slot is still active, so a tongue tip slides from /l/ to /d/
// instead of teleporting; an idle slot takes a new gesture at once.
struct ConstrictionSlotState final {
    CriticallyDampedTracker activation{};
    CriticallyDampedTracker place{};       // m from the lips, reference scale
    CriticallyDampedTracker target_area{}; // m^2, reference scale
    CriticallyDampedTracker width{};       // m
    CriticallyDampedTracker tip_down{};    // 0 apical .. 1 laminal
    double obstacle_from_lips{-1.0};

    void begin(const ConstrictionGesture& gesture, double activation_target, double angular_rate) noexcept;
    [[nodiscard]] ConstrictionGesture current() const noexcept;
    void advance(double step) noexcept;
};

struct ArticulatorState final {
    std::array<CriticallyDampedTracker, shape_parameter_count> shape{};
    CriticallyDampedTracker lip_area{};
    CriticallyDampedTracker lip_protrusion{};
    std::array<ConstrictionSlotState, constriction_slot_count> slot{};
    // The tongue body's lingering displacement after the blade dragged it
    // (cm, reference frame). It follows the pull with the body's own speed,
    // so after a /d/ release the body drifts back instead of snapping.
    CriticallyDampedTracker body_offset_x{};
    CriticallyDampedTracker body_offset_y{};
    CriticallyDampedTracker port_area{};          // m^2, reference scale
    CriticallyDampedTracker glottal_rest_area{};  // m^2
    CriticallyDampedTracker lung_pressure{};      // Pa
    CriticallyDampedTracker pitch{};              // Hz
    CriticallyDampedTracker effort{};             // lung pressure multiplier
    CriticallyDampedTracker loudness{};           // output gain

    void advance(double step) noexcept;
    void advance_geometry(double step) noexcept;
    void advance_physiology(double step) noexcept;
    void apply_speeds(const ArticulatorSpeeds& speeds) noexcept;
    void place_vowel(const VowelShape& vowel) noexcept;
    void target_vowel(const VowelShape& vowel, double angular_rate) noexcept;
};

// Turns articulator positions into a TractShape for one voice.
//
// The epilarynx is a fixed narrow tube. From its top to the lips the
// cross-sections come from the midsagittal tongue (tongue.hpp): body and jaw
// set the vowel; the tongue-tip slot places the tip (the body follows when
// the blade cannot reach) and the tongue-body slot presses the dorsum against
// the palate. The lips are their own tube, narrowed by the lips slot.
class TractModel final {
public:
    void prepare(const VoiceProfile& voice, std::size_t oral_sections, std::size_t nasal_sections);
    // Sizes a shape for this model. Allocates.
    void allocate(TractShape& shape) const;
    void build(const ArticulatorState& state, TractShape& shape) noexcept;
    void build_vowel(const VowelShape& vowel, TractShape& shape) noexcept;

    [[nodiscard]] std::size_t oral_sections() const noexcept;
    [[nodiscard]] std::size_t nasal_sections() const noexcept;
    [[nodiscard]] std::size_t port_section() const noexcept;
    [[nodiscard]] std::size_t epilarynx_sections() const noexcept;
    [[nodiscard]] std::size_t lip_sections() const noexcept;
    [[nodiscard]] LipProfile lip_profile(const ArticulatorState& state) const noexcept;
    [[nodiscard]] double length_scale() const noexcept; // voice length / reference length
    [[nodiscard]] double area_scale() const noexcept;   // length_scale^2
    // d ln(area_i) / d(posture parameter) for every oral section at a vowel
    // posture, by central differences of the geometry. Layout
    // [section * shape_parameter_count + parameter]. Uses internal scratch.
    void log_area_gradient(const VowelShape& vowel, std::span<double> output);
    [[nodiscard]] const TongueShape& tongue() const noexcept;
    // Body displacement the current tip gesture demands of the posture alone
    // (cm); the synthesizer feeds it to the body-offset trackers.
    [[nodiscard]] Point2 required_drag() const noexcept;
    // Distance of each section centre from the lip opening after the last build.
    [[nodiscard]] double distance_from_lips(std::size_t section) const noexcept;
    // Narrowest oral section beyond the epilarynx after the last build.
    [[nodiscard]] std::size_t narrowest_section() const noexcept;
    [[nodiscard]] double narrowest_area() const noexcept;
    // Section nearest to a reference distance from the lips, after the last build.
    [[nodiscard]] std::size_t section_at_distance(double reference_distance) const noexcept;
    [[nodiscard]] double oral_volume() const noexcept;
    [[nodiscard]] std::span<const SideCavity> side_cavities() const noexcept;
    // With a floor, the tongue is a body in the mouth: build() also writes
    // the space under it (TractShape::sublingual_area, sublingual_opening)
    // and splits the open column in front of the tip. Off by default.
    void set_floor(bool enabled) noexcept;
    void set_conservative_body(bool enabled) noexcept { conservative_body_=enabled; }
    [[nodiscard]] const char* geometry_error() const noexcept { return geometry_error_; }
    // The oral sections that can hold space under the tongue (between
    // MidsagittalTongue::floor_back_x and incisors_x).
    void floor_sections(std::size_t& first, std::size_t& count) const noexcept;

private:
    TongueShape pose(const TonguePosture& posture,const TongueGesture& body,const TongueGesture& tip) noexcept;
    void write_posture(const TonguePosture& posture, const TongueGesture& body, const TongueGesture& tip,
                       std::span<const double> residual, double lip_area, double lip_protrusion, bool vowel_floor,
                       TractShape& shape) noexcept;
    void apply_constriction(const ConstrictionGesture& gesture, double activation, TractShape& shape) noexcept;
    void finish(TractShape& shape) noexcept;

    std::size_t oral_sections_{};
    std::size_t nasal_sections_{};
    std::size_t port_section_{};
    std::size_t epilarynx_sections_{};
    std::size_t lip_sections_{};
    double length_scale_{1.0};
    double area_scale_{1.0};
    double section_length_{};
    std::vector<double> neutral_log_area_{}; // per oral section, reference areas (epilarynx)
    std::vector<double> arc_{};              // reference arc position of each oral section, cm
    std::vector<double> residual_basis_{};   // [section * anatomy_parameter_count + m]
    std::vector<double> gap_{};              // scratch gaps, cm
    bool floor_{false};
    bool conservative_body_{false};
    const char* geometry_error_{nullptr};
    ConservativeTongueSurface body_surface_{};
    std::vector<double> lower_gap_{};        // scratch, cm below the tongue
    std::vector<double> open_above_{};       // scratch, 1 where the space below meets the channel above
    TongueShape tongue_{};
    Point2 required_drag_{};
    TractShape scratch_{};
    std::vector<double> nasal_area_{};
    std::vector<double> nasal_length_{};
    std::vector<SideCavity> side_cavities_{};
    std::vector<double> side_area_{};
    std::vector<double> side_length_{};
    std::vector<double> distance_{};
    std::size_t narrowest_{};
    double narrowest_area_{};
    double oral_volume_{};
};

}
