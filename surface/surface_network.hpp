#pragma once
#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {

struct Vec3 final {
    double x{};
    double y{};
    double z{};
};

enum class SurfaceRole { wall, glottis, aperture };

// The material of a wall piece, per unit area: how heavy, how springy and how
// lossy the tissue is.
struct WallMaterial final {
    double mass{15.0};          // kg/m^2
    double resistance{1.0e4};   // Pa s/m
    double stiffness{5.3e5};    // Pa/m (resonance near 30 Hz with the default mass)
};

// One piece of the cavity surface. It collects the fluctuations arriving
// through the air from the pieces it sees, scatters them passively, and
// emits the result back along its paths. A yielding piece also moves under
// the pressure according to its material, absorbing some energy and
// re-emitting the rest with the delay of its own motion.
struct SurfacePiece final {
    Vec3 position{};    // m
    Vec3 normal{};      // unit, into the cavity
    SurfaceRole role{SurfaceRole::wall};
    std::size_t station{}; // flat station index (all tubes)
    std::size_t tube{};
    WallMaterial material{};
    bool yielding{true};
};

// A one-way air path. The fluctuation emitted by one piece arrives at the
// other after the travel time: an integer delay plus a first-order all-pass
// for the fraction, which keeps every frequency at full strength.
struct AirLink final {
    std::size_t from{};
    std::size_t to{};
    std::size_t whole_delay{};
    double allpass{};    // (1 - f) / (1 + f), fraction f in [0.5, 1.5)
    double admittance{}; // m^3/(Pa s): the air the path carries
    double scale{};      // admittance per unit crossing area
    std::size_t buffer_offset{};
};

inline constexpr std::size_t no_parent = static_cast<std::size_t>(-1);

enum class TubeEnd { closed, radiating };

// One tube of the cavity surface. Tube 0 is the airway from the glottis to
// the lips; every other tube opens from a station of an earlier tube (the
// nose from the velar port, the piriform fossae beside the larynx, the
// sinuses from the nose). Its stations run from that opening outward.
struct TubeGeometry final {
    std::size_t stations{};
    double length{};                 // m
    std::size_t ring_pieces{6};
    std::size_t parent{no_parent};
    std::size_t parent_station{};    // station index within the parent tube
    TubeEnd end{TubeEnd::closed};
    WallMaterial material{};
    bool yielding{true};
    // Wetted perimeter relative to a circle of the same area. Vocal-tract
    // cross-sections are flattened (ellipses, slits) and the nose is lined
    // with turbinates, so their walls rub more air per unit area than a
    // circle's (HYPOTHESIS: 1.5 for the airway).
    double perimeter_factor{1.5};
    // Share of that perimeter that yields (soft tissue rather than palate,
    // teeth or cartilage). 0.67 of the flattened perimeter is the yielding
    // area of a circle of the same section, which puts the closed-tract
    // wall resonance near 180 Hz as measured (MEASURED here: /u/ F1 310 Hz
    // rigid, 360 Hz yielding).
    double yielding_fraction{0.67};

    // The tube in space. A plain tube is straight with its rings stacked
    // along it. A solid tube places every piece of every ring on its
    // cross-section: an ellipse against the outer wall (the back of the
    // pharynx, the palate), as wide as the anatomy and as deep as the area
    // then requires, so a constriction is a slit along that wall. A bent
    // solid tube is the airway: it rises from the glottis, turns under the
    // palate and runs forward to the lips; its inner side is shorter.
    bool solid{false};
    bool bent{false};
    double width{0.0}; // m, a straight solid tube's width; 0: a circle of its area
    // Air paths across each ring, between every pair of its pieces that
    // face each other: what carries the cavity's cross modes.
    bool cross_paths{false};
    // With cross paths the air of each ring is shared: this fraction of its
    // compliance travels along the tube, the rest across it. Every path's
    // travel time is then shortened by the square root of its direction's
    // share, so sound still crosses the air at its own speed in every
    // direction (the scaling of a waveguide mesh). 1 without cross paths.
    double axial_share{1.0};
    // Travel time of the cross paths relative to the mesh scaling: with
    // boundary pieces only, the lumped share of the ring's air at the walls
    // lowers the cross modes (MEASURED: a 3.57 cm circle's first mode at
    // 4.95 kHz against 5.75 kHz analytic at 1.0).
    double cross_timing{1.0};
    // The material along the surface: neighbouring pieces pull on each
    // other through their shared edges (N/m, a membrane tension), so a
    // fluctuation also travels through the tissue.
    double edge_tension{0.0};
    // A parallel tube lies beside stations parent_station ... of its parent,
    // one of its stations beside each (the space under the tongue beside the
    // mouth above it). It has no opening at its first station; it meets the
    // parent wherever set_openings() opens the way between a station and the
    // one beside it.
    bool parallel{false};
};

struct CavityGeometry final {
    std::vector<TubeGeometry> tubes{};
    std::vector<double> station_area{}; // m^2, the stations of every tube in tube order
};

struct AirConstants final {
    double density{1.14};
    double sound_speed{350.0};
    double viscosity{1.86e-5};      // Pa s
    // Heat exchange with the wall adds (gamma - 1)/sqrt(Prandtl) of the
    // viscous boundary-layer loss.
    double heat_loss_factor{0.475};
};

// Offline, frozen, zero-mean-flow linearization. A directed pressure-wave path
// has transfer q^delay * (b0+b1*q+b2*q*q)/(1+a1*q+a2*q*q), q=z^-1.
// Power-normalized paths have the same pressure transfer at fixed admittance.
struct SurfaceLinearPath final {
    std::size_t from{};
    std::size_t to{};
    std::size_t delay{};
    double admittance{};
    std::array<double, 3> numerator{};
    std::array<double, 3> denominator{1.0, 0.0, 0.0};
};

struct SurfaceLinearPort final {
    std::size_t piece{};
    double admittance{};
    std::array<double, 3> numerator{}; // reflection filter
    std::array<double, 3> denominator{1.0, 0.0, 0.0};
    bool radiation{};
};

struct SurfaceLinearSnapshot final {
    double sample_rate{};
    std::size_t pieces{};
    std::size_t glottis{};
    std::size_t allocated_delay_samples{};
    std::vector<std::size_t> piece_tubes{}; // anatomical partition, including terminal pieces
    std::vector<SurfaceLinearPath> paths{}; // reciprocal paths remain adjacent
    std::vector<SurfaceLinearPort> ports{};
};

// A single airway tube with the given areas (no branches), open at the lips.
[[nodiscard]] CavityGeometry single_tube(std::span<const double> station_area, double length, std::size_t ring_pieces);

// The cavity as a passive scattering network of surface pieces linked
// through air, in physical units: pressures in Pa, flows in m^3/s,
// admittances in m^3/(Pa s).
//
// Per sample: gather() reads every arriving fluctuation and reports the
// pressure the glottis piece would have with no flow of its own, and its
// driving-point impedance; the caller solves the glottal flow against it;
// scatter() injects that flow and the turbulent pressures, forms every
// piece's pressure, and emits.
//
// Kinetic losses. Where the airway widens in the direction of the flow, the
// air leaves the narrower section as a jet and its kinetic energy is lost in
// the expansion (Borda-Carnot): a pressure drop
//   rho U |U| / 2 (1/A_narrow - 1/A_wide)^2
// in series at that boundary, and at an opening into the free air
// rho U |U| / (2 A^2). In wave terms the resistance is a small fraction of
// the boundary's characteristic admittance (about half the jet Mach number):
// negligible for sound, decisive for the steady flow, so the pressure behind
// a fricative or a closing stop builds as in the body. It varies smoothly
// with the areas, so a moving constriction never relocates a pressure drop
// in one step. The turbulence of the main jet enters as a series pressure
// at the boundary where it leaves the narrowest section.
class SurfaceNetwork final {
public:
    void build(const CavityGeometry& geometry, double sample_rate, const AirConstants& air);
    void reset() noexcept;
    // Reshapes in place: the pieces and paths stay; path weights follow the
    // new cross-sections, gliding there linearly over glide_samples samples
    // (a stepped change every block would itself be a tone at the block rate).
    void reshape(std::span<const double> station_area, std::size_t glide_samples = 1) noexcept;
    // Turbulence of the main jet enters where the air leaves this station
    // (flat index) downstream; no_parent turns it off. A moved source fades
    // in over fade_samples while the old place fades out: a noise pressure
    // that jumped places would step, and a step is a click.
    void set_turbulence_exit(std::size_t station, std::size_t fade_samples = 1) noexcept;
    // Openings (m^2, per flat station; only stations of parallel tubes are
    // read) between a parallel tube's station and the parent station beside
    // it. Takes effect at the next reshape, gliding with it.
    void set_openings(std::span<const double> opening_area) noexcept;

    void gather() noexcept;
    [[nodiscard]] double glottis_free_pressure() const noexcept; // Pa
    [[nodiscard]] double glottis_impedance() const noexcept;     // Pa s/m^3
    // Returns the microphone signal: the rate of change of the flow radiated
    // through every opening (m^3/s^2).
    double scatter(double glottal_flow, double glottal_turbulence, double jet_turbulence) noexcept;

    // Flow (m^3/s, last sample) leaving a station downstream: into the next
    // station of its tube, or out through its opening.
    [[nodiscard]] double exit_flow(std::size_t station) const noexcept;
    // Diagnostics: an opening's gliding area and the resistance R Y of its
    // kinetic loss.
    [[nodiscard]] double opening_area(std::size_t opening) const noexcept;
    [[nodiscard]] double opening_resistance(std::size_t opening) const noexcept;
    // R Y of the boundary into a station from its upstream station.
    [[nodiscard]] double entry_resistance(std::size_t station) const noexcept;
    // Volume flow of every yielding wall moving outward (m^3/s, last
    // sample): what the neck and cheeks pass to the air during a closure.
    [[nodiscard]] double wall_flow() const noexcept;
    [[nodiscard]] double station_pressure(std::size_t station) const noexcept; // Pa, mean over its ring
    [[nodiscard]] double radiated_flow(std::size_t opening) const noexcept;    // m^3/s, last sample
    [[nodiscard]] std::size_t opening_count() const noexcept;
    [[nodiscard]] std::size_t station_count() const noexcept;
    [[nodiscard]] std::size_t tube_station(std::size_t tube, std::size_t station) const noexcept;
    [[nodiscard]] std::size_t piece_count() const noexcept;
    [[nodiscard]] std::size_t link_count() const noexcept;
    // Diagnostics: a flow (m^3/s) injected into one piece at the next
    // scatter; one piece's pressure (Pa, last sample); a station's first piece.
    void add_flow(std::size_t piece, double flow) noexcept;
    [[nodiscard]] double piece_pressure(std::size_t piece) const noexcept;
    [[nodiscard]] std::size_t first_piece(std::size_t station) const noexcept;
    // Allocates offline. Rejects moving geometry, nonzero mean-flow kinetic
    // resistance, and membrane edge tension (not represented by these ports).
    [[nodiscard]] SurfaceLinearSnapshot linear_snapshot() const;

private:
    struct Station final {
        std::size_t tube{};
        std::size_t first_piece{};
        std::size_t pieces{};
        std::size_t upstream{no_parent};   // flat station it opens from; no_parent: the glottis or nothing
        std::size_t downstream{no_parent}; // next station in the same tube
        double spacing{};                  // m
        double perimeter_factor{1.0};
        double yielding_fraction{1.0};
        // Solid tubes: the outer wall's point and inward direction at this
        // station, the cross-section's half width, and the outer wall's
        // radius in a bend (0 when straight): a path at depth d from the
        // outer wall advances spacing (R - d) / (R - c) per station, c the
        // section's half depth, so the tube's length is measured along the
        // centres of its sections and the inner side is shorter.
        bool solid{false};
        bool cross_paths{false};
        double axial_share{1.0};
        double cross_timing{1.0};
        double edge_tension{0.0};
        Vec3 wall{};
        Vec3 inward{};
        double half_width{};
        double fixed_width{}; // m, 0: from the area as a circle
        double bend_radius{0.0}; // m
        bool parallel{false};
    };
    // A radiating opening: an aperture piece with its radiation port.
    struct Opening final {
        std::size_t piece{};
        std::size_t station{};
        double area{};      // m^2, the current geometry
        double admittance{};
        double b0{};
        double b1{};
        double a1{};
        double previous_into{};
        double previous_out{};
        double history{};
        double flow{};
        double previous_flow{};
    };

    void add_links(std::size_t first_piece, std::size_t first_count, std::size_t second_piece,
                   std::size_t second_count, double advance, double share_root, std::size_t lining_station,
                   double& longest);
    void write_wall_filter(std::size_t piece) noexcept;
    void write_radiation(Opening& opening, double area) noexcept;
    [[nodiscard]] std::array<double, 3> boundary_loss(std::size_t link, double crossing) const noexcept;
    [[nodiscard]] std::size_t exit_boundary(std::size_t station) const noexcept;
    void advance_glide() noexcept;
    void refresh_admittances() noexcept;
    [[nodiscard]] double crossing_area(std::size_t link, std::span<const double> station_area) const noexcept;

    AirConstants air_{};
    double sample_rate_{};
    std::vector<Station> stations_{};
    std::vector<std::size_t> tube_first_station_{};
    std::vector<SurfacePiece> pieces_{};
    std::vector<AirLink> links_{};
    std::vector<std::size_t> incoming_start_{};
    std::vector<std::size_t> incoming_links_{};
    std::vector<std::size_t> outgoing_start_{};
    std::vector<std::size_t> outgoing_links_{};
    std::vector<double> admittance_sum_{};      // per piece, including wall and radiation ports
    std::vector<double> arriving_{};            // per link, this sample
    std::vector<double> allpass_input_{};       // per link, previous integer-delayed sample
    std::vector<double> allpass_output_{};      // per link, previous arriving value
    std::vector<double> weighted_{};            // per piece, sum of admittance * arriving
    std::vector<double> pressure_{};            // per piece
    std::vector<double> history_{};             // one ring buffer per link
    std::size_t ring_size_{};
    std::size_t ring_mask_{};
    std::size_t write_index_{};
    std::size_t glottis_{};

    // Boundary-layer loss along each path (viscous friction and heat exchange
    // at the wall, rising as sqrt(f)): a first-order shelf whose loss is
    // fitted to sqrt(f) at and above 3 kHz, passive at every frequency.
    std::vector<double> link_length_{};                // m, axial advance
    std::vector<std::size_t> link_station_{};          // station whose wall lines the path
    std::vector<std::array<double, 3>> loss_filter_{}; // b0 b1 a1
    std::vector<std::array<double, 2>> loss_state_{};  // x1 y1

    // Wall ports: each yielding piece is a mass-spring-damper seen through a
    // port of reference admittance (wetted area)/(rho c); its return wave is
    // the bilinear-discretized reflection R(s) = (z(s) - 1)/(z(s) + 1),
    // z(s) = (m s + r + k/s)/(rho c). The tissue answers in the same
    // instant: the piece pressure and the wall's return are solved together.
    std::vector<double> wall_admittance_{};
    std::vector<std::array<double, 5>> wall_filter_{}; // b0 b1 b2 a1 a2
    std::vector<std::array<double, 4>> wall_state_{};  // x1 x2 y1 y2
    std::vector<double> wall_history_{};

    std::vector<Opening> openings_{};
    std::vector<std::size_t> opening_of_piece_{}; // per piece; no_parent if none

    // Solid tubes. Each wall piece sits at angle phi round its ring's
    // ellipse (0 at the outer wall).
    void add_paired_links(std::size_t first_piece, std::size_t second_piece, std::size_t count, double advance,
                          double share_root, std::size_t lining_station, double& longest);
    static void place_station(const TubeGeometry& shape, double along, Station& station) noexcept;
    void solid_targets(std::span<const double> station_area, std::size_t glide) noexcept;
    [[nodiscard]] double half_width_of(std::size_t station, double area) const noexcept;
    [[nodiscard]] double path_factor(const Station& station, double depth, double area, std::size_t flat) const noexcept;
    double piece_place(std::size_t piece, double area, Vec3& position, Vec3& inward) const noexcept;
    [[nodiscard]] double piece_arc(std::size_t piece, double area) const noexcept;
    std::vector<double> piece_cos_{};      // per piece, cos phi
    std::vector<double> piece_sin_{};      // per piece, sin phi
    std::vector<double> piece_wetted_{};   // per piece, m^2 of surface (solid tubes), current
    // Axial paths of solid tubes have a travel time that follows the depth
    // of their pieces (the inner side of the bend is shorter): their
    // all-pass glides within [whole + 1.05, whole + 3].
    std::vector<double> allpass_step_{};   // per link
    std::vector<double> axial_scale_{};    // per link: sqrt(axial share), 1 for plain tubes
    // Cross paths (links from first_cross_ on): their travel time follows
    // the chord between their pieces and is read by linear interpolation;
    // delay and admittance glide over each block.
    std::size_t first_cross_{};
    // Couplings of parallel tubes (links from first_coupling_ on, inside the
    // interpolated range): the opening between a station and the one beside
    // it, an air path of admittance A_open / (rho c) shared over the parent
    // ring's pieces, travel time from centre to centre.
    std::size_t first_coupling_{};      // first normalized path of a parallel tube (along it, then openings)
    std::size_t first_opening_link_{};  // first opening between a parallel station and its parent
    std::vector<double> opening_target_{};       // per flat station, m^2
    // Openings glide linearly per block and then pass the same two-pole
    // smoothing as the station areas.
    std::vector<double> opening_now_{};
    std::vector<double> opening_step_{};
    std::vector<double> opening_first_{};
    std::vector<double> opening_area_{};
    std::vector<std::size_t> beside_{};          // per flat station: parent station beside it, or no_parent
    std::vector<double> cross_delay_{};          // per cross link, samples (>= 1)
    std::vector<double> cross_delay_step_{};
    std::vector<double> cross_admittance_{};     // per cross link: chords m^3/(Pa s) per m^2 of ring area; others unused
    std::vector<double> cross_admittance_step_{};
    std::vector<double> root_admittance_{};      // per cross link, sqrt(admittance)
    std::vector<double> ring_weight_{};          // scratch per station
    // Edges: each wall piece of a solid tube with tension is joined to its
    // ring neighbours and its axial neighbours. The pull of the edges is a
    // pressure on the wall, e = T / S sum g (x_j - x_i), from the
    // displacements of the previous sample.
    std::vector<std::size_t> edge_start_{};   // per piece + 1
    std::vector<std::size_t> edge_piece_{};
    std::vector<double> edge_gain_{};         // T g / S: Pa per m of displacement difference
    std::vector<double> displacement_{};      // per piece, m outward
    std::vector<double> edge_pressure_{};     // per piece, Pa, this sample
    std::vector<double> injected_{};          // per piece, m^3/s for the next scatter

    // Gliding reshape: loss filters glide linearly; station areas glide
    // linearly and then pass a two-pole smoothing, and every admittance is
    // recomputed from the smoothed areas each sample.
    std::vector<std::array<double, 3>> loss_step_{}; // per link
    std::vector<double> geometry_first_{};           // first smoothing stage, per station
    std::vector<double> geometry_area_{};            // smoothed areas the network uses
    double motion_smoothing_{};
    double contact_smoothing_{};
    std::size_t glide_left_{};

    // Kinetic losses. Boundary b < stations is the passage into station b
    // from its upstream station; b = stations + k is the passage out
    // through opening k. Each path into a downstream piece across a boundary
    // carries the boundary's resistance R Y (from the previous sample's flow).
    std::vector<std::size_t> link_boundary_{};  // per link; no_parent if none
    std::vector<double> boundary_resistance_{}; // R Y per path
    std::vector<double> boundary_flow_{};       // m^3/s downstream, last sample
    // The kinetic resistance follows the boundary's flow through a 1 kHz
    // one-pole: updated from the raw previous-sample flow it chatters at the
    // sample rate where a nearly closed section holds trapped waves.
    std::vector<double> boundary_smoothed_{};
    double smoothing_{};
    std::vector<std::size_t> boundary_upstream_{};   // station before the boundary
    std::vector<std::size_t> boundary_downstream_{}; // station after it; no_parent: the free air
    // Station areas glide with the path weights, so the losses always see
    // the cross-sections the paths currently carry.
    std::vector<double> station_area_now_{};
    std::vector<double> station_area_step_{};
    std::vector<double> boundary_admittance_{}; // scratch: path admittance across each boundary
    std::vector<std::size_t> piece_boundary_{}; // per piece: the boundary its resisted paths cross
    std::vector<double> resisted_sum_{};        // per piece, admittance of its resisted paths
    std::size_t turbulence_boundary_{no_parent};
    std::size_t fading_boundary_{no_parent};  // the previous place, fading out
    double turbulence_weight_{1.0};           // of the current place
    double turbulence_fade_{};                // weight gained per sample
    double wall_flow_{};
};

}
