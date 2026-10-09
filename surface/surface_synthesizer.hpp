#pragma once
#include "surface_network.hpp"
#include "ducklips/air_path.hpp"
#include "ducklips/articulation.hpp"
#include "ducklips/cues.hpp"
#include "ducklips/sources.hpp"
#include "ducklips/voice.hpp"
#include "ducklips/vowel_library.hpp"
#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {

struct SurfaceSynthesisOptions final {
    SynthesisRate rate{};                  // output rate
    std::size_t oversample{4};             // the surface runs at rate * oversample
    std::size_t oral_sections{64};
    std::size_t nasal_sections{24};
    // Pieces per airway ring. With ring-to-ring paths only, every piece of
    // a ring carried the same plane wave (1 and 3 rendered identically,
    // MEASURED -188 dB). The solid airway's cross paths make them differ:
    // they carry the cross modes (a 3.57 cm circle's first mode within 1%
    // at 8 pieces; 12 pieces changed nothing, MEASURED).
    std::size_t ring_pieces{8};
    std::size_t branch_ring_pieces{1};     // pieces per ring of the nose and side cavities
    // The airway as a solid surface: pieces on bent elliptical
    // cross-sections with air paths across every ring (TubeGeometry::solid,
    // cross_paths). Needs ring_pieces of 6 or more.
    bool solid_airway{true};
    double axial_share{0.6};
    double cross_timing{0.86};
    double edge_tension{0.0};              // N/m
    // A floor under the tongue: the space below the lifted blade and tip
    // becomes a channel beside the mouth, open to it in front of the tip
    // (TractModel::set_floor, TubeGeometry::parallel).
    bool floor{true};
    double floor_opening_scale{1.0};       // diagnostics: scales the openings (0 isolates the space)
    double floor_minimum_area{0.0};        // diagnostics: m^2 floor on the space's sections
    std::size_t control_block{48};         // output samples between articulator updates
    std::size_t glottis_substeps{1};
    std::size_t cue_capacity{4096};
    // Fold surface (see FoldParameters): the narrowest point slides along a
    // rounded medial surface and the closing glottis seals gradually.
    // At 1 and 0.1 mm^2 the closure was so gradual that the source's rate
    // of change fell about 14 dB per octave above 1 kHz and every vowel was
    // dark above 3 kHz; at 0.25 and 0.03 mm^2 it falls about 6 dB per
    // octave, as a glottal pulse with an abrupt closure does (MEASURED on
    // a sustained /a/).
    double medial_roundness{0.25e-6};      // m^2
    double contact_softness{0.03e-6};      // m^2
    // Glottal turbulence: the voice's aspiration gain times this scale,
    // shaped by two poles at the corner (MEASURED against a human recording
    // for the baritone: 0.012 of jet dynamic pressure, 2 kHz).
    double aspiration_scale{4.8};
    double aspiration_corner{2000.0};      // Hz
    double whisper_aspiration_corner{4500.0}; // Hz, when the folds do not phonate
    // Frication: the voice's frication gain times this scale (MEASURED: a
    // human recording's sibilants sit 4.6 dB under its vowels).
    double frication_scale{40.0};
    // The jet's noise falls above about Strouhal 0.2 of its speed over its
    // width (HYPOTHESIS); the corner follows the jet each block, bounded here.
    double frication_strouhal{0.2};
    double lowest_frication_corner{1000.0};  // Hz
    double highest_frication_corner{12000.0}; // Hz
    // A jet that strikes the teeth is a stronger source than a free jet
    // (Shadle's obstacle experiments: 10-20 dB). The source then sits at the
    // teeth, so only the small space beyond them colours it.
    double burst_scale{0.35};              // a stop's release noise, relative to a fricative's
    double obstacle_boost{1.6};            // times the voice's obstacle gain
    // A reversed glottal jet (inhaling) points into the trachea; less of its
    // turbulence reaches the mouth.
    double inhaled_aspiration{0.3};
    // A jet that leaves the lips (f, v, a lip release) is turbulent in the
    // free air beyond them: its fluctuating force radiates directly, not
    // through the nearly closed opening. Modelled as an extra radiated flow
    // q A / (rho c) at the lips, times this (HYPOTHESIS).
    double free_jet_radiation{3.0};
    // Cycle-scale unsteadiness, at the voice's jitter of 0.006.
    double tremor{0.002};                  // fraction of tension, 40 Hz band
    // The voice's slow tension wander (12 Hz band) is scaled by this: at
    // the preset 0.6% rms it was a 25-40 Hz flutter of about 1.5% in F0,
    // heard as gargling in the child (MEASURED cycle by cycle).
    double slow_wander_scale{0.4};
    // Laryngeal pitch control: integral correction per second of the log
    // error between intended and measured cycle rate.
    double pitch_correction_rate{30.0};
    double adduction_wander{0.08};         // fraction of the rest opening, 20 Hz band
    // Diagnostics: which radiators reach the microphones (bit 0 lips,
    // bit 1 nostrils, bit 2 throat).
    unsigned radiator_mask{7u};
    // Diagnostics: per-internal-sample trace from trace_start (s) on.
    double trace_start{-1.0};
    std::size_t trace_samples{0};
    AirPathOptions air_path{};
    AirProperties air{};
};

struct SurfaceTraceSample final {
    double lip_area{};       // m^2, gliding
    double lip_resistance{}; // R Y of the lip opening's kinetic loss
    double lip_flow{};       // m^3/s
    double glottal_flow{};   // m^3/s
    double frication{};      // Pa
    double aspiration{};     // Pa
    double entry_last{};     // R Y into the last airway station
    double entry_before{};   // R Y into the station before it
    double lip_pressure{};   // Pa, last airway station
};

struct SurfaceProbe final {
    double glottal_flow{};       // m^3/s
    double glottal_area{};       // m^2
    double intraoral_pressure{}; // Pa, behind the narrowest section
    double lung_pressure{};      // Pa
    double intended_pitch{};     // Hz
    double measured_pitch{};     // Hz
    double narrowest_area{};     // m^2
    double jet_flow{};           // m^3/s
    std::size_t jet_section{};   // narrowest oral section
    double lip_area{};           // m^2
    double port_area{};          // m^2
};

// The DuckLips body with the surface engine as its acoustics.
//
// The behavioural layer is the modal synthesizer's: cues set targets, the
// articulators move with their own speeds and ramps, the tongue model turns
// postures into cross-sections. Every control block the cross-sections are
// handed to the surface network, which glides to them sample by sample. The
// nose opens from the velar port and the side cavities hang off the airway
// and the nose as more surface. The folds are driven by lung pressure minus
// the surface pressure above them; the jet leaves the narrowest section with
// its kinetic loss and turbulence. Lips, nostrils and the yielding walls
// (through the neck and cheeks) radiate to two microphones.
//
// prepare() allocates. render() does not allocate, lock or throw.
class SurfaceSynthesizer final {
public:
    void prepare(const VoiceProfile& voice, const SurfaceSynthesisOptions& options);
    void reset() noexcept;
    [[nodiscard]] bool enqueue(const Cue& cue) noexcept;
    [[nodiscard]] std::size_t queue_space() const noexcept;
    void render(std::span<float> stereo) noexcept;

    [[nodiscard]] bool finished() const noexcept;
    [[nodiscard]] SampleCount position() const noexcept;
    [[nodiscard]] const VowelLibrary& vowels() const noexcept;
    [[nodiscard]] const SurfaceProbe& probe() const noexcept;
    // Read-only offline diagnostics; callers may copy a prepared network for
    // independent frozen-geometry impulse experiments.
    [[nodiscard]] const SurfaceNetwork& acoustic_network() const noexcept { return network_; }
    [[nodiscard]] double unit_tension_pitch() const noexcept;
    [[nodiscard]] std::size_t piece_count() const noexcept;
    [[nodiscard]] std::size_t link_count() const noexcept;
    [[nodiscard]] const std::vector<SurfaceTraceSample>& trace() const noexcept;
    // Every station's target area (m^2) from the last control block: the
    // airway's oral sections first (glottis to lips), then the nose and the
    // side cavities.
    [[nodiscard]] std::span<const double> station_areas() const noexcept;

private:
    double calibrate_folds() const noexcept;
    void apply_cue(const Cue& cue) noexcept;
    void begin_block() noexcept;
    void write_station_areas() noexcept;
    void advance_internal(double fraction) noexcept;
    void track_pitch(double glottal_area) noexcept;
    [[nodiscard]] double decimated(std::size_t radiator) const noexcept;

    VoiceProfile voice_{};
    SurfaceSynthesisOptions options_{};
    AirProperties air_{};
    double internal_rate_{};
    double internal_period_{};
    VowelLibrary vowels_{};
    TractModel model_{};
    TractShape shape_{};
    SurfaceNetwork network_{};
    std::vector<double> station_area_{};
    std::size_t nose_tube_{};
    std::vector<std::size_t> side_tube_{};
    std::size_t floor_first_{};            // first oral section with space under the tongue
    std::size_t floor_count_{};            // 0 without a floor
    std::vector<double> opening_{};        // per flat station, m^2
    TwoMassGlottis glottis_{};
    CueQueue queue_{};
    ArticulatorState articulators_{};
    AirPath air_path_{};

    // Sources.
    NoiseGenerator buffet_noise_{};
    TurbulenceNoise aspiration_noise_{};
    TurbulenceNoise frication_noise_{};
    std::array<double, 3> aspiration_shaper_{};
    std::array<double, 3> frication_shaper_{};
    double aspiration_coefficient_{};
    double frication_coefficient_{};
    SlowNoise jitter_noise_{};
    SlowNoise shimmer_noise_{};
    SlowNoise tremor_noise_{};
    SlowNoise adduction_noise_{};
    SlowNoise roughness_noise_{};
    SlowNoise roughness_wander_{};

    // Decimation of the three radiators (lips, nostrils, throat) to the
    // output rate: a Blackman-windowed sinc evaluated at output instants.
    std::vector<double> kernel_{};
    std::vector<double> reversed_kernel_{};
    std::vector<double> radiator_history_{}; // [radiator * 2 * history_size + index], stored twice
    std::size_t history_size_{};
    std::size_t history_index_{};

    // Timing.
    SampleCount position_{};
    SampleCount finish_time_{-1};
    std::size_t block_remaining_{};
    std::size_t block_internal_{};
    std::size_t internal_in_block_{};
    double sample_period_{};

    // Block targets and their previous values (per-sample interpolation).
    double lung_pressure_{};
    double previous_lung_pressure_{};
    double rest_area_{};
    double previous_rest_area_{};
    double tension_{1.0};
    double previous_tension_{1.0};
    double output_gain_{1.0};
    double previous_output_gain_{1.0};
    double jet_area_{};
    double previous_jet_area_{};
    std::size_t jet_station_{};
    bool jet_at_lips_{};
    double free_jet_flow_{};
    double exit_flow_smoothed_{}; // smoothed pressure drop across the jet (Pa)
    double exit_speed_smoothed_{};  // smoothed flow over area through the jet (m/s)
    std::size_t jet_behind_{};
    std::size_t jet_beyond_{no_parent};
    double exit_smoothing_{};
    double frication_gain_{};
    double previous_frication_gain_{};
    double tremor_phase_{};

    // Sample state.
    double transglottal_{};
    double glottal_flow_{};
    double glottal_area_{};
    double previous_lip_flow_{};
    double previous_nose_flow_{};
    double previous_wall_flow_{};

    // Laryngeal pitch control from the measured fold cycle.
    double unit_pitch_{120.0};
    double pitch_integral_{};
    double last_area_{};
    double trough_area_{};
    bool rising_{};
    SampleCount last_peak_{-1};
    SampleCount last_voiced_{};
    SampleCount internal_position_{};
    double measured_pitch_{};

    SurfaceProbe probe_{};
    std::vector<SurfaceTraceSample> trace_{};
    SampleCount trace_begin_{-1};
};

}
