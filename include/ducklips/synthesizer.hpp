#pragma once
#include "ducklips/acoustic_network.hpp"
#include "ducklips/air_path.hpp"
#include "ducklips/articulation.hpp"
#include "ducklips/cues.hpp"
#include "ducklips/modal.hpp"
#include "ducklips/canonical_speech.hpp"
#include "ducklips/local_speech.hpp"
#include "ducklips/mean_airway.hpp"
#include "ducklips/sources.hpp"
#include "ducklips/tree_operator.hpp"
#include "ducklips/voice.hpp"
#include "ducklips/vowel_library.hpp"
#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {

// Where the output is taken from (diagnostics): the microphones, the lip
// radiation before the air path, or the glottal flow derivative alone.
enum class OutputTap { microphones, lips, glottal_source };

struct SynthesisOptions final {
    bool local_flow{false}; // bounded local exchanges; unified mean/acoustic air state
    std::size_t local_substeps{4};
    bool conservative_tongue{false}; // material body map; explicit research option
    bool canonical_motion{false}; // experimental full physical moving-state path
    bool spatial_pressure{false}; // experimental local mean-pressure/flow coupling; requires canonical motion
    std::size_t canonical_substeps{8};
    SynthesisRate rate{};
    std::size_t oral_sections{64};
    std::size_t nasal_sections{24};
    std::size_t mode_capacity{40};
    double maximum_mode_frequency{10000.0}; // Hz; modes above are not simulated
    std::size_t control_block{48};          // samples between articulator updates
    double geometry_step_limit{0.15};       // largest log-area change per geometry step
    std::size_t maximum_geometry_steps{8};  // geometry steps per control block
    std::size_t glottis_substeps{2};
    std::size_t cue_capacity{4096};
    ReprojectionInvariant invariant{ReprojectionInvariant::passive};
    AirPathOptions air_path{};
    AirProperties air{};
    bool diagnostics{false};                // measure basis orthogonality every block (costly)
    OutputTap tap{OutputTap::microphones};
    bool radiate_mean_flow{true};           // radiate d/dt of the aerodynamic outflow (release pop)
    double trace_start{-1.0};               // s; >= 0 records per-sample radiator values (diagnostics)
    std::size_t trace_samples{0};
};

struct SynthesisStatistics final {
    double geometry_seconds{}; // measured control-rate body construction, separate from audio updates
    std::size_t airflow_failures{};
    double airflow_failure_time{};
    double airflow_failure_residual{};
    std::size_t airflow_failure_iterations{};
    std::size_t geometry_failures{};
    const char* geometry_error{nullptr};
    double geometry_failure_time{};
    double maximum_body_wall_penetration{}; // cm; contact validity is separate from solver convergence
    double maximum_tip_target_error{}; // cm, unattainable tip target under body contact and blade reach
    double maximum_attachment_motion{}; // cm, local material displacement of blade attachment
    std::size_t nonfinite_samples{};
    std::size_t clipped_samples{};
    std::size_t blocks{};
    std::size_t warm_modes{};
    std::size_t bisected_modes{};
    std::size_t factorizations{};
    std::size_t fewest_modes{static_cast<std::size_t>(-1)};
    std::size_t most_modes{};
    std::size_t pitch_cycles{};
    std::size_t geometry_updates{};
    std::size_t most_geometry_steps{};
    double maximum_cavity_mass_residual{}; // m^3/s, moving mean-cavity balance
    double worst_orthogonality{};           // max |phi_i^T C phi_j - delta_ij| when diagnostics are on
    double worst_orthogonality_time{};      // s
};

// One sample of radiated flow derivatives (m^3/s^2) and the acoustic pressure
// above the glottis, for diagnostics.
struct TraceSample final {
    double lips{};
    double nostrils{};
    double throat{};
    double glottal_pressure{};
    double energy{};
    std::size_t modes{};
    double glottal_area{};  // m^2
    double glottal_flow{};  // m^3/s
    double tension{};
    // Geometry-event impulses, zero on samples without an event. Lip values
    // are radiation-flow derivatives (m^3/s^2), pressure values are Pa.
    double geometry_lip_metric_jump{};
    double geometry_lip_projection_jump{};
    double geometry_pressure_jump{};
    double geometry_lip_slope_jump{};
    double geometry_energy_ratio{};
    double geometry_parameter_change{};
    double frication_pressure{};
    std::size_t frication_branch{};
    double constriction_area{};
    double cavity_pressure{};
    double cavity_displacement_flow{};
    double cavity_mass_residual{};
    double oral_mean_flow{},nasal_mean_flow{},constriction_flow{};
};

// Values from the most recent control block, for diagnostics.
struct SynthesisProbe final {
    double glottal_flow{};       // m^3/s, mean over the block
    double glottal_area{};       // m^2, mean over the block
    double intraoral_pressure{}; // Pa
    double lung_pressure{};      // Pa
    double intended_pitch{};     // Hz
    double measured_pitch{};     // Hz, from the last fold cycle
    double tension{};            // Q
    double narrowest_area{};     // m^2
    double acoustic_energy{};    // J-scaled modal energy
    double transglottal{};       // Pa, last sample
    double glottal_pressure{};   // Pa, acoustic pressure above the glottis, last sample
    std::size_t modes{};
    std::array<double, 4> lowest_mode_hz{};
    // Actual geometry used by the acoustic solve; populated at every block.
    Point2 blade_root{};
    double attachment_motion{},tip_target_error{},body_wall_penetration{};
    // Geometry attribution, populated only with diagnostics enabled.
    double oral_volume{},floor_volume{},volume_without_body{},volume_without_tip{};
    Point2 tongue_centre{},tongue_tip{};
    double body_activation{},tip_activation{};
};

// The DuckLips engine: a simulated body that is always running.
//
// Every control block (1 ms by default) the articulators move toward their
// cue targets, the tract geometry is rebuilt, its acoustic eigenstates are
// updated by warm-started sparse solves, and the acoustic state is carried
// into the new eigenstates. Every sample the vocal folds, the airflow, the
// turbulence and the modes advance together, and the radiated flow is carried
// through the air to two microphones.
//
// prepare() allocates. enqueue() belongs to one control thread. render()
// belongs to one audio thread: it does not allocate, lock or throw, and fills
// interleaved L,R floats, which is the gui_forms::AudioGenerator contract.
class Synthesizer final {
public:
    void prepare(const VoiceProfile& voice, const SynthesisOptions& options);
    void reset() noexcept;
    [[nodiscard]] bool enqueue(const Cue& cue) noexcept;
    [[nodiscard]] std::size_t queue_space() const noexcept;
    void render(std::span<float> stereo) noexcept;

    [[nodiscard]] bool finished() const noexcept;
    [[nodiscard]] SampleCount position() const noexcept;
    [[nodiscard]] const SynthesisStatistics& statistics() const noexcept;
    [[nodiscard]] const SynthesisProbe& probe() const noexcept;
    [[nodiscard]] const VowelLibrary& vowels() const noexcept;
    [[nodiscard]] const VoiceProfile& voice() const noexcept;
    [[nodiscard]] const SynthesisOptions& options() const noexcept;
    [[nodiscard]] double unit_tension_pitch() const noexcept; // Hz at Q = 1, measured at prepare
    [[nodiscard]] const std::vector<TraceSample>& trace() const noexcept;

private:
    double calibrate_folds() const noexcept;
    void apply_cue(const Cue& cue) noexcept;
    void begin_block() noexcept;
    void begin_geometry_step() noexcept;
    void render_samples(std::span<float> stereo) noexcept;
    void render_local_samples(std::span<float> stereo) noexcept;
    void track_pitch(double glottal_area) noexcept;

    VoiceProfile voice_{};
    SynthesisOptions options_{};
    AirProperties air_{};
    VowelLibrary vowels_{};
    TractModel model_{};
    TractModel diagnostic_model_{};
    TractShape diagnostic_shape_{};
    TractShape shape_{};          // target geometry of the current block
    TractShape start_shape_{};    // geometry at the start of the block
    TractShape geometry_shape_{}; // geometry currently simulated
    TractShape local_tangent_shape_{};
    AcousticNetwork local_tangent_network_{};
    std::array<AcousticNetwork, 2> networks_{};
    std::array<ModalBasis, 2> bases_{};
    std::size_t current_{};
    TreeOperator matrix_{};
    ModalSolver solver_{};
    ModalResonator resonator_{};
    CanonicalSpeech canonical_{};
    LocalSpeech local_{};
    bool canonical_ready_{false};
    TwoMassGlottis glottis_{};
    IntraoralPressure intraoral_{};
    MeanAirway mean_airway_{};
    TurbulenceNoise aspiration_noise_{};
    TurbulenceNoise frication_noise_{};
    SlowNoise jitter_noise_{};
    SlowNoise shimmer_noise_{};
    AirPath air_path_{};
    CueQueue queue_{};
    ArticulatorState articulators_{};

    // Timing.
    SampleCount position_{};
    SampleCount finish_time_{-1};
    std::size_t block_remaining_{};
    std::size_t geometry_remaining_{};
    std::size_t geometry_steps_{1};
    std::size_t geometry_step_index_{};
    bool geometry_ready_{};
    std::size_t aspiration_branch_{};
    std::size_t frication_branch_{};
    double sample_period_{};

    // Block constants.
    double lung_pressure_{};
    double shimmer_multiplier_{1.0};
    double pitch_multiplier_{1.0};
    double pitch_correction_{1.0};
    double rest_area_{};
    double tension_{1.0};
    double oral_outlet_{};
    double nasal_outlet_{};
    double previous_oral_outlet_{};
    double previous_nasal_outlet_{};
    double cavity_compliance_{};
    double previous_cavity_compliance_{};
    double cavity_volume_{};
    double previous_cavity_volume_{};
    double constriction_area_{};
    double frication_scale_{};
    double tremor_phase_{};
    double output_gain_{1.0};

    // Sample state.
    double transglottal_{};
    double glottal_pressure_{};
    double glottal_flow_{};
    double tap_previous_flow_{};
    double glottal_area_{};
    double mean_coefficient_{};
    double mean_stage_{};
    double mean_flow_{};
    double previous_oral_outflow_{};
    double previous_nasal_outflow_{};
    double oral_stage_{};
    double nasal_stage_{};

    // Laryngeal pitch control from the measured fold cycle.
    double unit_pitch_{120.0};
    double pitch_integral_{};
    double last_area_{};
    double trough_area_{};
    bool rising_{};
    SampleCount last_peak_{-1};
    double measured_pitch_{};

    // Probe accumulation.
    double flow_sum_{};
    double area_sum_{};
    std::size_t probe_samples_{};

    SynthesisStatistics statistics_{};
    SynthesisProbe probe_{};
    std::vector<TraceSample> trace_{};
    double geometry_lip_metric_jump_{};
    double geometry_lip_projection_jump_{};
    double geometry_pressure_jump_{};
    double geometry_lip_slope_jump_{};
    double geometry_energy_ratio_{};
    double geometry_parameter_change_{};
    SampleCount trace_begin_{-1};
};

}
