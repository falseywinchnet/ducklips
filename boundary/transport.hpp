#pragma once
// Boundary-surface transport: air is not a resonator. Each patch's pressure
// is the time-domain Kirchhoff-Helmholtz collocation sum over the patches it
// sees, with propagation delay:
//
//   p_j(t) = 2 R_j sum_i A_i [ cos(b) p_i(t-r/c) / (4 pi r^2)
//                            + cos(b) dp_i/dt(t-r/c) / (4 pi c r)
//                            - rho a_i(t-r/c) / (4 pi r) ] + forcing_j(t)
//
// where b is the angle between the source normal (into its wall) and the
// direction from receiver to source, and a_i is the source wall's normal
// acceleration into the wall. The microphone uses the same sum with factor 1.
// Elastic surfaces carry their own modal decomposition (strips joined at
// edges); pressure plucks the modes and modal acceleration re-emits.
#include "geometry.hpp"
#include "ducklips/physics.hpp"
#include <cstddef>
#include <cstdint>
#include <vector>

namespace boundary {

struct TransportSettings final {
    double audio_rate{48000.0};
    std::size_t oversampling{4};
    bool membranes{true};
    double quadrature_resolution{3.0}; // subdivisions per (patch edge / distance), clamped to [2, 8]
    double reflectance_override{-1.0};  // when > 0, replaces every tissue reflectance
    // true: links only between patches that face each other with a clear line
    // of sight (the architect's rule). false: the full Kirchhoff kernel, with
    // signed cos(b) from behind and no occlusion (exact for a closed boundary).
    bool visibility_culling{true};
    // Two-point time average of each new wall pressure with the previous one
    // (a Rynne-Smith style MOT stabilizer; gain cos(pi f / rate), half-sample lag).
    bool temporal_averaging{false};
    ducklips::AirProperties air{};
};

struct CompileStatistics final {
    std::size_t patches{0};
    std::size_t linked_pairs{0};
    std::size_t pressure_taps{0};
    std::size_t acceleration_taps{0};
    std::size_t microphone_taps{0};
    std::size_t clamped_subpoints{0};   // delay below one internal sample, clamped to one
    std::size_t occluded_subpoints{0};
    std::size_t membrane_modes{0};
    double highest_membrane_hz{0.0};
    double largest_static_gain{0.0};    // max over receivers of the summed static (DC) pressure gain
    std::size_t ring_length{0};
    double compile_seconds{0.0};
};

class BoundaryTransport final {
public:
    // Allocates every buffer. The scene is not retained.
    void compile(const Scene& scene, const TransportSettings& settings);
    void reset() noexcept;
    // Replaces the scene's source patches by one wall-acceleration patch
    // (reciprocity probes). Call after compile; rebuilds the active tap lists (allocates).
    void inject_acceleration_at(std::size_t patch);
    // Advances one internal sample; source is the jet force pressure (Pa) or
    // the wall acceleration (m/s^2). Returns the microphone pressure (Pa).
    // Does not allocate.
    double step(double source) noexcept;
    [[nodiscard]] double patch_pressure(std::size_t patch) const noexcept;
    [[nodiscard]] double internal_rate() const noexcept;
    [[nodiscard]] const CompileStatistics& statistics() const noexcept;

private:
    struct Tap final {
        std::uint32_t delay{0};
        std::uint32_t source{0};
        double weight{0.0};
    };

    struct TapOrder final {
        bool operator()(const Tap& a, const Tap& b) const noexcept;
    };

    void compile_membranes(const Scene& scene);
    // Copies the acceleration taps whose source can move (an elastic patch
    // with membranes on, or an injected patch) into the active lists.
    void select_active_acceleration_taps();
    void accumulate_pair(const Scene& scene, const Vec3& receiver, const Vec3& receiver_normal, bool has_normal,
                         std::size_t source, double factor, std::vector<double>& pressure_by_delay,
                         std::vector<double>& acceleration_by_delay);
    void add_interpolated(std::vector<double>& by_delay, double delay, double weight);
    void flush_taps(std::vector<double>& by_delay, std::size_t source, std::vector<Tap>& taps);
    [[nodiscard]] bool segment_blocked(const Scene& scene, const Vec3& from, const Vec3& to) const;

    std::size_t patch_count_{0};
    std::size_t ring_length_{0};
    std::size_t ring_mask_{0};
    std::uint64_t step_index_{0};
    double rate_{0.0};
    double time_step_{0.0};
    ducklips::AirProperties air_{};
    double quadrature_resolution_{3.0};
    bool visibility_culling_{true};
    bool temporal_averaging_{false};
    CompileStatistics statistics_{};

    // Receiver j's pressure taps are [pressure_begin_[j], pressure_begin_[j + 1]).
    std::vector<std::size_t> pressure_begin_{};
    std::vector<Tap> pressure_taps_{};
    std::vector<std::size_t> acceleration_begin_{};
    std::vector<Tap> acceleration_taps_{};
    std::vector<Tap> microphone_pressure_taps_{};
    std::vector<Tap> microphone_acceleration_taps_{};
    // Every compiled acceleration tap; the lists above hold the active subset.
    std::vector<std::size_t> all_acceleration_begin_{};
    std::vector<Tap> all_acceleration_taps_{};
    std::vector<Tap> all_microphone_acceleration_taps_{};
    std::vector<bool> patch_moves_{};      // elastic with membranes on

    // Time-major rings: slot * patch_count_ + patch.
    std::vector<double> pressure_ring_{};
    std::vector<double> acceleration_ring_{};
    std::vector<double> forcing_sign_{};
    std::vector<double> injection_sign_{};
    std::vector<double> area_{};

    // Membrane modes. Mode k covers patches [mode_first_patch_[k], + mode_patch_count_[k]);
    // its mass-normalized shape is mode_shape_[mode_shape_offset_[k] + local].
    std::vector<std::size_t> mode_first_patch_{};
    std::vector<std::size_t> mode_patch_count_{};
    std::vector<std::size_t> mode_shape_offset_{};
    std::vector<double> mode_shape_{};
    std::vector<double> mode_omega_{};
    std::vector<double> mode_damping_{};
    std::vector<double> mode_position_{};
    std::vector<double> mode_velocity_{};
    std::vector<double> membrane_acceleration_{};

    std::size_t largest_delay_{0};
    bool acceleration_active_{false}; // membranes or an injected wall acceleration
};

}
