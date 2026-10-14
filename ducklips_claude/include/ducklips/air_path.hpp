#pragma once
#include "ducklips/physics.hpp"
#include <array>
#include <cstddef>
#include <vector>

namespace ducklips {

struct Position3 final {
    double x{}; // m, forward from the head centre
    double y{}; // m, to the speaker's left
    double z{}; // m, up
};

// Optional rectangular room. Images up to maximum_order reflections are kept
// while their pressure gain stays above cutoff_ratio of the direct path.
struct RoomGeometry final {
    bool enabled{false};
    double length{5.0};          // m, along x
    double width{4.0};           // m, along y
    double height{2.8};          // m, along z
    Position3 head{2.2, 1.9, 1.6};
    double wall_reflection{0.75}; // pressure reflection coefficient per bounce
    std::size_t maximum_order{3};
    double cutoff_ratio{0.02};
};

struct AirPathOptions final {
    Position3 lips{0.09, 0.0, -0.03};
    Position3 nostrils{0.10, 0.0, 0.0};
    Position3 throat{0.05, 0.0, -0.10};
    Position3 left_microphone{0.40, 0.085, 0.0};
    Position3 right_microphone{0.40, -0.085, 0.0};
    RoomGeometry room{};
    double full_scale_pressure{2.0}; // Pa mapped to sample value 1
    double dc_block_hz{30.0};
};

enum class Radiator : std::size_t { lips = 0, nostrils = 1, throat = 2 };
inline constexpr std::size_t radiator_count = 3;
inline constexpr std::size_t microphone_count = 2;

// Propagation from the three radiating openings to two microphones.
//
// The structure is the photonic transport field's: emitters, sparse links
// with per-link extinction, receivers, and unit responses compiled once for
// fixed geometry and scaled by the live source. Unlike light energy, sound at
// a microphone is coherent, so each link carries a delay and a signed pressure
// gain instead of an energy fraction. A radiator's volume-flow derivative
// dU/dt becomes far-field pressure rho/(4 pi d) dU/dt(t - d/c) per link
// (acoustic monopole). Reflected links pass a one-pole low-pass per link for
// frequency-dependent wall loss.
class AirPath final {
public:
    // Allocates and compiles the links.
    void prepare(const AirPathOptions& options, const AirProperties& air, SynthesisRate rate);
    void reset() noexcept;
    // One sample of radiator flow derivatives (m^3/s^2) in, microphone samples out.
    void process(const std::array<double, radiator_count>& flow_derivative,
                 std::array<double, microphone_count>& microphone) noexcept;
    [[nodiscard]] std::size_t link_count() const noexcept;

private:
    struct Link final {
        std::size_t radiator{};
        std::size_t microphone{};
        double delay_samples{};
        double gain{};
        double smoothing{1.0}; // one-pole coefficient; 1 passes unfiltered
        double filter_state{};
    };
    void add_links(const Position3& source, std::size_t radiator, const Position3& receiver, std::size_t microphone,
                   const AirPathOptions& options, const AirProperties& air, SynthesisRate rate);
    [[nodiscard]] double read_delayed(std::size_t radiator, double delay_samples) const noexcept;

    std::vector<Link> links_{};
    std::vector<double> history_{}; // [radiator * ring_size + index]
    std::size_t ring_size_{};
    std::size_t ring_mask_{};
    std::size_t write_index_{};
    double output_scale_{};
    double dc_coefficient_{};
    std::array<double, microphone_count> dc_input_{};
    std::array<double, microphone_count> dc_output_{};
};

}
