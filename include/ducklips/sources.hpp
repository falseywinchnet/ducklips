#pragma once
#include "ducklips/physics.hpp"
#include <array>
#include <cstddef>
#include <cstdint>

namespace ducklips {

// Two-mass vocal fold parameters (Ishizaka and Flanagan 1972, in the
// simplified form of Steinecke and Herzel 1995), converted to SI. The tension
// factor Q scales stiffness by Q and mass by 1/Q, so the natural frequencies
// scale by Q while damping ratios are unchanged.
struct FoldParameters final {
    double lower_mass{1.25e-4};          // kg
    double upper_mass{2.5e-5};           // kg
    double lower_stiffness{80.0};        // N/m
    double upper_stiffness{8.0};         // N/m
    double coupling_stiffness{25.0};     // N/m
    double collision_ratio{3.0};         // contact stiffness / mass stiffness
    double lower_damping_ratio{0.1};
    double upper_damping_ratio{0.6};
    double closed_damping_increase{1.0}; // added damping ratio while a mass is in contact
    double length{0.014};                // m, glottal length
    double lower_thickness{0.0025};      // m
    double upper_thickness{0.0005};      // m
    // Strips along the glottal length (at most maximum_fold_strips). With one
    // strip the folds are a rigid rectangle that closes along its whole length
    // at once. With several, each strip moves by the half-sine shape of the
    // fold vibration (the anterior commissure and the vocal process barely
    // move) and the rest opening is a posterior chink, so the glottis closes
    // progressively and a small leak stays open.
    std::size_t strips{1};
    // Flow area at which the closing glottis seals gradually rather than at
    // a corner (m^2; 0: hard corner). HYPOTHESIS: the curved, mucus-covered
    // medial surface seals over a fraction of a millisecond, which is what
    // ends a real voice's harmonics above about 4 kHz.
    double contact_softness{0.03e-6};
    // Shared by both engines: the modal path formerly left these at zero,
    // reintroducing the hard-minimum corner already diagnosed in the surface
    // experiments. These are explicit model parameters, not measured anatomy.
    // Area scale over which the narrowest point moves between the lower and
    // upper mass (m^2; used with contact_softness > 0). Real medial surfaces
    // are rounded, so the flow area has no corner at the exchange.
    double medial_roundness{0.25e-6};
    // Tension of the right fold relative to the left (1: symmetric). An
    // asymmetric larynx vibrates with subharmonics and irregular cycles
    // (Steinecke and Herzel 1995): the physical origin of a rough, gravelly
    // voice.
    double asymmetry{1.0};
    // Irregular fold tissue (thickened, mucus-coated): rms fractional
    // fluctuation of stiffness and adduction at the scale of single cycles
    // (150 Hz band). Measured rough and hoarse voices have several percent
    // jitter and shimmer; this is where a gravelly voice gets them. Applied
    // by the synthesizer.
    double roughness{0.0};
};

inline constexpr std::size_t maximum_fold_strips = 64;

struct FoldState final {
    std::array<double,2> lower_position{},lower_velocity{},upper_position{},upper_velocity{};
    double rest_area{};
};

// Self-oscillating vocal folds. Positions are the lateral displacement of
// each mass from its rest position; the glottal area at a mass is
// rest_area + 2 * length * displacement. Flow separates at the narrowest
// point, so only the lower mass is driven by the Bernoulli pressure.
class TwoMassGlottis final {
public:
    void reset(const FoldParameters& folds) noexcept;
    // Advances one audio sample in substeps of symplectic Euler. Pressure is
    // the transglottal drop (Pa); rest_area is the abduction (m^2); tension is
    // Q. surface_turbulence (Pa) is the fluctuating jet pressure on the lower
    // fold surface. Without it the symmetric rest posture is an exact, if
    // unstable, equilibrium and the folds never start.
    void advance(double transglottal_pressure, double rest_area, double tension, double sample_period,
                 std::size_t substeps, double surface_turbulence) noexcept;
    [[nodiscard]] double lower_area() const noexcept;
    [[nodiscard]] double upper_area() const noexcept;
    [[nodiscard]] double minimum_area() const noexcept;
    [[nodiscard]] FoldState state() const noexcept;
    // Restore physical state under the already configured tissue parameters.
    // Nonfinite input is rejected without changing the folds.
    [[nodiscard]] bool restore(const FoldState& state) noexcept;

private:
    // position_sum: the two folds' displacements added (left + right).
    double strip_gap(double rest_area, std::size_t strip, double position_sum) const noexcept;

    FoldParameters folds_{};
    std::array<double, maximum_fold_strips> shape_{};  // displacement per unit mass displacement, mean 1
    std::array<double, maximum_fold_strips> opening_{}; // share of the rest area, mean 1
    std::size_t strips_{1};
    // Per fold (0 left, 1 right): lower and upper mass positions and velocities.
    std::array<double, 2> lower_position_{};
    std::array<double, 2> lower_velocity_{};
    std::array<double, 2> upper_position_{};
    std::array<double, 2> upper_velocity_{};
    double rest_area_{};
};

// Quasi-steady Bernoulli flow through an orifice of area a with an extra
// series impedance Z between the free pressure drop and the orifice:
//   U = a sqrt(2 (drop - Z U) / rho), sign following the drop.
// Solved in closed form from the quadratic in U.
[[nodiscard]] double orifice_flow(double area, double pressure_drop, double series_impedance,
                                  const AirProperties& air) noexcept;

// Mean pressure behind the oral constriction. One lumped cavity with
// compliance C (air plus yielding walls) is filled by the glottal flow and
// vented through the oral and nasal outlets:
//   C dP/dt = U_in - A_out sqrt(2 |P| / rho) sign(P).
// The step is implicit in P and reduces to a quadratic in sqrt(|P|).
class IntraoralPressure final {
public:
    void reset() noexcept;
    double advance(double inflow, double outlet_area, double compliance, double sample_period,
                   const AirProperties& air) noexcept;
    // Linearized mass balance for moving geometry:
    // C_new P_new - C_old P_old = dt*(U_in-U_out) - (V_new-V_old).
    // Compliance includes yielding walls; volume_change is geometric air
    // volume displacement only. Volumes/flows are referenced to ambient.
    double advance_moving(double inflow,double outlet_area,double old_compliance,
                          double new_compliance,double volume_change,double sample_period,
                          const AirProperties& air) noexcept;
    [[nodiscard]] double pressure() const noexcept;

private:
    double pressure_{};
};

// Deterministic Gaussian noise (xorshift64* with Box-Muller) for turbulence
// and for slow physiological variation.
class NoiseGenerator final {
public:
    void seed(std::uint64_t value) noexcept;
    [[nodiscard]] double uniform() noexcept;   // [0, 1)
    [[nodiscard]] double gaussian() noexcept;  // unit variance

private:
    std::uint64_t state_{0x9E3779B97F4A7C15ULL};
    bool has_spare_{};
    double spare_{};
};

// Unit-variance noise low-passed by three cascaded poles: a slowly wandering
// physiological variable (pitch jitter, pressure shimmer). Muscles and lungs
// cannot move at audio rates; a single pole would leave a white residue
// (-6 dB/octave only) that reaches the air as hiss wherever the variable
// drives the flow directly. Three poles fall 18 dB/octave.
class SlowNoise final {
public:
    void prepare(std::uint64_t seed, double cutoff_hz, double update_rate_hz) noexcept;
    double next() noexcept;

private:
    NoiseGenerator source_{};
    double coefficient_{};
    double normalization_{};
    std::array<double, 3> stages_{};
};

// Turbulent source noise: Gaussian noise through a one-pole high-pass, which
// removes low-frequency rumble that a turbulent jet does not produce.
class TurbulenceNoise final {
public:
    void prepare(std::uint64_t seed, double highpass_hz, double sample_rate_hz) noexcept;
    double next() noexcept;

private:
    NoiseGenerator source_{};
    double coefficient_{};
    double previous_input_{};
    double previous_output_{};
};

}
