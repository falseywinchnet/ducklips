#pragma once
#include "ducklips/acoustic_network.hpp"
#include "ducklips/tree_operator.hpp"
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {

// The lowest acoustic eigenstates of an AcousticNetwork: K phi = w^2 C phi.
// Pressure shapes are C-orthonormal (phi_k^T C phi_j = delta_kj), so the
// network's acoustic energy is (1/2) sum_k (a_k^2 + b_k^2) in modal
// coordinates. Layout of shape arrays: [mode * node_count + node].
struct ModalBasis final {
    std::size_t node_count{};
    std::size_t capacity{};
    std::size_t mode_count{};
    std::vector<double> eigenvalue{};        // w^2, rad^2/s^2, ascending
    std::vector<double> angular_frequency{}; // rad/s
    std::vector<double> damping{};           // amplitude decay rate, 1/s
    std::vector<double> pressure_shape{};    // phi
    std::vector<double> unit_shape{};        // C^1/2 phi, unit Euclidean norm

    void prepare(std::size_t node_count, std::size_t capacity);
    [[nodiscard]] double pressure(std::size_t mode, std::size_t node) const noexcept;
};

struct ModalSolveReport final {
    std::size_t warm_modes{};       // accepted from the previous basis by Rayleigh refinement
    std::size_t bisected_modes{};   // placed by inertia bisection and inverse iteration
    std::size_t factorizations{};   // O(node_count) tree factorizations performed
};

// Finds every eigenstate below a frequency limit (up to the basis capacity).
//
// Each mode k is first refined from mode k of the previous basis by Rayleigh
// quotient iteration; the first quotient is the first-order perturbation
// estimate of the deformed eigenvalue. The inertia of S - shift I then
// verifies that the refined value is the k-th eigenvalue. A mode that fails
// verification is placed by inertia bisection and recovered by inverse
// iteration with deflation against nearly equal neighbours. No mode tracking
// across frames is needed for correctness: the resonator reprojects physical
// quantities.
class ModalSolver final {
public:
    void prepare(std::size_t node_count);
    // matrix must hold network.write_operator(). previous may be empty.
    ModalSolveReport solve(const AcousticNetwork& network, TreeOperator& matrix, double maximum_angular_frequency,
                           const ModalBasis& previous, ModalBasis& output) noexcept;
    // Amplitude decay of each mode from branch resistances (viscous, wall,
    // radiation) and node thermal conductance, evaluated at the mode frequency.
    static void compute_damping(const AcousticNetwork& network, ModalBasis& basis) noexcept;

private:
    bool refine_rayleigh(TreeOperator& matrix, std::span<double> vector, double& eigenvalue,
                         ModalSolveReport& report) noexcept;
    bool index_verified(TreeOperator& matrix, std::size_t index, double eigenvalue, ModalSolveReport& report) noexcept;
    double bisect_eigenvalue(TreeOperator& matrix, std::size_t index, double upper, ModalSolveReport& report) noexcept;
    void inverse_iterate(TreeOperator& matrix, double eigenvalue, std::span<double> vector, const ModalBasis& output,
                         std::size_t index, ModalSolveReport& report) noexcept;
    double rayleigh_quotient(const TreeOperator& matrix, std::span<const double> vector) noexcept;
    // Removes the components along every accepted mode below index.
    static void deflate(std::span<double> vector, const ModalBasis& output, std::size_t index) noexcept;

    std::vector<double> work_{};
    std::vector<double> product_{};
};

// Branch flow shape of a mode, M-orthonormal: psi_k(j) = g_j (phi_k(from) - phi_k(to)) / w_k.
[[nodiscard]] double branch_flow_shape(const AcousticNetwork& network, const ModalBasis& basis,
                                       std::size_t mode, std::size_t branch) noexcept;

// What a geometry change holds fixed when the state moves to new modes.
//   passive:         node value min(C_old, C_new) p and branch value
//                    min(g_old, g_new) M_old u. Each variable keeps whichever of
//                    its two invariants (pressure or charge, flow or momentum)
//                    gives the lower energy, so by Bessel's inequality in the
//                    new metrics the projected energy never exceeds the old.
//   charge_momentum: node charge C p and branch momentum M u (adiabatic walls).
//                    No energy bound: a fixed charge in a vanishing volume
//                    has unbounded pressure (HYPOTHESIS for exact closures;
//                    the measured lip closure lost energy, ratio 0.979).
//   pressure_flow:   node pressure and branch flow. No energy bound; the
//                    measured lip closure gained energy (ratio 1.12).
enum class ReprojectionInvariant { passive, charge_momentum, pressure_flow };

// Modal state and the per-block coefficients of the acoustic resonator.
//
// For mode k with z_k = a_k + i b_k:
//   dz/dt = (-alpha_k + i w_k) z + (phi_k^T s_node) + i (psi_k^T s_series)
// where s_node are volume-flow injections into nodes and s_series are
// pressure sources in series with branches. Over one sample with the inputs
// held, z <- propagator z + input_gain * input, exactly.
class ModalResonator final {
public:
    void prepare(std::size_t capacity, std::size_t node_count, std::size_t branch_count);
    void clear() noexcept;
    // Writes coefficients for basis. Does not change the state.
    void configure(const AcousticNetwork& network, const ModalBasis& basis, double sample_period,
                   std::size_t aspiration_branch, std::size_t frication_branch) noexcept;
    // Carries the state from the old network and basis to the new ones across
    // an instantaneous geometry change, then projects it onto the new modes.
    void reproject(const AcousticNetwork& old_network, const ModalBasis& old_basis,
                   const AcousticNetwork& new_network, const ModalBasis& new_basis,
                   ReprojectionInvariant invariant) noexcept;
    [[nodiscard]] double energy() const noexcept;

    std::size_t mode_count{};
    std::vector<std::complex<double>> state{};
    std::vector<std::complex<double>> propagator{};
    std::vector<std::complex<double>> input_gain{};
    std::vector<double> glottal_weight{};    // phi_k(0)
    std::vector<double> aspiration_weight{}; // psi_k(aspiration branch)
    std::vector<double> frication_weight{};  // psi_k(frication branch)
    std::vector<double> lip_output{};        // d(lip flow)/dt per unit a_k
    std::vector<double> nostril_output{};
    std::vector<double> wall_output{};       // summed oral wall flow derivative per unit a_k
    // Response of the next-sample glottal pressure to this sample's glottal
    // flow: sum_k phi_k(0)^2 Re(input_gain_k), Pa s/m^3.
    double glottal_impedance{};
    // Lip branch conductance and whether the frication source sits in that
    // branch; a source there also drives the radiated flow directly:
    // d(lip flow)/dt = g_lip (p_last + p_source).
    double lip_conductance{};
    bool frication_at_lips{};

private:
    std::vector<double> charge_{};   // per node
    std::vector<double> momentum_{}; // per branch
};

}
