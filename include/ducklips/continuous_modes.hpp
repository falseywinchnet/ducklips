#pragma once
#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {
// Instantaneous lossless moving-mode coefficients. Matrices are row-major.
// a_dot=-Omega*b+pressure*a; b_dot=Omega*a+flow*b.
// The owner derives transport from physical geometry, with a continuously
// oriented basis. This view does not solve eigenvectors or handle crossings.
struct ModalMotionSample final {
    std::span<const double> angular_frequency;
    std::span<const double> pressure_transport;
    std::span<const double> flow_transport;
    // If flow_transport is empty, use the canonical dual form
    // flow = -Omega pressure^T Omega^-1 + diag(flow_diagonal).
    // Physical lossless transport has flow_diagonal = omega_dot / omega.
    // This form requires strictly positive angular frequencies.
    std::span<const double> flow_diagonal{};
};

// Local interaction-coordinate RK4: integrate the diagonal acoustic phase
// analytically (Simpson for a quadratic frequency trajectory), then integrate
// only the moving transport in the rotating coordinates. A stationary mode
// advances by an exact sine/cosine rotation, without Gauss frequency bias.
// Reference kernel: no forcing, closure policy or physical loss realization.
class ContinuousModesReference final {
public:
    void prepare(std::size_t modes);
    // Coefficients at t, t+dt/2, t+dt in the SAME continuous modal convention.
    // Requires prepared positive mode count, spans n/n*n, finite coefficients,
    // nonnegative frequencies and finite positive dt. Allocation-free.
    void advance(const ModalMotionSample& first,const ModalMotionSample& middle,
                 const ModalMotionSample& last,double dt) noexcept;
    // [a_0...a_n-1, b_0...b_n-1], energy = dot(state,state)/2.
    std::vector<double> state{};
private:
    void coupling(const ModalMotionSample& motion,std::span<const double> input,
                  std::span<const double> cosine,std::span<const double> sine,
                  std::span<double> output) noexcept;
    std::size_t modes_{};
    std::array<std::vector<double>,4> k_{};
    std::vector<double> work_{},rotated_{},dual_{};
    std::vector<double> half_cos_{},half_sin_{},whole_cos_{},whole_sin_{};
};
}
