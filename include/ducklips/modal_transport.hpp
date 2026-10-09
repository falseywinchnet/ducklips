#pragma once
#include "ducklips/continuous_modes.hpp"
#include "ducklips/modal.hpp"

namespace ducklips {
// Derives continuous transport from a C-orthonormal pressure basis and physical
// parameter rates. Applies to ascending isolated positive eigenvalues with fixed mode
// count; the caller must keep modal signs/order continuous across samples.
// Closures and crossing/degenerate mode groups need a separate representation.
class ModalTransportReference final {
public:
    void prepare(std::size_t modes,std::size_t nodes,std::size_t branches);
    // Rates are dC/dt (m^3/(Pa*s)) and dG/dt (m^4/(kg*s)). No allocation.
    // Rejects nonfinite coefficients/rates, inconsistent sizes, nonpositive
    // frequencies, and eigenvalue gaps <= 1e-8 * max(omega^2). On rejection
    // the last valid output is retained. Does not solve the eigenproblem.
    [[nodiscard]] bool derive(const AcousticNetwork& network,const ModalBasis& basis,
                              std::span<const double> compliance_rate,
                              std::span<const double> conductance_rate) noexcept;
    [[nodiscard]] ModalMotionSample view() const noexcept;
    [[nodiscard]] ModalMotionSample canonical_view() const noexcept;
    [[nodiscard]] std::span<const double> frequency_rate() const noexcept;
private:
    std::size_t modes_{},nodes_{},branches_{};
    std::vector<double> flow_diagonal_{},next_diagonal_{};
    std::vector<double> omega_{},pressure_{},flow_{},frequency_rate_{};
    std::vector<double> weighted_pressure_{},weighted_differences_{};
    std::vector<double> differences_{},next_pressure_{},next_flow_{},next_rate_{};
};
}
