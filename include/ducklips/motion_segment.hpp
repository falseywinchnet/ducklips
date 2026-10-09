#pragma once
#include "ducklips/continuous_modes.hpp"

namespace ducklips {
struct CanonicalMotionKnot final {
    std::span<const double> angular_frequency,frequency_rate,pressure_transport;
};
struct CanonicalMotionFrame final {
    std::vector<double> angular_frequency,pressure_transport,flow_diagonal;
    void prepare(std::size_t modes);
    [[nodiscard]] ModalMotionSample view() const noexcept;
};
// Uniform control knots at -2,-1,0,1,2,3; evaluate interval [0,1]. Pressure
// uses cubic Hermite interpolation with shared fourth-order centered knot
// derivatives, so both pressure and its slope join continuously. Frequency uses values and physical rates at
// -1,0,1 in a quintic Hermite polynomial. Its own derivative supplies the
// canonical flow diagonal. This needs THREE future knots: at 1kHz, 3ms of
// lookahead. It is a segment evaluator, not an online command scheduler.
class CanonicalMotionSegment final {
public:
    void prepare(std::size_t modes);
    // Allocation-free; invalid preparation retains the last valid segment.
    [[nodiscard]] bool configure(const std::array<CanonicalMotionKnot,6>& knots,double interval) noexcept;
    // Frame must be prepared to the same size. Rejects out-of-interval time
    // or a nonfinite/nonpositive frequency. Failure leaves frame unchanged.
    [[nodiscard]] bool evaluate(double fraction,CanonicalMotionFrame& frame) noexcept;
private:
    std::size_t modes_{};
    bool configured_{};
    double interval_{};
    std::vector<double> pressure_{},next_pressure_{},frequency_{},next_frequency_{};
    std::vector<double> evaluated_frequency_{},evaluated_diagonal_{};
};
}
