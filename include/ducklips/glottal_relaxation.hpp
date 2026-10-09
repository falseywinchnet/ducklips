#pragma once
#include <array>
namespace ducklips {
// A = inverse inertia, R = linear resistance, A = L L^T. Return
// B = L beta(step L^T R L) L^T, beta(x) = coth(x/2)-2/x.
// Inactive coordinates are exactly eliminated. For PSD R, 0 <= B <= A.
// Units of step and A R must be reciprocal. No allocations; failure leaves B.
[[nodiscard]] bool glottal_relaxation_bias(const std::array<double,49>& inverse_metric,
    const std::array<double,49>& resistance,const std::array<bool,7>& active,
    double step,std::array<double,49>& bias) noexcept;
}
