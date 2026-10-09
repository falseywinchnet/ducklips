#pragma once
#include <complex>

namespace ducklips {
// Exact homogeneous evolution for a lossless-tube mode under uniform
// exponential area motion: C=s*C0, M=M0/s, s=exp(rate*t). The frozen frequency
// stays fixed. With a=sqrt(C)*p and b=sqrt(M)*u, stored energy is (a*a+b*b)/2:
//   a'=-alpha*a-omega*b-rate*a/2
//   b'= omega*a-alpha*b+rate*b/2.
// The rate terms exchange work with the moving wall. They are not artificial
// damping and cannot in general be represented by one complex multiplier.
// This is a reference building block, NOT a nonuniform moving-tract solver.
struct UniformAreaModeStep final {
    double aa{1.0}, ab{}, ba{}, bb{1.0};
    // omega in rad/s, alpha in 1/s, logarithmic area rate in 1/s, dt in s.
    // Rejects nonfinite inputs, negative omega/alpha/dt and overflow.
    void prepare(double omega, double alpha, double log_area_rate, double dt);
    [[nodiscard]] std::complex<double> apply(std::complex<double> state) const noexcept;
};
}
