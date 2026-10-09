#include "ducklips/moving_mode.hpp"
#include <cmath>
#include <stdexcept>

namespace ducklips {
void UniformAreaModeStep::prepare(double omega,double alpha,double rate,double dt) {
    if (!std::isfinite(omega) || !std::isfinite(alpha) || !std::isfinite(rate) || !std::isfinite(dt)
        || omega<0 || alpha<0 || dt<0) { throw std::invalid_argument("Invalid moving-mode parameter"); }
    const double half=rate/2;
    const double square=omega*omega-half*half;
    double cosine=1.0,sine_over_frequency=dt;
    if (square>0) {
        const double frequency=std::sqrt(square);
        cosine=std::cos(frequency*dt);
        sine_over_frequency=std::sin(frequency*dt)/frequency;
    } else if (square<0) {
        const double frequency=std::sqrt(-square);
        cosine=std::cosh(frequency*dt);
        sine_over_frequency=std::sinh(frequency*dt)/frequency;
    }
    const double decay=std::exp(-alpha*dt);
    UniformAreaModeStep step;
    step.aa=decay*(cosine-half*sine_over_frequency);
    step.ab=-decay*omega*sine_over_frequency;
    step.ba=-step.ab;
    step.bb=decay*(cosine+half*sine_over_frequency);
    if (!std::isfinite(square) || !std::isfinite(step.aa) || !std::isfinite(step.ab) || !std::isfinite(step.bb)) {
        throw std::invalid_argument("Moving-mode propagator overflow");
    }
    *this=step;
}
std::complex<double> UniformAreaModeStep::apply(std::complex<double> state) const noexcept {
    return {aa*state.real()+ab*state.imag(),ba*state.real()+bb*state.imag()};
}
}
