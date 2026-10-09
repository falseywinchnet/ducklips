#include "ducklips/continuous_modes.hpp"
#include "real_dot.hpp"
#include <cmath>
#include <stdexcept>

namespace ducklips {
void ContinuousModesReference::prepare(std::size_t modes) {
    if (modes==0 || modes>4096) { throw std::invalid_argument("invalid continuous mode count"); }
    modes_=modes;state.assign(2*modes,0.0);
    for (std::vector<double>& stage:k_) { stage.resize(2*modes); }
    work_.resize(2*modes);rotated_.resize(2*modes);dual_.resize(modes);
    half_cos_.resize(modes);half_sin_.resize(modes);whole_cos_.resize(modes);whole_sin_.resize(modes);
}
void ContinuousModesReference::coupling(const ModalMotionSample& motion,std::span<const double> input,
                                        std::span<const double> cosine,std::span<const double> sine,
                                        std::span<double> output) noexcept {
    if (!cosine.empty()) {
        for (std::size_t i=0;i<modes_;++i) {
            rotated_[i]=cosine[i]*input[i]-sine[i]*input[modes_+i];
            rotated_[modes_+i]=sine[i]*input[i]+cosine[i]*input[modes_+i];
        }
        input=rotated_;
    }
    if (motion.flow_transport.empty()) {
        for (std::size_t i=0;i<modes_;++i) {
            dual_[i]=motion.flow_diagonal[i]*input[modes_+i]/motion.angular_frequency[i];
        }
        for (std::size_t i=0;i<modes_;++i) {
            const double scale=-input[modes_+i]/motion.angular_frequency[i];
            output[i]=detail::dot_and_accumulate(motion.pressure_transport.data()+i*modes_,
                                                 input.data(),dual_.data(),scale,modes_);
        }
        for (std::size_t i=0;i<modes_;++i) {
            const double a=output[i],b=motion.angular_frequency[i]*dual_[i];
            if (cosine.empty()) { output[modes_+i]=b; }
            else {
                output[i]=cosine[i]*a+sine[i]*b;
                output[modes_+i]=-sine[i]*a+cosine[i]*b;
            }
        }
        return;
    }
    for (std::size_t i=0;i<modes_;++i) {
        const std::size_t row=i*modes_;
        const double a=detail::real_dot(motion.pressure_transport.data()+row,input.data(),modes_);
        const double b=detail::real_dot(motion.flow_transport.data()+row,input.data()+modes_,modes_);
        if (cosine.empty()) { output[i]=a;output[modes_+i]=b; }
        else {
            output[i]=cosine[i]*a+sine[i]*b;
            output[modes_+i]=-sine[i]*a+cosine[i]*b;
        }
    }
}
void ContinuousModesReference::advance(const ModalMotionSample& first,const ModalMotionSample& middle,
                                       const ModalMotionSample& last,double dt) noexcept {
    for (std::size_t i=0;i<modes_;++i) {
        const double half=dt*(5*first.angular_frequency[i]+8*middle.angular_frequency[i]-last.angular_frequency[i])/24;
        const double whole=dt*(first.angular_frequency[i]+4*middle.angular_frequency[i]+last.angular_frequency[i])/6;
        half_cos_[i]=std::cos(half);half_sin_[i]=std::sin(half);
        const double correction=dt*(last.angular_frequency[i]-first.angular_frequency[i])/4;
        if (std::abs(correction)<=1e-3) {
            // theta_full = 2*theta_half + correction. Slow metric motion
            // makes the correction tiny; reuse the expensive half rotation.
            // The omitted cosine term is <= 1e-18/720, below double epsilon.
            const double square=correction*correction;
            const double c=1-square/2+square*square/24;
            const double s=correction*(1-square/6+square*square/120);
            const double doubled_cos=half_cos_[i]*half_cos_[i]-half_sin_[i]*half_sin_[i];
            const double doubled_sin=2*half_cos_[i]*half_sin_[i];
            whole_cos_[i]=doubled_cos*c-doubled_sin*s;
            whole_sin_[i]=doubled_sin*c+doubled_cos*s;
        } else {
            whole_cos_[i]=std::cos(whole);whole_sin_[i]=std::sin(whole);
        }
    }
    coupling(first,state,{}, {},k_[0]);
    for (std::size_t i=0;i<state.size();++i) { work_[i]=state[i]+dt*.5*k_[0][i]; }
    coupling(middle,work_,half_cos_,half_sin_,k_[1]);
    for (std::size_t i=0;i<state.size();++i) { work_[i]=state[i]+dt*.5*k_[1][i]; }
    coupling(middle,work_,half_cos_,half_sin_,k_[2]);
    for (std::size_t i=0;i<state.size();++i) { work_[i]=state[i]+dt*k_[2][i]; }
    coupling(last,work_,whole_cos_,whole_sin_,k_[3]);
    for (std::size_t i=0;i<modes_;++i) {
        const std::size_t j=modes_+i;
        const double a=state[i]+dt*(k_[0][i]+2*k_[1][i]+2*k_[2][i]+k_[3][i])/6;
        const double b=state[j]+dt*(k_[0][j]+2*k_[1][j]+2*k_[2][j]+k_[3][j])/6;
        state[i]=whole_cos_[i]*a-whole_sin_[i]*b;
        state[j]=whole_sin_[i]*a+whole_cos_[i]*b;
    }
}
}
