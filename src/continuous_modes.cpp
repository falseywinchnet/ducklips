#include "ducklips/continuous_modes.hpp"
#include "real_dot.hpp"
#include <cmath>
#include <stdexcept>

namespace ducklips {
void ContinuousModesReference::prepare(std::size_t modes) {
    if (modes==0 || modes>4096) { throw std::invalid_argument("invalid continuous mode count"); }
    for (std::size_t i=0;i<3;++i) {
        GaussSample& sample=gauss_[i];sample.omega.resize(modes);sample.pressure.resize(modes*modes);
        sample.flow.resize(modes*modes);sample.diagonal.resize(modes);
        sample.cosine.resize(modes);sample.sine.resize(modes);
        gauss_stage_[i].resize(2*modes);gauss_next_[i].resize(2*modes);
    }
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
void ContinuousModesReference::advance_gauss(const ModalMotionSample& first,const ModalMotionSample& middle,
                                             const ModalMotionSample& last,double dt) noexcept {
    constexpr double root15=3.8729833462074168852;
    constexpr double nodes[3]={.5-root15/10,.5,.5+root15/10};
    constexpr double weights[3]={5.0/18,4.0/9,5.0/18};
    constexpr double tableau[3][3]={
        {5.0/36,2.0/9-root15/15,5.0/36-root15/30},
        {5.0/36+root15/24,2.0/9,5.0/36-root15/24},
        {5.0/36+root15/30,2.0/9+root15/15,5.0/36}};
    const bool canonical=first.flow_transport.empty();
    std::array<ModalMotionSample,3> views;
    for (std::size_t stage=0;stage<3;++stage) {
        GaussSample& sample=gauss_[stage];const double c=nodes[stage],square=c*c,cube=square*c;
        const double l0=2*square-3*c+1,lm=4*c-4*square,l1=2*square-c;
        for (std::size_t i=0;i<modes_;++i) {
            sample.omega[i]=l0*first.angular_frequency[i]+lm*middle.angular_frequency[i]+l1*last.angular_frequency[i];
            const double phase=dt*((2*cube/3-1.5*square+c)*first.angular_frequency[i]
                +(2*square-4*cube/3)*middle.angular_frequency[i]+(2*cube/3-.5*square)*last.angular_frequency[i]);
            sample.cosine[i]=std::cos(phase);sample.sine[i]=std::sin(phase);
            if (canonical) { sample.diagonal[i]=l0*first.flow_diagonal[i]+lm*middle.flow_diagonal[i]+l1*last.flow_diagonal[i]; }
        }
        // The center is already supplied. Only the two outer matrices need
        // interpolation; all three Picard passes reuse the same coefficients.
        if (stage==1) { views[stage]=middle; }
        else {
            for (std::size_t i=0;i<modes_*modes_;++i) {
                sample.pressure[i]=l0*first.pressure_transport[i]+lm*middle.pressure_transport[i]+l1*last.pressure_transport[i];
                if (!canonical) { sample.flow[i]=l0*first.flow_transport[i]+lm*middle.flow_transport[i]+l1*last.flow_transport[i]; }
            }
            views[stage]={sample.omega,sample.pressure,canonical ? std::span<const double>{} : std::span<const double>{sample.flow},sample.diagonal};
        }
        coupling(views[stage],state,sample.cosine,sample.sine,gauss_stage_[stage]);
    }
    for (std::size_t iteration=0;iteration<2;++iteration) {
        for (std::size_t stage=0;stage<3;++stage) {
            for (std::size_t i=0;i<state.size();++i) {
                work_[i]=state[i]+dt*(tableau[stage][0]*gauss_stage_[0][i]
                    +tableau[stage][1]*gauss_stage_[1][i]+tableau[stage][2]*gauss_stage_[2][i]);
            }
            coupling(views[stage],work_,gauss_[stage].cosine,gauss_[stage].sine,gauss_next_[stage]);
        }
        gauss_stage_.swap(gauss_next_);
    }
    for (std::size_t i=0;i<modes_;++i) {
        const std::size_t j=modes_+i;
        const double a=state[i]+dt*(weights[0]*gauss_stage_[0][i]+weights[1]*gauss_stage_[1][i]+weights[2]*gauss_stage_[2][i]);
        const double b=state[j]+dt*(weights[0]*gauss_stage_[0][j]+weights[1]*gauss_stage_[1][j]+weights[2]*gauss_stage_[2][j]);
        const double phase=dt*(first.angular_frequency[i]+4*middle.angular_frequency[i]+last.angular_frequency[i])/6;
        const double cosine=std::cos(phase),sine=std::sin(phase);
        state[i]=cosine*a-sine*b;state[j]=sine*a+cosine*b;
    }
}

}
