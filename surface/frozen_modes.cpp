#include "frozen_modes.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace ducklips {
void FrozenModes::load(const std::string& path) {
    std::ifstream file(path);
    std::string magic;
    std::size_t count=0;
    FrozenModes model;
    if (!(file >> magic >> model.rate_ >> count >> model.inputs_ >> model.outputs_) ||
        magic!="DUCKLIPS_FROZEN_MODES_V1" || !std::isfinite(model.rate_) || model.rate_<=0 ||
        count==0 || count>4096 || model.inputs_==0 || model.inputs_>16 ||
        model.outputs_!=model.inputs_+1) {
        throw std::invalid_argument("Invalid frozen modal header");
    }
    std::vector<double> poles_r(count),poles_i(count),drive_r(count*model.inputs_),drive_i(count*model.inputs_),
        observe_r(count*model.outputs_),observe_i(count*model.outputs_);
    model.direct_.resize(model.inputs_*model.outputs_);
    auto read = [&file](double& value) {
        if (!(file >> value) || !std::isfinite(value)) { throw std::invalid_argument("Invalid frozen modal coefficient"); }
    };
    for (double& value : model.direct_) { read(value); }
    for (std::size_t mode=0; mode<count; ++mode) {
        read(poles_r[mode]); read(poles_i[mode]);
        if (poles_i[mode]<0 || std::hypot(poles_r[mode],poles_i[mode])>1+1e-9) {
            throw std::invalid_argument("Unstable or unpaired frozen mode");
        }
        for (std::size_t input=0; input<model.inputs_; ++input) {
            read(drive_r[input*count+mode]); read(drive_i[input*count+mode]);
            if (poles_i[mode]==0 && drive_i[input*count+mode]!=0) { throw std::invalid_argument("Complex drive on real pole"); }
        }
        for (std::size_t output=0; output<model.outputs_; ++output) {
            read(observe_r[output*count+mode]); read(observe_i[output*count+mode]);
            if (poles_i[mode]==0 && observe_i[output*count+mode]!=0) { throw std::invalid_argument("Complex readout on real pole"); }
        }
        if (poles_i[mode]>0) { ++model.complex_count_; }
    }
    file >> std::ws;
    if (!file.eof()) { throw std::invalid_argument("Trailing frozen modal data"); }
    // Real coordinates: [Re(complex modes), Im(complex modes), real modes].
    // A real pole needs one state, and every readout becomes one real dot
    // product. No complex arithmetic or zero imaginary coefficients at runtime.
    const std::size_t states=count+model.complex_count_;
    model.pole_r_.resize(count);model.pole_i_.resize(model.complex_count_);
    model.state_.resize(states);model.drive_.resize(states*model.inputs_);model.observe_.resize(states*model.outputs_);
    std::size_t complex_index=0,real_index=0;
    for (std::size_t original=0;original<count;++original) {
        const bool paired=poles_i[original]>0;
        const std::size_t mode=paired ? complex_index++ : model.complex_count_+real_index++;
        const std::size_t state=paired ? mode : mode+model.complex_count_;
        model.pole_r_[mode]=poles_r[original];
        if (paired) { model.pole_i_[mode]=poles_i[original]; }
        for (std::size_t input=0;input<model.inputs_;++input) {
            model.drive_[input*states+state]=drive_r[input*count+original];
            if (paired) { model.drive_[input*states+state+model.complex_count_]=drive_i[input*count+original]; }
        }
        for (std::size_t output=0;output<model.outputs_;++output) {
            model.observe_[output*states+state]=observe_r[output*count+original];
            if (paired) { model.observe_[output*states+state+model.complex_count_]=-observe_i[output*count+original]; }
        }
    }
    *this=std::move(model);
}
void FrozenModes::reset() noexcept {
    std::fill(state_.begin(),state_.end(),0.0);
}
double FrozenModes::free_pressure(std::size_t output) const noexcept {
    const double* coefficients=observe_.data()+output*state_.size();
    const double* states=state_.data();
    const std::size_t count=state_.size();
    std::size_t state=0;
    double result=0;
#if defined(__aarch64__)
    // Four independent SIMD accumulators; no fast-math. Mode coordinates are
    // independent, so the readout is a real dot product with explicit rounding.
    float64x2_t sum0=vdupq_n_f64(0),sum1=sum0,sum2=sum0,sum3=sum0;
    for (;state+8<=count;state+=8) {
        sum0=vfmaq_f64(sum0,vld1q_f64(coefficients+state),vld1q_f64(states+state));
        sum1=vfmaq_f64(sum1,vld1q_f64(coefficients+state+2),vld1q_f64(states+state+2));
        sum2=vfmaq_f64(sum2,vld1q_f64(coefficients+state+4),vld1q_f64(states+state+4));
        sum3=vfmaq_f64(sum3,vld1q_f64(coefficients+state+6),vld1q_f64(states+state+6));
    }
    result=vaddvq_f64(vaddq_f64(vaddq_f64(sum0,sum1),vaddq_f64(sum2,sum3)));
#else
    double sums[4]{};
    for (;state+4<=count;state+=4) {
        for (std::size_t lane=0;lane<4;++lane) { sums[lane]+=coefficients[state+lane]*states[state+lane]; }
    }
    result=(sums[0]+sums[1])+(sums[2]+sums[3]);
#endif
    for (;state<count;++state) { result+=coefficients[state]*states[state]; }
    return result;
}
double FrozenModes::direct(std::size_t output,std::size_t input) const noexcept {
    return direct_[output*inputs_+input];
}
void FrozenModes::tick(std::span<const double> inputs,std::span<double> outputs) {
    if (inputs.size()!=inputs_ || outputs.size()!=outputs_) { throw std::invalid_argument("Frozen modal port count mismatch"); }
    for (std::size_t output=0;output<outputs_;++output) {
        double value=free_pressure(output);
        for (std::size_t input=0;input<inputs_;++input) { value+=direct(output,input)*inputs[input]; }
        outputs[output]=value;
    }
    for (std::size_t mode=0;mode<complex_count_;++mode) {
        const double re=state_[mode],im=state_[mode+complex_count_];
        state_[mode]=pole_r_[mode]*re-pole_i_[mode]*im;
        state_[mode+complex_count_]=pole_i_[mode]*re+pole_r_[mode]*im;
    }
    for (std::size_t mode=complex_count_;mode<modes();++mode) { state_[mode+complex_count_]*=pole_r_[mode]; }
    for (std::size_t input=0;input<inputs_;++input) {
        const double flow=inputs[input];
        if (flow==0.0) { continue; }
        const std::size_t offset=input*state_.size();
        for (std::size_t state=0;state<state_.size();++state) { state_[state]+=drive_[offset+state]*flow; }
    }
}
double FrozenModes::tick_glottis(double flow) noexcept {
    double radiation=free_pressure(outputs_-1)+direct(outputs_-1,0)*flow;
    for (std::size_t mode=0;mode<complex_count_;++mode) {
        const double re=state_[mode],im=state_[mode+complex_count_];
        state_[mode]=pole_r_[mode]*re-pole_i_[mode]*im+drive_[mode]*flow;
        state_[mode+complex_count_]=pole_i_[mode]*re+pole_r_[mode]*im+drive_[mode+complex_count_]*flow;
    }
    for (std::size_t mode=complex_count_;mode<modes();++mode) {
        const std::size_t state=mode+complex_count_;
        state_[state]=pole_r_[mode]*state_[state]+drive_[state]*flow;
    }
    return radiation;
}
}
