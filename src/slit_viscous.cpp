#include "ducklips/slit_viscous.hpp"
#include "ducklips/viscous_spectrum.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace ducklips {
bool SlitViscous::diagonalize(double ratio) noexcept {
    constexpr std::size_t n=dimension;std::array<double,n*n> a{};vectors_.fill(0.);a[0]=3*ratio;
    for(std::size_t i=0;i<n;++i){vectors_[n*i+i]=1.;}
    for(std::size_t j=0;j<modes;++j){const double w=viscous_spectrum::weight[j],p=viscous_spectrum::pole[j];
        a[0]+=ratio*w;a[n*(j+1)+j+1]=p;a[j+1]=a[n*(j+1)]=-std::sqrt(ratio*w*p);
    }
    constexpr double epsilon=std::numeric_limits<double>::epsilon();bool converged=false;
    for(unsigned sweep=0;sweep<80;++sweep){bool changed=false;
        for(std::size_t p=0;p<n;++p){for(std::size_t q=p+1;q<n;++q){
            const double off=a[n*p+q],pp=a[n*p+p],qq=a[n*q+q];
            if(off==0 || std::abs(off)<=4*epsilon*std::sqrt(pp)*std::sqrt(qq)){continue;}
            changed=true;const double tau=(qq-pp)/(2*off);
            const double t=std::copysign(1.,tau)/(std::abs(tau)+std::hypot(1.,tau));
            const double c=1/std::sqrt(1+t*t),s=t*c;
            for(std::size_t k=0;k<n;++k){if(k==p || k==q){continue;}const double kp=a[n*k+p],kq=a[n*k+q];
                a[n*k+p]=a[n*p+k]=c*kp-s*kq;a[n*k+q]=a[n*q+k]=s*kp+c*kq;
            }
            a[n*p+p]=pp-t*off;a[n*q+q]=qq+t*off;a[n*p+q]=a[n*q+p]=0;
            for(std::size_t k=0;k<n;++k){const double kp=vectors_[n*k+p],kq=vectors_[n*k+q];
                vectors_[n*k+p]=c*kp-s*kq;vectors_[n*k+q]=s*kp+c*kq;
            }
        }}
        if(!changed){converged=true;break;}
    }
    if(!converged){return false;}
    for(std::size_t i=0;i<n;++i){rate_[i]=a[n*i+i];if(!std::isfinite(rate_[i]) || rate_[i]<=0){return false;}}
    return true;
}
bool SlitViscous::prepare(const AirProperties& air,double inertance,double conductance,double gap) noexcept {
    if(!std::isfinite(air.viscosity) || air.viscosity<=0 || !std::isfinite(air.density) || air.density<=0){return false;}
    SlitViscous next;next.viscosity_=air.viscosity/air.density;double work;
    if(!next.change_geometry(inertance,conductance,gap,work)){return false;}*this=next;return true;
}
double SlitViscous::energy() const noexcept {
    if(inertance_==0){return 0;}double value=0;
    for(std::size_t i=0;i<modes;++i){value+=viscous_spectrum::pole[i]/viscous_spectrum::weight[i]*momentum[i]*momentum[i];}
    return .5*value/inertance_;
}
bool SlitViscous::change_geometry(double inertance,double conductance,double gap,double& work) noexcept {
    if(inertance==0 && conductance==0 && gap==0 && viscosity_>0){
        for(double z:momentum){if(!std::isfinite(z)){return false;}}
        const double before=energy();if(!std::isfinite(before)){return false;}
        momentum.fill(0.);inertance_=conductance_=half_gap_=0.;work=-before;return true;
    }
    if(!std::isfinite(inertance) || inertance<=0 || !std::isfinite(conductance) || conductance<0
        || !std::isfinite(gap) || gap<=0 || viscosity_<=0 || !std::isfinite(viscosity_/gap/gap)
        || !std::isfinite(inertance*conductance) || (conductance>0 && inertance*conductance<=0)){return false;}
    double quadratic=0;for(std::size_t i=0;i<modes;++i){
        if(!std::isfinite(momentum[i])){return false;}quadratic+=viscous_spectrum::pole[i]/viscous_spectrum::weight[i]*momentum[i]*momentum[i];
    }
    const double before=inertance_>0?.5*quadratic/inertance_:0.,after=.5*quadratic/inertance;
    if(!std::isfinite(before) || !std::isfinite(after)){return false;}
    if(conductance>0 && inertance*conductance!=inertance_*conductance_){
        SlitViscous next=*this;if(!next.diagonalize(inertance*conductance)){return false;}
        *this=next;
    }
    work=after-before;inertance_=inertance;conductance_=conductance;half_gap_=gap;return true;
}
bool SlitViscous::equilibrate(double flow) noexcept {
    if(conductance_==0 && viscosity_>0){if(flow!=0){return false;}momentum.fill(0.);return true;}
    if(!std::isfinite(flow) || inertance_<=0){return false;}std::array<double,modes> next;
    for(std::size_t i=0;i<modes;++i){next[i]=inertance_*viscous_spectrum::weight[i]/viscous_spectrum::pole[i]*flow;if(!std::isfinite(next[i])){return false;}}
    momentum=next;return true;
}
bool SlitViscous::advance(double& primary,double dt) noexcept {
    if(!std::isfinite(primary) || !std::isfinite(dt) || dt<=0 || viscosity_<=0){return false;}
    if(inertance_==0){
        for(double z:momentum){if(z!=0){return false;}}
        last={};return true;
    }
    constexpr std::size_t n=dimension;const double time=dt*viscosity_/half_gap_/half_gap_;
    if(!std::isfinite(time)){return false;}std::array<double,n> state{},increment{},metric{};
    if(conductance_==0){
        SlitViscousStep receipt;std::array<double,modes> next;
        for(std::size_t i=0;i<modes;++i){const double z=momentum[i],pole=viscous_spectrum::pole[i];
            if(!std::isfinite(z)){return false;}
            const double fraction=-std::expm1(-time*pole);
            next[i]=z*std::exp(-time*pole);receipt.blocked_impulse+=z*fraction;
            receipt.dissipated_energy+=.5*pole/(inertance_*viscous_spectrum::weight[i])*z*z*(-std::expm1(-2*time*pole));
        }
        const double before=energy();double after=0;
        for(std::size_t i=0;i<modes;++i){after+=.5*viscous_spectrum::pole[i]/(inertance_*viscous_spectrum::weight[i])*next[i]*next[i];}
        receipt.energy_change=after-before;receipt.energy_residual=receipt.energy_change+receipt.dissipated_energy;
        if(!std::isfinite(receipt.energy_residual) || !std::isfinite(receipt.blocked_impulse)){return false;}
        momentum=next;last=receipt;return true;
    }
    metric[0]=std::sqrt(conductance_);state[0]=metric[0]*primary;
    for(std::size_t i=0;i<modes;++i){metric[i+1]=std::sqrt(viscous_spectrum::pole[i]/(inertance_*viscous_spectrum::weight[i]));state[i+1]=metric[i+1]*momentum[i];}
    SlitViscousStep receipt;
    // Use increments for small changes and direct modal reconstruction once
    // all modes decay appreciably; avoid leaving cancellation residue after
    // an otherwise fully relaxed step.
    const bool direct=time*(*std::min_element(rate_.begin(),rate_.end()))>.5;
    for(std::size_t j=0;j<n;++j){double amplitude=0;for(std::size_t i=0;i<n;++i){amplitude+=vectors_[n*i+j]*state[i];}
        const double change=amplitude*(direct?std::exp(-time*rate_[j]):std::expm1(-time*rate_[j]));
        for(std::size_t i=0;i<n;++i){increment[i]+=vectors_[n*i+j]*change;}
        receipt.dissipated_energy+=.5*amplitude*amplitude*(-std::expm1(-2*time*rate_[j]));
    }
    std::array<double,n> next;double before=0,after=0;
    for(std::size_t i=0;i<n;++i){const double value=(direct?0.:state[i])+increment[i];next[i]=value/metric[i];if(!std::isfinite(next[i])){return false;}
        before+=.5*state[i]*state[i];after+=.5*value*value;
    }
    receipt.energy_change=after-before;receipt.energy_residual=receipt.energy_change+receipt.dissipated_energy;
    if(!std::isfinite(receipt.energy_residual) || receipt.dissipated_energy<0){return false;}
    primary=next[0];for(std::size_t i=0;i<modes;++i){momentum[i]=next[i+1];}last=receipt;return true;
}
}
