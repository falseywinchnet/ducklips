#include "ducklips/slit_thermal.hpp"
#include "ducklips/thermal_spectrum.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace ducklips {
bool SlitThermal::prepare(const AirProperties& air,double compliance,double half_gap) noexcept {
    if(!std::isfinite(air.heat_capacity_ratio) || air.heat_capacity_ratio<=1
        || !std::isfinite(air.thermal_diffusivity) || air.thermal_diffusivity<=0){return false;}
    SlitThermal next;next.diffusivity_=air.thermal_diffusivity;
    std::array<double,modes*modes> matrix{};
    for(std::size_t i=0;i<modes;++i){
        next.capacity_[i]=(air.heat_capacity_ratio-1)*thermal_spectrum::weight[i];
        next.capacity_sum_+=next.capacity_[i];
        next.scale_[i]=std::sqrt(next.capacity_[i]*thermal_spectrum::pole[i]);
        next.vectors_[modes*i+i]=1.;
    }
    // Eliminate the exactly conserved total charge. In y_j=(z_j-z_eq)/
    // sqrt(kappa_j*pole_j), the remaining generator is SPD:
    // S=diag(pole)+sqrt(kappa*pole)*sqrt(kappa*pole)^T.
    for(std::size_t i=0;i<modes;++i){for(std::size_t j=0;j<modes;++j){
        matrix[modes*i+j]=next.scale_[i]*next.scale_[j]+(i==j?thermal_spectrum::pole[i]:0.);
    }}
    constexpr double epsilon=std::numeric_limits<double>::epsilon();bool converged=false;
    for(unsigned sweep=0;sweep<80;++sweep){bool changed=false;
        for(std::size_t p=0;p<modes;++p){for(std::size_t q=p+1;q<modes;++q){
            const double off=matrix[modes*p+q],pp=matrix[modes*p+p],qq=matrix[modes*q+q];
            if(off==0 || std::abs(off)<=4*epsilon*std::sqrt(pp)*std::sqrt(qq)){continue;}
            changed=true;const double tau=(qq-pp)/(2*off);
            const double tangent=std::copysign(1.,tau)/(std::abs(tau)+std::hypot(1.,tau));
            const double cosine=1/std::sqrt(1+tangent*tangent),sine=tangent*cosine;
            for(std::size_t k=0;k<modes;++k){if(k==p || k==q){continue;}
                const double kp=matrix[modes*k+p],kq=matrix[modes*k+q];
                matrix[modes*k+p]=matrix[modes*p+k]=cosine*kp-sine*kq;
                matrix[modes*k+q]=matrix[modes*q+k]=sine*kp+cosine*kq;
            }
            matrix[modes*p+p]=pp-tangent*off;matrix[modes*q+q]=qq+tangent*off;
            matrix[modes*p+q]=matrix[modes*q+p]=0.;
            for(std::size_t k=0;k<modes;++k){const double kp=next.vectors_[modes*k+p],kq=next.vectors_[modes*k+q];
                next.vectors_[modes*k+p]=cosine*kp-sine*kq;next.vectors_[modes*k+q]=sine*kp+cosine*kq;
            }
        }}
        if(!changed){converged=true;break;}
    }
    if(!converged){return false;}
    for(std::size_t i=0;i<modes;++i){next.rate_[i]=matrix[modes*i+i];if(!std::isfinite(next.rate_[i]) || next.rate_[i]<=0){return false;}}
    double work;if(!next.change_geometry(compliance,half_gap,work)){return false;}*this=next;return true;
}
double SlitThermal::stored_air() const noexcept {double result=0;for(double z:charge){result+=z;}return result;}
double SlitThermal::energy() const noexcept {
    if(compliance_==0){return 0;}double result=0;
    for(std::size_t i=0;i<modes;++i){result+=charge[i]*charge[i]/capacity_[i];}
    return .5*result/compliance_;
}
bool SlitThermal::change_geometry(double compliance,double half_gap,double& work) noexcept {
    if(!std::isfinite(compliance) || compliance<=0 || !std::isfinite(half_gap) || half_gap<=0 || diffusivity_<=0){return false;}
    double quadratic=0;for(std::size_t i=0;i<modes;++i){
        if(!std::isfinite(charge[i]) || capacity_[i]<=0){return false;}quadratic+=charge[i]*charge[i]/capacity_[i];
    }
    const double before=compliance_>0?.5*quadratic/compliance_:0.,after=.5*quadratic/compliance;
    if(!std::isfinite(before) || !std::isfinite(after) || !std::isfinite(diffusivity_/half_gap/half_gap)){return false;}
    work=after-before;compliance_=compliance;half_gap_=half_gap;return true;
}
bool SlitThermal::equilibrate(double pressure) noexcept {
    if(!std::isfinite(pressure) || compliance_<=0){return false;}
    std::array<double,modes> next;for(std::size_t i=0;i<modes;++i){next[i]=capacity_[i]*compliance_*pressure;if(!std::isfinite(next[i])){return false;}}
    charge=next;return true;
}
bool SlitThermal::advance(double& primary,double dt) noexcept {
    if(!std::isfinite(primary) || !std::isfinite(dt) || dt<=0 || compliance_<=0){return false;}
    const double total=primary+stored_air(),equilibrium=total/(1+capacity_sum_),time=dt*diffusivity_/half_gap_/half_gap_;
    if(!std::isfinite(total) || !std::isfinite(time)){return false;}
    std::array<double,modes> deviation{},increment{},next{};
    for(std::size_t i=0;i<modes;++i){deviation[i]=(charge[i]-capacity_[i]*equilibrium)/scale_[i];}
    SlitThermalStep receipt;
    for(std::size_t j=0;j<modes;++j){double amplitude=0;
        for(std::size_t i=0;i<modes;++i){amplitude+=vectors_[modes*i+j]*deviation[i];}
        const double change=amplitude*std::expm1(-time*rate_[j]);
        for(std::size_t i=0;i<modes;++i){increment[i]+=vectors_[modes*i+j]*change;}
        receipt.dissipated_energy+=.5*rate_[j]*amplitude*amplitude/compliance_*(-std::expm1(-2*time*rate_[j]));
    }
    double transferred=0,after=0;
    for(std::size_t i=0;i<modes;++i){const double change=scale_[i]*increment[i];next[i]=charge[i]+change;transferred+=change;
        if(!std::isfinite(next[i])){return false;}after+=.5*next[i]*next[i]/capacity_[i]/compliance_;
    }
    const double next_primary=primary-transferred,before=.5*primary*primary/compliance_+energy();
    after+=.5*next_primary*next_primary/compliance_;
    receipt.energy_change=after-before;receipt.energy_residual=receipt.energy_change+receipt.dissipated_energy;
    double next_total=next_primary;for(double z:next){next_total+=z;}receipt.stored_air_residual=next_total-total;
    if(!std::isfinite(next_primary) || !std::isfinite(receipt.energy_residual) || !std::isfinite(receipt.stored_air_residual) || receipt.dissipated_energy<0){return false;}
    primary=next_primary;charge=next;last=receipt;return true;
}
}
