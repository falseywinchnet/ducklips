#include "ducklips/radiation_boundary.hpp"
#include <cmath>
#include <initializer_list>
namespace ducklips {
bool RadiationBoundary::prepare(const AirProperties& air,double area,double correction,double factor) noexcept {
    for(double v:{air.density,air.sound_speed,correction,factor}){if(!std::isfinite(v) || v<=0){return false;}}
    RadiationBoundary candidate;candidate.inverse_scale_=std::sqrt(pi)/(air.density*correction);
    candidate.coefficient_=factor*air.density/(4*pi*air.sound_speed);
    double work;if(!candidate.change_area(area,work)){return false;}*this=candidate;return true;
}
bool RadiationBoundary::change_area(double area,double& work) noexcept {
    if(!std::isfinite(area) || area<0 || inverse_scale_<=0 || coefficient_<=0 || !std::isfinite(momentum)){return false;}
    const double inverse=inverse_scale_*std::sqrt(area);
    const double before=energy(),after=.5*inverse*momentum*momentum;
    if(!std::isfinite(inverse) || !std::isfinite(before) || !std::isfinite(after)){return false;}
    work=after-before;area_=area;inverse_=inverse;if(area==0){momentum=0.;}return true;
}
double RadiationBoundary::resistance() const noexcept {
    return inverse_>0 ? 1/(coefficient_*inverse_*inverse_) : 0.;
}
bool RadiationBoundary::advance(double& duct,double g,double dt) noexcept {
    if(inverse_scale_<=0 || coefficient_<=0 || !std::isfinite(duct) || !std::isfinite(momentum)
        || !std::isfinite(g) || g<0 || !std::isfinite(dt) || dt<=0){return false;}
    RadiationBoundaryStep receipt;
    if(area_==0){if(g!=0 || momentum!=0){return false;}last=receipt;return true;}
    const double metric=g+inverse_,difference=g*duct-inverse_*momentum;
    const double rate=resistance()*metric;
    // Infinite positive rate is the limiting fully relaxed mode. expm1
    // also preserves the small-rate increment without subtracting two ones.
    if(!std::isfinite(metric) || !std::isfinite(difference) || std::isnan(rate) || rate<=0){return false;}
    const double exchanged=difference/metric*(-std::expm1(-dt*rate));
    const double next_duct=duct-exchanged,next_momentum=momentum+exchanged;
    receipt.radiated_energy=.5*difference*(difference/metric)*(-std::expm1(-2*dt*rate));
    const double before=.5*g*duct*duct+energy();
    const double after=.5*g*next_duct*next_duct+.5*inverse_*next_momentum*next_momentum;
    receipt.energy_change=after-before;receipt.energy_residual=receipt.energy_change+receipt.radiated_energy;
    if(!std::isfinite(next_duct) || !std::isfinite(next_momentum) || !std::isfinite(receipt.radiated_energy)
        || !std::isfinite(receipt.energy_residual) || receipt.radiated_energy<0){return false;}
    duct=next_duct;momentum=next_momentum;last=receipt;return true;
}
}
