#include "ducklips/gas_storage.hpp"
#include <cmath>

namespace ducklips {
namespace {
double log_sinhc(double x) noexcept {
    x=std::abs(x);
    if(x<.01){const double z=x*x;return z*(1./6+z*(-1./180+z*(1./2835-z/37800)));}
    if(x>20){return x-std::log(2*x)+std::log1p(-std::exp(-2*x));}
    return std::log(std::sinh(x)/x);
}
double log_density(double r) noexcept {return std::abs(r-1)<.5?std::log1p(r-1):std::log(r);}
bool positive(double x) noexcept {return std::isfinite(x) && x>0;}
}
bool GasStorage::configure(const AirProperties& air) noexcept {
    const double bulk=air.density*air.sound_speed*air.sound_speed;
    if(!positive(air.density) || !positive(air.sound_speed) || !positive(bulk)
        || !std::isfinite(air.heat_capacity_ratio) || air.heat_capacity_ratio<=1){return false;}
    bulk_=bulk;gamma_=air.heat_capacity_ratio;return true;
}
double GasStorage::density_energy(double r) const noexcept {
    const double d=r-1;
    if(std::abs(d)<.01){
        // Cancel the constant/linear terms analytically. First term is
        // exactly the ordinary quadratic acoustic compression energy.
        double term=.5*d*d,sum=term;
        for(int order=3;order<=12;++order){term*=d*(gamma_-double(order-1))/double(order);sum+=term;}
        return bulk_*sum;
    }
    return bulk_/(gamma_*(gamma_-1))*(std::expm1(gamma_*log_density(r))-gamma_*d);
}
bool GasStorage::evaluate(double mass,double volume,GasEnergy& result) const noexcept {
    if(!positive(bulk_) || !positive(mass) || !positive(volume)){return false;}
    const double r=mass/volume;if(!positive(r)){return false;}
    const double l=log_density(r);
    const GasEnergy candidate{volume*density_energy(r),bulk_/gamma_*std::expm1(gamma_*l),bulk_/(gamma_-1)*std::expm1((gamma_-1)*l)};
    if(!std::isfinite(candidate.energy) || candidate.energy<0 || !std::isfinite(candidate.pressure) || !std::isfinite(candidate.mass_effort)){return false;}
    result=candidate;return true;
}
bool GasStorage::mass_at_pressure(double volume,double pressure,double& mass) const noexcept {
    if(!positive(bulk_) || !positive(volume) || !std::isfinite(pressure)){return false;}
    const double excess=gamma_*pressure/bulk_;if(excess<=-1){return false;}
    const double candidate=volume*std::exp(std::log1p(excess)/gamma_);
    if(!positive(candidate)){return false;}mass=candidate;return true;
}
bool GasStorage::transfer(double mass0,double volume0,double mass1,double volume1,GasTransfer& result) const noexcept {
    if(!positive(bulk_) || !positive(mass0) || !positive(mass1) || !positive(volume0) || !positive(volume1)){return false;}
    const double r0=mass0/volume0,r1=mass1/volume1;if(!positive(r0) || !positive(r1)){return false;}
    const double l0=log_density(r0),l1=log_density(r1),h=.5*(l1-l0);
    // Secant of the density energy f(r). Homogeneity H=V f(M/V) gives
    // Delta H = mu_bar Delta M - P_bar Delta V exactly. A proportional
    // mass/volume change at ambient density has exactly zero effort.
    const double exponent=.5*(gamma_-1)*(l0+l1)+log_sinhc(gamma_*h)-log_sinhc(h);
    const double mu=bulk_/(gamma_-1)*std::expm1(exponent);
    const double pressure=.5*(r0+r1)*mu-.5*(density_energy(r0)+density_energy(r1));
    const double delta=(r1-r0)/r0;double secant_slope=0;
    if(std::abs(delta)<.001){
        double term=1.,sum=.5;
        for(int order=1;order<=10;++order){term*=delta*(gamma_-1-double(order))/double(order);sum+=term/double(order+2);}
        secant_slope=bulk_*std::exp((gamma_-2)*l0)*sum;
    }else{
        const double endpoint=bulk_/(gamma_-1)*std::expm1((gamma_-1)*l1);
        secant_slope=(endpoint-mu)/(r1-r0);
    }
    const double slope=-r0*r1/volume1*secant_slope;
    if(!std::isfinite(mu) || !std::isfinite(pressure) || !std::isfinite(slope) || slope>0){return false;}
    result={mu,pressure,slope};return true;
}
}
