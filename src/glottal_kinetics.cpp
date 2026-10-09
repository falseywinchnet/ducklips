#include "ducklips/glottal_kinetics.hpp"
#include <cmath>

namespace ducklips {
namespace {
bool positive(double x) noexcept{return std::isfinite(x) && x>0;}
}
bool GlottalKinetics::configure(const FoldParameters& f,const AirProperties& air,double outlet_loss) noexcept {
    if(!positive(f.length) || !positive(f.lower_thickness) || !positive(f.upper_thickness) || !positive(air.density)
        || !std::isfinite(air.viscosity) || air.viscosity<0 || !std::isfinite(outlet_loss) || outlet_loss<0){return false;}
    density_=air.density;viscosity_=air.viscosity;width_=f.length;thickness_={f.lower_thickness,f.upper_thickness};outlet_loss_=outlet_loss;return true;
}
bool GlottalKinetics::evaluate(const std::array<double,2>& mass,const std::array<double,3>& momentum,GlottalKineticState& result) const noexcept {
    if(!positive(density_) || !positive(mass[0]) || !positive(mass[1])){return false;}
    for(double p:momentum){if(!std::isfinite(p)){return false;}}
    const double a=thickness_[0]*thickness_[0],b=thickness_[1]*thickness_[1];
    const double l0=.5*density_*a/mass[0],l1=.5*density_*b/mass[1];
    GlottalKineticState s;s.inverse_inertance={1/l0,1/(l0+l1),1/l1};
    for(std::size_t j=0;j<3;++j){if(!positive(s.inverse_inertance[j])){return false;}s.energy+=.5*momentum[j]*momentum[j]*s.inverse_inertance[j];}
    // Each chamber contributes half its longitudinal kinetic energy at
    // each adjoining passage. The serial middle passage has common mass flow.
    const double qmid=momentum[1]*s.inverse_inertance[1];
    s.mass_effort[0]=momentum[0]*momentum[0]/(density_*a)+.25*density_*a*(qmid/mass[0])*(qmid/mass[0]);
    s.mass_effort[1]=momentum[2]*momentum[2]/(density_*b)+.25*density_*b*(qmid/mass[1])*(qmid/mass[1]);
    if(!std::isfinite(s.energy) || !std::isfinite(s.mass_effort[0]) || !std::isfinite(s.mass_effort[1])){return false;}
    result=s;return true;
}
bool GlottalKinetics::transfer(const std::array<double,2>& m0,const std::array<double,2>& m1,
    const std::array<double,3>& p0,const std::array<double,3>& p1,GlottalKineticTransfer& result) const noexcept {
    GlottalKineticState a,b;if(!evaluate(m0,p0,a) || !evaluate(m1,p1,b)){return false;}
    GlottalKineticTransfer s;
    for(std::size_t j=0;j<3;++j){s.flow[j]=.25*(a.inverse_inertance[j]+b.inverse_inertance[j])*(p0[j]+p1[j]);}
    const double h0=thickness_[0]*thickness_[0],h1=thickness_[1]*thickness_[1];
    // Symmetrized coordinate secants of G_mid = 2 M0 M1 /
    // (rho0 (h0^2 M1 + h1^2 M0)). These rational expressions avoid
    // subtracting nearly equal G values and give an exact finite chain rule.
    const double d00=h0*m0[1]+h1*m0[0],d10=h0*m0[1]+h1*m1[0];
    const double d01=h0*m1[1]+h1*m0[0],d11=h0*m1[1]+h1*m1[0];
    const double dg0=h0/density_*(m0[1]*m0[1]/(d00*d10)+m1[1]*m1[1]/(d01*d11));
    const double dg1=h1/density_*(m0[0]*m0[0]/(d00*d01)+m1[0]*m1[0]/(d10*d11));
    const double middle=.25*(p0[1]*p0[1]+p1[1]*p1[1]);
    s.mass_effort[0]=.5*(p0[0]*p0[0]+p1[0]*p1[0])/(density_*h0)+middle*dg0;
    s.mass_effort[1]=.5*(p0[2]*p0[2]+p1[2]*p1[2])/(density_*h1)+middle*dg1;
    s.energy_change=b.energy-a.energy;s.chain_residual=s.energy_change;
    for(std::size_t i=0;i<2;++i){s.chain_residual-=s.mass_effort[i]*(m1[i]-m0[i]);}
    for(std::size_t j=0;j<3;++j){s.chain_residual-=s.flow[j]*(p1[j]-p0[j]);}
    if(!std::isfinite(s.chain_residual)){return false;}result=s;return true;
}
bool GlottalKinetics::losses(const std::array<double,2>& mass,const std::array<double,2>& volume,
    std::array<double,3>& resistance,std::array<double,3>& quadratic) const noexcept {
    if(!positive(density_)){return false;}std::array<double,2> r;
    for(std::size_t i=0;i<2;++i){if(!positive(mass[i]) || !positive(volume[i])){return false;}
        const double area=volume[i]/thickness_[i],relative_density=mass[i]/volume[i];
        r[i]=6*viscosity_*thickness_[i]*width_*width_/(area*area*area*relative_density*relative_density);
        if(!std::isfinite(r[i]) || r[i]<0){return false;}
    }
    const double b=.5*outlet_loss_*density_*(thickness_[1]/mass[1])*(thickness_[1]/mass[1]);
    if(!std::isfinite(b) || !std::isfinite(r[0]+r[1])){return false;}
    resistance={r[0],r[0]+r[1],r[1]};quadratic={0.,0.,b};return true;
}
}
