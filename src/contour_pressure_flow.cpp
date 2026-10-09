#include "ducklips/contour_pressure_flow.hpp"
#include <algorithm>
#include <cmath>
namespace ducklips {
namespace {
constexpr double points8[]={.18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
constexpr double weights8[]={.36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
constexpr double points16[]={.09501250983763744019,.28160355077925891323,.45801677765722738634,.61787624440264374845,.75540440835500303390,.86563120238783174388,.94457502307323257608,.98940093499164993260};
constexpr double weights16[]={.18945061045506849629,.18260341504492358887,.16915651939500253819,.14959598881657673208,.12462897125553387205,.09515851168249278481,.06225352393864789286,.02715245941175409485};
}
bool ContourPressureKinetics::configure(const FoldParameters& f,const AirProperties& air) noexcept {
    ContourPressureStorage storage;if(!storage.configure(f,air) || !std::isfinite(f.length) || f.length<=0 || !std::isfinite(air.viscosity) || air.viscosity<0){return false;}
    storage_=storage;thickness_={f.lower_thickness,f.upper_thickness};density_=air.density;bulk_=air.density*air.sound_speed*air.sound_speed;gamma_=air.heat_capacity_ratio;viscosity_=air.viscosity;width_=f.length;return true;
}
bool ContourPressureKinetics::sample(const ContourPressureFlow& state,double z,ContourFlowSample& result) const noexcept {
    if(density_<=0 || !std::isfinite(z) || z<state.support[0] || z>state.support[1]){return false;}
    const double c0=.5*thickness_[0],c1=thickness_[0]+.5*thickness_[1],end=thickness_[0]+thickness_[1];
    const double t=std::clamp((z-c0)/(c1-c0),0.,1.),a=(1-t)*state.area[0]+t*state.area[1],mu=(1-t)*state.effort[0]+t*state.effort[1];
    if(a<=0 || 1+(gamma_-1)*mu/bulk_<=0){return false;}
    ContourFlowSample candidate;candidate.area=a;candidate.relative_density=std::exp(std::log1p((gamma_-1)*mu/bulk_)/(gamma_-1));
    const bool right=z>.5*(state.support[0]+state.support[1]);ContourPressureState partial;
    if(!storage_.interval(state.area,state.effort,right?z:0.,right?end:z,partial)){return false;}
    const double sign=right?1.:-1.;const std::array<double,2> c={partial.compliance[0]+partial.compliance[2],partial.compliance[1]+partial.compliance[3]},d={partial.mass_area[0]+partial.mass_area[2],partial.mass_area[1]+partial.mass_area[3]};
    for(std::size_t j=0;j<5;++j){candidate.flow[j]=(j==(right?2u:0u)?1.:0.)+sign*(c[0]*state.effort_rate[j]+c[1]*state.effort_rate[5+j]+(j>=3?d[j-3]:0.));if(!std::isfinite(candidate.flow[j])){return false;}}
    for(std::size_t j:{0u,2u}){if(state.area[j==0?0:1]<=0){candidate.flow[j]=0.;}}
    if(!std::isfinite(candidate.relative_density)){return false;}result=candidate;return true;
}
bool ContourPressureKinetics::integrate(const ContourPressureFlow& state,double lo,double hi,std::array<double,25>& result,int depth) const noexcept {
    const double middle=.5*(lo+hi),radius=.5*(hi-lo);std::array<double,25> low{},high{},absolute{};
    for(int rule=0;rule<2;++rule){const int count=rule==0?4:8;const double* nodes=rule==0?points8:points16;const double* weights=rule==0?weights8:weights16;auto& sum=rule==0?low:high;
        for(int n=0;n<count;++n){for(double side:{-1.,1.}){ContourFlowSample point;if(!sample(state,middle+side*radius*nodes[n],point)){return false;}
            const double factor=radius*weights[n]*density_/(point.relative_density*point.area);
            for(std::size_t j=0;j<5;++j){for(std::size_t k=0;k<5;++k){const double term=factor*point.flow[j]*point.flow[k];sum[5*j+k]+=term;if(rule==1){absolute[5*j+k]+=std::abs(term);}}}}}
    }
    // The wall-motion columns can vanish exactly in a uniform tube. Bound
    // roundoff with the other diagonal energies after dimensional scaling,
    // rather than demanding relative accuracy of a numerically zero entry.
    const double length=thickness_[0]+thickness_[1];double characteristic=0;
    for(std::size_t j=0;j<5;++j){const double scale=j<3?1.:length;characteristic=std::max(characteristic,std::abs(high[5*j+j])/(scale*scale));}
    bool converged=true;for(std::size_t j=0;j<25;++j){if(!std::isfinite(high[j])){return false;}const double scale=(j/5<3?1.:length)*(j%5<3?1.:length);
        if(std::abs(high[j]-low[j])>2e-10*absolute[j]+2e-14*characteristic*scale+1e-290){converged=false;}}
    if(converged){result=high;return true;}if(depth==12){return false;}
    std::array<double,25> a,b;if(!integrate(state,lo,middle,a,depth+1) || !integrate(state,middle,hi,b,depth+1)){return false;}
    for(std::size_t j=0;j<25;++j){result[j]=a[j]+b[j];}return true;
}
bool ContourPressureKinetics::prepare(const std::array<double,2>& area,const std::array<double,2>& effort,ContourPressureFlow& result) const noexcept {
    ContourPressureState gas;if(!storage_.evaluate(area,effort,gas)){return false;}
    const auto& c=gas.compliance;const double det=c[0]*c[3]-c[1]*c[2];if(det<=0 || !std::isfinite(det)){return false;}
    ContourPressureFlow candidate;candidate.area=area;candidate.effort=effort;
    for(std::size_t j=0;j<5;++j){double a=j==0?1.:j==1?-1.:0.,b=j==1?1.:j==2?-1.:0.;if(j>=3){a=-gas.mass_area[j-3];b=-gas.mass_area[2+j-3];}
        candidate.effort_rate[j]=(c[3]*a-c[1]*b)/det;candidate.effort_rate[5+j]=(-c[2]*a+c[0]*b)/det;}
    const double c0=.5*thickness_[0],c1=thickness_[0]+.5*thickness_[1],end=thickness_[0]+thickness_[1];candidate.support={0.,end};
    if(area[0]<=0){candidate.support[0]=c0+(c1-c0)*(-area[0])/(area[1]-area[0]);}
    if(area[1]<=0){candidate.support[1]=c0+(c1-c0)*area[0]/(area[0]-area[1]);}
    // Closed exterior flows are constrained to zero. Their formal basis
    // columns would carry finite flow into a zero-area endpoint and have
    // infinite energy; remove them from the admissible kinetic metric.
    for(std::size_t j:{0u,2u}){if(area[j==0?0:1]<=0){candidate.effort_rate[j]=candidate.effort_rate[5+j]=0.;}}
    if(area[0]<=0 || area[1]<=0){const double t=-area[0]/(area[1]-area[0]),slope=std::abs(area[1]-area[0])/(c1-c0),factor=12*viscosity_*width_*width_/(slope*slope*slope);candidate.contact_loss={factor*(1-t)*(1-t),factor*(1-t)*t,factor*(1-t)*t,factor*t*t};}
    const double nodes[]={candidate.support[0],std::clamp(c0,candidate.support[0],candidate.support[1]),std::clamp(c1,candidate.support[0],candidate.support[1]),candidate.support[1]};
    for(std::size_t n=0;n<3;++n){if(nodes[n+1]<=nodes[n]){continue;}std::array<double,25> part;if(!integrate(candidate,nodes[n],nodes[n+1],part)){return false;}
        for(std::size_t j=0;j<25;++j){candidate.inertia[j]+=part[j];}}
    for(std::size_t j:{0u,2u}){if(area[j==0?0:1]<=0){for(std::size_t k=0;k<5;++k){candidate.inertia[5*j+k]=candidate.inertia[5*k+j]=0.;}}}
    for(double v:candidate.contact_loss){if(!std::isfinite(v)){return false;}}
    for(double v:candidate.inertia){if(!std::isfinite(v)){return false;}}
    result=candidate;return true;
}
bool ContourPressureKinetics::energy(const ContourPressureFlow& state,const std::array<double,5>& rate,double& result) const noexcept {
    if((state.area[0]<=0 && rate[0]!=0) || (state.area[1]<=0 && rate[2]!=0)){return false;}
    double value=0;for(std::size_t j=0;j<5;++j){if(!std::isfinite(rate[j])){return false;}for(std::size_t k=0;k<5;++k){value+=.5*rate[j]*state.inertia[5*j+k]*rate[k];}}
    if(!std::isfinite(value)){return false;}result=value;return true;
}
}
