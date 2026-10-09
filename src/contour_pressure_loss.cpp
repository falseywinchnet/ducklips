#include "ducklips/contour_pressure_loss.hpp"
#include "contour_pressure_intervals.hpp"
#include <algorithm>
#include <cmath>
namespace ducklips {
namespace {
#include "slit_bgk_table.inc"
constexpr double points8[]={.18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
constexpr double weights8[]={.36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
constexpr double points16[]={.09501250983763744019,.28160355077925891323,.45801677765722738634,.61787624440264374845,.75540440835500303390,.86563120238783174388,.94457502307323257608,.98940093499164993260};
constexpr double weights16[]={.18945061045506849629,.18260341504492358887,.16915651939500253819,.14959598881657673208,.12462897125553387205,.09515851168249278481,.06225352393864789286,.02715245941175409485};
}
bool slit_bgk_conductance(double delta,double& result) noexcept {
    if(!std::isfinite(delta) || delta<=0){return false;}
    double value;
    if(delta<=1e-8){value=(-std::log(delta)+.6341765026477007)/1.772453850905516;}
    else if(delta>=100.){const double inverse=1/delta;value=delta/6+1.0162+inverse*(1.0653-2.1354*inverse);}
    else {const double coordinate=(std::log(delta)-bgk_log_min)/bgk_log_step;const int index=std::clamp(static_cast<int>(coordinate),0,199);const double t=coordinate-index;const double* c=bgk_coefficients[index];value=std::exp(((c[0]*t+c[1])*t+c[2])*t+c[3]);}
    if(!std::isfinite(value) || value<=0){return false;}result=value;return true;
}
bool ContourPressureDissipation::configure(const FoldParameters& f,const AirProperties& air,ContourPressureProfile profile) noexcept {
    ContourPressureKinetics kinetics;if(!kinetics.configure(f,air,profile)){return false;}
    kinetics_=kinetics;thickness_={f.lower_thickness,f.upper_thickness};basis_=profile==ContourPressureProfile::linear_span?std::array<double,2>{0.,thickness_[0]+thickness_[1]}:std::array<double,2>{.5*thickness_[0],thickness_[0]+.5*thickness_[1]};density_=air.density;gamma_=air.heat_capacity_ratio;thermal_speed_=air.sound_speed*std::sqrt(2/gamma_);viscosity_=air.viscosity;width_=f.length;return true;
}
bool ContourPressureDissipation::integrate(const ContourPressureFlow& state,double start,double end,int mapping,double lo,double hi,std::array<double,25>& result,const std::array<double,25>* budget,int depth) const noexcept {
    const double middle=.5*(lo+hi),radius=.5*(hi-lo),span=end-start;std::array<double,25> low{},high{},absolute{};
    for(int rule=0;rule<2;++rule){const int count=rule==0?4:8;const double* nodes=rule==0?points8:points16;const double* weights=rule==0?weights8:weights16;auto& sum=rule==0?low:high;
        for(int n=0;n<count;++n){for(double side:{-1.,1.}){
            const double t=middle+side*radius*nodes[n];const double z=mapping==0?start+span*t:mapping==1?start+span*t*t*t:end-span*t*t*t;
            const double jacobian=mapping==0?span:3*span*t*t;ContourFlowSample point;if(!kinetics_.sample(state,z,point)){return false;}
            const double thermal=thermal_speed_*std::pow(point.relative_density,.5*(gamma_-1));
            const double delta=(point.area/width_)*density_*point.relative_density*thermal/(2*viscosity_);double g;
            if(!slit_bgk_conductance(delta,g)){return false;}
            // Factor F/A before squaring. Admissible flow vanishes at a contact
            // tip, so F/A remains bounded although resistance per flow diverges.
            const double factor=radius*weights[n]*jacobian*density_*thermal*width_/(point.relative_density*g);
            for(std::size_t j=0;j<5;++j){for(std::size_t k=0;k<5;++k){const double term=factor*(point.flow[j]/point.area)*(point.flow[k]/point.area);sum[5*j+k]+=term;if(rule==1){absolute[5*j+k]+=std::abs(term);}}}
        }}
    }
    // Allocate absolute error from the entire segment to each subinterval.
    // A purely local relative test cannot converge at a logarithmic endpoint.
    const auto& whole=budget?*budget:absolute;
    const double length=thickness_[0]+thickness_[1];double characteristic=0;
    for(std::size_t j=0;j<5;++j){const double scale=j<3?1.:length;characteristic=std::max(characteristic,whole[5*j+j]/(scale*scale));}
    bool converged=true;for(std::size_t j=0;j<25;++j){if(!std::isfinite(high[j])){return false;}const double scale=(j/5<3?1.:length)*(j%5<3?1.:length);
        if(std::abs(high[j]-low[j])>2e-8*std::max(absolute[j],whole[j]*(hi-lo))+2e-12*characteristic*scale*(hi-lo)+1e-290){converged=false;}}
    if(converged){result=high;return true;}if(depth==16){return false;}
    std::array<double,25> a,b;if(!integrate(state,start,end,mapping,lo,middle,a,&whole,depth+1) || !integrate(state,start,end,mapping,middle,hi,b,&whole,depth+1)){return false;}
    for(std::size_t j=0;j<25;++j){result[j]=a[j]+b[j];}return true;
}
bool ContourPressureDissipation::prepare(const ContourPressureFlow& state,ContourPressureLoss& result) const noexcept {
    if(basis_[0]==0 && state.area[1]<state.area[0]){ContourPressureLoss reflected;if(!prepare(pressure_reflect_flow(state,basis_[1]),reflected)){return false;}result.resistance=pressure_reflect_matrix(reflected.resistance);return true;}
    if(density_<=0 || state.support[1]<=state.support[0]){return false;}ContourPressureLoss candidate;
    if(viscosity_==0){result=candidate;return true;}
    const double c0=basis_[0],c1=basis_[1];
    const double nodes[]={state.support[0],std::clamp(c0,state.support[0],state.support[1]),std::clamp(c1,state.support[0],state.support[1]),state.support[1]};
    for(std::size_t n=0;n<3;++n){if(nodes[n+1]<=nodes[n]){continue;}
        const int mapping=nodes[n]==state.support[0] && state.area[0]<=0?1:nodes[n+1]==state.support[1] && state.area[1]<=0?2:0;
        std::array<double,514> knots;const std::size_t count=n==1?pressure_interval_knots(nodes[n],nodes[n+1],state.area[0],state.area[1],knots):pressure_interval_knots(nodes[n],nodes[n+1],1.,1.,knots);if(count==0){return false;}
        for(std::size_t k=0;k+1<count;++k){std::array<double,25> part;if(!integrate(state,knots[k],knots[k+1],mapping,0.,1.,part)){return false;}for(std::size_t j=0;j<25;++j){candidate.resistance[j]+=part[j];}}}
    result=candidate;return true;
}
bool ContourPressureDissipation::power(const ContourPressureFlow& state,const ContourPressureLoss& loss,const std::array<double,5>& rate,double& result) const noexcept {
    if((state.area[0]<=0 && rate[0]!=0) || (state.area[1]<=0 && rate[2]!=0)){return false;}
    double value=0;for(std::size_t j=0;j<5;++j){if(!std::isfinite(rate[j])){return false;}for(std::size_t k=0;k<5;++k){value+=rate[j]*loss.resistance[5*j+k]*rate[k];}}
    if(!std::isfinite(value)){return false;}result=value;return true;
}
}
