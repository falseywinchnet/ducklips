#include "ducklips/contour_pressure_storage.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace ducklips {
namespace {
using Values=std::array<double,13>;
constexpr double points8[]={.18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
constexpr double weights8[]={.36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
constexpr double points16[]={.09501250983763744019,.28160355077925891323,.45801677765722738634,.61787624440264374845,.75540440835500303390,.86563120238783174388,.94457502307323257608,.98940093499164993260};
constexpr double weights16[]={.18945061045506849629,.18260341504492358887,.16915651939500253819,.14959598881657673208,.12462897125553387205,.09515851168249278481,.06225352393864789286,.02715245941175409485};
struct Integrand final {
    const GasStorage& gas;
    const std::array<double,2>& area;
    const std::array<double,2>& effort;
    double c0,c1,bulk,gamma;
    bool sample(double u,double t,double a,Values& f) const noexcept {
        const double mu=u*effort[0]+t*effort[1],scaled=(gamma-1)*mu/bulk;
        if(a<0 || scaled<=-1 || !std::isfinite(scaled)){return false;}
        double r;GasEnergy e;if(!gas.evaluate_effort(mu,e,r)){return false;}const double slope=std::pow(r,2-gamma)/bulk;
        f={a,a*r*u,a*r*t,a*e.energy,a*e.pressure,-e.pressure*u,-e.pressure*t,
            a*slope*u*u,a*slope*u*t,a*slope*t*t,r*u*u,r*u*t,r*t*t};
        for(double v:f){if(!std::isfinite(v)){return false;}}return true;
    }
    // Interpolate endpoint areas in the integration coordinate. Re-evaluating
    // signed area from rounded axial coordinates loses relative accuracy in
    // a short interval next to contact and can defeat adaptive quadrature.
    template<std::size_t N> bool quadrature(double lo,double hi,double alo,double ahi,const double (&nodes)[N],const double (&weights)[N],Values& sum,Values& absolute) const noexcept {
        const double radius=.5*(hi-lo),tlo=std::clamp((lo-c0)/(c1-c0),0.,1.),thi=std::clamp((hi-c0)/(c1-c0),0.,1.),ulo=std::clamp((c1-lo)/(c1-c0),0.,1.),uhi=std::clamp((c1-hi)/(c1-c0),0.,1.);sum.fill(0.);absolute.fill(0.);
        for(std::size_t j=0;j<N;++j){for(double sign:{-1.,1.}){const double fraction=.5*(1+sign*nodes[j]);Values f;if(!sample((1-fraction)*ulo+fraction*uhi,(1-fraction)*tlo+fraction*thi,(1-fraction)*alo+fraction*ahi,f)){return false;}
            for(std::size_t k=0;k<f.size();++k){const double term=radius*weights[j]*f[k];sum[k]+=term;absolute[k]+=std::abs(term);}}}return true;
    }
    bool integrate(double lo,double hi,double alo,double ahi,Values& output,int depth=0) const noexcept {
        Values low,high,absolute,unused;if(!quadrature(lo,hi,alo,ahi,points8,weights8,low,unused) || !quadrature(lo,hi,alo,ahi,points16,weights16,high,absolute)){return false;}
        bool converged=true;for(std::size_t k=0;k<high.size();++k){if(std::abs(high[k]-low[k])>2e-12*absolute[k]+1e-300){converged=false;}}
        if(converged){output=high;return true;}if(depth==12){return false;}
        Values left,right;const double middle=.5*(lo+hi);
        if(!integrate(lo,middle,alo,.5*(alo+ahi),left,depth+1) || !integrate(middle,hi,.5*(alo+ahi),ahi,right,depth+1)){return false;}
        for(std::size_t k=0;k<output.size();++k){output[k]=left[k]+right[k];}return true;
    }
};
}
double ContourPressureStorage::effort_difference(const std::array<double,2>& first_area,
    const std::array<double,2>& first_effort,const std::array<double,2>& last_area,
    const std::array<double,2>& last_effort) noexcept {
    constexpr double infinity=std::numeric_limits<double>::infinity();
    for(const auto* values:{&first_area,&first_effort,&last_area,&last_effort}){
        for(double v:*values){if(!std::isfinite(v)){return infinity;}}
    }
    const double d0=first_effort[0]-last_effort[0],d1=first_effort[1]-last_effort[1];
    if(!std::isfinite(d0) || !std::isfinite(d1)){return infinity;}
    double result=0;
    for(const auto* area:{&first_area,&last_area}){
        if((*area)[0]<=0 && (*area)[1]<=0){return infinity;}
        const double scale=std::max(std::abs((*area)[0]),std::abs((*area)[1]));
        const double a=(*area)[0]/scale,b=(*area)[1]/scale;
        // The difference is affine in the hat coordinate. Its absolute
        // maximum on each occupied interval occurs at an endpoint. Form both
        // weights directly, avoiding 1-t cancellation in a small air pocket.
        if(a>0){result=std::max(result,std::abs(d0));}
        if(b>0){result=std::max(result,std::abs(d1));}
        if(a<=0 || b<=0){
            const double u=b/(b-a),t=-a/(b-a);
            result=std::max(result,std::abs(std::fma(u,d0,t*d1)));
        }
    }
    return result;
}
bool ContourPressureStorage::configure(const FoldParameters& folds,const AirProperties& air,ContourPressureProfile profile) noexcept {
    if(!std::isfinite(folds.lower_thickness) || folds.lower_thickness<=0 || !std::isfinite(folds.upper_thickness) || folds.upper_thickness<=0){return false;}
    if(profile!=ContourPressureProfile::flat_caps && profile!=ContourPressureProfile::linear_span){return false;}
    GasStorage gas;if(!gas.configure(air)){return false;}
    gas_=gas;thickness_={folds.lower_thickness,folds.upper_thickness};bulk_=air.density*air.sound_speed*air.sound_speed;gamma_=air.heat_capacity_ratio;basis_=profile==ContourPressureProfile::linear_span?std::array<double,2>{0.,thickness_[0]+thickness_[1]}:std::array<double,2>{.5*thickness_[0],thickness_[0]+.5*thickness_[1]};return true;
}
bool ContourPressureStorage::evaluate(const std::array<double,2>& area,const std::array<double,2>& effort,ContourPressureState& result) const noexcept {
    ContourPressureState candidate;if(!interval(area,effort,0.,thickness_[0]+thickness_[1],candidate) || candidate.volume<=0){return false;}result=candidate;return true;
}
bool ContourPressureStorage::interval(const std::array<double,2>& area,const std::array<double,2>& effort,double first,double last,ContourPressureState& result) const noexcept {
    if(!std::isfinite(first) || !std::isfinite(last) || first<0 || last<first || last>thickness_[0]+thickness_[1]){return false;}
    if(bulk_<=0){return false;}for(double v:area){if(!std::isfinite(v)){return false;}}for(double v:effort){if(!std::isfinite(v)){return false;}}
    const double c0=basis_[0],c1=basis_[1],end=thickness_[0]+thickness_[1];
    const double nodes[]={0.,c0,c1,end};const Integrand integrand{gas_,area,effort,c0,c1,bulk_,gamma_};Values sum{};
    const auto opening=[&](double z){
        if(z<=c0){return area[0];}if(z>=c1){return area[1];}
        if((area[0]<=0 && area[1]>0) || (area[1]<=0 && area[0]>0)){const double root=area[0]==0?c0:area[1]==0?c1:c0+(c1-c0)*(-area[0])/(area[1]-area[0]);return (z-root)*(area[1]-area[0])/(c1-c0);}
        const double t=(z-c0)/(c1-c0);return (1-t)*area[0]+t*area[1];};
    for(std::size_t j=0;j<3;++j){double lo=std::max(first,nodes[j]),hi=std::min(last,nodes[j+1]);if(hi<=lo){continue;}const double a=opening(lo),b=opening(hi);if(a<=0 && b<=0){continue;}
        if(a<=0){lo=area[0]==0?c0:area[1]==0?c1:c0+(c1-c0)*(-area[0])/(area[1]-area[0]);}else if(b<=0){hi=area[0]==0?c0:area[1]==0?c1:c0+(c1-c0)*(-area[0])/(area[1]-area[0]);}
        // Density must be positive on the entire occupied interval, including
        // its contact endpoint, rather than only at quadrature nodes.
        for(double z:{lo,hi}){const double t=std::clamp((z-c0)/(c1-c0),0.,1.);if(1+(gamma_-1)*((1-t)*effort[0]+t*effort[1])/bulk_<=0){return false;}}
        Values part;if(!integrand.integrate(lo,hi,std::max(0.,a),std::max(0.,b),part)){return false;}for(std::size_t k=0;k<sum.size();++k){sum[k]+=part[k];}
    }
    ContourPressureState candidate;candidate.effort=effort;candidate.volume=sum[0];candidate.mass={sum[1],sum[2]};candidate.energy=sum[3];candidate.coenergy=sum[4];candidate.area_effort={sum[5],sum[6]};candidate.compliance={sum[7],sum[8],sum[8],sum[9]};candidate.mass_area={sum[10],sum[11],sum[11],sum[12]};
    for(double v:sum){if(!std::isfinite(v)){return false;}}result=candidate;return true;
}
bool ContourPressureStorage::recover(const std::array<double,2>& area,const std::array<double,2>& mass,const std::array<double,2>& guess,ContourPressureState& result) const noexcept {
    for(double m:mass){if(!std::isfinite(m) || m<=0){return false;}}
    const double scale=mass[0]+mass[1];if(!std::isfinite(scale)){return false;}ContourPressureState state;if(!evaluate(area,guess,state)){return false;}
    for(int iteration=0;iteration<48;++iteration){const std::array<double,2> residual={state.mass[0]-mass[0],state.mass[1]-mass[1]};const double norm=std::max(std::abs(residual[0]),std::abs(residual[1]));
        if(norm<=2e-13*scale){result=state;return true;}
        const auto& c=state.compliance;const double determinant=c[0]*c[3]-c[1]*c[2];if(determinant<=0 || !std::isfinite(determinant)){return false;}
        const std::array<double,2> delta={(-c[3]*residual[0]+c[1]*residual[1])/determinant,(c[2]*residual[0]-c[0]*residual[1])/determinant};
        bool accepted=false;double fraction=1.;for(int line=0;line<48;++line){const std::array<double,2> next={state.effort[0]+fraction*delta[0],state.effort[1]+fraction*delta[1]};ContourPressureState candidate;
            if(evaluate(area,next,candidate) && std::max(std::abs(candidate.mass[0]-mass[0]),std::abs(candidate.mass[1]-mass[1]))<norm){state=candidate;accepted=true;break;}fraction*=.5;}
        if(!accepted){return false;}
    }return false;
}
}
