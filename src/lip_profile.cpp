#include "ducklips/lip_profile.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace ducklips {
namespace {
constexpr double tail=3.0348542587702927017; // sqrt(-log(1e-4)), existing slot cutoff
constexpr std::array<double,4> gx{.18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
constexpr std::array<double,4> gw{.36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
constexpr std::array<double,2> hx{.33998104358485626480,.86113631159405257522};
constexpr std::array<double,2> hw{.65214515486254614263,.34785484513745385737};
double integrand(const LipProfile& profile,double x,int kind) noexcept {
    const double area=profile.area(x);
    return kind==2?1/area:area*(kind==1?x:1);
}
bool integral(const LipProfile& p,double a,double b,int kind,unsigned depth,double& result) noexcept {
    const double middle=a+(b-a)*.5,half=(b-a)*.5;
    double coarse=0,fine=0;
    for(std::size_t i=0;i<gx.size();++i){fine+=gw[i]*(integrand(p,middle-half*gx[i],kind)+integrand(p,middle+half*gx[i],kind));}
    for(std::size_t i=0;i<hx.size();++i){coarse+=hw[i]*(integrand(p,middle-half*hx[i],kind)+integrand(p,middle+half*hx[i],kind));}
    fine*=half;coarse*=half;
    if(!std::isfinite(fine) || !std::isfinite(coarse)){return false;}
    if(std::abs(fine-coarse)<=2e-12*std::abs(fine)){result=fine;return true;}
    if(depth==36 || middle==a || middle==b){return false;}
    double left=0,right=0;
    if(!integral(p,a,middle,kind,depth+1,left) || !integral(p,middle,b,kind,depth+1,right)){return false;}
    result=left+right;return std::isfinite(result);
}
bool pieces(const LipProfile& p,double a,double b,int kind,double& result) noexcept {
    std::array<double,6> cuts{a,b,p.place-p.core,p.place+p.core,p.place-p.core-tail*p.width,p.place+p.core+tail*p.width};
    for(double& c:cuts){c=std::clamp(c,a,b);}std::sort(cuts.begin(),cuts.end());result=0;
    for(std::size_t i=1;i<cuts.size();++i){
        const double first=cuts[i-1],last=cuts[i];if(first==last){continue;}double value=0;
        const double middle=first+(last-first)*.5,distance=std::abs(middle-p.place)-p.core;
        if(distance<=0 || distance>tail*p.width){
            const double area=distance<=0?p.minimum_area():p.base_area;
            value=(last-first)*(kind==2?1/area:area*(kind==1?middle:1));
        }else{
            // Integrate in local shoulder coordinates. Evaluating x-core in
            // world coordinates loses relative precision in a tiny opening,
            // making a valid inverse-area integral fail its error check.
            const double sign=middle>=p.place?1.:-1.,origin=p.place+sign*p.core;
            double lo=sign*(first-origin)/p.width,hi=sign*(last-origin)/p.width;
            if(lo>hi){std::swap(lo,hi);}
            LipProfile local=p;local.place=0;local.core=0;local.width=1;
            if(kind==1){
                double area_integral=0,first_moment=0;
                if(!integral(local,lo,hi,0,0,area_integral) || !integral(local,lo,hi,1,0,first_moment)){return false;}
                value=p.width*(origin*area_integral+sign*p.width*first_moment);
            }else{
                if(!integral(local,lo,hi,kind,0,value)){return false;}value*=p.width;
            }
        }
        result+=value;
    }
    return std::isfinite(result);
}
bool bounds(const LipProfile& p,double a,double b) noexcept {
    return p.valid() && std::isfinite(a) && std::isfinite(b) && a>=0 && b>a && b<=p.length;
}
}
bool LipProfile::valid() const noexcept {
    for(double v:{length,base_area,activation,place,core,width,target_area}){if(!std::isfinite(v)){return false;}}
    return length>0 && base_area>=0 && core>=0 && width>0 && target_area>=0;
}
double LipProfile::minimum_area() const noexcept {
    const double a=std::clamp(activation,0.,1.),remaining=1-a,target=std::min(base_area,target_area);
    // 1-smoothstep(a) = (1-a)^2(1+2a). Form the small opening directly;
    // subtracting a rounded smoothstep from one can invent premature contact.
    return target+(base_area-target)*remaining*remaining*(1+2*a);
}
double LipProfile::area(double x) const noexcept {
    const double minimum=minimum_area(),distance=std::abs(x-place)-core;
    if(distance<=0){return minimum;}
    const double u=distance/width;
    if(u>tail){return base_area;}
    // Stable at contact: subtracting exp(-u*u) from one would erase tiny
    // positive shoulder openings and could turn them into false contact.
    return minimum+(base_area-minimum)*(-std::expm1(-u*u));
}
bool LipProfile::measure(double a,double b,double& volume,double& centroid) const noexcept {
    if(!bounds(*this,a,b)){return false;}
    double v=0,moment=0;
    if(!pieces(*this,a,b,0,v) || !pieces(*this,a,b,1,moment)){return false;}
    volume=v;centroid=v>0?moment/v:a+(b-a)*.5;
    return std::isfinite(centroid) && centroid>=a && centroid<=b;
}
bool LipProfile::conductance(double a,double b,double density,double& result) const noexcept {
    if(!bounds(*this,a,b) || !std::isfinite(density) || density<=0){return false;}
    if(base_area==0 || (minimum_area()==0 && std::max(a,place-core)<=std::min(b,place+core))){result=0;return true;}
    double inverse=0;
    if(!pieces(*this,a,b,2,inverse) || inverse<=0){return false;}
    const double g=1/(density*inverse);if(!std::isfinite(g) || g<=0){return false;}result=g;return true;
}
}
