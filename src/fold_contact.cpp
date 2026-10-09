#include "ducklips/fold_contact.hpp"
#include <algorithm>
#include <cmath>
namespace ducklips {
bool FoldContact::configure(const FoldParameters& f,double tension) noexcept {
    for(double x:{f.length,f.lower_stiffness,f.upper_stiffness,f.asymmetry,tension}){if(!std::isfinite(x) || x<=0){return false;}}
    if(!std::isfinite(f.collision_ratio) || f.collision_ratio<0){return false;}FoldContact candidate;
    const std::array<double,2> stiffness={f.lower_stiffness,f.upper_stiffness};
    for(std::size_t i=0;i<2;++i){const double left=stiffness[i]*tension,right=left*f.asymmetry;
        const double pair=1/(1/left+1/right);candidate.area_stiffness_[i]=f.collision_ratio*pair/(f.length*f.length);
        if(!std::isfinite(candidate.area_stiffness_[i])){return false;}}
    candidate.configured_=true;*this=candidate;return true;
}
bool FoldContact::evaluate(const std::array<double,2>& area,FoldContactState& result) const noexcept {
    if(!configured_){return false;}FoldContactState candidate;
    for(std::size_t i=0;i<2;++i){if(!std::isfinite(area[i])){return false;}const double a=std::min(0.,area[i]);
        candidate.area_effort[i]=area_stiffness_[i]*a;candidate.energy+=.5*area_stiffness_[i]*a*a;}
    if(!std::isfinite(candidate.energy)){return false;}result=candidate;return true;
}
bool FoldContact::transfer(const std::array<double,2>& first,const std::array<double,2>& last,FoldContactTransfer& result) const noexcept {
    FoldContactState a,b;if(!evaluate(first,a) || !evaluate(last,b)){return false;}FoldContactTransfer candidate;candidate.energy_change=b.energy-a.energy;double work=0;
    for(std::size_t i=0;i<2;++i){const double u=std::min(0.,first[i]),v=std::min(0.,last[i]);
        const double slope=first[i]==last[i]?(first[i]<0?1.:0.):(v-u)/(last[i]-first[i]);
        candidate.area_effort[i]=.5*area_stiffness_[i]*(u+v)*slope;work+=candidate.area_effort[i]*(last[i]-first[i]);}
    candidate.chain_residual=candidate.energy_change-work;if(!std::isfinite(candidate.chain_residual)){return false;}result=candidate;return true;
}
}
