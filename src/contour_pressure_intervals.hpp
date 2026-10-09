#pragma once
#include "ducklips/contour_pressure_flow.hpp"
#include <array>
#include <algorithm>
namespace ducklips {
// Exact reversal of the full-span coordinate z -> H-z. Exterior and
// middle flows change sign; wall rates and enthalpy moments exchange ends.
inline std::array<double,25> pressure_reflect_matrix(const std::array<double,25>& matrix) noexcept {
    constexpr std::size_t permutation[]={2,1,0,4,3};constexpr double sign[]={-1,-1,-1,1,1};std::array<double,25> result{};
    for(std::size_t i=0;i<5;++i){for(std::size_t j=0;j<5;++j){result[5*i+j]=sign[i]*sign[j]*matrix[5*permutation[i]+permutation[j]];}}return result;
}
inline ContourPressureFlow pressure_reflect_flow(const ContourPressureFlow& state,double length) noexcept {
    ContourPressureFlow result;result.area={state.area[1],state.area[0]};result.effort={state.effort[1],state.effort[0]};result.support={0.,length};if(result.area[0]<=0){result.support[0]=result.area[0]==0?0.:length*(-result.area[0])/(result.area[1]-result.area[0]);}if(result.area[1]<=0){result.support[1]=result.area[1]==0?length:length*result.area[0]/(result.area[0]-result.area[1]);}
    constexpr std::size_t permutation[]={2,1,0,4,3};constexpr double sign[]={-1,-1,-1,1,1};
    for(std::size_t i=0;i<2;++i){for(std::size_t j=0;j<5;++j){result.effort_rate[5*i+j]=sign[j]*state.effort_rate[5*(1-i)+permutation[j]];}}
    result.inertia=pressure_reflect_matrix(state.inertia);result.contact_loss={state.contact_loss[3],state.contact_loss[2],state.contact_loss[1],state.contact_loss[0]};return result;
}
// Resolve a positive affine gap by equal ratios of area. This changes only
// quadrature partitioning, retaining the complete physical integration span.
inline std::size_t pressure_interval_knots(double lo,double hi,double a,double b,std::array<double,514>& knots) noexcept {
    knots[0]=lo;knots[1]=hi;if(a<=0 || b<=0 || std::max(a,b)<=4*std::min(a,b)){return 2;}
    const bool reverse=a>b;const double small=std::min(a,b),large=std::max(a,b),span=hi-lo;std::size_t count=1;knots[0]=reverse?hi:lo;
    for(double value=small*4;value<large;value*=4){if(count>=513){return 0;}const double fraction=(value-small)/(large-small),z=reverse?hi-span*fraction:lo+span*fraction;if(z==knots[count-1]){return 0;}knots[count++]=z;}
    knots[count++]=reverse?lo:hi;if(reverse){std::reverse(knots.begin(),knots.begin()+count);}return count;
}
}
