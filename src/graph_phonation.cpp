#include "ducklips/graph_phonation.hpp"
#include <algorithm>
#include <cmath>

namespace ducklips {
namespace {
bool copy_state(const GraphOrifice& from,GraphOrifice& to) noexcept {
    if(from.state.charge.size()!=to.state.charge.size() || from.state.momentum.size()!=to.state.momentum.size()
        || from.state.average_pressure.size()!=to.state.average_pressure.size() || from.state.average_flow.size()!=to.state.average_flow.size()){return false;}
    std::copy(from.state.charge.begin(),from.state.charge.end(),to.state.charge.begin());
    std::copy(from.state.momentum.begin(),from.state.momentum.end(),to.state.momentum.begin());
    std::copy(from.state.average_pressure.begin(),from.state.average_pressure.end(),to.state.average_pressure.begin());
    std::copy(from.state.average_flow.begin(),from.state.average_flow.end(),to.state.average_flow.begin());
    to.momentum=from.momentum;to.flow=from.flow;to.load_pressure=from.load_pressure;to.load_impedance=from.load_impedance;
    return true;
}
bool finite_state(const TwoMassGlottis& folds) noexcept {
    TwoMassGlottis scratch=folds;return scratch.restore(folds.state());
}
}
void GraphPhonation::prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,const FoldParameters& parameters){
    source.prepare(nodes,branches,source_node);trial_.prepare(nodes,branches,source_node);folds.reset(parameters);
    midpoint_area=0;transglottal_pressure=0;coupling_residual=0;iterations=0;
}
bool GraphPhonation::advance(std::span<const double> c1,std::span<const double> c2,
    std::span<const double> g1,std::span<const double> g2,double dt,
    double rest_area,double tension,double reservoir_pressure,double inlet_inertance,
    double surface_pressure,const AirProperties& air,
    std::span<const double> nodal_flow,std::span<const double> series_pressure) noexcept {
    if(!std::isfinite(rest_area) || !std::isfinite(tension) || tension<=0 || !std::isfinite(surface_pressure)
        || !std::isfinite(dt) || dt<=0 || !finite_state(folds) || !copy_state(source,trial_)){return false;}
    // A zero-flow trial supplies the free network pressure and unit-flow
    // response. It changes scratch state only. Retain the accepted inlet
    // momentum for the real source solve below.
    if(!trial_.advance(c1,c2,g1,g2,dt,0,reservoir_pressure,inlet_inertance,air,nodal_flow,series_pressure)){return false;}
    const double free_drop=reservoir_pressure-trial_.load_pressure;
    const double history=inlet_inertance>0?2*source.momentum/dt:0;
    const double effective_drop=free_drop+history;
    const double impedance=trial_.load_impedance+2*inlet_inertance/dt;
    if(!std::isfinite(effective_drop) || !std::isfinite(impedance)){return false;}
    double low=std::min({0.,free_drop,effective_drop}),high=std::max({0.,free_drop,effective_drop});
    const double tolerance=1e-10*(1+std::max(std::abs(low),std::abs(high)));
    double drop=0,area=0,residual=0;std::size_t count=0;
    TwoMassGlottis middle;
    for(;count<80;++count){
        drop=.5*low+.5*high;middle=folds;
        middle.advance(drop,rest_area,tension,.5*dt,1,surface_pressure);
        area=std::max(middle.minimum_area(),0.);
        if(!finite_state(middle) || !std::isfinite(area)){return false;}
        // A sealed source removes its boundary-path momentum. The retained
        // downstream air still applies its free pressure to the closed folds.
        double actual=free_drop;
        if(area>0){
            const double flow=orifice_flow(area,effective_drop,impedance,air);
            const double velocity=flow/area;actual=.5*air.density*velocity*std::abs(velocity);
        }
        residual=drop-actual;if(!std::isfinite(residual)){return false;}
        if(std::abs(residual)<=tolerance){break;}
        if(residual>0){high=drop;}else{low=drop;}
    }
    if(count==80){return false;}
    TwoMassGlottis next=middle;next.advance(drop,rest_area,tension,.5*dt,1,surface_pressure);
    if(!finite_state(next) || !copy_state(source,trial_)){return false;}
    if(!trial_.advance(c1,c2,g1,g2,dt,area,reservoir_pressure,inlet_inertance,air,nodal_flow,series_pressure)){return false;}
    // Verify against the committed source solve, including the exact seal.
    const double velocity=area>0?trial_.flow/area:0;
    const double actual=area>0?.5*air.density*velocity*std::abs(velocity):reservoir_pressure-trial_.load_pressure;
    if(!std::isfinite(actual) || std::abs(drop-actual)>tolerance){return false;}
    if(!copy_state(trial_,source)){return false;}
    folds=next;midpoint_area=area;transglottal_pressure=drop;coupling_residual=drop-actual;iterations=count+1;
    return true;
}
}
