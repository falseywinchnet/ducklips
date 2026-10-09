#include "ducklips/graph_orifice.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
void GraphOrifice::prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node){
    if(source_node>=nodes){throw std::invalid_argument("orifice source node outside graph");}
    state.prepare(nodes,branches);free_.prepare(nodes,branches);response_.prepare(nodes,branches);
    unit_flow_.assign(nodes,0.);unit_flow_[source_node]=1.;source_node_=source_node;
    momentum=0;flow=0;load_pressure=0;load_impedance=0;
}
bool GraphOrifice::advance(std::span<const double> c1,std::span<const double> c2,
    std::span<const double> g1,std::span<const double> g2,double dt,
    double area,double reservoir_pressure,double inlet_inertance,const AirProperties& air,
    std::span<const double> nodal_flow,std::span<const double> series_pressure) noexcept {
    if(unit_flow_.empty() || !std::isfinite(area) || area<0 || !std::isfinite(reservoir_pressure)
        || !std::isfinite(inlet_inertance) || inlet_inertance<0 || !std::isfinite(momentum)
        || !std::isfinite(air.density) || air.density<=0
        || state.charge.size()!=free_.charge.size() || state.momentum.size()!=free_.momentum.size()
        || state.average_pressure.size()!=free_.average_pressure.size() || state.average_flow.size()!=free_.average_flow.size()){return false;}
    std::copy(state.charge.begin(),state.charge.end(),free_.charge.begin());
    std::copy(state.momentum.begin(),state.momentum.end(),free_.momentum.begin());
    std::fill(response_.charge.begin(),response_.charge.end(),0.);
    std::fill(response_.momentum.begin(),response_.momentum.end(),0.);
    if(!free_.advance(c1,c2,g1,g2,dt,nodal_flow,series_pressure)
        || !response_.advance(c1,c2,g1,g2,dt,unit_flow_)){return false;}
    const double impedance=response_.average_pressure[source_node_];
    // A negative driving-point response is not a passive held-flow load.
    // Reject rather than silently choosing a Bernoulli root for that case.
    if(!std::isfinite(impedance) || impedance<0){return false;}
    const double history=inlet_inertance>0?2*momentum/dt:0;
    const double inertive_impedance=2*inlet_inertance/dt;
    const double effective_drop=reservoir_pressure-free_.average_pressure[source_node_]+history;
    const double effective_impedance=impedance+inertive_impedance;
    if(!std::isfinite(effective_drop) || !std::isfinite(effective_impedance)){return false;}
    const double incoming=orifice_flow(area,effective_drop,effective_impedance,air);
    const double next_momentum=(inlet_inertance>0 && area>0)?2*inlet_inertance*incoming-momentum:0;
    if(!std::isfinite(incoming) || !std::isfinite(next_momentum)){return false;}
    if(area>0){
        const double pressure=free_.average_pressure[source_node_]+impedance*incoming;
        const double velocity=incoming/area,loss=.5*air.density*velocity*std::abs(velocity);
        const double inertive_drop=inlet_inertance>0?(next_momentum-momentum)/dt:0;
        const double residual=reservoir_pressure-pressure-loss-inertive_drop;
        const double scale=std::abs(reservoir_pressure)+std::abs(pressure)+std::abs(loss)+std::abs(inertive_drop);
        if(!std::isfinite(residual) || !std::isfinite(scale) || std::abs(residual)>1e-10*scale){return false;}
    }
    const auto combine=[incoming](std::vector<double>& physical,const std::vector<double>& response){
        for(std::size_t i=0;i<physical.size();++i){physical[i]+=incoming*response[i];if(!std::isfinite(physical[i])){return false;}}
        return true;
    };
    if(!combine(free_.charge,response_.charge) || !combine(free_.momentum,response_.momentum)
        || !combine(free_.average_pressure,response_.average_pressure) || !combine(free_.average_flow,response_.average_flow)){return false;}
    std::copy(free_.charge.begin(),free_.charge.end(),state.charge.begin());
    std::copy(free_.momentum.begin(),free_.momentum.end(),state.momentum.begin());
    std::copy(free_.average_pressure.begin(),free_.average_pressure.end(),state.average_pressure.begin());
    std::copy(free_.average_flow.begin(),free_.average_flow.end(),state.average_flow.begin());
    momentum=next_momentum;flow=incoming;load_pressure=state.average_pressure[source_node_];load_impedance=impedance;
    return true;
}
}
