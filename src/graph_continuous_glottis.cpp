#include "ducklips/graph_continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace ducklips {
namespace {
bool combine(std::vector<double>& free,const std::vector<double>& response,double flow) noexcept {
    for(std::size_t i=0;i<free.size();++i){free[i]+=flow*response[i];if(!std::isfinite(free[i])){return false;}}
    return true;
}
}
void GraphContinuousGlottis::prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
    const FoldParameters& folds,const AirProperties& air,double tension,const FoldState& initial,
    const std::array<double,2>& effort,const std::array<double,3>& flow,double outlet_loss,ContourPressureProfile profile,double outlet_inertance){
    if(source_node>=nodes){throw std::invalid_argument("continuous source outside graph");}
    GraphContinuousGlottis candidate;
    if(!candidate.source.prepare(folds,air,tension,initial,effort,flow,outlet_loss,profile,outlet_inertance)){
        throw std::invalid_argument("invalid continuous source");
    }
    candidate.network.prepare(nodes,branches);candidate.free_.prepare(nodes,branches);candidate.response_.prepare(nodes,branches);
    candidate.unit_flow_.assign(nodes,0.);candidate.unit_flow_[source_node]=1.;candidate.source_node_=source_node;
    *this=std::move(candidate);
}
bool GraphContinuousGlottis::advance(std::span<const double> c1,std::span<const double> c2,
    std::span<const double> g1,std::span<const double> g2,double dt,double reservoir_pressure,
    const std::array<double,4>& external_force,std::span<const double> nodal_flow,
    std::span<const double> series_pressure,ContinuousGlottalUpdate update,double outlet_inertance) noexcept {
    if(unit_flow_.empty() || network.charge.size()!=free_.charge.size() || network.momentum.size()!=free_.momentum.size()
        || network.average_pressure.size()!=free_.average_pressure.size() || network.average_flow.size()!=free_.average_flow.size()){return false;}
    std::copy(network.charge.begin(),network.charge.end(),free_.charge.begin());
    std::copy(network.momentum.begin(),network.momentum.end(),free_.momentum.begin());
    std::fill(response_.charge.begin(),response_.charge.end(),0.);
    std::fill(response_.momentum.begin(),response_.momentum.end(),0.);
    if(!free_.advance(c1,c2,g1,g2,dt,nodal_flow,series_pressure)
        || !response_.advance(c1,c2,g1,g2,dt,unit_flow_)){return false;}
    GraphContinuousGlottalStep receipt;
    receipt.free_pressure=free_.average_pressure[source_node_];receipt.load_slope=response_.average_pressure[source_node_];
    ContinuousGlottis candidate=source;
    if(!candidate.advance(dt,reservoir_pressure,receipt.free_pressure,receipt.load_slope,external_force,update,outlet_inertance)){return false;}
    receipt.source_geometry_work=candidate.last.geometry_work;
    receipt.flow=candidate.last.flow[2];
    if(!combine(free_.charge,response_.charge,receipt.flow) || !combine(free_.momentum,response_.momentum,receipt.flow)
        || !combine(free_.average_pressure,response_.average_pressure,receipt.flow)
        || !combine(free_.average_flow,response_.average_flow,receipt.flow)){return false;}
    receipt.pressure=free_.average_pressure[source_node_];receipt.source_port_work=candidate.last.port_work;
    receipt.graph_port_work=dt*receipt.flow*receipt.pressure;receipt.port_residual=receipt.graph_port_work-receipt.source_port_work;
    for(std::size_t i=0;i<nodal_flow.size();++i){receipt.graph_external_work+=dt*nodal_flow[i]*free_.average_pressure[i];}
    for(std::size_t i=0;i<series_pressure.size();++i){receipt.graph_external_work+=dt*series_pressure[i]*free_.average_flow[i];}
    if(!std::isfinite(receipt.graph_port_work) || !std::isfinite(receipt.port_residual) || !std::isfinite(receipt.graph_external_work)){return false;}
    std::copy(free_.charge.begin(),free_.charge.end(),network.charge.begin());
    std::copy(free_.momentum.begin(),free_.momentum.end(),network.momentum.begin());
    std::copy(free_.average_pressure.begin(),free_.average_pressure.end(),network.average_pressure.begin());
    std::copy(free_.average_flow.begin(),free_.average_flow.end(),network.average_flow.begin());
    source=candidate;last=receipt;return true;
}
}
