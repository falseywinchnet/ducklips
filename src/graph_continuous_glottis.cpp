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
void GraphContinuousGlottis::prepare_thermal(const AirProperties& air,std::span<const ThermalSite> sites){
    std::vector<SlitThermal> candidate(sites.size());std::vector<std::size_t> nodes;
    for(std::size_t i=0;i<sites.size();++i){const auto& s=sites[i];
        if(s.node>=network.charge.size() || std::find(nodes.begin(),nodes.end(),s.node)!=nodes.end()
            || !candidate[i].prepare(air,s.compliance,s.half_gap)){throw std::invalid_argument("invalid thermal site");}
        nodes.push_back(s.node);
    }
    auto trial=candidate;thermal=std::move(candidate);thermal_trial_=std::move(trial);thermal_nodes_=std::move(nodes);
}
bool GraphContinuousGlottis::advance(std::span<const double> c1,std::span<const double> c2,
    std::span<const double> g1,std::span<const double> g2,double dt,double reservoir_pressure,
    const std::array<double,4>& external_force,std::span<const double> nodal_flow,
    std::span<const double> series_pressure,ContinuousGlottalUpdate update,double outlet_inertance,const RadiationGeometry* radiation_geometry,
    std::span<const ThermalGeometry> thermal_geometry) noexcept {
    if(unit_flow_.empty() || network.charge.size()!=free_.charge.size() || network.momentum.size()!=free_.momentum.size()
        || network.average_pressure.size()!=free_.average_pressure.size() || network.average_flow.size()!=free_.average_flow.size()
        || thermal.size()!=thermal_nodes_.size() || thermal.size()!=thermal_trial_.size() || thermal.size()!=thermal_geometry.size()){return false;}
    std::copy(network.charge.begin(),network.charge.end(),free_.charge.begin());
    std::copy(network.momentum.begin(),network.momentum.end(),free_.momentum.begin());
    RadiationBoundary next_radiation=radiation;GraphContinuousGlottalStep receipt;
    std::copy(thermal.begin(),thermal.end(),thermal_trial_.begin());
    for(std::size_t i=0;i<thermal_trial_.size();++i){auto& t=thermal_trial_[i];
        if(!t.advance(free_.charge[thermal_nodes_[i]],dt/2)){return false;}
        receipt.thermal_dissipation+=t.last.dissipated_energy;receipt.thermal_energy_residual+=t.last.energy_residual;
        receipt.thermal_air_residual+=t.last.stored_air_residual;
    }
    if(radiation_geometry){
        const auto& r=*radiation_geometry;double work;
        if(r.branch>=free_.momentum.size() || !next_radiation.change_area(r.midpoint_area,work)
            || !next_radiation.advance(free_.momentum[r.branch],r.first_conductance,dt/2)){return false;}
        receipt.radiation_geometry_work=work;receipt.radiated_energy=next_radiation.last.radiated_energy;
        receipt.radiation_energy_residual=next_radiation.last.energy_residual;
    }
    std::fill(response_.charge.begin(),response_.charge.end(),0.);
    std::fill(response_.momentum.begin(),response_.momentum.end(),0.);
    if(!free_.advance(c1,c2,g1,g2,dt,nodal_flow,series_pressure)
        || !response_.advance(c1,c2,g1,g2,dt,unit_flow_)){return false;}
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
    if(radiation_geometry){
        const auto& r=*radiation_geometry;double work;
        if(!next_radiation.advance(free_.momentum[r.branch],r.last_conductance,dt/2)
            || !next_radiation.change_area(r.last_area,work)){return false;}
        receipt.radiation_geometry_work+=work;receipt.radiated_energy+=next_radiation.last.radiated_energy;
        receipt.radiation_energy_residual+=next_radiation.last.energy_residual;
    }
    for(std::size_t i=0;i<thermal_trial_.size();++i){auto& t=thermal_trial_[i];const auto& geometry=thermal_geometry[i];double work;
        if(!t.change_geometry(geometry.compliance,geometry.half_gap,work)
            || !t.advance(free_.charge[thermal_nodes_[i]],dt/2)){return false;}
        receipt.thermal_geometry_work+=work;receipt.thermal_dissipation+=t.last.dissipated_energy;
        receipt.thermal_energy_residual+=t.last.energy_residual;receipt.thermal_air_residual+=t.last.stored_air_residual;
    }
    std::copy(free_.charge.begin(),free_.charge.end(),network.charge.begin());
    std::copy(free_.momentum.begin(),free_.momentum.end(),network.momentum.begin());
    std::copy(free_.average_pressure.begin(),free_.average_pressure.end(),network.average_pressure.begin());
    std::copy(free_.average_flow.begin(),free_.average_flow.end(),network.average_flow.begin());
    source=candidate;if(radiation_geometry){radiation=next_radiation;}
    std::copy(thermal_trial_.begin(),thermal_trial_.end(),thermal.begin());last=receipt;return true;
}
}
