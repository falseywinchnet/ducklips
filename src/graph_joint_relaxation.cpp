#include "ducklips/graph_continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace ducklips {
namespace {
bool partition(std::span<const std::size_t> owners,std::span<const ViscousGeometry> geometry){
    for(std::size_t i=0;i<owners.size();++i){const auto& x=geometry[i];double sum=0;
        if(!std::isfinite(x.segment_inertance) || x.segment_inertance<0 || !std::isfinite(x.branch_conductance) || x.branch_conductance<0){return false;}
        for(std::size_t j=0;j<owners.size();++j){if(owners[i]==owners[j]){if(x.branch_conductance!=geometry[j].branch_conductance){return false;}sum+=geometry[j].segment_inertance;}}
        if(!std::isfinite(sum) || sum*x.branch_conductance>1+1e-13){return false;}}
    return true;
}
bool metric_work(std::span<const double> q,std::span<const double> p,std::span<const double> ca,std::span<const double> cb,std::span<const double> ga,std::span<const double> gb,double& work){
    if(ca.size()!=q.size() || cb.size()!=q.size() || ga.size()!=p.size() || gb.size()!=p.size()){return false;}
    double value=0;
    for(std::size_t i=0;i<q.size();++i){if(!std::isfinite(ca[i]) || !std::isfinite(cb[i]) || ca[i]<=0 || cb[i]<=0){return false;}value+=.5*q[i]*q[i]*(1/cb[i]-1/ca[i]);}
    for(std::size_t i=0;i<p.size();++i){if(!std::isfinite(ga[i]) || !std::isfinite(gb[i]) || ga[i]<0 || gb[i]<0){return false;}value+=.5*p[i]*p[i]*(gb[i]-ga[i]);}
    if(!std::isfinite(value)){return false;}work=value;return true;
}
}
void GraphContinuousGlottis::prepare_joint(std::size_t radiation_branch,RelaxationOrder order){
    if(network.charge.empty() || (radiation_branch!=ground_node && radiation_branch>=joint_edges_.size())){throw std::invalid_argument("joint relaxation topology");}
    std::vector<std::size_t> no,bo;
    for(auto node:thermal_nodes_){for(std::size_t j=0;j<SlitThermal::modes;++j){no.push_back(node);}}
    for(auto branch:viscous_branches_){for(std::size_t j=0;j<SlitViscous::modes;++j){bo.push_back(branch);}}
    if(radiation_branch!=ground_node){bo.push_back(radiation_branch);}
    joint_free_.prepare(network.charge.size(),joint_edges_,no,bo,order);joint_response_.prepare(network.charge.size(),joint_edges_,no,bo,order);joint_final_.prepare(network.charge.size(),joint_edges_,no,bo,order);
    joint_node_metrics_.resize(no.size());joint_branch_metrics_.resize(bo.size());joint_resistance_.resize(joint_edges_.size());joint_drive_.resize(network.charge.size());
    joint_radiation_branch_=radiation_branch;joint_ready_=true;
}
bool GraphContinuousGlottis::advance_joint(const JointGraphGeometry& geometry,double dt,double reservoir,
    const std::array<double,4>& external,std::span<const double> nodal,std::span<const double> series,double outlet_inertance) noexcept {
    if(!joint_ready_ || network.charge.size()!=joint_free_.charge.size() || network.momentum.size()!=joint_free_.momentum.size()
        || network.average_pressure.size()!=network.charge.size() || network.average_flow.size()!=network.momentum.size()
        || geometry.last_compliance.size()!=network.charge.size() || geometry.last_conductance.size()!=network.momentum.size()
        || (!nodal.empty() && nodal.size()!=network.charge.size())
        || geometry.midpoint_thermal.size()!=thermal.size() || geometry.last_thermal.size()!=thermal.size()
        || geometry.midpoint_viscous.size()!=viscous.size() || geometry.last_viscous.size()!=viscous.size()
        || joint_node_metrics_.size()!=thermal.size()*SlitThermal::modes
        || joint_branch_metrics_.size()!=viscous.size()*SlitViscous::modes+(joint_radiation_branch_!=ground_node)
        || !std::isfinite(dt) || dt<=0){return false;}
    if(!partition(viscous_branches_,geometry.midpoint_viscous) || !partition(viscous_branches_,geometry.last_viscous)){return false;}
    GraphContinuousGlottalStep receipt;double work;
    if(!metric_work(network.charge,network.momentum,geometry.first_compliance,geometry.midpoint_compliance,geometry.first_conductance,geometry.midpoint_conductance,work)){return false;}
    receipt.graph_geometry_work=work;
    std::copy(thermal.begin(),thermal.end(),thermal_trial_.begin());std::copy(viscous.begin(),viscous.end(),viscous_trial_.begin());
    RadiationBoundary next_radiation=radiation;
    for(std::size_t i=0;i<thermal_trial_.size();++i){const auto& x=geometry.midpoint_thermal[i];
        if(!thermal_trial_[i].change_geometry(x.compliance,x.half_gap,work)){return false;}receipt.thermal_geometry_work+=work;
        for(std::size_t j=0;j<SlitThermal::modes;++j){joint_node_metrics_[i*SlitThermal::modes+j]={thermal_trial_[i].memory_storage(j),thermal_trial_[i].memory_rate(j)};}}
    std::fill(joint_resistance_.begin(),joint_resistance_.end(),0.);
    for(std::size_t i=0;i<viscous_trial_.size();++i){const auto& x=geometry.midpoint_viscous[i];const auto branch=viscous_branches_[i];
        if(x.branch_conductance!=geometry.midpoint_conductance[branch] || !viscous_trial_[i].change_geometry(x.segment_inertance,x.branch_conductance,x.half_gap,work)){return false;}
        receipt.viscous_geometry_work+=work;joint_resistance_[branch]+=viscous_trial_[i].steady_resistance();
        for(std::size_t j=0;j<SlitViscous::modes;++j){joint_branch_metrics_[i*SlitViscous::modes+j]={viscous_trial_[i].memory_storage(j),viscous_trial_[i].memory_rate(j)};}}
    if(joint_radiation_branch_!=ground_node){
        if(!next_radiation.change_area(geometry.midpoint_radiation_area,work)){return false;}receipt.radiation_geometry_work+=work;
        const double h=next_radiation.inverse_inertance();joint_branch_metrics_.back()=h>0?RelaxationMetric{1/h,next_radiation.resistance()*h}:RelaxationMetric{};}
    const auto load=[&](RelaxationGraphReference& graph){
        graph.charge=network.charge;graph.momentum=network.momentum;
        for(std::size_t i=0;i<thermal_trial_.size();++i){for(std::size_t j=0;j<SlitThermal::modes;++j){graph.node_memory[i*SlitThermal::modes+j]=thermal_trial_[i].charge[j];}}
        for(std::size_t i=0;i<viscous_trial_.size();++i){for(std::size_t j=0;j<SlitViscous::modes;++j){graph.branch_memory[i*SlitViscous::modes+j]=viscous_trial_[i].momentum[j];}}
        if(joint_radiation_branch_!=ground_node){graph.branch_memory.back()=next_radiation.momentum;}
    };
    load(joint_free_);
    if(!joint_free_.advance(geometry.midpoint_compliance,geometry.midpoint_conductance,joint_resistance_,joint_node_metrics_,joint_branch_metrics_,dt,nodal,series)){return false;}
    for(auto* v:{&joint_response_.charge,&joint_response_.momentum,&joint_response_.node_memory,&joint_response_.branch_memory}){std::fill(v->begin(),v->end(),0.);}
    if(!joint_response_.advance(geometry.midpoint_compliance,geometry.midpoint_conductance,joint_resistance_,joint_node_metrics_,joint_branch_metrics_,dt,unit_flow_)){return false;}
    receipt.free_pressure=joint_free_.average_pressure[source_node_];receipt.load_slope=joint_response_.average_pressure[source_node_];
    ContinuousGlottis candidate=source;
    if(!candidate.advance(dt,reservoir,receipt.free_pressure,receipt.load_slope,external,ContinuousGlottalUpdate::relaxation,outlet_inertance)){return false;}
    receipt.flow=candidate.last.flow[2];receipt.source_port_work=candidate.last.port_work;receipt.source_geometry_work=candidate.last.geometry_work;
    for(std::size_t i=0;i<joint_drive_.size();++i){joint_drive_[i]=nodal.empty()?0:nodal[i];}joint_drive_[source_node_]+=receipt.flow;
    // Re-evaluate the accepted linear combination so its quadratic work and
    // dissipation ledger includes the source response cross terms exactly.
    load(joint_final_);
    if(!joint_final_.advance(geometry.midpoint_compliance,geometry.midpoint_conductance,joint_resistance_,joint_node_metrics_,joint_branch_metrics_,dt,joint_drive_,series)){return false;}
    receipt.pressure=joint_final_.average_pressure[source_node_];receipt.graph_port_work=dt*receipt.flow*receipt.pressure;receipt.port_residual=receipt.graph_port_work-receipt.source_port_work;
    receipt.graph_external_work=joint_final_.last.external_work-receipt.graph_port_work;
    receipt.passive_numerical_dissipation=joint_final_.last.numerical_dissipation;receipt.passive_energy_residual=joint_final_.last.energy_residual;receipt.passive_air_residual=joint_final_.last.air_residual;
    for(std::size_t i=0;i<thermal_trial_.size();++i){
        for(std::size_t j=0;j<SlitThermal::modes;++j){const auto k=i*SlitThermal::modes+j;thermal_trial_[i].charge[j]=joint_final_.node_memory[k];receipt.thermal_dissipation+=joint_final_.node_dissipation[k];}
        const auto& x=geometry.last_thermal[i];if(!thermal_trial_[i].change_geometry(x.compliance,x.half_gap,work)){return false;}receipt.thermal_geometry_work+=work;}
    for(std::size_t i=0;i<viscous_trial_.size();++i){
        for(std::size_t j=0;j<SlitViscous::modes;++j){const auto k=i*SlitViscous::modes+j;viscous_trial_[i].momentum[j]=joint_final_.branch_memory[k];receipt.viscous_dissipation+=joint_final_.branch_dissipation[k];}
        const auto& x=geometry.last_viscous[i];if(x.branch_conductance!=geometry.last_conductance[viscous_branches_[i]] || !viscous_trial_[i].change_geometry(x.segment_inertance,x.branch_conductance,x.half_gap,work)){return false;}receipt.viscous_geometry_work+=work;}
    for(double heat:joint_final_.resistance_dissipation){receipt.viscous_dissipation+=heat;}
    if(joint_radiation_branch_!=ground_node){next_radiation.momentum=joint_final_.branch_memory.back();receipt.radiated_energy=joint_final_.branch_dissipation.back();
        if(!next_radiation.change_area(geometry.last_radiation_area,work)){return false;}receipt.radiation_geometry_work+=work;}
    if(!metric_work(joint_final_.charge,joint_final_.momentum,geometry.midpoint_compliance,geometry.last_compliance,geometry.midpoint_conductance,geometry.last_conductance,work)){return false;}
    receipt.graph_geometry_work+=work;receipt.passive_energy_change=joint_final_.last.energy_change+receipt.graph_geometry_work+receipt.thermal_geometry_work+receipt.viscous_geometry_work+receipt.radiation_geometry_work;
    if(!std::isfinite(receipt.passive_energy_change) || !std::isfinite(receipt.port_residual)){return false;}
    network.charge=joint_final_.charge;network.momentum=joint_final_.momentum;network.average_pressure=joint_final_.average_pressure;network.average_flow=joint_final_.average_flow;
    for(std::size_t i=0;i<network.momentum.size();++i){if(geometry.last_conductance[i]==0){network.momentum[i]=0.;}}
    source=candidate;radiation=next_radiation;
    std::copy(thermal_trial_.begin(),thermal_trial_.end(),thermal.begin());std::copy(viscous_trial_.begin(),viscous_trial_.end(),viscous.begin());last=receipt;return true;
}
}
