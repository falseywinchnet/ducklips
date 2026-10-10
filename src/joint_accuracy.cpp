#include "ducklips/joint_accuracy.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace ducklips {
namespace {
double scaled(double a,double b,double scale){return std::abs(a-b)/scale;}
constexpr double infinity=std::numeric_limits<double>::infinity();
}
void JointAccuracy::prepare(const GraphContinuousGlottis& model,std::size_t maximum_trials){
    if(!maximum_trials || maximum_trials>1000000 || model.network.charge.empty() || !model.joint_ready_){throw std::invalid_argument("joint accuracy preparation");}
    accepted_=full_=halves_=model;maximum_trials_=maximum_trials;
    const auto n=model.network.charge.size(),m=model.network.momentum.size();
    end_compliance_.resize(n);end_conductance_.resize(m);average_pressure.assign(n,0.);average_flow.assign(m,0.);
    pressure_sum_.resize(n);flow_sum_.resize(m);trial_pressure_.resize(n);trial_flow_.resize(m);
    end_times_.resize(2*maximum_trials);trial_times_.resize(2*maximum_trials);pending_ends_.resize(maximum_trials);last={};
}
bool JointAccuracy::compatible(const GraphContinuousGlottis& model) const noexcept {
    return model.joint_ready_ && model.source_node_==accepted_.source_node_
        && model.joint_radiation_branch_==accepted_.joint_radiation_branch_
        && model.joint_free_.order()==accepted_.joint_free_.order()
        && model.thermal_nodes_==accepted_.thermal_nodes_ && model.viscous_branches_==accepted_.viscous_branches_
        && model.joint_edges_.size()==accepted_.joint_edges_.size()
        && std::equal(model.joint_edges_.begin(),model.joint_edges_.end(),accepted_.joint_edges_.begin(),[](const TreeEdge& a,const TreeEdge& b){return a.first_node==b.first_node && a.second_node==b.second_node;})
        && model.network.charge.size()==average_pressure.size() && model.network.momentum.size()==average_flow.size()
        && model.network.charge.size()==accepted_.network.charge.size() && model.network.momentum.size()==accepted_.network.momentum.size()
        && model.thermal.size()==accepted_.thermal.size() && model.viscous.size()==accepted_.viscous.size();
}
double JointAccuracy::difference(const GraphContinuousGlottis& a,const GraphContinuousGlottis& b,const JointGraphGeometry& geometry) const noexcept {
    ContinuousGlottalState x,y;if(!a.source.state(x) || !b.source.state(y)){return infinity;}
    double error=0;
    for(std::size_t j=0;j<2;++j){
        error=std::max(error,scaled(x.effort[j],y.effort[j],1000.));error=std::max(error,scaled(x.area[j],y.area[j],1e-6));
        error=std::max(error,scaled(x.mass[j],y.mass[j],1e-8));
        error=std::max(error,scaled(x.folds.lower_position[j],y.folds.lower_position[j],1e-4));
        error=std::max(error,scaled(x.folds.upper_position[j],y.folds.upper_position[j],1e-4));
        error=std::max(error,scaled(x.folds.lower_velocity[j],y.folds.lower_velocity[j],.03571428571428571));
        error=std::max(error,scaled(x.folds.upper_velocity[j],y.folds.upper_velocity[j],.03571428571428571));
    }
    for(std::size_t j=0;j<3;++j){error=std::max(error,scaled(x.flow[j],y.flow[j],1e-4));}
    for(std::size_t j=0;j<a.network.charge.size();++j){error=std::max(error,scaled(a.network.charge[j],b.network.charge[j],1000.*geometry.last_compliance[j]));}
    for(std::size_t j=0;j<a.network.momentum.size();++j){error=std::max(error,std::abs(a.network.momentum[j]-b.network.momentum[j])*std::sqrt(geometry.last_conductance[j]/1e-6));}
    for(std::size_t i=0;i<a.thermal.size();++i){for(std::size_t j=0;j<SlitThermal::modes;++j){
        const double k=a.thermal[i].memory_storage(j);if(k!=b.thermal[i].memory_storage(j)){return infinity;}
        const double delta=std::abs(a.thermal[i].charge[j]-b.thermal[i].charge[j]);if(k==0){if(delta!=0){return infinity;}}else{error=std::max(error,delta/std::sqrt(k*1e-6));}
    }}
    for(std::size_t i=0;i<a.viscous.size();++i){for(std::size_t j=0;j<SlitViscous::modes;++j){
        const double k=a.viscous[i].memory_storage(j);if(k!=b.viscous[i].memory_storage(j)){return infinity;}
        const double delta=std::abs(a.viscous[i].momentum[j]-b.viscous[i].momentum[j]);
        if(k==0){if(delta!=0){return infinity;}}else{error=std::max(error,delta/std::sqrt(k*1e-6));}
    }}
    if(a.radiation.inverse_inertance()!=b.radiation.inverse_inertance()){return infinity;}
    error=std::max(error,std::abs(a.radiation.momentum-b.radiation.momentum)*std::sqrt(a.radiation.inverse_inertance()/1e-6));
    return std::isfinite(error)?error:infinity;
}
void JointAccuracy::accumulate(const GraphContinuousGlottis& model,double dt,JointAccuracyStep& receipt,
    std::vector<double>& pressure,std::vector<double>& flow) noexcept {
    const auto& s=model.source.last;const auto& g=model.last;
    receipt.energy_change+=s.energy_change+g.passive_energy_change;
    receipt.external_work+=s.reservoir_work+s.external_work+s.geometry_work+g.graph_external_work+g.graph_geometry_work+g.thermal_geometry_work+g.viscous_geometry_work+g.radiation_geometry_work;
    receipt.physical_loss+=s.fluid_loss+s.tissue_loss+s.outlet_loss+g.thermal_dissipation+g.viscous_dissipation+g.radiated_energy;
    receipt.numerical_loss+=s.numerical_loss+g.passive_numerical_dissipation;
    receipt.reservoir_volume+=dt*s.flow[0];receipt.source_port_work+=g.source_port_work;receipt.graph_port_work+=g.graph_port_work;
    for(std::size_t j=0;j<pressure.size();++j){pressure[j]+=dt*model.network.average_pressure[j];}
    for(std::size_t j=0;j<flow.size();++j){flow[j]+=dt*model.network.average_flow[j];}
}
bool JointAccuracy::advance(GraphContinuousGlottis& model,double first,double end,JointGeometryProvider provider,void* context,const JointAccuracyOptions& options) noexcept {
    failure=JointAccuracyFailure::configuration;attempt_time=first;attempt_step=end-first;attempt_error=infinity;
    if(!maximum_trials_ || !compatible(model) || !provider || !std::isfinite(first) || !std::isfinite(end) || end<=first
        || !std::isfinite(options.tolerance) || options.tolerance<=0 || !std::isfinite(options.maximum_step) || options.maximum_step<=0
        || !std::isfinite(options.minimum_step) || options.minimum_step<=0 || options.maximum_step<2*options.minimum_step
        || !options.maximum_trials || options.maximum_trials>maximum_trials_ || !std::isfinite(end-first)){return false;}
    accepted_=model;JointAccuracyStep receipt;receipt.minimum_step=infinity;
    std::fill(pressure_sum_.begin(),pressure_sum_.end(),0.);std::fill(flow_sum_.begin(),flow_sum_.end(),0.);
    const double block_count=std::ceil((end-first)/options.maximum_step);
    if(!std::isfinite(block_count) || block_count>static_cast<double>(options.maximum_trials)){return false;}
    const auto blocks=static_cast<std::size_t>(std::max(1.,block_count));
    // Preserve requested interval endpoints exactly. Repeated addition of a
    // small step can otherwise leave an unresolvable floating-point tail.
    for(std::size_t block=0;block<blocks;++block){
      double t=std::lerp(first,end,static_cast<double>(block)/static_cast<double>(blocks));
      double finish=std::lerp(first,end,static_cast<double>(block+1)/static_cast<double>(blocks));
      std::size_t pending=0;
      while(true){
        if(receipt.trials==options.maximum_trials){failure=JointAccuracyFailure::trial_limit;return false;}++receipt.trials;
        const double h=finish-t,middle=t+.5*h;
        attempt_time=t;attempt_step=h;
        if(middle==t || middle==finish || h/2<options.minimum_step){failure=JointAccuracyFailure::step_limit;return false;}
        full_=accepted_;halves_=accepted_;JointAccuracyStep trial;
        std::fill(trial_pressure_.begin(),trial_pressure_.end(),0.);std::fill(trial_flow_.begin(),trial_flow_.end(),0.);
        JointAdvanceInputs inputs;
        bool solved=provider(context,t,finish,inputs);
        if(solved){solved=inputs.geometry.last_compliance.size()==end_compliance_.size() && inputs.geometry.last_conductance.size()==end_conductance_.size();}
        if(solved){end_compliance_.assign(inputs.geometry.last_compliance.begin(),inputs.geometry.last_compliance.end());end_conductance_.assign(inputs.geometry.last_conductance.begin(),inputs.geometry.last_conductance.end());
            solved=full_.advance_joint(inputs.geometry,h,inputs.reservoir_pressure,inputs.external_force,inputs.nodal_flow,inputs.series_pressure,inputs.outlet_inertance);}
        if(solved){solved=provider(context,t,middle,inputs) && halves_.advance_joint(inputs.geometry,middle-t,inputs.reservoir_pressure,inputs.external_force,inputs.nodal_flow,inputs.series_pressure,inputs.outlet_inertance);}
        if(solved){accumulate(halves_,middle-t,trial,trial_pressure_,trial_flow_);
            solved=provider(context,middle,finish,inputs) && halves_.advance_joint(inputs.geometry,finish-middle,inputs.reservoir_pressure,inputs.external_force,inputs.nodal_flow,inputs.series_pressure,inputs.outlet_inertance);}
        double error=infinity;
        if(solved){
            solved=std::equal(end_compliance_.begin(),end_compliance_.end(),inputs.geometry.last_compliance.begin(),inputs.geometry.last_compliance.end())
                && std::equal(end_conductance_.begin(),end_conductance_.end(),inputs.geometry.last_conductance.begin(),inputs.geometry.last_conductance.end());
            if(solved){accumulate(halves_,finish-middle,trial,trial_pressure_,trial_flow_);error=difference(full_,halves_,inputs.geometry);}
        }
        attempt_error=error;
        if(solved && error<=options.tolerance){
            accepted_=halves_;receipt.energy_change+=trial.energy_change;receipt.external_work+=trial.external_work;receipt.physical_loss+=trial.physical_loss;receipt.numerical_loss+=trial.numerical_loss;
            receipt.reservoir_volume+=trial.reservoir_volume;receipt.source_port_work+=trial.source_port_work;receipt.graph_port_work+=trial.graph_port_work;
            receipt.maximum_accepted_error=std::max(receipt.maximum_accepted_error,error);receipt.minimum_step=std::min(receipt.minimum_step,std::min(middle-t,finish-middle));
            trial_times_[receipt.accepted_steps++]=middle;trial_times_[receipt.accepted_steps++]=finish;
            for(std::size_t j=0;j<pressure_sum_.size();++j){pressure_sum_[j]+=trial_pressure_[j];}for(std::size_t j=0;j<flow_sum_.size();++j){flow_sum_[j]+=trial_flow_[j];}
            t=finish;if(!pending){break;}finish=pending_ends_[--pending];
        }else{++receipt.rejected;pending_ends_[pending++]=finish;finish=middle;}
      }
    }
    receipt.energy_residual=receipt.energy_change+receipt.physical_loss+receipt.numerical_loss-receipt.external_work;
    if(!std::isfinite(receipt.energy_residual)){return false;}
    for(double p:pressure_sum_){if(!std::isfinite(p/(end-first))){return false;}}for(double q:flow_sum_){if(!std::isfinite(q/(end-first))){return false;}}
    model=accepted_;for(std::size_t j=0;j<average_pressure.size();++j){average_pressure[j]=pressure_sum_[j]/(end-first);}for(std::size_t j=0;j<average_flow.size();++j){average_flow[j]=flow_sum_[j]/(end-first);}
    std::copy_n(trial_times_.begin(),receipt.accepted_steps,end_times_.begin());last=receipt;failure=JointAccuracyFailure::none;return true;
}
}
