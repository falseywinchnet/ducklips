#include "ducklips/graph_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {
bool solve(std::array<double,9> a,std::array<double,3> b,std::array<double,3>& x) noexcept {
    for(std::size_t j=0;j<3;++j){
        std::size_t pivot=j;for(std::size_t i=j+1;i<3;++i){if(std::abs(a[3*i+j])>std::abs(a[3*pivot+j])){pivot=i;}}
        if(!std::isfinite(a[3*pivot+j]) || a[3*pivot+j]==0){return false;}
        for(std::size_t k=j;k<3;++k){std::swap(a[3*j+k],a[3*pivot+k]);}std::swap(b[j],b[pivot]);
        for(std::size_t i=j+1;i<3;++i){const double f=a[3*i+j]/a[3*j+j];
            for(std::size_t k=j+1;k<3;++k){a[3*i+k]-=f*a[3*j+k];}b[i]-=f*b[j];}
    }
    for(std::size_t j=3;j-->0;){double sum=b[j];for(std::size_t k=j+1;k<3;++k){sum-=a[3*j+k]*x[k];}x[j]=sum/a[3*j+j];if(!std::isfinite(x[j])){return false;}}
    return true;
}
}
void GraphGlottalStorage::prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
    const FoldParameters& tissue,const AirProperties& air,double tension,
    const FoldState& initial,const std::array<double,2>& pressure,const GlottalPassages& passages){
    if(source_node>=nodes){throw std::invalid_argument("glottal source outside graph");}
    for(std::size_t j=0;j<3;++j){
        if(!std::isfinite(passages.inertance[j]) || passages.inertance[j]<=0
            || !std::isfinite(passages.resistance[j]) || passages.resistance[j]<0
            || !std::isfinite(passages.quadratic_loss[j]) || passages.quadratic_loss[j]<0){throw std::invalid_argument("invalid fixed glottal passage");}
    }
    GraphGlottalStorage candidate;
    if(!candidate.gas_.configure(air) || !candidate.chambers.prepare(tissue,air,tension,initial,pressure)){throw std::invalid_argument("invalid glottal chambers");}
    candidate.network.prepare(nodes,branches);candidate.free_.prepare(nodes,branches);candidate.response_.prepare(nodes,branches);
    candidate.unit_flow_.assign(nodes,0.);candidate.unit_flow_[source_node]=1.;candidate.source_node_=source_node;candidate.passages_=passages;
    *this=std::move(candidate);
}
void GraphGlottalStorage::prepare_articulated(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
    const FoldParameters& tissue,const AirProperties& air,double tension,
    const FoldState& initial,const std::array<double,2>& pressure,double outlet_loss){
    if(source_node>=nodes){throw std::invalid_argument("glottal source outside graph");}
    GraphGlottalStorage candidate;
    if(!candidate.gas_.configure(air) || !candidate.chambers.prepare(tissue,air,tension,initial,pressure)
        || !candidate.kinetics_.configure(tissue,air,outlet_loss)){throw std::invalid_argument("invalid articulated glottal passages");}
    candidate.network.prepare(nodes,branches);candidate.free_.prepare(nodes,branches);candidate.response_.prepare(nodes,branches);
    candidate.unit_flow_.assign(nodes,0.);candidate.unit_flow_[source_node]=1.;candidate.source_node_=source_node;candidate.articulated_=true;
    *this=std::move(candidate);
}
bool GraphGlottalStorage::source_energy(double& result) const noexcept {
    double e=0;if(!chambers.energy(e)){return false;}
    if(articulated_){GlottalKineticState kinetic;if(!kinetics_.evaluate(chambers.state.mass,momentum,kinetic)){return false;}e+=kinetic.energy;}
    else{for(std::size_t j=0;j<3;++j){e+=.5*momentum[j]*momentum[j]/passages_.inertance[j];}}
    if(!std::isfinite(e)){return false;}result=e;return true;
}
bool GraphGlottalStorage::advance(std::span<const double> c1,std::span<const double> c2,
    std::span<const double> g1,std::span<const double> g2,double dt,double reservoir_pressure,
    double rest_area,const std::array<double,4>& external_force,
    std::span<const double> nodal_flow,std::span<const double> series_pressure) noexcept {
    if(unit_flow_.empty() || !std::isfinite(dt) || dt<=0
        || network.charge.size()!=free_.charge.size() || network.momentum.size()!=free_.momentum.size()
        || network.average_pressure.size()!=free_.average_pressure.size() || network.average_flow.size()!=free_.average_flow.size()){return false;}
    double mass=0,energy0=0;GasEnergy reservoir;
    if(!gas_.mass_at_pressure(1.,reservoir_pressure,mass) || !gas_.evaluate(mass,1.,reservoir) || !source_energy(energy0)){return false;}
    std::copy(network.charge.begin(),network.charge.end(),free_.charge.begin());
    std::copy(network.momentum.begin(),network.momentum.end(),free_.momentum.begin());
    std::fill(response_.charge.begin(),response_.charge.end(),0.);std::fill(response_.momentum.begin(),response_.momentum.end(),0.);
    if(!free_.advance(c1,c2,g1,g2,dt,nodal_flow,series_pressure) || !response_.advance(c1,c2,g1,g2,dt,unit_flow_)){return false;}
    const double pfree=free_.average_pressure[source_node_],z=response_.average_pressure[source_node_];
    if(!std::isfinite(z) || z<0){return false;}
    GlottalStorage trial;
    GlottalKineticState kinetic0;
    std::array<double,2> volume0{};
    if(articulated_ && (!kinetics_.evaluate(chambers.state.mass,momentum,kinetic0) || !chambers.volumes(volume0))){return false;}
    std::array<double,3> trial_momentum{},resistance{},quadratic{},kinetic_port{};
    GlottalKineticTransfer kinetic_transfer;
    std::array<double,3> q{},residual{};double norm=0,scale=0;
    // Zero flow is mass-feasible even if the old momenta point towards a
    // nearly empty cell. It is only an implicit-solve starting guess.
    const auto evaluate=[&](const std::array<double,3>& flow){
        trial=chambers;if(!trial.advance(flow,rest_area,external_force,dt)){return false;}
        resistance=passages_.resistance;quadratic=passages_.quadratic_loss;kinetic_port={};
        if(articulated_){
            GlottalKineticState kinetic1;std::array<double,2> volume1{},mid_mass{},mid_volume{};
            if(!kinetics_.evaluate(trial.state.mass,{},kinetic1) || !trial.volumes(volume1)){return false;}
            for(std::size_t j=0;j<3;++j){trial_momentum[j]=4*flow[j]/(kinetic0.inverse_inertance[j]+kinetic1.inverse_inertance[j])-momentum[j];}
            if(!kinetics_.transfer(chambers.state.mass,trial.state.mass,momentum,trial_momentum,kinetic_transfer)){return false;}
            for(std::size_t i=0;i<2;++i){mid_mass[i]=.5*(chambers.state.mass[i]+trial.state.mass[i]);mid_volume[i]=.5*(volume0[i]+volume1[i]);}
            if(!kinetics_.losses(mid_mass,mid_volume,resistance,quadratic)){return false;}
            const auto k=kinetic_transfer.mass_effort;kinetic_port={k[0],k[1]-k[0],-k[1]};
        }else{for(std::size_t j=0;j<3;++j){trial_momentum[j]=2*passages_.inertance[j]*flow[j]-momentum[j];}}
        norm=0;scale=1.;
        for(std::size_t j=0;j<3;++j){
            const double inertia=articulated_?(trial_momentum[j]-momentum[j])/dt:2*(passages_.inertance[j]*flow[j]-momentum[j])/dt;
            const double loss=resistance[j]*flow[j]+quadratic[j]*flow[j]*std::abs(flow[j]);
            const double drive=j==0?reservoir.mass_effort:j==2?-(pfree+z*flow[2]):0.;
            residual[j]=inertia+loss+trial.last.flow_effort[j]+kinetic_port[j]-drive;
            if(!std::isfinite(residual[j])){return false;}
            norm=std::max(norm,std::abs(residual[j]));scale=std::max(scale,std::abs(inertia)+std::abs(loss)+std::abs(trial.last.flow_effort[j])+std::abs(kinetic_port[j])+std::abs(drive));
        }
        return std::isfinite(scale);
    };
    std::size_t iteration=0;
    for(;iteration<48;++iteration){
        if(!evaluate(q)){return false;}if(norm<=1e-8+1e-12*scale){break;}
        const auto base=residual;const double old_norm=norm;std::array<double,9> jacobian{};
        for(std::size_t j=0;j<3;++j){
            double h=1e-5*std::max({std::abs(q[j]),.01*std::min(chambers.state.mass[0],chambers.state.mass[1])/dt,1e-8});
            bool differentiated=false;
            for(int attempt=0;attempt<12;++attempt){
                auto plus=q,minus=q;plus[j]+=h;minus[j]-=h;
                const bool plus_ok=evaluate(plus);const auto rp=residual;
                const bool minus_ok=evaluate(minus);const auto rm=residual;
                if(plus_ok && minus_ok){for(std::size_t i=0;i<3;++i){jacobian[3*i+j]=(rp[i]-rm[i])/(2*h);}differentiated=true;break;}
                h*=.5;
            }
            if(!differentiated){return false;}
        }
        auto rhs=base;for(double& v:rhs){v=-v;}std::array<double,3> delta{};if(!solve(jacobian,rhs,delta)){return false;}
        double fraction=1.;bool accepted=false;
        for(int line=0;line<48;++line){auto next=q;for(std::size_t j=0;j<3;++j){next[j]+=fraction*delta[j];}
            if(evaluate(next) && (norm<old_norm || norm<=1e-8+1e-12*scale)){q=next;accepted=true;break;}fraction*=.5;
        }
        if(!accepted){return false;}
    }
    if(iteration==48){return false;}
    GraphGlottalStep receipt;receipt.flow=q;receipt.law_residual=residual;receipt.iterations=iteration+1;
    receipt.load_effort=pfree+z*q[2];receipt.reservoir_effort=reservoir.mass_effort;
    receipt.reservoir_work=dt*reservoir.mass_effort*q[0];receipt.network_work=dt*receipt.load_effort*q[2];
    const auto next_momentum=trial_momentum;double energy1=0;if(!trial.energy(energy1)){return false;}
    if(articulated_){GlottalKineticState kinetic1;if(!kinetics_.evaluate(trial.state.mass,next_momentum,kinetic1)){return false;}
        energy1+=kinetic1.energy;receipt.kinetic_mass_effort=kinetic_transfer.mass_effort;
    }else{for(std::size_t j=0;j<3;++j){energy1+=.5*next_momentum[j]*next_momentum[j]/passages_.inertance[j];}}
    for(std::size_t j=0;j<3;++j){receipt.passage_loss+=dt*(resistance[j]+quadratic[j]*std::abs(q[j]))*q[j]*q[j];}
    receipt.energy_change=energy1-energy0;
    receipt.energy_residual=receipt.energy_change+trial.last.damping_loss+receipt.passage_loss+receipt.network_work
        -receipt.reservoir_work-trial.last.actuator_work-trial.last.external_work;
    if(!std::isfinite(receipt.energy_residual)){return false;}
    const auto combine=[&](std::vector<double>& a,const std::vector<double>& b){for(std::size_t i=0;i<a.size();++i){a[i]+=q[2]*b[i];if(!std::isfinite(a[i])){return false;}}return true;};
    if(!combine(free_.charge,response_.charge) || !combine(free_.momentum,response_.momentum)
        || !combine(free_.average_pressure,response_.average_pressure) || !combine(free_.average_flow,response_.average_flow)){return false;}
    std::copy(free_.charge.begin(),free_.charge.end(),network.charge.begin());std::copy(free_.momentum.begin(),free_.momentum.end(),network.momentum.begin());
    std::copy(free_.average_pressure.begin(),free_.average_pressure.end(),network.average_pressure.begin());std::copy(free_.average_flow.begin(),free_.average_flow.end(),network.average_flow.begin());
    chambers=trial;momentum=next_momentum;last=receipt;return true;
}
}
