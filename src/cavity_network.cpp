#include "ducklips/cavity_network.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace ducklips {
void CavityNetwork::prepare(std::size_t nodes,std::span<const CavityConnection> connections) {
    ready_=false;
    if (!nodes) { throw std::invalid_argument("cavity network needs at least one node"); }
    std::vector<TreeEdge> edges;
    connections_.assign(connections.begin(),connections.end());slots_.assign(connections.size(),cavity_reservoir);
    for (const CavityConnection& c:connections_) {
        if (c.from>=nodes || (c.to!=cavity_reservoir && (c.to>=nodes || c.to==c.from))) {
            throw std::invalid_argument("invalid cavity connection");
        }
        if (c.to!=cavity_reservoir) { edges.push_back({c.from,c.to}); }
    }
    matrix_.prepare(nodes,edges,0);std::size_t edge=0;
    for (std::size_t i=0;i<connections_.size();++i) {
        if (connections_[i].to!=cavity_reservoir) { slots_[i]=matrix_.edge_slot(edge++); }
    }
    pressure_.assign(nodes,0);volume_.resize(nodes);compliance_.resize(nodes);flow_.assign(connections.size(),0);
    candidate_.resize(nodes);trial_.resize(nodes);residual_.resize(nodes);gradient_.resize(nodes);delta_.resize(nodes);
    storage_.resize(nodes);rhs_.resize(nodes);roundoff_.resize(nodes);working_flow_.resize(connections.size());linearization_flow_.resize(connections.size());ready_=false;
}
bool CavityNetwork::reset(std::span<const double> volumes,std::span<const double> compliances,
                          std::span<const double> pressures) noexcept {
    const std::size_t n=pressure_.size();
    if (volumes.size()!=n || compliances.size()!=n || pressures.size()!=n || !n) { return false; }
    for (std::size_t i=0;i<n;++i) {
        if (!std::isfinite(volumes[i]) || volumes[i]<0 || !std::isfinite(compliances[i]) || compliances[i]<=0
            || !std::isfinite(pressures[i])) { return false; }
    }
    std::copy(volumes.begin(),volumes.end(),volume_.begin());
    std::copy(compliances.begin(),compliances.end(),compliance_.begin());
    std::copy(pressures.begin(),pressures.end(),pressure_.begin());
    std::fill(flow_.begin(),flow_.end(),0);ready_=true;return true;
}
double CavityNetwork::evaluate(std::span<const double> pressure,std::span<const CavityPassage> passages,
                              const AirProperties& air,bool jacobian) noexcept {
    constexpr double eps=std::numeric_limits<double>::epsilon();
    balance_residual_=0;
    if(jacobian){balance_roundoff_=0;}
    for (std::size_t i=0;i<pressure.size();++i) {
        residual_[i]=storage_[i]*pressure[i]-rhs_[i];
        balance_residual_+=residual_[i];
        if (jacobian) {
            matrix_.diagonal[i]=storage_[i];matrix_.parent_coupling[i]=0;
            roundoff_[i]=6*eps*std::abs(storage_[i]*pressure[i])+4*eps*std::abs(rhs_[i]);
            balance_roundoff_+=roundoff_[i];
        }
    }
    for (std::size_t e=0;e<connections_.size();++e) {
        const CavityConnection& c=connections_[e];const CavityPassage& passage=passages[e];
        const double downstream=c.to==cavity_reservoir?passage.reservoir_pressure:pressure[c.to];
        const double drop=pressure[c.from]-downstream;
        const double loss=drop>=0?passage.forward_loss:passage.reverse_loss;
        double flow=0;
        if(passage.area>0){flow=loss>0?orifice_flow(passage.area/std::sqrt(loss),drop,passage.resistance,air):drop/passage.resistance;}
        working_flow_[e]=flow;residual_[c.from]+=flow;
        if(jacobian){linearization_flow_[e]=flow;}
        if (c.to!=cavity_reservoir) { residual_[c.to]-=flow; }
        else {balance_residual_+=flow;}
        if (jacobian && passage.area>0) {
            const double slope=1/(passage.resistance+loss*air.density*std::abs(flow)/(passage.area*passage.area));
            // Nearly equal node pressures cannot encode arbitrarily small
            // drops. Bound the induced flow error by local sensitivity times
            // pressure representation error, in addition to arithmetic error.
            const double uncertainty=2*eps*slope*(std::abs(pressure[c.from])+std::abs(downstream))+4*eps*std::abs(flow);
            roundoff_[c.from]+=uncertainty;
            if(c.to==cavity_reservoir){balance_roundoff_+=uncertainty;}
            matrix_.diagonal[c.from]+=slope;
            if (c.to!=cavity_reservoir) {
                roundoff_[c.to]+=uncertainty;
                matrix_.diagonal[c.to]+=slope;matrix_.parent_coupling[slots_[e]]=-slope;
            }
        }
    }
    double norm=0;
    for (double value:residual_) {
        if (!std::isfinite(value)) { return INFINITY; }
        norm=std::max(norm,std::abs(value));
    }
    if(jacobian){std::copy(residual_.begin(),residual_.end(),gradient_.begin());}
    return norm;
}
bool CavityNetwork::balanced(double tolerance) const noexcept {
    for(std::size_t i=0;i<residual_.size();++i){
        if(!std::isfinite(residual_[i]) || !std::isfinite(roundoff_[i]) || std::abs(residual_[i])>tolerance+roundoff_[i]){return false;}
    }
    // Internal port errors cancel in total mass. Do not let their large local
    // pressure-quantization allowance excuse an error at external boundaries.
    return std::isfinite(balance_residual_) && std::isfinite(balance_roundoff_)
        && std::abs(balance_residual_)<=double(residual_.size())*tolerance+balance_roundoff_;
}
double CavityNetwork::potential_change(std::span<const CavityPassage> passages,const AirProperties& air) const noexcept {
    // The residual is the gradient of a strictly convex storage + passage
    // potential. Use that potential to globalize Newton, not the largest
    // individual mass residual, which can temporarily rise on a good step.
    // Use gradient dot displacement plus a nonnegative Bregman remainder.
    // This avoids cancelling large storage and passage work near convergence.
    double change=0;
    for(std::size_t i=0;i<candidate_.size();++i){
        const double d=trial_[i]-candidate_[i];
        change+=d*gradient_[i]+.5*storage_[i]*d*d;
    }
    for(std::size_t e=0;e<passages.size();++e){
        const CavityPassage& p=passages[e];if(p.area==0){continue;}
        const double q0=linearization_flow_[e],q1=working_flow_[e],a0=std::abs(q0),a1=std::abs(q1);
        const double loss0=q0>=0?p.forward_loss:p.reverse_loss,loss1=q1>=0?p.forward_loss:p.reverse_loss;
        const double dq=q1-q0,scale=air.density/(2*p.area*p.area);
        change+=.5*p.resistance*dq*dq;
        if(q0*q1>=0 && loss0==loss1){change+=scale*loss0*dq*dq*(a0+2*a1)/3;}
        else{change+=scale*(loss0*a0*a0*a0/3+loss1*a1*a1*(a0+2*a1/3));}
    }
    return change;
}
CavityStepReport CavityNetwork::advance(std::span<const double> volumes,std::span<const double> compliances,
    std::span<const CavityPassage> passages,std::span<const double> inflow,double dt,const AirProperties& air) noexcept {
    CavityStepReport report{};const std::size_t n=pressure_.size();
    if (!ready_ || volumes.size()!=n || compliances.size()!=n || inflow.size()!=n || passages.size()!=connections_.size()
        || !std::isfinite(dt) || dt<=0 || !std::isfinite(air.density) || air.density<=0) { return report; }
    double scale=0;
    for (std::size_t i=0;i<n;++i) {
        if (!std::isfinite(volumes[i]) || volumes[i]<0 || !std::isfinite(compliances[i]) || compliances[i]<=0
            || !std::isfinite(inflow[i])) { return report; }
        storage_[i]=compliances[i]/dt;
        rhs_[i]=(compliance_[i]*pressure_[i]-(volumes[i]-volume_[i]))/dt+inflow[i];
        if (!std::isfinite(storage_[i]) || !std::isfinite(rhs_[i])) { return report; }
        scale=std::max(scale,std::abs(rhs_[i]));candidate_[i]=pressure_[i];
    }
    for (const CavityPassage& p:passages) {
        if (!std::isfinite(p.area) || p.area<0 || !std::isfinite(p.resistance) || p.resistance<0
            || (p.area>0 && p.resistance==0) || !std::isfinite(p.reservoir_pressure)) { return report; }
        if(!std::isfinite(p.forward_loss) || p.forward_loss<0 || !std::isfinite(p.reverse_loss) || p.reverse_loss<0){return report;}
    }
    const double tolerance=1e-14+1e-11*scale;
    for (std::size_t iteration=0;iteration<64;++iteration) {
        const double norm=evaluate(candidate_,passages,air,true);
        report.maximum_mass_residual=norm;report.iterations=iteration;
        if (balanced(tolerance)) {
            std::copy(candidate_.begin(),candidate_.end(),pressure_.begin());
            std::copy(volumes.begin(),volumes.end(),volume_.begin());
            std::copy(compliances.begin(),compliances.end(),compliance_.begin());
            std::copy(working_flow_.begin(),working_flow_.end(),flow_.begin());
            report.converged=true;return report;
        }
        if (!std::isfinite(norm) || matrix_.factor(0)!=0) { return report; }
        for (std::size_t i=0;i<n;++i) { delta_[i]=-residual_[i]; }
        matrix_.solve(delta_);
        double descent=0;for(std::size_t i=0;i<n;++i){descent+=residual_[i]*delta_[i];}
        bool accepted=false;double fraction=1;
        for (std::size_t search=0;search<24;++search) {
            for (std::size_t i=0;i<n;++i) { trial_[i]=candidate_[i]+fraction*delta_[i]; }
            const double next=evaluate(trial_,passages,air,false);
            const double change=potential_change(passages,air);
            if (std::isfinite(next) && (balanced(tolerance) || (std::isfinite(change) && change<=1e-4*fraction*descent))) {
                std::copy(trial_.begin(),trial_.end(),candidate_.begin());accepted=true;break;
            }
            fraction*=.5;
        }
        if (!accepted) { return report; }
    }
    return report;
}
}
