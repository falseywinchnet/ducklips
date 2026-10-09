#include "ducklips/cavity_network.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <utility>
#include <type_traits>

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
    // Select a breadth-first spanning forest. Zero-coupling links connect its
    // roots only inside the preconditioner, never in the physical graph.
    std::vector<std::vector<std::pair<std::size_t,std::size_t>>> neighbours(nodes);
    for(std::size_t e=0;e<edges.size();++e){
        neighbours[edges[e].first_node].push_back({edges[e].second_node,e});
        neighbours[edges[e].second_node].push_back({edges[e].first_node,e});
    }
    std::vector<bool> visited(nodes,false),selected(edges.size(),false);
    std::vector<std::size_t> queue,roots;queue.reserve(nodes);
    for(std::size_t root=0;root<nodes;++root){
        if(visited[root]){continue;}roots.push_back(root);queue.clear();queue.push_back(root);visited[root]=true;
        for(std::size_t q=0;q<queue.size();++q){
            for(const std::pair<std::size_t,std::size_t>& link:neighbours[queue[q]]){
                if(!visited[link.first]){visited[link.first]=true;selected[link.second]=true;queue.push_back(link.first);}
            }
        }
    }
    direct_tree_=roots.size()==1 && edges.size()+1==nodes;
    std::vector<TreeEdge> forest;
    std::vector<std::size_t> selected_index(edges.size(),cavity_reservoir);
    for(std::size_t e=0;e<edges.size();++e){if(selected[e]){selected_index[e]=forest.size();forest.push_back(edges[e]);}}
    for(std::size_t i=1;i<roots.size();++i){forest.push_back({0,roots[i]});}
    matrix_.prepare(nodes,forest,0);std::size_t edge=0;
    for(std::size_t i=0;i<connections_.size();++i){
        if(connections_[i].to!=cavity_reservoir){
            if(selected_index[edge]!=cavity_reservoir){slots_[i]=matrix_.edge_slot(selected_index[edge]);}++edge;
        }
    }
    component_.resize(nodes);open_degree_.resize(nodes);anchored_.resize(nodes);pinned_.resize(nodes);
    slope_.resize(connections.size());boundary_diagonal_.resize(nodes);
    cg_residual_.resize(nodes);cg_direction_.resize(nodes);cg_image_.resize(nodes);cg_preconditioned_.resize(nodes);
    pressure_.assign(nodes,0);volume_.resize(nodes);compliance_.resize(nodes);flow_.assign(connections.size(),0);
    candidate_.resize(nodes);trial_.resize(nodes);residual_.resize(nodes);gradient_.resize(nodes);delta_.resize(nodes);
    storage_.resize(nodes);rhs_.resize(nodes);roundoff_.resize(nodes);working_flow_.resize(connections.size());linearization_flow_.resize(connections.size());ready_=false;
}
bool CavityNetwork::reconfigure(std::span<const CavityConnection> connections,
                                std::span<const std::size_t> old_node_for_new) {
    if(!ready_ || old_node_for_new.empty()){return false;}
    std::vector<bool> used(pressure_.size(),false);
    std::vector<double> volume(old_node_for_new.size(),0),compliance(old_node_for_new.size(),0),pressure(old_node_for_new.size(),0);
    for(std::size_t next=0;next<old_node_for_new.size();++next){
        const std::size_t old=old_node_for_new[next];if(old==cavity_reservoir){continue;}
        if(old>=pressure_.size() || used[old]){return false;}
        used[old]=true;volume[next]=volume_[old];compliance[next]=compliance_[old];pressure[next]=pressure_[old];
    }
    for(std::size_t old=0;old<pressure_.size();++old){
        if(!used[old] && (volume_[old]!=0 || compliance_[old]!=0)){return false;}
    }
    CavityNetwork rebuilt;
    try{rebuilt.prepare(old_node_for_new.size(),connections);}catch(const std::invalid_argument&){return false;}
    if(!rebuilt.reset(volume,compliance,pressure)){return false;}
    static_assert(std::is_nothrow_move_assignable_v<CavityNetwork>);
    *this=std::move(rebuilt);return true;
}
bool CavityNetwork::reset(std::span<const double> volumes,std::span<const double> compliances,
                          std::span<const double> pressures) noexcept {
    const std::size_t n=pressure_.size();
    if (volumes.size()!=n || compliances.size()!=n || pressures.size()!=n || !n) { return false; }
    for (std::size_t i=0;i<n;++i) {
        if (!std::isfinite(volumes[i]) || volumes[i]<0 || !std::isfinite(compliances[i]) || compliances[i]<0
            || (compliances[i]==0 && volumes[i]!=0)
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
            matrix_.diagonal[i]=pinned_[i]?1:storage_[i];matrix_.parent_coupling[i]=0;boundary_diagonal_[i]=matrix_.diagonal[i];
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
        if(jacobian){linearization_flow_[e]=flow;slope_[e]=0;}
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
            slope_[e]=slope;
            if(c.to==cavity_reservoir){matrix_.diagonal[c.from]+=slope;boundary_diagonal_[c.from]+=slope;}
            if (c.to!=cavity_reservoir) {
                roundoff_[c.to]+=uncertainty;
                if(slots_[e]!=cavity_reservoir){
                    matrix_.diagonal[c.from]+=slope;matrix_.diagonal[c.to]+=slope;matrix_.parent_coupling[slots_[e]]=-slope;
                }else{
                    // A zero-storage cell can lose its selected forest exits
                    // while retaining a cycle edge. Ground only its
                    // preconditioner with that omitted conductance; the actual
                    // Jacobian still contains the full conservative edge.
                    if(storage_[c.from]==0){matrix_.diagonal[c.from]+=slope;}
                    if(storage_[c.to]==0){matrix_.diagonal[c.to]+=slope;}
                }
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
std::size_t CavityNetwork::component_root(std::size_t node) noexcept {
    while(component_[node]!=node){component_[node]=component_[component_[node]];node=component_[node];}
    return node;
}
bool CavityNetwork::constrain_empty_cells(std::span<const CavityPassage> passages) noexcept {
    std::fill(pinned_.begin(),pinned_.end(),false);
    if(std::none_of(storage_.begin(),storage_.end(),[](double c){return c==0;})){return true;}
    for(std::size_t i=0;i<storage_.size();++i){component_[i]=i;open_degree_[i]=0;anchored_[i]=false;}
    for(std::size_t e=0;e<connections_.size();++e){
        if(passages[e].area==0){continue;}const CavityConnection& c=connections_[e];++open_degree_[c.from];
        if(c.to!=cavity_reservoir){++open_degree_[c.to];component_[component_root(c.from)]=component_root(c.to);}
    }
    for(std::size_t i=0;i<storage_.size();++i){if(storage_[i]>0){anchored_[component_root(i)]=true;}}
    for(std::size_t e=0;e<connections_.size();++e){
        if(passages[e].area>0 && connections_[e].to==cavity_reservoir){anchored_[component_root(connections_[e].from)]=true;}
    }
    for(std::size_t i=0;i<storage_.size();++i){
        if(storage_[i]==0 && open_degree_[i]==0){
            // This gauge has neither air nor a passage. No epsilon storage:
            // nonzero inherited mass or imposed inflow must cause failure.
            if(rhs_[i]!=0){return false;}
            pinned_[i]=true;candidate_[i]=0;anchored_[component_root(i)]=true;
        }
    }
    for(std::size_t i=0;i<storage_.size();++i){if(!anchored_[component_root(i)]){return false;}}
    return true;
}
void CavityNetwork::multiply_jacobian(std::span<const double> x,std::span<double> y) const noexcept {
    // Flux differences retain the small storage term without subtracting two
    // large diagonal/off-diagonal sums for nearly uniform pressure.
    for(std::size_t i=0;i<x.size();++i){y[i]=boundary_diagonal_[i]*x[i];}
    for(std::size_t e=0;e<connections_.size();++e){
        const CavityConnection& c=connections_[e];if(c.to==cavity_reservoir){continue;}
        const double flux=slope_[e]*(x[c.from]-x[c.to]);y[c.from]+=flux;y[c.to]-=flux;
    }
}
bool CavityNetwork::solve_jacobian() noexcept {
    if(matrix_.factor(0)!=0){return false;}
    if(direct_tree_){
        for(std::size_t i=0;i<delta_.size();++i){delta_[i]=-residual_[i];}
        matrix_.solve(delta_);return true;
    }
    // Preconditioned conjugate gradients solve the complete SPD Jacobian.
    // The forest is a preconditioner only; omitted cycle edges still contribute
    // to every multiply and the independently recomputed final residual.
    const std::size_t n=delta_.size();double scale=0;
    std::fill(delta_.begin(),delta_.end(),0);
    for(std::size_t i=0;i<n;++i){cg_residual_[i]=-residual_[i];scale=std::max(scale,std::abs(residual_[i]));}
    const double tolerance=1e-17+1e-11*scale;
    double previous=0;bool restart=true;
    for(std::size_t iteration=0;iteration<4*n+64;++iteration){
        std::copy(cg_residual_.begin(),cg_residual_.end(),cg_preconditioned_.begin());matrix_.solve(cg_preconditioned_);
        long double product=0;for(std::size_t i=0;i<n;++i){product+=static_cast<long double>(cg_residual_[i])*cg_preconditioned_[i];}
        const double current=static_cast<double>(product);
        if(!std::isfinite(current) || current<=0){return false;}
        const double beta=restart?0:current/previous;
        for(std::size_t i=0;i<n;++i){cg_direction_[i]=restart?cg_preconditioned_[i]:cg_preconditioned_[i]+beta*cg_direction_[i];}
        multiply_jacobian(cg_direction_,cg_image_);
        product=0;for(std::size_t i=0;i<n;++i){product+=static_cast<long double>(cg_direction_[i])*cg_image_[i];}
        const double denominator=static_cast<double>(product);if(!std::isfinite(denominator) || denominator<=0){return false;}
        const double alpha=current/denominator;double error=0;
        for(std::size_t i=0;i<n;++i){delta_[i]+=alpha*cg_direction_[i];cg_residual_[i]-=alpha*cg_image_[i];error=std::max(error,std::abs(cg_residual_[i]));}
        restart=error<=tolerance || (iteration+1)%64==0;
        if(restart){
            multiply_jacobian(delta_,cg_image_);error=0;
            for(std::size_t i=0;i<n;++i){cg_residual_[i]=-residual_[i]-cg_image_[i];error=std::max(error,std::abs(cg_residual_[i]));}
            if(error<=tolerance){return true;}
        }
        previous=current;
    }
    return false;
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
        if (!std::isfinite(volumes[i]) || volumes[i]<0 || !std::isfinite(compliances[i]) || compliances[i]<0
            || (compliances[i]==0 && volumes[i]!=0)
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
    if(!constrain_empty_cells(passages)){return report;}
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
        if (!std::isfinite(norm) || !solve_jacobian()) { return report; }
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
