#include "ducklips/relaxation_graph.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace ducklips {
namespace {
using Pair=std::array<double,2>;using Block=std::array<double,4>;
constexpr Block inverse_a{1.5,.5,-4.5,2.5};
Pair multiply(const Block& a,const Pair& x){return {a[0]*x[0]+a[1]*x[1],a[2]*x[0]+a[3]*x[1]};}
bool finite(std::span<const double> a){for(double v:a){if(!std::isfinite(v)){return false;}}return true;}
bool inverse(const Block& a,Block& result){
    const double scale=std::max({std::abs(a[0]),std::abs(a[1]),std::abs(a[2]),std::abs(a[3])});
    if(!std::isfinite(scale) || scale==0){return false;}
    Block b;for(std::size_t i=0;i<4;++i){b[i]=a[i]/scale;}
    const double det=b[0]*b[3]-b[1]*b[2];if(!std::isfinite(det) || det==0){return false;}
    result={b[3]/det/scale,-b[1]/det/scale,-b[2]/det/scale,b[0]/det/scale};return finite(result);
}
bool kernel(RelaxationMetric metric,double initial,double dt,Block& decay,Block& rate){
    const double k=metric.storage,lambda=metric.rate;
    if(!std::isfinite(k) || !std::isfinite(lambda) || !std::isfinite(initial) || k<0 || lambda<0){return false;}
    if(k==0 || lambda==0){if(k!=0 || lambda!=0 || initial!=0){return false;}decay={};rate={};return true;}
    const double x=dt*lambda;
    if(x<=1){
        const double denominator=1+(2./3)*x+x*x/6;
        decay={(1+x/4)/denominator,x/12/denominator,-.75*x/denominator,(1+5*x/12)/denominator};
        for(std::size_t i=0;i<4;++i){rate[i]=lambda*decay[i];}
    }else{
        // lambda*(I+dt*lambda*A)^-1, avoiding overflow and subtraction of
        // nearly equal equilibrated terms when the relaxation is stiff.
        const double v=1/x,denominator=v*v+(2./3)*v+1./6;
        const Block scaled{v+.25,1./12,-.75,v+5./12};
        for(std::size_t i=0;i<4;++i){rate[i]=scaled[i]/denominator/dt;decay[i]=rate[i]/lambda;}
    }
    return finite(decay) && finite(rate);
}
double mean(const Pair& a){return .75*a[0]+.25*a[1];}
double square(double x){return x*x;}
}
void RelaxationGraphReference::prepare(std::size_t nodes,std::span<const TreeEdge> edges,
    std::span<const std::size_t> node_owners,std::span<const std::size_t> branch_owners){
    if(!nodes || nodes>4096){throw std::invalid_argument("relaxation graph node count");}
    for(const auto& edge:edges){if(edge.first_node>=nodes || (edge.second_node!=ground_node && (edge.second_node>=nodes || edge.second_node==edge.first_node))){throw std::invalid_argument("relaxation graph edge");}}
    for(auto i:node_owners){if(i>=nodes){throw std::invalid_argument("relaxation node owner");}}
    for(auto i:branch_owners){if(i>=edges.size()){throw std::invalid_argument("relaxation branch owner");}}
    nodes_=nodes;edges_.assign(edges.begin(),edges.end());node_owners_.assign(node_owners.begin(),node_owners.end());branch_owners_.assign(branch_owners.begin(),branch_owners.end());
    const auto n=nodes,m=edges.size(),a=node_owners.size(),b=branch_owners.size();
    charge.assign(n,0);momentum.assign(m,0);node_memory.assign(a,0);branch_memory.assign(b,0);average_pressure.assign(n,0);average_flow.assign(m,0);
    node_dissipation.assign(a,0);branch_dissipation.assign(b,0);resistance_dissipation.assign(m,0);
    matrix_.resize(4*n*n);original_.resize(matrix_.size());right_.resize(2*n);original_right_.resize(2*n);pressure_.resize(2*n);
    next_q_.resize(n);next_p_.resize(m);next_node_.resize(a);next_branch_.resize(b);
    trial_node_heat_.resize(a);trial_branch_heat_.resize(b);trial_resistance_heat_.resize(m);
    node_kernel_.resize(a);node_decay_.resize(a);branch_kernel_.resize(b);branch_decay_.resize(b);admittance_.resize(m);
    free_flow_.resize(m);flow_.resize(m);dq_.resize(n);dp_.resize(m);node_drive_.resize(n);branch_drive_.resize(m);
}
bool RelaxationGraphReference::solve() noexcept {
    const auto n=right_.size();
    for(std::size_t i=0;i<n;++i){double scale=0;for(std::size_t j=0;j<n;++j){scale=std::max(scale,std::abs(matrix_[i*n+j]));}
        if(!std::isfinite(scale) || scale==0){return false;}right_[i]/=scale;for(std::size_t j=0;j<n;++j){matrix_[i*n+j]/=scale;}}
    original_=matrix_;original_right_=right_;
    for(std::size_t col=0;col<n;++col){std::size_t pivot=col;
        for(std::size_t row=col+1;row<n;++row){if(std::abs(matrix_[row*n+col])>std::abs(matrix_[pivot*n+col])){pivot=row;}}
        if(!std::isfinite(matrix_[pivot*n+col]) || matrix_[pivot*n+col]==0){return false;}
        if(pivot!=col){for(std::size_t j=col;j<n;++j){std::swap(matrix_[col*n+j],matrix_[pivot*n+j]);}std::swap(right_[col],right_[pivot]);}
        for(std::size_t row=col+1;row<n;++row){const double factor=matrix_[row*n+col]/matrix_[col*n+col];matrix_[row*n+col]=0;
            for(std::size_t j=col+1;j<n;++j){matrix_[row*n+j]-=factor*matrix_[col*n+j];}right_[row]-=factor*right_[col];}
    }
    for(std::size_t i=n;i-->0;){double value=right_[i];for(std::size_t j=i+1;j<n;++j){value-=matrix_[i*n+j]*pressure_[j];}pressure_[i]=value/matrix_[i*n+i];}
    if(!finite(pressure_)){return false;}
    for(std::size_t i=0;i<n;++i){double residual=-original_right_[i],scale=std::abs(original_right_[i]);
        for(std::size_t j=0;j<n;++j){const double x=original_[i*n+j]*pressure_[j];residual+=x;scale+=std::abs(x);}
        if(!std::isfinite(residual) || std::abs(residual)>1e-11*scale){return false;}}
    return true;
}
bool RelaxationGraphReference::advance(std::span<const double> c,std::span<const double> g,std::span<const double> resistance,
    std::span<const RelaxationMetric> node_metrics,std::span<const RelaxationMetric> branch_metrics,double dt,
    std::span<const double> nodal,std::span<const double> series) noexcept {
    const auto n=nodes_,m=edges_.size(),d=2*n;
    if(!n || c.size()!=n || g.size()!=m || resistance.size()!=m || node_metrics.size()!=node_owners_.size() || branch_metrics.size()!=branch_owners_.size()
        || charge.size()!=n || momentum.size()!=m || node_memory.size()!=node_metrics.size() || branch_memory.size()!=branch_metrics.size()
        || average_pressure.size()!=n || average_flow.size()!=m || node_dissipation.size()!=node_metrics.size() || branch_dissipation.size()!=branch_metrics.size() || resistance_dissipation.size()!=m
        || (!nodal.empty() && nodal.size()!=n) || (!series.empty() && series.size()!=m) || !std::isfinite(dt) || dt<=0
        || !finite(c) || !finite(g) || !finite(resistance) || !finite(charge) || !finite(momentum) || !finite(nodal) || !finite(series)){return false;}
    for(double value:c){if(value<=0){return false;}}for(double value:g){if(value<0){return false;}}for(double value:resistance){if(value<0){return false;}}
    std::fill(matrix_.begin(),matrix_.end(),0.);std::fill(node_drive_.begin(),node_drive_.end(),Pair{});std::fill(branch_drive_.begin(),branch_drive_.end(),Pair{});
    for(std::size_t i=0;i<n;++i){const double source=nodal.empty()?0:nodal[i];
        right_[2*i]=2*charge[i]/dt+source;right_[2*i+1]=-2*charge[i]/dt+source;
        for(std::size_t r=0;r<2;++r){for(std::size_t s=0;s<2;++s){matrix_[(2*i+r)*d+2*i+s]=c[i]/dt*inverse_a[2*r+s];}}}
    for(std::size_t i=0;i<node_metrics.size();++i){
        if(!kernel(node_metrics[i],node_memory[i],dt,node_decay_[i],node_kernel_[i])){return false;}const auto owner=node_owners_[i];
        const Pair drive=multiply(node_kernel_[i],{node_memory[i],node_memory[i]});
        for(std::size_t r=0;r<2;++r){right_[2*owner+r]+=drive[r];node_drive_[owner][r]+=drive[r];
            for(std::size_t s=0;s<2;++s){matrix_[(2*owner+r)*d+2*owner+s]+=node_metrics[i].storage*node_kernel_[i][2*r+s];}}
    }
    // Build each branch's two-stage impedance, accumulating all shared modes.
    for(std::size_t e=0;e<m;++e){admittance_[e]={resistance[e],0,0,resistance[e]};}
    for(std::size_t i=0;i<branch_metrics.size();++i){
        if(!kernel(branch_metrics[i],branch_memory[i],dt,branch_decay_[i],branch_kernel_[i])){return false;}const auto e=branch_owners_[i];
        const Pair drive=multiply(branch_kernel_[i],{branch_memory[i],branch_memory[i]});
        for(std::size_t r=0;r<2;++r){branch_drive_[e][r]+=drive[r];}
        for(std::size_t j=0;j<4;++j){admittance_[e][j]+=branch_metrics[i].storage*branch_kernel_[i][j];}
    }
    for(std::size_t e=0;e<m;++e){const auto edge=edges_[e];
        if(g[e]==0){admittance_[e]={};free_flow_[e]={};continue;}
        for(std::size_t j=0;j<4;++j){admittance_[e][j]+=inverse_a[j]/dt/g[e];}
        Block inverse_block;if(!inverse(admittance_[e],inverse_block)){return false;}admittance_[e]=inverse_block;
        const double source=series.empty()?0:series[e];
        free_flow_[e]=multiply(inverse_block,{branch_drive_[e][0]+source+2*momentum[e]/dt,branch_drive_[e][1]+source-2*momentum[e]/dt});
        for(std::size_t r=0;r<2;++r){right_[2*edge.first_node+r]-=free_flow_[e][r];if(edge.second_node!=ground_node){right_[2*edge.second_node+r]+=free_flow_[e][r];}
            for(std::size_t s=0;s<2;++s){const double value=inverse_block[2*r+s];matrix_[(2*edge.first_node+r)*d+2*edge.first_node+s]+=value;
                if(edge.second_node!=ground_node){matrix_[(2*edge.second_node+r)*d+2*edge.second_node+s]+=value;matrix_[(2*edge.first_node+r)*d+2*edge.second_node+s]-=value;matrix_[(2*edge.second_node+r)*d+2*edge.first_node+s]-=value;}}}
    }
    if(!finite(matrix_) || !finite(right_) || !solve()){return false;}
    RelaxationGraphStep receipt;double numeric=0;
    for(std::size_t i=0;i<n;++i){const double source=nodal.empty()?0:nodal[i];next_q_[i]=charge[i]+dt*source;dq_[i]={source,source};receipt.external_work+=dt*source*(.75*pressure_[2*i]+.25*pressure_[2*i+1]);}
    for(std::size_t e=0;e<m;++e){const auto edge=edges_[e];Pair drop{pressure_[2*edge.first_node],pressure_[2*edge.first_node+1]};
        if(edge.second_node!=ground_node){drop[0]-=pressure_[2*edge.second_node];drop[1]-=pressure_[2*edge.second_node+1];}
        const Pair response=multiply(admittance_[e],drop);flow_[e]={free_flow_[e][0]+response[0],free_flow_[e][1]+response[1]};
        const double transfer=dt*mean(flow_[e]),source=series.empty()?0:series[e];
        next_q_[edge.first_node]-=transfer;if(edge.second_node!=ground_node){next_q_[edge.second_node]+=transfer;}
        for(std::size_t r=0;r<2;++r){dq_[edge.first_node][r]-=flow_[e][r];if(edge.second_node!=ground_node){dq_[edge.second_node][r]+=flow_[e][r];}
            dp_[e][r]=g[e]>0?drop[r]+source-resistance[e]*flow_[e][r]:0;}
        next_p_[e]=g[e]>0?flow_[e][1]/g[e]:0;
        trial_resistance_heat_[e]=dt*resistance[e]*(.75*square(flow_[e][0])+.25*square(flow_[e][1]));receipt.physical_dissipation+=trial_resistance_heat_[e];receipt.external_work+=source*transfer;
    }
    for(unsigned kind=0;kind<2;++kind){const bool node=kind==0;
        const auto metrics=node?node_metrics:branch_metrics;const auto& owners=node?node_owners_:branch_owners_;
        const auto& memory=node?node_memory:branch_memory;auto& next=node?next_node_:next_branch_;auto& heat=node?trial_node_heat_:trial_branch_heat_;
        const auto& kernels=node?node_kernel_:branch_kernel_;const auto& decays=node?node_decay_:branch_decay_;
        for(std::size_t i=0;i<metrics.size();++i){const auto owner=owners[i];const double k=metrics[i].storage,lambda=metrics[i].rate;
            if(k==0){next[i]=0;heat[i]=0;continue;}
            const Pair drive=node?Pair{pressure_[2*owner],pressure_[2*owner+1]}:flow_[owner];
            const Pair initial=multiply(decays[i],{memory[i],memory[i]});
            // I-T = dt*A*(lambda*T), evaluated without subtracting two ones.
            const Pair response=multiply(kernels[i],drive);
            const Pair stage{initial[0]+dt*k*(5*response[0]-response[1])/12,initial[1]+dt*k*(3*response[0]+response[1])/4};
            const Pair derivative=multiply(kernels[i],{k*drive[0]-memory[i],k*drive[1]-memory[i]});
            next[i]=stage[1];heat[i]=dt*(.75*square(derivative[0])+.25*square(derivative[1]))/k/lambda;
            receipt.physical_dissipation+=heat[i];numeric+=square(dt*(derivative[0]-derivative[1]))/32/k;
            receipt.energy_change+=(next[i]-memory[i])*(next[i]+memory[i])/(2*k);
            if(node){next_q_[owner]-=next[i]-memory[i];for(std::size_t r=0;r<2;++r){dq_[owner][r]-=derivative[r];}}
            else if(g[owner]>0){for(std::size_t r=0;r<2;++r){dp_[owner][r]-=derivative[r];}}
        }
    }
    double balance=0;
    for(std::size_t i=0;i<n;++i){receipt.energy_change+=(next_q_[i]-charge[i])*(next_q_[i]+charge[i])/(2*c[i]);numeric+=square(dt*(dq_[i][0]-dq_[i][1]))/32/c[i];balance+=next_q_[i]-charge[i]-dt*(nodal.empty()?0:nodal[i]);}
    for(std::size_t i=0;i<node_memory.size();++i){balance+=next_node_[i]-node_memory[i];}
    for(std::size_t e=0;e<m;++e){receipt.energy_change+=.5*g[e]*(next_p_[e]-momentum[e])*(next_p_[e]+momentum[e]);numeric+=g[e]*square(dt*(dp_[e][0]-dp_[e][1]))/32;if(edges_[e].second_node==ground_node){balance+=dt*mean(flow_[e]);}}
    receipt.numerical_dissipation=numeric;receipt.air_residual=balance;receipt.energy_residual=receipt.energy_change+receipt.physical_dissipation+numeric-receipt.external_work;
    if(!finite(next_q_) || !finite(next_p_) || !finite(next_node_) || !finite(next_branch_) || !std::isfinite(receipt.energy_residual) || !std::isfinite(balance) || !std::isfinite(numeric) || !std::isfinite(receipt.physical_dissipation)){return false;}
    charge=next_q_;momentum=next_p_;node_memory=next_node_;branch_memory=next_branch_;
    node_dissipation=trial_node_heat_;branch_dissipation=trial_branch_heat_;resistance_dissipation=trial_resistance_heat_;
    for(std::size_t i=0;i<n;++i){average_pressure[i]=.75*pressure_[2*i]+.25*pressure_[2*i+1];}
    for(std::size_t e=0;e<m;++e){average_flow[e]=mean(flow_[e]);}last=receipt;return true;
}
}
