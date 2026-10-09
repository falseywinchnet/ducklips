#include "ducklips/moving_graph.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {
bool finite(std::span<const double> values) noexcept {
    return std::all_of(values.begin(),values.end(),[](double value){return std::isfinite(value);});
}
}
void MovingGraphReference::prepare(std::size_t nodes,std::span<const TreeEdge> branches){
    if(nodes==0 || nodes>4096){throw std::invalid_argument("moving graph requires 1..4096 reference nodes");}
    for(const TreeEdge& edge:branches){
        if(edge.first_node>=nodes || (edge.second_node!=ground_node && (edge.second_node>=nodes || edge.second_node==edge.first_node))){
            throw std::invalid_argument("invalid moving graph branch");
        }
    }
    nodes_=nodes;dimension_=2*nodes;branches_.assign(branches.begin(),branches.end());
    matrix_.resize(dimension_*dimension_);original_.resize(matrix_.size());
    right_.resize(dimension_);original_right_.resize(dimension_);stage_.resize(dimension_);
    charge.assign(nodes,0.);momentum.assign(branches.size(),0.);
    average_pressure.assign(nodes,0.);average_flow.assign(branches.size(),0.);
    candidate_charge_.resize(nodes);candidate_momentum_.resize(branches.size());
    candidate_pressure_.resize(nodes);candidate_flow_.resize(branches.size());
}
bool MovingGraphReference::solve() noexcept {
    const std::size_t n=dimension_;
    // Partial pivoting after physical row scaling. No absolute dimensional
    // pivot floor discards a small but valid compliance or passage.
    for(std::size_t column=0;column<n;++column){
        std::size_t pivot=column;
        for(std::size_t row=column+1;row<n;++row){if(std::abs(matrix_[row*n+column])>std::abs(matrix_[pivot*n+column])){pivot=row;}}
        if(!std::isfinite(matrix_[pivot*n+column]) || matrix_[pivot*n+column]==0){return false;}
        if(pivot!=column){for(std::size_t k=column;k<n;++k){std::swap(matrix_[column*n+k],matrix_[pivot*n+k]);}std::swap(right_[column],right_[pivot]);}
        for(std::size_t row=column+1;row<n;++row){
            const double multiplier=matrix_[row*n+column]/matrix_[column*n+column];matrix_[row*n+column]=0;
            for(std::size_t k=column+1;k<n;++k){matrix_[row*n+k]-=multiplier*matrix_[column*n+k];}
            right_[row]-=multiplier*right_[column];
        }
    }
    for(std::size_t row=n;row-->0;){
        double value=right_[row];for(std::size_t k=row+1;k<n;++k){value-=matrix_[row*n+k]*stage_[k];}
        stage_[row]=value/matrix_[row*n+row];if(!std::isfinite(stage_[row])){return false;}
    }
    // Check the original, scaled equations, not the eliminated triangular
    // system. The scale includes cancellation terms, with no pressure floor.
    for(std::size_t row=0;row<n;++row){
        double residual=-original_right_[row],scale=std::abs(original_right_[row]);
        for(std::size_t k=0;k<n;++k){const double value=original_[row*n+k]*stage_[k];residual+=value;scale+=std::abs(value);}
        if(!std::isfinite(residual) || !std::isfinite(scale) || std::abs(residual)>1e-11*scale){return false;}
    }
    return true;
}
bool MovingGraphReference::advance(std::span<const double> c1,std::span<const double> c2,
    std::span<const double> g1,std::span<const double> g2,double dt,
    std::span<const double> nodal_flow,std::span<const double> series_pressure) noexcept {
    const std::size_t n=nodes_,m=branches_.size(),d=dimension_;
    if(!n || c1.size()!=n || c2.size()!=n || g1.size()!=m || g2.size()!=m
        || charge.size()!=n || momentum.size()!=m || average_pressure.size()!=n || average_flow.size()!=m
        || (!nodal_flow.empty() && nodal_flow.size()!=n) || (!series_pressure.empty() && series_pressure.size()!=m)
        || !std::isfinite(dt) || dt<=0 || !finite(c1) || !finite(c2) || !finite(g1) || !finite(g2)
        || !finite(charge) || !finite(momentum) || !finite(nodal_flow) || !finite(series_pressure)){return false;}
    for(std::size_t i=0;i<n;++i){if(c1[i]<=0 || c2[i]<=0){return false;}}
    for(std::size_t i=0;i<m;++i){if(g1[i]<0 || g2[i]<0){return false;}}
    constexpr double a[2][2]={{.25,.25-.28867513459481288225},{.25+.28867513459481288225,.25}};
    std::fill(matrix_.begin(),matrix_.end(),0.);
    for(std::size_t i=0;i<n;++i){
        matrix_[(2*i)*d+2*i]=c1[i];matrix_[(2*i+1)*d+2*i+1]=c2[i];
        const double drive=nodal_flow.empty()?0:nodal_flow[i];
        for(std::size_t r=0;r<2;++r){right_[2*i+r]=charge[i]+dt*(a[r][0]+a[r][1])*drive;}
    }
    for(std::size_t e=0;e<m;++e){
        const TreeEdge edge=branches_[e];const double g[2]={g1[e],g2[e]},source=series_pressure.empty()?0:series_pressure[e];
        for(std::size_t r=0;r<2;++r){
            double load=0;for(std::size_t t=0;t<2;++t){load+=dt*a[r][t]*g[t]*(momentum[e]+dt*(a[t][0]+a[t][1])*source);}
            right_[2*edge.first_node+r]-=load;if(edge.second_node!=ground_node){right_[2*edge.second_node+r]+=load;}
            for(std::size_t s=0;s<2;++s){
                double block=0;for(std::size_t t=0;t<2;++t){block+=dt*dt*a[r][t]*g[t]*a[t][s];}
                matrix_[(2*edge.first_node+r)*d+2*edge.first_node+s]+=block;
                if(edge.second_node!=ground_node){
                    matrix_[(2*edge.second_node+r)*d+2*edge.second_node+s]+=block;
                    matrix_[(2*edge.first_node+r)*d+2*edge.second_node+s]-=block;
                    matrix_[(2*edge.second_node+r)*d+2*edge.first_node+s]-=block;
                }
            }
        }
    }
    for(std::size_t row=0;row<d;++row){
        const double c=row%2?c2[row/2]:c1[row/2];right_[row]/=c;
        for(std::size_t col=0;col<d;++col){matrix_[row*d+col]/=c;}
    }
    if(!finite(matrix_) || !finite(right_)){return false;}
    std::copy(matrix_.begin(),matrix_.end(),original_.begin());std::copy(right_.begin(),right_.end(),original_right_.begin());
    if(!solve()){return false;}
    std::copy(charge.begin(),charge.end(),candidate_charge_.begin());
    for(std::size_t i=0;i<n;++i){candidate_pressure_[i]=.5*(stage_[2*i]+stage_[2*i+1]);}
    for(std::size_t e=0;e<m;++e){
        const TreeEdge edge=branches_[e];double pressure[2];const double source=series_pressure.empty()?0:series_pressure[e];
        for(std::size_t r=0;r<2;++r){pressure[r]=stage_[2*edge.first_node+r]+source;
            if(edge.second_node!=ground_node){pressure[r]-=stage_[2*edge.second_node+r];}}
        const double first=momentum[e]+dt*(a[0][0]*pressure[0]+a[0][1]*pressure[1]);
        const double second=momentum[e]+dt*(a[1][0]*pressure[0]+a[1][1]*pressure[1]);
        candidate_flow_[e]=.5*(g1[e]*first+g2[e]*second);const double transfer=dt*candidate_flow_[e];
        candidate_charge_[edge.first_node]-=transfer;if(edge.second_node!=ground_node){candidate_charge_[edge.second_node]+=transfer;}
        candidate_momentum_[e]=momentum[e]+dt*.5*(pressure[0]+pressure[1]);
    }
    for(std::size_t i=0;i<n;++i){if(!nodal_flow.empty()){candidate_charge_[i]+=dt*nodal_flow[i];}}
    if(!finite(candidate_charge_) || !finite(candidate_momentum_) || !finite(candidate_pressure_) || !finite(candidate_flow_)){return false;}
    std::copy(candidate_charge_.begin(),candidate_charge_.end(),charge.begin());
    std::copy(candidate_momentum_.begin(),candidate_momentum_.end(),momentum.begin());
    std::copy(candidate_pressure_.begin(),candidate_pressure_.end(),average_pressure.begin());
    std::copy(candidate_flow_.begin(),candidate_flow_.end(),average_flow.begin());
    return true;
}
}
