#include "ducklips/moving_tree.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {
using Matrix = std::array<double,4>;
using Pair = std::array<double,2>;
Matrix product(const Matrix& a,const Matrix& b) noexcept {
    return {a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3]};
}
Pair product(const Matrix& a,const Pair& b) noexcept {
    return {a[0]*b[0]+a[1]*b[1],a[2]*b[0]+a[3]*b[1]};
}
Matrix inverse(const Matrix& a) noexcept {
    // Scale first: physical acoustic compliance is commonly near 1e-11.
    const double scale=std::max({std::abs(a[0]),std::abs(a[1]),std::abs(a[2]),std::abs(a[3])});
    const Matrix b{a[0]/scale,a[1]/scale,a[2]/scale,a[3]/scale};
    const double reciprocal=1/(scale*(b[0]*b[3]-b[1]*b[2]));
    return {b[3]*reciprocal,-b[1]*reciprocal,-b[2]*reciprocal,b[0]*reciprocal};
}
}
void MovingTreeReference::prepare(std::size_t nodes,std::span<const TreeEdge> branches) {
    if (nodes==0) { throw std::invalid_argument("moving tree requires nodes"); }
    std::vector<std::vector<std::size_t>> neighbours(nodes);
    for (const TreeEdge& edge:branches) {
        if (edge.first_node>=nodes || (edge.second_node!=ground_node &&
            (edge.second_node>=nodes || edge.first_node==edge.second_node))) {
            throw std::invalid_argument("invalid moving-tree branch");
        }
        if (edge.second_node!=ground_node) {
            neighbours[edge.first_node].push_back(edge.second_node);
            neighbours[edge.second_node].push_back(edge.first_node);
        }
    }
    for (std::vector<std::size_t>& adjacent:neighbours) {
        std::sort(adjacent.begin(),adjacent.end());
        adjacent.erase(std::unique(adjacent.begin(),adjacent.end()),adjacent.end());
    }
    std::vector<std::size_t> parent(nodes,ground_node),order{0};
    std::vector<bool> visited(nodes,false);visited[0]=true;
    for (std::size_t position=0;position<order.size();++position) {
        const std::size_t node=order[position];
        for (std::size_t other:neighbours[node]) {
            if (other==parent[node]) { continue; }
            if (visited[other]) { throw std::invalid_argument("moving network contains a cycle"); }
            visited[other]=true;parent[other]=node;order.push_back(other);
        }
    }
    if (order.size()!=nodes) { throw std::invalid_argument("moving network is disconnected"); }
    branches_.assign(branches.begin(),branches.end());parent_=std::move(parent);order_=std::move(order);
    diagonal_.resize(nodes);coupling_.resize(nodes);inverse_.resize(nodes);stage_.resize(nodes);
    charge.assign(nodes,0.0);momentum.assign(branches.size(),0.0);
}
void MovingTreeReference::advance(std::span<const double> c1,std::span<const double> c2,
                                  std::span<const double> g1,std::span<const double> g2,double dt) noexcept {
    constexpr double a11=.25,a12=.25-0.28867513459481288225;
    constexpr double a21=.25+0.28867513459481288225,a22=.25;
    for (std::size_t node=0;node<charge.size();++node) {
        diagonal_[node]={c1[node],0,0,c2[node]};coupling_[node]={};
        stage_[node]={charge[node],charge[node]};
    }
    for (std::size_t branch=0;branch<branches_.size();++branch) {
        const TreeEdge& edge=branches_[branch];
        const double first=g1[branch],second=g2[branch],h2=dt*dt;
        // h^2 A diag(G1,G2) A: 2x2 stage coupling along this branch.
        const Matrix block{h2*(a11*first*a11+a12*second*a21),h2*(a11*first*a12+a12*second*a22),
                           h2*(a21*first*a11+a22*second*a21),h2*(a21*first*a12+a22*second*a22)};
        const Pair load{dt*(a11*first+a12*second)*momentum[branch],
                        dt*(a21*first+a22*second)*momentum[branch]};
        for (std::size_t k=0;k<4;++k) { diagonal_[edge.first_node][k]+=block[k]; }
        for (std::size_t k=0;k<2;++k) { stage_[edge.first_node][k]-=load[k]; }
        if (edge.second_node!=ground_node) {
            const std::size_t child=parent_[edge.first_node]==edge.second_node ? edge.first_node : edge.second_node;
            for (std::size_t k=0;k<4;++k) {
                diagonal_[edge.second_node][k]+=block[k];coupling_[child][k]-=block[k];
            }
            for (std::size_t k=0;k<2;++k) { stage_[edge.second_node][k]+=load[k]; }
        }
    }
    // Block tree elimination costs O(nodes+branches), including new metrics.
    for (std::size_t position=order_.size();position>0;--position) {
        const std::size_t node=order_[position-1],up=parent_[node];
        inverse_[node]=inverse(diagonal_[node]);
        if (up!=ground_node) {
            const Matrix multiplier=product(coupling_[node],inverse_[node]);
            const Matrix correction=product(multiplier,coupling_[node]);
            const Pair load=product(multiplier,stage_[node]);
            for (std::size_t k=0;k<4;++k) { diagonal_[up][k]-=correction[k]; }
            for (std::size_t k=0;k<2;++k) { stage_[up][k]-=load[k]; }
        }
    }
    for (std::size_t node:order_) {
        const std::size_t up=parent_[node];
        if (up!=ground_node) {
            const Pair load=product(coupling_[node],stage_[up]);
            for (std::size_t k=0;k<2;++k) { stage_[node][k]-=load[k]; }
        }
        stage_[node]=product(inverse_[node],stage_[node]);
    }
    for (std::size_t branch=0;branch<branches_.size();++branch) {
        const TreeEdge& edge=branches_[branch];
        Pair pressure=stage_[edge.first_node];
        if (edge.second_node!=ground_node) {
            for (std::size_t k=0;k<2;++k) { pressure[k]-=stage_[edge.second_node][k]; }
        }
        const double first=momentum[branch]+dt*(a11*pressure[0]+a12*pressure[1]);
        const double second=momentum[branch]+dt*(a21*pressure[0]+a22*pressure[1]);
        const double transfer=dt*.5*(g1[branch]*first+g2[branch]*second);
        charge[edge.first_node]-=transfer;
        if (edge.second_node!=ground_node) { charge[edge.second_node]+=transfer; }
        momentum[branch]+=dt*.5*(pressure[0]+pressure[1]);
    }
}
}
