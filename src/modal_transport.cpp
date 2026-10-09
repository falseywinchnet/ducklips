#include "ducklips/modal_transport.hpp"
#include "real_dot.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
void ModalTransportReference::prepare(std::size_t modes,std::size_t nodes,std::size_t branches) {
    if (modes==0 || modes>nodes || nodes>4096 || branches>16384) {
        throw std::invalid_argument("invalid transport dimensions");
    }
    modes_=modes;nodes_=nodes;branches_=branches;
    omega_.assign(modes,0);frequency_rate_.assign(modes,0);next_rate_.resize(modes);
    pressure_.assign(modes*modes,0);flow_.assign(modes*modes,0);
    next_pressure_.resize(modes*modes);next_flow_.resize(modes*modes);
    differences_.resize(modes*branches);weighted_differences_.resize(modes*branches);
    weighted_pressure_.resize(modes*nodes);
}
ModalMotionSample ModalTransportReference::view() const noexcept { return {omega_,pressure_,flow_}; }
std::span<const double> ModalTransportReference::frequency_rate() const noexcept { return frequency_rate_; }
bool ModalTransportReference::derive(const AcousticNetwork& network,const ModalBasis& basis,
                                     std::span<const double> compliance_rate,
                                     std::span<const double> conductance_rate) noexcept {
    if (modes_==0 || basis.mode_count!=modes_ || basis.node_count!=nodes_ || network.nodes.size()!=nodes_ ||
        network.branches.size()!=branches_ || compliance_rate.size()!=nodes_ || conductance_rate.size()!=branches_ ||
        basis.angular_frequency.size()<modes_ || basis.pressure_shape.size()<modes_*nodes_) { return false; }
    for (double rate:compliance_rate) { if (!std::isfinite(rate)) { return false; } }
    for (double rate:conductance_rate) { if (!std::isfinite(rate)) { return false; } }
    double maximum=0;
    for (std::size_t i=0;i<modes_;++i) {
        const double omega=basis.angular_frequency[i];
        if (!std::isfinite(omega) || omega<=0 || !std::isfinite(omega*omega)) { return false; }
        maximum=std::max(maximum,omega*omega);
    }
    for (std::size_t i=0;i<modes_;++i) {
        for (std::size_t j=i+1;j<modes_;++j) {
            const double gap=basis.angular_frequency[j]*basis.angular_frequency[j]-basis.angular_frequency[i]*basis.angular_frequency[i];
            if (gap<=1e-8*maximum) { return false; }
        }
        const double* phi=basis.pressure_shape.data()+i*nodes_;
        for (std::size_t node=0;node<nodes_;++node) {
            if (!std::isfinite(phi[node])) { return false; }
            weighted_pressure_[i*nodes_+node]=phi[node]*compliance_rate[node];
        }
        for (std::size_t branch=0;branch<branches_;++branch) {
            const AcousticBranch& edge=network.branches[branch];
            if (edge.from>=nodes_ || (edge.to!=ground_node && edge.to>=nodes_)) { return false; }
            const double difference=phi[edge.from]-(edge.to==ground_node ? 0 : phi[edge.to]);
            differences_[i*branches_+branch]=difference;
            weighted_differences_[i*branches_+branch]=difference*conductance_rate[branch];
        }
    }
    for (std::size_t i=0;i<modes_;++i) {
        const double wi=basis.angular_frequency[i],li=wi*wi;
        for (std::size_t j=i;j<modes_;++j) {
            const double wj=basis.angular_frequency[j],lj=wj*wj;
            const double c=detail::real_dot(basis.pressure_shape.data()+i*nodes_,weighted_pressure_.data()+j*nodes_,nodes_);
            const double g=detail::real_dot(differences_.data()+i*branches_,weighted_differences_.data()+j*branches_,branches_);
            // D=U^T S_dot U, F=U^T diag(C_dot/C) U/2 + U^T U_dot.
            const double derivative=g-.5*(li+lj)*c;
            const double connection=i==j ? 0 : derivative/(lj-li);
            const double rate=i==j ? derivative/(2*li) : 0;
            const double pij=-.5*c-connection,pji=-.5*c+connection;
            const double fij=wi/wj*(.5*c-connection)+rate;
            const double fji=wj/wi*(.5*c+connection)+rate;
            if (!std::isfinite(pij) || !std::isfinite(pji) || !std::isfinite(fij) || !std::isfinite(fji)) { return false; }
            next_pressure_[i*modes_+j]=pij;next_pressure_[j*modes_+i]=pji;
            next_flow_[i*modes_+j]=fij;next_flow_[j*modes_+i]=fji;
            if (i==j) {
                next_rate_[i]=rate*wi;
                if (!std::isfinite(next_rate_[i])) { return false; }
            }
        }
    }
    std::copy_n(basis.angular_frequency.begin(),modes_,omega_.begin());
    pressure_.swap(next_pressure_);flow_.swap(next_flow_);frequency_rate_.swap(next_rate_);
    return true;
}
}
