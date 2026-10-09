#include "ducklips/canonical_speech.hpp"
#include <algorithm>
#include <cmath>
namespace ducklips {
void CanonicalSpeech::prepare(const AcousticNetwork& layout,const AirProperties& air) {
    density_=air.density;
    const double low=2*pi*1e-4,high=2*pi*1e10;
    const double spacing=std::log(high/low)/(2*loss_terms);
    for (std::size_t term=0;term<loss_terms;++term) {
        loss_poles_[term]=low*std::exp((double(term)+.5)*2*spacing);
        loss_weights_[term]=2*std::sqrt(2.0)/pi*std::sqrt(loss_poles_[term])*spacing;
    }
    const std::size_t channels=layout.node_count()+layout.branch_count();
    physical_memory_.assign(channels*(loss_terms+1),0);response_memory_.resize(physical_memory_.size());
    loss_scale_.resize(channels);static_resistance_.resize(channels);rad_pole_.resize(channels);rad_weight_.resize(channels);
    std::vector<TreeEdge> edges;
    for (const AcousticBranch& branch:layout.branches) { edges.push_back({branch.from,branch.to}); }
    physical_.prepare(layout.node_count(),edges);response_.prepare(layout.node_count(),edges);
    c1_.resize(layout.node_count());c2_.resize(layout.node_count());cend_.resize(layout.node_count());
    g1_.resize(edges.size());g2_.resize(edges.size());gend_.resize(edges.size());old_flow_.resize(edges.size());
    nodal_.assign(layout.node_count(),0);nodal_[0]=1;series_.resize(edges.size());reset();
}
void CanonicalSpeech::reset() noexcept {
    std::fill(physical_memory_.begin(),physical_memory_.end(),0);
    std::fill(physical_.charge.begin(),physical_.charge.end(),0);
    std::fill(physical_.momentum.begin(),physical_.momentum.end(),0);
    std::fill(cend_.begin(),cend_.end(),1);std::fill(gend_.begin(),gend_.end(),0);
}
void CanonicalSpeech::free_step(const AcousticNetwork& start,const AcousticNetwork& end,
                               double fraction,double increment,double dt,std::size_t substeps,
                               std::size_t aspiration_branch,std::size_t frication_branch,
                               double aspiration,double frication) noexcept {
    std::fill(response_memory_.begin(),response_memory_.end(),0);
    std::fill(response_.charge.begin(),response_.charge.end(),0);
    std::fill(response_.momentum.begin(),response_.momentum.end(),0);
    std::fill(series_.begin(),series_.end(),0);
    series_[aspiration_branch]+=aspiration;series_[frication_branch]+=frication;
    for (std::size_t b=0;b<g1_.size();++b) {
        const double gs=inertive_conductance(start,b),ge=inertive_conductance(end,b);
        const double g=gs+fraction*(ge-gs);
        old_flow_[b]=g*physical_.momentum[b];
    }
    load_pressure_=0;load_impedance_=0;
    const double h=dt/double(substeps);
    for (std::size_t step=0;step<substeps;++step) {
        const double t1=fraction+increment*(double(step)+.21132486540518711775)/double(substeps);
        const double t2=fraction+increment*(double(step)+.78867513459481288225)/double(substeps);
        const double tm=fraction+increment*(double(step)+.5)/double(substeps);
        for (std::size_t n=0;n<c1_.size();++n) {
            const double c=start.nodes[n].compliance,change=end.nodes[n].compliance-c;
            c1_[n]=c+t1*change;c2_[n]=c+t2*change;cend_[n]=c+tm*change;
            loss_scale_[n]=start.nodes[n].thermal_coefficient+tm*(end.nodes[n].thermal_coefficient-start.nodes[n].thermal_coefficient);
        }
        for (std::size_t b=0;b<g1_.size();++b) {
            const double g=inertive_conductance(start,b),change=inertive_conductance(end,b)-g;
            g1_[b]=g+t1*change;g2_[b]=g+t2*change;gend_[b]=g+tm*change;
            const std::size_t channel=c1_.size()+b;
            static_resistance_[channel]=start.branches[b].resistance+tm*(end.branches[b].resistance-start.branches[b].resistance);
            loss_scale_[channel]=start.branches[b].viscous_coefficient+tm*(end.branches[b].viscous_coefficient-start.branches[b].viscous_coefficient);
            const double radiation=start.branches[b].radiation_coefficient+tm*(end.branches[b].radiation_coefficient-start.branches[b].radiation_coefficient);
            rad_weight_[channel]=0;rad_pole_[channel]=1;
            if (radiation>0 && gend_[b]>0) {
                const std::size_t node=start.branches[b].from;
                const double area=start.nodes[node].area+tm*(end.nodes[node].area-start.nodes[node].area);
                const double mass=density_*end.losses.end_correction*std::sqrt(area/pi)/area;
                rad_pole_[channel]=mass/radiation;rad_weight_[channel]=mass*mass/radiation;
            }
        }
        losses(physical_,physical_memory_,h/2);losses(response_,response_memory_,h/2);
        physical_.advance(c1_,c2_,g1_,g2_,h,{},series_);
        response_.advance(c1_,c2_,g1_,g2_,h,nodal_);
        load_pressure_+=physical_.average_pressure[0]/double(substeps);
        load_impedance_+=response_.average_pressure[0]/double(substeps);
        losses(physical_,physical_memory_,h/2);losses(response_,response_memory_,h/2);
        close_passages(physical_,physical_memory_);close_passages(response_,response_memory_);
    }
    for (std::size_t n=0;n<cend_.size();++n) { cend_[n]=start.nodes[n].compliance+(fraction+increment)*(end.nodes[n].compliance-start.nodes[n].compliance); }
    for (std::size_t b=0;b<gend_.size();++b) { gend_[b]=inertive_conductance(start,b)+(fraction+increment)*(inertive_conductance(end,b)-inertive_conductance(start,b)); }
    close_passages(physical_,physical_memory_);close_passages(response_,response_memory_);
}
void CanonicalSpeech::close_passages(MovingTreeReference& state,std::vector<double>& memory) noexcept {
    // At zero area the air path no longer exists. Its inverse inertance is
    // zero, so a pressure difference must not integrate an unobservable
    // momentum and release it as flow on reopening. The branch's relaxation
    // states likewise belong to the absent air path, not to the surrounding
    // cavities. Keep all nodal pressure and all open-path history intact.
    for(std::size_t branch=0;branch<gend_.size();++branch){
        if(gend_[branch]!=0){continue;}
        state.momentum[branch]=0;
        const std::size_t offset=(cend_.size()+branch)*(loss_terms+1);
        std::fill_n(memory.begin()+static_cast<std::ptrdiff_t>(offset),loss_terms+1,0.);
    }
}
double CanonicalSpeech::inertive_conductance(const AcousticNetwork& network,std::size_t index) const noexcept {
    const AcousticBranch& branch=network.branches[index];
    if (branch.conductance==0 || branch.radiation_coefficient==0) { return branch.conductance; }
    const double area=network.nodes[branch.from].area;
    const double mass=density_*network.losses.end_correction*std::sqrt(area/pi)/area;
    return 1/(1/branch.conductance-mass);
}
void CanonicalSpeech::losses(MovingTreeReference& state,std::vector<double>& memory,double h) noexcept {
    // Each channel is an inertance/compliance connected to parallel passive
    // relaxation branches. Store their momenta/charges, not flow/pressure,
    // so changing coefficients do not reset their physical history.
    for (std::size_t channel=0;channel<loss_scale_.size();++channel) {
        const bool node=channel<cend_.size();
        const double metric=node ? 1/cend_[channel] : gend_[channel-cend_.size()];
        double& primary=node ? state.charge[channel] : state.momentum[channel-cend_.size()];
        const std::size_t offset=channel*(loss_terms+1);
        double resistance=static_resistance_[channel],load=0;
        for (std::size_t term=0;term<=loss_terms;++term) {
            const double pole=term==loss_terms ? rad_pole_[channel] : loss_poles_[term];
            const double weight=term==loss_terms ? rad_weight_[channel] : loss_scale_[channel]*loss_weights_[term];
            const double denominator=1+h*.5*pole;
            resistance+=weight/denominator;load+=pole*memory[offset+term]/denominator;
        }
        const double midpoint=(primary+h*.5*load)/(1+h*.5*metric*resistance);
        const double conjugate=metric*midpoint;
        for (std::size_t term=0;term<=loss_terms;++term) {
            const double pole=term==loss_terms ? rad_pole_[channel] : loss_poles_[term];
            const double weight=term==loss_terms ? rad_weight_[channel] : loss_scale_[channel]*loss_weights_[term];
            double& value=memory[offset+term];
            value+=h*(weight*conjugate-pole*value)/(1+h*.5*pole);
        }
        primary=2*midpoint-primary;
    }
}
void CanonicalSpeech::inject(double flow) noexcept {
    for (std::size_t i=0;i<physical_memory_.size();++i) { physical_memory_[i]+=flow*response_memory_[i]; }
    for (std::size_t n=0;n<cend_.size();++n) { physical_.charge[n]+=flow*response_.charge[n]; }
    for (std::size_t b=0;b<gend_.size();++b) { physical_.momentum[b]+=flow*response_.momentum[b]; }
}
double CanonicalSpeech::load_pressure() const noexcept { return load_pressure_; }
double CanonicalSpeech::load_impedance() const noexcept { return load_impedance_; }
double CanonicalSpeech::pressure() const noexcept { return physical_.charge[0]/cend_[0]; }
double CanonicalSpeech::flow_derivative(std::size_t branch,double dt) const noexcept {
    return (gend_[branch]*physical_.momentum[branch]-old_flow_[branch])/dt;
}
double CanonicalSpeech::energy() const noexcept {
    double total=0;
    for (std::size_t n=0;n<cend_.size();++n) { total+=physical_.charge[n]*physical_.charge[n]/cend_[n]; }
    for (std::size_t b=0;b<gend_.size();++b) { total+=gend_[b]*physical_.momentum[b]*physical_.momentum[b]; }
    for (std::size_t channel=0;channel<loss_scale_.size();++channel) {
        for (std::size_t term=0;term<=loss_terms;++term) {
            const double pole=term==loss_terms ? rad_pole_[channel] : loss_poles_[term];
            const double weight=term==loss_terms ? rad_weight_[channel] : loss_scale_[channel]*loss_weights_[term];
            const double value=physical_memory_[channel*(loss_terms+1)+term];
            if (weight>0) { total+=pole/weight*value*value; }
        }
    }
    return total/2;
}
}
