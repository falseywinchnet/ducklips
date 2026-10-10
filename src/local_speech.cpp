#include "ducklips/local_speech.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
void LocalSpeech::prepare(const AcousticNetwork& layout,const AirProperties& air,double wall_compliance,double glottal_length) {
    if(!std::isfinite(wall_compliance) || wall_compliance<=0){throw std::invalid_argument("local wall compliance must be positive");}
    if(!std::isfinite(glottal_length) || glottal_length<=0){throw std::invalid_argument("glottal length must be positive");}
    air_=air;wall_compliance_=wall_compliance;glottal_length_=glottal_length;
    branches_=layout.branches;
    const std::size_t n=layout.nodes.size(),b=branches_.size();
    charge.resize(n);inverse_c_.resize(n);volume_.resize(n);
    momentum.resize(b);wall_charge.resize(b);radiation_momentum.resize(b);
    g_.resize(b);wall_inverse_c_.resize(b);old_flow_.resize(b);
    transfer_p_.resize(b);transfer_d_.resize(b);momentum_p_.resize(b);momentum_d_.resize(b);
    radiation_h_.resize(b);radiation_decay_.resize(b);reset();
    jet_coefficient_.resize(b);half_step_.resize(b);
}
void LocalSpeech::reset() noexcept {
    for(std::vector<double>* state:{&charge,&momentum,&wall_charge,&radiation_momentum,&g_,&old_flow_,&volume_}){
        std::fill(state->begin(),state->end(),0.);
    }
    std::fill(inverse_c_.begin(),inverse_c_.end(),1.);
    std::fill(wall_inverse_c_.begin(),wall_inverse_c_.end(),0.);
    std::fill(radiation_h_.begin(),radiation_h_.end(),0.);
    glottal_flow_=0;glottal_momentum_=0;glottal_conductance_=0;ready_=false;
}
void LocalSpeech::exchange(std::size_t b,double series) noexcept {
    const AcousticBranch& edge=branches_[b];
    if(g_[b]==0){return;}
    // The exterior R || M load exchanges momentum locally. Its decaying
    // difference and conserved total are known directly, even for tiny ports.
    const double h=radiation_h_[b];
    if(h>0){
        const double sum=momentum[b]+radiation_momentum[b];
        const double equilibrium=h*sum/(g_[b]+h);
        momentum[b]=equilibrium+(momentum[b]-equilibrium)*radiation_decay_[b];
        radiation_momentum[b]=sum-momentum[b];
    }
    const double downstream=edge.to!=ground_node ? pressure(edge.to) : wall_inverse_c_[b]*wall_charge[b];
    const double difference=pressure(edge.from)-downstream+series;
    double transferred=transfer_p_[b]*momentum[b]+transfer_d_[b]*difference;
    if(jet_coefficient_[b]>0){
        // Vented jet / abrupt expansion: pressure loss beta Q |Q|. Include
        // it in the same local exchange, rather than applying a drag filter
        // after pressure has already pushed too much air through the opening.
        const double linear=transferred/half_step_[b];
        const double mid_flow=2*linear/(1+std::sqrt(1+4*jet_coefficient_[b]*std::abs(linear)));
        transferred=half_step_[b]*mid_flow;
        momentum[b]=2*mid_flow/g_[b]-momentum[b];
    }else{momentum[b]=momentum_p_[b]*momentum[b]+momentum_d_[b]*difference;}
    charge[edge.from]-=transferred;
    if(edge.to!=ground_node){charge[edge.to]+=transferred;}
    else if(edge.kind==BranchKind::wall){wall_charge[b]+=transferred;}
    if(h>0){
        const double sum=momentum[b]+radiation_momentum[b];
        const double equilibrium=h*sum/(g_[b]+h);
        momentum[b]=equilibrium+(momentum[b]-equilibrium)*radiation_decay_[b];
        radiation_momentum[b]=sum-momentum[b];
    }
}
void LocalSpeech::inject(double area,double lung,double h) noexcept {
    // The air crossing the folds has inertia as well as a separated-jet
    // pressure drop. Store its momentum instead of treating the opening as
    // an instantaneous resistor charging a tiny isolated inlet volume.
    if(area<=1e-12){glottal_momentum_=0;glottal_conductance_=0;return;}
    const double drop=lung-pressure();
    const double g=area/(air_.density*glottal_length_);
    glottal_conductance_=g;
    const double denominator=1+.25*h*h*g*inverse_c_[0];
    const double linear=g*(glottal_momentum_+.5*h*drop)/denominator;
    const double beta=.5*air_.density/(area*area);
    const double coefficient=.5*h*g*beta/denominator;
    const double flow=2*linear/(1+std::sqrt(1+4*coefficient*std::abs(linear)));
    glottal_momentum_=2*flow/g-glottal_momentum_;
    const double transferred=h*flow;
    charge[0]+=transferred;glottal_flow_+=transferred;
}
void LocalSpeech::advance(const AcousticNetwork& first,const AcousticNetwork& last,
                          double fraction,double increment,double dt,std::size_t substeps,
                          double area,double lung,std::size_t aspiration_branch,
                          std::size_t frication_branch,std::size_t jet_branch,double aspiration,double frication) noexcept {
    const double middle=fraction+.5*increment,end=fraction+increment;
    const double bulk=air_.density*air_.sound_speed*air_.sound_speed;
    const double h=dt/double(substeps),half=.5*h;
    double total_wall=0;
    for(std::size_t n=0;n<charge.size();++n){
        total_wall+=std::lerp(first.nodes[n].wall_area,last.nodes[n].wall_area,middle);
    }
    for(std::size_t b=0;b<branches_.size();++b){old_flow_[b]=flow(b);}
    for(std::size_t n=0;n<charge.size();++n){
        const double c=std::lerp(first.nodes[n].compliance,last.nodes[n].compliance,middle);
        const double v=bulk*std::lerp(first.nodes[n].compliance,last.nodes[n].compliance,end);
        // Charge is excess air volume: reducing the cavity displaces air,
        // while changing its compliance does not reset its pressure history.
        if(ready_){charge[n]-=.5*(v-volume_[n]);}
        inverse_c_[n]=1/c;
    }
    for(std::size_t b=0;b<branches_.size();++b){
        const AcousticBranch& a=first.branches[b];const AcousticBranch& z=last.branches[b];
        const double original_g=std::lerp(a.conductance,z.conductance,middle);
        double resistance=std::lerp(a.resistance,z.resistance,middle);
        const double node_area=std::lerp(first.nodes[a.from].area,last.nodes[a.from].area,middle);
        double inertia=original_g>0 ? 1/original_g : 0;
        radiation_h_[b]=0;wall_inverse_c_[b]=0;
        jet_coefficient_[b]=0;half_step_[b]=half;
        if(a.kind==BranchKind::wall){
            const double surface=std::lerp(first.nodes[a.from].wall_area,last.nodes[a.from].wall_area,middle);
            wall_inverse_c_[b]=total_wall/(wall_compliance_*surface);
        }else if(original_g>0){
            // Circular Poiseuille drag. Frequency-dependent boundary-layer
            // memory is deliberately not claimed by this first local engine.
            double aperture=node_area;
            if(a.to!=ground_node){aperture=std::min(aperture,std::lerp(first.nodes[a.to].area,last.nodes[a.to].area,middle));}
            resistance+=8*pi*air_.viscosity*inertia/(air_.density*aperture);
            // A smooth tube's changes of section exchange pressure and
            // momentum; they are not a succession of separated turbulent
            // jets. Separation loss belongs at the constriction and exits.
            if(a.to==ground_node || b==jet_branch){jet_coefficient_[b]=.5*air_.density/(aperture*aperture);}
            const double radiation=std::lerp(a.radiation_coefficient,z.radiation_coefficient,middle);
            if(radiation>0){
                const double exterior=air_.density*last.losses.end_correction/std::sqrt(pi*node_area);
                inertia=std::max(inertia-exterior,1e-12);
                radiation_h_[b]=1/exterior;
                const double r=exterior*exterior/radiation;
                radiation_decay_[b]=std::exp(-r*(1/inertia+1/exterior)*half*.5);
            }
        }
        g_[b]=original_g>0 ? 1/inertia : 0;
        if(g_[b]==0){momentum[b]=0;radiation_momentum[b]=0;continue;}
        const double inverse_sum=inverse_c_[a.from]+(a.to!=ground_node ? inverse_c_[a.to] : wall_inverse_c_[b]);
        // Exact algebra of this local midpoint exchange. All coefficients are
        // computed once per sample; each sweep then uses multiply-adds only.
        const double denominator=1+.5*half*resistance*g_[b]+.25*half*half*g_[b]*inverse_sum;
        jet_coefficient_[b]*=.5*half*g_[b]/denominator;
        transfer_p_[b]=half*g_[b]/denominator;
        transfer_d_[b]=.5*half*transfer_p_[b];
        momentum_p_[b]=2/denominator-1;
        momentum_d_[b]=half/denominator;
    }
    glottal_flow_=0;
    for(std::size_t step=0;step<substeps;++step){
        inject(area,lung,half);
        for(std::size_t b=0;b<branches_.size();++b){
            exchange(b,(b==aspiration_branch?aspiration:0)+(b==frication_branch?frication:0));
        }
        for(std::size_t b=branches_.size();b-->0;){
            exchange(b,(b==aspiration_branch?aspiration:0)+(b==frication_branch?frication:0));
        }
        inject(area,lung,half);
    }
    glottal_flow_/=dt;
    for(std::size_t n=0;n<charge.size();++n){
        const double c=std::lerp(first.nodes[n].compliance,last.nodes[n].compliance,end);
        const double v=bulk*c;
        if(ready_){charge[n]-=.5*(v-volume_[n]);}
        volume_[n]=v;inverse_c_[n]=1/c;
    }
    // Observe the physical endpoint, including exact closure. The midpoint
    // coefficients used above belong only to the time integration interval.
    for(std::size_t b=0;b<branches_.size();++b){
        const AcousticBranch& a=first.branches[b];const AcousticBranch& z=last.branches[b];
        const double conductance=std::lerp(a.conductance,z.conductance,end);
        if(conductance==0){g_[b]=0;momentum[b]=0;radiation_momentum[b]=0;radiation_h_[b]=0;continue;}
        double inertia=1/conductance;
        if(std::lerp(a.radiation_coefficient,z.radiation_coefficient,end)>0){
            const double aperture=std::lerp(first.nodes[a.from].area,last.nodes[a.from].area,end);
            const double exterior=air_.density*last.losses.end_correction/std::sqrt(pi*aperture);
            inertia=std::max(inertia-exterior,1e-12);radiation_h_[b]=1/exterior;
        }
        g_[b]=1/inertia;
    }
    ready_=true;
}
double LocalSpeech::pressure(std::size_t node) const noexcept {return charge[node]*inverse_c_[node];}
double LocalSpeech::flow(std::size_t branch) const noexcept {return g_[branch]*momentum[branch];}
double LocalSpeech::flow_derivative(std::size_t branch,double dt) const noexcept {return (flow(branch)-old_flow_[branch])/dt;}
double LocalSpeech::energy() const noexcept {
    double sum=0;
    for(std::size_t n=0;n<charge.size();++n){sum+=charge[n]*charge[n]*inverse_c_[n];}
    for(std::size_t b=0;b<momentum.size();++b){sum+=g_[b]*momentum[b]*momentum[b]+wall_inverse_c_[b]*wall_charge[b]*wall_charge[b]+radiation_h_[b]*radiation_momentum[b]*radiation_momentum[b];}
    return .5*sum+glottal_kinetic_energy();
}
}
