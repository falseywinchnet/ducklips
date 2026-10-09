#include "ducklips/glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ducklips {
namespace {
std::array<double,4> positions(const FoldState& s) noexcept {return {s.lower_position[0],s.lower_position[1],s.upper_position[0],s.upper_position[1]};}
std::array<double,4> velocities(const FoldState& s) noexcept {return {s.lower_velocity[0],s.lower_velocity[1],s.upper_velocity[0],s.upper_velocity[1]};}
bool finite(const std::array<double,4>& a) noexcept {for(double v:a){if(!std::isfinite(v)){return false;}}return true;}
bool solve(std::array<double,16> a,std::array<double,4> b,std::array<double,4>& x) noexcept {
    for(std::size_t j=0;j<4;++j){
        std::size_t pivot=j;for(std::size_t i=j+1;i<4;++i){if(std::abs(a[4*i+j])>std::abs(a[4*pivot+j])){pivot=i;}}
        if(a[4*pivot+j]==0 || !std::isfinite(a[4*pivot+j])){return false;}
        for(std::size_t k=j;k<4;++k){std::swap(a[4*j+k],a[4*pivot+k]);}std::swap(b[j],b[pivot]);
        for(std::size_t i=j+1;i<4;++i){const double f=a[4*i+j]/a[4*j+j];
            for(std::size_t k=j+1;k<4;++k){a[4*i+k]-=f*a[4*j+k];}b[i]-=f*b[j];}
    }
    for(std::size_t j=4;j-->0;){double sum=b[j];for(std::size_t k=j+1;k<4;++k){sum-=a[4*j+k]*x[k];}x[j]=sum/a[4*j+j];}
    return finite(x);
}
}
bool GlottalStorage::prepare(const FoldParameters& folds,const AirProperties& air,double tension,
    const FoldState& initial,const std::array<double,2>& pressure) noexcept {
    if(folds.strips!=1 || !std::isfinite(tension) || tension<=0 || !std::isfinite(folds.asymmetry) || folds.asymmetry<=0
        || !finite(positions(initial)) || !finite(velocities(initial)) || !std::isfinite(initial.rest_area)){return false;}
    for(double x:{folds.length,folds.lower_thickness,folds.upper_thickness,folds.lower_mass,folds.upper_mass,folds.lower_stiffness,folds.upper_stiffness}){
        if(!std::isfinite(x) || x<=0){return false;}
    }
    for(double x:{folds.coupling_stiffness,folds.lower_damping_ratio,folds.upper_damping_ratio}){if(!std::isfinite(x) || x<0){return false;}}
    GlottalStorage candidate;
    if(!candidate.gas_.configure(air)){return false;}
    candidate.length_=folds.length;candidate.thickness_={folds.lower_thickness,folds.upper_thickness};candidate.state.folds=initial;
    for(std::size_t i=0;i<2;++i){candidate.wall_area_[i]=candidate.length_*candidate.thickness_[i];}
    for(std::size_t side=0;side<2;++side){const double q=tension*(side?folds.asymmetry:1.);
        candidate.coupling_[side]=folds.coupling_stiffness*q;
        candidate.tissue_mass_[side]=folds.lower_mass/q;candidate.tissue_mass_[side+2]=folds.upper_mass/q;
        candidate.stiffness_[side]=folds.lower_stiffness*q;candidate.stiffness_[side+2]=folds.upper_stiffness*q;
        candidate.damping_[side]=2*folds.lower_damping_ratio*std::sqrt(candidate.tissue_mass_[side]*candidate.stiffness_[side]);
        candidate.damping_[side+2]=2*folds.upper_damping_ratio*std::sqrt(candidate.tissue_mass_[side+2]*candidate.stiffness_[side+2]);
    }
    if(!finite(candidate.tissue_mass_) || !finite(candidate.stiffness_) || !finite(candidate.damping_)){return false;}
    std::array<double,2> v;
    if(!candidate.volumes(v)){return false;}
    for(std::size_t i=0;i<2;++i){if(!candidate.gas_.mass_at_pressure(v[i],pressure[i],candidate.state.mass[i])){return false;}}
    double e=0;if(!candidate.energy(e)){return false;}*this=candidate;return true;
}
bool GlottalStorage::volumes(std::array<double,2>& result) const noexcept {
    const auto x=positions(state.folds);std::array<double,2> v;
    if(!finite(x) || !std::isfinite(state.folds.rest_area) || length_<=0){return false;}
    for(std::size_t i=0;i<2;++i){v[i]=thickness_[i]*(state.folds.rest_area+length_*(x[2*i]+x[2*i+1]));if(!std::isfinite(v[i]) || v[i]<=0){return false;}}
    result=v;return true;
}
double GlottalStorage::mechanical_energy(const std::array<double,4>& x,const std::array<double,4>& v) const noexcept {
    double sum=0;for(std::size_t i=0;i<4;++i){sum+=.5*(tissue_mass_[i]*v[i]*v[i]+stiffness_[i]*x[i]*x[i]);}
    for(std::size_t i=0;i<2;++i){const double d=x[i]-x[i+2];sum+=.5*coupling_[i]*d*d;}return sum;
}
bool GlottalStorage::energy(double& result) const noexcept {
    std::array<double,2> v;if(!volumes(v)){return false;}
    const auto velocity=velocities(state.folds);if(!finite(velocity)){return false;}
    double sum=mechanical_energy(positions(state.folds),velocity);
    for(std::size_t i=0;i<2;++i){GasEnergy e;if(!gas_.evaluate(state.mass[i],v[i],e)){return false;}sum+=e.energy;}
    if(!std::isfinite(sum)){return false;}result=sum;return true;
}
bool GlottalStorage::advance(const std::array<double,3>& flow,double rest_area,
    const std::array<double,4>& external_force,double dt) noexcept {
    if(!std::isfinite(dt) || dt<=0 || !std::isfinite(rest_area) || !finite(external_force)){return false;}
    for(double q:flow){if(!std::isfinite(q)){return false;}}
    std::array<double,2> v0;if(!volumes(v0)){return false;}double energy0=0;if(!energy(energy0)){return false;}
    const std::array<double,2> m1={state.mass[0]+dt*(flow[0]-flow[1]),state.mass[1]+dt*(flow[1]-flow[2])};
    for(double m:m1){if(!std::isfinite(m) || m<=0){return false;}}
    const auto x0=positions(state.folds),velocity0=velocities(state.folds);
    // Solve the small displacement directly, rather than recovering it by
    // subtracting two large positions at every Newton/velocity update.
    std::array<double,4> displacement{};
    // Feasible initial guess keeps the old volume even if the rest opening
    // moves. This is a Newton starting point, not a volume floor or a state reset.
    for(double& x:displacement){x-=.5*(rest_area-state.folds.rest_area)/length_;}
    std::array<GasTransfer,2> effort;std::array<double,2> v1;
    std::array<double,4> residual{},delta{};std::array<double,16> jacobian{};
    double norm=0,scale=0,roundoff=0;
    const auto evaluate=[&](const std::array<double,4>& step){
        auto x=x0;for(std::size_t i=0;i<4;++i){x[i]+=step[i];}
        if(!finite(step) || !finite(x)){return false;}norm=0;scale=0;roundoff=0;jacobian.fill(0.);
        for(std::size_t i=0;i<2;++i){
            v1[i]=thickness_[i]*(rest_area+length_*(x[2*i]+x[2*i+1]));
            if(!gas_.transfer(state.mass[i],v0[i],m1[i],v1[i],effort[i])){return false;}
        }
        for(std::size_t i=0;i<4;++i){
            const std::size_t cell=i/2,side=i%2,other=i<2?i+2:i-2;
            const double speed=step[i]/dt;
            const double inertial=2*tissue_mass_[i]*(speed-velocity0[i])/dt;
            const double spring=.5*stiffness_[i]*(x[i]+x0[i]);
            const double coupled=.5*coupling_[side]*(x[i]+x0[i]-x[other]-x0[other]);
            const double drag=damping_[i]*speed,force=wall_area_[cell]*effort[cell].pressure;
            residual[i]=inertial+spring+coupled+drag-force-external_force[i];
            norm=std::max(norm,std::abs(residual[i]));
            scale=std::max(scale,std::abs(inertial)+std::abs(spring)+std::abs(coupled)+std::abs(drag)+std::abs(force)+std::abs(external_force[i]));
            // Subtracting close velocities has a representable
            // force floor. Keep it separate from physical-force convergence:
            // scaling tolerance by both large opposing inertial terms lets
            // accepted force error grow as the time step shrinks.
            const double subtraction=(2*tissue_mass_[i]/(dt*dt)+damping_[i]/dt)*std::abs(step[i])
                +2*tissue_mass_[i]*std::abs(velocity0[i])/dt;
            roundoff=std::max(roundoff,8*std::numeric_limits<double>::epsilon()*subtraction);
            jacobian[4*i+i]=2*tissue_mass_[i]/(dt*dt)+.5*(stiffness_[i]+coupling_[side])+damping_[i]/dt;
            jacobian[4*i+other]=-.5*coupling_[side];
            for(std::size_t j=2*cell;j<2*cell+2;++j){jacobian[4*i+j]-=wall_area_[cell]*wall_area_[cell]*effort[cell].pressure_volume_slope;}
        }
        return finite(residual) && std::isfinite(norm) && std::isfinite(scale) && std::isfinite(roundoff);
    };
    std::size_t iteration=0;
    for(;iteration<64;++iteration){
        if(!evaluate(displacement)){return false;}
        if(norm<=1e-10*scale+roundoff){break;}
        auto rhs=residual;for(double& v:rhs){v=-v;}if(!solve(jacobian,rhs,delta)){return false;}
        const double old_norm=norm;double fraction=1.;bool accepted=false;
        for(int line=0;line<64;++line){auto trial=displacement;for(std::size_t i=0;i<4;++i){trial[i]+=fraction*delta[i];}
            if(evaluate(trial) && (norm<old_norm || norm<=1e-10*scale+roundoff)){displacement=trial;accepted=true;break;}fraction*=.5;
        }
        if(!accepted){return false;}
    }
    if(iteration==64){return false;}
    auto x1=x0;for(std::size_t i=0;i<4;++i){x1[i]+=displacement[i];}
    GlottalStorageState next=state;next.mass=m1;next.folds.rest_area=rest_area;
    std::array<double,4> velocity1{};GlottalStorageStep receipt;
    for(std::size_t i=0;i<4;++i){
        const double dx=displacement[i],speed=dx/dt;velocity1[i]=2*speed-velocity0[i];
        receipt.external_work+=external_force[i]*dx;receipt.damping_loss+=dt*damping_[i]*speed*speed;
        receipt.wall_work+=effort[i/2].pressure*wall_area_[i/2]*dx;
    }
    next.folds.lower_position={x1[0],x1[1]};next.folds.upper_position={x1[2],x1[3]};
    next.folds.lower_velocity={velocity1[0],velocity1[1]};next.folds.upper_velocity={velocity1[2],velocity1[3]};
    double energy1=mechanical_energy(x1,velocity1);
    for(std::size_t i=0;i<2;++i){GasEnergy e;if(!gas_.evaluate(m1[i],v1[i],e)){return false;}energy1+=e.energy;
        receipt.average_pressure[i]=effort[i].pressure;receipt.average_mass_effort[i]=effort[i].mass_effort;
        receipt.actuator_work-=effort[i].pressure*thickness_[i]*(rest_area-state.folds.rest_area);
    }
    receipt.flow_effort={effort[0].mass_effort,effort[1].mass_effort-effort[0].mass_effort,-effort[1].mass_effort};
    for(std::size_t i=0;i<3;++i){receipt.port_work+=dt*receipt.flow_effort[i]*flow[i];}
    receipt.energy_change=energy1-energy0;
    receipt.energy_residual=receipt.energy_change+receipt.damping_loss-receipt.port_work-receipt.actuator_work-receipt.external_work;
    receipt.iterations=iteration+1;
    if(!finite(velocity1) || !std::isfinite(energy1) || !std::isfinite(receipt.energy_residual)){return false;}
    state=next;last=receipt;return true;
}
}
