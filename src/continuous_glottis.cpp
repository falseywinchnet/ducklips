#include "ducklips/continuous_glottis.hpp"
#include "ducklips/glottal_relaxation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace ducklips {
namespace {
using Vector=std::array<double,7>;using Matrix=std::array<double,49>;
constexpr double energy_scale=1e-6,time_scale=1e-3,mass_scale=1e-8,area_scale=1e-6,position_scale=1e-4;
constexpr double coordinate_scale[]={mass_scale,mass_scale,mass_scale,position_scale,position_scale,position_scale,position_scale};
template<std::size_t N> bool finite(const std::array<double,N>& a){for(double v:a){if(!std::isfinite(v)){return false;}}return true;}
Vector multiply(const Matrix& a,const Vector& x){Vector y{};for(std::size_t i=0;i<7;++i){for(std::size_t j=0;j<7;++j){y[i]+=a[7*i+j]*x[j];}}return y;}
double dot(const Vector& a,const Vector& b){double s=0;for(std::size_t i=0;i<7;++i){s+=a[i]*b[i];}return s;}
bool solve(Matrix a,Vector b,Vector& x){
    for(std::size_t j=0;j<7;++j){std::size_t pivot=j;for(std::size_t i=j+1;i<7;++i){if(std::abs(a[7*i+j])>std::abs(a[7*pivot+j])){pivot=i;}}
        if(!std::isfinite(a[7*pivot+j]) || a[7*pivot+j]==0){return false;}
        for(std::size_t k=j;k<7;++k){std::swap(a[7*j+k],a[7*pivot+k]);}std::swap(b[j],b[pivot]);
        for(std::size_t i=j+1;i<7;++i){const double v=a[7*i+j]/a[7*j+j];for(std::size_t k=j+1;k<7;++k){a[7*i+k]-=v*a[7*j+k];}b[i]-=v*b[j];}}
    for(std::size_t j=7;j-->0;){double v=b[j];for(std::size_t k=j+1;k<7;++k){v-=a[7*j+k]*x[k];}x[j]=v/a[7*j+j];}return finite(x);
}
bool inverse(const Matrix& a,const std::array<bool,7>& active,Matrix& result){
    Matrix reduced=a;for(std::size_t j=0;j<7;++j){if(!active[j]){for(std::size_t k=0;k<7;++k){reduced[7*j+k]=reduced[7*k+j]=0;}reduced[7*j+j]=1;}}
    Matrix candidate{};for(std::size_t j=0;j<7;++j){if(!active[j]){continue;}Vector b{},x{};b[j]=1;if(!solve(reduced,b,x)){return false;}for(std::size_t i=0;i<7;++i){candidate[7*i+j]=active[i]?x[i]:0;}}
    for(std::size_t j=0;j<7;++j){for(std::size_t i=0;i<j;++i){candidate[7*i+j]=candidate[7*j+i]=.5*(candidate[7*i+j]+candidate[7*j+i]);}}
    result=candidate;return true;
}
Matrix transform(const std::array<double,25>& a,double width,double multiplier){
    Matrix result{};for(std::size_t i=0;i<7;++i){const std::size_t u=i<3?i:3+(i-3)/2;const double si=coordinate_scale[i]/time_scale*(i<3?1.:width);
        for(std::size_t j=0;j<7;++j){const std::size_t v=j<3?j:3+(j-3)/2;const double sj=coordinate_scale[j]/time_scale*(j<3?1.:width);result[7*i+j]=a[5*u+v]*si*sj*multiplier;}}
    return result;
}
}
std::array<double,4> ContinuousGlottis::coordinates(const Vector& q) const noexcept {
    return {initial_mass_[0]/mass_scale+q[0]-q[1],initial_mass_[1]/mass_scale+q[1]-q[2],
        (rest_area_+width_*position_scale*(q[3]+q[4]))/area_scale,(rest_area_+width_*position_scale*(q[5]+q[6]))/area_scale};
}
bool ContinuousGlottis::configuration(const std::array<double,4>& z,Configuration& result) const noexcept {
    if(!finite(z) || z[0]<=0 || z[1]<=0){return false;}
    Configuration candidate;candidate.z=z;candidate.active={z[2]>0,true,z[3]>0,true,true,true,true};const std::array<double,2> a={z[2]*area_scale,z[3]*area_scale},m={z[0]*mass_scale,z[1]*mass_scale};
    if(!storage_.recover(a,m,effort_guess_,candidate.gas) || !kinetics_.prepare(a,candidate.gas.effort,candidate.flow)){return false;}
    FoldContactState contact;if(!contact_.evaluate(a,contact)){return false;}
    candidate.potential=(candidate.gas.energy+contact.energy)/energy_scale;
    candidate.potential_gradient={candidate.gas.effort[0]*mass_scale/energy_scale,candidate.gas.effort[1]*mass_scale/energy_scale,
        (candidate.gas.area_effort[0]+contact.area_effort[0])*area_scale/energy_scale,(candidate.gas.area_effort[1]+contact.area_effort[1])*area_scale/energy_scale};
    candidate.metric=transform(candidate.flow.inertia,width_,1/energy_scale);
    // The same outlet transfer coordinate carries the attached half-column.
    // Its inertia belongs in the shared metric, not an independent flow state.
    candidate.metric[16]+=outlet_inertance_*mass_scale*mass_scale/(time_scale*time_scale*energy_scale);
    for(std::size_t i=0;i<4;++i){candidate.metric[7*(3+i)+3+i]+=mass_[i]*position_scale*position_scale/(time_scale*time_scale*energy_scale);}
    if(!inverse(candidate.metric,candidate.active,candidate.inverse)){return false;}result=candidate;return true;
}
bool ContinuousGlottis::resistance(const Configuration& c,Matrix& result) const noexcept {
    ContourPressureLoss loss;if(!loss_.prepare(c.flow,loss)){return false;}Matrix candidate=transform(loss.resistance,width_,time_scale/energy_scale);
    for(std::size_t i=0;i<4;++i){candidate[7*(3+i)+3+i]+=damping_[i]*position_scale*position_scale/(time_scale*energy_scale);}result=candidate;return true;
}
double ContinuousGlottis::mechanical(const Vector& q) const noexcept {
    double e=0;for(std::size_t i=0;i<4;++i){e+=.5*stiffness_[i]*q[3+i]*q[3+i];}
    for(std::size_t i=0;i<2;++i){const double d=q[3+i]-q[5+i];e+=.5*coupling_[i]*d*d;}return e*position_scale*position_scale/energy_scale;
}
bool ContinuousGlottis::prepare(const FoldParameters& f,const AirProperties& air,double tension,const FoldState& initial,const std::array<double,2>& effort,const std::array<double,3>& flow,double outlet_loss,ContourPressureProfile profile,double outlet_inertance) noexcept {
    for(double v:{f.lower_mass,f.upper_mass,f.lower_stiffness,f.upper_stiffness,f.length,f.asymmetry,tension}){if(!std::isfinite(v) || v<=0){return false;}}
    for(double v:{f.coupling_stiffness,f.lower_damping_ratio,f.upper_damping_ratio,outlet_loss,outlet_inertance}){if(!std::isfinite(v) || v<0){return false;}}
    if(!finite(effort) || !finite(flow) || !std::isfinite(initial.rest_area) || !finite(initial.lower_position) || !finite(initial.upper_position) || !finite(initial.lower_velocity) || !finite(initial.upper_velocity)){return false;}
    ContinuousGlottis candidate;if(!candidate.storage_.configure(f,air,profile) || !candidate.kinetics_.configure(f,air,profile) || !candidate.loss_.configure(f,air,profile) || !candidate.contact_.configure(f,tension) || !candidate.gas_.configure(air)){return false;}
    candidate.profile_=profile;candidate.rest_area_=initial.rest_area;candidate.width_=f.length;candidate.density_=air.density;candidate.bulk_=air.density*air.sound_speed*air.sound_speed;candidate.gamma_=air.heat_capacity_ratio;candidate.outlet_loss_=outlet_loss;candidate.outlet_inertance_=outlet_inertance;candidate.effort_guess_=effort;
    Vector velocity{};for(std::size_t j=0;j<3;++j){velocity[j]=flow[j]*time_scale/mass_scale;}
    for(std::size_t side=0;side<2;++side){const double q=tension*(side?f.asymmetry:1.);candidate.coupling_[side]=f.coupling_stiffness*q;
        candidate.mass_[side]=f.lower_mass/q;candidate.mass_[side+2]=f.upper_mass/q;candidate.stiffness_[side]=f.lower_stiffness*q;candidate.stiffness_[side+2]=f.upper_stiffness*q;
        candidate.damping_[side]=2*f.lower_damping_ratio*std::sqrt(candidate.mass_[side]*candidate.stiffness_[side]);candidate.damping_[side+2]=2*f.upper_damping_ratio*std::sqrt(candidate.mass_[side+2]*candidate.stiffness_[side+2]);
        candidate.position_[3+side]=initial.lower_position[side]/position_scale;candidate.position_[5+side]=initial.upper_position[side]/position_scale;
        velocity[3+side]=initial.lower_velocity[side]*time_scale/position_scale;velocity[5+side]=initial.upper_velocity[side]*time_scale/position_scale;}
    auto z=candidate.coordinates(candidate.position_);const std::array<double,2> a={z[2]*area_scale,z[3]*area_scale};ContourPressureState gas;if(!candidate.storage_.evaluate(a,effort,gas)){return false;}
    candidate.initial_mass_=gas.mass;candidate.active_={a[0]>0,true,a[1]>0,true,true,true,true};for(std::size_t j:{0u,2u}){if(!candidate.active_[j] && flow[j]!=0){return false;}}
    Configuration c;if(!candidate.configuration(candidate.coordinates(candidate.position_),c)){return false;}candidate.momentum_=multiply(c.metric,velocity);ContinuousGlottalState checked;if(!candidate.state(checked)){return false;}*this=candidate;return true;
}
bool ContinuousGlottis::state(ContinuousGlottalState& result) const noexcept {
    Configuration c;if(!configuration(coordinates(position_),c)){return false;}const Vector velocity=multiply(c.inverse,momentum_);ContinuousGlottalState candidate;candidate.folds.rest_area=rest_area_;candidate.mass=c.gas.mass;candidate.effort=c.gas.effort;candidate.area=c.flow.area;
    for(std::size_t j=0;j<3;++j){candidate.flow[j]=velocity[j]*mass_scale/time_scale;}
    for(std::size_t j=0;j<2;++j){candidate.folds.lower_position[j]=position_[3+j]*position_scale;candidate.folds.upper_position[j]=position_[5+j]*position_scale;candidate.folds.lower_velocity[j]=velocity[3+j]*position_scale/time_scale;candidate.folds.upper_velocity[j]=velocity[5+j]*position_scale/time_scale;}
    const Vector canonical=multiply(c.metric,velocity);for(std::size_t j=0;j<7;++j){candidate.momentum[j]=canonical[j]*energy_scale*time_scale/coordinate_scale[j];}
    candidate.outlet_inertance=outlet_inertance_;candidate.energy=energy_scale*(c.potential+mechanical(position_)+.5*dot(momentum_,velocity));if(!std::isfinite(candidate.energy)){return false;}result=candidate;return true;
}
bool ContinuousGlottis::advance(double dt,double reservoir_pressure,double free_load,double load_slope,const std::array<double,4>& force,ContinuousGlottalUpdate update,double outlet_inertance) noexcept {
    if(!std::isfinite(outlet_inertance) || (outlet_inertance<0 && outlet_inertance!=-1.)){return false;}
    if(outlet_inertance>=0 && outlet_inertance!=outlet_inertance_){
        // Prescribed geometry changes Hamiltonian coefficients, not canonical
        // state. Symmetric midpoint freezing is second order for smooth motion.
        ContinuousGlottis candidate=*this;ContinuousGlottalState before,middle,after;
        if(!state(before)){return false;}
        candidate.outlet_inertance_=.5*outlet_inertance_+.5*outlet_inertance;
        if(!candidate.state(middle)){return false;}
        double geometry_work=middle.energy-before.energy;
        if(!candidate.advance(dt,reservoir_pressure,free_load,load_slope,force,update)){return false;}
        if(!candidate.state(middle)){return false;}
        candidate.outlet_inertance_=outlet_inertance;
        if(!candidate.state(after)){return false;}
        geometry_work+=after.energy-middle.energy;
        candidate.last.geometry_work=geometry_work;candidate.last.energy_change=after.energy-before.energy;
        candidate.last.energy_residual=candidate.last.energy_change+candidate.last.numerical_loss+candidate.last.tissue_loss+candidate.last.fluid_loss+candidate.last.load_work-candidate.last.reservoir_work-candidate.last.external_work-geometry_work;
        if(!std::isfinite(candidate.last.energy_residual)){return false;}
        *this=candidate;return true;
    }
    if(!std::isfinite(dt) || dt<=0 || !std::isfinite(free_load) || !std::isfinite(load_slope) || load_slope<0 || !finite(force)){return false;}
    const bool dissipative=update==ContinuousGlottalUpdate::dissipative;
    const bool fitted=update==ContinuousGlottalUpdate::relaxation;
    if(update!=ContinuousGlottalUpdate::symmetric && !dissipative && !fitted){return false;}
    Configuration c0;if(!configuration(coordinates(position_),c0)){return false;}double reservoir_mass;GasEnergy reservoir;if(!gas_.mass_at_pressure(1.,reservoir_pressure,reservoir_mass) || !gas_.evaluate(reservoir_mass,1.,reservoir)){return false;}
    const double step=dt/time_scale;const Vector v0=multiply(c0.inverse,momentum_);Vector change{},residual{},p1{},q1{},speed{},external{};for(std::size_t j=0;j<7;++j){change[j]=step*v0[j];}
    external[0]=reservoir.mass_effort*mass_scale/energy_scale;external[2]=-free_load*mass_scale/energy_scale;for(std::size_t j=0;j<4;++j){external[3+j]=force[j]*position_scale/energy_scale;}
    if(dissipative && (c0.z[2]==0 || c0.z[3]==0) && dot(v0,v0)==0){
        // From rest at exact contact, the velocity predictor lies precisely
        // on a changing constraint. Predict displacement from the physical
        // pressure/spring acceleration to select a trial phase, not a gap floor.
        const auto& g=c0.potential_gradient;Vector rhs=external;rhs[0]-=g[0];rhs[1]-=g[1]-g[0];rhs[2]+=g[1];
        for(std::size_t i=0;i<4;++i){const std::size_t other=i<2?i+2:i-2;rhs[3+i]-=g[2+i/2]*width_*position_scale/area_scale+(stiffness_[i]*position_[3+i]+coupling_[i%2]*(position_[3+i]-position_[3+other]))*position_scale*position_scale/energy_scale;}
        const Vector acceleration=multiply(c0.inverse,rhs);for(std::size_t j=0;j<7;++j){change[j]+=.5*step*step*acceleration[j];}
    }
    Configuration c1,cm;Matrix drag{},relaxation{};double norm=0,load_loss=0,exit_loss=0,roundoff_scale=0;bool converged=false,roundoff_limited=false;
    const auto evaluate=[&](const Vector& trial){
        if(!finite(trial)){return false;}for(std::size_t j=0;j<7;++j){q1[j]=position_[j]+trial[j];speed[j]=trial[j]/step;}
        // Eliminate a closed coordinate exactly. Newton still solves its
        // algebraic residual, but a residual tolerance must not leak mass.
        const auto proposed=coordinates(q1);for(std::size_t j:{0u,2u}){if(proposed[j==0?2:3]<=0){q1[j]=position_[j];speed[j]=0;}}
        const auto z1=coordinates(q1);std::array<double,4> zm;for(std::size_t k=0;k<4;++k){zm[k]=.5*(c0.z[k]+z1[k]);}
        if(!configuration(z1,c1) || !configuration(zm,cm)){return false;}
        if((!dissipative || profile_==ContourPressureProfile::flat_caps) && (c1.active!=active_ || cm.active!=active_)){return false;}
        const Configuration& loss_configuration=dissipative?c1:cm;if(!resistance(loss_configuration,drag)){return false;}
        if(dissipative){p1=multiply(c1.metric,speed);}
        else {
        Matrix average;for(std::size_t j=0;j<49;++j){average[j]=.5*(c0.inverse[j]+c1.inverse[j]);}
        if(fitted){
            if(!glottal_relaxation_bias(average,drag,active_,step,relaxation)){return false;}
            Matrix implicit,explicit_part;for(std::size_t j=0;j<49;++j){implicit[j]=.5*(average[j]+relaxation[j]);explicit_part[j]=.5*(average[j]-relaxation[j]);}
            Vector right=multiply(explicit_part,momentum_);for(std::size_t j=0;j<7;++j){right[j]=speed[j]-right[j];if(!active_[j]){implicit[7*j+j]=1.;}}
            if(!solve(implicit,right,p1)){return false;}for(std::size_t j=0;j<7;++j){if(!active_[j]){p1[j]=0.;}}
        }else{
        for(std::size_t j=0;j<7;++j){if(!active_[j]){average[7*j+j]=1;}}
        Vector pavg;if(!solve(average,speed,pavg)){return false;}for(std::size_t j=0;j<7;++j){p1[j]=2*pavg[j]-momentum_[j];if(!active_[j]){p1[j]=0;}}
        }
        }
        // Symmetric discrete kinetic derivative in the four actual state
        // coordinates. Preserving this dependence retains the cyclic common
        // transferred-volume coordinate and its throughflow momentum law.
        std::array<double,4> gradient=cm.potential_gradient,dz;double norm2=0,work=0;
        for(std::size_t k=0;k<4;++k){
            // Gas moments can become much smaller than the nominal volume
            // scale as contact narrows the occupied support. A unit floor
            // then makes the capped difference perturb 10% of the remaining
            // gas, biasing the kinetic force. Differentiate positive gas
            // moments relatively; signed area keeps its contact-side rule.
            double h=2e-5*(k<2?std::abs(zm[k]):std::max(1.,std::abs(zm[k])));Matrix derivative;
            if(k>=2 && zm[k]==0){
                // At exact contact use the closed-side derivative. This is
                // also a defined Newton starting point for opening from rest.
                auto minus=zm,twice=zm;minus[k]-=h;twice[k]-=2*h;Configuration cn,c2;if(!configuration(minus,cn) || !configuration(twice,c2)){return false;}
                for(std::size_t j=0;j<49;++j){derivative[j]=(3*cm.inverse[j]-4*cn.inverse[j]+c2.inverse[j])/(2*h);}
            }else{
                h=std::min(h,.1*std::abs(zm[k]));if(h<=0){return false;}auto plus=zm,minus=zm;plus[k]+=h;minus[k]-=h;Configuration cp,cn;if(!configuration(plus,cp) || !configuration(minus,cn)){return false;}
                for(std::size_t j=0;j<49;++j){derivative[j]=(cp.inverse[j]-cn.inverse[j])/(2*h);}
            }
            gradient[k]+=dissipative?.5*dot(momentum_,multiply(derivative,momentum_)):.25*(dot(momentum_,multiply(derivative,momentum_))+dot(p1,multiply(derivative,p1)));dz[k]=z1[k]-c0.z[k];norm2+=dz[k]*dz[k];work+=gradient[k]*dz[k];}
        Matrix difference;for(std::size_t j=0;j<49;++j){difference[j]=c1.inverse[j]-c0.inverse[j];}
        const double delta=c1.potential-c0.potential+(dissipative?.5*dot(momentum_,multiply(difference,momentum_)):.25*(dot(momentum_,multiply(difference,momentum_))+dot(p1,multiply(difference,p1))));
        // Dividing cancellation in endpoint energies by a vanishing state
        // change creates a spurious force near steady throughflow. Correct
        // only a discrepancy resolved above the endpoint arithmetic scale.
        double arithmetic_scale=std::abs(c0.potential)+std::abs(c1.potential);
        for(std::size_t i=0;i<7;++i){for(std::size_t j=0;j<7;++j){const double inverse_scale=std::abs(c0.inverse[7*i+j])+std::abs(c1.inverse[7*i+j]);arithmetic_scale+=inverse_scale*(dissipative?.5*std::abs(momentum_[i]*momentum_[j]):.25*(std::abs(momentum_[i]*momentum_[j])+std::abs(p1[i]*p1[j])));}}
        const double energy_roundoff=64*std::numeric_limits<double>::epsilon()*arithmetic_scale;
        roundoff_scale=energy_roundoff;roundoff_limited=std::abs(delta-work)<=energy_roundoff;
        // Remove only the resolved excess above the arithmetic uncertainty.
        // A hard on/off correction jumps by energy_roundoff/|dz| at the
        // threshold and can leave Newton with no continuous residual to solve.
        // Soft thresholding leaves at most energy_roundoff uncorrected while
        // keeping this force continuous; the actual energy error is reported.
        if(norm2>1e-28 && !roundoff_limited){const double correction=std::copysign(std::abs(delta-work)-energy_roundoff,delta-work)/norm2;for(std::size_t k=0;k<4;++k){gradient[k]+=correction*dz[k];}}
        Vector generalized={gradient[0],gradient[1]-gradient[0],-gradient[1],0,0,0,0};
        for(std::size_t i=0;i<4;++i){const std::size_t other=i<2?i+2:i-2;const double x=.5*(position_[3+i]+q1[3+i]),y=.5*(position_[3+other]+q1[3+other]);
            generalized[3+i]=gradient[2+i/2]*width_*position_scale/area_scale+(stiffness_[i]*x+coupling_[i%2]*(x-y))*position_scale*position_scale/energy_scale;}
        Vector resistance_force=multiply(drag,speed);const double outflow=speed[2]*mass_scale/time_scale;double beta=0;
        if(loss_configuration.active[2] && outlet_loss_>0){const double r=std::pow(1+(gamma_-1)*loss_configuration.gas.effort[1]/bulk_,1/(gamma_-1));beta=outlet_loss_*density_/(2*r*r*loss_configuration.flow.area[1]*loss_configuration.flow.area[1]);}
        const double load=load_slope*outflow+beta*outflow*std::abs(outflow);resistance_force[2]+=load*mass_scale/energy_scale;load_loss=load*outflow*dt;
        exit_loss=beta*outflow*outflow*std::abs(outflow)*dt;
        norm=0;converged=true;
        for(std::size_t j=0;j<7;++j){
            residual[j]=(dissipative?c1.active[j]:active_[j])?(p1[j]-momentum_[j])/step+generalized[j]+resistance_force[j]-external[j]:trial[j]/step;
            norm=std::max(norm,std::abs(residual[j]));
            // Subtracting the two represented momenta limits the resolvable
            // acceleration as dt approaches zero. Do not demand a residual
            // below that arithmetic floor (the trajectory tolerance is
            // checked separately by coupled step doubling).
            const double momentum_roundoff=32*std::numeric_limits<double>::epsilon()*(std::abs(p1[j])+std::abs(momentum_[j]))/step;
            if(std::abs(residual[j])>2e-9+2e-11*(std::abs(generalized[j])+std::abs(resistance_force[j])+std::abs(external[j]))+momentum_roundoff){converged=false;}}
        return finite(residual) && std::isfinite(norm);
    };
    std::size_t iteration=0;for(;iteration<24;++iteration){if(!evaluate(change)){if(iteration!=0){return false;}change.fill(0);if(!evaluate(change)){return false;}}if(converged){break;}
        const Vector base=residual;const double old_norm=norm;const auto jacobian_active=c1.active;const auto jacobian_area=c1.flow.area;Matrix jacobian{};
        for(std::size_t j=0;j<7;++j){if(!jacobian_active[j]){jacobian[7*j+j]=1/step;continue;}double h=1e-5*std::max(std::abs(change[j]),1e-5);
            // A newly opening port may carry much less than the generic
            // flux perturbation. Scale its Jacobian step by the acoustic
            // volume that its actual opening can transport during this step.
            if(j==0 || j==2){h=std::min(h,1e-5*jacobian_area[j==0?0:1]*std::sqrt(bulk_/density_)*dt/mass_scale);}
            if(h<=0){return false;}bool differentiated=false;
            for(int attempt=0;attempt<12;++attempt){auto plus=change,minus=change;plus[j]+=h;minus[j]-=h;const bool okp=evaluate(plus);const Vector rp=residual;const bool okm=evaluate(minus);const Vector rm=residual;
                if(okp && okm){for(std::size_t i=0;i<7;++i){jacobian[7*i+j]=(rp[i]-rm[i])/(2*h);}differentiated=true;break;}h*=.5;}
            if(!differentiated){return false;}}
        Vector rhs=base,increment{};for(double& v:rhs){v=-v;}if(!solve(jacobian,rhs,increment)){return false;}double fraction=1;bool accepted=false;
        for(int line=0;line<24;++line){Vector next=change;for(std::size_t j=0;j<7;++j){next[j]+=fraction*increment[j];}if(evaluate(next) && (norm<old_norm || converged)){change=next;accepted=true;break;}fraction*=.5;}if(!accepted){return false;}
    }
    // The final permitted Newton correction may have converged in the
    // line search. Do not reject it solely because the loop budget ended.
    if(iteration==24 && !converged){return false;}for(std::size_t j=0;j<7;++j){if(!c1.active[j]){change[j]=0;}}ContinuousGlottis next=*this;next.position_=q1;next.momentum_=p1;next.effort_guess_=c1.gas.effort;next.active_=c1.active;
    // The three transferred volumes contain a cyclic throughflow coordinate.
    // Rebase their origins after each step, retaining the two actual gas
    // moments. Otherwise long throughflow subtracts large nearly equal
    // cumulative volumes and destroys pressure/Jacobian precision.
    next.initial_mass_={c1.z[0]*mass_scale,c1.z[1]*mass_scale};
    for(std::size_t j=0;j<3;++j){next.position_[j]=0.;}
    ContinuousGlottalStep receipt;
    for(std::size_t j=0;j<3;++j){receipt.flow[j]=speed[j]*mass_scale/time_scale;}
    receipt.reservoir_work=reservoir.mass_effort*change[0]*mass_scale;receipt.load_work=free_load*change[2]*mass_scale+load_loss;
    receipt.port_work=(free_load+load_slope*receipt.flow[2])*receipt.flow[2]*dt;receipt.outlet_loss=exit_loss;
    for(std::size_t j=0;j<4;++j){receipt.external_work+=force[j]*change[3+j]*position_scale;receipt.tissue_loss+=dt*damping_[j]*std::pow(speed[3+j]*position_scale/time_scale,2);}
    receipt.fluid_loss=step*energy_scale*dot(speed,multiply(drag,speed))-receipt.tissue_loss;
    const double e0=energy_scale*(c0.potential+mechanical(position_)+.5*dot(momentum_,v0));const double e1=energy_scale*(c1.potential+mechanical(q1)+.5*dot(p1,multiply(c1.inverse,p1)));
    if(dissipative){Vector delta_p;for(std::size_t j=0;j<7;++j){delta_p[j]=p1[j]-momentum_[j];}receipt.numerical_loss=.5*energy_scale*dot(delta_p,multiply(c1.inverse,delta_p));}
    else if(fitted){Vector delta_p;for(std::size_t j=0;j<7;++j){delta_p[j]=p1[j]-momentum_[j];}receipt.numerical_loss=.5*energy_scale*dot(delta_p,multiply(relaxation,delta_p));}
    receipt.energy_roundoff_scale=energy_scale*roundoff_scale;receipt.roundoff_limited_gradient=roundoff_limited;receipt.energy_change=e1-e0;receipt.energy_residual=receipt.energy_change+receipt.numerical_loss+receipt.tissue_loss+receipt.fluid_loss+receipt.load_work-receipt.reservoir_work-receipt.external_work;receipt.iterations=iteration+1;receipt.maximum_scaled_residual=norm;
    if(!std::isfinite(receipt.energy_residual)){return false;}next.last=receipt;*this=next;return true;
}
}
