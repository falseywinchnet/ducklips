#include "ducklips/contour_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace ducklips {
namespace {
std::array<double,4> positions(const FoldState& s){return {s.lower_position[0],s.lower_position[1],s.upper_position[0],s.upper_position[1]};}
std::array<double,4> velocities(const FoldState& s){return {s.lower_velocity[0],s.lower_velocity[1],s.upper_velocity[0],s.upper_velocity[1]};}
template<std::size_t N> bool finite(const std::array<double,N>& a){for(double x:a){if(!std::isfinite(x)){return false;}}return true;}
bool solve(std::array<double,49> a,std::array<double,7> b,std::array<double,7>& x){
    for(std::size_t j=0;j<7;++j){std::size_t p=j;for(std::size_t i=j+1;i<7;++i){if(std::abs(a[7*i+j])>std::abs(a[7*p+j])){p=i;}}
        if(!std::isfinite(a[7*p+j]) || a[7*p+j]==0){return false;}
        for(std::size_t k=j;k<7;++k){std::swap(a[7*j+k],a[7*p+k]);}std::swap(b[j],b[p]);
        for(std::size_t i=j+1;i<7;++i){const double f=a[7*i+j]/a[7*j+j];for(std::size_t k=j+1;k<7;++k){a[7*i+k]-=f*a[7*j+k];}b[i]-=f*b[j];}}
    for(std::size_t j=7;j-->0;){double v=b[j];for(std::size_t k=j+1;k<7;++k){v-=a[7*j+k]*x[k];}x[j]=v/a[7*j+j];}return finite(x);
}
}
bool ContourGlottalStorage::prepare(const FoldParameters& f,const AirProperties& air,double tension,
    const FoldState& initial,const std::array<double,2>& pressure,double outlet_loss) noexcept {
    // Validate body parameters with the existing contract. Its rectangular
    // positivity restriction does not apply to this retained wedge geometry.
    if(!finite(positions(initial)) || !finite(velocities(initial)) || !std::isfinite(initial.rest_area)){return false;}
    FoldState validation_state;validation_state.rest_area=3e-6;
    GlottalStorage validation;if(!validation.prepare(f,air,tension,validation_state,{})){return false;}
    ContourGlottalStorage candidate;if(!candidate.gas_.configure(air) || !candidate.kinetics_.configure(f,air,outlet_loss) || !candidate.contact_.configure(f,tension)){return false;}
    candidate.width_=f.length;candidate.state.folds=initial;
    for(std::size_t side=0;side<2;++side){const double q=tension*(side?f.asymmetry:1.);
        candidate.coupling_[side]=f.coupling_stiffness*q;
        candidate.tissue_mass_[side]=f.lower_mass/q;candidate.tissue_mass_[side+2]=f.upper_mass/q;
        candidate.stiffness_[side]=f.lower_stiffness*q;candidate.stiffness_[side+2]=f.upper_stiffness*q;
        candidate.damping_[side]=2*f.lower_damping_ratio*std::sqrt(candidate.tissue_mass_[side]*candidate.stiffness_[side]);
        candidate.damping_[side+2]=2*f.upper_damping_ratio*std::sqrt(candidate.tissue_mass_[side+2]*candidate.stiffness_[side+2]);}
    const auto x=positions(initial);const GlottalContour contour{f.lower_thickness,f.upper_thickness,
        initial.rest_area+f.length*(x[0]+x[1]),initial.rest_area+f.length*(x[2]+x[3])};GlottalContourGeometry geometry;
    if(!contour.geometry(geometry)){return false;}
    for(std::size_t i=0;i<2;++i){if(!candidate.gas_.mass_at_pressure(geometry.chamber[i].volume,pressure[i],candidate.state.mass[i])){return false;}}
    double e=0;if(!candidate.energy(e)){return false;}*this=candidate;return true;
}
bool ContourGlottalStorage::configuration(std::array<double,2>& area,ContourKineticState& result) const noexcept {
    const auto x=positions(state.folds);const std::array<double,2> a={state.folds.rest_area+width_*(x[0]+x[1]),state.folds.rest_area+width_*(x[2]+x[3])};
    if(!kinetics_.evaluate(state.mass,a,momentum,result)){return false;}area=a;return true;
}
double ContourGlottalStorage::mechanical_energy(const std::array<double,4>& x,const std::array<double,4>& v) const noexcept {
    double e=0;for(std::size_t i=0;i<4;++i){e+=.5*(tissue_mass_[i]*v[i]*v[i]+stiffness_[i]*x[i]*x[i]);}
    for(std::size_t i=0;i<2;++i){const double d=x[i]-x[i+2];e+=.5*coupling_[i]*d*d;}return e;
}
bool ContourGlottalStorage::energy(double& result) const noexcept {
    ContourKineticState kinetic;std::array<double,2> area;if(!configuration(area,kinetic)){return false;}
    FoldContactState contact;if(!contact_.evaluate(area,contact)){return false;}
    double sum=mechanical_energy(positions(state.folds),velocities(state.folds))+kinetic.energy+contact.energy;
    for(std::size_t i=0;i<2;++i){GasEnergy gas;if(!gas_.evaluate(state.mass[i],kinetic.geometry.chamber[i].volume,gas)){return false;}sum+=gas.energy;}
    if(!std::isfinite(sum)){return false;}result=sum;return true;
}
bool ContourGlottalStorage::advance(double dt,double reservoir_pressure,double free_load,double load_slope,
    double rest_area,const std::array<double,4>& external_force,ContourUpdate update) noexcept {
    if(!std::isfinite(dt) || dt<=0 || !std::isfinite(free_load) || !std::isfinite(load_slope) || load_slope<0
        || !std::isfinite(rest_area) || !finite(external_force)){return false;}
    double e0=0,reservoir_mass=0;GasEnergy reservoir;ContourKineticState k0;std::array<double,2> a0;
    if(!energy(e0) || !configuration(a0,k0) || !gas_.mass_at_pressure(1.,reservoir_pressure,reservoir_mass)
        || !gas_.evaluate(reservoir_mass,1.,reservoir)){return false;}
    const bool dissipative=update==ContourUpdate::dissipative;
    const auto x0=positions(state.folds),v0=velocities(state.folds);
    std::array<double,7> y{},residual{},tolerance{};
    for(std::size_t i=0;i<4;++i){y[3+i]=dissipative?dt*v0[i]:-.5*(rest_area-state.folds.rest_area)/width_;}
    std::array<double,2> m1{},a1{},area_effort{};std::array<double,3> p1{},flow{};std::array<GasTransfer,2> gas;
    FoldContactTransfer contact;std::array<bool,3> active{};std::array<double,3> ports{};
    ContourKineticState k1,mid;ContourKineticTransfer kinetic;double norm=0;bool converged=false;
    const auto evaluate=[&](const std::array<double,7>& trial){
        if(!finite(trial)){return false;}
        // Accumulate area displacement directly to avoid subtracting fold positions.
        for(std::size_t i=0;i<2;++i){a1[i]=a0[i]+rest_area-state.folds.rest_area+width_*(trial[3+2*i]+trial[4+2*i]);}
        active=dissipative?std::array<bool,3>{a1[0]>0,true,a1[1]>0}:std::array<bool,3>{k0.geometry.passage[0].open,k0.geometry.passage[1].open,k0.geometry.passage[2].open};
        for(std::size_t j=0;j<3;++j){flow[j]=active[j]?trial[j]:0.;}
        m1={state.mass[0]+dt*(flow[0]-flow[1]),state.mass[1]+dt*(flow[1]-flow[2])};
        if(!kinetics_.evaluate(m1,a1,{},k1)){return false;}
        for(std::size_t j=0;j<3;++j){
            if(!dissipative && k0.geometry.passage[j].open!=k1.geometry.passage[j].open){return false;}
            p1[j]=active[j]?(dissipative?flow[j]/k1.inverse_inertance[j]:4*flow[j]/(k0.inverse_inertance[j]+k1.inverse_inertance[j])-momentum[j]):0.;
        }
        if(!kinetics_.transfer(state.mass,a0,momentum,m1,a1,p1,kinetic,update)){return false;}
        std::array<double,2> mm,am;for(std::size_t i=0;i<2;++i){mm[i]=.5*(state.mass[i]+m1[i]);am[i]=.5*(a0[i]+a1[i]);
            if(!gas_.transfer(state.mass[i],k0.geometry.chamber[i].volume,m1[i],k1.geometry.chamber[i].volume,gas[i])){return false;}}
        if(dissipative){mid=k1;}else if(!kinetics_.evaluate(mm,am,{},mid)){return false;}
        if(!contact_.transfer(a0,a1,contact)){return false;}
        for(std::size_t k=0;k<2;++k){area_effort[k]=kinetic.effort[k+2]+contact.area_effort[k];
            // Contact clips the volume nonlinearly. Use the finite volume
            // gradient belonging to the same contour transfer as the inertia.
            for(std::size_t i=0;i<2;++i){area_effort[k]-=gas[i].pressure*kinetic.volume_gradient[i][k];}}
        const std::array<double,2> mu={gas[0].mass_effort+kinetic.effort[0],gas[1].mass_effort+kinetic.effort[1]};
        ports={mu[0]-reservoir.mass_effort,mu[1]-mu[0],free_load+load_slope*flow[2]-mu[1]};
        norm=0;converged=true;
        for(std::size_t j=0;j<3;++j){
            if(!active[j]){residual[j]=trial[j];tolerance[j]=1e-20;norm=std::max(norm,std::abs(trial[j])/1e-4);if(std::abs(trial[j])>tolerance[j]){converged=false;}continue;}
            const double inertia=(p1[j]-momentum[j])/dt,loss=(mid.resistance[j]+mid.quadratic_loss[j]*std::abs(flow[j]))*flow[j];
            residual[j]=inertia+loss+ports[j];tolerance[j]=1e-8+1e-11*(std::abs(inertia)+std::abs(loss)+std::abs(ports[j]));
            norm=std::max(norm,std::abs(residual[j])/1000.);if(std::abs(residual[j])>tolerance[j]){converged=false;}}
        for(std::size_t i=0;i<4;++i){const std::size_t side=i%2,other=i<2?i+2:i-2;const double speed=trial[3+i]/dt;
            const double inertia=2*tissue_mass_[i]*(speed-v0[i])/dt,spring=stiffness_[i]*(x0[i]+.5*trial[3+i]),
                coupled=coupling_[side]*(x0[i]-x0[other]+.5*(trial[3+i]-trial[3+other])),drag=damping_[i]*speed,
                fluid=width_*area_effort[i/2];
            residual[3+i]=inertia+spring+coupled+drag+fluid-external_force[i];
            const double subtraction=(2*tissue_mass_[i]/(dt*dt)+damping_[i]/dt)*std::abs(trial[3+i])+2*tissue_mass_[i]*std::abs(v0[i])/dt;
            tolerance[3+i]=1e-12+1e-10*(std::abs(inertia)+std::abs(spring)+std::abs(coupled)+std::abs(drag)+std::abs(fluid)+std::abs(external_force[i]))+8*std::numeric_limits<double>::epsilon()*subtraction;
            norm=std::max(norm,std::abs(residual[3+i])/.01);if(std::abs(residual[3+i])>tolerance[3+i]){converged=false;}}
        return finite(residual) && finite(tolerance) && std::isfinite(norm);
    };
    std::size_t iteration=0;
    for(;iteration<48;++iteration){if(!evaluate(y)){
            if(iteration!=0 || !dissipative){return false;}
            y.fill(0.);for(std::size_t i=0;i<4;++i){y[3+i]=-.5*(rest_area-state.folds.rest_area)/width_;}if(!evaluate(y)){return false;}
        }if(converged){break;}
        const auto base=residual;const double old_norm=norm;std::array<double,49> jacobian;
        for(std::size_t j=0;j<7;++j){double h=1e-5*(j<3?std::max({std::abs(y[j]),.01*std::min(state.mass[0],state.mass[1])/dt,1e-8}):std::max({std::abs(y[j]),std::abs(v0[j-3]*dt),1e-8}));
            bool differentiated=false;
            for(int attempt=0;attempt<16;++attempt){auto plus=y,minus=y;plus[j]+=h;minus[j]-=h;const bool okp=evaluate(plus);const auto rp=residual;const bool okm=evaluate(minus);const auto rm=residual;
                if(okp && okm){for(std::size_t i=0;i<7;++i){jacobian[7*i+j]=(rp[i]-rm[i])/(2*h);}differentiated=true;break;}h*=.5;}
            if(!differentiated){return false;}}
        auto rhs=base;for(double& r:rhs){r=-r;}std::array<double,7> delta;if(!solve(jacobian,rhs,delta)){return false;}
        bool accepted=false;double fraction=1.;for(int line=0;line<48;++line){auto next=y;for(std::size_t j=0;j<7;++j){next[j]+=fraction*delta[j];}
            if(evaluate(next) && (norm<old_norm || converged)){y=next;accepted=true;break;}fraction*=.5;}
        if(!accepted){return false;}}
    if(iteration==48){return false;}
    ContourGlottalStorage next=*this;next.state.mass=m1;next.momentum=p1;next.state.folds.rest_area=rest_area;auto x1=x0,v1=v0;ContourGlottalStep receipt;
    for(std::size_t i=0;i<4;++i){x1[i]+=y[3+i];v1[i]=2*y[3+i]/dt-v0[i];receipt.force_residual[i]=residual[3+i];
        receipt.external_work+=external_force[i]*y[3+i];receipt.damping_loss+=damping_[i]*y[3+i]*y[3+i]/dt;}
    next.state.folds.lower_position={x1[0],x1[1]};next.state.folds.upper_position={x1[2],x1[3]};next.state.folds.lower_velocity={v1[0],v1[1]};next.state.folds.upper_velocity={v1[2],v1[3]};
    for(std::size_t j=0;j<3;++j){receipt.flow[j]=flow[j];receipt.law_residual[j]=active[j]?residual[j]:0.;if(!active[j]){receipt.constraint_effort[j]=-(p1[j]-momentum[j])/dt-ports[j];}receipt.passage_loss+=dt*(mid.resistance[j]+mid.quadratic_loss[j]*std::abs(flow[j]))*flow[j]*flow[j];}
    for(std::size_t i=0;i<2;++i){receipt.pressure[i]=gas[i].pressure;receipt.actuator_work+=area_effort[i]*(rest_area-state.folds.rest_area);}
    receipt.kinetic_effort=kinetic.effort;receipt.area_effort=area_effort;receipt.contact_area_effort=contact.area_effort;receipt.reservoir_work=dt*reservoir.mass_effort*flow[0];receipt.load_work=dt*(free_load+load_slope*flow[2])*flow[2];
    double e1=0;if(!next.energy(e1)){return false;}receipt.energy_change=e1-e0;
    receipt.numerical_loss=kinetic.numerical_loss;
    receipt.energy_residual=receipt.energy_change+receipt.damping_loss+receipt.passage_loss+receipt.numerical_loss+receipt.load_work-receipt.reservoir_work-receipt.actuator_work-receipt.external_work;
    receipt.iterations=iteration+1;if(!std::isfinite(receipt.energy_residual)){return false;}next.last=receipt;*this=next;return true;
}
void GraphContourGlottal::prepare(std::size_t nodes,std::span<const TreeEdge> branches,std::size_t source_node,
    const FoldParameters& f,const AirProperties& air,double tension,const FoldState& initial,const std::array<double,2>& pressure,double loss){
    if(source_node>=nodes){throw std::invalid_argument("contour source outside graph");}GraphContourGlottal candidate;
    if(!candidate.source.prepare(f,air,tension,initial,pressure,loss)){throw std::invalid_argument("invalid contour source");}
    candidate.network.prepare(nodes,branches);candidate.free_.prepare(nodes,branches);candidate.response_.prepare(nodes,branches);
    candidate.unit_flow_.assign(nodes,0.);candidate.unit_flow_[source_node]=1.;candidate.source_node_=source_node;*this=std::move(candidate);
}
bool GraphContourGlottal::advance(std::span<const double> c1,std::span<const double> c2,std::span<const double> g1,std::span<const double> g2,
    double dt,double reservoir_pressure,double rest_area,const std::array<double,4>& external_force,std::span<const double> nodal_flow,std::span<const double> series_pressure,ContourUpdate update) noexcept {
    if(unit_flow_.empty() || network.charge.size()!=free_.charge.size() || network.momentum.size()!=free_.momentum.size()
        || network.average_pressure.size()!=free_.average_pressure.size() || network.average_flow.size()!=free_.average_flow.size()){return false;}
    std::copy(network.charge.begin(),network.charge.end(),free_.charge.begin());std::copy(network.momentum.begin(),network.momentum.end(),free_.momentum.begin());
    std::fill(response_.charge.begin(),response_.charge.end(),0.);std::fill(response_.momentum.begin(),response_.momentum.end(),0.);
    if(!free_.advance(c1,c2,g1,g2,dt,nodal_flow,series_pressure) || !response_.advance(c1,c2,g1,g2,dt,unit_flow_)){return false;}
    auto candidate=source;if(!candidate.advance(dt,reservoir_pressure,free_.average_pressure[source_node_],response_.average_pressure[source_node_],rest_area,external_force,update)){return false;}
    const double q=candidate.last.flow[2];
    const auto combine=[&](std::vector<double>& a,const std::vector<double>& b){for(std::size_t i=0;i<a.size();++i){a[i]+=q*b[i];if(!std::isfinite(a[i])){return false;}}return true;};
    if(!combine(free_.charge,response_.charge) || !combine(free_.momentum,response_.momentum) || !combine(free_.average_pressure,response_.average_pressure) || !combine(free_.average_flow,response_.average_flow)){return false;}
    std::copy(free_.charge.begin(),free_.charge.end(),network.charge.begin());std::copy(free_.momentum.begin(),free_.momentum.end(),network.momentum.begin());
    std::copy(free_.average_pressure.begin(),free_.average_pressure.end(),network.average_pressure.begin());std::copy(free_.average_flow.begin(),free_.average_flow.end(),network.average_flow.begin());source=candidate;return true;
}
}
