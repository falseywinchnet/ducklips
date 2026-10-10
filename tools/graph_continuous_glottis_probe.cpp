#include "ducklips/graph_continuous_glottis.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>

namespace {std::atomic<std::size_t> allocations{};}
void* operator new(std::size_t size){++allocations;if(void* result=std::malloc(size)){return result;}throw std::bad_alloc();}
void operator delete(void* pointer) noexcept {std::free(pointer);}
void operator delete(void* pointer,std::size_t) noexcept {std::free(pointer);}

namespace {
using ducklips::GraphContinuousGlottis;
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
double number(){double value;require(std::scanf("%lf",&value)==1,"fixture number");return value;}
void read(std::vector<double>& values){for(double& v:values){v=number();}}
void print(std::span<const double> values){for(double v:values){std::printf(" %.17g",v);}}
double graph_energy(const ducklips::MovingGraphReference& graph,std::span<const double> c,std::span<const double> g){
    double result=0.;for(std::size_t i=0;i<c.size();++i){result+=.5*graph.charge[i]*graph.charge[i]/c[i];}
    for(std::size_t i=0;i<g.size();++i){result+=.5*g[i]*graph.momentum[i]*graph.momentum[i];}return result;
}
ducklips::FoldState folds(const std::array<double,9>& state){
    ducklips::FoldParameters parameters;ducklips::FoldState result;result.rest_area=3e-6;
    result.lower_position.fill((state[2]-result.rest_area)/(2*parameters.length));
    result.upper_position.fill((state[3]-result.rest_area)/(2*parameters.length));
    result.lower_velocity.fill(state[7]/(2*parameters.length));result.upper_velocity.fill(state[8]/(2*parameters.length));return result;
}
std::array<double,9> physical(const ducklips::ContinuousGlottalState& s){
    return {s.effort[0],s.effort[1],s.area[0],s.area[1],s.flow[0],s.flow[1],s.flow[2],
        .014*(s.folds.lower_velocity[0]+s.folds.lower_velocity[1]),.014*(s.folds.upper_velocity[0]+s.folds.upper_velocity[1])};
}
constexpr std::array<double,9> initial={112.18465038080178,-3.138675468892969,
    3.6868231667933133e-6,4.149026342131548e-6,.00017466261521028343,
    .00017466261521028343,.00017466261521028343,0.,0.};
void prepare(GraphContinuousGlottis& model,std::size_t nodes,std::span<const ducklips::TreeEdge> edges,
    std::size_t port,const std::array<double,9>& state=initial,double loss=1.,double column=0.){
    model.prepare(nodes,edges,port,{}, {},1.,folds(state),{state[0],state[1]},
        {state[4],state[5],state[6]},loss,ducklips::ContourPressureProfile::linear_span,column);
}
void fixture(bool column=false,bool radiation=false,bool thermal=false,bool viscosity=false){
    unsigned n,m,steps,port;require(std::scanf("%u%u%u%u",&n,&m,&steps,&port)==4,"fixture dimensions");
    require(n>0 && n<=64 && m<=256 && steps<=1000000,"fixture size");
    std::vector<ducklips::TreeEdge> edges(m);
    for(auto& edge:edges){int first,second;require(std::scanf("%d%d",&first,&second)==2 && first>=0 && second>=-1,"fixture edge");edge={static_cast<std::size_t>(first),second<0?ducklips::ground_node:static_cast<std::size_t>(second)};}
    std::array<double,9> y;for(double& v:y){v=number();}const double loss=number(),inertance=column?number():0.;GraphContinuousGlottis model;prepare(model,n,edges,port,y,loss,inertance);
    if(radiation){const double area=number();require(model.radiation.prepare({},area),"radiation prepare");model.radiation.momentum=number();}
    read(model.network.charge);read(model.network.momentum);
    std::vector<ducklips::ThermalGeometry> thermal_geometry;
    if(thermal){std::vector<ducklips::ThermalSite> sites;
        for(std::size_t i=0;i<n;++i){const double c=number(),gap=number();sites.push_back({i,c,gap});thermal_geometry.push_back({c,gap});}
        model.prepare_thermal({},sites);for(auto& site:model.thermal){for(double& z:site.charge){z=number();}}
    }
    std::vector<ducklips::ViscousGeometry> viscous_geometry;
    if(viscosity){unsigned count;require(std::scanf("%u",&count)==1 && count<=256,"viscous count");std::vector<ducklips::ViscousSite> sites;
        for(unsigned i=0;i<count;++i){const auto branch=static_cast<std::size_t>(number());const double l=number(),g=number(),a=number();sites.push_back({branch,{l,g,a}});viscous_geometry.push_back({l,g,a});}
        model.prepare_viscous({},sites);for(auto& site:model.viscous){for(double& z:site.momentum){z=number();}}
    }
    std::vector<double> c1(n),c2(n),g1(m),g2(m),nodal(n),series(m);
    for(unsigned step=0;step<steps;++step){
        const double dt=number(),pressure=number();const int mode=static_cast<int>(number());const double next_inertance=column?number():-1.;
        ducklips::RadiationGeometry boundary;boundary.branch=m-1;
        if(radiation){boundary.midpoint_area=number();boundary.last_area=number();boundary.first_conductance=number();boundary.last_conductance=number();}
        if(thermal){for(auto& site:thermal_geometry){site.compliance=number();site.half_gap=number();}}
        if(viscosity){for(auto& site:viscous_geometry){site.segment_inertance=number();site.branch_conductance=number();site.half_gap=number();}}
        read(c1);read(c2);read(g1);read(g2);read(nodal);read(series);
        require(mode>=0 && mode<=2,"fixture update mode");
        const ducklips::ContinuousGlottalUpdate updates[]={ducklips::ContinuousGlottalUpdate::symmetric,ducklips::ContinuousGlottalUpdate::dissipative,ducklips::ContinuousGlottalUpdate::relaxation};
        if(!model.advance(c1,c2,g1,g2,dt,pressure,{},nodal,series,updates[mode],next_inertance,radiation?&boundary:nullptr,thermal_geometry,viscous_geometry)){
            std::fprintf(stderr,"coupled fixture failed at step %u\n",step);throw std::runtime_error("fixture advance");
        }
        ducklips::ContinuousGlottalState state;require(model.source.state(state),"fixture state");
        const auto y1=physical(state);print(y1);print(model.network.charge);print(model.network.momentum);
        print(model.network.average_pressure);print(model.network.average_flow);
        const auto& s=model.source.last;const auto& r=model.last;
        const double ledger[]={state.energy,state.mass[0],state.mass[1],s.reservoir_work,s.external_work,s.fluid_loss,s.tissue_loss,
            s.outlet_loss,s.numerical_loss,s.energy_residual,r.graph_port_work,r.graph_external_work,r.port_residual,r.flow,r.pressure,r.free_pressure,r.load_slope,
            s.flow[0],s.load_work-s.port_work-s.outlet_loss};print(ledger);if(column){const double extra[]={s.geometry_work,state.outlet_inertance};print(extra);}if(radiation){const double extra[]={model.radiation.momentum,model.radiation.energy(),r.radiated_energy,r.radiation_geometry_work,r.radiation_energy_residual};print(extra);}
        if(thermal){double energy=0,mass=0;for(const auto& site:model.thermal){print(site.charge);energy+=site.energy();mass+=site.stored_air();}
            const double extra[]={energy,mass,r.thermal_dissipation,r.thermal_geometry_work,r.thermal_energy_residual};print(extra);
        }if(viscosity){double energy=0;for(const auto& site:model.viscous){print(site.momentum);energy+=site.energy();}
            const double extra[]={energy,r.viscous_dissipation,r.viscous_geometry_work,r.viscous_energy_residual};print(extra);
        }std::puts("");
    }
}
void radiation_controls(bool thermal=false,bool viscosity=false){
    const std::array<ducklips::TreeEdge,4> edges={{{0,1},{1,2},{2,0},{2,ducklips::ground_node}}};
    const std::array<double,3> c={1e-11,1.4e-11,8e-12};const std::array<double,4> g={.002,.003,.001,.0025};
    GraphContinuousGlottis model;prepare(model,3,edges,1,initial,1.,80.);
    require(model.radiation.prepare({},1e-4),"radiation prepare");model.radiation.momentum=.004;
    model.network.charge={1e-9,-2e-9,5e-10};model.network.momentum={.01,-.02,.015,.005};
    std::array<ducklips::ThermalGeometry,3> thermal_geometry{};
    std::array<ducklips::ViscousGeometry,3> viscous_geometry={{{.4/g[1],g[1],.001},{.6/g[1],g[1],.002},{.8/g[3],g[3],.0005}}};
    if(viscosity){const std::array<ducklips::ViscousSite,3> sites={{{1,viscous_geometry[0]},{1,viscous_geometry[1]},{3,viscous_geometry[2]}}};model.prepare_viscous({},sites);}
    if(thermal){const std::array<ducklips::ThermalSite,3> sites={{{0,c[0],.001},{1,c[1],.002},{2,c[2],.0005}}};model.prepare_thermal({},sites);}
    const auto thermal_energy=[&](){double value=0;for(const auto& t:model.thermal){value+=t.energy();}for(const auto& v:model.viscous){value+=v.energy();}return value;};
    const auto stored_air=[&](){double value=0;for(double q:model.network.charge){value+=q;}for(const auto& t:model.thermal){value+=t.stored_air();}return value;};
    ducklips::ContinuousGlottalState state;require(model.source.state(state),"source state");
    const double initial_energy=state.energy+graph_energy(model.network,c,g)+model.radiation.energy();
    double work=0,loss=0,worst=0,radiated=0,heat=0,mass_error=0,transfer=0;const double initial_mass=state.mass[0]+state.mass[1]+stored_air();
    for(int step=0;step<300;++step){
        const double area=1e-4*(1+.2*std::sin(2*ducklips::pi*(step+1)/300.));
        const ducklips::RadiationGeometry geometry{3,g[3],g[3],.5*(area+model.radiation.area()),area};
        for(std::size_t i=0;i<3;++i){thermal_geometry[i]={c[i],.001*(1+.2*std::sin(.01*(step+1)+i))};}
        const auto count=allocations.load();
        require(model.advance(c,c,g,g,2.5e-6,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::relaxation,-1.,&geometry,thermal?std::span<const ducklips::ThermalGeometry>(thermal_geometry):std::span<const ducklips::ThermalGeometry>{},viscosity?std::span<const ducklips::ViscousGeometry>(viscous_geometry):std::span<const ducklips::ViscousGeometry>{}),"radiation coupled advance");
        require(allocations.load()==count,"radiation advance allocated");require(model.source.state(state),"accepted source");
        const auto& s=model.source.last;const auto& r=model.last;
        work+=s.reservoir_work+s.external_work+s.geometry_work+r.graph_external_work+r.radiation_geometry_work+r.thermal_geometry_work+r.viscous_geometry_work;
        loss+=s.fluid_loss+s.tissue_loss+s.outlet_loss+s.numerical_loss+r.radiated_energy+r.thermal_dissipation+r.viscous_dissipation;radiated+=r.radiated_energy;heat+=r.thermal_dissipation;
        transfer+=2.5e-6*(s.flow[0]-model.network.average_flow[3]);mass_error=std::max(mass_error,std::abs(state.mass[0]+state.mass[1]+stored_air()-initial_mass-transfer));
        require(r.radiated_energy>=0,"radiation passive");
        worst=std::max(worst,std::abs(state.energy+graph_energy(model.network,c,g)+model.radiation.energy()+thermal_energy()-initial_energy+loss-work));
    }
    require(worst<1e-13 && radiated>0 && mass_error<1e-20 && (!thermal || heat>0),"joint radiation/thermal work and air");
    const auto saved=model;const auto physical_before=physical(state);
    const ducklips::RadiationGeometry bad{3,g[3],g[3],model.radiation.area(),-1.};
    require(!model.advance(c,c,g,g,2.5e-6,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::relaxation,-1.,&bad,thermal?std::span<const ducklips::ThermalGeometry>(thermal_geometry):std::span<const ducklips::ThermalGeometry>{},viscosity?std::span<const ducklips::ViscousGeometry>(viscous_geometry):std::span<const ducklips::ViscousGeometry>{}),"reject endpoint radiation after source trial");
    require(model.source.state(state) && physical(state)==physical_before && model.network.charge==saved.network.charge
        && model.network.momentum==saved.network.momentum && model.radiation.momentum==saved.radiation.momentum
        && model.radiation.area()==saved.radiation.area() && model.last.radiated_energy==saved.last.radiated_energy,"radiation transaction");
    if(thermal){
        for(std::size_t i=0;i<3;++i){require(model.thermal[i].charge==saved.thermal[i].charge,"thermal rollback");}
        thermal_geometry.back().compliance=-1.;
        require(!model.advance(c,c,g,g,2.5e-6,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::relaxation,-1.,nullptr,thermal_geometry,viscosity?std::span<const ducklips::ViscousGeometry>(viscous_geometry):std::span<const ducklips::ViscousGeometry>{}),"late thermal geometry rejection");
        require(model.source.state(state) && physical(state)==physical_before && model.network.charge==saved.network.charge,"late thermal rollback");
        for(std::size_t i=0;i<3;++i){require(model.thermal[i].charge==saved.thermal[i].charge,"all thermal states preserved");}
        std::printf("PASS thermal coupling: dissipated %.17g J, stored-air residual %.17g m3; ",heat,mass_error);
    }
    if(viscosity){
        for(std::size_t i=0;i<3;++i){require(model.viscous[i].momentum==saved.viscous[i].momentum,"viscous rollback");}
        for(std::size_t i=0;i<3;++i){thermal_geometry[i].compliance=c[i];}
        viscous_geometry[0].half_gap=-1;
        require(!model.advance(c,c,g,g,2.5e-6,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::relaxation,-1.,nullptr,thermal_geometry,viscous_geometry),"late viscous rejection");
        require(model.network.momentum==saved.network.momentum && model.source.state(state) && physical(state)==physical_before,"viscous transaction");
        for(std::size_t i=0;i<3;++i){require(model.viscous[i].momentum==saved.viscous[i].momentum,"all viscous memories preserved");}
        const std::array<ducklips::ViscousSite,2> duplicate={{{1,{.6/g[1],g[1],.001}},{1,{.6/g[1],g[1],.001}}}};
        bool rejected=false;try{model.prepare_viscous({},duplicate);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected && model.viscous.size()==saved.viscous.size(),"no duplicated plug inertia");
        auto sealed_g=g;sealed_g[1]=0.;viscous_geometry[0].half_gap=.001;
        const auto open_geometry=viscous_geometry;
        for(std::size_t i=0;i<2;++i){double metric_work;
            require(model.viscous[i].equilibrate(1e-4),"seed surviving profile");
            viscous_geometry[i].branch_conductance=0.;
            require(model.viscous[i].change_geometry(viscous_geometry[i].segment_inertance,0.,viscous_geometry[i].half_gap,metric_work),"block local through-flow");
        }
        model.network.momentum[1]=0.;const double surviving=model.viscous[1].energy();
        const auto count=allocations.load();
        require(model.advance(c,c,sealed_g,sealed_g,2.5e-6,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::relaxation,-1.,nullptr,thermal_geometry,viscous_geometry),"sealed viscous graph");
        require(allocations.load()==count && model.network.average_flow[1]==0,"sealed flow and allocation");
        require(model.viscous[1].energy()>0 && model.viscous[1].energy()<surviving,"surviving air relaxes behind seal");
        viscous_geometry[0]={0,0,0};
        require(model.advance(c,c,sealed_g,sealed_g,2.5e-6,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::relaxation,-1.,nullptr,thermal_geometry,viscous_geometry),"retire absent segment transaction");
        require(model.viscous[0].energy()==0 && model.viscous[1].energy()>0,"retire only vanished air");
        for(double z:model.viscous[0].momentum){require(z==0,"no absent history");}
        require(model.advance(c,c,g,g,2.5e-6,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::relaxation,-1.,nullptr,thermal_geometry,open_geometry),"reopen viscous graph");
        std::printf("PASS shared-branch viscous coupling; ");
    }
    std::printf("PASS radiating coupled source: work residual %.17g J, radiated %.17g J; allocation-free and rollback\n",worst,radiated);
}
void controls(ducklips::ContinuousGlottalUpdate update=ducklips::ContinuousGlottalUpdate::symmetric,bool column=false){
    const std::array<ducklips::TreeEdge,4> edges={{{0,1},{1,2},{2,0},{2,ducklips::ground_node}}};
    const std::array<double,4> c={1e-11,1.4e-11,8e-12,2e-11},g={.002,.003,.001,.0025};
    const std::array<double,4> nodal={0.,1e-6,0.,-2e-6},series={0.,30.,0.,-10.};
    GraphContinuousGlottis model;prepare(model,4,edges,1,initial,1.,column?80.:0.);
    model.network.charge={1e-9,-2e-9,5e-10,3e-9};model.network.momentum={.01,-.02,.015,.005};
    ducklips::ContinuousGlottalState before;require(model.source.state(before),"initial source");
    if(column){
        GraphContinuousGlottis bare;prepare(bare,4,edges,1);ducklips::ContinuousGlottalState base;
        require(bare.source.state(base),"bare state");
        require(std::abs(before.energy-base.energy-.5*80.*initial[6]*initial[6])<1e-18,"column stored kinetic energy");
    }
    const double initial_energy=before.energy+graph_energy(model.network,c,g);double work=0,loss=0,worst_energy=0,worst_mass=0,worst_direct=0,worst_port=0;
    double mass=before.mass[0]+before.mass[1];for(double v:model.network.charge){mass+=v;}double transfer=0;
    for(int step=0;step<300;++step){
        ducklips::MovingGraphReference direct=model.network;const double dt=2.5e-6;
        const std::size_t count=allocations.load();
        require(model.advance(c,c,g,g,dt,1400.,{},nodal,series,update,column?80.*(1.+.2*std::sin(2*ducklips::pi*(step+1)/300.)):-1.),"coupled fixed advance");
        require(allocations.load()==count,"coupled advance allocated");
        auto drives=nodal;drives[1]+=model.last.flow;require(direct.advance(c,c,g,g,dt,drives,series),"direct graph advance");
        for(std::size_t i=0;i<c.size();++i){worst_direct=std::max(worst_direct,std::abs(direct.charge[i]-model.network.charge[i])/c[i]);worst_direct=std::max(worst_direct,std::abs(direct.average_pressure[i]-model.network.average_pressure[i]));}
        for(std::size_t i=0;i<g.size();++i){require(std::abs(direct.momentum[i]-model.network.momentum[i])<1e-12,"direct momentum");require(std::abs(direct.average_flow[i]-model.network.average_flow[i])<1e-14,"direct flow");}
        const auto& s=model.source.last;work+=s.reservoir_work+s.external_work+s.geometry_work+model.last.graph_external_work;
        require(s.geometry_work==model.last.source_geometry_work,"geometry receipt");
        loss+=s.fluid_loss+s.tissue_loss+s.outlet_loss+s.numerical_loss;
        require(s.outlet_loss>=0 && std::abs(s.load_work-s.port_work-s.outlet_loss)<1e-18,"exit loss split");
        ducklips::ContinuousGlottalState state;require(model.source.state(state),"accepted source");
        worst_energy=std::max(worst_energy,std::abs(state.energy+graph_energy(model.network,c,g)-initial_energy+loss-work));
        double current=state.mass[0]+state.mass[1];for(double v:model.network.charge){current+=v;}
        transfer+=dt*(s.flow[0]+nodal[0]+nodal[1]+nodal[2]+nodal[3]-model.network.average_flow[3]);
        worst_mass=std::max(worst_mass,std::abs(current-mass-transfer));worst_port=std::max(worst_port,std::abs(model.last.port_residual));
    }
    require(worst_energy<1e-14 && worst_mass<1e-20 && worst_direct<1e-8 && worst_port<1e-20,"coupled conservation");
    const auto saved=model;ducklips::ContinuousGlottalState saved_source;require(model.source.state(saved_source),"saved source");
    require(!model.advance(c,c,g,g,2.5e-6,1400.,{},{},{},update,-2.),"reject negative outlet column");
    require(!model.advance(c,c,g,g,2.5e-6,1400.,{},{},{},update,std::numeric_limits<double>::quiet_NaN()),"reject nonfinite outlet column");
    auto bad=c;bad[0]=0.;require(!model.advance(bad,c,g,g,2.5e-6,1400.),"reject graph metric");
    std::array<double,4> force{};force[0]=std::numeric_limits<double>::quiet_NaN();
    require(!model.advance(c,c,g,g,2.5e-6,1400.,force),"reject source after graph trials");
    ducklips::ContinuousGlottalState after;require(model.source.state(after),"rollback state");
    require(physical(after)==physical(saved_source) && after.momentum==saved_source.momentum && after.energy==saved_source.energy && after.outlet_inertance==saved_source.outlet_inertance
        && model.network.charge==saved.network.charge && model.network.momentum==saved.network.momentum
        && model.network.average_pressure==saved.network.average_pressure && model.network.average_flow==saved.network.average_flow
        && model.last.pressure==saved.last.pressure && model.last.graph_port_work==saved.last.graph_port_work
        && model.source.last.energy_residual==saved.source.last.energy_residual,"rollback all accepted state");
    bool threw=false;try{prepare(model,4,edges,4);}catch(const std::invalid_argument&){threw=true;}
    require(threw && model.network.charge==saved.network.charge,"prepare transaction");
    // Exact sealed outlet: source receives no graph work, graph evolves freely.
    auto sealed=initial;sealed[2]=3e-6;sealed[3]=-1e-7;sealed[0]=sealed[1]=0.;sealed[4]=sealed[5]=sealed[6]=0.;
    prepare(model,4,edges,1,sealed);model.network.charge[1]=1e-9;
    require(model.advance(c,c,g,g,1e-7,0.,{},{},{},update),"sealed outlet advance");
    require(model.last.flow==0. && model.last.graph_port_work==0. && model.source.last.outlet_loss==0.,"sealed outlet zero transport/work");
    auto reverse=initial;reverse[0]=reverse[1]=0.;reverse[4]=reverse[5]=reverse[6]=-1e-4;
    prepare(model,4,edges,1,reverse);model.network.charge[1]=500*c[1];
    require(model.advance(c,c,g,g,1e-7,0.) && model.last.flow<0. && model.last.graph_port_work<0.
        && model.source.last.outlet_loss>0.,"reverse flow work and dissipation");
    prepare(model,4,edges,1,initial,0.);
    require(model.advance(c,c,g,g,1e-7,1400.) && model.source.last.outlet_loss==0.,"loss-free outlet coupling");
    // The hard roundoff switch formerly made the source residual discontinuous
    // at step 118 of this coupled startup. Do not loosen its force tolerance.
    prepare(model,3,edges,1);
    const std::span<const double> startup_c(c.data(),3);
    const std::array<double,3> startup_drive={0.,1e-6,-2e-6};
    model.network.charge={120*c[0],-80*c[1],40*c[2]};model.network.momentum={.01,-.02,.015,.005};
    for(int step=0;step<140;++step){require(model.advance(startup_c,startup_c,g,g,3.125e-7,1400.,{},startup_drive,series),"roundoff threshold continuity");}
    std::printf("PASS continuous graph: work %.17g J, mass %.17g m3, direct graph %.17g Pa, port %.17g J; rollback and sealed outlet\n",worst_energy,worst_mass,worst_direct,worst_port);
}
}
int main(int argc,char** argv){try{
    if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}
    else if(argc==2 && std::string(argv[1])=="--column-fixture"){fixture(true);}
    else if(argc==2 && std::string(argv[1])=="--radiation-fixture"){fixture(true,true);}
    else if(argc==2 && std::string(argv[1])=="--thermal-fixture"){fixture(true,true,true);}
    else if(argc==2 && std::string(argv[1])=="--viscous-fixture"){fixture(true,true,true,true);}
    else if(argc==2 && std::string(argv[1])=="--radiation"){radiation_controls();}
    else if(argc==2 && std::string(argv[1])=="--thermal"){radiation_controls(true);}
    else if(argc==2 && std::string(argv[1])=="--viscous"){radiation_controls(true,true);}
    else if(argc==2 && std::string(argv[1])=="--column"){controls(ducklips::ContinuousGlottalUpdate::relaxation,true);}
    else if(argc==2 && std::string(argv[1])=="--relaxation"){controls(ducklips::ContinuousGlottalUpdate::relaxation);}
    else{require(argc==1,"arguments");controls();}return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
