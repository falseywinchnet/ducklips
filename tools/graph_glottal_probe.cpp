#include "ducklips/graph_glottal_storage.hpp"
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
void* operator new(std::size_t n){++allocations;if(void* p=std::malloc(n)){return p;}throw std::bad_alloc();}
void operator delete(void* p) noexcept{std::free(p);}void operator delete(void* p,std::size_t) noexcept{std::free(p);}
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
void values(const char* name,std::span<const double> v){std::printf("\"%s\":[",name);for(std::size_t i=0;i<v.size();++i){std::printf("%s%.17g",i?",":"",v[i]);}std::putchar(']');}
double energy(const ducklips::GraphGlottalStorage& p,std::span<const double> c,std::span<const double> g){
    double e=0;require(p.source_energy(e),"source energy");
    for(std::size_t i=0;i<c.size();++i){e+=.5*p.network.charge[i]*p.network.charge[i]/c[i];}
    for(std::size_t i=0;i<g.size();++i){e+=.5*g[i]*p.network.momentum[i]*p.network.momentum[i];}return e;
}
}
int main(int argc,char** argv){try{
    const bool articulated=argc>1 && std::string(argv[argc-1])=="--articulated";
    if(articulated){--argc;}
    const bool trajectory=argc==5 && std::string(argv[1])=="--trajectory";
    int steps=960;double duration=.005;std::string selected;
    if(trajectory){steps=std::stoi(argv[2]);duration=std::stod(argv[3]);selected=argv[4];require(steps>0 && steps<=3000000 && std::isfinite(duration) && duration>0,"invalid trajectory");}
    else{require(argc==1,"unknown arguments");}
    const ducklips::TreeEdge edges[]={{0,1},{1,2},{0,2},{0,ducklips::ground_node}};
    const double c[]={2e-10,3e-10,1e-10},g[]={.002,.003,.004,.002};
    const ducklips::AirProperties air{};const ducklips::FoldParameters tissue{};
    const ducklips::GlottalPassages links{{150.,200.,100.},{1e6,2e6,1e6},{6e10,4e10,6e10}};
    bool found=false;
    for(const std::string mode:{"forward","reverse","relaxation","rest_motion"}){
        if(trajectory && mode!=selected){continue;}found=true;
        ducklips::GraphGlottalStorage p;ducklips::FoldState initial;initial.rest_area=3e-6;
        const auto initial_pressure=mode=="relaxation"?std::array<double,2>{700.,300.}:std::array<double,2>{};
        if(articulated){p.prepare_articulated(3,edges,2,tissue,air,1.,initial,initial_pressure);}
        else{p.prepare(3,edges,2,tissue,air,1.,initial,initial_pressure,links);}
        if(mode=="relaxation"){p.momentum={.001,-.0005,.0004};p.network.charge={2e-8,-1e-8,3e-8};}
        const double initial_mass=p.chambers.state.mass[0]+p.chambers.state.mass[1]+p.network.charge[0]+p.network.charge[1]+p.network.charge[2];
        const double e0=energy(p,c,g),dt=duration/double(steps);double work=0,loss=0,actuator=0,entered=0,left=0,max_mass=0,max_law=0,max_balance=0,max_local=0,max_flow_difference=0,min_in=0,max_in=0;std::size_t max_iterations=0;
        const std::size_t before=allocations.load();
        for(int i=0;i<steps;++i){
            const double reservoir=mode=="relaxation"?0.:mode=="reverse"?-200.:700.;
            const double rest=mode=="rest_motion"?3e-6-5e-7*double(i+1)/double(steps):3e-6;
            if(!p.advance(c,c,g,g,dt,reservoir,rest)){std::fprintf(stderr,"%s failed at %d/%d\n",mode.c_str(),i,steps);return 1;}
            entered+=dt*p.last.flow[0];left+=dt*p.network.average_flow.back();work+=p.last.reservoir_work;
            loss+=p.last.passage_loss+p.chambers.last.damping_loss;actuator+=p.chambers.last.actuator_work;
            const double stored=p.chambers.state.mass[0]+p.chambers.state.mass[1]+p.network.charge[0]+p.network.charge[1]+p.network.charge[2];
            max_mass=std::max(max_mass,std::abs(stored-initial_mass+left-entered));
            max_balance=std::max(max_balance,std::abs(energy(p,c,g)-e0+loss-work-actuator));max_local=std::max(max_local,std::abs(p.last.energy_residual));
            for(double r:p.last.law_residual){max_law=std::max(max_law,std::abs(r));}
            max_iterations=std::max(max_iterations,p.last.iterations);max_flow_difference=std::max(max_flow_difference,std::abs(p.last.flow[0]-p.last.flow[2]));
            min_in=std::min(min_in,p.last.flow[0]);max_in=std::max(max_in,p.last.flow[0]);
        }
        require(allocations.load()==before,"coupled stepping allocated");require(max_mass<1e-17,"coupled mass balance");
        require(max_balance<1e-11 && max_local<1e-12,"coupled work balance");if(max_law>=1e-5){std::array<double,2> volume{};require(p.chambers.volumes(volume),"diagnostic volume");
            std::fprintf(stderr,"law %.17g Pa, work %.17g J, mass %.17g m3; volumes %.17g %.17g; masses %.17g %.17g; momenta %.17g %.17g %.17g\n",max_law,max_balance,max_mass,volume[0],volume[1],p.chambers.state.mass[0],p.chambers.state.mass[1],p.momentum[0],p.momentum[1],p.momentum[2]);}
        require(max_law<1e-5,"passage momentum balance");
        require(max_flow_difference>1e-9,"storage forced one through-flow");
        if(mode=="forward"){require(max_in>0,"no forward flow");}if(mode=="reverse"){require(min_in<0,"no reverse flow");}
        if(mode=="relaxation"){require(energy(p,c,g)<=e0,"undriven system gained energy");}
        if(trajectory){
            std::putchar('{');values("charge",p.network.charge);std::putchar(',');values("network_momentum",p.network.momentum);std::putchar(',');values("passage_momentum",p.momentum);
            std::putchar(',');values("mass",p.chambers.state.mass);const auto& f=p.chambers.state.folds;
            const double x[]={f.lower_position[0],f.lower_position[1],f.upper_position[0],f.upper_position[1]},v[]={f.lower_velocity[0],f.lower_velocity[1],f.upper_velocity[0],f.upper_velocity[1]};
            std::putchar(',');values("position",x);std::putchar(',');values("velocity",v);
            std::printf(",\"maximum_mass_balance\":%.17g,\"maximum_work_balance\":%.17g,\"maximum_local_balance\":%.17g,\"maximum_law_residual\":%.17g,\"maximum_flow_difference\":%.17g,\"maximum_iterations\":%zu}\n",max_mass,max_balance,max_local,max_law,max_flow_difference,max_iterations);
        }else{std::printf("PASS %s mass %.3g m3 work %.3g J law %.3g Pa; iterations %zu\n",mode.c_str(),max_mass,max_balance,max_law,max_iterations);}
        const auto saved=p;
        require(!p.advance(c,c,g,g,dt,std::numeric_limits<double>::quiet_NaN(),3e-6),"NaN pressure accepted");
        const double bad[]={2e-10,0,1e-10};require(!p.advance(bad,c,g,g,dt,700,3e-6),"zero compliance accepted");
        require(!p.advance(c,c,g,g,dt,700,std::numeric_limits<double>::quiet_NaN()),"NaN rest accepted");
        require(p.network.charge==saved.network.charge && p.network.momentum==saved.network.momentum && p.momentum==saved.momentum
            && p.chambers.state.mass==saved.chambers.state.mass && p.chambers.state.folds.lower_position==saved.chambers.state.folds.lower_position
            && p.chambers.state.folds.upper_position==saved.chambers.state.folds.upper_position && p.chambers.state.folds.lower_velocity==saved.chambers.state.folds.lower_velocity
            && p.chambers.state.folds.upper_velocity==saved.chambers.state.folds.upper_velocity && p.last.flow==saved.last.flow
            && p.last.energy_residual==saved.last.energy_residual,"rejected step changed state");
        ducklips::GraphGlottalStorage restored;
        if(articulated){restored.prepare_articulated(3,edges,2,tissue,air,1.,initial,{});}
        else{restored.prepare(3,edges,2,tissue,air,1.,initial,{},links);}
        restored.network.charge=p.network.charge;restored.network.momentum=p.network.momentum;restored.momentum=p.momentum;restored.chambers.state=p.chambers.state;
        require(p.advance(c,c,g,g,dt,700,3e-6) && restored.advance(c,c,g,g,dt,700,3e-6),"continuation failed");
        require(p.network.charge==restored.network.charge && p.network.momentum==restored.network.momentum && p.momentum==restored.momentum
            && p.chambers.state.mass==restored.chambers.state.mass && p.chambers.state.folds.lower_position==restored.chambers.state.folds.lower_position
            && p.chambers.state.folds.upper_position==restored.chambers.state.folds.upper_position,"restored physical state diverged");
    }
    require(found,"unknown mode");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
