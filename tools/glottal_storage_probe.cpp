#include "ducklips/glottal_storage.hpp"
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
void require(bool condition,const char* message){if(!condition){throw std::runtime_error(message);}}
void state(const ducklips::GlottalStorage& s){
    const auto& f=s.state.folds;
    std::printf("{\"mass\":[%.17g,%.17g],\"position\":[%.17g,%.17g,%.17g,%.17g],\"velocity\":[%.17g,%.17g,%.17g,%.17g],\"rest_area\":%.17g}",
        s.state.mass[0],s.state.mass[1],f.lower_position[0],f.lower_position[1],f.upper_position[0],f.upper_position[1],
        f.lower_velocity[0],f.lower_velocity[1],f.upper_velocity[0],f.upper_velocity[1],f.rest_area);
}
}
int main(int argc,char** argv){try{
    const ducklips::AirProperties air{};
    if(argc==2 && std::string(argv[1])=="--gas"){
        ducklips::GasStorage gas;require(gas.configure(air),"gas configure");double m0=0,v0=0,m1=0,v1=0;
        while(std::scanf("%lf %lf %lf %lf",&m0,&v0,&m1,&v1)==4){ducklips::GasEnergy a,b;ducklips::GasTransfer d;
            const bool ok=gas.evaluate(m0,v0,a) && gas.evaluate(m1,v1,b) && gas.transfer(m0,v0,m1,v1,d);
            if(!ok){std::puts("{\"accepted\":false}");continue;}
            std::printf("{\"accepted\":true,\"energy0\":%.17g,\"energy1\":%.17g,\"pressure0\":%.17g,\"pressure1\":%.17g,\"mass_effort\":%.17g,\"pressure\":%.17g,\"pressure_volume_slope\":%.17g}\n",a.energy,b.energy,a.pressure,b.pressure,d.mass_effort,d.pressure,d.pressure_volume_slope);
        }return 0;
    }
    const bool trajectory=argc==5 && std::string(argv[1])=="--trajectory";
    int steps=960;double duration=.005;std::string selected;
    if(trajectory){steps=std::stoi(argv[2]);duration=std::stod(argv[3]);selected=argv[4];require(steps>0 && steps<=1000000 && std::isfinite(duration) && duration>0,"trajectory parameters");}
    else{require(argc==1,"arguments");}
    for(const std::string mode:{"sealed","filling","rest_motion","asymmetric","narrow_upper","throughflow"}){
        if(trajectory && mode!=selected){continue;}
        ducklips::FoldParameters tissue;ducklips::FoldState initial;initial.rest_area=3e-6;
        if(mode=="asymmetric"){tissue.asymmetry=1.15;}
        if(mode=="narrow_upper"){initial.upper_position={-(3e-6-1e-8)/(.028),-(3e-6-1e-8)/(.028)};}
        ducklips::GlottalStorage s;const std::array<double,2> pressure=mode=="sealed"?std::array<double,2>{700.,0.}:std::array<double,2>{0.,0.};
        require(s.prepare(tissue,air,1.,initial,pressure),"prepare chambers");
        double e0=0;require(s.energy(e0),"initial energy");const auto initial_mass=s.state.mass;
        const double dt=duration/double(steps);double port=0,actuator=0,external=0,heat=0,max_residual=0,work=0,min_volume=1.;std::size_t maximum_iterations=0;
        const std::size_t before=allocations.load();
        for(int i=0;i<steps;++i){
            const std::array<double,3> flow=mode=="throughflow"?std::array<double,3>{1e-6,.6e-6,.2e-6}:(mode=="filling" || mode=="asymmetric" || mode=="narrow_upper")?std::array<double,3>{1e-6,0.,0.}:std::array<double,3>{0.,0.,0.};
            const double rest=mode=="rest_motion"?3e-6-1e-6*double(i+1)/double(steps):3e-6;
            const std::array<double,4> force=mode=="throughflow"?std::array<double,4>{1e-3,-2e-4,3e-4,-1e-4}:std::array<double,4>{};
            if(!s.advance(flow,rest,force,dt)){std::fprintf(stderr,"%s failed at %d/%d\n",mode.c_str(),i,steps);return 1;}
            std::array<double,2> volume;require(s.volumes(volume),"accepted volume");min_volume=std::min({min_volume,volume[0],volume[1]});
            port+=s.last.port_work;actuator+=s.last.actuator_work;external+=s.last.external_work;heat+=s.last.damping_loss;work+=s.last.wall_work;
            max_residual=std::max(max_residual,std::abs(s.last.energy_residual));maximum_iterations=std::max(maximum_iterations,s.last.iterations);
        }
        require(allocations.load()==before,"chamber stepping allocated");double e1=0;require(s.energy(e1),"final energy");
        const double balance=e1-e0+heat-port-actuator-external;
        require(std::abs(balance)<1e-13 && max_residual<1e-14,"complete work balance");
        const double expected=mode=="throughflow"?duration*.4e-6:mode=="filling" || mode=="asymmetric" || mode=="narrow_upper"?duration*1e-6:0;
        require(std::abs(s.state.mass[0]-initial_mass[0]-expected)<1e-20 && std::abs(s.state.mass[1]-initial_mass[1]-(mode=="throughflow"?duration*.4e-6:0.))<1e-20,"chamber mass balance");
        if(mode=="sealed"){require(e1<=e0,"sealed damped chambers gained energy");}
        if(trajectory){
            std::printf("{\"mode\":\"%s\",\"final\":",mode.c_str());state(s);
            std::printf(",\"energy0\":%.17g,\"energy1\":%.17g,\"port_work\":%.17g,\"actuator_work\":%.17g,\"wall_work\":%.17g,\"damping_loss\":%.17g,\"balance\":%.17g,\"maximum_step_balance\":%.17g,\"minimum_volume\":%.17g,\"maximum_iterations\":%zu}\n",e0,e1,port,actuator,work,heat,balance,max_residual,min_volume,maximum_iterations);
        }else{
            std::printf("PASS %s: energy balance %.17g J; wall work %.17g J; min volume %.17g m3; iterations %zu\n",mode.c_str(),balance,work,min_volume,maximum_iterations);
        }
        const auto saved=s.state;const double saved_receipt=s.last.energy_change;
        require(!s.advance({-1.,0.,0.},3e-6,{},dt),"negative remaining mass accepted");
        require(!s.advance({},std::numeric_limits<double>::quiet_NaN(),{},dt),"NaN rest opening accepted");
        require(s.state.mass==saved.mass && s.state.folds.lower_position==saved.folds.lower_position && s.state.folds.upper_position==saved.folds.upper_position
            && s.state.folds.lower_velocity==saved.folds.lower_velocity && s.state.folds.upper_velocity==saved.folds.upper_velocity
            && s.state.folds.rest_area==saved.folds.rest_area && s.last.energy_change==saved_receipt,"failed chamber step changed accepted state");
    }
    if(trajectory){require(selected=="sealed" || selected=="filling" || selected=="rest_motion" || selected=="asymmetric" || selected=="narrow_upper" || selected=="throughflow","unknown trajectory mode");}
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
