#include "ducklips/continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* text){if(!value){throw std::runtime_error(text);}}
template<std::size_t N> void array(const std::array<double,N>& x){std::printf("[");for(std::size_t j=0;j<N;++j){std::printf("%s%.17g",j?",":"",x[j]);}std::printf("]");}
}
int main(int argc,char** argv){try{
    if(argc==1){
        ducklips::ContinuousGlottis source;ducklips::FoldState folds;folds.rest_area=3e-6;require(source.prepare({},{},1.,folds,{}),"equilibrium prepare");ducklips::ContinuousGlottalState initial,after;require(source.state(initial),"equilibrium state");require(source.advance(1e-6,0.,0.,0.) && source.state(after),"equilibrium advance");require(after.energy==initial.energy && after.flow==initial.flow && after.mass==initial.mass,"equilibrium moved");require(!source.advance(-1.,0.,0.,0.) && source.state(after) && after.mass==initial.mass,"transactional rejection");std::printf("PASS continuous coupled equilibrium and invalid-step rejection\n");return 0;
    }
    const std::string mode=argv[1];const int steps=argc>2?std::stoi(argv[2]):100;const double duration=argc>3?std::stod(argv[3]):.0001;const double amplitude=argc>4?std::stod(argv[4]):1.;require(argc<=5 && steps>0 && duration>0 && amplitude>0,"arguments");require(mode=="free" || mode=="linear" || mode=="driven" || mode=="lower" || mode=="upper","mode");
    ducklips::FoldParameters f;ducklips::AirProperties air;ducklips::FoldState folds;folds.rest_area=3e-6;std::array<double,2> effort={100.,-30.};std::array<double,3> flow={1e-6,-2e-6,3e-6};double inlet=0.,outlet_loss=0.;
    folds.lower_position={1e-6,-2e-6};folds.upper_position={-1e-6,3e-6};folds.lower_velocity={.001,-.002};folds.upper_velocity={.003,-.001};f.asymmetry=1.17;
    if(mode=="free" || mode=="linear"){air.viscosity=0.;f.lower_damping_ratio=f.upper_damping_ratio=0.;}
    if(mode=="driven"){inlet=700.;outlet_loss=1.;effort={0.,0.};flow={};}
    if(mode=="lower"){folds.lower_position={-4e-6/(2*f.length),-4e-6/(2*f.length)};flow[0]=0;}
    if(mode=="upper"){folds.upper_position={-4e-6/(2*f.length),-4e-6/(2*f.length)};flow[2]=0;}
    if(mode=="linear"){for(double& v:effort){v*=amplitude;}for(double& v:flow){v*=amplitude;}for(auto* a:{&folds.lower_position,&folds.upper_position,&folds.lower_velocity,&folds.upper_velocity}){for(double& v:*a){v*=amplitude;}}}
    ducklips::ContinuousGlottis source;require(source.prepare(f,air,.9,folds,effort,flow,outlet_loss),"motion prepare");ducklips::ContinuousGlottalState initial,state;require(source.state(initial),"initial state");double loss=0,work=0,net=0,max_energy=0,max_mass=0,max_residual=0;std::size_t max_iterations=0;const double dt=duration/steps;
    for(int n=0;n<steps;++n){if(!source.advance(dt,inlet,0.,0.)){std::fprintf(stderr,"%s failed step %d/%d dt %.17g\n",mode.c_str(),n,steps,dt);return 1;}require(source.state(state),"state");const auto& s=source.last;loss+=s.fluid_loss+s.tissue_loss;work+=s.reservoir_work-s.load_work+s.external_work;net+=dt*(s.flow[0]-s.flow[2]);max_energy=std::max(max_energy,std::abs(state.energy-initial.energy+loss-work));max_mass=std::max(max_mass,std::abs(state.mass[0]+state.mass[1]-initial.mass[0]-initial.mass[1]-net));max_residual=std::max(max_residual,s.maximum_scaled_residual);max_iterations=std::max(max_iterations,s.iterations);require(s.fluid_loss>=-1e-20 && s.tissue_loss>=0,"nonnegative loss");}
    require(max_energy<1e-12 && max_mass<1e-20,"motion balance");std::printf("{\"mode\":\"%s\",\"steps\":%d,\"duration_s\":%.17g,\"energy0_J\":%.17g,\"energy1_J\":%.17g,\"physical_loss_J\":%.17g,\"net_work_J\":%.17g,\"maximum_energy_residual_J\":%.17g,\"maximum_mass_residual_m3\":%.17g,\"maximum_scaled_residual\":%.17g,\"maximum_iterations\":%zu,\"mass\":",mode.c_str(),steps,duration,initial.energy,state.energy,loss,work,max_energy,max_mass,max_residual,max_iterations);array(state.mass);std::printf(",\"effort\":");array(state.effort);std::printf(",\"area\":");array(state.area);std::printf(",\"flow\":");array(state.flow);std::printf(",\"position\":");array(std::array<double,4>{state.folds.lower_position[0],state.folds.lower_position[1],state.folds.upper_position[0],state.folds.upper_position[1]});std::printf(",\"velocity\":");array(std::array<double,4>{state.folds.lower_velocity[0],state.folds.lower_velocity[1],state.folds.upper_velocity[0],state.folds.upper_velocity[1]});std::printf(",\"momentum\":");array(state.momentum);std::printf("}\n");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
