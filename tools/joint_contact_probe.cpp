#include "ducklips/joint_accuracy.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <string>

namespace {std::atomic<std::size_t> allocations{};}
void* operator new(std::size_t size){++allocations;if(void* p=std::malloc(size)){return p;}throw std::bad_alloc();}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}

namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
double number(){double value;require(std::scanf("%lf",&value)==1 && std::isfinite(value),"restart number");return value;}
bool geometry(void* context,double,double,ducklips::JointAdvanceInputs& input) noexcept {
    input=*static_cast<ducklips::JointAdvanceInputs*>(context);return true;
}
double energy(const ducklips::GraphContinuousGlottis& m,const ducklips::ContinuousGlottalState& s,
              std::span<const double> c,std::span<const double> g){
    double result=s.energy+m.radiation.energy();
    for(std::size_t i=0;i<c.size();++i){result+=.5*m.network.charge[i]*m.network.charge[i]/c[i];}
    for(std::size_t i=0;i<g.size();++i){result+=.5*g[i]*m.network.momentum[i]*m.network.momentum[i];}
    for(const auto& t:m.thermal){result+=t.energy();}for(const auto& v:m.viscous){result+=v.energy();}return result;
}
double mass(const ducklips::GraphContinuousGlottis& m,const ducklips::ContinuousGlottalState& s){
    double result=s.mass[0]+s.mass[1];for(double q:m.network.charge){result+=q;}
    for(const auto& t:m.thermal){result+=t.stored_air();}return result;
}
void values(const char* key,std::span<const double> values){
    std::printf(",\"%s\":[",key);for(std::size_t i=0;i<values.size();++i){std::printf("%s%.17g",i?",":"",values[i]);}std::printf("]");
}
}
int main(int argc,char** argv){try{
    require(argc>=6,"mode duration macro_step tolerance contact_step [--initial-state] [--stop-after count] [--compare-step seconds]");
    bool restart=false;std::size_t stop_after=0;double compare_step=0;
    for(int i=6;i<argc;++i){const std::string option=argv[i];
        if(option=="--initial-state" && !restart){restart=true;}
        else if(option=="--compare-step" && compare_step==0 && i+1<argc){
            std::size_t used=0;const std::string value=argv[++i];compare_step=std::stod(value,&used);
            require(used==value.size() && std::isfinite(compare_step) && compare_step>0,"positive comparison step");
        }
        else if(option=="--stop-after" && !stop_after && i+1<argc){
            std::size_t used=0;const std::string value=argv[++i];stop_after=std::stoul(value,&used);
            require(used==value.size() && stop_after>0,"positive checkpoint index");
        }else{require(false,"unknown option");}
    }
    const std::string mode=argv[1];const double duration=std::stod(argv[2]),macro_step=std::stod(argv[3]),tolerance=std::stod(argv[4]),contact_step=std::stod(argv[5]);
    require(mode=="lower" || mode=="upper" || mode=="zero-lower" || mode=="zero-upper" || mode=="phonation","mode");
    for(double value:{duration,macro_step,tolerance,contact_step}){require(std::isfinite(value) && value>0,"positive finite settings");}
    const auto count=static_cast<std::size_t>(std::ceil(duration/macro_step));require(count>0 && count<=1000000,"macro interval count");
    require(stop_after<=count,"checkpoint outside trajectory");const std::size_t last_step=stop_after?stop_after:count;
    const bool zero=mode.starts_with("zero-"),lower=mode=="lower" || mode=="zero-lower",phonation=mode=="phonation";
    ducklips::AirProperties air;ducklips::FoldParameters f;ducklips::FoldState folds;folds.rest_area=3e-6;
    double pressure=0,reservoir=0;
    if(phonation){folds.lower_position={1e-6,1e-6};reservoir=1400.;}
    else if(zero){folds.rest_area=0;(lower?folds.upper_position:folds.lower_position).fill(3e-6/(2*f.length));pressure=reservoir=700.;}
    else{(lower?folds.lower_position:folds.upper_position).fill(-2.9e-6/(2*f.length));(lower?folds.lower_velocity:folds.upper_velocity).fill(-.05);}
    std::array<double,2> effort={pressure,pressure};std::array<double,3> initial_flow{};
    if(restart){
        folds.rest_area=number();for(double& v:effort){v=number();}
        for(auto* values:{&folds.lower_position,&folds.upper_position,&folds.lower_velocity,&folds.upper_velocity}){for(double& v:*values){v=number();}}
        for(double& v:initial_flow){v=number();}
    }
    constexpr std::size_t nodes=8;constexpr double area=2e-4,dx=.175/nodes,gap=area/(2*.03);
    std::vector<ducklips::TreeEdge> edges;for(std::size_t i=0;i<nodes;++i){edges.push_back({i,i+1<nodes?i+1:ducklips::ground_node});}
    std::vector<double> c(nodes,area*dx/(air.density*air.sound_speed*air.sound_speed)),g(nodes,area/(air.density*dx));g.back()*=2;
    ducklips::GraphContinuousGlottis model;
    model.prepare(nodes,edges,0,f,air,1.,folds,effort,initial_flow,1.,ducklips::ContourPressureProfile::linear_span,air.density*dx/(2*area));
    std::vector<ducklips::ThermalSite> thermal_sites;std::vector<ducklips::ThermalGeometry> thermal;
    std::vector<ducklips::ViscousSite> viscous_sites;std::vector<ducklips::ViscousGeometry> viscous;
    for(std::size_t i=0;i<nodes;++i){
        thermal_sites.push_back({i,c[i],gap});thermal.push_back({c[i],gap});
        viscous_sites.push_back({i,{1/g[i],g[i],gap}});viscous.push_back({1/g[i],g[i],gap});model.network.charge[i]=c[i]*pressure;
    }
    model.prepare_thermal(air,thermal_sites);model.prepare_viscous(air,viscous_sites);
    for(auto& site:model.thermal){for(std::size_t i=0;i<ducklips::SlitThermal::modes;++i){site.charge[i]=site.memory_storage(i)*pressure;}}
    require(model.radiation.prepare(air,area),"radiation");
    if(restart){
        for(double& v:model.network.charge){v=number();}for(double& v:model.network.momentum){v=number();}model.radiation.momentum=number();
        for(auto& site:model.thermal){for(double& v:site.charge){v=number();}}for(auto& site:model.viscous){for(double& v:site.momentum){v=number();}}
    }
    model.prepare_joint(nodes-1,ducklips::RelaxationOrder::fifth);
    ducklips::JointAdvanceInputs input;input.geometry={c,c,c,g,g,g,thermal,thermal,viscous,viscous,area,area};input.reservoir_pressure=reservoir;
    if(compare_step>0){
        auto full=model,halves=model;
        require(full.advance_joint(input.geometry,compare_step,reservoir),"comparison full step");
        require(halves.advance_joint(input.geometry,.5*compare_step,reservoir)
            && halves.advance_joint(input.geometry,.5*compare_step,reservoir),"comparison half steps");
        std::printf("{\"step_s\":%.17g,\"states\":[",compare_step);bool comma=false;
        for(const auto* m:{&full,&halves}){
            ducklips::ContinuousGlottalState s;require(m->source.state(s),"comparison source");
            std::printf("%s{\"energy\":%.17g",comma?",":"",s.energy);comma=true;
            values("area",s.area);values("effort",s.effort);values("mass",s.mass);values("flow",s.flow);
            values("lower_position",s.folds.lower_position);values("upper_position",s.folds.upper_position);
            values("lower_velocity",s.folds.lower_velocity);values("upper_velocity",s.folds.upper_velocity);
            values("charge",m->network.charge);values("momentum",m->network.momentum);std::putchar('}');
        }
        std::printf("]}\n");return 0;
    }
    ducklips::JointAccuracy controller;controller.prepare(model);
    ducklips::JointAccuracyOptions options{tolerance,macro_step,1e-13,4096,contact_step};
    ducklips::ContinuousGlottalState first,state;require(model.source.state(first),"initial source");
    const double e0=energy(model,first,c,g),m0=mass(model,first);double work=0,loss=0,numerical=0,net=0,worst_work=0,worst_mass=0;
    std::size_t closures=0,reopenings=0,substeps=0,rejected=0;
    for(std::size_t i=0;i<=last_step;++i){
        const double time=duration*static_cast<double>(i)/static_cast<double>(count);
        if(i){
            const double start=duration*static_cast<double>(i-1)/static_cast<double>(count);const auto allocation_count=allocations.load();
            if(!controller.advance(model,start,time,geometry,&input,options)){
                std::fprintf(stderr,"contact advance failed mode=%s macro=%zu reason=%d time=%.17g trial=%.17g error=%.17g\n",mode.c_str(),i,static_cast<int>(controller.failure),controller.attempt_time,controller.attempt_step,controller.attempt_error);return 1;
            }
            require(allocations.load()==allocation_count,"contact advance allocated");
            const auto& r=controller.last;work+=r.external_work;loss+=r.physical_loss;numerical+=r.numerical_loss;net+=r.reservoir_volume-(time-start)*controller.average_flow.back();
            closures+=r.source_closures;reopenings+=r.source_reopenings;substeps+=r.accepted_steps;rejected+=r.rejected;
            require(r.physical_loss>=0 && r.numerical_loss>=0 && r.maximum_accepted_error<=tolerance,"accepted accuracy and nonnegative losses");
        }
        require(model.source.state(state),"source state");
        require(state.mass[0]>0 && state.mass[1]>0,"positive source gas");
        for(std::size_t j=0;j<2;++j){if(state.area[j]<=0){require(state.flow[j==0?0:2]==0,"zero closed-port flow");}}
        worst_work=std::max(worst_work,std::abs(energy(model,state,c,g)-e0+loss+numerical-work));worst_mass=std::max(worst_mass,std::abs(mass(model,state)-m0-net));
        std::printf("{\"time\":%.17g,\"closures\":%zu,\"reopenings\":%zu,\"substeps\":%zu,\"rejected\":%zu,\"energy\":%.17g,\"work\":%.17g,\"physical_loss\":%.17g,\"numerical_loss\":%.17g,\"work_residual\":%.17g,\"mass_residual\":%.17g",time,closures,reopenings,substeps,rejected,energy(model,state,c,g),work,loss,numerical,worst_work,worst_mass);
        std::printf(",\"outlet_flow\":%.17g,\"radiation_flow\":%.17g,\"radiation_resistance\":%.17g",g.back()*model.network.momentum.back(),model.radiation.inverse_inertance()*model.radiation.momentum,model.radiation.resistance());
        values("area",state.area);values("gas_mass",state.mass);values("effort",state.effort);values("flow",state.flow);values("lower_position",state.folds.lower_position);values("upper_position",state.folds.upper_position);
        std::array<double,nodes> pressure_values;for(std::size_t j=0;j<nodes;++j){pressure_values[j]=model.network.charge[j]/c[j];}values("pressure",pressure_values);
        if(i==last_step || (restart && i==0)){
            // Physical-state continuation for this fixed tube/body control.
            // The numeric stdin order is source, charge, momentum, radiation,
            // thermal memories, viscous memories. It does not reset profiles.
            std::printf(",\"restart\":{\"version\":1");
            const std::array<double,14> source={state.folds.rest_area,state.effort[0],state.effort[1],state.folds.lower_position[0],state.folds.lower_position[1],state.folds.upper_position[0],state.folds.upper_position[1],state.folds.lower_velocity[0],state.folds.lower_velocity[1],state.folds.upper_velocity[0],state.folds.upper_velocity[1],state.flow[0],state.flow[1],state.flow[2]};
            values("source",source);values("charge",model.network.charge);values("momentum",model.network.momentum);std::printf(",\"radiation\":%.17g",model.radiation.momentum);
            std::vector<double> thermal_memory,viscous_memory;for(const auto& site:model.thermal){thermal_memory.insert(thermal_memory.end(),site.charge.begin(),site.charge.end());}for(const auto& site:model.viscous){viscous_memory.insert(viscous_memory.end(),site.momentum.begin(),site.momentum.end());}
            values("thermal",thermal_memory);values("viscous",viscous_memory);std::printf("}");
        }
        std::printf("}\n");std::fflush(stdout);
    }
    require(worst_work<1e-12 && worst_mass<1e-18,"whole coupled work and mass");
    if(!phonation && !restart){require(reopenings>0 && (zero || closures>0),"contact/reopening occurred");}
    std::fprintf(stderr,"PASS %s: closures=%zu reopenings=%zu accepted=%zu rejected=%zu work=%.17g mass=%.17g; no speech\n",mode.c_str(),closures,reopenings,substeps,rejected,worst_work,worst_mass);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
