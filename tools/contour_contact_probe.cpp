#include "ducklips/contour_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
int main(int argc,char** argv){try{
    const bool fitted=argc>1 && std::string(argv[argc-1])=="--relaxation";if(fitted){--argc;}
    const bool dissipative=argc>2 && std::string(argv[argc-1])=="--dissipative";if(dissipative){--argc;}
    const auto update=fitted?ducklips::ContourUpdate::relaxation:dissipative?ducklips::ContourUpdate::dissipative:ducklips::ContourUpdate::symmetric;
    if(argc!=2){throw std::invalid_argument("usage: contour_contact_probe STEPS_OVER_100_MS");}const int steps=std::stoi(argv[1]);if(steps<1000 || steps>1000000){throw std::invalid_argument("step range");}
    const ducklips::TreeEdge edges[]={{0,1},{1,2},{0,2},{0,ducklips::ground_node}};const double c[]={2e-10,3e-10,1e-10},g[]={.002,.003,.004,.002};
    const ducklips::FoldParameters f;const ducklips::AirProperties air;ducklips::FoldState initial;initial.rest_area=3e-6;
    ducklips::GraphContourGlottal engine;engine.prepare(3,edges,2,f,air,1.,initial,{});ducklips::GasStorage gas;if(!gas.configure(air)){return 1;}
    double maximum_density=1,min_area=3e-6,max_work=0,numerical_loss=0;std::array<int,2> closures{},reopenings{};std::array<bool,2> was_open={true,true};int accepted=0;const double dt=.1/double(steps);
    for(;accepted<steps;++accepted){if(!engine.advance(c,c,g,g,dt,700.,3e-6,{},{},{},update)){break;}
        std::array<double,2> a;ducklips::ContourKineticState k;if(!engine.source.configuration(a,k)){return 1;}
        for(std::size_t i=0;i<2;++i){maximum_density=std::max(maximum_density,engine.source.state.mass[i]/k.geometry.chamber[i].volume);min_area=std::min(min_area,a[i]);}
        for(std::size_t i=0;i<2;++i){const bool open=a[i]>0;if(was_open[i] && !open){++closures[i];}if(!was_open[i] && open){++reopenings[i];}was_open[i]=open;}
        numerical_loss+=engine.source.last.numerical_loss;max_work=std::max(max_work,std::abs(engine.source.last.energy_residual));}
    const auto& source=engine.source;const auto& folds=source.state.folds;std::array<double,2> area;ducklips::ContourKineticState k;if(!source.configuration(area,k)){return 1;}
    const double rates[]={f.length*(folds.lower_velocity[0]+folds.lower_velocity[1]),f.length*(folds.upper_velocity[0]+folds.upper_velocity[1])};
    std::printf("{\"status\":\"%s\",\"planned_steps\":%d,\"accepted_steps\":%d,\"time_s\":%.17g,\"dt_s\":%.17g,\"minimum_area_m2\":%.17g,\"maximum_relative_density\":%.17g,\"maximum_step_work_residual_J\":%.17g,\"numerical_loss_J\":%.17g,\"closures\":[%d,%d],\"reopenings\":[%d,%d],\"cells\":[",accepted==steps?"completed":"solver_stopped",steps,accepted,accepted*dt,dt,min_area,maximum_density,max_work,numerical_loss,closures[0],closures[1],reopenings[0],reopenings[1]);
    for(std::size_t i=0;i<2;++i){ducklips::GasEnergy e;if(!gas.evaluate(source.state.mass[i],k.geometry.chamber[i].volume,e)){return 1;}
        std::printf("%s{\"area_m2\":%.17g,\"area_rate_m2_s\":%.17g,\"volume_m3\":%.17g,\"mass_volume_m3\":%.17g,\"relative_density\":%.17g,\"pressure_Pa\":%.17g}",i?",":"",area[i],rates[i],k.geometry.chamber[i].volume,source.state.mass[i],source.state.mass[i]/k.geometry.chamber[i].volume,e.pressure);}
    std::printf("],\"momentum_Pa_s\":[%.17g,%.17g,%.17g],\"average_flow_m3_s\":[%.17g,%.17g,%.17g]}\n",source.momentum[0],source.momentum[1],source.momentum[2],source.last.flow[0],source.last.flow[1],source.last.flow[2]);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
