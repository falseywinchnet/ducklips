#include "ducklips/contour_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void snapshot(const ducklips::GraphContourGlottal& e,double time,double dt,const char* tag){
    std::array<double,2> a;ducklips::ContourKineticState k;
    if(!e.source.configuration(a,k)){throw std::runtime_error("configuration");}
    const auto& f=e.source.state.folds;
    std::printf("{\"tag\":\"%s\",\"time\":%.17g,\"dt\":%.17g,\"area\":[%.17g,%.17g],\"volume\":[%.17g,%.17g],\"state\":[",tag,time,dt,a[0],a[1],k.geometry.chamber[0].volume,k.geometry.chamber[1].volume);
    bool first=true;const auto number=[&](double x){std::printf("%s%.17g",first?"":",",x);first=false;};
    for(double x:e.source.state.mass){number(x);}for(double x:f.lower_position){number(x);}for(double x:f.upper_position){number(x);}for(double x:f.lower_velocity){number(x);}for(double x:f.upper_velocity){number(x);}for(double x:e.source.momentum){number(x);}
    const double c[]={2e-10,3e-10,1e-10},g[]={.002,.003,.004,.002};
    for(std::size_t i=0;i<3;++i){number(e.network.charge[i]/c[i]);}for(std::size_t i=0;i<4;++i){number(e.network.momentum[i]*g[i]);}
    std::printf("]}\n");std::fflush(stdout);
}
}
int main(int argc,char** argv){try{
    if(argc!=2){throw std::invalid_argument("usage: contour_pocket_probe STEPS_OVER_100_MS");}
    const int steps=std::stoi(argv[1]);if(steps<1000 || steps>1000000){throw std::invalid_argument("steps");}
    const ducklips::TreeEdge edges[]={{0,1},{1,2},{0,2},{0,ducklips::ground_node}};
    const double c[]={2e-10,3e-10,1e-10},g[]={.002,.003,.004,.002};
    ducklips::FoldState initial;initial.rest_area=3e-6;ducklips::GraphContourGlottal e;e.prepare(3,edges,2,{},{},1.,initial,{});
    const double base=.1/steps;double dt=base,time=0,maximum_work=0,loss=0;std::size_t accepted=0,rejected=0;bool checkpoint=false;double next_volume=1e-14,minimum_volume=1.,maximum_density=1.;
    const char* status="step_budget";
    for(std::size_t n=0;n<200000;++n){
        if(!checkpoint && time>=.01088){snapshot(e,time,dt,"checkpoint");checkpoint=true;}
        if(time>=.011){status="completed";break;}
        dt=std::min(dt,.011-time);
        if(!e.advance(c,c,g,g,dt,700.,3e-6,{},{},{},ducklips::ContourUpdate::relaxation)){
            ++rejected;if(rejected==1){snapshot(e,time,dt,"first_rejection");}
            dt*=.5;if(dt<1e-14){status="minimum_step";break;}continue;
        }
        time+=dt;++accepted;maximum_work=std::max(maximum_work,std::abs(e.source.last.energy_residual));loss+=e.source.last.numerical_loss;
        std::array<double,2> a;ducklips::ContourKineticState k;if(!e.source.configuration(a,k)){throw std::runtime_error("state");}
        minimum_volume=std::min(minimum_volume,k.geometry.chamber[0].volume);maximum_density=std::max(maximum_density,e.source.state.mass[0]/k.geometry.chamber[0].volume);
        if(k.geometry.chamber[0].volume<next_volume){snapshot(e,time,dt,"volume_decade");while(next_volume>k.geometry.chamber[0].volume){next_volume*=.1;}}
        dt=std::min(base,dt*1.2);
    }
    snapshot(e,time,dt,"final");std::printf("{\"status\":\"%s\",\"accepted\":%zu,\"rejected\":%zu,\"maximum_work_residual_J\":%.17g,\"numerical_loss_J\":%.17g,\"minimum_volume_m3\":%.17g,\"maximum_relative_density\":%.17g}\n",status,accepted,rejected,maximum_work,loss,minimum_volume,maximum_density);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
