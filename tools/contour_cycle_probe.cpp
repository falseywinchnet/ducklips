#include "ducklips/contour_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {void require(bool value,const char* text){if(!value){throw std::runtime_error(text);}}}
int main(int argc,char** argv){try{
    const bool fitted=argc>1 && std::string(argv[argc-1])=="--relaxation";if(fitted){--argc;}
    const auto update=fitted?ducklips::ContourUpdate::relaxation:ducklips::ContourUpdate::dissipative;
    const int steps=argc>1?std::stoi(argv[1]):2000;const std::string mode=argc>2?argv[2]:"lower";require(argc<=3 && steps>0 && steps<=1000000 && (mode=="lower" || mode=="upper"),"arguments");
    const bool lower=mode=="lower";ducklips::FoldParameters f;ducklips::FoldState initial;initial.rest_area=3e-6;
    (lower?initial.lower_position:initial.upper_position)={-2.9e-6/(2*f.length),-2.9e-6/(2*f.length)};
    (lower?initial.lower_velocity:initial.upper_velocity)={-.05,-.05};
    ducklips::ContourGlottalStorage source;require(source.prepare(f,{},1.,initial,{}),"cycle prepare");
    double e0=0;require(source.energy(e0),"initial energy");const double m0=source.state.mass[0]+source.state.mass[1],dt=.002/double(steps);
    double numerical=0,physical=0,net=0,max_balance=0,max_mass=0,min_area=1.,max_density=0,first_close=-1,first_open=-1,max_sealing_energy=0;
    std::array<double,2> previous_area;ducklips::ContourKineticState previous;require(source.configuration(previous_area,previous),"initial configuration");
    int closes=0,opens=0;bool was_open=true;
    for(int i=0;i<steps;++i){const auto old_momentum=source.momentum;if(!source.advance(dt,0,0,0,3e-6,{},update)){std::fprintf(stderr,"cycle %s failed at %d/%d\n",mode.c_str(),i,steps);return 1;}
        std::array<double,2> a;ducklips::ContourKineticState k;require(source.configuration(a,k),"configuration");const bool open=a[lower?0:1]>0;
        if(was_open && !open){const std::size_t j=lower?0:2;max_sealing_energy=std::max(max_sealing_energy,.5*previous.inverse_inertance[j]*old_momentum[j]*old_momentum[j]);++closes;if(first_close<0){first_close=double(i+1)*dt;}}if(!was_open && open){++opens;if(first_open<0){first_open=double(i+1)*dt;}}was_open=open;previous=k;
        for(std::size_t cell=0;cell<2;++cell){min_area=std::min(min_area,a[cell]);max_density=std::max(max_density,source.state.mass[cell]/k.geometry.chamber[cell].volume);
            if(a[cell]<=0){const std::size_t j=cell==0?0:2;require(source.momentum[j]==0 && source.last.flow[j]==0,"closed flux/momentum");}}
        numerical+=source.last.numerical_loss;physical+=source.last.passage_loss+source.last.damping_loss;net+=dt*(source.last.flow[0]-source.last.flow[2]);
        require(source.last.numerical_loss>=0,"negative numerical loss");double e=0;require(source.energy(e),"energy");max_balance=std::max(max_balance,std::abs(e-e0+numerical+physical));max_mass=std::max(max_mass,std::abs(source.state.mass[0]+source.state.mass[1]-m0-net));
    }
    require(closes>0 && opens>0,"no complete contact/reopening cycle");require(max_balance<1e-12 && max_mass<1e-18,"cycle balances");double e1=0;require(source.energy(e1) && e1<=e0,"passivity");
    const auto& s=source.state.folds;std::printf("{\"mode\":\"%s\",\"steps\":%d,\"duration_s\":0.002,\"closures\":%d,\"reopenings\":%d,\"first_closed_endpoint_s\":%.17g,\"first_reopened_endpoint_s\":%.17g,\"initial_energy_J\":%.17g,\"final_energy_J\":%.17g,\"maximum_energy_before_sealing_J\":%.17g,\"numerical_loss_J\":%.17g,\"physical_loss_J\":%.17g,\"maximum_work_residual_J\":%.17g,\"maximum_mass_residual_m3\":%.17g,\"minimum_area_m2\":%.17g,\"maximum_density_ratio\":%.17g,\"mass\":[%.17g,%.17g],\"position\":[%.17g,%.17g,%.17g,%.17g],\"velocity\":[%.17g,%.17g,%.17g,%.17g],\"momentum\":[%.17g,%.17g,%.17g]}\n",mode.c_str(),steps,closes,opens,first_close,first_open,e0,e1,max_sealing_energy,numerical,physical,max_balance,max_mass,min_area,max_density,source.state.mass[0],source.state.mass[1],s.lower_position[0],s.lower_position[1],s.upper_position[0],s.upper_position[1],s.lower_velocity[0],s.lower_velocity[1],s.upper_velocity[0],s.upper_velocity[1],source.momentum[0],source.momentum[1],source.momentum[2]);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
