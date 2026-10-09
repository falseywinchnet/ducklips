#include "ducklips/continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {void require(bool b,const char* m){if(!b){throw std::runtime_error(m);}}}
int main(int argc,char** argv){try{
 require(argc==2,"mode");const std::string mode=argv[1];require(mode=="lower" || mode=="upper","mode");const bool lower=mode=="lower";ducklips::FoldParameters f;ducklips::FoldState initial;initial.rest_area=3e-6;(lower?initial.lower_position:initial.upper_position)={-2.9e-6/(2*f.length),-2.9e-6/(2*f.length)};(lower?initial.lower_velocity:initial.upper_velocity)={-.05,-.05};ducklips::ContinuousGlottis source;
 require(source.prepare(f,{},1.,initial,{}, {},1.,ducklips::ContourPressureProfile::linear_span),"prepare");ducklips::ContinuousGlottalState first,state;require(source.state(first),"initial state");double loss=0,net=0,max_work=0,max_mass=0;int closes=0,opens=0;bool was_open=true;const double dt=1e-5;
 for(int n=0;n<200;++n){require(source.advance(dt,0.,0.,0.,{},ducklips::ContinuousGlottalUpdate::dissipative),"advance");require(source.state(state),"state");const bool open=state.area[lower?0:1]>0;if(was_open && !open){++closes;}if(!was_open && open){++opens;}was_open=open;require(state.mass[0]>0 && state.mass[1]>0,"retained air");if(!open){require(source.last.flow[lower?0:2]==0 && state.flow[lower?0:2]==0,"closed exterior flow");}require(source.last.numerical_loss>=-1e-20,"negative numerical loss");loss+=source.last.fluid_loss+source.last.tissue_loss+source.last.load_work+source.last.numerical_loss;net+=dt*(source.last.flow[0]-source.last.flow[2]);max_work=std::max(max_work,std::abs(state.energy-first.energy+loss));max_mass=std::max(max_mass,std::abs(state.mass[0]+state.mass[1]-first.mass[0]-first.mass[1]-net));}
 // Opening from an exactly zero cap and zero velocity must have a defined
 // Newton starting state, without a hidden positive opening floor.
 ducklips::FoldState resting;resting.rest_area=0.;(lower?resting.upper_position:resting.lower_position)={3e-6/(2*f.length),3e-6/(2*f.length)};ducklips::ContinuousGlottis onset;
 require(onset.prepare(f,{},1.,resting,{700.,700.}, {},1.,ducklips::ContourPressureProfile::linear_span),"zero-cap prepare");require(onset.advance(1e-7,700.,0.,0.,{},ducklips::ContinuousGlottalUpdate::dissipative),"zero-cap onset");ducklips::ContinuousGlottalState opened;require(onset.state(opened) && opened.area[lower?0:1]>0,"zero-cap opening");
 require(closes==1 && opens==1,"complete contact cycle");require(max_work<1e-12 && max_mass<1e-18 && state.energy<first.energy,"cycle balances and passivity");std::printf("PASS %s contact/reopening, work %.17g J mass %.17g m3\n",mode.c_str(),max_work,max_mass);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
