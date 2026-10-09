#include "ducklips/continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
int main(){try{
 ducklips::FoldParameters f;ducklips::FoldState initial;initial.rest_area=3e-6;initial.lower_position.fill((3.6868231667933133e-6-initial.rest_area)/(2*f.length)+1e-9);initial.upper_position.fill((4.149026342131548e-6-initial.rest_area)/(2*f.length));ducklips::ContinuousGlottis continuous,restarted;
 if(!continuous.prepare(f,{},1.,initial,{112.18465038080178,-3.138675468892969},{.00017466261521028343,.00017466261521028343,.00017466261521028343},1.,ducklips::ContourPressureProfile::linear_span)){throw std::runtime_error("prepare");}restarted=continuous;double worst=0,work=0;
 for(int n=0;n<1000;++n){ducklips::ContinuousGlottalState before,a,b;if(!restarted.state(before) || !restarted.prepare(f,{},1.,before.folds,before.effort,before.flow,1.,ducklips::ContourPressureProfile::linear_span)){throw std::runtime_error("restart");}
  if(!continuous.advance(1e-5,1400.,0.,0.,{},ducklips::ContinuousGlottalUpdate::dissipative) || !restarted.advance(1e-5,1400.,0.,0.,{},ducklips::ContinuousGlottalUpdate::dissipative) || !continuous.state(a) || !restarted.state(b)){std::fprintf(stderr,"failed step %d\n",n);throw std::runtime_error("long throughflow");}
  for(std::size_t j=0;j<2;++j){worst=std::max(worst,std::abs(a.area[j]-b.area[j])/1e-6);worst=std::max(worst,std::abs(a.effort[j]-b.effort[j])/1000);worst=std::max(worst,std::abs(a.mass[j]-b.mass[j])/1e-8);}for(std::size_t j=0;j<3;++j){worst=std::max(worst,std::abs(a.flow[j]-b.flow[j])/1e-4);}work=std::max(work,std::abs(continuous.last.energy_residual));
 }
 if(worst>1e-7 || work>1e-14){throw std::runtime_error("restart invariance");}std::printf("PASS continuous/restarted 10 ms throughflow: maximum scaled difference %.17g, step work residual %.17g J\n",worst,work);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
