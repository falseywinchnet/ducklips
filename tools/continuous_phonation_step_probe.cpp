#include "ducklips/continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
namespace {
void state(const ducklips::ContinuousGlottalState& s,double width){std::printf("[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",s.effort[0],s.effort[1],s.area[0],s.area[1],s.flow[0],s.flow[1],s.flow[2],width*(s.folds.lower_velocity[0]+s.folds.lower_velocity[1]),width*(s.folds.upper_velocity[0]+s.folds.upper_velocity[1]));}
}
int main(){int steps;double duration,pressure,rest,u0,u1,a,b,q0,q1,q2,v0,v1;
 while(std::scanf("%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&steps,&duration,&pressure,&rest,&u0,&u1,&a,&b,&q0,&q1,&q2,&v0,&v1)==13){ducklips::FoldParameters f;ducklips::FoldState initial;initial.rest_area=rest;initial.lower_position.fill((a-rest)/(2*f.length));initial.upper_position.fill((b-rest)/(2*f.length));initial.lower_velocity.fill(v0/(2*f.length));initial.upper_velocity.fill(v1/(2*f.length));ducklips::ContinuousGlottis model;ducklips::ContinuousGlottalState first,last;bool ok=steps>0 && duration>0 && model.prepare(f,{},1.,initial,{u0,u1},{q0,q1,q2},1.,ducklips::ContourPressureProfile::linear_span) && model.state(first);double residual=0,roundoff=0;int completed=0,limited=0;
 for(;ok && completed<steps;++completed){ok=model.advance(duration/steps,pressure,0.,0.,{},ducklips::ContinuousGlottalUpdate::dissipative);residual=std::max(residual,std::abs(model.last.energy_residual));roundoff=std::max(roundoff,model.last.energy_roundoff_scale);limited+=model.last.roundoff_limited_gradient?1:0;}ok=ok && model.state(last);if(!ok){std::printf("{\"accepted\":false,\"completed\":%d}\n",completed);continue;}std::printf("{\"accepted\":true,\"initial\":");state(first,f.length);std::printf(",\"final\":");state(last,f.length);std::printf(",\"maximum_step_work_residual\":%.17g,\"maximum_energy_roundoff_scale\":%.17g,\"roundoff_limited_steps\":%d}\n",residual,roundoff,limited);std::fflush(stdout);
 }return 0;}
