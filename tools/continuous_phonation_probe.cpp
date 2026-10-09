#include "ducklips/continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {
struct Ledger {
 double physical{},numerical{},work{},net{},minimum_dt{std::numeric_limits<double>::infinity()},first_closed{-1},first_reopened{-1};int segments{},transitions{},closures{},reopenings{};
 void add(const Ledger& s){physical+=s.physical;numerical+=s.numerical;work+=s.work;net+=s.net;minimum_dt=std::min(minimum_dt,s.minimum_dt);segments+=s.segments;transitions+=s.transitions;closures+=s.closures;reopenings+=s.reopenings;if(first_closed<0){first_closed=s.first_closed;}if(first_reopened<0){first_reopened=s.first_reopened;}}
 void step(const ducklips::ContinuousGlottis& source,double dt){physical+=source.last.fluid_loss+source.last.tissue_loss+source.last.load_work;numerical+=source.last.numerical_loss;work+=source.last.reservoir_work+source.last.external_work;net+=dt*(source.last.flow[0]-source.last.flow[2]);minimum_dt=std::min(minimum_dt,dt);++segments;}
};
// Trial copies make a rejected subdivision transactional. No physical force,
// opening, or time is clamped. Dissipation is used only on a verified change
// of exterior-port constraints; other rejected steps are bisected.
bool segment(ducklips::ContinuousGlottis& source,double dt,double time,double pressure,int remaining,Ledger& ledger){
 ducklips::ContinuousGlottis trial=source;
 if(trial.advance(dt,pressure,0.,0.,{},ducklips::ContinuousGlottalUpdate::symmetric)){ledger.step(trial,dt);source=trial;return true;}
 trial=source;ducklips::ContinuousGlottalState before,after;
 if(source.state(before) && trial.advance(dt,pressure,0.,0.,{},ducklips::ContinuousGlottalUpdate::dissipative) && trial.state(after) && ((before.area[0]>0)!=(after.area[0]>0) || (before.area[1]>0)!=(after.area[1]>0))){
  Ledger accepted;accepted.step(trial,dt);accepted.transitions=1;
  for(std::size_t j=0;j<2;++j){if(before.area[j]>0 && after.area[j]<=0){++accepted.closures;accepted.first_closed=time+dt;}if(before.area[j]<=0 && after.area[j]>0){++accepted.reopenings;accepted.first_reopened=time+dt;}}
  ledger.add(accepted);source=trial;return true;
 }
 if(remaining==0){return false;}trial=source;Ledger halves;
 if(!segment(trial,dt/2,time,pressure,remaining-1,halves) || !segment(trial,dt/2,time+dt/2,pressure,remaining-1,halves)){return false;}source=trial;ledger.add(halves);return true;
}
}
int main(int argc,char** argv){try{
 const int steps=argc>1?std::stoi(argv[1]):1000;const double duration=argc>2?std::stod(argv[2]):.01,pressure=argc>3?std::stod(argv[3]):700.,rest=argc>4?std::stod(argv[4]):3e-6;if(steps<=0 || !std::isfinite(duration) || duration<=0){throw std::runtime_error("arguments");}
 ducklips::FoldParameters f;ducklips::FoldState initial;initial.rest_area=rest;std::array<double,2> effort{};std::array<double,3> flow{};
 if(argc>5){double a,b,v0,v1;if(std::string(argv[5])!="--state" || std::scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf",&effort[0],&effort[1],&a,&b,&flow[0],&flow[1],&flow[2],&v0,&v1)!=9){throw std::runtime_error("initial state");}initial.lower_position.fill((a-rest)/(2*f.length));initial.upper_position.fill((b-rest)/(2*f.length));initial.lower_velocity.fill(v0/(2*f.length));initial.upper_velocity.fill(v1/(2*f.length));}
 const std::string mode=argc>6?argv[6]:"dissipative";if(mode!="dissipative" && mode!="symmetric" && mode!="transition" && mode!="adaptive"){throw std::runtime_error("mode");}
 ducklips::ContinuousGlottis source;if(!source.prepare(f,{},1.,initial,effort,flow,1.,ducklips::ContourPressureProfile::linear_span)){throw std::runtime_error("prepare");}ducklips::ContinuousGlottalState first,s;if(!source.state(first)){throw std::runtime_error("state");}const double dt=duration/steps;Ledger total;double max_work=0,max_mass=0;
 for(int n=0;n<=steps;++n){if(!source.state(s)){throw std::runtime_error("state");}max_work=std::max(max_work,std::abs(s.energy-first.energy+total.physical+total.numerical-total.work));max_mass=std::max(max_mass,std::abs(s.mass[0]+s.mass[1]-first.mass[0]-first.mass[1]-total.net));
  if(n%10==0 || n==steps){std::printf("{\"step\":%d,\"time\":%.17g,\"area\":[%.17g,%.17g],\"effort\":[%.17g,%.17g],\"flow\":[%.17g,%.17g,%.17g],\"area_rate\":[%.17g,%.17g],\"energy\":%.17g,\"physical_loss\":%.17g,\"numerical_loss\":%.17g,\"work\":%.17g}\n",n,n*dt,s.area[0],s.area[1],s.effort[0],s.effort[1],s.flow[0],s.flow[1],s.flow[2],f.length*(s.folds.lower_velocity[0]+s.folds.lower_velocity[1]),f.length*(s.folds.upper_velocity[0]+s.folds.upper_velocity[1]),s.energy,total.physical,total.numerical,total.work);std::fflush(stdout);}
  if(n==steps){std::fprintf(stderr,"{\"complete\":true,\"mode\":\"%s\",\"steps\":%d,\"duration\":%.17g,\"pressure\":%.17g,\"rest_area\":%.17g,\"maximum_work_residual\":%.17g,\"maximum_mass_residual\":%.17g,\"segments\":%d,\"minimum_step_s\":%.17g,\"dissipative_transitions\":%d,\"closures\":%d,\"reopenings\":%d,\"first_closed_s\":%.17g,\"first_reopened_s\":%.17g}\n",mode.c_str(),steps,duration,pressure,rest,max_work,max_mass,total.segments,total.minimum_dt,total.transitions,total.closures,total.reopenings,total.first_closed,total.first_reopened);return 0;}
  Ledger accepted;bool advanced;
  if(mode=="transition" || mode=="adaptive"){advanced=segment(source,dt,n*dt,pressure,mode=="adaptive"?12:0,accepted);}else{advanced=source.advance(dt,pressure,0.,0.,{},mode=="symmetric"?ducklips::ContinuousGlottalUpdate::symmetric:ducklips::ContinuousGlottalUpdate::dissipative);if(advanced){accepted.step(source,dt);}}
  if(!advanced){std::fprintf(stderr,"{\"complete\":false,\"mode\":\"%s\",\"failed_step\":%d,\"time\":%.17g,\"area\":[%.17g,%.17g],\"effort\":[%.17g,%.17g],\"maximum_work_residual\":%.17g,\"maximum_mass_residual\":%.17g}\n",mode.c_str(),n,n*dt,s.area[0],s.area[1],s.effort[0],s.effort[1],max_work,max_mass);return 2;}total.add(accepted);
 }
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
