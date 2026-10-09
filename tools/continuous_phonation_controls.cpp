#include "ducklips/continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
int main(int argc,char** argv){try{
 if(argc!=2){throw std::runtime_error("pressure");}const double pressure=std::stod(argv[1]);if(pressure!=700 && pressure!=1400){throw std::runtime_error("pressure");}
 // Independently solved full-span equilibria with a small modal perturbation.
 const std::array<double,9> cases[]={{119.58409897364274,-3.5075776989433929,3.7310808400520681e-06,4.2220140855856485e-06,0.0001185640176737039,0.00011856402658241322,0.00011856404618802033,1.0573390035818121e-09,-2.0203048066616277e-08},{112.19146094389872,-3.1389466651766678,3.6867231667931326e-06,4.1489365978136318e-06,0.00017465739767136106,0.00017465740705842998,0.00017465742312911634,-1.9755208388722954e-09,-1.5171180853447539e-08}};const auto& y=cases[pressure==700?0:1];ducklips::FoldParameters f;ducklips::FoldState initial;initial.rest_area=3e-6;initial.lower_position.fill((y[2]-initial.rest_area)/(2*f.length));initial.upper_position.fill((y[3]-initial.rest_area)/(2*f.length));initial.lower_velocity.fill(y[7]/(2*f.length));initial.upper_velocity.fill(y[8]/(2*f.length));ducklips::ContinuousGlottis source;
 if(!source.prepare(f,{},1.,initial,{y[0],y[1]},{y[4],y[5],y[6]},1.,ducklips::ContourPressureProfile::linear_span)){throw std::runtime_error("prepare");}double worst=0,scale=0;int limited=0;
 for(int n=0;n<40;++n){if(!source.advance(2.5e-6,pressure,0.,0.,{},ducklips::ContinuousGlottalUpdate::dissipative)){throw std::runtime_error("near-steady advance");}worst=std::max(worst,std::abs(source.last.energy_residual));scale=std::max(scale,source.last.energy_roundoff_scale);limited+=source.last.roundoff_limited_gradient?1:0;}
 ducklips::ContinuousGlottalState final;if(!source.state(final) || limited==0 || worst>1e-15 || scale<=0 || scale>1e-15){throw std::runtime_error("roundoff accounting");}std::printf("PASS %.0f Pa near-steady motion, %d roundoff-limited steps, residual %.17g J scale %.17g J\n",pressure,limited,worst,scale);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
