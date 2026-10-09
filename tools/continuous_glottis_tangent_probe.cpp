#include "ducklips/continuous_glottis.hpp"
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {void require(bool b,const char* m){if(!b){throw std::runtime_error(m);}}
template<std::size_t N> void array(const std::array<double,N>& a){std::printf("[");for(std::size_t i=0;i<N;++i){std::printf("%s%.17g",i?",":"",a[i]);}std::printf("]");}}
int main(int argc,char** argv){try{
 const double dt=argc>1?std::stod(argv[1]):1e-9;require(argc<=2 && dt>0,"dt");ducklips::FoldParameters f;f.asymmetry=1.17;f.lower_damping_ratio=f.upper_damping_ratio=0.;ducklips::AirProperties air;air.viscosity=0.;ducklips::FoldState initial;initial.rest_area=3e-6;initial.lower_position={-3e-5,-1e-6/f.length+3e-5};initial.upper_position={4e-5,1e-6/f.length-4e-5};initial.lower_velocity={.02,-.01};initial.upper_velocity={.03,-.01};
 ducklips::ContinuousGlottis model;require(model.prepare(f,air,1.,initial,{700.,-200.},{8e-5,5e-5,3e-5},0.),"prepare");ducklips::ContinuousGlottalState before,after;require(model.state(before) && model.advance(dt,0.,0.,0.) && model.state(after),"advance");
 std::printf("{\"dt\":%.17g,\"initial_mass\":",dt);array(before.mass);std::printf(",\"initial_momentum\":");array(before.momentum);std::printf(",\"final_momentum\":");array(after.momentum);std::printf(",\"energy_residual_J\":%.17g}\n",model.last.energy_residual);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
