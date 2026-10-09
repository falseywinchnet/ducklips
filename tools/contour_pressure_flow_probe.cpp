#include "ducklips/contour_pressure_flow.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool b,const char* text){if(!b){throw std::runtime_error(text);}}
template<std::size_t N> void array(const std::array<double,N>& x){std::printf("[");for(std::size_t j=0;j<N;++j){std::printf("%s%.17g",j?",":"",x[j]);}std::printf("]");}
}
int main(int argc,char** argv){try{
    ducklips::ContourPressureKinetics model;require(model.configure({},{}),"configure");
    if(argc==2 && std::string(argv[1])=="--batch"){
        std::array<double,2> a,u;while(std::scanf("%lf%lf%lf%lf",&a[0],&a[1],&u[0],&u[1])==4){ducklips::ContourPressureFlow s;require(model.prepare(a,u,s),"batch prepare");
            std::printf("{\"inertia\":");array(s.inertia);std::printf(",\"effort_rate\":");array(s.effort_rate);std::printf(",\"support\":");array(s.support);std::printf(",\"contact_loss\":");array(s.contact_loss);std::printf(",\"samples\":[");
            for(int j=1;j<=5;++j){const double z=s.support[0]+(s.support[1]-s.support[0])*j/6.;ducklips::ContourFlowSample p;require(model.sample(s,z,p),"sample");std::printf("%s{\"z\":%.17g,\"area\":%.17g,\"density\":%.17g,\"flow\":",j==1?"":",",z,p.area,p.relative_density);array(p.flow);std::printf("}");}std::printf("]}\n");}return 0;
    }
    require(argc==1,"arguments");ducklips::ContourPressureFlow state;require(model.prepare({3e-6,3e-6},{700.,700.},state),"uniform tube");
    const double r=std::pow(1+.4*700/139650.,2.5),q=1e-5;double energy=0;require(model.energy(state,{q,q,q,0.,0.},energy),"throughflow");const double expected=.5*1.14*.003/(r*3e-6)*q*q;require(std::abs(energy/expected-1)<1e-11,"uniform exact inertia");
    double maximum_flow=0;for(int j=0;j<=20;++j){ducklips::ContourFlowSample p;require(model.sample(state,.003*j/20.,p),"uniform field");maximum_flow=std::max(maximum_flow,std::abs(p.flow[0]+p.flow[1]+p.flow[2]-1));}require(maximum_flow<1e-12,"constant throughflow");
    for(const std::array<double,2> area:{std::array<double,2>{-5.2e-6,1e-6},std::array<double,2>{3e-6,-1e-6}}){require(model.prepare(area,{-100.,900.},state),"contact field");std::array<double,5> rate={1e-5,2e-5,3e-5,-.001,.002};rate[area[0]<=0?0:2]=0.;require(model.energy(state,rate,energy) && energy>0,"contact energy");rate[area[0]<=0?0:2]=1e-5;const double previous=energy;require(!model.energy(state,rate,energy) && energy==previous,"closed flow rejection");}
    std::printf("PASS uniform kinetic energy, local throughflow %.17g, both contact supports, closed-flow rejection\n",maximum_flow);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
