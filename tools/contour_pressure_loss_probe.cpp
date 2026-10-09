#include "ducklips/contour_pressure_loss.hpp"
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {
void require(bool b,const char* text){if(!b){throw std::runtime_error(text);}}
template<std::size_t N> void array(const std::array<double,N>& x){std::printf("[");for(std::size_t j=0;j<N;++j){std::printf("%s%.17g",j?",":"",x[j]);}std::printf("]");}
}
int main(int argc,char** argv){try{
    ducklips::ContourPressureKinetics kinetics;ducklips::ContourPressureDissipation model;require(kinetics.configure({},{}) && model.configure({},{}),"configure");
    if(argc==2 && std::string(argv[1])=="--conductance"){
        double d,g;while(std::scanf("%lf",&d)==1){require(ducklips::slit_bgk_conductance(d,g),"conductance");std::printf("%.17g\n",g);}return 0;
    }
    if(argc==2 && std::string(argv[1])=="--batch"){
        std::array<double,2> a,u;while(std::scanf("%lf%lf%lf%lf",&a[0],&a[1],&u[0],&u[1])==4){ducklips::ContourPressureFlow s;ducklips::ContourPressureLoss loss;require(kinetics.prepare(a,u,s),"batch flow");require(model.prepare(s,loss),"batch loss");std::printf("{\"resistance\":");array(loss.resistance);std::printf("}\n");}return 0;
    }
    require(argc==1,"arguments");double g=123;for(double d:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){require(!ducklips::slit_bgk_conductance(d,g) && g==123,"invalid conductance transactional");}
    require(ducklips::slit_bgk_conductance(1e-300,g) && g>0,"molecular limit");require(ducklips::slit_bgk_conductance(1e300,g) && std::abs(g/1e300-1./6)<1e-15,"continuum limit");
    for(double d:{1e-8,100.}){double l,r;require(ducklips::slit_bgk_conductance(d*(1-1e-10),l) && ducklips::slit_bgk_conductance(d*(1+1e-10),r) && std::abs(l/r-1)<3e-10,"continuous asymptotic joins");}
    ducklips::ContourPressureFlow state;ducklips::ContourPressureLoss loss;double power;
    for(const std::array<double,2> area:{std::array<double,2>{-5.2e-6,1e-6},std::array<double,2>{3e-6,-1e-6},std::array<double,2>{0.,1e-6},std::array<double,2>{1e-6,0.}}){require(kinetics.prepare(area,{700.,700.},state),"contact flow");require(model.prepare(state,loss),"contact loss");std::array<double,5> rate={1e-5,2e-5,3e-5,-.001,.002};rate[area[0]<=0?0:2]=0.;require(model.power(state,loss,rate,power) && power>0,"finite positive contact power");rate[area[0]<=0?0:2]=1e-5;const double previous=power;require(!model.power(state,loss,rate,power) && power==previous,"closed port rejection");}
    ducklips::AirProperties air;air.viscosity=1e-12;require(model.configure({},air),"continuum configure");require(kinetics.prepare({3e-6,3e-6},{0.,0.},state) && model.prepare(state,loss),"uniform tube");require(model.power(state,loss,{1e-5,1e-5,1e-5,0.,0.},power),"uniform power");
    const double expected=12*air.viscosity*.014*.014*.003/(3e-6*3e-6*3e-6)*1e-10;require(std::abs(power/expected-1)<1e-8,"classical Poiseuille limit");
    air.viscosity=0.;require(model.configure({},air) && model.prepare(state,loss) && model.power(state,loss,{1e-5,1e-5,1e-5,0.,0.},power) && power==0.,"inviscid limit");
    std::printf("PASS constitutive limits and joins, finite moving-contact loss, closed-port rejection, uniform Poiseuille and inviscid limits\n");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
