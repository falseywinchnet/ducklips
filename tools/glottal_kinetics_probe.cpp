#include "ducklips/glottal_kinetics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
bool close(double a,double b){return std::abs(a-b)<=1e-12*std::max(std::abs(a),std::abs(b));}
void values(const char* name,const auto& v){std::printf("\"%s\":[",name);for(std::size_t i=0;i<v.size();++i){std::printf("%s%.17g",i?",":"",v[i]);}std::putchar(']');}
}
int main(int argc,char** argv){try{
    ducklips::GlottalKinetics k;const ducklips::FoldParameters f;const ducklips::AirProperties air;
    require(k.configure(f,air,1.),"configure kinetics");
    if(argc==2 && std::string(argv[1])=="--pairs"){
        std::array<double,2> m0{},m1{};std::array<double,3> p0{},p1{};
        while(std::scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&m0[0],&m0[1],&m1[0],&m1[1],&p0[0],&p0[1],&p0[2],&p1[0],&p1[1],&p1[2])==10){
            ducklips::GlottalKineticState a,b;ducklips::GlottalKineticTransfer t;
            if(!k.evaluate(m0,p0,a) || !k.evaluate(m1,p1,b) || !k.transfer(m0,m1,p0,p1,t)){std::puts("{\"accepted\":false}");continue;}
            std::printf("{\"accepted\":true,\"energy0\":%.17g,\"energy1\":%.17g,\"chain_residual\":%.17g,",a.energy,b.energy,t.chain_residual);
            values("metric0",a.inverse_inertance);std::putchar(',');values("instant_effort0",a.mass_effort);std::putchar(',');values("flow",t.flow);std::putchar(',');values("mass_effort",t.mass_effort);std::puts("}");
        }return 0;
    }
    require(argc==1,"unknown arguments");
    const double area=3e-6,flow=5e-5;const std::array<double,2> mass={area*f.lower_thickness,area*f.upper_thickness};
    ducklips::GlottalKineticState zero,s;require(k.evaluate(mass,{},zero),"zero state");std::array<double,3> momentum{};
    for(std::size_t j=0;j<3;++j){momentum[j]=flow/zero.inverse_inertance[j];}
    require(k.evaluate(mass,momentum,s),"uniform state");
    const double velocity=flow/area,head=.5*air.density*velocity*velocity;
    require(close(s.mass_effort[0],head) && close(s.mass_effort[1],head),"uniform Bernoulli kinetic head");
    require(close(s.energy,.5*air.density*area*(f.lower_thickness+f.upper_thickness)*velocity*velocity),"uniform kinetic energy");
    std::array<double,3> r{},b{},r2{},b2{};require(k.losses(mass,mass,r,b),"uniform losses");
    const std::array<double,2> doubled={2*mass[0],2*mass[1]},squeezed={.5*mass[0],.5*mass[1]};
    require(k.losses(doubled,doubled,r2,b2),"larger losses");
    for(std::size_t j=0;j<3;++j){require(close(r2[j],r[j]/8),"parallel-plate area scaling");}require(close(b2[2],b[2]/4),"kinetic-head area scaling");
    require(k.losses(mass,squeezed,r2,b2),"compressed losses");
    for(std::size_t j=0;j<3;++j){require(close(r2[j],2*r[j]),"fixed-mass compression scaling");}
    const auto saved=s;require(!k.evaluate({0.,mass[1]},momentum,s) && s.energy==saved.energy,"invalid mass changed result");
    std::puts("PASS uniform kinetic energy/head, geometric and compression loss scaling, invalid mass rejection");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
