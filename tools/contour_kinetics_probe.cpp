#include "ducklips/contour_kinetics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <random>
using namespace ducklips;
namespace {void require(bool v,const char* text){if(!v){throw std::runtime_error(text);}}
void output(const ContourKineticState& k){std::printf("{\"energy\":%.17g,\"effort\":[",k.energy);for(int i=0;i<4;++i){std::printf("%s%.17g",i?",":"",k.effort[i]);}
    std::printf("],\"inverse\":[%.17g,%.17g,%.17g],\"resistance\":[%.17g,%.17g,%.17g]}\n",k.inverse_inertance[0],k.inverse_inertance[1],k.inverse_inertance[2],k.resistance[0],k.resistance[1],k.resistance[2]);}}
int main(int argc,char** argv){try{
    FoldParameters f;AirProperties air;ContourKinetics model;require(model.configure(f,air),"configure");
    if(argc==2 && std::string(argv[1])=="--batch"){std::array<double,2> m,a;std::array<double,3> p;while(std::scanf("%lf%lf%lf%lf%lf%lf%lf",&m[0],&m[1],&a[0],&a[1],&p[0],&p[1],&p[2])==7){ContourKineticState k;require(model.evaluate(m,a,p,k),"batch evaluate");output(k);}return 0;}
    if(argc==2 && std::string(argv[1])=="--pairs"){
        std::array<double,2> m0,a0,m1,a1;std::array<double,3> p0,p1;
        while(std::scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&m0[0],&m0[1],&a0[0],&a0[1],&p0[0],&p0[1],&p0[2],&m1[0],&m1[1],&a1[0],&a1[1],&p1[0],&p1[1],&p1[2])==14){
            ContourKineticTransfer t;require(model.transfer(m0,a0,p0,m1,a1,p1,t),"batch transfer");
            std::printf("{\"effort\":[%.17g,%.17g,%.17g,%.17g],\"flow\":[%.17g,%.17g,%.17g],\"chain_residual\":%.17g}\n",t.effort[0],t.effort[1],t.effort[2],t.effort[3],t.flow[0],t.flow[1],t.flow[2],t.chain_residual);
        }return 0;
    }
    require(argc==1,"arguments");std::mt19937_64 rng(7125);std::uniform_real_distribution<double> random(0.,1.);double max_chain=0,max_metric=0,max_limit=0,max_symmetry=0;
    for(int n=0;n<1000;++n){std::array<double,2> a0,a1,m0,m1;std::array<double,3> p0,p1;
        for(int i=0;i<2;++i){a0[i]=std::pow(10.,-8.+4*random(rng));a1[i]=std::pow(10.,-8.+4*random(rng));m0[i]=a0[i]*(i?f.upper_thickness:f.lower_thickness);m1[i]=a1[i]*(i?f.upper_thickness:f.lower_thickness);}
        if(n<8){a0={3e-6,3e-6};a1={3e-6*(1+double(n)*1e-10),3e-6*(1-double(n)*1e-10)};}
        for(int j=0;j<3;++j){p0[j]=.001*(2*random(rng)-1);p1[j]=.001*(2*random(rng)-1);}
        ContourKineticState k0,k1;require(model.evaluate(m0,a0,p0,k0) && model.evaluate(m1,a1,p1,k1),"evaluate");
        for(int j=0;j<3;++j){double inertia=0;for(int i=0;i<2;++i){inertia+=air.density*k0.geometry.chamber[i].volume/m0[i]*k0.geometry.passage[j].inverse_area[i];}
            max_metric=std::max(max_metric,std::abs(inertia*k0.inverse_inertance[j]-1));}
        ContourKineticTransfer t,r,limit;require(model.transfer(m0,a0,p0,m1,a1,p1,t) && model.transfer(m1,a1,p1,m0,a0,p0,r) && model.transfer(m0,a0,p0,m0,a0,p0,limit),"transfer");
        const double scale=k0.energy+k1.energy;max_chain=std::max(max_chain,std::abs(t.chain_residual)/scale);
        for(int i=0;i<4;++i){max_symmetry=std::max(max_symmetry,std::abs(t.effort[i]-r.effort[i])/std::max(1.,std::abs(t.effort[i])));max_limit=std::max(max_limit,std::abs(limit.effort[i]-k0.effort[i])/std::max(1.,std::abs(k0.effort[i])));}
    }
    std::printf("metric %.17g chain %.17g limit %.17g symmetry %.17g\n",max_metric,max_chain,max_limit,max_symmetry);
    require(max_metric<1e-10,"metric mismatch");require(max_chain<1e-10,"kinetic finite chain");require(max_limit<1e-12,"coincident derivative");require(max_symmetry<1e-10,"time symmetry");
    ContourKineticState sentinel;sentinel.energy=123;require(model.evaluate({1e-9,1e-9},{0,1e-6},{},sentinel) && sentinel.inverse_inertance[0]==0 && sentinel.geometry.chamber[0].volume>0,"sealed cap storage");sentinel.energy=123;
    require(!model.evaluate({1e-9,1e-9},{0,1e-6},{.001,0,0},sentinel) && sentinel.energy==123,"sealed momentum rejection");
    require(!model.evaluate({1e-9,1e-9},{0,0},{},sentinel) && sentinel.energy==123,"empty storage rejection");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
