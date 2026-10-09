#include "ducklips/contour_kinetics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
namespace {void require(bool v,const char* text){if(!v){throw std::runtime_error(text);}}}
int main(int argc,char** argv){try{
    ducklips::ContourKinetics model;require(model.configure({},{}),"configure");
    if(argc==2 && std::string(argv[1])=="--pairs"){
        std::array<double,2> m0,a0,m1,a1;std::array<double,3> p0,p1;
        while(std::scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&m0[0],&m0[1],&a0[0],&a0[1],&p0[0],&p0[1],&p0[2],&m1[0],&m1[1],&a1[0],&a1[1],&p1[0],&p1[1],&p1[2])==14){
            ducklips::ContourKineticTransfer k;require(model.transfer(m0,a0,p0,m1,a1,p1,k,ducklips::ContourUpdate::dissipative),"batch transfer");
            std::printf("{\"flow\":[%.17g,%.17g,%.17g],\"effort\":[%.17g,%.17g,%.17g,%.17g],\"numerical_loss\":%.17g,\"energy_change\":%.17g,\"chain_residual\":%.17g}\n",k.flow[0],k.flow[1],k.flow[2],k.effort[0],k.effort[1],k.effort[2],k.effort[3],k.numerical_loss,k.energy_change,k.chain_residual);
        }return 0;
    }
    require(argc==1,"arguments");const std::array<double,2> mass={7.5e-9,1.5e-9},area={3e-6,3e-6};double max_balance=0;
    for(std::size_t j=0;j<3;++j){std::array<double,3> p0{};p0[j]=.002;ducklips::ContourKineticState k;require(model.evaluate(mass,area,p0,k),"evaluate");
        for(double stiffness:{1e-8,1e-4,1.,100.,1e12}){const double dt=stiffness/(k.resistance[j]*k.inverse_inertance[j]);auto p1=p0;p1[j]/=1+stiffness;
            ducklips::ContourKineticTransfer t;require(model.transfer(mass,area,p0,mass,area,p1,t,ducklips::ContourUpdate::dissipative),"transfer");
            const double physical=dt*k.resistance[j]*t.flow[j]*t.flow[j];max_balance=std::max(max_balance,std::abs(t.energy_change+physical+t.numerical_loss)/k.energy);
            require(t.numerical_loss>=0,"negative discretization loss");if(stiffness==1e12){require(std::abs(p1[j]/p0[j])<1.01e-12,"stiff mode not extinguished");}
        }
    }
    require(max_balance<1e-14,"fixed-passage dissipation identity");
    const std::array<double,2> sealed={-1e-7,3e-6};ducklips::ContourKineticState before;require(model.evaluate(mass,area,{.002,0,0},before),"preseal");ducklips::ContourKineticTransfer seal;
    require(model.transfer(mass,area,{.002,0,0},mass,sealed,{},seal,ducklips::ContourUpdate::dissipative),"seal transfer");
    require(seal.flow[0]==0 && seal.numerical_loss==0 && std::abs(seal.chain_residual)<1e-20 && std::abs(seal.energy_change+before.energy)<1e-20,"sealing work not accounted");
    ducklips::ContourKineticTransfer reopen;require(model.transfer(mass,sealed,{},mass,area,{.002,0,0},reopen,ducklips::ContourUpdate::dissipative),"reopening transfer");
    require(reopen.flow[0]>0 && std::abs(reopen.numerical_loss-reopen.energy_change)<1e-20,"reopening energy/loss split");
    std::printf("PASS fixed-passage work %.17g relative; stiff decay, sealing geometry work, reopening numerical loss\n",max_balance);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
