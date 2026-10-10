#include "ducklips/slit_thermal.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
void fixture(){double c,a,dt,q;
    while(std::scanf("%lf %lf %lf %lf",&c,&a,&dt,&q)==4){ducklips::SlitThermal model;require(model.prepare({},c,a),"prepare");
        for(double& z:model.charge){require(std::scanf("%lf",&z)==1,"input");}
        require(model.advance(q,dt),"advance");std::printf("%.17g",q);for(double z:model.charge){std::printf(" %.17g",z);}
        std::printf(" %.17g %.17g %.17g\n",model.last.dissipated_energy,model.last.energy_residual,model.last.stored_air_residual);
    }
}
void controls(){double worst=0,mass=0;std::size_t count=0;
    for(double c:{1e-13,1e-10,1e-7})for(double a:{1e-6,1e-4,.01})for(double dt:{1e-10,1e-5,10.}){
        ducklips::SlitThermal model;require(model.prepare({},c,a),"prepare");double q=700*c;
        require(model.equilibrate(-200),"equilibrate");const double initial=.5*q*q/c+model.energy(),stored=q+model.stored_air();
        require(model.advance(q,dt),"advance");require(model.last.dissipated_energy>=0 && model.last.dissipated_energy<=initial*(1+1e-12),"passivity");
        worst=std::max(worst,std::abs(model.last.energy_residual)/initial);mass=std::max(mass,std::abs(q+model.stored_air()-stored)/(c*700));++count;
    }
    require(worst<1e-11 && mass<1e-14,"work and air");
    ducklips::SlitThermal model;require(model.prepare({},1e-10,.001),"sealed prepare");double q=700e-10;
    require(model.advance(q,1e5),"sealed cool");require(std::abs(q/1e-10-500)<1e-8,"isothermal pressure");
    require(std::abs(q+model.stored_air()-700e-10)<1e-21,"sealed conserved air");
    const double e=model.energy(),m=model.stored_air();double work=0;
    require(model.change_geometry(2e-10,.002,work) && model.stored_air()==m,"moving retained charge");
    require(std::abs(work+.5*e)<1e-20,"moving energy work");
    const auto saved=model;const double primary=q;
    require(!model.advance(q,-1) && q==primary && model.charge==saved.charge,"advance rollback");
    require(!model.change_geometry(0,.001,work) && model.compliance()==saved.compliance(),"no empty storage floor");
    require(!model.prepare({},1e-10,0) && model.charge==saved.charge,"prepare rollback");
    require(!model.equilibrate(std::numeric_limits<double>::infinity()) && model.charge==saved.charge,"equilibrium rollback");
    std::printf("PASS finite-gap thermal: %zu controls, relative energy %.17g, relative stored-air %.17g; sealed equilibrium and geometry work\n",count,worst,mass);
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{require(argc==1,"arguments");controls();}return 0;}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
