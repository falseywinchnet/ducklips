#include "ducklips/slit_viscous.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
void fixture(){double l,g,a,dt,p;
    while(std::scanf("%lf %lf %lf %lf %lf",&l,&g,&a,&dt,&p)==5){ducklips::SlitViscous model;require(model.prepare({},l,g,a),"prepare");
        for(double& z:model.momentum){require(std::scanf("%lf",&z)==1,"input");}
        require(model.advance(p,dt),"advance");std::printf("%.17g",p);for(double z:model.momentum){std::printf(" %.17g",z);}
        std::printf(" %.17g %.17g\n",model.last.dissipated_energy,model.last.energy_residual);
    }
}
void controls(){double worst=0;std::size_t count=0;
    for(double l:{1.,100.,1e5})for(double fraction:{1e-4,.5,1.})for(double a:{1e-6,.001,.01})for(double dt:{1e-10,1e-5,1e4}){
        ducklips::SlitViscous model;const double g=fraction/l;require(model.prepare({},l,g,a),"prepare");double p=.001/g;
        require(model.equilibrate(-.0003),"equilibrate");const double initial=.5*g*p*p+model.energy();
        require(model.advance(p,dt),"advance");require(model.last.dissipated_energy>=0 && model.last.dissipated_energy<=initial*(1+1e-11),"passivity");
        worst=std::max(worst,std::abs(model.last.energy_residual)/initial);++count;
    }
    require(worst<1e-11,"work identity");
    ducklips::SlitViscous model;require(model.prepare({},100.,.01,.001) && model.equilibrate(.001),"profile prepare");
    require(std::abs(model.energy()-.1*100*.001*.001)<1e-18,"profile inertia adds one fifth at DC");
    const double energy=model.energy();const auto memory=model.momentum;double work;
    require(model.change_geometry(200.,.005,.002,work) && model.momentum==memory,"retain momentum on geometry change");
    require(std::abs(work+.5*energy)<1e-18,"geometry work");
    const auto saved=model;double p=.02;
    require(!model.change_geometry(0,.005,.002,work) && model.momentum==saved.momentum && model.segment_inertance()==saved.segment_inertance(),"empty path requires explicit retirement");
    require(!model.advance(p,-1) && p==.02 && model.momentum==saved.momentum,"rollback");
    require(!model.prepare({},100.,-.01,.001) && model.momentum==saved.momentum,"prepare rollback");
    require(model.advance(p,1e8) && p==0,"long decay leaves no artificial residue");
    for(double z:model.momentum){require(z==0,"relaxed memory");}
    require(model.prepare({},100.,.01,.001) && model.equilibrate(.001),"closure profile");
    const auto open_memory=model.momentum;const double open_energy=model.energy();
    require(model.change_geometry(100.,0.,.001,work) && work==0 && model.momentum==open_memory,"blocked path retains local air");
    require(!model.equilibrate(.001),"blocked flow cannot equilibrate");
    p=0;require(model.advance(p,1e-4) && p==0 && model.last.blocked_impulse>0,"blocked reaction without flow");
    require(model.energy()>0 && model.energy()<open_energy && std::abs(model.last.energy_residual)<1e-18,"blocked profile decays physically");
    const double remaining=model.energy();
    require(model.change_geometry(0,0,0,work) && work==-remaining,"absent metric work");
    require(model.energy()==0 && model.advance(p,1.) && p==0,"absent path");
    for(double z:model.momentum){require(z==0,"absence retires history");}
    require(model.change_geometry(100.,.01,.001,work) && work==0 && model.advance(p,1e-4) && p==0,"reopening cannot resurrect history");
    require(model.equilibrate(.001),"open branch profile");
    const double before_absence=model.energy();
    require(model.change_geometry(0.,.01,0.,work) && work==-before_absence,"zero-length site on open branch");
    p=.1;require(model.advance(p,.001) && p==.1 && model.energy()==0 && model.branch_conductance()==.01,"absent site leaves branch momentum alone");
    require(model.equilibrate(.001),"open flow needs no absent profile");
    for(double z:model.momentum){require(z==0,"open absence retires history");}
    require(model.change_geometry(100.,.01,.001,work) && work==0,"new site starts without old energy");
    for(double z:model.momentum){require(z==0,"new shoulder has no resurrected profile");}
    std::printf("PASS finite-gap viscosity: %zu controls, relative energy %.17g; shared inertia, geometry work and rollback\n",count,worst);
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{require(argc==1,"arguments");controls();}return 0;}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
