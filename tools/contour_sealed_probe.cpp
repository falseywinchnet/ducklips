#include "ducklips/contour_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool v,const char* text){if(!v){throw std::runtime_error(text);}}
void vector(const char* label,const double* x,std::size_t n){std::printf("\"%s\":[",label);for(std::size_t i=0;i<n;++i){std::printf("%s%.17g",i?",":"",x[i]);}std::putchar(']');}
}
int main(int argc,char** argv){try{
    if(argc==2 && std::string(argv[1])=="--kernel"){
        ducklips::ContourKinetics model;require(model.configure({},{}),"configure");std::array<double,2> m,a;std::array<double,3> p;
        while(std::scanf("%lf%lf%lf%lf%lf%lf%lf",&m[0],&m[1],&a[0],&a[1],&p[0],&p[1],&p[2])==7){ducklips::ContourKineticState k;if(!model.evaluate(m,a,p,k)){std::puts("{\"accepted\":false}");continue;}
            std::printf("{\"accepted\":true,\"energy\":%.17g,",k.energy);vector("effort",k.effort.data(),4);std::putchar(',');vector("inverse",k.inverse_inertance.data(),3);std::putchar(',');vector("resistance",k.resistance.data(),3);
            const double volume[]={k.geometry.chamber[0].volume,k.geometry.chamber[1].volume};std::putchar(',');vector("volume",volume,2);std::puts("}");}return 0;
    }
    if(argc==2 && std::string(argv[1])=="--pairs"){
        ducklips::ContourKinetics model;require(model.configure({},{}),"configure");std::array<double,2> m0,a0,m1,a1;std::array<double,3> p0,p1;
        while(std::scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&m0[0],&m0[1],&a0[0],&a0[1],&p0[0],&p0[1],&p0[2],&m1[0],&m1[1],&a1[0],&a1[1],&p1[0],&p1[1],&p1[2])==14){
            ducklips::ContourKineticTransfer k;if(!model.transfer(m0,a0,p0,m1,a1,p1,k)){std::puts("{\"accepted\":false}");continue;}
            std::printf("{\"accepted\":true,\"chain_residual\":%.17g,",k.chain_residual);vector("effort",k.effort.data(),4);std::putchar(',');vector("flow",k.flow.data(),3);
            std::putchar(',');vector("volume_gradient_lower",k.volume_gradient[0].data(),2);std::putchar(',');vector("volume_gradient_upper",k.volume_gradient[1].data(),2);std::puts("}");}return 0;
    }
    int steps=100;std::string selected;const bool trajectory=argc==4 && std::string(argv[1])=="--trajectory";
    if(trajectory){steps=std::stoi(argv[2]);selected=argv[3];require(steps>0 && steps<=100000,"step count");}else{require(argc==1,"arguments");}
    if(argc==1){ducklips::FoldParameters f;f.asymmetry=1.17;ducklips::FoldContact contact;require(contact.configure(f,1.),"contact configure");
        for(double a:{-5e-6,-1e-12,0.,1e-12,2e-6}){for(double b:{-5e-6,-1e-12,0.,1e-12,2e-6}){
            const std::array<double,2> first={a,.2*a},last={b,.2*b};ducklips::FoldContactTransfer t,r;require(contact.transfer(first,last,t) && contact.transfer(last,first,r),"contact transfer");
            require(std::abs(t.chain_residual)<1e-20 && t.area_effort==r.area_effort,"contact chain/symmetry");
            ducklips::FoldContactState state;require(contact.evaluate(first,state) && state.energy>=0,"contact energy");if(a>=0){require(state.energy==0 && state.area_effort==std::array<double,2>{},"contact acted across gap");}
        }}
    }
    bool found=false;for(const std::string mode:{"lower","upper","lower_rest","upper_rest"}){if(trajectory && mode!=selected){continue;}found=true;
        const bool lower=mode.starts_with("lower");ducklips::FoldParameters f;f.asymmetry=1.17;ducklips::FoldState initial;initial.rest_area=3e-6;
        (lower?initial.lower_position:initial.upper_position)={-3.1e-6/(2*f.length),-3.1e-6/(2*f.length)};
        ducklips::ContourGlottalStorage source;require(source.prepare(f,{},1.,initial,{100.,50.}),"sealed prepare");source.momentum=lower?std::array<double,3>{0.,.0002,.0001}:std::array<double,3>{.0001,.0002,0.};
        const std::array<double,4> force=lower?std::array<double,4>{-.025,-.025,0.,0.}:std::array<double,4>{0.,0.,-.008,-.008};
        double e0=0;require(source.energy(e0),"energy0");const double mass0=source.state.mass[0]+source.state.mass[1];double work=0,loss=0,net=0,max_work=0,max_mass=0,max_contact=0;const double dt=.0002/double(steps);
        for(int i=0;i<steps;++i){const double rest=3e-6-(mode.ends_with("_rest")?2e-8*double(i+1)/double(steps):0.);require(source.advance(dt,700.,0.,0.,rest,force),"sealed advance");const auto& r=source.last;
            require(r.flow[lower?0:2]==0 && source.momentum[lower?0:2]==0,"sealed passage leaked");
            work+=r.reservoir_work-r.load_work+r.actuator_work+r.external_work;loss+=r.passage_loss+r.damping_loss;net+=dt*(r.flow[0]-r.flow[2]);
            double e=0;require(source.energy(e),"energy1");max_work=std::max(max_work,std::abs(e-e0+loss-work));max_mass=std::max(max_mass,std::abs(source.state.mass[0]+source.state.mass[1]-mass0-net));max_contact=std::max(max_contact,std::abs(r.contact_area_effort[lower?0:1]));}
        require(max_work<1e-13 && max_mass<1e-19 && max_contact>0,"sealed work/mass/contact");
        if(trajectory){std::putchar('{');vector("mass",source.state.mass.data(),2);std::putchar(',');vector("momentum",source.momentum.data(),3);const auto& s=source.state.folds;
            const double x[]={s.lower_position[0],s.lower_position[1],s.upper_position[0],s.upper_position[1]},v[]={s.lower_velocity[0],s.lower_velocity[1],s.upper_velocity[0],s.upper_velocity[1]};
            std::putchar(',');vector("position",x,4);std::putchar(',');vector("velocity",v,4);std::printf(",\"maximum_work_balance\":%.17g,\"maximum_mass_balance\":%.17g}\n",max_work,max_mass);
        }else{std::printf("PASS %s sealed drive: work %.17g J mass %.17g m3 contact %.17g Pa m\n",mode.c_str(),max_work,max_mass,max_contact);}
        const auto saved=source;source.momentum[lower?0:2]=.001;const auto damaged=source;require(!source.advance(dt,700.,0.,0.,3e-6,force),"stale sealed momentum accepted");require(source.state.mass==damaged.state.mass && source.momentum==damaged.momentum && source.last.energy_residual==damaged.last.energy_residual,"rejection mutated state");source=saved;
    }
    require(found,"unknown mode");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
