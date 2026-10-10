#include "ducklips/contour_pressure_storage.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
void array(const std::array<double,2>& x){std::printf("[%.17g,%.17g]",x[0],x[1]);}
void matrix(const std::array<double,4>& x){std::printf("[%.17g,%.17g,%.17g,%.17g]",x[0],x[1],x[2],x[3]);}
void output(const ducklips::ContourPressureState& s){
    std::printf("{\"volume\":%.17g,\"energy\":%.17g,\"coenergy\":%.17g,\"mass\":",s.volume,s.energy,s.coenergy);array(s.mass);
    std::printf(",\"effort\":");array(s.effort);std::printf(",\"area_effort\":");array(s.area_effort);std::printf(",\"compliance\":");matrix(s.compliance);std::printf(",\"mass_area\":");matrix(s.mass_area);std::printf("}\n");
}
}
int main(int argc,char** argv){try{
    ducklips::ContourPressureStorage model;require(model.configure({},{}),"configure");
    if(argc==2 && std::string(argv[1])=="--effort-difference"){
        std::array<double,2> a,u,b,v;
        while(std::scanf("%lf%lf%lf%lf%lf%lf%lf%lf",&a[0],&a[1],&u[0],&u[1],&b[0],&b[1],&v[0],&v[1])==8){
            std::printf("%.17g\n",ducklips::ContourPressureStorage::effort_difference(a,u,b,v));
        }
        return 0;
    }
    if(argc==2 && std::string(argv[1])=="--batch"){
        std::array<double,2> a,u;while(std::scanf("%lf%lf%lf%lf",&a[0],&a[1],&u[0],&u[1])==4){ducklips::ContourPressureState s;require(model.evaluate(a,u,s),"batch evaluate");output(s);}return 0;
    }
    if(argc==2 && std::string(argv[1])=="--recover"){
        std::array<double,2> a,m;while(std::scanf("%lf%lf%lf%lf",&a[0],&a[1],&m[0],&m[1])==4){ducklips::ContourPressureState s;require(model.recover(a,m,{},s),"batch recover");output(s);}return 0;
    }
    if(argc==3 && std::string(argv[1])=="--sweep"){
        const int steps=std::stoi(argv[2]);require(steps>=10 && steps<=100000,"sweep steps");
        ducklips::ContourPressureState initial,state;require(model.evaluate({-4.8e-6,1e-6},{700.,700.},initial),"sweep initial");state=initial;
        double work=0,maximum_mass_error=0,maximum_density=1.;const double da=-.4e-6/steps;
        for(int n=1;n<=steps;++n){const std::array<double,2> a={-4.8e-6+n*da,1e-6};ducklips::ContourPressureState next;
            require(model.recover(a,initial.mass,state.effort,next),"sweep recovery");work+=.5*(state.area_effort[0]+next.area_effort[0])*da;
            for(std::size_t j=0;j<2;++j){maximum_mass_error=std::max(maximum_mass_error,std::abs(next.mass[j]-initial.mass[j]));}
            const double t=-a[0]/(a[1]-a[0]);const double max_mu=std::max((1-t)*next.effort[0]+t*next.effort[1],next.effort[1]);
            maximum_density=std::max(maximum_density,std::pow(1+.4*max_mu/139650.,2.5));state=next;
        }
        const double error=state.energy-initial.energy-work;require(maximum_mass_error<1e-20,"sweep mass conservation");
        std::printf("{\"steps\":%d,\"maximum_mass_error_m3\":%.17g,\"maximum_relative_density\":%.17g,\"work_J\":%.17g,\"work_residual_J\":%.17g,\"initial\":",steps,maximum_mass_error,maximum_density,work,error);output(initial);std::printf(",\"final\":");output(state);std::printf("}\n");return 0;
    }
    require(argc==1,"arguments");double maximum_mass=0,maximum_legendre=0,maximum_force=0;
    const auto distance=ducklips::ContourPressureStorage::effort_difference;
    require(distance({1.,2.},{1000.,-7.},{3.,4.},{0.,2.})==1000.,"open support retains coefficient maximum");
    require(std::abs(distance({-99.,1.},{1000.,0.},{-99.,1.},{0.,0.})-10.)<1e-12,"narrow lower contact field");
    require(std::abs(distance({1.,-99.},{0.,1000.},{1.,-99.},{0.,0.})-10.)<1e-12,"narrow upper contact field");
    require(std::abs(distance({-99.,1.},{1000.,0.},{-9.,1.},{0.,0.})-100.)<1e-12,"union includes wider occupied support");
    require(distance({-99.,1.},{1000.,0.},{1.,-99.},{0.,0.})==1000.,"opposite supports keep both endpoints");
    require(distance({0.,1.},{1000.,0.},{0.,1.},{0.,0.})==1000.,"contact limit included");
    require(std::isinf(distance({-1.,0.},{0.,0.},{1.,1.},{0.,0.})),"empty effort support rejected");
    for(double area:{-5.2e-6,-5.000001e-6,-5e-6,-4.999999e-6,-4.8e-6,0.,3e-6}){
        const std::array<double,2> a={area,1e-6};ducklips::ContourPressureState s,recovered;require(model.evaluate(a,{-100.,900.},s),"cross seam");
        require(model.recover(a,s.mass,{},recovered),"inverse moments");
        require(s.compliance[0]*s.compliance[3]>s.compliance[1]*s.compliance[2],"positive susceptibility");
        maximum_legendre=std::max(maximum_legendre,std::abs(s.energy+s.coenergy-s.effort[0]*s.mass[0]-s.effort[1]*s.mass[1]));
        for(std::size_t j=0;j<2;++j){maximum_mass=std::max(maximum_mass,std::abs(recovered.mass[j]-s.mass[j])/(s.mass[0]+s.mass[1]));
            if(j==0 && area==0){continue;} // Zero-area plateau has distinct one-sided forces.
            auto plus=a,minus=a;const double h=1e-10;plus[j]+=h;minus[j]-=h;ducklips::ContourPressureState p,m;
            require(model.recover(plus,s.mass,s.effort,p) && model.recover(minus,s.mass,s.effort,m),"fixed moment shape perturbation");
            const double derivative=(p.energy-m.energy)/(2*h);maximum_force=std::max(maximum_force,std::abs(derivative-s.area_effort[j]));}
    }
    require(maximum_mass<3e-13 && maximum_legendre<1e-15 && maximum_force<2e-5,"thermodynamic work");
    ducklips::ContourPressureState sentinel;sentinel.energy=17.;require(!model.evaluate({-1e-6,-2e-6},{},sentinel) && sentinel.energy==17.,"empty support rejection");
    require(!model.recover({3e-6,3e-6},{-1.,1.},{},sentinel) && sentinel.energy==17.,"invalid mass rejection");
    require(!model.evaluate({3e-6,3e-6},{-1e9,0.},sentinel) && sentinel.energy==17.,"vacuum domain rejection");
    std::printf("PASS continuous storage across old seam; mass %.17g relative; Legendre %.17g J; fixed-moment shape force %.17g Pa m\n",maximum_mass,maximum_legendre,maximum_force);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
