#include "ducklips/radiation_boundary.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {
void require(bool v,const char* text){if(!v){throw std::runtime_error(text);}}
void fixture(){double a,b,g,h,p,q;
    while(std::scanf("%lf %lf %lf %lf %lf %lf",&a,&b,&g,&h,&p,&q)==6){
        ducklips::RadiationBoundary r;require(r.prepare({},a),"prepare");r.momentum=q;double work;
        require(r.change_area(b,work) && r.advance(p,g,h),"advance");
        std::printf("%.17g %.17g %.17g %.17g %.17g %.17g\n",p,r.momentum,r.last.radiated_energy,work,r.last.energy_residual,r.energy());
    }
}
void controls(){double worst=0.;std::size_t count=0;
    for(double a:{1e-12,1e-8,1e-4,.01})for(double g:{0.,1e-9,.002,1e4})for(double h:{1e-12,1e-7,.1}){
        ducklips::RadiationBoundary r;require(r.prepare({},a),"prepare");r.momentum=-.01;double p=.003;
        const double total=p+r.momentum,energy=.5*g*p*p+r.energy();
        require(r.advance(p,g,h),"finite advance");
        require(std::abs(p+r.momentum-total)<1e-17,"momentum transfer");
        require(r.last.radiated_energy>=0 && r.last.radiated_energy<=energy*(1+1e-14),"passivity");
        worst=std::max(worst,std::abs(r.last.energy_residual)/(energy+1e-300));++count;
    }
    require(worst<1e-14,"energy identity");
    ducklips::RadiationBoundary r;require(r.prepare({},1e-4),"prepare closure");r.momentum=.02;
    double p=.003,work=0.;const double energy=r.energy();require(r.change_area(0,work),"closure");
    require(r.momentum==0 && r.energy()==0 && work==-energy,"retire absent state");
    require(r.advance(p,0.,1e-5) && p==.003 && r.last.radiated_energy==0,"sealed boundary");
    require(!r.advance(p,.001,1e-5),"reject contradictory closed aperture");
    require(r.change_area(1e-4,work) && r.momentum==0 && work==0,"no reopening memory");
    const auto before=r;const double old=p;
    require(!r.advance(p,-1.,1e-5) && p==old && r.momentum==before.momentum,"rollback");
    require(!r.change_area(std::numeric_limits<double>::quiet_NaN(),work),"reject area");
    require(!r.prepare({},1e-4,-.7) && r.area()==before.area(),"prepare rollback");
    std::printf("PASS radiation: %zu exact-exchange controls, relative energy residual %.17g; closure and rollback\n",count,worst);
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{require(argc==1,"arguments");controls();}return 0;}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
