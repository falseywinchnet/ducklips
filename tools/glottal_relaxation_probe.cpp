#include "ducklips/glottal_relaxation.hpp"
#include "ducklips/contour_kinetics.hpp"
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool good,const char* message){if(!good){throw std::runtime_error(message);}}
void tests(){
    std::array<double,49> a{},r{},b{};std::array<bool,7> active;active.fill(true);
    for(std::size_t i=0;i<7;++i){a[7*i+i]=i+1.;r[7*i+i]=std::pow(10.,static_cast<double>(i)-2)/(i+1.);}
    require(ducklips::glottal_relaxation_bias(a,r,active,.01,b),"diagonal relaxation");
    for(std::size_t i=0;i<7;++i){for(std::size_t j=0;j<7;++j){const double expected=i==j?a[7*i+i]*ducklips::contour_relaxation_weight(.01*a[7*i+i]*r[7*i+i]):0.;require(std::abs(b[7*i+j]-expected)<1e-14,"diagonal weight");}}
    a.fill(0);r.fill(0);for(std::size_t i=0;i<7;++i){a[7*i+i]=1.;}
    const double c=std::cos(.37),s=std::sin(.37);r[0]=c*c+1000*s*s;r[8]=s*s+1000*c*c;r[1]=r[7]=(1-1000)*c*s;
    require(ducklips::glottal_relaxation_bias(a,r,active,.01,b),"rotated relaxation");
    const double low=ducklips::contour_relaxation_weight(.01),high=ducklips::contour_relaxation_weight(10.);
    require(std::abs(b[0]-low*c*c-high*s*s)<1e-14 && std::abs(b[8]-low*s*s-high*c*c)<1e-14
        && std::abs(b[1]-(low-high)*c*s)<1e-14,"rotated weights");
    active[0]=false;a[0]=0.;require(ducklips::glottal_relaxation_bias(a,r,active,.01,b),"eliminated coordinate");
    for(std::size_t i=0;i<7;++i){require(b[i]==0 && b[7*i]==0,"inactive bias");}
    const auto saved=b;r[8]=-1.;r[1]=r[7]=0.;require(!ducklips::glottal_relaxation_bias(a,r,active,.01,b) && b==saved,"negative resistance rollback");
    active.fill(true);a.fill(0);r.fill(0);for(std::size_t i=0;i<7;++i){a[7*i+i]=r[7*i+i]=1.;}
    r[1]=1e-16;require(ducklips::glottal_relaxation_bias(a,r,active,.01,b),"cancelled off-diagonal rounding");
    r[1]=.01;require(!ducklips::glottal_relaxation_bias(a,r,active,.01,b),"material asymmetry");
    std::puts("PASS glottal relaxation: diagonal and rotated modes, eliminated coordinates, invalid resistance rollback");
}
void fixture(){
    double step;while(std::scanf("%lf",&step)==1){
        std::array<double,49>a,r,b;std::array<bool,7>active;
        for(bool& value:active){int flag;require(std::scanf("%d",&flag)==1,"active");value=flag!=0;}
        for(double& value:a){require(std::scanf("%lf",&value)==1,"inverse inertia");}
        for(double& value:r){require(std::scanf("%lf",&value)==1,"resistance");}
        require(ducklips::glottal_relaxation_bias(a,r,active,step,b),"fixture relaxation");
        for(std::size_t j=0;j<49;++j){std::printf("%s%.17g",j?" ":"",b[j]);}std::puts("");
    }
}
}
int main(int argc,char**argv){try{if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{tests();}return 0;}catch(const std::exception&e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
