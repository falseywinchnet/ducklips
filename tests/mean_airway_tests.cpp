#include "ducklips/mean_airway.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {
std::size_t failures{};
void check(bool valid,const char* message,double value){std::printf("%s %s %.12g\n",valid?"PASS":"FAIL",message,value);if(!valid){++failures;}}
ducklips::TractShape tube(std::size_t n){
    ducklips::TractShape s;s.oral_area.assign(n,1e-4);s.oral_length.assign(n,.175/double(n));
    s.nasal_area.assign(2,1e-4);s.nasal_length.assign(2,.05);return s;
}
double steady(std::size_t n,bool narrow,bool reverse){
    const ducklips::AirProperties air{};ducklips::TractShape s=tube(n);
    if(narrow){for(std::size_t i=n/4;i<n/2;++i){s.oral_area[i]=1e-5;}}
    ducklips::AcousticNetwork layout;layout.prepare(n,2,1,{});layout.assemble(s,air,{});
    ducklips::MeanAirway mean;mean.prepare(layout,s,air,1e-8);
    bool valid=true;
    for(std::size_t i=0;i<1000;++i){valid=mean.advance(1,2.5e-5,reverse?-600:600,1e5,.001).converged&&valid;}
    double resistance=1e5;for(std::size_t i=0;i<n;++i){resistance+=8*ducklips::pi*air.viscosity*s.oral_length[i]/(s.oral_area[i]*s.oral_area[i]);}
    double inverse_area_squared=1/(2.5e-5*2.5e-5);
    if(narrow){inverse_area_squared+=1/(1e-5*1e-5)-1/(1e-4*1e-4);}
    const double expected=ducklips::orifice_flow(1/std::sqrt(inverse_area_squared),reverse?-600:600,resistance,air);
    const double error=std::abs(mean.glottal_flow()-expected);
    check(valid&&error<1e-11,"spatial steady flow matches summed viscous and contraction losses",error);
    check(std::abs(mean.lip_flow()-expected)<1e-11,"steady inlet and mouth flow agree",mean.lip_flow()-expected);
    return mean.glottal_flow();
}
}
int main(){
    const double coarse=steady(4,false,false),fine=steady(32,false,false);
    check(std::abs(coarse-fine)<1e-12,"subdividing a uniform tube adds no orifice loss",std::abs(coarse-fine));
    const double constricted=steady(4,true,false),refined=steady(32,true,false);
    check(std::abs(constricted-refined)<1e-12,"subdividing a constriction preserves its total head loss",std::abs(constricted-refined));
    const double reversed=steady(32,true,true);
    check(std::abs(reversed+refined)<1e-11,"symmetric airway reverses its steady flow",reversed+refined);
    const ducklips::AirProperties air{};ducklips::TractShape shape=tube(4);shape.oral_area[1]=0;
    ducklips::AcousticNetwork layout;layout.prepare(4,2,1,{});layout.assemble(shape,air,{});
    ducklips::MeanAirway moving,frozen;moving.prepare(layout,shape,air,1e-8);frozen.prepare(layout,shape,air,1e-8);
    bool valid=true;
    for(std::size_t i=0;i<100;++i){valid=moving.advance(1,2.5e-5,600,1e5,.001).converged&&valid;valid=frozen.advance(1,2.5e-5,600,1e5,.001).converged&&valid;}
    shape.oral_area[3]*=.5;moving.target(shape);double error=0,bridge=0,vented=0;
    for(std::size_t i=0;i<480;++i){
        valid=moving.advance(double(i+1)/480,2.5e-5,600,1e5,1./48000).converged&&valid;
        valid=frozen.advance(1,2.5e-5,600,1e5,1./48000).converged&&valid;
        error=std::max(error,std::abs(moving.pressure()-frozen.pressure()));
        bridge=std::max(bridge,std::abs(moving.oral_flow(1)));vented=std::max(vented,moving.lip_flow());
    }
    check(valid&&error<1e-6&&bridge==0,"spatial front compression cannot cross a sealed stop",error);
    check(vented>0,"front compression vents through the mouth",vented);
    shape=tube(4);shape.oral_area.back()=0;shape.port_area=1e-5;
    layout.assemble(shape,air,{});moving.prepare(layout,shape,air,1e-8);valid=true;
    for(std::size_t i=0;i<1000;++i){valid=moving.advance(1,2.5e-5,600,1e5,.001).converged&&valid;}
    check(valid && moving.lip_flow()==0 && moving.nasal_flow()>0
          && std::abs(moving.nasal_flow()-moving.glottal_flow())<1e-10,
          "open velum routes flow through the nose with a sealed mouth",moving.nasal_flow());
    return failures?1:0;
}
