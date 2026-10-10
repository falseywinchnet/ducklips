#include "ducklips/lip_viscous.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
void fixture(){unsigned count;while(std::scanf("%u",&count)==1){
    ducklips::LipProfile p;double first,last,width;require(std::scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&p.length,&p.base_area,&p.activation,&p.place,&p.core,&p.width,&p.target_area,&first,&last,&width)==10,"input");
    require(count<=32,"order");std::array<ducklips::LipViscousSegment,32> segments;ducklips::LipViscousMeasure measure;
    require(ducklips::lip_viscous_segments(p,first,last,width,1.14,std::span(segments).first(count),measure),"measure");
    std::printf("%.17g %.17g %.17g",measure.conductance,measure.quadrature_inertance,measure.relative_inertance_correction);
    for(std::size_t i=0;i<count;++i){const auto& s=segments[i];std::printf(" %.17g %.17g %.17g",s.position,s.inertance,s.half_gap);}std::puts("");
}}
void controls(){
    const ducklips::LipProfile p{.02,1e-4,1.,0.,.001,.004,0.};
    for(std::size_t n:{4,8,16,32}){
        std::array<ducklips::LipViscousSegment,32> s{};ducklips::LipViscousMeasure m;auto active=std::span(s).first(n);
        require(ducklips::lip_viscous_segments(p,.002,.008,.03,1.14,active,m),"open narrowing");
        double sum=0;for(const auto& a:active){sum+=a.inertance;require(a.position>.002 && a.position<.008 && a.half_gap>0,"positive open samples");}
        require(std::abs(sum*m.conductance-1)<1e-14,"same integrated inertia");
        const auto saved=s;const auto before=m;
        require(!ducklips::lip_viscous_segments(p,0,.008,.03,1.14,active,m),"exact seal rejected");
        require(m.conductance==before.conductance && s[0].inertance==saved[0].inertance,"closure rollback");
        require(!ducklips::lip_viscous_segments(p,.002,.008,0.,1.14,active,m),"invalid width");
    }
    std::printf("PASS lip viscous quadrature: four orders, shared integrated inertia, exact seal and rollback\n");
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{require(argc==1,"arguments");controls();}return 0;}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
