#include "ducklips/lip_viscous.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
void fixture(bool outlet=false){unsigned count;while(std::scanf("%u",&count)==1){
    ducklips::LipProfile p;double first,last,width;require(std::scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&p.length,&p.base_area,&p.activation,&p.place,&p.core,&p.width,&p.target_area,&first,&last,&width)==10,"input");
    require(count<=(outlet?33:32),"order");std::array<ducklips::LipViscousSegment,33> segments;ducklips::LipViscousMeasure measure;
    require(outlet?first==0 && ducklips::lip_outlet_viscous_segments(p,last,width,1.14,std::span(segments).first(count),measure):ducklips::lip_viscous_segments(p,first,last,width,1.14,std::span(segments).first(count),measure),"measure");
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
    std::array<ducklips::LipViscousSegment,17> outlet{};ducklips::LipViscousMeasure measure;
    require(ducklips::lip_outlet_viscous_segments(p,.003,.03,1.14,outlet,measure),"closed outlet");
    require(measure.conductance==0 && outlet[0].inertance==0 && outlet[0].half_gap==0,"core absent exactly");
    for(std::size_t i=1;i<outlet.size();++i){require(outlet[i].inertance>0 && outlet[i].half_gap>0,"shoulder air survives closure");}
    auto open=p;open.activation=.9;
    require(ducklips::lip_outlet_viscous_segments(open,.003,.03,1.14,outlet,measure),"reopened outlet");
    double sum=0;for(const auto& site:outlet){sum+=site.inertance;}
    require(std::abs(sum*measure.conductance-1)<1e-14,"outlet integrated inertia");
    for(double last:{.5*open.core,open.core}){
        require(ducklips::lip_outlet_viscous_segments(open,last,.03,1.14,outlet,measure),"core-only outlet");
        const double expected=1.14*last/open.minimum_area();
        require(std::abs(outlet[0].inertance/expected-1)<1e-14 && std::abs(measure.conductance*expected-1)<1e-14,"exact shorter core inertia");
        for(std::size_t i=1;i<outlet.size();++i){require(outlet[i].inertance==0 && outlet[i].half_gap==0,"zero-length shoulder absent");}
        require(ducklips::lip_outlet_viscous_segments(p,last,.03,1.14,outlet,measure),"closed core-only outlet");
        require(measure.conductance==0,"closed core-only conductance");
        for(const auto& s:outlet){require(s.inertance==0 && s.half_gap==0,"closed core-only air absent");}
    }
    std::printf("PASS lip viscous quadrature: four orders, shared integrated inertia, exact seal and rollback\n");
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else if(argc==2 && std::string(argv[1])=="--outlet-fixture"){fixture(true);}else{require(argc==1,"arguments");controls();}return 0;}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
