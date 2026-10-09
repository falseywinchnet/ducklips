#include "ducklips/graph_phonation.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
namespace {std::atomic<std::size_t> allocations{};}
void* operator new(std::size_t size){++allocations;if(void* p=std::malloc(size)){return p;}throw std::bad_alloc();}
void operator delete(void* p) noexcept{std::free(p);}
void operator delete(void* p,std::size_t) noexcept{std::free(p);}
namespace {
void require(bool result,const char* reason){if(!result){throw std::runtime_error(reason);}}
bool same(const ducklips::FoldState& a,const ducklips::FoldState& b){
    return a.lower_position==b.lower_position && a.lower_velocity==b.lower_velocity && a.upper_position==b.upper_position
        && a.upper_velocity==b.upper_velocity && a.rest_area==b.rest_area;
}
void values(const char* name,std::span<const double> v){
    std::printf("\"%s\":[",name);for(std::size_t i=0;i<v.size();++i){std::printf("%s%.17g",i?",":"",v[i]);}std::putchar(']');
}
void output(const ducklips::GraphPhonation& p){
    std::putchar('{');values("charge",p.source.state.charge);std::putchar(',');values("momentum",p.source.state.momentum);
    const auto s=p.folds.state();const double fs[]={s.lower_position[0],s.lower_position[1],s.lower_velocity[0],s.lower_velocity[1],
        s.upper_position[0],s.upper_position[1],s.upper_velocity[0],s.upper_velocity[1],s.rest_area};
    std::putchar(',');values("fold_state",fs);
    std::printf(",\"source_momentum\":%.17g,\"source_flow\":%.17g,\"midpoint_area\":%.17g,\"transglottal_pressure\":%.17g,\"residual\":%.17g,\"iterations\":%zu}",
        p.source.momentum,p.source.flow,p.midpoint_area,p.transglottal_pressure,p.coupling_residual,p.iterations);
}
}
int main(int argc,char** argv){try{
    ducklips::GraphPhonation p;const ducklips::AirProperties air{};const ducklips::FoldParameters tissue{};
    const ducklips::TreeEdge edges[]={{0,1},{1,2},{0,2},{0,ducklips::ground_node}};
    const double c[]={2e-10,3e-10,1e-10},g[]={.002,.003,.004,.002};
    p.prepare(3,edges,2,tissue);ducklips::FoldState initial;initial.lower_position={1e-6,1e-6};initial.rest_area=3e-6;
    require(p.folds.restore(initial),"initial fold state failed");
    int steps=2400;double duration=.025;
    const bool trajectory=argc==4 && std::string(argv[1])=="--trajectory";
    if(trajectory){steps=std::stoi(argv[2]);duration=std::stod(argv[3]);require(steps>0 && steps<=1000000 && std::isfinite(duration) && duration>0,"invalid trajectory");}
    else{require(argc==1,"unknown arguments");}
    const double dt=duration/double(steps);double inlet=0,outlet=0,mass_error=0,max_residual=0,min_area=1.,max_area=0;
    const std::size_t before=allocations.load();
    for(int i=0;i<steps;++i){
        if(!p.advance(c,c,g,g,dt,3e-6,1.,700.,25.,0,air)){
            std::fprintf(stderr,"phonation failed at step %d of %d (%.17g s)\n",i,steps,i*dt);output(p);return 1;
        }
        inlet+=dt*p.source.flow;outlet+=dt*p.source.state.average_flow.back();
        mass_error=std::max(mass_error,std::abs(p.source.state.charge[0]+p.source.state.charge[1]+p.source.state.charge[2]+outlet-inlet));
        max_residual=std::max(max_residual,std::abs(p.coupling_residual));min_area=std::min(min_area,p.midpoint_area);max_area=std::max(max_area,p.midpoint_area);
    }
    require(allocations.load()==before,"phonation step allocated");require(mass_error<1e-17,"phonation lost stored air");
    if(trajectory){std::printf("{\"final\":");output(p);std::printf(",\"mass_error\":%.17g,\"maximum_coupling_residual\":%.17g,\"aperture_range\":[%.17g,%.17g]}\n",mass_error,max_residual,min_area,max_area);return 0;}
    const auto q=p.source.state.charge,m=p.source.state.momentum,ap=p.source.state.average_pressure,af=p.source.state.average_flow;
    const auto folds=p.folds.state();const double sm=p.source.momentum,flow=p.source.flow,load=p.source.load_pressure,z=p.source.load_impedance;
    const double midpoint=p.midpoint_area,drop=p.transglottal_pressure,error=p.coupling_residual;const auto iterations=p.iterations;
    require(!p.advance(c,c,g,g,dt,3e-6,0,700,25,0,air),"zero tension accepted");
    const double bad[]={2e-10,0,1e-10};require(!p.advance(bad,c,g,g,dt,3e-6,1,700,25,0,air),"zero compliance accepted");
    require(!p.advance(c,c,g,g,dt,3e-6,1,700,1e308,0,air),"unrepresentable inlet accepted");
    require(q==p.source.state.charge && m==p.source.state.momentum && ap==p.source.state.average_pressure && af==p.source.state.average_flow
        && same(folds,p.folds.state()) && p.source.momentum==sm && p.source.flow==flow && p.source.load_pressure==load && p.source.load_impedance==z
        && p.midpoint_area==midpoint && p.transglottal_pressure==drop && p.coupling_residual==error && p.iterations==iterations,"failed coupling changed accepted state");
    ducklips::FoldState invalid=folds;invalid.lower_velocity[0]=std::numeric_limits<double>::quiet_NaN();
    require(!p.folds.restore(invalid) && same(folds,p.folds.state()),"nonfinite fold restoration changed state");
    // A fresh coupling instance continues every physical state exactly.
    ducklips::GraphPhonation continued;continued.prepare(3,edges,2,tissue);continued.source.state.charge=q;continued.source.state.momentum=m;
    continued.source.momentum=sm;require(continued.folds.restore(folds),"fold continuation failed");
    require(p.advance(c,c,g,g,dt,3e-6,1,700,25,0,air) && continued.advance(c,c,g,g,dt,3e-6,1,700,25,0,air),"continued advance failed");
    require(same(p.folds.state(),continued.folds.state()) && p.source.state.charge==continued.source.state.charge
        && p.source.state.momentum==continued.source.state.momentum && p.source.momentum==continued.source.momentum,"continuation changed the trajectory");
    // A genuinely sealed glottis removes inlet momentum, keeps downstream
    // air, and applies the static pressure drop rather than the stop impulse
    // to the already closed fold surfaces.
    ducklips::FoldParameters hard=tissue;hard.contact_softness=0;
    p.prepare(1,{},0,hard);const double sealed_c[]={1e-10};p.source.state.charge[0]=700*sealed_c[0];p.source.momentum=.001;
    require(p.advance(sealed_c,sealed_c,{},{},dt,-1e-6,1,900,25,0,air),"sealed folds failed");
    require(p.midpoint_area==0 && p.source.flow==0 && p.source.momentum==0,"sealed folds leaked or kept inlet momentum");
    require(p.source.state.charge[0]==700*sealed_c[0] && std::abs(p.transglottal_pressure-200)<1e-7,"sealed folds lost air or applied an inertial stop impulse as steady pressure");
    std::printf("PASS coupled moving folds, pressure residual, stored air, rollback, exact continuation, no allocation; mass %.17g m3 residual %.17g Pa aperture [%.17g, %.17g] m2\n",mass_error,max_residual,min_area,max_area);
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
