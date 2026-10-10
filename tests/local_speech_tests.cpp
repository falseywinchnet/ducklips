#include "ducklips/local_speech.hpp"
#include "ducklips/nudges.hpp"
#include "ducklips/tongue_surface.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <vector>

namespace {
std::size_t failures=0,allocations=0;
void check(bool value,const char* text){std::printf("%s %s\n",value?"PASS":"FAIL",text);if(!value){++failures;}}
ducklips::AcousticNetwork pair(double resistance=0){
    ducklips::AcousticNetwork n;n.prepare(4,2,1,{});
    for(ducklips::AcousticNode& node:n.nodes){node.compliance=1;node.area=1;node.wall_area=1;}
    for(ducklips::AcousticBranch& edge:n.branches){edge.conductance=0;}
    n.branches[0]={.from=0,.to=1,.conductance=1,.resistance=resistance};return n;
}
void step(ducklips::LocalSpeech& state,const ducklips::AcousticNetwork& n,double h){
    state.advance(n,n,0,0,h,1,0,0,ducklips::ground_node,ducklips::ground_node,ducklips::ground_node,0,0);
}
double oscillator(double dt){
    ducklips::AirProperties air;air.viscosity=0;
    ducklips::AcousticNetwork n=pair();ducklips::LocalSpeech s;s.prepare(n,air,1);step(s,n,dt);
    s.charge[0]=1;const double initial=s.energy();
    for(std::size_t i=0;i<std::size_t(std::lround(1/dt));++i){step(s,n,dt);}
    check(std::abs(s.energy()-initial)<2e-13,"lossless local exchange conserves energy");
    check(std::abs(s.charge[0]+s.charge[1]-1)<2e-13,"neighbor exchange conserves air");
    return std::abs(s.charge[0]-.5*(1+std::cos(std::sqrt(2.))));
}
double inlet_charge(double dt,double lung){
    ducklips::AirProperties air;air.viscosity=0;
    ducklips::AcousticNetwork n=pair();
    for(ducklips::AcousticBranch& edge:n.branches){edge.conductance=0;}
    for(ducklips::AcousticNode& node:n.nodes){node.compliance=1e-11;}
    ducklips::LocalSpeech s;s.prepare(n,air,1e-8,.003);step(s,n,dt);
    const double area=1e-5;
    double previous=s.energy();bool passive=true,conservative=true;
    // 1/3000 s is an exact integer number of steps at all rates below.
    for(std::size_t i=0;i<std::size_t(std::lround((1./3000)/dt));++i){
        const double before=s.charge[0];
        s.advance(n,n,0,0,dt,1,area,lung,ducklips::ground_node,ducklips::ground_node,ducklips::ground_node,0,0);
        const double air_added=s.charge[0]-before;
        conservative&=std::abs(air_added-dt*s.glottal_flow())<1e-20;
        const double energy=s.energy();
        passive&=energy-previous<=lung*air_added+1e-18;
        previous=energy;
    }
    check(passive,"inertive inlet never creates energy beyond reservoir work");
    check(conservative,"inertive inlet preserves supplied air including reversed breath");
    check(s.glottal_kinetic_energy()>0,"source energy includes moving air between the folds");
    const double pressure=s.pressure();
    s.advance(n,n,0,0,dt,1,0,lung,ducklips::ground_node,ducklips::ground_node,ducklips::ground_node,0,0);
    check(s.glottal_flow()==0 && s.glottal_kinetic_energy()==0 && s.pressure()==pressure,
          "sealed glottis stops inlet flow and retains the cavity air");
    return pressure;
}
}
void* operator new(std::size_t n){++allocations;void* p=std::malloc(std::max<std::size_t>(n,1));if(!p){throw std::bad_alloc();}return p;}
void* operator new[](std::size_t n){return ::operator new(n);}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete[](void* p) noexcept {std::free(p);}
int main(int argc,char** argv){
    if(argc==3){
        const std::size_t count=std::strtoul(argv[1],nullptr,10),substeps=std::strtoul(argv[2],nullptr,10);
        if(count<4 || count>128 || substeps<1 || substeps>16){return 2;}
        ducklips::AirProperties air;air.viscosity=0;
        ducklips::AcousticNetwork n;n.prepare(count,2,1,{});
        const double area=1e-4,length=.175/double(count),c=area*length/(air.density*air.sound_speed*air.sound_speed);
        for(ducklips::AcousticNode& node:n.nodes){node.compliance=c;node.area=area;node.wall_area=1;}
        for(ducklips::AcousticBranch& edge:n.branches){edge.conductance=0;}
        for(std::size_t b=0;b<count-1;++b){n.branches[b]={.from=b,.to=b+1,.conductance=area/(air.density*length)};}
        n.branches[count-1]={.from=count-1,.to=ducklips::ground_node,.kind=ducklips::BranchKind::lip_radiation,.conductance=2*area/(air.density*length)};
        ducklips::LocalSpeech s;s.prepare(n,air,1e-8);step(s,n,1./48000);s.charge[0]=1e-3*c;
        for(std::size_t i=0;i<24000;++i){
            std::printf("%.17g\n",s.pressure());
            s.advance(n,n,0,0,1./48000,substeps,0,0,ducklips::ground_node,ducklips::ground_node,ducklips::ground_node,0,0);
        }
        return 0;
    }
    const double coarse=oscillator(.05),fine=oscillator(.025);
    check(coarse/fine>3.99 && coarse/fine<4.01,"phase error converges quadratically against analytic oscillator");
    const double inlet_coarse=inlet_charge(1./48000,600);
    const double inlet_fine=inlet_charge(1./96000,600);
    const double inlet_reference=inlet_charge(1./384000,600);
    check(std::abs(inlet_coarse-inlet_reference)>3*std::abs(inlet_fine-inlet_reference),
          "inertive glottal exchange converges under timestep refinement");
    check(std::abs(inlet_charge(1./48000,-600)+inlet_coarse)<1e-10,
          "inhale and exhale obey the same pressure-flow law");
    ducklips::AirProperties air;air.viscosity=0;
    ducklips::AcousticNetwork n=pair(.7);ducklips::LocalSpeech s;s.prepare(n,air,1);step(s,n,.1);s.charge[0]=1;
    double previous=s.energy();bool passive=true;
    for(std::size_t i=0;i<100;++i){step(s,n,.1);const double energy=s.energy();passive&=energy<=previous+1e-14;previous=energy;}
    check(passive && previous<.251,"resistive exchange dissipates energy without leaking stored air");
    n.branches[0].conductance=0;const double first=s.charge[0],second=s.charge[1];
    step(s,n,.1);check(s.charge[0]==first && s.charge[1]==second && s.flow(0)==0,"sealed passage carries no flow and retains cavity air");
    // Wall pressure is supported by displacement rather than a DC air leak.
    n=pair();n.branches[0].conductance=0;
    n.branches[n.wall_branch(0)]={.from=0,.to=ducklips::ground_node,.kind=ducklips::BranchKind::wall,.conductance=1,.resistance=1};
    s.prepare(n,air,6);step(s,n,.05);s.charge[0]=1;
    for(std::size_t i=0;i<1000;++i){step(s,n,.05);}
    check(std::abs(s.charge[0]+s.wall_charge[n.wall_branch(0)]-1)<1e-12 && std::abs(s.pressure()-.5)<1e-8,"elastic wall reaches finite static pressure and preserves air");
    ducklips::ConservativeTongueSurface surface;
    check(surface.configure({5.,5.}, {.width=1.}),"prepare polygon body");
    double largest=0;
    for(std::size_t i=0;i<100;++i){
        const double x=3.+4.*double(i)/99;
        surface.set_polygon_queries(false);const double exact=surface.ray_gap({x,10.},{0.,-1.});
        surface.set_polygon_queries(true);const double polygon=surface.ray_gap({x,10.},{0.,-1.});
        if(std::isfinite(exact)&&std::isfinite(polygon)){largest=std::max(largest,std::abs(exact-polygon));}
    }
    check(largest<.002,"carried polygon gap agrees with smooth circle within 20 micrometres");
    ducklips::SynthesisOptions options;options.local_flow=true;options.local_substeps=1;options.oral_sections=24;options.nasal_sections=9;options.conservative_tongue=true;
    ducklips::NudgeStream a,b;a.prepare(ducklips::baritone_voice(),1,4096,options);b.prepare(ducklips::baritone_voice(),1,4096,options);
    std::string error;
    for(const char* cue:{"0 lungs 600 40","0 glottis .000025 170","0 vowel AA 75",".1 glottis .000003 170",".15 diphthong AY .2 75"}){
        check(a.submit(cue,error)&&b.submit(cue,error),"stream accepts articulatory nudge");
    }
    std::vector<float> x(48000),y(x.size());const std::size_t before=allocations;
    a.render(x);
    for(std::size_t i=0;i<24000;i+=17){b.render(std::span<float>(y).subspan(2*i,2*std::min<std::size_t>(17,24000-i)));}
    const std::size_t after=allocations;
    check(x==y,"new engine state is identical across arbitrary host chunks");
    check(before==after,"new engine renders without allocation");
    check(std::all_of(x.begin(),x.end(),[](float v){return std::isfinite(v);}),"streamed speech stays finite");
    check(a.engine().statistics().nonfinite_samples==0 && a.engine().statistics().clipped_samples==0,
          "finite output did not rely on invalid-sample replacement or clipping");
    return failures?1:0;
}
