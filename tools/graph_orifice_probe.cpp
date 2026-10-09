#include "ducklips/graph_orifice.hpp"
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
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
double energy(const ducklips::GraphOrifice& source,std::span<const double> c,std::span<const double> g,double inertance){
    double total=inertance>0?.5*source.momentum*source.momentum/inertance:0;
    for(std::size_t i=0;i<c.size();++i){total+=.5*source.state.charge[i]*source.state.charge[i]/c[i];}
    for(std::size_t i=0;i<g.size();++i){total+=.5*g[i]*source.state.momentum[i]*source.state.momentum[i];}
    return total;
}
}
int main(int argc,char** argv){try{
    ducklips::GraphOrifice source;const ducklips::AirProperties air{};
    const ducklips::TreeEdge edges[]={{0,1},{1,2},{0,2}};
    const double c[]={2e-10,3e-10,1e-10},g[]={.002,.003,.004};
    const double dt=1./96000,area=3e-6,inertance=25.;
    if(argc==3 && std::string(argv[1])=="--trajectory"){
        const int steps=std::stoi(argv[2]);require(steps>0 && steps<=1000000,"invalid trajectory step count");
        source.prepare(3,edges,2);
        for(int i=0;i<steps;++i){require(source.advance(c,c,g,g,.004/double(steps),area,700,inertance,air),"trajectory source failed");}
        std::printf("{\"charge\":[%.17g,%.17g,%.17g],\"momentum\":[%.17g,%.17g,%.17g],\"source_momentum\":%.17g}\n",
            source.state.charge[0],source.state.charge[1],source.state.charge[2],
            source.state.momentum[0],source.state.momentum[1],source.state.momentum[2],source.momentum);
        return 0;
    }
    require(argc==1,"unknown probe arguments");
    source.prepare(3,edges,2);double work=0,heat=0,injected=0,maximum_balance=0,maximum_law=0;
    const std::size_t before=allocations.load();
    for(int k=0;k<10000;++k){
        const double reservoir=k<5000?700.:-200.;const double old_momentum=source.momentum;
        require(source.advance(c,c,g,g,dt,area,reservoir,inertance,air),"source step failed");
        const double u=source.flow,p=source.load_pressure;
        const double loss=air.density*u*std::abs(u)/(2*area*area);
        maximum_law=std::max(maximum_law,std::abs(reservoir-p-loss-(source.momentum-old_momentum)/dt));
        work+=dt*reservoir*u;heat+=dt*loss*u;injected+=dt*u;
        const double stored=source.state.charge[0]+source.state.charge[1]+source.state.charge[2];
        maximum_balance=std::max(maximum_balance,std::abs(stored-injected));
        require(heat>=0,"orifice generated energy");
    }
    require(allocations.load()==before,"source stepping allocated");
    require(maximum_law<1e-8,"orifice momentum law mismatch");
    require(maximum_balance<1e-18,"total inlet volume not stored");
    require(std::abs(energy(source,c,g,inertance)+heat-work)<1e-13,"reservoir work not balanced by stored energy and orifice heat");
    // At pressure equilibrium a connected closed network accepts no flow.
    source.prepare(3,edges,2);for(std::size_t i=0;i<3;++i){source.state.charge[i]=c[i]*700.;}
    require(source.advance(c,c,g,g,dt,area,700,inertance,air),"equilibrium failed");
    require(std::abs(source.flow)<1e-16,"equilibrium injected flow");
    // The aperture's closure removes only the source passage momentum.
    source.momentum=.001;require(source.advance(c,c,g,g,dt,0,900,inertance,air),"closed aperture failed");
    require(source.flow==0 && source.momentum==0,"closed aperture leaked or retained latent source momentum");
    // Reject malformed geometry/source parameters without changing accepted outputs.
    const auto q=source.state.charge,m=source.state.momentum,ap=source.state.average_pressure,af=source.state.average_flow;
    const double saved_m=source.momentum,saved_p=source.load_pressure,saved_z=source.load_impedance,saved_u=source.flow;
    const double bad[]={2e-10,0,1e-10};
    require(!source.advance(bad,c,g,g,dt,area,700,inertance,air),"zero compliance accepted");
    require(!source.advance(c,c,g,g,dt,-area,700,inertance,air),"negative aperture accepted");
    require(!source.advance(c,c,g,g,dt,area,700,-inertance,air),"negative inlet inertia accepted");
    require(!source.advance(c,c,g,g,dt,area,700,1e308,air),"unrepresentable source impedance accepted");
    require(!source.advance(c,c,g,g,std::numeric_limits<double>::infinity(),area,700,inertance,air),"infinite time accepted");
    require(q==source.state.charge && m==source.state.momentum && ap==source.state.average_pressure && af==source.state.average_flow
        && source.momentum==saved_m && source.flow==saved_u && source.load_pressure==saved_p && source.load_impedance==saved_z,"failed source step changed accepted state");
    // With no inlet inertia the same port has the exact algebraic Bernoulli limit.
    source.prepare(1,{},0);const double lone_c[]={1e-10};
    require(source.advance(lone_c,lone_c,{},{},dt,area,700,0,air),"nodal limit failed");
    require(std::abs(source.load_pressure-.5*dt*source.flow/lone_c[0])<1e-10,"nodal port used endpoint pressure instead of average");
    std::printf("PASS source work/heat, total mass, reverse flow, equilibrium, closure, rollback, allocation-free; mass %.17g m3 law %.17g Pa\n",maximum_balance,maximum_law);
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
