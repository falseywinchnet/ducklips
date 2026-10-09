#include "ducklips/moving_graph.hpp"
#include "ducklips/moving_tree.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
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
double energy(const ducklips::MovingGraphReference& s,std::span<const double> c,std::span<const double> g){
    double result=0;for(std::size_t i=0;i<c.size();++i){result+=.5*s.charge[i]*s.charge[i]/c[i];}
    for(std::size_t i=0;i<g.size();++i){result+=.5*g[i]*s.momentum[i]*s.momentum[i];}return result;
}
void tests(){
    ducklips::MovingGraphReference s;
    const ducklips::TreeEdge ring[]={{0,1},{1,2},{2,3},{3,0}};
    const std::vector<double> c(4,1.),g(4,1.);s.prepare(4,ring);s.momentum.assign(4,.3);
    const std::size_t before=allocations.load();
    for(int k=0;k<100;++k){require(s.advance(c,c,g,g,.01),"cycle circulation solve failed");}
    require(allocations.load()==before,"graph stepping allocated");
    for(double q:s.charge){require(q==0,"cycle circulation changed nodal mass");}
    for(double p:s.momentum){require(p==.3,"a cycle path was removed");}
    s.charge={1,-1,1,-1};s.momentum.assign(4,0.);const double e0=energy(s,c,g);
    for(int k=0;k<5000;++k){require(s.advance(c,c,g,g,.001),"ring wave solve failed");}
    require(std::abs(s.charge[0]-std::cos(10.))<1e-10,"ring standing wave mismatch");
    require(std::abs(energy(s,c,g)/e0-1)<1e-12,"closed ring energy drift");
    // A shared outlet joins both mouth regions. Keep the direct mouth route,
    // both outlet entrances and the outlet's own pressure, creating a cycle.
    const ducklips::TreeEdge shared[]={{0,1},{0,2},{1,2},{2,ducklips::ground_node}};
    const double cc[]={1,2,3},gg[]={.2,.3,.4,.5},source[]={.01,0,0},series[]={.1,.2,-.1,.5};
    s.prepare(3,shared);double work=0;
    for(int k=0;k<100;++k){
        require(s.advance(cc,cc,gg,gg,.03,source,series),"shared outlet solve failed");
        for(std::size_t i=0;i<3;++i){work+=.03*source[i]*s.average_pressure[i];}
        for(std::size_t i=0;i<4;++i){work+=.03*series[i]*s.average_flow[i];}
    }
    require(std::abs(energy(s,cc,gg)-work)<1e-13,"graph source ports do not conserve work");
    const auto q=s.charge,p=s.momentum,ap=s.average_pressure,af=s.average_flow;
    const double bad_c[]={1,0,3};
    require(!s.advance(bad_c,cc,gg,gg,.03,source,series),"zero compliance was accepted");
    require(!s.advance(cc,cc,gg,gg,std::numeric_limits<double>::infinity()),"infinite timestep was accepted");
    require(!s.advance(cc,cc,gg,gg,1e308),"unrepresentable stage matrix was accepted");
    require(s.charge==q && s.momentum==p && s.average_pressure==ap && s.average_flow==af,"rejected step changed accepted state");
    const ducklips::TreeEdge isolated[]={{0,1},{2,3}};s.prepare(4,isolated);
    const double isolated_g[]={1,0},input[]={.1,0,0,0};
    for(int k=0;k<20;++k){require(s.advance(c,c,isolated_g,isolated_g,.01,input),"disconnected graph failed");}
    require(s.charge[2]==0 && s.charge[3]==0 && s.average_flow[1]==0,"disconnected graph transmitted sound");
    require(std::abs(s.charge[0]+s.charge[1]-.02)<1e-15,"disconnected component mass mismatch");
    const ducklips::TreeEdge tree_edges[]={{0,1},{1,2},{1,2},{2,ducklips::ground_node}};
    ducklips::MovingTreeReference tree;tree.prepare(3,tree_edges);s.prepare(3,tree_edges);
    s.charge={1e-8,-2e-9,3e-9};tree.charge=s.charge;s.momentum={.02,-.03,.01,.04};tree.momentum=s.momentum;
    double error=0;
    for(int k=0;k<1000;++k){
        const double c1[]={1e-10*(1+.1*std::sin(k*.01)),2e-10,3e-10},c2[]={c1[0]*1.0001,2.001e-10,3.001e-10};
        const double g1[]={.02,.03,.04,.05},g2[]={.0201,.0301,.0401,.0501},drive[]={1e-5,0,0},pressure[]={0,0,.3,0};
        require(s.advance(c1,c2,g1,g2,1./48000,drive,pressure),"moving graph tree limit failed");
        tree.advance(c1,c2,g1,g2,1./48000,drive,pressure);
        for(std::size_t i=0;i<3;++i){error=std::max(error,std::abs(s.charge[i]-tree.charge[i])/c2[i]);}
        for(std::size_t i=0;i<4;++i){require(std::abs(s.momentum[i]-tree.momentum[i])<1e-10,"tree momentum limit mismatch");}
    }
    require(error<1e-7,"tree pressure limit mismatch");
    // Changing volume/flow units together must not change the pressures or
    // momenta. This also checks that a dimensional pivot floor is absent.
    std::vector<double> reference_q,reference_m;
    for(double factor:{1.,1e-6,1e6}){
        s.prepare(3,shared);std::vector<double> scaled_c(3),scaled_g(4),scaled_source(3);
        for(std::size_t i=0;i<3;++i){scaled_c[i]=cc[i]*factor;scaled_source[i]=source[i]*factor;}
        for(std::size_t i=0;i<4;++i){scaled_g[i]=gg[i]*factor;}
        for(int k=0;k<20;++k){require(s.advance(scaled_c,scaled_c,scaled_g,scaled_g,.03,scaled_source,series),"graph scale control failed");}
        if(factor==1){reference_q=s.charge;reference_m=s.momentum;}
        else{
            for(std::size_t i=0;i<3;++i){require(std::abs(s.charge[i]/factor-reference_q[i])<1e-14,"graph pressure changed with common scale");}
            for(std::size_t i=0;i<4;++i){require(std::abs(s.momentum[i]-reference_m[i])<1e-14,"graph momentum changed with common scale");}
        }
    }
    std::printf("PASS cycle circulation, standing wave, energy/source work, shared outlet, isolation, rollback, allocation-free steps; tree pressure error %.17g Pa\n",error);
}
void values(const char* key,std::span<const double> v,bool comma){
    std::printf("%s\"%s\":[",comma?",":"",key);for(std::size_t i=0;i<v.size();++i){std::printf("%s%.17g",i?",":"",v[i]);}std::putchar(']');
}
void read(std::span<double> values){for(double& x:values){require(bool(std::cin>>x),"truncated fixture");}}
void fixture(){
    std::size_t n=0,m=0,steps=0;require(bool(std::cin>>n>>m>>steps) && n>0 && n<=256 && m<=4096 && steps<=100000,"invalid fixture sizes");
    std::vector<ducklips::TreeEdge> edges(m);
    for(auto& edge:edges){long long to=0;require(bool(std::cin>>edge.first_node>>to) && to>=-1,"invalid fixture edge");edge.second_node=to==-1?ducklips::ground_node:static_cast<std::size_t>(to);}
    ducklips::MovingGraphReference s;s.prepare(n,edges);read(s.charge);read(s.momentum);
    std::vector<double> c1(n),c2(n),g1(m),g2(m),nodal(n),series(m);
    for(std::size_t step=0;step<steps;++step){
        double dt=0;require(bool(std::cin>>dt),"missing step period");read(c1);read(c2);read(g1);read(g2);read(nodal);read(series);
        require(s.advance(c1,c2,g1,g2,dt,nodal,series),"graph fixture step failed");
        std::putchar('{');values("charge",s.charge,false);values("momentum",s.momentum,true);values("average_pressure",s.average_pressure,true);values("average_flow",s.average_flow,true);std::puts("}");
    }
}
}
int main(int argc,char** argv){
    try{if(argc==1){tests();}else if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{return 2;}}
    catch(const std::exception& error){std::fprintf(stderr,"%s\n",error.what());return 1;}
}
