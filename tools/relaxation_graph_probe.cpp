#include "ducklips/relaxation_graph.hpp"
#include "ducklips/viscous_spectrum.hpp"
#include "ducklips/physics.hpp"
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <string>
namespace {std::atomic<std::size_t> allocations{};}
void* operator new(std::size_t size){++allocations;if(void* result=std::malloc(size)){return result;}throw std::bad_alloc();}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
double number(){double x;require(std::scanf("%lf",&x)==1,"number");return x;}
void read(std::vector<double>& values){for(double& x:values){x=number();}}
void print(std::span<const double> values){for(double x:values){std::printf(" %.17g",x);}}
void fixture(){unsigned n,m,a,b,steps;double dt;
    while(std::scanf("%u %u %u %u %lf %u",&n,&m,&a,&b,&dt,&steps)==6){
        require(n>0 && n<100 && m<200 && a<2000 && b<2000 && steps>0 && steps<100000,"dimensions");
        std::vector<ducklips::TreeEdge> edges(m);for(auto& e:edges){e.first_node=static_cast<std::size_t>(number());e.second_node=static_cast<std::size_t>(number());if(e.second_node==n){e.second_node=ducklips::ground_node;}}
        std::vector<std::size_t> na(a),ba(b);for(auto& i:na){i=static_cast<std::size_t>(number());}for(auto& i:ba){i=static_cast<std::size_t>(number());}
        ducklips::RelaxationGraphReference model;model.prepare(n,edges,na,ba);
        std::vector<double> c(n),g(m),r(m),nodal(n),series(m);read(c);read(g);read(r);
        std::vector<ducklips::RelaxationMetric> nm(a),bm(b);for(auto& x:nm){x={number(),number()};}for(auto& x:bm){x={number(),number()};}
        read(model.charge);read(model.momentum);read(model.node_memory);read(model.branch_memory);read(nodal);read(series);
        for(unsigned i=0;i<steps;++i){require(model.advance(c,g,r,nm,bm,dt,nodal,series),"coupled step");}
        print(model.charge);print(model.momentum);print(model.node_memory);print(model.branch_memory);print(model.average_pressure);print(model.average_flow);
        print(model.node_dissipation);print(model.branch_dissipation);print(model.resistance_dissipation);
        const auto& s=model.last;const double receipt[]={s.energy_change,s.physical_dissipation,s.numerical_dissipation,s.external_work,s.energy_residual,s.air_residual};print(receipt);std::puts("");
    }
}
void controls(){
    ducklips::RelaxationGraphReference model;const std::array<ducklips::TreeEdge,4> edges={{{0,1},{1,2},{2,0},{2,ducklips::ground_node}}};
    const std::array<std::size_t,2> no{0,2};const std::array<std::size_t,3> bo{1,1,3};model.prepare(3,edges,no,bo);
    const std::array<double,3> c{1e-11,1.4e-11,8e-12},drive{0,1e-5,0};const std::array<double,4> g{.002,.003,.001,.0025},r{0,1e5,0,1e4};
    const std::array<ducklips::RelaxationMetric,2> nm={{{4e-12,1e5},{3e-12,2e7}}};std::array<ducklips::RelaxationMetric,3> bm={{{12,1e4},{18,3e7},{400,1e5}}};
    model.charge={1e-9,-2e-9,5e-10};model.momentum={.01,-.02,.015,.005};model.node_memory={2e-10,-1e-10};model.branch_memory={.001,-.002,.005};
    double worst=0;for(unsigned i=0;i<500;++i){const auto count=allocations.load();require(model.advance(c,g,r,nm,bm,2.5e-6,drive),"passive coupled step");require(allocations.load()==count,"step allocated");
        require(model.last.physical_dissipation>=0 && model.last.numerical_dissipation>=0 && std::abs(model.last.air_residual)<1e-20,"passivity and air");worst=std::max(worst,std::abs(model.last.energy_residual));}
    require(worst<1e-17,"work identity");const auto saved=model;bm.back().storage=-1;
    require(!model.advance(c,g,r,nm,bm,2.5e-6,drive) && model.charge==saved.charge && model.branch_memory==saved.branch_memory && model.last.energy_residual==saved.last.energy_residual,"transaction rollback");
    bm.back()={0,0};require(!model.advance(c,g,r,nm,bm,2.5e-6,drive),"absent memory cannot silently disappear");
    auto closed=g;closed[3]=0.;model.branch_memory.back()=0;
    require(model.advance(c,closed,r,nm,bm,2.5e-6,drive) && model.average_flow[3]==0 && model.momentum[3]==0,"exact closed branch");
    std::printf("PASS coupled relaxation graph: work residual %.17g J, air, shared modes, closure, rollback and no allocation\n",worst);
}
}
int main(int argc,char** argv){try{if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{require(argc==1,"arguments");controls();}return 0;}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
