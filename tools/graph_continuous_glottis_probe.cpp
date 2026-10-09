#include "ducklips/graph_continuous_glottis.hpp"
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
void* operator new(std::size_t size){++allocations;if(void* result=std::malloc(size)){return result;}throw std::bad_alloc();}
void operator delete(void* pointer) noexcept {std::free(pointer);}
void operator delete(void* pointer,std::size_t) noexcept {std::free(pointer);}

namespace {
using ducklips::GraphContinuousGlottis;
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
double number(){double value;require(std::scanf("%lf",&value)==1,"fixture number");return value;}
void read(std::vector<double>& values){for(double& v:values){v=number();}}
void print(std::span<const double> values){for(double v:values){std::printf(" %.17g",v);}}
double graph_energy(const ducklips::MovingGraphReference& graph,std::span<const double> c,std::span<const double> g){
    double result=0.;for(std::size_t i=0;i<c.size();++i){result+=.5*graph.charge[i]*graph.charge[i]/c[i];}
    for(std::size_t i=0;i<g.size();++i){result+=.5*g[i]*graph.momentum[i]*graph.momentum[i];}return result;
}
ducklips::FoldState folds(const std::array<double,9>& state){
    ducklips::FoldParameters parameters;ducklips::FoldState result;result.rest_area=3e-6;
    result.lower_position.fill((state[2]-result.rest_area)/(2*parameters.length));
    result.upper_position.fill((state[3]-result.rest_area)/(2*parameters.length));
    result.lower_velocity.fill(state[7]/(2*parameters.length));result.upper_velocity.fill(state[8]/(2*parameters.length));return result;
}
std::array<double,9> physical(const ducklips::ContinuousGlottalState& s){
    return {s.effort[0],s.effort[1],s.area[0],s.area[1],s.flow[0],s.flow[1],s.flow[2],
        .014*(s.folds.lower_velocity[0]+s.folds.lower_velocity[1]),.014*(s.folds.upper_velocity[0]+s.folds.upper_velocity[1])};
}
constexpr std::array<double,9> initial={112.18465038080178,-3.138675468892969,
    3.6868231667933133e-6,4.149026342131548e-6,.00017466261521028343,
    .00017466261521028343,.00017466261521028343,0.,0.};
void prepare(GraphContinuousGlottis& model,std::size_t nodes,std::span<const ducklips::TreeEdge> edges,
    std::size_t port,const std::array<double,9>& state=initial,double loss=1.){
    model.prepare(nodes,edges,port,{}, {},1.,folds(state),{state[0],state[1]},
        {state[4],state[5],state[6]},loss,ducklips::ContourPressureProfile::linear_span);
}
void fixture(){
    unsigned n,m,steps,port;require(std::scanf("%u%u%u%u",&n,&m,&steps,&port)==4,"fixture dimensions");
    require(n>0 && n<=64 && m<=256 && steps<=1000000,"fixture size");
    std::vector<ducklips::TreeEdge> edges(m);
    for(auto& edge:edges){int first,second;require(std::scanf("%d%d",&first,&second)==2 && first>=0 && second>=-1,"fixture edge");edge={static_cast<std::size_t>(first),second<0?ducklips::ground_node:static_cast<std::size_t>(second)};}
    std::array<double,9> y;for(double& v:y){v=number();}const double loss=number();GraphContinuousGlottis model;prepare(model,n,edges,port,y,loss);
    read(model.network.charge);read(model.network.momentum);
    std::vector<double> c1(n),c2(n),g1(m),g2(m),nodal(n),series(m);
    for(unsigned step=0;step<steps;++step){
        const double dt=number(),pressure=number();const int mode=static_cast<int>(number());
        read(c1);read(c2);read(g1);read(g2);read(nodal);read(series);
        if(!model.advance(c1,c2,g1,g2,dt,pressure,{},nodal,series,mode?ducklips::ContinuousGlottalUpdate::dissipative:ducklips::ContinuousGlottalUpdate::symmetric)){
            std::fprintf(stderr,"coupled fixture failed at step %u\n",step);throw std::runtime_error("fixture advance");
        }
        ducklips::ContinuousGlottalState state;require(model.source.state(state),"fixture state");
        const auto y1=physical(state);print(y1);print(model.network.charge);print(model.network.momentum);
        print(model.network.average_pressure);print(model.network.average_flow);
        const auto& s=model.source.last;const auto& r=model.last;
        const double ledger[]={state.energy,state.mass[0],state.mass[1],s.reservoir_work,s.external_work,s.fluid_loss,s.tissue_loss,
            s.outlet_loss,s.numerical_loss,s.energy_residual,r.graph_port_work,r.graph_external_work,r.port_residual,r.flow,r.pressure,r.free_pressure,r.load_slope,
            s.flow[0],s.load_work-s.port_work-s.outlet_loss};print(ledger);std::puts("");
    }
}
void controls(){
    const std::array<ducklips::TreeEdge,4> edges={{{0,1},{1,2},{2,0},{2,ducklips::ground_node}}};
    const std::array<double,4> c={1e-11,1.4e-11,8e-12,2e-11},g={.002,.003,.001,.0025};
    const std::array<double,4> nodal={0.,1e-6,0.,-2e-6},series={0.,30.,0.,-10.};
    GraphContinuousGlottis model;prepare(model,4,edges,1);
    model.network.charge={1e-9,-2e-9,5e-10,3e-9};model.network.momentum={.01,-.02,.015,.005};
    ducklips::ContinuousGlottalState before;require(model.source.state(before),"initial source");
    const double initial_energy=before.energy+graph_energy(model.network,c,g);double work=0,loss=0,worst_energy=0,worst_mass=0,worst_direct=0,worst_port=0;
    double mass=before.mass[0]+before.mass[1];for(double v:model.network.charge){mass+=v;}double transfer=0;
    for(int step=0;step<300;++step){
        ducklips::MovingGraphReference direct=model.network;const double dt=2.5e-6;
        const std::size_t count=allocations.load();
        require(model.advance(c,c,g,g,dt,1400.,{},nodal,series),"coupled fixed advance");
        require(allocations.load()==count,"coupled advance allocated");
        auto drives=nodal;drives[1]+=model.last.flow;require(direct.advance(c,c,g,g,dt,drives,series),"direct graph advance");
        for(std::size_t i=0;i<c.size();++i){worst_direct=std::max(worst_direct,std::abs(direct.charge[i]-model.network.charge[i])/c[i]);worst_direct=std::max(worst_direct,std::abs(direct.average_pressure[i]-model.network.average_pressure[i]));}
        for(std::size_t i=0;i<g.size();++i){require(std::abs(direct.momentum[i]-model.network.momentum[i])<1e-12,"direct momentum");require(std::abs(direct.average_flow[i]-model.network.average_flow[i])<1e-14,"direct flow");}
        const auto& s=model.source.last;work+=s.reservoir_work+s.external_work+model.last.graph_external_work;
        loss+=s.fluid_loss+s.tissue_loss+s.outlet_loss+s.numerical_loss;
        require(s.outlet_loss>=0 && std::abs(s.load_work-s.port_work-s.outlet_loss)<1e-18,"exit loss split");
        ducklips::ContinuousGlottalState state;require(model.source.state(state),"accepted source");
        worst_energy=std::max(worst_energy,std::abs(state.energy+graph_energy(model.network,c,g)-initial_energy+loss-work));
        double current=state.mass[0]+state.mass[1];for(double v:model.network.charge){current+=v;}
        transfer+=dt*(s.flow[0]+nodal[0]+nodal[1]+nodal[2]+nodal[3]-model.network.average_flow[3]);
        worst_mass=std::max(worst_mass,std::abs(current-mass-transfer));worst_port=std::max(worst_port,std::abs(model.last.port_residual));
    }
    require(worst_energy<1e-14 && worst_mass<1e-20 && worst_direct<1e-8 && worst_port<1e-20,"coupled conservation");
    const auto saved=model;ducklips::ContinuousGlottalState saved_source;require(model.source.state(saved_source),"saved source");
    auto bad=c;bad[0]=0.;require(!model.advance(bad,c,g,g,2.5e-6,1400.),"reject graph metric");
    std::array<double,4> force{};force[0]=std::numeric_limits<double>::quiet_NaN();
    require(!model.advance(c,c,g,g,2.5e-6,1400.,force),"reject source after graph trials");
    ducklips::ContinuousGlottalState after;require(model.source.state(after),"rollback state");
    require(physical(after)==physical(saved_source) && after.momentum==saved_source.momentum && after.energy==saved_source.energy
        && model.network.charge==saved.network.charge && model.network.momentum==saved.network.momentum
        && model.network.average_pressure==saved.network.average_pressure && model.network.average_flow==saved.network.average_flow
        && model.last.pressure==saved.last.pressure && model.last.graph_port_work==saved.last.graph_port_work
        && model.source.last.energy_residual==saved.source.last.energy_residual,"rollback all accepted state");
    bool threw=false;try{prepare(model,4,edges,4);}catch(const std::invalid_argument&){threw=true;}
    require(threw && model.network.charge==saved.network.charge,"prepare transaction");
    // Exact sealed outlet: source receives no graph work, graph evolves freely.
    auto sealed=initial;sealed[2]=3e-6;sealed[3]=-1e-7;sealed[0]=sealed[1]=0.;sealed[4]=sealed[5]=sealed[6]=0.;
    prepare(model,4,edges,1,sealed);model.network.charge[1]=1e-9;
    require(model.advance(c,c,g,g,1e-7,0.),"sealed outlet advance");
    require(model.last.flow==0. && model.last.graph_port_work==0. && model.source.last.outlet_loss==0.,"sealed outlet zero transport/work");
    auto reverse=initial;reverse[0]=reverse[1]=0.;reverse[4]=reverse[5]=reverse[6]=-1e-4;
    prepare(model,4,edges,1,reverse);model.network.charge[1]=500*c[1];
    require(model.advance(c,c,g,g,1e-7,0.) && model.last.flow<0. && model.last.graph_port_work<0.
        && model.source.last.outlet_loss>0.,"reverse flow work and dissipation");
    prepare(model,4,edges,1,initial,0.);
    require(model.advance(c,c,g,g,1e-7,1400.) && model.source.last.outlet_loss==0.,"loss-free outlet coupling");
    // The hard roundoff switch formerly made the source residual discontinuous
    // at step 118 of this coupled startup. Do not loosen its force tolerance.
    prepare(model,3,edges,1);
    const std::span<const double> startup_c(c.data(),3);
    const std::array<double,3> startup_drive={0.,1e-6,-2e-6};
    model.network.charge={120*c[0],-80*c[1],40*c[2]};model.network.momentum={.01,-.02,.015,.005};
    for(int step=0;step<140;++step){require(model.advance(startup_c,startup_c,g,g,3.125e-7,1400.,{},startup_drive,series),"roundoff threshold continuity");}
    std::printf("PASS continuous graph: work %.17g J, mass %.17g m3, direct graph %.17g Pa, port %.17g J; rollback and sealed outlet\n",worst_energy,worst_mass,worst_direct,worst_port);
}
}
int main(int argc,char** argv){try{
    if(argc==2 && std::string(argv[1])=="--fixture"){fixture();}else{require(argc==1,"arguments");controls();}return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
