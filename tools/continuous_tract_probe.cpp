#include "ducklips/graph_continuous_glottis.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
struct Ledger final {
    double work{},loss{},numerical{},net{},port{},minimum_dt{std::numeric_limits<double>::infinity()};
    std::size_t segments{},closures{},reopenings{};
    void add(const Ledger& other){work+=other.work;loss+=other.loss;numerical+=other.numerical;net+=other.net;port+=other.port;
        minimum_dt=std::min(minimum_dt,other.minimum_dt);segments+=other.segments;closures+=other.closures;reopenings+=other.reopenings;}
    void step(const ducklips::GraphContinuousGlottis& model,double dt){const auto& s=model.source.last;
        work+=s.reservoir_work+s.external_work+model.last.graph_external_work;loss+=s.fluid_loss+s.tissue_loss+s.outlet_loss;
        numerical+=s.numerical_loss;net+=dt*(s.flow[0]-model.network.average_flow.back());port+=model.last.graph_port_work;
        minimum_dt=std::min(minimum_dt,dt);++segments;}
};
// Only frozen tract metrics here. Each subdivision repeats the graph solve
// from its own accepted initial state and recomputes its affine load. Moving
// geometry would need new Gauss-stage metrics at the substep times.
bool segment(ducklips::GraphContinuousGlottis& model,std::span<const double> c,std::span<const double> g,
    double dt,double contact_step,ducklips::ContinuousGlottalUpdate update,int remaining,Ledger& ledger){
    auto trial=model;
    if(trial.advance(c,c,g,g,dt,1400.,{},{},{},update)){ledger.step(trial,dt);model=std::move(trial);return true;}
    trial=model;ducklips::ContinuousGlottalState before,after;
    if((contact_step==0 || dt<=contact_step) && model.source.state(before) && trial.advance(c,c,g,g,dt,1400.,{},{},{},ducklips::ContinuousGlottalUpdate::dissipative)
        && trial.source.state(after) && ((before.area[0]>0)!=(after.area[0]>0) || (before.area[1]>0)!=(after.area[1]>0))){
        Ledger accepted;accepted.step(trial,dt);
        for(std::size_t j=0;j<2;++j){accepted.closures+=before.area[j]>0 && after.area[j]<=0?1:0;accepted.reopenings+=before.area[j]<=0 && after.area[j]>0?1:0;}
        ledger.add(accepted);model=std::move(trial);return true;
    }
    if(!remaining){return false;}trial=model;Ledger halves;
    if(!segment(trial,c,g,dt/2,contact_step,update,remaining-1,halves) || !segment(trial,c,g,dt/2,contact_step,update,remaining-1,halves)){return false;}
    ledger.add(halves);model=std::move(trial);return true;
}
// Reference accuracy control, not a realtime path. Compare physical states
// on the same endpoint; accept the two half steps without extrapolation so
// their mass/work receipts remain valid. Scales are explicit absolute state
// scales, matching the independent ODE controls (not a relative-to-zero test).
double difference(const ducklips::GraphContinuousGlottis& a,const ducklips::GraphContinuousGlottis& b){
    ducklips::ContinuousGlottalState x,y;
    if(!a.source.state(x) || !b.source.state(y)){return std::numeric_limits<double>::infinity();}
    double error=0.;
    for(std::size_t j=0;j<2;++j){error=std::max(error,std::abs(x.effort[j]-y.effort[j])/1000.);error=std::max(error,std::abs(x.area[j]-y.area[j])/1e-6);}
    for(std::size_t j=0;j<3;++j){error=std::max(error,std::abs(x.flow[j]-y.flow[j])/1e-4);}
    for(std::size_t j=0;j<2;++j){
        error=std::max(error,std::abs(x.folds.lower_velocity[j]-y.folds.lower_velocity[j])/.03571428571428571);
        error=std::max(error,std::abs(x.folds.upper_velocity[j]-y.folds.upper_velocity[j])/.03571428571428571);
    }
    for(std::size_t j=0;j<a.network.charge.size();++j){error=std::max(error,std::abs(a.network.charge[j]-b.network.charge[j])/1e-8);}
    for(std::size_t j=0;j<a.network.momentum.size();++j){error=std::max(error,std::abs(a.network.momentum[j]-b.network.momentum[j])/.1);}
    return error;
}
bool accurate_segment(ducklips::GraphContinuousGlottis& model,std::span<const double> c,std::span<const double> g,
    double dt,double contact_step,ducklips::ContinuousGlottalUpdate update,double tolerance,int remaining,Ledger& ledger){
    auto full=model,halves=model;Ledger coarse,fine;
    const bool solved=segment(full,c,g,dt,contact_step,update,0,coarse) &&
       segment(halves,c,g,dt/2,contact_step,update,0,fine) &&
       segment(halves,c,g,dt/2,contact_step,update,0,fine);
    // The accuracy controller owns subdivision; a failed trial returns here
    // immediately instead of nesting another recursive integration inside it.
    // Use the raw difference conservatively: no assumption of second-order
    // error cancellation across a nonsmooth contact event.
    if(solved && difference(full,halves)<=tolerance){ledger.add(fine);model=std::move(halves);return true;}
    if(!remaining){return false;}
    auto refined=model;Ledger accepted;
    if(!accurate_segment(refined,c,g,dt/2,contact_step,update,tolerance,remaining-1,accepted) ||
       !accurate_segment(refined,c,g,dt/2,contact_step,update,tolerance,remaining-1,accepted)){return false;}
    ledger.add(accepted);model=std::move(refined);return true;
}
double graph_energy(const ducklips::MovingGraphReference& graph,std::span<const double> c,std::span<const double> g){
    double result=0;for(std::size_t j=0;j<c.size();++j){result+=.5*graph.charge[j]*graph.charge[j]/c[j];}
    for(std::size_t j=0;j<g.size();++j){result+=.5*g[j]*graph.momentum[j]*graph.momentum[j];}return result;
}
double mass(const ducklips::GraphContinuousGlottis& model,const ducklips::ContinuousGlottalState& state){
    double result=state.mass[0]+state.mass[1];for(double charge:model.network.charge){result+=charge;}return result;
}
}
int main(int argc,char** argv){try{
    if(argc<3){throw std::runtime_error("steps duration [--graph-state] [--contact-step seconds] [--relaxation] [--tolerance scaled-error]; stdin: 9 physical source values, optionally 8 charges and 8 momenta");}
    const int steps=std::stoi(argv[1]);const double duration=std::stod(argv[2]);
    bool graph_state=false;double contact_step=0.,tolerance=0.;ducklips::ContinuousGlottalUpdate update=ducklips::ContinuousGlottalUpdate::symmetric;
    for(int argument=3;argument<argc;++argument){
        const std::string option=argv[argument];
        if(option=="--graph-state"){graph_state=true;}
        else if(option=="--tolerance" && argument+1<argc){tolerance=std::stod(argv[++argument]);if(!std::isfinite(tolerance) || tolerance<=0){throw std::runtime_error("tolerance must be positive");}}
        else if(option=="--relaxation"){update=ducklips::ContinuousGlottalUpdate::relaxation;}
        else if(option=="--contact-step" && argument+1<argc){contact_step=std::stod(argv[++argument]);if(!std::isfinite(contact_step) || contact_step<=0){throw std::runtime_error("contact step must be positive");}}
        else{throw std::runtime_error("unknown option");}
    }
    if(steps<=0 || duration<=0 || !std::isfinite(duration)){throw std::runtime_error("invalid duration or steps");}
    std::array<double,9> y;for(double& v:y){if(std::scanf("%lf",&v)!=1){throw std::runtime_error("source seed");}}
    constexpr std::size_t nodes=8;constexpr double length=.175,area=2e-4,dx=length/nodes;
    ducklips::AirProperties air;ducklips::FoldParameters f;ducklips::FoldState initial;initial.rest_area=3e-6;
    initial.lower_position.fill((y[2]-initial.rest_area)/(2*f.length));initial.upper_position.fill((y[3]-initial.rest_area)/(2*f.length));
    initial.lower_velocity.fill(y[7]/(2*f.length));initial.upper_velocity.fill(y[8]/(2*f.length));
    std::vector<ducklips::TreeEdge> edges;for(std::size_t j=0;j<nodes;++j){edges.push_back({j,j+1<nodes?j+1:ducklips::ground_node});}
    std::vector<double> c(nodes,area*dx/(air.density*air.sound_speed*air.sound_speed)),g(nodes,area/(air.density*dx));
    g.back()*=2.; // Last cell center to the ideal ambient-pressure tube end.
    ducklips::GraphContinuousGlottis model;
    model.prepare(nodes,edges,0,f,air,1.,initial,{y[0],y[1]},{y[4],y[5],y[6]},1.,ducklips::ContourPressureProfile::linear_span);
    for(std::size_t j=0;j<nodes;++j){model.network.momentum[j]=y[6]/g[j];}
    if(graph_state){
        for(double& value:model.network.charge){if(std::scanf("%lf",&value)!=1 || !std::isfinite(value)){throw std::runtime_error("initial graph charge");}}
        for(double& value:model.network.momentum){if(std::scanf("%lf",&value)!=1 || !std::isfinite(value)){throw std::runtime_error("initial graph momentum");}}
    }
    ducklips::ContinuousGlottalState state;
    if(!model.source.state(state)){throw std::runtime_error("initial state");}
    const double initial_energy=state.energy+graph_energy(model.network,c,g),initial_mass=mass(model,state),dt=duration/steps;
    Ledger total;double maximum_work=0,maximum_mass=0,minimum_pressure=0,maximum_pressure=0;
    for(int step=0;step<=steps;++step){
        if(!model.source.state(state)){throw std::runtime_error("source state");}
        maximum_work=std::max(maximum_work,std::abs(state.energy+graph_energy(model.network,c,g)-initial_energy+total.loss+total.numerical-total.work));
        maximum_mass=std::max(maximum_mass,std::abs(mass(model,state)-initial_mass-total.net));
        for(std::size_t j=0;j<nodes;++j){const double pressure=model.network.charge[j]/c[j];minimum_pressure=std::min(minimum_pressure,pressure);maximum_pressure=std::max(maximum_pressure,pressure);}
        if(step%20==0 || step==steps || tolerance>0){
            std::printf("{\"step\":%d,\"time\":%.17g,\"source\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g],\"charge\":[",step,step*dt,
                state.effort[0],state.effort[1],state.area[0],state.area[1],state.flow[0],state.flow[1],state.flow[2],
                f.length*(state.folds.lower_velocity[0]+state.folds.lower_velocity[1]),f.length*(state.folds.upper_velocity[0]+state.folds.upper_velocity[1]));
            for(std::size_t j=0;j<nodes;++j){std::printf("%s%.17g",j?",":"",model.network.charge[j]);}std::printf("],\"momentum\":[");
            for(std::size_t j=0;j<nodes;++j){std::printf("%s%.17g",j?",":"",model.network.momentum[j]);}
            std::printf("],\"closures\":%zu,\"reopenings\":%zu,\"port_work\":%.17g,\"numerical_loss\":%.17g}\n",total.closures,total.reopenings,total.port,total.numerical);std::fflush(stdout);
        }
        if(step==steps){
            std::fprintf(stderr,"{\"complete\":true,\"steps\":%d,\"duration\":%.17g,\"nodes\":%zu,\"length_m\":%.17g,\"area_m2\":%.17g,\"segments\":%zu,\"minimum_dt\":%.17g,\"closures\":%zu,\"reopenings\":%zu,\"maximum_work_residual_J\":%.17g,\"maximum_mass_residual_m3\":%.17g,\"numerical_loss_J\":%.17g,\"port_work_J\":%.17g,\"pressure_range_Pa\":[%.17g,%.17g]}\n",steps,duration,nodes,length,area,total.segments,total.minimum_dt,total.closures,total.reopenings,maximum_work,maximum_mass,total.numerical,total.port,minimum_pressure,maximum_pressure);return 0;
        }
        Ledger accepted;const bool advanced=tolerance>0 ? accurate_segment(model,c,g,dt,contact_step,update,tolerance,12,accepted) : segment(model,c,g,dt,contact_step,update,12,accepted);
        if(!advanced){
            std::fprintf(stderr,"{\"complete\":false,\"failed_step\":%d,\"time\":%.17g,\"closures\":%zu,\"reopenings\":%zu}\n",step,step*dt,total.closures,total.reopenings);return 2;
        }total.add(accepted);
    }
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
