// Offline variable-graph mean-pressure reference. Reconfiguration allocates;
// this is not the speech callback. Geometry supplies the continuation map.
// stdin: initial_nodes steps; initial nodes V C p;
// each step: nodes edges dt; nodes: old_index(-1 birth) V C inflow;
// With --partition each node row is: old_index child_V child_C V C inflow.
// edges: from to(-1 reservoir) A R reservoir_pressure forward_loss reverse_loss.
#include "ducklips/cavity_network.hpp"
#include <algorithm>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>
int main(int argc,char** argv){
    const bool partition=argc==2 && std::string(argv[1])=="--partition";
    if(argc!=1 && !partition){std::fprintf(stderr,"usage: cavity_region_motion_probe [--partition] < input.txt\n");return 2;}

    std::size_t nodes=0,steps=0;
    if(std::scanf("%zu %zu",&nodes,&steps)!=2 || !nodes || nodes>65536 || !steps || steps>100000){return 2;}
    std::vector<double> volume(nodes),compliance(nodes),initial(nodes),inflow;
    for(std::size_t i=0;i<nodes;++i){if(std::scanf("%lf %lf %lf",&volume[i],&compliance[i],&initial[i])!=3){return 2;}}
    ducklips::CavityNetwork network;network.prepare(nodes,{});if(!network.reset(volume,compliance,initial)){return 2;}
    std::puts("step,node,pressure");for(std::size_t i=0;i<nodes;++i){std::printf("0,%zu,%.17g\n",i,network.pressures()[i]);}
    const ducklips::AirProperties air{};double residual=0;
    for(std::size_t step=1;step<=steps;++step){
        std::size_t count=0;double dt=0;
        if(std::scanf("%zu %zu %lf",&nodes,&count,&dt)!=3 || !nodes || nodes>65536 || count>262144){return 2;}
        std::vector<ducklips::CavityPartition> children(nodes);
        std::vector<std::size_t> map(nodes);volume.resize(nodes);compliance.resize(nodes);inflow.resize(nodes);
        for(std::size_t i=0;i<nodes;++i){
            long long old=0;
            if(partition){if(std::scanf("%lld %lf %lf %lf %lf %lf",&old,&children[i].volume,&children[i].compliance,&volume[i],&compliance[i],&inflow[i])!=6){return 2;}}
            else if(std::scanf("%lld %lf %lf %lf",&old,&volume[i],&compliance[i],&inflow[i])!=4){return 2;}
            if(old < -1){return 2;}
            map[i]=old==-1?ducklips::cavity_reservoir:static_cast<std::size_t>(old);children[i].old_node=map[i];
        }
        std::vector<ducklips::CavityConnection> links(count);std::vector<ducklips::CavityPassage> passages(count);
        for(std::size_t e=0;e<count;++e){
            long long to=0;auto& p=passages[e];
            if(std::scanf("%zu %lld %lf %lf %lf %lf %lf",&links[e].from,&to,&p.area,&p.resistance,&p.reservoir_pressure,&p.forward_loss,&p.reverse_loss)!=7 || to < -1){return 2;}
            links[e].to=to==-1?ducklips::cavity_reservoir:static_cast<std::size_t>(to);
        }
        try{if(!(partition?network.partition(links,children):network.reconfigure(links,map))){std::fprintf(stderr,"reconfiguration rejected at step %zu\n",step);return 1;}}
        catch(const std::exception& error){std::fprintf(stderr,"%s\n",error.what());return 2;}
        const auto result=network.advance(volume,compliance,passages,inflow,dt,air);
        if(!result.converged){std::fprintf(stderr,"pressure step %zu failed, residual %.17g\n",step,result.maximum_mass_residual);return 1;}
        residual=std::max(residual,result.maximum_mass_residual);
        for(std::size_t i=0;i<nodes;++i){std::printf("%zu,%zu,%.17g\n",step,i,network.pressures()[i]);}
    }
    std::fprintf(stderr,"steps=%zu max_mass_residual=%.17g\n",steps,residual);
}
