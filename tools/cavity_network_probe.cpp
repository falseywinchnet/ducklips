// Two communicating cavities: prescribed lung reservoir -> rear -> front -> ambient.
// This is a pressure reference, not a speech renderer. Input has no header:
// time,Vrear,Vfront,Crear,Cfront,Aglottis,Aneck,Aout,Rglottis,Rneck,Rout,Plung
#include "ducklips/cavity_network.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
int main(int argc,char** argv) {
    const bool lumped=argc==4 && std::string(argv[3])=="--lumped";
    if (argc!=3 && !lumped) { std::fprintf(stderr,"usage: cavity_network_probe initial_rear_pa initial_front_pa [--lumped] < input.csv > output.csv\n");return 2; }
    std::array<double,2> initial{};
    for(std::size_t i=0;i<initial.size();++i){
        char* end=nullptr;initial[i]=std::strtod(argv[i+1],&end);
        if(end==argv[i+1] || *end!='\0' || !std::isfinite(initial[i])){std::fprintf(stderr,"invalid initial pressure\n");return 2;}
    }
    ducklips::CavityNetwork network{};
    const std::array<ducklips::CavityConnection,3> links{{{0,ducklips::cavity_reservoir},{0,1},{1,ducklips::cavity_reservoir}}};
    const std::array<ducklips::CavityConnection,2> single_links{{{0,ducklips::cavity_reservoir},{0,ducklips::cavity_reservoir}}};
    if (lumped) { network.prepare(1,single_links); } else { network.prepare(2,links); }
    const ducklips::AirProperties air{};
    const std::array<double,2> inflow{0,0};std::array<double,2> volumes{},compliances{};
    std::array<ducklips::CavityPassage,3> ports{};double time=0,previous=0;bool ready=false;int scanned=0;
    std::puts("time,rear_pa,front_pa,glottal_inflow,neck_flow,mouth_outflow,mass_residual,iterations");
    while ((scanned=std::scanf("%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&time,
        &volumes[0],&volumes[1],&compliances[0],&compliances[1],&ports[0].area,&ports[1].area,&ports[2].area,
        &ports[0].resistance,&ports[1].resistance,&ports[2].resistance,&ports[0].reservoir_pressure))==12) {
        ducklips::CavityStepReport report{};
        if (!std::isfinite(time)) { return 2; }
        for(std::size_t i=0;i<2;++i){
            if(!std::isfinite(volumes[i]) || volumes[i]<0 || !std::isfinite(compliances[i]) || compliances[i]<=0){return 2;}
        }
        for(const ducklips::CavityPassage& port:ports){
            if(!std::isfinite(port.area) || port.area<0 || !std::isfinite(port.resistance) || port.resistance<0
                || (port.area>0 && port.resistance==0) || !std::isfinite(port.reservoir_pressure)){return 2;}
        }
        std::array<double,1> merged_volume{volumes[0]+volumes[1]},merged_compliance{compliances[0]+compliances[1]},merged_pressure{initial[0]},merged_inflow{0};
        std::array<ducklips::CavityPassage,2> merged_ports{ports[0],{0,ports[1].resistance+ports[2].resistance,0}};
        if (ports[1].area>0 && ports[2].area>0) {
            merged_ports[1].area=1/std::sqrt(1/(ports[1].area*ports[1].area)+1/(ports[2].area*ports[2].area));
        }
        if (!ready) {
            const bool okay=lumped?network.reset(merged_volume,merged_compliance,merged_pressure):network.reset(volumes,compliances,initial);
            if (!okay) { return 2; }
            ready=true;
        } else {
            report=lumped?network.advance(merged_volume,merged_compliance,merged_ports,merged_inflow,time-previous,air)
                :network.advance(volumes,compliances,ports,inflow,time-previous,air);
            if (!report.converged) { std::fprintf(stderr,"pressure solve failed at %.9g, residual %.9g\n",time,report.maximum_mass_residual);return 1; }
        }
        previous=time;
        std::printf("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu\n",time,network.pressures()[0],network.pressures()[lumped?0:1],
            -network.flows()[0],network.flows()[1],network.flows()[lumped?1:2],report.maximum_mass_residual,report.iterations);
    }
    return scanned==EOF && ready?0:2;
}
