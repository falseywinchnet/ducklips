#include "ducklips/graph_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {
struct CompressionPeak {double ratio{},time{},mass{},volume{},volume_rate{},mass_rate{},pressure{},resistance_in{},resistance_out{},relaxation{};std::array<double,3> flow{};std::array<double,2> areas{};};
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}}
int main(int argc,char** argv){try{
    require(argc==3,"usage: ducklips_glottal_geometry_motion_probe STEPS CSV");const int steps=std::stoi(argv[1]);require(steps>=1000 && steps<=1000000,"step count");
    const ducklips::TreeEdge edges[]={{0,1},{1,2},{0,2},{0,ducklips::ground_node}};
    const double c[]={2e-10,3e-10,1e-10},g[]={.002,.003,.004,.002};const ducklips::AirProperties air;const ducklips::FoldParameters tissue;
    ducklips::GraphGlottalStorage p;ducklips::FoldState initial;initial.rest_area=3e-6;p.prepare_articulated(3,edges,2,tissue,air,1.,initial,{});
    const double h[]={tissue.lower_thickness,tissue.upper_thickness};const double dt=.1/double(steps);
    std::FILE* csv=std::fopen(argv[2],"w");require(csv,"open motion CSV");
    std::fprintf(csv,"time_s,lower_area_m2,upper_area_m2,inlet_m3_s,between_m3_s,outlet_m3_s,lower_pressure_Pa,upper_pressure_Pa,load_effort_Pa\n");
    double min_area=1.,max_mach=0,max_reynolds=0,min_density=1.,max_density=1.,max_work=0;
    ducklips::GlottalKinetics kinetics;require(kinetics.configure(tissue,air,1.),"kinetic diagnostics");std::array<CompressionPeak,2> peaks{};
    std::array<double,2> tail_min={1.,1.},tail_max{};std::size_t samples=0;
    for(int i=0;i<steps;++i){const auto m0=p.chambers.state.mass;std::array<double,2> v0{},v1{};require(p.chambers.volumes(v0),"initial volume");
        require(p.advance(c,c,g,g,dt,700,3e-6),"motion advance");require(p.chambers.volumes(v1),"final volume");
        for(std::size_t j=0;j<2;++j){const double mid_mass=.5*(m0[j]+p.chambers.state.mass[j]),mid_volume=.5*(v0[j]+v1[j]);
            const double area=mid_volume/h[j],r=mid_mass/mid_volume;
            if(r>peaks[j].ratio){
                std::array<double,2> mm{},vv{};for(std::size_t cell=0;cell<2;++cell){mm[cell]=.5*(m0[cell]+p.chambers.state.mass[cell]);vv[cell]=.5*(v0[cell]+v1[cell]);}
                std::array<double,3> resistance{},quadratic{};require(kinetics.losses(mm,vv,resistance,quadratic),"peak losses");
                const double rin=resistance[j]+2*quadratic[j]*std::abs(p.last.flow[j]),rout=resistance[j+1]+2*quadratic[j+1]*std::abs(p.last.flow[j+1]);
                const double compliance=mid_volume/(air.density*air.sound_speed*air.sound_speed*std::pow(r,air.heat_capacity_ratio-2));
                peaks[j]={r,(i+.5)*dt,mid_mass,mid_volume,(v1[j]-v0[j])/dt,p.last.flow[j]-p.last.flow[j+1],p.chambers.last.average_pressure[j],rin,rout,compliance/(1/rin+1/rout),p.last.flow,{vv[0]/h[0],vv[1]/h[1]}};
            }
            min_area=std::min(min_area,v1[j]/h[j]);min_density=std::min(min_density,r);max_density=std::max(max_density,r);
            for(std::size_t port=j;port<=j+1;++port){const double speed=std::abs(p.last.flow[port])*h[j]/mid_mass;
                max_mach=std::max(max_mach,speed/air.sound_speed);
                const double gap=area/tissue.length,hydraulic_diameter=2*gap*tissue.length/(gap+tissue.length);
                max_reynolds=std::max(max_reynolds,air.density*r*speed*hydraulic_diameter/air.viscosity);
            }
            if((i+1)*dt>=.08){tail_min[j]=std::min(tail_min[j],v1[j]/h[j]);tail_max[j]=std::max(tail_max[j],v1[j]/h[j]);}
        }
        max_work=std::max(max_work,std::abs(p.last.energy_residual));
        if(i%std::max(1,steps/1000)==0 || i+1==steps){++samples;
            std::fprintf(csv,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",(i+1)*dt,v1[0]/h[0],v1[1]/h[1],p.last.flow[0],p.last.flow[1],p.last.flow[2],p.chambers.last.average_pressure[0],p.chambers.last.average_pressure[1],p.last.load_effort);
        }
    }
    require(std::fclose(csv)==0,"close motion CSV");
    std::printf("{\"steps\":%d,\"duration_s\":0.1,\"samples\":%zu,\"minimum_area_m2\":%.17g,\"maximum_mach\":%.17g,\"maximum_hydraulic_reynolds\":%.17g,\"density_ratio_range\":[%.17g,%.17g],\"maximum_step_work_residual_J\":%.17g,\"last_20ms_area_min\":[%.17g,%.17g],\"last_20ms_area_max\":[%.17g,%.17g]",steps,samples,min_area,max_mach,max_reynolds,min_density,max_density,max_work,tail_min[0],tail_min[1],tail_max[0],tail_max[1]);
    std::printf(",\"compression_peaks\":[");
    for(std::size_t j=0;j<2;++j){const auto& e=peaks[j];std::printf("%s{\"cell\":%zu,\"density_ratio\":%.17g,\"time_s\":%.17g,\"mass_m3\":%.17g,\"volume_m3\":%.17g,\"volume_rate_m3_s\":%.17g,\"mass_rate_m3_s\":%.17g,\"pressure_Pa\":%.17g,\"inlet_differential_resistance\":%.17g,\"outlet_differential_resistance\":%.17g,\"local_viscous_relaxation_s\":%.17g,\"flow\":[%.17g,%.17g,%.17g],\"areas_m2\":[%.17g,%.17g]}",j?",":"",j,e.ratio,e.time,e.mass,e.volume,e.volume_rate,e.mass_rate,e.pressure,e.resistance_in,e.resistance_out,e.relaxation,e.flow[0],e.flow[1],e.flow[2],e.areas[0],e.areas[1]);}
    std::puts("]}");
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
