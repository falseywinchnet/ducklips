#include "ducklips/tissue_source.hpp"
#include "ducklips/wav.hpp"
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc,char** argv){
    if(argc!=2 && argc!=3){return 2;}
    const std::filesystem::path root=argv[1];std::filesystem::create_directories(root);
    ducklips::TissueMaterial material;
    const bool apposed=argc==3 && std::string(argv[2])=="--apposed";
    if(apposed){material.rest_gap=0;}
    ducklips::TissueSource source;source.reset(material);
    constexpr double fs=48000,dt=1/fs;
    constexpr std::size_t frames=48000,substeps=2;
    const bool adducted=argc==3 && std::string(argv[2])=="--adducted";
    if(argc==3 && !adducted && !apposed){return 2;}
    if(adducted){
        source.set_adduction(-.0006);
        // Initialize the prephonatory posture under zero breath. This is
        // tissue equilibration before recording, not an audio fade or cut.
        for(std::size_t i=0;i<28800;++i){source.advance(0,0,dt/substeps);}
    }
    std::ofstream setup(root/"setup.json");
    setup<<"{\"adducted\":"<<(adducted?"true":"false")<<",\"apposed\":"<<(apposed?"true":"false")
        <<",\"central_rest_gap_m\":"<<(adducted?-.0006:material.rest_gap)
        <<",\"relaxation_fraction\":"<<material.relaxation_fraction
        <<",\"relaxation_seconds\":"<<material.relaxation_seconds
        <<",\"relaxation_log_rate_sd\":"<<material.relaxation_log_rate_sd
        <<",\"relaxation_seed\":"<<material.relaxation_seed
        <<",\"initialization_seconds\":"<<(adducted?.3:0)
        <<",\"initial_kinetic_j\":"<<source.kinetic_energy()
        <<",\"cover_shear_pa\":"<<material.cover_shear<<",\"cover_viscosity_pa_s\":"<<material.viscosity<<"}\n";
    std::ofstream trace(root/"trace.csv");trace<<std::setprecision(17)<<"time,lung,flow,area,kinetic,lower,center,upper\n";
    std::ofstream surface(root/"surface.csv");surface<<std::setprecision(10)<<"time,x,z,displacement,body,gap,pressure,contact\n";
    ducklips::AudioBuffer audio;audio.samples.resize(2*frames);
    const double scale=material.air_density/(4*3.14159265358979323846*.4*128);
    double previous=0,simulation_seconds=0;
    for(std::size_t frame=0;frame<frames;++frame){
        const double time=frame/fs;
        double flow=0,lung=0;
        const auto started=std::chrono::steady_clock::now();
        for(std::size_t step=0;step<substeps;++step){
            const double t=time+step*dt/substeps;
            // Breath pressure rises from rest, stays steady, then releases.
            // This controls the breath, never the vibration within a cycle.
            lung=600*(1-std::exp(-t/.012));
            if(t>.8){lung*=std::exp(-(t-.8)/.012);}
            flow+=source.advance(lung,0,dt/substeps)/substeps;
        }
        simulation_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        audio.samples[2*frame]=audio.samples[2*frame+1]=float(scale*(flow-previous)/dt);previous=flow;
        constexpr std::size_t center=ducklips::TissueSource::nx/2;
        trace<<time<<','<<lung<<','<<flow<<','<<source.minimum_area()<<','<<source.kinetic_energy()<<','
            <<source.surface()[center]<<','<<source.surface()[(ducklips::TissueSource::nz/2)*ducklips::TissueSource::nx+center]<<','
            <<source.surface()[(ducklips::TissueSource::nz-1)*ducklips::TissueSource::nx+center]<<'\n';
        if(frame%24==0){
            for(std::size_t z=0;z<ducklips::TissueSource::nz;++z){for(std::size_t x=0;x<ducklips::TissueSource::nx;++x){
                const std::size_t i=z*ducklips::TissueSource::nx+x;
                surface<<time<<','<<x<<','<<z<<','<<source.surface()[i]<<','<<source.body()[i]<<','<<source.gap()[i]<<','<<source.pressure()[i]<<','<<source.contact()[i]<<'\n';
            }}
        }
    }
    if(!ducklips::write_wav_float((root/"source-raw.wav").string(),audio)){return 1;}
    std::ofstream timing(root/"timing.json");
    timing<<std::setprecision(17)<<"{\"source_computation_seconds\":"<<simulation_seconds
        <<",\"audio_seconds\":1,\"realtime_factor\":"<<1/simulation_seconds
        <<",\"scope\":\"source only; excludes diagnostic file writes and posture initialization\"}\n";
    std::cout<<std::setprecision(10)<<"simulation_seconds="<<simulation_seconds<<" realtime_factor="<<1/simulation_seconds<<'\n';
}
