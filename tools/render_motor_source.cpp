#include "ducklips/tissue_source.hpp"
#include "ducklips/motor_hold.hpp"
#include "ducklips/wav.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc,char** argv){
    if(argc!=2 && argc!=3){return 2;}
    const std::filesystem::path root=argv[1];std::filesystem::create_directories(root);
    const std::uint64_t seed=argc==3?std::stoull(argv[2]):17;
    ducklips::TissueMaterial material;material.rest_gap=.0003;material.relaxation_seed=seed;
    ducklips::TissueSource source;source.reset(material);
    ducklips::MotorHold hold;hold.reset(material.rest_gap,.00025,source.posture_stiffness(),seed+1);
    constexpr double fs=48000,dt=1/fs,compliance=1e-6;
    constexpr std::size_t frames=67200,substeps=2;
    // A finite passive recoil supply: pressure falls with exhaled volume.
    // Initial volume and compliance are prototype body assumptions.
    double remaining_volume=compliance*600,expired=0;
    double previous_force=0,next_force=0;
    std::ofstream setup(root/"setup.json");
    setup<<"{\"renderer\":\"ducklips_motor_source_study\",\"motor_seed\":"<<seed
        <<",\"passive_gap_m\":0.0003,\"sensed_gap_target_m\":0.00025,\"sensor_filter_s\":0.02,\"feedback_delay_s\":0.03,"
        <<"\"fast_twitch_s\":0.015,\"slow_twitch_s\":0.035,\"motor_units\":64,\"interval_cv\":0.2,"
        <<"\"relaxation_fraction\":0.1,\"relaxation_seconds\":0.1,\"relaxation_log_rate_sd\":0.6,"
        <<"\"respiratory_compliance_m3_pa\":1e-6,\"initial_recoil_volume_m3\":0.0006,\"motor_release_s\":1.0,\"release_gap_target_m\":0.0012,\"antagonistic_pools\":true,"
        <<"\"scope\":\"Prototype posture control; no anatomical nerve-loop validation or accepted timbre\"}\n";
    std::ofstream trace(root/"trace.csv");trace<<std::setprecision(17)<<"time,lung,flow,area,kinetic,lower,center,upper,muscle_force,neural_drive,sensed_gap,spikes,expired_volume,muscle_work\n";
    std::ofstream surface(root/"surface.csv");surface<<std::setprecision(10)<<"time,x,z,displacement,body,gap,pressure,contact\n";
    ducklips::AudioBuffer audio;audio.samples.resize(2*frames);
    const double scale=material.air_density/(4*3.14159265358979323846*.4*128);
    double previous=0,simulation_seconds=0;
    constexpr std::size_t center=ducklips::TissueSource::nx/2;
    constexpr std::size_t upper=(ducklips::TissueSource::nz-1)*ducklips::TissueSource::nx+center;
    for(std::size_t frame=0;frame<frames;++frame){
        const double time=frame/fs;double flow=0,lung=0,muscle_force=0;
        const auto started=std::chrono::steady_clock::now();
        if(frame%48==0){previous_force=next_force;next_force=hold.advance(source.gap()[upper],time<1);}
        for(std::size_t step=0;step<substeps;++step){
            const double t=time+step*dt/substeps;
            muscle_force=previous_force+(next_force-previous_force)*(frame%48+double(step)/substeps)/48;
            source.set_holding_force(muscle_force);
            // Initial respiratory hold releases over 12 ms; subsequent
            // pressure is determined by remaining volume, not a fade curve.
            lung=std::max(0.,remaining_volume/compliance-600*std::exp(-t/.012));
            const double flux=source.advance(lung,0,dt/substeps);flow+=flux/substeps;
            const double volume=flux*dt/substeps;remaining_volume=std::max(0.,remaining_volume-volume);expired+=volume;
        }
        simulation_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        audio.samples[2*frame]=audio.samples[2*frame+1]=float(scale*(flow-previous)/dt);previous=flow;
        trace<<time<<','<<lung<<','<<flow<<','<<source.minimum_area()<<','<<source.kinetic_energy()<<','
            <<source.surface()[center]<<','<<source.surface()[(ducklips::TissueSource::nz/2)*ducklips::TissueSource::nx+center]<<','
            <<source.surface()[upper]<<','<<muscle_force<<','<<hold.drive()<<','<<hold.sensed_gap()<<','<<hold.spikes()<<','<<expired<<','<<source.holding_work()<<'\n';
        if(frame%24==0){for(std::size_t z=0;z<ducklips::TissueSource::nz;++z){for(std::size_t x=0;x<ducklips::TissueSource::nx;++x){
            const std::size_t i=z*ducklips::TissueSource::nx+x;
            surface<<time<<','<<x<<','<<z<<','<<source.surface()[i]<<','<<source.body()[i]<<','<<source.gap()[i]<<','<<source.pressure()[i]<<','<<source.contact()[i]<<'\n';
        }}}
    }
    if(!ducklips::write_wav_float((root/"source-raw.wav").string(),audio)){return 1;}
    std::ofstream timing(root/"timing.json");timing<<std::setprecision(17)
        <<"{\"source_computation_seconds\":"<<simulation_seconds<<",\"audio_seconds\":1.4,\"realtime_factor\":"<<1.4/simulation_seconds
        <<",\"scope\":\"source, relaxation, motor controller and finite breath; excludes diagnostic writes\"}\n";
    std::cout<<"simulation_seconds="<<simulation_seconds<<" realtime_factor="<<1.4/simulation_seconds<<'\n';
}
