// Diagnostic harness: current local engine, steady source, passive observation taps.
#include "ducklips/nudges.hpp"
#include "ducklips/body_config.hpp"
#include "ducklips/wav.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

constexpr double dt=1./48000.;
constexpr std::size_t frames=62400;
constexpr double pi=3.14159265358979323846;
void save(const std::filesystem::path& root,const std::string& name,const std::vector<float>& x){
    ducklips::AudioBuffer b;b.samples=x;
    for(float v:x){if(!std::isfinite(v)||std::abs(v)>1){throw std::runtime_error("invalid raw diagnostic audio");}}
    if(!ducklips::write_wav_float((root/(name+"-raw.wav")).string(),b)){throw std::runtime_error("write failed");}
}
int main(int argc,char** argv){
    if(argc!=3){return 2;}
    const std::filesystem::path root=argv[2];std::filesystem::create_directories(root);
    ducklips::VoiceProfile voice=ducklips::baritone_voice();std::string error;
    if(!ducklips::load_body_config(argv[1],voice,error)){throw std::runtime_error(error);}
    voice.jitter=0;voice.shimmer=0;voice.tremor_depth=0;voice.aspiration_gain=0;voice.frication_gain=0;
    ducklips::AirProperties air;
    const double scale=air.density/(4*pi*.4*128);
    // The same inertive glottal inlet, discharging into a constant-pressure
    // reservoir (infinite-compliance limit), with no resonant airway branches.
    ducklips::AcousticNetwork reservoir;reservoir.prepare(4,2,1,{});
    for(ducklips::AcousticNode& n:reservoir.nodes){n.compliance=1e100;n.area=1;n.wall_area=1;}
    for(ducklips::AcousticBranch& b:reservoir.branches){b.conductance=0;}
    ducklips::LocalSpeech inlet;inlet.prepare(reservoir,air,voice.static_wall_compliance,
        voice.folds.lower_thickness+voice.folds.upper_thickness);
    ducklips::TwoMassGlottis folds;folds.reset(voice.folds);
    std::vector<float> free_source(2*frames);double previous=0;
    std::ofstream free_trace(root/"free-source.csv");free_trace<<"time,area,flow,pressure\n";
    for(std::size_t i=0;i<frames;++i){
        folds.advance(voice.lung_pressure-inlet.pressure(),voice.adducted_area,1,dt,2,0,inlet.glottal_flow(),air.density,inlet.pressure());
        inlet.advance(reservoir,reservoir,0,0,dt,1,std::max(folds.minimum_area(),0.),voice.lung_pressure,
            ducklips::ground_node,ducklips::ground_node,ducklips::ground_node,0,0);
        free_source[2*i]=free_source[2*i+1]=float(scale*(inlet.glottal_flow()-previous)/dt);previous=inlet.glottal_flow();
        free_trace<<i*dt<<','<<folds.minimum_area()<<','<<previous<<','<<inlet.pressure()<<'\n';
    }
    save(root,"01-free-folds",free_source);
    for(int condition=0;condition<3;++condition){
        ducklips::SynthesisOptions options;options.local_flow=true;options.conservative_tongue=true;
        options.oral_sections=24;options.nasal_sections=9;options.local_substeps=1;
        options.diagnostic_steady_source=true;options.diagnostic_initial_vowel=condition==0?ducklips::Phone::AX:ducklips::Phone::AA;
        options.trace_start=0;options.trace_samples=frames;options.air_path.full_scale_pressure=128;
        ducklips::NudgeStream stream;stream.prepare(voice,1,4096,options);
        if(condition==2){
            for(const char* cue:{".35 constrict P 1.25 180",".50 release lips 180",".80 constrict P 1.25 180",".95 release lips 180"}){
                if(!stream.submit(cue,error)){throw std::runtime_error(error);}
            }
        }
        std::vector<float> final(2*frames);stream.render(final);
        const auto& stats=stream.engine().statistics();
        if(stats.nonfinite_samples||stats.clipped_samples||stats.geometry_failures||stats.airflow_failures){
            throw std::runtime_error("engine audit failed");
        }
        const auto& trace=stream.engine().trace();if(trace.size()!=frames){throw std::runtime_error("trace truncated");}
        const std::string prefix=condition==0?"neutral":condition==1?"ah":"closures";
        std::vector<float> source(final.size()),lip(final.size());
        std::array<std::vector<float>,3> inside,composite;
        std::array<ducklips::AirPath,3> paths;
        for(std::size_t j=0;j<3;++j){inside[j].resize(final.size());composite[j].resize(final.size());paths[j].prepare(options.air_path,air,{});}
        previous=0;double reconstruction_error=0;
        std::ofstream tr(root/(prefix+"-trace.csv"));tr<<"time,source_flow,source_area,p0,lip_derivative,nose_derivative,wall_derivative\n";
        for(std::size_t i=0;i<frames;++i){
            const auto& s=trace[i];
            source[2*i]=source[2*i+1]=float(scale*(s.glottal_flow-previous)/dt);previous=s.glottal_flow;
            lip[2*i]=lip[2*i+1]=float(scale*s.lips);
            for(std::size_t j=0;j<3;++j){
                inside[j][2*i]=inside[j][2*i+1]=float(scale*s.oral_stage_derivative[j]);
                std::array<double,2> mic{};paths[j].process({s.lips,j>=1?s.nostrils:0,j>=2?voice.tissue_transmission*s.throat:0},mic);
                for(std::size_t ch=0;ch<2;++ch){composite[j][2*i+ch]=float(mic[ch]);}
            }
            reconstruction_error=std::max(reconstruction_error,double(std::abs(composite[2][2*i]-final[2*i])));
            tr<<i*dt<<','<<s.glottal_flow<<','<<s.glottal_area<<','<<s.glottal_pressure<<','<<s.lips<<','<<s.nostrils<<','<<s.throat<<'\n';
        }
        save(root,prefix+"-source",source);save(root,prefix+"-lips",lip);save(root,prefix+"-final",final);
        if(condition==1){
            for(std::size_t j=0;j<3;++j){save(root,"ah-position-"+std::to_string(j+1),inside[j]);
                save(root,"ah-composite-"+std::to_string(j+1),composite[j]);}
        }
        std::cout<<prefix<<": full composite reconstruction error "<<reconstruction_error<<'\n';
        if(reconstruction_error>1e-9){throw std::runtime_error("composite mismatch");}
    }
}
