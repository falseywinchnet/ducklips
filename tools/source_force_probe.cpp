// Built by study_source_forces.py against the frozen source-stage engine.
// The observer reads actual forces at each tissue substep; it changes no state.
#include "ducklips/body_config.hpp"
#include "ducklips/local_speech.hpp"
#include "ducklips/sources.hpp"
#include "ducklips/wav.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace ducklips {
std::ostream* fold_force_trace=nullptr;
std::size_t fold_force_frame=0;
}

int main(int argc,char** argv){
    if(argc!=3){return 2;}
    const std::filesystem::path root=argv[2];
    std::filesystem::create_directories(root);
    ducklips::VoiceProfile voice=ducklips::baritone_voice();std::string error;
    if(!ducklips::load_body_config(argv[1],voice,error)){throw std::runtime_error(error);}
    if(voice.folds.strips!=1||voice.folds.asymmetry!=1){throw std::runtime_error("probe expects symmetric single-strip baseline");}
    ducklips::AirProperties air;
    ducklips::AcousticNetwork reservoir;reservoir.prepare(4,2,1,{});
    for(ducklips::AcousticNode& node:reservoir.nodes){node.compliance=1e100;node.area=1;node.wall_area=1;}
    for(ducklips::AcousticBranch& branch:reservoir.branches){branch.conductance=0;}
    ducklips::LocalSpeech inlet;inlet.prepare(reservoir,air,voice.static_wall_compliance,
        voice.folds.lower_thickness+voice.folds.upper_thickness);
    ducklips::TwoMassGlottis folds;folds.reset(voice.folds);
    std::ofstream forces(root/"forces.csv");forces<<std::setprecision(17);
    forces<<"frame,substep,lower_x,upper_x,lower_v,upper_v,lower_pressure,upper_pressure,lower_elastic,upper_elastic,lower_contact,upper_contact,lower_damping,upper_damping,lower_force,upper_force,lower_gap,upper_gap,lower_mass,upper_mass\n";
    ducklips::fold_force_trace=&forces;
    std::ofstream flow(root/"flow.csv");flow<<std::setprecision(17)<<"time,area,flow,derivative\n";
    constexpr double dt=1./48000.;constexpr std::size_t frames=62400;
    const double scale=air.density/(4*3.14159265358979323846*.4*128);
    ducklips::AudioBuffer audio;audio.samples.resize(2*frames);double previous=0;
    for(std::size_t i=0;i<frames;++i){
        ducklips::fold_force_frame=i;
        folds.advance(voice.lung_pressure-inlet.pressure(),voice.adducted_area,1,dt,2,0,inlet.glottal_flow(),air.density,inlet.pressure());
        inlet.advance(reservoir,reservoir,0,0,dt,1,std::max(folds.minimum_area(),0.),voice.lung_pressure,
            ducklips::ground_node,ducklips::ground_node,ducklips::ground_node,0,0);
        const double derivative=(inlet.glottal_flow()-previous)/dt;previous=inlet.glottal_flow();
        audio.samples[2*i]=audio.samples[2*i+1]=float(scale*derivative);
        flow<<i*dt<<','<<folds.minimum_area()<<','<<previous<<','<<derivative<<'\n';
    }
    ducklips::fold_force_trace=nullptr;
    if(!ducklips::write_wav_float((root/"free-folds-raw.wav").string(),audio)){return 1;}
}
