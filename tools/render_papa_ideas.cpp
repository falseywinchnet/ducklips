#include "ducklips/synthesizer.hpp"
#include "ducklips/planner.hpp"
#include "ducklips/wav.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char** argv){
    if(argc!=3){std::cerr<<"usage: ducklips_papa_ideas IDEA_0_TO_5 OUTPUT_DIRECTORY\n";return 2;}
    const int idea=std::stoi(argv[1]);if(idea<0||idea>5){return 2;}
    const char* names[]{"current system reference","shared material relaxation","contact recovery damping","energy-based strain hardening","shared axial glottal hydraulics","measured tissue posture feedback"};
    const std::filesystem::path root=argv[2];std::filesystem::create_directories(root);
    auto voice=ducklips::baritone_voice();voice.seed=41;voice.speaking_rate=.8;
    voice.jitter=voice.shimmer=voice.tremor_depth=0;voice.aspiration_gain=voice.frication_gain=0;voice.audible_breath=false;
    ducklips::SynthesisOptions options;options.local_flow=true;options.living_control=true;options.conservative_tongue=true;
    options.tissue_idea=static_cast<ducklips::TissueIdea>(idea);
    options.oral_sections=24;options.nasal_sections=9;options.local_substeps=1;options.air_path.full_scale_pressure=128;
    options.trace_start=0;options.trace_samples=48000*4;
    ducklips::Synthesizer engine;engine.prepare(voice,options);
    ducklips::UtterancePlanner planner;planner.prepare(voice,engine.vowels(),options.rate);
    const auto plan=planner.plan_phonemes("P AA0 P AA1");
    for(const auto& cue:plan.cues){if(!engine.enqueue(cue)){throw std::runtime_error("cue queue full");}}
    const std::size_t frames=std::size_t(std::ceil(plan.duration*48000/48))*48;
    ducklips::AudioBuffer audio;audio.samples.resize(frames*2);
    std::ofstream control(root/"control.csv");control<<std::setprecision(17)<<"time,coordinate,target,position,velocity,stress_memory,muscle,sensed,correction,load\n";
    double elapsed=0;
    for(std::size_t frame=0;frame<frames;frame+=48){
        const auto start=std::chrono::steady_clock::now();engine.render(std::span<float>(audio.samples).subspan(frame*2,96));
        elapsed+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const auto& observation=engine.somatic_control().observations();
        for(std::size_t i=0;i<observation.size();++i){const auto& o=observation[i];control<<double(frame)/48000<<','<<ducklips::SomaticControl::name(i)<<','<<o.target<<','<<o.position<<','<<o.velocity<<','<<o.stress_memory<<','<<o.muscle<<','<<o.sensed<<','<<o.correction<<','<<o.load<<'\n';}
    }
    if(!ducklips::write_wav_float((root/"voice-raw.wav").string(),audio)){return 1;}
    std::ofstream trace(root/"trace.csv");trace<<std::setprecision(17)<<"time,lung,flow,area,lip_derivative,pressure\n";
    for(std::size_t i=0;i<engine.trace().size();++i){const auto& t=engine.trace()[i];trace<<double(i)/48000<<','<<t.lung_pressure<<','<<t.glottal_flow<<','<<t.glottal_area<<','<<t.lips<<','<<t.glottal_pressure<<'\n';}
    std::ofstream timeline(root/"phones.csv");timeline<<"phone,start,end\n";
    for(const auto& phone:plan.phones){timeline<<ducklips::phone_spec(phone.phone).name<<','<<phone.start<<','<<phone.end<<'\n';}
    const auto& s=engine.statistics();
    std::ofstream report(root/"run.json");report<<std::setprecision(17)
        <<"{\"source\":\"elastic-volume TissueSource\",\"control\":\"SomaticControl, 29 coordinates\",\"experiment\":"<<idea<<",\"idea\":\""<<names[idea]<<"\",\"phones\":\"P AA0 P AA1\",\"seed\":41,\"noise\":0.08,\"renderer\":\"ducklips_papa_ideas\","
        <<"\"audio_seconds\":"<<double(frames)/48000<<",\"frames\":"<<frames<<",\"render_seconds\":"<<elapsed<<",\"realtime_factor\":"<<double(frames)/48000/elapsed
        <<",\"somatic_seconds\":"<<s.somatic_seconds<<",\"tissue_seconds\":"<<s.tissue_seconds<<",\"airway_seconds\":"<<s.airway_seconds<<",\"geometry_seconds\":"<<s.geometry_seconds
        <<",\"nonfinite\":"<<s.nonfinite_samples<<",\"clipped\":"<<s.clipped_samples<<",\"airflow_failures\":"<<s.airflow_failures<<",\"geometry_failures\":"<<s.geometry_failures
        <<",\"trace_frames\":"<<engine.trace().size()<<",\"status\":\"fresh Papa experiment; listening unassessed\"}\n";
    std::cout<<idea<<" "<<names[idea]<<": "<<double(frames)/48000<<" seconds, "<<double(frames)/48000/elapsed<<"x realtime; failures="<<s.airflow_failures<<','<<s.geometry_failures<<"\n";
    if(s.geometry_error){std::cerr<<s.geometry_error<<'\n';}
    return s.airflow_failures||s.geometry_failures||s.nonfinite_samples?1:0;
}
