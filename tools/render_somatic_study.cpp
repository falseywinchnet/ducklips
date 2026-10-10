#include "ducklips/synthesizer.hpp"
#include "ducklips/wav.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string_view>

int main(int argc,char** argv){
    if(argc!=2 && !(argc==3 && std::string_view(argv[2])=="--no-control-variation")){return 2;}
    const std::filesystem::path root=argv[1];std::filesystem::create_directories(root);
    auto voice=ducklips::baritone_voice();voice.seed=41;voice.jitter=voice.shimmer=voice.tremor_depth=0;
    voice.aspiration_gain=voice.frication_gain=0;voice.audible_breath=false;
    ducklips::SynthesisOptions options;options.local_flow=true;options.living_control=true;options.conservative_tongue=true;
    if(argc==3){options.somatic.noise=0;}
    options.oral_sections=24;options.nasal_sections=9;options.local_substeps=1;options.air_path.full_scale_pressure=128;
    options.trace_start=0;options.trace_samples=57600;
    ducklips::Synthesizer engine;engine.prepare(voice,options);
    auto send=[&engine](double t,ducklips::CueKind kind,double value,double rate){
        ducklips::Cue cue;cue.time=std::llround(t*48000);cue.kind=kind;cue.value=value;cue.angular_rate=rate;
        if(!engine.enqueue(cue)){throw std::runtime_error("cue queue full");}
    };
    ducklips::Cue vowel;vowel.kind=ducklips::CueKind::vowel;vowel.vowel=engine.vowels().shape(ducklips::Phone::AA);vowel.angular_rate=60;
    if(!engine.enqueue(vowel)){return 1;}
    send(0,ducklips::CueKind::glottis,voice.adducted_area,60);
    send(0,ducklips::CueKind::lungs,600,60);
    send(.75,ducklips::CueKind::lungs,0,60);
    send(.75,ducklips::CueKind::glottis,voice.abducted_area,60);
    ducklips::AudioBuffer audio;audio.samples.resize(57600*2);
    std::ofstream control(root/"control.csv");control<<std::setprecision(17)<<"time,coordinate,target,position,velocity,stress_memory,muscle,sensed,correction,load\n";
    double elapsed=0;
    for(std::size_t frame=0;frame<57600;frame+=48){
        const auto started=std::chrono::steady_clock::now();
        engine.render(std::span<float>(audio.samples).subspan(frame*2,96));
        elapsed+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        const auto& observation=engine.somatic_control().observations();
        for(std::size_t i=0;i<observation.size();++i){const auto& o=observation[i];
            control<<double(frame)/48000<<','<<ducklips::SomaticControl::name(i)<<','<<o.target<<','<<o.position<<','<<o.velocity<<','<<o.stress_memory<<','<<o.muscle<<','<<o.sensed<<','<<o.correction<<','<<o.load<<'\n';}
    }
    if(!ducklips::write_wav_float((root/"voice-raw.wav").string(),audio)){return 1;}
    std::ofstream trace(root/"trace.csv");trace<<std::setprecision(17)<<"time,lung,flow,area,lip_derivative,pressure\n";
    for(std::size_t i=0;i<engine.trace().size();++i){const auto& t=engine.trace()[i];trace<<double(i)/48000<<','<<t.lung_pressure<<','<<t.glottal_flow<<','<<t.glottal_area<<','<<t.lips<<','<<t.glottal_pressure<<'\n';}
    const auto& statistics=engine.statistics();
    std::ofstream report(root/"run.json");report<<std::setprecision(17)
        <<"{\"source\":\"elastic-volume TissueSource\",\"control\":\"SomaticControl, 29 coordinates\",\"seed\":41,\"noise\":"<<options.somatic.noise<<",\"opposing_tone\":0.02,"
        <<"\"audio_seconds\":1.2,\"render_seconds\":"<<elapsed<<",\"realtime_factor\":"<<1.2/elapsed
        <<",\"somatic_seconds\":"<<statistics.somatic_seconds<<",\"tissue_seconds\":"<<statistics.tissue_seconds<<",\"airway_seconds\":"<<statistics.airway_seconds
        <<",\"geometry_seconds\":"<<statistics.geometry_seconds<<",\"nonfinite\":"<<statistics.nonfinite_samples<<",\"clipped\":"<<statistics.clipped_samples
        <<",\"airflow_failures\":"<<statistics.airflow_failures<<",\"geometry_failures\":"<<statistics.geometry_failures
        <<",\"trace_frames\":"<<engine.trace().size()<<",\"remaining_breath_m3\":"<<engine.somatic_control().remaining_breath()
        <<",\"status\":\"fresh integrated experiment; perceptual quality unaccepted\"}\n";
    std::cout<<"render_seconds="<<elapsed<<" realtime_factor="<<1.2/elapsed<<" geometry="<<statistics.geometry_seconds
        <<" failures="<<statistics.airflow_failures<<','<<statistics.geometry_failures<<" nonfinite="<<statistics.nonfinite_samples<<'\n';
    return statistics.airflow_failures||statistics.geometry_failures||statistics.nonfinite_samples?1:0;
}
