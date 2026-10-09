// Frozen-geometry speed and impulse probe. This is not a speech renderer.
#include "frozen_modes.hpp"
#include "ducklips/sources.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc,char** argv) {
    try {
        if (argc!=4 && argc!=6) { throw std::invalid_argument("usage: surface_modal_probe model.modes seconds voice|ports|impulse [source out.f64]"); }
        ducklips::FrozenModes model;
        model.load(argv[1]);
        std::size_t used=0;
        const std::string duration_text=argv[2];
        const double seconds=std::stod(duration_text,&used);
        if (used!=duration_text.size() || !std::isfinite(seconds) || seconds<=0 || seconds>60) { throw std::invalid_argument("duration must be (0,60]"); }
        const bool voice=std::string(argv[3])=="voice";
        const bool ports=std::string(argv[3])=="ports";
        const bool impulse=std::string(argv[3])=="impulse";
        if ((!voice && !ports && !impulse) || (impulse ? argc!=6 : argc!=4)) { throw std::invalid_argument("invalid mode/arguments"); }
        std::ofstream file;
        std::size_t source=0;
        if (impulse) {
            const std::string source_text=argv[4];source=std::stoul(source_text,&used);
            if (used!=source_text.size() || source>=model.inputs()) { throw std::invalid_argument("invalid source"); }
            file.open(argv[5],std::ios::binary);
            if (!file) { throw std::runtime_error("cannot open impulse output"); }
        }
        const std::size_t samples=static_cast<std::size_t>(seconds*model.sample_rate());
        std::vector<double> input(model.inputs()),output(model.outputs());
        // Allocate impulse storage outside the measured loop, then write once.
        std::vector<double> captured(impulse ? samples*model.outputs() : 0);
        ducklips::TwoMassGlottis glottis;
        ducklips::FoldParameters folds;folds.strips=8;
        glottis.reset(folds);
        ducklips::NoiseGenerator noise;noise.seed(173);
        ducklips::AirProperties air;
        double previous_pressure=0,previous_radiation=0,energy=0,peak=0;
        const auto start=std::chrono::steady_clock::now();
        for (std::size_t sample=0;sample<samples;++sample) {
            double microphone=0;
            if (voice) {
                const double lung=800.0;
                glottis.advance(lung-previous_pressure,3e-6,1.0,1/model.sample_rate(),1,0.1*noise.gaussian());
                const double free=model.free_pressure(0);
                const double flow=ducklips::orifice_flow(glottis.minimum_area(),lung-free,model.direct(0,0),air);
                previous_pressure=free+model.direct(0,0)*flow;
                const double radiation=model.tick_glottis(flow);
                microphone=model.sample_rate()*(radiation-previous_radiation);previous_radiation=radiation;
            } else if (ports) {
                // Numerical throughput exercise: all five ports active, no
                // audio file or claim that these arbitrary drives are speech.
                for (std::size_t port=0;port<input.size();++port) { input[port]=1e-6*noise.gaussian(); }
                model.tick(input,output);
                microphone=model.sample_rate()*(output.back()-previous_radiation);previous_radiation=output.back();
            } else {
                input[source]=sample==0 ? 1e-10 : 0.0;
                model.tick(input,output);
                for (std::size_t channel=0;channel<model.outputs();++channel) { captured[sample*model.outputs()+channel]=output[channel]; }
                microphone=model.sample_rate()*(output.back()-previous_radiation);previous_radiation=output.back();
            }
            energy+=microphone*microphone;peak=std::max(peak,std::abs(microphone));
        }
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if (!std::isfinite(energy) || energy<=0) { throw std::runtime_error("silent or nonfinite probe"); }
        if (impulse) { file.write(reinterpret_cast<const char*>(captured.data()),static_cast<std::streamsize>(captured.size()*sizeof(double)));if (!file) { throw std::runtime_error("impulse write failed"); } }
        std::cout<<std::setprecision(17)<<"{\"mode\":\""<<(voice ? "coupled_frozen_voice" : (ports ? "all_ports_diagnostic" : "impulse"))<<"\",\"sample_rate\":"<<model.sample_rate()<<",\"samples\":"<<samples<<",\"stored_modes\":"<<model.modes()<<",\"render_seconds\":"<<elapsed<<",\"realtime_multiple\":"<<(static_cast<double>(samples)/model.sample_rate()/elapsed)<<",\"microphone_rms\":"<<std::sqrt(energy/samples)<<",\"microphone_peak\":"<<peak<<",\"accepted_as_speech_renderer\":false}\n";
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
