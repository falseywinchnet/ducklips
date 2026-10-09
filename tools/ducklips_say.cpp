// Command-line front end: text or ARPAbet in, a stereo 48 kHz float WAV out.
// Cues are streamed into the engine a short way ahead of the render position,
// the way an interactive host would feed it.
#include "ducklips/planner.hpp"
#include "ducklips/body_config.hpp"
#include "ducklips/synthesizer.hpp"
#include "ducklips/voice.hpp"
#include "ducklips/wav.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Arguments final {
    std::string voice{"baritone"};
    std::string body_path{};
    std::string text{};
    std::string phonemes{};
    std::string output{"ducklips.wav"};
    std::string dictionary{};
    std::string probe_path{};
    bool conservative_tongue{false};
    bool canonical_motion{false};
    std::size_t canonical_substeps{8};
    bool room{false};
    bool report{false};
    double rate{0.0};
    bool normalize{false};
    std::string invariant{"passive"};
    std::size_t control_block{48};
    std::size_t chunk_frames{480};
    bool diagnostics{false};
    bool no_pop{false};
    std::string tap{"microphones"};
    std::string trace_path{};
    double trace_start{-1.0};
    double frication_gain{-1.0};
    double aspiration_gain{-1.0};
    double lung_pressure{-1.0};
};

void print_usage() {
    std::printf("usage: ducklips_say [--voice NAME] [--body FILE] [--rate FACTOR] (--text TEXT | --phonemes ARPABET) [--out FILE.wav]\n"
                "                    [--room] [--normalize] [--report] [--dictionary cmudict.txt]\n"
                "                    [--probe FILE.csv] [--chunk FRAMES] [--control-block SAMPLES] [--diagnostics]\n"
                "                    [--no-pop] [--frication-gain G] [--aspiration-gain G] [--lung-pressure PA]\n"
                "                    [--trace START_S FILE.csv] [--tap microphones|lips|source]\n"
                "                    [--invariant passive|charge|pressure]\n"
                "[--canonical-motion] [--conservative-tongue] [--canonical-substeps 1..64]\n"
                "voices: baritone alto child gravel whisper\n");
}

bool parse_arguments(int count, char** values, Arguments& arguments) {
    for (int index = 1; index < count; ++index) {
        std::string flag = values[index];
        bool has_value = index + 1 < count;
        if (flag == "--conservative-tongue") {
            arguments.conservative_tongue=true;
        } else if (flag == "--canonical-motion") {
            arguments.canonical_motion=true;
        } else if (flag == "--canonical-substeps" && has_value) {
            std::istringstream parsed(values[++index]);std::string extra;
            if (!(parsed>>arguments.canonical_substeps) || (parsed>>extra) || arguments.canonical_substeps<1 || arguments.canonical_substeps>64) { return false; }
        } else if (flag == "--rate" && has_value) {
            std::istringstream parsed(values[++index]);
            std::string extra{};
            if (!(parsed >> arguments.rate) || (parsed >> extra) || !std::isfinite(arguments.rate)
                || arguments.rate < 0.25 || arguments.rate > 4.0) { return false; }
        } else if (flag == "--body" && has_value) {
            arguments.body_path = values[++index];
        } else if (flag == "--voice" && has_value) {
            arguments.voice = values[++index];
        } else if (flag == "--text" && has_value) {
            arguments.text = values[++index];
        } else if (flag == "--phonemes" && has_value) {
            arguments.phonemes = values[++index];
        } else if (flag == "--out" && has_value) {
            arguments.output = values[++index];
        } else if (flag == "--dictionary" && has_value) {
            arguments.dictionary = values[++index];
        } else if (flag == "--probe" && has_value) {
            arguments.probe_path = values[++index];
        } else if (flag == "--control-block" && has_value) {
            arguments.control_block = static_cast<std::size_t>(std::atoi(values[++index]));
        } else if (flag == "--chunk" && has_value) {
            arguments.chunk_frames = static_cast<std::size_t>(std::atoi(values[++index]));
        } else if (flag == "--trace" && index + 2 < count) {
            arguments.trace_start = std::atof(values[++index]);
            arguments.trace_path = values[++index];
        } else if (flag == "--tap" && has_value) {
            arguments.tap = values[++index];
        } else if (flag == "--no-pop") {
            arguments.no_pop = true;
        } else if (flag == "--frication-gain" && has_value) {
            arguments.frication_gain = std::atof(values[++index]);
        } else if (flag == "--lung-pressure" && has_value) {
            arguments.lung_pressure = std::atof(values[++index]);
        } else if (flag == "--aspiration-gain" && has_value) {
            arguments.aspiration_gain = std::atof(values[++index]);
        } else if (flag == "--diagnostics") {
            arguments.diagnostics = true;
        } else if (flag == "--room") {
            arguments.room = true;
        } else if (flag == "--report") {
            arguments.report = true;
        } else if (flag == "--normalize") {
            arguments.normalize = true;
        } else if (flag == "--invariant" && has_value) {
            arguments.invariant = values[++index];
        } else {
            return false;
        }
    }
    bool has_input = !arguments.text.empty() || !arguments.phonemes.empty();
    return has_input && (arguments.tap == "microphones" || arguments.tap == "lips" || arguments.tap == "source")
        && arguments.control_block >= 8 && arguments.control_block <= 480
        && arguments.chunk_frames > 0 && arguments.chunk_frames <= 48000;
}

}

int main(int argc, char** argv) {
    Arguments arguments{};
    if (!parse_arguments(argc, argv, arguments)) {
        print_usage();
        return 2;
    }
    ducklips::VoiceProfile voice{};
    if (!ducklips::voice_by_name(arguments.voice, voice)) {
        std::fprintf(stderr, "unknown voice %s\n", arguments.voice.c_str());
        return 2;
    }
    if (!arguments.body_path.empty()) {
        std::string error{};
        if (!ducklips::load_body_config(arguments.body_path, voice, error)) {
            std::fprintf(stderr, "%s\n", error.c_str());
            return 2;
        }
    }
    if (arguments.frication_gain >= 0.0) {
        voice.frication_gain = arguments.frication_gain;
    }
    if (arguments.lung_pressure >= 0.0) {
        voice.lung_pressure = arguments.lung_pressure;
    }
    if (arguments.aspiration_gain >= 0.0) {
        voice.aspiration_gain = arguments.aspiration_gain;
    }
    if (arguments.rate > 0.0) { voice.speaking_rate = arguments.rate; }
    ducklips::SynthesisOptions options{};
    options.conservative_tongue=arguments.conservative_tongue;
    if (options.conservative_tongue && !arguments.canonical_motion) { std::fprintf(stderr,"--conservative-tongue requires --canonical-motion\n");return 2; }
    options.canonical_motion=arguments.canonical_motion;
    options.canonical_substeps=arguments.canonical_substeps;
    std::printf("tongue body: %s\n",options.conservative_tongue?"conservative material map":"rigid circle");
    std::printf("acoustic path: %s\n",options.canonical_motion ? "canonical moving network (causal-loss integration reference)" : "compact modal");
    options.radiate_mean_flow = !arguments.no_pop;
    if (arguments.tap == "lips") {
        options.tap = ducklips::OutputTap::lips;
    } else if (arguments.tap == "source") {
        options.tap = ducklips::OutputTap::glottal_source;
    }
    if (!arguments.trace_path.empty()) {
        options.trace_start = arguments.trace_start;
        options.trace_samples = 48000 * 8;
    }
    options.air_path.room.enabled = arguments.room;
    options.control_block = arguments.control_block;
    options.diagnostics = arguments.diagnostics;
    if (arguments.invariant == "charge") {
        options.invariant = ducklips::ReprojectionInvariant::charge_momentum;
    } else if (arguments.invariant == "pressure") {
        options.invariant = ducklips::ReprojectionInvariant::pressure_flow;
    }

    std::chrono::steady_clock::time_point prepare_start = std::chrono::steady_clock::now();
    ducklips::Synthesizer synthesizer{};
    synthesizer.prepare(voice, options);
    std::chrono::steady_clock::time_point prepare_end = std::chrono::steady_clock::now();
    double prepare_seconds = std::chrono::duration<double>(prepare_end - prepare_start).count();
    std::printf("voice %s: tract %.1f cm, folds %.1f Hz at Q=1, prepared in %.3f s\n", voice.name.c_str(),
                voice.tract_length * 100.0, synthesizer.unit_tension_pitch(), prepare_seconds);
    if (arguments.report) {
        const std::vector<ducklips::FormantFit>& fits = synthesizer.vowels().fits();
        std::printf("vowel  target F1 F2 F3 (Hz)        posture F1 F2 F3 (Hz)        worst  iterations\n");
        for (std::size_t index = 0; index < fits.size(); ++index) {
            const ducklips::FormantFit& fit = fits[index];
            std::printf("%-5s  %5.0f %5.0f %5.0f            %5.0f %5.0f %5.0f            %5.1f%%  %zu\n",
                        fit.phone == ducklips::Phone::count ? "locus" : ducklips::phone_spec(fit.phone).name, fit.target[0], fit.target[1], fit.target[2],
                        fit.achieved[0], fit.achieved[1], fit.achieved[2], 100.0 * fit.largest_relative_error,
                        fit.iterations);
        }
    }

    ducklips::UtterancePlanner planner{};
    planner.prepare(voice, synthesizer.vowels(), options.rate);
    if (!arguments.dictionary.empty()) {
        std::size_t added = planner.load_dictionary(arguments.dictionary);
        std::printf("dictionary: %zu definitions installed from %s\n", added, arguments.dictionary.c_str());
    }
    ducklips::UtterancePlan plan{};
    if (!arguments.text.empty()) {
        plan = planner.plan_text(arguments.text);
    } else {
        plan = planner.plan_phonemes(arguments.phonemes);
    }
    std::printf("phones: %s\n", plan.transcription.c_str());
    for (std::size_t index = 0; index < plan.transitions.size(); ++index) {
        std::printf("transition %s\n", plan.transitions[index].c_str());
    }
    for (std::size_t index = 0; index < plan.guessed_words.size(); ++index) {
        std::printf("letter-to-sound guess: %s\n", plan.guessed_words[index].c_str());
    }
    std::printf("plan: %zu cues, %.2f s\n", plan.cues.size(), plan.duration);
    if (arguments.report) {
        std::printf("timeline:");
        for (std::size_t index = 0; index < plan.phones.size(); ++index) {
            std::printf(" %s@%.3f", ducklips::phone_spec(plan.phones[index].phone).name, plan.phones[index].start);
        }
        std::printf("\n");
    }

    std::FILE* probe = nullptr;
    if (!arguments.probe_path.empty()) {
        probe = std::fopen(arguments.probe_path.c_str(), "w");
        if (probe != nullptr) {
            std::fprintf(probe, "time,flow_cm3s,glottal_area_cm2,intraoral_pa,lung_pa,intended_hz,measured_hz,"
                                "tension,narrowest_cm2,modes,m1_hz,m2_hz,m3_hz,m4_hz,energy,transglottal_pa,p0_pa");
            if (arguments.diagnostics) {
                std::fprintf(probe,",oral_volume_m3,floor_volume_m3,volume_without_body_m3,volume_without_tip_m3,tongue_x_cm,tongue_y_cm,tip_x_cm,tip_y_cm,body_activation,tip_activation");
            }
            std::fprintf(probe,"\n");
        }
    }

    const std::size_t chunk_frames = std::max<std::size_t>(arguments.chunk_frames, 1);
    const ducklips::SampleCount lookahead = ducklips::samples_from_seconds(options.rate, 0.25);
    const ducklips::SampleCount limit = ducklips::samples_from_seconds(options.rate, plan.duration + 2.0);
    std::vector<float> chunk(2 * chunk_frames, 0.0f);
    ducklips::AudioBuffer audio{};
    audio.channels = 2;
    audio.sample_rate = options.rate.samples_per_second;
    std::size_t next_cue = 0;
    double render_seconds = 0.0;
    while (!synthesizer.finished() && synthesizer.position() < limit) {
        while (next_cue < plan.cues.size() && plan.cues[next_cue].time < synthesizer.position() + lookahead
               && synthesizer.queue_space() > 0) {
            if (!synthesizer.enqueue(plan.cues[next_cue])) {
                break;
            }
            next_cue += 1;
        }
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        synthesizer.render(std::span<float>(chunk));
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        render_seconds += std::chrono::duration<double>(end - start).count();
        audio.samples.insert(audio.samples.end(), chunk.begin(), chunk.end());
        if (probe != nullptr) {
            const ducklips::SynthesisProbe& state = synthesizer.probe();
            double time = ducklips::seconds_from_samples(options.rate, synthesizer.position());
            std::fprintf(probe, "%.4f,%.2f,%.4f,%.1f,%.1f,%.1f,%.1f,%.3f,%.4f,%zu,%.0f,%.0f,%.0f,%.0f,%.4e,%.1f,%.1f", time,
                         state.glottal_flow * 1.0e6, state.glottal_area * 1.0e4, state.intraoral_pressure,
                         state.lung_pressure, state.intended_pitch, state.measured_pitch, state.tension,
                         state.narrowest_area * 1.0e4, state.modes, state.lowest_mode_hz[0],
                         state.lowest_mode_hz[1], state.lowest_mode_hz[2], state.lowest_mode_hz[3],
                         state.acoustic_energy, state.transglottal, state.glottal_pressure);
            if (arguments.diagnostics) {
                std::fprintf(probe,",%.12e,%.12e,%.12e,%.12e,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f",
                             state.oral_volume,state.floor_volume,state.volume_without_body,state.volume_without_tip,
                             state.tongue_centre.x,state.tongue_centre.y,state.tongue_tip.x,state.tongue_tip.y,
                             state.body_activation,state.tip_activation);
            }
            std::fprintf(probe,"\n");
        }
    }
    if (probe != nullptr) {
        std::fclose(probe);
    }

    double peak = 0.0;
    double energy = 0.0;
    for (std::size_t index = 0; index < audio.samples.size(); ++index) {
        double value = static_cast<double>(audio.samples[index]);
        peak = std::max(peak, std::fabs(value));
        energy += value * value;
    }
    double audio_seconds = static_cast<double>(audio.samples.size() / 2) / options.rate.samples_per_second;
    double rms = std::sqrt(energy / std::max<double>(1.0, static_cast<double>(audio.samples.size())));
    if (arguments.normalize && peak > 0.0) {
        float gain = static_cast<float>(0.89 / peak);
        for (std::size_t index = 0; index < audio.samples.size(); ++index) {
            audio.samples[index] *= gain;
        }
    }
    const ducklips::SynthesisStatistics& statistics = synthesizer.statistics();
    std::printf("output audit: nonfinite %zu, clipped %zu\n", synthesizer.statistics().nonfinite_samples, synthesizer.statistics().clipped_samples);
    std::printf("rendered %.2f s in %.3f s: real-time factor %.3f (lower is faster)\n", audio_seconds, render_seconds,
                render_seconds / audio_seconds);
    std::printf("geometry steps %zu (at most %zu per block), modes %zu..%zu, warm-started %zu, bisected %zu, "
                "factorizations/step %.1f\n",
                statistics.blocks, statistics.most_geometry_steps, statistics.fewest_modes, statistics.most_modes, statistics.warm_modes,
                statistics.bisected_modes,
                static_cast<double>(statistics.factorizations) / std::max<double>(1.0, static_cast<double>(statistics.blocks)));
    if (options.conservative_tongue) {
        std::printf("geometry solve audit: failures %zu; sampled body/wall overlap %.9g cm\n",statistics.geometry_failures,statistics.maximum_body_wall_penetration);
        std::printf("maximum constrained tip target error: %.9g cm\n",statistics.maximum_tip_target_error);
    }
    if (options.canonical_motion) { std::printf("moving cavity maximum mass residual: %.9e m^3/s\n",statistics.maximum_cavity_mass_residual); }
    std::printf("fold cycles %zu, peak %.3f, rms %.4f (before normalization)\n", statistics.pitch_cycles, peak, rms);
    if (arguments.diagnostics && !options.canonical_motion) {
        std::printf("worst basis orthogonality defect %.3e at %.3f s\n", statistics.worst_orthogonality,
                    statistics.worst_orthogonality_time);
    }
    if (!arguments.trace_path.empty()) {
        std::FILE* trace = std::fopen(arguments.trace_path.c_str(), "w");
        if (trace != nullptr) {
            std::fprintf(trace, "time,lips,nostrils,throat,p0,energy,modes,area,flow,tension,geometry_lip_metric_jump,geometry_lip_projection_jump,geometry_pressure_jump,geometry_lip_slope_jump,geometry_energy_ratio,geometry_parameter_change,frication_pressure,frication_branch,constriction_area,cavity_pressure,cavity_displacement_flow,cavity_mass_residual\n");
            const std::vector<ducklips::TraceSample>& samples = synthesizer.trace();
            for (std::size_t index = 0; index < samples.size(); ++index) {
                double time = arguments.trace_start + static_cast<double>(index) / options.rate.samples_per_second;
                std::fprintf(trace, "%.6f,%.6e,%.6e,%.6e,%.3f,%.6e,%zu,%.6e,%.6e,%.5f,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%zu,%.9e,%.9e,%.9e,%.9e\n", time, samples[index].lips,
                             samples[index].nostrils, samples[index].throat, samples[index].glottal_pressure,
                             samples[index].energy, samples[index].modes, samples[index].glottal_area,
                             samples[index].glottal_flow, samples[index].tension, samples[index].geometry_lip_metric_jump,
                             samples[index].geometry_lip_projection_jump, samples[index].geometry_pressure_jump,
                             samples[index].geometry_lip_slope_jump, samples[index].geometry_energy_ratio,
                             samples[index].geometry_parameter_change, samples[index].frication_pressure, samples[index].frication_branch,
                             samples[index].constriction_area, samples[index].cavity_pressure,
                             samples[index].cavity_displacement_flow,samples[index].cavity_mass_residual);
            }
            std::fclose(trace);
        }
    }
    if (!ducklips::write_wav_float(arguments.output, audio)) {
        std::fprintf(stderr, "could not write %s\n", arguments.output.c_str());
        return 1;
    }
    std::printf("wrote %s\n", arguments.output.c_str());
    if (synthesizer.statistics().geometry_failures) {
        std::fprintf(stderr,"geometry audit failed at %.6f s: %s; diagnostic WAV retained\n",synthesizer.statistics().geometry_failure_time,synthesizer.statistics().geometry_error);return 1;
    }
    if (synthesizer.statistics().nonfinite_samples || synthesizer.statistics().clipped_samples) {
        std::fprintf(stderr,"output audit failed; diagnostic WAV retained\n");return 1;
    }
    return 0;
}
