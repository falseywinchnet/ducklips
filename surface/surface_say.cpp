// Speaks text or ARPAbet through the surface engine: the DuckLips planner and
// body, with the surface network as the acoustics.
#include "surface_synthesizer.hpp"
#include "ducklips/planner.hpp"
#include "ducklips/body_config.hpp"
#include "ducklips/voice.hpp"
#include "ducklips/wav.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Arguments final {
    std::string voice{"baritone"};
    std::string body_path{};
    std::string text{};
    std::string phonemes{};
    std::string output{"surface_say.wav"};
    std::string dictionary{};
    std::string probe_path{};
    std::string area_path{};
    bool normalize{false};
    bool list_cues{false};
    bool steady{false};
    double pitch_correction{30.0};
    double asymmetry{-1.0};
    double lungs{-1.0};
    double base_pitch{-1.0};
    double rate{-1.0};
    double accent{-1.0};
    double adducted{-1.0e9};
    double tremor{0.002};
    double adduction{0.08};
    std::size_t ring{8};
    bool solid{true};
    bool floor{true};
    double floor_opening{1.0};
    double floor_minimum{0.0};
    double roundness{-1.0};
    double softness{-1.0};
    double share{0.6};
    double timing{0.86};
    double tension{0.0};
    double frication_scale{40.0};
    double aspiration_scale{4.8};
    unsigned radiators{7u};
    double trace_start{-1.0};
    std::string trace_path{};
};

void print_usage() {
    std::printf("usage: surface_say [--voice NAME] [--body FILE] [--rate FACTOR] (--text TEXT | --phonemes ARPABET) [--out FILE.wav]\n"
                "                   [--normalize] [--dictionary cmudict.txt] [--ring PIECES] [--solid | --tube] [--no-floor] [--share A] [--timing T]\n"
                "                   [--tension N_PER_M] [--probe FILE.csv]\n"
                "                   [--frication-scale S] [--aspiration-scale S]\n"
                "voices: baritone alto child gravel whisper\n");
}

bool parse_arguments(int count, char** values, Arguments& arguments) {
    for (int index = 1; index < count; ++index) {
        std::string flag = values[index];
        bool has_value = index + 1 < count;
        if (flag == "--body" && has_value) {
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
        } else if (flag == "--areas" && has_value) {
            arguments.area_path = values[++index];
        } else if (flag == "--probe" && has_value) {
            arguments.probe_path = values[++index];
        } else if (flag == "--ring" && has_value) {
            arguments.ring = static_cast<std::size_t>(std::atoi(values[++index]));
        } else if (flag == "--roundness" && has_value) {
            arguments.roundness = std::atof(values[++index]) * 1.0e-6;
        } else if (flag == "--softness" && has_value) {
            arguments.softness = std::atof(values[++index]) * 1.0e-6;
        } else if (flag == "--solid") {
            arguments.solid = true;
        } else if (flag == "--floor-opening" && has_value) {
            arguments.floor_opening = std::atof(values[++index]);
        } else if (flag == "--floor-minimum" && has_value) {
            arguments.floor_minimum = std::atof(values[++index]) * 1.0e-4;
        } else if (flag == "--no-floor") {
            arguments.floor = false;
        } else if (flag == "--tube") {
            arguments.solid = false;
            arguments.ring = 1;
        } else if (flag == "--share" && has_value) {
            arguments.share = std::atof(values[++index]);
        } else if (flag == "--timing" && has_value) {
            arguments.timing = std::atof(values[++index]);
        } else if (flag == "--tension" && has_value) {
            arguments.tension = std::atof(values[++index]);
        } else if (flag == "--frication-scale" && has_value) {
            arguments.frication_scale = std::atof(values[++index]);
        } else if (flag == "--aspiration-scale" && has_value) {
            arguments.aspiration_scale = std::atof(values[++index]);
        } else if (flag == "--trace" && index + 2 < count) {
            arguments.trace_start = std::atof(values[++index]);
            arguments.trace_path = values[++index];
        } else if (flag == "--radiators" && has_value) {
            arguments.radiators = static_cast<unsigned>(std::atoi(values[++index]));
        } else if (flag == "--pitch-correction" && has_value) {
            arguments.pitch_correction = std::atof(values[++index]);
        } else if (flag == "--tremor" && has_value) {
            arguments.tremor = std::atof(values[++index]);
        } else if (flag == "--adduction" && has_value) {
            arguments.adduction = std::atof(values[++index]);
        } else if (flag == "--rate" && has_value) {
            arguments.rate = std::atof(values[++index]);
        } else if (flag == "--accent" && has_value) {
            arguments.accent = std::atof(values[++index]);
        } else if (flag == "--base-pitch" && has_value) {
            arguments.base_pitch = std::atof(values[++index]);
        } else if (flag == "--adducted" && has_value) {
            arguments.adducted = std::atof(values[++index]) * 1.0e-6;
        } else if (flag == "--lungs" && has_value) {
            arguments.lungs = std::atof(values[++index]);
        } else if (flag == "--asymmetry" && has_value) {
            arguments.asymmetry = std::atof(values[++index]);
        } else if (flag == "--steady") {
            arguments.steady = true;
        } else if (flag == "--cues") {
            arguments.list_cues = true;
        } else if (flag == "--normalize") {
            arguments.normalize = true;
        } else {
            return false;
        }
    }
    bool has_input = !arguments.text.empty() || !arguments.phonemes.empty();
    return has_input && arguments.ring >= 1;
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
    if (arguments.rate > 0.0) {
        voice.speaking_rate = arguments.rate;
    }
    if (arguments.accent >= 0.0) {
        voice.accent_rise = arguments.accent;
    }
    if (arguments.base_pitch > 0.0) {
        voice.base_pitch = arguments.base_pitch;
    }
    if (arguments.adducted > -1.0e8) {
        voice.adducted_area = arguments.adducted;
    }
    if (arguments.lungs > 0.0) {
        voice.lung_pressure = arguments.lungs;
    }
    if (arguments.asymmetry > 0.0) {
        voice.folds.asymmetry = arguments.asymmetry;
    }
    if (arguments.steady) {
        // Diagnostics: no physiological unsteadiness at all.
        voice.jitter = 0.0;
        voice.shimmer = 0.0;
        voice.tremor_depth = 0.0;
    }
    ducklips::SurfaceSynthesisOptions options{};
    options.ring_pieces = arguments.ring;
    options.solid_airway = arguments.solid;
    options.floor = arguments.floor;
    options.floor_opening_scale = arguments.floor_opening;
    options.floor_minimum_area = arguments.floor_minimum;
    if (!arguments.body_path.empty()) {
        options.medial_roundness = voice.folds.medial_roundness;
        options.contact_softness = voice.folds.contact_softness;
    }
    if (arguments.roundness >= 0.0) { options.medial_roundness = arguments.roundness; }
    if (arguments.softness >= 0.0) { options.contact_softness = arguments.softness; }
    options.axial_share = arguments.share;
    options.cross_timing = arguments.timing;
    options.edge_tension = arguments.tension;
    options.frication_scale = arguments.frication_scale;
    options.aspiration_scale = arguments.aspiration_scale;
    options.radiator_mask = arguments.radiators;
    options.pitch_correction_rate = arguments.pitch_correction;
    options.tremor = arguments.tremor;
    options.adduction_wander = arguments.adduction;
    if (!arguments.trace_path.empty()) {
        options.trace_start = arguments.trace_start;
        options.trace_samples = 192000 / 2; // 0.5 s at the internal rate
    }

    std::chrono::steady_clock::time_point prepare_start = std::chrono::steady_clock::now();
    ducklips::SurfaceSynthesizer synthesizer{};
    synthesizer.prepare(voice, options);
    std::chrono::steady_clock::time_point prepare_end = std::chrono::steady_clock::now();
    std::printf("voice %s: tract %.1f cm, folds %.1f Hz at Q=1, %zu pieces, %zu paths, prepared in %.2f s\n",
                voice.name.c_str(), voice.tract_length * 100.0, synthesizer.unit_tension_pitch(),
                synthesizer.piece_count(), synthesizer.link_count(),
                std::chrono::duration<double>(prepare_end - prepare_start).count());

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
    std::printf("plan: %zu cues, %.2f s\n", plan.cues.size(), plan.duration);
    if (arguments.list_cues) {
        const char* kinds[] = {"vowel", "constriction", "release", "velum", "glottis", "lungs", "pitch", "effort",
                               "loudness", "finish"};
        for (std::size_t index = 0; index < plan.cues.size(); ++index) {
            const ducklips::Cue& cue = plan.cues[index];
            if (cue.kind == ducklips::CueKind::pitch) {
                continue;
            }
            std::printf("cue %.3f %s slot %zu value %.4g rate %.0f\n",
                        ducklips::seconds_from_samples(options.rate, cue.time),
                        kinds[static_cast<std::size_t>(cue.kind)], static_cast<std::size_t>(cue.slot), cue.value,
                        cue.angular_rate);
        }
    }
    std::printf("timeline:");
    for (std::size_t index = 0; index < plan.phones.size(); ++index) {
        std::printf(" %s@%.3f", ducklips::phone_spec(plan.phones[index].phone).name, plan.phones[index].start);
    }
    std::printf("\n");

    std::FILE* probe = nullptr;
    if (!arguments.probe_path.empty()) {
        probe = std::fopen(arguments.probe_path.c_str(), "w");
        if (probe != nullptr) {
            std::fprintf(probe, "time,flow_cm3s,glottal_area_cm2,intraoral_pa,lung_pa,intended_hz,measured_hz,"
                                "narrowest_cm2,jet_cm3s,jet_section,lip_cm2,port_cm2\n");
        }
    }
    // With a probe, one chunk per control block so every block is recorded.
    const std::size_t chunk_frames = probe != nullptr ? options.control_block : 480;
    const ducklips::SampleCount lookahead = ducklips::samples_from_seconds(options.rate, 0.25);
    const ducklips::SampleCount limit = ducklips::samples_from_seconds(options.rate, plan.duration + 2.0);
    std::vector<float> chunk(2 * chunk_frames, 0.0f);
    ducklips::AudioBuffer audio{};
    audio.channels = 2;
    audio.sample_rate = options.rate.samples_per_second;
    std::size_t next_cue = 0;
    double render_seconds = 0.0;
    // --areas FILE: time (s), then the airway's areas (cm^2, glottis to
    // lips), every chunk (10 ms).
    std::FILE* area_file = nullptr;
    if (!arguments.area_path.empty()) {
        area_file = std::fopen(arguments.area_path.c_str(), "w");
    }
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
        if (area_file != nullptr) {
            std::span<const double> areas = synthesizer.station_areas();
            std::fprintf(area_file, "%.4f", ducklips::seconds_from_samples(options.rate, synthesizer.position()));
            for (std::size_t section = 0; section < options.oral_sections; ++section) {
                std::fprintf(area_file, " %.4f", areas[section] * 1.0e4);
            }
            // Then every other station (nose, side cavities, the space
            // under the tongue), after a separator.
            std::fprintf(area_file, " |");
            for (std::size_t station = options.oral_sections; station < areas.size(); ++station) {
                std::fprintf(area_file, " %.4f", areas[station] * 1.0e4);
            }
            std::fprintf(area_file, "\n");
        }
        if (probe != nullptr) {
            const ducklips::SurfaceProbe& state = synthesizer.probe();
            double time = ducklips::seconds_from_samples(options.rate, synthesizer.position());
            std::fprintf(probe, "%.4f,%.2f,%.4f,%.1f,%.1f,%.1f,%.1f,%.4f,%.2f,%zu,%.4f,%.4f\n", time,
                         state.glottal_flow * 1.0e6, state.glottal_area * 1.0e4, state.intraoral_pressure,
                         state.lung_pressure, state.intended_pitch, state.measured_pitch,
                         state.narrowest_area * 1.0e4, state.jet_flow * 1.0e6, state.jet_section,
                         state.lip_area * 1.0e4, state.port_area * 1.0e4);
        }
    }
    if (probe != nullptr) {
        std::fclose(probe);
    }
    double peak = 0.0;
    for (std::size_t index = 0; index < audio.samples.size(); ++index) {
        peak = std::max(peak, std::fabs(static_cast<double>(audio.samples[index])));
    }
    if (arguments.normalize && peak > 0.0) {
        float gain = static_cast<float>(0.89 / peak);
        for (std::size_t index = 0; index < audio.samples.size(); ++index) {
            audio.samples[index] *= gain;
        }
    }
    double audio_seconds = static_cast<double>(audio.samples.size() / 2) / options.rate.samples_per_second;
    std::printf("rendered %.2f s in %.2f s: real-time factor %.2f (lower is faster); peak %.3f\n", audio_seconds,
                render_seconds, render_seconds / audio_seconds, peak);
    if (!ducklips::write_wav_float(arguments.output, audio)) {
        std::fprintf(stderr, "cannot write %s\n", arguments.output.c_str());
        return 1;
    }
    std::printf("wrote %s\n", arguments.output.c_str());
    if (area_file != nullptr) {
        std::fclose(area_file);
    }
    if (!arguments.trace_path.empty()) {
        std::FILE* trace = std::fopen(arguments.trace_path.c_str(), "w");
        if (trace != nullptr) {
            std::fprintf(trace, "sample,lip_area_cm2,lip_resistance,lip_flow_cm3s,glottal_flow_cm3s,frication_pa,aspiration_pa,entry_last,entry_before,lip_pa\n");
            const std::vector<ducklips::SurfaceTraceSample>& samples = synthesizer.trace();
            for (std::size_t index = 0; index < samples.size(); ++index) {
                const ducklips::SurfaceTraceSample& sample = samples[index];
                std::fprintf(trace, "%zu,%.6f,%.5f,%.5f,%.3f,%.4f,%.4f,%.5f,%.5f,%.3f\n", index,
                             sample.lip_area * 1.0e4, sample.lip_resistance, sample.lip_flow * 1.0e6,
                             sample.glottal_flow * 1.0e6, sample.frication, sample.aspiration, sample.entry_last,
                             sample.entry_before, sample.lip_pressure);
            }
            std::fclose(trace);
        }
    }
    return 0;
}
