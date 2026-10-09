// Drives DuckLips through the GUI.Forms AudioEngine as a generator voice.
//   offline (default): renders through the engine's mixer without a device and
//                      compares the result with a direct Synthesizer render.
//   --play:            opens the audio device and plays the utterance.
#include "ducklips/gui_forms_voice.hpp"
#include "ducklips/planner.hpp"
#include "ducklips/wav.hpp"
#include "gui_forms/audio/audio.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace {

const char* status_name(gui_forms::AudioStatus status) {
    switch (status) {
    case gui_forms::AudioStatus::ok: return "ok";
    case gui_forms::AudioStatus::invalid_format: return "invalid_format";
    case gui_forms::AudioStatus::invalid_value: return "invalid_value";
    case gui_forms::AudioStatus::file_error: return "file_error";
    case gui_forms::AudioStatus::allocation_failed: return "allocation_failed";
    case gui_forms::AudioStatus::quota_exceeded: return "quota_exceeded";
    case gui_forms::AudioStatus::device_unavailable: return "device_unavailable";
    case gui_forms::AudioStatus::backend_error: return "backend_error";
    case gui_forms::AudioStatus::closed: return "closed";
    case gui_forms::AudioStatus::cancelled: return "cancelled";
    }
    return "unknown";
}

// Feeds cues that fall within the lookahead window of the render position.
void feed(ducklips::Synthesizer& synthesizer, const std::vector<ducklips::Cue>& cues, std::size_t& next,
          ducklips::SampleCount lookahead) {
    while (next < cues.size() && cues[next].time < synthesizer.position() + lookahead) {
        if (!synthesizer.enqueue(cues[next])) {
            break;
        }
        next += 1;
    }
}

}

int main(int argc, char** argv) {
    std::string text = "Hello from inside GUI dot Forms audio.";
    std::string voice_name = "alto";
    std::string output = "out/gui_forms_offline.wav";
    bool play = false;
    for (int index = 1; index < argc; ++index) {
        std::string flag = argv[index];
        if (flag == "--play") {
            play = true;
        } else if (flag == "--text" && index + 1 < argc) {
            text = argv[++index];
        } else if (flag == "--voice" && index + 1 < argc) {
            voice_name = argv[++index];
        } else if (flag == "--out" && index + 1 < argc) {
            output = argv[++index];
        }
    }
    ducklips::VoiceProfile voice{};
    if (!ducklips::voice_by_name(voice_name, voice)) {
        std::fprintf(stderr, "unknown voice\n");
        return 2;
    }
    ducklips::SynthesisOptions options{};

    // Reference: the synthesizer rendered directly.
    ducklips::Synthesizer reference{};
    reference.prepare(voice, options);
    ducklips::UtterancePlanner planner{};
    planner.prepare(voice, reference.vowels(), options.rate);
    ducklips::UtterancePlan plan = planner.plan_text(text);
    const ducklips::SampleCount lookahead = ducklips::samples_from_seconds(options.rate, 0.25);
    const std::size_t total_frames = static_cast<std::size_t>(ducklips::samples_from_seconds(options.rate, plan.duration));
    const std::size_t chunk_frames = 480;
    std::vector<float> direct(2 * total_frames, 0.0f);
    std::size_t next = 0;
    for (std::size_t frame = 0; frame < total_frames; frame += chunk_frames) {
        feed(reference, plan.cues, next, lookahead);
        std::size_t frames = std::min(chunk_frames, total_frames - frame);
        reference.render(std::span<float>(direct.data() + 2 * frame, 2 * frames));
    }

    // The same body as a GUI.Forms generator voice.
    std::shared_ptr<ducklips::GuiFormsVoice> source = std::make_shared<ducklips::GuiFormsVoice>();
    (*source).synthesizer().prepare(voice, options);
    gui_forms::AudioEngine engine{};
    gui_forms::AudioStatus opened = engine.open(!play);
    std::printf("AudioEngine::open(%s): %s\n", play ? "device" : "offline", status_name(opened));
    if (opened != gui_forms::AudioStatus::ok) {
        return 1;
    }
    gui_forms::AudioVoice handle{};
    gui_forms::AudioStatus attached = engine.generator(source, handle);
    gui_forms::AudioStatus started = handle.play();
    std::printf("generator voice: %s, play: %s\n", status_name(attached), status_name(started));
    if (attached != gui_forms::AudioStatus::ok || started != gui_forms::AudioStatus::ok) {
        return 1;
    }
    ducklips::Synthesizer& body = (*source).synthesizer();
    next = 0;
    if (play) {
        // Control thread: keep the cue queue ahead of the device callback.
        while (!body.finished()) {
            feed(body, plan.cues, next, lookahead);
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(300));
        engine.shutdown();
        std::printf("played %.2f s\n", plan.duration);
        return 0;
    }
    std::vector<float> mixed(2 * total_frames, 0.0f);
    for (std::size_t frame = 0; frame < total_frames; frame += chunk_frames) {
        feed(body, plan.cues, next, lookahead);
        std::size_t frames = std::min(chunk_frames, total_frames - frame);
        gui_forms::AudioStatus rendered = engine.render(std::span<float>(mixed.data() + 2 * frame, 2 * frames));
        if (rendered != gui_forms::AudioStatus::ok) {
            std::printf("AudioEngine::render: %s\n", status_name(rendered));
            return 1;
        }
    }
    engine.shutdown();
    double largest_difference = 0.0;
    double energy = 0.0;
    for (std::size_t index = 0; index < mixed.size(); ++index) {
        double difference = static_cast<double>(mixed[index]) - static_cast<double>(direct[index]);
        largest_difference = std::max(largest_difference, std::fabs(difference));
        energy += static_cast<double>(mixed[index]) * static_cast<double>(mixed[index]);
    }
    double rms = std::sqrt(energy / static_cast<double>(std::max<std::size_t>(mixed.size(), 1)));
    std::printf("offline mixer render %.2f s, rms %.4f, largest difference from direct render %.3e\n",
                static_cast<double>(total_frames) / options.rate.samples_per_second, rms, largest_difference);
    ducklips::AudioBuffer audio{};
    audio.samples = mixed;
    if (!ducklips::write_wav_float(output, audio)) {
        return 1;
    }
    std::printf("wrote %s\n", output.c_str());
    return 0;
}
