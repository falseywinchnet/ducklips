// Incremental stdin control / stdout little-endian float32 stereo at 48 kHz.
// 'render UNTIL' emits through that absolute score time. Stdout is audio only.
#include "ducklips/body_config.hpp"
#include "ducklips/nudges.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
    static_assert(std::endian::native == std::endian::little);
    ducklips::VoiceProfile voice = ducklips::baritone_voice();
    double speed = 1.0;
    ducklips::SynthesisOptions options{};
    std::string error{};
    for (int i = 1; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--conservative-tongue") {
            options.conservative_tongue=true;
        } else if (flag == "--local-flow") {
            options.local_flow=true;
            options.oral_sections=24;options.nasal_sections=9;options.local_substeps=1;
        } else if (flag == "--canonical-motion") {
            options.canonical_motion=true;
        } else if (flag == "--spatial-pressure") {
            options.spatial_pressure=true;
        } else if (flag == "--body" && i + 1 < argc) {
            if (!ducklips::load_body_config(argv[++i], voice, error)) { std::cerr << error << '\n'; return 2; }
        } else if (flag == "--speed" && i + 1 < argc) {
            std::istringstream value(argv[++i]);
            std::string extra{};
            if (!(value >> speed) || (value >> extra) || !std::isfinite(speed) || speed < .25 || speed > 4) { return 2; }
        } else {
            std::cerr << "usage: ducklips_stream [--body FILE] [--speed FACTOR] [--local-flow] [--canonical-motion] [--conservative-tongue] [--spatial-pressure] < nudges > stereo.f32\n";
            return 2;
        }
    }
    if(options.local_flow && (options.canonical_motion || options.spatial_pressure)){std::cerr<<"--local-flow selects a distinct engine\n";return 2;}
    if (options.conservative_tongue && !options.canonical_motion && !options.local_flow) { std::cerr << "--conservative-tongue requires --local-flow or --canonical-motion\n";return 2; }
    if(options.spatial_pressure && !options.canonical_motion){std::cerr<<"--spatial-pressure requires --canonical-motion\n";return 2;}
    std::cerr << "tongue body: " << (options.conservative_tongue?"conservative material map":"rigid circle") << '\n';
    ducklips::NudgeStream stream{};
    stream.prepare(voice, speed,4096,options);
    std::cerr << "acoustic path: " << (options.local_flow ? "local air exchanges" : options.canonical_motion ? "canonical moving network" : "compact modal") << '\n';
    std::cerr << "mean airflow: " << (options.spatial_pressure?"spatial pressure and local constriction flow":"single oral cavity") << '\n';
    std::array<float, 960> chunk{};
    std::string line{};
    std::size_t line_number = 0;
    double render_seconds = 0.0;
    while (std::getline(std::cin, line)) {
        ++line_number;
        std::istringstream command(line.substr(0, line.find('#')));
        std::string word{};
        command >> word;
        if (word != "render") {
            if (!stream.submit(line, error)) {
                std::cerr << "line " << line_number << ": " << error << '\n'; return 2;
            }
            continue;
        }
        double until = 0.0;
        std::string extra{};
        if (!(command >> until) || (command >> extra) || !std::isfinite(until) || until < 0 || until > 86400) {
            std::cerr << "line " << line_number << ": invalid render time\n"; return 2;
        }
        ducklips::SampleCount end = ducklips::samples_from_seconds({}, until / speed);
        if (end < stream.engine().position()) { std::cerr << "cannot render backwards\n"; return 2; }
        while (stream.engine().position() < end) {
            std::size_t count = static_cast<std::size_t>(std::min<ducklips::SampleCount>(480, end - stream.engine().position()));
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            stream.render(std::span<float>(chunk.data(), count * 2));
            if (stream.engine().statistics().geometry_failures) { std::cerr << "geometry audit failed: " << stream.engine().statistics().geometry_error << '\n';return 1; }
            if(stream.engine().statistics().airflow_failures){std::cerr<<"spatial airflow solve failed at "<<stream.engine().statistics().airflow_failure_time<<" s\n";return 1;}
            render_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (std::fwrite(chunk.data(), sizeof(float), count * 2, stdout) != count * 2) { return 1; }
        }
        if (std::fflush(stdout) != 0) { return 1; }
    }
    if (stream.pending() != 0) { std::cerr << "unrendered nudges remain; end with render UNTIL\n"; return 2; }
    const ducklips::SynthesisStatistics& statistics = stream.engine().statistics();
    if (statistics.nonfinite_samples || statistics.clipped_samples) {
        std::cerr << "render failed output audit: nonfinite " << statistics.nonfinite_samples
                  << ", clipped " << statistics.clipped_samples << '\n'; return 1;
    }
    double seconds = ducklips::seconds_from_samples({}, stream.engine().position());
    std::cerr << "audio " << seconds << " s, render " << render_seconds << " s, speed "
              << (render_seconds > 0 ? seconds / render_seconds : 0) << "x realtime\n";
    return 0;
}
