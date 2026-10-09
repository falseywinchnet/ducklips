// Offline frozen port-response probe. No speech source, normalization, or EQ.
#include "surface_network.hpp"
#include "surface_synthesizer.hpp"
#include "ducklips/body_config.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
void coefficients(std::ostream& out, const std::array<double, 3>& value) {
    out << '[' << value[0] << ',' << value[1] << ',' << value[2] << ']';
}
double number(const std::string& text) {
    std::size_t used = 0;
    double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value)) { throw std::invalid_argument("invalid finite number"); }
    return value;
}
}

int main(int argc, char** argv) {
    try {
        std::string fixture = "body";
        std::string prefix = "surface-response";
        std::string body;
        std::string source = "glottis";
        double duration = 0.25;
        double impulse = 1.0e-10;
        for (int index = 1; index < argc; index += 2) {
            if (index + 1 >= argc) { throw std::invalid_argument("every option needs a value"); }
            std::string option = argv[index], value = argv[index + 1];
            if (option == "--fixture") { fixture = value; }
            else if (option == "--out") { prefix = value; }
            else if (option == "--body") { body = value; }
            else if (option == "--source") { source = value; }
            else if (option == "--duration") { duration = number(value); }
            else if (option == "--impulse") { impulse = number(value); }
            else { throw std::invalid_argument("unknown option: " + option); }
        }
        if (!(duration > 0.0 && duration <= 10.0 && impulse > 0.0 && impulse <= 1.0e-3)) {
            throw std::invalid_argument("duration must be (0,10] seconds and impulse (0,1e-3] m^3/s");
        }
        if (source != "glottis" && source != "lateral") { throw std::invalid_argument("source must be glottis or lateral"); }
        const double rate = 192000.0;
        ducklips::SurfaceNetwork network{};
        std::size_t middle = 0;
        const Clock::time_point prepare_start = Clock::now();
        if (fixture == "body") {
            ducklips::VoiceProfile voice = ducklips::baritone_voice();
            if (!body.empty()) {
                std::string error;
                if (!ducklips::load_body_config(body, voice, error)) { throw std::invalid_argument(error); }
            }
            ducklips::SurfaceSynthesizer synth{};
            synth.prepare(voice, ducklips::SurfaceSynthesisOptions{});
            network = synth.acoustic_network();
            middle = network.first_piece(32);
        } else {
            if (fixture != "tube" && fixture != "solid" && fixture != "bent" && fixture != "floor") {
                throw std::invalid_argument("fixture must be body, tube, solid, bent, or floor");
            }
            std::vector<double> areas(50, 3.0e-4);
            ducklips::CavityGeometry geometry = ducklips::single_tube(areas, 0.175, fixture == "tube" ? 1 : 8);
            ducklips::TubeGeometry& tract = geometry.tubes[0];
            if (fixture != "tube") {
                tract.solid = true;
                tract.bent = fixture == "bent" || fixture == "floor";
                tract.cross_paths = true;
                tract.axial_share = 0.6;
                tract.cross_timing = 0.86;
                tract.perimeter_factor = 1.0;
            }
            if (fixture == "floor") {
                ducklips::TubeGeometry floor{};
                floor.stations = 8;
                floor.length = 8 * 0.175 / 50;
                floor.ring_pieces = 1;
                floor.parallel = true;
                floor.parent = 0;
                floor.parent_station = 35;
                geometry.tubes.push_back(floor);
                geometry.station_area.resize(58, 0.5e-4);
            }
            network.build(geometry, rate, ducklips::AirConstants{});
            if (fixture == "floor") {
                std::vector<double> openings(58, 0.0);
                // Two openings form a genuine loop, not a single side branch.
                openings[55] = 0.25e-4;
                openings[56] = 0.5e-4;
                network.set_openings(openings);
                network.reshape(geometry.station_area, 1);
            }
            middle = network.first_piece(25);
        }
        network.reset();
        const ducklips::SurfaceLinearSnapshot snapshot = network.linear_snapshot();
        std::size_t injection = source == "glottis" ? snapshot.glottis : middle;
        std::vector<std::size_t> observations{snapshot.glottis, middle};
        if (fixture != "tube") { observations.push_back(middle + 4); }
        for (const ducklips::SurfaceLinearPort& port : snapshot.ports) {
            if (port.radiation) { observations.push_back(port.piece); }
        }
        double prepare_seconds = std::chrono::duration<double>(Clock::now() - prepare_start).count();
        std::size_t samples = static_cast<std::size_t>(duration * rate);
        std::size_t channels = observations.size() + 1;
        std::vector<double> response(samples * channels);
        double gather_seconds = 0.0, scatter_seconds = 0.0;
        const Clock::time_point render_start = Clock::now();
        for (std::size_t sample = 0; sample < samples; ++sample) {
            Clock::time_point begin = Clock::now();
            network.gather();
            Clock::time_point gathered = Clock::now();
            if (sample == 0) { network.add_flow(injection, impulse); }
            double microphone = network.scatter(0.0, 0.0, 0.0);
            Clock::time_point scattered = Clock::now();
            gather_seconds += std::chrono::duration<double>(gathered - begin).count();
            scatter_seconds += std::chrono::duration<double>(scattered - gathered).count();
            for (std::size_t channel = 0; channel < observations.size(); ++channel) {
                response[sample * channels + channel] = network.piece_pressure(observations[channel]);
            }
            response[sample * channels + observations.size()] = microphone;
        }
        double render_seconds = std::chrono::duration<double>(Clock::now() - render_start).count();
        for (double value : response) {
            if (!std::isfinite(value)) { throw std::runtime_error("nonfinite observed response"); }
        }
        std::ofstream binary(prefix + ".f64", std::ios::binary);
        binary.write(reinterpret_cast<const char*>(response.data()), static_cast<std::streamsize>(response.size() * sizeof(double)));
        if (!binary) { throw std::runtime_error("cannot write response"); }
        std::ofstream out(prefix + ".json");
        out << std::setprecision(17);
        out << "{\n\"fixture\":\"" << fixture << "\",\"source\":\"" << source << "\",\"sample_rate\":" << rate
            << ",\"pieces\":" << snapshot.pieces << ",\"allocated_delay_samples\":" << snapshot.allocated_delay_samples
            << ",\"source_piece\":" << injection
            << ",\"impulse_m3_per_s\":" << impulse << ",\"samples\":" << samples << ",\"channels\":" << channels
            << ",\"pressure_pieces\":[";
        for (std::size_t i = 0; i < observations.size(); ++i) { if (i) { out << ','; } out << observations[i]; }
        out << "],\"piece_tubes\":[";
        for (std::size_t i = 0; i < snapshot.piece_tubes.size(); ++i) { if (i) { out << ','; } out << snapshot.piece_tubes[i]; }
        out << "],\"last_channel\":\"sum radiation flow derivative\",\"paths\":[\n";
        for (std::size_t i = 0; i < snapshot.paths.size(); ++i) {
            const ducklips::SurfaceLinearPath& path = snapshot.paths[i];
            if (i) { out << ",\n"; }
            out << '[' << path.from << ',' << path.to << ',' << path.admittance << ',' << path.delay << ',';
            coefficients(out, path.numerator); out << ','; coefficients(out, path.denominator); out << ']';
        }
        out << "\n],\"ports\":[\n";
        for (std::size_t i = 0; i < snapshot.ports.size(); ++i) {
            const ducklips::SurfaceLinearPort& port = snapshot.ports[i];
            if (i) { out << ",\n"; }
            out << '[' << port.piece << ',' << port.admittance << ',' << (port.radiation ? "true" : "false") << ',';
            coefficients(out, port.numerator); out << ','; coefficients(out, port.denominator); out << ']';
        }
        out << "\n],\"prepare_seconds\":" << prepare_seconds << ",\"render_seconds\":" << render_seconds
            << ",\"gather_seconds\":" << gather_seconds << ",\"scatter_seconds\":" << scatter_seconds << "}\n";
        if (!out) { throw std::runtime_error("cannot write snapshot"); }
        std::cout << fixture << ' ' << source << ": " << snapshot.pieces << " pieces, " << snapshot.paths.size()
                  << " paths; " << duration / render_seconds << "x; gather " << gather_seconds
                  << " s, scatter " << scatter_seconds << " s\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
