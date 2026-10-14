// Plucks the surface network with a glottal pulse train and records what a
// microphone hears. Shapes come from area files ("area_cm2 length_m" lines);
// two files make a glide from the first shape to the second.
#include "surface_network.hpp"
#include "ducklips/wav.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.14159265358979323846;

struct Profile final {
    std::vector<double> area{};   // m^2
    std::vector<double> length{}; // m
};

Profile read_profile(const std::string& path) {
    Profile profile{};
    std::FILE* file = std::fopen(path.c_str(), "r");
    double area = 0.0;
    double length = 0.0;
    while (file != nullptr && std::fscanf(file, "%lf %lf", &area, &length) == 2) {
        profile.area.push_back(area * 1.0e-4);
        profile.length.push_back(length);
    }
    if (file != nullptr) {
        std::fclose(file);
    }
    return profile;
}

double total_length(const Profile& profile) {
    double total = 0.0;
    for (std::size_t index = 0; index < profile.length.size(); ++index) {
        total += profile.length[index];
    }
    return total;
}

void resample(const Profile& profile, std::size_t stations, std::vector<double>& output) {
    double total = total_length(profile);
    output.assign(stations, 0.0);
    for (std::size_t station = 0; station < stations; ++station) {
        double begin = total * static_cast<double>(station) / static_cast<double>(stations);
        double end = total * static_cast<double>(station + 1) / static_cast<double>(stations);
        double running = 0.0;
        double covered = 0.0;
        double position = 0.0;
        for (std::size_t index = 0; index < profile.area.size(); ++index) {
            double low = std::max(begin, position);
            double high = std::min(end, position + profile.length[index]);
            if (high > low) {
                running += profile.area[index] * (high - low);
                covered += high - low;
            }
            position += profile.length[index];
        }
        output[station] = running / std::max(covered, 1.0e-12);
    }
}

// Two cascaded one-pole low-passes (unity gain at low frequency): the
// spectrum of turbulent jet pressure falls above its broadband peak.
struct NoiseShaper final {
    double coefficient{1.0};
    double first{};
    double second{};
};

void prepare_shaper(NoiseShaper& shaper, double corner_hz, double rate) {
    shaper.coefficient = corner_hz > 0.0 ? 1.0 - std::exp(-2.0 * pi * corner_hz / rate) : 1.0;
    shaper.first = 0.0;
    shaper.second = 0.0;
}

double shape_noise(NoiseShaper& shaper, double input) {
    shaper.first += shaper.coefficient * (input - shaper.first);
    shaper.second += shaper.coefficient * (shaper.first - shaper.second);
    return shaper.second;
}

// Decimates by the oversampling factor through a Blackman-windowed sinc
// low-pass: pass to 0.45 of the output rate, stop band about -74 dB from the
// output Nyquist rate on, so nothing above it folds back into the voice.
std::vector<double> decimate(const std::vector<double>& signal, std::size_t oversample) {
    const std::size_t taps = 64 * oversample + 1;
    std::vector<double> kernel(taps, 0.0);
    double kernel_sum = 0.0;
    for (std::size_t index = 0; index < taps; ++index) {
        double centre = static_cast<double>(index) - static_cast<double>(taps - 1) / 2.0;
        double cutoff = 0.47 / static_cast<double>(oversample);
        double sinc = centre == 0.0 ? 2.0 * cutoff : std::sin(2.0 * pi * cutoff * centre) / (pi * centre);
        double phase = 2.0 * pi * static_cast<double>(index) / static_cast<double>(taps - 1);
        double window = 0.42 - 0.5 * std::cos(phase) + 0.08 * std::cos(2.0 * phase);
        kernel[index] = sinc * window;
        kernel_sum += kernel[index];
    }
    std::vector<double> decimated{};
    for (std::size_t position = taps; position < signal.size(); position += oversample) {
        double value = 0.0;
        for (std::size_t index = 0; index < taps; ++index) {
            value += kernel[index] * signal[position - index];
        }
        decimated.push_back(value / kernel_sum);
    }
    return decimated;
}

bool write_normalized(const std::string& path, const std::vector<double>& signal) {
    double peak = 0.0;
    for (std::size_t index = 0; index < signal.size(); ++index) {
        peak = std::max(peak, std::fabs(signal[index]));
    }
    ducklips::AudioBuffer audio{};
    audio.channels = 2;
    audio.sample_rate = 48000.0;
    for (std::size_t index = 0; index < signal.size(); ++index) {
        float sample = static_cast<float>(0.8 * signal[index] / std::max(peak, 1.0e-30));
        audio.samples.push_back(sample);
        audio.samples.push_back(sample);
    }
    bool written = ducklips::write_wav_float(path, audio);
    return written;
}

}

int main(int argc, char** argv) {
    std::string first_file{};
    std::string second_file{};
    std::string output = "out/surface.wav";
    double seconds = 1.5;
    std::string flow_output{};
    double pitch = 110.0;
    double lung_pressure = 600.0;   // Pa, relaxed speech (800 sounded pressed)
    double rest_area = 4.0e-6;      // m^2, posterior chink
    double abducted_area = 30.0e-6; // m^2, folds apart before and after voicing
    double onset_seconds = 0.08;
    double turbulence_gain = 0.012; // breath noise (MEASURED against the architect's recording)
    double buffet_gain = 0.002;     // fold-surface buffeting, fraction of lung pressure
    double noise_corner = 2000.0;   // Hz, turbulence spectrum corner (0: unshaped white)
    double wobble_gain = 1.0;       // scale of all body unsteadiness
    double tremor_gain = 0.002;     // fast tension tremor, fraction of tension
    double adduction_gain = 0.08;   // adduction wander, fraction of the rest opening
    double unit_pitch = 136.0;      // Hz at tension 1 (MEASURED for these folds at 600 Pa)
    double medial_roundness = 1.0e-6; // m^2, smooth exchange of the narrowest point
    double contact_softness = 0.1e-6; // m^2, gradual sealing of the closing glottis
    std::size_t strips = 1;         // strips along the folds (more took energy out of 1-4 kHz)
    for (int index = 1; index + 1 < argc; index += 2) {
        std::string flag = argv[index];
        std::string value = argv[index + 1];
        if (flag == "--from") { first_file = value; }
        if (flag == "--to") { second_file = value; }
        if (flag == "--out") { output = value; }
        if (flag == "--seconds") { seconds = std::atof(value.c_str()); }
        if (flag == "--pitch") { pitch = std::atof(value.c_str()); }
        if (flag == "--flow-out") { flow_output = value; }
        if (flag == "--lungs") { lung_pressure = std::atof(value.c_str()); }
        if (flag == "--rest") { rest_area = std::atof(value.c_str()) * 1.0e-6; }
        if (flag == "--turbulence") { turbulence_gain = std::atof(value.c_str()); }
        if (flag == "--tremor") { tremor_gain = std::atof(value.c_str()); }
        if (flag == "--adduction") { adduction_gain = std::atof(value.c_str()); }
        if (flag == "--noise-corner") { noise_corner = std::atof(value.c_str()); }
        if (flag == "--wobble") { wobble_gain = std::atof(value.c_str()); }
        if (flag == "--buffet") { buffet_gain = std::atof(value.c_str()); }
        if (flag == "--unit-pitch") { unit_pitch = std::atof(value.c_str()); }
        if (flag == "--softness") { contact_softness = std::atof(value.c_str()) * 1.0e-6; }
        if (flag == "--roundness") { medial_roundness = std::atof(value.c_str()) * 1.0e-6; }
        if (flag == "--onset") { onset_seconds = std::atof(value.c_str()); }
        if (flag == "--strips") { strips = static_cast<std::size_t>(std::atoi(value.c_str())); }
    }
    Profile first = read_profile(first_file);
    Profile second = second_file.empty() ? first : read_profile(second_file);
    const std::size_t stations = 50;
    const std::size_t oversample = 4;
    const double rate = 48000.0 * static_cast<double>(oversample);
    std::vector<double> first_areas{};
    std::vector<double> second_areas{};
    resample(first, stations, first_areas);
    resample(second, stations, second_areas);
    ducklips::CavityGeometry geometry = ducklips::single_tube(std::span<const double>(first_areas), total_length(first), 6);
    ducklips::AirConstants air{};
    ducklips::SurfaceNetwork network{};
    network.build(geometry, rate, air);

    // The vocal folds: the two-mass folds, driven by the pressure difference
    // between the lungs (a silent reservoir) and the surface just above them,
    // with the flow solved against the surface's own response every sample.
    ducklips::FoldParameters folds{};
    folds.strips = strips;
    folds.contact_softness = contact_softness;
    folds.medial_roundness = medial_roundness;
    ducklips::TwoMassGlottis glottis{};
    glottis.reset(folds);
    ducklips::AirProperties fold_air{};
    ducklips::NoiseGenerator surface_noise{};
    surface_noise.seed(7);
    ducklips::TurbulenceNoise jet_noise{};
    jet_noise.prepare(11, 300.0, rate);
    NoiseShaper jet_shaper{};
    prepare_shaper(jet_shaper, noise_corner, rate);
    // Body wobble: slow tension wander, a little flutter, slow pressure wander.
    ducklips::SlowNoise tension_wander{};
    tension_wander.prepare(13, 6.0, rate);
    ducklips::SlowNoise flutter{};
    flutter.prepare(17, 12.0, rate);
    ducklips::SlowNoise pressure_wander{};
    pressure_wander.prepare(19, 4.0, rate);
    // Faster physiological unsteadiness, at the scale of single cycles:
    // muscle-tension tremor (pitch jitter) and arytenoid adduction wander
    // (cycle amplitude and closure, shimmer).
    ducklips::SlowNoise tremor{};
    tremor.prepare(23, 40.0, rate);
    ducklips::SlowNoise adduction{};
    adduction.prepare(29, 20.0, rate);
    const double critical_reynolds = 1800.0;

    const std::size_t total = static_cast<std::size_t>(seconds * rate);
    std::vector<double> heard(total, 0.0);
    std::vector<double> flows(total, 0.0);
    std::vector<double> areas(stations, 0.0);
    double transglottal = 0.0;
    double flow = 0.0;
    for (std::size_t sample = 0; sample < total; ++sample) {
        if (sample % 192 == 0) {
            double t = static_cast<double>(sample) / static_cast<double>(total);
            double ease = t * t * (3.0 - 2.0 * t);
            for (std::size_t station = 0; station < stations; ++station) {
                areas[station] = first_areas[station] + ease * (second_areas[station] - first_areas[station]);
            }
            network.reshape(areas, 192);
        }
        double time = static_cast<double>(sample) / rate;
        double breath = std::min(1.0, time / 0.05) * std::min(1.0, std::max(0.0, (seconds - time) / 0.05));
        double lungs = lung_pressure * breath * (1.0 + wobble_gain * 0.015 * pressure_wander.next());
        double wobble = 1.0 + wobble_gain * (0.004 * tension_wander.next() + 0.002 * flutter.next() + tremor_gain * tremor.next());
        // Breathy onset and offset: the folds start apart, close in while the
        // breath begins, and part again at the end (smooth, about 80 ms).
        double closing_in = std::clamp(time / onset_seconds, 0.0, 1.0);
        double parting = std::clamp((seconds - time) / onset_seconds, 0.0, 1.0);
        double adducted = closing_in * closing_in * (3.0 - 2.0 * closing_in) * parting * parting * (3.0 - 2.0 * parting);
        double posture = rest_area + (1.0 - adducted) * (abducted_area - rest_area);
        double opening = posture * std::max(0.0, 1.0 + wobble_gain * adduction_gain * adduction.next());
        double tension = pitch / unit_pitch * wobble;

        // Turbulence of the glottal jet: slight, in bursts with each opening.
        double area = std::max(glottis.minimum_area(), 0.0);
        double jet = area > 1.0e-10 ? std::fabs(flow) / area : 0.0;
        double slit = 2.0 * area / folds.length;
        double reynolds = jet * slit / (air.density > 0.0 ? 1.6e-5 : 1.6e-5);
        double gate = reynolds > critical_reynolds ? 1.0 - (critical_reynolds / reynolds) * (critical_reynolds / reynolds) : 0.0;
        double turbulence = turbulence_gain * 0.5 * air.density * jet * jet * gate * shape_noise(jet_shaper, jet_noise.next());
        double surface_buffet = buffet_gain * lungs * surface_noise.gaussian();

        glottis.advance(transglottal, opening, tension, 1.0 / rate, 1, surface_buffet + turbulence);
        area = std::max(glottis.minimum_area(), 0.0);
        network.gather();
        double above = network.glottis_free_pressure();
        double impedance = network.glottis_impedance();
        flow = ducklips::orifice_flow(area, lungs - above, impedance, fold_air);
        transglottal = lungs - (above + impedance * flow);
        heard[sample] = network.scatter(flow, turbulence, 0.0);
        flows[sample] = flow;
    }
    if (!write_normalized(output, decimate(heard, oversample))) {
        return 1;
    }
    if (!flow_output.empty() && !write_normalized(flow_output, decimate(flows, oversample))) {
        return 1;
    }
    std::printf("wrote %s (%zu pieces, %zu paths)\n", output.c_str(), network.piece_count(), network.link_count());
    return 0;
}
