// Boundary-surface sibilant prototype: experiments and measurements.
//
//   boundary_lab run OUT_DIR [--oversampling N] [--patch-cm X]
//   boundary_lab analyze FILE.wav FROM_S TO_S
#include "analysis.hpp"
#include "geometry.hpp"
#include "transport.hpp"
#include "tube.hpp"
#include "ducklips/wav.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

using boundary::BoundaryTransport;
using boundary::Decimator;
using boundary::Peak;
using boundary::Scene;
using boundary::SpectrumSummary;
using boundary::TransportSettings;
using boundary::TubeControl;

constexpr double audio_rate = 48000.0;
constexpr std::size_t not_found = static_cast<std::size_t>(-1);

struct RunOptions final {
    std::string out_dir{"out"};
    std::size_t oversampling{4};
    double patch_cm{0.4};
    double face_extent{0.0}; // m; > 0 adds the face (or the duct flange) as a baffle
    bool visibility_culling{true};
    double face_patch_scale{2.0};
    double tooth_cm{0.2};
};

double seconds_since(std::chrono::steady_clock::time_point started) {
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - started;
    double result = elapsed.count();
    return result;
}

// Runs the boundary model over an internal-rate source and fills audio at 48 kHz.
void run_boundary(BoundaryTransport& model, const std::vector<double>& source, std::size_t oversampling,
                  std::vector<double>& audio) {
    Decimator decimator{};
    decimator.prepare(oversampling, model.internal_rate(), 20000.0);
    audio.assign(source.size() / oversampling, 0.0);
    model.reset();
    std::size_t written = 0;
    double sample = 0.0;
    for (std::size_t i = 0; i < source.size(); ++i) {
        double microphone = model.step(source[i]);
        if (decimator.push(microphone, sample) && written < audio.size()) {
            audio[written] = sample;
            written += 1;
        }
    }
}

void run_tube(TubeControl& model, double internal_rate, const std::vector<double>& source, std::size_t oversampling,
              std::vector<double>& audio) {
    Decimator decimator{};
    decimator.prepare(oversampling, internal_rate, 20000.0);
    audio.assign(source.size() / oversampling, 0.0);
    model.reset();
    std::size_t written = 0;
    double sample = 0.0;
    for (std::size_t i = 0; i < source.size(); ++i) {
        double microphone = model.step(source[i]);
        if (decimator.push(microphone, sample) && written < audio.size()) {
            audio[written] = sample;
            written += 1;
        }
    }
}

void print_peaks(const char* label, const std::vector<Peak>& peaks) {
    std::printf("  %-28s", label);
    std::size_t shown = std::min<std::size_t>(peaks.size(), 6);
    for (std::size_t k = 0; k < shown; ++k) {
        std::printf(" %6.0f Hz (%+5.1f dB, B %4.0f)", peaks[k].hz, peaks[k].db, peaks[k].bandwidth_hz);
    }
    std::printf("\n");
}

void print_summary(const char* label, const SpectrumSummary& s, double real_time_factor) {
    std::printf("  %-34s rms %7.2f dB  peak %6.0f Hz  centroid %6.0f Hz  2.5-4k %7.2f  4-8k %7.2f  RTF %7.3f\n", label,
                s.rms_db, s.peak_hz, s.centroid_hz, s.band_2500_4000_db, s.band_4000_8000_db, real_time_factor);
}

void print_statistics(const Scene& scene, const boundary::CompileStatistics& s) {
    std::printf("  scene %-4s surfaces %zu patches %zu linked pairs %zu taps p/a/mic %zu/%zu/%zu ring %zu\n",
                scene.name.c_str(), scene.surfaces.size(), s.patches, s.linked_pairs, s.pressure_taps,
                s.acceleration_taps, s.microphone_taps, s.ring_length);
    std::printf("    clamped subpoints %zu occluded %zu largest static gain %.3f membrane modes %zu (highest %.0f Hz)"
                " compile %.2f s\n",
                s.clamped_subpoints, s.occluded_subpoints, s.largest_static_gain, s.membrane_modes,
                s.highest_membrane_hz, s.compile_seconds);
}

std::size_t find_surface(const Scene& scene, const std::string& name) {
    for (std::size_t s = 0; s < scene.surfaces.size(); ++s) {
        if (scene.surfaces[s].name == name) {
            return s;
        }
    }
    return not_found;
}

std::size_t middle_patch(const Scene& scene, std::size_t surface_index) {
    const boundary::Surface& surface = scene.surfaces[surface_index];
    std::size_t u = surface.cells_u / 2;
    std::size_t v = surface.cells_v / 2;
    std::size_t result = surface.first_patch + v * surface.cells_u + u;
    return result;
}

bool finite_signal(const std::vector<double>& signal) {
    for (std::size_t i = 0; i < signal.size(); ++i) {
        if (!std::isfinite(signal[i])) {
            return false;
        }
    }
    return true;
}

// Writes [gap, a, gap, b, gap] as stereo float WAV, both scaled by one gain
// so their relative level is preserved.
bool write_pair(const std::string& path, const std::vector<double>& a, const std::vector<double>& b) {
    double peak = 1.0e-30;
    for (std::size_t i = 0; i < a.size(); ++i) {
        peak = std::max(peak, std::fabs(a[i]));
    }
    for (std::size_t i = 0; i < b.size(); ++i) {
        peak = std::max(peak, std::fabs(b[i]));
    }
    double gain = 0.5 / peak;
    std::size_t gap = static_cast<std::size_t>(0.15 * audio_rate);
    std::size_t frames = 3 * gap + a.size() + b.size();
    ducklips::AudioBuffer buffer{};
    buffer.channels = 2;
    buffer.sample_rate = audio_rate;
    buffer.samples.assign(2 * frames, 0.0f);
    std::size_t position = gap;
    for (std::size_t i = 0; i < a.size(); ++i) {
        float value = static_cast<float>(gain * a[i]);
        buffer.samples[2 * (position + i)] = value;
        buffer.samples[2 * (position + i) + 1] = value;
    }
    position = 2 * gap + a.size();
    for (std::size_t i = 0; i < b.size(); ++i) {
        float value = static_cast<float>(gain * b[i]);
        buffer.samples[2 * (position + i)] = value;
        buffer.samples[2 * (position + i) + 1] = value;
    }
    bool written = ducklips::write_wav_float(path, buffer);
    return written;
}

struct PhoneResult final {
    std::vector<double> boundary_audio{};
    std::vector<double> boundary_rigid_audio{};
    std::vector<double> tube_audio{};
    SpectrumSummary boundary{};
    SpectrumSummary boundary_rigid{};
    SpectrumSummary tube{};
    double boundary_rtf{0.0};
    double boundary_rigid_rtf{0.0};
    double tube_rtf{0.0};
};

void duct_validation(const RunOptions& options) {
    std::printf("\n[1] Closed-open rigid duct, L = 2.0 cm, 1 x 1 cm, piston at the closed end\n");
    Scene scene{};
    boundary::build_duct_scene(0.02, 0.01, 0.0025, false, options.face_extent, scene);
    TransportSettings settings{};
    settings.oversampling = options.oversampling;
    settings.membranes = false;
    settings.visibility_culling = options.visibility_culling;
    BoundaryTransport model{};
    model.compile(scene, settings);
    print_statistics(scene, model.statistics());
    std::vector<double> source{};
    boundary::make_band_impulse(model.internal_rate(), static_cast<std::size_t>(0.1 * model.internal_rate()), source);
    std::vector<double> audio{};
    run_boundary(model, source, options.oversampling, audio);
    std::vector<Peak> peaks{};
    boundary::transfer_peaks(audio, audio_rate, 1000.0, 14000.0, 30.0, 60.0, peaks);
    print_peaks("boundary transport", peaks);
    TubeControl tube{};
    tube.compile(scene, model.internal_rate(), settings.air);
    std::vector<double> tube_audio{};
    run_tube(tube, model.internal_rate(), source, options.oversampling, tube_audio);
    boundary::transfer_peaks(tube_audio, audio_rate, 1000.0, 14000.0, 30.0, 60.0, peaks);
    print_peaks("1-D tube control", peaks);
    double radius = std::sqrt(1.0e-4 / ducklips::pi);
    double unflanged = settings.air.sound_speed / (4.0 * (0.02 + 0.61 * radius));
    double flanged = settings.air.sound_speed / (4.0 * (0.02 + 0.85 * radius));
    std::printf("  reference quarter-wave: %.0f Hz (unflanged end correction), %.0f Hz (flanged); 3rd harmonic x3\n",
                unflanged, flanged);
}

void reciprocity(const RunOptions& options, const Scene& scene) {
    std::printf("\n[3] Reciprocity in /%s/ (rigid): wall acceleration at A, pressure at B, and swapped\n",
                scene.name.c_str());
    std::size_t surface_a = find_surface(scene, "tongue_blade");
    std::size_t surface_b = find_surface(scene, "lower_lip");
    if (surface_a == not_found || surface_b == not_found) {
        std::printf("  surfaces missing\n");
        return;
    }
    std::size_t patch_a = middle_patch(scene, surface_a);
    std::size_t patch_b = middle_patch(scene, surface_b);
    TransportSettings settings{};
    settings.oversampling = options.oversampling;
    settings.visibility_culling = options.visibility_culling;
    settings.membranes = false;
    BoundaryTransport model{};
    model.compile(scene, settings);
    double rate = model.internal_rate();
    std::size_t count = 32768;
    std::vector<double> source{};
    boundary::make_band_impulse(rate, count, source);
    std::vector<double> at_b(count, 0.0);
    std::vector<double> at_a(count, 0.0);
    model.inject_acceleration_at(patch_a);
    model.reset();
    for (std::size_t i = 0; i < count; ++i) {
        model.step(source[i]);
        at_b[i] = model.patch_pressure(patch_b);
    }
    model.inject_acceleration_at(patch_b);
    model.reset();
    for (std::size_t i = 0; i < count; ++i) {
        model.step(source[i]);
        at_a[i] = model.patch_pressure(patch_a);
    }
    double area_a = scene.patches[patch_a].area;
    double area_b = scene.patches[patch_b].area;
    std::vector<std::complex<double>> forward(count);
    std::vector<std::complex<double>> backward(count);
    for (std::size_t i = 0; i < count; ++i) {
        forward[i] = std::complex<double>(at_b[i] / area_a, 0.0);
        backward[i] = std::complex<double>(at_a[i] / area_b, 0.0);
    }
    boundary::fft_in_place(forward);
    boundary::fft_in_place(backward);
    double bin_hz = rate / static_cast<double>(count);
    double magnitude_error = 0.0;
    double phase_error = 0.0;
    double worst_magnitude = 0.0;
    std::size_t bins = 0;
    for (std::size_t k = 0; k < count / 2; ++k) {
        double hz = static_cast<double>(k) * bin_hz;
        if (hz < 1000.0 || hz > 10000.0) {
            continue;
        }
        std::complex<double> ratio = forward[k] / backward[k];
        double db = 20.0 * std::log10(std::abs(ratio));
        double degrees = std::arg(ratio) * 180.0 / ducklips::pi;
        magnitude_error += std::fabs(db);
        phase_error += std::fabs(degrees);
        worst_magnitude = std::max(worst_magnitude, std::fabs(db));
        bins += 1;
    }
    std::printf("  A = %s centre, B = %s centre: mean |H_AB/H_BA| %.2f dB (worst %.2f), mean |phase| %.1f deg, 1-10 kHz\n",
                scene.surfaces[surface_a].name.c_str(), scene.surfaces[surface_b].name.c_str(),
                magnitude_error / static_cast<double>(bins), worst_magnitude, phase_error / static_cast<double>(bins));
}

void sibilant(const RunOptions& options, const Scene& scene, PhoneResult& result) {
    std::printf("\n[2] /%s/ geometry\n", scene.name.c_str());
    TransportSettings settings{};
    settings.oversampling = options.oversampling;
    settings.visibility_culling = options.visibility_culling;
    settings.membranes = true;
    BoundaryTransport model{};
    model.compile(scene, settings);
    print_statistics(scene, model.statistics());
    double rate = model.internal_rate();

    // Transfer function from the jet force to the microphone.
    std::vector<double> impulse{};
    boundary::make_band_impulse(rate, static_cast<std::size_t>(1.0 * rate), impulse);
    std::vector<double> response{};
    run_boundary(model, impulse, options.oversampling, response);
    std::size_t early_end = static_cast<std::size_t>(0.05 * audio_rate);
    std::size_t late_begin = static_cast<std::size_t>(0.8 * audio_rate);
    double early = boundary::rms(response, 0, early_end);
    double late = boundary::rms(response, late_begin, response.size());
    std::printf("  stability (impulse, 1 s): finite %s, late/early rms %.1f dB\n",
                finite_signal(response) ? "yes" : "NO", 20.0 * std::log10(std::max(late, 1.0e-300) / early));
    std::vector<double> head(response.begin(), response.begin() + static_cast<std::ptrdiff_t>(early_end * 2));
    std::vector<Peak> peaks{};
    boundary::transfer_peaks(head, audio_rate, 1000.0, 12000.0, 25.0, 100.0, peaks);
    print_peaks("boundary (membranes) peaks", peaks);
    TubeControl tube{};
    tube.compile(scene, rate, settings.air);
    std::vector<double> tube_response{};
    std::vector<double> short_impulse(impulse.begin(), impulse.begin() + static_cast<std::ptrdiff_t>(0.1 * rate));
    run_tube(tube, rate, short_impulse, options.oversampling, tube_response);
    boundary::transfer_peaks(tube_response, audio_rate, 1000.0, 12000.0, 25.0, 100.0, peaks);
    print_peaks("1-D tube control peaks", peaks);

    // The same jet noise drives every model.
    const double seconds = 0.5;
    std::vector<double> jet{};
    boundary::make_jet_noise(20261008, rate, static_cast<std::size_t>(seconds * rate), jet);
    std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    run_boundary(model, jet, options.oversampling, result.boundary_audio);
    result.boundary_rtf = seconds_since(started) / seconds;

    TransportSettings rigid_settings = settings;
    rigid_settings.membranes = false;
    BoundaryTransport rigid{};
    rigid.compile(scene, rigid_settings);
    started = std::chrono::steady_clock::now();
    run_boundary(rigid, jet, options.oversampling, result.boundary_rigid_audio);
    result.boundary_rigid_rtf = seconds_since(started) / seconds;

    started = std::chrono::steady_clock::now();
    run_tube(tube, rate, jet, options.oversampling, result.tube_audio);
    result.tube_rtf = seconds_since(started) / seconds;

    result.boundary = boundary::summarize(result.boundary_audio, audio_rate, 0.05, seconds);
    result.boundary_rigid = boundary::summarize(result.boundary_rigid_audio, audio_rate, 0.05, seconds);
    result.tube = boundary::summarize(result.tube_audio, audio_rate, 0.05, seconds);
    print_summary("boundary, membranes", result.boundary, result.boundary_rtf);
    print_summary("boundary, rigid walls", result.boundary_rigid, result.boundary_rigid_rtf);
    print_summary("1-D tube control", result.tube, result.tube_rtf);
    double difference = boundary::rms(result.boundary_audio, 0, result.boundary_audio.size());
    std::vector<double> delta(result.boundary_audio.size(), 0.0);
    for (std::size_t i = 0; i < delta.size(); ++i) {
        delta[i] = result.boundary_audio[i] - result.boundary_rigid_audio[i];
    }
    double residual = boundary::rms(delta, 0, delta.size());
    std::printf("  membrane contribution: (membranes - rigid) rms %.1f dB re signal\n",
                20.0 * std::log10(std::max(residual, 1.0e-300) / difference));
}

// Open duct: resonance frequency and bandwidth at the microphone against
// patch size and oversampling.
void duct_convergence(double patch_mm, std::size_t oversampling, double flange_cm, bool culling, bool averaging) {
    Scene scene{};
    boundary::build_duct_scene(0.02, 0.01, patch_mm * 1.0e-3, false, flange_cm * 0.01, scene);
    TransportSettings settings{};
    settings.oversampling = oversampling;
    settings.membranes = false;
    settings.visibility_culling = culling;
    settings.temporal_averaging = averaging;
    BoundaryTransport model{};
    model.compile(scene, settings);
    double rate = model.internal_rate();
    std::vector<double> source{};
    boundary::make_band_impulse(rate, static_cast<std::size_t>(0.1 * rate), source);
    std::vector<double> audio{};
    std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    run_boundary(model, source, oversampling, audio);
    double rtf = seconds_since(started) / 0.1;
    std::vector<Peak> peaks{};
    boundary::transfer_peaks(audio, audio_rate, 1000.0, 14000.0, 30.0, 60.0, peaks);
    const boundary::CompileStatistics& s = model.statistics();
    std::printf("  avg %d duct  %-7s flange %3.1f cm  patch %5.2f mm  OS %zu  patches %4zu  taps %8zu  static %.3f  RTF %8.2f  finite %s ",
                averaging ? 1 : 0, culling ? "culled" : "full", flange_cm, patch_mm, oversampling, s.patches, s.pressure_taps + s.acceleration_taps, s.largest_static_gain, rtf,
                finite_signal(audio) ? "yes" : "NO");
    std::size_t shown = std::min<std::size_t>(peaks.size(), 3);
    for (std::size_t k = 0; k < shown; ++k) {
        std::printf(" %6.0f Hz (B %5.0f, %+5.1f dB)", peaks[k].hz, peaks[k].bandwidth_hz, peaks[k].db);
    }
    std::printf("\n");
}

// Closed rigid box, reflectance 1: lossless in the continuum. Axial modes at
// c / (2 L) multiples; the DC mode is exactly neutral, so any quadrature
// excess shows up as growth.
void box_check(double patch_mm, std::size_t oversampling, double reflectance, bool averaging) {
    Scene scene{};
    boundary::build_duct_scene(0.04, 0.01, patch_mm * 1.0e-3, true, 0.0, scene);
    TransportSettings settings{};
    settings.oversampling = oversampling;
    settings.membranes = false;
    settings.reflectance_override = reflectance;
    settings.temporal_averaging = averaging;
    BoundaryTransport model{};
    model.compile(scene, settings);
    double rate = model.internal_rate();
    std::size_t count = static_cast<std::size_t>(0.15 * rate);
    std::vector<double> source{};
    boundary::make_band_impulse(rate, count, source);
    std::size_t probe = middle_patch(scene, find_surface(scene, "end_cap"));
    std::vector<double> response(count, 0.0);
    model.reset();
    for (std::size_t i = 0; i < count; ++i) {
        model.step(source[i]);
        response[i] = model.patch_pressure(probe);
    }
    // Remove the DC (static) part before looking for modes.
    double mean = 0.0;
    std::size_t tail_begin = count / 2;
    for (std::size_t i = tail_begin; i < count; ++i) {
        mean += response[i];
    }
    mean /= static_cast<double>(count - tail_begin);
    std::vector<double> centred(count, 0.0);
    for (std::size_t i = 0; i < count; ++i) {
        centred[i] = response[i] - mean;
    }
    double early = boundary::rms(centred, static_cast<std::size_t>(0.01 * rate), static_cast<std::size_t>(0.03 * rate));
    double late = boundary::rms(centred, static_cast<std::size_t>(0.12 * rate), count);
    std::vector<Peak> peaks{};
    boundary::transfer_peaks(centred, rate, 1000.0, 14000.0, 40.0, 20.0, peaks);
    const boundary::CompileStatistics& s = model.statistics();
    std::printf("  avg %d box   patch %5.2f mm  OS %zu  R %.3f  patches %4zu  static %.3f  finite %s  late/early %+7.1f dB  DC %.3g ",
                averaging ? 1 : 0, patch_mm, oversampling, reflectance, s.patches, s.largest_static_gain,
                finite_signal(response) ? "yes" : "NO", 20.0 * std::log10(std::max(late, 1.0e-300) / early), mean);
    std::size_t shown = std::min<std::size_t>(peaks.size(), 3);
    for (std::size_t k = 0; k < shown; ++k) {
        std::printf(" %6.0f Hz (B %5.0f)", peaks[k].hz, peaks[k].bandwidth_hz);
    }
    std::printf("\n");
}

// Impulse stability of one sibilant scene variant (rigid walls).
void stability_variant(const char* phone, double patch_cm, double face_cm, double face_scale, double tooth_cm,
                       bool culling, std::size_t oversampling, bool averaging) {
    boundary::SibilantPosture posture = boundary::s_posture();
    if (std::string(phone) == "sh") {
        posture = boundary::sh_posture();
    }
    posture.tooth_thickness = tooth_cm;
    Scene scene{};
    boundary::build_sibilant_scene(posture, patch_cm * 0.01, face_cm * 0.01, face_scale, scene);
    TransportSettings settings{};
    settings.oversampling = oversampling;
    settings.membranes = false;
    settings.visibility_culling = culling;
    settings.temporal_averaging = averaging;
    BoundaryTransport model{};
    model.compile(scene, settings);
    std::vector<double> impulse{};
    boundary::make_band_impulse(model.internal_rate(), static_cast<std::size_t>(0.3 * model.internal_rate()), impulse);
    std::vector<double> response{};
    run_boundary(model, impulse, oversampling, response);
    double early = boundary::rms(response, 0, static_cast<std::size_t>(0.03 * audio_rate));
    double late = boundary::rms(response, static_cast<std::size_t>(0.25 * audio_rate), response.size());
    std::vector<double> head(response.begin(), response.begin() + static_cast<std::ptrdiff_t>(0.1 * audio_rate));
    std::vector<Peak> peaks{};
    boundary::transfer_peaks(head, audio_rate, 1000.0, 12000.0, 25.0, 100.0, peaks);
    std::printf("  avg %d /%-2s/ %-6s patch %.2f face %.1f x%.1f tooth %.2f OS %zu  patches %3zu static %.3f finite %-3s late/early %+8.1f dB ",
                averaging ? 1 : 0, phone, culling ? "culled" : "full", patch_cm, face_cm, face_scale, tooth_cm, oversampling,
                model.statistics().patches, model.statistics().largest_static_gain,
                finite_signal(response) ? "yes" : "NO", 20.0 * std::log10(std::max(late, 1.0e-300) / early));
    std::size_t shown = std::min<std::size_t>(peaks.size(), 4);
    for (std::size_t k = 0; k < shown; ++k) {
        std::printf(" %5.0f Hz (B %4.0f, %+5.1f)", peaks[k].hz, peaks[k].bandwidth_hz, peaks[k].db);
    }
    std::printf("\n");
}

int stability_study() {
    stability_variant("s", 0.4, 2.0, 1.0, 0.2, false, 4, true);
    stability_variant("sh", 0.4, 2.0, 1.0, 0.2, false, 4, true);
    stability_variant("s", 0.4, 0.0, 1.0, 0.2, false, 4, true);
    stability_variant("s", 0.4, 2.0, 1.0, 0.2, true, 4, true);
    stability_variant("sh", 0.4, 2.0, 1.0, 0.2, true, 4, true);
    return 0;
}

int validate() {
    std::printf("Open duct L 2 cm, 1 x 1 cm (1-D control 3557 Hz; quarter-wave 3529-3733 Hz)\n");
    duct_convergence(2.5, 4, 0.0, true, false);
    duct_convergence(1.6667, 4, 0.0, true, false);
    duct_convergence(1.25, 4, 0.0, true, false);
    duct_convergence(2.5, 4, 2.0, true, false);
    duct_convergence(2.5, 4, 0.0, false, false);
    duct_convergence(2.5, 4, 2.0, false, false);
    duct_convergence(2.5, 4, 4.0, false, false);
    duct_convergence(1.6667, 4, 2.0, false, false);
    duct_convergence(2.5, 4, 2.0, false, true);
    duct_convergence(2.5, 8, 2.0, false, true);
    std::printf("Closed box 4 x 1 x 1 cm (axial modes 4375, 8750, 13125 Hz)\n");
    box_check(5.0, 4, 1.0, false);
    box_check(2.5, 4, 1.0, false);
    box_check(2.5, 8, 1.0, false);
    box_check(1.6667, 4, 1.0, false);
    box_check(5.0, 4, 1.0, true);
    box_check(2.5, 4, 1.0, true);
    box_check(2.5, 8, 1.0, true);
    return 0;
}

int run(const RunOptions& options) {
    std::printf("Boundary-surface sibilant prototype: oversampling %zu (internal %.0f Hz), patch %.2f cm, face %.1f cm,"
                " kernel %s\n",
                options.oversampling, audio_rate * static_cast<double>(options.oversampling), options.patch_cm,
                options.face_extent * 100.0, options.visibility_culling ? "culled" : "full");
    duct_validation(options);
    Scene s_scene{};
    Scene sh_scene{};
    boundary::build_sibilant_scene(boundary::s_posture(), options.patch_cm * 0.01, options.face_extent,
                                   options.face_patch_scale, s_scene);
    boundary::build_sibilant_scene(boundary::sh_posture(), options.patch_cm * 0.01, options.face_extent,
                                   options.face_patch_scale, sh_scene);
    PhoneResult s_result{};
    PhoneResult sh_result{};
    sibilant(options, s_scene, s_result);
    sibilant(options, sh_scene, sh_result);
    reciprocity(options, sh_scene);
    std::printf("\n[4] /sh/ - /s/ level (same jet source in every model)\n");
    std::printf("  boundary, membranes   %+6.2f dB\n", sh_result.boundary.rms_db - s_result.boundary.rms_db);
    std::printf("  boundary, rigid       %+6.2f dB\n", sh_result.boundary_rigid.rms_db - s_result.boundary_rigid.rms_db);
    std::printf("  1-D tube control      %+6.2f dB\n", sh_result.tube.rms_db - s_result.tube.rms_db);
    std::string suffix = "";
    if (options.face_extent > 0.0) {
        suffix = "_face";
    }
    if (!options.visibility_culling) {
        suffix += "_full";
    }
    std::string boundary_path = options.out_dir + "/boundary_s_sh" + suffix + ".wav";
    std::string tube_path = options.out_dir + "/tube_s_sh.wav";
    bool ok = write_pair(boundary_path, s_result.boundary_audio, sh_result.boundary_audio);
    ok = write_pair(tube_path, s_result.tube_audio, sh_result.tube_audio) && ok;
    std::printf("\nwrote %s and %s (%s)\n", boundary_path.c_str(), tube_path.c_str(), ok ? "ok" : "FAILED");
    return 0;
}

int analyze(const std::string& path, double from_s, double to_s) {
    ducklips::AudioBuffer buffer{};
    if (!ducklips::read_wav(path, buffer)) {
        std::printf("cannot read %s\n", path.c_str());
        return 1;
    }
    std::vector<double> mono(buffer.samples.size() / buffer.channels, 0.0);
    for (std::size_t i = 0; i < mono.size(); ++i) {
        mono[i] = static_cast<double>(buffer.samples[i * buffer.channels]);
    }
    SpectrumSummary summary = boundary::summarize(mono, buffer.sample_rate, from_s, to_s);
    print_summary(path.c_str(), summary, 0.0);
    return 0;
}

}

int main(int argc, char** argv) {
    if (argc >= 5 && std::string(argv[1]) == "analyze") {
        int status = analyze(argv[2], std::atof(argv[3]), std::atof(argv[4]));
        return status;
    }
    if (argc >= 2 && std::string(argv[1]) == "stability") {
        int status = stability_study();
        return status;
    }
    if (argc >= 2 && std::string(argv[1]) == "validate") {
        int status = validate();
        return status;
    }
    if (argc >= 3 && std::string(argv[1]) == "run") {
        RunOptions options{};
        options.out_dir = argv[2];
        for (int i = 3; i + 1 < argc; i += 2) {
            std::string flag = argv[i];
            if (flag == "--oversampling") {
                options.oversampling = static_cast<std::size_t>(std::atoi(argv[i + 1]));
            } else if (flag == "--full-kernel") {
                options.visibility_culling = std::atoi(argv[i + 1]) == 0;
            } else if (flag == "--face-scale") {
                options.face_patch_scale = std::atof(argv[i + 1]);
            } else if (flag == "--face-cm") {
                options.face_extent = std::atof(argv[i + 1]) * 0.01;
            } else if (flag == "--patch-cm") {
                options.patch_cm = std::atof(argv[i + 1]);
            }
        }
        int status = run(options);
        return status;
    }
    std::printf("usage: boundary_lab run OUT_DIR [--oversampling N] [--patch-cm X]\n"
                "       boundary_lab analyze FILE.wav FROM_S TO_S\n");
    return 2;
}
