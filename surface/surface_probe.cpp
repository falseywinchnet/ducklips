// Plucks a cavity of surface pieces with a flow impulse at the glottis and
// reports where it rings (peaks of the radiated spectrum).
#include "surface_network.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.14159265358979323846;

void spectrum(const std::vector<double>& signal, double rate, double step_hz, double top_hz, std::vector<double>& magnitude) {
    std::size_t bins = static_cast<std::size_t>(top_hz / step_hz);
    magnitude.assign(bins, 0.0);
    for (std::size_t bin = 1; bin < bins; ++bin) {
        double frequency = static_cast<double>(bin) * step_hz;
        std::complex<double> sum(0.0, 0.0);
        std::complex<double> rotation = std::polar(1.0, -2.0 * pi * frequency / rate);
        std::complex<double> phase(1.0, 0.0);
        for (std::size_t index = 0; index < signal.size(); ++index) {
            sum += signal[index] * phase;
            phase *= rotation;
        }
        magnitude[bin] = std::abs(sum);
    }
}

void read_profile(const std::string& path, std::size_t stations, ducklips::CavityGeometry& geometry) {
    std::vector<double> areas{};
    std::vector<double> lengths{};
    std::FILE* file = std::fopen(path.c_str(), "r");
    double area = 0.0;
    double length = 0.0;
    while (file != nullptr && std::fscanf(file, "%lf %lf", &area, &length) == 2) {
        areas.push_back(area * 1.0e-4);
        lengths.push_back(length);
    }
    if (file != nullptr) {
        std::fclose(file);
    }
    double total = 0.0;
    for (std::size_t index = 0; index < lengths.size(); ++index) {
        total += lengths[index];
    }
    if (areas.empty()) {
        return;
    }
    geometry.tubes[0].length = total;
    for (std::size_t station = 0; station < stations; ++station) {
        double begin = total * static_cast<double>(station) / static_cast<double>(stations);
        double end = total * static_cast<double>(station + 1) / static_cast<double>(stations);
        double running = 0.0;
        double covered = 0.0;
        double position = 0.0;
        for (std::size_t index = 0; index < areas.size(); ++index) {
            double low = std::max(begin, position);
            double high = std::min(end, position + lengths[index]);
            if (high > low) {
                running += areas[index] * (high - low);
                covered += high - low;
            }
            position += lengths[index];
        }
        geometry.station_area[station] = running / std::max(covered, 1.0e-12);
    }
}

}

int main(int argc, char** argv) {
    std::size_t stations = 50;
    std::size_t ring = 6;
    double oversample = 4.0;
    double area_cm2 = 3.0;
    bool rigid = false;
    double perimeter = 1.5;
    std::string area_file{};
    bool solid = false;
    bool bent = false;
    bool cross = false;
    double share = 1.0;
    double tension = 0.0;
    double width = 0.0;
    double length = 0.175;
    double top_hz = 6000.0;
    bool kick = false;
    double timing = 1.0;
    std::size_t beside = 0;      // stations of a parallel tube beside the middle of the airway
    double beside_area = 1.0;    // cm^2
    double opening_cm2 = 0.5;    // cm^2 at every station of it
    for (int index = 1; index + 1 < argc; index += 2) {
        std::string flag = argv[index];
        std::string value = argv[index + 1];
        if (flag == "--stations") { stations = static_cast<std::size_t>(std::atoi(value.c_str())); }
        if (flag == "--ring") { ring = static_cast<std::size_t>(std::atoi(value.c_str())); }
        if (flag == "--oversample") { oversample = std::atof(value.c_str()); }
        if (flag == "--area") { area_cm2 = std::atof(value.c_str()); }
        if (flag == "--rigid") { rigid = value == "1"; }
        if (flag == "--areas") { area_file = value; }
        if (flag == "--perimeter") { perimeter = std::atof(value.c_str()); }
        if (flag == "--solid") { solid = value == "1"; }
        if (flag == "--bent") { bent = value == "1"; }
        if (flag == "--cross") { cross = value == "1"; }
        if (flag == "--share") { share = std::atof(value.c_str()); }
        if (flag == "--tension") { tension = std::atof(value.c_str()); }
        if (flag == "--width") { width = std::atof(value.c_str()); }
        if (flag == "--length") { length = std::atof(value.c_str()); }
        if (flag == "--top") { top_hz = std::atof(value.c_str()); }
        if (flag == "--kick") { kick = value == "1"; }
        if (flag == "--timing") { timing = std::atof(value.c_str()); }
        if (flag == "--beside") { beside = static_cast<std::size_t>(std::atoi(value.c_str())); }
        if (flag == "--beside-area") { beside_area = std::atof(value.c_str()); }
        if (flag == "--opening") { opening_cm2 = std::atof(value.c_str()); }
    }
    const double rate = 48000.0 * oversample;
    std::vector<double> uniform(stations, area_cm2 * 1.0e-4);
    ducklips::CavityGeometry geometry = ducklips::single_tube(std::span<const double>(uniform), length, ring);
    geometry.tubes[0].perimeter_factor = perimeter;
    geometry.tubes[0].solid = solid;
    geometry.tubes[0].bent = bent;
    geometry.tubes[0].cross_paths = cross;
    geometry.tubes[0].axial_share = share;
    geometry.tubes[0].edge_tension = tension;
    geometry.tubes[0].width = width;
    geometry.tubes[0].cross_timing = timing;
    if (solid) {
        geometry.tubes[0].perimeter_factor = 1.0;
    }
    if (!area_file.empty()) {
        read_profile(area_file, stations, geometry);
    }
    if (rigid) {
        geometry.tubes[0].yielding = false;
    }
    if (beside > 0) {
        ducklips::TubeGeometry under{};
        under.stations = beside;
        under.length = geometry.tubes[0].length * static_cast<double>(beside) / static_cast<double>(stations);
        under.ring_pieces = 1;
        under.parent = 0;
        under.parent_station = stations / 2;
        under.parallel = true;
        under.end = ducklips::TubeEnd::closed;
        geometry.tubes.push_back(under);
        for (std::size_t index = 0; index < beside; ++index) {
            geometry.station_area.push_back(beside_area * 1.0e-4);
        }
    }
    ducklips::SurfaceNetwork network{};
    network.build(geometry, rate, ducklips::AirConstants{});
    if (beside > 0) {
        std::vector<double> openings(geometry.station_area.size(), 0.0);
        for (std::size_t index = stations; index < openings.size(); ++index) {
            openings[index] = opening_cm2 * 1.0e-4;
        }
        network.set_openings(std::span<const double>(openings));
        network.reshape(std::span<const double>(geometry.station_area), 1);
    }
    std::vector<double> response(static_cast<std::size_t>(0.25 * rate), 0.0);
    // --kick 1: a flow impulse into one piece of the middle ring, heard as
    // the pressure on the piece opposite (cross modes); otherwise the
    // glottis is plucked and the lips are heard.
    std::size_t middle = network.first_piece(stations / 2);
    std::size_t opposite = middle + ring / 2;
    for (std::size_t index = 0; index < response.size(); ++index) {
        network.gather();
        double pluck = (!kick && index == 0) ? 1.0e-3 : 0.0;
        if (kick && index == 0) {
            network.add_flow(middle, 1.0e-3);
        }
        double heard = network.scatter(pluck, 0.0, 0.0);
        response[index] = kick ? network.piece_pressure(opposite) - network.piece_pressure(middle + (ring / 4)) : heard;
    }
    double tail = 0.0;
    double whole = 0.0;
    for (std::size_t index = 0; index < response.size(); ++index) {
        whole = std::max(whole, std::fabs(response[index]));
        if (index >= response.size() - response.size() / 10) {
            tail = std::max(tail, std::fabs(response[index]));
        }
    }
    std::printf("  peak %.3e, last-tenth peak %.3e (%s)\n", whole, tail, tail < 1.0e-2 * whole ? "decays" : "DOES NOT DECAY");
    // Travel time: first sample where the radiated signal exceeds 10 % of its
    // first-pass peak (the pluck's arrival at the lips).
    double first_peak = 0.0;
    for (std::size_t index = 0; index < response.size() && index < static_cast<std::size_t>(0.0012 * rate); ++index) {
        first_peak = std::max(first_peak, std::fabs(response[index]));
    }
    std::size_t arrival = 0;
    while (arrival < response.size() && std::fabs(response[arrival]) < 0.1 * first_peak) {
        arrival += 1;
    }
    std::printf("  first arrival at the lips: %.1f samples = %.3f ms (length/c = %.3f ms)\n", static_cast<double>(arrival),
                1000.0 * static_cast<double>(arrival) / rate, 1000.0 * geometry.tubes[0].length / 350.0);
    std::vector<double> magnitude{};
    spectrum(response, rate, 10.0, top_hz, magnitude);
    std::printf("stations %zu ring %zu x%.0f %s perimeter %.2f: pieces %zu paths %zu\n", stations, ring,
                oversample, rigid ? "rigid" : "yielding", perimeter, network.piece_count(), network.link_count());
    std::printf("  peaks Hz(dB,width):");
    double top = *std::max_element(magnitude.begin(), magnitude.end());
    for (std::size_t bin = 2; bin + 2 < magnitude.size(); ++bin) {
        if (magnitude[bin] > magnitude[bin - 1] && magnitude[bin] >= magnitude[bin + 1] && magnitude[bin] > 0.03 * top) {
            std::size_t left = bin;
            while (left > 0 && magnitude[left] > magnitude[bin] / std::sqrt(2.0)) { left -= 1; }
            std::size_t right = bin;
            while (right + 1 < magnitude.size() && magnitude[right] > magnitude[bin] / std::sqrt(2.0)) { right += 1; }
            std::printf(" %zu(%.0f,%zu)", bin * 10, 20.0 * std::log10(magnitude[bin] / top), (right - left) * 10);
        }
    }
    std::printf("\n");
    return 0;
}
