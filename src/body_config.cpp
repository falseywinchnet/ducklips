#include "ducklips/body_config.hpp"
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <vector>

namespace ducklips {
namespace {
struct Field final { const char* name; double* value; double low; double high; };
std::string trim(const std::string& text) {
    std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) { return {}; }
    return text.substr(first, text.find_last_not_of(" \t\r\n") - first + 1);
}
}

bool read_body_config(std::istream& input, VoiceProfile& output, std::string& error) {
    VoiceProfile candidate = output;
    std::vector<Field> fields{
        {"tract_length", &candidate.tract_length, 0.1, 0.23},
        {"nasal_length", &candidate.nasal_length, 0.06, 0.17},
        {"port_position", &candidate.port_position, 0.1, 0.9},
        {"nostril_area", &candidate.nostril_area, 1e-06, 0.0008},
        {"tissue_transmission", &candidate.tissue_transmission, 0, 1},
        {"lung_pressure", &candidate.lung_pressure, 50, 2000},
        {"adducted_area", &candidate.adducted_area, -2e-06, 8e-05},
        {"abducted_area", &candidate.abducted_area, 1e-06, 0.0002},
        {"static_wall_compliance", &candidate.static_wall_compliance, 1e-11, 1e-06},
        {"base_pitch", &candidate.base_pitch, 60, 500},
        {"accent_rise", &candidate.accent_rise, 0, 1},
        {"declination", &candidate.declination, 0, 0.8},
        {"aspiration_gain", &candidate.aspiration_gain, 0, 0.05},
        {"frication_gain", &candidate.frication_gain, 0, 0.1},
        {"obstacle_gain", &candidate.obstacle_gain, 0, 10},
        {"jitter", &candidate.jitter, 0, 0.05},
        {"shimmer", &candidate.shimmer, 0, 0.2},
        {"tremor_rate", &candidate.tremor_rate, 0, 30},
        {"tremor_depth", &candidate.tremor_depth, 0, 0.1},
        {"breath_loudness", &candidate.breath_loudness, 0, 1},
        {"speaking_rate", &candidate.speaking_rate, 0.25, 4},
        {"folds.lower_mass", &candidate.folds.lower_mass, 1e-06, 0.001},
        {"folds.upper_mass", &candidate.folds.upper_mass, 1e-06, 0.001},
        {"folds.lower_stiffness", &candidate.folds.lower_stiffness, 1, 500},
        {"folds.upper_stiffness", &candidate.folds.upper_stiffness, 0.1, 200},
        {"folds.coupling_stiffness", &candidate.folds.coupling_stiffness, 0.1, 200},
        {"folds.length", &candidate.folds.length, 0.005, 0.025},
        {"folds.lower_thickness", &candidate.folds.lower_thickness, 0.0001, 0.006},
        {"folds.upper_thickness", &candidate.folds.upper_thickness, 0.0001, 0.006},
        {"folds.lower_damping_ratio", &candidate.folds.lower_damping_ratio, 0.01, 2},
        {"folds.upper_damping_ratio", &candidate.folds.upper_damping_ratio, 0.01, 2},
        {"folds.contact_softness", &candidate.folds.contact_softness, 0, 2e-06},
        {"folds.medial_roundness", &candidate.folds.medial_roundness, 0, 5e-06},
        {"folds.asymmetry", &candidate.folds.asymmetry, 0.5, 1.5},
        {"folds.roughness", &candidate.folds.roughness, 0, 0.2},
        {"speeds.tongue_tip", &candidate.speeds.tongue_tip, 5, 500},
        {"speeds.lips", &candidate.speeds.lips, 5, 500},
        {"speeds.tongue_body", &candidate.speeds.tongue_body, 5, 500},
        {"speeds.jaw", &candidate.speeds.jaw, 5, 500},
        {"speeds.velum", &candidate.speeds.velum, 5, 500},
        {"speeds.glottis", &candidate.speeds.glottis, 5, 500},
        {"speeds.lungs", &candidate.speeds.lungs, 5, 500},
        {"speeds.larynx_pitch", &candidate.speeds.larynx_pitch, 5, 500},
    };
    std::set<std::string> seen{};
    std::string line{};
    std::size_t number = 0;
    while (std::getline(input, line)) {
        ++number;
        line = trim(line.substr(0, line.find('#')));
        if (line.empty()) { continue; }
        std::size_t separator = line.find('=');
        std::string key = trim(line.substr(0, separator));
        std::string value = separator == std::string::npos ? "" : trim(line.substr(separator + 1));
        std::string reason{};
        if (separator == std::string::npos || value.empty()) {
            reason = "expected key=value";
        } else if (!seen.insert(key).second) {
            reason = "duplicate field " + key;
        } else if (key == "preset") {
            if (seen.size() != 1 || !voice_by_name(value, candidate)) {
                reason = "preset must be a known voice and the first setting";
            }
        } else if (key == "name") {
            if (value.size() > 80) { reason = "name exceeds 80 characters"; }
            else { candidate.name = value; }
        } else if (key == "folds.strips") {
            std::istringstream parser(value);
            std::size_t strips = 0;
            std::string extra{};
            if (!(parser >> strips) || (parser >> extra) || strips < 1 || strips > maximum_fold_strips) {
                reason = "folds.strips must be an integer in [1," + std::to_string(maximum_fold_strips) + "]";
            } else { candidate.folds.strips = strips; }
        } else if (key == "phonation" || key == "audible_breath" || key == "side_cavities") {
            if (value != "true" && value != "false") { reason = "expected true or false for " + key; }
            else if (key == "phonation") { candidate.phonation = value == "true"; }
            else if (key == "audible_breath") { candidate.audible_breath = value == "true"; }
            else { candidate.side_cavities = value == "true"; }
        } else {
            bool found = false;
            for (const Field& field : fields) {
                if (key != field.name) { continue; }
                found = true;
                std::istringstream parser(value);
                double number_value = 0.0;
                std::string extra{};
                if (!(parser >> number_value) || (parser >> extra) || !std::isfinite(number_value)
                    || number_value < field.low || number_value > field.high) {
                    reason = "invalid or out-of-range value for " + key + " [" + std::to_string(field.low)
                        + ", " + std::to_string(field.high) + "]";
                } else { *field.value = number_value; }
                break;
            }
            if (!found) { reason = "unknown field " + key; }
        }
        if (!reason.empty()) {
            error = "body config line " + std::to_string(number) + ": " + reason;
            return false;
        }
    }
    if (input.bad()) { error = "body config read failed"; return false; }
    if (candidate.abducted_area < candidate.adducted_area) {
        error = "abducted_area must be at least adducted_area"; return false;
    }
    output = candidate;
    error.clear();
    return true;
}

bool load_body_config(const std::string& path, VoiceProfile& output, std::string& error) {
    std::ifstream input(path);
    if (!input) { error = "cannot open body config: " + path; return false; }
    return read_body_config(input, output, error);
}
}
