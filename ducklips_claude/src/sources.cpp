#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>

namespace ducklips {

void TwoMassGlottis::reset(const FoldParameters& folds) noexcept {
    folds_ = folds;
    strips_ = std::clamp<std::size_t>(folds.strips, 1, maximum_fold_strips);
    shape_.fill(1.0);
    opening_.fill(1.0);
    if (strips_ > 1) {
        // Strip centres from the anterior commissure (0) to the vocal process
        // (1): half-sine vibration shape; rest opening growing toward the back.
        double shape_sum = 0.0;
        double opening_sum = 0.0;
        for (std::size_t strip = 0; strip < strips_; ++strip) {
            double along = (static_cast<double>(strip) + 0.5) / static_cast<double>(strips_);
            shape_[strip] = std::sin(3.14159265358979323846 * along);
            opening_[strip] = along * along;
            shape_sum += shape_[strip];
            opening_sum += opening_[strip];
        }
        for (std::size_t strip = 0; strip < strips_; ++strip) {
            shape_[strip] *= static_cast<double>(strips_) / shape_sum;
            opening_[strip] *= static_cast<double>(strips_) / opening_sum;
        }
    }
    lower_position_.fill(0.0);
    lower_velocity_.fill(0.0);
    upper_position_.fill(0.0);
    upper_velocity_.fill(0.0);
    rest_area_ = 0.0;
}

// The opening of one strip (m^2, negative when the folds press together),
// from the two folds' displacements added.
double TwoMassGlottis::strip_gap(double rest_area, std::size_t strip, double position_sum) const noexcept {
    double share = 1.0 / static_cast<double>(strips_);
    double gap = share * (rest_area * opening_[strip] + folds_.length * shape_[strip] * position_sum);
    return gap;
}

void TwoMassGlottis::advance(double transglottal_pressure, double rest_area, double tension, double sample_period,
                             std::size_t substeps, double surface_turbulence) noexcept {
    rest_area_ = rest_area;
    const double step = sample_period / static_cast<double>(substeps);
    const double length = folds_.length;
    const double share = 1.0 / static_cast<double>(strips_);
    // Each fold's tension: the left at Q, the right at Q times the asymmetry.
    std::array<double, 2> side_tension{tension, tension * folds_.asymmetry};
    std::array<double, 2> lower_mass{};
    std::array<double, 2> upper_mass{};
    std::array<double, 2> lower_stiffness{};
    std::array<double, 2> upper_stiffness{};
    std::array<double, 2> coupling{};
    std::array<double, 2> lower_critical{};
    std::array<double, 2> upper_critical{};
    for (std::size_t side = 0; side < 2; ++side) {
        lower_mass[side] = folds_.lower_mass / side_tension[side];
        upper_mass[side] = folds_.upper_mass / side_tension[side];
        lower_stiffness[side] = folds_.lower_stiffness * side_tension[side];
        upper_stiffness[side] = folds_.upper_stiffness * side_tension[side];
        coupling[side] = folds_.coupling_stiffness * side_tension[side];
        lower_critical[side] = 2.0 * std::sqrt(lower_mass[side] * lower_stiffness[side]);
        upper_critical[side] = 2.0 * std::sqrt(upper_mass[side] * upper_stiffness[side]);
    }
    for (std::size_t substep = 0; substep < substeps; ++substep) {
        // Each strip: Bernoulli pressure on the lower masses where both
        // levels are open, the full drop where only the upper level has
        // closed, and contact where the folds press together (each fold takes
        // half the overlap). Forces reach the masses through the strip's
        // vibration shape.
        double lower_sum = lower_position_[0] + lower_position_[1];
        double upper_sum = upper_position_[0] + upper_position_[1];
        double lower_drive = 0.0;
        double lower_overlap = 0.0; // shape-weighted half-overlap (m), negative when pressed
        double upper_overlap = 0.0;
        double lower_closed = 0.0;
        double upper_closed = 0.0;
        for (std::size_t strip = 0; strip < strips_; ++strip) {
            double lower_gap = strip_gap(rest_area, strip, lower_sum);
            double upper_gap = strip_gap(rest_area, strip, upper_sum);
            double pressure = 0.0;
            if (lower_gap > 0.0 && upper_gap > 0.0) {
                double ratio = std::min(lower_gap, upper_gap) / lower_gap;
                pressure = transglottal_pressure * (1.0 - ratio * ratio) + surface_turbulence;
            } else if (lower_gap > 0.0) {
                pressure = transglottal_pressure;
            }
            lower_drive += share * shape_[strip] * pressure;
            if (lower_gap < 0.0) {
                lower_overlap += share * shape_[strip] * lower_gap / (2.0 * length * share);
                lower_closed += share;
            }
            if (upper_gap < 0.0) {
                upper_overlap += share * shape_[strip] * upper_gap / (2.0 * length * share);
                upper_closed += share;
            }
        }
        for (std::size_t side = 0; side < 2; ++side) {
            double lower_contact = folds_.collision_ratio * lower_stiffness[side];
            double upper_contact = folds_.collision_ratio * upper_stiffness[side];
            double lower_damping = (folds_.lower_damping_ratio + folds_.closed_damping_increase * lower_closed)
                * lower_critical[side];
            double upper_damping = (folds_.upper_damping_ratio + folds_.closed_damping_increase * upper_closed)
                * upper_critical[side];
            double lower_force = length * folds_.lower_thickness * lower_drive
                - lower_stiffness[side] * lower_position_[side]
                - coupling[side] * (lower_position_[side] - upper_position_[side])
                - lower_contact * lower_overlap
                - lower_damping * lower_velocity_[side];
            double upper_force = -upper_stiffness[side] * upper_position_[side]
                - coupling[side] * (upper_position_[side] - lower_position_[side])
                - upper_contact * upper_overlap
                - upper_damping * upper_velocity_[side];
            lower_velocity_[side] += step * lower_force / lower_mass[side];
            upper_velocity_[side] += step * upper_force / upper_mass[side];
            lower_position_[side] += step * lower_velocity_[side];
            upper_position_[side] += step * upper_velocity_[side];
        }
    }
}

double TwoMassGlottis::lower_area() const noexcept {
    double area = rest_area_ + folds_.length * (lower_position_[0] + lower_position_[1]);
    return area;
}

double TwoMassGlottis::upper_area() const noexcept {
    double area = rest_area_ + folds_.length * (upper_position_[0] + upper_position_[1]);
    return area;
}

// The flow area: per strip, the narrower of the two openings. One strip
// keeps the signed value (negative while pressed shut).
double TwoMassGlottis::minimum_area() const noexcept {
    if (folds_.contact_softness > 0.0) {
        // Each strip seals smoothly: softness * log(1 + exp(gap / softness)).
        double softness = folds_.contact_softness / static_cast<double>(strips_);
        double roundness = std::max(folds_.medial_roundness, 1.0e-12) / static_cast<double>(strips_);
        double area = 0.0;
        for (std::size_t strip = 0; strip < strips_; ++strip) {
            // The narrowest point slides smoothly along the rounded medial
            // surface between the two masses: a smooth minimum, so the flow
            // area has no corner when the masses exchange the narrowest point.
            double lower_gap = strip_gap(rest_area_, strip, lower_position_[0] + lower_position_[1]);
            double upper_gap = strip_gap(rest_area_, strip, upper_position_[0] + upper_position_[1]);
            double low = std::min(lower_gap, upper_gap);
            double narrowest = low - roundness * std::log1p(std::exp(-std::fabs(lower_gap - upper_gap) / roundness));
            double ratio = narrowest / softness;
            double smooth = ratio > 30.0 ? ratio : std::log1p(std::exp(ratio));
            area += softness * smooth;
        }
        return area;
    }
    if (strips_ == 1) {
        double area = std::min(lower_area(), upper_area());
        return area;
    }
    double area = 0.0;
    for (std::size_t strip = 0; strip < strips_; ++strip) {
        double narrowest = std::min(strip_gap(rest_area_, strip, lower_position_[0] + lower_position_[1]),
                                    strip_gap(rest_area_, strip, upper_position_[0] + upper_position_[1]));
        area += std::max(narrowest, 0.0);
    }
    return area;
}

double orifice_flow(double area, double pressure_drop, double series_impedance, const AirProperties& air) noexcept {
    if (area <= 0.0 || pressure_drop == 0.0) {
        return 0.0;
    }
    double beta = 2.0 * area * area / air.density;
    double magnitude = std::fabs(pressure_drop);
    double impedance_term = beta * series_impedance;
    double root = std::sqrt(impedance_term * impedance_term + 4.0 * beta * magnitude);
    double flow = 2.0 * beta * magnitude / (impedance_term + root);
    if (pressure_drop < 0.0) {
        flow = -flow;
    }
    return flow;
}

void IntraoralPressure::reset() noexcept {
    pressure_ = 0.0;
}

double IntraoralPressure::advance(double inflow, double outlet_area, double compliance, double sample_period,
                                  const AirProperties& air) noexcept {
    // alpha s^2 + beta s - |gamma| = 0 with P = sign(gamma) s^2.
    double alpha = compliance / sample_period;
    double beta = std::max(outlet_area, 0.0) * std::sqrt(2.0 / air.density);
    double gamma = alpha * pressure_ + inflow;
    double magnitude = std::fabs(gamma);
    double root = 0.0;
    if (magnitude > 0.0) {
        root = 2.0 * magnitude / (beta + std::sqrt(beta * beta + 4.0 * alpha * magnitude));
    }
    double pressure = root * root;
    if (gamma < 0.0) {
        pressure = -pressure;
    }
    pressure_ = pressure;
    return pressure_;
}

double IntraoralPressure::pressure() const noexcept {
    return pressure_;
}

void NoiseGenerator::seed(std::uint64_t value) noexcept {
    // splitmix64 so that nearby seeds give unrelated streams; never zero.
    std::uint64_t mixed = value + 0x9E3779B97F4A7C15ULL;
    mixed = (mixed ^ (mixed >> 30)) * 0xBF58476D1CE4E5B9ULL;
    mixed = (mixed ^ (mixed >> 27)) * 0x94D049BB133111EBULL;
    mixed = mixed ^ (mixed >> 31);
    if (mixed == 0) {
        mixed = 0x2545F4914F6CDD1DULL;
    }
    state_ = mixed;
    has_spare_ = false;
    spare_ = 0.0;
}

double NoiseGenerator::uniform() noexcept {
    std::uint64_t x = state_;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    state_ = x;
    std::uint64_t scrambled = x * 0x2545F4914F6CDD1DULL;
    double value = static_cast<double>(scrambled >> 11) * (1.0 / 9007199254740992.0);
    return value;
}

double NoiseGenerator::gaussian() noexcept {
    if (has_spare_) {
        has_spare_ = false;
        return spare_;
    }
    double first = uniform();
    double second = uniform();
    double radius = std::sqrt(-2.0 * std::log(1.0 - first));
    double angle = 2.0 * pi * second;
    spare_ = radius * std::sin(angle);
    has_spare_ = true;
    double value = radius * std::cos(angle);
    return value;
}

void SlowNoise::prepare(std::uint64_t seed, double cutoff_hz, double update_rate_hz) noexcept {
    source_.seed(seed);
    // Three equal poles, each placed so the cascade is 3 dB down at the
    // cutoff (1 / sqrt(2^(1/3) - 1) = 1.961).
    coefficient_ = 1.0 - std::exp(-2.0 * pi * 1.961 * cutoff_hz / update_rate_hz);
    // Variance of the cascade driven by unit white noise, with x = (1 - c)^2:
    // c^6 sum C(n+2,2)^2 x^n = c^6 (1 + 4x + x^2) / (1 - x)^5.
    double pole = 1.0 - coefficient_;
    double x = pole * pole;
    double c3 = coefficient_ * coefficient_ * coefficient_;
    double variance = c3 * c3 * (1.0 + 4.0 * x + x * x) / std::pow(1.0 - x, 5.0);
    normalization_ = 1.0 / std::sqrt(variance);
    stages_.fill(0.0);
}

double SlowNoise::next() noexcept {
    double input = source_.gaussian();
    for (std::size_t stage = 0; stage < stages_.size(); ++stage) {
        stages_[stage] += coefficient_ * (input - stages_[stage]);
        input = stages_[stage];
    }
    double output = input * normalization_;
    return output;
}

void TurbulenceNoise::prepare(std::uint64_t seed, double highpass_hz, double sample_rate_hz) noexcept {
    source_.seed(seed);
    coefficient_ = std::exp(-2.0 * pi * highpass_hz / sample_rate_hz);
    previous_input_ = 0.0;
    previous_output_ = 0.0;
}

double TurbulenceNoise::next() noexcept {
    double input = source_.gaussian();
    double output = coefficient_ * (previous_output_ + input - previous_input_);
    previous_input_ = input;
    previous_output_ = output;
    return output;
}

}
