#include "ducklips/vowel_library.hpp"
#include <algorithm>
#include <cmath>

namespace ducklips {
namespace {

constexpr double formant_search_limit = 2.0 * pi * 6000.0; // rad/s
constexpr std::size_t inversion_mode_capacity = 12;
// Per-parameter regularization: light on the tongue, heavier on the
// anatomical residual so it only absorbs what the tongue cannot do.
constexpr std::array<double, shape_parameter_count> regularization_weight{
    0.0004, 0.0004, 0.0004, 0.0004, 0.0004, 0.02, 0.02, 0.02, 0.02};
// Neutral tongue posture the inversion is regularized toward, and the
// anatomical box it may not leave (body centre cm, jaw 0..1).
constexpr std::array<double, shape_parameter_count> neutral_posture{2.8, 6.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
constexpr std::array<double, shape_parameter_count> lowest_posture{1.4, 4.6, 0.0, -1.0, -1.0, -1.2, -1.2, -1.2, -1.2};
constexpr std::array<double, shape_parameter_count> highest_posture{4.6, 7.8, 1.0, 1.0, 1.0, 1.2, 1.2, 1.2, 1.2};
constexpr std::array<double, fitted_formant_count> formant_weight{1.0, 1.0, 0.6};

using SquareMatrix = std::array<std::array<double, shape_parameter_count>, shape_parameter_count>;
using ShapeVector = std::array<double, shape_parameter_count>;

// Solves a small symmetric positive definite system by Gaussian elimination
// with partial pivoting. Returns false for a singular system.
bool solve_dense(SquareMatrix matrix, ShapeVector right_side, ShapeVector& solution) noexcept {
    const std::size_t size = shape_parameter_count;
    for (std::size_t column = 0; column < size; ++column) {
        std::size_t pivot_row = column;
        for (std::size_t row = column + 1; row < size; ++row) {
            if (std::fabs(matrix[row][column]) > std::fabs(matrix[pivot_row][column])) {
                pivot_row = row;
            }
        }
        if (std::fabs(matrix[pivot_row][column]) < 1.0e-300) {
            return false;
        }
        std::swap(matrix[column], matrix[pivot_row]);
        std::swap(right_side[column], right_side[pivot_row]);
        for (std::size_t row = column + 1; row < size; ++row) {
            double factor = matrix[row][column] / matrix[column][column];
            for (std::size_t inner = column; inner < size; ++inner) {
                matrix[row][inner] -= factor * matrix[column][inner];
            }
            right_side[row] -= factor * right_side[column];
        }
    }
    for (std::size_t row = size; row > 0; --row) {
        std::size_t index = row - 1;
        double sum = right_side[index];
        for (std::size_t inner = index + 1; inner < size; ++inner) {
            sum -= matrix[index][inner] * solution[inner];
        }
        solution[index] = sum / matrix[index][index];
    }
    return true;
}

// A_x dg/dA_x for a tube branch through half of a section of area A_x and length l_x.
double half_section_area_derivative(double conductance, double length, double area, double density) noexcept {
    double value = conductance * conductance * density * length / (2.0 * area);
    return value;
}

}

void FormantInverter::prepare(const VoiceProfile& voice, std::size_t oral_sections, std::size_t nasal_sections,
                              const AirProperties& air) {
    air_ = air;
    model_.prepare(voice, oral_sections, nasal_sections);
    model_.allocate(tract_);
    network_.prepare(oral_sections, nasal_sections, model_.port_section(), model_.side_cavities());
    matrix_.prepare(network_.node_count(), network_.tree_edges(), model_.port_section());
    solver_.prepare(network_.node_count());
    basis_.prepare(network_.node_count(), inversion_mode_capacity);
    empty_.prepare(network_.node_count(), 0);
    area_gradient_.assign(oral_sections * shape_parameter_count, 0.0);
}

const TractModel& FormantInverter::model() const noexcept {
    return model_;
}

bool FormantInverter::solve_posture(const VowelShape& shape, std::array<std::size_t, fitted_formant_count>& oral_modes) {
    model_.build_vowel(shape, tract_);
    network_.assemble(tract_, air_, NetworkLosses{});
    network_.write_operator(matrix_);
    solver_.solve(network_, matrix_, formant_search_limit, empty_, basis_);
    // With the port closed the nasal branch is decoupled; its modes carry no
    // oral energy. Keep modes whose energy is mostly oral.
    std::size_t found = 0;
    const std::size_t oral = network_.oral_sections();
    for (std::size_t mode = 0; mode < basis_.mode_count && found < fitted_formant_count; ++mode) {
        double oral_energy = 0.0;
        for (std::size_t node = 0; node < oral; ++node) {
            double pressure = basis_.pressure(mode, node);
            oral_energy += network_.nodes[node].compliance * pressure * pressure;
        }
        if (oral_energy > 0.5) {
            oral_modes[found] = mode;
            found += 1;
        }
    }
    bool complete = found == fitted_formant_count;
    return complete;
}

bool FormantInverter::formants(const VowelShape& shape, FormantArray& output) {
    std::array<std::size_t, fitted_formant_count> modes{};
    if (!solve_posture(shape, modes)) {
        return false;
    }
    for (std::size_t formant = 0; formant < fitted_formant_count; ++formant) {
        output[formant] = basis_.angular_frequency[modes[formant]] / (2.0 * pi);
    }
    return true;
}

bool FormantInverter::formant_jacobian(const VowelShape& shape, FormantArray& frequencies, ShapeJacobian& jacobian) {
    // Chain rule: eigenvalue sensitivity to each section's log area, times
    // the tongue geometry's d ln(area) / d(posture parameter).
    model_.log_area_gradient(shape, std::span<double>(area_gradient_));
    std::array<std::size_t, fitted_formant_count> modes{};
    if (!solve_posture(shape, modes)) {
        return false;
    }
    const std::size_t oral = network_.oral_sections();
    for (std::size_t formant = 0; formant < fitted_formant_count; ++formant) {
        std::size_t mode = modes[formant];
        double eigenvalue = basis_.eigenvalue[mode];
        frequencies[formant] = basis_.angular_frequency[mode] / (2.0 * pi);
        ShapeVector gradient{};
        for (std::size_t section = 0; section < oral; ++section) {
            double area = tract_.oral_area[section];
            double length = tract_.oral_length[section];
            double pressure = basis_.pressure(mode, section);
            double sensitivity = -eigenvalue * network_.nodes[section].compliance * pressure * pressure;
            if (section > 0) {
                const AcousticBranch& below = network_.branches[network_.oral_branch(section - 1)];
                double difference = basis_.pressure(mode, section - 1) - pressure;
                sensitivity += half_section_area_derivative(below.conductance, length, area, air_.density)
                    * difference * difference;
            }
            if (section + 1 < oral) {
                const AcousticBranch& above = network_.branches[network_.oral_branch(section)];
                double difference = pressure - basis_.pressure(mode, section + 1);
                sensitivity += half_section_area_derivative(above.conductance, length, area, air_.density)
                    * difference * difference;
            }
            const AcousticBranch& wall = network_.branches[network_.wall_branch(section)];
            sensitivity += 0.5 * wall.conductance * pressure * pressure;
            for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
                gradient[parameter] += sensitivity * area_gradient_[section * shape_parameter_count + parameter];
            }
        }
        for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
            // d ln F = d lambda / (2 lambda)
            jacobian[formant][parameter] = gradient[parameter] / (2.0 * eigenvalue);
        }
    }
    return true;
}

FormantFit FormantInverter::fit(Phone phone, const FormantArray& target, VowelShape& shape) {
    FormantFit result{};
    result.phone = phone;
    result.target = target;
    FormantArray frequencies{};
    ShapeJacobian jacobian{};
    double damping = 1.0e-2;
    std::size_t iteration = 0;
    for (iteration = 0; iteration < 60; ++iteration) {
        if (!formant_jacobian(shape, frequencies, jacobian)) {
            break;
        }
        FormantArray residual{};
        double cost = 0.0;
        double largest = 0.0;
        for (std::size_t formant = 0; formant < fitted_formant_count; ++formant) {
            residual[formant] = formant_weight[formant] * std::log(frequencies[formant] / target[formant]);
            cost += residual[formant] * residual[formant];
            largest = std::max(largest, std::fabs(frequencies[formant] / target[formant] - 1.0));
        }
        for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
            double offset = shape.coefficients[parameter] - neutral_posture[parameter];
            cost += regularization_weight[parameter] * offset * offset;
        }
        if (largest < 0.003) {
            break;
        }
        bool improved = false;
        for (std::size_t attempt = 0; attempt < 8 && !improved; ++attempt) {
            SquareMatrix normal{};
            ShapeVector gradient{};
            for (std::size_t row = 0; row < shape_parameter_count; ++row) {
                for (std::size_t column = 0; column < shape_parameter_count; ++column) {
                    double sum = 0.0;
                    for (std::size_t formant = 0; formant < fitted_formant_count; ++formant) {
                        double weight = formant_weight[formant] * formant_weight[formant];
                        sum += weight * jacobian[formant][row] * jacobian[formant][column];
                    }
                    normal[row][column] = sum;
                }
                normal[row][row] += regularization_weight[row] + damping;
                double projected = 0.0;
                for (std::size_t formant = 0; formant < fitted_formant_count; ++formant) {
                    projected += formant_weight[formant] * jacobian[formant][row] * residual[formant];
                }
                gradient[row] = -(projected + regularization_weight[row] * (shape.coefficients[row] - neutral_posture[row]));
            }
            ShapeVector step{};
            if (!solve_dense(normal, gradient, step)) {
                damping *= 10.0;
                continue;
            }
            VowelShape trial = shape;
            for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
                double moved = trial.coefficients[parameter] + step[parameter];
                trial.coefficients[parameter] = std::clamp(moved, lowest_posture[parameter], highest_posture[parameter]);
            }
            FormantArray trial_frequencies{};
            if (!formants(trial, trial_frequencies)) {
                damping *= 4.0;
                continue;
            }
            double trial_cost = 0.0;
            for (std::size_t formant = 0; formant < fitted_formant_count; ++formant) {
                double value = formant_weight[formant] * std::log(trial_frequencies[formant] / target[formant]);
                trial_cost += value * value;
            }
            for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
                double offset = trial.coefficients[parameter] - neutral_posture[parameter];
                trial_cost += regularization_weight[parameter] * offset * offset;
            }
            if (trial_cost < cost) {
                shape = trial;
                damping = std::max(damping / 3.0, 1.0e-6);
                improved = true;
            } else {
                damping *= 4.0;
            }
        }
        if (!improved) {
            break;
        }
    }
    FormantArray achieved{};
    if (formants(shape, achieved)) {
        result.achieved = achieved;
    }
    double largest = 0.0;
    for (std::size_t formant = 0; formant < fitted_formant_count; ++formant) {
        largest = std::max(largest, std::fabs(result.achieved[formant] / target[formant] - 1.0));
    }
    result.iterations = iteration;
    result.largest_relative_error = largest;
    return result;
}

void VowelLibrary::build(const VoiceProfile& voice, std::size_t oral_sections, std::size_t nasal_sections,
                         const AirProperties& air) {
    FormantInverter inverter{};
    inverter.prepare(voice, oral_sections, nasal_sections, air);
    double frequency_scale = 1.0 / inverter.model().length_scale();
    fits_.clear();
    for (std::size_t index = 0; index < phone_count; ++index) {
        Phone phone = static_cast<Phone>(index);
        const PhoneSpec& spec = phone_spec(phone);
        if (spec.kind != PhoneClass::vowel) {
            continue;
        }
        VowelShape shape{};
        shape.lip_area = spec.lip_area_cm2 * square_centimetre;
        shape.lip_protrusion = spec.lip_protrusion_cm * centimetre;
        FormantArray target{spec.f1 * frequency_scale, spec.f2 * frequency_scale, spec.f3 * frequency_scale};
        FormantFit fit = inverter.fit(phone, target, shape);
        shapes_[index] = shape;
        fits_.push_back(fit);
    }
}

const VowelShape& VowelLibrary::shape(Phone phone) const noexcept {
    return shapes_[static_cast<std::size_t>(phone)];
}

const std::vector<FormantFit>& VowelLibrary::fits() const noexcept {
    return fits_;
}

}
