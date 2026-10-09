#include "ducklips/modal.hpp"
#include <algorithm>
#include <cmath>

namespace ducklips {
namespace {

double dot(std::span<const double> left, std::span<const double> right) noexcept {
    double sum = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        sum += left[index] * right[index];
    }
    return sum;
}

double normalize(std::span<double> vector) noexcept {
    double norm = std::sqrt(dot(vector, vector));
    if (norm > 0.0 && std::isfinite(norm)) {
        double inverse = 1.0 / norm;
        for (std::size_t index = 0; index < vector.size(); ++index) {
            vector[index] *= inverse;
        }
    }
    return norm;
}

double pressure_difference(const AcousticNetwork& network, const ModalBasis& basis, std::size_t mode,
                           std::size_t branch) noexcept {
    const AcousticBranch& link = network.branches[branch];
    double difference = basis.pressure(mode, link.from);
    if (link.to != ground_node) {
        difference -= basis.pressure(mode, link.to);
    }
    return difference;
}

}

void ModalBasis::prepare(std::size_t nodes, std::size_t modes) {
    node_count = nodes;
    capacity = modes;
    mode_count = 0;
    eigenvalue.assign(modes, 0.0);
    angular_frequency.assign(modes, 0.0);
    damping.assign(modes, 0.0);
    pressure_shape.assign(modes * nodes, 0.0);
    unit_shape.assign(modes * nodes, 0.0);
}

double ModalBasis::pressure(std::size_t mode, std::size_t node) const noexcept {
    return pressure_shape[mode * node_count + node];
}

double branch_flow_shape(const AcousticNetwork& network, const ModalBasis& basis, std::size_t mode,
                         std::size_t branch) noexcept {
    double difference = pressure_difference(network, basis, mode, branch);
    double flow = network.branches[branch].conductance * difference / basis.angular_frequency[mode];
    return flow;
}

void ModalSolver::prepare(std::size_t node_count) {
    work_.assign(node_count, 0.0);
    product_.assign(node_count, 0.0);
}

double ModalSolver::rayleigh_quotient(const TreeOperator& matrix, std::span<const double> vector) noexcept {
    matrix.multiply(vector, std::span<double>(product_));
    double quotient = dot(vector, std::span<const double>(product_));
    return quotient;
}

bool ModalSolver::refine_rayleigh(TreeOperator& matrix, std::span<double> vector, double& eigenvalue,
                                  ModalSolveReport& report) noexcept {
    const std::size_t count = vector.size();
    double shift = rayleigh_quotient(matrix, vector);
    double residual = 0.0;
    for (std::size_t iteration = 0; iteration < 4; ++iteration) {
        matrix.factor(shift);
        report.factorizations += 1;
        for (std::size_t index = 0; index < count; ++index) {
            work_[index] = vector[index];
        }
        matrix.solve(std::span<double>(work_));
        double norm = normalize(std::span<double>(work_));
        if (!(norm > 0.0) || !std::isfinite(norm)) {
            return false;
        }
        for (std::size_t index = 0; index < count; ++index) {
            vector[index] = work_[index];
        }
        shift = rayleigh_quotient(matrix, vector);
        double residual_square = 0.0;
        for (std::size_t index = 0; index < count; ++index) {
            double component = product_[index] - shift * vector[index];
            residual_square += component * component;
        }
        residual = std::sqrt(residual_square);
        if (residual <= 1.0e-10 * std::max(std::fabs(shift), 1.0)) {
            break;
        }
    }
    eigenvalue = shift;
    bool converged = residual <= 1.0e-7 * std::max(std::fabs(shift), 1.0);
    return converged;
}

bool ModalSolver::index_verified(TreeOperator& matrix, std::size_t index, double eigenvalue,
                                 ModalSolveReport& report) noexcept {
    std::size_t below = matrix.count_below(eigenvalue * (1.0 - 1.0e-8));
    std::size_t at_or_below = matrix.count_below(eigenvalue * (1.0 + 1.0e-8));
    report.factorizations += 2;
    bool verified = below <= index && at_or_below >= index + 1;
    return verified;
}

double ModalSolver::bisect_eigenvalue(TreeOperator& matrix, std::size_t index, double upper,
                                      ModalSolveReport& report) noexcept {
    double low = 0.0;
    double high = upper * (1.0 + 1.0e-9);
    for (std::size_t iteration = 0; iteration < 90; ++iteration) {
        if (high - low <= 1.0e-13 * high) {
            break;
        }
        double middle = 0.5 * (low + high);
        std::size_t below = matrix.count_below(middle);
        report.factorizations += 1;
        if (below >= index + 1) {
            high = middle;
        } else {
            low = middle;
        }
    }
    double eigenvalue = 0.5 * (low + high);
    return eigenvalue;
}

void ModalSolver::deflate(std::span<double> vector, const ModalBasis& output, std::size_t index) noexcept {
    const std::size_t count = vector.size();
    for (std::size_t earlier = 0; earlier < index; ++earlier) {
        std::span<const double> accepted(output.unit_shape.data() + earlier * count, count);
        double overlap = dot(accepted, vector);
        for (std::size_t node = 0; node < count; ++node) {
            vector[node] -= overlap * accepted[node];
        }
    }
}

void ModalSolver::inverse_iterate(TreeOperator& matrix, double eigenvalue, std::span<double> vector,
                                  const ModalBasis& output, std::size_t index, ModalSolveReport& report) noexcept {
    // The start vector may be an exact eigenvector of a neighbouring mode (a
    // failed warm start converges onto one). Inverse iteration cannot leave
    // such a vector quickly when the neighbour is close, so a fixed pattern is
    // mixed in, every accepted lower mode is deflated on every pass, and the
    // iteration runs until the residual has converged.
    const std::size_t count = vector.size();
    for (std::size_t node = 0; node < count; ++node) {
        double position = static_cast<double>(node) + 1.0;
        double pattern = std::sin(0.37 * static_cast<double>(index + 1) * position) + 0.5 * std::cos(1.3 * position);
        vector[node] += 0.05 * pattern;
    }
    deflate(vector, output, index);
    normalize(vector);
    matrix.factor(eigenvalue);
    report.factorizations += 1;
    for (std::size_t iteration = 0; iteration < 12; ++iteration) {
        matrix.solve(vector);
        deflate(vector, output, index);
        normalize(vector);
        double quotient = rayleigh_quotient(matrix, vector);
        double residual_square = 0.0;
        for (std::size_t node = 0; node < count; ++node) {
            double component = product_[node] - quotient * vector[node];
            residual_square += component * component;
        }
        if (std::sqrt(residual_square) <= 1.0e-9 * std::max(std::fabs(quotient), 1.0)) {
            break;
        }
    }
}

ModalSolveReport ModalSolver::solve(const AcousticNetwork& network, TreeOperator& matrix,
                                    double maximum_angular_frequency, const ModalBasis& previous,
                                    ModalBasis& output) noexcept {
    ModalSolveReport report{};
    const std::size_t count = network.node_count();
    const double upper = maximum_angular_frequency * maximum_angular_frequency;
    std::size_t available = matrix.count_below(upper);
    report.factorizations += 1;
    std::size_t modes = std::min(available, output.capacity);
    output.mode_count = modes;
    for (std::size_t mode = 0; mode < modes; ++mode) {
        std::span<double> vector(output.unit_shape.data() + mode * count, count);
        bool accepted = false;
        double eigenvalue = 0.0;
        bool warm = mode < previous.mode_count && previous.node_count == count;
        if (warm) {
            for (std::size_t node = 0; node < count; ++node) {
                double root_compliance = std::sqrt(network.nodes[node].compliance);
                vector[node] = root_compliance * previous.pressure_shape[mode * count + node];
            }
            double norm = normalize(vector);
            if (norm > 0.0 && std::isfinite(norm)) {
                accepted = refine_rayleigh(matrix, vector, eigenvalue, report);
            }
            if (accepted) {
                accepted = index_verified(matrix, mode, eigenvalue, report);
            }
            // A converged eigenvector is orthogonal to every other mode to
            // rounding; any visible overlap means the refinement landed on a
            // mode that is already in the basis.
            for (std::size_t earlier = 0; earlier < mode && accepted; ++earlier) {
                std::span<const double> lower(output.unit_shape.data() + earlier * count, count);
                if (std::fabs(dot(lower, vector)) > 1.0e-4) {
                    accepted = false;
                }
            }
            if (!accepted) {
                for (std::size_t node = 0; node < count; ++node) {
                    double root_compliance = std::sqrt(network.nodes[node].compliance);
                    vector[node] = root_compliance * previous.pressure_shape[mode * count + node];
                }
            }
        } else {
            for (std::size_t node = 0; node < count; ++node) {
                vector[node] = 0.0;
            }
        }
        if (accepted) {
            report.warm_modes += 1;
        } else {
            eigenvalue = bisect_eigenvalue(matrix, mode, upper, report);
            inverse_iterate(matrix, eigenvalue, vector, output, mode, report);
            report.bisected_modes += 1;
        }
        output.eigenvalue[mode] = eigenvalue;
        output.angular_frequency[mode] = std::sqrt(std::max(eigenvalue, 0.0));
        for (std::size_t node = 0; node < count; ++node) {
            double root_compliance = std::sqrt(network.nodes[node].compliance);
            output.pressure_shape[mode * count + node] = vector[node] / root_compliance;
        }
    }
    compute_damping(network, output);
    return report;
}

void ModalSolver::compute_damping(const AcousticNetwork& network, ModalBasis& basis) noexcept {
    const std::size_t count = network.node_count();
    for (std::size_t mode = 0; mode < basis.mode_count; ++mode) {
        double frequency = basis.angular_frequency[mode];
        double root_frequency = std::sqrt(frequency);
        double decay = 0.0;
        for (std::size_t branch = 0; branch < network.branches.size(); ++branch) {
            const AcousticBranch& link = network.branches[branch];
            if (link.conductance <= 0.0) {
                continue;
            }
            double flow = branch_flow_shape(network, basis, mode, branch);
            double resistance = link.viscous_coefficient * root_frequency + link.resistance
                + link.radiation_coefficient * frequency * frequency;
            decay += 0.5 * resistance * flow * flow;
        }
        for (std::size_t node = 0; node < count; ++node) {
            double pressure = basis.pressure(mode, node);
            double conductance = network.nodes[node].thermal_coefficient * root_frequency;
            decay += 0.5 * conductance * pressure * pressure;
        }
        double hertz = frequency / (2.0 * pi);
        double onset = network.losses.cross_mode_onset;
        if (hertz > onset) {
            double excess = (hertz - onset) / onset;
            // An amplitude decay alpha is a bandwidth of alpha / pi Hz.
            decay += pi * network.losses.cross_mode_bandwidth * excess * excess;
        }
        basis.damping[mode] = decay;
    }
}

void ModalResonator::prepare(std::size_t capacity, std::size_t node_count, std::size_t branch_count) {
    mode_count = 0;
    state.assign(capacity, std::complex<double>(0.0, 0.0));
    propagator.assign(capacity, std::complex<double>(0.0, 0.0));
    input_gain.assign(capacity, std::complex<double>(0.0, 0.0));
    glottal_weight.assign(capacity, 0.0);
    aspiration_weight.assign(capacity, 0.0);
    frication_weight.assign(capacity, 0.0);
    lip_output.assign(capacity, 0.0);
    nostril_output.assign(capacity, 0.0);
    wall_output.assign(capacity, 0.0);
    charge_.assign(node_count, 0.0);
    momentum_.assign(branch_count, 0.0);
    glottal_impedance = 0.0;
}

void ModalResonator::clear() noexcept {
    for (std::size_t mode = 0; mode < state.size(); ++mode) {
        state[mode] = std::complex<double>(0.0, 0.0);
    }
}

void ModalResonator::configure(const AcousticNetwork& network, const ModalBasis& basis, double sample_period,
                               std::size_t aspiration_branch, std::size_t frication_branch) noexcept {
    mode_count = basis.mode_count;
    const std::size_t lip_node = network.oral_sections() - 1;
    const std::size_t nostril_node = network.nostril_node();
    const double lips = network.branches[network.lip_branch()].conductance;
    const double nostril_conductance = network.branches[network.nostril_branch()].conductance;
    double impedance = 0.0;
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
        std::complex<double> pole(-basis.damping[mode], basis.angular_frequency[mode]);
        std::complex<double> step = std::exp(pole * sample_period);
        propagator[mode] = step;
        input_gain[mode] = (step - std::complex<double>(1.0, 0.0)) / pole;
        glottal_weight[mode] = basis.pressure(mode, 0);
        aspiration_weight[mode] = branch_flow_shape(network, basis, mode, aspiration_branch);
        frication_weight[mode] = branch_flow_shape(network, basis, mode, frication_branch);
        lip_output[mode] = basis.pressure(mode, lip_node) * lips;
        nostril_output[mode] = basis.pressure(mode, nostril_node) * nostril_conductance;
        double wall_sum = 0.0;
        for (std::size_t node = 0; node < network.oral_sections(); ++node) {
            wall_sum += basis.pressure(mode, node) * network.branches[network.wall_branch(node)].conductance;
        }
        wall_output[mode] = wall_sum;
        impedance += glottal_weight[mode] * glottal_weight[mode] * input_gain[mode].real();
    }
    glottal_impedance = impedance;
    lip_conductance = lips;
    frication_at_lips = frication_branch == network.lip_branch();
}

void ModalResonator::reproject(const AcousticNetwork& old_network, const ModalBasis& old_basis,
                               const AcousticNetwork& new_network, const ModalBasis& new_basis,
                               ReprojectionInvariant invariant) noexcept {
    const std::size_t node_total = old_network.node_count();
    const std::size_t branch_total = old_network.branch_count();
    // Node value v_i, so that a'_k = sum_i phi'_k(i) v_i.
    for (std::size_t node = 0; node < node_total; ++node) {
        double pressure = 0.0;
        for (std::size_t mode = 0; mode < old_basis.mode_count; ++mode) {
            pressure += old_basis.pressure(mode, node) * state[mode].real();
        }
        double old_compliance = old_network.nodes[node].compliance;
        double new_compliance = new_network.nodes[node].compliance;
        double weight = std::min(old_compliance, new_compliance);
        if (invariant == ReprojectionInvariant::charge_momentum) {
            weight = old_compliance;
        } else if (invariant == ReprojectionInvariant::pressure_flow) {
            weight = new_compliance;
        }
        charge_[node] = weight * pressure;
    }
    // Branch value w_j, so that b'_k = sum_j (dphi'_k(j) / w'_k) w_j. With the
    // old momentum lambda_j = M u = sum_k (dphi_k / w_k) b_k, the flow is
    // g_old lambda_j and the new mode's momentum weight is g_new lambda_j.
    for (std::size_t branch = 0; branch < branch_total; ++branch) {
        double momentum = 0.0;
        for (std::size_t mode = 0; mode < old_basis.mode_count; ++mode) {
            double difference = pressure_difference(old_network, old_basis, mode, branch);
            momentum += difference / old_basis.angular_frequency[mode] * state[mode].imag();
        }
        double old_conductance = old_network.branches[branch].conductance;
        double new_conductance = new_network.branches[branch].conductance;
        double weight = std::min(old_conductance, new_conductance);
        if (invariant == ReprojectionInvariant::charge_momentum) {
            weight = new_conductance;
        } else if (invariant == ReprojectionInvariant::pressure_flow) {
            weight = old_conductance;
        }
        momentum_[branch] = weight * momentum;
    }
    for (std::size_t mode = 0; mode < new_basis.mode_count; ++mode) {
        double real_part = 0.0;
        for (std::size_t node = 0; node < node_total; ++node) {
            real_part += new_basis.pressure(mode, node) * charge_[node];
        }
        double imaginary_part = 0.0;
        double inverse_frequency = 1.0 / new_basis.angular_frequency[mode];
        for (std::size_t branch = 0; branch < branch_total; ++branch) {
            double difference = pressure_difference(new_network, new_basis, mode, branch);
            imaginary_part += difference * inverse_frequency * momentum_[branch];
        }
        state[mode] = std::complex<double>(real_part, imaginary_part);
    }
    for (std::size_t mode = new_basis.mode_count; mode < state.size(); ++mode) {
        state[mode] = std::complex<double>(0.0, 0.0);
    }
}

double ModalResonator::energy() const noexcept {
    double total = 0.0;
    for (std::size_t mode = 0; mode < mode_count; ++mode) {
        total += 0.5 * std::norm(state[mode]);
    }
    return total;
}

}
