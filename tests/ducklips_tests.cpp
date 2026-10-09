// Numerical checks of the DuckLips building blocks against dense references
// and closed forms. Each check prints its measured error.
#include "ducklips/acoustic_network.hpp"
#include "ducklips/articulation.hpp"
#include "ducklips/cues.hpp"
#include "ducklips/modal.hpp"
#include "ducklips/sources.hpp"
#include "ducklips/tree_operator.hpp"
#include "ducklips/vowel_library.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::size_t failures = 0;

void check(bool condition, const char* name, double measured, double limit) {
    std::printf("%-58s %s  measured %.3e  limit %.1e\n", name, condition ? "pass" : "FAIL", measured, limit);
    if (!condition) {
        failures += 1;
    }
}

// Cyclic Jacobi eigenvalues of a dense symmetric matrix (row-major, n x n).
std::vector<double> dense_eigenvalues(std::vector<double> matrix, std::size_t n) {
    for (std::size_t sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t column = row + 1; column < n; ++column) {
                off += matrix[row * n + column] * matrix[row * n + column];
            }
        }
        if (off < 1.0e-26) {
            break;
        }
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                double apq = matrix[p * n + q];
                if (std::fabs(apq) < 1.0e-300) {
                    continue;
                }
                double app = matrix[p * n + p];
                double aqq = matrix[q * n + q];
                double theta = (aqq - app) / (2.0 * apq);
                double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0);
                double s = t * c;
                for (std::size_t k = 0; k < n; ++k) {
                    double akp = matrix[k * n + p];
                    double akq = matrix[k * n + q];
                    matrix[k * n + p] = c * akp - s * akq;
                    matrix[k * n + q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    double apk = matrix[p * n + k];
                    double aqk = matrix[q * n + k];
                    matrix[p * n + k] = c * apk - s * aqk;
                    matrix[q * n + k] = s * apk + c * aqk;
                }
            }
        }
    }
    std::vector<double> values(n, 0.0);
    for (std::size_t index = 0; index < n; ++index) {
        values[index] = matrix[index * n + index];
    }
    std::sort(values.begin(), values.end());
    return values;
}

std::vector<double> dense_from_tree(const ducklips::TreeOperator& matrix) {
    const std::size_t n = matrix.node_count();
    std::vector<double> dense(n * n, 0.0);
    for (std::size_t node = 0; node < n; ++node) {
        dense[node * n + node] = matrix.diagonal[node];
        std::size_t up = matrix.parent(node);
        if (up != ducklips::TreeOperator::no_parent) {
            dense[node * n + up] = matrix.parent_coupling[node];
            dense[up * n + node] = matrix.parent_coupling[node];
        }
    }
    return dense;
}

void test_tree_operator() {
    // A branched tree with deterministic values.
    const std::size_t n = 31;
    std::vector<ducklips::TreeEdge> edges{};
    for (std::size_t node = 1; node < n; ++node) {
        std::size_t parent = (node * 7 + 3) % node;
        edges.push_back(ducklips::TreeEdge{.first_node = parent, .second_node = node});
    }
    ducklips::TreeOperator matrix{};
    matrix.prepare(n, edges, 0);
    for (std::size_t node = 0; node < n; ++node) {
        matrix.diagonal[node] = 2.0 + std::sin(1.7 * static_cast<double>(node));
        matrix.parent_coupling[node] = 0.6 * std::cos(0.9 * static_cast<double>(node));
    }
    std::vector<double> eigenvalues = dense_eigenvalues(dense_from_tree(matrix), n);
    std::size_t worst = 0;
    for (std::size_t index = 0; index + 1 < n; ++index) {
        double shift = 0.5 * (eigenvalues[index] + eigenvalues[index + 1]);
        std::size_t counted = matrix.count_below(shift);
        std::size_t error = counted > index + 1 ? counted - index - 1 : index + 1 - counted;
        worst = std::max(worst, error);
    }
    check(worst == 0, "tree inertia counts match dense eigenvalues", static_cast<double>(worst), 0.0);

    std::vector<double> x(n, 0.0);
    std::vector<double> b(n, 0.0);
    for (std::size_t node = 0; node < n; ++node) {
        x[node] = std::cos(0.3 * static_cast<double>(node));
    }
    double shift = 0.37;
    matrix.multiply(x, b);
    for (std::size_t node = 0; node < n; ++node) {
        b[node] -= shift * x[node];
    }
    matrix.factor(shift);
    matrix.solve(b);
    double error = 0.0;
    for (std::size_t node = 0; node < n; ++node) {
        error = std::max(error, std::fabs(b[node] - x[node]));
    }
    check(error < 1.0e-10, "tree shifted solve reproduces x", error, 1.0e-10);
}

struct Fixture final {
    ducklips::VoiceProfile voice{ducklips::baritone_voice()};
    ducklips::TractModel model{};
    ducklips::TractShape shape{};
    ducklips::AcousticNetwork network{};
    ducklips::TreeOperator matrix{};
    ducklips::ModalSolver solver{};
    ducklips::ModalBasis basis{};
    ducklips::ModalBasis empty{};

    Fixture() : Fixture(ducklips::baritone_voice()) {}
    explicit Fixture(const ducklips::VoiceProfile& body) : voice(body) {
        model.prepare(voice, 64, 24);
        model.allocate(shape);
        network.prepare(64, 24, model.port_section(), model.side_cavities());
        matrix.prepare(network.node_count(), network.tree_edges(), model.port_section());
        solver.prepare(network.node_count());
        basis.prepare(network.node_count(), 40);
        empty.prepare(network.node_count(), 0);
    }
};

void test_modal_solver_against_dense() {
    Fixture fixture{};
    ducklips::VowelShape vowel{};
    vowel.coefficients = {2.6, 6.3, 0.4, 0.2, -0.1, 0.1, -0.1, 0.05, 0.0};
    fixture.model.build_vowel(vowel, fixture.shape);
    fixture.shape.port_area = 0.6 * ducklips::square_centimetre;
    fixture.network.assemble(fixture.shape, ducklips::AirProperties{}, ducklips::NetworkLosses{});
    fixture.network.write_operator(fixture.matrix);
    const double limit = 2.0 * ducklips::pi * 12000.0;
    fixture.solver.solve(fixture.network, fixture.matrix, limit, fixture.empty, fixture.basis);
    std::vector<double> dense = dense_eigenvalues(dense_from_tree(fixture.matrix), fixture.matrix.node_count());
    std::size_t expected = 0;
    while (expected < dense.size() && dense[expected] < limit * limit) {
        expected += 1;
    }
    check(fixture.basis.mode_count == expected, "mode count below 12 kHz equals dense count",
          static_cast<double>(fixture.basis.mode_count), static_cast<double>(expected));
    double worst = 0.0;
    for (std::size_t mode = 0; mode < std::min(expected, fixture.basis.mode_count); ++mode) {
        worst = std::max(worst, std::fabs(fixture.basis.eigenvalue[mode] / dense[mode] - 1.0));
    }
    check(worst < 1.0e-9, "cold eigenvalues match dense Jacobi (relative)", worst, 1.0e-9);

    // C-orthonormality of the computed pressure shapes.
    double orthogonality = 0.0;
    for (std::size_t first = 0; first < fixture.basis.mode_count; ++first) {
        for (std::size_t second = first; second < fixture.basis.mode_count; ++second) {
            double sum = 0.0;
            for (std::size_t node = 0; node < fixture.network.node_count(); ++node) {
                sum += fixture.network.nodes[node].compliance * fixture.basis.pressure(first, node)
                    * fixture.basis.pressure(second, node);
            }
            double expected_value = first == second ? 1.0 : 0.0;
            orthogonality = std::max(orthogonality, std::fabs(sum - expected_value));
        }
    }
    check(orthogonality < 1.0e-8, "pressure shapes are C-orthonormal", orthogonality, 1.0e-8);

    // Warm start after a small deformation equals a cold solve.
    ducklips::ModalBasis previous = fixture.basis;
    vowel.coefficients[1] += 0.02;
    fixture.model.build_vowel(vowel, fixture.shape);
    fixture.shape.port_area = 0.62 * ducklips::square_centimetre;
    fixture.network.assemble(fixture.shape, ducklips::AirProperties{}, ducklips::NetworkLosses{});
    fixture.network.write_operator(fixture.matrix);
    ducklips::ModalSolveReport warm = fixture.solver.solve(fixture.network, fixture.matrix, limit, previous, fixture.basis);
    ducklips::ModalBasis cold{};
    cold.prepare(fixture.network.node_count(), 40);
    fixture.solver.solve(fixture.network, fixture.matrix, limit, fixture.empty, cold);
    double difference = 0.0;
    for (std::size_t mode = 0; mode < cold.mode_count; ++mode) {
        difference = std::max(difference, std::fabs(fixture.basis.eigenvalue[mode] / cold.eigenvalue[mode] - 1.0));
    }
    check(difference < 1.0e-9 && fixture.basis.mode_count == cold.mode_count,
          "warm-started solve equals cold solve", difference, 1.0e-9);
    std::printf("    warm start accepted %zu of %zu modes with %zu factorizations\n", warm.warm_modes,
                fixture.basis.mode_count, warm.factorizations);
}

void test_uniform_tube() {
    // Closed at the glottis, open at the lips without end correction, walls
    // or side cavities: f_k = (2k - 1) c / (4 L).
    ducklips::VoiceProfile plain = ducklips::baritone_voice();
    plain.side_cavities = false;
    Fixture fixture(plain);
    ducklips::VowelShape vowel{};
    fixture.model.build_vowel(vowel, fixture.shape);
    double length = 0.0;
    for (std::size_t section = 0; section < fixture.shape.oral_area.size(); ++section) {
        fixture.shape.oral_area[section] = 3.0 * ducklips::square_centimetre;
        length += fixture.shape.oral_length[section];
    }
    ducklips::NetworkLosses losses{};
    losses.end_correction = 0.0;
    losses.wall_mass_per_area = 1.0e15;
    ducklips::AirProperties air{};
    fixture.network.assemble(fixture.shape, air, losses);
    fixture.network.write_operator(fixture.matrix);
    fixture.solver.solve(fixture.network, fixture.matrix, 2.0 * ducklips::pi * 6000.0, fixture.empty, fixture.basis);
    // Remove the decoupled nasal modes: keep the modes with oral energy.
    std::vector<double> oral_frequencies{};
    for (std::size_t mode = 0; mode < fixture.basis.mode_count; ++mode) {
        double oral_energy = 0.0;
        for (std::size_t node = 0; node < fixture.network.oral_sections(); ++node) {
            double pressure = fixture.basis.pressure(mode, node);
            oral_energy += fixture.network.nodes[node].compliance * pressure * pressure;
        }
        if (oral_energy > 0.5) {
            oral_frequencies.push_back(fixture.basis.angular_frequency[mode] / (2.0 * ducklips::pi));
        }
    }
    double worst = 0.0;
    for (std::size_t mode = 0; mode < 4 && mode < oral_frequencies.size(); ++mode) {
        double exact = static_cast<double>(2 * mode + 1) * air.sound_speed / (4.0 * length);
        worst = std::max(worst, std::fabs(oral_frequencies[mode] / exact - 1.0));
        std::printf("    tube mode %zu: %.1f Hz (quarter-wave %.1f Hz)\n", mode + 1, oral_frequencies[mode], exact);
    }
    check(worst < 5.0e-3 && oral_frequencies.size() >= 4, "uniform tube matches quarter-wave modes", worst, 5.0e-3);
}

void test_reprojection_identity() {
    Fixture fixture{};
    ducklips::VowelShape vowel{};
    fixture.model.build_vowel(vowel, fixture.shape);
    fixture.shape.port_area = 0.4 * ducklips::square_centimetre;
    fixture.network.assemble(fixture.shape, ducklips::AirProperties{}, ducklips::NetworkLosses{});
    fixture.network.write_operator(fixture.matrix);
    fixture.solver.solve(fixture.network, fixture.matrix, 2.0 * ducklips::pi * 12000.0, fixture.empty, fixture.basis);
    ducklips::ModalResonator resonator{};
    resonator.prepare(40, fixture.network.node_count(), fixture.network.branch_count());
    resonator.configure(fixture.network, fixture.basis, 1.0 / 48000.0, 0, 1);
    for (std::size_t mode = 0; mode < fixture.basis.mode_count; ++mode) {
        resonator.state[mode] = std::complex<double>(std::sin(1.0 + mode), std::cos(2.0 * mode));
    }
    std::vector<std::complex<double>> before(resonator.state.begin(), resonator.state.end());
    resonator.reproject(fixture.network, fixture.basis, fixture.network, fixture.basis, ducklips::ReprojectionInvariant::passive);
    double error = 0.0;
    for (std::size_t mode = 0; mode < fixture.basis.mode_count; ++mode) {
        error = std::max(error, std::abs(resonator.state[mode] - before[mode]));
    }
    check(error < 1.0e-9, "reprojection onto the same basis is the identity", error, 1.0e-9);
}

// Energy across a geometry change that closes the lips. Passive reprojection
// must not add energy; the two single-invariant modes are measured for the
// record.
double reprojected_energy(ducklips::ReprojectionInvariant invariant) {
    Fixture open{};
    Fixture closing{};
    ducklips::VowelShape vowel{};
    open.model.build_vowel(vowel, open.shape);
    vowel.lip_area = 0.004 * ducklips::square_centimetre;
    closing.model.build_vowel(vowel, closing.shape);
    for (std::size_t pass = 0; pass < 2; ++pass) {
        Fixture& fixture = pass == 0 ? open : closing;
        fixture.network.assemble(fixture.shape, ducklips::AirProperties{}, ducklips::NetworkLosses{});
        fixture.network.write_operator(fixture.matrix);
        fixture.solver.solve(fixture.network, fixture.matrix, 2.0 * ducklips::pi * 12000.0, fixture.empty, fixture.basis);
    }
    ducklips::ModalResonator resonator{};
    resonator.prepare(40, open.network.node_count(), open.network.branch_count());
    resonator.configure(open.network, open.basis, 1.0 / 48000.0, 0, 1);
    for (std::size_t mode = 0; mode < open.basis.mode_count; ++mode) {
        resonator.state[mode] = std::complex<double>(std::sin(3.0 + mode), std::cos(1.0 + 2.0 * mode));
    }
    double before = resonator.energy();
    resonator.reproject(open.network, open.basis, closing.network, closing.basis, invariant);
    resonator.configure(closing.network, closing.basis, 1.0 / 48000.0, 0, 1);
    double after = resonator.energy();
    double ratio = after / before;
    return ratio;
}

void test_passive_reprojection() {
    double passive = reprojected_energy(ducklips::ReprojectionInvariant::passive);
    double charge = reprojected_energy(ducklips::ReprojectionInvariant::charge_momentum);
    double pressure = reprojected_energy(ducklips::ReprojectionInvariant::pressure_flow);
    std::printf("    lip closure energy ratio: passive %.4f, charge/momentum %.4e, pressure/flow %.4e\n",
                passive, charge, pressure);
    check(passive <= 1.0 + 1.0e-12, "passive reprojection does not add energy at a closure", passive, 1.0);
}

void test_sensitivity() {
    ducklips::FormantInverter inverter{};
    inverter.prepare(ducklips::baritone_voice(), 64, 24, ducklips::AirProperties{});
    ducklips::VowelShape shape{};
    shape.coefficients = {2.9, 6.1, 0.35, 0.1, 0.2, 0.1, 0.0, -0.1, 0.05};
    ducklips::FormantArray frequencies{};
    ducklips::ShapeJacobian jacobian{};
    bool solved = inverter.formant_jacobian(shape, frequencies, jacobian);
    double worst = 0.0;
    const double step = 1.0e-5;
    for (std::size_t parameter = 0; parameter < ducklips::shape_parameter_count && solved; ++parameter) {
        ducklips::VowelShape plus = shape;
        ducklips::VowelShape minus = shape;
        plus.coefficients[parameter] += step;
        minus.coefficients[parameter] -= step;
        ducklips::FormantArray high{};
        ducklips::FormantArray low{};
        bool ok = inverter.formants(plus, high) && inverter.formants(minus, low);
        for (std::size_t formant = 0; formant < ducklips::fitted_formant_count && ok; ++formant) {
            double numeric = (std::log(high[formant]) - std::log(low[formant])) / (2.0 * step);
            double analytic = jacobian[formant][parameter];
            double scale = std::max(std::fabs(numeric), 0.05);
            worst = std::max(worst, std::fabs(analytic - numeric) / scale);
        }
    }
    check(solved && worst < 2.0e-3, "analytic formant sensitivities match finite differences", worst, 2.0e-3);
}

void test_sources() {
    ducklips::AirProperties air{};
    double area = 0.12 * ducklips::square_centimetre;
    double drop = 640.0;
    double impedance = 2.0e6;
    double flow = ducklips::orifice_flow(area, drop, impedance, air);
    double residual = flow - area * std::sqrt(2.0 * (drop - impedance * flow) / air.density);
    check(std::fabs(residual) < 1.0e-12, "orifice flow satisfies Bernoulli with series impedance", std::fabs(residual), 1.0e-12);

    ducklips::IntraoralPressure cavity{};
    double compliance = 1.0e-8;
    double step = 1.0 / 48000.0;
    double inflow = 3.0e-4;
    double outlet = 0.05 * ducklips::square_centimetre;
    double before = cavity.pressure();
    double after = cavity.advance(inflow, outlet, compliance, step, air);
    double balance = compliance * (after - before) / step - inflow + outlet * std::sqrt(2.0 * after / air.density);
    check(std::fabs(balance) < 1.0e-12, "implicit intraoral step satisfies its equation", std::fabs(balance), 1.0e-12);
    // Sealed moving cavity: ambient-equivalent air volume is conserved,
    // including elastic wall storage and changing air compliance.
    cavity.reset();const double bulk=air.density*air.sound_speed*air.sound_speed;
    const double initial_volume=5e-5,wall_compliance=1e-8;
    double old_volume=initial_volume,worst_mass_error=0;
    for (int i=1;i<=200;++i) {
        const double new_volume=initial_volume-1e-5*std::sin(ducklips::pi*i/200.0);
        const double old_c=old_volume/bulk+wall_compliance,new_c=new_volume/bulk+wall_compliance;
        const double pressure=cavity.advance_moving(0,0,old_c,new_c,new_volume-old_volume,step,air);
        worst_mass_error=std::max(worst_mass_error,std::abs(new_volume+new_c*pressure-initial_volume));
        if (i==100) { check(pressure>900,"tongue compression raises pressure in a sealed cavity",pressure,900); }
        old_volume=new_volume;
    }
    check(worst_mass_error<1e-16 && std::abs(cavity.pressure())<1e-9,
          "sealed compression/expansion conserves air and returns to ambient",worst_mass_error,1e-16);
    // Vented compression and expansion must satisfy the same mass balance,
    // including signed backflow during expansion.
    cavity.reset();double worst_flow_error=0;old_volume=initial_volume;
    for (int i=1;i<=200;++i) {
        const double new_volume=initial_volume+1e-5*std::sin(2*ducklips::pi*i/200.0);
        const double old_c=old_volume/bulk+wall_compliance,new_c=new_volume/bulk+wall_compliance;
        const double old_p=cavity.pressure();
        const double new_p=cavity.advance_moving(inflow,outlet,old_c,new_c,new_volume-old_volume,step,air);
        const double vent=std::copysign(outlet*std::sqrt(2*std::abs(new_p)/air.density),new_p);
        const double error=(new_c*new_p-old_c*old_p+new_volume-old_volume)/step-inflow+vent;
        worst_flow_error=std::max(worst_flow_error,std::abs(error));old_volume=new_volume;
    }
    check(worst_flow_error<1e-12,"moving vented cavity conserves signed mass flow",worst_flow_error,1e-12);
}

void test_cue_queue() {
    ducklips::CueQueue queue{};
    queue.prepare(4);
    std::size_t accepted = 0;
    for (std::size_t index = 0; index < 6; ++index) {
        ducklips::Cue cue{};
        cue.time = static_cast<ducklips::SampleCount>(index);
        if (queue.push(cue)) {
            accepted += 1;
        }
    }
    ducklips::Cue first{};
    bool has = queue.peek(first);
    check(accepted == 4 && has && first.time == 0, "cue queue holds its capacity in order", static_cast<double>(accepted), 4.0);
}

}

// Report (not a pass/fail check): eigenfrequencies and -3 dB bandwidths of the
// fitted /a/ posture, with the share of each loss mechanism.
void report_vowel_bandwidths() {
    ducklips::VoiceProfile voice = ducklips::baritone_voice();
    ducklips::VowelLibrary library{};
    library.build(voice, 64, 24, ducklips::AirProperties{});
    Fixture fixture{};
    fixture.model.build_vowel(library.shape(ducklips::Phone::AA), fixture.shape);
    fixture.network.assemble(fixture.shape, ducklips::AirProperties{}, voice.losses);
    fixture.network.write_operator(fixture.matrix);
    fixture.solver.solve(fixture.network, fixture.matrix, 2.0 * ducklips::pi * 12000.0, fixture.empty, fixture.basis);
    std::printf("    /a/ modes: frequency, bandwidth (viscous, wall, radiation, thermal) Hz\n");
    for (std::size_t mode = 0; mode < fixture.basis.mode_count; ++mode) {
        double oral_energy = 0.0;
        for (std::size_t node = 0; node < fixture.network.oral_sections(); ++node) {
            double pressure = fixture.basis.pressure(mode, node);
            oral_energy += fixture.network.nodes[node].compliance * pressure * pressure;
        }
        double w = fixture.basis.angular_frequency[mode];
        std::printf("      [oral energy %.2f, lip weight %.3e] ", oral_energy,
                    std::fabs(fixture.basis.pressure(mode, fixture.network.oral_sections() - 1)
                              * fixture.network.branches[fixture.network.lip_branch()].conductance));
        double viscous = 0.0;
        double wall = 0.0;
        double radiation = 0.0;
        double thermal = 0.0;
        for (std::size_t branch = 0; branch < fixture.network.branch_count(); ++branch) {
            const ducklips::AcousticBranch& link = fixture.network.branches[branch];
            if (link.conductance <= 0.0) {
                continue;
            }
            double flow = ducklips::branch_flow_shape(fixture.network, fixture.basis, mode, branch);
            viscous += 0.5 * link.viscous_coefficient * std::sqrt(w) * flow * flow;
            wall += 0.5 * link.resistance * flow * flow;
            radiation += 0.5 * link.radiation_coefficient * w * w * flow * flow;
        }
        for (std::size_t node = 0; node < fixture.network.node_count(); ++node) {
            double pressure = fixture.basis.pressure(mode, node);
            thermal += 0.5 * fixture.network.nodes[node].thermal_coefficient * std::sqrt(w) * pressure * pressure;
        }
        // An amplitude decay alpha gives a -3 dB bandwidth of alpha / pi Hz.
        std::printf("      %6.0f Hz  B %6.1f  (%5.1f, %5.1f, %5.1f, %5.1f)\n", w / (2.0 * ducklips::pi),
                    fixture.basis.damping[mode] / ducklips::pi, viscous / ducklips::pi, wall / ducklips::pi,
                    radiation / ducklips::pi, thermal / ducklips::pi);
    }
}

// Report: area functions (cm^2, every 4th section, glottis to lips) of fitted postures.
void report_area_functions() {
    ducklips::VoiceProfile voice = ducklips::baritone_voice();
    ducklips::VowelLibrary library{};
    library.build(voice, 64, 24, ducklips::AirProperties{});
    Fixture fixture{};
    const ducklips::Phone phones[3] = {ducklips::Phone::AX, ducklips::Phone::IY, ducklips::Phone::AA};
    const char* names[3] = {"AX", "IY", "AA"};
    for (std::size_t index = 0; index < 3; ++index) {
        const ducklips::VowelShape& shape = library.shape(phones[index]);
        fixture.model.build_vowel(shape, fixture.shape);
        std::printf("    %s posture (%.2f %.2f %.2f %.2f %.2f | %.2f %.2f %.2f %.2f):", names[index], shape.coefficients[0],
                    shape.coefficients[1], shape.coefficients[2], shape.coefficients[3], shape.coefficients[4],
                    shape.coefficients[5], shape.coefficients[6], shape.coefficients[7], shape.coefficients[8]);
        for (std::size_t section = 0; section < 64; section += 3) {
            std::printf(" %.1f", fixture.shape.oral_area[section] / ducklips::square_centimetre);
        }
        std::printf("\n");
    }
}

int main() {
    test_tree_operator();
    test_modal_solver_against_dense();
    test_uniform_tube();
    test_reprojection_identity();
    test_passive_reprojection();
    test_sensitivity();
    test_sources();
    test_cue_queue();
    report_vowel_bandwidths();
    report_area_functions();
    std::printf("%zu failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
