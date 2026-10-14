// Writes the area functions (cm^2, glottis to lips, one per line) of fitted
// vowel postures from the DuckLips vowel library, for the surface network.
#include "ducklips/articulation.hpp"
#include "ducklips/vowel_library.hpp"
#include "ducklips/voice.hpp"
#include "ducklips/acoustic_network.hpp"
#include "ducklips/modal.hpp"
#include "ducklips/tree_operator.hpp"
#include <cmath>
#include <cstdio>
#include <string>

namespace {

ducklips::TractShape shape_without_sides(const ducklips::TractShape& shape, ducklips::TractModel& model) {
    ducklips::TractShape copy{};
    model.allocate(copy);
    copy.oral_area = shape.oral_area;
    copy.oral_length = shape.oral_length;
    copy.nasal_area = shape.nasal_area;
    copy.nasal_length = shape.nasal_length;
    copy.port_area = 0.0;
    return copy;
}

}

int main() {
    ducklips::VoiceProfile voice = ducklips::baritone_voice();
    ducklips::VowelLibrary library{};
    library.build(voice, 64, 24, ducklips::AirProperties{});
    ducklips::TractModel model{};
    model.prepare(voice, 64, 24);
    ducklips::TractShape shape{};
    model.allocate(shape);
    const ducklips::Phone phones[3] = {ducklips::Phone::AA, ducklips::Phone::IY, ducklips::Phone::UW};
    const char* names[3] = {"AA", "IY", "UW"};
    for (std::size_t index = 0; index < 3; ++index) {
        model.build_vowel(library.shape(phones[index]), shape);
        std::string path = std::string("out/areas_") + names[index] + ".txt";
        std::FILE* file = std::fopen(path.c_str(), "w");
        for (std::size_t section = 0; section < shape.oral_area.size(); ++section) {
            std::fprintf(file, "%.4f %.6f\n", shape.oral_area[section] / ducklips::square_centimetre, shape.oral_length[section]);
        }
        std::fclose(file);
        // Reference: the same area function as a rigid-walled tube without
        // side cavities, solved by the eigen engine.
        ducklips::VoiceProfile rigid_voice = voice;
        rigid_voice.side_cavities = false;
        ducklips::TractModel rigid_model{};
        rigid_model.prepare(rigid_voice, 64, 24);
        ducklips::AcousticNetwork network{};
        network.prepare(64, 24, rigid_model.port_section(), rigid_model.side_cavities());
        ducklips::NetworkLosses losses{};
        losses.wall_mass_per_area = 1.0e15;
        network.assemble(shape_without_sides(shape, rigid_model), ducklips::AirProperties{}, losses);
        ducklips::TreeOperator matrix{};
        matrix.prepare(network.node_count(), network.tree_edges(), rigid_model.port_section());
        network.write_operator(matrix);
        ducklips::ModalSolver solver{};
        solver.prepare(network.node_count());
        ducklips::ModalBasis basis{};
        basis.prepare(network.node_count(), 40);
        ducklips::ModalBasis empty{};
        empty.prepare(network.node_count(), 0);
        solver.solve(network, matrix, 2.0 * ducklips::pi * 6000.0, empty, basis);
        std::printf("%s rigid tube, eigen engine:", names[index]);
        for (std::size_t mode = 0; mode < basis.mode_count; ++mode) {
            double oral = 0.0;
            for (std::size_t node = 0; node < network.oral_sections(); ++node) {
                oral += network.nodes[node].compliance * basis.pressure(mode, node) * basis.pressure(mode, node);
            }
            if (oral > 0.5) {
                std::printf(" %.0f", basis.angular_frequency[mode] / (2.0 * ducklips::pi));
            }
        }
        std::printf("\n");
    }
    return 0;
}
