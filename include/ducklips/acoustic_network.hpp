#pragma once
#include "ducklips/physics.hpp"
#include "ducklips/tree_operator.hpp"
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {

// Geometry of the airway at one instant. Sizes are fixed when the owner is
// prepared; every build writes all elements.
struct TractShape final {
    std::vector<double> oral_area{};    // m^2 per section; [0] at the glottis, back() at the lips
    std::vector<double> oral_length{};  // m per section
    std::vector<double> nasal_area{};   // m^2; [0] at the velar port, back() at the nostrils
    std::vector<double> nasal_length{}; // m
    double port_area{0.0};              // m^2 velopharyngeal opening; 0 closes the nasal branch
    double port_length{0.01};           // m
    // Opening at the mouth plane, distinct from the air volume immediately
    // behind it. Negative retains the uniform last-section convention.
    double lip_aperture{-1.0};         // m^2; zero is actual lip contact
    // Side cavities (piriform fossae, paranasal sinuses): every section of
    // every cavity, in the order of the network's SideCavity list.
    std::vector<double> side_area{};    // m^2
    std::vector<double> side_length{};  // m
    // The space under the tongue (TractModel::set_floor), per oral section:
    // its cross-section (m^2, 0 where the tongue rests on the floor) and the
    // face (m^2) through which it opens forward into the mouth, at the last
    // section under the tip. Zero without a floor; the modal engine ignores
    // them.
    std::vector<double> sublingual_area{};
    std::vector<double> sublingual_opening{};
};

// A closed side branch hanging off one node of the airway: a chain of
// sections whose last section is closed. It adds a leaf to the tree.
struct SideCavity final {
    std::size_t attach_node{};          // oral or nasal node it opens from
    std::size_t sections{};
    double perimeter_factor{1.0};

};

// Loss and boundary constants. Wall values follow Ishizaka, French and
// Flanagan's cheek measurement (1.5 g/cm^2, about 1000 dyn s/cm^3) converted to
// SI; the perimeter factors and the radiation factor are modelling choices,
// not measurements.
struct NetworkLosses final {
    double wall_mass_per_area{15.0};         // kg/m^2
    double wall_resistance_per_area{1.0e4};  // Pa s/m
    double oral_perimeter_factor{1.3};       // elliptical and irregular sections versus a circle
    double nasal_perimeter_factor{3.0};      // turbinate surface
    double radiation_factor{1.4};            // head baffling between a free monopole (1) and a baffle (2)
    double end_correction{0.7};              // radiation end correction / opening radius
    double minimum_area{1.0e-9};             // m^2 floor for section compliance
    // Above the first cross-mode cutoff the plane-wave tube no longer holds:
    // energy scatters into 3-D modes of the irregular tract. The 1-D model
    // would instead trap its highest modes in narrow sections as high-Q
    // whistles. Added bandwidth = cross_mode_bandwidth * ((f - onset)/onset)^2
    // above the onset. A modelling choice standing in for missing 3-D physics.
    double cross_mode_onset{4500.0};         // Hz
    double cross_mode_bandwidth{800.0};      // Hz at twice the onset
};

enum class BranchKind { tract, lip_radiation, velar_port, nasal, nostril_radiation, wall, side };

inline constexpr std::size_t ground_node = static_cast<std::size_t>(-1);

// An inertive air path between two nodes, or between a node and still air.
// conductance is the reciprocal inertance 1/M (m^4/kg); a closed path has 0.
// Resistance at angular frequency w is
//   viscous_coefficient * sqrt(w) + resistance + radiation_coefficient * w^2
// in Pa s/m^3.
struct AcousticBranch final {
    std::size_t from{};
    std::size_t to{ground_node};
    BranchKind kind{BranchKind::tract};
    double conductance{};
    double viscous_coefficient{};
    double resistance{};
    double radiation_coefficient{};
    double aperture{-1.0};             // m^2 at an exterior boundary, independent of node storage
};

// The compressible air volume of one tube section. Its shunt thermal
// conductance at w is thermal_coefficient * sqrt(w) (m^3/(Pa s)).
struct AcousticNode final {
    double compliance{};          // m^3/Pa
    double thermal_coefficient{};
    double area{};                // m^2 after the compliance floor
    double wall_area{};           // m^2 of yielding surface
    double compliance_rate{};     // m^3/(Pa s), tangent of moving local geometry
};

// The airway as a lumped acoustic network: oral sections 0..oral-1 from the
// glottis to the lips, nasal sections oral..oral+nasal-1 from the velar port
// to the nostrils, then the sections of each side cavity. Pressures live on
// nodes, volume flows on branches.
//
// Branch layout (fixed):
//   [0, oral-1)                     oral tract, section i to i+1
//   oral-1                          lips to still air
//   oral                            velar port, port section to first nasal section
//   [oral+1, oral+nasal)            nasal tract
//   oral+nasal                      nostrils to still air
//   oral+nasal+1+node               yielding wall of each node (all nodes)
//   after the walls                 side-cavity chains, one branch per section
// The tract, port, nasal and side branches form a tree; radiation and wall
// branches only add to the diagonal of the stiffness operator.
class AcousticNetwork final {
public:
    void prepare(std::size_t oral_sections, std::size_t nasal_sections, std::size_t port_section,
                 std::span<const SideCavity> side_cavities);
    // Writes every node and branch from the shape. No allocation.
    void assemble(const TractShape& shape, const AirProperties& air, const NetworkLosses& losses) noexcept;
    // Writes S = C^-1/2 K C^-1/2, K = G M^-1 G^T, into a TreeOperator prepared
    // with tree_edges() rooted at port_section().
    void write_operator(TreeOperator& matrix) const noexcept;

    [[nodiscard]] std::span<const TreeEdge> tree_edges() const noexcept;
    [[nodiscard]] std::size_t node_count() const noexcept;
    [[nodiscard]] std::size_t branch_count() const noexcept;
    [[nodiscard]] std::size_t oral_sections() const noexcept;
    [[nodiscard]] std::size_t nasal_sections() const noexcept;
    [[nodiscard]] std::size_t port_section() const noexcept;
    [[nodiscard]] std::size_t oral_branch(std::size_t section) const noexcept;
    [[nodiscard]] std::size_t lip_branch() const noexcept;
    [[nodiscard]] std::size_t port_branch() const noexcept;
    [[nodiscard]] std::size_t nostril_branch() const noexcept;
    [[nodiscard]] std::size_t nostril_node() const noexcept;
    [[nodiscard]] std::size_t side_node_count() const noexcept;
    [[nodiscard]] std::size_t wall_branch(std::size_t node) const noexcept;

    std::vector<AcousticNode> nodes{};
    std::vector<AcousticBranch> branches{};
    NetworkLosses losses{}; // as of the last assemble()

private:
    std::size_t oral_sections_{};
    std::size_t nasal_sections_{};
    std::size_t port_section_{};
    std::vector<SideCavity> side_cavities_{};
    std::size_t side_nodes_{};
    std::vector<TreeEdge> tree_edges_{};
    std::vector<std::size_t> edge_branch_{}; // branch index of each tree edge
};

}
