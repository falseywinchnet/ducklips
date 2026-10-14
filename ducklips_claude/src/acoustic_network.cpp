#include "ducklips/acoustic_network.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {

// Areas below this are treated as a closed passage.
constexpr double closed_area = 1.0e-12; // m^2

double circular_perimeter(double area) noexcept {
    double radius = std::sqrt(area / pi);
    double perimeter = 2.0 * pi * radius;
    return perimeter;
}

struct TubeHalf final {
    double area{};
    double length{};
    double perimeter_factor{};
};

// A branch through half of each adjoining section.
void write_tube_branch(AcousticBranch& branch, std::size_t from, std::size_t to, BranchKind kind,
                       const TubeHalf& first, const TubeHalf& second,
                       const AirProperties& air, double viscous_scale) noexcept {
    branch.from = from;
    branch.to = to;
    branch.kind = kind;
    branch.resistance = 0.0;
    branch.radiation_coefficient = 0.0;
    if (first.area <= closed_area || second.area <= closed_area) {
        branch.conductance = 0.0;
        branch.viscous_coefficient = 0.0;
        return;
    }
    double inertance = air.density * (first.length / (2.0 * first.area) + second.length / (2.0 * second.area));
    branch.conductance = 1.0 / inertance;
    double first_wall = first.perimeter_factor * circular_perimeter(first.area) * first.length / 2.0;
    double second_wall = second.perimeter_factor * circular_perimeter(second.area) * second.length / 2.0;
    double first_term = first_wall / (first.area * first.area);
    double second_term = second_wall / (second.area * second.area);
    branch.viscous_coefficient = viscous_scale * (first_term + second_term);
}

// A branch from the last half section through the radiation end correction
// into still air. Radiation resistance is the monopole value rho w^2/(4 pi c)
// scaled by the baffle factor.
void write_radiation_branch(AcousticBranch& branch, std::size_t from, BranchKind kind, const TubeHalf& half,
                            const AirProperties& air, const NetworkLosses& losses, double viscous_scale) noexcept {
    branch.from = from;
    branch.to = ground_node;
    branch.kind = kind;
    branch.resistance = 0.0;
    if (half.area <= closed_area) {
        branch.conductance = 0.0;
        branch.viscous_coefficient = 0.0;
        branch.radiation_coefficient = 0.0;
        return;
    }
    double radius = std::sqrt(half.area / pi);
    double end_length = losses.end_correction * radius;
    double inertance = air.density * (half.length / (2.0 * half.area) + end_length / half.area);
    branch.conductance = 1.0 / inertance;
    double wall = half.perimeter_factor * circular_perimeter(half.area) * half.length / 2.0;
    branch.viscous_coefficient = viscous_scale * wall / (half.area * half.area);
    branch.radiation_coefficient = losses.radiation_factor * air.density / (4.0 * pi * air.sound_speed);
}

void write_node(AcousticNode& node, double raw_area, double length, double perimeter_factor,
                const AirProperties& air, const NetworkLosses& losses) noexcept {
    const double bulk_stiffness = air.density * air.sound_speed * air.sound_speed;
    const double thermal_scale = (air.heat_capacity_ratio - 1.0) / bulk_stiffness * std::sqrt(air.thermal_diffusivity / 2.0);
    double area = std::max(raw_area, losses.minimum_area);
    double wall = perimeter_factor * circular_perimeter(area) * length;
    node.area = area;
    node.compliance = area * length / bulk_stiffness;
    node.wall_area = wall;
    node.thermal_coefficient = thermal_scale * wall;
}

}

void AcousticNetwork::prepare(std::size_t oral_sections, std::size_t nasal_sections, std::size_t port_section,
                              std::span<const SideCavity> side_cavities) {
    if (oral_sections < 4 || nasal_sections < 2 || port_section + 1 >= oral_sections) {
        throw std::invalid_argument("AcousticNetwork section counts are too small");
    }
    oral_sections_ = oral_sections;
    nasal_sections_ = nasal_sections;
    port_section_ = port_section;
    side_cavities_.assign(side_cavities.begin(), side_cavities.end());
    side_nodes_ = 0;
    for (std::size_t cavity = 0; cavity < side_cavities_.size(); ++cavity) {
        if (side_cavities_[cavity].sections == 0 || side_cavities_[cavity].attach_node >= oral_sections + nasal_sections) {
            throw std::invalid_argument("AcousticNetwork side cavity is invalid");
        }
        side_nodes_ += side_cavities_[cavity].sections;
    }
    std::size_t node_total = oral_sections + nasal_sections + side_nodes_;
    nodes.assign(node_total, AcousticNode{});
    branches.assign(2 * node_total + 1, AcousticBranch{});
    tree_edges_.clear();
    edge_branch_.clear();
    for (std::size_t section = 0; section + 1 < oral_sections; ++section) {
        tree_edges_.push_back(TreeEdge{.first_node = section, .second_node = section + 1});
        edge_branch_.push_back(oral_branch(section));
    }
    tree_edges_.push_back(TreeEdge{.first_node = port_section, .second_node = oral_sections});
    edge_branch_.push_back(port_branch());
    for (std::size_t section = 0; section + 1 < nasal_sections; ++section) {
        tree_edges_.push_back(TreeEdge{.first_node = oral_sections + section, .second_node = oral_sections + section + 1});
        edge_branch_.push_back(oral_sections + 1 + section);
    }
    // Side chains: their branches follow the walls; one branch enters each section.
    std::size_t node = oral_sections + nasal_sections;
    std::size_t branch = oral_sections + nasal_sections + 1 + node_total;
    for (std::size_t cavity = 0; cavity < side_cavities_.size(); ++cavity) {
        std::size_t upstream = side_cavities_[cavity].attach_node;
        for (std::size_t section = 0; section < side_cavities_[cavity].sections; ++section) {
            tree_edges_.push_back(TreeEdge{.first_node = upstream, .second_node = node});
            edge_branch_.push_back(branch);
            upstream = node;
            node += 1;
            branch += 1;
        }
    }
}

void AcousticNetwork::assemble(const TractShape& shape, const AirProperties& air, const NetworkLosses& wall_losses) noexcept {
    losses = wall_losses;
    const double viscous_scale = std::sqrt(air.density * air.viscosity / 2.0);
    for (std::size_t section = 0; section < oral_sections_; ++section) {
        write_node(nodes[section], shape.oral_area[section], shape.oral_length[section],
                   losses.oral_perimeter_factor, air, losses);
    }
    for (std::size_t section = 0; section < nasal_sections_; ++section) {
        write_node(nodes[oral_sections_ + section], shape.nasal_area[section], shape.nasal_length[section],
                   losses.nasal_perimeter_factor, air, losses);
    }
    for (std::size_t section = 0; section + 1 < oral_sections_; ++section) {
        TubeHalf first{.area = shape.oral_area[section], .length = shape.oral_length[section],
                       .perimeter_factor = losses.oral_perimeter_factor};
        TubeHalf second{.area = shape.oral_area[section + 1], .length = shape.oral_length[section + 1],
                        .perimeter_factor = losses.oral_perimeter_factor};
        write_tube_branch(branches[oral_branch(section)], section, section + 1, BranchKind::tract,
                          first, second, air, viscous_scale);
    }
    std::size_t lip_section = oral_sections_ - 1;
    TubeHalf lips{.area = shape.oral_area[lip_section], .length = shape.oral_length[lip_section],
                  .perimeter_factor = losses.oral_perimeter_factor};
    write_radiation_branch(branches[lip_branch()], lip_section, BranchKind::lip_radiation, lips, air, losses, viscous_scale);

    // Velar port: a short tube of its own plus half of each adjoining section.
    AcousticBranch& port = branches[port_branch()];
    port.from = port_section_;
    port.to = oral_sections_;
    port.kind = BranchKind::velar_port;
    port.resistance = 0.0;
    port.radiation_coefficient = 0.0;
    double velum_area = shape.oral_area[port_section_];
    double nasal_entry_area = shape.nasal_area[0];
    if (shape.port_area <= closed_area || velum_area <= closed_area || nasal_entry_area <= closed_area) {
        port.conductance = 0.0;
        port.viscous_coefficient = 0.0;
    } else {
        double port_term = shape.port_length / shape.port_area;
        double velum_term = shape.oral_length[port_section_] / (2.0 * velum_area);
        double nasal_term = shape.nasal_length[0] / (2.0 * nasal_entry_area);
        double inertance = air.density * (port_term + velum_term + nasal_term);
        port.conductance = 1.0 / inertance;
        double port_wall = circular_perimeter(shape.port_area) * shape.port_length;
        port.viscous_coefficient = viscous_scale * port_wall / (shape.port_area * shape.port_area);
    }

    for (std::size_t section = 0; section + 1 < nasal_sections_; ++section) {
        TubeHalf first{.area = shape.nasal_area[section], .length = shape.nasal_length[section],
                       .perimeter_factor = losses.nasal_perimeter_factor};
        TubeHalf second{.area = shape.nasal_area[section + 1], .length = shape.nasal_length[section + 1],
                        .perimeter_factor = losses.nasal_perimeter_factor};
        write_tube_branch(branches[oral_sections_ + 1 + section], oral_sections_ + section, oral_sections_ + section + 1,
                          BranchKind::nasal, first, second, air, viscous_scale);
    }
    std::size_t nostril_section = nasal_sections_ - 1;
    TubeHalf nostrils{.area = shape.nasal_area[nostril_section], .length = shape.nasal_length[nostril_section],
                      .perimeter_factor = losses.nasal_perimeter_factor};
    write_radiation_branch(branches[nostril_branch()], oral_sections_ + nostril_section, BranchKind::nostril_radiation,
                           nostrils, air, losses, viscous_scale);

    // Side cavities: a chain from the attach node; the first branch takes only
    // the half section inside the cavity (the opening is the attach node's).
    std::size_t side_node = oral_sections_ + nasal_sections_;
    std::size_t side_branch = oral_sections_ + nasal_sections_ + 1 + nodes.size();
    std::size_t side_section = 0;
    for (std::size_t cavity = 0; cavity < side_cavities_.size(); ++cavity) {
        const SideCavity& spec = side_cavities_[cavity];
        std::size_t upstream = spec.attach_node;
        TubeHalf previous{.area = nodes[spec.attach_node].area, .length = 0.0, .perimeter_factor = spec.perimeter_factor};
        for (std::size_t section = 0; section < spec.sections; ++section) {
            double area = shape.side_area[side_section];
            double length = shape.side_length[side_section];
            write_node(nodes[side_node], area, length, spec.perimeter_factor, air, losses);
            TubeHalf current{.area = area, .length = length, .perimeter_factor = spec.perimeter_factor};
            write_tube_branch(branches[side_branch], upstream, side_node, BranchKind::side, previous, current, air,
                              viscous_scale);
            previous = current;
            upstream = side_node;
            side_node += 1;
            side_branch += 1;
            side_section += 1;
        }
    }

    // Yielding walls: a mass and resistance per unit area, driven by the node pressure.
    for (std::size_t node_index = 0; node_index < nodes.size(); ++node_index) {
        AcousticBranch& wall = branches[wall_branch(node_index)];
        double wall_area = nodes[node_index].wall_area;
        wall.from = node_index;
        wall.to = ground_node;
        wall.kind = BranchKind::wall;
        wall.conductance = wall_area / losses.wall_mass_per_area;
        wall.viscous_coefficient = 0.0;
        wall.resistance = losses.wall_resistance_per_area / wall_area;
        wall.radiation_coefficient = 0.0;
    }
}

void AcousticNetwork::write_operator(TreeOperator& matrix) const noexcept {
    std::size_t node_total = nodes.size();
    for (std::size_t node_index = 0; node_index < node_total; ++node_index) {
        matrix.diagonal[node_index] = 0.0;
    }
    for (std::size_t branch_index = 0; branch_index < branches.size(); ++branch_index) {
        const AcousticBranch& branch = branches[branch_index];
        matrix.diagonal[branch.from] += branch.conductance;
        if (branch.to != ground_node) {
            matrix.diagonal[branch.to] += branch.conductance;
        }
    }
    for (std::size_t node_index = 0; node_index < node_total; ++node_index) {
        matrix.diagonal[node_index] /= nodes[node_index].compliance;
    }
    for (std::size_t edge = 0; edge < tree_edges_.size(); ++edge) {
        const AcousticBranch& branch = branches[edge_branch_[edge]];
        double scale = std::sqrt(nodes[branch.from].compliance * nodes[branch.to].compliance);
        matrix.parent_coupling[matrix.edge_slot(edge)] = -branch.conductance / scale;
    }
}

std::span<const TreeEdge> AcousticNetwork::tree_edges() const noexcept {
    return std::span<const TreeEdge>(tree_edges_);
}

std::size_t AcousticNetwork::node_count() const noexcept {
    return nodes.size();
}

std::size_t AcousticNetwork::branch_count() const noexcept {
    return branches.size();
}

std::size_t AcousticNetwork::oral_sections() const noexcept {
    return oral_sections_;
}

std::size_t AcousticNetwork::nasal_sections() const noexcept {
    return nasal_sections_;
}

std::size_t AcousticNetwork::port_section() const noexcept {
    return port_section_;
}

std::size_t AcousticNetwork::oral_branch(std::size_t section) const noexcept {
    return section;
}

std::size_t AcousticNetwork::lip_branch() const noexcept {
    return oral_sections_ - 1;
}

std::size_t AcousticNetwork::port_branch() const noexcept {
    return oral_sections_;
}

std::size_t AcousticNetwork::nostril_branch() const noexcept {
    return oral_sections_ + nasal_sections_;
}

std::size_t AcousticNetwork::nostril_node() const noexcept {
    return oral_sections_ + nasal_sections_ - 1;
}

std::size_t AcousticNetwork::side_node_count() const noexcept {
    return side_nodes_;
}

std::size_t AcousticNetwork::wall_branch(std::size_t node) const noexcept {
    return oral_sections_ + nasal_sections_ + 1 + node;
}

}
