#include "ducklips/tree_operator.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {

void TreeOperator::prepare(std::size_t node_count, std::span<const TreeEdge> edges, std::size_t root) {
    if (node_count == 0 || edges.size() + 1 != node_count || root >= node_count) {
        throw std::invalid_argument("TreeOperator requires a spanning tree");
    }
    std::vector<std::vector<std::size_t>> neighbours(node_count);
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        const TreeEdge& link = edges[edge];
        neighbours[link.first_node].push_back(link.second_node);
        neighbours[link.second_node].push_back(link.first_node);
    }
    parent_.assign(node_count, no_parent);
    std::vector<bool> visited(node_count, false);
    std::vector<std::size_t> breadth_first_order{};
    breadth_first_order.reserve(node_count);
    breadth_first_order.push_back(root);
    visited[root] = true;
    for (std::size_t position = 0; position < breadth_first_order.size(); ++position) {
        std::size_t node = breadth_first_order[position];
        for (std::size_t neighbour_index = 0; neighbour_index < neighbours[node].size(); ++neighbour_index) {
            std::size_t neighbour = neighbours[node][neighbour_index];
            if (!visited[neighbour]) {
                visited[neighbour] = true;
                parent_[neighbour] = node;
                breadth_first_order.push_back(neighbour);
            }
        }
    }
    if (breadth_first_order.size() != node_count) {
        throw std::invalid_argument("TreeOperator edges do not connect every node");
    }
    elimination_order_.assign(breadth_first_order.rbegin(), breadth_first_order.rend());
    edge_slot_.assign(edges.size(), 0);
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        const TreeEdge& link = edges[edge];
        if (parent_[link.first_node] == link.second_node) {
            edge_slot_[edge] = link.first_node;
        } else if (parent_[link.second_node] == link.first_node) {
            edge_slot_[edge] = link.second_node;
        } else {
            throw std::invalid_argument("TreeOperator edge is not a tree edge");
        }
    }
    diagonal.assign(node_count, 0.0);
    parent_coupling.assign(node_count, 0.0);
    pivot_.assign(node_count, 0.0);
    multiplier_.assign(node_count, 0.0);
}

std::size_t TreeOperator::node_count() const noexcept {
    return parent_.size();
}

std::size_t TreeOperator::edge_slot(std::size_t edge) const noexcept {
    return edge_slot_[edge];
}

std::size_t TreeOperator::parent(std::size_t node) const noexcept {
    return parent_[node];
}

std::size_t TreeOperator::factor(double shift) noexcept {
    const std::size_t count = parent_.size();
    double scale = 0.0;
    for (std::size_t node = 0; node < count; ++node) {
        pivot_[node] = diagonal[node] - shift;
        scale = std::max(scale, std::fabs(diagonal[node]));
    }
    // Replace only pivots that are zero to working precision of the largest entry.
    const double floor = scale * 1.0e-30 + 1.0e-300;
    std::size_t negative = 0;
    for (std::size_t position = 0; position < count; ++position) {
        std::size_t node = elimination_order_[position];
        double pivot = pivot_[node];
        if (std::fabs(pivot) < floor) {
            if (pivot < 0.0) {
                pivot = -floor;
            } else {
                pivot = floor;
            }
            pivot_[node] = pivot;
        }
        if (pivot < 0.0) {
            negative += 1;
        }
        std::size_t up = parent_[node];
        if (up != no_parent) {
            double coupling = parent_coupling[node];
            double multiplier = coupling / pivot;
            multiplier_[node] = multiplier;
            pivot_[up] -= multiplier * coupling;
        }
    }
    return negative;
}

std::size_t TreeOperator::count_below(double shift) noexcept {
    std::size_t negative = factor(shift);
    return negative;
}

void TreeOperator::solve(std::span<double> x) const noexcept {
    const std::size_t count = parent_.size();
    // L y = b, children first: a finished node's value is folded into its parent.
    for (std::size_t position = 0; position < count; ++position) {
        std::size_t node = elimination_order_[position];
        std::size_t up = parent_[node];
        if (up != no_parent) {
            x[up] -= multiplier_[node] * x[node];
        }
    }
    for (std::size_t node = 0; node < count; ++node) {
        x[node] /= pivot_[node];
    }
    // L^T x = z, root first.
    for (std::size_t position = count; position > 0; --position) {
        std::size_t node = elimination_order_[position - 1];
        std::size_t up = parent_[node];
        if (up != no_parent) {
            x[node] -= multiplier_[node] * x[up];
        }
    }
}

void TreeOperator::multiply(std::span<const double> x, std::span<double> y) const noexcept {
    const std::size_t count = parent_.size();
    for (std::size_t node = 0; node < count; ++node) {
        y[node] = diagonal[node] * x[node];
    }
    for (std::size_t node = 0; node < count; ++node) {
        std::size_t up = parent_[node];
        if (up != no_parent) {
            double coupling = parent_coupling[node];
            y[node] += coupling * x[up];
            y[up] += coupling * x[node];
        }
    }
}

}
