#pragma once
#include <cstddef>
#include <span>
#include <vector>

namespace ducklips {

struct TreeEdge final {
    std::size_t first_node{};
    std::size_t second_node{};
};

// Symmetric real matrix S whose off-diagonal sparsity pattern is a tree.
//
// Eliminating nodes leaves-first produces no fill, so a shifted LDL^T
// factorization, its inertia and one solve each cost O(node_count). By
// Sylvester's law of inertia the number of negative pivots of S - shift I
// equals the number of eigenvalues of S below the shift. The modal solver uses
// that count to place each eigenvalue by its index.
//
// This is the same elimination the wrench transport solver performs over a
// rigid-body contact graph: eliminating a node folds its load into its parent.
class TreeOperator final {
public:
    static constexpr std::size_t no_parent = static_cast<std::size_t>(-1);

    // Allocates. edges must form a spanning tree over node_count nodes.
    void prepare(std::size_t node_count, std::span<const TreeEdge> edges, std::size_t root);

    [[nodiscard]] std::size_t node_count() const noexcept;
    // The child node of an edge; parent_coupling[child] holds S(child, parent).
    [[nodiscard]] std::size_t edge_slot(std::size_t edge) const noexcept;
    [[nodiscard]] std::size_t parent(std::size_t node) const noexcept;

    // Matrix values. The owner writes them between factorizations; their sizes
    // are fixed by prepare(). diagonal[i] = S(i, i); parent_coupling[i] =
    // S(i, parent(i)); the root's coupling slot is unused.
    std::vector<double> diagonal{};
    std::vector<double> parent_coupling{};

    // Factors S - shift I and returns the number of negative pivots. A pivot
    // smaller than 1e-30 of the largest diagonal magnitude is replaced by that
    // floor with its sign kept, the usual Sturm-count convention.
    std::size_t factor(double shift) noexcept;
    [[nodiscard]] std::size_t count_below(double shift) noexcept;
    // Solves (S - shift I) x = b in place with the most recent factorization.
    void solve(std::span<double> right_side_in_solution_out) const noexcept;
    // y = S x. x and y must not overlap.
    void multiply(std::span<const double> x, std::span<double> y) const noexcept;

private:
    std::vector<std::size_t> elimination_order_{}; // children before parents, root last
    std::vector<std::size_t> parent_{};
    std::vector<std::size_t> edge_slot_{};
    std::vector<double> pivot_{};
    std::vector<double> multiplier_{}; // L(parent, node) of each non-root node
};

}
