#pragma once
#include "ducklips/physics.hpp"
#include "ducklips/tree_operator.hpp"
#include <span>
#include <vector>

namespace ducklips {
inline constexpr std::size_t cavity_reservoir=static_cast<std::size_t>(-1);
struct CavityConnection final {
    std::size_t from{},to{cavity_reservoir};
};
// Quasi-steady passage: dp = R Q + loss * rho Q |Q| / (2 A^2).
// Forward/reverse coefficients allow contraction losses without multiplying
// a complete orifice loss by the number of tube cells.
// Area zero seals the passage exactly. R must be positive for an open port.
// reservoir_pressure is used only when the connection ends at a reservoir.
struct CavityPassage final {
    double area{},resistance{},reservoir_pressure{};
    double forward_loss{1.0},reverse_loss{1.0}; // multiplier of rho Q |Q| / (2 A^2); zero is purely viscous
};
struct CavityStepReport final {
    bool converged{false};
    std::size_t iterations{};
    double maximum_mass_residual{}; // m^3/s, per node
};
// Spatial mean-pressure reference. Each air space stores its own pressure;
// a closed passage cannot transmit mass or deformation-induced pressure.
// Backward Euler conserves C p + V per node. This is not an acoustic solver.
class CavityNetwork final {
public:
    // Internal connections may form cycles or disconnected components. A tree
    // keeps the direct solve; other graphs use that elimination only as a
    // preconditioner, retaining every passage in the actual pressure equation.
    void prepare(std::size_t nodes,std::span<const CavityConnection> connections);
    bool reset(std::span<const double> volumes,std::span<const double> compliances,
               std::span<const double> pressures) noexcept;
    // Geometry and passages are the end-of-step values; inflow is m^3/s.
    // Failure leaves stored state and the last accepted port flows unchanged.
    [[nodiscard]] CavityStepReport advance(std::span<const double> volumes,std::span<const double> compliances,
        std::span<const CavityPassage> passages,std::span<const double> inflow,double dt,const AirProperties& air) noexcept;
    [[nodiscard]] std::span<const double> pressures() const noexcept { return pressure_; }
    [[nodiscard]] std::span<const double> flows() const noexcept { return flow_; }
private:
    double evaluate(std::span<const double> pressure,std::span<const CavityPassage> passages,
                    const AirProperties& air,bool jacobian) noexcept;
    double potential_change(std::span<const CavityPassage> passages,const AirProperties& air) const noexcept;
    bool balanced(double tolerance) const noexcept;
    bool solve_jacobian() noexcept;
    void multiply_jacobian(std::span<const double> x,std::span<double> y) const noexcept;
    TreeOperator matrix_{};
    std::vector<CavityConnection> connections_{};
    std::vector<std::size_t> slots_{};
    std::vector<double> pressure_{},volume_{},compliance_{},flow_{};
    std::vector<double> candidate_{},trial_{},residual_{},gradient_{},delta_{},storage_{},rhs_{},working_flow_{},linearization_flow_{},roundoff_{};
    std::vector<double> slope_{},boundary_diagonal_{},cg_residual_{},cg_direction_{},cg_image_{},cg_preconditioned_{};
    double balance_residual_{},balance_roundoff_{};
    bool ready_{false};
    bool direct_tree_{false};
};
}
