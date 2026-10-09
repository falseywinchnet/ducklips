#pragma once
#include "ducklips/tongue.hpp"
#include <array>
#include <span>

namespace ducklips {

// A material map, not a fitted speech trajectory. The stream function
// psi(x,y)=-integral(lift(x) dx)*gate(y) gives a divergence-free velocity.
// Implicit midpoint preserves its planar area; keeping z unchanged also
// preserves 3D volume. The command interval is dimensionless, not elapsed
// physiological time. This map supplies kinematics, not muscle dynamics.
struct TongueLiftKnot final {
    double x{},lift{},slope{}; // cm, cm per command unit, slope per command unit
};
struct TongueMaterialPoint final {
    Point2 position{};
    std::array<double,4> gradient{1,0,0,1}; // row-major material deformation gradient
    bool valid{false};
};
class ConservativeTongueMap final {
public:
    static constexpr std::size_t maximum_knots=513;
    // Knots must increase in x. The profile must meet zero lift and slope at
    // both ends. The caller is responsible for the wall/contact profile.
    // Lower tissue stays fixed at/below anchor_y; the gate is one
    // at/above free_y. The primitive is zero at centre_x. No allocation.
    bool configure(std::span<const TongueLiftKnot> knots,double anchor_y,double free_y,double centre_x) noexcept;
    [[nodiscard]] TongueMaterialPoint map(Point2 point,double command,std::size_t steps=32) const noexcept;
    [[nodiscard]] double lift(double x) const noexcept;
private:
    struct Piece final { double x{},width{},a{},b{},c{},d{},integral{}; };
    struct Field final { Point2 velocity{};std::array<double,4> derivative{}; };
    [[nodiscard]] Field field(Point2 point) const noexcept;
    void profile(double x,double& value,double& derivative,double& integral) const noexcept;
    std::array<Piece,maximum_knots-1> pieces_{};
    std::size_t count_{};
    double anchor_y_{},height_{},primitive_offset_{},total_integral_{};
};
}
