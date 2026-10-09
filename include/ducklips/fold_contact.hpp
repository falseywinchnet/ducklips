#pragma once
#include "ducklips/sources.hpp"
namespace ducklips {
struct FoldContactState final {
    double energy{};
    std::array<double,2> area_effort{};
};
struct FoldContactTransfer final {
    double energy_change{},chain_residual{};
    std::array<double,2> area_effort{};
};
// Two unilateral elastic contacts, one at each tissue control. The paired
// stiffness is the series stiffness of the opposing tissues, scaled by the
// body's collision_ratio. Contact forces on the two sides are equal. This
// is a stated lumped elastic law, not a measured collision/adhesion model.
class FoldContact final {
public:
    [[nodiscard]] bool configure(const FoldParameters&,double tension) noexcept;
    [[nodiscard]] bool evaluate(const std::array<double,2>& signed_area,FoldContactState&) const noexcept;
    [[nodiscard]] bool transfer(const std::array<double,2>& first,const std::array<double,2>& last,FoldContactTransfer&) const noexcept;
private:
    std::array<double,2> area_stiffness_{}; // paired N/m divided by width^2
    bool configured_{};
};
}
