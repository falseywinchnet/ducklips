#pragma once
#include "ducklips/glottal_contour.hpp"
#include "ducklips/sources.hpp"
namespace ducklips {
struct ContourKineticState final {
    std::array<double,3> inverse_inertance{},resistance{},quadratic_loss{};
    // Coordinates: lower mass, upper mass, lower area, upper area.
    std::array<double,4> effort{};
    GlottalContourGeometry geometry{};
    double energy{};
};
struct ContourKineticTransfer final {
    std::array<double,3> flow{};
    std::array<double,4> effort{};
    std::array<std::array<double,2>,2> volume_gradient{}; // chamber, control area
    double energy_change{},chain_residual{};
};
// Positive-storage contour reference. Density is uniform within each chamber; the three
// dual paths join the moving volume centroids. Mass and shape derivatives
// belong to the same kinetic energy. Sealed caps have zero inverse inertia
// and require zero momentum. Both chambers must retain positive volume.
// This constitutive kernel does not perform sealing/reopening state transport.
class ContourKinetics final {
public:
    [[nodiscard]] bool configure(const FoldParameters&,const AirProperties&,double outlet_loss=1.) noexcept;
    [[nodiscard]] bool evaluate(const std::array<double,2>& mass,const std::array<double,2>& area,
        const std::array<double,3>& momentum,ContourKineticState&) const noexcept;
    [[nodiscard]] bool transfer(const std::array<double,2>& mass0,const std::array<double,2>& area0,
        const std::array<double,3>& momentum0,const std::array<double,2>& mass1,const std::array<double,2>& area1,
        const std::array<double,3>& momentum1,ContourKineticTransfer&) const noexcept;
private:
    [[nodiscard]] bool metric(const std::array<double,4>& first,const std::array<double,4>& last,
        std::size_t coordinate,std::array<double,3>& value,std::array<double,3>& slope,
        std::array<double,2>* volume_slope=nullptr) const noexcept;
    double density_{},viscosity_{},width_{},outlet_loss_{};
    std::array<double,2> thickness_{};
};
}
