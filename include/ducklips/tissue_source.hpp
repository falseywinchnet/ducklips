#pragma once
#include <array>
#include <cstddef>

namespace ducklips {

// Experimental lateral-displacement continuum. Two elastic layers, attached
// at the anterior/posterior ends; free inferior/superior edges. Material
// constants and dimensions determine the motion. No prescribed period,
// waveform, temporal phase, pitch servo, or random excitation.
struct TissueMaterial final {
    double length{.015};
    double height{.004};
    double cover_depth{.001};
    double body_depth{.004};
    double density{1030};
    // A paired human in-vivo Voigt estimate (subject 2, PMC3961822).
    // This grounds the cover's starting material, not the rest of the body.
    double cover_shear{878};             // Pa, transverse shear modulus
    double body_shear{1200};
    double longitudinal_stress{6000};    // Pa, prestress
    double body_longitudinal_stress{10000};
    double viscosity{.52};              // Pa s, Kelvin-Voigt shear viscosity
    double body_viscosity{.16};
    double rest_gap{.00010};             // full central gap, m
    double convergence{.02};            // radians per medial surface
    double contact_depth{.0003};         // compliant contact layer depth, m
    double air_density{1.14};
    double air_viscosity{1.86e-5};
};

class TissueSource final {
public:
    static constexpr std::size_t nx=13,nz=13,count=nx*nz;
    using Field=std::array<double,count>;
    void reset(const TissueMaterial& material={});
    // Muscle posture moves the elastic equilibrium. Absolute surface/body
    // positions and their velocities remain continuous across the change.
    void set_adduction(double central_rest_gap);
    // Quasi-steady slit airflow and distributed surface pressure. The flow
    // approximation is explicit and provisional; tissue motion is dynamic.
    double advance(double lung_pressure,double outlet_pressure,double dt);
    [[nodiscard]] const Field& surface() const noexcept{return cover_;}
    [[nodiscard]] const Field& body() const noexcept{return body_;}
    [[nodiscard]] const Field& gap() const noexcept{return gap_;}
    [[nodiscard]] const Field& pressure() const noexcept{return pressure_;}
    [[nodiscard]] const Field& contact() const noexcept{return contact_;}
    [[nodiscard]] const Field& velocity() const noexcept{return cover_v_;}
    [[nodiscard]] const std::array<double,nz>& area() const noexcept{return area_;}
    [[nodiscard]] double flow() const noexcept{return flow_;}
    [[nodiscard]] double minimum_area() const noexcept{return minimum_area_;}
    [[nodiscard]] double kinetic_energy() const noexcept;
    [[nodiscard]] double elastic_energy() const noexcept;
private:
    void geometry();
    void aerodynamics(double lung,double outlet);
    TissueMaterial material_{};
    Field cover_{},body_{},cover_v_{},body_v_{},gap_{},rest_{},pressure_{},contact_{};
    Field cover_a_{},body_a_{};
    // Rotation about the longitudinal axis, located between vertical nodes.
    // The last z row is unused. Rotation is distinct from shear strain.
    Field rotation_{},rotation_v_{},rotation_a_{},shear_{};
    std::array<double,nz> area_{},resistance_{};
    double dx_{},dz_{},flow_{},minimum_area_{};
};
}
