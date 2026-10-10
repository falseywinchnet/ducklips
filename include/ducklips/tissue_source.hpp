#pragma once
#include <array>
#include <cstddef>
#include "ducklips/tissue_body_data.hpp"
#include "ducklips/tissue_relaxation.hpp"

namespace ducklips {

enum class TissueIdea { reference, shared_relaxation, contact_recovery, strain_hardening, axial_airway, measured_posture };

// Experimental reduced elastic volume. The body compiler derives coupled
// three-dimensional displacement from geometry, density and material law.
// Runtime pressure and contact evolve the deformation coordinates. No
// prescribed period, waveform, temporal phase, pitch servo or random kick.
// Recompile data/tissue_body.json to change the elastic body; rest posture,
// contact and airflow parameters remain runtime controls.
struct TissueMaterial final {
    double length{compiled_tissue::length};
    double height{compiled_tissue::height};
    double cover_depth{compiled_tissue::cover_depth};
    double body_depth{compiled_tissue::body_depth};
    double density{compiled_tissue::density};
    // A paired human in-vivo Voigt estimate (subject 2, PMC3961822).
    // This grounds the cover's starting material, not the rest of the body.
    double cover_shear{compiled_tissue::cover_shear}; // Pa
    double body_shear{compiled_tissue::body_shear};
    double longitudinal_stress{compiled_tissue::longitudinal_stress}; // Pa
    double body_longitudinal_stress{compiled_tissue::body_longitudinal_stress};
    double viscosity{compiled_tissue::cover_viscosity}; // Pa s
    double body_viscosity{compiled_tissue::body_viscosity};
    double rest_gap{.00010};             // full central gap, m
    double convergence{.02};            // radians per medial surface
    double contact_depth{.0003};         // compliant contact layer depth, m
    double relaxation_fraction{.1};
    double relaxation_seconds{.1};
    double relaxation_log_rate_sd{.6}; // experimental stochastic rate spread
    std::uint64_t relaxation_seed{1};
    bool shared_relaxation{false};
    bool contact_recovery{false};
    double hardening_strain{0}; // >0 adds a positive quartic material energy
    double air_density{1.14};
    double air_viscosity{1.86e-5};
};

class TissueSource final {
public:
    static constexpr std::size_t nx=compiled_tissue::nx,nz=compiled_tissue::nz,count=nx*nz;
    static constexpr std::size_t modes=compiled_tissue::modes;
    using Field=std::array<double,count>;
    void reset(const TissueMaterial& material={});
    // Muscle posture moves the elastic equilibrium. Absolute surface/body
    // positions and their velocities remain continuous across the change.
    void set_adduction(double central_rest_gap);
    // Active generalized closing force. Positive force pulls the passive
    // body inward; it supplies work through displacement, not teleportation.
    void set_holding_force(double newtons){holding_force_=newtons;}
    [[nodiscard]] static double posture_stiffness() noexcept;
    [[nodiscard]] double holding_work() const noexcept{return holding_work_;}
    // Quasi-steady slit airflow and distributed surface pressure. The flow
    // approximation is explicit and provisional; tissue motion is dynamic.
    double advance(double lung_pressure,double outlet_pressure,double dt);
    void set_tension_scale(double value) noexcept{tension_scale_=value;}
    [[nodiscard]] const Field& surface() const noexcept{return cover_;}
    [[nodiscard]] const Field& body() const noexcept{return body_;}
    [[nodiscard]] const Field& gap() const noexcept{return gap_;}
    [[nodiscard]] const Field& pressure() const noexcept{return pressure_;}
    [[nodiscard]] const Field& contact() const noexcept{return contact_;}
    [[nodiscard]] const Field& velocity() const noexcept{return cover_v_;}
    [[nodiscard]] const std::array<double,nz>& area() const noexcept{return area_;}
    [[nodiscard]] double flow() const noexcept{return flow_;}
    [[nodiscard]] double flow_resistance() const noexcept{return flow_resistance_;}
    [[nodiscard]] double flow_inertance() const noexcept{return flow_inertance_;}
    [[nodiscard]] double jet_area() const noexcept{return jet_area_;}
    [[nodiscard]] double mean_gap() const noexcept;
    [[nodiscard]] double minimum_area() const noexcept{return minimum_area_;}
    [[nodiscard]] double kinetic_energy() const noexcept;
    [[nodiscard]] double elastic_energy() const noexcept;
private:
    void geometry();
    void aerodynamics(double lung,double outlet);
    TissueMaterial material_{};
    TissueRelaxation relaxation_{};
    Field cover_{},body_{},cover_v_{},gap_{},rest_{},pressure_{},contact_{};
    std::array<double,modes> displacement_{},velocity_{},equilibrium_{},force_{};
    std::array<double,nz> area_{},resistance_{};
    double dx_{},dz_{},flow_{},minimum_area_{},reference_gap_{};
    double holding_force_{},holding_work_{},tension_scale_{1};
    std::array<double,modes*modes> geometric_stiffness_{};
    double hardening_energy_{},flow_resistance_{},flow_inertance_{},jet_area_{};
};
}
