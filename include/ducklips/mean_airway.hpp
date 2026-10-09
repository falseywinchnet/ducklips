#pragma once
#include "ducklips/acoustic_network.hpp"
#include "ducklips/cavity_network.hpp"

namespace ducklips {
// Spatial quasi-steady pressure storage on the acoustic airway topology.
// Viscous resistance is distributed by length. Bernoulli acceleration losses
// telescope through contractions; expansions assume zero pressure recovery.
// This deliberately simple closure is not calibrated human aerodynamics.
class MeanAirway final {
public:
    void prepare(const AcousticNetwork& layout,const TractShape& shape,
                 const AirProperties& air,double oral_wall_compliance);
    void reset() noexcept;
    void target(const TractShape& shape) noexcept;
    // fraction is the end of this sample in the control block. Inlet pressure
    // and resistance may include the acoustic load for implicit glottal flow.
    [[nodiscard]] CavityStepReport advance(double fraction,double glottal_area,
        double inlet_pressure,double inlet_resistance,double dt) noexcept;
    [[nodiscard]] double pressure() const noexcept { return network_.pressures()[0]; }
    [[nodiscard]] double glottal_flow() const noexcept { return -network_.flows()[0]; }
    [[nodiscard]] double glottal_load_pressure() const noexcept { return pressure()+inlet_viscosity_*glottal_flow(); }
    [[nodiscard]] double oral_flow(std::size_t section) const noexcept;
    [[nodiscard]] double oral_area(std::size_t section) const noexcept { return area_[section]; }
    [[nodiscard]] double lip_flow() const noexcept { return network_.flows()[lip_port_]; }
    [[nodiscard]] double nasal_flow() const noexcept { return network_.flows()[nasal_port_]; }
    [[nodiscard]] std::span<const double> pressures() const noexcept { return network_.pressures(); }
private:
    struct Link final { std::size_t from{},to{};bool port{false};bool side_entrance{false}; };
    void copy_shape(const TractShape& shape,std::vector<double>& area,std::vector<double>& length) noexcept;
    void build(double fraction) noexcept;
    AirProperties air_{};
    CavityNetwork network_{};
    std::vector<Link> links_{};
    std::vector<CavityPassage> passages_{};
    std::vector<std::size_t> oral_ports_{};
    std::vector<double> old_area_{},new_area_{},old_length_{},new_length_{},area_{},length_{},
        volume_{},compliance_{},wall_compliance_{},inflow_{},initial_pressure_{};
    std::size_t oral_{},nasal_{},lip_port_{},nasal_port_{};
    double old_port_area_{},new_port_area_{},old_port_length_{},new_port_length_{};
    double inlet_viscosity_{};
    bool ready_{false};
};
}
