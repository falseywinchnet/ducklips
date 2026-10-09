#pragma once
#include "ducklips/moving_tree.hpp"
namespace ducklips {
// End-to-end integration reference. Full physical branch history, continuously
// interpolated C/G, constant series sources per audio sample, implicit glottal
// feedback via the separately integrated unit volume-flow response.
// Losses: positive diffusive quadrature of sqrt(2s) for viscous/thermal
// channels; physical wall resistance; passive one-pole radiation load.
// Numerical reference, not a calibrated human body or a fast renderer.
class CanonicalSpeech final {
public:
    void prepare(const AcousticNetwork& layout,const AirProperties& air);
    void reset() noexcept;
    void free_step(const AcousticNetwork& start,const AcousticNetwork& end,
                   double fraction,double increment,double dt,std::size_t substeps,
                   std::size_t aspiration_branch,std::size_t frication_branch,
                   double aspiration,double frication) noexcept;
    void inject(double glottal_fluctuation) noexcept;
    // Work-conjugate average pressure and response to unit held nodal flow.
    [[nodiscard]] double load_pressure() const noexcept;
    [[nodiscard]] double load_impedance() const noexcept;
    [[nodiscard]] double pressure() const noexcept;
    [[nodiscard]] double energy() const noexcept;
    [[nodiscard]] double flow_derivative(std::size_t branch,double dt) const noexcept;
private:
    static constexpr std::size_t loss_terms=24;
    void losses(MovingTreeReference& state,std::vector<double>& memory,double h) noexcept;
    void close_passages(MovingTreeReference& state,std::vector<double>& memory) noexcept;
    double inertive_conductance(const AcousticNetwork& network,std::size_t branch) const noexcept;
    std::array<double,loss_terms> loss_poles_{},loss_weights_{};
    std::vector<double> physical_memory_{},response_memory_{},loss_scale_{},static_resistance_{},rad_pole_{},rad_weight_{};
    double density_{};
    double load_pressure_{},load_impedance_{};
    MovingTreeReference physical_{},response_{};
    std::vector<double> c1_{},c2_{},g1_{},g2_{},cend_{},gend_{},old_flow_{},nodal_{},series_{};
};
}
