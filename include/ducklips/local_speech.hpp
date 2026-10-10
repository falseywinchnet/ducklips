#pragma once
#include "ducklips/acoustic_network.hpp"
#include <vector>

namespace ducklips {
// Local air storage and momentum exchange. Fixed forward/backward sweeps;
// no global matrix, eigenbasis, nonlinear iteration or rejected timestep.
// Both mean and acoustic flow use the same physical state.
class LocalSpeech final {
public:
    void prepare(const AcousticNetwork& layout, const AirProperties& air, double wall_compliance,
                 double glottal_length=0.003);
    void reset() noexcept;
    void advance(const AcousticNetwork& first, const AcousticNetwork& last,
                 double fraction, double increment, double dt, std::size_t substeps,
                 double glottal_area, double lung_pressure, std::size_t aspiration_branch,
                 std::size_t frication_branch, std::size_t jet_branch, double aspiration, double frication) noexcept;
    [[nodiscard]] double pressure(std::size_t node=0) const noexcept;
    [[nodiscard]] double flow(std::size_t branch) const noexcept;
    [[nodiscard]] double flow_derivative(std::size_t branch,double dt) const noexcept;
    [[nodiscard]] double glottal_flow() const noexcept { return glottal_flow_; }
    [[nodiscard]] double glottal_endpoint_flow() const noexcept {
        return glottal_conductance_*glottal_momentum_;
    }
    [[nodiscard]] double glottal_kinetic_energy() const noexcept {
        return .5*glottal_conductance_*glottal_momentum_*glottal_momentum_;
    }
    [[nodiscard]] double energy() const noexcept;
    // Exposed physical state permits independent conservation/reference tests.
    std::vector<double> charge{}, momentum{}, wall_charge{}, radiation_momentum{};
private:
    void exchange(std::size_t branch, double series) noexcept;
    void inject(double area, double lung, double h) noexcept;
    AirProperties air_{};
    std::vector<AcousticBranch> branches_{};
    std::vector<double> inverse_c_{},g_{},wall_inverse_c_{},old_flow_{},volume_{};
    std::vector<double> transfer_p_{},transfer_d_{},momentum_p_{},momentum_d_{};
    std::vector<double> radiation_h_{},radiation_decay_{};
    std::vector<double> jet_coefficient_{},half_step_{};
    double wall_compliance_{},glottal_flow_{},glottal_length_{},glottal_momentum_{},glottal_conductance_{};
    bool ready_{};
};
}
