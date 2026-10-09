#pragma once
#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace ducklips {
// Experimental frozen linear acoustics only. Geometry and mean-flow losses
// cannot be changed by this representation. Preparation may allocate; ticks do not.
class FrozenModes final {
public:
    void load(const std::string& path);
    void reset() noexcept;
    void tick(std::span<const double> inputs, std::span<double> outputs);
    [[nodiscard]] double free_pressure(std::size_t output) const noexcept;
    [[nodiscard]] double direct(std::size_t output, std::size_t input) const noexcept;
    // Common voiced case: only the glottal flow drives the acoustic network.
    // Return the radiation flow and advance states; pressure was read above.
    double tick_glottis(double flow) noexcept;
    [[nodiscard]] double sample_rate() const noexcept { return rate_; }
    [[nodiscard]] std::size_t inputs() const noexcept { return inputs_; }
    [[nodiscard]] std::size_t outputs() const noexcept { return outputs_; }
    [[nodiscard]] std::size_t modes() const noexcept { return pole_r_.size(); }
private:
    double rate_{};
    std::size_t inputs_{}, outputs_{}, complex_count_{};
    std::vector<double> pole_r_, pole_i_, state_, drive_, observe_, direct_;
};
}
