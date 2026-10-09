#include "ducklips/physics.hpp"
#include <cmath>

namespace ducklips {

SampleCount samples_from_seconds(SynthesisRate rate, double seconds) noexcept {
    double exact = seconds * rate.samples_per_second;
    SampleCount rounded = static_cast<SampleCount>(std::llround(exact));
    return rounded;
}

double seconds_from_samples(SynthesisRate rate, SampleCount samples) noexcept {
    double seconds = static_cast<double>(samples) / rate.samples_per_second;
    return seconds;
}

}
