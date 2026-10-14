#pragma once
#include <cstdint>

namespace ducklips {

inline constexpr double pi = 3.14159265358979323846;

// SI units throughout DuckLips: metres, square metres, kilograms, seconds,
// pascals, cubic metres per second. Literature values quoted in cgs units are
// converted where they are declared.
inline constexpr double centimetre = 1.0e-2;        // m
inline constexpr double square_centimetre = 1.0e-4; // m^2

// Warm, saturated air inside the vocal tract.
struct AirProperties final {
    double density{1.14};               // kg/m^3
    double sound_speed{350.0};          // m/s
    double viscosity{1.86e-5};          // Pa s, dynamic
    double heat_capacity_ratio{1.4};
    double thermal_diffusivity{2.3e-5}; // m^2/s, conductivity / (density c_p)
};

// Audio time is counted in samples; planning time is in seconds. Conversion
// from seconds rounds half away from zero to the nearest sample.
using SampleCount = std::int64_t;

struct SynthesisRate final {
    double samples_per_second{48000.0};
};

[[nodiscard]] SampleCount samples_from_seconds(SynthesisRate rate, double seconds) noexcept;
[[nodiscard]] double seconds_from_samples(SynthesisRate rate, SampleCount samples) noexcept;

}
