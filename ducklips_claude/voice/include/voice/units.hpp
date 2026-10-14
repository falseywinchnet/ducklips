#pragma once

// Physical constants and sampling of the voice body. Units are CGS
// (centimetres, grams, seconds, dyn/cm^2) inside the body; the score
// language uses pascals and hertz for readability.

#include <cstdint>

namespace voice {

using Index = std::int32_t;

constexpr double kSoundSpeed = 35000.0;     // cm/s, warm humid air
constexpr double kAirDensity = 0.00114;     // g/cm^3
constexpr double kAirViscosity = 0.15;      // cm^2/s, kinematic
constexpr double kPi = 3.14159265358979323846;

constexpr double kAudioRate = 48000.0;      // output samples per second
constexpr Index kStepsPerSample = 2;        // wave steps per output sample
constexpr double kStepRate = kAudioRate * kStepsPerSample;
constexpr double kStepSeconds = 1.0 / kStepRate;
constexpr double kSectionLength = kSoundSpeed / kStepRate;  // 0.3646 cm

constexpr double kPascal = 10.0;            // dyn/cm^2 per pascal

// Characteristic impedance factor: Z = rho c / A.
constexpr double kRhoC = kAirDensity * kSoundSpeed;

}  // namespace voice
