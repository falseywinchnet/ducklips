#pragma once

// Units are CGS throughout (cm, g, s, dyn, dyn/cm^2 = 0.1 Pa), as in the
// voice/ engine and in VocalTractLab.

#include <cstddef>
#include <cstdint>

namespace tract {

using Index = std::int32_t;

constexpr double kPi = 3.14159265358979323846;
constexpr double kHalfPi = 0.5 * kPi;
constexpr double kTwoPi = 2.0 * kPi;

constexpr double kStepRate = 96000.0;              // state updates per second
constexpr double kStepSeconds = 1.0 / kStepRate;
constexpr double kAudioRate = 48000.0;             // written audio
// The control tick: muscles, tract shape, breath, pitch command, variation
// and jet placement update every kControlSteps state steps; the steps
// between interpolate. 12 steps = 0.125 ms (2026-10-15, the owner: timing
// resolution 8x; it was 96 steps, 1 ms).
constexpr long kControlSteps = 12;
constexpr double kControlSeconds = static_cast<double>(kControlSteps) * kStepSeconds;
constexpr double kAirDensity = 1.14e-3;            // g/cm^3, warm moist air
constexpr double kSoundSpeed = 35000.0;            // cm/s
constexpr double kRhoC = kAirDensity * kSoundSpeed;
constexpr double kAirViscosity = 1.86e-4;          // poise (dynamic)
// Air moves one node per step: nodes are c / step rate apart (0.3646 cm), so
// every axial air path is an exact one-step delay.
constexpr double kNodeLength = kSoundSpeed / kStepRate;
constexpr double kPascal = 10.0;                   // dyn/cm^2 per Pa
constexpr double kMinimumArea = 1.0e-4;            // cm^2, an air path never quite vanishes

}  // namespace tract
