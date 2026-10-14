#pragma once

// Holding under effort, ported from Astra's somatic control (web/tube
// commit 42c7337: motor_hold.hpp, somatic_control.cpp).
//
// A held coordinate is not a spring pulled to its target. Two opposing motor
// pools (64 units each) fire at irregular intervals; each discharge starts an
// overlapping twitch (15 ms, or 35 ms for every tenth unit). The coordinate
// is sensed through a 20 ms filter and a 30 ms transport delay, and the error
// changes the requested effort (feedforward, proportional, integral). The
// tissue keeps a strain memory that relaxes toward its current deformation at
// a fluctuating positive rate, so a held posture creeps and must be re-held.
// Nothing here is periodic: no oscillator, no prescribed vibration.
//
// HYPOTHESIS (Astra's, kept): 15/35 ms twitches and interval CV 0.2 follow
// Manriquez et al. (2017); unit count, thresholds, maximum force, the 0.08
// retained pool fluctuation, 0.02 opposing tone and 0.15 relaxation are model
// assumptions, not measured physiology.

#include <array>
#include <cstdint>

namespace tract {

// A finite pool of motor units: irregular discharges, overlapping twitches.
class MotorPool {
 public:
  static constexpr int kUnits = 64;

  void reset(std::uint64_t seed);
  // Mean activation the pool produces at a drive (0..1), without the
  // discharge irregularity.
  double expected_activation(double drive) const;
  double drive_for_activation(double activation) const;
  // Advance by dt seconds at a drive; returns the activation (0..~1).
  double advance(double drive, double dt);
  double activation() const { return fast_y_ + slow_y_; }

 private:
  static double firing_rate(double drive, int unit);
  double next_interval();
  double gauss();

  std::array<double, kUnits> phase_{};
  std::array<double, kUnits> interval_{};
  double fast_x_ = 0.0;
  double fast_y_ = 0.0;
  double slow_x_ = 0.0;
  double slow_y_ = 0.0;
  std::uint64_t random_ = 1U;
  double spare_ = 0.0;
  bool has_spare_ = false;
};

struct HoldSettings {
  double noise = 0.08;       // share of the finite pool's fluctuation retained
  double tone = 0.02;        // opposing activation, also at neutral posture
  double feedback = 1.0;
  double relaxation = 0.15;  // share of the posture carried by relaxing strain memory
  // Common drive: the slow (tau 0.15 s) synaptic fluctuation a motor pool's
  // units share, as a share of the effort. Astra's model kept 0.002 x noise.
  // HYPOTHESIS: force-steadiness studies put the coefficient of variation of
  // a held moderate contraction near 1-3%, mostly below ~2 Hz and mostly
  // from common input; 0.02 is in that range, not fitted.
  double common = 0.02;
  // The fastest a coordinate is moved, 1/s (Astra's 180; speech gestures
  // such as a 28 ms tip closure need more).
  double fastest = 180.0;
  // The burst: the request passes through ((1 + T s) / (1 + b s))^2, T the
  // twitch time (17 ms), so the force follows the intent with b instead of
  // T (an agonist burst, the inverse of the muscle's own lag, clipped by the
  // pools' capacity). 0: no burst (Astra's: half a move took 41 ms).
  double burst = 0.006;
};

// The coordinates held together, ticked at 1 kHz. Each coordinate belongs to
// a group whose members share a slow neural fluctuation (tau 0.15 s).
class SomaticHold {
 public:
  static constexpr int kMaxCoordinates = 48;
  static constexpr int kGroups = 6;
  static constexpr double kTick = 0.001;  // s

  void reset(std::uint64_t seed, const HoldSettings& settings);
  // A coordinate in physical units: its rest value (origin), the span that
  // is one unit of effort (scale), its limits and group. Returns its index.
  int add(double origin, double scale, double low, double high, int group);
  // The target and the speed of reaching it (1/s, clamped 8..fastest).
  void set_target(int index, double value, double rate);
  // An external load on a coordinate, in units of its scale.
  void set_load(int index, double load) { coordinates_[static_cast<std::size_t>(index)].load = load; }
  // Place every coordinate at its target, held and settled (the opening
  // posture of an utterance).
  void settle();
  // One 1 ms tick of every coordinate.
  void tick();
  double value(int index) const { return coordinates_[static_cast<std::size_t>(index)].value; }
  double effort(int index) const { return coordinates_[static_cast<std::size_t>(index)].muscle; }
  int count() const { return count_; }

 private:
  struct Coordinate {
    double origin = 0.0;
    double scale = 1.0;
    double low = 0.0;
    double high = 0.0;
    int group = 0;
    double target = 0.0;  // normalized
    double rate = 40.0;
    double load = 0.0;
    double q = 0.0;
    double v = 0.0;
    double memory = 0.0;
    double expected = 0.0;  // the internal model's memory, at the nominal rate
    double sense = 0.0;
    std::array<double, 30> delayed{};
    // Efference copy: the trajectory the command predicts, sensed through
    // the same filter and delay.
    double predicted_q = 0.0;
    double predicted_v = 0.0;
    double predicted_sense = 0.0;
    std::array<double, 30> predicted_delayed{};
    double predicted_force = 0.0;   // the predicted equilibrium, first stage
    double predicted_settle = 0.0;  // ... second stage
    double burst_first = 0.0;       // the burst's two low-pass states
    double burst_second = 0.0;
    int cursor = 0;
    double integral = 0.0;
    double rate_noise = 0.0;
    double twitch = 0.0;
    double twitch_input = 0.0;
    double mean_twitch = 0.0;
    double mean_input = 0.0;
    double muscle = 0.0;
    double value = 0.0;
    std::uint64_t random = 1U;
    MotorPool positive;
    MotorPool negative;
  };

  void build_recruitment();
  double drive(double activation) const;
  void tick_coordinate(Coordinate& c);
  static double noise(std::uint64_t& state);

  std::array<Coordinate, kMaxCoordinates> coordinates_{};
  int count_ = 0;
  std::array<double, kGroups> common_{};
  std::array<double, 257> recruitment_{};  // drive for activation i/256
  HoldSettings settings_;
  std::uint64_t seed_ = 1U;
  std::uint64_t random_ = 1U;
};

}  // namespace tract
