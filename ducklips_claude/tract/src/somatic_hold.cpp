#include "tract/somatic_hold.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

constexpr double kFastTwitch = 0.015;  // s
constexpr double kSlowTwitch = 0.035;  // s, every tenth unit
constexpr double kCapacity = 8.0;      // requested effort at full recruitment, normalized units
constexpr double kTwitch = 0.017;      // s, the mean-field twitch time

std::uint64_t next_random(std::uint64_t& state) {
  state ^= state >> 12U;
  state ^= state << 25U;
  state ^= state >> 27U;
  return state * 2685821657736338717ULL;
}

double uniform(std::uint64_t& state) {
  return static_cast<double>(next_random(state) >> 11U) * 0x1.0p-53;
}

}  // namespace

// --- MotorPool ---------------------------------------------------------------

void MotorPool::reset(std::uint64_t seed) {
  random_ = seed != 0U ? seed : 1U;
  has_spare_ = false;
  for (int i = 0; i < kUnits; ++i) {
    phase_[static_cast<std::size_t>(i)] = uniform(random_);
    interval_[static_cast<std::size_t>(i)] = next_interval();
  }
  fast_x_ = 0.0;
  fast_y_ = 0.0;
  slow_x_ = 0.0;
  slow_y_ = 0.0;
}

double MotorPool::gauss() {
  if (has_spare_) {
    has_spare_ = false;
    return spare_;
  }
  double u = std::max(uniform(random_), 1e-300);
  double w = uniform(random_);
  double radius = std::sqrt(-2.0 * std::log(u));
  spare_ = radius * std::sin(2.0 * 3.14159265358979323846 * w);
  has_spare_ = true;
  return radius * std::cos(2.0 * 3.14159265358979323846 * w);
}

double MotorPool::next_interval() {
  return std::max(0.2, 1.0 + 0.2 * gauss());
}

double MotorPool::firing_rate(double drive, int unit) {
  double threshold = 0.5 * static_cast<double>(unit) / kUnits;
  if (drive <= threshold) {
    return 0.0;
  }
  return 8.0 + 32.0 * std::clamp((drive - threshold) / (1.0 - threshold), 0.0, 1.0);
}

double MotorPool::expected_activation(double drive) const {
  double sum = 0.0;
  for (int i = 0; i < kUnits; ++i) {
    sum += firing_rate(drive, i) / 40.0;
  }
  return sum / kUnits;
}

double MotorPool::drive_for_activation(double activation) const {
  if (activation <= 0.0) {
    return 0.0;
  }
  if (activation >= 1.0) {
    return 1.0;
  }
  double low = 0.0;
  double high = 1.0;
  for (int i = 0; i < 12; ++i) {
    double middle = 0.5 * (low + high);
    if (expected_activation(middle) < activation) {
      low = middle;
    } else {
      high = middle;
    }
  }
  return 0.5 * (low + high);
}

double MotorPool::advance(double drive, double dt) {
  double fast_decay = std::exp(-dt / kFastTwitch);
  double slow_decay = std::exp(-dt / kSlowTwitch);
  fast_y_ = (fast_y_ + dt / kFastTwitch * fast_x_) * fast_decay;
  fast_x_ *= fast_decay;
  slow_y_ = (slow_y_ + dt / kSlowTwitch * slow_x_) * slow_decay;
  slow_x_ *= slow_decay;
  for (int i = 0; i < kUnits; ++i) {
    double rate = firing_rate(drive, i);
    if (rate == 0.0) {
      continue;
    }
    std::size_t u = static_cast<std::size_t>(i);
    double old = phase_[u];
    phase_[u] += rate * dt;
    if (phase_[u] >= interval_[u]) {
      // The discharge fell inside this step: its twitch starts there.
      double remaining = dt - (interval_[u] - old) / rate;
      bool slow = i % 10 == 0;
      double tau = slow ? kSlowTwitch : kFastTwitch;
      // Per-unit force: full drive gives a mean activation of 1.
      double jump = 1.0 / (kUnits * 40.0 * tau);
      double tail = jump * std::exp(-remaining / tau);
      if (slow) {
        slow_x_ += tail;
        slow_y_ += remaining / tau * tail;
      } else {
        fast_x_ += tail;
        fast_y_ += remaining / tau * tail;
      }
      phase_[u] -= interval_[u];
      interval_[u] = next_interval();
    }
  }
  return activation();
}

// --- SomaticHold -------------------------------------------------------------

double SomaticHold::noise(std::uint64_t& state) {
  // Unit-variance uniform noise (Astra's): sqrt(12) (u - 1/2).
  return std::sqrt(12.0) * (uniform(state) - 0.5);
}

void SomaticHold::build_recruitment() {
  MotorPool pool;
  for (std::size_t i = 0; i < recruitment_.size(); ++i) {
    recruitment_[i] = pool.drive_for_activation(static_cast<double>(i) / 256.0);
  }
}

double SomaticHold::drive(double activation) const {
  double x = std::clamp(activation, 0.0, 1.0) * 256.0;
  std::size_t i = std::min(static_cast<std::size_t>(x), static_cast<std::size_t>(255));
  double f = x - static_cast<double>(i);
  return recruitment_[i] + f * (recruitment_[i + 1] - recruitment_[i]);
}

void SomaticHold::reset(std::uint64_t seed, const HoldSettings& settings) {
  settings_ = settings;
  seed_ = seed != 0U ? seed : 1U;
  random_ = seed_;
  common_.fill(0.0);
  count_ = 0;
  build_recruitment();
}

int SomaticHold::add(double origin, double scale, double low, double high, int group) {
  if (count_ >= kMaxCoordinates) {
    return -1;
  }
  int index = count_++;
  Coordinate& c = coordinates_[static_cast<std::size_t>(index)];
  c = Coordinate{};
  c.origin = origin;
  c.scale = scale;
  c.low = low;
  c.high = high;
  c.group = std::clamp(group, 0, kGroups - 1);
  c.value = origin;
  c.random = seed_ + 1009U * static_cast<std::uint64_t>(index + 1);
  if (c.random == 0U) {
    c.random = 1U;
  }
  c.positive.reset(c.random);
  c.negative.reset(c.random + 503U);
  return index;
}

void SomaticHold::set_target(int index, double value, double rate) {
  Coordinate& c = coordinates_[static_cast<std::size_t>(index)];
  c.target = (std::clamp(value, c.low, c.high) - c.origin) / c.scale;
  c.rate = std::clamp(rate, 8.0, settings_.fastest);
}

void SomaticHold::settle() {
  for (int i = 0; i < count_; ++i) {
    Coordinate& c = coordinates_[static_cast<std::size_t>(i)];
    c.q = c.target;
    c.v = 0.0;
    c.memory = c.target;
    c.expected = c.target;
    c.sense = c.target;
    c.delayed.fill(c.target);
    c.predicted_q = c.target;
    c.predicted_v = 0.0;
    c.predicted_sense = c.target;
    c.predicted_delayed.fill(c.target);
    c.integral = 0.0;
    double request = (1.0 - settings_.relaxation) * c.target;
    c.predicted_force = c.target;
    c.predicted_settle = c.target;
    c.burst_first = request;
    c.burst_second = request;
    c.twitch = request;
    c.twitch_input = request;
    c.muscle = request;
    c.value = c.origin + c.scale * c.q;
  }
}

void SomaticHold::tick() {
  double common_decay = std::exp(-kTick / 0.15);
  double common_gain = std::sqrt(1.0 - common_decay * common_decay);
  for (double& value : common_) {
    value = common_decay * value + common_gain * noise(random_);
  }
  for (int i = 0; i < count_; ++i) {
    tick_coordinate(coordinates_[static_cast<std::size_t>(i)]);
  }
}

void SomaticHold::tick_coordinate(Coordinate& c) {
  const double dt = kTick;
  double common = common_[static_cast<std::size_t>(c.group)];
  // The strain memory relaxes toward the present deformation at a positive
  // rate that wanders (log-normal, 0.12 s local and 0.15 s shared).
  double local_decay = std::exp(-dt / 0.12);
  c.rate_noise = local_decay * c.rate_noise + std::sqrt(1.0 - local_decay * local_decay) * noise(c.random);
  double relaxation_rate = std::exp(0.12 * settings_.noise * (c.rate_noise + common)) / 0.2;
  c.memory += -std::expm1(-dt * relaxation_rate) * (c.q - c.memory);
  // The internal model expects the memory to relax at its nominal rate; the
  // muscle supplies what that memory lacks. (Astra's request, (1 - g) of the
  // target, assumed the memory already at the target: during a move it lags
  // 0.2 s, the integral wound up, and every posture overshot ~7% of its
  // travel and returned over seconds.) The true rate wanders; the
  // difference is the creep the feedback must catch.
  c.expected += -std::expm1(-dt / 0.2) * (c.q - c.expected);
  // Sensing: 20 ms filter, 30 ms transport delay.
  c.sense += -std::expm1(-dt / 0.02) * (c.q - c.sense);
  std::size_t cursor = static_cast<std::size_t>(c.cursor);
  double delayed = c.delayed[cursor];
  c.delayed[cursor] = c.sense;
  c.cursor = (c.cursor + 1) % static_cast<int>(c.delayed.size());
  // The command's own prediction, sensed alike: feedback acts only on what
  // the prediction did not foresee (creep, twitches, loads). Error taken
  // from the target wound the integral up during every move.
  // The predicted equilibrium follows the target through the force's lag
  // (two first-order stages: the burst's time, or the twitch's without it).
  double lag = settings_.burst > 0.0 ? settings_.burst : 0.017;
  double lag_gain = -std::expm1(-dt / lag);
  c.predicted_force += lag_gain * (c.target - c.predicted_force);
  c.predicted_settle += lag_gain * (c.predicted_force - c.predicted_settle);
  double predicted_displacement = c.predicted_q - c.predicted_settle;
  double predicted_slope = c.predicted_v + c.rate * predicted_displacement;
  double predicted_decay = std::exp(-c.rate * dt);
  c.predicted_q = c.predicted_settle + predicted_decay * (predicted_displacement + dt * predicted_slope);
  c.predicted_v = predicted_decay * (c.predicted_v - c.rate * dt * predicted_slope);
  c.predicted_sense += -std::expm1(-dt / 0.02) * (c.predicted_q - c.predicted_sense);
  double predicted_delayed = c.predicted_delayed[cursor];
  c.predicted_delayed[cursor] = c.predicted_sense;
  double error = predicted_delayed - delayed;
  double integral = c.integral + settings_.feedback * 0.15 * error * dt / 0.2;
  double intent = c.target - settings_.relaxation * c.expected + settings_.feedback * 0.15 * error + integral;
  if (std::abs(intent) < kCapacity || intent * error < 0.0) {
    c.integral = std::clamp(integral, -2.0, 2.0);
  }
  double requested = intent;
  if (settings_.burst > 0.0) {
    // (1 + T s)/(1 + b s) = T/b + (1 - T/b)/(1 + b s), twice.
    double ratio = kTwitch / settings_.burst;
    double gain = -std::expm1(-dt / settings_.burst);
    c.burst_first += gain * (intent - c.burst_first);
    double stage = ratio * intent + (1.0 - ratio) * c.burst_first;
    c.burst_second += gain * (stage - c.burst_second);
    requested = ratio * stage + (1.0 - ratio) * c.burst_second;
  }
  double limit = kCapacity * (1.0 - settings_.tone);
  requested = std::clamp(requested, -limit, limit);
  requested *= 1.0 + settings_.common * common;
  double activation = std::abs(requested) / kCapacity;
  // Opposing pools: their mean tone cancels, their twitches do not.
  double positive_drive = drive(settings_.tone + (requested > 0.0 ? activation : 0.0));
  double negative_drive = drive(settings_.tone + (requested < 0.0 ? activation : 0.0));
  double pool = kCapacity * (c.positive.advance(positive_drive, dt) - c.negative.advance(negative_drive, dt));
  double mean = kCapacity * (c.positive.expected_activation(positive_drive) -
                             c.negative.expected_activation(negative_drive));
  // Mean-field twitch of the request carries the population; the finite
  // pool adds a share of its fluctuation about its own mean twitch.
  double decay = std::exp(-dt / kTwitch);
  c.twitch = requested + decay * ((c.twitch - requested) + dt / kTwitch * (c.twitch_input - requested));
  c.twitch_input = requested + (c.twitch_input - requested) * decay;
  c.mean_twitch = mean + decay * ((c.mean_twitch - mean) + dt / kTwitch * (c.mean_input - mean));
  c.mean_input = mean + (c.mean_input - mean) * decay;
  c.muscle = c.twitch + settings_.noise * (pool - c.mean_twitch);
  // Passive memory + active effort + load: the exact critically damped
  // update toward that equilibrium.
  double equilibrium = c.muscle + settings_.relaxation * c.memory + c.load;
  double displacement = c.q - equilibrium;
  double slope = c.v + c.rate * displacement;
  double body_decay = std::exp(-c.rate * dt);
  c.q = equilibrium + body_decay * (displacement + dt * slope);
  c.v = body_decay * (c.v - c.rate * dt * slope);
  double position = c.origin + c.scale * c.q;
  if (position < c.low || position > c.high) {
    position = std::clamp(position, c.low, c.high);
    c.q = (position - c.origin) / c.scale;
    c.v = 0.0;
  }
  c.value = position;
}

}  // namespace tract
