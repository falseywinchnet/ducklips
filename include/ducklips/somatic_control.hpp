#pragma once
#include "ducklips/articulation.hpp"
#include "ducklips/motor_hold.hpp"
#include <array>
#include <cstdint>

namespace ducklips {
// Shared maintenance mechanism for the existing articulated coordinates.
// These are reduced coordinate dynamics, not measured anatomical muscles.
// Targets remain intents. Position, velocity, stress memory, twitch force,
// delayed sensing and adaptation persist across cue changes.
class SomaticControl final {
public:
    static constexpr std::size_t count=29;
    struct Settings final {
        double noise{.08}; // fraction of finite motor-pool fluctuation retained
        double tone{.02}; // opposing activation, including at neutral posture
        double feedback{1};
        double relaxation{.15};
        double respiratory_compliance{1e-6}; // m3/Pa, prototype body
        double breath_volume{.003}; // m3 available during this exhalation
    };
    static void validate(const Settings& settings);
    struct Observation final {
        double target{},position{},velocity{},stress_memory{},muscle{},sensed{},correction{},load{};
    };
    void reset(ArticulatorState& body,std::uint64_t seed,const Settings& settings);
    void advance(ArticulatorState& body,double seconds,double airway_pressure,double expired_flow);
    // An idle gesture acquires a contact reference, not a physical jump of
    // the tongue. Carry its small displacement and neural history with it.
    void relocate_idle_gesture(ArticulatorState& body,std::size_t slot) noexcept;
    // Smooth physiological coordinates across an acoustic control block.
    void physiology(ArticulatorState& body,double fraction) const noexcept;
    [[nodiscard]] const std::array<Observation,count>& observations() const noexcept{return observations_;}
    [[nodiscard]] static const char* name(std::size_t index) noexcept;
    [[nodiscard]] double remaining_breath() const noexcept{return breath_;}
    // Diagnostic generalized load in normalized coordinate units. The same
    // input is used for the tract-pressure and contact-load adapters.
    void set_load(std::size_t channel,double load) noexcept{loads_[channel]=load;}
private:
    struct Coordinate final {
        double origin{},scale{1},low{},high{},q{},v{},memory{},sense{},integral{};
        double twitch_input{},twitch{},mean_input{},mean_twitch{},previous{},current{};
        double rate_noise{},rate{40};
        std::array<double,30> delayed{};
        MotorPool positive{},negative{};
        std::uint64_t random{1};
        std::size_t group{},cursor{};
    };
    static std::array<CriticallyDampedTracker*,count> coordinates(ArticulatorState& body) noexcept;
    static double noise(std::uint64_t& state) noexcept;
    static double drive(double activation) noexcept;
    void tick(ArticulatorState& body,double pressure,double flow);
    Settings settings_{};
    std::array<Coordinate,count> state_{};
    std::array<Observation,count> observations_{};
    std::array<double,count> loads_{};
    std::array<double,5> common_{};
    std::uint64_t random_{1};
    double breath_{};
};
}
