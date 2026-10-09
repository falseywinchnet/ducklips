#pragma once
#include "ducklips/articulation.hpp"
#include "ducklips/physics.hpp"
#include <atomic>
#include <cstddef>
#include <vector>

namespace ducklips {

// One behavioural instruction to the body. Cues set targets; the articulator
// dynamics, the folds and the airflow decide what actually happens.
enum class CueKind {
    vowel,        // tongue/jaw/lip posture; angular_rate <= 0 places it at once
    constriction, // set a slot's gesture and drive its activation toward value
    release,      // drive a slot's activation toward value (normally 0)
    velum,        // velar port area target (m^2, reference scale)
    glottis,      // glottal rest area target (m^2)
    lungs,        // subglottal pressure target (Pa)
    pitch,        // intended F0 (Hz) for the laryngeal pitch controller
    effort,       // input energy: lung pressure multiplier, with firmer adduction above 1
    loudness,     // output modulation: gain applied to the microphone signal
    finish        // the utterance is over once this time is reached
};

struct Cue final {
    SampleCount time{};
    CueKind kind{CueKind::finish};
    ConstrictionSlot slot{ConstrictionSlot::lips};
    VowelShape vowel{};
    ConstrictionGesture gesture{};
    double value{};
    double angular_rate{};
};

// Orders cues by time, keeping the planner's order for equal times.
struct CueEarlier final {
    bool operator()(const Cue& left, const Cue& right) const noexcept;
};

// Single-producer single-consumer ring of cues. The control thread pushes; the
// audio thread peeks and pops. Capacity is fixed at prepare(); push never
// allocates and fails when full.
class CueQueue final {
public:
    void prepare(std::size_t capacity);
    [[nodiscard]] bool push(const Cue& cue) noexcept;
    [[nodiscard]] bool peek(Cue& output) const noexcept;
    void pop() noexcept;
    [[nodiscard]] std::size_t free_space() const noexcept;

private:
    std::vector<Cue> slots_{};
    std::atomic<std::size_t> write_{0};
    std::atomic<std::size_t> read_{0};
};

}
