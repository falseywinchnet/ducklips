#pragma once
#include "ducklips/synthesizer.hpp"
#include <string>
#include <vector>

namespace ducklips {
// Control-thread adapter. Nudges change targets on a running body; no text,
// recording, waveform splice, or engine reset is involved. Units: seconds,
// Hz, Pa, m^2, rad/s. Times are absolute since prepare, scaled by speed.
// submit() and render() must be externally serialized on one owner thread.
// For separate producer/audio threads use Synthesizer::enqueue directly.
// Pending gestures may overlap. Cues are sorted before being delivered at
// each 1 ms boundary. Queue/memory is bounded by pending_capacity.
class NudgeStream final {
public:
    void prepare(const VoiceProfile& body, double speed = 1.0, std::size_t pending_capacity = 4096);
    // "TIME vowel PHONE RATE", "TIME diphthong PHONE DURATION RATE",
    // "TIME constrict PHONE ACTIVATION RATE", "TIME release SLOT RATE",
    // "TIME pitch|lungs|glottis|velum|effort|loudness VALUE RATE".
    // SLOT = lips|tip|body. Empty/comment lines are accepted.
    // Errors leave pending cues untouched. Past/control-committed cues fail.
    [[nodiscard]] bool submit(const std::string& line, std::string& error);
    void render(std::span<float> stereo) noexcept;
    [[nodiscard]] const Synthesizer& engine() const noexcept;
    [[nodiscard]] std::size_t pending() const noexcept;
private:
    Synthesizer engine_{};
    std::vector<Cue> pending_{};
    std::size_t capacity_{};
    double speed_{1.0};
    SampleCount committed_until_{};
};
}
