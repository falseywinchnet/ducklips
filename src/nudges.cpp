#include "ducklips/nudges.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace ducklips {
void NudgeStream::prepare(const VoiceProfile& body, double speed, std::size_t pending_capacity) {
    if (!std::isfinite(speed) || speed < 0.25 || speed > 4.0 || pending_capacity < 4 || pending_capacity > 65536) {
        throw std::invalid_argument("speed must be in [0.25,4]; cue capacity in [4,65536]");
    }
    SynthesisOptions options{};
    options.cue_capacity = pending_capacity;
    engine_.prepare(body, options);
    speed_ = speed;
    capacity_ = pending_capacity;
    pending_.clear();
    pending_.reserve(capacity_);
    committed_until_ = 0;
}

bool NudgeStream::submit(const std::string& line, std::string& error) {
    std::istringstream input(line.substr(0, line.find('#')));
    input >> std::ws;
    if (input.eof()) { error.clear(); return true; }
    double time = 0.0;
    std::string kind{};
    if (!(input >> time >> kind) || !std::isfinite(time) || time < 0.0 || time > 86400.0) {
        error = "expected finite time in [0,86400] followed by nudge kind"; return false;
    }
    Cue cue{};
    cue.time = samples_from_seconds(engine_.options().rate, time / speed_);
    if (cue.time < committed_until_) {
        error = "nudge is behind the committed control boundary"; return false;
    }
    std::vector<Cue> additions{};
    double rate = 0.0;
    if (kind == "vowel" || kind == "diphthong" || kind == "constrict") {
        std::string symbol{};
        Phone phone{};
        int stress = -1;
        if (!(input >> symbol) || !parse_phone(symbol, phone, stress)) {
            error = "unknown ARPAbet phone"; return false;
        }
        const PhoneSpec& spec = phone_spec(phone);
        if (kind == "constrict") {
            if (!(input >> cue.value >> rate) || !std::isfinite(cue.value) || cue.value < 0.0 || cue.value > 1.5
                || !spec.primary.present || is_vowel(phone)) {
                error = "constrict requires a consonant gesture, activation [0,1.5], and rate"; return false;
            }
            cue.kind = CueKind::constriction;
            cue.slot = spec.primary.slot;
            cue.gesture = gesture_from_spec(spec.primary);
            additions.push_back(cue);
            if (spec.secondary.present) {
                cue.slot = spec.secondary.slot;
                cue.gesture = gesture_from_spec(spec.secondary);
                additions.push_back(cue);
            }
        } else {
            if (!is_vowel(phone)) { error = "vowel target requires a vowel"; return false; }
            double duration = 0.0;
            if (kind == "diphthong" && (!(input >> duration) || !std::isfinite(duration)
                || duration < 0.02 || duration > 10.0 || spec.glide_to == Phone::count)) {
                error = "diphthong requires AY/AW/EY/OW/OY and duration [0.02,10]"; return false;
            }
            if (!(input >> rate)) { error = "missing motion rate"; return false; }
            cue.kind = CueKind::vowel;
            cue.vowel = engine_.vowels().shape(phone);
            additions.push_back(cue);
            if (kind == "diphthong") {
                // Two articulatory targets. The tissue dynamics produce the
                // continuous trajectory; no prerecorded transition is used.
                cue.time += samples_from_seconds(engine_.options().rate, 0.55 * duration / speed_);
                cue.vowel = engine_.vowels().shape(spec.glide_to);
                additions.push_back(cue);
            }
        }
    } else if (kind == "release") {
        std::string slot{};
        if (!(input >> slot >> rate)) { error = "release requires slot and rate"; return false; }
        if (slot == "lips") { cue.slot = ConstrictionSlot::lips; }
        else if (slot == "tip") { cue.slot = ConstrictionSlot::tongue_tip; }
        else if (slot == "body") { cue.slot = ConstrictionSlot::tongue_body; }
        else { error = "slot must be lips, tip, or body"; return false; }
        cue.kind = CueKind::release;
        additions.push_back(cue);
    } else {
        double low = 0.0;
        double high = 0.0;
        if (kind == "pitch") { cue.kind = CueKind::pitch; low = 60; high = 500; }
        else if (kind == "lungs") { cue.kind = CueKind::lungs; low = -500; high = 2000; }
        else if (kind == "glottis") { cue.kind = CueKind::glottis; low = -2e-6; high = 2e-4; }
        else if (kind == "velum") { cue.kind = CueKind::velum; high = 3e-4; }
        else if (kind == "effort") { cue.kind = CueKind::effort; low = 0.3; high = 2; }
        else if (kind == "loudness") { cue.kind = CueKind::loudness; high = 1; }
        else { error = "unknown nudge kind: " + kind; return false; }
        if (!(input >> cue.value >> rate) || !std::isfinite(cue.value) || cue.value < low || cue.value > high) {
            error = "invalid value or missing motion rate for " + kind; return false;
        }
        additions.push_back(cue);
    }
    std::string extra{};
    if ((input >> extra) || !std::isfinite(rate) || rate < 1.0 || rate > 500.0) {
        error = "motion rate must be finite in [1,500], with no trailing fields"; return false;
    }
    if (pending_.size() + additions.size() > capacity_) { error = "pending nudge capacity exceeded"; return false; }
    for (Cue& addition : additions) {
        // Speed scales requested motion, but body speed limits still apply.
        addition.angular_rate = rate * speed_;
        pending_.push_back(addition);
    }
    std::stable_sort(pending_.begin(), pending_.end(), CueEarlier{});
    error.clear();
    return true;
}

void NudgeStream::render(std::span<float> stereo) noexcept {
    const SampleCount block = static_cast<SampleCount>(engine_.options().control_block);
    std::size_t done = 0;
    const std::size_t frames = stereo.size() / 2;
    while (done < frames) {
        SampleCount position = engine_.position();
        if (position % block == 0) {
            committed_until_ = position + block;
            std::size_t delivered = 0;
            while (delivered < pending_.size() && pending_[delivered].time < committed_until_) {
                // At most capacity_ pending, and all previously queued cues
                // were consumed in the previous control block.
                if (!engine_.enqueue(pending_[delivered])) { break; }
                ++delivered;
            }
            pending_.erase(pending_.begin(), pending_.begin() + static_cast<std::ptrdiff_t>(delivered));
        }
        std::size_t run = std::min(frames - done, static_cast<std::size_t>(block - position % block));
        engine_.render(stereo.subspan(2 * done, 2 * run));
        done += run;
    }
}
const Synthesizer& NudgeStream::engine() const noexcept { return engine_; }
std::size_t NudgeStream::pending() const noexcept { return pending_.size(); }
}
