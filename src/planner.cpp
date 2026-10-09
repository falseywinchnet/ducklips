#include "ducklips/planner.hpp"
#include "ducklips/connected_speech.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <sstream>

namespace ducklips {
namespace {

// Appends cues in seconds, converting to samples once.
class CueWriter final {
public:
    CueWriter(std::vector<Cue>& cues, SynthesisRate rate) : cues_(cues), rate_(rate) {}

    void scalar(double time, CueKind kind, double value, double angular_rate) {
        Cue cue{};
        cue.time = samples_from_seconds(rate_, std::max(time, 0.0));
        cue.kind = kind;
        cue.value = value;
        cue.angular_rate = angular_rate;
        cues_.push_back(cue);
    }

    void vowel(double time, const VowelShape& shape, double angular_rate) {
        Cue cue{};
        cue.time = samples_from_seconds(rate_, std::max(time, 0.0));
        cue.kind = CueKind::vowel;
        cue.vowel = shape;
        cue.angular_rate = angular_rate;
        cues_.push_back(cue);
    }

    void gesture(double time, const ConstrictionSpec& spec, double activation, double angular_rate) {
        settle_release(spec.slot, time);
        Cue cue{};
        cue.time = samples_from_seconds(rate_, std::max(time, 0.0));
        cue.kind = CueKind::constriction;
        cue.slot = spec.slot;
        cue.gesture = gesture_from_spec(spec);
        cue.value = activation;
        cue.angular_rate = angular_rate;
        cues_.push_back(cue);
    }

    // A release is held back until the articulator's next gesture is known:
    // if that gesture begins before the release time, the articulator moves
    // straight on (/z/ to /d/, /n/ to /d/) and the release is dropped.
    void release(double time, ConstrictionSlot slot, double angular_rate) {
        release_to(time, slot, 0.0, angular_rate);
    }

    // Drives a slot's activation toward a value (0 opens it fully).
    void release_to(double time, ConstrictionSlot slot, double value, double angular_rate) {
        std::size_t index = static_cast<std::size_t>(slot);
        flush_release(slot);
        Cue cue{};
        cue.time = samples_from_seconds(rate_, std::max(time, 0.0));
        cue.kind = CueKind::release;
        cue.slot = slot;
        cue.value = value;
        cue.angular_rate = angular_rate;
        pending_[index] = cue;
        has_pending_[index] = true;
    }

    void flush_all() {
        for (std::size_t index = 0; index < constriction_slot_count; ++index) {
            flush_release(static_cast<ConstrictionSlot>(index));
        }
    }

private:
    void settle_release(ConstrictionSlot slot, double onset) {
        std::size_t index = static_cast<std::size_t>(slot);
        if (!has_pending_[index]) {
            return;
        }
        // A release followed within 30 ms by the same articulator's next
        // gesture is dropped too: the articulator slides on (/z/ to /d/)
        // instead of opening and slamming shut again, which cut a whisper's
        // airflow off with a thump.
        SampleCount onset_sample = samples_from_seconds(rate_, std::max(onset - 0.030, 0.0));
        if (pending_[index].time >= onset_sample) {
            has_pending_[index] = false;
            return;
        }
        flush_release(slot);
    }

    void flush_release(ConstrictionSlot slot) {
        std::size_t index = static_cast<std::size_t>(slot);
        if (has_pending_[index]) {
            cues_.push_back(pending_[index]);
            has_pending_[index] = false;
        }
    }

    std::vector<Cue>& cues_;
    SynthesisRate rate_;
    std::array<Cue, constriction_slot_count> pending_{};
    std::array<bool, constriction_slot_count> has_pending_{};
};

// Deterministic variation in [-1, 1] for one phone of one voice.
double variation(std::uint64_t seed, std::size_t index) noexcept {
    std::uint64_t mixed = seed * 0x9E3779B97F4A7C15ULL + static_cast<std::uint64_t>(index) * 0xBF58476D1CE4E5B9ULL;
    mixed = (mixed ^ (mixed >> 31)) * 0x94D049BB133111EBULL;
    mixed = mixed ^ (mixed >> 29);
    double unit = static_cast<double>(mixed >> 11) * (1.0 / 9007199254740992.0);
    double value = 2.0 * unit - 1.0;
    return value;
}

bool voiced_phone(Phone phone) noexcept {
    bool voiced = phone_spec(phone).voiced;
    return voiced;
}

bool obstruent(PhoneClass kind) noexcept {
    bool result = kind == PhoneClass::stop || kind == PhoneClass::fricative || kind == PhoneClass::affricate
        || kind == PhoneClass::aspirate;
    return result;
}

double slot_rate(ConstrictionSlot slot) noexcept {
    switch (slot) {
    case ConstrictionSlot::lips: return 150.0;
    case ConstrictionSlot::tongue_tip: return 165.0;
    case ConstrictionSlot::tongue_body: return 110.0;
    }
    return 120.0;
}

struct Accent final {
    double peak{};   // s
    double height{}; // fractional rise
};

struct Intonation final {
    double base{};
    double declination{};
    double start{};
    double end{};
    bool question{};
    std::vector<Accent> accents{};

    double pitch(double time) const noexcept {
        double span = std::max(end - start, 0.05);
        double progress = std::clamp((time - start) / span, 0.0, 1.0);
        double value = base * (1.0 + declination * (0.5 - progress));
        // Falling accents (MEASURED in the architect's reading: every
        // stressed word starts high and falls within itself, "sea" 140 -> 105
        // Hz, "shore" 120 -> 80, resetting at the next accent). Each accent
        // rises quickly into its vowel onset and decays over about 130 ms;
        // the highest accent at any moment wins. A bell centred mid-vowel was
        // heard as an Irish lilt, and a hat with a plateau as monotone.
        double raise = 0.0;
        for (std::size_t index = 0; index < accents.size(); ++index) {
            double offset = time - accents[index].peak;
            double shape = offset < 0.0 ? std::exp(-(offset / 0.03) * (offset / 0.03)) : std::exp(-offset / 0.13);
            raise = std::max(raise, accents[index].height * shape);
        }
        value *= 1.0 + raise;
        if (question) {
            double ramp = std::clamp((time - (end - 0.30)) / 0.30, 0.0, 1.0);
            value *= 1.0 + 0.45 * ramp * ramp;
        } else {
            double ramp = std::clamp((time - (end - 0.25)) / 0.25, 0.0, 1.0);
            value *= 1.0 - 0.10 * ramp;
        }
        return value;
    }
};

const char* digit_name(char digit) noexcept {
    switch (digit) {
    case '0': return "zero";
    case '1': return "one";
    case '2': return "two";
    case '3': return "three";
    case '4': return "four";
    case '5': return "five";
    case '6': return "six";
    case '7': return "seven";
    case '8': return "eight";
    default: return "nine";
    }
}

}

void UtterancePlanner::prepare(const VoiceProfile& voice, const VowelLibrary& vowels, SynthesisRate rate) {
    voice_ = voice;
    vowels_ = vowels;
    rate_ = rate;
}

std::size_t UtterancePlanner::load_dictionary(const std::string& path) {
    std::size_t added = lexicon_.load_cmu_file(path);
    return added;
}

double UtterancePlanner::base_duration(const std::vector<PhoneToken>& tokens, std::size_t index) const {
    const PhoneToken& token = tokens[index];
    if (token.duration_override > 0.0) {
        return token.duration_override;
    }
    const PhoneSpec& spec = phone_spec(token.phone);
    bool stressed = token.stress == 1 || token.stress == 2;
    double duration = 0.07;
    switch (spec.kind) {
    case PhoneClass::vowel:
        if (token.phone == Phone::AX) {
            duration = 0.05;
        } else if (spec.glide_to != Phone::count) {
            duration = stressed ? 0.19 : 0.12;
        } else if (spec.tense) {
            duration = stressed ? 0.15 : 0.09;
        } else if (token.phone == Phone::AE) {
            duration = stressed ? 0.13 : 0.08;
        } else {
            duration = stressed ? 0.10 : 0.06;
        }
        if (token.stress == 2) {
            duration *= 0.85;
        }
        if (index + 1 < tokens.size() && !is_vowel(tokens[index + 1].phone)) {
            if (voiced_phone(tokens[index + 1].phone)) {
                duration *= 1.15;
            } else {
                duration *= 0.90;
            }
        }
        break;
    case PhoneClass::stop:
        duration = spec.voiced ? 0.060 : 0.075;
        break;
    case PhoneClass::fricative:
        if (spec.voiced) {
            duration = 0.070;
        } else if (token.phone == Phone::S || token.phone == Phone::SH) {
            duration = 0.110;
        } else {
            duration = 0.095;
        }
        break;
    case PhoneClass::affricate:
        duration = 0.110;
        break;
    case PhoneClass::nasal:
        duration = 0.065;
        break;
    case PhoneClass::liquid:
        duration = 0.060;
        break;
    case PhoneClass::glide:
        duration = 0.055;
        break;
    case PhoneClass::aspirate:
        duration = 0.060;
        break;
    }
    if (spec.kind != PhoneClass::vowel) {
        bool before = index > 0 && !is_vowel(tokens[index - 1].phone);
        bool after = index + 1 < tokens.size() && !is_vowel(tokens[index + 1].phone);
        if (before || after) {
            duration *= 0.8;
        }
        // A stop needs a closure long enough to build pressure and be heard.
        if (spec.kind == PhoneClass::stop || spec.kind == PhoneClass::affricate) {
            duration = std::max(duration, 0.058);
        }
    } else {
        // Short unstressed vowels still need a few voiced cycles.
        duration = std::max(duration, 0.055);
    }
    return duration;
}

UtterancePlan UtterancePlanner::plan_phrases(const std::vector<Phrase>& phrases) const {
    UtterancePlan plan{};
    CueWriter writer(plan.cues, rate_);
    const double speed = voice_.speaking_rate;
    const double lead = 0.6;
    // Shortest pause that holds an audible inhale (s).
    const double breath_room = 0.30;
    double time = lead;
    double previous_end = 0.0;

    // The tongue starts in the posture of the first vowel.
    Phone first_vowel = Phone::AX;
    bool found_vowel = false;
    for (std::size_t phrase = 0; phrase < phrases.size() && !found_vowel; ++phrase) {
        for (std::size_t index = 0; index < phrases[phrase].tokens.size(); ++index) {
            if (is_vowel(phrases[phrase].tokens[index].phone)) {
                first_vowel = phrases[phrase].tokens[index].phone;
                found_vowel = true;
                break;
            }
        }
    }
    writer.vowel(0.0, vowels_.shape(first_vowel), 0.0);
    writer.scalar(0.0, CueKind::pitch, voice_.base_pitch, 0.0);

    std::size_t global_index = 0;
    double last_pause = 0.0;
    for (std::size_t phrase_index = 0; phrase_index < phrases.size(); ++phrase_index) {
        const Phrase& phrase = phrases[phrase_index];
        const std::vector<PhoneToken>& tokens = phrase.tokens;
        const std::size_t count = tokens.size();
        if (count == 0) {
            continue;
        }
        std::vector<double> durations(count, 0.0);
        std::size_t last_vowel = count;
        for (std::size_t index = 0; index < count; ++index) {
            double wobble = 1.0 + 0.05 * variation(voice_.seed, global_index + index);
            durations[index] = base_duration(tokens, index) * tokens[index].duration_scale * wobble / speed;
            if (tokens[index].flap) {
                durations[index] = 0.024 / speed;
            }
            if (is_vowel(tokens[index].phone)) {
                last_vowel = index;
            }
        }
        for (std::size_t index = last_vowel; index < count; ++index) {
            if (tokens[index].duration_override <= 0.0) {
                durations[index] *= 1.4;
            }
        }
        const double phrase_start = time;
        double phrase_end = phrase_start;
        for (std::size_t index = 0; index < count; ++index) {
            phrase_end += durations[index];
        }

        // Breathing: inhale through nose and mouth, then pressurize. The
        // inhale must fit in the pause after the previous phrase; a short
        // pause (a comma) gets no breath, only the pressure back.
        const PhoneSpec& opening = phone_spec(tokens[0].phone);
        const double room = phrase_index == 0 ? lead : phrase_start - previous_end - 0.05;
        const bool inhale = voice_.audible_breath && room >= breath_room;
        if (inhale) {
            const double begin = std::min(0.50, room);
            writer.scalar(phrase_start - begin, CueKind::glottis, voice_.abducted_area * 1.6, 60.0);
            writer.scalar(phrase_start - begin, CueKind::velum, 0.8 * square_centimetre, 50.0);
            writer.scalar(phrase_start - begin + 0.04, CueKind::lungs, -110.0, 30.0);
            // Output modulation: the inhale is heard, but under the speech.
            writer.scalar(phrase_start - begin - 0.02, CueKind::loudness, voice_.breath_loudness, 60.0);
            writer.scalar(phrase_start - 0.12, CueKind::loudness, 1.0, 60.0);
            // Pressure is up before the first sound: an utterance-initial
            // voiced fricative ("the") needs voicing and airflow at once.
            // A voiced first sound: pressure rises early and the folds come
            // together soon after, so the rising pressure does not blow an
            // audible breath through the open glottis. A voiceless first
            // sound keeps the glottis open (closing it would voice a "heh"
            // before "hello"), so its pressure rises later.
            // For a voiced first sound the folds close as the pressure rises
            // (quickly, from 0.12 s before), so there is neither an audible
            // breath nor an "uh" before it.
            // A voiceless first sound gets its pressure just before it: earlier,
            // the open glottis lets out a long "hhh" (heard as hesitancy).
            const double pressurize = opening.voiced ? phrase_start - std::min(0.12, begin - 0.10)
                                                     : phrase_start - 0.05;
            writer.scalar(pressurize, CueKind::lungs, voice_.lung_pressure, 60.0);
            writer.scalar(pressurize, CueKind::glottis,
                          opening.voiced ? voice_.adducted_area : voice_.abducted_area, 60.0);
            // The velum seals with the pressure rise unless the phrase opens on
            // a nasal; voicing behind closed lips would otherwise hum.
            if (opening.kind != PhoneClass::nasal) {
                writer.scalar(pressurize, CueKind::velum, 0.0, 60.0);
            }
        } else {
            writer.scalar(phrase_start - std::min(0.20, std::max(room, 0.08)), CueKind::lungs, voice_.lung_pressure, 40.0);
        }
        if (opening.kind != PhoneClass::nasal) {
            writer.scalar(phrase_start - 0.10, CueKind::velum, 0.0, 45.0);
        }
        if (opening.voiced) {
            writer.scalar(phrase_start - 0.09, CueKind::glottis, voice_.adducted_area, 70.0);
        }

        Intonation intonation{};
        intonation.base = voice_.base_pitch;
        intonation.declination = voice_.declination;
        intonation.start = phrase_start;
        intonation.end = phrase_end;
        intonation.question = phrase.question;
        std::size_t primary_accents = 0;

        double start = phrase_start;
        double last_vowel_middle = phrase_start - 0.3;
        double last_vowel_end = phrase_start - 0.1;
        for (std::size_t index = 0; index < count; ++index) {
            const PhoneToken& token = tokens[index];
            const PhoneSpec& spec = phone_spec(token.phone);
            const double duration = durations[index];
            const double end = start + duration;
            plan.phones.push_back(PlannedPhone{.phone = token.phone, .stress = token.stress, .start = start, .end = end});
            const bool next_exists = index + 1 < count;
            const Phone next_phone = next_exists ? tokens[index + 1].phone : Phone::AX;
            const bool next_voiced = next_exists && voiced_phone(next_phone);

            if (spec.kind == PhoneClass::vowel) {
                // Vowel-to-vowel coarticulation: the tongue body heads for this
                // vowel partway through the preceding consonants. A tip or
                // dorsum consonant drags the body by geometry while it holds.
                double consonant_span = std::max(start - last_vowel_end, 0.0);
                double when = std::max(start - 0.055 - 0.5 * consonant_span, last_vowel_middle);
                // Before /r/ the American open-mid back vowel rises toward [o]
                // ("shore", "four": F1 about 500 Hz in the human reading).
                Phone posture = token.phone;
                if (token.phone == Phone::AO && next_exists && next_phone == Phone::R) {
                    posture = Phone::OW;
                }
                writer.vowel(when, vowels_.shape(posture), 42.0 * speed);
                if (spec.glide_to != Phone::count) {
                    writer.vowel(start + 0.45 * duration, vowels_.shape(spec.glide_to), 30.0 * speed);
                }
                if (token.stress == 1) {
                    // Vocal effort rises on an accented syllable: more lung
                    // pressure and firmer adduction, not a louder output.
                    writer.scalar(start - 0.05, CueKind::effort, 1.12, 40.0);
                    writer.scalar(end - 0.02, CueKind::effort, 1.0, 30.0);
                    // The peak sits at the stressed vowel's onset and falls
                    // through it (American H*); peaks at 55% of the vowel held
                    // the pitch up and fell late, a mark of Irish English.
                    double height = voice_.accent_rise * std::pow(0.92, static_cast<double>(primary_accents));
                    intonation.accents.push_back(Accent{.peak = start + 0.03 * duration, .height = height});
                    primary_accents += 1;
                } else if (token.stress == 2) {
                    intonation.accents.push_back(Accent{.peak = start + 0.03 * duration, .height = 0.4 * voice_.accent_rise});
                }
                last_vowel_middle = start + 0.5 * duration;
                last_vowel_end = end;
                start = end;
                continue;
            }

            if (obstruent(spec.kind) && !spec.voiced && !token.flap) {
                double abduction = voice_.abducted_area;
                if (spec.kind == PhoneClass::aspirate) {
                    abduction *= 0.8;
                }
                // Fricatives part the folds well before the noise so voicing
                // stops at its onset ("sea", not "zea"); stops devoice within
                // their closure.
                double abduct_lead = spec.kind == PhoneClass::fricative ? 0.035 : 0.012;
                writer.scalar(start - abduct_lead, CueKind::glottis, abduction, 110.0);
                if (next_voiced) {
                    // Voice onset time: aspiration after a stressed voiceless
                    // stop; otherwise the folds start closing at the release.
                    // After a fricative they close during its last 40 ms, so
                    // voicing starts under the constriction and the tongue then
                    // releases into the vowel (the human "shore"). Released
                    // first, the open tract passed a weak breath: a near-silent
                    // gap between "sh" and the vowel.
                    double lag = spec.kind == PhoneClass::fricative ? -0.040 : -0.010;
                    if (spec.kind == PhoneClass::stop && is_vowel(next_phone) && tokens[index + 1].stress == 1) {
                        lag = 0.025;
                    }
                    writer.scalar(end + lag, CueKind::glottis, voice_.adducted_area, 160.0);
                }
            }

            if (spec.kind == PhoneClass::stop && spec.voiced && !token.flap && next_exists
                && is_vowel(next_phone)) {
                // An English voiced stop before a vowel mostly devoices during
                // its closure (MEASURED in the architect's reading: no voice
                // bar in "duck" or "dog"): the folds part a little, pressure
                // behind the closure rises toward lung pressure, the release
                // bursts, and voicing returns within about 15 ms (short lag).
                double partial = voice_.adducted_area + 0.5 * (voice_.abducted_area - voice_.adducted_area);
                writer.scalar(start - 0.015, CueKind::glottis, partial, 160.0);
                writer.scalar(end - 0.012, CueKind::glottis, voice_.adducted_area, 160.0);
            }
            if ((spec.kind == PhoneClass::fricative || spec.kind == PhoneClass::affricate) && spec.voiced
                && !token.flap) {
                // A voiced fricative needs airflow through a voicing glottis:
                // the folds part a little, keep vibrating, and the jet through
                // the constriction becomes turbulent.
                double partial = voice_.adducted_area + 0.3 * (voice_.abducted_area - voice_.adducted_area);
                writer.scalar(start - 0.012, CueKind::glottis, partial, 110.0);
                // Before an obstruent the next sound sets the glottis itself;
                // closing it here overrode a following stop's opening ("is
                // duck" voiced straight through its /d/: "bucklips").
                const bool obstruent_next = next_exists && obstruent(phone_spec(next_phone).kind);
                if (!obstruent_next) {
                    writer.scalar(end - 0.010, CueKind::glottis, voice_.adducted_area, 140.0);
                }
            }

            const ConstrictionSpec& primary = spec.primary;
            const ConstrictionSpec& secondary = spec.secondary;
            // The first consonant of a phrase is formed before voicing starts.
            const double early = index == 0 ? 0.08 : 0.0;
            switch (spec.kind) {
            case PhoneClass::stop:
            case PhoneClass::nasal:
                if (token.flap) {
                    // A tap: the tongue tip strikes the ridge once without holding.
                    writer.gesture(start - 0.012, primary, 1.05, 320.0);
                    writer.release(end - 0.004, primary.slot, 260.0);
                    break;
                }
                // A stop opens at about 0.6 of the closing speed (0.3 sounded like /l/) (measured
                // release rates of tens of cm^2/s): the jet through the
                // widening gap makes a burst lasting several milliseconds. The
                // closure is centred on its slot: the slow release starts before
                // the slot ends, so the closure does not eat the next vowel.
                if (spec.kind == PhoneClass::stop) {
                    // Two-stage release: the pressed tissue first eases back to
                    // bare contact, then the closure opens from rest. Released
                    // in one motion, the articulator crossed contact at full
                    // speed and the gap opened in about a millisecond, leaving
                    // no burst (MEASURED: /d/ 0 -> 1.6 cm^2 in 6 ms).
                    // The closure approaches at 0.7 of the slot's speed: sealed in
                    // a few milliseconds it cut a whisper's airflow off at once,
                    // a thump louder than the vowels.
                    writer.gesture(start - 0.015 - early, primary, 1.25, 0.7 * slot_rate(primary.slot));
                    writer.release_to(end - 0.040, primary.slot, 1.0, slot_rate(primary.slot));
                    writer.release(end - 0.012, primary.slot, 0.35 * slot_rate(primary.slot));
                } else {
                    writer.gesture(start - 0.020 - early, primary, 1.25, slot_rate(primary.slot));
                    writer.release(end, primary.slot, slot_rate(primary.slot));
                }
                break;
            case PhoneClass::fricative:
                writer.gesture(start - 0.025 - early, primary, 1.0, 120.0);
                // The constriction holds until the folds have closed (voicing
                // starts under it), then releases into the vowel.
                writer.release(end - 0.005, primary.slot, 150.0);
                break;
            case PhoneClass::affricate: {
                ConstrictionSpec closure = primary;
                closure.area_cm2 = 0.0;
                closure.obstacle_cm = -1.0;
                writer.gesture(start - 0.020, closure, 1.25, slot_rate(primary.slot));
                writer.gesture(start + 0.5 * duration, primary, 1.0, 140.0);
                writer.release(end, primary.slot, 110.0);
                break;
            }
            case PhoneClass::liquid:
            case PhoneClass::glide:
                if (token.phone == Phone::R && index > 0 && is_vowel(tokens[index - 1].phone)) {
                    // American post-vocalic /r/ colours the vowel before it from
                    // about its middle (the human "shore": F3 falls steadily from
                    // 2.2 to 1.6 kHz across the word). Starting 80 ms early left
                    // a long flat vowel and a light /r/, heard as Irish.
                    writer.gesture(std::min(start - 0.080, last_vowel_middle - 0.020), primary, 1.0, 40.0);
                    // At the end of a phrase the /r/ is held while the voice
                    // fades; released during voicing it added an "uh" (F3 back
                    // up to 2.4 kHz): a released Scottish or Irish r.
                    writer.release(next_exists ? end + 0.005 : end + 0.250, primary.slot, 60.0);
                    break;
                }
                writer.gesture(start - 0.035 - early, primary, 1.0, 75.0);
                writer.release(end + 0.005, primary.slot, 60.0);
                break;
            case PhoneClass::aspirate:
            case PhoneClass::vowel:
                break;
            }
            if (secondary.present) {
                writer.gesture(start - 0.035, secondary, 1.0, 75.0);
                writer.release(end + 0.005, secondary.slot, 60.0);
            }
            if (spec.kind == PhoneClass::nasal) {
                writer.scalar(start - 0.050, CueKind::velum, 1.1 * square_centimetre, 55.0);
                const PhoneClass next_kind = next_exists ? phone_spec(next_phone).kind : PhoneClass::vowel;
                bool next_nasal = next_exists && next_kind == PhoneClass::nasal;
                bool next_obstruent = next_exists && obstruent(next_kind);
                if (next_obstruent) {
                    // The velum seals before a following obstruent so that it can
                    // build pressure: sealed by the nasal's end, or the stop's
                    // closure is a nasal and "dog" turns into "nog".
                    writer.scalar(std::max(start + 0.4 * duration, end - 0.050), CueKind::velum, 0.0, 220.0);
                } else if (!next_nasal) {
                    writer.scalar(end + 0.015, CueKind::velum, 0.0, 55.0);
                }
            }
            start = end;
        }
        global_index += count;

        writer.flush_all();
        writer.scalar(phrase_end - 0.25, CueKind::effort, 0.88, 20.0);
        writer.scalar(phrase_end + 0.20, CueKind::effort, 1.0, 20.0);
        // The phrase trails off: pressure falls, the folds part. Before a
        // short pause without a breath the pressure only eases.
        const bool last_phrase = phrase_index + 1 == phrases.size();
        const bool breath_follows = last_phrase || phrase.pause_after - 0.05 >= breath_room;
        (void)breath_follows; // the pressure falls at every pause; only the inhale needs room
        // The voice fades through the last syllable (the human reading falls
        // about 20 dB over its final 300 ms).
        writer.scalar(phrase_end - 0.08, CueKind::lungs, 0.0, 15.0);
        writer.scalar(phrase_end + 0.03, CueKind::glottis, voice_.abducted_area, 40.0);
        for (double moment = phrase_start - 0.10; moment <= phrase_end + 0.05; moment += 0.01) {
            writer.scalar(moment, CueKind::pitch, intonation.pitch(moment), 120.0);
        }
        time = phrase_end + phrase.pause_after;
        last_pause = phrase.pause_after;
        previous_end = phrase_end;

        std::ostringstream text{};
        if (!plan.transcription.empty()) {
            text << " | ";
        }
        for (std::size_t index = 0; index < count; ++index) {
            if (index > 0) {
                text << ' ';
            }
            text << phone_spec(tokens[index].phone).name;
            if (tokens[index].flap) {
                text << '~';
            }
            if (tokens[index].stress >= 0) {
                text << tokens[index].stress;
            }
        }
        plan.transcription += text.str();
    }
    double finish = time - last_pause + 0.35;
    writer.scalar(finish, CueKind::finish, 0.0, 0.0);
    plan.duration = finish;
    std::stable_sort(plan.cues.begin(), plan.cues.end(), CueEarlier{});
    return plan;
}

UtterancePlan UtterancePlanner::plan_text(const std::string& text) const {
    std::vector<Phrase> phrases{};
    Phrase current{};
    std::vector<std::string> guessed{};
    std::vector<std::string> events{};
    std::vector<SpokenWord> words{};
    std::string word{};
    std::size_t position = 0;
    while (position <= text.size()) {
        char letter = '\0';
        if (position < text.size()) {
            letter = text[position];
        }
        bool word_letter = std::isalnum(static_cast<unsigned char>(letter)) != 0 || letter == '\'';
        if (word_letter) {
            word.push_back(letter);
            position += 1;
            continue;
        }
        // Keep an explicitly known hyphenated compound as one lexical word.
        // Unknown compounds retain the ordinary word-boundary behavior.
        if (letter == '-' && !word.empty()) {
            std::size_t end = position;
            std::size_t known_end = position;
            while (end + 1 < text.size() && text[end] == '-'
                   && std::isalnum(static_cast<unsigned char>(text[end + 1]))) {
                end += 1;
                while (end < text.size()
                       && (std::isalnum(static_cast<unsigned char>(text[end])) || text[end] == '\'')) {
                    end += 1;
                }
                if (lexicon_.pronounce(word + text.substr(position, end - position)).from_dictionary) {
                    known_end = end;
                }
            }
            if (known_end > position) {
                word += text.substr(position, known_end - position);
                position = known_end;
                continue;
            }
        }
        if (!word.empty()) {
            bool numeric = true;
            for (std::size_t index = 0; index < word.size(); ++index) {
                if (!std::isdigit(static_cast<unsigned char>(word[index]))) {
                    numeric = false;
                }
            }
            std::vector<std::string> spoken{};
            if (numeric) {
                for (std::size_t index = 0; index < word.size(); ++index) {
                    spoken.push_back(digit_name(word[index]));
                }
            } else {
                spoken.push_back(word);
            }
            for (std::size_t index = 0; index < spoken.size(); ++index) {
                WordPronunciation pronunciation = lexicon_.pronounce(spoken[index]);
                if (!pronunciation.from_dictionary) {
                    guessed.push_back(pronunciation.word);
                }
                words.push_back(SpokenWord{.word = pronunciation.word, .phones = pronunciation.phones});
            }
            word.clear();
        }
        bool minor_break = letter == ',' || letter == ';' || letter == ':';
        bool major_break = letter == '.' || letter == '!' || letter == '?' || letter == '\n' || letter == '\0';
        if ((minor_break || major_break) && !words.empty()) {
            current.tokens = join_words(words, events);
            current.pause_after = minor_break ? 0.22 : 0.40;
            current.question = letter == '?';
            phrases.push_back(current);
            current = Phrase{};
            words.clear();
        }
        position += 1;
    }
    UtterancePlan plan = plan_phrases(phrases);
    plan.guessed_words = guessed;
    plan.transitions = events;
    return plan;
}

UtterancePlan UtterancePlanner::plan_phonemes(const std::string& phonemes) const {
    std::vector<Phrase> phrases{};
    Phrase current{};
    std::istringstream stream(phonemes);
    std::string symbol{};
    while (stream >> symbol) {
        if (symbol == "|") {
            if (!current.tokens.empty()) {
                current.pause_after = 0.3;
                phrases.push_back(current);
                current = Phrase{};
            }
            continue;
        }
        double override_seconds = 0.0;
        std::size_t colon = symbol.find(':');
        if (colon != std::string::npos) {
            override_seconds = std::atof(symbol.substr(colon + 1).c_str());
            symbol = symbol.substr(0, colon);
        }
        Phone phone = Phone::AX;
        int stress = -1;
        if (parse_phone(symbol, phone, stress)) {
            current.tokens.push_back(PhoneToken{.phone = phone, .stress = stress, .duration_override = override_seconds});
        }
    }
    if (!current.tokens.empty()) {
        phrases.push_back(current);
    }
    UtterancePlan plan = plan_phrases(phrases);
    return plan;
}

}
