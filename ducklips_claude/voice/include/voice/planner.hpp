#pragma once

// From phrased phones to the body language. Each phone becomes timed
// postures for the parts it needs. Parts it does not need keep moving toward
// their own next targets, so the tongue body travels toward the next vowel
// through the consonants before it (coarticulation is overlap of gestures,
// not interpolation of sounds).

#include "voice/score.hpp"
#include "voice/text.hpp"

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace voice {

// A phone's timing given from outside (a human reading's labels and pitch
// track), so the engine can be heard with natural prosody. Zero keeps the
// planner's own value.
struct PhoneTiming {
  double duration = 0.0;  // seconds
  double pitch = 0.0;     // Hz at the middle (vowels)
};

// A quantity learned from a human reading as an intercept plus named
// factors: vowel pitch in semitones re the speaker's median
// (tools/pitch_fit.py: stress:1, nuclear, phrase_progress:2, final:end,
// after:voiceless, ...) or phone duration in log seconds
// (tools/duration_fit.py: phone:aa, final_vowel:comma, cluster, ...).
struct FactorModel {
  bool loaded = false;
  double intercept = 0.0;
  std::map<std::string, double> factors;
};

// Reads "intercept V" (or "tempo V") and "factor NAME V" lines; false if
// the file is unreadable.
bool load_factor_model(const std::string& path, FactorModel& model);

// A posture of the articulatory body's organs (tract/articulation.hpp): jaw
// angle (degrees), tongue body centre and tip (cm, the jaw's frame), lip
// aperture, protrusion and width (cm), larynx drop (cm).
struct OrganPosture {
  double jaw = 6.0;
  double body_x = 4.0;
  double body_y = 6.6;
  double tip_x = 6.9;
  double tip_y = 7.4;
  double aperture = 1.0;
  double protrusion = 0.0;
  double width = 3.2;
  double larynx = 0.0;
};

struct PlanOptions {
  double rate = 1.0;          // >1 speaks faster
  // Token-to-token variation, as a fraction (0.1: about 10%). A speaker
  // never says a sound the same way twice: within each phone, bdl's and
  // rms's log durations spread 0.36 and 0.54, ours 0.11 (consonants
  // 0.06-0.09) - a metronome. It scales the per-phone draws below.
  double variation = 0.1;
  std::uint32_t seed = 1U;    // the variation's random sequence
  double base_pitch = 120.0;  // Hz at phrase start
  double pitch_floor = 100.0;  // Hz at a falling statement end (90: creak at 79 Hz under the final /l/ of "Golf.")
  double accent = 22.0;       // Hz added by a pitch accent
  double speech_pressure = 750.0;  // Pa at 125 Hz; +10 Pa per Hz
  double lead_in = 0.10;      // seconds before the first breath
  // Stop releases (tools/identify_tune):
  double release_channel = 0.6;   // cm^2 a stop first opens to
  double voiced_release_move = 0.05;  // seconds for a voiced stop to open (d/g bursts vs bdl+rms)
  double release_move = 0.03;     // seconds for a stop or affricate to open (burst spectra; was 0.022)
  double burst_scale = 1.0;       // scales the burst-to-end interval
  double release_peel_area = 0.0;  // cm^2 a stop first peels open to (0: no peel stage)
  double release_peel_time = 0.015;  // seconds of the peel stage
  double voice_lag = 0.02;        // voicing after a stop's segment, into stress
  // Unstressed syllables are driven less and voiced more laxly.
  double unstressed_pressure = 1.0;  // scales lung pressure
  double unstressed_gap = 0.0;       // added to the fold rest gap (cm)
  // Affricates: when the hiss phase starts (fraction of the segment), and
  // how far a voiced affricate's folds part through it (NATO "Juliett":
  // 2-4 kHz noise, voicing 26 dB down, suggested 0.04 cm).
  double postalveolar_lips = 3.0;      // cm^2 lip opening for /sh zh ch jh/ (fricatives vs bdl+rms; was 1.2)
  double postalveolar_protrusion = 0.0;  // cm (was 0.5)
  double affricate_hiss_at = 0.3;    // (jh/ch vs bdl+rms: 9.13 -> 8.10 dB; was 0.4)
  double affricate_hiss_gap = 0.04;
  double voiceless_pressure = 0.85;     // lung pressure through a voiceless consonant (x speech pressure)
  double fricative_abduct_lead = 0.03;  // seconds a voiceless fricative starts parting the folds early
  double abduction_move = 0.05;         // seconds for the folds to part for a voiceless sound
  double adduction_move = 0.02;         // seconds for the folds to close in after a voiceless sound
  double adduction_lead = 0.0;
  double onset_gap = -1.0;              // after a voiceless sound, adduct first to this gap (cm; below voicing_gap: off)
  double onset_hold = 0.02;             // seconds before settling to the modal gap
  double onset_settle = 0.03;           // seconds that settling takes          // seconds the adduction starts before voice onset
  double stress_press = 0.0;           // extra fold pressing on stressed vowels, cm (0.005 left the physiological folds unable to start: phrase-initial "Alfa" was silent)
  double voicing_gap = -0.006;          // modal fold rest gap (cm; negative presses). -0.008 on the physiological folds let voicing die and not restart (9 silent vowels in a0001-a0060; 2 at -0.006)
  double vowel_compression = 0.6;        // medial compression (0..1) of stressed vowels
  double unstressed_compression = 0.6;  // medial compression (0..1) of unstressed vowels
  // Glottal settings of the noisy voiced and breathy sounds.
  double aspirate_gap = 0.06;          // /h/ fold rest gap (cm). 0.09 before: a quiet voiceless hiss (/hh/ 22%); the owner's "hello"
                                       // opens breathy and partly voiced. 0.06: /hh/ 49%, overall 36.3 -> 36.8% (one seed)
  // /v z zh/ fold rest gap (cm). 0.04 (2026-10-15, articulatory body):
  // at 0.012 the folds closed fully every cycle, the mean flow was about
  // 50 cm^3/s and /z/ barely hissed (4-7.5 kHz 31 dB under its peak; 0.025:
  // 22, 0.04: 17, 0.06: 14; bdl's hiss about level with its voicing).
  double voiced_fricative_gap = 0.04;
  double vowel_reach = 0.35;           // longest tongue movement into a vowel (s). 0.15 (fitted to an F2-speed figure) left the
                                       // formants parked in blocks in the spectrogram; 0.35 with arrival 0.8 keeps them moving
  double vowel_arrival = 0.8;          // fraction of the vowel by which the tongue arrives (0.5 before)
  double labial_jaw = 1.0;             // jaw raise (0..1) for a labial consonant
  double jaw_lead = 0.04;              // seconds a consonant's jaw raise starts early
  double pause_pressure = 0.0;         // Pa held through a minor (comma) pause (was 300: the next closure clicked)
  double phrase_fade_tail = 0.06;      // seconds past the phrase end that breath support takes to fade (was 0.12)
  double phrase_fade_from = 0.7;       // share of the last vowel before the breath starts to fade (0.45 before)
  double phrase_open_delay = 0.03;     // seconds past the phrase end before the folds part (was 0.10)
  double fricative_body_area = 0.9;    // cm^2: the tongue body's least aperture during a tip fricative (0: off)
  double fricative_hold = 0.0;         // fraction of an onset tip fricative before the vowel's tongue moves
  double breath_lead = 0.04;           // seconds before a phrase that breath pressure starts to rise (NATO onsets; was 0.08)
  double breath_rise = 0.11;           // seconds that rise takes
  double vowel_breath_lead = 0.12;     // the same before a phrase that opens on a vowel
  double velar_front = 1.8;            // cm a velar moves forward beside a front vowel (was 0.85)
  double final_devoicing = 0.0;        // fold gap a voiced fricative opens to before no vowel (0: off)
  // Phrase-final lengthening, as a factor on the whole duration. rms, re
  // the medial median: before a comma pause ax 3.0x, iy 2.2, er 1.8, ow
  // 1.8; at the sentence end ax 2.5, iy 1.6, er 1.4, ow 1.5; consonants
  // 1.2-2.0x at either (n 2.0, t 1.5, d 1.7-1.9, s 1.2-1.6).
  double final_vowel_scale = 1.8;      // last vowel before a sentence end
  FactorModel pitch_model;             // learned vowel pitch (used when loaded)
  FactorModel duration_model;          // learned phone durations (used when loaded)
  double pitch_median = 122.0;         // Hz the learned model's semitones are relative to (bdl: median of
                                       // 120 sentence medians 122.1, IQR 118.5-125.6; was 110)
  // Least squares shrinks the learned factors toward the mean: their
  // held-out correlation with bdl is 0.55, so the predicted contour carries
  // about 55% of his spread (1.3 of 2.3 semitones). The factors are widened
  // by this; with the per-syllable variation (sd ~0.5 st) the spread is
  // about his again. (Enhanced fields: bdl's harmonics sweep with the
  // intonation where ours ran nearly flat.)
  double pitch_spread = 1.5;
  // A phrase-final voiceless stop keeps its breath until its release (true
  // for the segment-module body; the front end's own body faded it through
  // the vowel to avoid a 7x burst at "Foxtrot.").
  bool final_stop_breath = false;
  double lowest_pitch = 80.0;          // Hz: the folds break below about 78 (probe: 75 Hz jitters 14%)
  double final_fall = 0.9;             // a statement's last vowel falls to this share of its own pitch
  double comma_pitch = 0.93;           // pitch through the last vowel before a comma, re its own (was 1.06)
  // Sloped syllable pitch targets (0 peak: one level target per vowel).
  // Stressed: +peak at 40% into the vowel, then -fall by the next vowel;
  // unstressed: -drift. The flat staircase was the early-TTS look.
  double accent_peak = 0.05;  // (the planned pitch already holds the accent: +30% put stressed NATO vowels at 200-250 Hz)
  double accent_fall = 0.25;
  double unstressed_drift = 0.06;
  // Gliding intonation (replaces the sloped targets when on): the pitch
  // heads continuously from one vowel's anchor to the next, so it turns
  // only where the anchors turn. Stressed vowels, and each phrase's first
  // and last, keep their planned pitch; the vowels between take the line
  // between their neighbours, lowered by unstressed_dip semitones. (Sloped
  // targets rose and fell 5-7 semitones on nearly every syllable; the
  // owner heard a warble.)
  bool pitch_glide = true;
  // The articulatory body: emit organ positions (body, apex, mouth, spread,
  // larynx, angle, groove, lateral) from these vowel postures instead of
  // the area-profile parts (`say --organs FILE`).
  bool organs = false;
  std::map<std::string, OrganPosture> organ_postures;
  // /s z/ are made with the lips open and unrounded: at least this lip area
  // (cm^2) and at most this protrusion through them, and a labial released
  // into one opens to it. (Released straight to the next vowel's rounding,
  // "jumps over" put /s/ behind 0.3 cm^2 protruded lips: its peak fell
  // from 6.5 to 5 kHz and the owner heard "jumpsh".)
  double sibilant_lips = 1.5;
  double lip_press = 0.01;             // cm^2 (as a negative signed area) the lips press past touching in /p b m/
  double sibilant_protrusion = 0.1;
  double unstressed_dip = 0.5;
  double comma_vowel_scale = 2.4;      // last vowel before a minor boundary (continuation)
  double final_consonant_scale = 1.6;  // consonants after it, either boundary
  std::vector<PhoneTiming> timing;     // per phone, in order (empty: planned)
  std::vector<double> pauses;          // seconds after each phrase (empty: planned)
  double final_vowel_floor_ms = 110.0; // no last vowel shorter (rms final ax 110 ms, 2.5x medial)
  double sonorant_devoicing = 0.5;     // share of a /w r l j/ after an aspirated stop that is voiceless
  bool round_ahead = true;             // lips round through the consonant before /w r/
  double final_relax_gap = -0.008;     // fold gap the last vowel of a phrase relaxes to (0.006: the folds of "Golf." double-pulsed through the /l/)
};

struct TimedPhone {
  std::string name;
  std::string phoneme;        // before allophonic rules (a flap stays /t/ or /d/)
  std::int32_t stress = -1;
  double start = 0.0;
  double duration = 0.0;
  bool word_start = false;
  bool word_end = false;
  bool accent = false;
  bool nuclear = false;
  bool phrase_final = false;  // in the last syllable of its phrase
  double pitch = 0.0;         // intended Hz at the middle (vowels)
  std::string word;
  // Standard normal draws for this token (PlanOptions::variation scales them).
  double vary_duration = 0.0;
  double vary_pitch = 0.0;
  double vary_place = 0.0;
  double vary_area = 0.0;
  double vary_lips = 0.0;
  double vary_effort = 0.0;
};

struct Plan {
  std::vector<TimedPhone> phones;
  Score score;
};

Plan plan_speech(const std::vector<Phrase>& phrases, const PlanOptions& options);

// Overrides phone parameters from a text file, one phone per line:
//   NAME key value [key value ...]   keys: place area floor inherent minimum
//                                          vplace varea vlips vprotrusion
[[nodiscard]] bool load_phone_overrides(const std::string& path);
// Vowel postures for the articulatory body: lines "PHONE jaw=.. tx=.. ty=..
// tipx=.. tipy=.. ap=.. pro=.. lw=.. drop=..". Returns false if unreadable.
bool load_organ_postures(const std::string& path, PlanOptions& options);

}  // namespace voice
