#pragma once

// The speaker: reads the body language (voice::Score, written by the text
// front end: parsed phrases, breath marks, phonemes, learned timing and
// pitch) and deforms the segment-module body to satisfy it.
//
// Each part of the score drives a motor: a critically damped muscle pulling
// the part toward its latest target, whose command wanders slightly with the
// effort of holding a shape away from rest. Once per millisecond the motors
// set the support frames of the vocal tube (through the tract geometry of
// the front end's body: tongue body, tip, dorsum, lips, teeth, jaw), the
// arytenoid posture of the cord module, breath, velum, the opening under the
// tongue and the intended pitch. Sound itself comes only from the body.

#include "tract/articulation.hpp"
#include "tract/voice_body.hpp"
#include "voice/score.hpp"

#include <array>
#include <vector>
#include <cstdint>

namespace tract {

// The rest layout the speaker articulates: 48 nodes from the glottis to the
// lips (17.5 cm, the front end's oral grid), firm larynx walls for the first
// 2 cm, soft pharynx and mouth, lips over the last 1 cm.
BodyLayout speech_layout();
// The layout for the articulatory body: 60 nodes (21.9 cm, room for a
// lowered larynx and protruded lips), of which the open end uses as many as
// the geometry's length needs.
BodyLayout articulatory_layout();
// Does the score place organs (body, apex, mouth, spread, larynx, angle,
// groove)? Such a score is spoken by the articulatory body.
bool uses_organs(const voice::Score& score);

class Speaker {
 public:
  Speaker(const BodyLayout& layout, const voice::Score& score, std::uint32_t seed);

  // During a voiced closure the pharyngeal walls let go: the share of their
  // support stiffness released at most (0 turns it off).
  void set_closure_relaxation(double share) { relaxation_depth_ = share; }
  // ... and the pharyngeal walls are pulled outward by up to this (cm;
  // 0.25 lifts voiced stops from 26% to 38% voiced over 24 bdl sentences).
  void set_closure_pull(double offset) { pull_depth_ = offset; }
  // ... and the adductors press the folds together by up to this (cm of rest gap).
  void set_closure_press(double offset) { press_depth_ = offset; }
  // An opening constriction (lips, tip, dorsum) moves this many times faster
  // than it was asked to: tissue under pressure parts abruptly.
  void set_release_speed(double factor) { release_speed_ = factor; }
  void set_lung_follow(double seconds) { lung_follow_ = seconds; }
  // Scales every motor's wander (diagnostics; 1 as tuned, 0 holds steady).
  void set_wander_gain(double gain) { wander_gain_ = gain; }
  // The posterior chink while voicing (cm^2).
  void set_voicing_chink(double area) { voicing_chink_ = area; }
  // The flared open end: growth per node past the lips and its widest area.
  // Past the lip opening the articulatory body's tube is a near-closed
  // sliver of this area (cm^2).
  void set_sliver_area(double area) { sliver_area_ = area; }
  void set_flare(double growth, double widest) {
    flare_growth_ = growth;
    flare_widest_ = widest;
  }
  // The voicing posture: rest displacements of the folds' lower and upper
  // edges at the vocal processes (cm; VocalTractLab's modal 0.00495, 0.0004).
  void set_voicing_posture(double lower, double upper) {
    voicing_lower_ = lower;
    voicing_upper_ = upper;
  }

  // Streaming: a further event of the body language, arriving while the
  // body speaks. It is placed in time order among the events not yet
  // reached; an event whose time has already passed acts at once.
  void push(const voice::PartEvent& event);
  // The end of the utterance, when it is known only as the stream goes.
  void set_duration(double seconds) { score_.duration = seconds; }
  double time() const { return time_; }

  // One body step; returns the radiated sound at the step rate.
  double step();
  bool finished() const { return time_ >= score_.duration; }
  VoiceBody& body() { return body_; }
  bool articulatory() const { return articulatory_; }
  MidsagittalTract& anatomy() { return anatomy_; }
  const TractProfile& profile() const { return profile_; }
  // The organs where the motors hold them now.
  Articulators organs();

 private:
  struct Motor {
    double position = 0.0;
    double velocity = 0.0;
    double target = 0.0;
    double omega = 40.0;          // 1/s; reaches the target in about 4.5 / omega
    double rest = 0.0;            // posture held without effort
    double span = 1.0;            // distance from rest that costs one unit of effort
    double wander = 0.0;          // current deviation of the command
    double wander_scale = 0.0;    // wander per unit effort
    double wander_seconds = 0.2;  // correlation time of the wander
  };

  static constexpr int kMotorCount = 2 * static_cast<int>(voice::kPartCount);

  Motor& motor(voice::Part part, int coordinate);
  double value(voice::Part part, int coordinate);
  void set_motor(voice::Part part, int coordinate, double rest, double wander_scale, double wander_seconds,
                 double span);
  void apply(const voice::PartEvent& event);
  void advance_motors();
  void shape_body();
  // The tract's areas from the old area-profile parts (tongue, tip, dorsum,
  // lips, jaw); returns the narrowest area in the mouth (a seal).
  double shape_area_parts();
  // The tract from organ positions through the midsagittal model; returns
  // the narrowest area in the mouth.
  double shape_organs();
  void profile_to_nodes(std::vector<double>& areas, double length) const;
  double gauss();

  voice::Score score_;
  VoiceBody body_;
  std::array<Motor, kMotorCount> motors_{};
  std::size_t next_event_ = 0;
  long tick_ = 0;
  double time_ = 0.0;
  std::uint32_t random_ = 1U;
  double relaxation_depth_ = 0.0;
  double pull_depth_ = 0.25;
  double press_depth_ = 0.004;
  double release_speed_ = 1.0;
  double lung_follow_ = 0.006;  // s: the lung pressure follows its 1 kHz command over this
  // Firmer than VocalTractLab's modal (0.00495, 0.0004), the upper edge
  // pressed: over 24 bdl sentences H1-H2 in vowel middles +4.8 -> +0.1 dB
  // (bdl -3.4), jitter 2.7 -> 1.6% (bdl 1.5%), voiced stops 38 -> 47%
  // voiced (bdl 61%), vowels 82% (bdl 81%). (0.002, -0.004 was a first step.)
  // 2026-10-14: the lower edge 0.008 (from 0.004). At 0.004 the glottis
  // shut in only 60% of /iy/'s cycles (86-99% elsewhere) and the voice
  // jittered: over a0001-a0015 RAP 0.52% (p90 4.0%, 12 jumps per 100
  // cycles; the owner 0.23-0.42%, p90 1.5-2.0%, 4-5). At 0.008: /iy/ 83%
  // shut, RAP 0.32% (p90 1.5%, 3.8 jumps), H1-H2 and vowel tilt unchanged
  // (+1.7, -20.1 dB). Pressing the edges harder closed less, not more.
  double voicing_lower_ = 0.008;
  double voicing_upper_ = -0.012;  // (-0.008 before the folds' damping was raised; see Speaker::Speaker)
  // The articulatory body: geometry every kGeometryTicks control ticks, node
  // areas and the open end's length interpolated between frames.
  static constexpr int kGeometryTicks = 8;
  bool articulatory_ = false;
  double sliver_area_ = 0.02;
  MidsagittalTract anatomy_;
  TractProfile profile_;
  int geometry_tick_ = 0;
  std::vector<double> node_from_;
  std::vector<double> node_to_;
  std::vector<double> node_scratch_;
  double length_from_ = 17.5;
  double length_to_ = 17.5;
  double wander_gain_ = 1.0;
  double voicing_chink_ = 0.0002;
  double flare_growth_ = 1.6;
  double flare_widest_ = 30.0;
  double relaxation_ = 0.0;  // 0..1, how far the pharyngeal walls have let go
};

}  // namespace tract
