#pragma once

// The body language: a time-ordered stream of postures for named parts of
// the vocal system, plus behavioural markers. The planner writes it; the
// body reads it as it arrives and moves each part toward its latest target.
//
// Text form, one event per line (times and move durations in seconds):
//
//   0.512 tongue 12.80 0.350 0.080     place_cm area_cm2 move_s
//   0.512 tip    16.00 0.000 0.030     place_cm area_cm2 move_s
//   0.700 lungs  780 0 0.120           pressure_Pa unused move_s
//   0.700 pitch  118 0 0.150           intended_Hz unused move_s
//   0.300 mark breath 0.350            marker duration_s
//
// A constriction channel (tip, dorsum) with area >= kReleased exerts no
// constraint: the part is out of the way.

#include "voice/units.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace voice {

enum class Part : std::uint8_t {
  lungs,    // a: alveolar pressure (Pa)
  pitch,    // a: intended fundamental (Hz); the speaker's ear closes the loop
  glottis,  // a: arytenoid rest half-gap (cm), b: medial compression 0..1
  tongue,   // a: body constriction place from glottis (cm), b: area (cm^2)
  tip,      // a: tip/blade constriction place (cm), b: area (cm^2)
  dorsum,   // a: dorsum constriction place (cm), b: area (cm^2)
  lips,     // a: lip opening area (cm^2), b: protrusion (cm)
  velum,    // a: velopharyngeal port area (cm^2)
  floor,    // a: sublingual space volume (cm^3) under the raised tip
  jaw,      // a: gap between the incisors (cm^2); >= kReleased leaves the teeth apart
  lateral,  // a: opening of the pocket over the blade behind a lateral contact (cm^2); 0 closed
  count
};

constexpr Index kPartCount = static_cast<Index>(Part::count);
constexpr float kReleased = 8.0F;

enum class Mark : std::uint8_t {
  breath,   // take a breath for `a` seconds
  phrase,   // a prosodic phrase starts
  word,     // a word starts
  stop      // the utterance ends
};

struct PartEvent {
  double time = 0.0;
  Part part = Part::lungs;
  float a = 0.0F;
  float b = 0.0F;
  float move = 0.05F;  // seconds to (nearly) reach the target
};

struct MarkEvent {
  double time = 0.0;
  Mark mark = Mark::word;
  float a = 0.0F;
  std::string text;
};

struct Score {
  std::vector<PartEvent> parts;
  std::vector<MarkEvent> marks;
  double duration = 0.0;
};

const char* part_name(Part part);
const char* mark_name(Mark mark);
bool part_from_name(const std::string& name, Part& part);
bool mark_from_name(const std::string& name, Mark& mark);

void sort_score(Score& score);
std::string score_text(const Score& score);
bool parse_score_line(const std::string& line, Score& score);

}  // namespace voice
