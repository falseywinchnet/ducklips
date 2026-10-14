#include "voice/geometry.hpp"

#include <algorithm>
#include <cmath>

namespace voice {

namespace {

constexpr double kOralLength = kOralSections * kSectionLength;
// Epilarynx: the narrow tube from the folds to the aryepiglottic rim. Its
// own resonance near 3 kHz gathers F3 and F4 into the speaker's cluster.
// With its top at 3.0 cm (narrow for about 2.2 cm, then widening) the
// cluster sat near 4 kHz and vowels were 5-12 dB weak at 2.75-3.5 kHz
// against bdl; at 3.5 cm the vowel band balance matches within 3 dB.
constexpr double kEpilarynxTop = 3.5;
constexpr double kEpilarynxArea = 0.45;
constexpr double kJawNarrowing = 0.45;  // area loss at full jaw closure

double area_from_root(double root) {
  return root <= 0.0 ? 0.0 : root * root;
}

}  // namespace

TractLayout tract_areas(const TractShape& shape, std::array<double, kOralSections>& areas) {
  double place = shape.tongue_place;
  double tongue_area = std::max(shape.tongue_area, 0.05);
  double tip_place = shape.tip_place;
  double tip_root = shape.tip_root;
  double dorsum_place = shape.dorsum_place;
  double dorsum_root = shape.dorsum_root;
  double lip_root = shape.lip_root;
  double protrusion = std::clamp(shape.protrusion, 0.0, 1.2);
  double lip_length = 0.9 + protrusion;
  double lip_start = kOralLength - lip_length;
  double lip_area = area_from_root(lip_root);
  // Incisors: about 1 cm behind the lip opening, 0.5 cm deep.
  double teeth_gap = std::max(shape.teeth_gap, 0.02);
  // Jaw closure: a raised jaw lifts the tongue floor and narrows the whole
  // front of the mouth (from about 8 cm to the lips), lowering F1 toward
  // every consonant. Vowel postures are fitted with the jaw open (0).
  double jaw_closure = std::clamp(shape.jaw_closure, 0.0, 1.0);
  double teeth_start = lip_start - 0.5;
  double teeth_end = lip_start;
  double constriction_radius = std::sqrt(tongue_area / kPi);
  double curvature = 0.025 * std::max(1.2 - constriction_radius, 0.1);
  // The tip is the front of the tongue: once it is raised to a contact or
  // a channel, nothing of the body lies ahead of it, and the space between
  // the tip and the teeth (with the sublingual space below the blade) is
  // the front cavity that tunes a sibilant. Drawn as one parabola, the body
  // of the /i/ in "She" filled that space to 0.26-0.41 cm^2, and /sh/ came
  // out flat from 1.8 to 5 kHz (owner: "a spoon in the mouth"; the owner's
  // own /sh/ peaks at 2.7 kHz). A released tip (root 2.83) frees nothing.
  double tip_engaged = std::clamp((1.0 - shape.tip_root * std::abs(shape.tip_root)) / 0.7, 0.0, 1.0);
  double free_front = std::clamp(shape.free_front, 0.0, 1.0) * tip_engaged;

  for (Index index = 0; index < kOralSections; ++index) {
    double x = (index + 0.5) * kSectionLength;
    double distance = x - place;
    double radius = constriction_radius + curvature * distance * distance;
    double ceiling = 1.55;
    if (x > lip_start - 2.0) {
      ceiling = 1.4;
    }
    radius = std::min(radius, ceiling);
    double area = kPi * radius * radius;
    double ahead = std::clamp((x - tip_place) / kSectionLength, 0.0, 1.0);
    area += free_front * ahead * (kPi * ceiling * ceiling - area);
    double jaw_reach = std::clamp((x - 8.0) / 3.0, 0.0, 1.0);
    area *= 1.0 - kJawNarrowing * jaw_closure * jaw_reach;
    // Epilarynx: the narrow tube above the folds.
    if (x < kEpilarynxTop) {
      double blend = std::max(x - (kEpilarynxTop - 1.0), 0.0) / 1.0;
      area = kEpilarynxArea + blend * blend * (area - kEpilarynxArea);
    }
    // Lips: a boundary that covers part of a section blends it, so moving
    // the lips never snaps a whole section between shapes.
    double lip_cover = std::clamp((x - lip_start) / kSectionLength + 0.5, 0.0, 1.0);
    double into_lips = std::clamp((x - lip_start) / 0.5, 0.0, 1.0);
    area = area + lip_cover * into_lips * (lip_area - area);
    // Tip and dorsum constrictions: signed root areas; pressing (negative)
    // flattens tissue over a wider region. Each section takes the narrowest
    // point within its own span, not at its centre: sampled at centres, a
    // tip aimed between two sections opened by up to 0.13 cm^2, so a tip
    // sliding from /s/ to /t/ reopened mid-way and the jet burst (a click).
    double tip_distance = std::max(std::abs(x - tip_place) - 0.5 * kSectionLength, 0.0);
    double tip_signed = tip_root * std::abs(tip_root);
    double tip_area = tip_signed + 4.0 * tip_distance * tip_distance;
    double dorsum_distance = std::max(std::abs(x - dorsum_place) - 0.5 * kSectionLength, 0.0);
    double dorsum_signed = dorsum_root * std::abs(dorsum_root);
    double dorsum_area = dorsum_signed + 1.2 * dorsum_distance * dorsum_distance;
    area = std::min(area, std::max(tip_area, 0.0));
    area = std::min(area, std::max(dorsum_area, 0.0));
    area = area + lip_cover * (std::min(area, lip_area) - area);
    double teeth_cover = std::clamp((x - teeth_start) / kSectionLength + 0.5, 0.0, 1.0) *
                         std::clamp((teeth_end - x) / kSectionLength + 0.5, 0.0, 1.0);
    area = area + teeth_cover * (std::min(area, teeth_gap) - area);
    areas[static_cast<std::size_t>(index)] = std::max(area, 0.0);
  }

  TractLayout layout;
  layout.lip_start = lip_start;
  layout.teeth_start = teeth_start;
  layout.teeth_end = teeth_end;
  return layout;
}

}  // namespace voice
