#pragma once

// The midsagittal articulatory model: the shape of the vocal tract from the
// positions of its organs, in place of a formula for its area function.
//
// A fixed hull (posterior pharyngeal wall, velum, hard palate, alveolar
// ridge, upper incisors) faces the moving organs: the larynx (its height
// moves the glottis), the jaw (a rotation about the condyle that carries
// the lower incisors, the lower lip and the whole tongue), the tongue (one
// contour: a root rising from the epiglottis, a body circle, a blade and a
// tip) and the lips (opening and protrusion). A semipolar grid of cut lines
// (horizontal through the pharynx, a fan through the mouth, vertical at the
// front) crosses the airway; each cut gives the midsagittal width between
// the organs and the hull and the channel's midpoint. Widths become areas
// by region (A = alpha d^beta); the distances between successive midpoints
// are the sections' lengths, so the tract lengthens as the larynx drops or
// the lips protrude. Our own design, after the semipolar grids of
// Mermelstein (1973) and Maeda (1990) and the area rules of Heinz & Stevens
// (1965); fitted to Story's (2008) MRI area functions.
//
// Coordinates in cm: x forward, y up; the glottis at rest is (1.0, 0.0).

#include "tract/units.hpp"

#include <array>

namespace tract {

// Where each organ is. The tongue's points are in the jaw's frame (the jaw
// at rest, teeth together): opening the jaw carries them.
struct Articulators {
  double jaw_angle = 6.0;        // degrees, 0 teeth together .. about 18 wide open
  double tongue_x = 4.0;         // body circle centre (jaw frame)
  double tongue_y = 6.6;
  double tip_x = 6.9;            // tongue tip (jaw frame)
  double tip_y = 7.4;
  double lip_aperture = 1.0;     // cm between the lips (negative: pressed together)
  double lip_protrusion = 0.0;   // cm the lips move forward (0 .. 1.2)
  double lip_width = 4.0;        // cm, the lip opening's width (spread 4.5 .. rounded 2.0)
  double larynx_drop = 0.0;      // cm the larynx is lowered (0 .. 1.5)
  double lateral = 0.0;          // 0..1: side channels past a tip contact (/l/)
  double groove = 0.0;           // 0..1: a channel along the blade's midline (/s/)
};

// The region a cut belongs to (sets its area rule).
enum class CutRegion { larynx, pharynx, velar, palatal, alveolar, teeth, lips };

constexpr int kMaxCuts = 128;

struct TractProfile {
  int count = 0;                              // cuts in use, glottis first
  std::array<double, kMaxCuts> width{};       // midsagittal width, cm (0 closed)
  std::array<double, kMaxCuts> area{};        // cm^2
  std::array<double, kMaxCuts> arc{};         // cm from the glottis to the cut's midpoint
  std::array<double, kMaxCuts> mid_x{};       // the channel's midpoint
  std::array<double, kMaxCuts> mid_y{};
  std::array<CutRegion, kMaxCuts> region{};
  double length = 0.0;                        // cm, glottis to the lip opening
  double teeth_arc = 0.0;                     // cm from the glottis to the incisors' edge
  double lip_arc = 0.0;                       // where the lips begin
  double velum_arc = 0.0;                     // the nasal port
  double tip_arc = 0.0;                       // nearest cut to the tongue tip
};

// Area rule A = alpha d^beta for one region.
struct AreaRule {
  double alpha = 1.5;
  double beta = 1.5;
};

class MidsagittalTract {
 public:
  MidsagittalTract();

  // The shape for these organ positions.
  void compute(const Articulators& organs, TractProfile& profile);

  void set_area_rule(CutRegion region, const AreaRule& rule) { rules_[static_cast<std::size_t>(region)] = rule; }
  const AreaRule& area_rule(CutRegion region) const { return rules_[static_cast<std::size_t>(region)]; }
  // The speaker's size: every length of the anatomy times this.
  void set_scale(double scale) { scale_ = scale; }
  // How far the tongue root's base follows the body forward and back (the
  // genioglossus pulls the root forward with the body: /i/'s wide pharynx).
  void set_root_follow(double share) { root_follow_ = share; }
  // The epilarynx (the narrow tube above the folds) runs this far up (cm):
  // its resonance gathers F3 and F4 into the speaker's cluster.
  void set_larynx_height(double height) { larynx_height_ = height; }

  // The outlines of the last computed shape, for drawing: hull and organs.
  int hull_points() const { return hull_count_; }
  double hull_x(int index) const { return hull_x_[static_cast<std::size_t>(index)]; }
  double hull_y(int index) const { return hull_y_[static_cast<std::size_t>(index)]; }
  int organ_points() const { return organ_count_; }
  double organ_x(int index) const { return organ_x_[static_cast<std::size_t>(index)]; }
  double organ_y(int index) const { return organ_y_[static_cast<std::size_t>(index)]; }

 private:
  static constexpr int kMaxOutline = 160;

  void build_hull(const Articulators& organs);
  void build_organs(const Articulators& organs);
  void add_hull(double x, double y);
  void add_organ(double x, double y);
  void jaw_point(const Articulators& organs, double x, double y, double& out_x, double& out_y) const;
  // Nearest crossing of the ray (origin, direction) with an outline at a
  // ray parameter above `from`; false when it misses.
  bool cross(const std::array<double, kMaxOutline>& xs, const std::array<double, kMaxOutline>& ys, int count,
             bool closed, double origin_x, double origin_y, double direction_x, double direction_y, double from,
             bool farthest_below, double below, double& hit) const;
  void add_cut(const Articulators& organs, double origin_x, double origin_y, double direction_x, double direction_y,
               CutRegion region, TractProfile& profile);
  double area_of(const Articulators& organs, CutRegion region, double width, double mid_x) const;
  bool inside_organs(double x, double y) const;

  std::array<AreaRule, 7> rules_{};
  double scale_ = 1.0;
  double root_follow_ = 0.0;
  double larynx_height_ = 2.3;
  std::array<double, kMaxOutline> hull_x_{};
  std::array<double, kMaxOutline> hull_y_{};
  int hull_count_ = 0;
  std::array<double, kMaxOutline> organ_x_{};
  std::array<double, kMaxOutline> organ_y_{};
  int organ_count_ = 0;
  double glottis_y_ = 0.0;
  double lip_front_ = 9.0;       // x of the lip opening
  double lip_back_ = 8.1;        // x where the lips begin
  double tip_abs_x_ = 0.0;
  double tip_abs_y_ = 0.0;
};

}  // namespace tract
