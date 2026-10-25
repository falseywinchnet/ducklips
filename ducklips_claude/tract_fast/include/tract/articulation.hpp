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

// The tongue's upper surface as points from the root's base (by the
// epiglottis) over the dorsum to the tip, in the jaw's frame: the shape a
// physical tongue (tract/tongue.hpp) holds, in place of the body circle.
constexpr int kTongueContourPoints = 24;
struct TongueContour {
  int count = 0;
  std::array<double, kTongueContourPoints> x{};
  std::array<double, kTongueContourPoints> y{};
};

// A thick tongue's front below its upper surface (jaw frame): the rounded
// tip's cap (three points from its top to its bottom, past the front), the
// underside back from the tip, and the frenulum, where the underside meets
// the floor of the mouth.
constexpr int kTongueUndersidePoints = 10;
struct TongueUnderside {
  int count = 0;
  std::array<double, kTongueUndersidePoints> x{};
  std::array<double, kTongueUndersidePoints> y{};
};

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
  // When set, the tongue is this contour (its last point the tip), not the
  // body circle and tip above.
  const TongueContour* contour = nullptr;
  // With a contour, a thick tongue's cap, underside and frenulum (the space
  // between the underside and the floor is then air, under the tongue).
  const TongueUnderside* underside = nullptr;
  // With a contour: each cut's width is also bounded by the hull's distance
  // to the nearest organ surface (continuous as the tip slides past a cut).
  bool nearest_width = false;
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

// A tongue drawn as VocalTractLab 2.4 draws its midsagittal contour
// (VocalTract::calcTongue, GPL: mathematics only; owner 2026-10-18: "your
// base is basically.. a bad tongue. real tongues dont look like this", then
// chose a VTL-style tongue on our own hull and jaw): the root a rational
// quadratic Bezier (weights 1, 2, 1) from the root's base through a control
// point to its tangent on the body ellipse, the ellipse's arc over the top,
// the blade a Bezier from the ellipse through a control point to its tangent
// on the tip circle, and the tip circle's arc to its front. In the jaw's
// frame: the body ellipse centred on the body's position, the tip circle so
// its front lands on the planned tip. The blade's control point (VTL keeps
// one per shape) is placed `blade_along` of the way from the body centre to
// the tip and `blade_lift` cm above that line: JD2's 80 shapes put it at a
// median 0.61 and 1.15 cm. The root's control point is midway from the base
// to the ellipse and `root_bulge` cm behind that line (VTL fixes it by
// constants of its own frame).
struct VtlTongueShape {
  bool on = false;
  double body_rx = 1.8;      // cm, the body ellipse's half axes (JD2)
  double body_ry = 1.8;
  double tip_radius = 0.2;   // cm (JD2)
  double blade_along = 0.61;
  double blade_lift = 1.15;  // cm
  double root_bulge = 0.3;   // cm
};

class MidsagittalTract {
 public:
  MidsagittalTract();

  // The shape for these organ positions.
  void compute(const Articulators& organs, TractProfile& profile);
  // The geometric tongue (root base, body circle, blade, tip) for these
  // positions as kTongueContourPoints points evenly spaced along it, in the
  // jaw's frame: the shapes a physical tongue aims at and rests in.
  void tongue_outline(const Articulators& organs, TongueContour& out) const;
  // The VocalTractLab-style tongue (VtlTongueShape) for these positions,
  // likewise kTongueContourPoints points in the jaw's frame, root to tip.
  void vtl_tongue_outline(const Articulators& organs, TongueContour& out) const;
  void set_vtl_tongue(const VtlTongueShape& shape) { vtl_tongue_ = shape; }
  const VtlTongueShape& vtl_tongue() const { return vtl_tongue_; }
  // The root's base (skull frame) and the conversions between frames.
  void root_base(const Articulators& organs, double& x, double& y) const;
  void to_jaw_frame(const Articulators& organs, double x, double y, double& out_x, double& out_y) const;
  void from_jaw_frame(const Articulators& organs, double x, double y, double& out_x, double& out_y) const {
    jaw_point(organs, x, y, out_x, out_y);
  }

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
  // ... and its cross-section times this. (Every vowel of ours sat 7-16 dB
  // under bdl at 3.0-3.4 kHz and 3-6 dB over at 0.8-2 kHz, whatever the
  // folds' pulse: the narrower the epilarynx against the pharynx, the
  // stronger its resonance gathers F3 and F4, 2026-10-18.)
  void set_larynx_area_scale(double scale) { larynx_area_scale_ = scale; }
  // An /l/'s two side channels together (cm^2). (0.12, set for F1 near 290
  // Hz; against bdl, rms and slt our /l/ then sat 7-19 dB short at
  // 0.25-1.3 kHz, 2026-10-19.)
  void set_lateral_area(double area) { lateral_area_ = area; }
  // ... running this far behind the tip contact (cm). (At 2.0, three of ten
  // /l/ in "Please let all little lambs loll lazily" sealed 0.001 cm^2
  // further back, the blade on the palate: silent laterals.)
  void set_lateral_reach(double reach) { lateral_reach_ = reach; }

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
  void build_front_organs(const Articulators& organs, double under_x, double under_y);
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
  // The distance from a point to the nearest segment of the organ outline.
  double organ_distance(double x, double y) const;

  std::array<AreaRule, 7> rules_{};
  double scale_ = 1.0;
  double root_follow_ = 0.0;
  VtlTongueShape vtl_tongue_;
  double larynx_height_ = 2.3;
  double larynx_area_scale_ = 1.0;
  double lateral_area_ = 0.12;
  double lateral_reach_ = 2.0;
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
