#include "tract/articulation.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

// Anatomy (cm), an adult male, midsagittal. The glottis at rest is at
// (1.0, 0.0); the hull is fixed to the skull.
constexpr double kGlottisX = 1.0;
constexpr double kLarynxBack = 0.75;     // posterior wall of the larynx tube
constexpr double kLarynxFront = 1.25;    // its anterior wall (epilarynx 0.5 cm deep)
// (The epilarynx's length is MidsagittalTract::larynx_height_.)
constexpr double kJawPivotX = -0.5;      // condyle
constexpr double kJawPivotY = 9.5;
constexpr double kTongueRadius = 2.0;    // body circle
constexpr double kLipAxisY = 7.85;       // height of the lip opening's centre
constexpr double kLipBack = 8.1;         // where the lips begin, unprotruded
constexpr double kLipLength = 1.0;       // front to back, unprotruded
// The grid.
constexpr double kStep = 0.3;            // cm between parallel cuts
constexpr double kFanX = 3.4;            // centre of the fan through the mouth
constexpr double kFanY = 5.6;
constexpr double kFanFirst = 180.0;      // degrees: pointing back, continuing the pharynx cuts
constexpr double kFanLast = 38.0;        // pointing forward at the incisors
constexpr double kFanStep = 5.5;
constexpr double kFanFineBelow = 70.0;   // finer steps over the ridge, where the tip makes contact
constexpr double kFanFineStep = 2.5;
constexpr double kFrontFirstX = 7.35;    // the first vertical cut (the teeth onward)
constexpr double kTeethZoneX = 7.45;     // from here to the lips the passage is between the incisors

// Fixed hull from the top of the trachea's back wall forward to the incisors
// (the lips are added per shape): x, y.
constexpr int kHullFixed = 18;
constexpr double kHullFixedX[kHullFixed] = {0.0, 0.0, 0.15, 0.6, 1.4, 2.2, 2.9, 3.6, 4.4,
                                            5.2, 5.9, 6.5, 6.9, 7.3, 7.6, 7.85, 7.9, 8.0};
constexpr double kHullFixedY[kHullFixed] = {3.0, 8.0, 8.7, 9.2, 9.55, 9.85, 10.1, 10.35, 10.5,
                                            10.45, 10.2, 9.8, 9.3, 8.8, 8.3, 8.0, 8.1, 8.25};
constexpr double kVelumPortX = 0.5;      // the nasal port opens here, behind the velum
constexpr double kVelumPortY = 8.7;
constexpr double kDegree = kPi / 180.0;

double distance(double ax, double ay, double bx, double by) {
  return std::sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
}

}  // namespace

MidsagittalTract::MidsagittalTract() {
  // Fitted to Story's (2008) eleven MRI vowels together (tools: analysis/
  // geofit.py; rms log-area error 0.39 on the grid with the fine fan over
  // the ridge): area rules per region, the speaker's size and how far
  // the root follows the body.
  rules_[static_cast<std::size_t>(CutRegion::larynx)] = AreaRule{1.4313, 1.3997};
  rules_[static_cast<std::size_t>(CutRegion::pharynx)] = AreaRule{0.9890, 1.3612};
  rules_[static_cast<std::size_t>(CutRegion::velar)] = AreaRule{0.8807, 1.5268};
  rules_[static_cast<std::size_t>(CutRegion::palatal)] = AreaRule{0.6459, 2.0071};
  rules_[static_cast<std::size_t>(CutRegion::alveolar)] = AreaRule{0.8904, 1.7104};
  rules_[static_cast<std::size_t>(CutRegion::teeth)] = AreaRule{2.1730, 1.0033};
  rules_[static_cast<std::size_t>(CutRegion::lips)] = AreaRule{1.0, 1.0};
  // Story's speaker fitted 0.953; bdl's formants fit a smaller tract (mean
  // distance over six vowels 3.49 / 3.09 / 2.87 st at 0.953 / 0.92 / 0.89).
  scale_ = 0.90;
  root_follow_ = 0.2215;
}

void MidsagittalTract::add_hull(double x, double y) {
  if (hull_count_ < kMaxOutline) {
    hull_x_[static_cast<std::size_t>(hull_count_)] = x;
    hull_y_[static_cast<std::size_t>(hull_count_)] = y;
    ++hull_count_;
  }
}

void MidsagittalTract::add_organ(double x, double y) {
  if (organ_count_ < kMaxOutline) {
    organ_x_[static_cast<std::size_t>(organ_count_)] = x;
    organ_y_[static_cast<std::size_t>(organ_count_)] = y;
    ++organ_count_;
  }
}

// A point given in the jaw's frame (jaw at rest), carried by its rotation.
// Opening the jaw turns points ahead of the condyle down and back.
void MidsagittalTract::jaw_point(const Articulators& organs, double x, double y, double& out_x,
                                 double& out_y) const {
  double angle = -organs.jaw_angle * kDegree;
  double dx = x - kJawPivotX;
  double dy = y - kJawPivotY;
  out_x = kJawPivotX + dx * std::cos(angle) - dy * std::sin(angle);
  out_y = kJawPivotY + dx * std::sin(angle) + dy * std::cos(angle);
}

void MidsagittalTract::to_jaw_frame(const Articulators& organs, double x, double y, double& out_x,
                                    double& out_y) const {
  double angle = organs.jaw_angle * kDegree;
  double dx = x - kJawPivotX;
  double dy = y - kJawPivotY;
  out_x = kJawPivotX + dx * std::cos(angle) - dy * std::sin(angle);
  out_y = kJawPivotY + dx * std::sin(angle) + dy * std::cos(angle);
}

// The tongue root's base, beside the epiglottis: it follows the hyoid part
// of the way down, and (root_follow) the body forward and back.
void MidsagittalTract::root_base(const Articulators& organs, double& x, double& y) const {
  double drop = std::clamp(organs.larynx_drop, 0.0, 2.0);
  x = 1.35 + root_follow_ * (organs.tongue_x - 4.0);
  y = 2.6 - 0.7 * drop;
}

void MidsagittalTract::tongue_outline(const Articulators& organs, TongueContour& out) const {
  // The outline in the jaw's frame: the root's base (converted), the
  // tangent to the back of the body circle, the arc over its top, the
  // tangent to the tip; then even spacing along it.
  constexpr int kDense = 120;
  std::array<double, kDense + 4> xs{};
  std::array<double, kDense + 4> ys{};
  int dense = 0;
  double root_x = 0.0;
  double root_y = 0.0;
  root_base(organs, root_x, root_y);
  double jaw_root_x = 0.0;
  double jaw_root_y = 0.0;
  to_jaw_frame(organs, root_x, root_y, jaw_root_x, jaw_root_y);
  xs[static_cast<std::size_t>(dense)] = jaw_root_x;
  ys[static_cast<std::size_t>(dense)] = jaw_root_y;
  ++dense;
  double centre_x = organs.tongue_x;
  double centre_y = organs.tongue_y;
  double to_x = centre_x - jaw_root_x;
  double to_y = centre_y - jaw_root_y;
  double reach = std::sqrt(to_x * to_x + to_y * to_y);
  double back_angle = std::atan2(to_y, to_x) + kPi;
  if (reach > kTongueRadius) {
    back_angle = std::atan2(-to_y, -to_x) - std::acos(kTongueRadius / reach);
  }
  double from_x = organs.tip_x - centre_x;
  double from_y = organs.tip_y - centre_y;
  double tip_reach = std::sqrt(from_x * from_x + from_y * from_y);
  double front_angle = std::atan2(from_y, from_x);
  if (tip_reach > kTongueRadius) {
    front_angle += std::acos(kTongueRadius / tip_reach);
  }
  double start = back_angle;
  while (start < front_angle) {
    start += 2.0 * kPi;
  }
  constexpr int kArc = kDense - 2;
  for (int point = 0; point <= kArc; ++point) {
    double angle = start + (front_angle - start) * static_cast<double>(point) / kArc;
    xs[static_cast<std::size_t>(dense)] = centre_x + kTongueRadius * std::cos(angle);
    ys[static_cast<std::size_t>(dense)] = centre_y + kTongueRadius * std::sin(angle);
    ++dense;
  }
  xs[static_cast<std::size_t>(dense)] = organs.tip_x;
  ys[static_cast<std::size_t>(dense)] = organs.tip_y;
  ++dense;
  // Even spacing by arc length.
  std::array<double, kDense + 4> along{};
  for (int point = 1; point < dense; ++point) {
    std::size_t i = static_cast<std::size_t>(point);
    along[i] = along[i - 1] + distance(xs[i - 1], ys[i - 1], xs[i], ys[i]);
  }
  double total = along[static_cast<std::size_t>(dense - 1)];
  out.count = kTongueContourPoints;
  int segment = 1;
  for (int point = 0; point < kTongueContourPoints; ++point) {
    double wanted = total * static_cast<double>(point) / static_cast<double>(kTongueContourPoints - 1);
    while (segment < dense - 1 && along[static_cast<std::size_t>(segment)] < wanted) {
      ++segment;
    }
    std::size_t b = static_cast<std::size_t>(segment);
    double span = std::max(along[b] - along[b - 1], 1.0e-9);
    double share = std::clamp((wanted - along[b - 1]) / span, 0.0, 1.0);
    out.x[static_cast<std::size_t>(point)] = xs[b - 1] + share * (xs[b] - xs[b - 1]);
    out.y[static_cast<std::size_t>(point)] = ys[b - 1] + share * (ys[b] - ys[b - 1]);
  }
}

namespace {

// The angle on an ellipse (centre cx, cy; half axes a, b) of the point whose
// tangent passes through (hx, hy); `clockwise` picks which of the two
// (VocalTractLab's getEllipseTangent).
double ellipse_tangent(double hx, double hy, double cx, double cy, double a, double b, bool clockwise) {
  double rx = cx - hx;
  double ry = cy - hy;
  double radicand = std::max(rx * rx * b * b - b * b * a * a + ry * ry * a * a, 0.0);
  double root = std::sqrt(radicand);
  double denominator = rx * rx * b * b + ry * ry * a * a;
  if (std::abs(denominator) < 1.0e-6) {
    denominator = 1.0e-6;
  }
  double alpha0 = std::atan2(-b * (ry * a * a + rx * root) / denominator, -a * (rx * b * b - ry * root) / denominator);
  double alpha1 = std::atan2(-b * (ry * a * a - rx * root) / denominator, -a * (rx * b * b + ry * root) / denominator);
  double cosine = std::cos(alpha0);
  double sine = std::sin(alpha0);
  double v1x = cx + a * cosine - hx;
  double v1y = cy + b * sine - hy;
  double test = v1x * (-a * sine) + v1y * (b * cosine);
  if ((test >= 0.0 && clockwise) || (test < 0.0 && !clockwise)) {
    return alpha1;
  }
  return alpha0;
}

// The angle on a circle (centre cx, cy, radius r) of the tangent point from
// (hx, hy) (VocalTractLab's getCircleTangent).
double circle_tangent(double hx, double hy, double cx, double cy, double r, bool clockwise) {
  double p = hx - cx;
  double q = hy - cy;
  double radicand = p * p + q * q - r * r;
  if (radicand < 0.0) {
    return 0.0;
  }
  double n = std::sqrt(radicand);
  if (clockwise) {
    n = -n;
  }
  return std::atan2(p * n + q * r, p * r - q * n);
}

// A rational quadratic Bezier with weights 1, 2, 1 at t.
void bezier_point(double x0, double y0, double x1, double y1, double x2, double y2, double t, double& x, double& y) {
  double b0 = (1.0 - t) * (1.0 - t);
  double b1 = 2.0 * t * (1.0 - t) * 2.0;
  double b2 = t * t;
  double sum = b0 + b1 + b2;
  x = (b0 * x0 + b1 * x1 + b2 * x2) / sum;
  y = (b0 * y0 + b1 * y1 + b2 * y2) / sum;
}

}  // namespace

void MidsagittalTract::vtl_tongue_outline(const Articulators& organs, TongueContour& out) const {
  constexpr int kPieces = 32;
  constexpr int kDense = 3 * kPieces + 9;
  std::array<double, kDense> xs{};
  std::array<double, kDense> ys{};
  int dense = 0;
  const VtlTongueShape& shape = vtl_tongue_;
  double base_x = 0.0;
  double base_y = 0.0;
  double skull_x = 0.0;
  double skull_y = 0.0;
  root_base(organs, skull_x, skull_y);
  to_jaw_frame(organs, skull_x, skull_y, base_x, base_y);
  double cx = organs.tongue_x;
  double cy = organs.tongue_y;
  double a = shape.body_rx;
  double b = shape.body_ry;
  double r = shape.tip_radius;
  double tx = organs.tip_x - r;  // the tip circle's centre: its front on the planned tip
  double ty = organs.tip_y;
  // The root: base -> control -> tangent on the ellipse.
  double back = ellipse_tangent(base_x, base_y, cx, cy, a, b, true);
  double touch_x = cx + a * std::cos(back);
  double touch_y = cy + b * std::sin(back);
  double along_x = touch_x - base_x;
  double along_y = touch_y - base_y;
  double length = std::max(std::sqrt(along_x * along_x + along_y * along_y), 1.0e-6);
  // Behind the base-to-ellipse line: its right-hand normal (the line runs up
  // and forward; behind is toward the pharynx).
  double root_cx = base_x + 0.5 * along_x + shape.root_bulge * along_y / length;
  double root_cy = base_y + 0.5 * along_y - shape.root_bulge * along_x / length;
  double alpha0 = ellipse_tangent(root_cx, root_cy, cx, cy, a, b, true);
  double end_x = cx + a * std::cos(alpha0);
  double end_y = cy + b * std::sin(alpha0);
  for (int i = 0; i < kPieces; ++i) {
    bezier_point(base_x, base_y, root_cx, root_cy, end_x, end_y, static_cast<double>(i) / kPieces,
                 xs[static_cast<std::size_t>(dense)], ys[static_cast<std::size_t>(dense)]);
    ++dense;
  }
  // The blade's control point.
  double to_tip_x = tx - cx;
  double to_tip_y = ty - cy;
  double span = std::max(std::sqrt(to_tip_x * to_tip_x + to_tip_y * to_tip_y), 1.0e-6);
  double blade_x = cx + shape.blade_along * to_tip_x - shape.blade_lift * to_tip_y / span;
  double blade_y = cy + shape.blade_along * to_tip_y + shape.blade_lift * to_tip_x / span;
  double alpha1 = ellipse_tangent(blade_x, blade_y, cx, cy, a, b, false);
  if (alpha1 > alpha0 && alpha1 - alpha0 < 0.25 * kPi) {
    alpha1 = alpha0;
  }
  // The arc over the top, from the root's tangent to the blade's (angles falling).
  double first = alpha0 < 0.0 ? alpha0 + 2.0 * kPi : alpha0;
  double delta = alpha1 - first;
  if (delta > -1.0e-6) {
    delta = -1.0e-6;
  }
  for (int i = 0; i < kPieces; ++i) {
    double angle = first + delta * static_cast<double>(i) / kPieces;
    xs[static_cast<std::size_t>(dense)] = cx + a * std::cos(angle);
    ys[static_cast<std::size_t>(dense)] = cy + b * std::sin(angle);
    ++dense;
  }
  // The blade: ellipse -> control -> tangent on the tip circle.
  double start_x = cx + a * std::cos(alpha1);
  double start_y = cy + b * std::sin(alpha1);
  double alpha2 = circle_tangent(blade_x, blade_y, tx, ty, r, true);
  double tip_in_x = tx + r * std::cos(alpha2);
  double tip_in_y = ty + r * std::sin(alpha2);
  for (int i = 0; i < kPieces; ++i) {
    bezier_point(start_x, start_y, blade_x, blade_y, tip_in_x, tip_in_y, static_cast<double>(i) / kPieces,
                 xs[static_cast<std::size_t>(dense)], ys[static_cast<std::size_t>(dense)]);
    ++dense;
  }
  // The tip's arc to its front (angle 0).
  double from = alpha2 < 0.0 ? alpha2 + 2.0 * kPi : alpha2;
  double to = 0.0 - from;
  if (to > -1.0e-6) {
    to = -1.0e-6;
  }
  for (int i = 0; i < 8; ++i) {
    double angle = from + to * static_cast<double>(i) / 7.0;
    xs[static_cast<std::size_t>(dense)] = tx + r * std::cos(angle);
    ys[static_cast<std::size_t>(dense)] = ty + r * std::sin(angle);
    ++dense;
  }
  // Even spacing by arc length.
  std::array<double, kDense> along{};
  for (int point = 1; point < dense; ++point) {
    std::size_t i = static_cast<std::size_t>(point);
    along[i] = along[i - 1] + distance(xs[i - 1], ys[i - 1], xs[i], ys[i]);
  }
  double total = along[static_cast<std::size_t>(dense - 1)];
  out.count = kTongueContourPoints;
  int segment = 1;
  for (int point = 0; point < kTongueContourPoints; ++point) {
    double wanted = total * static_cast<double>(point) / static_cast<double>(kTongueContourPoints - 1);
    while (segment < dense - 1 && along[static_cast<std::size_t>(segment)] < wanted) {
      ++segment;
    }
    std::size_t s1 = static_cast<std::size_t>(segment);
    double piece = std::max(along[s1] - along[s1 - 1], 1.0e-9);
    double share = std::clamp((wanted - along[s1 - 1]) / piece, 0.0, 1.0);
    out.x[static_cast<std::size_t>(point)] = xs[s1 - 1] + share * (xs[s1] - xs[s1 - 1]);
    out.y[static_cast<std::size_t>(point)] = ys[s1 - 1] + share * (ys[s1] - ys[s1 - 1]);
  }
}

// The hull: the larynx's back wall (moving with the larynx), the posterior
// pharyngeal wall, velum, palate, ridge, upper incisors, then the upper lip.
void MidsagittalTract::build_hull(const Articulators& organs) {
  hull_count_ = 0;
  double drop = std::clamp(organs.larynx_drop, 0.0, 2.0);
  glottis_y_ = -drop;
  add_hull(kLarynxBack, glottis_y_ - 0.2);
  add_hull(kLarynxBack, glottis_y_ + larynx_height_);
  add_hull(0.0, std::max(glottis_y_ + larynx_height_ + 0.7, 2.0));
  // (The fixed wall's first point is replaced by the larynx-dependent one above.)
  for (int point = 1; point < kHullFixed; ++point) {
    add_hull(kHullFixedX[point], kHullFixedY[point]);
  }
  double protrusion = std::clamp(organs.lip_protrusion, 0.0, 1.5);
  double opening = std::max(organs.lip_aperture, 0.0);
  lip_back_ = kLipBack + 0.3 * protrusion;
  lip_front_ = kLipBack + kLipLength + protrusion;
  add_hull(lip_back_, kLipAxisY + 0.5 * opening + 0.15);
  add_hull(lip_back_ + 0.2, kLipAxisY + 0.5 * opening);
  add_hull(lip_front_, kLipAxisY + 0.5 * opening);
  add_hull(lip_front_ + 0.3, kLipAxisY + 0.5 * opening + 0.6);
}

// The organs as one closed outline: the larynx's front wall, the tongue
// root, the body circle, the blade, the tip, the floor behind the lower
// incisors, the lower incisors, the lower lip, and back under the chin.
void MidsagittalTract::build_organs(const Articulators& organs) {
  organ_count_ = 0;
  double centre_x = 0.0;
  double centre_y = 0.0;
  jaw_point(organs, organs.tongue_x, organs.tongue_y, centre_x, centre_y);
  double tip_x = 0.0;
  double tip_y = 0.0;
  const TongueContour* contour = organs.contour;
  if (contour != nullptr && contour->count >= 2 && organs.underside != nullptr && organs.underside->count > 1) {
    // A thick tongue's tip is its cap's front (the groove and the side
    // channels are placed from it; from the cap's top, the cap sealed where
    // the /s/ groove should run).
    jaw_point(organs, organs.underside->x[1], organs.underside->y[1], tip_x, tip_y);
  } else if (contour != nullptr && contour->count >= 2) {
    std::size_t last = static_cast<std::size_t>(contour->count - 1);
    jaw_point(organs, contour->x[last], contour->y[last], tip_x, tip_y);
  } else {
    jaw_point(organs, organs.tip_x, organs.tip_y, tip_x, tip_y);
  }
  tip_abs_x_ = tip_x;
  tip_abs_y_ = tip_y;
  // Larynx front wall, then the epiglottis base, which follows the hyoid
  // part of the way down.
  add_organ(kLarynxFront, glottis_y_ - 0.2);
  add_organ(kLarynxFront, glottis_y_ + larynx_height_);
  double root_x = 0.0;
  double root_y = 0.0;
  root_base(organs, root_x, root_y);
  if (contour != nullptr && contour->count >= 2) {
    // A physical tongue: its surface from the root's base to the tip.
    add_organ(root_x, root_y);
    for (int point = 1; point < contour->count; ++point) {
      double x = 0.0;
      double y = 0.0;
      jaw_point(organs, contour->x[static_cast<std::size_t>(point)], contour->y[static_cast<std::size_t>(point)], x, y);
      add_organ(x, y);
    }
    std::size_t last = static_cast<std::size_t>(contour->count - 1);
    double under_x = 0.0;
    double under_y = 0.0;
    if (organs.underside != nullptr && organs.underside->count > 0) {
      // A thick tongue: the cap, the underside back to the frenulum.
      for (int point = 0; point < organs.underside->count; ++point) {
        std::size_t p = static_cast<std::size_t>(point);
        jaw_point(organs, organs.underside->x[p], organs.underside->y[p], under_x, under_y);
        add_organ(under_x, under_y);
      }
      // The floor of the mouth forward from the frenulum (jaw frame), to
      // the floor behind the lower incisors.
      std::size_t frenulum = static_cast<std::size_t>(organs.underside->count - 1);
      constexpr int kFloorPoints = 3;
      constexpr double kFloorX[kFloorPoints] = {4.5, 6.0, 6.5};
      constexpr double kFloorY[kFloorPoints] = {5.0, 6.0, 6.3};
      for (int point = 0; point < kFloorPoints; ++point) {
        if (kFloorX[point] > organs.underside->x[frenulum] + 0.05) {
          jaw_point(organs, kFloorX[point], kFloorY[point], under_x, under_y);
          add_organ(under_x, under_y);
        }
      }
      build_front_organs(organs, under_x, under_y);
      return;
    }
    jaw_point(organs, contour->x[last] - 0.2, contour->y[last] - 0.6, under_x, under_y);
    add_organ(under_x, under_y);
    build_front_organs(organs, under_x, under_y);
    return;
  }
  add_organ(root_x, root_y);
  // Tangent from the root base to the back of the body circle.
  double to_x = centre_x - root_x;
  double to_y = centre_y - root_y;
  double reach = std::sqrt(to_x * to_x + to_y * to_y);
  double back_angle = std::atan2(to_y, to_x) + kPi;  // fallback: straight behind
  if (reach > kTongueRadius) {
    double spread = std::acos(kTongueRadius / reach);
    back_angle = std::atan2(-to_y, -to_x) - spread;  // the tangent on the back (posterior) side
  }
  // Tangent from the tip to the front of the circle.
  double from_x = tip_x - centre_x;
  double from_y = tip_y - centre_y;
  double tip_reach = std::sqrt(from_x * from_x + from_y * from_y);
  double front_angle = std::atan2(from_y, from_x);
  if (tip_reach > kTongueRadius) {
    front_angle += std::acos(kTongueRadius / tip_reach);
  }
  // Arc over the top from the back tangent to the front tangent (clockwise
  // in x-forward, y-up: angles decrease).
  double start = back_angle;
  double end = front_angle;
  while (start < end) {
    start += 2.0 * kPi;
  }
  constexpr int kArcPoints = 40;
  for (int point = 0; point <= kArcPoints; ++point) {
    double angle = start + (end - start) * static_cast<double>(point) / kArcPoints;
    add_organ(centre_x + kTongueRadius * std::cos(angle), centre_y + kTongueRadius * std::sin(angle));
  }
  // Blade to the tip, then the tip's underside down to the floor.
  add_organ(tip_x, tip_y);
  double under_x = 0.0;
  double under_y = 0.0;
  jaw_point(organs, organs.tip_x - 0.2, organs.tip_y - 0.6, under_x, under_y);
  add_organ(under_x, under_y);
  build_front_organs(organs, under_x, under_y);
}

// The organs' outline past the tongue: the floor behind the lower incisors,
// the incisors, the lower lip, the chin and the floor back to the larynx.
void MidsagittalTract::build_front_organs(const Articulators& organs, double under_x, double under_y) {
  double drop = std::clamp(organs.larynx_drop, 0.0, 2.0);
  // Floor behind the lower incisors, the incisors, the lower lip.
  double floor_x = 0.0;
  double floor_y = 0.0;
  jaw_point(organs, 7.1, 6.6, floor_x, floor_y);
  add_organ(std::max(floor_x, under_x + 0.05), std::min(floor_y, under_y));
  double lingual_x = 0.0;
  double lingual_y = 0.0;
  jaw_point(organs, 7.5, 7.3, lingual_x, lingual_y);
  add_organ(lingual_x, lingual_y);
  double edge_x = 0.0;
  double edge_y = 0.0;
  jaw_point(organs, 7.7, 7.75, edge_x, edge_y);
  add_organ(edge_x, edge_y);
  double opening = organs.lip_aperture;
  // A negative aperture presses the lips together: the lower lip's edge
  // rises past the upper's. (It once fell below it, and every /p b m/ was a
  // 0.1 cm slit, about 0.25 cm^2.)
  double lower_lip = kLipAxisY - 0.5 * std::max(opening, 0.0) - std::min(opening, 0.0);
  add_organ(lip_back_, std::min(lower_lip - 0.15, edge_y - 0.05));
  add_organ(lip_back_ + 0.2, lower_lip);
  add_organ(lip_front_, lower_lip);
  add_organ(lip_front_ + 0.3, lower_lip - 0.6);
  // Chin and the floor of the mouth back to the larynx (far from the airway).
  add_organ(lip_front_, 3.0);
  add_organ(3.0, 1.5 - drop);
  add_organ(kLarynxFront + 0.2, glottis_y_ - 0.3);
}

bool MidsagittalTract::cross(const std::array<double, kMaxOutline>& xs, const std::array<double, kMaxOutline>& ys,
                             int count, bool closed, double origin_x, double origin_y, double direction_x,
                             double direction_y, double from, bool farthest_below, double below,
                             double& hit) const {
  bool found = false;
  int segments = closed ? count : count - 1;
  for (int segment = 0; segment < segments; ++segment) {
    std::size_t a = static_cast<std::size_t>(segment);
    std::size_t b = static_cast<std::size_t>((segment + 1) % count);
    double ex = xs[b] - xs[a];
    double ey = ys[b] - ys[a];
    double denominator = direction_x * ey - direction_y * ex;
    if (std::abs(denominator) < 1.0e-12) {
      continue;
    }
    double wx = xs[a] - origin_x;
    double wy = ys[a] - origin_y;
    double t = (wx * ey - wy * ex) / denominator;
    double u = (wx * direction_y - wy * direction_x) / denominator;
    if (u < 0.0 || u > 1.0) {
      continue;
    }
    if (farthest_below) {
      if (t < below && (!found || t > hit)) {
        hit = t;
        found = true;
      }
    } else if (t > from && (!found || t < hit)) {
      hit = t;
      found = true;
    }
  }
  return found;
}

// Ray casting: is the point inside the closed organ outline?
bool MidsagittalTract::inside_organs(double x, double y) const {
  bool inside = false;
  for (int point = 0; point < organ_count_; ++point) {
    std::size_t a = static_cast<std::size_t>(point);
    std::size_t b = static_cast<std::size_t>((point + 1) % organ_count_);
    double ay = organ_y_[a];
    double by = organ_y_[b];
    if ((ay > y) != (by > y)) {
      double crossing = organ_x_[a] + (y - ay) * (organ_x_[b] - organ_x_[a]) / (by - ay);
      if (x < crossing) {
        inside = !inside;
      }
    }
  }
  return inside;
}

double MidsagittalTract::organ_distance(double x, double y) const {
  double nearest = 1.0e9;
  for (int point = 0; point < organ_count_; ++point) {
    std::size_t a = static_cast<std::size_t>(point);
    std::size_t b = static_cast<std::size_t>((point + 1) % organ_count_);
    double dx = organ_x_[b] - organ_x_[a];
    double dy = organ_y_[b] - organ_y_[a];
    double span = dx * dx + dy * dy;
    double t = span > 0.0 ? std::clamp(((x - organ_x_[a]) * dx + (y - organ_y_[a]) * dy) / span, 0.0, 1.0) : 0.0;
    nearest = std::min(nearest, distance(x, y, organ_x_[a] + t * dx, organ_y_[a] + t * dy));
  }
  return nearest;
}

double MidsagittalTract::area_of(const Articulators& organs, CutRegion region, double width, double mid_x) const {
  if (region == CutRegion::lips) {
    // An elliptical opening: width across, the aperture between the lips.
    return width > 0.0 ? 0.25 * kPi * std::max(organs.lip_width, 0.5) * width : 0.0;
  }
  const AreaRule& rule = rules_[static_cast<std::size_t>(region)];
  // (A closed midline still passes the side channels and groove below.)
  double area = width > 0.0 ? rule.alpha * std::pow(width, rule.beta) : 0.0;
  // The tongue root never seals the pharynx in English: it stays at least
  // 0.15 cm^2 open. (A body drawn back for /l/ before /ow/ closed it and the
  // voice stopped for 100 ms in "hello".)
  if (region == CutRegion::pharynx) {
    area = std::max(area, 0.15);
  }
  if (region == CutRegion::larynx && larynx_area_scale_ != 1.0) {
    area *= larynx_area_scale_;
  }
  // Side channels past a tip contact (a lateral) keep a passage open where
  // the midline closes at the tip. A sibilant's groove runs along the whole
  // blade, from about 2 cm behind the tip to the tip: with the jaw raised
  // and the body high (/s/ before /i/), the blade lay on the palate behind
  // the tip and closed the channel there, and "speak" lost its /s/.
  // The sides of the tongue are lowered along the whole blade, so the side
  // channels run from about 2 cm behind the tip contact. (Within 0.8 cm of
  // the tip only, the blade sealed 0.81 cm back, 0.006 cm^2, the flow
  // through the mouth stopped for 40 ms: "hello" heard as "hetlo".)
  bool palatal_or_front = region == CutRegion::alveolar || region == CutRegion::palatal || region == CutRegion::teeth;
  if (palatal_or_front && mid_x > tip_abs_x_ - lateral_reach_ && mid_x < tip_abs_x_ + 0.8) {
    // Both channels together 0.12 cm^2. (At 0.35, with the window now
    // reaching back along the blade, the /l/ was a vowel: F1 525 Hz, no
    // dip; at 0.175 still F1 526; at 0.12, F1 about 290, 9 dB down.)
    area = std::max(area, lateral_area_ * std::clamp(organs.lateral, 0.0, 1.0));
  }
  if (palatal_or_front && mid_x > tip_abs_x_ - 2.2 && mid_x < tip_abs_x_ + 0.5) {
    // Narrowest at the tip (0.06 cm^2), widening to 0.5 two centimetres
    // back: the hiss forms in a short channel just behind the teeth. (A
    // uniform groove made a 2 cm channel and the /s/ peaked at 2.5 kHz, a
    // /sh/; bdl's peaks near 4.5.)
    // (0.12 cm^2 at the tip let a voiced /z/'s 120 cm^3/s through at
    // 10 m/s, 28 Pa behind it: no hiss. 0.06 doubles the jet's speed.)
    double behind = std::clamp((tip_abs_x_ - mid_x) / 2.2, 0.0, 1.0);
    area = std::max(area, (0.06 + 0.34 * behind) * std::clamp(organs.groove, 0.0, 1.0));
  }
  return area;
}

// One cut: the hull crossing nearest along the ray, and the organ crossing
// nearest below it (the organs' surface facing the airway).
void MidsagittalTract::add_cut(const Articulators& organs, double origin_x, double origin_y, double direction_x,
                               double direction_y, CutRegion region, TractProfile& profile) {
  if (profile.count >= kMaxCuts) {
    return;
  }
  double outer = 0.0;
  if (!cross(hull_x_, hull_y_, hull_count_, false, origin_x, origin_y, direction_x, direction_y, 0.0, false, 0.0,
             outer)) {
    return;
  }
  double inner = 0.0;
  if (!cross(organ_x_, organ_y_, organ_count_, true, origin_x, origin_y, direction_x, direction_y, 0.0, true,
             outer + 1.0e-9, inner)) {
    inner = outer;
  }
  // Where an organ is pressed past the hull, the hull's point lies inside
  // the organ: the channel is closed there. (Taking the nearest organ
  // crossing below the hull instead found the tongue's underside and made an
  // 8 cm "gap" through the tongue at every palatal contact.)
  if (inside_organs(origin_x + direction_x * outer, origin_y + direction_y * outer)) {
    inner = outer;
  }
  double width = std::max(outer - inner, 0.0);
  double middle = width > 0.0 ? 0.5 * (inner + outer) : outer;
  if (organs.contour != nullptr && organs.nearest_width && region != CutRegion::lips && width > 0.0) {
    // With the physical tongue the passage is as wide as the hull's distance
    // to the nearest organ surface, which moves continuously with the
    // tongue. (Along the fixed ray alone, a sharp tip sliding past the ray
    // switched its crossing from the tongue's top to its underside: in /dh/
    // the narrowest area jumped 0.005-0.15 cm^2 from one millisecond to the
    // next while the tip moved smoothly.)
    width = std::min(width, organ_distance(origin_x + direction_x * outer, origin_y + direction_y * outer));
    middle = outer - 0.5 * width;
  }
  std::size_t index = static_cast<std::size_t>(profile.count);
  profile.width[index] = width;
  profile.mid_x[index] = origin_x + direction_x * middle;
  profile.mid_y[index] = origin_y + direction_y * middle;
  profile.region[index] = region;
  profile.area[index] = area_of(organs, region, width, profile.mid_x[index]);
  ++profile.count;
}

void MidsagittalTract::compute(const Articulators& organs, TractProfile& profile) {
  build_hull(organs);
  build_organs(organs);
  profile.count = 0;
  // Larynx and pharynx: horizontal cuts looking back at the posterior wall.
  for (double y = glottis_y_ + 0.15; y < kFanY; y += kStep) {
    CutRegion region = y < glottis_y_ + larynx_height_ ? CutRegion::larynx : CutRegion::pharynx;
    add_cut(organs, 9.5, y, -1.0, 0.0, region, profile);
  }
  // The fan through the mouth.
  // (Coarse over the back and the palate, fine over the ridge: with 5.5 deg
  // steps ending at 52 deg, a tip contact on the ridge fell between the
  // fan's last cut and the first vertical one, and /t d n/ never closed.)
  double angle = kFanFirst - 0.5 * kFanStep;
  while (angle > kFanLast) {
    double direction_x = std::cos(angle * kDegree);
    double direction_y = std::sin(angle * kDegree);
    CutRegion region = angle > 135.0  ? CutRegion::pharynx
                       : angle > 95.0 ? CutRegion::velar
                       : angle > 62.0 ? CutRegion::palatal
                                      : CutRegion::alveolar;
    add_cut(organs, kFanX, kFanY, direction_x, direction_y, region, profile);
    angle -= angle > kFanFineBelow ? kFanStep : kFanFineStep;
  }
  // The front: vertical cuts looking up, to the lip opening.
  for (double x = kFrontFirstX; x < lip_front_; x += kStep) {
    // Between the incisors and the lips the passage spans the dental arch,
    // not a midline slit: its own rule.
    CutRegion region = x >= lip_back_ ? CutRegion::lips : x >= kTeethZoneX ? CutRegion::teeth : CutRegion::alveolar;
    add_cut(organs, x, 2.0, 0.0, 1.0, region, profile);
  }
  // Lengths along the midline.
  double last_x = kGlottisX;
  double last_y = glottis_y_;
  double arc = 0.0;
  double best_tip = 1.0e9;
  double best_port = 1.0e9;
  profile.teeth_arc = 0.0;
  profile.lip_arc = 0.0;
  for (int cut = 0; cut < profile.count; ++cut) {
    std::size_t index = static_cast<std::size_t>(cut);
    arc += distance(last_x, last_y, profile.mid_x[index], profile.mid_y[index]);
    profile.arc[index] = arc;
    last_x = profile.mid_x[index];
    last_y = profile.mid_y[index];
    double to_tip = distance(last_x, last_y, tip_abs_x_, tip_abs_y_);
    if (to_tip < best_tip) {
      best_tip = to_tip;
      profile.tip_arc = arc;
    }
    double to_port = distance(last_x, last_y, kVelumPortX, kVelumPortY);
    if (to_port < best_port) {
      best_port = to_port;
      profile.velum_arc = arc;
    }
    if (profile.teeth_arc == 0.0 && last_x >= 7.8) {
      profile.teeth_arc = arc;
    }
    if (profile.lip_arc == 0.0 && profile.region[index] == CutRegion::lips) {
      profile.lip_arc = arc;
    }
  }
  profile.length = arc + distance(last_x, last_y, lip_front_, kLipAxisY);
  // A larger or smaller speaker: lengths and widths scale, areas follow
  // their rules on the scaled widths.
  if (scale_ != 1.0) {
    for (int cut = 0; cut < profile.count; ++cut) {
      std::size_t index = static_cast<std::size_t>(cut);
      profile.arc[index] *= scale_;
      profile.width[index] *= scale_;
      profile.area[index] = profile.region[index] == CutRegion::lips
                                ? profile.area[index] * scale_ * scale_
                                : area_of(organs, profile.region[index], profile.width[index], profile.mid_x[index]);
    }
    profile.length *= scale_;
    profile.teeth_arc *= scale_;
    profile.lip_arc *= scale_;
    profile.velum_arc *= scale_;
    profile.tip_arc *= scale_;
  }
}

}  // namespace tract
