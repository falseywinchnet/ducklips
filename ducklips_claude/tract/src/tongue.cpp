#include "tract/tongue.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

// Where the extrinsic muscles attach (cm): the genial tubercle inside the
// chin (jaw frame), the styloid process behind and above the pharynx and
// the hyoid (skull frame, the hyoid following the larynx down).
constexpr double kGenialX = 6.2;
constexpr double kGenialY = 5.3;
constexpr double kStyloidX = -1.0;
constexpr double kStyloidY = 8.8;
constexpr double kHyoidX = 2.0;
constexpr double kHyoidY = 2.2;
// The floor of the mouth under the tongue (jaw frame), closing the outline
// whose area the hydrostat keeps.
constexpr double kFloorFrontX = 6.2;
constexpr double kFloorFrontY = 5.3;
constexpr double kFloorBackX = 2.8;
constexpr double kFloorBackY = 3.4;
// The upper incisors' edge lies before this (cm): the hull the tongue meets ends there.
constexpr double kIncisorReach = 8.0;

// Each muscle's reach along the surface (0 the root, 1 the tip): a bell
// around `centre`; the origin it pulls toward (0 genial, 1 styloid, 2
// hyoid) or, for the intrinsic muscles, a fixed direction (jaw frame).
struct MuscleLayout {
  double centre;
  double width;
  int origin;            // -1: pulls along (direction_x, direction_y); -2: along the surface's normal;
                         // -3: along the surface toward the tip (elongation)
  double direction_x;    // (origin -2: +1 outward, -1 inward)
  double direction_y;
  double strength;       // times muscle_force
};

// (The longitudinal muscles shorten the upper or the lower surface: the
// blade curls up or down, along its own normal. As a fixed up-and-back pull
// the tip could not be lifted forward onto the ridge: no /t d n/ sealed.)
// The extrinsic muscles are split into their fibre bands, and the intrinsic
// transversus (narrowing: the surface rises there) and verticalis
// (flattening: it falls) act along the normal at five places: with seven
// muscles the best least-squares hold of a target shape left 0.8 cm of
// error, and the feedforward found almost nothing to do.
constexpr std::array<MuscleLayout, Tongue::kMuscles> kMuscleLayout = {{
    {0.15, 0.10, 0, 0.0, 0.0, 1.0},   // genioglossus, root fibres: the root forward
    {0.40, 0.10, 0, 0.0, 0.0, 1.0},   // genioglossus, middle fibres: the dorsum down and forward
    {0.75, 0.10, 0, 0.0, 0.0, 1.0},   // genioglossus, anterior fibres: the blade down
    {0.40, 0.12, 1, 0.0, 0.0, 1.0},   // styloglossus, back: up and back
    {0.62, 0.12, 1, 0.0, 0.0, 1.0},   // styloglossus, front: the dorsum up and back
    {0.35, 0.12, 2, 0.0, 0.0, 1.0},   // hyoglossus, back: down and back
    {0.55, 0.12, 2, 0.0, 0.0, 1.0},   // hyoglossus, front
    {0.92, 0.08, -2, 1.0, 0.0, 3.0},  // superior longitudinal: the blade curled up
    {0.92, 0.08, -2, -1.0, 0.0, 3.0}, // inferior longitudinal: the blade curled down
    {0.55, 0.25, -1, 0.0, 1.0, 1.0},  // the floor's lift (mylohyoid): the body up
    {0.20, 0.07, -2, 1.0, 0.0, 1.0},  // transversus: the surface up here
    {0.35, 0.07, -2, 1.0, 0.0, 1.0},
    {0.50, 0.07, -2, 1.0, 0.0, 1.0},
    {0.65, 0.07, -2, 1.0, 0.0, 1.0},
    {0.80, 0.07, -2, 1.0, 0.0, 1.0},
    {0.20, 0.07, -2, -1.0, 0.0, 1.0}, // verticalis: the surface down here
    {0.35, 0.07, -2, -1.0, 0.0, 1.0},
    {0.50, 0.07, -2, -1.0, 0.0, 1.0},
    {0.65, 0.07, -2, -1.0, 0.0, 1.0},
    {0.80, 0.07, -2, -1.0, 0.0, 1.0},
    // Elongation: squeezed across (transversus and verticalis together), the
    // hydrostat lengthens; the blade and tip are pushed along the surface
    // toward the tip (protrusion). Without it nothing could carry the tip
    // forward: it stopped 0.4 cm short of the incisors for /dh/.
    {0.88, 0.10, -3, 1.0, 0.0, 2.0},
}};

// The tissue grows more pliant toward the tip: the support and the bending
// stiffness of a section, times this (1 at the root, tip_pliancy at the tip).
double pliancy(int section, int sections, double tip_pliancy) {
  double along = static_cast<double>(section) / static_cast<double>(sections - 1);
  return 1.0 + (tip_pliancy - 1.0) * along * along;
}

// The pieces: root (sections up to 5), dorsum (to 14), blade (to 20), tip.
bool piece_physical(int section, const TongueSettings& settings) {
  if (section <= 5) {
    return settings.root_physical > 0.5;
  }
  if (section <= 14) {
    return settings.dorsum_physical > 0.5;
  }
  if (section <= 20) {
    return settings.blade_physical > 0.5;
  }
  return settings.tip_physical > 0.5;
}

double reach_weight(const MuscleLayout& layout, int section, int sections) {
  double along = static_cast<double>(section) / static_cast<double>(sections - 1);
  double offset = (along - layout.centre) / layout.width;
  return std::exp(-0.5 * offset * offset);
}

}  // namespace

void Tongue::reset(const TongueContour& rest, const TongueSettings& settings) {
  settings_ = settings;
  rest_ = rest;
  contour_ = rest;
  target_ = rest;
  vx_.fill(0.0);
  vy_.fill(0.0);
  activation_.fill(0.0);
  integral_.fill(0.0);
  for (int i = 0; i + 1 < kSections; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    rest_spacing_[a] = std::hypot(rest.x[a + 1] - rest.x[a], rest.y[a + 1] - rest.y[a]);
  }
  for (int i = 1; i + 1 < kSections; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    rest_bend_x_[a] = 2.0 * rest.x[a] - rest.x[a - 1] - rest.x[a + 1];
    rest_bend_y_[a] = 2.0 * rest.y[a] - rest.y[a - 1] - rest.y[a + 1];
  }
  rest_area_ = area_of(rest_);
  held_area_ = rest_area_;
  geometric_ = rest;
  anchor_ = rest;
  aim_ = rest;
  guard_ = rest;
  intrusion_.fill(0.0);
  penetration_.fill(0.0);
  guard_weight_.fill(0.0);
  guarding_ = false;
  shaping_.fill(1.0);
  geometric_vx_.fill(0.0);
  geometric_vy_.fill(0.0);
  geometric_set_ = false;
  std::array<double, 4> pieces{};
  pieces.fill(0.0);
  set_pins(pieces);
  impedance_.fill(0.0);
  impedance_damping_.fill(0.0);
  holding_ = false;
  uses_geometric_ = pinning_ || settings_.guard > 0.5 || (settings_.roles > 0.5 && (settings_.role_anchor > 0.5 || settings_.shaper_pinned > 0.0 ||
                                                           settings_.valve_pinned > 0.0 ||
                                                           settings_.nozzle_pinned > 0.0));
  feedforward_.fill(0.0);
  pressure_ = 0.0;
  measured_ = false;
  if (settings_.spine > 0.5) {
    spine_reset();
  }
  if (settings_.skeleton > 0.5) {
    skeleton_reset();
  }
  ready_ = true;
}

void Tongue::place(const TongueContour& shape) {
  contour_ = shape;
  vx_.fill(0.0);
  vy_.fill(0.0);
  if (settings_.spine > 0.5) {
    spine_place(shape);
  }
  if (settings_.skeleton > 0.5) {
    skeleton_place(shape);
  }
}

void Tongue::target_error(double& rms, double& largest) const {
  double sum = 0.0;
  largest = 0.0;
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double gap = std::hypot(target_.x[s] - contour_.x[s], target_.y[s] - contour_.y[s]);
    sum += gap * gap;
    largest = std::max(largest, gap);
  }
  rms = std::sqrt(sum / static_cast<double>(kSections - 1));
}

// An outline's area: the surface, the tip's underside, the floor.
double Tongue::outline_area(const TongueContour& shape) {
  std::array<double, kSections + 3> xs{};
  std::array<double, kSections + 3> ys{};
  for (int i = 0; i < kSections; ++i) {
    xs[static_cast<std::size_t>(i)] = shape.x[static_cast<std::size_t>(i)];
    ys[static_cast<std::size_t>(i)] = shape.y[static_cast<std::size_t>(i)];
  }
  std::size_t tip = static_cast<std::size_t>(kSections - 1);
  xs[tip + 1] = shape.x[tip] - 0.2;
  ys[tip + 1] = shape.y[tip] - 0.6;
  xs[tip + 2] = kFloorFrontX;
  ys[tip + 2] = kFloorFrontY;
  xs[tip + 3] = kFloorBackX;
  ys[tip + 3] = kFloorBackY;
  double twice = 0.0;
  std::size_t count = xs.size();
  for (std::size_t i = 0; i < count; ++i) {
    std::size_t j = (i + 1) % count;
    twice += xs[i] * ys[j] - xs[j] * ys[i];
  }
  return 0.5 * std::abs(twice);
}

// The floor of the mouth and the lower incisors' inner face as a solid
// below the tongue (jaw frame): from the floor's back, forward under the
// frenulum to the incisors' edge, closed far below.
constexpr int kFloorSolidPoints = 9;
constexpr double kFloorSolidX[kFloorSolidPoints] = {2.8, 4.5, 6.0, 6.5, 7.1, 7.5, 7.7, 7.7, 2.8};
constexpr double kFloorSolidY[kFloorSolidPoints] = {3.4, 5.0, 6.0, 6.3, 6.6, 7.3, 7.75, -20.0, -20.0};

double Tongue::area_of(const TongueContour& shape) const {
  return settings_.blade_body > 0.5 ? thick_area(shape) : outline_area(shape);
}

// The thick tongue's front from its upper surface. The tip's cap is a half
// circle of tip_radius below the last section (its top), oriented by the
// blade (sections 20 to 23), turning through the front to its bottom; the
// free blade's underside lies below sections 23 back to 20, thickening from
// twice tip_radius to blade_thickness; the frenulum drops from its back
// straight down (jaw frame) onto the floor of the mouth.
void Tongue::derive_underside(const TongueContour& shape, TongueUnderside& out) const {
  constexpr int kCap = 3;
  constexpr int kUnder = 4;
  double radius = settings_.tip_radius;
  out.count = kCap + kUnder + 1;
  std::size_t tip = static_cast<std::size_t>(kSections - 1);
  std::size_t blade = tip - static_cast<std::size_t>(kUnder - 1);
  double tx = shape.x[tip] - shape.x[blade];
  double ty = shape.y[tip] - shape.y[blade];
  double norm = std::max(std::hypot(tx, ty), 1.0e-9);
  tx /= norm;
  ty /= norm;
  double nx = -ty;  // outward (up): left of root-to-tip
  double ny = tx;
  double cx = shape.x[tip] - radius * nx;
  double cy = shape.y[tip] - radius * ny;
  for (int point = 0; point < kCap; ++point) {
    double angle = kPi * static_cast<double>(point + 1) / static_cast<double>(kCap + 1);
    out.x[static_cast<std::size_t>(point)] = cx + radius * (std::cos(angle) * nx + std::sin(angle) * tx);
    out.y[static_cast<std::size_t>(point)] = cy + radius * (std::cos(angle) * ny + std::sin(angle) * ty);
  }
  for (int under = 0; under < kUnder; ++under) {
    std::size_t m = tip - static_cast<std::size_t>(under);
    double thickness = 2.0 * radius + (settings_.blade_thickness - 2.0 * radius) * static_cast<double>(under) /
                                          static_cast<double>(kUnder - 1);
    std::size_t o = static_cast<std::size_t>(kCap + under);
    out.x[o] = shape.x[m] - thickness * nx;
    out.y[o] = shape.y[m] - thickness * ny;
  }
  // The frenulum: below the blade's back, on the floor (or the blade's back
  // itself where it already lies on or below the floor).
  std::size_t back = static_cast<std::size_t>(kCap + kUnder - 1);
  std::size_t last = back + 1;
  double floor_y = floor_height(out.x[back]);
  out.x[last] = out.x[back];
  out.y[last] = std::min(out.y[back], floor_y);
}

// The floor solid's upper boundary height at `x` (jaw frame).
double Tongue::floor_height(double x) const {
  for (int i = 0; i + 1 < 7; ++i) {
    double ax = kFloorSolidX[i];
    double bx = kFloorSolidX[i + 1];
    if (x >= ax && x <= bx && bx > ax) {
      return kFloorSolidY[i] + (x - ax) / (bx - ax) * (kFloorSolidY[i + 1] - kFloorSolidY[i]);
    }
  }
  return x < kFloorSolidX[0] ? kFloorSolidY[0] : kFloorSolidY[6];
}

// The section an underside point hangs from (the cap's from the tip; -1
// for the frenulum, fixed).
int Tongue::underside_parent(int point) {
  constexpr int kCap = 3;
  constexpr int kUnder = 4;
  if (point < kCap) {
    return kSections - 1;
  }
  if (point < kCap + kUnder) {
    return kSections - 1 - (point - kCap);
  }
  return -1;
}

double Tongue::thick_area(const TongueContour& shape) const {
  TongueUnderside under;
  derive_underside(shape, under);
  std::array<double, kSections + kTongueUndersidePoints + 1> xs{};
  std::array<double, kSections + kTongueUndersidePoints + 1> ys{};
  std::size_t count = 0;
  for (int i = 0; i < kSections; ++i) {
    xs[count] = shape.x[static_cast<std::size_t>(i)];
    ys[count] = shape.y[static_cast<std::size_t>(i)];
    ++count;
  }
  for (int i = 0; i < under.count; ++i) {
    xs[count] = under.x[static_cast<std::size_t>(i)];
    ys[count] = under.y[static_cast<std::size_t>(i)];
    ++count;
  }
  // From the frenulum straight to the floor's back. (Through the floor's
  // front corner too, the outline crossed itself whenever the blade's back
  // lay behind that corner, and the hydrostat chased a wrong area: the body
  // crumpled.)
  xs[count] = kFloorBackX;
  ys[count] = kFloorBackY;
  ++count;
  double twice = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    std::size_t j = (i + 1) % count;
    twice += xs[i] * ys[j] - xs[j] * ys[i];
  }
  return 0.5 * std::abs(twice);
}

bool Tongue::inside_floor(double px, double py) const {
  bool inside = false;
  for (int i = 0; i < kFloorSolidPoints; ++i) {
    int j = (i + 1) % kFloorSolidPoints;
    double ay = kFloorSolidY[i];
    double by = kFloorSolidY[j];
    if ((ay > py) != (by > py)) {
      double crossing = kFloorSolidX[i] + (py - ay) * (kFloorSolidX[j] - kFloorSolidX[i]) / (by - ay);
      if (px < crossing) {
        inside = !inside;
      }
    }
  }
  return inside;
}

// The nearest point on the floor solid's upper boundary (its first seven
// points, the floor and the incisors' inner face).
void Tongue::nearest_floor_point(double px, double py, double& qx, double& qy) const {
  double best = 1.0e9;
  for (int i = 0; i + 1 < 7; ++i) {
    double ax = kFloorSolidX[i];
    double ay = kFloorSolidY[i];
    double dx = kFloorSolidX[i + 1] - ax;
    double dy = kFloorSolidY[i + 1] - ay;
    double span = dx * dx + dy * dy;
    double t = span > 0.0 ? std::clamp(((px - ax) * dx + (py - ay) * dy) / span, 0.0, 1.0) : 0.0;
    double cx = ax + t * dx;
    double cy = ay + t * dy;
    double gap = std::hypot(px - cx, py - cy);
    if (gap < best) {
      best = gap;
      qx = cx;
      qy = cy;
    }
  }
}

// The cap and underside against the solids. A point inside the hull moves
// the section it hangs from by the way out, and that section loses its
// speed along it (the tip may not enter the teeth or the ridge). A point
// inside the floor stays on the floor: the blade resting there squashes
// (its thickness yields), and nothing above it is pushed. (Pushing the
// sections, the blade's thickness lifted a resting tip 1.3 cm every
// substep and the body crumpled.)
void Tongue::stop_underside() {
  derive_underside(contour_, underside_);
  for (int point = 0; point < underside_.count; ++point) {
    int parent = underside_parent(point);
    if (parent < 1) {
      continue;
    }
    std::size_t p = static_cast<std::size_t>(point);
    std::size_t s = static_cast<std::size_t>(parent);
    double px = underside_.x[p];
    double py = underside_.y[p];
    double qx = px;
    double qy = py;
    if (inside_floor(px, py)) {
      underside_.y[p] = std::max(underside_.y[p], floor_height(px));
      continue;
    }
    if (!inside_hull(px, py)) {
      continue;
    }
    double nx = 0.0;
    double ny = 0.0;
    nearest_hull_point(px, py, qx, qy, nx, ny);
    double ox = qx - px;
    double oy = qy - py;
    double out = std::hypot(ox, oy);
    if (out <= 1.0e-9) {
      continue;
    }
    contour_.x[s] += ox;
    contour_.y[s] += oy;
    double ux = ox / out;
    double uy = oy / out;
    double into = vx_[s] * ux + vy_[s] * uy;
    if (into < 0.0) {
      vx_[s] -= into * ux;
      vy_[s] -= into * uy;
    }
  }
  derive_underside(contour_, underside_);
  rest_on_floor();
}

// The underside's points held on the floor where the blade rests on it.
void Tongue::rest_on_floor() {
  for (int point = 0; point < underside_.count; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    if (inside_floor(underside_.x[p], underside_.y[p])) {
      underside_.y[p] = std::max(underside_.y[p], floor_height(underside_.x[p]));
    }
  }
}

// Muscle `muscle`'s pull on `section` per unit activation (dyn).
void Tongue::muscle_pull(int muscle, int section, const TongueContour& shape, const std::array<double, 6>& origins,
                         double& fx, double& fy) const {
  const MuscleLayout& layout = kMuscleLayout[static_cast<std::size_t>(muscle)];
  double weight = reach_weight(layout, section, kSections) * settings_.muscle_force * layout.strength;
  std::size_t s = static_cast<std::size_t>(section);
  if (layout.origin == -1) {
    fx = weight * layout.direction_x;
    fy = weight * layout.direction_y;
    return;
  }
  if (layout.origin == -3) {
    std::size_t before = s > 0 ? s - 1 : s;
    std::size_t after = s + 1 < static_cast<std::size_t>(kSections) ? s + 1 : s;
    double tx = shape.x[after] - shape.x[before];
    double ty = shape.y[after] - shape.y[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    fx = weight * layout.direction_x * tx / norm;
    fy = weight * layout.direction_x * ty / norm;
    return;
  }
  if (layout.origin == -2) {
    // Along the surface's outward normal (left of root-to-tip), or inward.
    std::size_t before = s > 0 ? s - 1 : s;
    std::size_t after = s + 1 < static_cast<std::size_t>(kSections) ? s + 1 : s;
    double tx = shape.x[after] - shape.x[before];
    double ty = shape.y[after] - shape.y[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    fx = weight * layout.direction_x * -ty / norm;
    fy = weight * layout.direction_x * tx / norm;
    return;
  }
  double dx = origins[static_cast<std::size_t>(2 * layout.origin)] - shape.x[s];
  double dy = origins[static_cast<std::size_t>(2 * layout.origin + 1)] - shape.y[s];
  double length = std::max(std::hypot(dx, dy), 1.0e-6);
  fx = weight * dx / length;
  fy = weight * dy / length;
}

// Feedback: each muscle's drive is the target's pull on the sections it
// reaches (the error, less a share of the velocity) projected on its own
// pull; proportional plus integral, then the activation follows.
void Tongue::drive_muscles(double seconds, const std::array<double, 6>& origins) {
  double follow = 1.0 - std::exp(-seconds / std::max(settings_.activation_time, 1.0e-4));
  for (int muscle = 0; muscle < kMuscles; ++muscle) {
    double projected = 0.0;
    double weights = 0.0;
    for (int section = 1; section < kSections; ++section) {
      double fx = 0.0;
      double fy = 0.0;
      muscle_pull(muscle, section, contour_, origins, fx, fy);
      double norm = std::hypot(fx, fy);
      if (norm <= 0.0) {
        continue;
      }
      std::size_t s = static_cast<std::size_t>(section);
      const TongueContour& aim = settings_.feedback_aim > 0.5 ? aim_ : target_;
      double ex = aim.x[s] - contour_.x[s] - settings_.damping_gain * vx_[s];
      double ey = aim.y[s] - contour_.y[s] - settings_.damping_gain * vy_[s];
      projected += fx * ex + fy * ey;  // dyn cm (the pull along the error)
      weights += norm;
    }
    std::size_t m = static_cast<std::size_t>(muscle);
    double error = weights > 0.0 ? projected / weights : 0.0;  // cm along the muscle's pull
    double drive = feedforward_[m] + settings_.gain * error + integral_[m];
    double wanted = std::clamp(drive, 0.0, 1.0);
    // The integral only while it can act (no wind-up against a limit).
    bool pushing_out = (drive >= 1.0 && error > 0.0) || (drive <= 0.0 && error < 0.0);
    if (!pushing_out) {
      integral_[m] = std::clamp(integral_[m] + settings_.integral_rate * settings_.gain * error * seconds, -1.0, 1.0);
    }
    activation_[m] += follow * (wanted - activation_[m]);
  }
}

// Moves a shape's surface (not its root nor its tip) along its normals
// until its outline has `wanted` area (two passes; the shift per pass is the
// missing area over the surface's length).
void Tongue::shift_to_area(TongueContour& shape, double wanted) const {
  for (int pass = 0; pass < 2; ++pass) {
    double missing = wanted - area_of(shape);
    double length = 0.0;
    for (int i = 0; i + 1 < kSections; ++i) {
      std::size_t a = static_cast<std::size_t>(i);
      length += std::hypot(shape.x[a + 1] - shape.x[a], shape.y[a + 1] - shape.y[a]);
    }
    double shift = missing / std::max(length, 1.0e-6);
    TongueContour before = shape;
    for (int i = 1; i + 1 < kSections; ++i) {
      std::size_t a = static_cast<std::size_t>(i);
      double tx = before.x[a + 1] - before.x[a - 1];
      double ty = before.y[a + 1] - before.y[a - 1];
      double norm = std::max(std::hypot(tx, ty), 1.0e-9);
      // Outward: to the left of the root-to-tip direction.
      shape.x[a] += shift * -ty / norm;
      shape.y[a] += shift * tx / norm;
    }
  }
}

void Tongue::set_target(const TongueContour& target) {
  target_ = target;
  planned_area_ = area_of(target_);
  if (ready_ && settings_.spine <= 0.5 && settings_.skeleton <= 0.5) {
    shift_to_area(target_, held_area_);
    if (measured_ && settings_.statics <= 0.5) {
      solve_feedforward();
    }
  }
}

// The hydrostat: the area restored by moving the surface along its normals.
void Tongue::hold_volume() {
  shift_to_area(contour_, held_area_);
}

// A piece's pin strength, by the piece's own setting (a piece that is not a
// body is pinned throughout) or the role schedule.
void Tongue::set_pins(const std::array<double, 4>& pieces) {
  pinning_ = false;
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    int piece = section <= 5 ? 0 : section <= 14 ? 1 : section <= 20 ? 2 : 3;
    double weight = pieces[static_cast<std::size_t>(piece)];
    if (!piece_physical(section, settings_)) {
      weight = 1.0;
    }
    pin_[s] = std::clamp(weight, 0.0, 1.0);
    pinning_ = pinning_ || pin_[s] > 0.0;
  }
}

// The guard (the inverse tongue): each section's guard point lies a share
// 1 - guard_keep of the way from the planned surface to the hull, along the
// gap; its normal points from the planned surface toward the hull.
void Tongue::build_guard() {
  guard_.count = geometric_.count;
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double px = geometric_.x[s];
    double py = geometric_.y[s];
    double qx = 0.0;
    double qy = 0.0;
    double hull_nx = 0.0;
    double hull_ny = 0.0;
    nearest_hull_point(px, py, qx, qy, hull_nx, hull_ny);
    double dx = qx - px;
    double dy = qy - py;
    double gap = std::hypot(dx, dy);
    // The planned surface at or past the hull: the guard on the hull, its
    // normal the hull's (into the airway, reversed: toward the hull).
    bool open = gap > 1.0e-6 && (dx * hull_nx + dy * hull_ny) < 0.0;
    double nx = open ? dx / gap : -hull_nx;
    double ny = open ? dy / gap : -hull_ny;
    double share = open ? 1.0 - settings_.guard_keep : 1.0;
    guard_.x[s] = open ? px + share * dx : qx;
    guard_.y[s] = open ? py + share * dy : qy;
    guard_nx_[s] = nx;
    guard_ny_[s] = ny;
    guard_weight_[s] = settings_.roles > 0.5 ? shaping_[s] : 1.0;
  }
  guard_.x[0] = geometric_.x[0];
  guard_.y[0] = geometric_.y[0];
}

void Tongue::guard_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const {
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double weight = guard_weight_[s];
    if (weight <= 0.0) {
      continue;
    }
    double into = (contour_.x[s] - guard_.x[s]) * guard_nx_[s] + (contour_.y[s] - guard_.y[s]) * guard_ny_[s];
    if (into <= 0.0) {
      continue;
    }
    double speed = std::max(vx_[s] * guard_nx_[s] + vy_[s] * guard_ny_[s], 0.0);
    double push = weight * (settings_.guard_stiffness * into + settings_.guard_loss * speed);
    fx[s] -= push * guard_nx_[s];
    fy[s] -= push * guard_ny_[s];
  }
}

void Tongue::measure_intrusion() {
  intrusion_[0] = 0.0;
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double into = (contour_.x[s] - guard_.x[s]) * guard_nx_[s] + (contour_.y[s] - guard_.y[s]) * guard_ny_[s];
    intrusion_[s] = guard_weight_[s] > 0.0 ? std::max(into, 0.0) * guard_weight_[s] : 0.0;
  }
}

void Tongue::set_shaping(const std::array<double, 4>& shares) {
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    std::size_t piece = static_cast<std::size_t>(section <= 5 ? 0 : section <= 14 ? 1 : section <= 20 ? 2 : 3);
    shaping_[s] = std::clamp(shares[piece], 0.0, 1.0);
  }
}

// The feedback's aim, section by section: the motors' posture now (at the
// tongue's own area) where the piece shapes a vowel, the look-ahead target
// where it strikes or holds a channel (feedback_aim 2); the posture
// everywhere (1).
void Tongue::build_aim() {
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double toward_posture = settings_.feedback_aim > 1.5 ? shaping_[s] : 1.0;
    aim_.x[s] = target_.x[s] + toward_posture * (anchor_.x[s] - target_.x[s]);
    aim_.y[s] = target_.y[s] + toward_posture * (anchor_.y[s] - target_.y[s]);
  }
  aim_.count = target_.count;
}

void Tongue::set_impedance(const std::array<double, 4>& stiffness, const std::array<double, 4>& damping) {
  holding_ = false;
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    std::size_t piece = static_cast<std::size_t>(section <= 5 ? 0 : section <= 14 ? 1 : section <= 20 ? 2 : 3);
    impedance_[s] = std::max(stiffness[piece], 0.0);
    impedance_damping_[s] = std::max(damping[piece], 0.0);
    holding_ = holding_ || impedance_[s] > 0.0;
  }
}

// Co-contraction: each section held toward its target posture by its
// piece's stiffness, with its damping (fraction of critical). The target is
// the one the feedback aims at; a target past the hull presses into it.
void Tongue::impedance_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const {
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double stiffness = impedance_[s];
    if (stiffness <= 0.0) {
      continue;
    }
    double along = static_cast<double>(section) / static_cast<double>(kSections - 1);
    double own_mass = settings_.section_mass * (1.0 + (settings_.tip_mass - 1.0) * along);
    double loss = 2.0 * impedance_damping_[s] * std::sqrt(stiffness * own_mass);
    if (settings_.role_anchor > 0.5) {
      fx[s] += stiffness * (anchor_.x[s] - contour_.x[s]) + loss * (geometric_vx_[s] - vx_[s]);
      fy[s] += stiffness * (anchor_.y[s] - contour_.y[s]) + loss * (geometric_vy_[s] - vy_[s]);
    } else {
      fx[s] += stiffness * (target_.x[s] - contour_.x[s]) - loss * vx_[s];
      fy[s] += stiffness * (target_.y[s] - contour_.y[s]) - loss * vy_[s];
    }
  }
}

void Tongue::set_geometric(const TongueContour& shape, double seconds) {
  for (int section = 0; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    geometric_vx_[s] = geometric_set_ ? (shape.x[s] - geometric_.x[s]) / seconds : 0.0;
    geometric_vy_[s] = geometric_set_ ? (shape.y[s] - geometric_.y[s]) / seconds : 0.0;
  }
  geometric_ = shape;
  geometric_set_ = true;
  if (settings_.role_anchor > 0.5) {
    anchor_ = shape;
    shift_to_area(anchor_, held_area_);
  }
}

// The pin: each section pulled toward the geometric outline (and its speed)
// by a critically damped spring, in proportion to its pin strength. Pinned,
// a piece follows the geometric tongue within about a millimetre; the bodies
// beside it feel it through stretch and bending.
void Tongue::pin_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const {
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double weight = pin_[s];
    if (weight <= 0.0) {
      continue;
    }
    double along = static_cast<double>(section) / static_cast<double>(kSections - 1);
    double own_mass = settings_.section_mass * (1.0 + (settings_.tip_mass - 1.0) * along);
    double stiffness = weight * settings_.pin_stiffness;
    double loss = 2.0 * std::sqrt(stiffness * own_mass);
    fx[s] += stiffness * (geometric_.x[s] - contour_.x[s]) + loss * (geometric_vx_[s] - vx_[s]);
    fy[s] += stiffness * (geometric_.y[s] - contour_.y[s]) + loss * (geometric_vy_[s] - vy_[s]);
  }
}

// The hydrostat with pinned pieces: only the bodies' surface moves (each
// section by its share of freedom, 1 - pin) to restore the rest area; the
// pinned pieces cannot give way, so the bodies bulge or yield for them.
void Tongue::hold_volume_weighted() {
  for (int pass = 0; pass < 2; ++pass) {
    double missing = held_area_ - area_of(contour_);
    double length = 0.0;
    for (int i = 1; i + 1 < kSections; ++i) {
      std::size_t a = static_cast<std::size_t>(i);
      double free_share = 1.0 - pin_[a];
      length += free_share * 0.5 * std::hypot(contour_.x[a + 1] - contour_.x[a - 1], contour_.y[a + 1] - contour_.y[a - 1]);
    }
    if (length <= 1.0e-6) {
      return;
    }
    double shift = missing / length;
    TongueContour before = contour_;
    for (int i = 1; i + 1 < kSections; ++i) {
      std::size_t a = static_cast<std::size_t>(i);
      double free_share = 1.0 - pin_[a];
      if (free_share <= 0.0) {
        continue;
      }
      double tx = before.x[a + 1] - before.x[a - 1];
      double ty = before.y[a + 1] - before.y[a - 1];
      double norm = std::max(std::hypot(tx, ty), 1.0e-9);
      contour_.x[a] += free_share * shift * -ty / norm;
      contour_.y[a] += free_share * shift * tx / norm;
    }
  }
}

// Contact: a section past the hull (beyond contact_depth) is returned to it
// and loses its speed into it. The hull runs from the larynx up the
// pharyngeal wall and forward over the palate to the upper incisors' edge
// (the lips past it are not a wall the tongue can reach from behind); the
// airway is on its right. At a corner the side is judged by the bisector of
// the two segments' normals: by the nearest segment's alone, a tip below and
// behind the incisors' edge read as inside the tooth (the edge-to-lip
// segment's normal) and was pushed back; no /dh/ could reach the teeth.
// The hull's nearest point to a section, and its normal toward the airway
// (at a corner, the bisector of the two segments' normals).
void Tongue::nearest_hull(int section, double& qx, double& qy, double& nx_out, double& ny_out) const {
  std::size_t s = static_cast<std::size_t>(section);
  nearest_hull_point(contour_.x[s], contour_.y[s], qx, qy, nx_out, ny_out);
}

void Tongue::nearest_hull_point(double px, double py, double& qx, double& qy, double& nx_out, double& ny_out) const {
  double best = 1.0e9;
  qx = 0.0;
  qy = 0.0;
  nx_out = 0.0;
  ny_out = 0.0;
  for (int i = 0; i + 1 < hull_count_; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    double ax = hull_x_[a];
    double ay = hull_y_[a];
    double dx = hull_x_[a + 1] - ax;
    double dy = hull_y_[a + 1] - ay;
    double span = dx * dx + dy * dy;
    if (span <= 0.0) {
      continue;
    }
    double t = std::clamp(((px - ax) * dx + (py - ay) * dy) / span, 0.0, 1.0);
    double cx = ax + t * dx;
    double cy = ay + t * dy;
    double gap = std::hypot(px - cx, py - cy);
    if (gap < best) {
      double length = std::sqrt(span);
      double nx = dy / length;  // the right normal: toward the airway
      double ny = -dx / length;
      std::size_t neighbour = a;
      if (t <= 0.0 && a > 0) {
        neighbour = a - 1;
      } else if (t >= 1.0 && a + 2 < static_cast<std::size_t>(hull_count_)) {
        neighbour = a + 1;
      }
      if (neighbour != a) {
        double ex = hull_x_[neighbour + 1] - hull_x_[neighbour];
        double ey = hull_y_[neighbour + 1] - hull_y_[neighbour];
        double other = std::max(std::hypot(ex, ey), 1.0e-9);
        nx += ey / other;
        ny += -ex / other;
        double norm = std::max(std::hypot(nx, ny), 1.0e-9);
        nx /= norm;
        ny /= norm;
      }
      best = gap;
      qx = cx;
      qy = cy;
      nx_out = nx;
      ny_out = ny;
    }
  }
  if (best >= 1.0e8) {
    nx_out = 0.0;
    ny_out = 0.0;
  }
}

// The closed hull's corners: its points from the larynx to the incisors'
// edge, then far above the edge, far above and behind, and behind the first.
void Tongue::hull_corner(int index, double& x, double& y) const {
  constexpr double kFar = 40.0;
  std::size_t first = 0;
  std::size_t last = static_cast<std::size_t>(hull_count_ - 1);
  if (index < hull_count_) {
    x = hull_x_[static_cast<std::size_t>(index)];
    y = hull_y_[static_cast<std::size_t>(index)];
  } else if (index == hull_count_) {
    x = hull_x_[last];
    y = kFar;
  } else if (index == hull_count_ + 1) {
    x = hull_x_[first] - kFar;
    y = kFar;
  } else {
    x = hull_x_[first] - kFar;
    y = hull_y_[first];
  }
}

// Ray casting on the closed hull (the skull's side).
bool Tongue::inside_hull(double px, double py) const {
  if (hull_count_ < 2) {
    return false;
  }
  int corners = hull_count_ + 3;
  bool inside = false;
  for (int i = 0; i < corners; ++i) {
    double ax = 0.0;
    double ay = 0.0;
    double bx = 0.0;
    double by = 0.0;
    hull_corner(i, ax, ay);
    hull_corner((i + 1) % corners, bx, by);
    if ((ay > py) != (by > py)) {
      double crossing = ax + (py - ay) * (bx - ax) / (by - ay);
      if (px < crossing) {
        inside = !inside;
      }
    }
  }
  return inside;
}

double Tongue::lift() const {
  double sum = 0.0;
  for (int section = 8; section <= 16; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    sum += contour_.y[s] - rest_.y[s];
  }
  return sum / 9.0;
}

void Tongue::measure_penetration() {
  penetration_[0] = 0.0;
  if (settings_.blade_body > 0.5) {
    derive_underside(contour_, underside_);
    rest_on_floor();
  }
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    penetration_[s] = 0.0;
    if (inside_hull(contour_.x[s], contour_.y[s])) {
      double qx = 0.0;
      double qy = 0.0;
      double nx = 0.0;
      double ny = 0.0;
      nearest_hull(section, qx, qy, nx, ny);
      penetration_[s] = std::hypot(contour_.x[s] - qx, contour_.y[s] - qy);
    }
  }
  bool thick = settings_.blade_body > 0.5 || settings_.spine > 0.5 || settings_.skeleton > 0.5;
  for (int point = 0; thick && point < underside_.count; ++point) {
    int parent = underside_parent(point);
    if (parent < 1) {
      continue;
    }
    std::size_t p = static_cast<std::size_t>(point);
    double qx = underside_.x[p];
    double qy = underside_.y[p];
    if (inside_hull(underside_.x[p], underside_.y[p])) {
      double nx = 0.0;
      double ny = 0.0;
      nearest_hull_point(underside_.x[p], underside_.y[p], qx, qy, nx, ny);
    }
    std::size_t s = static_cast<std::size_t>(parent);
    penetration_[s] = std::max(penetration_[s], std::hypot(underside_.x[p] - qx, underside_.y[p] - qy));
  }
}

void Tongue::stop_at_hull(double depth) {
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    if (pin_[s] > 0.5) {
      continue;
    }
    double qx = 0.0;
    double qy = 0.0;
    double nx = 0.0;
    double ny = 0.0;
    nearest_hull(section, qx, qy, nx, ny);
    if (settings_.solid_hull > 0.5) {
      // A section inside the solid goes back to its nearest hull point (less
      // `depth`), whichever face that is; it loses its speed into the solid.
      if (!inside_hull(contour_.x[s], contour_.y[s])) {
        continue;
      }
      double ox = qx - contour_.x[s];
      double oy = qy - contour_.y[s];
      double out = std::hypot(ox, oy);
      if (out <= depth || out <= 1.0e-9) {
        continue;
      }
      double ux = ox / out;  // out of the solid, toward the airway
      double uy = oy / out;
      contour_.x[s] = qx - depth * ux;
      contour_.y[s] = qy - depth * uy;
      double into = vx_[s] * ux + vy_[s] * uy;
      if (into < 0.0) {
        vx_[s] -= into * ux;
        vy_[s] -= into * uy;
      }
      continue;
    }
    if (nx == 0.0 && ny == 0.0) {
      continue;
    }
    double side = (contour_.x[s] - qx) * nx + (contour_.y[s] - qy) * ny;
    if (side < -depth) {
      contour_.x[s] = qx - depth * nx;
      contour_.y[s] = qy - depth * ny;
      double into = vx_[s] * nx + vy_[s] * ny;
      if (into < 0.0) {
        vx_[s] -= into * nx;
        vy_[s] -= into * ny;
      }
    }
  }
}

// Soft contact (contact_stiffness > 0): the mucosa between the tongue and
// the hull as a compliant, lossy layer. A section pressed p past the hull
// is pushed back by stiffness p^1.5 (Hertzian, as the folds' contact) and
// loses its normal speed through contact_loss while it presses. (With the
// hard stop alone, a tip held on the ridge was snapped back every substep
// and the seal chattered between 0.000 and 0.015 cm^2 from one millisecond
// to the next.)
void Tongue::contact_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const {
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    if (pin_[s] > 0.5) {
      continue;
    }
    double qx = 0.0;
    double qy = 0.0;
    double nx = 0.0;
    double ny = 0.0;
    nearest_hull(section, qx, qy, nx, ny);
    double pressed = -((contour_.x[s] - qx) * nx + (contour_.y[s] - qy) * ny);
    if (pressed <= 0.0) {
      continue;
    }
    double speed = vx_[s] * nx + vy_[s] * ny;
    double push = settings_.contact_stiffness * pressed * std::sqrt(pressed) - settings_.contact_loss * speed;
    push = std::max(push, 0.0);
    fx[s] += push * nx;
    fy[s] += push * ny;
  }
}

// The attachments (genial tubercle, styloid, hyoid) and the hull, in the
// jaw's frame; the root's base held beside the epiglottis.
void Tongue::prepare_frame(const Articulators& frame, const MidsagittalTract& anatomy) {
  origins_[0] = kGenialX;
  origins_[1] = kGenialY;
  anatomy.to_jaw_frame(frame, kStyloidX, kStyloidY, origins_[2], origins_[3]);
  double drop = std::clamp(frame.larynx_drop, 0.0, 2.0);
  anatomy.to_jaw_frame(frame, kHyoidX, kHyoidY - drop, origins_[4], origins_[5]);
  // The hull up to the upper incisors' edge (its lowest point in front).
  int count = std::min(anatomy.hull_points(), static_cast<int>(hull_x_.size()));
  int last = count - 1;
  double lowest = 1.0e9;
  for (int i = 0; i < count; ++i) {
    if (anatomy.hull_x(i) > 7.0 && anatomy.hull_x(i) < kIncisorReach && anatomy.hull_y(i) < lowest) {
      lowest = anatomy.hull_y(i);
      last = i;
    }
  }
  hull_count_ = last + 1;
  for (int i = 0; i < hull_count_; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    anatomy.to_jaw_frame(frame, anatomy.hull_x(i), anatomy.hull_y(i), hull_x_[a], hull_y_[a]);
  }
  double root_x = 0.0;
  double root_y = 0.0;
  anatomy.root_base(frame, root_x, root_y);
  anatomy.to_jaw_frame(frame, root_x, root_y, root_x_, root_y_);
}

// One substep of the body under the present activations: support, muscles,
// elasticity; then the hydrostat and (with `touching`) contact.
void Tongue::advance(double h, bool touching) {
  contour_.x[0] = root_x_;
  contour_.y[0] = root_y_;
  vx_[0] = 0.0;
  vy_[0] = 0.0;
  double mass = settings_.section_mass;
  double support_c = 2.0 * settings_.support_loss * std::sqrt(mass * settings_.support);
  double elastic_c = 2.0 * settings_.elastic_loss * std::sqrt(mass * settings_.elasticity);
  // The sections lighten toward the tip (tip_mass of the root's), linearly.
  std::array<double, kSections> fx{};
  std::array<double, kSections> fy{};
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    // Support: back toward the rest shape.
    double pliant = pliancy(section, kSections, settings_.tip_pliancy);
    fx[s] += pliant * (-settings_.support * (contour_.x[s] - rest_.x[s]) - support_c * vx_[s]);
    fy[s] += pliant * (-settings_.support * (contour_.y[s] - rest_.y[s]) - support_c * vy_[s]);
    // Muscles.
    for (int muscle = 0; muscle < kMuscles; ++muscle) {
      double mx = 0.0;
      double my = 0.0;
      muscle_pull(muscle, section, contour_, origins_, mx, my);
      double a = activation_[static_cast<std::size_t>(muscle)];
      fx[s] += a * mx;
      fy[s] += a * my;
    }
  }
  if (touching && settings_.contact_stiffness > 0.0) {
    contact_forces(fx, fy);
  }
  bool pinned = pinning_ && touching;
  if (pinned) {
    pin_forces(fx, fy);
  }
  if (holding_ && touching) {
    impedance_forces(fx, fy);
  }
  if (guarding_ && touching && settings_.guard_stiffness > 0.0) {
    guard_forces(fx, fy);
  }
  // Elasticity: stretch and its loss between neighbours, and bending.
  for (int i = 0; i + 1 < kSections; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    double dx = contour_.x[a + 1] - contour_.x[a];
    double dy = contour_.y[a + 1] - contour_.y[a];
    double length = std::max(std::hypot(dx, dy), 1.0e-9);
    double ux = dx / length;
    double uy = dy / length;
    double closing = (vx_[a + 1] - vx_[a]) * ux + (vy_[a + 1] - vy_[a]) * uy;
    // The blade stretches more freely than the body (it lengthens as the
    // tongue protrudes: held to its rest length, the tip could not reach the
    // incisors for /dh/).
    double stretch_pliancy = pliancy(i, kSections, settings_.tip_stretch);
    double tension = stretch_pliancy * (settings_.elasticity * (length - rest_spacing_[a]) + elastic_c * closing);
    fx[a] += tension * ux;
    fy[a] += tension * uy;
    fx[a + 1] -= tension * ux;
    fy[a + 1] -= tension * uy;
  }
  for (int i = 1; i + 1 < kSections; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    double bend = settings_.bending * pliancy(i, kSections, settings_.tip_pliancy);
    double bx = 2.0 * contour_.x[a] - contour_.x[a - 1] - contour_.x[a + 1] - rest_bend_x_[a];
    double by = 2.0 * contour_.y[a] - contour_.y[a - 1] - contour_.y[a + 1] - rest_bend_y_[a];
    fx[a] -= 2.0 * bend * bx;
    fy[a] -= 2.0 * bend * by;
    fx[a - 1] += bend * bx;
    fy[a - 1] += bend * by;
    fx[a + 1] += bend * bx;
    fy[a + 1] += bend * by;
  }
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double along = static_cast<double>(section) / static_cast<double>(kSections - 1);
    double own_mass = mass * (1.0 + (settings_.tip_mass - 1.0) * along);
    vx_[s] += h * fx[s] / own_mass;
    vy_[s] += h * fy[s] / own_mass;
    contour_.x[s] += h * vx_[s];
    contour_.y[s] += h * vy_[s];
  }
  if (pinned) {
    hold_volume_weighted();
  } else {
    hold_volume();
  }
  if (touching) {
    // Soft contact keeps a hard stop far behind it (no tunnelling through the hull).
    stop_at_hull(settings_.contact_stiffness > 0.0 ? 4.0 * settings_.contact_depth : settings_.contact_depth);
    if (settings_.blade_body > 0.5) {
      stop_underside();
    }
  }
}

// Each muscle's effect on the shape: switched on alone (activation 1, no
// contact, the frame as given) the body settles over kSettle seconds; its
// displacement from the settled rest is the muscle's influence. The
// feedforward then asks, for each target, which activations (0..1) hold the
// nearest shape, superposing the influences (least squares).
void Tongue::measure_influences() {
  constexpr double kSettle = 0.6;
  constexpr double kProbeActivation = 0.1;
  double h = 1.0e-3 / static_cast<double>(kSubsteps);
  int steps = static_cast<int>(kSettle / h);
  TongueContour placed = contour_;
  TongueContour settled = rest_;
  // The rest itself settled (the root may sit off the rest outline).
  contour_ = rest_;
  vx_.fill(0.0);
  vy_.fill(0.0);
  activation_.fill(0.0);
  for (int step = 0; step < steps; ++step) {
    advance(h, false);
  }
  settled = contour_;
  for (int muscle = 0; muscle < kMuscles; ++muscle) {
    contour_ = settled;
    vx_.fill(0.0);
    vy_.fill(0.0);
    activation_.fill(0.0);
    // (Measured small and scaled: at full activation a muscle alone carries
    // the tongue far, its pull turns and the hydrostat bends the response;
    // superposed, those shapes held 0.9 cm from a target they predicted to
    // 0.2.)
    activation_[static_cast<std::size_t>(muscle)] = kProbeActivation;
    for (int step = 0; step < steps; ++step) {
      advance(h, false);
    }
    for (int section = 0; section < kSections; ++section) {
      std::size_t s = static_cast<std::size_t>(section);
      influence_x_[static_cast<std::size_t>(muscle)][s] = (contour_.x[s] - settled.x[s]) / kProbeActivation;
      influence_y_[static_cast<std::size_t>(muscle)][s] = (contour_.y[s] - settled.y[s]) / kProbeActivation;
    }
  }
  settled_ = settled;
  contour_ = placed;
  vx_.fill(0.0);
  vy_.fill(0.0);
  activation_.fill(0.0);
  measured_ = true;
}

// The feedforward: activations in 0..1 minimising the distance between the
// settled rest plus their influences and the target (projected gradient on
// the kMuscles normal equations).
void Tongue::solve_feedforward() {
  std::array<std::array<double, kMuscles>, kMuscles> normal{};
  std::array<double, kMuscles> right{};
  for (int m = 0; m < kMuscles; ++m) {
    std::size_t a = static_cast<std::size_t>(m);
    for (int section = 1; section < kSections; ++section) {
      std::size_t s = static_cast<std::size_t>(section);
      double wx = target_.x[s] - settled_.x[s];
      double wy = target_.y[s] - settled_.y[s];
      right[a] += influence_x_[a][s] * wx + influence_y_[a][s] * wy;
      for (int n = 0; n < kMuscles; ++n) {
        std::size_t b = static_cast<std::size_t>(n);
        normal[a][b] += influence_x_[a][s] * influence_x_[b][s] + influence_y_[a][s] * influence_y_[b][s];
      }
    }
  }
  // Coordinate descent: each muscle in turn to its own best value given the
  // others, held in 0..1 (bounded least squares; converges for this
  // positive semidefinite system).
  for (int sweep = 0; sweep < 60; ++sweep) {
    for (int m = 0; m < kMuscles; ++m) {
      std::size_t a = static_cast<std::size_t>(m);
      if (normal[a][a] <= 1.0e-12) {
        continue;
      }
      double slope = -right[a];
      for (int n = 0; n < kMuscles; ++n) {
        slope += normal[a][static_cast<std::size_t>(n)] * feedforward_[static_cast<std::size_t>(n)];
      }
      feedforward_[a] = std::clamp(feedforward_[a] - slope / normal[a][a], 0.0, 1.0);
    }
  }
  // How well the muscles can hold the target at all (cm, rms).
  double sum = 0.0;
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double x = settled_.x[s];
    double y = settled_.y[s];
    for (int m = 0; m < kMuscles; ++m) {
      std::size_t a = static_cast<std::size_t>(m);
      x += feedforward_[a] * influence_x_[a][s];
      y += feedforward_[a] * influence_y_[a][s];
    }
    sum += (x - target_.x[s]) * (x - target_.x[s]) + (y - target_.y[s]) * (y - target_.y[s]);
  }
  reach_error_ = std::sqrt(sum / static_cast<double>(kSections - 1));
}

// The tissue's passive forces on a resting shape (no velocities): support
// toward the rest shape, stretch between neighbours and bending, as
// advance() applies them.
void Tongue::passive_forces(const TongueContour& shape, std::array<double, kSections>& fx,
                            std::array<double, kSections>& fy) const {
  fx.fill(0.0);
  fy.fill(0.0);
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double pliant = pliancy(section, kSections, settings_.tip_pliancy);
    fx[s] += pliant * -settings_.support * (shape.x[s] - rest_.x[s]);
    fy[s] += pliant * -settings_.support * (shape.y[s] - rest_.y[s]);
  }
  for (int i = 0; i + 1 < kSections; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    double dx = shape.x[a + 1] - shape.x[a];
    double dy = shape.y[a + 1] - shape.y[a];
    double length = std::max(std::hypot(dx, dy), 1.0e-9);
    double tension = pliancy(i, kSections, settings_.tip_stretch) * settings_.elasticity * (length - rest_spacing_[a]);
    fx[a] += tension * dx / length;
    fy[a] += tension * dy / length;
    fx[a + 1] -= tension * dx / length;
    fy[a + 1] -= tension * dy / length;
  }
  for (int i = 1; i + 1 < kSections; ++i) {
    std::size_t a = static_cast<std::size_t>(i);
    double bend = settings_.bending * pliancy(i, kSections, settings_.tip_pliancy);
    double bx = 2.0 * shape.x[a] - shape.x[a - 1] - shape.x[a + 1] - rest_bend_x_[a];
    double by = 2.0 * shape.y[a] - shape.y[a - 1] - shape.y[a + 1] - rest_bend_y_[a];
    fx[a] -= 2.0 * bend * bx;
    fy[a] -= 2.0 * bend * by;
    fx[a - 1] += bend * bx;
    fy[a - 1] += bend * by;
    fx[a + 1] += bend * bx;
    fy[a + 1] += bend * by;
  }
}

// The feedforward by inverse statics: at the target shape itself (not a
// linearization at rest), the activations in 0..1 whose pulls, with the
// hydrostat's pressure (a uniform push along the surface's outward normal,
// either sign), best balance the tissue's passive forces there; bounded
// least squares by coordinate descent, warm-started from the last target's.
// A target inside the hull is held as if the hull were absent: the muscles
// press into the contact.
void Tongue::solve_statics() {
  passive_forces(target_, statics_passive_x_, statics_passive_y_);
  for (int m = 0; m < kMuscles; ++m) {
    std::size_t a = static_cast<std::size_t>(m);
    for (int section = 1; section < kSections; ++section) {
      std::size_t s = static_cast<std::size_t>(section);
      muscle_pull(m, section, target_, origins_, statics_column_x_[a][s], statics_column_y_[a][s]);
    }
  }
  // The pressure's column: the outward normal times the section's share of
  // the surface, scaled by muscle_force (its unknown is then of order 1).
  std::size_t p = static_cast<std::size_t>(kMuscles);
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    std::size_t before = s - 1;
    std::size_t after = s + 1 < static_cast<std::size_t>(kSections) ? s + 1 : s;
    double tx = target_.x[after] - target_.x[before];
    double ty = target_.y[after] - target_.y[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    double share = 0.5 * norm * settings_.muscle_force;
    statics_column_x_[p][s] = share * -ty / norm;
    statics_column_y_[p][s] = share * tx / norm;
  }
  for (std::size_t a = 0; a <= p; ++a) {
    double right = 0.0;
    for (int section = 1; section < kSections; ++section) {
      std::size_t s = static_cast<std::size_t>(section);
      right -= statics_column_x_[a][s] * statics_passive_x_[s] + statics_column_y_[a][s] * statics_passive_y_[s];
    }
    statics_right_[a] = right;
    for (std::size_t b = 0; b <= p; ++b) {
      double sum = 0.0;
      for (int section = 1; section < kSections; ++section) {
        std::size_t s = static_cast<std::size_t>(section);
        sum += statics_column_x_[a][s] * statics_column_x_[b][s] + statics_column_y_[a][s] * statics_column_y_[b][s];
      }
      statics_normal_[a][b] = sum;
    }
  }
  for (int sweep = 0; sweep < 60; ++sweep) {
    for (std::size_t a = 0; a <= p; ++a) {
      if (statics_normal_[a][a] <= 1.0e-12) {
        continue;
      }
      double slope = -statics_right_[a];
      for (std::size_t b = 0; b < p; ++b) {
        slope += statics_normal_[a][b] * feedforward_[b];
      }
      slope += statics_normal_[a][p] * pressure_;
      if (a == p) {
        pressure_ = std::clamp(pressure_ - slope / statics_normal_[a][a], -4.0, 4.0);
      } else {
        feedforward_[a] = std::clamp(feedforward_[a] - slope / statics_normal_[a][a], 0.0, 1.0);
      }
    }
  }
  // The force left unbalanced, as a displacement against the support (cm rms).
  double sum = 0.0;
  for (int section = 1; section < kSections; ++section) {
    std::size_t s = static_cast<std::size_t>(section);
    double x = statics_passive_x_[s] + pressure_ * statics_column_x_[p][s];
    double y = statics_passive_y_[s] + pressure_ * statics_column_y_[p][s];
    for (std::size_t a = 0; a < p; ++a) {
      x += feedforward_[a] * statics_column_x_[a][s];
      y += feedforward_[a] * statics_column_y_[a][s];
    }
    sum += x * x + y * y;
  }
  reach_error_ = std::sqrt(sum / static_cast<double>(kSections - 1)) / settings_.support;
}

// ---- The spine tongue --------------------------------------------------

// The spine's distance below the upper surface at a station: thin at the
// root, body_thickness over the body, thinning over the blade to the tip's
// radius.
double Tongue::spine_up(int station) const {
  double along = static_cast<double>(station) / static_cast<double>(kSections - 1);
  double body = settings_.body_thickness;
  double tip = settings_.tip_radius;
  if (along < 0.15) {
    return 0.6 + (body - 0.6) * along / 0.15;
  }
  if (along < 0.65) {
    return body;
  }
  double share = (along - 0.65) / 0.35;
  return body + (tip - body) * share;
}

// Below the spine to the underside (the free blade only): the tip's radius
// at the tip, thickening backward.
double Tongue::spine_down(int station) const {
  int from_tip = kSections - 1 - station;
  return settings_.tip_radius + 0.12 * static_cast<double>(from_tip);
}

int Tongue::spine_base_segments() const {
  return std::clamp(static_cast<int>(std::lround(settings_.spine_base_end)), 2, kSections - 2);
}

int Tongue::spine_tip_root() const {
  return std::clamp(static_cast<int>(std::lround(settings_.spine_tip_start)), 1, spine_base_segments() - 1);
}

// A shape's beam bends: its spine (the surface less the spine's depth along
// the surface's normals), then the segments' headings; the base beam's
// bends from the clamp, the tip beam's from the base beam's segment at its
// root. (The spine's segment lengths are the tongue's own; only the
// headings are taken from the shape.)
void Tongue::spine_bends_of(const TongueContour& shape, std::array<double, kSections>& base,
                            std::array<double, kSections>& tip) const {
  std::array<double, kSections> sx{};
  std::array<double, kSections> sy{};
  for (int station = 0; station < kSections; ++station) {
    std::size_t i = static_cast<std::size_t>(station);
    std::size_t before = i > 0 ? i - 1 : i;
    std::size_t after = i + 1 < static_cast<std::size_t>(kSections) ? i + 1 : i;
    double tx = shape.x[after] - shape.x[before];
    double ty = shape.y[after] - shape.y[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    double depth = spine_up(station);
    sx[i] = shape.x[i] - depth * -ty / norm;
    sy[i] = shape.y[i] - depth * tx / norm;
  }
  std::array<double, kSections> heading{};
  for (int segment = 0; segment + 1 < kSections; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    heading[g] = std::atan2(sy[g + 1] - sy[g], sx[g + 1] - sx[g]);
  }
  int segments = spine_base_segments();
  int root = spine_tip_root();
  double previous = clamp_heading_;
  for (int segment = 0; segment < segments; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    base[g] = std::remainder(heading[g] - previous, 2.0 * kPi);
    previous = heading[g];
  }
  previous = heading[static_cast<std::size_t>(root)];
  for (int segment = root; segment + 1 < kSections; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    tip[g] = std::remainder(heading[g] - previous, 2.0 * kPi);
    previous = heading[g];
  }
}

// The beams' stations from their bends: the base beam from the anchor (the
// root base plus the rest offset), the tip beam from the base beam's station
// at its root; the spine is the base beam to the tip beam's root, the tip
// beam past the base beam's end, and between them a blend.
void Tongue::spine_positions(const std::array<double, kSections>& base, const std::array<double, kSections>& tip,
                             std::array<double, kSections>& sx, std::array<double, kSections>& sy,
                             std::array<double, kSections>& bx, std::array<double, kSections>& by,
                             std::array<double, kSections>& tx, std::array<double, kSections>& ty) const {
  int segments = spine_base_segments();
  int root = spine_tip_root();
  bx[0] = root_x_ + anchor_dx_;
  by[0] = root_y_ + anchor_dy_;
  double heading = clamp_heading_;
  std::array<double, kSections> base_heading{};
  for (int segment = 0; segment < segments; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    heading += base[g];
    base_heading[g] = heading;
    bx[g + 1] = bx[g] + spine_length_[g] * std::cos(heading);
    by[g + 1] = by[g] + spine_length_[g] * std::sin(heading);
  }
  std::size_t r = static_cast<std::size_t>(root);
  tx[r] = bx[r];
  ty[r] = by[r];
  heading = base_heading[r];
  for (int segment = root; segment + 1 < kSections; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    heading += tip[g];
    tx[g + 1] = tx[g] + spine_length_[g] * std::cos(heading);
    ty[g + 1] = ty[g] + spine_length_[g] * std::sin(heading);
  }
  for (int station = 0; station < kSections; ++station) {
    std::size_t i = static_cast<std::size_t>(station);
    if (station <= root) {
      sx[i] = bx[i];
      sy[i] = by[i];
    } else if (station <= segments) {
      double share = static_cast<double>(station - root) / static_cast<double>(segments - root);
      sx[i] = bx[i] + share * (tx[i] - bx[i]);
      sy[i] = by[i] + share * (ty[i] - by[i]);
    } else {
      sx[i] = tx[i];
      sy[i] = ty[i];
    }
  }
}

// A shape's spine: the surface less the spine's depth along its normals.
void Tongue::spine_points_of(const TongueContour& shape, std::array<double, kSections>& sx,
                             std::array<double, kSections>& sy) const {
  for (int station = 0; station < kSections; ++station) {
    std::size_t i = static_cast<std::size_t>(station);
    std::size_t before = i > 0 ? i - 1 : i;
    std::size_t after = i + 1 < static_cast<std::size_t>(kSections) ? i + 1 : i;
    double tx = shape.x[after] - shape.x[before];
    double ty = shape.y[after] - shape.y[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    double depth = spine_up(station);
    sx[i] = shape.x[i] - depth * -ty / norm;
    sy[i] = shape.y[i] - depth * tx / norm;
  }
}

// The goals' reach: three damped least-squares steps on the base beam's
// goal bends toward the shape's spine at the base beam's end, then on the
// tip beam's toward it at the tip (the turn of joint k moves a station by
// the perpendicular of its lever from the joint).
void Tongue::spine_reach_goals(const std::array<double, kSections>& goal_x,
                               const std::array<double, kSections>& goal_y) {
  double damping = settings_.spine_reach_damping;
  int steps = std::max(1, static_cast<int>(std::lround(settings_.spine_reach_steps)));
  int segments = spine_base_segments();
  int root = spine_tip_root();
  std::array<double, kSections> sx{};
  std::array<double, kSections> sy{};
  std::array<double, kSections> bx{};
  std::array<double, kSections> by{};
  std::array<double, kSections> px{};
  std::array<double, kSections> py{};
  for (int beam = settings_.spine_reach > 1.5 ? 1 : 0; beam < 2; ++beam) {
    bool tip_beam = beam == 1;
    std::size_t station = static_cast<std::size_t>(tip_beam ? kSections - 1 : segments);
    int first = tip_beam ? root : 0;
    int last = tip_beam ? kSections - 1 : segments;
    for (int iteration = 0; iteration < steps; ++iteration) {
      spine_positions(base_goal_, tip_goal_, sx, sy, bx, by, px, py);
      const std::array<double, kSections>& cx = tip_beam ? px : bx;
      const std::array<double, kSections>& cy = tip_beam ? py : by;
      double rx = goal_x[station] - cx[station];
      double ry = goal_y[station] - cy[station];
      // A = J J^T + damping, J's column k = (-(y_s - y_k), x_s - x_k).
      double a11 = damping;
      double a12 = 0.0;
      double a22 = damping;
      for (int joint = first; joint < last; ++joint) {
        std::size_t j = static_cast<std::size_t>(joint);
        double jx = -(cy[station] - cy[j]);
        double jy = cx[station] - cx[j];
        a11 += jx * jx;
        a12 += jx * jy;
        a22 += jy * jy;
      }
      double determinant = a11 * a22 - a12 * a12;
      if (std::abs(determinant) < 1.0e-12) {
        break;
      }
      double ux = (a22 * rx - a12 * ry) / determinant;
      double uy = (-a12 * rx + a11 * ry) / determinant;
      for (int joint = first; joint < last; ++joint) {
        std::size_t j = static_cast<std::size_t>(joint);
        double jx = -(cy[station] - cy[j]);
        double jy = cx[station] - cx[j];
        double turn = jx * ux + jy * uy;
        if (tip_beam) {
          tip_goal_[j] += turn;
        } else {
          base_goal_[j] += turn;
        }
      }
    }
  }
}

// A shape's spine segment lengths (its spine as in spine_bends_of).
void Tongue::spine_lengths_of(const TongueContour& shape, std::array<double, kSections>& lengths) const {
  std::array<double, kSections> sx{};
  std::array<double, kSections> sy{};
  for (int station = 0; station < kSections; ++station) {
    std::size_t i = static_cast<std::size_t>(station);
    std::size_t before = i > 0 ? i - 1 : i;
    std::size_t after = i + 1 < static_cast<std::size_t>(kSections) ? i + 1 : i;
    double tx = shape.x[after] - shape.x[before];
    double ty = shape.y[after] - shape.y[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    double depth = spine_up(station);
    sx[i] = shape.x[i] - depth * -ty / norm;
    sy[i] = shape.y[i] - depth * tx / norm;
  }
  for (int segment = 0; segment + 1 < kSections; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    lengths[g] = std::hypot(sx[g + 1] - sx[g], sy[g + 1] - sy[g]);
  }
}

void Tongue::spine_reset() {
  // The clamp's heading and the segment lengths from the rest shape's spine.
  std::array<double, kSections> sx{};
  std::array<double, kSections> sy{};
  for (int station = 0; station < kSections; ++station) {
    std::size_t i = static_cast<std::size_t>(station);
    std::size_t before = i > 0 ? i - 1 : i;
    std::size_t after = i + 1 < static_cast<std::size_t>(kSections) ? i + 1 : i;
    double tx = rest_.x[after] - rest_.x[before];
    double ty = rest_.y[after] - rest_.y[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    double depth = spine_up(station);
    sx[i] = rest_.x[i] - depth * -ty / norm;
    sy[i] = rest_.y[i] - depth * tx / norm;
  }
  for (int segment = 0; segment + 1 < kSections; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    spine_length_[g] = std::hypot(sx[g + 1] - sx[g], sy[g + 1] - sy[g]);
  }
  spine_rest_length_ = spine_length_;
  spine_goal_length_ = spine_length_;
  clamp_heading_ = std::atan2(sy[1] - sy[0], sx[1] - sx[0]);
  anchor_dx_ = sx[0] - rest_.x[0];
  anchor_dy_ = sy[0] - rest_.y[0];
  base_rest_bend_.fill(0.0);
  tip_rest_bend_.fill(0.0);
  spine_bends_of(rest_, base_rest_bend_, tip_rest_bend_);
  base_bend_ = base_rest_bend_;
  tip_bend_ = tip_rest_bend_;
  base_speed_.fill(0.0);
  tip_speed_.fill(0.0);
  base_goal_ = base_rest_bend_;
  tip_goal_ = tip_rest_bend_;
  root_x_ = rest_.x[0];
  root_y_ = rest_.y[0];
  spine_positions(base_bend_, tip_bend_, spine_x_, spine_y_, base_x_, base_y_, tip_x_, tip_y_);
  base_rest_y_ = base_y_;
}

void Tongue::spine_place(const TongueContour& shape) {
  spine_bends_of(shape, base_bend_, tip_bend_);
  spine_lengths_of(shape, spine_goal_length_);
  for (int segment = 0; segment + 1 < kSections; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    spine_length_[g] = std::clamp(spine_goal_length_[g], spine_rest_length_[g] * (1.0 - settings_.spine_stretch),
                                  spine_rest_length_[g] * (1.0 + settings_.spine_stretch));
  }
  base_speed_.fill(0.0);
  tip_speed_.fill(0.0);
}

// The surface around the spine: the upper surface at spine_up along the
// spine's normals (its first point the root base itself), the tip's cap a
// half circle about the spine's end, the free blade's underside at
// spine_down, the frenulum straight down onto the floor.
void Tongue::spine_surface() {
  spine_positions(base_bend_, tip_bend_, spine_x_, spine_y_, base_x_, base_y_, tip_x_, tip_y_);
  std::array<double, kSections> nx{};
  std::array<double, kSections> ny{};
  for (int station = 0; station < kSections; ++station) {
    std::size_t i = static_cast<std::size_t>(station);
    std::size_t before = i > 0 ? i - 1 : i;
    std::size_t after = i + 1 < static_cast<std::size_t>(kSections) ? i + 1 : i;
    double tx = spine_x_[after] - spine_x_[before];
    double ty = spine_y_[after] - spine_y_[before];
    double norm = std::max(std::hypot(tx, ty), 1.0e-9);
    nx[i] = -ty / norm;
    ny[i] = tx / norm;
    double depth = spine_up(station);
    contour_.x[i] = spine_x_[i] + depth * nx[i];
    contour_.y[i] = spine_y_[i] + depth * ny[i];
  }
  contour_.x[0] = root_x_;
  contour_.y[0] = root_y_;
  contour_.count = kSections;
  constexpr int kCap = 3;
  constexpr int kUnder = 4;
  std::size_t end = static_cast<std::size_t>(kSections - 1);
  double radius = spine_up(kSections - 1);
  double tx = ny[end];
  double ty = -nx[end];
  underside_.count = kCap + kUnder + 1;
  for (int point = 0; point < kCap; ++point) {
    double angle = kPi * static_cast<double>(point + 1) / static_cast<double>(kCap + 1);
    underside_.x[static_cast<std::size_t>(point)] = spine_x_[end] + radius * (std::cos(angle) * nx[end] + std::sin(angle) * tx);
    underside_.y[static_cast<std::size_t>(point)] = spine_y_[end] + radius * (std::cos(angle) * ny[end] + std::sin(angle) * ty);
  }
  for (int under = 0; under < kUnder; ++under) {
    std::size_t i = end - static_cast<std::size_t>(under);
    double depth = spine_down(kSections - 1 - under);
    std::size_t o = static_cast<std::size_t>(kCap + under);
    underside_.x[o] = spine_x_[i] - depth * nx[i];
    underside_.y[o] = spine_y_[i] - depth * ny[i];
  }
  std::size_t back = static_cast<std::size_t>(kCap + kUnder - 1);
  underside_.x[back + 1] = underside_.x[back];
  underside_.y[back + 1] = std::min(underside_.y[back], floor_height(underside_.x[back]));
  rest_on_floor();
}

// The solid hull against the upper surface and the cap: a point inside is
// returned to its nearest hull point by turning the joints behind it (the
// least turn that moves it so, over the beams that carry it), and those
// joints lose half their speed.
void Tongue::spine_contact() {
  int segments = spine_base_segments();
  int root = spine_tip_root();
  for (int pass = 0; pass < 3; ++pass) {
    bool moved = false;
    for (int point = 6; point < kSections + 3; ++point) {
      bool cap = point >= kSections;
      int station = cap ? kSections - 1 : point;
      double px = cap ? underside_.x[static_cast<std::size_t>(point - kSections)] : contour_.x[static_cast<std::size_t>(point)];
      double py = cap ? underside_.y[static_cast<std::size_t>(point - kSections)] : contour_.y[static_cast<std::size_t>(point)];
      if (!inside_hull(px, py)) {
        continue;
      }
      double qx = 0.0;
      double qy = 0.0;
      double hx = 0.0;
      double hy = 0.0;
      nearest_hull_point(px, py, qx, qy, hx, hy);
      double dx = qx - px;
      double dy = qy - py;
      // The joints that carry the station: the base beam's up to the station
      // (or its end), weighted 1 - share in the overlap; the tip beam's from
      // its root, weighted by the share.
      double share = station <= root ? 0.0 : station <= segments ? static_cast<double>(station - root) / static_cast<double>(segments - root) : 1.0;
      double norm = 0.0;
      for (int pass_beam = 0; pass_beam < 2; ++pass_beam) {
        bool tip_beam = pass_beam == 1;
        double weight = tip_beam ? share : 1.0 - share;
        int first = tip_beam ? root : 0;
        int last = tip_beam ? station : std::min(station, segments);
        for (int joint = first; weight > 0.0 && joint < last; ++joint) {
          std::size_t j = static_cast<std::size_t>(joint);
          double jx = tip_beam ? tip_x_[j] : base_x_[j];
          double jy = tip_beam ? tip_y_[j] : base_y_[j];
          double lx = -(py - jy);
          double ly = px - jx;
          norm += weight * weight * (lx * lx + ly * ly);
        }
      }
      if (norm <= 1.0e-9) {
        continue;
      }
      for (int pass_beam = 0; pass_beam < 2; ++pass_beam) {
        bool tip_beam = pass_beam == 1;
        double weight = tip_beam ? share : 1.0 - share;
        int first = tip_beam ? root : 0;
        int last = tip_beam ? station : std::min(station, segments);
        for (int joint = first; weight > 0.0 && joint < last; ++joint) {
          std::size_t j = static_cast<std::size_t>(joint);
          double jx = tip_beam ? tip_x_[j] : base_x_[j];
          double jy = tip_beam ? tip_y_[j] : base_y_[j];
          double turn = weight * (-(py - jy) * dx + (px - jx) * dy) / norm;
          if (tip_beam) {
            tip_bend_[j] += turn;
            tip_speed_[j] *= 0.5;
          } else {
            base_bend_[j] += turn;
            base_speed_[j] *= 0.5;
          }
        }
      }
      moved = true;
      spine_surface();
    }
    if (!moved) {
      break;
    }
  }
}

// One step of the spine tongue: each bend driven toward the shape's
// (base_goal_, tip_goal_) as a damped second-order system; the base beam's
// butt pulled back toward its rest bends (spine_butt, fading to nothing at
// its end) and its forward end held from rising more than spine_rise.
void Tongue::spine_step(double seconds) {
  spine_bends_of(target_, base_goal_, tip_goal_);
  spine_lengths_of(target_, spine_goal_length_);
  double follow = 1.0 - std::exp(-seconds / std::max(settings_.spine_stretch_time, 1.0e-4));
  for (int segment = 0; segment + 1 < kSections; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    double wanted = std::clamp(spine_goal_length_[g], spine_rest_length_[g] * (1.0 - settings_.spine_stretch),
                               spine_rest_length_[g] * (1.0 + settings_.spine_stretch));
    spine_length_[g] += follow * (wanted - spine_length_[g]);
  }
  if (settings_.spine_reach > 0.5) {
    std::array<double, kSections> goal_x{};
    std::array<double, kSections> goal_y{};
    spine_points_of(target_, goal_x, goal_y);
    spine_reach_goals(goal_x, goal_y);
  }
  int segments = spine_base_segments();
  int root = spine_tip_root();
  double h = seconds / static_cast<double>(kSubsteps);
  double base_omega = 2.0 * kPi * settings_.spine_base_hz;
  double tip_omega = 2.0 * kPi * settings_.spine_tip_hz;
  for (int substep = 0; substep < kSubsteps; ++substep) {
    spine_positions(base_bend_, tip_bend_, spine_x_, spine_y_, base_x_, base_y_, tip_x_, tip_y_);
    std::size_t front = static_cast<std::size_t>(segments);
    double rise = base_y_[front] - base_rest_y_[front] - settings_.spine_rise;
    for (int joint = 0; joint < segments; ++joint) {
      std::size_t j = static_cast<std::size_t>(joint);
      double butt = settings_.spine_butt * (1.0 - static_cast<double>(joint) / static_cast<double>(segments));
      double pull = (base_goal_[j] - base_bend_[j]) + butt * (base_rest_bend_[j] - base_bend_[j]);
      double accel = base_omega * base_omega * pull - 2.0 * settings_.spine_damping * base_omega * base_speed_[j];
      if (rise > 0.0) {
        // Turning joint j raises the forward end by its lever (the end's x
        // less the joint's); the limit turns it back in proportion.
        double lever = base_x_[front] - base_x_[j];
        accel -= 4.0 * base_omega * base_omega * rise * lever;
      }
      base_speed_[j] += h * accel;
      base_bend_[j] += h * base_speed_[j];
    }
    for (int joint = root; joint + 1 < kSections; ++joint) {
      std::size_t j = static_cast<std::size_t>(joint);
      double accel = tip_omega * tip_omega * (tip_goal_[j] - tip_bend_[j]) - 2.0 * settings_.spine_damping * tip_omega * tip_speed_[j];
      tip_speed_[j] += h * accel;
      tip_bend_[j] += h * tip_speed_[j];
    }
  }
  TongueContour before = contour_;
  spine_surface();
  spine_contact();
  for (int station = 0; station < kSections; ++station) {
    std::size_t i = static_cast<std::size_t>(station);
    vx_[i] = (contour_.x[i] - before.x[i]) / seconds;
    vy_[i] = (contour_.y[i] - before.y[i]) / seconds;
  }
}

void Tongue::step(double seconds, const Articulators& frame, const MidsagittalTract& anatomy) {
  prepare_frame(frame, anatomy);
  if (settings_.spine > 0.5) {
    spine_step(seconds);
    measure_penetration();
    return;
  }
  if (settings_.skeleton > 0.5) {
    skeleton_step(seconds);
    measure_penetration();
    return;
  }
  if (settings_.lateral_spread > 0.0 || settings_.lateral_narrow > 0.0) {
    // The area the plan asks for: of the motors' smooth posture when the
    // roles hold toward it (the look-ahead target steps as events enter it).
    double planned = settings_.lateral_from_posture > 0.5 && geometric_set_ ? area_of(geometric_) : planned_area_;
    double wanted = std::clamp(planned, rest_area_ * (1.0 - settings_.lateral_spread),
                               rest_area_ * (1.0 + settings_.lateral_narrow));
    held_area_ += (1.0 - std::exp(-seconds / std::max(settings_.lateral_time, 1.0e-4))) * (wanted - held_area_);
  }
  if (settings_.feedback_aim > 0.5) {
    build_aim();
  }
  guarding_ = settings_.guard > 0.5 && geometric_set_;
  if (guarding_) {
    build_guard();
  }
  if (settings_.statics > 0.5) {
    solve_statics();
  } else if (!measured_) {
    measure_influences();
  }
  double h = seconds / static_cast<double>(kSubsteps);
  for (int substep = 0; substep < kSubsteps; ++substep) {
    drive_muscles(h, origins_);
    advance(h, true);
  }
  if (guarding_) {
    measure_intrusion();
  }
  measure_penetration();
}

}  // namespace tract
