// The skeleton tongue: two internal supports, the tongue skinned around
// them (TongueSettings::skeleton; tongue.hpp).

#include "tract/tongue.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

constexpr int kResiduals = 2 * kTongueContourPoints + Tongue::kBoneParameters + 2 + 2 * (Tongue::kBoneSegments - 1);
constexpr std::size_t kAnchorX = static_cast<std::size_t>(2 * (Tongue::kBoneSegments + 1));
constexpr std::size_t kAnchorY = kAnchorX + 1;

double smooth_share(double value, double low, double high) {
  double share = std::clamp((value - low) / (high - low), 0.0, 1.0);
  return share * share * (3.0 - 2.0 * share);
}

// Solves the n x n system (row-major, n <= kBoneParameters) by Gaussian
// elimination with partial pivoting; returns false if it is singular.
bool solve_system(std::array<std::array<double, Tongue::kBoneParameters>, Tongue::kBoneParameters>& matrix,
                  std::array<double, Tongue::kBoneParameters>& right, int n) {
  for (int column = 0; column < n; ++column) {
    std::size_t c = static_cast<std::size_t>(column);
    std::size_t pivot = c;
    for (int row = column + 1; row < n; ++row) {
      if (std::abs(matrix[static_cast<std::size_t>(row)][c]) > std::abs(matrix[pivot][c])) {
        pivot = static_cast<std::size_t>(row);
      }
    }
    if (std::abs(matrix[pivot][c]) < 1.0e-14) {
      return false;
    }
    std::swap(matrix[pivot], matrix[c]);
    std::swap(right[pivot], right[c]);
    for (int row = column + 1; row < n; ++row) {
      std::size_t r = static_cast<std::size_t>(row);
      double factor = matrix[r][c] / matrix[c][c];
      for (int k = column; k < n; ++k) {
        matrix[r][static_cast<std::size_t>(k)] -= factor * matrix[c][static_cast<std::size_t>(k)];
      }
      right[r] -= factor * right[c];
    }
  }
  for (int row = n - 1; row >= 0; --row) {
    std::size_t r = static_cast<std::size_t>(row);
    double sum = right[r];
    for (int k = row + 1; k < n; ++k) {
      sum -= matrix[r][static_cast<std::size_t>(k)] * right[static_cast<std::size_t>(k)];
    }
    right[r] = sum / matrix[r][r];
  }
  return true;
}

// A point at share `share` (0..1) along a support's points.
void along_support(const std::array<double, Tongue::kBoneSegments + 1>& xs,
                   const std::array<double, Tongue::kBoneSegments + 1>& ys, double share, double& x, double& y) {
  double place = std::clamp(share, 0.0, 1.0) * static_cast<double>(Tongue::kBoneSegments);
  std::size_t segment = std::min(static_cast<std::size_t>(place), static_cast<std::size_t>(Tongue::kBoneSegments - 1));
  double fraction = place - static_cast<double>(segment);
  x = xs[segment] + fraction * (xs[segment + 1] - xs[segment]);
  y = ys[segment] + fraction * (ys[segment + 1] - ys[segment]);
}

// The unit tangent of a polyline at segment `segment`, share `share`: the
// segments' directions blended toward each neighbour across each vertex
// (continuous along the curve).
void curve_tangent(const Tongue::CurvePoints& cx, const Tongue::CurvePoints& cy, int segment, double share,
                   double& tx, double& ty) {
  std::size_t g = static_cast<std::size_t>(segment);
  double dx = cx[g + 1] - cx[g];
  double dy = cy[g + 1] - cy[g];
  double length = std::max(std::hypot(dx, dy), 1.0e-9);
  double ux = dx / length;
  double uy = dy / length;
  std::size_t neighbour = share < 0.5 ? (g > 0 ? g - 1 : g) : (g + 2 < static_cast<std::size_t>(Tongue::kCurvePoints) ? g + 1 : g);
  double nx = cx[neighbour + 1] - cx[neighbour];
  double ny = cy[neighbour + 1] - cy[neighbour];
  double other = std::max(std::hypot(nx, ny), 1.0e-9);
  double weight = share < 0.5 ? 0.5 - share : share - 0.5;
  tx = (1.0 - weight) * ux + weight * nx / other;
  ty = (1.0 - weight) * uy + weight * ny / other;
  double norm = std::max(std::hypot(tx, ty), 1.0e-9);
  tx /= norm;
  ty /= norm;
}

}  // namespace

// The central curve: the base support to the tip support's root (15
// points), the overlap where the base support's last quarter fades into the
// tip support's start (10 points, smoothstep weights), then the tip support
// to its end.
void Tongue::skeleton_curve_of(const BonePoints& xs, const BonePoints& ys, CurvePoints& cx, CurvePoints& cy) const {
  constexpr int kBase = 15;
  constexpr int kOverlap = 10;
  constexpr int kTip = kCurvePoints - kBase - kOverlap;
  double root = tip_root_share_;
  // The overlap's length on the tip support: as long as the base support's last quarter.
  double tip_overlap = std::clamp((1.0 - root) * bone_rest_length_[0] / std::max(bone_rest_length_[1], 1.0e-9), 0.05, 0.6);
  int index = 0;
  for (int point = 0; point < kBase; ++point) {
    double share = root * static_cast<double>(point) / static_cast<double>(kBase);
    along_support(xs[0], ys[0], share, cx[static_cast<std::size_t>(index)], cy[static_cast<std::size_t>(index)]);
    ++index;
  }
  for (int point = 0; point < kOverlap; ++point) {
    double u = static_cast<double>(point) / static_cast<double>(kOverlap);
    double weight = smooth_share(u, 0.0, 1.0);
    double bx = 0.0;
    double by = 0.0;
    double tx = 0.0;
    double ty = 0.0;
    along_support(xs[0], ys[0], root + u * (1.0 - root), bx, by);
    along_support(xs[1], ys[1], u * tip_overlap, tx, ty);
    cx[static_cast<std::size_t>(index)] = (1.0 - weight) * bx + weight * tx;
    cy[static_cast<std::size_t>(index)] = (1.0 - weight) * by + weight * ty;
    ++index;
  }
  for (int point = 0; point < kTip; ++point) {
    double share = tip_overlap + (1.0 - tip_overlap) * static_cast<double>(point) / static_cast<double>(kTip - 1);
    along_support(xs[1], ys[1], share, cx[static_cast<std::size_t>(index)], cy[static_cast<std::size_t>(index)]);
    ++index;
  }
}

// The rest supports and every skin point's attachment. The base support
// runs from a third of the way from the root base to the body's centre, to
// past the centre toward the tip; the tip support from the base support's
// last quarter to the tip cap's centre. A skin point attaches to each
// support at its nearest point there, by its offset along and across that
// segment.
void Tongue::skeleton_reset() {
  TongueUnderside under;
  derive_underside(rest_, under);
  skin_count_ = kSections + under.count - 1;  // the frenulum is not skinned
  std::array<double, kSkinPoints> sx{};
  std::array<double, kSkinPoints> sy{};
  double centre_x = 0.0;
  double centre_y = 0.0;
  for (int point = 0; point < skin_count_; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    bool upper = point < kSections;
    sx[p] = upper ? rest_.x[p] : under.x[static_cast<std::size_t>(point - kSections)];
    sy[p] = upper ? rest_.y[p] : under.y[static_cast<std::size_t>(point - kSections)];
    centre_x += sx[p];
    centre_y += sy[p];
  }
  centre_x /= static_cast<double>(skin_count_);
  centre_y /= static_cast<double>(skin_count_);
  // The tip cap's centre: tip_radius below the tip along the blade's normal.
  std::size_t tip = static_cast<std::size_t>(kSections - 1);
  std::size_t blade = tip - 3;
  double bx = rest_.x[tip] - rest_.x[blade];
  double by = rest_.y[tip] - rest_.y[blade];
  double norm = std::max(std::hypot(bx, by), 1.0e-9);
  double tip_x = rest_.x[tip] - settings_.tip_radius * -by / norm;
  double tip_y = rest_.y[tip] - settings_.tip_radius * bx / norm;
  double anchor_x = rest_.x[0] + 0.35 * (centre_x - rest_.x[0]);
  double anchor_y = rest_.y[0] + 0.35 * (centre_y - rest_.y[0]);
  double end_x = centre_x + 0.3 * (tip_x - centre_x);
  double end_y = centre_y + 0.3 * (tip_y - centre_y);
  bone_anchor_dx_ = anchor_x - rest_.x[0];
  bone_anchor_dy_ = anchor_y - rest_.y[0];
  double root_x = anchor_x + tip_root_share_ * (end_x - anchor_x);
  double root_y = anchor_y + tip_root_share_ * (end_y - anchor_y);
  double starts_x[2] = {anchor_x, root_x};
  double starts_y[2] = {anchor_y, root_y};
  double ends_x[2] = {end_x, tip_x};
  double ends_y[2] = {end_y, tip_y};
  for (int bone = 0; bone < 2; ++bone) {
    std::size_t b = static_cast<std::size_t>(bone);
    double heading = std::atan2(ends_y[b] - starts_y[b], ends_x[b] - starts_x[b]);
    bone_rest_length_[b] = std::hypot(ends_x[b] - starts_x[b], ends_y[b] - starts_y[b]);
    for (int point = 0; point <= kBoneSegments; ++point) {
      double share = static_cast<double>(point) / static_cast<double>(kBoneSegments);
      bone_rest_x_[b][static_cast<std::size_t>(point)] = starts_x[b] + share * (ends_x[b] - starts_x[b]);
      bone_rest_y_[b][static_cast<std::size_t>(point)] = starts_y[b] + share * (ends_y[b] - starts_y[b]);
    }
    for (int segment = 0; segment < kBoneSegments; ++segment) {
      bone_rest_heading_[b][static_cast<std::size_t>(segment)] = heading;
    }
    // Attachments: the nearest point of the support to each skin point.
    for (int point = 0; point < skin_count_; ++point) {
      std::size_t p = static_cast<std::size_t>(point);
      double best = 1.0e9;
      for (int segment = 0; segment < kBoneSegments; ++segment) {
        std::size_t g = static_cast<std::size_t>(segment);
        double ax = bone_rest_x_[b][g];
        double ay = bone_rest_y_[b][g];
        double dx = bone_rest_x_[b][g + 1] - ax;
        double dy = bone_rest_y_[b][g + 1] - ay;
        double span = dx * dx + dy * dy;
        double share = std::clamp(((sx[p] - ax) * dx + (sy[p] - ay) * dy) / span, 0.0, 1.0);
        double qx = ax + share * dx;
        double qy = ay + share * dy;
        double gap = std::hypot(sx[p] - qx, sy[p] - qy);
        if (gap < best) {
          best = gap;
          double length = std::sqrt(span);
          double tx = dx / length;
          double ty = dy / length;
          skin_segment_[b][p] = segment;
          skin_share_[b][p] = share;
          skin_along_[b][p] = (sx[p] - qx) * tx + (sy[p] - qy) * ty;
          skin_across_[b][p] = (sx[p] - qx) * -ty + (sy[p] - qy) * tx;
        }
      }
    }
  }
  for (int point = 0; point < skin_count_; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    skin_tip_weight_[p] = point >= kSections ? 1.0 : smooth_share(static_cast<double>(point), 11.0, 17.0);
  }
  bone_parameters_.fill(0.0);
  bone_goal_.fill(0.0);
  bone_speed_.fill(0.0);
  root_x_ = rest_.x[0];
  root_y_ = rest_.y[0];
  skeleton_pose(bone_parameters_, bone_x_, bone_y_);
  // Attachments to the rest central curve: its nearest point to each skin
  // point, and the offsets along and across it there.
  CurvePoints cx{};
  CurvePoints cy{};
  skeleton_curve_of(bone_rest_x_, bone_rest_y_, cx, cy);
  for (int point = 0; point < skin_count_; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    double best = 1.0e9;
    for (int segment = 0; segment + 1 < kCurvePoints; ++segment) {
      std::size_t g = static_cast<std::size_t>(segment);
      double dx = cx[g + 1] - cx[g];
      double dy = cy[g + 1] - cy[g];
      double span = std::max(dx * dx + dy * dy, 1.0e-12);
      double share = std::clamp(((sx[p] - cx[g]) * dx + (sy[p] - cy[g]) * dy) / span, 0.0, 1.0);
      double qx = cx[g] + share * dx;
      double qy = cy[g] + share * dy;
      double gap = std::hypot(sx[p] - qx, sy[p] - qy);
      if (gap < best) {
        best = gap;
        double tx = 0.0;
        double ty = 0.0;
        curve_tangent(cx, cy, segment, share, tx, ty);
        curve_segment_[p] = segment;
        curve_share_[p] = share;
        curve_along_[p] = (sx[p] - qx) * tx + (sy[p] - qy) * ty;
        curve_across_[p] = (sx[p] - qx) * -ty + (sy[p] - qy) * tx;
      }
    }
  }
}

// The supports' points for a set of parameters: each support's segments
// turn by its bends (cumulative from its start) and lengthen by its
// stretch; the base support starts at the anchor (with the root base), the
// tip support at its root on the base support, turned with it.
void Tongue::skeleton_pose(const BoneParameters& parameters, BonePoints& xs, BonePoints& ys) const {
  constexpr std::size_t kStride = static_cast<std::size_t>(kBoneSegments + 1);
  double move = settings_.skeleton_anchor_move;
  xs[0][0] = root_x_ + bone_anchor_dx_ + (move > 0.0 ? std::clamp(parameters[kAnchorX], -move, move) : 0.0);
  ys[0][0] = root_y_ + bone_anchor_dy_ + (move > 0.0 ? std::clamp(parameters[kAnchorY], -move, move) : 0.0);
  double base_stretch = std::exp(std::clamp(parameters[kBoneSegments], -settings_.skeleton_stretch, settings_.skeleton_stretch));
  double base_length = bone_rest_length_[0] * base_stretch / static_cast<double>(kBoneSegments);
  std::array<double, kBoneSegments> turned{};
  double turn = 0.0;
  for (int segment = 0; segment < kBoneSegments; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    turn += parameters[g];
    turned[g] = turn;
    double heading = bone_rest_heading_[0][g] + turn;
    xs[0][g + 1] = xs[0][g] + base_length * std::cos(heading);
    ys[0][g + 1] = ys[0][g] + base_length * std::sin(heading);
  }
  double place = tip_root_share_ * static_cast<double>(kBoneSegments);
  std::size_t root_segment = std::min(static_cast<std::size_t>(place), static_cast<std::size_t>(kBoneSegments - 1));
  double share = place - static_cast<double>(root_segment);
  xs[1][0] = xs[0][root_segment] + share * (xs[0][root_segment + 1] - xs[0][root_segment]);
  ys[1][0] = ys[0][root_segment] + share * (ys[0][root_segment + 1] - ys[0][root_segment]);
  double tip_stretch = std::exp(std::clamp(parameters[2 * kStride - 1], -settings_.skeleton_stretch, settings_.skeleton_stretch));
  double tip_length = bone_rest_length_[1] * tip_stretch / static_cast<double>(kBoneSegments);
  turn = turned[root_segment];
  for (int segment = 0; segment < kBoneSegments; ++segment) {
    std::size_t g = static_cast<std::size_t>(segment);
    turn += parameters[kStride + g];
    double heading = bone_rest_heading_[1][g] + turn;
    xs[1][g + 1] = xs[1][g] + tip_length * std::cos(heading);
    ys[1][g + 1] = ys[1][g] + tip_length * std::sin(heading);
  }
}

// The skin: each point on each support by its attachment (the posed
// segment's point at its share, plus its offsets along and across the
// segment), blended by its share on the tip support; the first point is
// the root base itself.
void Tongue::skeleton_skin(const BoneParameters& parameters, std::array<double, kSkinPoints>& xs,
                           std::array<double, kSkinPoints>& ys) const {
  BonePoints bx{};
  BonePoints by{};
  skeleton_pose(parameters, bx, by);
  if (settings_.skeleton_curve > 0.5) {
    CurvePoints cx{};
    CurvePoints cy{};
    skeleton_curve_of(bx, by, cx, cy);
    for (int point = 0; point < skin_count_; ++point) {
      std::size_t p = static_cast<std::size_t>(point);
      std::size_t g = static_cast<std::size_t>(curve_segment_[p]);
      double qx = cx[g] + curve_share_[p] * (cx[g + 1] - cx[g]);
      double qy = cy[g] + curve_share_[p] * (cy[g + 1] - cy[g]);
      double tx = 0.0;
      double ty = 0.0;
      curve_tangent(cx, cy, curve_segment_[p], curve_share_[p], tx, ty);
      xs[p] = qx + curve_along_[p] * tx + curve_across_[p] * -ty;
      ys[p] = qy + curve_along_[p] * ty + curve_across_[p] * tx;
    }
    xs[0] = root_x_;
    ys[0] = root_y_;
    return;
  }
  for (int point = 0; point < skin_count_; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    double on_x[2] = {0.0, 0.0};
    double on_y[2] = {0.0, 0.0};
    for (int bone = 0; bone < 2; ++bone) {
      std::size_t b = static_cast<std::size_t>(bone);
      std::size_t g = static_cast<std::size_t>(skin_segment_[b][p]);
      double dx = bx[b][g + 1] - bx[b][g];
      double dy = by[b][g + 1] - by[b][g];
      double length = std::max(std::hypot(dx, dy), 1.0e-9);
      double tx = dx / length;
      double ty = dy / length;
      double qx = bx[b][g] + skin_share_[b][p] * dx;
      double qy = by[b][g] + skin_share_[b][p] * dy;
      on_x[b] = qx + skin_along_[b][p] * tx + skin_across_[b][p] * -ty;
      on_y[b] = qy + skin_along_[b][p] * ty + skin_across_[b][p] * tx;
    }
    double weight = skin_tip_weight_[p];
    xs[p] = (1.0 - weight) * on_x[0] + weight * on_x[1];
    ys[p] = (1.0 - weight) * on_y[0] + weight * on_y[1];
  }
  xs[0] = root_x_;
  ys[0] = root_y_;
}

// A shape's points as the fit's goals: those inside the hull held on its
// nearest point (a closure presses to the hull, not into it).
void Tongue::skeleton_goals_of(const TongueContour& shape, std::array<double, kSections>& goal_x,
                               std::array<double, kSections>& goal_y) const {
  for (int point = 0; point < kSections; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    goal_x[p] = shape.x[p];
    goal_y[p] = shape.y[p];
    if (inside_hull(shape.x[p], shape.y[p])) {
      double nx = 0.0;
      double ny = 0.0;
      nearest_hull_point(shape.x[p], shape.y[p], goal_x[p], goal_y[p], nx, ny);
    }
  }
}

// Gauss-Newton on the parameters (forward-difference Jacobian): the skin's
// upper points toward the goals (each by its weight), every parameter
// pulled toward rest (skeleton_stiffness), and the base support's forward
// end held below its rest height plus skeleton_rise.
void Tongue::skeleton_fit(const std::array<double, kSections>& goal_x, const std::array<double, kSections>& goal_y,
                          const std::array<double, kSections>& weight, BoneParameters& parameters, int iterations,
                          const BoneParameters* around) const {
  // With `around` (minimal contact): each parameter held to its value there
  // instead of rest; no smoothness, rise or bend limit; the anchor fixed.
  constexpr double kHold = 0.2;
  constexpr double kStep = 1.0e-4;
  constexpr double kLevenberg = 1.0e-3;
  std::array<double, kResiduals> residual{};
  std::array<double, kResiduals> nudged{};
  std::array<std::array<double, kResiduals>, kBoneParameters> jacobian{};
  std::array<double, kSkinPoints> xs{};
  std::array<double, kSkinPoints> ys{};
  BonePoints bx{};
  BonePoints by{};
  double rest_end_y = bone_rest_y_[0][static_cast<std::size_t>(kBoneSegments)] - (bone_rest_y_[0][0] - (root_y_ + bone_anchor_dy_));
  for (int iteration = 0; iteration < iterations; ++iteration) {
    for (int column = -1; column < kBoneParameters; ++column) {
      if (around != nullptr && column >= 0 && static_cast<std::size_t>(column) >= kAnchorX) {
        for (int k = 0; k < kResiduals; ++k) {
          jacobian[static_cast<std::size_t>(column)][static_cast<std::size_t>(k)] = 0.0;
        }
        continue;
      }
      BoneParameters trial = parameters;
      if (column >= 0) {
        trial[static_cast<std::size_t>(column)] += kStep;
      }
      skeleton_skin(trial, xs, ys);
      skeleton_pose(trial, bx, by);
      std::array<double, kResiduals>& out = column < 0 ? residual : nudged;
      std::size_t r = 0;
      for (int point = 1; point < kSections; ++point) {
        std::size_t p = static_cast<std::size_t>(point);
        double root_weight = std::sqrt(weight[p]);
        out[r++] = root_weight * (xs[p] - goal_x[p]);
        out[r++] = root_weight * (ys[p] - goal_y[p]);
      }
      for (int k = 0; k < kBoneParameters; ++k) {
        std::size_t i = static_cast<std::size_t>(k);
        if (around != nullptr) {
          out[r++] = kHold * (trial[i] - (*around)[i]);
          continue;
        }
        double stiffness = i >= kAnchorX ? settings_.skeleton_anchor_stiffness : settings_.skeleton_stiffness;
        out[r++] = stiffness * trial[i];
      }
      for (int support = 0; around == nullptr && settings_.skeleton_smoothness > 0.0 && support < 2; ++support) {
        std::size_t first = static_cast<std::size_t>(support * (kBoneSegments + 1));
        for (int k = 1; k < kBoneSegments; ++k) {
          std::size_t i = first + static_cast<std::size_t>(k);
          out[r++] = settings_.skeleton_smoothness * (trial[i] - trial[i - 1]);
        }
      }
      double rise = by[0][static_cast<std::size_t>(kBoneSegments)] - rest_end_y - settings_.skeleton_rise;
      out[r++] = rise > 0.0 && around == nullptr ? 3.0 * rise : 0.0;
      if (settings_.skeleton_bend_limit > 0.0 && around == nullptr) {
        double turn = 0.0;
        for (int k = 0; k < kBoneSegments; ++k) {
          turn += trial[static_cast<std::size_t>(k)];
        }
        double beyond = std::abs(turn) - settings_.skeleton_bend_limit;
        out[r++] = beyond > 0.0 ? 3.0 * beyond : 0.0;
      }
      while (r < static_cast<std::size_t>(kResiduals)) {
        out[r++] = 0.0;
      }
      if (column >= 0) {
        for (int k = 0; k < kResiduals; ++k) {
          std::size_t i = static_cast<std::size_t>(k);
          jacobian[static_cast<std::size_t>(column)][i] = (nudged[i] - residual[i]) / kStep;
        }
      }
    }
    std::array<std::array<double, kBoneParameters>, kBoneParameters> normal{};
    std::array<double, kBoneParameters> right{};
    for (int a = 0; a < kBoneParameters; ++a) {
      std::size_t ia = static_cast<std::size_t>(a);
      for (int b = 0; b < kBoneParameters; ++b) {
        std::size_t ib = static_cast<std::size_t>(b);
        double sum = 0.0;
        for (int k = 0; k < kResiduals; ++k) {
          sum += jacobian[ia][static_cast<std::size_t>(k)] * jacobian[ib][static_cast<std::size_t>(k)];
        }
        normal[ia][ib] = sum + (a == b ? kLevenberg : 0.0);
      }
      double sum = 0.0;
      for (int k = 0; k < kResiduals; ++k) {
        sum += jacobian[ia][static_cast<std::size_t>(k)] * residual[static_cast<std::size_t>(k)];
      }
      right[ia] = -sum;
    }
    if (!solve_system(normal, right, kBoneParameters)) {
      return;
    }
    double limit = around != nullptr ? 0.1 : 0.3;
    for (int k = 0; k < kBoneParameters; ++k) {
      std::size_t i = static_cast<std::size_t>(k);
      parameters[i] += std::clamp(right[i], -limit, limit);
    }
    parameters[kBoneSegments] = std::clamp(parameters[kBoneSegments], -settings_.skeleton_stretch, settings_.skeleton_stretch);
    parameters[kAnchorX - 1] = std::clamp(parameters[kAnchorX - 1], -settings_.skeleton_stretch, settings_.skeleton_stretch);
  }
}

// The fit's cost of a set of parameters (its goal and stiffness terms).
double Tongue::skeleton_cost(const std::array<double, kSections>& goal_x, const std::array<double, kSections>& goal_y,
                             const std::array<double, kSections>& weight, const BoneParameters& parameters) const {
  std::array<double, kSkinPoints> xs{};
  std::array<double, kSkinPoints> ys{};
  skeleton_skin(parameters, xs, ys);
  double sum = 0.0;
  for (int point = 1; point < kSections; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    sum += weight[p] * ((xs[p] - goal_x[p]) * (xs[p] - goal_x[p]) + (ys[p] - goal_y[p]) * (ys[p] - goal_y[p]));
  }
  for (int k = 0; k < kBoneParameters; ++k) {
    std::size_t i = static_cast<std::size_t>(k);
    double stiffness = i >= kAnchorX ? settings_.skeleton_anchor_stiffness : settings_.skeleton_stiffness;
    sum += stiffness * stiffness * parameters[i] * parameters[i];
  }
  return sum;
}

void Tongue::skeleton_place(const TongueContour& shape) {
  std::array<double, kSections> goal_x{};
  std::array<double, kSections> goal_y{};
  std::array<double, kSections> weight{};
  weight.fill(1.0);
  skeleton_goals_of(shape, goal_x, goal_y);
  bone_goal_.fill(0.0);
  skeleton_fit(goal_x, goal_y, weight, bone_goal_, 8);
  bone_parameters_ = bone_goal_;
  bone_speed_.fill(0.0);
  skeleton_surface();
}

// The skin into the contour and the underside (its frenulum dropped onto
// the floor), and the supports' points.
void Tongue::skeleton_surface() {
  std::array<double, kSkinPoints> xs{};
  std::array<double, kSkinPoints> ys{};
  skeleton_skin(bone_parameters_, xs, ys);
  skeleton_pose(bone_parameters_, bone_x_, bone_y_);
  for (int point = 0; point < kSections; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    contour_.x[p] = xs[p];
    contour_.y[p] = ys[p];
  }
  contour_.count = kSections;
  int skinned_under = skin_count_ - kSections;
  underside_.count = skinned_under + 1;
  for (int point = 0; point < skinned_under; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    underside_.x[p] = xs[static_cast<std::size_t>(kSections + point)];
    underside_.y[p] = ys[static_cast<std::size_t>(kSections + point)];
  }
  std::size_t back = static_cast<std::size_t>(skinned_under - 1);
  underside_.x[back + 1] = underside_.x[back];
  underside_.y[back + 1] = std::min(underside_.y[back], floor_height(underside_.x[back]));
  rest_on_floor();
}

// The solid hull: upper points (past the root) and the cap inside it are
// sent to their nearest hull points by a Gauss-Newton step on the
// parameters themselves (the rest held where they are), up to twice; the
// speeds halve.
void Tongue::skeleton_contact() {
  for (int pass = 0; pass < 2; ++pass) {
    std::array<double, kSections> goal_x{};
    std::array<double, kSections> goal_y{};
    std::array<double, kSections> weight{};
    bool inside = false;
    for (int point = 0; point < kSections; ++point) {
      std::size_t p = static_cast<std::size_t>(point);
      goal_x[p] = contour_.x[p];
      goal_y[p] = contour_.y[p];
      weight[p] = 0.05;
      if (point >= 4 && inside_hull(contour_.x[p], contour_.y[p])) {
        double nx = 0.0;
        double ny = 0.0;
        nearest_hull_point(contour_.x[p], contour_.y[p], goal_x[p], goal_y[p], nx, ny);
        weight[p] = 1.0;
        inside = true;
      }
    }
    // The cap's front stands for the tip: inside the hull, the tip is sent
    // out by the cap's way out.
    std::size_t tip = static_cast<std::size_t>(kSections - 1);
    for (int point = 0; point < 3 && underside_.count > 3; ++point) {
      std::size_t p = static_cast<std::size_t>(point);
      if (inside_hull(underside_.x[p], underside_.y[p])) {
        double qx = 0.0;
        double qy = 0.0;
        double nx = 0.0;
        double ny = 0.0;
        nearest_hull_point(underside_.x[p], underside_.y[p], qx, qy, nx, ny);
        goal_x[tip] = contour_.x[tip] + (qx - underside_.x[p]);
        goal_y[tip] = contour_.y[tip] + (qy - underside_.y[p]);
        weight[tip] = 1.0;
        inside = true;
      }
    }
    if (!inside) {
      return;
    }
    // The bends answer the contact; the anchor (on the hyoid) is not shoved
    // by it. (Moved by contact, the anchor sprang back and drove the tongue
    // into the palate again: the /n/s bounced every 30 ms and clicked.)
    double anchor_x = bone_parameters_[kAnchorX];
    double anchor_y = bone_parameters_[kAnchorY];
    if (settings_.skeleton_minimal_contact > 0.5) {
      BoneParameters present = bone_parameters_;
      skeleton_fit(goal_x, goal_y, weight, bone_parameters_, 1, &present);
    } else {
      skeleton_fit(goal_x, goal_y, weight, bone_parameters_, 1);
    }
    bone_parameters_[kAnchorX] = anchor_x;
    bone_parameters_[kAnchorY] = anchor_y;
    for (double& speed : bone_speed_) {
      speed *= 0.5;
    }
    skeleton_surface();
  }
}

// One step: the goals refitted to the planner's shape (warm, twice), each
// parameter following its goal as a critically damped second-order system
// (the base support's at skeleton_base_hz, the tip support's at
// skeleton_tip_hz), then the skin and the hull.
void Tongue::skeleton_step(double seconds) {
  std::array<double, kSections> goal_x{};
  std::array<double, kSections> goal_y{};
  std::array<double, kSections> weight{};
  for (int point = 0; point < kSections; ++point) {
    weight[static_cast<std::size_t>(point)] = point >= 18 ? settings_.skeleton_tip_weight : 1.0;
  }
  skeleton_goals_of(target_, goal_x, goal_y);
  skeleton_fit(goal_x, goal_y, weight, bone_goal_, std::max(1, static_cast<int>(std::lround(settings_.skeleton_iterations))));
  since_restart_ += seconds;
  if (settings_.skeleton_restart > 0.0 && since_restart_ >= settings_.skeleton_restart) {
    since_restart_ = 0.0;
    BoneParameters fresh{};
    fresh.fill(0.0);
    skeleton_fit(goal_x, goal_y, weight, fresh,
                 std::max(1, static_cast<int>(std::lround(settings_.skeleton_restart_iterations))));
    if (skeleton_cost(goal_x, goal_y, weight, fresh) < skeleton_cost(goal_x, goal_y, weight, bone_goal_)) {
      bone_goal_ = fresh;
    }
  }
  {
    std::array<double, kSkinPoints> xs{};
    std::array<double, kSkinPoints> ys{};
    skeleton_skin(bone_goal_, xs, ys);
    std::size_t tip = static_cast<std::size_t>(kSections - 1);
    goal_tip_miss_ = std::hypot(xs[tip] - goal_x[tip], ys[tip] - goal_y[tip]);
    double sum = 0.0;
    for (int point = 1; point < kSections; ++point) {
      std::size_t p = static_cast<std::size_t>(point);
      sum += (xs[p] - goal_x[p]) * (xs[p] - goal_x[p]) + (ys[p] - goal_y[p]) * (ys[p] - goal_y[p]);
    }
    goal_miss_ = std::sqrt(sum / static_cast<double>(kSections - 1));
  }
  double h = seconds / static_cast<double>(kSubsteps);
  double base_omega = 2.0 * kPi * settings_.skeleton_base_hz;
  double tip_omega = 2.0 * kPi * settings_.skeleton_tip_hz;
  for (int substep = 0; substep < kSubsteps; ++substep) {
    for (int k = 0; k < kBoneParameters; ++k) {
      std::size_t i = static_cast<std::size_t>(k);
      double omega = k <= kBoneSegments || static_cast<std::size_t>(k) >= kAnchorX ? base_omega : tip_omega;
      double accel = omega * omega * (bone_goal_[i] - bone_parameters_[i]) - 2.0 * omega * bone_speed_[i];
      bone_speed_[i] += h * accel;
      bone_parameters_[i] += h * bone_speed_[i];
    }
  }
  TongueContour before = contour_;
  skeleton_surface();
  skeleton_contact();
  for (int point = 0; point < kSections; ++point) {
    std::size_t p = static_cast<std::size_t>(point);
    vx_[p] = (contour_.x[p] - before.x[p]) / seconds;
    vy_[p] = (contour_.y[p] - before.y[p]) / seconds;
  }
}

}  // namespace tract
