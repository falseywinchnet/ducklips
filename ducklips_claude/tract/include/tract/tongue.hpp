#pragma once

// The physical tongue: the tongue as a body, moved by its muscles.
//
// Its upper surface is kTongueContourPoints sections from the root's base
// (held at the hyoid, by the epiglottis) over the dorsum to the tip, in the
// jaw's frame (the jaw carries the tongue). Each section has the three
// pieces (the owner's spec):
//   mass       - its inertia;
//   elasticity - its coupling to its neighbours along the surface: stretch
//                (rest spacing) and bending (rest curvature);
//   support    - the passive tissue pulling it back toward the rest shape.
// Two constraints hold for the whole body:
//   volume  - a muscular hydrostat: its midsagittal area (the surface closed
//             through the floor of the mouth) is held at its rest value, so
//             squeezing one part bulges another;
//   contact - the hull (pharyngeal wall, velum, palate, ridge, teeth) stops
//             it; a contact may press contact_depth past the hull, so a
//             closure seals.
// Muscles pull the sections from where they attach:
//   genioglossus posterior (from the genial tubercle, inside the chin; on
//   the root and back of the dorsum: forward, the root out of the pharynx),
//   genioglossus anterior (same origin; on the blade: down),
//   styloglossus (from the styloid, behind and above; on the dorsum: up and
//   back), hyoglossus (from the hyoid; on the dorsum: down and back),
//   superior and inferior longitudinal (along the surface; on the tip: curled
//   up and back, or down), the floor's lift (mylohyoid; on the body: up).
// Each muscle's activation follows its drive through activation_time. The
// drive is feedback: the planner's target shape (the geometric tongue it
// aims for) against where the sections are, projected on each muscle's pull,
// proportional, derivative and integral. Energy enters only through the
// muscles answering that error; nothing is moved along a prescribed path.

#include "tract/articulation.hpp"

#include <array>

namespace tract {

struct TongueSettings {
  double section_mass = 1.0;      // g per section
  double support = 360.0;         // dyn/cm, each section to its rest place (about 3 Hz)
  double support_loss = 0.7;      // fraction of critical
  double elasticity = 20000.0;    // dyn/cm, stretch between neighbours
  double bending = 2000.0;        // dyn/cm, against curvature away from the rest shape
  double elastic_loss = 0.3;      // fraction of critical, between neighbours
  double muscle_force = 1500.0;   // dyn per section at full activation
  double activation_time = 0.03;  // s
  double gain = 1.0;              // drive per (muscle_force per cm of error)
  double damping_gain = 0.08;     // s: the drive's derivative part
  double integral_rate = 6.0;     // 1/s: its integral part
  double contact_depth = 0.04;    // cm a contact may press past the hull
  double tip_pliancy = 0.2;       // support and bending at the tip, times the root's
  double lead = 0.0;              // s: the targets due this far ahead are aimed at now
  double tip_stretch = 1.0;       // stretch stiffness at the tip, times the root's
  double tip_mass = 1.0;          // the tip section's mass, times section_mass (tapering from the root's 1)
  double contact_stiffness = 0.0; // dyn/cm^1.5: soft mucosal contact (0: the hard stop at contact_depth)
  double contact_loss = 600.0;    // dyn s/cm: its loss while pressed
  double statics = 0.0;           // 1: the feedforward by inverse statics at the target (0: influences at rest)
  double smooth_airway = 0.0;     // 1: passage widths by nearest distance (Articulators::nearest_width)
  // 1: the same for the VocalTractLab-style tongue. (Along a fixed cut ray
  // its tip switched the wall below from blade to floor: at jaw 5.699 ->
  // 5.700 degrees one cut went 2.455 -> 1.017 cm^2 and every arc past it
  // moved 0.093 cm; 0.4 cm^3 left the front cavity within 2 ms, a thump at
  // vowel level after speech ended, every ~145 ms as the wander crossed it.)
  double vtl_smooth_airway = 0.0;
  // The tongue's pieces, each a body (1) or pinned to the geometric tongue's
  // outline (0) throughout: root (sections 1-5), dorsum (6-14), blade (15-20),
  // tip (21-23).
  double root_physical = 1.0;
  double dorsum_physical = 1.0;
  double blade_physical = 1.0;
  double tip_physical = 1.0;
  // The pieces' roles by the phone sounding (tongue_roles.hpp), 1 on. Each
  // role is represented by the body, held to its target posture with the
  // role's stiffness (dyn/cm per section; co-contraction, the equilibrium-
  // point view of motor control) and damping (fraction of critical), or
  // pinned to the geometric outline (pinned 1).
  double roles = 0.0;
  double shaper_pinned = 0.0;
  double valve_pinned = 0.0;
  double nozzle_pinned = 0.0;
  double shaper_stiffness = 0.0;
  double valve_stiffness = 16000.0;   // about 20 Hz on a 1 g section
  double nozzle_stiffness = 16000.0;
  double shaper_damping = 1.0;
  double valve_damping = 1.0;
  double nozzle_damping = 1.5;
  // What the roles' stiffness holds toward: the feedback's target (0; it
  // steps as planner events enter the look-ahead), or the motors' smooth
  // posture (1: the geometric outline as the motors hold it, at the
  // tongue's own area), the equilibrium point moving continuously.
  double role_anchor = 0.0;
  // The tongue keeps its volume, not its midsagittal area: flattened it
  // spreads toward the teeth (the midline area falls), bunched it narrows
  // (the area rises). The midline area held follows the planner's shape
  // within rest x (1 - lateral_spread) .. rest x (1 + lateral_narrow), over
  // lateral_time. (0, 0: the midline area is held at rest.)
  double lateral_spread = 0.0;
  double lateral_narrow = 0.0;
  double lateral_time = 0.03;     // s
  double lateral_from_posture = 0.0;
  // What the muscle feedback aims at: 0 the look-ahead target (the next
  // consonant's shape during a vowel: the dorsum rose early and swelled the
  // vowel), 1 the motors' posture now (with role_anchor), 2 by role: the
  // pieces shaping a vowel aim at the posture, the valves and nozzles ahead.
  double feedback_aim = 0.0;
  // The airway's guard (the inverse tongue): the planned shape (the
  // geometric outline as the motors hold it) says where air should be. At
  // each section a share guard_keep of the planned gap between the planned
  // tongue and the hull is protected; a section pressing into it is pushed
  // back by guard_stiffness (dyn/cm) with guard_loss (dyn s/cm) on its
  // intruding speed. Where the plan closes the gap the guard lies on the
  // hull; with the roles, the guard yields on valves and nozzles (it is
  // weighted by the resonator's share). guard 1 builds and measures it;
  // guard_stiffness 0 only measures.
  // The hull as a solid (1): a section inside the closed hull (the skull's
  // side of the pharyngeal wall, velum, palate, ridge and upper teeth) is
  // returned to its nearest hull point, wherever that lies. (By the nearest
  // segment's normal alone, a tip pressed past the ridge's corner read the
  // incisors' face as the nearest and stayed inside the palate.)
  double solid_hull = 0.0;
  // A thick tongue (1): a fat, rounded tip (tip_radius) and a free blade
  // whose underside runs back from the tip, thickening to blade_thickness
  // three sections behind it, then the frenulum straight down to the floor
  // of the mouth; beneath the blade, air. (frenulum_x, frenulum_y: unused.) The underside
  // follows the upper surface as the lower face of a thick beam; what it
  // touches (the solid hull above, the floor and lower incisors below)
  // moves the sections above it. The hydrostat holds the whole body's area.
  double blade_body = 0.0;
  // The spine tongue (1): two overlapping beams inside the tongue carry its
  // shape, and the surface lies around them at set thicknesses (an
  // interpolant around a spine). The base beam is clamped at the root like
  // a fishing rod's butt: stiff there (spine_butt: its pull back toward the
  // rest bend, against the drive), more flexible forward, so it bends back
  // progressively; its forward end (station spine_base_end) may rise at most
  // spine_rise cm above its rest height. The tip beam, rooted on the base
  // beam at station spine_tip_start and overlapping it to spine_base_end
  // (the spine blends from one to the other there), curls and molds the
  // tip. Each beam's bends follow the planner's shape (its spine) as a
  // damped second-order system (spine_base_hz, spine_tip_hz,
  // spine_damping); the solid hull pushes back through the joints. With a
  // set thickness, bending keeps the body's area: no hydrostat projection.
  // The skeleton tongue (1): two internal supports, bent like rods, with
  // the tongue skinned around them. The base support rises from the
  // tongue's base (near the root, low in the body) up and forward through
  // its core to past its centre; the tip support springs from the base
  // support's last quarter and runs out to the tip. Each is
  // kBoneSegments segments that bend (and stretch) together; every surface
  // point rides the supports by its rest attachment (the back of the
  // tongue on the base support, the blade and tip on the tip support,
  // blended between: linear blend skinning). The bends that best skin the
  // planner's shape (its points inside the hull held on the hull) are the
  // goals, followed as damped second-order systems (skeleton_base_hz,
  // skeleton_tip_hz, critical); the base support's forward end may rise at
  // most skeleton_rise above its rest height; the solid hull pushes back
  // through the bends. The tip is fat (tip_radius) with a free blade.
  double skeleton = 0.0;
  double skeleton_base_hz = 10.0;
  double skeleton_tip_hz = 16.0;
  double skeleton_rise = 1.0;
  double skeleton_stiffness = 0.02;  // the fit's pull of each bend toward rest
  double skeleton_tip_weight = 2.0;  // the fit's weight on the blade and tip (sections 18-23)
  double skeleton_iterations = 2.0;  // Gauss-Newton iterations on the goals per control tick
  // Curve skinning (1): the two supports are blended as curves (over the
  // base support's last quarter it fades into the tip support's start),
  // and every skin point rides that one central curve by its rest offset
  // along and across it. (Blending the skinned points instead, a point in
  // the overlap was the average of two places when the supports diverged:
  // /oh/ spiked at section 16.) skeleton_bend_limit (radians): the base
  // support's whole turn beyond it is resisted in the fit.
  double skeleton_curve = 0.0;
  double skeleton_bend_limit = 0.0;  // 0: no limit
  // The anchor rides on the hyoid: it may move up to skeleton_anchor_move
  // (cm) each way, held by skeleton_anchor_stiffness (the fit's pull per cm).
  // skeleton_smoothness: the fit's pull against a bend differing from its
  // neighbour's (each support curves evenly instead of zigzagging).
  double skeleton_anchor_move = 0.0;
  double skeleton_anchor_stiffness = 0.3;
  double skeleton_smoothness = 0.0;
  // Minimal contact (1): the hull's correction changes the bends as little
  // as it can (held to their present values, not to rest, smoothness or the
  // rise limit), leaves the anchor, and turns each bend at most 0.1 rad.
  double skeleton_minimal_contact = 0.0;
  double skeleton_stretch = 0.35;  // each support's length within rest x exp(+-this)
  // Every skeleton_restart seconds the goals are also fitted afresh from
  // rest (skeleton_restart_iterations), and the better fit kept: the warm
  // fit stayed in a local solution set by where the tongue came from (/ow/
  // after a schwa: F1 555, F2 1060 Hz; from rest: 500, 980). 0: off.
  double skeleton_restart = 0.0;
  double skeleton_restart_iterations = 6.0;
  double spine = 0.0;
  double spine_base_end = 16.0;
  double spine_tip_start = 12.0;
  double spine_base_hz = 10.0;
  double spine_tip_hz = 18.0;
  double spine_damping = 1.0;
  double spine_butt = 2.0;
  double spine_rise = 0.5;
  double body_thickness = 1.6;    // cm, the spine to the upper surface over the body
  // The spine's segments lengthen or shorten (the intrinsic muscles: the
  // tongue protrudes and retracts) toward the shape's, within rest x (1 +-
  // spine_stretch), over spine_stretch_time.
  double spine_stretch = 0.25;
  double spine_stretch_time = 0.02;  // s
  // The goals reach (1): the beams' goal bends corrected so the base beam's
  // end and the tip land where the shape's spine has them (damped least
  // squares on the bends, as near as the beams' lengths allow).
  double spine_reach = 1.0;          // 1 both beams, 2 the tip beam only
  double spine_reach_damping = 0.05; // cm^2
  double spine_reach_steps = 3.0;
  double tip_radius = 0.4;        // cm: the tip is twice this thick
  double blade_thickness = 1.5;   // cm, at the back of the free blade
  double frenulum_x = 6.5;
  double frenulum_y = 6.3;
  double guard = 0.0;
  double guard_keep = 0.7;
  double guard_stiffness = 200000.0;
  double guard_loss = 400.0;  // 1: the plan's area from the motors' posture (with role_anchor), not the look-ahead target
  double role_lead = 0.02;        // s: a role is in force this long before its phone
  double role_hold = 0.01;        // s: and after it
  double role_ramp = 0.025;       // s: coming into and out of force
  double pin_stiffness = 1.0e6;   // dyn/cm per section, fully pinned (critically damped)
};

class Tongue {
 public:
  static constexpr int kMuscles = 21;

  // The rest shape (jaw frame); the sections start there.
  void reset(const TongueContour& rest, const TongueSettings& settings);
  bool ready() const { return ready_; }
  // Start the sections at a shape (the first target), at rest velocity.
  void place(const TongueContour& shape);
  // The target shape, given the tongue's own volume: inflated or deflated
  // along its normals to the rest area (the geometric tongue's area ran 16-30
  // cm^2 across a sentence; a hydrostat's stays near its own).
  void set_target(const TongueContour& target);
  // The guard's surface (jaw frame) and each section's intrusion into the
  // guarded airway (cm, 0 outside it).
  const TongueContour& guard_surface() const { return guard_; }
  double intrusion(int section) const { return intrusion_[static_cast<std::size_t>(section)]; }
  // How far the tongue's body (sections 8-16) stands above its rest, on
  // average (cm, jaw frame; up positive).
  double lift() const;
  // How far each section lies inside the hull (cm, 0 outside), after the step.
  double penetration(int section) const { return penetration_[static_cast<std::size_t>(section)]; }
  const TongueContour& geometric() const { return geometric_; }
  // Whether any piece may be pinned to the geometric outline.
  bool uses_geometric() const { return uses_geometric_; }
  // The geometric tongue's outline now (jaw frame), `seconds` after the last.
  void set_geometric(const TongueContour& shape, double seconds);
  // How strongly each piece (root, dorsum, blade, tip) is pinned (0..1).
  void set_pins(const std::array<double, 4>& pieces);
  // Each piece's stiffness toward its target posture and its damping ratio.
  void set_impedance(const std::array<double, 4>& stiffness, const std::array<double, 4>& damping);
  // Each piece's share in the resonator's role (0..1; feedback_aim 2).
  void set_shaping(const std::array<double, 4>& shares);
  // Advance `seconds` with the jaw, larynx and hull of `frame` (the hull as
  // the anatomy last computed it, skull frame).
  void step(double seconds, const Articulators& frame, const MidsagittalTract& anatomy);
  const TongueContour& contour() const { return contour_; }
  const TongueContour& target() const { return target_; }
  double activation(int muscle) const { return activation_[static_cast<std::size_t>(muscle)]; }
  double feedforward(int muscle) const { return feedforward_[static_cast<std::size_t>(muscle)]; }
  double integral(int muscle) const { return integral_[static_cast<std::size_t>(muscle)]; }
  // How far the sections are from the target (cm): rms and largest.
  void target_error(double& rms, double& largest) const;
  // The midsagittal areas (cm^2) of the target shape, the tongue now and at rest.
  double target_area() const { return area_of(target_); }
  // The planner's target before it is brought to the tongue's own area.
  double planned_area() const { return planned_area_; }
  double current_area() const { return area_of(contour_); }
  // The thick tongue's cap, underside and frenulum now (jaw frame; blade_body).
  const TongueUnderside& underside() const { return underside_; }
  double rest_area() const { return rest_area_; }
  double held_area() const { return held_area_; }
  // The feedforward's own estimate of how near the muscles can hold the target (cm rms).
  double reach_error() const { return reach_error_; }

 private:
  static constexpr int kSections = kTongueContourPoints;
  static constexpr int kSubsteps = 8;

  void muscle_pull(int muscle, int section, const TongueContour& shape, const std::array<double, 6>& origins,
                   double& fx, double& fy) const;
  void drive_muscles(double seconds, const std::array<double, 6>& origins);
  void hold_volume();
  void hold_volume_weighted();
  void build_aim();
  void nearest_hull_point(double px, double py, double& qx, double& qy, double& nx, double& ny) const;
  void build_guard();
  void guard_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const;
  void measure_intrusion();
  bool inside_hull(double px, double py) const;
  void hull_corner(int index, double& x, double& y) const;
  void measure_penetration();
  void pin_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const;
  void impedance_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const;
  void prepare_frame(const Articulators& frame, const MidsagittalTract& anatomy);
  void advance(double h, bool touching);
  void measure_influences();
  void solve_feedforward();
  void passive_forces(const TongueContour& shape, std::array<double, kSections>& fx,
                      std::array<double, kSections>& fy) const;
  void solve_statics();
  void shift_to_area(TongueContour& shape, double wanted) const;
  void nearest_hull(int section, double& qx, double& qy, double& nx, double& ny) const;
  void stop_at_hull(double depth);
  void contact_forces(std::array<double, kSections>& fx, std::array<double, kSections>& fy) const;
  static double outline_area(const TongueContour& shape);
  // The body's area: the outline's, or with blade_body the thick tongue's.
  double area_of(const TongueContour& shape) const;
  double thick_area(const TongueContour& shape) const;
  void derive_underside(const TongueContour& shape, TongueUnderside& out) const;
  static int underside_parent(int point);
  double floor_height(double x) const;
  bool inside_floor(double px, double py) const;
  void nearest_floor_point(double px, double py, double& qx, double& qy) const;
  void stop_underside();
  // The spine tongue.
  void spine_reset();
  void spine_place(const TongueContour& shape);
  void spine_bends_of(const TongueContour& shape, std::array<double, kSections>& base,
                      std::array<double, kSections>& tip) const;
  void spine_lengths_of(const TongueContour& shape, std::array<double, kSections>& lengths) const;
  void spine_points_of(const TongueContour& shape, std::array<double, kSections>& sx,
                       std::array<double, kSections>& sy) const;
  void spine_reach_goals(const std::array<double, kSections>& goal_x, const std::array<double, kSections>& goal_y);
  void spine_positions(const std::array<double, kSections>& base, const std::array<double, kSections>& tip,
                       std::array<double, kSections>& sx, std::array<double, kSections>& sy,
                       std::array<double, kSections>& bx, std::array<double, kSections>& by,
                       std::array<double, kSections>& tx, std::array<double, kSections>& ty) const;
  void spine_surface();
  void spine_step(double seconds);
  void spine_contact();
  double spine_up(int station) const;
  double spine_down(int station) const;
  int spine_base_segments() const;
  int spine_tip_root() const;
  // The skeleton tongue (tongue_skeleton.cpp).
 public:
  static constexpr int kBoneSegments = 4;
  // Each support's bends and its stretch, then the anchor's move (x, y).
  static constexpr int kBoneParameters = 2 * (kBoneSegments + 1) + 2;
  static constexpr int kSkinPoints = kSections + kTongueUndersidePoints;
  static constexpr int kCurvePoints = 46;  // the skeleton's central curve
  using CurvePoints = std::array<double, kCurvePoints>;
  // The supports' points now (jaw frame): support 0 the base, 1 the tip.
  double bone_x(int bone, int point) const { return bone_x_[static_cast<std::size_t>(bone)][static_cast<std::size_t>(point)]; }
  double bone_y(int bone, int point) const { return bone_y_[static_cast<std::size_t>(bone)][static_cast<std::size_t>(point)]; }
  // How far the goal pose's skin lies from the fit's goals (cm): at the tip, and rms.
  double goal_tip_miss() const { return goal_tip_miss_; }
  double goal_miss() const { return goal_miss_; }
 private:
  using BoneParameters = std::array<double, kBoneParameters>;
  using BonePoints = std::array<std::array<double, kBoneSegments + 1>, 2>;
  void skeleton_curve_of(const BonePoints& xs, const BonePoints& ys, CurvePoints& cx, CurvePoints& cy) const;
  void skeleton_reset();
  void skeleton_place(const TongueContour& shape);
  void skeleton_pose(const BoneParameters& parameters, BonePoints& xs, BonePoints& ys) const;
  void skeleton_skin(const BoneParameters& parameters, std::array<double, kSkinPoints>& xs,
                     std::array<double, kSkinPoints>& ys) const;
  void skeleton_fit(const std::array<double, kSections>& goal_x, const std::array<double, kSections>& goal_y,
                    const std::array<double, kSections>& weight, BoneParameters& parameters, int iterations,
                    const BoneParameters* around = nullptr) const;
  void skeleton_goals_of(const TongueContour& shape, std::array<double, kSections>& goal_x,
                         std::array<double, kSections>& goal_y) const;
  void skeleton_surface();
  double skeleton_cost(const std::array<double, kSections>& goal_x, const std::array<double, kSections>& goal_y,
                       const std::array<double, kSections>& weight, const BoneParameters& parameters) const;
  void skeleton_contact();
  void skeleton_step(double seconds);
  void rest_on_floor();

  TongueSettings settings_;
  bool ready_ = false;
  TongueContour rest_;
  TongueContour target_;
  TongueContour contour_;
  TongueContour geometric_;
  TongueContour anchor_;
  TongueContour aim_;
  TongueContour guard_;
  std::array<double, kSections> guard_nx_{};
  std::array<double, kSections> guard_ny_{};
  std::array<double, kSections> guard_weight_{};
  std::array<double, kSections> intrusion_{};
  std::array<double, kSections> penetration_{};
  TongueUnderside underside_;
  // The spine tongue's state: the base and tip beams' bends (radians, one
  // per segment: the first relative to the clamp or to the base beam's
  // segment where the tip beam roots), their speeds, the rest bends, the
  // segment lengths, the clamp's heading and the anchor's offset from the
  // root base; the spine's stations and the base beam's at rest.
  std::array<double, kSections> base_bend_{};
  std::array<double, kSections> base_speed_{};
  std::array<double, kSections> base_rest_bend_{};
  std::array<double, kSections> tip_bend_{};
  std::array<double, kSections> tip_speed_{};
  std::array<double, kSections> tip_rest_bend_{};
  std::array<double, kSections> spine_length_{};
  std::array<double, kSections> spine_rest_length_{};
  std::array<double, kSections> spine_goal_length_{};
  std::array<double, kSections> spine_x_{};
  std::array<double, kSections> spine_y_{};
  std::array<double, kSections> base_x_{};
  std::array<double, kSections> base_y_{};
  std::array<double, kSections> tip_x_{};
  std::array<double, kSections> tip_y_{};
  std::array<double, kSections> base_rest_y_{};
  std::array<double, kSections> base_goal_{};
  std::array<double, kSections> tip_goal_{};
  // The skeleton: each support's rest points, rest headings and segment
  // lengths; each skin point's attachment to each support (its segment, the
  // share along it, and its offset along and across it) and its share on
  // the tip support; the bends and stretches (parameters), their goals and
  // speeds; the supports' points now; the anchor's offset from the root base.
  BonePoints bone_rest_x_{};
  BonePoints bone_rest_y_{};
  std::array<std::array<double, kBoneSegments>, 2> bone_rest_heading_{};
  std::array<double, 2> bone_rest_length_{};
  std::array<std::array<int, kSkinPoints>, 2> skin_segment_{};
  std::array<std::array<double, kSkinPoints>, 2> skin_share_{};
  std::array<std::array<double, kSkinPoints>, 2> skin_along_{};
  std::array<std::array<double, kSkinPoints>, 2> skin_across_{};
  std::array<double, kSkinPoints> skin_tip_weight_{};
  int skin_count_ = 0;
  std::array<int, kSkinPoints> curve_segment_{};
  std::array<double, kSkinPoints> curve_share_{};
  std::array<double, kSkinPoints> curve_along_{};
  std::array<double, kSkinPoints> curve_across_{};
  double tip_root_share_ = 0.75;
  BoneParameters bone_parameters_{};
  BoneParameters bone_goal_{};
  BoneParameters bone_speed_{};
  BonePoints bone_x_{};
  BonePoints bone_y_{};
  double goal_tip_miss_ = 0.0;
  double since_restart_ = 0.0;
  double goal_miss_ = 0.0;
  double bone_anchor_dx_ = 0.0;
  double bone_anchor_dy_ = 0.0;
  double clamp_heading_ = 0.0;
  double anchor_dx_ = 0.0;
  double anchor_dy_ = 0.0;
  bool guarding_ = false;
  std::array<double, kSections> shaping_{};
  std::array<double, kSections> geometric_vx_{};
  std::array<double, kSections> geometric_vy_{};
  bool geometric_set_ = false;
  std::array<double, kSections> pin_{};
  bool pinning_ = false;
  std::array<double, kSections> impedance_{};
  std::array<double, kSections> impedance_damping_{};
  bool holding_ = false;
  bool uses_geometric_ = false;
  std::array<double, kSections> vx_{};
  std::array<double, kSections> vy_{};
  std::array<double, kSections> rest_spacing_{};
  std::array<double, kSections> rest_bend_x_{};
  std::array<double, kSections> rest_bend_y_{};
  double rest_area_ = 0.0;
  double planned_area_ = 0.0;
  double held_area_ = 0.0;        // the midline area the hydrostat holds now
  std::array<double, kMuscles> activation_{};
  std::array<double, kMuscles> feedforward_{};
  std::array<std::array<double, kSections>, kMuscles> influence_x_{};
  std::array<std::array<double, kSections>, kMuscles> influence_y_{};
  TongueContour settled_;
  bool measured_ = false;
  double reach_error_ = 0.0;
  // Inverse statics (solve_statics): the hydrostat's pressure and the
  // system's storage, kept here so no step allocates.
  double pressure_ = 0.0;
  std::array<double, kSections> statics_passive_x_{};
  std::array<double, kSections> statics_passive_y_{};
  std::array<std::array<double, kSections>, kMuscles + 1> statics_column_x_{};
  std::array<std::array<double, kSections>, kMuscles + 1> statics_column_y_{};
  std::array<std::array<double, kMuscles + 1>, kMuscles + 1> statics_normal_{};
  std::array<double, kMuscles + 1> statics_right_{};
  std::array<double, 6> origins_{};
  double root_x_ = 0.0;
  double root_y_ = 0.0;
  std::array<double, kMuscles> integral_{};
  // The hull in the jaw's frame, for contact (rebuilt each step).
  std::array<double, 160> hull_x_{};
  std::array<double, 160> hull_y_{};
  int hull_count_ = 0;
};

}  // namespace tract
