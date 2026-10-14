#pragma once

// The voice body. Air moves through connected tubes as pressure waves that
// cross one section per step and scatter at every junction (the
// instantaneous pressure field). Every tube wall, and each vocal fold, is
// built from three distinct pieces:
//
//   cover    - light and elastic; faces the air, joined to its neighbours
//              at its edges, so ripples run along it (edge eigenwaves)
//   mass     - heavy; carries the inertia under the cover
//   support  - cartilage, bone or muscle tone; holds the shape the speaker
//              asks for, and wanders slightly with the effort of holding it
//
// No piece is driven by a prescribed waveform. Sound starts because the
// lungs push air through the adducted folds: the folds open and close on
// the flow, and the flow is bent by the pressure waves that return from the
// trachea below, the mouth and nose above, the sinuses and the space under
// the tongue. Every coupling is a local closed-form update; there is no
// iterative solver.

#include "voice/score.hpp"
#include "voice/units.hpp"

#include <array>
#include <cstdint>

namespace voice {

constexpr Index kOralSections = 48;     // glottis to lips, 17.5 cm
constexpr Index kSubSections = 40;      // trachea and main bronchi, 14.6 cm
constexpr Index kNasalSections = 32;    // velum to nostrils, 11.7 cm
constexpr Index kMaxSections = 48;
constexpr Index kVelumJunction = 25;    // oral junction where the nose joins
constexpr Index kFloorJunction = 40;    // oral junction above the floor space
constexpr Index kPiriformJunction = 6;  // oral junction where the piriform recesses open
constexpr Index kPharynxFirst = 8;      // oral junctions of the expandable pharynx (about 3-9 cm)
constexpr Index kPharynxLast = 24;
constexpr Index kSamplesPerControl = 48;  // 1 kHz motor and shape updates
constexpr double kControlSeconds = kSamplesPerControl / kAudioRate;

// Three-piece wall, per unit area of wall.
struct WallMaterial {
  double cover_mass = 0.15;          // g/cm^2 (mucosa)
  double mass_mass = 1.35;           // g/cm^2 (muscle and fat)
  double cover_stiffness = 4.8e8;    // cover to mass, dyn/cm^3 (local mode 9 kHz)
  double cover_damping = 5000.0;     // dyn s/cm^3
  double support_stiffness = 2.5e5;  // mass to support, dyn/cm^3
  double support_damping = 1200.0;   // dyn s/cm^3
  double edge_stiffness = 2.0e6;     // cover to neighbouring cover
  double mass_edge_stiffness = 2.0e5;  // mass to neighbouring mass
  double wander = 2.0e-4;            // support wander scale, cm
  double coupling = 1.0;             // 0 makes this wall rigid (diagnostic)
};

struct WallState {
  std::array<double, kMaxSections> cover_y{};
  std::array<double, kMaxSections> cover_v{};
  std::array<double, kMaxSections> mass_y{};
  std::array<double, kMaxSections> mass_v{};
  std::array<double, kMaxSections> support_y{};
  std::array<double, kMaxSections> cover_free{};  // cover velocity before this step's pressure
  std::array<double, kMaxSections> previous_pressure{};
  std::array<double, kMaxSections> cover_force{};   // edge forces held over the step
  std::array<double, kMaxSections> mass_force{};
};

// Trapezoidal (bilinear) update of one wall section's linear state
// [cover y, cover v, mass y, mass v]: x' = step x + pressure (p0 + p1)
// + cover_input f_c + mass_input f_m. Passive walls stay passive.
struct WallStepper {
  std::array<std::array<double, 4>, 4> step{};
  std::array<double, 4> pressure{};
  std::array<double, 4> mass_input{};
};

// A tube of `count` sections. right[i] is the pressure wave arriving at the
// right end of section i; left[i] is the wave arriving at its left end.
struct Tube {
  Index count = 0;
  std::array<double, kMaxSections> area{};       // cm^2, this step
  std::array<double, kMaxSections> area_middle{};  // first smoothing stage of the area
  std::array<double, kMaxSections> area_rate_middle{};  // d(area)/dt, first 50 Hz smoothing stage
  std::array<double, kMaxSections> area_rate{};    // d(area)/dt smoothed below the 1 kHz control steps, cm^2/s
  std::array<double, kMaxSections> right{};
  std::array<double, kMaxSections> left{};
  std::array<double, kMaxSections> next_right{};
  std::array<double, kMaxSections> next_left{};
  std::array<double, kMaxSections> pressure{};   // junction pressure, junction i between i and i+1
  std::array<double, kMaxSections> flow{};       // volume velocity at junction i
  std::array<double, kMaxSections> noise{};      // series pressure source at junction i
  WallMaterial wall_material;
  WallState wall;
  WallStepper stepper;
  double loss = 0.0008;  // viscous attenuation scale, per sqrt(cm^2)
  // Per-step caches of quantities that depend only on area: the loss factor
  // of each section and the wall area at each junction (computed once per
  // step instead of once per use).
  std::array<double, kMaxSections> gain{};
  std::array<double, kMaxSections> wall_surface{};
};

// A side cavity reached through a short neck (Helmholtz resonator).
struct Cavity {
  double neck_area = 0.1;    // cm^2
  double neck_length = 0.5;  // cm
  double volume = 10.0;      // cm^3
  double damping = 0.15;     // fraction of critical
  double flow = 0.0;         // into the cavity, cm^3/s
  double pressure = 0.0;     // dyn/cm^2
  double sink_alpha = 0.0;   // this step's flow = alpha + beta * p
  double sink_beta = 0.0;
};

// Radiation from an opening into the room: inertance in parallel with
// resistance.
struct Radiation {
  double inertive_flow = 0.0;
  double radiated_flow = 0.0;
  double previous_flow = 0.0;
};

// One vocal fold (the pair is symmetric): two cover pieces joined at their
// edge, a body mass under them, and the thyroarytenoid/cartilage support.
struct Fold {
  double lower_x = 0.0;
  double lower_v = 0.0;
  double upper_x = 0.0;
  double upper_v = 0.0;
  double body_x = 0.0;
  double body_v = 0.0;
  double flow = 0.0;          // through the glottis, cm^3/s
  double area = 0.0;          // minimum glottal area, cm^2
  double sub_pressure = 0.0;  // dyn/cm^2 below the glottis
  double supra_pressure = 0.0;
  double contact = 0.0;       // contact depth (cm), for the ear
};

// Second-order movement of one part coordinate toward its target, with a
// slow random wander whose size grows with the effort of holding.
struct Motor {
  double position = 0.0;
  double velocity = 0.0;
  double target = 0.0;
  double omega = 60.0;
  double rest = 0.0;
  double wander = 0.0;
  double wander_scale = 0.0;   // per unit effort
  double wander_seconds = 0.2;
  double span = 1.0;           // distance from rest that counts as full effort
};

// A control value set once per millisecond and followed per step through
// two first-order stages, so neither it nor its rate of change jumps (a
// stepped shape would buzz at the control rate).
struct Glide {
  double goal = 0.0;
  double middle = 0.0;
  double value = 0.0;

  void reset(double now) {
    goal = now;
    middle = now;
    value = now;
  }
  void advance(double gain) {
    middle += gain * (goal - middle);
    value += gain * (middle - value);
  }
};

struct Ear {
  double mean_gap = 0.0;
  double previous_gap = 0.0;
  double last_crossing = -1.0;
  double heard_hz = 0.0;
  double correction = 0.0;  // learned log tension offset
  double voiced_seconds = 0.0;
};

// The vocal fold's three pieces: two cover pieces joined at their edge
// (lower and upper margin), the thyroarytenoid body under them, and the
// cartilage support set by adduction and cricothyroid tension.
// Values found by tools/fold_search.py (open quotient ~0.5, jitter ~0.25%
// across /a/ and /o/ at 105 and 145 Hz). The lower "thickness" is the depth
// over which subglottal pressure pushes the lower cover: the medial surface
// (~0.3 cm) plus the projection of the sloping conus elasticus beneath it
// (~0.5-0.8 cm). HYPOTHESIS: the fitted 1.03 cm is that effective depth.
// Physiological body-cover folds (Story & Titze 1995 layout; geometry
// after the Titze & Story 2002 rules for a modal male fold: length 1.6
// cm, thickness 0.3 cm split at the nodal point 0.125 cm, cover depth 0.3
// cm, tissue 1.04 g/cm^3). Cover stiffness puts each cover piece near
// 110 Hz on its own; coupling carries the mucosal wave between them.
// The previous materials came from a random search scored on open
// quotient, jitter and pitch only (lower piece 1.03 cm thick, upper 0.0275
// cm, stiffnesses 177x apart): the area opened fast and closed slowly or
// fired a 3 ms pressed spike with a rebound, open quotient 0.22-0.37, and
// the vibration doubled its period at 125-155 Hz (tools/fold_shape.py).
struct FoldMaterial {
  double lower_mass = 0.062;        // g
  double upper_mass = 0.087;        // g
  double body_mass = 0.1;           // g
  double lower_stiffness = 29600.0; // dyn/cm, cover to body, at tension 1
  double upper_stiffness = 41500.0;
  double shear_stiffness = 20000.0; // lower to upper cover edge
  double body_stiffness = 90000.0;  // body to support
  double lower_damping = 0.05;      // fraction of critical
  double upper_damping = 0.05;
  double body_damping = 0.05;
  double convergence = 0.0;         // lower margin rest gap above upper, cm
  double contact_scale = 1.0;       // contact stiffness / cover stiffness
  double contact_damping = 0.5;     // fraction of critical damping while in contact (1.0: pressed folds stopped at 150 Hz; speech B2 147 Hz, pops 2.3/s)
};

struct BodyOptions {
  FoldMaterial fold;
  double fold_length = 1.4;        // cm
  double lower_thickness = 0.125;  // cm (to the nodal point)
  double upper_thickness = 0.175;  // cm
  double abduction_stiffening = 40.0; // cover and body stiffen by this per cm of rest gap past abduction_onset
  double abduction_onset = 0.02;     // cm: the widest rest gap still set for voicing
  double abduction_damping = 0.0;    // cover damping grows by this factor per cm of rest gap past abduction_onset
  double pressed_flattening = 1.0;  // how far pressing flattens the zipper closure shape (0: always zips)
  double fold_mass_exponent = 1.0; // stretch thins the folds: mass / tension^e, stiffness x tension^(2-e) (0.5: amplitude fell 3x from 100 to 170 Hz)
  double chink = 0.04;             // posterior opening at full compression 0, cm^2 (0.02: soft folds settled still after nasals and /h/ and did not restart; 6 silent vowels in a0001-a0060, 2 at 0.04)
  double oral_loss = 0.0008;       // oral tract attenuation per step, per sqrt(cm^2)
  double strouhal = 0.2;           // frication noise corner = strouhal x velocity / slit height
  double frication_band = 3.0;     // the jet noise rises below corner / this (12 dB/oct)
  int frication_stages = 2;        // low-pass stages above the jet noise corner (2: -12 dB/oct; s/sh tops vs bdl+rms)
  double corner_ceiling = 5500.0;  // Hz, the highest that corner goes (fricatives vs bdl+rms; was 8000)
  double jet_growth = 0.4;         // cm over which a jet becomes turbulent
  double jet_decay = 1.2;          // cm over which its turbulence decays
  double free_front = 1.0;         // a raised tip frees the space ahead of it from the tongue body (0: one parabola)
  double labiodental = 0.6;        // jet weight where a lip channel strikes the upper incisors (0.3: /f/ -30 dB re its vowel, rms -18; 0.6: -25)
  double obstacle_reach = 0.8;     // cm a jet stays coherent enough to strike the teeth (bursts vs bdl+rms)
  double piriform_scale = 1.0;     // scales the piriform recesses' cross-section (diagnostic)
  double wave_diffusion = 0.06;    // boundary-layer loss: spatial smoothing of waves per step (0 off)
  bool jet_incoherent = false;     // each junction along the jet has its own turbulence (REJECTED as default: see NOTES)
  bool sink_closed_only = true;    // the pharynx sink acts only while the mouth is shut
  Index sink_last = kPharynxLast;  // last oral junction whose area change draws or pushes air
  bool flow_carry = true;          // keep a section's flow through its area changes (false: the old held waves)
  double glide_ms = 1.5;           // time constant of each of the two control-smoothing stages (was 0.6)
  double voiced_breath = 1.5;      // turbulent breath of vibrating folds (4.5 was set for the pressed source; on the physiological folds it doubled B2: 219 Hz vs human 105)
  double viscous_depth = 1.0;      // depth of the glottal slit for viscous loss, cm (per-pulse clicks vs bdl; was 0.3)
  double reference_hz = 90.0;      // fold natural pitch at tension 1
  double ear_gain = 0.06;          // per glottal cycle
  double hold_noise = 1.0;         // scales all wander
  double aspiration = 0.03;        // scales glottal turbulence (burst and /h/ spectra vs bdl; was 0.09)
  double frication = 0.6;          // scales tract turbulence (fricative levels between bdl and rms)
  double skin = 0.5;               // radiation from the neck wall (burst and /h/ spectra; was 0.3)
  double microphone_hz = 9000.0;   // listener off-axis at a distance: high-frequency roll-off
  double wall_coupling = 1.0;      // 0 makes every wall rigid (diagnostic)
  // The pocket over the blade behind a lateral (/l/) contact: the lateral
  // channels bypass it, so it hangs off the tract as a side resonator.
  double lateral_volume = 2.5;     // cm^3
  double lateral_neck = 0.5;       // cm
  Index lateral_junction = 44;     // oral junction it opens into (just behind 16.7 cm)
  std::uint32_t seed = 1U;
};

struct BodyTrace {
  double pitch_heard = 0.0;
  double glottal_area = 0.0;
  double glottal_flow = 0.0;
  double sub_pressure = 0.0;
  double tension = 0.0;
  double lip_flow = 0.0;
  double nose_flow = 0.0;
  double lip_area = 0.0;
  double velum_area = 0.0;
  double narrowest_area = 0.0;   // smallest oral section area, cm^2
  double narrowest_place = 0.0;  // its distance from the glottis, cm
};

// One jet-turbulence generator's filter state (see Body::step_noise).
struct JetNoise {
  double low = 0.0;     // below a third of the Strouhal corner (removed, two stages)
  double lower = 0.0;
  double band = 0.0;    // low-passed at the corner
  double second = 0.0;  // its second low-pass stage
};

class Body {
 public:
  explicit Body(const BodyOptions& options);

  void apply(const PartEvent& event);
  void render(float* output, Index sample_count);
  [[nodiscard]] BodyTrace trace() const;
  [[nodiscard]] const Tube& oral() const { return oral_; }

  // Adds a pressure pulse travelling up from the glottis (for tract probes).
  void kick(double pressure) { oral_.right[0] += pressure; }

  // Settle every motor at its current target (used before speech starts).
  void settle();

 private:
  void control_tick();
  void advance_glides();
  void shape_tract();
  void step();
  void step_tube_junctions(Tube& tube, Index first, Index last);
  double step_junction(Tube& tube, Index junction, double side_alpha, double side_beta);
  void step_glottis();
  void step_velum();
  void step_piriform();
  void step_cavity_junction(Cavity& cavity, Index junction);
  void diffuse_waves(Tube& tube);
  void prepare_walls(Tube& tube);
  void finish_walls(Tube& tube);
  void wall_sink(const Tube& tube, Index junction, double& alpha, double& beta) const;
  static void refresh_cache(Tube& tube, Index first, Index last);
  void step_fold();
  void step_noise();
  double shape_jet_noise(JetNoise& state, double floor_gain, double corner_gain);
  double step_lips();
  double step_nostrils();
  void wander_motor(Motor& motor, double effort);
  double random_uniform();
  double random_gauss();
  double motor_value(Part part, Index coordinate) const;

  BodyOptions options_;
  Tube oral_;
  Tube sub_;
  Tube nasal_;
  std::array<Tube, 2> piriform_{};  // blind recesses beside the larynx; index 0 at the opening
  std::array<Cavity, 3> sinuses_{};
  std::array<Index, 3> sinus_junction_{};
  Cavity floor_;
  Cavity pocket_;                  // lateral pocket (Part::lateral)
  std::array<double, kMaxSections> diffuse_scratch_{};       // waves before boundary-layer smoothing
  std::array<double, kMaxSections> diffuse_scratch_left_{};
  std::array<double, kMaxSections> diffuse_root_{};
  Radiation lips_;
  Radiation nostrils_;
  Fold fold_;
  Ear ear_;
  std::array<Motor, kPartCount * 2> motors_{};
  std::array<double, kOralSections> oral_target_{};
  std::array<double, kNasalSections> nasal_target_{};
  double tension_ = 1.0;          // muscle state, per control tick
  double tension_target_ = 1.0;
  Glide tension_glide_;
  Glide lung_glide_;
  Glide velum_glide_;
  Glide gap_glide_;
  Glide compression_glide_;
  Glide floor_glide_;
  Glide pocket_glide_;
  double glide_gain_ = 0.0;        // per-step gain of the control smoothing (from glide_ms)
  double lung_pressure_ = 0.0;
  double velum_area_ = 0.0;
  double fold_rest_gap_ = 0.05;
  double compression_ = 0.5;
  double tremor_phase_ = 0.0;
  double time_ = 0.0;
  double output_previous_ = 0.0;
  double highpass_state_ = 0.0;
  double microphone_state_ = 0.0;  // off-axis distance roll-off
  std::array<double, 4> alias_state_{};  // two biquads at the step rate (direct form II)
  std::array<JetNoise, kOralSections> jet_noise_{};  // turbulence filter state per junction (index 0: shared)
  double lung_return_low_ = 0.0;
  double airway_flow_ = 0.0;      // slow flow out of the lungs, cm^3/s  // slow part of the wave returning to the lungs
  double aspiration_low_ = 0.0;    // glottal turbulence below a third of its corner (removed,
  double aspiration_lower_ = 0.0;  // two stages)
  double aspiration_band_ = 0.0;   // above that, low-passed at the Strouhal corner
  Index control_countdown_ = 0;
  bool mouth_shut_ = false;              // a front section is shut (control rate)
  double sink_weight_ = 0.0;             // the pharynx sink's share, following mouth_shut_
  Index jet_narrowest_ = -1;                       // where the jet forms, or -1
  std::array<double, kOralSections> jet_weights_{};  // turbulence strength along its path
  std::uint64_t random_state_ = 1U;
};

}  // namespace voice
