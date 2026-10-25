#pragma once

// The cord module: the two vocal folds, a segment whose facing linings meet.
// It is self-oscillating: the flow through it pushes on its pieces and the
// returning pressures sustain the vibration (no prescribed motion).
//
// Pieces, as in every module: ELASTIC (each fold's lower and upper edge
// springs and the coupling between the edges), MASS (the lower and upper edge
// masses) and SUPPORT (the arytenoid frame that sets the rest gap and the
// cricothyroid stretch that sets tension).
//
// Mechanics follow VocalTractLab 2.4's triangular glottis (Birkholz et al.
// 2011): the rest gap is a wedge, zero at the anterior commissure and set by
// the frame at the vocal processes, so contact starts at the front and
// spreads backward smoothly (contact force quadratic at first touch);
// contact damping grows with the contact fraction; both masses are stepped
// implicitly (backward Euler); tension Q divides the masses and multiplies
// the springs while the edge coupling stays constant. Constants from
// TriangularGlottis.cpp and the M01 speaker file.

#include "tract/pinch_tube.hpp"
#include "tract/units.hpp"

#include <algorithm>
#include <cstdint>

namespace tract {

struct CordPosture {
  double lower_rest = 0.00495;  // cm, rest displacement of each fold's lower edge at the vocal processes
  double upper_rest = 0.0004;   // cm, upper edge
  double chink = 0.0002;        // cm^2, posterior (cartilaginous) opening in parallel
};

CordPosture modal_posture();
CordPosture breathy_posture();
CordPosture voiced_fricative_posture();
CordPosture voiceless_plosive_posture();
CordPosture aspirate_posture();
CordPosture voiceless_fricative_posture();

// The fold-tube mode's tissue (both folds alike; totals over the fold,
// shared equally among its slices; at tension Q = 1). Each slice along the
// flow is the body's three pieces: an ELASTIC cover (its own longitudinal
// restoring, and the sheet's vertical tension to the slices above and below,
// which carries the mucosal wave), a MASS body coupled under it, and the
// SUPPORT, the vocalis muscle, holding the body to the frame. Tension Q (the
// cricothyroid's stretch) divides masses and depth and multiplies stiffness
// and length, as in VocalTractLab's glottis.
struct FoldTubeSettings {
  double depth = 0.30;                 // cm, inferior to superior extent of the vibrating fold
  // Masses and stiffnesses three times the first values (same pitch): the
  // lighter fold was swung by the tract's load, its pitch stepping 4.0% at
  // nasal closures and releases (median over ARCTIC a0001-a0060; p90 7.1%)
  // against bdl's 2.4% (4.3%), seams through the harmonics the owner saw;
  // x2: 2.8% (5.8%); x3: 2.4% (5.4%), vowel middles 1.0% (bdl 1.4%).
  double cover_mass = 0.09;             // g
  double body_mass = 0.30;              // g
  double cover_stiffness = 36000.0;     // dyn/cm, the cover's own longitudinal restoring
  double coupling_stiffness = 90000.0;  // dyn/cm, cover on body
  double support_stiffness = 150000.0;  // dyn/cm, body on muscle (the vocalis)
  double vertical_stiffness = 60000.0;  // dyn/cm, between neighbouring cover slices
  double cover_loss = 0.1;              // fraction of critical, cover against body
  double support_loss = 0.4;            // fraction of critical, body against frame
  double contact_stiffness = 900000.0;  // dyn/cm, tissue pressed on tissue
  // The air's opening closes through a soft contact this wide (cm): tissue
  // compresses, so the last microns of the gap close smoothly. (With the
  // opening's corner at zero, every closure's kink spread a broadband spike:
  // 4-10 local clicks a second against 0-1 in the owner's and bdl's voices.)
  double contact_softness = 0.0005;
  // The one unsteadiness in the whole fold: the muscle's stiffness wanders
  // by this rms share, drawn every fold step (1.3 us), through one pole at
  // muscle_corner Hz (0: white at the fold step).
  // (Held /aa/, periods timed to 0.1 us: the physics alone repeats to a
  // RAP of 1.8 us; 0.01 through 100 Hz adds about 1.4 us, mean |dT| 4.3 ->
  // 6.2 us. White at the fold step it averages away within each cycle.)
  double muscle_fluctuation = 0.01;
  double muscle_corner = 100.0;
  // Open-loop tension for the intended pitch: Q = 1 + (pitch - natural) / per.
  // (Held /aa/ at 650 Pa with the elongation below; 217 / 217 gave 104.8 /
  // 116.8 / 136.3 Hz for 100 / 115 / 140.)
  double natural_pitch = 197.0;
  double pitch_per_tension = 171.0;
  // Elongation by the swing (Cord::set_elongation): a swinging fold is
  // stretched along its length, so more breath raises the pitch, as in
  // people (2-6 Hz per cmH2O). Without it the blown tube's pitch fell with
  // breath.
  // (Held /aa/, 400 -> 1000 Pa: 3.7 / 3.6 / 3.2 Hz per cmH2O at 100 / 115 /
  // 140 Hz; without it -0.8, from 116 to 111 Hz at 115.)
  double elongation = 70.0;              // 1/cm^2
  double elongation_reference = 0.0026;  // cm^2, the swing at which the stretch is 1
  double elongation_release = 0.01;   // s
  double elongation_floor = 0.9;
  // The muscle anticipates its own breath: the open-loop tension is divided
  // by the stretch the lungs' pressure will bring (the swing it expects,
  // swing_per_pascal (lungs - swing_onset), cm^2). What the tract does to
  // the pressure across the folds is not anticipated. (Held /aa/ at 115 Hz
  // the swing was 1.3e-3 cm^2 at 400 Pa and 4.3e-3 at 1000.)
  double swing_per_pascal = 4.9e-6;
  double swing_onset = 130.0;  // Pa
  // The air's opening along the fold (2026-10-16 experiment): 0, the
  // displacement uniform along the length; 1, a string's first mode.
  double string_mode = 0.0;
  // Strain stiffening of the cover, the coupling and the body's support
  // (1/cm^2; 0 linear springs): each spring pulls k x (1 + eta x^2). Story &
  // Titze (1995) give the body-cover model's springs eta = 100. Linear, the
  // tube rang as a damped sine with a corner at closure (flow H2 17 dB under
  // H1); the owner heard "the same mechanical tone from the vocal cords".
  double strain_stiffening = 0.0;
  // Parallel strips along the fold's length (1: one uniform fold; up to 3,
  // front first). Each is a column of the slices' three pieces with its share
  // of the length, its own stiffness and mass factors, and longitudinal
  // tension to its neighbours (dyn/cm across the fold, cover and body).
  // Under the same air they answer it differently and pull on each other;
  // the glottis is their sum. (The owner, of the sheet tube's sections: "the
  // parallel pieces allowed adding body because they responded differently,
  // they also responded to each other, and they cumulatively shaped the
  // resonant tone itself".)
  double strips = 1.0;
  double strip1_stiffness = 1.6;  // the front: stiff (the anterior macula)
  double strip1_mass = 1.0;
  double strip2_stiffness = 0.6;  // the middle of the membranous fold: elastic
  double strip2_mass = 0.9;
  double strip3_stiffness = 1.0;  // the back, toward the vocal process: adjusted
  double strip3_mass = 1.1;
  double strip_coupling = 20000.0;
  // Each slice's own pieces (front-to-back sections along the flow, slice 1
  // the inlet, 6 the outlet): factors on its mass piece (the body's mass),
  // its elastic piece (the cover's stiffness and its coupling to the body)
  // and its support (the body's spring to the muscle). All 1: six alike.
  // (The owner's sheet-tube sections: "the interplay gave good resonance".)
  double slice1_mass = 1.0;
  double slice1_elasticity = 1.0;
  double slice1_support = 1.0;
  double slice2_mass = 1.0;
  double slice2_elasticity = 1.0;
  double slice2_support = 1.0;
  double slice3_mass = 1.0;
  double slice3_elasticity = 1.0;
  double slice3_support = 1.0;
  double slice4_mass = 1.0;
  double slice4_elasticity = 1.0;
  double slice4_support = 1.0;
  double slice5_mass = 1.0;
  double slice5_elasticity = 1.0;
  double slice5_support = 1.0;
  double slice6_mass = 1.0;
  double slice6_elasticity = 1.0;
  double slice6_support = 1.0;
  // The elasticity itself wanders in each slice (rms share, through one
  // pole at elasticity_corner Hz), tamped by the slice's support: the
  // wander times exp(-elasticity_tamping (support - 1)). ("letting the
  // elasticity itself wander in a few places ... tamped by stiffness, gave
  // organic character", the owner.)
  double elasticity_wander = 0.0;
  double elasticity_corner = 20.0;
  double elasticity_tamping = 1.0;
  // ... and of the vertical coupling between slices (1/cm^2; 0 linear).
  double vertical_stiffening = 0.0;
  // The contact force's power of the penetration (1 linear; 1.5 Hertzian).
  double contact_exponent = 1.0;
};

class Cord {
 public:
  Cord();

  void set_target(const CordPosture& posture, double move_seconds);
  void set_pitch(double hz) { pitch_ = hz; }
  // The lungs' pressure (Pa) the muscle expects (FoldTubeSettings::swing_per_pascal).
  void set_expected_breath(double pascal) { expected_breath_ = pascal; }
  // Holding: the stretching muscle's effort scales tension (1 = intended)
  // and the arytenoid frame's hold shifts both rest gaps (cm).
  void set_holding(double tension_factor, double rest_offset) {
    tension_factor_ = tension_factor;
    rest_offset_target_ = rest_offset;
  }
  // Damping ratios of the lower and upper edge masses (VocalTractLab: 0.1, 0.6).
  void set_damping(double lower, double upper) {
    damping_ratio_[0] = lower;
    damping_ratio_[1] = upper;
  }
  // Elongation by deflection (1/cm^2): a fold swinging x from rest is
  // stretched along its length as a plucked string is, so its tension Q
  // grows as 1 + gamma <x^2>, with <x^2> the lower edge's squared
  // deflection averaged over about a cycle (5 ms). A harder breath swings
  // the folds wider and raises the pitch; without it this glottis's pitch
  // fell with subglottal pressure (150 Hz at 360 Pa, 122 at 950). The
  // factor is divided by its value at the reference deflection (cm^2), so
  // the intended pitch holds at speaking pressure. 0 is VocalTractLab's
  // glottis.
  void set_elongation(double gamma, double reference, double release, double floor) {
    elongation_ = gamma;
    elongation_reference_ = reference;
    elongation_release_ = release;
    elongation_floor_ = floor;
  }
  // Kinematic mode: the folds are moved, not blown. Each edge swings as a
  // cosine at the intended pitch (times the stretching muscle's holding
  // factor), the upper edge lagging by `phase_lag` rad, with an amplitude
  // of `amplitude_gain` (cm per sqrt(dyn/cm^2)) times the root of the
  // transglottal pressure (low-passed at 25 Hz), fading out between 200
  // and 100 Pa. The opening along the fold (the rest gap growing linearly
  // to the vocal process, plus the swing as a half sine along the length)
  // closes through a soft contact `contact_softness` cm wide; jitter is a
  // random walk of the phase, `jitter` (rad per sqrt(s)). After Titze's
  // (1989) four-parameter glottis, as VocalTractLab 2.4's default source.
  // The lungs stay the only energy; the flow is solved as before.
  void set_kinematic(bool kinematic, double amplitude_gain, double phase_lag, double contact_softness,
                     double jitter, std::uint32_t seed) {
    kinematic_ = kinematic;
    amplitude_gain_ = amplitude_gain;
    phase_lag_ = phase_lag;
    contact_softness_ = contact_softness;
    jitter_ = jitter;
    jitter_random_ = seed != 0U ? seed : 1U;
  }
  // VocalTractLab 2019's additions to the kinematic folds (Birkholz,
  // Drechsel & Stone 2019, Interspeech; the math as published, no code):
  //   skew s in [-1, 1]: each cycle's phase warped piecewise-linearly so the
  //     opening (trough to peak of the swing) takes pi (1 + s) of the 2 pi and
  //     the closing the rest; s > 0 opens slowly and shuts fast (their Fig. 3
  //     did not survive extraction; this knee is our reading of the text);
  //   double pulsing d in [0, 1]: over the phase's 4 pi, from trough to trough
  //     (1.5 pi to 3.5 pi) the swing is (1 - d/2)(sin + 1) - 1, every second
  //     cycle opening d less;
  //   flutter FL (Klatt & Klatt 1990, which they follow): the frequency plus
  //     (FL / 50)(f / 100)(sin 2 pi 12.7 t + sin 2 pi 7.1 t + sin 2 pi 4.7 t);
  //   length scaling: the swing times (rest length / length)^k (theirs: k = 1,
  //     x_max = c x_rel sqrt(P) L0 / L).
  // All 0: the kinematic folds as before.
  void set_kinematic_shape(double skew, double double_pulsing, double flutter, double length_scaling) {
    skew_ = std::clamp(skew, -0.9, 0.9);
    double_pulsing_ = std::clamp(double_pulsing, 0.0, 1.0);
    flutter_ = flutter;
    length_scaling_ = length_scaling;
  }
  // VocalTractLab 2.4's default folds, "Geometric glottis 2025" (math from
  // its GeometricGlottis2025.cpp, GPL: mathematics only, our own code):
  //   collision: each edge follows A cos(2 pi t) until it reaches where the
  //     folds begin to touch, A2 = -rest / pi, and through the closed phase
  //     the larger of two exponentials matched in slope to the cosine at t1
  //     and t4, both settling toward A3 = A2 - 0.2 A: the swing eases into
  //     closure and out of it instead of crossing the midline;
  //   pulse shape p in [-1, 1]: t4 = t3 + 0.1 p (p > 0, a later release) or
  //     t1 = t2 + 0.1 p (p < 0, an earlier one), t2, t3 where the cosine
  //     crosses A2;
  //   the swing from the commanded lung pressure less the pressure above the
  //     glottis through 25 Hz low-pass (two Butterworth sections; theirs a
  //     4th-order Chebyshev), times the relative amplitude;
  //   the area along the fold at 32 points, cut at the midline (no soft
  //     contact), plus the chink.
  // Their modal shape (JD2): rests 0.0403 / 0.03995 cm, chink 0.02 cm^2, lag
  // 1.222 rad, pulse shape 0.5, flutter 25%, 800 Pa; jitter 10% (0.176 rad
  // per sqrt(s)). Off: the kinematic folds as before.
  // `contact_softness` (cm, 0 = theirs, cut at the midline): the opening at
  // each of the 32 points passes through the midline softly, softness
  // log(1 + exp(x / softness)), as the older kinematic folds did. Measured
  // 2026-10-18 (NATO fit): the closures clicked where little else sounds
  // above 5 kHz (nasals, voicing onsets), 1.8/s against the reader's 0.2;
  // the cut leaves a kink in the area each time a point shuts.
  void set_geometric_2025(bool on, double pulse_shape, double relative_amplitude, double contact_softness) {
    geometric_2025_ = on;
    pulse_shape_ = std::clamp(pulse_shape, -1.0, 1.0);
    relative_amplitude_ = std::max(relative_amplitude, 0.0);
    softness_2025_ = std::max(contact_softness, 0.0);
  }
  // The 2025 folds' swing times this (0..1): the speaker fades it as the
  // arytenoids abduct (Speaker, BodyLayout::abduction_fade_*), as
  // VocalTractLab's voiceless shapes set a relative amplitude near 0.
  void set_abduction_scale(double scale) { abduction_scale_ = std::clamp(scale, 0.0, 1.0); }
  // The 2025 folds' rest length, thickness (cm) and the pitch they rest at
  // (Hz): JD2's 1.6, 0.45, 120 (male); W02's 1.2, 0.40, 220 (female).
  void set_geometry_2025(double length, double thickness, double pitch) {
    rest_length_2025_ = length;
    rest_thickness_2025_ = thickness;
    rest_pitch_2025_ = pitch;
  }
  // Cycle-to-cycle variation of the 2025 folds (0: none, as before): each
  // cycle draws offsets of the pulse shape, the upper edge's lag (rad), the
  // swing (share) and both rests (cm), these standard deviations, and the
  // values glide linearly across the cycle to them. Real folds never repeat
  // a pulse exactly (asymmetry, mucus, the jet's turbulence); measured
  // 2026-10-18, our harmonics' phase distortion scattered a half to a third
  // of bdl's and the NATO reader's, and phase distortion is blind to timing
  // jitter (it cancels a whole cycle's shift), so only the pulse's shape
  // (or noise in the harmonics) can scatter it.
  void set_cycle_variation(double shape, double lag, double amplitude, double rest) {
    cycle_sd_[0] = shape;
    cycle_sd_[1] = lag;
    cycle_sd_[2] = amplitude;
    cycle_sd_[3] = rest;
  }
  // Fold-tube mode: the folds are blown, never moved: a short tube of
  // kFoldSlices slices of the body's three pieces (FoldTubeSettings), stepped
  // eight times per step. The air's pressure on each slice is Bernoulli's
  // along the converging channel up to where the jet separates (past the
  // narrowest slice, once the channel widens by 20%), and the outlet's above
  // it; the lungs are the only energy.
  void set_fold_tube(bool on, const FoldTubeSettings& settings, std::uint32_t seed) {
    fold_tube_ = on;
    fold_ = settings;
    if (on) {
      set_elongation(settings.elongation, settings.elongation_reference, settings.elongation_release,
                     settings.elongation_floor);
    }
    if (seed != 0U) {
      jitter_random_ = seed;
    }
  }
  bool fold_tube() const { return fold_tube_; }
  // The pinch tube in the cord's place (pinch_tube.hpp): no fold, the
  // airway's pinched, elastic sections.
  void set_pinch(bool on, const PinchSettings& settings, std::uint32_t seed) {
    pinch_on_ = on;
    if (on) {
      pinch_.configure(settings, seed);
    }
  }
  bool pinch() const { return pinch_on_; }
  const PinchTube& pinch_tube() const { return pinch_; }
  // VocalTractLab 2.4's glottal flow (TdsModel.cpp, GPL: mathematics only):
  // the glottis is its two sections (each edge's thickness long), with
  //   viscous loss per section T max(0.03 A^-2.6, 0.0065 sqrt(2250) / A)
  //     (their empirical law, dyn s/cm^5, in place of the Poiseuille slit);
  //   the air's mass rho T / A per section, on areas of at least
  //     `inertance_floor` (theirs: their minimum area, 1e-4 cm^2), with no
  //     inlet or outlet lengths;
  //   the kinetic terms as before (full into a contraction, 20% recovered
  //     out of a glottal section), which are theirs.
  // Off: the flow as before.
  void set_vtl_flow(bool on, double inertance_floor) {
    vtl_flow_ = on;
    vtl_inertance_floor_ = std::max(inertance_floor, kMinimumArea);
  }
  // The air mass in the glottis is taken on areas of at least this (cm^2).
  void set_inertance_floor(double area) { inertance_floor_ = area; }
  // Turbulent pressure of the glottal jet acting on the folds' upper
  // surfaces (dyn/cm^2), applied from the next step.
  void set_surface_pressure(double pressure) { surface_pressure_ = pressure; }
  // Turbulent fluctuation of the volume flux through the glottis (cm^3/s),
  // added to the flow both sides see from the next step: a monopole source
  // that enters the tract exactly as the voiced flow does.
  void set_flux_noise(double flow) { flux_noise_ = flow; }

  // One step at the glottal junction. Waves arriving from below (top of the
  // trachea) and from above (bottom of the vocal tube) and the areas of
  // those two sections give the flow; the waves sent back are returned.
  void step(double wave_below, double area_below, double wave_above, double area_above, double& out_below,
            double& out_above);

  double flow() const { return flow_; }               // cm^3/s
  double area() const { return area_; }               // narrowest glottal area, cm^2
  double lower_area() const { return section_area_[0]; }
  double upper_area() const { return section_area_[1]; }
  // Contact fraction: of the lower edge, or in the fold tube of its most
  // closed slice.
  double contact() const;
  double tension() const { return tension_; }
  double pressure_below() const { return pressure_sub_; }   // dyn/cm^2, at the glottal inlet
  double pressure_above() const { return pressure_supra_; } // dyn/cm^2, at the glottal outlet

 private:
  void advance_posture();
  void compute_geometry();
  void advance_mechanics();
  void advance_kinematics();
  void advance_fold_tube();
  int strip_count() const;
  double strip_stiffness(int strip) const;
  double strip_mass(int strip) const;
  void slice_factors(int slice, double& mass, double& elasticity, double& support) const;
  void fold_geometry();
  void fold_air_areas();
  void fold_pressures(double flow, double p_below, double p_above, double area_below);
  double jitter_gauss();
  double shaped_swing(double cosine_phase) const;
  double collided_swing(double cosine_phase, double amplitude, double rest) const;
  double supraglottal_low(double pressure);
  void advance_geometric_2025(double frequency);

  // Posture (support piece), glided toward its target.
  CordPosture posture_;
  CordPosture guide_;   // the first of the two glide stages
  CordPosture target_;
  double glide_gain_ = 1.0;
  double pitch_ = 120.0;
  double expected_breath_ = 0.0;  // Pa
  double tension_factor_ = 1.0;
  double rest_offset_ = 0.0;         // glides toward its target (2 ms): holds arrive at 1 kHz
  double rest_offset_target_ = 0.0;
  double rest_offset_guide_ = 0.0;
  double surface_pressure_ = 0.0;
  double flux_noise_ = 0.0;
  double damping_ratio_[2] = {0.1, 0.6};
  double elongation_ = 0.0;  // 1/cm^2
  double elongation_reference_ = 0.0;  // cm^2
  double elongation_release_ = 0.01;   // s
  double elongation_floor_ = 0.0;      // lowest stretch factor
  double deflection_ = 0.0;  // cm^2, the lower edge's averaged squared deflection
  double tension_ = 1.0;
  double tension_guide_ = 1.0;     // the first of the tension's two stages
  // Mechanics.
  double x_[2] = {0.0, 0.0};
  double x_previous_[2] = {0.0, 0.0};
  double thickness_[2] = {0.24, 0.06};
  double length_ = 1.3;
  double open_length_[2] = {0.0, 0.0};
  double mean_width_[2] = {0.0, 0.0};
  double alpha_[2] = {0.0, 0.0};
  double rest_star_[2] = {0.0, 0.0};
  double section_area_[2] = {0.0, 0.0};
  // Kinematic mode.
  bool kinematic_ = false;
  double amplitude_gain_ = 0.0007;
  double phase_lag_ = 1.22;
  double contact_softness_ = 0.002;
  double jitter_ = 0.17;
  double phase_ = 0.0;
  double transglottal_low_ = 0.0;     // transglottal pressure through the first 25 Hz pole
  double transglottal_lower_ = 0.0;   // and through the second
  std::uint32_t jitter_random_ = 1U;
  double skew_ = 0.0;
  double double_pulsing_ = 0.0;
  double flutter_ = 0.0;
  double length_scaling_ = 0.0;
  bool second_cycle_ = false;   // which of the phase's two cycles (its 4 pi)
  double flutter_time_ = 0.0;   // s
  bool geometric_2025_ = false;
  bool vtl_flow_ = false;
  double vtl_inertance_floor_ = 1.0e-4;
  double pulse_shape_ = 0.0;
  double relative_amplitude_ = 1.0;
  double softness_2025_ = 0.0;
  double rest_length_2025_ = 1.6;
  double abduction_scale_ = 1.0;
  double rest_thickness_2025_ = 0.45;
  double rest_pitch_2025_ = 120.0;
  double cycle_sd_[4] = {0.0, 0.0, 0.0, 0.0};     // shape, lag, amplitude, rest
  double cycle_from_[4] = {0.0, 0.0, 0.0, 0.0};
  double cycle_to_[4] = {0.0, 0.0, 0.0, 0.0};
  double cycle_now_[4] = {0.0, 0.0, 0.0, 0.0};
  double supra_filter_[2][2] = {{0.0, 0.0}, {0.0, 0.0}};  // the two low-pass sections' states
  // Fold-tube mode.
  static constexpr int kFoldSlices = 6;
  static constexpr int kFoldSubsteps = 8;
  bool fold_tube_ = false;
  bool pinch_on_ = false;
  PinchTube pinch_;
  FoldTubeSettings fold_;
  static constexpr int kMaxStrips = 3;
  double cover_q_[kMaxStrips][kFoldSlices] = {};  // cm, each strip's slices' cover displacement (opening positive)
  double cover_v_[kMaxStrips][kFoldSlices] = {};
  double body_q_[kMaxStrips][kFoldSlices] = {};   // cm, its body
  double body_v_[kMaxStrips][kFoldSlices] = {};
  double strip_open_[kMaxStrips][kFoldSlices] = {};   // cm of the strip open at the slice
  double strip_alpha_[kMaxStrips][kFoldSlices] = {};  // closed share of the strip
  double strip_star_[kMaxStrips][kFoldSlices] = {};   // the rest gap's offset where it is in contact
  double slice_area_[kFoldSlices] = {};   // cm^2, both folds' opening at the slice
  double slice_open_[kFoldSlices] = {};   // cm of fold length open at the slice
  double slice_alpha_[kFoldSlices] = {};  // closed share of the length
  double slice_pressure_[kFoldSlices] = {};  // dyn/cm^2 on the slice's medial surface
  double muscle_ = 0.0;                   // the muscle's stiffness fluctuation (share)
  double elasticity_drift_[kFoldSlices] = {};  // each slice's elasticity wander (share)
  // Air.
  double inertance_floor_ = 0.01;
  double area_ = 0.0;
  double flow_ = 0.0;
  double inertial_pressure_ = 0.0;  // pressure across the glottal air mass, dyn/cm^2
  double pressure_sub_ = 0.0;
  double pressure_supra_ = 0.0;
  double pressure_lower_ = 0.0;
  double pressure_upper_ = 0.0;
};

}  // namespace tract
