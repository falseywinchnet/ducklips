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

#include "tract/units.hpp"

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
  void fold_geometry();
  void fold_air_areas();
  void fold_pressures(double flow, double p_below, double p_above, double area_below);
  double jitter_gauss();

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
  // Fold-tube mode.
  static constexpr int kFoldSlices = 6;
  static constexpr int kFoldSubsteps = 8;
  bool fold_tube_ = false;
  FoldTubeSettings fold_;
  double cover_q_[kFoldSlices] = {};      // cm, each slice's cover displacement (opening positive)
  double cover_v_[kFoldSlices] = {};
  double body_q_[kFoldSlices] = {};       // cm, its body
  double body_v_[kFoldSlices] = {};
  double slice_area_[kFoldSlices] = {};   // cm^2, both folds' opening at the slice
  double slice_open_[kFoldSlices] = {};   // cm of fold length open at the slice
  double slice_alpha_[kFoldSlices] = {};  // closed share of the length
  double slice_star_[kFoldSlices] = {};   // the rest gap's offset where the slice is in contact
  double slice_pressure_[kFoldSlices] = {};  // dyn/cm^2 on the slice's medial surface
  double muscle_ = 0.0;                   // the muscle's stiffness fluctuation (share)
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
