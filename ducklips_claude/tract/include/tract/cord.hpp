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

class Cord {
 public:
  Cord();

  void set_target(const CordPosture& posture, double move_seconds);
  void set_pitch(double hz) { pitch_ = hz; }
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
  double contact() const { return alpha_[0]; }        // contact fraction of the lower edge
  double tension() const { return tension_; }
  double pressure_below() const { return pressure_sub_; }   // dyn/cm^2, at the glottal inlet
  double pressure_above() const { return pressure_supra_; } // dyn/cm^2, at the glottal outlet

 private:
  void advance_posture();
  void compute_geometry();
  void advance_mechanics();

  // Posture (support piece), glided toward its target.
  CordPosture posture_;
  CordPosture guide_;   // the first of the two glide stages
  CordPosture target_;
  double glide_gain_ = 1.0;
  double pitch_ = 120.0;
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
  // Air.
  double area_ = 0.0;
  double flow_ = 0.0;
  double inertial_pressure_ = 0.0;  // pressure across the glottal air mass, dyn/cm^2
  double pressure_sub_ = 0.0;
  double pressure_supra_ = 0.0;
  double pressure_lower_ = 0.0;
  double pressure_upper_ = 0.0;
};

}  // namespace tract
