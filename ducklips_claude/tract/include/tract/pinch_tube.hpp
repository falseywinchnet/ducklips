#pragma once

// The pinch tube: the voice's source as a short stretch of the airway whose
// walls are sections of tissue, pinched nearly shut in the middle. Each
// section is the owner's three pieces:
//   mass       - its wall's inertia;
//   elasticity - its wall's pull on the neighbouring sections' walls;
//   support    - its wall's pull back to the rest shape.
// Air flows through from the trachea to the tract. Every section's wall
// feels the air's pressure there: Bernoulli's up to the narrowest section,
// the tract's past it (the jet separates). Where opposite walls meet they
// press as tissue. Nothing is a prescribed oscillator and there is no fold:
// the sound is the air's resonance (the trachea below, the tract above, the
// air's own mass in the pinch) interacting with the walls' elasticity.
//
// Geometry: each section is a slit `width` across with two opposite walls;
// its opening is its rest gap plus both walls' outward displacement. The
// rest gaps run wide, narrower, pinched (sections 3 and 4, set by the
// posture's rest positions), narrower, wide; the posture's chink is a small
// leak in parallel. Tension (from the intended pitch) stretches the walls:
// every section's elasticity and support times (pitch / natural_pitch)^2.
// CGS units.

#include <array>
#include <cstdint>

namespace tract {

struct PinchSettings {
  double on = 0.0;                 // 1: the pinch tube is the source
  double width = 1.3;              // cm across the walls
  double section_length = 0.4;     // cm along the flow, each section
  double surface_density = 0.1;    // g/cm^2 of wall, at mass factor 1
  double support = 30000.0;        // dyn/cm per wall, at support factor 1
  double elasticity = 30000.0;     // dyn/cm between neighbouring walls, at factor 1
  double loss = 0.2;               // fraction of critical damping
  double contact_stiffness = 300000.0;  // dyn/cm per wall per cm the walls overlap
  double contact_loss = 1.0;       // fraction of critical, while in contact
  double natural_pitch = 130.0;    // Hz at tension 1 (the pitch map's anchor)
  double outer_gap = 0.8;          // cm, rest gap of sections 1 and 6
  double inner_gap = 0.25;         // cm, rest gap of sections 2 and 5
  double wander = 0.05;            // rms share of each section's elasticity
  double wander_corner = 4.0;      // Hz
  double wander_tamping = 1.0;     // the wander times exp(-tamping (support - 1))
  // The separated jet's turbulence pressing on the walls past the pinch
  // (rms share of the jet's kinetic pressure, 1/2 rho v^2), each section's
  // own eddies (one pole at turbulence_corner Hz). 0: off. Measured
  // 2026-10-17: the pinch's cycles were nearly identical (jitter 0.2%,
  // shimmer 4%) against the owner's held vowels (9-14%, 12-19%, same
  // picker); a strictly periodic 90 Hz pulse train is the "oscillatory
  // rattle" the owner hears at low pitch.
  double jet_turbulence = 0.0;
  double turbulence_corner = 1000.0;
  // The harp (owner, 2026-10-18: "feedback and resonance loops. pretend this
  // is a musical instrument like a harp. 1/n2 type stuff"): each wall a
  // string across the slit, fixed at both ends, as string_modes standing
  // waves sin(n pi x / width). The string's tension along the slit gives mode
  // n the stiffness T (n pi / width)^2 (width / 2): harmonic frequencies, and
  // a pluck (the walls meeting) answered about 1/n^2 per mode. Every mode is
  // a resonance looped through the air: it shapes the opening (a lens), the
  // opening the flow, the flow the pressure that drives it. 0 modes: the
  // walls are pistons as before.
  double string_modes = 0.0;        // count, at most kMaxStringModes
  double string_tension = 15000.0;  // dyn along the slit, times the tension (pitch / natural_pitch)^2
  double string_loss_growth = 0.0;  // mode n's damping is loss n^growth
  double string_taper = 0.0;        // cm the rest gap widens toward one end (the posterior chink)
  // The owner's six sections (mass, elasticity, support factors), inlet first.
  double s1_mass = 1.0;
  double s1_elasticity = 2.0;
  double s1_support = 1.0;
  double s2_mass = 2.0;
  double s2_elasticity = 4.0;
  double s2_support = 0.1;
  double s3_mass = 2.0;
  double s3_elasticity = 1.0;
  double s3_support = 0.3;
  double s4_mass = 4.0;
  double s4_elasticity = 0.1;
  double s4_support = 3.0;
  double s5_mass = 0.1;
  double s5_elasticity = 3.0;
  double s5_support = 0.1;
  double s6_mass = 1.0;
  double s6_elasticity = 2.0;
  double s6_support = 1.0;
};

class PinchTube {
 public:
  static constexpr int kSections = 6;
  static constexpr int kMaxStringModes = 12;
  static constexpr int kStringPoints = 16;  // points across the slit for the opening and contact

  void configure(const PinchSettings& settings, std::uint32_t seed);
  // The posture: the pinch sections' rest gaps (cm, wall to wall) and the
  // leak beside them (cm^2).
  void set_rest(double lower_gap, double upper_gap, double leak);
  void set_pitch(double hz) { pitch_ = hz; }
  // One step at the tube rate, from the waves arriving below and above.
  void step(double wave_below, double area_below, double wave_above, double area_above, double& out_below,
            double& out_above);

  double flow() const { return flow_; }             // cm^3/s
  double area() const { return narrowest_ + leak_; }  // cm^2
  double contact() const { return contact_; }       // cm, the deepest overlap of opposite walls
  double tension() const { return tension_; }
  double pressure_below() const { return pressure_below_; }
  double pressure_above() const { return pressure_above_; }
  double opening(int section) const { return opening_[static_cast<std::size_t>(section)]; }

 private:
  void section_factors(int section, double& mass, double& elasticity, double& support) const;
  void areas();
  void pressures(double p_in, double p_out, double area_below);
  void advance_walls(double h);
  void advance_strings(double h);
  void string_openings();
  double gauss();

  PinchSettings settings_;
  std::array<double, kSections> rest_gap_{};
  std::array<double, kSections> y_{};      // cm, each wall's outward displacement
  std::array<double, kSections> v_{};
  std::array<double, kSections> drift_{};  // elasticity wander, share
  std::array<double, kSections> eddy_{};   // jet turbulence, unit rms
  // The strings: each section's mode displacements (cm, at the antinode) and
  // speeds, the displacement, speed and opening at each point across the
  // slit, and the mode shapes there.
  int modes_ = 0;
  std::array<std::array<double, kMaxStringModes>, kSections> mode_y_{};
  std::array<std::array<double, kMaxStringModes>, kSections> mode_v_{};
  std::array<std::array<double, kStringPoints>, kSections> point_y_{};
  std::array<std::array<double, kStringPoints>, kSections> point_v_{};
  std::array<std::array<double, kStringPoints>, kSections> point_opening_{};
  std::array<std::array<double, kStringPoints>, kMaxStringModes> shape_{};
  std::array<double, kMaxStringModes> pressure_weight_{};  // cm^2: uniform pressure to mode force
  std::array<double, kMaxStringModes> mode_loss_{};        // fraction of critical
  std::array<double, kSections> opening_{};  // cm, wall to wall
  std::array<double, kSections> area_{};     // cm^2
  std::array<double, kSections> pressure_{}; // dyn/cm^2 on each section's walls
  double leak_ = 0.0;
  double narrowest_ = 0.0;
  int narrowest_index_ = 0;
  double pitch_ = 120.0;
  double tension_ = 1.0;
  double flow_ = 0.0;
  double inertial_pressure_ = 0.0;
  double contact_ = 0.0;
  double pressure_below_ = 0.0;
  double pressure_above_ = 0.0;
  std::uint64_t random_ = 1U;
};

}  // namespace tract
