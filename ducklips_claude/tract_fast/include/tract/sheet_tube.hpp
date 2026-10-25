#pragma once

// The connected sheet tube, ported from Astra's web/tube (commit 42c7337,
// the owner-accepted below-mouth stage and its connected speech stages).
//
// The wall is one connected, periodic material sheet (length along the tube,
// circumference around it). It is NOT a stack of mass and spring layers:
// sections follow one another along the tube, each tuned for one property —
// a firm SUPPORT section (heavy ring support), a heavy MASS section, a light
// ELASTIC section — and all of them are part of the same sheet network:
// joined at their edges by membrane tension, carried on an elastic
// foundation, loaded by the air on every surface, and connected by delayed,
// lossy tissue paths. The cord is one such (elastic) section; nothing in the
// tube is a prescribed oscillator. The only sustained energy is the lung
// reservoir at the inlet.
//
// With ring supports and a ring-uniform inlet every angular sector moves
// alike, so one representative sector is computed (the axisymmetric family
// of Astra's eight-sector sheet, as her normal playback path). Its elastic
// step is the implicit midpoint of the whole supported sheet:
//   (M + h^2 K / 4) y = M (2 q + h v),  q' = y - q,  v' = (2 y - 4 q)/h - v,
// one tridiagonal sweep (the cached receptance of her engine, exactly).
// SI units throughout.

#include "tract/cord.hpp"
#include "tract/variation.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace tract {

struct SheetMaterial {
  std::string name;
  double density = 1040.0;   // kg/m^3
  double modulus = 8000.0;   // Pa
  double thickness = 0.002;  // m
  double loss = 0.12;        // fraction of critical damping (viscous)
  double fir[3] = {0.74, 0.17, 0.06};  // tissue-path attenuation taps
};

SheetMaterial sheet_material(const std::string& name);  // "soft", "fold", "firm"

// One section of the sheet: its geometry, material, the three property
// factors (mass, elasticity, ring support), and its ring support's rest
// offset and footprint (the articulation handle).
struct SheetBand {
  std::string name;
  std::string role;           // inlet, adapter, firm, soft, cord, outlet, tract
  double length = 0.01;       // m
  double radius = 0.005;      // m (rest radius before the support offset)
  std::string material = "soft";
  double mass = 1.0;
  double elasticity = 1.0;
  double offset = 0.0;        // m, ring support's radial offset
  double anchor = 0.5;        // 0..1, footprint position within the section
  double support = 1.0;
  int points = 8;
};

struct SheetDesign {
  std::vector<SheetBand> bands;
  double pressure = 0.0;          // Pa, reservoir target (source tube only)
  double wander = 0.0;            // share of slow reservoir wander
  double variation = 0.0;         // elasticity variation depth (<= 0.1)
  std::uint32_t seed = 73021U;
  double elasticity_scale = 1.0;  // whole-body elasticity multiplier
  // Acoustic node spacing: when true, each section gets as many nodes as fit
  // at c / rate, so every air path is exactly one sample long. (Astra's
  // engine kept its sections' 4-8 points at 24 kHz; the one-sample minimum
  // delay then made sound cross each 2-6 mm node in 1/24000 s, ~50-140 m/s:
  // an acoustically several-times-longer tube.)
  bool acoustic_spacing = false;
  bool loss_per_length = false;  // path losses per unit length, not per hop
  bool axial_air = false;        // axisymmetric air: one axial path, no ring self-loops
  bool radiating = false;        // the outlet radiates (R parallel L) instead of Astra's resistive end
  // The air through the cord as a quasi-steady orifice, not a linear wave
  // path: the flow loses its kinetic pressure where the jet separates
  // (Bernoulli), plus viscous (Poiseuille) loss and the air's inertia, and
  // the cord's surfaces feel the intraglottal pressure. Without it a narrow
  // section can only flutter acoustically (a balloon neck): the linear
  // waveguide has no flow-dependent loss, so the flow is never cut off.
  bool orifice = false;
  // Outward opening (m) at which the cord's stiffness has doubled (strain
  // stiffening, q^3); 0 leaves it linear.
  double cord_strain = 0.0;
  // With `orifice`: the cord as a slit between two folds of this length (m)
  // rather than a ring. The section's "radius" is then the half-gap; the
  // pressure acts on the folds' fixed medial surfaces (no hoop curvature);
  // the flow area is 2 L (half-gap - contact). A ring's area grows as r^2
  // and its pressed surface with r, so subglottal pressure could blow it
  // open like a balloon; a glottis is a slit. 0 keeps the ring.
  double cord_length = 0.0;
  // A slit cord's taper (m): the gap grows by this from the front of the
  // folds to the back (the arytenoids hold the posterior end open, as in
  // VocalTractLab's triangular glottis), so contact starts at the front and
  // spreads backward and the open area falls quadratically near closure
  // instead of stopping at once. 0: parallel folds.
  double cord_taper = 0.0;
  // The cord's own viscous loss (fraction of critical), or < 0 to keep its
  // material's. The phonation onset pressure grows with the folds' damping
  // (Titze 1988): with the fold material's 0.055 the cord stopped voicing
  // once held 0.3 mm open.
  double cord_loss = -1.0;
  // Scales every non-cord section's viscous loss (wall damping: the tract's
  // main loss at low frequencies; it widens the formants without bleeding
  // the steady flow as a lossy air path does).
  double wall_loss = 1.0;
  // Elastic fold contact (a slit cord): the folds may press into each other
  // past contact, resisted by a spring this many times the section's own
  // stiffness (cubic beyond 50 um) with contact damping, instead of being
  // stopped at contact. With the taper, contact starts at the front and
  // zips backward and the flow tapers off; a pressed rest (half-gap below
  // zero) is medial compression. 0: the clamp (the flow was cut at once:
  // MFDR 2000-3000 cm^3/s/ms with the lungs' branch, a human 200-600).
  double cord_contact = 0.0;
  // Jet losses along the air paths: where the flow leaves a narrowing into a
  // wider section it loses rho |U| U / 2 (1/A1 - 1/A2)^2 (Borda-Carnot), and
  // every path has Poiseuille loss. A series resistance at the expanding end
  // of each axial path, solved inside the junction (its value from the last
  // step's flow). Without it no pressure builds behind a constriction.
  bool jet_loss = false;
  // The air paths' attenuation taps (per hop; Astra's 0.985, 0.010, 0.003).
  // A lossier, spread set is a cavity whose walls take the highs (the nose:
  // narrow passages over a large mucosal surface).
  double air_taps[3] = {0.985, 0.010, 0.003};
  // The taps of the air paths below the cord (trachea and the bands under
  // the folds), when the tract above takes other taps. A first value below
  // 0 uses air_taps everywhere.
  double subglottal_taps[3] = {-1.0, 0.0, 0.0};
  // A one-pole low-pass (corner, Hz) on every air path's arriving waves:
  // a cavity whose walls take the highs over its length. 0: off.
  double air_lowpass = 0.0;
  // How often the sheet's pieces follow the material and the supports (Hz;
  // Astra's 200). A moving support's stiffness then steps at this rate: at
  // 200 Hz every moving wall chopped the sound at 200 Hz (lines at 200 and
  // 400 Hz in a held /f/ whose motor wander moved the lips).
  double control_rate = 200.0;
  // The lungs behind the reservoir inlet: a pressure source through the
  // airways' resistance (Pa s/m^3) and an inertance that opens the branch
  // above `lung_corner` Hz. At DC the inlet then holds the lung pressure
  // under flow (Astra's matched inlet alone lost 160 Pa per 100 cm^3/s, so
  // an open glottis drained the pressure below it); in the audio band the
  // inlet is still the matched one. 0: off.
  double lung_resistance = 0.0;
  double lung_corner = 20.0;
};

// Reads the text design format (data/sheet/*.txt). Returns false on error.
bool read_sheet_design(const std::string& path, SheetDesign& design);

// A port joining this tube's node to another tube (inlet, outlet or side).
struct SheetPort {
  int index = 0;
  double air_y = 0.0;       // admittance of the port's air path
  double air_in = 0.0;      // pressure wave arriving this step
  double air_out = 0.0;     // pressure wave leaving (after step)
  double solid_z = 0.0;     // impedance of the port's tissue path
  double solid_in = 0.0;    // velocity wave arriving this step
  double solid_out = 0.0;   // velocity wave leaving (after step)
};

class SheetTube {
 public:
  // open: Astra's resistive end (admittance A / (8 rho c) per sector, nearly
  // matched: little reflection); radiating: a piston's radiation load, R in
  // parallel with L, which reflects as a real mouth opening does.
  enum class End { reservoir, open, closed, connected, radiating };

  SheetTube(const SheetDesign& design, double rate);

  int count() const { return n_; }
  double rate() const { return rate_; }
  // Ends.
  void set_inlet(End end) { inlet_ = end; }
  void set_outlet(End end) { outlet_ = end; }
  // Reservoir pressure target, Pa.
  void set_pressure(double pascal) { pressure_target_ = pascal; }
  double drive() const { return drive_; }

  // Articulation: a section's ring support moves radially (its rest offset,
  // m). The rest radius glides to it over the next control period; the
  // sheet's pieces are recompiled at control rate.
  void set_band_offset(int band, double offset);
  // Muscle control of a section's elasticity (pitch intent stiffens the
  // cord): multiplies its varying elasticity factor; propagated to sheet,
  // damping, tissue impedance and delay as the variation is.
  void set_band_elasticity(int band, double factor) { band_factor_[static_cast<std::size_t>(band)] = factor; }
  // The tissue's strain memory (Astra's TissueRelaxation, per node): a share
  // of the foundation stiffness relaxes toward the present deformation at a
  // positive rate whose logarithm wanders (sd, correlation `seconds`, 1 kHz).
  // Below that rate the wall creeps and yields; at audio rates it does not.
  // Dissipative: it adds no force of its own.
  void set_strain_relaxation(double fraction, double seconds, double log_rate_sd, std::uint32_t seed);
  // Muscle control of a section's ring support (1 = its design's): the
  // pieces follow at the next control update.
  void set_band_support(int band, double factor) {
    std::size_t b = static_cast<std::size_t>(band);
    if (band_support_[b] != factor) {
      band_support_[b] = factor;
      geometry_dirty_ = true;
    }
  }
  int band_count() const { return static_cast<int>(bands_.size()); }
  int band_first(int band) const { return starts_[static_cast<std::size_t>(band)]; }
  int band_points(int band) const { return counts_[static_cast<std::size_t>(band)]; }
  int band_of(int node) const { return band_at_[static_cast<std::size_t>(node)]; }

  // Ports (connections to other tubes). Returns the port's index.
  int add_port(int node);
  SheetPort& port(int index) { return ports_[static_cast<std::size_t>(index)]; }

  // A volume flow (m^3/s, whole ring) into a node's air for the next step
  // only (a probe or turbulence source); several add.
  void inject_flow(int node, double flow) {
    injected_[static_cast<std::size_t>(node)] += flow;
    any_injected_ = true;
  }

  // Prepare the step: material variation and drive (call before step()).
  // Then step() advances air, tissue and sheet by one sample.
  bool step();

  double pressure(int node) const { return pressure_[static_cast<std::size_t>(node)]; }
  // Flow radiated from a radiating outlet this step, m^3/s (all sectors).
  double radiated_flow() const { return radiated_flow_; }
  // A pressure source (Pa) in series at the lower end of the air path into
  // `node` (from node - 1), for the next step only: a jet's turbulence.
  void add_series_pressure(int node, double pascal) { series_source_[static_cast<std::size_t>(node)] += pascal; }
  // The folds' motion from a triangular glottis (VocalTractLab-style, blown
  // into vibration by the flow it passes) instead of the sheet's own cord
  // section: the sheet's junctions below and above the cord section are its
  // incident waves and impedances; its flow enters them. It steps at its own
  // rate (96 kHz), twice per sheet step at 48 kHz. (The sheet's own slit
  // cord was intermittently chaotic in speech.) Needs `orifice`.
  void attach_glottis(Cord* cord) { glottis_ = cord; }
  // A pressure source in series with the cord orifice for the next step:
  // the glottal jet's turbulence (aspiration).
  void set_glottal_noise(double pascal) { glottal_noise_ = pascal; }
  // Flow from `node` to node + 1 along their axial air path (whole ring,
  // m^3/s; with `jet_loss`), and the narrowest radius between them.
  double axial_flow(int node) const;
  // Volume flow through the cord orifice, m^3/s (with `orifice`).
  double glottal_flow() const { return orifice_velocity_ * orifice_area_; }
  double glottal_velocity() const { return orifice_velocity_; }  // m/s
  double glottal_area() const { return orifice_area_; }          // m^2
  int cord_last() const { return cord_last_; }
  // Take every section's set offset and elasticity at once, without the
  // glide (the opening posture: gliding there moved the whole tract within
  // a few milliseconds and pushed air).
  void settle_geometry();
  double displacement(int node) const { return q_[static_cast<std::size_t>(node)]; }
  double radius(int node) const { return radii_[static_cast<std::size_t>(node)]; }
  double rest_radius(int node) const { return r_[static_cast<std::size_t>(node)]; }
  double node_dx(int node) const { return dx_[static_cast<std::size_t>(node)]; }
  double mass(int node) const { return mass_[static_cast<std::size_t>(node)]; }
  double axial_tie(int link) const { return axial_[static_cast<std::size_t>(link)]; }
  double impedance(int node) const { return impedance_now_[static_cast<std::size_t>(node)]; }
  const double* fir(int node) const { return fir_[static_cast<std::size_t>(node)].data(); }
  double factor(int node) const { return factor_now_[static_cast<std::size_t>(node)]; }
  long contacts() const { return contacts_; }
  bool fault() const { return fault_; }
  // Narrowest current radius over the cord section(s), m (2x is the gap).
  double cord_radius() const;

 private:
  struct Link {
    int a = 0;
    int b = 0;
    bool self = false;
    double weight = 1.0;
    double gain = 0.0;    // air: admittance scale; tissue: impedance
    double scale = 1.0;   // air: this step's sqrt(admittance)
    int delay = 1;
    double h[4] = {0.985, 0.010, 0.003, 0.0};
    double taps[3] = {0.985, 0.010, 0.003};
    double distance = 0.0;
    std::vector<double> ab;  // history of waves a -> b
    std::vector<double> ba;  // history of waves b -> a
    int mask = 15;
    int pos = 0;
    double ia = 0.0;  // wave arriving at a this step
    double ib = 0.0;  // wave arriving at b this step
    // Jet loss (air): series resistances at each end (per sector, Pa s/m^3),
    // this step's port admittances, and the flow a -> b (whole ring, m^3/s).
    double ra = 0.0;
    double rb = 0.0;
    double ya = 0.0;
    double yb = 0.0;
    double flow = 0.0;
    double low_a = 0.0;  // low-passed waves arriving at a and b (air_lowpass)
    double low_b = 0.0;
  };

  void compile();
  void compile_pieces();
  void recompile_pieces();
  void build_links();
  void configure_links();
  void set_filter(Link& link, double delay);
  void scale_taps(double* taps, double distance) const;
  void advance_control();
  void update_variable_material();
  void receive_air();
  void receive_solid();
  void send(std::vector<Link>& links, const std::vector<double>& field, bool air);
  void send_air();
  void update_jet_loss(Link& link);
  void advance_sheet();
  void solve_orifice();
  void solve_glottis();

  std::vector<SheetBand> bands_;
  std::vector<SheetMaterial> materials_;  // per band
  double rate_ = 24000.0;
  double h_ = 1.0 / 24000.0;
  int n_ = 0;
  double length_ = 0.0;
  // Geometry per node.
  std::vector<double> x_, dx_, r_, area0_;
  std::vector<int> band_at_, starts_, counts_;
  // Compiled sheet (nominal).
  std::vector<double> mass_, root_mass_, surface_, tension_, hoop_, foundation_, loss_, impedance_;
  std::vector<double> axial_;                          // n - 1
  std::vector<double> support_diagonal_, support_off_;  // n, n - 1
  std::vector<double> base_diagonal_;                  // foundation + axial ties + support
  std::vector<double> viscous_;                        // nominal viscous damping
  std::vector<std::vector<double>> fir_;
  // Material now (with variation, interpolated per sample).
  std::vector<double> factor_from_, factor_to_, factor_now_;
  std::vector<double> axial_now_, viscous_now_, impedance_now_;
  // Tridiagonal workspace.
  std::vector<double> diag_, off_, rhs_, scratch_;
  // State.
  std::vector<double> q_, v_, pressure_, sum_, admit_, solid_incoming_, solid_admit_, radii_, force_;
  std::vector<Link> air_;
  std::vector<Link> solid_;
  std::vector<SheetPort> ports_;
  // Variation and drive.
  Variation variation_;
  double variation_depth_ = 0.0;
  int control_period_ = 120;
  int control_phase_ = 0;
  long tick_ = 0;
  double pressure_target_ = 0.0;
  double wander_depth_ = 0.0;
  double wander_ = 0.0;
  double drive_ = 0.0;
  std::uint32_t random_ = 0x513fe2U;
  End inlet_ = End::reservoir;
  End outlet_ = End::open;
  bool geometry_dirty_ = false;
  std::vector<double> band_factor_;      // per band, muscle control of elasticity
  std::vector<double> band_support_;     // per band, muscle control of the ring support
  std::vector<double> radius_goal_;      // per node, the rest radius being glided to
  std::vector<double> radius_step_;      // per node, glide per sample
  int gliding_ = 0;                      // samples of glide left
  // Strain relaxation.
  double relax_fraction_ = 0.0;
  double relax_seconds_ = 0.1;
  double relax_sd_ = 0.6;
  int relax_period_ = 48;                // samples per 1 ms
  int relax_phase_ = 0;
  std::uint32_t relax_random_ = 1U;
  std::vector<double> strain_memory_, relax_log_rate_, relax_gain_;
  void advance_relaxation();
  double relax_gauss();
  bool loss_per_length_ = false;
  double radiation_inertance_flow_ = 0.0;  // per sector
  double radiation_previous_pressure_ = 0.0;
  double radiated_flow_ = 0.0;
  // The exit's kinetic loss (jet_loss): the jet leaving the mouth into the
  // room loses rho U^2 / 2A^2, in series with the radiation load. Without it
  // a lip channel built no pressure: /f/ at 334 cm^3/s through 0.088 cm^2
  // had 120 Pa behind it and a 14 m/s jet (the loss alone is about 870 Pa).
  double exit_resistance_ = 0.0;   // per sector, from the last step's flow
  double exit_admittance_ = 0.0;   // this step's port admittance (per sector)
  double exit_history_ = 0.0;      // the radiation load's history flow (per sector)
  double exit_load_ = 0.0;         // the radiation load's impedance (per sector)
  std::vector<double> injected_;  // per node, flow for the next step
  bool any_injected_ = false;
  bool axial_air_ = false;
  // Cord orifice.
  bool orifice_ = false;
  int cord_first_ = -1;
  int cord_last_ = -1;
  double cord_length_ = 0.0;
  double cord_strain_ = 0.0;
  double fold_length_ = 0.0;              // m, slit cord (0: ring)
  double fold_taper_ = 0.0;               // m, gap added from front to back
  double fold_contact_ = 0.0;             // contact stiffness multiple (0: clamp)
  bool soft_contact(std::size_t node) const { return fold_contact_ > 0.0 && slit(node); }
  double rest_of(std::size_t node) const;
  std::vector<double> arc_;               // per node, one sector's width of wall (m)
  bool slit(std::size_t node) const { return fold_length_ > 0.0 && in_cord_[node] != 0; }
  double cord_area(double radius) const;
  double angular_distance(std::size_t node) const;
  double orifice_velocity_ = 0.0;  // m/s through the narrowest ring
  double orifice_area_ = 0.0;      // m^2
  std::vector<double> glottal_pressure_;  // per node (cord nodes only)
  std::vector<char> in_cord_;
  bool jet_loss_ = false;
  std::vector<double> series_source_;  // per node, Pa, next step
  std::vector<int> axial_link_;         // per node, index in air_ of its path to node + 1 (-1: none)
  double glottal_noise_ = 0.0;
  Cord* glottis_ = nullptr;
  double lung_resistance_ = 0.0;
  double air_taps_[3] = {0.985, 0.010, 0.003};
  double subglottal_taps_[3] = {0.985, 0.010, 0.003};
  int first_cord_band_ = -1;  // bands before it lie below the folds
  double air_lowpass_gain_ = 0.0;  // 0: off
  double lung_gain_ = 0.0;     // g = 1 / (2 L / h + R)
  double lung_keep_ = 0.0;     // a = (2 L / h - R) / (2 L / h + R)
  double lung_flow_ = 0.0;     // whole ring, m^3/s
  double lung_drop_ = 0.0;     // drive - inlet pressure, last step
  long contacts_ = 0;
  bool fault_ = false;
};

// A joint between two tubes' nodes (a side branch: the nose at the velum,
// a sinus, the space under the tongue): delayed air and tissue paths as the
// connected system's, with an air opening the speaker may set (the velum's
// port). Prepare before the tubes step, send after.
class SheetJoint {
 public:
  void connect(SheetTube& a, int node_a, SheetTube& b, int node_b, double rate);
  // The air path's opening, m^2; negative: the narrower of the two rings.
  void set_opening(double area) { opening_ = area; }
  void prepare();
  void send();

 private:
  double filter(const std::vector<double>& history, double delay, const double* taps) const;

  SheetTube* a_ = nullptr;
  SheetTube* b_ = nullptr;
  int port_a_ = 0;
  int port_b_ = 0;
  double opening_ = -1.0;
  double air_delay_ = 1.0;
  double solid_delay_ = 1.0;
  double taps_[3] = {0.985, 0.010, 0.003};
  double air_scale_ = 0.0;
  double solid_scale_ = 0.0;
  std::vector<double> air_ab_, air_ba_, solid_ab_, solid_ba_;
  int mask_ = 255;
  int pos_ = 0;
};

// Connected tubes (Astra's version-4 speech system): bidirectional air and
// tissue paths between tube nodes, all read before any tube advances.
class SheetSystem {
 public:
  explicit SheetSystem(double rate) : rate_(rate) {}
  // Loads a .system file (tubes and connections) from a directory.
  bool load(const std::string& path);
  int tube_index(const std::string& id) const;
  SheetTube& tube(int index) { return tubes_[static_cast<std::size_t>(index)]; }
  int tube_count() const { return static_cast<int>(tubes_.size()); }
  int source() const { return source_; }
  int output() const { return output_; }
  // One step; returns the outlet's pressure (Pa) at the output tube's end.
  double step();
  bool fault() const { return fault_; }

 private:
  struct Coupling {
    int tube_a = 0;
    int port_a = 0;
    int tube_b = 0;
    int port_b = 0;
    double air_delay = 1.0;
    double solid_delay = 1.0;
    double taps[3] = {0.985, 0.010, 0.003};
    double air_scale = 0.0;
    double solid_scale = 0.0;
    std::vector<double> air_ab, air_ba, solid_ab, solid_ba;
    int mask = 255;
    int pos = 0;
  };
  void configure(Coupling& coupling);
  void prepare(Coupling& coupling);
  void send(Coupling& coupling);
  double filter(const Coupling& coupling, const std::vector<double>& history, double delay, const double* taps) const;

  double rate_;
  std::vector<std::string> ids_;
  std::vector<SheetTube> tubes_;
  std::vector<Coupling> couplings_;
  int source_ = 0;
  int output_ = 0;
  bool fault_ = false;
};

}  // namespace tract
