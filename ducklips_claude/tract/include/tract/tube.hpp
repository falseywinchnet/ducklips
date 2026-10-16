#pragma once

// A tube of three-piece segment modules with air flowing through its lumen.
//
// Each module is one length of the vocal system made of three distinct
// pieces (owner's brief): an ELASTIC piece (the cover sheet lining the lumen:
// light, held by axial and hoop tension, edge-joined to the neighbouring
// modules' sheets), a MASS piece (the tissue body behind it: heavy, coupled
// to the sheet by distributed springs) and a SUPPORT piece (the frame that
// anchors the body to the skeleton: stiff, positioned by articulation and
// held by effort). The sheet's axial coupling makes the whole tube's lining
// one elastic system; its implicit-midpoint step is the exact sum of all of
// its eigenmodes and costs one tridiagonal solve per step (no global solve).
//
// Air: one junction per node; waves travel between neighbouring nodes in one
// step (nodes are kNodeLength apart). Pressure pushes the sheet; the sheet's
// motion adds or removes air at its junction.

#include "tract/units.hpp"

#include <cmath>

#include <string>
#include <vector>

namespace tract {

// The three pieces of one module type. Per unit wall area unless noted.
struct SegmentMaterial {
  std::string name;
  // Elastic piece (cover sheet).
  // Defaults: soft wall. The body is mass-like across the speech band (about
  // 1.5 g/cm^2 on a soft, heavily damped support, near 40 Hz), as measured
  // tract walls are; the lining rides on it with its own resonances in the
  // kHz range above the formants. (First values: a 0.2 g/cm^2 lining on a 180 Hz spring floated
  // free of the body above 180 Hz and the walls swallowed F1: /i/ F1 650 Hz
  // where the rigid tube gives 298.)
  double cover_density = 1.04;      // g/cm^3
  double cover_thickness = 0.1;     // cm
  double axial_tension = 3000.0;    // dyn/cm, along the tube (edge coupling)
  double hoop_tension = 3000.0;     // dyn/cm, around the lumen
  double cover_loss = 0.15;         // fraction of critical damping of its local mode
  // Coupling between the elastic and mass pieces.
  double coupling_stiffness = 4.0e8;  // dyn/cm^3 (lining on body: about 10 kHz; at 3.5 kHz it absorbed F3-F5)
  // Mass piece (tissue body).
  double body_density = 1.04;       // g/cm^3
  double body_depth = 1.5;          // cm
  double body_loss = 0.2;           // (unused: lining-body friction is cover_loss)
  // Support piece (frame to skeleton).
  double support_stiffness = 1.0e5;  // dyn/cm^3 (body on support: about 40 Hz)
  double support_loss = 1.0;
};

// Viscous loss per cm of duct of the given area (VocalTractLab 2.4).
double viscous_resistance_per_cm(double area);

// A series element between two travelling waves on a line of impedance z:
// viscous resistance, the kinetic loss of a jet leaving a narrowing
// (Borda-Carnot, only where the air expands: narrow -> wide) and a series
// pressure source. `forward` heads toward the second side, `backward` toward
// the first; returns the flow (first -> second +) and sets the waves that
// leave toward each side.
double series_element(double forward, double backward, double z, double resistance, double source,
                      double first_area, double second_area, double& to_second, double& to_first);

// Built-in module types (prototype values, not measured anatomy).
SegmentMaterial soft_wall();        // pharynx, cheek, mouth floor
SegmentMaterial firm_wall();        // tongue surface braced, velum
SegmentMaterial hard_wall();        // palate, teeth region (bone support)
SegmentMaterial ring_wall();        // trachea (cartilage rings)
SegmentMaterial lip_wall();         // lips

struct TubeNode {
  double rest_area = 1.0;     // cm^2, set by the support frame (articulation)
  int material = 0;           // index into the tube's materials
};

class Tube {
 public:
  Tube(std::vector<SegmentMaterial> materials, std::vector<TubeNode> nodes);

  Index count() const { return static_cast<Index>(nodes_.size()); }
  // The moving open end (lips): only the first active_count() nodes carry
  // air; the radiation load sits on the last of them, behind an air mass
  // for the remaining fraction of a node, so the acoustic length follows
  // `length_cm` continuously. Nodes past the end keep their lining (and are
  // brought in quietly when the tube lengthens).
  void set_active_length(double length_cm);
  Index active_count() const { return active_count_; }
  double area(Index node) const { return area_[static_cast<std::size_t>(node)]; }
  double pressure(Index node) const { return pressure_[static_cast<std::size_t>(node)]; }
  double wall_flow(Index node) const { return wall_flow_[static_cast<std::size_t>(node)]; }
  double total_wall_flow() const { return total_wall_flow_; }
  double frame_rate(Index node) const { return radius_rate_[static_cast<std::size_t>(node)]; }
  // Air flow along link node -> node + 1 this step (cm^3/s, toward the lips +).
  double link_flow(Index link) const { return link_flow_[static_cast<std::size_t>(link)]; }

  // Articulation: the support frame sets each node's rest area. Applied at
  // control rate; the lining follows through its own mechanics.
  void set_rest_area(Index node, double area);
  // A one-pole low-pass on every hop (k; 0 off): a loss growing with frequency.
  void set_wave_smoothing(double k) { wave_smoothing_ = k; }
  double rest_area(Index node) const { return nodes_[static_cast<std::size_t>(node)].rest_area; }
  void set_support_offset(Index node, double offset_cm);
  // Holding effort varies the support piece's stiffness (1 = nominal).
  // The support's holding factor glides to this (2 ms), recompiling the
  // node's pieces every step on the way.
  void set_support_factor(Index node, double factor) { support_target_[static_cast<std::size_t>(node)] = factor; }
  // Take new rest areas and support factors (control rate). The frame does
  // not jump there: each node's rest radius follows its target as a
  // critically damped mass (about 1 ms), and the pieces' masses and
  // stiffnesses are recompiled every step while it moves. (Jumps at the
  // 1 kHz control rate rang the tube as a 1 kHz comb and clicked at every
  // movement.) A moving frame also pushes or draws air at its node.
  // `ramp_steps`: the steps until the next command; the target is ramped
  // linearly to the new radius over them, so it has no steps either.
  void update_rest_geometry(int ramp_steps);
  // Place the frame at the current rest areas at once, at rest (for the
  // body's starting posture, before any sound).
  void settle_geometry();
  // Nodes first..last are bounded by an incompressible body (the tongue):
  // frame motion there moves air from node to node but its net volume over
  // the span is zero. (last < first: no such span.)
  // Nodes from `node` on are outside air (a flaring open end): their area
  // changes push no air. (Default: none.)
  void set_quiet_from(Index node) { quiet_from_ = node; }
  void set_conserved_span(Index first, Index last) {
    conserved_first_ = first;
    conserved_last_ = last;
  }

  // End ports. A port couples an end node to something outside the tube
  // through one-step waves: `incoming` is the wave arriving from outside
  // this step; after step() the wave leaving toward outside is
  // outgoing_left() / outgoing_right(). port_area sets the port impedance.
  void set_left_port(double incoming, double port_area);
  void set_right_port(double incoming, double port_area);
  double outgoing_left() const { return outgoing_left_; }
  // Flow entering the first node through the left port this step, cm^3/s.
  double left_flow() const { return (left_incoming_ - outgoing_left_) * left_port_area_ / kRhoC; }
  double outgoing_right() const { return outgoing_right_; }
  // Side branches (nose at the velum, sinuses, piriform recesses, the
  // under-tongue pocket): a branch joins a node through a port of the given
  // area. Each step the branch's arriving wave is set; after step() the wave
  // leaving into the branch is side_outgoing(). Returns the port's index.
  int add_side_port(Index node);
  void set_side_port(int port, double incoming, double port_area);
  // Move a side port to another node (the articulatory body's velum and
  // tongue tip move along the tract).
  void move_side_port(int port, Index node) { side_node_[static_cast<std::size_t>(port)] = node; }
  double side_outgoing(int port) const { return side_outgoing_[static_cast<std::size_t>(port)]; }
  // A radiating right end (lips): parallel R-L radiation load.
  void set_radiating_right(bool radiating) { radiating_right_ = radiating; }
  // A radiating opening inside the tube (the lips of a body whose length
  // moves): the radiation load of an opening of `area` (R in parallel with L)
  // at `position` (nodes, continuous), shared between the two nodes around
  // it by its fraction, so it slides continuously along the tube. Commanded
  // once per control tick; the place and area glide there every step
  // (0.5 ms), and the node pair follows the glided place. (Stepped once per
  // tick, the load's weights and size put lines at 8 and 16 kHz under the
  // whole voice, the owner's "ringing"; gliding the share alone and clamping
  // it at each node change jumped the load: spikes at voicing onsets.)
  // position < 0: off.
  void set_opening(double position, double area);
  double radiated_flow() const { return radiated_flow_; }

  // A series pressure source (dipole turbulence) on the air path between
  // node and node + 1, for this step only.
  void add_series_source(Index node, double pressure);

  void step();

 private:
  void compile_materials();
  void compile_node(std::size_t i);
  void compile_link(std::size_t i);
  void refresh_geometry();
  void step_lining();
  void step_air();
  void follow_opening();
  void scatter_links();

  std::vector<SegmentMaterial> materials_;
  std::vector<TubeNode> nodes_;
  // Per node, refreshed with geometry.
  std::vector<double> rest_radius_;    // current (gliding) rest radius
  std::vector<double> radius_target_;  // where the frame is going (ramped over the control interval)
  std::vector<double> radius_goal_;    // the latest commanded radius
  std::vector<double> radius_ramp_;    // target change per step
  int ramp_left_ = 0;
  Index active_count_ = 0;
  double end_extra_ = 0.0;             // cm of tube past the last active node's centre
  std::vector<double> radius_rate_;    // cm/s, the frame's radial speed
  std::vector<double> frame_flow_;     // cm^3/s the frame draws from each node's air
  std::vector<double> forward_smooth_;   // each link's waves after the hop's low-pass
  std::vector<double> backward_smooth_;
  double wave_smoothing_ = 0.0;           // the hop low-pass's k (0: lossless hops)
  std::vector<double> frame_low_;      // the frame's flow through the first of two 250 Hz poles
  std::vector<double> frame_lower_;    // and through the second
  double frame_gain_ = 1.0 - std::exp(-2.0 * 3.14159265358979323846 * 250.0 * kStepSeconds);
  Index quiet_from_ = 1 << 30;
  Index conserved_first_ = 0;
  Index conserved_last_ = -1;
  std::vector<double> surface_;        // lining area of the node, cm^2
  std::vector<double> cover_mass_;     // g
  std::vector<double> body_mass_;      // g
  std::vector<double> hoop_k_;         // dyn/cm
  std::vector<double> coupling_k_;     // dyn/cm
  std::vector<double> support_k_;      // dyn/cm
  std::vector<double> axial_k_;        // dyn/cm, between node i and i + 1
  std::vector<double> cover_c_;        // dyn s/cm
  std::vector<double> body_c_;
  std::vector<double> support_offset_; // cm, articulated frame displacement
  std::vector<double> support_factor_; // holding effort on the support stiffness
  std::vector<double> support_target_; // where it is gliding
  // Lining state.
  std::vector<double> cover_q_;        // radial displacement, cm (outward +)
  std::vector<double> cover_v_;
  std::vector<double> body_q_;
  std::vector<double> body_v_;
  // Tridiagonal workspace (no allocation in the step).
  std::vector<double> diag_;
  std::vector<double> upper_;
  std::vector<double> rhs_;
  std::vector<double> scratch_;
  // Air.
  std::vector<double> area_;           // current node area, cm^2
  std::vector<double> link_area_;      // between node i and i + 1
  std::vector<double> forward_;        // wave leaving node i to the right (arrives at i + 1 next step)
  std::vector<double> backward_;       // wave leaving node i + 1 to the left (arrives at i next step)
  std::vector<double> forward_next_;
  std::vector<double> backward_next_;
  std::vector<double> series_;         // series pressure source on link i this step
  std::vector<double> pressure_;
  std::vector<double> wall_flow_;      // cm^3/s into the air from the lining
  std::vector<double> link_resistance_; // viscous series resistance of link i
  std::vector<double> link_flow_;       // flow through link i, cm^3/s
  int resistance_countdown_ = 0;        // steps until the viscous law is refreshed
  std::vector<double> resistance_step_; // viscous resistance glide per step
  std::vector<char> fast_node_;         // node's frame moving fast this step
  double total_wall_flow_ = 0.0;
  // Side ports.
  std::vector<Index> side_node_;
  std::vector<double> side_incoming_;
  std::vector<double> side_area_;
  std::vector<double> side_outgoing_;
  // Ports.
  double left_incoming_ = 0.0;
  double left_port_area_ = 0.0;
  double right_incoming_ = 0.0;
  double right_port_area_ = 0.0;
  double outgoing_left_ = 0.0;
  double outgoing_right_ = 0.0;
  bool radiating_right_ = false;
  double radiation_inertance_flow_ = 0.0;
  double radiation_previous_pressure_ = 0.0;
  double radiated_flow_ = 0.0;
  Index opening_node_ = -1;
  double opening_share_ = 0.0;
  double opening_area_ = 1.0;
  double opening_position_ = -1.0;       // nodes, the load's glided place
  double opening_position_target_ = -1.0;
  double opening_area_target_ = 1.0;
  double opening_flow_[2] = {0.0, 0.0};      // the full load's inertance flow at each of its two nodes
  double opening_pressure_[2] = {0.0, 0.0};  // and the pressure there last step
};

}  // namespace tract
