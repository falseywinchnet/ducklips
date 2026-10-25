#include "tract/voice_body.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

std::vector<SegmentMaterial> trachea_materials() {
  std::vector<SegmentMaterial> materials;
  materials.push_back(ring_wall());
  return materials;
}

// The subglottal tube, lungs first: the bronchial tree's sections (their
// total area growing geometrically from bronchial_area at the lungs to the
// trachea's lowest), then the trachea.
std::vector<double> subglottal_areas(const BodyLayout& layout) {
  std::vector<double> areas;
  Index bronchial = static_cast<Index>(std::lround(layout.bronchial_length / kNodeLength));
  double top = layout.trachea_areas.empty() ? 2.5 : layout.trachea_areas.front();
  for (Index node = 0; node < bronchial; ++node) {
    double share = (static_cast<double>(node) + 0.5) / static_cast<double>(bronchial);
    areas.push_back(layout.bronchial_area * std::pow(top / layout.bronchial_area, share));
  }
  areas.insert(areas.end(), layout.trachea_areas.begin(), layout.trachea_areas.end());
  return areas;
}

std::vector<TubeNode> trachea_nodes(const std::vector<double>& areas) {
  std::vector<TubeNode> nodes;
  for (double area : areas) {
    TubeNode node;
    node.rest_area = area;
    node.material = 0;
    nodes.push_back(node);
  }
  return nodes;
}

std::vector<TubeNode> vocal_nodes(const BodyLayout& layout) {
  std::vector<TubeNode> nodes;
  for (std::size_t index = 0; index < layout.vocal_areas.size(); ++index) {
    TubeNode node;
    node.rest_area = layout.vocal_areas[index];
    node.material = layout.vocal_materials[index];
    nodes.push_back(node);
  }
  return nodes;
}

// The vocal walls' materials, their support stiffness and loss scaled as
// the layout asks (BodyLayout::wall_stiffness_scale, wall_loss_scale).
std::vector<SegmentMaterial> vocal_materials(const BodyLayout& layout) {
  std::vector<SegmentMaterial> materials = layout.vocal_material_list;
  for (SegmentMaterial& material : materials) {
    material.support_stiffness *= layout.wall_stiffness_scale;
    material.support_loss *= layout.wall_loss_scale;
  }
  return materials;
}

}  // namespace

BodyLayout vowel_layout(const std::vector<double>& section_areas, double section_length) {
  BodyLayout layout;
  // Trachea: 12 cm of cartilage rings, wider just below the larynx
  // (VocalTractLab: 2.5 cm^2, the top sections 3.0 and 4.0).
  Index trachea_count = static_cast<Index>(std::lround(12.0 / kNodeLength));
  for (Index node = 0; node < trachea_count; ++node) {
    double from_top = (trachea_count - node - 0.5) * kNodeLength;
    double area = from_top < 1.0 ? 4.0 - from_top : 2.5;
    layout.trachea_areas.push_back(std::max(area, 2.5));
  }
  // Vocal tube: resample the profile at node centres.
  double total = section_length * static_cast<double>(section_areas.size());
  Index count = static_cast<Index>(std::lround(total / kNodeLength));
  layout.vocal_material_list.push_back(firm_wall());
  layout.vocal_material_list.push_back(soft_wall());
  layout.vocal_material_list.push_back(lip_wall());
  for (Index node = 0; node < count; ++node) {
    double x = (node + 0.5) * kNodeLength;
    double position = x / section_length - 0.5;
    Index below = std::clamp(static_cast<Index>(std::floor(position)), 0,
                             static_cast<Index>(section_areas.size()) - 1);
    Index above = std::min(below + 1, static_cast<Index>(section_areas.size()) - 1);
    double weight = std::clamp(position - below, 0.0, 1.0);
    double area = (1.0 - weight) * section_areas[static_cast<std::size_t>(below)] +
                  weight * section_areas[static_cast<std::size_t>(above)];
    layout.vocal_areas.push_back(area);
    int material = 1;
    if (x < 2.0) {
      material = 0;
    } else if (x > total - 1.0) {
      material = 2;
    }
    layout.vocal_materials.push_back(material);
  }
  return layout;
}

// Dang & Honda (1994) subject 4 nasal areas (cm^2), 19 sections of 0.57 cm
// from the velum to the nostrils (as used by VocalTractLab).
constexpr double kNoseAreas[19] = {1.30, 1.25, 1.30, 1.20, 1.60, 2.50, 3.60, 2.25, 2.50, 1.75,
                                   2.50, 2.10, 1.50, 1.50, 0.75, 1.25, 1.60, 1.65, 1.40};
constexpr double kNoseSection = 0.57;
// Paranasal sinuses: volume cm^3, neck area cm^2, nose section they open from.
constexpr double kSinusVolume[4] = {11.3, 6.8, 33.0, 6.2};
constexpr double kSinusNeck[4] = {0.185, 0.185, 0.145, 0.11};
constexpr int kSinusSection[4] = {8, 9, 11, 12};

SegmentMaterial nasal_wall() {
  SegmentMaterial material = hard_wall();
  material.name = "nasal";
  material.cover_thickness = 0.12;
  material.cover_loss = 0.5;  // mucosa over the turbinates: lossy
  return material;
}

// A closed pocket: a neck of `neck_nodes` nodes of area neck, then a cavity of
// the given volume spread over `cavity_nodes` nodes.
Tube pocket(const SegmentMaterial& material, double neck, int neck_nodes, double volume, int cavity_nodes) {
  std::vector<SegmentMaterial> materials;
  materials.push_back(material);
  std::vector<TubeNode> nodes;
  for (int node = 0; node < neck_nodes; ++node) {
    TubeNode piece;
    piece.rest_area = neck;
    nodes.push_back(piece);
  }
  for (int node = 0; node < cavity_nodes; ++node) {
    TubeNode piece;
    piece.rest_area = volume / (cavity_nodes * kNodeLength);
    nodes.push_back(piece);
  }
  return Tube(materials, nodes);
}

// Resistance of the bronchial tree to steady flow, dyn s/cm^5.
constexpr double kAirwayResistance = 1.5;

// First node above the epilarynx (3.5 cm from the glottis).
constexpr Index kFirstOralNode = 10;

// Variation channels: 0 breath, 1 tension, 2 fold posture, 3 shared wall
// effort, then one per vocal-tube node.
constexpr int kFirstWallChannel = 4;

VoiceBody::VoiceBody(const BodyLayout& layout)
    : trachea_(trachea_materials(), trachea_nodes(subglottal_areas(layout))),
      vocal_(vocal_materials(layout), vocal_nodes(layout)),
      skin_gain_(layout.skin_gain),
      skin_pole_gain_(1.0 - std::exp(-kTwoPi * layout.skin_corner * kStepSeconds)),
      settings_(layout),
      variation_(kFirstWallChannel + static_cast<int>(layout.vocal_areas.size()), layout.seed) {
  vocal_.set_radiating_right(true);
  vocal_.set_smoothing_by_area(layout.hop_loss_by_area > 0.5);
  vocal_.set_viscous_scale(layout.viscous_scale);
  if (layout.ear_latency > 0.0) {
    ear_delay_.assign(static_cast<std::size_t>(std::max(std::lround(layout.ear_latency / kControlSeconds), 1L)), 0.0);
  }
  vocal_.set_radiation_scales(layout.radiation_mass_scale, layout.radiation_resistance_scale);
  vocal_.set_support_log_glide(layout.support_log_glide);
  if (layout.frame_corner != 250.0) {
    vocal_.set_frame_corner(layout.frame_corner);
  }
  vocal_.set_area_limit(layout.vtl_area_limit > 0.5);
  cord_.set_kinematic(layout.kinematic_glottis, layout.kinematic_amplitude, layout.kinematic_phase_lag,
                      layout.kinematic_softness, layout.kinematic_jitter, layout.seed ^ 0x5bd1e995U);
  cord_.set_kinematic_shape(layout.kinematic_skew, layout.kinematic_double_pulsing, layout.kinematic_flutter,
                            layout.kinematic_length_scaling);
  cord_.set_geometric_2025(layout.geometric_2025 > 0.5, layout.pulse_shape, layout.relative_amplitude,
                           layout.contact_softness_2025);
  cord_.set_vtl_flow(layout.vtl_glottal_flow > 0.5, layout.vtl_inertance_floor);
  cord_.set_cycle_variation(layout.cycle_shape, layout.cycle_lag, layout.cycle_amplitude, layout.cycle_rest);
  cord_.set_geometry_2025(layout.fold_length_2025, layout.fold_thickness_2025, layout.fold_rest_pitch_2025);
  cord_.set_fold_tube(layout.fold_tube, layout.fold, layout.seed ^ 0x5bd1e995U);
  cord_.set_pinch(layout.pinch.on > 0.5, layout.pinch, layout.seed ^ 0x2545f491U);
  // Tremor band normalisation: the 2-15 Hz band (two poles at each corner)
  // of unit-variance noise at the control rate, measured once on a
  // separate generator (100,000 control ticks).
  {
    double slow = 0.0;
    double slower = 0.0;
    double fast = 0.0;
    double faster = 0.0;
    double slow_gain = 1.0 - std::exp(-kTwoPi * 2.0 * kControlSeconds);
    double fast_gain = 1.0 - std::exp(-kTwoPi * 15.0 * kControlSeconds);
    std::uint32_t random = 0x9E3779B9U;
    double sum = 0.0;
    constexpr int kWarm = 5000;
    constexpr int kTicks = 100000;
    for (int tick = 0; tick < kWarm + kTicks; ++tick) {
      double uniform_sum = 0.0;
      for (int draw = 0; draw < 4; ++draw) {
        random ^= random << 13U;
        random ^= random >> 17U;
        random ^= random << 5U;
        uniform_sum += static_cast<double>(random) / 4294967296.0;
      }
      double white = (uniform_sum - 2.0) * std::sqrt(3.0);
      slow += slow_gain * (white - slow);
      slower += slow_gain * (slow - slower);
      fast += fast_gain * (white - fast);
      faster += fast_gain * (fast - faster);
      if (tick >= kWarm) {
        sum += (faster - slower) * (faster - slower);
      }
    }
    tremor_scale_ = 1.0 / std::sqrt(sum / kTicks);
  }
  jet_weights_.assign(static_cast<std::size_t>(vocal_.count()), 0.0);
  jet_obstacle_weights_.assign(static_cast<std::size_t>(vocal_.count()), 0.0);
  aspiration_node_ = std::min(static_cast<Index>(std::lround(1.5 / kNodeLength)), vocal_.count() - 2);

  // Nose, opening from the pharynx at the velum (about 8.5 cm above the
  // glottis) and radiating at the nostrils.
  {
    std::vector<SegmentMaterial> materials;
    materials.push_back(nasal_wall());
    std::vector<TubeNode> nodes;
    double length = 19 * kNoseSection;
    Index count = static_cast<Index>(std::lround(length / kNodeLength));
    for (Index node = 0; node < count; ++node) {
      double x = (node + 0.5) * kNodeLength / kNoseSection;
      int section = std::min(static_cast<int>(x), 18);
      TubeNode piece;
      piece.rest_area = kNoseAreas[section];
      nodes.push_back(piece);
    }
    Tube nose(materials, nodes);
    nose.set_radiating_right(true);
    Index velum_node = std::min(static_cast<Index>(std::lround(8.5 / kNodeLength)), vocal_.count() - 2);
    add_branch(0, velum_node, -1.0, nose);
    nose_ = static_cast<int>(branch_tubes_.size()) - 1;
    branch_tubes_[static_cast<std::size_t>(nose_)].set_area_limit(layout.vtl_area_limit > 0.5);
  }
  // Sinuses on the nose: closed resonant cavities with bony walls.
  for (int sinus = 0; layout.sinuses && sinus < 4; ++sinus) {
    Index node = static_cast<Index>(std::lround((kSinusSection[sinus] + 0.5) * kNoseSection / kNodeLength));
    double neck = layout.sinus_neck_scale * kSinusNeck[sinus];
    add_branch(1, node, neck, pocket(nasal_wall(), neck, 1, kSinusVolume[sinus], 4));
  }
  // Piriform recesses beside the larynx, both sides together: 3 cm, 2 cm^3,
  // tapering 1.2 -> 0.13 cm^2, closed (VocalTractLab). They open beside the
  // laryngeal entrance into the hypopharynx, about 2.5 cm above the glottis,
  // through two openings of about 0.3 cm^2. (Joined at 1.5 cm through 1.2 cm^2,
  // inside the narrow epilarynx, they lowered /a/ F1 from 732 to 627 Hz.)
  if (layout.piriform) {
    std::vector<SegmentMaterial> materials;
    // Lined with lossy mucosa (as the nose): with the soft wall's own loss
    // its zero took 3.7 dB from 4-5 kHz over 24 bdl sentences.
    SegmentMaterial mucosa = soft_wall();
    mucosa.name = "piriform";
    mucosa.cover_loss = settings_.piriform_loss;
    materials.push_back(mucosa);
    std::vector<TubeNode> nodes;
    // Depth 2 cm, about 1.8 cm^3 (1.6 -> 0.2 cm^2): its quarter-wave zero
    // falls near 4.5 kHz, where Dang & Honda (1997) measured it. (3 cm
    // deep it cut every vowel off above 3.6 kHz.)
    Index count = static_cast<Index>(std::lround(2.0 / kNodeLength));
    for (Index node = 0; node < count; ++node) {
      TubeNode piece;
      piece.rest_area = settings_.piriform_scale * (1.6 * (1.0 - (node + 0.5) / count) + 0.2);
      nodes.push_back(piece);
    }
    Index entry = std::min(static_cast<Index>(std::lround(3.0 / kNodeLength)), vocal_.count() - 2);
    Tube recess(materials, nodes);
    recess.set_wave_smoothing(settings_.piriform_smoothing);
    add_branch(0, entry, settings_.piriform_port, recess);
  }
  // The cheek pockets (buccal vestibules): between the teeth and the cheeks,
  // on both sides, closed behind the last molars; they open into the mouth
  // along the teeth and move with them (group 3, set_teeth).
  if (layout.buccal_scale > 0.0) {
    constexpr double kBuccalDepth = 3.0;  // cm, canines to molars
    std::vector<SegmentMaterial> materials;
    SegmentMaterial cheek = soft_wall();
    cheek.name = "buccal";
    materials.push_back(cheek);
    std::vector<TubeNode> nodes;
    Index count = static_cast<Index>(std::lround(kBuccalDepth / kNodeLength));
    for (Index node = 0; node < count; ++node) {
      TubeNode piece;
      piece.rest_area = 0.8 * layout.buccal_scale;
      nodes.push_back(piece);
    }
    Index entry = std::clamp(static_cast<Index>(std::lround(teeth_start_ / kNodeLength - 0.5)), static_cast<Index>(1),
                             vocal_.count() - 2);
    add_branch(0, entry, 0.4 * layout.buccal_scale, Tube(materials, nodes));
    branches_.back().group = 3;
    tube_group_.back() = 3;
    buccal_ = true;
  }
  // The pocket under the tongue: the tongue is an object over the mouth
  // floor, and the space beneath its front is a closed side cavity.
  if (layout.under_tongue) {
    Index entry = std::max(vocal_.count() - static_cast<Index>(std::lround(3.0 / kNodeLength)), 1);
    add_branch(0, entry, -2.0, pocket(soft_wall(), 0.4, 1, 1.5, 3));
  }
  if (layout.lateral_pocket > 0.0) {
    Index entry = std::max(vocal_.count() - static_cast<Index>(std::lround(5.0 / kNodeLength)), 1);
    int cavity = std::max(static_cast<int>(std::lround(layout.lateral_pocket_length / kNodeLength)) - 1, 1);
    add_branch(0, entry, -4.0,
               pocket(soft_wall(), layout.lateral_pocket_area, 1,
                      layout.lateral_pocket_area * static_cast<double>(cavity) * kNodeLength, cavity));
  }
}

// Moves the port of every branch of `group` (1 the nose at the velum, 2 the
// pocket under the tongue) to the node at `place` cm from the glottis.
void VoiceBody::move_branch_group(int group, double place) {
  Index node = std::clamp(static_cast<Index>(std::lround(place / kNodeLength - 0.5)), static_cast<Index>(1),
                          vocal_.active_count() - 2);
  for (Branch& branch : branches_) {
    if (branch.parent == 0 && branch.group == group) {
      if (branch.shadow_port >= 0 && node != branch.node) {
        // The node it leaves keeps the area for now and hands it over.
        branch.shadow_node = branch.node;
        branch.shadow_share = 1.0;
        vocal_.move_side_port(branch.shadow_port, branch.shadow_node);
      }
      branch.node = node;
      vocal_.move_side_port(branch.port, node);
    }
  }
}

void VoiceBody::set_velum_place(double place) { move_branch_group(1, place); }

Index VoiceBody::nose_port_node() const {
  for (const Branch& branch : branches_) {
    if (branch.tube == nose_) {
      return branch.node;
    }
  }
  return -1;
}

void VoiceBody::set_teeth(double start, double end) {
  teeth_start_ = start;
  teeth_end_ = end;
  if (buccal_) {
    move_branch_group(3, start);
  }
}

void VoiceBody::set_under_tongue_place(double place) { move_branch_group(2, place); }

void VoiceBody::set_lateral_pocket_place(double place) { move_branch_group(4, place); }

void VoiceBody::add_branch(int parent, Index node, double port_area, Tube tube) {
  Tube& owner = parent == 0 ? vocal_ : branch_tubes_[static_cast<std::size_t>(nose_)];
  Branch branch;
  branch.parent = parent;
  branch.port = owner.add_side_port(node);
  branch.node = node;
  if (parent == 0 && settings_.branch_glide > 0.0) {
    branch.shadow_port = owner.add_side_port(node);
    branch.shadow_node = node;
  }
  branch.port_area = port_area;
  branch.tube = static_cast<int>(branch_tubes_.size());
  branch.group = parent == 1 || port_area == -1.0 ? 1 : port_area == -2.0 ? 2 : port_area == -4.0 ? 4 : 0;
  branch_tubes_.push_back(tube);
  tube_group_.push_back(branch.group);
  branches_.push_back(branch);
}

// Waves cross each side port with a one-step delay, as along any air path,
// and through the port's own series element: a port is a short duct (one
// node long) with viscous loss and the jet loss where it opens out. (As a
// bare lossless line, a nearly shut velum stored the pressure difference
// across it step after step in its two waves; when it opened, the store
// came out as a 1000 cm^3/s gush into the nose against the pressure, and
// a 1 kHz buzz.)
void VoiceBody::exchange_branches() {
  // Ports open and close continuously, not at the control rate: two
  // stages of 1 ms each, so the port's speed has no step at each tick
  // either (one 2 ms stage left a 1 kHz comb at 3-8 kHz through every
  // nasal, 10-20 dB over the spectrum between its teeth).
  double glide = 1.0 - std::exp(-kStepSeconds / 0.001);
  double velum_target = std::max(velum_area_, settings_.velum_leak * oral_open_);
  velum_guide_ += glide * (velum_target - velum_guide_);
  velum_now_ += glide * (velum_guide_ - velum_now_);
  under_tongue_guide_ += glide * (under_tongue_area_ - under_tongue_guide_);
  under_tongue_now_ += glide * (under_tongue_guide_ - under_tongue_now_);
  // A closing glide never reaches zero: about 1.4 s after a port shut, its
  // area underflowed, rho c / area overflowed and the port returned NaN (the
  // first long passage died 20.7 s in). Below a micro-area a port is shut.
  constexpr double kShutPort = 1.0e-6;  // cm^2
  if (velum_now_ < kShutPort) {
    velum_now_ = 0.0;
  }
  if (under_tongue_now_ < kShutPort) {
    under_tongue_now_ = 0.0;
  }
  lateral_pocket_guide_ += glide * (lateral_pocket_area_ - lateral_pocket_guide_);
  lateral_pocket_now_ += glide * (lateral_pocket_guide_ - lateral_pocket_now_);
  if (lateral_pocket_now_ < kShutPort) {
    lateral_pocket_now_ = 0.0;
  }
  for (Branch& branch : branches_) {
    Tube& child = branch_tubes_[static_cast<std::size_t>(branch.tube)];
    Tube& owner = branch.parent == 0 ? vocal_ : branch_tubes_[static_cast<std::size_t>(nose_)];
    double area = branch.port_area >= 0.0 ? branch.port_area
                  : branch.port_area > -1.5 ? velum_now_
                  : branch.port_area > -3.0 ? under_tongue_now_
                                            : lateral_pocket_now_;
    double into_child = 0.0;
    double into_owner = 0.0;
    if (branch.shadow_port >= 0) {
      // The two nodes as two ducts in parallel into the branch: each its
      // share of the area; the branch hears their waves weighted by area.
      branch.shadow_share *= std::exp(-kStepSeconds / settings_.branch_glide);
      if (branch.shadow_share < 1.0e-4) {
        branch.shadow_share = 0.0;
      }
      double main_area = area * (1.0 - branch.shadow_share);
      double shadow_area = area * branch.shadow_share;
      double main_child = 0.0;
      double main_owner = 0.0;
      double shadow_child = 0.0;
      double shadow_owner = 0.0;
      if (main_area >= 1.0e-6) {
        series_element(owner.side_outgoing(branch.port), child.outgoing_left(), kRhoC / main_area,
                       viscous_resistance_per_cm(main_area) * kNodeLength, 0.0, owner.area(branch.node), child.area(0),
                       main_child, main_owner);
      }
      if (shadow_area >= 1.0e-6) {
        series_element(owner.side_outgoing(branch.shadow_port), child.outgoing_left(), kRhoC / shadow_area,
                       viscous_resistance_per_cm(shadow_area) * kNodeLength, 0.0, owner.area(branch.shadow_node),
                       child.area(0), shadow_child, shadow_owner);
      }
      double open = area >= 1.0e-6 ? area : 0.0;
      owner.set_side_port(branch.port, main_owner, main_area >= 1.0e-6 ? main_area : 0.0);
      owner.set_side_port(branch.shadow_port, shadow_owner, shadow_area >= 1.0e-6 ? shadow_area : 0.0);
      double heard = open > 0.0 ? (main_area * main_child + shadow_area * shadow_child) / open : 0.0;
      child.set_left_port(heard, open);
      continue;
    }
    if (area >= 1.0e-6) {
      double scale = branch.tube == nose_ ? settings_.port_resistance_scale : 1.0;
      series_element(owner.side_outgoing(branch.port), child.outgoing_left(), kRhoC / area,
                     scale * viscous_resistance_per_cm(area) * kNodeLength, 0.0, owner.area(branch.node), child.area(0),
                     into_child, into_owner);
    }
    double open = area >= 1.0e-6 ? area : 0.0;
    owner.set_side_port(branch.port, into_owner, open);
    child.set_left_port(into_child, open);
  }
}

// Unit-variance white noise (sum of four uniforms) for the tremor.
double VoiceBody::tremor_white() {
  double sum = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    tremor_random_ ^= tremor_random_ << 13U;
    tremor_random_ ^= tremor_random_ >> 17U;
    tremor_random_ ^= tremor_random_ << 5U;
    sum += static_cast<double>(tremor_random_) / 4294967296.0;
  }
  return (sum - 2.0) * std::sqrt(3.0);
}

// Physiological tremor of the laryngeal muscles: motor units fire
// asynchronously, so the held force carries unit-rms noise in about 2-15 Hz.
// Two poles at each corner: with one, the band kept a 6 dB/oct skirt to
// 500 Hz, and the adductors' share of it moved the folds' rest gap from one
// glottal cycle to the next: cycle-to-cycle jitter (RAP) 0.82% where people
// and VocalTractLab's synthesis show 0.23-0.6% (0.32% without the tremor).
double VoiceBody::tremor_band(int channel) {
  double white = tremor_white();
  double slow_gain = 1.0 - std::exp(-kTwoPi * 2.0 * kControlSeconds);
  double fast_gain = 1.0 - std::exp(-kTwoPi * 15.0 * kControlSeconds);
  tremor_slow_[channel] += slow_gain * (white - tremor_slow_[channel]);
  tremor_slower_[channel] += slow_gain * (tremor_slow_[channel] - tremor_slower_[channel]);
  tremor_fast_[channel] += fast_gain * (white - tremor_fast_[channel]);
  tremor_faster_[channel] += fast_gain * (tremor_fast_[channel] - tremor_faster_[channel]);
  return tremor_scale_ * (tremor_faster_[channel] - tremor_slower_[channel]);
}

// Slow causes, once per millisecond: breath, the stretching muscle, the
// arytenoid hold, and every wall's support effort (shared plus local).
void VoiceBody::control_tick() {
  variation_.advance(kControlSeconds);
  if (settings_.ear_latency > 0.0) {
    // The listening ear's error reaches the muscle ear_latency late, and
    // only while there is voice to hear (the lungs pressing).
    std::size_t length = ear_delay_.size();
    double arriving = ear_delay_[ear_slot_];
    ear_delay_[ear_slot_] = lung_ > 2500.0 ? heard_error_ : 0.0;
    ear_slot_ = (ear_slot_ + 1) % length;
    pitch_correction_ += settings_.ear_rate * kControlSeconds * arriving;
    pitch_correction_ = std::clamp(pitch_correction_, -settings_.pitch_correction_limit, settings_.pitch_correction_limit);
  }
  for (int channel = 0; channel < kTremorChannels; ++channel) {
    tremor_value_[channel] = tremor_band(channel);
  }
  double shared = variation_.value(3);
  for (Index node = 0; node < vocal_.count(); ++node) {
    double local = variation_.value(kFirstWallChannel + node);
    double mix = settings_.wall_shared * shared + (1.0 - settings_.wall_shared) * local;
    double release = node >= release_first_ && node <= release_last_ ? 1.0 - release_share_ : 1.0;
    // Past the lips the flared open end is air, not tissue: rigid. (As soft
    // wall up to 30 cm^2, its huge compliance rang near 40 Hz at every change
    // of the flow through the lips: thumps, "tapping the microphone".)
    double outside = lip_end_ > 0 && node > lip_end_ + 1 ? 1000.0 : 1.0;
    vocal_.set_support_factor(node, outside * release * (1.0 + settings_.wall_variation * mix));
    vocal_.set_support_offset(node, node >= pull_first_ && node <= pull_last_ ? pull_offset_ : 0.0);
  }
  vocal_.update_rest_geometry(static_cast<int>(kControlSteps));
  place_jet();
  // The muscles' causes for the next millisecond; each step interpolates
  // from the last tick's values to these (a staircase at 1 kHz, even
  // glided over 2 ms, left a 1 kHz comb 15 dB above the floor in pauses).
  for (int cause = 0; cause < kCauses; ++cause) {
    cause_from_[cause] = cause_to_[cause];
  }
  cause_to_[0] = std::exp(pitch_correction_) *
                 (1.0 + settings_.tension_variation * variation_.value(1) + settings_.tremor * tremor_value_[0]);
  cause_to_[1] = settings_.posture_variation * variation_.value(2) + settings_.adductor_tremor * tremor_value_[1];
  cause_to_[2] = 1.0 + settings_.breath_variation * variation_.value(0) + settings_.breath_tremor * tremor_value_[2];
}

// Turbulence above the glottis. Strength after VocalTractLab 2.4: amplitude
// gain v^3 sqrt(A), low-passed at 500 Hz, rising with a 3 ms time constant on
// its square root and falling at once. The jet speed is the instantaneous
// flow over the instantaneous area, capped by what the pressure across the
// glottis allows. (Dividing the 500 Hz low-passed flow by the instantaneous
// area, as read, gave 500,000 cm/s during closure: HNR -2 dB, shimmer 39%.)
// Spectrum: the same jet turbulence as in the mouth, centred near a Strouhal
// number of 0.2 on the slit height (the glottis is a slit as long as the
// folds), rising 12 dB/oct below a third of that corner and falling
// 12 dB/oct above it; total power as VocalTractLab's. (VocalTractLab's
// 1 kHz low-pass left /h/ and aspiration 7-10 dB weak above 2 kHz against
// the front end's bdl-tuned source.)
double VoiceBody::aspiration_sample() {
  double area = std::max(cord_.area(), kMinimumArea);
  double drop = std::abs(trachea_.pressure(trachea_.count() - 1) - vocal_.pressure(0));
  double velocity = std::min(std::max(cord_.flow(), 0.0) / area, std::sqrt(2.0 * drop / kAirDensity));
  double instantaneous = 0.2e-6 * velocity * velocity * velocity * std::sqrt(area);
  double flow_gain = 1.0 - std::exp(-kTwoPi * 500.0 * kStepSeconds);
  aspiration_flow_ += flow_gain * (instantaneous - aspiration_flow_);
  double root = std::sqrt(std::max(aspiration_flow_, 0.0));
  if (root > aspiration_root_) {
    aspiration_root_ += (root - aspiration_root_) * (1.0 - std::exp(-kStepSeconds / 0.003));
  } else {
    aspiration_root_ = root;
  }
  // Corner from the jet's width: the slit height (area over the fold length)
  // or, when `glottal_round_jet`, its hydraulic diameter (the jet that
  // reaches the false folds and epiglottis has merged into a round jet).
  double slit_height = settings_.glottal_round_jet ? std::sqrt(4.0 * area / kPi) : std::max(area / 1.3, 0.005);
  double corner = std::clamp(settings_.strouhal * velocity / slit_height, 500.0, 12000.0);
  double corner_gain = 1.0 - std::exp(-kTwoPi * corner / kStepRate);
  double floor_gain = 1.0 - std::exp(-kTwoPi * corner / (settings_.frication_band * kStepRate));
  double pole = 1.0 - corner_gain;
  double gain_squared = corner_gain * corner_gain;
  double pole_term = 1.0 - pole * pole;
  double passed = gain_squared * gain_squared * (1.0 + pole * pole) / (pole_term * pole_term * pole_term);
  double white = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    noise_random_ ^= noise_random_ << 13U;
    noise_random_ ^= noise_random_ >> 17U;
    noise_random_ ^= noise_random_ << 5U;
    white += static_cast<double>(noise_random_) / 4294967296.0;
  }
  white = (white - 2.0) * std::sqrt(3.0);
  noise_low_ += floor_gain * (white - noise_low_);
  double once = white - noise_low_;
  noise_lower_ += floor_gain * (once - noise_lower_);
  noise_z1_ += corner_gain * ((once - noise_lower_) - noise_z1_);
  noise_z2_ += corner_gain * (noise_z1_ - noise_z2_);
  // VocalTractLab's source: noise of rms 0.29 through a 1 kHz second-order
  // low-pass at 48 kHz (about 0.152 of its amplitude), times 0.6.
  constexpr double kSourceRms = 0.2887 * 0.152 * 0.6;
  double unit = noise_z2_ / std::sqrt(passed);
  // Closing-phase instability (smoothed at 2 kHz): the deceleration over
  // 1 ms relative to the pulse's peak flow (held with a 10 ms decay), so it
  // measures how sharply each pulse is cut off. (Relative to the
  // instantaneous flow it grew without bound at closure: /i/ +22 dB.)
  double flow_now = cord_.flow();
  double rate = (flow_now - glottal_flow_previous_) * kStepRate;
  glottal_flow_previous_ = flow_now;
  deceleration_ += (1.0 - std::exp(-kTwoPi * 2000.0 * kStepSeconds)) * (std::max(-rate, 0.0) - deceleration_);
  flow_peak_ = std::max(flow_now, flow_peak_ * std::exp(-kStepSeconds / 0.01));
  double instability = std::min(1.0e-3 * deceleration_ / std::max(flow_peak_, 50.0), 3.0);
  double gain = 1.0 + settings_.closing_turbulence * instability;
  // The same turbulence as a fluctuation of the flux: rms share jet_flux of
  // the flow once the jet is turbulent (Reynolds number on the slit height
  // past about 1200).
  double reynolds = velocity * 2.0 * slit_height * kAirDensity / kAirViscosity;
  double onset = 1200.0 / std::max(reynolds, 1.0);
  double turbulent = std::max(1.0 - onset * onset, 0.0);
  double flux_unit = unit;
  if (settings_.flux_corner > 0.0) {
    double flux_gain = 1.0 - std::exp(-kTwoPi * settings_.flux_corner / kStepRate);
    flux_z1_ += flux_gain * (white - flux_z1_);
    flux_z2_ += flux_gain * (flux_z1_ - flux_z2_);
    flux_unit = flux_z2_;
  }
  // The flux belongs to the voice: it plays while the folds meet every
  // cycle (the glottal area's running minimum under 0.02 cm^2) and fades out
  // as they stay apart (0.06 cm^2; the aspiration source speaks for an open
  // glottis). (Following the whole flow, every voiceless consonant, whose
  // open glottis passes flow continuously, carried a loud low hiss: /hh/
  // 11.7 dB from bdl's bands, 3.7 without the flux; /f s sh/ +10-28 dB at
  // 0.5-2 kHz. Capping the flow at a voiced level changed nothing.)
  double area_now = cord_.area();
  area_floor_ = area_now < area_floor_ ? area_now
                                       : area_floor_ + (1.0 - std::exp(-kStepSeconds / 0.015)) * (area_now - area_floor_);
  double voiced_share = std::clamp((0.06 - area_floor_) / 0.04, 0.0, 1.0);
  cord_.set_flux_noise(settings_.jet_flux * voiced_share * turbulent * std::max(cord_.flow(), 0.0) * flux_unit);
  return settings_.aspiration * gain * aspiration_root_ * aspiration_root_ * kSourceRms * unit;
}

// VocalTractLab 2.4's aspiration (TdsModel.cpp, GPL: mathematics only): the
// glottal outflow through one pole at 500 Hz, over the upper glottal
// section's area, the jet speed v (outflow only); amplitude 0.2e-6 v^3
// sqrt(A) (Stevens), tapered to 0 from 0.5 to 1 cm^2; rising on its square
// root with a 3 ms time constant, falling at once, silent under 0.001; its
// noise of rms 1/sqrt(12) within +-1 (here a sum of four uniforms, which is
// that) through a critically damped second-order low-pass at 1 kHz, times
// 0.6 (and (1000 / f_c)^1.5 = 1). Their noise runs at 48 kHz: ours, at
// twice the rate, is raised by sqrt(2) to keep the same density in band.
double VoiceBody::vtl_aspiration_sample() {
  double pole = std::exp(-kTwoPi * 500.0 * kStepSeconds);
  vtl_noise_flow_ = (1.0 - pole) * cord_.flow() + pole * vtl_noise_flow_;
  double area = std::max(cord_.upper_area(), 1.0e-4);
  double velocity = std::max(vtl_noise_flow_, 0.0) / area;
  double gain = 0.2e-6 * settings_.vtl_aspiration_gain;
  if (area >= 0.5) {
    gain *= std::max(1.0 - (area - 0.5) / 0.5, 0.0);
  }
  double target = gain * velocity * velocity * velocity * std::sqrt(area);
  if (target > vtl_noise_amp_) {
    double hold = std::exp(-1.0 / (0.003 * kStepRate));
    double root = (1.0 - hold) * std::sqrt(target) + hold * std::sqrt(vtl_noise_amp_);
    vtl_noise_amp_ = root * root;
  } else {
    vtl_noise_amp_ = target;
  }
  if (vtl_noise_amp_ < 0.001) {
    vtl_noise_amp_ = 0.0;
    vtl_noise_state_[0] = 0.0;
    vtl_noise_state_[1] = 0.0;
    return 0.0;
  }
  double sum = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    noise_random_ ^= noise_random_ << 13U;
    noise_random_ ^= noise_random_ >> 17U;
    noise_random_ ^= noise_random_ << 5U;
    sum += static_cast<double>(noise_random_) / 4294967296.0;
  }
  double white = 0.5 * (sum - 2.0) * std::sqrt(kStepRate / 48000.0);
  // Second-order Butterworth low-pass at 1 kHz (bilinear, prewarped).
  double corner = std::clamp(settings_.vtl_aspiration_corner, 100.0, 0.45 * kStepRate);
  double k = std::tan(kPi * corner * kStepSeconds);
  double norm = 1.0 / (1.0 + std::sqrt(2.0) * k + k * k);
  double b0 = k * k * norm;
  double a1 = 2.0 * (k * k - 1.0) * norm;
  double a2 = (1.0 - std::sqrt(2.0) * k + k * k) * norm;
  double out = b0 * white + vtl_noise_state_[0];
  vtl_noise_state_[0] = 2.0 * b0 * white - a1 * out + vtl_noise_state_[1];
  vtl_noise_state_[1] = b0 * white - a2 * out;
  double factor = 1000.0 / corner;
  return out * vtl_noise_amp_ * factor * std::sqrt(factor) * 0.6;
}

// Where does a jet form? At the narrowest open section of the mouth (rest
// geometry, past the pharynx, narrower than 0.6 cm^2). It leaves the
// constriction laminar, turns turbulent over the next few millimetres and
// decays over about a centimetre; where it still holds together at the
// incisors it strikes them, several times stronger (an obstacle source).
// Ported from the front end's body, where these weights were tuned against
// bdl and the owner's sibilants and bursts.
void VoiceBody::place_jet() {
  // (Up to the lips: past them, a flaring body's open end is outside air.)
  Index count = lip_end_ > 0 ? std::min(lip_end_ + 1, vocal_.active_count()) : vocal_.active_count();
  Index narrowest = -1;
  double narrowest_area = 0.6;
  // (Search above the epilarynx, about 3.5 cm: its 0.45 cm^2 tube is not a
  // jet-forming constriction, and taken as one it hid every oral place
  // wider than itself.)
  for (Index node = kFirstOralNode; node < count; ++node) {
    double area = vocal_.rest_area(node);
    if (area < narrowest_area && area > 0.0) {
      narrowest_area = area;
      narrowest = node;
    }
  }
  jet_narrowest_ = narrowest;
  // How open the mouth is, for the velum's leak: shut at 0.05 cm^2, open
  // from 0.35 (no jet-forming constriction: open).
  oral_open_ = narrowest < 0 ? 1.0 : std::clamp((narrowest_area - 0.05) / 0.3, 0.0, 1.0);
  std::fill(jet_weights_.begin(), jet_weights_.end(), 0.0);
  std::fill(jet_obstacle_weights_.begin(), jet_obstacle_weights_.end(), 0.0);
  if (narrowest < 0) {
    return;
  }
  double exit_x = (narrowest + 1) * kNodeLength;
  for (Index link = narrowest; link <= count - 2; ++link) {
    double x = (link + 1) * kNodeLength;
    double distance = std::max(x - exit_x, 0.0);
    double teeth_cover = std::clamp((x - teeth_start_) / kNodeLength + 0.5, 0.0, 1.0) *
                         std::clamp((teeth_end_ - x) / kNodeLength + 0.5, 0.0, 1.0);
    double coherence = std::exp(-std::max(distance - settings_.obstacle_reach, 0.0) / 0.5);
    double growth = 1.0 - std::exp(-distance / settings_.jet_growth);
    double along = 0.35 * growth * std::exp(-distance / settings_.jet_decay);
    // Two sources: the free jet's turbulence (wide by now: jet_spread) and
    // the still-coherent jet striking the incisors (an obstacle source, as
    // narrow as the slit). A sibilant is mostly the second, a velar or
    // labial release only the first.
    jet_weights_[static_cast<std::size_t>(link)] = along;
    jet_obstacle_weights_[static_cast<std::size_t>(link)] = along * 3.0 * teeth_cover * coherence;
  }
  if (narrowest >= count - 3) {
    // Labiodental: the jet leaves between the upper incisors and the lip.
    jet_obstacle_weights_[static_cast<std::size_t>(count - 2)] += settings_.labiodental;
  }
}

double VoiceBody::jet_gauss() {
  double sum = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    jet_random_ ^= jet_random_ << 13U;
    jet_random_ ^= jet_random_ >> 17U;
    jet_random_ ^= jet_random_ << 5U;
    sum += static_cast<double>(jet_random_) / 4294967296.0;
  }
  return (sum - 2.0) * std::sqrt(3.0);
}

// The jet's turbulence, each step: dipole strength about v^3 sqrt(A) once
// the jet is turbulent (Reynolds number past about 1800 on its diameter),
// broadband with its energy near a Strouhal number of 0.2 on the slit
// height (oral constrictions are slits about 1.2 cm wide): rising
// 12 dB/oct below a third of that corner, falling 12 dB/oct above it. The
// cavities in front of the constriction then shape it.
void VoiceBody::sound_jet() {
  if (jet_narrowest_ < 0) {
    return;
  }
  Index count = lip_end_ > 0 ? std::min(lip_end_ + 1, vocal_.active_count()) : vocal_.active_count();
  Index narrowest = jet_narrowest_;
  double narrowest_area = vocal_.area(narrowest);
  for (Index node = std::max(narrowest - 3, kFirstOralNode); node <= std::min(narrowest + 3, count - 1);
       ++node) {
    if (vocal_.area(node) < narrowest_area) {
      narrowest_area = vocal_.area(node);
      narrowest = node;
    }
  }
  double area = std::max(narrowest_area, kMinimumArea);
  double flow = vocal_.link_flow(std::min(narrowest, count - 2));
  Index upstream = std::max(narrowest - 2, static_cast<Index>(0));
  Index downstream = std::min(narrowest + 2, count - 2);
  // At the lips the jet exits into the room, at atmospheric pressure.
  double outside = narrowest >= count - 3 ? 0.0 : vocal_.pressure(downstream);
  double drop = std::abs(vocal_.pressure(upstream) - outside) + 1.0;
  double velocity = std::min(std::abs(flow) / area, std::sqrt(2.0 * drop / kAirDensity));
  double diameter = std::sqrt(4.0 * area / kPi);
  double reynolds = velocity * diameter * kAirDensity / kAirViscosity;
  double onset = 1800.0 / std::max(reynolds, 1.0);
  double turbulent = std::max(1.0 - onset * onset, 0.0);
  double fade = std::clamp((0.6 - area) / 0.3, 0.0, 1.0);
  double dynamic = 0.5 * kAirDensity * velocity * velocity * turbulent * fade * fade * (velocity / 2000.0);
  // Free-jet turbulence forms a few millimetres downstream, where the jet
  // has widened (by about 0.2 of the distance): that width sets its corner.
  // At the incisors the jet is still about as narrow as the slit.
  // Tongue constrictions are slits about 1.2 cm wide; the lips (and the lip
  // against the upper incisors in /f v/) about 2.5 cm. (At 1.2 cm /f/ came
  // out low-heavy, 0/-5/-10/-14 dB in 0.5-8 kHz bands where bdl's rises to
  // 8 kHz, -17/-12/-8/0.)
  double slit_width = narrowest >= count - 3 ? settings_.lip_slit_width : 1.2;
  double slit = std::max(area / slit_width, 0.01);
  // With no spread the two are one coherent source (as fitted on the
  // front end's body; split, /s/ lost 4 dB and its 3-5 kHz peak).
  double free_noise = shaped_jet_noise(jet_free_, velocity, slit + settings_.jet_spread);
  double obstacle_noise = settings_.jet_spread > 0.0 ? shaped_jet_noise(jet_obstacle_, velocity, slit) : free_noise;
  double strength = 0.2 * settings_.frication * dynamic * std::sqrt(area);
  // Turbulence takes a few milliseconds to develop as a jet forms: the
  // noise builds with jet_attack and falls at once. (Following the velocity
  // instantly, every release rose 41-47 dB in 5 ms against bdl's 29-34, and
  // peaked 5-8 dB over bdl's bursts.)
  double attack_gain = 1.0 - std::exp(-kStepSeconds / std::max(settings_.jet_attack, kStepSeconds));
  if (strength > jet_strength_) {
    jet_strength_ += attack_gain * (strength - jet_strength_);
  } else {
    jet_strength_ = strength;
  }
  strength = jet_strength_;
  for (Index link = std::min(narrowest, count - 2); link <= count - 2; ++link) {
    double free_weight = jet_weights_[static_cast<std::size_t>(link)];
    double obstacle_weight = jet_obstacle_weights_[static_cast<std::size_t>(link)];
    if (free_weight > 0.0 || obstacle_weight > 0.0) {
      vocal_.add_series_source(link, strength * (free_weight * free_noise + obstacle_weight * obstacle_noise));
    }
  }
}

namespace {

// SplitMix64: a hash, so the reversible jet's noise for any link and step
// can be drawn again in any order (the reversible core draws it the same).
std::uint64_t split_mix(std::uint64_t key) {
  key += 0x9E3779B97F4A7C15ULL;
  key = (key ^ (key >> 30U)) * 0xBF58476D1CE4E5B9ULL;
  key = (key ^ (key >> 27U)) * 0x94D049BB133111EBULL;
  return key ^ (key >> 31U);
}

// Unit-variance noise for (link, step): four hashed uniforms summed.
double hashed_gauss(std::uint64_t seed, long link, long step) {
  double sum = 0.0;
  for (std::uint64_t draw = 0; draw < 4U; ++draw) {
    std::uint64_t key = seed ^ ((static_cast<std::uint64_t>(step) * 4096U + static_cast<std::uint64_t>(link)) * 4U + draw);
    sum += static_cast<double>(split_mix(key) >> 11U) * (1.0 / 9007199254740992.0);
  }
  return (sum - 2.0) * std::sqrt(3.0);
}

}  // namespace

// The reversible jet: from the narrowest node (area under 0.6 cm^2, above
// the epilarynx) to the open end, each link a turbulence source of strength
// gain * (rho/2) v^2 * turbulent * fade^2 * (v / 2000) * sqrt(A) * weight,
// v = |U_link| / A from the link's own flow of the step before (the jet
// passing it) and A the narrowest area; weight exp(-distance / 1 cm) from
// the narrowest node's exit. The noise is hashed per link and step and
// shaped per link as shaped_jet_noise shapes it. Everything that drives a
// source is known to a peel that has come from the lips to that link.
void VoiceBody::reversible_jet() {
  Index count = vocal_.active_count();
  if (link_jets_.size() < static_cast<std::size_t>(vocal_.count())) {
    link_jets_.resize(static_cast<std::size_t>(vocal_.count()));
  }
  long step = jet_step_++;
  Index narrowest = -1;
  double narrowest_area = 0.6;
  for (Index node = kFirstOralNode; node < count; ++node) {
    if (vocal_.area(node) < narrowest_area) {
      narrowest_area = vocal_.area(node);
      narrowest = node;
    }
  }
  double area = std::max(narrowest_area, kMinimumArea);
  double diameter = std::sqrt(4.0 * area / kPi);
  double fade = std::clamp((0.6 - area) / 0.3, 0.0, 1.0);
  if (settings_.reversible_jet_closing > 0.0) {
    double closed = 0.25 * settings_.reversible_jet_closing;
    fade *= std::clamp((area - closed) / (settings_.reversible_jet_closing - closed), 0.0, 1.0);
  }
  double slit = std::max(area / 1.2, 0.01);
  for (Index link = 0; link + 1 < count; ++link) {
    double velocity = std::abs(vocal_.link_flow(link)) / area;
    JetNoise& state = link_jets_[static_cast<std::size_t>(link)];
    double corner = std::clamp(settings_.strouhal * velocity / slit, 1500.0, settings_.jet_corner_ceiling);
    double corner_gain = 1.0 - std::exp(-kTwoPi * corner / kStepRate);
    double floor_gain = 1.0 - std::exp(-kTwoPi * corner / (settings_.frication_band * kStepRate));
    double pole = 1.0 - corner_gain;
    double gain_squared = corner_gain * corner_gain;
    double pole_term = 1.0 - pole * pole;
    double passed = gain_squared * gain_squared * (1.0 + pole * pole) / (pole_term * pole_term * pole_term);
    double white = hashed_gauss(settings_.seed, link, step);
    state.low += floor_gain * (white - state.low);
    double once = white - state.low;
    state.lower += floor_gain * (once - state.lower);
    state.band += corner_gain * ((once - state.lower) - state.band);
    state.second += corner_gain * (state.band - state.second);
    if (narrowest < 0 || link < narrowest) {
      continue;
    }
    double reynolds = velocity * diameter * kAirDensity / kAirViscosity;
    double onset = 1800.0 / std::max(reynolds, 1.0);
    double turbulent = std::max(1.0 - onset * onset, 0.0);
    double dynamic = 0.5 * kAirDensity * velocity * velocity * turbulent * fade * fade * (velocity / 2000.0);
    double distance = static_cast<double>(link - narrowest) * kNodeLength;
    double strength = settings_.reversible_jet * jet_gate_ * dynamic * std::sqrt(area) * std::exp(-distance / 1.0);
    vocal_.add_series_source(link, strength * state.second / std::sqrt(passed));
  }
}

// One jet-noise generator of unit variance: white noise with the band below
// a third of the corner removed (12 dB/oct) and two low-pass stages at the
// corner, the corner a Strouhal number of 0.2 on the jet's width.
double VoiceBody::shaped_jet_noise(JetNoise& state, double velocity, double width) {
  double corner = std::clamp(settings_.strouhal * velocity / width, 1500.0, settings_.jet_corner_ceiling);
  double corner_gain = 1.0 - std::exp(-kTwoPi * corner / kStepRate);
  double floor_gain = 1.0 - std::exp(-kTwoPi * corner / (settings_.frication_band * kStepRate));
  // Two identical poles a = 1 - g pass g^4 (1 + a^2) / (1 - a^2)^3 of white
  // noise's power.
  double pole = 1.0 - corner_gain;
  double gain_squared = corner_gain * corner_gain;
  double pole_term = 1.0 - pole * pole;
  double passed = gain_squared * gain_squared * (1.0 + pole * pole) / (pole_term * pole_term * pole_term);
  double white = jet_gauss();
  state.low += floor_gain * (white - state.low);
  double once = white - state.low;
  state.lower += floor_gain * (once - state.lower);
  state.band += corner_gain * ((once - state.lower) - state.band);
  state.second += corner_gain * (state.band - state.second);
  return state.second / std::sqrt(passed);
}

// The ear: at each closure the heard period trims the stretching muscle.
// Only cycles of a voice: lungs over 250 Pa and a period within half an
// octave of the intended one (folds flickering shut for a breath, or the
// first cycle after a stop, are not heard as pitch).
void VoiceBody::hear_pitch() {
  bool closed = cord_.contact() > 0.0;
  if (closed && !was_closed_) {
    if (last_closure_ >= 0) {
      double heard = kStepRate / static_cast<double>(tick_ - last_closure_);
      double ratio = heard / intended_pitch_;
      bool voiced = lung_ > 2500.0 && ratio > 0.707 && ratio < 1.414;
      if (voiced && settings_.ear_latency > 0.0) {
        // The listening ear: the error, smoothed by the period's share of
        // the smoothing time (applied at the control ticks, delayed).
        double period = static_cast<double>(tick_ - last_closure_) / kStepRate;
        double gain = 1.0 - std::exp(-period / std::max(settings_.ear_smoothing, 1.0e-3));
        heard_error_ += gain * (std::log(intended_pitch_ / heard) - heard_error_);
      } else if (voiced) {
        pitch_correction_ += settings_.pitch_feedback * std::log(intended_pitch_ / heard);
        pitch_correction_ = std::clamp(pitch_correction_, -settings_.pitch_correction_limit,
                                       settings_.pitch_correction_limit);
      }
    }
    last_closure_ = tick_;
  }
  was_closed_ = closed;
}

void VoiceBody::set_lung_pressure(double pascal, double move_seconds) {
  // The command is approached along a straight ramp over the next 1 ms
  // (commands arrive at the control rate; a staircase left a 1 kHz comb
  // in the breath noise of an inhale), then by the first-order follow.
  lung_ramp_from_ = lung_ramp_;
  lung_target_ = pascal * kPascal;
  lung_ramp_left_ = static_cast<int>(kControlSteps);
  double tau = std::max(move_seconds / 3.0, kStepSeconds);
  lung_gain_ = 1.0 - std::exp(-kStepSeconds / tau);
}

double VoiceBody::step() {
  if (tick_ % kControlSteps == 0) {
    control_tick();
  }
  double share = static_cast<double>(tick_ % kControlSteps) / static_cast<double>(kControlSteps);
  ++tick_;
  if (lung_ramp_left_ > 0) {
    lung_ramp_ += (lung_target_ - lung_ramp_from_) / static_cast<double>(kControlSteps);
    --lung_ramp_left_;
  } else {
    lung_ramp_ = lung_target_;
  }
  lung_ += lung_gain_ * (lung_ramp_ - lung_);
  cord_.set_holding(cause_from_[0] + share * (cause_to_[0] - cause_from_[0]),
                    cause_from_[1] + share * (cause_to_[1] - cause_from_[1]));
  // Breath variation and tremor, interpolated and glided (2 ms).
  double breath_target = cause_from_[2] + share * (cause_to_[2] - cause_from_[2]);
  breath_factor_ += (1.0 - std::exp(-kStepSeconds / 0.002)) * (breath_target - breath_factor_);
  double breath = lung_ * breath_factor_;
  if (settings_.vtl_aspiration > 0.5) {
    // VocalTractLab's: a dipole 1.5 cm above the glottis, shared between
    // the two sections it falls between; no other glottal noise.
    double at = 1.5 / kNodeLength;
    Index below = std::min(static_cast<Index>(at), vocal_.count() - 2);
    double share = std::clamp(at - static_cast<double>(below), 0.0, 1.0);
    double turbulence = vtl_aspiration_sample();
    vocal_.add_series_source(below, (1.0 - share) * turbulence);
    vocal_.add_series_source(below + 1, share * turbulence);
    cord_.set_flux_noise(0.0);
    cord_.set_surface_pressure(0.0);
  } else {
    double turbulence = aspiration_sample();
    vocal_.add_series_source(aspiration_node_, turbulence);
    cord_.set_surface_pressure(settings_.fold_turbulence * turbulence);
  }
  sound_jet();
  if (settings_.reversible_jet > 0.0) {
    reversible_jet();
  }
  // The lungs: alveolar pressure behind the bronchial tree. The tree's
  // branching absorbs sound (a matched end, no reflection) yet passes steady
  // flow almost freely (airway resistance about 1.5 dyn s/cm^5). As a plain
  // matched port the steady flow met rho c / A = 16 dyn s/cm^5, and an open
  // glottis (200 cm^3/s) pulled the pressure below the folds from 7500 to
  // 2300: weak bursts and weak vowels after every voiceless consonant. The
  // source is therefore raised by (Z - R_airway) times the flow below 20 Hz.
  double port_area = trachea_.area(0) * settings_.lung_port_scale;
  double impedance = kRhoC / port_area;
  double lung_flow = (0.5 * breath - trachea_.outgoing_left()) / impedance;  // previous step's flow into the trachea
  lung_flow_slow_ += (1.0 - std::exp(-kTwoPi * 20.0 * kStepSeconds)) * (lung_flow - lung_flow_slow_);
  // (Only while the port's impedance exceeds the airway's resistance: with a
  // wide, open lung end (subglottal.lung_port_scale) it fell to 0.35, the
  // correction turned negative with a loop gain above 1, and the pressure
  // below the folds drained until they stopped.)
  double compensation = std::max(impedance - kAirwayResistance, 0.0) * lung_flow_slow_;
  trachea_.set_left_port(0.5 * (breath + compensation), port_area);
  trachea_.set_right_port(from_cord_below_, trachea_.area(trachea_.count() - 1));
  vocal_.set_left_port(from_cord_above_, vocal_.area(0));
  exchange_branches();
  trachea_.step();
  vocal_.step();
  velum_quiet_ = velum_now_ > 0.0 ? 0 : velum_quiet_ + 1;
  under_quiet_ = under_tongue_now_ > 0.0 ? 0 : under_quiet_ + 1;
  lateral_quiet_ = lateral_pocket_now_ > 0.0 ? 0 : lateral_quiet_ + 1;
  constexpr long kQuietSteps = static_cast<long>(0.05 * kStepRate);
  for (std::size_t index = 0; index < branch_tubes_.size(); ++index) {
    int group = tube_group_[index];
    if ((group == 1 && velum_quiet_ > kQuietSteps) || (group == 2 && under_quiet_ > kQuietSteps) ||
        (group == 4 && lateral_quiet_ > kQuietSteps)) {
      continue;
    }
    branch_tubes_[index].step();
  }
  cord_.set_expected_breath(lung_ / kPascal);
  double returning = vocal_.outgoing_left();
  if (forcing_) {
    returning = forced_wave_;
    forced_flow_ = (from_cord_above_ - returning) * vocal_.area(0) / kRhoC;
    if (observe_gain_ > 0.0 && std::isfinite(observed_flow_)) {
      returning += observe_gain_ * kRhoC / vocal_.area(0) * (forced_flow_ - observed_flow_);
    }
  }
  cord_.step(trachea_.outgoing_right(), trachea_.area(trachea_.count() - 1), returning, vocal_.area(0),
             from_cord_below_, from_cord_above_);
  hear_pitch();
  // Sound: the rate of change of the flow leaving the lips, plus the
  // lining's own motion heard through the skin.
  // The lining's motion reaches the air through the neck and cheeks, which
  // pass only the lows: two poles at skin_corner. (Passed whole, a sealed
  // tract rang through the skin at every 1 kHz, its closed-closed tube's
  // resonances, all through /b/ closures: the owner heard "body" as
  // "bally" with a ringing. A real voice bar sits under 500 Hz.)
  skin_low_ += skin_pole_gain_ * (vocal_.total_wall_flow() - skin_low_);
  skin_lower_ += skin_pole_gain_ * (skin_low_ - skin_lower_);
  double flow = vocal_.radiated_flow() + branch_tubes_[static_cast<std::size_t>(nose_)].radiated_flow() +
                skin_gain_ * skin_lower_;
  double sound = (flow - previous_flow_) * kStepRate;
  previous_flow_ = flow;
  if (settings_.warmth > 0.0) {
    // The swing itself: the flow less its steady stream (two one-pole
    // high-passes), times 2 pi f so its fundamental matches the rate of
    // change's.
    double pole = 1.0 - std::exp(-kTwoPi * settings_.warmth_corner * kStepSeconds);
    warm_low_ += pole * (flow - warm_low_);
    double swing = flow - warm_low_;
    warm_lower_ += pole * (swing - warm_lower_);
    swing -= warm_lower_;
    warm_pitch_ += (1.0 - std::exp(-kStepSeconds / 0.02)) * (intended_pitch_ - warm_pitch_);
    sound = (1.0 - settings_.warmth) * sound + settings_.warmth * kTwoPi * std::max(warm_pitch_, 50.0) * swing;
  }
  return sound;
}

}  // namespace tract
