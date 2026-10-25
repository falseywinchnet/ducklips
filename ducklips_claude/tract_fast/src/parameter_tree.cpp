#include "tract/parameter_tree.hpp"

#include <cstdio>
#include <cstdlib>
#include <istream>
#include <ostream>

namespace tract {

namespace {

// Where a declared path lives.
enum class SlotKind { layout, fold, speaker, vocal_tube, nose_tube, tongue, pinch };

struct ParameterSlot {
  const char* path = "";
  SlotKind kind = SlotKind::layout;
  double BodyLayout::*layout_field = nullptr;
  double FoldTubeSettings::*fold_field = nullptr;
  double TongueSettings::*tongue_field = nullptr;
  double PinchSettings::*pinch_field = nullptr;
  void (Speaker::*speaker_set)(double) = nullptr;
  double (Speaker::*speaker_get)() const = nullptr;
};

ParameterSlot layout_slot(const char* path, double BodyLayout::*field) {
  ParameterSlot slot;
  slot.path = path;
  slot.kind = SlotKind::layout;
  slot.layout_field = field;
  return slot;
}

ParameterSlot fold_slot(const char* path, double FoldTubeSettings::*field) {
  ParameterSlot slot;
  slot.path = path;
  slot.kind = SlotKind::fold;
  slot.fold_field = field;
  return slot;
}

ParameterSlot pinch_slot(const char* path, double PinchSettings::*field) {
  ParameterSlot slot;
  slot.path = path;
  slot.kind = SlotKind::pinch;
  slot.pinch_field = field;
  return slot;
}

ParameterSlot tongue_slot(const char* path, double TongueSettings::*field) {
  ParameterSlot slot;
  slot.path = path;
  slot.kind = SlotKind::tongue;
  slot.tongue_field = field;
  return slot;
}

ParameterSlot speaker_slot(const char* path, void (Speaker::*set)(double), double (Speaker::*get)() const) {
  ParameterSlot slot;
  slot.path = path;
  slot.kind = SlotKind::speaker;
  slot.speaker_set = set;
  slot.speaker_get = get;
  return slot;
}

ParameterSlot tube_slot(const char* path, SlotKind kind) {
  ParameterSlot slot;
  slot.path = path;
  slot.kind = kind;
  return slot;
}

// The declared tree, block by block. The blocks group causes, not code:
// what the breath, the folds, the ear, the sources of noise, the walls,
// the side branches and the mouth each do.
const std::vector<ParameterSlot>& declared_slots() {
  static const std::vector<ParameterSlot> slots = {
      speaker_slot("breath.lung_follow", &Speaker::set_lung_follow, &Speaker::lung_follow),
      fold_slot("fold.depth", &FoldTubeSettings::depth),
      fold_slot("fold.cover_mass", &FoldTubeSettings::cover_mass),
      fold_slot("fold.body_mass", &FoldTubeSettings::body_mass),
      fold_slot("fold.cover_stiffness", &FoldTubeSettings::cover_stiffness),
      fold_slot("fold.coupling_stiffness", &FoldTubeSettings::coupling_stiffness),
      fold_slot("fold.support_stiffness", &FoldTubeSettings::support_stiffness),
      fold_slot("fold.vertical_stiffness", &FoldTubeSettings::vertical_stiffness),
      fold_slot("fold.cover_loss", &FoldTubeSettings::cover_loss),
      fold_slot("fold.support_loss", &FoldTubeSettings::support_loss),
      fold_slot("fold.contact_stiffness", &FoldTubeSettings::contact_stiffness),
      fold_slot("fold.contact_softness", &FoldTubeSettings::contact_softness),
      fold_slot("fold.muscle_fluctuation", &FoldTubeSettings::muscle_fluctuation),
      fold_slot("fold.muscle_corner", &FoldTubeSettings::muscle_corner),
      fold_slot("fold.natural_pitch", &FoldTubeSettings::natural_pitch),
      fold_slot("fold.pitch_per_tension", &FoldTubeSettings::pitch_per_tension),
      fold_slot("fold.elongation", &FoldTubeSettings::elongation),
      fold_slot("fold.elongation_reference", &FoldTubeSettings::elongation_reference),
      fold_slot("fold.elongation_release", &FoldTubeSettings::elongation_release),
      fold_slot("fold.elongation_floor", &FoldTubeSettings::elongation_floor),
      fold_slot("fold.swing_per_pascal", &FoldTubeSettings::swing_per_pascal),
      fold_slot("fold.swing_onset", &FoldTubeSettings::swing_onset),
      fold_slot("fold.string_mode", &FoldTubeSettings::string_mode),
      fold_slot("fold.strain_stiffening", &FoldTubeSettings::strain_stiffening),
      fold_slot("fold.vertical_stiffening", &FoldTubeSettings::vertical_stiffening),
      fold_slot("fold.contact_exponent", &FoldTubeSettings::contact_exponent),
      fold_slot("fold.strips", &FoldTubeSettings::strips),
      fold_slot("fold.strip1_stiffness", &FoldTubeSettings::strip1_stiffness),
      fold_slot("fold.strip1_mass", &FoldTubeSettings::strip1_mass),
      fold_slot("fold.strip2_stiffness", &FoldTubeSettings::strip2_stiffness),
      fold_slot("fold.strip2_mass", &FoldTubeSettings::strip2_mass),
      fold_slot("fold.strip3_stiffness", &FoldTubeSettings::strip3_stiffness),
      fold_slot("fold.strip3_mass", &FoldTubeSettings::strip3_mass),
      fold_slot("fold.strip_coupling", &FoldTubeSettings::strip_coupling),
      fold_slot("fold.slice1_mass", &FoldTubeSettings::slice1_mass),
      fold_slot("fold.slice1_elasticity", &FoldTubeSettings::slice1_elasticity),
      fold_slot("fold.slice1_support", &FoldTubeSettings::slice1_support),
      fold_slot("fold.slice2_mass", &FoldTubeSettings::slice2_mass),
      fold_slot("fold.slice2_elasticity", &FoldTubeSettings::slice2_elasticity),
      fold_slot("fold.slice2_support", &FoldTubeSettings::slice2_support),
      fold_slot("fold.slice3_mass", &FoldTubeSettings::slice3_mass),
      fold_slot("fold.slice3_elasticity", &FoldTubeSettings::slice3_elasticity),
      fold_slot("fold.slice3_support", &FoldTubeSettings::slice3_support),
      fold_slot("fold.slice4_mass", &FoldTubeSettings::slice4_mass),
      fold_slot("fold.slice4_elasticity", &FoldTubeSettings::slice4_elasticity),
      fold_slot("fold.slice4_support", &FoldTubeSettings::slice4_support),
      fold_slot("fold.slice5_mass", &FoldTubeSettings::slice5_mass),
      fold_slot("fold.slice5_elasticity", &FoldTubeSettings::slice5_elasticity),
      fold_slot("fold.slice5_support", &FoldTubeSettings::slice5_support),
      fold_slot("fold.slice6_mass", &FoldTubeSettings::slice6_mass),
      fold_slot("fold.slice6_elasticity", &FoldTubeSettings::slice6_elasticity),
      fold_slot("fold.slice6_support", &FoldTubeSettings::slice6_support),
      fold_slot("fold.elasticity_wander", &FoldTubeSettings::elasticity_wander),
      fold_slot("fold.elasticity_corner", &FoldTubeSettings::elasticity_corner),
      fold_slot("fold.elasticity_tamping", &FoldTubeSettings::elasticity_tamping),
      speaker_slot("posture.voicing_lower", &Speaker::set_voicing_lower, &Speaker::voicing_lower),
      speaker_slot("posture.voicing_upper", &Speaker::set_voicing_upper, &Speaker::voicing_upper),
      speaker_slot("posture.voicing_chink", &Speaker::set_voicing_chink, &Speaker::voicing_chink),
      speaker_slot("posture.closure_press", &Speaker::set_closure_press, &Speaker::closure_press),
      speaker_slot("larynx.tongue_pull_tension", &Speaker::set_tongue_pull_tension, &Speaker::tongue_pull_tension),
      speaker_slot("larynx.tongue_pull_adduction", &Speaker::set_tongue_pull_adduction, &Speaker::tongue_pull_adduction),
      speaker_slot("larynx.tongue_pull_raise", &Speaker::set_tongue_pull_raise, &Speaker::tongue_pull_raise),
      layout_slot("ear.pitch_feedback", &BodyLayout::pitch_feedback),
      layout_slot("ear.pitch_correction_limit", &BodyLayout::pitch_correction_limit),
      layout_slot("ear.latency", &BodyLayout::ear_latency),
      layout_slot("ear.smoothing", &BodyLayout::ear_smoothing),
      layout_slot("ear.rate", &BodyLayout::ear_rate),
      layout_slot("variation.breath", &BodyLayout::breath_variation),
      layout_slot("variation.tension", &BodyLayout::tension_variation),
      layout_slot("variation.posture", &BodyLayout::posture_variation),
      layout_slot("variation.wall", &BodyLayout::wall_variation),
      layout_slot("variation.wall_shared", &BodyLayout::wall_shared),
      layout_slot("variation.tremor", &BodyLayout::tremor),
      layout_slot("variation.adductor_tremor", &BodyLayout::adductor_tremor),
      layout_slot("variation.breath_tremor", &BodyLayout::breath_tremor),
      speaker_slot("variation.wander_gain", &Speaker::set_wander_gain, &Speaker::wander_gain),
      speaker_slot("motors.smoothing", &Speaker::set_motor_smoothing, &Speaker::motor_smoothing),
      speaker_slot("motors.move_scale", &Speaker::set_move_scale, &Speaker::move_scale),
      layout_slot("glottal_noise.aspiration", &BodyLayout::aspiration),
      layout_slot("glottal_noise.closing_turbulence", &BodyLayout::closing_turbulence),
      layout_slot("glottal_noise.jet_flux", &BodyLayout::jet_flux),
      layout_slot("glottal_noise.flux_corner", &BodyLayout::flux_corner),
      layout_slot("glottal_noise.fold_turbulence", &BodyLayout::fold_turbulence),
      layout_slot("oral_noise.frication", &BodyLayout::frication),
      layout_slot("oral_noise.jet_attack", &BodyLayout::jet_attack),
      layout_slot("oral_noise.strouhal", &BodyLayout::strouhal),
      layout_slot("oral_noise.frication_band", &BodyLayout::frication_band),
      layout_slot("oral_noise.jet_growth", &BodyLayout::jet_growth),
      layout_slot("oral_noise.jet_spread", &BodyLayout::jet_spread),
      layout_slot("oral_noise.jet_corner_ceiling", &BodyLayout::jet_corner_ceiling),
      layout_slot("oral_noise.lip_slit_width", &BodyLayout::lip_slit_width),
      layout_slot("oral_noise.jet_decay", &BodyLayout::jet_decay),
      layout_slot("oral_noise.obstacle_reach", &BodyLayout::obstacle_reach),
      layout_slot("oral_noise.labiodental", &BodyLayout::labiodental),
      layout_slot("walls.skin_gain", &BodyLayout::skin_gain),
      layout_slot("walls.skin_corner", &BodyLayout::skin_corner),
      layout_slot("recording.highpass", &BodyLayout::recording_highpass),
      layout_slot("recording.presence_db", &BodyLayout::recording_presence_db),
      layout_slot("recording.presence_hz", &BodyLayout::recording_presence_hz),
      layout_slot("recording.presence_q", &BodyLayout::recording_presence_q),
      layout_slot("walls.support_log_glide", &BodyLayout::support_log_glide),
      layout_slot("walls.frame_corner", &BodyLayout::frame_corner),
      layout_slot("branches.glide", &BodyLayout::branch_glide),
      layout_slot("harp.warmth", &BodyLayout::warmth),
      layout_slot("harp.warmth_corner", &BodyLayout::warmth_corner),
      layout_slot("walls.stiffness_scale", &BodyLayout::wall_stiffness_scale),
      layout_slot("walls.loss_scale", &BodyLayout::wall_loss_scale),
      speaker_slot("walls.closure_pull", &Speaker::set_closure_pull, &Speaker::closure_pull),
      speaker_slot("walls.closure_relaxation", &Speaker::set_closure_relaxation, &Speaker::closure_relaxation),
      layout_slot("branches.piriform_loss", &BodyLayout::piriform_loss),
      layout_slot("branches.piriform_port", &BodyLayout::piriform_port),
      layout_slot("branches.piriform_scale", &BodyLayout::piriform_scale),
      layout_slot("branches.piriform_smoothing", &BodyLayout::piriform_smoothing),
      layout_slot("branches.sinus_neck_scale", &BodyLayout::sinus_neck_scale),
      layout_slot("branches.nose_smoothing", &BodyLayout::nose_smoothing),
      layout_slot("tract.hop_loss", &BodyLayout::vocal_hop_loss),
      layout_slot("branches.buccal_scale", &BodyLayout::buccal_scale),
      layout_slot("branches.velum_leak", &BodyLayout::velum_leak),
      layout_slot("tract.viscous_scale", &BodyLayout::viscous_scale),
      pinch_slot("pinch.on", &PinchSettings::on),
      pinch_slot("pinch.width", &PinchSettings::width),
      pinch_slot("pinch.section_length", &PinchSettings::section_length),
      pinch_slot("pinch.surface_density", &PinchSettings::surface_density),
      pinch_slot("pinch.support", &PinchSettings::support),
      pinch_slot("pinch.elasticity", &PinchSettings::elasticity),
      pinch_slot("pinch.loss", &PinchSettings::loss),
      pinch_slot("pinch.contact_stiffness", &PinchSettings::contact_stiffness),
      pinch_slot("pinch.contact_loss", &PinchSettings::contact_loss),
      pinch_slot("pinch.natural_pitch", &PinchSettings::natural_pitch),
      pinch_slot("pinch.outer_gap", &PinchSettings::outer_gap),
      pinch_slot("pinch.inner_gap", &PinchSettings::inner_gap),
      pinch_slot("pinch.wander", &PinchSettings::wander),
      pinch_slot("pinch.wander_corner", &PinchSettings::wander_corner),
      pinch_slot("pinch.wander_tamping", &PinchSettings::wander_tamping),
      pinch_slot("pinch.jet_turbulence", &PinchSettings::jet_turbulence),
      pinch_slot("pinch.turbulence_corner", &PinchSettings::turbulence_corner),
      pinch_slot("pinch.string_modes", &PinchSettings::string_modes),
      pinch_slot("pinch.string_tension", &PinchSettings::string_tension),
      pinch_slot("pinch.string_loss_growth", &PinchSettings::string_loss_growth),
      pinch_slot("pinch.string_taper", &PinchSettings::string_taper),
      pinch_slot("pinch.s1_mass", &PinchSettings::s1_mass),
      pinch_slot("pinch.s1_elasticity", &PinchSettings::s1_elasticity),
      pinch_slot("pinch.s1_support", &PinchSettings::s1_support),
      pinch_slot("pinch.s2_mass", &PinchSettings::s2_mass),
      pinch_slot("pinch.s2_elasticity", &PinchSettings::s2_elasticity),
      pinch_slot("pinch.s2_support", &PinchSettings::s2_support),
      pinch_slot("pinch.s3_mass", &PinchSettings::s3_mass),
      pinch_slot("pinch.s3_elasticity", &PinchSettings::s3_elasticity),
      pinch_slot("pinch.s3_support", &PinchSettings::s3_support),
      pinch_slot("pinch.s4_mass", &PinchSettings::s4_mass),
      pinch_slot("pinch.s4_elasticity", &PinchSettings::s4_elasticity),
      pinch_slot("pinch.s4_support", &PinchSettings::s4_support),
      pinch_slot("pinch.s5_mass", &PinchSettings::s5_mass),
      pinch_slot("pinch.s5_elasticity", &PinchSettings::s5_elasticity),
      pinch_slot("pinch.s5_support", &PinchSettings::s5_support),
      pinch_slot("pinch.s6_mass", &PinchSettings::s6_mass),
      pinch_slot("pinch.s6_elasticity", &PinchSettings::s6_elasticity),
      pinch_slot("pinch.s6_support", &PinchSettings::s6_support),
      layout_slot("subglottal.bronchial_length", &BodyLayout::bronchial_length),
      layout_slot("subglottal.bronchial_area", &BodyLayout::bronchial_area),
      layout_slot("subglottal.lung_port_scale", &BodyLayout::lung_port_scale),
      tube_slot("tract.hop_loss", SlotKind::vocal_tube),
      tube_slot("tract.nose_hop_loss", SlotKind::nose_tube),
      layout_slot("tract.hop_loss_by_area", &BodyLayout::hop_loss_by_area),
      layout_slot("tract.radiation_mass_scale", &BodyLayout::radiation_mass_scale),
      layout_slot("tract.radiation_resistance_scale", &BodyLayout::radiation_resistance_scale),
      layout_slot("tongue.physical", &BodyLayout::physical_tongue),
      layout_slot("anatomy.scale", &BodyLayout::anatomy_scale),
      layout_slot("anatomy.epilarynx_height", &BodyLayout::epilarynx_height),
      layout_slot("anatomy.epilarynx_area", &BodyLayout::epilarynx_area),
      layout_slot("anatomy.lateral_area", &BodyLayout::lateral_area),
      layout_slot("anatomy.lateral_reach", &BodyLayout::lateral_reach),
      layout_slot("branches.lateral_pocket", &BodyLayout::lateral_pocket),
      layout_slot("branches.lateral_pocket_length", &BodyLayout::lateral_pocket_length),
      layout_slot("branches.lateral_pocket_area", &BodyLayout::lateral_pocket_area),
      layout_slot("branches.lateral_pocket_back", &BodyLayout::lateral_pocket_back),
      layout_slot("tongue.vtl", &BodyLayout::vtl_tongue),
      layout_slot("tongue.vtl_body_radius", &BodyLayout::vtl_body_radius),
      layout_slot("tongue.vtl_tip_radius", &BodyLayout::vtl_tip_radius),
      layout_slot("tongue.vtl_blade_along", &BodyLayout::vtl_blade_along),
      layout_slot("tongue.vtl_blade_lift", &BodyLayout::vtl_blade_lift),
      layout_slot("tongue.vtl_root_bulge", &BodyLayout::vtl_root_bulge),
      tongue_slot("tongue.section_mass", &TongueSettings::section_mass),
      tongue_slot("tongue.support", &TongueSettings::support),
      tongue_slot("tongue.support_loss", &TongueSettings::support_loss),
      tongue_slot("tongue.elasticity", &TongueSettings::elasticity),
      tongue_slot("tongue.bending", &TongueSettings::bending),
      tongue_slot("tongue.elastic_loss", &TongueSettings::elastic_loss),
      tongue_slot("tongue.muscle_force", &TongueSettings::muscle_force),
      tongue_slot("tongue.activation_time", &TongueSettings::activation_time),
      tongue_slot("tongue.gain", &TongueSettings::gain),
      tongue_slot("tongue.damping_gain", &TongueSettings::damping_gain),
      tongue_slot("tongue.integral_rate", &TongueSettings::integral_rate),
      tongue_slot("tongue.contact_depth", &TongueSettings::contact_depth),
      tongue_slot("tongue.tip_pliancy", &TongueSettings::tip_pliancy),
      tongue_slot("tongue.lead", &TongueSettings::lead),
      tongue_slot("tongue.tip_mass", &TongueSettings::tip_mass),
      tongue_slot("tongue.tip_stretch", &TongueSettings::tip_stretch),
      tongue_slot("tongue.statics", &TongueSettings::statics),
      tongue_slot("tongue.contact_stiffness", &TongueSettings::contact_stiffness),
      tongue_slot("tongue.contact_loss", &TongueSettings::contact_loss),
      tongue_slot("tongue.smooth_airway", &TongueSettings::smooth_airway),
      tongue_slot("tongue.vtl_smooth_airway", &TongueSettings::vtl_smooth_airway),
      tongue_slot("tongue.root_physical", &TongueSettings::root_physical),
      tongue_slot("tongue.dorsum_physical", &TongueSettings::dorsum_physical),
      tongue_slot("tongue.blade_physical", &TongueSettings::blade_physical),
      tongue_slot("tongue.tip_physical", &TongueSettings::tip_physical),
      tongue_slot("tongue.roles", &TongueSettings::roles),
      tongue_slot("tongue.shaper_pinned", &TongueSettings::shaper_pinned),
      tongue_slot("tongue.valve_pinned", &TongueSettings::valve_pinned),
      tongue_slot("tongue.nozzle_pinned", &TongueSettings::nozzle_pinned),
      tongue_slot("tongue.shaper_stiffness", &TongueSettings::shaper_stiffness),
      tongue_slot("tongue.valve_stiffness", &TongueSettings::valve_stiffness),
      tongue_slot("tongue.nozzle_stiffness", &TongueSettings::nozzle_stiffness),
      tongue_slot("tongue.shaper_damping", &TongueSettings::shaper_damping),
      tongue_slot("tongue.valve_damping", &TongueSettings::valve_damping),
      tongue_slot("tongue.nozzle_damping", &TongueSettings::nozzle_damping),
      tongue_slot("tongue.role_anchor", &TongueSettings::role_anchor),
      tongue_slot("tongue.lateral_spread", &TongueSettings::lateral_spread),
      tongue_slot("tongue.lateral_narrow", &TongueSettings::lateral_narrow),
      tongue_slot("tongue.lateral_time", &TongueSettings::lateral_time),
      tongue_slot("tongue.lateral_from_posture", &TongueSettings::lateral_from_posture),
      tongue_slot("tongue.feedback_aim", &TongueSettings::feedback_aim),
      tongue_slot("tongue.solid_hull", &TongueSettings::solid_hull),
      tongue_slot("tongue.blade_body", &TongueSettings::blade_body),
      tongue_slot("tongue.skeleton", &TongueSettings::skeleton),
      tongue_slot("tongue.skeleton_base_hz", &TongueSettings::skeleton_base_hz),
      tongue_slot("tongue.skeleton_tip_hz", &TongueSettings::skeleton_tip_hz),
      tongue_slot("tongue.skeleton_rise", &TongueSettings::skeleton_rise),
      tongue_slot("tongue.skeleton_stiffness", &TongueSettings::skeleton_stiffness),
      tongue_slot("tongue.skeleton_tip_weight", &TongueSettings::skeleton_tip_weight),
      tongue_slot("tongue.skeleton_iterations", &TongueSettings::skeleton_iterations),
      tongue_slot("tongue.skeleton_curve", &TongueSettings::skeleton_curve),
      tongue_slot("tongue.skeleton_bend_limit", &TongueSettings::skeleton_bend_limit),
      tongue_slot("tongue.skeleton_anchor_move", &TongueSettings::skeleton_anchor_move),
      tongue_slot("tongue.skeleton_anchor_stiffness", &TongueSettings::skeleton_anchor_stiffness),
      tongue_slot("tongue.skeleton_smoothness", &TongueSettings::skeleton_smoothness),
      tongue_slot("tongue.skeleton_minimal_contact", &TongueSettings::skeleton_minimal_contact),
      tongue_slot("tongue.skeleton_stretch", &TongueSettings::skeleton_stretch),
      tongue_slot("tongue.skeleton_restart", &TongueSettings::skeleton_restart),
      tongue_slot("tongue.skeleton_restart_iterations", &TongueSettings::skeleton_restart_iterations),
      tongue_slot("tongue.spine", &TongueSettings::spine),
      tongue_slot("tongue.spine_base_end", &TongueSettings::spine_base_end),
      tongue_slot("tongue.spine_tip_start", &TongueSettings::spine_tip_start),
      tongue_slot("tongue.spine_base_hz", &TongueSettings::spine_base_hz),
      tongue_slot("tongue.spine_tip_hz", &TongueSettings::spine_tip_hz),
      tongue_slot("tongue.spine_damping", &TongueSettings::spine_damping),
      tongue_slot("tongue.spine_butt", &TongueSettings::spine_butt),
      tongue_slot("tongue.spine_rise", &TongueSettings::spine_rise),
      tongue_slot("tongue.body_thickness", &TongueSettings::body_thickness),
      tongue_slot("tongue.spine_stretch", &TongueSettings::spine_stretch),
      tongue_slot("tongue.spine_stretch_time", &TongueSettings::spine_stretch_time),
      tongue_slot("tongue.spine_reach", &TongueSettings::spine_reach),
      tongue_slot("tongue.spine_reach_damping", &TongueSettings::spine_reach_damping),
      tongue_slot("tongue.spine_reach_steps", &TongueSettings::spine_reach_steps),
      tongue_slot("tongue.tip_radius", &TongueSettings::tip_radius),
      tongue_slot("tongue.blade_thickness", &TongueSettings::blade_thickness),
      tongue_slot("tongue.frenulum_x", &TongueSettings::frenulum_x),
      tongue_slot("tongue.frenulum_y", &TongueSettings::frenulum_y),
      tongue_slot("tongue.guard", &TongueSettings::guard),
      tongue_slot("tongue.guard_keep", &TongueSettings::guard_keep),
      tongue_slot("tongue.guard_stiffness", &TongueSettings::guard_stiffness),
      tongue_slot("tongue.guard_loss", &TongueSettings::guard_loss),
      tongue_slot("tongue.role_lead", &TongueSettings::role_lead),
      tongue_slot("tongue.role_hold", &TongueSettings::role_hold),
      tongue_slot("tongue.role_ramp", &TongueSettings::role_ramp),
      tongue_slot("tongue.pin_stiffness", &TongueSettings::pin_stiffness),
      speaker_slot("mouth.release_speed", &Speaker::set_release_speed, &Speaker::release_speed),
      speaker_slot("mouth.sliver_area", &Speaker::set_sliver_area, &Speaker::sliver_area),
      speaker_slot("mouth.continuous_opening", &Speaker::set_continuous_opening, &Speaker::continuous_opening),
      speaker_slot("mouth.continuous_cuts", &Speaker::set_continuous_cuts, &Speaker::continuous_cuts),
      speaker_slot("mouth.flare_growth", &Speaker::set_flare_growth, &Speaker::flare_growth),
      speaker_slot("mouth.flare_widest", &Speaker::set_flare_widest, &Speaker::flare_widest),
      layout_slot("kinematic.amplitude", &BodyLayout::kinematic_amplitude),
      layout_slot("kinematic.phase_lag", &BodyLayout::kinematic_phase_lag),
      layout_slot("kinematic.softness", &BodyLayout::kinematic_softness),
      layout_slot("kinematic.jitter", &BodyLayout::kinematic_jitter),
      layout_slot("kinematic.skew", &BodyLayout::kinematic_skew),
      layout_slot("kinematic.double_pulsing", &BodyLayout::kinematic_double_pulsing),
      layout_slot("kinematic.flutter", &BodyLayout::kinematic_flutter),
      layout_slot("kinematic.length_scaling", &BodyLayout::kinematic_length_scaling),
      layout_slot("kinematic.geometric_2025", &BodyLayout::geometric_2025),
      layout_slot("kinematic.pulse_shape", &BodyLayout::pulse_shape),
      layout_slot("kinematic.relative_amplitude", &BodyLayout::relative_amplitude),
      layout_slot("kinematic.contact_softness_2025", &BodyLayout::contact_softness_2025),
      layout_slot("kinematic.cycle_shape", &BodyLayout::cycle_shape),
      layout_slot("kinematic.cycle_lag", &BodyLayout::cycle_lag),
      layout_slot("kinematic.cycle_amplitude", &BodyLayout::cycle_amplitude),
      layout_slot("kinematic.cycle_rest", &BodyLayout::cycle_rest),
      layout_slot("kinematic.abduction_fade_start", &BodyLayout::abduction_fade_start),
      layout_slot("kinematic.abduction_fade_end", &BodyLayout::abduction_fade_end),
      layout_slot("kinematic.fold_length_2025", &BodyLayout::fold_length_2025),
      layout_slot("kinematic.fold_thickness_2025", &BodyLayout::fold_thickness_2025),
      layout_slot("kinematic.fold_rest_pitch_2025", &BodyLayout::fold_rest_pitch_2025),
      layout_slot("vtl.glottal_flow", &BodyLayout::vtl_glottal_flow),
      layout_slot("vtl.inertance_floor", &BodyLayout::vtl_inertance_floor),
      layout_slot("vtl.aspiration", &BodyLayout::vtl_aspiration),
      layout_slot("vtl.aspiration_gain", &BodyLayout::vtl_aspiration_gain),
      layout_slot("vtl.aspiration_corner", &BodyLayout::vtl_aspiration_corner),
      layout_slot("vtl.area_limit", &BodyLayout::vtl_area_limit),
  };
  return slots;
}

const ParameterSlot* find_slot(const std::string& path) {
  for (const ParameterSlot& slot : declared_slots()) {
    if (path == slot.path) {
      return &slot;
    }
  }
  return nullptr;
}

}  // namespace

bool apply_layout_parameters(const std::vector<ParameterEntry>& entries, BodyLayout& layout, std::string& error) {
  for (const ParameterEntry& entry : entries) {
    const ParameterSlot* slot = find_slot(entry.path);
    if (slot == nullptr) {
      error = "no such parameter: " + entry.path;
      return false;
    }
    if (slot->kind == SlotKind::layout) {
      layout.*(slot->layout_field) = entry.value;
    } else if (slot->kind == SlotKind::fold) {
      layout.fold.*(slot->fold_field) = entry.value;
    } else if (slot->kind == SlotKind::pinch) {
      layout.pinch.*(slot->pinch_field) = entry.value;
    } else if (slot->kind == SlotKind::tongue) {
      layout.tongue.*(slot->tongue_field) = entry.value;
    }
  }
  return true;
}

bool apply_speaker_parameters(const std::vector<ParameterEntry>& entries, Speaker& speaker, std::string& error) {
  for (const ParameterEntry& entry : entries) {
    const ParameterSlot* slot = find_slot(entry.path);
    if (slot == nullptr) {
      error = "no such parameter: " + entry.path;
      return false;
    }
    if (slot->kind == SlotKind::speaker) {
      (speaker.*(slot->speaker_set))(entry.value);
    } else if (slot->kind == SlotKind::vocal_tube) {
      speaker.body().vocal().set_wave_smoothing(entry.value);
    } else if (slot->kind == SlotKind::nose_tube) {
      speaker.body().nose().set_wave_smoothing(entry.value);
    }
  }
  return true;
}

void collect_parameters(const BodyLayout& layout, Speaker& speaker, std::vector<ParameterEntry>& entries) {
  entries.clear();
  for (const ParameterSlot& slot : declared_slots()) {
    ParameterEntry entry;
    entry.path = slot.path;
    if (slot.kind == SlotKind::layout) {
      entry.value = layout.*(slot.layout_field);
    } else if (slot.kind == SlotKind::fold) {
      entry.value = layout.fold.*(slot.fold_field);
    } else if (slot.kind == SlotKind::pinch) {
      entry.value = layout.pinch.*(slot.pinch_field);
    } else if (slot.kind == SlotKind::tongue) {
      entry.value = layout.tongue.*(slot.tongue_field);
    } else if (slot.kind == SlotKind::speaker) {
      entry.value = (speaker.*(slot.speaker_get))();
    } else if (slot.kind == SlotKind::vocal_tube) {
      entry.value = speaker.body().vocal().wave_smoothing();
    } else {
      entry.value = speaker.body().nose().wave_smoothing();
    }
    entries.push_back(entry);
  }
}

}  // namespace tract
