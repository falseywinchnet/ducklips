#include "tract/speaker.hpp"

#include "voice/geometry.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {


int index_of(voice::Part part, int coordinate) {
  return 2 * static_cast<int>(part) + coordinate;
}

double signed_root(double area) {
  return area >= 0.0 ? std::sqrt(area) : -std::sqrt(-area);
}

double area_of_root(double root) {
  return root <= 0.0 ? 0.0 : root * root;
}

}  // namespace

BodyLayout speech_layout() {
  std::vector<double> neutral(static_cast<std::size_t>(voice::kOralSections), 3.0);
  BodyLayout layout = vowel_layout(neutral, kNodeLength);
  return layout;
}

BodyLayout articulatory_layout() {
  constexpr std::size_t kNodes = 60;
  std::vector<double> neutral(kNodes, 3.0);
  BodyLayout layout = vowel_layout(neutral, kNodeLength);
  // The source: the kinematic folds (VocalTractLab 2.4's default kind; the
  // owner's A/B: the self-oscillating folds "still clicks a little").
  layout.kinematic_glottis = true;
  layout.pitch_feedback = 0.0;
  // The piriform recesses at a quarter of the speech body's cross-sections
  // (same depth: the 4.5 kHz zero stays). At full size their compliance
  // beside the narrow hypopharynx lowered /aa/'s F1 from 770 to 630 Hz
  // (an ideal transmission line agrees: 590); on Story's MRI shapes, which
  // already give real formants without them, that cannot be a real
  // speaker's, and the vowel tuning flattened every shape toward a uniform
  // tube to win F1 back (the owner's "talking through a tube").
  // Now full size again, the vowels retuned under the shape constraint: at
  // a quarter, voiced sound kept 7 dB too much above 4 kHz (harmonic and
  // noise alike) against the owner's NATO reading; at full size within 2
  // dB in every band (tract/analysis/nato_snr.py, /tmp/es/hn.py).
  layout.piriform_scale = 1.0;
  // The voice's breath: turbulent flux (rms 4% of the flow, gated to
  // voicing), a phase random walk of 0.5 rad/sqrt(s) and a firmer contact
  // (5 um). Without it the band HNR above 1 kHz stood 6-10 dB over bdl's
  // and the owner's ("Perfect Paul"); with the flux low-passed at 1.5 kHz,
  // 4-8 kHz fell 6-10 dB under the owner's recordings (m51b, m52a).
  layout.jet_flux = 0.04;
  layout.kinematic_jitter = 0.5;
  layout.kinematic_softness = 0.0005;
  // Frication 2.0 (the speech body's 0.8): /f/ 10.6 -> 3.6 dB from bdl's
  // bands, /s/ 8.5 -> 6.4, /sh/ 13.4 -> 6.6, with the flux gated (ARCTIC).
  layout.frication = 2.0;
  // Lip tissue where lips usually are (vowel_layout put it in the last cm).
  for (std::size_t node = 0; node < layout.vocal_materials.size(); ++node) {
    double x = (static_cast<double>(node) + 0.5) * kNodeLength;
    layout.vocal_materials[node] = x < 2.0 ? 0 : x > 16.6 ? 2 : 1;
  }
  return layout;
}

bool is_organ_part(voice::Part part) {
  return part == voice::Part::body || part == voice::Part::apex || part == voice::Part::mouth ||
         part == voice::Part::spread || part == voice::Part::larynx || part == voice::Part::angle ||
         part == voice::Part::groove;
}

bool uses_organs(const voice::Score& score) {
  for (const voice::PartEvent& event : score.parts) {
    if (is_organ_part(event.part)) {
      return true;
    }
  }
  return false;
}

Speaker::Speaker(const BodyLayout& layout, const voice::Score& score, std::uint32_t seed)
    : score_(score), body_(layout), random_(seed != 0U ? seed : 1U) {
  // The fold edges' damping. VocalTractLab's 0.1, 0.6 left the loudest
  // vowels 15-25 dB too strong at 1.5-4 kHz against bdl, the owner and
  // VocalTractLab's own synthesis (a flat, buzzy envelope; vowel tilt -13
  // dB). 0.3, 0.8 darkened it (-23) but rounded the pulse: H1-H2 +4.8 dB
  // (the owner -3.5, bdl -2.7), a hooty voice. 0.2, 0.8 with the upper edge
  // pressed to -0.012: tilt -20, H1-H2 +0.6; identification 38.9 -> 39.6%.
  body_.cord().set_damping(0.2, 0.8);
  // A loss growing with frequency on every hop of the vocal tube (k 0.05):
  // the upper resonances (5-7.5 kHz) rang as an even comb, a uniform
  // tube's, where a person's are damped away; costs about 5% on every
  // formant's frequency (an effective tract 5% longer).
  body_.vocal().set_wave_smoothing(0.05);
  // The nose's narrow, mucosa-lined passages take the highs much harder:
  // nasal murmurs sat 9-11 dB over bdl at 3-7 kHz; at 0.2 the top band is
  // within 4 dB of bdl. 0.4 with the narrowed sinus ostia (see
  // BodyLayout::sinus_neck_scale): at 0.2, 2-4 kHz stood 4-6 dB over bdl.
  body_.nose().set_wave_smoothing(layout.nose_smoothing);
  // Rest postures, wander per unit effort and wander time (the front end's
  // body: these are the speeds and unsteadiness it was tuned with).
  set_motor(voice::Part::lungs, 0, 0.0, 12.0, 0.35, 1000.0);
  set_motor(voice::Part::pitch, 0, 120.0, 0.8, 0.25, 100.0);
  set_motor(voice::Part::glottis, 0, 0.12, 0.0015, 0.12, 0.1);
  set_motor(voice::Part::glottis, 1, 0.5, 0.03, 0.2, 1.0);
  set_motor(voice::Part::tongue, 0, 9.5, 0.06, 0.18, 4.0);
  set_motor(voice::Part::tongue, 1, 1.6, 0.015, 0.15, 1.0);
  set_motor(voice::Part::tip, 0, 15.8, 0.04, 0.2, 1.0);
  set_motor(voice::Part::tip, 1, std::sqrt(static_cast<double>(voice::kReleased)), 0.012, 0.12, 1.0);
  set_motor(voice::Part::dorsum, 0, 12.0, 0.04, 0.2, 1.0);
  set_motor(voice::Part::dorsum, 1, std::sqrt(static_cast<double>(voice::kReleased)), 0.012, 0.12, 1.0);
  set_motor(voice::Part::lips, 0, 1.6, 0.015, 0.12, 1.0);
  set_motor(voice::Part::lips, 1, 0.0, 0.02, 0.2, 1.0);
  set_motor(voice::Part::velum, 0, 0.0, 0.01, 0.2, 1.0);
  set_motor(voice::Part::floor, 0, 0.5, 0.05, 0.2, 2.0);
  set_motor(voice::Part::jaw, 0, voice::kReleased, 0.01, 0.2, 8.0);
  set_motor(voice::Part::jaw, 1, 0.0, 0.0, 0.2, 1.0);
  set_motor(voice::Part::lateral, 0, 0.0, 0.0, 0.2, 0.5);
  // The articulatory body's organs (cm, degrees; tongue points in the jaw's frame).
  set_motor(voice::Part::body, 0, 4.0, 0.03, 0.18, 1.5);
  set_motor(voice::Part::body, 1, 6.6, 0.03, 0.18, 1.5);
  set_motor(voice::Part::apex, 0, 6.9, 0.02, 0.15, 1.0);
  set_motor(voice::Part::apex, 1, 7.4, 0.02, 0.15, 1.0);
  set_motor(voice::Part::mouth, 0, 1.0, 0.01, 0.12, 1.0);
  set_motor(voice::Part::mouth, 1, 0.0, 0.01, 0.2, 1.0);
  set_motor(voice::Part::spread, 0, 3.2, 0.03, 0.2, 1.5);
  set_motor(voice::Part::larynx, 0, 0.0, 0.01, 0.3, 1.0);
  set_motor(voice::Part::angle, 0, 6.0, 0.15, 0.2, 10.0);
  set_motor(voice::Part::groove, 0, 0.0, 0.0, 0.2, 1.0);
  articulatory_ = uses_organs(score_) && body_.vocal().count() >= 56;
  physical_tongue_ = layout.physical_tongue > 0.5;
  tongue_settings_ = layout.tongue;
  {
    VtlTongueShape shape;
    shape.on = layout.vtl_tongue > 0.5 && !physical_tongue_;
    shape.body_rx = layout.vtl_body_radius;
    shape.body_ry = layout.vtl_body_radius;
    shape.tip_radius = layout.vtl_tip_radius;
    shape.blade_along = layout.vtl_blade_along;
    shape.blade_lift = layout.vtl_blade_lift;
    shape.root_bulge = layout.vtl_root_bulge;
    anatomy_.set_vtl_tongue(shape);
  }
  lateral_pocket_ = layout.lateral_pocket;
  lateral_pocket_back_ = layout.lateral_pocket_back;
  abduction_fade_start_ = layout.abduction_fade_start;
  abduction_fade_end_ = layout.abduction_fade_end;
  if (layout.anatomy_scale != 1.0) {
    anatomy_.set_scale(layout.anatomy_scale);
  }
  if (layout.epilarynx_height >= 0.0) {
    anatomy_.set_larynx_height(layout.epilarynx_height);
  }
  if (layout.epilarynx_area != 1.0) {
    anatomy_.set_larynx_area_scale(layout.epilarynx_area);
  }
  if (layout.lateral_area != 0.12) {
    anatomy_.set_lateral_area(layout.lateral_area);
  }
  if (layout.lateral_reach != 2.0) {
    anatomy_.set_lateral_reach(layout.lateral_reach);
  }
  if (articulatory_) {
    // The articulatory body's voice (2026-10-15, against bdl's vowel bands
    // over a0001-a0060): folds resting at 0.011 / 0.006 cm and a lighter
    // hop loss. At 0.015 / 0.010 with loss 0.05 the vowels were 6-8 dB
    // short at 3-8 kHz (the owner: "muffled"); now within about 3 dB.
    voicing_lower_ = 0.011;
    voicing_upper_ = 0.006;
    // (The blown fold tube rests here too. Light, at these rests it closed
    // within 60 us (1.9e6 cm^3/s^2) and clicked every cycle; it was given
    // 0.030 / 0.015. With three times the tissue and a soft contact for the
    // air, it closes at 4.3e5 here (a held /aa/, open quotient 0.65); at
    // 0.030 / 0.015 it closed at 2.7e5, open 75%: the highs fell under the
    // noise, band HNR 4-7 kHz -13 dB, and 0-120 Hz rose to -12 dB.)
    body_.vocal().set_wave_smoothing(0.02);
    if (layout.vocal_hop_loss >= 0.0) {
      body_.vocal().set_wave_smoothing(layout.vocal_hop_loss);
    }
    // The lips radiate from wherever the opening is (set_opening), not from
    // the tube's last node.
    body_.vocal().set_radiating_right(false);
    std::size_t nodes = static_cast<std::size_t>(body_.vocal().count());
    node_from_.assign(nodes, 3.0);
    node_to_.assign(nodes, 3.0);
    node_scratch_.assign(nodes, 3.0);
  }
  // The score's opening postures are where the body starts.
  while (next_event_ < score_.parts.size() && score_.parts[next_event_].time <= 0.0) {
    apply(score_.parts[next_event_]);
    ++next_event_;
  }
  for (Motor& each : motors_) {
    each.position = each.target;
    each.velocity = 0.0;
  }
  shape_body();
  if (articulatory_) {
    // Start on the opening shape, not gliding into it (nor from a default
    // length: that popped at the first sample).
    length_to_ = profile_.length;
    profile_to_nodes(node_to_, length_to_);
    node_from_ = node_to_;
    length_from_ = length_to_;
    for (std::size_t node = 0; node < node_to_.size(); ++node) {
      body_.vocal().set_rest_area(static_cast<Index>(node), node_to_[node]);
    }
  }
  // The body starts in the opening posture, not gliding into it from the
  // neutral layout: that glide moved the whole tract within milliseconds and
  // the air it pushed was the loudest event of every utterance (a thump
  // 0.7 of full scale, which the normalisation then took as the peak).
  body_.vocal().settle_geometry();
  // The tongue spans the tract from the top of the epilarynx to the incisors.
  body_.vocal().set_conserved_span(static_cast<Index>(std::lround(3.5 / kNodeLength)),
                                   static_cast<Index>(std::lround(16.0 / kNodeLength)));
}

Speaker::Motor& Speaker::motor(voice::Part part, int coordinate) {
  return motors_[static_cast<std::size_t>(index_of(part, coordinate))];
}

double Speaker::value(voice::Part part, int coordinate) {
  return motor(part, coordinate).position;
}

void Speaker::set_motor(voice::Part part, int coordinate, double rest, double wander_scale, double wander_seconds,
                        double span) {
  Motor& each = motor(part, coordinate);
  each.rest = rest;
  each.position = rest;
  each.target = rest;
  each.aim = rest;
  each.wander_scale = wander_scale;
  each.wander_seconds = wander_seconds;
  each.span = span;
}

double Speaker::gauss() {
  double sum = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    random_ ^= random_ << 13U;
    random_ ^= random_ >> 17U;
    random_ ^= random_ << 5U;
    sum += static_cast<double>(random_) / 4294967296.0;
  }
  return (sum - 2.0) * std::sqrt(3.0);
}

// A part's new target. Constriction areas are carried as signed roots, so
// closing and pressing move smoothly through zero.
void Speaker::apply(const voice::PartEvent& event) {
  double move = std::max(static_cast<double>(event.move), 0.004);
  if (move_scale_ != 1.0 && event.part != voice::Part::lungs && event.part != voice::Part::pitch &&
      event.part != voice::Part::glottis) {
    move = std::max(move * move_scale_, 0.004);
  }
  double omega = 4.5 / move;
  Motor& first = motor(event.part, 0);
  Motor& second = motor(event.part, 1);
  double a = static_cast<double>(event.a);
  double b = static_cast<double>(event.b);
  switch (event.part) {
    case voice::Part::lungs:
    case voice::Part::pitch:
    case voice::Part::floor:
    case voice::Part::lateral:
      first.target = a;
      first.omega = omega;
      break;
    case voice::Part::jaw:
    case voice::Part::glottis:
      first.target = a;
      second.target = b;
      first.omega = omega;
      second.omega = omega;
      break;
    case voice::Part::tongue:
    case voice::Part::tip:
    case voice::Part::dorsum:
      first.target = a;
      second.target = signed_root(b);
      first.omega = omega;
      second.omega = omega;
      if (event.part != voice::Part::tongue && second.position < 0.3 && second.target > second.position) {
        second.omega = omega * release_speed_;  // a closure opening
      }
      break;
    case voice::Part::lips:
      first.target = signed_root(a);
      second.target = b;
      first.omega = omega;
      second.omega = omega;
      if (first.position < 0.3 && first.target > first.position) {
        first.omega = omega * release_speed_;
      }
      break;
    case voice::Part::velum:
      first.target = std::sqrt(std::max(a, 0.0));
      first.omega = omega;
      break;
    case voice::Part::body:
    case voice::Part::apex:
      first.target = a;
      second.target = b;
      first.omega = omega;
      second.omega = omega;
      break;
    case voice::Part::mouth:
      first.target = a;
      second.target = b;
      first.omega = omega;
      second.omega = omega;
      if (first.position < 0.05 && first.target > first.position) {
        first.omega = omega * release_speed_;  // the lips parting
      }
      break;
    case voice::Part::spread:
    case voice::Part::larynx:
    case voice::Part::angle:
    case voice::Part::groove:
      first.target = a;
      first.omega = omega;
      break;
    case voice::Part::count:
      break;
  }
}

// Muscles: each pulls its part toward the target plus its wander; the
// part's own mass smooths the wander. Holding a shape far from rest costs
// effort, and effort makes the hold less steady.
void Speaker::advance_motors() {
  for (Motor& each : motors_) {
    double goal = each.target;
    if (motor_smoothing_ > 0.0) {
      each.aim += (1.0 - std::exp(-kControlSeconds * motor_smoothing_ * each.omega)) * (each.target - each.aim);
      goal = each.aim;
    }
    double pull = each.omega * each.omega * (goal + each.wander - each.position);
    double drag = 2.0 * each.omega * each.velocity;
    each.velocity += kControlSeconds * (pull - drag);
    each.position += kControlSeconds * each.velocity;
    double effort = 0.4 + std::abs(each.target - each.rest) / each.span;
    double decay = std::exp(-kControlSeconds / each.wander_seconds);
    each.wander = each.wander * decay + wander_gain_ * each.wander_scale * effort * std::sqrt(1.0 - decay * decay) * gauss();
  }
}

double Speaker::shape_area_parts() {
  voice::TractShape shape;
  shape.tongue_place = value(voice::Part::tongue, 0);
  shape.tongue_area = std::max(area_of_root(value(voice::Part::tongue, 1)), 0.05);
  shape.tip_place = value(voice::Part::tip, 0);
  shape.tip_root = value(voice::Part::tip, 1);
  shape.dorsum_place = value(voice::Part::dorsum, 0);
  shape.dorsum_root = value(voice::Part::dorsum, 1);
  shape.lip_root = value(voice::Part::lips, 0);
  shape.protrusion = value(voice::Part::lips, 1);
  shape.teeth_gap = value(voice::Part::jaw, 0);
  shape.jaw_closure = value(voice::Part::jaw, 1);
  std::array<double, voice::kOralSections> areas{};
  voice::TractLayout layout = voice::tract_areas(shape, areas);
  Tube& vocal = body_.vocal();
  Index count = std::min(vocal.count(), static_cast<Index>(voice::kOralSections));
  double sealed = areas[static_cast<std::size_t>(voice::kOralSections - 1)];
  for (Index node = static_cast<Index>(std::lround(10.5 / kNodeLength)); node < count; ++node) {
    sealed = std::min(sealed, areas[static_cast<std::size_t>(node)]);
  }
  for (Index node = 0; node < count; ++node) {
    vocal.set_rest_area(node, areas[static_cast<std::size_t>(node)]);
  }
  vocal.set_active_length(static_cast<double>(count) * kNodeLength);
  body_.set_teeth(layout.teeth_start, layout.teeth_end);
  return sealed;
}

// The articulatory body: organ positions through the midsagittal model,
// once every kGeometryTicks control ticks; between frames the node areas
// and the open end's length glide linearly from the last frame to this one.
double Speaker::shape_organs() {
  Tube& vocal = body_.vocal();
  if (geometry_tick_ == 0) {
    Articulators held = organs();
    if (physical_tongue_) {
      move_tongue(held, static_cast<double>(kGeometryTicks) * kControlSeconds);
    } else if (anatomy_.vtl_tongue().on) {
      anatomy_.vtl_tongue_outline(held, vtl_contour_);
      held.contour = &vtl_contour_;
      held.nearest_width = tongue_settings_.vtl_smooth_airway > 0.5;
    }
    anatomy_.compute(held, profile_);
    node_from_.swap(node_to_);
    // (The glide starts from where the last frame's glide ended.) The
    // length follows the geometry through an 8 ms one-pole: a contact
    // regrouping the cuts shifted it by a node or more within a frame, the
    // lips' nodes swapped with the room's and pumped 3000 cm^3/s, a "tap on
    // the microphone".
    length_from_ = length_to_;
    length_to_ += (1.0 - std::exp(-1.0e-3 / 0.008)) * (profile_.length - length_to_);
    profile_to_nodes(node_to_, length_to_);
    body_.set_velum_place(profile_.velum_arc);
    body_.set_under_tongue_place(profile_.tip_arc);
    if (lateral_pocket_ > 0.0) {
      body_.set_lateral_pocket_place(profile_.tip_arc - lateral_pocket_back_);
    }
    body_.set_teeth(profile_.teeth_arc - 0.25, profile_.teeth_arc + 0.25);
  }
  double share = static_cast<double>(geometry_tick_ + 1) / static_cast<double>(kGeometryTicks);
  geometry_tick_ = (geometry_tick_ + 1) % kGeometryTicks;
  Index count = vocal.count();
  for (Index node = 0; node < count; ++node) {
    std::size_t index = static_cast<std::size_t>(node);
    vocal.set_rest_area(node, node_from_[index] + share * (node_to_[index] - node_from_[index]));
  }
  double length = length_from_ + share * (length_to_ - length_from_);
  // The whole tube carries air; the open end flares past the lips.
  // The radiation load sits at the lip opening, shared between the two
  // nodes around it, sliding as the length moves.
  double position = std::max(length / kNodeLength - 0.5, 1.0);
  Index lip_end = std::clamp(static_cast<Index>(std::floor(position)), static_cast<Index>(2), count - 3);
  double lip_area = profile_.count > 0 ? profile_.area[static_cast<std::size_t>(profile_.count - 1)] : 1.0;
  vocal.set_opening(std::clamp(position, 2.0, static_cast<double>(count - 2) - 1.0e-6), lip_area);
  body_.set_lip_end(lip_end);
  vocal.set_quiet_from(lip_end + 2);
  // The narrowest point of the mouth (past 10.5 cm, up to the lips).
  Index lips = std::min(static_cast<Index>(std::lround(length / kNodeLength)), count);
  double sealed = 10.0;
  for (Index node = static_cast<Index>(std::lround(10.5 / kNodeLength)); node < lips; ++node) {
    std::size_t index = static_cast<std::size_t>(node);
    sealed = std::min(sealed, node_from_[index] + share * (node_to_[index] - node_from_[index]));
  }
  return sealed;
}

void Speaker::move_tongue(Articulators& held, double seconds) {
  Articulators aimed = held;
  aimed.tongue_x = motor(voice::Part::body, 0).target;
  aimed.tongue_y = motor(voice::Part::body, 1).target;
  aimed.tip_x = motor(voice::Part::apex, 0).target;
  aimed.tip_y = motor(voice::Part::apex, 1).target;
  // The tongue plans ahead: a body with mass must start before its target is
  // due (the planner's leads were set for 28 ms motors), so it aims at the
  // next body and tip targets due within tongue.lead of now. (The next, not
  // the last: aiming at the last, a consonant's contact was skipped for the
  // vowel after it.)
  bool body_seen = false;
  bool apex_seen = false;
  for (std::size_t look = next_event_;
       look < score_.parts.size() && score_.parts[look].time <= time_ + tongue_settings_.lead; ++look) {
    const voice::PartEvent& coming = score_.parts[look];
    if (coming.part == voice::Part::body && !body_seen) {
      aimed.tongue_x = static_cast<double>(coming.a);
      aimed.tongue_y = static_cast<double>(coming.b);
      body_seen = true;
    } else if (coming.part == voice::Part::apex && !apex_seen) {
      aimed.tip_x = static_cast<double>(coming.a);
      aimed.tip_y = static_cast<double>(coming.b);
      apex_seen = true;
    }
  }
  anatomy_.tongue_outline(aimed, tongue_target_);
  if (!tongue_.ready()) {
    // At rest: the motors' rests (the tongue relaxed, jaw as now).
    Articulators rest = held;
    rest.tongue_x = motor(voice::Part::body, 0).rest;
    rest.tongue_y = motor(voice::Part::body, 1).rest;
    rest.tip_x = motor(voice::Part::apex, 0).rest;
    rest.tip_y = motor(voice::Part::apex, 1).rest;
    TongueContour relaxed;
    anatomy_.tongue_outline(rest, relaxed);
    tongue_.reset(relaxed, tongue_settings_);
    tongue_.place(tongue_target_);
  }
  tongue_.set_target(tongue_target_);
  if (tongue_settings_.roles > 0.5) {
    role_state(phones_, phone_cursor_, time_, tongue_settings_, tongue_roles_);
    tongue_.set_impedance(tongue_roles_.stiffness, tongue_roles_.damping);
    tongue_.set_shaping(tongue_roles_.shaping);
  }
  if (tongue_.uses_geometric()) {
    // Where pinned pieces go: the geometric tongue as the motors hold it now.
    anatomy_.tongue_outline(held, tongue_geometric_);
    tongue_.set_geometric(tongue_geometric_, seconds);
    if (tongue_settings_.roles > 0.5) {
      tongue_.set_pins(tongue_roles_.pins);
    }
  }
  tongue_.step(seconds, held, anatomy_);
  tongue_frame_ = held;
  tongue_frame_.contour = nullptr;
  double lifted = std::max(tongue_.lift(), 0.0);
  tongue_lift_ += (1.0 - std::exp(-seconds / 0.02)) * (lifted - tongue_lift_);
  held.contour = &tongue_.contour();
  bool thick = tongue_settings_.blade_body > 0.5 || tongue_settings_.spine > 0.5 || tongue_settings_.skeleton > 0.5;
  held.underside = thick ? &tongue_.underside() : nullptr;
  held.nearest_width = tongue_settings_.smooth_airway > 0.5;
}

Articulators Speaker::organs() {
  Articulators held;
  held.tongue_x = value(voice::Part::body, 0);
  held.tongue_y = value(voice::Part::body, 1);
  held.tip_x = value(voice::Part::apex, 0);
  held.tip_y = value(voice::Part::apex, 1);
  held.lip_aperture = value(voice::Part::mouth, 0);
  held.lip_protrusion = std::max(value(voice::Part::mouth, 1), 0.0);
  held.lip_width = value(voice::Part::spread, 0);
  held.larynx_drop = std::max(value(voice::Part::larynx, 0) - tongue_pull_raise_ * tongue_lift_, 0.0);
  held.jaw_angle = value(voice::Part::angle, 0);
  held.lateral = std::max(value(voice::Part::lateral, 0), 0.0);
  held.groove = std::max(value(voice::Part::groove, 0), 0.0);
  return held;
}

// Node areas from the profile: the area at the node's centre, or less if a
// narrower cut falls within the node's span (a closure is never stepped
// over); past the lip opening, the lips' area.
double Speaker::profile_area_at(double arc) const {
  int cuts = profile_.count;
  if (arc <= profile_.arc[0]) {
    return profile_.area[0];
  }
  for (int cut = 1; cut < cuts; ++cut) {
    std::size_t here = static_cast<std::size_t>(cut);
    if (profile_.arc[here] >= arc) {
      double span = std::max(profile_.arc[here] - profile_.arc[here - 1], 1.0e-6);
      double weight = std::clamp((arc - profile_.arc[here - 1]) / span, 0.0, 1.0);
      return (1.0 - weight) * profile_.area[here - 1] + weight * profile_.area[here];
    }
  }
  return profile_.area[static_cast<std::size_t>(cuts - 1)];
}

void Speaker::profile_to_nodes(std::vector<double>& areas, double length) const {
  int cuts = profile_.count;
  int cut = 0;
  for (std::size_t node = 0; node < areas.size(); ++node) {
    double low = static_cast<double>(node) * kNodeLength;
    double centre = low + 0.5 * kNodeLength;
    double high = low + kNodeLength;
    if (cuts == 0) {
      areas[node] = 3.0;
      continue;
    }
    if (centre >= length) {
      // Past the lips: a near-closed sliver (the radiation itself is a load
      // at the lip opening, Tube::set_opening). (A flaring "room" here acted
      // as part of the tract: /aa/'s F1 widened to 140-170 Hz and the 1-2 kHz
      // band stood 7-8 dB over bdl's, the owner's "merman gurgling".) The
      // node the opening crosses blends by how much of it lies inside.
      double lips = std::max(profile_.area[static_cast<std::size_t>(cuts - 1)], 0.0);
      double flared = sliver_area_;
      double inside = std::clamp((length - low) / kNodeLength, 0.0, 1.0);
      areas[node] = inside * lips + (1.0 - inside) * flared;
      continue;
    }
    while (cut + 1 < cuts && profile_.arc[static_cast<std::size_t>(cut + 1)] < centre) {
      ++cut;
    }
    std::size_t here = static_cast<std::size_t>(cut);
    double area = profile_.area[here];
    double inside = continuous_opening_ != 0.0 ? std::clamp((length - low) / kNodeLength, 0.0, 1.0) : 1.0;
    if (cut + 1 < cuts && centre > profile_.arc[here]) {
      std::size_t next = here + 1;
      double span = std::max(profile_.arc[next] - profile_.arc[here], 1.0e-6);
      double weight = std::clamp((centre - profile_.arc[here]) / span, 0.0, 1.0);
      area = (1.0 - weight) * profile_.area[here] + weight * profile_.area[next];
    }
    for (int each = std::max(cut - 2, 0); each < std::min(cut + 4, cuts); ++each) {
      std::size_t index = static_cast<std::size_t>(each);
      if (profile_.arc[index] >= low && profile_.arc[index] < high) {
        area = std::min(area, profile_.area[index]);
      }
    }
    if (continuous_cuts_ != 0.0) {
      area = std::min(area, std::min(profile_area_at(low), profile_area_at(std::min(high, length))));
    }
    areas[node] = inside * area + (1.0 - inside) * sliver_area_;
  }
}

void Speaker::shape_body() {
  Tube& vocal = body_.vocal();
  Index count = vocal.active_count();
  // Voiced closure: the mouth sealed above the pharynx while the folds are
  // set to vibrate. English speakers keep voicing through /b d g/ by letting
  // the cavity behind the closure grow; the pharyngeal walls relax and yield
  // to the air the folds pass, so the pressure across the folds holds for a
  // while. Their support stiffness is released toward (1 - depth) over about
  // 20 ms and restored over about 30 ms once the mouth opens or the folds
  // part. (Restored with a 40 ms time constant, the pull still held the
  // pharynx 30-40% wide 100 ms into the next vowel: after /b/, /aa/ F2 1420
  // against its 1240, "body" heard as "bady".)
  // (Prescribing the enlargement instead, +35% pharynx area over 60 ms,
  // sucked the cavity to -1.7 kPa, more than the folds could fill,
  // and the suction pulled them shut. bdl voices /b d g/ 61% of their
  // duration; ours was 36%.)
  double sealed = articulatory_ ? shape_organs() : shape_area_parts();
  bool voicing = value(voice::Part::glottis, 0) < 0.03;
  bool closure = sealed < 0.02 && voicing;
  double tau = closure ? 0.05 : 0.01;
  relaxation_ += (1.0 - std::exp(-kControlSeconds / tau)) * ((closure ? 1.0 : 0.0) - relaxation_);
  Index pharynx_first = static_cast<Index>(std::lround(3.5 / kNodeLength));
  Index pharynx_last = std::min(static_cast<Index>(std::lround(10.0 / kNodeLength)), count - 1);
  body_.set_wall_release(pharynx_first, pharynx_last, relaxation_depth_ * relaxation_);
  body_.set_wall_pull(pharynx_first, pharynx_last, pull_depth_ * relaxation_);

  // Arytenoid posture. The score's rest half-gap (cm) is about -0.006 for
  // modal voice, 0 while voicing through a closure, 0.14 for voiceless
  // sounds and 0.45 for a breath; the cord's rest displacements at the
  // vocal processes are 0.005 (modal) to 0.1 (voiceless fricative) in
  // VocalTractLab's terms. Every voicing setting (gap up to 0.002) is the
  // modal posture; wider settings open from there. Compression 0..1 narrows
  // the lower edge a little. (Mapping the closure-voicing gap 0 to a posture
  // twice as open as modal, the folds sat apart through every /b d g/
  // closure, leaking a steady 75 cm^3/s instead of vibrating: voiced stops
  // were voiced 27% of their duration, bdl's 61%.)
  double gap = value(voice::Part::glottis, 0);
  double compression = std::clamp(value(voice::Part::glottis, 1), 0.0, 1.0);
  double open = 0.7 * std::max(gap - 0.002, 0.0);
  CordPosture posture;
  // In a voiced closure the airway's pressure inflates the folds apart
  // (0.04 cm^2, a slow wobble that never collides: the voice bar was 32-40 dB
  // under the vowel, real ones 20-30); the adductors hold them together.
  double press = press_depth_ * relaxation_;
  double pulled = tongue_pull_adduction_ * tongue_lift_;
  posture.lower_rest = std::max(voicing_lower_ + open - 0.002 * (compression - 0.5) - press - pulled, -0.006);
  posture.upper_rest = std::max(voicing_upper_ + open - press - pulled, -0.008);
  posture.chink = voicing_chink_ + 0.12 * open;
  body_.cord().set_target(posture, 0.015);
  if (abduction_fade_end_ > abduction_fade_start_ && abduction_fade_end_ > 0.0) {
    body_.cord().set_abduction_scale((abduction_fade_end_ - gap) / (abduction_fade_end_ - abduction_fade_start_));
  }

  body_.set_lung_pressure(value(voice::Part::lungs, 0), lung_follow_);
  body_.set_pitch(std::max(value(voice::Part::pitch, 0), 40.0) * (1.0 + tongue_pull_tension_ * tongue_lift_));
  double velum_root = std::max(value(voice::Part::velum, 0), 0.0);
  body_.set_velum(velum_root * velum_root);
  // A raised tip (/r/, /sh/, /er/) lifts the blade off the floor of the
  // mouth and opens the space beneath it.
  double floor = value(voice::Part::floor, 0);
  body_.set_under_tongue(0.4 * std::clamp((floor - 0.5) / 1.5, 0.0, 1.0));
  if (lateral_pocket_ > 0.0) {
    body_.set_lateral_pocket(lateral_pocket_ * std::clamp(value(voice::Part::lateral, 0), 0.0, 1.0));
  }
}

void Speaker::push(const voice::PartEvent& event) {
  std::vector<voice::PartEvent>& parts = score_.parts;
  std::size_t place = parts.size();
  while (place > next_event_ && parts[place - 1].time > event.time) {
    --place;
  }
  parts.insert(parts.begin() + static_cast<std::ptrdiff_t>(place), event);
}

double Speaker::step() {
  if (tick_ % kControlSteps == 0) {
    while (next_event_ < score_.parts.size() && score_.parts[next_event_].time <= time_) {
      apply(score_.parts[next_event_]);
      ++next_event_;
    }
    advance_motors();
    shape_body();
  }
  ++tick_;
  time_ = static_cast<double>(tick_) * kStepSeconds;
  return body_.step();
}

}  // namespace tract
