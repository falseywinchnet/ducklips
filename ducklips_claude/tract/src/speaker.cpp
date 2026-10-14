#include "tract/speaker.hpp"

#include "voice/geometry.hpp"

#include <algorithm>
#include <cmath>

namespace tract {

namespace {

constexpr double kControlSeconds = 1.0e-3;
constexpr long kStepsPerControl = 96;

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
  body_.nose().set_wave_smoothing(0.4);
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
    case voice::Part::count:
      break;
  }
}

// Muscles: each pulls its part toward the target plus its wander; the
// part's own mass smooths the wander. Holding a shape far from rest costs
// effort, and effort makes the hold less steady.
void Speaker::advance_motors() {
  for (Motor& each : motors_) {
    double pull = each.omega * each.omega * (each.target + each.wander - each.position);
    double drag = 2.0 * each.omega * each.velocity;
    each.velocity += kControlSeconds * (pull - drag);
    each.position += kControlSeconds * each.velocity;
    double effort = 0.4 + std::abs(each.target - each.rest) / each.span;
    double decay = std::exp(-kControlSeconds / each.wander_seconds);
    each.wander = each.wander * decay + wander_gain_ * each.wander_scale * effort * std::sqrt(1.0 - decay * decay) * gauss();
  }
}

void Speaker::shape_body() {
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
  // Voiced closure: the mouth sealed above the pharynx while the folds are
  // set to vibrate. English speakers keep voicing through /b d g/ by letting
  // the cavity behind the closure grow; the pharyngeal walls relax and yield
  // to the air the folds pass, so the pressure across the folds holds for a
  // while. Their support stiffness is released toward (1 - depth) over about
  // 20 ms and restored over about 30 ms once the mouth opens or the folds
  // part. (Prescribing the enlargement instead, +35% pharynx area over 60 ms,
  // sucked the cavity to -1.7 kPa, more than the folds could fill,
  // and the suction pulled them shut. bdl voices /b d g/ 61% of their
  // duration; ours was 36%.)
  double sealed = areas[static_cast<std::size_t>(voice::kOralSections - 1)];
  for (Index node = static_cast<Index>(std::lround(10.5 / kNodeLength)); node < count; ++node) {
    sealed = std::min(sealed, areas[static_cast<std::size_t>(node)]);
  }
  bool voicing = value(voice::Part::glottis, 0) < 0.03;
  bool closure = sealed < 0.02 && voicing;
  double tau = closure ? 0.05 : 0.04;
  relaxation_ += (1.0 - std::exp(-kControlSeconds / tau)) * ((closure ? 1.0 : 0.0) - relaxation_);
  Index pharynx_first = static_cast<Index>(std::lround(3.5 / kNodeLength));
  Index pharynx_last = std::min(static_cast<Index>(std::lround(10.0 / kNodeLength)), count - 1);
  for (Index node = 0; node < count; ++node) {
    vocal.set_rest_area(node, areas[static_cast<std::size_t>(node)]);
  }
  body_.set_wall_release(pharynx_first, pharynx_last, relaxation_depth_ * relaxation_);
  body_.set_wall_pull(pharynx_first, pharynx_last, pull_depth_ * relaxation_);
  body_.set_teeth(layout.teeth_start, layout.teeth_end);

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
  posture.lower_rest = std::max(voicing_lower_ + open - 0.002 * (compression - 0.5) - press, -0.006);
  posture.upper_rest = std::max(voicing_upper_ + open - press, -0.008);
  posture.chink = 0.0002 + 0.12 * open;
  body_.cord().set_target(posture, 0.015);

  body_.set_lung_pressure(value(voice::Part::lungs, 0), lung_follow_);
  body_.set_pitch(std::max(value(voice::Part::pitch, 0), 40.0));
  double velum_root = std::max(value(voice::Part::velum, 0), 0.0);
  body_.set_velum(velum_root * velum_root);
  // A raised tip (/r/, /sh/, /er/) lifts the blade off the floor of the
  // mouth and opens the space beneath it.
  double floor = value(voice::Part::floor, 0);
  body_.set_under_tongue(0.4 * std::clamp((floor - 0.5) / 1.5, 0.0, 1.0));
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
  if (tick_ % kStepsPerControl == 0) {
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
