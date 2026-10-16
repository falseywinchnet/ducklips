#include "voice/body.hpp"
#include "voice/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace voice {

namespace {

constexpr double kTwoPi = 2.0 * kPi;
constexpr double kMinimumAdmittanceArea = 1.0e-7;
constexpr double kPressingRoot = -0.5;   // signed root area a closure aims for: pressing tissue spreads contact over ~0.5 cm

Index motor_index(Part part, Index coordinate) {
  return static_cast<Index>(part) * 2 + coordinate;
}

double signed_root_area(double area) {
  if (area <= 0.0) {
    return kPressingRoot;
  }
  return std::sqrt(area);
}

double area_from_root(double root) {
  if (root <= 0.0) {
    return 0.0;
  }
  return root * root;
}

// Viscous and thermal attenuation for a wave crossing one section.
double section_gain(double loss, double area) {
  double safe_area = std::max(area, 0.02);
  return 1.0 - loss / std::sqrt(safe_area);
}

constexpr double kClosedNeck = 1.0e-3;  // cm^2: a side cavity's neck counts as shut
constexpr double kWallContactArea = 0.01;
constexpr double kCarryArea = 0.01;
constexpr double kShutArea = 0.02;  // cm^2: a front section this narrow shuts the mouth for the pharynx sink  // cm^2: below this a section is treated as closing, not carrying flow  // cm^2: below this a junction's wall flow fades out

double wall_area(double tube_area) {
  // Flattened (elliptic) cross-section: perimeter about 1.25 x circular.
  double safe_area = std::max(tube_area, 0.05);
  return 1.25 * 2.0 * std::sqrt(kPi * safe_area) * kSectionLength;
}

void prepare_cavity(Cavity& cavity) {
  double inertance = kAirDensity * cavity.neck_length / cavity.neck_area;
  double compliance = cavity.volume / (kAirDensity * kSoundSpeed * kSoundSpeed);
  double resistance = 2.0 * cavity.damping * std::sqrt(inertance / compliance);
  double denominator = 1.0 + kStepSeconds * resistance / inertance +
                       kStepSeconds * kStepSeconds / (inertance * compliance);
  cavity.sink_alpha = (cavity.flow - kStepSeconds / inertance * cavity.pressure) / denominator;
  cavity.sink_beta = kStepSeconds / inertance / denominator;
}

void finish_cavity(Cavity& cavity, double junction_pressure) {
  double compliance = cavity.volume / (kAirDensity * kSoundSpeed * kSoundSpeed);
  cavity.flow = cavity.sink_alpha + cavity.sink_beta * junction_pressure;
  cavity.pressure += kStepSeconds * cavity.flow / compliance;
}

// Radiation load at an open end, solved jointly with the arriving wave
// (backward Euler on the inertive branch) and with the jet's exit loss: air
// leaving a narrow opening into the room gives up its kinetic pressure
// rho U|U| / 2 A^2 (inflow is a contraction and nearly lossless). Without it
// a lip constriction offered no resistance, and /v/ voiced as loudly as a
// vowel. Returns the reflected wave.
double radiate(Radiation& radiation, double arriving, double area) {
  double admittance = std::max(area, 0.0) / kRhoC;
  double load_area = std::max(area, 0.01);
  double radius = std::sqrt(load_area / kPi);
  double inertance = 8.0 * kAirDensity / (3.0 * kPi * kPi * radius);
  double conductance = 9.0 * kPi * kPi * load_area / (128.0 * kRhoC);
  double load = conductance + kStepSeconds / inertance;
  // Flow U = Y (2 arriving - p), p = p_load + k U|U|, U = inertive + load p_load:
  // U + c U|U| = drive, one quadratic.
  double drive = admittance * (2.0 * arriving + radiation.inertive_flow / load) / (1.0 + admittance / load);
  double exit_area = std::max(area, 1.0e-4);
  double kinetic = drive > 0.0 ? kAirDensity / (2.0 * exit_area * exit_area) : 0.0;
  double coefficient = admittance * kinetic / (1.0 + admittance / load);
  double magnitude = 2.0 * std::abs(drive) / (1.0 + std::sqrt(1.0 + 4.0 * coefficient * std::abs(drive)));
  double flow = drive >= 0.0 ? magnitude : -magnitude;
  double load_pressure = (flow - radiation.inertive_flow) / load;
  double pressure = load_pressure + kinetic * flow * std::abs(flow);
  radiation.inertive_flow += kStepSeconds * load_pressure / inertance;
  radiation.radiated_flow = radiation.inertive_flow + conductance * load_pressure;
  return pressure - arriving;
}

// Builds the trapezoidal update of a wall section for one step length.
WallStepper make_wall_stepper(const WallMaterial& material) {
  double mc = material.cover_mass;
  double mm = material.mass_mass;
  double k = material.cover_stiffness;
  double c = material.cover_damping;
  double ks = material.support_stiffness;
  double cs = material.support_damping;
  double system[4][4] = {{0.0, 1.0, 0.0, 0.0},
                         {-k / mc, -c / mc, k / mc, c / mc},
                         {0.0, 0.0, 0.0, 1.0},
                         {k / mm, c / mm, -(k + ks) / mm, -(c + cs) / mm}};
  double half = 0.5 * kStepSeconds;
  // left = I - h/2 A, right = I + h/2 A; invert left by Gauss-Jordan.
  double left[4][8] = {};
  double right[4][4] = {};
  for (Index row = 0; row < 4; ++row) {
    for (Index column = 0; column < 4; ++column) {
      double identity = row == column ? 1.0 : 0.0;
      left[row][column] = identity - half * system[row][column];
      right[row][column] = identity + half * system[row][column];
      left[row][column + 4] = identity;
    }
  }
  for (Index pivot = 0; pivot < 4; ++pivot) {
    Index best = pivot;
    for (Index row = pivot + 1; row < 4; ++row) {
      if (std::abs(left[row][pivot]) > std::abs(left[best][pivot])) {
        best = row;
      }
    }
    for (Index column = 0; column < 8; ++column) {
      std::swap(left[pivot][column], left[best][column]);
    }
    double scale = 1.0 / left[pivot][pivot];
    for (Index column = 0; column < 8; ++column) {
      left[pivot][column] *= scale;
    }
    for (Index row = 0; row < 4; ++row) {
      if (row == pivot) {
        continue;
      }
      double factor = left[row][pivot];
      for (Index column = 0; column < 8; ++column) {
        left[row][column] -= factor * left[pivot][column];
      }
    }
  }
  WallStepper stepper;
  for (Index row = 0; row < 4; ++row) {
    for (Index column = 0; column < 4; ++column) {
      double sum = 0.0;
      for (Index inner = 0; inner < 4; ++inner) {
        sum += left[row][inner + 4] * right[inner][column];
      }
      stepper.step[row][column] = sum;
    }
    stepper.pressure[row] = half * left[row][4 + 1] / mc;
    stepper.mass_input[row] = 2.0 * half * left[row][4 + 3] / mm;
  }
  return stepper;
}

// The fold edge along its length: it moves most at the middle and is held
// at the ends, so it closes like a zipper rather than all at once.
constexpr Index kFoldSlices = 8;
// Normalized to unit mean, so the mean swing equals the mass displacement.
constexpr double kFoldMode[kFoldSlices] = {0.30448, 0.86709, 1.29769, 1.53074, 1.53074, 1.29769, 0.86709, 0.30448};

struct FoldEdge {
  double open_width = 0.0;   // mean positive gap along the length, cm
  double penetration = 0.0;  // mode-weighted mean overlap where the edges press, cm
  double pressed = 0.0;      // fraction of the length in contact
};

FoldEdge fold_edge(double rest, double displacement, double flatness) {
  FoldEdge edge;
  for (Index slice = 0; slice < kFoldSlices; ++slice) {
    double shape = kFoldMode[slice] + flatness * (1.0 - kFoldMode[slice]);
    double gap = rest + displacement * shape;
    if (gap > 0.0) {
      edge.open_width += gap;
    } else {
      edge.penetration -= gap * shape;
      edge.pressed += 1.0;
    }
  }
  edge.open_width /= kFoldSlices;
  edge.penetration /= kFoldSlices;
  edge.pressed /= kFoldSlices;
  return edge;
}

// Pressed folds meet along most of their length at once: the vibration
// shape flattens as the arytenoids squeeze, so closure is abrupt.
double pressed_flatness(double rest_gap) {
  return std::clamp(-rest_gap / 0.03, 0.0, 1.0);
}

// Open width where air can pass: per slice, the narrower of the two margins.
double fold_passage(double lower_rest, double lower_x, double upper_rest, double upper_x, double flatness) {
  double sum = 0.0;
  for (Index slice = 0; slice < kFoldSlices; ++slice) {
    double shape = kFoldMode[slice] + flatness * (1.0 - kFoldMode[slice]);
    double lower = lower_rest + lower_x * shape;
    double upper = upper_rest + upper_x * shape;
    sum += std::max(std::min(lower, upper), 0.0);
  }
  return sum / kFoldSlices;
}

// 4th-order Butterworth low-pass at 10 kHz for a 96 kHz step rate: b0 b1 b2 a1 a2.
constexpr double kAliasFilter[2][5] = {
    {0.0661302228219, 0.132260445644, 0.0661302228219, -1.0155428255, 0.280063716783},
    {0.0838008420835, 0.167601684167, 0.0838008420835, -1.28690544681, 0.622108815145}};

}  // namespace

Body::Body(const BodyOptions& options) : options_(options) {
  random_state_ = 0x9E3779B97F4A7C15ULL ^ (static_cast<std::uint64_t>(options.seed) * 0x100000001B3ULL);
  oral_.count = kOralSections;
  sub_.count = kSubSections;
  nasal_.count = kNasalSections;

  oral_.loss = options_.oral_loss;

  // Trachea: stiff cartilage rings, lossy; bronchi below.
  sub_.wall_material.mass_mass = 0.8;
  sub_.wall_material.support_stiffness = 1.0e6;
  sub_.wall_material.support_damping = 2000.0;
  sub_.loss = 0.0012;
  for (Index index = 0; index < kSubSections; ++index) {
    double from_glottis = (kSubSections - index - 0.5) * kSectionLength;
    double area = 2.5;
    if (from_glottis < 1.5) {
      area = 1.6 + 0.6 * from_glottis;  // conus elasticus widening below
    }
    sub_.area[index] = area;
  }

  // Nose: bony walls under thick mucosa, high loss from the turbinates.
  nasal_.wall_material.mass_mass = 1.0;
  nasal_.wall_material.support_stiffness = 2.0e6;
  nasal_.wall_material.support_damping = 3000.0;
  nasal_.loss = 0.0035;
  for (Index index = 0; index < kNasalSections; ++index) {
    double fraction = (index + 0.5) / kNasalSections;
    double area = 2.2;
    if (fraction < 0.12) {
      area = 1.6;  // nasopharynx above the velum
    } else if (fraction < 0.75) {
      area = 2.6;  // main cavity between the turbinates
    } else {
      area = 2.6 - (fraction - 0.75) / 0.25 * 1.6;  // vestibule to nostrils
    }
    nasal_target_[index] = area;
    nasal_.area[index] = area;
  }
  nasal_.area[0] = 0.0;
  nasal_target_[0] = 0.0;

  // Sinuses: maxillary (both sides), sphenoid, frontal.
  sinuses_[0].neck_area = 0.15;
  sinuses_[0].neck_length = 0.6;
  sinuses_[0].volume = 30.0;
  sinus_junction_[0] = 14;
  sinuses_[1].neck_area = 0.08;
  sinuses_[1].neck_length = 0.5;
  sinuses_[1].volume = 7.0;
  sinus_junction_[1] = 6;
  sinuses_[2].neck_area = 0.1;
  sinuses_[2].neck_length = 1.2;
  sinuses_[2].volume = 5.0;
  sinus_junction_[2] = 20;

  floor_.neck_area = 0.8;
  floor_.neck_length = 0.6;
  floor_.volume = 0.5;
  floor_.damping = 0.25;

  // The lateral pocket starts closed; its junction stays between the floor
  // space and the lips.
  pocket_.neck_area = kClosedNeck;
  glide_gain_ = 1.0 - std::exp(-1.0 / (kStepRate * options_.glide_ms * 0.001));
  pocket_.neck_length = options_.lateral_neck;
  pocket_.volume = options_.lateral_volume;
  pocket_.damping = 0.2;
  options_.lateral_junction = std::clamp(options_.lateral_junction, kFloorJunction + 1, kOralSections - 2);

  // Motors: rest posture, wander per unit effort, wander time.
  Motor& lungs = motors_[motor_index(Part::lungs, 0)];
  lungs.wander_scale = 12.0;
  lungs.wander_seconds = 0.35;
  lungs.span = 1000.0;
  Motor& pitch = motors_[motor_index(Part::pitch, 0)];
  pitch.rest = 120.0;
  pitch.position = 120.0;
  pitch.target = 120.0;
  pitch.wander_scale = 0.8;
  pitch.wander_seconds = 0.25;
  pitch.span = 100.0;
  Motor& gap = motors_[motor_index(Part::glottis, 0)];
  gap.rest = 0.12;
  gap.position = 0.12;
  gap.target = 0.12;
  gap.wander_scale = 0.0015;
  gap.wander_seconds = 0.12;
  gap.span = 0.1;
  Motor& compression = motors_[motor_index(Part::glottis, 1)];
  compression.rest = 0.5;
  compression.position = 0.5;
  compression.target = 0.5;
  compression.wander_scale = 0.03;
  Motor& place = motors_[motor_index(Part::tongue, 0)];
  place.rest = 9.5;
  place.position = 9.5;
  place.target = 9.5;
  place.wander_scale = 0.06;
  place.wander_seconds = 0.18;
  place.span = 4.0;
  Motor& tongue_root = motors_[motor_index(Part::tongue, 1)];
  tongue_root.rest = 1.6;
  tongue_root.position = 1.6;
  tongue_root.target = 1.6;
  tongue_root.wander_scale = 0.015;
  tongue_root.wander_seconds = 0.15;
  Part constrictions[2] = {Part::tip, Part::dorsum};
  for (Index which = 0; which < 2; ++which) {
    Motor& where = motors_[motor_index(constrictions[which], 0)];
    where.rest = which == 0 ? 15.8 : 12.0;
    where.position = where.rest;
    where.target = where.rest;
    where.wander_scale = 0.04;
    Motor& root = motors_[motor_index(constrictions[which], 1)];
    root.rest = std::sqrt(static_cast<double>(kReleased));
    root.position = root.rest;
    root.target = root.rest;
    root.wander_scale = 0.012;
    root.wander_seconds = 0.12;
  }
  Motor& lip_root = motors_[motor_index(Part::lips, 0)];
  lip_root.rest = 1.6;
  lip_root.position = 1.6;
  lip_root.target = 1.6;
  lip_root.wander_scale = 0.015;
  lip_root.wander_seconds = 0.12;
  Motor& protrusion = motors_[motor_index(Part::lips, 1)];
  protrusion.wander_scale = 0.02;
  Motor& velum = motors_[motor_index(Part::velum, 0)];
  velum.wander_scale = 0.01;
  Motor& jaw = motors_[motor_index(Part::jaw, 0)];
  jaw.rest = kReleased;
  jaw.position = kReleased;
  jaw.target = kReleased;
  jaw.wander_scale = 0.01;
  jaw.span = 8.0;
  Motor& floor = motors_[motor_index(Part::floor, 0)];
  floor.rest = 0.5;
  floor.position = 0.5;
  floor.target = 0.5;
  floor.wander_scale = 0.05;
  floor.span = 2.0;
  // The pocket opening follows the tongue's shape; it is not a muscle held
  // against effort, and with wander it stood slightly ajar between laterals.
  Motor& lateral = motors_[motor_index(Part::lateral, 0)];
  lateral.wander_scale = 0.0;
  lateral.span = 0.5;
  for (Index index = 0; index < kPartCount * 2; ++index) {
    motors_[index].omega = 40.0;
  }

  // Piriform recesses: two blind pockets either side of the larynx, about
  // 2 cm deep; their quarter-wave zeros carve the spectrum near 4-5 kHz.
  Index piriform_lengths[2] = {5, 6};
  for (Index side = 0; side < 2; ++side) {
    Tube& pocket = piriform_[side];
    pocket.count = piriform_lengths[side];
    pocket.loss = 0.0015;
    for (Index index = 0; index < pocket.count; ++index) {
      double fraction = (index + 0.5) / pocket.count;
      pocket.area[index] = options_.piriform_scale * (0.9 - 0.5 * fraction);
    }
    pocket.stepper = make_wall_stepper(pocket.wall_material);
  }
  oral_.stepper = make_wall_stepper(oral_.wall_material);
  sub_.stepper = make_wall_stepper(sub_.wall_material);
  nasal_.stepper = make_wall_stepper(nasal_.wall_material);
  // Fixed-shape tubes: their caches never change.
  refresh_cache(sub_, 0, kSubSections - 1);
  refresh_cache(nasal_, 0, kNasalSections - 1);
  refresh_cache(piriform_[0], 0, piriform_[0].count - 1);
  refresh_cache(piriform_[1], 0, piriform_[1].count - 1);
  shape_tract();
  for (Index index = 0; index < kOralSections; ++index) {
    oral_.area[index] = oral_target_[index];
    oral_.area_middle[index] = oral_target_[index];
  }
  tension_glide_.reset(tension_);
  tremor_phase_ = random_uniform() * kTwoPi;
}

double Body::random_uniform() {
  random_state_ ^= random_state_ << 13U;
  random_state_ ^= random_state_ >> 7U;
  random_state_ ^= random_state_ << 17U;
  return static_cast<double>(random_state_ >> 11U) * (1.0 / 9007199254740992.0);
}

double Body::random_gauss() {
  double sum = random_uniform() + random_uniform() + random_uniform() + random_uniform();
  return (sum - 2.0) * 1.7320508075688772;
}

double Body::motor_value(Part part, Index coordinate) const {
  return motors_[motor_index(part, coordinate)].position;
}

void Body::apply(const PartEvent& event) {
  double move = std::max(static_cast<double>(event.move), 0.004);
  double omega = 4.5 / move;
  Index first = motor_index(event.part, 0);
  Index second = motor_index(event.part, 1);
  switch (event.part) {
    case Part::lungs:
    case Part::pitch:
    case Part::floor:
    case Part::lateral:
      motors_[first].target = event.a;
      motors_[first].omega = omega;
      break;
    case Part::jaw:
      motors_[first].target = event.a;
      motors_[second].target = event.b;
      motors_[first].omega = omega;
      motors_[second].omega = omega;
      break;
    case Part::glottis:
      motors_[first].target = event.a;
      motors_[second].target = event.b;
      motors_[first].omega = omega;
      motors_[second].omega = omega;
      break;
    case Part::tongue:
    case Part::tip:
    case Part::dorsum:
      motors_[first].target = event.a;
      motors_[second].target = signed_root_area(event.b);
      motors_[first].omega = omega;
      motors_[second].omega = omega;
      break;
    case Part::lips:
      motors_[first].target = signed_root_area(event.a);
      motors_[second].target = event.b;
      motors_[first].omega = omega;
      motors_[second].omega = omega;
      break;
    case Part::velum:
      motors_[first].target = std::sqrt(std::max(static_cast<double>(event.a), 0.0));
      motors_[first].omega = omega;
      break;
    case Part::body:
    case Part::apex:
    case Part::mouth:
    case Part::spread:
    case Part::larynx:
    case Part::angle:
    case Part::groove:
      // (organ positions belong to the articulatory body)
    case Part::count:
      break;
  }
}

void Body::settle() {
  for (Index index = 0; index < kPartCount * 2; ++index) {
    motors_[index].position = motors_[index].target;
    motors_[index].velocity = 0.0;
  }
  ear_.correction = 0.0;
  tension_ = motors_[motor_index(Part::pitch, 0)].target / options_.reference_hz;
  lung_pressure_ = motor_value(Part::lungs, 0) * kPascal;
  fold_rest_gap_ = motor_value(Part::glottis, 0);
  compression_ = std::clamp(motor_value(Part::glottis, 1), 0.0, 1.0);
  shape_tract();
  for (Index index = 0; index < kOralSections; ++index) {
    oral_.area[index] = oral_target_[index];
    oral_.area_middle[index] = oral_target_[index];
  }
  tension_glide_.reset(tension_);
  lung_glide_.reset(lung_pressure_);
  gap_glide_.reset(fold_rest_gap_);
  compression_glide_.reset(compression_);
  velum_glide_.reset(velum_glide_.goal);
  velum_area_ = velum_glide_.goal;
  floor_glide_.reset(floor_glide_.goal);
  floor_.volume = floor_glide_.goal;
  pocket_glide_.reset(pocket_glide_.goal);
  pocket_.neck_area = std::max(pocket_glide_.goal, kClosedNeck);
}

void Body::wander_motor(Motor& motor, double effort) {
  double decay = std::exp(-kControlSeconds / motor.wander_seconds);
  double spread = motor.wander_scale * options_.hold_noise * effort;
  motor.wander = motor.wander * decay + spread * std::sqrt(1.0 - decay * decay) * random_gauss();
}

void Body::control_tick() {
  for (Index index = 0; index < kPartCount * 2; ++index) {
    Motor& motor = motors_[index];
    // Wander lives in the muscle command; the part's own mass smooths it.
    double pull = motor.omega * motor.omega * (motor.target + motor.wander - motor.position);
    double drag = 2.0 * motor.omega * motor.velocity;
    motor.velocity += kControlSeconds * (pull - drag);
    motor.position += kControlSeconds * motor.velocity;
    double effort = 0.4 + std::abs(motor.target - motor.rest) / motor.span;
    wander_motor(motor, effort);
  }

  // Lungs: the alveolar pressure follows the respiratory muscles.
  lung_glide_.goal = motor_value(Part::lungs, 0) * kPascal;

  // Folds: the speaker sets cricothyroid tension from the intended pitch and
  // what the ear has learned; the muscle responds with a short delay.
  double intended = std::max(motor_value(Part::pitch, 0), 40.0);
  tension_target_ = intended / options_.reference_hz * std::exp(ear_.correction);
  double muscle = 1.0 - std::exp(-kControlSeconds / 0.03);
  tension_ += (tension_target_ - tension_) * muscle;
  tension_glide_.goal = tension_;
  gap_glide_.goal = motor_value(Part::glottis, 0);
  compression_glide_.goal = std::clamp(motor_value(Part::glottis, 1), 0.0, 1.0);

  shape_tract();
}

void Body::shape_tract() {
  TractShape shape;
  shape.tongue_place = motor_value(Part::tongue, 0);
  shape.tongue_area = std::max(area_from_root(motor_value(Part::tongue, 1)), 0.05);
  shape.tip_place = motor_value(Part::tip, 0);
  shape.tip_root = motor_value(Part::tip, 1);
  shape.dorsum_place = motor_value(Part::dorsum, 0);
  shape.dorsum_root = motor_value(Part::dorsum, 1);
  shape.lip_root = motor_value(Part::lips, 0);
  shape.protrusion = motor_value(Part::lips, 1);
  shape.teeth_gap = motor_value(Part::jaw, 0);
  shape.jaw_closure = motor_value(Part::jaw, 1);
  shape.free_front = options_.free_front;
  TractLayout layout = tract_areas(shape, oral_target_);
  double teeth_start = layout.teeth_start;
  double teeth_end = layout.teeth_end;

  double velum_root = std::max(motor_value(Part::velum, 0), 0.0);
  velum_glide_.goal = velum_root * velum_root;
  floor_glide_.goal = std::max(motor_value(Part::floor, 0), 0.2);
  pocket_glide_.goal = std::max(motor_value(Part::lateral, 0), 0.0);

  // Where does the jet form? Narrowest place in the mouth; the noise arises
  // where the jet strikes the wall downstream, or the teeth.
  Index narrowest = -1;
  double narrowest_area = 0.6;
  for (Index index = 6; index < kOralSections; ++index) {
    if (oral_target_[index] < narrowest_area && oral_target_[index] > 0.0) {
      narrowest_area = oral_target_[index];
      narrowest = index;
    }
  }
  // The jet leaves the constriction laminar and turns turbulent over the
  // next few millimetres (growth over 0.4 cm); its turbulence then decays
  // over about a centimetre, and where it strikes the incisors (an
  // obstacle) it is several times stronger. A source at the exit itself
  // drove the cavity behind the constriction and filled /s/ with 1-3 kHz.
  // Weights change continuously as the tongue moves.
  jet_narrowest_ = narrowest;
  // Is the mouth shut (a stop's closure)? Only then does the pharynx sink act.
  double least = oral_target_[kOralSections - 1];
  for (Index index = kPharynxLast; index < kOralSections; ++index) {
    least = std::min(least, oral_target_[static_cast<std::size_t>(index)]);
  }
  mouth_shut_ = least < kShutArea;
  for (Index index = 0; index < kOralSections; ++index) {
    jet_weights_[static_cast<std::size_t>(index)] = 0.0;
  }
  if (narrowest >= 0) {
    double exit_x = (narrowest + 1) * kSectionLength;
    for (Index junction = narrowest; junction <= kOralSections - 2; ++junction) {
      double x = (junction + 1) * kSectionLength;
      double distance = std::max(x - exit_x, 0.0);
      double teeth_cover = std::clamp((x - teeth_start) / kSectionLength + 0.5, 0.0, 1.0) *
                           std::clamp((teeth_end - x) / kSectionLength + 0.5, 0.0, 1.0);
      // The incisors are an obstacle only to a jet that still holds
      // together when it reaches them: a sibilant's jet a centimetre away
      // does; a velar's, 3 cm away, has spread and slowed (a free jet's
      // centreline decays within a few slit heights).
      double coherence = std::exp(-std::max(distance - options_.obstacle_reach, 0.0) / 0.5);
      double obstacle = 1.0 + 3.0 * teeth_cover * coherence;
      double growth = 1.0 - std::exp(-distance / options_.jet_growth);
      jet_weights_[static_cast<std::size_t>(junction)] = 0.35 * growth * std::exp(-distance / options_.jet_decay) * obstacle;
    }
    if (narrowest >= kOralSections - 3) {
      // Labiodental: the jet leaves between the upper incisors and the lip.
      jet_weights_[static_cast<std::size_t>(kOralSections - 2)] += options_.labiodental;
    }
  }
}

// One junction between sections a and b. The wall at the junction draws
// flow alpha + beta p_a; air expanding into the wider side loses the jet's
// kinetic pressure (Borda-Carnot, solved exactly); a turbulence source acts
// in series. All of it is one closed-form solve.
void Body::step_tube_junctions(Tube& tube, Index first, Index last) {
  for (Index junction = first; junction <= last; ++junction) {
    step_junction(tube, junction, 0.0, 0.0);
  }
}

// One scattering junction. A side branch at the junction (wall motion, a
// side cavity) takes flow alpha + beta p from the node on side a. Returns the
// pressure on side a.
double Body::step_junction(Tube& tube, Index junction, double side_alpha, double side_beta) {
  double area_a = tube.area[junction];
  double area_b = tube.area[junction + 1];
  double arriving_a = tube.right[junction] * tube.gain[junction];
  double arriving_b = tube.left[junction + 1] * tube.gain[junction + 1];
  double admittance_a = std::max(area_a, kMinimumAdmittanceArea) / kRhoC;
  double impedance_b = kRhoC / std::max(area_b, kMinimumAdmittanceArea);
  // The wall's flow reaches the junction only through the opening itself:
  // with the junction shut and the wall still coupled, the release of a
  // closure let the wall dump 35 cm^3/s through a 0.00005 cm^2 gap (whose
  // own limit is 0.16 cm^3/s), and every stop and affricate release
  // clicked (owner: "a lot of clicks").
  double wall_alpha = 0.0;
  double wall_beta = 0.0;
  wall_sink(tube, junction, wall_alpha, wall_beta);
  double opening = std::min(area_a, area_b);
  double reach = std::clamp(opening / kWallContactArea, 0.0, 1.0);
  double sink_alpha = side_alpha + reach * wall_alpha;
  double sink_beta = side_beta + reach * wall_beta;
  double node = admittance_a + sink_beta;
  // Kinetic loss k U|U| solved exactly with the linear part Z U = drive:
  // no lag, so a near-closed opening cannot ring numerically.
  double drive = (2.0 * admittance_a * arriving_a - sink_alpha) / node - 2.0 * arriving_b + tube.noise[junction];
  double linear = 1.0 / node + impedance_b;
  // The jet loses its kinetic energy once, where it expands into a wider
  // space (Borda-Carnot: rho/2 (1/A_up - 1/A_down)^2); a contraction or a
  // uniform channel is nearly lossless. Charging rho/(2 A^2) at every
  // narrow junction made a constriction N sections long cost N times the
  // pressure: the /i/ channel after a /k/ ate 6000 dyn/cm^2 and choked the
  // aspiration to 200 cm^3/s.
  double upstream = std::max(drive >= 0.0 ? area_a : area_b, 1.0e-4);
  double downstream = std::max(drive >= 0.0 ? area_b : area_a, 1.0e-4);
  double expansion = upstream < downstream ? 1.0 / upstream - 1.0 / downstream : 0.0;
  double kinetic = 0.5 * kAirDensity * expansion * expansion;
  double magnitude = 2.0 * std::abs(drive) / (linear + std::sqrt(linear * linear + 4.0 * kinetic * std::abs(drive)));
  double flow = drive >= 0.0 ? magnitude : -magnitude;
  if (area_a <= 0.0 || area_b <= 0.0) {
    flow = 0.0;
  }
  double pressure_a = (2.0 * admittance_a * arriving_a - sink_alpha - flow) / node;
  double pressure_b = 2.0 * arriving_b + impedance_b * flow;
  tube.next_left[junction] = pressure_a - arriving_a;
  tube.next_right[junction + 1] = pressure_b - arriving_b;
  tube.pressure[junction] = pressure_a;
  tube.flow[junction] = flow;
  return pressure_a;
}

void Body::step_glottis() {
  Fold& fold = fold_;
  double below = sub_.right[kSubSections - 1];
  double above = oral_.left[0];
  double sub_area = sub_.area[kSubSections - 1];
  double supra_area = std::max(oral_.area[0], 0.05);
  double chink = options_.chink * (1.0 - compression_) + 0.002;
  double membranous = 2.0 * options_.fold_length *
                      fold_passage(fold_rest_gap_ + options_.fold.convergence, fold.lower_x, fold_rest_gap_,
                                   fold.upper_x, options_.pressed_flattening * pressed_flatness(fold_rest_gap_));
  double glottal_area = membranous + chink;
  double effective = sub_area * supra_area / (sub_area + supra_area);
  double drive = 2.0 * (below - above);
  // Transglottal pressure = kinetic (Bernoulli) loss + viscous slit loss,
  // and the incident waves set it: kinetic U^2 + resistive U = drive.
  constexpr double kShearViscosity = 1.86e-4;  // poise
  double kinetic = kAirDensity / (2.0 * glottal_area * glottal_area);
  double slit_length = options_.fold_length;
  double slit_depth = options_.viscous_depth;
  double viscous = 12.0 * kShearViscosity * slit_depth * slit_length * slit_length /
                   (glottal_area * glottal_area * glottal_area);
  double resistive = kRhoC / effective + viscous;
  double magnitude = 2.0 * std::abs(drive) /
                     (resistive + std::sqrt(resistive * resistive + 4.0 * kinetic * std::abs(drive)));
  double flow = drive >= 0.0 ? magnitude : -magnitude;
  sub_.next_left[kSubSections - 1] = below - kRhoC / sub_area * flow;
  oral_.next_right[0] = above + kRhoC / supra_area * flow;
  fold.flow = flow;
  fold.area = glottal_area;
  fold.sub_pressure = below + sub_.next_left[kSubSections - 1];
  fold.supra_pressure = above + oral_.next_right[0];

  // Aspiration: turbulence just above the folds. Its strength follows the
  // jet's dynamic pressure, which the lungs bound, once the flow is
  // turbulent (Reynolds number above onset).
  // A jet cannot be faster than its pressure drop allows.
  double drop = std::abs(fold_.sub_pressure - fold_.supra_pressure) + 1.0;
  double velocity = std::min(std::abs(flow) / glottal_area, std::sqrt(2.0 * drop / kAirDensity));
  double diameter = std::sqrt(4.0 * glottal_area / kPi);
  double reynolds = velocity * diameter / kAirViscosity;
  double onset = 1200.0 / std::max(reynolds, 1.0);
  double turbulent = std::max(1.0 - onset * onset, 0.0);
  // Dipole turbulence: sound pressure grows about as velocity cubed, so a
  // slow jet (open folds behind an oral constriction) hisses far less than
  // a fast one. Normalized at 2000 cm/s.
  double dynamic = 0.5 * kAirDensity * velocity * velocity * turbulent * (velocity / 2000.0);
  // The same broadband turbulence as in the mouth (see step_noise): centred
  // near a Strouhal number of 0.2 on the slit height (area over the fold
  // length), rising 12 dB/octave below a third of that corner and falling
  // 6 dB/octave above. A near-white source made /h/ a 6-8 kHz hiss.
  double glottal_height = std::max(glottal_area / options_.fold_length, 0.01);
  double corner = std::clamp(0.2 * velocity / glottal_height, 500.0, 8000.0);
  double corner_gain = 1.0 - std::exp(-2.0 * kPi * corner / kStepRate);
  double floor_gain = 1.0 - std::exp(-2.0 * kPi * corner / (3.0 * kStepRate));
  double white = random_gauss();
  aspiration_low_ += floor_gain * (white - aspiration_low_);
  double once = white - aspiration_low_;
  aspiration_lower_ += floor_gain * (once - aspiration_lower_);
  aspiration_band_ += corner_gain * ((once - aspiration_lower_) - aspiration_band_);
  double unit = std::sqrt((2.0 - corner_gain) / corner_gain);
  // A wider jet carries a larger turbulent region: scale by its diameter.
  // Adducted, vibrating folds make a cleaner jet than parted ones (h,
  // aspiration, breath), but not a clean one: each open phase is turbulent.
  // Held at 0.06 here, bdl's vowels carried 12 dB more 4.5-7.5 kHz than ours
  // (-31.9 vs -44.2 dB re 0.3-4.5 kHz) - the owner's "speaking into a tube".
  double opening = std::clamp((fold_rest_gap_ - 0.01) / 0.08, 0.0, 1.0);
  double parted = std::max(opening, 0.06);
  double hiss = 0.10 * options_.aspiration * parted * (1.0 + 12.0 * parted * parted);
  // Its own term, so /h/ and the stop aspiration (fitted to bdl) keep their
  // level: fades out as the folds part.
  double breath = 0.10 * options_.voiced_breath * (1.0 - opening);
  oral_.noise[1] = (hiss + breath) * dynamic * std::sqrt(glottal_area) * unit * aspiration_band_;
}

void Body::step_velum() {
  Index junction = kVelumJunction;
  double area_a = oral_.area[junction];
  double area_b = oral_.area[junction + 1];
  double area_n = std::max(velum_area_, 0.0);
  double arriving_a = oral_.right[junction] * oral_.gain[junction];
  double arriving_b = oral_.left[junction + 1] * oral_.gain[junction + 1];
  double arriving_n = nasal_.left[0];
  double admittance_a = std::max(area_a, kMinimumAdmittanceArea) / kRhoC;
  double admittance_b = std::max(area_b, kMinimumAdmittanceArea) / kRhoC;
  double admittance_n = area_n / kRhoC;
  double sink_alpha = 0.0;
  double sink_beta = 0.0;
  wall_sink(oral_, junction, sink_alpha, sink_beta);
  double total = admittance_a + admittance_b + admittance_n + sink_beta;
  double pressure = (2.0 * (admittance_a * arriving_a + admittance_b * arriving_b + admittance_n * arriving_n) -
                     sink_alpha) / total;
  oral_.next_left[junction] = pressure - arriving_a;
  oral_.next_right[junction + 1] = pressure - arriving_b;
  nasal_.next_right[0] = pressure - arriving_n;
  oral_.pressure[junction] = pressure;
  oral_.flow[junction] = admittance_a * (2.0 * arriving_a - pressure);
}

// The hypopharynx junction: the tract continues, and both piriform
// recesses open into it.
void Body::step_piriform() {
  Index junction = kPiriformJunction;
  double area_a = oral_.area[junction];
  double area_b = oral_.area[junction + 1];
  double arriving_a = oral_.right[junction] * oral_.gain[junction];
  double arriving_b = oral_.left[junction + 1] * oral_.gain[junction + 1];
  double admittance_a = std::max(area_a, kMinimumAdmittanceArea) / kRhoC;
  double admittance_b = std::max(area_b, kMinimumAdmittanceArea) / kRhoC;
  double sink_alpha = 0.0;
  double sink_beta = 0.0;
  wall_sink(oral_, junction, sink_alpha, sink_beta);
  double sum = 2.0 * (admittance_a * arriving_a + admittance_b * arriving_b);
  double total = admittance_a + admittance_b + sink_beta;
  double arriving_side[2] = {0.0, 0.0};
  for (Index side = 0; side < 2; ++side) {
    const Tube& pocket = piriform_[side];
    double admittance = pocket.area[0] / kRhoC;
    arriving_side[side] = pocket.left[0];
    sum += 2.0 * admittance * arriving_side[side];
    total += admittance;
  }
  double pressure = (sum - sink_alpha) / total;
  oral_.next_left[junction] = pressure - arriving_a;
  oral_.next_right[junction + 1] = pressure - arriving_b;
  oral_.pressure[junction] = pressure;
  oral_.flow[junction] = admittance_a * (2.0 * arriving_a - pressure);
  for (Index side = 0; side < 2; ++side) {
    Tube& pocket = piriform_[side];
    pocket.next_right[0] = pressure - arriving_side[side];
    step_tube_junctions(pocket, 0, pocket.count - 2);
    Index last = pocket.count - 1;
    pocket.next_left[last] = pocket.right[last] * pocket.gain[last];
  }
}

// An oral junction with a side cavity (the floor space under the tongue, the
// lateral pocket). A cavity whose neck has closed holds no flow and follows
// the junction pressure, so it reopens without a jolt. (The floor junction
// had its own solve without the expansion loss; adding it changed nothing
// measurable: identification 37.8% both ways over three seeds.)
void Body::step_cavity_junction(Cavity& cavity, Index junction) {
  if (cavity.neck_area <= kClosedNeck) {
    cavity.flow = 0.0;
    cavity.pressure = step_junction(oral_, junction, 0.0, 0.0);
    return;
  }
  prepare_cavity(cavity);
  double pressure = step_junction(oral_, junction, cavity.sink_alpha, cavity.sink_beta);
  finish_cavity(cavity, pressure);
}

double Body::step_lips() {
  Index last = kOralSections - 1;
  double arriving = oral_.right[last] * oral_.gain[last];
  oral_.next_left[last] = radiate(lips_, arriving, oral_.area[last]);
  return lips_.radiated_flow;
}

double Body::step_nostrils() {
  Index last = kNasalSections - 1;
  double arriving = nasal_.right[last] * nasal_.gain[last];
  nasal_.next_left[last] = radiate(nostrils_, arriving, nasal_.area[last]);
  return nostrils_.radiated_flow;
}

// Wall pieces, in two halves around the junction solve. Before: the edge
// forces from neighbouring pieces and the support are frozen for the step,
// and the cover velocity is written as free part + admittance x pressure.
// The junction then solves pressure and wall flow together. After: the whole
// section state advances with the solved pressure.
void Body::prepare_walls(Tube& tube) {
  const WallMaterial& material = tube.wall_material;
  const WallStepper& stepper = tube.stepper;
  WallState& wall = tube.wall;
  Index last = tube.count - 2;
  for (Index junction = 0; junction <= last; ++junction) {
    double cover_left = junction > 0 ? wall.cover_y[junction - 1] : wall.cover_y[junction];
    double cover_right = junction < last ? wall.cover_y[junction + 1] : wall.cover_y[junction];
    double mass_left = junction > 0 ? wall.mass_y[junction - 1] : wall.mass_y[junction];
    double mass_right = junction < last ? wall.mass_y[junction + 1] : wall.mass_y[junction];
    wall.cover_force[junction] =
        -material.edge_stiffness * (2.0 * wall.cover_y[junction] - cover_left - cover_right);
    wall.mass_force[junction] =
        -material.mass_edge_stiffness * (2.0 * wall.mass_y[junction] - mass_left - mass_right) +
        material.support_stiffness * wall.support_y[junction];
    wall.previous_pressure[junction] = tube.pressure[junction];
    const std::array<double, 4>& row = stepper.step[1];
    double free_velocity = row[0] * wall.cover_y[junction] + row[1] * wall.cover_v[junction] +
                           row[2] * wall.mass_y[junction] + row[3] * wall.mass_v[junction] +
                           stepper.pressure[1] * (wall.previous_pressure[junction] + 2.0 * wall.cover_force[junction]) +
                           stepper.mass_input[1] * wall.mass_force[junction];
    wall.cover_free[junction] = free_velocity;
  }
}

void Body::wall_sink(const Tube& tube, Index junction, double& alpha, double& beta) const {
  // A pharynx that grows holds more air: it draws flow L dA/dt (half from
  // each section meeting at this junction), and a shrinking one pushes air
  // out. Without this an actively enlarged pharynx drew no air, and a
  // voiced stop's closure could not keep the folds vibrating. Only the wide
  // pharynx takes part: the parametric tongue does not conserve its own
  // volume elsewhere, and at a narrow release the smoothed rate lags the
  // area and forced flow into a closed section (a ringing burst).
  double narrow = std::min(tube.area[static_cast<std::size_t>(junction)], tube.area[static_cast<std::size_t>(junction + 1)]);
  // Only while the mouth is shut: that is when an enlarged pharynx has to
  // draw its air through the folds (the voiced stop's voice bar). With the
  // mouth open the parametric tongue's volume is not conserved, and the
  // sink turned every move of the body into a gust: the /i/ -> dark /l/ of
  // "Kilo" pushed 900-1000 cm^3/s through the moving channel (glottis
  // passing 10-60), which hissed and clicked (owner: "kilo is terrible").
  double weight = options_.sink_closed_only ? sink_weight_ : 1.0;
  if (&tube == &oral_ && weight > 0.0 && junction >= kPharynxFirst && junction <= options_.sink_last && narrow > 0.3) {
    alpha += weight * 0.5 * kSectionLength *
             (tube.area_rate[static_cast<std::size_t>(junction)] + tube.area_rate[static_cast<std::size_t>(junction + 1)]);
  }
  double coupling = options_.wall_coupling * tube.wall_material.coupling;
  if (coupling <= 0.0) {
    return;
  }
  double area = tube.wall_surface[junction] * coupling;
  alpha += area * tube.wall.cover_free[junction];
  beta += area * tube.stepper.pressure[1];
}

void Body::finish_walls(Tube& tube) {
  const WallStepper& stepper = tube.stepper;
  WallState& wall = tube.wall;
  Index last = tube.count - 2;
  for (Index junction = 0; junction <= last; ++junction) {
    double state[4] = {wall.cover_y[junction], wall.cover_v[junction], wall.mass_y[junction], wall.mass_v[junction]};
    double drive = wall.previous_pressure[junction] + tube.pressure[junction] + 2.0 * wall.cover_force[junction];
    double next[4] = {0.0, 0.0, 0.0, 0.0};
    for (Index row = 0; row < 4; ++row) {
      next[row] = stepper.step[row][0] * state[0] + stepper.step[row][1] * state[1] +
                  stepper.step[row][2] * state[2] + stepper.step[row][3] * state[3] +
                  stepper.pressure[row] * drive + stepper.mass_input[row] * wall.mass_force[junction];
    }
    wall.cover_y[junction] = next[0];
    wall.cover_v[junction] = next[1];
    wall.mass_y[junction] = next[2];
    wall.mass_v[junction] = next[3];
  }
}

void Body::step_fold() {
  Fold& fold = fold_;
  double tremor = 1.0 + 0.003 * std::sin(tremor_phase_);
  double tension = tension_glide_.value * tremor;
  // Cricothyroid stretch lengthens the fold: it becomes thinner and
  // lighter as well as stiffer (Ishizaka-Flanagan Q: mass / Q, stiffness x Q).
  const FoldMaterial& material = options_.fold;
  // Mass falls and stiffness rises with stretch; pitch follows tension.
  // Pitch goes as sqrt(stiffness / mass), so the two exponents sum to 2.
  double lightening = std::pow(tension, options_.fold_mass_exponent);
  const double kLowerMass = material.lower_mass / lightening;
  const double kUpperMass = material.upper_mass / lightening;
  const double kBodyMass = material.body_mass / lightening;
  double stretch_stiffening = tension * tension / lightening;
  // Parting the folds for a voiceless sound (posterior cricoarytenoid
  // pulling the arytenoids out and back) also stretches and tenses them,
  // so the vibration dies as they part instead of flapping on: soft
  // physiological folds otherwise swung 400 cm^3/s pulses through a
  // 0.2 cm^2 glottis between "is" and "Papa" and cut them off as the lips
  // closed (a pop).
  stretch_stiffening *= 1.0 + options_.abduction_stiffening * std::max(fold_rest_gap_ - options_.abduction_onset, 0.0);
  double lower_stiffness = material.lower_stiffness * stretch_stiffening;
  double upper_stiffness = material.upper_stiffness * stretch_stiffening;
  double shear_stiffness = material.shear_stiffness * stretch_stiffening;
  double body_stiffness = material.body_stiffness * stretch_stiffening;
  double thinning = 1.0 / std::sqrt(tension);
  // Cover loss is the same whether the folds rest together or apart. (A cut
  // of up to 85% for folds resting apart, meant to keep breathy voice alive,
  // let them swing chaotically while closing in after every voiceless
  // consonant: the next vowel had 35% shimmer against the speaker's 7%.)
  // Parted folds lose energy faster than they gain it from the flow: held
  // apart for a voiceless sound they do not flap. (Done by stiffening
  // instead, the flapping stopped but every adduction and abduction
  // re-tuned the folds and the pitch lurched.)
  double looseness = 1.0 + options_.abduction_damping * std::max(fold_rest_gap_ - options_.abduction_onset, 0.0);
  double lower_damping = 2.0 * material.lower_damping * looseness * std::sqrt(kLowerMass * lower_stiffness);
  double upper_damping = 2.0 * material.upper_damping * looseness * std::sqrt(kUpperMass * upper_stiffness);
  double body_damping = 2.0 * material.body_damping * std::sqrt(kBodyMass * body_stiffness);
  constexpr double kHardening = 100.0;

  double flatness = options_.pressed_flattening * pressed_flatness(fold_rest_gap_);
  FoldEdge lower_edge = fold_edge(fold_rest_gap_ + options_.fold.convergence, fold.lower_x, flatness);
  FoldEdge upper_edge = fold_edge(fold_rest_gap_, fold.upper_x, flatness);
  double lower_gap = lower_edge.open_width;
  double upper_gap = upper_edge.open_width;
  double sub_pressure = fold.sub_pressure;
  double supra_pressure = fold.supra_pressure;
  double lower_pressure = supra_pressure;
  double upper_pressure = supra_pressure;
  constexpr double kShut = 1.0e-5;
  if (lower_gap <= kShut) {
    lower_pressure = sub_pressure;
    upper_pressure = supra_pressure;
  } else if (upper_gap <= kShut) {
    lower_pressure = sub_pressure;
    upper_pressure = sub_pressure;
  } else if (lower_gap > upper_gap) {
    double ratio = upper_gap / lower_gap;
    lower_pressure = sub_pressure - (sub_pressure - supra_pressure) * ratio * ratio;
    upper_pressure = supra_pressure;
  }
  double lower_force = lower_pressure * options_.fold_length * options_.lower_thickness * thinning;
  double upper_force = upper_pressure * options_.fold_length * options_.upper_thickness * thinning;

  double lower_stretch = fold.lower_x - fold.body_x;
  double upper_stretch = fold.upper_x - fold.body_x;
  double lower_spring = lower_stiffness * lower_stretch * (1.0 + kHardening * lower_stretch * lower_stretch);
  double upper_spring = upper_stiffness * upper_stretch * (1.0 + kHardening * upper_stretch * upper_stretch);
  double lower_drag = lower_damping * (fold.lower_v - fold.body_v);
  double upper_drag = upper_damping * (fold.upper_v - fold.body_v);
  double shear = shear_stiffness * (fold.lower_x - fold.upper_x);

  // Contact, slice by slice along the length.
  double lower_contact = 0.0;
  if (lower_edge.pressed > 0.0) {
    double depth = lower_edge.penetration;
    lower_contact = material.contact_scale * lower_stiffness * depth * (1.0 + 500.0 * depth * depth) -
                    lower_edge.pressed * 2.0 * material.contact_damping * std::sqrt(kLowerMass * material.contact_scale * lower_stiffness) *
                        fold.lower_v;
  }
  double upper_contact = 0.0;
  if (upper_edge.pressed > 0.0) {
    double depth = upper_edge.penetration;
    upper_contact = material.contact_scale * upper_stiffness * depth * (1.0 + 500.0 * depth * depth) -
                    upper_edge.pressed * 2.0 * material.contact_damping * std::sqrt(kUpperMass * material.contact_scale * upper_stiffness) *
                        fold.upper_v;
  }
  fold.contact = std::max(lower_edge.penetration, upper_edge.penetration);

  double lower_total = lower_force - lower_spring - lower_drag - shear + lower_contact;
  double upper_total = upper_force - upper_spring - upper_drag + shear + upper_contact;
  double body_total = lower_spring + lower_drag + upper_spring + upper_drag - body_stiffness * fold.body_x -
                      body_damping * fold.body_v;
  fold.lower_v += kStepSeconds * lower_total / kLowerMass;
  fold.upper_v += kStepSeconds * upper_total / kUpperMass;
  fold.body_v += kStepSeconds * body_total / kBodyMass;
  fold.lower_x += kStepSeconds * fold.lower_v;
  fold.upper_x += kStepSeconds * fold.upper_v;
  fold.body_x += kStepSeconds * fold.body_v;
}

void Body::step_noise() {
  if (jet_narrowest_ < 0) {
    return;
  }
  Index narrowest = jet_narrowest_;
  double narrowest_area = oral_.area[static_cast<std::size_t>(narrowest)];
  for (Index index = std::max(narrowest - 3, 6); index <= std::min(narrowest + 3, kOralSections - 1); ++index) {
    if (oral_.area[static_cast<std::size_t>(index)] < narrowest_area) {
      narrowest_area = oral_.area[static_cast<std::size_t>(index)];
      narrowest = index;
    }
  }
  double area = std::max(narrowest_area, 1.0e-4);
  double flow = oral_.flow[static_cast<std::size_t>(std::min(narrowest, kOralSections - 2))];
  Index upstream = std::max(narrowest - 2, 0);
  Index downstream = std::min(narrowest + 2, kOralSections - 2);
  // At the lips the jet exits into the room, at atmospheric pressure.
  double outside = narrowest >= kOralSections - 3 ? 0.0 : oral_.pressure[static_cast<std::size_t>(downstream)];
  double drop = std::abs(oral_.pressure[static_cast<std::size_t>(upstream)] - outside) + 1.0;
  double velocity = std::min(std::abs(flow) / area, std::sqrt(2.0 * drop / kAirDensity));
  double diameter = std::sqrt(4.0 * area / kPi);
  double reynolds = velocity * diameter / kAirViscosity;
  double onset = 1800.0 / std::max(reynolds, 1.0);
  double turbulent = std::max(1.0 - onset * onset, 0.0);
  // Fade as the opening widens, so the source never switches abruptly.
  double fade = std::clamp((0.6 - area) / 0.3, 0.0, 1.0);
  // Dipole turbulence: pressure grows about as velocity cubed.
  double dynamic = 0.5 * kAirDensity * velocity * velocity * turbulent * fade * fade * (velocity / 2000.0);
  // Jet turbulence is broadband, its energy centred near a Strouhal number
  // of about 0.2 on the jet's thin dimension: oral constrictions are slits
  // about 1.2 cm wide, so the height (area / 1.2) sets it, not the
  // equivalent round diameter. As a dipole (the jet striking teeth and
  // wall) the source rises about 12 dB/octave up to a third of that corner,
  // is flat to the corner, and falls 6 dB/octave above. The cavities in front of the constriction then shape it: the /s/
  // peak near 5 kHz, /sh/ near 2.5 kHz, the velar burst near 1.5-2 kHz.
  // (A narrow band at 1.2 v/d, clamped at 7 kHz, had put every sibilant's
  // and burst's energy at 6-8 kHz.)
  double slit_height = std::max(area / 1.2, 0.01);
  double corner = std::clamp(options_.strouhal * velocity / slit_height, 1500.0, options_.corner_ceiling);
  double corner_gain = 1.0 - std::exp(-2.0 * kPi * corner / kStepRate);
  double floor_gain = 1.0 - std::exp(-2.0 * kPi * corner / (options_.frication_band * kStepRate));
  // Unit variance whatever the corner, so the v^3 law alone sets the power.
  // One pole a = 1 - g passes g / (2 - g) of white noise's power; two
  // identical poles pass g^4 (1 + a^2) / (1 - a^2)^3.
  double pole = 1.0 - corner_gain;
  double passed = options_.frication_stages >= 2
                      ? std::pow(corner_gain, 4.0) * (1.0 + pole * pole) / std::pow(1.0 - pole * pole, 3.0)
                      : corner_gain / (2.0 - corner_gain);
  double unit = 1.0 / std::sqrt(passed);
  double strength = 0.2 * options_.frication * dynamic * std::sqrt(area) * unit;
  // At the lips the narrowest section can be the last one, which has no
  // junction downstream; the jet still sounds at the lip junction. (Starting
  // the loop at the narrowest section injected nothing whenever the minimum
  // flickered onto the last section, and /f/ came out as sparse bursts.)
  Index first_junction = std::min(narrowest, kOralSections - 2);
  // The eddies along a jet are independent beyond a few millimetres (about
  // the slit height), so each junction sounds its own turbulence. One
  // signal copied to every junction from the exit to the teeth (1.5 cm for
  // /sh/) interfered with itself and smeared the front cavity's resonance:
  // driven by white noise, /sh/ came out flat within 4 dB from 1 to 8 kHz.
  // Independent sources add in power, so the strength is rescaled to keep
  // the coherent sum's level: sum(w) / sqrt(sum(w^2)).
  double weight_sum = 0.0;
  double weight_squares = 0.0;
  for (Index junction = first_junction; junction <= kOralSections - 2; ++junction) {
    double weight = jet_weights_[static_cast<std::size_t>(junction)];
    weight_sum += weight;
    weight_squares += weight * weight;
  }
  if (weight_squares <= 0.0) {
    return;
  }
  if (options_.jet_incoherent) {
    double level = strength * weight_sum / std::sqrt(weight_squares);
    for (Index junction = first_junction; junction <= kOralSections - 2; ++junction) {
      double weight = jet_weights_[static_cast<std::size_t>(junction)];
      if (weight > 0.0) {
        double shaped = shape_jet_noise(jet_noise_[static_cast<std::size_t>(junction)], floor_gain, corner_gain);
        oral_.noise[static_cast<std::size_t>(junction)] += weight * level * shaped;
      }
    }
    return;
  }
  double shared = shape_jet_noise(jet_noise_[0], floor_gain, corner_gain);
  for (Index junction = first_junction; junction <= kOralSections - 2; ++junction) {
    double weight = jet_weights_[static_cast<std::size_t>(junction)];
    if (weight > 0.0) {
      oral_.noise[static_cast<std::size_t>(junction)] += weight * strength * shared;
    }
  }
}

// One turbulence generator: white noise with the band below a third of the
// corner removed (two stages, 12 dB/octave) and low-passed at the corner.
// Above its peak, turbulence falls faster than one pole: with -6 dB/oct
// here and +6 from radiation, /s/ rose to 8 kHz and /sh/ stayed flat, no
// front-cavity peak standing out (owner: "all instances of s mangled",
// "she" heard as "fe"). A second stage gives -12 dB/oct.
double Body::shape_jet_noise(JetNoise& state, double floor_gain, double corner_gain) {
  double white = random_gauss();
  state.low += floor_gain * (white - state.low);
  double once = white - state.low;
  state.lower += floor_gain * (once - state.lower);
  state.band += corner_gain * ((once - state.lower) - state.band);
  state.second += corner_gain * (state.band - state.second);
  return options_.frication_stages >= 2 ? state.second : state.band;
}

// Boundary-layer loss. Air rubbing on the walls (viscosity) and trading
// heat with them damps a wave in proportion to about the square root of
// frequency; the per-section gain alone is frequency-independent, and it
// left F3-F5 of /a/ as strong as F1 with bandwidths of 80-110 Hz (human
// 150-400 Hz): sharp, strong high resonances that the owner heard as a
// "computer sound". A wave moves one section per step, so smoothing it
// against its neighbours in space is a symmetric filter in time: a loss
// growing as 2 e (1 - cos w) with no delay, so no formant moves. (A
// one-pole low-pass per section would have lengthened the tract by about
// 13% at the same loss.)
void Body::diffuse_waves(Tube& tube) {
  double strength = options_.wave_diffusion;
  if (strength <= 0.0) {
    return;
  }
  // Smooth the power-carrying wave (pressure x root area), and only where
  // the tube is smooth: across a closure or a step in area the mix would
  // move power between unlike sections (smoothing raw pressure diverged at
  // once).
  Index last = tube.count - 1;
  for (Index index = 0; index <= last; ++index) {
    std::size_t here = static_cast<std::size_t>(index);
    double root = std::sqrt(std::max(tube.area[here], kMinimumAdmittanceArea));
    diffuse_root_[here] = root;
    diffuse_scratch_[here] = tube.right[here] * root;
    diffuse_scratch_left_[here] = tube.left[here] * root;
  }
  for (Index index = 1; index < last; ++index) {
    std::size_t here = static_cast<std::size_t>(index);
    double low = std::min(std::min(tube.area[here - 1], tube.area[here]), tube.area[here + 1]);
    double high = std::max(std::max(tube.area[here - 1], tube.area[here]), tube.area[here + 1]);
    if (low < 0.1 || high > 1.5 * low) {
      continue;
    }
    double right = diffuse_scratch_[here] +
                   strength * (diffuse_scratch_[here - 1] - 2.0 * diffuse_scratch_[here] + diffuse_scratch_[here + 1]);
    double left = diffuse_scratch_left_[here] + strength * (diffuse_scratch_left_[here - 1] -
                                                              2.0 * diffuse_scratch_left_[here] +
                                                              diffuse_scratch_left_[here + 1]);
    tube.right[here] = right / diffuse_root_[here];
    tube.left[here] = left / diffuse_root_[here];
  }
}

void Body::refresh_cache(Tube& tube, Index first, Index last) {
  for (Index index = first; index <= last; ++index) {
    tube.gain[static_cast<std::size_t>(index)] = section_gain(tube.loss, tube.area[static_cast<std::size_t>(index)]);
  }
  for (Index junction = std::max(first - 1, 0); junction <= std::min(last, tube.count - 2); ++junction) {
    tube.wall_surface[static_cast<std::size_t>(junction)] = wall_area(
        0.5 * (tube.area[static_cast<std::size_t>(junction)] + tube.area[static_cast<std::size_t>(junction + 1)]));
  }
}

void Body::advance_glides() {
  // Two smoothing stages of BodyOptions::glide_ms each. At 0.6 ms, about 7% of every 1 ms control step survived as a
  // 1 kHz staircase on the areas, the fold gap, the tension and the lung
  // pressure - a ripple under every movement (owner: "sounds of a motor").
  double kGain = glide_gain_;
  for (Index index = 0; index < kOralSections; ++index) {
    double goal = oral_target_[index];
    double middle = oral_.area_middle[index] + kGain * (goal - oral_.area_middle[index]);
    double area = oral_.area[index] + kGain * (middle - oral_.area[index]);
    area = area < 1.0e-5 ? 0.0 : area;
    // Carry the flow, not the wave difference, through an area change: the
    // section's flow is A (right - left) / rho c, so held waves made any
    // steady flow grow with the area. A constriction sliding back from /i/
    // to dark /l/ ("Kilo") widened the sections in front of it from 0.2 to
    // about 3 cm^2 and turned 100 cm^3/s from the glottis into 900-1000 at
    // the lips, hissing through the moving channel. Pressure (right + left)
    // is kept as it is.
    double previous = oral_.area[index];
    if (options_.flow_carry && previous > kCarryArea && area > kCarryArea) {
      std::size_t here = static_cast<std::size_t>(index);
      double pressure = oral_.right[here] + oral_.left[here];
      double difference = (oral_.right[here] - oral_.left[here]) * previous / area;
      oral_.right[here] = 0.5 * (pressure + difference);
      oral_.left[here] = 0.5 * (pressure - difference);
    }
    // Articulators move below 20 Hz, but their targets step at the 1 kHz
    // control rate: two 50 Hz stages keep that ripple out of the volume
    // flow (it buzzed at 1, 2, 3 ... kHz) while the volume moved is kept.
    constexpr double kRateGain = 1.0 - 0.99673;  // exp(-2 pi 50 Hz / 96 kHz)
    double rate = (area - oral_.area[index]) / kStepSeconds;
    oral_.area_rate_middle[index] += kRateGain * (rate - oral_.area_rate_middle[index]);
    oral_.area_rate[index] += kRateGain * (oral_.area_rate_middle[index] - oral_.area_rate[index]);
    oral_.area_middle[index] = middle;
    oral_.area[index] = area;
  }
  tension_glide_.advance(kGain);
  lung_glide_.advance(kGain);
  velum_glide_.advance(kGain);
  gap_glide_.advance(kGain);
  compression_glide_.advance(kGain);
  floor_glide_.advance(kGain);
  pocket_glide_.advance(kGain);
  // The sink fades in and out over about 5 ms as the mouth shuts and opens.
  constexpr double kSinkGain = 1.0 - 0.99792;  // exp(-1 / (96 kHz * 5 ms))
  sink_weight_ += kSinkGain * ((mouth_shut_ ? 1.0 : 0.0) - sink_weight_);
  lung_pressure_ = lung_glide_.value;
  velum_area_ = std::max(velum_glide_.value, 0.0);
  fold_rest_gap_ = gap_glide_.value;
  compression_ = compression_glide_.value;
  floor_.volume = floor_glide_.value;
  pocket_.neck_area = std::max(pocket_glide_.value, kClosedNeck);
  refresh_cache(oral_, 0, kOralSections - 1);
}

void Body::step() {
  advance_glides();
  nasal_.area[0] = velum_area_;
  nasal_.area[1] = 0.5 * (velum_area_ + nasal_target_[2]);
  refresh_cache(nasal_, 0, 2);

  for (Index index = 0; index < kOralSections; ++index) {
    oral_.noise[index] = 0.0;
  }
  step_noise();
  prepare_walls(oral_);
  prepare_walls(sub_);
  prepare_walls(nasal_);
  prepare_walls(piriform_[0]);
  prepare_walls(piriform_[1]);

  step_glottis();
  // Oral junctions, with the floor space under the tongue at its junction.
  step_tube_junctions(oral_, 0, kPiriformJunction - 1);
  step_piriform();
  step_tube_junctions(oral_, kPiriformJunction + 1, kVelumJunction - 1);
  step_velum();
  step_tube_junctions(oral_, kVelumJunction + 1, kFloorJunction - 1);
  step_cavity_junction(floor_, kFloorJunction);
  Index lateral = options_.lateral_junction;
  step_tube_junctions(oral_, kFloorJunction + 1, lateral - 1);
  step_cavity_junction(pocket_, lateral);
  step_tube_junctions(oral_, lateral + 1, kOralSections - 2);
  double lip_flow = step_lips();

  // Below the folds: the lungs feed the bronchi and trachea.
  // The bronchial tree: an ideal pressure source for slow changes (the lungs
  // hold their pressure), a mostly absorbing end for sound (thousands of
  // branching airways), with a weak reflection that leaves the subglottal
  // resonances present but damped.
  double returning = sub_.left[0];
  lung_return_low_ += (returning - lung_return_low_) * (kTwoPi * 40.0 * kStepSeconds);
  // Airway resistance of the bronchial tree: heavy flow (open folds in a
  // voiceless sound) lowers the pressure that reaches the trachea.
  double out_flow = (sub_.right[0] - returning) * sub_.area[0] / kRhoC;
  airway_flow_ += (out_flow - airway_flow_) * (kTwoPi * 20.0 * kStepSeconds);
  constexpr double kAirwayResistance = 2.0;  // dyn s / cm^5
  double driving = lung_pressure_ - kAirwayResistance * std::max(airway_flow_, 0.0);
  sub_.next_right[0] = -0.35 * returning - 0.65 * lung_return_low_ + driving;
  step_tube_junctions(sub_, 0, kSubSections - 2);

  // Nose, with its sinuses hanging off three junctions.
  for (Index which = 0; which < 3; ++which) {
    prepare_cavity(sinuses_[which]);
  }
  for (Index junction = 0; junction <= kNasalSections - 2; ++junction) {
    double area_a = nasal_.area[junction];
    double area_b = nasal_.area[junction + 1];
    double arriving_a = nasal_.right[junction] * nasal_.gain[junction];
    double arriving_b = nasal_.left[junction + 1] * nasal_.gain[junction + 1];
    double admittance_a = std::max(area_a, kMinimumAdmittanceArea) / kRhoC;
    double admittance_b = std::max(area_b, kMinimumAdmittanceArea) / kRhoC;
    double sink_alpha = 0.0;
    double sink_beta = 0.0;
    wall_sink(nasal_, junction, sink_alpha, sink_beta);
    Index sinus = -1;
    for (Index which = 0; which < 3; ++which) {
      if (sinus_junction_[which] == junction) {
        sinus = which;
        sink_alpha += sinuses_[which].sink_alpha;
        sink_beta += sinuses_[which].sink_beta;
      }
    }
    double pressure = (2.0 * (admittance_a * arriving_a + admittance_b * arriving_b) - sink_alpha) /
                      (admittance_a + admittance_b + sink_beta);
    nasal_.next_left[junction] = pressure - arriving_a;
    nasal_.next_right[junction + 1] = pressure - arriving_b;
    nasal_.pressure[junction] = pressure;
    if (sinus >= 0) {
      finish_cavity(sinuses_[sinus], pressure);
    }
  }
  double nose_flow = step_nostrils();

  Tube* tubes[5] = {&oral_, &sub_, &nasal_, &piriform_[0], &piriform_[1]};
  for (Index which = 0; which < 5; ++which) {
    Tube& tube = *tubes[which];
    for (Index index = 0; index < tube.count; ++index) {
      tube.right[index] = tube.next_right[index];
      tube.left[index] = tube.next_left[index];
    }
    finish_walls(tube);
  }
  diffuse_waves(oral_);
  diffuse_waves(nasal_);
  // The piriform recesses are mucosa-lined too: undamped, their pole between
  // the two quarter-wave zeros rang as a steady line near 5.06 kHz through
  // every vowel (bandwidth 148 Hz, unmoved by lips, tongue or walls; it
  // vanished when the recesses were shrunk).
  diffuse_waves(piriform_[0]);
  diffuse_waves(piriform_[1]);
  step_fold();
  tremor_phase_ += kTwoPi * 5.3 * kStepSeconds;
  time_ += kStepSeconds;

  // The ear: one glottal cycle per upward crossing of the mean opening.
  double gap = fold_.area;
  ear_.mean_gap += (gap - ear_.mean_gap) * (kStepSeconds / 0.02);
  bool crossing = ear_.previous_gap < ear_.mean_gap && gap >= ear_.mean_gap;
  ear_.previous_gap = gap;
  if (crossing) {
    double period = time_ - ear_.last_crossing;
    ear_.last_crossing = time_;
    bool voicing = lung_pressure_ > 1500.0 && fold_rest_gap_ < 0.06;
    if (period > 1.0 / 450.0 && period < 1.0 / 50.0 && voicing) {
      ear_.heard_hz = 1.0 / period;
      double intended = std::max(motor_value(Part::pitch, 0), 40.0);
      double error = std::log(intended / ear_.heard_hz);
      // A percept far from the intention (a skipped or doubled cycle) is not
      // trusted; real corrections are small and slow.
      if (std::abs(error) < 0.3) {
        ear_.correction += options_.ear_gain * error;
        ear_.correction = std::clamp(ear_.correction, -0.4, 0.4);
      }
    }
  }

  (void)lip_flow;
  (void)nose_flow;
}

void Body::render(float* output, Index sample_count) {
  for (Index sample = 0; sample < sample_count; ++sample) {
    if (control_countdown_ <= 0) {
      control_tick();
      control_countdown_ = kSamplesPerControl;
    }
    --control_countdown_;
    // Decimate the radiated flow from the step rate: average the sub-steps
    // (a box filter) so sharp closures do not alias into the audio band.
    // Anti-alias the radiated flow at the step rate (two Butterworth
    // biquads, 10 kHz) before keeping every second step.
    double radiated = 0.0;
    for (Index sub_step = 0; sub_step < kStepsPerSample; ++sub_step) {
      step();
      double raw_skin = 0.0;
      for (Index junction = 0; junction < 20; ++junction) {
        raw_skin += oral_.wall_surface[junction] * oral_.wall.cover_v[junction];
      }
      double value = lips_.radiated_flow + nostrils_.radiated_flow + options_.skin * raw_skin;
      for (Index stage = 0; stage < 2; ++stage) {
        const double* coefficient = kAliasFilter[stage];
        double& first = alias_state_[static_cast<std::size_t>(stage * 2)];
        double& second = alias_state_[static_cast<std::size_t>(stage * 2 + 1)];
        double middle = value - coefficient[3] * first - coefficient[4] * second;
        value = coefficient[0] * middle + coefficient[1] * first + coefficient[2] * second;
        second = first;
        first = middle;
      }
      radiated = value;
    }
    double total = radiated;
    double pressure = (total - output_previous_) * kAudioRate;
    output_previous_ = total;
    highpass_state_ += (pressure - highpass_state_) * (kTwoPi * 30.0 / kAudioRate);
    double heard = pressure - highpass_state_;
    double roll = 1.0 - std::exp(-kTwoPi * options_.microphone_hz / kAudioRate);
    microphone_state_ += roll * (heard - microphone_state_);
    output[sample] = static_cast<float>(microphone_state_ * 1.0e-6);
  }
}

BodyTrace Body::trace() const {
  BodyTrace result;
  result.pitch_heard = ear_.heard_hz;
  result.glottal_area = fold_.area;
  result.glottal_flow = fold_.flow;
  result.sub_pressure = fold_.sub_pressure;
  result.tension = tension_;
  result.lip_flow = lips_.radiated_flow;
  result.nose_flow = nostrils_.radiated_flow;
  result.lip_area = oral_.area[kOralSections - 1];
  result.velum_area = velum_area_;
  result.narrowest_area = 1.0e9;
  for (Index index = 6; index < kOralSections; ++index) {
    if (oral_.area[static_cast<std::size_t>(index)] < result.narrowest_area) {
      result.narrowest_area = oral_.area[static_cast<std::size_t>(index)];
      result.narrowest_place = (index + 0.5) * kSectionLength;
    }
  }
  return result;
}

}  // namespace voice
