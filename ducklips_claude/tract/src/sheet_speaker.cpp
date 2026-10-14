#include "tract/sheet_speaker.hpp"

#include "voice/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace tract {

namespace {

constexpr double kAir = 1.204;             // kg/m^3
constexpr double kKinematic = 1.5e-5;      // m^2/s, air
constexpr double kTractCm = voice::kOralSections * voice::kSectionLength;  // 17.5 cm
constexpr double kNeutralArea = 3.0;       // cm^2
constexpr double kStrouhal = 0.2;

double signed_root(double area) {
  return area >= 0.0 ? std::sqrt(area) : -std::sqrt(-area);
}

double area_of_root(double root) {
  return root <= 0.0 ? 0.0 : root * root;
}

SheetBand make_band(const char* name, double length, double radius, const char* material, double mass,
                    double elasticity, double offset, double anchor, double support, int points, const char* role) {
  SheetBand band;
  band.name = name;
  band.length = length;
  band.radius = radius;
  band.material = material;
  band.mass = mass;
  band.elasticity = elasticity;
  band.offset = offset;
  band.anchor = anchor;
  band.support = support;
  band.points = points;
  band.role = role;
  return band;
}

}  // namespace

SheetDesign speech_sheet_design(const SheetVoice& voice) {
  int modules = voice.modules;
  SheetDesign design;
  design.cord_taper = voice.cord_taper;
  design.cord_contact = voice.cord_contact;
  design.cord_loss = voice.cord_loss;
  design.control_rate = voice.control_rate;
  design.wall_loss = voice.wall_loss;
  for (int k = 0; k < 3; ++k) {
    design.air_taps[k] = voice.oral_taps[k];
    design.subglottal_taps[k] = voice.subglottal_taps[k];
  }
  design.pressure = 0.0;
  design.wander = 0.0;
  design.variation = 0.0;
  design.acoustic_spacing = true;
  design.axial_air = true;
  design.radiating = true;
  design.orifice = true;
  design.cord_length = 0.014;
  design.jet_loss = voice.jet_loss;
  design.lung_resistance = voice.lung_resistance;  // the front end's airway resistance
  design.lung_corner = voice.lung_corner;
  // Astra's subglottal sections (42c7337). Their supports narrow the airway
  // to 0.25-0.39 cm^2 (a person's trachea: 2-3 cm^2); subglottal_narrowing
  // scales those offsets (1: Astra's, 0: the bands' own radii, 1.5-3.8 cm^2).
  design.bands.push_back(make_band("01_Inlet_expansion", 0.05, 0.009, "soft", 1.0, 1.0,
                                   voice.subglottal_narrowing * -0.0036271115257626475, 0.54,
                                   1.0, 8, "inlet"));
  design.bands.push_back(make_band("02_Lower_elastic_adapter", 0.0175, 0.008, "fold", 0.2, 3.0,
                                   voice.subglottal_narrowing * -0.005313555762881324, 0.53, 0.2, 4, "adapter"));
  design.bands.push_back(make_band("03_Lower_firm_support", 0.0225, 0.011, "firm", 4.7, 0.2,
                                   voice.subglottal_narrowing * -0.008313555762881324,
                                   0.54, 3.2, 4, "firm"));
  design.bands.push_back(make_band("04_Lower_heavy_tissue", 0.0175, 0.007, "soft", 3.5, 1.5,
                                   voice.subglottal_narrowing * -0.004313555762881324,
                                   0.50, 0.2, 4, "soft"));
  // The slit cord: 4 mm thick, 14 mm folds, rest half-gap 0.1 mm (found
  // steady on /a/ /i/ /u/ at 800 Pa with heavy walls: 101-104 Hz, closed
  // 18-27%, MFDR 350-900 cm^3/s/ms; the 3 mm cord of mass 4 sat at 127 Hz
  // and could not go lower).
  design.bands.push_back(make_band("05_Vocal_cord", voice.cord_thickness, 0.0001, "fold", voice.cord_mass,
                                   voice.cord_elasticity, -0.0001, 0.55, voice.cord_support, 8, "cord"));
  double neutral = std::sqrt(kNeutralArea / kPi) / 100.0;  // m
  double section = kTractCm / 100.0 / (3.0 * modules);      // m
  for (int m = 0; m < modules; ++m) {
    std::string index = std::to_string(m + 1);
    std::string elastic = "elastic_" + index;
    std::string support = "support_" + index;
    std::string mass = "mass_" + index;
    design.bands.push_back(make_band(elastic.c_str(), section, neutral, "fold", 30.0 * voice.tract_mass, 3.0, 0.0, 0.5, 0.2, 4,
                                     "elastic"));
    design.bands.push_back(make_band(support.c_str(), section, neutral, "firm", 21.0 * voice.tract_mass, 0.2, 0.0, 0.5, 3.2, 4,
                                     "support"));
    if (m == modules - 1) {
      design.bands.push_back(make_band("lips", section, neutral, "fold", 36.0 * voice.tract_mass, 1.5, 0.0, 0.5, 0.3, 4, "lips"));
    } else {
      design.bands.push_back(make_band(mass.c_str(), section, neutral, "soft", 45.0 * voice.tract_mass, 1.5, 0.0, 0.5, 0.2, 4,
                                       "mass"));
    }
  }
  return design;
}

SheetDesign nasal_sheet_design(const SheetVoice& voice, double area_share, double length_share) {
  int modules = 5;
  SheetDesign design;
  for (int k = 0; k < 3; ++k) {
    design.air_taps[k] = voice.nasal_taps[k];
  }
  design.air_lowpass = voice.nasal_lowpass;
  design.wall_loss = voice.nasal_wall_loss;
  design.acoustic_spacing = true;
  design.axial_air = true;
  design.radiating = true;
  design.jet_loss = true;
  double section = length_share * voice.nasal_length / (3.0 * modules);
  for (int m = 0; m < modules; ++m) {
    for (int k = 0; k < 3; ++k) {
      int index = 3 * m + k;
      double fraction = (index + 0.5) / (3.0 * modules);
      double area = 2.6;
      if (fraction < 0.12) {
        area = 1.6;
      } else if (fraction >= 0.75) {
        area = 2.6 - (fraction - 0.75) / 0.25 * 1.6;
      }
      double radius = std::sqrt(area_share * area / kPi) / 100.0;
      std::string name = "nose_" + std::to_string(index + 1);
      if (k == 0) {
        design.bands.push_back(make_band(name.c_str(), section, radius, "fold", 30.0 * voice.nasal_mass, 3.0, 0.0, 0.5, 0.2, 4, "elastic"));
      } else if (k == 1) {
        design.bands.push_back(make_band(name.c_str(), section, radius, "firm", 21.0 * voice.nasal_mass, 0.2, 0.0, 0.5, 3.2, 4, "support"));
      } else {
        design.bands.push_back(make_band(name.c_str(), section, radius, "soft", 45.0 * voice.nasal_mass, 1.5, 0.0, 0.5, 0.2, 4, "mass"));
      }
    }
  }
  return design;
}

SheetDesign sinus_sheet_design(const SheetVoice& voice, double radius_share) {
  SheetDesign design;
  design.acoustic_spacing = true;
  design.axial_air = true;
  design.bands.push_back(make_band("sinus_neck", 0.012, 0.0018, "firm", 1.0, 1.0, 0.0, 0.5, 2.0, 4, "tract"));
  design.bands.push_back(make_band("sinus_chamber", 0.035, radius_share * voice.sinus_radius, "firm", 2.0, 1.0, 0.0, 0.5, 3.0, 4, "tract"));
  return design;
}

SheetDesign piriform_sheet_design(double depth) {
  SheetDesign design;
  design.acoustic_spacing = true;
  design.axial_air = true;
  double front = std::sqrt(0.9 / kPi) / 100.0;
  double back = std::sqrt(0.4 / kPi) / 100.0;
  design.bands.push_back(make_band("piriform_mouth", 0.5 * depth, front, "soft", 15.0, 1.5, 0.0, 0.5, 0.5, 4, "tract"));
  design.bands.push_back(make_band("piriform_floor", 0.5 * depth, back, "soft", 15.0, 1.5, 0.0, 0.5, 0.5, 4, "tract"));
  return design;
}

SheetDesign lateral_sheet_design() {
  SheetDesign design;
  design.acoustic_spacing = true;
  design.axial_air = true;
  design.bands.push_back(make_band("lateral_neck", 0.005, std::sqrt(0.3 / kPi) / 100.0, "soft", 15.0, 1.5, 0.0, 0.5, 0.5, 4, "tract"));
  design.bands.push_back(make_band("lateral_chamber", 0.015, std::sqrt(1.67 / kPi) / 100.0, "soft", 15.0, 1.5, 0.0, 0.5, 0.5, 4, "tract"));
  return design;
}

SheetDesign under_tongue_sheet_design() {
  SheetDesign design;
  design.acoustic_spacing = true;
  design.axial_air = true;
  design.bands.push_back(make_band("under_entry", 0.012, 0.0025, "soft", 2.0, 1.0, 0.0, 0.5, 0.3, 4, "tract"));
  design.bands.push_back(make_band("floor_pocket", 0.025, 0.0055, "soft", 3.0, 1.5, 0.0, 0.5, 0.2, 4, "tract"));
  return design;
}

SheetSpeaker::SheetSpeaker(const voice::Score& score, std::uint32_t seed, const SheetVoice& voice)
    : score_(score),
      voice_(voice),
      design_(speech_sheet_design(voice)),
      tube_(design_, voice.rate),
      nose_(nasal_sheet_design(voice, voice.split_nose ? 0.55 : 1.0, 1.0), voice.rate),
      sinus_(sinus_sheet_design(voice, 1.0), voice.rate),
      nose_right_(nasal_sheet_design(voice, 0.45, 1.06), voice.rate),
      sinus_right_(sinus_sheet_design(voice, 1.25), voice.rate),
      under_(under_tongue_sheet_design(), voice.rate),
      piriform_left_(piriform_sheet_design(0.018), voice.rate),
      piriform_right_(piriform_sheet_design(0.022), voice.rate),
      lateral_(lateral_sheet_design(), voice.rate),
      random_(seed != 0U ? seed : 1U) {
  per_control_ = std::max(1, static_cast<int>(std::lround(voice_.rate * SomaticHold::kTick)));
  tube_.set_outlet(SheetTube::End::radiating);
  if (voice_.triangular_glottis) {
    tube_.attach_glottis(&glottis_);
    glottis_.set_elongation(voice_.fold_elongation, voice_.fold_elongation_reference,
                             voice_.fold_elongation_release, voice_.fold_elongation_floor);
  }
  if (voice_.tissue_relaxation > 0.0) {
    tube_.set_strain_relaxation(voice_.tissue_relaxation, 0.1, 0.6, seed + 17U);
  }
  for (std::size_t b = 0; b < design_.bands.size(); ++b) {
    if (design_.bands[b].role == "cord") {
      cord_band_ = static_cast<int>(b);
      first_tract_band_ = static_cast<int>(b) + 1;
    }
  }
  first_tract_node_ = tube_.cord_last() + 1;
  band_center_.assign(design_.bands.size(), 0.0);
  double start = 0.0;
  for (std::size_t b = static_cast<std::size_t>(first_tract_band_); b < design_.bands.size(); ++b) {
    band_center_[b] = 100.0 * (start + 0.5 * design_.bands[b].length);
    start += design_.bands[b].length;
  }
  node_center_.assign(static_cast<std::size_t>(tube_.count()), -1.0);
  jet_weights_.assign(static_cast<std::size_t>(tube_.count()), 0.0);
  for (std::size_t b = static_cast<std::size_t>(first_tract_band_); b < design_.bands.size(); ++b) {
    int first = tube_.band_first(static_cast<int>(b));
    int points = tube_.band_points(static_cast<int>(b));
    double length = 100.0 * design_.bands[b].length;
    for (int j = 0; j < points; ++j) {
      node_center_[static_cast<std::size_t>(first + j)] = band_center_[b] - 0.5 * length + (j + 0.5) * length / points;
    }
  }
  // Side branches: the nose at the velum (9.3 cm above the glottis, the
  // front end's junction 25), the sinus off the main nasal cavity, the space
  // under the tongue below the blade (13.5 cm).
  nose_.set_inlet(SheetTube::End::connected);
  nose_.set_outlet(SheetTube::End::radiating);
  sinus_.set_inlet(SheetTube::End::connected);
  sinus_.set_outlet(SheetTube::End::closed);
  nose_right_.set_inlet(SheetTube::End::connected);
  nose_right_.set_outlet(SheetTube::End::radiating);
  sinus_right_.set_inlet(SheetTube::End::connected);
  sinus_right_.set_outlet(SheetTube::End::closed);
  under_.set_inlet(SheetTube::End::connected);
  under_.set_outlet(SheetTube::End::closed);
  skin_gain_ = 1.0 - std::exp(-2.0 * kPi * voice_.skin_corner / voice_.rate);
  pharynx_node_ = tract_node_at(5.0);
  mouth_node_ = tract_node_at(11.0);
  if (voice_.side_tubes) {
    velum_joint_.connect(tube_, tract_node_at(9.3), nose_, 0, voice_.rate);
    sinus_joint_.connect(nose_, nose_.count() * 2 / 5, sinus_, 0, voice_.rate);
    if (voice_.split_nose) {
      velum_right_joint_.connect(tube_, tract_node_at(9.3), nose_right_, 0, voice_.rate);
      sinus_right_joint_.connect(nose_right_, nose_right_.count() / 2, sinus_right_, 0, voice_.rate);
      velum_right_joint_.set_opening(0.0);
    }
    under_joint_.connect(tube_, tract_node_at(voice_.under_tongue_at), under_, 0, voice_.rate);
  }
  velum_joint_.set_opening(0.0);
  under_joint_.set_opening(0.0);
  if (voice_.recesses) {
    for (SheetTube* recess : {&piriform_left_, &piriform_right_}) {
      recess->set_inlet(SheetTube::End::connected);
      recess->set_outlet(SheetTube::End::closed);
    }
    piriform_left_joint_.connect(tube_, tract_node_at(2.2), piriform_left_, 0, voice_.rate);
    piriform_right_joint_.connect(tube_, tract_node_at(2.2), piriform_right_, 0, voice_.rate);
  }
  if (voice_.recesses || voice_.lateral_pocket) {
    lateral_.set_inlet(SheetTube::End::connected);
    lateral_.set_outlet(SheetTube::End::closed);
    lateral_joint_.connect(tube_, tract_node_at(16.0), lateral_, 0, voice_.rate);
    lateral_joint_.set_opening(0.0);
  }

  HoldSettings settings;
  settings.noise = voice_.hold_noise;
  settings.common = voice_.hold_common;
  settings.fastest = voice_.fastest;
  hold_.reset(seed, settings);
  coordinates_.fill(-1);
  // Relaxed postures (origins), the span that costs one unit of effort (wide
  // enough that a full gesture's burst stays within the pools' capacity), the
  // limits and the neural groups (0 tongue, 1 lips and jaw, 2 velum, 3
  // larynx, 4 breath).
  add_coordinate(voice::Part::lungs, 0, 0.0, 600.0, 0.0, 2500.0, 4);
  add_coordinate(voice::Part::pitch, 0, 100.0, 60.0, 50.0, 400.0, 3);
  add_coordinate(voice::Part::glottis, 0, 0.45, 0.6, -0.02, 0.8, 3);
  add_coordinate(voice::Part::tongue, 0, 9.5, 4.0, 2.0, 17.0, 0);
  add_coordinate(voice::Part::tongue, 1, 1.6, 3.0, -1.5, 3.0, 0);
  add_coordinate(voice::Part::tip, 0, 15.8, 1.0, 6.0, 18.0, 0);
  add_coordinate(voice::Part::tip, 1, std::sqrt(static_cast<double>(voice::kReleased)), 3.0, -1.5, 3.0, 0);
  add_coordinate(voice::Part::dorsum, 0, 12.0, 1.0, 6.0, 17.0, 0);
  add_coordinate(voice::Part::dorsum, 1, std::sqrt(static_cast<double>(voice::kReleased)), 3.0, -1.5, 3.0, 0);
  add_coordinate(voice::Part::lips, 0, 1.6, 3.0, -1.5, 3.0, 1);
  add_coordinate(voice::Part::lips, 1, 0.0, 1.0, 0.0, 1.2, 1);
  add_coordinate(voice::Part::jaw, 0, static_cast<double>(voice::kReleased), 8.0, 0.02, 10.0, 1);
  add_coordinate(voice::Part::jaw, 1, 0.0, 1.0, 0.0, 1.0, 1);
  add_coordinate(voice::Part::velum, 0, 0.0, 1.0, 0.0, 2.0, 2);
  add_coordinate(voice::Part::floor, 0, 0.5, 2.0, 0.0, 4.0, 0);
  add_coordinate(voice::Part::lateral, 0, 0.0, 0.5, 0.0, 2.0, 0);
  applied_.assign(score_.parts.size() + 4096, 0);
  // The score's opening postures are where the body starts, held and settled.
  while (next_event_ < score_.parts.size() && score_.parts[next_event_].time <= 0.0) {
    apply(score_.parts[next_event_]);
    ++next_event_;
  }
  hold_.settle();
  for (int i = 0; i < hold_.count(); ++i) {
    previous_[static_cast<std::size_t>(i)] = hold_.value(i);
  }
  per_shape_ = std::max(1, static_cast<int>(std::lround(voice_.rate / voice_.shape_rate)));
  shape_body();
  tube_.settle_geometry();
}

// The first node of the tract section whose centre is nearest `cm` above the glottis.
int SheetSpeaker::tract_node_at(double cm) const {
  std::size_t best = static_cast<std::size_t>(first_tract_band_);
  for (std::size_t b = best; b < design_.bands.size(); ++b) {
    if (std::abs(band_center_[b] - cm) < std::abs(band_center_[best] - cm)) {
      best = b;
    }
  }
  return tube_.band_first(static_cast<int>(best));
}

void SheetSpeaker::add_coordinate(voice::Part part, int which, double origin, double scale, double low, double high,
                                  int group) {
  int index = hold_.add(origin, scale, low, high, group);
  hold_.set_target(index, origin, 40.0);
  coordinates_[static_cast<std::size_t>(2 * static_cast<int>(part) + which)] = index;
}

double SheetSpeaker::held(voice::Part part, int which) const {
  int index = coordinate(part, which);
  if (index < 0) {
    return 0.0;
  }
  double from = previous_[static_cast<std::size_t>(index)];
  return from + blend_ * (hold_.value(index) - from);
}

// A part's new target. Constriction areas are held as signed roots, so
// closing and pressing move smoothly through zero.
void SheetSpeaker::apply(const voice::PartEvent& event) {
  double move = std::max(static_cast<double>(event.move), 0.004);
  double rate = 4.5 / move;
  bool oral = event.part == voice::Part::tongue || event.part == voice::Part::tip ||
              event.part == voice::Part::dorsum || event.part == voice::Part::lips || event.part == voice::Part::jaw;
  if (oral) {
    rate *= voice_.articulator_speed;
  }
  if (event.part == voice::Part::glottis) {
    rate *= voice_.glottis_speed;
    if (event.a < 0.03 && held(voice::Part::glottis, 0) > 0.03) {
      rate *= voice_.voicing_speed;
    }
  }
  if (event.part == voice::Part::velum) {
    rate *= voice_.velum_speed;
  }
  if (event.part == voice::Part::lungs && static_cast<double>(event.a) < held(voice::Part::lungs, 0)) {
    rate *= voice_.breath_release_speed;
  }
  double a = static_cast<double>(event.a);
  double b = static_cast<double>(event.b);
  if ((event.part == voice::Part::tongue || event.part == voice::Part::tip || event.part == voice::Part::dorsum) &&
      b <= 0.0) {
    b = -voice_.closure_press * voice_.closure_press;
  }
  if (event.part == voice::Part::lips && a <= 0.0) {
    a = -voice_.closure_press * voice_.closure_press;
  }
  int first = coordinate(event.part, 0);
  int second = coordinate(event.part, 1);
  // A closure opening: the constriction's root rises from at or below contact.
  double opening_rate = rate;
  if (event.part == voice::Part::tongue || event.part == voice::Part::tip || event.part == voice::Part::dorsum) {
    if (held(event.part, 1) < 0.3 && signed_root(b) > held(event.part, 1)) {
      opening_rate = rate * voice_.release_speed;
    }
  }
  if (event.part == voice::Part::lips && held(event.part, 0) < 0.3 && signed_root(a) > held(event.part, 0)) {
    opening_rate = rate * voice_.release_speed;
  }
  switch (event.part) {
    case voice::Part::lungs:
    case voice::Part::pitch:
    case voice::Part::glottis:
    case voice::Part::floor:
    case voice::Part::lateral:
      hold_.set_target(first, a, rate);
      break;
    case voice::Part::velum:
      hold_.set_target(first, std::sqrt(std::max(a, 0.0)), rate);
      break;
    case voice::Part::jaw:
      hold_.set_target(first, a, rate);
      hold_.set_target(second, b, rate);
      break;
    case voice::Part::tongue:
    case voice::Part::tip:
    case voice::Part::dorsum:
      hold_.set_target(first, a, rate);
      hold_.set_target(second, signed_root(b), opening_rate);
      break;
    case voice::Part::lips:
      hold_.set_target(first, signed_root(a), opening_rate);
      hold_.set_target(second, b, rate);
      break;
    case voice::Part::count:
      break;
  }
}

void SheetSpeaker::shape_body() {
  voice::TractShape shape;
  shape.tongue_place = held(voice::Part::tongue, 0);
  shape.tongue_area = std::max(area_of_root(held(voice::Part::tongue, 1)), 0.05);
  shape.tip_place = held(voice::Part::tip, 0);
  shape.tip_root = held(voice::Part::tip, 1);
  shape.dorsum_place = held(voice::Part::dorsum, 0);
  shape.dorsum_root = held(voice::Part::dorsum, 1);
  shape.lip_root = held(voice::Part::lips, 0);
  shape.protrusion = held(voice::Part::lips, 1);
  shape.teeth_gap = held(voice::Part::jaw, 0);
  shape.jaw_closure = held(voice::Part::jaw, 1);
  std::array<double, voice::kOralSections> areas{};
  voice::TractLayout layout = voice::tract_areas(shape, areas);
  teeth_start_ = layout.teeth_start;
  teeth_end_ = layout.teeth_end;
  target_narrowest_ = 100.0;
  for (std::size_t k = 8; k < areas.size(); ++k) {
    target_narrowest_ = std::min(target_narrowest_, areas[k]);
  }
  // A voiced closure: the mouth sealed beyond 10.5 cm with the folds set to
  // vibrate. The pharyngeal walls let go toward it (tau 50 ms) and return
  // (tau 40 ms).
  double sealed = 100.0;
  for (std::size_t k = static_cast<std::size_t>(10.5 / voice::kSectionLength); k < areas.size(); ++k) {
    sealed = std::min(sealed, areas[k]);
  }
  bool voiced_closure = sealed < 0.02 && held(voice::Part::glottis, 0) < 0.03;
  double shape_seconds = static_cast<double>(per_shape_) / voice_.rate;
  // Behind a voiced oral closure with the nose shut, the throat and cheeks
  // let the trapped air expand them (the larynx lowers, the cheeks give):
  // a volume sink that keeps air crossing the folds, so voicing goes on
  // through the closure. It begins over 10 ms and ends over 15 ms.
  // The walls can only give so far (closure_volume, cm^3): once the cavity
  // has grown by that much, pressure builds again for the release; between
  // closures the walls return (tau 100 ms).
  bool expanding = voiced_closure && held(voice::Part::velum, 0) < 0.05 &&
                   expanded_volume_ < 1.0e-6 * voice_.closure_volume;
  expanded_volume_ += expansion_ * shape_seconds;
  if (!voiced_closure) {
    expanded_volume_ *= std::exp(-shape_seconds / 0.1);
  }
  double expansion_target = expanding ? 1.0e-6 * voice_.closure_expansion : 0.0;
  double expansion_tau = expanding ? 0.010 : 0.015;
  expansion_ += (1.0 - std::exp(-shape_seconds / expansion_tau)) * (expansion_target - expansion_);
  double tau = voiced_closure ? 0.05 : 0.04;
  closure_relaxation_ += (1.0 - std::exp(-shape_seconds / tau)) * ((voiced_closure ? 1.0 : 0.0) - closure_relaxation_);
  double pull = voice_.closure_pull * closure_relaxation_;
  // Every tract section's ring support moves to the narrowest area over its
  // span: a constriction shorter than a section (an alveolar channel is 3-4
  // mm, a section 7.3) is what the air passing it sees. (Taking the area at
  // the centre missed /s/'s channel: 0.21 cm^2, a 5 m/s jet, no hiss.)
  for (std::size_t b = static_cast<std::size_t>(first_tract_band_); b < design_.bands.size(); ++b) {
    double half = 50.0 * design_.bands[b].length;  // cm
    double from = (band_center_[b] - half) / voice::kSectionLength;
    double to = (band_center_[b] + half) / voice::kSectionLength;
    std::size_t first = static_cast<std::size_t>(std::clamp(std::floor(from), 0.0, static_cast<double>(voice::kOralSections - 1)));
    std::size_t last = static_cast<std::size_t>(std::clamp(std::ceil(to) - 1.0, 0.0, static_cast<double>(voice::kOralSections - 1)));
    double area = areas[first];
    for (std::size_t k = first; k <= last; ++k) {
      area = std::min(area, areas[k]);
    }
    double radius = std::sqrt(std::max(area, 0.0) / kPi) / 100.0;
    if (band_center_[b] >= 3.5 && band_center_[b] <= 10.0) {
      radius += pull;
      tube_.set_band_support(static_cast<int>(b), 1.0 - voice_.closure_release * closure_relaxation_);
    }
    tube_.set_band_offset(static_cast<int>(b), radius - design_.bands[b].radius);
  }
  // The arytenoids: every voicing setting (rest half-gap up to 0.002 cm)
  // presses the cord to its voicing adduction; wider settings open it.
  double gap = held(voice::Part::glottis, 0);
  double offset = voice_.pressed + 0.01 * std::max(gap - 0.002, 0.0);
  tube_.set_band_offset(cord_band_, offset);
  double pitch = std::max(held(voice::Part::pitch, 0), 40.0);
  if (voice_.triangular_glottis) {
    // The arytenoid posture as the segment-module speaker's: every voicing
    // setting (half-gap up to 0.002 cm) is the voicing posture; wider
    // settings open from there, with a posterior chink.
    // Settings narrower than a voiceless consonant's (voiced fricatives,
    // 0.012) open by voiced_open_gain: the folds part enough for the flow
    // that feeds the frication while they still vibrate.
    double open = (gap < 0.03 ? voice_.voiced_open_gain : 0.7) * std::max(gap - 0.002, 0.0);
    CordPosture posture;
    posture.lower_rest = std::max(voice_.voicing_lower + open, -0.006);
    posture.upper_rest = std::max(voice_.voicing_upper + open, -0.008);
    posture.chink = 0.0002 + 0.12 * open;
    glottis_.set_target(posture, 0.002);
    glottis_.set_pitch(pitch);
    // Devoicing stiffens the folds (cricothyroid with the abduction of a
    // voiceless consonant, not a breath): the stretch is there at once and
    // relaxes after the glottis closes, so the next vowel starts high and
    // falls. The ear expects it (efference), so it does not correct it away.
    bool devoicing = gap > 0.03 && gap < voice_.devoicing_widest;
    double devoicing_target = devoicing ? voice_.devoicing_raise : 0.0;
    double devoicing_tau = devoicing ? 0.002 : voice_.devoicing_relax;
    devoicing_ += (1.0 - std::exp(-shape_seconds / devoicing_tau)) * (devoicing_target - devoicing_);
    glottis_.set_holding(std::exp(pitch_correction_ + devoicing_), 0.0);
    intended_pitch_ = pitch * std::exp(devoicing_);
  } else {
    double stretch = std::pow(pitch / voice_.base_pitch, 1.0 / voice_.pitch_exponent);
    tube_.set_band_elasticity(cord_band_, std::clamp(stretch, voice_.stretch_low, voice_.stretch_high));
  }
  tube_.set_pressure(voice_.lung_scale * std::max(held(voice::Part::lungs, 0), 0.0));
  // The velum's port (cm^2 from its held root) and, under a raised tip, the
  // opening into the space beneath the blade (the front end's 0.4 cm^2 at
  // full lift).
  if (voice_.side_tubes) {
    double velum = std::max(held(voice::Part::velum, 0), 0.0);
    double port = 1e-4 * voice_.velum_gain * velum * velum;
    if (voice_.split_nose) {
      // The septum divides the passage behind the port: each side takes
      // its share of the opening.
      velum_joint_.set_opening(0.55 * port);
      velum_right_joint_.set_opening(0.45 * port);
    } else {
      velum_joint_.set_opening(port);
    }
    double floor = held(voice::Part::floor, 0);
    under_joint_.set_opening(1e-4 * 0.4 * std::clamp((floor - 0.5) / 1.5, 0.0, 1.0));
  }
  if (voice_.recesses || voice_.lateral_pocket) {
    lateral_joint_.set_opening(1e-4 * std::max(held(voice::Part::lateral, 0), 0.0));
  }
}

double SheetSpeaker::gauss() {
  double sum = 0.0;
  for (int draw = 0; draw < 4; ++draw) {
    random_ ^= random_ << 13U;
    random_ ^= random_ >> 17U;
    random_ ^= random_ << 5U;
    sum += static_cast<double>(random_) / 4294967296.0;
  }
  return (sum - 2.0) * std::sqrt(3.0);
}

// Unit-variance jet noise: white noise with the band below a third of the
// corner removed and two low-pass stages at the corner, the corner a
// Strouhal number of 0.2 on the jet's width (as the front end's jets).
double SheetSpeaker::turbulence(Turbulence& state, double velocity, double width, double lowest, double highest) {
  double step = 1.0 / voice_.rate;
  double corner = std::clamp(kStrouhal * velocity / std::max(width, 1e-5), lowest, highest);
  double corner_gain = 1.0 - std::exp(-2.0 * kPi * corner * step);
  double floor_gain = 1.0 - std::exp(-2.0 * kPi * corner / 3.0 * step);
  double pole = 1.0 - corner_gain;
  double gain_squared = corner_gain * corner_gain;
  double pole_term = 1.0 - pole * pole;
  double passed = gain_squared * gain_squared * (1.0 + pole * pole) / (pole_term * pole_term * pole_term);
  double white = gauss();
  state.low += floor_gain * (white - state.low);
  double once = white - state.low;
  state.lower += floor_gain * (once - state.lower);
  state.band += corner_gain * ((once - state.lower) - state.band);
  state.second += corner_gain * (state.band - state.second);
  return state.second / std::sqrt(passed);
}

// A jet's turbulence: a dipole of strength about rho v^3 sqrt(A) / 2 once
// its Reynolds number (on its diameter) passes ~1800. (The front end's
// cgs source, 0.2 frication (rho v^2/2)(v/2000 cm/s) sqrt(A cm^2), is this
// in SI.)
void SheetSpeaker::sound_glottis() {
  double area = tube_.glottal_area();
  double velocity = std::max(tube_.glottal_velocity(), 0.0);
  if (area <= 1e-10 || velocity <= 0.0) {
    return;
  }
  double diameter = std::sqrt(4.0 * area / kPi);
  double reynolds = velocity * diameter / kKinematic;
  double onset = 1800.0 / std::max(reynolds, 1.0);
  double turbulent = std::max(1.0 - onset * onset, 0.0);
  double strength = voice_.aspiration * 0.5 * kAir * velocity * velocity * velocity * std::sqrt(area) * turbulent;
  tube_.set_glottal_noise(strength * turbulence(glottal_noise_, velocity, diameter, 500.0, 12000.0));
}

void SheetSpeaker::sound_constriction() {
  int count = tube_.count();
  int narrowest = -1;
  double narrowest_radius = 1.0;
  // The epilarynx (below 2.2 cm) is not an oral constriction.
  for (int node = first_tract_node_ + 1; node < count; ++node) {
    if (node_center_[static_cast<std::size_t>(node)] > 2.2 && tube_.radius(node) < narrowest_radius) {
      narrowest_radius = tube_.radius(node);
      narrowest = node;
    }
  }
  jet_velocity_ = 0.0;
  jet_strength_ = 0.0;
  double area = kPi * narrowest_radius * narrowest_radius;
  double fade = std::clamp((0.6 - area * 1e4) / 0.3, 0.0, 1.0);
  if (narrowest < 0 || fade <= 0.0) {
    return;
  }
  // The constriction: the run of nodes narrower than twice the narrowest;
  // the jet's speed from the drop across all of it.
  double wide = 2.0 * narrowest_radius * narrowest_radius;
  int first = narrowest;
  while (first > first_tract_node_ + 1 && tube_.radius(first - 1) * tube_.radius(first - 1) < wide) {
    --first;
  }
  int last = narrowest;
  while (last < count - 1 && tube_.radius(last + 1) * tube_.radius(last + 1) < wide) {
    ++last;
  }
  double upstream = tube_.pressure(std::max(first - 1, first_tract_node_));
  double downstream = last >= count - 1 ? 0.0 : tube_.pressure(last + 1);
  double drop = std::max(upstream - downstream, 0.0) + 1.0;
  double flow = std::abs(tube_.axial_flow(std::min(narrowest, count - 2)));
  double velocity = std::min(flow / std::max(area, 1e-8), std::sqrt(2.0 * drop / kAir));
  double diameter = std::sqrt(4.0 * area / kPi);
  double reynolds = velocity * diameter / kKinematic;
  double onset = 1800.0 / std::max(reynolds, 1.0);
  double turbulent = std::max(1.0 - onset * onset, 0.0);
  // The jet's thin dimension: oral constrictions are slits about 1.2 cm wide.
  double slit = area / 0.012;
  double strength =
      voice_.frication * 0.5 * kAir * velocity * velocity * velocity * std::sqrt(area) * turbulent * fade * fade;
  double noise = turbulence(oral_noise_, velocity, slit, voice_.jet_corner_floor, voice_.jet_corner_ceiling);
  jet_velocity_ = velocity;
  jet_strength_ = strength;
  // Along the jet's path past the exit (the front end's weights). A source
  // at the exit itself drives the cavity behind the constriction and fills
  // /s/ with 1-3 kHz (ours: 6-7 dB too strong at 0.7-1.5 kHz, 16 dB weak at
  // 3-4.5 kHz, heard as /f/).
  double exit_x = node_center_[static_cast<std::size_t>(last)] + 0.5 * 100.0 * tube_.node_dx(last);
  for (int node = last + 1; node < count; ++node) {
    double x = node_center_[static_cast<std::size_t>(node)];
    double distance = std::max(x - exit_x, 0.0);
    double section = 100.0 * tube_.node_dx(node);
    double teeth_cover = std::clamp((x - teeth_start_) / section + 0.5, 0.0, 1.0) *
                         std::clamp((teeth_end_ - x) / section + 0.5, 0.0, 1.0);
    double coherence = std::exp(-std::max(distance - voice_.obstacle_reach, 0.0) / 0.5);
    double obstacle = 1.0 + voice_.obstacle_gain * teeth_cover * coherence;
    double growth = 1.0 - std::exp(-distance / voice_.jet_growth);
    double weight = 0.35 * growth * std::exp(-distance / voice_.jet_decay) * obstacle;
    if (node == count - 1 && last >= count - 3) {
      weight += voice_.labiodental;  // the jet leaves between the upper incisors and the lip
    }
    if (weight > 0.0) {
      if (voice_.jet_monopole) {
        double area_here = kPi * tube_.radius(node) * tube_.radius(node);
        tube_.inject_flow(node, weight * strength * noise * area_here / (kAir * 343.0));
      } else {
        tube_.add_series_pressure(node, weight * strength * noise);
      }
    }
  }
  if (last >= count - 1) {
    // The lips themselves are the narrowest: no section lies past them, and
    // the loop above sounded nothing (our /f/ was 7-9 dB weak at 0.7-7 kHz).
    // The jet still sounds at the lip junction, between the upper incisors
    // and the lip (the front end's).
    tube_.add_series_pressure(count - 1, voice_.labiodental * strength * noise);
  }
}

void SheetSpeaker::trace(std::FILE* out) const {
  int count = tube_.count();
  int narrowest = first_tract_node_ + 1;
  for (int node = first_tract_node_ + 1; node < count; ++node) {
    if (tube_.radius(node) < tube_.radius(narrowest)) {
      narrowest = node;
    }
  }
  double radius = tube_.radius(narrowest);
  double rest = tube_.rest_radius(narrowest);
  std::fprintf(out, "%.3f %6.0f %7.3f %5.0f %7.1f %6.3f %6.0f %6.3f %d target %.3f rest %.3f p %.0f jet %.1f %.1f ear %.3f velum %.3f below %.0f above %.0f\n", time_,
               held(voice::Part::lungs, 0), held(voice::Part::glottis, 0), held(voice::Part::pitch, 0),
               1e6 * tube_.glottal_flow(), 1e4 * tube_.glottal_area(), tube_.pressure(first_tract_node_ - 9),
               1e4 * kPi * radius * radius, narrowest, target_narrowest_, 1e4 * kPi * rest * rest,
               tube_.pressure(narrowest), jet_velocity_, jet_strength_, pitch_correction_,
               std::pow(std::max(held(voice::Part::velum, 0), 0.0), 2.0), 0.1 * glottis_.pressure_below(),
               0.1 * glottis_.pressure_above());
}

// The speaker hears its own pitch at each closure of the folds and corrects
// their tension toward the intended pitch.
void SheetSpeaker::hear_pitch() {
  // The ear listens only while the speaker phonates (glottis adducted,
  // breath up). Opening for a breath, the folds flickered into contact and
  // those "cycles" drove the correction from -0.16 to +0.54: the next
  // phrase began 3 semitones high ("vietnamese").
  bool phonating = held(voice::Part::glottis, 0) < 0.03 && held(voice::Part::lungs, 0) > 250.0;
  if (!phonating) {
    was_closed_ = false;
    last_closure_ = -1;
    return;
  }
  bool closed = glottis_.contact() > 0.0;
  if (closed && !was_closed_) {
    if (last_closure_ >= 0) {
      double heard = voice_.rate / static_cast<double>(tick_ - last_closure_);
      if (heard > 50.0 && heard < 500.0) {
        pitch_correction_ += voice_.pitch_feedback * std::log(intended_pitch_ / heard);
        pitch_correction_ = std::clamp(pitch_correction_, -voice_.pitch_correction_limit, voice_.pitch_correction_limit);
      }
    }
    last_closure_ = tick_;
  }
  was_closed_ = closed;
}

void SheetSpeaker::push(const voice::PartEvent& event) {
  std::vector<voice::PartEvent>& parts = score_.parts;
  std::size_t place = parts.size();
  while (place > next_event_ && parts[place - 1].time > event.time) {
    --place;
  }
  parts.insert(parts.begin() + static_cast<std::ptrdiff_t>(place), event);
  applied_.insert(applied_.begin() + static_cast<std::ptrdiff_t>(place), 0);
}

double SheetSpeaker::step() {
  if (tick_ % per_control_ == 0) {
    // Events in time order; the larynx's (breath, glottis, pitch) act
    // laryngeal_lead earlier. An early laryngeal event may pass ahead of
    // articulator events still waiting, so both cursors are scanned.
    while (next_event_ < score_.parts.size() && score_.parts[next_event_].time <= time_ + voice_.lead) {
      if (!applied_[next_event_]) {
        apply(score_.parts[next_event_]);
      }
      ++next_event_;
    }
    for (std::size_t k = next_event_; k < score_.parts.size(); ++k) {
      const voice::PartEvent& event = score_.parts[k];
      if (event.time > time_ + voice_.lead + voice_.laryngeal_lead) {
        break;
      }
      // The breath only: leading the glottis too moved voicing boundaries
      // against the consonants (identification 26.1 -> 25.0% at 30 ms).
      bool larynx = event.part == voice::Part::lungs;
      if (larynx && !applied_[k]) {
        apply(event);
        applied_[k] = 1;
      }
    }
    for (std::size_t k = next_event_; k < score_.parts.size(); ++k) {
      const voice::PartEvent& event = score_.parts[k];
      if (event.time > time_ + voice_.lead + voice_.voicing_lead) {
        break;
      }
      // A voicing onset: the glottis closing from open toward phonation.
      bool onset = event.part == voice::Part::glottis && event.a < 0.03 && held(voice::Part::glottis, 0) > 0.03;
      if (onset && !applied_[k]) {
        apply(event);
        applied_[k] = 1;
      }
    }
    for (int i = 0; i < hold_.count(); ++i) {
      previous_[static_cast<std::size_t>(i)] = hold_.value(i);
    }
    hold_.tick();
  }
  if (tick_ % per_shape_ == 0) {
    blend_ = static_cast<double>(tick_ % per_control_ + per_shape_) / per_control_;
    blend_ = std::min(blend_, 1.0);
    shape_body();
  }
  sound_glottis();
  sound_constriction();
  if (expansion_ > 0.0) {
    tube_.inject_flow(pharynx_node_, -0.5 * expansion_);
    tube_.inject_flow(mouth_node_, -0.5 * expansion_);
  }
  ++tick_;
  time_ = static_cast<double>(tick_) / voice_.rate;
  if (voice_.side_tubes) {
    velum_joint_.prepare();
    sinus_joint_.prepare();
    under_joint_.prepare();
    if (!nose_.step() || !sinus_.step() || !under_.step()) {
      return 0.0;
    }
    if (voice_.split_nose) {
      velum_right_joint_.prepare();
      sinus_right_joint_.prepare();
      if (!nose_right_.step() || !sinus_right_.step()) {
        return 0.0;
      }
    }
  }
  if (voice_.recesses) {
    piriform_left_joint_.prepare();
    piriform_right_joint_.prepare();
    if (!piriform_left_.step() || !piriform_right_.step()) {
      return 0.0;
    }
  }
  if (voice_.recesses || voice_.lateral_pocket) {
    lateral_joint_.prepare();
    if (!lateral_.step()) {
      return 0.0;
    }
  }
  if (!tube_.step()) {
    return 0.0;
  }
  if (voice_.recesses) {
    piriform_left_joint_.send();
    piriform_right_joint_.send();
  }
  if (voice_.recesses || voice_.lateral_pocket) {
    lateral_joint_.send();
  }
  if (voice_.triangular_glottis) {
    hear_pitch();
  }
  if (voice_.side_tubes) {
    velum_joint_.send();
    sinus_joint_.send();
    if (voice_.split_nose) {
      velum_right_joint_.send();
      sinus_right_joint_.send();
    }
    under_joint_.send();
  }
  // The mouth and the nostrils both radiate, and so do the throat's and
  // cheeks' walls: the tract's volume pulsing outward is a source like the
  // lips' flow (a voiced closure is heard only through them).
  double flow = tube_.radiated_flow();
  double nasal = nose_.radiated_flow() + (voice_.split_nose ? nose_right_.radiated_flow() : 0.0);
  double mouth = voice_.listen == 2 ? 0.0 : flow - previous_flow_;
  double nostrils = voice_.listen == 1 ? 0.0 : nasal - previous_nasal_flow_;
  double skin = 0.0;
  if (voice_.skin_radiation > 0.0) {
    double volume = 0.0;
    // Only the walls' vibration about their supported rest radii radiates;
    // the articulators' own movement (the supports gliding) is far below
    // hearing, and stepped at the shape rate it came out as clicks
    // ("Bravo"'s end: peak 19 dB over the speech, 11 without the skin).
    for (int node = first_tract_node_; node < tube_.count(); ++node) {
      double radius = tube_.radius(node);
      double rest = tube_.rest_radius(node);
      volume += kPi * (radius * radius - rest * rest) * tube_.node_dx(node);
    }
    // Skin and the tissue under it are heavy: they pass the walls' slow
    // pulsing (the voice bar) and not their fast ripples (two poles at
    // skin_corner).
    double raw_flow = (volume - previous_volume_) * voice_.rate;
    skin_low_[0] += skin_gain_ * (raw_flow - skin_low_[0]);
    skin_low_[1] += skin_gain_ * (skin_low_[0] - skin_low_[1]);
    double wall_flow = skin_low_[1];
    skin = voice_.skin_radiation * (wall_flow - previous_wall_flow_);
    previous_volume_ = volume;
    previous_wall_flow_ = wall_flow;
  }
  double sound = (mouth + nostrils + skin) * voice_.rate;
  previous_flow_ = flow;
  previous_nasal_flow_ = nasal;
  return sound;
}

}  // namespace tract
