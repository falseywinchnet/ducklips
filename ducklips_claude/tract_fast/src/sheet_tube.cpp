#include "tract/sheet_tube.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace tract {

namespace {

constexpr double kRho = 1.204;      // kg/m^3 (Astra's engine constants)
constexpr double kSound = 343.0;    // m/s
constexpr int kSectors = 8;         // circumferential samples of the sheet
constexpr double kAngle = 2.0 * 3.141592653589793 / kSectors;
constexpr double kPiValue = 3.141592653589793;
constexpr double kFloor = 0.000025;  // m: contact radius floor
constexpr double kAirTaps[3] = {0.985, 0.010, 0.003};
constexpr double kReferenceHop = 0.0044;  // m: Astra's mean node spacing (26.5 cm / 60)

int power_of_two_at_least(double value) {
  int size = 16;
  while (size < value) {
    size *= 2;
  }
  return size;
}

}  // namespace

SheetMaterial sheet_material(const std::string& name) {
  SheetMaterial material;
  material.name = name;
  if (name == "fold") {
    material.density = 1040.0;
    material.modulus = 22000.0;
    material.thickness = 0.001;
    material.loss = 0.055;
    material.fir[0] = 0.85;
    material.fir[1] = 0.10;
    material.fir[2] = 0.03;
  } else if (name == "firm") {
    material.density = 1100.0;
    material.modulus = 65000.0;
    material.thickness = 0.002;
    material.loss = 0.16;
    material.fir[0] = 0.61;
    material.fir[1] = 0.23;
    material.fir[2] = 0.10;
  }
  return material;
}

bool read_sheet_design(const std::string& path, SheetDesign& design) {
  std::ifstream input(path);
  if (!input) {
    return false;
  }
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream fields(line);
    std::string key;
    fields >> key;
    if (key == "pressure") {
      fields >> design.pressure;
    } else if (key == "wander") {
      fields >> design.wander;
    } else if (key == "variation") {
      fields >> design.variation;
    } else if (key == "seed") {
      fields >> design.seed;
    } else if (key == "elasticity_scale") {
      fields >> design.elasticity_scale;
    } else if (key == "radiating") {
      int flag = 0;
      fields >> flag;
      design.radiating = flag != 0;
    } else if (key == "control_rate") {
      fields >> design.control_rate;
    } else if (key == "air_lowpass") {
      fields >> design.air_lowpass;
    } else if (key == "air_taps") {
      fields >> design.air_taps[0] >> design.air_taps[1] >> design.air_taps[2];
    } else if (key == "subglottal_taps") {
      fields >> design.subglottal_taps[0] >> design.subglottal_taps[1] >> design.subglottal_taps[2];
    } else if (key == "lung_resistance") {
      fields >> design.lung_resistance;
    } else if (key == "lung_corner") {
      fields >> design.lung_corner;
    } else if (key == "jet_loss") {
      int flag = 0;
      fields >> flag;
      design.jet_loss = flag != 0;
    } else if (key == "cord_contact") {
      fields >> design.cord_contact;
    } else if (key == "wall_loss") {
      fields >> design.wall_loss;
    } else if (key == "cord_loss") {
      fields >> design.cord_loss;
    } else if (key == "cord_taper") {
      fields >> design.cord_taper;
    } else if (key == "cord_length") {
      fields >> design.cord_length;
    } else if (key == "cord_strain") {
      fields >> design.cord_strain;
    } else if (key == "orifice") {
      int flag = 0;
      fields >> flag;
      design.orifice = flag != 0;
    } else if (key == "axial_air") {
      int flag = 0;
      fields >> flag;
      design.axial_air = flag != 0;
    } else if (key == "loss_per_length") {
      int flag = 0;
      fields >> flag;
      design.loss_per_length = flag != 0;
    } else if (key == "acoustic_spacing") {
      int flag = 0;
      fields >> flag;
      design.acoustic_spacing = flag != 0;
    } else if (key == "band") {
      SheetBand band;
      fields >> band.name >> band.length >> band.radius >> band.material >> band.mass >> band.elasticity >>
          band.offset >> band.anchor >> band.support >> band.points >> band.role;
      if (!fields) {
        return false;
      }
      design.bands.push_back(band);
    }
  }
  return !design.bands.empty();
}

SheetTube::SheetTube(const SheetDesign& design, double rate)
    : bands_(design.bands),
      rate_(rate),
      h_(1.0 / rate),
      variation_(static_cast<int>(design.bands.size()), design.seed),
      variation_depth_(design.variation),
      pressure_target_(design.pressure),
      wander_depth_(design.wander),
      loss_per_length_(design.loss_per_length),
      axial_air_(design.axial_air),
      orifice_(design.orifice),
      cord_strain_(design.cord_strain),
      fold_length_(design.orifice ? design.cord_length : 0.0),
      fold_taper_(std::max(design.cord_taper, 0.0)),
      fold_contact_(design.orifice && design.cord_length > 0.0 ? std::max(design.cord_contact, 0.0) : 0.0),
      jet_loss_(design.jet_loss),
      lung_resistance_(design.lung_resistance) {
  for (int k = 0; k < 3; ++k) {
    air_taps_[k] = design.air_taps[k];
    subglottal_taps_[k] = design.subglottal_taps[0] >= 0.0 ? design.subglottal_taps[k] : design.air_taps[k];
  }
  for (std::size_t b = 0; b < design.bands.size(); ++b) {
    if (design.bands[b].role == "cord") {
      first_cord_band_ = static_cast<int>(b);
      break;
    }
  }
  if (design.air_lowpass > 0.0) {
    air_lowpass_gain_ = 1.0 - std::exp(-2.0 * kPiValue * design.air_lowpass / rate);
  }
  if (lung_resistance_ > 0.0) {
    double inertance = lung_resistance_ / (2.0 * kPiValue * std::max(design.lung_corner, 1.0));
    double twice = 2.0 * inertance / h_;
    lung_gain_ = 1.0 / (twice + lung_resistance_);
    lung_keep_ = (twice - lung_resistance_) / (twice + lung_resistance_);
  }
  if (design.radiating) {
    outlet_ = End::radiating;
  }
  for (SheetBand& band : bands_) {
    band.elasticity *= design.elasticity_scale;
    materials_.push_back(sheet_material(band.material));
    if (band.role == "cord" && design.cord_loss >= 0.0) {
      materials_.back().loss = design.cord_loss;
    }
    if (band.role != "cord") {
      materials_.back().loss *= design.wall_loss;
    }
    if (design.acoustic_spacing && !(design.orifice && band.role == "cord")) {
      double spacing = kSound / rate;
      // Rounded down (at least one node): every node is at least one sample
      // of sound long, the remainder a fractional delay. (Rounded, with at
      // least two, short sections got nodes shorter than a sample; clamped
      // to one sample, the air crossed them at ~170 m/s.)
      band.points = std::max(1, static_cast<int>(std::floor(band.length / spacing)));
    }
  }
  for (const SheetBand& band : bands_) {
    starts_.push_back(n_);
    counts_.push_back(band.points);
    n_ += band.points;
  }
  std::size_t n = static_cast<std::size_t>(n_);
  in_cord_.assign(n, 0);
  injected_.assign(n, 0.0);
  series_source_.assign(n, 0.0);
  axial_link_.assign(n, -1);
  glottal_pressure_.assign(n, 0.0);
  for (std::size_t b = 0; b < bands_.size(); ++b) {
    if (bands_[b].role != "cord") {
      continue;
    }
    for (int j = 0; j < counts_[b]; ++j) {
      in_cord_[static_cast<std::size_t>(starts_[b] + j)] = 1;
    }
    if (cord_first_ < 0) {
      cord_first_ = starts_[b];
    }
    cord_last_ = starts_[b] + counts_[b] - 1;
    cord_length_ += bands_[b].length;
  }
  if (cord_first_ <= 0 || cord_last_ + 1 >= n_) {
    orifice_ = false;  // the orifice needs air on both sides
  }
  for (std::vector<double>* field :
       {&x_, &dx_, &r_, &area0_, &arc_, &mass_, &root_mass_, &surface_, &tension_, &hoop_, &foundation_, &loss_,
        &impedance_, &support_diagonal_, &base_diagonal_, &viscous_, &factor_from_, &factor_to_, &factor_now_,
        &viscous_now_, &impedance_now_, &diag_, &rhs_, &scratch_, &q_, &v_, &pressure_, &sum_, &admit_,
        &solid_incoming_, &solid_admit_, &radii_, &force_}) {
    field->assign(n, 0.0);
  }
  for (std::vector<double>* field : {&axial_, &support_off_, &axial_now_, &off_}) {
    field->assign(n > 0 ? n - 1 : 0, 0.0);
  }
  band_at_.assign(n, 0);
  band_factor_.assign(bands_.size(), 1.0);
  band_support_.assign(bands_.size(), 1.0);
  radius_goal_.assign(n, 0.0);
  radius_step_.assign(n, 0.0);
  fir_.assign(n, std::vector<double>(3, 0.0));
  factor_from_.assign(n, 1.0);
  factor_to_.assign(n, 1.0);
  factor_now_.assign(n, 1.0);
  control_period_ = std::max(1, static_cast<int>(std::lround(rate_ / std::max(design.control_rate, 1.0))));
  compile();
  build_links();
  configure_links();
  viscous_now_ = viscous_;
  impedance_now_ = impedance_;
  for (int node = 0; node < n_; ++node) {
    radii_[static_cast<std::size_t>(node)] = r_[static_cast<std::size_t>(node)];
  }
}

// After the supports moved: the sheet's pieces at the new rest radii (the
// node positions stay; masses, tensions, foundation, supports, damping and
// tissue impedance follow the geometry).
void SheetTube::recompile_pieces() {
  for (std::size_t i = 0; i < static_cast<std::size_t>(n_); ++i) {
    area0_[i] = kPiValue * r_[i] * r_[i];
  }
  compile_pieces();
  update_variable_material();
}

// Geometry and the nominal supported sheet (Astra's compile + compileSheet,
// the axisymmetric family).
void SheetTube::compile() {
  double start = 0.0;
  int node = 0;
  for (std::size_t b = 0; b < bands_.size(); ++b) {
    const SheetBand& band = bands_[b];
    for (int j = 0; j < band.points; ++j) {
      std::size_t i = static_cast<std::size_t>(node);
      dx_[i] = band.length / band.points;
      x_[i] = start + (j + 0.5) * band.length / band.points;
      band_at_[i] = static_cast<int>(b);
      r_[i] = rest_of(i);
      area0_[i] = kPiValue * r_[i] * r_[i];
      ++node;
    }
    start += band.length;
  }
  length_ = start;
  compile_pieces();
}

void SheetTube::compile_pieces() {
  std::size_t n = static_cast<std::size_t>(n_);
  for (std::size_t i = 0; i < n; ++i) {
    const SheetBand& band = bands_[static_cast<std::size_t>(band_at_[i])];
    const SheetMaterial& material = materials_[static_cast<std::size_t>(band_at_[i])];
    bool flat = slit(i);
    // A slit cord's sector carries an eighth of the two folds' medial width.
    double circ = flat ? 2.0 * fold_length_ / kSectors : r_[i] * kAngle;
    arc_[i] = circ;
    double slope = (!flat && i > 0 && i + 1 < n) ? (r_[i + 1] - r_[i - 1]) / (x_[i + 1] - x_[i - 1]) : 0.0;
    double area = circ * dx_[i] * std::hypot(1.0, slope);
    surface_[i] = circ * dx_[i];
    mass_[i] = material.density * material.thickness * area * band.mass;
    root_mass_[i] = std::sqrt(mass_[i]);
    double stretch = flat ? 0.0 : std::abs(std::log(r_[i] / band.radius)) * 0.02;
    tension_[i] = material.modulus * material.thickness * band.elasticity * (0.035 + stretch);
    hoop_[i] = tension_[i] * dx_[i] / circ;
    double curvature = flat ? 0.0 : tension_[i] * surface_[i] / (r_[i] * r_[i]);
    foundation_[i] = curvature + material.modulus * band.elasticity * surface_[i] / 0.015;
    loss_[i] = material.loss;
    for (int k = 0; k < 3; ++k) {
      fir_[i][static_cast<std::size_t>(k)] = material.fir[k];
    }
    impedance_[i] = std::sqrt(material.density * material.thickness * band.mass * tension_[i]) * circ * 0.08;
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    axial_[i] = 2.0 * tension_[i] * tension_[i + 1] / (tension_[i] + tension_[i + 1]) * (arc_[i] + arc_[i + 1]) /
                2.0 / std::hypot(x_[i + 1] - x_[i], r_[i + 1] - r_[i]);
  }
  std::fill(support_diagonal_.begin(), support_diagonal_.end(), 0.0);
  std::fill(support_off_.begin(), support_off_.end(), 0.0);
  // Each section's ring support: a footprint (hat weights) within it.
  for (std::size_t b = 0; b < bands_.size(); ++b) {
    const SheetBand& band = bands_[b];
    int count = counts_[b];
    std::size_t first = static_cast<std::size_t>(starts_[b]);
    double at = band.anchor * (count - 1);
    double k = band.support * band_support_[b] * (tension_[first] * arc_[first] / dx_[first] + foundation_[first]);
    for (int i = 0; i < count; ++i) {
      double wi = std::max(0.0, 1.0 - std::abs(i - at));
      support_diagonal_[first + static_cast<std::size_t>(i)] += k * wi * wi;
      if (i + 1 < count) {
        double wj = std::max(0.0, 1.0 - std::abs(i + 1 - at));
        support_off_[first + static_cast<std::size_t>(i)] += k * wi * wj;
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    double diagonal = foundation_[i] + support_diagonal_[i];
    if (i > 0) {
      diagonal += axial_[i - 1];
    }
    if (i + 1 < n) {
      diagonal += axial_[i];
    }
    base_diagonal_[i] = diagonal;
    viscous_[i] = 2.0 * loss_[i] * std::sqrt(mass_[i] * diagonal);
  }
  axial_now_ = axial_;
}

void SheetTube::build_links() {
  // Air: axial, diagonal (to the next sector: the same one in the
  // axisymmetric quotient) and around each ring; tissue: axial and around.
  // With `axial_air` the axisymmetric quotient keeps only the axial air path
  // (its weight the axial plus diagonal ones, 0.6): with every sector alike,
  // the ring crossing takes a fraction of a sample and changes nothing, but
  // as a one-sample self-loop at every node it added compliance and slowed
  // the air (an /a/ tube resonated at 257/2098 Hz, not ~700/1100).
  for (int i = 0; i + 1 < n_; ++i) {
    bool through_cord = orifice_ && i + 1 >= cord_first_ && i <= cord_last_;
    if (through_cord) {
      // The orifice carries the air across the cord (solve_orifice).
      Link tissue;
      tissue.a = i;
      tissue.b = i + 1;
      solid_.push_back(tissue);
      continue;
    }
    Link axial;
    axial.a = i;
    axial.b = i + 1;
    axial_link_[static_cast<std::size_t>(i)] = static_cast<int>(air_.size());
    axial.weight = axial_air_ ? 0.60 : 0.36;
    air_.push_back(axial);
    if (!axial_air_) {
      Link diagonal;
      diagonal.a = i;
      diagonal.b = i + 1;
      diagonal.weight = 0.24;
      diagonal.distance = -1.0;  // marks the angular (diagonal) path
      air_.push_back(diagonal);
    }
    Link tissue;
    tissue.a = i;
    tissue.b = i + 1;
    solid_.push_back(tissue);
  }
  for (int i = 0; i < n_; ++i) {
    if (!axial_air_) {
      Link ring;
      ring.a = i;
      ring.b = i;
      ring.self = true;
      ring.weight = 0.10;
      ring.distance = -1.0;
      air_.push_back(ring);
    }
    Link tissue;
    tissue.a = i;
    tissue.b = i;
    tissue.self = true;
    tissue.distance = -1.0;
    solid_.push_back(tissue);
  }
}

// Astra's taps are a loss per hop of her grid (about 4.4 mm per node). With
// `loss_per_length` the loss becomes per unit length: a hop of length d
// keeps (sum of taps)^(d / 4.4 mm), same shape. (Per hop, a 96 kHz tube with
// four times the hops lost four times as much and its oscillation died.)
void SheetTube::scale_taps(double* taps, double distance) const {
  if (!loss_per_length_) {
    return;
  }
  double total = taps[0] + taps[1] + taps[2];
  double kept = std::pow(total, distance / kReferenceHop);
  for (int k = 0; k < 3; ++k) {
    taps[k] *= kept / total;
  }
}

void SheetTube::set_filter(Link& link, double delay) {
  double clamped = std::clamp(delay, 1.0, 4090.0);
  int size = power_of_two_at_least(clamped + 5.0);
  if (size > static_cast<int>(link.ab.size())) {
    std::vector<double> ab(static_cast<std::size_t>(size), 0.0);
    std::vector<double> ba(static_cast<std::size_t>(size), 0.0);
    int old = static_cast<int>(link.ab.size());
    for (int lag = 1; lag <= old; ++lag) {
      ab[static_cast<std::size_t>((size - lag) & (size - 1))] = link.ab[static_cast<std::size_t>((link.pos - lag + old) & link.mask)];
      ba[static_cast<std::size_t>((size - lag) & (size - 1))] = link.ba[static_cast<std::size_t>((link.pos - lag + old) & link.mask)];
    }
    link.ab = ab;
    link.ba = ba;
    link.mask = size - 1;
    link.pos = 0;
  }
  link.delay = static_cast<int>(std::floor(clamped));
  double f = clamped - link.delay;
  // Fractional delay and three-tap attenuation fused into one FIR.
  link.h[0] = (1.0 - f) * link.taps[0];
  link.h[1] = f * link.taps[0] + (1.0 - f) * link.taps[1];
  link.h[2] = f * link.taps[1] + (1.0 - f) * link.taps[2];
  link.h[3] = f * link.taps[2];
}

void SheetTube::configure_links() {
  for (Link& link : air_) {
    if (link.ab.empty()) {
      link.ab.assign(16, 0.0);
      link.ba.assign(16, 0.0);
    }
    bool angular = link.distance < 0.0;
    std::size_t a = static_cast<std::size_t>(link.a);
    std::size_t b = static_cast<std::size_t>(link.b);
    double theta = angular ? kAngle : 0.0;
    double distance = std::hypot(x_[a] - x_[b],
                                 0.5 * std::sqrt(std::max(r_[a] * r_[a] + r_[b] * r_[b] - 2.0 * r_[a] * r_[b] * std::cos(theta), 0.0)));
    bool below = first_cord_band_ >= 0 && band_at_[a] < first_cord_band_ && band_at_[b] < first_cord_band_;
    for (int k = 0; k < 3; ++k) {
      link.taps[k] = below ? subglottal_taps_[k] : air_taps_[k];
    }
    scale_taps(link.taps, std::max(distance, 1e-6));
    set_filter(link, distance / kSound * rate_);
    link.gain = std::sqrt(link.weight * kPiValue / (kSectors * kRho * kSound));
  }
  std::fill(solid_admit_.begin(), solid_admit_.end(), 0.0);
  for (Link& link : solid_) {
    if (link.ab.empty()) {
      link.ab.assign(16, 0.0);
      link.ba.assign(16, 0.0);
    }
    bool angular = link.distance < 0.0;
    std::size_t a = static_cast<std::size_t>(link.a);
    std::size_t b = static_cast<std::size_t>(link.b);
    link.gain = std::sqrt(impedance_[a] * impedance_[b]);
    solid_admit_[a] += link.gain;
    solid_admit_[b] += link.gain;
    double distance = angular ? angular_distance(a) : std::hypot(x_[a] - x_[b], r_[a] - r_[b]);
    std::size_t ia = std::min(a, static_cast<std::size_t>(n_ - 2));
    std::size_t ib = std::min(b, static_cast<std::size_t>(n_ - 2));
    double speed = std::sqrt((axial_[ia] * dx_[a] * dx_[a] / mass_[a] + axial_[ib] * dx_[b] * dx_[b] / mass_[b]) / 2.0);
    double delay = distance / std::max(0.1, speed) * rate_;
    for (int k = 0; k < 3; ++k) {
      link.taps[k] = std::sqrt(fir_[a][static_cast<std::size_t>(k)] * fir_[b][static_cast<std::size_t>(k)]);
    }
    scale_taps(link.taps, std::max(distance, 1e-6));
    // Reserve the slowest-delay history (elasticity at -10%) first.
    // (Muscle control may slacken a section to a fifth of its elasticity
    // and supports move: reserve four times the nominal delay.)
    set_filter(link, 4.0 * delay);
    set_filter(link, delay);
    if (angular) {
      link.distance = -1.0;
    }
  }
}

void SheetTube::set_strain_relaxation(double fraction, double seconds, double log_rate_sd, std::uint32_t seed) {
  std::size_t n = static_cast<std::size_t>(n_);
  relax_fraction_ = std::clamp(fraction, 0.0, 0.9);
  relax_seconds_ = std::max(1e-3, seconds);
  relax_sd_ = std::max(0.0, log_rate_sd);
  relax_period_ = std::max(1, static_cast<int>(std::lround(rate_ / 1000.0)));
  relax_phase_ = 0;
  relax_random_ = seed != 0U ? seed : 1U;
  strain_memory_.assign(n, 0.0);
  relax_log_rate_.assign(n, 0.0);
  relax_gain_.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    strain_memory_[i] = q_[i];
    relax_log_rate_[i] = relax_sd_ * relax_gauss();
  }
  advance_relaxation();
}

double SheetTube::relax_gauss() {
  relax_random_ ^= relax_random_ << 13U;
  relax_random_ ^= relax_random_ >> 17U;
  relax_random_ ^= relax_random_ << 5U;
  double u = (static_cast<double>(relax_random_) + 1.0) / 4294967297.0;
  relax_random_ ^= relax_random_ << 13U;
  relax_random_ ^= relax_random_ >> 17U;
  relax_random_ ^= relax_random_ << 5U;
  double w = static_cast<double>(relax_random_) / 4294967296.0;
  return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * kPiValue * w);
}

// Once per millisecond: the log rates take an innovation (an OU process of
// correlation relax_seconds_), and each node's per-sample gain follows.
void SheetTube::advance_relaxation() {
  double retain = std::exp(-0.001 / relax_seconds_);
  double innovation = relax_sd_ * std::sqrt(1.0 - retain * retain);
  for (std::size_t i = 0; i < relax_gain_.size(); ++i) {
    relax_log_rate_[i] = retain * relax_log_rate_[i] + innovation * relax_gauss();
    double rate = std::exp(relax_log_rate_[i] - 0.5 * relax_sd_ * relax_sd_) / relax_seconds_;
    relax_gain_[i] = -std::expm1(-h_ * rate);
  }
}

void SheetTube::settle_geometry() {
  std::size_t n = static_cast<std::size_t>(n_);
  for (std::size_t i = 0; i < n; ++i) {
    r_[i] = rest_of(i);
    radii_[i] = soft_contact(i) ? r_[i] + q_[i] : std::max(kFloor, r_[i] + q_[i]);
    double factor = band_factor_[static_cast<std::size_t>(band_at_[i])];
    factor_from_[i] = factor;
    factor_to_[i] = factor;
    factor_now_[i] = factor;
  }
  geometry_dirty_ = false;
  gliding_ = 0;
  recompile_pieces();
}

void SheetTube::set_band_offset(int band, double offset) {
  bands_[static_cast<std::size_t>(band)].offset = offset;
  geometry_dirty_ = true;
}

int SheetTube::add_port(int node) {
  SheetPort port;
  port.index = node;
  ports_.push_back(port);
  return static_cast<int>(ports_.size()) - 1;
}

double SheetTube::cord_radius() const {
  double narrowest = 1.0;
  for (int node = 0; node < n_; ++node) {
    if (bands_[static_cast<std::size_t>(band_at_[static_cast<std::size_t>(node)])].role == "cord") {
      narrowest = std::min(narrowest, radii_[static_cast<std::size_t>(node)]);
    }
  }
  return narrowest;
}

// Material variation (Astra's advanceVariation): at 200 Hz the sections'
// elasticity factors move; within each period the sheet's factors are
// interpolated per sample, and the damping, tissue-path impedances and
// delays follow at the update.
void SheetTube::advance_control() {
  if (control_phase_ == 0) {
    double dt = static_cast<double>(control_period_) / rate_;
    if (variation_depth_ > 0.0) {
      variation_.advance(dt);
    }
    for (std::size_t i = 0; i < static_cast<std::size_t>(n_); ++i) {
      factor_from_[i] = factor_to_[i];
      double varied = variation_depth_ > 0.0 ? 1.0 + variation_depth_ * variation_.value(band_at_[i]) : 1.0;
      factor_to_[i] = varied * band_factor_[static_cast<std::size_t>(band_at_[i])];
    }
    if (geometry_dirty_) {
      // The supports' new rest radii: glide there over this period.
      for (std::size_t i = 0; i < static_cast<std::size_t>(n_); ++i) {
        radius_goal_[i] = rest_of(i);
        radius_step_[i] = (radius_goal_[i] - r_[i]) / control_period_;
      }
      geometry_dirty_ = false;
      gliding_ = control_period_;
    }
    update_variable_material();
  }
  if (gliding_ > 0) {
    for (std::size_t i = 0; i < static_cast<std::size_t>(n_); ++i) {
      r_[i] += radius_step_[i];
    }
    if (--gliding_ == 0) {
      recompile_pieces();
    }
  }
  double mix = static_cast<double>(control_phase_ + 1) / control_period_;
  for (std::size_t i = 0; i < static_cast<std::size_t>(n_); ++i) {
    factor_now_[i] = factor_from_[i] + mix * (factor_to_[i] - factor_from_[i]);
  }
  control_phase_ = (control_phase_ + 1) % control_period_;
}

void SheetTube::update_variable_material() {
  std::size_t n = static_cast<std::size_t>(n_);
  std::vector<double>& f = factor_to_;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double a = tension_[i];
    double b = tension_[i + 1];
    axial_now_[i] = axial_[i] * f[i] * f[i + 1] * (a + b) / (a * f[i] + b * f[i + 1]);
  }
  for (std::size_t i = 0; i < n; ++i) {
    double k = (foundation_[i] + support_diagonal_[i]) * f[i] + (i > 0 ? axial_now_[i - 1] : 0.0) +
               (i + 1 < n ? axial_now_[i] : 0.0);
    viscous_now_[i] = viscous_[i] * std::sqrt(k / base_diagonal_[i]);
    impedance_now_[i] = impedance_[i] * std::sqrt(f[i]);
  }
  std::fill(solid_admit_.begin(), solid_admit_.end(), 0.0);
  for (Link& link : solid_) {
    std::size_t a = static_cast<std::size_t>(link.a);
    std::size_t b = static_cast<std::size_t>(link.b);
    link.gain = std::sqrt(impedance_now_[a] * impedance_now_[b]);
    solid_admit_[a] += link.gain;
    solid_admit_[b] += link.gain;
    bool angular = link.distance < 0.0;
    double distance = angular ? angular_distance(a) : std::hypot(x_[a] - x_[b], r_[a] - r_[b]);
    std::size_t ia = std::min(a, n - 2);
    std::size_t ib = std::min(b, n - 2);
    double speed = std::sqrt((axial_now_[ia] * dx_[a] * dx_[a] / mass_[a] + axial_now_[ib] * dx_[b] * dx_[b] / mass_[b]) / 2.0);
    double delay = std::clamp(distance / std::max(0.1, speed) * rate_, 1.0, 4090.0);
    int whole = static_cast<int>(std::floor(delay));
    if (whole + 3 > link.mask) {
      std::fprintf(stderr, "tissue delay %d over its history %d (link %d-%d)\n", whole, link.mask, link.a, link.b);
      fault_ = true;
      return;
    }
    double frac = delay - whole;
    link.delay = whole;
    link.h[0] = (1.0 - frac) * link.taps[0];
    link.h[1] = frac * link.taps[0] + (1.0 - frac) * link.taps[1];
    link.h[2] = frac * link.taps[1] + (1.0 - frac) * link.taps[2];
    link.h[3] = frac * link.taps[2];
  }
}

void SheetTube::receive_air() {
  for (Link& link : air_) {
    std::size_t a = static_cast<std::size_t>(link.a);
    std::size_t b = static_cast<std::size_t>(link.b);
    double s = link.gain * std::min(radii_[a], radii_[b]);
    double y = s * s;
    int m = link.mask;
    int p = (link.pos - link.delay) & m;
    std::size_t p0 = static_cast<std::size_t>(p);
    std::size_t p1 = static_cast<std::size_t>((p - 1) & m);
    std::size_t p2 = static_cast<std::size_t>((p - 2) & m);
    std::size_t p3 = static_cast<std::size_t>((p - 3) & m);
    double va = (link.h[0] * link.ba[p0] + link.h[1] * link.ba[p1] + link.h[2] * link.ba[p2] + link.h[3] * link.ba[p3]) / s;
    double vb = link.self ? va
                          : (link.h[0] * link.ab[p0] + link.h[1] * link.ab[p1] + link.h[2] * link.ab[p2] +
                             link.h[3] * link.ab[p3]) / s;
    if (air_lowpass_gain_ > 0.0) {
      link.low_a += air_lowpass_gain_ * (va - link.low_a);
      link.low_b += air_lowpass_gain_ * (vb - link.low_b);
      va = link.low_a;
      vb = link.self ? va : link.low_b;
    }
    link.scale = s;
    link.ia = va;
    link.ib = vb;
    if (jet_loss_ && !link.self) {
      // Series resistance and source at the ends: the line's Thevenin end
      // (2 v, 1/y) behind R (and e) is a Norton source y'(2 v + e), y' = y/(1 + R y).
      double ya = y / (1.0 + link.ra * y);
      double yb = y / (1.0 + link.rb * y);
      double eb = link.b == link.a + 1 ? series_source_[b] : 0.0;
      link.ya = ya;
      link.yb = yb;
      sum_[a] += ya * va;
      sum_[b] += yb * (vb + 0.5 * eb);
      admit_[a] += ya;
      admit_[b] += yb;
      continue;
    }
    sum_[a] += y * va;
    sum_[b] += y * vb;
    admit_[a] += y;
    admit_[b] += y;
  }
}

void SheetTube::receive_solid() {
  for (Link& link : solid_) {
    std::size_t a = static_cast<std::size_t>(link.a);
    std::size_t b = static_cast<std::size_t>(link.b);
    int m = link.mask;
    int p = (link.pos - link.delay) & m;
    std::size_t p0 = static_cast<std::size_t>(p);
    std::size_t p1 = static_cast<std::size_t>((p - 1) & m);
    std::size_t p2 = static_cast<std::size_t>((p - 2) & m);
    std::size_t p3 = static_cast<std::size_t>((p - 3) & m);
    double va = link.h[0] * link.ba[p0] + link.h[1] * link.ba[p1] + link.h[2] * link.ba[p2] + link.h[3] * link.ba[p3];
    double vb = link.self ? va
                          : link.h[0] * link.ab[p0] + link.h[1] * link.ab[p1] + link.h[2] * link.ab[p2] +
                                link.h[3] * link.ab[p3];
    link.ia = va;
    link.ib = vb;
    solid_incoming_[a] += link.gain * va;
    solid_incoming_[b] += link.gain * vb;
  }
}

void SheetTube::send(std::vector<Link>& links, const std::vector<double>& field, bool air) {
  for (Link& link : links) {
    double s = air ? link.scale : 1.0;
    std::size_t at = static_cast<std::size_t>(link.pos);
    link.ab[at] = (field[static_cast<std::size_t>(link.a)] - link.ia) * s;
    link.ba[at] = (field[static_cast<std::size_t>(link.b)] - link.ib) * s;
    link.pos = (link.pos + 1) & link.mask;
  }
}

// Air waves leave through the series elements: the line sees the port's
// pressure, the node's plus the resistance's drop less the source.
void SheetTube::send_air() {
  for (Link& link : air_) {
    std::size_t a = static_cast<std::size_t>(link.a);
    std::size_t b = static_cast<std::size_t>(link.b);
    std::size_t at = static_cast<std::size_t>(link.pos);
    if (link.self) {
      link.ab[at] = (pressure_[a] - link.ia) * link.scale;
      link.ba[at] = (pressure_[b] - link.ib) * link.scale;
      link.pos = (link.pos + 1) & link.mask;
      continue;
    }
    double eb = link.b == link.a + 1 ? series_source_[b] : 0.0;
    double into_a = link.ya * (2.0 * link.ia - pressure_[a]);
    double into_b = link.yb * (2.0 * link.ib + eb - pressure_[b]);
    double port_a = pressure_[a] + link.ra * into_a;
    double port_b = pressure_[b] + link.rb * into_b - eb;
    link.ab[at] = (port_a - link.ia) * link.scale;
    link.ba[at] = (port_b - link.ib) * link.scale;
    link.pos = (link.pos + 1) & link.mask;
    link.flow = kSectors * 0.5 * (into_b - into_a);
    update_jet_loss(link);
  }
  std::fill(series_source_.begin(), series_source_.end(), 0.0);
}

// The next step's series resistances from this step's flow: Borda-Carnot at
// the end where the flow expands, Poiseuille shared by both ends. Per sector
// a resistance is eight times the whole ring's.
void SheetTube::update_jet_loss(Link& link) {
  constexpr double kViscosity = 1.8e-5;  // Pa s
  std::size_t a = static_cast<std::size_t>(link.a);
  std::size_t b = static_cast<std::size_t>(link.b);
  double area_a = kPiValue * radii_[a] * radii_[a];
  double area_b = kPiValue * radii_[b] * radii_[b];
  double narrow = std::min(area_a, area_b);
  double length = std::abs(x_[b] - x_[a]);
  double viscous = 8.0 * kPiValue * kViscosity * length / (narrow * narrow);
  double flow = std::abs(link.flow);
  double expand_b = link.flow > 0.0 && area_b > area_a ? 1.0 / area_a - 1.0 / area_b : 0.0;
  double expand_a = link.flow < 0.0 && area_a > area_b ? 1.0 / area_b - 1.0 / area_a : 0.0;
  link.ra = kSectors * (0.5 * viscous + 0.5 * kRho * flow * expand_a * expand_a);
  link.rb = kSectors * (0.5 * viscous + 0.5 * kRho * flow * expand_b * expand_b);
}

double SheetTube::axial_flow(int node) const {
  int index = axial_link_[static_cast<std::size_t>(node)];
  return index >= 0 ? air_[static_cast<std::size_t>(index)].flow : 0.0;
}

// The supported sheet's implicit midpoint: (M + h^2 K / 4) y = M (2 q + h v).
void SheetTube::advance_sheet() {
  std::size_t n = static_cast<std::size_t>(n_);
  double quarter = 0.25 * h_ * h_;
  const std::vector<double>& f = factor_now_;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double a = tension_[i];
    double b = tension_[i + 1];
    axial_now_[i] = axial_[i] * f[i] * f[i + 1] * (a + b) / (a * f[i] + b * f[i + 1]);
  }
  for (std::size_t i = 0; i < n; ++i) {
    double k = (foundation_[i] + support_diagonal_[i]) * f[i] + (i > 0 ? axial_now_[i - 1] : 0.0) +
               (i + 1 < n ? axial_now_[i] : 0.0);
    diag_[i] = mass_[i] + quarter * k;
    if (i + 1 < n) {
      off_[i] = quarter * (support_off_[i] * f[i] - axial_now_[i]);
    }
    rhs_[i] = mass_[i] * (2.0 * q_[i] + h_ * v_[i]);
    if (relax_fraction_ > 0.0) {
      // The relaxed share of the foundation's stress, held at its memory:
      // a force through the midpoint step (y = q + q').
      rhs_[i] += 0.5 * h_ * h_ * relax_fraction_ * foundation_[i] * f[i] * strain_memory_[i];
    }
  }
  scratch_[0] = n > 1 ? off_[0] / diag_[0] : 0.0;
  rhs_[0] /= diag_[0];
  for (std::size_t i = 1; i < n; ++i) {
    double denominator = diag_[i] - off_[i - 1] * scratch_[i - 1];
    scratch_[i] = i + 1 < n ? off_[i] / denominator : 0.0;
    rhs_[i] = (rhs_[i] - off_[i - 1] * rhs_[i - 1]) / denominator;
  }
  for (std::size_t i = n - 1; i > 0; --i) {
    rhs_[i - 1] -= scratch_[i - 1] * rhs_[i];
  }
  for (std::size_t i = 0; i < n; ++i) {
    double old = q_[i];
    double response = rhs_[i];
    q_[i] = response - old;
    v_[i] = (2.0 * response - 4.0 * old) / h_ - v_[i];
  }
}

// The cord's air: a quasi-steady orifice between the junctions below (a) and
// above (b). With the narrowest ring's area A and the air's velocity there v,
//   rho L dv/dt = (Pa - Pb) - A v (Za + Zb) - rho v|v| / 2 - 8 pi mu L v / A,
// Pa, Pb the junctions' pressures without the orifice and Za, Zb their
// impedances (whole ring). Backward Euler in v: one quadratic, closed form.
// As A closes, the viscous term stops the flow (no inertial flow through a
// shut cord). Upstream of the narrowest ring the surfaces feel the Bernoulli
// pressure; past it the jet has separated and they feel the pressure above.
void SheetTube::solve_orifice() {
  constexpr double kViscosity = 1.8e-5;  // Pa s, air
  std::size_t a = static_cast<std::size_t>(cord_first_ - 1);
  std::size_t b = static_cast<std::size_t>(cord_last_ + 1);
  double za = 1.0 / (kSectors * admit_[a]);
  double zb = 1.0 / (kSectors * admit_[b]);
  double pa = 2.0 * sum_[a] / admit_[a];
  double pb = 2.0 * sum_[b] / admit_[b];
  std::size_t narrowest = static_cast<std::size_t>(cord_first_);
  for (std::size_t i = narrowest; i <= static_cast<std::size_t>(cord_last_); ++i) {
    if (radii_[i] < radii_[narrowest]) {
      narrowest = i;
    }
  }
  double area = std::max(cord_area(radii_[narrowest]), 1e-12);
  // Poiseuille: a ring 8 pi mu L / A^2, a slit 12 mu L W^2 / A^3 (W the
  // folds' length), here per unit of velocity.
  double viscous = fold_length_ > 0.0 ? 12.0 * kViscosity * cord_length_ * fold_length_ * fold_length_ / (area * area)
                                      : 8.0 * kPiValue * kViscosity * cord_length_ / area;
  double inertia = kRho * cord_length_ / h_;
  double linear = inertia + area * (za + zb) + viscous;
  double push = pa - pb + glottal_noise_ + inertia * orifice_velocity_;
  glottal_noise_ = 0.0;
  double kinetic = 0.5 * kRho;
  double magnitude = 2.0 * std::abs(push) / (linear + std::sqrt(linear * linear + 4.0 * kinetic * std::abs(push)));
  double velocity = push >= 0.0 ? magnitude : -magnitude;
  double flow = velocity * area;
  orifice_velocity_ = velocity;
  orifice_area_ = area;
  sum_[a] -= 0.5 * flow / kSectors;
  sum_[b] += 0.5 * flow / kSectors;
  double below = pa - flow * za;
  double above = pb + flow * zb;
  double area_below = kPiValue * radii_[a] * radii_[a];
  for (std::size_t i = static_cast<std::size_t>(cord_first_); i <= static_cast<std::size_t>(cord_last_); ++i) {
    if (i <= narrowest) {
      double here = std::max(cord_area(radii_[i]), 1e-12);
      double bernoulli = 0.5 * kRho * flow * flow * (1.0 / (here * here) - 1.0 / (area_below * area_below));
      glottal_pressure_[i] = std::max(above, below - bernoulli);
    } else {
      glottal_pressure_[i] = above;
    }
  }
}

// The tissue path around a ring, to the next sector. On a slit cord all
// sectors move alike in this quotient and the path carries nothing (a delayed
// self-loop there echoed: the slit's 3.5 mm arc at tissue speed broke its
// voicing): the shortest, one-sample path.
double SheetTube::angular_distance(std::size_t node) const {
  return slit(node) ? 0.0 : arc_[node];
}

// A section's rest radius: its ring support's offset from its radius, held
// off the axis; a softly contacting cord may rest pressed (below contact).
double SheetTube::rest_of(std::size_t node) const {
  const SheetBand& band = bands_[static_cast<std::size_t>(band_at_[node])];
  if (soft_contact(node)) {
    return band.radius + band.offset;
  }
  return std::max(kFloor, band.radius + band.offset);
}

double SheetTube::cord_area(double radius) const {
  if (fold_length_ > 0.0) {
    // The gap along the folds is (radius - contact) + taper z / L, z from
    // the front; the open area is twice its positive part, integrated.
    double front = radius - kFloor;
    if (fold_taper_ <= 0.0) {
      return 2.0 * fold_length_ * std::max(front, 0.0);
    }
    double back = front + fold_taper_;
    if (back <= 0.0) {
      return 0.0;
    }
    if (front >= 0.0) {
      return 2.0 * fold_length_ * (front + 0.5 * fold_taper_);
    }
    return fold_length_ * back * back / fold_taper_;
  }
  return kPiValue * radius * radius;
}

// The attached glottis between the junctions below (a) and above (b): each
// junction, seen from the glottis, is a source of open-circuit pressure P
// through an impedance Z, i.e. an incident wave P/2 on a line of impedance Z
// (CGS for the glottis: 10 dyn/cm^2 per Pa, Z in dyn s/cm^5 = 1e-5 Pa s/m^3).
void SheetTube::solve_glottis() {
  std::size_t a = static_cast<std::size_t>(cord_first_ - 1);
  std::size_t b = static_cast<std::size_t>(cord_last_ + 1);
  double za = 1.0 / (kSectors * admit_[a]);
  double zb = 1.0 / (kSectors * admit_[b]);
  double pa = 2.0 * sum_[a] / admit_[a] + glottal_noise_;
  double pb = 2.0 * sum_[b] / admit_[b];
  glottal_noise_ = 0.0;
  double wave_below = 0.5 * pa * 10.0;
  double wave_above = 0.5 * pb * 10.0;
  double area_below = kRhoC / (za * 1e-5);
  double area_above = kRhoC / (zb * 1e-5);
  double out_below = 0.0;
  double out_above = 0.0;
  int substeps = std::max(1, static_cast<int>(std::lround(kStepRate / rate_)));
  double flow_sum = 0.0;
  for (int k = 0; k < substeps; ++k) {
    glottis_->step(wave_below, area_below, wave_above, area_above, out_below, out_above);
    flow_sum += glottis_->flow();
  }
  double flow = 1e-6 * flow_sum / substeps;  // m^3/s
  orifice_area_ = 1e-4 * glottis_->area();
  orifice_velocity_ = flow / std::max(orifice_area_, 1e-10);
  sum_[a] -= 0.5 * flow / kSectors;
  sum_[b] += 0.5 * flow / kSectors;
  // The sheet's cord section is tissue only now: it feels the mean of the
  // pressures around it.
  double mean = 0.5 * ((pa - flow * za) + (pb + flow * zb));
  for (std::size_t i = static_cast<std::size_t>(cord_first_); i <= static_cast<std::size_t>(cord_last_); ++i) {
    glottal_pressure_[i] = mean;
  }
}

bool SheetTube::step() {
  if (fault_) {
    return false;
  }
  std::size_t n = static_cast<std::size_t>(n_);
  ++tick_;
  advance_control();
  if (relax_fraction_ > 0.0) {
    for (std::size_t i = 0; i < n; ++i) {
      strain_memory_[i] += relax_gain_[i] * (q_[i] - strain_memory_[i]);
    }
    if (++relax_phase_ >= relax_period_) {
      relax_phase_ = 0;
      advance_relaxation();
    }
  }
  if (tick_ % 96 == 0) {
    random_ ^= random_ << 13U;
    random_ ^= random_ >> 17U;
    random_ ^= random_ << 5U;
    wander_ += 0.018 * (static_cast<double>(random_) / 2147483648.0 - 1.0) - 0.008 * wander_;
  }
  drive_ += (pressure_target_ * (1.0 + wander_depth_ * wander_) - drive_) * h_ / 0.035;
  std::fill(sum_.begin(), sum_.end(), 0.0);
  std::fill(admit_.begin(), admit_.end(), 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    radii_[i] = soft_contact(i) ? r_[i] + q_[i] : std::max(kFloor, r_[i] + q_[i]);
  }
  receive_air();
  double port_y = 1.0 / (kSectors * kRho * kSound);
  if (inlet_ == End::reservoir) {
    double yi = area0_[0] * port_y;
    admit_[0] += yi;
    sum_[0] += yi * drive_ / 2.0;
    if (lung_resistance_ > 0.0) {
      // The airway branch (R + L, trapezoid) as a Norton source, per sector.
      double carried = lung_keep_ * lung_flow_ + lung_gain_ * lung_drop_;
      admit_[0] += lung_gain_ / kSectors;
      sum_[0] += 0.5 * (carried + lung_gain_ * drive_) / kSectors;
    }
  }
  if (outlet_ == End::open) {
    admit_[n - 1] += kPiValue * radii_[n - 1] * radii_[n - 1] * port_y;
  }
  double radiation_resistance = 0.0;
  double radiation_gain = 0.0;
  if (outlet_ == End::radiating) {
    // Per sector (an eighth of the ring): R and L eight times the whole's.
    double area = kPiValue * radii_[n - 1] * radii_[n - 1];
    radiation_resistance = kSectors * 128.0 * kRho * kSound / (9.0 * kPiValue * kPiValue * area);
    double inertance = kSectors * 8.0 * kRho / (3.0 * kPiValue * std::sqrt(kPiValue * area));
    radiation_gain = 0.5 * h_ / inertance;
    double load_admittance = 1.0 / radiation_resistance + radiation_gain;
    double history = radiation_inertance_flow_ + radiation_gain * radiation_previous_pressure_;
    if (jet_loss_) {
      // The load (outflow = Y p_r + history) behind the exit's resistance:
      // outflow = y (p + history Z), y = 1 / (R + Z), Z = 1 / Y.
      exit_load_ = 1.0 / load_admittance;
      exit_history_ = history;
      exit_admittance_ = 1.0 / (exit_resistance_ + exit_load_);
      admit_[n - 1] += exit_admittance_;
      sum_[n - 1] -= 0.5 * exit_admittance_ * history * exit_load_;
    } else {
      admit_[n - 1] += load_admittance;
      sum_[n - 1] -= 0.5 * history;
    }
  }
  std::fill(solid_incoming_.begin(), solid_incoming_.end(), 0.0);
  receive_solid();
  for (SheetPort& port : ports_) {
    std::size_t i = static_cast<std::size_t>(port.index);
    admit_[i] += port.air_y;
    sum_[i] += port.air_y * port.air_in;
    solid_admit_[i] += port.solid_z;
    solid_incoming_[i] += port.solid_z * port.solid_in;
  }
  if (orifice_) {
    if (glottis_ != nullptr) {
      solve_glottis();
    } else {
      solve_orifice();
    }
  }
  if (any_injected_) {
    // Flow sources into the junctions: per sector, an eighth of each.
    for (std::size_t i = 0; i < n; ++i) {
      sum_[i] += 0.5 * injected_[i] / kSectors;
      injected_[i] = 0.0;
    }
    any_injected_ = false;
  }
  advance_sheet();
  // Collocated pressure/volume and tissue force/velocity exchange: one
  // scalar port per node; the sheet's modes already propagated elastically.
  for (std::size_t i = 0; i < n; ++i) {
    double surface = 2.0 * kPiValue * radii_[i] * dx_[i] / kSectors;
    double mass = mass_[i];
    double v = v_[i];
    if (orifice_ && in_cord_[i] != 0) {
      // The cord's surface feels the intraglottal pressure; its own motion
      // does not load the air (the orifice is quasi-steady). A slit's
      // pressed surface is the folds' fixed medial face.
      if (slit(i)) {
        surface = surface_[i];
      }
      double resistance = solid_admit_[i] + viscous_now_[i];
      double drive = surface * glottal_pressure_[i] + 2.0 * solid_incoming_[i];
      double contact_resistance = 0.0;
      if (soft_contact(i)) {
        // Pressed past contact: a spring (cubic beyond 50 um) and damping.
        double pressed = kFloor - (r_[i] + q_[i]);
        if (pressed > 0.0) {
          double stiffness = fold_contact_ * (foundation_[i] + support_diagonal_[i]) * factor_now_[i];
          double depth = pressed / 50e-6;
          drive += stiffness * pressed * (1.0 + depth * depth);
          contact_resistance = 2.0 * 0.5 * std::sqrt(mass * stiffness);
        }
      }
      if (q_[i] > 0.0 && cord_strain_ > 0.0) {
        // Tissue stiffens with strain (collagen): opening outward, the
        // restoring force grows as q^3 (doubling the stiffness at
        // `cord_strain`).
        // Without it subglottal pressure blew a widening ring open like a
        // balloon: its pressed surface grows while its hoop stiffness falls.
        double strain = q_[i] / cord_strain_;
        drive -= foundation_[i] * factor_now_[i] * q_[i] * strain * strain;
      }
      resistance += contact_resistance;
      double next = ((mass - h_ * resistance / 2.0) * v + h_ * drive) / (mass + h_ * resistance / 2.0);
      pressure_[i] = glottal_pressure_[i];
      v_[i] = next;
      force_[i] = 0.5 * (v + next);
      if (!std::isfinite(next) || std::abs(next) > 60.0 || std::abs(q_[i]) > 0.03) {
        std::fprintf(stderr, "cord node %zu: v %.3g q %.3g r %.3g p %.4g\n", i, next, q_[i], r_[i], glottal_pressure_[i]);
        fault_ = true;
        return false;
      }
      if (r_[i] + q_[i] < kFloor) {
        ++contacts_;
        if (!soft_contact(i)) {
          q_[i] = kFloor - r_[i];
          double bounded = std::max(0.0, v_[i]);
          v_[i] = bounded;
          force_[i] = 0.5 * (v + bounded);
        }
      }
      continue;
    }
    double admit = admit_[i];
    double resistance = surface * surface / admit + solid_admit_[i] + viscous_now_[i];
    double drive = surface * 2.0 * sum_[i] / admit + 2.0 * solid_incoming_[i];
    double next = ((mass - h_ * resistance / 2.0) * v + h_ * drive) / (mass + h_ * resistance / 2.0);
    double mid = 0.5 * (v + next);
    pressure_[i] = (2.0 * sum_[i] - surface * mid) / admit;
    v_[i] = next;
    force_[i] = mid;
    if (!std::isfinite(next) || std::abs(next) > 60.0 || std::abs(q_[i]) > 0.03) {
      std::fprintf(stderr, "node %zu (band %d): v %.3g q %.3g r %.3g p %.4g\n", i, band_at_[i], next, q_[i], r_[i], pressure_[i]);
      fault_ = true;
      return false;
    }
    if (r_[i] + q_[i] < kFloor) {
      q_[i] = kFloor - r_[i];
      double bounded = std::max(0.0, v_[i]);
      v_[i] = bounded;
      force_[i] = 0.5 * (v + bounded);
      pressure_[i] = (2.0 * sum_[i] - surface * force_[i]) / admit;
      ++contacts_;
    }
  }
  if (inlet_ == End::reservoir && lung_resistance_ > 0.0) {
    double drop = drive_ - pressure_[0];
    lung_flow_ = lung_keep_ * lung_flow_ + lung_gain_ * (lung_drop_ + drop);
    lung_drop_ = drop;
  }
  if (outlet_ == End::radiating) {
    double p = pressure_[n - 1];
    if (jet_loss_) {
      // The pressure at the load, past the exit's loss; the next step's
      // resistance from this step's flow.
      double outflow = exit_admittance_ * (p + exit_history_ * exit_load_);
      p = p - exit_resistance_ * outflow;
      double area = kPiValue * radii_[n - 1] * radii_[n - 1];
      double whole = kSectors * std::abs(outflow);
      exit_resistance_ = kSectors * 0.5 * kRho * whole / (area * area);
    }
    radiation_inertance_flow_ += radiation_gain * (p + radiation_previous_pressure_);
    radiation_previous_pressure_ = p;
    radiated_flow_ = kSectors * (p / radiation_resistance + radiation_inertance_flow_);
  }
  if (jet_loss_) {
    send_air();
  } else {
    send(air_, pressure_, true);
  }
  send(solid_, force_, false);
  for (SheetPort& port : ports_) {
    std::size_t i = static_cast<std::size_t>(port.index);
    port.air_out = pressure_[i] - port.air_in;
    port.solid_out = force_[i] - port.solid_in;
    solid_admit_[i] -= port.solid_z;
  }
  return true;
}

// ---------------------------------------------------------------------------

bool SheetSystem::load(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    return false;
  }
  std::string directory = path.substr(0, path.find_last_of('/') + 1);
  double pressure = 0.0;
  double wander = 0.0;
  double variation = 0.0;
  double scale = 1.0;
  std::string source;
  std::string output;
  struct TubeLine {
    std::string id;
    std::string file;
    std::string end;
    std::uint32_t seed = 73021U;
  };
  std::vector<TubeLine> tube_lines;
  std::vector<std::pair<std::string, std::string>> connections;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream fields(line);
    std::string key;
    fields >> key;
    if (key == "pressure") {
      fields >> pressure;
    } else if (key == "wander") {
      fields >> wander;
    } else if (key == "variation") {
      fields >> variation;
    } else if (key == "elasticity_scale") {
      fields >> scale;
    } else if (key == "source") {
      fields >> source;
    } else if (key == "output") {
      fields >> output;
    } else if (key == "tube") {
      TubeLine tube;
      fields >> tube.id >> tube.file >> tube.end >> tube.seed;
      tube_lines.push_back(tube);
    } else if (key == "connect") {
      std::string a;
      std::string b;
      fields >> a >> b;
      connections.emplace_back(a, b);
    }
  }
  for (const TubeLine& line_item : tube_lines) {
    SheetDesign design;
    if (!read_sheet_design(directory + line_item.file, design)) {
      return false;
    }
    bool is_source = line_item.id == source;
    design.pressure = is_source ? pressure : 0.0;
    design.wander = is_source ? wander : 0.0;
    design.variation = variation;
    design.elasticity_scale = scale;
    design.seed = line_item.seed;
    ids_.push_back(line_item.id);
    tubes_.emplace_back(design, rate_);
    SheetTube& tube = tubes_.back();
    tube.set_inlet(is_source ? SheetTube::End::reservoir : SheetTube::End::closed);
    tube.set_outlet(line_item.end == "open" ? SheetTube::End::open : SheetTube::End::closed);
  }
  source_ = tube_index(source);
  output_ = tube_index(output);
  for (const std::pair<std::string, std::string>& connection : connections) {
    Coupling coupling;
    int* tube_fields[2] = {&coupling.tube_a, &coupling.tube_b};
    int* port_fields[2] = {&coupling.port_a, &coupling.port_b};
    const std::string* ends[2] = {&connection.first, &connection.second};
    for (int side = 0; side < 2; ++side) {
      std::string text = *ends[side];
      std::vector<std::string> parts;
      std::stringstream split(text);
      std::string part;
      while (std::getline(split, part, ':')) {
        parts.push_back(part);
      }
      int tube = tube_index(parts[0]);
      if (tube < 0) {
        return false;
      }
      SheetTube& sheet = tubes_[static_cast<std::size_t>(tube)];
      int node = 0;
      if (parts.size() == 2 && parts[1] == "inlet") {
        node = 0;
        sheet.set_inlet(SheetTube::End::connected);
      } else if (parts.size() == 2 && parts[1] == "outlet") {
        node = sheet.count() - 1;
        sheet.set_outlet(SheetTube::End::connected);
      } else {
        int section = std::stoi(parts[1]);
        double position = std::stod(parts[2]);
        node = sheet.band_first(section) + static_cast<int>(std::lround(position * (sheet.band_points(section) - 1)));
      }
      *tube_fields[side] = tube;
      *port_fields[side] = sheet.add_port(node);
    }
    configure(coupling);
    couplings_.push_back(coupling);
  }
  return true;
}

int SheetSystem::tube_index(const std::string& id) const {
  for (std::size_t index = 0; index < ids_.size(); ++index) {
    if (ids_[index] == id) {
      return static_cast<int>(index);
    }
  }
  return -1;
}

void SheetSystem::configure(Coupling& coupling) {
  SheetTube& a = tubes_[static_cast<std::size_t>(coupling.tube_a)];
  SheetTube& b = tubes_[static_cast<std::size_t>(coupling.tube_b)];
  int i = a.port(coupling.port_a).index;
  int j = b.port(coupling.port_b).index;
  double distance = 0.5 * (a.node_dx(i) + b.node_dx(j));
  coupling.air_delay = std::max(1.0, distance / kSound * rate_);
  double speed = std::sqrt((a.axial_tie(std::min(i, a.count() - 2)) * a.node_dx(i) * a.node_dx(i) / a.mass(i) +
                            b.axial_tie(std::min(j, b.count() - 2)) * b.node_dx(j) * b.node_dx(j) / b.mass(j)) /
                           2.0);
  coupling.solid_delay = std::max(1.0, distance / std::max(0.1, speed) * rate_);
  for (int k = 0; k < 3; ++k) {
    coupling.taps[k] = std::sqrt(a.fir(i)[k] * b.fir(j)[k]);
  }
  int capacity = power_of_two_at_least(std::max(coupling.air_delay, coupling.solid_delay / std::sqrt(0.9)) + 5.0);
  coupling.mask = std::max(capacity, 256) - 1;
  for (std::vector<double>* history : {&coupling.air_ab, &coupling.air_ba, &coupling.solid_ab, &coupling.solid_ba}) {
    history->assign(static_cast<std::size_t>(coupling.mask + 1), 0.0);
  }
}

double SheetSystem::filter(const Coupling& coupling, const std::vector<double>& history, double delay,
                           const double* taps) const {
  int whole = static_cast<int>(std::floor(delay));
  double f = delay - whole;
  int p = coupling.pos - whole;
  int m = coupling.mask;
  double h0 = history[static_cast<std::size_t>(p & m)];
  double h1 = history[static_cast<std::size_t>((p - 1) & m)];
  double h2 = history[static_cast<std::size_t>((p - 2) & m)];
  double h3 = history[static_cast<std::size_t>((p - 3) & m)];
  return (1.0 - f) * (taps[0] * h0 + taps[1] * h1 + taps[2] * h2) + f * (taps[0] * h1 + taps[1] * h2 + taps[2] * h3);
}

void SheetSystem::prepare(Coupling& coupling) {
  SheetTube& a = tubes_[static_cast<std::size_t>(coupling.tube_a)];
  SheetTube& b = tubes_[static_cast<std::size_t>(coupling.tube_b)];
  SheetPort& pa = a.port(coupling.port_a);
  SheetPort& pb = b.port(coupling.port_b);
  double radius = std::max(kFloor, std::min(a.rest_radius(pa.index) + a.displacement(pa.index),
                                            b.rest_radius(pb.index) + b.displacement(pb.index)));
  coupling.air_scale = std::sqrt(kPiValue / (kSectors * kRho * kSound)) * radius;
  pa.air_y = coupling.air_scale * coupling.air_scale;
  pb.air_y = pa.air_y;
  double factor = std::sqrt(a.factor(pa.index) * b.factor(pb.index));
  double z = std::sqrt(a.impedance(pa.index) / std::sqrt(a.factor(pa.index)) * b.impedance(pb.index) /
                       std::sqrt(b.factor(pb.index))) *
             std::sqrt(factor);
  coupling.solid_scale = std::sqrt(z);
  pa.solid_z = z;
  pb.solid_z = z;
  pa.air_in = filter(coupling, coupling.air_ba, coupling.air_delay, kAirTaps) / coupling.air_scale;
  pb.air_in = filter(coupling, coupling.air_ab, coupling.air_delay, kAirTaps) / coupling.air_scale;
  double delay = std::max(1.0, coupling.solid_delay / std::sqrt(factor));
  pa.solid_in = filter(coupling, coupling.solid_ba, delay, coupling.taps) / coupling.solid_scale;
  pb.solid_in = filter(coupling, coupling.solid_ab, delay, coupling.taps) / coupling.solid_scale;
}

void SheetSystem::send(Coupling& coupling) {
  SheetPort& pa = tubes_[static_cast<std::size_t>(coupling.tube_a)].port(coupling.port_a);
  SheetPort& pb = tubes_[static_cast<std::size_t>(coupling.tube_b)].port(coupling.port_b);
  std::size_t at = static_cast<std::size_t>(coupling.pos);
  coupling.air_ab[at] = pa.air_out * coupling.air_scale;
  coupling.air_ba[at] = pb.air_out * coupling.air_scale;
  coupling.solid_ab[at] = pa.solid_out * coupling.solid_scale;
  coupling.solid_ba[at] = pb.solid_out * coupling.solid_scale;
  coupling.pos = (coupling.pos + 1) & coupling.mask;
}

double SheetSystem::step() {
  if (fault_) {
    return 0.0;
  }
  for (Coupling& coupling : couplings_) {
    prepare(coupling);
  }
  for (SheetTube& tube : tubes_) {
    if (!tube.step()) {
      fault_ = true;
      return 0.0;
    }
  }
  for (Coupling& coupling : couplings_) {
    send(coupling);
  }
  SheetTube& out = tubes_[static_cast<std::size_t>(output_)];
  return out.pressure(out.count() - 1);
}

// --- SheetJoint ------------------------------------------------------------

void SheetJoint::connect(SheetTube& a, int node_a, SheetTube& b, int node_b, double rate) {
  a_ = &a;
  b_ = &b;
  port_a_ = a.add_port(node_a);
  port_b_ = b.add_port(node_b);
  double distance = 0.5 * (a.node_dx(node_a) + b.node_dx(node_b));
  air_delay_ = std::max(1.0, distance / kSound * rate);
  double speed = std::sqrt((a.axial_tie(std::min(node_a, a.count() - 2)) * a.node_dx(node_a) * a.node_dx(node_a) /
                                a.mass(node_a) +
                            b.axial_tie(std::min(node_b, b.count() - 2)) * b.node_dx(node_b) * b.node_dx(node_b) /
                                b.mass(node_b)) /
                           2.0);
  solid_delay_ = std::max(1.0, distance / std::max(0.1, speed) * rate);
  for (int k = 0; k < 3; ++k) {
    taps_[k] = std::sqrt(a.fir(node_a)[k] * b.fir(node_b)[k]);
  }
  // Room for the tissue delay at a fifth of the elasticity.
  int capacity = power_of_two_at_least(std::max(air_delay_, 2.5 * solid_delay_) + 5.0);
  mask_ = std::max(capacity, 256) - 1;
  for (std::vector<double>* history : {&air_ab_, &air_ba_, &solid_ab_, &solid_ba_}) {
    history->assign(static_cast<std::size_t>(mask_ + 1), 0.0);
  }
  pos_ = 0;
}

double SheetJoint::filter(const std::vector<double>& history, double delay, const double* taps) const {
  int whole = static_cast<int>(std::floor(delay));
  double f = delay - whole;
  int p = pos_ - whole;
  int m = mask_;
  double h0 = history[static_cast<std::size_t>(p & m)];
  double h1 = history[static_cast<std::size_t>((p - 1) & m)];
  double h2 = history[static_cast<std::size_t>((p - 2) & m)];
  double h3 = history[static_cast<std::size_t>((p - 3) & m)];
  return (1.0 - f) * (taps[0] * h0 + taps[1] * h1 + taps[2] * h2) + f * (taps[0] * h1 + taps[1] * h2 + taps[2] * h3);
}

void SheetJoint::prepare() {
  SheetPort& pa = a_->port(port_a_);
  SheetPort& pb = b_->port(port_b_);
  double radius = 0.0;
  if (opening_ < 0.0) {
    radius = std::max(kFloor, std::min(a_->rest_radius(pa.index) + a_->displacement(pa.index),
                                       b_->rest_radius(pb.index) + b_->displacement(pb.index)));
  } else {
    radius = std::sqrt(opening_ / kPiValue);
  }
  // A shut port passes nothing; its waves still decay through the histories.
  air_scale_ = std::sqrt(kPiValue / (kSectors * kRho * kSound)) * std::max(radius, 1e-7);
  pa.air_y = air_scale_ * air_scale_;
  pb.air_y = pa.air_y;
  double factor = std::sqrt(a_->factor(pa.index) * b_->factor(pb.index));
  double z = std::sqrt(a_->impedance(pa.index) / std::sqrt(a_->factor(pa.index)) * b_->impedance(pb.index) /
                       std::sqrt(b_->factor(pb.index))) *
             std::sqrt(factor);
  solid_scale_ = std::sqrt(z);
  pa.solid_z = z;
  pb.solid_z = z;
  pa.air_in = filter(air_ba_, air_delay_, kAirTaps) / air_scale_;
  pb.air_in = filter(air_ab_, air_delay_, kAirTaps) / air_scale_;
  double delay = std::min(std::max(1.0, solid_delay_ / std::sqrt(factor)), static_cast<double>(mask_ - 4));
  pa.solid_in = filter(solid_ba_, delay, taps_) / solid_scale_;
  pb.solid_in = filter(solid_ab_, delay, taps_) / solid_scale_;
}

void SheetJoint::send() {
  SheetPort& pa = a_->port(port_a_);
  SheetPort& pb = b_->port(port_b_);
  std::size_t at = static_cast<std::size_t>(pos_);
  air_ab_[at] = pa.air_out * air_scale_;
  air_ba_[at] = pb.air_out * air_scale_;
  solid_ab_[at] = pa.solid_out * solid_scale_;
  solid_ba_[at] = pb.solid_out * solid_scale_;
  pos_ = (pos_ + 1) & mask_;
}

}  // namespace tract
