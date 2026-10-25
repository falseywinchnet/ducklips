// The sustain tool's run (see tract/sustain_run.hpp). The arithmetic of
// the step loop is tools/sustain.cpp's as it stood (2026-10-20), moved here
// unchanged; this file is compiled with that tool's flags (no
// -ffp-contract=fast, no -mcpu=native: CMakeLists.txt, target tract_run),
// so the interpolations it does round as they did there.

#include "tract/sustain_run.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>

namespace tract {

DecimationSection decimation_lowpass(double corner, double rate, double q) {
  double w = std::tan(kPi * corner / rate);
  double norm = 1.0 / (1.0 + w / q + w * w);
  DecimationSection filter;
  filter.b0 = w * w * norm;
  filter.b1 = 2.0 * filter.b0;
  filter.b2 = filter.b0;
  filter.a1 = 2.0 * (w * w - 1.0) * norm;
  filter.a2 = (1.0 - w / q + w * w) * norm;
  return filter;
}

double run_section(DecimationSection& filter, double input) {
  double output = filter.b0 * input + filter.z1;
  filter.z1 = filter.b1 * input - filter.a1 * output + filter.z2;
  filter.z2 = filter.b2 * input - filter.a2 * output;
  return output;
}

void write_wave_file(const char* path, const std::vector<float>& samples, int rate) {
  std::ofstream out(path, std::ios::binary);
  std::uint32_t data_bytes = static_cast<std::uint32_t>(samples.size() * 2);
  out.write("RIFF", 4);
  std::uint32_t chunk = 36 + data_bytes;
  out.write(reinterpret_cast<const char*>(&chunk), 4);
  out.write("WAVEfmt ", 8);
  std::uint32_t format_size = 16;
  std::uint16_t format = 1;
  std::uint16_t channels = 1;
  std::uint32_t sample_rate = static_cast<std::uint32_t>(rate);
  std::uint32_t byte_rate = sample_rate * 2;
  std::uint16_t align = 2;
  std::uint16_t bits = 16;
  out.write(reinterpret_cast<const char*>(&format_size), 4);
  out.write(reinterpret_cast<const char*>(&format), 2);
  out.write(reinterpret_cast<const char*>(&channels), 2);
  out.write(reinterpret_cast<const char*>(&sample_rate), 4);
  out.write(reinterpret_cast<const char*>(&byte_rate), 4);
  out.write(reinterpret_cast<const char*>(&align), 2);
  out.write(reinterpret_cast<const char*>(&bits), 2);
  out.write("data", 4);
  out.write(reinterpret_cast<const char*>(&data_bytes), 4);
  for (float sample : samples) {
    double clipped = std::fmax(-1.0, std::fmin(1.0, static_cast<double>(sample)));
    std::int16_t value = static_cast<std::int16_t>(std::lround(clipped * 32767.0));
    out.write(reinterpret_cast<const char*>(&value), 2);
  }
}

const char* ProcessEnvironment::find(const char* name) const { return std::getenv(name); }

void SustainSink::checked(long, VoiceBody&) {}

namespace {

bool read_profile(const std::string& vowel, std::vector<double>& areas, double& section_length) {
  std::ifstream input("data/story_2008_areas.txt");
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream fields(line);
    std::string name;
    fields >> name;
    if (name != vowel) {
      continue;
    }
    fields >> section_length;
    double area = 0.0;
    while (fields >> area) {
      areas.push_back(area);
    }
    return !areas.empty();
  }
  return false;
}

// A variable's value as the command line reads it (strtod of the text),
// into `field` when the variable is set.
void take_number(const EnvironmentLookup& environment, const char* name, std::optional<double>& field) {
  const char* text = environment.find(name);
  if (text != nullptr) {
    field = std::strtod(text, nullptr);
  }
}

void take_number(const EnvironmentLookup& environment, const char* name, double& field) {
  const char* text = environment.find(name);
  if (text != nullptr) {
    field = std::strtod(text, nullptr);
  }
}

void take_path(const EnvironmentLookup& environment, const char* name, std::string& field) {
  const char* text = environment.find(name);
  if (text != nullptr) {
    field = text;
  }
}

// "first,second" as the command line reads it.
bool pair_of(const EnvironmentLookup& environment, const char* name, double& first, double& second) {
  const char* text = environment.find(name);
  if (text == nullptr) {
    return false;
  }
  first = std::strtod(text, nullptr);
  second = std::strtod(std::strchr(text, ',') + 1, nullptr);
  return true;
}

bool format_error(std::string& error, int& status, const char* text) {
  error = text;
  status = 1;
  return false;
}

}  // namespace

bool settings_from_environment(const std::vector<std::string>& arguments, const EnvironmentLookup& environment,
                               SustainSettings& settings, std::string& error, int& status) {
  if (arguments.size() < 2) {
    format_error(error, status, "usage: sustain VOWEL OUT.wav [pitch_hz] [lung_pa] [seconds]\n");
    status = 2;
    return false;
  }
  settings.vowel = arguments[0];
  settings.wave_path = arguments[1];
  settings.areas.clear();
  // AREAFILE=path: the tube's profile instead of a Story vowel: the section
  // length (cm) on the first line, then one area (cm^2) per line.
  if (environment.find("AREAFILE") != nullptr) {
    std::ifstream area_file(environment.find("AREAFILE"));
    double area = 0.0;
    settings.section_length = 0.0;
    area_file >> settings.section_length;
    while (area_file >> area) {
      settings.areas.push_back(area);
    }
    if (settings.areas.empty() || settings.section_length <= 0.0) {
      char text[1024];
      std::snprintf(text, sizeof(text), "cannot read AREAFILE %s\n", environment.find("AREAFILE"));
      return format_error(error, status, text);
    }
  } else if (!read_profile(settings.vowel, settings.areas, settings.section_length)) {
    char text[1024];
    std::snprintf(text, sizeof(text), "no profile %s in data/story_2008_areas.txt\n", settings.vowel.c_str());
    return format_error(error, status, text);
  }
  if (!options_from_environment(environment, settings, error)) {
    status = 1;
    return false;
  }
  settings.pitch = arguments.size() > 2 ? std::strtod(arguments[2].c_str(), nullptr) : 120.0;
  settings.lung = arguments.size() > 3 ? std::strtod(arguments[3].c_str(), nullptr) : 800.0;
  settings.seconds = arguments.size() > 4 ? std::strtod(arguments[4].c_str(), nullptr) : 1.5;
  status = 0;
  return true;
}

bool options_from_environment(const EnvironmentLookup& environment, SustainSettings& settings, std::string& error) {
  int status = 0;
  // AREAFILE2=path with GLIDE=start,duration (s).
  if (environment.find("AREAFILE2") != nullptr && environment.find("GLIDE") != nullptr) {
    std::ifstream area_file(environment.find("AREAFILE2"));
    double second_length = 0.0;
    double area = 0.0;
    area_file >> second_length;
    settings.areas_to.clear();
    while (area_file >> area) {
      settings.areas_to.push_back(area);
    }
    settings.glide = true;
    pair_of(environment, "GLIDE", settings.glide_start, settings.glide_duration);
    if (settings.areas_to.size() != settings.areas.size() || settings.glide_duration <= 0.0) {
      char text[1024];
      std::snprintf(text, sizeof(text), "AREAFILE2 needs %zu areas and GLIDE=start,duration\n", settings.areas.size());
      return format_error(error, status, text);
    }
  }
  // AREASEQ=path: keyframes, one per line: a time (s) then the node areas.
  if (environment.find("AREASEQ") != nullptr) {
    std::ifstream sequence_file(environment.find("AREASEQ"));
    std::string line;
    while (std::getline(sequence_file, line)) {
      std::istringstream fields(line);
      double time = 0.0;
      if (!(fields >> time)) {
        continue;
      }
      std::vector<double> frame;
      double area = 0.0;
      while (fields >> area) {
        frame.push_back(area);
      }
      if (frame.size() != settings.areas.size()) {
        char text[1024];
        std::snprintf(text, sizeof(text), "AREASEQ keyframe at %g s has %zu areas, not %zu\n", time, frame.size(),
                      settings.areas.size());
        return format_error(error, status, text);
      }
      settings.key_times.push_back(time);
      settings.key_areas.insert(settings.key_areas.end(), frame.begin(), frame.end());
    }
  }
  if (environment.find("VELUMSEQ") != nullptr) {
    std::ifstream velum_file(environment.find("VELUMSEQ"));
    double time = 0.0;
    double area = 0.0;
    while (velum_file >> time >> area) {
      settings.velum_times.push_back(time);
      settings.velum_areas.push_back(area);
    }
  }
  if (environment.find("FOLDSEQ") != nullptr) {
    std::ifstream fold_file(environment.find("FOLDSEQ"));
    double time = 0.0;
    CordPosture fold_posture;
    double move = 0.0;
    while (fold_file >> time >> fold_posture.lower_rest >> fold_posture.upper_rest >> fold_posture.chink >> move) {
      settings.fold_times.push_back(time);
      settings.fold_postures.push_back(fold_posture);
      settings.fold_moves.push_back(move);
    }
  }
  if (environment.find("JETSEQ") != nullptr) {
    std::ifstream jet_file(environment.find("JETSEQ"));
    double time = 0.0;
    double gate = 0.0;
    while (jet_file >> time >> gate) {
      settings.jet_times.push_back(time);
      settings.jet_gates.push_back(gate);
    }
  }
  // AREATICKS=path: one line a tick, the node areas; VELUMTICKS=path: one
  // area a line; FOLDTICKS=path: lower rest, upper rest, chink a line.
  if (environment.find("AREATICKS") != nullptr) {
    std::ifstream tick_file(environment.find("AREATICKS"));
    std::string line;
    while (std::getline(tick_file, line)) {
      std::istringstream fields(line);
      std::vector<double> frame;
      double area = 0.0;
      while (fields >> area) {
        frame.push_back(area);
      }
      if (frame.empty()) {
        continue;
      }
      if (frame.size() != settings.areas.size()) {
        char text[1024];
        std::snprintf(text, sizeof(text), "AREATICKS row %zu has %zu areas, not %zu\n",
                      settings.tick_areas.size() / settings.areas.size(), frame.size(), settings.areas.size());
        return format_error(error, status, text);
      }
      settings.tick_areas.insert(settings.tick_areas.end(), frame.begin(), frame.end());
    }
  }
  if (environment.find("VELUMTICKS") != nullptr) {
    std::ifstream tick_file(environment.find("VELUMTICKS"));
    double area = 0.0;
    while (tick_file >> area) {
      settings.tick_velum.push_back(area);
    }
  }
  if (environment.find("FOLDTICKS") != nullptr) {
    std::ifstream tick_file(environment.find("FOLDTICKS"));
    CordPosture tick_posture;
    while (tick_file >> tick_posture.lower_rest >> tick_posture.upper_rest >> tick_posture.chink) {
      settings.tick_folds.push_back(tick_posture);
    }
  }
  if (environment.find("FOLDKEYS") != nullptr) {
    std::ifstream key_file(environment.find("FOLDKEYS"));
    double time = 0.0;
    CordPosture key_posture;
    while (key_file >> time >> key_posture.lower_rest >> key_posture.upper_rest >> key_posture.chink) {
      settings.fold_key_times.push_back(time);
      settings.fold_key_postures.push_back(key_posture);
    }
  }
  take_number(environment, "PITCHEND", settings.pitch_end);
  take_number(environment, "RELEASE", settings.release);
  take_number(environment, "TRACE", settings.trace_ms);
  if (environment.find("BRANCHES") != nullptr) {
    settings.branches = std::string(environment.find("BRANCHES"));
  }
  take_number(environment, "ASPIRATION", settings.aspiration);
  take_number(environment, "FRICATION", settings.frication);
  take_number(environment, "FLUX", settings.flux);
  take_number(environment, "CLOSING", settings.closing);
  take_number(environment, "TREMOR", settings.tremor);
  take_number(environment, "ADDTREMOR", settings.adductor_tremor);
  take_number(environment, "BREATHTREMOR", settings.breath_tremor);
  take_number(environment, "PITCHFB", settings.pitch_feedback);
  take_number(environment, "PORTRESIST", settings.port_resistance);
  take_number(environment, "REVJET", settings.reversible_jet);
  take_number(environment, "REVJETCLOSE", settings.reversible_jet_closing);
  take_number(environment, "FOLDTURB", settings.fold_turbulence);
  take_number(environment, "SKIN", settings.skin);
  take_number(environment, "WALLS", settings.wall_scale);
  take_number(environment, "WALLMASS", settings.mass_scale);
  take_number(environment, "AXIAL", settings.axial_scale);
  take_number(environment, "WALLVAR", settings.wall_variation);
  take_number(environment, "VISCOUS", settings.viscous);
  take_number(environment, "ACTIVELENGTH", settings.active_length);
  take_number(environment, "STATICAREA", settings.static_area);
  take_number(environment, "REVNOSE", settings.reversible_nose);
  take_number(environment, "SEAL", settings.seal);
  take_number(environment, "SEALBLEND", settings.seal_blend);
  take_number(environment, "JETBOUND", settings.jet_bound);
  if (pair_of(environment, "NOSERAD", settings.nose_radiation_mass, settings.nose_radiation_resistance)) {
    settings.nose_radiation = true;
  }
  take_number(environment, "NOSEVISCOUS", settings.nose_viscous);
  take_number(environment, "NOSESHELF", settings.nose_shelf);
  take_number(environment, "NOSEHOP", settings.nose_hop);
  take_number(environment, "HOPLOSS", settings.hop_loss);
  if (pair_of(environment, "REST", settings.rest_lower, settings.rest_upper)) {
    settings.rest = true;
  }
  take_number(environment, "CHINK", settings.chink);
  if (pair_of(environment, "FOLDDAMP", settings.fold_damping_lower, settings.fold_damping_upper)) {
    settings.fold_damping = true;
  }
  take_number(environment, "VELUM", settings.velum);
  take_path(environment, "SOURCE", settings.source_path);
  take_path(environment, "RAW", settings.raw_path);
  if (environment.find("FORCEWAVE") != nullptr) {
    std::ifstream wave_file(environment.find("FORCEWAVE"), std::ios::binary);
    double wave = 0.0;
    while (wave_file.read(reinterpret_cast<char*>(&wave), sizeof(wave))) {
      settings.forced_waves.push_back(wave);
    }
  }
  take_number(environment, "OBSERVE", settings.observe_gain);
  if (environment.find("FORCEFLOW") != nullptr) {
    std::ifstream flow_file(environment.find("FORCEFLOW"), std::ios::binary);
    double flow = 0.0;
    while (flow_file.read(reinterpret_cast<char*>(&flow), sizeof(flow))) {
      settings.observed_flows.push_back(flow);
    }
  }
  if (environment.find("RAWDOUBLE") != nullptr) {
    settings.raw_double = std::strtod(environment.find("RAWDOUBLE"), nullptr) != 0.0;
  }
  if (environment.find("RAWCORD") != nullptr) {
    settings.raw_cord = true;
  }
  take_path(environment, "PORTWAVES", settings.port_waves_path);
  take_path(environment, "PORTDUMP", settings.port_dump_path);
  take_path(environment, "GEODUMP", settings.geometry_path);
  take_path(environment, "AREAS", settings.areas_path);
  take_path(environment, "WALLDUMP", settings.wall_dump_path);
  take_path(environment, "NOSEDUMP", settings.nose_dump_path);
  return status == 0;
}

bool validate_settings(const SustainSettings& settings, std::string& error) {
  if (settings.areas.empty() || settings.section_length <= 0.0) {
    error = "the tube needs areas and a section length\n";
    return false;
  }
  std::size_t count = settings.areas.size();
  if (settings.key_areas.size() != settings.key_times.size() * count) {
    error = "AREASEQ keyframes need the tube's count of areas each\n";
    return false;
  }
  if (settings.tick_areas.size() % count != 0) {
    error = "AREATICKS rows need the tube's count of areas each\n";
    return false;
  }
  if (settings.glide && (settings.areas_to.size() != count || settings.glide_duration <= 0.0)) {
    error = "AREAFILE2 needs the tube's count of areas and GLIDE=start,duration\n";
    return false;
  }
  if (settings.velum_times.size() != settings.velum_areas.size() ||
      settings.jet_times.size() != settings.jet_gates.size() ||
      settings.fold_times.size() != settings.fold_postures.size() ||
      settings.fold_times.size() != settings.fold_moves.size() ||
      settings.fold_key_times.size() != settings.fold_key_postures.size()) {
    error = "keyframe times and values differ in count\n";
    return false;
  }
  return true;
}

void port_wave_row(VoiceBody& body, double row[6]) {
  Index port_node = body.nose_port_node();
  row[0] = body.nose().left_incoming();
  row[1] = body.vocal().side_incoming(0);
  row[2] = body.vocal().pressure(port_node);
  row[3] = body.vocal().backward_wave(port_node);
  row[4] = body.nose().radiated_flow();
  row[5] = body.vocal().radiated_flow();
}

void geometry_row(Tube& vocal, Index node, double row[12]) {
  for (int column = 0; column < 12; ++column) {
    row[column] = 0.0;
  }
  row[0] = vocal.area(node);
  row[1] = vocal.frame_flow(node);
  row[2] = node + 1 < vocal.count() ? vocal.link_resistance(node) : 0.0;
  vocal.lining(node, row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11]);
}

StreamFiles::StreamFiles(const SustainSettings& settings, bool trace) : raw_double_(settings.raw_double) {
  if (!settings.raw_path.empty()) {
    raw_ = std::fopen(settings.raw_path.c_str(), "wb");
  }
  if (!settings.port_waves_path.empty()) {
    port_waves_ = std::fopen(settings.port_waves_path.c_str(), "wb");
  }
  if (!settings.port_dump_path.empty()) {
    port_dump_ = std::fopen(settings.port_dump_path.c_str(), "wb");
  }
  if (!settings.geometry_path.empty()) {
    geometry_ = std::fopen(settings.geometry_path.c_str(), "wb");
  }
  trace_every_ = trace && settings.trace_ms > 0.0 ? static_cast<long>(settings.trace_ms * 1e-3 * kStepRate) : 0;
}

StreamFiles::~StreamFiles() {
  for (std::FILE* file : {raw_, port_waves_, port_dump_, geometry_}) {
    if (file != nullptr) {
      std::fclose(file);
    }
  }
}

void StreamFiles::record(long, VoiceBody& body, const StepFlows& flows) {
  // GEODUMP=path: per step, each node's area, the air the frame draws
  // from it, its link's viscous resistance (0 for the last) and its
  // lining (WALLDUMP's nine), float64, for a replica of a moving tube.
  if (geometry_ != nullptr) {
    Tube& vocal = body.vocal();
    for (Index node = 0; node < vocal.count(); ++node) {
      double row[12];
      geometry_row(vocal, node, row);
      std::fwrite(row, sizeof(row), 1, geometry_);
    }
  }
  // PORTWAVES=path: per step (float64) the wave the port delivered into the
  // nose, the wave it sent into the vocal node, that node's pressure, the
  // wave arriving at it from its outer link next step, the nostrils' and
  // the lips' flows.
  if (port_waves_ != nullptr) {
    double row[6];
    port_wave_row(body, row);
    std::fwrite(row, sizeof(row), 1, port_waves_);
  }
  // PORTDUMP=path: the nasal port's area each step (float64).
  if (port_dump_ != nullptr) {
    double port = body.velum_port();
    std::fwrite(&port, sizeof(port), 1, port_dump_);
  }
  if (raw_ != nullptr) {
    if (raw_double_) {
      // RAWDOUBLE=1: float64 pairs (a search's differences of 1e-3 in a
      // log area drowned in float32's rounding near -70 dB).
      double pair[2] = {flows.entering, flows.radiated};
      std::fwrite(pair, sizeof(pair), 1, raw_);
    } else {
      float pair[2] = {static_cast<float>(flows.entering), static_cast<float>(flows.radiated)};
      std::fwrite(pair, sizeof(pair), 1, raw_);
    }
  }
}

void StreamFiles::checked(long step, VoiceBody& body) {
  if (trace_every_ > 0 && step % trace_every_ == 0) {
    Tube& trachea = body.trachea();
    Tube& vocal = body.vocal();
    std::printf("%.4f area %.4f flow %7.1f contact %.2f sub %7.0f supra %7.0f lips %8.1f wall %8.2f nose_p %8.1f nose %8.1f\n",
                step / kStepRate, body.cord().area(), body.cord().flow(), body.cord().contact(),
                trachea.pressure(trachea.count() - 1), vocal.pressure(0), vocal.radiated_flow(),
                vocal.total_wall_flow(), body.nose().pressure(0), body.nose().radiated_flow());
  }
}

namespace {

BodyLayout layout_of(const SustainSettings& settings) {
  BodyLayout layout = vowel_layout(settings.areas, settings.section_length);
  if (settings.branches.has_value()) {
    // BRANCHES=pus: keep piriform (p), under-tongue (u), sinuses (s).
    const std::string& keep = *settings.branches;
    layout.piriform = keep.find('p') != std::string::npos;
    layout.under_tongue = keep.find('u') != std::string::npos;
    layout.sinuses = keep.find('s') != std::string::npos;
  }
  if (settings.aspiration.has_value()) {
    layout.aspiration = *settings.aspiration;
  }
  if (settings.frication.has_value()) {
    layout.frication = *settings.frication;
  }
  if (settings.flux.has_value()) {
    layout.jet_flux = *settings.flux;
  }
  if (settings.closing.has_value()) {
    layout.closing_turbulence = *settings.closing;
  }
  if (settings.tremor.has_value()) {
    layout.tremor = *settings.tremor;
  }
  if (settings.adductor_tremor.has_value()) {
    layout.adductor_tremor = *settings.adductor_tremor;
  }
  if (settings.breath_tremor.has_value()) {
    layout.breath_tremor = *settings.breath_tremor;
  }
  if (settings.pitch_feedback.has_value()) {
    layout.pitch_feedback = *settings.pitch_feedback;
  }
  if (settings.port_resistance.has_value()) {
    layout.port_resistance_scale = *settings.port_resistance;
  }
  if (settings.reversible_jet.has_value()) {
    layout.reversible_jet = *settings.reversible_jet;
  }
  if (settings.reversible_jet_closing.has_value()) {
    layout.reversible_jet_closing = *settings.reversible_jet_closing;
  }
  if (settings.fold_turbulence.has_value()) {
    layout.fold_turbulence = *settings.fold_turbulence;
  }
  if (settings.skin.has_value()) {
    layout.skin_gain = *settings.skin;
  }
  // WALLS, WALLMASS and AXIAL scale the walls (1: as they are).
  for (SegmentMaterial& material : layout.vocal_material_list) {
    material.axial_tension *= settings.axial_scale;
    material.axial_tension *= settings.wall_scale;
    material.hoop_tension *= settings.wall_scale;
    material.coupling_stiffness *= settings.wall_scale;
    material.support_stiffness *= settings.wall_scale;
    material.cover_thickness *= settings.mass_scale;
    material.body_depth *= settings.mass_scale;
  }
  if (settings.wall_variation.has_value()) {
    layout.wall_variation = *settings.wall_variation;
  }
  if (settings.viscous.has_value()) {
    layout.viscous_scale = *settings.viscous;
  }
  return layout;
}

// The keyframe interpolation the command line does at each 1 ms tick: the
// keyframes around `now` (held before the first and after the last) and
// the share between them.
void keyframe_share(const std::vector<double>& times, double now, std::size_t& from, std::size_t& to, double& share) {
  std::size_t after = 0;
  while (after < times.size() && times[after] <= now) {
    ++after;
  }
  from = after == 0 ? 0 : after - 1;
  to = after < times.size() ? after : times.size() - 1;
  share = to == from ? 0.0 : std::clamp((now - times[from]) / (times[to] - times[from]), 0.0, 1.0);
}

}  // namespace

SustainRun::SustainRun(const SustainSettings& settings) : body_(layout_of(settings)) {
  // ACTIVELENGTH=cm: the vocal tube's length (Tube::set_active_length).
  if (settings.active_length.has_value()) {
    body_.vocal().set_active_length(*settings.active_length);
  }
  // STATICAREA=1: the air sees the rest geometry only.
  if (settings.static_area.has_value()) {
    body_.vocal().set_static_admittance(*settings.static_area != 0.0);
  }
  // REVNOSE=1: the nose as the reversible engine needs it.
  if (settings.reversible_nose.has_value() && *settings.reversible_nose != 0.0) {
    body_.nose().scale_axial(0.0);
    body_.nose().set_static_admittance(true);
  }
  if (settings.seal.has_value()) {
    body_.vocal().set_seal_area(*settings.seal);
  }
  if (settings.seal_blend.has_value() && *settings.seal_blend != 0.0) {
    body_.vocal().set_seal_blend(true);
  }
  if (settings.jet_bound.has_value() && *settings.jet_bound != 0.0) {
    body_.set_jet_bound(std::sqrt(2.0 * settings.lung * kPascal / kAirDensity));
  }
  if (settings.nose_radiation) {
    body_.nose().set_radiation_scales(settings.nose_radiation_mass, settings.nose_radiation_resistance);
  }
  if (settings.nose_viscous.has_value()) {
    body_.nose().set_viscous_scale(*settings.nose_viscous);
  }
  if (settings.nose_shelf.has_value()) {
    body_.nose().set_wave_shelf(*settings.nose_shelf);
  }
  if (settings.nose_hop.has_value()) {
    body_.nose().set_wave_smoothing(*settings.nose_hop);
  }
  if (settings.hop_loss.has_value()) {
    body_.vocal().set_wave_smoothing(*settings.hop_loss);
  }
  body_.set_pitch(settings.pitch);
  CordPosture posture = modal_posture();
  if (settings.rest) {
    posture.lower_rest = settings.rest_lower;
    posture.upper_rest = settings.rest_upper;
  }
  if (settings.chink.has_value()) {
    posture.chink = *settings.chink;
  }
  body_.cord().set_target(posture, 0.001);
  if (settings.fold_damping) {
    body_.cord().set_damping(settings.fold_damping_lower, settings.fold_damping_upper);
  }
  body_.set_lung_pressure(settings.lung, 0.08);
  if (settings.velum.has_value()) {
    body_.set_velum(*settings.velum);
  }
  first_ = decimation_lowpass(18000.0, kStepRate, 0.5412);
  second_ = decimation_lowpass(18000.0, kStepRate, 1.3066);
  source_first_ = decimation_lowpass(18000.0, kStepRate, 0.5412);
  source_second_ = decimation_lowpass(18000.0, kStepRate, 1.3066);
  total_ = total_steps(settings);
  release_step_ = release_step(settings);
  keep_audio_ = !settings.wave_path.empty();
  keep_source_ = !settings.source_path.empty();
  if (keep_audio_) {
    audio_.reserve(static_cast<std::size_t>(total_ / 2 + 1));
  }
}

long trigger_step(double time) { return std::lround(time * kStepRate); }

long total_steps(const SustainSettings& settings) { return static_cast<long>(settings.seconds * kStepRate); }

long release_step(const SustainSettings& settings) {
  long total = total_steps(settings);
  return settings.release > 0.0 ? total - static_cast<long>(settings.release * kStepRate) : total;
}

void tick_controls(const SustainSettings& settings, long step, long total, TickControls& controls) {
  double pitch = settings.pitch;
  double pitch_end = settings.pitch_end.has_value() ? *settings.pitch_end : pitch;
  double progress = static_cast<double>(step) / static_cast<double>(total);
  controls.pitch = pitch + (pitch_end - pitch) * progress;
  double now = step / kStepRate;
  std::size_t from = 0;
  std::size_t to = 0;
  double share = 0.0;
  // The field's row at this tick (the last row after the field ends).
  std::size_t tick = static_cast<std::size_t>(step / 96);
  controls.fold_key = !settings.fold_key_times.empty() || !settings.tick_folds.empty();
  if (!settings.tick_folds.empty()) {
    controls.fold_posture = settings.tick_folds[std::min(tick, settings.tick_folds.size() - 1)];
  } else if (controls.fold_key) {
    keyframe_share(settings.fold_key_times, now, from, to, share);
    const CordPosture& start = settings.fold_key_postures[from];
    const CordPosture& end = settings.fold_key_postures[to];
    controls.fold_posture.lower_rest = start.lower_rest + share * (end.lower_rest - start.lower_rest);
    controls.fold_posture.upper_rest = start.upper_rest + share * (end.upper_rest - start.upper_rest);
    controls.fold_posture.chink = start.chink + share * (end.chink - start.chink);
  }
  controls.jet = !settings.jet_times.empty();
  if (controls.jet) {
    keyframe_share(settings.jet_times, now, from, to, share);
    controls.jet_gate = settings.jet_gates[from] + share * (settings.jet_gates[to] - settings.jet_gates[from]);
  }
  controls.velum = !settings.velum_times.empty() || !settings.tick_velum.empty();
  if (!settings.tick_velum.empty()) {
    controls.velum_area = settings.tick_velum[std::min(tick, settings.tick_velum.size() - 1)];
  } else if (controls.velum) {
    keyframe_share(settings.velum_times, now, from, to, share);
    controls.velum_area = settings.velum_areas[from] + share * (settings.velum_areas[to] - settings.velum_areas[from]);
  }
  std::size_t count = settings.areas.size();
  controls.rest_areas.resize(count);
  controls.area_source = 0;
  if (!settings.tick_areas.empty()) {
    controls.area_source = 3;
    std::size_t rows = settings.tick_areas.size() / count;
    const double* row = settings.tick_areas.data() + std::min(tick, rows - 1) * count;
    std::copy(row, row + count, controls.rest_areas.begin());
  } else if (!settings.key_times.empty()) {
    controls.area_source = 1;
    keyframe_share(settings.key_times, now, from, to, share);
    const double* start = settings.key_areas.data() + from * count;
    const double* end = settings.key_areas.data() + to * count;
    for (std::size_t node = 0; node < count; ++node) {
      controls.rest_areas[node] = start[node] + share * (end[node] - start[node]);
    }
  } else if (settings.glide) {
    controls.area_source = 2;
    double glide_share = std::clamp((step / kStepRate - settings.glide_start) / settings.glide_duration, 0.0, 1.0);
    for (std::size_t node = 0; node < count; ++node) {
      controls.rest_areas[node] = settings.areas[node] + glide_share * (settings.areas_to[node] - settings.areas[node]);
    }
  }
}

// The inputs of one step before the body steps (the command line's loop
// head, unchanged in order and arithmetic).
void SustainRun::control(const SustainSettings& settings, long step) {
  bool tick = step % 96 == 0;
  if (tick) {
    tick_controls(settings, step, total_, tick_);
    body_.set_pitch(tick_.pitch);
  }
  if (step == release_step_) {
    body_.set_lung_pressure(0.0, settings.release);
    body_.cord().set_target(breathy_posture(), 0.5 * settings.release);
  }
  for (std::size_t key = 0; key < settings.fold_times.size(); ++key) {
    if (trigger_step(settings.fold_times[key]) == step) {
      body_.cord().set_target(settings.fold_postures[key], settings.fold_moves[key]);
    }
  }
  if (!tick) {
    return;
  }
  if (tick_.fold_key) {
    body_.cord().set_target(tick_.fold_posture, 0.001);
  }
  if (tick_.jet) {
    body_.set_jet_gate(tick_.jet_gate);
  }
  if (tick_.velum) {
    body_.set_velum(tick_.velum_area);
  }
  if (tick_.area_source != 0) {
    std::size_t count = tick_.rest_areas.size();
    for (std::size_t node = 0; node < count && static_cast<Index>(node) < body_.vocal().count(); ++node) {
      body_.vocal().set_rest_area(static_cast<Index>(node), tick_.rest_areas[node]);
    }
  }
}

namespace {

bool same_bits(double first, double second) { return std::memcmp(&first, &second, sizeof(double)) == 0; }

bool same_bits(const std::vector<double>& first, const std::vector<double>& second) {
  return first.size() == second.size() &&
         (first.empty() || std::memcmp(first.data(), second.data(), first.size() * sizeof(double)) == 0);
}

bool same_bits(const std::optional<double>& first, const std::optional<double>& second) {
  return first.has_value() == second.has_value() && (!first.has_value() || same_bits(*first, *second));
}

bool same_bits(const CordPosture& first, const CordPosture& second) {
  return same_bits(first.lower_rest, second.lower_rest) && same_bits(first.upper_rest, second.upper_rest) &&
         same_bits(first.chink, second.chink);
}

bool same_ticks(const TickControls& first, const TickControls& second) {
  return same_bits(first.pitch, second.pitch) && first.fold_key == second.fold_key &&
         (!first.fold_key || same_bits(first.fold_posture, second.fold_posture)) && first.jet == second.jet &&
         (!first.jet || same_bits(first.jet_gate, second.jet_gate)) && first.velum == second.velum &&
         (!first.velum || same_bits(first.velum_area, second.velum_area)) && first.area_source == second.area_source &&
         (first.area_source == 0 || same_bits(first.rest_areas, second.rest_areas));
}

// What a step's forced wave does: whether the folds are forced, the wave,
// and whether RAW reports the forced flow.
struct ForcedStep {
  bool forced = false;
  double wave = 0.0;
  bool reported = false;
};

ForcedStep forced_step(const SustainSettings& settings, long step) {
  ForcedStep out;
  if (settings.forced_waves.empty()) {
    return out;
  }
  long count = static_cast<long>(settings.forced_waves.size());
  double wave = step < count ? settings.forced_waves[static_cast<std::size_t>(step)] : 0.0;
  out.forced = std::isfinite(wave);
  out.wave = out.forced ? wave : 0.0;
  out.reported = step < count && out.forced;
  return out;
}

// What a step's observed flow does (a NaN, any NaN, does nothing).
struct ObservedStep {
  bool observing = false;
  bool finite = false;
  double flow = 0.0;
};

ObservedStep observed_step(const SustainSettings& settings, long step) {
  ObservedStep out;
  out.observing = !settings.observed_flows.empty() && settings.observe_gain > 0.0;
  if (!out.observing) {
    return out;
  }
  long count = static_cast<long>(settings.observed_flows.size());
  double flow = step < count ? settings.observed_flows[static_cast<std::size_t>(step)] : std::nan("");
  out.finite = std::isfinite(flow);
  out.flow = out.finite ? flow : 0.0;
  return out;
}

}  // namespace

bool same_construction(const SustainSettings& a, const SustainSettings& b) {
  return same_bits(a.section_length, b.section_length) && same_bits(a.areas, b.areas) && same_bits(a.lung, b.lung) &&
         a.branches == b.branches && same_bits(a.aspiration, b.aspiration) && same_bits(a.frication, b.frication) &&
         same_bits(a.flux, b.flux) && same_bits(a.closing, b.closing) && same_bits(a.tremor, b.tremor) &&
         same_bits(a.adductor_tremor, b.adductor_tremor) && same_bits(a.breath_tremor, b.breath_tremor) &&
         same_bits(a.pitch_feedback, b.pitch_feedback) && same_bits(a.port_resistance, b.port_resistance) &&
         same_bits(a.reversible_jet, b.reversible_jet) &&
         same_bits(a.reversible_jet_closing, b.reversible_jet_closing) &&
         same_bits(a.fold_turbulence, b.fold_turbulence) && same_bits(a.skin, b.skin) &&
         same_bits(a.wall_scale, b.wall_scale) && same_bits(a.mass_scale, b.mass_scale) &&
         same_bits(a.axial_scale, b.axial_scale) && same_bits(a.wall_variation, b.wall_variation) &&
         same_bits(a.viscous, b.viscous) && same_bits(a.active_length, b.active_length) &&
         same_bits(a.static_area, b.static_area) && same_bits(a.reversible_nose, b.reversible_nose) &&
         same_bits(a.seal, b.seal) && same_bits(a.seal_blend, b.seal_blend) && same_bits(a.jet_bound, b.jet_bound) && a.nose_radiation == b.nose_radiation &&
         (!a.nose_radiation || (same_bits(a.nose_radiation_mass, b.nose_radiation_mass) &&
                                same_bits(a.nose_radiation_resistance, b.nose_radiation_resistance))) &&
         same_bits(a.nose_viscous, b.nose_viscous) && same_bits(a.nose_shelf, b.nose_shelf) &&
         same_bits(a.nose_hop, b.nose_hop) && same_bits(a.hop_loss, b.hop_loss) && a.rest == b.rest &&
         (!a.rest || (same_bits(a.rest_lower, b.rest_lower) && same_bits(a.rest_upper, b.rest_upper))) &&
         same_bits(a.chink, b.chink) && a.fold_damping == b.fold_damping &&
         (!a.fold_damping ||
          (same_bits(a.fold_damping_lower, b.fold_damping_lower) && same_bits(a.fold_damping_upper, b.fold_damping_upper))) &&
         same_bits(a.velum, b.velum) && a.raw_cord == b.raw_cord && a.raw_double == b.raw_double &&
         a.wave_path.empty() == b.wave_path.empty() && a.source_path.empty() == b.source_path.empty() &&
         a.raw_path == b.raw_path && a.port_waves_path == b.port_waves_path && a.port_dump_path == b.port_dump_path &&
         a.geometry_path == b.geometry_path && a.areas_path.empty() && b.areas_path.empty() &&
         a.wall_dump_path.empty() && b.wall_dump_path.empty() && a.nose_dump_path.empty() && b.nose_dump_path.empty();
}

long first_divergent_step(const SustainSettings& a, const SustainSettings& b) {
  if (!same_construction(a, b)) {
    return -1;
  }
  long total_a = total_steps(a);
  long total_b = total_steps(b);
  long limit = std::min(total_a, total_b);
  // RELEASE: its step and its glide.
  long release_a = release_step(a);
  long release_b = release_step(b);
  if (release_a != release_b || !same_bits(a.release, b.release)) {
    long earliest = std::min(release_a, release_b);
    if (earliest < limit) {
      limit = std::max(earliest, 0L);
    }
  }
  // FOLDSEQ triggers: the steps they fire at and what each sets, in order.
  for (std::size_t key = 0; key < std::max(a.fold_times.size(), b.fold_times.size()); ++key) {
    bool in_a = key < a.fold_times.size();
    bool in_b = key < b.fold_times.size();
    if (in_a && in_b && trigger_step(a.fold_times[key]) == trigger_step(b.fold_times[key]) &&
        same_bits(a.fold_postures[key], b.fold_postures[key]) && same_bits(a.fold_moves[key], b.fold_moves[key])) {
      continue;
    }
    // A trigger that differs: the earlier of its two firings (or the one).
    long earliest = limit;
    if (in_a) {
      earliest = std::min(earliest, trigger_step(a.fold_times[key]));
    }
    if (in_b) {
      earliest = std::min(earliest, trigger_step(b.fold_times[key]));
    }
    limit = std::max(std::min(limit, earliest), 0L);
  }
  // The forced waves and observed flows, step by step.
  for (long step = 0; step < limit; ++step) {
    ForcedStep forced_a = forced_step(a, step);
    ForcedStep forced_b = forced_step(b, step);
    ObservedStep observed_a = observed_step(a, step);
    ObservedStep observed_b = observed_step(b, step);
    if (forced_a.forced != forced_b.forced || !same_bits(forced_a.wave, forced_b.wave) ||
        forced_a.reported != forced_b.reported || observed_a.observing != observed_b.observing ||
        (observed_a.observing && (!same_bits(a.observe_gain, b.observe_gain) || observed_a.finite != observed_b.finite ||
                                  !same_bits(observed_a.flow, observed_b.flow)))) {
      limit = step;
      break;
    }
  }
  // The ticks' controls.
  TickControls tick_a;
  TickControls tick_b;
  for (long step = 0; step < limit; step += 96) {
    tick_controls(a, step, total_a, tick_a);
    tick_controls(b, step, total_b, tick_b);
    if (!same_ticks(tick_a, tick_b)) {
      return step;
    }
  }
  return limit;
}

// AREAS, WALLDUMP and NOSEDUMP: text files after the first step.
void SustainRun::first_step_dumps(const SustainSettings& settings) {
  // AREAS=path: the vocal tube's node areas after the first step (cm^2).
  if (!settings.areas_path.empty()) {
    std::FILE* areas_out = std::fopen(settings.areas_path.c_str(), "w");
    if (areas_out != nullptr) {
      for (Index node = 0; node < body_.vocal().count(); ++node) {
        std::fprintf(areas_out, "%.9g\n", body_.vocal().area(node));
      }
      std::fclose(areas_out);
    }
  }
  // WALLDUMP=path: each node's lining after the first step.
  if (!settings.wall_dump_path.empty()) {
    std::FILE* walls_out = std::fopen(settings.wall_dump_path.c_str(), "w");
    if (walls_out != nullptr) {
      for (Index node = 0; node < body_.vocal().count(); ++node) {
        double values[9] = {};
        body_.vocal().lining(node, values[0], values[1], values[2], values[3], values[4], values[5], values[6],
                             values[7], values[8]);
        for (double value : values) {
          std::fprintf(walls_out, "%.12g ", value);
        }
        std::fprintf(walls_out, "\n");
      }
      std::fclose(walls_out);
    }
  }
  // NOSEDUMP=path: the nose's port node, node count, hop low-pass and
  // viscous scale on the first line, then per node its area, its link's
  // viscous resistance (0 for the last) and its lining.
  if (!settings.nose_dump_path.empty()) {
    std::FILE* nose_out = std::fopen(settings.nose_dump_path.c_str(), "w");
    if (nose_out != nullptr) {
      Tube& nose = body_.nose();
      std::fprintf(nose_out, "%ld %ld %.17g %.17g\n", static_cast<long>(body_.nose_port_node()),
                   static_cast<long>(nose.count()), nose.wave_smoothing(), nose.viscous_scale());
      for (Index node = 0; node < nose.count(); ++node) {
        double row[12] = {};
        row[0] = nose.area(node);
        row[2] = node + 1 < nose.count() ? nose.link_resistance(node) : 0.0;
        nose.lining(node, row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11]);
        for (int column = 0; column < 12; ++column) {
          std::fprintf(nose_out, "%.17g ", row[column]);
        }
        std::fprintf(nose_out, "\n");
      }
      std::fclose(nose_out);
    }
  }
}

bool SustainRun::advance(const SustainSettings& settings, long end, SustainSink& sink) {
  // The run's length and release are the settings' it is advanced with (a
  // resumed run may be longer or shorter than the one it was taken from).
  total_ = total_steps(settings);
  release_step_ = release_step(settings);
  long stop = std::min(end, total_);
  const std::vector<double>& forced_waves = settings.forced_waves;
  const std::vector<double>& observed_flows = settings.observed_flows;
  long forced_count = static_cast<long>(forced_waves.size());
  long observed_count = static_cast<long>(observed_flows.size());
  bool observing = !observed_flows.empty() && settings.observe_gain > 0.0;
  while (!failed_ && step_ < stop) {
    long step = step_;
    control(settings, step);
    // FORCEWAVE: the folds driven by these returning waves (a NaN wave:
    // the folds read their own tube that step, running free).
    if (!forced_waves.empty()) {
      double wave = step < forced_count ? forced_waves[static_cast<std::size_t>(step)] : 0.0;
      body_.force_return(std::isfinite(wave), std::isfinite(wave) ? wave : 0.0);
    }
    if (observing) {
      double flow = step < observed_count ? observed_flows[static_cast<std::size_t>(step)] : std::nan("");
      body_.observe_flow(settings.observe_gain, flow);
    }
    double sound = body_.step();
    if (step == 0) {
      first_step_dumps(settings);
    }
    StepFlows flows;
    flows.entering = settings.raw_cord ? body_.cord().flow() : body_.vocal().left_flow();
    if (!forced_waves.empty() && step < forced_count && std::isfinite(forced_waves[static_cast<std::size_t>(step)])) {
      flows.entering = body_.forced_flow();
    }
    flows.radiated = body_.vocal().radiated_flow() + body_.nose().radiated_flow();
    sink.record(step, body_, flows);
    ++step_;
    double filtered = run_section(second_, run_section(first_, sound));
    if (!std::isfinite(filtered)) {
      std::fprintf(stderr, "non-finite output at %.4f s\n", step / kStepRate);
      failed_ = true;
      return false;
    }
    if (keep_source_) {
      double flow = body_.cord().flow();
      double flow_rate = run_section(source_second_, run_section(source_first_, (flow - previous_flow_) * kStepRate));
      previous_flow_ = flow;
      if (step % 2 == 0) {
        source_.push_back(flow_rate);
        source_peak_ = std::fmax(source_peak_, std::fabs(flow_rate));
      }
    }
    if (step % 2 == 0) {
      peak_ = std::fmax(peak_, std::fabs(filtered));
      if (keep_audio_) {
        audio_.push_back(filtered);
      }
    }
    sink.checked(step, body_);
  }
  return !failed_;
}

void SustainRun::finish(const SustainSettings& settings) const {
  if (failed_) {
    return;
  }
  if (!settings.wave_path.empty()) {
    std::vector<float> samples;
    samples.reserve(audio_.size());
    for (double value : audio_) {
      samples.push_back(static_cast<float>(peak_ > 0.0 ? 0.7 * value / peak_ : 0.0));
    }
    write_wave_file(settings.wave_path.c_str(), samples, static_cast<int>(kAudioRate));
  }
  if (!settings.source_path.empty()) {
    std::vector<float> source_samples;
    for (double value : source_) {
      source_samples.push_back(static_cast<float>(source_peak_ > 0.0 ? 0.7 * value / source_peak_ : 0.0));
    }
    write_wave_file(settings.source_path.c_str(), source_samples, static_cast<int>(kAudioRate));
  }
}

}  // namespace tract
