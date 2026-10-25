// Sustain one vowel through the segment-module tract.
// usage: sustain VOWEL OUT.wav [pitch_hz] [lung_pa] [seconds]
//   VOWEL is a row of data/story_2008_areas.txt (iy ih ey eh ae ah aa ao ow uh uw).
// Environment: TRACE=ms prints the cord and pressures every ms milliseconds;
// SKIN=gain sets how much lining motion is heard through the skin;
// PITCHEND=hz glides the intended pitch to that value by the end; RELEASE=s
// ends the vowel as a speaker does: breath let go and folds parted over the
// last s seconds.
// TREMOR=share (rms of cord tension) and FOLDTURB=share (jet pressure on
// the folds) set the fast causes; SOURCE=path also writes the glottal flow derivative (the source alone).

#include "tract/voice_body.hpp"
#include "wave_out.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

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

}  // namespace

int main(int count, char** values) {
  if (count < 3) {
    std::fprintf(stderr, "usage: sustain VOWEL OUT.wav [pitch_hz] [lung_pa] [seconds]\n");
    return 2;
  }
  std::vector<double> areas;
  double section_length = 0.0;
  // AREAFILE=path: the tube's profile instead of a Story vowel: the section
  // length (cm) on the first line, then one area (cm^2) per line.
  if (std::getenv("AREAFILE") != nullptr) {
    std::ifstream area_file(std::getenv("AREAFILE"));
    double area = 0.0;
    area_file >> section_length;
    while (area_file >> area) {
      areas.push_back(area);
    }
    if (areas.empty() || section_length <= 0.0) {
      std::fprintf(stderr, "cannot read AREAFILE %s\n", std::getenv("AREAFILE"));
      return 1;
    }
  } else if (!read_profile(values[1], areas, section_length)) {
    std::fprintf(stderr, "no profile %s in data/story_2008_areas.txt\n", values[1]);
    return 1;
  }
  // AREAFILE2=path with GLIDE=start,duration (s): the tube glides from
  // AREAFILE's areas to these (the same count), linearly at each 1 ms
  // control tick, starting at `start`.
  std::vector<double> areas_to;
  double glide_start = 0.0;
  double glide_duration = 0.0;
  if (std::getenv("AREAFILE2") != nullptr && std::getenv("GLIDE") != nullptr) {
    std::ifstream area_file(std::getenv("AREAFILE2"));
    double second_length = 0.0;
    double area = 0.0;
    area_file >> second_length;
    while (area_file >> area) {
      areas_to.push_back(area);
    }
    glide_start = std::strtod(std::getenv("GLIDE"), nullptr);
    glide_duration = std::strtod(std::strchr(std::getenv("GLIDE"), ',') + 1, nullptr);
    if (areas_to.size() != areas.size() || glide_duration <= 0.0) {
      std::fprintf(stderr, "AREAFILE2 needs %zu areas and GLIDE=start,duration\n", areas.size());
      return 1;
    }
  }
  // AREASEQ=path: keyframes, one per line: a time (s) then the node areas
  // (AREAFILE's count); at each 1 ms tick the rest areas are interpolated
  // linearly between the keyframes around it (held before the first and
  // after the last). Overrides AREAFILE2 and GLIDE.
  std::vector<double> key_times;
  std::vector<std::vector<double>> key_areas;
  if (std::getenv("AREASEQ") != nullptr) {
    std::ifstream sequence_file(std::getenv("AREASEQ"));
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
      if (frame.size() != areas.size()) {
        std::fprintf(stderr, "AREASEQ keyframe at %g s has %zu areas, not %zu\n", time, frame.size(), areas.size());
        return 1;
      }
      key_times.push_back(time);
      key_areas.push_back(frame);
    }
  }
  // VELUMSEQ=path: the velum's port area (cm^2) over time, one "time area"
  // per line, linear between, set at each 1 ms tick.
  std::vector<double> velum_times;
  std::vector<double> velum_areas;
  if (std::getenv("VELUMSEQ") != nullptr) {
    std::ifstream velum_file(std::getenv("VELUMSEQ"));
    double time = 0.0;
    double area = 0.0;
    while (velum_file >> time >> area) {
      velum_times.push_back(time);
      velum_areas.push_back(area);
    }
  }
  // FOLDSEQ=path: the folds' posture over time, one "time lower_rest
  // upper_rest chink move_seconds" per line (cm, cm, cm^2, s): at that time
  // the folds start toward that posture (Cord::set_target), as an
  // abduction for a voiceless phone and the adduction after it.
  std::vector<long> fold_steps;
  std::vector<tract::CordPosture> fold_postures;
  std::vector<double> fold_moves;
  if (std::getenv("FOLDSEQ") != nullptr) {
    std::ifstream fold_file(std::getenv("FOLDSEQ"));
    double time = 0.0;
    tract::CordPosture fold_posture;
    double move = 0.0;
    while (fold_file >> time >> fold_posture.lower_rest >> fold_posture.upper_rest >> fold_posture.chink >> move) {
      fold_steps.push_back(std::lround(time * tract::kStepRate));
      fold_postures.push_back(fold_posture);
      fold_moves.push_back(move);
    }
  }
  // JETSEQ=path: the reversible jet's gate over time, one "time gate" per
  // line, linear between, set at each 1 ms tick (VoiceBody::set_jet_gate).
  std::vector<double> jet_times;
  std::vector<double> jet_gates;
  if (std::getenv("JETSEQ") != nullptr) {
    std::ifstream jet_file(std::getenv("JETSEQ"));
    double time = 0.0;
    double gate = 0.0;
    while (jet_file >> time >> gate) {
      jet_times.push_back(time);
      jet_gates.push_back(gate);
    }
  }
  // FOLDKEYS=path: the folds' posture as keyframes, one "time lower_rest
  // upper_rest chink" per line, linear between and set at each 1 ms tick
  // (Cord::set_target over 1 ms): unlike FOLDSEQ's triggers, smooth in the
  // keyframes' times (a trigger moved by less than a step did nothing,
  // then jumped, and a search's time derivatives were spikes).
  std::vector<double> fold_key_times;
  std::vector<tract::CordPosture> fold_key_postures;
  if (std::getenv("FOLDKEYS") != nullptr) {
    std::ifstream key_file(std::getenv("FOLDKEYS"));
    double time = 0.0;
    tract::CordPosture key_posture;
    while (key_file >> time >> key_posture.lower_rest >> key_posture.upper_rest >> key_posture.chink) {
      fold_key_times.push_back(time);
      fold_key_postures.push_back(key_posture);
    }
  }
  double pitch = count > 3 ? std::strtod(values[3], nullptr) : 120.0;
  double lung = count > 4 ? std::strtod(values[4], nullptr) : 800.0;
  double seconds = count > 5 ? std::strtod(values[5], nullptr) : 1.5;
  double pitch_end = std::getenv("PITCHEND") != nullptr ? std::strtod(std::getenv("PITCHEND"), nullptr) : pitch;
  double release = std::getenv("RELEASE") != nullptr ? std::strtod(std::getenv("RELEASE"), nullptr) : 0.0;
  double trace_ms = std::getenv("TRACE") != nullptr ? std::strtod(std::getenv("TRACE"), nullptr) : 0.0;

  tract::BodyLayout layout = tract::vowel_layout(areas, section_length);
  if (std::getenv("BRANCHES") != nullptr) {
    // BRANCHES=pus: keep piriform (p), under-tongue (u), sinuses (s).
    std::string keep = std::getenv("BRANCHES");
    layout.piriform = keep.find('p') != std::string::npos;
    layout.under_tongue = keep.find('u') != std::string::npos;
    layout.sinuses = keep.find('s') != std::string::npos;
  }
  if (std::getenv("ASPIRATION") != nullptr) {
    layout.aspiration = std::strtod(std::getenv("ASPIRATION"), nullptr);
  }
  if (std::getenv("FRICATION") != nullptr) {
    layout.frication = std::strtod(std::getenv("FRICATION"), nullptr);
  }
  if (std::getenv("FLUX") != nullptr) {
    layout.jet_flux = std::strtod(std::getenv("FLUX"), nullptr);
  }
  if (std::getenv("CLOSING") != nullptr) {
    layout.closing_turbulence = std::strtod(std::getenv("CLOSING"), nullptr);
  }
  if (std::getenv("TREMOR") != nullptr) {
    layout.tremor = std::strtod(std::getenv("TREMOR"), nullptr);
  }
  if (std::getenv("ADDTREMOR") != nullptr) {
    layout.adductor_tremor = std::strtod(std::getenv("ADDTREMOR"), nullptr);
  }
  if (std::getenv("BREATHTREMOR") != nullptr) {
    layout.breath_tremor = std::strtod(std::getenv("BREATHTREMOR"), nullptr);
  }
  if (std::getenv("PITCHFB") != nullptr) {
    layout.pitch_feedback = std::strtod(std::getenv("PITCHFB"), nullptr);
  }
  // PORTRESIST=s: a factor on the velum port's viscous resistance.
  if (std::getenv("PORTRESIST") != nullptr) {
    layout.port_resistance_scale = std::strtod(std::getenv("PORTRESIST"), nullptr);
  }
  // REVJET=gain: the reversible jet (VoiceBody::reversible_jet).
  if (std::getenv("REVJET") != nullptr) {
    layout.reversible_jet = std::strtod(std::getenv("REVJET"), nullptr);
  }
  // REVJETCLOSE=cm^2: the reversible jet fades out as its narrowest node
  // closes (none at a quarter of it, full from it).
  if (std::getenv("REVJETCLOSE") != nullptr) {
    layout.reversible_jet_closing = std::strtod(std::getenv("REVJETCLOSE"), nullptr);
  }
  if (std::getenv("FOLDTURB") != nullptr) {
    layout.fold_turbulence = std::strtod(std::getenv("FOLDTURB"), nullptr);
  }
  if (std::getenv("SKIN") != nullptr) {
    layout.skin_gain = std::strtod(std::getenv("SKIN"), nullptr);
  }
  // WALLS=s scales every wall stiffness and tension (diagnostic: large s makes
  // the walls rigid); WALLMASS=m scales the wall masses.
  double wall_scale = std::getenv("WALLS") != nullptr ? std::strtod(std::getenv("WALLS"), nullptr) : 1.0;
  double mass_scale = std::getenv("WALLMASS") != nullptr ? std::strtod(std::getenv("WALLMASS"), nullptr) : 1.0;
  // AXIAL=s scales only the walls' axial tension (0: each node's lining moves
  // on its own, which lets the tube be driven backward node by node).
  double axial_scale = std::getenv("AXIAL") != nullptr ? std::strtod(std::getenv("AXIAL"), nullptr) : 1.0;
  for (tract::SegmentMaterial& material : layout.vocal_material_list) {
    material.axial_tension *= axial_scale;
    material.axial_tension *= wall_scale;
    material.hoop_tension *= wall_scale;
    material.coupling_stiffness *= wall_scale;
    material.support_stiffness *= wall_scale;
    material.cover_thickness *= mass_scale;
    material.body_depth *= mass_scale;
  }
  // WALLVAR=share: the walls' support variation (0: steady walls).
  if (std::getenv("WALLVAR") != nullptr) {
    layout.wall_variation = std::strtod(std::getenv("WALLVAR"), nullptr);
  }
  // VISCOUS=s scales the links' viscous resistance (0: none).
  if (std::getenv("VISCOUS") != nullptr) {
    layout.viscous_scale = std::strtod(std::getenv("VISCOUS"), nullptr);
  }
  tract::VoiceBody body(layout);
  // ACTIVELENGTH=cm: the vocal tube's length, whole nodes and the rest as
  // the open end's air mass (Tube::set_active_length), so a tube's length
  // moves continuously.
  if (std::getenv("ACTIVELENGTH") != nullptr) {
    body.vocal().set_active_length(std::strtod(std::getenv("ACTIVELENGTH"), nullptr));
  }
  // STATICAREA=1: the air sees the rest geometry only (Tube::set_static_admittance).
  if (std::getenv("STATICAREA") != nullptr) {
    body.vocal().set_static_admittance(std::strtod(std::getenv("STATICAREA"), nullptr) != 0.0);
  }
  // REVNOSE=1: the nose as the reversible engine needs it: linings without
  // axial tension, the air on the rest geometry.
  if (std::getenv("REVNOSE") != nullptr && std::strtod(std::getenv("REVNOSE"), nullptr) != 0.0) {
    body.nose().scale_axial(0.0);
    body.nose().set_static_admittance(true);
  }
  // SEAL=area: the vocal tube's nodes at or under this area seal (Tube::set_seal_area).
  if (std::getenv("SEAL") != nullptr) {
    body.vocal().set_seal_area(std::strtod(std::getenv("SEAL"), nullptr));
  }
  // NOSERAD=mass,resistance: scales on the nostrils' radiation load (default 1,1).
  if (std::getenv("NOSERAD") != nullptr) {
    body.nose().set_radiation_scales(std::strtod(std::getenv("NOSERAD"), nullptr),
                                     std::strtod(std::strchr(std::getenv("NOSERAD"), ',') + 1, nullptr));
  }
  // NOSEVISCOUS=s: the nose's viscous resistance scale (default 1).
  if (std::getenv("NOSEVISCOUS") != nullptr) {
    body.nose().set_viscous_scale(std::strtod(std::getenv("NOSEVISCOUS"), nullptr));
  }
  // NOSESHELF=b: the nose's hop filter as a shelf (Tube::set_wave_shelf; with NOSEHOP its pole).
  if (std::getenv("NOSESHELF") != nullptr) {
    body.nose().set_wave_shelf(std::strtod(std::getenv("NOSESHELF"), nullptr));
  }
  // NOSEHOP=k: the nose's per-hop low-pass (its default 0: lossless hops).
  if (std::getenv("NOSEHOP") != nullptr) {
    body.nose().set_wave_smoothing(std::strtod(std::getenv("NOSEHOP"), nullptr));
  }
  // HOPLOSS=k: the vocal tube's per-hop low-pass (0: lossless hops).
  if (std::getenv("HOPLOSS") != nullptr) {
    body.vocal().set_wave_smoothing(std::strtod(std::getenv("HOPLOSS"), nullptr));
  }
  body.set_pitch(pitch);
  tract::CordPosture posture = tract::modal_posture();
  if (std::getenv("REST") != nullptr) {
    // REST=lower,upper (cm): the folds' rest displacements at the vocal processes.
    posture.lower_rest = std::strtod(std::getenv("REST"), nullptr);
    posture.upper_rest = std::strtod(std::strchr(std::getenv("REST"), ',') + 1, nullptr);
  }
  // CHINK=cm^2: the folds' always-open chink (0: a full closure).
  if (std::getenv("CHINK") != nullptr) {
    posture.chink = std::strtod(std::getenv("CHINK"), nullptr);
  }
  body.cord().set_target(posture, 0.001);
  if (std::getenv("FOLDDAMP") != nullptr) {
    // FOLDDAMP=lower,upper: damping ratios of the fold edges.
    body.cord().set_damping(std::strtod(std::getenv("FOLDDAMP"), nullptr),
                            std::strtod(std::strchr(std::getenv("FOLDDAMP"), ',') + 1, nullptr));
  }
  body.set_lung_pressure(lung, 0.08);
  if (std::getenv("VELUM") != nullptr) {
    body.set_velum(std::strtod(std::getenv("VELUM"), nullptr));
  }

  std::vector<float> samples;
  Biquad first = lowpass(18000.0, tract::kStepRate, 0.5412);
  Biquad second = lowpass(18000.0, tract::kStepRate, 1.3066);
  long total = static_cast<long>(seconds * tract::kStepRate);
  long trace_every = trace_ms > 0.0 ? static_cast<long>(trace_ms * 1e-3 * tract::kStepRate) : 0;
  double peak = 0.0;
  std::vector<double> raw;
  raw.reserve(static_cast<std::size_t>(total / 2 + 1));
  const char* source_path = std::getenv("SOURCE");
  Biquad source_first = lowpass(18000.0, tract::kStepRate, 0.5412);
  Biquad source_second = lowpass(18000.0, tract::kStepRate, 1.3066);
  std::vector<double> source;
  double source_peak = 0.0;
  double previous_flow = 0.0;
  // RAW=path: per step, the glottal flow and the radiated flows (lips + nose
  // + skin as heard), float32 pairs at the step rate, for transfer estimates.
  std::ofstream raw_out;
  if (std::getenv("RAW") != nullptr) {
    raw_out.open(std::getenv("RAW"), std::ios::binary);
  }
  std::vector<double> forced_waves;
  if (std::getenv("FORCEWAVE") != nullptr) {
    std::ifstream wave_file(std::getenv("FORCEWAVE"), std::ios::binary);
    double wave = 0.0;
    while (wave_file.read(reinterpret_cast<char*>(&wave), sizeof(wave))) {
      forced_waves.push_back(wave);
    }
  }
  // FORCEFLOW=path (float64 per step, NaN: none) with OBSERVE=gain: the
  // forced folds observed against these flows (VoiceBody::observe_flow).
  std::vector<double> observed_flows;
  double observe_gain = std::getenv("OBSERVE") != nullptr ? std::strtod(std::getenv("OBSERVE"), nullptr) : 0.0;
  if (std::getenv("FORCEFLOW") != nullptr) {
    std::ifstream flow_file(std::getenv("FORCEFLOW"), std::ios::binary);
    double flow = 0.0;
    while (flow_file.read(reinterpret_cast<char*>(&flow), sizeof(flow))) {
      observed_flows.push_back(flow);
    }
  }
  bool raw_double = std::getenv("RAWDOUBLE") != nullptr && std::strtod(std::getenv("RAWDOUBLE"), nullptr) != 0.0;
  std::ofstream port_waves_out;
  if (std::getenv("PORTWAVES") != nullptr) {
    port_waves_out.open(std::getenv("PORTWAVES"), std::ios::binary);
  }
  std::ofstream port_out;
  if (std::getenv("PORTDUMP") != nullptr) {
    port_out.open(std::getenv("PORTDUMP"), std::ios::binary);
  }
  std::ofstream geometry_out;
  if (std::getenv("GEODUMP") != nullptr) {
    geometry_out.open(std::getenv("GEODUMP"), std::ios::binary);
  }
  long release_step = release > 0.0 ? total - static_cast<long>(release * tract::kStepRate) : total;
  for (long step = 0; step < total; ++step) {
    if (step % 96 == 0) {
      double progress = static_cast<double>(step) / static_cast<double>(total);
      body.set_pitch(pitch + (pitch_end - pitch) * progress);
    }
    if (step == release_step) {
      body.set_lung_pressure(0.0, release);
      body.cord().set_target(tract::breathy_posture(), 0.5 * release);
    }
    for (std::size_t key = 0; key < fold_steps.size(); ++key) {
      if (fold_steps[key] == step) {
        body.cord().set_target(fold_postures[key], fold_moves[key]);
      }
    }
    if (!fold_key_times.empty() && step % 96 == 0) {
      double now = step / tract::kStepRate;
      std::size_t after = 0;
      while (after < fold_key_times.size() && fold_key_times[after] <= now) {
        ++after;
      }
      std::size_t from = after == 0 ? 0 : after - 1;
      std::size_t to = after < fold_key_times.size() ? after : fold_key_times.size() - 1;
      double share = to == from ? 0.0
                                : std::clamp((now - fold_key_times[from]) / (fold_key_times[to] - fold_key_times[from]),
                                             0.0, 1.0);
      const tract::CordPosture& start = fold_key_postures[from];
      const tract::CordPosture& end = fold_key_postures[to];
      tract::CordPosture blended;
      blended.lower_rest = start.lower_rest + share * (end.lower_rest - start.lower_rest);
      blended.upper_rest = start.upper_rest + share * (end.upper_rest - start.upper_rest);
      blended.chink = start.chink + share * (end.chink - start.chink);
      body.cord().set_target(blended, 0.001);
    }
    if (!jet_times.empty() && step % 96 == 0) {
      double now = step / tract::kStepRate;
      std::size_t after = 0;
      while (after < jet_times.size() && jet_times[after] <= now) {
        ++after;
      }
      std::size_t from = after == 0 ? 0 : after - 1;
      std::size_t to = after < jet_times.size() ? after : jet_times.size() - 1;
      double share = to == from ? 0.0 : std::clamp((now - jet_times[from]) / (jet_times[to] - jet_times[from]), 0.0, 1.0);
      body.set_jet_gate(jet_gates[from] + share * (jet_gates[to] - jet_gates[from]));
    }
    if (!velum_times.empty() && step % 96 == 0) {
      double now = step / tract::kStepRate;
      std::size_t after = 0;
      while (after < velum_times.size() && velum_times[after] <= now) {
        ++after;
      }
      std::size_t from = after == 0 ? 0 : after - 1;
      std::size_t to = after < velum_times.size() ? after : velum_times.size() - 1;
      double share = to == from ? 0.0 : std::clamp((now - velum_times[from]) / (velum_times[to] - velum_times[from]), 0.0, 1.0);
      body.set_velum(velum_areas[from] + share * (velum_areas[to] - velum_areas[from]));
    }
    if (!key_times.empty() && step % 96 == 0) {
      double now = step / tract::kStepRate;
      std::size_t after = 0;
      while (after < key_times.size() && key_times[after] <= now) {
        ++after;
      }
      std::size_t from = after == 0 ? 0 : after - 1;
      std::size_t to = after < key_times.size() ? after : key_times.size() - 1;
      double share = to == from ? 0.0 : std::clamp((now - key_times[from]) / (key_times[to] - key_times[from]), 0.0, 1.0);
      for (std::size_t node = 0; node < areas.size() && static_cast<tract::Index>(node) < body.vocal().count(); ++node) {
        body.vocal().set_rest_area(static_cast<tract::Index>(node),
                                   key_areas[from][node] + share * (key_areas[to][node] - key_areas[from][node]));
      }
    } else if (!areas_to.empty() && step % 96 == 0) {
      double share = std::clamp((step / tract::kStepRate - glide_start) / glide_duration, 0.0, 1.0);
      for (std::size_t node = 0; node < areas.size() && static_cast<tract::Index>(node) < body.vocal().count(); ++node) {
        body.vocal().set_rest_area(static_cast<tract::Index>(node), areas[node] + share * (areas_to[node] - areas[node]));
      }
    }
    // FORCEWAVE=path (float64 per step): the folds driven by these returning
    // waves (VoiceBody::force_return); RAW's first column is then the flow
    // they pass (VoiceBody::forced_flow).
    // (A NaN wave: the folds read their own tube that step, running free;
    // for the stretches the sound cannot reveal.)
    if (!forced_waves.empty()) {
      double wave = step < static_cast<long>(forced_waves.size()) ? forced_waves[static_cast<std::size_t>(step)] : 0.0;
      body.force_return(std::isfinite(wave), std::isfinite(wave) ? wave : 0.0);
    }
    if (!observed_flows.empty() && observe_gain > 0.0) {
      double flow = step < static_cast<long>(observed_flows.size()) ? observed_flows[static_cast<std::size_t>(step)]
                                                                    : std::nan("");
      body.observe_flow(observe_gain, flow);
    }
    double sound = body.step();
    // GEODUMP=path: per step, each node's area, the air the frame draws
    // from it, its link's viscous resistance (0 for the last) and its
    // lining (WALLDUMP's nine), float64, for a replica of a moving tube.
    if (geometry_out.is_open()) {
      tract::Tube& vocal = body.vocal();
      for (tract::Index node = 0; node < vocal.count(); ++node) {
        double row[12] = {};
        row[0] = vocal.area(node);
        row[1] = vocal.frame_flow(node);
        row[2] = node + 1 < vocal.count() ? vocal.link_resistance(node) : 0.0;
        vocal.lining(node, row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11]);
        geometry_out.write(reinterpret_cast<const char*>(row), sizeof(row));
      }
    }
    // AREAS=path: the vocal tube's node areas after the first step (cm^2),
    // one per line, for a replica of the tube (the reversible core).
    if (step == 0 && std::getenv("AREAS") != nullptr) {
      std::FILE* areas_out = std::fopen(std::getenv("AREAS"), "w");
      if (areas_out != nullptr) {
        for (tract::Index node = 0; node < body.vocal().count(); ++node) {
          std::fprintf(areas_out, "%.9g\n", body.vocal().area(node));
        }
        std::fclose(areas_out);
      }
    }
    // WALLDUMP=path: each node's lining after the first step (surface,
    // cover and body masses, hoop, coupling and support stiffnesses, cover
    // and body damping, rest radius), one node per line.
    if (step == 0 && std::getenv("WALLDUMP") != nullptr) {
      std::FILE* walls_out = std::fopen(std::getenv("WALLDUMP"), "w");
      if (walls_out != nullptr) {
        for (tract::Index node = 0; node < body.vocal().count(); ++node) {
          double values[9] = {};
          body.vocal().lining(node, values[0], values[1], values[2], values[3], values[4], values[5], values[6],
                              values[7], values[8]);
          for (double value : values) {
            std::fprintf(walls_out, "%.12g ", value);
          }
          std::fprintf(walls_out, "\n");
        }
        std::fclose(walls_out);
      }
    }
    // PORTWAVES=path: per step (float64) the wave the port delivered into the
    // nose, the wave it sent into the vocal node, that node's pressure, the
    // wave arriving at it from its outer link next step, the nostrils' and
    // the lips' flows: the reversible engine's window at a nasal release.
    if (port_waves_out.is_open()) {
      tract::Index port_node = body.nose_port_node();
      double row[6] = {body.nose().left_incoming(), body.vocal().side_incoming(0), body.vocal().pressure(port_node),
                       body.vocal().backward_wave(port_node), body.nose().radiated_flow(), body.vocal().radiated_flow()};
      port_waves_out.write(reinterpret_cast<const char*>(row), sizeof(row));
    }
    // PORTDUMP=path: the nasal port's area each step (float64).
    if (port_out.is_open()) {
      double port = body.velum_port();
      port_out.write(reinterpret_cast<const char*>(&port), sizeof(port));
    }
    // NOSEDUMP=path: after the first step, the nose: its port node, node
    // count, hop low-pass and viscous scale on the first line, then per node
    // its area, its link's viscous resistance (0 for the last) and its
    // lining (WALLDUMP's nine).
    if (step == 0 && std::getenv("NOSEDUMP") != nullptr) {
      std::FILE* nose_out = std::fopen(std::getenv("NOSEDUMP"), "w");
      if (nose_out != nullptr) {
        tract::Tube& nose = body.nose();
        std::fprintf(nose_out, "%ld %ld %.17g %.17g\n", static_cast<long>(body.nose_port_node()),
                     static_cast<long>(nose.count()), nose.wave_smoothing(), nose.viscous_scale());
        for (tract::Index node = 0; node < nose.count(); ++node) {
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
    if (raw_out.is_open()) {
      double entering = std::getenv("RAWCORD") != nullptr ? body.cord().flow() : body.vocal().left_flow();
      if (!forced_waves.empty() && step < static_cast<long>(forced_waves.size()) &&
          std::isfinite(forced_waves[static_cast<std::size_t>(step)])) {
        entering = body.forced_flow();
      }
      double radiated = body.vocal().radiated_flow() + body.nose().radiated_flow();
      if (raw_double) {
        // RAWDOUBLE=1: float64 pairs (a search's differences of 1e-3 in a
        // log area drowned in float32's rounding near -70 dB).
        double pair[2] = {entering, radiated};
        raw_out.write(reinterpret_cast<const char*>(pair), sizeof(pair));
      } else {
        float pair[2] = {static_cast<float>(entering), static_cast<float>(radiated)};
        raw_out.write(reinterpret_cast<const char*>(pair), sizeof(pair));
      }
    }
    double filtered = run(second, run(first, sound));
    if (!std::isfinite(filtered)) {
      std::fprintf(stderr, "non-finite output at %.4f s\n", step / tract::kStepRate);
      return 1;
    }
    double flow = body.cord().flow();
    double flow_rate = run(source_second, run(source_first, (flow - previous_flow) * tract::kStepRate));
    previous_flow = flow;
    if (step % 2 == 0) {
      raw.push_back(filtered);
      peak = std::fmax(peak, std::fabs(filtered));
      if (source_path != nullptr) {
        source.push_back(flow_rate);
        source_peak = std::fmax(source_peak, std::fabs(flow_rate));
      }
    }
    if (trace_every > 0 && step % trace_every == 0) {
      tract::Tube& trachea = body.trachea();
      tract::Tube& vocal = body.vocal();
      std::printf("%.4f area %.4f flow %7.1f contact %.2f sub %7.0f supra %7.0f lips %8.1f wall %8.2f nose_p %8.1f nose %8.1f\n",
                  step / tract::kStepRate, body.cord().area(), body.cord().flow(), body.cord().contact(),
                  trachea.pressure(trachea.count() - 1), vocal.pressure(0), vocal.radiated_flow(),
                  vocal.total_wall_flow(), body.nose().pressure(0), body.nose().radiated_flow());
    }
  }
  for (double value : raw) {
    samples.push_back(static_cast<float>(peak > 0.0 ? 0.7 * value / peak : 0.0));
  }
  write_wave(values[2], samples, static_cast<int>(tract::kAudioRate));
  if (source_path != nullptr) {
    std::vector<float> source_samples;
    for (double value : source) {
      source_samples.push_back(static_cast<float>(source_peak > 0.0 ? 0.7 * value / source_peak : 0.0));
    }
    write_wave(source_path, source_samples, static_cast<int>(tract::kAudioRate));
  }
  std::fprintf(stderr, "%s: %.2f s, peak %.4g (normalized to 0.7)\n", values[2], seconds, peak);
  return 0;
}
