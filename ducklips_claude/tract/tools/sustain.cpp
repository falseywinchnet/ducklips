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
  if (!read_profile(values[1], areas, section_length)) {
    std::fprintf(stderr, "no profile %s in data/story_2008_areas.txt\n", values[1]);
    return 1;
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
  for (tract::SegmentMaterial& material : layout.vocal_material_list) {
    material.axial_tension *= wall_scale;
    material.hoop_tension *= wall_scale;
    material.coupling_stiffness *= wall_scale;
    material.support_stiffness *= wall_scale;
    material.cover_thickness *= mass_scale;
    material.body_depth *= mass_scale;
  }
  tract::VoiceBody body(layout);
  body.set_pitch(pitch);
  tract::CordPosture posture = tract::modal_posture();
  if (std::getenv("REST") != nullptr) {
    // REST=lower,upper (cm): the folds' rest displacements at the vocal processes.
    posture.lower_rest = std::strtod(std::getenv("REST"), nullptr);
    posture.upper_rest = std::strtod(std::strchr(std::getenv("REST"), ',') + 1, nullptr);
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
    double sound = body.step();
    if (raw_out.is_open()) {
      float pair[2] = {static_cast<float>(std::getenv("RAWCORD") != nullptr ? body.cord().flow()
                                                                            : body.vocal().left_flow()),
                       static_cast<float>(body.vocal().radiated_flow() + body.nose().radiated_flow())};
      raw_out.write(reinterpret_cast<const char*>(pair), sizeof(pair));
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
