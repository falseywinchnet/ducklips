// Transfer function of a vocal tube alone: closed glottis end, one pressure
// impulse entering at the glottis, radiated flow recorded at the lips.
// Prints resonance peaks (frequency, level, -3 dB bandwidth) and band levels.
// usage: transfer VOWEL [rigid]
//   rigid: walls 1000x stiffer (the lining barely moves)

#include "tract/voice_body.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
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

double magnitude_db(const std::vector<double>& response, double hz) {
  std::complex<double> sum(0.0, 0.0);
  double omega = tract::kTwoPi * hz / tract::kStepRate;
  for (std::size_t n = 0; n < response.size(); ++n) {
    sum += response[n] * std::polar(1.0, -omega * static_cast<double>(n));
  }
  return 20.0 * std::log10(std::abs(sum) + 1e-30);
}

}  // namespace

int main(int count, char** values) {
  if (count < 2) {
    std::fprintf(stderr, "usage: transfer VOWEL [rigid]\n");
    return 2;
  }
  std::vector<double> areas;
  double section_length = 0.0;
  if (!read_profile(values[1], areas, section_length)) {
    std::fprintf(stderr, "no profile %s\n", values[1]);
    return 1;
  }
  bool rigid = count > 2 && std::string(values[2]) == "rigid";
  tract::BodyLayout layout = tract::vowel_layout(areas, section_length);
  std::vector<tract::TubeNode> nodes;
  for (std::size_t index = 0; index < layout.vocal_areas.size(); ++index) {
    tract::TubeNode node;
    node.rest_area = layout.vocal_areas[index];
    node.material = layout.vocal_materials[index];
    nodes.push_back(node);
  }
  if (rigid) {
    for (tract::SegmentMaterial& material : layout.vocal_material_list) {
      material.axial_tension *= 1000.0;
      material.hoop_tension *= 1000.0;
      material.coupling_stiffness *= 1000.0;
      material.support_stiffness *= 1000.0;
    }
  }
  tract::Tube tube(layout.vocal_material_list, nodes);
  tube.set_radiating_right(true);
  // VISCOUS=x: a factor on the links' viscous resistance; HOPLOSS=k: the hop low-pass.
  if (std::getenv("VISCOUS") != nullptr) {
    tube.set_viscous_scale(std::strtod(std::getenv("VISCOUS"), nullptr));
  }
  if (std::getenv("HOPLOSS") != nullptr) {
    tube.set_wave_smoothing(std::strtod(std::getenv("HOPLOSS"), nullptr));
  }
  // 200 ms of response: an impulse wave at a nearly closed glottis.
  long steps = static_cast<long>(0.2 * tract::kStepRate);
  std::vector<double> response;
  double previous = 0.0;
  // DCFLOW=w: a steady driving wave under the impulse (air flowing, as in
  // phonation); settle 100 ms first, then subtract the steady response.
  double steady = std::getenv("DCFLOW") != nullptr ? std::strtod(std::getenv("DCFLOW"), nullptr) : 0.0;
  // GLOTTIS=a: the glottal end as an absorbing port of area a (cm^2) instead
  // of the nearly closed 0.001.
  double glottis = std::getenv("GLOTTIS") != nullptr ? std::strtod(std::getenv("GLOTTIS"), nullptr) : 1.0e-3;
  for (long step = 0; step < static_cast<long>(0.1 * tract::kStepRate) && steady != 0.0; ++step) {
    tube.set_left_port(steady, glottis);
    tube.step();
  }
  // The response is the difference from an identical twin without the
  // impulse, so any slow settling of the steady flow cancels exactly.
  tract::Tube twin = tube;
  double twin_previous = twin.radiated_flow();
  previous = tube.radiated_flow();
  if (steady != 0.0) {
    std::printf("steady flow %.1f cm3/s\n", tube.radiated_flow());
  }
  // PULSE=1: the steady drive pulses at 120 Hz (half-wave sine, as glottal
  // pulses), so the losses see phonation-sized flow swings.
  bool pulse = std::getenv("PULSE") != nullptr;
  for (long step = 0; step < steps; ++step) {
    double drive = steady;
    if (pulse) {
      drive = 2.0 * steady * std::max(std::sin(tract::kTwoPi * 120.0 * step / tract::kStepRate), 0.0);
    }
    tube.set_left_port(drive + (step == 0 ? 1.0 : 0.0), glottis);
    twin.set_left_port(drive, glottis);
    tube.step();
    twin.step();
    double flow = tube.radiated_flow();
    double twin_flow = twin.radiated_flow();
    response.push_back(((flow - previous) - (twin_flow - twin_previous)) * tract::kStepRate);
    previous = flow;
    twin_previous = twin_flow;
  }
  // BODY=1: the same impulse as a series source on the first link, once in
  // the bare tube and once in the whole body at rest (no breath, folds
  // apart as at rest, branches off), to see what the body adds.
  if (std::getenv("BODY") != nullptr) {
    tract::BodyLayout quiet = tract::vowel_layout(areas, section_length);
    quiet.piriform = false;
    quiet.under_tongue = false;
    quiet.sinuses = false;
    quiet.breath_variation = 0.0;
    quiet.tremor = 0.0;
    quiet.adductor_tremor = 0.0;
    quiet.breath_tremor = 0.0;
    quiet.aspiration = 0.0;
    quiet.frication = 0.0;
    tract::VoiceBody body(quiet);
    // BODY=voice: the body phonates (/a/ at 120 Hz, 800 Pa) for 0.5 s first,
    // and the impulse response is taken against an identical twin body.
    bool voicing = std::string(std::getenv("BODY")) == "voice";
    if (voicing) {
      body.set_pitch(120.0);
      tract::CordPosture posture = tract::modal_posture();
      if (std::getenv("REST") != nullptr) {
        posture.lower_rest = std::strtod(std::getenv("REST"), nullptr);
        posture.upper_rest = std::strtod(std::strchr(std::getenv("REST"), ',') + 1, nullptr);
      }
      body.cord().set_target(posture, 0.001);
      double lung = std::getenv("LUNG") != nullptr ? std::strtod(std::getenv("LUNG"), nullptr) : 800.0;
      body.set_lung_pressure(lung, 0.08);
      for (long step = 0; step < static_cast<long>(0.5 * tract::kStepRate); ++step) {
        body.step();
      }
    }
    tract::VoiceBody twin_body = body;
    tract::Tube alone(layout.vocal_material_list, nodes);
    alone.set_radiating_right(true);
    response.clear();
    std::vector<double> body_response;
    double alone_previous = 0.0;
    double body_previous = 0.0;
    for (long step = 0; step < steps; ++step) {
      if (step == 0) {
        alone.add_series_source(0, 1.0);
        body.vocal().add_series_source(0, 1.0);
      }
      alone.step();
      body.step();
      twin_body.step();
      double alone_flow = alone.radiated_flow();
      double body_flow = body.vocal().radiated_flow() - twin_body.vocal().radiated_flow();
      response.push_back((alone_flow - alone_previous) * tract::kStepRate);
      body_response.push_back((body_flow - body_previous) * tract::kStepRate);
      alone_previous = alone_flow;
      body_previous = body_flow;
    }
    std::printf("alone:");
    for (double hz = 500.0; hz <= 8000.0; hz += 500.0) {
      std::printf(" %.0f:%+.0f", hz, magnitude_db(response, hz) - magnitude_db(response, 1000.0));
    }
    std::printf("\nbody: ");
    for (double hz = 500.0; hz <= 8000.0; hz += 500.0) {
      std::printf(" %.0f:%+.0f", hz, magnitude_db(body_response, hz) - magnitude_db(body_response, 1000.0));
    }
    std::printf("\n");
    return 0;
  }
  // Scan 100 Hz - 6 kHz for peaks.
  std::vector<double> grid;
  for (double hz = 100.0; hz <= 8000.0; hz += 10.0) {
    grid.push_back(magnitude_db(response, hz));
  }
  double top = -1e9;
  for (double level : grid) {
    top = std::fmax(top, level);
  }
  for (std::size_t i = 1; i + 1 < grid.size(); ++i) {
    if (grid[i] > grid[i - 1] && grid[i] >= grid[i + 1] && grid[i] > top - 45.0) {
      std::size_t low = i;
      while (low > 0 && grid[low] > grid[i] - 3.0) {
        --low;
      }
      std::size_t high = i;
      while (high + 1 < grid.size() && grid[high] > grid[i] - 3.0) {
        ++high;
      }
      std::printf("peak %5.0f Hz  %+6.1f dB  bandwidth %4.0f Hz\n", 100.0 + 10.0 * i, grid[i] - top,
                  10.0 * static_cast<double>(high - low));
    }
  }
  // Level every 500 Hz (dB re the strongest peak).
  std::printf("levels:");
  for (double hz = 500.0; hz <= 8000.0; hz += 500.0) {
    std::printf(" %.0f:%+.0f", hz, magnitude_db(response, hz) - top);
  }
  std::printf("\n");
  return 0;
}
