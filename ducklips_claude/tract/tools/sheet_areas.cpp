// The sheet tube's cross-sections while a score is spoken: steps the speaker
// to a time and prints every node's distance above the glottis, band, radius
// and area, so a posture's real area function can be checked.
// usage: sheet_areas SCORE SECONDS
//   Columns: node, cm from the cord's top, band name, radius (mm), area (cm^2).

#include "tract/sheet_speaker.hpp"
#include "voice/score.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace {

bool read_score(const char* path, voice::Score& score) {
  std::ifstream input(path);
  if (!input) {
    return false;
  }
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (!voice::parse_score_line(line, score)) {
      return false;
    }
  }
  voice::sort_score(score);
  return true;
}

}  // namespace

int main(int count, char** values) {
  if (count < 3) {
    std::fprintf(stderr, "usage: sheet_areas SCORE SECONDS\n");
    return 2;
  }
  voice::Score score;
  if (!read_score(values[1], score)) {
    std::fprintf(stderr, "cannot read %s\n", values[1]);
    return 1;
  }
  double until = std::strtod(values[2], nullptr);
  tract::SheetVoice voice;
  tract::SheetSpeaker speaker(score, 1U, voice);
  while (speaker.time() < until) {
    speaker.step();
  }
  const tract::SheetTube& tube = speaker.tube();
  const tract::SheetDesign& design = speaker.design();
  constexpr double kPi = 3.14159265358979323846;
  double position = 0.0;
  bool above_cord = false;
  for (int node = 0; node < tube.count(); ++node) {
    const tract::SheetBand& band = design.bands[static_cast<std::size_t>(tube.band_of(node))];
    double spacing = band.length / static_cast<double>(band.points);
    if (above_cord) {
      position += spacing;
    }
    if (band.role == "cord") {
      above_cord = true;
      position = 0.0;
    }
    double radius = tube.radius(node);
    std::printf("%3d %6.2f %-26s %6.2f %7.3f\n", node, 100.0 * position, band.name.c_str(), 1000.0 * radius,
                1.0e4 * kPi * radius * radius);
  }
  return 0;
}
