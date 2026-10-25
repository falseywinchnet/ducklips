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
// Every variable and file the tool reads is documented where it is read:
// tract::settings_from_environment (src/sustain_run.cpp). The run itself is
// tract::SustainRun, shared with the in-process engine (bindings/), so the
// two cannot drift.

#include "tract/sustain_run.hpp"

#include <cstdio>
#include <string>
#include <vector>

int main(int count, char** values) {
  std::vector<std::string> arguments;
  for (int index = 1; index < count; ++index) {
    arguments.push_back(values[index]);
  }
  tract::SustainSettings settings;
  tract::ProcessEnvironment environment;
  std::string error;
  int status = 0;
  if (!tract::settings_from_environment(arguments, environment, settings, error, status)) {
    std::fprintf(stderr, "%s", error.c_str());
    return status;
  }
  tract::SustainRun run(settings);
  {
    tract::StreamFiles files(settings, true);
    if (!run.advance(settings, run.total(), files)) {
      return 1;
    }
  }
  run.finish(settings);
  std::fprintf(stderr, "%s: %.2f s, peak %.4g (normalized to 0.7)\n", values[2], settings.seconds, run.peak());
  return 0;
}
