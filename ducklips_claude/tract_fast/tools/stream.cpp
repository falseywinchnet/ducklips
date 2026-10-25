// The body language as a stream: score lines arrive on standard input while
// the body speaks, and the sound leaves on standard output as 16-bit 48 kHz
// mono samples, block by block (10 ms). The body never waits for the whole
// utterance: it speaks up to the time of the latest event it has heard.
// usage: say "text" x.wav --score /dev/stdout --plan-only | stream > out.raw
//   (play: sox -t raw -r 48000 -e signed -b 16 -c 1 out.raw -d)

#include "tract/speaker.hpp"
#include "voice/score.hpp"
#include "wave_out.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr long kBlockSteps = 960;  // 10 ms at the step rate

void flush_subnormals() {
#if defined(__aarch64__)
  std::uint64_t control = 0;
  __asm__ volatile("mrs %0, fpcr" : "=r"(control));
  control |= (1ULL << 24);
  __asm__ volatile("msr fpcr, %0" : : "r"(control));
#endif
}

}  // namespace

int main() {
  flush_subnormals();
  // The opening postures (events at time 0) come first; everything after is
  // pushed as it arrives.
  voice::Score opening;
  std::vector<voice::PartEvent> later;
  std::string line;
  bool started = false;
  double heard_until = 0.0;
  Biquad first = lowpass(18000.0, tract::kStepRate, 0.5412);
  Biquad second = lowpass(18000.0, tract::kStepRate, 1.3066);
  std::vector<std::int16_t> block;
  block.reserve(static_cast<std::size_t>(kBlockSteps / 2));
  tract::BodyLayout layout = tract::speech_layout();
  // The speaker is built at the first non-opening event.
  std::vector<tract::Speaker> holder;
  holder.reserve(1);
  // Output scale, fixed (a stream cannot normalise afterwards): radiated
  // sound peaks near 1e6 in body units for speech at 600-800 Pa.
  double level = 5.0e-7;
  while (std::getline(std::cin, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    voice::Score parsed;
    if (!voice::parse_score_line(line, parsed) || parsed.parts.empty()) {
      continue;  // marks and the duration line carry no posture
    }
    voice::PartEvent event = parsed.parts.front();
    if (!started && event.time <= 0.0) {
      opening.parts.push_back(event);
      continue;
    }
    if (!started) {
      opening.duration = 1.0e9;
      holder.emplace_back(layout, opening, layout.seed);
      started = true;
    }
    holder.front().push(event);
    heard_until = std::max(heard_until, static_cast<double>(event.time));
    // Speak everything up to the latest event heard.
    tract::Speaker& speaker = holder.front();
    while (speaker.time() + kBlockSteps * tract::kStepSeconds <= heard_until) {
      block.clear();
      for (long step = 0; step < kBlockSteps; ++step) {
        double filtered = run(second, run(first, speaker.step()));
        if (step % 2 == 0) {
          double clipped = std::fmax(-1.0, std::fmin(1.0, filtered * level));
          block.push_back(static_cast<std::int16_t>(std::lround(clipped * 32767.0)));
        }
      }
      std::fwrite(block.data(), sizeof(std::int16_t), block.size(), stdout);
    }
  }
  // The stream has ended: let the last postures sound for half a second.
  if (started) {
    tract::Speaker& speaker = holder.front();
    double end = heard_until + 0.5;
    while (speaker.time() < end) {
      block.clear();
      for (long step = 0; step < kBlockSteps; ++step) {
        double filtered = run(second, run(first, speaker.step()));
        if (step % 2 == 0) {
          double clipped = std::fmax(-1.0, std::fmin(1.0, filtered * level));
          block.push_back(static_cast<std::int16_t>(std::lround(clipped * 32767.0)));
        }
      }
      std::fwrite(block.data(), sizeof(std::int16_t), block.size(), stdout);
    }
  }
  return 0;
}
