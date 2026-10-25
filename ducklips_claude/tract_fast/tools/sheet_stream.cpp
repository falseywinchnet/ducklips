// The body language as a stream, spoken by the sheet tube: score lines
// arrive on standard input while the tube speaks, and the sound leaves on
// standard output as 16-bit 48 kHz mono samples, block by block (10 ms). The
// speaker never waits for the whole utterance: it speaks up to the time of
// the latest event it has heard.
// usage: say "text" x.wav --score /dev/stdout --plan-only --pitch 120 \
//          --phones ../tract/data/phones_sheet.txt | sheet_stream > out.raw
//   (play: sox -t raw -r 48000 -e signed -b 16 -c 1 out.raw -d)
// Environment: LEVEL scales the output (default 0.6: speech at 600-800 Pa
// peaks near 0.5 of full scale; 0.8 touched full scale).

#include "tract/sheet_speaker.hpp"
#include "voice/score.hpp"
#include "wave_out.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void flush_subnormals() {
#if defined(__aarch64__)
  std::uint64_t control = 0;
  __asm__ volatile("mrs %0, fpcr" : "=r"(control));
  control |= (1ULL << 24);
  __asm__ volatile("msr fpcr, %0" : : "r"(control));
#endif
}

// The output chain of one block: down to 48 kHz (two Butterworth biquads at
// 18 kHz, every other sample) and a 40 Hz high-pass (the walls' slow
// pumping), then the fixed level.
struct Output {
  Biquad first;
  Biquad second;
  double b0 = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
  double x1 = 0.0;
  double x2 = 0.0;
  double y1 = 0.0;
  double y2 = 0.0;
  double level = 0.6;
};

void prepare_output(Output& output, double rate) {
  output.first = lowpass(18000.0, rate, 0.5412);
  output.second = lowpass(18000.0, rate, 1.3066);
  double corner = std::tan(3.14159265358979323846 * 40.0 / 48000.0);
  double norm = 1.0 / (1.0 + std::sqrt(2.0) * corner + corner * corner);
  output.b0 = norm;
  output.a1 = 2.0 * (corner * corner - 1.0) * norm;
  output.a2 = (1.0 - std::sqrt(2.0) * corner + corner * corner) * norm;
}

void speak_block(tract::SheetSpeaker& speaker, Output& output, long steps, int decimate,
                 std::vector<std::int16_t>& block) {
  block.clear();
  for (long step = 0; step < steps; ++step) {
    double filtered = run(output.second, run(output.first, speaker.step()));
    if (step % decimate != 0) {
      continue;
    }
    double y0 = output.b0 * (filtered - 2.0 * output.x1 + output.x2) - output.a1 * output.y1 - output.a2 * output.y2;
    output.x2 = output.x1;
    output.x1 = filtered;
    output.y2 = output.y1;
    output.y1 = y0;
    double clipped = std::fmax(-1.0, std::fmin(1.0, y0 * output.level));
    block.push_back(static_cast<std::int16_t>(std::lround(clipped * 32767.0)));
  }
  std::fwrite(block.data(), sizeof(std::int16_t), block.size(), stdout);
  std::fflush(stdout);
}

}  // namespace

int main() {
  flush_subnormals();
  tract::SheetVoice voice;
  const long block_steps = static_cast<long>(std::lround(voice.rate * 0.01));
  const int decimate = std::max(1, static_cast<int>(std::lround(voice.rate / 48000.0)));
  Output output;
  prepare_output(output, voice.rate);
  const char* level_text = std::getenv("LEVEL");
  if (level_text != nullptr) {
    output.level = std::strtod(level_text, nullptr);
  }
  // The opening postures (events at time 0) come first; everything after is
  // pushed as it arrives.
  voice::Score opening;
  std::string line;
  bool started = false;
  double heard_until = 0.0;
  std::vector<std::int16_t> block;
  block.reserve(static_cast<std::size_t>(block_steps / decimate + 1));
  std::vector<tract::SheetSpeaker> holder;  // built at the first non-opening event
  holder.reserve(1);
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
      holder.emplace_back(opening, 1U, voice);
      started = true;
    }
    tract::SheetSpeaker& speaker = holder.front();
    speaker.push(event);
    heard_until = std::max(heard_until, static_cast<double>(event.time));
    // Speak everything up to the latest event heard.
    while (speaker.time() + static_cast<double>(block_steps) / voice.rate <= heard_until) {
      speak_block(speaker, output, block_steps, decimate, block);
    }
  }
  // The stream has ended: let the last postures sound for half a second.
  if (started) {
    tract::SheetSpeaker& speaker = holder.front();
    double end = heard_until + 0.5;
    while (speaker.time() < end) {
      speak_block(speaker, output, block_steps, decimate, block);
    }
  }
  return 0;
}
