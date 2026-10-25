#pragma once

// The space between the mouth and the microphone. A real recording is never
// dry: a small room returns a short reverberant tail and a faint background,
// so closures and pauses are never digital silence. The references (CMU
// ARCTIC, the NATO recordings) are close-miked: at -16 dB the tail of each
// vowel buried the weak consonants after it (/th/ read -21 dB against the
// speaker's -31.5; dry it reads -28.7), so the room stays faint.

#include <cstdint>
#include <vector>

namespace voice {

struct RoomOptions {
  double reverb_seconds = 0.25;   // RT60
  double reverb_level_db = -32.0; // reverberant energy relative to direct
  double noise_level_db = -58.0;  // background relative to the speech peak
  // A close microphone colours the voice: the proximity effect lifts the
  // low harmonics (low shelf), and a speech microphone's presence rise
  // lifts 2.5-5 kHz. Against the ARCTIC long-term spectrum ours was 3-4 dB
  // short below 315 Hz and 4 dB short at 2.5-5 kHz. Off since the folds
  // became physiological: the colouring had made up for the pressed
  // source, and on the new one put the fundamental 6 dB and 2.5-8 kHz
  // 3.5-6.5 dB above bdl; flat, every band is within about 4 dB.
  double proximity_db = 0.0;      // low shelf gain below proximity_hz (was 4)
  double proximity_hz = 250.0;
  double presence_db = 0.0;       // peak gain at presence_hz (about two octaves wide; was 6)
  double presence_hz = 3500.0;
  double presence_q = 0.7;        // 0.7: about two octaves between half-gain points
  std::uint32_t seed = 7U;
};

// Adds a feedback-delay-network room tail and a background floor in place.
void apply_room(std::vector<float>& samples, double rate, const RoomOptions& options);

}  // namespace voice
