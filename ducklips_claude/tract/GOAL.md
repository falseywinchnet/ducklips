# Goal (set 2026-10-13 at the owner's direction)

Work continuously without checking in with the owner. Stop only for a
decision the rules below do not cover, or for something that needs the
owner's hands.

## What to build

English text goes in. A parser finds natural spacing and breaths and
inserts those behavioural marks. Words become phonemes; phonemes become
configurations of the vocal system, with marks and timing. That language
streams to an engine that deforms the vocal system to satisfy it, in real
time and for several voices.

The body, as the owner specified it:

- Every section combines three different pieces: one of distinguishing
  elasticity, one of distinguishing mass, one of distinguishing support.
- All of them take part in edge-based eigenwave amplitude shifts and in
  instantaneous pressure waves.
- The pieces wander slightly under the effort of holding a shape.
- Feedback exists below the cords, below the tongue and in the sinuses.
- Energy comes through feedback, never from direct mass-spring
  stimulation: that always sounds mechanical, and mass-spring parts become
  oscillators.
- No expensive solvers.

## The bar

The voice must sound like a person to the owner. VocalTractLab sounds
natural; what exists so far does not beat Microsoft's Whistler TTS.

- **m2i** (`listen/m2i_new_engine.wav`, the segment-module engine) is the
  owner's best so far.
- The sheet engine (m22–m28) "warbles like a vibrator on their throat",
  and the voice is lost.
- **m30** (the segment-module engine at 750 Pa) "warbles less" but "still
  sounds like a mid early 2000s tts model".

## The lesson about measurement

The comparison and analysis measures used so far "do absolutely nothing
useful to identify what is wrong with the sound". They are long-term
averages, band levels, an identification classifier and slow pitch
statistics. Each one called a version fine that the ear rejected.

New transforms are needed, built to see what a listener hears:

- the fine time structure inside each glottal cycle;
- how the frequency bands move together;
- noise against periodicity, cycle by cycle;
- the shape of transitions.

Each one must be validated: a candidate transform counts only if it
separates real speech (bdl, the owner, NATO) from our engines and from
known-bad versions, as the owner's ear does. A measure that cannot tell m28
from bdl is not used to steer.

## Standing rules

- Heavy compute runs on the Mini.
- No commits.
- Downloads need explicit permission.
- No Apple recogniser and no whisper.cpp.
- VocalTractLab is GPL: study its mathematics only; do not run or copy its
  code. Its download holds Windows binaries and no synthesised audio.
- House style for first-party code.
- Use the epistemic labels.
- Delete only my own scratch files.
- Watch the Mini's disk.
- Keep `NOTES.md` current.
- Send listening files when a version is meaningfully better by the
  validated transforms.

## State (2026-10-13 evening)

- **Main line:** the segment-module engine (`tools/speak`, `stream`, `voices`; `Speaker` / `VoiceBody`), planned by `voice/build/say --plan-only --pressure 750 --phones tract/data/phones_tract.txt`.
  - It meets the spec: per-node cover / body / support pieces, the triangular glottis blown by the flow, the trachea, nose, sinus and under-tongue branches, and motor wander. It runs in real time for several voices.
  - The sheet engine (`sheet_speak`) is the owner-rejected branch ("vibrator").
- **Listening files sent today:** m31–m39 (`listen/`). m39 is the current best: from text, flowing articulation, sloped intonation, damped folds, a lossy tube and nose, tuned vowels, softer bursts, a breathy /h/, the room.
- **Identification** (four-view, three seeds): 37.5% (human ceiling 54.9).
- **How we look now:** four-view enhanced images (`voice/build/tf_analyze` + `voice/tools/tf_view.py`) and zoomed wideband (`voice/tools/spectrogram2.py`), beside bdl, VTL's synthesis and the owner's own reading of the test sentences (`ref/user_three.wav`). Then targeted measures for what the eye found (`/tmp/slice.py`, `/tmp/tilt.py`, `/tmp/bursts2.py`, `/tmp/phone_spectrum.py`).

## State (2026-10-14)

- Owner on today's `talk_remote.sh`: more human, "still sickly, warbling a bit"; "evening" and "human" not clearly legible.
- **Measured and fixed:**
  - the pitch plan swung on almost every syllable; gliding intonation now turns only at stressed vowels (contour reversals 2.95 → 1.21 per voiced second; people 0.8–1.4);
  - the adductor tremor leaked to cycle rate (one-pole band) and doubled the cycle jitter; the band now has two poles at each corner;
  - the glottis shut in only 60% of /iy/'s cycles; lower-edge rest 0.008. RAP 0.32%, p90 1.5% (owner 0.23%, 1.5%);
  - nasal murmurs: sinus ostia ×0.25, nose loss 0.4 (within 2 dB of bdl in five bands).
- **Validated transform:** `analysis/wmatch.py` (cycle periods from audio by waveform matching; agrees with the engine's true closures).
- `listen/m43_steady_voice.wav` sent. Identification 38.0%.
