# Speech iteration: render, depict, inspect, change

The user's required loop is: **every round generates a sound and depicts that
sound with the full enhanced analysis.** Inspect the audible/visible result
before selecting the next model change. Numerical correctness supports the
work; it does not establish speech quality and cannot substitute for this loop.

Use `tools/render_visual_round.py` with the scientific Python runtime. It
freezes the current local-engine executable, source, body and dictionary;
renders fresh complete speech; writes a listening copy with one declared fixed
gain; computes the preserved enhanced field; and checks that the depiction's
source hash equals the delivered WAV's hash. A round with a missing or failed
field remains incomplete. The command does not award an intelligibility score.

```sh
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python \
  tools/render_visual_round.py --text 'Papa.' --out out/visual-round-NNN
```

The complete enhanced field retains original magma colors, all frequency rows
and all computed time columns. Current round settings are 48 kHz, 128-sample
hop, 0–8 kHz, apertures 1024/2048/4096 and the original seven offsets. The hop is
a genuine analysis spacing, not image interpolation. Native-cell views and
integer magnification preserve the field; no ordinary-STFT companion,
grayscale conversion, excerpt substitution or display downsampling is used.

Compare the actual pattern of voiced regions, harmonic continuity, formant
motion, closure/release and noise with the human reference. A visual match alone
is not proof of intelligibility or realism. Log what is actually visible, the
next chosen change and the subsequent audio/depiction. Do not declare a word
working because its plan contains the expected phones or its samples are finite.

Round 001 is the freshly rendered current candidate with its unfinished motion
changes, not a new travelling-wave engine or a recovered working model. Exact
provenance is in `out/visual-round-001/run.json`; the full paired fields are
`papa-enhanced.npz` and `human-papa-enhanced.npz` in that directory. The human
reference remains valeatory's CC BY-SA 3.0 recording (see `samples/manifest.json`).

## Inspected rounds

- **001:** The full depiction shows narrow persistent bands and noise in the
  lead-in and between the syllables. The human reference has a brief, weak
  first voiced region and a stronger second one. Merely lengthening the first
  planned vowel is not justified by that comparison.
- **002:** Added local thermal/friction relaxation memories to test whether
  they removed the visible ringing. The enhanced depiction remained almost
  unchanged at the targeted bands. Rejected: the experiment was removed from
  the working engine, and its exact sources remain in the round snapshot.
- **003:** Returned to the 001 engine and used `bodies/word-audition.body` to
  omit the automatic preceding inhalation. The earliest diffuse noise is gone,
  but a broad noisy lead-in and persistent narrow bands remain. This is a
  configuration control, not a recovered word or an established improvement
  in intelligibility. The two voiced regions still do not resemble the human
  reference clearly. The next investigation should account for that excess
  noise and weak first voiced region before adding acoustic complexity.

Round 002 and 003 files are `speech.wav`, `enhanced.npz`, `enhanced.json` and
`round.json` within their respective `out/visual-round-NNN` directories. Every
recording was freshly generated during this iteration; no earlier recording
is substituted. The WAV-to-depiction hash is checked by the round runner.

The combined comparison contains the human reference and all three complete
fields, on the same 70 dB magma scale and display reference. PNG palette storage
is lossless: decoded RGBA pixels are required to equal the original field
colors exactly. The browser view was inspected with all 684 rows and 514
synthetic time columns present at integer pixel sizes; switching rounds was
verified. Full-duration human reference has 305 time columns. Different clip
lengths retain their real time axes, rather than being time-warped to match.

Focused local speech tests and their AddressSanitizer build passed for the
retained motion changes. That is an implementation integrity result only.
No round has received listening acceptance; no word is counted as working.

## Papa continuation after rejection of rounds 001–003

The user heard all three as the same and asked that the images remain internal
working tools. Subsequent rounds therefore keep complete audio and enhanced
depictions locally, without presenting a comparison viewer to the user.

- **004:** Coordinated lung pressurization with the first stop and anticipated
  glottal motor delay. A short harmonic first region appeared, but noise still
  occupied the intended closures. The physical trace showed that the lip slot
  never sealed: a section-centre sample missed contact at the mouth plane.
- **005:** Represented mouth-plane aperture separately from retained air
  storage. Fixed volume quadrature retains the partially open lip section.
  The prolonged noisy leak vanished, leaving brief transients and the final
  vowel. The first vowel was absent: the folds stayed at equilibrium.
- **006:** Preserved tissue position and momentum as the adduction equilibrium
  moved. This corrects instantaneous tissue motion, but alone did not restore
  the first vowel in the enhanced field.
- **007:** Used actual inertive glottal flow in the simplified lower-fold
  pressure force. Previously the force assumed instantaneous steady flow,
  inconsistent with the newly stored air momentum. First-vowel harmonics
  appeared after release, followed by a quiet closure and the second vowel.
- **008:** Changed the NATO entry to `P AA0 P AA1`, preserving the first open
  vowel in [NATO's pah-pah pronunciation guide](https://www.nato.int/content/dam/nato/webready/news/2010-2019/2017/12/21/20180111_nato-alphabet-sign-signal.pdf).
  Two open-vowel regions appeared, but a preparatory lip-closing transient
  preceded the word and the first region was short.
- **009:** Initialized the isolated utterance in its first stop closure before
  acoustic storage is initialized. Subsequent closures remain physical motion.
  Used the existing speaking-rate input at 0.8, rather than editing/cutting the
  audio. Its inspected full field has two harmonic vowel regions separated
  by closure, with the unrelated earlier transient absent. Synthetic timbre
  and narrow upper-frequency ridges remain. This is the current Papa candidate;
  human listening acceptance remains pending. Exact hashes and observations
  are recorded in `papa_visual_rounds.json`.

The local engine still uses fixed local air exchanges. It has not been replaced
by the abandoned coupled solver or an older speech engine. No human recording,
neural voice, waveform splice, or post-render repair supplies the synthesis.
The lower-fold force remains a reduced separated-jet model with its existing
non-suction assumption; this is not a complete fluid/tissue energy proof.

Unprompted transcription is a secondary diagnostic, never an acceptance score.
Whisper base.en returned “Saha” (004), “Ahh”/“Ah” (005/006), “the pa” (007), and
“Buh-bye!” (008/009). A second, larger recognizer, Whisper small.en, returned
“Papa” for 009 and for the human control. The recognizers were given neither
the intended word nor a NATO vocabulary prompt. The disagreement remains
recorded; no human listening approval is inferred from either result.
