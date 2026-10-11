# NATO reference work — 2026-10-11

The user rejected single-airway papa/mama take 039 as equally bad or worse and
asked for NATO codewords compared visually against human enhanced spectrograms.
Do not ask for more listening judgments. **No identical reference match has been
established.** A finite waveform, recognizable word, conservative junction, or
passing test is not evidence of natural speech.

The browser's previously recognized (but rattling) take-036 design remains the
rollback default. All work here is explicitly experimental. No recording from
that older engine is presented as a newly improved result.

## Reproduce the current experiments

From the repository root, render a fresh Papa with the current sheet engine:

```sh
node web/tube/render-speech.mjs papa out/new-papa \
  docs/experiments/nato-reference/body-053.json reference \
  --plan docs/experiments/nato-reference/papa-plan-052.json
```

`body-043.json` retains the 16 mm fold-flow channel with a larger front mouth.
`body-051.json` tests a 4 mm flow channel and 15 mm span across the slit.
`body-053.json` retains the 16 mm channel and reduces front-cavity radius from
16.5 to 12.5 mm to test the overly low first resonance. These are experiments,
not anatomy measurements or accepted speech presets.

Render all 26 NATO words to a new directory:

```sh
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python \
  tools/nato_visual_batch.py --out out/new-nato \
  --design docs/experiments/nato-reference/body-051.json \
  --plan-options docs/experiments/nato-reference/batch-options-052.json \
  --render-only
```

The batch uses the general reference score, including its Papa timing; the
separate Papa plan above contains the additional manually inspected timing
experiment. They are not interchangeable takes. Each output folder stores its
exact plan and design.

Then generate the complete enhanced fields and comparisons:

```sh
OPENBLAS_NUM_THREADS=1 /Users/joshuahkuttenkuler/miniforge3/bin/python \
  tools/nato_visual_batch.py --out out/new-nato --fields-only \
  --reference-cache out/nato-visual-040 --workers 2
```

A reference cache is accepted only when its WAV hash, channel, sample rate,
start/end samples, hop, and row count agree. The manifest updates after each
finished job, supports resuming fields, and distinguishes failed renders from
completed fields and unverified matches. Do not edit renderer source during a
render; its source hash guard deliberately rejects that take.

## Visual contract

Every generated candidate receives the preserved enhanced depiction. Use the
complete field at 48 kHz, hop 128, 684 frequency rows through approximately
8 kHz, magma, 70 dB. `compare_nato_fields.py` places both complete fields at one
pixel per cell with physical axes. No time warping, image resizing, grayscale,
crop, or ordinary-STFT companion is used. The NPZ retains original numerical
fields for inspection.

Display reference is 0.17808298734486477. Synthetic display gain is the fixed
0.1614252377244924 calibration established with take 040, kept across subsequent
takes and words. This is only a display calibration; raw PCM is unchanged.
Do not renormalize each candidate's brightness to conceal lost vocal energy.

## Findings and rejected candidates

- 038: projecting contact before force exchange did not remove irregular sound.
- 039: one serial air pressure field removed the old duplicate fold-pressure
  forcing, but the user explicitly rejected the resulting voice.
- 040: 22 of 26 original NATO renders completed; Charlie, Delta, Mike, and
  Whiskey exceeded the radiation guard. Those four later completed after
  localizing tongue constrictions and including actual constrictions between
  acoustic sampling points. Do not count 040 as 26 successful words.
- 041–044: shorter/wider upper tract, localized articulation and correct initial
  silence still left major formant, pitch and voicing errors. Charlie 042 is
  much longer than its reference; the later reference-score timing fixes that
  duration discrepancy, without proving a phonetic match.
- 045: varying elasticity across the source body produced a descending
  harmonic spacing close to the reference during the stressed vowel. The
  strongest synthetic band remained near 650 Hz instead of 800–900 Hz, and the
  first vowel was weak. Neither matching harmonic spacing nor internal motion
  established acceptable sound.
- 046: reciprocal moving-fold volume coupling satisfied its discrete volume
  balance test but produced a worse enhanced field, with broad excess energy
  and a long tail. It is retained behind `foldVolume`, disabled by default in
  the reference score. The remaining source junction is therefore incomplete;
  this is not a validated energy-conserving fluid/structure model.
- 047: batching the native material calls is bit-identical to 046. Its WAV hash
  is identical, so 046's exact enhanced field applies. Batching is not a sound
  improvement.
- 048: widening the lower throat did not repair the dominant resonance.
- 049–050: earlier first gesture, fold release during silence, and separate
  unvoiced abduction improve control semantics but still visibly fail the
  reference. Flow traces expose late phonation after fold adduction, low mean
  glottal flow, and ringing. The first vowel remains too weak.
- 051: separating the short flow-channel dimension from the slit span raises
  flow and changes the harmonic structure, but weakens the first vowel and
  shifts pitch above the reference. It is not a proven improvement.
- 052: fold adduction begins before the lip release. The first vowel becomes
  more visible than in 051, but its shape and the stressed-vowel spectrum remain
  substantially different from the human reference.

All 26 fresh renders in `out/nato-visual-052` completed without a model fault.
Alpha's first attempt was rejected by the source-change guard after a test file
changed; it was rerendered with frozen source. All 26 full-field pairs were inspected; none matches. The per-word observations
and exact WAV/comparison hashes are in `visual-audit-052.json`. The old 040 queue was stopped
after its first twelve reference fields to concentrate analysis on fresh takes;
already-running child field jobs were allowed to finish. Every old artifact is
preserved.

## Mechanism changes and limits

`CoupledAir` has one serial airway spanning the source body, mouth and lips,
with bidirectional nasal and under-tongue branches. Every non-fold material
patch owns one acoustic pressure junction. Slit flow integrates inertia and
viscous resistance along the actual channel; its force profile remains an
approximate attached/separated-flow model. It is not a replacement voice
waveform or a sampled reference source.

Local closures operate on a small part of the moving sheet instead of collapsing
an entire mouth segment. Acoustic edges include all intervening sheet apertures.
Optional fluctuating jet losses remain positive and are driven by existing flow;
the unforced system stays exactly silent with them enabled. They do not supply
arbitrary noise independent of airflow.

Muscle effort changes material elasticity continuously, with a 35 ms response,
and the existing bounded noise/sawtooth/bubble variation remains active. Effort
and fold-abduction gestures are independent from lip release. These gestures
require further reference calibration; their presence is not a quality result.

The native batching representation preserves every material coefficient and
transport history and has an exact sample-parity test. Current renders remain
below the 10x realtime goal (roughly 3–5x with enhanced analyses also running).
Do not claim the efficiency requirement is met. The 55-test suite includes the
original take-036 exact hash, native/reference mechanics, transport checks, and
new single-airway and radiation checks. Those checks establish implementation properties,
not realistic speech.

## Physical dimension references

The distinction between anterior-posterior length and vertical medial-surface
thickness motivates the short-channel experiment, rather than fitting those two
physical dimensions with one length. The parameter choice is still ours:

- [Measured medial-surface thickness methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC12014861/)
- [Vocal-fold geometry and layer variation](https://pmc.ncbi.nlm.nih.gov/articles/PMC6993787/)

These papers do not validate this implementation or its sound.

## Further mechanism trials: 053–066

No candidate in these trials has passed the human-reference visual comparison.
Exact render hashes and source hashes are in `receipts-053-066.json`; raw takes,
full fields, and comparisons are retained in the corresponding `out/nato-*`
directories. The tracked plans reproduce the parameter choices with current
source; the source frozen beside each original WAV reproduces that exact take.

- 053: reduced front-cavity radius; the dominant low band remained misplaced.
- 054: AA instead of AX in Papa's first vowel, consistent with the working
  transcription and the similar formant regions visible in both reference
  vowels. The first vowel still did not match.
- 055: closed the under-tongue branch during AA. This did not repair the first
  resonance. A separate terminal-flow trace found negligible nasal radiation
  during this oral vowel, rejecting nasal cancellation as the main cause.
- 056: established Alpha's initial posture with zero pressure during the
  preparation period. This reduced the initial setup transient but did not
  fix vowel structure or timing.
- 057–058: narrower posterior geometry moved the low resonance, but also
  changed source pitch through feedback. Lowering effort restored some pitch
  range, while leaving the wrong dominant bands. Neither trial is accepted.
- 059: localized fricative areas of 2.5/4/5 mm² for lip/tip/blade constrictions.
  The desired narrow opening exists, but fricative energy remains inadequate.
- 060: a causal aperture-reflection FIR reduced some upper ringing in Papa.
  Its main vowel bands still differ substantially from the human reference.
- 061: Foxtrot with preformed posture and aperture radiation. The initial F
  opening is in place before pressure arrives, but mean flow rises too slowly
  and the enhanced field still lacks the reference's fricative structure.
- 062: a measured 44-area human AA profile spanning the whole post-fold tract.
  It changes the resonances but produces too much upper energy, the wrong
  dominant low bands, and an extended decay. A measured tract from another
  speaker is not the NATO speaker's tract, nor validation of our propagation.
- 063–064: a propagation FIR whose coefficients sum to one preserves steady
  pressure instead of attenuating DC at every acoustic step. The stationary
  pressure test passes. Foxtrot 063 still has missing fricative energy and large
  release lines; this correction is insufficient to produce the desired voice.
- 065: continuous geometry-dependent oral loss replaces the 1.5 mm on/off
  threshold. The reduced loss is proportional to squared relative expansion,
  with the existing discharge factor retained. It vanishes in a uniform tube.
  This remains a quasi-steady approximation, not a resolved turbulent jet.
- 066: prescribed inlet pressure replaces the matched acoustic source
  impedance. In 065's initial F at 0.155 s, approximately 482 Pa requested at
  the source yielded only 177 Pa at the inlet. Holding the applied pressure
  raises flow through the F opening from about 30 to 50 ml/s at 0.16 s.
  This changes the inlet reflection as well as the steady flow; it is a
  distinct boundary hypothesis, not evidence that a human lung is an ideal
  pressure boundary.

The new options remain opt-in (`apertureRadiation`, `vowelAreas`, `airLossFIR`,
`geometricJetLoss`, `pressureBoundary`). The old browser default remains
unchanged. A stationary-pressure test, exact native batching parity, or higher
flow is never used to declare the speech natural.

### Geometry and radiation provenance

`story-2008-areas.json` contains numerical cross-sectional areas from Table I of
[Story 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2377017/), in cm², with the
original sampling lengths. The trial interpolates areas onto moving material
patches over the existing tract length. It uses geometry only, not reference
PCM or a replacement excitation. The coarse air edges still take the minimum
aperture among their material patches, so narrow shapes may occupy too much
acoustic length. This remains an unresolved representation error.

The aperture FIR follows the gamma-shaped approximation in equations 12–15 of
[Silva et al.](https://arxiv.org/abs/0811.3625). Fixed-radius coefficients are
nonnegative before the reflection sign, sum to one, and have the specified end
correction. The tests cover their fixed-radius passivity and frequency response.
They do not prove passivity under moving-radius modulation or validate mouth
radiation in this engine.

The latest full suite passes 55 tests, including the original take-036 hash.
Two replays also verify that precomputing acoustic neighbors leaves the entire
065 and 066 WAVs byte-identical. Current rendering is still below 10x realtime;
recent trials with enhanced analyses running measure about 3–6x. Native sheet
batching reduces crossings without changing material state, but the efficiency
requirement remains unmet.

### Completed visual review

All 26 words in batch 052 fail the reference comparison. The repeated defects
are wrong vowel-band locations and transitions, missing sustained fricative
energy, sharp artificial releases, and low-frequency decay beyond the human
reference. The audit records a separate observation and artifact hashes for
every word. No listening assessment was requested from the user.

Take 064 retains 062's wrong vowel structure despite DC-preserving propagation.
Take 065 visibly reduces several large vertical release spikes relative to
063. Its weak fricatives, vowel resonances, and timing still fail the human
Foxtrot. This is a localized artifact reduction, not a successful word model.
The earlier 043 Mike and Whiskey full comparisons were also reviewed and fail
on onset/timing, resonance, release spikes, and decay.

Take 066's complete field is also rejected. Prescribed inlet pressure increases
upper-band energy and fricative-region texture, but the sustained reference S
field is still missing, the vowel bands remain misplaced, and the final portion
runs late. Increased airflow is not a sufficient fix. All 053–066 candidates
have now been depicted and visually reviewed; none is promoted as a speech
quality result. The 065 reduction in release spikes is the clearest localized
improvement among the final transport trials.

LAN service was restarted with the updated source-freeze list and hop-128
renderer. HTTP responses at `http://192.168.1.18:8765/` return the exact local
coupled-air, articulation, radiation, and WASM module bytes. This verifies
serving over the host's LAN address; it is not a new test from the MacBook.
