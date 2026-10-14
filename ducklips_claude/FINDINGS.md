# DuckLips findings

Environment for every **MEASURED** value below, unless stated otherwise:

- Machine: Apple M4 Mac mini, macOS 26.5, one thread.
- Build: Apple clang 21, `-O2` Release, C++20, no fast-math.
- Engine: 48 kHz, 64 oral and 24 nasal sections, modes below 12 kHz (at most
  40), 48-sample control blocks, passive reprojection.
- Date: 2026-10-08.
- Reproduce with `./build/ducklips_tests`, `ducklips_say` and
  `ducklips_analyze` as listed.

## Numerical core

| Check | Result |
|---|---|
| Tree inertia counts vs dense Jacobi eigenvalues (31-node branched tree) | exact at every midpoint |
| Cold eigenvalues vs dense Jacobi, real network with velar port open | 3.0e-13 relative; 20 of 20 modes below 12 kHz |
| C-orthonormality of computed shapes | 4.4e-16 |
| Warm start after a small deformation vs cold solve | 3.7e-14; all 20 modes accepted warm, 80 factorizations |
| Uniform 17.5 cm closed-open tube vs quarter-wave modes | 500.0, 1499.7, 2498.4, 3495.7 Hz (0.12% worst) |
| Analytic formant sensitivities vs central differences | 8.8e-8 relative |
| Passive reprojection across a lip closure | energy ratio 0.979 (bounded by 1) |
| Bernoulli flow and implicit cavity step | residual 3e-20 |

**/a/ bandwidths (MEASURED).** Mode frequency, then bandwidth B split into
viscous, wall, radiation and thermal shares, all in Hz:

| Mode | B | Viscous | Wall | Radiation | Thermal |
|---|---|---|---|---|---|
| 733 | 32 | 6 | 16 | 5 | 5 |
| 1357 | 39 | 12 | 4 | 17 | 6 |
| 2536 | 64 | 16 | 1 | 37 | 9 |
| 3702 | 118 | 21 | 1 | 85 | 12 |

Walls dominate F1 and radiation dominates above F2, as in the acoustic
phonetics literature.

## Speech

**Speed (MEASURED).** "The quick brown fox…" (6.7–7.5 s of speech) renders
14–20× faster than real time:

| Voice | Real-time factor | Modes | Factorizations per geometry step |
|---|---|---|---|
| baritone | 0.070 | 20–24 | 80 |
| alto | 0.057 | 17–19 | 66 |
| child | 0.049 | 14–17 | 56 |
| gravel | 0.071 | 21–24 | 82 |
| whisper | 0.062 | 18–21 | 72 |

About 0.4% of mode solves fall back to bisection. Vowel inversion at prepare
takes 0.15–0.28 s per voice.

**Vowel postures (MEASURED).** The posture's own eigenfrequencies reach all
16 targets within 1.4–4.8%. The inversion is regularized toward the neutral
tract and stops at that tolerance.

**Vowels in the audio (MEASURED).** LPC order 12 on a 12 kHz copy, median over
voiced frames, baritone, 1.0 s sustained vowels:

| Vowel | Measured F1/F2/F3 | Target | Δ |
|---|---|---|---|
| IY | 310 / 2280 / 2840 | 342 / 2322 / 3000 | F1 −9% |
| IH | 350 / 2025 / 2710 | 427 / 2034 / 2684 | F1 −18% |
| EH | 550 / 1775 / 2560 | 580 / 1799 / 2605 | |
| AE | 545 / 1860 / 2510 | 588 / 1952 / 2601 | |
| AA | 680 / 1315 / 2525 | 768 / 1333 / 2522 | F1 −11% |
| AO | 595 / 970 / 2545 | 652 / 997 / 2538 | |
| UH | 395 / 1130 / 2475 | 469 / 1122 / 2434 | F1 −16% |
| UW | 310 / 1010 / 2365 | 378 / 997 / 2343 | F1 −18% |
| AH | 610 / 1185 / 2595 | 623 / 1200 / 2550 | |
| ER | 410 / 1370 / 1765 | 474 / 1379 / 1710 | low F3 reached |

F2 and F3 land within about 5%. F1 reads 5–18% low; see open problems.

**Pitch (MEASURED).** The feedback loop holds intent: on a sustained /a/,
measured 114.8 Hz against intended 114.7 Hz. Sentence medians are 117 Hz
(baritone, base 110), 219 Hz (alto, 205), 279 Hz (child, 265) and 93 Hz
(gravel, 88). Question endings rise visibly in the spectrogram.

**Turbulence levels (MEASURED, rms relative to /a/; calibrated).**

| Phone | Level | Natural speech, approx. |
|---|---|---|
| /s/ | −7.8 dB | −10 |
| /ʃ/ | −13 dB | −8 |
| /f/ | −29 dB | −22 |
| /h/ | −13 dB | −15 to −20 |

**GUI.Forms Audio (MEASURED).** As a generator voice in `AudioEngine`, an
offline render of 3.35 s differs from the direct render by 0.0 (bit-identical).

**Crest factor (MEASURED).** On the "Hello world…" sentence, peak/rms is
0.31/0.022, about 23 dB. Natural speech is roughly 18–20 dB. The largest
event is the /k/ release burst.

## Failures kept on the record

Each was found by measurement, diagnosed, and fixed or rejected.

1. **The folds never started (REJECTED as a model; fixed).** With equal rest
   areas, the two-mass equilibrium is an exact fixed point. It is unstable
   above threshold, but nothing perturbed it: zero cycles, and the "120 Hz"
   calibration was its fallback. Fix: glottal-jet turbulence acts on the
   lower fold surface. This is physical, and it is the same noise as
   aspiration.
2. **Charge-conserving reprojection blamed for blow-up (REJECTED
   diagnosis).** A sentence diverged at about 2.5 s with energy doubling
   every millisecond. Charge conservation (pressure q/C as C→0) was the first
   suspect, and a passive invariant with a proof was introduced. The actual
   cause was item 3. After that fix, charge/momentum and pressure/flow also
   run the sentence stably (peaks 0.30 and 0.31). Passive stays the default
   because it is the only invariant with a bound. Pressure/flow measurably
   adds energy at a lip closure (×1.12).
3. **Duplicate eigenvectors (fixed).** A failed warm refinement leaves the
   vector on a neighbouring mode's exact eigenvector. Three inverse-iteration
   passes could not leave it when the nasal and oral eigenvalues were about
   1e-5 apart. Measured basis orthogonality defect: 1.0. Fix: perturbed
   restart, deflation against all accepted modes, iteration to convergence,
   and an all-modes overlap test. Defect is now 2.6e-7, and bisections fell
   from 13,052 to 783.
4. **Labial frication never radiated (fixed).** A series source in the lip
   branch did not reach the radiated flow (`dU/dt = (p + p_s)/M`). /f/ was
   45 dB below /s/.
5. **Gesture teleport clicks (fixed).** A new gesture in an active slot
   (L→D, Z→DH) moved place and area within one 1 ms block. Slots now glide.
6. **Closure kink clicks (reduced).** A linear activation-to-area map closes
   with a kink that splatters energy into the highest modes. A smoothstep
   contact cut the peak from 1.0 to 0.63.
7. **Adaptive geometry sub-stepping as the click cure (REJECTED
   hypothesis).** Hypothesis: closures scatter energy because one projection
   per block is non-adiabatic. Up to 8 interpolated steps per block (19% more
   solves) left the peak unchanged at 0.64. The mechanism stays in, because
   it bounds the geometry change per step, but it did not cure the clicks.
8. **Cavity overfill by a lagging mean flow (fixed).** Feeding the cavity
   the low-passed glottal flow overfilled it during /k/ closures (1018 Pa
   against 800 Pa lung pressure). Glottal flow reversed, and the acoustic
   cavity charged and released a broadband click. The cavity now receives
   the true flow; only the radiated outflow derivative is smoothed. Peak fell
   from 0.64 to 0.31.

## Round 2 (2026-10-08, after the architect listened)

**GIVEN** (listening report): voiced output often sounded clipped. "DuckLips"
was heard as "octlips", and voiced passages were mostly unintelligible. The
whisper voice was intelligible ("I am a physical model of a talking person").

**Diagnosis (MEASURED, per-phone probe of the baritone demo).** Whisper
working pointed at the voicing path, not the tract.

- 7 of 39 vowels never closed the glottis. After a voiceless stop the folds
  re-adducted too slowly, so short unstressed syllables were breathy and
  swallowed: physi-*cal*, *ta*lking.
- The /d/ of "DuckLips" never sealed: its peak intraoral pressure was 9 Pa.
  The /z/ before it shares the tongue-tip slot. The gliding tip changed degree
  too slowly, and /z/'s release cue arrived after /d/'s closure cue and
  cancelled it.

**Fixes.**

- Degree changes in place are fast (450 rad/s) while place glides stay at
  160 rad/s.
- A release is now dropped when the same articulator moves straight into its
  next gesture.
- The velum seals before an obstruent that follows a nasal.
- Stops have a minimum closure; short vowels have a minimum length.
- Voice onset is −10 ms after a release, or +25 ms aspiration before a
  stressed vowel.

Result (MEASURED): 1 of 39 vowels breathy, no unsealed stops. The
baritone's 90th-percentile F0 went from 500 Hz (garbage frames) to 140 Hz.

**Connected speech (added).** Weak forms, nasal place assimilation,
palatalization, cluster loss, geminate merging and flapping. Example:
"Did you want to put the water in the cup, and then go to bed?" becomes
`D IH JH UW0 W AA N T AX0 P AH T DH AX0 W AO T~ ER0 … AX0 N DH EH NG G OW T~ AX0 B EH D`.

**Energy modulation (added).** The `effort` cue scales the input energy:
lung pressure, with firmer adduction above 1. The planner raises it on
accented syllables and lowers it phrase-finally. The `loudness` cue
modulates the output waveform; the planner uses it to keep inhalations
under the speech.

**Side cavities (added).**

- Piriform fossae, as a 2-cm leaf beside the epilarynx.
- Maxillary sinuses, as a Helmholtz leaf (ostium 0.14 cm² into 30 cm³) off
  the nasal tract.

Both are leaves of the same tree, so the solver is unchanged. MEASURED: all
tests pass, including the dense cross-check on the larger tree. The /a/
posture gains a fossa-localized mode at 4285 Hz that does not radiate,
consistent with the known 4–5 kHz piriform antiresonance.

Cost: the fossae lower formants. Open vowels then underfit (REJECTED fix: an
8-mode shape basis changed nothing). Widening the open-vowel lip apertures
(/ɑ/ 6 cm², /æ/ 5 cm²) and loosening regularization gives /ɑ/ 6.3%, /æ/
6.7%, front vowels' F3 about 5% low, and the other ten vowels within 1.5%.

## Round 3 (2026-10-08): the tongue, motion, and the whistle

**GIVEN** (listening): "DuckLips" was heard as "bucklips", /n/ was not heard,
there was a persistent high-frequency sound that does not belong in speech,
and "without a tongue you can't do a /d/". The architect also specified that
mouth parts move at different speeds, with a small ramp-up and ramp-down.
Speech recognizers (Apple's or whisper.cpp) are not to be used.

**Missing tongue (MEASURED, then fixed).** The old tip "constriction" was a
pinch over a vowel tube; /d/ showed its alveolar F2 for about 10 ms, close to
/b/. It is replaced by a midsagittal tongue (`tongue.hpp`):

- A fixed outer wall, and a tongue body (circle) with a 2.2 cm blade and a
  rounded tip pad. Jaw, blade raising and tongue root are posture
  parameters.
- Areas come from the wall-to-tongue gap, using smooth minima so the tissue
  has no corners.
- /t d n/ place the tip pad on the ridge, and a blade too short to reach
  drags the body.
- The drag feeds a body-offset tracker with tongue-body inertia, so after
  release the body drifts back.
- /k g/ press the dorsum.

The cosine vowel basis is gone. A penalized 4-term log-area residual stands
in for anatomy the coarse geometry lacks.

Result (MEASURED):

- /d/ seals, with intraoral pressure up to 640 Pa.
- At release, F2 falls 1455 → 1265 Hz over about 50 ms; /b/ stays flat
  around 1140 Hz.
- Vowel fits: worst 8.6% (/ɔ/), and /ʌ ə ʊ u ɝ/ within 3%. The /d/ locus
  starts lower than the textbook 1800 Hz.

**REJECTED along the way:**

- Body and jaw alone could not reach the front vowels: /i/ F2 was 1932
  against 2322 Hz.
- A flat 5–6 cm² hypopharynx chamber lowered F3 in every vowel.
- A thin blade line could never seal, because it crossed the ridge over less
  than one section.
- A purely geometric drag snapped the body back within 1 ms after release.

**The whistle (MEASURED, then fixed).** Above about 6 kHz the 1-D tube's
highest modes were trapped in narrow throat sections and stopped radiating.
Their bandwidths fell to 130–190 Hz at 9.0, 9.9, 10.8 and 11.3 kHz: narrow,
persistent whistles on every glottal pulse. The fix is a cross-mode loss
above 4.5 kHz, standing in for the 3-D scattering the 1-D column lacks (a
modelling choice), plus a 10 kHz top mode. The 8–11 kHz band dropped from
−51 to −63 dB.

**Motion.** Every articulator follows a third-order critically damped
trajectory: zero velocity and zero acceleration at the start, no overshoot.
Each part has its own speed limit from the body (tongue tip 230, lips 180,
glottis 170, tongue body 75, velum 65, jaw 45, lungs 40 rad/s). A cue cannot
make tissue move faster than its limit.

**Nasals (MEASURED).** /n/ and /m/ are 4–5 dB below the vowels, which is
normal; their absence to the ear is not a level problem. After the tongue
change, /n/ shows the same falling F2 transition as /d/.

**Boundary-surface probe (separate worker; boundary/NOTES.md).**

- Delayed surface exchange reproduces closed-box modes within 1%.
- Line-of-sight-only links cannot make the open-end reflection (REJECTED).
- The full kernel is unstable when marched explicitly in time.
- Tissue membrane modes lie below 610 Hz and do not matter for sibilants.
- 40–150× the cost of the 1-D voice.

## Open problems

- **Intelligibility is unmeasured.** No speech recognizer is installed on
  either machine, and judging by ear is the architect's call. The `say`
  control spectrograms are structurally comparable but are not a
  perception test.
- **F1 reads low (HYPOTHESIS).** The 5–18% may be LPC bias at a 115 Hz F0
  (the posture's own F1 is within 5%), or a source-tract interaction lowering
  F1 during the open phase. It needs a formant measure robust to harmonics,
  for example analysis-by-synthesis against the known modes.
- **/ʃ/ is too weak relative to /s/.** Both inject at the same teeth position.
  A per-gesture obstacle strength or a better jet model is needed.
- **Lexicon coverage.** The core dictionary is about 250 words, and the
  letter-to-sound rules are deliberately crude (they flag their guesses).
  The CMU dictionary loads with `--dictionary`, but no file is vendored;
  fetching one is an architect decision.
- **Partial physics.** No lateral channel for /l/, no paranasal sinuses
  (nasal antiresonances), no subglottal tract, no diffuse late reverberation,
  and no mouth directivity.
- **Calibrated constants.** Turbulence gains, static wall compliance, port
  geometry and gesture timings were tuned against levels and spectrograms,
  not measured.

## House style

The author has reviewed every file under `research/ducklips` against
`planning/PROGRAMMING_HOUSE_STYLE.md`; nobody has audited it independently.
It uses no `auto`, lambdas, structured bindings, ranges pipelines or
fast-math. Pointer access is written `(*p).m`. The real-time path does not
allocate.

Known deviations:

- Some conditional expressions use ternaries.
- `TreeOperator::diagonal`, `parent_coupling` and the network's
  `nodes`/`branches` are public data with documented fixed sizes.
- The planner's `CueWriter` keeps a reference member.
- The analysis tool's FFT uses a compound `for` header.
