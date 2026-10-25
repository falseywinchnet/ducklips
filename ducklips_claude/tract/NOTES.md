> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. For the past week fit was substituted for understanding. Progress comes only from structure: the physics, the direct inverse, and knowing why something works or breaks. Search scaffolding added since 2026-10-19 (arrival grids, time sweeps, joint guards, sane clamps, block polishes, vowel revisits) is labelled REJECTED as a method: it is fit, not understanding. Round-trip dB numbers from it show only that the model can be inverted; they are not progress.

# Segment-module tract — notes

Design: `../learned/segment_design.md`. Sources of the ideas: `../learned/`.
Built on the M4 Mini (`~/ducklips/tract`, `cmake -S . -B build`), tools run
from this folder (they read `data/`).

## Milestone 1 — sustained vowels (2026-10-13)

**Engine** (`include/tract/`, `src/`):
- Tube of three-piece modules on a 0.3646 cm node grid (one-step air delay at
  96 kHz). Elastic piece: lining sheet with axial and hoop tension; its
  implicit-midpoint step is one tridiagonal sweep per step (the exact sum of
  all of its modes). Mass piece: tissue body on distributed springs. Support
  piece: frame with stiffness, articulation offset and holding effort. Air
  junction per node with an implicit lining–air exchange (Astra's port solve);
  VocalTractLab viscous loss per link; R‖L lip radiation; output = d/dt of lip
  flow (+ optional skin emission of the lining).
- Cord module after VocalTractLab's triangular glottis: wedge-shaped rest gap
  pinned at the anterior commissure, continuous contact with growing contact
  damping, backward-Euler masses, tension Q (masses /Q, springs ×Q, coupling
  constant), Bernoulli in / 20% recovery out, plus glottal air inertance.
- Lungs as a matched pressure reservoir under 12 cm of ring-module trachea.
- Variation of causes (Astra's drift / sawtooth / bubble): breath ±3%,
  tension ±1.5%, fold rest gap ±0.0015 cm, wall support ±10% (30% shared).
- Glottal turbulence after VTL (v³√A dipole 1.5 cm above the glottis).

**Findings** (MEASURED with `tools/sustain`, Story 2008 shapes, 120 Hz, 800 Pa):
- Without glottal air inertance the folds were blown to a static opening
  (area 0.034 cm², steady 92 cm³/s): no oscillation. With it they
  self-oscillate at the intended pitch: period 8.1 ms, flow rising over about
  3 ms to 242 cm³/s and cut off in about 2 ms, full contact about 30% of the
  cycle; 31 consecutive cycles within 0.1% in peak flow.
- First wall values (lining 0.2 g/cm² on a 180 Hz spring, body on 126 Hz)
  swallowed F1 (/i/ F1 650–700 Hz, /u/ 640–860, pitch pulled to 135 Hz).
  Rigid walls give the tube's own formants: /a/ 684/977, /i/ 298/2152,
  /u/ 309/739 Hz. Recalibrated (body about 1.5 g/cm² on a soft, heavily damped
  40 Hz support; lining tied to it with its own modes in the kHz range):
  /a/ 690/953 (B1 98), /i/ 296/2078 (B1 144), /u/ 315/740 (B1 150).
- Aspiration as first ported (low-passed flow ÷ instantaneous area) blew up
  during closure (v ≈ 500,000 cm/s; HNR −2 dB, shimmer 39%). Fixed with the
  instantaneous jet speed capped by √(2Δp/ρ): HNR 21–28 dB.
- Voice quality: jitter 0.16–0.20%, shimmer 0.17–0.22% — more regular than a
  human sustained vowel. Pitch on /i/ 134.5 Hz for 120 intended (open-loop
  tension; source–tract interaction). The "level wobble" at F0/4–F0/5 in the
  checker is an artefact of 10 ms envelope windows beating with the period.
- PENDING: owner's listening verdict on `listen/m1_vowels.wav`.
- Owner verdict on `m1_vowels.wav`: "sounds like a vocoder monk."

## Milestone 1b — after "vocoder monk" (2026-10-13)

Diagnosis (MEASURED, `tools/transfer`, `tools/sustain`, band and cycle checks):
- Wall friction acted on absolute lining velocity (B1 360 Hz); now relative
  to the body (B1 90). The lining's own resonance near 3.5 kHz absorbed
  F3–F5; coupling raised to 4e8 (about 10 kHz): F3/F4/F5 bandwidths 60–180 Hz.
- Source alone (`SOURCE=` flow derivative): −4 dB at 250, −12 at 1 kHz, −24
  at 2 kHz re H1 — about −6 dB/oct, a normal modal slope. The 2 kHz valley in
  /a/ output (−43 dB) comes from the Story /a/ shape (F2 983, F3 2964), not
  the source.
- Side branches: the piriform recess (2 cm³) lowers /a/ F1 732 → 633 Hz
  wherever it is joined (1.5 or 2.5 cm, 1.2 or 0.6 cm²). Story's MRI areas
  exclude it, so the lowering is kept as real (HYPOTHESIS: correct in sign;
  size unverified). The under-tongue pocket lowered /i/ F2 2164 → 1874 and
  F3 3058 → 2380 while open; it now opens only with a raised tongue tip
  (`set_under_tongue`), closed for vowels. Sinuses have no effect with the
  velum shut.
- Too regular: jitter 0.16%, shimmer 0.16%. Variation of causes was only slow
  (≥ 80 ms). Added physiological tremor, 2–15 Hz random force (unit-rms band of
  two one-pole low-passes) on three muscles: stretching (0.5% of tension),
  adductors (0.0007 cm rest gap) and breathing (1.5% of lung pressure).
  Klatt's three-sine flutter removed (a repeating pattern). Result /a/ /i/ /u/:
  jitter 0.45–0.65%, shimmer 1.6–2.2% (normal human ranges).
  Tremor alone: jitter 0.31%, shimmer 0.5%; breath tremor drives shimmer.
- Jet turbulence on the folds' upper surfaces (VTL aspiration pressure):
  no measurable effect at 1× or 4× — REJECTED as a jitter source at this level.
- Auditory pitch feedback: at each closure the heard period trims tension
  (8% of the log error per cycle, bounded ±30%). /i/ 134 → 120 Hz, /u/
  133 → 120, /a/ 123 → 119 for 120 intended.
- `listen/m1b_vowels.wav`: /a/ /i/ /u/, 1.3 s each, intended pitch gliding
  128 → 104 Hz, released by letting breath go and parting the folds over the
  last 180 ms. PENDING: owner verdict.
- Owner verdict on `m1b_vowels.wav`: "sounds more human but still machined.
  about 2x as good."

## Milestone 2 — text to speech through the segment-module body (2026-10-13)

**Pipeline** (all MEASURED working on "Hello, my name is Papa", NATO words,
"The quick brown fox", arctic_a0003): text → the front end (`../voice`
`say --plan-only --score`: parser, breath and phrase marks, phonemes,
learned bdl durations and pitch, articulator score) → `tools/speak` →
`Speaker` (`src/speaker.cpp`): one critically damped motor per score part
(pulls toward its target; its command wanders with the effort of holding a
posture away from rest, the front end's tuned rates), once per ms setting the
support frames of the 48-node vocal tube through the front end's tract
geometry, the cord's arytenoid posture, breath, velum, under-tongue
opening and intended pitch (the ear closes the pitch loop). 3.3–3.5× real
time on the Mini (one voice).

**New physics needed to speak** (each found by a measured defect):
- Links carry a true series element (viscous + Borda-Carnot jet loss on
  expansion + series sources, closed form): pressure builds behind closures
  (5200 dyn/cm² behind /p/, released over ~6 ms). Wave attenuation alone
  never reflects. Vowel formants unchanged.
- Oral jet turbulence ported from the front end's body (narrowest place past
  the epilarynx, growth/decay along the jet, incisor obstacle, labiodental,
  Strouhal-centred band): /z/ /s/ hiss in place.
- The frame moves, it does not jump: rest radius follows its ramped target
  (critically damped, 160 Hz), masses recompiled while moving, and a moving
  frame pushes air (surface × frame speed). Jumps at the 1 ms control rate
  rang the tube as a 1 kHz comb and clicked at every movement. Tension,
  lung pressure and port areas likewise glide per step.
- Side ports (velum, sinus necks, piriform, under-tongue) are short ducts
  with the same series element. As bare lossless one-step lines a nearly
  shut velum integrated the pressure difference across it in its two waves;
  on opening it released ~1000 cm³/s into the nose against the pressure
  (nose 3000 dyn/cm², pharynx −800) and buzzed at 1 kHz. Fixed: comb +55 dB
  → +14 dB (noise floor) at that moment.
- Lungs: the bronchial tree absorbs sound (matched end) but passes steady
  flow with ~1.5 dyn s/cm⁵; as a plain matched port (16 dyn s/cm⁵) an open
  glottis dropped subglottal pressure 7500 → 2300 and every vowel after a
  voiceless consonant was 4–8 dB weak. Source raised by (Z − R) × flow below
  20 Hz.
- Piriform recess 2 cm deep (zero near 4.5 kHz, Dang & Honda); 3 cm deep it
  cut vowels off above 3.6 kHz.
- Glottal aspiration spectrum: Strouhal-centred like the jet (VTL's 1 kHz
  low-pass left /h/ 7 dB weak above 2 kHz).

**Comparisons** (MEASURED):
- NATO band-contour error vs owner (26 words, DTW): old engine 5.12/4.45/5.22
  (low/mid/high), new 5.78/5.59/6.00. New is ~1 dB worse in every band; the
  old engine was fitted to this metric, the new one has not been tuned.
- Long-term spectrum ("The quick brown fox", re 0.1–1 kHz): new 4–5 kHz −29,
  5–6 kHz −26; old −22/−22; bdl a0003 −22/−22; owner −18/−16. The new body
  is 5–10 dB dark above 4 kHz. Not walls (rigid: same), not branches (≤3 dB),
  not fold damping (≤4 dB). The source derivative falls ~10 dB/oct above
  1 kHz.
- bdl vowels: harmonic-to-between ratio 0–5 dB from 1 to 6 kHz; ours 11–25
  (partly bdl's pitch movement inside the window). HYPOTHESIS: too-clean
  upper harmonics are part of "machined".
- REJECTED (this form): more VTL aspiration (×30 matches vowel noise but
  /h/ 30 dB too loud); closing-phase turbulence (deceleration-driven jet
  instability; per peak flow): uneven — /i/ +18 dB at 4–6 kHz where /a/ /ae/
  barely move. Left in code, gain 0.
- `listen/m2_new_engine.wav` (four sentences, new engine) and
  `listen/m2_old_then_new.wav` (each sentence old then new, then bdl
  arctic_a0003). PENDING: owner verdict.

### Milestone 2, continued — where the brightness and continuity go (2026-10-13)

- Glottal air mass now stepped by the trapezoid rule (lossless; 0.95 decay
  of its pressure so a shut glottis cannot hold a Nyquist ripple). Backward
  Euler had made it a resistance growing as f² (≈42 dyn s/cm⁵ at 5 kHz).
  Effect on the voiced high band: small (≤2 dB). Kept as the correct step.
- Probe fixes: the transfer probe's plain DFT leaked a settling steady flow
  across the spectrum (a false "flow flattens the formants"); it now subtracts
  a twin run. Steady or 120 Hz pulsing flow changes the bare tube's transfer
  by < 1 dB (Borda-Carnot losses are 0.1–0.25% of Z at /a/'s narrowings).
- MEASURED: the whole body at rest has the bare tube's transfer; the same
  body while phonating (twin-body impulse difference) is 20–30 dB lower at
  4–7.5 kHz. Not aspiration, not jet noise, not walls, not the branches,
  not the delivered-vs-cord flow, not an absorbing port of 0.03–0.1 cm².
  HYPOTHESIS: the open glottis (kinetic resistance ≈ 45 + its air mass,
  close to the tract's ≈ 80 dyn s/cm⁵) absorbs the upper formants during
  each open phase; the high band rises with lung pressure (600 → 1100 Pa:
  5–6 kHz −50 → −35 dB re 0.1–1 kHz), as a voice brightens when louder.
  Same-sentence LTAS (arctic_a0003): new −28/−31 dB at 4–5/5–6 kHz, bdl
  −22/−22, old engine −22/−25; 6–8 kHz matches.
- An LPC (order 18) inverse filter whitens every source alike (bdl's and
  ours): useless for comparing source tilt. Needs IAIF-style separation.
- Continuity: bdl voices straight through phrases (voice bar through
  /dh v n m/, formants gliding); ours (and the old engine, same front end)
  comes in blocks separated by silent gaps and clicks. Per-phone voiced
  share, arctic_a0003, ours vs bdl's labels: /ay/ "time" 56% vs 94%,
  /iy/ "twentieth" 70% vs 100%, /ax/ 40% vs 100%, /th/ 0% vs 67%, /v/ 0% vs
  25%. Partly label conventions (bdl's aspiration sits inside the stop),
  partly late re-adduction after voiceless stops (planner adducts ~45 ms after
  release, plus voicing build-up).
- A swept closure pumps air: /th/ → /t/ ("twentieth time") slides the tip
  contact 17.0 → 15.0 cm in 28 ms while closing; with frames that displace air
  this compresses the mouth to 9400 dyn/cm² (> lungs) and drives flow back
  into the trachea. Physically right for that motion; the motion is not human
  (a tip releases the dental contact and makes a new alveolar one). Front-end
  score issue, not yet changed.
- Open-phase absorption confirmed (MEASURED, twin-body impulse, /a/): the
  phonating body's high band depends on where in the cycle the impulse lands;
  landing in a closed phase (pressed posture) it is the bare tube's (+27 dB at
  5 kHz re 1 kHz), in an open phase 20–30 dB lower. At 5 kHz the open glottis
  (≈ 39 + j129 dyn s/cm⁵: jet resistance, trachea, air mass) against the
  epilarynx (≈ 85) returns only ≈ 58% of the power per reflection.
- Closed quotient from glottal area (< 0.003 cm²): 40–42% for modal at 100–
  120 Hz (human EGG range 0.4–0.5); lower-edge contact alone reads 27%.
  So the folds do not stay open too long.
- Voicing continuity over arctic_a0003: bdl 60% of the span voiced (27 runs,
  median 30 ms); new 49% (24, 35 ms); old 43% (23, 40 ms). Less of a gap than
  the spectrograms suggested.
- Open question: bdl's 4–6 kHz band in vowels is largely aperiodic
  (harmonic/between 2–5 dB; ours 14). The two ways tried to add voiced
  turbulence failed (uniform aspiration: /h/ 30 dB too loud; deceleration-
  driven: /i/ overshoots). Waiting on the owner's ear before the next step.
- bdl /hh/ (51 tokens): −14.5 dB overall, −27 dB above 2 kHz re the
  utterance's loudest vowel. Ours at VTL aspiration 1: −23/−32. Aspiration
  default now 2 (/h/ −22/−26). Pressure-source aspiration reaches /i/ but
  not /a/ (×5: /i/ 4–6 kHz +6 dB, /a/ unchanged).
- REJECTED: glottal flux noise (a monopole fluctuation of the glottal flux,
  rms share of the flow): 3% already makes /h/ as loud as a vowel (−3 dB)
  and floods /i/ above 4 kHz; broadband flux noise radiates rising 6 dB/oct.
  Left in code at 0.
- Opening thump: the tube began at the neutral 3 cm² layout and glided to the
  score's first posture within milliseconds; the air it pushed was the loudest
  event of every utterance (0.7 of full scale) and every m2 file was
  normalized by it. The speaker now settles the frame in the opening posture
  before any sound (`Tube::settle_geometry`). `listen/m2b_*`: the same four
  sentences re-rendered (also aspiration 2, trapezoidal glottal air mass),
  levelled by speech loudness instead of peak.
- Enhanced field (voice/enhanced_specgram, 0–2 kHz, arctic_a0003 at 16 kHz):
  bdl's harmonics sweep in long curves with the intonation; ours run flatter
  with small fast wiggles (tremor). Median F0 ours 114 Hz, bdl 129.
- Intonation (front end, shared by both engines): the engine tracks the
  intended pitch within ±0.2 st per vowel; the flatness was the plan's. bdl
  pitch refit (bdl_full): held-out correlation 0.55 (rms 1.95 st vs flat
  2.33). Planner now widens the learned factors by `pitch_spread` 1.5 and
  takes bdl's sentence median, 122 Hz (was 110). Intended per-sentence
  vowel range 3.8–5.0 → 5.4–7.5 st; realized 6.0–7.8. Side effect: with
  pitch, speech pressure rises (+10 Pa/Hz): 2–4 kHz now 4–7 dB above bdl on
  arctic_a0003 while 4–6 kHz stays 5–9 dB low. NATO vs owner: mean 5.79 →
  5.88 dB (single words carry little intonation). `listen/m2c_new_engine.wav`.

### Voicing through consonants (2026-10-13)

Yardstick: `/tmp/corpus_cmp.py N first` (on the Mini) renders N bdl ARCTIC
texts through front end + this engine and compares with bdl's recordings
(long-term bands re 0.1–1 kHz; voiced share of the speech span);
`/tmp/vclass.py` splits voiced share by phone class (ours: planner phones;
bdl: labels). 24 sentences (a0101–a0124), MEASURED:

| step | 4–5k | 5–6k | voiced diff | vowel | voiced stop | voiced fric |
|---|---|---|---|---|---|---|
| bdl | | | | 81% | 61% | 42% |
| start (m2c) | −1.7 | −2.5 | −0.12 | 74% | 36% | 32% |
| prescribed pharynx widening +35% | | | −0.12 | 73% | 44% | 31% |
| tongue volume conserved | −3.5 | −4.3 | −0.11 | 78% | 27% | 31% |
| voicing gaps → modal posture | −3.6 | −4.6 | −0.10 | 79% | 26% | 32% |
| + wall pull 0.25 cm | | | −0.09 | 79% | 38% | 35% |

- The parabolic tract geometry does not conserve volume: /aa/ → /iy/
  through a /d/ closure widened the pharynx with nothing shrinking, and the
  frame's air displacement sucked the cavity behind the closure to −1.2 kPa,
  pulling the folds shut. The tongue is incompressible: over 3.5–16 cm the
  frame now moves air between nodes with zero net volume
  (`Tube::set_conserved_span`).
- Prescribed widening (frame) sucked −1.7 kPa: REJECTED. Wall release
  (support stiffness down 80–95%): < 3 points. The active part, muscles
  pulling the pharyngeal walls outward through the support anchor
  (`VoiceBody::set_wall_pull`, the tissue follows through its own mass
  against the air), holds ~3 kPa across the folds through a /d/ closure.
- The score's closure-voicing gap (0) had mapped to a posture twice as open
  as modal: the folds were inflated apart by the pressurized airway and
  leaked a steady 75 cm³/s. All voicing gaps (≤ 0.002) now map to modal.
- Remaining: voiced stops 38% vs 61%, voiced fricatives 35% vs 42%, voiced
  share −0.09; the 4–6 kHz deficit grew to −3.6/−4.6 dB over these
  changes. A faint 1 kHz comb remains in breath pauses (low level).
- Breath and adductor holds now glide per step (2 ms) instead of stepping at
  1 kHz. The faint 1 kHz comb in pauses is ~50 dB under the speech and comes
  mostly from articulator wander moving the frame; left.
- Through a /d/ closure the folds sat inflated apart (0.04 cm², a slow
  non-colliding wobble): the voice bar 32–40 dB under the vowel. Extra
  adduction during voiced closures (`set_closure_press`, up to 0.008 cm):
  no change in voiced share — kept at 0.004 (harmless), effect unproven.
- Intervocalic /b d g/ closure level re the louder neighbouring vowel, 16
  each: ours −28 dB, bdl −11.8. Flaps: the planner's flap rule fires
  (/t d/ between a vowel and an unstressed vowel → DX, a 28 ms tap), but the
  learned duration model had no "dx" factor (bdl's labels write flaps as t/d)
  so every flap took a generic stop's 80–95 ms. Flaps now take their own
  28 ms. Voiced share diff −0.09 → −0.08; vowels 80% (bdl 81%).
- Piriform and the 4–6 kHz deficit (24 sentences): without the recess
  4–5 kHz is +0.2 dB re bdl (with it −3.5), 5–6 kHz −3.1 (−4.5), but 1–2 kHz
  +3.5 (+2.2). A lossier lining (cover loss 0.15 → 0.5 → 1.0) changes nothing:
  the zero is the air column's and the lining barely moves. Size and opening
  are the remaining handles; not fitted without an anatomical reason.
  Lining loss left at 0.5 (mucosa).

### Stability and speed (2026-10-13)

- A 21 s passage died at 20.7 s (non-finite output). A shut port's area glides
  toward zero exponentially and never reaches it; ~1.4 s after the
  under-tongue space last opened ("forest") the area underflowed, ρc/area
  overflowed and the port returned NaN. Ports below 1e-6 cm² are now shut.
- Speed on the M4 Mini, one voice, 21.11 s passage (MEASURED; spectra and
  voicing statistics unchanged to ≤1 dB / 1 run at each step):
  | step | time | real-time factor |
  |---|---|---|
  | start | 6.91 s | 3.1× |
  | viscous law every 16 steps (pow was 21%), noise pows inlined, getenv out of the loop, subnormals flushed | 4.77 s | 4.4× |
  | nose+sinuses / under-tongue not stepped 50 ms after their port shuts; moving nodes recompiled every 8 steps | 3.65 s | 5.8× |
  | -mcpu=native -ffp-contract=fast | 3.5–4.0 s | 5.3–6.1× (noise) |
  One voice per core: several voices in real time on the Mini's 10 cores.
  Remaining cost: Tube::step_air, step_lining, refresh_geometry.
- Multiple voices (`tools/voices`, one thread per voice, own seed, pitch
  moved 3 semitones per further voice), M4 Mini, 4.28 s sentence:
  4 voices together 1.02 s wall (4.2× real time), 8 voices 1.58 s (2.7×).
  (4 performance + 6 efficiency cores.)
- Streaming (`Speaker::push`, `tools/stream`): score lines on stdin while the
  body speaks, 10 ms blocks of 16-bit 48 kHz on stdout; the body speaks up to
  the latest event heard. Same long-term spectrum as the file render within
  1–2 dB; fixed output scale 5e-7 (peaks 18.8k–27.0k of 32767 on the four
  test sentences, none clipped).
  `say "text" x.wav --score /dev/stdout --plan-only | stream > out.raw`

### Voice quality in running speech (2026-10-13)

`/tmp/vjit.py` (jitter, shimmer, H1–H2, HNR in 80 ms vowel middles) and
`/tmp/vhnr.py` (harmonic-to-between per band, pitch-stable vowel middles),
24 bdl sentences, ours (planner phones) vs bdl (labels). MEASURED:

| | jitter | shimmer | H1–H2 | HNR | 1–2k | 2–3k | 3–4k | 4–6k |
|---|---|---|---|---|---|---|---|---|
| bdl | 1.47% | 6.72% | −3.4 | 6.5 | 15.7 | 10.6 | 8.5 | 4.6 |
| VTL modal posture | 2.70% | 6.90% | +4.8 | 5.8 | 20.6 | 14.5 | 11.3 | 8.8 |
| 0.002, −0.004 (now) | 1.71% | 6.76% | +1.0 | 4.9 | 18.5 | 11.5 | 9.3 | 6.3 |
| 0.003, −0.006 | 1.54% | 6.93% | +0.4 | 5.5 | 19.1 | 14.3 | 11.1 | 8.3 |

- An 8 dB too-strong first harmonic (H1–H2 +4.8 vs −3.4) is a soft, hooty,
  "monk" voice; the voicing posture (upper edge pressed) takes 3.8 dB of it.
  HYPOTHESIS: part of bdl's negative H1–H2 is his recording chain's low cut.
- More aspiration (3–4×) barely moves 1–4 kHz HNR: the clean harmonics are
  not missing noise.
- Cost: 1–4 kHz long-term bands +2.6..+3.0 dB re bdl (was +2.0..+2.1).
- NATO vs owner (26 words, mean of band errors): 5.88 dB after the pitch
  change; 6.67 dB now. Of the separable changes since, each helps or is
  neutral on its own (posture VTL → firmer: 6.82 → 6.67; wall pull and
  press on vs off: 6.95 → 6.82; tongue conservation off: 6.56). The ~0.8 dB
  regression is unexplained (the owner's voice is higher than bdl's, whom
  the prosody now follows; contour timing and fricative levels dominate the
  worst words). Recorded, not chased.
- Further along the same line (pressing the upper edge lowers H1–H2, the
  lower edge raises it): 0.004/−0.008 → H1–H2 +0.1, jitter 1.58%, shimmer
  7.1%, HNR 1–2/2–3/3–4/4–6 kHz 18.8/12.0/10.7/6.9, voiced share −0.07,
  vowels 82% (bdl 81%), voiced stops 47% (bdl 61%). 0.006/−0.01: H1–H2 +0.4,
  voiced −0.06, 4–6 kHz HNR 8.2. Default now 0.004/−0.008.
  NATO vs owner with 0.004/−0.008: 6.60 dB (6.67 at 0.002/−0.004). `listen/m2f_new_engine.wav`.

### Breath at phrase ends (2026-10-13; front end, shared)

- NATO words the new body loses most against the old (owner reference):
  Quebec +4.2, X-ray +3.9, Papa +3.4, Foxtrot +3.4, Mike +2.9, Oscar +2.8,
  Charlie +2.8 dB — stops. Spectrograms vs the owner: the final /k/ of
  "Mike" was followed by 250 ms of flat 2–8 kHz hiss (the owner: a compact
  burst and ~50 ms of aspiration); the final /k/ of "Quebec" released no
  burst at all.
- Cause 1: the plan held lung pressure through every phrase's pause, so an
  open glottis hissed through it. Now the breath is let go 50 ms after a
  phrase's last sound (each phrase sets its own breath at its onset).
- Cause 2: the plan fades the breath through the vowel before a final
  voiceless stop (a fix for the front end's body, whose "Foxtrot." burst at
  7× the level around it), but whether it took effect depended on event
  order: "Mike" kept full pressure (a later voiceless-pressure event undid
  the fade), "Quebec" released at 1400 dyn/cm² (no burst). Option
  `--final-stop-breath` (used for this body: nato_compare ENGINE=tract,
  corpus_cmp, listening renders) keeps the breath to the release; the
  release then drains it in milliseconds.
- Result: final /k/ and /t/ burst (short, broadband; the owner's is 60–80 ms
  and concentrated at 1.5–4 kHz). NATO vs owner 6.60 → 6.56 dB.
  Remaining worst: Foxtrot, X-ray, Mike — mid band ~10 dB; the owner's /m/
  murmur in "Mike" is long and strong, ours barely shows.
- Phrase-initial /m/ of "Mike": voicing starts ~30 ms in while the breath is
  still rising (pressure starts 40 ms before the phrase, 110 ms rise).
  Longer breath lead, NATO vs owner: 0.04 s 6.56, 0.08 s 6.62, 0.12 s 7.49
  dB (an early rise reads as an early onset): kept 0.04.
  `listen/m2g_new_engine.wav`: with the phrase-end breath and final-stop fixes.

### Glottal noise spectrum (2026-10-13)

- The final /k/ burst of "Mike" (80 ms after release), band levels re the
  strongest, 0.5–1.5 / 1.5–3 / 3–5 / 5–8 kHz: owner −5/0/−4/−18; ours
  −19/−9/−4/0 — and it was nearly all glottal noise (frication off: no
  change; aspiration off: −10 dB at 3–8 kHz). The velar jet exists (×5
  gain: +7 dB at 1.5–3 kHz) but lay ~10 dB under the glottal noise.
- /hh/ shape (dB re strongest band, 0.3–1/1–2/2–3/3–5/5–7.8 kHz), bdl (267)
  0/−13/−7/−4/−11; ours (slit-height corner) 0/−16/−20/−9/−6. The glottal
  noise corner, 0.2 v / width, used the slit height (area / 1.3 cm): corner
  3–12 kHz. On the jet's hydraulic diameter √(4A/π) (the jet reaching the
  false folds and epiglottis is round): corner ~1–2 kHz; with gain 6 the
  shape is 0/−13/−13/−2/−9; gain 8 puts /h/ above 2 kHz at −29 dB re the
  loudest vowel (bdl −27). Default now: round jet, gain 8.
- 24 sentences: 4–5 kHz −3.5 → −2.8 dB re bdl; voice quality unchanged
  (H1–H2 +0.1, jitter 1.56%, shimmer 6.86%). NATO vs owner 6.56 → 6.10 dB.
- Still: the burst's 5–8 kHz stays strong (the velar jet's own corner uses a
  slit height too); /hh/ 2–3 kHz −6 dB re bdl.
  `listen/m2h_new_engine.wav`.
- Final /k/ of "Mike" with the round-jet glottal noise: ours 0/−6/−5/−3
  (0.5–1.5/1.5–3/3–5/5–8 kHz), owner −5/0/−4/−18. Still mostly glottal
  noise; the velar jet (the owner's 1.5–3 kHz peak: the front cavity rung by
  the release jet) sits ~12 dB under it. NEXT: velar burst strength — the
  jet model's weights were fitted to sibilants (incisor obstacle) in the
  front end's body; a velar jet has no obstacle within reach.

### Fricative and burst levels (2026-10-13)

`/tmp/bursts.py` (40 bdl sentences; levels re the utterance's loudest
vowel, MEASURED):
- Fricatives (middle 60%): at frication 0.6 /s/ /sh/ −22 dB vs bdl −13;
  at 1.5: /s/ −14, /sh/ −14, /f/ −34 (−32), /th/ −34 (−29), /z/ −16 (−13),
  /v/ −12 (−16). Frication default now 1.5. Long-term 4–7 kHz re bdl now
  +1..+3 dB (was −3..−5).
- Stop releases (50 ms from 10 ms before the stop's end): ours /p/ −22,
  /t/ −16, /k/ −22; bdl −12, −13, −14. Jet gain barely moves them; a faster
  opening of the constriction (×2, ×4) neither. With an alignment-free
  measure (loudest 10 ms within ±50 ms of the stop's end) the gaps are /p/ 8,
  /k/ 5, /t/ 0 dB — but that window reaches the next vowel. Unresolved:
  needs release-aligned windows on both sides.
- NATO vs owner with frication 1.5: 6.10 → 6.02 dB. `listen/m2i_new_engine.wav`.

### Burst spectra (2026-10-13)

`/tmp/bursts2.py`: releases found the same way on both sides (sharpest
2.5 ms rise of 2–8 kHz energy within 40 ms of the stop's end), 15 ms burst.
40 sentences, MEASURED. Band shape 0.5–1.5/1.5–3/3–5/5–8 kHz re strongest,
burst level re loudest vowel:

| | bdl | ours |
|---|---|---|
| /p/ | −1/−5/−6/−9, −21 dB | −7/−9/−10/−1, −22 dB |
| /t/ | −7/−5/−3/−6, −18 dB | −17/−7/−1/−4, −16 dB |
| /k/ | −2/−6/−6/−14, −20 dB | −17/−13/−4/−1, −22 dB |

- Burst LEVELS match within 2 dB (the earlier "10 dB weak" was the window
  missing the release). Shapes do not: ours put the energy at 5–8 kHz,
  bdl's /p/ and /k/ are low- and mid-heavy. Onsets: ours /p/ rises 50 dB
  per 5 ms, bdl 34.
- Cause: the jet-noise corner, 0.2 v / slit height, with the slit near zero
  at a release → corner at the 12 kHz clamp. Widening the jet where the
  turbulence forms (`jet_spread`, cm added to the slit height): 0.08 moves
  /p/ to −3/−3/−12/−3 and /k/ to −12/−10/−4/−3, but /s/ goes from
  −31/−15/−7/0 (bdl −28/−22/−1/0) to −19/−8/−5/0 and /sh/ loses its
  3–5 kHz. Splitting free-jet and obstacle (incisor) noise with the spread
  only on the free jet: /s/ still degraded (−23/−11/−7/0, level −18 vs −13),
  stops no better. REJECTED both; default spread 0 (one coherent source,
  as before).
- /f/ shape wrong in every variant: ours 0/−5/−10/−14, bdl −17/−12/−8/0.
  HYPOTHESIS: glottal noise (gain 8, round-jet corner ~1 kHz) through the
  wide-open glottis dominates /f/. NEXT.

### /f/ (2026-10-13)

- Fine bands (0–7.8 kHz) of /f/ middles: ours −3/−4/−3/−8/−17/…/−23 (energy
  below 1.5 kHz), bdl −2/−16/−17/−15/−13/−12/−11/−11/−6/−2/−2 (rising to
  8 kHz). Not glottal noise (aspiration off: still low-heavy), not the lip
  slit width (1.2/2.5/4 cm: no change). Trace of "a fee": lip opening 0.135
  cm², 300 cm³/s, only ~70 Pa behind the lips (a real /f/ builds 400–600
  Pa) — the jet at the lips ~1100 cm/s, its turbulence weak. The front end's
  tuned F area (0.18) was fitted on its own body.
- `data/phones_tract.txt` (the front end's table with this body's values;
  used by speak renders, nato_compare ENGINE=tract, corpus_cmp): F area
  0.06 → shape −18/−9/−3/0 (bdl −17/−12/−8/0); labiodental source 0.6 → 1.2
  → /f/ −33 dB (bdl −32), /th/ −31 (−29). Lip slit width 2.5 cm kept
  (physical; no measurable effect).
- 40 sentences now: fricatives /s/ −14 (−13), /sh/ −14 (−13), /f/ −33 (−32),
  /th/ −31 (−29), /z/ −17 (−13), /v/ −12 (−16); voiced share −0.06; voice
  quality unchanged; NATO vs owner 6.03 dB.

### The ringing (2026-10-13)

- Owner (on m2h/m2i): "some kind of high pitched faint ringing about twice a
  second."
- MEASURED on "Papa, papa" with both noise sources off: a 12.0 kHz ring
  decaying over ~15 ms at every lip closure and release, and a weaker one
  at 5.98 kHz. 12 kHz = 96 kHz / 8 and 6 kHz = 96 kHz / 16: the speed work
  had moving nodes recompiled every 8 steps and the viscous law evaluated
  every 16, so while a constriction moved fast its pieces and losses jumped
  periodically. Fixed: nodes whose frame moves faster than 0.5 cm/s are
  recompiled every step, links next to them get the exact viscous law every
  step, and elsewhere the viscous law glides between its 16-step
  evaluations. High band at the closures: ~44 dB → −5 dB. Speed 5.8× →
  5.1× real time.
- The "too bright" stop bursts were largely this ring. Release shapes
  (0.5–1.5/1.5–3/3–5/5–8 kHz): /p/ 0/−6/−9/−9 (bdl −1/−5/−6/−9), /t/
  −13/−5/−2/−8 (−7/−5/−3/−6), /k/ −8/−8/−3/−5 (−2/−6/−6/−14). Onsets still
  steeper (39–46 dB per 5 ms vs 29–34). NATO vs owner 6.11 dB (6.03 before;
  within run noise). Jet spread and corner-ceiling experiments above were
  measuring through this ring; their rejections stand only as measured.
- `listen/m2j_new_engine.wav`.
- 1 kHz comb sweep (`/tmp/artifacts.py`, `/tmp/combscore.py`, 21 s
  passage): in speech the comb at exact kHz multiples averages +2.4 dB
  (chance harmonic alignment; nothing to fix); in quiet windows it was +15 dB
  — the muscle causes (tension, rest offset, breath) arrived as a 1 kHz
  staircase. They are now interpolated across each tick (`cause_from_/to_`)
  and the lung command is ramped across the tick: pauses between phrases are
  silent; the opening inhale keeps +8.8 dB of comb (from +12.1), from
  articulator wander moving the frame, ~50 dB under speech. Spectra
  unchanged.

### "Microsoft speech engine" / "belly" (2026-10-13)

- Owner: "he still sounds like a really bad microsoft speech engine";
  "hello my name is papa and i speak with a belly."
- Breath pressure 750 → 600 Pa for this body (`say --pressure 600` in all
  its tools): 40 sentences, every band 1–7 kHz within ±1 dB of bdl (was
  +2.5..+3.7 at 1–4 kHz); voiced stops 52% (bdl 56%); H1–H2 −0.4.
- "body" → "belly": the flap never touched (tip closest 0.21 cm²): a stop's
  28 ms approach and 50 ms voiced release launched 16 ms apart. Flaps now
  tap ballistically (12 ms close, held to the phone's end, 15 ms release):
  < 0.01 cm² for ~13 ms. Planner (shared).
- Vowel chart (hVd, LPC, vs Hillenbrand men): /æ/ 711/1547 (588/1952) — the
  front end's tuned table set AE to place 7.5/1.5 (fitted on its own body);
  refitted on this body in `data/phones_tract.txt`: 9.0/0.75 → 624/1736 in
  "had". /ʊ/ F1 579 → 500 (lips 0.4, 9.75/0.5). /aa/ 756/1328 (768/1333).
  /u/ and /ʌ/ sit at the front end's fronted/connected targets by design.
- The Mini's disk is full (189 GB used on Data; ~/Developer 34 GB,
  ~/ducklips 20 GB, /private/var/folders 14 GB). Only this session's scratch
  renders in /tmp/es were removed (84 MB) to finish m2k.
- `listen/m2k_new_engine.wav`.
- NEXT (the owner's verdict points here): split test — the new body driven by
  bdl's own timing and pitch (prosody transplant) against our planner's, to
  locate the "Microsoft" quality between the front end and the body.

## Milestone 3 — the sheet tube (2026-10-13)

Owner: "your engine doesnt really sound any different than your old engine
... the elastic, mass, and support are not elastic, mass, and support. that
would be a mass and spring model which is horribly confused. the sheet tube
has THREE SECTIONS, one of each of which is tuned for support, mass, and
elastic properties, all of which are part of the tube's networks ... air
acts on the surface of all. the experiments where we worked on relaxation
dynamic behavior are ... necessary"; "mass-spring systems sound terrible.
they always wind up becoming oscillators." Base chosen by the owner:
correct acoustics.

- `src/sheet_tube.cpp`, `include/tract/sheet_tube.hpp`: Astra's connected
  sheet tube (web/tube, 42c7337) ported: one supported periodic sheet
  (axisymmetric family; implicit midpoint = her cached receptance, one
  tridiagonal sweep), sections in sequence with their own material and
  mass/elasticity/support factors and ring supports, air waveguide network,
  tissue velocity-wave network, collocated scalar ports, 25 µm contact,
  reservoir with 35 ms response and wander, per-section elasticity variation
  propagated to damping, impedance and tissue delay; connected systems.
  MEASURED parity with her JS engine: identical F0 (312/304 Hz), level,
  level variation and band shape; waveform correlation 0.9999 over the
  first 50–200 ms (later chaotic divergence). 10–17× real time.
- Her air was slow: every hop at least one 24 kHz sample over 2–6 mm nodes
  (~50–140 m/s). With nodes at c/rate (`acoustic_spacing`, rounded down)
  her cord section still self-oscillates (46 of 108 settings, 87–290 Hz).
  In the axisymmetric quotient the ring/diagonal air paths become one-sample
  self-loops that add compliance: `axial_air` keeps one axial path of the
  same total weight. Her outlet was near-matched (25% reflection):
  `radiating` adds the piston load (R ∥ L, trapezoid); the sound is d/dt of
  the radiated flow.
- Her sections are light (elastic 0.2, soft ~7, firm ~10 kg/m²): the sheet's
  own modes filled 180–700 Hz and swamped the formants (transfer probe
  `tools/sheet_transfer`). Tract modules now: elastic (fold, ~10 kg/m²,
  elasticity 3), support (firm, ~15 kg/m², support 3.2), mass (soft,
  ~31 kg/m²) in sequence, sized to Story's area functions (`make_vowel_tube`).
  Eight modules at 48 kHz: /a/ 850/1050/2980 (Story 730/1090/2440), /i/
  390/2160/2950 (270/2290/3010), /u/ 410/810/2380 (300/870/2240).
- Voiced (cord 16 mm, elasticity 0.6, mass 2.0, offset +0.3 mm, support 0.2,
  1200 Pa): /a/ F0 120, formants 836/1070/2970, jitter 0.18%, shimmer
  0.78%, H1–H2 −2.7, HNR 19; /i/ 118 Hz, 371/2161/2953, jitter 0.15%,
  shimmer 1.6%; /u/ formants right but voicing unstable (shimmer 23%).
- `listen/m3_sheet_vowels.wav`: /a/ /i/ /u/ from the sheet tube, then
  Astra's original below-mouth stage.
- NEXT: Astra's somatic holding control (stress relaxation with memory,
  opposing twitching motor pools, delayed PI correction) for the cord's and
  sections' holding; /u/ stability; articulation of ring supports from the
  front end's score.

### Holding, and the balloon (2026-10-13)

Owner, on the sheet tube (m3/m4): "it sounds like a rubber balloon now. or a
chicken toy."

- OBSERVED cause: the sheet cord never closed. Its narrowest diameter swung
  1.27–1.67 mm (±13%) around an open 0.64 mm ring: a balloon neck. The sheet
  tube's air is a linear waveguide, with no flow-dependent loss at a
  narrowing (no Bernoulli, no jet separation), so the cord could only
  flutter acoustically and never cut the flow off.
- `orifice 1` (design key): the air across the cord section is a
  quasi-steady orifice between the junctions below and above:
  rho L dv/dt = (Pa − Pb) − A v (Za + Zb) − rho v|v|/2 − 8 pi mu L v / A.
  This is one closed-form quadratic per sample (backward Euler in the
  velocity; as A closes, the viscous term stops the flow). The cord's
  surfaces feel the Bernoulli pressure up to the narrowest ring and the
  pressure above it past that point. With the orifice on, the cord section
  keeps its own node count, 8, not the acoustic spacing. Astra parity is
  unchanged with the key off (correlation 1.0000).
- Robust closing cords (`/tmp/robust2.py`; widest adduction ranges with
  steady voicing): length 3–5 mm, radius 1.3 mm, mass 1.5–3, elasticity
  0.4–0.6, support 0.1–0.2. Each voices across 400–600 µm of adduction, and
  F0 is set by mass and elasticity.
  - At adduction −0.35 mm the flow is pulsed, peaks at about 180 cm³/s and
    closes each cycle.
  - At −0.45 mm the pulse takes the textbook skew (slow rise, fast closing).
  - Static /a/ at −0.35 mm: jitter 0.24%, shimmer 1.8%, H1–H2 +3.4,
    HNR 25.
- Before the orifice, a robustness search on the old acoustic cord found
  that its 0.3 mm candidate voiced only within ±20 µm. Mass 4 with support
  0.8 voiced over 360 µm.
- Holding (`tract/somatic_hold.hpp`, ported from Astra's motor_hold and
  somatic_control):
  - Two opposing 64-unit pools with irregular discharges and 15/35 ms
    twitches.
  - Sensing with a 20 ms filter and 30 ms delay, then PI correction.
  - A relaxing strain memory with a wandering log-rate.
  - Effort is measured from a relaxed posture: an open cord (+1 mm), a
    slack stretch (0.6), a neutral 3 cm² tract and no breath. Holding at a
    coordinate's own rest costs nothing and wanders nothing (first attempt).
  - Added `common` (shared slow drive, a share of effort; Astra's model
    used 0.002 × noise). HYPOTHESIS: force-steadiness CV of 1–3% puts it
    near 0.02–0.03.
  - With noise 0.2 and common 0.03, the held vowels have F0 sd 8–15 cents,
    jitter 0.2–0.8% and shimmer 1–7%.
- Tissue strain relaxation inside the sheet
  (`SheetTube::set_strain_relaxation`; Astra's TissueRelaxation, per
  node): 10% of the foundation relaxes toward the present deformation at a
  log-normal wandering rate. It is dissipative and adds no force of its own.
- Muscle control of the sheet: `set_band_elasticity` (stretch, propagated
  as the variation is) and gliding `set_band_offset` (the rest radius moves
  over one control period; the pieces are recompiled after it).
- Astra's controller overshoots about 7% of the travel at onset and settles
  over about 1.3 s (integral rate). It is kept as written for now.
- `listen/m5_closing_cord.wav`: held /a/ /i/ /u/ with the closing cord at
  three settings: −0.35 mm about 126 Hz, −0.25 mm about 125 Hz, then a
  deeper cord about 97–110 Hz.
- NEXT: the owner's verdict on m5; articulation of the ring supports from
  the score (a sheet-tube Speaker); jet noise at the orifice (aspiration);
  the holding overshoot.

### The slit cord, heavy walls, and the sheet speaker (2026-10-13)

Owner, on m5: "sounds more like a anemic chinese boy now."

What was wrong, and what changed. Measurements are /a/ unless noted.

**Child-sized formants.**
- OBSERVED: F1 was 830 Hz, against the area function's 668 Hz (lossless chain-matrix model, `/tmp/klres.py`).
- Story's area function gives F3 ≈ 2950 Hz at any resolution. The 2440 Hz I compared against belonged to another speaker.
- The extra F1 came from walls too light to carry the tissue behind them. At three times the vowel tubes' wall masses (elastic 30, support 21, mass 45, lips 36), the transfer probe gives /a/ 720/1020/2960, /i/ 310/2140/2930 and /u/ 300/770/2360 Hz.
- At ten times the formant structure broke.

**Anemic: the source.**
- The ring cord stayed open about 85% of each cycle.
- Its stiffer, more closed settings blew out: subglottal pressure inflated a widening ring 3 cm, a balloon (its pressed surface grows with r while its hoop stiffness falls).
- `cord_length` makes the cord a slit between two 14 mm folds:
  - the half-gap replaces the radius;
  - the pressure acts on the fixed medial faces, with no hoop curvature;
  - the flow area is 2L(gap − contact), with slit Poiseuille loss.
- The probe shows a mucosal wave: the folds open from below first, the opening travels up, and closure also starts from below.
- `cord_strain` (cubic outward stiffening) is kept as an option; at 1 mm it stopped the blowout but raised F0 to 145 Hz.
- The slit's around-the-ring tissue path is one sample. With its arc (3.5 mm at tissue speed), a delayed self-loop at every cord node echoed, and the cord stopped voicing.

**A cord for every vowel.**
- `/tmp/allvowel_search.py` required steady voicing on /a/, /i/ and /u/ with heavy walls at 800 Pa. 11 of 360 cords passed.
- Chosen: 3 mm thick, mass 4, elasticity 0.8, support 0.1, adduction −0.1 mm.
  - F0 127–129 Hz on all three vowels.
  - Closed 26–29% of the cycle.
  - MFDR 190–320 cm³/s/ms.
- Pitch comes from stretch: 1.0–2.0 gives 127–172 Hz (F0 ∝ stretch^0.43), independent of pressure from 600 to 1000 Pa. Below 1.0 voicing is not steady.

**Holding.**
- *Overshoot.* Astra's request, (1 − g)·target, assumed the relaxing memory already sat at the target, so every move overshot about 7% and returned over seconds.
  - An internal model now asks for target − g·(memory expected at the nominal rate).
  - An efference copy predicts the commanded trajectory through the same 20 ms sensing filter and 30 ms delay. Feedback acts only on what the prediction missed.
  - Breath now overshoots 3%.
- *Speed.* The twitch filter (second order, 17 ms) made half a move take 41 ms.
  - A burst passes the request through ((1 + Ts)/(1 + bs))² with b = 6 ms: the inverse of the muscle's own lag, clipped by the pools' capacity.
  - Half a move now takes 23 ms.
  - Effort spans were widened (3 root units for constrictions, 0.6 cm for the glottis) so that a full gesture's burst fits the capacity.
- *Timing.* Commands act 20 ms before their events (EMG precedes movement). Oral articulators move at twice the score's rate.
  - At the score's rate /s/ reached its 0.08 cm² channel only as it ended. Area goes as the square of the root, so the last percent of travel sets a channel.

**Air.**
- `jet_loss`: Borda–Carnot loss, ρ|U|U/2·(1/A₁ − 1/A₂)², where a jet expands, plus Poiseuille loss. It is a series resistance at each axial path's end, solved inside the junction, with its value from the last step's flow.
- Series pressure sources place a jet's turbulence. The glottal one enters the orifice equation.
- `lung_resistance` (1.5e5 Pa·s/m³, a 20 Hz corner): a lung branch in parallel with Astra's matched inlet.
  - Her inlet alone lost 160 Pa per 100 cm³/s, so an open glottis drained the pressure below it.
  - Now /s/ holds about 400 Pa below the glottis at 330–440 cm³/s.

**The sheet speaker** (`tract/sheet_speaker.hpp`, `tools/sheet_speak`).
- Reads the front end's score.
- Every part is a held coordinate of a SomaticHold.
- Each millisecond:
  - the front end's geometry turns the held posture into an area function;
  - every tract section's ring moves to the narrowest area over its span (taking the area at the centre missed /s/'s 3–4 mm channel);
  - the glottis part sets the slit cord's adduction, pitch its stretch, and the lungs the breath.
- Side tubes join through `SheetJoint`, using delayed air and tissue paths as Astra's system does:
  - the nose (11.7 cm, the front end's nasal areas) at the velum (9.3 cm), opened by the velum part, radiating at the nostrils;
  - Astra's maxillary sinus off the nasal cavity;
  - Astra's under-tongue pocket at 13.5 cm, opened by the floor part.
- Before the nose, every nasal was a silent stop.
- Turbulence uses the front end's source in SI units: frication × ½ρv³√A once Re > 1800, shaped at a Strouhal corner of 0.2.
- Runs at about 8× real time with the four tubes at 48 kHz.
- `listen/m7_sheet_speech.wav`: "Hello, my name is Papa. I speak with my body." and "The quick brown fox jumps over the lazy dog."

### The goat (2026-10-13)

Owner, on m7: "this is illegible goat making goat sounds ... way worse than
the boy. the boy was in turn better than the balloon." The ranking so far:
balloon (m3/m4) < boy (m5, held vowels) > goat (m7, speech).

- OBSERVED: the cord had been tuned on sustained vowels only, and nothing measured speech.
  - `/tmp/vowel_voicing.py`: the share of 10 ms frames in the score's vowels that are periodic (autocorrelation ≥ 0.35).
  - `/tmp/ltas_diff.py`: the long-term spectrum of voiced frames against the owner's `ref/user_nato.wav`, as mean |dB| over 160 Hz–10 kHz.
- Baseline:
  - The front end's own engine: 89% voiced, 2.0 dB.
  - The sheet speaker as m7 was rendered: 67–72% voiced, 4.4 dB.
  - The sheet voice was 9.5 dB too strong at 1.6–2.5 kHz, 7–10 dB too strong above 4 kHz, and 4 dB weak at 160–300 Hz: thin and buzzing.
- Cause: the source closed far too abruptly.
  - With the lung branch the clamp-contact slit cord reached an MFDR of 2000–3000 cm³/s/ms (a human's is 200–600), with ±2.5 kPa water-hammer spikes below the glottis.
  - In "Papa" the /p/ devoicing held the glottis 0.13 cm open well into the vowel. That is a whispered vowel, and the pitch tracker read 400 Hz from its noise. A faster glottis did not help, because the slit cord stops voicing once held 0.3 mm open.
- New design keys:
  - `cord_taper`: the gap grows toward the back, so the open area is integrated along the folds, and closure zips from the front with the area falling quadratically.
  - `cord_contact`: elastic fold contact (a spring cubic beyond 50 µm, plus contact damping) in place of the clamp. A pressed rest is now real medial compression.
  - `cord_loss`: the cord's own damping. Lower damping widened the voicing range from 0.1 to 0.3 mm but roughened it, so it stays at the material's value.
- With a 0.3 mm taper and contact ×1: 70–78% voiced, 3.1 dB.
  - Aspiration 0.3 changed little.
  - Remaining: +3 to +6 dB above 4 kHz and −3 dB at 160–600 Hz.
- NATO contours: 7.16 dB. The front end's own engine scored 5.4–5.8; on these measures it is still ahead.
- Deeper cord adopted: 4 mm, mass 3, elasticity 0.4, support 0.2. It voices at 101–104 Hz on every vowel, and stretch 0.9–2.2 gives 98–145 Hz.
- Turbulence sources sit on the first path past the constriction. Inside it, in series with the expansion loss, a tenfold gain changed nothing. Oscar's /s/ at 2.5–7 kHz now measures −16 dB at frication 1.5, against NATO's −14 to −21.
- `listen/m8_sheet_speech.wav`: the same two sentences as m7, with these changes.

### Still a goat: the cord is intermittently chaotic in speech (2026-10-13)

Owner, on m8: "still a goat" ... "this is just really really realy really
terrible. its just just a buzz, its like he's talking through a fan".

- OBSERVED: a held /a/ through the full speaker (`/tmp/es/held_a.score`) measures jitter 2–3%, shimmer 11–15% and HNR 0–2 dB.
- The flow trace (wander off) is a regular 8.5 ms oscillation that closes. Every 0.1–0.2 s closure fails for 40–60 ms, with the minimum flow climbing to 130 cm³/s, and then recovers. The pressure below the folds stays steady. This is intermittency: the chop the owner hears.
- No single setting removes it. With the held /a/ through the speaker:
  - Side tubes, tissue relaxation or jet loss switched off: HNR 0.1–1.8 dB.
  - Motor-pool wander off: HNR 2.6 dB.
  - Parallel clamped folds: HNR 1.0 dB (8.2 dB only at a halved 62 Hz).
  - No lung branch: HNR 1.2–2.0 dB.
  - Cord damping 0.15–0.6: the folds stop, or jump to 322 Hz.
  - Pressing the folds 0.3–1.0 mm: worse, with period doubling and 81–322 Hz regimes.
- The calibration had warned of this. On /a/ this cord was steady at stretch 0.9, 1.0, 1.5 and 2.2, and broke up at 1.2 and 1.8. Speech glides through the islands.
- The earlier "steady" cords were tested in vowel tubes without the lung branch and without the speaker's tract. That result did not transfer.
- OPEN: the direction is the owner's call. Either keep searching the sheet cord over many speaking conditions at once, or excite the sheet tube with a flow-driven fold model of proven regularity.

### The triangular glottis in the sheet (2026-10-13)

Owner's choice: "Proven folds in sheet". The sheet tube keeps:
- its tract, holding, nose, sinus and under-tongue feedback;
- the lungs, jet losses and turbulence.

Only the folds' motion comes from the triangular glottis (`tract/Cord`, VocalTractLab-style: blown into vibration by the flow it passes, with contact zipping from the front).

- `SheetTube::attach_glottis`: the junctions below and above the cord section are the glottis's incident waves (P/2) on lines of their impedances, converted to CGS. It steps twice per 48 kHz sheet step (its own 96 kHz), and its mean flow enters the junctions. The sheet's cord section stays as tissue, under the mean of the pressures around it.
- The ear (`SheetSpeaker::hear_pitch`, the segment-module body's): at each closure the folds' tension is corrected by 0.08·log(intended/heard), within ±0.3. Open-loop, this glottis gave 141 Hz for an intended 120.
- Held /a/ through the speaker:

  | Wander (noise / common) | Jitter | Shimmer | HNR | F0 sd |
  |---|---|---|---|---|
  | off | 0.10% | 0.39% | 25 dB | 2 cents |
  | 0.1 / 0.015 | 0.25% | 1.0% | 18 dB | 7 cents |
  | 0.2 / 0.03 | 0.84% | 3.6% | 8 dB | 22 cents |

  The sheet's own slit cord gave 2–3%, 11–15% and 0–2 dB. Adopted: 0.1 / 0.02.
- Speech: 79–85% of vowel frames voiced (front end engine 89%; m8 70–78%). Long-term spectrum 2.3 dB from the owner's (front end engine 2.0; m8 3.1). Remaining: −5 dB at 160–300 Hz, +4 dB at 1.6–2.5 kHz.
- NATO contours: 7.61 dB (m8 7.16). The new loss is the low band in voiceless consonants (Oscar, X-ray, Victor, Delta): 10–20 dB too much below 700 Hz while the glottis is open. NEXT.
- `listen/m9_sheet_glottis.wav`: the same two sentences.

### Intelligibility, the fine grid and the nose (2026-10-13)

- `voice/tools/render_arctic.py` and `nato_compare.py` take `ENGINE=sheet`, with `SHEET` holding the sheet speaker's arguments.
- Four-view phone identification (ARCTIC a1–a60 against bdl; bdl's held-out sentences 54.9%):
  - front end's engine: 35.7%;
  - sheet speaker (m9 settings): 17.2%.
- **96 kHz, 16 modules.** The tract's sections become the front end's 3.6 mm oral grid; the triangular glottis runs natively. Voicing 92–96% (front 89%), long-term spectrum 1.6 dB from the owner's (front 2.0), 1.8× real time. Identification unchanged at 17.0%.
- **Command lead.** Identification scores the phones on the score's timeline:

  | Lead | 0.07 s | 0.04 s | 0.02 s | 0 | −0.02 s |
  |---|---|---|---|---|---|
  | Identification | 2.4% | 8.4% | 18.8% | 23.4% | 15.3% |

  Adopted 0. Articulator speed 1–3 changes little (23.2–24.3%); 1.5 adopted.
- **Closures pressed past contact** (root −0.6 for a score area of 0). Before, the tip of /n/ sealed only as /n/ ended, and it read as /d/. Now /n/ seals at once, voicing continues through it and no pressure builds behind the seal, but /n/ → /d/ still leads the confusions (77).
- **Phone spectra** (`/tmp/phone_spectrum.py`, the middle half of each instance, re the utterance):
  - bdl's /n m/ peak at 200–400 Hz (−4.8 dB). Ours are 3–14 dB under bdl through 200–1500 Hz with the 2 kHz nasal low-pass, and 11–17 dB too bright at 400–1500 Hz without it. `air_lowpass` (a per-tube one-pole on the air paths) is the knob. A three-tap filter at 96 kHz could not low-pass at speech frequencies.
  - A wider velum port changes nothing.
  - Our ARCTIC /ɑ/ is 17 dB too strong at 1.5–2.2 kHz. The planner's /ɑ/ area function gives 745/1385/2580 Hz in an ideal tube, so the sheet sits F2 a little high. OPEN.
- Current: identification 24–25%, voiced 85–99%, spectrum 1.8 dB.
- `listen/m10_sheet_fine.wav`: the same two sentences.

### m10 verdict, and the wall chop (2026-10-13)

Owner, on m10: "a bit of ringing during papa. he does sound more human.
there's some warbling. hello my name is papa. i speak through(garbled) my
body. the fit brown fox jumps(warbles) over the lazy dog(warbles)".

- **The chop.** OBSERVED in a held /f/ with frication and aspiration off: periodic sound with lines at 200 and 400 Hz, 27 dB under a held vowel.
  - The motor wander moves the rings every millisecond, but the sheet rebuilt its pieces (mass, tension, foundation, supports, damping) only at its 200 Hz control rate. Every moving wall stepped its stiffness every 5 ms, and lines at 1 kHz came from the motor ticks.
  - Fix: `control_rate` (design key; Astra's 200, the speaker's 8000) and `SheetSpeaker::shape_rate` (8000), with the held posture interpolated between the motor pools' 1 ms ticks.
  - Every audio band of the noiseless held /f/ is now at or below 20 dB (it was 29–34 at 150–1500 Hz). What is left is wall pumping below 50 Hz, cut by a 40 Hz high-pass in `sheet_speak`, as a recording chain does.
- **The exit's kinetic loss.** Under `jet_loss` there is a series resistance ρ|U|/2A² ahead of the radiation load at the outlet. Before it, a lip channel built no pressure: /f/ had 120 Pa behind it and a 14 m/s jet, where 870 Pa was due. Now about 300 Pa and 22 m/s.
- **Jet placement** is the front end's: turbulence grows over 0.4 cm and decays over 1.2 cm past the exit, ×4 where a coherent jet strikes the incisors, with a labiodental weight.
- Nasal low-pass 4.5 kHz. Frication 3.
- Identification 26.6%, voiced 87–96%, spectrum 1.7 dB.
- OPEN, "Papa":
  - The /p/ release bursts at 85–88 dB, as loud as the vowel (human: 10–20 dB under).
  - The vowel after it sits at about 72 dB, breathy while the glottis returns.
- `listen/m11_sheet_smooth.wav`.
- "Papa" ringing:
  - The m10 spectrogram (`/tmp/es/papa_spec.png`) shows horizontal lines at 2, 3, 4 and 5 kHz through the aspirated "pa": harmonics of the motor pools' 1 ms tick, stepping the walls. m11 interpolates between ticks and the lines are gone.
  - Burst levels (`/tmp/burst_level.py`: the loudest 10 ms in a stop, re the following vowel) match bdl within about 2 dB for every stop (/p/ −5.3 against −4.2).
  - The breathy "pa" is the score's own: the glottis stays open (0.14) until 20 ms into the vowel. Ours then adducts in about 11 ms. Glottis speed ×3–6 changes nothing.
- /k p/ middles: ours are the silent closure, while bdl's labels include burst and aspiration. Aspiration 1–2.5 does not change them (identification 26.8 / 25.8%). This is the planner's timeline more than the sheet's sound; the front end's engine scores /k/ 13%.
- Nasal murmur: a bigger sinus chamber (radius 12–16 mm, against Astra's 8) changes nothing (26.2–26.6%). Ours stays 5–9 dB too strong at 400–1000 Hz and weak at 200–400 Hz. OPEN.
- Against the owner's NATO reading (`REF=owner`): sheet 10.44 dB, front end engine 5.16.
  - Many of our words end in flat DTW tails. Our audio is trimmed at the last phone's end plus 30 ms while the voice is still fading. The lungs fall over the score's 125 ms and voicing stops about 60 ms past the last phone (in "Echo").
  - The owner's words are 630–830 ms; ours are 370–630.
  - Not resolved: the front end shares the plan and scores 5.16.
- Word ends (`/tmp/es`, all 26 NATO words): ours reached −30 dB re peak a median 49 ms after the last phone; the front end's engine reaches it at the phone's end (median −5 ms).
  - `breath_release_speed` (a falling breath target moves ×3) brings "Echo" to its end.
  - Against the owner: 10.44 → 10.20 (×6: 10.16, but voicing 80%). Adopted ×3.
  - The tails were not the main cause; the gap to the front end's 5.16 is spread over all three bands. OPEN.
- REJECTED as default (kept as `--recesses 1`): the front end's piriform recesses (blind pockets 1.8 and 2.2 cm deep, 0.9 → 0.4 cm², at 2.2 cm) and its lateral pocket (a 0.5 cm neck into 2.5 cm³ at 16 cm, opened by the `lateral` part), as sheet tubes on joints.
  - Spectrum 1.7 → 2.7 dB: −10 dB at 2.5–4 kHz, the 2.2 cm pocket's zero cutting near 3.9 kHz too deep.
  - /l/ 16% → 7%; identification 26.1 → 24.5%.
- Word onsets (first 10 ms within 20 dB of the peak, re the first phone): ours median +65 ms, the front end's +30 ms. The voice started late; vowel-initial words +30–40 ms against 0–10.
  - `laryngeal_lead` commands the breath ahead of the articulators.

  | Lead | Against the owner | Identification |
  |---|---|---|
  | none | 10.20 | 26.1% |
  | breath 30 ms | 9.30 | 25.5% |
  | breath 60 ms | 9.55 | 24.6% |
  | breath, glottis and pitch 30 ms | 9.12 | 25.0% |
  | breath, glottis and pitch 50 ms | 9.05 | 21.3% |

  Leading the glottis moves voicing boundaries against the consonants. Adopted: breath 30 ms.
- /sh/ and /ch/ with the current defaults: 33% and 19% (near 0 before the exit loss and the front end's jet placement). Our /sh/ still peaks too high (3–7 kHz, −3 to −5 dB re the utterance) and is 10 dB weak at 1.5–3 kHz.
  - REJECTED: opening the under-tongue pocket ahead of the tip (15.2 or 15.8 cm instead of 13.5). The pocket resonated near 1.2 kHz; /sh/ 22–28%, /ch/ 0–6%.
- **Nasal murmur against an ideal model** (`/tmp/nasal_model.py`: lossless chain matrix with the pharynx to the velum, the closed mouth branch to the seal, and the nose to the radiating nostrils; /n/ posture from the score). Peaks at 350 and 1015 Hz (−8), a zero near 700 Hz.
  - Held /n/ through the sheet (harmonic envelope): with the 4.5 kHz nasal low-pass, one broad peak at 275 Hz sliding down with no zero, and everything above 1.5 kHz below −48 dB. Without it, peak 360–480, dip 720–960, peak at 1080 Hz (−12): the ideal's structure.
  - Adopted: flat losses on the nose's paths (taps 0.96, 0.005, 0.002), no low-pass. /n/ 5% → 16% and its bands above 700 Hz within 2 dB of bdl's.
  - The murmur is about 13 dB too quiet at 200–400 Hz, and its peak sits near 500 Hz against bdl's 250–300. /m/ still 0%. OPEN.

### m12 verdict, pitch, the tube effect (2026-10-13)

Owner, on m12: "he now sounds vietnamese. i speak in my body is warbling.
the big brown fox jumps over the lazy dog. there is a bit of a tube
effect". Sentence two was heard nearly whole ("quick" → "big").

- **"Vietnamese".** Vowels missed the plan by a median 157 cents, with 113 cents of swing inside a vowel (`/tmp/pitch_error.py`); the front end's engine 79 and 80.
  - Open loop this glottis sits about 20% high in the sheet (145 Hz for 120), and its pitch falls with pressure (148 Hz at 450 Pa, 122 at 1050).
  - The ear was too weak: gain 0.08, range ±0.3.
  - Opening for a breath, the folds flickered into contact and the ear took those flickers for cycles: its correction swung from −0.16 to +0.54 and entered the next phrase at +0.14. "I" started at 178 Hz against a planned 121.
  - Now the ear listens only while phonating (glottis under 0.03 cm, lungs over 250 Pa), with gain 0.2 and range ±1.2: 74 and 45 cents, swing 82–85.
- **Nose.** Flat losses (taps 0.96) replace the low-pass; see above. Nose wall mass ×4–15 changes nothing.
- **Tube effect.** Lossier oral paths bled the steady flow: voicing 58–89%, spectrum 3.5–10 dB off. REJECTED. Wall damping ×4 (`wall_loss`) gives voiced 82–99% and spectrum 1.6 dB; LPC B1 stays 85–90 Hz, already human, so narrow bandwidths are probably not the tube effect. OPEN: what the owner hears as "tube".
- `listen/m14_sheet.wav`.
- **The room.** `sheet_speak` applies the front end's `voice::apply_room`: a faint reverberant tail (−32 dB, RT60 0.25 s) and a background (−58 dB), as real recordings are never dry (`--room 0` leaves it dry). Against the owner's NATO reading: 9.30 → 6.12 dB (front end engine 5.16); our digitally silent closures and pauses had counted against every word. Identification 25.8%.
- Nose length 14 or 17 cm (against 11.7): /n/ 13% / 9%, /m/ 0%. No gain.
  - A held /n/ (sealed, voiced) already peaks on its first harmonics (140–280 Hz) and slopes down, as bdl's murmur does.
  - The 400–700 Hz excess in the ARCTIC averages is therefore from running speech: probably short unstressed /n/s ("in", "and") whose tip is not sealed for most of the phone, so the middle is a nasalized vowel. NEXT: the seal's timing in short nasals.
- **Nasals, found.** The oral seal is made for the whole middle of almost every /n m/ (10 ARCTIC sentences). Two other faults:
  - The velum, at the score's rate (100 ms moves, starting about 27 ms ahead), was 0–0.1 cm² open as each nasal began. `velum_speed` ×3.
  - The flat nose losses (taps summing to 0.967) bled the steady flow. With the mouth sealed, pressure built and the folds stopped (glottal flow 0 through a 129 ms /m/); the murmur read 15 dB under bdl's.
  - Now: spread taps summing to about 1 (0.4, 0.35, 0.248) and nose wall damping ×8 (`nasal_wall_loss`). /n/ 26% (5–17 before), /m/ 9% (0), 200–400 Hz at bdl's level. Still about 8–16 dB too strong at 400–1000 Hz. Identification 26.0%.
- `listen/m16_sheet.wav`.
- **/f/, a port bug.** When the lips themselves are the narrowest section (the tube's last), the jet loop started past the end and sounded nothing, labiodental weight included. Now the jet sounds at the lip junction (the front end's). /f/ 27% → 92% (front 92%); identification 27.6%.
- **Oral jet corner within 1.5–5.5 kHz** (the front end's `corner_ceiling`, fitted to bdl's fricatives; ours allowed 500 Hz–12 kHz): identification 28.9%, /s/ 79%, /f/ 94%, /sh/ 33%, /ch/ 6%.
  - /sh/ is still too strong at 3–7 kHz (−1 to −4 dB re the utterance; bdl −10 to −12) and 10–14 dB weak at 1.5–3 kHz. That is the front cavity resonating too high, not the source. OPEN, with /ch/.
- /ch/ 6% (front 69%): its spectrum is /sh/'s fault, stronger. Ideal lossless response of the /sh/ area function (`/tmp/front_cavity.py`; posture from the score, source past the constriction): peaks near 1.1 and 4.8 kHz, none at the 2.5–3 kHz of bdl's /sh/. The shared mouth geometry gives /sh ch/ no long front cavity; the front end's body compensates with its sublingual space. OPEN: a /sh/ posture whose front cavity resonates near 2.7 kHz.
- /sh/ in the ideal model, with the lips extended by 0–3 sections (the score's protrusion lengthens a real tract; our fixed grid instead moves the lips inward): no clean 2.5–3 kHz peak. A 1.1 kHz resonance dominates every variant, because a source in series at a junction drives the back cavity as strongly as the front. A real sibilant mostly excites the front cavity. OPEN: the noise's coupling (perhaps a monopole at the teeth, or a source impedance that decouples the back cavity), and lip protrusion as real lengthening.
- **Releases.** Pressed past contact (`closure_press`), a closure then opened at the score's rate. /ch/ in "chair" released 30 ms after the score did, and its frication was cut to half the phone. `release_speed` ×3: an opening constriction (root rising from ≤ 0.3) moves three times faster, as the old speaker's did (tissue under pressure parts abruptly).
  - Identification: 28.9 → 29.4% at ×2 (/ch/ 6 → 25%), 29.9% at ×4.
- **/l/.** bdl's /l/ has a deep valley at 1.5–2.2 kHz (−41.5 dB re the utterance); ours had −27.
  - The planner emits lateral-pocket events only when a phone sets `pocket`, and /L/ does not. A pocket on the tuned /L/ made it worse (8–12%).
  - The planner's built-in dark /l/ (tongue body 8.25 cm / 0.25 cm², tip contact at 16.7 cm) gives /l/ 37% (17% with the front end's tuned L, which was fitted to its own body) and /r/ 43%.
  - `tract/data/phones_sheet.txt`: the front end's table without its L line. `render_arctic.py` and `nato_compare.py` pass it for `ENGINE=sheet`.
  - Identification 31.0%; NATO 5.43 dB (front end 5.4–5.8); against the owner 5.75 (front 5.16).
- `lateral_pocket` (`--lateral 1`) splits the lateral pocket from the piriform recesses; unused by default.
- `listen/m19_sheet.wav`: the two sentences planned with `phones_sheet.txt`.
- /w/ 22–29% (front 63%). The lips reach their 0.14–0.2 cm². The tongue has about 42 ms toward /u/ before the next vowel's command. Faster articulators (cap 800/s, ×2) do not help: identification 30.0–31.0%, /w/ 22–25%. OPEN: our /w/ is 8–9 dB too strong at 0.4–1.5 kHz against bdl.
- Identification over three seeds (current defaults): sheet 31.0 / 29.5 / 28.9 (mean 29.8%); front end engine 35.7 / 35.8 / 35.0 (mean 35.5%). The single-seed figures above were the favourable end. Single-seed differences under about 1.5 points are noise.
- Three-seed search (`/tmp/sheet_search.sh`) against the defaults' 29.8%:

  | Change | Mean |
  |---|---|
  | wall loss 1 | 29.1 |
  | frication 2 | 28.6 |
  | release ×2 | 29.6 |
  | press 0.3 | 29.6 |
  | aspiration 1 | 29.9 |

  All are within the noise; the defaults stand. Intelligibility is on a plateau for single-parameter changes. The remaining gap to the front end (35.5%) lies in specific phones: /sh ch/ geometry, /w/, /d k/, close vowels.
- **Voiced closures.** Voicing during /b d g/ closures (glottal-flow swings over 30 cm³/s, ARCTIC a1–a8): /b/ 12%, /d/ 60%, /g/ 38%; bdl about 61%.
  - `closure_pull` (pharyngeal rings drawn outward during a voiced closure, the segment-module speaker's): 3 mm gives 47 → 49% overall closure voicing and identification 30.2% (three seeds; noise). The heavy, supported walls barely change the cavity's volume.
  - Off by default. NEXT if pursued: release the pharyngeal walls' support during the closure (a `set_band_support` control), as the old speaker's `wall_release` did.
- `closure_release` (`SheetTube::set_band_support`: the pharyngeal rings' support released during a voiced closure) changes nothing: the foundation term (the tissue's modulus) carries the walls' stiffness, not the support. Off by default.
- The /b/ figure (12%) rests on 2 instances. One is phrase-initial, with the glottis still closing from the breath (0.027–0.040 cm), devoiced as English often devoices it. Not a sound basis for tuning.
- `tools/sheet_stream`: score lines on stdin while the tube speaks, 16-bit 48 kHz on stdout in 10 ms blocks, spoken up to the latest event heard (fixed level 0.6; 0.8 touched full scale). `say … --score /dev/stdout --plan-only --pitch 120 --phones ../tract/data/phones_sheet.txt | sheet_stream`: 4.15 s of speech in 2.85 s of wall time, planner included; 93% of vowel frames voiced.
- `tools/sheet_voices`: several sheet voices on threads, each with its own seed and pitch (PITCH_STEP semitones apart). Four together: 4.00 s of speech in 2.78 s of wall time (1.4× real time) on the Mini.

### How the brief is met now (2026-10-13)

- **Tokens in.** `voice/` parses English text into phrases, inserts breath and phrase marks, breaks words into phonemes, and plans each phoneme as postures of named parts with timing: the score. It streams (`--score /dev/stdout`).
- **The engine.** `tract/SheetSpeaker` deforms the sheet tube to satisfy the score:
  - every part is held by Astra's motor pools, so pieces wander under the effort of holding;
  - ring supports move to the area function;
  - glottis, pitch, lungs, velum, floor and lateral are driven from the score.
- **Sections.** The tract is 16 modules of three sections in sequence (elastic, support, mass) on one sheet network. They are joined by membrane tension and delayed tissue paths, and air acts on every surface.
- **Feedback.**
  - Below the cords: Astra's subglottal sections with the lungs' branch.
  - Below the tongue: the under-tongue pocket.
  - In the sinuses: Astra's maxillary sinus on the nasal tube.
  - The nose itself: radiating at the nostrils, joined at the velum.
- **The folds.** Owner's choice: the triangular glottis, blown into vibration by the flow, not driven. There is no direct mass-spring stimulation anywhere: the lungs are the only energy.
- **Real time and voices.** 1.7× real time for one voice, 1.4× for four together at 96 kHz.
- **Open:** /sh ch/ geometry, /w/, /d k/, close vowels, and the owner's verdict.
- /dh/: the tuned table closes it (area 0; the front end did this because "with" was heard right only once /dh/ closed). Pressed and released fast in the sheet, it is a full stop: /dh/ 4%. The planner's built-in /dh/ (a 0.18 cm² channel) gives /dh/ 14% but no overall change (29.5 against 29.8, three seeds). Kept closed pending the owner's ear on "speak with".
- **The sibilants' source.** The oral jet's turbulence as a local flow (monopole: the pressure strength over the local ρc/A) instead of a pressure in series. A series source drives the back cavity as strongly as the front: in the ideal model /sh/'s 1.1 kHz back resonance dominated. Adopted with frication 1 (`jet_monopole`; the lip junction stays series).

  | Source | Identification (three seeds) |
  |---|---|
  | series | 29.8 |
  | monopole, frication 0.4 | 28.3 |
  | monopole, frication 1 | 29.8 |
  | monopole, frication 3 | 30.1 |

  At frication 1: /sh/ 50% (33), /ch/ 44% (about 20), /s/ 74%; /s/ within 2–3 dB of bdl at 3–7 kHz. /sh/ is right at 2.2–3 kHz but 14 dB loud at 3–4.5 kHz.
- `obstacle_gain` (the incisor obstacle's ×(1 + g); the front end's 3): 1 or 0 barely moves /sh/ at 3–4.5 kHz (+2.1 / +1.8 dB against +2.8). The excess is the front cavity's own resonance, a little high in this geometry. Kept at 3.
- **/sh/ posture** (`phones_sheet.txt`): area 0.13 at 14.0 cm (the front end's 0.1 at 14.4).

  | Place, area | /sh/ | /ch/ | Identification |
  |---|---|---|---|
  | 14.4, 0.1 | 50% | 44% | — |
  | 13.8, 0.1 | 56% | 38% | — |
  | 14.1, 0.1 | 56% | 44% | 30.9% |
  | 14.0, 0.13 | 61% | 38% | 30.8% |

  At 14.0 / 0.13 the peak sits at 2.2–3 kHz as bdl's. The top band (4.5–7 kHz) is 9 dB weak, and /sh/ is about 10 dB loud overall.
- Three seeds with the monopole sibilants and the /sh/ posture: 30.8 / 29.5 / 30.1, mean 30.1% (29.8 before; front end 35.5).
- **/w/** (`phones_sheet.txt`): lip opening 0.07 cm² (the planner's 0.20). /w/ was 8–9 dB strong at 0.4–1.5 kHz on bdl, a lip aperture too open for its low formants. Now within 1–4 dB over 0.4–4.5 kHz: /w/ 25% → 37% (0.12) → 39% (0.07). Three seeds: 31.4 / 30.0 / 30.1, mean 30.5% (30.1 before).
- /k g/ place 12.05 or 12.6 (the table's 11.2 / 11.65): /k/ stays 5%, /g/ 0–3%. The velar burst's place is not the problem; the timeline labels are. Not adopted.
- Rounded back vowels with smaller lips (/uw/ 0.3 or 0.15, /ow/ 0.4 or 0.25, /ao/ 0.6): /uw ow ao uh/ stay at 0–9%. Not adopted. The front end's engine is also weak there (/uw/ 12%, /ow/ 0%).
- `listen/m21_sheet.wav`: the two sentences with the current table (/sh/, /w/).
- **/iy/** (`phones_sheet.txt`): channel 0.45 cm² (the table's 0.179), place and lips kept. The sheet's kinetic losses at a narrowing left the narrow channel 5–9 dB quiet on bdl.
  - /iy/ 18% → 46%; the planner's built-in /iy/ gives 37%.
  - Three seeds 31.6 / 31.3 / 31.1, mean 31.3% (30.5 before).
- /er r/ constriction 0.4 cm² (0.25): three seeds mean 31.6% against 31.3 (noise); /er r/ themselves unchanged. Not adopted.
- Wider channels elsewhere (one seed; baseline 31.6): /ao/ tongue 0.5 (32.1, /ao/ 0%), /v/ 0.2 (31.8, 0%), /z/ 0.08 (31.2, 22%), /th/ 0.12 (31.7, 0%). No phone improves; /iy/'s gain does not generalise. Not adopted.

### m21 verdict: "regression, slightly. some sounds are too faint." (2026-10-13)

- Word levels against m12 (the owner's "more human"):
  - "I" −3.0, "speak" −3.4, "body" −2.4, the second "the" −1.9 dB.
  - "speak" was the narrow /iy/, fixed after m21 (−12.0 against m12's −12.1).
  - /w/ peaked 34–42 dB under the loudest sound in "with" and "quick", in m12 too: nearly silent.
- **/w/.** After /k/ in "speak with" the glottis stayed 0.14 cm open through all of the 92 ms /w/: steady flow, no voicing, a whisper. The score's voicing onset came 46 ms into the /w/ (timed for the old engine's quick folds), and ours took 30 ms more to close.
  - `voicing_lead` 0.03 s and `voicing_speed` ×2: a glottis event closing toward phonation acts 30 ms early and twice as fast; openings keep their time (leading every glottis event had cost the consonants).
  - Identification, three seeds: 33.8 / 31.8 / 32.5, mean 32.7% (31.3 before; lead 0.05: 31.1).
  - Words against m12 now within 1.7 dB or louder ("I" −1.2, "body" −1.1, "the" −1.7; "with" +2.8, "name" +5.8).
- `listen/m22_sheet.wav`.
- Confusions with the voicing lead: /d/ → /t/ 54, /k/ → /t/ 35, /ah/ → /eh r ih er ae/ 20–33 each.
  - The planner's built-in vowel postures for /ah/, or for /ah eh ae aa/ (one seed; baseline 33.8): 33.7 / 33.9, /ah/ 28–29%. Vowel postures are not the problem. Not adopted.
- /d/ between voiced sounds (ARCTIC a1–a20, n = 11): closure voiced a median 71%, voicing resumes at the release (median 0 ms). /t/ in the same position (n = 25) is flapped by the planner: closure fully voiced, onset 0 ms. The /d/ → /t/ confusion is not voicing; intervocalic /d/ and the flap are near-identical (in bdl too). Largely a labelling matter for the test.
- **Warble check, m14 → m22** (sentence 1):

  | | Pitch error | Swing inside vowels | Loudness wobble |
  |---|---|---|---|
  | m14 | 74 cents | 82 cents | 1.28 dB per 10 ms |
  | m22 | 57 cents | 142 cents | 1.63 dB per 10 ms |

  Sentence 2: wobble 0.97 → 1.37, swing 85 → 78.
  - Not caused by the voicing lead, release speed, monopole or room.
  - The old nose (taps 0.96, wall ×1, velum ×1) together gives swing 103, wobble 1.53. No single one of the three reproduces it: taps alone 131 and 1.99, wall alone 143 and 1.59.
  - The nose couples to pitch steadiness, but on one sentence this measure is too noisy to isolate. Kept pending the owner's ear on m22.
- Per vowel (sentence 1), pitch swing grew most next to nasals:

  | Vowel | Neighbours | m14 | m22 |
  |---|---|---|---|
  | /ow/ "hello" | before /m/ | 196 | 498 cents |
  | /ey/ "name" | between /n/ and /m/ | 52 | 142 |
  | /ih/ "is" | after /m/ | 52 | 130 |
  | /ay/ "my" | after /m/ | 103 | 181 |
  | /iy/ "speak" | — | 101 | 239 |

  The livelier nose (the spread taps that keep the flow, giving real murmurs) couples back onto the folds while the velum is open. HYPOTHESIS, untested: the ear could correct more gently while the velum is open.
- **The warble's amplifier is the ear.** With it off, the near-nasal swings drop: "hello" 498 → 144 cents, "name" 142 → 95, "is" 130 → 65. Gain 0.2 had been chosen before the ear was gated to phonation.

  | Gain | Sentence 1 error / swing | Sentence 2 error / swing |
  |---|---|---|
  | 0.2 | 57 / 142 | 35 / 78 |
  | 0.12 | 57 / 106 | 57 / 85 |
  | 0.08 | 76 / 98 | 55 / 80 |
  | 0.05 | 99 / 93 | 52 / 66 |

  Adopted 0.12. `listen/m23_sheet.wav`.

### Searching for the warble and the tube (2026-10-13, after m23)

- Ear gain 0.12 (m23): identification 32.4% over three seeds, against 32.7 at 0.2. No cost.
- **Fixed colorations: REJECTED as the tube effect.**
  - Method: the long-term spectrum's fine structure after removing its 1/3-octave trend (`/tmp/ripple.py`: voiced loud frames, 300 Hz–7 kHz). A resonance or notch that stays put whatever the phone survives as a ripple.
  - Over 29 ARCTIC sentences: ours 0.81 dB RMS, bdl 0.86. Two sentences: ours 1.77, the front end 1.26, bdl 1.35.
  - No standing notch or peak sets us apart. The side tubes (`--side-tubes 0`) change nothing.
- **Too-clean harmonics: REJECTED.** Harmonic-to-noise ratio by band in steady voiced frames (`/tmp/hnr.py`):

  | Band | bdl | ours | the owner |
  |---|---|---|---|
  | < 1 kHz | 21.6 | 20.8 | 23.8 |
  | 1–2 kHz | 9.3 | 10.9 | 11.9 |
  | 2–3 kHz | 6.2 | 7.3 | 8.2 |
  | 3–4 kHz | 4.4 | 4.7 | 6.5 |
  | 4–7 kHz | 1.0 | 0.9 | 1.7 |

- **Spectral motion.** Spectral flux and band-envelope modulation (`/tmp/flux.py`, 29 sentences) are near bdl's. The one difference: ours has 0.13 of its modulation energy in 8–16 Hz, against 0.10.
- **Wobble in steady vowels: near human.** Inside vowels of at least 120 ms, after a quadratic detrend (`/tmp/vw.py`):
  - pitch 6.0 cents (three seeds 5.9–6.2), bdl 5.4;
  - level 0.32 dB, bdl 0.24.
  - So the warble is not a tremor of the hold. Hold noise 0, common drive 0 and the ear off do not reduce it.
- **The warble's source: the glottis's pitch falls with breath pressure.**
  - Held /a/, ear off, dry (`/tmp/fp.sh`): 150 Hz at 390 Pa, 145 at 650, 130 at 845, 122 at 1040. Held /a/ and held /n/ show the same.
  - A ramp shows the fall continuous from about 280 Pa (171 Hz) to about 830 (125 Hz), with the closed share rising from 0 to 0.44, then flat at the commanded pitch. The fall lies in speech's own range.
  - People rise by about 2–6 Hz per cmH2O. VocalTractLab's triangular glottis has no pressure term in its tension; this fall is its own.
  - In running speech this gives the scoop seen in the trace (`/tmp/es/ptrace_v1.png`). Every phrase onset starts sharp at low pressure and slides down: "my" starts at 160 Hz and falls to 125 within 70 ms. Consonants that cut the transglottal pressure push pitch up.
- REJECTED: lateral strain stiffening of the edge springs, k(1 + eta x^2) (Ishizaka–Flanagan's cubic). Eta 100–10000 raises pitch at every pressure, but the fall persists (eta 1000: 180 → 189 → 174 Hz), and eta 3000 and above period-doubles.
- **Elongation by deflection** (`Cord::set_elongation`). The fold's tension Q (in this model pitch ∝ Q) is multiplied by (1 + gamma d)/(1 + gamma d_ref), floored at 0.9. Here d is the lower edge's peak squared deflection, held and released over 10 ms, and d_ref = 2e-3 cm^2 (about 650 Pa). A swinging fold is stretched along its length, as a plucked string is.
  - The cycle-averaged deflection (5 ms mean of x^2) failed: when the folds close for 44% of the cycle, the average falls although the swing is wider, and pitch dropped at 845 Pa.
  - Peak follower, held /a/ from 390 to 1040 Pa:
    - gamma 200: 121 → 140 Hz;
    - gamma 300: 114 → 146 Hz (4.7 Hz per cmH2O);
    - gamma 400: 109 → 148 Hz.
  - A 30 ms release: pitch error 80 cents against 70. Without the floor, phrase-final vowels lost voicing ("hello" /ow/ 80 → 45% voiced) as the dying breath let the factor fall to about 0.63.
  - Over 29 ARCTIC sentences (`/tmp/pe29.py`: plan against heard, every 15 ms in the vowels):

    | Setting | Voiced | Pitch error | Swing | Onset error |
    |---|---|---|---|---|
    | none (three seeds) | 86–88% | 77–78 | 118–127 | 138–139 |
    | gamma 200 | 89.2% | 69 | 113 | 111 |
    | gamma 300 | 89.9% | 70 | 106 | 112 |

    Errors are in cents.
  - Identification: 32.3% (gamma 200) and 31.2–31.4% (gamma 300), against 32.4. Adopted: gamma 200.
  - The trace with elongation (`/tmp/es/ptrace_v1g.png`): "my" starts at 107 Hz and rises to 115 Hz. The onset scoop and the breath blip are gone.
  - The aggregate 8–20 Hz pitch wobble over whole voiced runs did not move (about 120 cents against 114). That measure mixes in consonant effects; the plan-referenced figures above are the cleaner ones.
- **Pitch after voiceless consonants** (`/tmp/onset.py`: the first 30 ms of a vowel after a pause or voiceless consonant, against the vowel's own median).
  - bdl: first frame +218 cents, 0–15 ms +125, 15–30 ms +20. This is the known raised F0 after voiceless obstruents: the folds are stiffened for devoicing and relax afterwards.
  - Ours: −6 / −10 without elongation, −32 / −19 with it.
  - `devoicing_raise`: tension × exp(raise) while the glottis is open for a voiceless consonant (0.03–0.3; breaths excluded), relaxing with `devoicing_relax` once it closes. The ear's intended pitch carries the same factor (efference), so the ear does not correct the raise away.

    | Raise, relax | First frame | 0–15 ms | 15–30 ms | Pitch error |
    |---|---|---|---|---|
    | 0.08, 20 ms | −2 | −5 | −10 | 71 |
    | 0.15, 20 ms | +57 | +20 | −9 | 72 |
    | 0.25, 20 ms | +118 | +54 | +1 | 72 |
    | 0.15, 40 ms | +56 | +24 | +2 | 74 |

    The first three columns are in cents. Voicing rises to 91–91.5%.
- Identification over three seeds:
  - gamma 200 with the floor: 31.5%;
  - with the devoicing raise 0.15 / 40 ms: 31.7%;
  - with 0.25 / 40 ms: 31.8%.
  - The current defaults before these changes: 32.4%.
  - The elongation costs about a point. It is on every seed and borderline against the noise. It is kept for the onsets.
- Adopted: devoicing raise 0.25, relax 40 ms. Onsets after voiceless sounds now read +131 / +82 / +15 cents (bdl +218 / +125 / +20). The plan-referenced pitch error rises to 83 cents because the plan has no microprosody. `listen/m24_sheet.wav`.
- **Faint sounds** (`/tmp/phone_level.py`: each phone's middle 60%, dB re the sentence's loudest 20 ms, median, 29 sentences, against bdl):

  | Phones | Too quiet by |
  |---|---|
  | /b/ | 22 dB |
  | /th/ (5 instances) | 21 dB |
  | /w/ | 15 dB |
  | /g/ | 13 dB |
  | /f/ | 9 dB |
  | /iy/ | 8.6 dB |
  | /dh p y/ | about 8 dB |
  | /hh/ | 7 dB |
  | /r ey ih/ | about 6 dB |

  /aa ae ao sh ch/ are +0 to +2.4 dB.
  - The jet loss (`--jet-loss 0`) takes /w/ from −14.5 to −8.6 dB; nothing else moves.
- **Close vowels.**
  - Running-speech /iy/ formants match bdl's (305 / 2205 / 2936 Hz against 322 / 2325 / 2802), so the posture is right; the energy is not.
  - H1−H2: bdl −3.2 dB in open vowels and −6.6 in close vowels (F1 lifts H2); ours +3.5 and +3.3. Our close vowels' F1 is too weak.
  - (`/tmp/es/held_iy.score` with tongue 13.75 / 0.45 is not the planner's /iy/: F2 1337 Hz. With 0.3 it gives F2 2279 Hz.)
- **The sheet's transfer against the ideal.** `tools/sheet_areas` prints the sheet's cross-sections at a moment of a score. `/tmp/sheet_transfer.py` compares the radiated output over the glottal flow, at the harmonics, with a plane-wave chain over the same areas.
  - /aa/ matches the ideal within 1 dB below 1 kHz.
  - /iy/ is 19 dB short at F1 and 20 dB short at F2/F3, relative to the valley between them: the valley is filled.
  - Not the cause: nostrils muted (`--listen 1`), side tubes off, jet loss off, adduction (`--voicing-lower` 0.002 to −0.002), wall loss 1–32, wall mass ×100 (`--tract-mass`).
  - **The cause is the oral air taps' steady loss:** Astra's 0.985 / 0.010 / 0.003 sum to 0.998, about 0.2% per 4.4 mm hop. With taps summing to 1 (0.987 / 0.010 / 0.003):
    - /iy/ comes within 1–5 dB of the ideal; with wall loss 1–8, within 1–3 dB.
    - Held /iy/ F1 moves 220 → 262 Hz, B1 173 → 29 Hz; /aa/ B1 91 → 47 Hz. With wall loss 8: /iy/ B1 49 Hz, /aa/ 43 Hz.
- REJECTED (for now): oral taps summing to 1 in running speech.
  - The same taps made the trachea lossless. Psub swung ±1500 Pa, with a −6000 Pa spike at the first adduction, and vowel voicing fell to 73–75%.
  - `SheetDesign::subglottal_taps` / `--subglottal-taps` now lets the airways below the folds keep their own taps. With lossier ones (0.95 / 0.03 / 0.01 to 0.8 / 0.1 / 0.05):
    - voicing 82–92%;
    - the spectrum went dark above 2.5 kHz (−8 to −14 dB) with the heavier ones;
    - /iy/ re the sentence moved only −8.6 → −7.1 to −9.7 dB, as open vowels gained as much.
  - Half the steady loss (oral sums 0.999 and 0.9995, trachea unchanged): voicing 80–84%, swing 127–137 cents, /iy/ −8.5. No gain.
  - With the sheet's own glottal spectrum (about −13 dB/octave from 120 to 800 Hz), even the ideal tube puts held /iy/ about 10 dB under /aa/. People show about 5 in running speech (bdl 4.8). OPEN: the remaining close-vowel deficit: the source's low harmonics against people's, and B1.
- **Voiced closures are silent because nothing radiates, not only because voicing stops.**
  - In "big" (/b/ at 3.98 s), the folds vibrate about 35 ms into the closure, then stop as the mouth's pressure climbs to 285 Pa.
  - `closure_expansion` (`--expansion`, cm^3/s): behind a voiced oral closure with the velum shut, a volume sink at 5 and 11 cm (the larynx lowering, the cheeks giving). It begins over 10 ms and ends over 15 ms.
    - At 60 cm^3/s the folds vibrate through the whole closure.
    - Phone levels did not move: /b/ −23.6 → −23.1 to −24.2 dB. The sealed tube's only outlets, mouth and nostrils, are shut. /d/ fell to −12 to −15 dB as the burst lost pressure.
  - `skin_radiation` (`--skin`): the tract's volume velocity above the folds (from the walls' moving radii) radiated like the lips' flow, through a two-pole 400 Hz low-pass (`skin_corner`: skin and tissue are heavy).
    - Without the low-pass the walls' fast ripples put +16 and +27 dB into 4–10 kHz.
  - Results over 29 sentences (dB re bdl; spectrum = mean |LTAS difference|):

    | Setting | /b/ | /g/ | /d/ | /w/ | /iy/ | Spectrum | Vowels voiced |
    |---|---|---|---|---|---|---|---|
    | none | −23.6 | −13.2 | −5.2 | −14.0 | −8.3 | 2.2 dB | 90% |
    | skin 1 | −4.6 | −8.6 | −2.7 | −9.2 | −6.1 | 2.1 dB (160–300 Hz −7 → −5) | 85.5% |
    | skin 1, expansion 60 | −4.0 | −2.0 | −5.2 | −9.1 | −7.3 | 2.0 dB | 86% |
- Identification over three seeds:

  | Setting | Mean |
  |---|---|
  | skin 1 | 31.5% |
  | skin 1, expansion 30 | 31.3% |
  | skin 1, expansion 60 | 29.8% |

  The m24 defaults: 31.8%. Adopted: skin 1 (400 Hz corner), expansion off.
  - The skin component alone (renders with and without at fixed gain differ only by it, since it does not feed back) sits a median 9 dB under the total in voiced frames (90th percentile −1 dB).
  - Real neck radiation during vowels is probably weaker. The sheet's walls pulse strongly at H1–H2, and that compliance may be part of the rubbery "tube" sound. OPEN.
  - Skin 0.5 with a 200 Hz corner: 13.8 dB under, but /b/ only −11 dB.
- `listen/m25_sheet.wav`.
- **Formant bandwidths in running speech** (`/tmp/vformants.py`, LPC in the middle 40%, 29 sentences): ours are as wide as bdl's or wider. /iy/ B1 106 Hz against 42, /ih/ 141 against 82, /uw/ 102 against 52. Narrow, ringing resonances are not the tube effect. REJECTED.
- **Vowel targets were off.** With the front end's postures, against bdl in running speech (F1 / F2, Hz):
  - /ah/ 667 / 1838 against 577 / 1248 (nearly /eh/);
  - /ao/ 668 / 1079 against 433 / 850;
  - /aa/ 812 / 1511 against 691 / 1215;
  - /ow/, /uw/ similar errors.
  - Mean F1/F2 distance 4.07 semitones (`/tmp/vtune.py`: plans and renders a0001–a0060 with a given table, 34 s a trial).
- `/tmp/vopt.py`, four rounds of coordinate search: 1.77 semitones. Identification 33.0% (33.5 / 32.7 / 32.7), against 31.5. The new postures are in `data/phones_sheet.txt`. `voice/tools/render_arctic.py` takes `PHONES=` for another table.
- Second pass of `/tmp/vopt.py` (finer steps; IH EY UH ER added): 1.49 semitones. AH 2.8 → 0.3, AO 2.1 → 0.6, IH 1.2 → 0.3, ER 1.6 → 0.5; OW stays 2.9 (a diphthong); UH 6.6 (2 instances).
  - Identification 33.3% (33.7 / 32.6 / 33.5).
  - AH's lips are now 0.6 cm^2, AA's 1.0: narrower than a person's. In this tract they stand in for a longer front cavity.
  - `listen/m26_sheet.wav`.
- **Wall compliance as the tube effect (HYPOTHESIS, under test).** The skin's share in voiced frames at tract wall mass ×1 / ×3 / ×10: −9.4 / −12.3 / −14.3 dB re the mouth (a0001).
  - The 29 sentences re-planned with the new vowel table (`/tmp/es/rs`; the old set is in `rs_old`):

    | Wall mass | H1−H2 (bdl −3.2) | Spectrum | /b/ |
    |---|---|---|---|
    | ×1 | +4.0 dB | 1.9 dB | −7.1 dB |
    | ×3 | +2.9 dB | 1.6 dB | −11.2 dB |
    | ×10 | +2.5 dB | 1.7 dB | −11.3 dB |

  - Identification 33.3% (×1) against 32.7% (×3): within noise.
  - For the owner's ear: `listen/m26_heavy_walls.wav` (×3) beside `listen/m26_sheet.wav`.
- **Where identification fails now** (seed 1, 33.7%; `ROWS=` in `voice/tools/phone_test.py` prints what given phones were heard as):

  | Phone | Accuracy | Heard mostly as |
  |---|---|---|
  | /dh/ | 2% | /b/ (33) |
  | /v/ | 4% | /n/, /r/ |
  | /z/ | 14% | /s/ |
  | /ng/ | 0% | /r/ |
  | /m/ | 18% | /r/, /n/ |
  | /n/ | 31% | /r/, /m/, /w/ |
  | /d/ | 18% | /t/ (44) |
  | /ow/ | 0% | /aa/ |
  | /l/ | 18% | /w/ |

  - Band spectra (`/tmp/phone_spectrum.py`): our /v n m dh/ are +7 to +19 dB at 400–1000 Hz against bdl. A vowel-like F1 sounds through these consonants.
  - Mid-/n/ the mouth is sealed (0.000 cm^2) and the nostrils carry it all, but 700–1000 Hz is still about 8 dB strong.
- REJECTED, all on one seed against 33.7% overall:

  | Change | Overall | Phones |
  |---|---|---|
  | DH area 0 → 0.15, V 0.12 → 0.2, Z 0.05 → 0.1 | 32.4% | dh 4, v 0, z 13 |
  | `--voiced-open` 1.5 (new: the opening gain below a 0.03 setting; 0.7 before) | 32.9% | — |
  | `--voiced-open` 2.5 | 33.4% | — |
  | nasal wall loss 24 | 33.4% | — |
  | nasal taps 0.3 / 0.4 / 0.3 | 32.8% | n 20 |
  | both nasal changes | 33.2% | — |

  The voiced-open gains changed /z v dh/ by nothing.
- **Prosody** (`/tmp/prosody.py`, 60 sentences):
  - F0 SD per sentence 2.77 semitones against bdl's 2.33; 5–95% range 6.7 against 5.7. Pitch movement is human-sized.
  - Periodic share inside phones (autocorrelation at least 0.5, 30 ms):

    | | Vowels | Nasals | Liquids | Voiced stops | Voiced fricatives | Voiceless |
    |---|---|---|---|---|---|---|
    | bdl | 0.76 | 0.81 | 0.68 | 0.48 | 0.34 | 0.10 |
    | ours | 0.75 | 0.73 | 0.61 | 0.13 | 0.27 | 0.03 |

  - Vowel-middle periodicity: 0.63 against 0.69.
  - The voicing gap is in the voiced stops (hence /d/ → /t/).
- REJECTED: the closure expansion with a volume budget (`--closure-volume`, cm^3 per closure; the walls return with tau 100 ms between closures).
  - One seed: 60 / 1.5: 32.7%; 100 / 2.5: 32.9% (d 18 → 34%, b 37 → 50%, g 17 → 7%); 60 / 3: 32.3%.
  - Three seeds, 100 / 2.5: 31.9% against 33.3. Voicing the closure helps /d b/ and costs the other stops.
- The 400–1000 Hz excess in /m n v/ is not the skin: without it (`--skin 0`) the band spectra are the same within 1–2 dB.
- Nor is it the shape (`/tmp/midareas.py`: `sheet_areas` at phone middles). /v/'s lips are at 0.09–0.2 cm^2 and /m/ is sealed (0.000) at every middle checked. The excess is acoustic (the sonorant /v/; the nose's response above 700 Hz). OPEN.
- REJECTED: narrower /v/ and /dh/.
  - One seed: V 0.05 with DH 0.05: 33.6% (dh 2 → 16%, th 22 → 11%, v 4%). V 0.03: 33.7% (v 7%).
  - Three seeds, DH 0.05 alone: 32.7% against 33.3.
- Single changes to the phone table now sit within noise of 33.3%. The remaining consonant gaps (/v z dh/ voiced frication, nasal murmur above 700 Hz, closure voicing) look structural, not posture settings.
- **Nasals.**
  - The oral seal is in place: /m/ sealed at 25 / 50 / 75% of all 20 instances checked; /n/ at 50 and 75% of 30, and at 25% of 63% (`/tmp/sealfrac.py`). The 400–1000 Hz excess is the nose's response.
  - `--nasal-lowpass 15000` (per hop): /n/ 700–1000 Hz −13.2 → −18.1 dB (bdl −28.7). 400–700 Hz unchanged at +7 dB. Identification 33.6% one seed (n 32%, m 20%) against 33.7. `--nasal-lowpass 25000`: 33.1%.
  - Sinus radius 11 or 14 mm (the neck's Helmholtz zero near 600 Hz moving lower): /n m/ band spectra within 1 dB of the 8 mm sinus. Its coupling, through the 0.1 cm^2 neck at 2/5 of the nose, barely touches the murmur.
  - OPEN: the murmur's 400–700 Hz excess. A real nose's paired maxillary and other sinuses, and its divided, mucosa-lined passages, damp N2 heavily; ours is one passage with one small sinus.
- **NATO with the m26 defaults:** against the owner 5.12 dB (5.75 at m19; front end 5.16); against the NATO reference 6.02 (5.43 at m19).
  - Against the owner, the high band carries the most (6–9 dB in Quebec, Golf, Foxtrot, Victor, Mike).
  - Our isolated words are a median of about 540 ms against the owner's 730.
  - `say --rate` 0.85 / 0.75 lengthens them (Papa 540 → 620 / 700 ms) but scores 5.26 / 5.35. The DTW-aligned measure does not reward duration. Whether the owner hears our words as rushed is OPEN.
- **Clicks from the skin.** In the high-resolution spectrogram of the NATO words beside the owner's (`/tmp/es/quebec_cmp.png`, `voice/tools/spectrogram2.py`), ours showed one-sample clicks at word ends.
  - At the end of "Bravo" the peak was 11× the vowel's level, 19 dB over speech; 11 dB without the skin (`/tmp/clicks.sh`).
  - The skin radiated the tract's whole volume change, including the articulators' movement stepped at the shape rate.
  - Now only the walls' vibration about their supported rest radii (`SheetTube::rest_radius`, which glides with the articulation) radiates. Bravo 10.9 dB; the skin is −10.2 dB re the mouth in voiced frames (90th percentile −5.6, against −0.5).
  - Over 29 sentences: /b/ −10.4 dB, /g/ −8.8, /iy/ −5.7, spectrum 1.8 dB.
  - **Identification 34.3%** (35.0 / 33.3 / 34.6), against 33.3. The front end's engine is at 35.5.
  - NATO: against the owner 5.40 (5.12), against the reference 5.44 (6.02).
  - Silence between words: −56 to −58 dB re speech; the owner's −47 to −55.
- The spectrogram also shows our words far shorter than the owner's ("Bravo" about 0.35 s of sound against 0.75). Our vowels are also dark above about 4.5 kHz, where the owner's carry formants to 7 kHz. OPEN.
- `listen/m27_sheet.wav`.
- Vowel band spectra against bdl (`/tmp/phone_spectrum.py`, 60 sentences, dB re the utterance):
  - /aa/ is 7 dB short at 4.5–7 kHz and 5 dB at 3–4.5 kHz.
  - /ow/ is 9 dB short at 4.5–7 kHz and 11.5 dB at 3–4.5 kHz.
  - /ae eh iy/ are within ±3 dB above 4.5 kHz.
  - The dark highs are in the vowels whose lips the search narrowed (AA 1.0, OW 0.36 cm^2). The trade is formant position against brightness. OPEN.
  - All vowels are 2–4 dB short at 200–400 Hz and 3–4 dB strong at 80–200.
- **H1−H2, a candidate for the hollow, hooty "tube" quality.** Loud voiced frames:

  | Source | H1−H2 |
  |---|---|
  | the owner (`ref/user_nato.wav`) | −3.9 dB |
  | bdl | −2.7 dB |
  | the NATO reference | −11.5 dB |
  | ours (NATO words) | +3.4 dB |

  Our first harmonic is about 7 dB too strong against the second, in vowels of every kind.
  - Not the skin (`--skin 0`: +3.5 open / +2.1 close), nor the walls (×10 mass as well: +2.8). The ear and lossier subglottal taps make it worse.
- **The glottal pulse is the cause** (`/tmp/pulse.sh`, held /a/, 650 Pa):
  - open quotient 0.65; flow H1−H2 +13.1 dB; output +6.5.
  - The folds open and close in a clean triangle (peak 0.105 cm^2, closed 35%).
  - The flow lags into a flat top of about 150 cm^3/s (`/tmp/es/pulse.png`).
  - Bernoulli through 0.1 cm^2 at 600 Pa would pass about 320 cm^3/s. The plateau implies about 130 Pa left across the glottis at peak opening: the pressure collapses below or above it. Psub traces swing 300–1000 Pa.
- Glottis settings against H1−H2:

  | Setting | Open quotient | Output H1−H2 | F0 |
  |---|---|---|---|
  | default | 0.65 | +6.5 dB | 120 Hz |
  | voicing lower 0 | 0.69 | +7.7 dB | — |
  | voicing lower −0.004 | 0.99 (never under 0.002 cm^2) | +1.9 dB | 142 Hz |
  | voicing upper −0.016 | — | +6.3 dB | — |
  | `--contact 2` | — | +6.5 dB | — |
  | lung resistance 0 | — | +7.0 dB | — |
  | subglottal taps lossier | — | +7.5 dB | — |
  | subglottal taps lossless | — | +6.1 dB | — |
  | lungs ×1.5 | 0.56 | +3.7 dB | — |

- REJECTED: louder breath in running speech. Lungs ×1.3: H1−H2 +2.9, spectrum 4.3 dB (highs +5 to +7 dB), voicing 84.6%. Lungs ×1.5: H1−H2 +2.6, spectrum 5.7 dB (highs +8 to +11), voicing 80%.
- OPEN, NEXT: where the transglottal pressure goes during the open phase. Candidates: the narrow epilarynx (0.455 cm^2 over 2.5 cm above the folds: a large inertance), the subglottal sheet's swing. A trace of the pressures just below and above the folds over a cycle would show it.
- Mini disk: 330 MB free. `voice/out` holds 5.2 GB (arctic 2.7 GB rewritten in place each run; `evaluate_*`, `arctic_rms` and `four_view_tune` are older outputs, not deleted).
- **The flat-topped pulse comes from a choked airway below the folds.** Trace of the glottis's own inlet and outlet pressures (new trace columns `below`, `above`; `Cord::pressure_below/above`), held /a/ at 636 Pa (`/tmp/es/press.png`):
  - while the folds are open, the pressure below sags to about 300 Pa;
  - at closure it spikes to about 1400 Pa while the pressure above dips to −750.
  - Astra's subglottal bands carry support offsets that narrow the airway to 0.25–0.39 cm^2 (`sheet_areas`); a person's trachea is 2–3 cm^2.
- `subglottal_narrowing` (`--subglottal-narrowing`) scales those offsets (1: Astra's; 0: the bands' own 1.5–3.8 cm^2). Held /a/:

  | Narrowing | Open quotient | Flow H1−H2 | Output H1−H2 |
  |---|---|---|---|
  | 1 | 0.65 | +13.1 dB | +6.5 dB |
  | 0.7 | 0.56 | +8.9 dB | +2.3 dB |
  | 0.4 | 0.54 | — | +1.8 dB |
  | 0 | 0.54 | +8.3 dB | +1.7 dB |

- Over 29 sentences:

  | Setting | H1−H2 (open) | Pitch error | Spectrum | Identification (three seeds) |
  |---|---|---|---|---|
  | defaults | +4.2 dB | 89 cents | 1.8 dB | 34.3% |
  | narrowing 0 | +1.1 dB | 70 cents | 3.6 dB | — |
  | narrowing 0.4 | +1.3 dB | 68 cents | 3.2 dB (highs +3 to +6) | 31.9% |
  | narrowing 0.4, lungs ×0.8 | +1.7 dB | 91 cents | **1.5 dB** | 32.4% |

  Narrowing 0 and 0.4 give onset 56–57 cents. The defaults stay. The ear decides: `listen/m27_wide_airway.wav` (narrowing 0.4, lungs ×0.8) beside `listen/m27_sheet.wav`.
- Narrowing 0.7 with lungs ×0.9: H1−H2 +2.0 dB (open), +2.2 (close); spectrum 1.6 dB; voicing 86.8%; pitch error 86 cents.
  - Identification 33.6% (34.7 / 33.2 / 32.8) against 34.3: within noise.
  - Adopted as the defaults (`subglottal_narrowing` 0.7, `lung_scale` 0.9). `listen/m28_sheet.wav`.
- Elongation reference after the airway change, over 29 sentences:

  | Reference | Pitch error | Voicing |
  |---|---|---|
  | 2.0e-3 | 86 cents | 86.8% |
  | 1.6e-3 | 80 cents | 88.7% |
  | 1.3e-3 | 79 cents | — |
  | 1.0e-3 | 83 cents | — |

  Identification 33.4% at 1.6e-3, against 33.6. Adopted 1.6e-3.
  - Pitch swing within vowels is 140–148 cents since the airway change (106–116 before). OPEN.
- **Steady-vowel wobble** (`/tmp/vw.py`) on the re-planned set (`/tmp/es/rs`, planned by `vtune.py` without a seed):

  | Setting | Wobble |
  |---|---|
  | current | 10.8 cents |
  | elongation 0 | 15.3 |
  | devoicing raise 0 | 13.3 |
  | skin 0 | 10.9 |
  | all of today's changes off, old airway | 10.7 |
  | the original vowel table (`/tmp/es/phones_orig.txt`) | 11.4 |
  | the old airway | 12.5 |

  The earlier 5.9–6.2 cents came from the previous set, planned by `render_arctic` with its seed. The rise is in the material: the planner's seed and contours, with only 22–35 long vowels in a run. Today's engine changes did not cause it; the elongation lowers it.
- `split_nose` (`--split-nose 1`): two passages behind the port (55% / 45% of the cross-sections, the second 6% longer), each with its own sinus (the second's chamber ×1.25). The port's opening is shared 55 / 45.
  - Running-speech /n/ 700–1000 Hz: −13.2 → −15.2 dB (bdl −28.7); /m/: −13.4 → −15.6 (bdl −26.0).
  - 400–700 Hz unchanged (−8.3 / −8.2 against bdl −15.6 / −17.1).
  - With heavy oral walls (×10): /n/ 400–700 −10.0. With lossless oral taps: −11.7.
  - The held /n/ already has the ideal model's peak, dip and peak (above). The shallow fall above the murmur's peak belongs to running speech's short nasals. Off by default. OPEN.
- REJECTED: `--fastest 800` (articulator rate cap ×2) for the late /n/ seal. One seed: 33.8% against 34.9; n 38% (39), m 20% (18); the /n/ band spectrum is unchanged.

### Owner: "m2i was ironically our best so far ... everything on your new engine warbles like im using a vibrator on their throat" (2026-10-13)

- m27, m28 and m27_wide_airway: "all three of these are warbling and sound exactly the same". `listen/m2i_new_engine.wav` (the segment-module engine, frication 1.5, before the 12 kHz ring fix) is the owner's best. The sheet has fewer defects but has lost the voice.
- **10–40 Hz envelope modulation** (`/tmp/vibrator.py`, band power over loud stretches):

  | Voice | 10–20 Hz | 20–40 Hz |
  |---|---|---|
  | m28 | −1.4 | −5.8 |
  | m2i | −8.8 | −11.0 |
  | bdl | −9.6 | −12.3 |
  | owner | −8.4 | −11.4 |

  - Confounded by stretch length (m28 has a few long stretches that carry articulation). Median windowed 10–40 Hz modulation (100 ms windows): m28 0.117, m2i 0.109, bdl 0.107, owner 0.081, m12 0.142. Close.
- It falls with the pressure under the folds:
  - lung resistance 0: −6.7 / −8.3; restored to the same pressure (lungs ×1.3): −1.2 / −2.9 again;
  - lungs ×0.7: −5.9 / −7.6;
  - the lung corner (20–1000 Hz) changes nothing.
- A held /a/ is not modulated (−8 to −12 dB at 10–20 Hz).
- True fold irregularity from the trace (`/tmp/truejit.py`, closures from the glottal area): jitter 2.1%, shimmer 0.42 dB, no alternation (period doubling). Lungs ×0.6: 1.2% / 0.28. Holds and tissue off: 1.8% / 0.24. Mildly above a person; not obviously a vibrator.
- No measure here clearly finds what the owner hears. A/B set for the owner's ear (`listen/vibrator_ab/`, sentence 1): current; walls silent; stiff walls; no muscle wander; softer breath; no ear; no lung branch; old airway; flat pitch command; everything off.
- **The segment-module engine's flutter came with m2k.** 20–40 / 40–70 Hz envelope modulation over the m2 sentence set:

  | Version | 20–40 Hz | 40–70 Hz |
  |---|---|---|
  | m2c | −8.9 | −11.9 |
  | m2e | −11.5 | −16.1 |
  | m2g | −9.6 | −14.6 |
  | m2h | −11.0 | −15.2 |
  | m2i | −11.0 | −15.2 |
  | m2j | −11.7 | −15.8 |
  | **m2k** | **−4.0** | **−6.8** |
  | `speak` today on m2k's scores | −4.0 | −6.8 |

  - m2k lowered the breath from 750 to 600 Pa.
  - The same scores with the lungs ×1.25 (750 Pa): −10.9 / −14.6, m2i's level.
- `listen/m30_segment_750.wav`: today's segment-module engine at 750 Pa on m2i's four sentences. It is m2i with the ring and comb fixes, the flap taps and /ae/; the candidate return to the owner's preferred voice.
- Owner on m30: "its warbling less", then: "it still sounds like a mid early 2000s tts model. not even a good one. not even the voice labs we were working to improve from".
- Segment engine with lungs ×1.25 / ×1.4 / ×1.55: 20–40 Hz modulation −10.9 / −10.7 / −9.9 dB. No gain above 750 Pa.
- The old front-end engine (`/tmp/es/listen/old_1..4`) against m30 (segment, 750 Pa) and m28 (sheet), the same four sentences:

  | | old | m30 | m28 |
  |---|---|---|---|
  | 20–40 / 40–70 Hz modulation | −12.9 / −19.4 | −10.9 / −14.6 | −5.8 / −11.2 |
  | spectrum against bdl | 3.4 dB | 2.6 dB | 2.6 dB |
  | spectrum against the owner | 1.6 dB | 2.4 dB | 1.0 dB |

  - m30 is bright: 2.2–3 kHz +5 dB and 4.5–7 kHz +6 to +7 against both references.
  - HNR at 3–4 kHz: m30 8.5 dB, old 5.0, m28 4.8, bdl 3.8. m30 is too clean and buzzy in the upper bands, a classic synthetic trait (HYPOTHESIS for the "2000s TTS").
- The F0 tracker (`/tmp/prosody.py`) puts 29–30% of every m2-series file at its 302 Hz floor (the shortest lag). This is a tracking artifact on these files: the F0 SD and range figures for them are invalid.

## Milestone 4 — see what is wrong (2026-10-13)

Owner: VocalTractLab "sounds natural. what you have made doesnt even beat microsoft Whistler TTS". Then: keep going without checking in. The measures used so far "do absolutely nothing useful to identify what is wrong with the sound ... you may have to write wholly new transforms". The goal is in `GOAL.md`.

- **VTL's examples**, downloaded with permission to `scratchpad/vtl_dl/audio`, mirrored at `/tmp/es/vtl`:
  - copy-synthesis pairs `vtl2.3-example0[1-6]-orig/-synth` (German, 22 kHz);
  - expressive variants 07–09;
  - speakers `vtl2.4-example-{JD2,M01,W02}-synth`;
  - fricative sequences.
  - The download in `~/Downloads/VTL2.4` holds only Windows binaries and is not run.
- **Cycle-level transforms, first try** (`/tmp/hear.py`: band-envelope periodicity, cycle repeatability, 3–7 kHz crest, cycle-domain sidebands). REJECTED as steering: VTL's originals, bdl and the owner differ among themselves as much as from our engines (recording, speaker).
- **Motion statistics** (`/tmp/motion.py`: LPC F0/F1/F2 speeds every 5 ms) separate only weakly: our F2 curvature is 1800–2300 st/s^2 against 2400–4000 in every real or natural source. LPC tracking noise dominates.
- **What the eye sees** (four-view enhanced time-frequency images, `voice/build/tf_analyze` + `voice/tools/tf_view.py`; `/tmp/es/tfv/a.png`, `b.png`):
  - VTL's synthesis and bdl: harmonics curve as pitch glides; formants sweep in continuous arcs from sound to sound.
  - Ours (m30 segment, m28 sheet): blocks. Flat harmonics (one pitch level per syllable) and formant bands held still through each sound, then jumping. Steady states joined by quick transitions is the look of early-2000s TTS.
  - m28's harmonics also look kinked.
- The scores say why:
  - pitch is one target per syllable with 50–100 ms moves (a staircase);
  - the tongue moves in 70–150 ms and then holds;
  - bdl's F0 slopes within every syllable (accents rise to a peak and fall 20–40 Hz over 100–200 ms; `/tmp/es/f0cmp.png`).
- Prototypes, score post-processes for the segment engine at 750 Pa:
  - `/tmp/stretch.py`: tongue body, jaw and the lips' non-closing moves take ×2.5 longer, centred. Closures and tip, dorsum, velum, glottis and breath keep their timing. The formants now arc continuously through "my name is" (`/tmp/es/tfv/c.png`). `listen/m31_flowing.wav`.
  - `/tmp/contour.py`: sloped syllable pitch targets (after Xu's target approximation). Stressed: rise from 100 ms before the vowel to +30% at 40% into it, then fall to −30% by the next vowel. Unstressed: −6%. The harmonics now arc as VTL's do (`/tmp/es/tfv/d.png`). `listen/m32_intonation.wav` (with the stretch).
- Still visible: thin steady lines at 3–7 kHz in ours (clean stationary high harmonics, buzz) where VTL's high band is diffuse. NEXT.
- **The planner, from text.** `PlanOptions::accent_peak` / `accent_fall` / `unstressed_drift` (`say --accent-peak` ...) emit the sloped syllable pitch targets. 0 keeps the old single target. Defaults: peak 0.3, fall 0.3, drift 0.06.
  - `vowel_reach` 0.15 → 0.35 s and `vowel_arrival` 0.5 → 0.8: the tongue keeps travelling through each vowel.
  - Both are visible as arcs in `/tmp/es/tfv/e.png` (p0 / p1). The planner is shared, so the sheet engine gets them too.
  - `/tmp/four.sh` plans the four listening sentences from text (pressure 750) and speaks them through the segment engine.
- **Brightness: slices of the loudest vowels** (`/tmp/slice.py`: the five loudest 60 ms frames, narrowband, dB re the peak below 1.5 kHz; `/tmp/es/slices.png`):
  - bdl, the owner, VTL's original and VTL's synthesis: formant peaks with deep valleys; 1.5–2.5 kHz peaks −30 to −44 dB; 4–8 kHz peaks −41 to −55, floors −60 to −99.
  - Ours (segment, p1): a flat envelope from 0.5 to 3.7 kHz (peaks −9 to −12) and a 4–8 kHz floor near −48.
  - Turning off aspiration, jet flux, closing turbulence, frication or wall variation changes nothing: the brightness is harmonic. The source is too sharp.
  - Fold postures (VOICEPOST 0.008/0.0, 0.012/0.004) barely move it.
- **The fold edges' damping** (`speak` env `CORDDAMP=lower,upper`; VocalTractLab's 0.1, 0.6). Peaks:

  | Damping | 1.5–2.5 kHz | 2.5–4 kHz | 4–6 kHz | 6–8 kHz |
  |---|---|---|---|---|
  | 0.1 / 0.6 | −12 | −9 | −35 | −37 |
  | 0.2 / 0.7 | −13 | −14 | −32 | −43 |
  | 0.3 / 0.6 | −11 | −23 | −46 | −49 |
  | **0.3 / 0.8** | **−41** | **−35** | **−49** | **−45** |
  | 0.4 / 0.8 | −32 | −24 | −47 | −40 |

  - Per-vowel voicing unchanged; F0 5–10% higher at 0.3 / 0.8.
  - The formants stand apart with dark gaps between them (`/tmp/es/tfv/f.png`).
  - Default in `Speaker` now 0.3 / 0.8. `listen/m34_softer_folds.wav`: from text, the planner's new motion and intonation, these folds.
  - The sheet engine's Cord still has 0.1 / 0.6. NEXT: check it too.
- **The comb.** Zoomed wideband images (`voice/tools/spectrogram2.py`; `/tmp/es/zoom_cons.png`): the vowel of our "hello" shows evenly spaced horizontal bands (about 2.5, 3.3, 4.2, 5, 6 kHz), a uniform tube's resonances ringing through.
  - In bdl ("fighter in my life") the upper formants fade into diffuse texture, and the formant bands stand apart with black holes between them. Ours was a mush from 1 to 5 kHz.
  - The segment tube's viscous loss is frequency-flat (the boundary layer's value at 2250 Hz), and it has no heat conduction.
  - `Tube::set_wave_smoothing(k)`: a one-pole low-pass on every hop, y = (1 − k) x + k y'. No loss at DC, a loss growing with frequency.
  - k 0.03 / 0.06 / 0.1 fade the 5–7.5 kHz bands progressively (`/tmp/es/zoom_hf.png`). It delays each hop by about k samples: every formant about k lower.
  - Default in `Speaker` 0.05 (`speak` env `HFLOSS`). `listen/m35_damped_tube.wav`.
- The frame's displaced air (`Tube::frame_flow_`) now passes two poles at 250 Hz. An 18 dB haze around 1 kHz in a pause went; the absolute line levels did not change.
- **Cord glides in two stages** (posture, rest offset, tension; each stage half the old time constant): command velocity no longer jumps at every 1 ms tick. A pause's small tones went (+10 → +3.7 dB).
  - The 1/2/3 kHz "lines" at 1.34–1.42 s were the /m/ murmur's 9th, 18th and 27th harmonics (F0 about 111 Hz), not an artifact. `/tmp/tones.py` cannot separate the two in voiced audio.
- Open in the images: the /ow/ of "hello" holds steady formants (no diphthong glide); sonorant consonants (/l/, /m n/) leave weak landmarks.
- **The owner's own reading of sentence 1** is `ref/user_three.wav` (0.74–1.24 "hello", 1.50–2.98 "my name is papa", 3.86–5.78 "and I speak with a body"; 21.7 s, eleven phrases). Against it (`/tmp/es/own_vs_ours.png`):
  - ours is about 35% faster;
  - ours is brighter in vowels;
  - our stop bursts are loud full-band columns where the owner's are soft.
  - `/tmp/tilt.py` (vowel 2.5–5 kHz re 0–1.2 kHz; burst 2 ms peaks re the median voiced frame):

    | | Vowel tilt | Bursts |
    |---|---|---|
    | owner | −22 dB | −20 to −21 dB |
    | bdl | −17 to −19 | −18 to −22 |
    | VTL synth | −24 | −24 |
    | m35 | −18 to −20 | −14 to −18 |
    | m2i | −14 | −12 |

    m2i, the owner's favourite, is the brightest and burstiest, so neither decides "the voice".
- `PlanOptions::phrase_fade_from` (`--phrase-fade-from`): the breath fades from this share of the last vowel (0.45 → 0.7). "hello," /ow/ voiced 68–72 → 79–82%.
- Accent peak +30% on top of the planned stressed pitch put NATO's stressed vowels at 200–250 Hz (the plan already holds the accent). Now peak 0.05, fall 0.25: pitch 98–158 Hz, "Foxtrot" /aa/ voiced 45 → 90%. The arcs are smaller (`/tmp/es/tfv/j.png`).
- `listen/m36_from_text.wav`: all of the above, from text.
- Open: short unstressed vowels voiced 0% ("hello"'s /ax/, "papa"'s first /ah/, "and"'s /ax/).
- `speak` env `RELEASE` (closure opening speed ×1, ×0.5, ×0.3): fewer bursts reach the detector (8.3 → 6.2 /s), but their 2 ms peaks stay about −16 dB re voiced frames. Not the lever for loud bursts (bdl −18 to −22, the owner −20).
- The owner's "hello" begins with a breathy, partly voiced /h/ (first 60 ms −16 dB re the word, centroid 500 Hz); ours is a quiet voiceless hiss (−28 dB, centroid 1060 Hz, 12% above 3 kHz). Minor. OPEN.
- `voice/tools/render_arctic.py` `ENGINE=tract`: the segment engine (say at 750 Pa, `speak` with the seed; `PHONES` optional). `/tmp/seg_search.sh`: its three-seed identification per say-option set.
- **Identification, segment engine with the new planner** (`/tmp/seg_search.sh`, three seeds): **35.7%** (34.9 / 36.4 / 35.8), the best of our engines (the front end's 35.5). The same engine with the old planner settings (accent 0, reach 0.15, arrival 0.5, fade 0.45): 33.7%.
- Timing against the owner's sentence 1:

  | | Ours | Owner |
  |---|---|---|
  | "hello" | 0.51 s | 0.50 s |
  | "my name is papa" | 1.12 s | 1.48 s |
  | "and I speak with a body" | 1.26 s | 1.92 s |
  | pause before the last phrase | 0.22 s | 0.88 s |

  Our durations follow bdl's read-speech pace; the owner reads more slowly. Left as is pending the owner's ear.
- **Fricatives on the segment engine** (`/tmp/phone_spectrum.py`; `ENGINE=tract`, 60 sentences): at frication 1.5, /s/ was +7 dB at 4.5–7 kHz and +10 at 2.2–3 kHz against bdl; /sh/ +3 to +13 at 1.5–3 kHz; /f/ +13 at its low end.
  - At 0.8: /s/ −6.6 (bdl −9) at 4.5–7 kHz, /sh/ within about 3 dB, /f/'s top bands −37 (−32 to −37).
  - `BodyLayout::frication` default 1.5 → 0.8 (`speak` env `FRICATION`). /f/'s low bands are still about +13. OPEN.
  - `render_arctic.py` takes `SPEAK_ENV="VAR=value ..."` for speak.
- **A table mix-up.** `say` loads `data/phones_en.txt` only when run from `voice/`. The four listening sentences (`/tmp/four.sh`, run from `/tmp/es`) used the planner's built-in vowels; the ARCTIC renders used `phones_en.txt`. Both now use `tract/data/phones_tract.txt` (render_arctic's `ENGINE=tract` default; four.sh `--phones`).
- Segment engine vowels with `phones_tract.txt` (`ENGINE=tract /tmp/vtune.py`): mean 3.99 semitones.
  - /ao/ 7.0 (F1 639 / bdl 433)
  - /iy/ 6.2 (F1 226 / 323)
  - /ah/ 5.0 (F2 1658 / 1248)
  - /uw/ 4.5
  - /ow/ 4.1
  - Search running (`START=tract ENGINE=tract /tmp/vopt.py`).
- Segment vowel search: 3.99 → 1.67 semitones (AE's lips capped at 10 cm^2; the search drove them to 17.6). Identification, three seeds:
  - `phones_tract.txt` as it was, with frication 0.8: 36.3%;
  - tuned: **36.7%** (36.3 / 37.4 / 36.4).
  - Adopted into `data/phones_tract.txt`. `listen/m37_vowels.wav`.
- Segment engine phone levels re bdl (`/tmp/phone_level.py`, 60 sentences, tuned table): nearly all within ±5 dB (the sheet: /b/ −22, /w/ −15). Quiet: /th/ −7.6, /hh/ −5.8, /iy/ −4.3. Loud: /v/ +5.2, /t/ +4.5, /k/ +4.4, /dh/ +3.8.
- `PlanOptions::aspirate_gap` 0.09 → 0.06 (the /h/ fold gap). One seed:

  | Gap | Overall | /hh/ |
  |---|---|---|
  | 0.09 | 36.3% | 22% |
  | 0.06 | 36.8% | 49% |
  | 0.035 | 36.7% | 53% (/iy/ 27%) |

  The owner's /h/ is breathy and partly voiced.
- **Bursts on today's engine** (`/tmp/bursts2.py` on the 40-sentence hold-out corpus re-spoken, `/tmp/burstvar.sh`): P / T / K peak −13 / −13 / −15 dB re the loudest vowel (bdl −21 / −18 / −20), rising 41–47 dB per 5 ms (bdl 29–34). /k/ is flat across bands (−7 −7 −7 −8; bdl −2 −6 −6 −14): heard as /t/.
  - RELEASE ×0.5 or frication 0.4 move them partly, but frication is shared with the fricatives.
  - `BodyLayout::jet_attack` (`speak` env `JETATTACK`): the jet noise builds with this time constant and falls at once (turbulence takes milliseconds to develop).

    | Attack | Rise P / T / K (dB per 5 ms) | Peak P / T / K | /k/ shape |
    |---|---|---|---|
    | 0 | 41 / 47 / 41 | −13 / −13 / −15 | −7 −7 −7 −8 |
    | 2 ms | 36 / 45 / 35 | — | — |
    | 4 ms | 34 / 43 / 32 | — | — |
    | 8 ms | 31 / 43 / 31 | −16 / −16 / −16 | 0 −6 −8 −9 |

  - Adopted 8 ms. /t/ still rises 43: the release's pressure transient, not the noise. OPEN.
  - Fricatives with it (60 sentences): /s/ and /sh/ within about 3 dB of bdl; /f/'s low bands no longer +13 (now −4); /z/ 4–6 dB weak above 3 kHz.
- `listen/m38_softer_bursts.wav`.
- Identification with all of the above (three seeds): **37.5%** (37.3 / 37.4 / 37.9). The segment engine is now the main line; the sheet engine has not had today's fixes (frication, the jet attack, the vowel table). The planner's changes reach both.
- Against the owner's reading in four views (`/tmp/es/tfv/k.png`):
  - the owner's pitch is steadier and lower through "my name is papa";
  - our /p/ releases still show as bright full-height lines;
  - the owner's pauses sit on a room and microphone floor where ours fell to digital silence.
  - `speak` now applies `voice::apply_room` (the sheet tool's faint tail and background; `ROOM=0` dry): pauses −49.7 dB re speech (the owner's −47 to −55). `listen/m39_room.wav`.
- Real time with today's engine (Mini):
  - `stream` (`say ... --score /dev/stdout | stream`): 6.94 s of speech in 1.69 s, planner included (4.1×);
  - `voices`: four voices together 4.3× real time, eight 2.7×.
  - `stream` does not apply the room (apply_room works on a whole buffer). OPEN.
- Segment nasals against bdl (60 sentences): /m n ng/ are about 9 dB strong at 80–200 Hz, weak at 1–2.2 kHz, and 9–11 dB strong at 3–7 kHz. /m/ and /n/ have nearly the same spectra (bdl's differ by 4–5 dB at 1–1.5 and 2.2–3 kHz).
  - Identification /m/ 2%, /ng/ 0–5%, /n/ 59%. The same under the old vowel reach (`seg_try`), so not the new motion.
  - `speak` env `NOSEHF`: the hop low-pass on the nose (k 0.1 / 0.2 / 0.3: 4.5–7 kHz −40 / −42 / −44; bdl −46).
  - Default in `Speaker`: nose 0.2. 2.2–3 kHz stays 3–5 dB strong (not the nose). The low-band excess and /m/ against /n/ are OPEN.
- Nose hop loss 0.2: identification **38.9%** (38.9 / 39.0 / 38.9), from 37.5.
- **H1−H2 again.** m39 H1−H2 +5.0 dB (the owner −3.5, m2i −0.7, m28 +1.7): the damping that removed the buzz rounded the pulse. `/tmp/srcsweep.sh` (the four sentences; median H1−H2 and vowel tilt):

  | Damping (lower, upper) | Posture | H1−H2 | Tilt |
  |---|---|---|---|
  | 0.1, 0.6 | — | −0.3 | −13.2 |
  | 0.3, 0.8 | — | +4.8 | −23.3 |
  | 0.2, 0.8 | — | +1.0 | −20.1 |
  | 0.25, 0.8 | — | +2.8 | −22.0 |
  | 0.15, 0.8 | — | +0.1 | −16.5 |
  | 0.2, 1.0 | — | +3.3 | −19.4 |
  | 0.2, 0.8 | upper −0.012 | +0.6 | −20.0 |

  The lower edge's damping trades brightness against H1. The owner: −3.5 / −22; bdl about −2.7 / −18.
  - Adopted 0.2 / 0.8 with `voicing_upper_` −0.012. Identification **39.6%** (40.0 / 38.6 / 40.1). `listen/m40_balanced_source.wav`.
- m40 against the owner: spectrum 1.6 dB (bdl 2.2). Vowels voiced 0%: "hello"'s /ax/, "papa"'s first /ah/, "and"'s /ax/, and "Echo"'s stressed /eh/.
  - "Echo": the breath started 40 ms before the phrase (`breath_lead`, tuned for consonant-initial NATO words) and rose over 110 ms while the folds were shut. Vibration began weak at 3.47 s and reached full flow only as the folds opened for /k/.
  - `PlanOptions::vowel_breath_lead` 0.12 s for phrases that open on a vowel (`--vowel-breath-lead`): "and" 0 → 100%, "Echo" 0 → 100%, "Alfa" 53 → 90%.
  - Identification 39.8% (40.1 / 38.7 / 40.5). `listen/m41_vowel_onsets.wav`.
  - Still unvoiced: "hello"'s /ax/ after the breathy /h/, and "papa"'s first vowel between aspirated /p/s (the planner stresses "papa" pa-PA).
- Per-phone (seed 1, 40.1%): v / ng / uw / ow 0%, g 3%, uh 9%, er 10%, th 11%, m 14%, k 15%, l 16%. /uw/ and /ow/ are heard as /ah/; k → t 28; l → r / w; m → n 25.
  - "papa" is P AH2 P AA1 by the built-in lexicon (NATO's pa-PAH) though cmudict has P AA1 P AH2. Deliberate for NATO; left.
- **/uw/ back and rounded.** The table's fronted /uw/ (13.28 / 3.25 / 0.6) gave F2 1456 (bdl 1239). F2 by posture:

  | Place, area, lips | F1 / F2 | Distance to bdl |
  |---|---|---|
  | 8 / 0.4 / 0.3 | 466 / 974 | 7.0 st |
  | 9 / 0.5 / 0.35 | — / 1097 | — |
  | 10 / 0.6 / 0.3 | 379 / 1242 | 2.1 st |
  | 10.5 / 0.45 / 0.25 | — | — |

  At 10.5 / 0.45 / 0.25: /uw/ 18%, overall 40.6% (one seed). Adopted.
- /ow/ glides to a rounded /U/ (gplace 9.5, garea 0.4, glips 0.25; the built-in kUH has lips 1.0): /ow/ 0 → 4%, 40.8% (one seed). Its mid-vowel F1/F2 already match bdl (1.1 st). Its trouble is elsewhere. OPEN.
- Three seeds with both: **40.0%** (40.8 / 38.9 / 40.2). `listen/m42_rounded_vowels.wav`.
- REJECTED (one seed; baseline 40.8%):
  - the planner's built-in dark /l/ (`L` line removed): 40.5%, /l/ 8%;
  - velars further back (K 12.2, G 12.4): 40.5%, /k/ 15%.
- `tools/talk.sh "text" OUT.wav [seed]`: the main line end to end (say at 750 Pa with `phones_tract.txt`, then speak with its room). 5.65 s of speech in 1.35 s on the Mini. `listen/m42_new_text.wav`: "Good evening. I hope the voice sounds a little more human today."
- NATO (`nato_compare.py`, `ENGINE=tract`, `--pressure 750`): against the owner **5.59** dB, against the reference **6.01** (the front end 5.16; the sheet 5.12–5.40 earlier). The DTW contour score has disagreed with the owner's ear before (the m22–m28 sheet scored well and "warbles"). Recorded, not steered by.
- `tools/talk_remote.sh "text" [OUT.wav] [seed]`: from the laptop, runs `talk.sh` on the Mini over ssh (`MINI` host, default m4mini), copies the WAV back and plays it with afplay (`PLAY=0` skips). `~/ducklips` exists only on the Mini.

### Owner, on `talk_remote.sh` "i hope i sound a little more human today" (2026-10-14)

"he does. still sickly, warbling a bit. evening and human came through not necessarily legible. the warbling is the biggest problem to tame".

- **Diagnostics added to `speak`.**
  - `CYCLES=path`: one line per glottal closure, with time, period, peak glottal flow, the sound's peak-to-peak, the intended pitch, the ear's correction and the lung pressure.
  - Switches: `BREATHVAR`, `TENSIONVAR`, `POSTUREVAR`, `PITCHFB`, `BREATHTREMOR`, `WANDER` (`Speaker::set_wander_gain`), `SKIN`, `WALLK` / `WALLLOSS` (vocal wall support), `SINUSNECK`.
- **A validated cycle transform.** `analysis/wmatch.py` finds cycle periods from audio by waveform matching (pyin voicing, then the next cycle maximising the normalized correlation of one period with the next). It reports:
  - RAP: the 3-cycle relative average perturbation, which removes intonation slopes;
  - jumps: cycles with RAP > 3.3% per 100 cycles.
  - Validated against the engine's true closures (`CYCLES`) on three renders: RAP 0.77 / 0.40 / 0.31% against 0.82 / 0.40 / 0.31%, jumps 13.0 / 9.3 / 7.7 against 13.0 / 12.2 / 7.8.
- **MEASURED (RAP median %, p90 %, jumps / 100):**

  | Voice | RAP median | RAP p90 | Jumps |
  |---|---|---|---|
  | owner, three sentences | 0.42 | 2.0 | 4.6 |
  | owner, NATO | 0.23 | 1.5 | 3.9 |
  | bdl a0001–a0008 | 0.30–0.59 | 3.1–4.2 | 9–15 |
  | VTL originals | 0.36–1.04 | 1.6–4.6 | 3–16 |
  | VTL synthesis | 0.40–0.45 | 1.2–3.4 | 3–11 |
  | m2i | 0.39 | 3.5 | 10.8 |
  | m30 | 0.37 | 3.1 | 9.5 |
  | m28 (sheet) | 0.40 | 4.1 | 11.7 |
  | **today, owner's sentence** | **0.82** | **4.1** | **13.0** |

  - The engine since m30 doubled the cycle-to-cycle jitter: "sickly". RAP does not separate m28's "vibrator" from m2i, so it is not the warble's measure.
- **Cause: the adductor tremor.** Single ablations (RAP median):
  - tremor off (all three channels): 0.32%;
  - breath tremor off: 0.76%;
  - tension, posture, breath or wall variation off: 0.83–0.84%;
  - wander off: 0.80%;
  - tension tremor alone: 0.43%.
  - The 2–15 Hz band was the difference of two one-pole low-passes. Above 15 Hz it kept a 6 dB/oct skirt to 500 Hz, so the adductors moved the folds' rest gap from one cycle to the next.
  - Now two poles at each corner, with the scale measured once at construction. RAP **0.52%**, p90 2.8%, jumps 7.7 (from 0.82 / 4.1 / 13.0).
- **The pitch plan warbled.** The sloped targets (`accent_peak` / `accent_fall` / `unstressed_drift`) with per-vowel learned pitches swung the intended pitch on almost every syllable. "a little more human today": 117 → 138 → 119 → 114 → 111 → 106 → 89 → 120 → 140 → 106 Hz, 7.5 semitones up within 150 ms. The owner's reading of the three sentences moves in one rise and fall per phrase.
  - `PlanOptions::pitch_glide` (default on; `--pitch-glide 0` restores the sloped targets) with `unstressed_dip` 0.5 semitone (`--unstressed-dip`). Stressed vowels and each phrase's first and last vowel keep their planned pitch. The vowels between lie on the line (semitones, by time) between their kept neighbours, 0.5 st lower. The pitch heads continuously from one anchor to the next, so it turns only at stressed vowels.
  - Same sentence: 119 → 131 → 141 → 134 → 113 → 110 → 136 → 128 → 104 Hz.
  - Pitch-contour reversals (`analysis/contour.py`, turns ≥ 1 st per voiced second): owner 0.84 / 1.11, bdl 0–2.3, VTL synthesis 0.96–1.22, m2i 0.67, m30 0.50, m28 2.09. Noisy per sentence; not a validated steer.
- **Nasals ("evening", "human").** Wideband "evening" against bdl a0003: our /n/ and /ng/ keep vowel-like energy to 4 kHz, with the tube's 1 kHz-spaced upper resonances (5, 6, 7 kHz; 17.5 cm) standing as stripes.
  - A synchronous fold on 1 ms showed only +4.6 dB of tick-locked content: not a control-rate comb.
  - Nasal murmur bands (`analysis/nasalband.py`; middle 50% of n / m / ng re the loudest vowel; 0–0.5 / 0.5–1 / 1–2 / 2–4 / 4–8 kHz) over 111 nasals in a0001–a0030:

    | Setting | 0–0.5 | 0.5–1 | 1–2 | 2–4 | 4–8 | Distance |
    |---|---|---|---|---|---|---|
    | bdl (148, a0001–a0040) | −10.3 | −29.6 | −38.0 | −36.9 | −50.7 | — |
    | ours | −8.6 | −34.2 | −45.1 | −32.8 | −48.6 | 19.6 |
    | no sinuses | −9.2 | −21.7 | −32.5 | −30.5 | −48.2 | — |
    | nose loss 0.35 | −8.7 | −37.5 | −46.1 | −36.5 | −50.6 | 18.1 |
    | no sinuses, nose 0.4 | −9.4 | −22.7 | −39.4 | −33.9 | −51.4 | 12.9 |
    | **sinus necks ×0.25, nose 0.4** | −10.1 | −30.7 | −40.1 | −35.5 | −51.5 | **5.6** |

  - The sinuses' side resonances were far too deep. Adopted `BodyLayout::sinus_neck_scale` 0.25 (ostia 0.03–0.05 cm²) and nose hop loss 0.4.
  - The skin (`SKIN=0`) and stiffer walls (`WALLK=10`, which raised everything above 1 kHz by 15 dB) are not the cause.
- The velum and under-tongue ports now glide in two 1 ms stages (they were one 2 ms stage) so the port's speed has no step at the control tick. No measured change; kept for consistency with the cord's glide.
- Identification (seed 1) with the tremor fix: 39.4% with or without `pitch_glide`; seeds 2 and 3: 39.1 and 38.8% (40.0% mean before, 40.8 / 38.9 / 40.2).
- **Incomplete closure in /iy/.** `CYCLES` now times each cycle from the glottis opening (area rising through 0.02 cm² after falling below 0.01) and logs the cycle's smallest area. The lower edge's contact alone missed cycles whose upper edge shut, and every cycle of a pressed posture.
  - `analysis/closure.py`: the share of voiced cycles (audio waveform matching) in which the glottis shut (smallest area < 0.001 cm²; the chink is 0.0002).
  - a0001–a0015, posture (lower, upper):

    | Posture | /iy/ shut | Others shut | RAP median | RAP p90 | Jumps / 100 |
    |---|---|---|---|---|---|
    | 0.004, −0.012 | 60% | 86–99% | 0.52% | 4.0% | 12.0 |
    | 0.002, −0.012 | 30% | — | — | — | — |
    | 0.0, −0.012 | 15% | — | — | — | — |
    | −0.002, −0.014 | 7% | — | — | — | — |
    | 0.006, −0.012 | 79% | — | 0.42% | 2.5% | 6.9 |
    | **0.008, −0.012** | **83%** | 92–99% | **0.32%** | **1.5%** | **3.8** |

  - Pressing the edges together closed the glottis less, not more.
  - At 0.008 the miss from the intended pitch is 39 cents (64 at 0.004). H1−H2 / vowel tilt on the four m39 scores: +1.3 / −19.2 at 0.004, +1.7 / −20.1 at 0.008.
  - The owner's NATO reading: RAP 0.23%, p90 1.5%, jumps 3.9. Adopted 0.008.
  - A wider /iy/ constriction also helps (varea 0.2: 38%, 0.3: 60% by the old lower-edge count) but moves its formants; left at 0.138.
- Identification with all of today's changes (glide, tremor band, nasal, posture 0.008): **38.0%** (38.0 / 37.9 / 38.1), against 40.0% yesterday. Steps on seed 1: 40.8 → 39.4 (tremor band; glide on or off alike) → 39.1 (nasal) → 38.0 (posture). Each step is within about a seed's spread; together they cost 2 points. Recorded; today's steer is the owner's warble.
  - Per phone (seed 3): g 0, ow 0, v 4, uh 9, ng 10, th 11, l 15, y 15 … n 68, f 71, s 73%.
- `listen/m43_steady_voice.wav`: "I hope I sound a little more human today." / "Good evening. I hope the voice sounds a little more human today." / "Hello, my name is Papa, and I speak with a body. The quick brown fox jumps over the lazy dog." (`talk.sh`, seed 73021).

### Owner, on m43 (2026-10-14)

"hello my name is papa and i speak with a bally(ringing) the quip brown fox jumpsh over the lazy dog". No mention of warble.

- **"bally (ringing)".** Wideband spectrogram of "with a body": through the /b/ closure (and the /dh/), steady stripes at every 1 kHz to 7 kHz. With the lips shut, the 17.5 cm tract is sealed at both ends and rings at n·1 kHz. It reached the air through the skin path (the lining's flow, heard whole); `SKIN=0` removes the stripes but leaves the closure silent.
  - Now `BodyLayout::skin_corner` 250 Hz, two poles (`SKINCORNER`): the neck and cheeks pass only the lows. The /b/ closure keeps a voice bar under 300 Hz and the stripes are gone.
  - The flap of "body" does make contact (0.001 cm² for about 6 ms, under 0.02 for about 20 ms).
- **"jumpsh".** The /p/ released its lips to the next vowel's posture: "over"'s /ow/, 0.29 cm² and protruded, so the /s/ was spoken through rounded lips. Its spectrum: peak 5.0 kHz, 1.5–3 kHz at −16 to −20 dB (the /s/ of "fox": 6.5 kHz, −25 to −30 dB).
  - `PlanOptions::sibilant_lips` 1.5 cm² and `sibilant_protrusion` 0.1:
    - a labial released into /s z/ opens to at least these;
    - an /s z/ not after a labial holds them through its closure.
  - Now "jumps": 1.5–3 kHz at −23 to −27 dB, peak 5 kHz. bdl's /s/ (a0001): peak 4.5–5.5 kHz, 1.5–3 kHz at −22 to −29 dB.
- **"quip".** The /k/ of "quick" before "brown" is unreleased (the "uplep prawn" fix), so its place shows only in the vowel's transition. `velar_place` took its neighbour from the next vowel (/aU/, back), so the tongue slid back 0.85 cm from /I/ into the closure.
  - Formant transitions into closures (`analysis/pinch.py`; LPC F2 and F3, the vowel's last 25 ms against its middle; bdl 60 sentences against ours):

    | | Velar dF2 | Velar d(F3−F2) | Labial dF2 | Labial dF3 |
    |---|---|---|---|---|
    | bdl | +65 | −46 (the pinch) | −24 | −22 |
    | ours | +28 | +25 | +66 | +62 |

  - A velar right after its vowel and not right before one now takes that vowel as its neighbour.
- With all three fixes, seed 1: identification **39.3%** (38.0 before).
  - Per phone: b 23 → 50%, s 73 → 78%, p 43 → 37%, k 22 → 15%, m 27 → 21% (against seed 3's row).
  - Transitions: velar d(F3−F2) +25 → +3 Hz (bdl −46); labial dF2 +81, dF3 +205 (bdl −24, −22).
  - The measure may be confounded where our acoustic closure does not fall at the score's vowel end. OPEN: why our F2 and F3 rise into labial closures.
- `listen/m44_body_quick_jumps.wav`: the "Hello…/fox" sentences first, then the two "more human today" sentences.

### Owner (2026-10-15): "he still sounds sickly. body sounds like baldy. human sounds odd."

- **Sickly: measured, not found.**
  - Within-vowel pitch and level wobble (`analysis/vowel_wobble.py`, the ARCTIC set): bdl 4.0 cents / 0.26 dB, ours 5.2 / 0.25.
  - Harmonics-to-noise by band (`analysis/bandhnr.py`, cycle-to-cycle coherence on the waveform-matched cycles). m44 is 20.9 / 9.1 / 4.3 / −7.0 dB at 0.08–1 / 1–2 / 2–4 / 4–7 kHz. Inside the range of the owner (16.9–18.8 / 3.9–8.6 / 3.4–6.4 / −5.3–0.2), bdl and VocalTractLab (originals and synthesis).
  - H1−H2 stays near +2 dB (owner −3.5, bdl −2.7) under every source change tried:

    | Change | H1−H2 | Vowel tilt |
    |---|---|---|
    | upper edge −0.016 / −0.020 | +2.2 | — |
    | lower 0.010, upper −0.020 | +2.7 | — |
    | damping 0.15, 0.8 | +1.1 | −16 |
    | upper damping 0.6 | +2.5 | — |
    | no skin | +2.1 | — |
    | skin corner 1 kHz | +1.6 | — |
    | 950 Pa | +1.9 | −14 |
    | 1150 Pa | +2.1 | −13 |

  - Each cycle is open about 60% of the period (open 5.2 ms of 8.7), peak area 0.105 cm², and the flow is skewed (peak late, sharp cut-off). A normal male open quotient. OPEN.
  - A/B for the owner's ear (`listen/sickly_ab/`, `analysis/sickly_ab.sh`): the two sentences as heard (hand rules), with the learned models, firmer folds (damping 0.15), stronger breath (950 Pa), perfectly steady (all variation, tremor and wander off), and higher pitch (median 140 Hz).
- **`talk.sh` never loaded the learned pitch and duration models.** `say` reads `data/pitch_en.txt`, `data/durations_en.txt` and `data/phones_en.txt` relative to its working directory. Run from anywhere else it silently falls back to the hand rules. Every listening file since `talk.sh` (m42–m44, and the `/tmp/es/abl` renders) used them, while identification (`render_arctic.py`, run from `voice/`) used the learned models. `talk.sh` now runs `say` from `voice/` and passes extra `say` options. Contour reversals: hand rules 1.10, learned 1.49 per voiced second.
- **"baldy".** LPC tracks through "a body": F1 fell from 730 to 490 Hz over the last 15–20 ms of the /a/, under F2 1400 and F3 2400 (a dark /l/). Three causes in the planner's flap (DX):
  - the jaw rose 40 ms ahead, as for a stop;
  - the tongue body went to the alveolar locus and was then widened for a voice bar (area ×1.8), interrupting /a/ → /i/;
  - the next vowel's tongue began only 20 ms before the tap, so the body was still back when the tip rose, and F2 fell into the tap (1240 → 1090) as for a labial.
  - Now a flap keeps the jaw open, skips the locus and the voice-bar widening, and the next vowel's tongue starts 60 ms ahead. F2 rises into the tap (1320 → 1480), the closure reaches −25 dB (from −15), and /i/ reaches F2 2200 within 40 ms of the release.
- **"human".** The lips never sealed for /m/: aimed only to touch, they left 0.001–0.004 cm² open (the tip's /n/ sealed fully). `PlanOptions::lip_press` 0.01: /p b m/ aim at a slightly pressed seal (0.000 throughout). The /m/'s spectrum barely changed: its peak is still near 2 kHz with strong 1–2.5 kHz. The word's timing (/uw/ 113 ms, /m/ 111, /ax/ 40, /n/ 85) gives a long nasal against a clipped /ax/. OPEN, and now planned by the learned durations.
- `say --pitch-median` (the learned model's reference pitch).
- Identification with the flap, jaw and lip-press changes (seed 1): **38.6%** (39.3% before them). /b/ 57%, /m/ 27% (21), /p/ 43% (37); /d/ 30% (33); /dh/ 4% (6). Within a seed's spread; recorded.
- Owner on the A/B: "they sound exactly as sickly. all of them. they all sound just like a 1998 tts model". None of the five changes (learned models, firmer folds, stronger breath, all variation off, higher pitch) moved it.
- Control tick 1 ms → 0.125 ms (`kControlSteps` 12 in `units.hpp`; the owner: "timing resolution 8x", meaning the intermodulation behind the warble). Same speed (4.2× real time).
  - Inter-harmonic energy (`analysis/interharm.py`: cycles time-warped to 256 samples, 8 per FFT, harmonics on every 8th bin, ±2 bins counted harmonic), nine sentences:
    - noise and variation off: 1 ms −24.3 / −16.0 / −13.3 dB at harmonics 1–8 / 8–24 / 24–48; 0.125 ms −23.0 / −14.0 / −12.3;
    - full: 1 ms −23.4 / −11.7 / −8.0; 0.125 ms −23.5 / −14.2 / −9.0;
    - people and VTL: −17 to −28 / −9 to −19 / −1 to −11.
  - A held /a/ is 40–50 dB clean between its low harmonics. The tick was not the intermodulation; kept.
- Owner: "we need better articulation, and more consistent articulation ... our engine is possibly a tiny demo and not a real engine. we need to study their engine as well to figure out how many finite elements they use and think about their budget compared to ours". Study in `learned/vtl_budget_vs_ours.md`:
  - VTL 2.4 (Dec 2025): its articulation geometry is 12,605 lines and about 90% of its compute (ESTIMATE): 3D surfaces, a 407-vertex tongue, 129 cut planes × 96 lateral samples, rebuilt at 400 Hz, with section lengths from the anatomy. Ours is a 109-line parabola with two notches, a fixed 17.5 cm tract and about 1% of our compute. The budgets are inverted.
  - Proposal: a midsagittal articulatory model (one tongue, the jaw as a carrier, larynx height, lips with length), a variable-length tube (B1 waves with a fractional lip end, or B2 a lumped line), and gestures in physical parameter space with context consonant prototypes.
  - Redone on **VTL 2.4** (December 2025, the owner's choice; it was VTL3D 1.0.2 first). Its speakers select the **"Geometric glottis 2025"**, a kinematic valve:
    - cosine edges at exactly the commanded F0, the upper lagging 1.22 rad;
    - amplitude ∝ √(transglottal pressure);
    - exponential contact and a pulse-shape skew;
    - a 0.03 cm² chink;
    - Klatt flutter and phase-random-walk jitter.
  - It is not the self-oscillating triangular glottis we ported. Owner's decision: keep our cord or try a kinematic valve (A/B).

### Kinematic folds, an A/B (2026-10-15; owner: "yes and lets try b1")

- `Cord::set_kinematic` (`BodyLayout::kinematic_glottis`; `speak` `KINEMATIC=1`, `KAMP`, `KLAG`, `KSOFT`, `KJITTER`, and `CHINK` for the voicing chink). Our own design after Titze's four-parameter glottis, as VTL 2.4's default source; nothing copied:
  - each edge a cosine at the intended pitch × the holding factor, the upper lagging 1.22 rad;
  - amplitude 0.0007 cm per √(dyn/cm²) × √(transglottal pressure, two 25 Hz poles), fading out between 200 and 100 Pa;
  - the opening integrated over 16 points along the fold (rest gap linear to the vocal process, swing a half sine), with a softplus contact 0.002 cm wide;
  - jitter a phase random walk, 0.17 rad per √s;
  - the flow, trachea and tract unchanged; the ear's correction off.
- Owner's sentence: |deviation from the plan| 8 cents (self-oscillating 33); RAP 0.26% by audio (0.32; owner 0.23–0.42); 276 voiced cycles (242).
- H1−H2 and vowel tilt (the four m39 scores):

  | Source | H1−H2 | Tilt |
  |---|---|---|
  | self-oscillating | +2.2 | −20.3 |
  | kinematic | −5.0 | −14.5 |
  | kinematic, chink 0.03 (VTL's default) | −1.3 | −20.8 |
  | kinematic, chink 0.01 | −4.1 | −24.2 |
  | kinematic, chink 0.005 | −4.7 | −19.3 |
  | kinematic, softness 0.003 | −4.9 | −18.1 |
  | kinematic, rest 0.015 / 0.01 | **−4.0** | **−19.0** |
  | kinematic, rest 0.02 / 0.02 | −3.2 | −24.4 |
  | owner | −3.5 | −22 |
  | bdl | −2.7 | −18 |

  - H1−H2, stuck near +2 under every self-oscillating setting, reaches the human range at once.
- `listen/source_ab/` (`analysis/source_ab.sh`; "Hello … body", the fox, "I hope … today"; learned models, room):
  - A: self-oscillating;
  - B: kinematic, rest 0.015 / 0.01;
  - C: B with chink 0.03.
  PENDING the owner's ear.

### The articulatory body: A (midsagittal model), B1 (moving open end), C (organ planner) (2026-10-15)

Owner: "yes and lets try b1".

- **A. `tract/articulation.{hpp,cpp}` (`MidsagittalTract`).** Our own model after Mermelstein's and Maeda's semipolar grids and Heinz & Stevens's area rules. Unscaled coordinates in cm, x forward, y up, the glottis at (1, 0).
  - A fixed hull: larynx back wall, posterior pharyngeal wall, velum, palate, alveolar ridge, upper incisors; the upper lip added per shape.
  - The organs as one closed outline: the larynx front wall; the root, from the epiglottis base (following the hyoid and, 17%, the tongue body) tangent to a 2 cm body circle; the blade tangent to the tip; the tip's underside; the floor; the lower incisors and lower lip (carried by the jaw's rotation about the condyle); the chin.
  - Cuts: horizontal through the larynx and pharynx (every 0.3 cm), a fan from (3.4, 5.6) through the mouth (5.5°), vertical at the front (every 0.3 cm).
    - Width = the hull crossing minus the nearest organ crossing below it. Zero where the hull point lies inside the organs (a contact).
    - Area = α·d^β by region (larynx, pharynx, velar, palatal, alveolar, teeth); the lips an ellipse, width × aperture.
    - Lengths = distances between successive midpoints, so the tract lengthens with a lowered larynx and protruded lips.
  - Floors: side channels past a tip contact (lateral, 0.2 cm²) and a groove along the blade (0.12 cm² at the tip, widening to 0.5 two cm back), applied even where the midline is closed.
  - Fitted to Story's (2008) eleven MRI vowels together (`analysis/geofit.py`, `geofitplot.py`): rms log-area error 0.36. Lengths within 0.3 cm except /ae/ +0.9 and /iy/ −0.5 (Story 16.1–19.6 cm). Fitted speaker scale 0.953, root follow 0.17.
  - Bugs fixed on the way:
    - an organ pressed past the hull made a fan cut take the tongue's underside: an 8 cm "gap", the length jumping to 25 cm and the tube blowing up at the first /k/;
    - a zero-width cut skipped the groove floor, so /s/ became a stop.
  - `tools/geometry_probe` (key=value lines → cuts and outlines).
- **Vowel postures** (`tract/data/organs.txt`): the MRI postures tuned to bdl's F1–F3 by coordinate search (`analysis/organ_vopt.py`, held vowels through `speak`, kinematic source, penalised 0.5 semitone per cm of tongue movement). Before → after:
  - /aa/ 5.96 → 0.53 st (717 / 1127 / 2506 Hz);
  - /ah/ 6.65 → 0.43;
  - /uw/ 9.54 → 1.94 (340 / 1273);
  - /ey/ 2.81 → 0.64;
  - /ih/ 0.76; /eh/ 0.98; /ae/ 0.84; /ow/ 1.73; /uh/ 2.44; /iy/ 2.10 (F2 1975 against bdl's 2213); /ao/ 2.63 (stuck).
- **B1. `Tube::set_active_length`.** Only the first `active_count()` nodes carry air. The radiation load sits on the last, behind an extra air mass ρ·(remaining length)/A, so the acoustic length is continuous. Lengthening waits a quarter node past each boundary (hysteresis), and newly joined links start silent.
  - Side ports move (`Tube::move_side_port`; `VoiceBody::set_velum_place`, `set_under_tongue_place`).
  - The jets use the moving lip end.
  - `articulatory_layout()`: 60 nodes (21.9 cm).
- **Speaker.** A score holding organ parts (`body`, `apex`, `mouth`, `spread`, `larynx`, `angle`, `groove`, with `lateral` as 0..1) is spoken by the articulatory body (`uses_organs`).
  - The geometry runs every 8 control ticks (1 ms); node areas (`profile_to_nodes`: the centre's interpolated area, or a narrower cut inside the node) and the length glide linearly between frames.
  - 3.2–3.6× real time.
  - `PROFILE=t` also prints the held organs and cuts.
- **C. The organ planner** (`say --organs FILE`; `PlanOptions::organs`; `gesture_organs`):
  - Vowels set every organ, from their onset's start, arriving late in the vowel; diphthongs glide (AY = AA→IH, AW = AA→UH, OY = AO→IH, EY→IH, OW→UH); reduced vowels open the jaw and lips less.
  - A consonant constrains only its own articulator, placed against the hull; everything else follows the vowel it heads to (or leaves):
    - /p b m/ press the lips;
    - /f v/ a 0.06 cm lip slit;
    - /th dh/ the tip at the incisors' edge;
    - /t d n/ (and the flap) the tip 0.12 cm past the ridge;
    - /s z/ the tip 0.1 cm short of the ridge, with groove;
    - /sh zh ch jh/ further back, with groove and lips rounded;
    - /l/ a tip contact with side channels and the body drawn back;
    - /r/ bunched;
    - /k g ng/ the body raised against the roof above where the vowel holds it (palatal beside front vowels, velar beside back ones);
    - /w/ lips rounded, the body toward /uw/;
    - /y/ the body toward /iy/.
  - The area-profile parts are no longer written.
- **First tuning of the articulatory body.** Identification, one seed (own folds / moved folds):
  - first planner: 22.1 / 18.1%.
    - Almost every failing phone was heard as /r/: vowel F3 150–830 Hz under bdl (`analysis/f123.py`, running-speech medians). The vowel tuning had weighted F3 a quarter, against rough targets.
    - Retuned to bdl's measured F1–F3 with full weight at speaker scale 0.90 (Story's 0.953 is a larger speaker: mean distance over six vowels 3.49 / 3.09 / 2.87 st at 0.953 / 0.92 / 0.89). Identification 21.9 / 18.6%, nearly unchanged.
  - **Leaking closures** (`analysis/phoneband.py`, closure middles re the loudest vowel): /b d m z v/ carried −23 to −31 dB at 0.5–4 kHz where bdl's are −42 to −70.
    - A pressed lip aperture put the lower lip *below* the upper (a 0.1 cm slit, about 0.25 cm²): fixed.
    - The planner's vowel events re-set every organ just before the onset consonant's closure, and a consonant's hand-back reopened the next consonant's closure ("is Papa", "and"). Now each consonant owns its organs (`owned_organs`): a vowel leaves them to its onset consonants, a consonant hands back only what the consonants after it do not own, and the flap aims for the jaw midway between its vowels.
    - No cut passed through a tip contact: the fan ended at 52° (5.5° steps) and the vertical cuts began at 6.7 cm. Now the fan runs to 38° in 2.5° steps over the ridge, and the vertical cuts begin at the teeth. /t d n/ and the flap close; /s z/ hold their 0.12 cm² groove.
    - The area rules refitted to the MRI on the new grid (rms 0.39, root follow 0.22) and the vowels retuned: most within 0.2–1.5 st of bdl; /uw/ 2.7, /uh/ 2.6, /ao/ 3.3, /iy/ 2.0.
  - Now **25.2 / 24.7%** (old body 38–39%). Per phone (moved folds): p 48, t 44, b 33, d 22, m 23; still near 0: th v dh jh z y ng uw uh ow g w. /w/ is heard as /b/ (its lips too closed), /z/ as /n/, /v/ and /uw/ as /r/.
- `listen/m45_new_body_{moved,own}_folds.wav`: the three test sentences through the articulatory body (organ planner, learned models, room).
- Owner on m45: "he warbles less, but his voice is far weaker and there's some clicking as if going into overload".
  - `analysis/clicks.py` (largest samples re the 99th percentile, with the phone): the new body's peaks stood **20.8 dB** over the speech (the old body's 7.1). Every click fell exactly when the tube **joined a node** to its open end (`speak` `ENDLOG=1`: 0.7531, 1.8865, 2.1411, 3.4784, 4.6475, 6.0230 s …); dropping nodes was quiet.
  - The joined link started silent, so the old end's pressure collapsed for a step. The output normalisation then scaled the speech about 14 dB down: the "far weaker" voice.
  - Now a joined link starts with the waves that carry the old end's pressure and radiated flow, the new nodes at that pressure. The tube does not lengthen while its end is shut (lips sealed: the huge impedance made a click in "jumps"'s /p/); the length waits on the end's air mass (up to 2.5 nodes).
  - Peak re 99th: moved folds 20.8 → 8.0 dB; own folds 11.2 (a real /jh/ release against 600 Pa, probably too loud). Listening renders (room): 7.5 / 8.4 dB (old body 7.9).
- `listen/m46_new_body_{moved,own}_folds.wav`: m45 without the clicks.

### Owner on m46 (2026-10-15): "he sounds like he has a lisp or a fork in his mouth, but he isnt warbling. hello hiccups. papa clips. i speak with a body has a lolipop in his mouth. the quick brown fox ... is barely articulable. but this is progress. lets move forward."

- `analysis/phonebands.py`: every phone's median spectrum in 7 bands (0–0.5 … 6–8 kHz) re the loudest vowel, ours minus bdl over a0001–a0060, sorted by mean |difference|. Tuning on the moved folds.
- **Lisp.** /s/ was 25 dB short at 4–8 kHz, /z/ and /jh/ 34–41. The constriction sat at the ridge's crest, 1.3 cm behind the incisors, past the jet's coherence (0.8 cm).
  - /s z/ now 1 cm behind the incisors' edge, jaw 1.5°; /s/ 4–8 kHz short by 9–11 dB.
  - The groove is narrowest at the tip (0.06 cm², widening to 0.4 two cm back); 0.12 let /z/'s air through too slowly to hiss.
  - `voiced_fricative_gap` 0.012 → **0.04**: at 0.012 the folds shut every cycle, the mean flow was about 50 cm³/s and the pressure behind the /z/ channel about 28 Pa. /z/ 4–7.5 kHz re its peak: −31 / −22 / −17 / −14 dB at 0.012 / 0.025 / 0.04 / 0.06.
  - /sh/ moved forward to 6.65 cm fell from 33 to 0% identified; kept at 6.3.
- /th dh/ against the incisors' back face (aimed at their edge, the tip stayed 0.15 cm under: 0.35 cm²); /dh/ pressed lightly (0.04). /f v/ slit 0.03 cm. /w/ aperture 0.45.
- **"hello hiccups".** /l/ drew the body 0.5 cm back from /ow/'s already-back posture; the body reached its limit and closed the pharynx from 6 to 7.7 cm, and the voice stopped for 100 ms.
  - /l/ now retracts 0.2 only from a forward body (x > 3.4).
  - The geometry keeps the pharynx region at least 0.15 cm² open (English never seals it).
  - Side channels 0.35 cm².
- **"lollipop".** Vowels averaged +4.5 dB at 1–2 kHz and −6 to −7 at 3–8 kHz against bdl (the old body: +2.4 at 1–2, −1.7 at 3–4, −11 at 4–6).
  - Not the teeth passage (doubled: no change). Not the epilarynx (2.9 / 3.4 cm: worse). Not bandwidths: B1–B4 ours 149 / 106 / 116 / 126 Hz, bdl 99 / 112 / 159 / 142.
  - F3 is low in the back vowels (running speech: /ah/ −199, /ow/ −204, /uh/ −402, /uw/ −302 Hz) and /er/ had no r-colouring (+396).
  - Powell's method (`analysis/organ_vopt2.py`) where coordinate search stuck:
    - /er/ 417 / 1433 / 1923 Hz (bdl 419 / 1438 / 1853), with the pocket under the tongue open;
    - /uh/ F3 2615;
    - /uw/ 335 / 1237 / 2468 (an odd posture: open jaw, low body);
    - /ah/ 2633.
- Identification (moved folds): 24.7 → 27.9 (sibilants, voiced fricatives) → 27.1 → **26.5%**.
  - s 69, p 59, aa 56, ay 75, t 40, b 33;
  - still 0: th v dh jh y uw uh ow (the rounded vowels too, though their held formants now match).
- `listen/m47_new_body.wav`: the three sentences with these fixes (moved folds).
- Owner on m47: "he's stil clicking like audio overload."
  - `analysis/clicks2.py`: local clicks, samples of the 5 kHz high-pass above K× its rms over the surrounding 20 ms (K 8), per second. `clickphones.py`: per phone.
    - Owner 5.6/s; bdl a0001 5.6/s (K 5); VTL synthesis 10.3/s.
    - m47 12.3/s (K 8), many of them exactly at the open end's node changes, removals as well as joins.
  - **The open end no longer moves node by node.** All 60 nodes always carry air. Past the lips the area flares 1.6× per node from the lip opening to 30 cm² (the room), and the node the opening crosses blends by how much of it lies inside. The lips move without any switching; `VoiceBody::set_lip_end` tells the jets where the lips are.
  - The remaining clicks recur every glottal cycle, densest through /dh/ (84/s), /ao/, /n/, /d/, /iy/: closures and narrow tracts, where nothing else is high.
    - A per-node pressure dump (`speak WAVEDUMP=from,to`) traced one from node 0 (the glottis, at closure) up the tube one node per step: the closure's own excitation.
    - Two changes on the way, neither the cause of that excitation:
      - the lining's contact with the opposite wall is now a critically damped 1.5 kHz spring, not a hard stop (no change in the count);
      - the glottal air mass is taken on areas of at least 0.01 cm² (rho l / A blew up at closure in the trapezoidal step), which brings the own folds to 5.6/s.
  - m48: moved folds **10.1/s** (VTL synthesis 10.3), own folds **4.6/s** (owner 5.6). Peak re 99th 8.9 / 10.4 dB.
- `listen/m48_new_body_{moved,own}_folds.wav`.
- Owner on m48: "sounds like he's tapping the microphone while speaking. both sound the same."
  - `analysis/thumps.py` (bumps of the 0–120 Hz envelope 6 dB over its local median, and the band's share of energy):
    - m45–m47: none, −19 dB; m48: 3.5/s, **−8.8 dB**;
    - old body 0.2/s, −18.3; owner 0.0/s, −17.4; VTL synthesis −9.4.
  - Bisected over the saved renders: it came with the flared open end. Not its walls (made rigid: no change), not its shape (FLARE=growth,widest: no change), not the noise sources or branches.
  - A step dump (`speak STEPS`) showed the frame pumping −3300 cm³/s within 2 ms in /z/: a contact regrouped the cuts, the length jumped by a node or more within a frame, and the lips' nodes swapped with the room's, whose large area changes counted as moving walls.
  - Now:
    - the room's nodes push no air (`Tube::set_quiet_from`, masked before the frame's two 250 Hz poles; masked after them, each switch was a two-sample spike);
    - the length follows the geometry through an 8 ms one-pole (and starts at the opening posture's length: a default 17.5 cm popped at the first sample);
    - the room's areas no longer follow the lips (2 cm² growing 1.6× per node to 30).
  - /k g/ squeezed 12 cm³ out of the sealed space in 6 ms (body 0.15 cm into the dome, 40 ms in, 30 out): the pressure behind passed the lungs' and the burst stood 17 dB over the speech. Now 0.05 cm, 60 ms in, 50 ms out.
  - Result (owner's sentence set, moved folds): thumps 0/s, 0–120 Hz share −21 dB; peak re 99th 8.8 dB, in vowels; local clicks 12.4/s.
- `listen/m49_new_body_{moved,own}_folds.wav`.
- Owner on m49: "sounds a bit muffled. own folds still clicks a little. both voices seem to sound like a classic TTS model."
  - **Muffled.** Vowel bands (mean over 10 vowels, ours − bdl, 0–0.5 … 6–8 kHz), moved folds:

    | Setting | 0–0.5 | 0.5–1 | 1–2 | 2–3 | 3–4 | 4–6 | 6–8 |
    |---|---|---|---|---|---|---|---|
    | rest 0.015 / 0.010, hop loss 0.05 | −3.9 | +1.3 | +7.6 | −2.1 | −6.8 | −8.1 | −3.8 |
    | loss 0.02 | −4.2 | +1.6 | +7.5 | −2.0 | −6.4 | −5.8 | −2.0 |
    | loss 0 | −4.0 | +1.7 | +7.6 | −1.6 | −6.3 | −4.4 | −0.7 |
    | walls ×3 / ×10 support, or loss 0.4 | (no change) | | | | | | |
    | own folds | −2.8 | +0.9 | +5.7 | −1.9 | −3.7 | +0.6 | +11.9 |
    | rest 0.008 / 0.004 | −4.2 | +1.5 | +8.0 | +1.0 | −1.0 | 0.0 | +5.7 |
    | **rest 0.011 / 0.006, loss 0.02** | −4.2 | +1.3 | +7.9 | +0.2 | **−2.8** | **−2.2** | **+2.0** |

  - Adopted the last as the articulatory body's defaults (`articulatory_layout()` and `Speaker`): the kinematic folds, rests 0.011 / 0.006, hop loss 0.02. `speak KINEMATIC=0` restores the self-oscillating folds.
  - The 1–2 kHz excess (+6 to +8, with 0–0.5 kHz −3 to −4) is in every source variant: the body's vowels, not the source. B1 146–149 Hz against bdl's 99. OPEN.
  - **Classic TTS.** `listen/prosody_ab_new_body.wav`: for a0001, a0003, a0008 in turn — bdl himself, ours as planned, ours with bdl's durations, vowel pitches and pauses (`tools/prosody_transplant.py`, `ENGINE=tract`, organs). PENDING the owner's ear: if ours with bdl's prosody still sounds like TTS, the fault is the body; if not, the planner.
- `listen/m50_new_body.wav` (new defaults).
- Owner on m50 and the prosody A/B: "this thing STILL sounds like perfect paul. in both cases" (bdl's prosody did not help: the fault is the body), then "it sounds like im listening to a merman, gurgling his words".
  - **Perfect Paul.** Band HNR (`analysis/bandhnr.py`, dB, 0.08–1 / 1–2 / 2–4 / 4–7 kHz):
    - m50: 20 / 15 / 8.0 / +3.5;
    - bdl: 5.5–9 at 1–2, 3–4.5 at 2–4, −1 to −4.6 at 4–7;
    - owner: 3.9 / 3.4 / −5.3 at 1–2 / 2–4 / 4–7.
    - Ours is too periodic above 1 kHz.
  - **Merman.** LTAS (`analysis/ltas.py`) of the same sentences: ours +7–10 dB at 1–2 kHz and −4 to −7 dB below 300 Hz against bdl. Held /aa/ F1 654, F2 1097; FLARE=1.0,1.5 dropped F1 to 380. The flared exterior acted as part of the tract.
  - **The lips radiate from the opening** (`Tube::set_opening`): a Flanagan R‖L load (R = 128ρc/9π²A, L = 8ρ/3π√(πA), A the last cut's area) shared between the two nodes around the opening by its fractional position. Each slot keeps the full load's inertance flow and the node draws its share, so the load slides with no switching. Past the opening the tube is a 0.02 cm² sliver; the right end no longer radiates.
    - Held /aa/: F1 704–716, F2 1211. LPC B1 reads 170–186, but with the kinematic pulse (open quotient 0.32) that estimate follows the harmonic spacing (KINEMATIC=0, open quotient 0.80, reads 94). Not a valid damping measure here.
    - k1 LTAS 0–300 / 300–600 Hz: m50 −12.5 / −7.5, now −9.4 / −6.4 (bdl −4 to −7 / −3 to −4.6).
    - No thumps (0–120 Hz share −19 dB). Local clicks 4.5/s, mostly at /dh/ contact.
    - ARCTIC F1–F3 per vowel now within about 100 Hz of bdl's (F3 too), ao/uh/ow F1 +50–90.
      - (Beware stale renders: the voice/out/arctic set from before this had every front vowel's F2 near 1100 Hz.)
  - **Breath noise.** `jet_flux` (FLUX) adds the noise the folds lack, but shaped like the jet noise it is nearly white at the glottis.
    - FLUX 0.04: vowels' 6–8 kHz +16 dB over bdl.
    - New `flux_corner` (FLUXCORNER, Hz): two low-pass stages of white noise at unity gain below the corner. (Keeping the rms instead piled the noise into the passband: 4–6 kHz rose.)
    - Vowel means (ours − bdl, the seven bands):
      - FLUX 0.08, corner 1500, KJITTER 0.5: −1.5 +2.4 +4.7 +0.7 −5.0 −0.2 +6.3;
      - plain: −1.6 +2.5 +4.2 +0.2 −5.7 −5.4 −2.4.
    - k1 band HNR with it: 18.9 / 10.9 / 5.6 / −1.3 (m51b); plain 20.9 / 13.4 / 8.6 / 4.0 (m51a).
  - **Fixed spectral structure** in every vowel (60 ARCTIC sentences, fixed scores, pooled main vowels, 250 Hz bands, `/tmp/es/ub/ab.sh`):
    - a notch at 3.0–3.5 kHz (−6 to −8) that moves with the epilarynx length: 1.8 cm deepens it to −12 at 3.0–3.25; 3.0 cm moves it to 3.5–3.75 and sinks 4.5–5.5 kHz by 6–12;
    - the piriform zero at 5.25–5.75 kHz (−7 to −8; without the piriform branch, 4.25–5 kHz stands +10–12);
    - the broad +5 at 0.75–2.25 kHz, with 0–250 Hz −3, in every source variant (KINEMATIC=0 included).
    - The under-tongue branch is not involved (371 vowels, median 0.0 dB at 3.25–3.5). An earlier render_arctic.py variant table that said it was is wrong; it is not reproducible.
    - OPEN: the 1–2 kHz excess and the missing 3–3.5 kHz (F4 region).
- `listen/m51_lip_radiation.wav` (lip-opening radiation, defaults) and `listen/m51_lip_radiation_breath.wav` (plus FLUX 0.08, FLUXCORNER 1500, KJITTER 0.5), owner's sentence.
- Owner on m51: "the second one sounds a little less, but perfect paul/the merman is still audible in both."
  - **Reference check.** Full-band LTAS (48 kHz) against the owner's own recordings (user_three, user_nato, 44.1 kHz):
    - m51b 4–6 kHz −29 against −19.7, 6–8 kHz −33 against −26 to −29 (bdl's 16 kHz recordings are dull there and hid it);
    - 1–2 kHz −6.7 against −9 to −14.
    - FLUX 0.04 with no corner gives 4–6 kHz −23, 6–8 kHz −25: the corner removed the breath the owner needs. FLUXCORNER stays available, default off.
  - **Decimation** (speak): now an eighth-order Butterworth at 16 kHz (was fourth-order at 18). The 8–20 kHz content did not change: it is the fricatives' own (a held /aa/ falls smoothly to −64 dB at 16–18 kHz), not aliasing.
  - **Measured equal to bdl, so not the cause:**
    - H1−H2: −3.6 against −3.3 (`/tmp/h1h2.py`);
    - shimmer: 6.6% against 5.5% (`analysis/shim.py`);
    - formant flutter: median |second difference| F1–F3 5 / 8 / 12 Hz against 4 / 7 / 10 (`analysis/fflut.py`);
    - vowel level spread: 4.8 against 4.2 dB (`analysis/vlev.py`).
  - **Formant amplitudes** (`analysis/famp.py`; mean ours − bdl: H1−A1, A1−A2, A2−A3, A3−A4 | B1–B4):
    - old (FLUX 0.08, corner 1500, rest 0.015 / 0.010): −2.3 −1.3 +4.3 +3.4 | +34 −25 −16 +25. Front vowels' F2 too strong and sharp (B2 60–80 against 85–145); F3 and F4 weak; F4 200–270 Hz high.
    - **FLUX 0.04, no corner, KJITTER 0.5, KSOFT 0.0005, rests 0.011 / 0.006:** −3.6 −3.9 +1.3 −0.7 | +22 −20 −13 −35. The 3.0–3.5 kHz hole fills (pooled bands −1.7 / +1.6, was −8 / −6). The hole was missing noise, not geometry.
    - Plus KLAG 0.5: +0.5 −4.6 +1.4 −0.6 | +58 −14 −10 −18.
      - The edges' 1.22 rad lag is what keeps the open quotient at 0.32: the opening is the narrower edge's, and two cosines 70° apart are both open about 30% of the cycle. KLAG 0.8 gives 0.38 and 0.5 gives 0.43 (mean flow 43 → 73 cm³/s). KAMP and VOICEPOST hardly move it.
    - KAMP up only adds mid and high energy.
  - **Piriform.**
    - Lossless, its zero cuts 17–23 dB at 4.6–5.0 kHz (noise-excited transfer, folds still). Without it, 4–6 kHz matches the owner, but over ARCTIC 4.25–5 kHz stands +10–12 over bdl.
    - New PIRIPORT (default 0.6 cm²): at 0.3 the zero falls to about 3.6 kHz and digs −8 to −14 into the F4 region. REJECTED.
    - New PIRISMOOTH (hop loss in the recesses, default 0): it shifted the pole-zero pair and dug −7 to −9 at 3.0–3.5 kHz. REJECTED.
  - **Epilarynx.**
    - 2.6–3.0 cm puts F4 within 50–75 Hz of bdl's, but deepens 4.5–5.5 kHz by 2–12 dB.
    - Not adopted; with the full breath noise the 3 kHz region no longer needs it.
  - **Still OPEN:** the 0.75–2.25 kHz excess (+5 to +8 against bdl, +3 to +8 against the owner). It persists under every source variant tried (lag, amplitude, rest, contact, KINEMATIC=0) and with the walls ×3 or lossier.
  - **Identification** (phone_test, FLUX 0.04, KJITTER 0.5, KSOFT 0.0005): 21.2%, human ceiling 52.1%.
    - Liquids and glides 10%, nasals 19%.
    - Frequent n→r, ah→r and l→r: an r-colouring, consistent with the mid excess. Possibly the "merman".
- `listen/m52_breath_full.wav` (FLUX 0.04 no corner, KJITTER 0.5, KSOFT 0.0005) and `listen/m52_breath_full_longer_pulse.wav` (plus KLAG 0.5).
- Owner on m52: "full sounds a bit more like talking through a tube. full longer pulse - a bit more like through a tube AND like he has a bit of a lisp. but its the same merman, except you cant hear his roboty-ness behind the lisp and tuby".
  - Breath noise excites every resonance the tract has, so a fixed "tube" colouring became more audible.
  - Ruled out:
    - the trachea: its lung end is matched, so it does not ring;
    - the sliver past the lips (new SLIVER switch, 0.02 / 0.002 / 0.0002 cm²): no change.
  - **The tube.** Our tuned vowel shapes are compressed toward a uniform tube against Story's MRI area functions (10 positions, glottis to lips):
    - /aa/: pharynx 0.58–0.82 cm² against 0.25–0.29; mouth 5.6 against 9.3;
    - /ah/: mouth 2.3 against 6.6;
    - /iy/: pharynx 1–3 against 3.3–5.2.
    - The model can reach Story's shapes (the geofit postures: /aa/ mouth 7.4, pharynx 0.15). The formant tuning flattened them.
  - **Why** (`analysis/tlcheck.py`: an ideal transmission line, chain matrices with piston radiation, through the model's own profile, against the engine speaking the same posture):
    - the engine's formants came out lower: /aa/ F1 632 against 770, /ah/ 552 against 650, front vowels' F2 6–10% low;
    - without the piriform branch the engine matches the line (/aa/ 773 / 1038 against 770 / 1020);
    - the ideal line with the same recess agrees (`analysis/tlpiri.py`: /aa/ F1 590), so the engine is right and the recess is too big.
    - At full size (1.8 cm³ beside a back cavity of a few cm³) it pulled F1 down 18%. Story's MRI shapes, which have no recesses, already give real formants, so that cannot be a real speaker's. The tuner flattened every shape to win F1 back.
  - **Now.**
    - Piriform cross-sections ×0.25 in the articulatory layout (new `piriform_scale`, PIRISCALE): the depth and 4.5 kHz zero stay; /aa/ F1 −9%.
    - Vowels retuned by `analysis/organ_vopt3.py`: cost = formant distance to bdl (semitones) + 10 × rms log-area error against Story at 16 positions, Powell from the geofit postures.
  - A first run's monitor waited on an output file while the script had crashed at import (it ran organ_vopt.py's main); five hours lost. The tuner now carries its own formant routine and the wait watches the process.
  - Mini disk reached 320 MB free. Freed my scratch, then, on the owner's direction, deleted `~/ducklips/out` (13 GB of the October 9–11 DuckLips outputs, including snapshot copies of experiment code). Now 14 GB free.
  - **Result of the retune.**
    - Formant error 0.2–0.9 semitones for iy ih ey eh ae aa uw; ah 1.35, ow 1.6, uh 0.26 (but its shape is rms 0.82 from Story's), ao 4.1 (F2 784 against 918).
    - Odd features: AA's lips protrude 0.69 cm; AH's jaw is 13°. To watch.
    - ARCTIC with FLUX 0.04, KJITTER 0.5, KSOFT 0.0005:
      - F1–F3 within about 100 Hz except ao and ow F1 (+170, +120) and er F3 (+228);
      - famp −4.1 −4.6 +1.8 +0.5 | +24 −14 −30 −4;
      - pooled bands still +4 to +7 at 0.75–2.25 kHz;
      - 4–5 kHz now +13 to +17 over bdl (smaller recesses, shallower zero), but owner's sentence 4–6 kHz −18.7 against his −19.7;
      - identification 21.5% (unchanged); ah→eh rose from 45 to 115 (the AH posture also serves unstressed vowels).
- `listen/m53_story_shapes.wav` (vowel shapes near Story's, piriform ×0.25, m52a's breath settings).
- Owner on m53: "it has a faint ringing. and he has a bad lisp."
  - Then: "a lisp doesnt mean that in this context. it just means he sounds like he has a serious speech impediment".
  - **Tongue tip.** The free tuner had lifted the tip toward the ridge in vowels (AH 9.0, IY 8.4, IH 8.0) and protruded AA's lips 0.69 cm.
    - organ_vopt3 now bounds the tip behind the lower incisors (x 6.7–7.4, y 6.6–7.4, jaw frame) and unrounded vowels' protrusion to 0.15 cm.
    - Retuned from the m53 postures: most vowels within about 1 semitone; ao 4.4.
  - **Schwa.** AX and IX used AH's posture (jaw 10–13°), so every reduced syllable was a full open vowel. bdl's schwa ("ax", his most frequent phone, 158 in a0001–a0060) is at 549 / 1449 / 2660.
    - New AX posture tuned on formants only, jaw 3–8°, lips 0.6–1.3 cm: 556 / 1434 / 2648 alone, 567 / 1487 / 2624 in ARCTIC.
    - The planner gives IX the AX posture.
  - **Breath noise gated to voicing.** It followed the whole glottal flow, so every voiceless consonant's continuous open-glottis flow carried a loud low hiss:
    - FLUX 0.04 against 0, mean |d| from bdl's bands: /hh/ 11.7 against 3.7; /f/ +16–28 dB at 0–2 kHz;
    - capping the flow at 150–250 cm³/s changed nothing.
    - Now it follows the glottal area's running minimum (15 ms release): full under 0.02 cm², none over 0.06. Vowel HNR unchanged.
  - **Frication 2.0** (was 0.8). Bands (mean |d|), before → after: /f/ 10.6 → 3.1, /s/ 8.5 → 6.9, /sh/ 13.4 → 6.5, /v/ 8.7 → 4.4.
  - Made articulatory defaults (`articulatory_layout()`): jet_flux 0.04, kinematic_jitter 0.5, kinematic_softness 0.0005, frication 2.0, piriform_scale 0.25.
  - **/dh/** closes fully for about 30 ms (a dental stop). Its +20 dB at 0.5–4 kHz is the release and voicing. Not pursued.
  - **Identification** (phone_test): 21.5 → **28.6%**. Fricatives 27 → 37%, nasals 18 → 30%, stops 25 → 33%, vowels 21 → 28%. Still 0%: ao, jh, y, w, uw, uh, ow.
  - m54 against the owner's recording at 8–11 kHz: −19.4 against −29 (stronger frication; his microphone may roll off). Watch for hiss.
  - The "faint ringing" is not yet located.
- `listen/m54_tip_at_rest_schwa.wav`.
- Owner on m54: "he sounds like perfect paul exactly", then: "it literally does not matter how you change the mouth. how are you generating the sound to begin with".
  - The answer: the kinematic glottis, a cosine wound at the score's pitch. Structurally DECtalk's source.
  - Owner: "the problem is diagnosed. and you dont know how to make an organic sounding system". Then: "the tube idea that i talked about. repeatedly. ... that has to be the vocal cord", the owner's own spec in GOAL.md (three pieces per section, air on every surface, energy only through feedback, no driven mass-springs).
  - Then: "it's meant to have a small random fluctuation in exactly one component of this — the sort that it shifts the fundamental lines. that adds the organic quality. the shift is muscle-related. on the scale of microseconds and a fraction."
  - **The fold tube** (`Cord::set_fold_tube`, `FoldTubeSettings`; speak `FOLDTUBE=1`, `FOLD=key=value,...`):
    - **Pieces.** Six slices along the flow, each the body's three pieces: an elastic cover (own longitudinal restoring, plus vertical tension to its neighbours: the mucosal wave), a mass body coupled under it, and the support, the vocalis.
    - **Muscles set material and rest, never motion.** Tension Q divides masses and depth and multiplies stiffness and length. The rest gap is a wedge (zero at the anterior commissure) from the lower edge's rest to the upper's.
    - **Stepping.** Eight steps per engine step (1.3 µs), symplectic Euler.
    - **Air.** Bernoulli pressure on each slice up to where the jet separates (the first slice beyond the narrowest that is 20% wider), the outlet's pressure above that; suction bounded at one transglottal drop. The flow is solved with the trachea and tract loads as before, with Poiseuille along every slice and the slices' air mass.
    - **The one fluctuation.** The muscle's support stiffness, rms 1%, drawn every fold step through one pole at 100 Hz.
    - In fold mode the tension, posture and breath variations, the three tremors and the breath flux are off.
  - **Measured** (held /aa/, 750 Pa; periods from interpolated flow crossings, to about 0.1 µs; `analysis/periods.py`):
    - It self-oscillates from the lungs alone, closes fully every cycle, peak flow about 400 cm³/s.
    - Pitch is 198 Q Hz, so the open-loop map is natural 198, per 198.
    - The physics alone repeats to a RAP of 1.8 µs. The muscle term (0.01, 100 Hz) adds about 1.4 µs (mean |dT| 4.3 → 6.2 µs). White at the fold step, even at 30% rms, it averages away.
    - The breath flux shook the folds (RAP 42 µs), hence off.
  - **Closure.** At the kinematic body's rests (0.011 / 0.006 cm) the flow fell to zero within 60 µs (MFDR 1.9e6 cm³/s²): a click every cycle in the vowels, 56/s held. At 0.030 / 0.015 the closure zips along the wedge: MFDR 8.6e5, 1.5/s held. Now the fold tube's voicing rests.
  - **Owner's sentence.** F0 median 124 against an intended 119 (range 98–145). Clicks 4.7/s, no thumps. Band HNR 20.9 / 14.4 / 8.9 / 5.5 (no breath flux).
- `listen/m55_blown_fold_tube.wav`.
- Owner on m55: "well, he does sound different. not so much perfect paul. there is still that bit of ringing, and some warbling - but not as much as before."
  - **The ringing.** A time-resolved narrow-line detector (`analysis/ring2.py`: 21 ms frames, bins above 3 kHz standing 12 dB over ±400 Hz) found steady lines at 8 kHz and 16 kHz. They were in m54 and m55, in about 40% of all frames, silences included; none in the owner's recordings.
    - With a 10-step control tick they moved to 9.6 and 19.2 kHz: a control-tick artifact (96 kHz / 12).
    - The cause: the lip radiation load's place and area were set once per tick and held (Speaker → `Tube::set_opening`).
    - Gliding the share alone and clamping it at each node change jumped the load: spikes at voicing onsets (max second difference 0.45 against 0.14).
    - Now `set_opening(position, area)`: the continuous place and the area glide every step (0.5 ms), and the node pair follows the glided place.
    - Also: wall pieces are recompiled every step while a node or its support moves (they were recompiled every 8 steps and once per tick, and held between), and the support's holding factor glides (2 ms).
    - Result (owner's sentence, fold tube): the 16 kHz line in 68 frames (was 662), 8 kHz gone from the top lines. A faint remainder (dry: 16 kHz 95, 8 kHz 59 of 1322 frames) from other per-tick setters.
    - Onset spikes 0.071 (m55 0.135), largest spike as m55. Clicks 6.0/s, no thumps.
  - **The warble.** Not located.
    - Measured: the fold vibrates down to 200 Pa lung pressure and at rest gaps up to 0.07 cm. In the vowels where a 5 ms autocorrelation pitch tracker lost the voice ("with", "quick", "jumps", "over the"), the cycle log shows regular cycles. The tracker's failure was not the voice's.
    - Fast pitch wobble (3–8 Hz, cents) was not discriminating between m54, bdl and the owner.
- `listen/m56_fold_tube_no_ringing.wav`.
- Owner, with a spectrogram of ours beside the owner's own: "see how in the word hello, on the left, yours has 5 vertical lines? figure out what causes them."
  - The lines were in "my name is": five internal boundaries (m|ay, ay|n, n|ey, ey|m, m|ih), every one a nasal closure or release, where our harmonics broke.
  - **Ruled out:**
    - lip-section crossings: no dips (mean +0.7 dB) with the outside rigidity or the quiet mask disabled;
    - dropouts: the 1 ms energy envelope moves at most 3 dB;
    - overall spectral change at nasal boundaries: 1.94 against bdl's 1.81 (`analysis/seams.py`).
  - **The cause: pitch steps.** The blown fold's pitch jumped 2.5–3% within one or two cycles at each nasal closure and release, as the tract's load changed. Every harmonic moved at once: a seam.
    - ARCTIC, at nasal boundaries, the largest change between frames 10 ms apart (`analysis/f0seam.py`), median / p90:

      | Voice | Median | p90 |
      |---|---|---|
      | bdl | 2.4% | 4.3% (vowel middles 1.4%) |
      | fold tube | 4.0% | 7.1% |
      | kinematic | 1.4% | 3.8% (vowel middles 0.97%) |

  - **Heavier tissue** (masses and stiffnesses scaled together, same pitch):
    - ×2: 2.8% (5.8%); ×2 with cover loss 0.2: 3.8% (8.3%); ×2 with vertical ×2: 3.6%;
    - **×3: 2.4% (5.4%), vowel middles 1.0%.** Adopted: cover 0.09 g, body 0.30 g, kc 36000, kcb 90000, ks 150000, kv 60000, contact 9e5 dyn/cm.
    - Pitch map recalibrated to 217 (99.5 / 115.0 / 139.2 Hz for 100 / 115 / 140).
  - **Closure.** At the wide rests (0.030 / 0.015) the heavier fold closed gently (2.7e5 cm³/s², open 75%): weak highs, band HNR 4–7 kHz −13 dB.
    - Back to 0.011 / 0.006: 4.3e5, open 0.65. Its kinked closures spread spikes, so the air's opening now closes through a soft contact 5 µm wide (`contact_softness`; the pieces' contact forces keep the plain geometry). Held /aa/ clicks 10 → 6/s, no bounce.
    - On the clicks2 detector, real voices give 0–1.1/s (owner 0.5–1.1, bdl 0–0.6); m54 (kinematic) 3.8.
  - **m58** (owner's sentence):
    - boundary steps 1.1 / 1.2 / 2.8 / 1.8%; F0 median 117 against intended 119; clicks 5.4/s (3.1 with aspiration 0); no thumps; ringing lines faint (16 kHz in 56 frames of 1551);
    - band HNR 20.1 / 12.4 / 6.9 / 0.7;
    - 6–11 kHz 10–15 dB over the owner's recording, from frication 2.0 (FRICATION=0 gives −32: but the owner's recording is another sentence with fewer sibilants);
    - aspiration 0 lifts 4–7 kHz HNR to 4.5 and halves the clicks.
- `listen/m58_heavier_fold.wav`.
- Owner on m58: "he's essentially irish. sounds irish. sounds like he has a terrible accent, or a speech impedement."
  - **Pitch against the plan, per vowel** (owner's sentence, CYCLES; `/tmp/es/pv.py`, `/tmp/es/reg.py`):
    - kinematic: mean +5 cents, sd 18;
    - m58: mean −39, sd 42. Stressed vowels −47 against unstressed −24, so stress lost 23 cents of its rise. The high accents fell short: IY1 −77, IH1 −66, AH1 −92 cents.
  - **Cause: the blown tube's pitch fell with breath.** Held /aa/ at 115 Hz, ear off: 116.0 Hz at 400 Pa, 110.7 at 1000 (−0.8 Hz per cmH2O; people rise 2–6). Stressed syllables carry more breath, so their pitch was pulled down: the melody's stress pattern was inverted. The elongation by swing (`Cord::set_elongation`, which fixed the same fault in the two-mass folds) was never set for the articulatory body. The calibration had also drifted to −45 cents at 750 Pa since the rests and soft contact changed.
  - **The fixes:**
    - `FoldTubeSettings::elongation` 70 /cm², reference 0.0026 cm². Held /aa/, 400 → 1000 Pa: 3.7 / 3.6 / 3.2 Hz per cmH2O at 100 / 115 / 140 Hz. The swing (cover's lowest slice, peak squared) was 1.3e-3 cm² at 400 Pa and 4.3e-3 at 1000.
    - The pitch map was recalibrated at 650 Pa: natural 197, per 171 (99.9 / 115.1 / 140.4 Hz).
    - **The muscle anticipates its breath** (`swing_per_pascal` 4.9e-6, `swing_onset` 130 Pa). The open-loop tension is divided by the stretch the lungs' pressure (`lung_`, the reservoir, not the pressure across the folds) will bring. Without this the planned contour, which already contains the stress rise, overshot: mean +5, sd 58, stressed +17 over unstressed, error against lungs r = +0.71. With it, held /aa/ at 115 Hz runs 114.4 → 120.7 Hz from 400 to 1000 Pa, the remainder being what the tract does to the pressure.
    - **The ear in fold mode** (`speak` FOLDTUBE sets `pitch_feedback` 0.08). `Cord::contact()` now reports the fold tube's most closed slice; it had read the unused two-mass edge. The ear hears only lungs over 250 Pa and periods within half an octave of the intended one.
  - **m59** (owner's sentence): pitch error mean −2, sd 40; stressed +1 against unstressed −8 cents; F0 median 119 (intended 119). Band HNR 19.7 / 12.9 / 6.7 / 0.9; clicks 5.4/s (as m58).
    - Nasal-boundary steps 5.9 / 1.2 / 2.2 / 3.0 / 2.0%. The first is the voicing onset after the pause (107 → 111 → 118 Hz as the swing builds).
    - Not addressed: the "speech impediment" part (articulation; heard since m54's freed tongue tip).
- `listen/m59_breath_raises_pitch.wav`.
- Owner on m59: "okay. lets fix body. because he sounds like he's saying bady".
  - **"body" /aa/ F2** (LPC, `analysis/ftrack.py`):
    - bdl, five ARCTIC tokens (a0178, a0256, b0104, a0260, b0195): F1 about 700–780 Hz; F2 1100–1250 from the /b/ release, rising toward 1300–1450 only in the last third.
    - m59: F1 745 Hz; F2 1500 at the release, 1413 at best, then 1870 by the /d/. Toward /ae/.
    - The tuned /aa/ posture held alone gives F2 1243: the posture is right, the motion is not.
  - **Cause 1: the closure's pharynx pull lingered.** A held /aa/ with one lip closure and nothing else moving: F2 1420 after the release, back to 1240 only 150 ms later. PROFILE: 30 ms after the release the 2–10 cm areas were 1.5–2.3× their rest; 100 ms after, still +30–40%.
    - The source is `Speaker::shape_body`'s voiced-closure pull (0.25 cm outward on 3.5–10 cm), restored with a 40 ms time constant (its comment said about 30 ms).
    - Now 10 ms. Held test: F2 1311 at 45 ms (was 1377); in "body", F2 at best 1346.
  - **Cause 2: the jaw under the /b/.** Capped at 5 degrees (the /aa/ jaw is 12.8), it lifted the tongue: F2 near 1500 at the release. The lips seal at any jaw angle (tested at 12.8). Cap now 8 (`kLabialJaw`).
  - **Cause 3: across the flap the next vowel began 60 ms ahead.** Here that was 27 ms into the 87 ms /aa/: jaw, lips and body turned toward /iy/ before the /aa/ arrived. A vowel now holds for 60% of its length before the tongue heads for the next (`kVowelHold`); the flap's own lead still applies.
  - **m60** (score regenerated: `say ... --organs ../tract/data/organs.txt`, identical to k4 apart from these timings):
    - "body" /aa/ F2 1465 → 1274 → 1600 (bdl 1100–1250 → 1300–1450). The release is still about 200 Hz high.
    - clicks 5.0/s; band HNR 19.7 / 12.8 / 6.9 / 1.0; F0 median 119; boundary steps 4.8 / 1.2 / 2.4 / 1.0 / 1.5%.
- `listen/m60_body.wav`.
- Owner on m60: "the warble is greatly diminished. lets work on hello. it doesnt sound like hello. it sounds like halting- hetlo".
  - **The /l/ was a stop.** TRACE: 0.80–0.845 s the narrowest area 0.006–0.015 cm² at 13 cm, flow through the lips near 0; the level fell 38–43 dB.
    - PROFILE at 0.815: the side channels (0.35 cm² × lateral) covered only cuts within 0.8 cm of the tip. The blade 0.81 cm behind it lay 0.09 cm from the palate: 0.0057 cm².
    - Now the channels run from 2 cm behind the tip contact to 0.8 cm ahead (`MidsagittalTract::area_of`).
  - **Then the /l/ was a vowel** (F1 525 Hz, no dip, F2 1600).
    - bdl's /l/ between vowels (a0001–a0100, 46 tokens, `/tmp/es/lstat.py`): F1 414, F2 957, F3 2784 Hz, 2.8 dB under the next vowel (p25 −7.7).
    - A grid at "hello"'s /l/ (side channels × body), F at 0.80–0.85 s:
      - channels 0.35, 0.21 and 0.12 cm² (lateral 1, 0.6, 0.35);
      - body at /ow/'s (2.97, 7.29): F2 1610–1690;
      - body at (2.65, 6.48): F2 1250–1270;
      - F1 jumps from about 290 to 525–600 between 0.12 and 0.21 cm².
    - Adopted: channels 0.12 cm² (`kLateralArea`). Beside a central or back vowel the planner's /l/ body goes to (2.65, 6.5) (`kDarkLateralX/Y`); beside front vowels it keeps the 0.2 cm retraction.
    - Delaying /ow/'s lip rounding to the /l/'s end lifted F3 from 2140 to 2310; not adopted.
  - **m62:**
    - "hello" /l/: 9 dB under the vowel, F1 260–310, F2 1450 → 1185, F3 2400 → 1900 as the rounding arrives. F3 sits well under bdl's: an /r/ colour is possible.
    - "lazy" /l/ (front context): 8–11 dB down (was 5, a vowel), F2 about 2000–2100 (clear; bdl's 957 is dark).
    - clicks 5.6/s; band HNR 19.8 / 13.5 / 5.9 / 0.2.
- `listen/m62_hello.wav`.
- Owner on m62: "between fox and jumps there's a weird vocalization".
  - **The cause: a voiced vocoid between /s/ and /dZ/.** The spectrogram shows voicing pulses with vowel formants at 5.69–5.715 s, after the /s/ hiss and before /dZ/'s closure.
    - The planner voices 15 ms before a fricative's end (right before a vowel). Here the glottis closed at 5.687 while the tip only sealed for /dZ/ at about 5.716: 25 ms of voice through an open mouth.
  - **The fix.** After a fricative, a following closure gets the voice only once it has sealed (`kSealSettle` 25 ms after the fricative's end). Before an affricate the fricative emits no voicing at all: a voiced /dZ/ opens for its hiss and re-closes the folds at its own end, so it is devoiced across "fox jumps" as in English.
  - In this sentence only that one glottis event changed. The glottis now stays at 0.19–0.21 cm² from the /s/ through the /dZ/ closure.
- `listen/m63_fox_jumps.wav`.
- Owner: "you can commit these changes. okay. lets do the SNR comparison to NATO and work to refine." Committed and pushed 0f9b6fe (`claude/segment-voice`, `ducklips_claude/`).
- **Band signal-to-noise against NATO** (`analysis/nato_snr.py`). The 26 NATO words, rendered by the main line (say --organs, speak FOLDTUBE=1), against the NATO reading (`ref/nato/`) and the owner's two passes (`ref/user_nato.wav`).
  - **The method.** Every glottal cycle in the loud stretches (within 20 dB of the word's peak) is found by waveform matching, with none dropped for irregularity. Each is compared with the next band by band: HNR = 10 log10(r / (1 − r)). Medians per word.
    - Gated to clean cycles (match > 0.9, as `bandhnr.py`), ours scored 2–10 dB over both references: the gate kept only our best cycles. pyin-based cycle finding (wmatch) found no voice at all in 6 of our words (Delta, Kilo ...), though the fold log shows cycles.
  - **Recording floors.** The owner's recording is 63 dB under the speech in every band but 4–7 kHz (49): not limiting.
  - **Jitter.** Over clean cycles, RAP median 0.34% (NATO 0.29, owner 0.26), p90 1.78% (0.97, 0.84); shimmer 0.56 dB (0.56, 0.60). Our tail is onsets: about 3% for the first two cycles after each voice onset, under 1% from the third.
  - **Energy-weighted parts** (`/tmp/es/hn.py`): per band, harmonic r·E and noise (1 − r)·E summed over loud cycles, dB re the lowest band's total, median over words. At 1–4 kHz ours matched within 0.5 dB. At 4–7 kHz ours had 7 dB too much of both parts:

    | 4–7 kHz | Harmonic | Noise |
    |---|---|---|
    | ours | −23.9 | −22.0 |
    | owner | −30.9 | −28.2 |
    | NATO | −33.4 | −27.6 |

  - **Levers, at 4–7 kHz:**
    - closure softness (FOLD soft 0.002–0.004) tilts all bands: 2–4 kHz under the owner before 4–7 matches;
    - aspiration moves only the noise;
    - cover loss 0.3 is worse;
    - **piriform recesses at full size** (their 4.5 kHz zero): harmonic −29.6, noise −26.4.
  - **Adopted: piriform ×1, vowels retuned** (organ_vopt3, W 10, PIRISCALE=1, then a second pass for iy uw ao ah, MAXFEV 700). Formant errors at full size before retuning:
    - most F1 and F2 down 0.5–2 st; /aa/ F2 −2.0, /ao/ F2 −4.5.
  - Retuned: within 0.5 st: ih eh ae aa uh ey ax; ah 1.4, ow 1.2, iy 1.4 (F2 2150, stuck), uw 2.0 (F3 2244), ao 2.9.
    - Shape errors: aa 0.50 → 0.57, uh 0.42 → 0.69, others within about 0.05.
  - **Now, against the owner** (per-cycle medians, ours minus owner): +1.8 / +3.1 / +0.6 / +0.6 dB (was +2.1 / +4.0 / +2.0 / +1.5). Clean cycles 61% (owner 66). Energy-weighted harmonic within 1.7 dB in every band; noise at 4–7 kHz −24.5 against −28.2.
  - **Aspiration kept at 8.** The breath noise is concentrated in the loudest, highest-flow cycles (v³), the owner's more evenly; lowering the gain cannot fix that.
    - At 5: energy-weighted noise within 0.7 dB in every band, but the typical cycle too clean (per-cycle +3.3 / +4.1 / +3.4 / +3.2).
    - /hh/ above 2 kHz re the loudest 20 ms (`/tmp/es/hhlev.py`): bdl −30.5 (186 tokens; p25 −35.9). Ours ("hello"): m63 −28.0; piriform ×1 −31.3 at aspiration 8, −33.7 at 6, −35.3 at 5.
  - **m64** (owner's sentence): "body" /aa/ F2 now reaches 1172 Hz (m60 1274; bdl 1100–1250). Band HNR 19.4 / 12.8 / 6.4 / −0.9. Clicks 6.9/s (m62 5.6): OPEN, possibly the weaker steady sound above 4 kHz exposing transients.
  - **Next:** the 1–2 kHz per-cycle excess (+3 dB); the noise's dependence on flow (more even through the cycle and from cycle to cycle).
- `listen/m64_full_piriform.wav`.
- Owner: "lets try turning this into an blockwise AST system and finetuning to match the human voice generations, ie, a gradient descent behavioral finetune ... we want to dampen shifting or make it so massive shifts take a long time, too, so it doesnt explode."
- **The parameter tree** (`include/tract/parameter_tree.hpp`, `src/parameter_tree.cpp`). The engine's tunable settings are declared once as a tree of named blocks: breath, fold, posture, ear, variation, glottal_noise, oral_noise, walls, branches, tract, mouth, kinematic (74 leaves; 42 are free in the fit).
  - Each leaf is bound to a field of the body's layout (applied before the body is built) or a setting of the speaker (applied after). An undeclared path is an error.
  - `speak` reads a block file with `PARAMS=file` (over the defaults and FOLDTUBE, under the other switches); `PARAMS_OUT=file` writes the tree as spoken.
  - New layout fields `wall_stiffness_scale` and `wall_loss_scale` scale every vocal wall's support. The speaker gained single-value setters and getters.
  - `data/voice.params` is m64's tree. **MEASURED:** dumped and read back, Alfa renders bit-identical.
- **The fit** (`analysis/fit_blocks.py`). Block-coordinate SPSA in normalized units (log or linear per leaf, bounded): two renders of all 26 NATO words per gradient (same seed), Adam-like moments.
  - **Damping.** Each step is capped by its block's trust radius (≤ 0.1 units). A spring pulls every leaf back toward its start (ANCHOR). A step is rendered before it is kept: if the loss does not fall (or a render fails, goes non-finite or loses its voice) it is undone, the radius halved and the momentum cleared.
  - A second seed checks for overfitting every 10 iterations. One evaluation (26 words, rendered and measured) takes about 4 s on the Mini.
- **The behavioural loss failed validation; nothing was fitted.** The loss: per word, band HNR medians and spreads, RAP and shimmer, voiced band balance, pitch span, clicks and DTW level contours, each z-scored by how much the three human readings differ, through a Huber function. The test: the human readings scored against the other two (leave one out) must beat our sources.

  | Scored against | Loss |
  |---|---|
  | humans, leave one out | owner 1: 0.77, owner 2: 0.87, NATO: 0.76 |
  | m64 fold tube, same references | 1.17–1.23 |
  | no aspiration, same references | 1.09–1.15 |
  | **kinematic source (Perfect Paul, rejected)**, same references | **0.72–0.78** |

  - Against all three readings: kinematic 0.66, m64 1.09.
  - The rejected formula source scores like a person, better on hnr, perturbation, contour and span. Those statistics were matched onto the kinematic source in m51–m52. Scoring distance from the humans' mean also rewards averageness: a voice at the mean of three readings beats any one real reading (span: humans 1.07–1.32, kinematic 0.58–0.87).
  - **REJECTED as a steering loss.** It would steer toward Perfect Paul.
- **Exploration** (13 per-word candidates, scratch `explore.py`): high-band level against low-band level per cycle (slopes 1–2, 2–4, 4–7 kHz), pitch–level correlation and slope, lag-1 correlation of period and amplitude changes, onset cycles and onset/offset pitch, cycle-to-cycle tilt flicker, level and tilt ranges.
  - For every one, a human reading falls outside the other two readings' range (widened by half) at least as often (6–18 of 26 words) as the kinematic source falls outside all three (2–9). No summary statistic of 26 words separates people from the formula source.
- **OPEN:** a steering signal that does separate them. CANDIDATE: a pretrained audio critic (AudioSet's Audio Spectrogram Transformer has "Speech synthesizer" beside "Male speech"); it must first rate the owner, NATO and bdl human and the kinematic source synthetic. Needs a model download (owner's permission).
- Owner: "alternating least squares, is what i meant. blockwise." (not an abstract syntax tree; the parameter tree stays as the fit's address space).
- **The ALS fit** (`analysis/fit_blocks.py`; the SPSA version is kept as `analysis/fit_spsa_summary.py`, REJECTED; block files in `analysis/parameter_files.py`).
  - **Residuals, paired.** For each NATO word, each of the owner's loud frames (10 ms hop, 25 ms window, within 35 dB of the word's loudest) has eight features: four band levels (dB re the loudest frame) and four band periodicities (HNR one period apart).
    - Our frames are aligned to the owner's by DTW on the band levels. Each owner frame is met by the mean of our frames aligned to it.
    - Residual = (ours − owner) / scale, plus one per word: the median |cents| between heard and intended pitch, divided by 20.
  - **Blocks, alternated.** The alignment is its own block: solved by DTW at the start of each parameter block's turn, then held, so the residual vector keeps its length.
    - Each parameter block in turn: a forward-difference Jacobian (k + 1 renders of 26 words), then a Levenberg–Marquardt step on mean(r²) + ANCHOR·Σu².
  - **Dampers:**
    - λ per block (×4 per rejected step, ½ per kept one, floor 0.05);
    - each |d| capped at the block's radius (≤ 0.1 units);
    - the anchor's spring toward the start;
    - render-before-keep;
    - a second seed is scored after each sweep.
  - **Smoke test** (one sweep, stopped after two blocks): breath 1.1152 → 1.1137 kept. Fold: the first step (λ 1) raised the cost to 1.1208 and was undone; at λ 4, radius 0.025, 1.1153 kept. The mechanics work; gains of 0.1% are at noise level.
- **Data fault: the second half of `ref/user_nato.wav` is not the alphabet in order.**
  - `nato_compare.OWNER_SPANS[1]` (31.9–62.3 s) labels its 26 segments Alfa…Zulu, and the segment boundaries are right. But matched by MFCC DTW against the owner's own pass 1, only "Charlie" lands on its label.
  - The same matcher puts the NATO reader's words (a different speaker) on the right owner word 12 of 26 times.
  - The recording has 39 segments after 31.9 s. Before this fault was found, per-frame band levels differed 10–17 dB rms between the two "passes" in the loud frames.
  - The m64 "owner" SNR column (`nato_snr.py`) averaged both passes, so half of it is other words. OPEN: what the owner read there.
- **Validation, against pass 1 only** (scales from the faulty pass 2, which scales every row alike; mean squared residual):

  | Candidate | Loss |
  |---|---|
  | NATO reader (another person) | 1.024 |
  | no aspiration | 1.009 |
  | m64 | 0.856 |
  | kinematic (Perfect Paul) | 0.802 |

  - The rejected source again lands nearest the owner, on periodicity: 0.857 against m64's 0.960. The fit would chase the same per-frame noise statistics that made Perfect Paul.
  - **Not run beyond the smoke test.**
- Owner: "listen. what you're doing is tuning the model to make samples that are closely matched to the real samples. okay?" then "no. not my voice. real reference samples."
  - **GIVEN:** the objective is closeness to real reference recordings; no human-versus-rejected-source gate.
- **The fit retargeted to the NATO reference reader** (`ref/nato`, 26 words; the owner's recording is no longer a target).
  - **Frame features (17).** Spectral shape: twelve log-spaced bands, 200 Hz–7 kHz, re the frame's own level, floored at −40 dB. The frame's level re the word's loudest, floored at −50. Periodicity of four bands.
  - **Why the floors.** The reference is a processed recording: noise-gated to −120 dB in pauses, with quiet bands emptied to −100 dB as a lossy coder does (Alfa's "fa": 3.9–5.2 kHz at −94 to −113 dB). Unfloored, those frames dominated: 15–20 dB rms per band.
  - **Why the trim.** Both sides are trimmed to their loud span (35 dB, 3-frame margin) before DTW. Our renders carry about 0.6 s of silence the reference words do not, and DTW, forced end to end, paired reference speech with it.
  - Scale 3 dB for every feature.
  - **m64 against the reference (MEASURED):** ours minus reference, mean per band 200 Hz → 7 kHz:
    - +12.6 / +11.2 dB at 200–360 Hz;
    - −3 to +4 dB at 0.5–2 kHz;
    - +6.3 / +9.5 / +8.0 dB at 2.1–5.2 kHz;
    - per-frame scatter about 10 dB per band (vowel formants of two speakers).
  - **Costs:** m64 14.96 (shape 18.0, level 12.4, periodicity 6.7, pitch 1.9); kinematic 15.50.
- **Render noise against the residual** (MEASURED, segment-pooled):
  - m64's rms residual against the reference: 2.93 units (3 dB each);
  - the same tree at another seed changes it by 0.54;
  - fold.cover_mass × 1.001 changes it by 0.12.
  - Frame by frame (before pooling) every step of the first ALS run was undone: a nudged parameter re-draws the folds' cycle-to-cycle course, so per-frame differences were mostly chance. Residuals are now averaged over each third of a word's loud frames, with a 0.3-unit difference step.
- **Trust radius widened, 0.05 → 0.2 units (max 0.4).** At 0.05 the gain per step (about 0.005) sat under the chance in the cost (0.01–0.03), and the dampers stalled the fit. With unit 0.2 on a log scale, 0.2 units is 4% per turn.
  - The linear model then predicted the kept steps: fold 8.511 predicted against 8.560 rendered; posture 8.530 against 8.518.
- Owner: "can try adam instead if desired. blockwise ASL is just my 'safe bet' when working with unknown problems". ALS kept: it converged.
- **The fit: four sweeps (MEASURED).** The training cost fell 8.578 → 8.048 over four sweeps; the untrained seed 51117 was scored after each.

  | Cost | m64 | Fitted |
  |---|---|---|
  | seed 73021 (fitted on) | 8.578 | 7.990 (−6.9%) |
  | seed 51117 (never fitted) | 8.905 | 8.246 (−7.4%) |

  - Groups at 73021: shape 11.19 → 10.36, level 4.94 → 4.63, periodicity 2.22 → 2.25, pitch 1.93 → 1.70.
  - The branches block gave the largest gain (−0.09 per sweep).
- **What moved** (`data/fitted_m65.params`):
  - piriform port 0.6 → 0.42 cm², piriform scale 1.0 → 1.11, piriform loss 0.5 → 0.43, sinus necks 0.25 → 0.20;
  - aspiration 8 → 6.5, fold turbulence 1.0 → 0.87;
  - tract hop loss 0.02 → 0.015;
  - walls: support loss × 0.81, stiffness × 0.94, skin corner 250 → 324 Hz;
  - fold tissue within 1–2% of its start (λ climbed to 2000 there: its steps did not pay).
- **The PARAMS round trip** is bit-identical on the sentence (data/voice.params against no PARAMS).
- **m65, the owner's sentence:**
  - band HNR 19.6 / 12.2 / 6.5 / −2.7 (m64 now: 18.8 / 12.3 / 6.3 / −1.3);
  - clicks 7.1/s (6.2).
- `listen/m65_als_fit.wav`; `listen/m65_before_m64.wav` (m64 rendered by today's code, for A/B).
- **Next:** more sweeps, with the organs (vowel postures, organs.txt) as blocks. Shape, where the remaining cost lies, is mostly formants, which only the organs can move. Then sentences as well as words.
- Owner on m65: "honestly? it sounds completely the same."
  - The fit moved few parameters, by a few percent, on a coarse measure. The fold tissue never moved: its λ climbed to 2000.
- Owner: "build a simple neural network where the neurons track parameters. all of them. and we design it so that it doesnt explode too fast anywhere. then we can try ALS or adamw on it. we're not training a neural network- we're running a policy optimizer to finetune an engine."
- **The engine network** (`analysis/engine_net.py`). 186 neurons in 24 blocks: 69 tree leaves, plus 13 vowels × 9 posture coordinates.
  - **Left out (5):** the four kinematic.* leaves, which act only on the kinematic source, and glottal_noise.flux_corner, whose 0 is a switch.
  - **Activation.** Each neuron has a latent z (0 at the start) and a bounded, saturating activation:
    - log or linear: sigmoid between bounds, biased so z = 0 gives the start;
    - "onset": high·tanh(max(z, 0)), for causes that start off (the variations, closing turbulence, jet flux ...).
  - **Bounds:**
    - tree leaves ×/÷2.5 by default, explicit ranges for 22;
    - vowels: jaw ±4°, tongue body and tip ±0.6 cm, aperture ×/÷2, protrusion −0.1 to +0.5, width ±1, larynx −0.2 to +0.5.
  - Every candidate is re-planned by `say` with its own postures (under 10 ms a word).
  - **MEASURED:** at z = 0 the network reproduces m64's cost exactly (8.5776).
- **The policy optimizer** (`analysis/policy_fit.py`).
  - **adamw:** antithetic ES over all neurons. 6 pairs per step, σ 0.2; the same seed within a pair, rotating over 8 seeds across pairs. Rank-shaped costs feed AdamW (rate 0.03, β 0.9 / 0.99).
    - Dampers: each latent's step clipped at 0.05; decoupled decay 0.01 toward the start; saturating activations; the policy rendered at two seeds every 5 steps, rolled back to the best if 3% worse (rate halved, moments cleared).
  - **als:** the blockwise Levenberg–Marquardt of fit_blocks.py in latent space, over the 24 blocks.
  - **First adamw steps:** population mean 9.85 → 9.30 in 5 steps. The policy (mean of seeds 73021 and 51117) went 8.741 → 8.546, not yet at the four-sweep ALS fit on the body (8.118 on the same two seeds).
- **adamw, 80 steps (MEASURED):** policy cost (mean of seeds 73021 and 51117) 8.741 → 7.219 (best, step 70; −17%). The body-only ALS fit reached 8.118. Still falling: population 9.85 → 7.43. The rms latent 0.32; no rollbacks.
  - **Largest moves:**
    - fold rest gaps lower 0.011 → 0.018 cm (+67%), upper 0.006 → 0.010 (+70%);
    - fold body mass +49%; cover loss −27%; coupling and vertical stiffness −25%; elongation −27%;
    - motor wander −27%; skin gain −21%; pitch feedback +25%;
    - tension variation switched on (0.013); piriform smoothing on (0.13);
    - flare widest +29%; ER jaw +26%; vowel tongue shifts up to 0.2 cm.
  - **m66** (owner's sentence, `listen/m66_adamw_policy.wav`; `data/policy_m66.params`, `data/policy_m66.organs`):
    - band HNR 20.7 / 12.5 / 5.0 / −7.9 (m64: 18.8 / 12.3 / 6.3 / −1.3);
    - clicks 9.8/s (6.2);
    - pitch against intended: median +16 cents, sd 56 (m59: −2, 40).
    - The fit's measure does not see clicks; they rose.
  - Resumed for 120 more steps.
- **adamw continued** (killed at step 95 for the listener work): policy 7.061 at step 95.
- Owner on m66: "funny. he sounds about the same."
  - The cost fell 17%, but it is twelve band levels, a level and four periodicities per word third: not what the ear hears.
  - Owner approved downloading listener models ("yes").
- **Pretrained listeners, validated before steering** (`analysis/listeners.py`, `analysis/listener_check.py`).
  - The AASIST clone was refused by the session's safety policy (it would bring in and run the repository's code). Only stock-transformers classifiers were taken: config plus safetensors, in `~/models` on the Mini.
  - `transformers` must be kept off the Mini's broken jax (USE_FLAX=0).

  | Listener | Result (REJECTED) |
  |---|---|
  | AST ASVspoof5, AST ASVspoof2019 (MattyB95) | P(synthetic) 1.000 for every set, people included (NATO reader, owner, bdl) |
  | wav2vec2 mo-thecreator | the NATO reader and bdl "fake" (1.000); kinematic words "real" (0.000) |
  | wav2vec2 MelodyMachine V2 | ARCTIC sentences: bdl −11 log-odds, ours above (AUC 1.00), kinematic +7.7, m64 −4.1. But six of seven owner phrases come out "fake" (+10.6; the owner's whole sentence −11), one rms sentence +10.7, the NATO reader's words "fake" |

  - Each is keyed to its training recordings: ASVspoof and ARCTIC conditions, not people.
- **WavLM paired distance** (`analysis/wavlm_objective.py`; microsoft/wavlm-base-plus). Our ARCTIC sentence against bdl's real recording of the same sentence: trimmed, normalized, layers 6 and 9 (unit frames), DTW on cosine distance, mean along the path.
  - **MEASURED (a0001–a0010), layer 6 / 9:**

    | Candidate | Distance to bdl |
    |---|---|
    | rms (another person) | 0.283 / 0.156 |
    | m64 | 0.588 / 0.489 |
    | kinematic | 0.564 / 0.469 |
    | m66 | 0.571 / 0.473 |

  - rms is nearer bdl than each of ours in 10 of 10 sentences. Kinematic is nearer than m64 in 9 of 10: the measure sees real speech against ours, not the source specifically.
  - Adopted as the policy's cost (GIVEN: match real samples): `policy_fit.py WORKDIR wavlm`. Sentences a0001–a0010 steer; a0011–a0020 are held out. About 35 s a step.
  - Start: 0.5385 fit, 0.5539 held.
- **WavLM policy, 60 steps (MEASURED).** Distance to bdl: fit sentences 0.5385 → 0.4935 (−8.4%); held-out a0011–a0020 0.5539 → 0.5224 (−5.7%). Still falling. A person (rms) sits at about 0.22.
  - **Largest moves:**
    - jet decay +52%;
    - fold cover loss 0.10 → 0.068; vertical stiffness −23%;
    - turbulent glottal flux switched on (jet_flux 0.056; off since m50, when it shook the blown folds);
    - tremor 0.005, tension variation 0.015, posture variation 0.0014 switched on;
    - wander −20%;
    - ER tongue body −0.24 cm and forward 0.2; UW aperture +25%.
  - `listen/m67_wavlm_policy.wav`; `data/policy_m67.params`, `data/policy_m67.organs`.
  - **m67 (owner's sentence):** band HNR 20.8 / 1.5 / −1.3 / −13.9 (m64: 18.8 / 12.3 / 6.3 / −1.3): much noisier above 1 kHz, the turbulent flux. Clicks 2.2/s (m64 6.2).
- **WavLM's m67, for the owner's sentence:** "dog sounds as irish as ever". The owner asked for the real comparison, NATO against NATO, interleaved (`listen/nato_interleaved_m67.wav`). The owner flagged ten words: Charlie, Foxtrot, Golf, Juliett, Lima, Sierra, Uniform, Whiskey, X-ray, Zulu (`listen/nato_focus_m67.wav`).
  - In the owner's words: "ch-charlie (stutters). foxtrot warble. golf is most damaged of all. juliet stutters. llima. sierra sounds shierra. uniform reads oohnihform. whiskey reads wuiskey. xray reads like a boxer said it ie ekh-rahy".
- **Spectrograms against the reference** (`/tmp/es/specpairs.py`):
  - m67's vowels are filled with noise to 7 kHz (the policy's jet flux), where the reference's harmonics are clean. The work below runs on m64's tree.
  - Our starting consonants are 0.15–0.2 s; the reference's 0.05–0.1.
  - "Golf" opens with 0.25 s of near-silence: the /g/ closure with prevoicing blips.
- **Fixed: affricates never sealed** (planner.cpp, the organ branch for SH ZH CH JH).
  - The groove (0.6) was set at the closure's start; the area rule's groove floor (0.06 × 0.6 = 0.036 cm²) kept the tip open. The hiss leaked from the start, and the release burst inside it: a transient at 0.74 s as loud as the vowel in "Juliett".
  - Now the groove forms with the hiss (at 35% of the phone). **MEASURED:** the tip is sealed (0.000 cm²) 0.675–0.725 s, then releases once into the hiss.
- **Fixed: /k/ held through /s/** (planner.cpp, the hand-back).
  - A consonant handed its organs back to the vowel only if no later consonant before the vowel owned them. In "X-ray" (K S R) and "Foxtrot" (K S T R) the /r/ owns the body, so the /k/ never released: the tract was shut 0.805–0.885 s through the whole /s/.
  - Now a closure (stop, nasal, affricate) is kept only by the consonant right after it. **MEASURED:** frames of 4–8 kHz hiss within 20 dB of the word's peak: X-ray 0 → 4, Foxtrot 0 → 10.
- **OPEN: /s/ too low ("shierra").** Peak against the NATO reference:

  | Word | NATO | Ours |
  |---|---|---|
  | Sierra | 7.5 kHz | 3.1 kHz |
  | Whiskey | 6.7 kHz | 3.6 kHz |
  | X-ray | 5.8 kHz | 3.9 kHz |
  | Foxtrot | 5.3 kHz | 7.3 kHz |

  - **Ruled out:**
    - lip aperture floor 0.5 → 1.4 cm (centroid +1 kHz, peak unmoved);
    - groove taper 2.2 → 0.7 cm;
    - body lowered 0.4–1.2 cm during /s/;
    - constriction-to-teeth distance (1.4 cm in both Sierra and Foxtrot).
  - `speak` PROFILE now prints the teeth arc.
- `listen/nato_focus_m68.wav`: the ten words, NATO then ours (m64 tree, both fixes).
- Owner on m68: "they have a gravely seal. its still not fixed really".
  - **The seal itself is clean** (TRACE 0.5 ms): the tip closes smoothly to 0.000 cm², with no chatter.
  - **The fault is after it.** At the hiss-to-vowel transition one spike, the word's largest sample, then full-strength voice at once. The NATO reader's "Charlie": burst, about 40 ms of hiss, voice swelling over about 70 ms; "Juliett" voiced through the hiss.
  - **Cause** (TRACE 0.25 ms, Juliett 0.826–0.833 s): the folds swing from apart (0.08 cm²) to shut in about 6 ms under 200 cm³/s while the hiss channel opens. The oral pressure collapses 430 Pa → 0, a hard glottal attack.
  - Phrase-initial voiced closures wait with the folds apart (kWaitingGap), and the affricate parted them further for its hiss (affricate_hiss_gap 0.04), so the whole breath stood behind the channel.
- **Fixed (planner.cpp):**
  - An affricate's hiss lasts at most 60 ms (`kAffricateHiss`, `affricate_hiss_start()`): hiss 125 → 65 ms (Charlie), 90 → 45 (Juliett); NATO 40–60.
  - A voiced affricate closes the folds in from its hiss over 60 ms (`affricate_voicing_move`, replacing `affricate_hiss_gap`; say `--affricate-voicing-move`). Juliett's largest sample went from 1.30× the vowel's following peaks to 1.13×; crest 3.9 → 3.3 (NATO 2.3).
- **Not fixed: Charlie's hard onset after voiceless /tS/** (crest 4.1, NATO 2.0). Re-tested on the fold tube, crest per word:
  - adduction move 0.05 s: worse or equal in 7 of 8 words;
  - two-stage onset (`--onset-gap` 0.015): worse or equal;
  - both: worse.
  - Our crest is 3.1–4.1 in every word with a voiceless onset (NATO 2.0–2.8). OPEN, the old "hard voice onsets after /p t k/".
- `listen/nato_focus_m69.wav`.
- Owner: "let's do a timing and articulation gradient descent pass, just as we did a body shape finetune."
- **The planner's parameter tree** (`voice/include/voice/plan_parameters.hpp`, `voice/src/plan_parameters.cpp`; `say --plan-params FILE`, `--plan-params-out FILE`). 167 leaves:
  - timing (22), glottis (17), breath (10), release (2), pitch (12), mouth (8, the area planner's);
  - articulation (55): `ConsonantShapes`, the organ planner's consonant shapes, until now literals in `gesture_organs`;
  - duration (41): a factor on each phone's planned length, new: `PlanOptions::duration_scale`.
  - The constants kSealLead, kVowelHold, kSealSettle, kLabialJaw, kDarkLateral*, kAffricateHiss became options and shapes.
  - The block-file parser moved to `voice/block_file.hpp`, shared by `say` and `speak`.
  - **MEASURED:** after the refactor, the ten flagged words' scores and audio are byte-identical to m69. A dumped plan tree read back gives an identical score. `voice/data/plan.params` is the start.
- **engine_net.py, groups.** The network is built over the groups being fitted: body (69), organs (117), plan (144; left out: pitch.*, mouth.*, onset_gap, lead_in, pause_pressure). `GROUPS=organs,plan` fits timing and articulation with the body at m64.
  - Plan ranges: durations ×/÷2; tongue positions ±0.4 cm; presses ±0.15; jaw ceilings ±2°; moves ×/÷2.
- **The cost gained timing.** DTW forgives timing, so each word adds its length (log of our loud span over the reference's, / 0.1) and its inner timing (at 5 reference frames, the aligned position against the reference's, s / 0.03). Frame features are pooled over 6 runs per word.
  - **Start (MEASURED):** cost 10.40; shape 13.5, level 7.2, periodicity 3.4, duration 9.6, timing 4.5, pitch 1.9.
  - Word lengths, ours over NATO: Golf 2.56, Victor 1.73, Kilo 1.58, November 1.55, Charlie 1.50, Uniform 1.43, Mike 1.40, Quebec 1.39; the rest 0.89–1.32.
- AdamW, organs + plan (261 neurons), 80 steps: running.
- **Timing and articulation pass, adamw, 80 steps (MEASURED).** Policy cost (two seeds) 10.503 → 9.815 at best (step 40, −6.6%); flat within the seeds' chance (±0.3) after.
  - **Largest moves** (rms latent 0.30): durations /uw/ ×1.28, /w/ ×1.21, /ah/ ×1.13; /w/ body share of /uw/ 0.5 → 0.62; vowel postures within 0.2 cm (EY aperture +25%, ER jaw +0.7°); adduction lead +4.6 ms.
  - **Word lengths, trimmed** (NATO / m69 / m70):

    | Word | NATO | m69 | m70 |
    |---|---|---|---|
    | Charlie | 0.45 | 0.70 | 0.68 |
    | Golf | 0.39 | 0.74 | 0.71 |
    | Lima | 0.44 | 0.66 | 0.51 |
    | Sierra | 0.53 | 0.60 | 0.50 |
    | Juliett | 0.61 | 0.76 | 0.89 |
    | Uniform | 0.59 | 0.98 | 0.97 |
    | Zulu | 0.46 | 0.69 | 0.73 |

  - **The systematic part, about 1.4× too long, barely moved.** With 261 neurons, 6 pairs per step give a noisy gradient. With the step cap (0.05 latent) and decay, a 40% change in rate (latent ≈ 1.5) needs dozens of consistent steps.
  - `listen/nato_focus_m70_timing_pass.wav`; `voice/data/policy_m70.plan`, `tract/data/policy_m70.organs`.
- Owner on m70: "sounds like he has a fork in his mouth." Then: "timing and shape. we'll alternate. B-ALS on each, to converge. thus by alternation we will gradually move into shape. run it."
- **Alternating block ALS** (`policy_fit.py als`, `ORDER=alternate`, `START_FROM=` m70's best).
  - A timing block (duration, then timing, in turn) before each shape block (consonant shapes, each vowel's posture).
  - **Views.** Re-aligning every render, all residuals together, the linear model failed (duration: 9.68 predicted, 9.90 rendered): changed durations make the DTW path, and so the frame residuals, jump.
    - A timing block is now solved against the timing residuals alone (re-aligned each render).
    - A shape block is solved against the frame residuals and pitch, the alignment held at its start.
    - A step is kept if its view's cost falls and the whole cost rises at most 1%.
  - **First turns:** the timing view 5.80 → 5.50 (duration) → 5.15 (timing). The vowel blocks predict gains under 0.01, not realized: undone.
- **Alternating ALS, 3 sweeps (MEASURED).** Policy cost (seeds 73021 / 51117): start 9.745 / 9.885; after sweeps 9.649 / 10.085, 9.618 / 10.329, 9.606 / 10.299.
  - The fitted seed improved and the held-out seed worsened: the shape steps learned one seed's chance. By the two-seed mean, m70 stays best.
  - **Timing converged:** duration residual 7.9 → 5.2 (words about 28% → 23% too long); rate 1.00 → 1.16.
  - **Durations:** /uw/ 1.30, /ah/ 1.19, /w/ 1.18, /ow/ 1.16; /l/ /f/ /z/ 0.89–0.90.
  - The spectral shape barely moved (12.66 → 12.59). The consonant-shape block kept one step in three sweeps.
  - **Word lengths (s; NATO / m70 / m71):** Charlie 0.45 / 0.68 / 0.61; Golf 0.39 / 0.71 / 0.64; Juliett 0.61 / 0.89 / 0.81; Uniform 0.59 / 0.97 / 0.89; X-ray 0.51 / 0.69 / 0.63; Zulu 0.46 / 0.73 / 0.66; Lima 0.44 / 0.51 / 0.47; Sierra 0.53 / 0.50 / 0.46.
  - `listen/nato_focus_m71_alternating_als.wav` (the last state, not the two-seed best); `voice/data/policy_m71.plan`, `tract/data/policy_m71.organs`.
- Owner on m71: "the stuttering did not improve. lets work on charlie."
- **Charlie's stutter: two /tS/.** Envelope (m71): noise 0.680–0.710 s, silence, the real hiss 0.760–0.810, a crack at 0.815, a puff, voice from 0.845.
  - TRACE: the breath rises from 0.63 s, but the tip seals only at 0.73. It is 0.35 → 0.02 cm² open while 100–190 cm³/s hiss through: the first "ch".
  - **Cause:** a phrase's first consonant aimed its tip 120 ms early (close_time), but its jaw (and lips, body) only `shape.lead` (40 ms) early. The tip rides on the jaw and could not seal until the jaw rose.
  - **Fixed (planner.cpp, gesture_organs, `shape_time`):** a phrase-initial consonant's whole shape starts with its closure. **MEASURED:** before the release Charlie is silent (rms ≤ 0.001, was 0.01–0.02).
- **The crack at the hiss-to-vowel opening** (TRACE 0.5 ms): the pressure behind the channel fell 125 Pa within 1 ms as the tip left the palate, the folds still apart (0.20 cm²), closing in 20 ms.
  - **Fixed:** after a voiceless affricate the folds close over 80 ms (`glottis.affricate_adduction_move`), and its hiss channel and groove open over 0.1 s (`articulation.affricate_release_move`). Both are new tree leaves.
  - **MEASURED (three seeds; crack = the loudest 5 ms at the opening):** against the hiss 8.6 → 4.1, against the vowel 1.23 → 0.65.
  - Slower still (0.12 / 0.15): no further gain. Earlier fold closing (adduction_lead 0.03) or a wider channel (hiss press −0.2): within the seeds' spread.
  - (A first test, editing the plan with sed, showed 0.9× the hiss. m71's plan file had no lines for the new leaves, and the global tip_release_move it used also moved the release earlier. Re-tested on a full tree.)
- `listen/nato_focus_m72_charlie.wav`: the ten words on m71's plan and organs with both fixes. Lima's and Whiskey's first 150 ms dropped from 0.39 and 0.45 of the word's peak to 0.07 and 0.04 (their phrase-initial /l/ /w/ now in place before the air). Not yet heard.
- Owner on m72: "ooniform."
  - **LPC F2 of the /y/:** ours 1550 Hz, falling to /uw/'s 1270–1340. The NATO reader's "you": F2 2300–2700 throughout (a fronted /uw/).
  - **Organs during the /y/** (PROFILE): the body at /iy/'s (5.10, 7.57), but /uw/'s jaw (5.9° against /iy/'s 1.0; the body rides on it) and full rounding (aperture 0.18, protrusion 0.92). The /y/ owned only body and tip.
  - **Fixed (planner.cpp):** /y/ owns jaw, mouth and spread too. It takes /iy/'s jaw, and lips that keep `y_round_share` (0.3, new leaf) of the next vowel's shape; it hands them back at its end.
  - **MEASURED:** Uniform's /y/ F2 2150–2200 Hz for about 50 ms, gliding to 1300 by 0.79 s. Yankee's /y/ 2170–2340.
- `listen/nato_focus_m73_uniform.wav`.
- Owner on m73: "they all sound like the nato-words, but still more like a old tts. make four principled experimental changes to the vocal tube and read 'there is always money in the banana stand' with all four candidates."
- **Four tube experiments** (body tree leaves; defaults off, **MEASURED** bit-identical with them off):
  - **A. Boundary-layer loss by cross-section** (`tract.hop_loss_by_area`; Tube::set_smoothing_by_area). Friction and heat conduction at the walls go per length as perimeter over area (2√(π/A) for a round duct). Each link's hop k = k·√(3 cm² / A), at most 0.6: narrow passages lossier, wide ones less. Until now one k for every link.
  - **B. Cheek pockets** (`branches.buccal_scale`): the buccal vestibules between the teeth and cheeks, both sides as one closed side branch. 3 cm deep (canines to molars), 0.8 cm², port 0.4 cm², soft wall; it moves with the teeth (branch group 3, VoiceBody::set_teeth).
  - **C. Velopharyngeal leak** (`branches.velum_leak`, 0.03 cm²): the velum's least opening while the mouth is open. It fades as a constriction narrows from 0.35 to 0.05 cm² (the velopharynx seals under an obstruent's oral pressure).
  - **D. Head radiation** (`tract.radiation_mass_scale` 0.72, `radiation_resistance_scale` 1.04): the lip load of an unflanged end (end correction 0.61a against the baffled 0.85a; low-frequency radiation resistance about half), where a head stands between a baffle and none. Until now a piston in an infinite baffle.
- **"There is always money in the banana stand."** (m71 plan and organs with today's planner; body m64 plus each change alone; `listen/m74_tube/banana_*.wav`):

  | Variant | LTAS 0.5–1.5 / 1.5–3 / 3–5 / 5–7.5 kHz (dB re 80–500 Hz) | Band HNR (dB) | Clicks/s |
  |---|---|---|---|
  | base | −1.8 / −6.5 / −7.1 / −9.2 | 21.0 / 11.3 / 7.1 / −0.8 | 7.3 |
  | A boundary | −1.9 / −6.4 / −7.1 / −9.2 | 21.1 / 10.9 / 6.8 / −0.8 | 8.0 |
  | B cheeks | −1.5 / −7.2 / −9.4 / −9.2 | 21.9 / 11.1 / 7.3 / −0.6 | 8.0 |
  | C velum leak | −2.9 / −7.0 / −7.1 / −9.2 | 21.6 / 10.0 / 6.7 / +0.3 | 5.4 |
  | D head | −1.7 / −5.7 / −4.6 / −8.4 | 20.9 / 11.1 / 6.7 / −0.1 | 8.3 |
- Owner on the tube variants: "c and d maybe slightly better? but ALL still with the problem. dig deeper into changing how it makes the vocal sound".
- **The glottal pulse** (`/tmp/es/pulse.py` on STEPS at 96 kHz, "stand" /ae/, 9 cycles). Open quotient 0.65, speed quotient 1.96, peak flow 231 cm³/s, least area 0.0002 cm².
  - **Return phase 0.06 ms**: from the flow's steepest fall to closure, six samples. A human modal pulse: about 0.2–1 ms.
  - The flow itself is smooth (H2 −17 dB re H1; 2 kHz −56, 4 kHz −65). Almost all excitation above 1 kHz comes from one sharp corner per cycle: an impulse train's character (old formant synthesis), and the 5–9 local clicks per second (people 0–1).
- **Source experiments** (body tree; `listen/m74_tube/banana_e*.wav`; return phase, least area, clicks/s, band HNR):

  | Experiment | Return phase | Least area (cm²) | Clicks/s | Band HNR (dB) |
  |---|---|---|---|---|
  | base | 0.06 ms | 0.0002 | 7.3 | 21.0 / 11.3 / 7.1 / −0.8 |
  | **E1** string mode (`fold.string_mode`, new: displacement along the fold as a string's first mode, (π/2) sin(πz)) | 0.04 ms | 0.0009 | 6.7 | |
  | **E2** posterior chink 0.03 cm² | **1.0 ms** | 0.030 | **2.9** | 15.8 / −1.0 / −8.4 / −30 |
  | E6 chink 0.01 cm² | 0.48 ms | 0.010 | 4.8 | 21.3 / 6.4 / −3.5 / −12.9 |
  | E3 contact softness 20 µm | 0.09 ms | | 8.9 | |
  | E3 contact softness 50 µm | 0.36 ms | | 7.0 | 20.9 / 7.0 / −1.0 / −14.4 (2 kHz −70: dull) |
  | E5 chink 0.02 + softness 20 µm | 0.73 ms | | 2.2 | 17.7 / −7.5 / −4.4 / −12.9 |
  | E4 chink 0.03 + velum leak + head radiation | 0.97 ms | | 3.5 | 20.0 / −3.4 / −5.5 / −30 |

  - E1: closure is abrupt because the slices reach contact fast (a 0.01 cm rest gap swept in about 0.1 ms), not because of the shape along the length.
  - E2: with a chink, flow never stops and the pulse's end rounds off over about 1 ms. But the aspiration source (∝ v³√A) hisses continuously through it: much breathier.
- Owner on the source experiments: "softer, but the same mechanical tone underneath from the vocal cords."
  - The flow's H2 is 17 dB under H1: near a sine (human pulses about 6–12 dB).
  - Every spring of the fold tube is linear (cover, coupling, support, vertical). Real fold tissue stiffens with strain; body-cover models give the springs k(1 + ηx²), η = 100 cm⁻² (Story & Titze 1995).
- **E7 strain stiffening** (`fold.strain_stiffening`, new; cover, coupling and support springs; 0 bit-identical):

  | η (cm⁻²) | Flow H2 (dB re H1) | OQ | Peak flow (cm³/s) | Clicks/s | F0 (Hz) |
  |---|---|---|---|---|---|
  | 0 | −17.4 | 0.65 | 231 | 7.3 | 107 |
  | 100 | −14.6 | 0.60 | 221 | 7.3 | 106 |
  | 300 | −10.9 | 0.56 | 200 | 6.4 | 106 |
  | 1000 | −8.0 | 0.49 | 156 | 5.8 | 109 |

  - With chink 0.01: η 300 return phase 0.59 ms, H2 −12.4, clicks 4.8/s; η 1000 return phase 1.0 ms, H2 −9.0, clicks 3.5/s.
  - `listen/m74_tube/banana_e7_eta300.wav`, `_e7_eta1000`, `_e8_eta300_chink10`, `_e9_eta1000_chink10`.
- Owner on η 1000: "a little bit more human".
  - **Pitch: stiffening raises it** (amplitude-dependent stiffness): η 1000 +56 cents, 2000 +200, 3000 +295 uncalibrated; spread 42 → 85 cents.
  - Calibrated (`/tmp/es/caleta.sh`: natural_pitch and pitch_per_tension ×1.056 / ×1.159 / ×1.233), H2 re H1 −6.3 / −4.7 / −4.8 dB; OQ 0.48 / 0.45 / 0.45; peak flow 156 / 124 / 107 cm³/s.
  - Also new: `fold.vertical_stiffening` (vertical shear springs k(1 + ηx²)) and `fold.contact_exponent` (1.5: Hertzian contact).
- Owner: "for the vocal tube experiment, there were multiple segments. one was stiff, one was elastic, and one was i forget, and slightly adjusted. we never got beyond a more organic sounding 'aaaahh' sound, but we did get something like that. the parallel pieces allowed adding body because they responded differently, they also responded to each other, and they culmunatively shaped the resonant tone itself".
- **Parallel strips along the fold** (`fold.strips`, 1–3; Cord::fold_geometry, fold_air_areas, advance_fold_tube rewritten over [strip][slice]).
  - Each strip is a column of the slices' cover, body and support. It has its share of the length (so of the mass, stiffness and air force) times its own factors:
    - front (stiff): k ×1.6, m ×1.0;
    - middle (elastic): k ×0.6, m ×0.9;
    - back (adjusted): k ×1.0, m ×1.1.
  - Neighbours are joined by longitudinal tension (`fold.strip_coupling`, 20000 dyn/cm, cover and body). Every strip feels the same slice pressures through its own open length; the glottis is their sum.
  - **MEASURED:** one strip renders the sentence bit-identical to before.
- **Strip variants** (banana sentence, pitch-calibrated; `listen/m74_tube/banana_p*.wav`):

  | Variant | Return phase | OQ | Flow H2 (dB re H1) | Clicks/s | Pitch sd |
  |---|---|---|---|---|---|
  | P1 three strips | 0.24 ms | 0.64 | −18.0 | 6.1 | 44 cents |
  | P2 + η 1000 | 0.46 ms | 0.58 | −9.0 | 5.1 | 80 |
  | P3 + η 1000, coupling 5000 | 1.0 ms | 0.75 | −12.3 | 2.9 | 77 |
  | P4 + η 1000, contrast 2.5 / 0.4 / 1.2 | 1.25 ms | 0.84 | −16.8 | 3.2 | 64 |

  - The strips close at different times, which gives the return phase a chink gave, without a standing gap.
- Owner: "not parallel. sequential. but yea. parallel in a sense." Then the tube config:
  - section 1: mass 1.0, elasticity 2, stiffness 1;
  - section 2: extremely low stiffness, very high elasticity, mass half the elasticity;
  - section 3: mass twice the elasticity, low stiffness;
  - section 4: mass 4, almost no elasticity, stiffness 3;
  - section 5: almost no mass or stiffness, high elasticity;
  - outlet: mass 1, elasticity 2, stiffness 1.
  - Also: "letting the elasticity itself wander in a few places proportionally or maybe exponentially tamped by stiffness, gave organic character".
- **Per-slice pieces** (`fold.slice{1..6}_{mass,elasticity,support}`, new). The fold tube's six slices run inlet to outlet: the sequential sections.
  - Mass scales the body (the mass piece), elasticity the cover's stiffness and its coupling (the elastic piece), support the body's spring to the muscle. Damping follows each slice's pieces at a fixed critical share.
  - **Elasticity wander** (`fold.elasticity_wander`, rms share; `elasticity_corner` 20 Hz; `elasticity_tamping`): per slice, through one pole, times exp(−tamping (support − 1)).
  - **MEASURED:** all 1 and no wander, bit-identical.
- **The owner's config as factors** (m/e/s): 1/2/1, 2/4/0.1, 2/1/0.3, 4/0.1/3, 0.1/3/0.1, 1/2/1. Banana sentence, pitch-calibrated (`listen/m74_tube/banana_t*.wav`):

  | Variant | OQ | Return phase | Flow re H1: H2 / 1 / 2 / 4 / 6 kHz (dB) | Clicks/s | RAP | Pitch sd |
  |---|---|---|---|---|---|---|
  | base | 0.65 | 0.06 ms | −17.4 / −40.8 / −55.7 / −65.3 / −82.7 | 7.3 | 0.24% | 42 cents |
  | T1 sections | 0.67 | 0.04 ms | −23.9 / −37.7 / −47.4 / −57.4 / −73.1 | 4.8 | 0.30% | 43 |
  | T2 + wander 0.05 | 0.65 | 0.04 ms | −24.3 / −38.3 / −48.6 / −61.7 / −72.9 | 7.0 | 0.91% | 50 |
  | T3 + wander 0.05, η 1000 | 0.57 | 0.06 ms | −10.5 / −35.5 / −48.8 / −60.4 / −70.6 | 7.3 | 0.79% | 84 |
  | T4 + wander 0.15 | 0.63 | 0.04 ms | −20.8 / −38.1 / −51.2 / −61.7 / −74.6 | 9.6 | 2.59% | 97 |

  - The sections put 8–10 dB more into the flow pulse itself at 2–6 kHz: structure inside each pulse, not only its corner. The wander roughens the period (RAP 0.9–2.6%; the owner 0.26%, NATO 0.29%).
- Owner: "banana_t3_sections_wander_eta sounds a bit more organic. he has an organic sense of body, at least. it warbles on top of that, but that's mouth warble."
  - T3 saved as `data/voice_t3.params`: m64, the six sections, elasticity wander 0.05, η 1000, pitch map ×1.015.
- **Wobble in voiced stretches** (`/tmp/es/mwobble.py`: rms about a 120 ms average):

  | Recording | Level | Spectral balance (1–3 kHz re total) | Pitch, audio tracker |
  |---|---|---|---|
  | owner (user_three) | 1.13 dB | 1.86 dB | 16 cents |
  | bdl a0001–3 | 0.61–1.22 | 1.73–2.52 | 11–94 |
  | base | 1.53 | 2.57 | 77 |
  | T3 | 1.63 | 1.56 | 253 |

  - The audio tracker's pitch jumps octaves. From the engine's CYCLES, cycle-level pitch wobble is base 15, T1 14, T3 20 cents; cycles jumping >5%: 6, 11, 11%; no period doubling (sign alternation 0.36–0.53).
  - Motor wander off, wall variation off, both: level 1.49 / 1.65 / 1.75 dB (T3 1.63). On one sentence (2–4 stretches) not separable. `listen/m74_tube/banana_t3_w{1,2,3}_*.wav` for the owner's ear.
- Owner on the motor and wall variants: "they all sound the same. lets remove the tongue entirely and see how he sounds".
  - **Frozen** (`banana_t3_tongue_frozen.wav`): every tongue command after t = 0 removed (body, apex, groove, lateral); the tongue holds the opening schwa (AX: body 3.71, 5.96; tip 6.64, 6.85). Jaw, lips, velum, larynx, glottis and breath still move.
  - **Lowered** (`banana_t3_tongue_lowered.wav`): the same, with the body 1 cm and the tip 0.5 cm lower (PROFILE: ty 4.93, tipy 6.34).
  - Wobble (level / balance): T3 1.63 / 1.56 dB; frozen 1.48 / 1.14; lowered 1.27 / 1.22.
- Owner: "lowered is best, and no warble. almost human. frozen, a bit of clicky/wet sounds, but no warbling".
  - **The warble is the tongue's movement.** The voice without it is close: the T3 body.
- **The area function as articulated** (`speak AREALOG=path`, new: every millisecond, every vocal node's rest area to the lip opening; `/tmp/es/arearough.py`).
  - Fast fluctuation (rms log area about a 10 ms average, oral nodes): moving 0.117, frozen 0.045, lowered 0.033. Worst at 12.6–14.0 cm (the blade under the ridge).
  - But that is the consonants' closures. Inside the vowels the areas glide smoothly ("stand" /ae/: 12.9 cm 1.50 → 2.27 cm² over 40 ms, then slowly back): no geometric jitter.
  - The tongue never settles: moves up to `vowel_reach` 0.35 s, arriving `vowel_arrival` 0.88 into the vowel (from the timing fit). The formants slide through every vowel.
- **Timing variants on T3** (`listen/m74_tube/banana_t3_v*.wav`):
  - V1: vowel_arrival 0.5;
  - V2: also vowel_reach 0.12;
  - V3: also vowel_hold 0.8.
  - The wobble measure has too few voiced stretches here (1–2); for the owner's ear.
- Owner: "banana_t3_v3_arrive50_fast_hold subjectively maybe?"
- **What couples tongue movement into the sound.** Frozen and lowered both lost the warble, so it is the movement, not the posture. On V3 (`listen/m74_tube/banana_v3_x*.wav`):
  - X1: the ear off (`ear.pitch_feedback` 0): moving formants shift the heard period it corrects each cycle;
  - X2: walls ×3 stiffer (`walls.stiffness_scale` 3);
  - X3: no pharyngeal pull in voiced closures (`walls.closure_pull` 0);
  - X4: no motor wander (`variation.wander_gain` 0).
- Owner: "i think we need to re-articulate hte tongue itself maybe.. noear is better".
  - The ear off raised the pitch 260 cents: T3's map had been calibrated with the cycle-by-cycle ear running.
  - **A listening ear** (`ear.latency`, `ear.smoothing`, `ear.rate`, new): the log error smoothed, delayed, integrated; latency 0 keeps the old ear, bit-identical.
- **Calibrated on the V3 banana score** (`/tmp/es/earcal.sh`):

  | Ear | Pitch error (median) | Spread | Pitch map |
  |---|---|---|---|
  | per cycle | +3 cents | 64 | ×1.021 |
  | none | +4 | **51** | ×1.151 |
  | listening (latency 0.12 s, rate 3/s) | +22 | 93 (it hunts) | ×1.049 |

  - Without the ear the pitch is steadier than with either.
  - `data/voice_t3_noear.params`; `listen/m74_tube/banana_v3_no_ear_calibrated.wav`.
- Owner: "build it." (the physical tongue)
- **The physical tongue** (`include/tract/tongue.hpp`, `src/tongue.cpp`; `tongue.physical` 1; `speak TONGUELOG=path`).
  - **Body.** 24 sections of the tongue's surface from the root's base (held beside the epiglottis) over the dorsum to the tip, in the jaw's frame. Each has mass, elasticity (stretch and bending to its neighbours) and support (to the rest shape). Support and bending fall toward the tip (`tip_pliancy` 0.2).
  - **Constraints.** A hydrostat (the midsagittal area, surface closed through the floor, held at rest); contact (the hull stops it, `contact_depth` 0.04 cm past).
  - **Muscles.** Genioglossus posterior and anterior from the genial tubercle; styloglossus from the styloid; hyoglossus from the hyoid. Superior and inferior longitudinal curl the blade along its own normal (×3); the floor's lift.
  - **Control.** Activations follow their drive over 15 ms. The drive is feedback: the planner's raw targets for body and tip as a target shape, given the tongue's own volume, the error less 0.08 s of velocity, projected on each muscle's pull, plus an integral.
  - **Geometry.** `MidsagittalTract` takes an explicit tongue contour (`Articulators::contour`); `tongue_outline()` samples the geometric tongue as targets.
  - **MEASURED:** off, the sentence is bit-identical.
- **Findings on the way:**
  - **The geometric tongue does not conserve volume:** its area ran 16.4–30.0 cm² over the sentence (rest 18.65). Targets are now shifted along their normals to the rest area; median shape error 0.59 → 0.52 cm.
  - **The error does not fall with gain** (×1 to ×30): seven muscles cannot draw every circle-and-tip shape. The vowels are to be refit in the new body.
  - **Closures** (least oral area in the phone; geometric → physical, f4k):
    - DH (phrase-initial) 0 → 1.40 cm²; DH 0 → 0.008;
    - N 0 → 0.01–0.18; T 0.001 → 0.078; D 0 → 0;
    - S 0.055 → 0.060; Z 0.055 → 0.11–0.17;
    - L 0.09 → 0.35; R 0.45 → 0.10.
    - With uniform stiffness and a fixed up-and-back longitudinal pull, no /t d n/ sealed (0.16–0.97).
- `listen/m74_tube/banana_physical_tongue_1.wav` (T3, no ear, V3 timing; `data/voice_physical_tongue_1.params`).
- Owner on physical tongue 1: "it's golden ... This thing is garbage. It sounds terrible ... like he's got a very serious speech impediment, like his tongue has been mangled. But here's the most magical, wonderful thing about it. He sounds like he has a real physical tongue."
- **Control and reach of the physical tongue** (`tools/tongue_probe`: the tongue alone, holding rest then a target):
  - **Feedforward.** Each muscle's influence is measured at reset (alone, activation 0.1, scaled; settled 0.6 s, no contact). For each target, bounded least squares by coordinate descent gives the activations (0..1); feedback works on top.
    - The first solver (step 1/trace) barely moved: activations ≤ 0.33.
    - Influences measured at full activation superposed badly: held 0.9 cm from a target predicted to 0.2.
  - **Reach: 20 muscles.** Genioglossus 3 bands, styloglossus 2, hyoglossus 2, longitudinals 2, floor lift, transversus and verticalis at 5 places each, along the normal.
    - With 7 the least-squares hold of sentence targets left about 0.8 cm; with 20, median 0.13–0.20 cm.
  - **Probe (an /aa/-like target):** feedforward alone 0.46 → 0.09 cm in 260 ms; with feedback (gain 3) 0.18 at 80 ms, 0.10 at 140, 0.06 at 260.
  - **Speed.** Support, elasticity, bending and muscle force ×4 (support about 6 Hz): sentence error median 0.46 → 0.40 cm (×9: 0.35, closures worse).
  - **Look-ahead** (`tongue.lead`): it aims at the next body and tip targets due within the lead. Aiming at the last one skipped consonants for the following vowel.
- **Closures, lead 0.03** (least oral area in the phone, cm²):

  | Consonant | Geometric | Physical |
  |---|---|---|
  | N | 0 | 0, 0, 0, 0.06, 0.28 |
  | D | 0 | 0 |
  | DH (second) | 0 | 0.001 |
  | DH (first) | 0 | 0.52 |
  | T | 0.001 | 0.018 |
  | S | 0.055 | 0.044 |
  | Z | 0.055 | 0.32–0.72 |
  | L | 0.09 | 0.96 |
  | R | 0.45 | 0.82 |

  - Shape error median 0.40 cm.
- `listen/m74_tube/banana_physical_tongue_2.wav`; `data/voice_physical_tongue_2.params`.
- Owner on physical tongue 2: "He does still sound like his tongue has been injured. But we can work".
- **Consonants of the physical tongue** (`tools/tongue_probe`: hold a consonant target alone, the narrowest oral opening over time):
  - **Speed.** /l/ and /t/ closed only after about 180 ms (a consonant lasts 50–80 ms).
    - Lighter or tapered sections (`section_mass`, `tip_mass`) changed nothing.
    - The controller's velocity term (`damping_gain` 0.08 s: 1.6 cm of "error" at 20 cm/s) braked every move. At 0.02: /l/ 0.020 cm² at 40 ms, /t/ 0.017 and /z/ 0.058 at 80 ms.
    - Higher gains overshot and sprang open again.
  - **Contact at corners.** The tongue's hull now ends at the upper incisors' edge, and a vertex's side is judged by the bisecting normal. A point below and behind the edge had read as inside the tooth.
  - **Protrusion.** Nothing could push the tip forward: the genioglossus pulls toward the chin, the longitudinals along the normal, the support back, the hydrostat along normals. The tip stopped at x 7.34, 0.4 cm short of the incisors for /dh/, whatever the stretch stiffness (`tip_stretch`, new).
    - New 21st actuator: **elongation** (transversus and verticalis squeezing the hydrostat lengthen it), on the blade along the surface toward the tip, ×2. /dh/ in the probe 1.20 → 0.34 cm², now at the teeth.
  - The geometric /dh/ target lies inside the upper incisors (x 7.78, y 8.51; the back face at y ≈ 8.1 there). The old tongue had no contact and "closed" through the tooth.
- **Sentence, physical tongue 3** (damping_gain 0.02, no lead; `data/voice_physical_tongue_3.params`, `listen/m74_tube/banana_physical_tongue_3.wav`):
  - shape error median 0.25 cm (2: 0.29–0.40);
  - least oral area (cm²): N 0 / 0.0001 / 0 / 0.0002 / 0; DH 0 (second), 0.40 (first); T 0.005; D 0; S 0.027; L 0.12; Z 0.10 and 0.87 (the first, "is"); R 0.51 (geometric 0.45).
- Owner on physical tongue 3: "Massive regression. He's now warbling pretty badly,".
- **A warble measure that follows the owner** (`analysis/tongue_ringing.py`): from AREALOG, the oral nodes' log areas inside vowels, the rms of their 5–30 Hz part, median over vowels.

  | Version | Ringing | Owner |
  |---|---|---|
  | geometric tongue | 0.111 | warble |
  | physical tongue 2 | 0.041 | "injured", no warble |
  | physical tongue 3 | 0.116 | "warbling pretty badly" |

  - **The cause:** the controller's velocity braking cut to 0.02 s (for speed), with the 30 ms look-ahead dropped. The tongue overshot and rang around its targets (probe: the tip 0.5 cm past at 40 ms).

  | Variant (protrusion and contact fixes kept) | Ringing |
  |---|---|
  | braking 0.08, no lead | 0.058 |
  | braking 0.04 | 0.064 |
  | tissue viscosity (support_loss 2, elastic_loss 1), braking 0.02 | 0.095 |
  | **braking 0.08, lead 0.03** | **0.041** |
  | braking 0.06, lead 0.03 | 0.051 |
  | braking 0.04, lead 0.03 | 0.078 |
  | braking 0.08, lead 0.05 | 0.039, closures worse |

- **Physical tongue 4** = braking 0.08, lead 0.03, with elongation and the contact fix (`data/voice_physical_tongue_4.params`, `listen/m74_tube/banana_physical_tongue_4.wav`).
  - Closures (cm²): N 0.0001 / 0 / 0.029 / 0.0015 / 0; T 0.016; D 0; S 0.046; L 0.10 (tongue 2: 0.96); R 0.40 (0.82); Z 0.28, 1.10; DH 0.62, 0.18.
- **Rule from here:** every tongue change is checked against the ringing (≤ 0.04) before the owner hears it.
- Owner on tongue 4: "doesn't actually sound better than version 2. And we could probably refine it from this point with a bit more ML. But you be the judge of what to do next, and refine towards what you want to".
  - Judgement: the vowels (postures fit for the geometric tongue), not the last consonant closures.
- **Vowels refit through the physical tongue** (organ_vopt3, SPEAK_ENV FOLDTUBE=1 PARAMS=tongue 4, W 10, MAXFEV 260; from m71's postures; `data/organs_physical.txt`).
  - Formant distance to bdl (st), before → after: uh 7.00 → 1.44, ao 5.60 → 1.62, ah 2.98 → 1.41, iy 2.58 → 1.20, ow 2.57 → 0.18, ih 2.01 → 0.26, ax 1.96 → 0.15, uw 1.64 → 0.72, eh 1.07 → 0.15, aa 0.94 → 0.82, ae 0.88 → 0.35, ey 0.80 → 0.04.
  - Shape rms 0.31–0.48.
  - The jaw rose to 14° for aa, ow (the bound).
- **Physical tongue 5** (`listen/m74_tube/banana_physical_tongue_5_vowels.wav`): ringing median 0.046 (tongue 4: 0.041), largest vowel 0.42.
  - Closures (cm²): N 0 / 0 / 0.03 / 0.37 / 0.05; DH 0.63, 0; T 0.08; D 0.016; S 0.09; L 0.12; Z 0.42, 0.25; R 0.36.
- **Prepared for the ML refinement** (engine_net.py):
  - `ENGINE_TREE`, `ENGINE_ORGANS`, `ENGINE_PLAN` set the start;
  - group "tongue" (the physical tongue's settings, ×/÷2; not its switch);
  - the NATO residual gains each word's vowel ringing as a hinge (free to 0.035, then 1 per 0.002).
  - `data/voice_tongue4.params`: tongue 4's canonical tree.
- Owner on tongue 5: "it trades a nonphysical garble at 2.5 seconds for a click and a more physically plausible garble. thats positive."
- **ML refinement** (policy_fit.py adamw, GROUPS tongue,organs,plan, from tongue 5; /tmp/es/tml): policy 17.47 → 13.82 by step 15. MEASURED: nearly all of it is the ringing hinge (624 → 206); shape 14.17 → 14.10, duration 5.8 → 7.1. The tongue settings moved little (braking 0.08 → 0.091, muscle force −9%). It quiets the warble by slowing, not by better articulation.
- **Two faults behind the garble** (AREALOG plus the tip in TONGUELOG):
  - "the" (1.70 s): the tip moved smoothly while the narrowest area jumped 0.005 ↔ 0.15 cm² from one millisecond to the next. The cause was the fixed ray cuts: as the sharp tip slid past a ray, its crossing switched from the tongue's top to its underside.
  - "stand" (2.47–2.52 s): the tip rose toward the ridge at about 5 cm/s and never sealed /n d/ (0.016 cm²).
- **Why the tongue is slow** (tongue_probe, /æ/ → /n/ target):
  - It sealed at 160 ms; at braking 0.04 it seals at 82 ms; at 0.02 the tip overshot past the ridge.
  - With feedback off, the influence-based feedforward barely moved it (tip 1.25 cm off). The feedback did all the work, and the braking capped its speed at error / 0.08 s.
- **Inverse statics feedforward** (`tongue.statics 1`): at the target shape, the activations (0..1) plus the hydrostat's pressure (free sign) that balance the tissue's passive forces there. Bounded least squares by coordinate descent, warm-started, solved every control tick.
  - Feedforward alone now seals in 96 ms with the tip on its aim (6.75, 10.25).
  - Off is bit-identical.
- **Continuous passage width** (physical tongue only): each cut's width is min(ray width, the hull point's distance to the nearest organ segment); the midpoint follows.
  - Large flicker events (second difference of log narrowest area > 0.5) fell 46 → 20.
  - The geometric tongue is bit-identical.
- **Soft mucosal contact** (`tongue.contact_stiffness`, dyn/cm^1.5; `contact_loss` 600 dyn s/cm): Hertzian p^1.5 push, lossy while pressed, with a hard stop at 4× contact_depth.
  - It did not change the remaining flicker (C 40, G 39). That flicker is the tip sliding forward and down the ridge, sealed, as /n/ releases into /iy/ (about 80 cm/s; the contact crosses tube sections).
  - Off is bit-identical.
- Variants on the banana sentence (seed 73021; ringing median/max; closures as sealed ms / smallest cm²):

  | | statics | braking | soft | ringing | vowels closed | final N D |
  |---|---|---|---|---|---|---|
  | tongue 5 | 0 | 0.08 | – | 0.046 / 0.42 | IY, IH of "money in" (0.00) | 0 / 0.016 |
  | 6 = width fix | 0 | 0.08 | – | 0.050 / 0.26 | IY, IH (0.00) | 0 / 0.013 |
  | C | 1 | 0.04 | – | 0.057 / 0.19 | none | 73 / 28 ms |
  | G | 1 | 0.04 | 600000 | 0.056 / 0.25 | IH of "is" (0.00, two bursts) | 71 / 26 ms |
  | **7 = H** | 1 | 0.06 | 600000 | 0.056 / **0.094** | none (≥ 0.16) | 36 ms / 0.007 |

  - The statics feedforward reaches its targets, exposing targets the volume shift inflates into the palate (G's "is").
  - Tongue 5 already sealed the tract in the vowels of "money in".
- Sent: `listen/m74_tube/banana_physical_tongue_6_smooth_airway.wav`, `banana_physical_tongue_7_statics_soft.wav` (`data/voice_physical_tongue_7.params`).
- Owner on 6 and 7: "6 doesnt click, but is slightly muffled where we clicked before and is more tongue-tied and 7 does not sound better. none of these do. lets walk it back. we're trying too hard."
  - Walked back: statics, soft contact and the smooth airway are switches, all off by default (`tongue.smooth_airway`); tongue 5 and the geometric tongue reproduce bit-identically.
- Owner: "critically think about the mouth as a musical instrument, identify the role of the tongue in each sound ... and enable or disable or transition between states that we believe represent these roles."
- **The tongue's roles** (`include/tract/tongue_roles.hpp`, `src/tongue_roles.cpp`; speak reads the score's phone comments):
  - **shaper**, the resonator (vowels, glides, /r/, /h/, rest): it sets the cavities, slowly and stably;
  - **valve**, the occluder (/t d n l th dh/ and the flap with blade and tip; /k g ng/ with the dorsum): it strikes, seals and releases;
  - **nozzle**, the jet's aperture (/s z sh zh ch jh/, blade and tip): it holds a still, precise channel.
  - Each piece (root 1–5, dorsum 6–14, blade 15–20, tip 21–23) takes the role of the phone sounding. A role comes into force over role_ramp (25 ms), from role_lead (20 ms) before its phone to role_hold (10 ms) after.
- **Representation 1: pinned to the geometric outline** (spring 1e6 dyn/cm, critically damped; the hydrostat moves only the free sections). REJECTED for now:
  - R1 (valves and nozzles pinned) sealed every closure exactly as the geometric tongue does.
  - But the body wobbled back and forth inside the vowels (node 20 in /ao/, every 10 ms: 1.75 2.16 1.76 1.98 2.25 1.55 2.08). At each handover the pinned shape and the body's own disagree, and the body springs back.
- **The ringing measure is confounded.** Turning off the motors' wander left geometric's 5–30 Hz value unchanged (0.113 → 0.116). Inside the vowels the geometric tongue sweeps large and monotonically (/ao/ node 20: 0.84 → 0.19 → 1.58), with one reversal per node.
  - A 50 ms vowel between consonants is itself a 10 Hz movement. The band partly scores how much the tongue moves, so the ≤ 0.04 gate favoured a sluggish tongue.
  - Added `/tmp/es/bands.py`: slow (5–30 Hz), mid (30–100 Hz) and jitter (> 100 Hz) motion of the oral log areas in vowels. Wobble shows in mid and jitter.
- **Representation 2: impedance** (co-contraction; the equilibrium-point view). The role sets the piece's stiffness toward its posture (dyn/cm per section) and its damping ratio.
  - Anchored on the feedback's target (role_anchor 0), the stiff valves jittered (40000: jitter 0.0155). The target steps as planner events enter the look-ahead.
  - Anchored on the motors' smooth posture (role_anchor 1: the geometric outline as the motors hold it, shifted to the tongue's own area):

  | | valve / nozzle stiffness | shaper | slow | mid | jitter | sealed ms (N N N N N, final D) | /s t/ |
  |---|---|---|---|---|---|---|---|
  | tongue 5 | – | – | 0.046 | 0.016 | 0.0065 | 44 84 0 0 0, 0 | 0.08 |
  | geometric | – | – | 0.113 | 0.017 | 0.0053 | 55 50 67 51 90, 58 | sealed |
  | **A1** | 16000 | free | 0.082 | 0.011 | 0.0035 | 40 48 59 50 55, 13 | 4, 8 ms |
  | A2 | 40000 | free | 0.087 | 0.011 | 0.0035 | 42 43 41 45 15, 13 | 0.013 |
  | A4 | 80000 | free | 0.092 | 0.015 | 0.0059 | 36 7 25 34 6, 13 | 12, 13 ms |
  | **A5** | 40000 | 10000, ζ 1.5 | 0.081 | 0.010 | 0.0028 | 39 44 39 38 12, 14 | 0.011 |

  - A1 and A5 open the vowels of "money in" (/iy/ 0.04–0.07 and /ih/ 0.56–0.64; tongue 5 sealed both).
  - The sentence-initial /dh/ and the /z/ of "is" stay open.
- Sent: `listen/m74_tube/banana_roles_1_valves.wav` (A1, `data/voice_roles_A1.params`) and `banana_roles_2_valves_soft_vowels.wav` (A5, `data/voice_roles_A5.params`).
- Owner on roles 1 and 2: "still sounds like he got stung by a bee."
- **The swollen tongue, measured:**
  - The hydrostat held the midline area at rest (18.5 cm²) while the planner's shapes asked 14.4–29.6 cm² (TONGUELOG now logs the planned area). "the banana", AH and AE wanted 14.5–17.7 and were inflated.
  - The schwas left 8.7–11.7 cm³ of oral cavity against the geometric tongue's 13–15.
  - Formants (LPC at mid-vowel, `/tmp/es/fmt.py`, against the fitted targets): back and central vowels fronted. AO F2 1203 (target 924), AH 1383 (1158), schwas 1650–1740 (1452). Front vowels near target.
  - All-pinned R5 gave AH 1054 (geometric 1033): the outline path is right, and the fixed area is what fronts AH.
  - Stiffening the resonator alone (S1–S4, up to 60000) left AO at 1142 and AH at 1358, so speed was not the cause.
- **Volume kept by width** (`tongue.lateral_spread`, `lateral_narrow`, `lateral_time`; `lateral_from_posture`): a tongue keeps its volume, not its midline area. Flattened it spreads toward the teeth; bunched it narrows. The midline area held follows the plan's within rest × (1 − spread) .. (1 + narrow) over lateral_time. Off is bit-identical.
- Results (tongue 5 score, seed 73021):

  | | F1 / F2 median error (st) | AO, AH F2 | mid / jitter | sealed N N N N N, D |
  |---|---|---|---|---|
  | geometric | 1.36 / 1.11 | –, 1033 | 0.017 / 0.0053 | 55 50 67 51 90, 58 |
  | tongue 5 | 2.04 / 2.36 | 1203, 1383 | 0.016 / 0.0065 | 44 84 0 0 0, 0 |
  | roles A1 | 1.90 / 1.54 | 1214, 1377 | 0.011 / 0.0035 | 40 48 59 50 55, 13 |
  | **L3** (roles 40000/30000, spread 0.3, narrow 0.6) | 1.09 / 1.29 | 1025, 1268 | 0.024 / 0.0077 | 39 35 47 45 75, 40 |
  | **L3m** (L3, area from the motors' posture) | 1.54 / 0.96 | 967, 1230 | 0.016 / 0.0057 | 51 18 17 14 9, 8 |

  - L3 has a click at 1.889 s (the /n/ of "banana"; transient ratio 42); L3m has none.
- Sent: `listen/m74_tube/banana_roles_3_unswollen.wav` (L3, `data/voice_roles_L3.params`) and `banana_roles_4_unswollen_smooth.wav` (L3m: the same file plus `lateral_from_posture 1`).
- Owner on roles 3 and 4: "he still has a bee sting, but it's coming down in inflation."
- What was left (L3m against the geometric tongue, oral cavity split pharynx / mid 10.5–14 cm / front):
  - Mid-mouth volume was short in every non-high vowel: AO 10.2/13.7, AH 11.3/14.0, schwas 8.4–12.3/11.9–14.6 cm³.
  - IH and IY were not raised enough (F2 1770–2070 against 2073–2177).
  - The body's range was compressed.
- A pliable surface (P1–P5: stretch 80000 → 8000–20000, bending 8000 → 800–2000, support ÷4) changed nothing (F2 0.96–1.01 st): the shape follows its anchor.
- **The inflator was the muscle feedback aiming at the look-ahead target** (the next consonant's raised tongue):

  | | F1 / F2 error (st) | AH, IY, IH F2 | slow / mid / jitter | sealed N N N N N, D |
  |---|---|---|---|---|
  | L3m | 1.54 / 0.96 | 1230, 2070, 1980 | 0.075 / 0.016 / 0.0057 | 51 18 17 14 9, 8 |
  | lead 0 | 1.80 / 0.99 | 1142, 2046, 2020 | 0.080 / 0.016 / 0.0071 | 31 18 26 11 16, 0 |
  | **P7**, gain 1 | 1.40 / **0.81** | 1220, 2104, 1997 | 0.076 / 0.013 / 0.0051 | 48 17 23 23 10, 11 |
  | **P8**, feedback off | **0.95** / 1.00 | 1067, 2225, 2129 | 0.097 / 0.018 / 0.0058 | 52 47 42 47 86, 57 |
  | geometric | 1.36 / 1.11 | 1033, 2205, 2122 | 0.113 / 0.017 / 0.0053 | 55 50 67 51 90, 58 |

  - P8: the body moves by its roles' stiffness toward the motors' posture (with passive tissue, volume by width, contact). Its vowels are the nearest to target, its closures like the geometric tongue's, and its slow motion approaches the geometric tongue's. Its 2.566 s transient is the final /d/'s release burst.
- Sent: `listen/m74_tube/banana_roles_5_less_anticipation.wav` (P7, `data/voice_roles_P7.params`) and `banana_roles_6_no_feedback.wav` (P8, `data/voice_roles_P8.params`).
- Owner on roles 5 and 6: "5 is less swollen, but more of the robot gets through."
  - Reading: the feedback carries the organic part; its look-ahead aim carried the swelling.
- **Feedback aimed by role** (`tongue.feedback_aim`: 0 look-ahead, 1 the motors' posture, 2 by role, from the role schedule's shaper share per piece): full feedback (gain 3), but vowel-shaping pieces correct toward the posture now while valves and nozzles aim ahead.

  | | F1 / F2 (st) | AH, IY F2 | slow / mid / jitter | sealed N N N N N, D |
  |---|---|---|---|---|
  | L3m (aim ahead) | 1.54 / 0.96 | 1230, 2070 | 0.075 / 0.016 / 0.0057 | 51 18 17 14 9, 8 |
  | P7 (gain 1) | 1.40 / 0.81 | 1220, 2104 | 0.076 / 0.013 / 0.0051 | 48 17 23 23 10, 11 |
  | **G1** (aim 2) | 1.33 / 0.87 | 1192, 1975 | 0.075 / 0.016 / 0.0058 | 32 23 19 13 12, 7 |
  | G2 (aim 1) | 2.50 / 1.28 | 1198, 1989 | 0.079 / 0.016 / 0.0077 | 38 21 32 15 22, 1 |
  | **G3** (aim 2, no resonator stiffness) | 1.66 / 0.87 | 1182, 1951 | 0.078 / 0.012 / 0.0045 | 32 20 23 10 15, 3 |
  | G4 (aim 2, resonator 10000) | 1.48 / 0.84 | 1188, 1954 | 0.076 / 0.013 / 0.0048 | 33 21 20 11 15, 3 |

- Sent: `listen/m74_tube/banana_roles_7_aim_by_role.wav` (G1) and `banana_roles_8_aim_by_role_muscle_vowels.wav` (G3); `data/voice_roles_G1.params`, `voice_roles_G3.params`.
- Owner: "this is endless tweaking but still a nonphysical mouth that doesnt actually work. can you turn this into a video i can see where a physical 2d mouth articulates equivalently to the model while the model talks?"
- **The mouth video:**
  - `speak` GEOLOG=path (GEOEVERY ms, default 20) logs each frame's midsagittal geometry (hull, organ outline, the tongue's section count) and the vocal node areas.
  - `analysis/mouth_video.py` draws it beside the tube areas and the waveform with phone labels, and muxes the audio (ffmpeg; SLOW stretches the audio, pitch kept).
  - Sent `listen/video/mouth_realtime.mp4` and `mouth_slow4.mp4` (roles 8 = G3, banana sentence, seed 73021; the audio is identical to the render's).
- OBSERVED in frames:
  - At /ao/ (0.98 s) the tongue fills most of the mouth and its root narrows the pharynx to about 0.2 cm² at 4–6 cm.
  - At /n/ (1.50 s) the tip presses past the hull into the incisors' corner.
  - The jaw, lower teeth and lips are coarse polygons.
- Owner on the video: "that's not a tongue. that's a lump of flesh that is mostly attached to the base. thats our first problem ... we can add guards that epistemically protect the space for airflow we know should be present and guard them from contracting too far. imagine- like an inverse tongue."
- **The airway guard** (`tongue.guard`, `guard_keep`, `guard_stiffness`, `guard_loss`; tongue.cpp build_guard, guard_forces, measure_intrusion):
  - The planned shape (the geometric outline as the motors hold it) says where air should be. At each section, guard_keep (0.7) of the planned gap to the hull is protected.
  - A section pressing into it is pushed back (200000 dyn/cm, loss 400 on its intruding speed). Where the plan closes the gap, the guard lies on the hull. With the roles it yields on valves and nozzles (weighted by the resonator's share).
  - guard_stiffness 0 only measures. Off is bit-identical.
  - GEOLOG gains the planned surface, the guard surface and the intrusions. mouth_video.py draws the guard (blue), the plan (dashed) and intruding sections (yellow), and prints intrusion statistics. `/tmp/es/intrusion.py` lists them by phone.
- Intrusion into the guarded airway (5 ms frames; deepest per frame):

  | | median | 90th | largest | over 0.1 cm |
  |---|---|---|---|---|
  | tongue 5, measured only | 0 | 0.215 | **1.60** | IY 1.29, IH 1.60, AX 0.87, B 0.56, AE 0.57, R 0.39, N 0.37, W 0.30, L 0.23 |
  | roles 8, measured only | 0 | 0.060 | 0.32 | R 0.32, AO 0.31, L 0.25, IH 0.20 … |
  | roles 8 + guard | 0 | 0.006 | 0.066 | none |
  | tongue 5 + guard (no roles) | 0 | 0.001 | 0.145 | AX, AE; but its closures fail (N 0.40, 0.11, final 0.15): without roles the guard does not yield for valves |

- OBSERVED in the videos:
  - Tongue 5 at /iy/ of "money" (1.57 s): the tip stays stuck up on the ridge after /n/ while the plan has it down behind the lower teeth (intrusion 1.29 cm). The vowels of "money in" are near-silent in the waveform.
  - Roles 8 + guard at the same moment: the vowel is open, but the tip curls back under itself (the last three sections hook down and back).
  - Roles 8 + guard sounds nearly unchanged (formants F2 0.89 against 0.87 st): roles 8 mostly kept to the domain already.
- Sent `listen/video/guard_tongue5_detect_slow4.mp4` and `guard_roles8_enforced_slow4.mp4` (`data/voice_roles_G3_guard.params`).
- Owner on the guard video: "does not protect. tongue still penetrates top of mouth."
  - The guard yields on valves (the /n/), so it never stood between the tip and the palate.
  - MEASURED (`/tmp/es/penetration.py`, GEOLOG field 8 = each section's depth inside the closed hull): sections never went deeper than contact_depth (0.04 cm).
  - What showed was the tip as a needle. At 1.65 s the last section stood 0.4 cm straight up off the blade (sections 24 → 25) and jabbed the ridge. The planner aims the tip about 1 cm inside the palate for /n/, and the tip was pliant (0.2 of the root's support and bending).
- **The solid hull** (`tongue.solid_hull`): ray-casting inside test on the closed hull (to the incisors' edge, closed far above). A section inside returns to its nearest hull point less contact_depth, whichever face that is, and loses its speed into the solid. Off is bit-identical.
- Tip bend during tip consonants (`/tmp/es/spike.py`: degrees between the blade, sections 20 → 22, and the last segment 22 → 23) and the deepest section inside:

  | | contact_depth | tip_pliancy | tip bend (N N N N N, D) | inside | sealed N N N N N, D |
  |---|---|---|---|---|---|
  | roles 8 + guard | 0.04 | 0.2 | 45 73 60 69 71, 68° | 0.04 | 32 20 21 10 15, 4 |
  | H2 | 0 | 0.2 | 61 80 67 54 82, 77° | 0 | 23 16 21 19 11, 4 |
  | H3 | 0 | 0.5 | 38 72 59 54 72, 71° | 0 | 24 13 15 18 12, 3 |
  | **H4** | 0 | 1.0 | 19 49 45 39 27, 26° | 0 | 25 4 9 10 64, 38 |

  - Formants unchanged (H4: F1 1.49, F2 0.83 st).
- Sent `listen/video/solid_hull_slow4.mp4` (H4, `data/voice_roles_solid_hull.params`). At /n/ (1.65 s) the blade meets the ridge at the hull with no spike.
- Still open: the tongue's front is a slab dropping from the tip to the floor (no free blade, no sublingual space).
- Owner on the solid hull: "It's a terrible tongue, and the geometry of it is not very good. We should rework it. Be more like a physical tongue. Obviously, it needs a fat tip. But his n are at least physical now."
- **The thick tongue** (`tongue.blade_body`, tip_radius, blade_thickness):
  - The underside follows the upper surface as the lower face of a thick beam: a cap of three points, a four-section underside, and a frenulum dropping onto the floor. Articulators::underside puts it in the organ outline, so the air beneath the blade is outside the organs, and the floor is drawn along its line.
  - Contact: hull contact on the cap moves the tip. The floor solid (floor plus the lower incisors' inner face, jaw frame) only squashes the blade. Pushing sections from the floor lifted a resting tip 1.3 cm per substep and the body crumpled.
  - The body's area closes from the frenulum straight to the floor's back; through the floor's front corner the outline crossed itself.
  - The profile's tip (groove and side channels) is the cap's front.
  - K1 (radius 0.4, thickness 1.5): N sealed 61 52 83 64 97 ms, D 58; nothing inside the palate; F1 / F2 1.91 / 1.00 st. K2 (radius 0.5) stopped the /z/ of "is" and clicked.
  - Sent `listen/video/thick_tongue_slow4.mp4`; `data/voice_thick_tongue_K1.params`.
- Owner: "we need to flex a support beam from within. several that overlap- one for the tip curling/molding, one for the base. its an interpolant problem with a spine around which the rest of the tongue moves. the base one can bend back like a fishing pole while keeping its forward end from raising up too much."
- **The spine tongue** (`tongue.spine` and spine_*; tongue.cpp spine_*):
  - Two beams of bends. The base beam is clamped at the root (anchor = root base plus its rest offset), stations 0 to spine_base_end (16), with butt pull toward the rest bends (spine_butt, fading forward) and a forward-end rise limit (spine_rise). The tip beam roots on the base beam at spine_tip_start (12); the spine blends linearly between the beams over 12–16.
  - Bends follow the shape's spine (its surface less spine_up along its normals) as damped second-order systems (base 10 Hz, tip 18 Hz, critical).
  - Segment lengths follow the shape's within ±spine_stretch over 20 ms.
  - spine_reach: damped least squares on the goal bends toward the shape's spine at the base end and at the tip.
  - The surface is the spine plus spine_up: 0.6 at the root, body_thickness 1.6 over the body, thinning to tip_radius. The cap is a half circle about the spine's end; the underside is spine_down below the blade.
  - The solid hull turns the joints behind a point inside (least-norm over the beams carrying it).
  - Off is bit-identical.
- Results (H4 base, tongue 5 score):

  | | F1 / F2 (st) | AO, AH F2 | mid / jitter | sealed N N N N N, D | target error rms |
  |---|---|---|---|---|---|
  | K1 thick | 2.15 / 0.98 | 924, 1159 | 0.014 / 0.0048 | 61 52 83 64 97, 58 | 0.535 |
  | spine, butt 2, rise 0.5 | 1.45 / 1.91 | 1421, 1542 | 0.024 / 0.0089 | none | – |
  | spine, butt 0, rise 1.5, stretch 0.4 (SR0) | 1.31 / 1.16 | 1220, 1219 | 0.025 / 0.0090 | 0 0 36 0 67, 31 | 0.616 |
  | reach both beams, 3 steps | broken (EH F2 3439; slow 0.22–0.30) | | | | 0.73 |
  | **SR5g**: reach tip only, damping 0.5, 1 step | 1.53 / 1.07 | 1214, 1215 | 0.028 / 0.0106 | 0.11 0.09 57 0.23 75, 36 | 0.625 |

  - The fishing-pole butt (spine_butt 2) held the root to a third of the plan's turn: back vowels fronted and no seals.
  - SR5g's /s/ is open (0.018); its /t/ after it is not sealed (0.018). Nothing is inside the palate (0.062 largest, 2 frames).
- Sent `listen/video/spine_tongue_slow4.mp4`, `spine_tongue_realtime.mp4`, `listen/m74_tube/banana_spine_tongue.wav`; `data/voice_spine_tongue.params`.
  - OBSERVED: the shapes are smooth with a fat tip and a free blade. At /n/ of "money" the spine tongue stays below the palate (the plan's dashed outline is above it); the body sits lower than the plan in /iy/.
- Owner on the spine tongue: "lol. terrible! show me your two articulating supports." Their sketch: a base support rising from the tongue's base low in the body up and forward through its core; a tip support from its end out to the tip; the base flexing back, the tip curling.
  - My "spine" had been the surface moved inward by a thickness, parallel to the dorsum.
- **The skeleton tongue** (`tongue.skeleton`, skeleton_*; `src/tongue_skeleton.cpp`):
  - Rest supports: the base from 0.35 of the way from the root base to the body's centroid, to the centroid plus 0.3 toward the tip cap's centre. The tip support from the base support's last quarter to the tip cap's centre. Four segments each.
  - Parameters: each support's four bends (cumulative) and a stretch (log, ±0.35). The tip support turns with the base support at its root.
  - Skin: linear blend skinning of the rest upper surface and the thick tongue's rest underside. Each point attaches to its nearest point on each support (segment, share, offset along and across). The tip share rises with smoothstep from section 11 to 17 (underside 1). The first point is the root base.
  - Goals: Gauss–Newton (forward-difference Jacobian, Levenberg 1e-3, steps clamped ±0.3) of the skin's upper points onto the planner's shape (points inside the hull held on it; tip weight 2), plus each parameter pulled toward rest (0.02) and the base end's rise limit. Two warm iterations per tick.
  - Parameters follow their goals critically damped (base 10 Hz, tip 16 Hz). Hull contact refits the parameters themselves toward the nearest hull points.
  - GEOLOG field 9 holds the supports; mouth_video draws them (base red, tip green).
  - Off is bit-identical.
- SK1 (H4 base): target error rms 0.444 cm (K1 0.535, spine 0.62). Nothing inside the palate. F1 / F2 1.58 / 1.64 st (AO 1327, AH 1446: back vowels fronted). N sealed 10, 0 (0.027), 17, 0 (0.012), 64 ms; D 22. Mid / jitter 0.030 / 0.0082. The root surface kinks in some frames.
- Sent `listen/video/skeleton_tongue_slow4.mp4` (`data/voice_skeleton_tongue.params`).
- Owner on the skeleton tongue: "its better. it needs careful finetuning. but it's better. in the is positioned especially wrong."
- Rise limit and tip weight (SK2–SK5: rise 2.5, tip weight 6, stiffness 0.005) changed little.
  - Measured against the motors' posture (GEOLOG "planned"), the tip missed about 0.9 cm in every variant. Faster supports (40/60 Hz) made it worse (1.29), and more fit iterations changed nothing.
  - The goal fit itself is good: the goal pose's tip lies 0.18 cm from its goals and 0.30 cm rms (TONGUELOG now logs both, after the planned area).
- **The cause in "in the": the targets run a phone early.** With the tongue's look-ahead (lead 0.03) on top of the planner's own leads (set for 28 ms motors):
  - the /n/'s apex target (6.75, 10.27) was held only from 1.575 to 1.625 s (the /n/ sounds 1.613–1.668);
  - the /dh/'s (7.80, 8.69) from 1.625 (mid-/n/) to 1.680;
  - the schwa's from 1.680 (11 ms into the /dh/).
- **SK9 = lead 0:**
  - F1 / F2 1.18 / 1.41 st (AO 1131, AH 1279).
  - Mid / jitter 0.013 / 0.0046 (SK5: 0.039 / 0.0138).
  - N sealed 46, 1 (0.010), 56, 45, 83 ms; S 12, T 19, D 52.
  - Frames in "in the": /ih/ tip curled down behind the lower teeth; /n/ blade to the teeth and ridge (dental before /dh/); /dh/ fat tip at the upper incisors' back; schwa lowered.
  - Slower supports (7/10 Hz, SK10) smoother in the slow band (0.046) but sealed less. Lead 0.01 (SK11) in between.
- Sent `listen/video/skeleton_tongue_2_slow4.mp4`, `_realtime.mp4`, `listen/m74_tube/banana_skeleton_tongue_2.wav` (`data/voice_skeleton_tongue_2.params`).
- Owner on skeleton tongue 2: "issue just after one second. harsh click."
  - At 1.049 s (the /l/ of "always"; transient ratio 13.9, slope 0.20) the tip moved 0.01 cm and the narrowest area held at the /l/'s side-channel floor (0.119 cm²).
  - Between 1.044 and 1.046 s the cavity in front of the contact jumped (one node 2.97 → 5.66 cm², the lip-side nodes 0.69 → 2.97). The vertical cuts in front of the raised fat tip switched their organ crossing from the cap to the floor beneath the blade as the cap's front passed them.
- **SK12 = SK9 + smooth_airway 1** (width bounded by the hull point's distance to the nearest organ, continuous as the cap passes):
  - Transients over 8: SK9 1.049 (13.9), 1.478, 1.950, 2.265 (22.3), 2.444. SK12: only 2.436 (9.0).
  - F1 / F2 1.61 / 1.59 st (AO 1250, AH 1331; SK9 1131, 1279).
  - N sealed 46, 4 (0.006), 57, 48, 83 ms; D 53.
  - Mid / jitter 0.0135 / 0.0045.
- Sent `listen/m74_tube/banana_skeleton_tongue_3_no_click.wav`, `listen/video/skeleton_tongue_3_realtime.mp4` (`data/voice_skeleton_tongue_3.params`).
- Owner: "playing faster than realtime?"
  - The video and audio are real time (3.12 s at 50 fps; 3.13 s at 48 kHz).
  - The speech is fast: the plan (v3_arrive50_fast_hold) speaks the banana sentence at 6.3 syllables/s, 65 ms per phone (read speech about 4–5).
  - With `say --rate 0.75`: 4.1 syllables/s (speech 2.93 s). Skeleton tongue 3 then seals every closure: N 83 58 100 62 122 ms, D 85, T 54, DH 12. Transients over 8: 1.655 (8.7), 2.533, 2.779, 3.150 (12.7), 3.413.
  - Sent `listen/m74_tube/banana_skeleton_tongue_3_rate75.wav`, `listen/video/skeleton_tongue_3_rate75.mp4` (score /tmp/es/banana/pt5_r75.score).
- Owner on the 0.75-rate banana: "its our best result so far."
- **Sustained vowels** (owner: "sustain and hold aaaaaa, eeeee,,,ooooo"):
  - `say "ah." / "ee." / "oh."` (rate 0.75) with the plan's duration scales AA 5.4, IY 6.4, OW 5.0 (/tmp/es/sustain/hold.plan): each about 1.6 s. Spoken by skeleton tongue 3 (SK12.params).
  - LPC every 100 ms through the hold:
    - ah: F1 704–753 (target 732), F2 1396–1423 (target 1188).
    - ee: F1 363–381 (321), F2 2044–2074 (2177).
    - oh: F1 454–515 (465), F2 877 → 1031 (983; glides up over the hold).
  - Pitch 137 → 127 Hz. No wobble in the formant tracks.
  - Sent `listen/m74_tube/sustained_vowels.wav` (the three, 2 s each) and `listen/video/sustained_vowels.mp4`.
- Owner on the sustained vowels, with a spectrogram of their own aaa / eee / ooo: "notice how it is harder to hold eee. also- for all three. he sounds like a sick ship. this is a vocal cord optimization problem."
  - Their spectrogram: strong low harmonics; formant bands above them separated by dark valleys; fine noise speckle; the /iy/ harmonics waver.
  - Ours at the same scale: energy nearly even to 8 kHz, valleys filled, a fine regular ripple through the upper harmonics.
  - (`voice/ref/user_three.wav` is not those three vowels.)
- **Source measures** (`analysis/source_measures.py`; ours on the held vowels, bdl on 400 ARCTIC files' stressed aa / iy / ow, middle 60% of tokens ≥ 120 ms):

  | | tilt dB/oct | valley dB | H1–H2 dB | HNR dB | jitter % |
  |---|---|---|---|---|---|
  | ours aa / iy / ow | −9.1 / +0.3 / −7.8 | 20 / 15 / 25 | +5.3 / +3.0 / +5.4 | 30 / 30 / 27 | 0.21–0.27 |
  | bdl aa / iy / ow | −8.9 / +9.5 / −8.8 | 21 / 22 / 25 | −2.7 / −6.9 / −3.3 | 19 / 16 / 20 | 0.37–0.65 |

- **Source fit** (`analysis/sourcefit.py`, /tmp/es/srcfit): blockwise Levenberg–Marquardt (damping, step radius, kept only if the cost falls) over four blocks: posture (voicing_lower / upper / chink), fold cover (mass, stiffness, loss, contact_softness), glottal noise (aspiration, fold_turbulence, closing_turbulence), variation (tension, posture). 1 s holds of /aa iy ow/; residuals scaled per measure against bdl's medians.
  - Cost 7.84 → 5.12 in 6 sweeps. HNR 29 → 23–26; H1–H2 stuck at +2.5 to +4.5; jitter about 0.27.
  - The fit chose aspiration 5.3, fold_turbulence 1.9, closing_turbulence 0.6, cover_mass 0.071, cover_stiffness 44000, contact_softness 0.00052, voicing_lower 0.0085, upper 0.0052, chink 0.00018.
- The folds close every cycle (minimum area = the 0.0002 cm² chink). They vibrate near 138 Hz against an intended 118.6 (the ear is off in this tree).
  - Glottal flow (WAVEDUMP now prints the flow and area): open quotient 0.54. The flow rises steadily, sags slightly and cuts off abruptly. Flow H1–H2 = 10.6 dB (ordinary); the radiated +4.5 follows.
  - The flow already includes the glottis's air mass and the tract's wave impedance (cord.cpp). H1–H2 may partly reflect bdl's recording chain; men's H1–H2 commonly runs 0 to +5.
- Candidates (1 s holds; HNR and jitter on two windows, noisy):
  - aspiration ×3–6 and fold_turbulence ×3–6 lowered HNR only about 3 dB (25.5–28). The glottal noise is weak.
  - variation.tension 0.06 with posture 0.05 brought jitter to 0.33–0.47% (bdl's range) but raised pitch to about 157 Hz and H1–H2 to 6–11.5.
- Sent `listen/m74_tube/sustain_0_before.wav`, `sustain_1_source_fit.wav`, `sustain_2_breathier.wav` (aspiration 32, fold_turbulence 12), `sustain_3_breathier_irregular.wav` (plus tension 0.06, posture 0.05); params `data/voice_source_best/B2/B4.params`.
- Owner on /ee/: "holding this is actually bending the lower lip down a little bit, closing the mouth partially, so that the tongue has a narrower gap above it but is not itself tilted up. the effort here is in pushing the tongue up exercising the lower neck muscles, which in turn increases the feedback to the vocal cords."
  - OBSERVED (held /ee/, skeleton tongue 3): jaw 0.77°, tongue body high under the palate (narrowest 0.21 cm²), tip curled down behind the lower incisors. The articulation matches; no coupling to the larynx existed.
  - Also OBSERVED: the held /oh/ is malformed. The base support stands nearly vertical and the front of the tongue kinks.
- **The tongue's pull on the larynx** (speaker; `larynx.tongue_pull_tension`, `_adduction`, `_raise`): Tongue::lift() = mean height of sections 8–16 above rest (jaw frame), followed over 20 ms. Per cm of lift: the intended pitch × (1 + tension), the folds' rest gaps − adduction, larynx_drop − raise (not above rest). Off is bit-identical.
  - Lift in the held vowels: ah 0.00, ee 0.78, oh 0.52 cm.
  - TP3 (tension 0.08, adduction 0.01, raise 0.4), 1 s holds: pitch ah / ee / oh 152 / 159 / 156 Hz (without: 152 / 150 / 149); H1–H2 +4.5 / −1.4 / +0.8 (without: +4.5 / +2.5 / +3.2). bdl's aa / iy / ow: −2.7 / −6.9 / −3.3, the same ordering.
  - The larynx block must be appended as a block (cand.sh edits only existing blocks; the first TP renders were identical to the source fit).
- Sent `listen/m74_tube/sustained_tongue_pull.wav`, `listen/video/sustained_tongue_pull.mp4`, `listen/m74_tube/banana_tongue_pull_rate75.wav` (`data/voice_tongue_pull.params`). Banana transients over 8: 2.533, 2.653, 3.366, 3.404, 3.599 (22.9, the final /d/'s release).
- Owner: "fix the oh."
  - OBSERVED (held /oh/, TP3): the base support stood nearly vertical and section 16 spiked out of the dorsum, with a notch in front. Linear blend skinning: a point in the overlap (tip share rising over sections 11–17) averaged two places when the supports diverged (base up, tip forward and down).
- **Curve skinning** (`tongue.skeleton_curve`): the supports are blended as curves. The central curve is the base support to the tip support's root (15 points), an overlap where the base support's last quarter fades into the tip support's start (10 points, smoothstep), then the tip support (21 points). Each skin point rides that one curve by its rest offset along and across it, with tangents blended across vertices.
  - `skeleton_bend_limit`: a fit residual on the base support's whole turn beyond the limit (no effect at 0.5–0.8 here).
  - Smooth shapes, no spike. But /oh/ could not retract its back (F2 4.54 st; TP3's 0.31 came from its distorted shape).
- **The anchor rides the hyoid** (`skeleton_anchor_move` cm each way, `skeleton_anchor_stiffness` per cm). Two more parameters move the base support's anchor; they follow at the base support's rate.
  - `skeleton_smoothness`: fit residuals against a bend differing from its neighbour's (the tip support zigzagged).
  - Off, both are bit-identical (zero Jacobian columns and zero rows).
  - AM2 (move 1.0, smoothness 0.1), held vowels: shapes hug the plan in /ah iy ow/. F2 AA 1247 (1188), IY 2170 (2177), OW 1099 (983). F1 errors IY 3.2 st and OW 4.1 st; /ee/ closes to 0.07 cm².
  - The vowel postures were fitted through tongue 5, whose swelling is baked into them. Refitting them through this tongue: `organ_vopt3.py` (now takes SPEAK from the environment) with build_x speak and AM2.params, W 10, MAXFEV 260, all twelve vowels → /tmp/es/vskel/organs_skeleton.txt.
- **Vowels refit through the skeleton tongue** (`data/organs_skeleton.txt`): formant error before → after:
  - aa 1.34 → 0.07; iy 10.4 → 1.20; ow 2.94 → 1.79 (F 514 990); ae 0.85 → 0.13; ah 3.75 → 0.61; ao 4.89 → 0.80;
  - eh 0.60 → 0.22; ey 3.80 → 0.50; ih 2.90 → 0.11; uh 6.05 → 0.92; uw 2.97 → 1.61; ax 2.88 → 0.10.
- **Regression found and fixed: contact bounce.**
  - With the moving anchor and the new vowels, each /n/ sealed, sprang off the palate to 0.7–0.8 cm² and came back about every 30 ms (N sealed only 9–23 ms; transients 39–53 at 2.6–3.1 s).
  - At contact (1.924 → 1.928 s) the base support's end dropped 1.2 cm in one millisecond: the contact refit carried the full regularization (rest pull, smoothness, rise limit) against a weak hold.
  - Keeping the anchor out of the contact refit alone did not help.
  - **`skeleton_minimal_contact`**: the correction holds each parameter to its present value (0.2), drops smoothness, rise and bend limit, leaves the anchor, and turns each bend at most 0.1 rad.
  - Banana: N sealed 78 115 103 64 113 ms, D 86, T 55, DH 69. One transient over 12 (2.534, ratio 14). Nothing inside the palate.
- **The held /oh/:**
  - The planner's OW glide (toward UH) took 60% of the vowel from 35% in: a 1.6 s /oh/ slid for a second. Now `timing.diphthong_glide_max` (0.18 s) caps it, ending at 95% of the vowel. Diphthongs under 0.3 s are unchanged (the banana plan is bit-identical).
  - Path dependence: /ow/ held from rest gave F 491–525 / 956–999 Hz; arriving from a schwa, 529–558 / 977–1062; the planned word, 587 / 1150. **`skeleton_restart`**: every 10 ms the goals are also fitted afresh from rest (6 iterations), keeping the lower cost.
- **Skeleton tongue 4** (`data/voice_skeleton_tongue_4.params`, organs_skeleton), held vowels (F1-F2 every 0.3 s):
  - ah 719–742 / 1219–1235; ee 332–338 / 2073–2096; oh 517–559 / 981–1059, gliding only at its end.
  - Motion in the holds: slow 0.003–0.013, mid 0.001–0.006.
  - Banana: F2 0.81 st; N 78 115 106 70 120 ms, D 86; transients over 12 at 2.532 (23.9) and 3.602 (the final /d/).
- Sent `listen/m74_tube/sustained_fixed_oh.wav`, `listen/video/sustained_fixed_oh.mp4`, `listen/m74_tube/banana_skeleton_refit_vowels.wav`.
- Owner on skeleton tongue 4: "he still sounds somewhat like a sick sheep. ooo does not need to pull the tongue back. the banana stand has completely collapsed." Then: "stand is the only word that actually articulates in this recording."
- **The collapse, traced** (banana at rate 0.75, the "banana stand" stretch):

  | | voiced | median f0 | 10th-percentile level |
  |---|---|---|---|
  | best take's tree (L0) | 77% | 131 Hz | −38.5 dB |
  | + source fit | 76% | 146 Hz | −39.4 |
  | + tongue pull | 73% | 148 Hz | −43.2 |
  | + tongue changes | 75% | 151 Hz | −39.0 |
  | + new vowels (sent) | 67% | 149 Hz | −49.2 |

  - Closures and formants did not show it.
  - Per phone, the sent take was 2–4 dB quieter in every vowel and voiced consonant while /s/ held its level: only "stand" stood out.
  - The source fit raised f0 about 15 Hz (the fit held no pitch; the ear is off).
  - L0's render no longer matches SK12r bit for bit (diverges at 1.24 s; not traced).
- **Skeleton tongue 5** = the best take's source (SK12) + the skeleton tongue 4 mechanics + the refit vowels (`data/voice_skeleton_tongue_5.params`). Banana: levels back to L0's; F1 / F2 1.19 / 0.62 st; N 75 109 104 71 121 ms, D 86, T 54, DH 66; transients over 12 at 3.152 (23) and 3.599 (the final /d/, 75).
- **"Ooo" with the tongue forward:** organ_vopt3 bounds for uw and ow (owner's rule): tx ≥ 3.6, jaw ≤ 8°, protrusion ≥ 0.6, aperture ≤ 0.6.
  - UW: jaw 4.26, tx 4.14, pro 1.00, ap 0.17 → F 338 / 1222 / 2269 (bdl 337 / 1238 / 2496).
  - OW cannot reach its F2 without retraction (1193 vs 974, 4.0 st) and keeps its earlier posture.
  - `data/organs_skeleton.txt` now holds the new UW.
  - Held (F1-F2): ah 715–728 / 1201–1211; ee 335–362 / 2079–2103; oo 326–349 / 1207–1246. Motion in the holds below 0.003.
- **The sheep:** not found yet. Ruled out on the held vowels:
  - slow tremor: pitch 0.6–0.8% rms (6 Hz), loudness 0.4 dB (`analysis/tremor.py`);
  - roughness: pulse-picked from the audio the same way (`analysis/pulses.py`), ours jitter 0.9–1.1%, shimmer 2–4%; bdl's vowels 1.1–1.5%, 5–7%;
  - nasality: the velum closes (velum_leak 0);
  - period doubling: successive period changes correlate −0.50, as independent jitter does; the half-pitch subharmonic is 23–31 dB down.
- Sent `listen/m74_tube/sustained_ah_ee_oo.wav`, `listen/video/sustained_ah_ee_oo.mp4`, `listen/m74_tube/banana_skeleton_5.wav`, `listen/video/banana_skeleton_5.mp4`.
- Owner on skeleton tongue 5: "still a bit of a gurgle/death rattle. but better. the mouth is at least correct, its vocal cords that need tuning for the sustained. for the banana stand: stand sounds clear. there is always sounds like he has something in his mouth. the z especially at 1.62 is unable to reach his proper position correctly."
- **The /z/s:** neither narrowed (core narrowest 0.43 and 0.61 cm²); the goal fit itself missed the tip by 0.68 and 1.00 cm.
  - `skeleton_smoothness` 0.1 resisted the tip support's curl. Off: the second /z/'s miss fell to 0.23 cm and its core to 0.060.
  - The first /z/ ("is") needed the tongue longer: `tongue.skeleton_stretch` (was the constant ±0.35) at 0.5 → misses 0.09 / 0.10 cm, cores 0.060 / 0.062 (the groove). At 0.7, 0.05 / 0.11.
- **Skeleton tongue 6** (`data/voice_skeleton_tongue_6.params`; smoothness 0, stretch 0.5):
  - Banana F1 / F2 1.14 / 0.47 st; DH (sentence start) 0.035 cm²; Z 0.054 cm²; N 76 109 105 68 123 ms; D 86; nothing inside the palate.
  - Transients over 12 at 3.151 (19.6) and 3.598 (the final /d/, 90).
  - Held vowels unchanged (ah 718–727 / 1200–1218, ee 333–369 / 2097–2116, oo 330–344 / 1194–1219).
- Sent `listen/m74_tube/banana_skeleton_6_z.wav`, `listen/video/banana_skeleton_6.mp4`.
- **The rattle** (`analysis/rattle.py`: 10–60 Hz flutter of 1 ms log envelopes, whole band and 2–5 kHz, with pulse jitter and shimmer), held /ah/, one source off at a time:
  - all on: 0.42 / 0.72 dB, jitter 0.96%, shimmer 4.1%;
  - walls off: 0.41 / 0.73; muscle fluctuation off: 0.31 / 0.65; glottal noise off: 0.40 / 0.67;
  - **elasticity wander off: 0.13 / 0.46, jitter 0.16%, shimmer 0.7%**; everything off: 0.10 / 0.14.
  - The wander was noise low-passed at fold.elasticity_corner = 20 Hz, so its energy sat in the rattle band.
  - Corner 2 / 4 / 8 Hz (wander 0.05): flutter 0.20 / 0.21 / 0.34 dB, jitter 0.38 / 0.56 / 0.67%, shimmer 1.7 / 2.2 / 3.1%. Wander 0.1 at 4 Hz: 0.37 dB, 1.09%, 4.6%.
  - Adopted corner 4 Hz (pitch tremor 0.68%, 6 Hz).
- Owner: "always has a weird behavior still. and after that- we can finetune on nato. this is getting better. on teh vocal cords: adding more tube sections below to add resonant feedback from the lungs may be important here."
  - OBSERVED: in "always" the /w/ sealed the mouth (oral minimum 0.00–0.02 cm², 1.33–1.42 s); the plan itself had the dorsum on the palate.
  - The planner gives /w/ the /uw/ body height, and banana.score had been planned with the earlier UW (ty 8.59, jaw 12.5).
  - Re-planned with the refit UW (`banana6.score`): /w/ 0.259 cm², /l/ 0.087, AO 0.43, EY 0.43. F1 / F2 1.07 / 0.43 st. Closures unchanged; nothing inside the palate.
- **Skeleton tongue 7** (`data/voice_skeleton_tongue_7.params`, organs_skeleton with the new UW). Sent `listen/m74_tube/banana_skeleton_7.wav`, `listen/video/banana_skeleton_7.mp4`, `listen/m74_tube/sustained_slow_wander.wav`.
- **The subglottal system** (owner's suggestion), `subglottal.bronchial_length`, `bronchial_area`, `lung_port_scale` (BodyLayout):
  - The trachea (12 cm, 2.5 cm², wider at the top) gains a bronchial tube below it whose total area grows geometrically to bronchial_area at the lungs. The lungs' end is a port lung_port_scale times wider than the last section.
  - Before: a matched (absorbing) end, so no subglottal resonances.
  - Transmission-line estimate of the impedance peaks at the glottis (/tmp/es/ladder): trachea alone with an open end 568 / 2144 Hz; bronchi 6 cm to 40 cm², ×3: 496 / 1728; 10 cm: 456 / 1414 / 2058. Adult men measure about 550–650 / 1350–1550 / 2200; our first sits low with the 12 cm trachea.
  - With the lung port ×3 at 40 cm² the voice died within 0.9 s. The steady-flow correction (impedance − airway resistance) × flow went negative (port impedance 0.35 < 1.5) with a loop gain above 1. Now clamped at zero; the existing setups (port impedance about 16) are bit-identical.
  - SG8 (8 cm, 40 cm², ×3), banana: F1 / F2 0.98 / 0.52 st; closures unchanged; the final /d/'s transient 98 → 31.
  - Held /ah/: H1–H2 5.5 → 4.7, HNR 35.5 → 31.8, rattle flutter 0.21 → 0.26 dB, jitter 0.56 → 0.57%, shimmer 2.2 → 5.6% (bdl 5–7%).
- **Skeleton tongue 8** (`data/voice_skeleton_tongue_8.params`). Sent `listen/m74_tube/sustained_subglottal.wav`, `banana_subglottal.wav`.
- Owner: "always still has a sort of wrong sound to it. the nato finetunes should help ferret it out im sure."
- Owner on the subglottal renders: "sounds the same- still with a death rattle."
- **Rattle isolation** (held /ah/, SG8): `listen/m74_tube/rattle_isolation.wav`, four in a row:

  | | flutter 10–60 Hz (whole / 2–5 kHz) | jitter | shimmer |
  |---|---|---|---|
  | 1 current | 0.26 / 0.76 dB | 0.57% | 5.6% |
  | 2 steady (wall, wander, muscle off) | 0.12 / 0.62 | 0.16% | 1.0% |
  | 3 quiet (aspiration, fold and closing turbulence, frication off) | 0.26 / 0.48 | 0.61% | 5.5% |
  | 4 both | 0.09 / 0.16 | 0.16% | 0.3% |

  Asked the owner which still rattle.
- **NATO fine-tune launched** (/tmp/es/nato): policy_fit als, 2 sweeps, ORDER=alternate. Groups organs, plan, tongue (`ENGINE_TONGUE_KEYS` "tongue.skeleton_,tongue.tip_radius"), `ENGINE_RINGING=0` (the hinge is now switchable). From SG8.params, organs_final (= organs_skeleton), the v3 plan. 269 neurons in 20 blocks; start cost 12.65.
  - The main tract/build was rebuilt to current source and reproduces SG8 bit for bit.
- Owner on the four rattle takes: "it smoothly shifts from raspy zombie robot to squeezed chicken toy, with no sign of a human at any stage."
  - So the variation sources supply the rasp; the bare fold oscillation is the squeezed toy.
- Third-octave spectra (dB re 500 Hz) against bdl's tokens (`analysis/thirds.py`):
  - /ah/ within ±4 dB to 4 kHz (no twang peak at 2.5–4 kHz; 5 kHz −10).
  - /ee/ 630–2000 Hz +4 to +13 dB (the F1–F2 valley too full).
  - /oo/ 1–2 kHz +5 to +12.
  - (The 314 Hz band falls between our steady harmonics.)
  - With every side branch off (`BRANCHES`), the excess stays (+3 to +11); the piriform shapes only 5 kHz.
- **NAQ by inverse filtering** (`analysis/naq.py`: IAIF at 16 kHz; per cycle (flow peak-to-peak / largest closing rate) / period; the closed share):

  | | ours | bdl (quartiles) | closed share ours / bdl |
  |---|---|---|---|
  | ah | 0.080 (both off 0.079) | 0.098 (0.089–0.103) | 0.15 / 0.15 |
  | ee | 0.105 | 0.125 (0.110–0.142) | 0.11 / 0.11 |
  | oo | 0.133 | 0.123 (0.110–0.145) | 0.11 / 0.10 |

  - The direct-flow estimate (0.04) is not comparable: it is not inverse-filtered.
- Less pressed (rest gaps voicing_lower / upper 0.011 / 0.006 →):
  - 0.015 / 0.009: NAQ 0.085 / 0.118 / 0.140.
  - **0.019 / 0.012 (open2): NAQ 0.092 / 0.132 / 0.149**; /ah/ flutter 0.29 / 1.09 dB, jitter 0.64%, shimmer 4.2%.
  - **open2g** (+ wall 0.05, muscle fluctuation 0.005, elasticity wander 0.025): NAQ the same; jitter 0.33%, shimmer 2.2%.
- Sent `listen/m74_tube/sustained_less_pressed.wav` (open2) and `sustained_less_pressed_gentle.wav` (open2g); params `data/voice_source_less_pressed*.params`.
- NATO fine-tune: block 2 of 20 after about 20 minutes (plan.duration 12.81 → 12.34; organs.AA unchanged).
- Owner on the less-pressed takes: "the gentle sounds a faint bit more human. but its still massively missing what makes a voice human. which is that it is way more resonant with a lot more layered sound."
- **Formant bandwidths** (`analysis/bandwidths.py`: LPC at 10 kHz, order 12, 25 ms frames; F1–F3 the three lowest poles above 150 Hz with bandwidth under 700):

  | | ours (open2g) | bdl |
  |---|---|---|
  | ah | 119 / 91 / 97 Hz | 153 / 100 / 120 |
  | ee | **185** / 74 / 120 | **40** / 73 / 174 |
  | oo | **215** / 54 / 96 | **59** / 180 / 123 |

  - Ruled out as the cause of the broad F1:
    - pitch (110 / 120 / 150 Hz: 195 / 185 / 218);
    - walls.loss_scale 0.1–0.3 and tract.hop_loss 0.005;
    - wall stiffness ×3 to ×100 (lowers F1, bandwidth 164–227);
    - subglottal off (204) and skin off (171);
    - the oral jet (frication, jet_growth: identical).
  - The pressed source narrowed it a little (151 / 197).
- **Tube alone** (`tools/transfer`, closed glottis, Story's shapes; now VISCOUS and HOPLOSS):
  - iy / uw / aa F1 bandwidth 110 / 140 / 120 Hz; rigid walls 50 / 60 / 60.
  - The links' viscous resistance (`viscous_resistance_per_cm`: max(0.03 A^−2.6, 0.31 / A) per cm) is a steady-flow law. A boundary layer at a formant frequency in 0.2 cm² gives about 0.4 per cm against its 1.97.
  - VISCOUS 0.2: rigid F1 bandwidth 20–30 Hz (the radiation and boundary-layer floor); with walls 90 / 110 / 70 Hz. Walls still add 60–80 Hz.
- **`tract.viscous_scale`** (BodyLayout → Tube::set_viscous_scale; default 1, bit-identical).
  - V1 = open2g + 0.2: held F1 / F2 / F3 bandwidths ee 130 / 45 / 84, oo 170 / 52 / 53, ah 91 / 87 / 77. Wall loss ×0.5 on top changed nothing.
  - V3 = the pressed source + 0.2: ee 113, oo 146, ah 97.
  - Banana (V1): F1 / F2 0.58 / 0.52 st; closures unchanged; one transient over 12 (3.598, the final /d/).
- Sent `listen/m74_tube/sustained_resonant_gentle.wav` (V1), `sustained_resonant_pressed.wav` (V3), `banana_resonant.wav` (V1); params `data/voice_resonant_gentle.params`, `voice_resonant_pressed.params`.
- Owner on the resonance takes: "its about the same. a bleating sheep." Then, interrupting a three-strip fold trial: "stop trying to make a fold oscillator. make a tube that has a pinched section, with multiple resonant sections before and after that follow our design pattern. the fold doesnt exist. there is no fold. there is only resonance interacting with elasticity."
- History: the sheet tube's elastic section (milestone 3) never closed, a balloon neck at ±13% of an open ring. The slit cord and heavy walls followed, then folds inside the sheet.
- **The pinch tube** (`include/tract/pinch_tube.hpp`, `src/pinch_tube.cpp`; `pinch.*`; Cord delegates to it when pinch.on; off is bit-identical):
  - Six wall sections (slit 1.3 cm across, 0.4 cm each) in the owner's pattern (m / e / s: 1/2/1, 2/4/0.1, 2/1/0.3, 4/0.1/3, 0.1/3/0.1, 1/2/1). Rest gaps 0.8, 0.25, then 2 × the posture's lower and upper rest (the pinch), then 0.25, 0.8; the chink a parallel leak.
  - Each wall: mass, support to rest, elasticity to its neighbours (both times the tension, (pitch / natural_pitch)²; elasticity wanders, corner 4 Hz, tamped by support), loss, contact where the walls cross.
  - Air: one flow through the chain from the arriving waves (the sections' air mass, slit friction, the kinetic pressure into the narrowest, 20% recovered past it). Each wall feels Bernoulli's pressure up to the narrowest section and the outlet's past it.
- First render (support 30000, elasticity 30000, density 0.1): blown open, a balloon again: area 0.12–0.17 cm², flow 365–475 cm³/s, a slow 70 Hz wobble.
- Sweeps (`analysis/pinch_sweep.py`; held /ah/, WAVEDUMP 1.00–1.10 s):
  - Support 100000, elasticity 30000, density 0.1: self-oscillation at 189 Hz, flow 43–194, strength 0.94.
  - Support 60000, elasticity 10000, density 0.15: 128 Hz, area 0.014–0.091, flow 15–280, not closing. Through the hold the flow was nearly sinusoidal (24–264 cm³/s), never shut.
  - Pressing the rest gaps (posture 0.004/0 down to −0.004/−0.008) shrank the swing and raised the pitch.
  - **Loss 0.2 → 0.1 closes it:** rest 0.004/0, area down to 0.0002 (the chink), flow 0–308 cm³/s, 129 Hz, strength 0.91. Loss 0.05 overshoots (negative flow).
- **Pinch tube 1** (`data/voice_pinch_1.params`):
  - Held vowels voiced 100%: f0 119 / 116 / 109 Hz; F1-F2 ah 726–732 / 1206–1211, ee 393–408 / 2094–2111, oo 382–391 / 1215–1240.
  - NAQ 0.119 / 0.150 / 0.138 (bdl 0.098 / 0.125 / 0.123); H1–H2 10.0 / 6.8 / 6.5; HNR 22 / 25 / 23 dB (the fold 32–35; bdl 16–20).
  - Banana: voiced 80%, f0 131 (102–148), F1 / F2 1.44 / 0.76 st.
- Sent `listen/m74_tube/sustained_pinch.wav`, `banana_pinch.wav`.
- NATO fine-tune stopped (organs.ER, sweep 0, cost 12.06) once the source was to be replaced.
- Owner on pinch tube 1: "so sickly. maybe a bit better. give me some sliders! lemme finetune for you."
- **Sliders** (`analysis/slider_server.py`, on the Mini, bound to 127.0.0.1:8765; the laptop reaches it through `ssh -N -L 8765:127.0.0.1:8765 m4mini`; opened in the app's browser pane):
  - One slider per pinch setting (walls, geometry, the six sections' three pieces, wander), the pinch's rest posture and chink, breath noise, turbulence on the walls, tract friction and tract wall loss. Log scales where the range spans decades.
  - Utterances: held ah / ee / oo and the banana sentence (`data/slider_scores/`). Render on release, about 1 s; the last six takes are kept to A/B.
  - The base is the engine's whole tree for `data/voice_pinch_1.params` (speak PARAMS_OUT), so the sliders start at the engine's own values. Its default render is bit-identical to pinch tube 1.
  - Presets: `tract/data/sliders/NAME.params` (full tree) and `NAME.json` (the values) on the Mini. Every render's values are appended to `render_log.jsonl` there.
- Owner, after 216 slider renders (all held ah): "saved ee. disregard section settings. no matter what i do i cant disentangle oscillatory rattle from frequency- higher frequency voices, no rattle. and yes, real voices have a sort of rattle, but its more complicated than what i hear here".
- **OBSERVED — the saved ee does not phonate** (`analysis/rattle_pitch.py`; flow and opening dumped over 1.0–2.0 s):
  - The pinch is blown open: 0.97 cm², 900 cm³/s, no cycles. Along the render log (every 12th render re-rendered) it phonated with closure every cycle through render ~132. The last section moves (renders 144, 156, 215) blew it open.
  - ee's non-section values with the default sections ("ee-ns") phonate with closure.
- **MEASURED — the rattle against pitch** (ee-ns, held ah, tension anchor varied):
  - F0 75 / 93 / 111 / 128 Hz: 10–60 Hz flow flutter 1.27 / 0.49 / 0.25 / 0.13 dB, between-pulse drop 34 / 24 / 22 / 14 dB.
  - Flow jitter stays 0.2–0.3%, shimmer 0.5–0.7%.
  - Pitch: anchor 120 gives 93 Hz against the intended 119. Anchor 92 gives 119.5.
- **MEASURED — real voices drop as deep between pulses but vary far more** (`analysis/envelope_depth.py`, the 80 Hz–6 kHz Hilbert envelope; pulses.py jitter and shimmer):
  - bdl, 120 files: median drop 21–23 dB at F0 95–130 Hz. Owner (`voice/ref/user_three.wav`): 17–26 dB. Ours at matched F0: 22–24 dB. Depth alone is not the rattle.
  - Ours is too uniform: median 22, 90th percentile 24 dB. bdl: 23 and 28.
  - Wav level, owner's held vowels at 96–116 Hz: jitter 9–14%, shimmer 12–19%, 1 ms envelope flutter 0.6–1.4 dB. ee-ns at 93–128 Hz: 0.4–1.2%, 3.6–5.3%, 0.2–0.4 dB.
  - HYPOTHESIS: the "oscillatory rattle" is a strictly periodic pulse train. Below about 100 Hz it is heard as a mechanical buzz; higher, it fuses into a tone. Real voices are rattly too, but irregularly from cycle to cycle ("more complicated").
- **`pinch.jet_turbulence`, `pinch.turbulence_corner`** (default 0, bit-identical, verified with cmp):
  - The separated jet's eddies press on the walls past the narrowest section: rms share × ½ρv² (v = flow / pinch area), one-pole noise per section at the corner (1 kHz).
  - The walls are shaken hardest while the jet runs fastest and not at all while shut: a body-driven cycle-to-cycle variation, not a prescribed jitter.
  - Sweep, ee-ns at anchor 92 (F0 119.5): share 0 / 0.1 / 0.2 / 0.4 / 0.8 gives flow jitter 0.18 / 0.37 / 0.58 / 1.05 / 1.92%. Wav level at 0 / 0.2 / 0.4: jitter 1.2 / 1.6 / 2.9%, shimmer 3.7 / 7.1 / 19.5%, flutter 0.24 / 0.37 / 0.95 dB.
  - HYPOTHESIS, unheard. A share of 0.3–0.4 is at the high end of plausible near-jet pressure fluctuation.
- Sliders restarted on `data/sliders/ee_ns.params` (ee, default sections, anchor 92, jet eddies 0.3), with a new Jet group (eddies, eddy speed). "turbulence on the walls" relabelled "turbulence noise in the sound": it is the acoustic glottal noise, not a force on the walls.
- Owner: "feedback and resonance loops. pretend this is a musical instrument like a harp. 1/n2 type stuff."
- **String walls** (`pinch.string_modes`, `string_tension`, `string_loss_growth`, `string_taper`; 0 modes = pistons, bit-identical, verified with cmp):
  - Each section's wall is a string across the slit, fixed at both ends, as K standing waves sin(nπx/W), K ≤ 12. The opening and contact are evaluated at 16 points; the opening is a lens.
  - Per mode: half the wall's mass, support and neighbour elasticity, plus the string's T(nπ/W)²(W/2), growing as n² (harmonic frequencies). T is multiplied by the pitch tension.
  - Uniform pressure drives the odd modes. Contact at each point drives every mode by its shape there (the pluck, about 1/n² displacement). Taper (an additive wedge in the rest gap) breaks the symmetry so even modes can sound.
  - Render about 2.4 s for the held ah (pistons 0.5 s).
- **MEASURED** (ee-ns, jet 0.2, 8 modes; held ah):
  - The phonation window is narrow: T 1000 → static (opening 0.073); 3000 → 81 Hz, irregular, rarely closing; 10000 → 177 Hz at anchor 92; 30000 and above → shut. Anchors for about 119 Hz: T 6000 / 110, 10000 / 137, 15000 / 160, 20000 / 185.
  - Flow harmonic slope (n = 1..16): pistons −15.9 dB/oct; T 6000 −14.1, 10000 −14.3, 20000 −12.9 (1/n² = −12). Flow H2: −17 → −16 / −13 / −9 dB.
  - Sound H2..H8 come up 3–15 dB (more layered). Flow jitter 0.4–1.1%, subharmonic −40 to −27 dB.
  - OBSERVED, unexplained: ear.pitch_feedback 0.08 changed nothing (output identical).
- Sliders restarted on `data/sliders/ee_harp.params` (ee_ns, 8 modes, T 10000, anchor 137, jet 0.2), with a Harp strings group. A/B with pistons: modes 0 and anchor 92.
- Owner on the string walls: "its the same. when i adjust them, it adjusts the entire voice. every slider does the same thing. throw out the entire tube. the harp point was that vocal cords are a string instrument - it is the output that needs harmonics. not an extra component. figure out what vocallab does for their vocal cords". The pinch tube is abandoned (code kept behind pinch.on = 0).
- **OBSERVED — VocalTractLab's vocal folds** (Birkholz, Drechsel & Stone 2019, Interspeech, vocaltractlab.de/publications/birkholz-2019-interspeech.pdf; VTL 2.4 manual §2.1 and the glottis dialog; math only, no code read):
  - Three models: "Geometric glottis 2025" (the default; unpublished; the 2019 model plus a jitter parameter after Fraj, Schoentgen & Grenez 2012), "Geometric glottis 2019", and the self-oscillating "Triangular glottis" (Birkholz, Kröger & Neuschaefer-Rube 2011). The classical two-mass model was removed. The "aspiration strength" parameter was removed as unnecessary.
  - The geometric glottis is kinematic: the motion is prescribed. The lower (i = 1) and upper (i = 2) edges: x_i(y, t) = x_rest,i (1 − y / L) + x_max sin(πy / L) sin(φ_i(t)), y along the fold. The paper's words: like "ideal strings that are fixed at their two ends", but only the first mode in space and in time.
  - x_max = c · x_rel · sqrt(P_trans) · L0 / L, with c set so x_max = 1 mm at 800 Pa.
  - Phases α1 = 2π∫f0 dt, α2 = α1 − φ (upper edge lag), taken mod 4π. A piecewise-linear warp (s ∈ [−1, 1]) skews the pulse. Double pulsing d ∈ [0, 1] shrinks every second cycle's opening: x_dyn = x_max sin(πy / L)[(1 − d / 2)(sin φ + 1) − 1] for 1.5π < φ < 3.5π.
  - A_i = 2∫max(x_i, 0) dy + A_chink. The glottis is two tube sections of length T / 2 with areas A1 and A2. The flow comes from the tube's aerodynamics, so the harmonics come from the clipped area pulse and the tract's air mass, not from the folds.
  - f0 carries small quasi-random flutter after Klatt & Klatt 1990. Eqs. 1–2 (L and T from f0, f_nat) did not survive my text extraction of page 1.
  - Most natural male setting (30 listeners, 4050 stimuli): x_rest,1 0.3 mm, x_rest,2 0.6 mm (slightly divergent), φ 70°, s 0, d 0.05, chink 2 · x_rest,1 · 4 mm = 2.4 mm², L0 16 mm, T0 4.5 mm, f_nat 120 Hz, P_sub 800 Pa. It scored 2.2 of 3: "rather natural", not natural.
  - Ours: `Cord::set_kinematic` (2026-10-15) is this model without skew, double pulsing or flutter; the owner called it "Perfect Paul". The Cord's self-oscillating mechanics are a port of VTL's triangular glottis.
- Owner: "study what we did for waterfall over on the playsuite /games/ repo, specifically what sibling instance called waterfall did for voice b. when you make this voice, add those effects." Then: "its the 1/n² bridge force".
- **OBSERVED — the Canopy Falls harp** (`~/Developer/playsuite-games-falls`, `src/music/styles/heavenly_harp.cpp`, commit 912b87f; the "waterfall" session's transcript, 2026-10-18):
  - Tuning B (D4–B5) beat A, but the owner found it "brighter, so its more audible, but its also colder". The session traced the cold to the output: bridge force, with harmonics falling as 1/n. W2: 85% the strings' own displacement (1/n², round and warm) and 15% bridge force, plus a fuller soundboard (lowest resonance 200 Hz instead of 230). Owner: "w2 is nice. really nice. that's a commit."
  - Its other parts: Karplus-Strong waveguides with an all-pass for the fractional period; two planes 0.022–0.045% apart, the second ringing 0.6× as long; a finger-pad pluck at 0.40–0.50 of the length (harmonics to about the 8th); a shared passive bridge (coupling 0.0005) for sympathetic ringing; a 1.5 ms slip-off onset; soundboard band-passes 200/480/1050 Hz (widths 110/160/420 Hz, gains 0.9/0.45/0.14); an 8-line FDN hall, 3.8 s; a soft saturator x/(1 + 0.3|x|).
- **`harp.warmth`, `harp.warmth_corner`** (BodyLayout → VoiceBody::step; default 0, gated):
  - Sound = (1 − w) · d(flow)/dt + w · 2πf · (flow less its steady stream). f is the intended pitch followed over 20 ms, so the two match at the fundamental. The steady stream is removed by two one-pole high-passes at 50 Hz.
  - V1 (`voice_resonant_gentle`), held ah / ee / oo, harmonics 1–24: at warmth 0 / 0.5 / 0.85 / 1, slope ah −7.6 / −8.2 / −10.2 / −13.4 dB/oct, ee −5.3 / −5.9 / −7.9 / −11.1, oo −13.0 / −13.6 / −15.6 / −18.8. H1–H2 ah 7.9 → 12.7 at 0.85. Banana energy below 60 Hz −34.8 → −27.6 dB of the total.
  - Sent `listen/harp_warmth/sustained_w0_before.wav`, `sustained_w0.85_harp.wav`, `sustained_w1.wav`, `banana_w0_before.wav`, `banana_w0.85_harp.wav`, `banana_w1.wav`.
- Owner on the output warmth: "he sounds better but muffled. like he's speaking behind a wall. it muffles out the warbling. that bridge is only supposed to work on the vocal cord, which needs harmonics." Output warmth stays at 0.
- **MEASURED — the folds' own harmonic series** (V1, cord flow by WAVEDUMP over 1.0–1.5 s, harmonics 1–40):
  - Slope −14.9 / −15.9 / −15.8 dB/oct (ah / ee / oo).
  - Against 1/n² at n = 2, 4, 8, 16, 32: ah −2 +1 +3 −6 −9; ee 0 +1 −1 −1 −10; oo −1 +1 −1 −6 −8. On the line to the 8th harmonic, short above it.
  - Harmonics above the floor between them at n = 2, 4, 8, 16, 32: 57–59, 45–48, 32–40, 16–23, 9–15 dB. The upper series is buried.
- **`harp.bridge`** (BodyLayout → Cord::set_bridge; default 0, off path verified identical with cmp):
  - The flux the airways receive becomes flow + b · (d(flow)/dt / 2πf − swing). Swing is the flow less its steady stream (one pole, 50 Hz); f is the intended pitch over 20 ms; the rate of change passes two poles at 5 kHz.
  - Unbounded, it blew up at share 0.05 (held oo, 2.468 s, glottis open): its gain rose toward Nyquist and fed back through the tract's echo into the flow solve. With the 5 kHz poles, 0.05 and 0.1 are stable; 0.15 (the harp's W2) still blows up (held oo, 2.216 s).
  - Source against 1/n² (the bridged flow rebuilt offline from the dumped flow), n = 2, 4, 8, 16, 32:
    - 0.05: ah −2 +1 +4 −4 −3; ee 0 +1 0 −1 −5; oo −2 +2 0 −3 −4.
    - 0.1: ah −3 +1 +7 −2 +12; ee 0 +3 +3 +1 +5; oo −2 +2 +4 0 +11.
  - Sent `listen/harp_bridge/sustained_before.wav`, `sustained_bridge_0.05.wav`, `banana_before.wav`, `banana_bridge_0.05.wav`.
- Owner on the bridge: "terrible noise between words. bridge is broken, inaduble, no words. retract. implement more of vocallab math". **The bridge is retracted** (code, slot and member removed; the flow line is back to `flow + flux_noise_`). REJECTED: it amplified the flow's rate of change, breath noise included, between words.
- **VTL 2019 additions on the kinematic folds** (`Cord::set_kinematic_shape`; `kinematic.skew`, `double_pulsing`, `flutter`, `length_scaling`; all 0 → identical, verified with cmp):
  - The skew knee is our reading of the text: the opening takes π(1 + s) of each 2π. Area pulses checked: s 0.5 opens slowly and shuts fast; d 0.5 makes every second pulse 0.16 against 0.25 cm².
  - `data/voice_vtl2019_male.params` (their most natural male on V1's body, KINEMATIC=1): F0 114.5, jitter 1.06%, shimmer 3–5%. Sent `listen/vtl2019/`.
- Owner: "the provided samples from their 2019 model do not sound natural. only their 2.4 version sounds natural, so if they updated that we will need to figure it out."
- **OBSERVED — VTL 2.4's "Geometric glottis 2025"** (`~/Downloads/VTL2.4/API/Developer/Sources/GeometricGlottis2025.cpp` on the laptop; the speakers JD2 and M01 select it; GPL, math only):
  - L = L0 √(f0 / f_rest), T = T0 L0 / L (L0 1.6, T0 0.45 cm, f_rest 120 Hz).
  - Amplitude 0.0007 (L0 / 1.6) √P_trans. P_trans is the commanded lung pressure less the supraglottal pressure through a 4th-order Chebyshev at 25 Hz. Faded to 0 between 2000 and 1000 dPa, times the relative amplitude. The 2019 L0 / L scaling is commented out.
  - Edges: cosine phases, the upper lagging 1.22 rad.
  - Collision: A2 = −rest / π, A3 = A2 − 0.2 A. t2 = acos(A2 / A)/2π and t3 = 1 − t2. Pulse shape p moves t4 = t3 + 0.1p (p > 0) or t1 = t2 + 0.1p (p < 0). Between t1 and t4 the displacement is the largest of the cosine and two exponentials toward A3, slope-matched at t1 and t4.
  - Area: 32 points along the fold, x = rest · z + d · sin(πz) cut at 0, summed as 2∫x plus the chink.
  - Flutter: Klatt, 25%. Jitter: Schoentgen model I, phase += 2π · 0.28 · (jitter / 100) · (±1) · √dt; at 10% that is 0.176 rad/√s.
  - The glottis is two tube sections, T/2 each.
  - Modal shape: JD2 rests 0.0403 / 0.03995 cm, chink 0.02 cm², pulse shape 0.5, 800 Pa; M01 0.054 / 0.0718, chink 0, pulse shape 0.502.
- **`kinematic.geometric_2025`, `pulse_shape`, `relative_amplitude`** (`Cord::set_geometric_2025`; off → kinematic identical, verified with cmp):
  - Our code from that math. The supraglottal low-pass is two Butterworth sections (theirs Chebyshev). The commanded pressure is our lung pressure (`expected_breath_`).
  - `data/voice_vtl2025_jd2.params` (JD2 modal on V1's body, flutter 25, jitter 0.176): area 0.044–0.188 cm², flow 128–447 cm³/s, smooth closure; held ah F0 114, jitter 0.53%, shimmer 4.0%.
  - Sent `listen/vtl2025/sustained_vtl2025_jd2.wav`, `banana_vtl2025_jd2.wav`.
- Owner: "implement their glottal flow and aspiration too".
- **OBSERVED — VTL 2.4's glottal flow** (TdsModel.cpp prepareTimeStep; Constants.h):
  - Each section has L = ρ (len / 2) / A per half and R per half = (len / 2) · max(0.03 · A^−2.6, 0.0065 · A^−1 · √2250). Areas are floored at 1e-4. The glottis sections have no soft walls.
  - Kinetic terms are linearized on the previous flow: k = 1 into a contraction, 0.2 out of a glottal section expanding, 0.5 elsewhere.
  - The glottis is two sections of T/2. Their rate is 48 kHz, Crank–Nicolson (θ 0.5).
- **OBSERVED — VTL 2.4's aspiration** (calcNoiseSources, calcNoiseSample):
  - The constriction is the upper glottis section. Flow is its outgoing current through one pole at 500 Hz, outflow only; v = flow / A.
  - Gain 0.2e-6, f_c 1000 Hz, tapered between 0.5 and 1 cm². Amplitude gain · |v|v² · √A.
  - The rise is a 3 ms one-pole on the square root and the fall is instant; the source is silent under 0.001.
  - Noise: N(0, 1/√12) cut at ±1, through a 2nd-order low-pass (Q 1/√2) at f_c, times amplitude · (1000 / f_c)^1.5 · 0.6.
  - A dipole (series pressure source) at 1.5 cm above the glottis, split linearly between the section there and the next.
- **`vtl.glottal_flow`, `vtl.inertance_floor`** (Cord::set_vtl_flow; off → identical, verified with cmp for the kinematic, fold-tube and 2025 sources):
  - Their viscous law per edge section and their inertance (no inlet or outlet, floor 1e-4). Our kinetic terms already were theirs. Our quadratic solve is kept (they linearize).
- **`vtl.aspiration`, `vtl.aspiration_gain`** (VoiceBody::vtl_aspiration_sample):
  - Their source as above, in place of ours, the flux noise and the fold turbulence.
  - Bilinear Butterworth at 1 kHz. Noise raised √2 for our 96 kHz.
  - It is about 12× quieter than V1's aspiration (gain 12.12 on the same normalization).
- **MEASURED** (2025 JD2 on V1's body, held ah 1–2 s):
  - 2025 alone: jitter 0.53%, shimmer 4.0%, fuzz −14.7 dB.
  - + their flow: flow 110–433 cm³/s (min was 128), jitter 0.49%, shimmer 3.6%, fuzz −16.1 dB.
  - + their aspiration: jitter 0.43%, shimmer 1.9%, fuzz −19.2 dB, 2–5 kHz flutter 1.49 → 0.71 dB.
  - Banana clicks (clicks2): 2025 5, + flow 7, + aspiration 23 (V1 29). The added clicks spread through the voiced stretches: the amplitude follows the 500 Hz flow over the instantaneous area and falls at once, so it jumps within each cycle. That is their formula.
  - `data/voice_vtl2025_flow.params` (flow; our aspiration still on at V1's gain) and `voice_vtl2025_full.params` (flow and their aspiration).
  - Sent `listen/vtl2025/sustained|banana_vtl2025_flow.wav` and `sustained|banana_vtl2025_flow_aspiration.wav`.
- Owner: "lets finetune what can be finetuned against our human reference recordings for pure vocal cord sounds? and lets finetune against taking the angle on the specgram". Then: "For 'the angle on the spectrogram' i meant and now do taking np.angle ie the phase. if you want a different phase representation you can take irfft twice in time which unwinds the representation and places everything in the real domain."
- **`analysis/cord_features.py`**:
  - tilt (0.3–5 kHz, dB/oct), H1–H2, harmonic peaks over troughs in 0–1 / 1–2 / 2–4 / 4–6 kHz;
  - pitch drift / tremor / flutter (cents, 0.5–3 / 3–8 / 8–20 Hz) and level wobble;
  - **phase from np.angle**: per 5 ms frame, three periods (Hann), each harmonic's np.angle at k·f0; phase distortion PD_k = φ_(k+1) − φ_k − φ_1 (Degottex); per harmonic over time the deviation √(−2 ln R) (scatter, rad) and the mean angle, medians per band.
  - The double-irfft representation is not built yet.
- **Mistake, corrected by the owner**: "my recordings are not aaa eee or ooh". `voice/ref/user_three.wav` is 21.75 s of connected speech, "hello my name is ducklips, the fox jumped over the laazy dog. she sells seashells by the seashore", then the same whispered (the normal pass is 0–9.8 s). I had treated its voiced runs as held vowels. It also carries steady 400 and 800 Hz lines through the silences (hum).
- **The held-vowel fit is withdrawn** (`analysis/cord_fit.py`; `data/voice_vtl2025_fit.params` cost 0.324, `fit2` 0.175). Its owner-derived targets (drift, tremor, wobble) came from speech. Its glottis went fully adducted (rests 0) and the banana clicked 107 times (25.7/s). Not from the aspiration or the inertance floor: gain 1, gain 1e-4 and floor 0.01 all left 108–137.
- **MEASURED, same words** (voiced runs of at least 0.12 s, medians):
  - Owner's normal pass (14 runs): tilt −11.1, H1–H2 −0.4, HNR 17.5 / 8.8 / 7.3 / 6.8, phase scatter 0.88 / 1.16 / 1.09, PD mean +0.8 to +1.0, flutter 9.1, clicks 1.4/s.
  - V1 saying the same sentences: HNR 29.9 / 20.2 / 11.5 / 6.8, scatter 0.58 / 0.74 / 0.70, PD mean +2.6 to +2.9.
  - VTL 2025: HNR 29.1 / 18.8 / 9.7 / 6.2, scatter 0.29 / 0.49 / 0.62, PD mean −2.6 to −2.7.
  - bdl's vowel nuclei: PD mean +2.6 to +3.0. The PD mean differs between owner and bdl (speaker or recording), so it is not a target.
- Owner: "id rather you fit to the nato voice".
- **MEASURED, NATO words** (the reader's 26 clips against `say` of the 26 words, `/tmp/es/ref/nato.score`):
  - Reader (28 runs): tilt −13.4, H1–H2 −9.6, HNR 26.4 / 12.9 / 8.6 / 6.0, scatter 0.89 / 1.43 / 1.59, flutter 12.8, clicks 0.2/s.
  - V1: H1–H2 +6.0, HNR 22.8 / 12.9 / 7.9 / 6.8, scatter 0.50 / 0.65 / 0.77, flutter 8.1, clicks 5.0/s.
  - VTL 2025: H1–H2 +3.4, HNR 19.8 / 9.1 / 7.3 / 6.3, scatter 0.37 / 0.57 / 0.74, flutter 7.2, clicks 4.9/s.
  - Both humans' phase scatter is about twice ours in every band.
- **`analysis/cord_fit_sentences.py`**:
  - Fits the 13 settings to a human reading of the same words: tilt, H1–H2, HNR per band, phase scatter per band and flutter, plus a squared penalty for clicks above the reference's rate.
  - Workers get the base through a pool initializer. macOS spawns them fresh; the first held-vowel run rendered over engine defaults until fixed.
  - Running against NATO (30 generations, from `voice_vtl2025_full`): start cost 12.56.
- Owner: "double-irfft means instead of specgram which is abs of fft, take the irfft, then slice the top half. you wind up with a phase unrolled representation."
- **`analysis/real_spectrogram.py`**: per Hann frame (1024 samples, hop 128), numpy.fft.irfft of the frame, top half (Nyquist → 0, flipped). That is the frame's cosine transform: magnitude × cos(phase), signed. Shown on a signed log scale.
  - `listen/reference/real_alfa_bravo.png`: the NATO reader against VTL 2025 and V1 on "alfa" and "bravo".
  - OBSERVED by eye: the reader's harmonics bend with the falling pitch, with tilted fringes; ours run flat, with near-vertical fringes. Ours is a dense saturated mass under 1 kHz, the reader's troughs stay clean. The reader's words are about 0.4 s against our 0.65 (rate 0.75). Not yet a measure.
- **NATO fit** (`data/voice_vtl2025_nato.params`, 30 generations, cost 12.56 → 3.13):
  - Settings: rests 0.0086 / 0.0041 cm, chink 0, lag 0.90 rad, pulse shape −0.49, amplitude 0.00038, flutter 30, jitter 0.56, aspiration gain 13.4, variation breath 0.010, tension 0.019, tremor 0.0025, breath tremor 0.11.
  - Features against the reader: scatter 0.61 / 1.13 / 1.31 (0.89 / 1.43 / 1.59); H1–H2 −2.5 (−9.6); HNR 17.5 / 8.5 / 7.0 / 6.3 (26.4 / 12.9 / 8.6 / 6.0); flutter 10.8 (12.8); tilt −13.2 (−13.4); clicks 1.8/s (0.2). Banana: 13 clicks (3.1/s).
  - HYPOTHESIS: the reader's harmonics are clean under 1 kHz yet their phases scatter. Ours buy scatter only with noise (phase jitter, 13× aspiration), which fouls the low band: HNR 0–1 kHz fell from 19.8 to 17.5 against a target of 26.4. The largest remaining terms are HNR 0–1 kHz 8.7, H1–H2 8.2 and scatter 3.5–4.1.
  - Sent `listen/nato_fit/nato_words_fitted.wav`, `nato_words_before.wav`, `banana_fitted.wav`, `banana_before.wav`.
- Owner on the NATO fit: "much more human. some clicks, tho."
- **Click hunt** (fitted voice; `clicks2` K 8: NATO 51 = 1.8/s, banana 13 = 3.1/s; reader 0.2/s):
  - REJECTED as causes (counts unchanged or higher): VTL aspiration (gain 1e-4 or 1), our glottal flow instead of theirs, no phase jitter, no body variation, no oral frication, no skin, jet attack 0.02.
  - Pulse shape matters for the K 8 count: 0 → 6.1/s, −0.49 (the fit) → 1.8, −1.0 → 1.6. Phase lag 1.3–1.8 → 2.2–2.3.
  - These K 8 events are glottal closures ringing the formants above 5 kHz wherever little else sounds there (nasals, voicing onsets).
  - REJECTED: soft contact in the 2025 glottis (`kinematic.contact_softness_2025`, 0.0002–0.004 cm). Counts 1.8–2.1/s NATO, 3.4–4.1 banana; features unchanged.
  - Strong clicks (K 16, `/tmp/es/strong.py`): NATO 11.026, 11.053 (juliett's /t/), 11.128, 14.299 (pauses), 22.323 (/uw/ → /n/), 22.919 (/m/ → /v/); banana 1.890 (/ah/ → /n/). The same instants (to 0.1 ms and 0.1 dB) in V1, VTL 2025 and the fit, so they are not from the folds. Each is a one-sample spike with brief ringing.
  - Three coincide within 1–2 ms with the lip end moving a node (ENDLOG). Two follow score commands within a millisecond: 14.2986 `glottis`, 22.9187 `lungs 0`. The commands drive second-order motors, so nothing steps.
  - REJECTED for the strong clicks:
    - logarithmic support glide (`walls.support_log_glide` 5–20 ms; the outside-the-lips node's support jumps ×1000);
    - VTL's area-rate limit on the vocal tube and nose (`vtl.area_limit`: 400 cm²/s, and the t⁴ curve, c 480e6 opening and 40e6 closing; TdsModel::setTube);
    - the room (ROOM=0).
  - All new switches default off and leave the outputs identical (verified with cmp).
  - Sent `listen/nato_fit/click_candidates.wav` (the seven strong candidates, 0.6 s each) to learn which clicks the owner hears.
- Owner on the candidates: "5 clicks, the rest dont seem to actually be clicking". Candidate 5 is NATO 22.323 s, "uniform" /uw/ → /n/.
- **FOUND — click 5**:
  - Every-step pressures (WAVEDUMP) put its source at node 20–21 (7.5 cm), at 22.322875 s, spreading one node per step both ways.
  - The node's link flow jumped from 12.4 to 41.0 cm³/s in one step while a geometry update arrived. The velum's side port follows the articulation (`set_velum_place` → `move_branch_group`, the nearest node) and stepped from node 21 to 20 as the velum lowered for /n/. The nasal flow (about 25 cm³/s) landed on another node at once.
  - It is in both VTL voices but not V1: the kinematic folds drive more flow through the nose there.
  - REJECTED on the way: the trace's "narrowest" jump at 22.3243 (the lip end moving, a reporting artefact), the wall motion kink (STEPS sums the walls' motion before its two 250 Hz poles), and the under-tongue port waking (it opened smoothly at 22.3112).
- **`branches.glide`** (BodyLayout::branch_glide; 0 = the port jumps, as before; outputs identical, verified with cmp for the NATO fit and V1):
  - A moving side port on the vocal tube (velum, under the tongue, cheeks) keeps a shadow port at the node it left, holding a share of the area that fades with this time constant.
  - The two act as parallel ducts. The branch hears their waves weighted by area.
  - Strong clicks (K 16): fit 6 NATO, 1 banana → glide 3 ms: 0 and 0; 6 ms: 0 and 0; 12 ms: 1 ("yankee" 26.903) and 0. K 8 counts: NATO 51 → 49, banana 13 → 11 at 6 ms.
  - NATO features unchanged: H1–H2 −2.6, HNR 17.8 / 8.4 / 6.9 / 6.6, scatter 0.63 / 1.16 / 1.31, flutter 10.8.
  - `data/voice_vtl2025_nato_glide.params` (glide 6 ms). Sent `listen/nato_fit/nato_words_fitted_glide.wav`, `banana_fitted_glide.wav`.
- `speak PORTLOG=from,to`: every step's velum and under-tongue port areas (`VoiceBody::under_tongue_port()` added).
- Owner on the glide takes: "charlie, foxtrot, xray have harsh sounds. alpha, golf, mike, qeuebec, yankee sound wrong or are essentially misshapen. the rest are closer to reasonable than they have been before."
- **MEASURED — harsh frication** (`/tmp/es/fric.py`: unvoiced frames with 2–10 kHz energy within 20 dB of the word's highest; level re the word's loudest frame; third octaves 0.5–10 kHz):
  - Reader: charlie −15.4 dB, peak 3.15 kHz; foxtrot −13.3, peak 5 kHz; x-ray −11.7, peak 6.3 kHz. All fall off by 8–10 kHz.
  - Ours (fit + glide): −6.9 / −4.7 / −4.7 dB, peaking at 10 / 10 / 5–6.3 kHz, 8–10 kHz near the top. V1 is −9.1 / −7.2 / −7.8.
  - The VTL aspiration gain (1, 4, 13) changes nothing here; the oral frication does it.
  - Sweep over `oral_noise.frication` and `jet_corner_ceiling`. At 0.5 / 6000: −11.4 / −14.5 / −13.5 dB; charlie peaks at 3.15 kHz; foxtrot still flat to 10 kHz (its /f/ seems not to obey the ceiling; OPEN).
  - `data/voice_vtl2025_nato_fric.params`. Strong clicks: 2, both in pauses after a final /t/ (6.669 and 11.128 s; the second is candidate 3, which the owner did not hear as a click).
  - Sent `listen/nato_fit/nato_words_frication_fixed.wav`, `banana_frication_fixed.wav`.
- **MEASURED — misshapen words** (`/tmp/es/wordform.py`: LPC F1–F3, voiced frames, medians per fifth of the voiced stretch):
  - Our open vowels are too closed. F1, reader against ours: alfa ~1000 vs 650; golf 740–875 vs 450–600; delta ~840 vs 390–560; bravo 700–870 vs 580–700.
  - The postures were fitted to bdl; the NATO reader opens wider.
  - The reader's "mike" and "quebec" gave no LPC rows (OPEN: voicing or level gate).
- **NATO fine-tune relaunched** (/tmp/es/nato2): `policy_fit.py als 2`, ORDER=alternate, GROUPS organs,plan,tongue, ENGINE_TONGUE_KEYS "tongue.skeleton_,tongue.tip_radius", ENGINE_RINGING=0, from `voice_vtl2025_nato_fric.params` and organs_skeleton.
  - New switch `ENGINE_SOURCE=kinematic` in policy_fit.py and engine_net.py renders with KINEMATIC=1 instead of FOLDTUBE=1.
  - 267 neurons in 20 blocks; start cost 11.97 (the fold-tube run started at 12.65).
- Owner on the frication fix: "hasnt fixed them. lets show the tongue behavior."
- Sent `listen/video/problem_words_slow2.mp4` (alfa, golf, mike, quebec, yankee, charlie, foxtrot, x-ray at half speed; `voice_vtl2025_nato_fric.params`) and `problem_words_frames.png` (stills at the key vowels).
- OBSERVED in the stills:
  - The jaw barely moves between /ae/, /eh/, /y/ and /ng/.
  - For /ae/ the tongue body fills the oral cavity under a nearly closed jaw.
  - /ao/ (golf) holds a 0.2 cm² pharyngeal constriction over about 4.5 cm.
  - organs_skeleton jaw angles: EH 1.05, IY 1.95, EY 3.2, UW 4.3, AH 5.0, AE 6.5, AO 10.7, AA 12.8. The bdl formant fit reached /ae/ (F1 637) and /eh/ (595) with an almost shut jaw and tongue shape. The NATO reader's /ae/ F1 is about 1000.
  - "quebec" is planned K EH B EH K (no /w/), which matches the ICAO "keh-BECK".
- Owner on refitting the vowels with a working jaw: "yes, do that."
- NATO fine-tune (/tmp/es/nato2) stopped before it finished its first block; it started from the shut-jaw postures.
- **`organ_vopt3.py` jaw minimums** (`JAW_MIN`, degrees): aa 14, ae 13, ao 11, ah 9, eh 8, ey 6, ow 6, ih 5, uh 5, ax 4, iy 3, uw 2. The jaw's ceiling is raised from 14 to 20. The uw/ow rules (jaw ≤ 8, lips rounded, body forward) and the tip-behind-the-incisors rule stay.
  - Running: START = CURRENT = organs_skeleton, W 10, MAXFEV 260, SPEAK_ENV KINEMATIC=1 with `voice_vtl2025_nato_fric.params`, out `data/organs_jaw.txt`. ER keeps its posture (no bdl target).
- **Jaw refit** (`data/organs_jaw.txt`, bdl targets): jaws open (AA 14.2, AE 13.05, AO 14.9, AH 9.2, EH 8.2, IH 9.8, OW 8.0, UH 8.4, UW 4.7, AX 6.1, IY 3.0). Formants within 0.35–2.3 st of bdl's (OW 3.7).
  - It changes almost nothing in sound. The fit raised the tongue to hold bdl's F1 (AE ty 7.32 → 7.71, F1 632): the tube is almost as before and the tongue still fills the mouth (stills `/tmp/es/ref/frames_jaw_sheet.png`).
  - The target, not the jaw, held the vowels closed: bdl's formants are connected speech; the reader speaks stressed single words.
- **The reader's stressed vowels** (`/tmp/es/readervowels.py`: our score's phone timing mapped proportionally onto each clip's speech span, LPC over the middle half; 1–7 frames a word, noisy, e.g. charlie 1284 / 2554 / 4287):
  - Medians: aa 887 / 1457 / 3044, eh 881 / 1930 / 2900, ih 459 / 2707 / 3138, iy 472 / 2831 / 3303, uw 496 / 2202 / 2884 (fronted), ow 656 / 1413 / 2731, ae 605 / 2668 / 3372 (only "tango" and "yankee", before nasals, raised).
  - Scaled by 1 / 1.12 for his shorter tract (the F3s) into `data/reader_vowel_targets.json` (aa, eh, ih, iy, uw, ow; ae, ao, ah, uh, ax, ey stay bdl). `organ_vopt3.py` takes `TARGETS=file.json`.
- **Lexicon** (voice/src/lexicon.cpp): alfa AE1 → AA1 L F AH0 (the reader's "AHL-fah", F1 990); golf G AO1 → G AA1 L F (as CMUdict).
- **Reader-vowel refit** (`data/organs_reader.txt`, W 10, from organs_jaw, then `organs_reader2.txt`: eh and iy again at W 4, MAXFEV 400):
  - aa jaw 14.2, F1 763; eh jaw 11.7, F 708 / 1621 / 2501 (target 787 / 1723 / 2589); ih 385 / 2068; uw 443 / 1940 (jaw 7.8, fronted); ow 569 / 1264.
  - iy did not move at either weight (F2 2078 against 2528).
- **MEASURED (NATO words, `/tmp/es/ref/nato3.score`)**: F1 alfa 723–805 (was 546–660; reader 723–1001), golf 757 (was 449–609; reader 736–875), quebec /eh/ 600–660 (was 406–579), delta partly.
  - Strong clicks 4: 6.675 (pause after foxtrot), 10.461 (juliett /jh/), 18.959 (romeo /ow/ → /m/), 22.819 (uniform /r/ → /m/, −10.5 dB re peak). The last three are new closures; HYPOTHESIS: the wider jaw makes the lips travel and seal harder.
  - Sent `listen/nato_fit/nato_words_reader_vowels.wav`, `banana_reader_vowels.wav`, `listen/video/problem_words_reader_slow2.mp4`.
- Owner on the reader-vowel takes: "most are wrong. alpha is down up lips. golf is up down lips. mike is lips back front upper. charlie is up down back up down. foxtrot is lips up down up. xray is up down up-ish. xray and alpha are kind of kay." (The articulators' sequence each word needs; the video showed them missing.)
- **MEASURED — gestures** (`/tmp/es/gestures.py`: AREALOG every ms; per phone, the tightest node past 4 cm over the phone and when):
  - With organs_reader2, the velars no longer seal: mike /k/ 0.11 cm² (old postures 0.00), quebec 0.03–0.09 (0.00), x-ray 0.11 (0.00), yankee /ng/ 0.21 (0.01), golf /g/ 0.08 (0.05), foxtrot /k/ 0.13 (0.10).
  - alfa's /l/ was narrowest in the pharynx (4.6 cm). charlie's /r/ and /l/ made no oral constriction, with the old postures as well.
- **FOUND in the planner** (voice/src/planner.cpp):
  1. /k g ng/ send only a body target; the jaw stays at the vowel's (now 13–14 degrees) and the body aimed from it cannot reach the roof.
  2. /l/ and /r/ send no jaw command either.
  3. /r/ hands the tip back to the vowel at finish − 0.03 s, 10 ms after charlie's /l/ aims it (shape_time = start − lead), and undoes it.
  - New ConsonantShapes (defaults leave the plans identical; verified with cmp on the problem-word score): `articulation.velar_jaw`, `lateral_jaw`, `rhotic_jaw` (default 90, no ceiling), `articulation.liquid_yields_tip` (default 0). At 1, /l/ and /r/ skip the tip hand-back when the next phone aims the tip.
  - `voice/data/plan_jawcaps.params`: velar_jaw 6, lateral_jaw 6, rhotic_jaw 8, liquid_yields_tip 1.
- **Tongue reach**: TONGUELOG through alfa's /l/: the tip target (6.67, 10.61) in the jaw frame, the tip reaching (6.12, 9.72), goal_tip_miss 0.7–0.9 cm. `tongue.skeleton_stretch` 0.5 → 0.8 lets it reach; 1.2 changes nothing more. `data/voice_vtl2025_nato_stretch.params`.
- Velar timing: dorsum_close_move 0.06 → 0.04 and seal_lead 0.4 → 0.8 changed the score (moves 40 ms, 8 ms earlier) but not when the tongue arrives. The skeleton tongue's own dynamics set the pace (OPEN).
- **MEASURED (plan_jawcaps + stretch 0.8)**:
  - Tightest per phone: alfa /l/ 14.4 cm, 0.093; golf /g/ 0.011; mike /k/ 0.014 at 90%; quebec /k/ 0.011 at 82%; foxtrot /k/ 0.002 at 77%; x-ray /k/ 0.000 at 87%; charlie /l/ 13.7 cm, 0.079, /r/ 0.109.
  - yankee /ng/ 0.114 at the phone's end, unchanged (OPEN).
  - NATO words strong clicks 1 (x-ray /r/ 22.241), K 8 42.
  - Sent `listen/video/problem_words_gestures_slow2.mp4`, `listen/nato_fit/nato_words_gestures.wav`.
- Owner: "this need a lot more work. your base is basically.. a bad tongue. real tongues dont look like this." Asked how to rebuild, the owner chose "VTL-style tongue only" (keep our hull and jaw).
- **OBSERVED — VTL 2.4's tongue contour** (VocalTract::calcTongue 3806–3990; Geometry.cpp getEllipseTangent, getCircleTangent; GPL, math only):
  - The root is a rational quadratic Bezier (weights 1, 2, 1) from the hyoid through (TRX, TRY) to its tangent on the body ellipse (JD2: 1.8 × 1.8 cm).
  - Then the ellipse's arc (angles falling from the root's tangent to the blade's). The blade is a Bezier from the ellipse through (TBX, TBY) to its tangent on the tip circle (r 0.2). Then the tip arc to angle 0.
  - Automatic root: TRX = 0.938 |C − H| − 5.11, TRY = 0.831 TCX − 3.03 (constants in VTL's own frame).
  - JD2's 80 shapes put the blade control point at a median 0.61 of the way from body centre to tip (0.39–0.94), 1.15 cm above that line (0.35–2.15).
- **`tongue.vtl`** (+ vtl_body_radius 1.8, vtl_tip_radius 0.2, vtl_blade_along 0.61, vtl_blade_lift 1.15, vtl_root_bulge 0.3; `MidsagittalTract::vtl_tongue_outline`, VtlTongueShape):
  - With tongue.physical 0 the speaker hands this contour to the geometry (jaw frame, 24 points root to tip). The body ellipse is centred on the body's position; the tip circle's front lands on the planned tip.
  - Root control point: midway from our root base to the ellipse, 0.3 cm behind that line. Blade control point by the JD2 medians.
  - Off by default (the current voice is identical, verified with cmp). `geometry_probe` takes `anat.vtl=1`; `organ_vopt3.py` passes `PROBE_EXTRA`.
  - `data/voice_vtl_tongue.params` (stretch voice, physical 0, vtl 1).
  - OBSERVED with organs_reader2 (fitted for the old tongue): a smooth dome on a root, a sloping blade, a small tip. /g/ and /k/ press the dome on the roof (area 0 at 9–11 cm). /l/ raises the tip to the ridge.
  - Sent `listen/video/vtl_tongue_frames.png`, `vtl_tongue_slow2.mp4`.
  - Vowel refit for it running: reader targets, jaw minimums, W 6, MAXFEV 300 → `data/organs_vtl_tongue.txt`.
- **New-tongue vowel refit** (`data/organs_vtl_tongue.txt`: PROBE_EXTRA anat.vtl=1, reader targets, jaw minimums, W 6, MAXFEV 300). Formants within: aa 0.68 st, ah 0.12, ax 0.37, ow 0.49, ae 1.05, uw 1.08, uh 1.17, eh 1.77, ao 2.43, ey 3.19, ih 3.27, iy 3.61.
- Owner (2026-10-18): "we now have multiple measures- phase specgrams, wavelets, other transforms- and ample reference material- and ample design material. set a goal for yourself and seek to work carefully until you get to 50 db. timing and behavior may have to be adjusted ... if you need to download some more reference audio or another small model or two, that's fine ... don't ask me to check audio- use technical measures only ... use profiles and then to refine the profiles ... so that you can fit genders and depths ... go far and go long."
- **`analysis/fidelity.py`**: fidelity of a reading to a reference of the same text, in dB.
  - Phones from labels (ARCTIC), from our score, or carried by DTW on MFCCs (NATO).
  - Groups per phone class: spectrum (24 ERB bands, middle 60%, dB re the utterance's speech level); dynamics (std of six octave-band levels within a token); phase (np.angle phase-distortion deviation in 3 bands, at frame-wise pitch, pooled over front / central / back vowels); pitch (F0 percentiles 5 / 25 / 50 / 75 / 95, semitones); timing (mean ln duration).
  - Per group: 10 log10(spread of the reference's features / our squared error). Combined: 10 log10(groups / Σ error / spread).
  - Per-vowel pitch and per-vowel phase were tried first and failed calibration (bdl against himself −6.1 and 4.0 dB).
- **Calibration** (a0001–a0030; the first row against a0031–a0060):

  | | combined | spectrum | dynamics | phase | pitch | timing |
  |---|---|---|---|---|---|---|
  | bdl vs bdl (other sentences) | 8.5 | 13.2 | 10.1 | 4.5 | 20.4 | 7.2 |
  | bdl vs rms (same sentences) | −0.4 | 5.7 | 3.9 | −2.5 | −4.6 | 7.6 |
  | bdl vs slt | −5.0 | 7.5 | 6.0 | −5.8 | −10.6 | 10.2 |
  | ours (VTL tongue) vs bdl | −3.2 | 4.5 | 1.9 | −9.6 | 11.7 | 8.5 |
  | ours (skeleton) vs bdl, fixed-pitch phase | −2.4 | 4.4 | 2.1 | −8.7 | 11.8 | 8.5 |

  - bdl a0001–a0100 against a0101–a0200 (fixed-pitch phase version): combined 11.8, spectrum 18.6. The same-speaker ceiling grows with data.
  - 50 dB is far past what these statistics resolve. Working goals per profile: above the different-speaker floor (−0.4 dB), then toward the same-speaker figure (8.5–12 dB), then as far as the data allows.
- `analysis/render_set.py` (parallel ARCTIC renders) and `analysis/profile_fit.py` (diagonal ES over 24 tree leaves in groups source, noise, variation, tract, tongue; train and held-out ARCTIC ids; reference features cached).
- Downloaded slt (CMU ARCTIC female, festvox `cmu_arctic/cmu_arctic/packed/cmu_us_slt_arctic-0.95-release.tar.bz2`, 120 MB) to `voice/ref/arctic_slt_dl/`, linked as `voice/ref/arctic_slt/{wav,lab}`.
- New slot `anatomy.scale` (BodyLayout::anatomy_scale → MidsagittalTract::set_scale; 1 leaves it alone). Not yet built on the Mini: the bdl fit runs on the current binary.
- bdl profile fit running (/tmp/es/profile/bdl): train a0001–a0020, held out a0021–a0030, 40 generations, from voice_vtl_tongue + organs_vtl_tongue + plan_jawcaps.
- **bdl tree fit, first run** (/tmp/es/profile/bdl): train −2.8 → 5.0 dB (spectrum 3.8, dynamics 1.6, phase 5.6, pitch 14.9, timing 8.3), held out −0.2 → 2.2 (spectrum 4.3, dynamics 1.7, phase −2.0, pitch 7.0, timing 7.1). Plateau from generation 1 to 5; stopped.
- **Phase is speaker-specific**: in ARCTIC sentences bdl's phase scatter is low (0.33–0.77 rad by band and vowel group). Ours sat higher at 1–4 kHz (0.95–1.11) and lower under 1 kHz. The NATO reader's is high (0.9–1.6).
  - Cycle-to-cycle variation (`kinematic.cycle_shape/lag/amplitude/rest`, new) made bdl's phase worse (−6.5 → −9.7 dB, 10 sentences). It is kept for scattered voices.
  - `vtl.aspiration_corner` (new; their filter corner and gain (1000 / f_c)^1.5; 1000 leaves it alone). At 3000 it changed little for bdl.
- **Vowel targets** (`analysis/vowel_targets.py`, f123.py's method: 10 kHz, order 12, a0001–a0060): bdl reproduces f123's table exactly; rms deeper (F3 about 2400); slt higher (F3 about 2950).
  - slt /iy/ (F2 670, a picking error) and /uh/ (too few tokens) are replaced by Peterson & Barney's women's 310 / 2790 / 3310 and 470 / 1160 / 2680.
  - `data/targets_bdl.json`, `targets_rms.json`, `targets_slt.json`. The first version at 16 kHz and order 12 had too few poles (bdl /aa/ F2 2825); discarded.
- F0 percentiles 5 / 50 / 95 (a0001–a0020): bdl 98 / 122 / 149, rms 74 / 95 / 123, slt 167 / 188 / 223 Hz.
- New slots (dev build; defaults leave outputs identical, verified with cmp): `kinematic.fold_length_2025`, `fold_thickness_2025`, `fold_rest_pitch_2025` (JD2 1.6 / 0.45 / 120; VTL W02 1.2 / 0.40 / 220); `anatomy.scale`.
  - `build_dev/` holds the newest speak while fits run on `build/`. `render_set.py` takes SPEAK_BIN; `organ_vopt3.py` takes VOWEL_PITCH.
- Profiles: `data/profile_rms.{params,plan}` (anatomy 1.05, pitch about 95 Hz); `data/profile_slt.{params,plan}` (anatomy 0.87, W02 folds and modal shape: rests 0.0515, chink 0.01, lag 0.52, pulse shape −0.2; pitch about 188 Hz).
- **`analysis/profile_cycle.sh`**: rounds of tree (source, noise, variation, tract, tongue, cycle) → vowels to the speaker's targets → plan (pitch, timing, consonants). Single-threaded workers, build_dev.
  - bdl cycle running (/tmp/es/profile/bdlc.log): from the first run's best, train a0001–a0016, held out a0017–a0024, 3 rounds × 8 generations.
- **Spectral error by class** (bdl first-run best, a0001–a0010; bias low / mid / high thirds of the ERB bands):
  - f 14% of the error, +27.7 / +16.4 / −5.3 dB; hh 8.7%, +17.9 / +16.9 / +11.1; jh 7.3%, high −22.6; y and w +15 to +18 mid; ow +16.5 mid; z, th high −15 to −16; r −7 to −11.
  - HYPOTHESIS 1, the folds swinging through voiceless sounds (the 2025 swing follows the commanded breath, which an open glottis does not lower). New `kinematic.abduction_fade_start/end` (the planner's gap; off by default, identical with cmp). Fade 0.04 → 0.10: spectrum 3.3 → 3.4, not the cause.
  - HYPOTHESIS 2 (open): the constrictions form late (alfa /f/ tightest only at 53%), so a consonant's middle still holds the vowel. `articulation.lead`, `tip_close_move`, `dorsum_close_move`, `timing.seal_lead` added to the plan_consonants group; the fade to the source group.
- **Dynamics** (a0001–a0010): ours moves 1.7 dB less within a token than bdl (4.56 against 6.28 dB). Most in the low bands (−2.5 to −3.0 at 125–250 Hz) and in consonants (f −5.7, jh −5.3, w −5.0, hh −4.9, sh −4.6): smooth where bdl switches sharply.
- Profile starts (8 sentences): rms start against rms 2.3 dB combined (the bdl voice against rms −1.0); slt start against slt 0.5 (pitch 11.4).
- `profile_fit.py` NATO mode (REFDIR "nato"): the 26 words as one render, the reader's clips cut by our phones through DTW each evaluation; even words train, odd words are held out.
- **`analysis/duration_calibrate.py`** (closed-form timing): plan the reference's sentences (say --plan-only), compare each phone's planned durations with the labels, scale `duration.X` by the geometric mean ratio, 3–4 passes. a0001–a0100, organs_vtl_tongue:
  - bdl: rms log ratio 0.133 → 0.001 over 38 phones (factors AA 0.67, AE 0.67, AO 0.70, UH 0.73, AH 0.85, Y 0.85, T 1.10, D 1.12).
  - rms: 0.226 → 0.001 (AA 0.56, AO 0.54, UH 0.52, HH 0.59, G 0.75).
  - slt: 0.211 → 0.001 (AA 0.60, AE 0.64, UH 0.55, Y 0.65).
  - At render_set's rate 0.75 the open vowels ran about 1.5× long.
  - `data/profile_{bdl,rms,slt}_durations.plan`.
- `analysis/profile_cycle2.sh`: profile_cycle.sh plus a duration stage after each plan fit (DURATION_IDS, default a0001–a0100). Kept separate because bash reads a running script as it goes.
- **bdl cycle round 1** (bdlc; train a0001–a0016, held out a0017–a0024):
  - Tree: train 4.4 → 5.46. Held out 2.5 → 2.0 at the best.
  - Vowels to bdl's targets: all within 0.01–2.25 st (ow 2.25, uw 1.58, iy 1.36).
  - Plan: train 4.43, held out 0.5 (phase −4.8).
  - Bound-pinned: kinematic.amplitude 0.0003 (floor), tract.viscous_scale 2 (ceiling), radiation_mass_scale 0.5 (floor). Widened to 0.0001, 5, 0.2.
- **Held-out noise floor**: bdl's training references (16) against his held-out ones (8) score combined 6.1 (spectrum 11.9, dynamics 8.1, phase 3.1, pitch 9.6, timing 3.8). With 8 held-out sentences the measure resolves about 6 dB.
  - Sets enlarged: train a0001–a0030, held out a0101–a0120 (bdld, profile_cycle2: 4 rounds × 8 generations, with the duration stage).
- **Slip, fixed**: widening the bounds put each comment before its tuple's comma. That is valid syntax (a tuple called as a function), so ast.parse passed and the fitter crashed at import. The cycle's stderr went to /dev/null, so the first bdld run silently skipped every fit stage. Check by importing, not parsing.
- **bdld round 1** (train a0001–a0030, held out a0101–a0120): the tree stage did not improve on the start (train 5.46, held out 3.9: spectrum 6.3, dynamics 4.3, phase −0.3, pitch 10.5, timing 5.7).
  - bdl's own ceiling on these sets (train references against held-out references): combined 7.3 (spectrum 12.6, dynamics 8.9, phase 3.3, pitch 16.7, timing 5.7). Our timing equals it.
- **Held-out error by class** (bdld start, a0101–a0120): the voiceless obstruents carry +8 to +16 dB of low bands.
  - f +15.5 low; th +16.0; voiceless stops +14.8; ch +11.4; jh +11.0; sh +8.4. Also high-band deficits: th −7.9, ch −10.9, jh −9.9, sh −11.1.
  - Spectra (dB re speech level, 144 / 317 Hz): bdl's s −41 / −44 against ours −30 / −36; voiceless stops −38 / −46 against −20 / −26; f −47 / −50 against −26 / −39.
  - REJECTED as the cause (10 held-out sentences): VTL aspiration gain 1 (no change); fold fade 0.03 → 0.08 (no change); ROOM=0 (low excess +18 → +14 only, everything else worse: combined 4.9 → −0.1, phase 1.5 → −5.4; the room stays).
  - Not mainly late closures either: constrictions first reached at 2–16% of the phone, held 54–98% (closure_timing.py over 6 sentences: CH 0.92, F 0.84, K 0.54, P 0.66, S 0.98, SH 0.98, T 0.83).
  - FOUND in one /t/ (a0101, 1.519 s): the seal came 28% into the phone; the glottis stayed at 0.03–0.16 cm² and kept vibrating through the closure. The planner's voiceless gap (0.14, a constant) opens it too little for the 2025 folds.
- `PlanOptions::voiceless_gap` (`glottis.voiceless_gap`, default 0.14, plans identical with cmp) in voice/build_dev. New plan_glottis group in profile_fit: voiceless, aspirate and voiced-fricative gaps, abduction lead and moves, adduction lead.
  - render_set / duration_calibrate take SAY_BIN. `profile_cycle3.sh`: dev say and speak, plan_glottis in the plan stage, stage errors kept.
- **REJECTED (10 held-out sentences, mean-phase measure; combined / phase)**: the start voice 4.9 / 1.5.
  - `motors.smoothing` 1, 2, 4: 2.8, 1.9, 0.7 (phase −2.2, −3.6, −5.3).
  - `walls.frame_corner` 60 Hz: −0.6 (phase −6.4). It stays at 250.
  - voiceless gap 0.14 alone, raised: no gain.
- **Vowel refit can degrade the whole voice.** organ_vopt3 against bdl's targets brought every vowel within 2.3 st, yet held-out combined fell 4.9 → 0.7. So profile_cycle3 now checks the vowel stage: it renders the held-out set before and after and keeps the refit only if the score holds.
- **Phase group, recalibrated** (OBSERVED; the reference against itself):
  - With the mean of the PD deviations, the phase group scored 3.3 dB for bdl train vs held out, but −1.3 for held out vs a second held-out set (a0121–a0140). It was noise, and it dominated the combined score.
  - With the median over frames: 2.4 and −0.4. Not more stable.
  - The cause is the signal term. The phase features are absolute scatters (0.3–0.8 rad); their spread about their own mean is tiny. Taking the signal as their sum of squares (fidelity_v3, deployed as fidelity.py on the Mini 2026-10-18 05:40 EDT; fidelity_v1.py keeps the old one):

    | pair (bdl reference sets) | spectrum | dynamics | phase | pitch | timing | combined |
    |---|---|---|---|---|---|---|
    | bdl30 vs bdl20 | 12.6 | 8.9 | 14.8 | 16.7 | 5.7 | 9.9 |
    | bdl20 vs bdl20b | 13.2 | 8.7 | 11.8 | 14.5 | 6.9 | 10.1 |
    | bdl20 vs rms20 | 5.3 | 3.8 | 0.7 | −5.6 | 6.5 | −0.3 |
    | bdl20 vs slt20 | 7.1 | 6.4 | 2.6 | −10.4 | 8.1 | −3.9 |
    | ours (start) vs bdl a0101–10 | 6.1 | 4.0 | 10.5 | 13.9 | 6.2 | 6.9 |

  - The same speaker now sits at about 10 dB, and another speaker at 0 dB or below. The ranking of the rejected variants holds (frame_corner 60: 4.5; smoothing 1, 2, 4: 6.5, 6.9, 6.6).
  - The goal is still 50 dB. At these set sizes the measure's same-speaker ceiling is about 10 dB, so beyond 10 the score stops distinguishing a better voice from bdl himself. More dB needs larger sets or a finer measure.
- bdle (profile_cycle3, 4 rounds × 8 generations, train a0001–a0030, held out a0101–a0120) restarted with fidelity_v3.
- **Phase pooled over all vowels** (fidelity v4, deployed 2026-10-18 06:05 EDT; `FIDELITY_PHASE_GROUPS=1` restores front/central/back).
  - Our own voice, rendered with four seeds (a0101–a0120), front/central/back: phase 5.0–11.2 dB, back-vowel scatter 0.41–0.93 rad. Too few tokens per group, so the fitter chased seed noise.
  - The running bdle cycle's held-out phase went −0.8, 2.9, −3.3 over near-identical candidates; even profile_fit's `%.10g` rewrite of the start voice moved it 11 dB.
  - Pooled (OBSERVED): bdl30 vs bdl20 phase 19.7 (combined 10.1); bdl20 vs bdl20b 26.7 (10.7); vs rms 0.9 (−0.3); vs slt 3.2 (−3.8). Ours over the 4 seeds: phase 11.3–16.0, combined 6.1–6.7.
  - The phase scatter in the 0–1 / 1–2 / 2–4 kHz bands (rad): bdl 0.37 / 0.49 / 0.57; rms 0.55 / 0.86 / 1.20; slt 0.39 / 0.73 / 1.09; ours 0.35–0.50 / 0.42–0.51 / 0.68–0.75 (ours is noisier at 2–4 kHz).
  - bdle stopped.
- **Low thumps in obstruents: the geometry, not the source** (a0101 "You, you would not keep the truth from me.", TRACE / AREALOG / PROFILE):
  - Inside sealed /t k/ (1.52–1.72 s) the 90–350 Hz band jumped to vowel level (−1 to +1 dB re speech) in single 10 ms frames, with nothing above 2 kHz. bdl's closures sit near −40.
  - Not the skin (the profile's skin_gain is 0.037), not voicing (periodicity 0.0–0.3), not the geometry frame (1 ms; the control tick is 0.125 ms).
  - CAUSE 1, the lip opening. `profile_to_nodes` gave the node the opening crosses the full area until its centre passed the opening, then half: 2.32 → 1.17 cm² within a frame. Its wall closed at −95, an air pulse to 153 cm³/s out of the lips (1.589 s). New `mouth.continuous_opening` (Speaker; default 0, output identical with cmp) blends by the share inside over the whole node: peak 59 cm³/s, that frame −1 → −16 dB.
  - CAUSE 2, cuts sliding between nodes. A node's area was the min of its centre and the cuts inside its span, so a narrow cut leaving the span jumped it. The front cavity's volume moved in 1 ms steps of −98, −136, +50 cm³/s while the organs moved smoothly. New `mouth.continuous_cuts` (default 0, identical) adds the profile at the span's two ends to the min (continuous in the cuts' positions; closures still never stepped over): the same window's volume rate runs −13 … +5 smoothly, some ±20 jitter left later in /k/.
  - The TRACE line's "narrowest 0.020" in every fricative was the sliver past the lips (TRACE scans two nodes past them). `analysis/constriction_probe.py` reads AREALOG up to the lips.
- **Narrowest areas per phone** (constriction_probe, a0101–a0120, bdl_start), median cm²: TH 0.000 and DH 0.000 (sealed: stops, not fricatives); IY 0.036; Y 0.022; W 0.063; F 0.082; V 0.089; S 0.060; Z 0.060; SH 0.046; R 0.150; L 0.120; UW 0.100; ER 0.187.
  - Typical values in the literature (approximate, HYPOTHESIS as targets): fricatives 0.05–0.2; glides about 0.15–0.3; /i/ about 0.2–0.3.
- **Variants against bdl a0101–a0120** (fidelity v4, two render seeds pooled), rms class error change in dB:

  | variant | spectrum | dynamics | phase | combined | moved |
  |---|---|---|---|---|---|
  | bdl_start | 5.3 | 3.8 | 14.0 | 6.5 | worst classes: f 14.1, iy 14.0, th 11.5, voiceless stops 11.1, y 10.7 |
  | + continuous_opening + continuous_cuts | 5.4 | 4.0 | 8.9 | 6.3 | voiceless stops −4.1, sh −3.6, th −3.4, ch −3.3, dh −2.7; m n w +1.6–2.2 |
  | + IY ty 8.282 → 8.12 | 5.9 | 4.9 | 9.8 | 6.8 | iy −7.3, y −5.2; hh +4.2 |
  | IY ty 8.06 | 5.8 | 4.6 | 9.6 | 6.6 | iy −8.5; hh +6.6 |
  | y_raise −0.1 | — | — | — | — | y −4.9 alone, nothing with IY lowered |
  | dental_press −0.05 / −0.10 (one seed) | 5.1 / 4.9 | | | 6.4 / 6.3 | dh +5.9 / +7.7: bdl's /dh/ is a near-closure |

  - The phase loss with continuous cuts is not added noise. The vowels became too regular: scatter at 0–1 / 1–2 / 2–4 kHz 0.435 / 0.476 / 0.721 → 0.125 / 0.325 / 0.596 rad (bdl 0.373 / 0.490 / 0.567). Part of our low-band irregularity was the geometry jolting the flow. The source's variation is to be refitted.
  - New start `data/profile_bdl_s2.*` (both switches, IY ty 8.12): held out 6.8 (spectrum 5.9, dynamics 5.1, phase 10.0, pitch 10.0, timing 5.7).
- **`analysis/class_fit.py`** (separable fit): the 13 vowel postures' organ fields and 26 consonant plan shapes (143 settings), each ranked by its own phone classes' cost.
  - A class's cost is its spectrum, dynamics and timing squared errors over each group's reference spread, so the classes' costs add up to the combined objective's error.
  - Antithetic pairs (10 renders per generation, a new seed each generation); the mean is scored on held-out every 3 generations.
  - Run bdl_c1 from profile_bdl_s2 (18 generations).
- **Where the class errors sit** (profile_bdl_s2, a0101–a0120, two seeds; ours − bdl, dB re speech level, 24 ERB bands 144 Hz – 6.7 kHz):
  - f: +10 to +18 in every band. A level error: bdl's /f/ sits 15–20 dB under his /s/, ours about level with our /s/. The tree's `oral_noise.frication` had been raised to 2.5 (from 0.5) to lift weak sibilants (still −3 to −8 at 4.5–6.7 kHz), and with `labiodental` 1.2 it carries /f/ up too. class_fit now also fits `oral_noise.labiodental`, `frication` and `obstacle_reach` by class.
  - th: low bands +11 to +15, high −5 to −14. voiceless stops +5 to +11 throughout; hh +5 to +17.
  - iy: −7 to −12 at 144–317 Hz and −10 to −13 at 2.3–3.4 kHz (weak F1 and F2/F3 regions). ae +8 to +16 at 0.4–2.6 kHz.
  - m: +8 to +16 at 0.4–3 kHz. w: −12 to −20 below 400 Hz, +11 to +13 at 1.5–2 kHz.
- **Pitch** (same renders, semitones re 100 Hz, p5/25/50/75/95): bdl −0.14 / 2.26 / 3.56 / 4.86 / 6.26; ours −1.64 / 1.96 / 3.56 / 4.66 / 6.46. The median is right; our low tail sinks 1.5 st (91 Hz against 99), the final falls.
- `analysis/profile_cycle4.sh`: class_fit (best held-out mean) → profile_fit on source, noise, variation, cycle, plan_pitch (kept only if its held-out beats its start's) → duration_calibrate. Starts `data/profile_{bdl,rms,slt}_s2.*`; rms and slt carry both continuous switches, organs_vtl_tongue and their duration plans.
- **class_fit bdl_c2** (146 settings, 18 generations, from profile_bdl_s2): held-out combined 6.8 → 7.0 (spectrum 5.9 → 7.0, dynamics 5.1 → 5.0, phase 10.0 → 10.6, pitch 10.0 → 9.4). Training class cost 1.09 → 0.66 median. `data/profile_bdl_c2.*`.
  - Largest moves: IY drop 0.02 → 0.43; dental_jaw 3 → 0.4; labiodental_jaw 4 → 0.8; postalveolar aperture 0.6 → 1.23; rhotic_press −0.35 → −0.01, rhotic_up 0.3 → 0.03; w_aperture 0.45 → 0.05 (its floor); frication 2.5 → 1.5.
  - The labiodental obstacle weight wandered (0.84–1.41) rather than falling.
- **/f/ diagnosed** ("Before the fifth of February, a few fine fish were found far off."; a0101):
  - In the class-fitted voice /f/ seals (labiodental jaw 0.8°): oral flow 0–2 cm³/s, mouth pressure at lung pressure. Its low band is release and squeeze transients; its highs are silent (−61 dB). The fit preferred a sealed /f/ to an open one, because:
  - Open (profile_bdl_s2, 0.06–0.14 cm²), two faults:
    - (1) The folds keep swinging: glottal area 0.17–0.22 cm² at f0 with about 2200 dPa across them (the 2025 amplitude fades only from 2000 to 1000 dPa). That gives voicing at 100–500 Hz. The broadband periodicity check (0.05–0.45) missed it under the frication. `kinematic.abduction_fade_start/end` 0.04/0.10 stops the swing.
    - (2) The bulk flow (300 cm³/s) ramps. With every noise source off, the /f/ band levels are 0–50 Hz −18 to −22, 50–100 −20 to −28, 100–200 −24 to −52 dB re speech level. bdl: −42 to −47, −43 to −67, −41 to −64, with his energy at 4–8 kHz (−27 to −40). The glottis closes 0.17 → 0.03 cm² in 16 ms after /f/ (adduction_move 0.02 s): flow 330 → 62 cm³/s.
  - The highs of an open /f/ are the labiodental obstacle source: −12 to −17 dB above 3 kHz, −33 to −44 without it. bdl's are about 10–15 dB under ours, so roughly labiodental 0.3. With it off, the VTL aspiration (gain 12.7) is what remains up there. Jet frication contributes almost nothing to /f/.
  - Long-term spectrum (a0101–a0120) re the 150–300 Hz band: bdl 10–90 Hz −14 to −19 dB, ours 0 to −5.
- New, all off by default:
  - `recording.highpass` (BodyLayout; speak applies a two-pole Butterworth after the room; speak output identical with cmp): at 60/80/100 Hz the /f/ 0–50 Hz band falls 10/14/18 dB, 100+ Hz unchanged.
  - Plan `breath.syllable_rise` and `syllable_decay` (the vowel's breath pulse, 0.35 / 0.88; plans identical with cmp).
  - profile_fit: groups plan_breath (unstressed and voiceless pressure, syllable rise and decay) and recording.highpass in noise. class_fit TREE: the abduction fade (credited to the voiceless classes) and sinus_neck_scale (nasals).
- **Dynamics** (profile_bdl_s2, two seeds; std of octave levels within a token, bdl / ours dB): nasals 4–7 / 1.5–3.3; l 5–8 / 3.5–4.5; ae 3.6–5.5 / 1.1–2.4; uh 2.5–6.8 / 0.5–3.0; z low bands 11–14 / 6.7–7.8. Ours holds too steady within tokens.
- Chain `bdl4 → rms4 → slt4` (profile_cycle4: 3 rounds, class 9 generations, tree 10) started 06:45 EDT.
- **Vowel loudness envelope** (total level over each vowel token, median re its peak, at 0, 0.1 … 1 of the token; a0101–a0120):
  - bdl: −2.0 −1.2 −1.0 −0.9 −0.8 −0.9 −0.8 −1.2 −1.5 −2.3 −3.2 (flat-topped).
  - ours: −2.9 −2.6 −2.0 −1.5 −0.9 −0.7 −1.3 −1.8 −2.5 −3.6 −4.4 (already more peaked).
  - So the dynamics shortfall is spectral change within tokens (formant transitions, articulator speed), not the breath pulse. New profile_fit group plan_motion: vowel_hold, release_move, jaw_lead, tip and dorsum release moves, liquid_move, adduction_move.
  - `profile_cycle5.sh` = cycle 4 + plan_motion + line-buffered logs. rms5 and slt5 follow bdl4 (`/tmp/es/run_profiles5.sh`).
- class_fit: NATO mode (REFDIR "nato": the 26 words in one render, the reader's clips re-cut per render by DTW, even words train, odd held out).
- **bdl4 round 1** (class 9 generations, then tree 10): held-out 7.0 → 7.2 (class) → 7.3 (tree, kept; training 8.30 → 8.58, one improvement in 10 generations; pulse_shape at its bound 1). Snapshot /tmp/es/snap/bdl4r1.*: two seeds, held-out spectrum 7.9, dynamics 5.2, phase 10.0, pitch 9.6, timing 5.6, combined 7.2.
  - Class costs (sum 0.745) led by timing: uh 0.089, ae 0.057, ow 0.017, y 0.016.
- **Our phones and bdl's labels are different words.**
  - Our planner gives function words weak forms (AX) where ARCTIC's labels keep the citation vowel: 16 of 319 phones in a0101–a0110 ("that", "and", "would", "the"). bdl's held-out /ae/ class held 14 tokens, ours 6.
  - So the pooled duration calibration shortened our stressed /aa ae/ by 0.67 to match bdl's mix of stressed vowels and short function words (held out after it: ae 0.63, uh 0.57, ow 0.78 of bdl).
  - bdl does reduce them (his weak-form /ae ah/ last 45–65 ms, our AX 44); our full vowels were the short ones (aa 0.068 against 0.110 s, ae 0.063 against 0.095).
  - FIX (measure v5, deployed 07:30 EDT): `fidelity.aligned_names` renames our phones to the reference's labels where the sequences align one to one (difflib; equal-length replacements position by position). `ours_set(dir, ids, ref_dir)`; profile_fit, class_fit and variant_test pass it. On the bdl4r1 renders: timing 5.6 → 4.6, combined 7.2 → 6.9 (same words now compared).
  - `duration_calibrate.py` now pairs aligned tokens (keyed by our phone). On bdl4r1 it starts at rms log ratio 0.098 and settles: AA 0.95, AE 0.94, AH 1.07, UH 0.74, Y 0.85, T 1.10, P 0.89 (the pooled method had AA and AE at 0.67).
- **Nasal boundaries** (500–2000 Hz level re the utterance's loud frames, medians around the labelled boundary, bdl4r1 renders against bdl a0101–a0120):
  - Into a nasal: bdl's vowel +2 to +3 dB, falling to −7 within 25 ms; ours 0 → −5 within 5 ms (a smaller drop, the murmur too bright).
  - Nasal → vowel: bdl −8 → +8 to +9 within 10–15 ms; ours −5 → +3 to +5 (the vowel stays 4–5 dB weak).
  - Nasal → stop: bdl −6 → −19; ours dips to −30 at +30 ms.
  - The planner closes the velum from 10 ms before a nasal's end over 70 ms (open to 0.45, port 0.2 cm²), so the vowel after it is nasalised. New plan options `nasal.velum`, `open_lead`, `open_move`, `close_lead`, `close_move`, `anticipation`, `anticipation_move` (defaults the old constants; plans identical with cmp) and a profile_fit group plan_nasal.
  - class_fit TREE: `vtl.aspiration_gain` (hh +13/+6/+7/+15 dB by quarters of the band; credited to hh and the voiceless stops).
- `profile_cycle5.sh` splits the shared fit in two (source: source, noise, variation, cycle; plan: plan_pitch, plan_breath, plan_motion, plan_nasal). 49 settings in one ES of 9 candidates was too thin.
- **Stops** (bdl4r1 renders, level re speech over each token at 0, 0.1 … 1, medians):
  - p t k under 500 Hz: bdl −6 −10 −16 −28 −34 −40 −40 −34 −32 −30 −28; ours −7 −11 −17 −22 −26 −29 −27 −24 −17 −15 −13 (the closure 11–13 dB too loud). Above 2 kHz mid-closure: bdl −58 to −64, ours −53 to −56.
  - b d g under 500 Hz: bdl −4 … −16 (a strong voice bar all through); ours falls to −29 (voiced closures devoice).
  - class_fit TREE adds `walls.closure_pull`, `walls.closure_relaxation` and `posture.closure_press`, credited to the voiced stops and (half) the voiced fricatives.
- **bdl4 round 2**: class held-out flat at 7.2 (v4). The tree stage (v5 measure; start held-out 7.0, timing 4.6) found nothing in 10 generations, all candidates 6.8–8.0 against the start's 8.27. The ES is not elitist and its first step (0.15 in unit space over about 40 settings) is too wide near an optimum. `PROFILE_SIGMA` (profile_fit; default 0.15), set to 0.08 in cycle 5.
  - The aligned duration calibration then moved the plan (0.098 → 0.001).
- **Voiced closures** ("A big dog dug a bigger bag, and Bob added a good deed.", bdl4r1):
  - The folds swing all through a /b/ closure (glottal area 0.027–0.052 cm² at 0.25 ms resolution; a 4 ms trace had aliased it to a dying swing).
  - The mean mouth pressure goes negative (−1500 to −2200 dPa) under closure_pull 0.25 while the lips are sealed.
  - The voice bar radiates only through the skin, and the tree fit had cut `walls.skin_gain` to 0.037. Under 500 Hz at mid-closure, median re speech: skin 0.037 −27.9; 0.3 −21.3; 1.0 −12.1; 3.0 −3.9 (bdl −9 to −16).
  - Held out (two seeds): base 6.9; skin 0.3 6.8; skin 1.0 6.7 (stop_voiced not better; affricate_voiced +3.3, dh +1.8). The skin radiates every wall motion up to its 250 Hz corner, not just the voice bar.
  - closure_pull 0.7, closure_relaxation 0.6 and closure_press 0.012, singly or together: no change to the voice bar (−34 … −36 for /b/).
  - class_fit TREE: skin_gain and skin_corner, credited to the voiced closures.
- **bdl4 done** (3 rounds; aligned measure from round 2's tree stage on): final held-out 7.1–7.2.
  - Round 3: the class stage peaked at its start (its held-out fell to 6.7 by generation 9 while its training cost kept falling: 147 settings over 30 sentences fit noise in rare classes; uh has 3 tokens). The tree stage found nothing.
  - Normalized error by group (two seeds): dynamics 0.33, timing 0.29, spectrum 0.16, pitch 0.11, phase 0.09.
  - Phase scatter: bdl 0.37 / 0.49 / 0.57, ours 0.19 / 0.34 / 0.67. Pitch: our p5/25/50/75/95 run 0.3–1.3 st low.
  - bdl5 (cycle 5, from bdl4) will train on a0001–a0060.
- **Weak forms' length**: `timing.function_vowel_scale` (PlanOptions, default 1; TimedPhone::function_word from Word::function; plans identical with cmp). duration_calibrate fits it from our AX aligned to bdl's full-vowel labels (key "weak forms"): bdl4 → 1.23 (rms log ratio 0.019 → 0.011 in 3 passes). Held out (two seeds): timing 5.3 → 6.1 (bdl against himself 5.7–6.9), combined 7.1 → 7.5. `data/profile_bdl4.*`.
- **`analysis/pitch_calibrate.py`** (closed form): render the training sentences, compare F0 percentiles. pitch_median moves by the median's gap, pitch_spread by the interquartile ratio, lowest_pitch by the rest of the 5th percentile's gap.
  - bdl4 on a0001–a0030: gaps (st, p5/25/50/75/95) −0.72 +0.40 +0.20 −0.10 −0.80 → −0.10 +0.00 +0.10 +0.00 −0.40 in 2 passes (rms 0.52 → 0.19). pitch_median 122 → 117.5, spread 1.5 → 1.76, lowest_pitch 80 → 94.7.
  - Held out (two seeds): pitch 10.2 → 11.3, combined 7.5 → 7.5 (phase 11.7 → 10.4, inside its seed range).
  - `profile_cycle6.sh` = cycle 5 + this stage after the durations (not for NATO). slt5, nato5 and bdl5 run through it (same-length in-place launcher edit). bdl5 starts from bdl4 with the weak-form and pitch plans (`data/profile_bdl4.*`).
- `branches.nose_smoothing` (the nose's hop loss, fixed at 0.4 in Speaker) and `tract.hop_loss` (the vocal tube's; below 0 keeps the Speaker's 0.02) as layout settings; speak output identical with cmp. class_fit TREE: nose_smoothing credited to the nasals. profile_fit tract group: hop_loss. Cycle 6's source stage now includes the tract group.
- **/z/ devoicing** (under 400 Hz re speech at 0, 0.2 … 1 of the token; /z/ before no vowel): bdl −5 −5 −13 −22 −33 −34 (14 tokens); ours −7 −13 −17 −19 −19 −18 (voicing held). The plan's `glottis.final_devoicing` (open the folds from 40% of the token) was 0. class_fit consonants now include it (credited to /z v/).
- **Closure energy, still open** ("keep" in a0101, bdl4w): inside the /k/ closure the oral flow swings −11 … +42 cm³/s with a period of about 10 ms while the velum seal holds; under 500 Hz −15 to −19 dB re speech (bdl about −40). The start of /t/ /p/ closures also carries 20–30 ms of voicing tail (−4 → −24 dB).
  - Not motor step targets: motors.smoothing 1 / 2 (two seeds) 7.5 / 7.6 against 7.5.
  - Not motor wander: wander_gain 0 / 0.3 → 7.6 / 7.5 (voiceless stops −0.5 dB).
  - Suspected: the yielding walls ringing after the pressure step.
- Cycle 6 calibrates (durations, pitch) once before round 1 too: rms5 sat at held-out timing 1.1 through a whole round on old pooled durations.
  - walls.loss_scale 2 / 4 (the fit had it at 0.40): 7.5 / 7.6 against 7.5. Not the walls either.
  - What the /k/ shows: sealed only 33 of its 91 ms (closing 15 ms, the contact sliding from node 29 to 33 toward /iy/, released at 60%). The token's last 40% is burst and aspiration with the release flow rising (to 245 cm³/s). bdl's stops hold closure to about 70%, and their aspiration has almost nothing under 500 Hz.
  - The planner releases a stop `burst_scale × post_burst_interval` (capped at 40% of the phone) before its end. class_fit consonants now include burst_scale, seal_lead, aspirate_gap, and the tip and dorsum close moves (credited to the stops, nasals and affricates).
  - My "under 500 Hz" probes include sub-100 Hz flow ramps, which the fidelity bands (from 100 Hz) do not count. The class errors place the voiceless stops' excess at 144 Hz–1.3 kHz (+9, +9 dB by quarters).
- **Level normalisation exposes out-of-phone sound** (bdl4p renders against bdl a0101–a0120, 10 ms frames within 40 dB of each utterance's loudest):
  - Frames outside every phone: ours 16% of the window at −13.5 dB re speech level; bdl 1% at −22. They pull our speech level about 2.5 dB down, so nearly every class reads too loud: vowels +1.7 dB on average, consonants +3.1 (ae +7.3, hh +9.4).
  - Silence before the first phone: bdl −58.6 dB re the loudest frame, ours −41.2.
  - The room is not the cause (dry ROOM=0: lead-in −53.5 / −43.1 against −50.8 / −42.8).
  - In a0101: bursts up to +9 dB re speech level just before the first phone, and after the last, single-millisecond spikes at +3 to +4 every ~145 ms.
- **The spikes: a geometry flip at the tongue tip.** At 3.062 s nodes 37–40 lost 0.4 cm³ in 2 ms (node 38 1.758 → 1.367 cm²), and at 3.206 s got it back. The organs moved by jaw 5.699 → 5.700°, tongue y 8.010 → 8.011.
  - The cut at x 7.35 switched its lower wall: arc 14.079 → 13.730, area 2.455 → 1.017 cm², every later arc −0.093 cm (teeth arc 14.644 → 14.551). The motor wander pushes the tip back and forth across it.
  - The physical tongue already avoided this with `smooth_airway` (Articulators::nearest_width: width = the hull point's distance to the nearest organ). The VTL-tongue branch never set it.
  - New `tongue.vtl_smooth_airway` (default 0; speak identical with cmp). On: no spike after speech (all ≤ −27 dB re speech level), the +9 dB pre-speech burst gone (+3 left at the breath's onset).
  - Held out, unrefitted: 7.5 → 7.0 (spectrum 7.9 → 7.2, dynamics 4.9 → 4.0; affricate_voiced −1.5, th −1.4, dh −1.3, f −1.0; s +4.5, sh +2.7, w +1.6: nearer-distance widths narrow the sibilant channels the fit was tuned for). Testing with a class refit (bdl_sa, chain paused with SIGSTOP).
  - Still open: the 0.67 s lead-in carries inhalation and breath-onset noise at −13 to −25 dB re speech level (bdl's lead-in is silent).
  - Refit (bdl_sa: class_fit 9 generations from bdl4 with the switch on; one seed): held-out 7.1 → 7.3 → 7.4 → 7.3 (spectrum 7.2 → 8.3, dynamics 4.2 → 4.4). Without the switch, bdl4 scores 7.5 (two seeds), so the measure is neutral to slightly negative.
  - Dynamics fell with continuity again (4.9 → 4.4): the jolts had padded it. Adopted anyway (vowel-loud clicks outside phones are real artifacts the measure barely sees): bdl5 starts from bdl_sa (`data/profile_bdl_sa.*`); slt5 and nato5 starts carry the switch. rms5 (running) does not.
- `breath.inhale_pressure` (plan; default −110 Pa, the old constant; plans identical with cmp), in profile_fit's plan_breath group (−150 … 0).
- **Phase scatter fixed by per-cycle lag and less aspiration** (bdl4w with bdl4p plan, held-out a0101–a0120):
  - Scatter in the 0–1 / 1–2 / 2–4 kHz bands, rad (bdl 0.37 / 0.49 / 0.57), one seed:
    - base 0.201 / 0.354 / 0.675
    - cycle_shape 0.15: 0.215 / 0.351 / 0.689; 0.35: 0.235 / 0.364 / 0.681
    - cycle_lag 0.15: 0.250 / 0.416 / 0.745; 0.3: 0.269 / 0.426 / 0.703; 0.6: 0.299 / 0.457 / 0.760
    - lag 0.6 + shape 0.35: 0.348 / 0.453 / 0.771
    - + aspiration gain 6 (from 11.8): 0.347 / 0.455 / 0.634; gain 3: 0.349 / 0.457 / 0.612
  - Two seeds: base phase 10.4, combined 7.5. Lag 0.6, shape 0.35, aspiration 3: phase 24.6, spectrum 8.3, combined 8.0 (hh −2.6 dB, ih −1.8). Lag 0.7, shape 0.4, aspiration 2: phase 23.5, combined 8.0.
  - bdl against himself: phase 19.7–26.7. Ours is now inside that range.
  - The tree ES never moved cycle_lag from 0, and its bound (0.4) was below what bdl needed. Bounds widened (lag 1.2, shape 0.8).
  - `analysis/source_calibrate.py`: pitch as pitch_calibrate, plus cycle_lag toward the 0–2 kHz scatter and aspiration_gain toward the 2–4 kHz scatter. Replaces the pitch stage in cycle 6.
  - bdl5's start carries lag 0.6, shape 0.35, aspiration 3.
- **100-sentence sets** (measure v5):
  - bdl a0201–a0300 vs a0301–a0400: spectrum 18.7, dynamics 13.6, phase 21.9, pitch 19.2, timing 13.8, combined 16.3. Against rms: −0.6; slt: 0.2. Larger sets raise the ceiling from about 10 to 16 dB.
  - Ours (the bdl5 start: bdl_sa with lag 0.6, shape 0.35, aspiration 3) on a0201–a0300, one seed: spectrum 7.8, dynamics 4.8, phase 14.3, pitch 5.3, timing 12.4, combined 7.5.
  - Pitch fell there because bdl's F0 moves between sessions: a0201–a0300 p5/25/50/75/95 = 0.86 / 2.76 / 4.26 / 5.86 / 7.86 st, a0101–a0120 −0.14 / 2.26 / 3.56 / 4.86 / 6.26; ours −0.54 / 1.66 / 3.16 / 4.46 / 6.36.
  - id ranges take a step ("a0001-a0593/20"; fidelity and render_set). Cycle 6 calibrates durations on a0001-a0593/3 and the source on a0021-a0593/25 (which misses the held-out a0101–a0120).
- **Dynamics is the sound, not the boundaries** (bdl5 start against bdl a0101–a0120, two seeds):
  - Transitions through tokens, a0201–a0300, median level re speech level at 0, 1/8 … 1 of the token:
    - /s/ 250–500 Hz: bdl −6 −23 −34 −38 −40; ours −10 −18 −25 −31 −34. 4–8 kHz: bdl −20 −11 −5 −3; ours −17 −14 −9 −5.
    - /f/ 250–500 Hz: bdl down to −51 by 0.4; ours −25 … −35. 4–8 kHz: bdl flat about −25; ours −34 → −16.
    - Nasals: our murmur 5 dB high at 250–500 Hz; its highs start vowel-like (−21 against −29) and end too dark (−43 against −37).
  - The abduction fade the class fit had left inactive (start 0.036 > end 0.008). With 0.03 / 0.08: f −2.3, sh −2.0, s −1.2 dB, but dynamics 4.5 → 4.2 and combined 7.8 → 7.7.
  - `motors.move_scale` (Speaker; every articulator's move time, not lungs, pitch or glottis; default 1, identical with cmp): 0.6 → combined 7.5 (dynamics 3.7); 0.8 → 7.9 (dynamics 4.4).
  - `FIDELITY_BOUNDARIES=mapped` (bdl's labels carried onto our audio by MFCC DTW): base 7.8 → 6.6, dynamics 4.5 → 4.4, timing 6.1 → 3.7, phase 21.3 → 14.0. The cross-voice DTW jitters the boundaries, and dynamics stays the same, so planned boundaries stay.
  - New fit settings: profile_fit plan_motion `motors.move_scale` (0.3–1.5); class_fit `glottis.abduction_move` and `fricative_abduct_lead` (credited to the voiceless obstruents).
- **The vowels' upper spectrum** (bdl5 start against bdl, ERB bands, mean over vowel classes, ours − bdl, dB):
  - Every second band from 144 Hz: +5 +1 +2 0 +6 +5 +5 +2 −4 −10 +4 0. At 2.0 / 2.3 / 2.6 / 3.0 / 3.4 / 3.9 / 4.5 / 5.1 / 5.8 kHz: +2 0 −4 −8 −10 −3 +4 0 0.
  - Fine spectra (a0201–a0300, vowel middles): bdl carries an F3–F4 cluster at 2.3–3.6 kHz (−1 to −4 dB re his peak) and falls off above 3.7 kHz. Ours is 10–15 dB weaker at 1.8–3.7 kHz, with peaks at 2.65–2.85 and 4.15–4.6 kHz.
  - Not the source: pulse_shape 0.5 / 0, phase_lag 0.8 / 2.0 leave the 3.4 kHz band at −10 or −11.
  - Not the losses: viscous_scale 1 / 0.5 and hop_loss 0 / 0.005 (−11). radiation_resistance_scale 0.7 changed nothing at all (to check: the slot may not reach the articulatory body).
  - Piriform sinuses (the zero near 4–5 kHz): scale 0.6 → 4.5 kHz +11; 1.4 → 3.4 kHz −13; smoothing 0.3 → −17.
  - New `anatomy.epilarynx_height` and `anatomy.epilarynx_area` (Speaker → MidsagittalTract; defaults identical with cmp):
    - area 0.6 / 0.4 → 3.4 kHz −14 / −17, 4.5 kHz +7 / +6
    - height 2.8 / 3.3 / 3.8 → 3.4 kHz −10 / −10 / −11, 3.9 kHz −8 / −7 / −7, 4.5 kHz −1 / −2 / −6 (combined 7.7 / 7.6 / 6.3)
  - New `recording.presence_db`, `presence_hz`, `presence_q` (speak's room; RoomOptions::presence_q; default 0 dB, identical). 8 dB at 3.3 kHz, Q 1.5: 2.6–3.9 kHz +0 −2 −4 +2, but 4.5 kHz +8; combined 7.6 against 7.7.
  - These go to a profile_fit group "shape" (presence dB, Hz, Q; epilarynx height and area), in cycle 6's source stage.
  - The 0.8–2 kHz excess (+5) is common to every vowel too.
- Held-out checks in class_fit and profile_fit now average each group's dB over two render seeds (73021, 11); class_fit renders them in parallel. One seed's phase group alone ranged 11.3–16.0 dB over a fixed voice, and the stages' acceptance compares these numbers.
- rms5 (cycle 5): round 1 class 2.5 → 4.5; tree and plan rejected (the plan stage had held-out pitch 10.3 → 15.2 at one point, but combined did not rise); the aligned durations took timing 1.1 → 4.0, held-out 5.7. Round 2's class and tree stages: no gain.
- **Status, 2026-10-18 11:15 EDT** (measure v5: aligned labels, pooled absolute phase):
  - Ceilings, the speaker against himself or herself: 20-sentence sets about 10 dB; 100-sentence sets bdl 16.3, rms 13.6, slt 13.0. Another speaker scores about 0.
  - bdl (the bdl5 start: VTL tongue with smooth airway, continuous geometry, class-fitted organs and consonants, aligned durations with weak-form scale, calibrated pitch, per-cycle lag 0.6 / shape 0.35, aspiration 3):
    - held out a0101–a0120, two seeds: about 7.7–8.0 (spectrum 8.0, dynamics 4.4, phase 21, pitch 11.6, timing 6.1)
    - a0201–a0300, one seed: 7.5 (spectrum 7.8, dynamics 4.8, phase 14.3, pitch 5.3 (session drift), timing 12.4)
  - rms (rms5 after 2 rounds): held-out 5.8.
  - The largest remaining gaps against the 100-sentence ceiling: dynamics (4.8 / 13.6), spectrum (7.8 / 18.7: vowels short at 2.6–3.9 kHz and over at 0.8–2 kHz; voiceless obstruents +12 to +20 dB at 144–400 Hz), pitch across sessions.
- **Affricates** (bdl5 start, a0201–a0300; level re speech at 0, 1/8 … 1 of the token):
  - ch 250–500 Hz: bdl −11 −31 −47 −44 −41 −41 −34 −33 −28; ours −20 −13 −29 −36 −36 −37 −37 −38 −35. 2–4 kHz: bdl −31 −40 −38 −18 −9 −4 −4 −7 −10; ours −20 −13 −28 −34 −37 −31 −16 −14 −19. bdl hisses from about 40%; ours from 75%, and 10 dB weaker.
  - jh 250–500 Hz: bdl −3 … −40 (devoiced through the hiss); ours −2 −5 −19 −31 −21 −14 −13 −14 −12 (re-voiced). 2–4 kHz: bdl −15 to −7 in the hiss; ours −28 to −19.
  - b d g: bdl's voice bar −8 … −22; ours −2 … −26 (the closure's middle deeper than his).
  - The hiss starts at max(start + share × duration, end − affricate_hiss_max), and affricate_hiss_max was 0.06 s, so a long /ch/ hissed only its last 60 ms. New plan option `glottis.affricate_hiss_gap` (below 0: the voicing gap; plans identical with cmp).
  - class_fit consonants now include affricate_hiss_max, affricate_closure and affricate_hiss_gap.
- **Nasal murmurs** (mid-phone spectra, a0201–a0300, dB re bdl's peak):
  - m: ours 0–450 Hz +3 to +8 over bdl. 0.6–2 kHz 5–8 dB under: bdl has a second resonance near 0.9–1.1 kHz (−23/−24), ours −31 to −34. 2.3–3.9 kHz 5–10 dB over.
  - n: low +5 and highs +2 to +10 alike.
  - class_fit covers the highs (nose_smoothing) and the sinuses. The missing 1 kHz resonance is the nasal tube's own geometry. Open.
- **rms5 done** (3 rounds, cycle 5): held-out (two seeds) spectrum 4.1, dynamics 3.5, phase 11.8, pitch 20.5, timing 3.8, combined 5.8 (start 2.5). `data/profile_rms5.*`. A source calibration afterwards gave 5.7 (rejected). rms6 (cycle 6, 2 rounds, from rms5) queued last.
- **slt's phase**: our scatter 0.63 / 1.38 / 1.10 against her 0.39 / 0.73 / 1.09.
  - jitter 0.56 → 0.2 and flutter 30 → 10: no change.
  - Aspiration 13.4 → 3: 0.35 / 1.04 / 1.21; phase 5.9 → 12.0, spectrum 3.5 → 5.5, combined 3.5 → 4.2. Her aspiration was too loud across 1–2 kHz.
  - source_calibrate now sets aspiration by the 1–4 kHz bands together and lag by 0–1 kHz alone. Lag steps halved, pitch steps damped 0.6: on rms, lag went 0 → 0.24 → 0.08 and the median overshot.
  - slt5 restarted with the fixed rule (`/tmp/es/run_profiles6.sh`: slt5, nato5, bdl5, rms6).
- Samples sent to the owner (12:34 EDT, `out/samples_2026-10-18/`): bdl_morning (profile bdl_start), bdl_now (the bdl5 start), rms_now (rms5), slt_now (the slt start after round-0 calibration only), three sentences each.
- **slt5** (cycle 6): round 0 calibration (durations 0.144 → 0.030; aspiration 13.4 → 3.3) took held-out 3.5 → 5.2. Round 1: class 6.1, tree 6.5 (kept), plan 6.9 (kept; phase 20.4). Held-out pitch stays about 5.5: her a0101–a0120 block's F0 differs from the spread set (her own 20 against 100 sentences: pitch 1.3).
- **Fitting speed** (the owner, 14:20: "optimize that code if you didnt plan for it to take an hour"):
  - Per sentence, speak 2.25 s for 3.6 s of speech (1.6× real time); say --plan-only under 0.01 s; describe 0.5 s (pyin 0.48).
  - speak's profile (sample, 4 s): Tube::step_air 562, step_lining 475, refresh_geometry 393, pow 125 (the viscous law, each step next to fast nodes), add_cut 48, the 2025 folds 34. The physics at 96 kHz is the cost.
  - Compiler flags: -mcpu=native and -flto both give identical output (cmp), but CPU time is within noise (6.0–6.8 s against 6.0–6.5). Not adopted.
  - The waste was in the pipeline. profile_fit's start evaluation and held-out checks (30 and 40 renders) ran one by one in one pool worker while the others waited, about 100 s each and several times a stage. class_fit's two-seed held-out used 2 of 10 cores.
  - `render_set.render_many(jobs, threads)` (threaded subprocess renders): profile_fit renders those at WORKERS (10) at once, class_fit's held-out seeds at 5 each. Takes effect from each next stage.
- slt5 round 2: class 6.9; tree kept 7.3; plan rejected. Source calibration (spread set): pitch gap rms 0.18 st; phase 0.339 / 0.738 / 1.039 against 0.363 / 0.616 / 1.044.
- **slt5 done** (3 rounds, cycle 6): held-out 3.5 → 7.1 (round 3's plan stage 6.1 → 7.1). `data/profile_slt5.*`. Her held-out pitch swings (3.3–7.4): the source calibration matches the spread set, and her a0101–a0120 block's F0 differs from it.
- **The NATO reader is in the female range**: F0 p5/25/50/75/95 = 8.96 / 10.86 / 12.76 / 14.16 / 16.06 st re 100 Hz (median 209 Hz); ours, from the male VTL start, −3.14 … 6.86 (held-out pitch −12.7, combined −6.3). Phase scatter 0.66 / 1.21 / 1.30 (breathy).
  - Formants between bdl's and slt's: iy 421 / 2528 / 2949, aa 792 / 1301 / 2718, uw 443 / 1966 / 2575, eh 787 / 1723 / 2589.
  - nato5 restarted from slt5 with anatomy scale 0.94 (slt 0.87) and pitch +2.3 st (lowest_pitch 155): held-out at its first class generation 2.7 (pitch 4.7). `data/profile_nato_s3.*`. Chain `/tmp/es/run_profiles7.sh`: nato5, bdl5, rms6.
- **nato5 done** (3 rounds from the female start; ~15 min a round, one render per candidate): held-out about 3.7–3.8 (start 2.7 on the female base, −6.3 on the male). Round 1 tree and plan kept; later stages rejected; spectrum about 2, dynamics about 0. The measure is limited by 13-word halves (one take of 26 words, DTW-cut). `data/profile_nato5.*`.
- bdl5 started 16:12 (cycle 6, train a0001–a0060, 2 rounds, from bdl_sa with lag 0.6, shape 0.35, aspiration 3).
- **bdl5** round 0 calibration (durations a0001-a0593/3: 0.074 → 0.011; source on a0021-a0593/25: pitch gap rms 0.46 st, phase 0.377 / 0.471 / 0.567 against 0.362 / 0.542 / 0.569): held-out (two seeds) 8.4 (timing 8.3, phase 20.0, pitch 11.9). Round 1 class stage (60 training sentences): 8.7 (dynamics 5.5).
- **bdl5 round 1**: tree kept, held-out 8.7 → 9.3 (spectrum 8.2, dynamics 6.4, phase 22.0, pitch 12.6, timing 8.3). Plan rejected (9.2). `data/profile_bdl5r1.*`.
  - Settings moved: presence 0.5 dB at 4.3 kHz (unused); epilarynx height 2.3 → 1.99, area 1.20; highpass 18 → 51 Hz; hop_loss 0.02 → 0.004; viscous 2 → 1.26; wall loss 0.40 → 0.60; skin 0.15 → 0.50 (the voice bar); aspiration 3 → 0.3; jitter 0.45 → 0.88; flutter 46 → 60 (its bound, widened to 120); pulse_shape 1 → 0.41; frication 1.8 → 1.0; fade 0.046 / 0.045 (still inactive).
  - Then source calibration: phase 0.360 / 0.577 / 0.643 against 0.362 / 0.542 / 0.569.
- **bdl5 done** (round 2: tree 9.2 and plan 9.3 both rejected; final source calibration pitch rms 0.26 st, phase 0.392 / 0.537 / 0.558 against 0.362 / 0.542 / 0.569): held-out 9.3 (two seeds). `data/profile_bdl5.*`.
  - **a0201–a0300 (100 sentences, one seed): spectrum 9.1, dynamics 6.5, phase 18.2, pitch 12.9, timing 14.8, combined 10.5.** The bdl5 start had 7.5. Ceiling (bdl against himself): 18.7 / 13.6 / 21.9 / 19.2 / 13.8 / 16.3. Timing is above the ceiling.
  - Worst classes now: z 8.2, uw 7.8, ch 7.7, l 7.1, v 6.9, b d g 6.9 (rms dB).
  - bdl6 queued after rms6: cycle 6, 2 rounds, training on a0001-a0593/10, held out a0006-a0593/30 (60 and 20 sentences spread across the sessions, disjoint; every step of 20 or less hits a0101–a0120).
- **rms6 done** (cycle 6, 2 rounds from rms5): held-out 6.4 (rms5 5.8). `data/profile_rms6.*`. a0201–a0300 (one seed): spectrum 4.8, dynamics 4.3, phase 23.2, pitch 15.8, timing 6.5, combined 7.2 (his ceiling 13.6).
  - Worst: z 17.9 dB. rms's /z/ hisses like his /s/ (−5 to −10 dB re speech at 3.4–5.8 kHz); ours is 21–30 dB short there (no hiss: voiced_fricative_gap 0.04). jh likewise −17 to −22 above 2.6 kHz. class_fit consonants now include `glottis.voiced_fricative_gap`.
- **slt5 on a0201–a0300**: spectrum 7.2, dynamics 5.7, phase 17.4, pitch 4.1 (her F0 drifts between sessions), timing 12.4, combined 7.3 (her ceiling 13.0).
- rms7 and slt6 queued after bdl6 (cycle 6, 2 rounds, train a0001-a0593/10, held out a0006-a0593/30).
- **bdl6 done** (train a0001-a0593/10, held out a0006-a0593/30, from bdl5): held-out stayed about 9.5 (pitch 24–26 on the spread set: the spread calibration generalises). No stage improved: bdl has plateaued on the present structure.
- **What is left for bdl** (bdl5, a0201–a0300): normalized error spectrum 0.124, dynamics 0.221, timing 0.033.
  - Top class costs, dynamics first: dh 0.028, th 0.025, ch 0.019, m 0.019, z 0.018, v 0.016, n 0.015, oy 0.015, ng 0.015.
  - Octave spreads, bdl / ours (dB):
    - dh 11.2 14.2 14.5 13.1 13.1 11.1 / 6.6 9.3 9.8 9.5 11.8 10.3 (bdl's /dh/ is a near stop)
    - th 14.3 15.3 12.1 11.5 10.3 8.4 / 11.8 10.4 8.2 8.1 9.6 11.5
    - m 3.2 4.2 7.2 6.4 7.1 4.5 / 3.0 3.1 3.5 6.9 4.3 8.8
    - oy 1.7 2.3 5.6 8.0 6.1 4.4 / 5.7 6.5 8.8 8.9 6.6 5.4
  - Mean spread error by octave −0.23 −0.52 −0.74 −0.84 −0.94 +0.46.
  - New plan option `articulation.dh_press` (above −1 overrides dental_press for /dh/; plans identical with cmp). class_fit: dh_press credited to dh, dental_press to th only.
  - bdl7 queued after slt6 (from bdl5, spread sets).
- **rms7 done** (cycle 6, spread sets, from rms6): spread held-out 6.1 → 6.9 (round 2's tree kept). **a0201–a0300: spectrum 6.3, dynamics 5.8, phase 20.8, pitch 16.9, timing 6.5, combined 8.2** (rms6 7.2; his ceiling 13.6). /z/ is no longer among the worst (v 11.0, l 10.5, f 9.2, w 8.7, y 8.7, ng 8.4). `data/profile_rms7.*`.
- **slt6 done** (cycle 6, spread sets, from slt5): spread held-out 7.6 → 8.5 (both trees kept). **a0201–a0300: spectrum 8.0, dynamics 6.7, phase 22.7, pitch 13.5, timing 12.4, combined 10.2** (slt5 7.3; her ceiling 13.0). Training on the spread set fixed her session pitch (4.1 → 13.5). Worst: l 10.1, v 8.4, hh 6.9. `data/profile_slt6.*`.
- Summary at 06:5x EDT 2026-10-19, a0201–a0300 combined against each speaker's own 100-sentence ceiling: bdl5 10.5 / 16.3; rms7 8.2 / 13.6; slt6 10.2 / 13.0; nato5 about 3.8 (13-word halves).
- **/l/** (a0201–a0300, ours − reference, every second ERB band): bdl −2 −7 −6 −12 −10 −8 +3 +11 −1 +3 −3 −2; rms −10 −9 −12 −19 −18 −15 −1 +5 +5 +6 −5 +3; slt −5 −7 −8 −16 −13 −17 −13 −3 +4 +1 −5 −8. Our /l/ is weak from 250 Hz to 1.3 kHz for all three.
  - The lateral channels' area was a constant (0.12 cm²). New `anatomy.lateral_area` and `anatomy.lateral_reach` (MidsagittalTract; defaults identical with cmp). 0.2 / 0.3: /l/ −0.3 / +0.1 dB; held-out unchanged. Reach 3 / 4.5 cm: no change.
  - In "Please let all little lambs loll lazily" three of ten /l/ sealed (0.001 cm²): an /l/ released its side channels 20 ms before its end, after the next /l/ had opened them. New plan switch `articulation.lateral_hold` (default 0, plans identical): all ten open. No held-out effect (a0101–a0120 has no /l/ pairs).
  - So the weakness is elsewhere. Open. class_fit TREE: lateral_area, credited to /l/.
- **/l/ spectrum** (mid-phone, a0201–a0300, dB re bdl's peak):
  - bdl prevocalic: strong below 600 Hz (−1 to −4), a resonance at 0.8–1.25 kHz (−13 to −23), a deep trough at 1.4–2.5 kHz (−39 to −49), F3 2.75–3.4 kHz.
  - Ours: −7 to −10 below 600 Hz, −18 to −32 at 0.6–1.25 kHz, no trough (−28 to −35).
  - New **lateral pocket**: VoiceBody branch group 4, a closed tube (`branches.lateral_pocket` port cm² at full lateral, `_length` 4 cm, `_area` 0.6 cm², `_back` 3 cm behind the tip contact). Port = lateral flag × lateral_pocket; it glides and is skipped when shut, like the under-tongue pocket. Off by default; output identical with cmp.
  - In "Lily will feel well, little Ella, all alone in the hall." (band levels dB): port 0.3: 1.6–2.0 kHz 6.6 → −9.3; 0.6: the trough at 2–2.5 kHz; 0.6 and 5 cm: 1.3–2.0 kHz 7.7 / 6.6 → 0.3 / −12.3.
  - Held out (two seeds): /l/ −0.3 to −0.5 dB rms, combined 9.2 → 9.2–9.3. The low band's weakness remains (open).
  - class_fit TREE: lateral_pocket (0–0.8) and its length, credited to /l/.
  - /l/ against its preceding vowel, 100–600 Hz (a0201–a0300 medians): bdl −0.1 dB (79), ours −5.6 dB (84).
  - The source is not the cause ("Lily will feel well…", TRACE 0.25 ms): glottal swing 0.06–0.08 cm² and flow AC 31–59 cm³/s in /l/, the same as the vowels beside it; mouth pressure 520–910 dPa in /l/ against −34 to 970 in the vowels.
  - So the tract damps or misplaces /l/'s F1. A 0.3 cm² channel did not lift it. Open.
- Samples sent 2026-10-19 (out/samples_2026-10-19/): profiles bdl5, rms7, slt6, nato5, three sentences each.
- **Owner's verdict on the 2026-10-19 samples**: "really terrible"; "half the words cant even be understood"; "worse than it was before"; only slt gets most words right; "You, you would not keep the truth from me" is illegible in every voice but slt.
  - The fidelity measure rose (bdl 7.5 → 10.5, rms 7.2 → 8.2, slt 7.3 → 10.2) while intelligibility fell.
  - The measure compares phone-class averages (mean spectra, spreads, durations, pitch percentiles). It never asks whether a word can be recognised, so the class fit was free to trade articulation for average spectra (for example labiodental_jaw 4 → 0.8, w_aperture to its floor, rhotic_up 0.3 → 0.03).
  - This repeats [[organic-comes-from-bodies]]: fits and statistics did not move the owner.
  - Chain stopped 2026-10-19.

## Reversibility (2026-10-19, owner's /goal: "start with the known-answer round trip. go long until you have it going both ways")
- The owner: an epistemic utility understands the inside of the box. The engine must be **reversible end to end**: excited by sound, it performs the mouth shapes (hinted by the NATO codewords); driven forward by our folds, it gives the sound back, in = out.
- **`analysis/reverse_core.py`** (the reversible acoustic core):
  - Forward, a replica of Tube::step_air: waves hop one node per 96 kHz step; junctions at the nodes; link areas the mean of their nodes; the glottis end driven by the flow entering; the lips loaded R ‖ L; each link's Borda-Carnot kinetic loss (series_element with no viscous resistance).
  - Backward, `peel`: the tube driven by the sound, node by node from the lips to the glottis. At the lips the load gives the pressure step by step. On each link the arriving wave and the wave that left the far node fix the link's flow (U = (f' − b)/Z); from it, the entering wave (a = b + Z U + K U|U|/2) and so the near node's pressure; that junction gives the next link's arriving wave. One step per node.
  - Tool support: sustain has `VISCOUS=` (viscous scale), `HOPLOSS=` (the per-hop low-pass) and `AREAS=path` (node areas after the first step). Its `RAW=` gives the flow entering the vocal tube and the radiated flow, per step.
- **Known answer 1, a static tube** (sustain aa, Story 2008; walls ×1e6, no branches, no noise sources, viscous 0, hop loss 0; 120 Hz, 800 Pa, 0.6 s; 47 nodes, 17.14 cm):
  - Forward replica against the engine's lip flow: 29.2 dB linear; **47.3 dB with the kinetic loss** (band SNRs 44–98 dB).
  - Backward by inverting the transfer with an FFT: REJECTED. Synthetic linear round trip 4–36 dB by band: the nearly lossless tube rings past the signal's end, and the truncation smears the inverse.
  - **Backward by peeling: synthetic round trip 297.6 dB linear, 296.8 dB with the kinetic loss (exact). The engine's glottal flow recovered from its lip flow: 71.2 dB (swing 68.6 dB).**
- **Known answer 2: posture from sound** (`analysis/posture_recover.py`; geometry_probe with the VTL tongue → node areas, the narrowest of the cut profile over each node; sustain `AREAFILE=` renders exactly that tube; rigid lossless tube):
  - Choosing among the 11 Story vowels from the sound alone: by the peeled flow's energy 6 of 11 (REJECTED); twice-differenced 1 of 11 (REJECTED); **by ringing in the closed phases, 11 of 11**.
  - Continuous postures (true bdl5 /aa/, hint organs_vtl_tongue /aa/). Ringing alone found a posture equally still but wrong, a tube 0.73 cm² rms off: the true source is never fully still (its stillest 35% reads −17 dB, chink or not). Spectral ripple rated the impostor (4.79 dB) above the truth (5.80): REJECTED.
  - **Consistency with our own folds** (the owner's loop: the glottal flow a tube implies, the sound peeled back through it, against the flow our folds make driving that tube forward): truth −73 dB waveform mismatch, hint +4 dB, the ringing impostor 9.5 dB worse by harmonic levels.
  - The rigid lossless tube's landscape is a narrow well: 90% of the way from the hint −15 dB, the truth −73. Global search (differential evolution, 4428 evaluations, then Powell) stopped at −29.8 dB in another basin (tube 0.87 cm² rms off; re-rendered sound 15.9 dB waveform, 4.85 dB harmonic levels). The rigid tube's resonances are needle-sharp (radiation-damped only), so any F1 error rings strongly.
  - Tract length now continuous (whole nodes plus the rest as the open end's air mass, Tube::set_active_length; sustain `ACTIVELENGTH=`). The period finder takes the shortest lag within 90% of the peak (it had chosen 60 Hz for 120).
- **The reversible core with real losses**:
  - Viscous resistance (2(a − b) = (2Z + R)U + K U|U|) and the per-hop one-pole low-pass, both undone in the peel. Synthetic round trips 293–297 dB. Engine (rigid walls, viscous 1, hop 0.02): forward 55.7 dB, peel 74.9 dB.
  - Yielding walls matter: against the rigid replica the engine's lip flow sat at 10.7 dB (5 dB at 300–1000 Hz).
  - The walls' axial tension couples neighbouring nodes (an implicit tridiagonal solve), which a node-by-node peel cannot follow. Without it (sustain `AXIAL=0`) the sound changes by 31.9 dB (lip flow; glottal 42.3). For the reversible engine, local linings.
  - The walls' random support variation (wall_variation 0.10) is hidden state: `WALLVAR=0`.
  - numba versions `forward_tube` / `peel_tube`: each node's lining (cover on hoop and coupling springs, body on coupling and support) driven by that node's pressure; in reverse, by its known pressure history. Tube::lining(...) exposes each node's compiled lining; sustain `WALLDUMP=` writes it.
    - Synthetic walled round trip: 280.6 dB. Engine: forward 40.9 dB, peel 57.5 dB.
    - The wall-displaced area changed both links' admittances. New `Tube::set_static_admittance` (the air sees the rest geometry; the linings act through their flow; default off; sustain `STATICAREA=1`).
  - **With it, the engine and the reversible core are the same machine: forward 147.7 dB, peel 141.4 dB** (float32 logging). Local walls, viscous loss, kinetic loss, hop low-pass, radiation, continuous length: the tube runs both ways.
- **Known answer 3: posture from sound through the walled engine** (local steady linings, static admittance, viscous 1, hop 0.02; true bdl5 /aa/, hint organs_vtl_tongue /aa/; consistency over 15 periods from 0.1 s):
  - Truth −138.9 dB, hint +5.5 dB. Along the straight line from the hint to the truth the landscape falls smoothly (90% of the way: −14.2 dB).
  - An impostor posture scores −30.2 dB consistency and −21.4 dB lip-waveform mismatch. By the lip waveform alone the landscape is a needle (about −3 dB everywhere but at the truth, −386): REJECTED as a search measure.
  - Powell from the hint: −29.3 dB (754 evaluations), tube 1.64 cm² rms off, sound round trip 21.3 dB.
  - Twenty random starts of Powell over the nine organs: best −34.6 dB; ten land near −25 to −35 dB, ten near 0 dB. None reached the truth. The posture landscape has many shallow wrong basins: it is the wrong space to search in first.
  - Fault found: `render` left its temporary folder behind; 36236 folders filled the Mini's disk. It now cleans up.
- **Known answer 4: the tube from sound** (`analysis/tube_recover.py`; Levenberg–Marquardt on the residual driven − implied glottal flow, parameters the node log areas and the length past the last whole node, columns by finite differences in parallel):
  - From the true tube with 2% log-area noise: −95.0 dB consistency, **sound round trip 88.6 dB**. 10% noise: −81.0 dB. 30% noise: −24.2 dB (a wrong basin).
  - The tube it reached differs from the truth by 0.37 cm² rms, all of it a node-to-node zigzag in the wide front cavity (7.10 8.20 6.57 8.93 6.17 against 7.5 flat), left over from the starting noise. A ripple at the grid's own spacing moves the sound only near the half step rate: the grid's null space, not a second mouth. OBSERVED.
  - So sound → tube closes locally; what is open is reaching the true basin from the hint. The first run from the hint (16 smooth control points, bounded trust region) stopped at −27.9 dB with 47 nodes where the truth has 46: the residual jumps when the count of whole nodes changes. Now the count is fixed per run and several counts are tried.
  - From the posture basins, Levenberg–Marquardt over the nodes (trust-region with bounds stalled where it began; LM, the shape clipping the log areas, moved): seed 14 −33.0 → **−119.8 dB, tube 0.057 cm² rms off**; seed 2 stopped at −45.7, seed 9 at −34.8. Starting smooth (16 control points) from the hint, LM blew the areas up to the clip (tube 13–14 cm² rms off): REJECTED.
  - **Pipeline** (`tube_recover.recover`): 512 Sobol postures about the hint, each scored by consistency; the best 6 → LM over the nine organs (band 1.5 kHz, then whole); the best 3 → LM over every node (1.5 kHz, 3 kHz, whole ×2).
  - **bdl5 /aa/ from the hint alone: −94.2 dB consistency, tube 0.059 cm² rms off (46 of 46 nodes), sound → tube → sound round trip 91.4 dB.** Survey best −28.8; its posture search −40.8 (tube 0.77 off); the other two carried basins stopped at −40.6 and −51.3. MEASURED (one vowel; the other twelve running).
- **Known answer 5: a moving tube, both ways.** sustain `AREAFILE2=` with `GLIDE=start,duration` glides the rest areas linearly at each 1 ms tick; `GEODUMP=` writes per step and node the area, the air the frame draws (`Tube::frame_flow`), the link's viscous resistance (`Tube::link_resistance`) and the lining (float64). New accessors only; nothing else changes.
  - `reverse_core.forward_moving` / `peel_moving`: the walled core with per-step geometry, in Tube::step's order (the linings step on the last step's constants, the geometry recompiles, the air steps on the new one; each link scatters on the geometry of the step that made its waves; the frame's drawn air enters every junction).
  - bdl5 /aa/ gliding to /iy/ (resampled to /aa/'s 46 nodes) over 0.12–0.22 s, 0.35 s: **forward 150.0 dB, peel 137.8 dB; during the glide 150.9 / 136.7 dB.** The static peel through /aa/'s tube over the same stretch: −1.6 dB. MEASURED.
  - Frame motion (rest radius at 160 Hz critically damped, the displaced air through two poles at 250 Hz, linings recompiled every step, the viscous law gliding every 16 steps) depends only on the geometry's targets, so the peel takes it as known: the moving engine is reversible given its trajectory.
- **Robustness of the tube from sound** (bdl5 vowels, hints organs_vtl_tongue):
  - First pipeline on /ae/: −32.4 dB, round trip 21.8 dB; /ah/ −35.1. Every /ae/ posture basin had 47–50 nodes where the truth has 44, and the node search keeps its count. Several truths lie outside the hint spans (/ey/ tipx −1.49 spans, /ih/ drop +1.71): survey widened 1.5×.
  - Over 1024 surveyed postures, rank correlation with the tube's distance from the truth: closed-phase ringing 0.15–0.38, waveform consistency 0.02–0.14. Both are exact only near the truth. **Our own sound's harmonic levels against the sound's** (rms dB): 0.70 / 0.82 / 0.40 (/ae/ /ah/ /aa/) against the log-area shape's distance. MEASURED.
  - Why the levels are the same test: the peel is the forward pass's exact inverse, so our folds' flow equals the flow the sound implies exactly when our sound equals the sound. Per harmonic, log(ours / the sound's) is that test: the real part the levels, the imaginary part the phases. The waveform form is exact but rugged; the spectral form varies smoothly with geometry.
  - Our folds settle at a pitch the tube's load moves (period 761 against 769 steps, 0.05 spans from the true /ae/). Read at the sound's period, our upper harmonics were misread; each sound is now read at its own period, phases relative to the fundamental. The same drift over 15 periods is part of why the waveform test is rugged away from the truth.
  - **The node count was the cliff**: free, the posture → tube map jumps when the length crosses a whole node. From 0.05 spans off the true /ae/, a start rounding to 43 nodes (truth 44) never came back (+3.0 → +0.5 dB). `posture_recover.nodes(posture, count)` holds the count (the length past it carried by the open end's air mass, 0–2.5 nodes, as the engine clamps it). At a fixed count, the posture ladder (levels; spectrum with phase to 5, then 10 kHz; waveform) reached −52 to −72 dB from 5 of 6 starts within 0.1 spans; the posture itself ends 0.3–0.4 spans off at −70 dB (many postures, one tube).
  - A count-sweep by stretching the profile: REJECTED (it never chose another count; a basin 0.36 cm² rms off over the true 48 nodes had its 2 extra nodes at the lips, which a stretch cannot remove).
  - Spectrum rungs then waveform rungs over the nodes made near-exact posture tubes worse (/ae/ −63.3 → −44.7 dB). The node stage now runs the exact test only and keeps the best tube: **/ae/ −86.2 dB, tube 0.020 cm² rms off, sound round trip 82.0 dB** (first pipeline: 21.8 dB). /ah/ 40.1 dB, /aa/ 30.7 dB (stuck node searches; next items).
  - A posture tube with the length clamped at the most air mass started the node search on the fraction's bound: its difference column read zero, the length never moved. Started inside the bounds now.
  - **float32 was the floor**: from the true /aa/ tube with 10% log-area noise, the node search stopped at −70 dB on `xtol` (rejected steps): sustain's RAW is float32, and a 1e-3 step in a log area left ~1e-4 relative noise in each difference column, the size of the residual itself. sustain `RAWDOUBLE=1` writes float64 (default unchanged); the same search then fell −71.7 → −92.2 → −100.4 → −103.9 dB over four runs of 30. MEASURED.
- **Output error against equation error** (the decisive change):
  - With float64, all four free-running (output-error) results: /aa/ **123.2 dB** round trip (tube 0.058 cm² rms off); /ae/ 42.3, /ah/ 44.0, /ao/ 32.6. The last three stalled over the nodes from tubes 0.15–0.65 cm² off.
  - Cause: the free test runs our folds against their own candidate tube; a slightly wrong tube moves their pitch, and over 15 periods the flow drifts in phase: an output-error landscape, rugged.
  - **Equation error**: drive our folds with the tract the sound implies. The peel gives the first node's pressure and flow, hence the wave the tube returns toward the glottis, (p − U ρc / A0) / 2. `VoiceBody::force_return(on, wave)` has the folds read it in place of their tube's (`forced_flow()` the flow they pass; off by default); sustain `FORCEWAVE=path` (float64 per step). Default output bit-identical (cmp, two renders). The residual is the forced flow less the implied flow, step by step.
  - /ae/ at the truth: equation −247.2 dB, waveform −230.1. From the true tube with 10% log-area noise: **equation −44.3 → −113.9 → −133.1 → −142.4 dB, tube 0.005 cm² rms off**; the waveform test from the same start never moved (−19.4). From 30% / 50% jagged noise the equation error also stalls (−38 dB): it needs a smooth start, which the posture stage gives. MEASURED.
  - The equation error is now the test of the posture stage's last rung and of the node stage.
  - **With the equation error, from the hint alone: /ae/ −247.5 dB, tube exact (0.000 cm² rms), sound round trip 216.6 dB; /ah/ −129.5 dB, tube 0.015 cm² rms off, round trip 119.2 dB.** All three carried basins of each converged to the same tube (from −40 to −70 dB posture tubes). Before it: 42.3 and 44.0 dB. MEASURED.
- **All 13 bdl5 vowels from the hint alone** (equation error; survey 1024, 6 posture basins, 3 node searches): round trip /aa/ 224.5, /ae/ 216.6, /eh/ 232.2, /ey/ 246.3, /ax/ 188.1, /ih/ 132.2, /iy/ 127.7, /ah/ 119.2 dB: **8 of 13 exact** (tubes ≤ 0.019 cm² rms off). Not closed: /ao/ 47.4, /uw/ 44.7 (equation −79.3, tube 0.085 off, stalled), /uh/ 33.5, /er/ 12.4, /ow/ 2.1 dB. MEASURED.
  - All five are rounded or rhotic: a small lip opening (0.2–0.7 cm² at the last nodes), where the open end's air mass trades against the lip area. The count was wrong in the carried basins (/ow/: 45–50 nodes, none at the true 48; /er/: 46–47 against 44). Not an area clip (smallest true area 0.15 cm²).
  - The node stage now changes its count when no count is given: past 1.2 nodes of air mass it becomes nodes (the last area continued), and pinned at none, a node is given back; it runs until two runs gain under 1 dB (at most 10).
  - With the adaptive count: /ow/ 47 → 48 nodes (the truth's) and −33.1 → −97.1 dB over ten runs, still gaining ~1 dB a run; **round trip 2.1 → 67.7 dB**. /er/ unchanged (12.7 dB): its carried basins were at 46–47 nodes with 0.45–0.8 nodes of air mass, never pinned; the 43-node basin 0.17 cm² off was left behind (lower posture score). Next: carry all six basins.
  - Jacobian of the equation error at the true tube: the weakest directions are node-to-node zigzags in the wide regions (the grid's null space) for /aa/, /ow/ and /er/ alike (smallest singular values 0.015–0.022). /ow/'s strongest is 15× stronger than /aa/'s (condition 9.5e5 against 6.0e4): the small lip end, its resistance and kinetic terms going as 1/A and 1/A². A curved, steep valley, not a missing direction. Column scaling (x_scale "jac") from a smooth 15% perturbation of /ow/: −63.8 against −52.7 dB unscaled; neither converged. MEASURED.
- **/uw/'s floor**: from the true tube with 2% or 5% noise, the node search stops at −79.6 dB both times (truth −249.6). Near the truth the equation error is exactly linear in a random direction (|r| ∝ t from 1e-7 to 1e-2). The stopping point differs from the truth by a node-to-node zigzag (the near-null subspace), and along the line from it to the truth the error is bumpy at −80 dB (−79.6, −81.2, −77.6, −79.6, −85.8) before −246.9 at the truth. With the profile's bend penalized in the first runs (weights 1, 0.1, 0.01 on the second difference of log area, then none): −112.6 dB, tube 0.025 cm² off (3, 1, 0.3, 0.1: −106.0). A floor remains. The node stage now runs that schedule and carries all six posture basins. MEASURED.
  - Rerun (bend schedule, six basins): **/uh/ round trip 223.5 dB (tube exact), /er/ 71.3, /ao/ 63.6**; /ow/ 53.4 (67.7 in the previous run), /uw/ 44.7 (stopped at −79.3 again from its basin). Each winner was the 4th or 5th basin carried. **Vowels from the hint alone: 11 of 13 past 60 dB round trip, 9 exact.** MEASURED.
- **Known answer 6: a sequence of held tubes** (`analysis/sequence_recover.py`; sustain `AREASEQ=` keyframes, each phone held 70 ms, transitions 50–80 ms; `reverse_core.keyframe_schedule`). Four vowels keyframed (/aa/ /iy/ /uw/ /ae/, 0.45 s): through the geometry replica, forward 161.0 dB, peel 163.6 dB. The truth's equation error, /aa/ /iy/ /ae/ /uh/: −138.7 dB; one arrival 1 ms late: −21.3 dB.
  - First run: /aa/ (the first phone) exact (−247.7 dB); every later phone failed (−13 dB), the joint search with it (round trip 1.1 dB). Searched alone, a later phone's candidate is peeled from the sound's start through one static tube, and our forced folds (an oscillator with memory) reach the hold carrying that wrong history. Now each later phone is searched in its sequence (`PrefixContext`: the phones before as recovered, the candidate held from its hinted arrival); `tube_recover.recover(..., context=)` hears every rung that way.
  - Second run, each later phone heard in its sequence but at the hinted times (8 ms off): /iy/ still −12.9 dB. **Our folds are an oscillator whose phase is neutral**: an earlier mismatch (here the transition's timing) leaves them out of phase for good, so the equation error after it never recovers. In the first run the timing grid had found transition 1 within 1 ms with /iy/'s tube still wrong. Now: each transition's times first (with the hint's tube for the phone it reaches), then the phone's tube in its sequence, then the times again.
  - Timing first, then the tube in its sequence: transition 1 found at 103.0 / 147.6 ms (true 100 / 150) with the hint's /iy/ tube; /iy/ then −42.9 dB (0.275 cm² off); with it the times came to 100.2 / 150.0 ms. **Then rounds of times and the tube again (`tube_recover.node_search` in context): /iy/ −77.0, then −96.4 dB, tube 0.061 cm² off, transition 1 exact on the 0.2 ms grid**; transition 2 within 0.4 ms. MEASURED.
  - Phone by phone the errors piled up through the folds' phase: /aa/ exact, /iy/ −96.4, /ae/ −57.4 (its polish rounds went to −52.9), /uh/ −23.5 and its transition drifted (319 / 446 ms against 360 / 420). Now after each phone, a joint search over every phone so far and the times between (`prefix_joint`; bend weights 1, 0.1, then none), so the next phone is heard after a near-exact history.
- Owner, 2026-10-19 evening: changing the geometry or the scaffold (more pieces, bifurcations) is fair when a behaviour is painfully hard to reach; the test is whether the tuned behaviour looks nonphysical. Ideas to carry: per-piece viscosity in motion so that only the triggers to switch are sent (targets and switch times: what the sequence search already estimates); "picking up a voice" by replaying a person speaking. Wants to hear the results after the batch.
  - **With the joint after each phone: /aa/ /iy/ /ae/ /uh/ recovered from the sound alone (hints: organs_vtl_tongue postures, times 8 ms off): all six transition times exact (≤ 0.01 ms), tubes 0.28 / 0.05 / 0.14 / 0.04 cm² rms off; sound round trip 72.1 dB over the radiated flow (0.02–0.53 s).** Joints: 2 phones −99.2 dB, 3 phones −84.6, 4 phones −84.5 (against the sound's swing); the final joint over all could not go further (−79.8 dB, four runs unchanged). The engine's own wave files (each normalised to its peak, onset included, the rate of change of the flow) differ by 51.0 dB. MEASURED.
  - Listening files: out/reversible/ (original, recovered, and both with their difference ×100). Sent to the owner 2026-10-19.
- **Known answer 7: closures.** /aa/ → shut → /aa/ (keyframes 0.10/0.13/0.20/0.23 s). The forward replica is exact through every closure (155–164 dB).
  - Three shut nodes at the lips: peel 83.0 dB (3e-3 cm²), 56.8 (2e-3), −1.1 (1e-3); at 5e-4 and 1e-4 it overflows (NaN). Behind a long closure the links themselves are narrow, and the peel divides by each narrow link's admittance at every node behind it (about 1e4 a node, times the viscous resistance's ratio): three nodes amplify the last rounding by ~1e17. The sound carries almost nothing about a sealed cavity.
  - **One shut node** (inside the tract): peel 189–243 dB at 1e-3, 1e-4 and 1e-5 cm², the lip flow down 40–57 dB during the closure. A link's area is the mean of its two nodes, so the links beside a single shut node stay half open. **A closure one node long is a full closure and exactly reversible**: the reversible geometry closes in one node. MEASURED.
- **Known answer 8: frication, reversible.** The engine's jet (sound_jet) drives each source from the narrowest link's flow and the pressure two nodes upstream: in a peel from the lips those come after the source links. New `VoiceBody::reversible_jet` (layout `reversible_jet` = gain, default 0; sustain `REVJET=`): from the narrowest node (under 0.6 cm², from node 10) to the open end, each link's source = gain · (ρ/2)v² · turbulent · fade² · (v/2000) · √A · exp(−distance/1 cm), v = |U_link| / A from **the link's own flow the step before** (the jet passing it), the noise hashed per link and step (SplitMix64, so any link and step can be drawn again) and shaped per link as shaped_jet_noise. Default output bit-identical (cmp).
  - `reverse_core`: the same source in the forward replica and in the peel (each link's wave entered less S/2, S from that link's flow one step earlier, known by then).
  - /aa/ → one constricted node → /aa/, gain 0.8: near the teeth (0.08 cm²) noise −23.2 dB of the voice, **forward 201.8 dB, peel 206.8 dB (203.4 during the constriction)**; at the lips (0.15 cm²) noise −30.1 dB, peel 211.0 (213.3). A peel ignoring the jet: 19.5 / 29.5 dB. A mid-tract constriction (node 30, 0.1 cm²) stopped the render on a non-finite sample (open). MEASURED.
- **Known answer 9: the nose.** sustain `REVNOSE=1` (the nose's linings without axial tension, its air on the rest geometry; `Tube::scale_axial`), `VELUMSEQ=` (port area over time), `NOSEDUMP=` (port node, the nose's constants), `PORTDUMP=` (port area per step), `NOSEHOP=` (the nose's hop low-pass, default 0). Defaults bit-identical (cmp, velum shut and open).
  - `reverse_core.forward_branched`: the port a one-hop series element each step (VoiceBody::exchange_branches), the vocal tube with it as a side branch, then the nose; the sound lips + nostrils. Against the engine (/aa/, velum 0.01 → 0.4 → 0.01 cm²): **285.4 dB** (ignoring the nose 17.5).
  - The split: one microphone, two openings. A junction condition per step (the oral peel's pressure at the port node against the waves the port and the nasal peel imply) is met at the true nostril flow to 1.8e-11, and the glottal flow peeled through it comes back at 288.6 dB. MEASURED.
  - Solving it: Newton–GMRES on the whole condition (matrix-free, two peels per product): the residual stalled at 0.92 after 204 products, an hour and a half without a result (REJECTED: the condition is a causal shift of the share). Marching it, solving each nostril sample from the condition m steps earlier: matched the truth for 200 steps (2e-11), then diverged (rounding grew ~5% a step: REJECTED).
  - **`marching.py`**: the oral path peeled time-major (one diagonal per sample) to the port node, the port stepped *forward*, the nose run forward skewed (node j at time k + j, so the nostrils' flow is known m − 1 steps ahead of when the oral peel needs it). No solves. Each piece exact against the replica's internal signals (skewed nose 374.5 dB, time-major oral peel 269.2 dB). The coupled loop was unstable near 4.9 kHz once the port passed ~0.28 cm² (×300 per 256 steps): the oral inverse inside the physical loop. The engine's nose had lossless hops; with a nasal hop low-pass of 0.2 (about 0.985 a hop at 4.9 kHz, −4 dB a pass through the nose at 5 kHz, little below 1 kHz) **the march is stable: split 283.4 / 284.0 dB, glottal flow 289.1 / 289.2 dB (ports 0.4 / 0.8 cm²), 0.07 s for 0.4 s of sound.** 0.05 and 0.1: still unstable. A scaffold choice for reversibility (nasal passages are lossy); their bandwidths against real nasals: open. MEASURED.
  - **/ama/** (lips shut in one node 0.14–0.24 s, velum 0.5 cm²; nose hop 0.2): replica 182.2 dB; during the /m/ the nostrils carry the sound (lips −61.1 dB). The forward-nose march (B) held through the lips' closing and failed only once they were shut (143.3 ms): it inverts the oral path, which a closure makes weak. Peeling the nose and solving each nostril sample (A) fails at once even on a held /m/ (2.2 ms): REJECTED in principle. A per-step hybrid (A while shut) failed at the switch, A's own instability.
  - The rule: invert the strong path, run the weak one forward. Next: a true seal (no leak: the engine's 1e-4 cm² floor leaks −61 dB). Sealed, the lips' share is the isolated lip node's decay (forward), the nose is the single peeled path, the cavity behind the seal is a closed branch driven forward from the junction's pressure; for a stop with the velum shut, the region behind the seal runs forward with the folds until the release. A time-major engine that peels what is observable and runs the rest forward.
  - **True seal + modes** (`marching._march_sealed`): the oral path peeled while no seal falls in its reach; past that the nostrils carry the sound and the lips' own flow comes from the part in front of the seal run on alone (an isolated copy; sealed, it has no input). The front part runs forward from the port node's pressure throughout (`_oral_forward`, exact against the replica: 180.3 dB); at the release the oral peel's state is rebuilt from it (`_peel_from_front`). Then `_peel_below` takes the glottal flow from the port node down. Seal 0.01 cm²: **before the /m/ 222.2 dB, during the /m/ 128.7 dB** (the sealed-mode delivered wave reads ~2e-4 off: peeling through the nasal hop loss amplifies the high band ~(1.2/0.8)^30). Nasal viscous loss in place of the hop loss: ×5 and ×10 unstable while open, ×20 gives 104.6 dB in the /m/.
  - The release: once the nostril sample being taken passes the release, the lips radiate again and their flow depends on the port node's pressure up to 9 steps past what is known; holding the last value was wrong there and, fed back through the ahead copy, exploded within 5 steps. Next: a small Newton solve over the ~30 nostril samples after each release (the window run with them; the lips' flow = sound − share against the front model's).
  - **The nasal hop loss in tension**: the open-mouth march needs the nose damped near 4.9 kHz, the sealed one needs the nose peel well conditioned; the one-pole hop 0.2 cost ~106 dB at the half step rate over the nose, which the peel must undo. Stabilizers tried for the open-mouth march with a lossless nasal hop: nostril radiation scales (0.25–2): unstable; velum port resistance ×3 to ×100: unstable (REJECTED). **A shelf hop filter** g (x − b x′) + a y′ (`Tube::set_wave_shelf`, sustain `NOSESHELF=`; defaults bit-identical): the open-mouth march is stable once the shelf loses about −0.2 dB a hop at 5 kHz (a 0.6 b 0.57: stable; a 0.5 b 0.46, −0.16 dB: unstable); **a 0.8 b 0.79 loses −0.32 dB a hop at 5 kHz and only −0.47 at the half step rate (−14 dB over the nose)**, march glottal flow 289.4 dB. The shelf's inverse is a stable recursion in both nose peels.
  - /ama/ with the shelf nose and the seal: **before the /m/ 223.4 dB, during it 167.3 dB** (128.7 with the one-pole). MEASURED.
  - **The release** (`march_split`, Newton over the port node's pressures across the window: the lips hear them through the part in front of the seal before the march reaches them): consistent (residual 7.5e-4 at the true pressures, scale 2300); from a flat start Newton diverged, from the pressures one pitch period earlier (damped) it converged to 1e-10, **to another solution** (89.6 off the truth near the release), and the march on from it reproduces the sound with a glottal flow at 49.1 dB (61.8 after the /m/). A periodic prior (Tikhonov 1e-3–1e-1): 48.3 / 60.9 dB. Just after the release the lips are barely open and the port node's pressure is faintly heard: from the sound alone a family of glottal flows gives the same sound there. An information limit of sound-only inversion at a nasal release (~50–60 dB); in the search the folds' consistency (the equation error) is what decides among them. OPEN.
- **Known answer 10: a stop with the velum shut — peel where the sound reveals the inside, run forward with the folds where it cannot** (`analysis/observer.py`; sustain FORCEWAVE NaN = the folds read their own tube that step, free; defaults bit-identical). Blocks: the node-major peel over the whole sound records every node's state just before the closure reaches node 0's lookahead (the seal node's index steps before the seal); render A (folds forced before, free from there) gives the glottal flow through the closure; the replica runs forward on it (`_forward_steps`: the reversible core's step with the shelf and the seal) to one step after the last sealed step (the peel's first iteration inverts the step before's scatter: ending at the last sealed step left the flow after it at 26.3 dB); the node-major peel resumes from the forward state (`_peel_states`); render B gives the equation error outside the window.
  - Blocks against the engine (four keyframed vowels): forward steps 162.9 dB, peel 164.3, recorded state = forward state to ~1e-8 relative, forward from it 158.5, peel from a later state 164.3.
  - /aa/ → lips sealed (one node, seal 0.01 cm², velum shut) → /aa/: **implied glottal flow 223.6 dB before the closure, 185.4 after; through it the folds' own flow 250.9 dB and the lips 298.7 dB; equation error at the truth −186.0 dB.** MEASURED. The same principle should settle a nasal release: there too the sound cannot decide, the folds can.
- **Known answer 11: /ama/ by the same principle** (`observer.implied_branched`): windows the single-path peel cannot see = each seal and each stretch with the velum's port open, from when it reaches node 0's lookahead to one step past its end (`hidden_windows`); inside them the vocal tube, the port and the nose run forward on our folds' flow (`_forward_branched_steps`); outside, the shut-off nose rings on alone (`_nose_steps_isolated`) and its flow is taken off the sound before the peel. /aa/ → /m/ (lips sealed in one node, velum 0 → 0.5 → 0 cm², shelf nose 0.8 / 0.79): **implied glottal flow 231.1 dB where revealed, the folds' own 231.2 dB in the window, equation error at the truth −225.4 dB.** MEASURED. The release's ambiguity is gone: the folds decide what the sound cannot. The cost: inside the window the search has no equation error (only the forward sound against the sound, an output error, would test it); the split marches (`marching`) stay for nasalized vowels with the mouth open, if those need the equation error.
- **Words with closures in the search** (`sequence_recover.py`: phones `NAME@place` for a stop, `NAME@place^area` for a nasal; a stop's or nasal's tube is the phone before's with one node sealed; the residual from `observer.implied` / `implied_branched`, equation error where revealed, the forward sound against the sound inside the windows; velum keyframes tied to the phone times, `reverse_core.velum_port` = the engine's port to 5.6e-17). At the truth: /aa b aa/ −194.0 dB (closure 1 ms late −21.3), /aa m aa/ −172.6 dB (velum 0.45 for 0.5: −28.8; closure leaving 1 ms late −37.6).
  - /aa b aa/ first run: /aa/ exact, the closure's times exact; the release's 7 ms off and the vowel after it lost (−25.7 dB, round trip 1.3 dB): the window gave the search nothing (its residual zeroed), so the release was constrained only after the window. Now the stop's window also compares the forward sound with the sound. Rerunning; /aa m aa/ running.
  - /aa m aa/ with the whole velum-open stretch as one forward window: round trip 25.8 dB (velum 0.48 for 0.50, release 2.4 ms off, joint stalled at −26 dB): ~230 ms of free-running folds is an output error again, rugged.
  - **Nasal words now use the split march as the inverse, with only the release decided by the folds** (`sequence_recover.march_residual`, `marching.march_windows`): the march up to each release; a render forced with the waves implied so far and free from the window on dumps per step the port's waves, the port node's pressure, its outer wave and both flows (sustain `PORTWAVES=`; new `Tube::backward_wave`, `side_incoming`, `left_incoming`; defaults bit-identical); those drive every part of the march through the window; the march resumes.
    - Two faults found by testing against the engine's own dump: the window's first step reads the port pressure and outer wave of the steps before it (copied from the window on, they read zero: 44.4 dB after a window where a run without one gave 289.6); and the peel below the port reads port_node steps ahead, so its last waves before a window came from zeros and the render forced with them lost the window (−35.5 dB at the truth): the folds now run free from port_node + 2 steps before each window.
    - **/aa m aa/ at the truth: glottal flow 223.3 / 142.1 / 167.8 / 150.3 dB (before, during, the release window, after); the word's error −144.3 dB** (velum 0.45 for 0.5: −24.7; closure leaving 1 ms late: −28.4; last vowel 5% wider: −13.2). MEASURED.
  - /aa b aa/ second run: the release still 7 ms late. The seal broke when the lips, opening toward the next vowel's area, crossed the seal area: the moment depended on that area (the next vowel still at its hint) as much as on the leaving time. Now a sealed phone opens to 0.05 cm² 1 ms after its leaving (a keyframe), then moves to the next phone: the release's moment is the leaving time's. Rerunning both words.
  - Third run (release keyframe): the closure's arrival came out 23 ms early, and the vowel after it got a window too short for the node search (the LM stopped on fewer residuals than variables); /aa m aa/ likewise (closure 28 ms early, velum 0.01). Cause: a third of the coarse timing grid scored NaN (hypotheses closing early: the peel through the wrong narrowing lips overflowed) and `numpy.argmin` returns the first NaN. Non-finite residual steps now count as large (`tube_recover.fixed_clip`, as a short render's missing steps), and the grids take non-finite costs as infinite. The coarse grid then bottoms out at the truth (2.2e-15 against 2.9e5 at the next point). The second run's 7 ms may have had the same cause; this run tests the release keyframe and the guard together. MEASURED.
  - **A caveat on every sequence result so far**: the hinted times were 8 ms off the truth, a multiple of the timing grids' steps (4, 1, 0.2 ms), so the grids held the exact truth. "Transition times exact (≤ 0.01 ms)" in known answer 6 and the closure times above owe that to the lattice. The hints are now 7.3 ms off (0.1 ms off the finest lattice); the joint LM must finish the times. OPEN until rerun.
  - **Closure times near the truth were a spike** (whole word, the closure's arrival / the release late by 0.05, 0.1, 0.2, 0.5, 1, 2 ms): /aa d aa/ (alveolar seal) arrival −23.6 … −21.2 dB and release −7.4, 6.6, 27.1, 14.0, 29.0, 17.3 dB (truth −209.5): no basin, only the exact point. The peel through a narrowing but open node amplifies a wrong hypothesis's mismatch, and the folds, forced with it just before the window, carried it on. **Hidden from 0.1 cm² rather than from the seal** (`observer.closures(..., hide)`, `sequence_recover.HIDE_AREA`): /aa d aa/ arrival −40.6 … −27.0, release −42.0, −36.6, −30.6, −20.8, −11.0, −9.7 dB; /aa b aa/ alike; at 0.3 the same, at 1.0 worse (the release −19 dB at 0.05 ms). Truths unchanged: /aa b aa/ −193.0, /aa d aa/ −209.5, /aa g aa/ (velar, 18 nodes behind the open end) −208.6 dB. /aa n aa/ (alveolar nasal) at the truth: −98.6 dB. MEASURED.
  - **Nasal words had it worse**: /aa m aa/ (march residual) closure late 0.05 / 0.1 / 0.2 / 0.5 / 1 / 2 ms: −19.0, −14.6, +27.9, +1580, +47.3, +368.5 dB; release late 1 ms: inf (the march running away on the oral peel through the narrowing lips). Now with HIDE_AREA the folds also decide every stretch where a node in front of the port is narrower than 0.1 cm² but not shut: the closing (from when the narrowing enters the oral peel's reach to when the seal does; the sealed mode, which peels the nose, then takes over; `marching._march_run` now starts its ahead copy also when entering the sealed mode from a window) and the opening (from the release's window to past the narrowing). Truth −148.3 dB (−144.3 before). Closure late as above: −38.9, −36.2, −34.2, −32.1, −27.6, −27.3; release −42.3, −39.5, −36.1, −22.7, −20.5, −19.7. /aa n aa/ (alveolar) alike (−38.9 … −24.0; −43.3 … −18.7). MEASURED.
- **Fricatives in words** (`sequence_recover`: a phone `NAME@place~area`, the phone before's tube with that node narrowed, not sealed; the reversible jet on, gain 0.8, carried by every Target (`tube_recover.DEFAULT_JET`, `Target.jet`) and context to the peel (`reverse_core.peel_moving_glottis(jet=)`) and the renders; with the jet the first phone too is searched in its sequence, the static peel having no jet). At the truth, closure arrival / leaving late 0.05 … 2 ms: /aa s aa/ (alveolar, 0.08 cm²) −209.9 dB; −60.7 … −30.5; −56.5 … −21.9. /aa f aa/ (lips, 0.15 cm²) −195.9; −68.2 … −37.0; −74.8 … −40.6. Smooth (~6 dB a doubling) through the plain peel: a narrowing above ~0.05 cm² is no seal's ill-conditioning. A word with both a seal and a fricative needs the jet in the observer (not yet). MEASURED.
- **Voiceless phones** (sustain `FOLDSEQ=`: "time lower_rest upper_rest chink move_seconds" per line, `Cord::set_target` at that step; defaults bit-identical, cmp, plain and with REVJET/SEAL). A phone written with a trailing `!` holds the folds apart: `fold_keyframes` abducts them ABDUCT_LEAD (10 ms) before its arrival (cord.cpp's voiceless plosive or fricative posture) and adducts them from its leaving (+ VOICE_ONSET, 0), each move 20 ms: triggers tied to the times the search already estimates, no new parameters. The folds are only ever run forward, so the equation error is unchanged. At the truth (closure arrival / leaving late 0.05 … 2 ms): /aa p aa/ −189.9 dB (voicing stops in the closure: glottal flow rms 8, 5, 3 against ~100), −39.1 … −18.5, −40.6 … −10.3. MEASURED.
  - Voiceless fricatives with the jet at 0.8: through 0.08 cm² with the folds apart the render ran away 20 ms after the abduction (the reversible jet's source is the flow of the step before: an explicit feedback, unstable when strong); at 0.2 and 0.4 it runs (the noise at the lips ~10 against the vowel's ~92 rms at 0.2). **At 0.4: /aa s! aa/ (0.08 cm²) −205.0 dB, /aa f! aa/ (lips 0.15 cm²) −196.1, /aa s aa/ −210.1**, each smooth (~6 dB a doubling). JET_GAIN now 0.4. MEASURED.
- **Seals and fricatives in one word.** The observer now carries the reversible jet (`observer._forward_steps` / `_peel_states`: per link its flow at the step, its source's four filter states and the source it carried, STATE_COLUMNS 19; at a resumed step the source is the state's, that scatter being already in it). Against the engine (/aa s aa t! aa/, /aa k! aa s! aa/, by 20 ms): forward steps 176–290 dB, observer glottal flow 185–293 dB.
  - The engine's reversible jet ran away at a voiceless closure (the render stopped 2 ms into /k!/'s): each link's velocity is its own flow over the *narrowest node's* area, and with that node sealed (1e-4 cm²) the cavity in front's small flows read as huge jets. New `reversible_jet_closing` (sustain `REVJETCLOSE=`; `reverse_core._jet_source(closing=)`): the jet fades out as the narrowest node closes, none at a quarter of it, full from it; `posture_recover.JET_CLOSING` = 0.04 cm² on every render with a jet and a seal. Default and plain REVJET output bit-identical (cmp); known answer 8's constrictions (≥ 0.08 cm²) unaffected.
  - `observer.closures(hide=)` now hides only narrowings that reach the seal (a fricative's 0.08 cm² was hidden by the 0.1 threshold, and the folds ran free through it; it peels well).
  - At the truth, closure arrival / leaving late 0.05 … 2 ms: **/aa s aa t! aa/ −204.2 dB** (−63.3 … −33.0; −57.6 … −26.2), **/aa k! aa s! aa/ −198.4 dB** (−36.4 … −22.9; −46.0 … −17.5). MEASURED.
- **NATO codewords as known answers (the truth's error; vowels only in the profile: /r/ = ER, /w/ = UW, /j/ = IY, /l/ an alveolar narrowing of 0.25 cm² without a lateral branch).** A consonant's tube is now the last vowel's with its node changed (`vowel_before`; the phone before's carried a cluster's first constriction on), a word's leading consonants the first vowel's (searched with it: `PrefixContext(leading=)`). The jet's gain in words is 0.2 (0.3 ran away at /k!/'s release into /er/: the explicit jet's loop grows with the flow; a stable formulation is OPEN). A seal held across a velum closing (/n d/) goes wholly to the folds (the sealed march mode needs the port open).
  - Echo /eh k! ow/ −166.5 dB; Oscar /aa s! k! er/ −164.9; Alfa /ae l f! ah/ −170.4; India /ih n d iy ah/ −159.5; Papa /p! aa p! ah/ −171.9; Kilo /k! ih l ow/ −148.8; Delta /d eh l t! ah/ −165.3; Tango /t! ae ŋ g ow/ −155.3; Golf /g aa l f!/ −196.5. Each closure's or release's time late 0.05 … 2 ms rises monotonically from about −35 … −50 dB. MEASURED. These are the equation errors at the truth; the searches from hints are the test.
- **What a release's time does to the error after it** (/aa b aa/, true tubes, release scanned ±1.5 ms by 0.25): the truth 1.5e-11, ±0.25 ms 5e5 / 8.6e6, then a plateau rising slowly to 1.7e7 at ±1.5 ms; almost all of it after the window. With the next vowel at its hint the minimum moves to +0.5 ms, and the window's own part (the forward sound against the sound, release burst included) moves with it (the burst is shaped by the opening toward that vowel): scanning the window alone does not unbias the release. MEASURED.
  - Not fold phase: in a held /aa/ the forced folds are entrained (a nudge of 1e-3 to the forced wave decays to the 1e-11 floor within 10 ms). After the release the implied flow is right (against the truth's: 2.4e3 at the window's end, 9 by 55 ms), and the error between our folds and it starts small (4.5e4, 2.3e3 in the first 10 ms) and grows ~500× by 25–35 ms: **the voicing onset after a closure amplifies the folds' state difference** (an oscillation regrowing from small amplitude is an instability). Voicing onset encodes the release time sharply but in a narrow basin (~0.1–0.25 ms). MEASURED.
  - An observer on the forced folds (sustain `FORCEFLOW=` + `OBSERVE=gain`: the folds feel the forced wave plus gain (ρc/A0)(their flow − the implied flow); zero at the truth; defaults bit-identical): gain 0.1 changed nothing, 0.5 made the error ~10× larger on both sides. REJECTED as tried (`sequence_recover.OBSERVE_GAIN` stays 0).
- **A release with its own opening** (`sequence_recover.RELEASE_OPEN` 4 ms, `RELEASE_WIDE` 0.3 cm²: the seal opens to 0.05 cm² 1 ms after the leaving, to 0.3 cm² (or the next phone's area there if narrower) 4 ms after it, then on to the next phone): the folds' window (narrower than 0.1 cm²) now depends on the release's time alone. /aa b aa/ with the true vowel after: ±0.25 / ±0.5 ms cost 6.3e5–1.6e6 / 1.5e6–2.4e6 against ~1e7 before (a wider basin). With the vowel after at its hint, **the window's part alone has its minimum at the true release to 0.05 ms** (4.97e5 there, 5.0e5–5.1e5 at +0.05 … +0.15, 6.9e5 at −0.05), where the whole residual's is still at +0.5 ms. MEASURED.
  - So a transition out of a seal now takes the release's time first by the window alone (`Target.hidden_only`; ±4 ms by 0.25 ms, then ±0.3 by 0.025), then the next phone's arrival by the whole residual (4, 1, 0.2, 0.04 ms rungs). The batch is restarted with it (/aa b aa/, Echo, /aa m aa/, Papa).
  - First try in the search: the release came out 11 ms early. The window-only cost is not comparable between candidates: the window runs from the closure to the release, so an earlier release had a shorter window and fewer terms (and the scan reached only ±4 ms of a hint 7.3 ms off). REJECTED.
  - **The release's time by our whole engine run forward** (`sequence_recover.output_residual`, `Target.output_from`): our folds forced by the plain peel's waves up to 5 ms before the closure arrived (nothing hides the inside there), free from there; our sound against the sound over one stretch for every candidate (to 15 ms past the hinted release). /aa b aa/, release late by −10 … +2 ms: with the true vowel after, 0 at the truth, 5.7e4 / 2.3e5 at ±0.1 ms, 5.9e5 / 2.4e6 at ±0.5, 3.1e6 / 6.0e6 at ±1, rising to 6.5e7 at −10 (a shallow local dip at −7); **with the vowel after at its hint: the minimum at the truth (1986), 5.4e4 / 2.5e5 at ±0.1 ms**, the rest as with the true vowel. Nasals keep the window-only scan (the plain peel does not hold once the velum opens; the march's release window keeps one length). Batch restarted. MEASURED.
- **The jet in the nasal march, and frication as a gesture.**
  - `marching`: the reversible jet in the time-major oral peel (`_diagonal(jet=)`, per-link columns 15–20), the front's forward run (`_oral_forward(jet=)`, columns 14–19), the hand-over (`_peel_from_front`) and the peel below the port (`_peel_below(jet=)`); `jet_table` per step (narrowest node, its area, the gate). A first run divided by zero inside numba: the front runs ahead of the sound and read past the jet table (numba does not bound-check); its rows are now clamped as the geometry's. /aa m aa f! aa/ at the truth −145.3 dB (−145.2 without the jet).
  - Ungated, the jet ran the march away at the truth next to an open velum in /ao/ and /er/ (their own narrowings, /ao/'s pharynx 0.15 cm², /er/'s lips 0.225, hiss as jets): /aa f er m/ 32 dB, /aa f ao m/ 134 dB; /aa/ (0.37, a weak jet) and /ih/ (0.84, none) held. **The jet is now gated by the phone times** (sustain `JETSEQ=` "time gate", set at each 1 ms tick, `VoiceBody::set_jet_gate`; `reverse_core.jet_gate_steps`; `jet_gate` through reverse_core, observer and marching; `sequence_recover.jet_keyframes`: on through each fricative, ramped over the transitions into and out of it): triggers only, no new parameters. Defaults bit-identical (cmp). A fricative is written `~area`, a narrowing without noise (/l/) `=area`. JET_GAIN back to 0.4.
  - Uniform still failed without any jet: /uw n ih ao/ 341 dB at the truth. With the velum open and the mouth narrow in front of the port (/uw/'s lips 0.333, or 0.518 at node 26 just in front of the port at /iy/'s node count) the mouth is the weak path, and the march (which inverts the oral path inside the port–nose loop) ran away. **Such stretches are now the folds' too** (`WEAK_ORAL` 0.7 cm² with the port open; 0.5 still let /iy uw n/ run away off the truth).
  - At the truth with all of it: **Uniform /iy uw n ih f! ao er m/ −150.4 dB**, closure late 0.05 … 2 ms −54.3 … −34.2; /iy uw n ih/ −172.9; /uw n ih ao/ −140.4; Oscar −168.6; Alfa (with /l/ as `=0.25`) −170.4; /aa m aa/ −148.8. MEASURED.
  - (A debugging print inserted between an `if` and its `else` attached the `else` to itself and put the old release windows back on top of the new; caught by a truth that changed with an environment variable; fixed before any search loaded it.)
- **/aa b aa/ from hints with the release found by the forward run: round trip 25.0 dB.** The closure and release exact (0.02 / 0.00 ms off); the vowel after, 22.6 ms late and 0.43 cm² off; the prefix joint and the final joint (LM) did not move it at all. MEASURED.
  - Why: with the hint's tube for that vowel, its arrival's cost is flat (~1e7 from −16 to +24 ms; with the true tube a clean minimum at the truth), so the first grid placed it 17.5 ms late, and the vowel's window (from arrival + 20 ms) shrank to 15 ms. A late arrival spoils a window even at the hold's end (0.32–0.35 s, true tube: +5 ms −28 dB, +15 ms −13 dB) while an early one barely does (−5 ms −66 dB, −10 ms −59 dB). The polish round's coarse rung (±16 ms) then jumped the arrival 18 ms on a tube still 0.43 off.
  - Now after a release the next phone's first arrival is its hint less 5 ms (`EARLY_START`, no grid), and polish rounds re-time it with the found tube on ±3 ms rungs only. Rerunning.
  - The joints' finite-difference step was tube_recover.STEP (1e-3) for every parameter: 1 ms for the times, wider than a release's whole basin, so the LM's time columns were noise and the joint never moved. Times now step 2e-5 s (`SequenceShape.steps`, `TIME_STEP`).
- **Ten codewords at the truth with everything above** (release opening, gated jet, weak-mouth windows, /l/ as `=0.25`): Echo −152.4, Oscar −168.6, Alfa −170.4, India −143.6, Papa −165.9, Kilo −148.4, Delta −163.5, Tango −151.9, Golf −196.6, Uniform −150.4 dB; each first closure's or release's offset rises monotonically from −30 … −73 dB at 0.05 ms. The forward model and its inverse agree on all of them: the round trip from the truth is exact for the vocabulary's sounds (stops voiced and voiceless at three places, nasals at three, fricatives voiced and voiceless, /l/, clusters, leading and final consonants). The open problem is finding the truth from hints. MEASURED.
- **/aa b aa/ from hints (hint times 7.3 ms off, the vowels' hints from organs_vtl_tongue): round trip 134.4 dB.** Closure and release exact; the vowel after the release: first search −59.8 dB (0.232 cm² rms), arrival re-timed 0.28894 (true 0.29000), polish −81.3 (0.081), arrival 0.29022, polish −82.7; the prefix joint −140.0 dB against the sound's swing, the final joints −142.6, −147.8, −149.7, −151.8 dB, tubes 0.000 / 0.000 / 0.008 cm² rms, every time 0.00 ms off. MEASURED. **The first word with a stop carried back from its sound alone, end to end.** (Before the early start and the joints' time step: 25.0 dB.)
- **Echo from hints: round trip 38.5 dB.** Closure and release exact (≤ 0.01 ms); /ow/ after the release 0.33 cm² rms off (−54.4 dB after two polish rounds; /ow/ was one of the hard vowels alone, 53.4 dB); the prefix joint −68.8 dB, and the final joints did not move at all (288 evaluations each, nothing changed). Suspect: the voiceless folds' trigger (FOLDSEQ starts the move at a whole step, so a time moved by less than a step does nothing and then jumps: spikes in the joints' time columns). New sustain `FOLDKEYS=` (posture keyframes, linear between, set at each 1 ms tick; defaults bit-identical) and `fold_keyframes` as keyframes (modal → the phone's posture over FOLD_MOVE from ABDUCT_LEAD before its arrival, back over FOLD_MOVE from its leaving): smooth in the phone times. Truths: Echo −150.4, Papa −165.9, Oscar −164.4 dB. Echo queued again. MEASURED / OPEN.
- **/aa m aa/ from hints: round trip 27.3 dB.** Closure, release and the vowel after's arrival right (≤ 0.02, 0.7 ms); the velum grid chose 0.84 cm² for 0.50 (its window ran past the release, where the phones after were still the sealed /m/'s tube), and the prefix joint traded every tube 0.1 cm² (the first vowel's too, exact before) for the velum (0.529). Scanned with the vowel after at its hint: over the window to 30 ms past the release the cost was flat from 0.45 to 0.55 (6.1e5); **over the /m/'s hold alone (the march peels the nose there), an exact minimum at 0.50 (2.8e-17)**, 0.45 / 0.55 at 4.6e3 / 8.0e3. `velum_grid` now scores the hold, with a third rung (0.0005 cm²). The test's 0.5 also sat on the grid's 0.05 lattice: rerun with 0.47. MEASURED.
- Papa crashed in its first release scan: for a word-leading closure the forward run's stretch began 5 ms before time 0. It now starts from rest at 0 (exact there too). Requeued after Echo.
- **Echo /eh k! ow/ from hints, smooth voiceless folds: round trip 66.4 dB** (38.5 with the trigger folds). Every time exact (≤ 0.00 ms but the closure's leaving, found 1.7 ms early and corrected by the joint); /eh/ and the /k/ 0.006 cm² rms off, /ow/ 0.129; prefix joint −85.2 dB, the final joints −83.4 (the first moved the /ow/ tube, the rest stalled). /ow/ was the hardest vowel alone too (53.4–67.7 dB, the zigzag near-null space): what is left is that vowel's floor, not the consonant. MEASURED.
- **/aa m aa/ (velum 0.47) from hints: round trip 63.1 dB** (27.3 before the hold-only velum). Velum 0.474 from the grid, 0.470 after the joints (exact); every time exact; tubes 0.044 / 0.044 / 0.097 cm² rms (the prefix joint moved the first vowel, exact before, to buy the rest). MEASURED.
  - **Why the final joints never move** (here and in Echo): at the recovered point a Gauss-Newton step (and 0.1, 0.01 of it) raised the cost (0.46 → 5.1e7, 4.8e3, 2.9), only 0.001 of it lowered it (0.4595); LM with column scaling (`x_scale='jac'`, 60 evaluations a rung) did not move either. Along one bad column (the vowel after's node 11, log area) the residual is linear from −1e-2 to +1e-4 and then **jumps** (|Δr| 0.054 at +1e-4, 1.63 at +3e-4): not curvature, a discontinuity. The two residuals first differ at 0.2259 s, inside the release's window (0.2220–0.2362 s), where our folds run free: in a nasal word that window stays open ~14 ms past the release, until the lips pass WEAK_ORAL, and the free folds' response to the geometry there is discontinuous (and their phase, once off, is not pulled back after it). /aa b aa/'s window closes ~1.5 ms after its release, and its joints converged (134.4 dB). OPEN: shorter release windows for nasals, or folds that the revealed flow re-locks.
- **Papa from hints**: the leading /p!/'s release converged to the truth over the polish rounds (0.10050 → 0.10015 → 0.10000 s); the first vowel only to 1.14 cm² rms (−47.6 dB; its arrival 0.14562 for 0.15, the early start 12.3 ms before the truth). Then a crash: the prefix joint fits the first phones only, while the fold and jet schedules looked up the whole word's voiceless / fricative phones (index 2). Both schedules now skip phones beyond the arrivals they are given. The leading consonant's vowel is OPEN. Queue restarted: Golf, Alfa, Oscar, India, Papa.
- **2026-10-20: the parallel toolchain.** Owner: start children on each piece of the optimization, hunt algorithmic equivalences, build a parallel toolchain that is the same but far better (novel voices will need on-the-fly re-hunts). research/ducklips/tract_fast (Mini ~/ducklips/tract_fast; its rebuilt engine byte-identical to this tree's): four children, each owning separate files and gated on bit-identical (or reported-deviation) equivalence against this frozen reference: the engine in-process with snapshots and per-step trimming (NOTES_engine.md), the numba kernels with profiles and resume-from-step (NOTES_kernels.md), cheap exact Jacobians by AD (NOTES_gradients.md), the search driver's waste, stop rules and locality interface (NOTES_driver.md). This tree stays the reference; its only changes now are bug fixes like the one above.
- **Two corrections from the gradients child (both re-measured here).**
  - **The replica's viscous glide was a step out of phase**: Tube's constructor has already taken one glide step (a sixteenth of the law's resistance, countdown 15) when sustain's first step runs; reverse_core started it at 0 with countdown 0, so its link resistances ran a step behind GEODUMP's (up to 0.019, 29% of values in /aa iy ae/). In phase: 8.9e-16. **The truth's equation error: /aa iy ae/ −160.2 → −244.3 dB, /aa b aa/ −201.6 → −243.9, Echo −150.4 → −244.7**; /aa m aa/ unchanged at −147.8 (the nose's replica likely has its own; OPEN). `reverse_core.ENGINE_GLIDE_PHASE` (default 1; env 0 restores the old replica), in this tree and in tract_fast. MEASURED.
  - **The stalled joints' discontinuity is the pitch ear, not the voicing onset.** The engine's ear (`hear_pitch`, pitch_feedback 0.08 per closure, closures timed in whole steps) moves the folds' tension in quanta. Along /aa m aa/'s bad column (the vowel after's node 11) with PITCHFB=0 in the truth and every render, the residual is linear and symmetric from −1e-2 to +1e-2 (|Δr| 0.0057 per 1e-5) where with the ear it jumped between +1e-4 and +3e-4. The earlier attribution (the free folds in the nasal release's window) is withdrawn: REJECTED. Truths with the ear off: /aa b aa/ −242.3, Echo −245.6, /aa m aa/ −160.6 dB. The gradients child also found two smaller engine discontinuities: a triangular glottis whose rest gap goes negative (contact fraction 0 → 1) and the moving frame's fast-node switch. MEASURED.
  - Decision pending the engine child's equivalence gate (it drives this tree's render): the reversible renders will set PITCHFB=0 (the ear stays in the engine as an option; the reversible configuration already switches off what cannot be inverted or differentiated).
- **The gradients child** (tract_fast/ad, analysis/gradients.py, NOTES_gradients.md): forward-mode exact Jacobians (Dual<N> blocks, cord-only forced folds, per-node and per-phone sparsity) for the vowel-sequence residual: 143 columns in 4.5 s CPU against 43.6 s for the search's forward differences (0.30 s a reference residual), agreeing with central differences to the differences' own floor; one LM descent −19.5 → −81.6 dB with it where forward differences stalled at −19.8. Seals, the nose march and the jet: assessed, not built. MEASURED.
- **The engine child** (tract_fast: SustainRun shared by the CLI and a C library, bindings/, analysis/engine.py, NOTES_engine.md): an in-process engine with posture_recover.render's exact signature, snapshot/resume, and an Incremental renderer resuming from the last snapshot before the first changed input; per-step work trimmed. Re-run here: its gate passed (13 cases byte-identical to this tree's CLI and render; snapshots at 6 steps and 12 incremental renders identical). Renders 1.40–1.49× faster (89.7 → 61.8 ms for a 0.41 s stop word); Incremental 2.3–2.4× on a search's pattern; the subprocess and files were only 10–15 ms of a render. Its hot tube loops are now compiled with -ffp-contract=off and the reference's fused multiply-adds written as explicit std::fma (otherwise harmless edits moved every render ~1e-9). MEASURED.
- **PITCH_EAR = 0 in this tree's renders** (posture_recover: PITCHFB=0 in every render, truth and hypotheses; env PITCH_EAR restores a gain). Truths: /aa b aa/ −242.3, /aa m aa/ −160.6 dB. The paused Golf run loaded the old code (ear on) and keeps it consistently; the words after it take the ear off.
- **The driver child** (tract_fast/analysis: sequence_recover, tube_recover, posture_recover, driver_ledger, column_locality, driver_bench; NOTES_driver.md). Exact changes, each checked bitwise: leastsq in place of least_squares (scipy's wrapper computed an extra Jacobian a stage), a residual/Jacobian memo, inert columns (a stop's own areas never reach the sound: 143 → 98 evaluations a Jacobian), skipping exact repeats, renders and peels stopped at the compared stretch, Jacobian columns resumed from the base point's record, a first-vowel cache, the in-process engine (`RENDERER=incremental`; CPU per residual 0.24 s → 0.10 s). /aa b aa/ on the old model: every printed number equal to this tree's log, 50,983 evaluations in place of 63,866; 58% of the time in node-search rungs. Algorithmic options (`FAST_ALGO=halving,jointstop,stall1,jointmore`). **On the corrected model (glide in phase, ear off): /aa b aa/ 140.9 dB (8 joint rounds; tubes 0.0029 cm² rms), /aa m aa/ 131.2 dB (velum 0.4700, tubes 0.0037; the joints move now), Echo 0.8 dB (failed at its closure).** MEASURED.
  - Echo's closure on the corrected model (2-D map, leaving × arrival ±8 ms by 2): the arrival is a canyon: 1e-17 at the true pair, a plateau of 10^6.4–7 at ±2 ms of arrival whatever the leaving, a clean V along the leaving at the true arrival; with the leaving 8 ms off the arrival still has a funnel ±2.5 ms wide around the truth (10^5.3 at the bottom). The 2-D coarse grid (4 ms; the hint 7.3 ms off) never sat in it. The forward output error has the same shape. Now in tract_fast: into a seal, the arrival alone (±12 ms by 0.5), then the leaving alone, then the fine 2-D rungs. MEASURED; the rerun is going.
  - From here the working line is tract_fast (POOL=8 RENDERER=incremental FAST_ALGO=halving,jointstop,stall1,jointmore); this tree stays the yardstick. Queue: Echo, Golf, Alfa, Oscar, India, Papa.
- **Echo on tract_fast (corrected model): 48.4 dB.** The new closing scan put the /k!/'s arrival at 0.1497 s (true 0.1500); the joints took every time within 0.01 ms; /ow/ 0.143 cm² rms; the final joint did not move. At that point the times' Jacobian columns were ~2000× the tubes' (1.2e7, 8.2e6, 1.3e6 against ~5e3) and nonlinear inside their 2e-5 s step (the voiceless closure's abduction ends the voicing: a bifurcation, like the onset after a release), while the tubes' columns were linear (median disagreement 6e-4): even 0.1 of a Gauss-Newton step raised the cost (27 → 395). The tubes alone, times held: /ow/ 0.143 → 0.087 cm² but the error stayed at −70.9 dB (round trip 49.4): the times limit it. **Alternating the times alone (FD step 1e-6 s, ~25 evaluations a round) with the tubes alone: round trip 48.4 → 53.2, 56.2, 59.3 dB in three rounds, still climbing** (−70.5 → −81.5 dB). MEASURED. Now `FAST_ALGO blocks` (sequence_recover.block_rounds, after the final joint, to 12 rounds while a round gains 0.3 dB) and `polish_word.py` (block rounds from a saved answer).
- **Golf /g aa l f!/ on tract_fast: 85.7 dB.** Word-leading /g/'s release exact (0.10000); its single vowel only 0.81 cm² off after its own search, then the final joint −28.5 → −71.8 → −105.8 → −107.0 dB, tubes 0.001 cm² rms, every time exact. MEASURED. Echo and Golf being polished by block rounds.
- **Echo polished by block rounds (polish_word.py): 48.4 → 75.8 dB.** Ten rounds: −75.3, −78.4, −81.5, −84.5, −87.6, −90.6, −93.7, −95.7, −96.3, −96.4 dB; the times to below a microsecond (the /eh/'s leaving −30.2 → −0.01 µs, the release −6.0 → −0.62 µs); /eh/ and the /k!/ 0.0008 cm² rms; /ow/ 0.081 (its own floor: the zigzag near-null directions, as when it was searched alone). MEASURED.
- **Golf's polish gained nothing at a 1 µs time step** (−107.0 dB, 85.7 dB; truth −245.9). Its error is even (~83 dB by 20 ms) from the /g/'s release at 0.10 s on, ~120 dB in the final /f!/. At the point: a Gauss-Newton step, even 0.01 of it, raised the cost (0.0137 → 0.0179); the least linear columns are the /g/'s release (disagreement 0.88 between its step and a quarter of it, norm 2.8e6) and the /f!/'s leaving (norm 9.4e6); along the straight line to the truth the cost fell 0.0137 → 0.0096 over 75% of the way, then 5.8e-5 at 90%, 1.8e-16 at 100%: the last tenth is ~6 ns of time (the step to the truth: tubes 0.0017 in log area, times 60 ns at most). A word-leading release 60 ns off, and the voicing's start from rest turns it into an 83 dB floor. The block rounds' time step now shrinks with the times' moves (a tenth of the last round's largest, 1 ns to 1 µs). Re-polishing. MEASURED.
  - LM's model could not see the times' canyon at any step it was given (at 1 µs it moved them 16 ps). The block rounds' times are now a 1-D sweep (`time_sweep`: each time ±h at 1 µs, 0.3, 0.1 µs … 1 ns, every candidate of a scale evaluated together, kept when the cost falls). **Golf polished: 85.7 → 91.5 dB**, its /g/'s release 60 → 20 ns off; then no single time at any scale lowered the cost with the tube 0.0011 cm² off: what is left is coupled (time and tube together along the canyon). The floor is physical in kind: voicing that starts from rest after a word-leading release is hypersensitive to that release's moment. MEASURED / OPEN.
- **Alfa /ae l f! ah/ on tract_fast: 41.0 dB.** /ae/, /l/, /f!/ and every time right (times within a few µs); the last vowel /ah/ 0.89 cm² rms off: its first search (context: hinted times, a prefix not yet joined) reached −58.2 dB, its polishes −53.3 and −64.1, the joints −75.3. At the recovered times: the found answer −74.9 dB in /ah/'s window, −75.3 whole; the true /ah/ with the other tubes as found −83.7 in the window but −55.1 whole (the prefix joints had moved the earlier tubes 0.027 cm² off to buy the rest); at the true times the found answer −74.9: a joint local minimum, the first vowel search's failure. `revisit_vowel.py`: a vowel searched again from its hint's survey in its recovered sequence, then the joints and block rounds. Running on Alfa. MEASURED.
  - **Alfa revisited: 41.0 → 77.5 dB.** /ah/ searched again in its recovered sequence: −103.6 dB (0.123 cm² off), joints −73.8 → −100.2, −104.4, −104.5 dB; tubes 0.0015 / 0.0015 / 0.0015 / 0.087 cm² rms; times within 0.04 µs. A second search of each vowel once the word around it is known is now chained after the queue (Oscar 3, India 3,4, Papa 1,3, Echo 2). MEASURED.
- Oscar's first run crashed in /er/'s posture survey (a candidate's render gave a stretch with no period; `period_of` found no candidates). The spectrum residual now scores a silent, short or non-finite stretch as far as can be (1e3). Oscar requeued after the chain (with FAST_ALGO blocks), then its revisit. India running. MEASURED.
- **On "both ways"**: every round trip here is both directions: the sound carried back to the mouth (tubes, times, velum), then that recovered mouth driven forward by the same engine and its sound compared with the original's. The kept wave files (`_original`, `_recovered`, `_recovered_polished`, `_recovered_revisited`) are the forward direction from recovered parameters.
- **India's first pass: 5.0 dB.** /ih/ exact, /n/'s closure and velum exact (0.500), then /n/ → /d/ (one alveolar seal held while the velum closes) was scanned as a release: its leaving 17.6 ms early, and every later phone was lost. A seal-to-seal transition has no release; its times show only in the nasal murmur's fading inside the folds' window, which runs from the first seal's closing to the second's release whatever these times are, so that window's own sound is a fair comparison. With true tubes: −16.4 at the truth; the arrival a funnel (5.2–5.3 at ±1 ms up to 7.2 at ±12 ms, log10), the leaving nearly flat beyond ±1 ms (5.9 → 6.4). Now: the arrival alone (±16 ms by 1), the leaving alone, then 2-D rungs to 10 µs. India and Oscar rerun after Papa (with blocks), then the vowel revisits of Papa, Echo, India, Oscar. MEASURED.
- **Papa's first pass on tract_fast: 5.3 dB.** The first vowel after the leading /p!/ right after its joint (0.026 cm², times exact); the second /ah/ after the second /p!/ searched to −19.0 dB (2.98 cm² off); then the four-phone prefix joint, whose first rounds penalize the tubes' bends, left every tube 4.2 cm² off and a time 76 ms off (−5.9 dB) and the word was lost. A prefix joint whose result is worse on the unbent error than its start is now not taken. Papa rerun after India and Oscar; the revisits after. MEASURED.
- **India's second run: 2.4 dB.** The seal-to-seal transition now came out 2.75 / 0.25 ms off (leaving / arrival; the leaving's landscape is weak), then /d/'s release 3 ms off: its forward scan forced our folds by the plain peel to 5 ms before the /d/'s own arrival, where the /n/'s velum was still open and the plain peel does not hold. The release scan now starts before the whole chain of sealed phones the release ends, and before the velum begins to open if one of them is nasal. With true tubes: −15.3 at the true release, monotone to 5.5 / 6.1 at ±0.2 ms, 7.0 at ±1 ms, 8.0 at ±8 ms. India's third run after Papa's. MEASURED.
- **Oscar /aa s! k! er/ on tract_fast with block rounds: 72.4 dB.** Block rounds from the final joint: −67.7 (round 4) … −84.8 dB (round 11, the cap of 12), still gaining ~1.7 dB a round; /aa/, /s!/, /k!/ 0.0006 cm² rms, /er/ 0.078; the /er/'s arrival 1.6 µs off and closing. MEASURED.
- **Papa's second run: 4.6 dB.** The same as the first: the second /ah/'s survey basins all near −10 dB, then a four-phone joint that moved the second /p!/'s closure 76 ms (it scored better than its start, so the guard let it through). The renders of that search went non-finite from 0.38 s on (just after the second /p!/'s release). Free renders of the same hypotheses were finite; the forced ones were not: the posture candidates can put a node at 0.02 cm², the peel through it overflows (the ill-conditioning closures had before they were hidden), and one overflowed forced wave made the rest of the render non-finite, so a whole candidate scored 1e3. `posture_recover.sane_forcing`: a forced wave's non-finite or overflowed steps (|w| > 1e7 dyn/cm²; the lungs push ~8e3) become free (NaN), the residual's steps there count as large, the rest of the candidate stands; a wave with none is passed through untouched (bit-identical renders). Papa's third run after India's. MEASURED.
- **India's third run: 2.5 dB.** /d/'s release now 0.33 ms off; the seal-to-seal leaving still 2.75 ms; /iy/ searched to −12.3 dB (2.2 cm² off; this run predates sane_forcing); then the prefix and final joints relocated the word (every tube 23 cm² off, times 99–139 ms off) while scoring a little better than their bad start. A joint is now also not taken if it moves any time over 3 ms or any phone's tube over 1 cm² rms (`joint_sane`, prefix and final joints). India's fourth run after Papa's third, then the revisits. MEASURED.
- **Papa's third run reproduced its second exactly** (the second /ah/ −19.0 dB, 4.6 overall): non-finite forcing was not what sank /ah/. With the earlier phones true, /ah/'s window error: true tube at the true arrival −245.7 dB, hinted tube −18.9; **at the early start (12.3 ms before the truth: the hint 7.3 ms early, less EARLY_START) even the true tube −11.5 dB**: after a voiceless stop the voicing restarts there, and a wrong arrival spoils the window whatever the tube (the early start helped /aa b aa/, whose hint was late). A spectral-change estimate of the arrival from the sound failed (the held voice's frame-to-frame change never fell under any share of the transition's: every estimate ~62 ms late; REJECTED). **A quick posture fit at each candidate arrival** (9 organ parameters from the hint, 20 evaluations, on the window's equation error in the sequence): −12.7, −22.9, −24.3, −34.9, **−45.0**, −33.4, −43.7, −28.6 dB at −12, −8, −4, −2, 0, +2, +4, +8 ms from the true arrival. `arrival_by_posture`: for every vowel after a release, a 2 ms grid over ±12 ms about the hint, the deepest fit's arrival, then the vowel's search. Papa and India rerun with it, then the revisits. MEASURED.
