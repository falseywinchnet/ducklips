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
