# Surface-propagation engine (prototype, 2026-10-08)

**GIVEN** (architect): model the cavity as a deformed tube surface made of
pieces. Fluctuations travel between pieces through the air like light,
carrying amplitude and phase. Each piece receives, absorbs or blocks part
of what arrives and emits the rest according to its material. Energy is only
fed forward in time; resonance and feedback emerge. The glottis end is
plucked, and a microphone receives.

## Mechanism

- **Pieces.** Rings of pieces on a tube surface whose radius follows the area
  function, plus a glottis piece and an aperture piece.
- **Air paths.** Delay lines between pieces in neighbouring rings. A path's
  weight is the cross-section it passes through (harmonic mean of the two
  half-sections) shared over the ring-to-ring paths. Its delay is the axial
  advance of the wavefront divided by one propagation constant.
- **Scattering at a piece.** Pressure p = 2 Σ Y_k a_k / Σ Y_k from the
  arriving fluctuations a_k. The piece emits R (p − a_k) back along each
  path, where R is the material reflection (0.9995).
- **Stability.** Every piece is passive (R ≤ 1), so the network cannot gain
  energy. This is the stability the earlier boundary probe lacked.
- **Aperture.** A radiation port matched to the tube channel, returning the
  piston radiation reflection of the lip opening: −1 at low frequency,
  partial at high, with the end-correction delay. The microphone hears the
  rate of change of the radiated flow.
- **Reshaping.** Only the path weights change; the pieces and paths stay.
  It is continuous, with no remeshing and no mode switching.

## Measured (M4 Mac mini; 50 rings × 6 pieces; 192 kHz internal rate)

Resonances in Hz, compared with the exact transfer-matrix solution (rigid
walls, no side cavities):

| Shape | Exact | Surface network |
|---|---|---|
| uniform 3 cm² | 469 1407 2345 3283 | 470 1400 2330 3280 |
| /ɑ/ | 688 1359 2463 3580 | 660 1360 2520 3580 |
| /i/ | 252 2335 2757 3522 | 270 2310 2650 3600 |
| /u/ | 312 1080 2309 3741 | 310 1150 2470 3980 |

- Bandwidths are 20–190 Hz, rising with frequency through radiation.
- Rendering runs 2.5× slower than real time (unoptimized scalar code).

## Kept failures

- **Equal path weights:** shape had no effect; every shape rang as a uniform
  tube.
- **Full Euclidean delays with dense rings:** the wavefront zig-zags at about
  c/3, with resonances every 310 Hz.
- **A frequency-flat release piece at the lips:** bandwidths 230–430 Hz, far
  too damped.
- **An open-air flare of pieces instead of a radiation property:** it formed
  a second cavity. Open air has no surfaces to return energy; this agrees
  with the boundary probe's finding.
- **Paths straight across a ring:** delayed loops slowed low frequencies
  (first resonance 280 Hz instead of ~480). They are disabled below the
  cross-mode range.
- **End pieces with unit weight:** an impedance step inside each end pulled
  the first resonance down. They are now weighted by ring size.

## Not yet built

- Material propagation through the surface itself: piece eigenmodes, yielding
  walls, wall emission.
- Side cavities (sinuses) as more surface.
- Coupling to a real fold source.
- Fricative sources.
- Planned motion.
- Real-time optimization.

## Version 2 (2026-10-08): delays, wall material, coupled folds

**Delays.** Linear interpolation is replaced by integer delay plus a
first-order all-pass fraction (full strength at every frequency). The old
"speed factor" of 1.07–1.12 turned out to be the interpolator's smearing,
not mesh dispersion; with the all-pass no correction is needed. A pluck
reaches the lips in exactly length/c (MEASURED: 96 samples = 0.500 ms).

**Exact agreement (MEASURED, rigid walls, no calibration).** Resonances in
Hz, compared with the exact transfer-matrix solution:

| Shape | Exact | Surface network |
|---|---|---|
| /ɑ/ | 688 1359 2463 3580 | 680 1350 2450 3540 |
| /i/ | 252 2335 2757 3522 | 250 2330 2720 3490 |
| /u/ | 312 1080 2309 3741 | 320 1090 2320 3760 |

**Wall material.** Each piece is a mass-spring-damper per unit area
(default soft tissue: 15 kg/m², 1e4 Pa·s/m, 5.3e5 Pa/m) seen through a port
of admittance area/ρc. It is solved in the same instant as the piece
pressure.

- REJECTED: a one-sample-late wall answer. It hides a short dead-end pipe
  behind every piece and detuned the cavity.
- The radiation port is solved the same way.
- MEASURED: walls raise F1 (uniform tube 470 → 520 Hz; /ɑ/ 680 → 720 Hz).
- OPEN: on /i/ and /u/ the walls may over-damp F1.

**Folds.**

- **Coupling:** the two-mass folds are driven by lung pressure (a silent
  reservoir) minus the surface pressure just above them. Flow is solved
  against the surface's driving-point impedance every sample, so the
  coupling is two-way.
- **Turbulence:** slight glottal-jet turbulence, gated by Reynolds number and
  scaled by jet dynamic pressure, enters the paths leaving the glottis and
  buffets the folds.
- **Wobble from the body:** tension wander (0.4%, 6 Hz), flutter (0.2%,
  12 Hz) and lung-pressure wander (1.5%, 4 Hz).

MEASURED on /ɑ/: median F0 112.4 Hz for 110 asked (10th–90th percentile
111.9–112.7), LPC formants 650/1335/2490/3575 Hz. It renders about 2× slower
than real time.

## Version 3 (2026-10-08): softer closure, less noise, natural unsteadiness

Architect's verdict on v2: better, but forced, mechanical, buzzy and not
rich, with too much audible noise. Measured with `voice_quality.py` (band
levels against 0–1 kHz, harmonic-to-noise ratio per band, glottal open and
closing quotients, flow harmonic slope, jitter and shimmer).

- **Buzz: whole-length closure (MEASURED).** v2's 4–8 kHz band sat at
  −16 dB, where natural male vowels are nearer −30 dB. The two-mass folds were
  a rigid slab that shut along the whole glottal length at once. The folds now
  have 8 strips along their length that move by the half-sine vibration shape,
  and the rest opening is a posterior chink. They close progressively (zipper)
  and keep a small leak. 4–8 kHz on /ɑ/ falls to −26 dB; the flow slope goes
  from −13 to −17 dB/oct. One strip reproduces the old folds exactly.
- **Forced:** lung pressure 800 → 600 Pa (relaxed speech).
- **Noise, part 1 (MEASURED):** jet turbulence was the audible hiss.
  Buffeting of the fold surface contributed nothing measurable. Softer closure
  removed 9 dB of high harmonics that had masked it. The gain goes from 0.006
  to 0.001, which sits below the voice's own irregularity.
- **Noise, part 2 (MEASURED, a real defect).** Physiological wander (lung
  pressure, adduction) was one-pole-filtered white noise. It keeps a white
  residue at audio rates (−6 dB/oct only). Where the variable drives the flow
  directly, the microphone's d/dt turns that residue into flat hiss: adduction
  wander alone took /ɑ/ to an HNR near 0 dB in every band. `SlowNoise` is now
  three cascaded poles (−18 dB/oct), normalized in closed form.
- **Mechanical (MEASURED):** v2 cycles were too regular (jitter 0.24%,
  shimmer 0.7%; healthy voices are about 0.5–1% and 2–4%). New sources:
  muscle-tension tremor (0.4%, 40 Hz band) and arytenoid adduction wander
  (15% of the rest opening, 20 Hz band). Result: jitter 0.47%, shimmer 3.2%.
- **Kept caution: the HNR measure.** The 100 ms-frame HNR counts slow pitch
  wobble and jitter as noise. With all variation off it reads 43 dB in every
  band, so it is not a hiss meter on varied voices; compare at equal wobble.
- **Unchanged:** F0 109.6 Hz for 110 asked (`unit_pitch` 139 for 8 strips at
  600 Pa).
- **OPEN:** on /u/ the flow-peak jitter reads 1.6%, probably F1 ripple on the
  flow peaks rather than real period jitter.

## Version 4 (2026-10-08): first comparison with a human recording

The architect compared a spectrogram of their own "ah-ee" with v3's. Theirs
has steep high-frequency roll-off, deep valleys between harmonics and formants,
a breathy onset and a moving pitch; v3 was bright up to the top of the band.

- MEASURED: v3's band levels stopped falling above 8 kHz (−31…−35 dB per
  band out to 24 kHz). Of that, the jet noise was white, and the microphone
  derivative tilted it upward. Shaping the turbulence with two poles at 2 kHz
  takes 20–24 kHz from −35 to −44 dB, which lets breath noise come back
  (gain 0.004) without hiss.
- MEASURED: fold strips beyond 8 (16, 64) change the high band by less than
  1 dB, so strip kinks are not what keeps the voiced high end up.
- Not yet compared numerically against the human recording: the screenshot
  carries no level scale (v3 was peak-normalized to −2 dBFS). The next gate is
  the same measures on the architect's own recording.
- Absent by construction so far: breathy onset and offset, pitch contour, and
  any check of formant bandwidths against real ones.

## Human reference (2026-10-08): the architect's recording

Recording: 79 s, 44.1 kHz mono. NATO alphabet, digits, letters, and one
sentence. Measured with `voice_compare.py` over strongly voiced frames
(40 ms, periodicity > 0.6, within 30 dB of the loudest).

| Band (Hz) | Human | v4 ah-ee | sharp folds, ah+ee+oo | soft seal 0.6 mm² |
|---|---|---|---|---|
| 1–2k | −8.9 | −18.1 | −16.7 | −22.7 |
| 2–4k | −11.5 | −24.7 | −14.7 | −26.3 |
| 4–6k | −29.7 | −29.7 | −20.4 | −32.4 |
| 6–8k | −33.9 | −33.0 | −28.2 | −36.0 |
| 8–12k | −42.0 | −34.2 | −29.0 | −37.7 |
| 12–16k | −47.5 | −37.9 | −31.9 | −41.3 |
| 16–22k | −60.7 | −38.2 | −33.7 | −42.6 |

- **Human shape:** the voice is strong to 4 kHz, then a cliff of about 18 dB
  into 4–6 kHz, falling steadily to −61 dB. The recording chain is not the
  cause: /s/ frames peak at 4–6 kHz and keep −14 dB at 8–12 kHz.
- **Other human measures:** LPC bandwidths B1 84, B2 132, B3 247 Hz. F0
  90–142 Hz, moving 0.135 semitone per 10 ms (renders: 0.03). Of 81 voiced
  runs, 69 begin with noise before periodicity (renders: none).
- **v3/v4 were the wrong cure for buzz (MEASURED).** Progressive closure took
  energy out of 1–4 kHz, where the human voice is strong. v2-like sharp
  closure matches the mid band better. Contact softness lowers everything
  above 1 kHz together and makes no cliff.
- **The surface network is faithful (MEASURED).** Network transfer from
  glottal flow to microphone against the exact lossless plane-wave tube with
  Flanagan radiation on /ɑ/: within 1 dB to 6 kHz, and 3–19 dB *below* exact
  above 8 kHz. The exact 1-D tube itself rises toward high frequency (with
  the radiation derivative). So the missing cliff is physics the plane-wave
  tube lacks, not mesh error.
- **The floor is not noise (MEASURED).** Breath noise, buffeting and wobble
  change nothing above 8 kHz; the glottal flow itself falls to −95 dB.
- **HYPOTHESES for the cliff, untested:**
  - hypopharyngeal side cavities (piriform fossae: an antiresonance near
    4–5 kHz);
  - three-dimensional cross-modes and losses above about 4 kHz (no
    across-ring propagation exists yet);
  - a glottal source that is weaker above 4 kHz than the two-mass model.

## Version 5 (2026-10-08): the high-frequency floor was a corner in the flow area

- **Piriform fossae: REJECTED as the cause of the cliff.** Two closed
  2 cm side branches on the exact tube move 4–22 kHz by only a few dB.
  Under a −12 dB/oct source, the plain exact tube falls only about 5 dB by
  16–22 kHz, where the human falls about 60.
- **MEASURED cause.** The flow area was the hard minimum of the lower-mass
  and upper-mass openings. Each cycle, as the narrowest point passes from one
  mass to the other, the area has a corner, which sprays harmonics to the
  top of the band. Contact smoothing could not touch it. On a real fold the
  narrowest point slides along a rounded medial surface.
- **Fix.** A smooth minimum over `medial_roundness` (1 mm²), plus slight
  sealing softness (0.1 mm²). With wobble and noise off, 16–22 kHz goes from
  −43 dB to −76…−82 dB, and 0–2 kHz is unchanged.
- **Breath noise.** Shaped at 2 kHz with gain 0.012, matched to the
  recording's top bands.
- **Irregularity.** Tremor (0.2%) and adduction wander (8%) were halved to
  match the human mid-band HNR.
- **Breathy onset and offset:** the folds start apart (30 mm²) and close in
  over 80 ms.
- **Strips:** back to 1, because progressive closure removed mid-band energy.

MEASURED, /ɑ/ + /i/ + /u/ pooled, against the human recording:

| Band (Hz) | Level, human | Level, v5 | HNR, human | HNR, v5 |
|---|---|---|---|---|
| 0–1k | 0.0 | 0.0 | 15.1 | 15.7 |
| 1–2k | −8.9 | −16.4 | 10.4 | 9.9 |
| 2–4k | −11.5 | −18.5 | 6.3 | 5.7 |
| 4–6k | −29.7 | −28.8 | 2.6 | 2.6 |
| 6–8k | −33.9 | −41.8 | 1.4 | 0.9 |
| 8–12k | −42.0 | −48.4 | 0.2 | 0.3 |
| 12–16k | −47.5 | −53.2 | −0.1 | 0.1 |
| 16–22k | −60.7 | −54.0 | 0.1 | 0.0 |

- /ɑ/ flow: open quotient 0.62, jitter 0.60%, shimmer 2.3%, F0 110 Hz.
- **OPEN:**
  - The mid band is still about 7 dB weak, though the human set includes
    vowels the renders do not, so a vowel-matched reference is needed.
  - Pitch motion is a quarter of the human's (planner work).
  - The top band is 7 dB high: the noise shaping needs a third pole.

## Version 6 (2026-10-08): tones, branches, speech

**Two "modem" tones, both MEASURED.**

- A 1 kHz comb during glides: the cavity was reshaped in steps every 1 ms.
  Path weights, wall admittances, station areas, loss filters and radiation
  filters now all glide every sample.
- A 3.8 kHz whistle on /u/: the exact lossless /u/ tube has resonances
  3–7 Hz wide. Boundary-layer losses (viscous and thermal, ∝ √f, wetted
  perimeter 1.5× a circle) now act on every path as a passive first-order
  shelf fitted to √f above 3 kHz. They widen /u/ F4 from 20 to 50 Hz.
- Walls yield over 0.67 of that perimeter (palate, teeth and cartilage do
  not). This keeps F1 where the circle-wall model put it (/u/ 360 Hz).

**Tubes.** The surface is now a set of tubes:

- the airway;
- the nose, opening from the velar port, radiating at the nostrils, with
  turbinate perimeter 3× and 30% yielding;
- the piriform fossae, soft-walled, beside the larynx;
- the maxillary sinuses, bony, from the nose.

Ring to ring paths only; one radiating opening per open tube.

**Kinetic losses (Borda–Carnot).** At every boundary where the flow widens,
a series resistance ρ|U|(1/A_n − 1/A_w)²/2 sits on the paths. At an opening
it is ρ|U|/(2A²). In wave units it is about half the jet Mach number of the
path admittance. Pressure now builds behind closures as in the body
(MEASURED in "duck lips": /d/ 570 Pa, /k/ 870 Pa, /p/ 740 Pa).

- REJECTED: one lumped jet-exit resistance at the narrowest section. When
  the narrowest point moved, a 50 Pa drop relocated within one sample: a
  click.
- Also fixed (MEASURED by per-sample trace):
  - The losses read target areas while the paths were still gliding, so a
    closing lip's resistance jumped 45× at a block edge (lip pressure
    −92 → +1390 Pa in one sample). Station areas now glide.
  - A vanishing opening's radiation filter tends to a pure short, where the
    instantaneous solve divides by 1 + b0 → 0. Its radius is floored at
    k ≥ 1 (0.34 mm).
  - The jet speed used by the losses is bounded at 60 m/s and smoothed at
    1 kHz.

**Speech (`surface_synthesizer.*`, `surface_say`).**

- The modal synthesizer's behaviour layer (cues, articulators, tongue,
  planner, pitch control) drives the surface at 4× oversampling, with one
  articulator block per 1 ms gliding per sample.
- 225 pieces and 1216 paths (airway rings of 3, branch rings of 1). Real-time
  factor 1.2–1.3 on the M4 for one voice.

**Frication, matched to the architect's recording.** Human sibilants sit
4.6 dB under vowels, with 68% of their energy at 4–6 kHz.

- Noise is a series pressure where the jet leaves the narrowest section, or
  at the teeth for sibilants (obstacle).
- The source corner follows Strouhal 0.2 of jet speed over width
  (HYPOTHESIS), with three poles (two left 16–22 kHz 19 dB too strong).
- /s/ place moved from 2.3 to 1.8 cm from the lips: the front-cavity peak
  went from 3.3 kHz to 4–6 kHz.
- MEASURED, baritone sentences: sibilants −6.7 dB re vowels, 61% of their
  energy at 4–6 kHz.
- A reversed (inhaled) glottal jet radiates 0.3 as much turbulence.

**Pitch controller.** It counted slow wobbles of the open glottis during
breathing as cycles near 66 Hz and wound fold tension to its limit, and the
alto never phonated. It now accepts only cycles within 0.6–1.6 of intended
and relaxes after 30 ms without voicing.

**MEASURED, baritone connected speech against the human recording.**

| Measure | Human | Baritone |
|---|---|---|
| Level 1–2 kHz | −8.9 dB | −12.4 dB |
| Level 2–4 kHz | −11.5 dB | −15.5 dB |
| Level 16–22 kHz | −60.7 dB | −46 dB |
| LPC B2 / B3 | 132 / 247 Hz | 75 / 114 Hz |
| F0 motion | 0.135 st/10 ms | 0.075 st/10 ms |

**OPEN:**

- The planner abducts the glottis 30 ms into a phrase-final vowel ("hello").
  The modal engine does the same.
- Fixed station spacing ignores section-length changes (lip protrusion).
- Alto and child sibilants are 6–8 dB weaker than the baritone's.

## Version 7 (2026-10-08): architect's notes on the first sentences

**Architect's notes on the v6 checkpoint:**

- "my name is ucklips": the /d/ was broken.
- The middle sentence was unintelligible.
- "sea??ore".
- The child gargles.
- Whisper clips.

**/z/ and /d/ in "is duck" (MEASURED, probe and zoomed spectrogram).**

- /z/ made no noise. With the glottis fully adducted the jet through the
  tongue gap stayed near 10 m/s, under the turbulence threshold. The planner
  now parts the folds by 0.3 of the way to abduction for voiced fricatives
  (jet 150 cm³/s, 300 Pa behind it).
- /d/ opened from closure to 0.5 cm² in under 4 ms, so the burst lasted
  under 1 ms. Stops now release at 0.3 of the closing rate.

**Breath (planner).** The next phrase's inhale was timed 0.50 s before it,
whatever the pause. After a comma (0.22 s) the inhale cues landed inside
the previous phrase's final vowel ("hello"): the folds blew open and the
output dimmed. An inhale is now scheduled only when it fits (≥ 0.30 s);
otherwise pressure falls and returns. The phrase-end pressure fall is faster
(rate 30), so whisper noise does not trail into pauses.

**Lung pressure.** The presets drop by about a quarter (baritone 800 →
600 Pa). Conversational values are 500–700 Pa; 800 sounded pressed.
Accent rises: baritone ×1.5 (to 0.27); alto 0.30, child 0.25. Baritone
pitch motion is now 0.10 st/10 ms (human 0.135).

**Pitch controller.**

- It ignores cycles off by more than 0.6–1.6× from intended and relaxes
  after 30 ms without voicing. It had been winding up to 3.7× tension during
  breathing, and the alto never phonated.
- Correction is now per second (30/s) rather than per cycle.

**The child's gargle (MEASURED cycle by cycle).**

- With all unsteadiness off, cycles are perfectly regular; the earlier
  half-harmonic reading was intonation smear.
- Every preset had adult-length folds (14 mm), so the child's peak flow was
  about 260 cm³/s, roughly twice a child's. Fold length now scales with the
  body (alto 0.79 → 11 mm; child 0.57 → 8 mm), with stiffness scaled alike
  so threshold pressure is kept. Child peak flow is now 134 cm³/s; its
  resting gap is halved.
- Accent rises of 0.37 made a fast up-and-down glide on every syllable;
  the child is back at 0.25.
- The slow tension wander is scaled to 0.4 of the preset value.

**Whisper.**

- The /d/ release spike (3.6× the 99.9th percentile) is gone with the
  slower release.
- Small transients at contacts remain: a 0.3 cm³/s kink at /l/ onset is
  physical, but without voicing to mask it, it stands out.

**Fricatives, by place.**

- /f/ was silent: a series pressure in a nearly closed lip opening drives
  almost no flow. A jet leaving the lips now also radiates its turbulence
  directly into free air, as q·A/(ρc) × 3 (HYPOTHESIS, calibrated).
- MEASURED re vowel: /s/ −7.4 dB, /th/ −26 dB, /f/ −31 dB before the
  free-jet ×3.
- The frication source scales with body size (area relative to 0.1 cm² ×
  the body's area scale). Alto and child sibilants are now −11 dB against
  the baritone's −8.6.

**Lexicon.** "lazy" came out as L AE Z IY and "jumped" as JH AH M P EH D
(letter-to-sound guesses). Added the test words, the NATO alphabet and the
digits, so the architect's recording can be compared word for word.

**OPEN:**

- An utterance-initial voiced fricative (the "th" of "the") is lost.
- Unstressed ER does not reach its F3 dip in time.
- /r/ colouring overall is untested.

**Speed (MEASURED, M4, one voice).** One piece per ring renders identically
to three (difference −188 dB). With ring-to-ring paths only, every piece of a
ring carries the same plane wave. It is 5× faster: 6.1 s of speech in 1.5 s,
a real-time factor of 0.25 (97 pieces, 192 paths). The default is now 1;
more pieces return with across-ring paths. Downsampling and compiler flags
were not the bottleneck (no change).

**Utterance onset (planner).** Lungs pressurize 0.22 s before the phrase
and the glottis adducts 0.06 s later, so the rising pressure does not blow
1,100 cm³/s of audible breath through an open glottis. The first consonant
of a phrase is formed 0.08 s early, which made the "th" of "the" voiced and
present.

**Consonant–vowel joins (planner, from labelled spectrograms).**

- Voicing began behind the lips with the velum still open after the
  inhale, a hum before "by". The velum now seals with the pressure rise.
- After a voiceless fricative the glottis started closing only at the
  fricative's end and took about 30 ms. The tongue had already opened, so
  600 cm³/s of breath passed: "s-h-ee". Closing now starts 30 ms before the
  end. The fricative releases 20 ms before its end at 170 rad/s (was 5 ms
  at 110, which left about 50 ms of muffled voicing).
- MEASURED: "sea shore" now runs noise → voice with no gap.

## Version 8 (2026-10-08): the architect's reading of the three lines

**Reference.** `out/ref/three.wav` (from `~/Desktop/Untitled_three.wav`) is
the architect reading the three test lines, spoken (0–9.6 s) and whispered
(11–21.3 s).

**Architect's verdict on v7.**

- Baritone: "??? my name is ??lips. ?? ?? ??? dog. she sells seashells by
  the seashore."
- Child: "heh-hello. my name is ducklips. … ze fox flumped over the sleeping
  dog. she sells zeashells by the zeashore."
- Whisper: severe transients.

**Whisper transients (MEASURED).**

- Human whisper: 2 transients in 10.3 s. Ours: 56 in 8.7 s (17 in sentence
  1).
- They survived removing frication and aspiration, so the noise sources
  were not the cause. They sat at every closure and release, in the lip
  radiation.
- Cause: station areas glided linearly between 1 ms articulator updates,
  so the rate of closure jumped at every block edge. The microphone hears
  the rate of change of flow, and whisper's large steady airflow turned
  each jump into a click.
- Fix: areas glide linearly, then pass a two-pole smoothing at 600 Hz. Every
  admittance (paths, walls, ports, radiation) is recomputed from the
  smoothed areas each sample, so motion has continuous velocity. Loss
  filters still glide linearly.
- Result: sentence 1 went from 17 transients to 1; the full set has 2 in
  8.7 s. Vowel resonances are unchanged.

**Planner timing, from labelled spectrograms against the reading.**

- "Heh-hello": the early pressure rise also closed the glottis before a
  voiceless first sound, so a voiced blip preceded the /h/.
  - Now: pressure rises 0.12 s before the phrase.
  - Voiced first sound: the glottis closes with the rise (lung rate 60), so
    there is neither "uh" nor an audible breath.
  - Voiceless first sound: the glottis stays open.
- "Zeashells": voicing ran into /s/ after a voiced sound. Voiceless
  fricatives now part the folds 35 ms before their slot; stops still part
  them 12 ms before.
- "Flumped": the voiced affricate /j/ had a fully closed glottis and no
  noise. Voiced affricates now part the folds like voiced fricatives.
- Choppiness: stop closures ran 70–90 ms and ate the vowels (the gesture
  started 20 ms early and the slow release ran past the slot). Closures now
  start 10 ms before the slot and release 12 ms before its end.
- The reading says "the lazy fox jumped over the sleeping dog" in about
  1.9 s; the planner takes 2.6 s. Speaking rate is OPEN.

**OPEN:**

- Baritone brightness: the 3–5 kHz formant bands are as strong as F1–F2,
  and the resonances are narrower than the human's. The plane-wave tube
  lacks the losses of 3-D cross-modes.
- Speaking rate.
- /r/.

## Version 9 (2026-10-08): bursts, fricative edges, onsets

**Architect's verdict on v8.**

- Baritone more legible, but: "hello" hesitant; the /d/ of "duck" not
  enunciated; "the" not properly said; transients in "jumped", "seashells"
  and "seashore"; "dog" sounds like "lock".
- The child sounds like a woman.
- Whisper has no /d/ at all.

**Fricative edges (MEASURED).** Every /s/ switched its noise on and off like
a rectangle. The source jumped to the teeth, with a gain step, once the
tongue was 30% of the way into the gesture. Now the source sits at the teeth
for the whole gesture, the obstacle gain blends in with progress, and gain
and jet area are interpolated per sample. /s/ rises in about 30 ms, and the
three sentences show no transients.

**Jet speed.** It now comes from the pressure drop across the constriction,
√(2Δp/ρ) between sections two either side (smoothed at 1 kHz).

- Flow over area read 5 m/s at a /d/ release with 280 Pa behind it: through
  a gap that has just opened, the flow is still building against inertia.
- For a closing gap the smaller of the two is used, or every closure
  flashed like a release.
- The obstacle (teeth) gain doubled (`obstacle_boost` 1.6) to keep
  sibilants at about −7 dB re vowels. /f/ −19 dB, /th/ −21.5 dB.

**/d/ and /t/ bursts.** Compared with the human reading: the human /d/ has
no voice bar, then a sharp burst strongest at 4–6 kHz and 20–40 ms of noise.
Ours had a voice bar and no burst.

- Alveolar stops take the teeth as an obstacle, like /s/.
- Two-stage release: the pressed tissue first eases to bare contact, then
  the gap opens from rest at 0.35 of the closing rate. Released in one
  motion, the tongue crossed contact at full speed and opened to 1.6 cm² in
  6 ms. A slow single-stage release (0.3) made /d/ sound like /l/.
- A voiced stop before a vowel half-opens the glottis from 15 ms before its
  closure and closes it at release (English short-lag /b d g/). MEASURED:
  pressure behind "dog"'s /d/ is 410 Pa (was about 250 with voicing).
- A nasal before an obstruent seals the velum by its end. In "sleeping
  dog" the /d/ closure had been a nasal.
- MEASURED burst noise: "duck" peaks near 20 Pa over about 10 ms. "Dog" is
  still weaker (6 Pa, 6 ms) — OPEN.

**Onsets.**

- A voiceless first sound gets pressure 50 ms ahead with the glottis open;
  120 ms had made a long "hhh" (heard as hesitancy).
- A voiced first sound closes the glottis as pressure rises, 120 ms ahead.
- /th/ and /dh/ narrow to 0.09 cm².

**Child.** The tract went from 12.2 to 10.8 cm (nose 7.2 cm), about 27%
shorter than the alto's; at 12.2 cm it was heard as a woman.

## Version 10 (2026-10-08): "they sound Irish"; whisper, gravel

**Architect's verdict on v9.**

- They sound Irish.
- Whisper still clips "seashells" and "seashore".
- Gravel says "bucklips", and sounds queasy rather than gravelly.

**Intonation ("Irish").** Every stressed syllable carried its own Gaussian
rise and fall (σ 75 ms): a word-by-word lilt. Replaced with the "hat"
pattern: pitch rises into the first accent, holds across the following
accents (heights interpolated), and falls after the last, nuclear accent.
MEASURED baritone contour: rise to about 150 Hz, declining plateau, final
fall to about 90 Hz.

**Whisper "clipping" (MEASURED against the whispered reading).**

- Human whispered vowels sit 12–15 dB under the fricatives, with noise up to
  about 6 kHz. Ours sat 25–30 dB under and were dark, so each fricative was a
  loud slab between near-silences.
- A non-phonating glottis now shapes its turbulence at 4.5 kHz (was 2 kHz);
  the vowel gap is 12.0 dB.
- Fricatives also fell off a cliff at their end. They now release 25 ms
  before the end at 120 rad/s, so the noise fades as the gap widens.

**Gravel.**

- The "queasy" sound was the slow wander (2.5%) and tremor (1%): a slow
  wobble. These are reduced to 0.8% and 0.3%.
- The folds can now differ: the right fold's tension relative to the left
  (`FoldParameters::asymmetry`).
  - REJECTED for roughness: at asymmetries 0.35–0.7 and 400–1000 Pa the two
    folds stayed entrained (no subharmonics; jitter under 1%). Pressed
    adduction (negative rest area) stops phonation; this two-mass model
    cannot make vocal fry.
- Added `FoldParameters::roughness`: irregular fold tissue as fractional
  stiffness and adduction fluctuation in a 150 Hz band (cycle scale).
  Gravel 0.03 gives MEASURED jitter 2.8% and shimmer 3.3% (rough-voice
  range).

**"Bucklips" (cue conflict).** A voiced fricative closed the glottis 10 ms
before its end, ready for a vowel. That cue landed after the following
voiced stop's opening cue and cancelled it, so "is duck" voiced straight
through its /d/ with no pressure and no burst. Before an obstruent the
fricative now leaves the glottis alone. MEASURED gravel /d/: devoiced
closure at 375–480 Pa, then a burst. Baritone and child now burst too.

## Version 11 (2026-10-08): the click, the Irish /r/

**Architect's verdict on v10.** "ducklips" has a click; the voices still
sound Irish, "shore" especially.

**The click (MEASURED).**

- The /d/ burst peaked at 1.3× the vowel; in the human reading it is about
  0.2–0.45×, and in rms 6–10 dB under the vowel. With frication off it fell
  17×, so it was the stop's noise at the sibilant's obstacle gain.
- Stop releases now use `burst_scale` 0.35 of the obstacle noise. A gesture
  counts as a closure below 0.01 cm² target area, because a gliding target
  never reaches exactly 0.
- Voiced voices: /d/ 5 ms rms is now 0.2–0.3× the vowel.

**Whisper /d/ thump.**

- Between /z/ and /d/ the tongue tip released and closed again within
  10 ms, cutting a whisper's 480 cm³/s airflow off at once (3.3× the
  vowel).
- The planner now drops a release when the same articulator's next gesture
  starts within 30 ms (the tongue slides on).
- Stop closures approach at 0.7 of the slot's speed.
- The surface smooths motion near contact: from 600 Hz down to 150 Hz as a
  section closes below 0.2 cm² (soft tissue meets gradually).
- The whisper /d/ is now silence, burst, then vowel, at 2.3× the whispered
  vowel (human 1.1–1.3×).

**/r/ (Irish).**

- The "shore" /r/ came late and only reached F3 2.1 kHz.
- After a vowel, /r/ now starts 80 ms early at 110 rad/s, colouring the
  vowel as American post-vocalic /r/ does.
- Spec swept (MEASURED, minimum F3 in "shore" / "aa r aa"):

| Place, area, lips | "shore" | "aa r aa" |
|---|---|---|
| 3.6 cm, 0.30 cm², 1.1 cm² (old) | 2065 | 1670 |
| 3.6, 0.20, 1.1 | 1922 | 1610 |
| 3.6, 0.20, 0.6 | 1811 | 1576 |
| 4.5, 0.20, 0.6 (chosen) | 1592 | 1622 |
| 3.0, 0.15, 0.6 | 1965 | 1564 |

- The inhale is not voiced: the flow is a smooth −420 to −640 cm³/s.
  Breath noise had fooled the periodicity detector.

**"Seashore" still Irish (v11 verdict): compared with the human "by the
seashore" (MEASURED, LPC tracks).**

- Human "shore": F1 about 500 Hz throughout. F2 dips to about 900 Hz and
  rises toward 1300. F3 falls steadily from about 2200 to 1600 Hz across the
  word, and pitch falls within the word.
- Ours: a static open "aw" (F1 600, F2 1100, F3 2300 flat for 250 ms), then
  a light /r/ at the end. The /r/ gesture began only 80 ms before its slot.
- Post-vocalic /r/ now begins at the middle of the preceding vowel (rate
  60), and "aw" before /r/ takes the "o" posture (American raising).
- Now: F1 480–540, F3 2400 → 1570 by mid-vowel, reaching about 1470.
- OPEN: F2 drifts down to 890 Hz where the human's rises to 1300. Probably
  too much lip rounding for /r/.

**"Still Irish; your shore is exactly the same" (v12 verdict).** The formant
change was in the file (F3 fell to 1500 Hz) but did not change what was
heard, so formants were not the cue.

Compared with the human reading's pitch track (MEASURED): every stressed
word starts high and falls within itself, with a reset at the next stressed
word. "sells" 128 → 105 Hz, "seashells" 128 → 95, "sea" 140 → 105, "shore"
120 → 80. Ours was a smooth, nearly level slope (the hat had overcorrected
into monotone); before that, bell-shaped accents peaked at 55% of the vowel
(late peaks, a mark of Irish English).

- Accents are now falling (H*): a quick rise (30 ms) into the stressed
  vowel's onset, then decay over about 130 ms; the highest accent wins.
  Accent heights fade 0.92 per accent (was 0.75); the pitch tracker runs at
  120 rad/s (was 60).
- Compounds: SEASHORE, SEASHELLS, SEASHELL with stress on the first part
  (the human said "SEAshore"). The test text now uses the compounds.
- The post-vocalic /r/ approaches at 40 rad/s. At 60 its F3 dropped 800 Hz
  in about 40 ms, more like a tap than the human's 250 ms glide.
- MEASURED: the F0 contour now shows a fall on each stressed word and a
  range of 85–150 Hz (human 85–140).

**"Still the same; you're still getting shore wrong" (v13 verdict).**
Segment-by-segment comparison with the human "shore" (MEASURED):

- A gap between "sh" and the vowel. The tongue released before the glottis
  closed, and the open tract passed a weak breath for 15–20 ms (−28 dB):
  "sh-'ore". The glottis now closes during the fricative's last 40 ms and
  the tongue releases 5 ms before its end (rate 150), so voicing starts under
  the constriction as in the human recording.
- An "uh" after the /r/. The phrase-final /r/ released while voicing went
  on (F3 back to 2.4 kHz), like a released Scottish or Irish r. A
  phrase-final /r/ now holds until after the voice ends.
- The voice faded late and steeply. Lungs now fall from 80 ms before the
  phrase end at rate 15.
- The "sh" was dark: noise peak 1.9 kHz against the human 3.1 kHz (centroid
  2.1 against 3.3). Place sweep: 3.1 cm → 1.9 kHz, 2.7 → 2.4 kHz,
  2.4 → 2.9 kHz. Postalveolars (/sh zh ch jh/) now sit at 2.4 cm.
- Process failure, caught before sending: one render this round came from a
  binary older than the source (the build output was filtered away). The
  render script now prints the build tail and both timestamps, and every
  send is preceded by a measurement on the sent file.
- OPEN: the human "shore" vowel has almost nothing above 2.5 kHz; ours keeps
  bright upper formants at 3.5–5 kHz (the baritone brightness item).

## Version 12 (2026-10-08): the three-dimensional tube (first build)

**Architect.** "These all have the same problem … you do not understand the
spectrogram." Then: "The cause: you didn't build the 3d tube."

**The common problem (MEASURED in every render, against the human
spectrograms).**

- Fixed, sharp, bright bands at about 3.5, 4.2 and 5 kHz, the same in every
  vowel and as strong as F1–F2. The human vowels are dim and diffuse above
  about 3 kHz.
- The "surface" so far was a 1-D tube: rings of pieces with only
  ring-to-ring paths carry plane waves, and it matches exact 1-D theory.
  The 3-D cavity was never built. Every per-sound fix since was treating
  symptoms of this.

**`tube3d.*`: the air as a 3-D staggered FDTD grid.**

- Pressure lives at cell centres and velocity on faces. Cells are 3.2 mm at
  192 kHz (c·dt·√3 ≤ cell).
- The tube is bent: a vertical pharynx (0.37 L), a quarter-circle bend
  (outer radius 4 cm, length measured at 3 cm), then the horizontal mouth.
- Cross-sections are ellipses lying against the outer wall (back of the
  pharynx, the palate). Width follows an anatomical profile (1.2 cm at the
  larynx, 3 cm pharynx, 3.4 cm mouth, 4 cm lips); height is set by the area.
  Constrictions are slits along the outer wall.
- Cells are partly open from a smooth distance to the ellipse, scaled per
  slice to the exact area. Faces open as the less open neighbour. Open
  fractions glide over each block.
- Yielding tissue (mass, damping, stiffness per area) on each cell's
  exposed face area.
- The glottis drives the first slice; the lips radiate through an implicit
  lumped Flanagan port over the last slice.
- Only active cells are stepped. `tube3d_probe` measures the impulse
  response; `tube3d_voice` sustains a vowel with the two-mass folds.

**MEASURED, first build.**

- Stable; impulse responses decay.
- Uniform 3 cm² tube: 520/1380/2190/3060 Hz (exact 470/1407/2345/3283). The
  higher formants run 7% low: the length mapping through the bend still
  needs calibrating.
- /ɑ/: 690/1250/2280/3220/4120 Hz (1-D surface 710/1360/2450/3540/4380).
- The bend makes cross-mode resonances at 5.9, 7.0 and 7.7 kHz, but they
  are about 30 Hz wide (whistles). Friction on lateral motion (2π·300/s)
  did not damp them, so they are not lateral modes.
- Sustained vowels show a hard cliff near 4.5–5.3 kHz with little above
  (4–6 kHz −36 dB, 6–8 kHz −62 dB, against the 1-D −28/−41). Too abrupt to
  be physical. Suspected cause: at 3.2 mm with partly open cells the grid
  cannot carry sound that high through narrow passages. Published 3-D tract
  FDTD uses 1–2 mm cells.
- Voicing starts after about 0.25 s, against about 0.05 s in 1-D.
  Suspected cause: the glottal coupling through ρc/A and a one-step-late
  pressure.
- Speed: 8.6× slower than real time (1 s of vowel in 8.6 s on one M4 core).
  Halving the cell costs 16× (cells 8×, time step 2×).

**OPEN:**

- Validate the grid against analytic duct modes; find the cliff.
- Fix the whistle modes and the slow onset.
- Calibrate the length.
- Add the nose and side cavities.
- Couple the tube to the synthesizer.
- Decide on resolution against speed (GPU?).

**REJECTED (2026-10-09), by the architect:** a grid of air cells is not
the design. The design is a surface of pieces. The pieces interact at their
edges (the material) and across the air between them (delayed
fluctuations); the air itself holds no state. The grid code was moved out
of the tree.

## Version 13 (2026-10-09): the solid surface (pieces, edges, air paths)

The airway is now the architect's surface. Each ring holds 8 wall pieces
lying on the section's ellipse. The ellipse sits against the outer wall,
which is the back of the pharynx and then the palate. Each section has its
anatomical width, and its depth is whatever the area then requires, so a
constriction becomes a slit along that wall. The tube rises from the glottis,
turns under the palate (outer radius 4 cm) and runs forward to the lips.

Connections between the pieces:

- **Axial air paths:** piece k of one ring links to piece k of the next.
  Each path's travel time follows its depth, so the inner side of the bend
  is shorter. The length is measured along the centres of the sections.
- **Cross air paths:** every pair of pieces in a ring links directly across
  the air. A path's travel time is its chord. Its share of the ring's air is
  the measure of the rays between the two pieces (s_i s_j cos_i cos_j).
- **Waveguide-mesh scaling:** the axial paths carry 0.6 of each ring's
  compliance. Each path's time is then scaled by the square root of its
  direction's share (×0.86 for the cross paths, calibrated below), so sound
  crosses the air at c in every direction.
- **Edges:** neighbouring pieces pull on each other with a membrane tension
  (`edge_tension`), as a pressure on the tissue port. The default is 0,
  because its effect has not been measured yet.
- **Walls:** each piece keeps its own mass-spring-damper port. The wetted
  area and the boundary-layer loss now come from the ellipse's own perimeter.

MEASURED (surface_probe, rigid, 192 kHz):

- **Plane waves:** with cross paths on, a uniform 3 cm² tube rings at
  460/1380/2320/3250/4200/5140/6080/7010 Hz. Without them it rings at
  470/1400/2350/3300/4260/5220/6180/7150 Hz, so dispersion is −1.5 to −2%
  (share 0.6). At share 0.5 it is −3.5% at 7 kHz. At share 0.3 the axial
  delays clamp and the result is wrong.
- **Cross modes:** a 3.57 cm circle, kicked off-centre, has its first
  cross mode at 4.95 kHz at timing 1.0 and 5.68 kHz at timing 0.86. The
  analytic value is 5.745 kHz. 12 pieces per ring give the same result as 8.
  Without cross paths there are no cross modes at all.
- **Independent strips are wrong:** without cross paths each strip of the
  bent airway is its own tube, and its spectrum fills with spurious peaks
  (REJECTED as a configuration).
- **Bent vowels:** /ɑ/ 670/1320/2400/3470/4280 Hz against
  710/1360/2440/3520/4350 Hz for the 1-D tube; /i/ F2 2240 Hz against 2320 Hz.
  Above 5 kHz the shape differs from the 1-D tube (split peaks near
  5.8–6.0 kHz).
- **Cost:** the full voice runs at a real-time factor of 7.5, against 0.58
  for the 1-D tube.

Global faults found by separating source and tract in the per-sample trace
(lip flow over glottal flow, sustained /ɑ/ and /i/):

- **The fixed 5.8 kHz notch.** It appeared in every vowel and was the hard
  ceiling near 5.5 kHz in every spectrogram.
  - Cause: a closed side tube ended at its last node, half a section short,
    so the 2 cm piriform branch behaved as 1.6 cm.
  - Fixed with a rigid end piece half a section beyond the last ring. The
    notch is now at 4.5–4.8 kHz, where the piriform fossae put it.
- **The source was too dull.** With medial roundness 1 mm² and contact
  softness 0.1 mm², the glottal flow derivative fell about 14 dB per octave
  above 1 kHz (−46 dB at 3–4 kHz, −57 dB at 5–6 kHz relative to 0–1 kHz).
  At 0.25 / 0.03 mm² it falls about 6 dB per octave (−36.5, −43). These are
  now the defaults.
- **Result:** voiced-frame band balance (relative dB) for the baritone on
  the three test lines, against the human's spoken half of three.wav:

  | Band | Human | Before | Now |
  |---|---|---|---|
  | 3–4 kHz | −21 | −25 | −24.5 |
  | 4–5 kHz | −37 | −26 | −28.6 |
  | 5–6 kHz | −37 | −40 | −33.8 |
  | 6–8 kHz | −41 | −39 | −33 |

  The hard ceiling is gone.
- **Transients:** the new airway renders show 4–10 transients; the 1-D tube
  with the same source shows 127 (pitch pulses).

OPEN:

- **A fixed bright line near 4.3–4.4 kHz** appears through every vowel in
  the spectrogram: the tract's F5, which barely moves with the vowel.
  Removing the epilarynx leaves it in place but 5 dB weaker; a 3 cm
  epilarynx moves the boosted resonance to 3.35 kHz (MEASURED, 1-D probe).
  The human recording shows a steady band near 3.5 kHz instead.
- **The cross modes the airway can carry are not excited.** The sections
  are symmetric left to right and the glottis drives every piece of the
  first ring alike, so below about 8 kHz only the bend breaks the plane wave.
  Still missing: asymmetric geometry (the epiglottis and vallecula, the
  sublingual space, the teeth as intercepting obstacles), oblique air paths
  between rings with ray-traced visibility, and edge tension.
- **Speed:** a real-time factor of 7.5.

## Understanding the mouth (2026-10-09): what "shore" asks of it

Verdict on version 13: "about the same; shore still sounds like shure". The
architect's direction: understand the mouth better, not just fix the sound.

**The formants.** LPC tracks (MEASURED, formant_track.py):

- Human "shore", 9.28–9.53 s:
  - F1 rises from 380 to 530 Hz.
  - F2 falls from 1550 to about 900 Hz (out of the "sh"), then rises to
    1270 Hz.
  - F3 falls from 2270 to 1600 Hz.
- Ours: F1 falls to 400 Hz and then 300 Hz. F2 sinks to 830 Hz and stays
  there for 150 ms. F3 drops to 1420 Hz.
  - A low F1 with a low F2 is a close rounded vowel, so we produce "shure".

**Tool.** `mouth_sensitivity.py` does three things:

- solves our own area function (dumped every 10 ms by `surface_say --areas`)
  as a plane-wave tube;
- reports, for each region, what narrowing there does to F1–F3 (Fant's
  sensitivity functions);
- finds the smoothest change of log area that takes the model's F1–F3 to a
  target.

The answer is the smoothest of many shapes that give the same three
formants, so it names regions, not exact anatomy.

**What the human /r/ of "shore" asks for**, against our /r/ (model 430 / 936
/ 1464 Hz → human 527 / 1244 / 1594 Hz):

| Region | Needed / ours |
|---|---|
| Lips | ×2.6 wider (1.07 → 2.7 cm²) |
| Lower pharynx | ×0.69 (narrowest ×0.28) |
| Teeth / alveolar | ×1.37 wider |
| Hard palate | the palatal constriction stays (0.2 cm²) |

In words: the human /r/ is barely rounded. It narrows the pharynx (tongue
root drawn back) and the palate (bunched tongue) and leaves a large space in
front of the constriction. Ours does the work with the lips.

That matches the published account of American /r/: constrictions at the
lips (slight), the palate and the pharynx, each near a velocity maximum of
F3, with a large front cavity and the sublingual space below the tongue.
That account comes from the literature; the inversion here agrees with it.

**The vowel of "shore".** Human 474 / 924 / 2240 Hz, ours 482 / 1152 / 2416
Hz (from the solved area function):

- soft palate narrower (×0.74)
- upper pharynx narrower (×0.84)
- teeth region wider (×1.23)
- lips unchanged

A back-raised [o] made with the tongue, not the lips.

**Sensitivities of our tract** (narrowing effect):

- **The lips lower F1 hardest** (−0.24 in the vowel) and lower F2 and F3
  as well. A lip-heavy /r/ therefore turns the vowel into "oo".
- **The pharynx** raises F1.
- **The hard palate** lowers F1 and F3 and raises F2 during the /r/.

OPEN:

- **The model's articulators can't produce the human /r/ as it is made.**
  /r/ is a tip narrowing plus lips; there is no tongue-root (pharyngeal)
  gesture and no sublingual space.
- **Next:** invert the whole human recording frame by frame into a map of
  what the human mouth does through the three lines. Compare it with what
  our planner and tongue do, and look for systematic gaps, not
  sound-by-sound fixes.

## Version 14 (2026-10-09): a floor under the tongue

GIVEN (architect): "the tongue is an articulating object in the mouth, not
the bottom. so we do need a floor."

**Geometry** (`MidsagittalTongue::floor_space`):

- The floor of the mouth lies under the tongue and is lowered by the jaw;
  the lower incisors close the space in front.
- The blade is a slab: its underside lies one thickness below its upper
  surface, 1.2 cm at the root and the pad's diameter at the tip.
- The frenulum holds the underside no higher than the line from the tip
  pad down to the floor at x ≤ 5.8 cm.
- Each section's space is weighted by how much of it the tongue covers, and
  opens forward into the mouth through its front face.

**Acoustics.** The space is a parallel surface tube (TubeGeometry::parallel)
beside oral sections 4.0–8.0 cm. It is closed at the incisors and meets the
mouth only through the face openings (`set_openings`). It is on by default
(`--no-floor` turns it off).

REJECTED on the way (all MEASURED):

- **A paper-thin blade** (underside drawn as a line from the pad to the
  frenulum): up to 4 cm² under the tongue, and the /sh/ peak fell from
  2.7 to 1.5 kHz.
- **Splitting the open column in front of the tip** into two channels
  joined across by 1.6 cm paths: the lower half rang as a separate pocket
  near 1.6–1.8 kHz.
- **The space appearing whole as the tip passed a section centre:** a click
  at the "j" of "jumped". Replaced by coverage weighting.
- **Pressure-wave paths whose admittance swings to zero.** A stored wave in
  a path that no longer acts on its ends gains the pressure difference of
  those ends on every round trip; the first opening under the tongue blew
  the voice up. The cross paths, openings and parallel-tube paths now carry
  power-normalized waves (√Y times the pressure wave). Calibration is
  unchanged: cross mode 5680 Hz, plane wave 460/1380/2320 Hz.
- **Linear per-block glides of those paths:** a slope change at every block
  edge, so tip onsets were 7–13× louder than without the floor. They now
  follow the two-pole-smoothed geometry every sample: ratio 0.9–1.5 and
  global peak ratio 1.00. Solid render: 4 transients, against 0 without the
  floor.

Results (MEASURED, solid baritone, three test lines):

- **Space under the tongue:** 0.5–1.5 cm³ during tip consonants, present
  31% of the time, at most 5.2 cm³.
- **"shore" /r/:** F3 1580–1600 Hz with the floor, against 1445–1470 Hz
  without it; the human value is 1600. F2 stays at 830–860 Hz against the
  human 1250, which is the lip rounding identified by the inversion. F1 is
  430–480 Hz against 500–530 Hz.
- **/sh/:** noise peak 2.3–2.4 kHz with the floor against 2.0–2.7 kHz
  without. The pocket under /sh/ is the known physical source of its low
  spectrum.
- **/s/ regressed:** "sells" 3.1 kHz (3.9 without the floor), "sea" 1.6 kHz
  (4.1 without).
  - Our /s/ is apical: the tip is raised to the ridge, leaving a pocket of
    about 1.25 cm³ under it.
  - English /s/ is mostly laminal: the tip is down behind the lower
    incisors, the blade grooves against the ridge, and there is almost no
    pocket. That missing pocket is what separates /s/ from /sh/.
- **Laminal /s/ (`tip_down` = 1) is REJECTED for now.** Our blade is a
  straight segment, so with the tip down it can't bulge up behind the
  teeth. It leaves a large front cavity instead, and /s/ falls to 2.5 kHz
  with or without the floor. The mechanism stays in the code, unused.

OPEN (what the floor exposed about the mouth):

- **The blade must be flexible:** able to bulge to the ridge while the tip
  stays down (laminal), or curl up and back (retroflex).
- **The jaw must move with consonants:** sibilants raise it until the teeth
  nearly meet, which also raises the floor.
- **/r/ needs its tongue-root gesture and light lips** (see "Understanding
  the mouth").
