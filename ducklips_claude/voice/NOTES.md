# voice — text to a feedback-driven speaking body (2026-10-11)

New engine, written from scratch. Imports nothing from the earlier engines.

## Pipeline (GIVEN by the architect)
text → `parse_text` (numbers, phrases, breath decisions, pitch accents)
→ `pronounce` (lexicon, then NRL-style letter-to-sound rules, then a stress heuristic)
→ `plan_speech` (Klatt-style durations, intonation, overlapping gestures)
→ body language (`score.hpp`: timed part postures plus markers; text form via `say --score`)
→ `Body` receives each event when its time arrives and deforms toward it.

## Body (`body.cpp`)
- **Air:** pressure waves cross one 0.365 cm section per step at 96 kHz and scatter at every junction.
  - Tubes: trachea (40 sections), mouth (48), nose (32), two piriform recesses.
  - Helmholtz cavities: three sinuses and the floor of the mouth under the tongue.
  - The bronchial end is an ideal pressure source for slow changes and mostly absorbing for sound.
- **Walls:** three pieces at every junction.
  - cover: elastic, edge-coupled
  - mass: edge-coupled
  - support: holds the posture and wanders
  - Each section is advanced with the trapezoidal rule and solved jointly with its junction pressure (precomputed 4×4, no iteration). It is passive, so a sealed mouth stays stable.
- **Folds:** lower and upper cover pieces joined at an edge, a body mass, and a cartilage support.
  - They oscillate only from lung pressure acting through the flow, and the flow is shaped by the waves returning from below and above. There is no driven waveform.
  - The flow solve combines Bernoulli and viscous slit losses with the incident waves, in closed form.
- **Turbulence:** the jet's dynamic pressure, bounded by the pressure drop, with Reynolds onset; scaled by jet size. Sources sit at the glottis and downstream of the narrowest place, or at the incisors.
- **Motors:** each part moves second-order toward its target. Holding wander is an OU process in the muscle command, scaled by effort.
- **Ear:** the speaker hears each glottal cycle and slowly corrects cricothyroid tension toward the intended pitch.

## MEASURED (Neo, Apple A18 Pro; tools in `tools/`)
- Real-time factor: 6.6–9× (single thread).
- Closed-glottis schwa: F 565/1524/2544/3593 Hz, bandwidths 43/73/106/125 Hz.
- Vowel postures fitted on the body itself (`vowel_fit`): most within 2% of Hillenbrand men; EH and AE about 5–7% low in F2.
- Rhotic posture: F3 = 1630 Hz.
- Sustained /ɑ/: open quotient 0.64, speed quotient 2.2, jitter 0.23%, shimmer 0.8%.
- Band tilt vs the owner's recorded vowel (dB re 0–1.5 kHz):

  | band | ours | owner |
  |---|---|---|
  | 1.5–3 kHz | −19 | −16 |
  | 3–4.5 kHz | −17 | −15 |
  | 4.5–6 kHz | −37 | −36 |

## REJECTED along the way (kept for the reasons)
- Explicit wall coupling: unstable with a sealed mouth.
- A 450 Hz wall-cover mode: absorbed F1.
- Turbulence ∝ Re²: exploded at releases.
- Area stepped at the 1 kHz control rate: a 1 kHz buzz.
- Pure-Bernoulli glottis: leaked the trachea through the closed folds.
- Lung source with 0.6 kPa source impedance.

## Second pass (MEASURED)
- **Jet resistance.** Constrictions now lose the jet's kinetic pressure (ρU|U|/2A²) in one closed-form junction solve that includes the wall and the noise source. Before this, /f/ passed 830 cm³/s through a 0.14 cm² gap; now it passes about 150.
- **Click hunting** (`tools/clicks.py`) found three causes:
  - Lip and teeth region edges snapped a whole section; they are now soft partial coverage.
  - The output had no anti-alias averaging of sub-steps.
  - The lower-cover damping was raised (0.30).
- **Voiced stops** enlarge the pharynx during closure. The "November" voicing gap fell from about 100 ms to about 20 ms.
- **Phrase-initial nasal:** the velum opens with the lips, which removes the nasal pop.
- **Soft voice onset:** the folds adduct before the pressure rises.
- **Intonation:**
  - Statements fall through their final vowel to about 85 Hz; questions rise.
  - Pitch is set before the phrase begins.
  - The ear's prior was recalibrated: about 100 Hz per unit tension, measured.
  - Contours now follow the owner's "Alfa"/"Bravo": level stressed syllable, final fall to about 85–92 Hz.
- **Frication** is band-passed around 1.2 × jet speed ÷ width (capped at 7 kHz). /s/ bands, dB relative to 4–12 kHz:

  | band | ours | reference |
  |---|---|---|
  | 0–1 kHz | −30 | −57 |
  | 1–2.5 kHz | −32 | −54 |
  | 2.5–4 kHz | −27 | −40 |
  | 4–6 kHz | −20 | −16 |
  | 6–9 kHz | −6 | −1 |

  /ʃ/ sits about 15 dB lower in frequency than /s/.
- **Listener roll-off:** first order at 9 kHz.

## OPEN
- /s/ low-band rejection is still about 20 dB short of the reference; it is too strong above 9 kHz.
- Nasal murmur strength.
- EH/AE F2.
- The one-constriction tongue body.
- Prosody has been compared only on isolated words.
- Not yet judged by a listener.

Build: `cmake -S . -B build && cmake --build build`, then `build/say "text" out.wav [--score f]`.

## REJECTED (third pass)
- **A second front-dorsum shaping point for EH/AE.** F2 rose only to 1730 against 1840 Hz (EH) and stayed at 1590 against 1720 Hz (AE). The limit is the tongue model, not the posture.
- **A longer tongue-tip channel for /s/.** Low-band rejection did not improve (0–1 kHz went from −30 to −20 dB, and 1–2.5 kHz stayed near −33). The low-band /s/ energy is not back-cavity leakage through the tip.

## Third pass: what the spectrograms actually showed (2026-10-11)
The owner judged every sample garbage. Looking side by side, and not at band totals, showed five structural faults:
- **The source was a flat comb.** Every harmonic to about 4 kHz was equally bright, about 30 dB too much above 3.6 kHz. The cause was instantaneous closure along the whole fold length.
  - Fix: an 8-slice zipper closure (a sine vibration shape along the length, unit mean).
- **Each period never went quiet.** A flow ramp and aspiration fizz ran through the whole cycle.
  - Fix: aspiration now scales with how far the folds rest apart, so it is about 6% during modal voicing.
- **Vowels were held blocks.** Fix: the tongue now arrives late in each vowel, and vowels drift 30% toward schwa as they end.
- **The ear misheard octave errors** and drove tension up until voicing died. Fix: it now ignores percepts more than ±30% from the intention.
- **The folds blew open** under full pressure at low pitch, and on rounded or high vowels.
  - Fix: stretch scaling (mass / √T, stiffness × T^1.5, thickness / √T).
  - Fix: lung pressure now follows pitch (+5 Pa per Hz around 750 Pa at 125 Hz, clamped 560–1000). This matches the physiological correlation between subglottal pressure and F0.
- **Other changes:**
  - Epilarynx 3.0 cm × 0.45 cm² (F3/F4 cluster).
  - Wider /l/ side channel.
  - Rate 0.85.
  - Stronger phrase-final lengthening (×2 on the last vowel).
  - Pitch floor 90 Hz.

**MEASURED** on the Alfa vowel against the owner's:
- 4.5–6 kHz: −38 dB, owner −35.
- Stable voicing 100–155 Hz for /i ɑ æ u oʊ/ at speech pressures, jitter about 0.5%.

**OPEN:**
- Open quotient is about 0.7, with no pressed or modal range: adduction stops the oscillation instead.
- F2/F3 voiced bandwidths are wider than the owner's.

## Fourth pass: firm source by search (MEASURED)
`tools/fold_search.py` ran a random search over the 13 fold materials, scored on open quotient near 0.5, jitter and pitch accuracy for /ɑ/ and /oʊ/ at 105 and 145 Hz.
- Result: open quotient 0.46–0.54, jitter 0.2–0.3%, the old one-regime slackness (about 0.7) gone.
- The /ɑ/ long-term spectrum now tracks the owner's from 0 to 8 kHz, except the owner's F3/F4 cluster at 2.7–3.4 kHz, where ours is about 15 dB weaker.
- CANDIDATE: the fitted lower "thickness" (1.03 cm) is a pressure-to-stiffness ratio, about 3× anatomical; it is not a measured dimension.
- Rate 0.95.
- The lower pressure depth is reinterpreted as medial surface plus conus elasticus projection (HYPOTHESIS, plausible at about 1 cm).
- REJECTED: matching the owner's 2.7–3.4 kHz cluster with the epilarynx alone. Settings from 3.0 to 4.0 cm long and 0.25 to 0.45 cm² never matched the 1.5–3 and 3–4.5 kHz bands together. The cluster likely needs pharynx and piriform geometry as well.
- REJECTED: the piriform recesses as the source of the owner's 3–4.5 kHz energy. Depths of 1.5–2.9 cm and openings of 0.9–1.6 cm² only lowered that band (−24 to −32 dB against the owner's −15).
  - NEXT HYPOTHESIS: source energy at 3–4 kHz (closure shape), or a narrower F4 bandwidth.
- MEASURED: a sustained /æ/ at 125 Hz and 750 Pa gives −18 dB at 3–4.5 kHz, against the owner's −15. The 9 dB shortfall seen in the word "Alfa" is therefore prosodic (pitch and pressure path), not the source or the geometry.
  - Lower cover damping overshoots (−10 dB); the current setting is kept.
- REJECTED: one material set covering breathy (gap +0.03), modal and pressed (gap −0.02) voice. In a 90-trial search, pressed voice oscillated (open quotient about 0.5) only with a pressure depth of 3–5 cm, which is unphysical, and breathy never oscillated.
  - HYPOTHESIS: the effective pressure depth (the medial and conus surface exposed to subglottal pressure) must change with adduction. Pressed folds expose less of their lower surface; breathy folds expose more.
  - The engine keeps the modal set.
- MEASURED (by audio; `tools/glottal.py` open-quotient and cycle counts proved unreliable at these settings):
  - gap +0.02: no oscillation (periodicity 0.05).
  - gap 0: modal (periodicity 0.94, H1–H2 −3 dB).
  - gap −0.02 to −0.04: steady, but upper harmonics fall 10–28 dB. That is duller, the opposite of real pressed voice.
  - OPEN: adduction does not map onto the breathy–modal–pressed continuum. The planner uses only gap 0, so speech is unaffected.
- Breathy continuum restored: cover damping now falls as the folds rest apart (×0.15 at gap 0.02), so compressed tissue is lossy and parted tissue carries the mucosal wave.
  - MEASURED: gap +0.01 voiced with steeper tilt; +0.02 voiced and noisy (periodicity 0.67); +0.03 whisper; gap 0 modal unchanged.
  - Pressed voice still dulls instead of brightening (OPEN).
- Pressed voice fixed in direction: the zipper vibration shape flattens with adduction (flat at gap −0.03), so pressed folds close along most of their length at once.
  - MEASURED at 120 Hz, 750 Pa: gap −0.015 gives H6 +16 dB re H1; gap −0.03 gives H5 +11, H6 +23, louder. Modal unchanged.
  - Adduction now spans breathy → modal → pressed in the correct spectral directions.
  - At breathy gap +0.01 the harmonics fall but the 1.5–4.5 kHz band totals rise. Probably hiss; not separated (OPEN).
- The planner uses voice quality:
  - Stressed vowels get slightly pressed folds (gap −0.008).
  - Phrase-final vowels relax to breathy (gap +0.012) through their second half.
  - Vowels after voiceless consonants keep their aspiration timing.
- REJECTED: closing the Alfa 3–4.5 kHz gap (now 7 dB) with a stronger stressed press. Gaps of −0.015 and −0.022 gave −25 and −21 dB against −22, with no trend.

## Fifth pass: "sounds like it hits the microphone limit" (owner)
- **Aliasing.** Two-sample averaging left energy rising again at 16–20 kHz (−38 dB, against the owner's −61). Replaced by a 4th-order 10 kHz Butterworth at 96 kHz before decimation; now −57 to −65 dB.
- **Gated envelopes.** Vowels sat on a flat plateau and words stopped in about 20 ms, like a limiter.
  - Each syllable is now a breath pulse: pressure peaks early in the vowel and decays to 72%.
  - Pressure eases (×0.65) during voiceless consonants.
  - Phrases run down: pressure decays from mid final vowel while voicing continues, and the folds part 0.1 s after the end.
- **Exact kinetic jet loss.** The quadratic is solved per junction, which removed the lagged-resistance ringing at near-closed openings.
- **Airway resistance** at the lung end: 2 dyn s/cm⁵.
- **Frication** gain ×8. MEASURED against references, dB re vowel: /s/ −17.4 (ref −16.9), /f/ −17.3 (ref −13.2).
  - Fixed: the lip-jet drop is taken to the room.
  - Labiodental obstacle gain 10.
- Vowel breath decay deepened to 55%. Loudness of voiced frames over Alfa/Bravo/Charlie: mean −11.5 dB (owner −13.6), lower quartile −16.6 (owner −20.1).

## Sixth pass: four time-frequency views and a human reference corpus (2026-10-12)
**Owner's hearing of the samples:** "elpha bravo farley delfa efo fa-fa/oooo. my mame is papa and i speap with a baeelee...". Stops and affricates are heard as /f/, /n/ as /m/, and some vowels shift.

**Transforms** (`analysis/time_frequency.cpp`, `build/tf_analyze`): all four run in about 0.5 s on 3.2 s of 16 kHz audio.
- Garcés constant-Q Gabor atoms: Q_N = 1/(2^(1/2N) − 2^(−1/2N)) and M_N = 2√(ln 2)·Q_N, computed per band by a small baseband inverse FFT.
- Lindeberg time-causal Gabor: limit kernel with 8 recursive stages, c = √2, mean-delay compensated.
- Smoothed pseudo-Wigner–Ville, enhanced by the minimum over three window shapes (cross-term consensus).
- STFT.
- Each is enhanced like the owner's STFT: fusion across three resolutions plus synchrosqueezed texture.

**Reference corpus:** CMU ARCTIC bdl, utterances a0001–a0060 with phone labels (free licence, `ref/arctic_bdl/COPYING`).

**Comparison tool** (`tools/compare_arctic.py`): renders each transcript and aligns phones, then reports level and level-normalized shape per phone and class in all four views. It also has a one-line summary for tuning.

**MEASURED** on a0001–a0020 (shape = RMS dB difference, cq view):
- Vowels and nasals about 8 dB; stops and fricatives about 12–16 dB.
- /t/, /z/, /k/, /p/, /dh/, /sh/ worst.
- /hh/ was 17.6 dB too quiet (now −11); /s/ is 12 dB quiet.

**Changes:**
- Fast tip and lip gestures (28 ms) that aim past their target.
- Aspiration rises steeply with fold parting.
- Voiced fricatives part the folds slightly (0.012).
- American flapping of /t, d/ before unstressed vowels.
- Vowel breath decay back to 72%.
- Frication 1.5 and aspiration 2 (best of a 3×2 scan; the total moved only 1%).

**OPEN** (from the a0003 four-view comparison):
- Speaking rate is 40% slower than bdl.
- Pitch is flatter.
- Formant sweeps are weaker.
- /s/ lacks the compact 4–7.5 kHz front-cavity peak.
- Click-like vertical lines remain at some consonants.

## Seventh pass: tuned against ARCTIC with objective scores
The summary score is fricative |level| + stop-release |level| + consonant shape + vowel shape, in dB, over utterances a0001–a0020. It went from 35.3 to 31.7.

- **Vowel durations:** Klatt's inherent durations ran vowels 2.2× too long for read speech. Scaled by 0.5 at rate 1.0. Duration ratio ours/human: utterance 0.98, vowels 1.05, consonants 0.87.
- **Pitch:** range and movement already matched (4.5 vs 5.0 semitones; 17.9 vs 17.2 semitones/s). Base pitch raised from 112 to 120 Hz (human median 122).
- **Vowel postures:** IY IH EY EH AE AA refitted to the bdl medians (`tools/vowel_survey.py`); AX and UW given connected-speech targets. LPC is unreliable for back vowels there, so those keep their fits.
- **Formant movement** is still about 35% slower than the human's (F1 ~1000 vs 1500 Hz/s; F2 ~1970 vs 2820 Hz/s) (OPEN).
- **Jet turbulence:**
  - Pressure now grows as v³ (dipole), so the glottis behind an oral constriction no longer out-hisses the jet.
  - The source is distributed downstream of the constriction (1.2 cm decay, ×4 at the incisors), with no if/else source placement.
  - Fricative levels vs human: /s/ −2.4, /z/ 0, /ʃ/ +1.9, /f/ +1.3 dB.
  - /s/ constriction at 15.6 cm.
- **Aspiration growth** is 1 + 12·parted². /h/ went from −17.6 to −3.7 dB.
- REJECTED: leading stop and nasal closures by 1.0–1.4 movement times (total 31.8–32.0, against 31.7). The remaining stop "steady" shape error (about 15 dB) lies inside closures, where the human recording carries faint voicing and room sound. It may matter little perceptually.
- Vowel drift toward schwa removed (total 31.66 → 31.43).
- REJECTED: a syllable-rate jaw cycle realised as incisor-gap narrowing (gap 1.5 or 0.8 cm²) on oral consonants. F1 movement was unchanged (about 910 Hz/s against the human's 1508).
  - OPEN: a real jaw lowers F1 by narrowing the whole front oral cavity. It needs to enter the area function, and the vowel postures must then be refitted with the jaw as a parameter.

## Eighth pass: closures never sealed (found from the owner's transcription)
**Owner heard:** "delpa bravo sharlie delsa echo foxtraha / broccoli ??? ??? the azy og ?? easily".

**MEASURED** with the new narrowest-constriction trace (`CONSTRICTION=1 build/say ...`):
- /d/ in "dog" only narrowed to 0.065–0.098 cm². Pressing contact was a parabola reaching zero within ±0.08 cm, narrower than one 0.36 cm section, so no tongue closure ever sealed. Stops became fricatives: t → s, d → nothing.
- Fix: the pressing root goes from −0.15 to −0.5, spreading contact over about 0.5 cm. Every stop, nasal and affricate in a 20-closure test now seals to 0.000.
- /l/ side channel set back to 0.20 cm².
- The ARCTIC spectral score barely moves (31.7 → 31.6): it does not measure place cues, and silent closures differ in shape from the human recording's room-tone closures. Release shapes improved (/t/ 11.9 → 10.4).
- Formant-movement deficit and jaw: see the seventh pass (OPEN).

**Routing:** the Mini is reached over Wi-Fi (`ssh m4mini`); AWDL works only nearby.
- The jaw is now in the area function: closure narrows the mouth by up to 45% from 8 cm to the lips. Oral consonants raise it and vowels drop it.
  - Score 31.6 → 31.8 (neutral); F1 movement unchanged (883 Hz/s).
  - HYPOTHESIS: the human F1 "movement" of 1500 Hz/s (15 Hz per 10 ms frame) is partly LPC tracking jitter on a real recording. That statistic is no longer used as a target.

## Ninth pass: the owner's verdicts against the phone-test tuning
**Owner (pass/fail on the tuned-inventory samples):** echo fail; foxtrot fail but less; other NATO words better. In the fox sentence only "lazy" was legible. "Is it raining today?" was barely legible. "Hello my name is..." then nothing. "Sea shore" came out as "si-shooor".

**MEASURED** (ARCTIC a0001–a0060, `tools/render_arctic.py`, per-phone mean duration ours/human):
- /iy/ 0.77, /sh/ 0.52, /ao/ 1.43, /aa/ 1.56, /ch/ 0.45, /jh/ 0.64. This is exactly "si-shooor".
- The kNN tuner had moved durations and back vowels where its fingerprint, not a listener, liked them.
- Inherent durations were rescaled by the ratio. Every phone is now within 4% of the human.

**Back vowels and AH:** taken from the engine impulse codebook (`out/inversion/engine`), matched to Hillenbrand 1995 men. Jaw open, lips constrained by rounding class. UW aims at the fronted American /u/ (F2 about 1450).
- Running-speech OW went from 416/1497 to 537/952 Hz.
- AH went from 419/1709 to 603/1530 Hz.

**Rhotics** (`vowel_fit rhotic`, wider grid):
- The tuned ER measured 400/1280/2880 Hz in isolation: no r-colour.
- The fit gives 480/1370/1720 Hz: dorsum bunched at 13.0 cm (0.40 cm²), root 5.0 cm (2.0 cm²), lips 1.2 cm², protrusion 0.3, floor 1.0.
- Applied to ER and R. Running-speech ER F3 went from 2960 to 2378 Hz (undershoot remains).

**Phone test:** 30.0 → 29.0% (flat). REJECTED as the sole guide: it scored the changes behind "si-shooor" as improvements.

**Later in the ninth pass: three engine faults found from the "Echo"/"Foxtrot" spectrograms against `ref/nato16`.**

1. **The vocal folds glided shut instead of snapping shut.** MEASURED with `vowel_probe` TRACE and `tools/source_slope.py`:
   - At rest gap 0 the folds met the midline at the end of their swing. Flow decayed over about 1 ms (MFDR/peak flow 1.1 per ms), and the opening was the sharp corner.
   - That is the reverse of a human pulse. The flow derivative fell about 10 dB/octave, and vowels in running speech lacked about 8 dB at 2–4 kHz against bdl.
   - Fix: modal voicing at rest gap −0.015 (`kVoicingGap`), stressed −0.02. This gives gradual opening, abrupt closure, a closed phase, and MFDR/peak 1.4–1.8 per ms (human range).
   - The glottal viscous loss now uses a 0.3 cm channel depth (`viscous_depth`) instead of the 1.06 cm fitted thickness (HYPOTHESIS: only the narrowest channel's depth matters).
   - Result: the 2–4 kHz band went from −22 to −18 dB (human −15). LPC now finds F3 where the human has it (EH 2832/2851 Hz, AE 2827/2786, IY 3129/3073).
2. **Lip openings were too wide.** The bdl-refit postures used 8.6–14 cm².
   - MEASURED (`tools/transfer_probe`): at 9 cm², F3 and F4 bandwidths are about 250 Hz. At 4 cm² they are about 130 Hz.
   - All vowels were refit (`vowel_fit bounded`) to Hillenbrand men within anthropometric lip ranges. The glottis was 0.01 apart during the fit: fully closed, the subglottal tract couples through the chink and splits F1 into a peak plus a broad bump the fitter mistook for F1.
   - Rhotic refit: dorsum 13.0/0.25, root 6.5/2.0, lips 1.2, protrusion 0.7 gives 470/1380/1750 Hz.
   - OPEN: front F2 tops out near 1700 Hz for EH/AE in this geometry.
3. **Tip and dorsum constrictions aliased on the section grid.** The tip parabola (4·d²) was sampled at section centres 0.36 cm apart, so a tip aimed between them opened by up to 0.13 cm².
   - "Foxtrot" /s/ → /t/ reopened mid-slide (0.081 → 0.010 → 0.117 → 0 cm²). The jet burst into an 8 kHz click, and the loudness normaliser then turned whole utterances down (NATO sentence gain 0.29).
   - Fix: each section takes the narrowest point within its span. Fricative and lateral areas were raised so their steady openings stay unchanged (S 0.133, Z 0.122, SH 0.120, ZH 0.178, TH 0.45, DH 0.276, L 0.248).
   - NATO gain is now 1.05. The remaining crest peaks (about 10× local RMS) are /t/ releases into /l/ and /ð/.

Phone test: 28.8%. New tools: `render_arctic.py`, `hvd_test.py` (LPC misses F1 on some synthetic vowels; unreliable), `transfer_probe`, `area_formants.py`, `source_slope.py`; AREAS snapshots in `say` and `vowel_probe`.

**Ninth pass, fricatives** (band shares in `/tmp`-style surveys; per-phone dB share by band, ARCTIC a0001–a0060):
- **Jet noise source.** It was a narrow band-pass at 1.2·v/d, clamped at 7 kHz, which put every sibilant and burst at 6–8 kHz (/s/ 6–8 kHz −1.3 dB share; human peaks at 4–6 kHz).
  - Now broadband: the Strouhal corner is 0.2·v/h on the slit height h = area / 1.2 cm, clamped to 1.5–8 kHz. The low skirt rises 12 dB/octave up to a third of the corner (a dipole), and the response is flat to the corner with −6 dB/octave above.
  - Unit variance at every corner, so the v³ law sets the power.
- **Turbulence grows downstream** (1 − e^(−d/0.4 cm)) instead of sounding at the constriction exit. A source there drove the back cavity and filled /s/ with 1–3 kHz.
- **/s/ and /z/ place 15.45.** At 15.15 the shape was best but the level was 10 dB too loud (front-cavity alignment). At 15.6 the level was 6 dB quiet and the peak too high.
  - MEASURED /s/ shares, 0–8 kHz: −25 −34 −26 −25 −19 −11 −6 −2 (human −32 −34 −29 −28 −15 −6 −5 −7). Level −4.9 dB (human −6.0).
  - /ʃ/ now peaks at 2–4 kHz like the human.
- **Bug: /f/ was intermittent.** The noise loop started at the narrowest section. When the minimum flickered onto the last lip section (47, no downstream junction), nothing was injected, so /f/ came out as sparse crackle.
  - The loop now always covers the lip junction. /f/ is continuous broadband noise (shares about −15 dB flat across 1–8 kHz, human −13 to −15).
  - Phone-test /f/ went from 33% to 65%.
- **Voice onset time** (`tools/vot.py`, from the body trace): about 60–75 ms before stress and about 30 ms unstressed, human-like. Aspiration is still weak against the release puff (OPEN).
- Phone test 29.9%.

**Ninth pass, airflow** (found from Kilo, Papa and Tango against `ref/nato16`):
- **Bug: kinetic loss charged at every narrow junction.** The junction solve applied ρU|U|/(2A²) on the narrower side at every junction, so a constriction N sections long cost N times the pressure.
  - MEASURED: in "Kilo" the /i/ channel after the /k/ release ate about 6000 dyn/cm². Post-release flow was 200 cm³/s with only 50–100 dyn/cm² across the open glottis, so there was almost no aspiration (about 80 ms of near-silence).
  - Fix: Borda–Carnot. The loss ρ/2·(1/A_up − 1/A_down)² is charged only where flow expands; contractions and uniform channels are lossless. Post-release flow is now about 690 cm³/s, with a strong aspiration burst like the human.
- **Lip exit loss.** With Borda–Carnot inside the tube, air leaving through the lips had no loss at all, and /v/ voiced as loudly as a vowel. `radiate` now solves the exit loss ρU|U|/(2A²) jointly with the radiation load (inflow lossless).
- **Fold flutter at abduction.** With realistic flow, the 85% cover-damping cut for folds resting apart let them flutter without bound at phrase-end abduction (glottal area 0.01 → 0.47 cm² → NaN).
  - The cut now applies only in the breathy range: it fades out between rest gaps 0.04 and 0.08.
- **Room.** At −16 dB its tail buried weak consonants: /θ/ read −21 dB against the speaker's −31.5; dry it reads −28.7. The reverb level is now −32 dB (the references are close-miked).
- **Calibration** (dB re utterance, ARCTIC a0001–a0060):
  - Aspiration 2.0 → 0.5; frication 1.5 → 0.75; /h/ rest gap 0.035 → 0.09.
  - Phone settings: S/Z place 15.6; Z and V area 0.05; F area 0.18; TH area 0.04 at 16.7; DH area 0.01 at 16.7.
  - Ours vs human: s −4.6/−6.0, z −8.0/−8.1, f −27.1/−27.7, θ −29.3/−31.5, v −10.1/−10.1.
  - OPEN: ʃ −10.0/−4.7 (shape right, 5 dB quiet), ð −12.7/−18.9, h −3.6/−11.2 (too much voicing).
  - Stop-release windows are not comparable: the label conventions differ.
- **Phone test 30.9%** (f 77%, ð 38%, z 44%, b 62%, iy 55%). θ, jh, uw, ng, uh remain at 0%.

**Ninth pass, liquids** (fine-band envelopes, dB re the strongest band, ARCTIC a0001–a0060):
- **/l/ was heard as /r/** (top confusion, 33 cases).
  - The human /l/ has a deep F2–F3 valley (1.25–2.25 kHz at −35 to −45) and a strong 2.5–3.5 kHz region. Ours had a shallow valley and a 2250 Hz peak.
  - Cause: the planner gave /l/ a tongue-body posture only in codas (hard-coded), and its table posture was the neutral one with lips at 1.6 cm².
  - `tools/lateral_fit` (chain matrix around a fixed tip contact, jaw raised) gives a dark /l/ at tongue 8.25/0.25, lips 1.0, tip 0.18 at 16.7 cm: 380/920/2730 Hz against bdl's 380/950/2850.
  - Now `kDarkL`, used in every position. Ours: valley −35 to −44, 2.5 kHz −29.
- **/r/** now takes its table lips (1.0, protrusion 0.7) and floor instead of hard-coded 1.3/0.4/2.5. Floor 2.5 (sublingual cavity) matches the human /r/ envelope best (2 kHz −33 vs −37).
- **Phone test 32.1–32.6%**, the best so far (l 19%, iy 61%).
- Checked and left alone: the peak-to-median-vowel ratio is 8.0 (human 5.9). The peaks are stressed nuclear vowels with pressed voice, not artefacts.

**Ninth pass, affricates and onsets:**
- **/dʒ/ had no frication** ("Juliett" went from its burst straight into the vowel). The affricate release opened to 0.165 cm², more than twice the /ʃ/ fricative aim.
  - Now: /tʃ/ opens to 0.072 (the /ʃ/ aim), and its mid-frication spectrum matches /ʃ/ (2–4 kHz).
  - /dʒ/ opens to 0.04 with the folds parted (0.012) like the voiced fricatives. 0.035 hissed /s/-like; 0.05 was masked by voicing.
- **The voicing-onset "click"** is the F1 ringing of open and back vowels. Lip flow swings ±200 cm³/s while glottal flow is 30–200, and subglottal pressure peaks near 1000 Pa at closure.
  - The peak-to-median-vowel ratio is 9.3 (human 5.9). Peaks fall mostly in /ɔ aɪ oʊ ɑ/ (e.g. "for").
  - REJECTED (no gain): two-stage adduction at voicing onset (9.5); lower speech pressure 600 Pa (9.4); stressed vowels unpressed (8.6–8.8, phone test −0.2).
  - OPEN: F1 bandwidth during voicing may now be narrower than human. The removed per-junction loss had been damping it under DC flow.
- Phone test 31.7–32.6% across these runs (±0.5 is run-to-run noise from the changed inventory).
- **Vowel envelopes** (fine bands, ARCTIC): /iː ɪ/ match the human closely. All others are 5–12 dB weak at 2.75–3.5 kHz (the F3/F4 cluster). The /ʌ/ "F2 2100 Hz" from LPC was an F1/F2 merge: in situ the area function gives 710/1270 Hz.
- **Fold source scan** (`vowel_probe` EH at −0.015, output 2–3 / 3–4 kHz re strongest harmonic; baseline −33/−32):
  - Lower thickness 0.6: −24/−27. Convergence +0.01: −27/−29. Shear ×3: −25/−27. Shear ×0.3 or upper stiffness ×10: worse.
  - REJECTED in running speech, lower thickness 0.6 / 0.5: the vowel 2–4 kHz balance approached the human (−16/−20 at 0.6), but the crest exploded (19.6 / 27.8) and the phone test fell to 26.4 / 25.2%.
  - OPEN: a thinner pressure-bearing depth brightens the source but destabilises running phonation. Find what bursts before revisiting.
- **Bursts at vowel → /t/ closures.** Every voiceless consonant abducted the glottis 30 ms before its start. For stops this let /h/-like flow (650 cm³/s) run through the open tract, and the tongue tip then closed into that stream with a turbulence burst (the largest peaks with thinner folds, and "ah>t"/"ae>t" peaks at normal thickness).
  - Stops now abduct at closure; fricatives and /h/ still lead by 30 ms. Crest 9.3 → 8.8; voice onset times unchanged; phone test 31.7%.
  - Thinner folds retried after this fix and still REJECTED (0.8: crest 9.4, test 30.0%; 0.6: 15.3, 26.3%).
- **Phrase-final stops.**
  - Moving the stop's lung command to closure put a final /t/'s pressure command after the phrase-end lung release, so the lungs re-pressurised. The lung command is back at start − 30 ms; only the glottis waits for closure.
  - The planner also re-adducted the folds for voicing after a phrase-final stop even with nothing following (flapping noise at the end of "Foxtrot"). Voicing now resumes only if a voiced phone follows in the phrase.
  - Crest 8.4 (human 5.9); phone test 31.7%.
- **Epilarynx lengthened** (top 3.0 → 3.5 cm; narrow tube about 2.5–3 cm). At 3.0 the F3/F4 cluster sat near 4 kHz, and all vowels but /iː ɪ/ were 5–12 dB weak at 2.75–3.5 kHz.
  - MEASURED vowel band balance (re 0–1 kHz) at 3.5: −14.4 −15.9 −18.0 −33.5 −40.9 −43.5 (human −14.6 −15.0 −15.2 −32.4 −37.3 −39.4).
  - Crest 8.2; phone test 32.0%. A top of 4.0 overshoots (phone test 29.2%).
  - REJECTED: refitting the vowel, /ɝ/ and dark-/l/ postures on the new geometry. The impulse fits improved (AE F2 1990 Hz), but running speech scored worse (phone test 30.0%, crest 9.7, /l/ 7%). The earlier postures are kept.
- **Aspiration source.** The glottal noise now uses the same broadband Strouhal model as the jets in the mouth: corner 0.2·v/h on the glottal slit height (area / fold length), 12 dB/octave below a third of the corner, −6 dB/octave above. The near-white source had made /h/ a 6–8 kHz hiss.
  - Recalibrated: aspiration 0.5 → 0.1; /h/ rest gap 0.09.
  - /h/ shares 0–8 kHz: −5 −22 −14 −6 −8 −19 −18 −14 (human −1 −15 −21 −17 −17 −22 −26 −26); level −5.9 (human −11.2).
  - Phone test 32.4% (/h/ 31%). Gap 0.07 matched the spectrum best but voiced too loudly (level −2.2; /h/ heard as a vowel, 4%).
- **Formant movement** (`tools/formant_stats.py`): |dF2/dt| is 2251 Hz/s (human 2820); |dF1/dt| 1116 (human 1508).
  - REJECTED: vowel tongue arrival at 35% or 25% of the vowel (was 50%). F2 velocity rose only to about 2300, phone test 31.8–31.9%; the articulator motors' own dynamics set the pace.
- **"Echo" against the owner's own recording** (`ref/user_nato.wav`, about 6.2–7.0 s) and `nato16/nato_e`:
  - Vowel durations are close (E about 120–130 ms; owner's O 330 ms, NATO's 230, ours 240).
  - OPEN: our /k/ release is a single transient, then about 20 ms of weak noise and a strong first voicing pulse. Both references show a dense burst-frication of 20–30 ms.
  - OPEN: a steady band near 5.1 kHz runs through our vowels that the owner's recording lacks (the ARCTIC average at 5–6 kHz is not high, so it may be speaker-specific).
- **/k/ aspiration** (the owner failed "Echo"). Hiss above 1 kHz in the first 45 ms after release, re the following vowel: NATO "Kilo" −10.5 dB; ours −13.2 dB at aspiration 0.1, of which the glottal source is most (frication alone −20.6).
  - Aspiration 0.13 matches. /h/ level −4.2 (human −11.2); phone test 31.7% (32.4% at 0.1, within noise; /h/ 22% vs 31%).
  - REJECTED: a slower dorsal release (35–50 ms instead of 22) to lengthen the velar burst. No visible gain in the spectrograms; the test's /k/ fell from 17% to 8–10% and /g/ from 7% to 3%.
- REJECTED, longer jet-turbulence decay (2.0 and 3.0 cm instead of 1.2) to lift /ʃ/ (−10 dB, human −4.7). /ʃ/ did not move (−10.3 / −10.0), /s/ rose to −3.5, phone test 31.3–31.5%. The /ʃ/ deficit lies elsewhere (OPEN).

## Tenth pass: the owner's transcription of the ninth-pass samples
**Owner heard:** "hello my name is bekhar and i speak for everybody" / "elfa, bra?o *echo *foxtrot" / "uplep prawn fox fun cover the azy bag" / "is it raining today has a lot of pops".

- **Pops** (`tools/pops.py`: 2 ms bursts more than 2.2× the 30 ms level): 1.61/s against bdl's 1.18.
  - Cause: after a voiceless consonant the folds were brought together in 40 ms. The first pulses swung twice as wide as a vowel's (glottal area 0.19 cm², flow 425 cm³/s) and slammed shut.
  - Adduction now takes 70 ms (`kAdductionMove`), giving 1.36/s. Voice onset times are 84/41/68 ms for Kilo/Papa/Tango.
  - Trade-off: short vowels after voiceless consonants are breathier through their middle (crest median 11.3). At 40 ms: 1.78/s and crest 8.6.
- **No voice bars.** Human voiced stops (NATO Bravo/Delta) carry 20–30 ms of voicing in the closure. Ours were silent: a sealed tract filled to subglottal pressure.
  - HYPOTHESIS confirmed: the Kelly–Lochbaum tube had no volume term, so actively enlarging a section drew no air.
  - Added `area_rate` (dA/dt, smoothed by two 50 Hz stages so the 1 kHz control steps do not buzz) as a volume sink L·dA/dt in `wall_sink`.
  - The sink is restricted to the wide pharynx (junctions 8–24, both sections > 0.3 cm²). Everywhere else the parametric tongue does not conserve its own volume (tongue moves pushed out 500 cm³/s), and at narrow releases the lagged rate forced flow into closed sections (ringing bursts).
  - The voiced-stop enlargement now runs through the last 80 ms before release, and the folds slacken to touching (gap 0.0) meanwhile. Held pressed, the enlargement sucked the cavity below atmospheric and /g/ releases imploded (lip inflow −1400 cm³/s; peak 29.7× the vowel level).
  - The enlargement relaxes over 150 ms after the release. "Delta" now shows a 30 ms voice bar.
- **Voiceless stop before a voiced obstruent:** no aspiration lag ("quick brown": the /b/ closure had gone unvoiced, heard /p/).
- **Unreleased stops:** a stop before a stop or affricate at another place opens 20 ms after the next phone starts, inside its closure. "quick brown" had burst into silence.
- MEASURED: pops 1.36/s (human 1.18), worst spike 14.8× vowel rms, vowel band balance unchanged, phone test 32.1%.
- OPEN: lazy → "azy" (onset /l/ weak); dog → "bag" (place cue of /d/); jumps → "fun"; Papa → "bekhar".
- **lazy → "azy":** /l/ now releases in 18 ms (the tip snaps; F1 jumps): /l/ 9–13 → 14–16%.
  - REJECTED: a clear onset /l/ halfway to the next vowel. F2 did not move (the lips set it) and the envelope moved away from bdl's.
- **dog → "bag":** the /d/ burst level matched (3–7 kHz re vowel: ours −5.3, NATO −5.2 dB); the place cue was the transition. F2 at release into /ɑ/ was 1570 Hz for "dog" and 1430 for "bog" (human about 1700 vs 1000).
  - Alveolar stops and nasals now draw the tongue body to a fronted locus posture (`kAlveolarLocus`, /ɪ/-like) during the closure, and return after the release.
  - Labial stops and nasals open their lips over 40 ms (was 22).
  - The voiced-stop enlargement widens whichever closure posture applies and relaxes over 100 ms (was 150).
  - "dog" F2 now falls 1945 → 1755 Hz; "bog" rises 1360 → 1525.
  - Phone test 33.7% (best; d 29%, n 60%); pops 1.29/s.
- **Papa → "bekhar":** a lexicon error. "papa" was P AH0 P AA1 (stress on the second syllable), so the first /p/ came out unstressed and unaspirated (heard /b/).
  - Now P AA1 P AH0; the first /p/ has VOT 70 ms. Also fixed: bravo B R AA1 V OW0 (was secondary/primary reversed); "x" EH1 K S (X-ray had no vowel).
- **jumps → "fun":** /dʒ/ had a silent closure and almost no frication.
  - Voiced affricates now get the voiced-stop pharynx enlargement.
  - The frication phase starts at 40% of the affricate (was 50%), with the folds parting at 35%.
  - About 50 ms of voiced frication now follows the release. OPEN: the noise is flatter than the NATO 2–4 kHz peak, and the closure still lacks a voice bar.
  - Phone test 33.7% (jh 0 → 12%, ch 6 → 50%); pops 1.28/s.
- **/dʒ/ frication:** the NATO "Juliett" /dʒ/ frication is 2–4 kHz noise with voicing 26 dB down (English /dʒ/ largely devoices through its hiss). Ours was voicing-dominated.
  - Bug: the voiced-stop block, extended to affricates, reset the glottis to touching after the hiss-phase parting. That reset now applies to stops only.
  - /dʒ/ now hisses through the /ʃ/ channel (0.072) with the folds parted to 0.04.
  - Phone test 33.9% (best; jh 25%); pops 1.29/s.

## Eleventh pass: whole-voice quality (owner: "it sounds the same")
The per-consonant fixes of the tenth pass did not change what the owner heard, so this pass measured the whole voice (`tools/voice_quality.py`: jitter, shimmer, H1−H2 and HNR over ARCTIC vowels, middle 60 ms).
- Overall: ours jitter 2.07% / shimmer 9.6% / HNR 2.9 dB; bdl 1.86 / 7.3 / 6.4. Sustained vowels are clean (HNR 21 dB, shimmer 0.8%).
  - Ruled out: aspiration, frication, articulator wander, room, the pitch "ear" gain, tract loss.
- Split by context, the roughness is in **vowels after voiceless consonants**: jitter 5.8%, shimmer 35%, HNR 0.7 dB (bdl 2.1 / 7.1 / 7.9). That is a rough, crackly start to roughly a third of all syllables.
  - Cause 1: the 85% cover-damping cut for folds resting apart. Closing in from a voiceless consonant, the folds swung chaotically. Removed: cover loss is now the same apart or together.
  - Cause 2: the 70 ms adduction adopted against pops in the tenth pass. Adduction is now 20 ms, starting at voice onset (VOT 56–64 ms stressed, about 33 unstressed).
  - Post-voiceless vowels now: jitter 1.9%, shimmer 9.2%. Overall: jitter 1.81, shimmer 8.5, HNR 3.8.
- **Pops, re-examined.** Most events at the 2.2× threshold were mild first pulses. Strong pops (over 3× the 30 ms level) were 0.24/s against bdl's 0.05/s, all at /t k/ releases.
  - Cause: stops released straight towards "fully released", so the tongue crossed contact at full speed and the opening jumped.
  - Stops now open to a 0.6 cm² channel over the release time, then clear over 50 ms. Strong pops 0.05/s (bdl 0.05).
- REJECTED for pops: two-stage adduction (touch, then press). Pops 1.20/s, but post-voiceless shimmer 16–18%: the folds vibrate irregularly while passing through "touching". Lower lung pressure in voiceless consonants: no effect.
- Phone test 33.4–34.2%.

## Twelfth pass: the four-view baseline (owner: "check the four representative forms")
Two measures, both over all four views (enhanced STFT, constant-Q, time-causal Gabor, Wigner–Ville):
- **Four-view identification** (`tools/phone_test.py --views stft,cq,causal,wigner`; higher is better). Every one of our phones is named by a nearest-neighbour vote over human exemplars (bdl a0001–a0040). The human's own a0041–a0060 phones score 54.9%. Owner's target: 40%.
- **Four-view distance** (`tools/compare_arctic.py`, `four_view_summary`; lower is better). Per-phone level-normalised shape error in each view, plus fricative and stop-release level errors.
- `tools/evaluate.py` runs both (about 34 s on the Mini). `tools/four_view_tune.py` tunes on the distance; `tools/identify_tune.py` tunes on identification (objective a0001–a0040, check a0041–a0060).

Steps:
- MEASURED: the start of the pass was identification 33.3%, distance 29.60 (a1–20) / 29.11 (a21–40).
- `four_view_tune` (SH, global noise, vowels IY IH AH EH AE AA): distance 28.62 / hold-out 28.34. Adopted, with aspiration 0.09 and frication 0.6 now the body defaults.
- **Stops burst inside their segment.** In the bdl labels, the burst sits 24/18/43 ms (p/t/k) and 10/8/28 ms (b/d/g) before the label's end (`tools/burst_position.py`); the next phone begins at voice onset. Ours burst at the end, and our /k t/ slices at 75% were 10–20 dB short of the human aspiration (identified as /p/).
  - A stop into a vowel, liquid or glide now opens at end − `post_burst_interval` (place and voicing), and voicing begins at the end (+20 ms into stress).
  - Identification 33.3 → 35.6%. Distance a1–20 unchanged; a21–40 28.34 → 29.06, because the release windows went from too quiet to too loud (b +10, p +6 dB).
- REJECTED: a tissue low-pass on the neck-wall radiation. The /b p/ closures read 6–18 dB louder than bdl's at 0.4–1.2 kHz, and the skin radiation is the only closure path (the seal is complete; reverb was ruled out). Low-pass at 150/300/600 Hz: identification 35.1/35.7/35.9 against 35.8, distance 29.08/28.79/28.63 against 28.62. No gain.
- `identify_tune` pass 1 (W, L, R, ER, Y postures): 36.95 → 39.04% on a1–40, check 32.88 → 33.63.
  - R varea 3 was dropped: +0.14 tuned, −0.9 check.
  - The tuner now refuses any step that costs the check sentences more than 0.3%.
- `identify_tune` pass 2 (vowel and diphthong postures; diphthong end targets are now overridable as `gplace garea glips gprotrusion`): only AW garea, AO vlips and AH vplace moved. 39.25%, check 35.59. The vowel postures are at a local optimum for this test.
- MEASURED: full identification (a1–60) 38.1%; distance 28.72 / 28.73.
- **Weak forms.** Aligned against the bdl labels (vowels matched to vowels), the human said a schwa for 49 of our 71 AH0, 19 of our 27 unstressed UW ("to"), and 12 each of unstressed AE and EH ("and", "them").
  - Function words now take weak forms (`kWeakForms`: a, the, of, to, and, for, from, was, them…), and unstressed AH is AX everywhere.
  - "the" and "to" keep a tense vowel before a vowel; a function word ending its phrase keeps its full form.
  - Identification 38.1 → 38.6%; distance 28.62 / 28.52 (both better).
- `identify_tune` pass 3: AX vlips 1 (a rounder, darker schwa: ah → eh was the largest confusion, 58), D place 15.2. 40.69% on a1–40, check 36.84.
- MEASURED: full identification 39.6%; distance 28.67 / 28.77.
- `identify_tune` pass 4 (global noise, finer vowel steps, stop/fricative/liquid durations, NG place, fricative areas): nothing moved. The phone table is at a local optimum.
- The stop-release mechanics are now plan options (`--release-channel`, `--release-move`, `--burst-scale`, `--voice-lag`). Pass 5 kept every default; the measured burst positions are confirmed (burst-scale 0.6/1.4/1.8: 39.8/38.8/38.5 vs 41.0). Skin 0.15/0.5 and pressure 680/820 were also worse.
- Headroom ranking, (human − ours) × count: /l/ 3.0% of the total (human 84%, ours 22%), schwa 2.4, /d/ 1.6, /t/ 1.4, /z/ 1.2, /h/ 1.1.
- /l/ pass 6: L varea 0.18 (41.13%, check 37.14). Full identification 39.8%.
- REJECTED: unstressed effort. `--unstressed-pressure` 0.6–0.9 gave 27.9–40.0% against 41.1 (phonation continues, with a smaller swing); `--unstressed-gap` 0.005–0.03 did not improve either. Our schwa profile is 5–10 dB above bdl's relative to the utterance, but quieter unstressed syllables lose identification.
- /l/ band-profile fit (`tools/lateral_profile`: closed-form chain matrices, bdl /l/ − /ah/, 250 Hz–5 kHz). Best rms without a side branch 5.70 dB (tongue root at 4 cm, area 0.18, tip contact 14.6, lateral area 0.7, lips 1.7). With a blind 1.5 cm × 1.0 cm² branch behind the contact, 4.84 dB. Both fits are limited by the model's narrow F2 peak against a token-averaged target.
- Joint /l/ postures from the profile fit (pass 8) and single steps around them: nothing beat the current /l/ by more than noise. A single tube cannot make bdl's /l/. OPEN: a lateral side branch behind the tip contact (a new body part).
- **Glottal settings of noisy sounds** are now plan options: `--aspirate-gap` (/h/, 0.09), `--voiced-fricative-gap`, `--final-devoicing` (off).
  - /h/ gap 0.02–0.12: no gain. Final devoicing 0.03–0.14: no gain.
  - Voiced-fricative gap 0.012 → 0.006: 41.27%, check 38.05 (+0.9 on the check sentences). Adopted as the default.
- **MEASURED, end of pass:** four-view identification **40.2%** (start 33.3; human 54.9); distance 28.79 / 28.75 (start 29.60 / 29.11); pops 1.45/s (bdl 1.18); jitter 1.89%, shimmer 8.6%, HNR 4.0 dB (bdl 1.86 / 7.3 / 6.4).

## Thirteenth pass: from the owner's ear (2026-10-12)
Owner, at 40.2%: "helllo my name is bipple and i speap really badly. foxtrot is better. alpha is still elfa. the weap bron fox fence over the lazy dog. is it raining today. sthe stheshells by the seasthore. better." Next check-in: at 60. The human's own ceiling on this test is 54.9%, so the working check-in is 55.
- **Seed honesty.** The tuning ran at one seed. Adding a body part shifted the random sequence, and 40.2% read 37.3–38.4%. Identification is now the mean over three seeds (`tools/evaluate.py`), and `identify_tune` scores two seeds per candidate. MEASURED baseline after the shift: 37.8%.
- **Lateral pocket** (`Part::lateral`, a Helmholtz side cavity at oral junction 44, opened by /l/ while the tip holds contact): REJECTED for now. Openings 0.1–1.0 cm² gave 38.5–38.9% against 39.2 (two seeds). The part stays, closed by default.
  - Bug found on the way: the pocket motor had effort wander, so the "closed" pocket stood slightly ajar. Its wander is now 0.
- The floor junction now uses the generic junction solve, with the jet expansion loss its inline solve lacked. Measured: no difference (37.8% both ways).
- **Lisp ("sthe shells").** bdl /s/ rises from 3 kHz and peaks at 4–5 kHz; ours was about 10 dB short at 3–5 kHz and peaked near 6 kHz. bdl /sh/ is strong from 1.9 kHz; ours started at 1.2 kHz and was flat.
  - The identification test cannot see this (our /s/ scores 84%), so `identify_tune` gained a spectral-shape objective (`PROFILE=s,sh`: RMS of the middle-slice band profile, 1–7 kHz, against bdl).
  - S place 15.6 → 15.3 and SH area 0.12 → 0.10: shape error 5.70 → 4.23 dB (check 7.34 → 5.85). Identification is unchanged (37.9%); the distance is 0.2–0.3 worse.
- **"speak" → "speap", "quick" → "weap".**
  - NATO "Kilo": after the burst, 3.4–6.8 kHz noise sits at −26 to −30 dB re the loudest vowel band for 35 ms. Ours was 15–25 dB lower, behind a brief broadband click at 0.8–2.4 kHz (diffuse, as for a labial).
  - Glottal aspiration ×4 changed nothing above 3 kHz; that noise is the jet through the /i/ channel the dorsum opens into.
  - Velar place next to front vowels was centre + 0.85 = 12.05 cm while /i/ narrows at 13.75 cm. Leaving /i/ the tongue slid back 1.7 cm into the closure and F2 fell as for /p/. `velar_front` is now 1.8 (identification neutral: +0.07 / check −0.15).
- **"Alfa" → "elfa".** Our /ae/ (560/1700/2520 Hz) and /eh/ (580/1760/2570) were one vowel, and /ae/'s F2 was even below /eh/'s; bdl has 581/1885/2759 against 594/1785/2883. The body alone tops out near F2 1700.
  - Vowels can now carry a front-dorsum target (`dplace darea`). `vowel_fit front`: AE dorsum 11.5/0.40 + tongue 9.5/0.80 + lips 7.5 → 590/1910 Hz; EH 12.0/1.0 + 9.5/1.4 + lips 6 → 570/1780.
  - Two-point AE and EH adopted: identification 39.47 vs 39.36 (check 36.39 vs 36.02, two seeds).
- **Bursts against bdl.** A burst-slice objective (`PROFILE=k@4,t@4,p@4,k@5,t@5,p@5`: the 75% and end slices, level kept, 1–7 kHz) showed every stop carrying the same +8 to +14 dB at 0.8–2.4 kHz. A component shared by all places is not frication at the constriction; it masked the place cues.
  - It was glottal aspiration. 0.09 → 0.02 took the burst error from 7.9 to 5.4 dB, steadily. /h/ agrees: shape-and-level error 7.5 → 3.4 dB at 0.03.
  - Identification is flat across 0.04–0.09 (37.9–38.3%). The four-view tuner had set 0.09 to fit levels; the ear and the shapes want less.
  - Adopted: aspiration 0.03; stop release 30 ms (was 22); frication 0.45 (was 0.6); neck-wall skin 0.5 (was 0.3). Burst-and-/h/ error 4.6 dB (check 6.0); identification 38.1%; distance 29.21 / 28.83 (from 29.53 / 29.19).
  - OPEN: /t/ and /p/ still carry +9 to +14 dB at 0.8–1.2 kHz at the burst. /h/ lacks bdl's low voiced energy (−9 to −12 dB at 300–600 Hz), and a narrower /h/ gap made the profile worse.
- **"Papa" → "bipple": planner bug.** Each monophthong re-issued its own tongue posture at 65% of its length (a leftover "relax" with zero drift). After a vowel shorter than about 57 ms (every schwa), that event landed after the next vowel's tongue target and pulled the tongue back. The /a/ of "Papa", after the 56 ms /I/ of "is", sat at F1 431 / F2 2104 Hz: an /I/. Removed. Identification 38.1 → 39.1% (three seeds).
  - The tongue's arrival time (`--vowel-reach`, `--vowel-arrival`) and the jaw for labials (`--labial-jaw`, `--jaw-lead`) are now options. Faster vowel arrival is worse (38.6–39.9 vs 40.4, two seeds); the jaw options were neutral.
- REJECTED: more oral-tract loss (`--oral-loss` 0.0016/0.0028: 38.2/36.5%). Unstressed medial compression 0.15–0.8: no gain. Microphone colouring (`--proximity-db`, `--presence-db`; long-term spectrum ours 3–4 dB short below 315 Hz and 4 dB short at 2.5–5 kHz): no gain; it stays at 0 dB.
- Dynamic range ours vs bdl (frame levels re the 95th percentile): median −9.8 vs −9.0 dB. Not compressed.
- **Pronunciation.**
  - "brown" came out B R OW N (letter rules). Lexicon entries added for the -own/-ow /aw/ words.
  - Bug in the letter-to-sound matcher: a greedy `:` in a left context swallowed the consonant a following `^` needed. "#^:" never matched, so final -y gave /I/ (lazy, body, happy: bdl /i/), and the voiceless -ed, syllabic -le and plural rules failed too.
  - Doubled consonant letters gave doubled phones ("happy" P P).
  - Our flapped /d/ was scored as /t/ (DX → t), while bdl's labels keep /d/ (they never use dx). The score now reports the phoneme.
  - Unstressed medial /I/ in rule-derived words is now schwa (bdl ax: 24 → 19 mismatches).
  - Identification 38.6–39.4% across these (noise level); the distance improved (28.49 on a21–40).
  - OPEN: 42% of ARCTIC word tokens (1561 of 1831 distinct words) go through the letter rules, which misplace stress ("delicate") and vowels ("Philip" /faIlIp/). A real pronunciation dictionary is needed (CMUdict; download needs the owner's permission).
  - OPEN: the 54.9% ceiling is bdl against bdl. A second ARCTIC male would measure what a different real speaker scores (download needs permission).
- **CMUdict** (ref/cmudict, fetched with the owner's permission; sha256 81917843…) is now consulted after our lexicon and before the letter rules. Pronunciation mismatches against bdl's labels fell from about 150 to 64. What remains is mostly our weak forms against his citation-form labels ("and" ae, "was" aa). Identification 38.7%; distance 28.91 / 28.52.
- **The ceiling, calibrated.** CMU ARCTIC rms (a second US male, ref/arctic_rms, fetched with permission; sha256 86fbf03c…), four-view identification:
  - bdl a41–60 against bdl exemplars: 54.9%. rms a41–60 against rms: 51.1%.
  - rms against bdl: 28.7%. bdl against rms: 28.8% (a different real speaker).
  - Ours a41–60: against bdl 38.0%, against rms 29.3%, against both 36.9%.
  - MEASURED conclusion: on the reference speaker we did not tune to, our voice scores as a different human does. Above about 29–30% the test rewards resemblance to bdl, not humanness or intelligibility. 55–60 would mean cloning bdl.
- **A less pressed voice.** Voice quality and the long-term spectrum (ours against bdl):

  | modal rest gap | shimmer % | H1−H2 dB | HNR dB | LTAS 80–160 / 160–315 / 0.6–1.25k / 1.25–2.5k / 2.5–5k (re 315–630) | identification % |
  |---|---|---|---|---|---|
  | bdl | 7.25 | −3.7 | 6.4 | −6.5 / −3.0 / −5.0 / −9.5 / −8.7 | — |
  | −0.025 | 11.3 | −5.6 | 2.1 | −12.0 / −7.5 / −2.5 / −6.5 / −9.5 | — |
  | −0.015 (was) | 8.4 | −4.6 | 4.3 | −10.2 / −5.8 / −4.5 / −10.4 / −12.9 | 38.7 |
  | **−0.008** | 7.4 | −3.8 | 4.6 | −8.1 / −4.1 / −5.1 / −10.3 / −11.4 | 38.9 |
  | −0.004 | 7.25 | −3.5 | 4.7 | −7.4 / −3.5 / −4.3 / −8.8 / −10.4 | 36.2 |
  | 0.0 | 6.9 | −3.3 | 4.3 | −6.6 / −3.1 / −2.1 / −5.8 / −6.7 | 34.6 |

  - Adopted −0.008 (distance 28.71 / 28.37). The modal gap is now `--voicing-gap`.
  - OPEN: HNR stays near 4.5 dB against bdl's 6.4 at every setting.
- **HNR (ours 4.6 dB, bdl 6.4).** Unchanged with frication 0, aspiration 0, wander 0 or 0.5, or earlier vowel arrival. There are no stationary spurious tones in either recording (fine long-term spectrum, 2 Hz bins).
  - The deficit sits after voiceless consonants: median 3.5 vs 6.4 dB (p10 −2.3 vs 1.0); after voiced sounds 5.4 vs 6.7.
  - Adduction timing is now an option (`--adduction-move`, `--adduction-lead`). A faster or earlier adduction lifts post-voiceless HNR to 4.7 and overall to 5.0–5.1 dB, and pops to 0.82–0.94/s. REJECTED all the same: identification 37.7–38.3% vs 38.9, distance 29.2–30.4 vs 28.7 (the vowel onsets after voiceless stops move away from bdl's).
- /l/ shape (`PROFILE=l,l@2,l@4`): the lateral pocket is REJECTED again at openings 0.5–1.5 cm², volumes 0.6–1.6 cm³ and necks 0.25–0.8 cm. Contact 16.7 → 16.2 cm and lateral area 0.45 → 0.3 cm²: shape error 5.85 → 4.94 dB, identification 38.9% (unchanged), distance 28.71 / 28.39. Adopted.
- /t/ burst source isolation ("a tea", band levels re the loudest vowel band). 5 ms after the burst the 0.8–1.2 kHz band reads −36 dB with or without frication, aspiration or neck-wall radiation. The excess over bdl is the deterministic release transient, the tract ringing as the trapped pressure is let go, not a noise source. HYPOTHESIS: the opening's flow onset is too abrupt (the tip's signed root area passes contact at speed). OPEN; the release speed and channel were already swept against the burst profile (30 ms and 0.6 cm² are best).
  - The release speed drives that ring: at 60 ms instead of 30 the 0.8–1.2 kHz band falls about 5 dB, but a uniformly slow release smears the burst.
  - REJECTED: a two-stage release, a slow peel to 0.02–0.1 cm² first (`--release-peel-area`, `--release-peel-time`; tip and dorsum stops). The burst profile error rose from 5.95 to 6.4–7.3 dB. The option stays, off.
- **Within the human spread.** Per phone (cq, all seven slices, 150 Hz–7 kHz, level kept), the distance from our mean fingerprint to the mean of bdl and rms, against the distance bdl–rms. A third human would sit about 0.87× the human spread from the pair's mean.
  - Ours is closer than that for most phones (t 4.3 vs spread 8.0, p 3.8 vs 6.9, iy 4.3 vs 5.5).
  - At or beyond it: z 9.5 (8.3), th 10.1 (9.6), ao 7.5 (5.5), y 7.5 (6.6), sh 7.9 (7.6), ng 7.3 (6.4), m 6.6 (5.5).
  - HYPOTHESIS: what still sounds non-human is not in the average phone spectra but in what they average away: temporal dynamics, token-to-token variation, voice texture, prosody. Only the owner's ear can confirm it.
- **Token-to-token variation** (cq fingerprints, rms distance of each token from its phone's mean, median over phones): bdl 11.0 dB, rms 10.0, ours 9.2. By class (bdl / rms / ours): vowels 10.8 / 9.8 / 9.1; stops 13.6 / 12.3 / 11.6; fricatives 11.6 / 9.7 / 10.2; nasals 9.3 / 9.1 / 7.1; liquids 11.7 / 10.1 / 9.0.
  - Part of the human spread is label noise (forced alignment); our boundaries are exact.
  - HYPOTHESIS: our nasals are more uniform than a speaker's, and the rest is within label noise of human variation.
- **Outliers fitted to both speakers** (`identify_tune` with `BOTH=1`: reference = mean of bdl and rms profiles). Z place 15.0, TH area 0.07 and place 17.0, AO vplace 5 / varea 0.3, M inherent 105. Error 9.09 → 7.18 dB (check 7.49 → 6.01). Identification 38.7%, distance 28.60 / 28.49 (unchanged within noise). Adopted.
- REJECTED: the same two-speaker fit for Y, SH, NG and N. 6.44 → 6.22 dB on a1–40, but the check sentences worsened (6.75 → 6.89). Those phones stay as they were.

## Fourteenth pass: the owner's ear on the fresh samples (2026-10-12)
Owner: "hello my name is papa and i speep ??? baldy / aelpha, bravua, charlie, delta echoa, fox(computer sound)trot / the kluip fox femps over the lazy doll. is is raining today (computer sounding-ish) / she ??error zeashells by the seashore (a bit of computer sounds here and there)". "Papa", "seashore", "lazy", "raining today", "Charlie" and "Delta" now pass.
- **"Computer sound" = under-damped high formants.** /a/ through the closed-glottis probe: F3–F5 (2430, 3270, 3770 Hz) within 3 dB of F1, bandwidths 79–110 Hz; human F4/F5 bandwidths are 150–400 Hz. The per-section loss was frequency-independent: no viscous or thermal boundary-layer loss.
  - Spectrograms of "Foxtrot" against NATO show a steady line near 5 kHz through our vowels and harmonics bright to 8 kHz.
  - New: `diffuse_waves`, a spatial smoothing of the travelling waves each step. One section per step makes it a symmetric (zero-delay) time filter, with loss 2ε(1 − cos ω), so no formant moves.
  - Smoothing raw pressure diverged at once. It now smooths pressure × √area, and only where three neighbouring sections are within 1.5× in area and above 0.1 cm².
  - ε 0.06: F4 bandwidth 110 → 164 Hz and 5 dB down, the 3.8 kHz resonance damped out, F1–F3 unchanged.
  - Measured at ε 0.04 / 0.06 / 0.08 / 0.15: identification 39.1 / 38.8 / 38.6 / 37.6% (was 38.7); distance 28.47 / 28.63 / 28.90 / 29.87 on a1–20 and 28.12 / 28.14 / 28.33 / 29.52 on a21–40 (was 28.60 / 28.49); HNR 5.0 / — / 5.2 / 5.4 dB (was 4.6).
  - Adopted 0.06. OPEN: a 5.05 kHz resonance (−21 dB, 150 Hz bandwidth) survives; it lives where the areas step (teeth, lips) and the smoothing is skipped.
- **"bravua", "echoa": the phrase-final tail.** After a comma the voice ran on 120 ms past the last vowel at nearly full strength (breath decayed to 300 Pa, folds parting at +100 ms), then cut off. Heard as a second vowel.
  - The fade now ends 60 ms past the phrase end (was 120), and the folds part at +30 ms (was +100).
  - The comma pause holds no breath pressure (was 300 Pa). With pressure and open folds, the next phrase's first closure clicked onto moving air.
  - Options: `--phrase-fade-tail`, `--phrase-open-delay`, `--pause-pressure`. Identification 38.7%, distance 28.97 / 28.40 (from 28.63 / 28.14: phrase-final phones changed). Kept for the owner's ear.
- **The 5 kHz line.** A resonance at 5.06 kHz (bandwidth 148 Hz) stayed put with lips, tongue and rigid walls, and vanished when the piriform recesses were shrunk. It is the recesses' pole between their two quarter-wave zeros. The recesses now get the same boundary-layer loss: peak −24 → −31 dB, broad. Identification 38.5%, distance 29.01 / 28.41, HNR 5.2 dB. `piriform_scale` (diagnostic) and the probe's `PIRIFORM` were added for this.
- **"sells sea" → "sells zea": override bug.** The voiced fricative /z/ (and the voiced affricate /dZ/) re-closed the folds 10 ms before its end. That event came after a following voiceless sound had opened them (30 ms early), so the /s/ of "sea" was voiced through. Both now re-close only before a voiced sound. Identification 38.8%, distance 28.43 / 28.41.
- A click inside /z/ before /s/ is the frication surge as the folds part behind the narrow channel (dense noise, not an impulse). Slower parting is worse (`--abduction-move` 0.08/0.12) and so is faster (0.015–0.035, pops 1.02–1.11/s vs 0.96). The default 50 ms stays.
- **"dog" → "doll".** A phrase ending in a stop faded its breath through the vowel, so the final /g/ released no burst and its voiced closure sounded like /l/. Breath now holds to a final stop's release.
  - A phrase-final voiced stop devoices at its release; voiced through, "dog." ended in a puff of schwa.
- **Velar bursts and the teeth.** The jet's teeth boost (×4) now fades past `obstacle_reach` (0.8 cm; e-fold 0.5 cm): a jet 3 cm from the teeth has spread. Bursts and sibilants against bdl+rms: 4.14 → 3.91 dB. Identification 39.0%, distance 28.14 / 28.03 (best so far).
  - OPEN: NATO "Kilo" has −26 to −30 dB of 3.4–6.8 kHz noise for 20 ms after the burst; ours −37 to −50. bdl's /k/ is quieter there than NATO's, and the fits follow bdl.

## Fifteenth pass: variation, and the fricatives again (owner: "still sounds very synthetic"; "all instances of s mangled"; "fe sells"; "speaking into a tube")
- **Variation.** Within each phone, the spread of log duration: bdl 0.36, rms 0.54, ours 0.11 (vowels 0.21; stops 0.09, fricatives 0.06, nasals 0.07). Part of the human spread is systematic context and label coarseness, but our consonants were metronomic.
  - New: `PlanOptions::variation` (default 0.1) with its own seeded generator (`--seed` sets it with the body's). Every phone token draws standard normals for duration (×e^(0.1g)), pitch (accent excursion ×(1 + 0.2g), other vowels ×(1 + 0.03g)), vowel tongue place (+0.3g cm), area (×e^(0.1g)), lips (×e^(0.07g)) and breath effort (×(1 + 0.06g)).
  - Duration log-sd 0.11 → 0.15; fingerprint dispersion 9.2 → 9.3 dB.
  - Identification over three seeds: 39.2% without variation, 38.3% with. Distance 27.90 / 27.98 vs 27.94 / 28.36. Both measures average spectra, so randomness costs a little there; the owner's ear judges it.
- **Affricates.** The /dZ/ hiss used fixed numbers copied from an old /sh/ (14.8 cm, 0.072 cm²); it now takes /sh/'s own place and 0.6 × area. Its start moved from 40% to 30% of the segment (`--affricate-hiss-at`): jh/ch against bdl+rms 9.13 → 8.10 dB (check 6.87 → 5.42).
- **Sibilants.** /s/ rose to 8 kHz with no 4–5 kHz peak, and /sh/ was flat from 1.6 to 8 kHz. Neither the boundary-layer loss nor the teeth reach caused it.
  - The jet-noise corner is now capped at 5.5 kHz (`--corner-ceiling`, was 8 kHz): fricatives against bdl+rms 4.66 → 4.11 dB (check 4.14 → 3.56).
  - REJECTED: a second low-pass stage above the corner (`--frication-stages 2`), and Strouhal 0.14 / obstacle reach 1.2 (marginal).
  - Postalveolar lips 1.2 cm² / protrusion 0.5 → 3.0 / 0: 3.73 dB (check 3.26).
  - Frication 0.6 suits the fricatives (3.59) but not the bursts; it stays 0.45.
  - V area 0.05 → 0.12, voiced-fricative gap 0.006 → 0.012: z/v against bdl+rms 7.51 → 5.70 (check 6.04 → 4.73). Identification 38.5%, distance 28.03 / 27.94.
- Hard onsets after voiceless stops (first 30 ms peak / mid-vowel peak): ours median 0.47 (spikes over 1.3 in 10%), bdl 0.90 (8%). The spike heard in "Papa" is occasional, not systematic.
- **"Speaking into a tube."** Long-term spectrum (re 315–630 Hz), ours against bdl: 80–160 Hz −8.0 vs −6.5; 2.5–5 kHz −11.6 vs −8.7. The F1 region dominated: hollow. Halving the boundary-layer loss changed little.
  - A close-microphone colouring of +2 dB low shelf (250 Hz) and +3 dB presence (3.5 kHz) now brings every band up to 5 kHz within 0.6 dB (total deviation 8.3 → 4.2 dB). 5–8 kHz ends 2.8 dB high.
  - Identification 39.2% (from 38.5); distance 28.38 / 28.17 (from 28.03 / 27.94). Adopted for the owner's complaint.

## Sixteenth pass (owner: "papak", "eofut … a lot of clicks", "over (glitchy)", "sick breathing through tube", "weird r emphasis"; "speak", "jumps" and "sells seashells" now pass)
- **Clicks: wall flow through a closed junction.** Per-sample trace (`say` `FINE="t0:t1"`) of the /tS/ release in "Charlie": at an opening of 0.00005 cm² the junction carried −35 cm³/s, where the gap's limit is 0.16 cm³/s. The wall at the closed section stayed coupled to the junction and dumped its motion through the gap; the output rang at 4.8 kHz from silence within 60 µs.
  - Wall flow at a junction now fades with the opening below 0.01 cm² (`kWallContactArea`). The /tS/ release went from −1.0 (normalised) to −0.02.
  - Identification 39.3%, distance 28.28 / 28.30, pops 0.98/s.
  - Automated click hunts that did not find it: per-pulse high-band impulses (61.5/s of vowel vs bdl 24, unmoved by aspiration, compression, microphone or gap); isolated strong pulses (bdl more than ours); vowel crest factor (ours 3.13 vs bdl 3.10).
- **Final voiceless stops** no longer hold breath to their release: the /t/ of "Foxtrot." burst at full pressure ("mechanical sound at end"). Voiced final stops still hold ("dog").
- **"Alfa" → "elpha".** Our /ae/ in "Alfa" measured F1 398 / F2 1640 Hz (cepstral envelope); NATO's 1037 / 1681. The two-point fit had raised the front of the tongue toward Hillenbrand's raised /ae/ (588 / 1952). The fitter claimed 750 Hz for its open refit, but the rendered vowel still read about 400; only the rendered output counts.
  - Transfer-probe grid: tongue 7.5 cm / 1.5 cm², lips 9, no front-dorsum point → rendered 632 / 1546 Hz.
  - Against bdl+rms it is worse (7.77 vs 5.49 dB: their /ae/ is the raised Midwestern one), and identification/distance are unchanged (39.4%, 28.34 / 28.27). Adopted for the owner's ear.
- The liquids no longer raise the jaw (the tongue makes /l r/); no measurable change.
- **Hard voice onsets after /p t k/.** In 23% of the vowels after a voiceless stop the first pulse exceeds 1.5× the vowel's later peak; bdl none. Trace ("Papa"): aspiration carries about 940 cm³/s through a 0.46 cm² glottis, and adduction drives the folds shut so the first closure cuts about 320 cm³/s within 0.75 ms.
  - Not fixed by adduction speed or lead, lower voiceless breath pressure, or a two-stage onset (`--onset-gap` 0.01–0.03: 17–30%). Those options stay, off. OPEN.
- **"Clicks, sounds of a motor"** (owner, after the wall fix). Every glottal closure stood out above 4 kHz: 61.5 high-band impulses per second of vowel vs bdl 24.0.
  - The viscous throttling of the closing slit sets this. A Poiseuille term ∝ depth / A³ makes the flow taper as the gap closes, if the slit is deep enough. Depth 0.3 → 0.6 / 1.2 / 2.4 cm: 62 → 31 / 23 / 31 per second.
  - At 1.2 cm: jitter 1.87% (bdl 1.86), shimmer 7.9 (7.3), H1−H2 −3.7 (−3.7), HNR 5.7 (6.4). It darkens the long-term spectrum and costs identification (37.4%).
  - Adopted depth 1.0 with microphone colouring +4 dB bass / +6 dB presence: clicks 22.2/s, identification 38.6% (from 39.4), distance 28.69 / 28.44.
  - With the deeper slit, hard onsets after /p t k/ fell 23% → 15% (bdl 0%); pops 0.84/s (bdl 1.18).
  - Re-tested: adduction 35 ms or a two-stage onset gives 11%, but jitter 2.8–3.2% (bdl 1.9) and HNR 4.2–4.9. Not adopted. OPEN.
- **Vowel brightness ("speaking into a tube").** Side-by-side spectrograms of arctic_a0002: bdl's vowels carry soft energy to 7.5 kHz; ours went dark above about 4.5 kHz. Vowel middles, 4.5–7.5 kHz re 0.3–4.5 kHz: bdl −31.9 dB, ours −44.2. (Near-silent frames: ours 7.4%, bdl 5.6%.)
  - The aspiration term held vibrating folds at a 0.06 floor, about 46 dB under parted-fold hiss: a clean jet.
  - New `voiced_breath`: an open-phase turbulence term for adducted folds that fades as they part, so /h/ and stop aspiration keep their bdl-fitted level.
  - 4.5 → −31.6 dB; HNR 5.2–5.4; identification 38.3%; distance 28.39 / 28.25 (from 28.69 / 28.44). Adopted.
- **"speak" → "speet".** /p t k/ after /s/ in the same word were aspirated like word-initial stops; the /i/ of "speak" came out weak and breathy. They are now unaspirated (voice onset at the release), and the /i/ is voiced from the release (F2 ≈ 2300 Hz) before the /k/.
- **"with a".** /dh/ had a 0.01 cm² channel and burst like a stop. Against bdl+rms: area 0.06, place 16.2 → 8.60 → 6.28 dB (check 9.11 → 6.82). Adopted.
- "lazy": the click at the /i/ → /d/ closure is closing frication (300–400 cm³/s through the shrinking tip gap as the oral pressure rises to subglottal), not a numerical fault.

## Seventeenth pass (owner: "work on fixing she, quick, ppaat, really, and lazy. seashore also sounds inauthentic")
- `tools/word_compare.py`: the same word in the same ARCTIC sentence, reference above ours (`SPEAKER=rms` uses all 1132 rms sentences; only a0001–a0060 of bdl are here).
  - "sea": rms /s/ lies at 3.5–7.5 kHz with almost nothing below; ours ran from 1 to 8 kHz.
  - "she", "shore": rms /sh/ is loud at 2–7.5 kHz with a band at 2.5 kHz; ours was weak and smeared down to 0.5 kHz.
  - "quick", "body", "people": our bursts are thin, bright lines; rms's are softer.
- Sibilant source isolation: aspiration and voiced breath contribute nothing; it is frication.
  - Narrowing the jet-noise band (`--frication-band`, rise below corner/3 → corner/1.2) cuts /s/ 1.6–3 kHz by up to 7 dB, but the peak stays at 6.3–8 kHz.
  - A joint fit against bdl+rms (band, S place and area, teeth reach, corner ceiling, SH place and area) sits at a local optimum (3.58 → 3.50 dB, check worse). Nothing adopted. OPEN: what the averaged profiles miss and the spectrograms show.
- Pre-pause vowel length: ours median 109 ms, bdl 110, rms 125.
- **The 1 kHz staircase.** Per-sample traces of a /p/ release show the lip area's growth jumping at every 1 ms control tick, with the output bumping each time. Control targets (areas, fold gap, tension, lungs) step at 1 kHz behind two 0.6 ms smoothing stages, which pass about 7% of each step as a 1 kHz ripple under every movement.
  - Smoothing is now `glide_ms` 1.5 per stage. Identification 38.2% (38.3), distance 28.01 / 28.13 (from 28.57 / 28.30), hard onsets 15 → 11%, pops 0.78/s.
- **Sibilant tops.** Mean spectra over the 60 sentences (re peak): human /s/ falls 3–7 dB by 6.3–7.5 kHz; ours stayed flat (−0.4 / −1.4). Two low-pass stages above the jet corner (`frication_stages` 2): /s/ −1.3 / −3.1, /sh/ between bdl and rms. Identification 38.1%, distance 28.00 / 28.08. Adopted (the earlier rejection used averaged-profile shape only).
- **Fricative level re the adjacent vowel** (median over 60 sentences): /s/ bdl −7.6, rms −2.5, ours −6.8; /sh/ −9.3 / −4.0 / −8.6; /f/ −30.2 / −17.6 / −31.7. bdl is a soft-fricative speaker, and ours had been calibrated to him.
  - The utterance-initial "She" read −23 dB (rms −8). The first vowel's shape is set during the breath, so the /i/ body channel (0.18 cm² over 3 cm) stood behind the /sh/ channel and choked the flow to about 200 cm³/s.
  - New `fricative_body_area` (0.9 cm²): during a tip fricative the body is no narrower than this, and it takes the vowel's own aperture as the fricative ends. "She": −23 → −13.8 dB. (Holding the vowel's onset instead, `--fricative-hold`, changed nothing.)
  - Frication 0.45 → 0.6, putting the levels between the two speakers (/s/ −5.7, /sh/ −8.1). Distance against bdl 29.20 / 28.78; identification 37.7%. Adopted for the owner's "s mangled" / "she → fe".
- **Velar place beside front vowels** now follows the vowel: max(centre + 0.85, min(vowel place − 0.3, centre + velar_front)). A fixed +1.8 cm had put the /k/ of "quick" ahead of its /I/ (11.25 cm).
- **Burst levels** (peak 5 ms frame from 30% of the stop to 10 ms after, re the next vowel's 90th-percentile frame): /p t k/ match the humans (p −11.0 vs −10.0/−10.1; k −8.4 vs −7.4/−6.7). /b d g/ were 1.5–6.5 dB over (b −0.3 vs −6.5/−7.1).
  - Voiced stops now open over 50 ms (`voiced_release_move`; the labial rule keeps the slower time): d −3.1, g −3.2. /b/ stays near −1 dB at any release speed; its loudest frame is the vowel's first pulses (OPEN).
  - Identification 37.8%, distance 29.24 / 28.88.
- **"with a" → "really"** (owner: "the worst damage is now the really before badly"; "your visual decomposition methods aren't diverse enough").
  - Four-view (STFT, constant-Q, time-causal Gabor, Wigner–Ville) and band-level contours (80–600 Hz, 0.6–4 kHz, 10 ms) of rms "with a" (a0137, a0153, a0481) against ours. In all three rms tokens /dh/ is a near-closure: the low band falls 20–39 dB for 30–40 ms, the mid band 16–28 dB.
  - Ours (area 0.06 from the averaged-profile fit) dipped 5–8 dB low and 10–15 mid, so "with a" ran as one vowel. Even 0.01 cm² reached its channel too late in the 57 ms /dh/.
  - DH now closes (area 0) for 75 ms: low −13 to −19 dB, mid −30 to −41 dB for about 60 ms, no clicks at closure or release. Identification 38.1%, distance 28.82 / 28.58 (from 29.24 / 28.88).
  - The averaged-profile fit had pulled /dh/ open; a level contour across the word showed what it misses.

## Eighteenth pass: the whole NATO alphabet (owner: "do the same thing and look over other sounds. extend beyond alpha-foxtrot too")
- `tools/nato_compare.py`: all 26 NATO words against the reference recordings. Level contours in three bands (0.2–0.7, 0.7–2.5, 2.5–7 kHz; 10 ms; dB re the word's peak). Ours is trimmed to its phone timeline, the reference to 35 dB of its peak; DTW on the contours, mean |difference| per band, with a grid figure. (The first version trimmed ours by level and kept the audible breath; the low band started at 80 Hz, but NATO's microphone chain has little below 200 Hz.) Four-view images of the worst words via `tf_analyze` / `tf_view.py`.
- Intervocalic consonant dips (min level re the louder neighbouring vowel, 80–600 Hz / 0.6–4 kHz, 60 sentences, bdl / rms / ours): no consonant is more than 6 dB shallower than both humans; /dh/ −43/−41 now (bdl −30/−43). Mildly shallow: /v/ −12 low (−17/−16), /b/ −19 (−29/−23).
- **Whiskey/Oscar voicing through /s/.** NATO stops voicing at the /s/ onset; ours fades over 30–40 ms. Measured over 60 sentences (low band re vowel at +0/10/20/30 ms): bdl −1/−5/−10/−17, rms −11/−19/−28/−31, ours −9/−16/−27/−37. Ours is already rms-like; NATO's crisp offset is his style. Earlier or faster abduction (`--fricative-abduct-lead`, `--abduction-move`) cost 2 points of identification. REJECTED; the options stay at their defaults.
- **Phrase-initial /b d g dZ/ are no longer prevoiced.** The folds wait slightly apart (0.04 cm) and voice at the release; "Juliett" had opened with 90 ms of hum. Juliett 9.7 → 8.7, Bravo 4.3 → 3.9, Delta 6.1 → 6.5; identification 38.0%.
- **Breath onset** (`breath_lead` 0.08 → 0.04 s before the phrase): voiced openings no longer hum before their phone (Mike 8.2 → 6.2, Victor 9.4 → 7.7, Sierra 8.0 → 6.1, Oscar 11.4 → 9.6). NATO total 6.51 → 5.77 dB; identification 37.4% (seed-noise range), distance 28.79 / 28.58.
- **"kilo is terrible": the pharynx sink pumped air.** Four-view images of "Kilo" showed a broadband click at the /i/ → /l/ transition. Per-sample trace: 900–1000 cm³/s through the constriction as it slid back (sections 34 → 22), while the glottis passed 3–120 cm³/s.
  - Not the walls (rigid: same). Not held waves: `flow_carry`, which rescales the wave difference so a section's flow survives its area change, changed nothing. It stays on, as the physically right update.
  - It was the pharynx volume sink. As dark /l/'s body (7 cm) narrowed the pharynx, L·dA/dt pushed air out with nothing taking it back; the parametric tongue does not conserve its volume.
  - The sink exists for the voiced-stop voice bar, so it now acts only while the mouth is shut (`sink_closed_only`; a front section below 0.02 cm²; fading over 5 ms). The constriction flow now follows the glottis (10–120 cm³/s).
  - The /b/ voice bar is unchanged (−14 to −16 dB). NATO total 5.77 → 5.58 (Lima 7.4 → 4.9). Identification 38.1% (37.4), distance 28.01 / 28.17 (28.79 / 28.58).
- Oscar after the fixes: the low band is dominated by our strong fundamental against NATO's bass-light chain. Our /s/ after /a/ carries 0.5–2 kHz energy NATO's lacks.
  - Averaged over 60 sentences our /s/ lies within the human range below 2.5 kHz (bdl more, rms less).
  - The jet's growth and decay lengths are now options (`--jet-growth` 0.4, `--jet-decay` 1.2 cm). A shorter decay raises the /s/ low band (the source nearer the exit drives the back cavity); 2.0 cm lowers it about 2 dB. Left at default.
- NATO says "OSS-kah" (ICAO spelling, non-rhotic: final F1 about 940 Hz). Our final /er/ reads 480 / 1083 / 2273 Hz. Its F3 is high for an American /er/ (about 1600), so it is weakly r-coloured. OPEN, together with the owner's earlier "weird r emphasis" in "seashore".

## Nineteenth pass (owner: "papat", "alpha bravo charlie echo foxtrot" with Delta lost, "clip brown", "clipek", "gal (sickly)", "she sells … with a spoon in mouth")
- Heard correctly now: "with a", "speak", Hotel through Zulu except Papa, Quebec, Golf. The test sentence is "I speak with a body"; the owner wrote "really badly" for it until /dh/ closed.
- **The owner's own reading as a reference.** `tools/nato_compare.py` with `REF=owner` uses `ref/user_nato.wav` (two passes through the alphabet; word spans found by level, in the script; `PASS=2` for the second). `tools/nato_views.py`: NATO, owner and ours as spectrograms, one row per word. Totals (OBSERVED, before this pass): NATO 5.58 dB, owner 4.90; Papa worst against the owner (7.0; 810 ms against our 420).
- **"papat": phrase-final lengthening applied to the whole duration.** Before a pause, rms vowels run 1.8–2.1× their medial length (AX 40 → 85 ms, IY 100 → 178) and consonants 1.4–1.9×; ours ran 1.25–1.4× (final AX of "Papa," 48 ms, then cut). Klatt's factor stretched only the part above the minimum. Now `final_vowel_scale` 2.0 and `final_consonant_scale` 1.6 multiply the whole duration (`minor_final_share` 0.6 of that before a comma).
- **"papa", "quebec", "golf" in the built-in lexicon** now follow the owner's and NATO's readings: P AH2 P AA1 (both stress the second syllable: the owner holds it 550 ms), K EH2 B EH1 K (no /w/ in either), G AO1 L F (CMU's second variant: the owner's F1 600, F2 below 1100 Hz; our /a/ posture gives F2 1330).
- **"clip" / "clipek": /kw/.** rms's /w/ after /k/ is voiceless and dark (almost nothing above 1.2–1.5 kHz) and voices about 80 ms after the burst. Ours rounded only as the /k/ let go and voiced from the stop's end: a bright voiced glide with bands at 2.5–4.5 kHz.
  - `round_ahead`: the lips take a following /w r/'s rounding at the preceding consonant's closure (non-labials).
  - `sonorant_devoicing` 0.5: after an aspirated stop, voicing waits through half of a following /w r l j/.
  - Identification without either 37.9% against 37.7% with both (seed noise).
- **"spoon in mouth": the space ahead of the tongue tip.** The owner's "She sells…" (`ref/user_three.wav`, 18.2–21.3 s): /sh/ rises steeply at 2.4 kHz to a peak at 2.7 kHz and falls 12–18 dB above 4.5 kHz; /s/ peaks near 5 kHz with a sharp edge near 4 kHz. Ours were both flat white noise above the source filter's lower corner.
  - With a near-white source (corner ceiling 20 kHz, one stage, band 1) our /sh/ was still flat within ±4 dB from 1 to 8 kHz. The tube imposed no front-cavity resonance, and no loss option changed that (wave diffusion, oral loss, wall coupling, skin, viscous depth).
  - The area function during the /sh/ of "She": constriction 0.06 cm² at sections 39–40, then 0.26–0.41 cm² to the teeth. The tongue body, already at the /i/ posture and drawn as one parabola, filled the space ahead of the tip.
  - Geometry: a raised tip frees the space ahead of it from the body (`free_front`, by how engaged the tip is: none when released). The /sh/ front is now 0.65, 2.3, 3.4 cm², then the teeth.
  - RMS difference from the owner's spectra over 1.5–7.5 kHz (mean removed): /sh/ 4.6 → 4.1 dB, /s/ 8.8 → 7.0. /s/ now peaks at 4.75–5.75 kHz and falls above. Identification unchanged (36.2 vs 36.3%). /sh/ level re vowel −7.3 → −4.5 dB (rms −4.0, bdl −9.3).
  - /sh/ is still broad (about 1.9–5 kHz). Not fixed by: constriction place 13.6–14.4, sublingual volume 0.5–8, rounding (lips 1.0–1.5, protrusion 0.6–1.0), channel area 0.1–0.35 (wider falls off above 3 kHz but keeps the low edge), frication band 1.5–2, corner ceiling 4–6.5 kHz, jet decay, obstacle reach. A white source concentrated at the constriction exit (jet decay 0.6) gives a 2.25 kHz peak about 5 dB proud; the incisor source radiates almost flat. OPEN.
- **/f/ level** (`labiodental`, the jet weight where a lip channel strikes the incisors; was a fixed 0.3). /f/ re its vowel: 0.3 → −30.4 dB, 0.6 → −25.2, 1.0 → −20.9 (bdl −30.2, rms −17.6; the owner's /f/ in "Golf" is long and strong). Identification 37.7 / 37.2 / 36.2% (the kNN matches bdl's soft fricatives). Adopted 0.6.
- **"sickly": the end of "Golf."** Cycle-by-cycle glottal flow: the /a/ is steady (8.7 ms), then through the /l/ the flow double-pulses (3.6 + 7 ms: F1 ripple of the narrow /l/ tract on loosened folds) for about 80 ms, then creaks at 79 Hz. The cause was the phrase-final relaxation (folds eased to 0.006 cm, pitch falling to a 90 Hz floor) meeting the /l/ constriction. With the folds held at the modal −0.008 and a 100 Hz floor the cycles stay at 9–10.8 ms. Golf jitter 21% → 7%; Lima unchanged (1.7%).
  - `final_relax_gap` −0.008 (was 0.006), `pitch_floor` 100 (was 90). Identification 37.3%, NATO 5.43 → 5.32, owner 4.79 → 4.93.
- After the pass: identification 37.3% (38.1 before; the louder /f/ costs about 1 point), distance 28.01 / 28.71, pops 0.90/s (human 1.03), fricatives re vowel s −4.8, sh −4.5, f −25.2, th −26.1.
- REJECTED as default, kept as `jet_incoherent` / `--jet-incoherent`: an independent turbulence generator per junction along the jet (eddies are independent beyond a few millimetres; one signal copied over 1.5 cm interferes with itself), rescaled by sum(w)/sqrt(sum(w^2)) to keep the coherent level. Against the owner /s/ 7.0 → 6.3 dB, /sh/ 4.1 → 4.4. White-source /sh/ still shows no front-cavity peak (it tilts up instead). Identification 37.4% (37.3), but distance 28.57 / 29.02 (28.01 / 28.71) and NATO 5.32 → 5.52: the rescale raises the high band where the coherent copies had partly cancelled.

## Twentieth pass: owner, "you have not made a more human sounding model, nor holistically fixed systemic issues that make this voice sound revolting"
- The per-word passes improve meters, not the voice. Meters checked for a systemic fault this pass, all near human (ARCTIC a0001–a0060, vowel middles):
  - 3–7 kHz envelope crest within a pitch period: ours 12.4 dB, bdl 12.5, rms 14.5 (the vertical stripes in our spectrograms are not more impulsive pulses).
  - 3–7 kHz re 0.1–1 kHz: ours −20.3, bdl −15.4, rms −25.4.
  - Waveform crest per 30 ms vowel frame: ours 9.8 dB, bdl 9.4, rms 9.2; all active frames p99 17.5 / 16.0 / 16.6.
  - Yet peak-normalised, our sentences are about 8 dB quieter on average than rms's (a0003: RMS 0.037 vs 0.094): a few outlying peaks, not general peakiness.
- `tools/cross_synth.py`: a listening test that splits the problem by layer. LPC (order 18, 25 ms / 5 ms) separates each signal into excitation (residual: glottal source, noise, timing, pitch) and filter (tract, articulation, radiation). DTW on log spectra aligns ours with rms. Hybrids: our excitation through rms's filters, and rms's excitation through ours. Identity resynthesis error −291 dB. Files in `listen/cross/` (a0001, a0003). PENDING the owner's ear: which hybrid still sounds wrong.
- Owner on the swap (a0003): "humansource is.. closer than oursource. still sounds BAD." The excitation layer (our glottal source, its noise, our timing and pitch) carries most of the damage; the tract and articulation also fail on their own.
- Prosody split. `PlanOptions::timing` and `pauses` (`say --timing FILE`: "duration pitch" per phone, "pause s" per phrase break) replace the planned durations, vowel pitches and pauses after time_phrase and intone_phrase, before the gestures. `tools/prosody_transplant.py` aligns our phones to a human's labels and takes each vowel's F0 from the recording. The rms labels are unusable for this (stressed /eh/ 15 ms, /ao/ 15, /r/ 14 in a0003: automatic-label fragments), so it uses bdl (`SPEAKER=bdl`, a0001–a0060). Files: `listen/cross/bdl_*_{human,planned,human_prosody}.wav`. PENDING the owner's ear: if ours with bdl's prosody still sounds bad, the glottal source itself is the fault.
- Owner: "bdl_a0001_human_prosody is better." Prosody is one systemic fault. Durations-only and pitch-only splits (`listen/cross/bdl_a0001_{durations,pitch}_only.wav`) PENDING the owner's ear.
- **The source: the folds moved wrongly.** Inverse filtering (LPC 18 on the vowel middle, residual integrated) shows bdl's flow rising slowly over most of the cycle and dropping fast, the textbook pulse. Ours sat closed about 60% of each cycle and fired a short, nearly symmetric spike. Over a0001–a0060 vowel middles, same method: open quotient bdl 0.83, rms 0.94, ours 0.73; speed quotient (rise/fall) bdl 3.36, rms 2.38, ours 1.33.
  - `tools/fold_shape.py`: sustained vowels (vowel_probe TRACE) across pitch and adduction, judged on the simulator's own glottal area and flow: open quotient, area and flow skew, closure snap, rebound (a second area hump after closure), jitter, mean and AC flow.
  - The old materials (random search scored on open quotient, jitter and pitch: lower cover 1.03 cm thick and 179,621 dyn/cm, upper 0.0275 cm and 1013): open quotient 0.22–0.37, a 12–16% rebound after each closure, and period doubling at 125 and 155 Hz (jitter 27–36%). At gap +0.02 the area opened fast and closed slowly: the human pulse backwards.
  - New physiological body-cover folds (geometry after Titze & Story 2002 for a modal male: 0.062 / 0.087 / 0.1 g; 0.125 / 0.175 cm; cover stiffness for about 110 Hz per piece; coupling 20,000; damping 0.05; collision 1× cover stiffness). Stretch now thins as much as it stiffens (`fold_mass_exponent` 1.0; at 0.5 the amplitude fell 3× from 100 to 170 Hz and 170 Hz at the modal gap did not vibrate). Results: open quotient 0.54–0.62, flow skew 1.5–2.3, no rebound, jitter 0.3–0.6% (about 3% at 100 Hz), mean flow 76–168 cm³/s, AC 234–403. A −0.02 gap no longer vibrates; −0.013 fails at 140 Hz on the probe.
  - In speech: inverse-filtered open quotient 0.76, speed quotient 1.50 (from 0.73 / 1.33); jitter 3.8% (human 1.86), H1−H2 +4.7 (−3.7), HNR 4.1 (6.4), pops 2.05/s (from 0.90; human 1.03). The pops are at /p/ closures and vowel onsets; vowel-onset envelopes after stops match rms. Faster adduction, a two-stage onset and looser folds all made pops and jitter worse. `abduction_stiffening` (folds tensing as they part) changed nothing measurable; kept at 40 as physiology.
  - The microphone colouring (+4 dB below 250 Hz, +6 dB at 3.5 kHz) had compensated for the pressed source. On the new one it put the fundamental 6 dB and 2.5–8 kHz 3.5–6.5 dB above bdl; now off, every band is within about 4 dB.
  - Identification 37.1 → 33.5%, NATO 5.37 → 6.57 dB. The phone and tract fits were tuned on the old source. PENDING the owner's ear before any refit.
- **Prosody, measured while the owner listens.**
  - `say --plan-only --features FILE`: the planned phones with stress, accent, word and phrase position, without synthesis (0.1 s per sentence).
  - `tools/duration_fit.py`: log duration as phone identity plus context factors and a per-sentence tempo, Huber-weighted least squares against the rms labels (25,874 phones; labels under 20 ms dropped; one sentence in five held out). Held-out log error after removing tempo: learned 0.412 (median |e| 0.236), our hand rules 0.458 (0.279). The rms labels are too noisy for more; the learned lengthening (×1.4–1.55 at phrase ends, stress ×1.14 over unstressed) is close to the hand rules. Not wired into the planner.
  - Speaking rate: ours 11.4 phones per second of speech, bdl 12.1, rms 12.2. F0 spread per sentence: ours 2.05 semitones, bdl 2.05, rms 2.43.
  - **Commas rose; humans fall.** a0001's F0 against bdl: at "trail," and "Steels," bdl falls 140 → 100 Hz while ours rose to 155–160. Over 106 comma pauses in rms, the last vowel falls a median 1.3 semitones (52% fall, 13% rise); sentence ends fall 2.9 (82%). `comma_pitch` 0.93 (was a fixed 1.06 rise).
- **The mouth layer ("humansource still sounds BAD"): F2 was twice as broad as human.** LPC bandwidths of open-vowel middles (a0001–a0060): bdl B1/B2/B3 78/105/185 Hz, rms 148/104/241, ours 83/219/260.
  - The tract itself is not the cause: the closed-glottis transfer of /a/ gives B2 113 Hz, and a sustained probe /a/ through the new folds gives 130. In speech, `voiced_breath` (turbulence of vibrating folds, 4.5, set for the pressed source) filled the spectrum around F2: at 0, B2 148; at 1.5, B2 147, 3–7 kHz re 0.1–1 kHz −26.9 dB (rms −25.4, bdl −15.4), HNR 4.8 (from 4.5; human 6.4). Variation off does not narrow it (236).
  - Adopted 1.5. Identification 31.9% (33.5), NATO 6.57, owner reference 6.54. The identification test's kNN was built on the old source; the owner's ear decides.
- **Fold collision.** Human voice quality (vowel middles a0001–a0060; jitter %, shimmer %, H1−H2 dB, HNR dB): bdl 1.86 / 7.25 / −3.7 / 6.4; rms 1.74 / 4.98 / −3.1 / 7.7. Ours (new folds, breath 1.5): 2.81 / 9.43 / +3.1 / 4.8.
  - The fundamental is about 6 dB too strong. Pressing the folds lowers it but makes them irregular. Planner modal gap −0.002 / −0.005 / −0.008 / −0.011: H1−H2 6.4 / 5.2 / 3.1 / 1.4; HNR 2.9 / 4.0 / 4.8 / 3.1.
  - Collision damping (new `FoldMaterial::contact_damping`; was a fixed critical 1.0). Probe at gap −0.011: at 1.0 vibration fails at 150 Hz and jitters 4.4% at 125; at 0.5 it is regular (jitter 0.2–0.8%, shimmer 0.5–1.7%).
  - In speech, damping 0.5: B2 147 → 115 Hz, B3 237 → 164 (humans about 105 / 185–241), pops 2.31 → 1.86/s, but H1−H2 3.1 → 5.8 and HNR 4.8 → 4.0. Adopted for vowel clarity. Collision stiffness 3× (skew 2.0–2.4 on the probe) raised pops to 2.67/s and B3 to 271; kept at 1×. OPEN: pressed, regular phonation (H1−H2 −3 with HNR 6–8) is out of reach of this fold model.
- **Is pressed phonation the gap? Probably not.** The new folds' open quotient stays at 0.54–0.60 (probe, 130 Hz) whatever the rest gap, lower-piece thickness (0.125–0.22 cm), divergence (−0.01 to 0) or body stiffness (to 180,000). That is normal for male speech at comfortable loudness (flow open quotient about 0.5–0.6; Holmberg et al. 1988).
  - Our H1−H2 excess (+4 to +6 dB against −3) may partly be the human recordings' chain: ARCTIC and NATO carry little below 150–200 Hz (bdl's 80–160 Hz band 3 dB under ours with our microphone flat), which alone lowers a 100–130 Hz fundamental by several dB. HYPOTHESIS, untested.
  - `pressed_flattening` (0 removes the flattened, all-at-once closure of pressed folds): on the probe the folds then vibrate to gap −0.015 (was −0.011). In speech it is worse: pops 1.86 → 2.11/s, B2 115 → 132 Hz, H1−H2 5.8 → 6.4. REJECTED; default 1.
- **Pitch learned from the human reading (`tools/pitch_fit.py`, `data/pitch_en.txt`).** Each aligned vowel's F0 at its middle (rms, 1131 sentences, semitones re the sentence median) is fitted as a sum of factors the planner knows. Held-out vowels:

  | predictor | median error | correlation |
  |---|---|---|
  | learned | 1.15 st | 0.62 |
  | hand rules | 1.50 st | 0.40 |
  | a flat line | 1.38 st | 0 |

  - The hand rules did worse than monotone. Their accents added 22–26 Hz (3–3.5 st) where the speaker's add 0.35–0.43 st, and their declination ran by time rather than through the phrase.
  - Learned structure: phrase start +1.1 → end −1.2 st (by quarters of the phrase's vowels); stressed +0.6, unstressed −0.4; after a voiceless consonant +0.6, after a voiced one −0.16; comma −1.1; sentence end −3.0.
  - `intone_phrase_learned` uses it when `data/pitch_en.txt` exists (`say --rule-pitch` for the old rules); planner held-out correlation 0.60.
  - A statement's last vowel now falls to 0.9× its own pitch, not to a fixed 100 Hz floor (the learned final vowel sits near 78–90 Hz, so a fixed floor would have risen). No pitch below 80 Hz: on the probe the new folds are regular at 80–90 Hz (jitter 0.3%) and break at 75 (14%).
  - Identification 31.8 → 33.0%, distance 28.72 / 29.11 → 28.09 / 28.95, jitter 2.81 → 2.17% (human 1.86); NATO 6.38 → 6.48 (single words).
- **Durations learned (`data/durations_en.txt` from `tools/duration_fit.py`; `learned_seconds` in the planner mirrors its factors; `say --rule-durations` for the hand rules).** The fit's tempo came out 0 (the ridge split the level between phone factors and sentence offsets), which spoke 9.8 phones per second. Set to −0.219 to match rms's 12.2.
  - Identification 33.0 → 33.6%, distance 28.09 / 28.95 → 28.29 / 28.22, NATO 6.48 → 6.11, owner reference 6.55 → 5.97.
  - The hand-rule lengthening options (final_vowel_scale and so on) apply only with --rule-durations.
- **Articulation speed.** Median formant velocity (a0001–a0020): bdl F1 1508, F2 2820 Hz/s; ours 1290 / 2269. With `vowel_reach` (the longest tongue movement into a vowel) at 0.15 s instead of 0.20: 1403 / 2429, identification 34.1%, distance 28.27 / 28.27. 0.11 s gains nothing more (1436 / 2400, 33.7%). Adopted 0.15. The remaining 14% on F2 is OPEN.
- **Pops (2.20/s, human 1.03) now gather in /l/** (after ax 17, s 12, ih 11, eh 8, k 8, p 7, iy 6), often several per /l/ about 25 ms apart. The glottal flow through the /l/ of "Steels" is chaotic: periods 3.7–12 ms, peaks 34–268 cm³/s.
  - Reproduced on the probe with the /l/ body posture (7 cm, 0.18 cm², lips 1): jitter 14% at the current damping (0.05 tissue, 0.5 collision). It is regular at 0.12 / 1.0 (2.2%).
  - In speech, 0.12 / 1.0: pops 2.28/s (no gain), HNR 5.6 → 8.0, shimmer 12.9 → 9.3%, H1−H2 6.6 → 3.5, but B2 111 → 197 Hz, identification 34.1 → 32.0%, distance 28.27 / 28.27 → 29.98 / 29.82. 0.1 / 0.75: pops 2.12, B2 128, identification 33.0%.
  - Neither adopted. The /l/ pops are not the chaotic phonation alone; their source is OPEN.
- **The dark /l/ body was a pharyngeal constriction.** Through the /l/ of "Steels" the narrowest section jumps from the palate (37) to the pharynx (20, 0.2 cm²): `L varea 0.18`, from the band-profile fit. With the tip channel (0.3) in series it starved the folds, which pulsed irregularly (the waveform shows sparse large pulses with ringing between).
  - L body area 0.18 / 0.4 / 0.7 cm²: pops 2.20 / 1.97 / 1.93 per second; identification 34.1 / 34.3 / 34.2%; distance 28.27 / 28.27, 28.43 / 28.36, 28.22 / 28.13.
  - Adopted 0.7: the tongue retracts for a dark /l/ but does not close the pharynx. Frication off left pops unchanged (2.13), so the /l/ pops were not turbulence.
- **Remaining pops (1.93/s) are voice onsets after voiceless sounds** (hh>iy 12, s>iy 9, f>ih 9, t>ay 7, ...) and ih>ng (11). "He turned" (a0009): the first cycles after /h/ run 9.6, 8.9, 7.3, 6.5, 7.2 ms with peaks alternating 115–281 cm³/s, and some pulses carry sharp high-frequency bursts. Not the glottal noise sources (aspiration and voiced breath off: 1.94/s). OPEN: how the soft folds close in from apart (adduction onset dynamics).
  - Onset dynamics do not move it (pops/s; current 1.93): adduction move 0.035 → 2.23, 0.012 → 1.94; abduction stiffening 0 → 1.89, 120 → 2.10; two-stage onset (gap 0.004) → 2.25. REJECTED; defaults kept.
- **REGRESSION, owner: "alpha is gone now, there isnt even a word there."** Since the physiological folds, the planner's extra 0.005 cm pressing of stressed vowels (to −0.013) left them unable to start from rest. Every phrase-initial stressed vowel was silent: Alfa, Echo, India, X-ray (first-vowel level 0.0003–0.0005 against 0.12–0.18). Every batch sent since the fold change carried this. The probe checks had used sustained vowels, never a cold start.
  - `stress_press` (new option) is now 0. Vowels more than 26 dB under their sentence's peak (a0001–a0060, of 797): 21 → 9.
  - The other 9: voicing died during a nasal or /h/ and did not restart. "no idea" (a0057): through both vowels the glottis held a steady 0.016–0.023 cm² gap passing 5–18 cm³/s, with no oscillation. At gap −0.008 the soft folds have a stable still state beside the vibrating one.
  - Modal `voicing_gap` −0.008 / −0.006 / −0.004: silent vowels 9 / 2 / 2; identification 34.0 / 34.4 / 34.1%; pops 1.96 / 2.21 / 2.36 per second; jitter 2.05 / 2.55 / 3.82%. Adopted −0.006: a missing word outweighs a click.
  - NEW GUARD for every fold change: `tools/silent_vowels.py` (vowels under −26 dB re sentence peak) and the vowel-initial NATO words.
- **Full CMU ARCTIC bdl** (owner approved the download): `ref/arctic_bdl_full` (lab and wav linked into `ref/arctic_bdl_full_dl`, release 0.95, 1132 sentences). Labels under 20 ms: bdl 0.0%, rms 4.5%.
  - Refit on bdl, held out: durations log error 0.327 (rms-fitted model on bdl 0.369); pitch correlation 0.55 (0.52), median error 1.00 st (flat 1.26).
  - `data/durations_en.txt` and `data/pitch_en.txt` are now the bdl fits. Tempo needs no correction: 12.1 phones per second, as bdl.
  - Identification 34.6%, distance 28.44 / 28.07, NATO 5.54, owner reference 5.26; silent vowels 6.
- **Voicing that dies and stays dead.** "I have no idea" (a0057): the folds vibrate through "have", then stop about 15 ms into the /n/, as they re-close from the voiced-fricative gap (0.012) to modal (−0.006) under the nasal load, and settle into a still state.
  - Silent vowels by lever: collision damping 0.3 → 3, tissue damping 0.03 → 7, collision stiffness 0.5× → 5, posterior chink 0.03 → 3, chink 0.04 → 2.
  - Chink 0.04 cm² (was 0.02): identification 34.6 → 36.3%, distance 28.44 / 28.07 → 27.83 / 27.53, NATO 5.54 → 5.39, owner reference 5.26 → 4.93. Costs: pops 2.20 → 2.43 per second, H1−H2 8.0 → 8.9 dB, B2 116 → 132 Hz. Adopted (a small posterior chink is common in male voices).
