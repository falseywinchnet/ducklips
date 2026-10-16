# VocalTractLab against our engine: elements, budget, topology (2026-10-15)

The owner, after the six-way "sickly" A/B all sounded "exactly as sickly ...
just like a 1998 tts model": the issue is not that VocalTractLab copied a real
speaker; our engine may be "a tiny demo and not a real engine". We need better
and more consistent articulation, and to study how many elements VTL uses,
its budget against ours, and how much effort goes into the topology.

Sources:

- **VocalTractLab 2.4** (released December 2025; the owner: "VTL wasnt any good until 2025"), `~/Downloads/VTL2.4/API/Developer/Sources` (GPL). Read for structure and numbers only; nothing run or copied.
- VTL3D 1.0.2 (2022) where noted.
- Our earlier VTL 2.4 notes (`vtl_geometry.md`, `vtl_acoustics.md`, `vtl_gestures.md`, `vtl_glottis.md`).
- Our code (`voice/src/geometry.cpp`, `voice/src/planner.cpp`, `tract/src/*`).

Labels:

- OBSERVED: read in code.
- DERIVED: arithmetic on observed numbers.
- ESTIMATE: an order-of-magnitude count, not measured.

## 1. Element counts

### Articulation (the shape of the tract)

| | VTL 2.4 (OBSERVED, `VocalTract.h`) | Ours (`voice/src/geometry.cpp`) |
|---|---|---|
| Model | 3D anatomical surfaces, mirrored | One 1D area profile from a formula |
| Parameters | 19: larynx HX HY; jaw JX JA; lips LP LD; velum VS VO; tongue body TCX TCY, tip TTX TTY, blade TBX TBY, root TRX TRY (automatic from the body), sides TS1–TS3 | 10: tongue place / area, tip place / root, dorsum place / root, lip root, protrusion, teeth gap, jaw closure (+ `free_front`) |
| Tongue | 33 dynamic + 4 static ribs × 11 points = 407 vertices; ellipse body + tip circle + Bézier blade and root, limited against the hull | A parabola (radius = constriction radius + curvature·distance²) |
| Tip and dorsum | Part of the same tongue surface | Two extra parabolic notches, combined with the body by `min()` |
| Other surfaces | Upper cover 23 ribs × 6; lower cover 17 × 5; teeth 25 × 5; lips 9 × 10; epiglottis 4 × 5; uvula 4 × 5; fill 12 × 4; radiation 13 × 6. About 1,000 vertices, about 2,000 triangles per side (DERIVED) | None |
| Centre line | 129 points; a cut plane at each, both profiles sampled at 96 lateral points | None (straight 17.5 cm axis) |
| Area | Integrated from the cross-section profile; contact = tongue within 1 mm of the hull; lateral channels kept | Formula value per section |
| Tube sections | 18 pharynx + 32 mouth (`Tube.h`; VTL3D 1.0.2: 16 + 24), **length of each from the geometry**: the tract lengthens with a lowered larynx or protruded lips | 48 × 0.365 cm, **fixed 17.5 cm**; protrusion only moves the lip start inward |
| Update rate | Geometry rebuilt every 120 samples at 48 kHz (2.5 ms, 400 Hz); tube interpolated per sample (`Synthesizer.cpp`) | 8 kHz (0.125 ms since today) |
| Code | `VocalTract.cpp` 7,776 lines + Surface 1,239 + Splines 1,138 + Geometry 1,186 + AnatomyParams 1,266 = **12,605** | **109** |

### Gestures (how the shape moves)

| | VTL | Ours |
|---|---|---|
| Score | Tiers (vowel, lip, tip, body, velic, glottal, F0, pressure). Each parameter follows a **5th-order critically damped target approximation**, C⁴-continuous, τ 5–40 ms. | Per-phone emitted events (`planner.cpp`, 1,514 lines of rules and exceptions). Each articulator a **2nd-order critically damped motor** with random wander. |
| Consonants | Each consonant has three prototypes (in /a/, /i/, /u/ context), blended by the current vowel's position in the a–i–u plane. Place dominance decides which consonant owns which parameter. | One posture per phone. Special cases added one by one: flap, coda velar, sibilant lips, round-ahead, voiced-stop widening, jaw lead, … |
| Timing | Closure intervals solved against the geometry (gesture edges moved in 2.5 ms steps until the forward model closes exactly when intended); empirical timing laws (Kleiner, Fietkau, VOT polynomials, Kohler τ) | Hand-set leads and moves per class |
| Code | `GesturalScore.cpp` 5,423 + `SegmentSequence.cpp` 937 (VTL3D: 4,439; new in 2.4: pitch gestures, closing gestures, velic time constants, step-response timing, parameter limits) | `planner.cpp` 1,514 + `speaker.cpp` 288 |

### Acoustics and source

| | VTL | Ours |
|---|---|---|
| Network | 23 trachea + 2 glottis + 50 pharynx / mouth + 19 nose + 5 piriform + 4 sinus = **103** sections, 107 currents (`Tube.h`, `TdsModel.h`) | About 140 nodes: trachea 33, vocal 48, nose 30, sinuses 20, piriform 5, under-tongue 4 (DERIVED from `voice_body.cpp`) |
| Solver | Lumped transmission line, θ = 0.5, one banded (skyline) Cholesky per sample: about 3k operations per sample (DERIVED, 2.4 notes) | Travelling waves, one node per step; per node a three-piece wall (cover, body, support) with a tridiagonal lining solve per tube |
| Rate | 48 kHz | 96 kHz |
| Walls | One lumped mass-spring-damper per section (25 Hz, overdamped) | Cover / body / support pieces per node, wander, per-node variation |
| Glottis | **Selected by default: "Geometric glottis 2025"**, kinematic (section 1b). The triangular glottis we ported is a non-default option. | The triangular glottis, self-oscillating |
| Code | TdsModel 2,156 + Tube 617 + glottis models 2,329 = 5,102 | tube 527 + voice_body 669 + cord 220 = 1,416 (+ the sheet engine 2,317) |

### 1b. The 2025 glottis (`GeometricGlottis2025.cpp`, the speakers' default)

A kinematic model after Titze's, "extended and improved by Birkholz". The folds are moved, not blown.

- Each edge (lower, upper) is a cosine at the **commanded F0**; the upper lags by 1.22 rad.
- Amplitude = 0.0007 · (rest length / 1.6) · √(transglottal pressure). The supraglottal pressure is low-passed at 25 Hz. The amplitude fades to 0 between 2000 and 1000 dPa.
- Fold length ∝ √(F0 / rest F0), thickness ∝ 1/length; rest length 1.6 cm, thickness 0.45 cm.
- Contact: below the abduction-set touching amplitude the displacement follows exponentials toward a floor 20% of the amplitude below contact (no hard clip). `pulse_shape` (−1…1) shifts the closing or opening knee by up to 10% of the period.
- Area: integrated over 32 points along the fold, the rest abduction linear from the vocal process, plus a **chink of 0.03 cm²** (default; ours 0.0002 cm²).
- Flutter: Klatt & Klatt (1990) three sines (12.7, 7.1, 4.7 Hz) at 25% by default.
- Jitter: a random walk in phase (Schoentgen 2001 model I), b = 0.28 · 10% per √s.
- The acoustics then load this area (two glottal tube sections) like any constriction.

So in the voice the owner calls natural, F0 is exactly the planned contour, and the pulse is the same shape every cycle except for designed flutter and jitter. Amplitude follows pressure smoothly. None of the self-oscillator's behaviour is there:

- onset F0 overshoot (+10–18% in ours);
- closure failures (/iy/ 60% before today);
- pressure-dependent pitch;
- H1−H2 that no setting moves.

Our cord is the owner's specification ("energy through feedback ... mass-spring parts become oscillators"). A kinematic valve is neither a driven mass-spring nor a feedback oscillator: the lungs remain the only energy source. Whether to adopt one is the owner's decision.

## 2. Budget

ESTIMATE per second of speech:

| | VTL | Ours |
|---|---|---|
| Articulation geometry | 400 frames/s × (129 cut planes × about 2,000 triangles, plus surfaces and limits) ≈ 10⁶–10⁷ operations per frame ≈ **10⁹ per second** | 8,000 ticks/s × 48 sections × about 40 operations ≈ **1.5 × 10⁷ per second** |
| Acoustics | 48,000 × about 3,000 ≈ **1.4 × 10⁸ per second** | 96,000 × 140 nodes × about 150 ≈ **2 × 10⁹ per second** (runs 4.2× real time on one M4 core) |
| Share on articulation | **about 90%** | **about 1%** |

**The budgets are inverted.** VTL spends most of its work, and most of its
code, deciding *what shape the tract is* at each moment. Its acoustics is a
plain lumped line at 44.1 kHz. We spend almost everything on *how the air and
tissue respond* (96 kHz waves, three-piece walls, tissue variation). The shape
they respond to comes from a 109-line formula. A formula cannot keep the
organs consistent with each other, so every word's failure turned into a
special case in the planner:

- "speap", "quip": the tongue slid back as the dorsum closed;
- "baldy": the tip rose while the body stayed back;
- "jumpsh": the lips took the next vowel's rounding;
- "bally": the sealed tube rang;
- /uw/ and /ow/: rounding cannot lengthen a fixed-length tract.

Each is patched now, but the class remains. VTL does not have this class
because its tongue is one object, its jaw carries teeth, lips and tongue,
its tract length follows the anatomy, and its consonants are shapes in
context, not rules.

## 3. What our topology lacks (structural, not tuning)

1. **Variable tract length.** A lowered larynx and protruded lips lengthen the tract (VTL: section lengths from the geometry). Ours is fixed at 17.5 cm. Rounded vowels can only narrow the lips, not lengthen the tube: /uw/ 18%, /ow/ 4% identification.
2. **One tongue.** Body, blade, tip and root as one contour with tangent continuity and a fixed anatomy; the root follows the body (pharyngeal constriction for back vowels comes free). Ours: three independent profile pieces joined by `min()`, the pharynx the same parabola as the mouth.
3. **The jaw as a carrier.** It rotates the lower teeth, lower lip (VTL: lips at half the jaw angle) and tongue base. Ours: a multiplicative narrowing of the front 8–11 cm.
4. **Contact from geometry.** Closure where the tongue is within 1 mm of the hull, lateral channels for /l/, a fricative groove floor. Ours: signed-root areas aimed at zero.
5. **Context consonant shapes and dominance** in place of per-phone rules.
6. **Higher-order, C⁴-continuous target approximation**, with closure timing solved against the forward geometry.

What we have that VTL does not (the owner's specification, kept):

- three-piece tissue per module, with edge coupling and wander;
- feedback cavities (under the tongue, sinuses, subglottal);
- streaming and several voices in real time.

These shape the timbre around a tract that must first be shaped well.

## 4. Proposal: re-assert the topology

**0. The source (owner's decision).** Keep the self-oscillating cord, or put a kinematic glottis in its place (our own, after Titze's four-parameter model, as VTL 2.4 does). Cost: trivial. It would make F0 and pulse shape follow the plan exactly. An A/B on the same sentences would show by ear how much of "sickly" is the self-oscillator.

CANDIDATE, for the owner's decision. Our own design from anatomy and the phonetics literature, not VTL's code.

**A. A midsagittal articulatory model**, about 1,500–3,000 lines.

- Fixed hull from anatomy: palate, alveolar ridge, upper incisors, posterior pharyngeal wall, velum. It can be fitted to the MRI vowel areas we already hold (`tract/data/story_2008_areas.*`) and to bdl's formants.
- The jaw rotates the lower incisors, lower lip and tongue base.
- The tongue is one contour: body circle or ellipse, blade, tip, and a root following the body. Hyoid and larynx height move the glottis end.
- The lips have protrusion and opening.
- A centre line through the channel. At each of about 100 cut lines, the midsagittal distance d becomes an area by region-specific A = α·d^β (published for pharynx, velar, palatal, alveolar and labial regions). Lateral channels and a sibilant groove are explicit features.
- Output per frame: area **and length** per section.
- Cost: about 10⁶ operations per frame. At 400 Hz that is under 1% of one core: affordable.

**B. A variable-length tube.** Our waveguide ties node length to c·Δt. Two ways:

- **B1:** keep the waves. Resample the area function by arc length onto the fixed grid, and move the lip end continuously with a fractional-delay terminal section (the tract's length changes without whole-node jumps).
- **B2:** a lumped line like VTL's (θ-method, one banded solve per sample, O(N), about 3k operations). The three-piece walls become per-section wall impedances.

B1 keeps the owner's wave topology; B2 is the proven path.

**C. Gestures in a physical parameter space.**

- Consonant prototypes in /a/, /i/, /u/ context with place dominance.
- 5th-order (or our own C⁴) target approximation per parameter.
- Closure and release instants solved against the geometry forward model.

The planner's rules shrink to timing and prosody.

**Validation, before the owner listens:**

- formant trajectories through CV and VC transitions against bdl (velar pinch, labial fall, alveolar locus);
- vowel formants against bdl and the MRI areas;
- identification.

Then the owner's ear on the same sentences.
