# Segment-module tract — design (2026-10-13)

Status: CANDIDATE design for owner review. Sources: owner's brief (2026-10-08/09
and 2026-10-13), Astra's sheet tube at 42c7337 (`astra_lessons.md`), VTL 2.4
(`vtl_*.md`), and the current `voice/` engine (front end, learned prosody,
diagnostics).

## 1. The unit: a three-piece segment module

A module is one length of the vocal system (ℓ ≈ 0.3–1.5 cm, set by anatomy, not
by a wave-step grid). Internally it has three distinct pieces, each with a
distinguishing property; no piece plays two roles.

| Piece | Property | Physical counterpart | Representation |
|---|---|---|---|
| Elastic | elasticity (tension) | mucosa / cover | thin periodic sheet: axial × angular samples on the lumen surface, radial DOF; stiffness from axial tension, hoop tension and the coupling to the body; light; its modes (edge-based eigenwaves) are computed per module |
| Mass | inertia | muscle / soft tissue bulk | one lumped mass per angular sector behind the sheet, coupled to the sheet by distributed springs and dampers; low own stiffness |
| Support | structure | cartilage, bone, ring support, jaw | a frame per module that anchors the mass piece to the skeleton with stiffness k_s; its rest position and stiffness are what articulation and muscle effort set |

- Air pressure acts on the elastic piece; its velocity moves air (wall volume
  velocity), and its summed motion also radiates through the skin (voice bars,
  nasal murmur).
- Edge coupling: neighbouring modules' elastic sheets join along their shared
  ring edge (delayed tissue-wave ports, Astra), and their mass pieces share
  shear springs. Material changes update stiffness, damping, port impedance
  and delay together.
- Holding: the support piece is driven toward its target by muscle effort with
  relaxation memory (stress relaxes under held deformation) and delayed,
  imperfect correction (overlapping twitches), per Astra's somatic control.
  Variation acts on causes (support stiffness, sheet elasticity, effort) with
  temporal shape: drift / smoothed sawtooth / bubble, bounded, a small shared
  (breath/effort) component plus local deviations.
- Module types differ only in their three pieces' parameters: trachea ring,
  adapter, firm tissue, heavy soft tissue, cord, pharynx wall, velum, tongue
  surface, palate (hard support), cheek, lip, nasal wall, sinus wall.

## 2. Solving without expensive solvers
- Per module, the elastic-piece modes come from a small eigenproblem
  (axial × angular samples; angular families separate when a ring is uniform;
  otherwise the module solves its own small dense problem). The one-step
  implicit-midpoint response `R = Φ diag(1/(1 + h²ω²/4)) Φᵀ` is cached
  (Astra) and refreshed at control rate when geometry or effort changes (LDL
  recurrence for tridiagonal families). Audio rate is matrix application only.
- Mass and support pieces: small implicit (backward-Euler) updates per sector.
- Air: per-module air volumes joined by delayed pressure waves (axial,
  diagonal and around the ring; admittance-weighted junctions, histories as
  p·√Y) — the owner's photon-like propagation — with VTL's losses: viscous
  `max(0.03·A^−2.6, 0.0065·A^−1·√2250)` per cm and kinetic junction terms
  (k = 1 into contractions, 0.5 out, 0.2 out of the glottis).
- Glottal coupling staggered one sample (VTL), everything O(N) per sample.

## 3. The cord module (source)
The cord is a module whose two facing elastic pieces collide. It is
feedback-driven (no prescribed motion), following VTL's triangular glottis:
- elastic piece split into lower and upper edges across the fold thickness;
  rest gap triangular along the length, pinned at zero at the anterior
  commissure, abduction sets the posterior gap; posterior chink in parallel;
- continuous contact: fraction α and penetration at the closed part's
  midpoint (quadratic onset `F = k_c·x²/(2·g_post)`), contact stiffness ≈ 3×
  the edge spring, contact damping +1·critical·α;
- upper edge light, soft, heavily damped and tightly coupled (mode 2 ≥ 1.5×
  mode 1); implicit mechanical step;
- tension Q: masses /Q, springs ×Q, coupling constant, length √Q, thickness
  1/√Q; dF0/dQ calibrated from simulation;
- shapes: modal convergent and nearly adducted (posterior ≈ 0.005–0.007 cm);
  voiced fricative ≈ 0.02; voiceless plosive 0.05; /h/ 0.08; voiceless
  fricative 0.10 — vibration only at ≤ 0.02, the 0.03–0.06 range crossed
  quickly;
- pressures from the neighbouring air volumes; aspiration as a dipole 1.5 cm
  downstream, `v³√A`, flow low-passed at 500 Hz, 3 ms rise, instant fall.
- Variation through Q (flutter/drift) and holding effort, not phase.

## 4. Assembling the tract
- Below the cords: lung reservoir → trachea modules → Astra's accepted
  nine-section below-mouth stage (mixed module types around the cord).
- Above: pharynx modules (with piriform pockets), velum junction → nose
  modules (19 sections, Dang & Honda areas, quadratic port ramp) with four
  sinus pocket modules (closed cavities with walls), mouth modules over a
  floor with the tongue as an object above it (owner: "the tongue is an
  articulating object in the mouth, not the bottom"), under-tongue pocket,
  teeth/lip modules, radiation (in-network R‖L at lips and nostrils, then an
  output filter).
- Rest geometry per articulation from VTL's anatomy: run VTL's own geometry
  (GPL) offline as a tool to export, for each vowel and consonant prototype
  (C(a), C(i), C(u)) of a speaker file, the cross-section areas and
  perimeters along the centre line; our modules read these as support targets.
  Minimum-area floors (fricative groove 0.15, lateral 0.20, incisor zone
  0.15 cm²) and closure by contact.

## 5. Articulation and the cue stream
- Front end (existing `voice/`): text → phrases, breaths, words → phonemes →
  learned durations and pitch (bdl, 1132 sentences).
- Gestures (VTL): vowel tier plus lip / tongue-tip / tongue-body / velic /
  glottal / F0 / pressure tiers; consonant targets interpolated from three
  context prototypes by the vowel's (α, β), place dominance, return gestures
  with their own τ; VTL's timing laws (Kleiner vowel onsets, Fietkau
  diphthongs, VOT polynomials, glottal leads, pressure ramps) as priors.
- Motion: 5th-order critically damped target approximation (τ = move/10),
  driving the support pieces' targets; the holding control then maintains
  them imperfectly.

## 6. Milestones (owner listening at each)
1. **One sustained vowel**: lungs → below-mouth stage with the new cord module
   → a fixed /a/ tract of modules. Question for the owner: organic, no rattle?
2. **Vowel set** from VTL geometry; compare resonances with VTL's own transfer
   functions and bdl.
3. **CV syllables** with gestures; intent-vs-outcome audit; cross-synthesis
   and prosody-transplant checks.
4. **Full text pipeline** (hello, NATO, fox, seashore) at ≥ 10× real time.

Not adopted: VTL's kinematic glottis (prescribed motion — the owner's brief
excludes direct drive).
