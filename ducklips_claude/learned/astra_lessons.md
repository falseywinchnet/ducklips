# What to take from Astra's work (owner direction, 2026-10-13)

## Base to take
Not Astra's final tube workshop. The earlier sheet tube, right after we found
that each stage needs **multiple types of segments in combination**, with the
3D periodic wall (terms, response, air + tissue). Snapshot exported read-only
(Astra's repo `~/ducklips` on the Mini untouched):

- `~/ducklips_base/tube-42c7337/web/tube/` — commit 42c7337 "Connect the
  accepted vocal tube to mouth, feedback pockets, and lips" (contains the
  owner-accepted below-mouth stage f83173c). Decline began at e7c7394
  (articulated words: slit valve, coarse oral air graph) and 15797b5 (NATO
  corrections).
- `~/ducklips_base/tube-42c7337/src`, `include`, `reports/TISSUE_SOURCE_STUDY.md`,
  `reports/SOMATIC_CONTROL.md` — commit 36d670d, the C++ tissue source and
  somatic (holding) control.

## The segment, as the owner put it
"A segment of your model is a multiple-piece module that internally is a bit
more sophisticated and therefore physical — three piece. Then you combine these
segments into your total tract and articulate them." Each segment combines a
piece with a distinguishing elasticity, a piece with a distinguishing mass and
a piece with a distinguishing structural support (never the same piece); all
take part in edge-based modal (eigen) amplitude shifts and instantaneous
pressure waves; pieces wander slightly under the effort of holding.

## Sheet tube (42c7337) — what it is
- Wall = 2D periodic sheet (axial × circumferential, seam joined), radial
  motion; lumped M and symmetric K from axial tension, circumferential tension,
  hoop stiffness, an elastic foundation and ring supports. Each angular family
  separates exactly (DFT basis); all modes kept.
- Below-mouth stage: nine sections of different material types — inlet
  expansion, light adapter, firm tissue, heavy soft tissue, cord, heavy soft
  tissue, firm tissue, light adapter, throat expansion (mass/elasticity/support
  per section, e.g. heavy tissue 3.5/1.5/0.2, firm 4.7/0.2/3.2, cord elasticity
  1.5).
- Runtime: implicit-midpoint modal response cached as a sparse spatial
  receptance `R = Φ diag(1/(1 + h²ω²/4)) Φᵀ`; update
  `response = R(2·disp + z)`, `next_disp = response − disp`,
  `next_z = 2·response − 4·disp − z` — no solve at audio time.
- Air: delayed pressure-wave connections axially, diagonally and around rings,
  each with fractional-delay + attenuation FIR; admittance-weighted junctions;
  histories store p·√Y. Tissue: separate delayed velocity-wave connections.
- Variation: per-section seeded drift + smoothed sawtooth + bubble, bounded
  elasticity multiplier (≤ ±10%), 200 Hz control, propagated consistently to
  damping, port impedance and delay.
- Energy only from the lung reservoir; no oscillator. 10–12× real time (WASM).
- Known limits: contact = 25 µm radial floor (no airtight collision/jet),
  uncalibrated materials, no global energy ledger.

## Astra's lessons on organic sound (owner-supplied commentary)
1. Perturb physical causes (relaxation, recruitment, holding effort, pressure)
   and let the sound inherit the variation; no separate jitter/shimmer/tremor.
2. Relaxation needs memory (stress relaxing under held deformation), and
   holding posture needs continuing, delayed, imperfect muscle correction
   (overlapping twitches). New cues change intention but keep deformation,
   velocity, stress memory and twitch history.
3. Perturbations need temporal shape: drift (tens–hundreds of ms) 40%,
   smoothed sawtooth (≈0.35–1.25 Hz) 25%, bubbles (≈160–560 ms) 35%, smoothed,
   bounded (±10% is a maximum, not a constant).
4. A material change must propagate consistently (sheet response, damping,
   tissue-port impedance and delay together).
5. More irregularity ≠ more human. The "sickly"/rattling defect persisted with
   holding variation removed; recurring groups of unequal flow pulses are a
   *structured instability*. Distinguish small variation around a stable
   oscillatory regime from repeated disruption of it.
6. Spatial correlation unresolved: try a small shared (effort/breath)
   component plus local deviations.
