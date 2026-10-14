# VocalTractLab 2.4 — glottis models (study notes)

Source: VTL 2.4 `Sources/` (GPL; mathematics only): TriangularGlottis.cpp/.h,
Glottis.cpp, GeometricGlottis2019/2025, TdsModel coupling. Shape values from
the speaker files in `~/Downloads/VTL2.4/API/{M01,JD2,W02}.speaker`. CGS.
OBSERVED = code; DERIVED = arithmetic from constants; HYPOTHESIS = ours.

## Triangular glottis (self-oscillating two-mass) — the model to learn from

- Two lateral DOF, lower and upper mass, rigid translation along the length.
  Rest position `restX` is massless control.
- **Triangular rest gap**: gap = restX at the vocal processes (posterior),
  **0 at the anterior commissure**; posterior width `restX + relX`, anterior
  `relX` (restX ≥ 0); parallel when restX < 0.
- Constants (`.cpp:59–77`): lower T 0.24 cm, m 0.12 g, ζ 0.1, k 80 000,
  contact k 240 000; upper T 0.06 cm, m 0.03 g, ζ 0.6, k 8 000, contact 24 000;
  coupling kc 25 000 dyn/cm; L0 1.3 cm; inlet 0.05, outlet 0.01 cm. All springs
  linear — the only nonlinearity is contact geometry.
- Tension `Q = 1 + (f0 − natural_f0)/dF0dQ` (≥ 0.05); masses /Q, springs and
  contact springs ×Q (frequency ∝ Q, m·k constant, absolute damping
  `r = 2ζ√(mk)` independent of F0); length L0√Q, thickness T0/√Q. **Coupling
  kc kept constant** (Q² scaling failed above 400 Hz). natural_f0 ≈ 131 Hz,
  dF0/dQ ≈ 132 (speaker files), fitted from simulated oscillation, open loop.
  Flutter modulates Q: `F0 += 0.25(F0/100)(sin 2π12.7t + sin 2π7.1t + sin 2π4.7t)`.
- DERIVED in-vacuo modes at Q = 1: 118.7 and 189.6 Hz (ratio 1.6, mode 2
  strongly damped); kc/k1 0.31 at Q = 1, 0.51 at 80 Hz.
- Geometry per mass: both ends open → open length L, mean width
  `backX + frontX`; partly closed → apex `zApex = L(1 + relX/restX)`, open
  length zApex, contact length L − zApex. Section area
  `openLength·meanWidth + chinkArea` (arytenoid chink in parallel on both).
- **Contact**: fraction α = contactLength/L; spring on penetration at the
  midpoint of the closed part, `restXStar = restX(1 − meanContactZ/L)`;
  force `−k_contact·α·(relX + restXStar)`. For abducted folds this is
  `F = k_contact·relX²/(2·restX)` — zero force and slope at first touch,
  quadratic as contact spreads from the commissure; C¹ at both transitions.
- **Damping** `ζ = ζ_static + α·1.0` (up to +1 critical with contact).
- **Pressure forces** from the acoustic network's section pressures (no
  closed-form Ishizaka–Flanagan): `F0 = P_lower·openLength0·T0 +
  ¼(P_sub + P_lower)·inlet·L`, `F1 = P_upper·openLength1·T1 +
  ¼(P_supra + P_upper)·outlet·L`. Pressure acts on the open length only.
- **Integrator**: both masses implicit backward Euler (springs, coupling,
  contact, damping), pressure explicit; one 2×2 Cramer solve per sample at
  48 kHz.
- Shapes (rest_disp lower / upper cm, ary_area cm²): modal 0.00495 / 0.0004 /
  0.0002 (JD2 0.0074 / 0.0067 / 0.0062); stop −0.0297 / −0.0301; h 0.0802 /
  0.07985 / 0.0104; voiced-fricative 0.02035 / 0.01965; voiceless-plosive
  0.0501 / 0.0501; voiceless-fricative 0.10225 / 0.10225. Folds vibrate only
  at ≤ 0.02 cm posterior gap; modal is convergent, upper edge almost adducted.
- Why it does not flap and restarts (OBSERVED mechanisms): commissure pinned
  at zero gap (no whole-length grazing impact), contact damping grows with
  contact fraction, implicit mechanics, a light/soft/heavily damped upper
  mass tightly coupled (robust 1:1 entrainment), convergent prephonatory shape
  with near-closed upper edge so pressure builds and blows the heavy lower mass
  open (no kick needed), constant flutter perturbation, F0-invariant damping.

## Glottal flow (TdsModel)
- The glottis is two sections in the network: L = ρ·0.5·len/A, viscous
  `max(0.03·A^−2.6, 0.0065·A^−1·√2250)`, no soft walls, no area-rate limit,
  floor 1e-4 cm². Kinetic term linearized from previous flow; k = 1 into
  contractions (incl. trachea→glottis), k = 0.2 for expansions out of glottal
  sections (divergent step, glottal exit), 0.5 elsewhere — the asymmetry
  produces the intraglottal pressure difference that drives the masses.
- Staggered coupling: geometry → set sections → previous-step pressures →
  mechanics → acoustics. Interaction (ripple, skew, subglottal coupling)
  emerges from the network.

## Aspiration
- Dipole 1.5 cm downstream of the glottal exit, split between two inlets;
  `gain·v³·√A` (gain 0.2e-6, fc 1000 Hz), v from the upper-glottis outflow
  low-passed at 500 Hz (so modulated by the pulses), forward flow only; 3 ms
  √-domain rise limit, instant fall; 2nd-order low-pass at fc,
  ×(1000/fc)^1.5 × 0.6.

## Kinematic models (not to be adopted — the owner requires feedback-driven sound)
Geometric 2019/2025 prescribe edge motion (half-sine mode on a triangular
rest shape, amplitude 0.0007(L0/1.6)√Ptrans cm, phase lag ≈ 70°, skew by phase
warping, Schoentgen jitter). Recorded for completeness only.

## HYPOTHESES for our three-mass folds
1. Triangular rest gap along the slices (commissure pinned), chink separate.
2. Continuous within-slice contact (α, restXStar), contact stiffness ≈ 3× the
   tissue spring scaled with tension, contact damping +1·critical·α.
3. Implicit mechanical step (3×3 solve per sample), pressures explicit.
4. Upper mass light, soft, heavily damped (ζ ≈ 0.4–0.6), tightly coupled; mode
   2 ≥ 1.5× mode 1. Ours has the upper cover heavier than the lower — reversed.
5. Tension: masses /Q, springs ×Q, coupling constant, length √Q, thickness
   1/√Q; calibrate dF0/dQ from simulation.
6. Asymmetric junction recovery (k = 1 in, 0.2 out of the glottis) so the
   driving asymmetry does not depend on tract inertance (nasals).
7. VTL-like modal rest shape (convergent, nearly adducted); pass through
   0.03–0.06 cm quickly (that is VTL's voiceless region).
8. Flutter/jitter as modulation of tension Q, not phase; phrase ends by
   lowering lung pressure (to ≈ 200 Pa with a 20 ms time constant).
