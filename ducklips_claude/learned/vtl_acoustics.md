# VocalTractLab 2.4 — time-domain acoustics (study notes)

Source: VTL 2.4 `Sources/` (GPL; mathematics only): TdsModel.cpp/.h,
TlModel.cpp, Tube.cpp/.h, Synthesizer.cpp, Constants. CGS units.
OBSERVED = read in code; DERIVED = arithmetic from its constants.

## Discretization and solver
- 103 sections: trachea 23 (1 cm; 4.0, 3.0, then 2.5 cm²), glottis 2 (half
  fold thickness each, 0.3 cm default), pharynx 18 (glottis → nasal port),
  mouth 32 (port → lips); lengths vary with articulation (≈ 0.3–0.6 cm), not
  tied to c·T. Nose 19 × 0.57 cm from the last pharynx section (three-way
  junction); piriform fossa 5 sections from pharynx section 4; 4 sinuses on
  nose sections 8, 9, 11, 12. 107 branch currents. 48 kHz, 12 kHz bandwidth.
- Network: pressures at section centres, flows at boundaries (symmetric T):
  `L = ρ(l/2)/A`, `C = V/(ρc²)`, half-section resistances;
  `(L_A + L_B) du/dt + (R_A + R_B) u = p_A − p_B + S`.
- θ-method, θ = 0.5 (trapezoidal); pressures eliminated
  (`p_{n+1} = D + E(U_in − U_out)`, `E = θΔt/(C + α)`), leaving a symmetric
  107×107 system in currents, re-filled every sample and solved by envelope
  (skyline) Cholesky (≈ 1.9k multiply-adds + 0.7k substitution + 107 sqrt per
  sample, DERIVED). Minimum area 1e-4 cm².
- Per sample: interpolate tube/glottis → glottis uses previous-step pressures →
  set glottal areas → resistances/kinetic/noise from previous flows → solve.

## Losses
- Viscous (frequency independent): `R' = max(0.03·A^−2.6, 0.0065·A^−1·√2250)`
  per cm; `R0 = R1 = 0.5·l·R'`. Crossover A ≈ 0.23 cm² (DERIVED). No
  heat-conduction loss.
- Kinetic at every junction (and side-branch entrances), previous-flow
  linearized: `R1_i −= k·½ρu₁/A₁²`, `R0_{i+1} += k·½ρu₂/A₂²`; k = 1 into a
  contraction or with side branches; into an expansion k = 0.5 (tract),
  0.2 (glottis exit). Compare Borda–Carnot k_BC = 2r/(1 + r).
- Glottal flow emerges from the network (contraction k = 1, glottal viscous,
  exit recovery 0.2) — no separate flow law.
- No dV/dt pumping term. Area rate limits ("softContacts", not at the glottis):
  |ΔA| ≤ 400 cm²/s plus a quartic law `ΔA_max = K((A_old/K)^{1/4} + Δt)^4 −
  A_old`, K = 40e6 closing, 480e6 opening (closing from 0.1 cm² ≈ 7 ms,
  opening to 0.1 ≈ 3.8 ms, DERIVED).

## Walls
- Uniform per unit area: mass 2.0 g/cm², damping 1400 dyn·s/cm³, stiffness
  5e4 dyn/cm³ (resonance ≈ 25 Hz, overdamped ≈ 2.2, DERIVED). Surface
  `S = 2√(πA)·l` (sinus: sphere). `Lw = M/S, Rw = R/S, Cw = S/K`; θ-discretized
  `Uw = α·ṗ + β` — an added compliance. All sections except the glottis.
- **Skin radiation**: wall volume velocities of pharynx, mouth and nose
  sections are summed into the radiated flow (voice bars, nasal murmur).

## Noise sources (calcNoiseSources 911–1474)
- Dipoles: series pressure sources at section inlets. `P_rms = gain·v³·√A`,
  v from the narrowest section's outflow low-passed at 500 Hz; outflow only;
  no Reynolds threshold; gain tapers to 0 between 0.5 and 1.0 cm².
- Placement / gain / fc:
  - glottal (always): obstacle 1.5 cm above the glottis; 0.2e-6 / 1000 Hz
  - nasal (always): inlet of nose section 1; 0.3e-6 / 500 Hz
  - tongue, obstacle at incisors when the jet exit is within 2.5 cm (up to two
    constrictions): fricative 1.0e-6·e^{−d/3}, fc 6000·e^{−d/3}; plosive
    0.6e-6·e^{−d/3}, fc 1500·e^{−d/3}; blended by TS3; zero for laterals
  - tongue, wall source (/x, k/): middle of the section after the jet exit;
    0.5e-6 / 1500 Hz
  - lips: middle of the lip constriction; /f/ 0.7e-6 / 8 kHz, /p/ 0.2e-6 /
    500 Hz, blended by `exp(−(A_teeth − 0.15)/0.10)`
- Each amplitude split linearly between two adjacent inlets with independent
  draws. Gaussian noise → 2nd-order Butterworth low-pass at fc → gain
  ×(1000/fc)^1.5 × 0.6; attack one-pole on √amplitude τ = 3 ms, instantaneous
  release.

## Side cavities, velum
- Sinuses (V, neck l, neck A, f_H DERIVED): 11.3, 0.3, 0.185, ≈1300 Hz;
  6.8, 0.3, 0.185, ≈1680; 33.0, 0.45, 0.145, ≈550; 6.2, 1.0, 0.11, ≈740.
  Neck L = ρl/A, R = 15·8πμl/A² (×15 Hagen–Poiseuille), soft sphere wall.
- Piriform: 3.0 cm, 2.0 cm³ total, 5 × 0.6 cm tapering 1.2 → 0.13 cm², closed.
- Velum: first 4 nose sections ramp quadratically from the port area; no
  transvelar coupling.

## Radiation and output
- In-network load at lips and nostrils: `R = 128ρc/(9π²A)` ‖ `L = 8ρ/(3π√(πA))`.
- Radiated flow = mouth + nostrils + skin flow → output FIR: head–torso
  (default, 180 taps, 48 kHz-specific), piston-in-sphere (64 taps), or the
  plain first difference; then 8-pole Chebyshev low-pass at 12 kHz. Head–torso
  re first difference (DERIVED): +0.6 dB 250 Hz, +1.4 500 Hz, +0.2 1 kHz,
  +2.2 2 kHz, +3.6–3.7 4–6 kHz, +2.3–3.1 10–12 kHz.

## Subglottal
- Lung pressure is an ideal series source at trachea section 0; no lung
  compliance or bronchi. Subglottal F1 ≈ 550, F2 ≈ 1270 Hz (comment).
  Glottis coupling lagged one sample.

## CANDIDATE ideas for our engine
1. Glottis as two tube sections inside the network, with contraction k = 1 and
   exit recovery 0.2, so the flow and its loading emerge.
2. Test VTL's fixed-k recovery against our Borda–Carnot.
3. Per-unit-area walls everywhere except the glottis, and **skin radiation**
   of the summed wall motion.
4. Dipole noise at obstacles (teeth within 2.5 cm, else the next section),
   place-dependent fc/gain, outflow gating, 3 ms attack / instant release.
5. Quartic area-rate limiting near closure (keep our pumping term).
6. In-network R‖L radiation plus a head–torso output filter (redesigned for
   our sample rate).
