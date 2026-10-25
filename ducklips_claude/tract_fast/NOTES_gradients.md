> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. For the past week fit was substituted for understanding. Progress comes only from structure: the physics, the direct inverse, and knowing why something works or breaks. Search scaffolding added since 2026-10-19 (arrival grids, time sweeps, joint guards, sane clamps, block polishes, vowel revisits) is labelled REJECTED as a method: it is fit, not understanding. Round-trip dB numbers from it show only that the model can be inverted; they are not progress.

# Exact Jacobians for the speech-inversion search (2026-10-20)

Labels: MEASURED (reproduced by a named command, environment and data
below), OBSERVED (read in a named source), HYPOTHESIS, OPEN, REJECTED.

Environment for every measurement: the M4 Mac mini, `~/miniforge3/bin/python3`
(numba 0.63.1, numpy 2.3.5, scipy 1.16.3), Apple clang 21, shared with two
sibling projects (load average 23 to 34 on 10 cores throughout), so times
are CPU seconds (this process plus sustain's) unless marked wall. Data: the
known-answer test's truth (`sequence_recover`'s plan, bdl5 profile vowels,
the sound rendered by sustain) and a hypothesis near it (log areas +2% rms,
times +0.3 ms rms, seed 11): `analysis/bench_gradients.py`.

## Summary

- MEASURED: the exact Jacobian of the vowel-sequence equation error
  (`sequence_residual`, AA,IY,UW: 36,480 rows x 143 columns) costs 4.5 s CPU
  (median of 5), against 43.6 s CPU for the search's forward differences on
  the reference pipeline: 10x less CPU, 33 in-process residuals' worth.
  Wall 5.3 / 3.7 / 3.1 s with 1 / 2 / 3 native threads.
- MEASURED: it agrees with central differences to the differences' own
  error floor (about 1e-6 relative); the search's own forward differences
  (steps 1e-3 in log area, 2e-5 s in time) are off by 0.1-0.2% (areas) and
  4-8% (times) on this word, and by 65-890% on columns of a four-vowel word
  where a jump lies within the step.
- MEASURED: Levenberg-Marquardt with the exact Jacobian went from -19.5 dB
  to -81.6 dB in 30 residual evaluations; with the search's forward
  differences of the same residual it stalled at -19.8 dB (one seed).
- OBSERVED, then MEASURED: forced folds depend on the hypothesis only
  through the forced wave and the first node's area. A copy of lungs,
  trachea and cord alone (no vocal tube, nose or subprocess) replaces the
  sustain FORCEWAVE render: the in-process residual costs 0.134 s CPU
  against 0.303 s and equals the reference to 5e-12 relative.
- MEASURED: three sources of true discontinuities in the residual (the
  ear's whole-step closure timing; the triangular glottis's contact
  fraction when a rest gap is negative; the moving frame's fast-node
  switch). The exact derivative is the smooth piece's; finite differences
  straddle the jumps.
- OBSERVED, MEASURED: `reverse_core._moving_geometry` glides the links'
  viscous resistance one step behind the engine. Gliding as the engine
  does, the truth's equation error falls from -160.1 to -244.8 dB.

## 1. The residual's computational graph (OBSERVED)

Parameters theta: each phone's node log areas (K x n), the transitions'
times (2K - 1), (nasals: velum ports). Then:

1. `SequenceShape.split`: exp of clipped log areas; times with GAP clamps
   (`max`, so piecewise linear).
2. `keyframes` -> `reverse_core.keyframe_schedule`: rest areas per step,
   sampled at each 1 ms tick, linear between keyframes. Smooth in the
   keyframe areas; in the keyframe times piecewise linear (a kink each time
   a tick crosses a keyframe time, every 1 ms of motion).
3. `_moving_geometry`: per node, the rest radius followed by a critically
   damped 160 Hz frame (a sleep branch below 1e-9 cm / 1e-6 cm/s), its
   pushed air through two 250 Hz poles, its lining's nine constants from
   the radius; per link the viscous resistance, glided over 16 steps, or
   exact while either node moves faster than 0.5 cm/s.
4. `_peel_moving` (no side, no jet): the lips' load, then node by node to
   the glottis; each layer is causal in time except that its pressure reads
   the next step's wave (one step anti-causal per layer); the lining and the
   hop low-pass are the only recursions.
5. w = (p0 - U rho c / A0) / 2.
6. sustain FORCEWAVE (src/voice_body.cpp, VoiceBody::step): with forcing on,
   `returning = forced_wave_` and `cord_.step(trachea.outgoing_right(),
   trachea.area(last), returning, vocal_.area(0), ...)`; forced_flow_ =
   (from_cord_above_ - returning) A0 / rho c. Nothing else of the vocal
   tube reaches the folds. The folds' other inputs (lungs, trachea, body
   variation, tremor, the ear) do not depend on the hypothesis, except the
   ear, which reads the closure instants (whole steps).
7. residual = window(D - U) less its mean, / scale (`fixed_clip`); bend rows.

Structure used for the derivatives:

- a log area of node i in phone k moves only node i's rest schedule, and
  only from the tick its phone's keyframe interval begins: only node i's
  geometry (and links i - 1, i) moves;
- peel layer j reads nodes j - 1 .. j + 1, so node i enters at layer i + 1;
- each layer reads one step ahead: a tangent starting at step s reaches the
  glottis from step s - (layers) on;
- the forced folds are causal: their tangent is zero before their inputs'.

## 2. Strategy and costs

Notation: P parameters (143), T steps (39,360), n nodes (46), m residual
rows (36,480). A residual: geometry + peel 0.10 s, forced folds 0.03 s
(in-process) or 0.2 s (sustain).

| strategy | cost of one Jacobian | memory | status |
|---|---|---|---|
| forward differences (current) | (P + 1) residuals: 43.6 s CPU (reference), 19.3 s (in-process) | small | MEASURED; inaccurate across jumps |
| forward mode, P tangents | P x (tangent work per lane) | T x block per array | chosen |
| reverse mode / adjoint, full J | m sweeps (one per row) | n T states | REJECTED (m = 36,480) |
| reverse mode, gradient only | 3-5 residuals per gradient | ~0.2 GB (n x T x 6 peel states + 300 x T folds states), or checkpoints | OPEN: suits gradient methods, which the ill-conditioned problem punishes |
| Gauss-Newton by CG, matrix-free (J v forward, J^T u reverse) | ~5 residuals per CG step x 20-60 steps | as reverse | OPEN: wins only for P well above ~500 (per-ms trajectories, many phones jointly) |

Forward mode, as built (HYPOTHESIS for the estimate, then MEASURED): each
native pass carries a block of N directions as Dual<N> numbers. Work per
lane relative to one primal: geometry ~0 (only the block's own nodes move),
peel ~0.15 of a primal peel (layer and time starts below), folds ~0.45 of a
primal folds (time starts; the 33-node trachea dominates its state). The
measured split for AA,IY,UW (CPU, three runs): schedule tangent 0.01-0.02
s, geometry + peel tangents 2.1-2.3 s, folds tangents 1.6-1.9 s, the rest
0.04 s.

The hybrid's pieces (all implemented):

1. The forced render as lungs + trachea + cord only (section 1, item 6).
2. Node-local geometry: a block computes its own nodes' frames as duals and
   reads every other node's primal frame as plain doubles; layers none of
   whose nodes the block moves run with double geometry constants
   (`GeometryRead<S, false>`), no tangent arithmetic on them.
3. Layer starts: a block starts the peel at the outermost layer its nodes
   reach, from the primal layer output stored there (the primal pass keeps
   every layer's output: 2 n T doubles, 29 MB).
4. Time starts: a block starts the peel and the folds at its first moving
   tick (less one step per layer), its running state there replayed in
   double (peel) or promoted from a primal checkpoint every 1024 steps
   (folds).
5. Columns ordered times, then each phone's nodes from the lips, in blocks
   that never mix phones (a mixed block would start at the earliest phone);
   blocks of 8 for the peel (narrower blocks start nearer the glottis) and
   16 for the folds.

MEASURED block widths (one thread, CPU s, the phone-major order with the
time starts): peel 2.37 / 2.70 / 3.59 and folds 2.29 / 1.97 / 2.10 at 8 /
16 / 32. Before the time starts and the phone-major order (node-major blocks
of 16, wall at load ~20): peel 3.2 s, folds 2.8 s. MEASURED threads: the blocks are independent; 1 / 2 / 3 threads gave
5.3 / 3.7 / 3.1 s wall, results identical to the bit.

## 3. The prototype

Files (all new):

- `ad/dual.hpp`: Dual<N> (value and N derivatives), elementary functions,
  max/min/clamp choosing by value as std:: does.
- `ad/forced_folds.hpp`: templated copies of tract::Variation, the
  trachea (Tube for a static ring-wall tube with end ports), tract::Cord's
  triangular glottis, and VoiceBody::step's forcing path (lung ramp and
  compensation, causes, tremor, the ear); checkpoint promotion; diagnostic
  traces.
- `ad/forced_folds_api.cpp`: `ff_primal`, `ff_tangent` (block plan,
  checkpoints, threads), diagnostics `ff_primal_feedback`, `ff_trace`.
  Built with the engine's flags (`-O3 -mcpu=native -ffp-contract=fast`).
- `ad/geometry_peel.hpp`, `ad/peel_api.cpp`: templated `_moving_geometry`
  (per node and per link) and `_peel_moving`; `sj_peel` (primal, layer
  outputs, tangent blocks), `sj_geometry`. Built with `-ffp-contract=off`,
  as numba compiles.
- `ad/Makefile`: `make -C ad` builds `build_ad/libtract_ad.dylib`.
- `analysis/gradients.py`: the ctypes wrappers; `keyframe_tangent` (the
  schedule's tangent for any keyframe parametrization: d times, d frames);
  `sequence_keyframe_jacobians`; `ColumnOrder`; `SequenceJacobian`
  (`residual`, `jacobian`), a drop-in for `tube_recover.Problem`'s
  Jacobian; `sustain` (posture_recover.render's run, scratch under
  /tmp/esf/grad, never /tmp/es).
- `analysis/bench_gradients.py`: `check`, `jacobian`, `sweep`, `timing`,
  `descend`, `reference` (each command's output quoted below).

## 4. Validation

Primal (`bench_gradients.py check`, AA,IY,UW), MEASURED:

- native geometry vs `reverse_core.moving_geometry`: identical in all
  1,810,560 values of area, pushed air and link resistance;
- native peel vs numba `peel_moving_glottis`: identical in all 39,360 steps
  of implied flow and first-node pressure;
- native forced folds vs sustain FORCEWAVE, sustain's own first-node area
  and the same waves: not bit-identical. The first difference is 1 ulp at
  step 6 (values near 1e-150 as the lungs start); the oscillator carries it
  to at most 6.4e-10 cm^3/s (2.3e-12 of the peak) by the end. OPEN: which
  operation's rounding differs (`-ffp-contract=on` and `off` differ sooner,
  at steps 3 and 2);
- the residual, native vs reference: 5.2e-12 relative, both -24.36 dB;
  AE,UW,IY,AO: 1.2e-11, both -27.34 dB.

Semantics differences (OBSERVED): sustain lets the folds run free where a
forced wave is not finite (`std::isfinite`); the cord-only copy cannot (no
vocal tube). Only a hypothesis whose peel overflows reaches it, and there
`fixed_clip` walls the residual at 1e3 either way.

Jacobian against central differences of the same residual, MEASURED:

- `jacobian` (all columns, steps 1e-6 log area, 1e-8 s): AA,IY,UW median
  relative error 5.5e-6, 90% 1.6e-5, max 1.2e-4; time columns 3e-7 to
  1.2e-5 (the last leaving's column exactly 0 in both: nothing follows it);
  AE,UW,IY,AO (183 columns) median 6.6e-6, two columns 1.3-1.6%.
- `sweep` (one column, the step varied): the error falls as h^2 and then
  rises as 1/h, the exact column at the floor: log areas 1.3e-6, 1.2e-6,
  7.8e-7 (h 1e-4 to 1e-5); times 3.9e-8 and 1.2e-6 (h 1e-7 s). The two
  four-vowel columns: 2.4e-3 / 7.5e-3 / 1.6e-2 at h 1e-4 / 1e-5 / 1e-6,
  then 8.5e-5 at 1e-7: a jump 6.4e-7 away (section 6, source 3), not the
  derivative.
- `reference`, the search's forward differences on the reference pipeline
  against the exact columns: log areas at step 1e-3: 8.5e-4, 1.8e-3, 1.3e-3
  (at 1e-5: 7.5e-6 to 1.8e-5); times at 2e-5 s: 0.040 and 0.079 (at 2e-7 s:
  5.9e-5, 1.4e-4).

## 5. Speed and use in a search

`bench_gradients.py timing` (AA,IY,UW, medians of 5, CPU s): reference
residual 0.303 (numba peel, sustain render), native residual 0.134, exact
Jacobian 4.50 = 14.9 reference residuals = 33.5 native residuals. The
search's forward-difference Jacobian: 43.6 s CPU on the reference residual,
19.3 s on the native one. AE,UW,IY,AO (183 columns, 51,840 steps): 7.15 s
CPU.

`bench_gradients.py descend` (tube_recover.stage's settings: `lm`, 30
residuals, xtol/ftol 1e-12, gtol 1e-14; the native residual; one seed):

- exact Jacobian: 30 residuals, 21 Jacobians, 110 s CPU; -19.5 -> -81.6 dB;
  times at most 0.146 ms off; areas 0.064 cm^2 rms off (no bend rows here:
  the residual's near-null zigzags, which the search's bend rows hold);
- forward differences at the search's steps: 30 residuals, 17 Jacobians,
  358 s CPU; stalled at -19.8 dB. HYPOTHESIS: its time columns, made of
  jumps (section 6), point the steps wrong; tried once, OPEN for more
  seeds and for the staged search's own start points.

Use: `SequenceJacobian(shape, target).jacobian(x)` in place of
`tube_recover.Problem.jacobian` (same column order as the parameters), or
`least_squares(model.residual, x, jac=model.jacobian, method="lm")`.

## 6. Discontinuities: where the derivative is undefined

Each found by scanning the residual along one parameter (400 points) and
then stage by stage (schedule, geometry, peel, forced wave, folds) across
the largest step; MEASURED on AA,IY,UW unless noted.

1. The ear (VoiceBody::hear_pitch, `pitch_feedback` 0.08 per closure): the
   heard period is a whole number of steps, so a closure instant crossing a
   step boundary moves the tension by a fixed quantum. Along the third
   phone's arrival time, 10 jumps within +-20 us, residual-norm jumps 16-25
   and two of 1.37e3 (a smooth 0.1 us step changes it by 0.048). With the
   ear's gain 0 in the copy (`ff_primal_feedback`), 1 jump remained.
2. The triangular glottis with a negative rest gap (src/cord.cpp,
   compute_geometry): the upper edge's rest plus the body's posture
   variation reached -0.00042 cm at 0.293 s; with a negative rest the whole
   edge closes at once, so the contact fraction jumps 0 -> 1 as back =
   rest + x crosses 0, and the contact damping with it (0.6 -> 1.6 of
   critical). The folds' flow jumped 1.70 cm^3/s, and the jump stayed 1.70
   as the parameter interval was narrowed to 1.25e-11 s (`ff_trace`).
3. The moving frame's fast-node switch (`_moving_geometry`, |rate| > 0.5
   cm/s: the viscous link exact instead of glided): AE,UW,IY,AO node 7, a
   resistance jump of 2.7e-5 at link 6, residual jump 2.3e-5 at log-area
   distance 6.4e-7.
4. Kinks, not jumps: the schedule in the times (each 1 ms tick crossing),
   GAP clamps, `max(area, 1e-4)`, the sleep threshold (its tangent freezes
   at the converged value, as the value does).

The exact derivative is the derivative of the piece the trajectory is on.
At the search's steps a time column straddles several jumps of type 1
(spacing ~4 us against a 20 us step) and a log-area column can straddle
type 2 and 3 jumps; that is the 4-8% and 65-890% errors above.

Handling (proposals, not implemented):

- In the engine (the owner's choice; it changes the forward model and the
  known answers): time the ear's closures between steps (the contact's
  linear crossing), or turn the ear off (PITCHFB=0) in both the truth's and
  the forced renders of an inversion; give the triangular glottis a soft
  contact fraction for a negative rest (over ~1 um of penetration); blend
  the fast-node switch continuously.
- In the search: exact Jacobians with the trust region LM already has
  (a jump is then noise in the accept/reject test, not a wrong direction);
  no differences across jumps; grids only to find a basin.

## 7. A replica mismatch found on the way (OBSERVED, MEASURED)

`reverse_core._moving_geometry` starts the links' viscous glide at step 0
with the countdown at 0. The engine's Tube constructor already called
refresh_geometry once (src/tube.cpp: one glide step taken, countdown 15),
so the engine glides one step ahead and refreshes the law at steps 15, 31,
... Against sustain's GEODUMP the replica's link resistance differs by up to
0.0187 (the first glide) and 2.3e-5 relative later, in 37% of the values.
The copy here has the engine's phase as an option (`GRADIENT_ENGINE_PHASE=1`,
`link_track(..., engine_phase)`): it then matches GEODUMP to 9.4e-16, and the
truth's equation error falls from -160.1 to -244.8 dB (AA,IY,UW, glottal
reference). The fix belongs in reverse_core.py (its owner's file): resistance
target/16 and countdown 15 before step 0. With it adopted, the copy's
default should follow.

## 8. Extending to the full residual

A. Seals (`sealed_residual`, `observer.implied`), OPEN, estimated:
   - pieces: the peel segments with their recorded state (the layer peel
     already holds it), the forward tube through each window
     (`observer._forward_steps`: ~100 lines, as `_forward_moving` with the
     seal), and the folds running free inside the windows. Free folds read
     the engine's own vocal tube, whose state there is the history of
     sustain's whole run, so this needs a templated vocal Tube (static
     admittance, local linings, seal, moving geometry from the node tracks)
     stepped over the whole word for every block: about one more full peel
     without its layer sparsity, roughly +3 to +4 s CPU at P = 143;
   - window edges are whole steps set by area thresholds (seal, HIDE_AREA):
     the derivative holds the windows fixed; their motion is a jump. Smooth
     it by tapering: blend the peeled and the forward-modelled flow over a
     few steps by a continuous function of the narrowest area's distance
     from the threshold.
B. The nose march (`march_residual`), OPEN: the time-major diagonals with a
   Newton solve of the nose share per step (differentiate at the converged
   share by the implicit function theorem: one scalar division per step,
   not through the iterations), `_peel_below`, the nose forward and the
   port's series element; the windows rendered by the whole engine (vocal
   tube, port, nose, folds: the pieces of A plus the nose). Discrete: the
   port's 1e-6 cm^2 shut threshold, the nose's 50 ms quiet skip, the window
   edges; and inside a release's window the folds run free, so their
   neutral phase carries tangents (and jumps) to the word's end.
C. The reversible jet, OPEN: `_jet_source` is differentiable given its fixed
   hashed noise (clamps and the turbulence onset are kinks); the narrowest
   node (`_narrowest`, an argmin over nodes 10..n-1 under 0.6 cm^2) moves
   the source discretely when two areas cross. In the peel the source is an
   explicit one-step-delayed term of the link's own flow: plain Dual
   arithmetic. The gate keyframes go through keyframe_tangent's form.
D. Voiceless folds (FOLDKEYS): `ff_tangent` takes the targets' tangents
   already (`targets_dot`, untested); the Python side needs fold_keyframes'
   tangent in the times (the same piecewise-linear form as the schedule's).
E. The air mass past the last node (tube_recover.Shape's fraction, the
   peel's end_extra): one more dual input to the lips' load.

## 9. Open items

- The folds copy's first 1-ulp difference (section 4).
- Speed: the Dual<N> structs are arrays of structures; a structure-of-arrays
  layout (lanes contiguous across a node's state) might double the
  throughput (HYPOTHESIS: by a count of its operations the peel now runs
  near 8 Gflop/s on one core).
- Reverse mode / Gauss-Newton-CG for large parameter counts (section 2).
- A, B, C above; the engine-side smoothing of section 6 is the owner's call.
- The LM comparison of section 5 on more seeds and on the staged search's
  own start points.
