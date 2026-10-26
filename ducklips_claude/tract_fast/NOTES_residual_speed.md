> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. Progress comes only from structure.

# One equation-error residual, cheaper (2026-10-24/25)

Labels as AGENTS.md: GIVEN, OBSERVED, MEASURED, HYPOTHESIS, CANDIDATE,
REJECTED, DECIDED.

GIVEN (owner, 2026-10-24): "we need to fix that cost. dramatically." Brief:
make one canonical equation-error residual (`adjoint.ResidualCache`,
canonical mode, the voice solve's residual) at least 5x cheaper, measured,
bit-identical or with any difference labelled and measured; no constant
tuned; engine changes behind a switch or byte-identical by default.

Environment of every measurement: the M4 Mac mini (4 performance + 6
efficiency cores), `~/miniforge3/bin/python3` (numpy 2.3.5, numba 0.63.1),
Apple clang, Release (-O3). Development tree `/tmp/esf/rs_dev/tract_fast`
(this tree's files), the main tree `~/ducklips/tract_fast` as the baseline
build. Load 1.2-3 unless said (no other job ran; the main session was idle).
Times are process CPU seconds; a "residual" is the main session's
`residual_cost.py` workload generalized: each word at its hinted start,
then every articulator group's lags moved +-1 ms in turn, all events
+-0.4 ms, each event +-1 ms, every residual resumed in one cache from the
one before (the voice Jacobian's pattern). JETBOUND=1 SEALBLEND=1.

## Result (MEASURED)

| word | steps x nodes | before: CPU a residual | after | x |
|---|---|---|---|---|
| Delta `D@alveolar,EH,L@alveolar=0.25,T@alveolar!,AH` | 66,240 x 44 | 0.1551 s | 0.0880 s | 1.76 |
| Papa `P@lips!,AA,P@lips!,AH` | 51,840 x 46 | 0.1260 s | 0.0740 s | 1.70 |
| Alfa `AE,L@alveolar=0.25,F@lips~0.15!,AH` (jet, no seal) | 51,840 x 44 | 0.1504 s | 0.0689 s | 2.18 |
| Echo `EH,K@velar!,OW` | 39,360 x 44 | 0.0933 s | 0.0498 s | 1.87 |
| Golf `G@velar,AA,L@alveolar=0.25,F@lips~0.15!` (jet) | 51,840 x 46 | 0.1613 s | 0.0813 s | 1.98 |
| Oscar `AA,S@alveolar~0.08!,K@velar!,ER` (jet) | 51,840 x 46 | 0.1628 s | 0.1015 s | 1.60 |
| `IY,AA,UW` (no seal) | 39,360 x 43 | 0.0931 s | 0.0390 s | 2.39 |
| `IY,T@alveolar!,AA` | 39,360 x 43 | 0.0925 s | 0.0478 s | 1.94 |
| `UW,D@alveolar,IY` | 39,360 x 46 | 0.0977 s | 0.0522 s | 1.87 |
| `AA,P@lips!,IY` | 39,360 x 46 | 0.0968 s | 0.0548 s | 1.77 |

Every residual of every point identical to the bit to the canonical
residual as it stands (0 of 199 differ; relative difference 0, dB
difference 0), also in exact (non-canonical) mode on Delta and Golf (44
points, 0 differ). Every default engine render byte-identical to the main
tree's build (`bindings/engine_identity.py`: the 12 renders of
`bindings/equivalence.py`'s matrix and the ten words' canonical final
renders, 22 of 22 identical with `cmp`).

One Jacobian of the nine-word voice problem (`voice.jacobian_by_word`,
traverse, groups, cached, secant 1 ms, 9 workers, 392,640 x 35, one task a
word): identical to the bit (max |difference| 0). Wall: the first (cold
caches) 10.82 -> 7.84 s (1.38x), the second (caches warm from the first, a
solve's steady state) 6.98 -> 4.05 s (1.72x); the workers' CPU for both
Jacobians and building the nine words 127.7 -> 86.1 s. The wall is the
longest word's task (Delta's), so it follows Delta's residual.

The 5x target was not reached: 1.6-2.4x a residual, 1.87x on the ten
words' sum (1.229 -> 0.657 s). Section "Why not 5x" gives the measured floor.

## 1. Where a residual's time went (MEASURED, before)

Delta, canonical, resumed (first change at step 9,504 of 66,240; the
engine runs 57,600 steps a residual): 0.171-0.175 s CPU. cProfile tottime
per residual: engine `advance_into` 0.072 (42%), the numba peel
(`StatesPeel.advance`) 0.049 (29%), `_implied_run` own 0.017 (almost all
the hidden window's forward run `_forward_steps`, a numba call cProfile
does not separate), the moving geometry 0.015, `least_area`/`closures`
0.012 (two full passes a residual), other 0.006. Seal-less words: the peel
is `reverse_core.SidePeel` (Alfa 0.081 with its jet).

The engine (the C++ render, Delta's canonical final settings, whole render,
`/usr/bin/time -l` difference of 40 renders and none): 1,240-1,290 ns a step,
34,240 instructions and 5,465 cycles a step, 6.3 instructions a cycle. It
is instruction-bound, not stalled. Section timers (a scratch build with
`mach_absolute_time` around each part, ~14 ns overhead each, ns a step):

| part | ns/step |
|---|---|
| trachea (33 nodes): lining / geometry / air+scatter | 171 / 136 / 264 |
| vocal tube (44 nodes): lining / geometry / air+scatter | 115 / 209 / 296 |
| cord and the ear | 49 |
| jets, control | 10, 7 |

So the brief's "30 ns a node-step" is 16 ns a node-step over the 77 nodes
the engine steps (the trachea's 33 included), at 445 instructions each:
the losses, lining, viscous law, seal blend, kinetic loss, contact and
port handling of the general Tube step. There is no per-step allocation;
ctypes crossings are 7-8 a residual (one per 960-step chunk), negligible.

Divisions (MEASURED, a micro-benchmark): latency 2.3-2.5 ns, throughput
0.30 ns scalar and 0.14 ns a lane in vectors. Series chains of them cost;
independent ones do not.

## 2. What was built (each bit-identical, MEASURED)

### 2.1 The peel, native and time-major (`bindings/peel_core.cpp`, `analysis/peel_native.py`)

`peel_native.NativeStatesPeel` (observer.StatesPeel's interface: advance,
reach, resumable, resumed, boundary) and `NativeSidePeel` (SidePeel without
side port or stop node, `peel_native.peel_moving_record`). The numba peels
run node by node in 512-step chunks; each node's lining recurrence is a
chain of 3-4 divisions a step, latency-bound (26 ns a layer-step). The
native kernel runs all layers as a skewed wavefront: layer j's link runs
d_j = 1 + 2 (n - 2 - j) wavefronts behind the open end and its node one
step behind its link, so within a wavefront no layer reads another's
output. Each expression is the numba kernel's (compiled -ffp-contract=off:
numba fuses no multiply-add; same libm exp and pow; min/max with Python's
selection). The jet's narrowest node per row is kept in a ring filled as
the open end reaches each row (no `_narrowest_steps` pass).

Steps, Delta-like rows (IY,AA,UW's 39,360 x 43 for the bench), ns a
layer-step: numba 26; native time-major scalar 21.6 (275 instructions,
92 cycles, 3.0 a cycle: the reorder window held about two layer-steps'
100-cycle chains); two layers a NEON vector, structure-of-arrays state
12.4; links and nodes in two loops a wavefront (each loop's chain half a
layer-step's) 11.3; wave rings slotted by the writing wavefront (a pair's
ring values one vector access) 10.8; the geometry row 32 wavefronts ahead
prefetched (the layers sweep each row from the last node to the first as
~90 strided streams the hardware does not follow; in-cache 5.6, from
memory 10.9 before) 6.0. The kinetic term per pair: both directions'
coefficients from the areas, the one the flow's sign selects (the engine's
scatter_links keeps both likewise) -0.8 ns. The jet: a pair's two filters
in one vector where both coefficient caches hold their corner (Alfa
0.0745 -> 0.0675 s a residual).

### 2.2 The geometry and its notes (`bindings/geometry_core.cpp`)

`reverse_core._geometry_steps` to the bit, gathering per step what the
callers read in further passes: least area (closures), first area (the
forced waves' impedance), the jet's narrowest node. It reads the command's
tick rows directly (`numpy.repeat`'s 23 MB schedule not made) and keeps
the viscous law by link area (a pure function, as the engine's
`Tube::link_viscous`). `observer.closures_of` takes the kept least areas;
`implied_record_with` / `implied_resume(known=...)` take the closures and
first areas; the runs of sealed steps found with numpy (2 ms as a Python
loop).

MEASURED, Delta's schedule, the kernel alone: 5.64 ns a node-step writing
the full array, 5.50 into a 64-row ring: compute-bound, not
memory-bound (so fusing the geometry into the peel would save little;
section 4). With `pow` replaced (a diagnostic, not kept) 4.73.

### 2.3 The hidden windows' forward run (`bindings/forward_core.cpp`)

`observer._forward_steps` to the bit (per-link jet caches: a pure function
of the corner). Nodes two to a vector where both have linings, links two
to a vector where both seal shares are exactly 1. Bench (Delta's whole
word): 8.3 -> 7.1 ns a node-step; in the residual 16 ms numba -> ~7 ms.

### 2.4 The engine: ENTERINGONLY (a switch; `SustainSettings::entering_from`, `Tube::set_air_frozen`)

In the canonical residual the folds are forced from the last hidden
window's end to the word's end, and only the flow entering is read; the
folds then read nothing of the vocal tube but its first area. With
ENTERINGONLY=step the engine does not step the vocal tube's air and lining
from that step on and refreshes only node 0's geometry
(`Tube::refresh_first_area`, the fused build's follower written with its
two fmas). `adjoint.SharedFieldFlow` sets it to the step after the last
free (NaN) fold step; `RenderCache` never resumes a state frozen at F for
a render that reads the air past F (`RenderCache.serves`), and the engine
refuses that too. `validate_settings` refuses ENTERINGONLY unless
ASPIRATION, FLUX and FOLDTURB are 0: the aspiration's pressure drop (which
reads the vocal tube's first pressure) reaches the folds only through
them. Share of a residual's engine steps frozen: Delta ~30%, Papa ~26%,
Golf ~75%, seal-less words 100%.

Difference, labelled (OBSERVED, not MEASURED to occur): a render whose
vocal tube diverges (non-finite sound) inside a frozen stretch stops the
reference engine there (its later flows zero), while the frozen one goes
on. No such render occurred in the 199 points or the Jacobians.

### 2.5 The engine: two exact fast paths, default on (`src/tube_step.cpp`)

- `Tube::frame_settled` / `refresh_settled_geometry`: when no node's frame,
  support or pushed air moves (the trachea always), refresh_geometry's
  frame part leaves everything as it is; only the areas and links follow
  the linings. The fused build computes the viscous glide's step as one
  multiply-subtract (`fnmsub`, read from its object code): written as
  `std::fma`. (A first version without it moved every render; the identity
  gate caught it.) -3%.
- `Tube::plain_air` / `step_air_plain` / `scatter_links_plain`: a tube with
  its two end ports open and nothing else (the trachea): step_air and
  scatter_links two nodes to a NEON vector, each lane the pinned
  expressions (a vector fma is the scalar std::fma per lane), the drawn
  flows summed in node order after, scalar where a node touches the
  opposite wall. Trachea air 240 -> 150 ns a step; the engine 1,300 ->
  1,170 ns a step.

Both default on: 22 of 22 renders byte-identical (`engine_identity.py`).

## 3. Why not 5x (MEASURED floor)

After all the above, Delta's residual (0.083-0.088 s) is: engine ~0.048 (the
trachea ~0.023 of it: its lining is a tridiagonal elimination, 32 divisions
in series a step, ~125 ns of pure latency; the vocal tube's geometry
0.012, its air and lining while not frozen 0.015; cord and ear 0.004), the
geometry 0.0135, the peel 0.010, the windows 0.007, Python ~0.004. The
engine's work a residual is the physics of every step after the first
changed tick (the folds are a self-oscillator; their state never
re-converges to the base run's), so it cannot be shared across a
Jacobian's columns; bit-identity fixes each operation. A 5x residual
(0.034 s) would need the engine at ~260 ns a step, below what its
trachea chain and the vocal tube's divisions allow without changing the
model.

## 4. REJECTED and negative results (MEASURED)

- Native peel node-major-equivalent time-major scalar: 21.6 ns a
  layer-step (1.2x numba): window-bound (section 2.1).
- The geometry kernel two nodes to a vector with every node's lining
  compiled (exact: compiling an unchanged radius gives the copied values):
  8.7 against 6.8 ns a node-step under the same load: copying a static
  node is cheaper than its two square roots and divisions.
- Interleaving the trachea's and the vocal tube's phases in VoiceBody::step
  (so the trachea's division chain overlaps the vocal lining): exact, 1,177
  against 1,170 ns a step; reverted.
- Fusing the geometry into the peel (a ring of rows, no 280 MB array): not
  built; the geometry is compute-bound (5.50 vs 5.64 ns a node-step) and
  the peel with the row prefetch runs within 7% of in-cache.

## 5. CANDIDATE (not built; for the owner)

- The engine reading the replica's vocal geometry (both compute the same
  frame, linings, links and viscous law; the engine's vocal geometry is
  0.012 s of Delta's residual, the replica 0.0135). NOT bit-identical:
  MEASURED, the engine's GEODUMP rows against the replica's differ from
  step 0 by up to 1.3e-15 relative (area, linings), 3.5e-15 (link
  resistance) and 3e-10 (pushed air), from their initialisations (the
  engine starts from the word's first tube, the replica from the command's
  first tick). A difference of rounding size, which the folds would carry
  as the canonical form's 2.6e-12 is carried; to be measured if wanted.
- The vocal tube's air two nodes to a vector (as the trachea's: seals,
  blend, radiation, a closed velum port, contact): ~5% of a seal word's
  residual by the timers.
- A lean, vectorized vocal geometry refresh in the engine (moving nodes:
  the follower, compile_node, the viscous law): mirrors the fused build's
  multiply-adds node by node.

## 6. How to use (DECIDED within this worker's scope: switches, default off)

Environment, read at import (or set the module flags):
`NATIVE_PEEL=1` (observer.NATIVE_PEEL: the native peels, geometry,
windows; adjoint.ResidualCache reads it when made) and `ENTERING_ONLY=1`
(adjoint.ENTERING_ONLY). Build: `cmake --build build_engine --target
tract_engine peel_core` (`libpeel_core.dylib`, override with PEEL_CORE_LIB).
Off, every caller's results and code path are as before (the engine's two
fast paths in 2.5 are on by default and byte-identical). Reversal: unset
the variables; the engine paths revert by deleting their two `if`s.
Adopting the switches as defaults is the main session's call (`puppet.py`,
`voice.py` not edited here).

Verification commands (Mini): `bash /tmp/esf/rs_identity.sh TREE tag`
(engine renders, `cmp`), the per-word bit check `rs_check_fast.py` (in
`/tmp/esf`, scratch), `rs_voice_jacobian.py` (both Jacobians, compared).

## 7. House-style review (`planning/PROGRAMMING_HOUSE_STYLE.md`)

Scope reviewed: everything written here: `bindings/peel_core.h`,
`peel_core.cpp`, `geometry_core.cpp`, `forward_core.cpp`, `jet_core.hpp`,
`engine_identity.py`; `analysis/peel_native.py`; the added parts of
`src/tube_step.cpp` (frame_settled, refresh_settled_geometry,
refresh_first_area, the plain air), `src/tube.cpp` (three branches),
`src/sustain_run.cpp` (ENTERINGONLY), `include/tract/tube.hpp`,
`sustain_run.hpp`, `bindings/tract_engine.cpp` (one table row);
`analysis/adjoint.py` (RenderCache, SharedFieldFlow, ChainedGeometry,
ResidualCache changes) and `analysis/observer.py` (closures_of,
implied_record_with, the `known` argument, the switch). Checked: types
spelled (no auto), no lambdas or comprehensions, named functions for each
stage, retained state as named fields (peel state layout, rings, caches,
`RenderCache.frozen`), fixed inputs copied once a call (`Peel`), borrowed
pointers' lifetimes stated (ctypes records valid while the Python fields
hold the arrays), operation order and every fused operation written
explicitly, comments citing the measurements.

Recorded deviations: the vector paths evaluate both alternatives of a few
element branches (the kinetic coefficient's two directions, the follower's
the forward run's kinetic coefficient likewise, the jet pair's corner
test falling back whole) where the values are finite and the selected result is
the scalar one; justified by the measurements above, against the style's
rule. `[[gnu::always_inline]]` forces expansion in peel_core.cpp (a
specific measured reason: a call taking the local Peel's address would let
stores alias its fields). Legacy code (`observer.py`, `reverse_core.py`,
`engine.py`, the engine outside the added functions) read, not reviewed,
not claimed compliant.
