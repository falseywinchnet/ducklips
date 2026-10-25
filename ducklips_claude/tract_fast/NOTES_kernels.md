> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. For the past week fit was substituted for understanding. Progress comes only from structure: the physics, the direct inverse, and knowing why something works or breaks. Search scaffolding added since 2026-10-19 (arrival grids, time sweeps, joint guards, sane clamps, block polishes, vowel revisits) is labelled REJECTED as a method: it is fit, not understanding. Round-trip dB numbers from it show only that the model can be inverted; they are not progress.

# Inverse kernels: profile, exact speedups, resume (tract_fast, 2026-10-20)

Scope: `analysis/reverse_core.py`, `observer.py`, `marching.py` (edited) and the
new `analysis/kernels_jet.py`, `analysis/bench_kernels.py`. The frozen reference
(`~/ducklips/tract`, `research/ducklips/tract`) was never edited.
`sequence_recover.py`, `tube_recover.py`, `posture_recover.py`, `engine.py` were
not touched; every speedup below reaches them through the kernels' unchanged
entry points, except the resume APIs and `_peel_below(until=)`, which they
must call to use (recipes at the end).

## Method

- **Harness.** `analysis/bench_kernels.py` rebuilds each residual's kernel calls
  with sequence_recover's own helpers (keyframes, plan, velum keyframes, jet and
  fold keyframes, the march's window rule) for five words:
  a `AA,B@lips,AA` (sealed_residual, observer path), b `AA,M@lips^0.47,AA`
  (march_residual), c `AA,S@alveolar~0.08!,AA` (plain peel with the gated jet),
  d `IY,UW,N@alveolar^0.5,IH,F@lips~0.15!,AO,ER,M@lips^0.5` (Uniform: march with
  the jet), e `G@velar,AA,L@alveolar=0.25,F@lips~0.15!` (Golf: observer with the
  jet). Each at the truth, with the consonant's leaving +0.3 ms ("time"), and
  with the last phone's middle node +1e-3 in log area ("node").
- **Both trees in one process.** The reference's modules load as `ref_*` (their
  sibling imports pointed at each other while loading); the two kernel sets are
  evaluated alternately (5 repeats, order alternating), so machine load (10-17
  on 10 cores, 4 P + 6 E cores; numbers below are on that loaded Mini) hits
  both alike. Engine renders are memoized on disk by the bytes of every input:
  identical inputs get identical engine output, and kernel time is measured apart.
- **Equivalence gate.** Every kernel output (schedule, geometry digest, implied
  flow and pressure, revealed mask, modelled sound, nostril share, returned wave,
  port pressure, outer wave, mode, windows, each window render's forced input,
  the final residual) compared with `numpy.array_equal` (NaN equal to NaN).
- numba caches of bench runs go to `/tmp/esf/numba_cache`.

## Profile of the reference (MEASURED, per residual, CPU s, final A/B run)

| word (steps) | kernels | where |
|---|---|---|
| a /aa b aa/ (39 360) | 0.31 | `_peel_states` 0.23 (2 calls), geometry 0.05, schedule 0.011, closures 0.01, window forward 0.011 |
| b /aa m aa/ (39 360) | 0.50 | `_march_run` 0.24 (5), `_peel_below` 0.19 (3), geometry 0.05, velum_port 0.015, schedule 0.014 |
| c /aa s! aa/ (39 360) | 0.30 | `peel_moving_glottis` 0.23, geometry 0.05, schedule 0.015 |
| d Uniform (103 680) | 3.2 | `_peel_below` 2.5 (5 calls, jet), `_march_run` 0.72 (9), geometry 0.15, jet_table 0.10, schedule 0.05, velum_port 0.035 |
| e Golf (51 840) | 0.80 | `_peel_states` 0.67 (2, jet), geometry 0.06, window forward 0.034 |

Where the time goes, by mechanism (MEASURED):
- The node-major peels cost 55-63 ns per node-step without the jet and ~110 with
  it (a, c). Calibration on the same Mini: a dependent divide+add chain costs
  4.3 ns per link, independent divides 0.55 ns: the lining recurrence (four to
  five divisions in series per step) is latency-bound, and every per-step
  `geometry[t, :, 3:]` view cost numba two atomic reference counts.
- The jet: four splitmix hashes, two `exp`, three `sqrt` per link and step,
  although the gate is shut for most of a word (Uniform: open 47 136-66 335 of
  103 680) and the corner sits on its 1500 Hz floor almost always.
- The observer's first stretch was peeled to the end of the sound and thrown
  away past the closure (a: 39 360 + 17 830 steps peeled for 14 411 + 17 830 used).
- march_residual's window renders each re-peel below the port from step 0 to the
  end (Uniform: 5 whole peels with the jet).
- `_march_run`: a serial chain through ~105 nodes per step (oral peel 23, front
  22, nose forward 30, nose peel 30). Dropping each part in turn (b, CPU s):
  all 0.218; without front 0.198, without nose peel 0.170, without oral peel
  0.194, without nose forward 0.162; without all four 0.042.

## Changes, each bit-identical (MEASURED on a-e x truth/time/node: 15 of 15
combinations, every output array equal)

1. **Jet (kernels_jet.py).** White noise from a per-process table filled by the
   same hash (the hypothesis never changes it); the filter's coefficients cached
   by corner frequency; with the gain exactly 0 and a finite velocity the source
   is `(gain * 0.0) * state[3] / sqrt(passed)` (the reference's value to the bit)
   and the Reynolds number, two square roots and an `exp` are skipped. Filter
   state read in place (no views).
2. **Skewed chunks for the node-major peels** (`reverse_core.SidePeel`,
   `observer.StatesPeel`, `marching.BelowPeel`). Node j at t needs node j+1 to
   t+1 only, so after the open end has taken the sound to T, link j is carried
   to T-(n-2-j) and node j is final to T-(n-1-j). 512-step chunks keep the rows
   in cache across nodes. Each value has the same expression and inputs. Chunk
   sweep on c (CPU s): whole sound 0.197, 4096: 0.159, 1024: 0.150, 512: 0.149,
   256: 0.157, 128: 0.157; on a: whole sound 0.141, 512: 0.077-0.088.
3. **No per-step array views** in the peels, `_diagonal`, `_oral_forward` (lining
   and junction read in place; geometry's linings read from column 3).
   Peel on a: 0.089 -> 0.072 s.
4. **The observer peels only what it uses**: each closure's stretch to its record
   and hidden_from plus a tube length (`StatesPeel.reach`), not to the end.
5. **march_windows reads the linings from the geometry itself** (no 75%-of-geometry
   copy: 130 MB for b, 350 MB for Uniform per call) and counts the steps with a
   seal in front of the port once (prefix sums) instead of a reach x nodes scan
   per step.
6. **The peel below the port resumes its last call** (`marching._peel_below`):
   on the very same read-only geometry (moving_geometry now returns read-only
   arrays; nobody writes into them, grep), the next call starts from the chunk
   before its inputs first differ. The memo holds the geometry weakly only.
   (A first version never hit: `numpy.ascontiguousarray(g, float)` of a
   numba-made array returns a new view object every call; fixed.) Uniform's five
   calls: 0.19, 0.14, 0.13, 0.03, 0.0006 s.
7. **Schedules once per tick**: `keyframe_schedule` and `jet_gate_steps` interpolate
   once per 96-step tick and repeat (numpy.interp is per point, the same values);
   `velum_port`'s glide jitted; `jet_table`'s narrowing jitted;
   `closures` reduces the areas in one jitted pass (NaN as numpy's minimum).
8. **moving_geometry** allocates without zeroing (every entry written; the last
   node's link column set to 0 as the zeros had it) and copies a node's lining
   row from the step before when its radius did not move (same inputs, same
   numbers). It is memory-bound (12 doubles per node-step): 1.4-1.6x.

### Result (MEASURED, final A/B, CPU s per residual, medians of 5)

| word | reference kernels | fast kernels | x | notes |
|---|---|---|---|---|
| a /aa b aa/ | 0.31 | 0.13 | 2.4 | observer total 0.25 -> 0.094 |
| b /aa m aa/ | 0.50 | 0.36 | 1.4 | `_peel_below` 0.19 -> 0.094; `_march_run` 0.24 -> 0.23 |
| c /aa s! aa/ | 0.29 | 0.17 | 1.7 | peel 0.23 -> 0.14 |
| d Uniform | 3.2 | 1.4 | 2.3 | `_peel_below` 2.47 -> 0.53, `_march_run` 0.71 -> 0.68, jet_table 0.10 -> 0.005 |
| e Golf | 0.80 | 0.22 | 3.6 | observer 0.72 -> 0.18 |

"Kernels" = schedule + geometry + velum + gates + peel / observer / march_windows
(the march's last peel below the port is outside it: b 0.06, d ~0.0006 resumed).
Engine renders are the other agent's; at these speeds a residual's kernels are
now about one to three engine renders' worth.

## REJECTED (MEASURED)

- **fastmath** (all flags, or contraction only) on the peel: not bit-identical,
  max relative deviation of the implied flow 2.2e-14 (a) / 1.5e-13 (c) with all
  flags, 6.3e-15 / 1.6e-13 with contraction only; gain 0-27%. At a truth whose
  equation error is ~1e-15 relative it would raise the floor; no switch kept.
- **Time-major staggered peel** (all nodes per step, skew 2, so nodes in a step
  are independent): exact, but 0.077 s against the chunked 0.066 on a, 0.168
  against 0.126 on c.
- **Fused node groups** (2-6 nodes in one loop, lagging two steps; the three
  passes fused per step): exact, 2x slower than the separate-pass chunks (a 0.15-0.19
  against 0.079). The separate pass-1 and pass-2 loops pipeline; fusing lengthens
  the serial body.
- **Nose constants hoisted** out of the march (the nose is fixed: impedances,
  admittances, Borda-Carnot coefficients, radiation, shelf gain once): exact, no
  gain (`_march_run` 0.211-0.227 s against 0.220-0.227); latency-bound. Reverted.
- **prange**: the peel's nodes form a chain (node j needs node j+1), the march is
  a chain in time, and the search already runs Pool(10) on 10 cores; nothing
  was parallelized.

## OPEN

- `_march_run` (b 0.23 s, Uniform 0.68 s) is the largest remaining kernel. Its
  per-step chain could overlap the front and the nose peel (leaves) with the
  next step's oral peel by interleaving their node loops: a 1.5-2x estimate, not
  tried (the march's modes and hand-overs make it a large change).
- The geometry's 12 doubles per node-step (172 MB for 0.41 s, 464 MB for
  Uniform) bound geometry construction and resume by memory traffic; a compact
  form (linings by radius) would need every kernel and caller to change.
- `waves` from `march_windows_resume`: see the resume gate.

## Resume (locality): APIs and gate

Contract for every API: the new hypothesis's inputs equal the record's before
step `first_changed` (the schedule's rows, the sound, the jet gate, the port;
for engine callbacks also the fold and jet keyframes' effect: the engine must be
causal, which it is). The bench's `first_changed` is the minimum of the first
differing schedule row, gate, port, and the tick of the fold keyframe before
the first differing one.

- `reverse_core.moving_geometry_record(rest)` -> (geometry, GeometryRecord), state
  every 960 steps; `moving_geometry_resume(rest, record, base, first_changed)`
  (copies the base rows: memory-bound); `GeometryWorkspace(base, record).resume(rest, k)`
  writes only the rows from the boundary into one buffer (valid until its next
  resume; writable, so the `_peel_below` memo never trusts it).
- `reverse_core.peel_moving_record(...)` -> (flow, pressure, SidePeel);
  `peel_moving_resume(record, geometry, flow_out, first_changed, gate)`.
- `observer.implied_record(...)` -> (flow, pressure, revealed, modelled,
  ImpliedRecord); `implied_resume(record, geometry, sound, free_flow,
  first_changed, jet_gate)`: closures whose peel and window lie before the
  change are taken from the record, **free_flow (an engine render) is not called
  for them**; the stretch the change falls in resumes its peel from the last
  chunk boundary; the rest runs again.
- `marching.march_windows_record(...)` / `march_windows_resume(record, ...,
  first_changed=k)`: snapshots of the small per-node states every 512 steps;
  resumes from the last one a lookahead (tube + nose nodes) before k; a window
  whose render could not have changed (k past its end plus a nose length) is
  taken from the record, **render_window is not called for it**.
- `marching.BelowPeel(..., record=True)` and `.resumed(port_pressure, from_front,
  returned, k, geometry, port_y, jet)` for the peel below the port across
  hypotheses; `_peel_below(..., until=s)` stops a chunk past what the flow
  before s needs (entries after s are then not the whole peel's).

Gate (MEASURED, `bench_kernels.py resume`): each word's time and node variants
evaluated whole and resumed from the truth's record; geometry (both ways),
peel / observer outputs / all march outputs and the peel below compared:
**0 failures over 10 cases**. One documented exception: `march_windows`'s
`waves` output on rows *between* windows (a render's free run, never read by the
march, not used by march_residual) keeps a reused window's render: 106 rows
differ in b-time and d-time; all rows the march reads are identical.

| case | first change (left) | kernels whole -> resumed (CPU s) | engine calls |
|---|---|---|---|
| c time | 15 264 (61%) | 0.140 -> 0.082 | - |
| c node | 21 216 (46%) | 0.154 -> 0.071 | - |
| a time | 21 216 (46%) | 0.086 -> 0.063 | free_flow 1 -> 1 |
| a node | 21 600 (45%) | 0.090 -> 0.044 | free_flow 1 -> 0 |
| e time | 9 696 (81%) | 0.177 -> 0.172 | 1 -> 1 |
| e node | 34 656 (33%) | 0.177 -> 0.079 | 1 -> 0 |
| b time | 21 216 (46%) | 0.436 -> 0.252 | window renders 2 -> 1 |
| b node | 21 600 (45%) | 0.420 -> 0.218 | 2 -> 1 |
| d time | 34 656 (67%) | 1.54 -> 1.08 | 4 -> 3 |
| d node | 86 496 (17%) | 1.59 -> 0.56 | 4 -> 2 |

Geometry through the workspace: 0.006-0.08 s against 0.03-0.13 whole (the
copying resume was 0.5-0.7 s twice under memory pressure).

## Recipes for the residual owners (not applied: sequence_recover is theirs)

- march_residual's `render_window`: pass `until=max(begin - port_node - 2, 0)` to
  `marching._peel_below`; only that prefix is read (forced is NaN after it).
- Finite-difference columns: per worker, evaluate the base hypothesis once with
  the `*_record` APIs (and `moving_geometry_record`), then for each column compute
  `first_changed` as above and call the `*_resume` APIs with a
  `GeometryWorkspace` geometry. Order columns by decreasing first change to keep
  the workspace's copies to nothing. The final engine render of each residual
  still runs whole until the engine itself can resume.

## Caveats

- `moving_geometry` (and its record/resume) now returns **read-only** arrays.
  No caller in analysis/ writes into a geometry (grep); one that does will get
  numpy's error.
- The jet tuple from `marching.jet_table` has two more entries (noise table,
  coefficient cache); only jet_table builds it.
- `observer._forward_steps` takes two more arguments (noise, cache; None keeps
  the reference's jet); `marching._march_run` takes the seal counts and the
  linings' column offset.
- Disclosure: before caches were redirected, the first bench runs on the
  reference tree let numba add compiled specializations to
  `~/ducklips/tract/analysis/__pycache__` (index files of `_peel_moving_side`,
  `_peel_states`, `_forward_steps`, `_march_run`, 10:55-11:24). Source was not
  touched; numba's cache entries are per signature and do not change results.
  Left in place (the reference's own batch may have written some of them).
