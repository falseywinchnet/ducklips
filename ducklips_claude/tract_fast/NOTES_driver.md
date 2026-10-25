> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. For the past week fit was substituted for understanding. Progress comes only from structure: the physics, the direct inverse, and knowing why something works or breaks. Search scaffolding added since 2026-10-19 (arrival grids, time sweeps, joint guards, sane clamps, block polishes, vowel revisits) is labelled REJECTED as a method: it is fit, not understanding. Round-trip dB numbers from it show only that the model can be inverted; they are not progress.

# The word driver, made faster (tract_fast, 2026-10-20)

Scope: `analysis/sequence_recover.py`, `tube_recover.py`, `posture_recover.py`
(scratch, render switch, pitch-ear switch) and the new `driver_ledger.py`
(instrumentation), `column_locality.py` (the locality interface) and
`driver_bench.py` (equivalence checks and costs). The frozen reference
(`~/ducklips/tract`, `/tmp/es`) was read, never written. Labels: MEASURED
(reproduced here; run and file named), REJECTED, OPEN. All runs on the Mac
mini under the sibling projects' load (load average 15-32 on 10 cores), so
wall times are noisy; evaluation counts are exact and CPU per evaluation is
the steadier cost.

## 0. Switches and housekeeping

| env | default | effect |
|---|---|---|
| `POOL` | 10 (the reference's) | workers |
| `CHUNK` | 1 | tasks handed out one at a time (pool.map's default chunks made a 145-column Jacobian on 4 workers 15 chunks of 10: 40 tasks on the busiest worker against 37) |
| `FAST_EXACT` | 1 | every exact shortcut of section 2 (0: the reference's path, for A/B) |
| `RESUME` | 1 | finite-difference columns resumed from the base point's record (2.6) |
| `FIRST_VOWEL_CACHE` | 1 | the first vowel's search kept on disk (2.7) |
| `FAST_ALGO` | none | comma list of path-changing options: `halving`, `jointstop`, `stall1`, `reuse` (section 3) |
| `RENDERER` | `sustain` | `sustain`: the tool as a subprocess (the reference's); `engine`: `engine.render` in-process; `incremental`: `engine.Incremental` per process and render kind (2.8) |
| `PITCH_EAR` | `0.0` | PITCHFB for every render (model change, section 5); `engine`: the engine's own (the reference's model before today) |
| `LEDGER_DEPTH`, `LEDGER_TAG` | 4, "" | the time table's depth; the ledger file's suffix |

- Scratch: `/tmp/esf` (`posture_recover.SCRATCH`; the nose constants and
  port waves in `sequence_recover` too).
- Ported the coordinator's fix: `fold_keyframes` / `jet_keyframes` only use
  phones that exist in the arrivals they are given (prefix joints).
- `sequence_recover.main` is now `build_word` (the known answer),
  `configure_targets`, `recover_word` (the search) and the reports; the
  search itself is line for line the reference's, with stages named for the
  ledger.

## 1. Where the time goes (MEASURED)

Instrument: `driver_ledger` (a stage tree; evaluations counted where they
are computed, CPU of the pool's workers measured in the workers including the
engine subprocesses they wait for, main-process CPU by os.times; "saved" =
evaluations an exact shortcut made unnecessary, i.e. the reference would have
computed them). The ledger is printed at the end of every run and kept as
`~/ducklips/tract_fast/out/reversible/ledger_<word><tag>.json`.

### /aa b aa/, exact shortcuts on, sustain subprocess, POOL=4

Run: `bench/exact_AA_B_lips_AA.log` (Mini, `~/ducklips/tract_fast/bench/`).
Every number it printed (survey, every posture and node rung, transitions,
prefix and final joints, round trip 134.4 dB) is identical to the
reference's log `/tmp/es/word_AA_B_lips_AA.log`; only the joints'
evaluation counts differ (1970/2552/1873/1970 against 3054/3918/2910/3054).

Computed 50 983 evaluations, saved 12 883 (the reference path: 63 866).
Wall 105.5 min, CPU 285.8 min (0.34 CPU s per evaluation on average; workers
on the old kernels, renders by subprocess).

| stage | evaluations | saved | wall min | CPU min | share of wall |
|---|---|---|---|---|---|
| first /aa/ (no context): survey 1024, 6 posture basins, 1 node search (7 rungs, -245.5 dB) | 7 376 | 969 | 10.9 | 17.4 | 10% |
| transition 1 (into /b/): 4 grids, 290 | 290 | 2 | 0.5 | 1.5 | 0.5% |
| transition 2 (release scans) | 103 | 3 | 0.2 | 0.5 | 0.2% |
| second /aa/ recover: survey | 1 024 | 0 | 0.7 | 2.1 | 0.7% |
| second /aa/ recover: 6 posture basins | 2 903 | 563 | 9.5 | 10.2 | 9% |
| second /aa/ recover: 6 node searches | 18 539 | 1 890 | 48.8 | 128.3 | **46%** |
| second /aa/ polish x2 (transition 126 + node search) | 4 697 | 627 | 7.7 | 26.1 | 7% |
| prefix joint (4 rounds) | 7 686 | 4 250 | 12.1 | 44.2 | 11% |
| final joint (4 rounds) | 8 365 | 4 579 | 15.1 | 55.6 | 14% |
| **word** | **50 983** | **12 883** | **105.5** | **285.8** | |

By kind (exclusive): node-search rungs 26 532 evaluations, 60.8 wall min
(58%); final joint 15.1; posture stages 14.3; prefix joint 12.1; surveys
1.9; every timing grid together 0.9.

Reading: the cost is the node searches (every basin of every later vowel,
five to seven Levenberg-Marquardt runs of 30 over 47 parameters each), then
the joints. The timing grids the task suspected are under 1%. The node
search that won (-59.8 dB) was the third of six; the first alone took 23.6
wall min and stopped at -59.0.

### /aa m aa/, exact shortcuts on, the reference's model, POOL=4

Run: `bench/exact_AA_M_lips_0_47_AA.log`, run with `ENGINE_GLIDE_PHASE=0
PITCH_EAR=engine` (the reference's model) on the kernels' owner's new
kernels. Its first vowel is line for line the reference's; it then departs
(transition 1 leaving 0.0998 against the reference's 0.0999, velum grid
0.444 against 0.474) and ends at a round trip of 54.8 dB (reference 63.1).
OPEN, not the driver's (MEASURED): every driver change on the march path was
checked bitwise on this model (bench `truncate` and `offsets`: 5 targets at
the hint and 30 points of the transition-1 grid, all bitwise), and the
/aa b aa/ run, made before the kernels were swapped, matched the reference
to the last digit. Suspects: the new kernels with ENGINE_GLIDE_PHASE=0 not
reproducing the old replica bit for bit on the nasal path, or the engine
build the fast tree links (`build_dev -> build_ref`, relinked 10:55)
differing from the one that made the reference's log at 09:53.

Computed 27 006 evaluations, saved 5 957; wall 43.7 min, CPU 136.7 min
(0.30 CPU s per evaluation; load had fallen to ~5).

| stage | evaluations | saved | wall min | CPU min |
|---|---|---|---|---|
| first /aa/ (as in /aa b aa/) | 7 376 | 969 | 3.9 | 11.4 |
| transition 1 + velum grid (58) | 348 | 2 | 0.5 | 1.9 |
| second /aa/ recover: survey, 6 postures, 6 node searches | 9 448 | 2 235 | 18.4 | 52.4 |
| second /aa/ polish x2 | 7 768 | 787 | 16.4 | 56.1 |
| prefix joint | 1 843 | 1 407 | 3.9 | 13.2 |
| final joint: round 0 did not move, rounds 1-3 skipped as exact repeats | 120 | 554 | 0.4 | 0.9 |
| **word** | **27 006** | **5 957** | **43.7** | **136.7** |

The reference's four unmoving final rounds cost 4 x 314 = 1 256 evaluations;
here 120 (the discarded report Jacobians and the three repeats removed).
Node-search rungs are again the bulk (17 139 evaluations, 32 of 44 min).

## 2. Exact changes (the same numbers)

Each is checked bitwise by `driver_bench.py` (MEASURED, outputs quoted) and,
for the whole search, by the line-for-line identical /aa b aa/ log above.

2.1 **leastsq instead of least_squares (MINPACK lmder either way).**
scipy 1.16's `least_squares(method='lm')` evaluates f and J at x0 before
MINPACK asks for them again (cached by its VectorFunction) and, after
MINPACK returns, computes `J = jac(x)` at the result for its report: a whole
Jacobian per stage that the driver discarded. `leastsq` calls the same
`_minpack._lmder` with the same arguments (factor 100, diag None, the same
tolerances and maxfev). With a memo of the last residuals and the last
Jacobian by the parameters' bytes, its own shape checks cost nothing.
Bench `stage` (a posture stage, 9 parameters, 15 nfev): 75 -> 60
evaluations, parameters bitwise equal. The stalled joints of /aa m aa/ and
Echo spent 314 evaluations a round: 145 at x0, ~24 trial points, 145 for
the discarded report.

2.2 **Memoized residuals.** MINPACK evaluates f(x) then asks J(x): the
Jacobian's base point is the point just evaluated (one evaluation per
Jacobian), and least_squares' f(x0) twice.

2.3 **Inert columns.** A stop's (and a word's leading consonant's) own log
areas are parameters that `SequenceShape.split` replaces by its vowel's: 46
of the 143 columns of /aa b aa/'s joints are the base residual (bend 0) or
the base with only the bend rows moved. Built without a render
(`SequenceShape.inert_columns` / `inert_residual`). Bench `inert`: 144 -> 98
evaluations per Jacobian, Jacobian bitwise equal at bend 0 and bend 1.

2.4 **Repeated stages skipped.** A deterministic stage started from the
same bytes on the same target repeats itself: `StillRounds` skips a joint
round after one that did not move (the reference's /aa m aa/ and Echo:
rounds 1-3 of the final joint, 3 x 314 evaluations); the node search keeps
a memo of (start bytes, bend) -> result per search (a rung that did not
improve leaves the next unbent rung on the same start: /aa b aa/'s polish
rung 4, 77 evaluations). With a single node count (a word's) the posture
search's levels evaluation that decides nothing is skipped.

2.5 **Renders and peels only as far as the compared stretch.** The engine
is causal and its output does not depend on the run's length (a 0.41 s and
a 0.2 s render agreed bitwise over the 0.2 s), and the peels read a bounded
reach ahead (a 46-node tube's peel was bitwise the same on the sound cut 50
steps past the stretch, not at 25): every residual whose target is not
low-passed (filtfilt reads both ways) renders to its stretch's end and peels
4 node counts past it (`sequence_recover.horizon`, `reach`). A closure's
free render stops at the last closure's hidden end (observer reads its folds
only over [hidden_from, hidden_to)); a march window's render at its end
plus the nose's lag (the windowed step reads `share` at t + m - 1); the
peel below the port for a window render with `until=` (NOTES_kernels.md's
recipe). Bench `truncate` on /aa b aa/, /aa m aa/, Echo, five target kinds
each (whole, timing grid exact and low-passed, prefix, release scan): 15 of
15 bitwise. Single-call wall on /aa b aa/: grid exact 0.66 -> 0.24 s, whole
0.93 -> 0.58 s, release scan 0.22 -> 0.15 s. For a word of N phones the
early phones' searches no longer render the whole word (they did: every
evaluation of phone k rendered all N phones).

2.6 **Columns resumed from the base's record (the locality of section 4,
now live for sealed words).** `Problem.jacobian` hands each worker one
batch of columns (`column_batch`): the base point once with
`reverse_core.moving_geometry_record` and `observer.implied_record`, then
its columns in decreasing order of their first input divergence with a
`GeometryWorkspace` and `observer.implied_resume` (closures wholly before
the change keep their free render). Columns are dealt round robin by resume
step so every batch holds early and late ones. Bench `resume` on /aa b
aa/'s joint Jacobian (98 computed columns, POOL=4): whole 29.6 CPU s,
resumed 23.6 (sustain) / 17.1 (incremental engine); Jacobian bitwise equal
at bend 0 and bend 1. Found on the way (MEASURED, reported to the kernels'
owner): `GeometryWorkspace.resume` from boundary 0 starts from the base's
own start state (made from its first rest row); a first phone's column
changes that row, and 0.0016-0.024 came out off. The driver computes such
columns whole (`first_changed < record.every`).

2.7 **The first vowel's search kept.** With no context (a word starting with
a vowel and no jet) the first vowel's search depends only on its hint, its
window, the node count, the sound to the window's end plus the reach, the
search's constants, the model switches and the engine's bytes: /aa b aa/
and /aa m aa/ ran it identically line for line in the reference. Kept as
`out/cache/first_vowel_<sha>.npz` under that key. Bench `first`: the
equation error and the spectrum on the sound cut at the key's reach are
bitwise the whole sound's. Saves 10.9 wall min (7 376 evaluations) for every
further word with that opening.

2.8 **The in-process engine behind one switch** (`posture_recover.RENDERER`;
the render's signature unchanged). Bench `renderer` on /aa b aa/ and /aa m
aa/, five target kinds: sustain, engine and incremental bitwise equal; CPU s
per residual (medians of 3) whole 0.24 / 0.19 / 0.10, grid exact
0.13 / 0.10 / 0.06, prefix 0.21 / 0.16 / 0.09 (incremental's figure is its
best case, the same hypothesis again). The incremental engine is kept per
process and per render kind (a closure's free render and the driven render
alternate in one residual).

2.9 Smaller: timing grids memoize each rung's centre (the rung before's
best, 1 evaluation a rung); pool tasks one at a time.

## 3. Algorithmic changes (a different path; judged by the final result)

Options (`FAST_ALGO`), each changing the search's path:

- `halving`: every carried basin gets its node search's first rung (bend
  1.0), then only the deepest `HALVING_KEEP` = 2 go on (their first rung
  reused from the memo). Evidence for it: in the reference's /aa m aa/
  second /aa/ the winning basin (-79.0 dB) was second by its posture value
  but first after one rung; on the corrected model (section 5) the first
  rungs of all six basins came out equal (/aa b aa/ first vowel: -88.6 dB x
  6; second /aa/: -65.9 x 4, -37.2, -35.4; /aa m aa/ second /aa/: -72.2 x 6),
  i.e. the basins merge after one bent rung and the rest of their searches
  were near-duplicates.
- `stall1`: a node search ends after one unbent rung without 1 dB of gain
  (reference: two).
- `jointstop`: an unbent joint round gaining less than 0.5 dB ends its
  rounds (on its own residual, no extra evaluation).
- `jointmore`: the final joint goes on past 4 rounds while a round gains at
  least 1 dB, to 8 (MEASURED below: the 4-round result had not converged).
- `reuse`: a later vowel equal to an earlier one starts a node search from
  the earlier tube and skips the survey when it reaches -60 dB. Not used in
  the validation: in these known-answer words the repeated vowel's tube is
  identical by construction, so its benefit would be overstated here
  (OPEN: to be judged on words whose repeats differ, or on real speech).

Weak point (OPEN): with `halving` a first vowel whose first basin would
have reached EXACT (-150 dB) and ended the search still pays every basin's
first rung and the second carried basin's start (/aa b aa/ first vowel:
10 158 evaluations against the reference's 7 376). Run the posture-best
basin whole first, then halve the rest, when it does not reach EXACT.

## 4. Locality: the interface (goal 3)

`column_locality.py` holds the driver's side; section 2.6 uses it.

- `SequenceShape.input_tracks(parameters, target)`: the keyframe tracks the
  engine and the peels read (areas, velum, folds, jet).
- `first_divergence(base_tracks, moved_tracks)`: the first time any track
  can differ (linear tracks: a change at keyframe i bends the track from
  keyframe i - 1, or from the earlier of the two times when only the end of
  a held stretch moved).
- `resume_step(divergence, nodes)`: less the peels' reach
  (`LOOKAHEAD_NODES` = 4 node counts; the kernels' resume APIs take
  `first_changed` and handle their own reach, so the driver passes the first
  differing schedule row and track step instead).
- `ResumableResidual` protocol: `run_with_snapshots(parameters, target,
  steps) -> (residual, snapshots)`, `resume(parameters, target, snapshot)`,
  bitwise a whole computation; `FullRecompute` the stub (today's whole
  columns). The sealed residual implements it concretely as
  `sealed_core(..., keep=True)` / `sealed_core(..., base=record)`.
- A detail any resumed residual must keep: `fixed_clip` removes the
  stretch's mean, so a column's change moves every row; a resume carries the
  base's raw difference before its step and clips after.
- Validation (bench `locality`, /aa b aa/, 9 columns): the predicted resume
  step was at or before the first row that departs from the base's constant
  shift in every column (e.g. the third phone's nodes: resume 0.2148 s, first
  departing row 0.2169 s; the last leaving: no row departs).
- Projected share of column steps that need computing with resume
  (MEASURED over every Jacobian of the /aa b aa/ run): 47.0% (23 309
  columns). Today's realised gain is smaller (2.6: 20-42% CPU) because the
  driven engine render runs whole and first-phone columns are computed whole.

OPEN: the march path (nasal words) has `marching.march_windows_record` /
`march_windows_resume`; not yet wired into `march_residual` (same pattern
as `sealed_core`). The driven render of each residual is still whole on the
subprocess and engine renderers; the incremental renderer resumes it by
itself.

## 5. Model changes (not efficiency; kept apart)

- `PITCH_EAR = 0.0` (PITCHFB=0 in every render, truth and hypotheses; env
  `PITCH_EAR=engine` restores the engine's own): the coordinator's finding
  that the pitch ear's per-closure, whole-step tension correction is the
  residual's discontinuity behind the stalled joints.
- `reverse_core.ENGINE_GLIDE_PHASE` (the kernels' owner; default 1): the
  links' viscous glide in the engine's phase.
- Both are part of the first vowel's cache key.

## 6. Validation

All on the corrected model (section 5: PITCH_EAR 0, glide in phase),
`RENDERER=incremental`, exact shortcuts and resume on, `FAST_ALGO=halving,
jointstop,stall1` (Echo also `jointmore`), POOL=4, load ~5. Logs:
`~/ducklips/tract_fast/bench/fast_*.log`, `continue_*.log`; ledgers
`out/reversible/ledger_*_fast.json`. The truths' errors on this model:
/aa b aa/ -242.3, /aa m aa/ -160.6, Echo -245.6 dB (the coordinator's
figures, reproduced).

| word | reference round trip | fast run (4 joint rounds) | + jointmore | tubes rms max | times max | velum | wall (POOL 4) | CPU | evaluations computed (saved) | CPU s / eval |
|---|---|---|---|---|---|---|---|---|---|---|
| /aa b aa/ | 134.4 dB | 130.0 dB | **140.9 dB** (8 rounds; 141.2 after a 9th that gained 0.29) | 0.0077 -> 0.0029 cm^2 | 0.000 ms | - | 24.9 min (+8.5 for rounds 5-8, run beside another word) | 86.3 min (+30.4) | 50 739 (14 238) (+10 220) | 0.10 (0.18 beside another run) |
| /aa m aa/ | 63.1 dB | **131.2 dB** | 131.3 dB (a 5th round gained 0.40 and would end it) | 0.0037 | 0.000 ms | 0.4700 (true 0.470) | 50.4 min | 179.7 min | 39 330 (14 110) | 0.27 |
| Echo /eh k! ow/ | 66.4 dB | **0.8 dB: FAILED** | (jointstop ended the joint after 1 round) | 0.661 | 0.70 ms | - | 16.0 min | 53.5 min | 34 978 (8 050) | 0.09 |

Against the reference runs (POOL 10, 44-73 min a word, 0.61 / 1.15 CPU s per
residual for stop / nasal words): /aa b aa/ 24.9 + 8.5 min on 4 workers;
CPU per evaluation 0.10 s (/aa b aa/; 0.34 in the exact run on the old
kernels with subprocess renders), 0.27 s (/aa m aa/; the march path is not
yet resumed). The exact run of /aa b aa/ against this one: CPU 285.8 ->
86.3 min (3.3x) for about the same number of evaluations; the reference
path computed 63 866 for it.

Readings (MEASURED):

- /aa b aa/ with four joint rounds came out 4.4 dB under the reference's
  round trip (130.0 against 134.4) on the new model and path: its final
  joint had not converged (2.1-3.0 dB a round). `jointmore` (on to 8 rounds
  while a round gains 1 dB) brings it to 140.9 dB and the third tube from
  0.0077 to 0.0029 cm^2 (joint_continue.py from the 4-round result: the
  start rounded through the saved .npz, -150.9 dB and 130.0 dB re-measured
  identical). The validated configuration is therefore
  `halving,jointstop,stall1,jointmore`.
- /aa m aa/: the final joint now moves (2 676, 2 676, 2 480, 2 088
  evaluations; -147.2 -> -159.2 dB) where the reference's four rounds did
  not; round trip 63.1 -> 131.2 dB, velum exact to 4 digits. That gain is
  the model change's (the pitch ear's discontinuity gone), not the driver's.
- Echo FAILED, at its first transition, before any algorithmic change acts
  (the first vowel was exact, -248.4 dB; the timing grids are untouched by
  `FAST_ALGO`): the /eh/ -> /k!/ grids chose leaving 0.0943 s (true 0.1000),
  cost 5.9e4 against ~1e-17 at the true times, and nothing after recovered
  (the prefix joint went to -25.9 dB, /eh/'s tube 0.219 off). Isolated with
  `transition_probe.py` (the transition alone from the true first tube):
  PITCH_EAR=engine with the old glide: 0.1017 / 0.1499, the reference's own
  (1.7 ms early: the reference recovered from it); PITCH_EAR=0 with either
  glide: 0.0943 / 0.1500. Along the leaving with the arrival exactly true
  the cost is one clean V (807 at -0.5 ms, 1.7e3 at +0.5;
  `bench/scan.py`), so the 2-D grids' coarse rungs pick a wrong diagonal
  pair. REJECTED (MEASURED): a 1-D leaving scan (+-8 ms, 0.5 ms) with the
  arrival held, then the fine rungs again: 0.0962, cost 3.0e4 (/aa b aa/
  and /aa m aa/ unchanged at 0.1000 / 0.0999). OPEN for the model's owner:
  the voiceless stop's closing time on the pitch-ear-free model.

## 7. OPEN

- The march path's columns (nasal words) are not resumed
  (`marching.march_windows_record` / `_resume` exist; wire them as
  `sealed_core` is). /aa m aa/ spends 0.27 CPU s an evaluation against 0.10
  for /aa b aa/.
- /aa m aa/ on the reference's model in the fast tree departs from the
  reference's log at transition 1 (section 1); the driver's march-path
  changes are bitwise on that model; the kernels or the engine build are
  the suspects.
- Echo's first transition on the corrected model (section 6).
- An exact-only run on the corrected model, to split the validation's
  differences between the model and `FAST_ALGO`.
- `halving` on a first vowel whose best basin would reach EXACT (section 3).
- `reuse` judged on words whose repeated vowels differ.
- Node searches remain the largest cost (58% of the exact /aa b aa/ run's
  wall). With the pitch ear off the basins merge after one bent rung;
  deduplicating equal rung results (bitwise equal -> exact) or ending a
  basin whose first rung lands within 0.1 dB of another's are the next
  measurable cuts.
- The engine's driven render per residual: incremental resumes it within a
  worker by itself; with the subprocess or plain engine it runs whole.
- `GeometryWorkspace.resume` from boundary 0 (section 2.6), reported to the
  kernels' owner.

## 8. Files

- `analysis/sequence_recover.py`, `analysis/tube_recover.py`,
  `analysis/posture_recover.py` (owned; changed as above).
- `analysis/driver_ledger.py`, `analysis/column_locality.py`,
  `analysis/driver_bench.py` (tests: cost, truncate, offsets, first,
  renderer, stage, inert, locality, resume), `analysis/joint_continue.py`,
  `analysis/transition_probe.py` (new).
- `bench/*.log` here: the runs' logs (warnings filtered).
- Mini: `~/ducklips/tract_fast/bench/` (logs, chain scripts, scan.py,
  smoke scripts), `out/reversible/ledger_*.json`, `out/cache/`.
