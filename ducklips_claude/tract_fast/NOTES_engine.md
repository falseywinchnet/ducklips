> **OWNER RULE (2026-10-20): NO CALIBRATION, EVER.** No fine-tuning, fitting, per-word patches, per-speaker profile fits, or tuned constants. For the past week fit was substituted for understanding. Progress comes only from structure: the physics, the direct inverse, and knowing why something works or breaks. Search scaffolding added since 2026-10-19 (arrival grids, time sweeps, joint guards, sane clamps, block polishes, vowel revisits) is labelled REJECTED as a method: it is fit, not understanding. Round-trip dB numbers from it show only that the model can be inverted; they are not progress.

# The engine, in-process and faster (tract_fast, 2026-10-20)

Scope: tract_fast/src, include, tools, CMakeLists.txt, bindings/, and
analysis/engine.py. The frozen reference (~/ducklips/tract) was only read and
run, never edited. All measurements were taken on the M4 Mac mini, which
other jobs shared (load 5 to 31).

Labels: MEASURED (reproduced here, with the named script), REJECTED (tried,
failed its gate), OPEN (not done or not verified).

## What was built

- `include/tract/sustain_run.hpp`, `src/sustain_run.cpp`: the sustain tool's
  run as a reusable value. `SustainSettings` holds every input the tool takes
  from its arguments, its environment and the files those name, as numbers.
  `settings_from_environment` is the one parser of the environment's text.
  `SustainRun` builds the body and steps it; `StreamFiles` writes the
  per-step files. The step loop is the tool's loop moved over unchanged. The
  file is compiled with the tool's old flags (target `tract_run`: no
  `-mcpu=native`, no `-ffp-contract=fast`), so its keyframe interpolations
  round as before.
- `tools/sustain.cpp`: now a thin command line over that run. The command
  line and the engine cannot drift.
- `bindings/tract_engine.h`, `bindings/tract_engine.cpp`: a plain C ABI that
  builds `libtract_engine.dylib`. It takes settings named by the tool's
  environment variables (numbers, pairs, text, and arrays named by the files
  they replace), renders into caller buffers (RAW pair, PORTWAVES row, GEODUMP
  rows), and offers runs that can be snapshotted (`tract_run_copy`) and
  resumed with changed inputs, plus `tract_settings_divergence`.
- `analysis/engine.py` (ctypes, numpy only):
  - `render(...)`: posture_recover.render's exact signature and return
    values. It writes `port_waves`, `nose_dump` and `wave` exactly as the
    reference does, and resolves relative paths against ROOT as the
    reference's cwd does.
  - `prepare`, `Settings`, `Run` (`advance`, `snapshot`, `advance_into`,
    `nose`).
  - `Incremental`: render() that re-runs only from the first step whose
    inputs changed since its last call.
- Harness and benchmarks in `bindings/`: `equivalence.py` (the gate),
  `bench.py` (reference against engine per render), `bench_engine.py`
  (builds interleaved, in-process), `bench_instructions.py` (instructions and
  cycles by /usr/bin/time -l), `bench_incremental.py` (a search's pattern).
- Build: `cmake -S . -B build_engine -DCMAKE_BUILD_TYPE=Release && cmake
  --build build_engine -j 2`. Every target builds without warnings.

### Exactness details of render()

- The reference passes the pitch and the seconds with `%g` (6 significant
  digits: 117.123456789 Hz reaches the tool as 117.123). engine.render applies
  the same rounding, `float("%g" % x)`. MEASURED identical (gate case 2).
- Keyframe files are read by the tool as one stream of numbers taken 4 or 5 at
  a time (FOLDKEYS/FOLDSEQ) or 2 at a time (VELUMSEQ/JETSEQ). engine.py groups
  them the same way. AREASEQ is read line by line, as before.
- The reference's subprocess inherits the caller's environment. Any sustain
  variable render does not set itself is passed through and parsed by the
  tool's own parser. TRACE is dropped: it only prints, and the reference
  discards its output. MEASURED identical with TREMOR and PITCHEND set.
- Non-finite values written as text (NaN in a keyframe) stop the tool's
  reader. engine.py refuses them with ValueError rather than imitating that.
  Binary inputs (forced waves, observed flows) take NaN as the reference does.

## Equivalence (gate: bindings/equivalence.py, byte-for-byte via tobytes)

MEASURED, final build. In-process engine and this tree's CLI against the
frozen reference posture_recover.render driving ~/ducklips/tract/build_dev/
sustain. Every case is identical in pairs and written files:

| case | steps |
|---|---|
| plain /aa/ 0.3 s | 28800 |
| /iy/, pitch 117.123456789, 0.4123456 s (%g rounding) | 39585 |
| /aa b aa/ closure with seal, 0.41 s | 39360 |
| nasal /aa m aa/: velum keys, PORTWAVES file, NOSEDUMP file, 0.54 s | 51840 |
| jet /aa s iy/ with gate keys and closing (seal), 0.41 s | 39360 |
| jet without gate or seal | 28800 |
| voiceless /aa p aa/, FOLDKEYS | 39360 |
| FOLDSEQ triggers | 28800 |
| forced waves with NaN stretches | 28800 |
| observer flows on forced folds | 28800 |
| AREAFILE2 glide, with the WAV file | 28800 |
| everything: nasal, jet, folds, seal, 0.54 s | 51840 |
| caller's environment passed through (TREMOR, PITCHEND) | 28800 |

- Snapshot gate, MEASURED identical: render(0..N) equals render(0..k) +
  resume(k..N) for k = 1, 12, 96, 5000, 19213 and 30000. Both the copy and
  the original resume identically. With forced waves and the last keyframe
  changed from k = 18000 on, the resumed run equals a whole render of the
  changed inputs.
- Incremental gate, MEASURED identical against the reference. Twelve renders
  of changed inputs: the same again, a later frame, a keyframe 1.3 ms later, a
  longer render, a jet key moved, forced waves added and cut, back to the
  first, a nasal word, a velum key, a new initial tube. It reused 0 to 50880
  steps; a different body or a step-0 dump reuses nothing.
- The command line with its defaults, MEASURED identical against the
  reference binary (RAW, PORTDUMP, WAV bytes). This covers all side branches,
  aspiration, frication and wall variation, with and without VELUM, TREMOR,
  RELEASE, FRICATION and PITCHEND.

## Speed

All timings are MEASURED at load about 5. The reference is
posture_recover.render with the reference binary; the engine is
engine.render. Medians of 9 interleaved repeats (bindings/bench.py). Every
engine result was checked byte for byte against the reference.

| word | reference wall / cpu ms | engine wall / cpu ms | speedup |
|---|---|---|---|
| 0.41 s stop (seal) | 89.7 / 89.3 | 61.8 / 61.7 | 1.45x |
| 0.41 s nasal (nose) | 106.0 / 105.2 | 75.3 / 75.1 | 1.41x |
| 0.41 s fricative (jet) | 108.4 / 107.5 | 77.5 / 77.2 | 1.40x |
| 0.54 s stop (seal) | 125.8 / 118.7 | 84.5 / 84.3 | 1.49x |
| 0.54 s nasal (nose) | 138.5 / 134.8 | 96.4 / 95.8 | 1.44x |
| 0.54 s fricative (jet) | 144.2 / 143.4 | 97.8 / 97.5 | 1.47x |

### Where the speed comes from

- **Process and file overhead removed.** The subprocess, keyframe files, RAW
  file and temporary folder cost only 10 to 15 ms per render, not most of
  the 0.1 s (engine before any step change against reference: 1.15 to 1.2x).
  The stepping is the cost.
- **Per step.** bindings/bench_engine.py, minimum over interleaved rounds
  against the unoptimized in-process build: 0.75 to 0.82 of its time. Stop
  word 0.41 s: 67 to 52 ms; 0.54 s jet word: 108 to 81 ms.
- **Incremental, a search's pattern** (bindings/bench_incremental.py; 0.54 s
  nasal word; one later keyframe's time and one area moved per render):

  | keyframe moved | steps reused | reference ms | engine ms | incremental ms |
  |---|---|---|---|---|
  | the last | 31% | 132.9 | 95.8 | 54.7 (2.4x the reference) |
  | the third | 19% | not run | 96.7 | 73.2 |
  | the second, third and last in turn | 22% | 140.7 | 89.2 | 62.3 (2.3x the reference) |

  The later the phone a search moves, the more is reused. A search over the
  last phone of a long word, or over forced waves past a point, re-runs only
  the tail.

### The changes, each byte-identical under the gate

1. **Jet gain memo** (voice_body: reversible_jet, shaped_jet_noise).
   - What: a filter's corner gains are taken again only when its corner
     changes. The corner sits clamped at 1500 Hz or the ceiling most of the
     time. The filter state evolves exactly as before, so a later gate
     opening is unaffected.
   - MEASURED: 9 to 11% fewer instructions on jet words.
2. **Tremor normalization once per process.** A constant that was measured
   with 105,000 iterations at every body construction (about 0.3 ms per
   render).
3. **Viscous-law memo per link** (pow only when a link's area changes).
4. **Pinned tube arithmetic** (src/tube_step.cpp, compiled with
   `-ffp-contract=off`): see the REJECTED item below, which made it
   necessary.
   - The lining step, its factorization, the air junctions and the series
     links now carry each multiply-add the reference build fused as an
     explicit `std::fma`, read from the reference object code.
   - Their roundings are now fixed by the source, not by the compiler's
     fusion choices.
5. **Geometry-only quantities kept while their inputs hold** (compared with
   `==`; same inputs give the same values):
   - per link: A/rho c, rho c/A, and the kinetic coefficients in both
     directions;
   - per node: the junction's lining resistance and its step's factors;
   - the radiation load;
   - the lining's matrix, elimination and body denominators, now refactored
     only for recompiled nodes, propagating while the elimination factor
     changes.
6. **Side-port chains** built once a step, replacing a scan of every port at
   every node, twice. MEASURED: step_air plus scatter went from 10060 to 8178
   instructions per vocal step.
7. **Diagonal lining solve.** With no axial coupling (AXIAL=0's vocal tube,
   REVNOSE's nose) every elimination factor is a zero, and each node's
   division stands alone. All divisions are taken at once, then checked: a
   nonzero right side and a finite neighbour make fma(+-0, c, r) equal r. If
   any check fails, the reference elimination runs.
   - MEASURED: the vocal tube's lining went from about 1000 to 390 cycles
     per step (it ran at 2.1 instructions per cycle, latency-bound); the
     vocal tube's whole step from about 2900 to 2230 cycles.
8. **Links recompiled only beside a moved node.** A link's coupling depends
   on its two nodes' radii only.

## REJECTED

- **Caching inside the reference's `-ffp-contract=fast` code.**
  - What happened: keeping one quantity of step_air in a cache unpaired the
    reference's two-lane vector divide (drive and resistance together). The
    compiler then fused the drive's add into an FMA. Every render moved by
    up to 6e-8 absolute (about 1e-9 relative) from step 62.
  - `#pragma clang fp contract(off)` does not stop the backend under
    `-ffp-contract=fast`. MEASURED: the fused instruction is still emitted.
  - Fix: the separate file with `-ffp-contract=off` and explicit FMAs
    (item 4).
- **Timing under load.** Wall and CPU times swing by 2x with the core a
  process lands on (4 performance and 6 efficiency cores, load 15 to 31).
  Instructions retired are steady but blind to division latency. The
  comparisons above were therefore taken at load about 5, interleaved, and
  minimums were used for build-to-build ratios. Under `taskpolicy -b`
  (efficiency cores only), CPU time repeated within 1%, but runs starved
  under load.
- **Skipping the zero-gain noise sources** (sound_jet with FRICATION=0,
  aspiration with ASPIRATION=0) was not done. A non-finite input would
  inject NaN at a different step than the reference, so a diverging render
  could stop at a different step. The gain would have been about 2 to 3%.

## OPEN

- The equivalence gate covers render()'s settings and the CLI's defaults.
  The other tools (speak, voices, the sheet tools) use paths the gate does
  not reach and were not compared:
  - the articulatory body's opening slot, branch glide, area limit and
    conserved span;
  - hop loss by area, and a nose without REVNOSE (its full tridiagonal
    solve).
  The pinned step_air carries those paths as transcribed from the reference
  object code. Gate them before using tract_fast's speak.
- tract_fast/build_dev still points at build_ref (the reference binary).
  analysis/engine.py loads build_engine/libtract_engine.dylib. Pointing
  posture_recover at the engine is for the analysis owners:
  `posture_recover.render = engine.render`, or one `engine.Incremental()`
  per worker process.
- Further per-step work not done:
  - refresh_geometry is still in the fast-contract file; it could be
    pinned and given a static-tube path;
  - the trachea's tridiagonal chain (33 nodes) remains serial;
  - step_air could peel its boundary nodes.
- Snapshot size was not measured. Incremental keeps one every 960 steps
  (10 ms): about 50 per 0.5 s word.
