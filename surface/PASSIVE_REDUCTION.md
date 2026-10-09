# Passive frozen surface reduction

This experiment retains the detailed surface model's cavities, wall filters,
fractional delays and transverse/floor paths in a normalized-energy state
realization, then projects those states before changing to modal coordinates.
It is a frozen, zero-mean-flow acoustic prototype, not the completed moving
speech engine. The existing speech renderers remain in place.

## What changed in the representation

The full-body snapshot allocates 319,232 delay slots, but its exact dynamic
realization needs 20,220 energy states. Each integer delay stores only its
actual length. Allpass, interpolation, loss-shelf, radiation and mechanical
wall states use energy coordinates. A component's scattering step is
contractive; lossless junction scattering conserves power. Their connected
network therefore cannot create energy when its geometry is fixed.

External acoustic wave ports obey

```
x_next = A x + B u
v      = C x + D u
||x_next||² + ||v||² <= ||x||² + ||u||²
```

Here `u` and `v` are power-normalized incident and reflected waves, not PCM.
Prescribed volume flow closes these wave ports algebraically and recovers the
original pressure/flow network. The unreduced dynamic formulation agrees with
the independent path-elimination operator within 7.69e-11 relative complex
error across all five drives and six observations at eight sampled frequencies
from 48 Hz to 12 kHz. The largest component transfer error is 1.24e-13.
Random energy checks are supplementary diagnostics, not a proof by themselves;
the component contractions and junction power identity establish the energy
construction.

An orthonormal basis `V` gives the reduced scattering step

```
Ar = Vᵀ A V, Br = Vᵀ B, Cr = C V, Dr = D
```

Orthogonal projection cannot increase stored energy, so this operation
preserves contractivity. Volume-flow closure comes **after** projection.
Complex eigenmodes are merely coordinates for this reduced system: their
physical energy metric is transformed with them and is not an unweighted sum
of squared modal amplitudes. No poles are clamped and no wall damping is added.

## Why anatomical blocks matter

A single global response basis introduced enormous relative leakage through
the nearly closed nose, even when the strong oral responses looked accurate.
Tightening numerical rank tolerance did not cure it. Partitioning the basis
by anatomical cavity, with separate groups for connecting channels, preserves
the missing block connections. This prevents the projection from inventing a
strong direct mouth-to-nose connection. Snapshot metadata now identifies the
anatomical tube of terminal pieces as well as ordinary pieces.

A 48-anchor block basis had 473 states and looked good under glottal drive
(maximum magnitude discrepancy about 0.162 dB on the reported level mask).
Testing all five drives exposed a lateral response error of 1.97 dB. A
96-anchor basis increased the order to only 501 states and resolved that
poorly sampled resonance. Increasing order blindly was not the diagnosis;
the independent drive revealed which response was missing from the evidence.

The final 2,401-frequency sweep covers 50–12,000 Hz, all five inputs, and all
six outputs. It has:

- Maximum normalized RMS complex error: 1.04e-5.
- Maximum relative complex error with a 1%-of-each-transfer-peak floor: 3.49e-5.
- Maximum magnitude discrepancy above the same level threshold: 0.000147 dB.
- Maximum phase discrepancy above that threshold: 1.88e-5 radians.
- Reduced scattering-step singular value: 1 + 3.3e-15 (floating-point roundoff).

These are sampled frequency-response measurements, not a continuous-frequency
error bound or perceptual validation. Each transfer uses its own peak: pressure
and radiation flow are never scaled against each other. A 241-point sweep has
a different sampled peak and hence a different level mask; its maximum reported
magnitude error is 0.000393 dB. Both receipts are preserved.

## Native update and the speed failure that remains

The 501 real states comprise 200 complex conjugate pairs and 101 real modes
(301 stored mode records). The C++ runtime stores only the necessary real
coordinates. Real-coordinate readouts use an ARM64 vector dot product, with a
portable scalar fallback; source/coupling equations do not use fast-math.
Preparation allocates; sample updates operate on prepared vectors.

Native impulses for every input agree with the dense reduced equations within
2.49e-6 relative complex error using each transfer's 1% peak floor. This test
uses 0.0625-second impulses and the same 400/s exponential contour on both
sides, so the finite tail is explicitly suppressed. The modal export also
checks undamped responses at 241 frequencies. Weak nasal transfers expose
several parts per million of eigencoordinate roundoff; the report retains raw
errors as well as the floored metric.

Five consecutive 10-second probes on this M4, at 192 kHz:

| Probe | Measured speed | Gate |
|---|---:|---|
| Coupled eight-strip two-mass folds, glottal acoustic drive and radiation | 19.68–19.97× | Pass for this frozen case |
| All five source ports active, all six outputs, including Gaussian source generation | 8.73–8.80× | **Fail** |

The simultaneous-source case is a numerical throughput exercise, not a human
speech sample. It produces no audio file. All reference recordings remain
human speech. Timings exclude model preparation and file I/O. The earlier
full surface impulse timing of about 0.205× includes gather/scatter clock
instrumentation and uses a different workload, so it is not an exact speedup
ratio against the coupled-fold benchmark.

The general path initially needed scalar complex accumulations. Real states,
contiguous source updates and vector readouts improved its representation.
Additional blocked source passes and paired readouts did not improve measured
throughput and were removed. Isolated generation of 9.6 million Gaussian
values took about 0.097 s; it does not explain away the remaining failure.
`verify_surface_modes.py` deliberately exits nonzero on that speed result.
The new representation is **not accepted as the speech renderer**.

The next reduction question is how many of these energy coordinates are both
reachable from acoustic source ports and observable at the necessary pressure
and radiation ports. Any further reduction must retain anatomical isolation
and passivity. Dropping low-amplitude modal residues without carrying their
energy metric would not meet that requirement.

## Reproduce

Python requires NumPy and SciPy. The current local runtime is
`/Users/joshuahkuttenkuler/miniforge3/bin/python`.

```sh
cmake -S surface -B build/surface -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/surface
ctest --test-dir build/surface --output-on-failure
python tools/test_passive_surface.py
mkdir -p out/port_response
build/surface/surface_response --fixture body --duration 0.25 --out out/port_response/body
python tools/passive_surface.py out/port_response/body.json --out out/port_response/body-energy.json
python tools/reduce_surface.py out/port_response/body.json --anchors 96 --anatomical \
  --validation-points 2401 --out out/port_response/body-reduced-fine
python tools/export_surface_modes.py out/port_response/body-reduced-fine.npz \
  --out out/port_response/body-fine.modes
python tools/verify_surface_modes.py out/port_response/body-reduced-fine.npz \
  out/port_response/body-fine.modes --out out/port_response/body-fine-native.json
```

The last command currently returns failure because the all-port speed is below
10×. Full evidence, binary/model/source hashes, rejected candidates and timing
runs are in [`reports/passive_surface_reduction.json`](../reports/passive_surface_reduction.json).
The Git head there precedes the uncommitted-at-measurement work; hashes identify
the actual binary, model and mode file. Source hashes are taken at receipt
assembly, with a subsequent input-metadata guard identified in the receipt.
Generated arrays and mode files stay outside Git.

Release and ASan/UBSan builds each pass all four CTest suites. Three Python
component tests pass. The modal tests exercise real/conjugate recurrence,
source coupling, reset, invalid ports and transactional rejection of bad files.

## Unresolved scope

This does not yet cover other bodies and articulatory poses, moving geometry,
mean-flow kinetic resistance, source-location changes, or membrane edge
coupling. Offline basis preparation is much too slow for a per-nudge solve.
A moving implementation must account for articulator work and state transport
between energy coordinates, rather than crossfading unrelated frozen voices.
Human realism, flexible tongue-blade/root behavior and body matching still
require the experiments in [`RESEARCH.md`](../RESEARCH.md). Spectrogram
inspection uses the full enhanced field alone in its original colors/range.
