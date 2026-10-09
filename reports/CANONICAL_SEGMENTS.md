# Compact canonical transport with native control segments

This extends [CONTINUOUS_MODES.md](CONTINUOUS_MODES.md) from coefficients
supplied every half audio sample to physical control knots at 1 kHz. The
physical-rate builder, coefficient interpolator and acoustic step are all
native. The speech renderer has not been switched to this path; its remaining
geometry-clock lines are not claimed removed by this reference experiment.

## One transport matrix determines the other

Writing pressure transport as `P=-F`, the verified canonical equations give

```
a_dot = -Omega b + P a
b_dot =  Omega a - Omega P^T Omega^-1 b + diag(omega_dot/omega) b
```

Thus one dense matrix, frequencies and a diagonal connection replace two
dense matrices. The compact update forms `h=b/omega`, applies `P` to pressure
coordinates and `-P^T` to `h`, then restores the frequency scaling. The native
row traversal reuses each loaded matrix row for both products. It preserves
the existing interaction-coordinate integration; it is not a sparse threshold
or removal of small physical couplings.

`ModalMotionSample` retains the full-matrix form as a reference. An empty
flow matrix selects the canonical form with an explicit diagonal. The
physical builder exposes it through `canonical_view()`. The compact form
requires strictly positive frequencies and the canonical dual relation;
arbitrary independent damping/forcing matrices cannot be silently substituted.
Dense analytic rotation tests cover both representations, SIMD blocks and
odd tails, and confirm zero heap allocation during acoustic stepping.

## Reuse frequency derivatives; join pressure slopes

The physical builder already computes `omega_dot`. A quintic Hermite
frequency polynomial can use the values and rates at three knots, `-1,0,1`,
instead of requiring six independent frequency values. Its own derivative,
divided by its frequency, supplies the diagonal connection. Pressure and flow
are no longer interpolated as unrelated matrices.

An initial four-knot cubic pressure interpolant passed the highest-mode test,
but failed the mixed-state gate at 1.47e-6 maximum relative state error. The
1e-6 gate was retained. Its rejected receipt is
[canonical_segment_cubic_rejected.json](canonical_segment_cubic_rejected.json).

The final pressure segment is cubic Hermite. Shared knot slopes use the
fourth-order centered difference, in units of one control interval:

```
slope_i = (P_(i-2) - 8 P_(i-1) + 8 P_(i+1) - P_(i+2)) / 12
```

Both adjacent segments use the same value and slope at a knot. Pressure
transport is therefore C1 across control intervals. Frequency values and
rates also agree at shared knots. Tests check polynomial reproduction,
pressure value/slope continuity, and transactional rejection of invalid
intervals, times, and nonpositive evaluated frequencies.

[CanonicalMotionSegment](../include/ducklips/motion_segment.hpp) consumes six
uniformly spaced knots `-2,-1,0,1,2,3` and evaluates the central interval.
Only a cubic matrix polynomial is evaluated at audio rate; the quintic and
its derivative operate on one value per mode. At 1 kHz the segment requires
three future knots, or **3 ms lookahead**. The diagnostic supplies them.
A command-to-audio scheduler implementing that delay is not yet present.
A delayed audio path could make the lookahead causal, but input commands,
physiology controls, geometry and source feedback must be aligned and tested.

## Runtime representation and validation

The first native interpolation-plus-step prototype measured about 8.7x at
32 modes. It repeated finite-value checks and copied a full matrix at every
half sample. During segment configuration, an absolute coefficient-sum bound
now certifies the entire pressure polynomial on `[0,1]`, with overflow margin.
After frequency validation, matrix evaluation can write directly to the
output frame. This removes repeated work without removing the finite checks
at the point where new control data enter the segment.

The two acoustic phase angles also share information:

```
theta_full = 2 theta_half + dt*(omega_end-omega_start)/4
```

For correction magnitude at most 1e-3 radians, the implementation reuses the
half rotation and evaluates the tiny correction with sine/cosine Taylor
polynomials. The first omitted cosine term is bounded by `1e-18/720`, below
double precision rounding. Larger corrections use direct trigonometric
functions. Tests cover both sides of the threshold and the fallback against
the independently integrated phase.

## Native accuracy at the measured control cadence

The verifier starts from the five-node physical network and its independent
frozen eigensystem at each 1 kHz knot. The native builder derives frequency
rates and transport from `C_dot`, `G_dot` and pressure shapes. These native
outputs feed the native segment evaluator and compact acoustic step at
48 kHz. Physical pressure/flow states are then compared against independent
adaptive integration of the full canonical network over 40 ms.

| Initial condition | Maximum relative physical state error |
|---|---:|
| Lowest nominal mode | 3.50e-8 |
| Highest nominal mode | 4.91e-7 |
| Deterministic mixed state | 3.47e-7 |

All three pass the retained 1e-6 gate. These are positive, isolated modes in a
lossless network with gradient branch momentum. This fixture does not exercise
mode crossings, closures, viscothermal memory, or empirical human realism.
Its highest instantaneous mode is about 4 kHz; dynamic accuracy across the
full speech bandwidth still needs a broader reference.
Series-pressure sources can create nongradient momentum, so integrating
frication requires an additional physical treatment; simply wiring a noise
value into this restricted reference is not sufficient.

The receipt is [canonical_segment_reference.json](canonical_segment_reference.json).
It includes source/binary hashes, all errors and timing runs. The preliminary
frequency-derivative study, before the mixed-state pressure correction, is
[canonical_motion_reference.json](canonical_motion_reference.json).

## Timing scope

The interpolation-plus-step benchmark runs five simulated seconds using a
bounded ring of six control knots. Timing includes new knot filling, segment
configuration, two half-sample evaluations per audio sample, phase calculation,
and acoustic evolution. It asserts zero heap allocations throughout that
loop. Physical-rate construction is measured separately at 1 kHz on 96 nodes
and 191 branches. The verifier retains three runs each at 20, 24 and 32 modes
and also reports the sum of the separate worst costs as a budget estimate.

The final measured interpolation-plus-step ranges are 27.78–28.22x at
20 modes, 21.18–21.52x at 24, and 13.68–13.93x at 32. All nine runs record
zero heap allocations. Adding the separate worst coefficient-construction
costs gives budget estimates of 22.22x, 16.69x and 10.65x, respectively.

Neither benchmark includes geometry construction, eigenbasis solving,
continuously orienting/tracking modes, a real command scheduler, the vocal
fold source, full physical losses, or microphone radiation. The summed cost
is not a measured complete speech-rendering rate. A narrow margin at 32 modes
cannot be used to claim the full engine meets 10x. The 10x requirement remains
an acceptance gate for that integration.

## Reproduce and continue

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/canonical_motion_reference.py --out out/geometry_clock/canonical-motion.json
python tools/verify_motion_segment.py --out out/geometry_clock/canonical-segment-c1
```

The next work is the real geometry/eigenbasis producer and correctly delayed
source coupling, including nongradient forcing and crossing/closure handling.
Production PCM and the preserved enhanced spectrogram are unchanged here.

Validation: core Release 7/7, core ASan/UBSan 7/7, and surface Release 9/9
CTest suites pass. The native physical-knot/segment/acoustic verifier passes
its retained accuracy and component-cost gates.
