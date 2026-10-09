# Native continuous modal transport

The previous direct moving-tree reference converged accurately, but its useful
substepped configurations ran below 10x realtime. This implementation separates
rapid acoustic oscillation from the much slower physical geometry transport.
It remains a reference component: production speech PCM and its remaining
geometry-clock lines are unchanged.

## Local interaction coordinates

The independently verified moving-mode equations have the form

```
x_dot = (J Omega(t) + M(t)) x
x = [a, b], H = dot(x,x)/2
M = block_diag(-F, Omega_dot/Omega + Omega F^T Omega^-1)
```

For each output sample, set `x(t)=R(theta(t)) y(t)`, where `R` contains one
sine/cosine rotation per mode and `theta` integrates its angular frequency.
Then `y_dot=R^-1 M R y`. RK4 now integrates this geometry coupling; it does
not approximate the rapid diagonal oscillator with a polynomial step.

Frequency samples at the beginning, midpoint and end define a quadratic.
Its half- and whole-sample phase integrals are

```
theta_half = dt * (5 omega_start + 8 omega_middle - omega_end) / 24
theta_full = dt * (omega_start + 4 omega_middle + omega_end) / 6
```

This is exact for a quadratic frequency trajectory, not for every possible
motion. The four RK stages rotate into physical modal coordinates, apply the
two real transport matrices, and rotate the derivative back. The final state
is rotated by the whole-sample phase. A stationary oscillator therefore has
an exact acoustic phase step rather than the previous Gauss frequency bias.

Native API: [continuous_modes.hpp](../include/ducklips/continuous_modes.hpp).
After preparation, advance uses fixed storage. Tests verify one second of
stationary 12 kHz phase, the independently known integral of a quadratic
chirp, the exact uniform-area-motion reference, and dense spatial rotations
with an analytic solution at 8, 9, 20 and 32 modes. The dense tests exercise
SIMD blocks and odd tails and count heap allocations during advance: zero.

## Deriving the matrices from physical rates

[ModalTransportReference](../include/ducklips/modal_transport.hpp) accepts the
existing `AcousticNetwork`, its C-orthonormal `ModalBasis`, and physical
`dC/dt` and `dG/dt` arrays (`G` is reciprocal inertance). It computes

```
V_ij = sum_nodes phi_i * C_dot * phi_j
D_ij = sum_branches delta_phi_i * G_dot * delta_phi_j
       - (omega_i^2 + omega_j^2) V_ij/2
Q_ij = D_ij / (omega_j^2 - omega_i^2), i != j; Q_ii = 0
F = V/2 + Q
omega_dot_i = D_ii / (2 omega_i)
```

The public output is the two transport matrices and angular frequency/rate.
No geometry-specific empirical gain is fitted. The builder rejects nonfinite
input/results, nonpositive frequencies, wrong dimensions, and eigenvalue gaps
at or below 1e-8 of the largest squared frequency. Rejection retains the last
valid output. The caller must supply a physically valid basis with continuous
signs/order; this builder does not solve the eigenproblem or handle crossings.

The first implementation repeatedly multiplied parameter rates inside each
pairwise sum and accessed modal scalars through an out-of-line accessor.
At 32 modes its measured coefficient construction cost was 8.66x realtime for
1,000 constructions per simulated second, below the floor. The replacement
stores contiguous pressure/branch-difference vectors with rates absorbed once
per feature. Pairwise projection becomes two dot products. Both construction
and application share explicit independent accumulators, ARM64 SIMD where
available, and a scalar fallback; fast-math is not enabled.

## Independent accuracy checks

The native builder processes every half-sample of the five-node nonuniform
fixture. Those native-derived coefficients feed the native acoustic step.
The full resulting physical pressure/flow state is compared with independent
adaptive integration of the canonical network, not just the modal equations.

At one step per 48 kHz output sample, maximum relative state errors over
40 ms are approximately:

| Initial state | Native continuous transport | Prior direct Gauss step |
|---|---:|---:|
| Lowest nominal mode | 1.34e-8 | 9.11e-6 |
| Highest nominal mode | 4.42e-7 | 0.0630 |
| Deterministic mixed modal state | 2.75e-7 | Not measured in the prior receipt |

Native physical transport coefficients match the independent matrix-derivative
construction within 7.47e-14 using a unit-floor relative metric. Native and
Python interaction trajectories agree within about 3.4e-15 absolute state
norm. These are isolated, positive, lossless modes with gradient momentum;
they do not validate arbitrary branch forcing, viscothermal memory, closures,
mode crossings, or a complete human airway.

See [continuous_modes_reference.json](continuous_modes_reference.json) for
source/binary hashes, complete errors and timing receipts. The Python
subdivision experiment is [interaction_modal_reference.json](interaction_modal_reference.json).

## Timing scope and the remaining coefficient-sampling problem

Kernel timings include three coefficient fills per audio sample, phase
calculation, four transport products and the state update. They exclude
physical transport construction, which is timed separately on 96 nodes and
191 branches at 1,000 constructions per simulated second. That construction
benchmark uses representative dense basis/rate arrays to measure cost; its
input is not a fitted speech trajectory. Eigenbasis solving is excluded from
both measurements. Receipts retain three runs at each of 20, 24 and 32 modes.

The final receipt measures approximately 28x, 23x and 14.5x for the
kernel at those respective counts, and 110x, 78x and 47x for construction at
1 kHz. Summing the separate worst costs gives a budget estimate of 22.1x,
17.6x and 11.0x. This sum is **not** an end-to-end speech benchmark or evidence
that a 1 kHz coefficient stream has sufficient accuracy.

The accurate native fixture above supplies physical coefficients at every
half audio sample. A separate sampling experiment explicitly tests the gap
between that accuracy check and the cheaper control-rate timing:

| Local interpolation | Control rate | Highest-mode relative state error |
|---|---:|---:|
| Cubic, four knots | 1 kHz | 2.51e-4 |
| Cubic, four knots | 2 kHz | 1.60e-5 |
| Quintic, six knots | 1 kHz | 4.57e-6 |
| Quintic, six knots | 2 kHz | 4.30e-7 |

The six-knot experiment requires three future knots, hence 1.5 ms of lookahead
at 2 kHz. It is not an online predictor, and no causal coefficient scheduler
or native polynomial interpolator has been implemented here. Doubling
construction cadence also doubles its cost; at 32 modes, the summed kernel
plus 2 kHz construction estimate already falls below 10x before interpolation,
eigenbasis, source and microphone work. No whole-renderer speed acceptance is
claimed. The smaller current speech modal counts have more budget, but need
those same omitted costs measured.

Sampling receipt: [interpolated_motion_reference.json](interpolated_motion_reference.json).
The useful change is an accurate, efficient native acoustic transport and a
native physical coefficient builder. The next integration work must provide
continuous, causally scheduled coefficients and preserve the 10x floor with
actual articulation, source loading, losses and radiation.

## Reproduce

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/interaction_modal_reference.py --out out/geometry_clock/interaction.json
python tools/verify_continuous_modes.py --out out/geometry_clock/continuous-modes-native-transport
python tools/interpolate_motion_reference.py --out out/geometry_clock/interpolated-motion.json
```

The enhanced spectrogram, original human reference and production speech
rendering path are not modified by this reference implementation.

Validation on this revision: core Release 6/6, core ASan/UBSan 6/6, and surface
Release 8/8 CTest suites pass. The independent native physical-rate and acoustic
reference verifier passes its stated component gates. A fresh production Alfa
render is byte-identical to `out/geometry_clock/alfa.wav`, has zero nonfinite or
clipped samples, and renders 1.31 seconds in 0.079 seconds in this single-run
regression check. This does not measure the new components inside speech.
