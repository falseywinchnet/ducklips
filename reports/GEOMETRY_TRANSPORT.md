# Geometry-clock diagnosis and a moving-energy reference

The remaining horizontal lines continue to follow the geometry/control clock
after the physiological-force fix. In the 0.48–0.58 s part of synthetic Alfa,
changing the interval from 48 to 64 samples moves the strongest high bands from
near 4/7 kHz to near 5.25/4.5/3.75 kHz. These are multiples of the respective
1 kHz and 750 Hz clocks. This is measured from PCM; no second spectrogram or
changed enhanced display scale is involved.

## What the geometry event changes

The optional sample trace now separates the instantaneous lip-radiation
readout jump into two parts:

1. The expected result of the selected physical metric invariant and the new
   lip-branch conductance, using the old physical lip pressure.
2. The remaining error after projecting that pressure onto the new modal basis.

It also records the glottal pressure jump, lip-readout slope jump, acoustic
energy ratio, and largest relative compliance/conductance change. The trace
calculation does not alter the generated PCM (bit-identical to the previous
continuous-physiology render).

In the problematic 0.48–0.58 s window the direct lip metric contribution is
zero: the lips themselves have not changed there. The lip projection jumps
are small (largest about 9.04e-6 m³/s²), while their free-evolution slope jumps
reach 119.3 m³/s³. Interior geometry changes still act on the stored field.
This is not explained by a large direct switch of the lip radiation coefficient.
The diagnostic reports describe boundary jumps, not the complete subsequent
forced response.

Changing only the state-transfer invariant supplies a stronger causal test:

| Geometry rule | 4 kHz component, 0.48–0.58 s | Interpretation |
|---|---:|---|
| Default `passive` | about -57.2 dBFS | Chooses the lower-energy invariant at each metric change |
| `charge_momentum` | about -57.2 dBFS | Conserved charge/momentum applied as abrupt geometry steps |
| `pressure_flow` | about -99.1 dBFS | Holds pressure/flow through the geometry boundary |

The pressure/flow case is approximately 42 dB cleaner at this line, but its
largest measured geometry event increases acoustic energy by 1.014%. Such an
increase may represent wall work or an inconsistent mapping; this experiment
alone does not distinguish them. The production default was not switched based
on its appearance. A temporary correction of the radiation readout's damping
term did not materially affect the line and was removed.

Measured traces and spectra are in [geometry_clock_events.json](geometry_clock_events.json)
and [geometry_clock_lines.json](geometry_clock_lines.json). `--trace` now adds
six named geometry fields to the existing CSV. Geometry work is evaluated only
when tracing is enabled; ordinary rendering does not run the diagnostic sums.

## Exact reference: uniform area motion

A tract that changes shape nonuniformly involves changing modal shapes and
coupling between modes. First isolate the geometric energy terms in a case
with no changing eigenvectors: a lossless tube whose every area is multiplied
by the same positive factor `s(t)`, with lengths fixed. Compliance and reciprocal
inertance scale together:

```
C(t) = s(t) C0       G(t) = s(t) G0       M(t) = M0/s(t)
q = C p             lambda = M u
H = q²/(2C) + G lambda²/2
```

For the explicitly assumed canonical lumped-network equations
`q_dot = -B u`, `lambda_dot = B^T p`, the source-free power balance is

```
H_dot = -p² C_dot/2 + lambda² G_dot/2.
```

This is a reference for the acoustic network's metric transport. It does not
include mean wall-displacement volume flow, convective terms, losses, closures,
or all mechanics of a real moving airway. Those additions require their own
power accounting; the reference is not empirical human validation.

For each normalized mode set `a=q/sqrt(C)`, `b=sqrt(G)*lambda`, so that
`H=(a²+b²)/2`. For exponential area motion `s=exp(r*t)`, the frozen frequency
`omega` stays fixed and the exact homogeneous equation is

```
a_dot = -r*a/2 - omega*b
b_dot =  omega*a + r*b/2
H_dot = r*(b²-a²)/2.
```

The generator `A` satisfies `A²=-(omega²-r²/4) I`. Thus its propagator requires
only a 2×2 real map, with a closed-form sine/cosine expression (or its
hyperbolic counterpart). With constant modal damping `alpha`, multiply it by
`exp(-alpha*dt)`. In complex notation the new term couples `z` to its conjugate;
one ordinary complex multiplier cannot represent that motion.

The native implementation is [moving_mode.cpp](../src/moving_mode.cpp), with
its explicit restricted contract in [moving_mode.hpp](../include/ducklips/moving_mode.hpp).
Preparation calculates four coefficients; applying a prepared step needs four
real multiplications and two additions. It is a reference building block,
**not integrated as the general speech renderer**.

## Why finer passive stepping is insufficient

The reference uses a 1,234 Hz free mode and a fourfold area change over 0.1 s,
with initial normalized state `(0.3, 0.7)`. The exact solution and integrated
wall work agree to better than 2.3e-16 energy units. An independently computed
matrix exponential verifies the closed form.

Under uniform expansion, the old `passive` mapping scales both quadratures by
`1/sqrt(s_new/s_old)`. Under contraction it scales both by
`sqrt(s_new/s_old)`. Consequently, a fourfold expansion **or** contraction
leaves only a quarter of the starting energy, independent of subdivision.
The exact moving reference does not have that 75% artificial loss. Its RMS
state discrepancy from the stepped passive trajectory is roughly 31%, even
at one geometry update per audio sample. This proves that raising the update
rate alone cannot repair that model of geometric energy transfer.

Charge/momentum stepping instead converges toward the specified continuous
canonical model: its RMS state error falls from about 0.49% at 1 kHz geometry
updates to 0.010% at 48 kHz updates. Those finite-step errors still justify a
continuous treatment. Pressure/flow transport has the opposite metric signs
and does not converge to this particular canonical reference.

The report [moving_metric_reference.json](moving_metric_reference.json) retains
all three rules at six update intervals, including a stationary control. The
new tests check sample-subdivision invariance, the canonical determinant,
independently integrated wall work, the existing stationary modal limit, and
transactional rejection of invalid parameters.

## Reproduce and remaining work

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/moving_metric_reference.py --out out/geometry_clock/moving-reference.json
build/ducklips_say --dictionary samples/nato.dict --text 'Alfa.' \
  --trace 0 out/geometry_clock/events.csv --out out/geometry_clock/alfa.wav
build/ducklips_say --dictionary samples/nato.dict --text 'Alfa.' \
  --invariant pressure --trace 0 out/geometry_clock/pressure-events.csv \
  --out out/geometry_clock/pressure-invariant.wav
python tools/analyze_control_lines.py out/geometry_clock/alfa.wav \
  out/geometry_clock/pressure-invariant.wav --out out/geometry_clock/lines.json
```

The next generalization must include nonuniform compliance/inertance rates,
changing modal shapes, source/radiation coupling and articulator work, and
must be checked against this uniform limit. The current production render
retains the previously fixed continuous physiology and still has the remaining
geometry-clock defect. No claim of completed removal or perfect speech is made.

Validation: Release and ASan/UBSan each pass all three core CTest suites; the
surface build passes its five registered suites. A fresh seven-case, one-run-
per-case compact benchmark passes the 10× floor (see
[geometry_clock_performance.json](geometry_clock_performance.json)); this is
separate from the prior 21-run receipt. The new reference propagator is not
being counted as a complete speech-rendering benchmark.

The subsequent nonuniform work is documented in
[NONUNIFORM_MOTION.md](NONUNIFORM_MOTION.md): canonical-coordinate compression,
a direct moving-tree reference whose accurate substeps fail the 10x floor,
and a verified continuous modal-transport oracle for the next representation.
