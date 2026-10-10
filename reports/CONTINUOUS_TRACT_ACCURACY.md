# Separate work balance from motion accuracy

The coupled continuous source approaches an independent open-passage trajectory
under step refinement. This establishes accuracy over the tested interval
before the first closure; it does not establish accuracy through all later
closures, a human voice, or a speech rendering. Performance work remains parked.

## Independent open motion

`continuous_tract_open_reference.py` integrates the continuous-field source and
the same eight-cell lossless tube with SciPy Radau. It starts from the complete
8 ms checkpoint used in the contact controls and stops 425 microseconds later,
before the first seal. The smallest area is about 0.00166243 mm2. No positive
area floor is applied. The source's local gas effort and the graph's linearized
pressure remain the work-conjugate interface used in the earlier coupling test.

Two independent quadrature/tolerance settings use 96 points with `rtol=1e-8`
and 192 points with `rtol=2e-9`. Their largest normalized difference over 1,361
common samples is 5.20e-8. These are fresh ODE integrations. The
[receipt](continuous_tract_open_reference.json) records their inputs, hashes,
normalization, output prefixes and comparisons.

The resistance-fitted native update gives the following errors on **the same
18 sample times**, with no interpolation:

| Maximum native step | Maximum normalized state error |
|---|---:|
| 1.25 microseconds | 0.0529954 |
| 0.625 microseconds | 0.0130076 |
| 0.3125 microseconds | 0.00321070 |

Successive error ratios are 4.074 and 4.051. The largest error in each of these
comparisons is the upper gas effort, scaled by 1,000 Pa. On their respective
output grids, maximum tract pressure errors are 0.002524, 0.0006311 and
0.0001585 Pa. These small tract errors coexist with larger internal pressure
errors. The source is not interchangeable with a slowly varying flow source
just because its external tube pressure looks settled.

The broader 16 ms contact controls are a separate test. Prefix snapshots are
immutable and checksummed even when the longer job is still running. A passed
prefix is not a completed long trajectory.

## Audit the force derivative before changing it

`continuous_force_probe.cpp` evaluates the native fixed-mass, fixed-momentum
kinetic derivative in the symmetric tissue subspace. The independent Python
check differentiates the field inertia with complex steps, then transforms the
result from effort/area coordinates to mass/area coordinates. It does not use
the native finite difference as its reference.

At the four saved open states from 8.000 to 8.425 ms, the production difference
rule's largest errors are 2.43e-6 Pa for the gas mass derivatives and 4.79e-9
J/m2 for the area derivatives. Quarter and sixteenth difference steps are also
recorded. This finding does not support changing the derivative step to repair
the much larger long-trajectory error. It is not a validation of closed states,
contact derivatives, or asymmetric motion. See the
[force receipt](continuous_force_reference.json).

## Request accuracy explicitly

The tube research driver now accepts `--tolerance`. For each proposed interval,
it compares one full evolution with two half evolutions, each recomputing the
implicit source/tract load. If the physical state difference exceeds the
requested normalized tolerance, it subdivides and tries again. The accepted
state is the two-half-step result, with their original mass and work receipts;
there is no extrapolation of state or deletion of stored gas.

The norm uses absolute scales: 1,000 Pa gas effort, 1 mm2 aperture, 1e-4 m3/s
flow, 0.001 m2/s symmetric aperture rate, 1e-8 m3 graph charge, and 0.1 Pa s
branch momentum. Those scales are for this configured symmetric test body.
The raw full/half difference is used, without dividing by three across contact.
A local tolerance is not a bound on the accumulated 16 ms error. The default
fixed-step driver and production source remain unchanged.

This is an accuracy reference with copied trial states, not a realtime engine
or an efficiency result. Its purpose is to distinguish a successful nonlinear
solve from an accurate step before choosing the source's next time representation.

Fresh smooth checks cover the first 100 microseconds from the common checkpoint,
with a requested maximum step of 5 microseconds and 21 common output times:

| Requested local tolerance | Maximum sampled normalized error | Accepted substeps |
|---|---:|---:|
| 1e-4 | 3.744e-5 | 332 |
| 2.5e-5 | 1.131e-5 | 512 |
| 6.25e-6 | 4.229e-6 | 874 |

The error decreases when tolerance is tightened. This is a short smooth test;
it does not prove the same error bound over the full compression pulse or
through repeated contact. The controller owns subdivision: failed single-step
trials return directly to it, avoiding a recursive solve inside another
recursive solve. Earlier prototype jobs with nested subdivision were stopped
and are not counted as completed validations.

A fresh 100 microsecond closure control, starting from the saved pre-seal native
state, also completes under ASAN/UBSAN with one seal. It accepts 324 substeps,
reaching a minimum step of 9.765625 ns. Its maximum accumulated work residual
is 4.85e-20 J and mass residual is 1.19e-22 m3; its separately recorded numerical
loss is 4.26e-13 J. These are measured solver/ledger properties, not physical
calibration. The [adaptive receipt](continuous_tract_adaptive_reference.json)
records the commands, independent comparison, binaries and outputs.

## Completed longer refinements

The frozen eight-cell tube now has four completed resistance-fitted runs from
the same full checkpoint, each covering another 16 ms and recording three
closures and three reopenings. Maximum open steps are 2.5, 1.25, 0.625 and
0.3125 microseconds; crossing caps are 78.125, 78.125, 39.0625 and 19.53125 ns.
Thus the last two refinements reduce both open-step and crossing-step error.

| Adjacent maximum steps, microseconds | Maximum sampled tract-pressure difference | Final tract-pressure difference | Final aperture difference |
|---|---:|---:|---:|
| 2.5 / 1.25 | 658.15 Pa | 206.70 Pa | 5.0551 mm2 |
| 1.25 / 0.625 | 136.90 Pa | 34.199 Pa | 1.0787 mm2 |
| 0.625 / 0.3125 | 29.077 Pa | 8.4682 Pa | 0.30692 mm2 |

At 6 ms, the final tract-pressure differences across those pairs are much
smaller: 1.595, 0.3891 and 0.09293 Pa. The larger differences accumulate later,
through subsequent motion. Each comparison uses aligned recorded times, with
no interpolation. Their sample counts differ, so endpoint differences are the
cleaner common-time refinement comparison. These are adjacent native-run
differences, not errors against an independent contact oracle.

The finest run's maximum accumulated work residual is 2.61e-15 J, mass residual
is 1.52e-20 m3, and separately recorded numerical loss is 5.11e-10 J. The
[completed refinement receipt](continuous_tract_refinement.json) records the
snapshot binary, input, outputs, commands and early/late comparisons.

This is evidence that the coupled contact motion is approaching a common
trajectory. The remaining aperture difference and the lack of an independent
through-contact oracle prevent treating the finest result as exact. Short
open-trajectory agreement, repeated-contact refinement and human listening
remain distinct forms of evidence. None of these tests is a new NATO recording.
