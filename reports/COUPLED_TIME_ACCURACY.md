# Coupled source and airway time accuracy

The passive three-stage update reduced the lip-mode damping error, but did not
resolve the source's fast pressure evolution. The preceding fixed-source
comparison still had a maximum sampled internal pressure error of 3.994 Pa,
and the finest moving closing refinements differed by 9.97 Pa near the source.
Those errors concern the trajectory, even when its work balance closes.

`JointAccuracy` now applies step doubling to the complete continuous source,
pressure/flow graph, thermal memory, viscous memory and radiation state. It is
an accuracy reference for advancing the developing engine. It has not rendered
an utterance or established a speech improvement. The last listening evidence
still belongs to the spatial-pressure recording: November improved, Foxtrot's
/tr/ remained jammed.

## Accepted motion and geometry

Each trial compares one full step with two half steps. Only the two smaller
steps are committed; there is no extrapolation of the physical state. The
accepted work, physical heat, numerical loss, reservoir transfer and internal
port work are the sums of those steps' actual receipts. Pressure and flow
means are duration-weighted over the entire caller interval. Multiplying
these two means does not recover the port work, which is recorded separately.

The geometry provider receives the actual start and end times of each trial.
It supplies old, midpoint and endpoint geometry, moving-volume flow, external
drives and the source outlet inertia. It can be called for a rejected trial,
so it must not commit body state. The moving-mouth adapter recomputes the
actual contour and profile geometry at those times; it does not interpolate
the old endpoint acoustic metrics to invent intermediate openings.

The error is the largest normalized physical-state difference. Its scales
are 1000 Pa for source and tract pressure, 1e-6 m² for source opening,
1e-8 m³ for source gas moments, 1e-4 m for fold displacement, 1e-4 m³/s for
source flow, and 0.0357142857 m/s for individual fold velocities. Graph
momentum and relaxation memories use their endpoint energy metric with a
1e-6 J reference. These are explicit research scales, not a perceptual model
or a claim that a local tolerance bounds accumulated error.

Prepared scratch graphs and bounded work arrays avoid allocation during
advancement. Topology, memory ownership and passive stage count must match
the prepared layout. Failure leaves the caller's physical state, accepted
receipt, means and accepted-time trace unchanged. Attempt diagnostics identify
configuration rejection, exhausted trial budget or an unresolved minimum step.

Subdivision preserves the requested endpoints using a stack of pending
interval ends. Repeated floating-point addition previously left a tiny final
interval that the minimum-step check could not resolve. This correction
preserves elapsed time; it does not silently drop that remainder.

## Small-step source arithmetic

The independent fixed-geometry fixture begins with nonzero graph pressure and
momentum but zero thermal and viscous memories. This deliberately excites a
very fast diffusion startup. At tighter trajectory tolerance, the source
Newton solve failed after refinement reached subnanosecond steps.

The residual contains `(p1 - p0) / step`. Its absolute rounding uncertainty
increases as the dimensionless step decreases, even when the represented
momenta are accurate. The original fixed force threshold eventually requested
precision below that subtraction's resolution. At one failed 38.1 ps trial,
the remaining scaled residual was 5.88e-9, already near the arithmetic floor.

The convergence threshold now includes
`32 * epsilon * (abs(p1) + abs(p0)) / step` for each coordinate. This is an
arithmetic allowance for the represented momentum calculation, not a change
to the constitutive forces or trajectory tolerance. Actual residuals and work
errors remain reported. Some intermediate trials still fail nonlinear
convergence and are subdivided; the controller does not accept those trials.

## Independent continuous-source comparison

The native trajectories use the existing 200 requested times over 0.5 ms and
the same fixed cyclic tract, source, outlet column, thermal, viscous and
radiation states as the saved independent 164-state ODE. The oracle and its
generators are hash checked before reuse. The fresh native calculations are
recorded in [joint_accuracy_reference.json](joint_accuracy_reference.json).

| Local tolerance | Accepted substeps | Largest sampled source-pressure error | Final source-pressure error | Largest tract-pressure error |
|---|---:|---:|---:|---:|
| 1e-4 | 1156 | 3.82563 Pa | 0.08175 Pa | 0.003735 Pa |
| 2.5e-5 | 1880 | 1.31842 Pa | 0.10818 Pa | 0.001631 Pa |
| 6.25e-6 | 3046 | 0.55225 Pa | 0.19522 Pa | 0.000693 Pa |

The largest sampled error decreases, but the endpoint errors are not monotone.
An endpoint alone is therefore a poor measure of this oscillatory trajectory.
The table does not establish a global error bound or a converged utterance.
The tight runs reach a minimum accepted step of 1.192 ps during the artificial
nonequilibrium startup; this is not a proposed speech clock.

Across these runs, the direct cumulative work residual is at most 4.10e-16 J,
air residual at most 3.36e-21 m³, and every accepted local comparison meets its
requested tolerance. Physical and numerical dissipation are nonnegative.

All 15 focused source, contact, coupling and controller checks pass in Release
and ASAN/UBSAN builds. The controller check verifies allocation-free advance,
exact replay of accepted substeps, rejection of a different stage layout,
and rollback after an accepted prefix followed by trial-budget exhaustion.
The tight cold-start interval also passes ASAN/UBSAN: its 175 output columns
agree with Release except for three entries differing by at most two ulps.

## Moving-airway checks and remaining integration

The moving adapter supports `--joint-tolerance` and `--joint-max-step` with the
joint passive update. `tools/joint_accuracy_moving_reference.py` checks direct
work and air balance, actual accepted times, and state continuation over two
closing and two reopening intervals. It verifies that its source files and
executable remain unchanged during the run. The independent fixed-source
result above does not substitute for these moving checks.

The current executable has completed all four moving intervals at local
tolerance 1e-4. The [current moving receipt](joint_accuracy_moving_current.json)
records 3732, 572, 1040 and 592 accepted substeps respectively. Direct work
residual is at most 1.42e-16 J and air residual at most 3.82e-20 m³. Continuation
preserves graph and relaxation states exactly and checks source physical-state
round trips. The source files and executable were unchanged throughout this
run. All 31 historical JSON regressions also remain exact with this executable.
The full first closing interval also passes ASAN/UBSAN. It takes the same 3732
accepted substeps as Release, with a maximum common-time pressure difference
of 8.88e-9 Pa. Its independently accumulated work and air residuals are
1.42e-16 J and 1.96e-20 m³. The [completed validation receipt](joint_accuracy_validation.json)
binds the checked files to these results, including the unchanged speech WAV.

The earlier, separately identified executable completed all twelve moving
intervals across the three tolerances. Its four coarse JSON outputs match the
current executable exactly. Its incomplete source snapshot prevents treating
later-file source hashes as build provenance. The preserved binary, inputs,
outputs and this limitation are recorded in the
[exploratory refinement receipt](joint_accuracy_moving_exploration.json).
These finer results are not current-build validation or an independent ODE.

| Compared local tolerances | Closing maximum sampled pressure change | Opening maximum sampled pressure change |
|---|---:|---:|
| 1e-4 → 2.5e-5 | 0.07956 Pa | 0.22047 Pa |
| 2.5e-5 → 6.25e-6 | 0.02951 Pa | 0.05390 Pa |

The closing endpoint differences decrease from 0.04056 to 0.00806 Pa. Opening
endpoint differences increase from 0.01298 to 0.02460 Pa despite the decreasing
maximum sampled discrepancy. This is another reason to inspect complete
trajectories rather than infer accuracy from a final state alone.

The source-contact transition rule and complete utterance integration remain
unfinished. This controller provides an accuracy reference, not a production
speed result. Performance optimization remains parked by user instruction;
the eventual engine still must satisfy the original 10× realtime requirement.
