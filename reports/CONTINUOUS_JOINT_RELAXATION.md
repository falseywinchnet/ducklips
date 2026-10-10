# Common pressure, airflow and relaxation update

The developing continuous-source path now advances tract pressure, transported
air, viscous profile memory, thermal storage and aperture radiation in one
passive solve. The constant-pressure slit control reaches the correct steady
flow across the tested stiffness range. Source motion and moving geometry still
have measurable time error. This is a mechanics change, not a new speech
recording or a change to the spatial-pressure listening score.

## Why the update changed

The prior split update could dissipate exactly the energy it received while
transporting the wrong amount of air. With `x = dt R0 G`, its steady-flow ratio
was about 1.135 at x=1, 5.002 at x=10, 50 at x=100 and 5000 at x=10000.
Pressure accelerated the air in a separate substep; subsequent drag removed
endpoint momentum without correcting the already integrated transport.
The [prior control](viscous_forced_flow_reference.json) remains available.

`RelaxationGraphReference` advances the following fixed-metric equations:

- Node pressure is `P=q/C`. A thermal memory obeys
  `z'=lambda (k P-z)`, and its derivative is subtracted from node charge rate.
- Branch flow is `Q=G p`. A profile or radiation memory obeys
  `z'=lambda (k Q-z)`, and its derivative is subtracted from branch momentum rate.
- A branch can also carry a DC resistance `R0`. Multiple local profile segments
  share one flow and the branch's existing total inertia.
- Each memory stores `z^2/(2k)`. Physical dissipation is
  `R0 Q^2 + sum(z'^2/(k lambda))`, with node and branch memories included.

The update uses two-stage Radau IIA, with stage matrix
`[[5/12,-1/12],[3/4,1/4]]` and weights `[3/4,1/4]`.
Its fixed-metric energy identity includes a separate nonnegative numerical loss
`dt^2/32 * ||f1-f2||_H^2`. Numerical damping is never reported as viscosity,
thermal loss or radiation. The method and stiff stability properties follow
[Hairer and Wanner's Radau treatment](https://www.unige.ch/~hairer/preprints/coimbra.pdf);
the [energy-stability discussion by Hairer and Lubich](https://www.unige.ch/~hairer/preprints/gradientflow.pdf)
provides further context. The model-specific work identity is checked directly.

All memory stages and branch-flow stages are eliminated locally. The native
solve contains only two pressure unknowns per node. The tested actual airway
has 14 nodes, so this is a 28-unknown solve even with 192 thermal memories,
1128 viscous memories and one radiation memory. This preserves all those
states; it does not prune them. No rendering-speed claim follows from this
representation. Performance optimization remains parked.

## Integration and geometry

`GraphContinuousGlottis::prepare_joint` prepares the fixed ownership of every
memory after loss sites have been configured. `advance_joint` obtains the free
tract response and the response to unit source flow. Their work-conjugate mean
pressure loads the existing implicit continuous gas/tissue update. An accepted
passive solve then includes the source flow and all dissipation cross terms.
The source has not been converted to the passive Radau stage equations.

Prescribed geometry is frozen at the midpoint during the passive step. Old-to-
midpoint and midpoint-to-end metric changes account for geometry work in the
graph, thermal memories, viscous memories and radiation state. This midpoint
approximation does not inherit the third-order fixed-metric result for arbitrary
moving geometry. Exact closure retains the previously established distinction
between absent core air and surviving shoulder profile states. Zero conductance
means zero branch flow; surviving memories can still dissipate their energy.

The adapter exposes this path through `--joint-relaxation`. The old path remains
available for reproducible comparisons. These are fixed-topology Bravo controls
with state continuation across two intervals, not complete utterance evolution.

## Independent checks

The [generic reference receipt](relaxation_graph_reference.json) contains:

- 80 random fixed-metric cases checked against an independent, unreduced full-
  state Radau system. Maximum scaled state difference is 9.96e-13; maximum mean
  pressure difference is 3.87e-12 Pa.
- A transient checked against both a matrix exponential and an independent
  continuous ODE. Those oracles agree within 6.81e-17 in energy coordinates.
  Native errors at 16, 32, 64, 128 and 256 steps are 1.205e-5, 1.517e-6,
  1.889e-7, 2.354e-8 and 2.937e-9: approximately eightfold reduction per halving.
- Driven slit tests from rest and from steady flow at x=0.01, 0.1, 1, 10, 100
  and 10000. All measured steady-flow ratios are within 2e-14 of one.
  The pressure boundary uses a very large compliant reservoir; its largest
  pressure deviation is 2.34e-23 Pa. It is not a mathematically eliminated clamp.

The [continuous source comparison](continuous_joint_reference.json) uses fresh
native trajectories against the saved independent 164-state ODE. The saved
trajectory and all its independent generator files are hash-verified against
the original receipt; the ODE was not integrated again for this change. Its
quadrature/tolerance control differs by at most 5.49e-9 in the stated scales.

| Steps in 0.5 ms | Maximum sampled scaled error | Endpoint scaled error |
|---:|---:|---:|
| 200 | 0.0596675 | 0.00219647 |
| 400 | 0.0370847 | 0.000556797 |
| 800 | 0.0147976 | 0.000213685 |
| 1600 | 0.00399405 | 0.000182917 |

The finest maximum is the lower glottal gas effort near 82.19 microseconds,
about 3.994 Pa on its 1000 Pa normalization. The endpoint maximum is also a
source effort. The prior split path's finest maximum was 0.00399451 and its
endpoint error was 0.000175406. The common passive solve therefore does not
establish a general improvement in coupled-source trajectory accuracy.
The finest passive numerical loss is 2.49e-10 J, reported separately from
physical losses. Combined source/tract work residual stays below 4.33e-16 J.

The [moving-airway audit](planar_joint_reference.json) reconstructs geometry
independently and solves the complete sparse stage equations without the
native memory elimination. It propagates its own passive state across the
256-subdivision intervals, including geometry work and closure. It uses the
recorded coupled source flux, so it is not an independent nonlinear source
trajectory. Higher subdivisions check mass, work and refinement without
repeating that full stage reconstruction. A [finer run](planar_joint_fine_reference.json)
extends the comparison to 2048 subdivisions per millisecond.

The independent pressure comparison is within 1.46e-7 Pa across 2054 geometry
snapshots. Maximum energy-coordinate difference is 2.54e-12; physical
loss differs by at most 9.03e-18 J and geometry work by 7.91e-19 J.
The 256/512/1024 runs contain 7171 accepted steps, with combined source/tract
air residual at most 6.21e-20 m3.

| Successive subdivisions per millisecond | Closing endpoint difference | Opening endpoint difference |
|---|---:|---:|
| 256 to 512 | 0.918697 Pa | 1.51040 Pa |
| 512 to 1024 | 0.447036 Pa | 2.44321 Pa |
| 1024 to 2048 | 0.103139 Pa | 1.16791 Pa |

These are differences between native trajectories, not errors against a known
continuous moving-geometry solution. Closing improves at the final refinement;
opening is nonmonotonic and is not claimed converged. Its final 1.16791 Pa
difference lies in the last lip storage section, while the largest mouth-node
difference is 3.26e-5 Pa and the largest throat-node difference is 0.00265 Pa.
Physical loss and the work ledger alone cannot settle that local time error.

## Validation and remaining work

The [collected receipt](continuous_joint_validation.json) verifies artifact
hashes, 20 focused Release and 20 AddressSanitizer tests, the source-control
rerun after adding the fixture protocol, four actual moving sanitizer controls,
and 31 exact historical JSON regressions. The moving sanitizer outputs are
compared to Release within stated tolerances, not claimed bit-identical.
Native controls include no allocation after preparation, rollback after late
invalid geometry, shared memories, work exchange and exact closure.

The next accuracy work is source gas-pressure evolution and moving-geometry
trajectory error, with the independent references retained as controls.
Yielding walls, remaining unmodeled loss geometry, sustained source contact
under the full moving tract, remeshing and whole-utterance integration remain
unfinished. No word is newly scored as intelligible. The latest listening
artifact remains `out/nato-spatial-pressure.wav`, SHA-256
`8284f366f2df5f3496d3139f7cbbaf4a8947ad839737d6a7e40b78635083b307`.
