# Connect the continuous source to the moving mouth

The shared-inertia continuous gas/tissue source now drives the actual moving
Bravo mouth, integrated lips and epilarynx controls. The same source flow enters
the tract, and the tract's pressure acts back on the gas and tissue during the
step. This replaces the inherited algebraic source **in an explicit research
mode**. It is not yet a complete utterance renderer or a new recording.

The first required connection was the glottis-to-first-throat-node half-section.
Its air inertia cannot be silently omitted or represented by an unrelated flow
state. The source now includes it in the same kinetic metric as its three gas
transfer coordinates and four tissue coordinates.

## Attached column and changing geometry

For a throat half-section of length `l/2` and area `A`, the inviscid linearized
column inertia is `L = rho * l / (2 A)`, in Pa s2/m3. Its kinetic energy is
`L Q_out^2 / 2`. The same outlet coordinate carries source flow and column flow;
adding `L` to that coordinate's metric also changes the coupled canonical
momentum and inverse inertia. The column's storage remains in the throat node,
so no air volume is duplicated.

`ContinuousGlottis::prepare` accepts an optional initial outlet inertia. An
optional endpoint value in `advance` prescribes its change. The default zero
inertia and unchanged-parameter advance preserve the preceding API behavior.
A changing parameter is handled symmetrically:

1. Retain canonical state and change the coefficient to its endpoint midpoint.
2. Advance the coupled source with that held midpoint coefficient.
3. Retain the accepted canonical state and change to the endpoint coefficient.

The two coefficient changes have measured Hamiltonian energy changes, recorded
as `geometry_work`. They are not assigned to viscosity or numerical damping.
The source balance becomes

```
energy_change + physical_losses + numerical_loss + load_work
    = reservoir_work + external_force_work + geometry_work
```

The graph wrapper carries the same receipt. Preparation allocates; advancing a
prepared source/graph with changing column inertia passes the allocation-free
control. Invalid endpoint inertia and failed solves preserve accepted state.
A sealed outlet carries zero column flow and zero geometry work in the tested
coefficient change. Reverse flow retains positive exit dissipation.

The [independent column receipt](continuous_outlet_column_reference.json)
uses a separate continuous-field ODE coupled to a cyclic three-node tract.
The variable-column ODE includes both the additional inertia and its time
change in the momentum equation. Its independent geometric power is
`-Q_out^2 * dL/dt / 2`, evaluated at fixed canonical state. It is not obtained
from the native step-work receipt.

At the last step halving, endpoint-error ratios are 4.006 for a constant column
and 4.006 for a changing column. The changing-column work-error ratio is 3.998;
its finest work discrepancy is 3.18e-12 J against an accumulated 5.16e-7 J.
Independent quadrature/tolerance refinements differ by less than 1.96e-9 in the
declared normalization. Coarse startup pressure transients remain underresolved:
at a 1.25 microsecond step, maximum sampled state errors are about 0.0371 and
0.0343, although endpoint errors are about 5.33e-5. Endpoint agreement does not
establish accuracy throughout those transients.

## Actual articulated geometry

The explicit native mode is

```
ducklips_planar_inertive_motion_probe --lip-profile --epilarynx \
    --continuous-folds bodies/reference.body
```

It uses the existing two mouth regions, four integrated lip sections and eight
throat sections, including their changing volumes and certified region/cap
attachments. The final throat half-section's inertia uses its actual endpoint
area and length. The source loads the body's fold parameters and retains the
explicit full-span glottal profile used in the recent contact work.

After the existing geometry input, this mode reads rest area, reservoir pressure
and tension. Initial gas effort matches the initial throat effort; initial flow
is zero, with a one-time 1 micrometer lower-fold displacement. This initializes
an experiment; it does not prescribe a source waveform or an oscillation phase.
Rest area and tension remain held within these controls. External tissue force
is zero. The source uses the resistance-fitted update; source contact requiring
subdivision is not automatically resolved by this adapter.

With `--initial-state`, graph charges and branch momenta are followed by two
gas efforts, three instantaneous flows, and the lower positions, lower velocities,
upper positions and upper velocities (two values each). The adapter reconstructs
the source at the identical boundary geometry. This serialization has measured
round-trip error; it is not claimed to be bit-exact canonical restoration.

## Moving-control evidence

The [moving-airway receipt](planar_continuous_reference.json) covers Bravo
closing at 1.284–1.286 s and reopening at 1.374–1.376 s. These are separately
initialized controls, each carrying source and tract state across two 1 ms
intervals. Closing begins at 0 Pa, reopening at 700 Pa; both have a 700 Pa
reservoir, 3e-6 m2 rest area and tension 1. The geometry comes from the saved
native articulator export, not from the new source's own full-utterance run.

Twelve runs use 256, 512 and 1,024 subdivisions per millisecond, totaling
7,171 accepted steps. Independent Boolean geometry, lip quadrature and the
full graph stage solve check 3,079 geometry snapshots at the first refinement.
The independent graph receives the **recorded coupled source flux**. This tests
geometry and the connection; it is not an independent ODE for the entire moving
source/mouth system. The separate column ODE above supplies that component's
continuous-motion check.

| Measured maximum | Value |
|---|---:|
| Independent graph pressure difference | 3.13e-9 Pa |
| Local source/network stored-air residual | 9.44e-21 m3 |
| Per-interval combined stored-air residual | 8.22e-20 m3 |
| Source step work residual | 4.54e-17 J |
| Source/tract port-work mismatch | 5.30e-23 J |
| Source mass round-trip at interval boundary | 2.49e-24 m3 |
| Source canonical-momentum round-trip | 4.17e-17, in each coordinate's SI momentum unit |
| Source energy round-trip | 3.39e-21 J |

Combined mass includes the two source gas moments and the tract's actual volume
plus acoustic charge. Its only external transfers are reservoir inflow and lip
outflow. Fully sealed lip intervals carry exactly zero exterior flow. A step
ending at closure can still transport air before the endpoint; its average
flow must not be confused with endpoint flow.

The temporal accuracy limitation remains visible:

| Refinement, subdivisions/ms | Closing endpoint change | Opening endpoint change |
|---|---:|---:|
| 256 to 512 | 51.93 Pa | 39.01 Pa |
| 512 to 1,024 | 39.73 Pa | 4.08 Pa |

These moving controls are not temporally converged. The mouth/throat graph is
still inviscid with an ideal pressure outlet: wall/thermal losses, radiation,
general changing-region transport, continuous physiological actuation and
whole-utterance integration remain required. The tests do not calibrate a
human body or establish intelligibility. Performance work remains parked.

## Validation and reproduction

Eleven focused CTests pass in Release and eleven pass under ASAN/UBSAN,
covering the source, contact, onset, transfer rebasing, graph work, fitted decay
and changing-column coupling. Additional sealed-column and reverse-flow fixtures
pass in both builds. Seven historical planar adapter outputs reproduce their
entire JSON exactly, and the legacy equilibrium/closed-source controls pass.
All four actual moving-source intervals also pass under ASAN/UBSAN at 256
subdivisions/ms. The [validation receipt](continuous_moving_airway_validation.json)
pins the checked binaries, outputs, source inputs and test logs.

Use the existing scientific Python environment with NumPy, SciPy and Shapely:

```sh
cmake --build build --target ducklips_graph_continuous_glottis_probe ducklips_planar_inertive_motion_probe
python tools/continuous_outlet_column_reference.py
python tools/continuous_outlet_column_controls.py
python tools/planar_continuous_reference.py
```

The input trajectory/profile exports and their hashes are recorded in the
receipts. Outputs live under `out/continuous-outlet-column`,
`out/continuous-outlet-column-controls` and `out/planar-continuous`.
