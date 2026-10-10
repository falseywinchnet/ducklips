# Lip outlet viscosity through closure and reopening

The developing continuous-source airway now includes local viscous profile
states on the actual exterior lip path. Exact closure distinguishes a vanished
flat core from air that remains in the widening shoulder. The new controls
also expose a limitation of the current numerical coupling: alternating exact
drag decay with a pressure step can strongly overpredict driven throughput in
a stiff slit. Pressure and drag need a common update before this path is ready
for a full speech claim.

## Physical state at a seal

The existing positive viscous realization uses segment inertia L, full-branch
inverse inertia G, local half-gap a, and 24 profile momenta z_j. The branch owns
the plug-flow inertia. Local profile energy is

`E = sum(pole_j/weight_j * z_j^2)/(2 L)`.

When another part of the path seals, G becomes zero but local air can remain.
Through-flow is constrained to zero. Each surviving profile obeys
`dz_j/dt = -nu pole_j z_j/a^2`, dissipating its stored energy. Its pressure
impulse is supported by the closed-flow constraint, with zero flow work.
`SlitViscousStep.blocked_impulse` records that contribution. The object leaves
the branch's unobservable canonical momentum to the graph's closure handling;
the actual adapter resets it before reopening.

A vanished local gap is a different state. The explicit geometry tuple
`(segment_inertance, branch_conductance, half_gap) = (0,0,0)` denotes absence,
not a finite air inertia. Entering this state returns `-E_before` as the change
of the profile energy metric, then clears the local memory. Opening again
starts that site without resurrecting its previous profile. The metric change
is separate from viscous heat. In a resolved continuous seal the strong
pre-contact drag should already relax the core; abrupt endpoint retirement
alone is not proof of an accurate closing trajectory.

The shared graph now accepts surviving finite sites on zero-conductance
branches and explicit absent sites. Late failure remains transactional. Native
controls exercise blocked decay, partial retirement, reopening, and preservation
of the other site's remaining energy without allocation during graph advance.
Twenty-four independent stiff Radau ODE controls check blocked diffusion;
maximum scaled state error is 2.61e-14 and relative heat error 3.44e-15.

## Geometry of the exterior path

The explicit `--outlet-viscous` adapter option requires the existing
`--lip-viscous-width`. For the actual profile with place=0 and a positive flat
core, it uses one core site and 4, 8, 16 or 32 fixed-reference shoulder sites.
The shoulder extends from the core edge to the last storage cell's actual
volume centroid. When open, all weighted inertias are normalized to the
existing integrated path inertia; when the core is closed, its site is absent
and the finite shoulder sites remain. No aperture or area floor is introduced.

The independent spatial oracle uses an angle transformation to resolve the
shrinking shoulder scale, while the native state sites retain their fixed
reference abscissae. Twenty-five profiles cover the actual closing/opening
intervals plus activation values 1-10^-k, k=1..7: 100 native fixtures. Across
129 frequencies from 20 Hz to 200 kHz, maximum relative spatial error in the
viscous correction is 1.186% at four shoulder points, 0.1661% at eight,
0.03372% at sixteen and 0.01437% at thirty-two. These are spatial errors;
the separate finite relaxation-spectrum approximation remains as previously
reported. The geometry comparison is within 9.43e-13 relative.

This remains a local stationary-wall slit approximation. Narrow-region acoustic
models assume wavelengths long compared with cross-sectional scales; they do
not replace a full transient moving-wall flow calculation. The broader model
limits are described in [COMSOL's narrow-region documentation](https://doc.comsol.com/6.4/doc/com.comsol.help.aco/aco_ug_pressure.05.186.html).
The 0.03 m slit width is an explicit body approximation, not anatomical
calibration. The current helper intentionally rejects other outlet shapes;
general interior contact and remeshing need separate state transport.

## Actual moving controls

Each case is two continued one-millisecond Bravo geometry intervals coupled
to the continuous gas/tissue source, thermal throat storage, internal throat
and lip drag, and aperture radiation. Closing and opening have separate initial
states. These are not an utterance or a new recording.

| Case | Internal lip drag only, 512 to 1024 subdivisions/ms | With outlet drag, 256 to 512 | With outlet drag, 512 to 1024 |
| --- | ---: | ---: | ---: |
| Closing | 0.11336 Pa | 0.41984 Pa | 0.11364 Pa |
| Opening | 0.64207 Pa | 1.26300 Pa | 0.27300 Pa |

These values are maximum endpoint pressure differences across all 14 nodes.
Doubling spatial quadrature to sixteen changes the finest-run endpoint pressures
by at most 0.00004639 Pa in closing and 0.00021680 Pa in opening.
The lower opening peak difference does not establish uniformly better accuracy:
in the frozen inviscid tract diagnostic its projected error energy rises from
1.93e-12 to 9.39e-12 J. That basis excludes source/loss-memory feedback and
circulation invisible to node divergence; it is not the full damped spectrum.

The eight-point runs contain 7,171 accepted steps. Their independent replay
uses the recorded coupled-source flux and separately reconstructs geometry,
graph, radiation, thermal and viscous updates. It is not an independent
whole-source oracle. Maximum discrepancies are 8.03e-9 Pa in graph pressure,
8.81e-15 Pa s in viscous memory, 4.94e-19 J in step viscous dissipation and
5.32e-22 J in profile geometry work. Maximum interval air-balance residual is
6.06e-20 m³. Sixteen-point and sanitizer runs also pass the unchanged limits.

## Driven-flow limitation: the next coupling requirement

A new native fixed-slit diagnostic applies the current composition under a
constant pressure difference: half a viscous step, a pressure/transport step,
then half a viscous step. It starts from the correct continuous steady state.
An independent full-matrix fixed-point solve predicts the resulting numerical
steady throughput. This directly checks forcing, rather than only free decay.

The physical steady flow is `Q = pressure/R0`, with
`R0 = 3 L nu/a^2`. Define stiffness `x = dt R0 G`. The observed mean-throughput
ratios against that continuous solution are:

| x | Separated numerical flow / continuous steady flow |
| ---: | ---: |
| 0.01 | 1.000135 |
| 0.1 | 1.004283 |
| 1 | 1.135244 |
| 10 | 5.002173 |
| 100 | 50.000000 |
| 10000 | 5000.000000 |

The native time march and independent fixed-point solution agree within 1e-10
relative. The step work ledger closes even while the mean flow is wrong. At
large stiffness the endpoint velocity is almost zero, but the middle pressure
step has already transported too much air. Even the scalar resistance-only
case shows this: its ratio is `(x/2) coth(x/2)` and tends to `x/2`, not one.
The finite profile model adds its own frequency dependence at smaller x.

The finest actual closing control has 510 open-endpoint steps, all with the
endpoint DC decay rate times the step above one. Opening has 951 such steps
out of 2,048 open-endpoint steps. These endpoint diagnostics do not measure the
whole evolving step's error, but confirm that the stiff regime is encountered
in the current geometry rather than only in an unrelated synthetic example.

This diagnostic is not a numerical error estimate for the entire moving airway.
It establishes that the present separated scheme lacks the correct stiff
forced-flow limit. Further temporal refinement or more profile quadrature
cannot replace a pressure/drag update with the correct steady throughput.
That common update is the next integration requirement, before treating these
controls as evidence for a new speech renderer. No physical coefficient was
tuned to compensate for this numerical defect.

## Validation and reproduction

Eighteen focused Release tests and eighteen ASAN/UBSAN tests pass. All four
actual moving outlet intervals pass under sanitizers at 256 subdivisions/ms.
Twenty-seven historical adapter outputs reproduce their full JSON exactly,
including the preceding internal-lip path.

```sh
cmake --build build --target ducklips_slit_viscous_probe \
    ducklips_lip_viscous_probe ducklips_graph_continuous_glottis_probe \
    ducklips_planar_inertive_motion_probe ducklips_viscous_forced_flow_probe
python tools/blocked_viscous_reference.py
python tools/lip_outlet_viscous_reference.py
python tools/planar_continuous_reference.py --radiation \
    --thermal-slit-width .03 --viscous-slit-width .03 \
    --lip-viscous-width .03 --outlet-viscous
python tools/viscous_forced_flow_reference.py
python tools/radiation_regressions.py --output out/outlet-viscous-regressions \
    --include-radiation --include-thermal --include-viscous --include-lip-viscous
```

The [validation receipt](continuous_outlet_viscous_validation.json) verifies
nested artifact hashes, exact closure-state handling, quadrature comparisons
and the unchanged listening recording. Performance optimization remains parked.
Whole-utterance integration, other boundary/attachment paths, moving-wall
physics and fresh human listening remain unfinished.
