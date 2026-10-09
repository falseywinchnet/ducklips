# Continuous vocal-fold source with acoustic back pressure

`GraphContinuousGlottis` connects the developing shared-inertia, continuous-
pressure source to `MovingGraphReference`. The tract receives the source's
accepted outlet flux, and the source responds to the tract's pressure during
the same step. This is the new continuous source, not `GraphPhonation`,
`GraphContourGlottal`, or the earlier spatial-pressure speech renderer.

The connection is a small-network integration reference. It does not yet
replace the articulated speech path. The latest listening evidence remains
November improved and Foxtrot jammed at /tr/ in the
[spatial-pressure recording](LISTENING_SPATIAL_2026_10_09.md).

## Eliminate the tract at its source port

For prescribed positive stage compliances and nonnegative inverse inertances,
the acoustic graph is linear in a nodal source held over the step. Two graph
solves give its free evolution and response to unit source flow. At the
source node,

```
average outlet effort = free effort + load slope * average source outflow.
```

The continuous source includes this affine relation in its nonlinear momentum
solve. Once that solve accepts a flux, superposition recovers the graph's
endpoint charge, branch momentum, average pressure, and average branch flow.
There is no previous-step pressure lag or prescribed glottal waveform.

The shared flux is reference-volume mass flux. The graph's pressure variable
is its **linearized** mass effort, compatible with the source's nonlinear
enthalpy effort to acoustic order. This is not an exact finite-density
volume-flow/thermodynamic-pressure interface. Large tract compression or
nonlinear convective acoustics requires a different graph constitutive model.
The source still requires a nonnegative affine load slope; a stage combination
that produces a negative slope is rejected, not clipped.

Prepared graph topology stays fixed. Cycles, parallel branches and disconnected
components remain in the graph. Each accepted advance allocates no memory.
Failed graph solves, rejected source solves and invalid preparations preserve
all previously accepted source/network state and receipts.

## Separate transferred work from exit loss

The source's historical `load_work` combines affine-load work and its kinetic
exit loss. The receipt now exposes those components separately as `port_work`
and `outlet_loss`, preserving their sum. This matters because energy sent to
the tract must not be counted as dissipated energy.

```
W_port = dt * average outlet flux * average outlet effort
load_work = W_port + D_exit

Delta H_source + D_fluid + D_tissue + D_exit + D_numerical
    = W_reservoir + W_tissue - W_port + source residual

Delta H_graph = W_port + W_nodal + W_series + graph residual
```

The graph energy equation above applies to **fixed metrics**. Prescribed
moving metrics also do geometric work. The wrapper exposes port and external
source work for moving metrics but does not label their sum as the graph's
energy change. Its graph evolution is checked against separate continuous
moving-metric equations.

Reverse flow reverses the transferred work while retaining nonnegative exit
loss. A sealed source outlet has exactly zero flux, port work and exit loss;
the attached graph can continue its own motion.

## A roundoff discontinuity found by the coupling

The [previous source correction](CONTINUOUS_PHONATION_ONSET.md) suppressed
discrete-gradient corrections below an endpoint-arithmetic uncertainty scale.
That removed amplified cancellation near steady throughflow. Its hard switch
nevertheless introduced a force jump at the threshold. A preloaded-tract
startup at 0.3125 microseconds per step stopped at step 118: the tissue force
residual alternated around the switch and stalled above the unchanged Newton
tolerance. The passage was open; this was not a contact event.

Write `e` for the endpoint energy discrepancy, `epsE` for the existing
arithmetic scale, and `dz` for the configuration change. The correction now
uses only the resolved excess:

```
resolved discrepancy = sign(e) * max(abs(e) - epsE, 0)
gradient correction  = resolved discrepancy * dz / dot(dz, dz)
```

This is continuous at the threshold. When the correction is applied, it leaves
at most `epsE` of that energy discrepancy uncorrected. The existing guard for
vanishing configuration increments remains. It does not relax the force tolerance. The arithmetic
scale remains an estimate, not a certified error bound; the actual energy
residual is still reported. A native regression reproduces the formerly
failing startup through 140 steps.

## Independent checks and contact experiment

The independent checker joins the previously separate open continuous-field
ODE to a charge-and-momentum graph ODE and integrates it with Radau. Native
source/graph trajectories are compared under temporal refinement with both
fixed and smoothly varying metrics. A second independent calculation solves
the full graph collocation system, carrying its own charge/momentum state
and using each accepted native source flux. This checks the affine elimination
independently of the native pressure-only solve.

The preloaded initial graph excites fast source pressure motion. The original
10, 5 and 2.5 microsecond trials had nonmonotonic maximum trajectory error,
despite small work residuals. They are not an adequate accuracy demonstration.
The refinement receipt retains the coarse results alongside the finer ones.

The [independent receipt](graph_continuous_glottis_reference.json) contains
twelve 2 ms trajectories. The following errors are maximum absolute errors in
the declared mixed-state normalization, not relative percentages: source
effort 1000 Pa, area 1e-6 m2, flow 1e-4 m3/s, area rate 0.001 m2/s, graph charge
1e-8 m3, and branch momentum 0.1 Pa s.

| Maximum step | Fixed metrics: whole trajectory | Moving metrics: whole trajectory |
|---|---:|---:|
| 10 microseconds | 0.0623590 | 0.0583421 |
| 5 microseconds | 0.1056100 | 0.0993640 |
| 2.5 microseconds | 0.0697023 | 0.0652036 |
| 1.25 microseconds | 0.0372977 | 0.0348187 |
| 0.625 microseconds | 0.0115004 | 0.0107631 |
| 0.3125 microseconds | 0.00293191 | 0.00274383 |

The final whole-trajectory refinement ratios are 3.92248 and 3.92264. Endpoint
errors at 2 ms are 3.28914e-6 and 2.86435e-6, with final refinement ratios
4.00141 and 4.00049. The independent 32/48-point quadrature and tolerance
refinement changes the sampled reference states by at most 9.67e-10 in this
normalization. The source's final state changes by 2.35/3.01 normalized units
when back pressure is omitted in the independent fixed/moving cases: the
connection materially feeds back into the tissue/gas motion.

Across these runs, the independent full graph collocation differs by at most
3.44e-11 Pa in pressure-equivalent charge. Maximum coupled mass residual is
4.87e-21 m3. For fixed metrics, the largest cumulative joint work residual is
4.99e-15 J. These small balance residuals coexist with the coarse accuracy
failures in the table; they are separate checks.

The [local tangent diagnostic](continuous_tract_timescales.json) identifies
source/graph modes near 63,975 Hz and 144,193 Hz at that nonstationary initial
state. Their local decay rates are about 22,117/s and 37,339/s. Halving the
Jacobian difference step changes those frequencies by less than 0.000001 Hz.
These are modeled internal pressure timescales, not measured human sound
frequencies or an equilibrium-stability result. A 10 microsecond step cannot
resolve a 6.94 microsecond oscillation. This explains why conservative coarse
steps can still misrepresent the initial pressure exchange. It does not
establish the cause of artifacts in the older speech recording.

A separate experiment attaches the source to an eight-cell tube, 17.5 cm long
with uniform 2 cm2 area, and an ideal zero-effort outlet. The graph starts with
zero compression and the initial source outflow throughout its branches.
The source starts from the saved pre-contact state at 1400 Pa reservoir pressure.
This tube is a declared control geometry, not a measured vowel or the actual
articulated Foxtrot airway. It has no wall loss or radiation resistance.

The experiment selects a dissipative source step only when an exterior-port
constraint changes. A rejected step bisects both source and graph; each
substep computes a fresh affine load from its own initial graph state. Failed
recursive trials are transactional. This frozen-geometry research driver
copies scratch storage; its allocation behavior is not that of the single-
step public wrapper. Moving-geometry subdivision would require new stage
metrics at the actual substep times.

Numerical damping, physical losses, port work, total mass and total work are
recorded separately. Repeated closure alone does not establish a converged
trajectory or a settled phonation cycle. The open-domain oracle does not
validate the contact constitutive law.

All four [24 ms coupled contact runs](continuous_tract_motion.json) complete,
but **their contact trajectories are not yet converged**:

| Maximum step | Port closures / reopenings | Added numerical loss |
|---|---:|---:|
| 10 microseconds | 3 / 3 | 1.43679e-6 J |
| 5 microseconds | 3 / 3 | 2.50936e-7 J |
| 2.5 microseconds | 4 / 4 | 5.99554e-9 J |
| 1.25 microseconds | 3 / 3 | 2.84846e-9 J |

These count individual port-constraint changes, not complete phonation
periods. In the 2.5 microsecond run the lower port seals again near 20.7 ms;
that extra closure is absent on the finest grid. Across adjacent grids the
maximum endpoint tract pressure differences are 666.48, 228.78 and 107.35 Pa.
The final pair's upper areas differ by 2.14 mm2, with endpoint source outflows
of about 199 versus 83 mL/s. These are material motion errors, even though
the work/mass checks hold. Small port and energy residuals do not justify
calling this a working voice.

Maximum cumulative joint work residual across the four runs is 5.74e-14 J;
maximum coupled mass residual is 1.18e-20 m3. On the finest grid these are
1.29e-15 J and 1.18e-20 m3. Tract pressures range from roughly -1024 to 1017 Pa
on that grid. Resolving contact-driven pressure motion under tract feedback
is the next accuracy requirement, before connecting this source through the
moving articulated airway and assessing pronunciation.

## Reproduce

```sh
cmake -S . -B build
cmake --build build --target ducklips_graph_continuous_glottis_probe ducklips_continuous_tract_probe
ctest --test-dir build --output-on-failure -R ducklips_graph_continuous_glottis
OPENBLAS_NUM_THREADS=1 python tools/graph_continuous_glottis_reference.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_tract_experiment.py
OPENBLAS_NUM_THREADS=1 python tools/continuous_tract_timescales.py
```

The scripts record source, binary, seed and output hashes. Prior onset/contact
receipts remain historical measurements of their recorded source versions.
Eight focused native CTests pass in both Release and AddressSanitizer/
UndefinedBehaviorSanitizer builds. They cover the coupled port, direct graph
evolution, allocation-free stepping, rollback, sealed/reverse flow, exit-loss
separation, the new roundoff-switch regression, and existing source/contact/
throughflow controls. The roundoff regression also fails when compiled against
the preceding hard-switch source, confirming it detects the repaired defect.
This milestone does not establish pronunciation, perceptual realism or the
requested rendering throughput. Performance optimization remains parked.
