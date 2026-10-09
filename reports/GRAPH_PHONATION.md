# Moving folds coupled to the graph; glottal work accounting remains open

The new `GraphPhonation` library component connects the existing two-mass
vocal folds to the total-flow throat/mouth/lip network. Fold opening and
pressure feedback now agree within the same integration step. Fold positions,
velocities and inlet momentum can continue across saved geometry intervals.
The planar adapter loads the fold tissue parameters from the body file.

The numerical coupling works in the tested controls. A separate work audit
also exposes a physical gap in the inherited fold-force model: during parts
of closure it gives the tissue more aerodynamic work than the algebraic
orifice supplies in that interval. There is no internal glottal air-storage
or displacement state to account for that difference. This remains a
reference integration, not a completed physical speech engine. No new
recording or listening improvement is claimed.

## Coupling pressure at the correct location

The pressure at the first throat storage node is not the pressure immediately
above the glottal constriction. Their difference accelerates the air in the
inlet half-section. The source port already retains that path's momentum.
For an open glottis the pressure applied to the fold-force law is the drop
across the orifice itself:

```
D = P_lung - P_throat_average - (lambda_new - lambda_old) / dt
  = rho U |U| / (2 A_midpoint^2)
```

For a trial `D`, the existing fold mechanics advance half a step and supply
`A_midpoint`. The network's free pressure and unit-flow response determine
the resulting inlet flow and actual orifice drop. A bounded scalar solve
makes the trial and actual drops agree. The second fold half-step then
uses that same pressure. The network receives the resulting total flow.

The solve uses the pressure bound provided by the free network load and
inlet momentum. It checks the residual before accepting the result; no
claim of global root uniqueness is made for arbitrary configurations.
An unresolved root or invalid numerical state returns failure without
changing the accepted folds, air state or diagnostics. Stepping allocates
no memory after preparation.

The existing symplectic Euler fold mechanics are unchanged. There are two
mechanical half-steps per network step, so the coupled method is first
order overall; implicit pressure feedback does not make the inherited
mechanical integrator second order.

An exactly closed glottis passes zero flow and removes the source-path
momentum. The remaining static pressure difference acts on the closed
geometry; the numerical stop impulse is not interpreted as a sustained
fold pressure. The downstream air remains stored.

`FoldState` contains both sides' lower/upper positions and velocities plus
the current rest area. Restoration requires finite values and preserves
state on failure. Tissue parameters remain those of the configured body.
The stream supplies rest area, tension, lung pressure and optional surface
pressure for each step; it does not prescribe a waveform or fold phase.

## Native mouth integration

`ducklips_planar_inertive_motion_probe --lip-profile --epilarynx
--vocal-folds bodies/reference.body` uses the same 14-node geometry as the
[held-aperture source control](PLANAR_SOURCE.md). The source input becomes:

```
rest_area reservoir_pressure
fold_tension surface_pressure
```

`--initial-state` adds the nine fold-state values after node charges,
network momenta and inlet momentum. Output records midpoint aperture,
transglottal pressure, coupling residual and the complete fold state.

The [moving-geometry receipt](planar_phonation_reference.json) covers actual
Bravo closure at 1.284–1.286 s and reopening at 1.374–1.376 s. Each is a
separately initialized two-interval control. Both use rest area 3e-6 m²,
tension 1 and a 700 Pa reservoir. A one-time 1 µm displacement of each lower
fold seeds motion; surface-pressure forcing is zero. The mouth geometry is
from the saved native trajectory, but these laryngeal controls are not a
replay of the original recording's complete physiology.

All physical state continues exactly across the shared interval boundary.
Independent Boolean geometry, lip quadrature, full network stages and a
separate nonlinear fold/source solve audit 12 runs, 451 steps and 1,365
geometry snapshots. The oracle carries its own air and fold states.

| Maximum discrepancy | Value |
|---|---:|
| Network pressure | 1.12e-8 Pa |
| Source flow | 2.61e-16 m³/s |
| Fold position | 4.16e-16 m |
| Fold velocity | 3.78e-13 m/s |
| Transglottal pressure | 2.19e-7 Pa |
| Local stored-air balance | 6.44e-21 m³ |
| Complete two-interval stored-air balance | 2.72e-20 m³ |

The oracle agreement is a discretization check. The moving inviscid controls
are not temporally converged: the 32-to-64-subdivision endpoint change is
44.01 Pa during closing and 4.86 Pa during opening. Their undamped lip modes
and coarse mouth partition retain the limitations of the earlier controls.

## Independent fold/source checks

The [source receipt](graph_phonation_reference.json) includes a separate
three-node network with a cyclic path and an outlet. A 25 ms independently
evolved discrete reference covers opening and near contact. It agrees with
the native endpoint pressure within 2.11e-8 Pa. The smallest midpoint
aperture is about 5.22e-40 m² under the existing smooth contact law; no new
aperture floor was added.

An independent continuous ODE in pressure, flow and mechanical coordinates
checks the smooth first 3 ms. Its final three refinement ratios are about
1.997, 1.998 and 1.999, approaching first-order convergence. At 4,608 steps
over that prefix, endpoint errors are 0.00142 Pa, 2.16e-9 m³/s and 4.84 nm
in fold position. This explicit ODE comparison does not claim to resolve
the singular fluid relaxation of near contact later in the control.

Native checks cover exact state continuation, finite state restoration,
rejection rollback, zero-flow hard contact, pressure feedback and stepping
without allocation. AddressSanitizer/UndefinedBehaviorSanitizer pass the
native coupling checks and all four moving Bravo intervals at 16 subdivisions.
Seven focused CTests, including the main native test suite, pass. All 24
previous held-aperture source controls reproduce their historical JSON
output exactly. Recorded source, binary, input and output hashes were
verified after the runs.

## What the work audit changes

Matching the equations is not sufficient to establish a physically complete
source. The [work receipt](fold_port_work_reference.json) compares, over each
step of the same 25 ms control:

- Orifice port work: `dt * D * U`.
- Aerodynamic work on both folds: the inherited lower-surface pressure
  force times each accepted half-step displacement.

At 2,400, 4,800 and 9,600 steps, the peak amount by which fold power exceeds
orifice power is 8.025, 8.032 and 8.034 mW. The integrated positive shortfall
is 10.648, 10.684 and 10.698 µJ. Refinement preserves this interval mismatch.
Near the peak in the finest run, the folds receive about 8.15 mW while
through-flow supplies 0.116 mW, at a midpoint opening of 5.73e-9 m².
Each independent work trajectory is checked against the native fold/air
endpoint. The finest native pressure disagreement is below 1.51e-8 Pa.

This is not net energy creation over the whole control: the total orifice
work is about 2.442 mJ, whereas net aerodynamic fold work is about 47.72 µJ.
It is a missing account of how energy is exchanged during parts of closure.
Calling the full orifice work irreversible heat while also assigning
mechanical work to the folds would also double count that energy transfer.
The earlier held-aperture source energy identity remains valid for its
fixed orifice; it does not prove conservation of the coupled fold system.

The missing internal glottal air volume is a concrete next modelling target.
When the upper opening closes while the lower fold moves, the airflow that
fills the changing space need not equal the flow leaving it. The current
algebraic orifice uses one through-flow and has no state for that difference.
A physical replacement needs pressure forces and displaced volume that do
matching work, including any fluid storage or inertial transfer. That is
an implementation requirement suggested by this audit, not a completed
replacement or a uniquely established explanation of every residual.

This gate takes precedence over treating the new source as ready for speech.
Tract losses/radiation, full-utterance topology transport and anatomical
calibration also remain unfinished. Performance optimization stays parked.

Follow-up: [finite glottal storage](GLOTTAL_STORAGE.md) now independently
checks mass storage and matching pressure-volume work on the tissue. It
uses supplied flow ports and has not replaced this coupled source. The
work gap measured here therefore remains an open full-system check.

## Reproduction

```sh
cmake --build build --target ducklips_graph_phonation_probe ducklips_planar_inertive_motion_probe
python tools/graph_phonation_reference.py
python tools/planar_phonation_reference.py
python tools/fold_port_work_reference.py
```

The moving controls reuse the saved native poses, airway/lip exports and
initial network states documented in the held-aperture report. Receipts hash
those initial references. Outputs and work traces are under
`out/graph-phonation`, `out/planar-phonation` and `out/fold-port-work`.
The sanitizer receipt is `out/planar-phonation-asan/receipt.json`.
