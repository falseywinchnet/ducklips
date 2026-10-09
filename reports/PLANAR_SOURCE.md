# Connect the throat air column and an implicit glottal inlet

The moving planar mouth and integrated lip reference now connect to the
body's eight existing epilarynx sections. A reusable `GraphOrifice` source
injects total volume flow into this same pressure network and receives its
pressure feedback. Source-path momentum continues across saved intervals
along with the mouth, lip and throat states.

This is a held-aperture integration control. It does not yet move the vocal
folds, replace the production speech path or produce a new recording. The
latest listening result remains November improved, Foxtrot jammed at /tr/.

Subsequent [moving-fold integration](GRAPH_PHONATION.md) adds dynamic tissue
feedback and identifies a remaining glottal work-accounting gap. The fixed-
aperture results below retain their original scope.

## Source law and stored air

The source solves against the network's average pressure over the step,
which is conjugate to the held injected flow. The network is advanced twice:
once with its existing state, geometry displacement and any series sources,
and once from zero with unit inlet flow. Their linear superposition gives

```
P_load = P_free + Z_network U
```

The glottis-to-first-section-centre air path has inertance `L`. It carries
its own canonical momentum `lambda`, advanced by implicit midpoint:

```
U = (lambda_new + lambda_old) / (2 L)
(lambda_new - lambda_old) / dt = P_lung - P_load - rho U |U| / (2 A_glottis^2)
```

Substitution yields the existing signed Bernoulli quadratic with effective
impedance `Z_network + 2 L/dt` and pressure history `2 lambda_old/dt`.
The zero-inertance limit is an ideal nodal orifice. An exactly closed glottal
aperture has zero flow and discards only that absent source passage's
momentum. It does not discard the downstream stored air.

There is no separate mean-pressure cavity in this control. The same node
charge accounts for all entering air. Moving walls contribute `-dV/dt` to
that charge, and the inlet contributes `U`. Thus the tested global balance is

```
final sum(V + q) - initial sum(V + q) = integral(U_in - U_lips) dt
```

`GraphOrifice` allocates during preparation and steps without allocation.
Malformed source/geometry data or a failed solve leaves its accepted state
and source outputs unchanged. A negative network driving-point response is
rejected rather than treated as a passive held-flow load. General topology
changes and their state transport remain the caller's responsibility.

For fixed geometry, the network's source work is `dt * P_load * U`. Together
with the midpoint inlet inertia, this gives the discrete energy identity

```
change in stored acoustic + inlet kinetic energy
  = dt * P_lung * U - dt * rho |U|^3 / (2 A_glottis^2)
```

The source aperture and inlet inertance are held during a step. Network
geometry uses its existing two Gauss stages; the overall source coupling
is second order. This fixed-geometry energy identity is not a claim that
moving walls do zero work or that vocal-fold mechanics have been verified.

## Throat connection

The adapter's `--epilarynx` input takes section areas and lengths from the
unchanged native NATO airway export. The eight sections are ordered from
the mouth seam back toward the glottis. Their total length is 21.875 mm,
exactly the core's inlet seam coordinate for this body. Their volume is
included once, outside the planar core.

Core-to-throat attachment inertia includes the measured exposed cap and
the throat section's half length. Subsequent links include both adjacent
half lengths. The remaining half of the first section at the glottis is
included in the source inertance. When several cap regions feed the same
throat tube, entry shares use cap area; they do not create duplicate tubes.
The actual controls have one throat attachment and two lip attachments.

The fixed graph contains 14 storage nodes and 15 network branches, plus the
inlet's momentum. Complete lip closure retains its inner shoulder volume.
The throat sections require positive storage; disappearing throat cells are
not silently regularized. Both caps must remain open with certified region
ownership. The independent audit checks this for the tested intervals.

After the existing lip input, the epilarynx input is:

```
section_count attachment_count
[area0 length0 area1 length1] per section, seam toward glottis
[initial pressures] per section
[domain_edge region core_distance] per attachment
```

`--glottal-orifice` appends held aperture and reservoir pressure. With
`--initial-state`, canonical node charges, network momenta and then inlet
momentum follow. This continuation requires the same endpoint geometry and
node/branch identities; it is not a remeshing operation.

## Evidence and its limits

The [moving-geometry receipt](planar_source_reference.json) uses the actual
Bravo controls already established by the
[profile adapter](PLANAR_LIP_PROFILE.md): closing at 1.284–1.286 s and opening
at 1.374–1.376 s. Each carries state through two 1 ms intervals. The closing
network starts at 0 Pa; the opening network starts at 700 Pa. Both use a
held 3e-6 m² glottal aperture and a 700 Pa reservoir. These are separately
initialized controls, not a full simulation from the beginning of Bravo.

Twelve runs at 16, 32 and 64 subdivisions independently audit 451 steps and
1,365 geometry snapshots. The oracle uses independent Boolean geometry,
axial quadrature and a full charge/momentum stage solve, with a bracketed
root of the nonlinear source law. It carries its own state between intervals.
Native state continuation reproduces volume, compliance, conductance, charge,
network momentum and source momentum exactly at each shared boundary.

| Maximum independent discrepancy | Value |
|---|---:|
| Network pressure | 1.63e-8 Pa |
| Source load pressure | 1.42e-9 Pa |
| Source flow | 2.14e-16 m³/s |
| Network branch momentum | 1.46e-13 Pa s |
| Source momentum | 1.59e-14 Pa s |
| Local stored-air balance | 6.37e-21 m³ |
| Per-interval global stored-air balance | 5.43e-20 m³ |

The largest balance residual across a complete two-interval control is
4.07e-20 m³ over the tested refinements. Fully closed lip passages carry
exactly zero exterior flow.

Additional native-only runs refine through 512 subdivisions per millisecond.
**The moving-control endpoint pressures are not yet converged at that rate.**
The last refinement changes the closing endpoint by 51.96 Pa and the opening
endpoint by 38.97 Pa. Agreement with an independent discretization does not
remove this timestep error. The undamped lip modes identified in the earlier
profile work remain present; these controls do not yet include tract losses
or physical radiation. No accuracy or speed claim for a speech renderer is
based on these coarse endpoints.

The separate [continuous source receipt](graph_orifice_reference.json)
compares a static three-node cyclic network with an independently integrated
nonlinear ODE in pressure/flow coordinates. The final three refinement ratios
are 3.85, 3.96 and 3.99, approaching the expected factor four for second order.
At 16,384 steps over 4 ms, maximum endpoint differences are 6.97e-5 Pa and
2.28e-11 m³/s. This is a convergence control, not a production timestep.

Native controls also check reservoir work against stored energy and orifice
heat, total injected volume, reverse flow, equilibrium, exact source closure,
rollback and allocation-free stepping. A stationary, closed-lip version of
the actual 14-node geometry holds uniform 700 Pa within 8.27e-11 Pa, with
source flow below 1.36e-18 m³/s. Sealing the source gives exactly zero inlet
flow even with a different reservoir pressure.

Seven previous controls reproduce their entire historical JSON output
exactly, including the old tree, uniform lip and integrated profile modes.
AddressSanitizer/UndefinedBehaviorSanitizer pass all four moving source
intervals at 16 subdivisions and the native source controls. Five focused
CTest targets pass for the source, graph, lip profile and planar geometry.
Recorded source, binary, input and output hashes were rechecked after the
completed runs.

## Reproduce and continue

```sh
cmake --build build --target ducklips_graph_orifice_probe ducklips_planar_inertive_motion_probe
python tools/graph_orifice_reference.py
python tools/planar_source_reference.py
python tools/planar_source_controls.py
```

Inputs are the existing native airway/profile exports and the six Bravo
poses from the previous profile report. Diagnostic outputs are under
`out/graph-orifice`, `out/planar-source` and `out/planar-source-controls`.
The sanitizer receipt is `out/planar-source-asan/receipt.json`.

The source port is now available in the library. Dynamic vocal-fold motion
and its mechanical coupling still need integration with this port. Physical
losses/radiation, adequate spatial resolution and full-utterance topology
transport remain necessary before replacing the speech route. The current
coarse mouth grid is not anatomically calibrated. Performance optimization
remains parked; no new audio or enhanced spectrogram is attributed to these
controls.
