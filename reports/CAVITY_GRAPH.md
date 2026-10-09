# Pressure in passages that separate and rejoin

The [solid-airspace reference](PLANAR_AIRSPACE.md) retains air above and below
the tongue. Its passages can form loops, and a tongue closure can leave
separate fluid components. The previous `CavityNetwork::prepare` rejected
both cyclic graphs and disconnected graphs. It could not consume these
geometries while retaining their actual connections.

`CavityNetwork` now accepts both. This extends the native mean-pressure
engine; it does not connect the planar geometry to speech rendering yet.
The [latest listening assessment](LISTENING_SPATIAL_2026_10_09.md) remains:
November improved, Foxtrot still jammed around /tr/. There is no new listening
result, spectrogram, or claim of improved speech in this change.

## Preserved physical equation

Each node still stores `V + C p`. Backward Euler solves

```
(C_new p_new - C_old p_old + V_new - V_old)/dt
    + sum(outgoing passage flows) = prescribed inflow.
```

Every open passage retains the existing passive pressure law
`dp = R Q + loss rho Q |Q| / (2 A²)`, with independently specified forward
and reverse losses. Zero area means exactly zero flow. Parallel connections
are separate physical passages; their conductances both participate.
Reservoirs retain their specified pressure. Separate components do not gain
an artificial fluid connection.

The nonlinear mass residual and convex-potential Newton globalization are
unchanged. A tree still uses direct elimination. For other graphs, Newton's
positive Jacobian is solved with conjugate gradients, preconditioned by a
spanning forest plus storage and external-port conductances. The forest
helps solve the system; it does not replace the physical graph. All cycle
edges participate in every Jacobian multiplication, nonlinear residual,
passage potential and accepted flow. Jacobian products use pressure
differences to avoid cancelling large diagonal and off-diagonal terms near
uniform pressure.

Disconnected forest roots have zero-coupling bookkeeping links inside the
tree factorization only. Conjugate-gradient convergence is checked against
a freshly recomputed full Jacobian residual. Iteration limits return failure,
leaving accepted pressure and flow untouched. Scratch memory is allocated
in preparation, with no allocation during stepping. This is support for
needed physical topology, not a performance-optimization milestone.

## Verification

Native regressions exercise:

- A four-node ring with an independently known antisymmetric pressure mode.
  Maximum pressure error after one step: 6.83e−13 Pa. Removing a route would
  change this mode.
- Closing both contacts around one node and compressing that isolated
  chamber. Sealed flows remain exactly zero; other pressure changes stay
  below 1.43e−14 Pa.
- Reopening the loop: pressure equalizes, stored mass remains constant and
  energy decreases. Maximum charge discrepancy: 8.61e−23 m³.
- Disconnected chambers under compression, two parallel passages, and no
  allocations during cyclic stepping.

The full Release suite passed 14/14. The affected streaming, cavity-network
and mean-airway tests passed ASan/UBSan, 3/3. All four Foxtrot graph inputs also
passed under ASan/UBSan with the existing `-O1 -g` sanitizer configuration;
maximum difference from the Release pressures was 5.69e−14 Pa. The local
sanitizer receipt is `out/cavity-graph/sanitizer.json`.

The existing tree regressions remain in place, including pressure reversal,
nonlinear loss, closure/reopening, changing compliance, mass balance,
rollback, nasal routing, and streamed synthesis.

`tools/cavity_graph_reference.py` independently builds and solves 24 seeded
nonlinear graphs with loops, closed contacts, directional quadratic losses,
reservoirs, and sources. SciPy's root solver receives an independently
assembled residual and Jacobian. Maximum native/reference pressure difference
is 2.65e−9 Pa; maximum total mass discrepancy is 1.01e−18 m³ per step.
These fixtures use fixed geometry, not changing cell ownership.

Four saved Foxtrot poses are also converted into graphs using every strictly
positive cell and internal opening, without an aperture threshold:

| Pose | Air cells | Connections | Fluid components | Independent loops | Max pressure difference |
|---|---:|---:|---:|---:|---:|
| 5.025 s | 1,270 | 2,377 | 1 | 1,108 | 1.50e−8 Pa |
| 5.040 s | 1,270 | 2,378 | 2 | 1,110 | 3.15e−8 Pa |
| 5.071 s | 1,282 | 2,407 | 1 | 1,126 | 1.55e−8 Pa |
| 5.073 s | 1,286 | 2,416 | 1 | 1,131 | 1.20e−8 Pa |

Each is compared with an independently assembled sparse linear solve.
Maximum mass discrepancy within an individual fluid component is 3.47e−24 m³.
The closed /t/ snapshot starts with uniform rear pressure and nonuniform front
pressure: front equalization leaves the rear exactly unchanged in this run.
The other poses each have one connected air region.

For these closed, fixed, linear fixtures the energy identity is also checked:

```
E_new - E_old + 1/2 sum(C (p_new - p_old)²) + dt sum(R Q_new²) = 0.
```

The second term is backward Euler's numerical dissipation and the third is
passage dissipation. Maximum discrepancy is 1.13e−21 J. The graph pressure
solver therefore preserves the assessed geometry's connectivity and storage
under these tests, rather than joining all air into one pressure compartment.

Commands, input/output and source/binary hashes, reference-library versions,
and all individual results are in
[cavity_graph_reference.json](cavity_graph_reference.json). To reproduce:

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/cavity_graph_reference.py
```

The audit consumes the four saved snapshots produced by
`tools/planar_airspace_reference.py`. Native C++ adds no external dependency;
SciPy is used only by the offline independent reference.

## Limits before speech integration

The snapshot pressure test assigns existing circular-tube resistance to each
opening using the aperture and distance between cell centres. That is a
numerical fixture, **not physically justified planar viscosity**. Interior
grid faces are not tube walls; copying that resistance into the speech engine
would introduce artificial wall friction. Physical mean-flow resistance,
inertance and acoustic boundary terms must come from actual tissue geometry.

`CavityNetwork` still requires positive compliance and persistent node
identity. It does not transport stored air when a cut cell disappears,
splits, or merges. Passing zero volume with an arbitrary positive compliance
would merely hide that missing moving-contact treatment. The planar geometry
also still needs general resolution of cells containing disconnected air
pieces. The four audited snapshots do not contain such cells.

The production spatial-pressure speech route still uses section geometry,
and the acoustic representation remains a tree. Nasal flow, side cavities,
glottal coupling and radiation must survive its eventual replacement.
Conservative moving-state transport, tissue-volume preservation, transverse
anatomy, wall work, and new listening evidence remain required. Neither the
speech-realism goal nor the eventual speed gate is satisfied by this audit.
