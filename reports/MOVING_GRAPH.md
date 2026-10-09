# Canonical pressure and momentum on a general graph

`MovingGraphReference` extends the canonical moving-network equations to
cycles, parallel branches and disconnected components. This is needed to
connect several planar mouth regions to one shared lip tube: those regions
can already communicate inside the mouth, so the common outlet creates a
cycle. Separate ideal pressure loads cannot represent its shared storage
and pressure.

This change is a small-network integration reference. It does not replace the
production speech path or the existing tree solver, attach the actual lips,
or produce a new recording. The latest listening result remains November
improved and Foxtrot jammed around /tr/ in the spatial-pressure recording.

## Equations and state

The state remains nodal compression volume `q` and branch momentum `lambda`:

```
q'      = -B G(t) lambda + nodal_flow
lambda' = B^T C(t)^-1 q + series_pressure
p       = C(t)^-1 q
Q       = G(t) lambda
```

Every branch retains its momentum, including circulation around a loop that
has zero net nodal flow. The two-stage Gauss method evaluates positive nodal
compliance and nonnegative branch conductance at both stage times. Eliminating
the stage momenta gives a dense system of `2 * node_count` stage pressures;
it does not eliminate physical branch states. Compliance scaling and partial
pivoting precede a residual check against the original scaled equations.
There is no dimensional absolute pivot or aperture floor.

The caller can supply nodal volume flow and series pressure held over the
step. The outputs include average pressure and average branch flow, allowing
source-work and stored-air balance to be checked with the same discretization.
All stepping scratch is allocated in `prepare`. Invalid inputs, an
unrepresentable matrix, a failed solve or a nonfinite candidate leave the
accepted state and average outputs unchanged.

Zero `G` stops branch transport. As with the existing tree reference, physical
contact handling belongs to the caller: this integrator does not erase the
latent momentum or decide what happens to its energy. It requires positive
`C` and fixed topology. A disappearing air volume therefore needs separate
state transport; a zero-compliance closed lip is not supported by pretending
it has a small positive volume.

## Verification

The native probe checks preserved ring circulation, a ring standing wave
against its analytic solution, frozen energy, source work through a shared
outlet, isolated components, rollback, scale changes and allocation-free
steps. Over 1,000 steps on the same moving tree, its maximum
pressure-equivalent charge difference from `MovingTreeReference` is
5.33e-12 Pa.

The [independent receipt](moving_graph_reference.json) contains 24 seeded
graph cases with 20 steps each. These include cycles, parallel and ground
branches, disconnected components, moving nonuniform coefficients,
stagewise zero conductance and both kinds of source. The Python oracle solves
the full charge-and-momentum collocation system in normalized coordinates,
rather than the native pressure-only system. It carries its own state from
the initial condition throughout each case.

| Maximum discrepancy over 480 steps | Value |
|---|---:|
| Charge difference divided by second-stage compliance | 6.17e-13 Pa |
| Branch momentum | 8.33e-17 Pa s |
| Average pressure | 5.69e-13 Pa |
| Average flow | 3.26e-19 m³/s |
| Local stored-air balance | 1.54e-23 m³ |
| Global balance including boundary/source flow | 1.50e-23 m³ |

The receipt's `pressure_Pa` field means the first row: a pressure-equivalent
charge difference using `C2`. These arbitrary stage fixtures do not declare
endpoint compliance, so this is not an endpoint-pressure measurement.

A separate four-node shared-outlet fixture uses continuously varying
compliances and conductances. An independent DOP853 integration supplies the
continuous reference; tighter tolerances and a halved maximum step change
its endpoint by 1.67e-15 relative to the initial normalized state norm.

| Gauss subdivisions over 0.4 ms | Maximum relative state error | Refinement ratio |
|---|---:|---:|
| 32 | 1.18219e-4 | — |
| 64 | 7.47099e-6 | 15.824 |
| 128 | 4.68222e-7 | 15.956 |

This establishes fourth-order convergence on the declared smooth graph
fixture. Its coefficients are test inputs, not measured lip anatomy.

Fresh Release tree and graph CTests passed (2/2), together with the complete
independent graph audit. AddressSanitizer/UndefinedBehaviorSanitizer passed
the native graph checks and an eight-case discrete audit plus all three
continuous refinements. Its receipt is
`out/moving-graph-asan/receipt.json`. All 88 Release and 40 sanitizer recorded
source, binary and fixture artifact hashes were verified after the runs.
No performance benchmark was run.

## Lip and larynx connection still required

The existing planar export includes the entire 17.5 cm reference tract.
Appending the current articulator's lip and epilarynx sections would count
their resting lengths twice. For the present 64-section reference body,
the existing model allocates eight epilaryngeal sections and four lip
sections: 2.1875 cm and 1.09375 cm of resting length respectively. The planar
core must end at explicit replacement seams before those tubes are connected.
These section counts are existing discretization-rounded body fractions,
not new anatomical measurements.

The lip tube must use the final assembled section areas and lengths. Nominal
`lip_area` alone omits the separately applied lip constriction; protrusion
also changes section length. Complete closure introduces the zero-storage
problem described above. At the inlet, the actual fold aperture comes from
the two-mass dynamics and transglottal pressure feedback, not the planar cap
or a geometry-only replay.

The dense solver uses quadratic storage and cubic factorization work. It is
suited to the small coupled controls tested here, not a claim of full-grid
speech-rendering efficiency. Physical seams, closure-state transport,
fold/source coupling and radiation remain necessary before this path can
be assessed by listening.

## Reproduce

```sh
cmake --build build --target ducklips_moving_graph_probe ducklips_moving_tree_probe
ctest --test-dir build --output-on-failure -R 'ducklips_moving_graph_reference|ducklips_moving_tree_reference'
python tools/moving_graph_reference.py
```
