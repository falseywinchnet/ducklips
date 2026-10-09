# Moving pressure through the Foxtrot contact intervals

The native reference now composes moving polygonal airspace, time-integrated
openings and conservative pressure storage. It advances both complete 1 ms
intervals containing the previously located tongue closure and release. The
two runs start independently; they do not form a continuous Foxtrot trajectory
and produce no new audio. The user's latest listening result remains:
November improved, Foxtrot still jammed around /tr/.

## Geometry and state carried through the step

The inputs are the actual saved full-mouth contours at 5.031–5.032 s and
5.062–5.063 s. Corresponding vertices move linearly between the original
endpoints. The native probe rebuilds volumes and openings from those contours;
it does not blend endpoint volumes or areas.

Each experiment inserts a vertical cut through its contact point and a
horizontal cut at y = 0.09 m. This creates five persistent pressure regions
and four possible connecting ports. The horizontal cut separates the upper
neck from a lower opening on the same vertical line. Without this separation,
a whole-face aperture would incorrectly combine openings belonging to
different region pairs.

Initial pressures are 600 Pa in the two rear regions and 0 Pa in the three
front regions. The outer boundary is sealed. There is no glottal source or
lip radiation. Every region retains its own stored air, including after
reopening; geometric connectivity does not cause pressure averaging.

`PlanarFaceMotion` detects zero-to-positive and positive-to-zero opening
changes. These times become pressure-step boundaries in addition to uniform
time subdivisions. Within each step, the native face integrator supplies the
average opening along the actual contour motion. The fixture uses

```
C = V / B, B = 139650 Pa
G = mobility * average_open_area / grid_centre_distance
mobility = 1e-7 m²/(Pa s)
M_new - M_old + dt * sum(outflow) = 0
M = V + C p
```

The mobility is a prescribed linear numerical control, not calibrated airway
friction. Interval-averaged conductance with implicit endpoint pressure is a
time discretization; it is not an exact solution for continuously varying
pressure and geometry.

The native tool requires a declared persistent graph. Counts and cell IDs
are checked internally, but do not prove geometric identity. The independent
audit establishes five matching persistent-clear witnesses at every step and
checks their interior points at quarter times. These witnesses do not certify
the absence of every possible unsampled internal topology event. General
region creation, disappearance, graph changes and complete contact scheduling
remain separate work.

## Independent checks and results

The release audit runs 1, 2, 4, 8, 16, 32, 64 and 128 uniform subdivisions for
each event, with the contact time inserted: 16 trajectories, 526 pressure
steps. At every stored state, GEOS independently checks region volumes,
region-pair openings and whole-fluid connectivity. Independent temporal
quadrature checks each interval aperture. A separate dense pressure solve
evolves from the initial pressures, without resetting to native intermediate
pressures.

| Maximum discrepancy across all release runs | Result |
|---|---:|
| Region volume against independent Boolean geometry | 4.14e-19 m³ |
| Endpoint port area | 5.21e-18 m² |
| Interval-average port area | 1.80e-13 m² |
| Pressure against independently evolved oracle | 1.14e-8 Pa |
| Whole-network stored-air balance | 2.72e-20 m³ |
| Individual-region stored-air balance | 6.42e-21 m³ |
| Sealed rear/front group stored-air balance | 6.78e-21 m³ |
| Flow through fully closed throat in checked intervals | exactly zero |

Closed-interval group checks exclude the contact-time roundoff neighbourhood;
they do not interpret a floating-point zero-gap classification as exact real
geometry. The geometry audit uses the actual native sampled contour coordinates
at contact. Reconstructing that pose with a differently rounded interpolation
can change the topology at the final bit.

Three independent quadrature attempts reported warnings. They were retained
in the audit files and retried with finer uniform integration partitions;
every final attempt completed without a warning. Independent quadrature does
not use the native integrator's internal event roots.

Time refinement reduces the final pressure differences once the coarse
transient is resolved. The 32-to-64 and 64-to-128 subdivision differences are
0.849 and 0.432 Pa for closure, and 5.305 and 2.657 Pa for release. The first
coarse differences are not monotonic. These are successive-resolution
differences, not certified error bounds.

## What this exposes, and what it does not establish

The release fixture reaches a minimum relative density of 0.880886, a roughly
12% reduction; closure reaches 0.975745. The release result persists under
time refinement. It exceeds what this linear storage control can establish
as a credible small-perturbation human-airway calculation. The sealed exterior,
coarse pressure partition, prescribed mobility, body motion and storage law
all remain physical assumptions to examine. No coefficient was adjusted to
hide this result.

The experiment establishes conservative numerical transport through these
moving contact cases. It does not validate the pressure magnitudes, spatial
resolution, physical neck-flow law, acoustic coupling or speech realism.
It does not change the production speech renderer or resolve the audible /tr/
problem. Performance optimization remains deferred.

## Validation and reproduction

```
cmake --build build --target ducklips_planar_pressure_motion_probe
python tools/planar_pressure_motion_reference.py
```

The input pose files come from the preceding planar geometry audit. The
[release receipt](planar_pressure_motion_reference.json) records their hashes,
the native binary and source hashes, all 16 results and artifact hashes.
Detailed native states and independent audits are in
`out/planar-pressure-motion/`.

Fresh AddressSanitizer/UndefinedBehaviorSanitizer checks passed four trajectories
(1 and 16 subdivisions per event, 38 pressure steps), including their full
independent audits. The receipt is
`out/planar-pressure-motion-asan/receipt.json`.

The shared correspondence checker now uses an explicit coordinate-roundoff
allowance proportional to boundary length, coordinate magnitude and depth
when checking whether a Boolean witness is a subset of its endpoint overlap.
This replaces an absolute volume tolerance that was inappropriate across
thin fragments and whole-mouth regions. It changes no geometry or witness
classification and is not a certified GEOS error bound. The original 520
intervals and controls passed again, with identical per-interval output
hashes, in `out/planar-pressure-witness-regression/receipt.json`.

Source, binary, input and per-run artifact hashes in both new pressure receipts
were verified after completion. The library implementation is unchanged in
this checkpoint. The prior 15/15 full release suite at `b358959` is historical
baseline evidence; it was not rerun for this new reference tool.
