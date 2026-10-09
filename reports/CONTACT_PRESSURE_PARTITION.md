# Retain pressure states across a small opening

The two principal contact events in the saved Foxtrot interval are now located
geometrically. A pressure partition through each contact retains the narrowing
opening, instead of collapsing both sides into one pressure as soon as their
air becomes geometrically connected. A native zero-time partition operation
introduces that refinement without a pressure pulse or artificial mixing.

The experiment also exposed and fixed a native conservation stopping error:
a fixed absolute residual threshold could report flow through a tiny opening
while leaving both chambers' stored air unchanged.

These are local geometry/pressure controls using actual saved contours. They
are not the complete moving Foxtrot trajectory or new rendered speech. The
audible /tr/ problem remains unresolved.

## Located contact and release

Vertex-edge contact equations are solved within the previously identified
cells and time brackets. Corresponding vertices follow the same declared
piecewise linear paths used by the prior motion audit. Independent GEOS
component-count bisection checks the selected event.

| Event | Saved-pose time | Contact point, metres | Native contour identities |
|---|---:|---|---|
| Closure, cell 3760 | 5.0313321308 s | (0.07365453, 0.09543106) | Pad vertex 125, domain edge 234 |
| Release, cell 3758 | 5.0624035949 s | (0.07003643, 0.09607271) | Pad vertex 105, domain edge 228 |

The whole fluid domain changes from one component to two at closure and from
two to one at release. Algebraic and independent bisection fractions agree
within 1.48e−14 for these two events. This precision concerns the supplied
polygon model; it is not the accuracy of human anatomy, articulator motion,
or recorded speech timing.

The candidate calculation includes vertex-edge contacts and the local cell
boundary. It is a selected-event reference, not a complete continuous-time
scheduler for arbitrary polygon arrangements. Other event types and all
remaining Foxtrot split/merge relations still require integration.

## Keep a partition through the contact

For each event, clip the actual mouth domain to its affected original grid
cell and insert a vertical cut through the contact point. The native planar
geometry and independent Boolean geometry then measure the two pressure
volumes and their connecting opening. The artificial outer boundary is sealed
for the local pressure control.

On the open side of both events, the unrefined cell contains one connected
region. The refined cell retains two volumes joined by an opening whose area
tends to zero at contact. On the closed side the two volumes have no connecting
port. The cut should be introduced before the neck becomes restrictive and
retained through reopening; refinement cannot reconstruct a pressure gradient
already lost by an earlier coarse representation.

The frozen-geometry control starts the two local volumes at 600 Pa and 0 Pa.
It uses a declared linear mobility, `G = mobility * A / distance`, with
mobility 1e−7 m²/(Pa s), and advances 0.1 ms. This is a numerical flow fixture,
not calibrated vocal-tract friction. It tests the contact limit, not the
correct human pressure magnitude or a moving-wall trajectory.

As the sampled geometry approaches contact by factors of ten:

| Case | Opening area | Receiving pressure after the frozen step | Pressure if immediately combined |
|---|---:|---:|---:|
| Closure, first sample | 2.04717e−9 m² | 11.1914 Pa | 544.176 Pa |
| Closure, smallest opening | 2.04720e−14 m² | 0.000114870 Pa | 544.347 Pa |
| Release, first sample | 5.25030e−9 m² | 9.70218 Pa | 403.845 Pa |
| Release, smallest opening | 5.25032e−14 m² | 0.000100627 Pa | 405.282 Pa |

The smallest offsets are mathematical limiting tests of the polygon model,
not a claim that continuum human-airway physics is resolved at those scales.
Both sealed controls retain the initial 600 Pa / 0 Pa pressures exactly.

The instantaneous combined pressures follow directly from assigning one
uniform pressure to the cell while preserving its total stored air. They are
not measurements from a new speech recording. Their finite jump as opening
area tends to zero demonstrates why connectivity alone is insufficient for
choosing a pressure degree of freedom.

## Partition storage without advancing time

`CavityNetwork::partition` receives actual child volumes and compliances at
the same geometry and time as the parent. Children retain their parent's
pressure; their volumes and compliances must each sum to the parent's stored
values within 64 machine epsilons, relative, with no absolute volume floor.
The caller establishes that the children are an actual disjoint geometric
partition. Matching sums alone cannot establish that geometry.

For the declared linear storage model, this preserves both

```
M = V + C p                 stored air in reference-volume units
E = 1/2 C p²               quadratic compression energy
```

within the stated roundoff allowance. Different old parents cannot be merged
or pressure-averaged by this API. Unparented nodes must start empty. Physical
volume changes still go through the conservative time-step equation.
Invalid partitions leave the original state intact. Like graph reconfiguration,
partitioning currently allocates and is an offline/control operation; ordinary
pressure stepping remains allocation-free. Port momentum and acoustic state
are outside this mean-pressure operation.

Across the twelve actual-contact refinements, the pressure change is exactly
zero. Maximum reference-volume discrepancy is 1.70e−23 m³ and maximum
compression-energy discrepancy is 2.16e−23 J. These small nonzero discrepancies
reflect independently integrated parent/child geometry within the explicit
roundoff allowance; they are not hidden physical wall motion.

## Fix conservation accuracy for tiny regions

The initial smallest-opening control failed for a useful reason. The native
solver returned a flow of about 1.43e−15 m³/s but left pressures at 600 Pa and
0 Pa. Its old absolute convergence floor, 1e−14 m³/s, accepted the local
imbalance before the stored air was updated. Equal and opposite internal flows
still cancel in global mass, so checking only total mass would miss this.

Convergence now uses each node equation's own scale:

```
scale_i = |(C_new/dt) p_i| + |rhs_i| + sum |incident flows|
|residual_i| <= 1e−11 scale_i + explicit floating-point roundoff allowance
```

The separate whole-network balance retains its own scale and roundoff budget;
internal flow cancellation does not excuse a boundary imbalance. The inner
linear solve also uses a relative residual target without a dimensional floor.
The existing pressure-representation allowances remain in place.

The regression scales volumes, compliance and conductance from 10⁻⁶ to 10⁶
of the fixture size while retaining its pressure solution. Maximum relative
local mass-balance error is 1.67e−16. The actual-contact controls agree with the
independent two-pressure analytical solution within 1.14e−13 Pa, with maximum
local reference-volume balance error 1.51e−27 m³.

This change affects the native mean-pressure solver used by the experimental
speech path. No new speech has yet been rendered or assessed after it.

## Validation and evidence

- Full Release suite: 15/15 passed.
- Cavity-network and mean-airway tests under ASan/UBSan: 2/2 passed.
- Twelve open-neck cases, two sealed cases and twelve zero-time partitions
  passed in Release and under ASan/UBSan. All 26 local native geometries were
  independently checked against GEOS. Maximum refined-region volume
  discrepancy was 8.34e−23 m³.
- Existing independent regressions passed: 24 nonlinear graph fixtures and
  four fixed Foxtrot graphs. Maximum pressure differences were 2.65e−9 Pa
  and 3.15e−8 Pa respectively. These retain their previously declared numerical
  flow laws; they do not validate human viscosity or moving speech.

[contact_partition_reference.json](contact_partition_reference.json) records
contact identities, independent brackets, geometry/pressure/partition checks,
commands and hashes. Detailed outputs are in `out/contact-partition/`.
Sanitizer evidence is in `out/contact-partition-asan/receipt.json`; the broader
pressure regression receipt and exact driver are in
`out/contact-partition-regression/`. Recorded source, binary, input and output
hashes were checked against final files. Streaming was not rerun under
sanitizers in this change.

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/planar_contact_partition_reference.py
```

## Remaining integration

The pressure mesh can now be refined conservatively and retain a finite
connection through contact/reopening. The complete moving-region evolution
still needs an event scheduler, disjoint spatial partitions and state
correspondence across all relevant events, interval-consistent flow, and
verification of wall work. A single pair of contact cuts is not a complete
mesh for the whole articulated trajectory.

The eventual speech path must also retain the nasal and side branches,
glottal/radiation coupling and acoustic state. Transverse anatomy, tissue-volume
conservation, physical flux laws and new human listening remain unverified.
The latest recording is still the assessed spatial-pressure alphabet.
Performance work remains parked; no realism or throughput gate is declared
satisfied by these controls.
