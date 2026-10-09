# Keep stored air with the continuing region

The native mean-pressure network can now rebuild its graph while retaining
volume, compliance and pressure for each continuing region. A separate
geometric reference establishes conservative continuation witnesses from the
moving contours. Together they preserve the contents of two chambers behind
a moving seal, even when region indices are deliberately reversed at every
pose.

This is progress toward moving-airway state integration. It does not yet
advance the complete Foxtrot pressure trajectory or change the speech
recording. Splits, merges, transient topology, missing correspondence and
acoustics remain unresolved integration requirements.

## Endpoint overlap can cross a sealed wall

A wall moving through a cell leaves some space formerly on its right on its
left. Intersecting old and new fluid fragments can therefore associate air
across a seal. Normalizing these overlap volumes into mass-transfer weights
would make a geometric intersection do the work of an unproved transport law.

The reference instead constructs a conservative enclosure of the moving
boundaries. For an edge with linearly moving endpoints, every swept point is
a convex combination of its four endpoint corners. Their convex hull contains
the edge's entire sweep. Subtracting all such hulls from an old/new fragment
intersection leaves space that remains clear throughout the declared motion.
A positive area supplies a geometric continuation witness.

The hull can overestimate the swept area of a rotating or deforming edge.
Consequently, the absence of a witness is inconclusive. Zero-area sweeps do
not remove positive-area witnesses and are omitted from the area subtraction;
intermediate point-in-fluid checks additionally verify the recorded witness
points. All Boolean geometry uses floating-point GEOS operations, not certified
exact predicates. No invalid polygon is automatically repaired.

This is geometric continuation, not tracking individual air parcels. Witness
volume is not a fraction of old mass to copy. Old contents swept out of a
region still require the conservative flow equation.

## A native operation for stored state

`CavityNetwork::reconfigure` accepts a new graph and an explicit old-node index
for each continuing node. It carries the old **V, C and p together**, so the
stored quantity `M = V + C p` survives indexing changes. A newly introduced
node starts with zero volume, compliance and pressure; its subsequent filling
must be handled by `advance`.

The operation rejects duplicate old-node references and omission of any old
node with nonzero volume or compliance. To remove a disappearing region, retain
its node through the flow step with zero final storage, expel its contents
through valid interval-open passages, and only then discard the empty slot.
Invalid graph indices and invalid maps leave the original graph, stored state
and last accepted port flows unchanged.

Reconfiguration allocates and is currently an offline research/control
operation. Ordinary pressure stepping remains allocation-free. Successful
reconfiguration clears the diagnostic port flows; this quasi-steady solver
has no stored port momentum. This API is not an acoustic-state remapper or
a completed streaming implementation. It does not infer geometry or authorize
an arbitrary split/merge weighting.

## Coupled moving-seal control

The geometric test uses two sealed chambers separated by a 0.1 mm wall in a
4 cm by 2 cm planar domain with 2.2 cm extrusion. The wall moves from zero
shift to +0.2 mm, then −0.1 mm, then back to zero. The two regions swap their
output ordering at every pose. Initial chamber pressures are 600 Pa and 0 Pa.
There is no interchamber passage.

Persistent-clear witnesses provide the index map. The native graph operation
preserves the old storage, and `advance` applies the new geometric volumes.
For each chamber, the independent solution of the declared linear storage
model is:

```
p_new = (M_initial - V_new) / (V_new / B),  B = rho c².
```

At 1, 2, 4 and 8 subdivisions per motion segment:

- Maximum pressure discrepancy is 1.60e−12 Pa.
- Maximum individual-chamber reference-volume discrepancy is 1.70e−21 m³.
- The pressures return to 600 Pa and 0 Pa, within 5.2e−13 Pa.
- Model density ratios remain between 0.99433 and 1.01013.

The coarsest motion contains an endpoint overlap of 0.088 mL across the moving
seal. The witness method excludes that cross-seal association. The test is a
controlled mean-pressure geometry, not a measured vocal tract, a moving-wall
energy validation, or rendered human speech.

Native tests also verify state preservation on failed reconfiguration, new
empty storage, rejection of duplicated/deleted air, conservative draining,
and removal only after the region is empty.

## Correspondence in the saved Foxtrot motion

The reference independently reconstructs every native region at all 521 saved
125 µs poses from 5.020 to 5.085 s. Maximum region-volume discrepancy remains
4.50e−21 m³. Across the 520 adjacent intervals it records:

| Relation | Count |
|---|---:|
| One-to-one witness pairs | 663,429 |
| All positive witness pairs | 663,467 |
| One-to-many split candidates | 10 |
| Many-to-one merge candidates | 9, in eight intervals |
| Old regions without a witness | 209 |
| New regions without a witness | 227 |
| Positive endpoint overlaps excluded by the sweep enclosure | 10 |

These counts are relations between local computational regions, not counts of
anatomical chambers or proven air-pocket births/deaths. Excluded endpoint
overlaps are unsupported by this conservative witness construction; they are
not all demonstrated physical cross-seal transfers.

All recorded witness points were independently tested in the fluid geometry
at one quarter, one half and three quarters of each interval: 1,990,401 checks.
The finite checks supplement the sweep construction; they do not by themselves
certify continuous-time topology.

Examples of unresolved events include:

- Cell 2056 splits between 5.030000 and 5.030125 s.
- Cell 3760 splits across the established closure bracket,
  5.031250–5.031375 s.
- Cells 3758 and 3631 each have a merge relation across the reopening bracket,
  5.062375–5.062500 s.

Repeating correspondence only at the original 1 ms endpoints finds two split
relations and no merge relations. It misses the other sampled intermediate
relations; endpoint-only state indexing is therefore inadequate.

## What a witness does not establish

Two explicit controls prevent overinterpreting this construction:

1. A small pocket can translate farther than its own width. It has a real
   continuation but no stationary-clear witness. Missing correspondence must
   trigger further geometric treatment, not deletion of its contents.
2. A shutter can cross a narrow joining corridor and leave it open at both
   endpoints. The large chambers retain clear witness points, so the endpoint
   relation is one-to-one, yet the middle of the interval has two disconnected
   fluid components. A one-to-one witness does not certify unchanged topology
   throughout a timestep.

Thus even the one-to-one Foxtrot relations are not a certificate that an entire
pressure step may be taken with one unchanged region graph. Internal contact
events and passage chronology still need to delimit the pressure evolution.
When disconnected pieces reconnect through a tiny internal neck, geometric
connectivity alone also does not justify instant pressure equalization.
The coupling across that neck must be resolved physically.

## Verification and reproducibility

The full Release suite passed 15/15. The cavity-network and mean-airway tests
passed ASan/UBSan, 2/2. All four moving-seal trajectories also passed under
ASan/UBSan with the independent geometric correspondence and analytical
pressure checks. Streaming was not rerun under sanitizers in this change.

- [planar_region_motion_reference.json](planar_region_motion_reference.json)
  records all 520 interval summaries, input hashes and full-correspondence
  output hashes. Detailed records are in `out/planar-region-motion/`.
- [cavity_region_motion_reference.json](cavity_region_motion_reference.json)
  records the four native coupled trajectories, source/binary hashes and
  analytical discrepancies. The sanitizer receipt is
  `out/cavity-region-motion-asan/receipt.json`.
- The coarser comparison is retained in
  `out/planar-region-motion-1ms/receipt.json`.

```sh
cmake --build build
ctest --test-dir build --output-on-failure
python tools/planar_region_motion_reference.py
python tools/cavity_region_motion_reference.py
```

The next integration must resolve split/merge events and unsupported
correspondence, consume [opening histories](PLANAR_FACE_MOTION.md) in order,
and preserve mass without assuming instantaneous mixing through narrow
connections. Physical flux laws, acoustic state and moving-boundary work,
transverse anatomy, tissue-volume conservation, nasal/side branches and
new human listening evidence remain requirements. No speech realism or
throughput gate has been declared satisfied. Performance work remains parked.
