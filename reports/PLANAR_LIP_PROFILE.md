# Preserve lip contact storage in the coupled mouth network

The native planar reference now consumes the integrated axial lip profile.
Its mouth regions and lip sections form one pressure/momentum network. The
lip sections retain their measured air volume while the exterior passage
closes, and accepted state continues across successive saved intervals.

This connects the [contact-storage geometry](LIP_CONTACT_STORAGE.md) to the
[planar mouth adapter](PLANAR_LIP_COUPLING.md). It is still a coarse inviscid
reference, not a production speech path or a new recording. The latest
listening evidence remains November improved and Foxtrot jammed around /tr/.

## Coupling and interval continuation

`ducklips_planar_inertive_motion_probe --lip-profile` accepts two endpoint
profiles exported from the actual articulator state. Each gives length,
base area, activation, place, core width, shoulder width and target area.
At each acoustic stage the adapter integrates section volumes and volume
centroids, then integrates inverse area along the paths between centroids.
It also integrates the inlet and outlet portions of the lip tube.

For each mouth-region attachment, the branch inertance is the sum of:

1. The core path to the seam, using the measured exposed cap area.
2. The profile's entry-path inertance, allocated by that cap's share of the
   complete mouth opening.

Both mouth regions connect to the same first lip storage node. The entry
path is counted once across those attachments. Lip narrowing and protrusion
affect both storage and inertia; a zero outlet opening does not erase the
inner shoulder's storage. Fully absent storage and unsupported interior lip
seals are rejected explicitly.

The profile input includes initial lip pressures. The optional
`--initial-state` suffix carries canonical node charge and branch momentum
into the next interval. The caller must certify identical boundary geometry
and node/branch ownership. This is continuation on an established graph,
not a general remeshing operation. The tested continuation reproduces every
accepted volume, compliance, conductance, charge and momentum exactly at
the shared interval boundary. Closed passages retain no latent momentum.

Activation crossings into or out of contact are inserted as integration
boundaries when the target aperture is zero. The tested lip place and flat
core extent are fixed. General moving-contact-boundary event scheduling is
not established by these controls.

## Actual Bravo geometry

The experiment uses six newly reconstructed native poses, retaining the
same body and articulator trajectory as the saved NATO replay. The core
partition ends at the previously established lip and epilarynx seams. A
coarse horizontal division gives two communicating mouth regions; both
meet the lip seam and feed the common four-section lip tube.

Two separate 2 ms controls each carry their state through two successive
1 ms intervals:

- Closing, 1.284–1.286 s: all six storage nodes initially at 0 Pa.
- Opening, 1.374–1.376 s: all six storage nodes initially at 700 Pa.

The inlet seam is sealed, and the lip exterior is at ideal ambient pressure.
There is no prescribed glottal source in these controls. They include the
actual moving core volumes as well as the lip motion; they are not an
isolated lip tube, nor a continuous simulation from the start of Bravo.

Reconstruction exposed a replay-export issue: an idle dorsal tracker stores
raw width zero, while the engine's effective gesture uses its existing
1e-4 m minimum width. The exporter now applies that same conversion before
reconstructing the surface, and reports the failing geometry operation if
export fails. No articulator or body parameter was retuned. Two historical
Foxtrot exports reproduce their previous geometry exactly after this fix;
six sanitizer Bravo exports match the Release exports exactly.

## Independent evidence

The [receipt](planar_lip_profile_reference.json) records the exact inputs,
source and binary hashes, moving geometry, carried state and refinements.
The full independent audit covers both intervals of both controls at 16,
32 and 64 subdivisions: 12 runs, 451 steps and 1,365 endpoint/stage geometries.

Independent Boolean geometry checks each core volume, internal aperture and
cap attachment. Independent axial quadrature supplies lip volumes,
centroids and path inertances. A full charge-and-momentum collocation solve
carries its own state across both intervals. It never resets itself to the
native solution at an interval boundary.

Geometric continuation witnesses preserve both core-region identities.
Persistent clear segments certify that the internal mouth connection stays
open throughout each complete interval, and conservative solid-sweep bounds
keep the mouth caps clear. The lip contact time is checked independently
against the activation crossing. Fully closed lip intervals pass exactly
zero outlet flow, and stored-air balance includes the actual outlet flux.

| Maximum discrepancy in the full independent audit | Value |
|---|---:|
| Pressure | 7.97e-9 Pa |
| Branch momentum | 2.88e-13 Pa s |
| Average branch flow | 3.26e-15 m³/s |
| Local stored-air balance | 6.63e-21 m³ |
| Per-interval global balance including outlet flow | 2.72e-20 m³ |

Checking the initial and final total air across both carried intervals,
together with their summed outlet flow, gives a maximum residual of
1.36e-20 m³ across the tested refinements. The mouth caps remain at least
19.55 mm from the conservative solid sweep in these intervals. This is a
geometric clearance certificate, not an anatomical calibration.

Further native-only runs check temporal refinement through 8,192
subdivisions per saved millisecond. They are separate from the independent
geometry/state audit. The short inviscid lip network has fast undamped modes,
so the coarse endpoint pressures are not monotonically converged. These
subdivisions are accuracy controls, not a proposed production rendering
rate or a demonstration of rendering efficiency.

The final 4,096-to-8,192 subdivision comparison changes the endpoint pressure
by 0.009450 Pa in the closing control and 0.001405 Pa in the opening control.
These are the maximum differences across all six storage nodes after the
complete two-interval continuation.

For fine diagnostics, `--compact-geometry` omits repeated contour coordinates
from the log. It changes no geometry calculation or acoustic step. A paired
run verifies exact equality of all retained metrics and states with the full
log. The saved input still contains both complete endpoint contours.

A separate static test closes the entire exterior of the coupled mouth/lip
network and initializes it uniformly at 700 Pa. It holds that pressure with
maximum error below 9e-11 Pa and zero lip leakage. The previous three uniform-
lip controls and two sealed tree controls reproduce their historical output
exactly, preserving their earlier interpretation.

AddressSanitizer/UndefinedBehaviorSanitizer pass all four moving intervals
at 16 subdivisions, including continuation. Four focused native CTests pass
for planar airspace, face motion, the graph solver and the lip profile.
All 148 Release and 40 sanitizer recorded source, binary, input and fixture
hashes were verified after the runs.

## Reproduction and remaining speech integration

```sh
cmake --build build --target ducklips_airspace_geometry_probe ducklips_planar_inertive_motion_probe
mkdir -p out/planar-lip-profile/poses
build/ducklips_airspace_geometry_probe out/nato-airspace-geometry.csv 1.284 1.286 512 out/planar-lip-profile/closure-geometry.csv --frames out/planar-lip-profile/poses --core bodies/reference.body
build/ducklips_airspace_geometry_probe out/nato-airspace-geometry.csv 1.374 1.376 512 out/planar-lip-profile/opening-geometry.csv --frames out/planar-lip-profile/poses --core bodies/reference.body
python tools/planar_lip_profile_reference.py
```

The profile CSV comes from the native replay described in
[LIP_CONTACT_STORAGE.md](LIP_CONTACT_STORAGE.md). Detailed controls and
regression artifacts are under `out/planar-lip-profile/`; the sanitizer
receipt is `out/planar-lip-profile-asan/receipt.json`.

The native adapter now handles the tested complete terminal lip closure
without a compliance floor or discarded air. It still needs the epilarynx
and dynamic vocal-fold source, physical losses and radiation, and continuous
state transport through the full utterance's changing region topology.
Finer partitions with genuinely disappearing pressure cells remain outside
this fixed-graph continuation. Numerical agreement here does not establish
human realism, pronunciation or the requested speed. Performance optimization
remains parked.
