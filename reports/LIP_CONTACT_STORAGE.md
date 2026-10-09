# Separate lip storage from the contact aperture

The lip closure investigation found a geometric cause of the apparent
zero-storage problem. In the reference body, the final section is 2.734375 mm
long, while the lip gesture's fully closed flat core extends 1.5 mm inward
from the exterior. Its midpoint lies inside contact. Assigning that midpoint
area to the entire section therefore erases the still-open inner shoulder.

The new `LipProfile` integrates the existing axial constriction law to obtain
air volume, volume centroid and air-column inertia separately. At complete
exterior contact, the last section retains its inner air volume and remains
connected to the preceding section; its outlet conductance is exactly zero.
No minimum area, minimum volume, pressure reset or new lip anatomy is used.

This is implemented as reusable native geometry plus a coupled lip-network
reference. It has not yet replaced the planar adapter's piecewise-uniform
lip input or the production speech network. No new speech recording or
listening improvement is claimed.

## Integrate the field already specified by the articulator

`TractModel::lip_profile(state)` exports the current physical tube length,
base area, activation, place, flat-core half-width, Gaussian shoulder width
and target area. The body supplies its existing length and area scales;
protrusion supplies additional length. These are the same parameters used
by the existing lip-slot rule. Its small-weight cutoff is retained.

Let `x` measure distance inward from the lip opening. With base area `A0`,
target area `At`, activation `a` clamped to `[0,1]`, core half-width `h`,
place `p` and shoulder width `w`, the minimum area is

```
Amin = min(A0, At) + [A0 - min(A0, At)] (1-a)^2 (1+2a)
u    = max(|x-p|-h, 0) / w
A(x) = Amin + (A0-Amin) [1-exp(-u²)]
```

Outside the existing `exp(-u²) >= 1e-4` shoulder support, the rule returns
the base area. The factored complement of smoothstep preserves tiny positive
openings near `a=1`; subtracting a rounded smoothstep from one can invent
premature contact. `-expm1(-u²)` similarly preserves small shoulder openings.

For an axial storage interval:

```
V = integral A(x) dx
x_centre = integral x A(x) dx / V
C = V / (rho c²)
G = 1 / [rho integral_path dx / A(x)]
```

Branch paths run between these volume centroids, or from a centroid to a tube
end. A path intersecting contact has `G=0`. Its neighboring storage need not
be zero. Quadrature splits at the known core and shoulder boundaries and
checks error adaptively. Failure is explicit; there is no positive-aperture
threshold or fallback metric.
Shoulder integration uses local normalized coordinates. Subtracting the
contact location repeatedly in global coordinates loses relative precision
at tiny gaps and can prevent the inverse-area integral from meeting its
error check. The local-coordinate calculation retains a positive outlet
conductance even in the tested `activation = 1 - 1e-9` case.

For the completely closed, zero-target reference lip, the remaining volume
in its terminal interval has an independent analytic expression:

```
u = (section_length - core_length) / shoulder_width
V = A0 shoulder_width [u - sqrt(pi)/2 erf(u)]
```

At the saved Bravo pose at 1.285 s, that volume is 9.081920 microlitres. The
midpoint area is zero. These statements describe different geometric
quantities and are both true.

## Whole NATO geometry evidence

The extended native replay exports the profile through `--lip-profile CSV`.
Its assembled airway CSV and tongue-geometry CSV reproduce the previous
23,239-frame replay byte for byte. All four lip midpoint areas computed from
the exported field match the existing assembled sections to within
2.17e-19 m².

Across the alphabet, 562 frames have a zero final midpoint area. Every one
retains positive final-section air volume, a positive connection to the
preceding section, and a sealed exterior. The smallest retained terminal
volume is 8.199331 microlitres. Thus the disappearance of the entire terminal
volume in these frames was a midpoint representation artifact under the
existing law, not a required consequence of lip closure.

Independent SciPy quadrature checks selected poses, including closure and
release boundaries. Dividing the same closed profile into 1, 2, 4, 8, 16 and
32 sections preserves its total air volume. The finer partitions correctly
contain wholly empty terminal intervals. Geometry can represent them;
the positive-storage acoustic probe rejects those intervals pending state
transport rather than giving them a small compliance.

## Coupled closure and reopening controls

The native lip-network probe connects four integrated storage nodes with
their integrated inertances. The inlet is a declared ideal pressure boundary;
the exterior is at 0 Pa. Wall-displacement flow is the negative volume change
divided by each step duration. The graph solver evolves charge and momentum.

Two separate 2 ms controls use the actual saved profile parameters:

- Closing: 1.284–1.286 s, initially at 0 Pa, inlet at 0 Pa.
- Opening: 1.374–1.376 s, initially at 700 Pa, inlet at 700 Pa.

These controls are separately initialized, not a continuous replay of Bravo
from the start of the utterance. Parameters are interpolated between native
control frames and the area law is evaluated at the acoustic stages. The
closing step is split where activation reaches contact: fraction
0.4987278312279807 of its first millisecond under this interpolation.
Already closed passages do not accumulate latent branch momentum.

An independent full-state collocation solve uses independently integrated
metrics and carries its own volumes and state forward. Seven audited runs
(three subdivisions for each moving control and one static closed control)
check 1,545 stage/endpoint geometries. All tested fully closed intervals have
exactly zero outlet flow. Stored-air balance includes the actual inlet and
outlet fluxes.

| Maximum discrepancy in the independent dynamic audit | Value |
|---|---:|
| Pressure | 1.07e-8 Pa |
| Branch momentum | 2.93e-14 Pa s |
| Average branch flow | 5.42e-16 m³/s |
| Local stored-air balance | 2.10e-22 m³ |
| Global balance including external flow | 8.48e-22 m³ |

The static closed control starts at 700 Pa with the inlet also at 700 Pa.
All four pressures remain within 1e-8 Pa of that value, with positive storage
and no outlet leakage. Closure does not erase the stored air or force its
pressure to zero.

The isolated short-tube reference has fast undamped modes. Coarse time
stepping is not converged even when its linear solve is accurate. Further
native-only refinement is kept separate from the independent geometry/state
audit. The last 4,096-to-8,192 subdivision comparison changes endpoint pressure
by less than 0.01 Pa in closing and less than 0.002 Pa in opening. This is
numerical evidence for the declared inviscid control, not a required
production oversampling rate or a demonstration of realistic speech.

The [receipt](lip_profile_reference.json) records exact errors, commands,
source/binary/input hashes, the complete geometry run, the independent
checks and the native refinement results. AddressSanitizer and
UndefinedBehaviorSanitizer pass selected geometry frames, both moving
controls at 16 subdivisions, and the static closed control; their receipt is
`out/lip-profile-asan/receipt.json`. Fresh focused native CTests pass for the
core engine, graph solver and lip profile (3/3).
All 86 Release and 32 sanitizer recorded file hashes were verified after
the runs. Additional direct CLI checks confirm rejection of wholly empty
pressure sections and interior seals, and positive outlet conductance at
activation `1 - 1e-9`; their local record is
`out/lip-profile/edge-case-checks.json`.

## Integration limits and reproduction

This corrects the representation of the existing lip law; it does not
validate that law against measured lip tissue. Viscous effects, radiation,
mean-flow losses and tissue-force feedback are absent from these controls.
The standalone probe supports the connected remaining air behind terminal
contact. An interior seal that would split storage into separate regions is
rejected, as is any wholly empty pressure section. General region birth,
disappearance and changing partitions still require conservative transport.

The next integration step is to replace the planar lip adapter's midpoint
volumes and half-section inertances with this profile's separate metrics,
then couple source and losses on the resulting geometry. The current speech
path and its latest listening assessment remain unchanged.

```sh
cmake --build build --target ducklips_tongue_reach_audit ducklips_lip_profile_probe
mkdir -p out/lip-profile
build/ducklips_tongue_reach_audit bodies/reference.body samples/nato.dict 'Alfa. Bravo. Charlie. Delta. Echo. Foxtrot. Golf. Hotel. India. Juliett. Kilo. Lima. Mike. November. Oscar. Papa. Quebec. Romeo. Sierra. Tango. Uniform. Victor. Whiskey. X-ray. Yankee. Zulu.' --lip-profile out/lip-profile/nato.csv --airway out/lip-profile/airway.csv > out/lip-profile/geometry.csv
build/ducklips_lip_profile_probe
python tools/lip_profile_reference.py
```
