# Preserve boundary openings and their owning air regions

The planar geometry now retains selected domain boundary edges as explicit
air connections. The moving inertive reference can attach an ideal external
pressure to each resulting region/edge pair. This removes the forced sealed
exterior from the controlled experiment without inventing an internal route
through the tongue.

Opening the mouth-end cap to ambient pressure reduces the release experiment's
worst sampled density reduction from 3.7% to 0.65%. It confirms that the sealed
front-mouth boundary was a major driver of the preceding pressure drop.
These are pressure/flow controls, not a new speech recording.

## Geometric boundary ports

`PlanarAirspace::build_regions` and `build_between` accept an optional span of
unique domain-edge indices. `boundary_ports` returns the edge identity, owning
air region, exposed area, area centroid and outward unit normal. The default
empty selection preserves the existing closed-boundary behavior.

The extraction uses the same polygon strips and region identities as the
interior volume calculation. A sloping edge keeps its true length. A vertical
edge uses the strip on the domain's interior side. Solid intervals on the
boundary's supporting line are subtracted, including collinear tissue on the
outside of the domain. Thus exterior tissue contact can seal a cap even when
the interior cell still contains air. No positive-aperture threshold is used.

Different air regions meeting the same cap retain separate attachments.
Disjoint openings belonging to the same edge/region pair are aggregated by
area and first moment. Their centroid may lie in the tissue between them; it
is not a point guaranteed to be inside an opening. This aggregation does not
establish a single physical jet or a common radiation impedance for those
separate openings.

Output and scratch storage are reserved during preparation. Boundary
extraction and moving-boundary sampling add no allocations after that step.
Capacity failures, invalid edge indices and duplicate edge selections fail
explicitly. An ordinary geometry build clears old boundary attachments.

## End-cap provenance and the anatomical distinction

The native pose exporter now records two edge identities:

- `supraglottal_entry`: the closing edge between the inner and outer walls at
  the start of the reconstructed tongue/floor domain.
- `mouth_end`: the cap joining the outer and inner walls at its oral endpoint.

Across all 66 saved Foxtrot poses, the new exports match every historical
geometry field exactly; only these labels are added. The cap areas are about
264.55 mm² and 428.085 mm² respectively under the fixed 22 mm extrusion.

The cap labels are not a claim that these are complete anatomical openings.
The midsagittal reconstruction currently includes jaw/root and tongue geometry,
but not the separately controlled lip narrowing or protrusion used elsewhere
in the speech model. The supraglottal entry also is not the time-varying space
between the vocal folds. Attaching a full glottis or lip model requires those
remaining geometries and states; naming the caps does not supply them.

## Open-boundary pressure controls

The moving inertive probe's optional `--boundary-loads` input declares each
domain-edge/region attachment, an inertive path length, and a held external
pressure. It appends a branch to the external pressure boundary:

```
G_boundary = measured_open_area / (rho * declared_length)
Q_outward = G_boundary * lambda
lambda' = p_region - p_external
```

The control uses the normal distance between the coarse grid centre and the
cap as its path length. It does not add an empirical resistance or clamp the
interior region's pressure. Stored-air balance includes the actual integrated
external flow. The existing canonical series-pressure input supplies the
external pressure with the correct sign.

The two selected 1 ms intervals run separately with the same initial pressures
as the previous sealed control: 600 Pa in the rear regions and 0 Pa in front.
Two boundary configurations are compared:

1. Mouth-end cap open to 0 Pa, with the inlet still sealed.
2. Mouth-end cap at 0 Pa and the supraglottal entry at 600 Pa.

These are ideal pressure boundaries with air-column inertia. They have no
radiation loss or external acoustic field, and the prescribed inlet pressure
is not an oscillating vocal-fold source. The resulting pressure oscillations
must not be interpreted as a calibrated mouth response.

| Release result at 128 subdivisions | Sealed exterior | Mouth cap open | Both caps connected |
|---|---:|---:|---:|
| Worst sampled relative density | 0.962888 | 0.993450 | 0.993450 |
| Final lower-front-left pressure | −5182.8 Pa | +718.2 Pa | +716.5 Pa |
| Final lower-front-right pressure | −4646.6 Pa | +273.3 Pa | +271.5 Pa |
| Final upper-front-right pressure | −4681.9 Pa | +390.4 Pa | +382.1 Pa |
| Net ambient air entering the mouth cap | 0 | 0.55551 mL | 0.55939 mL |

Mass is expressed in reference-volume units under the linear storage model.
Positive final pressures do not mean that no negative transient occurred:
the quoted minimum density covers the sampled trajectory. The two open-cap
cases reach nearly the same minimum because their early front-mouth response
is isolated from the inlet by the tongue closure. Flow through that closed
contact remains exactly zero in the checked closed intervals.

The probe schedules internal contact events. It does not yet schedule arbitrary
boundary occlusion events. For these particular end caps, independent convex
enclosures of every moving solid-edge sweep, together with the endpoint solid
interiors, prove positive separation from the complete moving cap. Minimum
clearance is 16.49 mm at the mouth cap and 32.94 mm at the inlet over the two
cases. The tested caps therefore remain unobstructed; this does not license
using the same time stepping through an unexamined moving lip closure.

## Independent evidence

The [boundary geometry receipt](planar_boundary_reference.json) records 66
poses and 1,386 boundary region attachments. Independent Boolean geometry
checks every interior region and port, then intersects each selected cap with
its owning region while subtracting the solid union. Maximum area discrepancy
is 6.78e-20 m², centroid discrepancy 1.39e-17 m, and normal discrepancy zero.
Native controls additionally cover separate chambers on one cap, reversed
winding, sloping caps, exterior partial/full contact, a 1e-12 positive opening,
moving cap identity, invalid selections, and allocation-free execution.

The [open-boundary receipt](planar_open_boundary_reference.json) records sixteen
moving trajectories, 976 integration steps and 2,944 independently checked
endpoint/stage geometries. A separate dense collocation solve evolves from
the initial state using independent geometric metrics, with external pressure
and flow included. Its maximum endpoint-pressure discrepancy is 3.01e-9 Pa.
Maximum local stored-air discrepancy is 6.64e-21 m³; whole-domain balance,
including boundary flux, is within 2.72e-20 m³. In the receipt, the historical
`sealed_group_mass_error_m3` field also includes external flux when
`group_balance_includes_external_ports` is true: the tongue-sealed groups may
exchange air through their own external caps.

The maximum 64-to-128 subdivision endpoint difference is 0.227 Pa. This is
time refinement of the declared coarse equations, not a spatial or physical
validation. Both original sealed 128-subdivision trajectories are rerun and
reproduce every saved charge, momentum, pressure, flow and average flow
exactly, despite the added optional boundary support.

Analytical controls check a closed-inlet/open-outlet quarter-wave duct at
three spatial and three temporal resolutions, and a single compliant region
driven through an inertive cap by positive external pressure. The latter
checks the exact forced-oscillator pressure, inward flow direction and source
work against stored energy. An ideal pressure-release boundary preserves the
undamped model's energy; it must not be confused with radiative absorption.

Fresh native geometry tests and independent audits passed in Release and
under AddressSanitizer/UndefinedBehaviorSanitizer. The sanitizer geometry
audit covers seven saved poses; its dynamic audit covers all twelve analytical
controls and all four event/boundary combinations at 16 subdivisions. Their
receipts are `out/planar-boundary-asan/receipt.json` and
`out/planar-open-boundary-asan/receipt.json`.

## Reproduction and remaining integration

```
cmake --build build --target ducklips_planar_airspace_tests ducklips_planar_regions_probe ducklips_airspace_geometry_probe ducklips_planar_inertive_motion_probe
build/ducklips_airspace_geometry_probe out/nato-airspace-geometry.csv 5.020 5.085 512 out/planar-boundary/geometry.csv --frames out/planar-boundary/poses
build/ducklips_planar_airspace_tests
python tools/planar_boundary_reference.py
python tools/planar_open_boundary_reference.py
```

The boundary, input, source, binary and per-run output hashes are preserved in
the receipts. Historical artifacts are retained separately from these runs.
The next integration must supply actual lip and fold geometry and a physical
radiation/source load. Viscous and jet losses, adequate spatial inertance,
nasal paths and continuous full-utterance state transport remain unfinished.
The latest listening assessment still belongs to the spatial-pressure
recording: November improved, Foxtrot's /tr/ still jammed. No new pronunciation
or realism result, spectrogram, benchmark or full-suite pass is asserted here.
