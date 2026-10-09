# A closed stop must separate pressure storage

Historical component experiment: the subsequent opt-in speech integration is
documented in [SPATIAL_PRESSURE.md](SPATIAL_PRESSURE.md). The measurements below
remain the original two-chamber reference, not results from that integration.

The production renderer still assigns one mean pressure to the entire oral
volume. Its outlet area is the narrowest oral area (plus the nasal path).
When a tongue stop closes, motion anywhere in the mouth therefore changes
pressure on both sides of the closure. This gives front-chamber compression
an unphysical route to the larynx.

This change adds and tests a native `CavityNetwork` that stores pressure
separately in connected air spaces. **It is not yet integrated into speech.**
The existing audio still has the single-cavity defect. This report does not
claim a new pronunciation result or identify the cause of a heard click.
The current listening result is in [LISTENING_2026_10_09.md](LISTENING_2026_10_09.md).

## Evidence in the actual moving body

During the /k/ closure in Foxtrot, from 0.870 to 0.880 s, a fixed spatial cut
at oral section 42 separates the following changes:

- Rear air volume change: zero.
- Front air volume change: −9.984458765 mL.
- The neck stays sealed throughout the interval.

This is a geometry replay of the actual conservative-material body. Attachment
coordinates, attachment motion and tip error match the actual speech probe
at all 1,653 control steps to every printed digit. The original source trace
has cavity pressure reaching 1,674.8 Pa within that 10 ms interval, despite
lung pressure around 600–650 Pa. The wider /k/ interval reaches 1,769.9 Pa.
Pressure above lung pressure is possible under compression; the defect is
that front compression contributes to rear pressure through a sealed stop.

Conserving tongue tissue does not establish the accuracy of this air-volume
change. The current empirical midsagittal-gap-to-area conversion and omitted
sublingual/lateral spaces remain limitations. The following experiment tests
where a prescribed volume change acts, not whether 9.98 mL is physiological.

## Native pressure model

Each node stores its own volume V, compliance C and pressure p. Backward Euler
enforces the linearized mass balance

```
(C_new p_new − C_old p_old + V_new − V_old) / dt
    = prescribed_inflow − sum(outward passage flows).
```

Each passage obeys the quasi-steady passive pressure law

```
Δp = R Q + ρ Q |Q| / (2 A²).
```

Zero area seals a passage exactly. Open passages require positive R. Internal
connections form a tree, with any number of fixed-pressure reservoir ports.
Newton iteration uses the positive passage derivative
`dQ/dΔp = 1 / (R + ρ |Q| / A²)` and a tree solve, with backtracking. A failed
step leaves accepted pressures, volumes, compliances and flows unchanged.
`prepare` allocates; stepping does not allocate. This is a mean-pressure
component, not a replacement acoustic wave solver.

The mass residual tolerance is `1e−14 + 1e−11 max(abs(rhs))` m³/s. A sealed
connection transmits exactly zero flow. Disconnected pressures can differ
slightly between separate runs because the global convergence criterion can
stop Newton at different iterations; the independent native test measures
1.46e−8 Pa of such difference.

## Controlled comparison

The native geometry supplies volumes, neck area, front minimum area and a
one-cell neck length. The actual source trace supplies block-mean glottal
area; the block probe supplies lung pressure. The probe's glottal area lags
one block, so it is deliberately not used for that input.

For this reference only, wall compliance is split equally between two
chambers (5e−9 m³/Pa each), with air compliance added locally. The port's
linear resistance is circular Poiseuille resistance; the front passage uses
the neck's one-cell length as an explicit illustrative proxy. Neither that
length nor the wall-compliance split is a fit to human anatomy. Air constants
are the native defaults: ρ = 1.14 kg/m³, c = 350 m/s, μ = 1.86e−5 Pa s.

There are four runs: separate or lumped pressure storage, each with moving
or frozen front volume. All other inputs match. Initial separate pressures
are rear 600 Pa and front 0 Pa; the lumped pressure is 600 Pa. Therefore the
causal comparison is moving minus frozen **within each topology**, not an
absolute comparison of the two initial states.

| At 0.880 s | Separate chambers | Lumped chamber |
|---|---:|---:|
| Rear pressure, front moving | 625.603 Pa | 1,058.514 Pa |
| Rear pressure, front frozen | 625.603 Pa | 625.915 Pa |
| Rear rise caused by front motion | 0 Pa at emitted precision | 432.599 Pa |
| Mouth outflow, front moving | 782.115 mL/s | 0 mL/s |
| Flow across sealed neck | 0 mL/s | Outlet sealed |

Maximum mass residual across all four runs is 6.72e−14 m³/s. The lumped
experiment reproduces the topological mistake under the same passage laws;
it does not reproduce the full production pressure trace, which additionally
includes the coupled folds and acoustic feedback. These are prescribed-input
pressure experiments, not newly rendered speech.

Native checks additionally cover exact sealed compression, changing
compliance, independent series steady-flow balance, reverse flow, passive
equalization, branched topology, invalid-input rollback and no allocation.
The Release suite passed 11/11. After the final prepare-state invalidation
guard, the new component passed again in Release and ASan/UBSan (1/1).
No fresh full sanitizer suite is claimed.

## Remaining integration

The production path needs geometry-derived air spaces and passages that
remain conservative when closures form, move and reopen. It also needs:

- Glottal feedback from the rear pressure and the acoustic response.
- Frication driven by the local pressure drop across each constriction.
- Mouth and nasal mean flow routed to their own radiators.
- A consistent separation between mean and acoustic storage, so neither
  volume displacement nor glottal-flow fluctuations are counted twice.

Applying an orifice loss independently at every mesh cell is not an acceptable
shortcut: subdividing the same passage would add repeated nonlinear loss.
The diagnostic section 42 is not a production compartment rule or a Foxtrot
special case in the engine.

## Reproduce

```sh
cmake --build build
ctest --test-dir build --output-on-failure
./build/ducklips_say --canonical-motion --conservative-tongue \
  --body bodies/reference.body --dictionary samples/nato.dict --text 'Foxtrot.' \
  --chunk 48 --out out/foxtrot-source-audit.wav \
  --probe out/foxtrot-source-audit-probe.csv --trace 0.5 out/foxtrot-source-audit.csv
python3 tools/cavity_partition_reference.py \
  --probe out/foxtrot-source-audit-probe.csv --trace out/foxtrot-source-audit.csv
```

The stored receipt used the already-existing `foxtrot-local-material.csv`
probe. The independently traced audio is byte-identical to that prior speech
render. CSV inputs, outputs and native source/binary hashes are recorded in
[cavity_partition_reference.json](cavity_partition_reference.json). The full
enhanced Foxtrot fields from the preceding speech experiment remain valid;
there is no new audio or spectrogram attributed to this pressure component.
