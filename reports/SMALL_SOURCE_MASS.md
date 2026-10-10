# Resolving force derivatives as the air support shrinks

The sustained-source failure exposed a force-derivative error in a partly
closed contour. The source stores two weighted reference-air moments. As
contact pushes the remaining air toward one end, one of those moments can
become much smaller than the nominal 1e-8 m3 coordinate scale.

The previous kinetic derivative used

```
h = min(2e-5 * max(1, abs(z)), 0.1 * abs(z))
```

for each dimensionless state coordinate. For the first gas moment at the
last emitted state before the fixed-tube failure, this was a 10% perturbation.
Halving the requested finite-difference factor twice still hit that cap:
those ostensibly finer calculations did not refine that derivative at all.

Positive gas moments now use `h = 2e-5 * abs(z)`. Signed opening areas retain
their existing contact-side treatment. This changes the calculation of the
existing Hamiltonian force; it adds no source waveform, leakage, gap floor,
air reset, damping term or pressure adjustment.

## Independent force check

The [reference calculation](small_source_mass_reference.json) uses separate
quadrature, a nonlinear solve for the prescribed gas moments, and the
fixed-velocity fluid metric identity

```
d kinetic_energy / d gas_moment at fixed momentum = -vᵀ (dG/d gas_moment) v / 2
```

It compares that result with the native inverse-metric derivative at fixed
momentum. The fixture uses the observed gas flows and held walls; it does not
invent the tissue velocities absent from the failed trace. Two independent
relative perturbations check the reference derivative.

The old first-moment force differs by −3.52979 Pa. With relative moment
perturbations, the discrepancy is below 0.001 Pa over the checked refinements
(about 1.37e-8 Pa for the default factor). That establishes the derivative
defect at this configuration. It does not establish that the sustained-source
failure has been repaired.

The 606 kPa reported source effort is a profile coefficient extrapolated to
the closed end, not the pressure throughout the actual air pocket. The
occupied support is approximately 0.09915 mm long. Reconstructing the field
there gives physical endpoint pressures of +17.219 and −3.862 kPa. Small
weighted moments and large extrapolated coefficients can therefore coexist
with much smaller pressures in the air. They require numerical care rather
than interpreting the coefficient alone as a physical pressure spike.

## Trajectory evidence and remaining work

Six Release and six ASAN/UBSAN controls pass after the change: continuous
source, joint third/fifth-order passive coupling, coupled accuracy, and both
exact-zero source-opening controls. Logs and executable/source identities
are bound by the [validation receipt](small_source_mass_validation.json).

The old failing fixed-tube path is being replayed with its original physical
code and time grid to save the full state at 28.95 ms, immediately before
the reported failure. The new `joint_contact_probe --stop-after N` option
stops at an existing observation index without changing earlier timesteps.
The completed prefix matches the preserved old observations exactly.

A separately frozen corrected executable is queued to continue that state
for 50 microseconds. This is a targeted check of the actual failure location,
not a fresh onset with different initial conditions. Both jobs retain their
native files under `out/source-closure-checkpoint/` and
`out/source-small-mass-retry/`. The retry result remains pending in this
record. No sustained phonation, speech render or intelligibility improvement
is claimed yet.
