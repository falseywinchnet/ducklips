# An open lip outlet can end inside its flat core

Connecting the current pressure-driven source to an open actual-mouth pose
exposed a rejected valid geometry. The outlet-loss quadrature assumed the last
lip cell's volume centroid lay beyond the profile's flat core. That is true
for the preceding near-seal controls, but not for the wider opening at 1.500 s
in the saved articulation trace. Its outlet path ends near 1.371 mm, inside
the 1.500 mm core. The earlier executable exits during initialization when
outlet viscosity is enabled.

The outlet path now handles both cases with its existing state layout:

- Beyond the core, it retains the exact core element and the existing shoulder
  quadrature. The previous closing continuation produces identical full JSON.
- Inside the core, the whole path is a constant slit of length equal to the
  centroid distance. Its air inertia is exactly `rho * length / area`; the
  shoulder has zero length and contributes no inertia or loss memory.

This required distinguishing an absent *local loss element* from a closed
*whole flow branch*. `SlitViscous` now accepts zero local inertia and zero gap
on an open branch. It keeps the branch conductance, contributes no drag, and
leaves branch momentum unchanged. Retiring an element clears its local
memories and records their energy change as geometry work. Reappearance
starts without resurrecting those memories. An element with zero inertia and
a positive gap remains invalid. This extends the existing moving-metric
model; it does not establish a full moving-wall Stokes/advection solution.

## Evidence

The [independent geometry receipt](lip_outlet_core_reference.json) checks six
actual lip profiles at four quadrature orders, including the earlier sealing
and reopening poses and the two new open poses. Maximum relative geometry
discrepancy is 3.44e-13. For the core-only paths, the complete distributed
impedance reduces to the single constant-slit element, with eight absent
shoulder slots at the production reference's order.

Five focused Release and five ASAN/UBSAN checks cover the lip quadrature,
local viscous profile, joint solve and accuracy controller. A 10 µs actual
open-mouth integration with all existing thermal, viscous and radiation
connections passes in both builds. The historical closing interval remains
exactly equal to its recorded complete JSON. These checks are collected in
the [validation receipt](lip_outlet_core_validation.json).

## Actual articulation with the continuous source

`tools/actual_coupled_phonation.py` now drives the existing actual-geometry
adapter with a declared lung pressure. The reference helper's default remains
700 Pa, including its historical input bytes. The new run uses 1400 Pa and
the body's source parameters, with the current continuous gas/tissue engine
and joint thermal, viscous and radiation update.

The geometry is reconstructed from the saved physiological/articulatory
trace at 1.500 and 1.501 s. The sound is not taken from the old speech
recording. The run first holds the initial mouth for 1 ms, continues through
the actual 1 ms body/lip movement, then holds the final mouth while preserving
all source and airway memories. The first held interval and motion interval
complete, with work residuals below 4.07e-16 J and 1.24e-16 J, and air residuals
below 1.55e-20 and 2.87e-20 m³ respectively.

The longer run uses a frozen executable and an exact 120-file source snapshot
under `out/actual-coupled-phonation-v2`. It precedes one additional finite-value
guard for unrepresentably small positive gaps; the current build's short
open-mouth control separately checks the valid geometry. The original failed
attempt remains under `out/actual-coupled-phonation`.

This is progress toward voiced articulation, not a new NATO pronunciation
result. The mouth still uses two oral pressure regions, with four lip and
eight throat sections. General changing topology, further spatial resolution,
physiological control, speech realism and production efficiency remain open.

The separate fixed-tube reference has now completed 20 ms with five closings
and reopenings. Its work residual stays below 2.25e-15 J and its air residual
below 4.46e-21 m³ in the completed 15 ms continuation. Its 250 ms continuation
is running from the preserved endpoint. That reference and the actual-mouth
run have separate executables and records; neither is presented as an
intelligible speech recording.
