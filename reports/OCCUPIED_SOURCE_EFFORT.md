# Measure source error on the occupied gas

The source stores two coefficients of an affine specific-enthalpy effort
field. Their units are Pa, but they are not necessarily the pressure at an
occupied location. When one end closes and the remaining air occupies a short
wedge, the coefficient at the closed end extrapolates the field into tissue.

The old coupled error controller compared both coefficients directly. At the
saved 28.954 ms fixed-tube source state, a 2 ns full-step/two-half-step comparison gives
0.132547 Pa coefficient disagreement but only 0.00102687 Pa disagreement on
the occupied interval. Across five successful trial sizes, the ratio is
128–130. Larger trials also fail their nonlinear solve, so the error norm is
not the only remaining obstacle to closure.

## Changed measure

`ContourPressureStorage::effort_difference` returns the exact supremum of
the affine effort difference over the union of the two occupied supports,
including their contact limits. The union matters when the candidate tissue
positions differ: neither candidate may hide the other's occupied air.
An affine difference reaches its largest absolute value at an interval
endpoint, so at most four endpoint evaluations suffice. Contact weights are
formed directly from the two areas, avoiding subtraction of a number close
to one in a very short pocket.

This measure applies to both supported axial profiles because their pressure
and area use the same hat-coordinate interval. When both ends are open, the
new measure is exactly the old maximum coefficient difference. Empty support
is rejected. The controller still uses the same 1000 Pa effort scale and all
its existing area, gas-moment, tissue, flow and passive-memory checks. It still
commits two half steps without extrapolation and retains the same work and
contact rules. No physical force or loss law changed.

## Evidence

- An independent 80-digit Decimal oracle checks 2,403 field pairs, including
  1,200 swapped pairs. Maximum relative discrepancy is 2.22e-16; swapped pairs
  give identical distances. Analytic tests cover narrowing at either end,
  different occupied intervals, opposite contact sides and empty support.
- A fully open 10 microsecond trajectory is exactly identical with the old
  and new checks, including its complete saved state and receipts.
- From the same saved near-contact state, the 0.1 microsecond trajectory at
  tolerance 1e-4 uses 32 accepted leaves instead of 128. Its occupied-field
  endpoint difference from the old fine trajectory is 0.279 Pa. Tightening
  the new tolerance to 2.5e-5 reduces that difference to 0.0557 Pa. At
  1.5625e-6 it uses 128 leaves and reproduces the old trajectory endpoint
  exactly. Tolerances control local disagreement, not a global error bound.
- A fresh continuation advances 2 microseconds from the saved state to a
  total 28.956 ms. It retains positive gas, zero flow through the sealed lower
  port, complete passive memory and checked work/air balance. The maximum
  work residual is 8.20e-16 J; mass residual is 5.85e-22 m3. The upper opening
  remains positive at 1.339e-8 m2 and is still closing. This is not completion
  of a closure/reopening cycle.
- The storage test plus five coupled/contact regression tests pass in Release
  and sanitizer builds. The standalone Decimal oracle uses the final native
  field-distance implementation.

`tools/joint_contact_probe.cpp --compare-step SECONDS` exposes the raw full
and half-step source states for this diagnosis. It performs an unadapted
comparison and reports failure if either candidate cannot be solved; it does
not loosen the controller. `tools/occupied_effort_reference.py` reproduces the
independent field-distance calculation.

Frozen sources, executables, input states and traces are under
`out/source-occupied-effort` and `out/source-effort-error`. The committed
`occupied_source_effort.json` binds the completed evidence. A separate
`last-pocket-run.json` tracks the still-further continuation; it is not part
of the completed two-microsecond result above. No new speech recording or
intelligibility result is claimed.
