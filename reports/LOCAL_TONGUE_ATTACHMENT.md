# Local material motion at the blade attachment

The previous conservative-body integration left a 1.37 cm tip-target miss.
The cause is now localized: in the full NATO sequence, 47 control steps
(47 ms) around the /k/ → /s/ transition in Foxtrot required the blade to reach
forward while the dorsum was still raised. At the worst step, 4.954 s, body
activation was 0.9383 and tip activation 0.9064. Dorsal lift had pulled the
attachment backward, and the only response available to the fixed-length
blade was to drag the whole body through the closed dorsal contact.

This change adds a local, volume-preserving material motion around the
attachment. The blade retains its 2.2 cm length; nearby body tissue deforms
instead of translating the entire closed tongue. It is connected to the same
`--canonical-motion --conservative-tongue` speech path. No word-specific
waveform, gesture timing change, or output gain change supplies this result.

## Material field

Let the attachment start at c₀ and need displacement d. Its prescribed path
is c(t) = c₀ + t d, with t from zero to one. In coordinates z = x − c(t), use

```
w(z) = (1 − |z|²/R²)^4 inside radius R, zero outside
ψ(z) = (d_x z_y − d_y z_x) w(z)
v = (∂ψ/∂y, −∂ψ/∂x)
```

At the moving centre, v = d, so the attachment follows its requested path.
The field is divergence-free, and its value and first derivatives vanish at
the support edge. Material outside the swept support remains exactly fixed.
As with the dorsal map, extending this planar map with unchanged transverse
coordinate conserves three-dimensional volume locally, but does not supply
measured transverse anatomy.

The code uses implicit midpoint with 32 steps, propagates the deformation
gradient, and rejects a failed Newton solve. Each step has a trace-zero field
derivative; its Cayley derivative therefore has determinant one. Reverse
mapping traverses the same moving centres in reverse order. Airway ray
intersections and inside/outside classification use the composed dorsal and
attachment maps, including the composed inverse and gradient.

The support radius is derived from available geometry, not a fitted word
parameter. Bounds on distance to the oral wall, the lower anchor, and the
posterior edge of the oral region are evaluated along the attachment path.
Each bound is 1-Lipschitz; subtracting half the path sample spacing bounds the
whole path between samples. A small numerical margin keeps the support off
the fixed boundaries. A path with no admissible support fails explicitly.
The continuous field leaves those boundaries fixed. Discrete contact is
checked in the native fixtures and speech audits; this is not a global
collision certificate for arbitrary anatomy.

The conservative path now uses the same attachment mechanics at zero dorsal
activation. Previously that special case used the old rigid-body motion.
A regression checks continuity through zero with an active tip gesture.

## What is and is not established

The local field is a kinematic model. It does not minimize elastic energy,
compute muscle forces, or establish physiological timing. The stream still
sets targets; this field supplies a deformation that meets a formerly missed
geometric constraint while conserving material and retaining the lower anchor.
Large strains remain a concern. An exploratory worst-transition field with
2.2 cm support reached a maximum incremental principal stretch of 3.17; that
is not a justified human strain bound or a measurement of the final automatic
support rule.

The reason to permit deformation instead of treating the tongue as a rigid
link is consistent with continuum tongue modeling: intrinsic tissue can
elongate, shorten and bend while coordinating transverse deformation.
[Stavness, Lloyd and Fels (2012)](https://www.cs.usask.ca/faculty/stavness/papers/stavness2012-automatic-prediction-of-tongue-muscle-activations.pdf)
modeled these actions with a muscle-driven finite-element system and compared
selected predictions with experimental data. That paper supports the need
for deformable tissue; it does **not** validate this compact field or its
chosen attachment path.

The blade surface, pad contact, physiological source timing, and empirical
midsagittal-gap-to-airway-area law still need work. No NATO word is counted as
human-verified correct from these measurements.

## Native and speech evidence

- The local map's native grid test has maximum determinant error 2.89e-15,
  forward/inverse error 1.41e-15 cm, and independent finite-difference gradient
  discrepancy 6.45e-10. Material outside the support is unchanged exactly.
- Six combined surface cases (three dorsal activation levels, each with and
  without the tip gesture) agree with an independently intersected 8,192-point
  contour on hit classification; the largest gap difference is 9.35e-6 cm.
- The full alphabet geometry replay has maximum tip shortfall 1.83e-14 cm
  and maximum sampled body/wall overlap 2.51e-13 cm. The replay intentionally
  omits acoustic and physiological updates; it is geometry evidence.
- For a separate actual Foxtrot render, the replay's attachment coordinates,
  attachment-motion amount and tip miss match the speech probe at all 1,653
  control steps, to every printed digit. This connects the geometry audit to
  the geometry used in the actual acoustic solve.
- Actual Foxtrot has zero geometry failures, nonfinite samples or clipping,
  zero sampled body overlap, and maximum tip shortfall 1.83e-14 cm. The largest
  local attachment motion is 1.44082 cm. Its peak is 0.467 without normalization.
  The prior commit's matching word render has a 1.36951 cm tip miss and peak
  0.607. Those levels are output measurements, not intelligibility scores.

The previous-commit reference was rebuilt from `f39e834`, with the same body,
dictionary and canonical acoustic flags. It used 480-frame output chunks;
the new probe used 48-frame chunks, so the old file contains 7 ms of extra
terminal padding. Phone cues and timing are identical. The streaming tests
separately check host-chunk independence.

Final test scopes, alphabet audio results, hashes and source provenance are
in [local_tongue_attachment_reference.json](local_tongue_attachment_reference.json).

The actual full alphabet render completed 23,239 control blocks (23.239 s)
with zero geometry failures, nonfinite samples or clipping. Maximum tip
shortfall was 1.8208e-14 cm; maximum sampled body overlap was 2.5047e-13 cm.
Maximum local attachment motion was 1.82127 cm. Peak output was 0.530 without
normalization, and maximum moving-cavity mass residual was 5.5855e-16 m³/s.
The render's native source and binary hashes match the final implementation.
Its audio is `out/nato-local-material.wav`.

The final Release suite passed 10/10. ASan/UBSan passed the affected streaming,
material-flow and surface tests (3/3), followed by a rebuilt surface test
including the added active-tip continuity regression (1/1). No full sanitizer
suite beyond those affected tests is claimed for this change.

## Inspect and reproduce

```sh
cmake --build build
ctest --test-dir build --output-on-failure
./build/ducklips_tongue_reach_audit bodies/reference.body samples/nato.dict \
  'Foxtrot.' > out/foxtrot-geometry-local.csv
./build/ducklips_say --canonical-motion --conservative-tongue \
  --body bodies/reference.body --dictionary samples/nato.dict --text 'Foxtrot.' \
  --chunk 48 --out out/foxtrot-local-material.wav --probe out/foxtrot-local-material.csv
python tools/enhanced_specgram/render.py out/foxtrot-local-material.wav \
  --out out/foxtrot-local-material-enhanced --channel 0 --max-hz 8000 \
  --display-reference .015974152821882408
```

`out/foxtrot-local-material-enhanced.html` contains the full human and synthetic
fields. Human: 684×104 cells, samples 0–53,040. Synthetic: 684×155 cells,
samples 0–79,344. Each PNG preserves every row and column in magma at 70 dB;
its pixels were compared exactly with the saved linear field. The display
references are retained from the prior human/synthetic inspections and differ
between recordings, so brightness is not a calibrated relative level. There
is no conventional companion plot, grayscale conversion or shortened excerpt.
The human reference is valeatory's Foxtrot, credited in
[samples/README.md](../samples/README.md).

Visual inspection still shows a long synthetic pre-speech interval and
substantial differences in voicing and spectral distribution. Solving the
attachment reach does not resolve those remaining speech failures.
