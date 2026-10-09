# Tongue displacement and mean-cavity mass balance

The canonical speech renderer previously interpolated acoustic compliance and
branch inertance through each control interval, but its separate aerodynamic
cavity retained `C dP/dt = Q_in - Q_out` while geometry changed. It omitted the
volume swept by the moving tongue/lips and the change in stored air associated
with compliance. A sealed cavity could therefore shrink without developing
pressure. This affects the pressure drop driving folds, stops and frication.

The moving update conserves the existing model's linearized air quantity:

```
V_new + C_new P_new - (V_old + C_old P_old) = dt (Q_in - Q_out)
C = V / (rho c^2) + C_wall
Q_out = A_out sign(P_new) sqrt(2 abs(P_new) / rho)
```

`V` is geometric oral air volume and `C_wall` is the configured wall
compliance. Substitution leaves the same implicit quadratic in
`sqrt(abs(P_new))`, so expanding cavities and reverse flow are handled without
clamping pressure to zero. The new method accepts both endpoint compliances
and the geometric volume change. The canonical renderer interpolates these
at each audio sample, aligned with its acoustic geometry interval. On first
initialization both endpoints use the assembled initial body, avoiding an
artificial empty-to-full cavity transition.

This is conservation within a constant-bulk-modulus, small-pressure,
single-cavity approximation. It is not a distributed compressible flow solver
or a complete thermodynamic/moving-wall model. The acoustic perturbation and
mean pressure remain separate; this change does not claim a complete total
energy identity for their coupling.

## Evidence from speech and independent limits

- A sealed 50 mL cavity is compressed by 10 mL and expanded back. Its stored
  air is conserved within 1e-16 m^3; it develops over 900 Pa at maximum
  compression and returns to ambient within 1e-9 Pa.
- Vented compression/expansion verifies signed mass flow, including inward
  flow during expansion, within 1e-12 m^3/s.
- Streamed diphthongs/frication retain the mass-balance gate, no heap
  allocation, host-buffer identity and finite unclipped output.
- The updated canonical `Sierra` and `Papa. Quebec.` renders are finite and
  unclipped. The audited two-word render's maximum mass residual is
  2.4234e-16 m^3/s. The extra diagnostics do not change its PCM bytes.

The 2.21-second Papa/Quebec field retains all 684 frequency rows and 208 time
columns, original magma colors, 70 dB range and the existing synthetic display
reference. No ordinary-STFT companion or cropped field is delivered.

## What the audit exposes next

The Papa/Quebec gesture near 1.338 seconds reduces modeled oral volume at a
peak rate of about 3.167 L/s. Its cavity pressure remains bounded, but the
large modeled displacement needs an anatomical explanation. It must not be
silently capped or hidden by gain changes. The original coarse tract lacks
the detailed surface model's under-tongue airspace; that is one hypothesis to
check, not a proven explanation. Tongue deformation, passage geometry and
timing remain open parts of the requested physical body.

The frication trace also shows 19 source-branch changes in the baseline
Sierra word. Most occur at weak excitation. Source snapping and the existing
0.3 obstacle-activation switch remain unchanged pending a physical source
placement treatment; this work does not claim to fix them. The mass-balance
error was corrected first because it affects the actual source pressure.

No word-recognition or human-realism score is assigned. Performance work
remains deferred by the user; the integrated renderer is still an offline
reference and has not met the final 10x requirement.

Reproduce the speech-level audit:

```sh
./build/ducklips_say --canonical-motion --body bodies/reference.body \
  --dictionary samples/nato.dict --text 'Papa. Quebec.' \
  --out out/canonical-stops-volume-audit.wav \
  --trace 0 out/canonical-stops-volume-audit.csv
```

The source hashes, render logs, trace measurements and artifact hashes are
recorded in [moving_cavity_reference.json](moving_cavity_reference.json).

Follow-up [geometry attribution](TONGUE_GEOMETRY.md) identifies the rigid
tongue-body gesture as the source of the large Quebec displacement. It also
documents a separate corrected discontinuity at the onset of tip activation.
