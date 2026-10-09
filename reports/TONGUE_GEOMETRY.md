# Geometry attribution and continuous tip activation

## The large Quebec displacement

The new optional geometry diagnostics evaluate the actual control pose and
three diagnostic alternatives: include the existing floor-space estimate,
remove the body gesture, or remove the tip gesture. A separately prepared
TractModel owns these calculations. They do not enter the acoustic state.
`--diagnostics --probe FILE --chunk 48` records control poses at 1 kHz; without
`--diagnostics`, geometry columns are omitted rather than reported as zero.

In the Papa/Quebec render, the body gesture lifts the rigid tongue circle
from about 5.32 cm to 7.46 cm vertically. The oral air volume falls from about
81.19 mL to 32.27 mL through this gesture and the following vowel transition.
At the maximum displacement rate (control endpoint 1.339 s), volume changes
by -3.167 mL over one millisecond. Removing the body gesture leaves volume at
81.19 mL in this region; removing the tip gesture gives the original result.
Tip activation is zero. The existing floor-space estimate is also zero at
that peak and only grows to about 1.9 mL later.

Thus the existing under-tongue estimate does not explain away this change.
It is the programmed rigid-body motion, not a numerical cavity-balance
failure. This does not establish that the rigid-body model is anatomically
correct. Local dorsum deformation, contact constraints and tissue-volume
redistribution still need a physical treatment. No displacement cap, gain
change or empirically chosen speed reduction was applied.

The diagnostic render's common PCM prefix is byte-identical to the earlier
render. Their final lengths differ only because the CLI finishes on the
requested host chunk boundary (48 versus 480 frames). Streaming tests now
compare geometry diagnostics enabled against disabled while also changing
host buffer size, and assert no allocation during diagnostic rendering.

## A separate discontinuity in the tip

The old pose code blended the resting pad centre toward a target surface,
then subtracted the entire pad radius whenever activation was nonzero. The
zero-activation pose used the resting centre directly. Consequently even an
activation of 1e-9 moved the pad by 3.5 mm and dragged the body by 1.429 mm
in the measured neutral apical fixture. Smooth activation tracking could not
remove that geometric discontinuity.

Both interpolation endpoints now denote pad centres. The contact target is
`wall - normal * (gap + pad_radius)`, and the activation blend moves from the
resting centre to that target. Pressing above full activation retains its
existing definition. No interpolation rate or target gap was changed.

The native before/after probe reports:

| Activation | Old tip displacement (cm) | New tip displacement (cm) |
|---|---:|---:|
| 0 | 0 | 0 |
| 1e-9 | 0.35 | 0 |
| 1e-7 | 0.35 | 8.98e-14 |
| 0.001 | 0.349990 | 8.99e-6 |
| 0.01 | 0.349004 | 8.93e-4 |
| 0.1 | 0.256419 | 0.083932 |
| 1 | 2.997566 | 2.997566 |

Regression tests exercise apical, laminal and intermediate tip postures at
three places. They verify the limit of both tip and body positions at zero
activation and the specified surface gap at full activation. Frozen vowel
fitting has no active tip gesture, so this correction acts on gestured speech.
It does not solve the rigid-circle dorsum limitation identified above.

## Actual speech and reference view

Sierra is rendered through the canonical engine after the correction. It has
zero nonfinite and clipped samples, peak 0.044 and maximum moving-cavity mass
residual 2.04e-16 m^3/s. That is not a human pronunciation score.

The enhanced viewer contains complete Sierra fields: the attributed valeatory
human word (684 rows, 84 columns, 0.888 s), and generated speech (684 rows,
132 columns, 1.40 s). It retains the human's own saved amplitude reference,
the existing synthetic display reference, magma colors and 70 dB range.
Both full time axes are retained, with no waveform splice or time alignment.
There is no ordinary-STFT companion or grayscale/cropped view.

Reproduce:

```sh
./build/ducklips_say --canonical-motion --body bodies/reference.body \
  --dictionary samples/nato.dict --text 'Papa. Quebec.' --diagnostics \
  --chunk 48 --probe out/canonical-geometry-attribution.csv \
  --out out/canonical-geometry-attribution.wav
./build/ducklips_say --canonical-motion --body bodies/reference.body \
  --dictionary samples/nato.dict --text 'Sierra.' \
  --out out/canonical-sierra-continuous-tip.wav
```

Source hashes, the native limit comparison, geometry measurements and render
logs are recorded in [tongue_geometry_reference.json](tongue_geometry_reference.json).
Performance work remains deferred by user instruction. Anatomical fidelity,
perceptual pronunciation and the final throughput goal remain unfinished.
