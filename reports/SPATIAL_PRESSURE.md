# Spatial airflow exposes a geometric sampling discontinuity

The canonical moving-network engine now has an opt-in `--spatial-pressure`
path, including the conservative material tongue. It connects the separately
tested pressure network to actual streamed speech. This is a physical
integration experiment, not an accepted pronunciation improvement. It reduces
the remaining Foxtrot /k/→/s/ waveform transient but introduces a larger one
later in the same word. The user's [14 legible-but-mechanical words](LISTENING_2026_10_09.md)
remain the only human score, attached to the earlier specified recording.

## Mean airflow on the actual airway tree

Each oral, nasal and side-cavity node stores its own `C p + V`. Closed
passages carry exactly zero mass. The glottal port solves against local rear
pressure and the canonical acoustic impedance; the vocal folds receive that
local load. Frication uses local constriction flow divided by its local area.
Mouth and nasal mean flows go to their respective existing radiators.
Reverse constriction flow does not receive the forward teeth-obstacle gain.
This does not introduce a recorded source, word-specific waveform, source
gain correction, limiter, or click-repair filter.

The existing control-block geometry endpoints define linearly interpolated
node volumes. Air compliance is `V/(rho c²)`. The body's existing static oral
wall compliance is distributed by reference section length. That distribution
is an explicit assumption; it is not measured anatomy. It also leaves wall
compliance at a zero-area oral node, which remains a contact-model limitation.

Passages obey `dp = R Q + lambda rho Q |Q|/(2 Aref²)`. Circular Poiseuille
resistance is distributed by length. A contraction from area `a` to `b < a`
contributes `lambda = (Aref/b)² - (Aref/a)²`; an expansion assumes zero
recovery. Forward and reverse coefficients are independent. Subdividing a
uniform section therefore adds no spurious orifice loss, and consecutive
monotone contractions telescope. Outlet static pressure is ambient; the
lip/nostril contraction is not counted again at the outlet.

These choices are deliberately uncalibrated. Measurements and a fitted
contraction/expansion model are reported by
[Birkholz and Häsner, Interspeech 2024](https://www.isca-archive.org/interspeech_2024/birkholz24_interspeech.html).
Their abstract reports partial pressure recovery and a fitted viscous area
exponent. This implementation does not adopt those fitted parameters or claim
agreement with their measurements. It first tests topology and conservation.

## Nonlinear solve and verification

Newton uses a positive tree Jacobian. Globalization follows the convex
storage-plus-passage potential, whose pressure gradient is the mass residual.
Direct subtraction of that potential lost precision near convergence; the
implementation evaluates gradient work plus a nonnegative Bregman remainder.
This preserves the same physical equations. Failure preserves the previous
accepted pressure/flow state and explicitly stops speech, with time,
residual and iteration count reported; it never substitutes the old cavity.

Each node must meet its own residual tolerance plus a bound from pressure
representation and flow sensitivity. Total mass has a separate residual
computed from storage and external ports, because internal flow uncertainty
cancels globally. A global maximum uncertainty alone was too permissive:
the repeated closure/reversal regression's mass error was 1.51e−13 m³.
The corrected local-plus-global criterion reduces it to 1.45e−17 m³.

Native checks cover independently summed steady viscous and nonlinear losses,
4/32-section refinement, reverse flow, sealed front/rear isolation, nasal
routing with a sealed mouth, passive pressure bounds, changing compliance,
rollback and allocation-free stepping. Streamed body/tip gestures produce
identical samples with whole buffers or 17-frame buffers and with diagnostics
enabled or disabled. The complete Release suite passed 13/13.
The affected streaming, cavity-network and mean-airway tests also passed
ASan/UBSan, 3/3, using `-O1 -g -fsanitize=address,undefined
-fno-omit-frame-pointer`. The initial unoptimized sanitizer run was stopped
before completion and is not counted as a pass. CLI checks additionally
verified a 481-frame prefix ending inside a host chunk and rejected invalid
prefix values or spatial mode without canonical mode.

The prior single-cavity canonical route was separately rerendered through
5.65 s and matches all 271,200 preceding stereo frames exactly. The new
`--until` option preserves the complete original alphabet plan rather than
replanning an isolated word. No performance optimization was attempted.

## What the speech experiment revealed

The complete new alphabet is `out/nato-spatial-pressure.wav`, rendered with
all three flags: `--canonical-motion --conservative-tongue --spatial-pressure`.
The 23.239-second run has zero geometry failures, airflow solve failures,
nonfinite samples and clipping. Maximum nodal mass residual is
7.797645398e−13 m³/s. Those results validate execution, not realistic speech.
The exact commands and hashes are in
[spatial_pressure_reference.json](spatial_pressure_reference.json).

| Same full-sequence interval | Prior closed-contact render | Spatial-pressure experiment |
|---|---:|---:|
| Foxtrot /k/→/s/, largest left sample jump | 0.237927 | 0.054910 |
| Foxtrot 5.023–5.029 s, largest left sample jump | 0.103355 | 0.805382 |
| Whole Foxtrot peak | 0.140136 | 0.624117 |
| Whole November, largest left sample jump | 0.021027 | 0.020938 |

These are unnormalized waveform measurements, not perceived click ratings.
The larger later Foxtrot transient prevents treating this as a successful
click fix. No pronunciation improvement is asserted for any codeword.

Complete enhanced fields for Foxtrot (4.55–5.65 s) and November
(11.70–12.65 s) are saved in
`out/{foxtrot,november}-spatial-pressure-enhanced.{npz,png,json}`.
They retain all 684 frequency rows and all 104/90 time columns, magma,
the saved reference 0.015974152821882408 and 70 dB. PNGs exactly match the
saved matrices' color mapping at one pixel per row and four per column.
No ordinary STFT companion is displayed. Exact full-context audio excerpts
and visual metadata are recorded in [spatial_pressure_media.json](spatial_pressure_media.json).

The geometry replay agrees with all 23,239 speech-probe attachment positions,
attachment motions and target errors to every printed digit. Its full airway
CSV identifies a particular discontinuity at 5.025 s:

- Section 58's measurement line is at x = 7.854501096 cm.
- The tip centre moves from x = 7.522187706 to 7.502406814 cm.
- The rounded pad radius is 0.35 cm: the line changes from hitting the pad
  to missing it.
- The measured section area jumps from 0.330394 to 8.979708 cm².
- That one section gains 2.365047 mL; the complete oral volume gains
  2.291309 mL in the 1 ms block.

This is the midpoint ray switching from the tip surface to a deeper boundary,
not the material body suddenly gaining volume. Similar switches occur at
5.071 s in section 52 and 5.073 s in section 51. The spatial airflow model
reacts strongly to the prescribed false volume changes.

An independent controlled rounded-pad/floor calculation isolates the sampling
problem. It uses native reference dimensions but omits the body, blade,
anatomy residual and smooth minimum; it is not a full speech geometry fit.

| Pad movement across a section's measurement line | Midpoint volume jump | Integrated section volume change |
|---|---:|---:|
| 0.02 cm | 2.362422 mL | 0.178188 mL |
| 0.002 cm | 2.336010 mL | 0.017819 mL |
| 0.0002 cm | 2.327180 mL | 0.001782 mL |
| 0.00002 cm | 2.324346 mL | 0.000178 mL |

The midpoint jump remains finite as motion tends to zero. The integrated
volume change tends to zero. Merely increasing temporal resolution therefore
cannot remove this spatial sampling defect. The integrated quantity is still
only the upper-gap proxy, not total physical air volume. A real pad has air
below it as well as above it; when its silhouette passes a ray, these regions
join. Total air space must not appear or disappear at that join. The repository's
older `floor_space` diagnostic already recognizes partial tip coverage, but
neither current acoustic assembly nor this pressure integration consumes its
sublingual fields. This experiment localizes the consequence in the current
speech engine; it does not establish that a cell average alone restores anatomy.

Earlier source-isolation replays
also showed that turning off frication suppresses the /k/→/s/ transient, while
doubling acoustic substeps does not. Muting frication is not a speech fix;
the doubled-substep render also differed overall and does not establish
numerical convergence.

## Remaining work

The next geometric correction must distinguish cell air volume from the
apertures that connect cells, integrate partial occupation as the tongue
crosses a cell, account for connected space above and below the tongue, and
preserve exact closure. Averaging an aperture indiscriminately
could reopen a real stop; smoothing the output would conceal the defect.
The current opt-in experiment deliberately retains the observed geometry so
that this failure and its consequences remain reproducible.

Further limitations include the empirical gap-to-area law, omitted lateral
and sublingual connectivity, kinematic rather than muscle-driven tongue
motion, one selected turbulent jet, block-rate source selection, unvalidated
mean/acoustic separation, and near-contact acoustic loss-memory transport.
Neither finite sample output nor a cleaner spectrogram establishes NATO
pronunciation or human realism.

Reproduce the measurement after running the saved alphabet command:

```sh
./build/ducklips_tongue_reach_audit bodies/reference.body samples/nato.dict \
  'Alfa. Bravo. Charlie. Delta. Echo. Foxtrot. Golf. Hotel. India. Juliett. Kilo. Lima. Mike. November. Oscar. Papa. Quebec. Romeo. Sierra. Tango. Uniform. Victor. Whiskey. X-ray. Yankee. Zulu.' \
  --airway out/nato-spatial-airway.csv > out/nato-spatial-geometry.csv
python tools/spatial_pressure_reference.py
```

The measurement script verifies source/binary launch hashes, replay identity,
exact waveform windows and the independent controlled integration. It does
not render or modify audio, and it does not assign intelligibility scores.
