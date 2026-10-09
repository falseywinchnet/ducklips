# Horizontal lines in the synthetic enhanced field

The full enhanced field for `out/nato_alfa.wav` shows conspicuous horizontal
components in regions where the human reference is mostly empty. They are
real components of the generated PCM. In the quiet 0.16–0.28 s interval their
peaks lie at almost exactly 1, 2, 3, 5, 6 and 7 kHz. The compact engine's default
control interval is 48 samples at 48 kHz: a 1 kHz clock.

Changing only that interval to 64 samples moves the comb to multiples of
750 Hz; 96 samples produces multiples of 500 Hz. These are clock-related
artifacts, not a set of human formants or a visualization-only effect. The
unmodified 48-sample reproduction is sample-identical to the PCM saved inside
`out/enhanced_specgram/nato_generated.npz`.

## Isolated mechanisms

**Physiological forces:** the body trackers had smooth targets, but the engine
advanced them once per control block and held their resulting forces constant
for every sample in that block. A smooth lung-pressure target became a staircase
at the glottal flow equation. A temporary intervention that only interpolated
lung pressure reduced the strongest quiet-section tones by 39–41 dB.

The retained change separates geometry advancement from physiological
advancement. Lung pressure, glottal resting aperture, intended pitch, effort,
and loudness now advance through their existing third-order tracker dynamics
at each audio sample. The geometry solver retains its current cadence. The
pitch feedback multiplier is recalculated when its controller changes, rather
than calculating an exponential per sample. No output EQ, notch filter or gate
was added. The temporary pressure ramp was removed.

Measured with a Hann-windowed sinusoid projection of channel 0, reporting
peak amplitude relative to full-scale amplitude 1:

| Interval | Frequency | Before | Continuous physiology | Reduction |
|---|---:|---:|---:|---:|
| 0.16–0.28 s | 1 kHz | -76.06 dBFS | -117.95 dBFS | 41.89 dB |
| 0.16–0.28 s | 2 kHz | -69.17 dBFS | -108.80 dBFS | 39.63 dB |
| 1.05–1.20 s | 1 kHz | -84.26 dBFS | -115.15 dBFS | about 31 dB |
| 1.05–1.20 s | 2 kHz | -77.87 dBFS | -107.95 dBFS | about 30 dB |

**Moving tract:** strong components remain during the preparatory tract motion
around 0.48–0.58 s. With continuous physiology, the 4 kHz component is about
-57.2 dBFS. A diagnostic that freezes the tract after the first geometry solve
reduces it to about -123.3 dBFS. That isolates the geometry-update path as
another contributor. It does not make freezing geometry an acceptable fix.
The diagnostic was removed, and articulation remains enabled.

Complete removal is **not claimed**. The retained changes address the
physiological forcing staircase. Continuous geometry/state transport remains
an open requirement, including the projection and changing source/radiation
weights. Slowly varying shimmer/jitter multipliers also still update at
control rate and should be included in the remaining modulation audit.

## Evidence and validation

- Raw measurements and WAV hashes: [control_clock_lines.json](control_clock_lines.json).
- Baseline code: commit `b755253d7fbab3cbdf154396628739449a89f36c`.
- Reusable PCM tone measurement: `tools/analyze_control_lines.py`.
- A rendered-Alfa regression checks the quiet 1/2 kHz clock components below
  1e-5 peak full-scale amplitude, while requiring audible voiced output. The
  old 2 kHz component is about 3.48e-4, so the stored baseline fails that bound.
- Existing streaming checks still establish allocation-free rendering and
  bit-identical output for whole-buffer versus 17-frame host chunks.
- Release and ASan/UBSan each pass the numerical and streaming suites.
- All seven compact-engine performance cases pass three repeated runs:
  minimum speeds span 12.48–19.50×. This covers the five presets, configured
  reference body and double speech rate. See
  [control_clock_performance.json](control_clock_performance.json).
- CLI tap names are now validated. An unrecognized name previously silently
  selected microphone output, which could invalidate a source-isolation test.

The new full enhanced field is
`out/control_lines/continuous-enhanced.png` (684 rows × 123 time columns,
columns displayed at 4 pixels each). It uses the original magma map, 70 dB
range and the **same saved amplitude reference** as the old synthetic field.
No rows or time columns were discarded. The original views and audio remain
available for comparison. The remaining transition-region lines are visible
in the corrected field and have not been hidden by changing its display scale.

The human comparison uses the samples stored with the original enhanced-field
NPZ. The current derived human WAV no longer matches that older metadata's
file hash; it was not silently substituted in this investigation.
