# Enhanced-field spectrogram for DuckLips visual analysis

This is the preserved Cleanup2 enhanced double-IRFFT speech field. Use it to
visually examine DuckLips output: harmonic continuity, formant trajectories,
onsets, bursts, frication, and transitions between articulatory gestures.
The requested deliverable is an offline visualization tool for a sibling agent
examining a generative speech system. It does not modify DuckLips audio.

## Run on the Mac mini

The existing Mini runtime has NumPy, SciPy, and Matplotlib; no installation is
needed. From `/Users/joshuahkuttenkuler/ducklips`:

```sh
/Users/joshuahkuttenkuler/miniforge3/bin/python tools/enhanced_specgram/render.py \
  /absolute/path/to/generated-speech.wav \
  --channel 0 --max-hz 8000 --out out/enhanced_specgram/speech
```

Outputs: `speech.png` (the full enhanced field alone, original magma colors), `speech.npz` (linear
fields and exact axes), `speech.json` (settings, timing, source hash, display
references). Existing filenames are overwritten if the same output prefix is
reused; choose a different prefix for each experiment. The default includes
the entire WAV. Use `--start` and `--duration` only for an explicitly requested
excerpt; the JSON records its exact sample bounds.

The PNG preserves every computed frequency row and time column. Each row is
one image pixel high and each column is four pixels wide; there is no subplot
layout or image reduction. Use the full-field viewer in `../native_fields.py`
for labelled axes at native cell size. Narrow displays show the complete field
vertically rather than cropping it. Keep the original magma map, saved display
reference, and dB range. Ordinary STFT and intermediate arrays remain available
in the NPZ for targeted numerical diagnostics; they are not companion panels
in the speech inspection view.

The default preserves the original reference geometry: apertures 1024/2048/4096,
hop 512 samples, target N=2048, and 256 rows. At 48 kHz it covers 0–2989.741 Hz.
Use `--max-hz 8000` for wider speech inspection at a suitable sample rate, e.g.
fricatives and higher formants. That changes the field crop and registration
context, so record it as an explicit configuration change. Larger crops and long
clips cost more, but the renderer does not silently shorten them.

Use `--channel 0` or `--channel 1` to inspect individual microphones. By default
stereo is averaged, which can cancel phase-opposed components. PCM is converted
to floating point without per-clip loudness normalization. The input is not
resampled: window durations and frequency spacing depend on its actual rate.
Use the same sample rate and parameters for controlled comparisons.

`--hop` explicitly changes the spacing of computed time centers, while keeping
the apertures and frequency rows unchanged. For short source-motion diagnostics,
`--hop 48` at 48 kHz gives 1 ms center spacing. These are new field observations,
not an interpolated larger image; the JSON and NPZ retain the chosen spacing.
Denser centers do not shorten the 21.3–85.3 ms windows or establish finer
physical resolving power. A 50 ms diagnostic therefore still has substantial
window-padding effects and is unsuitable for claims of steady harmonic quality.

## What to call directly

The self-contained algorithm files are `oracle/uncertainty_fusion.py`,
`oracle/fourier_unrolled_densified.py`, and `oracle/circles.py`. They require
NumPy and SciPy. Matplotlib is only needed by the renderer.

```python
# With tools/enhanced_specgram on sys.path:
from oracle.fourier_unrolled_densified import (
    fourier_unrolled_densified, reassigned_texture_baseline,
)

# samples: mono float waveform; preserve actual sample_rate separately.
fusion = fourier_unrolled_densified(samples)
result = reassigned_texture_baseline(fusion.fused, fusion.support_union)
enhanced = result.merged                   # row by time, nonnegative amplitude
times = fusion.centers / sample_rate
hz = np.arange(enhanced.shape[0]) * sample_rate / (2 * (2048 - 1))
```

Import NumPy as `np` for the axis expression. The renderer provides WAV loading,
channel selection, excerpting, saving, and plotting around these exact calls.

## The exact construction

1. Hann-window a centered frame and take its real FFT. For each complex STFT
   column, `recovered = irfft(spectrum, n=N)` recovers the real windowed frame.
2. Apply another `irfft(recovered)` with its **default output length**
   `L = 2*(N-1)`, then take the absolute value and low-row crop. For N=2048,
   L=4094, not 4096. See `double_irfft_low_rows` in `uncertainty_fusion.py`.
   This second operation treats the recovered samples as real spectral
   coefficients: it yields a cosine observation, not a second audio recovery.
3. Repeat at center offsets `(-379, -241, -64, 0, 83, 214, 427)` samples. A
   cosine null in one shifted view may coincide with a ridge in another.
   Normalize each view by its 99.5th percentile into the central view's scale;
   estimate local displacement with the Fourier-circle registration routines;
   warp the **linear amplitude** fields into reference coordinates and take
   their per-pixel maximum. Logs prepare registration views only.
4. Perform that construction separately for N=1024, 2048, and 4096. Map their
   rows using the exact N-1 frequency relationship; register and take a maximum
   across apertures. Three-FFT reassigned power supplies registration detail.
   `fourier_unrolled_densified(...).fused` is this multiscale field;
   `.support_union` is the corresponding maximum of transported power support.
5. The **final enhanced reference** is `reassigned_texture_baseline(...)`:

   ```text
   cartoon = 3-by-3 uniform average(multiscale field)
   texture = sqrt(reassigned support power)
   scale   = Q0.995(multiscale field) / max(Q0.995(texture), 1e-30)
   enhanced = cartoon + scale * texture
   ```

The final `.merged` field is the requested full reference. Merely calling
`fourier_unrolled_densified` omits the final average-plus-texture step. A
quadrature envelope `hypot(cosine, quadrature)` was an earlier distinct variant;
it is not a replacement for this construction.

## Axes and interpretation for the sibling

- Row r maps to `r * sample_rate / [2*(target_N-1)]`; time column t uses the
  actual frame center, `centers[t]/sample_rate`. The renderer adds the excerpt
  start time. Saved NPZ arrays `hz` and `seconds` are the authoritative axes.
- For 48 kHz and N=2048, rows are 11.72447 Hz apart and hop 512 is 10.66667 ms.
  These are sampling-grid spacings, not measured two-source resolving power.
  The windows span about 21.3, 42.7, and 85.3 ms.
- This is a deterministic transformed field. Its apparent fine structure is
  affected by cosine phase, registration, maximum fusion, reassignment, and
  the explicit 3-by-3 average. The saved intermediate arrays allow targeted
  numerical checks when needed; enhanced-field inspection uses the enhanced
  field itself without automatically adding conventional comparison panels.
- The final field has no complex synthesis phase, is not calibrated PSD or
  energy, and is not an invertible replacement for a waveform. Inspect the
  existing speech system through it; do not feed it into that system as audio.
- The renderer maps each enhanced cell directly to its original magma RGBA
  color and saves the entire matrix. It performs no spatial interpolation or
  reduction. The original historical five-row comparison used independent log
  normalization per representation; that layout is not the inspection default.
- For comparisons across renders, pass the same `--display-reference VALUE`
  from a saved report and the same `--db-range`. This fixes plotting contrast;
  internal quantile matching remains part of the original algorithm, so image
  intensity is not an absolute acoustic-level measurement.
- Excerpts are zero-padded at their own boundaries. At 48 kHz the largest
  half-window plus offset reaches about 52 ms beyond a center. Include extra
  context and inspect the interior; registration and quantile estimates also
  depend on the selected excerpt. Use identical excerpts for comparisons.
- NPZ keys separate `literal`, `registered_subhop`, `multiscale`,
  `support_power`, `cartoon`, `texture`, and `enhanced`. Use those arrays for
  additional zooms, without confusing a denser display with new measurements.

## Provenance and the live-filter distinction

Recovered 2026-10-09 from the Mini's preserved tree:
`/Users/joshuahkuttenkuler/Developer/CodexBuilds/cleanup2-a07b457b1144/research/true_superresolution/oracle`.
The original source was `bfft/experiments/ostensibly_frontend` and
`bfft/personal_deblurrer/circles.py`, associated with chat
`01a03b18-6686-79a3-87b1-ddee0b5a7b14` ("Reconstruct 3D scene from image").
Recorded source commit: `6e64386f2da658a90493733c08b01dea6425ecdd`.
That source had local changes, so the file hashes are the precise identity.

`source-provenance.json` preserves the upstream manifest. The three algorithm
files and LICENSE in this package match its snapshot SHA-256 values exactly.
The manifest also names VAD files, which are intentionally not needed for this
visualization package. `render.py` and this README are new integration aids.
The upstream LICENSE is preserved under `oracle/LICENSE`.

Cleanup2's `docs/true-superresolution.md` explicitly records the full enhanced
field as the selected comparison reference, with registered single-aperture
sub-hop magnitude as a cheaper live candidate. Its later `Filter.py` defaults
to `GUIDE_REGISTRATION=False`, retaining only the central double-inverse
observation. Thus the live default does **not** reproduce the full enhanced
visual reference. This package deliberately calls the full reference pipeline.

Historical radio-fixture output is retained in Cleanup2's
`build/true-superresolution/probe.npz` under `trace_field`. A reproducibility
check and example render are recorded in `verification.json` alongside this
guide. Those checks establish recovery/parity, not a speech-quality score.
