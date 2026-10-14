#!/usr/bin/env python3
"""Render the preserved Cleanup2 enhanced speech field for visual inspection."""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np
from scipy.io import wavfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from display import enhanced_rgba

from oracle.fourier_unrolled_densified import (
    DEFAULT_APERTURES, DEFAULT_OFFSETS,
    fourier_unrolled_densified, reassigned_texture_baseline,
)
from oracle.uncertainty_fusion import centered_stft_lattice, double_irfft_low_rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("wav", type=Path)
    p.add_argument("--out", type=Path, required=True, help="Output filename prefix")
    p.add_argument("--start", type=float, default=0.0, help="Start in source seconds")
    p.add_argument("--duration", type=float, help="Explicit excerpt duration; default includes the rest of the WAV")
    p.add_argument("--max-hz", type=float, help="Highest row to include; default original 256 rows")
    p.add_argument("--channel", type=int, help="Zero-based channel; default arithmetic mono mix")
    p.add_argument("--hop", type=int, default=512, help="Time-center spacing in samples; default 512, no image resampling")
    p.add_argument("--db-range", type=float, default=70.0)
    p.add_argument("--display-reference", type=float, help="Fixed positive amplitude reference for enhanced-family panels")
    a = p.parse_args()
    if a.start < 0 or (a.duration is not None and a.duration <= 0) or a.db_range <= 0 or a.hop <= 0:
        p.error("start must be nonnegative; duration, db-range and hop must be positive")
    if a.display_reference is not None and a.display_reference <= 0:
        p.error("display-reference must be positive")
    fs, pcm = wavfile.read(a.wav)
    first = round(a.start * fs)
    last = len(pcm) if a.duration is None else min(len(pcm), first + round(a.duration * fs))
    if last <= first:
        p.error("excerpt has no samples")
    x = pcm[first:last]
    if x.dtype == np.uint8:
        x = (x.astype(np.float64) - 128.0) / 128.0
    elif np.issubdtype(x.dtype, np.integer):
        info = np.iinfo(x.dtype)
        x = x.astype(np.float64) / max(abs(info.min), info.max)
    else:
        x = x.astype(np.float64)
    if a.channel is not None:
        count = 1 if x.ndim == 1 else x.shape[1]
        if not 0 <= a.channel < count:
            p.error("channel index is outside this WAV")
        if x.ndim == 2:
            x = x[:, a.channel]
    elif x.ndim == 2:
        x = x.mean(axis=1)
    if not np.all(np.isfinite(x)):
        p.error("input contains nonfinite samples")
    n, hop = 2048, a.hop
    df = fs / (2 * (n - 1))
    rows = 256
    if a.max_hz is not None:
        if not 0 < a.max_hz <= fs / 2:
            p.error("max-hz must be positive and no greater than Nyquist")
        rows = min(n, int(np.ceil(a.max_hz / df)) + 1)
    if rows < 8 or len(x) < 4 * hop:
        p.error(f"use at least 8 frequency rows and an excerpt of at least {4 * hop} samples")
    started = time.perf_counter()
    fusion = fourier_unrolled_densified(x, target_crop_rows=rows, hop_length=hop)
    texture = reassigned_texture_baseline(fusion.fused, fusion.support_union)
    spectrum = centered_stft_lattice(x, n_fft=n, hop_length=hop)
    literal = double_irfft_low_rows(spectrum, n_fft=n, crop_rows=rows)
    registered = fusion.aperture_fields[DEFAULT_APERTURES.index(n)]
    ordinary = np.abs(spectrum)
    seconds = first / fs + fusion.centers / fs
    hz = np.arange(rows) * df
    reference = a.display_reference or max(float(np.max(texture.merged)), 1e-30)
    ordinary_reference = max(float(np.max(ordinary)), 1e-30)
    fields = dict(literal=literal, registered_subhop=registered,
                  multiscale=fusion.fused, support_power=fusion.support_union,
                  cartoon=texture.cartoon, texture=texture.texture,
                  enhanced=texture.merged, ordinary_stft_magnitude=ordinary)
    if any(not np.all(np.isfinite(v)) for v in fields.values()):
        raise RuntimeError("nonfinite output field")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(str(a.out) + ".npz", **fields, seconds=seconds, hz=hz,
                        sample_rate=fs, samples=x, hop=hop, target_n_fft=n,
                        display_reference=reference,
                        ordinary_hz=np.fft.rfftfreq(n, 1 / fs))
    # The delivered image is ONLY the full enhanced field. Figure layouts can
    # silently squeeze hundreds of rows into too few pixels even with nearest
    # interpolation. Map cells directly to RGBA, repeat columns by an integer,
    # and save without a plotting resampler. Exact axes remain in NPZ/JSON and
    # the full-field viewer adds axes without scaling away any image cells.
    rgba = enhanced_rgba(texture.merged, reference, a.db_range)
    plt.imsave(str(a.out) + "-cells.png", rgba)
    plt.imsave(str(a.out) + ".png", np.repeat(rgba, 4, axis=1))
    report = dict(source=str(a.wav.resolve()),
                  source_sha256=hashlib.sha256(a.wav.read_bytes()).hexdigest(),
                  sample_rate=int(fs), start_sample=first, end_sample=last,
                  channel=a.channel if a.channel is not None else "mean",
                  apertures=DEFAULT_APERTURES, offsets=DEFAULT_OFFSETS,
                  target_n_fft=n, hop=hop, rows=rows, row_spacing_hz=df,
                  highest_row_hz=float(hz[-1]), shape=list(texture.merged.shape),
                  display_reference=reference, ordinary_display_reference=ordinary_reference,
                  db_range=a.db_range, texture_amplitude_scale=texture.texture_amplitude_scale,
                  png=dict(field="enhanced", colormap="magma", row_pixels=1,
                           column_pixels=4, full_field=True),
                  elapsed_seconds=time.perf_counter()-started,
                  use="offline visual inspection of generated speech; enhanced is not a complex synthesis spectrum")
    Path(str(a.out) + ".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
