#!/usr/bin/env python3
"""Embed every spatial cell below 8 kHz in a crop/pan viewer; never resize fields."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import numpy as np


def encode(path, label):
    with np.load(path) as data:
        seconds = data['seconds']
        if len(seconds) < 2 or not np.allclose(np.diff(seconds), np.diff(seconds)[0]):
            raise ValueError('Viewer requires at least two uniformly spaced time columns')
        result = {'label': label, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                  'cols': len(data['seconds']), 't0': float(data['seconds'][0]),
                  'dt': float(np.diff(data['seconds'])[0])}
        for name, key, axis in [('ordinary', 'ordinary_stft_magnitude', 'ordinary_hz'),
                                ('enhanced', 'enhanced', 'hz')]:
            raw = data[key]
            field = raw[data[axis] <= 8010]
            hz = data[axis]
            if (raw.shape != (len(hz), len(seconds)) or len(field) < 256
                    or not np.allclose(np.diff(hz), np.diff(hz)[0])
                    or hz[0] != 0 or not np.all(np.isfinite(raw))
                    or np.any(raw < -8*np.finfo(float).eps*max(float(raw.max()), 1e-30))):
                raise ValueError(f'{name}: expected a finite nonnegative uniform field with >=256 rows')
            # Display quantization only: <0.00046 dB error above the -60 dB floor.
            # Full precision linear fields remain in the original NPZ; no spatial
            # averaging, interpolation, decimation, or row/column omission here.
            db = 20 * np.log10(np.maximum(field / max(float(raw.max()), 1e-30), 1e-3))
            packed = np.round(np.clip((db + 60) / 60, 0, 1) * 65535).astype('<u2')
            result[name] = {'rows': field.shape[0], 'df': float(np.diff(data[axis])[0]),
                            'data': base64.b64encode(packed.tobytes()).decode('ascii')}
        return result


def native_crops(path, prefix):
    """Small PNGs with explicit pixel-to-cell indexing, no plotting resampler."""
    from PIL import Image
    receipt = {'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'width': 384, 'height': 256, 'column_pixels': 8, 'enhanced_row_pixels': 2,
               'display': 'grayscale, white=-60 dB, black=0 dB, independent field peaks'}
    with np.load(path) as data:
        df = float(np.diff(data['hz'])[0])
        first_row, first_column = 90, 66
        high = (first_row + 128 - .5) * df
        low = (first_row - .5) * df
        seconds = data['seconds']
        if len(seconds) < first_column + 48 or len(data['hz']) < first_row + 128:
            raise ValueError('Static Alfa crop is outside this input; omit --crop-prefix')
        dt = float(np.diff(seconds)[0])
        receipt.update(frequency_edges_hz=[low, high],
                       time_edges_s=[float(seconds[first_column] - dt/2),
                                     float(seconds[first_column+47] + dt/2)])
        for name, key, axis in [('ordinary', 'ordinary_stft_magnitude', 'ordinary_hz'),
                                ('enhanced', 'enhanced', 'hz')]:
            raw, hz = data[key], data[axis]
            # Every output pixel samples exactly one cell. Enhanced rows repeat
            # twice; ordinary rows occupy approximately four pixels at the SAME
            # physical frequency scale (the two grids are not precisely 2:1).
            frequencies = high - (np.arange(256) + .5) * df / 2
            rows = np.floor(frequencies / (hz[1] - hz[0]) + .5).astype(int)
            cols = first_column + np.arange(384) // 8
            field = raw[rows[:, None], cols[None, :]]
            db = 20 * np.log10(np.maximum(field / max(float(raw.max()), 1e-30), 1e-3))
            gray = np.round(np.clip(-db / 60, 0, 1) * 255).astype('uint8')
            target = prefix.with_name(prefix.name + '-' + name + '.png')
            target.parent.mkdir(parents=True, exist_ok=True)
            Image.fromarray(gray).save(target)
            receipt[name] = {'path': str(target), 'distinct_rows': len(np.unique(rows)),
                             'hz_per_row': float(hz[1] - hz[0])}
    prefix.with_suffix('.json').write_text(json.dumps(receipt, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('human', type=Path)
    parser.add_argument('generated', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--crop-prefix', type=Path,
                        help='Also export fixed native 384x256 human Alfa crops and an axis receipt')
    args = parser.parse_args()
    records = [encode(args.human, 'Human: valeatory'), encode(args.generated, 'DuckLips: reference body')]
    template = Path(__file__).with_name('field_viewer.html').read_text()
    result = template.replace('__FIELD_DATA__', json.dumps(records, separators=(',', ':')))
    if len(result.encode()) >= 1_000_000:
        raise ValueError('Viewer exceeds 1 MB; select shorter input excerpts, never downsample fields')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(result)
    if args.crop_prefix:
        native_crops(args.human, args.crop_prefix)
    print(json.dumps({'path': str(args.out), 'bytes': len(result.encode()),
                      'fields': [{k: v for k, v in r.items() if k not in ('ordinary', 'enhanced')} for r in records]}, indent=2))

if __name__ == '__main__':
    main()
