#!/usr/bin/env python3
"""Show the entire enhanced field in its original colors, without spatial reduction."""
import argparse
import base64
import hashlib
import io
import json
from pathlib import Path
import numpy as np
from PIL import Image
from enhanced_specgram.display import enhanced_rgba


def png_bytes(rgba):
    buffer = io.BytesIO()
    Image.fromarray(rgba).save(buffer, format='PNG')
    return buffer.getvalue()


def encode(path, label, png_prefix=None):
    with np.load(path) as data:
        field, seconds, hz = data['enhanced'], data['seconds'], data['hz']
        if (field.shape != (len(hz), len(seconds)) or len(seconds) < 2 or len(hz) < 2
                or not np.allclose(np.diff(seconds), np.diff(seconds)[0])
                or not np.allclose(np.diff(hz), np.diff(hz)[0])):
            raise ValueError('Expected complete uniformly spaced frequency/time axes')
        # Use the exact existing display settings, including explicit references.
        receipt = json.loads(path.with_suffix('.json').read_text())
        reference = float(data['display_reference'])
        db_range = float(receipt['db_range'])
        rgba = enhanced_rgba(field, reference, db_range)
        dt, df = float(seconds[1]-seconds[0]), float(hz[1]-hz[0])
        result = {'label': label, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                  'cols': len(seconds), 'rows': len(hz), 'db_range': db_range,
                  'display_reference': reference,
                  'time_edges': [float(seconds[0]-dt/2), float(seconds[-1]+dt/2)],
                  'frequency_edges': [float(hz[0]-df/2), float(hz[-1]+df/2)],
                  'png': 'data:image/png;base64,' + base64.b64encode(png_bytes(rgba)).decode('ascii')}
        if png_prefix:
            # Integer replication only. The full matrix occupies the PNG.
            target = png_prefix.with_name(png_prefix.name + '-' + label.split(':')[0].lower() + '.png')
            target.parent.mkdir(parents=True, exist_ok=True)
            Image.fromarray(np.repeat(rgba, 4, axis=1)).save(target)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('human', type=Path)
    parser.add_argument('generated', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--png-prefix', type=Path, help='Also export both full-field color PNGs')
    args = parser.parse_args()
    records = [encode(args.human, 'Human: valeatory', args.png_prefix),
               encode(args.generated, 'DuckLips: reference body', args.png_prefix)]
    # Exact same magma color map as the field, left to right from floor to peak.
    import matplotlib
    ramp = matplotlib.colormaps['magma'](np.linspace(0, 1, 160), bytes=True)[None, :, :]
    legend = 'data:image/png;base64,' + base64.b64encode(png_bytes(np.repeat(ramp, 12, axis=0))).decode('ascii')
    template = Path(__file__).with_name('field_viewer.html').read_text()
    result = template.replace('__FIELD_DATA__', json.dumps({'records': records, 'legend': legend}, separators=(',', ':')))
    if len(result.encode()) >= 1_000_000:
        raise ValueError('Viewer exceeds 1 MB; use a standalone full-field image instead of reducing the data')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(result)
    print(json.dumps({'path': str(args.out), 'bytes': len(result.encode()),
                      'fields': [{k: v for k, v in r.items() if k != 'png'} for r in records]}, indent=2))

if __name__ == '__main__':
    main()
