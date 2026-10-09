#!/usr/bin/env python3
"""Fetch the small human-speech reference set, with strict byte/hash limits."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare', action='store_true', help='Also make 48 kHz mono analysis WAVs (requires numpy, scipy, soundfile)')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'samples/manifest.json').read_text())
    for item in manifest['recordings']:
        target = ROOT / item['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            if shutil.disk_usage(ROOT).free < 1024**3:
                raise RuntimeError('Keep at least 1 GiB free; no sample was downloaded')
            request = urllib.request.Request(item['url'], headers={'User-Agent': 'DuckLips/0.1 speech research'})
            with urllib.request.urlopen(request, timeout=45) as response:
                data = response.read(item['max_bytes'] + 1)
            if len(data) > item['max_bytes']:
                raise RuntimeError(f"Download exceeded limit: {item['id']}")
            if hashlib.sha256(data).hexdigest() != item['sha256']:
                raise RuntimeError(f"Source changed: {item['id']}; inspect before updating the manifest")
            temporary = target.with_suffix(target.suffix + '.part')
            temporary.write_bytes(data)
            temporary.replace(target)
        if hashlib.sha256(target.read_bytes()).hexdigest() != item['sha256']:
            raise RuntimeError(f'Hash mismatch: {target}')
        print(f"verified {item['id']}: {target.stat().st_size} bytes")
        if args.prepare:
            import numpy as np
            import soundfile as sf
            from scipy.signal import resample_poly
            from math import gcd
            signal, rate = sf.read(target, always_2d=True)
            # Select the first microphone; no stereo phase cancellation or loudness normalization.
            signal = signal[:, 0]
            divisor = gcd(rate, 48000)
            converted = resample_poly(signal, 48000 // divisor, rate // divisor)
            output = ROOT / 'samples/derived' / (item['id'] + '_48k.wav')
            output.parent.mkdir(parents=True, exist_ok=True)
            sf.write(output, converted, 48000, subtype='FLOAT')
            receipt = dict(source_sha256=item['sha256'], output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                           source_rate=rate, output_rate=48000, selected_channel=0, normalized=False,
                           seconds=len(converted)/48000, peak=float(np.max(np.abs(converted))),
                           modifications='First channel, scipy.signal.resample_poly to 48 kHz, float32 WAV',
                           license=item['license'], author=item['author'], source_page=item['source_page'])
            output.with_suffix('.json').write_text(json.dumps(receipt, indent=2) + '\n')
            print(f'prepared {output}: {receipt["seconds"]:.3f} s')
            if item['id'] == 'nato_valeatory':
                segments = json.loads((ROOT / 'samples/nato_segments.json').read_text())
                if segments['sourceSHA256'] != item['sha256']:
                    raise RuntimeError('Word boundaries belong to a different recording')
                clips = []
                for clip in segments['clips']:
                    start, end = round(clip['start'] * 48000), round(clip['end'] * 48000)
                    if not 0 <= start < end <= len(converted):
                        raise RuntimeError('Invalid word boundary')
                    path = output.parent / ('nato_' + clip['letter'].lower() + '.wav')
                    sf.write(path, converted[start:end], 48000, subtype='FLOAT')
                    clips.append(dict(word=clip['word'], start_sample=start, end_sample=end,
                                      path=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
                (output.parent / 'nato_clips.json').write_text(json.dumps(dict(
                    **receipt, segmentation_author='emonster / SkyStamp',
                    segmentation_source='https://skystamp.app/license/', clips=clips), indent=2) + '\n')
                print(f'prepared {len(clips)} labelled human word samples')


if __name__ == '__main__':
    main()
