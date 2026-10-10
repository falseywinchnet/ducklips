#!/usr/bin/env python3
"""Freeze the local engine and render fresh, auditable listening candidates."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
from scipy.io import wavfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--alphabet', action='store_true')
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for source, name in [('build/ducklips_say', 'ducklips_say'), ('bodies/reference.body', 'body.body'),
                         ('samples/nato.dict', 'nato.dict')]:
        shutil.copy2(ROOT/source, out/name)
    paths = [ROOT/'CMakeLists.txt', *sorted((ROOT/'src').glob('*.cpp')),
             *sorted((ROOT/'include/ducklips').glob('*.hpp')), Path(__file__),
             ROOT/'tools/ducklips_say.cpp', ROOT/'tools/ducklips_stream.cpp', ROOT/'tests/local_speech_tests.cpp']
    for path in paths:
        destination = out/'source-snapshot'/path.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    report = dict(scope='New local-engine speech candidates. Intelligibility unassessed; '
                        'excessive consonant transients remain unresolved.',
                  binary_sha256=sha(out/'ducklips_say'),
                  source_hashes={str(p.relative_to(ROOT)):sha(p) for p in paths},
                  body_sha256=sha(out/'body.body'), dictionary_sha256=sha(out/'nato.dict'), runs=[])
    words = dict(papa='Papa.', november='November.', foxtrot='Foxtrot.')
    if args.alphabet:
        words['alphabet'] = ('Alfa. Bravo. Charlie. Delta. Echo. Foxtrot. Golf. Hotel. India. Juliett. '
                             'Kilo. Lima. Mike. November. Oscar. Papa. Quebec. Romeo. Sierra. Tango. '
                             'Uniform. Victor. Whiskey. X-ray. Yankee. Zulu.')
    for name, text in words.items():
        command = [str(out/'ducklips_say'), '--local-flow', '--conservative-tongue',
                   '--body', str(out/'body.body'), '--dictionary', str(out/'nato.dict'),
                   '--text', text, '--full-scale-pa', '128', '--out', str(out/f'{name}-raw.wav')]
        started = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True)
        elapsed = time.perf_counter()-started
        log = result.stdout+result.stderr
        (out/f'{name}.log').write_text(log)
        entry = dict(name=name, text=text, command=command, returncode=result.returncode,
                     process_wall_seconds=elapsed, log_sha256=sha(out/f'{name}.log'))
        match = re.search(r'rendered ([0-9.]+) s in ([0-9.]+) s', log)
        if match:
            duration, rendering = map(float, match.groups())
            entry.update(duration_s=duration, render_seconds=rendering, realtime_multiple=duration/rendering)
        match = re.search(r'body construction ([0-9.]+) s', log)
        if match:
            entry['body_construction_seconds'] = float(match[1])
        if result.returncode == 0:
            rate, raw = wavfile.read(out/f'{name}-raw.wav')
            peak = float(np.max(abs(raw)))
            if not np.isfinite(raw).all() or not 0 < peak < 1:
                raise ValueError('invalid or clipped native waveform')
            gain = .8/peak
            pcm = np.rint(32767*gain*raw).astype(np.int16)
            wavfile.write(out/f'{name}.wav', rate, pcm)
            entry.update(sample_rate=rate, frames=len(raw), raw_peak=peak, full_scale_pressure_Pa=128,
                         fixed_playback_gain=gain, raw_sha256=sha(out/f'{name}-raw.wav'),
                         playback_sha256=sha(out/f'{name}.wav'))
        report['runs'].append(entry)
        (out/'run.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(entry), flush=True)
    if any(r['returncode'] for r in report['runs']):
        raise SystemExit('One or more candidates failed: no listening copy made for failed runs.')


if __name__ == '__main__':
    main()
