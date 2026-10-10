#!/usr/bin/env python3
"""Export a completed fresh word and its full enhanced field, optionally waiting."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from PIL import Image
from scipy.io import wavfile

from actual_coupled_audio import export
from coupled_phonation import ROOT, save, sha
from native_fields import enhanced_rgba


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--wait', action='store_true')
    args = parser.parse_args()
    run_path, out = args.run.resolve(), args.out.resolve()
    if out.exists():
        raise ValueError('use a new listening-output directory')
    while True:
        manifest = json.loads(run_path.read_text())
        if not manifest.get('word'):
            raise ValueError('this exporter requires an articulated word run')
        if manifest['status'] == 'complete':
            break
        if manifest['status'] != 'running' or not args.wait:
            raise ValueError(f'word render is {manifest["status"]}: {manifest.get("error", "not complete")}')
        time.sleep(10)
    # No audio or listening directory is created for a running/failed prefix.
    audio = export(run_path, out/'speech', peak=.8)
    command = [sys.executable, str(ROOT/'tools/enhanced_specgram/render.py'),
               str(out/'speech.wav'), '--channel', '0', '--max-hz', '8000',
               '--hop', '48', '--out', str(out/'enhanced')]
    subprocess.run(command, check=True)
    visual = json.loads((out/'enhanced.json').read_text())
    with np.load(out/'enhanced.npz') as fields:
        pixels = np.repeat(enhanced_rgba(fields['enhanced'], visual['display_reference'],
                                        visual['db_range']), 4, axis=1)
    actual = np.asarray(Image.open(out/'enhanced.png').convert('RGBA'))
    if not np.array_equal(actual, pixels):
        raise ValueError('enhanced image lost native field pixels or colors')
    rate, pcm = wavfile.read(out/'speech.wav')
    if (rate != 48000 or len(pcm)/rate != audio['audio_duration_s']
            or np.max(abs(pcm.astype(np.int32))) >= 32767):
        raise ValueError('invalid or clipped listening waveform')
    with np.load(out/'speech.npz') as native:
        expected = np.rint(32767 * audio['fixed_playback_gain_per_sqrt_W']
                           * native['resampled_power_wave_sqrt_W']).astype(np.int16)
        if not np.array_equal(pcm, expected):
            raise ValueError('listening waveform differs from physical export')
    if (visual['start_sample'] != 0 or visual['end_sample'] != len(pcm)
            or visual['source_sha256'] != sha(out/'speech.wav')):
        raise ValueError('enhanced field does not cover this complete waveform')
    report = dict(status='ready_for_listening', word=manifest['word'],
                  scope=audio['scope'], duration_s=audio['duration_s'],
                  intelligibility='unassessed; requires listening',
                  enhanced_pixels_exact=True, enhanced_shape=list(actual.shape),
                  enhanced_command=command,
                  hashes={str(p.relative_to(ROOT)):sha(p) for p in
                          [Path(__file__), ROOT/'tools/actual_coupled_audio.py', run_path,
                           *sorted(out.iterdir())]})
    save(out/'verification.json', report)
    print(json.dumps(dict(status=report['status'], word=report['word'],
                          audio=str(out/'speech.wav'), enhanced=str(out/'enhanced.png'))), flush=True)


if __name__ == '__main__':
    main()
