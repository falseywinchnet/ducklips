#!/usr/bin/env python3
"""One complete speech round: fresh audio, full enhanced field, and provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.io import wavfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--text', required=True)
    p.add_argument('--body', type=Path, default=ROOT/'bodies/reference.body')
    p.add_argument('--dictionary', type=Path, default=ROOT/'samples/nato.dict')
    p.add_argument('--rate', type=float, default=1)
    p.add_argument('--tap', choices=['microphones', 'lips', 'source'], default='microphones',
                   help='Observe the same engine at a diagnostic output; does not disconnect feedback')
    p.add_argument('--playback-gain', type=float,
                   help='Use a declared common gain for controlled comparisons instead of peak normalization')
    p.add_argument('--level-match', type=Path,
                   help='Match one constant playback gain to this reference WAV in the declared windows')
    p.add_argument('--level-window', type=float, nargs=2, action='append', metavar=('START', 'END'))
    p.add_argument('--hop', type=int, default=128)
    p.add_argument('--display-reference', type=float, default=.17808298734486477)
    p.add_argument('--reference-field', type=Path, help='Full enhanced human reference NPZ')
    p.add_argument('--viewer', type=Path, help='Comparison view; requires --reference-field')
    args = p.parse_args()
    if args.playback_gain is not None and (not np.isfinite(args.playback_gain) or args.playback_gain <= 0):
        p.error('--playback-gain must be finite and positive')
    if args.viewer and not args.reference_field:
        p.error('--viewer requires --reference-field')
    if args.level_match and (not args.level_window or args.playback_gain is None):
        p.error('--level-match requires --level-window and the reference --playback-gain')
    if args.level_window and not args.level_match:
        p.error('--level-window requires --level-match')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for source, name in [(ROOT/'build/ducklips_say', 'ducklips_say'),
                         (args.body, 'body.body'), (args.dictionary, 'nato.dict')]:
        shutil.copy2(source, out/name)
    sources = [ROOT/'CMakeLists.txt', *sorted((ROOT/'src').glob('*.cpp')),
               *sorted((ROOT/'include/ducklips').glob('*.hpp')), ROOT/'tools/ducklips_say.cpp', Path(__file__)]
    for source in sources:
        dest = out/'source-snapshot'/source.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), text=args.text,
                  engine='local-flow; conservative tongue; exact source snapshot retained',
                  status='rendering', commands=[],
                  files={str(x.relative_to(out)):digest(x) for x in out.rglob('*') if x.is_file()})

    def save():
        (out/'round.json').write_text(json.dumps(record, indent=2)+'\n')

    def run(command, log):
        record['commands'].append(command)
        save()
        with (out/log).open('w') as f:
            subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True, cwd=ROOT)

    try:
        run([str(out/'ducklips_say'), '--local-flow', '--conservative-tongue',
             '--body', str(out/'body.body'), '--dictionary', str(out/'nato.dict'),
             '--text', args.text, '--rate', str(args.rate), '--tap', args.tap, '--full-scale-pa', '128',
             '--trace', '0', str(out/'trace.csv'), '--probe', str(out/'probe.csv'),
             '--report', '--out', str(out/'speech-raw.wav')], 'speech.log')
        fs, raw = wavfile.read(out/'speech-raw.wav')
        peak = float(np.max(abs(raw)))
        if not np.isfinite(raw).all() or not 0 < peak < 1:
            raise ValueError('Native speech is invalid, silent or clipped')
        gain = args.playback_gain if args.playback_gain is not None else .8/peak
        if args.level_match:
            reference_fs, reference = wavfile.read(args.level_match)
            if reference_fs != fs or reference.shape != raw.shape or reference.dtype != np.int16:
                raise ValueError('Level reference must be a same-shape, same-rate PCM16 recording')
            selected = np.zeros(len(raw), dtype=bool)
            for start, end in args.level_window:
                if not np.isfinite([start, end]).all() or not 0 <= start < end <= len(raw)/fs:
                    raise ValueError('Level windows must be finite, ordered and inside the recording')
                selected[round(start*fs):round(end*fs)] = True
            reference_rms = float(np.sqrt(np.mean((reference[selected].astype(float)/32767)**2)))
            candidate_rms = float(np.sqrt(np.mean(raw[selected].astype(float)**2)))
            if not reference_rms > 0 or not candidate_rms > 0:
                raise ValueError('Level windows must contain nonzero audio')
            gain = reference_rms/candidate_rms
            args.display_reference *= gain/args.playback_gain
            record['level_match'] = dict(reference=str(args.level_match.resolve()),
                reference_sha256=digest(args.level_match), windows_seconds=args.level_window,
                reference_rms=reference_rms, candidate_raw_rms=candidate_rms,
                reference_playback_gain=args.playback_gain, display_reference=args.display_reference,
                method='one stereo RMS gain over union of windows; not perceptual loudness matching')
        if peak*gain > 1:
            raise ValueError('Declared playback gain would clip; choose a lower common gain')
        wavfile.write(out/'speech.wav', fs, np.rint(raw*gain*32767).astype(np.int16))
        record['playback'] = dict(fixed_gain=gain, gain_selection='voiced RMS match' if args.level_match else ('declared' if args.playback_gain is not None else 'peak 0.8'),
                                  tap=args.tap, frames=len(raw), sample_rate=fs,
                                  processing='one fixed gain over the full recording; no cuts or filters')
        record['status'] = 'audio rendered; enhanced field pending'
        run([sys.executable, str(ROOT/'tools/enhanced_specgram/render.py'), str(out/'speech.wav'),
             '--channel', '0', '--max-hz', '8000', '--hop', str(args.hop),
             '--display-reference', str(args.display_reference), '--out', str(out/'enhanced')], 'enhanced.log')
        receipt = json.loads((out/'enhanced.json').read_text())
        if receipt['source_sha256'] != digest(out/'speech.wav'):
            raise ValueError('Depiction does not belong to this audio')
        record['enhanced'] = receipt
        if args.viewer:
            run([sys.executable, str(ROOT/'tools/native_fields.py'), str(args.reference_field.resolve()),
                 str(out/'enhanced.npz'), '--out', str(args.viewer.resolve()),
                 '--human-label', 'Human reference',
                 '--generated-label', 'Fresh synthesis: '+args.text], 'viewer.log')
            viewer = args.viewer.resolve()
            record['viewer'] = dict(path=str(viewer), sha256=digest(viewer))
        record['status'] = 'audio and full enhanced field ready; visual inspection and listening assessment required'
    except Exception as error:
        record['status'] = 'incomplete'
        record['error'] = str(error)
        save()
        raise
    record['files'].update({x.name:digest(x) for x in out.iterdir() if x.is_file() and x.name != 'round.json'})
    save()
    print(json.dumps({'audio': str(out/'speech.wav'), 'depiction': str(out/'enhanced.png'),
                      'status':record['status']}, indent=2))


if __name__ == '__main__':
    main()
