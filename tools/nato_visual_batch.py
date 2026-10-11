#!/usr/bin/env python3
"""Render current-sheet NATO words and their complete enhanced fields.

Every take freezes its source, geometry, gestures, and WAV hash. Reference
fields can be reused only when their source hash and full-field settings match.
No time warping, image reduction, output filtering, or reference PCM synthesis.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from scipy.io import wavfile

ROOT = Path(__file__).resolve().parents[1]
WORDS = 'alpha bravo charlie delta echo foxtrot golf hotel india juliet kilo lima mike november oscar papa quebec romeo sierra tango uniform victor whiskey x-ray yankee zulu'.split()
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--design', type=Path)
p.add_argument('--workers', type=int, default=3)
p.add_argument('--words', nargs='+', choices=WORDS, default=WORDS)
p.add_argument('--reference-cache', type=Path)
p.add_argument('--plan-options', type=Path, help='JSON experiment flags merged onto each reference score')
mode = p.add_mutually_exclusive_group()
mode.add_argument('--render-only', action='store_true')
mode.add_argument('--fields-only', action='store_true')
a = p.parse_args()
if not 1 <= a.workers <= 4:
    p.error('Use one to four enhanced-field workers.')
if not a.fields_only and a.design is None:
    p.error('--design is required when rendering speech')
out = a.out.resolve()
out.mkdir(parents=True, exist_ok=True)
manifest_path = out / 'manifest.json'
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')

def save(manifest):
    temp = manifest_path.with_suffix('.tmp')
    temp.write_text(json.dumps(manifest, indent=2) + '\n')
    temp.replace(manifest_path)

if a.fields_only:
    manifest = json.loads(manifest_path.read_text())
else:
    if manifest_path.exists():
        p.error('This take already has a manifest; choose a fresh directory or --fields-only.')
    manifest = []
    for word in a.words:
        i = WORDS.index(word)
        folder = out / word
        folder.mkdir(parents=True, exist_ok=True)
        ref = ROOT / 'samples/derived' / f'nato_{chr(97+i)}.wav'
        cmd = ['node', 'web/tube/render-speech.mjs', word, str(folder), str(a.design.resolve()), 'reference']
        if a.plan_options:
            score = subprocess.check_output(['node', '--input-type=module', '-e',
                "import {referenceScore} from './web/tube/nato-score.mjs'; console.log(JSON.stringify(referenceScore(process.argv[1])));", word], cwd=ROOT)
            plan = json.loads(score)
            plan.update(json.loads(a.plan_options.read_text()))
            plan_path = folder / 'requested-plan.json'
            plan_path.write_text(json.dumps(plan, indent=2))
            cmd += ['--plan', str(plan_path)]
        with (out / (word + '-render.log')).open('w') as log:
            result = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        manifest.append(dict(word=word, letter=chr(65+i), reference=str(ref),
            reference_sha256=hashlib.sha256(ref.read_bytes()).hexdigest(), synthetic=str(folder/'outlet.wav'),
            status='rendered; enhanced pending' if result.returncode == 0 else 'render failed',
            render_ok=result.returncode == 0, recognition='Unverified; visual matching in progress'))
        save(manifest)
        print(word, manifest[-1]['status'], flush=True)
if a.render_only:
    sys.exit(0)

def valid_field(prefix, wav):
    meta = prefix.with_suffix('.json')
    if not meta.exists() or not prefix.with_suffix('.npz').exists():
        return False
    try:
        receipt = json.loads(meta.read_text())
        return (receipt['source_sha256'] == hashlib.sha256(wav.read_bytes()).hexdigest()
                and receipt['hop'] == 128 and receipt['rows'] == 684
                and receipt['start_sample'] == 0 and receipt['channel'] == 0
                and receipt['sample_rate'] == 48000
                and receipt['end_sample'] == len(wavfile.read(wav)[1]))
    except (KeyError, ValueError):
        return False

def field(job):
    i, kind = job
    item = manifest[i]
    folder = out / item['word']
    prefix = folder / kind
    wav = Path(item['reference'] if kind == 'reference' else item['synthetic'])
    if valid_field(prefix, wav):
        return i, kind, 'complete', str(prefix)
    if kind == 'reference' and a.reference_cache:
        cached = a.reference_cache.resolve() / item['word'] / 'reference'
        if valid_field(cached, wav):
            return i, kind, 'complete', str(cached)
    cmd = [sys.executable, str(ROOT/'tools/enhanced_specgram/render.py'), str(wav),
           '--channel', '0', '--max-hz', '8000', '--hop', '128',
           '--display-reference', '.17808298734486477', '--out', str(prefix)]
    with (folder / f'{kind}.log').open('w') as log:
        result = subprocess.run(cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    return i, kind, 'complete' if result.returncode == 0 else 'failed', str(prefix)

jobs = [(i, k) for i, item in enumerate(manifest) if item['word'] in a.words
        for k in ['reference', 'synthetic'] if k == 'reference' or item['render_ok']]
with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
    futures = [pool.submit(field, job) for job in jobs]
    for future in concurrent.futures.as_completed(futures):
        i, kind, result, prefix = future.result()
        item = manifest[i]
        item[kind+'_field'] = result
        item[kind+'_prefix'] = prefix
        if all(item.get(k+'_field') == 'complete' for k in ['reference', 'synthetic']):
            comparison = out / item['word'] / 'comparison.png'
            cmd = [sys.executable, str(ROOT/'tools/compare_nato_fields.py'),
                item['reference_prefix']+'.npz', item['synthetic_prefix']+'.npz',
                '--out', str(comparison), '--gain', '.1614252377244924', '--label', item['word'].title()+' · current sheet take']
            subprocess.run(cmd, cwd=ROOT, env=env, check=True)
            item['comparison'] = str(comparison)
            item['status'] = 'ready for full-field inspection; not a verified match'
        save(manifest)
        print(item['word'], kind, result, flush=True)
