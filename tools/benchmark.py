#!/usr/bin/env python3
"""Repeatable 10x realtime gate. Rendering, startup and file I/O are reported separately."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
TEXT = 'Alfa. Bravo. Charlie. Delta. Echo. Foxtrot. Golf. Hotel. India. Juliett. Kilo. Lima. Mike. November. Oscar. Papa. Quebec. Romeo. Sierra. Tango. Uniform. Victor. Whiskey. Xray. Yankee. Zulu.'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, default=ROOT / 'build/ducklips_say')
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--out', type=Path, default=ROOT / 'reports/performance.json')
    p.add_argument('--minimum-speed', type=float, default=10.0)
    args = p.parse_args()
    if args.repeats < 1 or args.minimum_speed <= 0:
        p.error('repeats and minimum-speed must be positive')
    output = ROOT / 'out/benchmark'
    output.mkdir(parents=True, exist_ok=True)
    cases = [(v, ['--voice', v]) for v in ['baritone','alto','child','gravel','whisper']]
    cases += [('body', ['--body', str(ROOT / 'bodies/reference.body')]),
              ('fast', ['--voice','baritone','--rate','2'])]
    rows = []
    for name, options in cases:
        runs = []
        for repeat in range(args.repeats):
            command = [str(args.binary.resolve()), *options, '--dictionary', str(ROOT / 'samples/nato.dict'),
                       '--text', TEXT, '--out', str(output / (name + '.wav'))]
            start = time.perf_counter()
            result = subprocess.run(command, text=True, capture_output=True, check=True, timeout=180)
            wall = time.perf_counter() - start
            (output / f'{name}-{repeat}.log').write_text(result.stdout + result.stderr)
            match = re.search(r'rendered ([\d.]+) s in ([\d.]+) s', result.stdout)
            if not match:
                raise RuntimeError('Missing timing report')
            duration, rendering = map(float, match.groups())
            audit = re.search(r'output audit: nonfinite (\d+), clipped (\d+)', result.stdout)
            if not audit:
                raise RuntimeError('Missing pre-clamp output audit')
            prepare = re.search(r'prepared in ([\d.]+) s', result.stdout)
            runs.append(dict(audio_seconds=duration,render_seconds=rendering,wall_seconds=wall,
                             prepare_seconds=float(prepare[1]), speed=duration/rendering,
                             end_to_end_speed=duration/wall, nonfinite=int(audit[1]),clipped=int(audit[2])))
        minimum = min(r['speed'] for r in runs)
        passed = minimum >= args.minimum_speed and all(r['nonfinite']==0 and r['clipped']==0 for r in runs)
        row = dict(case=name,command=command,runs=runs,minimum_speed=minimum,
                   median_speed=statistics.median(r['speed'] for r in runs),passed=passed)
        rows.append(row)
        print(f'{name}: minimum {minimum:.2f}x, median {row["median_speed"]:.2f}x, {"PASS" if passed else "FAIL"}', flush=True)
    report = dict(utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(), machine=platform.machine(),
                  cpu=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip() if platform.system()=='Darwin' else platform.processor(),
                  binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                  git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  working_tree_changes=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).splitlines(),
                  criterion='Every timed render >= minimum_speed and no nonfinite or clipped samples before output limiting',
                  minimum_speed=args.minimum_speed,includes_silent_intervals=True,
                  includes_prepare_in_render_time=False, excludes='Disk I/O and offline enhanced-field visualization',
                  cases=rows,passed=all(r['passed'] for r in rows))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0 if report['passed'] else 1)


if __name__=='__main__':
    main()
