#!/usr/bin/env python3
"""Run pressure-driven tip contact with a coupled, frozen acoustic tube."""
import argparse
import concurrent.futures
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
from continuous_phonation_reference import YSCALE

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--collect', action='store_true', help='collect completed local runs instead of rerunning them')
    args = parser.parse_args()
    out = ROOT/'out/graph-continuous-glottis'
    out.mkdir(parents=True, exist_ok=True)
    seed_path = ROOT/'reports/continuous_phonation_motion.json'
    seed = json.loads(seed_path.read_text())['seeds']['contact_state']
    seed_text = ' '.join(map(str, seed))+'\n'
    (out/'tube-seed.input').write_text(seed_text)
    binary = ROOT/'build/ducklips_continuous_tract_probe'
    counts = [2400, 4800, 9600, 19200]
    def run(n):
        result = subprocess.run([str(binary), str(n), '.024'], input=seed_text, capture_output=True, text=True)
        (out/f'tube-{n}.jsonl').write_text(result.stdout)
        (out/f'tube-{n}.stderr').write_text(result.stderr)
        assert result.returncode == 0, result.stderr
    if not args.collect:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(run, counts))
    runs = []
    c = np.full(8, 2e-4*(.175/8)/139650.)
    g = np.full(8, 2e-4/(1.14*(.175/8)))
    g[-1] *= 2
    for n in counts:
        receipt = json.loads((out/f'tube-{n}.stderr').read_text())
        rows = [json.loads(line) for line in (out/f'tube-{n}.jsonl').read_text().splitlines()]
        assert receipt['complete'] and rows[-1]['step'] == n
        assert receipt['closures'] > 0 and receipt['closures'] == receipt['reopenings']
        assert receipt['maximum_work_residual_J'] < 1e-12 and receipt['maximum_mass_residual_m3'] < 1e-19
        assert np.max(np.abs((np.array(rows[0]['source'])-seed)/YSCALE)) < 1e-8
        receipt['final'] = rows[-1]
        runs.append(receipt)
    differences = []
    for a, b in zip(runs, runs[1:]):
        x, y = a['final'], b['final']
        differences.append(dict(steps=[a['steps'], b['steps']],
            source_maximum_scaled_difference=float(np.max(np.abs((np.array(x['source'])-y['source'])/YSCALE))),
            tract_maximum_pressure_difference_Pa=float(np.max(np.abs(np.array(x['charge'])-y['charge'])/c)),
            tract_maximum_flow_difference_m3_s=float(np.max(np.abs(np.array(x['momentum'])-y['momentum'])*g))))
    paths = [Path(__file__), seed_path, ROOT/'tools/continuous_tract_probe.cpp', binary,
        ROOT/'src/graph_continuous_glottis.cpp', ROOT/'include/ducklips/graph_continuous_glottis.hpp',
        ROOT/'src/continuous_glottis.cpp', ROOT/'include/ducklips/continuous_glottis.hpp',
        ROOT/'src/moving_graph.cpp', out/'tube-seed.input',
        *[out/f'tube-{n}.{ext}' for n in counts for ext in ['jsonl', 'stderr']]]
    report = dict(scope='Pressure-driven source plus frozen eight-cell ideal tube, transient contact and feedback; not speech, anatomy calibration, radiation or settled phonation',
        seed=seed, pressure_Pa=1400., tube=dict(length_m=.175, area_m2=2e-4, nodes=8, boundary='ideal zero effort',
            initial_charge='zero', initial_flow='source initial outlet flow in every branch',
            viscosity=False, radiation_loss=False),
        method='symmetric source step except explicit port-constraint changes; rejected steps bisect both source and graph transactionally',
        temporal_convergence_established=False,
        convergence_note='Balance checks establish conservation, not trajectory accuracy. Phase counts and adjacent endpoint differences remain part of the accuracy assessment.',
        native=runs, adjacent_endpoint_differences=differences,
        hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/continuous_tract_motion.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(runs=[{k:v for k,v in r.items() if k!='final'} for r in runs], endpoint_differences=differences), indent=2))

if __name__ == '__main__':
    main()
