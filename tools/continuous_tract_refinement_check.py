#!/usr/bin/env python3
"""Collect completed fine contact runs; do not infer completion from prefixes."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'out/continuous-tract-contact'
SCALE = np.array([1000.,1000.,1e-6,1e-6,1e-4,1e-4,1e-4,.001,.001])

def main():
    cases = ['fitted-6400','fitted-12800','fitted-25600','fitted-51200']
    runs, rows, paths = {}, {}, [Path(__file__), OUT/'fitted-probe', OUT/'seed.input']
    original = json.loads((ROOT/'reports/continuous_tract_contact_reference.json').read_text())
    for name in cases:
        n = int(name.split('-')[1])
        path = OUT/(name+'.jsonl')
        data = list(map(json.loads, path.read_text().splitlines()))
        summary = json.loads((OUT/(name+'.stderr')).read_text())
        assert summary['complete'] and data[-1]['step'] == n
        assert summary['maximum_work_residual_J'] < 1e-12
        assert summary['maximum_mass_residual_m3'] < 1e-19
        for key in ['source','charge','momentum']:
            assert np.allclose(data[0][key],original['seed'][key],rtol=1e-13,atol=1e-20)
        cap = {6400:7.8125001e-8,12800:7.8125001e-8,25600:3.9062501e-8,51200:1.9531251e-8}[n]
        runs[name] = dict(summary=summary, final=data[-1], maximum_step_s=.016/n, contact_step_cap_s=cap,
            reproduction_command=['out/continuous-tract-contact/fitted-probe',str(n),'.016','--graph-state','--relaxation','--contact-step',str(cap)])
        rows[name] = {round(r['time'],12):r for r in data}
        paths += [path,OUT/(name+'.stderr')]
    comparisons = []
    for a,b in zip(cases,cases[1:]):
        common = sorted(rows[a].keys() & rows[b].keys())
        windows = []
        for end in [.000425,.006,.016]:
            selected = [t for t in common if t<=end]
            x = np.array([rows[a][t]['source'] for t in selected])
            y = np.array([rows[b][t]['source'] for t in selected])
            qx = np.array([rows[a][t]['charge'] for t in selected])
            qy = np.array([rows[b][t]['charge'] for t in selected])
            d = x-y
            windows.append(dict(end_s=end,last_sample_s=selected[-1],samples=len(selected),
                maximum_scaled_source_difference=float(np.max(abs(d/SCALE))),
                final_scaled_source_difference=float(np.max(abs(d[-1]/SCALE))),
                maximum_area_difference_mm2=float(np.max(abs(d[:,2:4]))*1e6),
                final_area_difference_mm2=float(np.max(abs(d[-1,2:4]))*1e6),
                maximum_tract_pressure_difference_Pa=float(np.max(abs(qx-qy))/(2e-4*.175/8/139650.)),
                final_tract_pressure_difference_Pa=float(np.max(abs(qx[-1]-qy[-1]))/(2e-4*.175/8/139650.))))
        comparisons.append(dict(cases=[a,b],windows=windows))
    report = dict(scope='Common-state 16 ms contact refinement with both open step and crossing cap refined on the finest pair. Adjacent-run differences are not errors against an independent contact oracle.',
        original_source_commit='224a4b12ce883fd52b89e7f93315e9a58e2a1c0e',
        native=runs,comparisons=comparisons,
        fine_contact_temporal_convergence_established=False,
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/continuous_tract_refinement.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(comparisons,indent=2))

if __name__ == '__main__':
    main()
