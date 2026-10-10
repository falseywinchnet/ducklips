#!/usr/bin/env python3
"""Compare saved native prefixes with independently refined open ODE motion."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'out/continuous-tract-contact'
NAMES = ['local-6400', 'local-12800', 'fitted-3200', 'fitted-6400',
    'fitted-12800', 'fitted-25600', 'fitted-51200']

def main():
    oracle = []
    paths = [Path(__file__)]
    for name in ['open-96', 'open-192']:
        receipt = json.loads((OUT/(name+'.json')).read_text())
        for namepath, digest in receipt['hashes'].items():
            assert hashlib.sha256((ROOT/namepath).read_bytes()).hexdigest() == digest, namepath
        oracle.append(np.load(OUT/(name+'.npz')))
        paths += [OUT/(name+'.json'), OUT/(name+'.npz')]
    a, b = oracle
    assert np.array_equal(a['times'], b['times'])
    ref_error = np.max(abs(a['scaled_states']-b['scaled_states']), axis=0)
    assert max(ref_error) < 1e-7
    scale = b['scale']
    records, runs = {}, {}
    for name in NAMES:
        # Save an immutable prefix, even when a longer native run is still
        # in flight. A later complete row proves this interval is finished.
        lines = (OUT/(name+'.jsonl')).read_text().splitlines()
        rows = [json.loads(line) for line in lines if line.endswith('}')]
        assert rows[-1]['time'] > b['times'][-1]
        rows = [r for r in rows if r['time'] <= b['times'][-1]+1e-15]
        prefix = OUT/(name+'-open.jsonl')
        prefix.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        paths.append(prefix)
        times = np.array([r['time'] for r in rows])
        indices = np.rint(times/(b['times'][1]-b['times'][0])).astype(int)
        assert np.max(abs(times-b['times'][indices])) < 1e-14
        values = np.array([r['source']+r['charge']+r['momentum'] for r in rows])/scale
        assert np.max(abs(values[0]-b['scaled_states'][0])) < 1e-8
        errors = values-b['scaled_states'][indices]
        runs[name] = {round(t, 12): e for t, e in zip(times, errors)}
        records[name] = dict(samples=len(rows), final_time_s=float(times[-1]),
            maximum_scaled_error=float(np.max(abs(errors))),
            maximum_scaled_error_by_component=np.max(abs(errors), axis=0).tolist(),
            final_scaled_error=errors[-1].tolist(),
            maximum_tract_pressure_error_Pa=float(np.max(abs(errors[:,9:17]))*1e-8/(2e-4*.175/8/139650.)))
    selected = ['fitted-12800', 'fitted-25600', 'fitted-51200']
    common = sorted(set.intersection(*(set(runs[n]) for n in selected)))
    common_errors = [float(np.max(abs(np.array([runs[n][t] for t in common])))) for n in selected]
    ratios = [x/y for x,y in zip(common_errors, common_errors[1:])]
    assert all(3.7 < r < 4.5 for r in ratios), ratios
    paths += [OUT/'fitted-probe', OUT/'seed.input']
    report = dict(scope='Open motion only, from absolute time 8 ms until 8.425 ms; before first seal. No claim about later contact convergence or speech.',
        oracle_refinement_maximum_scaled_error=float(max(ref_error)),
        oracle_refinement_by_component=ref_error.tolist(),
        normalization=scale.tolist(), native=records,
        common_grid=dict(cases=selected, samples=len(common), times_s=common,
            maximum_scaled_errors=common_errors, refinement_ratios=ratios),
        sampling_note='Per-run maxima use different output grids. Ratios above use the same sample times without interpolation.',
        native_provenance='Local controls use archived 3a9eb8b control-probe; fitted runs use the relaxation implementation committed in 224a4b1. Long-run completion is not claimed by prefix checks.',
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/continuous_tract_open_reference.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(oracle_error=max(ref_error), common_grid=report['common_grid']), indent=2))

if __name__ == '__main__':
    main()
