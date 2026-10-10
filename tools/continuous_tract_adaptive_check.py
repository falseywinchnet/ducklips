#!/usr/bin/env python3
"""Check the opt-in step-doubling control against an independent open ODE."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'out/continuous-tract-contact'

def main():
    oracle_receipt = json.loads((OUT/'open-192.json').read_text())
    for path, digest in oracle_receipt['hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    oracle = np.load(OUT/'open-192.npz')
    records, paths = [], [Path(__file__), ROOT/'tools/continuous_tract_probe.cpp',
        ROOT/'build/ducklips_continuous_tract_probe', ROOT/'.build/asan/ducklips_continuous_tract_probe',
        OUT/'open-192.npz', OUT/'open-192.json', OUT/'seed.input', OUT/'adaptive-contact.input']
    for label, tolerance in [('1e4',1e-4),('25e6',2.5e-5),('625e8',6.25e-6)]:
        name = 'checked-adaptive-open-'+label
        rows = list(map(json.loads,(OUT/(name+'.jsonl')).read_text().splitlines()))
        summary = json.loads((OUT/(name+'.stderr')).read_text())
        assert summary['complete'] and rows[-1]['step']==20 and len(rows)==21
        times = np.array([r['time'] for r in rows])
        index = np.rint(times/(oracle['times'][1]-oracle['times'][0])).astype(int)
        assert np.max(abs(times-oracle['times'][index])) < 1e-14
        values = np.array([r['source']+r['charge']+r['momentum'] for r in rows])/oracle['scale']
        error = values-oracle['scaled_states'][index]
        assert summary['maximum_work_residual_J'] < 1e-12
        assert summary['maximum_mass_residual_m3'] < 1e-19
        records.append(dict(tolerance=tolerance,summary=summary,
            maximum_scaled_error=float(np.max(abs(error))),
            final_scaled_error=float(np.max(abs(error[-1]))),
            maximum_scaled_error_by_component=np.max(abs(error),axis=0).tolist(),
            reproduction_command=['build/ducklips_continuous_tract_probe','20','.0001','--graph-state','--relaxation','--tolerance',str(tolerance),'--contact-step','1.9531251e-8']))
        paths += [OUT/(name+'.jsonl'),OUT/(name+'.stderr')]
    assert records[2]['maximum_scaled_error'] < records[1]['maximum_scaled_error'] < records[0]['maximum_scaled_error'] < 1e-3
    contact = json.loads((OUT/'checked-adaptive-contact.stderr').read_text())
    assert contact['complete'] and contact['closures']==1 and contact['reopenings']==0
    assert contact['maximum_work_residual_J'] < 1e-12 and contact['maximum_mass_residual_m3'] < 1e-19
    paths += [OUT/'checked-adaptive-contact.jsonl',OUT/'checked-adaptive-contact.stderr',
        OUT/'force-asan.jsonl',OUT/'adaptive-asan.jsonl',OUT/'adaptive-asan.stderr']
    report = dict(scope='Fresh short smooth controls against the independent ODE at three tolerances, plus one closure under ASAN/UBSAN. Local tolerance is not a global trajectory error bound.',
        duration_s=.0001,maximum_requested_step_s=5e-6,open=records,
        sanitizer_contact=contact,
        contact_reproduction_command=['.build/asan/ducklips_continuous_tract_probe','20','.0001','--graph-state','--relaxation','--tolerance','1e-4','--contact-step','1.9531251e-8'],
        unclaimed='The 425 microsecond compression interval and full 16 ms contact evolution are not validated by these short controls. Earlier recursively nested prototype jobs were stopped, not counted as passes.',
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (ROOT/'reports/continuous_tract_adaptive_reference.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(errors=[r['maximum_scaled_error'] for r in records],segments=[r['summary']['segments'] for r in records],contact=contact),indent=2))

if __name__ == '__main__':
    main()
