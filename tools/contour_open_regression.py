#!/usr/bin/env python3
"""Check current open motion against the preserved independent open oracle.

The oracle is reused from its committed receipt, not integrated again here.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    commit='5bc0d1ace9273b158346921fd1a652edbb8c508f';reference=ROOT/'reports/contour_glottal_reference.json'
    original=subprocess.check_output(['git','show',commit+':reports/contour_glottal_reference.json'],cwd=ROOT);assert original==reference.read_bytes()
    # The independent equations, gas law and default body/air must match the
    # original oracle; changed source solvers are what this regression tests.
    for path in ['tools/contour_glottal_reference.py','src/gas_storage.cpp','include/ducklips/physics.hpp','include/ducklips/sources.hpp']:
        assert subprocess.check_output(['git','show',commit+':'+path],cwd=ROOT)==(ROOT/path).read_bytes(),path
    prior=json.loads(original);binary=ROOT/'build/ducklips_contour_glottal_probe';out=ROOT/'out/contour-sealed';out.mkdir(parents=True,exist_ok=True)
    scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3,[700.]*3,[1e-4]*4];c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004,.002]);runs=[];acceptance=[]
    for mode in ['forward','reverse','relaxation','rest_motion']:
        oracle=next(x for x in prior['oracles'] if x['mode']==mode);expected=np.array(oracle['final_state']);previous=None
        for steps in ([3840,7680,15360] if mode=='relaxation' else [960,1920,3840]):
            command=[str(binary),'--trajectory',str(steps),str(prior['duration_s']),mode];run=subprocess.run(command,capture_output=True,text=True);assert run.returncode==0,run.stderr
            path=out/f'open-{mode}-{steps}.json';path.write_text(run.stdout);s=json.loads(run.stdout)
            actual=np.r_[s['mass'],s['position'],s['velocity'],s['passage_momentum'],np.array(s['charge'])/c,g*np.array(s['network_momentum'])];error=abs(actual-expected);norm=float(max(error/scale))
            runs.append(dict(mode=mode,steps=steps,normalized_error=norm,refinement_ratio=None if previous is None else previous/norm,maximum_network_pressure_error_Pa=float(max(error[13:16])),maximum_work_balance_J=s['maximum_work_balance'],command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path)));previous=norm
        acceptance.append(dict(mode=mode,second_order=all(3.5<r['refinement_ratio']<4.5 for r in runs[-2:]),endpoint_accuracy=previous<1e-4))
    sources=['src/contour_kinetics.cpp','src/contour_glottal_storage.cpp','src/fold_contact.cpp','src/glottal_contour.cpp','src/gas_storage.cpp','src/moving_graph.cpp','include/ducklips/contour_kinetics.hpp','include/ducklips/contour_glottal_storage.hpp','include/ducklips/fold_contact.hpp','tools/contour_open_regression.py','tools/contour_glottal_probe.cpp']
    report=dict(scope='current open dynamics compared with reused independent open oracle',oracle_commit=commit,reference=str(reference.relative_to(ROOT)),reference_sha256=sha(reference),runs=runs,acceptance=acceptance,hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),binary]})
    (ROOT/'reports/contour_open_regression.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(acceptance));assert all(x['second_order'] and x['endpoint_accuracy'] for x in acceptance)
if __name__=='__main__':main()
