#!/usr/bin/env python3
"""Refine the relaxation trajectory against the already verified ODE endpoint.

This explicitly reuses the continuous oracle; it is not a fresh ODE solve.
The original report retains its failed coarse-grid convergence criterion.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    reference=ROOT/'reports/contour_glottal_reference.json';prior=json.loads(reference.read_text())
    for name,digest in prior['hashes'].items():assert sha(ROOT/name)==digest,name
    out=ROOT/'out/contour-glottal';binary=ROOT/'build/ducklips_contour_glottal_probe'
    oracle=next(x for x in prior['oracles'] if x['mode']=='relaxation');expected=np.array(oracle['final_state'])
    previous=next(x['normalized_error'] for x in prior['runs'] if x['mode']=='relaxation' and x['steps']==3840)
    scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3,[700.]*3,[1e-4]*4];c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004,.002]);runs=[]
    for steps in [7680,15360]:
        command=[str(binary),'--trajectory',str(steps),str(prior['duration_s']),'relaxation'];run=subprocess.run(command,capture_output=True,text=True);assert run.returncode==0,run.stderr
        path=out/f'relaxation-{steps}.json';path.write_text(run.stdout);native=json.loads(run.stdout)
        actual=np.r_[native['mass'],native['position'],native['velocity'],native['passage_momentum'],np.array(native['charge'])/c,g*np.array(native['network_momentum'])]
        error=abs(actual-expected);norm=float(max(error/scale))
        runs.append(dict(steps=steps,normalized_error=norm,refinement_ratio=previous/norm,maximum_network_pressure_error_Pa=float(max(error[13:16])),maximum_position_error_m=float(max(error[2:6])),maximum_mass_balance_m3=native['maximum_mass_balance'],maximum_work_balance_J=native['maximum_work_balance'],command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path)));previous=norm
    acceptance=[a for a in prior['acceptance'] if a['mode']!='relaxation']+[dict(mode='relaxation',second_order=all(3.5<r['refinement_ratio']<4.5 for r in runs),endpoint_accuracy=previous<1e-4)]
    report=dict(scope='two additional native refinements compared with reused independent continuous relaxation endpoint',reference=str(reference.relative_to(ROOT)),reference_sha256=sha(reference),oracle_refinement_difference=oracle['normalized_refinement_difference'],runs=runs,combined_acceptance=acceptance,hashes={str(p.relative_to(ROOT)):sha(p) for p in [binary,Path(__file__)]})
    path=ROOT/'reports/contour_glottal_endpoint_reference.json';path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report['combined_acceptance']));assert all(x['second_order'] and x['endpoint_accuracy'] for x in acceptance)
if __name__=='__main__':main()
