#!/usr/bin/env python3
"""Additional lower-contact refinement using the preserved continuous oracle."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    reference=ROOT/'reports/contour_sealed_reference.json';prior=json.loads(reference.read_text())
    for name,digest in prior['hashes'].items():assert sha(ROOT/name)==digest,name
    binary=ROOT/'build/ducklips_contour_sealed_probe';out=ROOT/'out/contour-sealed';scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3];runs=[];acceptance=[]
    for mode in ['lower','lower_rest']:
        previous=next(r for r in prior['runs'] if r['mode']==mode and r['steps']==6400);oracle=next(r for r in prior['oracles'] if r['mode']==mode);command=[str(binary),'--trajectory','12800',mode]
        run=subprocess.run(command,capture_output=True,text=True);assert run.returncode==0,run.stderr
        path=out/f'{mode}-12800.json';path.write_text(run.stdout);native=json.loads(run.stdout);actual=np.r_[native['mass'],native['position'],native['velocity'],native['momentum']];error=abs(actual-oracle['endpoint']);norm=float(max(error/scale));ratio=previous['normalized_error']/norm
        runs.append(dict(mode=mode,steps=12800,normalized_error=norm,refinement_ratio=ratio,previous_refinement_ratio=previous['refinement_ratio'],maximum_position_error_m=float(max(error[2:6])),maximum_work_balance_J=native['maximum_work_balance'],maximum_mass_balance_m3=native['maximum_mass_balance'],command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path)))
        acceptance.append(dict(mode=mode,second_order=3.5<ratio<4.5 and 3.5<previous['refinement_ratio']<4.5,endpoint_accuracy=norm<1e-4))
    acceptance += [r for r in prior['acceptance'] if r['mode'] not in ['lower','lower_rest']]
    report=dict(scope='additional native lower-contact refinements; continuous oracles reused from the verified initial report',reference=str(reference.relative_to(ROOT)),reference_sha256=sha(reference),runs=runs,combined_acceptance=acceptance,hashes={str(p.relative_to(ROOT)):sha(p) for p in [binary,Path(__file__)]})
    (ROOT/'reports/contour_sealed_endpoint_reference.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(acceptance));assert all(r['second_order'] and r['endpoint_accuracy'] for r in acceptance)
if __name__=='__main__':main()
