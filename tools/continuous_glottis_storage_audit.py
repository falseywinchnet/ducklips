#!/usr/bin/env python3
"""Fresh 80-digit storage and near-ambient effort audit after coupling."""
import json,subprocess,hashlib
from pathlib import Path
import mpmath as mp
from contour_pressure_reference import oracle,unpack
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-glottis'
def main():
 rows=[r['input'] for r in json.loads((ROOT/'reports/contour_pressure_reference.json').read_text())['cases_detail']]
 rows += [[3e-6,3e-6,u,-u] for u in [1e-12,1e-9,1e-6,.001,1.,10.]]
 inp=''.join(' '.join(map(str,r))+'\n' for r in rows);(OUT/'storage-input.txt').write_text(inp);binary=ROOT/'build/ducklips_contour_pressure_probe';p=subprocess.run([str(binary),'--batch'],input=inp,capture_output=True,text=True);assert p.returncode==0,p.stderr;(OUT/'storage-output.jsonl').write_text(p.stdout);states=[json.loads(s) for s in p.stdout.splitlines()];mp.mp.dps=80;records=[]
 for i,(row,state) in enumerate(zip(rows,states)):
  expected=oracle(row);v=float(expected[0]);scale=[v,v,v,139650*v,139650*v,139650*.003,139650*.003,v/139650,v/139650,v/139650,.003,.003,.003];error=float(max(abs(mp.mpf(a)-b)/c for a,b,c in zip(unpack(state),expected,scale)));energy_relative=float(abs(mp.mpf(state['energy'])/expected[3]-1)) if expected[3] else 0.
  records.append(dict(input=row,maximum_scaled_error=error,energy_relative_error=energy_relative));assert error<1e-10
  if i>=72:assert energy_relative<1e-12,(row,energy_relative)
  if i%12==11:print('checked',i+1,flush=True)
 recover='3e-6 3e-6 6.0002e-9 3e-9\n';p=subprocess.run([str(binary),'--recover'],input=recover,capture_output=True,text=True);assert p.returncode==0,p.stderr;recovered=json.loads(p.stdout)
 paths=[binary,Path(__file__),ROOT/'src/gas_storage.cpp',ROOT/'include/ducklips/gas_storage.hpp',ROOT/'src/contour_pressure_storage.cpp',OUT/'storage-input.txt',OUT/'storage-output.jsonl']
 report=dict(scope='fresh gas-storage regression including signed effort as small as 1e-12 Pa',precision_digits=80,cases=records,maximum_scaled_error=max(r['maximum_scaled_error'] for r in records),maximum_near_ambient_energy_relative_error=max(r['energy_relative_error'] for r in records[72:]),previously_failing_recovery=dict(input=recover,result=recovered),hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 (ROOT/'reports/continuous_glottis_storage_audit.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS 78 independent storage cases and near-ambient inverse recovery')
if __name__=='__main__':main()
