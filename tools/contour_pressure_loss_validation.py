#!/usr/bin/env python3
"""Fresh regressions and sanitizer provenance for the slit loss component."""
import json,subprocess,hashlib
from pathlib import Path
import numpy as np
from contour_pressure_flow_reference import Reference,SCALE
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/slit-kinetic'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(args,name,text=None):
 p=subprocess.run(args,cwd=ROOT,input=text,capture_output=True,text=True);(OUT/(name+'.stdout')).write_text(p.stdout);(OUT/(name+'.stderr')).write_text(p.stderr);assert p.returncode==0,(args,p.stderr,p.stdout);return p.stdout

def main():
 ctest=['ctest','--test-dir','build','-R','^ducklips_contour_pressure_(storage|flow|loss)$','--output-on-failure'];command(ctest,'ctest')
 sanitizer=['ctest','--test-dir','.build/asan','-R','^ducklips_contour_pressure_(storage|flow|loss)$','--output-on-failure'];command(sanitizer,'asan-ctest')
 for mode in ['loss','conductance']:
  args=['--batch'] if mode=='loss' else ['--conductance'];inputs=(OUT/(mode+'-input.txt')).read_text();a=command(['build/ducklips_contour_pressure_loss_probe',*args],mode+'-release',inputs);b=command(['.build/asan/ducklips_contour_pressure_loss_probe',*args],mode+'-asan',inputs);assert a==b
 rows=[r['input'] for r in json.loads((ROOT/'reports/contour_pressure_flow_reference.json').read_text())['cases']];inputs=''.join(' '.join(map(str,r))+'\n' for r in rows);states=[json.loads(s) for s in command(['build/ducklips_contour_pressure_flow_probe','--batch'],'flow-regression',inputs).splitlines()];flow=[]
 for row,state in zip(rows,states):
  ref=Reference(row);m=ref.matrix(1e-11);actual=np.array(state['inertia']).reshape(5,5);error=float(np.max(abs(actual-m)/(ref.base*np.outer(SCALE,SCALE))));field_error=max(float(max(abs(np.array(s['flow'])-ref.flow(s['z']))/SCALE)) for s in state['samples']);assert error<2e-8 and field_error<1e-9;flow.append(dict(input=row,matrix_scaled_error=error,field_scaled_error=field_error))
 old=json.loads((ROOT/'reports/contour_pressure_reference.json').read_text());rows=[r['input'] for r in old['cases_detail']];inputs=(ROOT/'out/contour-pressure/recover-input.txt').read_text();recovered=[json.loads(s) for s in command(['build/ducklips_contour_pressure_probe','--recover'],'inverse-regression',inputs).splitlines()];effort_error=max(max(abs(np.array(state['effort'])-row[2:])) for row,state in zip(rows,recovered));assert effort_error<1e-4
 sweeps=[]
 for n in [100,200,400,800]:
  state=json.loads(command(['build/ducklips_contour_pressure_probe','--sweep',str(n)],f'sweep-{n}'));assert state['maximum_mass_error_m3']<1e-20 and state['maximum_relative_density']<1.32;sweeps.append(state)
 ratios=[abs(a['work_residual_J']/b['work_residual_J']) for a,b in zip(sweeps,sweeps[1:])];assert all(3.5<r<4.5 for r in ratios[-2:]) and abs(sweeps[-1]['work_residual_J'])<1e-14
 paths=[Path(__file__),ROOT/'CMakeLists.txt',*ROOT.glob('include/ducklips/contour_pressure_*.hpp'),*ROOT.glob('src/contour_pressure_*.cpp'),ROOT/'src/slit_bgk_table.inc',ROOT/'tools/slit_kinetic_reference.py',ROOT/'tools/contour_pressure_loss_probe.cpp',ROOT/'tools/contour_pressure_loss_reference.py',ROOT/'reports/slit_kinetic_reference.json',ROOT/'reports/slit_kinetic_initial_grid_failure.json',ROOT/'reports/contour_pressure_loss_reference.json',*OUT.glob('*.stdout'),*OUT.glob('*.stderr'),*OUT.glob('*-input.txt'),*OUT.glob('*-output.txt')]
 for folder in ['build','.build/asan']:
  paths.extend(ROOT/folder/name for name in ['ducklips_contour_pressure_probe','ducklips_contour_pressure_flow_probe','ducklips_contour_pressure_loss_probe'])
 report=dict(scope='instantaneous constitutive component, no momentum integration or new speech',base_commit='50d1a2c3a69aa4f0afcd53739bae5542598a6b14',ctest_command=ctest,focused_ctests_passed=3,sanitizer_command=sanitizer,sanitizer_ctests_passed=3,sanitizers='address,undefined',sanitizer_batch_and_conductance_byte_identical=True,flow_regression=flow,inverse_maximum_effort_error_Pa=float(effort_error),prescribed_geometry_sweeps=sweeps,work_refinement_ratios=ratios,hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))})
 (ROOT/'reports/contour_pressure_loss_validation.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS three focused tests and sanitizers, identical sanitized outputs, flow/inverse/seam regressions')
if __name__=='__main__':main()
