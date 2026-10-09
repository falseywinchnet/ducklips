#!/usr/bin/env python3
"""Focused source, sanitizer and reused constitutive-oracle receipts."""
import json,subprocess,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-glottis'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(args,name,inp=None):
 p=subprocess.run(args,cwd=ROOT,input=inp,text=True,capture_output=True);(OUT/(name+'.stdout')).write_text(p.stdout);(OUT/(name+'.stderr')).write_text(p.stderr);assert p.returncode==0,(args,p.stdout,p.stderr);return p.stdout

def main():
 pattern='^ducklips_(continuous_glottis|contour_pressure_(storage|flow|loss)|glottal_storage_reference)$'
 for folder,name in [('build','ctest'),('.build/asan','asan-ctest')]:command(['ctest','--test-dir',folder,'-R',pattern,'--output-on-failure'],name)
 for mode in ['free','driven','lower','upper']:
  release=command(['build/ducklips_continuous_glottis_probe',mode,'10','.000001'],mode+'-release-smoke');sanitized=command(['.build/asan/ducklips_continuous_glottis_probe',mode,'10','.000001'],mode+'-asan-smoke');assert release==sanitized
 a=command(['build/ducklips_continuous_glottis_tangent_probe'],'tangent-release-smoke');b=command(['.build/asan/ducklips_continuous_glottis_tangent_probe'],'tangent-asan-smoke');assert a==b
 # Reuse the independently checked constitutive fixture from 6426640,
 # with its saved input/output hashes verified before comparison.
 previous=json.loads((ROOT/'reports/contour_pressure_loss_validation.json').read_text());inp_path=ROOT/'out/slit-kinetic/loss-input.txt';expected_path=ROOT/'out/slit-kinetic/loss-output.txt'
 for p in [inp_path,expected_path]:assert sha(p)==previous['hashes'][str(p.relative_to(ROOT))]
 actual=command(['build/ducklips_contour_pressure_loss_probe','--batch'],'loss-regression',inp_path.read_text());old=[json.loads(s)['resistance'] for s in expected_path.read_text().splitlines()];new=[json.loads(s)['resistance'] for s in actual.splitlines()];errors=[]
 for a,b in zip(old,new):
  scale=max(abs(np.array(a)));error=float(max(abs(np.array(a)-b))/scale);assert error<1e-10;errors.append(error)
 rows=json.loads((ROOT/'reports/contour_pressure_reference.json').read_text())['cases_detail'];recover=(ROOT/'out/contour-pressure/recover-input.txt').read_text();states=[json.loads(s) for s in command(['build/ducklips_contour_pressure_probe','--recover'],'inverse-regression',recover).splitlines()];error=max(max(abs(np.array(s['effort'])-row['input'][2:])) for s,row in zip(states,rows));assert error<1e-4
 sweeps=[json.loads(command(['build/ducklips_contour_pressure_probe','--sweep',str(n)],f'seam-{n}')) for n in [100,200,400,800]];ratios=[abs(a['work_residual_J']/b['work_residual_J']) for a,b in zip(sweeps,sweeps[1:])];assert all(3.5<v<4.5 for v in ratios[-2:]) and abs(sweeps[-1]['work_residual_J'])<1e-14
 paths=[ROOT/'CMakeLists.txt',Path(__file__),*ROOT.glob('tools/continuous_glottis*'),ROOT/'src/continuous_glottis.cpp',ROOT/'include/ducklips/continuous_glottis.hpp',ROOT/'src/gas_storage.cpp',ROOT/'include/ducklips/gas_storage.hpp',ROOT/'src/contour_pressure_storage.cpp',*ROOT.glob('reports/continuous_glottis_*reference.json'),ROOT/'reports/continuous_glottis_storage_audit.json',*OUT.glob('*.stdout'),*OUT.glob('*.stderr')]
 for folder in ['build','.build/asan']:
  paths.extend(ROOT/folder/name for name in ['ducklips_continuous_glottis_probe','ducklips_continuous_glottis_tangent_probe','ducklips_glottal_storage_probe','ducklips_contour_pressure_probe','ducklips_contour_pressure_flow_probe','ducklips_contour_pressure_loss_probe'])
 report=dict(base_commit='6426640d3ffb116cd90549533c4f0d1c177fd3a1',scope='coupled motion within fixed topology; not a production speech renderer',focused_ctests_passed=5,sanitizer_ctests_passed=5,sanitizers='address,undefined',sanitized_four_motion_controls_and_tangent_byte_identical=True,constitutive_fixture_reused_from='6426640',constitutive_maximum_scaled_change=max(errors),inverse_maximum_effort_error_Pa=float(error),seam_sweeps=sweeps,seam_refinement_ratios=ratios,hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))})
 (ROOT/'reports/continuous_glottis_validation.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS five focused tests/sanitizers, identical motion/tangent sanitizer outputs, constitutive/inverse/seam regressions')
if __name__=='__main__':main()
