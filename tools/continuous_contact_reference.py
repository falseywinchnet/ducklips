#!/usr/bin/env python3
"""Check native contact refinement and collect exact validation provenance."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-contact'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 experiment=json.loads((OUT/'experiment.json').read_text());assert len(experiment)==6 and all(r['exit']==0 for r in experiment);records=[];paths=[]
 for mode in ['lower','upper']:
  runs=[];vectors=[]
  for n in [2000,4000,8000]:
   trace=OUT/f'cycle-{mode}-{n}.jsonl';summary=OUT/f'cycle-{mode}-{n}.stderr';paths.extend([trace,summary]);rows=[json.loads(s) for s in trace.read_text().splitlines()];s=json.loads(summary.read_text());assert s['closures']==s['reopenings']==1 and s['minimum_area']<0;assert rows[0]['time']==0 and rows[-1]['time']==.002 and rows[-1]['step']==n
   mass=max(abs(r['mass_balance']) for r in rows);work=max(abs(r['balance']) for r in rows);assert mass<1e-18 and max(work,s['maximum_balance_J'])<1e-12;assert s['numerical_loss_J']>0 and s['physical_loss_J']>0
   for r in rows:
    for area,flow in [(r['area'][0],r['flow'][0]),(r['area'][1],r['flow'][2])]:
     if area<=0:assert flow==0
   final=rows[-1];vectors.append(np.r_[np.array(final['area'])/1e-6,np.array(final['effort'])/1000,np.array(final['flow'])/1e-5,np.array(final['area_rate'])/.001]);runs.append(dict(steps=n,maximum_sampled_mass_residual_m3=mass,final=final,**s))
  differences=[float(np.max(np.abs(a-b))) for a,b in zip(vectors,vectors[1:])];assert 0<differences[1]<differences[0];assert all(a['numerical_loss_J']>b['numerical_loss_J'] for a,b in zip(runs,runs[1:]));events={}
  for key in ['first_closed_s','first_reopened_s']:
   spread=abs(runs[-1][key]-runs[-2][key]);assert spread<=.002/4000*1.01;events[key+'_fine_grid_change']=spread
  records.append(dict(mode=mode,runs=runs,successive_scaled_endpoint_differences=differences,endpoint_refinement_ratio=differences[0]/differences[1],**events))
 for name in ['final-ctest.log','verified-asan-ctest.log']:
  assert '100% tests passed, 0 tests failed out of 7' in (OUT/name).read_text();paths.append(OUT/name)
 for name in ['final-profile-reference.log','final-storage-reference.log','final-regression.log']:
  assert (OUT/name).read_text().startswith('PASS');paths.append(OUT/name)
 paths.extend([OUT/'experiment.json',OUT/'regression.json',ROOT/'CMakeLists.txt',*ROOT.glob('tools/continuous_contact*'),*ROOT.glob('tools/contour_profile_limit*'),ROOT/'tools/linear_contour_storage_reference.py',ROOT/'reports/contour_profile_limits.json',ROOT/'reports/linear_contour_storage_reference.json',*ROOT.glob('src/contour_pressure*'),ROOT/'src/continuous_glottis.cpp',*ROOT.glob('include/ducklips/contour_pressure*'),ROOT/'include/ducklips/continuous_glottis.hpp',OUT/'final-build.log',OUT/'final-asan-build.log'])
 for folder in ['build','.build/asan']:
  paths.extend(ROOT/folder/name for name in ['ducklips_continuous_contact_controls','ducklips_continuous_glottis_probe','ducklips_contour_pressure_probe','ducklips_contour_pressure_flow_probe','ducklips_contour_pressure_loss_probe','ducklips_glottal_storage_probe'])
 paths.extend([ROOT/'build/ducklips_continuous_contact_probe',ROOT/'build/ducklips_contour_profile_limit_probe'])
 report=dict(base_commit='dbc880fb3334d61e896217507c7c1f56e883be6f',scope='two free full-span tip-contact cycles and two exact-zero pressure-driven onset controls; not sustained phonation or speech',endpoint_scales=dict(area_m2=1e-6,effort_Pa=1000,flow_m3_s=1e-5,area_rate_m2_s=.001),cycles=records,release_ctests_passed=7,address_undefined_sanitizer_ctests_passed=7,default_regressions=json.loads((OUT/'regression.json').read_text()),hashes={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))});(ROOT/'reports/continuous_contact_reference.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS two contact refinements, decreasing numerical loss, mass/work ledger, and seven release/sanitizer controls')
if __name__=='__main__':main()
