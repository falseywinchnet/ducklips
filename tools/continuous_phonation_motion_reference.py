#!/usr/bin/env python3
"""Collect nonlinear trajectories without treating a failed run as a pass."""
import hashlib,json
from pathlib import Path
import numpy as np
from continuous_phonation_reference import YSCALE
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-phonation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def state(r):return np.r_[r['effort'],r['area'],r['flow'],r['area_rate']]
def trace(name):return [json.loads(l) for l in (OUT/(name+'.jsonl')).read_text().splitlines()]
def summary(name):return json.loads((OUT/(name+'.stderr')).read_text())
def main():
 report=dict(base_commit='dc04d2abbf49ba076774c0602d54144cd2bfc2bc',scope='pressure-driven onset and transient contact cycles; not a settled voice or new speech',seeds={name:list(map(float,(OUT/(file+'.input')).read_text().split())) for name,file in [('mature_state','mature'),('contact_state','contact-seed')]});paths=[]
 coarse=json.loads((OUT/'open-ode.json').read_text());fine=json.loads((OUT/'open-ode-fine.json').read_text());assert coarse['success'] and fine['success'];assert coarse['near_contact_event'] and fine['near_contact_event'];event_change=abs(coarse['final_time']-fine['final_time']);assert event_change<1e-7
 i=min(range(len(fine['times'])),key=lambda i:abs(fine['times'][i]-.04));target=np.array(fine['states'][i]);error=[];symmetric=[]
 for n in [4000,8000]:
  name=f'symmetric-{n}';r=trace(name)[-1];s=summary(name);assert s['complete'];e=float(np.max(np.abs(state(r)-target)/YSCALE));error.append(e);symmetric.append(dict(steps=n,dt_s=.04/n,maximum_scaled_endpoint_error=e,receipt=s));paths.extend([OUT/(name+'.jsonl'),OUT/(name+'.stderr')])
 assert error[1]<3e-4 and 3.5<error[0]/error[1]<4.5
 report['independent_open_reference']=dict(coarse={k:v for k,v in coarse.items() if k not in ['times','states']},fine={k:v for k,v in fine.items() if k not in ['times','states']},near_contact_time_change_s=event_change,symmetric=symmetric,symmetric_endpoint_refinement_ratio=error[0]/error[1],source_stage='after roundoff fix, before coordinate rebase; independently verified equations unchanged by rebase')
 diagnostics=[]
 for name,stage in [('pressure1400-15000','before both fixes'),('pressure1400-30000','before both fixes'),('mature-10000','roundoff fix, all-step dissipative'),('mature-20000','roundoff fix, all-step dissipative'),('transition-8000','roundoff fix, before coordinate rebase'),('transition-16000','roundoff fix, before coordinate rebase'),('rebased-8000','both fixes, no subdivision'),('rebased-16000','both fixes, no subdivision')]:
  r=trace(name);s=summary(name);diagnostics.append(dict(name=name,stage=stage,receipt=s,last_saved_state=r[-1]));paths.extend([OUT/(name+'.jsonl'),OUT/(name+'.stderr')])
 report['developmental_runs']=diagnostics
 runs=[];vectors=[]
 for n in [2400,4800,9600,19200]:
  name=f'adaptive-{n}';rows=trace(name);s=summary(name);assert s['complete'] and s['closures']==s['reopenings']==3 and s['dissipative_transitions']==6;assert s['maximum_work_residual']<1e-11 and s['maximum_mass_residual']<1e-18;assert s['segments']>=n and s['minimum_step_s']>0
  for r in rows:
   for a,q in [(r['area'][0],r['flow'][0]),(r['area'][1],r['flow'][2])]:
    if a<=0:assert q==0
  vectors.append(state(rows[-1])/YSCALE);runs.append(dict(receipt=s,initial_state=rows[0],final_state=rows[-1]));paths.extend([OUT/(name+'.jsonl'),OUT/(name+'.stderr')])
 differences=[float(np.max(np.abs(a-b))) for a,b in zip(vectors,vectors[1:])];assert all(b<a for a,b in zip(differences,differences[1:]));assert all(a['final_state']['numerical_loss']>b['final_state']['numerical_loss'] for a,b in zip(runs,runs[1:]))
 report['driven_contact_refinement']=dict(duration_s=.024,pressure_Pa=1400,runs=runs,successive_scaled_endpoint_differences=differences,endpoint_refinement_ratios=[a/b for a,b in zip(differences,differences[1:])],scale=YSCALE.tolist(),scope='three transient closures and reopenings, not proof of a settled periodic limit cycle')
 contacts=[]
 for mode in ['lower','upper']:
  name='contact-'+mode;s=summary(name);assert s['closures']==s['reopenings']==1;assert s['maximum_balance_J']<1e-12;current=trace(name)[-1];previous=json.loads((ROOT/f'out/continuous-contact/cycle-{mode}-2000.jsonl').read_text().splitlines()[-1]);difference=float(np.max(np.abs(state(current)-state(previous))/YSCALE));assert difference<1e-6;contacts.append(dict(mode=mode,receipt=s,maximum_scaled_endpoint_change=difference));paths.extend([OUT/(name+'.jsonl'),OUT/(name+'.stderr')])
 report['free_contact_regressions']=contacts;report['default_regressions']=json.loads((OUT/'regression.json').read_text());report['release_ctests_passed']=10;report['address_undefined_sanitizer_ctests_passed']=10
 for name in ['final-ctest.log','final-asan-ctest.log']:
  assert '100% tests passed, 0 tests failed out of 10' in (OUT/name).read_text();paths.append(OUT/name)
 assert (OUT/'rebase-control.log').read_text().startswith('PASS');report['restart_control']=(OUT/'rebase-control.log').read_text().strip()
 paths.extend([Path(__file__),*ROOT.glob('tools/continuous_phonation*'),ROOT/'tools/continuous_glottis_rebase_controls.cpp',ROOT/'CMakeLists.txt',ROOT/'src/continuous_glottis.cpp',ROOT/'include/ducklips/continuous_glottis.hpp',OUT/'open-ode.json',OUT/'open-ode-fine.json',OUT/'contact-seed.input',OUT/'mature.input',OUT/'regression.json',OUT/'rebase-control.log',OUT/'probe-before-subdivision',ROOT/'build/ducklips_continuous_phonation_probe',ROOT/'build/ducklips_continuous_phonation_controls',ROOT/'build/ducklips_continuous_glottis_rebase_controls',ROOT/'.build/asan/ducklips_continuous_phonation_controls',ROOT/'.build/asan/ducklips_continuous_glottis_rebase_controls'])
 report['hashes']={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))};(ROOT/'reports/continuous_phonation_motion.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS native pressure-driven contact refinement, independent open trajectory, restart invariance and regressions')
if __name__=='__main__':main()
