#!/usr/bin/env python3
"""Verify joint-update evidence and locate the remaining coupled-source error."""
import hashlib,json
from pathlib import Path
import numpy as np
from continuous_viscous_reference import YSCALE,XSCALE,ZSCALE,VZSCALE,C,DURATION
ROOT=Path(__file__).resolve().parents[1]
def sha(path):
 with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def main():
 verified={}
 def verify(value):
  if isinstance(value,dict):
   for key,item in value.items():
    if key in ['hashes','files','input_poses'] and isinstance(item,dict):
     for name,digest in item.items():
      assert sha(ROOT/name)==digest,('stale artifact',name);verified[name]=digest
    else:verify(item)
  elif isinstance(value,list):
   for item in value:verify(item)
 reports={}
 for name in ['relaxation_graph_reference','continuous_joint_reference','planar_joint_reference','planar_joint_fine_reference','planar_joint_asan_reference']:
  path=ROOT/f'reports/{name}.json';reports[name]=json.loads(path.read_text());verify(reports[name]);verified[str(path.relative_to(ROOT))]=sha(path)
 reg_path=ROOT/'out/joint-regressions/receipt.json';reg=json.loads(reg_path.read_text());verify(reg);assert len(reg['cases'])==31 and all(c['exact'] for c in reg['cases']);verified[str(reg_path.relative_to(ROOT))]=sha(reg_path)
 tests=[]
 for folder,log,follow in [('build','out/joint-ctest.log','out/joint-fixture-ctest.log'),('.build/asan','out/joint-asan-ctest.log','out/joint-fixture-asan-ctest.log')]:
  assert '100% tests passed, 0 tests failed out of 20' in (ROOT/log).read_text()
  assert '100% tests passed, 0 tests failed out of 1' in (ROOT/follow).read_text()
  for name in [log,follow]:verified[name]=sha(ROOT/name)
  tests.append(dict(build=folder,passed=20,failed=0,source_control_after_fixture_protocol_change=1))
 assert len(reports['planar_joint_asan_reference']['runs'])==4
 assert len(reports['planar_joint_reference']['runs'])==12 and len(reports['planar_joint_fine_reference']['runs'])==4
 assert len(reports['planar_joint_fine_reference']['comparisons'])==2
 generic=reports['relaxation_graph_reference'];assert generic['random_cases']==80
 assert max(abs(r['average_flow_ratio']-1) for r in generic['forced_flow'])<1e-11
 oracle=np.load(ROOT/'out/continuous-viscous/continuous.npz');full=oracle['fine'].T;location=[]
 for n in [200,400,800,1600]:
  data=np.loadtxt(ROOT/f'out/continuous-joint/fixed-{n}.txt');reference=full[np.arange(1,n+1)*(3200//n)]
  native=np.c_[data[:,:9]/YSCALE,data[:,9:16]/XSCALE,data[:,44]/.1,data[:,49:121]/ZSCALE.ravel(),data[:,126:198]/VZSCALE.ravel()]
  error=abs(native-np.c_[reference[:,:89],reference[:,91:163]]);peak=np.unravel_index(error.argmax(),error.shape);endpoint=int(error[-1].argmax())
  assert peak[1]<2 and endpoint<2
  groups={name:float(error[:,a:b].max()) for name,a,b in [('source',0,9),('tract',9,16),('radiation',16,17),('thermal',17,89),('viscous',89,161)]}
  location.append(dict(steps=n,maximum_sampled_scaled_error_by_group=groups,peak_source_effort_index=int(peak[1]),peak_time_s=float((peak[0]+1)*DURATION/n),peak_source_effort_error_Pa=float(error[peak]*1000),endpoint_source_effort_index=endpoint,endpoint_source_effort_error_Pa=float(error[-1,endpoint]*1000),maximum_tract_pressure_error_Pa=float(np.max(error[:,9:12]*XSCALE[:3]/C))))
 closure=[]
 for case in ['closing','opening']:
  closed_steps=retired=0;surviving=0.;numeric=physical=0.;energy_residual=air_residual=0.
  for interval in [0,1]:
   data=json.loads((ROOT/f'out/planar-joint-fine/{case}-{interval}-2048.json').read_text());previous=data['initial']
   for step in data['steps']:
    sites=[s for s in step['viscous_sites'] if s['branch']==6];prior=[s for s in previous['viscous_sites'] if s['branch']==6];assert len(sites)==9
    if sites[0]['segment_inertance']==0:
     closed_steps+=1;assert sites[0]['half_gap']==0 and sites[0]['branch_conductance']==0 and all(z==0 for z in sites[0]['momentum']);assert step['momentum'][6]==0
     if prior[0]['segment_inertance']>0:retired+=1
     if previous['conductance'][6]==0:assert step['average_flow'][6]==0
     surviving=max(surviving,sum(s['energy'] for s in sites[1:]))
    r=step['joint_passive'];numeric+=r['numerical_dissipation'];physical+=sum(r[k] for k in ['thermal_dissipation','viscous_dissipation','radiated_energy']);energy_residual=max(energy_residual,abs(r['energy_residual']));air_residual=max(air_residual,abs(r['air_residual']));previous=step
  assert (case=='closing' and retired==1 and closed_steps>0 and surviving>0) or (case=='opening' and retired==0)
  closure.append(dict(case=case,closed_endpoint_steps=closed_steps,core_retirements=retired,maximum_surviving_blocked_profile_energy_J=surviving,physical_dissipation_J=physical,passive_numerical_dissipation_J=numeric,maximum_passive_work_residual_J=energy_residual,maximum_passive_air_residual_m3=air_residual))
 paths=[Path(__file__),ROOT/'README.md',ROOT/'RESEARCH.md',ROOT/'reports/CONTINUOUS_JOINT_RELAXATION.md',ROOT/'CMakeLists.txt',ROOT/'out/nato-spatial-pressure.wav',ROOT/'tools/relaxation_graph_probe.cpp',ROOT/'.build/asan/ducklips_graph_continuous_glottis_probe']
 for path in paths:verified[str(path.relative_to(ROOT))]=sha(path)
 assert verified['out/nato-spatial-pressure.wav']=='8284f366f2df5f3496d3139f7cbbaf4a8947ad839737d6a7e40b78635083b307'
 result=dict(scope='Common passive update, continuous-source comparison, actual moving geometry and exact closure. No new speech or speed claim. Saved independent source ODE reused with generator and artifact hashes verified.',ctest=tests,historical_exact_json_regressions=31,moving_sanitizer_intervals=4,source_error_location=location,finest_moving_closure=closure,verified_artifact_count=len(verified),hashes=dict(sorted(verified.items())))
 (ROOT/'reports/continuous_joint_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(verified_artifact_count=len(verified),source_error_location=location,finest_moving_closure=closure),indent=2))
if __name__=='__main__':main()
