#!/usr/bin/env python3
"""Collect current fifth-order evidence and compare complete sampled controls."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
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
 for name in ['radau5_reference','lip_relaxation_mode','radau5_compatibility','continuous_radau5_reference','planar_radau5_reference','planar_radau5_asan_reference']:
  p=ROOT/f'reports/{name}.json';r=json.loads(p.read_text());verify(r);reports[name]=r;verified[str(p.relative_to(ROOT))]=sha(p)
 p=ROOT/'out/radau5-regressions/receipt.json';r=json.loads(p.read_text());verify(r);assert len(r['cases'])==31;verified[str(p.relative_to(ROOT))]=sha(p)
 for log in ['out/radau5-ctest.log','out/radau5-asan-ctest.log']:
  assert '100% tests passed, 0 tests failed out of 4' in (ROOT/log).read_text();verified[log]=sha(ROOT/log)
 assert len(reports['radau5_compatibility']['runs'])==18 and len(reports['planar_radau5_asan_reference']['runs'])==4
 assert reports['radau5_reference']['random_cases']==80
 # Historical data hashes are verified, without treating the old source code
 # digests as current. These are saved third-order trajectories, not reruns.
 baseline=json.loads((ROOT/'reports/planar_joint_reference.json').read_text());baseline_files={}
 for run in baseline['runs']:baseline_files.update(run['files'])
 comparisons=[]
 for order,folder in [(3,'planar-joint'),(5,'planar-radau5')]:
  for case in ['closing','opening']:
   for low,high in [(256,512),(512,1024)]:
    maxima=np.zeros(14);samples=0;endpoint=None;peak={};peak_value=-1.
    for interval in [0,1]:
     paths=[ROOT/f'out/{folder}/{case}-{interval}-{n}.json' for n in [low,high]]
     for path in paths:
      name=str(path.relative_to(ROOT));digest=sha(path)
      if order==3:assert digest==baseline_files[name],name
      verified[name]=digest
     coarse,fine=[json.loads(path.read_text()) for path in paths];hi={round(s['geometry']['fraction'],13):s for s in [fine['initial'],*fine['steps']]}
     for a in [coarse['initial'],*coarse['steps']]:
      key=round(a['geometry']['fraction'],13);assert key in hi,(case,interval,key);b=hi[key];assert abs(a['geometry']['fraction']-b['geometry']['fraction'])<1e-13
      delta=abs(np.array(b['pressure'])-a['pressure']);maxima=np.maximum(maxima,delta);samples+=1;endpoint=delta
      if max(delta)>peak_value:
       node=int(np.argmax(delta));peak_value=float(delta[node]);peak=dict(node=node,control_time_s=float((interval+a['geometry']['fraction'])*.001),coarse_pressure_Pa=a['pressure'][node],fine_pressure_Pa=b['pressure'][node])
    comparisons.append(dict(order=order,case=case,coarse_subdivisions=low,fine_subdivisions=high,common_samples=samples,peak_sample=peak,maximum_sampled_pressure_difference_Pa=float(max(maxima)),maximum_mouth_difference_Pa=float(max(maxima[:2])),maximum_lip_difference_Pa=float(max(maxima[2:6])),maximum_throat_difference_Pa=float(max(maxima[6:])),endpoint_pressure_difference_Pa=float(max(endpoint)),sampled_pressure_difference_by_node_Pa=maxima.tolist()))
 paths=[Path(__file__),ROOT/'README.md',ROOT/'RESEARCH.md',ROOT/'reports/LIP_TRANSIENT_ACCURACY.md',ROOT/'CMakeLists.txt',ROOT/'out/nato-spatial-pressure.wav',ROOT/'.build/asan/ducklips_relaxation_graph_probe',ROOT/'.build/asan/ducklips_graph_continuous_glottis_probe',ROOT/'tools/graph_continuous_glottis_probe.cpp',ROOT/'tools/relaxation_graph_probe.cpp']
 for p in paths:verified[str(p.relative_to(ROOT))]=sha(p)
 assert verified['out/nato-spatial-pressure.wav']=='8284f366f2df5f3496d3139f7cbbaf4a8947ad839737d6a7e40b78635083b307'
 result=dict(scope='Fifth-order passive time update in the same physical engine. Current independent references, sanitized controls, saved third-order comparison, and source error remaining. No new speech or speed claim.',release_focused_tests=4,asan_focused_tests=4,sanitized_moving_intervals=4,historical_exact_json_regressions=31,preserved_third_order_fixtures=18,sampled_refinement=comparisons,verified_artifact_count=len(verified),hashes=dict(sorted(verified.items())))
 (ROOT/'reports/radau5_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(verified_artifact_count=len(verified),sampled_refinement=comparisons),indent=2))
if __name__=='__main__':main()
