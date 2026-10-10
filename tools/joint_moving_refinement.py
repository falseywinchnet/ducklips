#!/usr/bin/env python3
"""One finer moving control, with ledger checks and comparison to the audited run."""
import csv,json
import numpy as np
from planar_continuous_reference import ROOT,setup,vector,run,sha
from planar_joint_reference import audit_joint
OUT=ROOT/'out/planar-joint-fine'
def main():
 OUT.mkdir(parents=True,exist_ok=True);rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));airway=list(csv.DictReader((ROOT/'out/nato-spatial-airway.csv').open()))
 runs=[];comparisons=[];inputs={};binary=ROOT/'build/ducklips_planar_inertive_motion_probe'
 for name,start,pressure in [('closing',1284,0.),('opening',1374,700.)]:
  previous=None;mass0=None;net=0.
  for interval in [0,1]:
   index=start+interval;p,q,paths=setup(index/1000);profiles=np.array([vector(rows[i]) for i in [index,index+1]])
   epi=np.array([[[float(airway[t][f'{field}_{i}']) for i in reversed(range(8))] for field in ['area','length']] for t in [index,index+1]])
   data,receipt=run(f'{name}-{interval}',p,q,profiles,epi,2048,pressure,previous,binary,OUT,True,.03,.03,.03,8,True,True)
   receipt.update(audit_joint(data,p,q,profiles,epi,False,True,.03,.03,.03,8,True));runs.append(receipt)
   if mass0 is None:mass0=receipt['initial_total_mass_m3']
   net+=receipt['net_external_transfer_m3'];previous=data['steps'][-1];inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
  residual=abs(receipt['final_total_mass_m3']-mass0-net);assert residual<1e-17
  old=ROOT/f'out/planar-joint/{name}-1-1024.json'
  if old.exists():
   prior=json.loads(old.read_text())['steps'][-1];delta=np.array(previous['pressure'])-prior['pressure'];inputs[str(old.relative_to(ROOT))]=sha(old)
   comparisons.append(dict(case=name,subdivisions=2048,endpoint_pressure_difference_Pa=float(max(abs(delta))),pressure_difference_Pa=delta.tolist(),whole_control_mass_residual_m3=residual))
  print(json.dumps(comparisons[-1] if comparisons and comparisons[-1]['case']==name else dict(case=name,comparison_pending=True)),flush=True)
 paths=[Path(__file__),binary,ROOT/'tools/planar_continuous_reference.py',ROOT/'tools/planar_joint_reference.py',ROOT/'src/graph_joint_relaxation.cpp',ROOT/'src/relaxation_graph.cpp',ROOT/'tools/planar_inertive_motion_probe.cpp']
 report=dict(scope='Fresh 2048-subdivision controls; same moving geometry and work/mass checks. No new independent stage reconstruction at this resolution; see planar_joint_reference.json for the independent 256-subdivision audit. Endpoint comparisons against 1024-subdivision native output.',runs=runs,comparisons=comparisons,input_poses=inputs,hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
 (ROOT/'reports/planar_joint_fine_reference.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':
 from pathlib import Path
 main()
