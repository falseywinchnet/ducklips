#!/usr/bin/env python3
"""Fifth-order passive update in the actual moving continuous-source controls."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from planar_continuous_reference import ROOT,setup,vector,run,sha
from planar_joint_reference import audit_joint

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--subdivisions',nargs='+',type=int,default=[64,256,512,1024]);parser.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_planar_inertive_motion_probe');parser.add_argument('--output',type=Path,default=ROOT/'out/planar-radau5');parser.add_argument('--report',type=Path,default=ROOT/'reports/planar_radau5_reference.json');options=parser.parse_args()
 out=options.output.resolve();out.mkdir(parents=True,exist_ok=True);binary=options.binary.resolve()
 rows=list(csv.DictReader((ROOT/'out/lip-profile/nato.csv').open()));airway=list(csv.DictReader((ROOT/'out/nato-spatial-airway.csv').open()));runs=[];comparisons=[];inputs={}
 for name,start,pressure in [('closing',1284,0.),('opening',1374,700.)]:
  last=None
  for n in options.subdivisions:
   previous=None;mass0=None;net=0.
   for interval in [0,1]:
    index=start+interval;p,q,paths=setup(index/1000);profiles=np.array([vector(rows[i]) for i in [index,index+1]])
    epi=np.array([[[float(airway[t][f'{field}_{i}']) for i in reversed(range(8))] for field in ['area','length']] for t in [index,index+1]])
    data,receipt=run(f'{name}-{interval}',p,q,profiles,epi,n,pressure,previous,binary,out,True,.03,.03,.03,8,True,True,True)
    receipt.update(audit_joint(data,p,q,profiles,epi,n==options.subdivisions[0],True,.03,.03,.03,8,True,3));runs.append(receipt)
    if mass0 is None:mass0=receipt['initial_total_mass_m3']
    net+=receipt['net_external_transfer_m3'];previous=data['steps'][-1];inputs.update({str(path.relative_to(ROOT)):sha(path) for path in paths})
   residual=abs(receipt['final_total_mass_m3']-mass0-net);assert residual<1e-17
   pressure_vector=np.array(previous['pressure']);record=dict(case=name,subdivisions=n,endpoint_pressure_difference_Pa=None if last is None else float(max(abs(pressure_vector-last))),pressure=pressure_vector.tolist(),whole_control_mass_residual_m3=residual)
   comparisons.append(record);last=pressure_vector;print(json.dumps(record),flush=True)
 paths=[Path(__file__),binary,ROOT/'tools/planar_continuous_reference.py',ROOT/'tools/planar_joint_reference.py',ROOT/'src/graph_joint_relaxation.cpp',ROOT/'src/relaxation_graph.cpp',ROOT/'include/ducklips/relaxation_graph.hpp',ROOT/'include/ducklips/graph_continuous_glottis.hpp',ROOT/'tools/planar_inertive_motion_probe.cpp',ROOT/'bodies/reference.body',ROOT/'out/lip-profile/nato.csv',ROOT/'out/nato-spatial-airway.csv']
 report=dict(scope='Three-stage passive Radau update with midpoint geometry and existing implicit gas/tissue source. Independent full-stage audit at the first resolution, using native coupled source flux; higher resolutions check work/mass and refinement. The entire moving/source method is not fifth order. No new speech.',runs=runs,comparisons=comparisons,input_poses=inputs,hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
 options.report.resolve().write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
