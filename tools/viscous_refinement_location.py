#!/usr/bin/env python3
"""Locate current endpoint refinement error in the frozen tract's wave modes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from planar_source_reference import B
ROOT=Path(__file__).resolve().parents[1]
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=ROOT/'out/planar-viscous');parser.add_argument('--report',type=Path,default=ROOT/'reports/viscous_refinement_location.json');args=parser.parse_args()
 records=[];paths=[]
 for name in ['closing','opening']:
  pair=[args.input.resolve()/f'{name}-1-{n}.json' for n in [512,1024]];paths+=pair
  a,b=[json.loads(p.read_text())['steps'][-1] for p in pair]
  assert a['volume']==b['volume'] and a['compliance']==b['compliance'] and a['conductance']==b['conductance']
  delta=np.array(b['pressure'])-a['pressure'];i=int(abs(delta).argmax())
  c=np.array(b['compliance']);g=np.array(b['conductance']);sqrtc=np.sqrt(c)
  lam,v=eigh(((B*g)@B.T)/sqrtc[:,None]/sqrtc[None,:])
  charge=np.array(b['charge'])-a['charge'];momentum=np.array(b['momentum'])-a['momentum']
  amplitude=v.T@(charge/sqrtc);speed=v.T@(-B@(g*momentum)/sqrtc)
  energy=.5*amplitude**2;valid=lam>1.;energy[valid]+=.5*speed[valid]**2/lam[valid];j=int(energy.argmax())
  diagnostic=dict(scope='Frozen inviscid endpoint tract wave basis; excludes source, thermal, viscous and radiation memory feedback, and branch circulation not visible to node divergence. Diagnostic basis, not the lossy coupled spectrum.',dominant_error_mode_frequency_Hz=float(np.sqrt(lam[j])/(2*np.pi)),dominant_mode_error_energy_fraction=float(energy[j]/energy.sum()),dominant_mode_last_lip_energy_coordinate_fraction=float(v[5,j]**2),dominant_mode_pressure_error_Pa=float(v[5,j]*amplitude[j]/sqrtc[5]),total_modal_error_energy_J=float(energy.sum()))
  record=dict(case=name,coarse_subdivisions=512,fine_subdivisions=1024,endpoint_pressure_difference_Pa=delta.tolist(),largest_node=i,largest_region='last lip storage section' if i==5 else 'other',maximum_mouth_difference_Pa=float(max(abs(delta[:2]))),maximum_throat_difference_Pa=float(max(abs(delta[6:]))),largest_node_fine_volume_m3=b['volume'][i],largest_node_fine_pressure_Pa=b['pressure'][i],largest_node_fine_compliance_m3_Pa=b['compliance'][i],lip_internal_conductance=b['conductance'][5],lip_outlet_conductance=b['conductance'][6],frozen_inviscid_tract_diagnostic=diagnostic)
  records.append(record);print(json.dumps(dict(case=name,**diagnostic)),flush=True)
 paths += [Path(__file__),ROOT/'tools/planar_inertive_motion_probe.cpp',ROOT/'tools/planar_source_reference.py']
 report=dict(scope='Localize the current endpoint timestep discrepancy. Nodes: mouth 0-1, lips 2-5, throat 6-13. A numerical-error location, not proof of its physical cause or audible consequence.',cases=records,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 args.report.write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
