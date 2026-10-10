#!/usr/bin/env python3
"""Continuous-source comparison with the three-stage passive solve."""
import json
from pathlib import Path
import numpy as np
from continuous_joint_reference import ROOT,sha,SavedOracle,run,DURATION
OUT=ROOT/'out/continuous-radau5'
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 receipt_path=ROOT/'reports/continuous_viscous_reference.json';saved_receipt=json.loads(receipt_path.read_text())
 oracle_path=ROOT/'out/continuous-viscous/continuous.npz'
 reused=['tools/continuous_viscous_reference.py','tools/continuous_outlet_column_reference.py','tools/continuous_phonation_reference.py','tools/graph_continuous_glottis_reference.py','tools/radiation_boundary_reference.py','tools/slit_thermal_reference.py','tools/slit_viscous_reference.py','out/continuous-viscous/continuous.npz']
 for name in reused:assert sha(ROOT/name)==saved_receipt['hashes'][name],name
 oracle=SavedOracle(oracle_path);delta=float(np.max(abs(oracle.y-oracle.coarse)));assert delta==saved_receipt['reference_refinement_scaled_error']
 runs=[]
 for n in [200,400,800,1600]:
  record=run(n,oracle,fixture="--joint-fifth-fixture",out=OUT);runs.append(record);print(json.dumps(record),flush=True)
 ratios=[a['endpoint_scaled_error']/b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])]
 assert all(a['endpoint_scaled_error']>b['endpoint_scaled_error'] for a,b in zip(runs,runs[1:])),ratios
 assert all(a['maximum_sampled_scaled_error']>b['maximum_sampled_scaled_error'] for a,b in zip(runs,runs[1:]))
 paths=[Path(__file__),ROOT/'tools/continuous_joint_reference.py',receipt_path,*[ROOT/name for name in reused],ROOT/'src/relaxation_graph.cpp',ROOT/'include/ducklips/relaxation_graph.hpp',ROOT/'src/graph_joint_relaxation.cpp',ROOT/'src/graph_continuous_glottis.cpp',ROOT/'include/ducklips/graph_continuous_glottis.hpp',ROOT/'tools/graph_continuous_glottis_probe.cpp',ROOT/'build/ducklips_graph_continuous_glottis_probe',*sorted(OUT.iterdir())]
 report=dict(scope='Fresh native three-stage fifth-order passive Radau update coupled to the continuous gas/tissue source, fixed cyclic tract, aperture radiation, thermal storage and three viscous segments. Compared to a reused, hash-verified independent continuous ODE trajectory. No new ODE integration was run; all independent generator files and the saved trajectory are verified against their original receipt. Tests trajectory accuracy, not speech.',reference_refinement_scaled_error=delta,duration_s=DURATION,runs=runs,endpoint_refinement_ratios=ratios,hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})
 (ROOT/'reports/continuous_radau5_reference.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
