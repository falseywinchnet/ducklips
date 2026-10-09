#!/usr/bin/env python3
"""Small perturbations of pressure equilibria, compared with independent modes."""
import concurrent.futures,json,subprocess,time
from pathlib import Path
import numpy as np
from continuous_phonation_reference import Reference,YSCALE,REST
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-phonation'
def run(case):
 pressure,n,y,direction=case;duration=.002;epsilon=1e-4;rows=[' '.join(map(str,[n,duration,pressure,REST,*((y+sign*epsilon*direction)*YSCALE)])) for sign in [-1,1]];inp='\n'.join(rows)+'\n';name=f'mode-{pressure}-{n}';(OUT/(name+'.input')).write_text(inp);start=time.monotonic();p=subprocess.run([str(ROOT/'build/ducklips_continuous_phonation_step_probe')],input=inp,text=True,capture_output=True);(OUT/(name+'.jsonl')).write_text(p.stdout);(OUT/(name+'.stderr')).write_text(p.stderr);r=dict(pressure=pressure,steps=n,duration=duration,epsilon=epsilon,exit=p.returncode,seconds=time.monotonic()-start);print(json.dumps(r),flush=True);return r

def main():
 ref=Reference();cases=[];modes=[]
 for pressure in [700,1400]:
  x,y,error=ref.equilibrium(pressure,[.12,-.003,3.7,4.2,1.2]);A=ref.jacobian(y,pressure,.0005);eigen,vectors=np.linalg.eig(A);j=max((j for j in range(9) if eigen[j].imag>0),key=lambda j:eigen[j].real);vector=vectors[:,j];vector*=np.exp(-1j*np.angle(vector[2]));direction=vector.real;direction/=np.max(np.abs(direction));modes.append(dict(pressure=pressure,equilibrium=y.tolist(),direction=direction.tolist(),jacobian=A.tolist(),eigenvalue=[float(eigen[j].real),float(eigen[j].imag)]));cases.extend((pressure,n,y,direction) for n in [200,400,800])
 (OUT/'modal-inputs.json').write_text(json.dumps(modes,indent=2)+'\n')
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(run,cases))
 assert all(r['exit']==0 for r in records);assert all(json.loads(line)['accepted'] for r in records for line in (OUT/f"mode-{r['pressure']}-{r['steps']}.jsonl").read_text().splitlines());(OUT/'modal-experiment.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
