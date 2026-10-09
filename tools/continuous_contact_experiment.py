#!/usr/bin/env python3
"""Run both linear-span sealing/reopening controls on three time grids."""
import concurrent.futures,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-contact'
def run(case):
 mode,n=case;command=[str(ROOT/'build/ducklips_continuous_contact_probe'),mode,str(n),'.002','linear','dissipative'];start=time.monotonic();p=subprocess.run(command,capture_output=True,text=True);(OUT/f'cycle-{mode}-{n}.jsonl').write_text(p.stdout);(OUT/f'cycle-{mode}-{n}.stderr').write_text(p.stderr);return dict(command=command,mode=mode,steps=n,seconds=time.monotonic()-start,exit=p.returncode,error=p.stderr)
def main():
 OUT.mkdir(parents=True,exist_ok=True);records=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for r in pool.map(run,[(m,n) for m in ['lower','upper'] for n in [2000,4000,8000]]):records.append(r);print(json.dumps(r),flush=True)
 (OUT/'experiment.json').write_text(json.dumps(records,indent=2)+'\n');assert all(r['exit']==0 for r in records)
if __name__=='__main__':main()
