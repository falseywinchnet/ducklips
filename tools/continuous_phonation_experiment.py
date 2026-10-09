#!/usr/bin/env python3
"""Reproduce four driven-contact grids from the recorded physical seed."""
import concurrent.futures,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-phonation'
def run(n):
 command=[str(ROOT/'build/ducklips_continuous_phonation_probe'),str(n),'.024','1400','3e-6','--state','adaptive'];start=time.monotonic();p=subprocess.run(command,input=(OUT/'contact-seed.input').read_text(),text=True,capture_output=True);(OUT/f'adaptive-{n}.jsonl').write_text(p.stdout);(OUT/f'adaptive-{n}.stderr').write_text(p.stderr);r=dict(command=command,exit=p.returncode,seconds=time.monotonic()-start);print(json.dumps(r),flush=True);return r

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if not (OUT/'contact-seed.input').exists():
  seed=json.loads((ROOT/'reports/continuous_phonation_motion.json').read_text())['seeds']['contact_state'];(OUT/'contact-seed.input').write_text(' '.join(map(str,seed))+'\n')
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(run,[2400,4800,9600,19200]))
 (OUT/'contact-experiment.json').write_text(json.dumps(records,indent=2)+'\n');assert all(r['exit']==0 for r in records)
if __name__=='__main__':main()
