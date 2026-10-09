#!/usr/bin/env python3
"""Generate all coupled-motion cases consumed by continuous_glottis_reference."""
import concurrent.futures,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-glottis'
def run(case):
 mode,n,duration,amplitude=case;name=f'{mode}-{n}'+(f'-{amplitude}' if mode=='linear' else '')+('-long' if duration>.0001 else '')
 command=[str(ROOT/'build/ducklips_continuous_glottis_probe'),mode,str(n),str(duration),str(amplitude)];start=time.monotonic();p=subprocess.run(command,capture_output=True,text=True);(OUT/(name+'.json')).write_text(p.stdout);(OUT/(name+'.stderr')).write_text(p.stderr)
 return dict(command=command,returncode=p.returncode,seconds=time.monotonic()-start,name=name,error=p.stderr)
def main():
 OUT.mkdir(parents=True,exist_ok=True);cases=[(m,n,.0001,1.) for m in ['free','driven','lower','upper'] for n in [100,200,400,800,1600,3200]]
 cases += [('lower',6400,.0001,1.)]
 cases += [('linear',n,.00001,a) for n in [100,200,400,800,1600] for a in [.001,.0001]]+[('driven',1000,.001,1.)];records=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for record in pool.map(run,cases):records.append(record);print(json.dumps(record),flush=True)
 (OUT/'experiment.json').write_text(json.dumps(records,indent=2)+'\n');assert all(r['returncode']==0 for r in records)
if __name__=='__main__':main()
