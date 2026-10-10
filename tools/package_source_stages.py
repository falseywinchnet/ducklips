#!/usr/bin/env python3
"""Package fresh source-isolation runs at one common gain, then depict all unique signals."""
from pathlib import Path
import hashlib,json,sys,subprocess,concurrent.futures
import numpy as np
from scipy.io import wavfile
root=Path(sys.argv[1]).resolve()
raw={p.name.removesuffix('-raw.wav'):wavfile.read(p)[1] for p in root.glob('*-raw.wav')}
gain=.7/max(float(np.max(abs(x))) for x in raw.values())
record={'sample_rate':48000,'duration_each':1.3,'common_gain':gain,'processing':'one common gain; no equalization, cuts, fades or per-clip normalization',
 'source_controls':{'lung_pa':600,'rest_area_m2':3e-6,'tension':1,'tissue_substeps':2,'air_substeps':1,'jitter':0,'shimmer':0,'tremor':0,'aspiration_noise':0,'frication_noise':0,'pitch_servo':'does not control fixed tension'},
 'conditions':{'free':'same inertive inlet discharges into constant-pressure reservoir; all airway branches absent','neutral':'connected AX tract, source controls held constant','ah':'connected AA tract, source controls held constant','closures':'AA plus lip-only closure/release nudges at .35/.50/.80/.95 s; no glottal suppression or word planner'},'files':{},'sequences':{}}
for name,x in raw.items():
 assert np.isfinite(x).all()
 p=root/(name+'.wav');pcm=np.rint(x*gain*32767).astype(np.int16);wavfile.write(p,48000,pcm)
 record['files'][name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'peak':float(np.max(abs(pcm.astype(float)))/32767),'raw_sha256':hashlib.sha256((root/(name+'-raw.wav')).read_bytes()).hexdigest()}
sequences={
 'progression':['01-free-folds','neutral-source','neutral-lips','ah-source','ah-lips','closures-lips'],
 'airway-and-composite':['ah-source','ah-position-1','ah-position-2','ah-position-3','ah-lips','ah-composite-1','ah-composite-2','ah-composite-3']}
for name,names in sequences.items():
 parts=[];timeline=[];at=0
 for i,n in enumerate(names):
  x=wavfile.read(root/(n+'.wav'))[1];parts.append(x);timeline.append({'start':at/48000,'end':(at+len(x))/48000,'signal':n});at+=len(x)
  if i+1<len(names):parts.append(np.zeros((16800,2),dtype=np.int16));at+=16800
 wavfile.write(root/(name+'.wav'),48000,np.concatenate(parts));record['sequences'][name]=timeline
(root/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
print('AUDIO READY',root,'gain',gain,flush=True)
# Exact duplicate signals share a depiction; every unique complete recording
# is analyzed at the preserved full 128-sample lattice and original colors.
seen={};jobs=[]
for name,item in record['files'].items():
 digest=item['sha256']
 if digest not in seen:seen[digest]=name;jobs.append(name)
 item['enhanced_signal']=seen[digest]
(root/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
renderer=Path(__file__).resolve().parent/'enhanced_specgram/render.py'
def depict(name):
 command=[sys.executable,str(renderer),str(root/(name+'.wav')),'--channel','0','--max-hz','8000','--hop','128','--display-reference',str(.17808298734486477*gain/882.4414167632915),'--out',str(root/(name+'-enhanced'))]
 with (root/(name+'-enhanced.log')).open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
 receipt=json.loads((root/(name+'-enhanced.json')).read_text());assert receipt['source_sha256']==record['files'][name]['sha256'];assert receipt['start_sample']==0 and receipt['end_sample']==62400
 print('DEPICTED',name,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(depict,jobs))
record['all_unique_depictions_complete']=True;(root/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
