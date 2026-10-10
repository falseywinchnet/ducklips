#!/usr/bin/env python3
"""Preserve a tissue-source experiment, depict its whole sound and deformation."""
from pathlib import Path
import hashlib,json,subprocess,sys,shutil
import numpy as np
from scipy.io import wavfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
root=Path(sys.argv[1]).resolve()
fs,audio=wavfile.read(root/'source-raw.wav')
assert fs==48000 and np.isfinite(audio).all()
reference_gain=315.98259903736243
gain=min(reference_gain,.7/max(float(abs(audio).max()),1e-30))
wavfile.write(root/'source.wav',fs,np.rint(audio*gain*32767).astype(np.int16))
record={'status':'unaccepted source experiment','sample_rate':fs,'frames':len(audio),
 'gain':gain,'processing':'one constant gain, no EQ, noise, cuts, or fades',
 'raw_sha256':hashlib.sha256((root/'source-raw.wav').read_bytes()).hexdigest(),
 'playback_sha256':hashlib.sha256((root/'source.wav').read_bytes()).hexdigest(),
 'sources':{},'enhanced_inspected':False}
if (root/'setup.json').exists():record['setup']=json.loads((root/'setup.json').read_text())
if (root/'timing.json').exists():record['timing']=json.loads((root/'timing.json').read_text())
snapshot=root/'snapshot';snapshot.mkdir(exist_ok=True)
executable=root/'render'
if not executable.exists():
    candidate=ROOT/'build/ducklips_tissue_source_study'
    if candidate.exists():shutil.copy2(candidate,executable)
if executable.exists():record['executable_sha256']=hashlib.sha256(executable.read_bytes()).hexdigest()
for path in ['include/ducklips/tissue_source.hpp','src/tissue_source.cpp','tools/render_tissue_source.cpp','tools/package_tissue_round.py',
             'include/ducklips/tissue_body_data.hpp','data/tissue_body.json','tools/compile_tissue_body.py']:
    data=(ROOT/path).read_bytes();(snapshot/Path(path).name).write_bytes(data)
    record['sources'][path]=hashlib.sha256(data).hexdigest()
(root/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
a=np.genfromtxt(root/'trace.csv',delimiter=',',names=True)
fig,ax=plt.subplots(4,1,figsize=(13,10),sharex=True)
ax[0].plot(a['time'],a['lung']);ax[0].set_ylabel('breath pressure, Pa')
for name in ['lower','center','upper']:ax[1].plot(a['time'],a[name]*1e3,label=name)
ax[1].legend();ax[1].set_ylabel('surface motion, mm')
ax[2].plot(a['time'],a['flow']*1e6);ax[2].set_ylabel('flow, cm³/s')
ax[3].plot(a['time'],audio[:,0]*128);ax[3].set_ylabel('source audification, Pa');ax[3].set_xlabel('s')
fig.suptitle(root.name+': complete pressure-driven tissue experiment');fig.tight_layout();fig.savefig(root/'trace-full.png',dpi=110)
for axis in ax:axis.set_xlim(.5,.54)
fig.suptitle(root.name+': supplementary cycle detail');fig.savefig(root/'trace-detail.png',dpi=110)
s=np.genfromtxt(root/'surface.csv',delimiter=',',names=True)
times=np.unique(s['time']);target=np.linspace(.5,.516,8)
fig,ax=plt.subplots(2,4,figsize=(14,7))
for axis,t in zip(ax.ravel(),target):
    actual=times[np.argmin(abs(times-t))];rows=s[s['time']==actual]
    shape=(int(rows['z'].max()+1),int(rows['x'].max()+1))
    image=axis.imshow(rows['gap'].reshape(shape)*1e3,origin='lower',extent=[0,15,0,4],aspect='auto',cmap='coolwarm',vmin=-.2,vmax=1.5,interpolation='nearest')
    axis.set_title(f'{actual:.4f} s');axis.set_xlabel('length, mm');axis.set_ylabel('height, mm')
fig.suptitle('Full medial surface at successive times; signed gap in mm')
fig.tight_layout();fig.savefig(root/'surface-motion.png',dpi=110)
print('AUDIO AND PHYSICAL DEPICTIONS READY',root,flush=True)
subprocess.run([sys.executable,str(ROOT/'tools/enhanced_specgram/render.py'),str(root/'source.wav'),'--channel','0','--max-hz','8000','--hop','128','--display-reference',str(.06376754775627493*gain/reference_gain),'--out',str(root/'enhanced')],check=True)
receipt=json.loads((root/'enhanced.json').read_text())
assert receipt['source_sha256']==record['playback_sha256'] and receipt['start_sample']==0 and receipt['end_sample']==len(audio)
record['enhanced_complete']=True;(root/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
