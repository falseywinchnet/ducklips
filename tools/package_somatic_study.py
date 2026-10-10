#!/usr/bin/env python3
"""Freeze and depict the full freshly rendered controlled-engine recording."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys
import numpy as np
from scipy.io import wavfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];root=Path(sys.argv[1]).resolve()
fs,raw=wavfile.read(root/'voice-raw.wav');assert fs==48000 and np.isfinite(raw).all()
gain=min(315.98259903736243,.7/max(float(abs(raw).max()),1e-30))
wavfile.write(root/'voice.wav',fs,np.rint(raw*gain*32767).astype(np.int16))
receipt=json.loads((root/'run.json').read_text());receipt.update(gain=gain,processing='one constant playback gain; no EQ, fade, crop or added output noise',enhanced_inspected=False)
receipt['audio_sha256']=hashlib.sha256((root/'voice.wav').read_bytes()).hexdigest()
receipt['raw_sha256']=hashlib.sha256((root/'voice-raw.wav').read_bytes()).hexdigest()
shutil.copy2(ROOT/'build'/receipt.get('renderer','ducklips_somatic_study'),root/'render')
receipt['executable_sha256']=hashlib.sha256((root/'render').read_bytes()).hexdigest()
paths=[ROOT/'CMakeLists.txt',*sorted((ROOT/'include/ducklips').glob('*.hpp')),*sorted((ROOT/'src').glob('*.cpp')),ROOT/'tools/render_somatic_study.cpp',ROOT/'tools/render_papa_ideas.cpp',Path(__file__)]
receipt['sources']={}
for p in paths:
 d=root/'snapshot'/p.relative_to(ROOT);d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,d);receipt['sources'][str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
(root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
a=np.genfromtxt(root/'trace.csv',names=True,delimiter=',')
fig,ax=plt.subplots(4,1,figsize=(13,10),sharex=True)
ax[0].plot(a['time'],a['lung']);ax[0].set_ylabel('breath, Pa')
ax[1].plot(a['time'],a['flow']*1e6);ax[1].set_ylabel('glottal flow, cm³/s')
ax[2].plot(a['time'],a['area']*1e6);ax[2].set_ylabel('glottal area, mm²')
ax[3].plot(np.arange(len(raw))/fs,raw[:,0]*128);ax[3].set_ylabel('microphone pressure, Pa');ax[3].set_xlabel('s')
fig.tight_layout();fig.savefig(root/'physical-trace.png',dpi=110)
b=np.genfromtxt(root/'control.csv',names=True,delimiter=',',dtype=None,encoding='utf8')
fig,ax=plt.subplots(5,2,figsize=(14,13),sharex=True)
for axis,key in zip(ax.ravel(),['tongue_x','tongue_y','jaw','blade','lip_area','lip_protrusion','velum','glottal_area','breath_pressure','laryngeal_pitch_intent']):
 rows=b[b['coordinate']==key];axis.plot(rows['time'],rows['target'],label='intent');axis.plot(rows['time'],rows['position'],label='actual');axis.plot(rows['time'],rows['sensed'],label='sensed',alpha=.6);axis.set_title(key)
ax[0,0].legend();fig.tight_layout();fig.savefig(root/'body-control.png',dpi=110)
print('Fresh audio and physical depictions ready:',root,flush=True)
subprocess.run([sys.executable,str(ROOT/'tools/enhanced_specgram/render.py'),str(root/'voice.wav'),'--channel','0','--max-hz','8000','--hop','128','--display-reference',str(.06376754775627493*gain/315.98259903736243),'--out',str(root/'enhanced')],check=True)
r=json.loads((root/'enhanced.json').read_text());assert r['source_sha256']==receipt['audio_sha256'] and r['end_sample']==len(raw) and r['start_sample']==0
receipt['enhanced_complete']=True;(root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
