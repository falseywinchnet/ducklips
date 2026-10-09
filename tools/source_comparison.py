#!/usr/bin/env python3
"""Compare fold geometry laws at the glottal-flow-derivative tap; no output EQ."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import soundfile as sf
from scipy.signal import welch

ROOT = Path(__file__).resolve().parents[1]
output = ROOT / 'out/source_comparison'
output.mkdir(parents=True, exist_ok=True)
rows = []
for name, body in [('sharp','sharp-contact-control'),('rounded','reference')]:
    path = output / (name + '.wav')
    command = [str(ROOT/'build/ducklips_say'), '--body', str(ROOT/'bodies'/f'{body}.body'),
               '--phonemes','AA1:1.5','--tap','source','--out',str(path)]
    result = subprocess.run(command,capture_output=True,text=True,check=True)
    (output / (name+'.log')).write_text(result.stdout)
    x, fs = sf.read(path)
    x = x[round(.9*fs):round(1.7*fs),0]
    f,p = welch(x,fs,nperseg=8192)
    reference = p[(f>=100)&(f<1000)].sum()
    bands = {f'{a}-{b}':float(10*np.log10(p[(f>=a)&(f<b)].sum()/reference))
             for a,b in [(1000,4000),(4000,8000),(8000,16000),(16000,22000)]}
    rows.append(dict(name=name,body=body,bands_db=bands,wav_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
report = dict(method='Welch spectrum, 8192-sample segments, 0.9-1.7 s excerpt, channel 0, relative to 100-1000 Hz power',
              measure='glottal flow derivative, before tract radiation; same body, cues and seed except aperture geometry',
              not_a='human realism, intelligibility, or calibrated sound-level score',cases=rows,
              change_db={key:rows[1]['bands_db'][key]-rows[0]['bands_db'][key] for key in rows[0]['bands_db']})
(ROOT/'reports/source_geometry.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
