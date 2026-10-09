#!/usr/bin/env python3
"""Measure the control-clock tones identified in the full enhanced speech field.

Uses windowed sinusoid projections of the saved PCM, not a companion image.
No filtering, normalization, or synthesis changes are applied to the recordings.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.io import wavfile
from scipy.signal import find_peaks


def analyze(path):
    fs,pcm=wavfile.read(path)
    if pcm.ndim==2:pcm=pcm[:,0]
    if np.issubdtype(pcm.dtype,np.integer):
        if pcm.dtype==np.uint8:pcm=(pcm.astype(float)-128)/128
        else:pcm=pcm.astype(float)/max(abs(np.iinfo(pcm.dtype).min),np.iinfo(pcm.dtype).max)
    pcm=pcm.astype(float)
    if not np.isfinite(pcm).all():raise ValueError('Nonfinite PCM')
    windows=[]
    for start,end in [(.16,.28),(.4,.52),(.48,.58),(1.05,1.2)]:
        x=pcm[round(start*fs):round(end*fs)]
        if len(x)!=round(end*fs)-round(start*fs):raise ValueError('Audio too short')
        w=np.hanning(len(x));times=np.arange(len(x))/fs
        amplitudes=[float(2*abs(np.sum(x*w*np.exp(-2j*np.pi*f*times)))/w.sum()) for f in range(1000,8001,1000)]
        spectrum=abs(np.fft.rfft(x*w,131072));hz=np.fft.rfftfreq(131072,1/fs)
        peaks=find_peaks(spectrum)[0];peaks=peaks[(hz[peaks]>500)&(hz[peaks]<8000)]
        peaks=peaks[np.argsort(spectrum[peaks])[-8:][::-1]]
        windows.append({'seconds':[start,end],'tone_hz':list(range(1000,8001,1000)),
                        'tone_peak_amplitudes':amplitudes,
                        'tone_dbfs':[float(20*np.log10(max(v,1e-30))) for v in amplitudes],
                        'strongest_peak_hz':hz[peaks].tolist(),'rms':float(np.sqrt(np.mean(x*x)))})
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'sample_rate':fs,'channel':0,'windows':windows}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('wav',type=Path,nargs='+');p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    report={'method':'Hann-windowed sinusoid projection; peak amplitude 2*abs(sum(w*x*exp(-jwt)))/sum(w), dB relative to full-scale amplitude 1',
            'observations':[analyze(path) for path in a.wav]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2)+'\n')
    for r in report['observations']:print(r['path'],r['windows'][0]['tone_dbfs'][:2])

if __name__=='__main__':main()
