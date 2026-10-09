#!/usr/bin/env python3
"""Coarse timing overview only. Use native_fields.py to inspect field detail.

Nearest interpolation does not prevent row loss when the axes or the displayed
figure are smaller than the input raster. This overview is deliberately not an
ordinary-versus-enhanced resolution comparison.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('human', type=Path)
p.add_argument('generated', type=Path)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
fig, axes = plt.subplots(3, 2, figsize=(14, 9), layout='constrained')
for column, (path, label) in enumerate([(a.human,'Human: valeatory — Alfa'), (a.generated,'DuckLips: reference body — Alfa')]):
    with np.load(path) as data:
        fs = int(data['sample_rate'])
        t = data['seconds']
        samples = data['samples']
        peak = max(np.max(np.abs(samples)),1e-30)
        axes[0,column].plot(t[0]+np.arange(len(samples))/fs,samples/peak,lw=.45,color='#245e85')
        axes[0,column].set(title=label,ylabel='Waveform / own peak',ylim=(-1.1,1.1))
        for row,key,hz_key in [(1,'ordinary_stft_magnitude','ordinary_hz'),(2,'enhanced','hz')]:
            field = data[key]
            hz = data[hz_key]
            db = 20*np.log10(np.maximum(field/max(np.max(field),1e-30),1e-3))
            axes[row,column].imshow(db,origin='lower',aspect='auto',interpolation='nearest',cmap='magma',
                                  vmin=-60,vmax=0,extent=[t[0],t[-1],hz[0],hz[-1]])
            axes[row,column].set(ylim=(0,8000),ylabel='Hz',title='Ordinary STFT' if row==1 else 'Full enhanced field')
        for ax in axes[:,column]:
            ax.set_xlim(0,1.5)
            ax.set_xlabel('Seconds from recording start (no time warping)')
fig.suptitle('Timing overview only — spatial detail reduced\nUse native_fields.py for ordinary/enhanced detail; independent display peaks',fontsize=15)
a.out.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(a.out,dpi=150)
