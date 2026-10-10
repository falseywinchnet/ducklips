#!/usr/bin/env python3
"""Supplement a complete enhanced speech depiction with its source-state trace.

Does not edit sound or replace/crop the full enhanced field. The detail view is
explicitly a waveform diagnostic; the whole utterance is always saved too.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('round', type=Path)
    parser.add_argument('--compare', type=Path)
    parser.add_argument('--detail', nargs=2, type=float, default=[1.0, 1.035])
    args = parser.parse_args()
    a = np.genfromtxt(args.round/'trace.csv', delimiter=',', names=True)
    t = a['time']
    # Trace timestamps are rounded for CSV output. Sample spacing is defined
    # by the audio receipt, not by differences of rounded timestamps.
    import json
    fs = json.loads((args.round/'round.json').read_text())['playback']['sample_rate']
    fig, ax = plt.subplots(4, 1, figsize=(13, 9), sharex=True)
    ax[0].plot(t, a['flow']*1e6, label='transferred mean Q')
    ax[0].plot(t, a['glottal_endpoint_flow']*1e6, label='endpoint Q', alpha=.7)
    ax[0].set_ylabel('cm³/s')
    for k in ['lower_fold_area', 'upper_fold_area', 'area']:
        ax[1].plot(t, a[k]*1e6, label=k)
    ax[1].set_ylabel('mm²')
    ax[2].plot(t, a['p0'], label='inlet pressure')
    ax[2].plot(t, a['lung_pressure'], label='lung')
    ax[2].set_ylabel('Pa')
    ax[3].plot(t, np.gradient(a['flow'])*fs, label='Q derivative')
    ax[3].plot(t, a['lips'], label='lip derivative', alpha=.65)
    ax[3].set_ylabel('m³/s²')
    ax[3].set_xlabel('s')
    if args.compare:
        b = np.genfromtxt(args.compare/'trace.csv', delimiter=',', names=True)
        ax[0].plot(b['time'], b['flow']*1e6, ':', color='black', label=args.compare.name+' Q')
    for axis in ax:
        axis.legend(loc='upper left', fontsize=8)
    fig.suptitle(args.round.name+': complete source trace')
    fig.tight_layout()
    fig.savefig(args.round/'source-trace-full.png', dpi=100)
    for axis in ax:
        axis.set_xlim(*args.detail)
    fig.suptitle(args.round.name+': supplementary cycle detail; complete enhanced field retained')
    fig.savefig(args.round/'source-cycle-detail.png', dpi=100)


if __name__ == '__main__':
    main()
