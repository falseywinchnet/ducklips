#!/usr/bin/env python3
"""Compare native frozen-modal impulses to dense equations; time coupled folds."""
import argparse
import hashlib
import json
import platform
import subprocess
import tempfile
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('model',type=Path);p.add_argument('modes',type=Path)
    p.add_argument('--binary',type=Path,default=Path('build/surface/surface_modal_probe'))
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    with np.load(args.model) as data:
        a,b,c,d,o,ou=[data[k] for k in ['A','B','C','D','O','OU']];fs=float(data['fs'])
    outputs=np.vstack([c,o]);direct=np.vstack([d,ou])
    duration=.0625;damping=400.0;samples=int(duration*fs)
    bins=np.unique(np.round(np.linspace(50,12000,241)*samples/fs).astype(int))
    hz=bins*fs/samples
    expected=np.array([outputs@np.linalg.solve(np.exp((damping+2j*np.pi*f)/fs)*np.eye(len(a))-a,b)+direct for f in hz])
    observed=np.empty_like(expected)
    impulse_runs=[]
    with tempfile.TemporaryDirectory(prefix='ducklips-modal-') as temp:
        for source in range(b.shape[1]):
            path=Path(temp)/f'{source}.f64'
            receipt=json.loads(subprocess.check_output([str(args.binary),str(args.modes),str(duration),'impulse',str(source),str(path)],text=True))
            impulse_runs.append(receipt)
            data=np.fromfile(path,dtype=np.float64).reshape(samples,outputs.shape[0])
            spectrum=np.fft.rfft(data*np.exp(-damping*np.arange(samples)/fs)[:,None],axis=0)/1e-10
            observed[:,:,source]=spectrum[bins]
    peak=abs(expected).max(axis=0)
    error=abs(observed-expected)/np.maximum(abs(expected),np.maximum(peak*.01,1e-30))
    runs=[json.loads(subprocess.check_output([str(args.binary),str(args.modes),'10','voice'],text=True)) for unused in range(5)]
    port_runs=[json.loads(subprocess.check_output([str(args.binary),str(args.modes),'10','ports'],text=True)) for unused in range(5)]
    report={'schema':1,'platform':platform.platform(),
            'binary_sha256':hashlib.sha256(args.binary.read_bytes()).hexdigest(),
            'model_sha256':hashlib.sha256(args.model.read_bytes()).hexdigest(),
            'modes_sha256':hashlib.sha256(args.modes.read_bytes()).hexdigest(),
            'git_head_at_measurement':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            'working_tree_modified':bool(subprocess.check_output(['git','status','--porcelain'],text=True).strip()),
            'impulse_duration_seconds':duration,'laplace_damping_per_second':damping,
            'frequency_count':len(hz),'frequency_range_hz':[float(hz.min()),float(hz.max())],
            'max_native_error_with_1pct_per_transfer_floor':float(error.max()),
            'per_transfer_native_error_with_1pct_floor':error.max(axis=0).tolist(),
            'impulse_runs':impulse_runs,'coupled_voice_runs':runs,'all_port_runs':port_runs,
            'minimum_all_ports_realtime_multiple':min(r['realtime_multiple'] for r in port_runs),
            'minimum_coupled_voice_realtime_multiple':min(r['realtime_multiple'] for r in runs),
            'accepted_as_speech_renderer':False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+'\n')
    print('Native transfer error:',report['max_native_error_with_1pct_per_transfer_floor'])
    print('Coupled voice speed:',[r['realtime_multiple'] for r in runs])
    if not np.isfinite(error).all() or error.max()>1e-5:raise SystemExit('Native coordinate verification failed')
    print('All port speed:',[r['realtime_multiple'] for r in port_runs])
    if min(report['minimum_coupled_voice_realtime_multiple'],report['minimum_all_ports_realtime_multiple'])<10:raise SystemExit('Below 10x: investigate representation before acceptance')

if __name__=='__main__':main()
