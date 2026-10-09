#!/usr/bin/env python3
"""Compare the native moving tree with a directly integrated nonuniform network.

Fixture exchange uses little-endian uint64 headers/edges and float64 samples.
The diagnostic executable is local to the same host; no audio is produced.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import numpy as np
from scipy.linalg import null_space
from nonuniform_metric_reference import MovingNetwork, integrate


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_fixture(binary, folder, network, subdivisions, label):
    n=network.n;m=network.m;dt=1/(48000*subdivisions)
    frames=round(network.duration/dt)
    edges=[(0,1),(1,2),(2,3),(2,4),(3,2**64-1)]+[(i,2**64-1) for i in range(n)]
    input_path=folder/f'{label}-{subdivisions}.fixture'
    output_path=folder/f'{label}-{subdivisions}.f64'
    with input_path.open('wb') as f:
        f.write(np.array([n,m,frames],dtype='<u8').tobytes());f.write(np.array([dt],dtype='<f8').tobytes())
        f.write(np.array(edges,dtype='<u8').tobytes())
        initial=np.r_[network.initial[:n]*np.sqrt(network.c0),network.initial[n:]/np.sqrt(network.g0)]
        f.write(initial.astype('<f8').tobytes())
        for frame in range(frames):
            c1,g1,_,_=network.parameters((frame+.5-np.sqrt(3)/6)*dt)
            c2,g2,_,_=network.parameters((frame+.5+np.sqrt(3)/6)*dt)
            f.write(np.r_[network.c0*c1,network.c0*c2,network.g0*g1,network.g0*g2].astype('<f8').tobytes())
    subprocess.run([str(binary),'--fixture',str(input_path),str(output_path)],check=True)
    raw=np.fromfile(output_path,dtype='<f8').reshape(frames,n+m)[subdivisions-1::subdivisions]
    result=[network.initial]
    for i,state in enumerate(raw):
        c,g,_,_=network.parameters((i+1)/48000)
        result.append(np.r_[state[:n]/np.sqrt(network.c0*c),state[n:]*np.sqrt(network.g0*g)])
    return np.array(result)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=Path('build/ducklips_moving_tree_probe'))
    p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True);binary=args.binary.resolve()
    comparisons=[]
    for label in ['gradient','circulation','highest_mode']:
        network=MovingNetwork()
        if label=='highest_mode':network.initial=np.r_[.3*network.pressure_modes[:,-1],.7*network.flow_modes[:,-1]]
        if label=='circulation':network.initial[network.n:]+=.2*null_space(network.k0)[:,0]
        times,reference,report=integrate(network)
        previous=None;cases=[]
        for subdivisions in [1,2,4]:
            native=run_fixture(binary,args.out,network,subdivisions,label)
            error=float(np.max(np.linalg.norm(native-reference,axis=1))/np.linalg.norm(network.initial))
            cases.append({'substeps':subdivisions,'max_relative_state_error':error,
                          'refinement_ratio':None if previous is None else previous/error})
            previous=error
        comparisons.append({'initial_condition':label,'reference':report,'cases':cases})
    # Do not overlap these timings with builds or numerical reference solves.
    timings=[json.loads(subprocess.check_output([str(binary),'--benchmark',str(substeps)]))
             for substeps in [1,2,4] for _ in range(3)]
    sources=['include/ducklips/moving_tree.hpp','src/moving_tree.cpp','tools/moving_tree_probe.cpp',
             'tools/nonuniform_metric_reference.py','tools/verify_moving_tree.py']
    phase_bias=[]
    for hz in [651,1234,3264,8000,12000]:
        z=2j*np.pi*hz/48000
        response=(1+z/2+z*z/12)/(1-z/2+z*z/12)
        actual=float(np.angle(response)*48000/(2*np.pi))
        phase_bias.append({'frequency_hz':hz,'numerical_frequency_hz':actual,'relative_frequency_error':actual/hz-1})
    result={'schema':1,'scope':'Lossless moving-tree kernel only; not a complete speech-rendering speed claim',
            'platform':platform.platform(),'binary_sha256':digest(binary),
            'source_sha256':{s:digest(s) for s in sources},'comparisons':comparisons,'stationary_phase_bias_at_48khz':phase_bias,'timings':timings}
    (args.out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    for comparison in comparisons:
        cases=comparison['cases']
        if any(c['refinement_ratio']<14 for c in cases[1:]):
            raise SystemExit('Moving-tree fourth-order convergence gate failed')
        if comparison['initial_condition']!='highest_mode' and cases[-1]['max_relative_state_error']>1e-6:
            raise SystemExit('Low-mode moving-tree accuracy gate failed')
    if min(t['realtime_multiple'] for t in timings)<10:
        raise SystemExit('Moving-tree speed below 10x: investigate representation')

if __name__=='__main__':main()
