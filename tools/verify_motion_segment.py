#!/usr/bin/env python3
"""Verify native physical knots -> canonical interpolation -> acoustic state."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import numpy as np
from moving_modal_reference import operator
from nonuniform_metric_reference import MovingNetwork,integrate


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--binary',type=Path,default=Path('build/ducklips_motion_segment_probe'))
    p.add_argument('--transport-binary',type=Path,default=Path('build/ducklips_modal_transport_probe'))
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    binary=args.binary.resolve();transport=args.transport_binary.resolve()
    network=MovingNetwork();n=network.n;interval=.001;cells=round(network.duration/interval);samples=48
    knot_times=np.arange(-2,cells+3)*interval
    operators=[operator(network,t) for t in knot_times]
    source=args.out/'physical-knots.fixture';output=args.out/'physical-knots.f64'
    edges=[(0,1),(1,2),(2,3),(2,4),(3,2**64-1)]+[(i,2**64-1) for i in range(n)]
    with source.open('wb') as f:
        f.write(np.array([n,n,network.m,len(operators)],dtype='<u8').tobytes());f.write(np.array(edges,dtype='<u8').tobytes())
        for t,(unused,v,m,w) in zip(knot_times,operators):
            c,g,rc,rg=network.parameters(t);phi=(v[:n,:n]/np.sqrt(network.c0*c)[:,None]).T
            f.write(np.r_[w,phi.ravel(),network.c0*c*rc,network.g0*g*rg].astype('<f8').tobytes())
    subprocess.run([str(transport),'--fixture',str(source),str(output)],check=True)
    native=np.fromfile(output,dtype='<f8').reshape(len(operators),n+2*n*n)
    records=[np.r_[entry[3],derived[:n],derived[n:n+n*n]].astype('<f8').tobytes() for entry,derived in zip(operators,native)]
    comparisons=[]
    for label in ['low','high','mixed']:
        initial=np.zeros(2*n);mode=n-1 if label=='high' else 0;initial[mode]=.3;initial[n+mode]=.7
        if label=='mixed':initial=np.random.default_rng(7109).normal(0,.2,2*n)
        network.initial=operators[2][1]@initial
        times,reference,reference_report=integrate(network)
        source=args.out/f'{label}.fixture';output=args.out/f'{label}.f64'
        with source.open('wb') as f:
            f.write(np.array([n,cells,samples],dtype='<u8').tobytes());f.write(np.array([interval],dtype='<f8').tobytes())
            f.write(initial.astype('<f8').tobytes())
            for record in records:f.write(record)
        subprocess.run([str(binary),'--fixture',str(source),str(output)],check=True)
        state=np.vstack([initial,np.fromfile(output,dtype='<f8').reshape(cells*samples,2*n)])
        physical=np.array([operator(network,t)[1]@x for t,x in zip(times,state)])
        error=float(np.max(np.linalg.norm(physical-reference,axis=1))/np.linalg.norm(network.initial))
        comparisons.append({'initial_condition':label,'max_relative_state_error':error,'reference':reference_report})
    timings=[json.loads(subprocess.check_output([str(binary),'--benchmark',str(n)])) for n in [20,24,32] for _ in range(3)]
    construction=[json.loads(subprocess.check_output([str(transport),'--benchmark',str(n)])) for n in [20,24,32] for _ in range(3)]
    budget=[]
    for n in [20,24,32]:
        seconds=max(t['wall_seconds'] for t in timings if t['modes']==n)+max(t['wall_seconds'] for t in construction if t['modes']==n)
        budget.append({'modes':n,'sum_of_separate_worst_wall_seconds':seconds,'equivalent_realtime_multiple':5/seconds})
    files=['include/ducklips/continuous_modes.hpp','src/continuous_modes.cpp','src/real_dot.hpp',
           'include/ducklips/modal_transport.hpp','src/modal_transport.cpp','include/ducklips/motion_segment.hpp',
           'src/motion_segment.cpp','tools/motion_segment_probe.cpp','tools/modal_transport_probe.cpp',
           'tools/verify_motion_segment.py','tools/canonical_motion_reference.py']
    report={'schema':1,'scope':'Native physical-rate construction at 1kHz, compact coefficient interpolation, and unforced acoustic stepping; isolated gradient modes. Timing excludes geometry/eigenbasis solving, source, physical losses and microphone.',
            'lookahead_seconds':.003,'command_scheduler':'Not implemented; diagnostic supplies future knots',
            'platform':platform.platform(),'binary_sha256':sha(binary),'transport_binary_sha256':sha(transport),
            'source_sha256':{f:sha(f) for f in files},'comparisons':comparisons,
            'interpolation_and_step_timings':timings,'physical_rate_construction_timings':construction,
            'summed_cost_not_complete_speech_benchmark':budget}
    (args.out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if any(c['max_relative_state_error']>1e-6 for c in comparisons):raise SystemExit('Native motion segment accuracy gate failed')
    if min(t['realtime_multiple'] for t in timings)<10 or min(t['equivalent_realtime_multiple'] for t in budget)<10:
        raise SystemExit('Native motion segment cost below 10x: investigate representation')

if __name__=='__main__':main()
