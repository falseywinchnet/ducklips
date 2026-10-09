#!/usr/bin/env python3
"""Verify native continuous-mode integration against the canonical network."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import numpy as np
from moving_modal_reference import operator
from nonuniform_metric_reference import MovingNetwork,integrate
from interaction_modal_reference import step


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=Path('build/ducklips_continuous_modes_probe'))
    p.add_argument('--transport-binary',type=Path,default=Path('build/ducklips_modal_transport_probe'))
    p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True);binary=args.binary.resolve()
    network=MovingNetwork();dt=1/48000;frames=round(network.duration/dt);n=network.n
    operators=[operator(network,i*dt/2) for i in range(2*frames+1)]
    # Derive every half-sample coefficient with the native physical-rate
    # builder, then feed THOSE coefficients to the native time integrator.
    transport_binary=args.transport_binary.resolve()
    transport_input=args.out/'transport.fixture';transport_output=args.out/'transport.f64'
    edges=[(0,1),(1,2),(2,3),(2,4),(3,2**64-1)]+[(i,2**64-1) for i in range(n)]
    with transport_input.open('wb') as f:
        f.write(np.array([n,n,network.m,len(operators)],dtype='<u8').tobytes())
        f.write(np.array(edges,dtype='<u8').tobytes())
        for i,(unused,v,m,w) in enumerate(operators):
            c,g,rc,rg=network.parameters(i*dt/2)
            phi=(v[:n,:n]/np.sqrt(network.c0*c)[:,None]).T
            f.write(np.r_[w,phi.ravel(),network.c0*c*rc,network.g0*g*rg].astype('<f8').tobytes())
    subprocess.run([str(transport_binary),'--fixture',str(transport_input),str(transport_output)],check=True)
    derived=np.fromfile(transport_output,dtype='<f8').reshape(len(operators),n+2*n*n)
    expected=np.array([np.r_[w*(np.diag(m[:n,:n])+np.diag(m[n:,n:])),m[:n,:n].ravel(),m[n:,n:].ravel()]
                       for unused,unused2,m,w in operators])
    coefficient_error=float(np.max(abs(derived-expected)/np.maximum(1,abs(expected))))
    packed=[]
    for i in range(frames):
        record=[]
        for j in range(2*i,2*i+3):
            record.extend([operators[j][3],derived[j,n:n+n*n],derived[j,n+n*n:]])
        packed.append(np.concatenate(record).astype('<f8').tobytes())
    comparisons=[]
    for label in ['low','high','mixed']:
        initial=np.zeros(2*n)
        mode=n-1 if label=='high' else 0;initial[mode]=.3;initial[n+mode]=.7
        if label=='mixed':initial=np.random.default_rng(7109).normal(0,.2,2*n)
        network.initial=operators[0][1]@initial
        times,reference,reference_report=integrate(network)
        source=args.out/f'{label}.fixture';output=args.out/f'{label}.f64'
        with source.open('wb') as f:
            f.write(np.array([n,frames],dtype='<u8').tobytes());f.write(np.array([dt],dtype='<f8').tobytes())
            f.write(initial.astype('<f8').tobytes())
            for record in packed:f.write(record)
        subprocess.run([str(binary),'--fixture',str(source),str(output)],check=True)
        modal=np.vstack([initial,np.fromfile(output,dtype='<f8').reshape(frames,2*n)])
        physical=np.array([operators[2*i][1]@x for i,x in enumerate(modal)])
        python=[initial];x=initial.copy()
        for i in range(frames):
            samples=[(entry[2],entry[3]) for entry in operators[2*i:2*i+3]]
            x=step(x,*samples,dt);python.append(x)
        error=float(np.max(np.linalg.norm(physical-reference,axis=1))/np.linalg.norm(network.initial))
        python_error=float(np.max(np.linalg.norm(modal-np.array(python),axis=1)))
        comparisons.append({'initial_condition':label,'max_relative_state_error':error,
                            'max_absolute_native_vs_python_error':python_error,'reference':reference_report})
    timings=[json.loads(subprocess.check_output([str(binary),'--benchmark',str(n)])) for n in [20,24,32] for _ in range(3)]
    construction_timings=[json.loads(subprocess.check_output([str(transport_binary),'--benchmark',str(n)])) for n in [20,24,32] for _ in range(3)]
    summed_cost=[]
    for modes in [20,24,32]:
        kernel=max(t['wall_seconds'] for t in timings if t['modes']==modes)
        construction=max(t['wall_seconds'] for t in construction_timings if t['modes']==modes)
        summed_cost.append({'modes':modes,'sum_of_separate_worst_wall_seconds':kernel+construction,
                            'equivalent_realtime_multiple':5/(kernel+construction)})
    sources=['include/ducklips/modal_transport.hpp','src/modal_transport.cpp','src/real_dot.hpp','tools/modal_transport_probe.cpp',
             'include/ducklips/continuous_modes.hpp','src/continuous_modes.cpp','tools/continuous_modes_probe.cpp',
             'tools/interaction_modal_reference.py','tools/verify_continuous_modes.py','tools/moving_modal_reference.py']
    report={'schema':1,'scope':'Unforced continuous modal transport kernel; includes three coefficient fills and phase evaluation per sample; excludes eigenbasis/physical transport construction, losses, source and microphone',
            'platform':platform.platform(),'binary_sha256':sha(binary),'source_sha256':{s:sha(s) for s in sources},
            'native_transport_max_component_relative_error_with_unit_floor':coefficient_error,
            'transport_binary_sha256':sha(transport_binary),'comparisons':comparisons,'timings':timings,
            'transport_construction_timings':construction_timings,
            'summed_cost_not_complete_speech_benchmark':summed_cost}
    (args.out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    if coefficient_error>1e-10:raise SystemExit('Native physical transport coefficient gate failed')
    if any(c['max_relative_state_error']>1e-6 or c['max_absolute_native_vs_python_error']>1e-11 for c in comparisons):
        raise SystemExit('Continuous modes accuracy gate failed')
    if min(t['realtime_multiple'] for t in timings)<10 or min(t['equivalent_realtime_multiple'] for t in summed_cost)<10:
        raise SystemExit('Continuous modes speed below 10x: investigate representation')

if __name__=='__main__':main()
