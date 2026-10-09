#!/usr/bin/env python3
"""Independent canonical graph audit, including cycle flow and shared outlets.

No performance benchmark, body fit or speech render is performed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
A=np.array([[.25,.25-np.sqrt(3)/6],[.25+np.sqrt(3)/6,.25]])


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def incidence(n,edges):
    b=np.zeros((n,len(edges)))
    for k,(i,j) in enumerate(edges):
        b[i,k]=1
        if j>=0:b[j,k]=-1
    return b


def scaled_operator(b,c,g,scale):
    n,m=b.shape
    physical=np.block([[np.zeros((n,n)),-b*g],[b.T/c,np.zeros((m,m))]])
    return physical*scale[np.newaxis,:]/scale[:,np.newaxis]


def collocation(z,b,c1,c2,g1,g2,dt,nodal,series,scale):
    size=len(z);identity=np.eye(size);source=np.r_[nodal,series]/scale
    f1,f2=[scaled_operator(b,c,g,scale) for c,g in [(c1,g1),(c2,g2)]]
    matrix=np.block([[identity-dt*A[0,0]*f1,-dt*A[0,1]*f2],[-dt*A[1,0]*f1,identity-dt*A[1,1]*f2]])
    right=np.concatenate([z+dt*sum(row)*source for row in A]);stages=np.linalg.solve(matrix,right).reshape(2,size)
    result=z+dt*(.5*(f1@stages[0]+f2@stages[1])+source)
    physical=stages*scale;n=b.shape[0]
    return result,.5*(physical[0,:n]/c1+physical[1,:n]/c2),.5*(g1*physical[0,n:]+g2*physical[1,n:])


def run(name,binary,out,n,edges,initial,steps):
    source=out/f'{name}.txt';destination=source.with_suffix('.jsonl');log=source.with_suffix('.log')
    lines=[f'{n} {len(edges)} {len(steps)}',*[' '.join(map(str,e)) for e in edges],
           ' '.join(map(str,initial[:n])),' '.join(map(str,initial[n:]))]
    for dt,c1,c2,g1,g2,nodal,series in steps:
        lines.append(str(dt));lines.extend(' '.join(map(str,v)) for v in [c1,c2,g1,g2,nodal,series])
    source.write_text('\n'.join(lines)+'\n')
    command=[str(binary),'--fixture']
    with source.open() as f,destination.open('w') as g:result=subprocess.run(command,stdin=f,stdout=g,stderr=subprocess.PIPE,text=True)
    log.write_text(result.stderr);assert result.returncode==0,(name,result.stderr)
    states=[json.loads(line) for line in destination.read_text().splitlines()];assert len(states)==len(steps)
    receipt=dict(name=name,nodes=n,branches=len(edges),steps=len(steps),command=command,
                 **{k:str(v.relative_to(ROOT)) for k,v in [('input',source),('output',destination),('log',log)]},
                 **{k+'_sha256':sha(v) for k,v in [('input',source),('output',destination),('log',log)]})
    return states,receipt


def discrete_cases(binary,out,count):
    rng=np.random.default_rng(581093);receipts=[]
    for case in range(count):
        n=3+case%8
        if case%4==0:
            edges=[(i,i+1) for i in range(n-1) if i!=n//2-1]
        else:edges=[(i,(i+1)%n) for i in range(n)]
        if case%3==0:edges.append(edges[0])
        if case%4!=0:edges.extend([(0,n-1),(n-1,-1)])
        m=len(edges);b=incidence(n,edges)
        c0=10**rng.uniform(-11,-9,n);g0=10**rng.uniform(-4,-1.5,m);scale=np.r_[np.sqrt(c0),1/np.sqrt(g0)]
        initial=scale*rng.normal(0,.001,n+m);z=initial/scale
        steps=[]
        for k in range(20):
            c1=c0*np.exp(rng.uniform(-.2,.2,n));c2=c0*np.exp(rng.uniform(-.2,.2,n))
            g1=g0*np.exp(rng.uniform(-.2,.2,m));g2=g0*np.exp(rng.uniform(-.2,.2,m))
            if k%5==0:g1[0]=g2[0]=0
            steps.append((1./(48000*(1+case%3)),c1,c2,g1,g2,rng.normal(0,1e-5,n),rng.normal(0,10,m)))
        states,receipt=run(f'graph-{case:02d}',binary,out,n,edges,initial,steps)
        maxima=dict(pressure_Pa=0.,momentum_Pa_s=0.,average_pressure_Pa=0.,average_flow_m3_s=0.,local_mass_m3=0.,global_mass_m3=0.)
        old=initial.copy()
        for state,step in zip(states,steps):
            dt,c1,c2,g1,g2,nodal,series=step
            z,pavg,qavg=collocation(z,b,c1,c2,g1,g2,dt,nodal,series,scale);reference=z*scale
            actual=np.r_[state['charge'],state['momentum']]
            maxima['pressure_Pa']=max(maxima['pressure_Pa'],float(max(abs((actual[:n]-reference[:n])/c2))))
            maxima['momentum_Pa_s']=max(maxima['momentum_Pa_s'],float(max(abs(actual[n:]-reference[n:]))))
            maxima['average_pressure_Pa']=max(maxima['average_pressure_Pa'],float(max(abs(np.array(state['average_pressure'])-pavg))))
            maxima['average_flow_m3_s']=max(maxima['average_flow_m3_s'],float(max(abs(np.array(state['average_flow'])-qavg))))
            mass=actual[:n]-old[:n]+dt*(b@np.array(state['average_flow'])-nodal)
            maxima['local_mass_m3']=max(maxima['local_mass_m3'],float(max(abs(mass))))
            external=sum(state['average_flow'][k] for k,(_,j) in enumerate(edges) if j<0)
            maxima['global_mass_m3']=max(maxima['global_mass_m3'],abs(sum(actual[:n]-old[:n])+dt*(external-sum(nodal))))
            old=actual
        assert maxima['pressure_Pa']<1e-6 and maxima['momentum_Pa_s']<1e-8 and maxima['average_pressure_Pa']<1e-6,maxima
        assert maxima['average_flow_m3_s']<1e-10 and maxima['local_mass_m3']<1e-18 and maxima['global_mass_m3']<1e-18,maxima
        receipt['maxima']=maxima;receipts.append(receipt)
    return receipts


def continuous_cases(binary,out):
    # Two mouth regions connect both to each other and to the same outlet
    # volume. A common pressure has a different response from separate loads.
    n=4;edges=[(0,1),(0,2),(1,2),(2,3),(3,-1)]
    b=incidence(n,edges);c0=np.array([2e-10,3e-10,1e-11,2e-11]);g0=np.array([.01,.005,.007,.003,.002])
    scale=np.r_[np.sqrt(c0),1/np.sqrt(g0)];initial=scale*np.array([.001,-.0004,0,0,.0001,-.0002,.0003,0,0])
    nodal=np.array([2e-5,0,0,0]);series=np.array([0,0,0,0,-20.]);duration=.0004
    def parameters(t):
        return c0*np.exp(.15*np.sin(2*np.pi*t/duration+np.arange(n)*.3)),g0*np.exp(.2*np.sin(2*np.pi*t/duration+np.arange(len(edges))*.2))
    def derivative(t,z):
        c,g=parameters(t);return scaled_operator(b,c,g,scale)@z+np.r_[nodal,series]/scale
    exact=solve_ivp(derivative,[0,duration],initial/scale,method='DOP853',rtol=2e-12,atol=1e-15,dense_output=True,max_step=duration/256)
    assert exact.success
    # A tighter independent integration bounds numerical reference sensitivity.
    tighter=solve_ivp(derivative,[0,duration],initial/scale,method='DOP853',rtol=3e-13,atol=1e-16,max_step=duration/512)
    assert tighter.success
    oracle_delta=float(np.linalg.norm(exact.y[:,-1]-tighter.y[:,-1])/np.linalg.norm(initial/scale))
    assert oracle_delta<1e-10
    results=[];previous=None
    for subdivisions in [32,64,128]:
        dt=duration/subdivisions;steps=[]
        for k in range(subdivisions):
            c1,g1=parameters((k+.5-np.sqrt(3)/6)*dt);c2,g2=parameters((k+.5+np.sqrt(3)/6)*dt)
            steps.append((dt,c1,c2,g1,g2,nodal,series))
        states,receipt=run(f'shared-outlet-{subdivisions}',binary,out,n,edges,initial,steps)
        error=max(float(np.linalg.norm(np.r_[s['charge'],s['momentum']]/scale-exact.sol((i+1)*dt))/np.linalg.norm(initial/scale)) for i,s in enumerate(states))
        if previous is not None:assert previous/error>14,(previous,error)
        receipt.update(maximum_relative_state_error=error,refinement_ratio=None if previous is None else previous/error,
                       continuous_oracle_endpoint_sensitivity=oracle_delta)
        previous=error;results.append(receipt)
    assert previous<1e-6
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=ROOT/'build/ducklips_moving_graph_probe')
    p.add_argument('--out',type=Path,default=ROOT/'out/moving-graph')
    p.add_argument('--report',type=Path,default=ROOT/'reports/moving_graph_reference.json')
    p.add_argument('--sanitizer',action='store_true');args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    discrete=discrete_cases(args.binary,args.out,8 if args.sanitizer else 24)
    continuous=continuous_cases(args.binary,args.out)
    sources=['include/ducklips/moving_graph.hpp','src/moving_graph.cpp','tools/moving_graph_probe.cpp','tools/moving_graph_reference.py',
             'include/ducklips/moving_tree.hpp','src/moving_tree.cpp']
    report=dict(use='general-graph canonical reference including a common outlet; not a complete articulated lip or speech render',
                command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],discrete=discrete,continuous=continuous,
                hashes={str(path.relative_to(ROOT)):sha(path) for path in [*(ROOT/s for s in sources),args.binary]},
                versions=dict(numpy=np.__version__,scipy=scipy.__version__))
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(discrete_maxima={k:max(r['maxima'][k] for r in discrete) for k in discrete[0]['maxima']},continuous=[{k:r[k] for k in ['name','maximum_relative_state_error','refinement_ratio']} for r in continuous]),indent=2))


if __name__=='__main__':main()
