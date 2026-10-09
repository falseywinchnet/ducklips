#!/usr/bin/env python3
"""Diffuse-wall, linearized isothermal BGK plane-channel conductance.

Cell-averaged integral transport, with exact spatial exponential integrals.
This is an offline constitutive reference, not a molecular speech simulation.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import solve_toeplitz,matmul_toeplitz
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/slit-kinetic'

def velocities(order=16):
    x,w=leggauss(order);edges=np.arange(-42.,3.,1.);cs=[];ws=[]
    for a,b in zip(edges,edges[1:]):
        c=np.exp((a+b)/2+(b-a)/2*x);weight=w*(b-a)/2*c*np.exp(-c*c)/np.sqrt(np.pi);cs.extend(c);ws.extend(weight)
    return np.array(cs),np.array(ws)

def conductance(delta,cells=512,order=16,details=False):
    c,w=velocities(order);t=delta/(cells*c);one=-np.expm1(-t)
    diagonal=np.empty_like(t);small=t<.01;z=t[small];diagonal[small]=z*(.5+z*(-1/6+z*(1/24+z*(-1/120+z*(1/720-z/5040)))))
    diagonal[~small]=1-one[~small]/t[~small]
    column=np.empty(cells);column[0]=np.sum(2*w*diagonal)
    column[1:]=(w*one*one/t)@np.exp(-np.outer(t,np.arange(cells-1)))
    # K is the cell-average transport kernel. Positive velocity V=-2*u_P
    # solves (I-K)V=K*1/delta, and G is its cross-sectional mean.
    cumulative=np.cumsum(column);row_sum=cumulative+cumulative[::-1]-column[0]
    matrix=-column;matrix[0]+=1
    v=solve_toeplitz(matrix,row_sum/delta,check_finite=True);g=float(np.mean(v))
    if not details:return g
    residual=matmul_toeplitz(matrix,v)-row_sum/delta
    return dict(delta=delta,cells=cells,velocity_order=order,conductance=g,minimum_velocity=float(min(v)),symmetry_error=float(max(abs(v-v[::-1]))/max(v)),linear_residual=float(max(abs(residual))/max(abs(row_sum/delta))),gaussian_weight_sum=float(sum(w)))

def low_limit(delta):
    return (-np.log(delta)+1.5*(1-np.euler_gamma))/np.sqrt(np.pi)

def high_limit(delta):
    # Sharipov & Seleznev (1998), Eq. 3.88, p.691, BGK diffuse wall.
    return delta/6+1.0162+1.0653/delta-2.1354/delta**2

def refined(delta):
    runs=[dict(cells=n,conductance=conductance(delta,n)) for n in [512,1024,2048]]
    while abs(runs[-1]['conductance']/runs[-2]['conductance']-1)>1e-5:
        n=runs[-1]['cells']*2
        if n>65536:raise RuntimeError('spatial refinement budget exceeded')
        runs.append(dict(cells=n,conductance=conductance(delta,n)))
    g=runs[-1]['conductance'];velocity=conductance(delta,runs[-1]['cells'],24)
    return dict(delta=delta,conductance=g,runs=runs,spatial_relative_change=abs(g/runs[-2]['conductance']-1),velocity_relative_change=abs(velocity/g-1))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='reports/slit_kinetic_reference.json');args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    published={.1:2.0327,.5:1.6018,1.:1.5386,2.:1.5948,5.:1.9907,10.:2.7686}
    reference=[]
    for delta,target in published.items():
        runs=[conductance(delta,n,details=True) for n in [256,512,1024,2048]];reference.append(dict(delta=delta,published=target,runs=runs,relative_error=abs(runs[-1]['conductance']/target-1)))
    deltas=np.geomspace(1e-8,100.,201);table=[]
    for delta in deltas:
        table.append(refined(float(delta)))
        if len(table)%20==0:print('table',len(table),table[-1],flush=True)
    # Positive C2 interpolant in log(G), with analytic tail value/slope joins.
    values=np.log([row['conductance'] for row in table]);step=np.log(deltas[1]/deltas[0])
    values[0]=np.log(low_limit(deltas[0]));values[-1]=np.log(high_limit(deltas[-1]))
    left_slope=-step/(np.sqrt(np.pi)*low_limit(deltas[0]));d=deltas[-1]
    right_slope=step*(d/6-1.0653/d+4.2708/d**2)/high_limit(d)
    spline=CubicSpline(np.arange(len(table)),values,bc_type=((1,left_slope),(1,right_slope)))
    checks=[]
    for i in range(len(table)-1):
        for fraction in [.25,.5,.75]:
            delta=float(np.exp(np.log(deltas[0])+(i+fraction)*step));r=refined(delta);actual=r['conductance'];value=float(np.exp(spline(i+fraction)))
            checks.append(dict(delta=delta,reference=actual,interpolated=value,relative_error=abs(value/actual-1),cells=r['runs'][-1]['cells'],spatial_relative_change=r['spatial_relative_change']))
        if (i+1)%40==0:print('interpolation',i+1,flush=True)
    asymptote=[dict(delta=float(d),actual=conductance(float(d),1024),leading=low_limit(d)) for d in [1e-4,1e-6,1e-8]]
    high=[]
    for d in [100.,200.,500.]:
        runs=[dict(cells=n,conductance=conductance(d,n)) for n in [8192,16384,32768]]
        extrapolated=(4*runs[-1]['conductance']-runs[-2]['conductance'])/3
        high.append(dict(delta=d,runs=runs,richardson=extrapolated,asymptote=high_limit(d),relative_error=abs(high_limit(d)/extrapolated-1)))
    report=dict(scope='linearized isothermal BGK gas, stationary diffuse parallel walls, small longitudinal pressure gradient; not a validated human-air surface model',published_source='https://srd.nist.gov/jpcrdreprint/1.556019.pdf',published_location='Table 1 column f page 679; high-delta Eq 3.88 page 691',initial_failure='slit_kinetic_initial_grid_failure.json: original 2048-cell limit failed 3e-4 spatial change at delta 68.13 and 100; retained, followed by refinement to stricter 1e-5',benchmarks=reference,table=table,interpolation_checks=checks,free_molecular=asymptote,hydrodynamic=high,versions=dict(numpy=np.__version__,scipy=scipy.__version__),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/args.out).write_text(json.dumps(report,indent=2)+'\n')
    assert max(r['relative_error'] for r in reference)<2e-4
    assert max(r['spatial_relative_change'] for r in table)<1e-5
    assert max(r['velocity_relative_change'] for r in table)<1e-9
    assert max(r['relative_error'] for r in checks)<2e-5
    assert max(r['relative_error'] for r in high)<2e-5
    assert abs(asymptote[-1]['actual']/asymptote[-1]['leading']-1)<1e-6
    rows=['// Generated by tools/slit_kinetic_reference.py; log(G), 200 cubic intervals.',
          '// BGK diffuse walls; independently checked at three interior points per interval.',
          'constexpr double bgk_log_min = '+format(float(np.log(deltas[0])),'.17g')+';',
          'constexpr double bgk_log_step = '+format(float(step),'.17g')+';',
          'constexpr double bgk_coefficients[200][4] = {']
    rows += ['    {'+', '.join(format(float(v),'.17g') for v in spline.c[:,i])+'},' for i in range(200)]
    rows += ['};','']
    (ROOT/'src/slit_bgk_table.inc').write_text('\n'.join(rows))
    print('PASS benchmarks, spatial/velocity refinement, 600 withheld interpolation points, asymptotic checks',flush=True)
if __name__=='__main__':main()
