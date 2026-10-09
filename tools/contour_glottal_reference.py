#!/usr/bin/env python3
"""Independent open-contour energy derivatives and continuous coupled motion."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import scipy
import mpmath as mp
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

# Independent interval integration, also usable with arbitrary precision or
# complex-step differentiation. The native paired-number code is not called.
def geometry(a,log):
    h0,h1=map(type(a[0]),['.0025','.0005']) if isinstance(a[0],mp.mpf) else (.0025,.0005)
    c0=h0/2;c1=h0+h1/2;total=h0+h1;knots=[0,c0,h0,c1,total]
    def area(z):
        if z.real<=c0:return a[0]
        if z.real>=c1:return a[1]
        return a[0]*(c1-z)/(c1-c0)+a[1]*(z-c0)/(c1-c0)
    volumes=[];centres=[]
    for lo,hi in [(0,h0),(h0,total)]:
        points=[lo]+[z for z in knots if lo<z<hi]+[hi];v=0;moment=0
        for x,y in zip(points,points[1:]):
            ax,ay=area(x),area(y);length=y-x
            v+=length*(ax+ay)/2;moment+=length*((2*x+y)*ax+(x+2*y)*ay)/6
        volumes.append(v);centres.append(moment/v)
    i1=[[0,0] for _ in range(3)];i3=[[0,0] for _ in range(3)]
    for p,(lo,hi) in enumerate(zip([0,*centres],[*centres,total])):
        points=[lo]+[z for z in knots if lo.real<z<hi.real]+[hi]
        for x,y in zip(points,points[1:]):
            ax,ay=area(x),area(y);d=y-x;cell=0 if (x.real+y.real)/2<h0 else 1
            r=(ay-ax)/ax
            if abs(r)<.05:
                series=0;power=1
                for k in range(30):series+=power/(k+1);power*=-r
                integral=d/ax*series
            else:integral=d*log(ay/ax)/(ay-ax)
            i1[p][cell]+=integral;i3[p][cell]+=d*(ax+ay)/(2*ax*ax*ay*ay)
    return volumes,i1,i3

def metrics(x,log=np.log):
    v,i1,i3=geometry(x[2:],log);rho=mp.mpf('1.14') if isinstance(x[0],mp.mpf) else 1.14
    inv=[1/(rho*sum(v[i]/x[i]*i1[j][i] for i in range(2))) for j in range(3)]
    return v,inv,i3

def audit(out):
    binary=ROOT/'build/ducklips_contour_kinetics_probe';rng=np.random.default_rng(741);rows=[];pairs=[]
    for n in range(180):
        area=10**rng.uniform(-8,-4,2)
        if n<20:area=np.array([3e-6,3e-6*(1+(n-10)*1e-10)])
        v,_,_=geometry(area,np.log);mass=np.array(v)*rng.uniform(.6,2,2);p=rng.uniform(-.005,.005,3)
        rows.append([*mass,*area,*p])
        if n<60:
            scale=1+rng.uniform(-1e-10,1e-10,4) if n<30 else 10**rng.uniform(-1,1,4)
            other=np.r_[mass,area]*scale;p1=rng.uniform(-.005,.005,3)
            pairs.append([*mass,*area,*p,*other,*p1])
    maxima=dict(energy_scaled=0.,metric_relative=0.,effort_scaled=0.,transfer_scaled=0.,flow_scaled=0.,resistance_relative=0.)
    mp.mp.dps=80
    def e(x,p):return sum(pi*pi*gi/2 for pi,gi in zip(p,metrics(x,mp.log)[1]))
    for tag,flag,data in [('kinetics','--batch',rows),('transfers','--pairs',pairs)]:
        inp=out/(tag+'-input.txt');inp.write_text(''.join(' '.join(map(str,r))+'\n' for r in data))
        run=subprocess.run([str(binary),flag],input=inp.read_text(),text=True,capture_output=True);assert run.returncode==0,run.stderr
        (out/(tag+'-output.jsonl')).write_text(run.stdout);results=[json.loads(s) for s in run.stdout.splitlines()];assert len(results)==len(data)
        for row,result in zip(data,results):
            row=list(map(mp.mpf,row));x=row[:4];p=row[4:7];energy=e(x,p)
            if tag=='kinetics':
                vol,inv,i3=metrics(x,mp.log)
                maxima['energy_scaled']=max(maxima['energy_scaled'],float(abs(result['energy']-energy)/energy))
                for j in range(3):
                    maxima['metric_relative']=max(maxima['metric_relative'],float(abs(result['inverse'][j]/inv[j]-1)))
                    resistance=12*mp.mpf('1.86e-5')*mp.mpf('.014')**2*sum((vol[i]/x[i])**2*i3[j][i] for i in range(2))
                    maxima['resistance_relative']=max(maxima['resistance_relative'],float(abs(result['resistance'][j]/resistance-1)))
                for k in range(4):
                    derivative=mp.diff(lambda z:e([z if j==k else x[j] for j in range(4)],p),x[k])
                    maxima['effort_scaled']=max(maxima['effort_scaled'],float(abs(result['effort'][k]-derivative)*x[k]/energy))
            else:
                y=row[7:11];p1=row[11:];scale=energy+e(y,p1);expected=[mp.mpf(0)]*4
                for order in [range(4),range(3,-1,-1)]:
                    current=x.copy()
                    for k in order:
                        nxt=current.copy();nxt[k]=y[k]
                        if nxt[k]==current[k]:
                            d=sum(mp.diff(lambda z:e([z if j==k else current[j] for j in range(4)],q),current[k]) for q in [p,p1])/4
                        else:d=sum(e(nxt,q)-e(current,q) for q in [p,p1])/(4*(nxt[k]-current[k]))
                        expected[k]+=d;current=nxt
                for k in range(4):maxima['transfer_scaled']=max(maxima['transfer_scaled'],float(abs(result['effort'][k]-expected[k])*max(x[k],y[k])/scale))
                for j,(g0,g1) in enumerate(zip(metrics(x,mp.log)[1],metrics(y,mp.log)[1])):
                    q=(g0+g1)*(p[j]+p1[j])/4;maxima['flow_scaled']=max(maxima['flow_scaled'],float(abs(result['flow'][j]-q)/((g0+g1)*(abs(p[j])+abs(p1[j]))/4)))
    assert max(maxima.values())<1e-9,maxima
    return dict(instant_cases=len(rows),finite_transfer_cases=len(pairs),digits=80,maxima=maxima)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--steps',type=int,nargs='+',default=[120,240,480,960,1920,3840]);parser.add_argument('--modes',nargs='+',default=['forward','reverse','relaxation','rest_motion']);args=parser.parse_args()
    out=ROOT/'out/contour-glottal';out.mkdir(parents=True,exist_ok=True);kernel=audit(out);print('kernel',kernel,flush=True)
    binary=ROOT/'build/ducklips_contour_glottal_probe';bulk=1.14*350**2;gamma=1.4;duration=.005;width=.014
    tissue_mass=np.array([1.25e-4]*2+[2.5e-5]*2);stiff=np.array([80.]*2+[8.]*2);damping=2*np.array([.1]*2+[.6]*2)*np.sqrt(tissue_mass*stiff)
    c=np.array([2e-10,3e-10,1e-10]);g=np.array([.002,.003,.004,.002]);b=np.array([[1,0,1,1],[-1,1,0,0],[0,-1,-1,0]])
    # Open geometry volume is linear in the controls. Get its independent
    # gradients by integrating each basis surface, not native force areas.
    volume_gradient=np.array([geometry(np.array([1.,0.])+1e-30,np.log)[0],geometry(np.array([0.,1.])+1e-30,np.log)[0]]).T
    runs=[];oracles=[];acceptance=[]
    for mode in args.modes:
        reservoir=0. if mode=='relaxation' else -200. if mode=='reverse' else 700.
        reservoir_mu=bulk/(gamma-1)*np.expm1((gamma-1)/gamma*np.log1p(gamma*reservoir/bulk))
        pressure=np.array([700.,300.]) if mode=='relaxation' else np.zeros(2)
        m0=np.array([.0025,.0005])*3e-6*np.power(1+gamma*pressure/bulk,1/gamma);y0=np.r_[m0,np.zeros(18)]
        if mode=='relaxation':y0[10:13]=[.001,-.0005,.0004];y0[13:16]=np.array([2e-8,-1e-8,3e-8])/c
        def derivative(t,y):
            m=y[:2];x=y[2:6];velocity=y[6:10];momentum=y[10:13];pressure=y[13:16];flow=y[16:]
            rest=3e-6-5e-7*t/duration if mode=='rest_motion' else 3e-6
            area=rest+width*np.array([x[0]+x[1],x[2]+x[3]]);assert min(area)>0
            coordinates=np.r_[m,area];vol,inv,i3=metrics(coordinates);vol=np.array(vol);inv=np.array(inv);q=inv*momentum
            # Complex-step differentiation of the independent scalar kinetic
            # energy. No native secants or solver Jacobian are reused.
            effort=np.empty(4)
            for k in range(4):
                z=coordinates.astype(complex);z[k]+=1e-25j
                effort[k]=sum(momentum[j]**2*gj/2 for j,gj in enumerate(metrics(z)[1])).imag/1e-25
            log_r=np.log(m/vol);p=bulk/gamma*np.expm1(gamma*log_r);mu=bulk/(gamma-1)*np.expm1((gamma-1)*log_r)
            area_effort=effort[2:]-p@volume_gradient
            force=-width*np.repeat(area_effort,2)-stiff*x-damping*velocity;coupled=25*(x[:2]-x[2:]);force[:2]-=coupled;force[2:]+=coupled
            resistance=12*1.86e-5*width**2*(np.array(i3)@((vol/m)**2));quadratic=np.array([0.,0.,.5*1.14*(vol[1]/(m[1]*area[1]))**2]);mu+=effort[:2]
            pdot=np.array([reservoir_mu-mu[0],mu[0]-mu[1],mu[1]-pressure[2]])-resistance*q-quadratic*q*abs(q)
            drive=-b@flow;drive[2]+=q[2]
            return np.r_[q[0]-q[1],q[1]-q[2],velocity,force/tissue_mass,pdot,drive/c,g*(b.T@pressure)]
        options=dict(method='DOP853',rtol=2e-11,atol=np.r_[[1e-23]*2,[1e-16]*4,[1e-10]*4,[1e-12]*3,[1e-8]*3,[1e-15]*4])
        ode=solve_ivp(derivative,(0,duration),y0,max_step=2e-6,**options);assert ode.success
        refined=solve_ivp(derivative,(0,duration),y0,max_step=1e-6,**options);assert refined.success
        scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3,[700.]*3,[1e-4]*4];difference=abs(ode.y[:,-1]-refined.y[:,-1]);assert max(difference/scale)<1e-8,difference
        oracles.append(dict(mode=mode,final_state=refined.y[:,-1].tolist(),normalized_refinement_difference=float(max(difference/scale))))
        previous=None
        for steps in args.steps:
            command=[str(binary),'--trajectory',str(steps),str(duration),mode];run=subprocess.run(command,capture_output=True,text=True);assert run.returncode==0,run.stderr
            path=out/f'{mode}-{steps}.json';path.write_text(run.stdout);native=json.loads(run.stdout)
            y=np.r_[native['mass'],native['position'],native['velocity'],native['passage_momentum'],np.array(native['charge'])/c,g*np.array(native['network_momentum'])]
            error=abs(y-refined.y[:,-1]);norm=float(max(error/scale))
            runs.append(dict(mode=mode,steps=steps,normalized_error=norm,refinement_ratio=None if previous is None else previous/norm,
                maximum_network_pressure_error_Pa=float(max(error[13:16])),maximum_position_error_m=float(max(error[2:6])),
                maximum_work_balance_J=native['maximum_work_balance'],maximum_mass_balance_m3=native['maximum_mass_balance'],
                command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path)));previous=norm
        print(json.dumps(runs[-1]),flush=True);acceptance.append(dict(mode=mode,second_order=all(3.5<r['refinement_ratio']<4.5 for r in runs[-2:]),endpoint_accuracy=norm<1e-4))
    sources=['include/ducklips/contour_kinetics.hpp','src/contour_kinetics.cpp','include/ducklips/contour_glottal_storage.hpp','src/contour_glottal_storage.cpp','include/ducklips/glottal_contour.hpp','src/glottal_contour.cpp','src/gas_storage.cpp','src/glottal_storage.cpp','src/moving_graph.cpp','tools/contour_kinetics_probe.cpp','tools/contour_glottal_probe.cpp','tools/contour_glottal_reference.py']
    report=dict(scope='open contour coupled to gas/tissue/passage momentum and linear tract; no contact transport or speech claim',duration_s=duration,kernel=kernel,runs=runs,oracles=oracles,acceptance=acceptance,hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),binary,ROOT/'build/ducklips_contour_kinetics_probe']},versions=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__))
    (ROOT/'reports/contour_glottal_reference.json').write_text(json.dumps(report,indent=2)+'\n');assert all(x['second_order'] and x['endpoint_accuracy'] for x in acceptance),acceptance
if __name__=='__main__':main()
