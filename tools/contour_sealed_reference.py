#!/usr/bin/env python3
"""Independent signed-contour constitutive and closed-passage motion checks."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'out/contour-sealed'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def geometry(area,log):
    precise=isinstance(area[0],mp.mpf);number=mp.mpf if precise else float
    h0,h1=map(number,['.0025','.0005']);c0=h0/2;c1=h0+h1/2;total=h0+h1;knots=[0,c0,h0,c1,total]
    def opening(z):
        if z.real<=c0:return area[0]
        if z.real>=c1:return area[1]
        return area[0]*(c1-z)/(c1-c0)+area[1]*(z-c0)/(c1-c0)
    volumes=[];centres=[]
    for start,end in [(0,h0),(h0,total)]:
        points=[start]+[z for z in knots if start<z<end]+[end];v=0;moment=0
        for lo,hi in zip(points,points[1:]):
            a,b=opening(lo),opening(hi)
            if a.real<=0 and b.real<=0:continue
            if a.real<=0:lo=lo+(hi-lo)*(-a)/(b-a);a=0
            elif b.real<=0:hi=lo+(hi-lo)*a/(a-b);b=0
            d=hi-lo;v+=d*(a+b)/2;moment+=d*((2*lo+hi)*a+(lo+2*hi)*b)/6
        assert v.real>0
        volumes.append(v);centres.append(moment/v)
    i1=[[0,0] for _ in range(3)];i3=[[0,0] for _ in range(3)];open_path=[area[0].real>0,True,area[1].real>0]
    for j,(first,last) in enumerate(zip([0,*centres],[*centres,total])):
        if not open_path[j]:continue
        points=[first]+[z for z in knots if first.real<z<last.real]+[last]
        for lo,hi in zip(points,points[1:]):
            a,b=opening(lo),opening(hi);assert a.real>0 and b.real>0
            d=hi-lo;cell=0 if (lo.real+hi.real)/2<h0 else 1;r=(b-a)/a
            if abs(r)<.05:
                series=0;power=1
                for k in range(30):series+=power/(k+1);power*=-r
                inverse=d/a*series
            else:inverse=d*log(b/a)/(b-a)
            i1[j][cell]+=inverse;i3[j][cell]+=d*(a+b)/(2*a*a*b*b)
    return volumes,i1,i3,open_path

def metrics(x,log=np.log):
    v,i1,i3,opened=geometry(x[2:],log);rho=mp.mpf('1.14') if isinstance(x[0],mp.mpf) else 1.14
    g=[1/(rho*sum(v[i]*i1[j][i]/x[i] for i in range(2))) if opened[j] else 0 for j in range(3)]
    return v,g,i3

def kernel(binary):
    rng=np.random.default_rng(9041);rows=[];pairs=[]
    for n in range(160):
        positive=10**rng.uniform(-7,-4);area=np.array([-positive*rng.uniform(.01,3),positive]) if n%2==0 else np.array([positive,-positive*rng.uniform(.01,.07)])
        if n<4:area[n%2]=0.
        volume=metrics(np.r_[1e-9,1e-9,area])[0];mass=np.array(volume)*rng.uniform(.7,1.5,2);p=rng.uniform(-.003,.003,3);p[0 if n%2==0 else 2]=0
        rows.append([*mass,*area,*p])
        if n>=4 and n<84:
            x=np.r_[mass,area];y=x*(1+rng.uniform(-1e-10,1e-10,4)) if n<44 else x*rng.uniform(.8,1.2,4)
            p1=rng.uniform(-.003,.003,3);p1[0 if n%2==0 else 2]=0;pairs.append([*x,*p,*y,*p1])
    # Constitutive endpoint transitions are checked independently of the
    # dynamic solver, which still rejects an unowned topology change.
    for n in range(12):
        area=np.array([3e-6,3e-6]);target=area.copy();target[n%2]=-1e-7
        v=metrics(np.r_[1e-9,1e-9,area])[0];w=metrics(np.r_[1e-9,1e-9,target])[0]
        p=rng.uniform(-.003,.003,3);p1=rng.uniform(-.003,.003,3);p1[0 if n%2==0 else 2]=0;pairs.append([*v,*area,*p,*w,*target,*p1])
    mp.mp.dps=80;maxima=dict(energy_scaled=0.,effort_scaled=0.,metric_relative=0.,resistance_relative=0.,transfer_scaled=0.,volume_gradient_scaled=0.,chain_scaled=0.)
    def energy(x,p):return sum(pi*pi*g/2 for pi,g in zip(p,metrics(x,mp.log)[1]))
    for name,flag,data in [('kernel','--kernel',rows),('pairs','--pairs',pairs)]:
        inp=OUT/(name+'-input.txt');inp.write_text(''.join(' '.join(map(str,r))+'\n' for r in data));run=subprocess.run([str(binary),flag],input=inp.read_text(),capture_output=True,text=True);assert run.returncode==0,run.stderr
        (OUT/(name+'-output.jsonl')).write_text(run.stdout);results=[json.loads(s) for s in run.stdout.splitlines()];assert len(results)==len(data)
        for row,s in zip(data,results):
            assert s['accepted'],row;row=list(map(mp.mpf,row));x=row[:4];p=row[4:7];e=energy(x,p)
            if name=='kernel':
                volume,inv,i3=metrics(x,mp.log);maxima['energy_scaled']=max(maxima['energy_scaled'],float(abs(s['energy']-e)/e))
                for j in range(3):
                    if inv[j]==0:assert s['inverse'][j]==s['resistance'][j]==0
                    else:
                        maxima['metric_relative']=max(maxima['metric_relative'],float(abs(s['inverse'][j]/inv[j]-1)))
                        r=12*mp.mpf('1.86e-5')*mp.mpf('.014')**2*sum((volume[i]/x[i])**2*i3[j][i] for i in range(2));maxima['resistance_relative']=max(maxima['resistance_relative'],float(abs(s['resistance'][j]/r-1)))
                for k in range(4):
                    derivative=mp.diff(lambda z:energy([z if j==k else x[j] for j in range(4)],p),x[k],direction=-1 if k>=2 and x[k]==0 else 0)
                    scale=abs(x[k]) if x[k]!=0 else max(abs(x[2]),abs(x[3]));maxima['effort_scaled']=max(maxima['effort_scaled'],float(abs(s['effort'][k]-derivative)*scale/e))
            else:
                y=row[7:11];p1=row[11:];scale=e+energy(y,p1);expected=[mp.mpf(0)]*4;dv=[[mp.mpf(0)]*2 for _ in range(2)]
                for order in [range(4),range(3,-1,-1)]:
                    current=x.copy()
                    for k in order:
                        nxt=current.copy();nxt[k]=y[k];delta=nxt[k]-current[k]
                        if delta==0:derivative=sum(mp.diff(lambda z:energy([z if j==k else current[j] for j in range(4)],q),current[k]) for q in [p,p1])/4
                        else:derivative=sum(energy(nxt,q)-energy(current,q) for q in [p,p1])/(4*delta)
                        expected[k]+=derivative
                        if k>=2:
                            v=metrics(current,mp.log)[0];w=metrics(nxt,mp.log)[0]
                            for i in range(2):
                                d=mp.diff(lambda z:metrics([z if j==k else current[j] for j in range(4)],mp.log)[0][i],current[k]) if delta==0 else (w[i]-v[i])/delta
                                dv[i][k-2]+=d/2
                        current=nxt
                for k in range(4):maxima['transfer_scaled']=max(maxima['transfer_scaled'],float(abs(s['effort'][k]-expected[k])*max(abs(x[k]),abs(y[k]))/scale))
                for i,label in enumerate(['volume_gradient_lower','volume_gradient_upper']):
                    for k in range(2):maxima['volume_gradient_scaled']=max(maxima['volume_gradient_scaled'],float(abs(s[label][k]-dv[i][k])/mp.mpf('.003')))
                maxima['chain_scaled']=max(maxima['chain_scaled'],float(abs(s['chain_residual'])/scale))
    assert max(maxima.values())<1e-9,maxima
    return dict(instant_cases=len(rows),transfer_cases=len(pairs),digits=80,maxima=maxima)

def motions(binary):
    width=.014;bulk=1.14*350**2;gamma=1.4;asym=1.17;q=np.array([1,asym,1,asym]);mass=np.array([1.25e-4]*2+[2.5e-5]*2)/q;spring=np.array([80.]*2+[8.]*2)*q;damping=2*np.array([.1]*2+[.6]*2)*np.sqrt(mass*spring);coupling=25*np.array([1,asym])
    pair_stiff=3/(1/spring[[0,2]]+1/spring[[1,3]]);runs=[];oracles=[]
    reservoir_mu=bulk/(gamma-1)*np.expm1((gamma-1)/gamma*np.log1p(gamma*700/bulk))
    for mode in ['lower','upper','lower_rest','upper_rest']:
        lower=mode.startswith('lower');closed=0 if lower else 2;x=np.zeros(4)
        if lower:x[:2]=-3.1e-6/(2*width)
        else:x[2:]=-3.1e-6/(2*width)
        area=3e-6+width*np.array([sum(x[:2]),sum(x[2:])]);volume=metrics(np.r_[1e-9,1e-9,area])[0]
        m0=np.array(volume)*np.power(1+gamma*np.array([100.,50.])/bulk,1/gamma);momentum=np.array([0.,.0002,.0001]) if lower else np.array([.0001,.0002,0.]);y0=np.r_[m0,x,np.zeros(4),momentum]
        force=np.array([-.025,-.025,0.,0.]) if lower else np.array([0.,0.,-.008,-.008])
        def derivative(t,y):
            m=y[:2];x=y[2:6];v=y[6:10];p=y[10:];rest=3e-6-(2e-8*t/.0002 if mode.endswith('_rest') else 0.);a=rest+width*np.array([sum(x[:2]),sum(x[2:])]);assert (a[0]<0<a[1]) if lower else (a[1]<0<a[0])
            coords=np.r_[m,a];volume,g,i3=metrics(coords);volume=np.array(volume);flow=np.array(g)*p;effort=np.zeros(4);volume_gradient=np.zeros((2,2))
            for k in range(4):
                z=coords.astype(complex);z[k]+=1e-25j;vz,gz,_=metrics(z);effort[k]=sum(p[j]*p[j]*gj/2 for j,gj in enumerate(gz)).imag/1e-25
                if k>=2:volume_gradient[:,k-2]=np.imag(vz)/1e-25
            log_r=np.log(m/volume);pressure=bulk/gamma*np.expm1(gamma*log_r);mu=bulk/(gamma-1)*np.expm1((gamma-1)*log_r)+effort[:2]
            area_effort=effort[2:]-pressure@volume_gradient+pair_stiff*np.minimum(a,0)/width**2
            acceleration=force-width*np.repeat(area_effort,2)-spring*x-damping*v;coupled=coupling*(x[:2]-x[2:]);acceleration[:2]-=coupled;acceleration[2:]+=coupled
            resistance=12*1.86e-5*width**2*(np.array(i3)@((volume/m)**2));quadratic=np.zeros(3)
            if a[1]>0:quadratic[2]=.5*1.14*(volume[1]/(m[1]*a[1]))**2
            pdot=np.array([reservoir_mu-mu[0],mu[0]-mu[1],mu[1]])-resistance*flow-quadratic*flow*abs(flow);pdot[closed]=0
            return np.r_[flow[0]-flow[1],flow[1]-flow[2],v,acceleration/mass,pdot]
        options=dict(method='DOP853',rtol=2e-11,atol=np.r_[[1e-23]*2,[1e-16]*4,[1e-10]*4,[1e-12]*3]);coarse=solve_ivp(derivative,(0,.0002),y0,max_step=2e-7,**options);fine=solve_ivp(derivative,(0,.0002),y0,max_step=1e-7,**options);assert coarse.success and fine.success
        scale=np.r_[[1e-8]*2,[1e-4]*4,[1.]*4,[.05]*3];refinement=float(max(abs(coarse.y[:,-1]-fine.y[:,-1])/scale));assert refinement<1e-8
        oracles.append(dict(mode=mode,endpoint=fine.y[:,-1].tolist(),normalized_refinement_difference=refinement));previous=None
        for steps in [100,200,400,800,1600,3200,6400]:
            command=[str(binary),'--trajectory',str(steps),mode];run=subprocess.run(command,capture_output=True,text=True);assert run.returncode==0,run.stderr
            path=OUT/f'{mode}-{steps}.json';path.write_text(run.stdout);s=json.loads(run.stdout);state=np.r_[s['mass'],s['position'],s['velocity'],s['momentum']];error=abs(state-fine.y[:,-1]);norm=float(max(error/scale))
            runs.append(dict(mode=mode,steps=steps,normalized_error=norm,refinement_ratio=None if previous is None else previous/norm,maximum_position_error_m=float(max(error[2:6])),maximum_mass_balance_m3=s['maximum_mass_balance'],maximum_work_balance_J=s['maximum_work_balance'],command=command,output=str(path.relative_to(ROOT)),output_sha256=sha(path)));previous=norm
        print(json.dumps(runs[-1]),flush=True)
    return runs,oracles

def main():
    OUT.mkdir(parents=True,exist_ok=True);binary=ROOT/'build/ducklips_contour_sealed_probe';kernels=kernel(binary);print('kernel',kernels,flush=True);runs,oracles=motions(binary)
    acceptance=[dict(mode=mode,second_order=all(3.5<r['refinement_ratio']<4.5 for r in [r for r in runs if r['mode']==mode][-2:]),endpoint_accuracy=[r for r in runs if r['mode']==mode][-1]['normalized_error']<1e-4) for mode in ['lower','upper','lower_rest','upper_rest']]
    sources=['include/ducklips/contour_kinetics.hpp','src/contour_kinetics.cpp','include/ducklips/contour_glottal_storage.hpp','src/contour_glottal_storage.cpp','include/ducklips/fold_contact.hpp','src/fold_contact.cpp','src/glottal_contour.cpp','src/gas_storage.cpp','src/glottal_storage.cpp','src/moving_graph.cpp','tools/contour_sealed_probe.cpp','tools/contour_sealed_reference.py']
    report=dict(scope='positive retained chambers with a fixed sealed cap and elastic tissue contact; no automatic sealing or reopening',kernel=kernels,runs=runs,oracles=oracles,acceptance=acceptance,hashes={str(p.relative_to(ROOT)):sha(p) for p in [*(ROOT/s for s in sources),binary]},versions=dict(numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__))
    (ROOT/'reports/contour_sealed_reference.json').write_text(json.dumps(report,indent=2)+'\n');assert all(x['second_order'] and x['endpoint_accuracy'] for x in acceptance),acceptance
if __name__=='__main__':main()
