#!/usr/bin/env python3
"""Independent continuous continuation of a sealed lower cap near its chamber seam."""
import argparse
from functools import partial
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[1]
WIDTH=.014; BULK=1.14*350**2; GAMMA=1.4
C=np.array([2e-10,3e-10,1e-10]); G=np.array([.002,.003,.004,.002])
B=np.array([[1,0,1,1],[-1,1,0,0],[0,-1,-1,0]])
MASS=np.array([1.25e-4,2.5e-5]); SPRING=np.array([80.,8.]); DAMPING=2*np.array([.1,.6])*np.sqrt(MASS*SPRING)
SCALE=np.r_[1.,1.,1e-6,1.,.01,.01,.001,.001,[1000.]*3,[1e-4]*4,1e-6]

def metric(c):
    """Coordinates M0,M1,lower opening, seam opening; lower contact, upper open.

    Integrals use local areas, avoiding subtraction of nearly coincident axial
    coordinates as the contact point approaches the fixed chamber seam.
    """
    m0,m1,a,s=c;b=(6*s-a)/5;slope=(b-a)/.0015;d=.00025
    v0=s*s/(2*slope);v1=d*(s+3*b)/2
    upper_offset=d*(s+11*b)/(3*(s+3*b));plateau=upper_offset-d;outlet=2*d-upper_offset
    i1=np.array([[np.log(1.5)/slope,np.log(b/s)/slope+plateau/b],[0.,outlet/b]])
    i3=np.array([[5/(8*slope*s*s),(1/(s*s)-1/(b*b))/(2*slope)+plateau/b**3],[0.,outlet/b**3]])
    volumes=np.array([v0,v1]);mass=np.array([m0,m1]);g=1/(1.14*(i1@(volumes/mass)))
    return volumes,g,i3,b

def decode(z):
    z=z*SCALE;m=np.exp(z[:2]);a=z[2];s=np.exp(z[3]);b=(6*s-a)/5
    return m,a,s,b,z[4:6],z[6:8],z[8:11],z[11:15]

def density_factor(r0,r1):
    # Integral of r dp gives Phi = B/(gamma+1) r**(gamma+1).
    # R = 12 mu width**2 I3 * (mu0-mu1)/(Phi0-Phi1).
    hi=max(r0,r1);lo=min(r0,r1);t=np.log(lo/hi)
    if t==0:return 1/(hi*hi)
    return (GAMMA+1)/(GAMMA-1)/(hi*hi)*np.expm1((GAMMA-1)*t)/np.expm1((GAMMA+1)*t)

def derivative(t,z,pressure_loss=False):
    m,a,s,b,rates,p,pressure,flow=decode(z);coords=np.r_[m,a,s];vol,g,i3,_=metric(coords);q=g*p
    eff=np.empty(4);dv=np.empty((2,2))
    for k in range(4):
        h=abs(coords[k])*1e-20;c=coords.astype(complex);c[k]+=h*1j;vz,gz,_,_=metric(c)
        eff[k]=np.imag(np.sum(p*p*gz/2))/h
        if k>=2:dv[:,k-2]=np.imag(vz)/h
    logr=np.log(m/vol);pg=BULK/GAMMA*np.expm1(GAMMA*logr);mu=BULK/(GAMMA-1)*np.expm1((GAMMA-1)*logr)+eff[:2]
    shape=eff[2:]-pg@dv
    # Convert (a,seam) shape derivatives to the original (a,b) controls.
    ae=np.array([shape[0]+shape[1]/6,5*shape[1]/6])+3*SPRING/2*np.minimum([a,b],0)/WIDTH**2
    x=(np.array([a,b])-3e-6)/(2*WIDTH);v=np.array([rates[0],(6*rates[1]-rates[0])/5])/(2*WIDTH)
    force=-WIDTH*ae-SPRING*x-DAMPING*v;coupled=25*(x[0]-x[1]);force+=np.array([-coupled,coupled])
    acc=2*WIDTH*force/MASS;rate_acc=np.array([acc[0],(acc[0]+5*acc[1])/6])
    resistance=12*1.86e-5*WIDTH**2*(i3@((vol/m)**2));quad=.5*1.14*(vol[1]/(m[1]*b))**2
    if pressure_loss:
        r=m/vol;outlet_density=(1+GAMMA*pressure[2]/BULK)**(1/GAMMA)
        resistance=12*1.86e-5*WIDTH**2*np.sum(i3,axis=1)*np.array([density_factor(r[0],r[1]),density_factor(r[1],outlet_density)])
    pdot=np.array([mu[0]-mu[1],mu[1]-pressure[2]])-resistance*q;pdot[1]-=quad*q[1]*abs(q[1])
    drive=-B@flow;drive[2]+=q[1]
    dissipation=np.sum(resistance*q*q)+quad*abs(q[1])**3+2*np.sum(DAMPING*v*v)
    return np.r_[-q[0]/m[0],(q[0]-q[1])/m[1],rates[0],rates[1]/s,rate_acc,pdot,drive/C,G*(B.T@pressure),dissipation]/SCALE

def energy(z):
    m,a,s,b,rates,p,pressure,flow=decode(z);vol,g,_,_=metric(np.r_[m,a,s]);r=m/vol
    gas=np.sum(vol*BULK/(GAMMA*(GAMMA-1))*(np.expm1(GAMMA*np.log(r))-GAMMA*(r-1)))
    x=(np.array([a,b])-3e-6)/(2*WIDTH);v=np.array([rates[0],(6*rates[1]-rates[0])/5])/(2*WIDTH)
    tissue=np.sum(MASS*v*v+SPRING*x*x)+25*(x[0]-x[1])**2
    contact=np.sum(.5*(3*SPRING/2)*(np.minimum([a,b],0)/WIDTH)**2)
    network=.5*np.sum(C*pressure*pressure)+.5*np.sum(flow*flow/G)
    return gas+.5*np.sum(g*p*p)+tissue+contact+network

def rebound(t,z):return z[5]
rebound.direction=1

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='out/contour-pocket/native-76800.jsonl');ap.add_argument('--pressure-loss',action='store_true');ap.add_argument('--stop-seam',type=float,default=1e-13);ap.add_argument('--rtol',type=float,default=1e-8);ap.add_argument('--duration',type=float,default=4e-5);ap.add_argument('--out',default='out/contour-pocket/continuous.json');args=ap.parse_args()
    path=ROOT/args.input;checkpoint=next(json.loads(s) for s in path.read_text().splitlines() if json.loads(s).get('tag')=='checkpoint');y=np.array(checkpoint['state']);a,b=checkpoint['area'];s=(a+5*b)/6
    assert np.array_equal(y[2:10:2],y[3:10:2]),'symmetric body required'
    rates=WIDTH*np.array([y[6]+y[7],y[8]+y[9]]);z0=np.r_[np.log(y[:2]),a,np.log(s),rates[0],(rates[0]+5*rates[1])/6,y[11:13],y[13:],0.]/SCALE
    def boundary(t,z):return z[3]-np.log(args.stop_seam)
    boundary.terminal=True;boundary.direction=-1
    sol=solve_ivp(partial(derivative,pressure_loss=args.pressure_loss),(0,args.duration),z0,method='Radau',rtol=args.rtol,atol=args.rtol*.01,max_step=1e-7,events=[rebound,boundary] if args.pressure_loss else rebound)
    events=[]
    samples=[("start",0.,sol.y[:,0])]+[("rebound",t,z) for t,z in zip(sol.t_events[0],sol.y_events[0])]+[("final",sol.t[-1],sol.y[:,-1])]
    for tag,t,z in samples:
        m,a,s,b,rates,p,pressure,flow=decode(z);vol,g,_,_=metric(np.r_[m,a,s]);logr=np.log(m/vol)
        events.append(dict(tag=tag,time=float(checkpoint['time']+t),area=[a,b],seam=s,volume=vol.tolist(),mass=m.tolist(),relative_density=(m/vol).tolist(),pressure=(BULK/GAMMA*np.expm1(GAMMA*logr)).tolist(),area_rate=[rates[0],(6*rates[1]-rates[0])/5],momentum=p.tolist(),energy_J=energy(z),dissipation_J=float(z[-1]*SCALE[-1]),work_residual_J=energy(z)+float(z[-1]*SCALE[-1])-energy(z0)))
    result=dict(loss_law="continuous_pressure" if args.pressure_loss else "piecewise_cell_density",stop_seam=args.stop_seam if args.pressure_loss else None,success=sol.success,message=sol.message,checkpoint=checkpoint,rtol=args.rtol,duration=args.duration,steps=len(sol.t)-1,nfev=sol.nfev,events=events,final_scaled_state=sol.y[:,-1].tolist(),hashes={args.input:hashlib.sha256(path.read_bytes()).hexdigest(),'tools/contour_pocket_reference.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    (ROOT/args.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert sol.success
if __name__=='__main__':main()
