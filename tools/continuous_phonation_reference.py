#!/usr/bin/env python3
"""Independent open full-span continuous equations and pressure-onset audit."""
import json,hashlib
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import root
from contour_pressure_loss_reference import spline,low,high
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-phonation'
H=.003;WIDTH=.014;RHO=1.14;BULK=139650.;GAMMA=1.4;ETA=1.86e-5;REST=3e-6
MASS=np.array([1.25e-4,2.5e-5])/(2*WIDTH**2);STIFF=np.array([[105.,-25.],[-25.,33.]])/(2*WIDTH**2);DAMP=2*np.array([.1,.6])*np.sqrt(np.array([1.25e-4,2.5e-5])*[80.,8.])/(2*WIDTH**2)
YSCALE=np.array([1000.,1000.,1e-6,1e-6,1e-4,1e-4,1e-4,.001,.001]);VSCALE=YSCALE[4:]
class Reference:
 def __init__(self,order=32,derivative=1e-20):
  x,w=leggauss(order);self.t=(x+1)/2;self.weight=w/2;self.derivative=derivative
 def fields(self,t,area,effort):
  phi=np.stack([1-t,t],axis=-1);a=phi@area;u=phi@effort;lr=np.log1p((GAMMA-1)*u/BULK)/(GAMMA-1);r=np.exp(lr);p=BULK/GAMMA*np.expm1(GAMMA*lr);return phi,a,r,p
 def configuration(self,area,effort,metric_only=False):
  assert np.min(area)>0
  t=self.t;w=H*self.weight;phi,a,r,p=self.fields(t,area,effort);C=np.einsum('n,ni,nj->ij',w*a*r**(2-GAMMA)/BULK,phi,phi);D=np.einsum('n,ni,nj->ij',w*r,phi,phi);E=np.linalg.solve(C,np.c_[np.array([[1.,-1.,0.],[0.,1.,-1.]]),-D])
  # Independently integrate continuity from the inlet at each quadrature node.
  tt=t[:,None]*t[None,:];pp,aa,rr,_=self.fields(tt,area,effort);cw=H*t[:,None]*self.weight[None,:];c=np.einsum('nk,nki->ni',cw*aa*rr**(2-GAMMA)/BULK,pp);d=np.einsum('nk,nki->ni',cw*rr,pp);F=-c@E;F[:,0]+=1;F[:,3:]-=d
  L=np.einsum('n,ni,nj->ij',w*RHO/(r*a),F,F);G=L.copy();G[3:,3:]+=np.diag(MASS)
  if metric_only:return G
  thermal=350*np.sqrt(2/GAMMA)*r**.2;delta=a/WIDTH*RHO*r*thermal/(2*ETA);conductance=np.exp(spline(np.log(delta)));conductance=np.where(delta>=100,high(delta),conductance);conductance=np.where(delta<=1e-8,low(delta),conductance)
  R=np.einsum('n,ni,nj->ij',w*RHO*thermal*WIDTH/(r*a*a*conductance),F,F);R[3:,3:]+=np.diag(DAMP)
  return dict(mass=phi.T@(w*a*r),C=C,D=D,area_effort=-phi.T@(w*p),G=G,R=R,L=L,phi=phi,a=a,r=r,w=w)
 def stationary(self,x,pressure):
  u=x[:2]*1000;area=x[2:4]*1e-6;Q=x[4]*1e-4;s=self.configuration(area,u);phi=s['phi'];a=s['a'];r=s['r'];w=s['w'];ku=-.5*RHO*Q*Q/BULK*(phi.T@(w*r**(-GAMMA)/a));ka=-.5*RHO*Q*Q*(phi.T@(w/(r*a*a)));km=np.linalg.solve(s['C'].T,ku);eta=u-km;ga=s['area_effort']-ka+s['D'].T@km+STIFF@(area-REST);v=np.r_[np.full(3,Q),0.,0.];force=np.r_[eta[0],eta[1]-eta[0],-eta[1],ga]+s['R']@v;force[0]-=BULK/(GAMMA-1)*np.expm1((GAMMA-1)/GAMMA*np.log1p(GAMMA*pressure/BULK));rout=(1+(GAMMA-1)*u[1]/BULK)**(1/(GAMMA-1));force[2]+=RHO/(2*rout**2*area[1]**2)*Q*abs(Q)
  return force/np.array([1000.,1000.,1000.,1.,1.])
 def rhs(self,y,pressure):
  physical=y*YSCALE;u=physical[:2];area=physical[2:4];v=physical[4:];s=self.configuration(area,u);udot=np.linalg.solve(s['C'],np.array([v[0]-v[1],v[1]-v[2]])-s['D']@v[3:]);derivatives=[]
  z=physical[:4]
  for j in range(4):
   h=self.derivative*YSCALE[j];perturbed=z.astype(complex);perturbed[j]+=1j*h;derivatives.append(self.configuration(perturbed[2:],perturbed[:2],metric_only=True).imag/h)
  kinetic=np.array([.5*v@d@v for d in derivatives]);km=np.linalg.solve(s['C'].T,kinetic[:2]);eta=u-km;ga=s['area_effort']-kinetic[2:]+s['D'].T@km+STIFF@(area-REST);force=-np.r_[eta[0],eta[1]-eta[0],-eta[1],ga]-s['R']@v
  force[0]+=BULK/(GAMMA-1)*np.expm1((GAMMA-1)/GAMMA*np.log1p(GAMMA*pressure/BULK));rout=(1+(GAMMA-1)*u[1]/BULK)**(1/(GAMMA-1));force[2]-=RHO/(2*rout**2*area[1]**2)*v[2]*abs(v[2]);gdot=sum(d*rate for d,rate in zip(derivatives,np.r_[udot,v[3:]]));accel=VSCALE*np.linalg.solve(s['G']*VSCALE[:,None]*VSCALE[None,:],VSCALE*(force-gdot@v));return np.r_[udot,v[3:],accel]/YSCALE
 def equilibrium(self,pressure,guess):
  sol=root(lambda x:self.stationary(x,pressure),guess,tol=1e-9);error=float(np.max(np.abs(self.stationary(sol.x,pressure))));assert error<1e-8,(pressure,error,sol.message);x=sol.x;y=np.r_[x[:4],np.full(3,x[4]),0.,0.];return x,y,error
 def jacobian(self,y,pressure,h=1e-3):
  columns=[]
  for j in range(9):
   plus=y.copy();minus=y.copy();plus[j]+=h;minus[j]-=h;columns.append((self.rhs(plus,pressure)-self.rhs(minus,pressure))/(2*h))
  return np.array(columns).T

def main():
 OUT.mkdir(parents=True,exist_ok=True);ref=Reference();guess=np.array([.1,0.,3.7,4.2,1.2]);records=[]
 for pressure in [0.,350.,700.,1050.,1400.]:
  x,y,error=ref.equilibrium(pressure,guess);guess=x;eig=np.linalg.eigvals(ref.jacobian(y,pressure));modes=sorted([dict(growth_s=float(z.real),frequency_Hz=float(abs(z.imag)/(2*np.pi))) for z in eig if z.imag>=0],key=lambda r:r['growth_s'],reverse=True);r=dict(pressure_Pa=pressure,equilibrium_scaled=y.tolist(),stationary_residual=error,rhs_residual=float(np.max(np.abs(ref.rhs(y,pressure)))),modes=modes);records.append(r);print(json.dumps(r),flush=True)
 (OUT/'onset-reference.json').write_text(json.dumps(dict(scope='independent smooth open full-span ODE, symmetric tissue subspace; not anatomy validation',cases=records),indent=2)+'\n')
if __name__=='__main__':main()
