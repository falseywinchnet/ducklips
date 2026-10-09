#!/usr/bin/env python3
"""Independent nonlinear Hamiltonian force check, beyond the linear limit."""
import json,subprocess,hashlib
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import root
from contour_pressure_flow_reference import Reference,C0,C1,H,BULK,GAMMA
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'out/continuous-glottis';WIDTH=.014
X,W=leggauss(32)
def gas(area,effort):
 m=np.zeros(2);C=np.zeros((2,2));force=np.zeros(2)
 for lo,hi in zip([0.,C0,C1],[C0,C1,H]):
  z=(lo+hi)/2+(hi-lo)/2*X;t=np.clip((z-C0)/(C1-C0),0,1);phi=np.array([1-t,t]);a=np.array(area)@phi;mu=np.array(effort)@phi;r=np.exp(np.log1p((GAMMA-1)*mu/BULK)/(GAMMA-1));p=BULK/GAMMA*np.expm1(GAMMA*np.log1p((GAMMA-1)*mu/BULK)/(GAMMA-1));weight=W*(hi-lo)/2
  m+=(phi*(a*r*weight)).sum(axis=1);C+=(phi*(a*r**(2-GAMMA)/BULK*weight))@phi.T;force-=(phi*(p*weight)).sum(axis=1)
 return m,C,force
T=np.zeros((5,7));T[:3,:3]=np.eye(3);T[3,3:5]=WIDTH;T[4,5:]=WIDTH
mass=np.array([1.25e-4,1.25e-4/1.17,2.5e-5,2.5e-5/1.17]);stiff=np.array([80.,80*1.17,8.,8*1.17]);coupled=[25.,25*1.17]
area=np.array([2e-6,4e-6]);effort=np.array([700.,-200.]);moments,C,gas_force=gas(area,effort);v=np.array([8e-5,5e-5,3e-5,.02,-.01,.03,-.01]);x=np.array([-3e-5,-1e-6/WIDTH+3e-5,4e-5,1e-6/WIDTH-4e-5]);z=np.r_[moments,area]
def metric(z):
 m,a=z[:2],z[2:];sol=root(lambda u:(gas(a,u)[0]-m)/1e-9,effort,jac=lambda u:gas(a,u)[1]/1e-9,tol=1e-10);assert max(abs(gas(a,sol.x)[0]-m))<1e-20,sol.message
 ref=Reference([*a,*sol.x]);G=T.T@ref.matrix(1e-11)@T;G[3:,3:]+=np.diag(mass);return G

def main():
 derivatives=[]
 for relative in [2e-5,1e-5]:
  kinetic=[]
  for j in range(4):
   h=relative*max(abs(z[j]),1e-9 if j<2 else 1e-6);plus=z.copy();minus=z.copy();plus[j]+=h;minus[j]-=h;derivative=(metric(plus)-metric(minus))/(2*h);kinetic.append(float(-.5*v@derivative@v))
  g=np.r_[effort,gas_force]+kinetic;general=np.r_[g[0],g[1]-g[0],-g[1],np.repeat(g[2:],2)*WIDTH+stiff*x]
  for side,k in enumerate(coupled):general[3+side]+=k*(x[side]-x[side+2]);general[5+side]+=k*(x[side+2]-x[side])
  derivatives.append(dict(relative_step=relative,kinetic_gradient=kinetic,total_momentum_rate=(-general).tolist()))
 scale=np.array([1000.]*3+[.01]*4);target=np.array(derivatives[-1]['total_momentum_rate']);refinement=float(max(abs(target-np.array(derivatives[0]['total_momentum_rate']))/scale));assert refinement<1e-6
 G=metric(z);p0=G@v;records=[];binary=ROOT/'build/ducklips_continuous_glottis_tangent_probe'
 for dt in [1e-9,5e-10,2.5e-10]:
  p=subprocess.run([str(binary),str(dt)],capture_output=True,text=True);assert p.returncode==0,p.stderr;r=json.loads(p.stdout);(OUT/f'tangent-{dt}.json').write_text(p.stdout);rate=(np.array(r['final_momentum'])-np.array(r['initial_momentum']))/dt;error=float(max(abs(rate-target)/scale));assert max(abs(np.array(r['initial_momentum'])-p0)/np.maximum(abs(p0),1e-12))<1e-9
  records.append(dict(**r,momentum_rate=rate.tolist(),maximum_scaled_rate_error=error));print(records[-1],flush=True)
 extrapolated=2*np.array(records[-1]['momentum_rate'])-np.array(records[-2]['momentum_rate']);extrapolation_error=float(max(abs(extrapolated-target)/scale));assert extrapolation_error<1e-6
 report=dict(scope='independent nonlinear canonical momentum derivative with substantial flow and moving tissue',reference='separate Gauss quadrature and moment recovery; differentiate fluid metric at fixed masses/areas and use -v^T dG v/2',derivatives=derivatives,reference_refinement_scaled_error=refinement,steps=records,richardson_scaled_error=extrapolation_error,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'tools/continuous_glottis_tangent_probe.cpp',ROOT/'tools/contour_pressure_flow_reference.py',binary,ROOT/'src/continuous_glottis.cpp']})
 (ROOT/'reports/continuous_glottis_tangent_reference.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS independent nonlinear momentum forces')
if __name__=='__main__':main()
