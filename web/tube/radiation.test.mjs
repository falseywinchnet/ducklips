import test from 'node:test';
import assert from 'node:assert/strict';
import {ApertureRadiation} from './radiation.mjs';
test('aperture FIR is causal, passive at fixed radius, and retains the end-correction delay',()=>{
 const model=new ApertureRadiation(24000);
 for(const radius of [.000025,.001,.006,.0125,.0165,.025,.06]){
  const c=model.coefficients(radius),sum=c.reduce((s,v)=>s+v,0),mean=c.reduce((s,v,k)=>s+k*v,0),expected=2*.6133*radius/343*24000;
  assert.ok(c.every(v=>Number.isFinite(v)&&v>=0));assert.ok(Math.abs(sum-1)<1e-12);assert.ok(Math.abs(mean-expected)<.001*expected+1e-8);
  for(let j=0;j<=64;j++){const w=Math.PI*j/64;let re=0,im=0;for(let k=0;k<c.length;k++){re+=c[k]*Math.cos(w*k);im+=c[k]*Math.sin(w*k);}assert.ok(Math.hypot(re,im)<=1+1e-12);}
 }
});
test('finite radiation approaches the published causal response as the time lattice becomes finer',()=>{
 const model=new ApertureRadiation(96000),r=.0165,c=model.coefficients(r),alpha=1.2266,shape=4*.6133**2;
 for(const f of [200,500,1000,3000]){
  const w=2*Math.PI*f/model.rate,ka=2*Math.PI*f*r/343;let re=0,im=0;for(let k=0;k<c.length;k++){re+=c[k]*Math.cos(w*k);im-=c[k]*Math.sin(w*k);}
  const expected=(1+(ka/alpha)**2)**(-shape/2),phase=-shape*Math.atan(ka/alpha);
  assert.ok(Math.abs(Math.hypot(re,im)-expected)<.004);assert.ok(Math.abs(Math.atan2(im,re)-phase)<.003);
 }
});
