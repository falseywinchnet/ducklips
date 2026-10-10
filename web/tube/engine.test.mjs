import test from 'node:test';import assert from 'node:assert/strict';
import {initialDesign,validateDesign,compile,TubeEngine,eigen,modeShape,N,P,C} from './engine.mjs';
import {receptance} from './sheet.mjs';
const near=(a,b,t=1e-9)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b} (tolerance ${t})`);
const design=initialDesign(),compiled=compile(design);
test('volume blocks add their actual enclosed volume',()=>{const d=initialDesign();d.blocks.push({x:.42,width:.14,volume:1.5});const b=compile(d);near(b.area.reduce((s,v,i)=>s+(v-compiled.area[i])*b.dx[i],0),1.5e-6,1e-14);});
test('symmetric eigendecomposition returns orthonormal eigenpairs',()=>{const K=[3,-1,0,-1,4,-2,0,-2,5],e=eigen(K,3);for(let m=0;m<3;m++)for(let i=0;i<3;i++){let value=0;for(let j=0;j<3;j++)value+=K[i*3+j]*e.vectors[j*3+m];near(value,e.values[m]*e.vectors[i*3+m],1e-8);}});
test('384 modes are sorted and satisfy the complete periodic sheet equation',()=>{
 const s=compiled.sheet;assert.equal(s.modes.length,N*C);assert.ok(s.modes.length>=40);assert.equal(s.channels[0].omega2.length,48);
 for(let m=1;m<s.modes.length;m++)assert.ok(s.modes[m].frequency>=s.modes[m-1].frequency);
 for(const mode of s.modes){const ch=s.channels[mode.channel],shape=modeShape(s,s.modes.indexOf(mode));let residual=0,scale=0;
  for(let t=0;t<C;t++)for(let i=0;i<N;i++){
   let f=0;for(let j=0;j<N;j++)f+=s.base[i*N+j]*shape[t*N+j];
   f+=s.hoop[i]*(2*shape[t*N+i]-shape[((t+1)%C)*N+i]-shape[((t+C-1)%C)*N+i]);
   const inertial=s.mass[i]*ch.omega2[mode.axial]*shape[t*N+i];residual=Math.max(residual,Math.abs(f-inertial));scale=Math.max(scale,Math.abs(inertial));
  }assert.ok(residual/scale<1e-8,`mode residual ${residual/scale}`);
 }
 // The material is joined across every section boundary and around the seam.
 for(let b=1;b<6;b++)assert.ok(s.base[(b*P-1)*N+b*P]<0);assert.ok(s.hoop.every(k=>k>0));
});
test('whole-body mass and elasticity scale the frequencies physically',()=>{
 const d=initialDesign();for(const b of d.bands)b.mass=4;const mass=compile(d);for(let i=0;i<384;i++)near(mass.sheet.modes[i].frequency*2,compiled.sheet.modes[i].frequency,1e-7);
 for(const b of d.bands){b.mass=1;b.elasticity=4;}const stiff=compile(d);for(let i=0;i<384;i++)near(stiff.sheet.modes[i].frequency,compiled.sheet.modes[i].frequency*2,1e-7);
});
test('moving one ring support and adding volume change global sheet modes',()=>{
 const d=initialDesign();d.bands[1].anchor=.15;const moved=compile(d);assert.notDeepEqual(moved.sheet.modes,compiled.sheet.modes);
 d.blocks.push({x:.55,width:.12,volume:2});assert.notDeepEqual(compile(d).sheet.modes,moved.sheet.modes);
});
test('cached receptance agrees with explicit sums of all 384 modes',()=>{
 const h=1/24000;
 for(const ch of compiled.sheet.channels){const k=receptance(ch,N,24000),q=Float64Array.from({length:N},(_,i)=>Math.sin(i*.7)*1e-7),z=Float64Array.from({length:N},(_,i)=>Math.cos(i*.3)*1e-8),expectedQ=new Float64Array(N),expectedZ=new Float64Array(N);
  for(let m=0;m<N;m++){
   let qm=0,zm=0;for(let i=0;i<N;i++){qm+=ch.phi[i*N+m]*q[i];zm+=ch.phi[i*N+m]*z[i];}
   const w=ch.omega2[m],zn=((1-h*h*w/4)*zm-h*h*w*qm)/(1+h*h*w/4),qn=qm+(zm+zn)/2;
   for(let i=0;i<N;i++){expectedQ[i]+=ch.phi[i*N+m]*qn;expectedZ[i]+=ch.phi[i*N+m]*zn;}
  }
  for(let i=0;i<N;i++){let r=0;for(let j=k.rows[i];j<k.rows[i+1];j++)r+=k.values[j]*(2*q[k.columns[j]]+z[k.columns[j]]);near(r-q[i],expectedQ[i],1e-18);near(2*r-4*q[i]-z[i],expectedZ[i],1e-18);}
  assert.ok(k.droppedRowBound<1e-12);
 }
});
test('rotational compression matches all eight sectors, including contact',()=>{
 const a=new TubeEngine(compiled),b=new TubeEngine(compiled,24000,{symmetry:false});let error=0;
 for(let i=0;i<5000;i++)error=Math.max(error,Math.abs(a.step()-b.step()));assert.ok(error<1e-7,`PCM error ${error}`);assert.equal(a.fault,null);assert.equal(b.fault,null);assert.ok(a.contactCount>0);
 for(let t=0;t<C;t++)for(let i=0;i<N;i++)near(a.displacement[i],b.displacement[t*N+i],1e-9);
});
test('circumferential modes really move and couple to air in the full sheet',()=>{
 const d=initialDesign();d.pressure=0;const c=compile(d),e=new TubeEngine(c,24000,{symmetry:false});
 const mode=c.sheet.modes.findIndex(m=>m.k===2&&m.kind==='sin');e.displaceMode(mode,1e-7);const start=Array.from(e.displacement);
 for(let i=0;i<200;i++)e.step();assert.equal(e.fault,null);assert.ok(e.pressure.some(p=>Math.abs(p)>1e-6));assert.notDeepEqual(Array.from(e.displacement),start);
 assert.ok(e.displacement.some((v,i)=>Math.abs(v-e.displacement[i%N])>1e-10));
 assert.throws(()=>new TubeEngine(c).displaceMode(mode,1e-7),/full sheet/);
});
test('without reservoir pressure, exact rest remains silent',()=>{const d=initialDesign();d.pressure=0;const e=new TubeEngine(compile(d));for(let i=0;i<12000;i++)near(e.step(),0);assert.equal(e.fault,null);assert.ok(e.displacement.every(v=>v===0));});
test('support recompilation preserves displacement and velocity',()=>{const d=initialDesign(),e=new TubeEngine(compiled);for(let i=0;i<1000;i++)e.step();const q=Array.from(e.displacement),v=Array.from(e.velocity);d.bands[1].anchor=.2;e.configure(compile(d));assert.deepEqual(Array.from(e.displacement),q);assert.deepEqual(Array.from(e.velocity),v);for(let i=0;i<1000;i++)e.step();assert.equal(e.fault,null);});
test('same design replays deterministically and produces finite pressure-driven output',()=>{const a=new TubeEngine(compiled),b=new TubeEngine(compiled);let power=0;for(let i=0;i<24000;i++){const x=a.step(),y=b.step();assert.equal(x,y);power+=x*x;}assert.ok(power>0);assert.equal(a.fault,null);});
test('legacy two-wall designs migrate to one ring support without mutating input',()=>{const d=initialDesign();d.version=1;for(const b of d.bands){delete b.offset;delete b.anchor;b.upper=.002;b.lower=.004;b.upperAnchor=.2;b.lowerAnchor=.8;}const next=validateDesign(d);assert.equal(d.version,1);assert.equal(next.version,2);near(next.bands[0].offset,.003);near(next.bands[0].anchor,.5);assert.equal(next.bands[0].upper,undefined);});
test('invalid imported physical parameters are rejected',()=>{const d=initialDesign();d.bands[0].mass=NaN;assert.throws(()=>compile(d));});
test('removing pressure permits ringdown in the same sheet',()=>{const e=new TubeEngine(compiled);let before=0,after=0;for(let i=0;i<24000;i++){const v=e.step();if(i>=18000)before+=v*v;}e.pressureTarget=0;for(let i=0;i<72000;i++){const v=e.step();if(i>=66000)after+=v*v;}assert.equal(e.fault,null);assert.ok(after<before*1e-6,`${after} vs ${before}`);});
test('reconfiguring unchanged geometry preserves every traveling-wave history',()=>{
 const a=new TubeEngine(compiled),b=new TubeEngine(compiled);for(let i=0;i<1500;i++){a.step();b.step();}b.configure(compiled);
 for(let i=0;i<3000;i++)near(a.step(),b.step(),1e-12);
});
