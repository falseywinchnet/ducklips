import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {initialDesign,compile,TubeEngine,MATERIALS,C} from './engine.mjs';
import {referenceDesign} from './presets.mjs';
import {receptance,VariableReceptance,compileSheet as compilePhysicalSheet} from './sheet.mjs';
import {ElasticityVariation} from './variation.mjs';
import {TubeAudioRenderer} from './audio-renderer.mjs';
const near=(a,b,t=1e-10)=>assert.ok(Math.abs(a-b)<=t,`${a} vs ${b}; tolerance ${t}`);

test('saved six-section reference matches the preserved user file exactly',()=>{
 assert.deepEqual(referenceDesign(),JSON.parse(fs.readFileSync(new URL('./presets/user-vocal-tube-2026-10-10.json',import.meta.url))));
});
test('below-mouth geometry preserves diameter and material constraints in a connected 480-mode sheet',()=>{
 const d=initialDesign(),c=compile(d),r=b=>b.radius+b.offset;
 assert.equal(d.bands.length,9);assert.equal(c.n,60);assert.equal(c.sheet.modes.length,480);assert.equal(c.sheet.channels[0].omega2.length,60);
 assert.equal(c.cordIndex,4);assert.equal(d.bands[4].elasticity,1.5);assert.equal(d.blocks.length,0);
 for(const b of [d.bands[0],d.bands.at(-1)]){near(r(b),2*r(d.bands[5]));assert.equal(b.mass,1);assert.equal(b.elasticity,1);assert.equal(b.support,1);}
 for(const b of d.bands.filter(b=>b.role==='adapter')){assert.equal(b.mass,.2);assert.equal(b.support,.2);assert.equal(b.elasticity,3);near(r(b),r(d.bands[5]));}
 for(const i of c.starts.slice(1))assert.ok(c.sheet.base[(i-1)*c.n+i]<0);
 const e=new TubeEngine(c);near(e.snapshot().gap,2*r(d.bands[4]));
});
test('independent slow material states are bounded, seeded, and continuous through bubble resets',()=>{
 const a=new ElasticityVariation(9,1),b=new ElasticityVariation(9,1),other=new ElasticityVariation(9,2);let maxStep=0,previous=new Float64Array(9).fill(1);
 for(let i=0;i<20000;i++){
  const f=a.advance(.005,.1);assert.deepEqual(f,b.advance(.005,.1));other.advance(.005,.1);
  for(let j=0;j<9;j++){assert.ok(f[j]>=.9&&f[j]<=1.1);maxStep=Math.max(maxStep,Math.abs(f[j]-previous[j]));}previous.set(f);
 }
 assert.notDeepEqual(a.factors,other.factors);assert.ok(new Set(a.factors).size===9);assert.ok(maxStep<.01);
 a.advance(.005,0);assert.ok(a.factors.every(x=>x===1));
});
test('varying cached response agrees with independent complete eigenanalysis at both elasticity limits',()=>{
 for(const master of [1,4]){
  const d=initialDesign();d.elasticityScale=master;const c=compile(d),n=c.n;
  for(const factors of [d.bands.map((_,i)=>i%2?1.1:.9),d.bands.map((_,i)=>i%2?.9:1.1)]){
   const changed=structuredClone(d);changed.bands.forEach((b,i)=>b.elasticity*=factors[i]);
   // A local value of .2 is a UI limit; the live 10% variation can go to .18.

   // Build independently from physical values, bypassing only UI range validation.
   const target=compilePhysicalSheet(c,changed.bands.map(b=>({...b,elasticity:b.elasticity*master})),MATERIALS,8);
   for(const ch of [c.sheet.channels[0],c.sheet.channels.at(-1)]){
    const base=receptance(ch,n,24000),v=new VariableReceptance(c.sheet,ch,24000,base);v.update(factors,0);
    const expected=receptance(target.channels.find(x=>x.k===ch.k&&x.kind===ch.kind),n,24000,{dense:true});
    for(let i=0;i<n*n;i++)near(v.matrix[i],expected.matrix[i],3e-12);
    for(let i=0;i<n;i++)near(v.viscous[i],target.viscous[i],1e-10);
    assert.ok(v.droppedRowBound<1e-12);
   }
  }
 }
});
test('packed air histories reproduce the original transport through reconfiguration and fractional delays',()=>{
 const d=initialDesign(),c=compile(d),a=new TubeEngine(c),b=new TubeEngine(c,24000,{compactAir:false});
 for(let phase=0;phase<3;phase++){
  for(let i=0;i<6000;i++)assert.equal(a.step(),b.step());
  if(phase===0){d.bands.forEach(x=>{x.length=.1;});}
  else {d.bands=initialDesign().bands;}
  const next=compile(d);a.configure(next);b.configure(next);
 }
 assert.equal(a.fault,null);assert.equal(b.fault,null);
});
test('unchanged and pressure-only edits preserve the evolving material state and wave histories',()=>{
 const d=initialDesign(),c=compile(d),a=new TubeEngine(c),b=new TubeEngine(c);
 for(let i=0;i<1501;i++){a.step();b.step();}b.configure(compile(d));
 for(let i=0;i<1200;i++)assert.equal(a.step(),b.step());
 d.pressure=700;b.configure(compile(d));a.pressureTarget=700;
 for(let i=0;i<1200;i++)assert.equal(a.step(),b.step());
 d.noiseSeed++;b.configure(compile(d));assert.equal(b.variationSeed,d.noiseSeed);
});
test('whole-body elasticity covers 1 through 4, with doubled dry frequencies and finite running output',()=>{
 let baseline;
 for(const scale of [1,2,3,4]){
  const d=initialDesign();d.elasticityScale=scale;const c=compile(d),e=new TubeEngine(c);
  if(!baseline)baseline=c.sheet.modes;
  for(let i=0;i<c.sheet.modes.length;i++)near(c.sheet.modes[i].frequency,baseline[i].frequency*Math.sqrt(scale),1e-7);
  for(let i=0;i<48000;i++)e.step();assert.equal(e.fault,null);
 }
});
test('expanded varying sheet retains rotational equivalence',()=>{
 const c=compile(initialDesign()),a=new TubeEngine(c),b=new TubeEngine(c,24000,{symmetry:false});let error=0;
 for(let i=0;i<4000;i++)error=Math.max(error,Math.abs(a.step()-b.step()));assert.ok(error<1e-7);assert.equal(a.fault,null);assert.equal(b.fault,null);
 for(let t=0;t<C;t++)for(let i=0;i<c.n;i++)near(a.displacement[i],b.displacement[t*c.n+i],1e-9);
});
test('changing section count resets audio interpolation and retains the requested running state',()=>{
 const r=new TubeAudioRenderer(48000,()=>{});r.handle({type:'configure',compiled:compile(referenceDesign())});r.handle({type:'run',value:true});r.render(new Float32Array(501));
 assert.equal(r.phase,1);r.handle({type:'configure',compiled:compile(initialDesign())});assert.equal(r.phase,0);assert.equal(r.a,0);assert.equal(r.b,0);assert.equal(r.engine.n,60);assert.equal(r.engine.pressureTarget,initialDesign().pressure);
 r.render(new Float32Array(1000));assert.equal(r.engine.fault,null);
});
test('material variation cannot excite exact rest, and output releases when reservoir pressure stops',()=>{
 const quiet=initialDesign();quiet.pressure=0;const rest=new TubeEngine(compile(quiet));
 for(let i=0;i<6000;i++)assert.equal(rest.step(),0);
 const e=new TubeEngine(compile(initialDesign()));let before=0,after=0;
 for(let i=0;i<48000;i++){const x=e.step();if(i>=42000)before+=x*x;}
 e.pressureTarget=0;for(let i=0;i<96000;i++){const x=e.step();if(i>=90000)after+=x*x;}
 assert.equal(e.fault,null);assert.ok(before>0);assert.ok(after<before*1e-6);
});
