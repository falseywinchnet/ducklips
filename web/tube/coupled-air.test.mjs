import test from 'node:test';
import assert from 'node:assert/strict';
import {compileSystem,SpeechSystem} from './system.mjs';
import {speechDesign} from './system-presets.mjs';
import {CoupledAir,expansionLoss} from './coupled-air.mjs';
import {SpeechSession} from './speech-session.mjs';
import {referenceScore} from './nato-score.mjs';
const compiled=compileSystem(speechDesign());
test('the serial airway has one wall owner and no additional glottal takeoff',()=>{
 const s=new SpeechSystem(compiled),ports=s.engines[0].ports.length,a=new CoupledAir(s);assert.equal(a.body.ports.length,ports);
 for(const e of s.engines)for(let i=0;i<e.n;i++){
  const count=a.nodes.reduce((sum,n)=>sum+n.walls.filter(w=>w.e===e&&w.i===i).length,0),fold=e===a.body&&i>=a.begin&&i<a.end;
  assert.equal(count,fold?0:1,'each patch belongs to either its airway junction or the fold element');
 }
});
test('unforced serial sheets remain at exact rest',()=>{
 const s=new SpeechSystem(compiled),a=new CoupledAir(s);a.flowNoise=true;s.pressureTarget=0;
 for(let i=0;i<2400;i++)assert.equal(a.step(),0);assert.equal(s.fault,null);
});
test('batching all sheets preserves every audible sample',()=>{
 const p={...referenceScore('papa'),foldVolume:true},a=new SpeechSession(new SpeechSystem(compiled),p),b=new SpeechSession(new SpeechSystem(compiled),{...p,nativeBatch:false});
 while(!a.finished){assert.equal(a.step(),b.step());assert.equal(a.engine.fault,null);assert.equal(b.engine.fault,null);}
 assert.ok(a.speaker.air.maxVolumeError<1e-12,`fold volume residual ${a.speaker.air.maxVolumeError}`);
});
test('a vowel-initial reference word does not sound during its initial silence',()=>{
 const p=referenceScore('alpha'),s=new SpeechSession(new SpeechSystem(compiled),p);let audible=0;
 while(audible<Math.floor(p.events[0].start*s.engine.rate)-2){const sample=s.step();if(s.audible){assert.equal(sample,0);audible++;}}
});
test('air attenuation preserves a stationary pressure field instead of leaking DC',()=>{
 const s=new SpeechSystem(compiled),a=new CoupledAir(s),pressure=500;a.airLossFIR=true;s.pressureTarget=pressure;a.body.drive=pressure;
 for(const n of a.nodes){n.open=false;n.walls=[];}
 for(const e of a.edges){const radius=Math.min(e.wall.e.c.r[e.wall.i],...e.apertures.map(w=>w.e.c.r[w.i]),e.aperture?e.aperture.e.c.r[e.aperture.i]:Infinity),wave=pressure*a.rootY*radius/2;e.r=e.l=e.lossR1=e.lossR2=e.lossL1=e.lossL2=wave;}
 a.step();assert.equal(s.fault,null);
 for(const e of a.edges){const expected=pressure*e.root/2;assert.ok(Math.abs(e.r-expected)<1e-12);assert.ok(Math.abs(e.l-expected)<1e-12);}
});

test('geometric oral loss vanishes in a uniform pipe and has no aperture cutoff',()=>{
 assert.equal(expansionLoss(.001,.001),0);assert.equal(expansionLoss(.01,.01),0);
 for(const r of [.0001,.0005,.001,.0015,.002,.005]){
  const k=expansionLoss(r,.008);assert.ok(k>=0&&k<=1);
  assert.ok(Math.abs(expansionLoss(r-1e-10,.008)-expansionLoss(r+1e-10,.008))<1e-6);
 }
});

test('a prescribed inlet holds the applied pressure under wall and flow loading',()=>{
 const s=new SpeechSystem(compiled),a=new CoupledAir(s);a.pressureBoundary=true;s.pressureTarget=200;
 for(let i=0;i<240;i++){a.step();assert.equal(a.nodes[a.inlet].p,a.body.drive);assert.equal(s.fault,null);}
});
