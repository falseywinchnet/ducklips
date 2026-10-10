import test from 'node:test';import assert from 'node:assert/strict';
import {initialDesign,compile,TubeEngine} from './engine.mjs';
import {systemDesign} from './system-presets.mjs';
import {compileSystem,SpeechSystem,validateSystem} from './system.mjs';
import {TubeAudioRenderer} from './audio-renderer.mjs';

test('compiled sample kernel is identical to JS through variation, contact, edits and fractional transport',()=>{
 const d=initialDesign(),a=new TubeEngine(compile(d)),b=new TubeEngine(compile(d),24000,{native:false});
 for(let phase=0;phase<3;phase++){
  for(let i=0;i<12000;i++)assert.equal(a.step(),b.step());
  assert.equal(a.fault,null);assert.deepEqual(a.snapshot(),b.snapshot());
  d.bands[4].offset=-.00055;d.bands[0].length=phase===0?.1:.05;
  a.configure(compile(d));b.configure(compile(d));
 }
 assert.ok(a.contactCount>0);
});
test('network kernel and JS agree sample for sample, including feedback and state-preserving edits',()=>{
 const d=systemDesign(),c=compileSystem(d),a=new SpeechSystem(c),b=new SpeechSystem(c,24000,{native:false});
 for(let i=0;i<12000;i++)assert.equal(a.step(),b.step());
 const before=a.engines.map(e=>Array.from(e.displacement));d.tubes[1].design.bands[1].offset=-.001;const edited=compileSystem(d);a.configure(edited);b.configure(edited);
 assert.deepEqual(a.engines.map(e=>Array.from(e.displacement)),before);
 for(let i=0;i<12000;i++)assert.equal(a.step(),b.step());assert.deepEqual(a.snapshot(),b.snapshot());assert.equal(a.fault,null);
});
test('connections convey downstream geometry changes back into the original source sheet',()=>{
 const base=systemDesign();
 for(const id of ['mouth','lips','sinus','undertongue']){
  const changed=structuredClone(base);changed.tubes.find(t=>t.id===id).design.bands[0].radius*=.5;
  const a=new SpeechSystem(compileSystem(base)),b=new SpeechSystem(compileSystem(changed));let difference=0;
  for(let i=0;i<16000;i++){a.step();b.step();if(i>8000)difference+=Math.abs(a.engines[0].pressure[32]-b.engines[0].pressure[32]);}
  assert.equal(a.fault,null);assert.equal(b.fault,null);assert.ok(difference>1e-3,id+' returns no source feedback');
 }
});
test('junction update is independent of child evaluation order',()=>{
 const d=systemDesign(),reverse=structuredClone(d);reverse.tubes.reverse();const a=new SpeechSystem(compileSystem(d)),b=new SpeechSystem(compileSystem(reverse));
 for(let i=0;i<12000;i++)assert.equal(a.step(),b.step());
});
test('all arrangements retain the accepted source and sorted complete child modes',()=>{
 for(const name of ['open','balanced','rounded']){const d=systemDesign(name),c=compileSystem(d);assert.deepEqual(d.tubes[0].design,initialDesign());assert.equal(new Set(d.tubes.map(t=>t.design.noiseSeed)).size,5);for(const t of c.tubes){assert.equal(t.compiled.sheet.modes.length,t.compiled.n*8);assert.ok(t.compiled.sheet.modes.length>=40);}}
 assert.equal(systemDesign('open').tubes[2].design.bands.length,3);assert.equal(systemDesign('rounded').tubes[1].design.bands.length,3);
});
test('unforced network stays at exact rest; pressure-off network releases',()=>{
 const d=systemDesign();d.pressure=0;const quiet=new SpeechSystem(compileSystem(d));for(let i=0;i<6000;i++)assert.equal(quiet.step(),0);
 d.pressure=850;const e=new SpeechSystem(compileSystem(d));let before=0,after=0;
 for(let i=0;i<48000;i++){const x=e.step();if(i>42000)before+=x*x;}
 e.pressureTarget=0;for(let i=0;i<144000;i++){const x=e.step();if(i>138000)after+=x*x;}
 assert.equal(e.fault,null);assert.ok(before>0);assert.ok(after<before*1e-5,`${after}/${before}`);
});
test('invalid network terminations and missing branch connections are rejected',()=>{
 for(const mutate of [d=>d.connections.pop(),d=>d.tubes[0].end='open',d=>d.connections[0].a.end='inlet',d=>d.connections.push(d.connections[0]),d=>d.tubes[2].end='closed']){const d=systemDesign();mutate(d);assert.throws(()=>validateSystem(d));}
});
test('network audio is buffer-size independent for LAN and worklet, and topology resets cleanly',()=>{
 function render(size){let capture;const r=new TubeAudioRenderer(48000,m=>{if(m.type==='recording')capture=m.samples;});r.handle({type:'configure',compiled:compileSystem(systemDesign())});r.handle({type:'run',value:true});r.handle({type:'capture',seconds:.2});for(let i=0;i<9600;i+=size)r.render(new Float32Array(Math.min(size,9600-i)));return {r,capture};}
 const a=render(128),b=render(2048);assert.deepEqual(a.capture,b.capture);const r=a.r;r.handle({type:'configure',compiled:compileSystem(systemDesign('open'))});assert.equal(r.phase,0);assert.equal(r.engine.engines[2].n,12);r.render(new Float32Array(4000));assert.equal(r.engine.fault,null);
});
