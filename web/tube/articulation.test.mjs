import test from 'node:test';
import {createHash} from 'node:crypto';
import {wav} from './engine.mjs';
import assert from 'node:assert/strict';
import {systemDesign,speechDesign} from './system-presets.mjs';
import {compileSystem,SpeechSystem} from './system.mjs';
import {Articulation,planSpeech} from './articulation.mjs';
import {SpeechSession} from './speech-session.mjs';
import {TubeAudioRenderer} from './audio-renderer.mjs';
const compiled=compileSystem(systemDesign('open'));
function spoken(quantum,level=0){
 const messages=[],renderer=new TubeAudioRenderer(48000,m=>{if(m.type!=='state')messages.push(m);});
 renderer.handle({type:'level',value:level});renderer.handle({type:'speak',compiled,plan:planSpeech('papa')});
 const monitor=new Float32Array(48000*3);for(let i=0;i<monitor.length;i+=quantum)renderer.render(monitor.subarray(i,Math.min(i+quantum,monitor.length)));
 return {messages,monitor,renderer};
}
test('spoken captures and completion agree between worklet and LAN buffer sizes',()=>{
 const a=spoken(128),b=spoken(2048,.05),raw=a.messages.find(m=>m.type==='speech-recording');
 assert.ok(raw);assert.equal(raw.samples.length,Math.ceil(planSpeech('papa').seconds*48000));assert.ok(raw.samples.some(x=>Math.abs(x)>.001));
 assert.deepEqual(raw.samples,b.messages.find(m=>m.type==='speech-recording').samples);
 assert.equal(a.messages.filter(m=>m.type==='speech-finished').length,1);assert.ok(a.monitor.every(x=>x===0));assert.equal(a.renderer.engine.fault,null);
});
test('cancelling before the articulation exists releases pressure without a null dereference',()=>{
 const e=new SpeechSystem(compiled),s=new SpeechSession(e,planSpeech('papa'));s.step();s.cancel();for(let i=0;i<2400;i++)assert.ok(Number.isFinite(s.step()));assert.equal(s.phase,'stopped');assert.ok(e.drive<1);
});
test('lip stop closes the physical aperture and the vowel opens it again',()=>{
 const e=new SpeechSystem(compiled),a=new Articulation(e,planSpeech('papa')),lips=a.parts.lips;
 const gap=()=>Math.min(...Array.from({length:lips.c.counts.at(-1)},(_,j)=>{const i=lips.c.starts.at(-1)+j;return lips.restRadii[i]+lips.displacement[i];}));
 for(let i=0;i<4800;i++)a.step();assert.ok(gap()<.0001,`closed radius ${gap()}`);
 for(let i=0;i<3600;i++)a.step();assert.ok(gap()>.004,`vowel radius ${gap()}`);assert.equal(e.fault,null);
});
test('articulator commands preserve positions before the next physical step',()=>{
 const e=new SpeechSystem(compiled),part=e.engines[1],before=Array.from(part.restRadii,x=>x);part.setSectionRadius(1,.002);
 assert.deepEqual(Array.from(part.restRadii),before);assert.ok(part.radiusTargets[part.c.starts[1]]===.002);
});
test('nasal gestures open a radiating branch while oral stops close it',()=>{
 const a=new Articulation(new SpeechSystem(compiled),planSpeech('mama')),b=new Articulation(new SpeechSystem(compiled),planSpeech('papa'));
 for(let i=0;i<4000;i++){a.step();b.step();}const gap=x=>Math.min(...Array.from({length:x.parts.sinus.c.counts[0]},(_,i)=>x.parts.sinus.restRadii[i]+x.parts.sinus.displacement[i]));assert.ok(gap(a)>.0005);assert.ok(gap(b)<.0001);assert.equal(a.system.fault,null);
});
test('quiet oral air graph is exactly silent without flow or moving targets',()=>{
 const d=systemDesign('open');d.pressure=0;const e=new SpeechSystem(compileSystem(d)),a=new Articulation(e,planSpeech('papa'));
 for(const part of e.engines)part.radiusTargets.set(part.restRadii);
 for(let i=0;i<2400;i++)assert.equal(a.air.step(),0);
});
test('mouth feedback changes source flow through the new bidirectional connection',()=>{
 const make=()=>{const e=new SpeechSystem(compiled);for(let i=0;i<24000;i++)e.step();return new Articulation(e,planSpeech('/aa/'));};
 const a=make(),b=make();b.parts.mouth.setSectionRadius(1,.001);let difference=0;
 for(let i=0;i<6000;i++){a.air.step();b.air.step();difference=Math.max(difference,Math.abs(a.air.glottalFlow-b.air.glottalFlow));}
 assert.ok(difference>1e-9);assert.equal(a.system.fault,null);assert.equal(b.system.fault,null);
});
test('unknown text is rejected rather than silently omitted',()=>{
 assert.throws(()=>planSpeech('abracadabra'),/No pronunciation/);assert.throws(()=>planSpeech('/invalid/'),/Unknown speech sound/);assert.throws(()=>planSpeech(''),/Enter a word/);
});
test('a slit retains its pressure surface as its air opening approaches closure',()=>{
 const d=systemDesign('open'),band=d.tubes[0].design.bands.find(b=>b.role==='cord');band.slotWidth=.01;
 const c=compileSystem(d).tubes[0].compiled,i=c.starts[c.cordIndex];
 assert.equal(c.area[i],2*band.slotWidth*c.r[i]);assert.equal(c.sheet.surface[i],2*band.slotWidth*c.dx[i]/8);
 assert.ok(c.sheet.modes.every(m=>Number.isFinite(m.frequency)&&m.frequency>0));
 const a=new SpeechSystem(compileSystem(d)),b=new SpeechSystem(compileSystem(d),24000,{native:false});
 for(let j=0;j<2400;j++)assert.equal(a.step(),b.step());
});
test('material edits preserve an articulated wall position and ongoing muscle motion',()=>{
 const e=new SpeechSystem(compiled),a=new Articulation(e,planSpeech('papa'));
 for(let i=0;i<1200;i++)a.step();const before=e.engines.map(p=>Array.from(p.restRadii,(r,i)=>r+p.displacement[i]));
 const d=systemDesign('open');d.tubes[0].design.bands[4].mass+=.1;e.configure(compileSystem(d));
 for(let j=0;j<e.engines.length;j++)assert.deepEqual(Array.from(e.engines[j].restRadii,(r,i)=>r+e.engines[j].displacement[i]),before[j]);
});

test('recognized two-word geometry renders with all modes and an unclipped monitor',()=>{
 const c=compileSystem(speechDesign()),plan=planSpeech('papa mama');let capture;
 const renderer=new TubeAudioRenderer(48000,m=>{if(m.type==='speech-recording')capture=m.samples;});
 renderer.handle({type:'level',value:.025});renderer.handle({type:'speak',compiled:c,plan});
 const block=new Float32Array(128);let peak=0;
 for(let n=0;n<1600&&!capture;n++){renderer.render(block);for(const value of block)peak=Math.max(peak,Math.abs(value));}
 assert.ok(capture);assert.equal(createHash('sha256').update(Buffer.from(wav(capture))).digest('hex'),'a5402335395dc68ca5c584410198ec1e744b049b0d97e54a1ed5acda064f4478','preserve the listener-recognized baseline');assert.equal(renderer.engine.fault,null);assert.ok(peak>0&&peak<.95);assert.equal(c.tubes.reduce((n,t)=>n+t.compiled.sheet.modes.length,0),832);
});
