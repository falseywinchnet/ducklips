import test from 'node:test';
import assert from 'node:assert/strict';
import {TubeAudioRenderer} from './audio-renderer.mjs';
import {initialDesign,compile} from './engine.mjs';
function render(quantum,level){
 let recording;
 const renderer=new TubeAudioRenderer(48000,message=>{if(message.type==='recording')recording=message.samples;});
 renderer.handle({type:'configure',compiled:compile(initialDesign())});
 renderer.handle({type:'level',value:level});renderer.handle({type:'run',value:true});
 renderer.handle({type:'capture',seconds:1});
 const monitor=new Float32Array(48000);
 for(let i=0;i<monitor.length;i+=quantum)renderer.render(monitor.subarray(i,Math.min(i+quantum,monitor.length)));
 return {monitor,recording};
}
test('worklet and LAN buffer sizes produce identical PCM and raw captures',()=>{
 const worklet=render(128,.25),lan=render(2048,.25);
 assert.deepEqual(lan.monitor,worklet.monitor);assert.deepEqual(lan.recording,worklet.recording);
});
test('LAN raw capture retains engine output while monitor is muted',()=>{
 const result=render(2048,0);assert.ok(result.monitor.every(x=>x===0));
 assert.ok(result.recording.some(x=>Math.abs(x)>.001));assert.equal(result.recording.length,48000);
});
