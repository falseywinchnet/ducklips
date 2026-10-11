import fs from 'node:fs';
import {createHash} from 'node:crypto';
import {speechDesign} from './system-presets.mjs';
import {compileSystem} from './system.mjs';
import {referenceScore} from './nato-score.mjs';
import {planSpeech} from './articulation.mjs';
import {TubeAudioRenderer} from './audio-renderer.mjs';
import {wav} from './engine.mjs';
const text=process.argv[2]||'papa',dir=process.argv[3]||'out/words-papa';
const design=process.argv[4]?JSON.parse(fs.readFileSync(process.argv[4])):speechDesign(),planFile=process.argv.indexOf('--plan'),plan=planFile>=0?JSON.parse(fs.readFileSync(process.argv[planFile+1])):process.argv[5]==='reference'?referenceScore(text):planSpeech(text);
if(process.argv[5]==='smooth-separation')plan.smoothSeparation=true;
if(process.argv[5]==='contact')plan.contactMode=true;
if(process.argv[5]==='coupled')plan.coupledAir=true;
if(process.argv.includes('--fold-volume'))plan.foldVolume=true;
fs.mkdirSync(dir+'/source/native',{recursive:true});const hashes={};
for(const file of [...fs.readdirSync('web/tube').filter(f=>f.endsWith('.mjs')),'native/core.c','native/build.sh']){const bytes=fs.readFileSync('web/tube/'+file);fs.writeFileSync(dir+'/source/'+file,bytes);hashes[file]=createHash('sha256').update(bytes).digest('hex');}
let samples;const renderer=new TubeAudioRenderer(48000,message=>{if(message.type==='speech-recording')samples=message.samples;});
renderer.handle({type:'level',value:0});renderer.handle({type:'speak',compiled:compileSystem(design),plan});
const block=new Float32Array(128),start=performance.now();
while(!samples){renderer.render(block);if(renderer.engine.fault)throw Error(renderer.engine.fault);if(performance.now()-start>30000)throw Error('Word render did not complete.');}
const elapsed=(performance.now()-start)/1000,peak=samples.reduce((peak,x)=>Math.max(peak,Math.abs(x)),0),listenGain=peak?Math.min(100,.7/peak):1,raw=Buffer.from(wav(samples));
fs.mkdirSync(dir,{recursive:true});fs.writeFileSync(dir+'/outlet.wav',raw);fs.writeFileSync(dir+'/listen.wav',Buffer.from(wav(Float32Array.from(samples,x=>x*listenGain))));
fs.writeFileSync(dir+'/plan.json',JSON.stringify(plan,null,2));fs.writeFileSync(dir+'/design.json',JSON.stringify(design,null,2));
for(const [file,hash] of Object.entries(hashes))if(createHash('sha256').update(fs.readFileSync('web/tube/'+file)).digest('hex')!==hash)throw Error('Source changed during render: '+file);
const receipt={text,seconds:samples.length/48000,preparationSeconds:1.6,wallSeconds:elapsed,totalSimulationRealtime:(plan.seconds+1.6)/elapsed,peak,listenGain,sha256:createHash('sha256').update(raw).digest('hex'),sourceSha256:hashes,recognition:'Not established by rendering or numerical tests.'};
fs.writeFileSync(dir+'/receipt.json',JSON.stringify(receipt,null,2));console.log({text,seconds:receipt.seconds,realtime:receipt.totalSimulationRealtime,peak,sha256:receipt.sha256});
