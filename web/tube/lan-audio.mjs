// HTTP LAN transport: run the identical renderer in a Worker and schedule PCM
// buffers through standard Web Audio. No insecure-origin flag or legacy DSP API.
export function createLanAudio(context){
 const worker=new Worker('./audio-worker.mjs',{type:'module'});
 let destination=null,nextTime=0,connected=false;
 const port={onmessage:null,postMessage(message){worker.postMessage(message);}};
 worker.onmessage=({data})=>{
  if(data.type==='audio'){
   const buffer=context.createBuffer(1,data.samples.length,context.sampleRate);
   buffer.copyToChannel(data.samples,0);
   const source=context.createBufferSource();source.buffer=buffer;source.connect(destination);
   const start=Math.max(nextTime,context.currentTime+.015);nextTime=start+buffer.duration;
   source.onended=()=>{source.disconnect();if(connected)worker.postMessage({type:'pull'});};
   source.start(start);
  }else if(data.type==='error'){
   port.onmessage?.({data:{type:'state',pressure:[],displacement:[],velocity:[],scope:new Float32Array(512),scopeIndex:0,drive:0,gap:0,fault:data.message}});
  }else port.onmessage?.({data});
 };
 worker.postMessage({type:'init',rate:context.sampleRate});
 return {port,connect(target){
  destination=target;connected=true;
  // Queue only after the caller's configuration has entered the worker mailbox.
  queueMicrotask(()=>{for(let i=0;i<4;i++)worker.postMessage({type:'pull'});});
 }};
}
