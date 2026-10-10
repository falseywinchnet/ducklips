import {TubeAudioRenderer} from './audio-renderer.mjs';
let renderer;
self.onmessage=({data})=>{
 try{
  if(data.type==='init')renderer=new TubeAudioRenderer(data.rate,(message,transfer)=>self.postMessage(message,transfer));
  else if(data.type==='pull'){
   const samples=new Float32Array(2048);
   renderer.render(samples);
   self.postMessage({type:'audio',samples},[samples.buffer]);
  }else renderer.handle(data);
 }catch(error){self.postMessage({type:'error',message:error.message});}
};
