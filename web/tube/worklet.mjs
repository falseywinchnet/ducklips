import {TubeAudioRenderer} from './audio-renderer.mjs';
class TubeProcessor extends AudioWorkletProcessor {
  constructor(){
    super();
    this.renderer=new TubeAudioRenderer(sampleRate,(data,transfer)=>this.port.postMessage(data,transfer));
    this.port.onmessage=({data})=>this.renderer.handle(data);
  }
  process(inputs,outputs){return this.renderer.render(outputs[0][0]);}
}
registerProcessor('flat-tube',TubeProcessor);
