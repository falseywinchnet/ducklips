import {TubeEngine} from './engine.mjs';
class TubeProcessor extends AudioWorkletProcessor{
 constructor(){super();this.engine=null;this.running=false;this.level=.25;this.clock=0;this.phase=0;this.a=0;this.b=0;this.capture=null;this.scope=new Float32Array(512);this.si=0;this.port.onmessage=({data:d})=>{
 if(d.type==='configure'){if(this.engine)this.engine.configure(d.compiled);else this.engine=new TubeEngine(d.compiled,sampleRate/2);this.engine.pressureTarget=this.running?d.compiled.pressure:0;}
 if(d.type==='run'){this.running=d.value;if(this.engine)this.engine.pressureTarget=this.running?this.engine.c.pressure:0;}
 if(d.type==='level')this.level=d.value;
 if(d.type==='capture'){this.capture=new Float32Array(Math.round(sampleRate*d.seconds));this.ci=0;}
 if(d.type==='reset'){this.engine=new TubeEngine(d.compiled,sampleRate/2);this.engine.pressureTarget=this.running?d.compiled.pressure:0;this.a=0;this.b=0;this.phase=0;}
 };}
 process(inputs,outputs){const out=outputs[0][0];for(let i=0;i<out.length;i++){
 let v=0;if(this.engine){if(this.phase===0){this.a=this.b;this.b=this.engine.step();}v=this.phase===0?this.a:(this.a+this.b)/2;this.phase=1-this.phase;}
 out[i]=Math.max(-.95,Math.min(.95,v*this.level)); // Monitor protection only; raw recording is untouched.
 if(this.capture){this.capture[this.ci++]=v;if(this.ci===this.capture.length){this.port.postMessage({type:'recording',samples:this.capture,rate:sampleRate},[this.capture.buffer]);this.capture=null;}}
 if(this.clock++%8===0)this.scope[this.si++%512]=v;
 }if(this.engine&&this.clock>=sampleRate/24){this.clock=0;this.port.postMessage({type:'state',...this.engine.snapshot(),scope:this.scope,scopeIndex:this.si%512});}return true;}
}
registerProcessor('flat-tube',TubeProcessor);
