import {makeEngine,sameTopology} from './system.mjs';
export class TubeAudioRenderer{
 constructor(rate,emit){this.rate=rate;this.emit=emit;this.engine=null;this.running=false;this.level=.25;this.clock=0;this.phase=0;this.a=0;this.b=0;this.capture=null;this.scope=new Float32Array(512);this.si=0;}
 handle(d){
 if(d.type==='configure'){if(this.engine&&sameTopology(this.engine,d.compiled))this.engine.configure(d.compiled);else {this.engine=makeEngine(d.compiled,this.rate/2);this.a=0;this.b=0;this.phase=0;}this.engine.pressureTarget=this.running?d.compiled.pressure:0;}
 if(d.type==='run'){this.running=d.value;if(this.engine)this.engine.pressureTarget=this.running?this.engine.c.pressure:0;}
 if(d.type==='level')this.level=d.value;
 if(d.type==='capture'){this.capture=new Float32Array(Math.round(this.rate*d.seconds));this.ci=0;}
 if(d.type==='reset'){this.engine=makeEngine(d.compiled,this.rate/2);this.engine.pressureTarget=this.running?d.compiled.pressure:0;this.a=0;this.b=0;this.phase=0;}
 }
 render(out){for(let i=0;i<out.length;i++){
 let v=0;if(this.engine){if(this.phase===0){this.a=this.b;this.b=this.engine.step();}v=this.phase===0?this.a:(this.a+this.b)/2;this.phase=1-this.phase;}
 out[i]=Math.max(-.95,Math.min(.95,v*this.level)); // Monitor protection only; raw recording is untouched.
 if(this.capture){this.capture[this.ci++]=v;if(this.ci===this.capture.length){this.emit({type:'recording',samples:this.capture,rate:this.rate},[this.capture.buffer]);this.capture=null;}}
 if(this.clock++%8===0)this.scope[this.si++%512]=v;
 }if(this.engine&&this.clock>=this.rate/24){this.clock=0;this.emit({type:'state',...this.engine.snapshot(),scope:this.scope,scopeIndex:this.si%512});}return true;}
}

