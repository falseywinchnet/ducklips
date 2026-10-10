import {SpeechSession} from './speech-session.mjs';
import {makeEngine,sameTopology} from './system.mjs';
export class TubeAudioRenderer{
 constructor(rate,emit){this.rate=rate;this.emit=emit;this.engine=null;this.speech=null;this.running=false;this.level=.25;this.clock=0;this.phase=0;this.a=0;this.b=0;this.capture=null;this.speechCapture=null;this.scope=new Float32Array(512);this.si=0;}
 handle(d){
 if(d.type==='configure'){if(this.engine&&sameTopology(this.engine,d.compiled)&&(!this.speech||sameAirGeometry(this.engine.c,d.compiled)))this.engine.configure(d.compiled);else {this.speech=null;this.speechCapture=null;this.speechComplete=true;this.engine=makeEngine(d.compiled,this.rate/2);this.a=0;this.b=0;this.phase=0;}this.engine.pressureTarget=this.running?d.compiled.pressure:0;if(this.speech?.speaker)this.speech.speaker.air.body.outletMode='open';}
 if(d.type==='speak'){this.engine=makeEngine(d.compiled,this.rate/2);this.speech=new SpeechSession(this.engine,d.plan);this.running=true;this.a=0;this.b=0;this.phase=0;this.capture=null;this.speechCapture=new Float32Array(Math.ceil(d.plan.seconds*this.rate/2)*2);this.speechCI=0;this.speechPlan=d.plan;this.speechComplete=false;}
 if(d.type==='stop-speech'){if(this.speech){this.speech.cancel();this.speechCapture=null;this.speechComplete=true;this.running=false;}}
 if(d.type==='run'){this.speechCapture=null;if(this.speech){this.engine=makeEngine(this.engine.c,this.rate/2);this.speech=null;this.a=this.b=this.phase=0;}this.running=d.value;if(this.engine)this.engine.pressureTarget=this.running?this.engine.c.pressure:0;}
 if(d.type==='level')this.level=d.value;
 if(d.type==='capture'){this.capture=new Float32Array(Math.round(this.rate*d.seconds));this.ci=0;}
 if(d.type==='reset'){this.speechCapture=null;this.speech=null;this.engine=makeEngine(d.compiled,this.rate/2);this.engine.pressureTarget=this.running?d.compiled.pressure:0;this.a=0;this.b=0;this.phase=0;}
 }
 render(out){for(let i=0;i<out.length;i++){
 let v=0;if(this.engine){if(this.phase===0){this.a=this.b;this.b=this.speech?this.speech.step():this.engine.step();if(this.speech?.finished&&!this.speechComplete){this.speechComplete=true;this.running=false;this.emit({type:'speech-finished'});}}v=this.phase===0?this.a:(this.a+this.b)/2;this.phase=1-this.phase;}
 if(this.speechCapture&&this.speech?.audible){this.speechCapture[this.speechCI++]=v;if(this.speechCI===this.speechCapture.length){this.emit({type:'speech-recording',samples:this.speechCapture,rate:this.rate,plan:this.speechPlan},[this.speechCapture.buffer]);this.speechCapture=null;}}
 out[i]=Math.max(-.95,Math.min(.95,v*this.level)); // Monitor protection only; raw recording is untouched.
 if(this.capture){this.capture[this.ci++]=v;if(this.ci===this.capture.length){this.emit({type:'recording',samples:this.capture,rate:this.rate},[this.capture.buffer]);this.capture=null;}}
 if(this.clock++%8===0)this.scope[this.si++%512]=v;
 }if(this.engine&&this.clock>=this.rate/24){this.clock=0;this.emit({type:'state',...this.engine.snapshot(),speech:this.speech?.snapshot(),scope:this.scope,scopeIndex:this.si%512});}return true;}
}


function sameAirGeometry(a,b){return a.tubes.every((t,i)=>{const x=t.compiled.x,y=b.tubes[i].compiled.x;return x.length===y.length&&x.every((v,j)=>v===y[j]);});}
