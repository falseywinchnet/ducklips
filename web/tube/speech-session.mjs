import {Articulation} from './articulation.mjs';
export class SpeechSession {
 constructor(engine,plan){if(engine.c.kind!=='network')throw Error('Choose a connected mouth arrangement to speak.');this.engine=engine;this.plan=plan;this.count=0;this.phase='preparing';this.audible=false;this.finished=false;this.speaker=null;engine.pressureTarget=engine.c.pressure;}
 cancel(){this.cancelled=true;this.finished=true;this.phase='stopped';this.engine.pressureTarget=0;}
 step(){
  if(this.cancelled)return this.speaker?this.speaker.air.step():this.engine.step();
  if(this.phase==='preparing'){
   this.engine.step();if(++this.count>=this.engine.rate){this.speaker=new Articulation(this.engine,this.plan);this.count=0;this.phase='forming closure';}return 0;
  }
  if(this.phase==='forming closure'){
   this.speaker.air.step();if(++this.count>=this.engine.rate*.6){this.phase='speaking';}return 0;
  }
  this.audible=true;const value=this.speaker.step();if(this.speaker.finished){this.finished=true;this.phase='finished';}return value;
 }
 snapshot(){return {...this.speaker?.snapshot(),text:this.plan.text,phase:this.phase,finished:this.finished};}
}
