import {CoupledAir} from './coupled-air.mjs';
import {OralAir} from './oral-air.mjs';
// A stream of articulatory targets, in metres and seconds. The moving supports
// supply work to the existing pressure-driven sheets; no PCM is stored here.
export const PHONES={
 aa:{mouth:[.0054,.011,.011,.008],lip:.006,voice:1,duration:.28},
 ae:{mouth:[.006,.009,.007,.009],lip:.009,voice:1,duration:.23},
 eh:{mouth:[.008,.008,.004,.008],lip:.007,voice:1,duration:.20},
 ih:{mouth:[.010,.009,.003,.006],lip:.006,voice:1,duration:.17},
 iy:{mouth:[.012,.010,.0018,.005],lip:.006,voice:1,duration:.25},
 uh:{mouth:[.005,.004,.008,.007],lip:.004,voice:1,duration:.18},
 uw:{mouth:[.010,.004,.008,.006],lip:.0025,voice:1,duration:.28},
 ao:{mouth:[.004,.007,.010,.007],lip:.004,voice:1,duration:.26},
 ax:{mouth:[.006,.006,.006,.006],lip:.006,voice:1,duration:.16},
 p:{closure:'lip',voice:0,duration:.12},b:{closure:'lip',voice:1,duration:.10},
 t:{closure:'tip',voice:0,duration:.10},d:{closure:'tip',voice:1,duration:.09},
 k:{closure:'back',voice:0,duration:.12},g:{closure:'back',voice:1,duration:.10},
 m:{closure:'lip',voice:1,nasal:1,duration:.14},n:{closure:'tip',voice:1,nasal:1,duration:.13},
 h:{voice:0,aspiration:1,duration:.08},
 f:{constriction:'lip',voice:0,duration:.17},v:{constriction:'lip',voice:1,duration:.14},
 s:{constriction:'tip',voice:0,duration:.18},z:{constriction:'tip',voice:1,duration:.15},
 sh:{constriction:'blade',voice:0,duration:.18},zh:{constriction:'blade',voice:1,duration:.16},ng:{closure:'back',voice:1,nasal:1,duration:.13},
 l:{mouth:[.007,.008,.003,.003],lip:.006,voice:1,duration:.13},
 r:{mouth:[.006,.009,.003,.004],lip:.004,voice:1,duration:.15},
 w:{mouth:[.010,.004,.008,.006],lip:.0025,voice:1,duration:.12},
 y:{mouth:[.012,.010,.0018,.005],lip:.006,voice:1,duration:.10},
 _: {voice:0,silence:true,duration:.16}
};
export const WORDS={papa:'p aa p aa',mama:'m aa m aa',baba:'b aa b aa',baby:'b ey b iy',poppy:'p aa p iy',puppy:'p uh p iy',paper:'p ey p er',people:'p iy p ax l',hello:'h eh l ow',up:'uh p',tea:'t iy',two:'t uw',go:'g ow',no:'n ow',yes:'y eh s',we:'w iy',you:'y uw',me:'m iy',pop:'p aa p',mom:'m aa m',dad:'d ae d',day:'d ey',bye:'b ay',wow:'w aw',
 alpha:'ae l f aa',bravo:'b r aa v ow',charlie:'t sh aa r l iy',delta:'d eh l t aa',echo:'eh k ow',foxtrot:'f aa k s t r aa t',golf:'g ao l f',hotel:'h ow t eh l',india:'ih n d iy aa',juliet:'d z uw l iy eh t',kilo:'k iy l ow',lima:'l iy m aa',mike:'m ay k',november:'n ow v eh m b er',oscar:'aa s k er',quebec:'k ax b eh k',romeo:'r ow m iy ow',sierra:'s iy eh r aa',tango:'t ae n g ow',uniform:'y uw n ih f ao r m',victor:'v ih k t er',whiskey:'w ih s k iy','x-ray':'eh k s r ey',yankee:'y ae n k iy',zulu:'z uw l uw'};
const diphthongs={ey:['eh','iy'],ay:['aa','iy'],ow:['ao','uw'],aw:['aa','uw'],oy:['ao','iy'],er:['ax','r']};
export function planSpeech(text,{speed=1}={}){
 if(typeof text!=='string'||text.length>500||!Number.isFinite(speed)||speed<.5||speed>2)throw Error('Use at most 500 characters and speed 0.5–2.');
 const words=text.toLowerCase().trim().split(/\s+/).filter(Boolean);if(!words.length)throw Error('Enter a word.');
 const events=[];let at=.08;
 for(const word of words){const phones=word.startsWith('/')&&word.endsWith('/')?word.slice(1,-1).split('-'):WORDS[word.replace(/[.,!?]/g,'')];if(!phones)throw Error(`No pronunciation for “${word}” yet. Use a listed word or /p-aa-p-aa/.`);
  for(const token of Array.isArray(phones)?phones:phones.split(' ')){const parts=diphthongs[token]||[token];for(const phone of parts){const spec=PHONES[phone];if(!spec)throw Error(`Unknown speech sound: ${phone}`);const duration=(diphthongs[token]?.length? .15:spec.duration)/speed;events.push({phone,start:at,end:at+duration,word});at+=duration;}}
  events.push({phone:'_',start:at,end:at+.18/speed,word});at+=.18/speed;
 }
 return {text,speed,events,seconds:at+.20};
}
export class Articulation {
 constructor(system,plan){if(system.c.kind!=='network'||!['mouth','lips','sinus','undertongue'].every(id=>system.c.tubes.some(t=>t.id===id)))throw Error('Speech requires mouth, lips, nasal, and under-tongue stages.');this.air=plan.coupledAir?new CoupledAir(system,{batch:plan.nativeBatch!==false}):new OralAir(system);this.air.smoothSeparation=!!plan.smoothSeparation;this.air.foldVolume=!!plan.foldVolume;this.air.flowNoise=!!plan.flowNoise;this.air.apertureRadiation=!!plan.apertureRadiation;this.air.airLossFIR=!!plan.airLossFIR;this.air.geometricJetLoss=!!plan.geometricJetLoss;this.air.pressureBoundary=!!plan.pressureBoundary;this.system=system;this.plan=plan;this.time=0;this.index=0;this.finished=false;this.control=0;this.phone='_';this.voice=0;
  this.parts=Object.fromEntries(system.c.tubes.map((t,i)=>[t.id,system.engines[i]]));
  if(!this.parts.mouth||!this.parts.lips)throw Error('Speech requires a connected mouth and lips.');
  this.targets=system.engines.map(e=>e.c.starts.map(i=>e.c.r[i]));this.current=this.targets.map(a=>a.slice());
  this.update(0);
 }
 step(){
  const s=this.system,t=this.time;if(this.control++%48===0)this.update(t);
  const v=this.air.step();this.time+=1/s.rate;if(this.time>=this.plan.seconds)this.finished=true;return v;
 }
 update(t){
  const {plan,system:s}=this;while(this.index<plan.events.length&&t>=plan.events[this.index].end)this.index++;
  const event=plan.events[this.index],prelude=plan.honorOnsets&&event&&t<event.start,spec=event&&!prelude?PHONES[event.phone]:PHONES._;this.phone=prelude?'_':event?.phone??'_';this.air.nasalTarget=spec.nasal?1:0;
  if(spec.silence){
   s.pressureTarget=0;
   if(plan.releaseOnSilence&&!(prelude&&plan.preparePosture)){const body=s.engines[s.c.source];body.setSectionRadius(body.c.cordIndex,.0015);}
   return;
  }
  let vowel=spec.mouth?spec:null;
  if(!vowel)for(let i=this.index+1;i<Math.min(this.index+5,plan.events.length);i++){const next=PHONES[plan.events[i].phone];if(next.mouth){vowel=next;break;}}
  vowel??=PHONES.ax;
  const mi=s.c.tubes.findIndex(x=>x.id==='mouth'),li=s.c.tubes.findIndex(x=>x.id==='lips'),bi=s.c.source;
  for(let j=0;j<this.targets[mi].length;j++)this.targets[mi][j]=s.engines[mi].c.r[s.engines[mi].c.starts[j]]*vowel.mouth[Math.min(j,3)]/PHONES.aa.mouth[Math.min(j,3)];
  this.targets[li]=s.engines[li].c.starts.map(i=>s.engines[li].c.r[i]);this.targets[li][this.targets[li].length-1]*=vowel.lip/PHONES.aa.lip;const lip=this.targets[li].length-1;
  const lipClosure=spec.closure==='lip';
  if(!plan.localClosures&&spec.closure==='tip')this.targets[mi][this.targets[mi].length-1]=-.001;
  if(!plan.localClosures&&spec.closure==='back')this.targets[mi][1]=-.001;
  if(spec.constriction==='lip'&&!plan.oralJetAreas)this.targets[li][lip]=.00035;
  if(!plan.localClosures&&spec.constriction==='tip')this.targets[mi][this.targets[mi].length-1]=.0003;
  if(!plan.localClosures&&spec.constriction==='blade')this.targets[mi][2]=.0004;
  const body=s.engines[bi],cord=body.c.cordIndex;
  const nextSpec=plan.events[this.index+1]&&PHONES[plan.events[this.index+1].phone];
  const preparingVoice=plan.prephonationLead&&(spec.closure||spec.constriction||spec.aspiration)&&nextSpec?.voice&&event.end-t<=plan.prephonationLead;
  this.air.breathy=!(spec.voice||preparingVoice);
  this.targets[bi][cord]=this.air.breathy?(plan.unvoicedGap??body.c.r[body.c.starts[cord]]+.00015):(plan.foldRadius??(body.c.slotWidth[body.c.starts[cord]]?0:.0003));
  this.air.constriction=spec.constriction??null;
  const ni=s.c.tubes.findIndex(x=>x.id==='sinus');this.targets[ni][0]=spec.nasal?s.engines[ni].c.r[0]:-.0003;
  if(plan.articulateFloor){
   const ui=s.c.tubes.findIndex(x=>x.id==='undertongue');
   this.targets[ui][0]=this.phone==='r'?s.engines[ui].c.r[0]:-.0003;
  }
  const progress=Math.max(0,Math.min(1,(t-event.start)/(event.end-event.start)));
  s.pressureTarget=s.c.pressure*((event.pressure??1)+((event.pressureEnd??event.pressure??1)-(event.pressure??1))*progress);
  if(event.effort!==undefined){const effort=event.effort+((event.effortEnd??event.effort)-event.effort)*progress;if(plan.effortScope==='body')for(let j=0;j<body.c.counts.length;j++)body.setSectionElasticity(j,effort);else body.setSectionElasticity(cord,effort);}
  const blend=1;
  for(let i=0;i<s.engines.length;i++)for(let j=0;j<this.current[i].length;j++){this.current[i][j]+=blend*(this.targets[i][j]-this.current[i][j]);s.engines[i].setSectionRadius(j,this.current[i][j]);}
  const vowelKey=plan.vowelAreas&&Object.keys(PHONES).find(key=>PHONES[key]===vowel);
  if(vowelKey&&plan.vowelAreas[vowelKey])this.applyAreaProfile(plan.vowelAreas[vowelKey]);
  if(plan.localClosures&&(spec.closure==='tip'||spec.closure==='back'||spec.constriction==='tip'||spec.constriction==='blade')){
   const e=this.parts.mouth,j=spec.closure==='back'?1:spec.constriction==='blade'?2:e.c.starts.length-1;
   const centerIndex=spec.closure==='back'?e.c.starts[j]+Math.floor(e.c.counts[j]/2):e.c.starts[j]+e.c.counts[j]-1,center=e.c.x[centerIndex],base=this.targets[mi][j],gap=spec.closure?-.0003:plan.oralJetAreas?Math.sqrt(plan.oralJetAreas[spec.constriction]*1e-6/Math.PI):spec.constriction==='blade'?.0004:.0003;
   for(let i=e.c.starts[j];i<e.c.starts[j]+e.c.counts[j];i++){const local=plan.vowelAreas?e.radiusTargets[i]:base;e.radiusTargets[i]=local-(local-gap)*Math.exp(-.5*((e.c.x[i]-center)/.003)**2);}
  }
  if(lipClosure||(spec.constriction==='lip'&&plan.oralJetAreas)){const e=this.parts.lips,center=e.c.x[e.n-1],base=this.targets[li][lip],gap=lipClosure?-.0003:Math.sqrt(plan.oralJetAreas.lip*1e-6/Math.PI);for(let i=e.c.starts[lip];i<e.n;i++){const local=plan.vowelAreas?e.radiusTargets[i]:base;e.radiusTargets[i]=local-(local-gap)*Math.exp(-.5*((e.c.x[i]-center)/.002)**2);}}
 }
 applyAreaProfile(areas){
  if(!Array.isArray(areas)||areas.length<2||areas.length>128||areas.some(a=>!Number.isFinite(a)||a<=0||a>20))throw Error('Vowel areas must be 2–128 positive cross-sections, at most 20 cm².');
  const b=this.parts.body??this.system.engines[this.system.c.source],m=this.parts.mouth,l=this.parts.lips,begin=b.c.starts[b.c.cordIndex]+b.c.counts[b.c.cordIndex],origin=b.c.x[begin]-b.c.dx[begin]/2,bodyLength=b.c.length-origin,total=bodyLength+m.c.length+l.c.length;
  for(const [e,start,offset] of [[b,begin,-origin],[m,0,bodyLength],[l,0,bodyLength+m.c.length]])for(let i=start;i<e.n;i++){
   const position=Math.max(0,Math.min(areas.length-1,(e.c.x[i]+offset)/total*areas.length-.5)),k=Math.min(areas.length-2,Math.floor(position)),f=position-k,area=areas[k]*(1-f)+areas[k+1]*f;
   e.radiusTargets[i]=Math.sqrt(area*1e-4/Math.PI);
  }
 }
 snapshot(){return {text:this.plan.text,phone:this.phone,time:this.time,seconds:this.plan.seconds,finished:this.finished};}
}
