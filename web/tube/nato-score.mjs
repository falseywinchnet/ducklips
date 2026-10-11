// Reference-guided articulator timing, not sampled excitation or output audio.
// Initial Papa score follows the visible reduced first vowel and stressed second
// vowel in the preserved valeatory human reference. Other words remain untuned.
import {NATO_TARGETS} from './nato-targets.mjs';
import {planSpeech} from './articulation.mjs';
export const NATO=['alpha','bravo','charlie','delta','echo','foxtrot','golf','hotel','india','juliet','kilo','lima','mike','november','oscar','papa','quebec','romeo','sierra','tango','uniform','victor','whiskey','x-ray','yankee','zulu'];
export function referenceScore(word){
 const target=NATO_TARGETS[word];if(!target)throw Error('No NATO reference target: '+word);
 const plan=planSpeech(word);plan.coupledAir=true;plan.foldVolume=false;plan.localClosures=true;plan.honorOnsets=true;plan.effortScope='body';
 const mapping={CH:'t-sh',JH:'d-zh',NG:'ng',HH:'h'},events=[];
 for(const token of target.phones){const key=token.replace(/[012]/g,''),stress=token.endsWith('1');for(const event of planSpeech('/'+(mapping[key]??key.toLowerCase())+'/').events.filter(e=>e.phone!=='_'))events.push({...event,word,weight:(event.end-event.start)*(stress?1.25:1),stress});}
 const total=events.reduce((s,e)=>s+e.weight,0),duration=target.seconds-.23;let at=.08;
 for(const event of events){event.start=at;event.end=at+duration*event.weight/total;at=event.end;event.pressure=event.stress?1:.65;if(event.stress){event.effort=2.5;event.effortEnd=1.5;}delete event.weight;}
 events.push({phone:'_',word,start:at,end:target.seconds});plan.events=events;plan.seconds=target.seconds;plan.referenceTiming='Clip duration and working stress transcription; individual boundaries unverified.';
 if(word==='papa'){
  plan.events=[
   {phone:'p',start:0,end:.125,word,pressure:.3,effort:1.2},
   {phone:'ax',start:.125,end:.20,word,pressure:.6,effort:1.2},
   {phone:'p',start:.20,end:.245,word,effort:3.2},
   {phone:'aa',start:.245,end:.555,word,pressure:1,pressureEnd:.75,effort:3.2,effortEnd:1.1},
   {phone:'_',start:.555,end:.811,word}
  ];plan.seconds=.811;plan.referenceTiming='Visual onset/stress estimate from human Papa; not an exact alignment.';
 }
 return plan;
}
