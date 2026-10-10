import {belowMouthDesign} from './presets.mjs';
function segment(name,length,radius,material,mass,elasticity,support){return {name,role:'tract',length,radius,material,mass,elasticity,support,offset:0,anchor:.5,points:4};}
function tube(id,name,bands,end='closed'){return {id,name,end,design:{version:3,pressure:0,wander:0,monitor:.25,elasticityScale:1,elasticityVariation:.1,noiseSeed:[...id].reduce((n,c)=>(n*31+c.charCodeAt(0))>>>0,73021),bands,blocks:[]}};}
export function systemDesign(arrangement='balanced'){
 const body=belowMouthDesign();
 const mouth=tube('mouth','Mouth',[
  segment('01 · Rear mouth',.025,.0054,'soft',2,1,1),
  segment('02 · Tongue body',.030,.008,'soft',3.5,1.5,.2),
  segment('03 · Palate chamber',.030,.011,'firm',3,.8,2),
  segment('04 · Front tongue',.020,.0055,'soft',2.5,1.5,.3)
 ],'connected');
 const lips=tube('lips','Teeth and lips',[
  segment('01 · Teeth gap',.010,.004,'firm',4,1,4),
  segment('02 · Lip tissue',.012,.0035,'fold',.8,1.5,.3)
 ],'open');
 const sinus=tube('sinus','Sinus pocket',[
  segment('01 · Sinus neck',.012,.0018,'firm',1,1,2),
  segment('02 · Sinus chamber',.035,.008,'firm',2,1,3)
 ]);
 const tongue=tube('undertongue','Under-tongue pocket',[
  segment('01 · Under-tongue entry',.012,.0025,'soft',2,1,.3),
  segment('02 · Floor pocket',.025,.0055,'soft',3,1.5,.2)
 ]);
 if(arrangement==='open'){mouth.design.bands[1].radius=.011;mouth.design.bands[3].radius=.008;lips.design.bands[0].radius=.006;lips.design.bands[1].radius=.006;sinus.design.bands[0].radius=.001;lips.design.bands.unshift(segment('01 · Lip vestibule',.008,.008,'soft',1.5,1,.3));lips.design.bands.forEach((b,i)=>b.name=String(i+1).padStart(2,'0')+b.name.slice(2));}
 if(arrangement==='rounded'){mouth.design.bands[1].radius=.005;mouth.design.bands[2].radius=.010;mouth.design.bands[3].radius=.007;lips.design.bands[1].length=.020;lips.design.bands[1].radius=.0025;tongue.design.bands[0].radius=.0035;mouth.design.bands.splice(3,1);Object.assign(mouth.design.bands[2],{name:'03 · Soft front chamber',length:.050,material:'soft',mass:3.5,elasticity:1.5,support:.2});}
 if(!['balanced','open','rounded'].includes(arrangement))throw Error('Unknown tube arrangement.');
 return {version:4,system:'connected-speech',arrangement,pressure:body.pressure,wander:body.wander,elasticityScale:1,elasticityVariation:.1,monitor:.25,
  tubes:[{id:'body',name:'Below-mouth',end:'connected',design:body},mouth,lips,sinus,tongue],
  connections:[
   {a:{tube:'body',end:'outlet'},b:{tube:'mouth',end:'inlet'}},
   {a:{tube:'mouth',end:'outlet'},b:{tube:'lips',end:'inlet'}},
   {a:{tube:'mouth',section:0,position:.65},b:{tube:'sinus',end:'inlet'}},
   {a:{tube:'mouth',section:2,position:.35},b:{tube:'undertongue',end:'inlet'}}
  ],source:'body',output:'lips'};
}
