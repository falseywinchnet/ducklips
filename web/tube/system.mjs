import {compile,validateDesign,TubeEngine,RHO,SOUND,C} from './engine.mjs';
const bounded=(x,a,b)=>Number.isFinite(x)&&x>=a&&x<=b;
export function validateSystem(d){
 if(d?.version!==4)return validateDesign(d);
 if(d.system!=='connected-speech'||!Array.isArray(d.tubes)||d.tubes.length<2||d.tubes.length>8||!Array.isArray(d.connections)||d.connections.length>16)throw Error('Invalid connected tube system.');
 if(!bounded(d.pressure,0,1600)||!bounded(d.wander,0,.1)||!bounded(d.elasticityScale,1,4)||!bounded(d.elasticityVariation,0,.1))throw Error('Invalid system pressure or elasticity.');
 const ids=new Set();for(const t of d.tubes){if(typeof t.id!=='string'||!/^[a-z][a-z0-9-]{0,23}$/.test(t.id)||ids.has(t.id)||typeof t.name!=='string'||t.name.length>60||!['open','closed','connected'].includes(t.end))throw Error('Invalid tube identity or termination.');ids.add(t.id);validateDesign(t.design);}
 if(!ids.has(d.source)||!ids.has(d.output)||d.source===d.output)throw Error('Invalid source or output tube.');
 for(const e of d.connections)for(const p of [e.a,e.b]){
  const t=d.tubes.find(t=>t.id===p?.tube);if(!t)throw Error('A connection refers to a missing tube.');
  if(p.end!==undefined){if(!['inlet','outlet'].includes(p.end))throw Error('Invalid tube end.');}
  else if(!Number.isInteger(p.section)||p.section<0||p.section>=t.design.bands.length||!bounded(p.position,0,1))throw Error('Invalid side connection.');
 }
 const ends=new Set(),adj=new Map(d.tubes.map(t=>[t.id,[]]));
 for(const e of d.connections){
  if(e.a.tube===e.b.tube)throw Error('A tube cannot connect to itself.');
  adj.get(e.a.tube).push(e.b.tube);adj.get(e.b.tube).push(e.a.tube);
  for(const p of [e.a,e.b])if(p.end){const key=p.tube+':'+p.end;if(ends.has(key))throw Error('An end can have only one connection; branch from a section.');ends.add(key);}
 }
 if(ends.has(d.source+':inlet'))throw Error('The source inlet must remain a pressure reservoir.');
 for(const t of d.tubes){
  if(t.id!==d.source&&!ends.has(t.id+':inlet'))throw Error('Each downstream inlet needs a connection.');
  if(ends.has(t.id+':outlet')!==(t.end==='connected'))throw Error('Outlet termination does not match its connection.');
 }
 if(d.tubes.find(t=>t.id===d.output).end!=='open')throw Error('The sampled outlet must be open.');
 const seen=new Set(),visit=id=>{if(seen.has(id))return;seen.add(id);adj.get(id).forEach(visit);};visit(d.source);
 if(seen.size!==d.tubes.length)throw Error('Every tube must connect to the source.');
 return d;
}
export function compileSystem(input){
 const d=validateSystem(input);if(d.version!==4)return compile(d);
 const tubes=d.tubes.map(t=>({...t,compiled:compile({...t.design,pressure:t.id===d.source?d.pressure:0,wander:t.id===d.source?d.wander:0,elasticityScale:d.elasticityScale,elasticityVariation:d.elasticityVariation})}));
 // One system elasticity control multiplies every section's local set value.
 const endpoint=p=>{const tube=tubes.findIndex(t=>t.id===p.tube),c=tubes[tube].compiled;return {tube,end:p.end,index:p.end==='inlet'?0:p.end==='outlet'?c.n-1:c.starts[p.section]+Math.round(p.position*(c.counts[p.section]-1))};};
 return {kind:'network',n:tubes.reduce((sum,t)=>sum+t.compiled.n,0),pressure:d.pressure,tubes,source:tubes.findIndex(t=>t.id===d.source),output:tubes.findIndex(t=>t.id===d.output),connections:d.connections.map(e=>({a:endpoint(e.a),b:endpoint(e.b)})),topology:JSON.stringify({ids:tubes.map(t=>[t.id,t.compiled.n,t.end]),connections:d.connections,source:d.source,output:d.output})};
}
// A causal bidirectional pressure/tissue connection between complete sheets.
// Each end enters its local collocated junction as another physical port.
class Coupling {
 constructor(a,b,rate){
  this.a=a;this.b=b;this.rate=rate;this.pos=0;this.mask=255;
  for(const key of ['airAB','airBA','solidAB','solidBA'])this[key]=new Float64Array(256);
  this.pa={index:a.index,airY:0,airIn:0,airOut:0,solidZ:0,solidIn:0,solidOut:0};this.pb={...this.pa,index:b.index};a.engine.ports.push(this.pa);b.engine.ports.push(this.pb);this.configure();
 }
 configure(){
  const {a,b,rate}=this,ca=a.engine.c,cb=b.engine.c,i=a.index,j=b.index;
  const distance=(ca.dx[i]+cb.dx[j])/2;
  this.airDelay=Math.max(1,distance/SOUND*rate);
  const sa=ca.sheet,sb=cb.sheet;
  const speed=Math.sqrt((sa.axial[Math.min(i,ca.n-2)]*ca.dx[i]**2/sa.mass[i]+sb.axial[Math.min(j,cb.n-2)]*cb.dx[j]**2/sb.mass[j])/2);
  this.solidDelay=Math.max(1,distance/Math.max(.1,speed)*rate);
  this.taps=sa.fir[i].map((v,k)=>Math.sqrt(v*sb.fir[j][k]));
  const capacity=2**Math.ceil(Math.log2(Math.max(this.airDelay,this.solidDelay/Math.sqrt(.45))+5));
  if(capacity>16384)throw Error('Connection exceeds available tissue history.');
  if(capacity>this.mask+1){for(const key of ['airAB','airBA','solidAB','solidBA']){const old=this[key],next=new Float64Array(capacity);for(let lag=1;lag<=old.length;lag++)next[(capacity-lag)&(capacity-1)]=old[(this.pos-lag)&this.mask];this[key]=next;}this.pos=0;this.mask=capacity-1;}
 }
 filter(history,delay,taps){const whole=Math.floor(delay),f=delay-whole,p=this.pos-whole,m=this.mask;return (1-f)*(taps[0]*history[p&m]+taps[1]*history[(p-1)&m]+taps[2]*history[(p-2)&m])+f*(taps[0]*history[(p-1)&m]+taps[1]*history[(p-2)&m]+taps[2]*history[(p-3)&m]);}
 prepare(solidOnly=false){
  const {a,b,pa,pb}=this,ca=a.engine.c,cb=b.engine.c;
  if(!solidOnly){const radius=Math.max(.000025,Math.min(a.engine.restRadii[a.index]+a.engine.displacement[a.index],b.engine.restRadii[b.index]+b.engine.displacement[b.index]));
  this.airScale=Math.sqrt(Math.PI/(C*RHO*SOUND))*radius;pa.airY=pb.airY=this.airScale*this.airScale;}else pa.airY=pb.airY=0;
  const fa=a.engine.variation.factors[ca.bandAt[a.index]],fb=b.engine.variation.factors[cb.bandAt[b.index]],factor=Math.sqrt(fa*fb);
  const z=Math.sqrt(ca.sheet.impedance[a.index]*cb.sheet.impedance[b.index])*Math.sqrt(factor);this.solidScale=Math.sqrt(z);pa.solidZ=pb.solidZ=z;
  if(!solidOnly){pa.airIn=this.filter(this.airBA,this.airDelay,AIR_TAPS)/this.airScale;pb.airIn=this.filter(this.airAB,this.airDelay,AIR_TAPS)/this.airScale;}
  const delay=Math.max(1,this.solidDelay/Math.sqrt(factor));pa.solidIn=this.filter(this.solidBA,delay,this.taps)/this.solidScale;pb.solidIn=this.filter(this.solidAB,delay,this.taps)/this.solidScale;
 }
 send(solidOnly=false){if(!solidOnly){this.airAB[this.pos]=this.pa.airOut*this.airScale;this.airBA[this.pos]=this.pb.airOut*this.airScale;}this.solidAB[this.pos]=this.pa.solidOut*this.solidScale;this.solidBA[this.pos]=this.pb.solidOut*this.solidScale;this.pos=(this.pos+1)&this.mask;}
}
const AIR_TAPS=[.985,.010,.003];
export class SpeechSystem {
 constructor(c,rate=24000,options={}){if(options.symmetry===false)throw Error('Network ports currently require ring-symmetric loading.');this.rate=rate;this.n=c.n;this.c=c;this.fault=null;this.engines=c.tubes.map(t=>new TubeEngine(t.compiled,rate,options));this.configureTerminations();this.links=c.connections.map(e=>new Coupling({...e.a,engine:this.engines[e.a.tube]},{...e.b,engine:this.engines[e.b.tube]},rate));}
 configureTerminations(){
  for(let i=0;i<this.engines.length;i++){const e=this.engines[i];e.inletMode=i===this.c.source?'reservoir':'closed';e.outletMode=this.c.tubes[i].end==='open'?'open':'closed';}
  for(const link of this.c.connections)for(const p of [link.a,link.b]){if(p.end==='inlet')this.engines[p.tube].inletMode='connected';if(p.end==='outlet')this.engines[p.tube].outletMode='connected';}
 }
 configure(c){if(c.topology!==this.c.topology)throw Error('A changed network topology needs a fresh engine.');this.c=c;for(let i=0;i<this.engines.length;i++)this.engines[i].configure(c.tubes[i].compiled);this.configureTerminations();for(const link of this.links)link.configure();}
 get pressureTarget(){return this.engines[this.c.source].pressureTarget;}
 set pressureTarget(value){this.engines[this.c.source].pressureTarget=value;}
 get drive(){return this.engines[this.c.source].drive;}
 get previousOutlet(){return this.engines[this.c.output].previousOutlet;}
 step(){if(this.fault)return 0;for(const link of this.links)link.prepare();let output=0;for(let i=0;i<this.engines.length;i++){const value=this.engines[i].step();if(i===this.c.output)output=value;if(this.engines[i].fault){this.fault=this.c.tubes[i].name+': '+this.engines[i].fault;return 0;}}for(const link of this.links)link.send();return output;}
 snapshot(){const tubes=this.engines.map(e=>e.snapshot());return {...tubes[this.c.source],tubes,modeCount:tubes.reduce((sum,s)=>sum+s.modeCount,0),fault:this.fault};}
}
export function makeEngine(c,rate=24000){return c.kind==='network'?new SpeechSystem(c,rate):new TubeEngine(c,rate);}
export function sameTopology(engine,c){return engine.c.kind===c.kind&&(c.kind==='network'?engine.c.topology===c.topology:engine.n===c.n);}
