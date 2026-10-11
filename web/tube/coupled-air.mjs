import {NativeSheets} from './native-core.mjs';
import {ApertureRadiation,RADIATION_TAPS} from './radiation.mjs';
import {OralAir} from './oral-air.mjs';
import {C,RHO,SOUND} from './engine.mjs';
// Sudden-expansion loss tends continuously to zero when neighboring areas
// match. This is a reduced quasi-steady loss, not a resolved turbulent jet.
export function expansionLoss(radius,neighborRadius){
 const ratio=Math.min(1,(radius/neighborRadius)**2);return (1-ratio)**2;
}
// Experimental single serial airway. Every wall has one pressure source;
// the fold is an inertial flow element between the two adjacent cavities.
export class CoupledAir extends OralAir {
 constructor(system,{batch=true}={}){
  super(system);const e=this.body,c=e.c;this.unified=true;
  e.ports.splice(e.ports.indexOf(this.port),1);this.port=null;
  this.begin=c.starts[c.cordIndex];this.end=this.begin+c.counts[c.cordIndex];
  if(!this.slit)throw Error('The coupled-air trial requires a slit cord.');
  const node=()=>{const i=this.nodes.length;this.nodes.push({ports:[],walls:[],p:0,radiation:new Float64Array(5)});return i;};
  const before=c.x[this.begin]-c.dx[this.begin]/2,after=c.x[this.end-1]+c.dx[this.end-1]/2;
  const upstream=Array.from({length:Math.max(2,Math.round(before*this.rate/SOUND))+1},node);
  this.inlet=upstream[0];this.upstream=upstream.at(-1);
  this.buildPath(upstream,Array.from({length:this.begin},(_,i)=>({e,i,x:c.x[i]})),before);
  const downstream=Array.from({length:Math.max(2,Math.round((c.length-after)*this.rate/SOUND))},node);downstream.push(this.main[0]);this.downstream=downstream[0];
  this.buildPath(downstream,Array.from({length:e.n-this.end},(_,j)=>({e,i:this.end+j,x:c.x[this.end+j]-after})),c.length-after);
  for(const edge of this.edges)edge.neighbors=Int32Array.from(new Set([...this.nodes[edge.a].ports,...this.nodes[edge.b].ports].map(p=>p.edge)));
  this.inletY=c.area[0]/(C*RHO*SOUND);this.foldArea=new Float64Array(this.end-this.begin);this.foldInertia=new Float64Array(this.foldArea.length);this.foldPressure=new Float64Array(this.foldArea.length);this.foldUp=new Float64Array(this.foldArea.length);this.foldOffset=new Float64Array(this.foldArea.length);this.foldSurface=new Float64Array(this.foldArea.length);this.maxVolumeError=0;
  this.glottalMomentum=0;this.glottalFlow=0;this.oral=system.engines;for(const part of this.oral)part.contactMode=true;this.batch=batch?new NativeSheets(this.oral):null;
  this.radiationModel=new ApertureRadiation(this.rate);for(const n of this.nodes)if(n.open){n.radiation=new Float64Array(RADIATION_TAPS);n.radiationCoefficients=new Float64Array(RADIATION_TAPS);}
 }
 buildPath(path,locations,length){
  const begin=this.edges.length;super.buildPath(path,locations,length);
  for(let k=0;k<path.length-1;k++){const edge=this.edges[begin+k];edge.apertures=locations.filter(w=>w.x>=k/(path.length-1)*length&&w.x<=(k+1)/(path.length-1)*length);edge.lossR1=edge.lossR2=edge.lossL1=edge.lossL2=0;}
 }
 lossVariation(flow,area,diameter,state,key,strength){
  if(!this.flowNoise)return 1;
  this.random=(1664525*this.random+1013904223)>>>0;
  state[key]=.4*state[key]+.6*(this.random/2147483648-1);
  const reynolds=RHO*Math.abs(flow)*diameter/(area*.0000181);
  // Fluctuating jet loss remains strictly positive. It cannot sound without
  // pressure-driven flow, and it exchanges the same flow with both cavities.
  return 1+strength*Math.max(0,Math.min(1,(reynolds-1200)/2000))*state[key];
 }
 step(){
  const s=this.s,h=1/this.rate;if(s.fault)return 0;
  for(const link of s.links)link.prepare(true);
  for(const edge of this.edges){const {e,i}=edge.wall;let r=Math.max(.000025,e.restRadii[i]+e.displacement[i]);if(edge.aperture){const a=edge.aperture;r=Math.min(r,Math.max(.000025,a.e.restRadii[a.i]+a.e.displacement[a.i]));}for(const w of edge.apertures)r=Math.min(r,Math.max(.000025,w.e.restRadii[w.i]+w.e.displacement[w.i]));edge.radius=r;edge.root=this.rootY*r;edge.y=edge.root*edge.root;}
  if(this.batch&&!this.batch.current())this.batch=new NativeSheets(this.oral);
  for(const e of this.oral){e.advanceControl();if(!this.batch)e.native.prepareWalls(h);}if(this.batch)this.batch.prepare(h);
  // First gather the incident waves and free wall motion. The valve exchanges
  // equal volume flow at both ends within these same junction equations.
  for(let k=0;k<this.nodes.length;k++){
   const n=this.nodes[k];let sum=0,admit=0,wallFlow=0,wallAdmit=0;
   for(const p of n.ports){const e=this.edges[p.edge];sum+=e.root*(p.side?e.r:e.l);admit+=e.y;}
   for(const w of n.walls){wallFlow+=w.e.sum[w.i];wallAdmit+=w.e.admit[w.i];}
   n.rhs=2*sum-wallFlow;n.admit=admit+wallAdmit;n.airY=admit;
   if(k===this.inlet){n.rhs+=this.inletY*this.body.drive;n.admit+=this.inletY;}
  }
  const upper=this.nodes[this.upstream],lower=this.nodes[this.downstream],body=this.body;
  let inertance=0,area=Infinity,viscous=0;
  for(let i=this.begin;i<this.end;i++){
   const k=i-this.begin,height=Math.max(.0000001,body.restRadii[i]+body.displacement[i]),a=2*body.c.slotWidth[i]*height/C;
   this.foldArea[k]=a;this.foldInertia[k]=RHO*body.c.dx[i]/a;inertance+=this.foldInertia[k];area=Math.min(area,a);
   // Laminar slit resistance, integrated along the actual opening.
   viscous+=12*.0000181*body.c.dx[i]/(a*(2*height)**2);
  }
  const p0=upper.rhs/upper.admit,p1=lower.rhs/lower.admit,oldMomentum=this.glottalMomentum;
  const resistance=RHO/(2*.7*.7*area*area)*this.lossVariation(this.glottalFlow,area,area*C/body.c.slotWidth[this.begin],this,'jetNoise',this.breathy?.12:.025),linear=inertance/h+viscous+1/upper.admit+1/lower.admit,drive=p0-p1+oldMomentum/h;
  const q=2*drive/(linear+Math.sqrt(linear*linear+4*resistance*Math.abs(drive)));this.glottalFlow=q;this.glottalMomentum=inertance*q;
  const sub=(upper.rhs-q)/upper.admit,oral=(lower.rhs+q)/lower.admit;
  const acceleration=(this.glottalMomentum-oldMomentum)/h;let accumulated=0,minimum=Infinity,separated=false,loss=0;
  for(let i=this.begin;i<this.end;i++){
   const k=i-this.begin,a=this.foldArea[k],height=a*C/(2*body.c.slotWidth[i]);minimum=Math.min(minimum,a);if(a>1.2*minimum)separated=true;
   accumulated+=this.foldInertia[k]/2;loss+=6*.0000181*body.c.dx[i]/(a*(2*height)**2);
   // One absolute force pressure, with the same fluid inertia and viscosity
   // used in the end-to-end flow equation. No mean-subtracted second pressure.
   const pressure=separated?oral:sub-acceleration*accumulated/inertance-loss*q-.5*RHO*(q/a)**2;
   body.pressure[i]=this.foldPressure[k]=pressure;this.foldUp[k]=separated?0:1-accumulated/inertance;
   accumulated+=this.foldInertia[k]/2;loss+=6*.0000181*body.c.dx[i]/(a*(2*height)**2);
  }
  if(this.foldVolume)this.exchangeFoldVolume(upper,lower,h,inertance,viscous,resistance,oldMomentum,sub,oral);else{upper.rhs-=q;lower.rhs+=q;}
  let flow=0;
  for(const n of this.nodes){
   let radiation=0,history=0,direct=0,loadY=0;
   if(n.open){const edge=this.edges[n.ports[0].edge];
    if(this.apertureRadiation){this.radiationModel.gather(n,edge.radius);direct=n.radiationDirect;history=n.radiationHistory;}
    else{const delay=Math.min(3.99,1.22*edge.radius*this.rate/SOUND),whole=Math.floor(delay),fraction=delay-whole;if(whole===0)direct=-.995*(1-fraction);else history=-.995*(1-fraction)*n.radiation[whole-1];history-=.995*fraction*n.radiation[whole];}
    loadY=n.airY;radiation=loadY*(1-direct)/(1+direct);history/=1+direct;
   }
   // A prescribed reservoir pressure is an explicit alternative to the
   // matched acoustic source impedance, which also drops steady pressure.
   const pressure=this.pressureBoundary&&n===this.nodes[this.inlet]?body.drive:(n.rhs+2*loadY*history)/(n.admit+radiation);n.p=pressure;
   if(n.open){const incoming=direct/(1+direct)*pressure+history;n.radiation.copyWithin(1,0,-1);n.radiation[0]=pressure-incoming;n.outletFlow=loadY*(pressure-2*incoming);flow+=n.outletFlow;}
   for(const w of n.walls)w.e.pressure[w.i]=pressure;
   for(const p of n.ports){const e=this.edges[p.edge],incoming=p.side?e.r:e.l,out=(e.root*pressure-incoming)*(this.airLossFIR?1:.997);if(p.side)e.nextL=out;else e.nextR=out;}
  }
  if(this.batch){if(this.batch.finish()<0){s.fault='Coupled sheet motion exceeded its domain.';return 0;}}else for(const e of this.oral){const contacts=e.native.finishWalls();if(contacts<0){s.fault='Coupled sheet motion exceeded its domain.';return 0;}e.contactCount+=contacts;}
  if(this.foldVolume){let actual=0;for(let i=this.begin;i<this.end;i++)actual+=this.foldSurface[i-this.begin]*body.force[i];this.volumeError=this.foldInflow-this.foldOutflow-actual;this.maxVolumeError=Math.max(this.maxVolumeError,Math.abs(this.volumeError));}
  for(const edge of this.edges){
   const r=edge.radius;let jetLoss=1;
   if(this.geometricJetLoss){let neighbor=r;for(const index of edge.neighbors)neighbor=Math.max(neighbor,this.edges[index].radius);jetLoss=expansionLoss(r,neighbor);}
   if(this.geometricJetLoss||r<.0015){const a=edge.nextR/edge.root,b=edge.nextL/edge.root,z=1/edge.y,area=Math.PI*r*r/C,resistance=jetLoss*RHO/(2*.7*.7*area*area)*this.lossVariation(edge.flow,area,2*r,edge,'jet',.12),drive=2*(a-b),linear=2*z,q=2*drive/(linear+Math.sqrt(linear*linear+4*resistance*Math.abs(drive)));edge.flow=q;edge.r=(b+q*z)*edge.root;edge.l=(a-q*z)*edge.root;}
   else{edge.r=edge.nextR;edge.l=edge.nextL;}
   if(this.airLossFIR){const r=edge.r,l=edge.l;edge.r=.985*r+.012*edge.lossR1+.003*edge.lossR2;edge.l=.985*l+.012*edge.lossL1+.003*edge.lossL2;edge.lossR2=edge.lossR1;edge.lossR1=r;edge.lossL2=edge.lossL1;edge.lossL1=l;}
  }
  for(const link of s.links)link.send(true);
  this.lips.previousOutlet=this.nodes[this.outlet].p;const output=(flow-this.lastFlow)/this.referenceY*this.rate/48000/10;this.lastFlow=flow;
  if(!Number.isFinite(output)||Math.abs(output)>20){s.fault='Coupled radiation exceeded its domain.';return 0;}return output;
 }
 exchangeFoldVolume(up,down,h,L,R,loss,momentum,predictUp,predictDown){
  const e=this.body,count=this.end-this.begin;let mask=0,pressureUp=predictUp,pressureDown=predictDown,q=this.glottalFlow;
  for(let k=0;k<count;k++){
   const i=this.begin+k,w=this.foldUp[k];this.foldOffset[k]=this.foldPressure[k]-w*predictUp-(1-w)*predictDown;
   this.foldSurface[k]=2*e.c.slotWidth[i]*e.c.dx[i]/C;
  }
  // The reciprocal two-port wall mobility is a 2x2 junction, independent of
  // the number of sheet modes. Contact only removes inward wall mobility.
  for(let pass=0;pass<=count;pass++){
   let m00=up.admit,m01=0,m11=down.admit,r0=up.rhs,r1=down.rhs;
   for(let k=0;k<count;k++){
    const i=this.begin+k,w=this.foldUp[k],v=1-w,closed=mask&(1<<k),Y=closed?0:e.admit[i],free=closed?this.foldSurface[k]*e.velocity[i]/2:e.sum[i]+Y*this.foldOffset[k];
    m00+=w*w*Y;m01+=w*v*Y;m11+=v*v*Y;r0-=w*free;r1-=v*free;
   }
   const determinant=m00*m11-m01*m01,p0=(r0*m11-r1*m01)/determinant,p1=(r1*m00-r0*m01)/determinant;
   const z=(m00+m11+2*m01)/determinant,linear=L/h+R+z,drive=p0-p1+momentum/h;
   q=2*drive/(linear+Math.sqrt(linear*linear+4*loss*Math.abs(drive)));
   pressureUp=p0-q*(m11+m01)/determinant;pressureDown=p1+q*(m00+m01)/determinant;
   let nextMask=0;
   for(let k=0;k<count;k++){
    const i=this.begin+k,w=this.foldUp[k],p=w*pressureUp+(1-w)*pressureDown+this.foldOffset[k];
    if(e.restRadii[i]+e.displacement[i]<=.0000001+1e-15&&e.sum[i]+e.admit[i]*p<this.foldSurface[k]*e.velocity[i]/2)nextMask|=1<<k;
   }
   if(nextMask===mask)break;if(pass===count)throw Error('Fold contact junction did not settle.');mask=nextMask;
  }
  let displacedUp=0,displacedDown=0;
  for(let k=0;k<count;k++){
   const i=this.begin+k,w=this.foldUp[k],p=w*pressureUp+(1-w)*pressureDown+this.foldOffset[k];e.pressure[i]=this.foldPressure[k]=p;
   const flow=mask&(1<<k)?this.foldSurface[k]*e.velocity[i]/2:e.sum[i]+e.admit[i]*p;displacedUp+=w*flow;displacedDown+=(1-w)*flow;
  }
  this.foldInflow=q+displacedUp;this.foldOutflow=q-displacedDown;this.glottalFlow=q;this.glottalMomentum=L*q;
  up.rhs-=this.foldInflow;down.rhs+=this.foldOutflow;
 }

}
