// Experimental bidirectional air graph; aerodynamic/contact splitting remains approximate. Material nodes project their
// volume velocity onto junctions spaced c/rate apart; their sheet modes remain
// unchanged. Area changes scatter normalized traveling waves without gain.
import {C,RHO,SOUND} from './engine.mjs';
export class OralAir {
 constructor(system){
  this.s=system;this.rate=system.rate;this.rootY=Math.sqrt(Math.PI/(C*RHO*SOUND));this.nodes=[];this.edges=[];this.walls=[];this.nasal=0;this.nasalTarget=0;this.lastFlow=0;this.glottalFlow=0;this.glottalMomentum=0;this.returnWave=0;this.random=873219;this.jetNoise=0;
  const node=()=>{const i=this.nodes.length;this.nodes.push({ports:[],walls:[],p:0,radiation:new Float64Array(5)});return i;};
  const stages=system.c.tubes.map((t,i)=>({id:t.id,e:system.engines[i]}));
  const part=id=>stages.find(t=>t.id===id).e;
  const mouth=part('mouth'),lips=part('lips'),body=system.engines[system.c.source];this.body=body;this.lips=lips;
  const main=Array.from({length:Math.max(3,Math.round((mouth.c.length+lips.c.length)*this.rate/SOUND))+1},node);this.main=main;
  const locations=[];for(const [e,offset] of [[mouth,0],[lips,mouth.c.length]])for(let i=0;i<e.n;i++)locations.push({e,i,x:offset+e.c.x[i]});
  const length=mouth.c.length+lips.c.length;
  this.buildPath(main,locations,length);this.edges.at(-1).aperture={e:lips,i:lips.n-1};
  for(const id of ['sinus','undertongue']){
   const tubeIndex=system.c.tubes.findIndex(t=>t.id===id),mouthIndex=system.c.tubes.findIndex(t=>t.id==='mouth'),connection=system.c.connections.find(link=>(link.a.tube===tubeIndex&&link.b.tube===mouthIndex)||(link.b.tube===tubeIndex&&link.a.tube===mouthIndex));
   if(!connection)throw Error('The oral branches must connect to the mouth.');const endpoint=connection.a.tube===mouthIndex?connection.a:connection.b,attach=mouth.c.x[endpoint.index]/length;
   const e=part(id),branch=[main[Math.max(1,Math.round(attach*(main.length-1)))]];
   for(let i=0;i<Math.max(2,Math.round(e.c.length*this.rate/SOUND));i++)branch.push(node());
   this.buildPath(branch,Array.from({length:e.n},(_,i)=>({e,i,x:e.c.x[i]})),e.c.length);
   if(id==='sinus'){this.nose=branch.at(-1);this.noseEdge=this.edges.at(-1);this.nodes[this.nose].open=true;}else this.nodes[branch.at(-1)].closed=true;
  }
  this.outlet=main.at(-1);this.nodes[this.outlet].open=true;
  this.slit=body.c.slotWidth[body.c.starts[body.c.cordIndex]]>0;this.meanPressure=new Float64Array(body.n);
  this.port={index:this.slit?body.c.starts[body.c.cordIndex]-1:body.c.starts[body.c.cordIndex]+body.c.counts[body.c.cordIndex]-1,airY:0,airIn:0,airOut:0,solidZ:0,solidIn:0,solidOut:0};body.ports.push(this.port);body.outletMode='open';
  this.oral=system.engines.filter(e=>e!==body);this.referenceY=Math.PI*.006**2/(C*RHO*SOUND);
  for(const e of this.oral){if(!e.native)throw Error('Articulated oral sheets require the compiled kernel.');}
 }
 buildPath(path,locations,length){
  for(let k=0;k<path.length-1;k++){
   const x=(k+.5)/(path.length-1)*length,wall=locations.reduce((a,b)=>Math.abs(a.x-x)<Math.abs(b.x-x)?a:b);
   const edge={a:path[k],b:path[k+1],wall,r:0,l:0,nextR:0,nextL:0,y:0,root:0,jet:0,flow:0};
   const index=this.edges.length;this.edges.push(edge);this.nodes[edge.a].ports.push({edge:index,side:0});this.nodes[edge.b].ports.push({edge:index,side:1});
  }
  for(const wall of locations){const k=Math.max(1,Math.min(path.length-1,Math.round(wall.x/length*(path.length-1))));this.nodes[path[k]].walls.push(wall);this.walls.push(wall);}
 }
 step(){
  const s=this.s,h=1/this.rate;if(s.fault)return 0;this.nasal+=(this.nasalTarget-this.nasal)/(this.rate*.02);
  for(const link of s.links){link.prepare(true);}
  for(const edge of this.edges){const {e,i}=edge.wall;let r=Math.max(.000025,e.restRadii[i]+e.displacement[i]);if(edge.aperture){const a=edge.aperture;r=Math.min(r,Math.max(.000025,a.e.restRadii[a.i]+a.e.displacement[a.i]));}edge.radius=r;edge.root=this.rootY*r;edge.y=edge.root*edge.root;}
  const source=this.body,fold=source.c.cordIndex,begin=source.c.starts[fold],end=begin+source.c.counts[fold];
  if(!this.slit)for(let i=begin;i<end;i++){const slot=source.c.slotWidth[i],r=Math.max(slot?.0000001:.000025,source.restRadii[i]+source.displacement[i]),area=(slot?2*slot*r:Math.PI*r*r)/C,u=this.glottalFlow/area;source.wallPressureOffset[i]=-(this.aeroScale??1)*Math.min(1.5*source.drive,.5*RHO*u*u);}
  const first=this.edges[0];this.port.airY=this.body.c.area[this.port.index]/(C*RHO*SOUND);this.port.airIn=this.returnWave;this.body.step();if(this.body.fault){s.fault=this.body.fault;return 0;}
  let gap=Infinity;const cord=this.body.c.cordIndex;for(let i=this.body.c.starts[cord];i<this.body.c.starts[cord]+this.body.c.counts[cord];i++)gap=Math.min(gap,Math.max(this.body.c.slotWidth[i]?.0000001:.000025,this.body.restRadii[i]+this.body.displacement[i]));
  this.random=(1664525*this.random+1013904223)>>>0;this.jetNoise=.65*this.jetNoise+.35*(this.random/2147483648-1);
  const slit=this.body.c.slotWidth[this.body.c.starts[cord]],area=(slit?2*slit*gap:Math.PI*gap*gap)/C,reynolds=RHO*Math.abs(this.glottalFlow)/area*2*gap/.0000181,turbulence=Math.max(0,Math.min(1,(reynolds-1200)/2000));
  const loss=RHO/(2*.7*.7*area*area)*(1+(this.breathy?.12:.015)*turbulence*this.jetNoise),inertance=RHO*this.body.c.counts[cord]*this.body.c.dx[this.body.c.starts[cord]]/area;
  const a=this.port.airOut,b=first.l/first.root,linear=inertance/h+1/this.port.airY+1/first.y,drive=2*(a-b)+this.glottalMomentum/h;
  this.glottalFlow=2*drive/(linear+Math.sqrt(linear*linear+4*loss*Math.abs(drive)));
  this.glottalMomentum=inertance*this.glottalFlow;
  if(this.slit){
   const sub=2*a-this.glottalFlow/this.port.airY,oral=2*b+this.glottalFlow/first.y;let separation=0,minimum=Infinity;
   for(let i=begin;i<end;i++){const height=Math.max(.0000001,source.restRadii[i]+source.displacement[i]),localArea=2*source.c.slotWidth[i]*height/C;minimum=Math.min(minimum,localArea);const opening=Math.max(0,Math.min(1,(localArea/minimum-1)/.4));separation=Math.max(separation,(this.smoothSeparation?opening*opening*(3-2*opening):+(localArea>1.2*minimum)));
    const speed=this.glottalFlow/localArea,attached=Math.max(oral-2*Math.abs(sub-oral),sub-.5*RHO*speed*speed),target=this.smoothSeparation?attached+(oral-attached)*separation:separation?oral:attached;
    this.meanPressure[i]+=(source.pressure[i]-this.meanPressure[i])/(this.rate*.02);source.wallPressureOffset[i]=target-this.meanPressure[i];
   }
  }
  this.returnWave=a-this.glottalFlow/this.port.airY;
  first.nextR=(b+this.glottalFlow/first.y)*first.root;
  for(const e of this.oral){e.tick++;if(e.variable)e.advanceVariation();e.native.prepareWalls(h);}
  let flow=0;
  for(let k=1;k<this.nodes.length;k++){
   const node=this.nodes[k];let sum=0,admit=0,wallFlow=0,wallAdmit=0;
   for(const p of node.ports){const e=this.edges[p.edge],incoming=p.side?e.r:e.l;sum+=e.root*incoming;admit+=e.y;}
   for(const w of node.walls){wallFlow+=w.e.sum[w.i];wallAdmit+=w.e.admit[w.i];}
   const open=node.open;
   // Both lips and nostrils radiate. Nasal coupling is controlled upstream
   // by the moving velum neck, not by closing the nostril termination.
   // The open-end reflection is a passive fractional-delay FIR. Its round
   // trip delay is the unflanged end correction (0.61 aperture radii each way).
   // Solve its direct feedthrough together with the wall/air junction.
   let radiation=0,history=0,direct=0,loadY=0;
   if(open){const edge=this.edges[node.ports[0].edge],r=edge.radius,delay=Math.min(3.99,1.22*r*this.rate/SOUND),whole=Math.floor(delay),fraction=delay-whole;
    if(whole===0)direct=-.995*(1-fraction);else history=-.995*(1-fraction)*node.radiation[whole-1];history-=.995*fraction*node.radiation[whole];
    loadY=admit;radiation=loadY*(1-direct)/(1+direct);history/=1+direct;
   }
   const pressure=(2*sum+2*loadY*history-wallFlow)/(admit+wallAdmit+radiation);node.p=pressure;
   if(open){const incoming=direct/(1+direct)*pressure+history;for(let j=4;j>0;j--)node.radiation[j]=node.radiation[j-1];node.radiation[0]=pressure-incoming;flow+=loadY*(pressure-2*incoming);}
   for(const w of node.walls)w.e.pressure[w.i]=pressure;
   for(const p of node.ports){const e=this.edges[p.edge],incoming=p.side?e.r:e.l,out=(e.root*pressure-incoming)*.997;if(p.side)e.nextL=out;else e.nextR=out;}

  }
  for(const e of this.oral){const contacts=e.native.finishWalls();if(contacts<0){s.fault=e.fault='Articulated sheet motion exceeded its domain.';return 0;}else e.contactCount+=contacts;}
  for(const edge of this.edges){
   const r=edge.radius;
   if(r<.0015){
    // The narrow aperture dissipates mean flow and reflects its obstruction
    // back into both cavities. A bounded fluctuation of that loss produces
    // jet noise; without pressure/flow this element is exactly silent.
    const area=Math.PI*r*r/C,z=1/edge.y;
    const a=edge.nextR/edge.root,b=edge.nextL/edge.root,drive=2*(a-b),linear=2*z;
    this.random=(1664525*this.random+1013904223)>>>0;edge.jet=.4*edge.jet+.6*(this.random/2147483648-1);
    const reynolds=RHO*Math.abs(edge.flow)/area*2*r/.0000181;
    const resistance=RHO/(2*.7*.7*area*area)*(1+.08*Math.max(0,Math.min(1,(reynolds-1200)/2000))*edge.jet);
    const q=2*drive/(linear+Math.sqrt(linear*linear+4*resistance*Math.abs(drive)));edge.flow=q;
    edge.r=(b+q*z)*edge.root;edge.l=(a-q*z)*edge.root;
   }else{edge.r=edge.nextR;edge.l=edge.nextL;}
  }
  for(const link of s.links)link.send(true);
  this.lips.previousOutlet=this.nodes[this.outlet].p;
  const output=(flow-this.lastFlow)/this.referenceY*this.rate/48000/10;this.lastFlow=flow;
  if(!Number.isFinite(output)||Math.abs(output)>20){s.fault='Oral radiation exceeded its domain.';return 0;}return output;
 }
}
