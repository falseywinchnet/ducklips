// Connected periodic sheet + pressure and tissue transport. No voice source.
import {compileSheet,CIRCUMFERENCE,modeShape,receptance} from './sheet.mjs';
export {eigen,modeShape} from './sheet.mjs';
export const RHO=1.204,SOUND=343,BANDS=6,P=8,N=BANDS*P,C=CIRCUMFERENCE;
export const MATERIALS={
 soft:{name:'Soft tissue',density:1040,modulus:8000,thickness:.002,loss:.12,fir:[.74,.17,.06],color:'#72c4b7'},
 fold:{name:'Elastic band',density:1040,modulus:22000,thickness:.001,loss:.055,fir:[.85,.10,.03],color:'#eab97a'},
 firm:{name:'Firm tissue',density:1100,modulus:65000,thickness:.002,loss:.16,fir:[.61,.23,.10],color:'#959dd9'}
};
export const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
export function initialDesign(){return {version:2,pressure:650,wander:.012,monitor:.25,blocks:[],bands:[
 {name:'01 · Return chamber',length:.050,radius:.009,material:'soft'},
 {name:'02 · Vocal band',length:.016,radius:.00032,material:'fold'},
 {name:'03 · Lower tract',length:.035,radius:.007,material:'soft'},
 {name:'04 · Oral body',length:.045,radius:.011,material:'soft'},
 {name:'05 · Front cavity',length:.035,radius:.008,material:'firm'},
 {name:'06 · Outlet',length:.025,radius:.006,material:'soft'}
 ].map(b=>({...b,mass:1,elasticity:1,offset:0,anchor:.5,support:1}))};}
export function validateDesign(d){
 if(!d||![1,2].includes(d.version)||!Array.isArray(d.bands)||d.bands.length!==BANDS||!Array.isArray(d.blocks)||d.blocks.length>128)throw Error('Expected a six-section tube design.');
 const finite=(v,a,b)=>Number.isFinite(v)&&v>=a&&v<=b;
 if(!finite(d.pressure,0,1600)||!finite(d.wander,0,.1))throw Error('Invalid pressure.');
 if(d.version===1){
  for(const b of d.bands)if(!finite(b.upper,-.02,.02)||!finite(b.lower,-.02,.02)||!finite(b.upperAnchor??.5,0,1)||!finite(b.lowerAnchor??.5,0,1))throw Error('Invalid legacy support.');
  d={...d,version:2,bands:d.bands.map(({upper,lower,upperAnchor=.5,lowerAnchor=.5,...b})=>({...b,offset:(upper+lower)/2,anchor:(upperAnchor+lowerAnchor)/2}))};
 }
 for(const b of d.bands)if(typeof b.name!=='string'||b.name.length>80||!finite(b.anchor,0,1)||!MATERIALS[b.material]||!finite(b.length,.008,.1)||!finite(b.radius,.00005,.025)||!finite(b.mass,.2,5)||!finite(b.elasticity,.2,5)||!finite(b.offset,-.02,.02)||!finite(b.support,.2,5))throw Error('Invalid section properties.');
 for(const b of d.blocks)if(!finite(b.x,0,1)||!finite(b.volume,.01,8)||!finite(b.width,.02,.5))throw Error('Invalid volume block.');return d;
}
export function compile(input){
 const d=validateDesign(input),x=[],dx=[],r=[],area=[];let start=0;
 const length=d.bands.reduce((s,b)=>s+b.length,0);
 for(const band of d.bands){for(let j=0;j<P;j++){dx.push(band.length/P);x.push(start+(j+.5)*band.length/P);r.push(Math.max(.000025,band.radius+band.offset));}start+=band.length;}
 const added=new Float64Array(N);
 for(const block of d.blocks){const weights=x.map((v,i)=>Math.max(0,1-Math.abs(v/length-block.x)/(block.width/2))*dx[i]);let sum=weights.reduce((s,v)=>s+v,0);
  if(sum===0){let nearest=0;for(let i=1;i<N;i++)if(Math.abs(x[i]/length-block.x)<Math.abs(x[nearest]/length-block.x))nearest=i;weights[nearest]=1;sum=1;}
  for(let i=0;i<N;i++)added[i]+=block.volume*1e-6*weights[i]/sum/dx[i];
 }
 for(let i=0;i<N;i++){area[i]=Math.PI*r[i]*r[i]+added[i];r[i]=Math.sqrt(area[i]/Math.PI);}
 const c={length,x,dx,r,upper:r,lower:r,area,pressure:d.pressure,wander:d.wander};c.sheet=compileSheet(c,d.bands,MATERIALS,P);return c;
}
class Link {
  constructor(a,b,delay,taps,weight=1){this.a=a;this.b=b;this.i=a%N;this.j=b%N;this.weight=weight;this.self=a===b;this.angular=false;this.z=0;this.y=0;this.scale=0;this.admittanceScale=0;this.ab=new Float64Array(16);this.ba=new Float64Array(16);this.mask=15;this.pos=0;this.ia=0;this.ib=0;this.setFilter(delay,taps);}
  setFilter(delay,taps){
    this.delay=clamp(delay,1,4090);this.taps=taps;
    const size=2**Math.ceil(Math.log2(this.delay+5));
    if(size>this.ab.length){const ab=new Float64Array(size),ba=new Float64Array(size);for(let i=1;i<=this.ab.length;i++){ab[(size-i)&(size-1)]=this.ab[(this.pos-i+this.ab.length)&this.mask];ba[(size-i)&(size-1)]=this.ba[(this.pos-i+this.ba.length)&this.mask];}this.ab=ab;this.ba=ba;this.mask=size-1;this.pos=0;}this.d=Math.floor(this.delay);const f=this.delay-this.d;
    // Fuse fractional propagation delay and three-tap attenuation into one FIR.
    this.h0=(1-f)*taps[0];this.h1=f*taps[0]+(1-f)*taps[1];this.h2=f*taps[1]+(1-f)*taps[2];this.h3=f*taps[2];
  }
  receive(){const i=(this.pos-this.d)&this.mask,j=(i-1)&this.mask,k=(i-2)&this.mask,l=(i-3)&this.mask;
    this.ia=this.h0*this.ba[i]+this.h1*this.ba[j]+this.h2*this.ba[k]+this.h3*this.ba[l];
    this.ib=this.self?this.ia:this.h0*this.ab[i]+this.h1*this.ab[j]+this.h2*this.ab[k]+this.h3*this.ab[l];
  }
  send(a,b){this.ab[this.pos]=a;this.ba[this.pos]=b;this.pos=(this.pos+1)&this.mask;}
}

// Contiguous delay histories and connection parameters keep the transport
// network in one representation instead of hundreds of independently indexed
// objects. Link objects remain the editable topology; this is its audio state.
class LinkBank {
 constructor(links){
  this.links=links;this.length=links.length;
  for(const name of ['a','b','i','j','offset','mask','pos','delay'])this[name]=new Int32Array(this.length);
  for(const name of ['h0','h1','h2','h3','gain','scale','ia','ib'])this[name]=new Float64Array(this.length);
  const size=links.reduce((sum,e)=>sum+e.ab.length,0);this.ab=new Float64Array(size);this.ba=new Float64Array(size);let offset=0;
  for(let k=0;k<links.length;k++){
   const e=links[k];for(const key of ['a','b','i','j','mask','pos','h0','h1','h2','h3'])this[key][k]=e[key];
   this.delay[k]=e.d;this.offset[k]=offset;this.gain[k]=e.z||e.admittanceScale;
   this.ab.set(e.ab,offset);this.ba.set(e.ba,offset);
   e.ab=this.ab.subarray(offset,offset+e.ab.length);e.ba=this.ba.subarray(offset,offset+e.ba.length);offset+=e.ab.length;
  }
 }
 sync(){for(let k=0;k<this.length;k++)this.links[k].pos=this.pos[k];}
 receiveAir(sum,admit,radii,areas,symmetry){
  const {a,b,i,j,offset,mask,pos,delay,h0,h1,h2,h3,ab,ba,gain,scale,ia,ib}=this;
  for(let k=0;k<this.length;k++){
   const o=offset[k],m=mask[k],p=(pos[k]-delay[k])&m,p0=o+p,p1=o+((p-1)&m),p2=o+((p-2)&m),p3=o+((p-3)&m);
   const s=gain[k]*(symmetry?Math.min(radii[i[k]],radii[j[k]]):Math.sqrt(Math.min(areas[i[k]],areas[j[k]])/Math.PI)),y=s*s;
   const va=(h0[k]*ba[p0]+h1[k]*ba[p1]+h2[k]*ba[p2]+h3[k]*ba[p3])/s;
   const vb=a[k]===b[k]?va:(h0[k]*ab[p0]+h1[k]*ab[p1]+h2[k]*ab[p2]+h3[k]*ab[p3])/s;
   scale[k]=s;ia[k]=va;ib[k]=vb;sum[a[k]]+=y*va;sum[b[k]]+=y*vb;admit[a[k]]+=y;admit[b[k]]+=y;
  }
 }
 receiveSolid(sum){
  const {a,b,offset,mask,pos,delay,h0,h1,h2,h3,ab,ba,gain,ia,ib}=this;
  for(let k=0;k<this.length;k++){
   const o=offset[k],m=mask[k],p=(pos[k]-delay[k])&m,p0=o+p,p1=o+((p-1)&m),p2=o+((p-2)&m),p3=o+((p-3)&m);
   const va=h0[k]*ba[p0]+h1[k]*ba[p1]+h2[k]*ba[p2]+h3[k]*ba[p3];
   const vb=a[k]===b[k]?va:h0[k]*ab[p0]+h1[k]*ab[p1]+h2[k]*ab[p2]+h3[k]*ab[p3];
   ia[k]=va;ib[k]=vb;sum[a[k]]+=gain[k]*va;sum[b[k]]+=gain[k]*vb;
  }
 }
 send(field,air){
  const {a,b,offset,mask,pos,ab,ba,scale,ia,ib}=this;
  for(let k=0;k<this.length;k++){
   const at=offset[k]+pos[k],s=air?scale[k]:1;
   ab[at]=(field[a[k]]-ia[k])*s;ba[at]=(field[b[k]]-ib[k])*s;pos[k]=(pos[k]+1)&mask[k];
  }
 }
}

export class TubeEngine {
 constructor(compiled,rate=24000,{symmetry=true}={}){
  this.rate=rate;this.tick=0;this.random=0x513fe2;this.wander=0;this.drive=0;this.fault=null;this.contactCount=0;
  // Ring controls and uniform inlet preserve rotation symmetry exactly. This
  // quotient stores one of eight identical sectors. The full-sheet path is
  // used for nonuniform loads/initial states and checked against the quotient.
  this.symmetry=symmetry;this.sectors=symmetry?1:C;this.size=N*this.sectors;
  for(const key of ['pressure','sum','admit','displacement','velocity','force','solidIncoming','solidAdmit','radii','dv'])this[key]=new Float64Array(this.size);
  this.areas=new Float64Array(N);this.air=[];this.solid=[];this.previousOutlet=0;this.configure(compiled);
 }
 configure(c){
  this.airBank?.sync();this.solidBank?.sync();this.c=c;this.pressureTarget=c.pressure;const sheet=c.sheet,h=1/this.rate;
  this.channels=sheet.channels.slice(0,this.symmetry?1:C).map(ch=>({...ch,kernel:receptance(ch,N,this.rate),nodalQ:new Float64Array(N),nodalV:new Float64Array(N),input:new Float64Array(N)}));
  if(!this.air.length){
   const next=t=>(t+1)%this.sectors,prev=t=>(t+this.sectors-1)%this.sectors;
   for(let t=0;t<this.sectors;t++){
    for(let i=0;i<N-1;i++){
     this.air.push(new Link(t*N+i,t*N+i+1,1,[.985,.010,.003],.36));
     for(const neighbour of this.symmetry?[next(t)]:[next(t),prev(t)]){const link=new Link(t*N+i,neighbour*N+i+1,1,[.985,.010,.003],this.symmetry?.24:.12);link.angular=true;this.air.push(link);}
     this.solid.push(new Link(t*N+i,t*N+i+1,1,sheet.fir[i]));
    }
    for(let i=0;i<N;i++){
     const air=new Link(t*N+i,next(t)*N+i,1,[.985,.010,.003],.10);air.angular=true;this.air.push(air);
     const solid=new Link(t*N+i,next(t)*N+i,1,sheet.fir[i]);solid.angular=true;this.solid.push(solid);
    }
   }
  }
  this.solidAdmit.fill(0);
  for(const e of this.air){
   const a=e.a%N,b=e.b%N,theta=e.angular?2*Math.PI/C:0;
   const distance=Math.hypot(c.x[a]-c.x[b],.5*Math.sqrt(c.r[a]**2+c.r[b]**2-2*c.r[a]*c.r[b]*Math.cos(theta)));
   e.setFilter(distance/SOUND*this.rate,e.taps);e.admittanceScale=Math.sqrt(e.weight*Math.PI/(C*RHO*SOUND));
  }
  for(const e of this.solid){
   const a=e.i,b=e.j;
   e.z=Math.sqrt(sheet.impedance[a]*sheet.impedance[b]);
   this.solidAdmit[e.a]+=e.z;this.solidAdmit[e.b]+=e.z;
   const distance=e.angular?c.r[a]*2*Math.PI/C:Math.hypot(c.x[a]-c.x[b],c.r[a]-c.r[b]);
   const speed=Math.sqrt((sheet.axial[Math.min(a,N-2)]*c.dx[a]**2/sheet.mass[a]+sheet.axial[Math.min(b,N-2)]*c.dx[b]**2/sheet.mass[b])/2);
   e.setFilter(distance/Math.max(.1,speed)*this.rate,sheet.fir[a].map((v,k)=>Math.sqrt(v*sheet.fir[b][k])));
  }
  this.airBank=new LinkBank(this.air);this.solidBank=new LinkBank(this.solid);
 }
 // A diagnostic initial displacement can excite any of the computed sheet
 // modes; no modal waveform is ever imposed on the running pressure source.
 displaceMode(index,metres){
  const shape=modeShape(this.c.sheet,index);
  if(this.symmetry&&this.c.sheet.modes[index].k!==0)throw Error('A circumferential mode requires the full sheet state.');
  for(let i=0;i<this.size;i++)this.displacement[i]+=metres*shape[i];this.configure(this.c);
 }
 advanceSheet(){
  const sheet=this.c.sheet,h=1/this.rate;
  if(this.symmetry){
   const ch=this.channels[0],{rows,columns,values}=ch.kernel,q=this.displacement,v=this.velocity,root=sheet.rootMass,input=ch.input;
   for(let i=0;i<N;i++)input[i]=(2*q[i]+h*v[i])*root[i];
   for(let i=0;i<N;i++){
    let response=0;for(let k=rows[i];k<rows[i+1];k++)response+=values[k]*input[columns[k]];
    const old=q[i];response/=root[i];q[i]=response-old;v[i]=(2*response-4*old)/h-v[i];
   }
   return;
  }
  // The cached receptance is the sum of every supported sheet mode. Each
  // step applies that response once; geometry edits rebuild it, not audio.
  for(const ch of this.channels){
   ch.nodalQ.fill(0);ch.nodalV.fill(0);
   for(let i=0;i<N;i++)for(let t=0;t<this.sectors;t++){
    const factor=(this.symmetry?1:ch.angular[t])*sheet.rootMass[i];
    ch.nodalQ[i]+=factor*this.displacement[t*N+i];ch.nodalV[i]+=factor*h*this.velocity[t*N+i];
   }
   for(let i=0;i<N;i++)ch.input[i]=2*ch.nodalQ[i]+ch.nodalV[i];
   const {rows,columns,values}=ch.kernel;
   for(let i=0;i<N;i++){
    let response=0;for(let k=rows[i];k<rows[i+1];k++)response+=values[k]*ch.input[columns[k]];
    const old=ch.nodalQ[i];ch.nodalQ[i]=(response-old)/sheet.rootMass[i];
    ch.nodalV[i]=(2*response-4*old-ch.nodalV[i])/(h*sheet.rootMass[i]);
   }
  }
  if(this.symmetry){this.displacement.set(this.channels[0].nodalQ);this.velocity.set(this.channels[0].nodalV);}
  else{
   this.displacement.fill(0);this.velocity.fill(0);
   for(const ch of this.channels)for(let t=0;t<C;t++)for(let i=0;i<N;i++){this.displacement[t*N+i]+=ch.angular[t]*ch.nodalQ[i];this.velocity[t*N+i]+=ch.angular[t]*ch.nodalV[i];}
  }
 }
 step(){
  if(this.fault)return 0;
  const h=1/this.rate,c=this.c,sheet=c.sheet;
  this.tick++;
  if(this.tick%96===0){this.random^=this.random<<13;this.random^=this.random>>>17;this.random^=this.random<<5;this.wander+=.018*((this.random>>>0)/2147483648-1)-.008*this.wander;}
  this.drive+=(this.pressureTarget*(1+c.wander*this.wander)-this.drive)*h/.035;
  this.sum.fill(0);this.admit.fill(0);this.areas.fill(0);
  for(let t=0;t<this.sectors;t++)for(let i=0;i<N;i++){
   const r=Math.max(.000025,c.r[i]+this.displacement[t*N+i]);this.radii[t*N+i]=r;this.areas[i]+=Math.PI*r*r/this.sectors;
  }
  this.airBank.receiveAir(this.sum,this.admit,this.radii,this.areas,this.symmetry);
  for(let t=0;t<this.sectors;t++){
   const first=t*N,last=first+N-1,yi=c.area[0]/(C*RHO*SOUND),yo=this.areas[N-1]/(C*RHO*SOUND);
   this.admit[first]+=yi;this.sum[first]+=yi*this.drive/2;this.admit[last]+=yo;
  }
  this.solidIncoming.fill(0);
  this.solidBank.receiveSolid(this.solidIncoming);
  this.advanceSheet();

  for(let t=0;t<this.sectors;t++)for(let i=0;i<N;i++){
   const j=t*N+i,surface=2*Math.PI*this.radii[j]*c.dx[i]/C,mass=sheet.mass[i],v=this.velocity[j],admit=this.admit[j];
   // Passive, collocated pressure-volume and tissue-force/velocity exchange.
   // Solve this scalar port only; the sheet's modes already performed its
   // elastic propagation. No global per-sample matrix solve is needed.
   const resistance=surface*surface/admit+this.solidAdmit[j]+sheet.viscous[i];
   const drive=surface*2*this.sum[j]/admit+2*this.solidIncoming[j];
   const next=((mass-h*resistance/2)*v+h*drive)/(mass+h*resistance/2),mid=(v+next)/2;
   this.pressure[j]=(2*this.sum[j]-surface*mid)/admit;
   this.dv[j]=next-v;this.velocity[j]=next;this.force[j]=mid;
   if(!Number.isFinite(next)||Math.abs(next)>60||Math.abs(this.displacement[j])>.03){this.fault='Sheet motion exceeded the model domain. Stop and adjust support or pressure.';return 0;}
   if(c.r[i]+this.displacement[j]<.000025){this.displacement[j]=.000025-c.r[i];const bounded=Math.max(0,this.velocity[j]);this.dv[j]+=bounded-this.velocity[j];this.velocity[j]=bounded;this.force[j]=(v+bounded)/2;this.pressure[j]=(2*this.sum[j]-surface*this.force[j])/admit;this.contactCount++;}
  }
  this.airBank.send(this.pressure,true);
  this.solidBank.send(this.force,false);
  let outlet=0;for(let t=0;t<this.sectors;t++)outlet+=this.pressure[t*N+N-1]/this.sectors;
  const output=(outlet-this.previousOutlet)*this.rate/48000/10;this.previousOutlet=outlet;
  if(!Number.isFinite(output)||Math.abs(output)>20){this.fault='Pressure exceeded the model domain. Stop and adjust the tube.';return 0;}return output;
 }
 snapshot(){
  const displacement=new Float64Array(N*2),pressure=new Float64Array(N*2),velocity=new Float64Array(N*2);let gap=Infinity;
  for(let i=0;i<N;i++)for(let side=0;side<2;side++){
   const from=(this.symmetry?0:side*C/2)*N+i,to=side*N+i;displacement[to]=this.displacement[from];pressure[to]=this.pressure[from];velocity[to]=this.velocity[from];
  }
  for(let i=P;i<2*P;i++)gap=Math.min(gap,2*this.c.r[i]+displacement[i]+displacement[N+i]);
  return {pressure:Array.from(pressure),displacement:Array.from(displacement),velocity:Array.from(velocity),gap,drive:this.drive,fault:this.fault,contacts:this.contactCount,modeCount:this.c.sheet.modes.length,rotationalSymmetry:this.symmetry};
 }
}

export function wav(samples,rate=48000){const data=new ArrayBuffer(44+samples.length*4),v=new DataView(data);const text=(i,s)=>{for(let j=0;j<s.length;j++)v.setUint8(i+j,s.charCodeAt(j));};text(0,'RIFF');v.setUint32(4,36+samples.length*4,true);text(8,'WAVEfmt ');v.setUint32(16,16,true);v.setUint16(20,3,true);v.setUint16(22,1,true);v.setUint32(24,rate,true);v.setUint32(28,rate*4,true);v.setUint16(32,4,true);v.setUint16(34,32,true);text(36,'data');v.setUint32(40,samples.length*4,true);for(let i=0;i<samples.length;i++)v.setFloat32(44+i*4,samples[i],true);return data;}
