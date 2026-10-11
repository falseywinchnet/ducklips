import {coreModule} from './core-wasm.mjs';
// A memory layout for the unchanged sample kernel. JS retains configuration,
// control-rate elasticity, provenance and UI; all views share the kernel state.
export function nativeArena(){const wasm=new WebAssembly.Instance(coreModule).exports;return {wasm,cursor:wasm.base()};}
export class NativeCore {
 constructor(engine,arena=nativeArena()){
  this.engine=engine;this.wasm=arena.wasm;this.memory=this.wasm.memory;this.cursor=arena.cursor;
  const copy=source=>{if(Array.isArray(source))source=Float64Array.from(source);const Type=source.constructor;this.cursor=(this.cursor+7)&~7;const view=new Type(this.memory.buffer,this.cursor,source.length);this.cursor+=source.byteLength;view.set(source);return view;};
  const bind=(object,key)=>{object[key]=copy(object[key]);return object[key].byteOffset;};
  const pointer=source=>copy(source).byteOffset;
  const bank=b=>{
   const fields=[b.length,+!!b.unitAir,b.unitPos||0,b.unitMask||0];
   for(const key of ['a','b','i','j','offset','mask','pos','delay','h0','h1','h2','h3','gain','scale','ia','ib','ab','ba'])fields.push(bind(b,key));
   fields.push(b.unitAB?bind(b,'unitAB'):0,b.unitBA?bind(b,'unitBA'):0);
   for(let k=0;k<b.length;k++){const e=b.links[k],start=b.offset[k],length=b.mask[k]+1;e.ab=b.ab.subarray(start,start+length);e.ba=b.ba.subarray(start,start+length);}
   return copy(Uint32Array.from(fields));
  };
  const e=engine,c=e.c,s=c.sheet,ch=e.channels[0],response=e.variable?ch.variable:ch.kernel;
  const fields=[e.n,bind(e,'restRadii'),...[c.dx,s.mass,s.rootMass,c.area].map(pointer)];
  for(const key of ['pressure','sum','admit','displacement','velocity','solidIncoming','solidAdmit','radii','areas','dv','force','viscous'])fields.push(bind(e,key));
  fields.push(bind(ch,'input'),pointer(response.rows),pointer(response.columns),bind(response,'values'),e.variable?bind(response,'delta'):0);
  this.air=bank(e.airBank);this.solid=bank(e.solidBank);fields.push(this.air.byteOffset,this.solid.byteOffset);
  this.ports=copy(new Float64Array(32*7));fields.push(0,this.ports.byteOffset,bind(e,'radiusTargets'),bind(e,'radiusVelocity'),bind(e,'wallPressureOffset'),pointer(c.slotWidth),bind(e,'airRadii'));this.core=copy(Uint32Array.from(fields));arena.cursor=this.cursor;
 }
 packPorts(){
  const e=this.engine;if(e.ports.length>32)throw Error('Too many sheet ports.');this.core[25]=e.ports.length;
  for(let k=0;k<e.ports.length;k++){const p=e.ports[k],i=k*7;this.ports[i]=p.index;this.ports[i+1]=p.airY;this.ports[i+2]=p.airIn;this.ports[i+3]=p.solidZ;this.ports[i+4]=p.solidIn;}
 }
 unpackPorts(){const e=this.engine;for(let k=0;k<e.ports.length;k++){e.ports[k].airOut=this.ports[k*7+5];e.ports[k].solidOut=this.ports[k*7+6];}}
 prepareWalls(h){this.packPorts();const e=this.engine;this.wasm.walls_prepare(this.core.byteOffset,h,e.controlMix,+e.variable,+!!e.contactMode);}
 finishWalls(){const contacts=this.wasm.walls_finish(this.core.byteOffset,+!!this.engine.contactMode);this.unpackPorts();return contacts;}
 step(h){
  this.packPorts();const e=this.engine;
  const contacts=this.wasm.step(this.core.byteOffset,h,e.drive,e.controlMix,+(e.inletMode==='reservoir'),+(e.outletMode==='open'),+e.variable,+!!e.contactMode);
  e.airBank.unitPos=this.air[2];
  for(let k=0;k<e.ports.length;k++){e.ports[k].airOut=this.ports[k*7+5];e.ports[k].solidOut=this.ports[k*7+6];}
  return contacts;
 }
}

// One shared material arena avoids ten JS/WASM crossings per physical sample.
// Each sheet keeps every coefficient, mode response, and transport history.
export class NativeSheets {
 constructor(engines){
  this.engines=engines;const arena=nativeArena();this.wasm=arena.wasm;
  for(const e of engines){e.airBank.sync();e.solidBank.sync();e.native=new NativeCore(e,arena);}
  this.natives=engines.map(e=>e.native);
  const copy=source=>{arena.cursor=(arena.cursor+7)&~7;const a=new source.constructor(this.wasm.memory.buffer,arena.cursor,source.length);a.set(source);arena.cursor+=source.byteLength;return a;};
  this.cores=copy(Uint32Array.from(this.natives,n=>n.core.byteOffset));this.mixes=copy(new Float64Array(engines.length));this.flags=copy(new Uint32Array(engines.length));this.contacts=copy(new Int32Array(engines.length));
 }
 current(){return this.engines.every((e,i)=>e.native===this.natives[i]);}
 prepare(h){
  for(let i=0;i<this.engines.length;i++){const e=this.engines[i];e.native.packPorts();this.mixes[i]=e.controlMix;this.flags[i]=+e.variable|(+!!e.contactMode<<1);}
  this.wasm.walls_prepare_batch(this.cores.byteOffset,this.cores.length,h,this.mixes.byteOffset,this.flags.byteOffset);
 }
 finish(){
  const result=this.wasm.walls_finish_batch(this.cores.byteOffset,this.cores.length,this.flags.byteOffset,this.contacts.byteOffset);
  for(let i=0;i<this.engines.length;i++){this.natives[i].unpackPorts();if(this.contacts[i]>0)this.engines[i].contactCount+=this.contacts[i];}
  return result;
 }
}
