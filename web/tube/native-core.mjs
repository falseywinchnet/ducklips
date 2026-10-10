import {coreModule} from './core-wasm.mjs';
// A memory layout for the unchanged sample kernel. JS retains configuration,
// control-rate elasticity, provenance and UI; all views share the kernel state.
export class NativeCore {
 constructor(engine){
  this.engine=engine;this.wasm=new WebAssembly.Instance(coreModule).exports;this.memory=this.wasm.memory;this.cursor=this.wasm.base();
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
  const fields=[e.n,...[c.r,c.dx,s.mass,s.rootMass,c.area].map(pointer)];
  for(const key of ['pressure','sum','admit','displacement','velocity','solidIncoming','solidAdmit','radii','areas','dv','force','viscous'])fields.push(bind(e,key));
  fields.push(bind(ch,'input'),pointer(response.rows),pointer(response.columns),bind(response,'values'),e.variable?bind(response,'delta'):0);
  this.air=bank(e.airBank);this.solid=bank(e.solidBank);fields.push(this.air.byteOffset,this.solid.byteOffset);
  this.ports=copy(new Float64Array(32*7));fields.push(0,this.ports.byteOffset);this.core=copy(Uint32Array.from(fields));
 }
 step(h){
  const e=this.engine;if(e.ports.length>32)throw Error('Too many sheet ports.');this.core[25]=e.ports.length;
  for(let k=0;k<e.ports.length;k++){const p=e.ports[k],i=k*7;this.ports[i]=p.index;this.ports[i+1]=p.airY;this.ports[i+2]=p.airIn;this.ports[i+3]=p.solidZ;this.ports[i+4]=p.solidIn;}
  const contacts=this.wasm.step(this.core.byteOffset,h,e.drive,e.controlMix,+(e.inletMode==='reservoir'),+(e.outletMode==='open'),+e.variable);
  e.airBank.unitPos=this.air[2];
  for(let k=0;k<e.ports.length;k++){e.ports[k].airOut=this.ports[k*7+5];e.ports[k].solidOut=this.ports[k*7+6];}
  return contacts;
 }
}
