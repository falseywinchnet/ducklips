// Causal unflanged-pipe reflection, Silva et al., equations 12–15:
// https://arxiv.org/abs/0811.3625 . The positive gamma-shaped impulse is
// integrated into a finite FIR; its negative is the pressure reflection.
// This is an aperture model, not a measured human-mouth radiation model.
const cache=new Map();
export const RADIATION_TAPS=64;
export class ApertureRadiation {
 constructor(rate,sound=343){
  this.rate=rate;this.sound=sound;this.maximumRadius=Math.min(.06,(RADIATION_TAPS-2)*sound/(12*rate));this.bins=512;
  const key=rate+':'+sound;if(cache.has(key)){this.table=cache.get(key);return;}
  const table=new Float64Array((this.bins+1)*RADIATION_TAPS),eta=.6133,beta=.5,alpha=eta/beta,nu=2*eta*eta/beta-1,steps=2048;
  const u=Float64Array.from({length:steps},(_,i)=>(i+.5)*20/steps),weight=Float64Array.from(u,x=>x**nu*Math.exp(-alpha*x));
  for(let b=0;b<=this.bins;b++){
   const base=b*RADIATION_TAPS,scale=b/this.bins*this.maximumRadius*rate/sound;
   for(let i=0;i<steps;i++){const lag=u[i]*scale,k=Math.floor(lag),f=lag-k;if(k<RADIATION_TAPS)table[base+k]+=weight[i]*(1-f);if(k+1<RADIATION_TAPS)table[base+k+1]+=weight[i]*f;}
   let sum=0;for(let k=0;k<RADIATION_TAPS;k++)sum+=table[base+k];
   for(let k=0;k<RADIATION_TAPS;k++)table[base+k]/=sum;
  }
  this.table=table;cache.set(key,table);
 }
 coefficients(radius,out=new Float64Array(RADIATION_TAPS)){
  if(!Number.isFinite(radius)||radius<0||radius>this.maximumRadius)throw Error('Aperture radius exceeds the radiation table.');
  const x=radius/this.maximumRadius*this.bins,b=Math.min(this.bins-1,Math.floor(x)),f=x-b,lo=b*RADIATION_TAPS,hi=lo+RADIATION_TAPS;
  for(let k=0;k<RADIATION_TAPS;k++)out[k]=this.table[lo+k]*(1-f)+this.table[hi+k]*f;
  return out;
 }
 gather(node,radius){
  this.coefficients(radius,node.radiationCoefficients);const c=node.radiationCoefficients;
  node.radiationDirect=-c[0];let history=0;
  for(let k=1;k<RADIATION_TAPS;k++)history-=c[k]*node.radiation[k-1];
  node.radiationHistory=history;
 }
}
