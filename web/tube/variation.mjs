// Slow material variation only. These values change elasticity, never the
// inlet waveform or rendered PCM. Every component and their sum are bounded.
export class ElasticityVariation {
 constructor(count,seed=73021){
  this.random=seed|0||1;this.factors=new Float64Array(count).fill(1);
  this.states=Array.from({length:count},()=>({phase:this.uniform(),frequency:.35+this.uniform()*.9,noise:0,target:0,hold:0,saw:0,bubble:0,age:2,duration:.3,sign:1,wait:.2+this.uniform(),value:0}));
 }
 uniform(){let x=this.random;x^=x<<13;x^=x>>>17;x^=x<<5;this.random=x;return (x>>>0)/4294967296;}
 advance(dt,depth){
  for(let i=0;i<this.states.length;i++){
   const s=this.states[i];s.hold-=dt;if(s.hold<=0){s.target=2*this.uniform()-1;s.hold=.08+.32*this.uniform();}
   s.noise+=(s.target-s.noise)*(1-Math.exp(-dt/.08));
   s.phase=(s.phase+dt*s.frequency)%1;s.saw+=(2*s.phase-1-s.saw)*(1-Math.exp(-dt/.045));
   s.wait-=dt;s.age+=dt;if(s.wait<=0){s.age=0;s.duration=.16+.4*this.uniform();s.sign=this.uniform()<.5?-1:1;s.wait=s.duration+.15+this.uniform()*1.2;}
   const x=s.age/s.duration;s.bubble=x<1?s.sign*(256/27)*x*(1-x)**3:0;
   const mixed=.4*s.noise+.25*s.saw+.35*s.bubble;s.value+=(mixed-s.value)*(1-Math.exp(-dt/.025));
   this.factors[i]=1+depth*s.value;
  }return this.factors;
 }
}
