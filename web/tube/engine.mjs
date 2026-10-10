// A new geometry-driven experiment. No imports from the previous speech engine.
export const RHO = 1.204, SOUND = 343, N = 24, BANDS = 6, P = 4;
export const MATERIALS = {
  soft: {name:'Soft tissue', density:1040, modulus:8000, thickness:.002, loss:.12, fir:[.74,.17,.06], color:'#72c4b7'},
  fold: {name:'Elastic band', density:1040, modulus:22000, thickness:.001, loss:.055, fir:[.85,.10,.03], color:'#eab97a'},
  firm: {name:'Firm tissue', density:1100, modulus:65000, thickness:.002, loss:.16, fir:[.61,.23,.10], color:'#959dd9'}
};
export const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
export function initialDesign(){
  return {version:1, pressure:650, wander:.012, monitor:.25, blocks:[], bands:[
    {name:'01 · Return chamber',length:.050,radius:.009,material:'soft'},
    {name:'02 · Vocal band',length:.016,radius:.00032,material:'fold'},
    {name:'03 · Lower tract',length:.035,radius:.007,material:'soft'},
    {name:'04 · Oral body',length:.045,radius:.011,material:'soft'},
    {name:'05 · Front cavity',length:.035,radius:.008,material:'firm'},
    {name:'06 · Outlet',length:.025,radius:.006,material:'soft'}
  ].map(b=>({...b,mass:1,elasticity:1,upper:0,lower:0,upperAnchor:.5,lowerAnchor:.5,support:1}))};
}
export function validateDesign(d){
  if(!d || d.version!==1 || !Array.isArray(d.bands)||d.bands.length!==6||!Array.isArray(d.blocks)||d.blocks.length>128)throw Error('Expected a six-band tube design.');
  const finite=(v,a,b)=>Number.isFinite(v)&&v>=a&&v<=b;
  if(!finite(d.pressure,0,1600)||!finite(d.wander,0,.1))throw Error('Invalid pressure.');
  for(const b of d.bands)if(typeof b.name!=='string'||b.name.length>80||!finite(b.upperAnchor??.5,0,1)||!finite(b.lowerAnchor??.5,0,1)||!MATERIALS[b.material]||!finite(b.length,.008,.1)||!finite(b.radius,.00005,.025)||!finite(b.mass,.2,5)||!finite(b.elasticity,.2,5)||!finite(b.upper,-.02,.02)||!finite(b.lower,-.02,.02)||!finite(b.support,.2,5))throw Error('Invalid band properties.');
  for(const b of d.blocks)if(!finite(b.x,0,1)||!finite(b.volume,.01,8)||!finite(b.width,.02,.5))throw Error('Invalid volume block.');
  return d;
}
// Jacobi decomposition of the small, symmetric mass-normalized material matrix.
export function eigen(matrix,n){
  const a=Float64Array.from(matrix), v=new Float64Array(n*n);for(let i=0;i<n;i++)v[i*n+i]=1;
  for(let k=0;k<80;k++){
    let p=0,q=1,max=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(a[i*n+j])>max){max=Math.abs(a[i*n+j]);p=i;q=j;}
    if(max<1e-9)break;
    const angle=.5*Math.atan2(2*a[p*n+q],a[q*n+q]-a[p*n+p]), c=Math.cos(angle),s=Math.sin(angle);
    for(let i=0;i<n;i++){let x=a[i*n+p],y=a[i*n+q];a[i*n+p]=c*x-s*y;a[i*n+q]=s*x+c*y;}
    for(let j=0;j<n;j++){let x=a[p*n+j],y=a[q*n+j];a[p*n+j]=c*x-s*y;a[q*n+j]=s*x+c*y;}
    for(let i=0;i<n;i++){let x=v[i*n+p],y=v[i*n+q];v[i*n+p]=c*x-s*y;v[i*n+q]=s*x+c*y;}
  }
  const order=Array.from({length:n},(_,i)=>i).sort((i,j)=>a[i*n+i]-a[j*n+j]);
  return {values:order.map(i=>a[i*n+i]),vectors:Array.from({length:n*n},(_,k)=>v[Math.floor(k/n)*n+order[k%n]])};
}
export function compile(d){
  validateDesign(d); const x=[],dx=[],r=[], upper=[],lower=[],area=[],patches=[];let start=0;
  const length=d.bands.reduce((s,b)=>s+b.length,0);
  for(let b=0;b<6;b++){
    const band=d.bands[b];for(let j=0;j<P;j++){
      const i=b*P+j;dx[i]=band.length/P;x[i]=start+(j+.5)*dx[i];r[i]=band.radius;
      upper[i]=r[i]+band.upper;lower[i]=r[i]+band.lower;
    }start+=band.length;
  }
  // A block contributes its stated volume exactly on this finite-volume grid.
  const added=new Float64Array(N);
  for(const block of d.blocks){
    const weights=x.map((v,i)=>Math.max(0,1-Math.abs(v/length-block.x)/(block.width/2))*dx[i]);
    let sum=weights.reduce((s,v)=>s+v,0);
    if(sum===0){let nearest=0;for(let i=1;i<N;i++)if(Math.abs(x[i]/length-block.x)<Math.abs(x[nearest]/length-block.x))nearest=i;weights[nearest]=1;sum=1;}
    for(let i=0;i<N;i++)added[i]+=block.volume*1e-6*weights[i]/sum/dx[i];
  }
  for(let i=0;i<N;i++){
    const base=Math.max(.000025,(upper[i]+lower[i])/2), next=Math.sqrt(base*base+added[i]/Math.PI), expansion=next-base;
    upper[i]+=expansion;lower[i]+=expansion;r[i]=Math.max(.000025,(upper[i]+lower[i])/2);area[i]=Math.PI*r[i]*r[i];
  }
  for(let side=0;side<2;side++)for(let b=0;b<6;b++){
    const band=d.bands[b],mat=MATERIALS[band.material],step=band.length/P;
    const width=.015, mass=mat.density*mat.thickness*width*step*band.mass;
    // Anchor displacement changes rest shape and support span/prestrain, hence K.
    const offset=side?band.lower:band.upper, wall=side?lower:upper;
    let arc=step;for(let j=1;j<P;j++)arc+=Math.hypot(step,wall[b*P+j]-wall[b*P+j-1]);
    const stretch=arc/band.length*Math.sqrt(1+(offset/band.length)**2);
    const tension=mat.modulus*mat.thickness*width*.035*band.elasticity*stretch;
    const spring=tension/step, K=new Float64Array(P*P);
    for(let j=0;j<P;j++){
      K[j*P+j]=2*spring+spring*.08+mat.modulus*width*step/(.015)*band.elasticity;
      if(j)K[j*P+j-1]=-spring;if(j<P-1)K[j*P+j+1]=-spring;
    }
    K[0]+=spring*(band.support-1);K[P*P-1]+=spring*(band.support-1);
    const anchor=clamp(side?(band.lowerAnchor??.5):(band.upperAnchor??.5),0,1)*(P-1);
    const shape=Array.from({length:P},(_,j)=>Math.max(0,1-Math.abs(j-anchor)));
    for(let i=0;i<P;i++)for(let j=0;j<P;j++)K[i*P+j]+=spring*band.support*shape[i]*shape[j];
    const modes=eigen(Array.from(K,v=>v/mass),P);
    patches.push({side,band:b,mass,width,dx:step,omega2:modes.values,phi:modes.vectors,loss:mat.loss,fir:mat.fir,z:.006, frequencies:modes.values.map(w=>Math.sqrt(w)/(2*Math.PI))});
  }
  return {length,x,dx,r,upper,lower,area,patches,pressure:d.pressure,wander:d.wander};
}
class Link {
  constructor(a,b,delay,taps,weight=1){this.a=a;this.b=b;this.weight=weight;this.ab=new Float64Array(16);this.ba=new Float64Array(16);this.mask=15;this.pos=0;this.ia=0;this.ib=0;this.setFilter(delay,taps);}
  setFilter(delay,taps){
    this.delay=clamp(delay,1,4090);this.taps=taps;
    const size=2**Math.ceil(Math.log2(this.delay+5));
    if(size>this.ab.length){const ab=new Float64Array(size),ba=new Float64Array(size);for(let i=1;i<=this.ab.length;i++){ab[(size-i)&(size-1)]=this.ab[(this.pos-i+this.ab.length)&this.mask];ba[(size-i)&(size-1)]=this.ba[(this.pos-i+this.ba.length)&this.mask];}this.ab=ab;this.ba=ba;this.mask=size-1;this.pos=0;}this.d=Math.floor(this.delay);const f=this.delay-this.d;
    // Fuse fractional propagation delay and three-tap attenuation into one FIR.
    this.h0=(1-f)*taps[0];this.h1=f*taps[0]+(1-f)*taps[1];this.h2=f*taps[1]+(1-f)*taps[2];this.h3=f*taps[2];
  }
  receive(){const i=(this.pos-this.d)&this.mask,j=(i-1)&this.mask,k=(i-2)&this.mask,l=(i-3)&this.mask;
    this.ia=this.h0*this.ba[i]+this.h1*this.ba[j]+this.h2*this.ba[k]+this.h3*this.ba[l];
    this.ib=this.h0*this.ab[i]+this.h1*this.ab[j]+this.h2*this.ab[k]+this.h3*this.ab[l];
  }
  send(a,b){this.ab[this.pos]=a;this.ba[this.pos]=b;this.pos=(this.pos+1)&this.mask;}
}
export class TubeEngine {
  constructor(compiled,rate=24000){
    this.rate=rate;this.tick=0;this.random=0x513fe2;this.wander=0;this.drive=0;this.fault=null;this.contactCount=0;
    this.pressure=new Float64Array(N*2);this.sum=new Float64Array(N*2);this.admit=new Float64Array(N*2);this.airFlux=new Float64Array(N*2);
    this.displacement=new Float64Array(N*2);this.velocity=new Float64Array(N*2);this.force=new Float64Array(N*2);this.solidIncoming=new Float64Array(N*2);
    this.radii=new Float64Array(N);this.areas=new Float64Array(N);this.air=[];this.solid=[];this.previousOutlet=0;
    this.configure(compiled);
  }
  configure(c){
    this.c=c;this.pressureTarget=c.pressure;this.patches=c.patches.map((p,index)=>{
      const q=new Float64Array(P),v=new Float64Array(P),base=p.side*N+p.band*P;
      for(let m=0;m<P;m++)for(let j=0;j<P;j++){q[m]+=p.phi[j*P+m]*this.displacement[base+j];v[m]+=p.phi[j*P+m]*this.velocity[base+j];}
      const h=1/this.rate;
      const damp=p.omega2.map(w=>2*p.loss*Math.sqrt(w)+2*p.z/p.mass);
      return {...p,q,v,damp,den:p.omega2.map((w,m)=>1+h*damp[m]/2+h*h*w/4)};
    });
    if(!this.air.length){
      for(let i=0;i<N-1;i++)for(let side=0;side<2;side++){
        this.air.push(new Link(side*N+i,side*N+i+1,1,[.985,.010,.003],.42));
        this.air.push(new Link(side*N+i,(1-side)*N+i+1,1,[.985,.010,.003],.08));
      }
      for(let i=0;i<N;i++)this.air.push(new Link(i,N+i,1,[.985,.010,.003],.16));
      for(let side=0;side<2;side++)for(let b=0;b<5;b++)for(let j=0;j<P;j++){
        const p=c.patches[side*6+b],next=c.patches[side*6+b+1];
        this.solid.push(new Link(side*N+b*P+j,side*N+(b+1)*P+j,1,p.fir.map((v,k)=>Math.sqrt(v*next.fir[k]))));
      }
    }
    for(const e of this.air){
      const a=e.a%N,b=e.b%N, y1=(e.a<N?1:-1)*c.r[a]/2,y2=(e.b<N?1:-1)*c.r[b]/2;
      e.setFilter(Math.hypot(c.x[a]-c.x[b],y1-y2)/SOUND*this.rate,e.taps);
    }
    for(const e of this.solid){
      const a=e.a%N,b=e.b%N,pa=c.patches[Math.floor(e.a/N)*6+Math.floor(a/P)],pb=c.patches[Math.floor(e.b/N)*6+Math.floor(b/P)];
      const speed=Math.sqrt((pa.omega2[0]*pa.dx**2+pb.omega2[0]*pb.dx**2)/2);
      e.setFilter(Math.abs(c.x[a]-c.x[b])/Math.max(1,speed)*this.rate,pa.fir.map((v,k)=>Math.sqrt(v*pb.fir[k])));
    }
  }
  step(){
    if(this.fault)return 0;const h=1/this.rate,c=this.c;
    this.tick++;
    if(this.tick%96===0){this.random^=this.random<<13;this.random^=this.random>>>17;this.random^=this.random<<5;this.wander+=.018*((this.random>>>0)/2147483648-1)-.008*this.wander;}
    // Reservoir has finite response; the only stochastic input is slow pressure wander.
    this.drive+=(this.pressureTarget*(1+c.wander*this.wander)-this.drive)*h/.035;
    this.sum.fill(0);this.admit.fill(0);this.airFlux.fill(0);
    for(let i=0;i<N;i++){this.radii[i]=Math.max(.000025,(c.upper[i]+c.lower[i]+this.displacement[i]+this.displacement[N+i])/2);this.areas[i]=Math.PI*this.radii[i]**2;}
    for(const e of this.air){
      e.receive();const a=e.a%N,b=e.b%N;e.y=e.weight*Math.min(this.areas[a],this.areas[b])/(RHO*SOUND);
      e.scale=Math.sqrt(e.y);e.ia/=e.scale;e.ib/=e.scale;
      this.sum[e.a]+=e.y*e.ia;this.sum[e.b]+=e.y*e.ib;this.admit[e.a]+=e.y;this.admit[e.b]+=e.y;
    }
    // Two ends: a resistive pressure reservoir and a radiating exterior load.
    for(let s=0;s<2;s++){
      const first=s*N,last=first+N-1,yi=c.area[0]/(2*RHO*SOUND),yo=this.areas[N-1]/(2*RHO*SOUND);
      this.admit[first]+=yi;this.sum[first]+=yi*this.drive/2;
      this.admit[last]+=yo;
    }
    this.solidIncoming.fill(0);
    for(const e of this.solid){e.receive();this.solidIncoming[e.a]+=e.ia;this.solidIncoming[e.b]+=e.ib;}
    // Advance each supported material's modal response to its solid ports.
    for(const p of this.patches){
      const base=p.side*N+p.band*P;
      for(let j=0;j<P;j++){
        const i=base+j,a=i%N,gap=c.upper[a]+c.lower[a]+this.displacement[a]+this.displacement[N+a];
        this.force[i]=2*p.z*this.solidIncoming[i];
        if(gap<.00005){this.force[i]+=1000*(.00005-gap)+.03*Math.max(0,-this.velocity[a]-this.velocity[N+a]);this.contactCount++;}
      }
      for(let m=0;m<P;m++){
        let f=0;for(let j=0;j<P;j++)f+=p.phi[j*P+m]*this.force[base+j];
        const old=p.v[m],w=p.omega2[m];
        const v=((1-h*p.damp[m]/2-h*h*w/4)*old-h*w*p.q[m]+h*f/p.mass)/p.den[m];
        p.q[m]+=h*(old+v)/2;p.v[m]=v;
      }
      for(let j=0;j<P;j++){
        const i=base+j,a=i%N;let q=0,v=0;
        for(let m=0;m<P;m++){q+=p.phi[j*P+m]*p.q[m];v+=p.phi[j*P+m]*p.v[m];}
        // A local matched pressure/motion exchange: the same surface converts
        // pressure to force and velocity to volume. No explicit stiff air kick.
        const surface=Math.PI*this.radii[a]*c.dx[a];
        const pressure=(2*this.sum[i]-surface*v)/(this.admit[i]+h*surface*surface/(2*p.mass));
        const dv=h*surface*pressure/p.mass;
        this.pressure[i]=pressure;this.displacement[i]=q+h*dv/2;this.velocity[i]=v+dv;
        if(!Number.isFinite(q)||Math.abs(q)>.03||Math.abs(v+dv)>60){this.fault='Motion exceeded the model domain. Stop and adjust support or pressure.';return 0;}
      }
      for(let m=0;m<P;m++){
        let q=0,v=0;for(let j=0;j<P;j++){q+=p.phi[j*P+m]*this.displacement[base+j];v+=p.phi[j*P+m]*this.velocity[base+j];}p.q[m]=q;p.v[m]=v;
      }
    }
    for(const e of this.air){
      const oa=this.pressure[e.a]-e.ia,ob=this.pressure[e.b]-e.ib;
      this.airFlux[e.a]+=e.y*(oa-e.ia);this.airFlux[e.b]+=e.y*(ob-e.ib);e.send(oa*e.scale,ob*e.scale);
    }
    for(const e of this.solid)e.send(this.velocity[e.a]-e.ia,this.velocity[e.b]-e.ib);
    // Outlet pressure difference is a fixed two-tap radiation observation, not excitation.
    const outlet=(this.pressure[N-1]+this.pressure[2*N-1])/2;
    const output=(outlet-this.previousOutlet)*this.rate/48000/10;this.previousOutlet=outlet;
    if(!Number.isFinite(output)||Math.abs(output)>20){this.fault='Pressure exceeded the model domain. Stop and adjust the tube.';return 0;}
    return output;
  }
  snapshot(){return {pressure:Array.from(this.pressure),displacement:Array.from(this.displacement),velocity:Array.from(this.velocity),gap:2*Math.min(...this.radii.slice(4,8)),drive:this.drive,fault:this.fault,contacts:this.contactCount};}
}
export function wav(samples,rate=48000){const data=new ArrayBuffer(44+samples.length*4),v=new DataView(data);const text=(i,s)=>{for(let j=0;j<s.length;j++)v.setUint8(i+j,s.charCodeAt(j));};text(0,'RIFF');v.setUint32(4,36+samples.length*4,true);text(8,'WAVEfmt ');v.setUint32(16,16,true);v.setUint16(20,3,true);v.setUint16(22,1,true);v.setUint32(24,rate,true);v.setUint32(28,rate*4,true);v.setUint16(32,4,true);v.setUint16(34,32,true);text(36,'data');v.setUint32(40,samples.length*4,true);for(let i=0;i<samples.length;i++)v.setFloat32(44+i*4,samples[i],true);return data;}
