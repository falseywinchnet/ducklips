// Radial displacement of a connected, periodic material sheet. The angular
// Fourier factorization is exact for ring supports and circumferentially
// uniform materials; it is not a pair of side-wall strips.
export const CIRCUMFERENCE = 8;
export function eigen(matrix,n){
 const a=Float64Array.from(matrix),v=new Float64Array(n*n);for(let i=0;i<n;i++)v[i*n+i]=1;
 const scale=Math.max(...Array.from({length:n},(_,i)=>Math.abs(a[i*n+i])),1);
 for(let sweep=0;sweep<60;sweep++){
  let max=0;
  for(let p=0;p<n-1;p++)for(let q=p+1;q<n;q++){
   const off=a[p*n+q];max=Math.max(max,Math.abs(off));if(Math.abs(off)<scale*1e-14)continue;
   const tau=(a[q*n+q]-a[p*n+p])/(2*off),t=(tau>=0?1:-1)/(Math.abs(tau)+Math.sqrt(1+tau*tau)),c=1/Math.sqrt(1+t*t),s=t*c;
   a[p*n+p]-=t*off;a[q*n+q]+=t*off;a[p*n+q]=a[q*n+p]=0;
   for(let i=0;i<n;i++)if(i!==p&&i!==q){const x=a[i*n+p],y=a[i*n+q];a[i*n+p]=a[p*n+i]=c*x-s*y;a[i*n+q]=a[q*n+i]=s*x+c*y;}
   for(let i=0;i<n;i++){const x=v[i*n+p],y=v[i*n+q];v[i*n+p]=c*x-s*y;v[i*n+q]=s*x+c*y;}
  }
  if(max<scale*1e-12)break;
  if(sweep===59)throw Error('Sheet eigendecomposition did not converge.');
 }
 const order=Array.from({length:n},(_,i)=>i).sort((i,j)=>a[i*n+i]-a[j*n+j]);
 return {values:Float64Array.from(order,i=>a[i*n+i]),vectors:Float64Array.from({length:n*n},(_,k)=>v[Math.floor(k/n)*n+order[k%n]])};
}
export function compileSheet(c,bands,materials,pointsPerBand){
 const n=c.x.length,C=CIRCUMFERENCE,angle=2*Math.PI/C,mass=new Float64Array(n),rootMass=new Float64Array(n),surface=new Float64Array(n),tension=new Float64Array(n),hoop=new Float64Array(n),foundation=new Float64Array(n),loss=new Float64Array(n),axial=new Float64Array(n-1),fir=[],impedance=new Float64Array(n);
 for(let i=0;i<n;i++){
  const b=bands[c.bandAt[i]],m=materials[b.material],circ=c.r[i]*angle;
  const slope=i&&i<n-1?(c.r[i+1]-c.r[i-1])/(c.x[i+1]-c.x[i-1]):0;
  const area=circ*c.dx[i]*Math.hypot(1,slope);
  surface[i]=circ*c.dx[i];mass[i]=m.density*m.thickness*area*b.mass;rootMass[i]=Math.sqrt(mass[i]);
  tension[i]=m.modulus*m.thickness*b.elasticity*(.035+Math.abs(Math.log(c.r[i]/b.radius))*.02);
  hoop[i]=tension[i]*c.dx[i]/circ;
  foundation[i]=tension[i]*surface[i]/(c.r[i]*c.r[i])+m.modulus*b.elasticity*surface[i]/.015;
  loss[i]=m.loss;fir.push(m.fir);
  impedance[i]=Math.sqrt(m.density*m.thickness*b.mass*tension[i])*circ*.08;
 }
 for(let i=0;i<n-1;i++)axial[i]=2*tension[i]*tension[i+1]/(tension[i]+tension[i+1]) * angle*(c.r[i]+c.r[i+1])/2 / Math.hypot(c.x[i+1]-c.x[i],c.r[i+1]-c.r[i]);
 const base=new Float64Array(n*n),supportDiagonal=new Float64Array(n),supportOff=new Float64Array(n-1);
 for(let i=0;i<n;i++)base[i*n+i]=foundation[i];
 for(let i=0;i<n-1;i++){const k=axial[i];base[i*n+i]+=k;base[(i+1)*n+i+1]+=k;base[i*n+i+1]-=k;base[(i+1)*n+i]-=k;}
 // Each handle is a ring support. Its footprint moves continuously within a band.
 for(let b=0;b<bands.length;b++){
  const band=bands[b],count=c.counts[b],at=band.anchor*(count-1),weights=Array.from({length:count},(_,j)=>Math.max(0,1-Math.abs(j-at)));
  const first=c.starts[b],k=band.support*(tension[first]*c.r[first]*angle/c.dx[first]+foundation[first]);
  for(let i=0;i<count;i++)for(let j=0;j<count;j++){const value=k*weights[i]*weights[j];base[(first+i)*n+first+j]+=value;if(i===j)supportDiagonal[first+i]+=value;else if(j===i+1)supportOff[first+i]+=value;}
 }
 const channels=[],modes=[],families=[];
 for(let k=0;k<=C/2;k++){
  const angularEigenvalue=4*Math.sin(Math.PI*k/C)**2,K=Float64Array.from(base);
  for(let i=0;i<n;i++)K[i*n+i]+=hoop[i]*angularEigenvalue;
  const normalized=Float64Array.from(K,(value,index)=>value/(rootMass[Math.floor(index/n)]*rootMass[index%n]));
  const e=eigen(normalized,n),damping=Float64Array.from(e.values,(w,m)=>{let l=0;for(let i=0;i<n;i++)l+=loss[i]*e.vectors[i*n+m]**2;return 2*l*Math.sqrt(w);});
  families.push({k,K,omega2:e.values,phi:e.vectors,damping});
  for(const kind of k===0||k===C/2?['cos']:['cos','sin']){
   const angular=Float64Array.from({length:C},(_,j)=>(k===0||k===C/2?1/Math.sqrt(C):Math.sqrt(2/C))*(kind==='cos'?Math.cos(k*j*angle):Math.sin(k*j*angle)));
   const channel=channels.length;channels.push({k,kind,angular,omega2:e.values,phi:e.vectors,damping});
   for(let m=0;m<n;m++)modes.push({channel,axial:m,k,kind,frequency:Math.sqrt(e.values[m])/(2*Math.PI)});
  }
 }
 modes.sort((a,b)=>a.frequency-b.frequency||a.channel-b.channel||a.axial-b.axial);
 const viscous=Float64Array.from(mass,(m,i)=>2*loss[i]*Math.sqrt(m*base[i*n+i]));
 return {bandAt:c.bandAt,tension,supportDiagonal,supportOff,n,circumference:C,mass,rootMass,surface,axial,hoop,foundation,base,families,channels,modes,impedance,fir,viscous};
}
export function modeShape(sheet,index){
 const mode=sheet.modes[index],ch=sheet.channels[mode.channel],shape=new Float64Array(sheet.n*sheet.circumference);let peak=0;
 for(let t=0;t<sheet.circumference;t++)for(let i=0;i<sheet.n;i++){const value=ch.angular[t]*ch.phi[i*sheet.n+mode.axial]/sheet.rootMass[i];shape[t*sheet.n+i]=value;peak=Math.max(peak,Math.abs(value));}
 for(let i=0;i<shape.length;i++)shape[i]/=peak;return shape;
}

// Compile the sum of ALL modes into its one-step receptance. In normalized
// coordinates R = Phi diag(1 / (1 + h² omega² / 4)) Phiᵀ. Midpoint advancement
// needs one application of R, not modal projections at every audio sample.
// Tiny spatial tails have an explicit maximum absolute row-error receipt.
export function receptance(channel,n,rate,{dense=false}={}){
 const h=1/rate,rows=[],columns=[],values=[],matrix=new Float64Array(n*n);let droppedRowBound=0;
 for(let i=0;i<n;i++){
  rows.push(values.length);let dropped=0;
  for(let j=0;j<n;j++){
   let value=0;for(let m=0;m<n;m++)value+=channel.phi[i*n+m]*channel.phi[j*n+m]/(1+h*h*channel.omega2[m]/4);
   matrix[i*n+j]=value;
   if(dense||Math.abs(value)>1e-14){columns.push(j);values.push(value);}else dropped+=Math.abs(value);
  }droppedRowBound=Math.max(droppedRowBound,dropped);
 }
 rows.push(values.length);
 if(droppedRowBound>1e-12)throw Error('Receptance compression exceeded its error bound.');
 return {rows:Uint32Array.from(rows),columns:Uint32Array.from(columns),values:Float64Array.from(values),matrix,droppedRowBound};
}

// The axial matrices are tridiagonal in each angular Fourier family. During
// slow elasticity variation their exact modal sum can be updated directly as
// the inverse of I + h² K_normalized / 4. This O(n²) control-rate recurrence
// avoids a new eigendecomposition; audio still applies one cached response.
export class VariableReceptance {
 constructor(sheet,channel,rate,baseKernel){
  this.sheet=sheet;this.channel=channel;this.h2=1/(4*rate*rate);const n=sheet.n;
  for(const key of ['pivot','lower','diagonal','off','inverseDiagonal','factor','axial','viscous','impedance'])this[key]=new Float64Array(n);
  this.matrix=new Float64Array(n*n);this.kernel=baseKernel;
  // Include a conservative spatial neighborhood and all remaining significant
  // tails. Verify the omitted row norm on every control-rate update.
  const rows=[],columns=[],values=[];
  for(let i=0;i<n;i++){rows.push(values.length);for(let j=0;j<n;j++)if(Math.abs(i-j)<=2||Math.abs(baseKernel.matrix[i*n+j])>1e-16){columns.push(j);values.push(baseKernel.matrix[i*n+j]);}}
  rows.push(values.length);this.rows=Uint32Array.from(rows);this.columns=Uint32Array.from(columns);this.values=Float64Array.from(values);this.next=Float64Array.from(values);this.delta=new Float64Array(values.length);this.droppedRowBound=0;
 }
 update(factors,mix=1){
  const s=this.sheet,n=s.n,h2=this.h2,angular=4*Math.sin(Math.PI*this.channel.k/s.circumference)**2;
  for(let i=0;i<n;i++)this.factor[i]=factors[s.bandAt[i]];
  for(let i=0;i<n-1;i++){
   const a=s.tension[i],b=s.tension[i+1],fa=this.factor[i],fb=this.factor[i+1];
   this.axial[i]=s.axial[i]*fa*fb*(a+b)/(a*fa+b*fb);
  }
  for(let i=0;i<n;i++){
   const f=this.factor[i],k=(s.foundation[i]+s.supportDiagonal[i]+s.hoop[i]*angular)*f+(i?this.axial[i-1]:0)+(i<n-1?this.axial[i]:0);
   this.diagonal[i]=1+h2*k/s.mass[i];
   if(i<n-1)this.off[i]=h2*(s.supportOff[i]*f-this.axial[i])/(s.rootMass[i]*s.rootMass[i+1]);
   this.pivot[i]=this.diagonal[i]-(i?this.off[i-1]*this.lower[i-1]:0);if(!(this.pivot[i]>0))throw Error('Elastic response lost positive definiteness.');
   this.lower[i]=i<n-1?this.off[i]/this.pivot[i]:0;
   this.viscous[i]=s.viscous[i]*Math.sqrt((k-s.hoop[i]*angular*f)/s.base[i*n+i]);this.impedance[i]=s.impedance[i]*Math.sqrt(f);
  }
  this.inverseDiagonal[n-1]=1/this.pivot[n-1];
  for(let i=n-2;i>=0;i--)this.inverseDiagonal[i]=1/this.pivot[i]+this.lower[i]**2*this.inverseDiagonal[i+1];
  for(let j=0;j<n;j++){
   let value=this.inverseDiagonal[j];this.matrix[j*n+j]=value;
   for(let i=j-1;i>=0;i--){value*=-this.lower[i];this.matrix[i*n+j]=this.matrix[j*n+i]=value;}
  }
  let bound=0;
  for(let i=0;i<n;i++){
   let at=this.rows[i],dropped=0;
   for(let j=0;j<n;j++){
    if(at<this.rows[i+1]&&this.columns[at]===j){const start=this.values[at]+mix*this.delta[at];this.values[at]=start;this.next[at]=this.matrix[i*n+j];this.delta[at]=this.next[at]-start;at++;}
    else dropped+=Math.abs(this.matrix[i*n+j]);
   }bound=Math.max(bound,dropped);
  }
  if(bound>1e-12)throw Error('Changing elasticity exceeded the material response compression bound.');this.droppedRowBound=Math.max(this.droppedRowBound,bound);
 }
}
