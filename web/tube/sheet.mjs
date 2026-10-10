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
  const b=bands[Math.floor(i/pointsPerBand)],m=materials[b.material],circ=c.r[i]*angle;
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
 const base=new Float64Array(n*n);
 for(let i=0;i<n;i++)base[i*n+i]=foundation[i];
 for(let i=0;i<n-1;i++){const k=axial[i];base[i*n+i]+=k;base[(i+1)*n+i+1]+=k;base[i*n+i+1]-=k;base[(i+1)*n+i]-=k;}
 // Each handle is a ring support. Its footprint moves continuously within a band.
 for(let b=0;b<bands.length;b++){
  const band=bands[b],at=band.anchor*(pointsPerBand-1),weights=Array.from({length:pointsPerBand},(_,j)=>Math.max(0,1-Math.abs(j-at)));
  const first=b*pointsPerBand,k=band.support*(tension[first]*c.r[first]*angle/c.dx[first]+foundation[first]);
  for(let i=0;i<pointsPerBand;i++)for(let j=0;j<pointsPerBand;j++)base[(first+i)*n+first+j]+=k*weights[i]*weights[j];
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
 return {n,circumference:C,mass,rootMass,surface,axial,hoop,foundation,base,families,channels,modes,impedance,fir,viscous};
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
