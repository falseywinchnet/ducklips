// The same per-sample sheet/transport operations as engine.mjs, in contiguous
// WebAssembly memory. Double precision, no fast-math, no modal truncation.
typedef unsigned int u32;
#define PI 3.141592653589793
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define ABS(a) __builtin_fabs(a)
#define SQRT(a) __builtin_sqrt(a)
typedef struct {
 u32 count,unit,unitPos,unitMask;
 int *a,*b,*i,*j,*offset,*mask,*pos,*delay;
 double *h0,*h1,*h2,*h3,*gain,*scale,*ia,*ib,*ab,*ba,*unitAB,*unitBA;
} Bank;
typedef struct {
 u32 n;
 double *radius,*dx,*mass,*root,*area0;
 double *pressure,*sum,*admit,*q,*v,*solidIncoming,*solidAdmit,*radii,*areas,*dv,*force,*viscous,*input;
 u32 *rows,*columns;double *values,*delta;
 Bank *air,*solid;
 u32 portCount;double *ports;double *targets,*muscleVelocity,*wallPressureOffset,*slotWidth,*airRadii;
} Core;
static __attribute__((always_inline)) inline void receive(Bank *b,double *sum,double *admit,double *radii,int air){
 u32 n=b->count;
 if(air&&b->unit){
  const u32 p0=((b->unitPos-1)&b->unitMask)*n,p1=((b->unitPos-2)&b->unitMask)*n,p2=((b->unitPos-3)&b->unitMask)*n;
  const double h0=b->h0[0],h1=b->h1[0],h2=b->h2[0];
  const double *restrict ab=b->unitAB,*restrict ba=b->unitBA;
  for(u32 k=0;k<n;k++){
   int a=b->a[k],other=b->b[k];double s=b->gain[k]*MIN(radii[b->i[k]],radii[b->j[k]]),y=s*s;
   double va=(h0*ba[p0+k]+h1*ba[p1+k]+h2*ba[p2+k])/s;
   double vb=a==other?va:(h0*ab[p0+k]+h1*ab[p1+k]+h2*ab[p2+k])/s;
   b->scale[k]=s;b->ia[k]=va;b->ib[k]=vb;sum[a]+=y*va;sum[other]+=y*vb;admit[a]+=y;admit[other]+=y;
  }
  return;
 }
 for(u32 k=0;k<n;k++){
  double s=1,y=b->gain[k];int a=b->a[k],other=b->b[k];
  if(air){s=b->gain[k]*MIN(radii[b->i[k]],radii[b->j[k]]);y=s*s;}
  int o=b->offset[k],m=b->mask[k],p=(b->pos[k]-b->delay[k])&m;
  int j0=o+p,j1=o+((p-1)&m),j2=o+((p-2)&m),j3=o+((p-3)&m);
  double va=(b->h0[k]*b->ba[j0]+b->h1[k]*b->ba[j1]+b->h2[k]*b->ba[j2]+b->h3[k]*b->ba[j3])/s;
  double vb=a==other?va:(b->h0[k]*b->ab[j0]+b->h1[k]*b->ab[j1]+b->h2[k]*b->ab[j2]+b->h3[k]*b->ab[j3])/s;
  b->scale[k]=s;b->ia[k]=va;b->ib[k]=vb;sum[a]+=y*va;sum[other]+=y*vb;
  if(air){admit[a]+=y;admit[other]+=y;}
 }
}
static __attribute__((always_inline)) inline void send(Bank *b,double *field,int air){
 u32 n=b->count,o=b->unitPos*n;
 for(u32 k=0;k<n;k++){
  double s=air?b->scale[k]:1;
  double a=(field[b->a[k]]-b->ia[k])*s,v=(field[b->b[k]]-b->ib[k])*s;
  if(air&&b->unit){b->unitAB[o+k]=a;b->unitBA[o+k]=v;}
  else {int at=b->offset[k]+b->pos[k];b->ab[at]=a;b->ba[at]=v;b->pos[k]=(b->pos[k]+1)&b->mask[k];}
 }
 if(air&&b->unit)b->unitPos=(b->unitPos+1)&b->unitMask;
}
__attribute__((export_name("step"))) int step(Core *c,double h,double drive,double mix,int inlet,int outlet,int variable){
 u32 n=c->n;int contacts=0;
 for(u32 i=0;i<n;i++){c->muscleVelocity[i]+=h*(6400*(c->targets[i]-c->radius[i])-160*c->muscleVelocity[i]);double dr=h*c->muscleVelocity[i];c->radius[i]+=dr;c->q[i]-=dr;c->sum[i]=0;c->admit[i]=0;c->solidIncoming[i]=0;c->radii[i]=MAX(c->slotWidth[i]>.0?.0000001:.000025,c->radius[i]+c->q[i]);c->areas[i]=c->slotWidth[i]>.0?2*c->slotWidth[i]*c->radii[i]:PI*c->radii[i]*c->radii[i];c->airRadii[i]=c->slotWidth[i]>.0?SQRT(c->areas[i]/PI):c->radii[i];}
 receive(c->air,c->sum,c->admit,c->airRadii,1);
 if(inlet){double yi=c->area0[0]/(8*1.204*343);c->admit[0]+=yi;c->sum[0]+=yi*drive/2;}
 if(outlet)c->admit[n-1]+=c->areas[n-1]/(8*1.204*343);
 receive(c->solid,c->solidIncoming,0,0,0);
 for(u32 p=0;p<c->portCount;p++){double *port=c->ports+7*p;int i=(int)port[0];c->admit[i]+=port[1];c->sum[i]+=port[1]*port[2];c->solidAdmit[i]+=port[3];c->solidIncoming[i]+=port[3]*port[4];}
 for(u32 i=0;i<n;i++)c->input[i]=(2*c->q[i]+h*c->v[i])*c->root[i];
 for(u32 i=0;i<n;i++){
  double response=0;for(u32 k=c->rows[i];k<c->rows[i+1];k++)response+=(c->values[k]+(variable?mix*c->delta[k]:0))*c->input[c->columns[k]];
  double old=c->q[i];response/=c->root[i];c->q[i]=response-old;c->v[i]=(2*response-4*old)/h-c->v[i];
 }
 for(u32 i=0;i<n;i++){
  double surface=(c->slotWidth[i]>.0?2*c->slotWidth[i]:2*PI*c->radii[i])*c->dx[i]/8,mass=c->mass[i],v=c->v[i],admit=c->admit[i];
  double resistance=surface*surface/admit+c->solidAdmit[i]+c->viscous[i];
  double force=surface*2*c->sum[i]/admit+2*c->solidIncoming[i]+surface*c->wallPressureOffset[i];
  double next=((mass-h*resistance/2)*v+h*force)/(mass+h*resistance/2),mid=(v+next)/2;
  c->pressure[i]=(2*c->sum[i]-surface*mid)/admit;c->dv[i]=next-v;c->v[i]=next;c->force[i]=mid;
  if(next!=next||ABS(next)>60||ABS(c->q[i])>.03)return -1;
  double floor=c->slotWidth[i]>.0?.0000001:.000025;if(c->radius[i]+c->q[i]<floor){c->q[i]=floor-c->radius[i];double bounded=MAX(0,c->v[i]);c->dv[i]+=bounded-c->v[i];c->v[i]=bounded;c->force[i]=(v+bounded)/2;c->pressure[i]=(2*c->sum[i]-surface*c->force[i])/admit;contacts++;}
 }
 send(c->air,c->pressure,1);send(c->solid,c->force,0);
 for(u32 p=0;p<c->portCount;p++){double *port=c->ports+7*p;int i=(int)port[0];port[5]=c->pressure[i]-port[2];port[6]=c->force[i]-port[4];c->solidAdmit[i]-=port[3];}
 return contacts;
}
extern unsigned char __heap_base;
__attribute__((export_name("base"))) u32 base(void){return (u32)&__heap_base;}
// Oral sheets share pressure junctions on a coarser, physically spaced air grid.
// First advance their material response; the air junction then solves its scalar
// pressure from the summed wall-flow admittance, and finish applies that pressure.
__attribute__((export_name("walls_prepare"))) void walls_prepare(Core *c,double h,double mix,int variable){
 u32 n=c->n;
 for(u32 i=0;i<n;i++){
  c->muscleVelocity[i]+=h*(6400*(c->targets[i]-c->radius[i])-160*c->muscleVelocity[i]);double dr=h*c->muscleVelocity[i];c->radius[i]+=dr;c->q[i]-=dr;
  c->radii[i]=MAX(c->slotWidth[i]>.0?.0000001:.000025,c->radius[i]+c->q[i]);c->areas[i]=c->slotWidth[i]>.0?2*c->slotWidth[i]*c->radii[i]:PI*c->radii[i]*c->radii[i];c->airRadii[i]=c->slotWidth[i]>.0?SQRT(c->areas[i]/PI):c->radii[i];c->solidIncoming[i]=0;
 }
 receive(c->solid,c->solidIncoming,0,0,0);
 for(u32 p=0;p<c->portCount;p++){double *port=c->ports+7*p;int i=(int)port[0];c->solidAdmit[i]+=port[3];c->solidIncoming[i]+=port[3]*port[4];}
 for(u32 i=0;i<n;i++)c->input[i]=(2*c->q[i]+h*c->v[i])*c->root[i];
 for(u32 i=0;i<n;i++){
  double response=0;for(u32 k=c->rows[i];k<c->rows[i+1];k++)response+=(c->values[k]+(variable?mix*c->delta[k]:0))*c->input[c->columns[k]];
  double old=c->q[i];response/=c->root[i];c->q[i]=response-old;c->v[i]=(2*response-4*old)/h-c->v[i];
  double surface=(c->slotWidth[i]>.0?2*c->slotWidth[i]:2*PI*c->radii[i])*c->dx[i]/8,den=c->mass[i]+h*(c->solidAdmit[i]+c->viscous[i])/2;
  c->force[i]=(c->mass[i]*c->v[i]+h*c->solidIncoming[i])/den;
  c->dv[i]=h*surface/(2*den);
  c->sum[i]=surface*c->force[i];c->admit[i]=surface*c->dv[i];
 }
}
__attribute__((export_name("walls_finish"))) int walls_finish(Core *c){
 int contacts=0;for(u32 i=0;i<c->n;i++){
  double mid=c->force[i]+c->dv[i]*c->pressure[i],next=2*mid-c->v[i];
  double floor=c->slotWidth[i]>.0?.0000001:.000025;if(c->radius[i]+c->q[i]<floor){c->q[i]=floor-c->radius[i];next=MAX(0,next);mid=(c->v[i]+next)/2;contacts++;}
  c->v[i]=next;c->force[i]=mid;if(next!=next||ABS(next)>60||ABS(c->q[i])>.03)return -1;
 }
 send(c->solid,c->force,0);
 for(u32 p=0;p<c->portCount;p++){double *port=c->ports+7*p;int i=(int)port[0];port[6]=c->force[i]-port[4];c->solidAdmit[i]-=port[3];}
 return contacts;
}
