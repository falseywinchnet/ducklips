#include "ducklips/contour_kinetics.hpp"
#include <algorithm>
#include <cmath>
namespace ducklips {
namespace {
// A scalar at two coordinate endpoints and its divided difference. Unlike
// subtracting two computed energies, these product/quotient rules retain
// the derivative when the endpoints coincide. The symmetric product rule
// is an exact finite chain rule, not a finite-difference perturbation.
struct Pair {double a{},b{},d{};Pair()=default;Pair(double x):a(x),b(x){} Pair(double x,double y,double s):a(x),b(y),d(s){} };
Pair operator+(Pair x,Pair y){return {x.a+y.a,x.b+y.b,x.d+y.d};}
Pair operator-(Pair x,Pair y){return {x.a-y.a,x.b-y.b,x.d-y.d};}
Pair operator*(Pair x,Pair y){return {x.a*y.a,x.b*y.b,.5*((x.a+x.b)*y.d+(y.a+y.b)*x.d)};}
Pair operator/(Pair x,Pair y){return {x.a/y.a,x.b/y.b,(x.d*.5*(y.a+y.b)-y.d*.5*(x.a+x.b))/(y.a*y.b)};}
Pair logarithm(Pair x){
    const double delta=x.b-x.a,relative=delta/x.a;
    const double secant=delta==0?1/x.a:(std::abs(relative)<.5?std::log1p(relative):std::log(x.b)-std::log(x.a))/delta;
    return {std::log(x.a),std::log(x.b),x.d*secant};
}
Pair clamp(Pair x,double low,double high){
    const double a=std::clamp(x.a,low,high),b=std::clamp(x.b,low,high);
    if(x.a==x.b){return {a,b,x.a>low && x.a<high?x.d:(x.a==low || x.a==high?.5*x.d:0.)};}
    return {a,b,x.d*(b-a)/(x.b-x.a)};
}
double phi(double t){
    if(std::abs(t)<.25){double sum=1,p=1;for(int k=1;k<=18;++k){p*=t*t;sum+=p/double(2*k+1);}return sum;}
    return .5*(std::log1p(t)-std::log1p(-t))/t;
}
Pair phi(Pair t){
    if(std::max(std::abs(t.a),std::abs(t.b))<.25){Pair sum(1),p(1),square=t*t;for(int k=1;k<=18;++k){p=p*square;sum=sum+p/Pair(double(2*k+1));}return sum;}
    if(std::min(std::abs(t.a),std::abs(t.b))<.01){return {phi(t.a),phi(t.b),t.d*(phi(t.b)-phi(t.a))/(t.b-t.a)};}
    return Pair(.5)*(logarithm(Pair(1)+t)-logarithm(Pair(1)-t))/t;
}
// Integral of 1/A for a positive affine segment, stable at equal areas.
Pair reciprocal_integral(Pair length,Pair a,Pair b){return Pair(2)*length*phi((b-a)/(b+a))/(a+b);}
Pair positive_part(Pair x){
    const double a=std::max(0.,x.a),b=std::max(0.,x.b);
    return {a,b,x.a==x.b?(x.a>0?x.d:0.):x.d*(b-a)/(x.b-x.a)};
}
bool valid_coordinates(const std::array<double,4>& x){for(double v:x){if(!std::isfinite(v)){return false;}}return x[0]>0 && x[1]>0;}
}
bool ContourKinetics::configure(const FoldParameters& f,const AirProperties& air,double loss) noexcept {
    for(double v:{f.lower_thickness,f.upper_thickness,f.length,air.density}){if(!std::isfinite(v) || v<=0){return false;}}
    if(!std::isfinite(air.viscosity) || air.viscosity<0 || !std::isfinite(loss) || loss<0){return false;}
    thickness_={f.lower_thickness,f.upper_thickness};width_=f.length;density_=air.density;viscosity_=air.viscosity;outlet_loss_=loss;return true;
}
bool ContourKinetics::metric(const std::array<double,4>& first,const std::array<double,4>& last,
    std::size_t coordinate,std::array<double,3>& value,std::array<double,3>& slope,std::array<double,2>* volume_slope) const noexcept {
    if(density_<=0 || !valid_coordinates(first) || !valid_coordinates(last)){return false;}
    // Widely separated coordinates, including sign changes, have a stable
    // endpoint secant. This also keeps piecewise contact branches separate.
    if(coordinate<4 && std::abs(last[coordinate]-first[coordinate])>1e-4*std::max(std::abs(first[coordinate]),std::abs(last[coordinate]))){
        std::array<double,3> start,unused;if(!metric(first,first,4,start,unused) || !metric(last,last,4,value,unused)){return false;}
        const double delta=last[coordinate]-first[coordinate];for(std::size_t j=0;j<3;++j){slope[j]=(value[j]-start[j])/delta;if(!std::isfinite(slope[j])){return false;}}
        if(volume_slope){volume_slope->fill(0.);if(coordinate>=2){GlottalContourGeometry a,b;
            if(!GlottalContour{thickness_[0],thickness_[1],first[2],first[3]}.geometry(a) || !GlottalContour{thickness_[0],thickness_[1],last[2],last[3]}.geometry(b)){return false;}
            for(std::size_t i=0;i<2;++i){(*volume_slope)[i]=(b.chamber[i].volume-a.chamber[i].volume)/delta;if(!std::isfinite((*volume_slope)[i])){return false;}}}}
        return true;
    }
    std::array<Pair,4> x;for(std::size_t i=0;i<4;++i){x[i]={first[i],last[i],i==coordinate?1.:0.};}
    const double h=thickness_[0],total=h+thickness_[1],c0=.5*h,c1=h+.5*thickness_[1];
    const std::array<double,5> nodes={0.,c0,h,c1,total};
    const auto area=[&](Pair z,std::size_t segment){
        if(segment==0){return x[2];}if(segment==3){return x[3];}
        return x[2]+(x[3]-x[2])*(z-Pair(c0))/Pair(c1-c0);
    };
    std::array<Pair,2> volume{},moment{},centroid{};
    for(std::size_t j=0;j<4;++j){const std::size_t cell=j<2?0:1;
        const Pair lo(nodes[j]),hi(nodes[j+1]),length=hi-lo,a=area(lo,j),b=area(hi,j);
        Pair v,local_moment;
        if(j==0 || j==3){const Pair opening=positive_part(a);v=length*opening;local_moment=Pair(.5)*length*length*opening;}
        else if(std::min({a.a,a.b,b.a,b.b})>=0){v=Pair(.5)*length*(a+b);local_moment=length*length*(a+Pair(2)*b)/Pair(6);}
        else if(std::max({a.a,a.b,b.a,b.b})<=0){continue;}
        else{const Pair pa=positive_part(a),pb=positive_part(b),difference=b-a,square=pb*pb-pa*pa;
            v=length*square/(Pair(2)*difference);
            local_moment=length*length*((pb*pb*pb-pa*pa*pa)/Pair(3)-a*square/Pair(2))/(difference*difference);}
        volume[cell]=volume[cell]+v;moment[cell]=moment[cell]+lo*v+local_moment;
    }
    for(std::size_t i=0;i<2;++i){if(!std::isfinite(volume[i].a) || !std::isfinite(volume[i].b) || volume[i].a<=0 || volume[i].b<=0){return false;}
        centroid[i]=moment[i]/volume[i];if(!std::isfinite(volume[i].d)){return false;}if(volume_slope){(*volume_slope)[i]=volume[i].d;}}
    const std::array<Pair,4> ends={Pair(0),centroid[0],centroid[1],Pair(total)};
    for(std::size_t p=0;p<3;++p){
        if((p==0 && first[2]<=0 && last[2]<=0) || (p==2 && first[3]<=0 && last[3]<=0)){value[p]=0.;slope[p]=0.;continue;}
        Pair inertia;
        for(std::size_t j=0;j<4;++j){const std::size_t cell=j<2?0:1;
            // Clamping both path endpoints lets a moving centroid cross a
            // profile knot without changing the finite work identity.
            const Pair lo=clamp(ends[p],nodes[j],nodes[j+1]),hi=clamp(ends[p+1],nodes[j],nodes[j+1]);
            if(lo.a==hi.a && lo.b==hi.b && lo.d==hi.d){continue;}
            inertia=inertia+Pair(density_)*volume[cell]/x[cell]*reciprocal_integral(hi-lo,area(lo,j),area(hi,j));
        }
        const Pair inverse=Pair(1)/inertia;value[p]=inverse.b;slope[p]=inverse.d;
        if(!std::isfinite(value[p]) || value[p]<=0 || !std::isfinite(slope[p])){return false;}
    }
    return true;
}
bool ContourKinetics::evaluate(const std::array<double,2>& mass,const std::array<double,2>& area,
    const std::array<double,3>& momentum,ContourKineticState& result) const noexcept {
    const std::array<double,4> x={mass[0],mass[1],area[0],area[1]};ContourKineticState candidate;
    if(!valid_coordinates(x)){return false;}for(double p:momentum){if(!std::isfinite(p)){return false;}}
    const GlottalContour profile{thickness_[0],thickness_[1],area[0],area[1]};if(!profile.geometry(candidate.geometry)){return false;}
    for(std::size_t j=0;j<3;++j){if(!candidate.geometry.passage[j].open && momentum[j]!=0){return false;}}
    for(std::size_t k=0;k<4;++k){std::array<double,3> gradient;
        if(!metric(x,x,k,candidate.inverse_inertance,gradient)){return false;}
        for(std::size_t j=0;j<3;++j){candidate.effort[k]+=.5*momentum[j]*momentum[j]*gradient[j];}if(!std::isfinite(candidate.effort[k])){return false;}
    }
    for(std::size_t j=0;j<3;++j){candidate.energy+=.5*momentum[j]*momentum[j]*candidate.inverse_inertance[j];
        for(std::size_t i=0;i<2;++i){const double inverse_density=candidate.geometry.chamber[i].volume/mass[i];
            candidate.resistance[j]+=12*viscosity_*width_*width_*inverse_density*inverse_density*candidate.geometry.passage[j].inverse_area_cubed[i];}
        if(!std::isfinite(candidate.resistance[j]) || !std::isfinite(candidate.energy)){return false;}
    }
    if(area[1]>0){const double outlet_speed=candidate.geometry.chamber[1].volume/(mass[1]*area[1]);
        candidate.quadratic_loss[2]=.5*outlet_loss_*density_*outlet_speed*outlet_speed;}
    if(!std::isfinite(candidate.quadratic_loss[2])){return false;}result=candidate;return true;
}
bool ContourKinetics::transfer(const std::array<double,2>& mass0,const std::array<double,2>& area0,const std::array<double,3>& p0,
    const std::array<double,2>& mass1,const std::array<double,2>& area1,const std::array<double,3>& p1,ContourKineticTransfer& result) const noexcept {
    const std::array<double,4> a={mass0[0],mass0[1],area0[0],area0[1]},b={mass1[0],mass1[1],area1[0],area1[1]};
    std::array<double,3> g0,g1,unused;if(!metric(a,a,4,g0,unused) || !metric(b,b,4,g1,unused)){return false;}
    ContourKineticTransfer candidate;double work=0;
    for(std::size_t j=0;j<3;++j){if(!std::isfinite(p0[j]) || !std::isfinite(p1[j]) || (g0[j]==0 && p0[j]!=0) || (g1[j]==0 && p1[j]!=0)){return false;}
        candidate.flow[j]=.25*(g0[j]+g1[j])*(p0[j]+p1[j]);work+=candidate.flow[j]*(p1[j]-p0[j]);
        candidate.energy_change+=.5*(p1[j]*p1[j]*g1[j]-p0[j]*p0[j]*g0[j]);}
    // Average the forward and reverse coordinate paths. This gives a
    // time-symmetric discrete gradient and its correct coincident limit.
    for(int direction:{1,-1}){auto current=a;
        for(std::size_t step=0;step<4;++step){const std::size_t k=direction==1?step:3-step;auto next=current;next[k]=b[k];std::array<double,3> gradient;
            std::array<double,2> volume_slope;if(!metric(current,next,k,unused,gradient,&volume_slope)){return false;}
            if(k>=2){for(std::size_t i=0;i<2;++i){candidate.volume_gradient[i][k-2]+=.5*volume_slope[i];}}
            for(std::size_t j=0;j<3;++j){candidate.effort[k]+=.125*(p0[j]*p0[j]+p1[j]*p1[j])*gradient[j];}current=next;
        }
    }
    for(std::size_t k=0;k<4;++k){work+=candidate.effort[k]*(b[k]-a[k]);}
    for(const auto& row:candidate.volume_gradient){for(double v:row){if(!std::isfinite(v)){return false;}}}
    candidate.chain_residual=candidate.energy_change-work;
    for(double v:candidate.effort){if(!std::isfinite(v)){return false;}}if(!std::isfinite(candidate.chain_residual)){return false;}result=candidate;return true;
}
}
