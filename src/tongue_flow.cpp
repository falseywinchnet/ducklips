#include "ducklips/tongue_flow.hpp"
#include <algorithm>
#include <cmath>

namespace ducklips {
bool ConservativeTongueMap::configure(std::span<const TongueLiftKnot> knots,double anchor_y,double free_y,double centre_x) noexcept {
    count_=0;
    if (knots.size()<2 || knots.size()>maximum_knots || !std::isfinite(anchor_y)
        || !std::isfinite(free_y) || !std::isfinite(centre_x) || free_y<=anchor_y
        || !std::isfinite(free_y-anchor_y)) { return false; }
    if (knots.front().lift!=0 || knots.front().slope!=0 || knots.back().lift!=0 || knots.back().slope!=0) { return false; }
    double integral=0;
    for (std::size_t i=0;i+1<knots.size();++i) {
        const TongueLiftKnot& left=knots[i];const TongueLiftKnot& right=knots[i+1];
        const double width=right.x-left.x;
        if (!std::isfinite(left.x) || !std::isfinite(left.lift) || !std::isfinite(left.slope)
            || !std::isfinite(right.x) || !std::isfinite(right.lift) || !std::isfinite(right.slope) || width<=0) { return false; }
        Piece& p=pieces_[i];p.x=left.x;p.width=width;p.a=left.lift;p.b=width*left.slope;
        p.c=3*(right.lift-left.lift)-width*(2*left.slope+right.slope);
        p.d=2*(left.lift-right.lift)+width*(left.slope+right.slope);p.integral=integral;
        integral+=width*(p.a+p.b/2+p.c/3+p.d/4);
        if (!std::isfinite(integral+p.a+p.b+p.c+p.d)) { return false; }
    }
    count_=knots.size()-1;anchor_y_=anchor_y;height_=free_y-anchor_y;total_integral_=integral;
    double value=0,derivative=0;profile(centre_x,value,derivative,primitive_offset_);
    return true;
}
void ConservativeTongueMap::profile(double x,double& value,double& derivative,double& integral) const noexcept {
    value=0;derivative=0;integral=0;
    if (!count_ || x<=pieces_[0].x) { return; }
    const Piece& last=pieces_[count_-1];
    if (x>=last.x+last.width) { integral=total_integral_;return; }
    std::size_t low=0,high=count_;
    while (low+1<high) { const std::size_t mid=(low+high)/2;if (pieces_[mid].x<=x) { low=mid; } else { high=mid; } }
    const Piece& p=pieces_[low];const double t=(x-p.x)/p.width;
    value=p.a+t*(p.b+t*(p.c+t*p.d));derivative=(p.b+t*(2*p.c+3*t*p.d))/p.width;
    integral=p.integral+p.width*t*(p.a+t*(p.b/2+t*(p.c/3+t*p.d/4)));
}
double ConservativeTongueMap::lift(double x) const noexcept {
    double value=0,derivative=0,integral=0;profile(x,value,derivative,integral);return value;
}
ConservativeTongueMap::Field ConservativeTongueMap::field(Point2 point) const noexcept {
    double value=0,derivative=0,integral=0;profile(point.x,value,derivative,integral);integral-=primitive_offset_;
    const double t=std::clamp((point.y-anchor_y_)/height_,0.0,1.0);
    const double gate=t*t*t*(10+t*(-15+6*t));
    const double first=30*t*t*(1-t)*(1-t)/height_;
    const double second=60*t*(1-t)*(1-2*t)/(height_*height_);
    return Field{.velocity={-integral*first,value*gate},.derivative={-value*first,-integral*second,derivative*gate,value*first}};
}
TongueMaterialPoint ConservativeTongueMap::map(Point2 point,double command,std::size_t steps) const noexcept {
    TongueMaterialPoint result{.position=point};
    if (!count_ || !steps || !std::isfinite(command) || !std::isfinite(point.x) || !std::isfinite(point.y)) { return result; }
    // The anchored region does not move. Above the gate, x is fixed and y
    // follows a vertical shear exactly. Integrating either region repeatedly
    // adds no mechanics; retain its exact material map and unit determinant.
    if(command==0 || point.y<=anchor_y_){result.valid=true;return result;}
    if(point.y>=anchor_y_+height_){
        double value=0,derivative=0,integral=0;profile(point.x,value,derivative,integral);
        const double end=point.y+command*value;
        if(std::isfinite(end) && end>=anchor_y_+height_){
            result.position.y=end;result.gradient[2]=command*derivative;result.valid=true;return result;
        }
    }
    const double h=command/static_cast<double>(steps),half=h/2;
    for (std::size_t step=0;step<steps;++step) {
        const Point2 initial=result.position;const Field start=field(initial);
        Point2 next{initial.x+h*start.velocity.x,initial.y+h*start.velocity.y};bool converged=false;
        for (std::size_t iteration=0;iteration<24;++iteration) {
            const Field mid=field({(initial.x+next.x)/2,(initial.y+next.y)/2});
            const std::array<double,4>& a=mid.derivative;
            const double b00=1-half*a[0],b01=-half*a[1],b10=-half*a[2],b11=1-half*a[3];
            const double determinant=b00*b11-b01*b10;
            if (!std::isfinite(determinant) || std::abs(determinant)<1e-12) { return result; }
            const double rx=next.x-initial.x-h*mid.velocity.x,ry=next.y-initial.y-h*mid.velocity.y;
            const double dx=(b11*rx-b01*ry)/determinant,dy=(-b10*rx+b00*ry)/determinant;
            next.x-=dx;next.y-=dy;
            if (!std::isfinite(next.x) || !std::isfinite(next.y)) { return result; }
            if (std::max(std::abs(dx),std::abs(dy))<1e-12) { converged=true;break; }
        }
        if (!converged) { return result; }
        const Field mid=field({(initial.x+next.x)/2,(initial.y+next.y)/2});const std::array<double,4>& a=mid.derivative;
        const double b00=1-half*a[0],b01=-half*a[1],b10=-half*a[2],b11=1-half*a[3];
        const double determinant=b00*b11-b01*b10;
        if (!std::isfinite(determinant) || std::abs(determinant)<1e-12) { return result; }
        const std::array<double,4> m{
            (b11*(1+half*a[0])-b01*half*a[2])/determinant,
            (b11*half*a[1]-b01*(1+half*a[3]))/determinant,
            (-b10*(1+half*a[0])+b00*half*a[2])/determinant,
            (-b10*half*a[1]+b00*(1+half*a[3]))/determinant};
        const std::array<double,4> old=result.gradient;
        result.gradient={m[0]*old[0]+m[1]*old[2],m[0]*old[1]+m[1]*old[3],m[2]*old[0]+m[3]*old[2],m[2]*old[1]+m[3]*old[3]};
        result.position=next;
        for (double entry:result.gradient) { if (!std::isfinite(entry)) { return result; } }
    }
    result.valid=true;return result;
}

bool TongueMaterialTranslation::configure(Point2 start,Point2 displacement,double radius) noexcept {
    radius_squared_=0;
    if (!std::isfinite(start.x) || !std::isfinite(start.y) || !std::isfinite(displacement.x)
        || !std::isfinite(displacement.y) || !std::isfinite(radius) || radius<=0
        || !std::isfinite(radius*radius) || radius*radius==0) { return false; }
    start_=start;displacement_=displacement;radius_squared_=radius*radius;return true;
}
TongueMaterialTranslation::Field TongueMaterialTranslation::field(Point2 point,double time) const noexcept {
    const double x=point.x-start_.x-time*displacement_.x,y=point.y-start_.y-time*displacement_.y;
    const double q=(x*x+y*y)/radius_squared_;
    if (q>=1) { return {}; }
    // psi=(vx*y-vy*x)*(1-r^2/R^2)^4, velocity=(psi_y,-psi_x).
    // The value and its first three derivatives vanish at the support edge.
    const double s=1-q,w=s*s*s*s,first=-8*s*s*s/radius_squared_;
    const double second=48*s*s/(radius_squared_*radius_squared_);
    const double wx=first*x,wy=first*y,wxx=first+second*x*x,wyy=first+second*y*y,wxy=second*x*y;
    const double vx=displacement_.x,vy=displacement_.y,linear=vx*y-vy*x;
    const double a00=vx*wx-vy*wy+linear*wxy;
    return Field{.velocity={vx*w+linear*wy,vy*w-linear*wx},
                 .derivative={a00,2*vx*wy+linear*wyy,2*vy*wx-linear*wxx,-a00}};
}
TongueMaterialPoint TongueMaterialTranslation::map(Point2 point,bool inverse,std::size_t steps) const noexcept {
    TongueMaterialPoint result{.position=point};
    if (radius_squared_<=0 || !steps || !std::isfinite(point.x) || !std::isfinite(point.y)) { return result; }
    // Outside the entire swept support, the field and its gradient vanish
    // for the whole command interval. A stationary point cannot enter it.
    const double length_squared=displacement_.x*displacement_.x+displacement_.y*displacement_.y;
    const double fraction=length_squared>0 ? std::clamp(((point.x-start_.x)*displacement_.x+(point.y-start_.y)*displacement_.y)/length_squared,0.,1.) : 0.;
    const double outside_x=point.x-start_.x-fraction*displacement_.x,outside_y=point.y-start_.y-fraction*displacement_.y;
    if(length_squared==0 || outside_x*outside_x+outside_y*outside_y>=radius_squared_){result.valid=true;return result;}
    const double h=(inverse?-1.:1.)/static_cast<double>(steps),half=h/2;
    for (std::size_t step=0;step<steps;++step) {
        const double time=(inverse?1.:0.)+(static_cast<double>(step)+.5)*h;
        const Point2 initial=result.position;const Field start=field(initial,time);
        Point2 next{initial.x+h*start.velocity.x,initial.y+h*start.velocity.y};bool converged=false;
        for (std::size_t iteration=0;iteration<24;++iteration) {
            const Field mid=field({(initial.x+next.x)/2,(initial.y+next.y)/2},time);
            const std::array<double,4>& a=mid.derivative;
            const double b00=1-half*a[0],b01=-half*a[1],b10=-half*a[2],b11=1-half*a[3],det=b00*b11-b01*b10;
            if (!std::isfinite(det) || std::abs(det)<1e-12) { return result; }
            const double rx=next.x-initial.x-h*mid.velocity.x,ry=next.y-initial.y-h*mid.velocity.y;
            const double dx=(b11*rx-b01*ry)/det,dy=(-b10*rx+b00*ry)/det;
            next.x-=dx;next.y-=dy;
            if (!std::isfinite(next.x) || !std::isfinite(next.y)) { return result; }
            if (std::max(std::abs(dx),std::abs(dy))<1e-12) { converged=true;break; }
        }
        if (!converged) { return result; }
        const Field mid=field({(initial.x+next.x)/2,(initial.y+next.y)/2},time);const std::array<double,4>& a=mid.derivative;
        const double b00=1-half*a[0],b01=-half*a[1],b10=-half*a[2],b11=1-half*a[3],det=b00*b11-b01*b10;
        if (!std::isfinite(det) || std::abs(det)<1e-12) { return result; }
        const std::array<double,4> m{(b11*(1+half*a[0])-b01*half*a[2])/det,(b11*half*a[1]-b01*(1+half*a[3]))/det,
            (-b10*(1+half*a[0])+b00*half*a[2])/det,(-b10*half*a[1]+b00*(1+half*a[3]))/det};
        const std::array<double,4> old=result.gradient;
        result.gradient={m[0]*old[0]+m[1]*old[2],m[0]*old[1]+m[1]*old[3],m[2]*old[0]+m[3]*old[2],m[2]*old[1]+m[3]*old[3]};
        for (double entry:result.gradient) { if (!std::isfinite(entry)) { return result; } }
        result.position=next;
    }
    result.valid=true;return result;
}
}
