#include "ducklips/planar_face_motion.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ducklips {
namespace {
long double evaluate(const long double* c,int degree,long double t) noexcept {
    long double value=c[degree];for(int i=degree-1;i>=0;--i){value=value*t+c[i];}return value;
}
// Derivative partitioning isolates every sign-changing root of degree <= 3.
// Near-zero stationary points are retained too (tangencies may have no sign
// change). Extra cuts are harmless; no distance tolerance merges nearby roots.
std::size_t polynomial_roots(const long double* c,int degree,std::array<long double,8>& output) noexcept {
    while(degree && c[degree]==0){--degree;}if(!degree){return 0;}
    if(degree==1){const long double t=-c[0]/c[1];if(t>0 && t<1){output[0]=t;return 1;}return 0;}
    std::array<long double,4> derivative{};for(int i=1;i<=degree;++i){derivative[std::size_t(i-1)]=i*c[i];}
    std::array<long double,8> stationary{};const std::size_t count=polynomial_roots(derivative.data(),degree-1,stationary);
    std::array<long double,10> cuts{};cuts[0]=0;for(std::size_t i=0;i<count;++i){cuts[i+1]=stationary[i];}cuts[count+1]=1;
    long double scale=0;for(int i=0;i<=degree;++i){scale+=std::abs(c[i]);}
    const long double roundoff=64*std::numeric_limits<long double>::epsilon()*scale;
    std::size_t found=0;
    for(std::size_t i=1;i<=count;++i){if(std::abs(evaluate(c,degree,cuts[i]))<=roundoff){output[found++]=cuts[i];}}
    for(std::size_t i=0;i<=count;++i){
        long double lo=cuts[i],hi=cuts[i+1],a=evaluate(c,degree,lo),b=evaluate(c,degree,hi);
        if(a==0 || b==0 || std::signbit(a)==std::signbit(b)){continue;}
        for(unsigned iteration=0;iteration<80;++iteration){
            const long double mid=lo+(hi-lo)/2;if(mid==lo || mid==hi){break;}
            const long double value=evaluate(c,degree,mid);if(value==0){lo=hi=mid;break;}
            if(std::signbit(value)==std::signbit(a)){lo=mid;a=value;}else{hi=mid;}
        }
        output[found++]=lo+(hi-lo)/2;
    }
    std::sort(output.begin(),output.begin()+static_cast<std::ptrdiff_t>(found));
    return static_cast<std::size_t>(std::unique(output.begin(),output.begin()+static_cast<std::ptrdiff_t>(found))-output.begin());
}
Point2 subtract(Point2 a,Point2 b) noexcept {return {a.x-b.x,a.y-b.y};}
Point2 local(Point2 p,PlanarFace f) noexcept {return f.vertical?Point2{p.x-f.coordinate,p.y}:Point2{p.y-f.coordinate,p.x};}
}
void PlanarFaceMotion::prepare(std::size_t maximum_edges,std::size_t maximum_events){
    if(maximum_edges<3 || maximum_edges>4096 || maximum_events<2){throw std::invalid_argument("invalid face motion capacity");}
    maximum_edges_=maximum_edges;maximum_events_=maximum_events;
    edges_.reserve(maximum_edges);crossings_.reserve(maximum_edges);events_.reserve(maximum_events);windows.reserve(maximum_events-1);
    hits_.reserve(maximum_edges);domain_left_.reserve(maximum_edges);domain_right_.reserve(maximum_edges);solids_.reserve(2*maximum_edges);
    windows.clear();average_area=error_estimate=0;
}
bool PlanarFaceMotion::add_event(long double time) noexcept {
    if(time<=0 || time>=1){return true;}if(events_.size()==maximum_events_){return false;}
    events_.push_back(static_cast<double>(time));return true;
}
bool PlanarFaceMotion::roots(const long double* polynomial,int degree) noexcept {
    std::array<long double,4> normalized{};long double scale=0;
    for(int i=0;i<=degree;++i){if(!std::isfinite(polynomial[i])){return false;}scale=std::max(scale,std::abs(polynomial[i]));}
    if(scale==0){return true;}for(int i=0;i<=degree;++i){normalized[std::size_t(i)]=polynomial[i]/scale;}
    std::array<long double,8> result{};const std::size_t count=polynomial_roots(normalized.data(),degree,result);
    for(std::size_t i=0;i<count;++i){if(!add_event(result[i])){return false;}}return true;
}
void PlanarFaceMotion::intervals(std::size_t contour,double t,bool positive,std::vector<Interval>& output) noexcept {
    hits_.clear();
    for(const Edge& e:edges_){
        if(e.contour!=contour){continue;}
        const Point2 a{e.a.x+t*e.da.x,e.a.y+t*e.da.y},b{e.b.x+t*e.db.x,e.b.y+t*e.db.y};
        // Limits from both sides of the face. Their intersection for domain,
        // union for solids seals a contour lying exactly along the face.
        const bool crossing=positive?((a.x<=0 && b.x>0)||(b.x<=0 && a.x>0))
                                    :((a.x<0 && b.x>=0)||(b.x<0 && a.x>=0));
        if(crossing){hits_.push_back(a.y-a.x*(b.y-a.y)/(b.x-a.x));}
    }
    std::sort(hits_.begin(),hits_.end());
    for(std::size_t i=0;i+1<hits_.size();i+=2){
        const double lo=std::max(face_.low,hits_[i]),hi=std::min(face_.high,hits_[i+1]);
        if(hi>lo){output.push_back({lo,hi});}
    }
}
double PlanarFaceMotion::area(double fraction) noexcept {
    domain_left_.clear();domain_right_.clear();solids_.clear();
    intervals(0,fraction,false,domain_left_);intervals(0,fraction,true,domain_right_);
    for(std::size_t contour=1;contour<contours_;++contour){intervals(contour,fraction,false,solids_);intervals(contour,fraction,true,solids_);}
    std::sort(solids_.begin(),solids_.end(),[](Interval a,Interval b){return a.low<b.low;});
    double length=0;
    for(Interval left:domain_left_){for(Interval right:domain_right_){
        double start=std::max(left.low,right.low);const double end=std::min(left.high,right.high);if(end<=start){continue;}
        for(Interval solid:solids_){
            if(solid.high<=start){continue;}if(solid.low>=end){break;}
            if(solid.low>start){length+=solid.low-start;}start=std::max(start,solid.high);if(start>=end){break;}
        }
        if(start<end){length+=end-start;}
    }}
    return depth_*length;
}
bool PlanarFaceMotion::integrate(double lo,double hi,double tolerance,unsigned level,double& integral,double& error) noexcept {
    // Independent 4/8 point Gauss estimates, used only BETWEEN event roots.
    constexpr std::array<double,4> x8{.18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
    constexpr std::array<double,4> w8{.36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
    constexpr std::array<double,2> x4{.33998104358485626480,.86113631159405257522};
    constexpr std::array<double,2> w4{.65214515486254614263,.34785484513745385737};
    const double half=(hi-lo)/2,mid=lo+half;double high=0,low=0;
    for(std::size_t i=0;i<x8.size();++i){high+=w8[i]*(area(mid-half*x8[i])+area(mid+half*x8[i]));}
    for(std::size_t i=0;i<x4.size();++i){low+=w4[i]*(area(mid-half*x4[i])+area(mid+half*x4[i]));}
    high*=half;low*=half;const double difference=std::abs(high-low);
    if(!std::isfinite(high) || !std::isfinite(low)){return false;}
    if(difference<=tolerance){integral=high;error=difference;return true;}
    if(level==20 || mid==lo || mid==hi){return false;}
    double first=0,last=0,a=0,b=0;
    if(!integrate(lo,mid,tolerance/2,level+1,first,a) || !integrate(mid,hi,tolerance/2,level+1,last,b)){return false;}
    integral=first+last;error=a+b;return true;
}
bool PlanarFaceMotion::build(PlanarContour first_domain,std::span<const PlanarContour> first_solids,
    PlanarContour last_domain,std::span<const PlanarContour> last_solids,PlanarFace face,double depth,double tolerance) noexcept {
    windows.clear();edges_.clear();crossings_.clear();events_.clear();average_area=error_estimate=0;
    if(!maximum_edges_ || !std::isfinite(face.coordinate) || !std::isfinite(face.low) || !std::isfinite(face.high) || face.high<=face.low
        || !std::isfinite(depth) || depth<=0 || !std::isfinite(tolerance) || tolerance<=0 || first_solids.size()!=last_solids.size()){return false;}
    face_=face;depth_=depth;contours_=first_solids.size()+1;events_.push_back(0);events_.push_back(1);
    for(std::size_t contour=0;contour<contours_;++contour){
        const PlanarContour first=contour?first_solids[contour-1]:first_domain,last=contour?last_solids[contour-1]:last_domain;
        if(first.size()<3 || first.size()!=last.size() || first.size()>maximum_edges_-edges_.size()){return false;}
        for(std::size_t i=0;i<first.size();++i){
            for(Point2 p:{first[i],last[i]}){if(!std::isfinite(p.x) || !std::isfinite(p.y)){return false;}}
        }
        for(std::size_t i=0;i<first.size();++i){
            const std::size_t j=(i+1)%first.size();
            const Point2 a=local(first[i],face),b=local(first[j],face),end_a=local(last[i],face),end_b=local(last[j],face);
            const Edge e{a,b,subtract(end_a,a),subtract(end_b,b),contour};
            for(Point2 p:{e.a,e.b,e.da,e.db}){if(!std::isfinite(p.x) || !std::isfinite(p.y)){return false;}}
            edges_.push_back(e);
            // Bounding boxes over linear endpoint paths safely discard edges
            // which never intersect the face's supporting line.
            if(std::min({a.x,b.x,end_a.x,end_b.x})>0 || std::max({a.x,b.x,end_a.x,end_b.x})<0
                || std::min({a.y,b.y,end_a.y,end_b.y})>face.high || std::max({a.y,b.y,end_a.y,end_b.y})<face.low){continue;}
            const long double vertex[2]{a.x,e.da.x};if(!roots(vertex,1)){return false;}
            const long double other_vertex[2]{b.x,e.db.x};if(!roots(other_vertex,1)){return false;}
            Crossing c;
            c.d[0]=static_cast<long double>(b.x)-a.x;c.d[1]=static_cast<long double>(e.db.x)-e.da.x;
            // Intersection's tangent coordinate = N(t)/D(t).
            c.n[0]=static_cast<long double>(a.y)*b.x-static_cast<long double>(b.y)*a.x;
            c.n[1]=static_cast<long double>(a.y)*e.db.x+static_cast<long double>(e.da.y)*b.x
                  -static_cast<long double>(b.y)*e.da.x-static_cast<long double>(e.db.y)*a.x;
            c.n[2]=static_cast<long double>(e.da.y)*e.db.x-static_cast<long double>(e.db.y)*e.da.x;
            if(!roots(c.d,1)){return false;}
            for(double bound:{face.low,face.high}){
                const long double polynomial[3]{c.n[0]-bound*c.d[0],c.n[1]-bound*c.d[1],c.n[2]};
                if(!roots(polynomial,2)){return false;}
            }
            crossings_.push_back(c);
        }
    }
    // An interval's endpoint can change owner when two edge intersections
    // exchange order. N1*D2-N2*D1 has degree <= 3. Include all roots, even
    // those outside the segments, to avoid a tolerance-based event filter.
    for(std::size_t i=0;i<crossings_.size();++i){for(std::size_t j=0;j<i;++j){
        const Crossing& a=crossings_[i];const Crossing& b=crossings_[j];std::array<long double,4> p{};
        for(std::size_t n=0;n<3;++n){for(std::size_t d=0;d<2;++d){p[n+d]+=a.n[n]*b.d[d]-b.n[n]*a.d[d];}}
        if(!roots(p.data(),3)){return false;}
    }}
    std::sort(events_.begin(),events_.end());events_.erase(std::unique(events_.begin(),events_.end()),events_.end());
    for(std::size_t i=1;i<events_.size();++i){
        const double first=events_[i-1],last=events_[i],width=last-first;double value=0,error=0;
        if(!integrate(first,last,tolerance*width,0,value,error)){return false;}
        windows.push_back({first,last,value/width,error/width});average_area+=value;error_estimate+=error;
    }
    return true;
}
}
