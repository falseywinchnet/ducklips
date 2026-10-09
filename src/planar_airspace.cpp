#include "ducklips/planar_airspace.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ducklips {
namespace {
double cross(Point2 a,Point2 b) noexcept {return a.x*b.y-a.y*b.x;}
Point2 minus(Point2 a,Point2 b) noexcept {return {a.x-b.x,a.y-b.y};}
bool increasing(std::span<const double> values) {
    if(values.size()<2){return false;}
    for(std::size_t i=0;i<values.size();++i){if(!std::isfinite(values[i]) || (i && values[i]<=values[i-1])){return false;}}
    return true;
}
}
void PlanarAirspace::prepare(std::span<const double> x,std::span<const double> y,std::size_t maximum_edges){
    if(!increasing(x) || !increasing(y) || maximum_edges<3 || maximum_edges>4096){throw std::invalid_argument("invalid planar grid or edge capacity");}
    x_.assign(x.begin(),x.end());y_.assign(y.begin(),y.end());maximum_edges_=maximum_edges;
    volume.resize(columns()*rows());vertical_face.resize(x_.size()*rows());horizontal_face.resize(columns()*y_.size());
    edges_.reserve(maximum_edges);hits_.reserve(maximum_edges);
    cuts_.reserve(x.size()+maximum_edges*(y.size()+1)+maximum_edges*(maximum_edges-1)/2);
    domain_.reserve(maximum_edges);solid_.reserve(maximum_edges);free_.reserve(2*maximum_edges);left_.reserve(2*maximum_edges);
}
void PlanarAirspace::polygon_intervals(PlanarContour polygon,double x,std::vector<Interval>& out) noexcept {
    hits_.clear();
    for(std::size_t i=0;i<polygon.size();++i){
        const Point2 a=polygon[i],b=polygon[(i+1)%polygon.size()];
        // Half-open crossing avoids double-counting polygon vertices.
        if((a.x<=x && x<b.x) || (b.x<=x && x<a.x)){
            hits_.push_back(a.y+(x-a.x)*(b.y-a.y)/(b.x-a.x));
        }
    }
    std::sort(hits_.begin(),hits_.end());
    for(std::size_t i=0;i+1<hits_.size();i+=2){if(hits_[i+1]>hits_[i]){out.push_back({hits_[i],hits_[i+1]});}}
}
void PlanarAirspace::free_intervals(PlanarContour domain,std::span<const PlanarContour> solids,double x) noexcept {
    domain_.clear();solid_.clear();free_.clear();polygon_intervals(domain,x,domain_);
    for(PlanarContour solid:solids){polygon_intervals(solid,x,solid_);}
    std::sort(solid_.begin(),solid_.end(),[](Interval a,Interval b){return a.low<b.low;});
    // Subtract the solid union from every connected slice of the domain.
    for(Interval air:domain_){
        double start=air.low;
        for(Interval solid:solid_){
            if(solid.high<=start){continue;}
            if(solid.low>=air.high){break;}
            if(solid.low>start){free_.push_back({start,std::min(solid.low,air.high)});}
            start=std::max(start,solid.high);if(start>=air.high){break;}
        }
        if(start<air.high){free_.push_back({start,air.high});}
    }
}
bool PlanarAirspace::build(PlanarContour domain,std::span<const PlanarContour> solids,double depth) noexcept {
    if(x_.empty() || !std::isfinite(depth) || depth<=0 || domain.size()<3){return false;}
    std::size_t count=domain.size();
    for(PlanarContour p:solids){if(p.size()<3 || p.size()>maximum_edges_-std::min(count,maximum_edges_)){return false;}count+=p.size();}
    if(count>maximum_edges_){return false;}
    edges_.clear();cuts_.clear();
    for(double value:x_){cuts_.push_back(value);}
    for(std::size_t contour=0;contour<=solids.size();++contour){
        const PlanarContour p=contour==0?domain:solids[contour-1];
        for(std::size_t i=0;i<p.size();++i){
            const Point2 a=p[i],b=p[(i+1)%p.size()];
            if(!std::isfinite(a.x) || !std::isfinite(a.y)){return false;}
            edges_.push_back({a,b});
            if(a.x>x_.front() && a.x<x_.back()){cuts_.push_back(a.x);}
            if(a.y!=b.y){for(double y:y_){
                const double t=(y-a.y)/(b.y-a.y);
                if(t>0 && t<1){const double x=a.x+t*(b.x-a.x);if(x>x_.front() && x<x_.back()){cuts_.push_back(x);}}
            }}
        }
    }
    // Union ordering can change where edges of different contours intersect.
    // Include all intersections; shared endpoints and collinear overlaps are
    // already split by the vertex cuts. No geometric raster sampling is used.
    for(std::size_t i=0;i<edges_.size();++i){for(std::size_t j=0;j<i;++j){
        const Segment a=edges_[i],b=edges_[j];const Point2 u=minus(a.b,a.a),v=minus(b.b,b.a),w=minus(b.a,a.a);
        const double determinant=cross(u,v);if(determinant==0){continue;}
        const double t=cross(w,v)/determinant,s=cross(w,u)/determinant;
        if(t>0 && t<1 && s>0 && s<1){const double x=a.a.x+t*u.x;if(x>x_.front() && x<x_.back()){cuts_.push_back(x);}}
    }}
    std::sort(cuts_.begin(),cuts_.end());cuts_.erase(std::unique(cuts_.begin(),cuts_.end()),cuts_.end());
    std::fill(volume.begin(),volume.end(),0);std::fill(vertical_face.begin(),vertical_face.end(),0);std::fill(horizontal_face.begin(),horizontal_face.end(),0);
    for(std::size_t cut=1;cut<cuts_.size();++cut){
        const double lo=cuts_[cut-1],hi=cuts_[cut],x=lo+(hi-lo)/2;if(x<=lo || x>=hi){continue;}
        const std::size_t column=static_cast<std::size_t>(std::upper_bound(x_.begin(),x_.end(),x)-x_.begin()-1);
        free_intervals(domain,solids,x);
        for(Interval air:free_){
            for(std::size_t row=0;row<rows();++row){
                const double height=std::max(0.,std::min(air.high,y_[row+1])-std::max(air.low,y_[row]));
                volume[row*columns()+column]+=depth*(hi-lo)*height;
            }
            for(std::size_t face=1;face<rows();++face){
                if(y_[face]>air.low && y_[face]<air.high){horizontal_face[face*columns()+column]+=depth*(hi-lo);}
            }
        }
    }
    // A shared face is open only where fluid reaches it from BOTH sides.
    // This seals a solid edge exactly even when it lies on the grid face.
    for(std::size_t face=1;face<columns();++face){
        const double x=x_[face];
        free_intervals(domain,solids,std::nextafter(x,-std::numeric_limits<double>::infinity()));left_=free_;
        free_intervals(domain,solids,std::nextafter(x,std::numeric_limits<double>::infinity()));
        for(Interval a:left_){for(Interval b:free_){
            const double low=std::max(a.low,b.low),high=std::min(a.high,b.high);if(high<=low){continue;}
            for(std::size_t row=0;row<rows();++row){
                vertical_face[row*x_.size()+face]+=depth*std::max(0.,std::min(high,y_[row+1])-std::max(low,y_[row]));
            }
        }}
    }
    return true;
}
}
