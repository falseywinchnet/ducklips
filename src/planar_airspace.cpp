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
    piece_capacity_=link_capacity_=0;regions.clear();region_links.clear();
    x_.assign(x.begin(),x.end());y_.assign(y.begin(),y.end());maximum_edges_=maximum_edges;
    volume.resize(columns()*rows());vertical_face.resize(x_.size()*rows());horizontal_face.resize(columns()*y_.size());
    edges_.reserve(maximum_edges);hits_.reserve(maximum_edges);
    cuts_.reserve(x.size()+maximum_edges*(y.size()+1)+maximum_edges*(maximum_edges-1)/2);
    domain_.reserve(maximum_edges);solid_.reserve(maximum_edges);free_.reserve(2*maximum_edges);left_.reserve(2*maximum_edges);
}
void PlanarAirspace::polygon_intervals(PlanarContour polygon,std::size_t first_edge,double x,std::vector<Interval>& out) noexcept {
    hits_.clear();
    for(std::size_t i=0;i<polygon.size();++i){
        const Point2 a=polygon[i],b=polygon[(i+1)%polygon.size()];
        // Half-open crossing avoids double-counting polygon vertices.
        if((a.x<=x && x<b.x) || (b.x<=x && x<a.x)){
            hits_.push_back({a.y+(x-a.x)*(b.y-a.y)/(b.x-a.x),first_edge+i});
        }
    }
    std::sort(hits_.begin(),hits_.end(),[](Hit a,Hit b){return a.y<b.y;});
    for(std::size_t i=0;i+1<hits_.size();i+=2){if(hits_[i+1].y>hits_[i].y){out.push_back({hits_[i].y,hits_[i+1].y,hits_[i].edge,hits_[i+1].edge});}}
}
void PlanarAirspace::free_intervals(PlanarContour domain,std::span<const PlanarContour> solids,double x) noexcept {
    domain_.clear();solid_.clear();free_.clear();polygon_intervals(domain,0,x,domain_);
    std::size_t offset=domain.size();
    for(PlanarContour solid:solids){polygon_intervals(solid,offset,x,solid_);offset+=solid.size();}
    std::sort(solid_.begin(),solid_.end(),[](Interval a,Interval b){return a.low<b.low;});
    // Subtract the solid union from every connected slice of the domain.
    for(Interval air:domain_){
        double start=air.low;std::size_t start_edge=air.low_edge;
        for(Interval solid:solid_){
            if(solid.high<=start){continue;}
            if(solid.low>=air.high){break;}
            if(solid.low>start){free_.push_back({start,solid.low,start_edge,solid.low_edge});}
            if(solid.high>start){start=solid.high;start_edge=solid.high_edge;}if(start>=air.high){break;}
        }
        if(start<air.high){free_.push_back({start,air.high,start_edge,air.high_edge});}
    }
}
bool PlanarAirspace::build(PlanarContour domain,std::span<const PlanarContour> solids,double depth) noexcept {
    if(x_.empty() || !std::isfinite(depth) || depth<=0 || domain.size()<3){return false;}
    std::size_t count=domain.size();
    for(PlanarContour p:solids){if(p.size()<3 || p.size()>maximum_edges_-std::min(count,maximum_edges_)){return false;}count+=p.size();}
    if(count>maximum_edges_){return false;}
    regions.clear();region_links.clear();
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
void PlanarAirspace::prepare_regions(std::size_t piece_capacity,std::size_t link_capacity){
    if(x_.empty() || piece_capacity==0 || link_capacity==0){throw std::invalid_argument("prepare the grid and positive region capacities first");}
    piece_capacity_=piece_capacity;link_capacity_=link_capacity;
    pieces_.reserve(piece_capacity);regions.reserve(piece_capacity);region_links.reserve(link_capacity);
    previous_rows_.resize(rows());current_rows_.resize(rows());
}
std::size_t PlanarAirspace::piece_root(std::size_t piece) noexcept {
    while(pieces_[piece].parent!=piece){pieces_[piece].parent=pieces_[pieces_[piece].parent].parent;piece=pieces_[piece].parent;}
    return piece;
}
double PlanarAirspace::boundary_y(std::size_t edge,double x) const noexcept {
    const Segment& s=edges_[edge];
    if(x==s.a.x){return s.a.y;}if(x==s.b.x){return s.b.y;}
    return s.a.y+(x-s.a.x)*(s.b.y-s.a.y)/(s.b.x-s.a.x);
}
bool PlanarAirspace::build_regions(PlanarContour domain,std::span<const PlanarContour> solids,double depth) noexcept {
    if(!piece_capacity_ || !build(domain,solids,depth)){return false;}
    pieces_.clear();regions.clear();region_links.clear();
    std::fill(previous_rows_.begin(),previous_rows_.end(),Range{});
    std::size_t previous_column=0;bool have_previous=false;
    for(std::size_t cut=1;cut<cuts_.size();++cut){
        const double lo=cuts_[cut-1],hi=cuts_[cut],x=lo+(hi-lo)/2;if(x<=lo || x>=hi){continue;}
        const std::size_t column=static_cast<std::size_t>(std::upper_bound(x_.begin(),x_.end(),x)-x_.begin()-1);
        free_intervals(domain,solids,x);
        for(std::size_t row=0;row<rows();++row){
            current_rows_[row].begin=pieces_.size();
            for(const Interval& air:free_){
                const double low=std::max(air.low,y_[row]),high=std::min(air.high,y_[row+1]);if(high<=low){continue;}
                if(pieces_.size()==piece_capacity_){return false;}
                const double ll=std::max(boundary_y(air.low_edge,lo),y_[row]),hl=std::min(boundary_y(air.high_edge,lo),y_[row+1]);
                const double lr=std::max(boundary_y(air.low_edge,hi),y_[row]),hr=std::min(boundary_y(air.high_edge,hi),y_[row+1]);
                const std::size_t slice=static_cast<std::size_t>(&air-free_.data());
                const std::size_t id=pieces_.size();pieces_.push_back({row*columns()+column,id,0,slice,ll,hl,lr,hr,depth*(hi-lo)*(high-low),{x,.5*(low+high)}});
                // Adjacent strips in one cell form a region only through an
                // opening of positive length. Across a grid face, retain a port.
                if(have_previous){for(std::size_t other=previous_rows_[row].begin;other<previous_rows_[row].end;++other){
                    const Piece& a=pieces_[other];const double opening=std::min(a.high_right,hl)-std::max(a.low_right,ll);
                    if(opening<=0){continue;}
                    if(column==previous_column){pieces_[piece_root(id)].parent=piece_root(other);}
                    else{
                        if(region_links.size()==link_capacity_){return false;}
                        region_links.push_back({other,id,row*x_.size()+column,true,depth*opening});
                    }
                }}
                // The same free slice crosses a horizontal cell boundary iff
                // its open interval strictly contains that boundary.
                if(row && air.low<y_[row] && air.high>y_[row]){
                    for(std::size_t other=current_rows_[row-1].begin;other<current_rows_[row-1].end;++other){
                        const Piece& a=pieces_[other];
                        if(a.slice==slice){
                            if(region_links.size()==link_capacity_){return false;}
                            region_links.push_back({other,id,row*columns()+column,false,depth*(hi-lo)});
                        }
                    }
                }
            }
            current_rows_[row].end=pieces_.size();
        }
        previous_rows_.swap(current_rows_);previous_column=column;have_previous=true;
    }
    constexpr std::size_t missing=static_cast<std::size_t>(-1);
    for(Piece& piece:pieces_){piece.region=missing;}
    for(std::size_t i=0;i<pieces_.size();++i){
        Piece& root=pieces_[piece_root(i)];
        if(root.region==missing){root.region=regions.size();regions.push_back({root.cell,0,root.interior});}
        regions[root.region].volume+=pieces_[i].volume;
    }
    for(PlanarRegionLink& link:region_links){
        link.first=pieces_[piece_root(link.first)].region;link.second=pieces_[piece_root(link.second)].region;
        if(link.first>link.second){std::swap(link.first,link.second);}
    }
    std::sort(region_links.begin(),region_links.end(),[](const PlanarRegionLink& a,const PlanarRegionLink& b){
        if(a.first!=b.first){return a.first<b.first;}if(a.second!=b.second){return a.second<b.second;}
        if(a.vertical!=b.vertical){return a.vertical<b.vertical;}return a.face<b.face;
    });
    std::size_t count=0;
    for(const PlanarRegionLink& link:region_links){
        if(count && region_links[count-1].first==link.first && region_links[count-1].second==link.second
            && region_links[count-1].vertical==link.vertical && region_links[count-1].face==link.face){region_links[count-1].area+=link.area;}
        else{region_links[count++]=link;}
    }
    region_links.resize(count);return true;
}

}
