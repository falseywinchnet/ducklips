// stdin: nx ny solid_count depth, x edges, y edges, then domain and each
// solid as: vertex_count followed by x y pairs. All coordinates are metres.
#include "ducklips/planar_airspace.hpp"
#include "ducklips/planar_face_motion.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <exception>
#include <vector>
namespace {
void contour(ducklips::PlanarContour points){
    std::putchar('[');for(std::size_t i=0;i<points.size();++i){std::printf("%s[%.17g,%.17g]",i?",":"",points[i].x,points[i].y);}std::putchar(']');
}
}
int main(int argc,char** argv){
    const bool ports_mode=(argc==2 && std::string(argv[1])=="--ports") || (argc==4 && std::string(argv[3])=="--ports");
    const bool face_mode=argc==3 && (std::string(argv[1])=="--face-v" || std::string(argv[1])=="--face-h");
    const bool moving=(argc==3 || (argc==4 && ports_mode)) && (std::string(argv[1])=="--between" || face_mode);
    if(argc!=1 && !moving && !(argc==2 && ports_mode)){std::fprintf(stderr,"usage: planar_regions_probe [--between FRACTION [--ports] | --ports | --face-v INDEX | --face-h INDEX] < contours.txt\nports input appends count and domain-edge indices\n");return 2;}
    double fraction=0;
    if(moving){char* end=nullptr;fraction=std::strtod(argv[2],&end);if(end==argv[2] || *end){return 2;}}
    std::size_t nx=0,ny=0,count=0;double depth=0;
    if(std::scanf("%zu %zu %zu %lf",&nx,&ny,&count,&depth)!=4 || nx==0 || ny==0 || nx>256 || ny>256 || count>32){return 2;}
    std::vector<double> x(nx+1),y(ny+1);for(double& v:x){if(std::scanf("%lf",&v)!=1){return 2;}}for(double& v:y){if(std::scanf("%lf",&v)!=1){return 2;}}
    std::vector<std::vector<ducklips::Point2>> polygons(count+1);std::size_t edges=0;
    for(std::vector<ducklips::Point2>& polygon:polygons){
        std::size_t size=0;if(std::scanf("%zu",&size)!=1 || size<3 || size>4096){return 2;}edges+=size;if(edges>4096){return 2;}polygon.resize(size);
        for(ducklips::Point2& p:polygon){if(std::scanf("%lf %lf",&p.x,&p.y)!=2){return 2;}}
    }
    std::vector<ducklips::PlanarContour> solids;for(std::size_t i=1;i<polygons.size();++i){solids.emplace_back(polygons[i]);}
    std::vector<std::vector<ducklips::Point2>> last(count+1);std::vector<ducklips::PlanarContour> last_solids;
    if(moving){
        for(std::size_t i=0;i<last.size();++i){
            std::size_t size=0;if(std::scanf("%zu",&size)!=1 || size!=polygons[i].size()){return 2;}
            last[i].resize(size);for(ducklips::Point2& p:last[i]){if(std::scanf("%lf %lf",&p.x,&p.y)!=2){return 2;}}
            if(i){last_solids.emplace_back(last[i]);}
        }
    }
    std::vector<std::size_t> selected;
    if(ports_mode){
        std::size_t size=0;if(std::scanf("%zu",&size)!=1 || size>polygons[0].size()){return 2;}
        selected.resize(size);for(std::size_t& edge:selected){if(std::scanf("%zu",&edge)!=1){return 2;}}
    }
    if(face_mode){
        char* end=nullptr;const unsigned long index=std::strtoul(argv[2],&end,10);if(end==argv[2] || *end){return 2;}
        const bool vertical=std::string(argv[1])=="--face-v";ducklips::PlanarFace face;face.vertical=vertical;
        if(vertical){if(index>=(nx+1)*ny || index%(nx+1)==0 || index%(nx+1)==nx){return 2;}
            face.coordinate=x[index%(nx+1)];face.low=y[index/(nx+1)];face.high=y[index/(nx+1)+1];}
        else{if(index>=nx*(ny+1) || index/nx==0 || index/nx==ny){return 2;}
            face.coordinate=y[index/nx];face.low=x[index%nx];face.high=x[index%nx+1];}
        ducklips::PlanarFaceMotion motion;motion.prepare(edges,262144);
        if(!motion.build(polygons[0],solids,last[0],last_solids,face,depth)){std::fprintf(stderr,"face motion failed\n");return 1;}
        std::printf("{\"average_area\":%.17g,\"error_estimate\":%.17g,\"windows\":[",motion.average_area,motion.error_estimate);
        for(std::size_t i=0;i<motion.windows.size();++i){const auto& w=motion.windows[i];
            std::printf("%s[%.17g,%.17g,%.17g,%.17g]",i?",":"",w.first,w.last,w.average_area,w.error_estimate);}
        std::puts("]}");return 0;
    }
    ducklips::PlanarAirspace space;
    try{space.prepare(x,y,edges);space.prepare_regions(262144,524288);}catch(const std::exception& error){std::fprintf(stderr,"%s\n",error.what());return 2;}
    const bool okay=moving?space.build_between(polygons[0],solids,last[0],last_solids,fraction,depth,selected):space.build_regions(polygons[0],solids,depth,selected);
    if(!okay){std::fprintf(stderr,"region decomposition failed\n");return 1;}
    std::printf("{\"regions\":[");
    for(std::size_t i=0;i<space.regions.size();++i){const ducklips::PlanarRegion& r=space.regions[i];std::printf("%s[%zu,%.17g,%.17g,%.17g]",i?",":"",r.cell,r.volume,r.interior.x,r.interior.y);}
    std::printf("],\"links\":[");
    for(std::size_t i=0;i<space.region_links.size();++i){const ducklips::PlanarRegionLink& r=space.region_links[i];std::printf("%s[%zu,%zu,%zu,%d,%.17g]",i?",":"",r.first,r.second,r.face,int(r.vertical),r.area);}
    std::putchar(']');
    if(ports_mode){
        std::printf(",\"boundary_ports\":[");
        for(std::size_t i=0;i<space.boundary_ports.size();++i){const auto& p=space.boundary_ports[i];
            std::printf("%s[%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g]",i?",":"",p.edge,p.region,p.area,p.centre.x,p.centre.y,p.outward_normal.x,p.outward_normal.y);}
        std::putchar(']');
    }
    if(moving){
        const auto sampled=space.interpolated_contours();std::printf(",\"fraction\":%.17g,\"domain\":",fraction);contour(sampled[0]);
        std::printf(",\"solids\":[");for(std::size_t i=1;i<sampled.size();++i){if(i>1){std::putchar(',');}contour(sampled[i]);}std::putchar(']');
    }
    std::puts("}");
}
