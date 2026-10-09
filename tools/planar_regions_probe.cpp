// stdin: nx ny solid_count depth, x edges, y edges, then domain and each
// solid as: vertex_count followed by x y pairs. All coordinates are metres.
#include "ducklips/planar_airspace.hpp"
#include <cstdio>
#include <exception>
#include <vector>
int main(){
    std::size_t nx=0,ny=0,count=0;double depth=0;
    if(std::scanf("%zu %zu %zu %lf",&nx,&ny,&count,&depth)!=4 || nx==0 || ny==0 || nx>256 || ny>256 || count>32){return 2;}
    std::vector<double> x(nx+1),y(ny+1);for(double& v:x){if(std::scanf("%lf",&v)!=1){return 2;}}for(double& v:y){if(std::scanf("%lf",&v)!=1){return 2;}}
    std::vector<std::vector<ducklips::Point2>> polygons(count+1);std::size_t edges=0;
    for(std::vector<ducklips::Point2>& polygon:polygons){
        std::size_t size=0;if(std::scanf("%zu",&size)!=1 || size<3 || size>4096){return 2;}edges+=size;if(edges>4096){return 2;}polygon.resize(size);
        for(ducklips::Point2& p:polygon){if(std::scanf("%lf %lf",&p.x,&p.y)!=2){return 2;}}
    }
    std::vector<ducklips::PlanarContour> solids;for(std::size_t i=1;i<polygons.size();++i){solids.emplace_back(polygons[i]);}
    ducklips::PlanarAirspace space;
    try{space.prepare(x,y,edges);space.prepare_regions(262144,524288);}catch(const std::exception& error){std::fprintf(stderr,"%s\n",error.what());return 2;}
    if(!space.build_regions(polygons[0],solids,depth)){std::fprintf(stderr,"region decomposition failed\n");return 1;}
    std::printf("{\"regions\":[");
    for(std::size_t i=0;i<space.regions.size();++i){const ducklips::PlanarRegion& r=space.regions[i];std::printf("%s[%zu,%.17g,%.17g,%.17g]",i?",":"",r.cell,r.volume,r.interior.x,r.interior.y);}
    std::printf("],\"links\":[");
    for(std::size_t i=0;i<space.region_links.size();++i){const ducklips::PlanarRegionLink& r=space.region_links[i];std::printf("%s[%zu,%zu,%zu,%d,%.17g]",i?",":"",r.first,r.second,r.face,int(r.vertical),r.area);}
    std::puts("]}");
}
