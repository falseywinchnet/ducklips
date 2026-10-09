// Native moving-geometry/interval-opening/pressure reference on a declared
// persistent region graph. The caller must independently establish identity;
// index/cell agreement alone is not a geometric continuation certificate.
// stdin: nx ny solids depth duration subdivisions; x edges; y edges;
// first then last contours (domain first; vertex_count, x y ...);
// node_count edge_count; initial pressures; edges: from to face vertical(0/1).
#include "ducklips/planar_airspace.hpp"
#include "ducklips/planar_face_motion.hpp"
#include "ducklips/cavity_network.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>
namespace {
using ducklips::Point2;using ducklips::PlanarContour;
using Polygons=std::vector<std::vector<Point2>>;
struct Port {ducklips::CavityConnection connection{};ducklips::PlanarFace face{};std::size_t index{};double distance{};};
void sample(const Polygons& first,const Polygons& last,double t,Polygons& result){
    for(std::size_t c=0;c<first.size();++c){for(std::size_t v=0;v<first[c].size();++v){
        result[c][v]={std::lerp(first[c][v].x,last[c][v].x,t),std::lerp(first[c][v].y,last[c][v].y,t)};
    }}
}
std::vector<PlanarContour> views(const Polygons& p){std::vector<PlanarContour> r;for(const auto& c:p){r.emplace_back(c);}return r;}
void values(const char* name,std::span<const double> a){
    std::printf(",\"%s\":[",name);for(std::size_t i=0;i<a.size();++i){std::printf("%s%.17g",i?",":"",a[i]);}std::putchar(']');
}
bool compatible(const ducklips::PlanarAirspace& geometry,std::span<const std::size_t> cells,std::span<const Port> ports){
    if(geometry.regions.size()!=cells.size()){return false;}
    for(std::size_t i=0;i<cells.size();++i){if(geometry.regions[i].cell!=cells[i]){return false;}}
    for(const auto& link:geometry.region_links){
        std::size_t count=0;for(const auto& p:ports){
            if(link.first==p.connection.from && link.second==p.connection.to && link.face==p.index && link.vertical==p.face.vertical){++count;}
        }
        if(count!=1){return false;}
    }
    // This reference permits at most one region pair on each supplied face.
    // Otherwise a whole-face integral cannot be assigned to that pair.
    for(std::size_t i=0;i<ports.size();++i){for(std::size_t j=0;j<i;++j){
        if(ports[i].index==ports[j].index && ports[i].face.vertical==ports[j].face.vertical){return false;}
    }}
    return true;
}
void contour(PlanarContour points){
    std::putchar('[');for(std::size_t i=0;i<points.size();++i){std::printf("%s[%.17g,%.17g]",i?",":"",points[i].x,points[i].y);}std::putchar(']');
}
void state(double fraction,const ducklips::PlanarAirspace& geometry,const ducklips::CavityNetwork& pressure,
           std::span<const PlanarContour> contours,std::span<const double> area,double residual,bool comma){
    std::printf("%s{\"fraction\":%.17g,\"regions\":[",comma?",":"",fraction);
    for(std::size_t i=0;i<geometry.regions.size();++i){const auto& r=geometry.regions[i];
        std::printf("%s[%zu,%.17g,%.17g,%.17g]",i?",":"",r.cell,r.volume,r.interior.x,r.interior.y);
    }
    std::printf("],\"links\":[");
    for(std::size_t i=0;i<geometry.region_links.size();++i){const auto& r=geometry.region_links[i];
        std::printf("%s[%zu,%zu,%zu,%d,%.17g]",i?",":"",r.first,r.second,r.face,int(r.vertical),r.area);
    }
    std::putchar(']');values("pressure",pressure.pressures());values("flow",pressure.flows());values("interval_area",area);
    std::printf(",\"domain\":");contour(contours[0]);std::printf(",\"solids\":[");
    for(std::size_t i=1;i<contours.size();++i){if(i>1){std::putchar(',');}contour(contours[i]);}std::putchar(']');
    std::printf(",\"maximum_mass_residual\":%.17g}",residual);
}
}
int main(){
    std::size_t nx=0,ny=0,solid_count=0,subdivisions=0;double depth=0,duration=0;
    if(std::scanf("%zu %zu %zu %lf %lf %zu",&nx,&ny,&solid_count,&depth,&duration,&subdivisions)!=6
        || !nx || !ny || nx>256 || ny>256 || solid_count>32 || !subdivisions || subdivisions>2048
        || !std::isfinite(duration) || duration<=0){return 2;}
    std::vector<double> x(nx+1),y(ny+1);
    for(double& v:x){if(std::scanf("%lf",&v)!=1){return 2;}}for(double& v:y){if(std::scanf("%lf",&v)!=1){return 2;}}
    Polygons first(solid_count+1),last(solid_count+1);std::size_t edges=0;
    for(std::size_t endpoint=0;endpoint<2;++endpoint){Polygons& polygons=endpoint?last:first;
        for(std::size_t c=0;c<polygons.size();++c){
            std::size_t count=0;if(std::scanf("%zu",&count)!=1 || count<3 || count>4096){return 2;}
            if(endpoint){if(count!=first[c].size()){return 2;}}else{edges+=count;if(edges>4096){return 2;}}
            polygons[c].resize(count);for(Point2& p:polygons[c]){if(std::scanf("%lf %lf",&p.x,&p.y)!=2){return 2;}}
        }
    }
    std::size_t nodes=0,port_count=0;if(std::scanf("%zu %zu",&nodes,&port_count)!=2 || !nodes || nodes>65536 || port_count>262144){return 2;}
    std::vector<double> initial(nodes),volume(nodes),compliance(nodes),inflow(nodes,0),areas(port_count);
    for(double& p:initial){if(std::scanf("%lf",&p)!=1){return 2;}}
    std::vector<Port> ports(port_count);std::vector<ducklips::CavityConnection> connections(port_count);
    for(std::size_t i=0;i<port_count;++i){auto& p=ports[i];unsigned vertical=0;
        if(std::scanf("%zu %zu %zu %u",&p.connection.from,&p.connection.to,&p.index,&vertical)!=4 || vertical>1
            || p.connection.from>=nodes || p.connection.to>=nodes || p.connection.from>=p.connection.to){return 2;}
        p.face.vertical=vertical!=0;
        if(vertical){const std::size_t column=p.index%(nx+1),row=p.index/(nx+1);if(!column || column>=nx || row>=ny){return 2;}
            p.face.coordinate=x[column];p.face.low=y[row];p.face.high=y[row+1];p.distance=.5*(x[column+1]-x[column-1]);}
        else{const std::size_t column=p.index%nx,row=p.index/nx;if(!row || row>=ny){return 2;}
            p.face.coordinate=y[row];p.face.low=x[column];p.face.high=x[column+1];p.distance=.5*(y[row+1]-y[row-1]);}
        connections[i]=p.connection;
    }
    ducklips::PlanarAirspace geometry;ducklips::PlanarFaceMotion face_motion;ducklips::CavityNetwork pressure;
    try{geometry.prepare(x,y,edges);geometry.prepare_regions(262144,524288);face_motion.prepare(edges,262144);pressure.prepare(nodes,connections);}
    catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 2;}
    const auto original_first=views(first),original_last=views(last);
    if(!geometry.build_regions(first[0],std::span<const PlanarContour>(original_first).subspan(1),depth) || geometry.regions.size()!=nodes){return 1;}
    const ducklips::AirProperties air{};const double bulk=air.density*air.sound_speed*air.sound_speed;constexpr double mobility=1e-7;
    std::vector<std::size_t> cells(nodes);
    for(std::size_t i=0;i<nodes;++i){cells[i]=geometry.regions[i].cell;volume[i]=geometry.regions[i].volume;compliance[i]=volume[i]/bulk;}
    if(!compatible(geometry,cells,ports) || !pressure.reset(volume,compliance,initial)){return 1;}
    // Pressure steps must end at actual opening/closure, not merely use an
    // average aperture spanning a change in connection. Other smooth face
    // events remain inside the geometric quadrature; time refinement controls
    // the pressure approximation between connection changes.
    std::vector<double> cuts;for(std::size_t i=0;i<=subdivisions;++i){cuts.push_back(double(i)/double(subdivisions));}
    std::vector<double> contact_cuts;
    for(const Port& p:ports){
        if(!face_motion.build(first[0],std::span<const PlanarContour>(original_first).subspan(1),last[0],std::span<const PlanarContour>(original_last).subspan(1),p.face,depth)){return 1;}
        bool first_window=true,was_open=false;
        for(const auto& w:face_motion.windows){
            const bool open=w.average_area>0;
            if(!first_window && open!=was_open){cuts.push_back(w.first);contact_cuts.push_back(w.first);}
            first_window=false;was_open=open;
        }
    }
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    std::printf("{\"contact_cuts\":[");for(std::size_t i=0;i<contact_cuts.size();++i){std::printf("%s%.17g",i?",":"",contact_cuts[i]);}
    std::printf("],\"states\":[");state(0,geometry,pressure,original_first,areas,0,false);
    Polygons a=first,b=first;const auto va=views(a),vb=views(b);
    std::vector<ducklips::CavityPassage> passages(port_count);
    double maximum_residual=0;
    for(std::size_t step=1;step<cuts.size();++step){
        const double t0=cuts[step-1],t1=cuts[step];sample(first,last,t0,a);sample(first,last,t1,b);
        for(std::size_t i=0;i<port_count;++i){
            if(!face_motion.build(a[0],std::span<const PlanarContour>(va).subspan(1),b[0],std::span<const PlanarContour>(vb).subspan(1),ports[i].face,depth)){return 1;}
            areas[i]=face_motion.average_area;const double g=mobility*areas[i]/ports[i].distance;
            passages[i]={areas[i],g>0?1/g:0,0,0,0};
        }
        if(!geometry.build_regions(b[0],std::span<const PlanarContour>(vb).subspan(1),depth) || !compatible(geometry,cells,ports)){
            std::fprintf(stderr,"declared region graph incompatible at fraction %.17g\n",t1);return 1;
        }
        for(std::size_t i=0;i<nodes;++i){volume[i]=geometry.regions[i].volume;compliance[i]=volume[i]/bulk;}
        const auto result=pressure.advance(volume,compliance,passages,inflow,duration*(t1-t0),air);
        if(!result.converged){std::fprintf(stderr,"pressure failed at fraction %.17g residual %.17g\n",t1,result.maximum_mass_residual);return 1;}
        maximum_residual=std::max(maximum_residual,result.maximum_mass_residual);
        state(t1,geometry,pressure,vb,areas,result.maximum_mass_residual,true);
    }
    std::puts("]}");std::fprintf(stderr,"steps=%zu max_mass_residual=%.17g\n",cuts.size()-1,maximum_residual);
}
