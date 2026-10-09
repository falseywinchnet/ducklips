// Native moving-geometry/inertive-pressure reference on a declared
// persistent region graph. The caller must independently establish identity;
// index/cell agreement alone is not a geometric continuation certificate.
// stdin: nx ny solids depth duration subdivisions; x edges; y edges;
// first then last contours (domain first; vertex_count, x y ...);
// node_count edge_count; initial pressures; edges: from to face vertical(0/1).
#include "ducklips/planar_airspace.hpp"
#include "ducklips/planar_face_motion.hpp"
#include "ducklips/cavity_network.hpp"
#include "ducklips/moving_tree.hpp"
#include "ducklips/moving_graph.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>
namespace {
using ducklips::Point2;using ducklips::PlanarContour;
using Polygons=std::vector<std::vector<Point2>>;
struct BoundaryLoad {std::size_t edge{},region{};double distance{},pressure{};};
struct TubeSection {double first_area{},first_length{},last_area{},last_length{};};
struct LipTube {std::vector<TubeSection> sections;std::vector<BoundaryLoad> attachments;};
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
bool compatible(const ducklips::PlanarAirspace& geometry,std::span<const std::size_t> cells,std::span<const Port> ports,std::span<const BoundaryLoad> loads,const LipTube& lip){
    if(geometry.regions.size()!=cells.size()){return false;}
    for(std::size_t i=0;i<cells.size();++i){if(geometry.regions[i].cell!=cells[i]){return false;}}
    for(const auto& link:geometry.region_links){
        std::size_t count=0;for(const auto& p:ports){
            if(link.first==p.connection.from && link.second==p.connection.to && link.face==p.index && link.vertical==p.face.vertical){++count;}
        }
        if(count!=1){return false;}
    }
    // The original tree controls permit one pair per face. The lip graph
    // can retain several pairs, but its caller must certify those apertures
    // stay open: whole-face contact events cannot detect a single pair's
    // closure while another pair on that face remains open.
    for(std::size_t i=0;i<ports.size();++i){for(std::size_t j=0;j<i;++j){
        if(lip.sections.empty() && ports[i].index==ports[j].index && ports[i].face.vertical==ports[j].face.vertical){return false;}
    }}
    for(const auto& port:geometry.boundary_ports){
        std::size_t count=0;for(const auto& load:loads){count+=load.edge==port.edge && load.region==port.region;}
        for(const auto& attachment:lip.attachments){count+=attachment.edge==port.edge && attachment.region==port.region;}
        if(count!=1){return false;}
    }
    return true;
}
void contour(PlanarContour points){
    std::putchar('[');for(std::size_t i=0;i<points.size();++i){std::printf("%s[%.17g,%.17g]",i?",":"",points[i].x,points[i].y);}std::putchar(']');
}
void geometry_json(double fraction,const ducklips::PlanarAirspace& geometry,
                   std::span<const PlanarContour> contours){
    std::printf("{\"fraction\":%.17g,\"regions\":[",fraction);
    for(std::size_t i=0;i<geometry.regions.size();++i){const auto& r=geometry.regions[i];
        std::printf("%s[%zu,%.17g,%.17g,%.17g]",i?",":"",r.cell,r.volume,r.interior.x,r.interior.y);
    }
    std::printf("],\"links\":[");
    for(std::size_t i=0;i<geometry.region_links.size();++i){const auto& r=geometry.region_links[i];
        std::printf("%s[%zu,%zu,%zu,%d,%.17g]",i?",":"",r.first,r.second,r.face,int(r.vertical),r.area);
    }
    std::putchar(']');
    std::printf(",\"boundary_ports\":[");
    for(std::size_t i=0;i<geometry.boundary_ports.size();++i){const auto& p=geometry.boundary_ports[i];
        std::printf("%s[%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g]",i?",":"",p.edge,p.region,p.area,p.centre.x,p.centre.y,p.outward_normal.x,p.outward_normal.y);}
    std::putchar(']');
    if(!contours.empty()){
        std::printf(",\"domain\":");contour(contours[0]);std::printf(",\"solids\":[");
        for(std::size_t i=1;i<contours.size();++i){if(i>1){std::putchar(',');}contour(contours[i]);}std::putchar(']');
    }
    std::putchar('}');
}
void metrics(const ducklips::PlanarAirspace& geometry,std::span<const Port> ports,std::span<const BoundaryLoad> loads,
             std::span<double> v,std::span<double> c,std::span<double> g,const LipTube& lip,double fraction){
    const ducklips::AirProperties air{};const double bulk=air.density*air.sound_speed*air.sound_speed;
    for(std::size_t i=0;i<geometry.regions.size();++i){v[i]=geometry.regions[i].volume;c[i]=v[i]/bulk;}
    std::fill(g.begin(),g.end(),0.);
    for(const auto& link:geometry.region_links){for(std::size_t i=0;i<ports.size();++i){
        if(link.first==ports[i].connection.from && link.second==ports[i].connection.to
            && link.face==ports[i].index && link.vertical==ports[i].face.vertical){g[i]+=link.area/(air.density*ports[i].distance);}
    }}
    for(const auto& port:geometry.boundary_ports){for(std::size_t i=0;i<loads.size();++i){
        if(port.edge==loads[i].edge && port.region==loads[i].region){g[ports.size()+i]+=port.area/(air.density*loads[i].distance);}
    }}
    if(lip.sections.empty()){return;}
    const std::size_t base=ports.size()+loads.size(),core=geometry.regions.size();
    double total_cap_area=0;
    for(const auto& port:geometry.boundary_ports){if(port.edge==lip.attachments[0].edge){total_cap_area+=port.area;}}
    const auto area=[&](std::size_t i){return std::lerp(lip.sections[i].first_area,lip.sections[i].last_area,fraction);};
    const auto length=[&](std::size_t i){return std::lerp(lip.sections[i].first_length,lip.sections[i].last_length,fraction);};
    for(std::size_t i=0;i<lip.sections.size();++i){
        v[core+i]=area(i)*length(i);c[core+i]=v[core+i]/bulk;
    }
    for(const auto& port:geometry.boundary_ports){for(std::size_t i=0;i<lip.attachments.size();++i){
        const auto& attachment=lip.attachments[i];
        if(port.edge==attachment.edge && port.region==attachment.region){
            // Allocate the common tube's first half-section across incoming
            // branches in proportion to cap area. The shares sum to one;
            // assigning the full tube area to each would duplicate its flow.
            const double share=port.area/total_cap_area;
            g[base+i]=1/(air.density*(attachment.distance/port.area+.5*length(0)/(area(0)*share)));
        }
    }}
    for(std::size_t i=1;i<lip.sections.size();++i){
        g[base+lip.attachments.size()+i-1]=1/(air.density*(.5*length(i-1)/area(i-1)+.5*length(i)/area(i)));
    }
    // Ideal ambient pressure at the end of the last section. Its half-length
    // is included once; no radiation resistance or end correction is implied.
    const std::size_t last=lip.sections.size()-1;
    g.back()=area(last)/(air.density*.5*length(last));
}
template<class Acoustic>
void state(const Acoustic& acoustic,std::span<const double> v,
           std::span<const double> c,std::span<const double> g){
    std::vector<double> pressure(c.size()),flow(g.size());double energy=0;
    for(std::size_t i=0;i<c.size();++i){pressure[i]=acoustic.charge[i]/c[i];energy+=.5*acoustic.charge[i]*pressure[i];}
    for(std::size_t i=0;i<g.size();++i){flow[i]=g[i]*acoustic.momentum[i];energy+=.5*g[i]*acoustic.momentum[i]*acoustic.momentum[i];}
    values("volume",v);values("compliance",c);values("conductance",g);
    values("charge",acoustic.charge);values("momentum",acoustic.momentum);
    values("pressure",pressure);values("flow",flow);values("average_flow",acoustic.average_flow);
    values("average_pressure",acoustic.average_pressure);std::printf(",\"energy\":%.17g",energy);
}

}
int main(int argc,char** argv){
    bool boundary_mode=false,lip_mode=false;
    for(int i=1;i<argc;++i){const std::string arg=argv[i];
        if(arg=="--boundary-loads" && !boundary_mode){boundary_mode=true;}
        else if(arg=="--lip-tube" && !lip_mode){lip_mode=true;}
        else{return 2;}
    }
    // Boundary loads append: count, then domain_edge region distance pressure.
    // These controls require independently established continuously open caps;
    // only the internal face-motion events are scheduled by this tool.

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
    std::vector<double> initial(nodes),volume(nodes),compliance(nodes),inflow(nodes,0);
    for(double& p:initial){if(std::scanf("%lf",&p)!=1 || !std::isfinite(p)){return 2;}}
    std::vector<Port> ports(port_count);
    for(std::size_t i=0;i<port_count;++i){auto& p=ports[i];unsigned vertical=0;
        if(std::scanf("%zu %zu %zu %u",&p.connection.from,&p.connection.to,&p.index,&vertical)!=4 || vertical>1
            || p.connection.from>=nodes || p.connection.to>=nodes || p.connection.from>=p.connection.to){return 2;}
        p.face.vertical=vertical!=0;
        if(vertical){const std::size_t column=p.index%(nx+1),row=p.index/(nx+1);if(!column || column>=nx || row>=ny){return 2;}
            p.face.coordinate=x[column];p.face.low=y[row];p.face.high=y[row+1];p.distance=.5*(x[column+1]-x[column-1]);}
        else{const std::size_t column=p.index%nx,row=p.index/nx;if(!row || row>=ny){return 2;}
            p.face.coordinate=y[row];p.face.low=x[column];p.face.high=x[column+1];p.distance=.5*(y[row+1]-y[row-1]);}
    }
    std::vector<BoundaryLoad> loads;std::vector<std::size_t> selected;
    if(boundary_mode){
        std::size_t count=0;if(std::scanf("%zu",&count)!=1 || count>262144){return 2;}loads.resize(count);
        for(std::size_t i=0;i<count;++i){auto& b=loads[i];
            if(std::scanf("%zu %zu %lf %lf",&b.edge,&b.region,&b.distance,&b.pressure)!=4 || b.edge>=first[0].size()
                || b.region>=nodes || !std::isfinite(b.distance) || b.distance<=0 || !std::isfinite(b.pressure)){return 2;}
            for(std::size_t j=0;j<i;++j){if(loads[j].edge==b.edge && loads[j].region==b.region){return 2;}}
            if(std::find(selected.begin(),selected.end(),b.edge)==selected.end()){selected.push_back(b.edge);}
        }
    }
    LipTube lip;
    // After optional pressure loads: section_count attachment_count;
    // rows area0 length0 area1 length1; rows domain_edge region core_distance.
    // End values are final assembled articulator sections, linearly sampled
    // within this declared geometry interval. Zero-storage closure is rejected.
    if(lip_mode){
        std::size_t sections=0,attachments=0;
        if(std::scanf("%zu %zu",&sections,&attachments)!=2 || !sections || sections>256 || !attachments || attachments>nodes){return 2;}
        lip.sections.resize(sections);lip.attachments.resize(attachments);
        for(auto& s:lip.sections){
            if(std::scanf("%lf %lf %lf %lf",&s.first_area,&s.first_length,&s.last_area,&s.last_length)!=4){return 2;}
            for(double value:{s.first_area,s.first_length,s.last_area,s.last_length}){
                if(!std::isfinite(value) || value<=0){std::fprintf(stderr,"lip reference requires positive section area and length; closed-lip storage transport is not implemented\n");return 2;}
            }
        }
        for(std::size_t i=0;i<attachments;++i){auto& b=lip.attachments[i];
            if(std::scanf("%zu %zu %lf",&b.edge,&b.region,&b.distance)!=3 || b.edge>=first[0].size() || b.region>=nodes
                || !std::isfinite(b.distance) || b.distance<=0){return 2;}
            if(i && b.edge!=lip.attachments[0].edge){return 2;}
            for(std::size_t j=0;j<i;++j){if(lip.attachments[j].region==b.region){return 2;}}
            for(const auto& load:loads){if(load.edge==b.edge){return 2;}}
            if(std::find(selected.begin(),selected.end(),b.edge)==selected.end()){selected.push_back(b.edge);}
        }
    }
    const std::size_t total_nodes=nodes+lip.sections.size();
    if(lip_mode && total_nodes>256){return 2;}
    volume.resize(total_nodes);compliance.resize(total_nodes);inflow.resize(total_nodes);
    const std::size_t channels=ports.size()+loads.size()+lip.attachments.size()+lip.sections.size();std::vector<double> series(channels,0.);
    for(std::size_t i=0;i<loads.size();++i){series[ports.size()+i]=-loads[i].pressure;}
    ducklips::PlanarAirspace geometry;ducklips::PlanarFaceMotion face_motion;
    ducklips::MovingTreeReference acoustic;ducklips::MovingGraphReference graph;
    try{
        geometry.prepare(x,y,edges);geometry.prepare_regions(262144,524288);face_motion.prepare(edges,262144);
        std::vector<ducklips::TreeEdge> branches;
        for(const auto& p:ports){branches.push_back({p.connection.from,p.connection.to});}
        for(const auto& b:loads){branches.push_back({b.region,ducklips::ground_node});}
        for(const auto& b:lip.attachments){branches.push_back({b.region,nodes});}
        for(std::size_t i=1;i<lip.sections.size();++i){branches.push_back({nodes+i-1,nodes+i});}
        if(lip_mode){branches.push_back({total_nodes-1,ducklips::ground_node});graph.prepare(total_nodes,branches);}
        else{acoustic.prepare(nodes,branches);}
    }
    catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 2;}
    const auto original_first=views(first),original_last=views(last);
    if(!geometry.build_regions(first[0],std::span<const PlanarContour>(original_first).subspan(1),depth,selected) || geometry.regions.size()!=nodes){return 1;}
    std::vector<std::size_t> cells(nodes);
    std::vector<double> gend(channels),old_volume(total_nodes),c1(total_nodes),c2(total_nodes),g1(channels),g2(channels),scratch(total_nodes);
    for(std::size_t i=0;i<nodes;++i){cells[i]=geometry.regions[i].cell;}
    if(!compatible(geometry,cells,ports,loads,lip)){return 1;}
    metrics(geometry,ports,loads,volume,compliance,gend,lip,0);
    auto& charge=lip_mode?graph.charge:acoustic.charge;
    auto& momentum=lip_mode?graph.momentum:acoustic.momentum;
    for(std::size_t i=0;i<nodes;++i){charge[i]=compliance[i]*initial[i];}
    const auto snapshot=[&](){if(lip_mode){state(graph,volume,compliance,gend);}else{state(acoustic,volume,compliance,gend);}};
    // Contact times split the integration. No dissipative mobility is assigned
    // to artificial grid faces. Inverse inertance is A/(rho*centre_distance).
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
    std::printf("],\"initial\":{\"geometry\":");geometry_json(0,geometry,original_first);
    snapshot();std::printf("},\"steps\":[");
    Polygons a=first;const auto va=views(a);
    constexpr double offset=.28867513459481288225;
    for(std::size_t step=1;step<cuts.size();++step){
        const double t0=cuts[step-1],t1=cuts[step],width=t1-t0,dt=duration*width;
        const double first_stage=t0+(.5-offset)*width,second_stage=t0+(.5+offset)*width;
        old_volume=volume;
        std::printf("%s{\"dt\":%.17g,\"stage1\":",step>1?",":"",dt);
        sample(first,last,first_stage,a);
        if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip)){return 1;}
        metrics(geometry,ports,loads,scratch,c1,g1,lip,first_stage);geometry_json(first_stage,geometry,{});
        std::printf(",\"stage2\":");sample(first,last,second_stage,a);
        if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip)){return 1;}
        metrics(geometry,ports,loads,scratch,c2,g2,lip,second_stage);geometry_json(second_stage,geometry,{});
        std::printf(",\"geometry\":");sample(first,last,t1,a);
        if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip)){return 1;}
        metrics(geometry,ports,loads,volume,compliance,gend,lip,t1);geometry_json(t1,geometry,va);
        for(std::size_t i=0;i<total_nodes;++i){inflow[i]=-(volume[i]-old_volume[i])/dt;}
        // A closed passage has no kinetic state to accumulate pressure history.
        // The zero-metric momentum is discarded, as in CanonicalSpeech.
        for(std::size_t i=0;i<channels;++i){if(g1[i]==0 && g2[i]==0){momentum[i]=0;}}
        if(lip_mode){if(!graph.advance(c1,c2,g1,g2,dt,inflow,series)){std::fprintf(stderr,"shared-lip solve failed\n");return 1;}}
        else{acoustic.advance(c1,c2,g1,g2,dt,inflow,series);}
        for(std::size_t i=0;i<channels;++i){if(gend[i]==0){momentum[i]=0;}}
        values("wall_flow",inflow);snapshot();std::putchar('}');
    }
    std::puts("]}");std::fprintf(stderr,"steps=%zu; inviscid canonical moving-volume reference; no speech\n",cuts.size()-1);
}
