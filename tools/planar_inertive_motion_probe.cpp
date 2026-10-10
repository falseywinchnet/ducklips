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
#include "ducklips/lip_profile.hpp"
#include "ducklips/lip_viscous.hpp"
#include "ducklips/graph_orifice.hpp"
#include "ducklips/graph_phonation.hpp"
#include "ducklips/graph_continuous_glottis.hpp"
#include "ducklips/joint_accuracy.hpp"
#include "ducklips/body_config.hpp"
#include <algorithm>
#include <array>
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
struct LipTube {
    std::vector<TubeSection> sections;std::vector<BoundaryLoad> attachments;
    bool integrated{};ducklips::LipProfile first_profile{},last_profile{};
    ducklips::LipProfile profile(double t)const{
        const auto& a=first_profile;const auto& b=last_profile;
        return {std::lerp(a.length,b.length,t),std::lerp(a.base_area,b.base_area,t),std::lerp(a.activation,b.activation,t),
                std::lerp(a.place,b.place,t),std::lerp(a.core,b.core,t),std::lerp(a.width,b.width,t),std::lerp(a.target_area,b.target_area,t)};
    }
};
struct Epilarynx {std::vector<TubeSection> sections;std::vector<BoundaryLoad> attachments;};
struct Port {ducklips::CavityConnection connection{};ducklips::PlanarFace face{};std::size_t index{};double distance{};};
void sample(const Polygons& first,const Polygons& last,double t,Polygons& result){
    for(std::size_t c=0;c<first.size();++c){for(std::size_t v=0;v<first[c].size();++v){
        result[c][v]={std::lerp(first[c][v].x,last[c][v].x,t),std::lerp(first[c][v].y,last[c][v].y,t)};
    }}
}
std::vector<PlanarContour> views(const Polygons& p){std::vector<PlanarContour> r;for(const auto& c:p){r.emplace_back(c);}return r;}
void fold_state(const ducklips::TwoMassGlottis& folds){
    const auto s=folds.state();
    std::printf(",\"fold_state\":[%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g]",
        s.lower_position[0],s.lower_position[1],s.lower_velocity[0],s.lower_velocity[1],
        s.upper_position[0],s.upper_position[1],s.upper_velocity[0],s.upper_velocity[1],s.rest_area);
}
void values(const char* name,std::span<const double> a){
    std::printf(",\"%s\":[",name);for(std::size_t i=0;i<a.size();++i){std::printf("%s%.17g",i?",":"",a[i]);}std::putchar(']');
}
bool compatible(const ducklips::PlanarAirspace& geometry,std::span<const std::size_t> cells,std::span<const Port> ports,std::span<const BoundaryLoad> loads,const LipTube& lip,const Epilarynx& epi){
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
        for(const auto& attachment:epi.attachments){count+=attachment.edge==port.edge && attachment.region==port.region;}
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
bool lip_metrics(const ducklips::PlanarAirspace& geometry,std::span<const Port> ports,std::span<const BoundaryLoad> loads,
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
    if(lip.sections.empty()){return true;}
    const std::size_t base=ports.size()+loads.size(),core=geometry.regions.size();
    double total_cap_area=0;
    for(const auto& port:geometry.boundary_ports){if(port.edge==lip.attachments[0].edge){total_cap_area+=port.area;}}
    if(lip.integrated){
        const auto profile=lip.profile(fraction);const std::size_t count=lip.sections.size();
        if(profile.minimum_area()==0 && profile.place-profile.core>0 && profile.place+profile.core<profile.length){
            std::fprintf(stderr,"interior lip contact requires separate storage regions\n");return false;
        }
        std::array<double,256> centre{};
        for(std::size_t i=0;i<count;++i){
            const double a=profile.length*double(count-i-1)/double(count),b=profile.length*double(count-i)/double(count);
            if(!profile.measure(a,b,v[core+i],centre[i])){std::fprintf(stderr,"lip volume integral failed\n");return false;}
            if(v[core+i]<=0){std::fprintf(stderr,"empty lip storage requires state transport\n");return false;}
            c[core+i]=v[core+i]/bulk;
        }
        double entry=0;
        if(!profile.conductance(centre[0],profile.length,air.density,entry)){return false;}
        for(const auto& port:geometry.boundary_ports){for(std::size_t i=0;i<lip.attachments.size();++i){
            const auto& attachment=lip.attachments[i];
            if(port.edge==attachment.edge && port.region==attachment.region && entry>0){
                const double share=port.area/total_cap_area;
                g[base+i]=1/(air.density*attachment.distance/port.area+1/(entry*share));
            }
        }}
        for(std::size_t i=1;i<count;++i){
            if(!profile.conductance(centre[i],centre[i-1],air.density,g[base+lip.attachments.size()+i-1])){return false;}
        }
        return profile.conductance(0,centre[count-1],air.density,g[base+lip.attachments.size()+count-1]);
    }
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
    g[base+lip.attachments.size()+last]=area(last)/(air.density*.5*length(last));
    return true;
}
bool metrics(const ducklips::PlanarAirspace& geometry,std::span<const Port> ports,std::span<const BoundaryLoad> loads,
             std::span<double> v,std::span<double> c,std::span<double> g,const LipTube& lip,const Epilarynx& epi,double fraction){
    if(!lip_metrics(geometry,ports,loads,v,c,g,lip,fraction)){return false;}
    if(epi.sections.empty()){return true;}
    const ducklips::AirProperties air{};const double bulk=air.density*air.sound_speed*air.sound_speed;
    const std::size_t node=geometry.regions.size()+lip.sections.size();
    const std::size_t branch=ports.size()+loads.size()+lip.attachments.size()+lip.sections.size();
    const auto area=[&](std::size_t i){return std::lerp(epi.sections[i].first_area,epi.sections[i].last_area,fraction);};
    const auto length=[&](std::size_t i){return std::lerp(epi.sections[i].first_length,epi.sections[i].last_length,fraction);};
    for(std::size_t i=0;i<epi.sections.size();++i){v[node+i]=area(i)*length(i);c[node+i]=v[node+i]/bulk;}
    double total=0;for(const auto& port:geometry.boundary_ports){if(port.edge==epi.attachments[0].edge){total+=port.area;}}
    for(const auto& port:geometry.boundary_ports){for(std::size_t i=0;i<epi.attachments.size();++i){
        const auto& a=epi.attachments[i];if(port.edge==a.edge && port.region==a.region){
            const double share=port.area/total;
            g[branch+i]=1/(air.density*(a.distance/port.area+.5*length(0)/(area(0)*share)));
        }
    }}
    for(std::size_t i=1;i<epi.sections.size();++i){
        g[branch+epi.attachments.size()+i-1]=1/(air.density*(.5*length(i-1)/area(i-1)+.5*length(i)/area(i)));
    }
    // The remaining glottis-to-centre half-section is handled by the
    // source port's own canonical momentum, not another storage cell.
    return true;
}
bool continuous_state(const ducklips::GraphContinuousGlottis& model){
    ducklips::ContinuousGlottalState s;if(!model.source.state(s)){return false;}
    std::printf(",\"continuous_source\":{\"energy\":%.17g,\"rest_area\":%.17g,\"outlet_inertance\":%.17g",s.energy,s.folds.rest_area,s.outlet_inertance);
    values("mass",s.mass);values("effort",s.effort);values("flow",s.flow);values("area",s.area);
    values("lower_position",s.folds.lower_position);values("lower_velocity",s.folds.lower_velocity);
    values("upper_position",s.folds.upper_position);values("upper_velocity",s.folds.upper_velocity);
    values("canonical_momentum",s.momentum);std::putchar('}');return true;
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
    bool boundary_mode=false,lip_mode=false,profile_mode=false,initial_state=false,compact_geometry=false,epi_mode=false,orifice_mode=false,fold_mode=false,continuous_mode=false,radiation_mode=false;
    std::string fold_body;
    double thermal_width=0,joint_tolerance=0,joint_maximum_step=1e-5;
    bool joint_fifth=false,joint_relaxation=false,outlet_viscous=false;double viscous_width=0,lip_viscous_width=0;std::size_t lip_viscous_order=8;
    for(int i=1;i<argc;++i){const std::string arg=argv[i];
        if(arg=="--boundary-loads" && !boundary_mode){boundary_mode=true;}
        else if(arg=="--lip-tube" && !lip_mode){lip_mode=true;}
        else if(arg=="--lip-profile" && !lip_mode){lip_mode=true;profile_mode=true;}
        else if(arg=="--epilarynx" && !epi_mode){epi_mode=true;}
        else if(arg=="--glottal-orifice" && !orifice_mode){orifice_mode=true;}
        else if(arg=="--continuous-folds" && !orifice_mode && i+1<argc){orifice_mode=true;continuous_mode=true;fold_body=argv[++i];}
        else if(arg=="--vocal-folds" && !orifice_mode && i+1<argc){orifice_mode=true;fold_mode=true;fold_body=argv[++i];}
        else if(arg=="--radiation" && !radiation_mode){radiation_mode=true;}
        else if(arg=="--thermal-slit-width" && thermal_width==0 && i+1<argc){
            try{std::size_t used=0;const std::string value=argv[++i];thermal_width=std::stod(value,&used);if(used!=value.size()){return 2;}}
            catch(...){return 2;}
            if(!std::isfinite(thermal_width) || thermal_width<=0){return 2;}
        }
        else if(arg=="--viscous-slit-width" && viscous_width==0 && i+1<argc){
            try{std::size_t used=0;const std::string value=argv[++i];viscous_width=std::stod(value,&used);if(used!=value.size()){return 2;}}
            catch(...){return 2;}
            if(!std::isfinite(viscous_width) || viscous_width<=0){return 2;}
        }
        else if(arg=="--lip-viscous-width" && lip_viscous_width==0 && i+1<argc){
            try{std::size_t used=0;const std::string value=argv[++i];lip_viscous_width=std::stod(value,&used);if(used!=value.size()){return 2;}}
            catch(...){return 2;}
            if(!std::isfinite(lip_viscous_width) || lip_viscous_width<=0){return 2;}
        }
        else if((arg=="--joint-tolerance" || arg=="--joint-max-step") && i+1<argc){
            double value;try{std::size_t used=0;const std::string token=argv[++i];value=std::stod(token,&used);if(used!=token.size()){return 2;}}catch(...){return 2;}
            if(!std::isfinite(value) || value<=0){return 2;}if(arg=="--joint-tolerance"){joint_tolerance=value;}else{joint_maximum_step=value;}
        }
        else if(arg=="--joint-fifth" && !joint_fifth){joint_fifth=true;}
        else if(arg=="--joint-relaxation" && !joint_relaxation){joint_relaxation=true;}
        else if(arg=="--outlet-viscous" && !outlet_viscous){outlet_viscous=true;}
        else if(arg=="--lip-viscous-order" && i+1<argc){
            try{std::size_t used=0;const std::string value=argv[++i];lip_viscous_order=std::stoul(value,&used);if(used!=value.size()){return 2;}}
            catch(...){return 2;}
            if(lip_viscous_order!=4 && lip_viscous_order!=8 && lip_viscous_order!=16 && lip_viscous_order!=32){return 2;}
        }
        else if(arg=="--initial-state" && !initial_state){initial_state=true;}
        else if(arg=="--compact-geometry" && !compact_geometry){compact_geometry=true;}
        else{return 2;}
    }
    const bool viscous_mode=viscous_width>0 || lip_viscous_width>0;
    if((joint_tolerance>0 && !joint_relaxation) || (joint_fifth && !joint_relaxation) || (joint_relaxation && !continuous_mode) || (outlet_viscous && lip_viscous_width<=0) || (lip_viscous_width>0 && !profile_mode) || ((radiation_mode || thermal_width>0 || viscous_mode) && !continuous_mode) || (epi_mode && !lip_mode) || (orifice_mode && !epi_mode)){return 2;}
    // Boundary loads append: count, then domain_edge region distance pressure.
    // These controls require independently established continuously open caps;
    // only the internal face-motion events are scheduled by this tool.

    std::size_t nx=0,ny=0,solid_count=0,subdivisions=0;double depth=0,duration=0;
    if(std::scanf("%zu %zu %zu %lf %lf %zu",&nx,&ny,&solid_count,&depth,&duration,&subdivisions)!=6
        || !nx || !ny || nx>256 || ny>256 || solid_count>32 || !subdivisions || subdivisions>(profile_mode?8192u:2048u)
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
        lip.sections.resize(sections);lip.attachments.resize(attachments);lip.integrated=profile_mode;
        // Profile mode replaces the section rows with two endpoint profiles:
        // length base_area activation place core width target_area (SI),
        // followed by the initial pressure of each lip storage node.
        if(profile_mode){
            for(auto* p:{&lip.first_profile,&lip.last_profile}){
                if(std::scanf("%lf %lf %lf %lf %lf %lf %lf",&p->length,&p->base_area,&p->activation,&p->place,&p->core,&p->width,&p->target_area)!=7 || !p->valid()){return 2;}
            }
            initial.resize(nodes+sections);
            for(std::size_t i=nodes;i<initial.size();++i){if(std::scanf("%lf",&initial[i])!=1 || !std::isfinite(initial[i])){return 2;}}
        }else for(auto& s:lip.sections){
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
    Epilarynx epi;
    // Epilarynx rows are ordered from mouth seam toward the glottis:
    // count attachments; area0 length0 area1 length1 rows; initial pressures;
    // domain_edge region core_distance attachment rows. Positive storage only.
    if(epi_mode){
        std::size_t count=0,attachments=0;
        if(std::scanf("%zu %zu",&count,&attachments)!=2 || !count || count>256 || !attachments || attachments>nodes){return 2;}
        epi.sections.resize(count);epi.attachments.resize(attachments);
        for(auto& e:epi.sections){
            if(std::scanf("%lf %lf %lf %lf",&e.first_area,&e.first_length,&e.last_area,&e.last_length)!=4){return 2;}
            for(double value:{e.first_area,e.first_length,e.last_area,e.last_length}){if(!std::isfinite(value) || value<=0){return 2;}}
        }
        const std::size_t first_epi=nodes+lip.sections.size();initial.resize(first_epi+count,0.);
        for(std::size_t i=first_epi;i<initial.size();++i){if(std::scanf("%lf",&initial[i])!=1 || !std::isfinite(initial[i])){return 2;}}
        for(std::size_t i=0;i<attachments;++i){auto& b=epi.attachments[i];
            if(std::scanf("%zu %zu %lf",&b.edge,&b.region,&b.distance)!=3 || b.edge>=first[0].size() || b.region>=nodes
                || !std::isfinite(b.distance) || b.distance<=0){return 2;}
            if(i && b.edge!=epi.attachments[0].edge){return 2;}
            for(std::size_t j=0;j<i;++j){if(epi.attachments[j].region==b.region){return 2;}}
            for(const auto& load:loads){if(load.edge==b.edge){return 2;}}
            for(const auto& attachment:lip.attachments){if(attachment.edge==b.edge){return 2;}}
            if(std::find(selected.begin(),selected.end(),b.edge)==selected.end()){selected.push_back(b.edge);}
        }
    }
    double glottal_area=0,reservoir_pressure=0,tension=1,surface_pressure=0;
    if(orifice_mode && (std::scanf("%lf %lf",&glottal_area,&reservoir_pressure)!=2 || !std::isfinite(glottal_area)
        || (!fold_mode && glottal_area<0) || !std::isfinite(reservoir_pressure))){return 2;}
    if(fold_mode && (std::scanf("%lf %lf",&tension,&surface_pressure)!=2 || !std::isfinite(tension) || tension<=0 || !std::isfinite(surface_pressure))){return 2;}
    if(continuous_mode && (std::scanf("%lf",&tension)!=1 || !std::isfinite(tension) || tension<=0)){return 2;}
    const std::size_t total_nodes=nodes+lip.sections.size()+epi.sections.size();
    initial.resize(total_nodes,0.);
    if(lip_mode && total_nodes>256){return 2;}
    volume.resize(total_nodes);compliance.resize(total_nodes);inflow.resize(total_nodes);
    const std::size_t channels=ports.size()+loads.size()+lip.attachments.size()+lip.sections.size()
        +(epi_mode?epi.attachments.size()+epi.sections.size()-1:0);std::vector<double> series(channels,0.);
    for(std::size_t i=0;i<loads.size();++i){series[ports.size()+i]=-loads[i].pressure;}
    ducklips::PlanarAirspace geometry;ducklips::PlanarFaceMotion face_motion;
    ducklips::MovingTreeReference acoustic;ducklips::MovingGraphReference passive;ducklips::GraphOrifice orifice;
    ducklips::GraphContinuousGlottis continuous;ducklips::VoiceProfile body;
    ducklips::GraphPhonation phonation;auto& source=fold_mode?phonation.source:orifice;
    auto& graph=continuous_mode?continuous.network:(orifice_mode?source.state:passive);
    try{
        geometry.prepare(x,y,edges);geometry.prepare_regions(262144,524288);face_motion.prepare(edges,262144);
        std::vector<ducklips::TreeEdge> branches;
        for(const auto& p:ports){branches.push_back({p.connection.from,p.connection.to});}
        for(const auto& b:loads){branches.push_back({b.region,ducklips::ground_node});}
        for(const auto& b:lip.attachments){branches.push_back({b.region,nodes});}
        for(std::size_t i=1;i<lip.sections.size();++i){branches.push_back({nodes+i-1,nodes+i});}
        if(lip_mode){
            branches.push_back({nodes+lip.sections.size()-1,ducklips::ground_node});
            const std::size_t base=nodes+lip.sections.size();
            for(const auto& b:epi.attachments){branches.push_back({b.region,base});}
            for(std::size_t i=1;i<epi.sections.size();++i){branches.push_back({base+i-1,base+i});}
            if(fold_mode || continuous_mode){
                std::string error;
                if(!ducklips::load_body_config(fold_body,body,error)){std::fprintf(stderr,"%s\n",error.c_str());return 2;}
                if(continuous_mode){
                    ducklips::FoldState seed;seed.rest_area=glottal_area;seed.lower_position.fill(1e-6);
                    const auto& e=epi.sections.back();const double column=ducklips::AirProperties{}.density*.5*e.first_length/e.first_area;
                    continuous.prepare(total_nodes,branches,total_nodes-1,body.folds,{},tension,seed,{initial.back(),initial.back()},{},1.,ducklips::ContourPressureProfile::linear_span,column);
                }else{phonation.prepare(total_nodes,branches,total_nodes-1,body.folds);}
            }else if(orifice_mode){source.prepare(total_nodes,branches,total_nodes-1);}else{graph.prepare(total_nodes,branches);}
        }
        else{acoustic.prepare(nodes,branches);}
    }
    catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 2;}
    const auto radiation_area=[&](double fraction){return lip.integrated?lip.profile(fraction).area(0.):std::lerp(lip.sections.back().first_area,lip.sections.back().last_area,fraction);};
    const std::size_t radiation_branch=ports.size()+loads.size()+lip.attachments.size()+lip.sections.size()-1;
    if(radiation_mode){
        const ducklips::NetworkLosses parameters;
        if(!continuous.radiation.prepare({},radiation_area(0.),parameters.end_correction,parameters.radiation_factor)){return 2;}
    }
    const auto original_first=views(first),original_last=views(last);
    if(!geometry.build_regions(first[0],std::span<const PlanarContour>(original_first).subspan(1),depth,selected) || geometry.regions.size()!=nodes){return 1;}
    std::vector<std::size_t> cells(nodes);
    std::vector<double> gend(channels),old_volume(total_nodes),old_compliance(total_nodes),old_conductance(channels),c1(total_nodes),c2(total_nodes),g1(channels),g2(channels),scratch(total_nodes);
    for(std::size_t i=0;i<nodes;++i){cells[i]=geometry.regions[i].cell;}
    if(!compatible(geometry,cells,ports,loads,lip,epi)){return 1;}
    if(!metrics(geometry,ports,loads,volume,compliance,gend,lip,epi,0)){return 1;}
    auto& charge=lip_mode?graph.charge:acoustic.charge;
    auto& momentum=lip_mode?graph.momentum:acoustic.momentum;
    for(std::size_t i=0;i<total_nodes;++i){charge[i]=compliance[i]*initial[i];}
    const std::size_t thermal_base=nodes+lip.sections.size();
    const std::size_t viscous_base=ports.size()+loads.size()+lip.attachments.size()+lip.sections.size()+epi.attachments.size();
    const std::size_t lip_viscous_offset=viscous_width>0?2*(epi.sections.size()-1):0;
    const std::size_t lip_viscous_branch=ports.size()+loads.size()+lip.attachments.size();
    std::vector<ducklips::ViscousGeometry> viscous_geometry(lip_viscous_offset+(lip_viscous_width>0?(lip.sections.size()-1)*lip_viscous_order+(outlet_viscous?lip_viscous_order+1:0):0));
    std::vector<std::size_t> viscous_branches(viscous_geometry.size());std::vector<double> viscous_positions(viscous_geometry.size());
    const auto set_viscous_geometry=[&](double fraction,std::span<const double> conductance){
        if(viscous_width>0){for(std::size_t i=1;i<epi.sections.size();++i){for(std::size_t half=0;half<2;++half){const auto& s=epi.sections[i-1+half];
            const double area=std::lerp(s.first_area,s.last_area,fraction),length=std::lerp(s.first_length,s.last_length,fraction);
            const auto index=2*(i-1)+half;const auto branch=viscous_base+i-1;viscous_branches[index]=branch;
            viscous_geometry[index]={ducklips::AirProperties{}.density*.5*length/area,conductance[branch],area/(2*viscous_width)};
        }}}
        if(lip_viscous_width>0){
            const auto profile=lip.profile(fraction);std::array<double,256> centre{};
            for(std::size_t i=0;i<lip.sections.size();++i){double unused;
                if(!profile.measure(profile.length*double(lip.sections.size()-i-1)/double(lip.sections.size()),profile.length*double(lip.sections.size()-i)/double(lip.sections.size()),unused,centre[i])){return false;}
            }
            for(std::size_t i=1;i<lip.sections.size();++i){std::array<ducklips::LipViscousSegment,32> segments{};ducklips::LipViscousMeasure measure;
                if(!ducklips::lip_viscous_segments(profile,centre[i],centre[i-1],lip_viscous_width,ducklips::AirProperties{}.density,std::span(segments).first(lip_viscous_order),measure)){return false;}
                const auto branch=lip_viscous_branch+i-1;
                for(std::size_t j=0;j<lip_viscous_order;++j){const auto index=lip_viscous_offset+(i-1)*lip_viscous_order+j;const auto& segment=segments[j];
                    viscous_branches[index]=branch;viscous_positions[index]=segment.position;
                    viscous_geometry[index]={segment.inertance,conductance[branch],segment.half_gap};
                }
            }
            if(outlet_viscous){
                std::array<ducklips::LipViscousSegment,33> segments{};ducklips::LipViscousMeasure measure;
                if(!ducklips::lip_outlet_viscous_segments(profile,centre[lip.sections.size()-1],lip_viscous_width,ducklips::AirProperties{}.density,std::span(segments).first(lip_viscous_order+1),measure)){return false;}
                const auto branch=lip_viscous_branch+lip.sections.size()-1;
                const auto offset=lip_viscous_offset+(lip.sections.size()-1)*lip_viscous_order;
                for(std::size_t j=0;j<=lip_viscous_order;++j){const auto index=offset+j;const auto& segment=segments[j];
                    viscous_branches[index]=branch;viscous_positions[index]=segment.position;
                    viscous_geometry[index]={segment.inertance,conductance[branch],segment.half_gap};
                }
            }
        }
        return true;
    };
    if(viscous_mode){
        if(!set_viscous_geometry(0,gend)){return 2;}std::vector<ducklips::ViscousSite> sites;
        for(std::size_t i=0;i<viscous_geometry.size();++i){sites.push_back({viscous_branches[i],viscous_geometry[i]});}
        try{continuous.prepare_viscous({},sites);}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 2;}
    }
    std::vector<ducklips::ThermalGeometry> thermal_geometry;
    if(thermal_width>0){
        std::vector<ducklips::ThermalSite> sites;thermal_geometry.resize(epi.sections.size());
        for(std::size_t i=0;i<epi.sections.size();++i){sites.push_back({thermal_base+i,compliance[thermal_base+i],epi.sections[i].first_area/(2*thermal_width)});}
        try{continuous.prepare_thermal({},sites);}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 2;}
        for(std::size_t i=0;i<sites.size();++i){if(!continuous.thermal[i].equilibrate(initial[thermal_base+i])){return 2;}}
    }
    std::vector<ducklips::ThermalGeometry> midpoint_thermal(thermal_geometry.size());
    std::vector<ducklips::ViscousGeometry> midpoint_viscous(viscous_geometry.size());
    if(joint_relaxation){try{continuous.prepare_joint(radiation_mode?radiation_branch:ducklips::ground_node,joint_fifth?ducklips::RelaxationOrder::fifth:ducklips::RelaxationOrder::third);}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 2;}}
    // Optional continuation appends canonical charge then branch momentum.
    // The caller must certify identical endpoint geometry and node/branch
    // ownership across calls; this is not a general remeshing operation.
    if(initial_state){
        for(double& q:charge){if(std::scanf("%lf",&q)!=1 || !std::isfinite(q)){return 2;}}
        for(double& m:momentum){if(std::scanf("%lf",&m)!=1 || !std::isfinite(m)){return 2;}}
        if(orifice_mode && !continuous_mode && (std::scanf("%lf",&source.momentum)!=1 || !std::isfinite(source.momentum))){return 2;}
        if(continuous_mode){
            std::array<double,2> effort;std::array<double,3> flow;ducklips::FoldState restored;restored.rest_area=glottal_area;
            for(double& v:effort){if(std::scanf("%lf",&v)!=1){return 2;}}
            for(double& v:flow){if(std::scanf("%lf",&v)!=1){return 2;}}
            for(auto* field:{&restored.lower_position,&restored.lower_velocity,&restored.upper_position,&restored.upper_velocity}){
                for(double& v:*field){if(std::scanf("%lf",&v)!=1){return 2;}}
            }
            const auto& e=epi.sections.back();const double column=ducklips::AirProperties{}.density*.5*e.first_length/e.first_area;
            if(!continuous.source.prepare(body.folds,{},tension,restored,effort,flow,1.,ducklips::ContourPressureProfile::linear_span,column)){return 2;}
        }
        if(radiation_mode){
            double value;if(std::scanf("%lf",&value)!=1 || !std::isfinite(value) || (continuous.radiation.area()==0 && value!=0)){return 2;}
            continuous.radiation.momentum=value;
        }
        if(thermal_width>0){for(auto& site:continuous.thermal){for(double& z:site.charge){if(std::scanf("%lf",&z)!=1 || !std::isfinite(z)){return 2;}}}}
        if(viscous_mode){for(auto& site:continuous.viscous){for(double& z:site.momentum){if(std::scanf("%lf",&z)!=1 || !std::isfinite(z)){return 2;}}}}
        if(fold_mode){
            ducklips::FoldState s;
            if(std::scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf",&s.lower_position[0],&s.lower_position[1],&s.lower_velocity[0],&s.lower_velocity[1],
                &s.upper_position[0],&s.upper_position[1],&s.upper_velocity[0],&s.upper_velocity[1],&s.rest_area)!=9 || !phonation.folds.restore(s)){return 2;}
        }
        for(std::size_t i=0;i<channels;++i){if(gend[i]==0){momentum[i]=0;}}
    }
    const auto snapshot=[&](){if(radiation_mode){
        std::printf(",\"radiation_state\":{\"area\":%.17g,\"momentum\":%.17g,\"inverse_inertance\":%.17g,\"energy\":%.17g}",continuous.radiation.area(),continuous.radiation.momentum,continuous.radiation.inverse_inertance(),continuous.radiation.energy());
    }if(thermal_width>0){
        std::printf(",\"thermal_sites\":[");for(std::size_t i=0;i<continuous.thermal.size();++i){const auto& site=continuous.thermal[i];
            std::printf("%s{\"node\":%zu,\"compliance\":%.17g,\"half_gap\":%.17g,\"stored_air\":%.17g,\"energy\":%.17g",i?",":"",thermal_base+i,site.compliance(),site.half_gap(),site.stored_air(),site.energy());
            values("charge",site.charge);std::putchar('}');}std::putchar(']');
    }if(viscous_mode){
        std::printf(",\"viscous_sites\":[");for(std::size_t i=0;i<continuous.viscous.size();++i){const auto& site=continuous.viscous[i];
            std::printf("%s{\"branch\":%zu,\"segment_inertance\":%.17g,\"branch_conductance\":%.17g,\"half_gap\":%.17g,\"energy\":%.17g",i?",":"",viscous_branches[i],site.segment_inertance(),site.branch_conductance(),site.half_gap(),site.energy());
            if(lip_viscous_width>0 && i>=lip_viscous_offset){std::printf(",\"position\":%.17g",viscous_positions[i]);}
            values("momentum",site.momentum);std::putchar('}');}std::putchar(']');
    }if(lip_mode){state(graph,volume,compliance,gend);}else{state(acoustic,volume,compliance,gend);}};
    // Contact times split the integration. No dissipative mobility is assigned
    // to artificial grid faces. Inverse inertance is A/(rho*centre_distance).
    std::vector<double> cuts;for(std::size_t i=0;i<=subdivisions;++i){cuts.push_back(double(i)/double(subdivisions));}
    std::vector<double> contact_cuts;
    std::vector<double> lip_contact_cuts;
    if(profile_mode){
        const auto& p=lip.first_profile;const auto& q=lip.last_profile;
        if(p.target_area==0 && q.target_area==0 && (p.activation<1)!=(q.activation<1)){
            const double contact=(1-p.activation)/(q.activation-p.activation);
            if(contact>0 && contact<1){cuts.push_back(contact);lip_contact_cuts.push_back(contact);}
        }
    }
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
    std::putchar(']');
    if(profile_mode){values("lip_contact_cuts",lip_contact_cuts);}
    std::printf(",\"initial\":{\"geometry\":");geometry_json(0,geometry,original_first);
    snapshot();if(continuous_mode && !continuous_state(continuous)){return 1;}if(orifice_mode && !continuous_mode){std::printf(",\"source_momentum\":%.17g",source.momentum);}if(fold_mode){fold_state(phonation.folds);}std::printf("},\"steps\":[");
    Polygons a=first;const auto va=views(a);
    if(joint_tolerance>0){
        ducklips::JointAccuracy controller;controller.prepare(continuous);
        const ducklips::JointAccuracyOptions options{joint_tolerance,joint_maximum_step,1e-11,4096};
        // Provider scratch is overwritten on every trial. Geometry is a
        // function of the requested times, never the last attempted step.
        auto provider=[&](double begin,double finish,ducklips::JointAdvanceInputs& request) noexcept {
            const double t0=begin/duration,t1=finish/duration,mid=.5*(t0+t1),dt=finish-begin;
            sample(first,last,t0,a);
            if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip,epi) || !metrics(geometry,ports,loads,old_volume,old_compliance,old_conductance,lip,epi,t0)){return false;}
            sample(first,last,mid,a);
            if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip,epi) || !metrics(geometry,ports,loads,scratch,c1,g1,lip,epi,mid)){return false;}
            for(std::size_t j=0;j<midpoint_thermal.size();++j){const auto& section=epi.sections[j];midpoint_thermal[j]={c1[thermal_base+j],std::lerp(section.first_area,section.last_area,mid)/(2*thermal_width)};}
            if(viscous_mode && !set_viscous_geometry(mid,g1)){return false;}midpoint_viscous=viscous_geometry;
            sample(first,last,t1,a);
            if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip,epi) || !metrics(geometry,ports,loads,volume,compliance,gend,lip,epi,t1)){return false;}
            for(std::size_t j=0;j<total_nodes;++j){inflow[j]=-(volume[j]-old_volume[j])/dt;}
            for(std::size_t j=0;j<thermal_geometry.size();++j){const auto& section=epi.sections[j];thermal_geometry[j]={compliance[thermal_base+j],std::lerp(section.first_area,section.last_area,t1)/(2*thermal_width)};}
            if(viscous_mode && !set_viscous_geometry(t1,gend)){return false;}
            request.geometry={old_compliance,c1,compliance,old_conductance,g1,gend,midpoint_thermal,thermal_geometry,midpoint_viscous,viscous_geometry,radiation_mode?radiation_area(mid):0.,radiation_mode?radiation_area(t1):0.};
            const auto& section=epi.sections.back();request.outlet_inertance=ducklips::AirProperties{}.density*.5*std::lerp(section.first_length,section.last_length,t1)/std::lerp(section.first_area,section.last_area,t1);
            request.reservoir_pressure=reservoir_pressure;request.nodal_flow=inflow;request.series_pressure=series;return true;
        };
        const auto callback=[](void* context,double begin,double finish,ducklips::JointAdvanceInputs& request) noexcept {return (*static_cast<decltype(provider)*>(context))(begin,finish,request);};
        std::size_t accepted=0,rejected=0;
        for(std::size_t step=1;step<cuts.size();++step){
            const double begin=cuts[step-1]*duration,finish=cuts[step]*duration;
            if(!controller.advance(continuous,begin,finish,callback,&provider,options)){std::fprintf(stderr,"joint accuracy failed at interval %zu reason=%d time=%.17g trial=%.17g difference=%.17g\n",step,static_cast<int>(controller.failure),controller.attempt_time,controller.attempt_step,controller.attempt_error);return 1;}
            const auto& r=controller.last;accepted+=r.accepted_steps;rejected+=r.rejected;
            std::printf("%s{\"dt\":%.17g,\"geometry\":",step>1?",":"",finish-begin);geometry_json(cuts[step],geometry,compact_geometry?std::span<const PlanarContour>{}:va);snapshot();if(!continuous_state(continuous)){return 1;}
            std::printf(",\"accuracy\":{\"energy_change\":%.17g,\"external_work\":%.17g,\"physical_loss\":%.17g,\"numerical_loss\":%.17g,\"energy_residual\":%.17g,\"source_port_work\":%.17g,\"graph_port_work\":%.17g,\"reservoir_volume\":%.17g,\"maximum_accepted_error\":%.17g,\"minimum_step\":%.17g,\"trials\":%zu,\"rejected\":%zu,\"accepted_steps\":%zu",r.energy_change,r.external_work,r.physical_loss,r.numerical_loss,r.energy_residual,r.source_port_work,r.graph_port_work,r.reservoir_volume,r.maximum_accepted_error,r.minimum_step,r.trials,r.rejected,r.accepted_steps);
            values("average_pressure",controller.average_pressure);values("average_flow",controller.average_flow);values("accepted_end_times",controller.accepted_end_times());std::printf("}}");
        }
        std::printf("]}\n");std::fprintf(stderr,"accepted_steps=%zu rejected=%zu; coupled accuracy control; no speech\n",accepted,rejected);return 0;
    }
    const double offset=joint_relaxation?0.:.28867513459481288225;
    for(std::size_t step=1;step<cuts.size();++step){
        const double t0=cuts[step-1],t1=cuts[step],width=t1-t0,dt=duration*width;
        const double first_stage=t0+(.5-offset)*width,second_stage=t0+(.5+offset)*width;
        const double radiation_first_g=radiation_mode?gend[radiation_branch]:0.;
        old_volume=volume;old_compliance=compliance;old_conductance=gend;
        std::printf("%s{\"dt\":%.17g,\"stage1\":",step>1?",":"",dt);
        sample(first,last,first_stage,a);
        if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip,epi)){return 1;}
        if(!metrics(geometry,ports,loads,scratch,c1,g1,lip,epi,first_stage)){return 1;}geometry_json(first_stage,geometry,{});
        std::printf(",\"stage2\":");sample(first,last,second_stage,a);
        if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip,epi)){return 1;}
        if(!metrics(geometry,ports,loads,scratch,c2,g2,lip,epi,second_stage)){return 1;}geometry_json(second_stage,geometry,{});
        std::printf(",\"geometry\":");sample(first,last,t1,a);
        if(!geometry.build_regions(a[0],std::span<const PlanarContour>(va).subspan(1),depth,selected) || !compatible(geometry,cells,ports,loads,lip,epi)){return 1;}
        if(!metrics(geometry,ports,loads,volume,compliance,gend,lip,epi,t1)){return 1;}geometry_json(t1,geometry,compact_geometry?std::span<const PlanarContour>{}:va);
        for(std::size_t i=0;i<total_nodes;++i){inflow[i]=-(volume[i]-old_volume[i])/dt;}
        // A closed passage has no kinetic state to accumulate pressure history.
        // The zero-metric momentum is discarded, as in CanonicalSpeech.
        for(std::size_t i=0;i<channels;++i){if(g1[i]==0 && g2[i]==0){momentum[i]=0;}}
        if(continuous_mode){
            const auto& e=epi.sections.back();
            const double column=ducklips::AirProperties{}.density*.5*std::lerp(e.first_length,e.last_length,t1)/std::lerp(e.first_area,e.last_area,t1);
            ducklips::RadiationGeometry radiation_geometry;
            if(radiation_mode){radiation_geometry={radiation_branch,radiation_first_g,gend[radiation_branch],radiation_area(.5*(t0+t1)),radiation_area(t1)};}
            for(std::size_t i=0;i<thermal_geometry.size();++i){const auto& section=epi.sections[i];thermal_geometry[i]={compliance[thermal_base+i],std::lerp(section.first_area,section.last_area,t1)/(2*thermal_width)};}
            if(joint_relaxation){
                const double midpoint=.5*(t0+t1);
                for(std::size_t i=0;i<midpoint_thermal.size();++i){const auto& section=epi.sections[i];midpoint_thermal[i]={c1[thermal_base+i],std::lerp(section.first_area,section.last_area,midpoint)/(2*thermal_width)};}
                if(viscous_mode && !set_viscous_geometry(midpoint,g1)){return 1;}midpoint_viscous=viscous_geometry;
            }
            if(viscous_mode && !set_viscous_geometry(t1,gend)){std::fprintf(stderr,"lip viscous path disappeared or failed at step %zu\n",step);return 1;}
            const ducklips::JointGraphGeometry whole{old_compliance,c1,compliance,old_conductance,g1,gend,midpoint_thermal,thermal_geometry,midpoint_viscous,viscous_geometry,radiation_geometry.midpoint_area,radiation_geometry.last_area};
            const bool success=joint_relaxation?continuous.advance_joint(whole,dt,reservoir_pressure,{},inflow,series,column):
                continuous.advance(c1,c2,g1,g2,dt,reservoir_pressure,{},inflow,series,ducklips::ContinuousGlottalUpdate::relaxation,column,radiation_mode?&radiation_geometry:nullptr,thermal_geometry,viscous_geometry);
            if(!success){
                std::fprintf(stderr,"continuous source solve failed at step %zu, dt %.17g\n",step,dt);return 1;
            }
            if(!continuous_state(continuous)){return 1;}
            if(joint_relaxation){const auto& r=continuous.last;
                std::printf(",\"joint_passive\":{\"energy_change\":%.17g,\"energy_residual\":%.17g,\"air_residual\":%.17g,\"numerical_dissipation\":%.17g,\"graph_geometry_work\":%.17g,\"external_work\":%.17g,\"viscous_dissipation\":%.17g,\"thermal_dissipation\":%.17g,\"radiated_energy\":%.17g,\"viscous_geometry_work\":%.17g,\"thermal_geometry_work\":%.17g,\"radiation_geometry_work\":%.17g}",r.passive_energy_change,r.passive_energy_residual,r.passive_air_residual,r.passive_numerical_dissipation,r.graph_geometry_work,r.graph_external_work,r.viscous_dissipation,r.thermal_dissipation,r.radiated_energy,r.viscous_geometry_work,r.thermal_geometry_work,r.radiation_geometry_work);
            }else{
            if(viscous_mode){const auto& r=continuous.last;
                std::printf(",\"viscous_dissipation\":%.17g,\"viscous_geometry_work\":%.17g,\"viscous_work_residual\":%.17g",r.viscous_dissipation,r.viscous_geometry_work,r.viscous_energy_residual);
            }
            if(thermal_width>0){const auto& r=continuous.last;
                std::printf(",\"thermal_dissipation\":%.17g,\"thermal_geometry_work\":%.17g,\"thermal_work_residual\":%.17g,\"thermal_air_residual\":%.17g",r.thermal_dissipation,r.thermal_geometry_work,r.thermal_energy_residual,r.thermal_air_residual);
            }
            if(radiation_mode){const auto& r=continuous.last;
                std::printf(",\"radiated_energy\":%.17g,\"radiation_geometry_work\":%.17g,\"radiation_work_residual\":%.17g",r.radiated_energy,r.radiation_geometry_work,r.radiation_energy_residual);
            }
            }
            const auto& r=continuous.source.last;
            std::printf(",\"source_flow\":%.17g,\"source_pressure\":%.17g,\"source_impedance\":%.17g",continuous.last.flow,continuous.last.pressure,continuous.last.load_slope);
            std::printf(",\"source_reservoir_flow\":%.17g,\"source_port_work\":%.17g,\"source_geometry_work\":%.17g,\"source_work_residual\":%.17g,\"source_port_residual\":%.17g",r.flow[0],r.port_work,r.geometry_work,r.energy_residual,continuous.last.port_residual);
        }
        else if(orifice_mode){
            const double midpoint=.5*(t0+t1);const auto& e=epi.sections.back();
            const double inlet_inertance=ducklips::AirProperties{}.density*.5*std::lerp(e.first_length,e.last_length,midpoint)/std::lerp(e.first_area,e.last_area,midpoint);
            const bool advanced=fold_mode?
                phonation.advance(c1,c2,g1,g2,dt,glottal_area,tension,reservoir_pressure,inlet_inertance,surface_pressure,ducklips::AirProperties{},inflow,series):
                source.advance(c1,c2,g1,g2,dt,glottal_area,reservoir_pressure,inlet_inertance,ducklips::AirProperties{},inflow,series);
            if(!advanced){
                std::fprintf(stderr,"glottal source solve failed\n");return 1;
            }
            std::printf(",\"source_flow\":%.17g,\"source_pressure\":%.17g,\"source_impedance\":%.17g",source.flow,source.load_pressure,source.load_impedance);
            std::printf(",\"source_momentum\":%.17g,\"source_inertance\":%.17g",source.momentum,inlet_inertance);
            if(fold_mode){
                fold_state(phonation.folds);
                std::printf(",\"glottal_midpoint_area\":%.17g,\"transglottal_pressure\":%.17g,\"coupling_residual\":%.17g,\"coupling_iterations\":%zu",
                    phonation.midpoint_area,phonation.transglottal_pressure,phonation.coupling_residual,phonation.iterations);
            }
        }
        else if(lip_mode){if(!graph.advance(c1,c2,g1,g2,dt,inflow,series)){std::fprintf(stderr,"shared-lip solve failed\n");return 1;}}
        else{acoustic.advance(c1,c2,g1,g2,dt,inflow,series);}
        for(std::size_t i=0;i<channels;++i){if(gend[i]==0){momentum[i]=0;}}
        values("wall_flow",inflow);snapshot();std::putchar('}');
    }
    std::puts("]}");std::fprintf(stderr,"steps=%zu; inviscid canonical moving-volume reference; no speech\n",cuts.size()-1);
}
