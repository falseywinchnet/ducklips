// Reconstruct selected native speech poses and measure an explicit extruded
// planar air domain. This is a geometry reference, not a speech renderer.
#include "ducklips/planar_airspace.hpp"
#include "ducklips/tongue_surface.hpp"
#include "ducklips/physics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>

namespace {
void points(std::ostream& out,ducklips::PlanarContour p){out<<'[';for(std::size_t i=0;i<p.size();++i){if(i){out<<',';}out<<'['<<p[i].x<<','<<p[i].y<<']';}out<<']';}
void numbers(std::ostream& out,std::span<const double> p){out<<'[';for(std::size_t i=0;i<p.size();++i){if(i){out<<',';}out<<p[i];}out<<']';}
}

int main(int argc,char** argv){
    if(argc!=6 && argc!=7 && argc!=8){std::fprintf(stderr,"usage: ducklips_airspace_geometry_probe GEOMETRY.csv START END CONTOUR_VERTICES OUTPUT.csv [LAST_POSE.json | --frames DIRECTORY]\n");return 2;}
    const bool frames=argc==8;
    if(frames && std::string(argv[6])!="--frames"){return 2;}
    if(frames){std::filesystem::create_directories(argv[7]);}
    double start=0,end=0;std::size_t vertices=0;
    try{start=std::stod(argv[2]);end=std::stod(argv[3]);vertices=std::stoul(argv[4]);}catch(...){return 2;}
    if(!std::isfinite(start) || !std::isfinite(end) || start<0 || end<start || vertices<32 || vertices>1024){return 2;}
    std::ifstream input(argv[1]);std::ofstream output(argv[5]);if(!input || !output){return 2;}
    constexpr std::size_t nx=64,ny=64;
    std::array<double,nx+1> xs{};std::array<double,ny+1> ys{};
    for(std::size_t i=0;i<=nx;++i){xs[i]=-.01+.11*double(i)/nx;}
    for(std::size_t i=0;i<=ny;++i){ys[i]=.105*double(i)/ny;}
    ducklips::PlanarAirspace cells;cells.prepare(xs,ys,2*vertices+600);
    std::vector<ducklips::Point2> body(vertices),pad(vertices),domain;domain.reserve(600);
    std::array<ducklips::Point2,4> blade{};
    const std::array<ducklips::PlanarContour,3> solids{{body,blade,pad}};
    std::vector<double> arcs;for(std::size_t i=0;i<=256;++i){arcs.push_back(17.5*double(i)/256);}
    for(double arc:{4.,6.5,7.,7.+ducklips::pi,7.+ducklips::pi+.5,7.+ducklips::pi+2.5,
                    7.+ducklips::pi+4.,7.+ducklips::pi+5.2,7.+ducklips::pi+6.,7.+ducklips::pi+6.4}){arcs.push_back(arc);}
    std::sort(arcs.begin(),arcs.end());arcs.erase(std::unique(arcs.begin(),arcs.end()),arcs.end());
    std::string line;std::getline(input,line);
    output<<std::setprecision(17)<<"time,air_volume_m3,body_polygon_area_m2,maximum_pose_error_cm,maximum_cell_volume_m3,open_vertical_faces,open_horizontal_faces\n";
    std::size_t rows=0;
    while(std::getline(input,line)){
        std::stringstream stream(line);std::string item;std::vector<double> values;std::size_t column=0;
        while(std::getline(stream,item,',')){values.push_back(column++==1?0:std::stod(item));}
        if(values.size()<26){std::fprintf(stderr,"geometry CSV must include posture_root\n");return 2;}const double time=values[0];if(time<start || time>end){continue;}
        const ducklips::TonguePosture posture{.body_x=values[10],.body_y=values[11],.jaw=values[12],.blade=values[13],
                                               .root=values[25],.offset_x=values[14],.offset_y=values[15]};
        const ducklips::TongueGesture dorsal{.activation=values[2],.place_from_lips=values[3],.gap=values[4],.width=values[5]};
        const ducklips::TongueGesture apical{.activation=values[6],.place_from_lips=values[7],.gap=values[8],.tip_down=values[9]};
        const ducklips::Point2 centre{posture.body_x+posture.offset_x,posture.body_y-.9*std::clamp(posture.jaw,0.,1.)+posture.offset_y};
        ducklips::ConservativeTongueSurface surface;if(!surface.configure(centre,dorsal)){return 1;}
        const ducklips::TongueShape tongue=ducklips::MidsagittalTongue::pose(posture,dorsal,apical,&surface);if(!tongue.valid){return 1;}
        const double error=std::max({std::abs(tongue.blade_root.x-values[16]),std::abs(tongue.blade_root.y-values[17]),
                                    std::abs(tongue.tip.x-values[18]),std::abs(tongue.tip.y-values[19])});
        if(error>1e-8){std::fprintf(stderr,"reconstructed pose differs at %.6f: %.12g cm\n",time,error);return 1;}
        for(std::size_t i=0;i<vertices;++i){
            const double angle=2*ducklips::pi*double(i)/double(vertices),c=std::cos(angle),s=std::sin(angle);
            const ducklips::TongueMaterialPoint point=surface.material({centre.x+2.7*c,centre.y+2.7*s});if(!point.valid){return 1;}
            body[i]={.01*point.position.x,.01*point.position.y};pad[i]={.01*(tongue.tip.x+.35*c),.01*(tongue.tip.y+.35*s)};
        }
        // Explicit diagnostic solidification of the former zero-thickness
        // blade. Its lower attachment is an interior BODY MATERIAL point,
        // transported by the same map as the upper attachment. Extending a
        // thickness normal outward from the upper root only touches the body
        // at one point, leaving an unphysical air route through the attachment.
        // This is not a muscle/constant-volume model of the combined blade.
        const double dx=tongue.tip.x-tongue.blade_root.x,dy=tongue.tip.y-tongue.blade_root.y,length=std::hypot(dx,dy);
        if(length<=0){return 1;}const ducklips::Point2 down{dy/length,-dx/length};
        const ducklips::TongueShape rest=ducklips::MidsagittalTongue::pose(posture,{},{});
        const ducklips::TongueMaterialPoint lower=surface.material({rest.blade_root.x,rest.blade_root.y-ducklips::MidsagittalTongue::blade_root_thickness});
        if(!lower.valid){return 1;}
        blade={{{.01*tongue.blade_root.x,.01*tongue.blade_root.y},{.01*tongue.tip.x,.01*tongue.tip.y},
                 {.01*(tongue.tip.x+.35*down.x),.01*(tongue.tip.y+.35*down.y)},
                 {.01*lower.position.x,.01*lower.position.y}}};
        domain.clear();
        for(double arc:arcs){const auto m=ducklips::MidsagittalTongue::midline(arc),n=ducklips::MidsagittalTongue::normal(arc);const double h=ducklips::MidsagittalTongue::outer_offset(arc);domain.push_back({.01*(m.x+h*n.x),.01*(m.y+h*n.y)});}
        for(std::size_t i=arcs.size();i-->0;){
            const double arc=arcs[i];const auto m=ducklips::MidsagittalTongue::midline(arc),n=ducklips::MidsagittalTongue::normal(arc);
            if(arc==7.+ducklips::pi){domain.push_back({.02,.01*(9-1-1.4*std::clamp(posture.jaw,0.,1.))});}
            const double h=ducklips::MidsagittalTongue::floor_offset(arc,posture);domain.push_back({.01*(m.x-h*n.x),.01*(m.y-h*n.y)});
            // Preserve the native vertical lower-incisor step instead of
            // replacing it with a diagonal across one wall sample interval.
            if(arc==7.+ducklips::pi+6.){domain.push_back({.08,.01*(9-1-1.4*std::clamp(posture.jaw,0.,1.))});}
        }
        if(!cells.build(domain,solids,.022)){return 1;}
        double area=0;for(std::size_t i=0;i<vertices;++i){const auto a=body[i],b=body[(i+1)%vertices];area+=a.x*b.y-a.y*b.x;}
        std::size_t vertical=0,horizontal=0;for(double v:cells.vertical_face){vertical+=v>0;}for(double v:cells.horizontal_face){horizontal+=v>0;}
        output<<time<<','<<std::accumulate(cells.volume.begin(),cells.volume.end(),0.)<<','<<std::abs(area)/2<<','<<error<<','<<*std::max_element(cells.volume.begin(),cells.volume.end())<<','<<vertical<<','<<horizontal<<'\n';++rows;
        if(argc>=7){
            std::ostringstream name;
            if(frames){name<<argv[7]<<"/pose-"<<std::fixed<<std::setprecision(3)<<time<<".json";}else{name<<argv[6];}
            std::ofstream snapshot(name.str());if(!snapshot){return 2;}snapshot<<std::setprecision(17)<<"{\"time\":"<<time<<",\"depth\":0.022,\"domain\":";points(snapshot,domain);
            // These are the end caps of the reconstructed tongue/floor domain.
            // Lip shaping and the time-varying fold aperture are separate body
            // controls, so neither cap is itself a complete anatomical port.
            snapshot<<",\"boundary_edges\":{\"supraglottal_entry\":"<<domain.size()-1
                    <<",\"mouth_end\":"<<arcs.size()-1<<'}';
            snapshot<<",\"solids\":[";for(std::size_t i=0;i<solids.size();++i){if(i){snapshot<<',';}points(snapshot,solids[i]);}
            snapshot<<"],\"x_edges\":";numbers(snapshot,xs);snapshot<<",\"y_edges\":";numbers(snapshot,ys);
            snapshot<<",\"volume\":";numbers(snapshot,cells.volume);snapshot<<",\"vertical_face\":";numbers(snapshot,cells.vertical_face);
            snapshot<<",\"horizontal_face\":";numbers(snapshot,cells.horizontal_face);snapshot<<"}\n";if(!snapshot){return 1;}
        }
    }
    if(!rows || !output){return 1;}std::fprintf(stderr,"measured %zu native poses with %zu contour vertices\n",rows,vertices);return 0;
}
