// Export the current native geometry for the constrained-tissue experiment.
// Arguments: circle centre x/y (cm), gesture place from lips (cm), target gap (cm).
#include "ducklips/tongue.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv) {
    if (argc!=5) { std::fprintf(stderr,"usage: tongue_geometry_probe centre_x centre_y place gap\n");return 2; }
    const double x=std::strtod(argv[1],nullptr),y=std::strtod(argv[2],nullptr);
    const double place=std::strtod(argv[3],nullptr),gap=std::strtod(argv[4],nullptr);
    if (!std::isfinite(x+y+place+gap) || place<=0 || place>=17.5 || gap<0) { return 2; }
    ducklips::TonguePosture posture{};posture.body_x=x;posture.body_y=y+.9*posture.jaw;
    constexpr std::size_t count=351;
    std::array<double,count> arcs{},gaps{};
    for (std::size_t i=0;i<count;++i) { arcs[i]=17.5*static_cast<double>(i)/(count-1); }
    const double arc=ducklips::MidsagittalTongue::length-place;
    const ducklips::Point2 mid=ducklips::MidsagittalTongue::midline(arc),normal=ducklips::MidsagittalTongue::normal(arc);
    const double offset=ducklips::MidsagittalTongue::outer_offset(arc);
    std::printf("{\"centre_cm\":[%.12g,%.12g],\"radius_cm\":%.12g,\"place_cm\":%.12g,\"gap_cm\":%.12g,\"wall_target_cm\":[%.12g,%.12g],\"normal\":[%.12g,%.12g],\"sections\":[",x,y,ducklips::MidsagittalTongue::body_radius,place,gap,mid.x+normal.x*offset,mid.y+normal.y*offset,normal.x,normal.y);
    for (std::size_t i=0;i<count;++i) {
        const ducklips::Point2 m=ducklips::MidsagittalTongue::midline(arcs[i]),n=ducklips::MidsagittalTongue::normal(arcs[i]);
        const double o=ducklips::MidsagittalTongue::outer_offset(arcs[i]);
        std::printf("%s[%.12g,%.12g,%.12g,%.12g,%.12g,%.12g]",i?",":"",arcs[i],m.x+n.x*o,m.y+n.y*o,n.x,n.y,o+ducklips::MidsagittalTongue::floor_offset(arcs[i],posture));
    }
    std::printf("],\"poses\":[");
    for (std::size_t step=0;step<=10;++step) {
        const double activation=static_cast<double>(step)/10;
        const ducklips::TongueGesture body{.activation=activation,.place_from_lips=place,.gap=gap};
        const ducklips::TongueShape s=ducklips::MidsagittalTongue::pose(posture,body,{});
        ducklips::MidsagittalTongue::gaps(s,posture,arcs,gaps);
        std::printf("%s{\"activation\":%.12g,\"centre_cm\":[%.12g,%.12g],\"blade_root_cm\":[%.12g,%.12g],\"tip_cm\":[%.12g,%.12g],\"gaps_cm\":[",step?",":"",activation,s.body_centre.x,s.body_centre.y,s.blade_root.x,s.blade_root.y,s.tip.x,s.tip.y);
        for (std::size_t i=0;i<count;++i) { std::printf("%s%.12g",i?",":"",gaps[i]); }
        std::printf("]}");
    }
    std::printf("]}\n");
}
