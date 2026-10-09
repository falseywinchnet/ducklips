#include "ducklips/glottal_contour.hpp"
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
namespace {
void require(bool value,const char* message){if(!value){throw std::runtime_error(message);}}
void output(const ducklips::GlottalContourGeometry& g){
    std::printf("{\"accepted\":true,\"chambers\":[");
    for(std::size_t i=0;i<2;++i){const auto& c=g.chamber[i];std::printf("%s{\"volume\":%.17g,\"centroid\":%.17g,\"derivative\":[%.17g,%.17g]}",i?",":"",c.volume,c.centroid,c.area_derivative[0],c.area_derivative[1]);}
    std::printf("],\"aperture\":[%.17g,%.17g,%.17g],\"passages\":[",g.aperture[0],g.aperture[1],g.aperture[2]);
    for(std::size_t i=0;i<3;++i){const auto& p=g.passage[i];std::printf("%s{\"open\":%s,\"inverse_area\":[%.17g,%.17g],\"inverse_area_cubed\":[%.17g,%.17g]}",i?",":"",p.open?"true":"false",p.inverse_area[0],p.inverse_area[1],p.inverse_area_cubed[0],p.inverse_area_cubed[1]);}std::puts("]}");
}
}
int main(int argc,char** argv){try{
    if(argc==2 && std::string(argv[1])=="--batch"){
        ducklips::GlottalContour p;
        while(std::scanf("%lf %lf %lf %lf",&p.lower_thickness,&p.upper_thickness,&p.lower_area,&p.upper_area)==4){ducklips::GlottalContourGeometry g;if(p.geometry(g)){output(g);}else{std::puts("{\"accepted\":false}");}}return 0;
    }
    require(argc==1,"unknown arguments");
    ducklips::GlottalContour p{.0025,.0005,3e-6,3e-6};ducklips::GlottalContourGeometry g;
    require(p.geometry(g),"uniform geometry");require(std::abs(g.chamber[0].volume-7.5e-9)<1e-23 && std::abs(g.chamber[1].volume-1.5e-9)<1e-23,"uniform volume");
    require(g.passage[0].open && g.passage[1].open && g.passage[2].open,"uniform openings");
    p.upper_area=0.;require(p.geometry(g),"upper contact geometry");
    require(g.chamber[1].volume>0 && std::abs(g.chamber[1].volume-6.25e-11)<1e-24,"contact erased wedge storage");
    require(!g.passage[2].open && g.aperture[2]==0 && g.passage[1].open,"upper contact connectivity");
    p.upper_area=-1e-6;require(p.geometry(g),"fully occupied upper cell");
    require(g.chamber[1].volume==0 && !g.passage[1].open && !g.passage[2].open && g.chamber[0].volume>0,"empty-cell geometry");
    p.lower_area=0;p.upper_area=3e-6;ducklips::GlottalContourMeasure plateau;
    require(p.measure(0.,.5*p.lower_thickness,plateau) && plateau.volume==0 && plateau.area_derivative[0]==0 && plateau.area_derivative[1]==0,"zero plateau acquired pressure-force area");
    p.lower_area=0;p.upper_area=0;require(p.geometry(g) && g.chamber[0].volume==0 && g.chamber[1].volume==0,"fully closed profile");
    p.upper_thickness=std::numeric_limits<double>::quiet_NaN();require(!p.geometry(g) && g.chamber[0].volume==0,"invalid profile changed result");
    std::puts("PASS uniform storage, wedge retained at contact, exact seal, empty-cell geometry, rejection");return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
