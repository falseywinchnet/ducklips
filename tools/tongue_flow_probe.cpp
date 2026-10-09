// CSV interface for independent geometry checks, not a speech renderer.
#include "ducklips/tongue_flow.hpp"
#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv) {
    if (argc!=7) { std::fprintf(stderr,"usage: tongue_flow_probe profile.csv anchor_y free_y centre_x command steps < points.csv\n");return 2; }
    FILE* file=std::fopen(argv[1],"r");if (!file) { return 2; }
    std::array<ducklips::TongueLiftKnot,ducklips::ConservativeTongueMap::maximum_knots> knots{};std::size_t count=0;
    while (count<knots.size() && std::fscanf(file,"%lf,%lf,%lf",&knots[count].x,&knots[count].lift,&knots[count].slope)==3) { ++count; }
    int extra=std::fgetc(file);
    while (extra!=EOF && std::isspace(static_cast<unsigned char>(extra))) { extra=std::fgetc(file); }
    std::fclose(file);if (extra!=EOF) { return 2; }
    ducklips::ConservativeTongueMap map{};
    if (!map.configure(std::span<const ducklips::TongueLiftKnot>(knots.data(),count),std::strtod(argv[2],nullptr),std::strtod(argv[3],nullptr),std::strtod(argv[4],nullptr))) { return 2; }
    const double command=std::strtod(argv[5],nullptr);const std::size_t steps=std::strtoul(argv[6],nullptr,10);
    std::puts("x,y,X,Y,Fxx,Fxy,Fyx,Fyy,valid");ducklips::Point2 point{};
    int scanned=0;
    while ((scanned=std::scanf("%lf,%lf",&point.x,&point.y))==2) {
        const ducklips::TongueMaterialPoint result=map.map(point,command,steps);
        std::printf("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",point.x,point.y,result.position.x,result.position.y,result.gradient[0],result.gradient[1],result.gradient[2],result.gradient[3],result.valid?1:0);
        if (!result.valid) { return 1; }
    }
    return scanned==EOF?0:2;
}
