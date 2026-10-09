#include "ducklips/planar_face_motion.hpp"
#include "ducklips/cavity_network.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
namespace {std::atomic<std::size_t> allocations{};unsigned failures=0;}
void* operator new(std::size_t size){++allocations;if(void* p=std::malloc(size)){return p;}throw std::bad_alloc();}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
namespace {
using ducklips::Point2;using ducklips::PlanarContour;
std::array<Point2,4> box(double x0,double x1,double y0,double y1){return {{{x0,y0},{x1,y0},{x1,y1},{x0,y1}}};}
void check(bool value,const char* name,double error=0){std::printf("%s %s %.12g\n",value?"PASS":"FAIL",name,error);if(!value){++failures;}}
double cubic_primitive(double t){return .5*t*t-4.5*t-3.24*std::log(t+1)+15.4*std::log(t+2);}
double closed_duration(const ducklips::PlanarFaceMotion& motion){double result=0;for(const auto& w:motion.windows){if(w.average_area==0){result+=w.last-w.first;}}return result;}
}
int main(){
    const auto domain=box(-3,4,-2,3);const ducklips::PlanarFace face{true,0,0,1};
    ducklips::PlanarFaceMotion motion;motion.prepare(64,4096);
    const auto a=box(-2,-1,-1,2),b=box(1,2,-1,2);const std::array<PlanarContour,1> first{a},last{b};
    const std::size_t before=allocations.load();bool valid=motion.build(domain,first,domain,last,face,1);const std::size_t after=allocations.load();
    check(valid && std::abs(motion.average_area-2./3)<1e-14 && std::abs(closed_duration(motion)-1./3)<1e-14,
          "passage open at both endpoints retains its intervening closure",motion.average_area-2./3);
    check(before==after,"moving face partition and integration allocate no storage");
    const auto narrow0=box(-.37,-.37+1e-7,-1,2),narrow1=box(.63,.63+1e-7,-1,2);
    const std::array<PlanarContour,1> n0{narrow0},n1{narrow1};valid=motion.build(domain,n0,domain,n1,face,1);
    check(valid && std::abs(closed_duration(motion)-1e-7)<1e-14,"events retain a closure much shorter than a uniform sampling step",closed_duration(motion)-1e-7);
    // N/D intersection: y=2/(2+t), with solid below, free face [y,2].
    const std::array<Point2,4> r0{{{-1,-1},{2,-1},{1,2},{-1,0}}},r1{{{-1,-1},{3,-1},{2,2},{-1,0}}};
    const std::array<PlanarContour,1> rational0{r0},rational1{r1};const double exact=2-2*std::log(1.5);
    valid=motion.build(domain,rational0,domain,rational1,{true,0,0,2},1,1e-13);
    check(valid && std::abs(motion.average_area-exact)<1e-13,"changing edge direction integrates the rational aperture",motion.average_area-exact);
    valid=motion.build(domain,rational1,domain,rational0,{true,0,0,2},1,1e-13);
    check(valid && std::abs(motion.average_area-exact)<1e-13,"time reversal preserves the aperture integral",motion.average_area-exact);
    // Two solids' edges exchange ownership of the union envelope at t=1/2.
    const auto up0=box(-1,1,-1,.2),up1=box(-1,1,-1,.8),down0=box(-1,1,.8,2),down1=box(-1,1,.2,2);
    const std::array<PlanarContour,2> overlap0{up0,down0},overlap1{up1,down1};
    valid=motion.build(domain,overlap0,domain,overlap1,face,1);
    check(valid && std::abs(motion.average_area-.15)<1e-14 && std::abs(closed_duration(motion)-.5)<1e-14,
          "overlapping solids close at their edge-order event",motion.average_area-.15);
    // A pressure source starts after this shutter closes. Whole-step average
    // conductance would invent transfer; consuming the windows in order must
    // leave the receiving chamber unchanged.
    valid=motion.build(domain,overlap0,domain,overlap1,face,1);
    ducklips::CavityNetwork pressure;
    const std::array<ducklips::CavityConnection,1> connection{{{0,1}}};pressure.prepare(2,connection);
    const ducklips::AirProperties air{};const double bulk=air.density*air.sound_speed*air.sound_speed;
    const std::array<double,2> volume{1e-6,1e-6},compliance{1e-6/bulk,1e-6/bulk},zero{};
    valid=pressure.reset(volume,compliance,zero)&&valid;double input_volume=0,transferred=0;
    for(const auto& w:motion.windows){
        const double dt=(w.last-w.first)*.001;const std::array<double,2> input{w.average_area==0?1e-8:0,0};
        const double conductance=w.average_area*1e-7;
        const std::array<ducklips::CavityPassage,1> passage{{{w.average_area,conductance>0?1/conductance:0,0,0,0}}};
        valid=pressure.advance(volume,compliance,passage,input,dt,air).converged&&valid;
        input_volume+=dt*input[0];transferred+=dt*pressure.flows()[0];
    }
    check(valid && pressure.pressures()[0]>0 && pressure.pressures()[1]==0 && transferred==0,
          "chronological opening windows block a pressure source after closure",transferred);
    check(std::abs(compliance[0]*pressure.pressures()[0]-input_volume)<1e-24,
          "closed chamber retains all subsequently supplied air",compliance[0]*pressure.pressures()[0]-input_volume);
    pressure.reset(volume,compliance,zero);
    const double mean_g=motion.average_area*1e-7;
    const std::array<ducklips::CavityPassage,1> averaged{{{motion.average_area,1/mean_g,0,0,0}}};
    const std::array<double,2> averaged_input{input_volume/.001,0};
    valid=pressure.advance(volume,compliance,averaged,averaged_input,.001,air).converged;
    check(valid && pressure.pressures()[1]>0 && pressure.flows()[0]>0,
          "whole-step averaging demonstrably transfers pressure through later contact",pressure.pressures()[1]);
    // Two line intersections have difference
    // (t-.2)(t-.5)(t-.8)/((t+1)(t+2)). All three cubic roots matter;
    // the opening is positive only on (.2,.5) and (.8,1).
    const auto cubic_domain=box(-3,4,-11,3);
    const std::array<Point2,4> upper0{{{-1,0},{0,-4.24},{0,2},{-1,2}}},upper1{{{-1,1},{1,-4.24},{1,2},{-1,2}}};
    const std::array<Point2,4> lower0{{{-1,-10},{1,-10},{1,-8.4},{-1,0}}},lower1{{{-1,-10},{2,-10},{2,-4.9},{-1,0}}};
    const std::array<PlanarContour,2> cubic0{upper0,lower0},cubic1{upper1,lower1};
    valid=motion.build(cubic_domain,cubic0,cubic_domain,cubic1,{true,0,-5,0},1);
    const double cubic_exact=cubic_primitive(.5)-cubic_primitive(.2)+cubic_primitive(1)-cubic_primitive(.8);
    check(valid && std::abs(motion.average_area-cubic_exact)<1e-13 && std::abs(closed_duration(motion)-.5)<1e-12,
          "three cubic edge-order events retain both separate open intervals",motion.average_area-cubic_exact);
    auto rd=domain,ra=a,rb=b;for(auto* p:{&rd,&ra,&rb}){for(Point2& v:*p){std::swap(v.x,v.y);}}
    const std::array<PlanarContour,1> rot0{ra},rot1{rb};valid=motion.build(rd,rot0,rd,rot1,{false,0,0,1},1);
    check(valid && std::abs(motion.average_area-2./3)<1e-14,"horizontal and vertical faces share the same event geometry");
    const auto wall=box(0,1,-1,2);const std::array<PlanarContour,1> seal{wall};
    valid=motion.build(domain,seal,domain,seal,face,1);
    check(valid && motion.average_area==0,"a solid boundary exactly on the face remains sealed");
    valid=motion.build(wall,{},wall,{},face,1);
    check(valid && motion.average_area==0,"a domain boundary has no two-sided fluid opening");
    auto reverse_a=a,reverse_b=b;std::reverse(reverse_a.begin(),reverse_a.end());std::reverse(reverse_b.begin(),reverse_b.end());
    const std::array<PlanarContour,1> rev0{reverse_a},rev1{reverse_b};valid=motion.build(domain,rev0,domain,rev1,face,1);
    check(valid && std::abs(motion.average_area-2./3)<1e-14,"contour winding does not change interval openings");
    check(!motion.build(domain,first,domain,{},face,1),"mismatched contours fail");
    check(!motion.build(domain,first,domain,last,{true,0,1,0},1),"invalid face fails");
    check(!motion.build(domain,first,domain,last,face,1,std::numeric_limits<double>::quiet_NaN()),"invalid tolerance fails");
    ducklips::PlanarFaceMotion bounded;bounded.prepare(64,2);
    check(!bounded.build(domain,first,domain,last,face,1),"event capacity exhaustion fails explicitly");
    bounded.prepare(3,16);check(!bounded.build(domain,{},domain,{},face,1),"edge capacity exhaustion fails explicitly");
    return failures?1:0;
}
