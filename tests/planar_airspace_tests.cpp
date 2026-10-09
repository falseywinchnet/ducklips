#include "ducklips/planar_airspace.hpp"
#include "ducklips/physics.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <numeric>
#include <vector>

namespace {std::atomic<std::size_t> allocations{};std::size_t failures{};}
void* operator new(std::size_t size){allocations++;if(void* p=std::malloc(size)){return p;}throw std::bad_alloc();}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
namespace {
void check(bool value,const char* name,double error=0){std::printf("%s %s %.12g\n",value?"PASS":"FAIL",name,error);if(!value){++failures;}}
std::vector<double> grid(std::size_t n,double first,double last){std::vector<double> values(n+1);for(std::size_t i=0;i<=n;++i){values[i]=first+(last-first)*double(i)/double(n);}return values;}
double sum(const std::vector<double>& values){return std::accumulate(values.begin(),values.end(),0.);}
std::array<ducklips::Point2,4> box(double x0,double x1,double y0,double y1){return {{{x0,y0},{x1,y0},{x1,y1},{x0,y1}}};}
double polygon_area(std::span<const ducklips::Point2> p){double result=0;for(std::size_t i=0;i<p.size();++i){const auto a=p[i],b=p[(i+1)%p.size()];result+=a.x*b.y-a.y*b.x;}return std::abs(result)/2;}
}
int main(){
    const auto domain=box(0,.04,0,.04),a=box(.01,.025,.01,.03),b=box(.02,.03,.005,.025);
    const std::array<ducklips::PlanarContour,2> solids{{a,b}};
    ducklips::PlanarAirspace space;space.prepare(grid(8,0,.04),grid(8,0,.04),520);
    const std::size_t before=allocations.load();const bool valid=space.build(domain,solids,.022);
    const std::size_t after=allocations.load();const double expected=(.04*.04-(.015*.02+.01*.02-.005*.015))*.022;
    check(valid && std::abs(sum(space.volume)-expected)<1e-18,"overlapping body parts exclude their union once",sum(space.volume)-expected);
    check(before==after,"planar cell and face integration allocates no render storage");
    ducklips::PlanarAirspace fine;fine.prepare(grid(17,0,.04),grid(13,0,.04),520);fine.build(domain,solids,.022);
    check(std::abs(sum(space.volume)-sum(fine.volume))<1e-18,"grid refinement preserves total polygonal air volume",sum(space.volume)-sum(fine.volume));
    // Polygon intersections inside cells change the union envelope. The
    // triangle pair overlaps in a known smaller triangle (area .0001).
    const std::array<ducklips::Point2,3> t1{{{0,0},{.03,0},{0,.03}}},t2{{{.01,0},{.04,0},{.04,.03}}};
    const std::array<ducklips::PlanarContour,2> triangles{{t1,t2}};
    space.build(domain,triangles,.022);
    const double triangle_expected=(.0016-(.00045+.00045-.0001))*.022;
    check(std::abs(sum(space.volume)-triangle_expected)<1e-18,"crossing solid edges are integrated at their union switch",sum(space.volume)-triangle_expected);
    // A solid slab spans the domain and lies exactly on grid faces. It must
    // isolate left/right fluid without deleting either chamber's volume.
    const auto seal=box(.015,.02,-.01,.05);const std::array<ducklips::PlanarContour,1> slab{{seal}};
    space.build(domain,slab,.022);double leak=0;
    for(std::size_t row=0;row<space.rows();++row){leak=std::max(leak,space.vertical_face[row*9+3]);leak=std::max(leak,space.vertical_face[row*9+4]);}
    check(leak==0 && std::abs(sum(space.volume)-.035*.04*.022)<1e-18,"closed contact separates full-volume chambers exactly",leak);
    const auto horizontal=box(-.01,.05,.015,.02);const std::array<ducklips::PlanarContour,1> flat{{horizontal}};
    space.build(domain,flat,.022);leak=0;for(std::size_t col=0;col<8;++col){leak=std::max({leak,space.horizontal_face[3*8+col],space.horizontal_face[4*8+col]});}
    check(leak==0,"contact lying on a horizontal face seals exactly",leak);
    const std::array<ducklips::Point2,8> concave{{{0,0},{.04,0},{.04,.04},{.03,.04},{.03,.01},{.01,.01},{.01,.04},{0,.04}}};
    space.build(concave,{},.022);
    check(std::abs(sum(space.volume)-.001*.022)<1e-18,"concave mouth domain retains every free interval",sum(space.volume)-.001*.022);
    auto reversed=concave;std::reverse(reversed.begin(),reversed.end());const std::vector<double> concave_volume=space.volume;
    space.build(reversed,{},.022);double winding_error=0;for(std::size_t i=0;i<space.volume.size();++i){winding_error=std::max(winding_error,std::abs(space.volume[i]-concave_volume[i]));}
    check(winding_error<1e-20,"polygon winding does not change air volume",winding_error);
    // Translate a rounded pad across a measurement line while staying inside
    // the mouth. Total free air must remain constant; individual cell changes
    // must tend to zero with the displacement, on both sides of the pad.
    std::vector<ducklips::Point2> pad(256);const double radius=.0035;
    std::vector<double> old(space.volume.size());double worst_total=0,previous_change=1,ratio=0;
    bool continuous=true,lower_present=true;
    for(double displacement:{2e-4,2e-5,2e-6,2e-7}){
        double change=0;
        for(int direction:{-1,1}){
            const double centre=.02-radius+double(direction)*displacement/2;
            for(std::size_t i=0;i<pad.size();++i){const double angle=2*ducklips::pi*double(i)/double(pad.size());pad[i]={centre+radius*std::cos(angle),.025+radius*std::sin(angle)};}
            const std::array<ducklips::PlanarContour,1> obstacle{{pad}};space.build(domain,obstacle,.022);
            worst_total=std::max(worst_total,std::abs(sum(space.volume)-(.0016-polygon_area(pad))*.022));
            lower_present=lower_present && space.volume[2*8+3]>.0049*.0049*.022;
            if(direction<0){old=space.volume;}else{for(std::size_t i=0;i<old.size();++i){change=std::max(change,std::abs(old[i]-space.volume[i]));}}
        }
        ratio=change/previous_change;if(previous_change<1 && ratio>.12){continuous=false;}previous_change=change;
    }
    check(worst_total<1e-18 && lower_present,"moving rounded pad preserves total air above and below it",worst_total);
    check(continuous && previous_change<1e-10,"partial-cell volume change vanishes with tip displacement",ratio);
    // Polygonization itself must converge to the analytical circle area.
    double old_error=1;
    for(std::size_t count:{32,64,128,256}){
        pad.resize(count);for(std::size_t i=0;i<count;++i){const double angle=2*ducklips::pi*double(i)/double(count);pad[i]={.02+radius*std::cos(angle),.02+radius*std::sin(angle)};}
        const std::array<ducklips::PlanarContour,1> obstacle{{pad}};space.build(domain,obstacle,.022);
        const double error=std::abs(sum(space.volume)-(.0016-ducklips::pi*radius*radius)*.022);
        check(error<old_error*.26,"polygon refinement approaches analytical circle air volume",error);old_error=error;
    }
    // Grid cells are control volumes, not an assertion of fluid connectivity.
    // A wall crossing the interior of a cell must create separate regions.
    ducklips::PlanarAirspace split;split.prepare(grid(1,0,1),grid(1,0,1),32);split.prepare_regions(2048,4096);
    const auto unit=box(0,1,0,1),interior_wall=box(.45,.55,-1,2);
    const std::array<ducklips::PlanarContour,1> internal{{interior_wall}};
    const std::size_t regions_before=allocations.load();const bool split_valid=split.build_regions(unit,internal,1);
    const std::size_t regions_after=allocations.load();
    check(split_valid && split.regions.size()==2 && split.region_links.empty()
        && std::abs(split.regions[0].volume-.45)<1e-14 && std::abs(split.regions[1].volume-.45)<1e-14,
        "one grid cell retains two separated air regions");
    check(regions_before==regions_after,"component decomposition performs no allocation after preparation");
    // Two pieces of tissue meet only at one vertex. Fluid cannot pass through
    // that point, even though they do not overlap in positive area.
    const std::array<ducklips::Point2,3> upper{{{.2,1.1},{.8,1.1},{.431,.517}}};
    const std::array<ducklips::Point2,3> lower{{{.2,-.1},{.8,-.1},{.431,.517}}};
    const std::array<ducklips::PlanarContour,2> contact{{upper,lower}};
    const bool contact_valid=split.build_regions(unit,contact,1);
    check(contact_valid && split.regions.size()==2 && split.region_links.empty(),"point contact closes a within-cell air route exactly",double(split.regions.size()));
    const auto narrow_wall=box(.5,.5+1e-12,-1,2);const std::array<ducklips::PlanarContour,1> thin{{narrow_wall}};
    check(split.build_regions(unit,thin,1) && split.regions.size()==2,"thin positive tissue is not removed by a contact threshold");
    // A small positive gap at the contact reconnects the air; all pieces in
    // this one cell now belong to one region, not one pressure per strip.
    auto raised=upper;for(ducklips::Point2& point:raised){point.y+=1e-10;}
    const std::array<ducklips::PlanarContour,2> opened{{raised,lower}};
    check(split.build_regions(unit,opened,1) && split.regions.size()==1,"small positive opening reconnects the cell without an aperture floor");
    ducklips::PlanarAirspace limited;limited.prepare(grid(1,0,1),grid(1,0,1),32);limited.prepare_regions(1,1);
    check(!limited.build_regions(unit,contact,1),"insufficient component capacity is explicit failure");
    for(bool vertical_wall:{false,true}){
        ducklips::PlanarAirspace ports;
        ports.prepare(grid(vertical_wall?1:2,0,1),grid(vertical_wall?2:1,0,1),32);ports.prepare_regions(2048,4096);
        const auto barrier=vertical_wall?box(.45,.55,-1,2):box(-1,2,.45,.55);
        const std::array<ducklips::PlanarContour,1> wall{{barrier}};
        bool ports_valid=ports.build_regions(unit,wall,1) && ports.regions.size()==4 && ports.region_links.size()==2;
        for(const ducklips::PlanarRegionLink& link:ports.region_links){
            const auto a=ports.regions[link.first].interior,b=ports.regions[link.second].interior;
            const double side_a=vertical_wall?a.x-.5:a.y-.5,side_b=vertical_wall?b.x-.5:b.y-.5;
            ports_valid=ports_valid && side_a*side_b>0 && link.vertical!=vertical_wall && std::abs(link.area-.45)<1e-14;
        }
        check(ports_valid,"two openings on one grid face attach to their own air regions");
    }
    split.prepare(grid(2,0,1),grid(3,0,1),32);
    check(!split.build_regions(unit,{},1),"changing the grid invalidates old component scratch dimensions");
    split.prepare_regions(2048,4096);
    check(split.build_regions(unit,{},1) && split.regions.size()==6,"component preparation can be repeated for a changed grid");
    split.build(unit,{},1);
    check(split.regions.empty() && split.region_links.empty(),"ordinary geometry build does not expose stale component outputs");
    return failures?1:0;
}
