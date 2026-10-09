#include "ducklips/planar_airspace.hpp"
#include "ducklips/physics.hpp"
#include "ducklips/cavity_network.hpp"
#include <limits>
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
    // Domain ports retain the one-sided owning region, including a boundary
    // lying inside a grid cell. An internal grid face cannot represent this.
    ducklips::PlanarAirspace boundary;boundary.prepare(grid(1,-.1,1.1),grid(1,-.1,1.1),32);boundary.prepare_regions(2048,4096);
    const std::array<std::size_t,4> all_edges{0,1,2,3};
    const std::size_t port_before=allocations.load();
    const bool port_valid=boundary.build_regions(unit,{},.02,all_edges);
    const std::size_t port_after=allocations.load();
    bool square_ports=port_valid && boundary.boundary_ports.size()==4;
    const std::array<ducklips::Point2,4> normals{{{0,-1},{1,0},{0,1},{-1,0}}};
    for(const auto& port:boundary.boundary_ports){square_ports=square_ports && port.region==0 && std::abs(port.area-.02)<1e-16
        && port.outward_normal.x==normals[port.edge].x && port.outward_normal.y==normals[port.edge].y;}
    check(square_ports,"domain caps expose exact area and outward direction inside a grid cell");
    check(port_before==port_after,"boundary extraction allocates nothing after preparation");
    const auto horizontal_barrier=box(-.1,1.1,.45,.55);const std::array<ducklips::PlanarContour,1> cut_mouth{{horizontal_barrier}};
    const std::array<std::size_t,1> right_cap{1};
    bool split_ports=boundary.build_regions(unit,cut_mouth,1,right_cap) && boundary.boundary_ports.size()==2;
    for(const auto& port:boundary.boundary_ports){split_ports=split_ports && std::abs(port.area-.45)<1e-14
        && ((port.centre.y<.5)==(boundary.regions[port.region].interior.y<.5));}
    check(split_ports,"two disconnected chambers on one cap retain separate source attachments");
    const auto exterior_contact=box(1,1.2,.2,.8);const std::array<ducklips::PlanarContour,1> exterior{{exterior_contact}};
    check(boundary.build_regions(unit,exterior,1,right_cap) && boundary.regions.size()==1 && boundary.boundary_ports.size()==1
        && std::abs(boundary.boundary_ports[0].area-.4)<1e-14 && std::abs(boundary.boundary_ports[0].centre.y-.5)<1e-14,
        "coincident exterior tissue closes its part of a boundary aperture");
    auto reversed_unit=unit;std::reverse(reversed_unit.begin(),reversed_unit.end());
    check(boundary.build_regions(reversed_unit,exterior,1,right_cap) && boundary.boundary_ports.size()==1
        && std::abs(boundary.boundary_ports[0].area-.4)<1e-14 && boundary.boundary_ports[0].outward_normal.x==1,
        "boundary winding does not reverse physical flow direction");
    const std::array<ducklips::Point2,3> triangular_mouth{{{0,0},{1,0},{0,1}}};
    check(boundary.build_regions(triangular_mouth,{},1,right_cap) && boundary.boundary_ports.size()==1
        && std::abs(boundary.boundary_ports[0].area-std::sqrt(2.))<1e-14
        && std::abs(boundary.boundary_ports[0].outward_normal.x-1/std::sqrt(2.))<1e-14
        && std::abs(boundary.boundary_ports[0].centre.x-.5)<1e-14,
        "sloping cap has geometric length normal and area centroid");
    const auto seal_cap=box(1,1.1,0,1);const std::array<ducklips::PlanarContour,1> exterior_seal{{seal_cap}};
    check(boundary.build_regions(unit,exterior_seal,1,right_cap) && boundary.boundary_ports.empty(),"complete external contact seals a boundary exactly");
    const auto nearly_sealed=box(1,1.1,0,1-1e-12);const std::array<ducklips::PlanarContour,1> narrow_cap{{nearly_sealed}};
    check(boundary.build_regions(unit,narrow_cap,1,right_cap) && boundary.boundary_ports.size()==1 && boundary.boundary_ports[0].area>0
        && boundary.boundary_ports[0].area<1.01e-12,"small positive boundary opening survives without an aperture threshold");
    const auto moved_cap=box(0,1.05,0,1);const std::size_t moving_port_before=allocations.load();
    bool moving_port_valid=true;
    for(double fraction:{0.,.25,.5,.75,1.}){
        moving_port_valid=boundary.build_between(unit,{},moved_cap,{},fraction,1,right_cap) && moving_port_valid;
        moving_port_valid=moving_port_valid && boundary.boundary_ports.size()==1 && std::abs(boundary.boundary_ports[0].area-1)<1e-14
            && std::abs(boundary.boundary_ports[0].centre.x-(1+.05*fraction))<1e-14;
    }
    check(moving_port_valid && allocations.load()==moving_port_before,"moving cap identity and area persist without new allocations");
    const std::array<std::size_t,2> duplicate_cap{1,1};
    const std::array<std::size_t,1> invalid_cap{4};
    check(!boundary.build_regions(unit,{},1,duplicate_cap) && !boundary.build_regions(unit,{},1,invalid_cap),"duplicate and invalid boundary identities are rejected");
    boundary.build_regions(unit,{},1,right_cap);boundary.build(unit,{},1);
    check(boundary.boundary_ports.empty(),"ordinary geometry clears prior boundary attachments");
    // A sloping tissue boundary sweeps a triangular air wedge into the right
    // cell. Its volume is h^2/2 and its shared aperture is h, not two linear
    // interpolants. All contour vertices have fixed correspondence.
    ducklips::PlanarAirspace motion;const std::array<double,3> mx{-1,0,1};const std::array<double,2> my{0,1};
    motion.prepare(mx,my,32);motion.prepare_regions(2048,4096);
    const auto mouth=box(-1,1,0,1);
    const std::array<ducklips::Point2,5> first{{{-2,2},{2,-2},{3,-2},{3,3},{-2,3}}};
    const std::array<ducklips::Point2,5> last{{{-2,2.4},{2.4,-2},{3,-2},{3,3},{-2,3}}};
    const std::array<ducklips::PlanarContour,1> from{{first}},to{{last}};
    const std::size_t motion_before=allocations.load();bool motion_valid=true;double volume_error=0,aperture_error=0;
    for(double fraction:{0.,.125,.25,.5,.75,1.}){
        motion_valid=motion.build_between(mouth,from,mouth,to,fraction,1)&&motion_valid;
        const double height=.4*fraction;
        volume_error=std::max(volume_error,std::abs(motion.volume[1]-.5*height*height));
        aperture_error=std::max(aperture_error,std::abs(motion.vertical_face[1]-height));
    }
    const std::size_t motion_after=allocations.load();
    check(motion_valid && volume_error<1e-14 && aperture_error<1e-14,"moving contours preserve quadratic wedge volume and linear aperture",volume_error);
    check(motion_before==motion_after,"sampling moving contour geometry allocates nothing after preparation");
    motion.build_between(mouth,from,mouth,to,.5,1);
    check(std::abs(motion.volume[1]-.02)<1e-14 && std::abs(motion.volume[1]-.04)>.019,
        "midpoint geometry does not substitute interpolated endpoint volumes",motion.volume[1]);
    check(!motion.build_between(mouth,from,mouth,to,std::numeric_limits<double>::quiet_NaN(),1)
        && !motion.build_between(mouth,from,mouth,to,1.01,1),"invalid motion fractions fail explicitly");
    const std::array<ducklips::PlanarContour,1> wrong{{mouth}};
    check(!motion.build_between(mouth,from,mouth,wrong,.5,1),"motion rejects mismatched vertex correspondence");
    // Couple the freshly measured wedge to the actual native pressure solver.
    // Midpoint aperture is the exact interval average for this linear motion.
    const ducklips::AirProperties air{};const double bulk=air.density*air.sound_speed*air.sound_speed;
    bool coupled_valid=true;double pressure_error=0,mass_error=0;
    for(std::size_t steps:{1u,2u,4u,8u}){for(bool closing:{false,true}){
        ducklips::CavityNetwork pressure;const std::array<ducklips::CavityConnection,1> reservoir{{{0,ducklips::cavity_reservoir}}};pressure.prepare(1,reservoir);
        const double expected_pressure=closing?.16/(.04-.16/bulk):-.16/(.04+.16/bulk);
        std::array<double,1> v{closing?.08:0},c{v[0]/bulk},p{closing?expected_pressure:0},in{};
        coupled_valid=pressure.reset(v,c,p)&&coupled_valid;const double initial=v[0]+c[0]*p[0];double discharged=0;
        for(std::size_t step=1;step<=steps;++step){
            const double t0=double(step-1)/double(steps),t1=double(step)/double(steps),middle=.5*(t0+t1);
            coupled_valid=motion.build_between(mouth,from,mouth,to,closing?1-middle:middle,1)&&coupled_valid;
            const double conductance=.1*motion.vertical_face[1];
            const std::array<ducklips::CavityPassage,1> opening{{{motion.vertical_face[1],1/conductance,0,0,0}}};
            coupled_valid=motion.build_between(mouth,from,mouth,to,closing?1-t1:t1,1)&&coupled_valid;
            v[0]=motion.volume[1];c[0]=v[0]/bulk;
            const auto result=pressure.advance(v,c,opening,in,t1-t0,air);coupled_valid=result.converged&&coupled_valid;
            pressure_error=std::max(pressure_error,std::abs(pressure.pressures()[0]-expected_pressure));
            discharged+=(t1-t0)*pressure.flows()[0];mass_error=std::max(mass_error,std::abs(v[0]+c[0]*pressure.pressures()[0]+discharged-initial));
        }
    }}
    check(coupled_valid && pressure_error<1e-10,"native moving contour measurements give exact finite wedge pressure under refinement",pressure_error);
    check(mass_error<1e-14,"native moving wedge exchanges its stored air conservatively",mass_error);
    return failures?1:0;
}
