#include "ducklips/tongue_surface.hpp"
#include "ducklips/articulation.hpp"
#include "ducklips/physics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>
namespace {
std::size_t failures=0;
void check(bool valid,const char* name,double error=0) { std::printf("%s: %s (%.9e)\n",name,valid?"pass":"FAIL",error);if (!valid) { ++failures; } }
double cross(ducklips::Point2 a,ducklips::Point2 b) { return a.x*b.y-a.y*b.x; }
double polygon_gap(const std::vector<ducklips::Point2>& p,ducklips::Point2 wall,ducklips::Point2 direction) {
    bool inside=false;double best=std::numeric_limits<double>::infinity();
    for (std::size_t i=0;i<p.size();++i) {
        const ducklips::Point2 a=p[i],b=p[(i+1)%p.size()],edge{b.x-a.x,b.y-a.y},offset{a.x-wall.x,a.y-wall.y};
        if ((a.y>wall.y)!=(b.y>wall.y)) { if (wall.x<a.x+(wall.y-a.y)*(b.x-a.x)/(b.y-a.y)) { inside=!inside; } }
        const double denominator=cross(direction,edge);
        if (std::abs(denominator)<1e-14) { continue; }
        const double t=cross(offset,edge)/denominator,u=cross(offset,direction)/denominator;
        if (t>=0 && u>=0 && u<=1) { best=std::min(best,t); }
    }
    return inside?0:best;
}
}
int main() {
    const ducklips::Point2 centre{3.059093051,5.322166367};
    const ducklips::TonguePosture posture{.body_x=centre.x,.body_y=centre.y,.jaw=0};
    const ducklips::TongueShape rest=ducklips::MidsagittalTongue::pose(posture,{},{});
    ducklips::ConservativeTongueSurface surface{};
    ducklips::TongueGesture body{.activation=1,.place_from_lips=6,.gap=0,.width=.9};
    check(surface.configure(centre,body),"dorsal surface configured");
    const ducklips::TongueShape posed=ducklips::MidsagittalTongue::pose(posture,body,{},&surface);
    const ducklips::TongueMaterialPoint root=surface.material(rest.blade_root);
    check(posed.valid && root.valid && std::hypot(posed.blade_root.x-root.position.x,posed.blade_root.y-root.position.y)<1e-12,"blade root is the mapped body material point");
    const ducklips::TongueGesture tip{.activation=1,.place_from_lips=2.2,.gap=.1};
    const ducklips::TongueShape attached=ducklips::MidsagittalTongue::pose(posture,body,tip,&surface);
    const ducklips::TongueMaterialPoint transported_root=surface.material(rest.blade_root);
    const double attachment=std::hypot(attached.blade_root.x-transported_root.position.x,attached.blade_root.y-transported_root.position.y);
    check(attached.valid && attachment<1e-12,"blade follows locally transported body material",attachment);
    const double length=std::hypot(attached.tip.x-attached.blade_root.x,attached.tip.y-attached.blade_root.y);
    check(length<=ducklips::MidsagittalTongue::blade_length+1e-12,"combined gesture preserves blade reach constraint",length);
    check(attached.body_wall_penetration<2e-10,"blade pull respects body wall contact",attached.body_wall_penetration);
    check(attached.tip_target_error<1e-10 && attached.attachment_motion>0 && attached.drag==0,
          "anterior deformation reaches combined target without rigid body drag",attached.tip_target_error);
    const ducklips::Point2 anchor{centre.x,centre.y-1.35};
    const ducklips::TongueMaterialPoint fixed=surface.material(anchor);
    check(fixed.valid && fixed.position.x==anchor.x && fixed.position.y==anchor.y,"attachment transport leaves lower anchor fixed");
    std::array<double,351> arcs{},old_gap{},new_gap{};
    for (std::size_t i=0;i<arcs.size();++i) { arcs[i]=17.5*static_cast<double>(i)/(arcs.size()-1); }
    ducklips::MidsagittalTongue::gaps(rest,posture,arcs,old_gap);
    body.activation=1e-7;check(surface.configure(centre,body),"vanishing gesture configured");
    const ducklips::TongueShape tiny=ducklips::MidsagittalTongue::pose(posture,body,{},&surface);
    ducklips::MidsagittalTongue::gaps(tiny,posture,arcs,new_gap);
    double onset=0;for (std::size_t i=0;i<arcs.size();++i) { onset=std::max(onset,std::abs(new_gap[i]-old_gap[i])); }
    check(onset<1e-8,"airway approaches the original circular body continuously",onset);
    body.activation=0;check(surface.configure(centre,body),"idle body accepts attachment deformation");
    const ducklips::TongueShape idle_tip=ducklips::MidsagittalTongue::pose(posture,body,tip,&surface);
    ducklips::MidsagittalTongue::gaps(idle_tip,posture,arcs,old_gap);
    body.activation=1e-7;check(surface.configure(centre,body),"near-idle body accepts same attachment deformation");
    const ducklips::TongueShape tiny_tip=ducklips::MidsagittalTongue::pose(posture,body,tip,&surface);
    ducklips::MidsagittalTongue::gaps(tiny_tip,posture,arcs,new_gap);
    onset=0;for (std::size_t i=0;i<arcs.size();++i) { onset=std::max(onset,std::abs(new_gap[i]-old_gap[i])); }
    check(idle_tip.valid && tiny_tip.valid && onset<1e-8,"body activation does not switch attachment mechanics at zero",onset);
    double worst=0,closure=0;bool hits_match=true,finite=true;
    std::vector<ducklips::Point2> contour(8192);
    for (double activation:{0.0,.5,1.0}) { for (bool move_tip:{false,true}) {
        body.activation=activation;check(surface.configure(centre,body),"surface sweep configured");
        if (move_tip) {
            const ducklips::TongueShape combined=ducklips::MidsagittalTongue::pose(posture,body,tip,&surface);
            check(combined.valid && combined.tip_target_error<1e-10,"combined surface sweep reaches tip target");
        }
        for (std::size_t i=0;i<contour.size();++i) {
            const double angle=2*ducklips::pi*static_cast<double>(i)/contour.size();
            const ducklips::TongueMaterialPoint p=surface.material({centre.x+2.7*std::cos(angle),centre.y+2.7*std::sin(angle)});finite=finite&&p.valid;contour[i]=p.position;
        }
        for (double arc:arcs) {
            const ducklips::Point2 n=ducklips::MidsagittalTongue::normal(arc),m=ducklips::MidsagittalTongue::midline(arc);
            const double offset=ducklips::MidsagittalTongue::outer_offset(arc);
            const ducklips::Point2 wall{m.x+n.x*offset,m.y+n.y*offset},direction{-n.x,-n.y};
            const double exact=surface.ray_gap(wall,direction),reference=polygon_gap(contour,wall,direction);
            finite=finite&&!std::isnan(exact);hits_match=hits_match&&(std::isfinite(exact)==std::isfinite(reference));
            if (std::isfinite(exact) && std::isfinite(reference)) { worst=std::max(worst,std::abs(exact-reference)); }
            if (activation==1 && std::abs(arc-11.5)<.14) { closure=std::max(closure,exact); }
        }
    } }
    check(finite && hits_match,"mapped intersections agree with independent dense contour hit topology");
    check(worst<5e-5,"mapped intersections agree with dense contour distances",worst);
    check(closure<1e-8,"native finite contact core closes at airway rays",closure);
    body.place_from_lips=10;check(!surface.configure(centre,body) && surface.error()!=nullptr,"unsupported wall target explicitly rejected");
    ducklips::TractModel model{};model.prepare(ducklips::baritone_voice(),64,24);model.set_conservative_body(true);
    ducklips::TractShape shape{};model.allocate(shape);ducklips::ArticulatorState state{};state.place_vowel({});
    ducklips::ConstrictionSlotState& slot=state.slot[static_cast<std::size_t>(ducklips::ConstrictionSlot::tongue_body)];
    slot.activation.place(1);slot.place.place(.10);slot.width.place(.009);slot.target_area.place(0);model.build(state,shape);
    check(model.geometry_error()!=nullptr,"tract model reports unsupported geometry without rigid fallback");
    return failures?1:0;
}
