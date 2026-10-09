#include "ducklips/tongue_surface.hpp"
#include "ducklips/physics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace ducklips {
namespace {
constexpr double radius=MidsagittalTongue::body_radius;
constexpr double patch_half_width=.15; // same contact core as TractModel
constexpr double missing=std::numeric_limits<double>::infinity();
double cross(Point2 a,Point2 b) noexcept { return a.x*b.y-a.y*b.x; }
double roof(double x) noexcept {
    if (x<2) { return 7+std::sqrt(std::max(2.9*2.9-(x-2)*(x-2),0.0)); }
    const double arc=7+pi+x-2;
    return 9+MidsagittalTongue::outer_offset(arc);
}
double roof_slope(double x) noexcept {
    // One-sided at the native wall knots; the profile is continuous there.
    constexpr double h=1e-6;
    return (roof(x+h)-roof(x-h))/(2*h);
}
}
bool ConservativeTongueSurface::configure(Point2 centre,const TongueGesture& gesture) noexcept {
    error_=nullptr;centre_=centre;
    if (!std::isfinite(centre.x) || !std::isfinite(centre.y) || !std::isfinite(gesture.activation)
        || !std::isfinite(gesture.place_from_lips) || !std::isfinite(gesture.gap)
        || !std::isfinite(gesture.width) || gesture.gap<0 || gesture.width<=0) {
        error_="invalid conservative body gesture";return false;
    }
    const double arc=MidsagittalTongue::length-gesture.place_from_lips;
    const Point2 normal=MidsagittalTongue::normal(arc),mid=MidsagittalTongue::midline(arc);
    if (std::abs(normal.x)>1e-12 || normal.y<.999999) { error_="conservative dorsum requires an oral-wall target";return false; }
    const double x0=mid.x,half=patch_half_width;
    const double left_room=x0-half-(centre.x-radius),right_room=centre.x+radius-(x0+half);
    if (left_room<=1e-5 || right_room<=1e-5) { error_="dorsal contact patch lies outside resting body";return false; }
    // Preserve the slot's width as a shoulder length, not a new contact width.
    const double shoulder=std::max(gesture.width,1e-4);
    const double left=std::min(shoulder,left_room*.95),right=std::min(shoulder,right_room*.95);
    const double start=x0-half-left,end=x0+half+right;
    std::array<double,ConservativeTongueMap::maximum_knots> positions{};
    std::size_t count=257;
    for (std::size_t i=0;i<count;++i) { positions[i]=start+(end-start)*static_cast<double>(i)/256; }
    for (double x:{x0-half,x0,x0+half,2.,4.5,6.,7.2,8.4,9.4}) {
        if (x>start && x<end) { positions[count++]=x; }
    }
    std::sort(positions.begin(),positions.begin()+static_cast<std::ptrdiff_t>(count));
    count=static_cast<std::size_t>(std::unique(positions.begin(),positions.begin()+static_cast<std::ptrdiff_t>(count))-positions.begin());
    std::array<TongueLiftKnot,ConservativeTongueMap::maximum_knots> knots{};
    for (std::size_t i=0;i<count;++i) {
        const double x=positions[i],dx=x-centre.x,radial=std::sqrt(std::max(radius*radius-dx*dx,0.0));
        const double span=x<x0?left:right;
        const double q=std::clamp((std::abs(x-x0)-half)/span,0.0,1.0);
        const double weight=1-q*q*q*(10+q*(-15+6*q));
        const double dw=-30*q*q*(1-q)*(1-q)*(x<x0?-1.:1.)/span;
        const double lift=roof(x)-gesture.gap-centre.y-radial;
        knots[i]={x,lift*weight,(roof_slope(x)+dx/radial)*weight+lift*dw};
    }
    knots[0].lift=knots[0].slope=knots[count-1].lift=knots[count-1].slope=0;
    const double far=std::max(std::abs(x0-half-centre.x),std::abs(x0+half-centre.x));
    const double free=centre.y+std::sqrt(radius*radius-far*far);
    const double anchor=centre.y-radius/2;
    if (!map_.configure(std::span<const TongueLiftKnot>(knots.data(),count),anchor,free,x0)) { error_="invalid conservative lift profile";return false; }
    const double reach=std::clamp(gesture.activation,0.0,1.0);command_=reach*reach*(3-2*reach);
    // Activation above contact cannot translate tissue through the palate.
    // Contact force is not modeled by this kinematic component.
    for (std::size_t i=0;i<contour_count;++i) {
        const double angle=2*pi*static_cast<double>(i)/contour_count;
        const TongueMaterialPoint point=material({centre.x+radius*std::cos(angle),centre.y+radius*std::sin(angle)});
        if (!point.valid) { error_="conservative contour mapping did not converge";return false; }
        contour_[i]=point.position;
    }
    return true;
}
TongueMaterialPoint ConservativeTongueSurface::material(Point2 point) const noexcept { return map_.map(point,command_); }
double ConservativeTongueSurface::wall_penetration(Point2 translation) const noexcept {
    double maximum=0;
    for (Point2 point:contour_) {
        point.x+=translation.x;point.y+=translation.y;
        maximum=std::max(maximum,std::max(-.9-point.x,point.y-roof(point.x)));
    }
    return maximum;
}
double ConservativeTongueSurface::admissible_pull(Point2 translation,Point2 pull) const noexcept {
    double fraction=1;
    for (Point2 point:contour_) {
        point.x+=translation.x;point.y+=translation.y;
        if (point.x<-.9-1e-10 || point.y>roof(point.x)+1e-10) { return 0; }
        std::array<double,10> events{};std::size_t count=0;
        if (std::abs(pull.x)>1e-14) {
            for (double knot:{-.9,2.,4.5,6.,7.2,8.4,9.4}) {
                const double time=(knot-point.x)/pull.x;
                if (time>0 && time<fraction) { events[count++]=time; }
            }
        }
        events[count++]=fraction;std::sort(events.begin(),events.begin()+static_cast<std::ptrdiff_t>(count));
        double low=0;
        for (std::size_t i=0;i<count;++i) {
            double high=events[i];const double x=point.x+high*pull.x,y=point.y+high*pull.y;
            if (x>=-.9 && y<=roof(x)+1e-10) { low=high;continue; }
            // On each linear roof segment, or the concave bend, admissible
            // translations form an interval. Visiting every wall knot
            // prevents a valid endpoint from tunneling through a corner.
            for (std::size_t iteration=0;iteration<40;++iteration) {
                const double mid=(low+high)/2,mx=point.x+mid*pull.x,my=point.y+mid*pull.y;
                if (mx>=-.9 && my<=roof(mx)+1e-10) { low=mid; } else { high=mid; }
            }
            fraction=std::min(fraction,low);break;
        }
    }
    return fraction;
}
double ConservativeTongueSurface::ray_gap(Point2 wall,Point2 inward,Point2 translation) const noexcept {
    wall.x-=translation.x;wall.y-=translation.y;
    const TongueMaterialPoint inverse=map_.map(wall,-command_);
    if (!inverse.valid) { return std::numeric_limits<double>::quiet_NaN(); }
    const double dx=inverse.position.x-centre_.x,dy=inverse.position.y-centre_.y;
    if (dx*dx+dy*dy<=radius*radius+1e-11) { return 0; }
    double best=missing;
    for (std::size_t i=0;i<contour_count;++i) {
        const Point2 a{contour_[i].x-wall.x,contour_[i].y-wall.y};
        const Point2 b{contour_[(i+1)%contour_count].x-wall.x,contour_[(i+1)%contour_count].y-wall.y};
        const double ga=cross(a,inward),gb=cross(b,inward);
        if ((ga>0 && gb>0) || (ga<0 && gb<0)) { continue; }
        double low=2*pi*static_cast<double>(i)/contour_count,high=2*pi*static_cast<double>(i+1)/contour_count;
        double g_low=ga,angle=(low+high)/2;bool converged=false;
        for (std::size_t iteration=0;iteration<40;++iteration) {
            const double cosine=std::cos(angle),sine=std::sin(angle);
            const TongueMaterialPoint p=material({centre_.x+radius*cosine,centre_.y+radius*sine});
            if (!p.valid) { return std::numeric_limits<double>::quiet_NaN(); }
            const Point2 offset{p.position.x-wall.x,p.position.y-wall.y};const double g=cross(offset,inward);
            if (std::abs(g)<1e-10) { const double t=offset.x*inward.x+offset.y*inward.y;if (t>=-1e-10) { best=std::min(best,std::max(t,0.0)); }converged=true;break; }
            if ((g_low<=0 && g>=0) || (g_low>=0 && g<=0)) { high=angle; } else { low=angle;g_low=g; }
            const Point2 tangent{radius*(-p.gradient[0]*sine+p.gradient[1]*cosine),radius*(-p.gradient[2]*sine+p.gradient[3]*cosine)};
            const double derivative=cross(tangent,inward);const double candidate=angle-g/derivative;
            angle=std::isfinite(candidate) && candidate>low && candidate<high?candidate:(low+high)/2;
        }
        if (!converged) { return std::numeric_limits<double>::quiet_NaN(); }
    }
    return best;
}
}
