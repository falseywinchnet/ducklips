#include "ducklips/glottal_contour.hpp"
#include <algorithm>
#include <cmath>
namespace ducklips {
namespace {
bool interval(const GlottalContour& p,double first,double last) noexcept {
    return p.valid() && std::isfinite(first) && std::isfinite(last) && first>=0 && last>first && last<=p.lower_thickness+p.upper_thickness;
}
std::size_t partition(const GlottalContour& p,double first,double last,std::array<double,5>& points) noexcept {
    std::size_t n=0;points[n++]=first;
    for(double z:{.5*p.lower_thickness,p.lower_thickness,p.lower_thickness+.5*p.upper_thickness}){if(z>first && z<last){points[n++]=z;}}
    points[n++]=last;return n;
}
double lower_weight(const GlottalContour& p,double z) noexcept {
    if(z<=.5*p.lower_thickness){return 1.;}
    if(z>=p.lower_thickness+.5*p.upper_thickness){return 0.;}
    return (p.lower_thickness+.5*p.upper_thickness-z)/(p.lower_thickness+.5*p.upper_thickness-.5*p.lower_thickness);
}
}
bool GlottalContour::valid() const noexcept {
    return std::isfinite(lower_thickness) && lower_thickness>0 && std::isfinite(upper_thickness) && upper_thickness>0
        && std::isfinite(lower_thickness+upper_thickness) && std::isfinite(lower_area) && std::isfinite(upper_area);
}
double GlottalContour::signed_area(double z) const noexcept {
    const double w=lower_weight(*this,z);return w*lower_area+(1-w)*upper_area;
}
bool GlottalContour::measure(double first,double last,GlottalContourMeasure& result) const noexcept {
    if(!interval(*this,first,last)){return false;}
    GlottalContourMeasure candidate;candidate.centroid=.5*(first+last);double moment=0;
    std::array<double,5> points{};const auto n=partition(*this,first,last,points);
    for(std::size_t k=1;k<n;++k){double a=points[k-1],b=points[k],aa=signed_area(a),ab=signed_area(b);
        if(!std::isfinite(aa) || !std::isfinite(ab)){return false;}if(aa<=0 && ab<=0){continue;}
        if(aa<=0){a+=(b-a)*(-aa)/(ab-aa);aa=0.;}
        else if(ab<=0){b=a+(b-a)*aa/(aa-ab);ab=0.;}
        const double length=b-a,sum=aa+ab,volume=.5*length*sum;
        if(!std::isfinite(volume) || volume<0){return false;}
        candidate.volume+=volume;moment+=volume*(a+length*(aa+2*ab)/(3*sum));
        // The moving contact endpoint contributes no boundary term because
        // its air area is zero. At a zero-area plateau this chooses the
        // closed-side derivative; the event has no unique two-sided derivative.
        const double da=.5*length*(lower_weight(*this,a)+lower_weight(*this,b));
        candidate.area_derivative[0]+=da;candidate.area_derivative[1]+=length-da;
    }
    if(candidate.volume>0){candidate.centroid=moment/candidate.volume;}
    if(!std::isfinite(candidate.volume) || !std::isfinite(candidate.centroid)){return false;}result=candidate;return true;
}
bool GlottalContour::passage(double first,double last,GlottalContourPassage& result) const noexcept {
    if(!interval(*this,first,last)){return false;}
    GlottalContourPassage candidate;
    // The profile is monotone. A zero or occupied endpoint seals the path,
    // including the logarithmic resistance singularity at exact contact.
    if(signed_area(first)<=0 || signed_area(last)<=0){result=candidate;return true;}
    candidate.open=true;std::array<double,5> points{};const auto n=partition(*this,first,last,points);
    for(std::size_t k=1;k<n;++k){const double a=points[k-1],b=points[k],aa=signed_area(a),ab=signed_area(b),length=b-a;
        if(aa<=0 || ab<=0 || !std::isfinite(aa) || !std::isfinite(ab)){return false;}
        const double relative=(ab-aa)/aa;
        const double logarithm=std::abs(relative)<.5?std::log1p(relative):std::log(ab)-std::log(aa);
        const double inverse=relative==0?length/aa:length*logarithm/(ab-aa);
        const double cubed=length*(aa+ab)/(2*aa*aa*ab*ab);
        if(!std::isfinite(inverse) || inverse<=0 || !std::isfinite(cubed) || cubed<=0){return false;}
        const std::size_t cell=a>=lower_thickness?1:0;
        candidate.inverse_area[cell]+=inverse;candidate.inverse_area_cubed[cell]+=cubed;
    }
    for(std::size_t i=0;i<2;++i){if(!std::isfinite(candidate.inverse_area[i]) || !std::isfinite(candidate.inverse_area_cubed[i])){return false;}}
    result=candidate;return true;
}
bool GlottalContour::geometry(GlottalContourGeometry& result) const noexcept {
    if(!valid()){return false;}GlottalContourGeometry candidate;const double total=lower_thickness+upper_thickness;
    if(!measure(0.,lower_thickness,candidate.chamber[0]) || !measure(lower_thickness,total,candidate.chamber[1])){return false;}
    candidate.aperture={std::max(0.,lower_area),std::max(0.,signed_area(lower_thickness)),std::max(0.,upper_area)};
    if(candidate.chamber[0].volume>0 && !passage(0.,candidate.chamber[0].centroid,candidate.passage[0])){return false;}
    if(candidate.chamber[0].volume>0 && candidate.chamber[1].volume>0
        && !passage(candidate.chamber[0].centroid,candidate.chamber[1].centroid,candidate.passage[1])){return false;}
    if(candidate.chamber[1].volume>0 && !passage(candidate.chamber[1].centroid,total,candidate.passage[2])){return false;}
    result=candidate;return true;
}
}
