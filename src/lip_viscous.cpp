#include "ducklips/lip_viscous.hpp"
#include <array>
#include <cmath>
namespace ducklips {
namespace {
constexpr std::array<double,2> x4{0.33998104358485626,0.86113631159405257};
constexpr std::array<double,2> w4{0.65214515486254654,0.34785484513745346};
constexpr std::array<double,4> x8{0.18343464249564984,0.52553240991632899,0.79666647741362673,0.96028985649753618};
constexpr std::array<double,4> w8{0.36268378337836182,0.31370664587788716,0.22238103445337445,0.10122853629037652};
constexpr std::array<double,8> x16{0.095012509837637441,0.28160355077925892,0.45801677765722737,0.61787624440264377,0.75540440835500311,0.86563120238783164,0.9445750230732326,0.98940093499164994};
constexpr std::array<double,8> w16{0.18945061045506859,0.18260341504492361,0.16915651939500262,0.14959598881657665,0.12462897125553395,0.095158511682492994,0.062253523938647609,0.027152459411754055};
constexpr std::array<double,16> x32{0.048307665687738317,0.14447196158279649,0.23928736225213706,0.33186860228212767,0.42135127613063533,0.50689990893222936,0.5877157572407623,0.66304426693021523,0.73218211874028971,0.79448379596794239,0.84936761373256997,0.89632115576605209,0.93490607593773967,0.96476225558750639,0.98561151154526827,0.99726386184948157};
constexpr std::array<double,16> w32{0.096540088514727576,0.095638720079274597,0.093844399080804386,0.09117387869576371,0.087652093004403686,0.083311924226946624,0.078193895787070117,0.072345794108848199,0.0658222227763616,0.058684093478535461,0.05099805926237596,0.042835898022226926,0.034273862913021633,0.025392065309262728,0.016274394730906486,0.0070186100094703618};
}
bool lip_viscous_segments(const LipProfile& profile,double first,double last,double width,double density,
    std::span<LipViscousSegment> segments,LipViscousMeasure& measure) noexcept {
    std::span<const double> points,weights;
    switch(segments.size()){
    case 4:points=x4;weights=w4;break;case 8:points=x8;weights=w8;break;
    case 16:points=x16;weights=w16;break;case 32:points=x32;weights=w32;break;default:return false;
    }
    if(!std::isfinite(width) || width<=0){return false;}
    double g;if(!profile.conductance(first,last,density,g) || g<=0){return false;}
    const double mid=first+(last-first)*.5,half=(last-first)*.5;
    std::array<LipViscousSegment,32> candidate{};double approximate=0;
    for(std::size_t i=0;i<points.size();++i){for(std::size_t side=0;side<2;++side){
        const auto j=2*i+side;const double x=mid+(side?1.:-1.)*half*points[i],area=profile.area(x);
        const double inertia=density*half*weights[i]/area,gap=area/(2*width);
        if(!std::isfinite(inertia) || inertia<=0 || !std::isfinite(gap) || gap<=0){return false;}
        candidate[j]={x,inertia,gap};approximate+=inertia;
    }}
    const double exact=1/g,normalization=exact/approximate;
    if(!std::isfinite(exact) || !std::isfinite(normalization) || normalization<=0){return false;}
    for(std::size_t i=0;i<segments.size();++i){candidate[i].inertance*=normalization;if(!std::isfinite(candidate[i].inertance) || candidate[i].inertance<=0){return false;}}
    for(std::size_t i=0;i<segments.size();++i){segments[i]=candidate[i];}
    measure={g,approximate,normalization-1};return true;
}
bool lip_outlet_viscous_segments(const LipProfile& profile,double last,double width,double density,
    std::span<LipViscousSegment> segments,LipViscousMeasure& measure) noexcept {
    if(!profile.valid() || profile.place!=0 || profile.core<=0 || profile.base_area<=0
        || !std::isfinite(last) || last<=0 || last>profile.length
        || !std::isfinite(width) || width<=0 || !std::isfinite(density) || density<=0){return false;}
    std::span<const double> points,weights;
    switch(segments.size()){
    case 5:points=x4;weights=w4;break;case 9:points=x8;weights=w8;break;
    case 17:points=x16;weights=w16;break;case 33:points=x32;weights=w32;break;default:return false;
    }
    double g;if(!profile.conductance(0,last,density,g)){return false;}
    const double minimum=profile.minimum_area();std::array<LipViscousSegment,33> candidate{};
    if(last<=profile.core){
        // A wide opening can place the outlet cell centroid inside the flat
        // core. Its path then has no shoulder: keep the existing site slots
        // but make those zero-length loss elements explicitly absent.
        const double inertia=minimum>0?density*last/minimum:0.;
        const double gap=minimum/(2*width);
        if(!std::isfinite(inertia) || !std::isfinite(gap)
            || (minimum>0 && (inertia<=0 || gap<=0))){return false;}
        candidate[0]={last*.5,inertia,gap};
        for(std::size_t i=1;i<segments.size();++i){candidate[i]={last,0.,0.};}
        for(std::size_t i=0;i<segments.size();++i){segments[i]=candidate[i];}
        measure={g,inertia,0.};return true;
    }
    candidate[0]={profile.core*.5,minimum>0?density*profile.core/minimum:0.,minimum/(2*width)};
    double approximate=candidate[0].inertance;
    const double half=(last-profile.core)*.5,middle=profile.core+half;
    for(std::size_t i=0;i<points.size();++i){for(std::size_t side=0;side<2;++side){
        const auto j=1+2*i+side;const double x=middle+(side?1.:-1.)*half*points[i],area=profile.area(x);
        candidate[j]={x,density*half*weights[i]/area,area/(2*width)};approximate+=candidate[j].inertance;
    }}
    const double normalization=g>0?1/(g*approximate):1.;
    if(!std::isfinite(approximate) || !std::isfinite(normalization) || normalization<=0){return false;}
    for(std::size_t i=0;i<segments.size();++i){auto& s=candidate[i];s.inertance*=normalization;
        if(!std::isfinite(s.inertance) || !std::isfinite(s.half_gap) || s.inertance<0 || s.half_gap<0
            || (i>0 && (s.inertance==0 || s.half_gap==0))){return false;}
    }
    for(std::size_t i=0;i<segments.size();++i){segments[i]=candidate[i];}
    measure={g,approximate,normalization-1};return true;
}
}
