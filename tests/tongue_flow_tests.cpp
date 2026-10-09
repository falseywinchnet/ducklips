#include "ducklips/tongue_flow.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace {
std::size_t failures=0;
void check(bool valid,const char* name,double error=0) {
    std::printf("%s: %s (%.9e)\n",name,valid?"pass":"FAIL",error);if (!valid) { ++failures; }
}
}
int main() {
    ducklips::ConservativeTongueMap flow{};
    const std::array<ducklips::TongueLiftKnot,4> knots{{{-2,0,0},{-1,2,0},{1,2,0},{2,0,0}}};
    check(flow.configure(knots,-1,1,0),"finite compact lift accepted");
    double determinant_error=0,reversal_error=0,gradient_error=0,anchor_error=0,shear_error=0;bool valid=true;
    for (std::size_t ix=0;ix<17;++ix) {
        for (std::size_t iy=0;iy<17;++iy) {
            const ducklips::Point2 p{-2.1+.27*static_cast<double>(ix),-1.2+.21*static_cast<double>(iy)};
            const ducklips::TongueMaterialPoint q=flow.map(p,1,32),back=flow.map(q.position,-1,32);valid=valid&&q.valid&&back.valid;
            const std::array<double,4>& f=q.gradient;
            determinant_error=std::max(determinant_error,std::abs(f[0]*f[3]-f[1]*f[2]-1));
            reversal_error=std::max(reversal_error,std::hypot(back.position.x-p.x,back.position.y-p.y));
            if (p.y<=-1) { anchor_error=std::max(anchor_error,std::hypot(q.position.x-p.x,q.position.y-p.y)); }
            if (p.y>=1) { shear_error=std::max(shear_error,std::hypot(q.position.x-p.x,q.position.y-p.y-flow.lift(p.x))); }
            const double h=1e-5;
            for (std::size_t coordinate=0;coordinate<2;++coordinate) {
                ducklips::Point2 plus=p,minus=p;
                if (coordinate==0) { plus.x+=h;minus.x-=h; } else { plus.y+=h;minus.y-=h; }
                const ducklips::TongueMaterialPoint a=flow.map(plus,1,32),b=flow.map(minus,1,32);valid=valid&&a.valid&&b.valid;
                gradient_error=std::max(gradient_error,std::abs((a.position.x-b.position.x)/(2*h)-f[coordinate]));
                gradient_error=std::max(gradient_error,std::abs((a.position.y-b.position.y)/(2*h)-f[2+coordinate]));
            }
        }
    }
    check(valid,"all native forward, inverse and difference probes converge");
    check(determinant_error<1e-11,"local area conservation",determinant_error);
    check(reversal_error<1e-10,"forward/inverse material recovery",reversal_error);
    check(gradient_error<2e-6,"gradient agrees with independently differenced positions",gradient_error);
    check(anchor_error==0,"anchored tissue remains exactly fixed",anchor_error);
    check(shear_error<1e-12,"upper tissue follows exact prescribed shear",shear_error);
    check(!flow.map({0,0},1,0).valid,"zero substep count rejected");
    check(!flow.map({0,0},std::numeric_limits<double>::quiet_NaN()).valid,"nonfinite command rejected");
    check(!flow.configure(knots,1,-1,0),"reversed anchor interval rejected");
    check(!flow.map({0,0},1).valid,"failed configuration disables map");
    ducklips::TongueMaterialTranslation local{};
    const ducklips::Point2 start{.3,.2},shift{1.0,.9};
    check(local.configure(start,shift,2.2),"local attachment field configured");
    const ducklips::TongueMaterialPoint target=local.map(start);
    check(target.valid && std::hypot(target.position.x-start.x-shift.x,target.position.y-start.y-shift.y)<1e-12,
          "attachment follows its exact prescribed displacement");
    determinant_error=0;gradient_error=0;reversal_error=0;valid=true;
    for (int ix=-12;ix<=12;++ix) { for (int iy=-12;iy<=12;++iy) {
        const ducklips::Point2 p{.25*ix,.25*iy};
        const ducklips::TongueMaterialPoint q=local.map(p),back=local.map(q.position,true);
        valid=valid&&q.valid&&back.valid;const std::array<double,4>& f=q.gradient;
        determinant_error=std::max(determinant_error,std::abs(f[0]*f[3]-f[1]*f[2]-1));
        reversal_error=std::max(reversal_error,std::hypot(back.position.x-p.x,back.position.y-p.y));
        const double h=1e-5;
        for (std::size_t coordinate=0;coordinate<2;++coordinate) {
            ducklips::Point2 plus=p,minus=p;
            if (coordinate==0) { plus.x+=h;minus.x-=h; } else { plus.y+=h;minus.y-=h; }
            const ducklips::TongueMaterialPoint a=local.map(plus),b=local.map(minus);valid=valid&&a.valid&&b.valid;
            gradient_error=std::max(gradient_error,std::abs((a.position.x-b.position.x)/(2*h)-f[coordinate]));
            gradient_error=std::max(gradient_error,std::abs((a.position.y-b.position.y)/(2*h)-f[2+coordinate]));
        }
    } }
    check(valid && determinant_error<1e-11,"moving local field preserves material area",determinant_error);
    check(reversal_error<1e-10,"moving local field reverses in material coordinates",reversal_error);
    check(gradient_error<2e-6,"local field gradient agrees with independent differences",gradient_error);
    const ducklips::TongueMaterialPoint fixed=local.map({-4,-4});
    check(fixed.valid && fixed.position.x==-4 && fixed.position.y==-4,"material outside swept support is exactly fixed");
    check(!local.map(start,false,0).valid,"local field rejects zero steps");
    check(!local.configure(start,shift,0) && !local.map(start).valid,"failed local field configuration disables mapping");
    return failures?1:0;
}
