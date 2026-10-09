#include "ducklips/glottal_relaxation.hpp"
#include "ducklips/contour_kinetics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace ducklips {
bool glottal_relaxation_bias(const std::array<double,49>& a,const std::array<double,49>& r,
    const std::array<bool,7>& active,double step,std::array<double,49>& bias) noexcept {
    if(!std::isfinite(step) || step<0){return false;}
    for(std::size_t i=0;i<49;++i){if(!std::isfinite(a[i]) || !std::isfinite(r[i])){return false;}}
    constexpr double epsilon=std::numeric_limits<double>::epsilon();
    for(std::size_t i=0;i<7;++i){for(std::size_t j=0;j<i;++j){
        // An off-diagonal can cross zero by cancellation of signed products.
        // Its rounding scale is set by the diagonal quadratic forms, not by
        // the vanishing result of that cancellation.
        const double as=std::max({std::abs(a[7*i+j]),std::abs(a[7*j+i]),std::sqrt(std::abs(a[7*i+i]))*std::sqrt(std::abs(a[7*j+j]))});
        const double rs=std::max({std::abs(r[7*i+j]),std::abs(r[7*j+i]),std::sqrt(std::abs(r[7*i+i]))*std::sqrt(std::abs(r[7*j+j]))});
        if(std::abs(a[7*i+j]-a[7*j+i])>64*epsilon*as || std::abs(r[7*i+j]-r[7*j+i])>64*epsilon*rs){return false;}
    }}
    std::array<double,49> l{},s{},v{},rl{},lv{},candidate{};
    for(std::size_t i=0;i<7;++i){if(!active[i]){continue;}for(std::size_t j=0;j<=i;++j){if(!active[j]){continue;}
        double value=a[7*i+j];for(std::size_t k=0;k<j;++k){value-=l[7*i+k]*l[7*j+k];}
        if(i==j){if(value<=0 || !std::isfinite(value)){return false;}l[7*i+j]=std::sqrt(value);}
        else{l[7*i+j]=value/l[7*j+j];}
    }}
    for(std::size_t i=0;i<7;++i){v[7*i+i]=1.;for(std::size_t j=0;j<7;++j){for(std::size_t k=0;k<7;++k){rl[7*i+j]+=.5*(r[7*i+k]+r[7*k+i])*l[7*k+j];}}}
    for(std::size_t i=0;i<7;++i){for(std::size_t j=0;j<=i;++j){
        double value=0;for(std::size_t k=0;k<7;++k){value+=l[7*k+i]*rl[7*k+j];}
        if(!std::isfinite(value)){return false;}s[7*i+j]=s[7*j+i]=value;
    }}
    // Cyclic Jacobi rotations retain the orthogonal resistance modes. Test
    // each off-diagonal against its own diagonal scale, not the fastest mode.
    bool converged=false;
    for(int sweep=0;sweep<64;++sweep){bool changed=false;
        for(std::size_t p=0;p<7;++p){for(std::size_t q=p+1;q<7;++q){
            const double off=s[7*p+q],pp=s[7*p+p],qq=s[7*q+q];
            if(off==0 || std::abs(off)<=8*epsilon*std::sqrt(std::abs(pp))*std::sqrt(std::abs(qq))){continue;}
            changed=true;const double tau=(qq-pp)/(2*off);
            const double t=std::copysign(1.,tau)/(std::abs(tau)+std::hypot(1.,tau));
            const double c=1/std::sqrt(1+t*t),sn=t*c;
            for(std::size_t k=0;k<7;++k){if(k==p || k==q){continue;}
                const double kp=s[7*k+p],kq=s[7*k+q];
                s[7*k+p]=s[7*p+k]=c*kp-sn*kq;s[7*k+q]=s[7*q+k]=sn*kp+c*kq;
            }
            s[7*p+p]=pp-t*off;s[7*q+q]=qq+t*off;s[7*p+q]=s[7*q+p]=0.;
            for(std::size_t k=0;k<7;++k){const double kp=v[7*k+p],kq=v[7*k+q];v[7*k+p]=c*kp-sn*kq;v[7*k+q]=sn*kp+c*kq;}
        }}
        if(!changed){converged=true;break;}
    }
    if(!converged){return false;}
    double scale=0;for(std::size_t i=0;i<7;++i){scale=std::max(scale,std::abs(s[7*i+i]));}
    std::array<double,7> weight{};
    for(std::size_t i=0;i<7;++i){const double rate=s[7*i+i];if(!std::isfinite(rate) || rate < -128*epsilon*scale){return false;}
        weight[i]=contour_relaxation_weight(step*std::max(0.,rate));}
    for(std::size_t i=0;i<7;++i){for(std::size_t j=0;j<7;++j){for(std::size_t k=0;k<7;++k){lv[7*i+j]+=l[7*i+k]*v[7*k+j];}}}
    for(std::size_t i=0;i<7;++i){for(std::size_t j=0;j<=i;++j){double value=0;
        for(std::size_t k=0;k<7;++k){value+=lv[7*i+k]*weight[k]*lv[7*j+k];}
        if(!std::isfinite(value)){return false;}candidate[7*i+j]=candidate[7*j+i]=value;
    }}
    bias=candidate;return true;
}
}
