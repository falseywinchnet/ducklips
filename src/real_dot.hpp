#pragma once
#include <cstddef>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace ducklips::detail {
// Explicit independent accumulators, no fast-math. Both modal transport and
// its application use contiguous real feature vectors. Roundoff may differ
// across architectures; physical accuracy is checked against the oracle.
inline double real_dot(const double* first,const double* second,std::size_t count) noexcept {
    std::size_t i=0;double result=0;
#if defined(__aarch64__)
    float64x2_t s0=vdupq_n_f64(0),s1=s0,s2=s0,s3=s0;
    for (;i+8<=count;i+=8) {
        s0=vfmaq_f64(s0,vld1q_f64(first+i),vld1q_f64(second+i));
        s1=vfmaq_f64(s1,vld1q_f64(first+i+2),vld1q_f64(second+i+2));
        s2=vfmaq_f64(s2,vld1q_f64(first+i+4),vld1q_f64(second+i+4));
        s3=vfmaq_f64(s3,vld1q_f64(first+i+6),vld1q_f64(second+i+6));
    }
    result=vaddvq_f64(vaddq_f64(vaddq_f64(s0,s1),vaddq_f64(s2,s3)));
#else
    double s0=0,s1=0,s2=0,s3=0;
    for (;i+4<=count;i+=4) {
        s0+=first[i]*second[i];s1+=first[i+1]*second[i+1];
        s2+=first[i+2]*second[i+2];s3+=first[i+3]*second[i+3];
    }
    result=(s0+s1)+(s2+s3);
#endif
    for (;i<count;++i) { result+=first[i]*second[i]; }
    return result;
}
}
