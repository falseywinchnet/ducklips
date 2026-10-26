// Two-lane double arithmetic for the engine's and the native peel's vector
// paths. On ARM this is <arm_neon.h> itself, so the Mac mini's builds are
// unchanged byte for byte. Elsewhere (the x86-64 compute pod, 2026-10-26) the
// same names are defined on SSE2, each lane the IEEE operation NEON performs:
// add, subtract, multiply, divide and square root correctly rounded, the
// comparisons as all-ones or all-zero lane masks, the selects bitwise, and
// vfmaq_f64 fused (FMA3 where the compiler targets it, std::fma otherwise),
// so a lane's result equals the scalar code's on either machine.
#pragma once

#if defined(__aarch64__) || defined(__ARM_NEON)

#include <arm_neon.h>

#else

#include <cmath>
#include <cstdint>
#include <immintrin.h>

typedef __m128d float64x2_t;
typedef __m128i uint64x2_t;

inline float64x2_t vld1q_f64(const double* source) { return _mm_loadu_pd(source); }
inline void vst1q_f64(double* destination, float64x2_t value) { _mm_storeu_pd(destination, value); }
inline float64x2_t vdupq_n_f64(double value) { return _mm_set1_pd(value); }

inline float64x2_t vaddq_f64(float64x2_t left, float64x2_t right) { return _mm_add_pd(left, right); }
inline float64x2_t vsubq_f64(float64x2_t left, float64x2_t right) { return _mm_sub_pd(left, right); }
inline float64x2_t vmulq_f64(float64x2_t left, float64x2_t right) { return _mm_mul_pd(left, right); }
inline float64x2_t vdivq_f64(float64x2_t left, float64x2_t right) { return _mm_div_pd(left, right); }
inline float64x2_t vsqrtq_f64(float64x2_t value) { return _mm_sqrt_pd(value); }
inline float64x2_t vnegq_f64(float64x2_t value) { return _mm_xor_pd(value, _mm_set1_pd(-0.0)); }
inline float64x2_t vabsq_f64(float64x2_t value) { return _mm_andnot_pd(_mm_set1_pd(-0.0), value); }

// NEON's vfmaq_f64(a, b, c) is a + b c with one rounding.
inline float64x2_t vfmaq_f64(float64x2_t addend, float64x2_t left, float64x2_t right) {
#if defined(__FMA__)
  return _mm_fmadd_pd(left, right, addend);
#else
  alignas(16) double a[2];
  alignas(16) double b[2];
  alignas(16) double c[2];
  _mm_store_pd(a, addend);
  _mm_store_pd(b, left);
  _mm_store_pd(c, right);
  return _mm_set_pd(std::fma(b[1], c[1], a[1]), std::fma(b[0], c[0], a[0]));
#endif
}

inline uint64x2_t vceqq_f64(float64x2_t left, float64x2_t right) { return _mm_castpd_si128(_mm_cmpeq_pd(left, right)); }
inline uint64x2_t vcltq_f64(float64x2_t left, float64x2_t right) { return _mm_castpd_si128(_mm_cmplt_pd(left, right)); }
inline uint64x2_t vcleq_f64(float64x2_t left, float64x2_t right) { return _mm_castpd_si128(_mm_cmple_pd(left, right)); }
inline uint64x2_t vcgtq_f64(float64x2_t left, float64x2_t right) { return _mm_castpd_si128(_mm_cmpgt_pd(left, right)); }
inline uint64x2_t vcgeq_f64(float64x2_t left, float64x2_t right) { return _mm_castpd_si128(_mm_cmpge_pd(left, right)); }

// The selects take each bit from the first operand where the mask is set.
inline float64x2_t vbslq_f64(uint64x2_t mask, float64x2_t chosen, float64x2_t otherwise) {
  __m128d bits = _mm_castsi128_pd(mask);
  return _mm_or_pd(_mm_and_pd(bits, chosen), _mm_andnot_pd(bits, otherwise));
}
inline uint64x2_t vbslq_u64(uint64x2_t mask, uint64x2_t chosen, uint64x2_t otherwise) {
  return _mm_or_si128(_mm_and_si128(mask, chosen), _mm_andnot_si128(mask, otherwise));
}
inline uint64x2_t vreinterpretq_u64_f64(float64x2_t value) { return _mm_castpd_si128(value); }

// 32-bit lane views of a mask (the peel's complement of a select mask): the
// same 128 bits, so the casts are the identity and NOT flips every bit.
typedef __m128i uint32x4_t;
inline uint32x4_t vreinterpretq_u32_u64(uint64x2_t value) { return value; }
inline uint64x2_t vreinterpretq_u64_u32(uint32x4_t value) { return value; }
inline uint32x4_t vmvnq_u32(uint32x4_t value) { return _mm_xor_si128(value, _mm_set1_epi32(-1)); }

// Lane access with the lane a compile-time constant, as NEON requires.
#define vgetq_lane_f64(value, lane) (tract_lane_f64<lane>(value))
#define vsetq_lane_f64(scalar, value, lane) (tract_set_lane_f64<lane>((scalar), (value)))
#define vgetq_lane_u64(value, lane) (tract_lane_u64<lane>(value))

template <int Lane>
inline double tract_lane_f64(float64x2_t value) {
  alignas(16) double out[2];
  _mm_store_pd(out, value);
  return out[Lane];
}
template <int Lane>
inline float64x2_t tract_set_lane_f64(double scalar, float64x2_t value) {
  alignas(16) double out[2];
  _mm_store_pd(out, value);
  out[Lane] = scalar;
  return _mm_load_pd(out);
}
template <int Lane>
inline std::uint64_t tract_lane_u64(uint64x2_t value) {
  alignas(16) std::uint64_t out[2];
  _mm_store_si128(reinterpret_cast<__m128i*>(out), value);
  return out[Lane];
}

#endif
