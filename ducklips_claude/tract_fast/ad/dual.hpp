#pragma once

// Forward-mode dual numbers carrying a fixed block of N tangents through one
// pass: value v and derivatives d[0..N-1] with respect to N parameter
// directions. The templated copies in this directory are written once over a
// scalar type S; S = double is the primal (to be compared bit for bit with
// the reference engine), S = Dual<N> carries N Jacobian columns.
//
// Branches in the copied code test values only (value_of), so a dual takes
// the branch its value takes: the derivative is that of the piece the
// trajectory is on (undefined only on the measure-zero switching set).

#include <algorithm>
#include <cmath>

namespace tract_ad {

template <int N>
struct Dual {
  double v;
  double d[N];
};

template <int N>
inline Dual<N> make_constant(double value) {
  Dual<N> out;
  out.v = value;
  for (int k = 0; k < N; ++k) {
    out.d[k] = 0.0;
  }
  return out;
}

inline double value_of(double x) { return x; }

template <int N>
inline double value_of(const Dual<N>& x) {
  return x.v;
}

// A scalar type's constant (double stays double).
template <class S>
struct Constant;

template <>
struct Constant<double> {
  static double make(double value) { return value; }
};

template <int N>
struct Constant<Dual<N>> {
  static Dual<N> make(double value) { return make_constant<N>(value); }
};

template <int N>
inline Dual<N> operator+(const Dual<N>& a, const Dual<N>& b) {
  Dual<N> out;
  out.v = a.v + b.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = a.d[k] + b.d[k];
  }
  return out;
}

template <int N>
inline Dual<N> operator+(const Dual<N>& a, double b) {
  Dual<N> out = a;
  out.v = a.v + b;
  return out;
}

template <int N>
inline Dual<N> operator+(double a, const Dual<N>& b) {
  Dual<N> out = b;
  out.v = a + b.v;
  return out;
}

template <int N>
inline Dual<N> operator-(const Dual<N>& a, const Dual<N>& b) {
  Dual<N> out;
  out.v = a.v - b.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = a.d[k] - b.d[k];
  }
  return out;
}

template <int N>
inline Dual<N> operator-(const Dual<N>& a, double b) {
  Dual<N> out = a;
  out.v = a.v - b;
  return out;
}

template <int N>
inline Dual<N> operator-(double a, const Dual<N>& b) {
  Dual<N> out;
  out.v = a - b.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = -b.d[k];
  }
  return out;
}

template <int N>
inline Dual<N> operator-(const Dual<N>& a) {
  Dual<N> out;
  out.v = -a.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = -a.d[k];
  }
  return out;
}

template <int N>
inline Dual<N> operator*(const Dual<N>& a, const Dual<N>& b) {
  Dual<N> out;
  out.v = a.v * b.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = a.d[k] * b.v + a.v * b.d[k];
  }
  return out;
}

template <int N>
inline Dual<N> operator*(const Dual<N>& a, double b) {
  Dual<N> out;
  out.v = a.v * b;
  for (int k = 0; k < N; ++k) {
    out.d[k] = a.d[k] * b;
  }
  return out;
}

template <int N>
inline Dual<N> operator*(double a, const Dual<N>& b) {
  Dual<N> out;
  out.v = a * b.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = a * b.d[k];
  }
  return out;
}

template <int N>
inline Dual<N> operator/(const Dual<N>& a, const Dual<N>& b) {
  Dual<N> out;
  out.v = a.v / b.v;
  double inverse = 1.0 / b.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = (a.d[k] - out.v * b.d[k]) * inverse;
  }
  return out;
}

template <int N>
inline Dual<N> operator/(const Dual<N>& a, double b) {
  Dual<N> out;
  out.v = a.v / b;
  double inverse = 1.0 / b;
  for (int k = 0; k < N; ++k) {
    out.d[k] = a.d[k] * inverse;
  }
  return out;
}

template <int N>
inline Dual<N> operator/(double a, const Dual<N>& b) {
  Dual<N> out;
  out.v = a / b.v;
  double slope = -out.v / b.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = slope * b.d[k];
  }
  return out;
}

template <int N>
inline Dual<N>& operator+=(Dual<N>& a, const Dual<N>& b) {
  a = a + b;
  return a;
}

template <int N>
inline Dual<N>& operator-=(Dual<N>& a, const Dual<N>& b) {
  a = a - b;
  return a;
}

// Elementary functions: value by the library function, derivative by its
// closed form at that value.
inline double ad_sqrt(double x) { return std::sqrt(x); }

template <int N>
inline Dual<N> ad_sqrt(const Dual<N>& x) {
  Dual<N> out;
  out.v = std::sqrt(x.v);
  double slope = 0.5 / out.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = slope * x.d[k];
  }
  return out;
}

inline double ad_abs(double x) { return std::abs(x); }

template <int N>
inline Dual<N> ad_abs(const Dual<N>& x) {
  if (x.v >= 0.0) {
    Dual<N> out = x;
    out.v = std::abs(x.v);
    return out;
  }
  return -x;
}

inline double ad_exp(double x) { return std::exp(x); }

template <int N>
inline Dual<N> ad_exp(const Dual<N>& x) {
  Dual<N> out;
  out.v = std::exp(x.v);
  for (int k = 0; k < N; ++k) {
    out.d[k] = out.v * x.d[k];
  }
  return out;
}

inline double ad_pow(double x, double exponent) { return std::pow(x, exponent); }

template <int N>
inline Dual<N> ad_pow(const Dual<N>& x, double exponent) {
  Dual<N> out;
  out.v = std::pow(x.v, exponent);
  double slope = exponent * out.v / x.v;
  for (int k = 0; k < N; ++k) {
    out.d[k] = slope * x.d[k];
  }
  return out;
}

// max/min choose the argument by value, as std::max/std::min do (the first
// on a tie for max (std::max returns a unless a < b), the first for min).
inline double ad_max(double a, double b) { return std::max(a, b); }
inline double ad_min(double a, double b) { return std::min(a, b); }

template <int N>
inline Dual<N> ad_max(const Dual<N>& a, const Dual<N>& b) {
  return (a.v < b.v) ? b : a;
}

template <int N>
inline Dual<N> ad_max(const Dual<N>& a, double b) {
  return (a.v < b) ? make_constant<N>(b) : a;
}

template <int N>
inline Dual<N> ad_min(const Dual<N>& a, const Dual<N>& b) {
  return (b.v < a.v) ? b : a;
}

template <int N>
inline Dual<N> ad_min(const Dual<N>& a, double b) {
  return (b < a.v) ? make_constant<N>(b) : a;
}

// std::clamp(x, low, high): low if x < low, else high if high < x, else x.
inline double ad_clamp(double x, double low, double high) { return std::clamp(x, low, high); }

template <int N>
inline Dual<N> ad_clamp(const Dual<N>& x, double low, double high) {
  if (x.v < low) {
    return make_constant<N>(low);
  }
  if (high < x.v) {
    return make_constant<N>(high);
  }
  return x;
}

}  // namespace tract_ad
