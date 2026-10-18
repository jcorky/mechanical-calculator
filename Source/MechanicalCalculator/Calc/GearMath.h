// Elementary and special functions for the gear-driven calculator, computed the way an operator
// at a stepped-drum machine would: only the machine's five procedures (add, subtract, multiply,
// divide, square root) applied in sequence, by series and argument reduction. The constants pi,
// ln 2 and ln 10 are themselves computed on the gears the first time they are needed (Machin's
// formula and atanh series), not typed in.
//
// The names mirror <cmath> so library code reads the same: gm::sqrt, gm::exp, gm::atan2, ...
#pragma once

#include "GearReal.h"

namespace mcalc
{
namespace gm
{
const Real& Pi();
const Real& Ln2();
const Real& Ln10();
inline Real Inf() { return Real::Infinity(false); }
inline Real NaN() { return Real::NaN(); }

inline bool isnan(const Real& X) { return X.IsNaN(); }
inline bool isinf(const Real& X) { return X.IsInf(); }
inline bool isfinite(const Real& X) { return X.IsFinite(); }
inline bool signbit(const Real& X) { return X.IsNegative(); }
inline Real abs(const Real& X) { return Real::Abs(X); }
inline Real fabs(const Real& X) { return Real::Abs(X); }
inline Real floor(const Real& X) { return Real::Floor(X); }
inline Real ceil(const Real& X) { return Real::Ceil(X); }
inline Real trunc(const Real& X) { return Real::Trunc(X); }
inline Real round(const Real& X) { return Real::Round(X); }
inline long long lround(const Real& X) { return Real::Round(X).ToLong(); }
inline long long llround(const Real& X) { return Real::Round(X).ToLong(); }
inline Real sqrt(const Real& X) { return Real::Sqrt(X); }
inline Real copysign(const Real& X, const Real& S) { return S.IsNegative() ? -Real::Abs(X) : Real::Abs(X); }
inline Real min(const Real& A, const Real& B) { return B < A ? B : A; }
inline Real max(const Real& A, const Real& B) { return A < B ? B : A; }
inline Real fmin(const Real& A, const Real& B) { return A.IsNaN() ? B : (B.IsNaN() ? A : min(A, B)); }
inline Real fmax(const Real& A, const Real& B) { return A.IsNaN() ? B : (B.IsNaN() ? A : max(A, B)); }

Real fmod(const Real& A, const Real& B);   // truncated remainder, sign of A
Real hypot(const Real& X, const Real& Y);
Real hypot(const Real& X, const Real& Y, const Real& Z);
Real cbrt(const Real& X);
Real exp(const Real& X);
Real expm1(const Real& X);
Real log(const Real& X);                   // natural
Real log1p(const Real& X);
Real log10(const Real& X);
Real log2(const Real& X);
Real pow(const Real& X, const Real& Y);
Real powi(const Real& X, long long N);     // repeated multiplication
Real sin(const Real& X);
Real cos(const Real& X);
Real tan(const Real& X);
Real asin(const Real& X);
Real acos(const Real& X);
Real atan(const Real& X);
Real atan2(const Real& Y, const Real& X);
Real sinh(const Real& X);
Real cosh(const Real& X);
Real tanh(const Real& X);
Real asinh(const Real& X);
Real acosh(const Real& X);
Real atanh(const Real& X);
Real erf(const Real& X);
Real erfc(const Real& X);
Real lgamma(const Real& X);                // log |Gamma(x)|
Real tgamma(const Real& X);
Real ldexp(const Real& X, int E);          // X * 2^E
} // namespace gm
} // namespace mcalc
