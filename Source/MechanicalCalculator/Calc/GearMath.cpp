#include "GearMath.h"

#include <string>

namespace mcalc
{
namespace gm
{
namespace gearmath_detail
{
// Below this a series term no longer reaches the 18th digit of the sum.
const Real& Tiny()
{
	static const Real T = Real::FromParts(false, 100000000000000000ull, -21);   // 1e-21
	return T;
}

bool Negligible(const Real& Term, const Real& Sum)
{
	return Real::Abs(Term) <= Real::Abs(Sum) * Tiny() || Term.IsZero();
}

// A long fixed-point register for the constants: 9 integer digits and 63 decimals in base-1e9
// limbs. Terms of the series are formed by the machine's exact division by small whole numbers,
// so the constants come out correct to far more than 18 digits; Hi + Lo splits them for use.
struct Fixed
{
	static constexpr int N = 8;          // D[N-1] is the integer limb
	uint64_t D[N] = {};

	static Fixed One()
	{
		Fixed F;
		F.D[N - 1] = 1;
		return F;
	}
	bool IsZero() const
	{
		for (uint64_t X : D)
		{
			if (X != 0)
			{
				return false;
			}
		}
		return true;
	}
	void DivSmall(uint64_t Div)
	{
		uint64_t Rem = 0;
		for (int i = N - 1; i >= 0; --i)
		{
			const uint64_t X = D[i] + Rem * 1000000000ull;
			D[i] = X / Div;
			Rem = X % Div;
		}
	}
	void MulSmall(uint64_t M)
	{
		uint64_t Carry = 0;
		for (int i = 0; i < N; ++i)
		{
			const uint64_t X = D[i] * M + Carry;
			D[i] = X % 1000000000ull;
			Carry = X / 1000000000ull;
		}
	}
	void Add(const Fixed& B)
	{
		uint64_t Carry = 0;
		for (int i = 0; i < N; ++i)
		{
			const uint64_t X = D[i] + B.D[i] + Carry;
			D[i] = X % 1000000000ull;
			Carry = X / 1000000000ull;
		}
	}
	void Sub(const Fixed& B)   // this >= B
	{
		uint64_t Borrow = 0;
		for (int i = 0; i < N; ++i)
		{
			const uint64_t Need = B.D[i] + Borrow;
			if (D[i] >= Need)
			{
				D[i] -= Need;
				Borrow = 0;
			}
			else
			{
				D[i] = D[i] + 1000000000ull - Need;
				Borrow = 1;
			}
		}
	}
	// Writes the value as Hi (first 18 significant digits, truncated) + Lo (the next 18).
	void Split(Real* Hi, Real* Lo, Real* Rounded) const
	{
		std::string Digits = std::to_string(D[N - 1]);
		const int IntDigits = D[N - 1] == 0 ? 0 : static_cast<int>(Digits.size());
		if (IntDigits == 0)
		{
			Digits.clear();
		}
		for (int i = N - 2; i >= 0; --i)
		{
			std::string Limb = std::to_string(D[i]);
			Digits += std::string(9 - Limb.size(), '0') + Limb;
		}
		// Digits holds IntDigits integer digits followed by 63 decimals.
		size_t First = Digits.find_first_not_of('0');
		const long long PointPos = IntDigits;   // decimal point after this many characters
		const std::string Sig = Digits.substr(First);
		const long long TopExp = PointPos - static_cast<long long>(First) - 1;   // power of the first digit
		auto Make = [](const std::string& Part, long long Exp) {
			const std::string Text = Part + "e" + std::to_string(Exp - static_cast<long long>(Part.size()) + 1);
			size_t Used = 0;
			return Real::Parse(Text.c_str(), &Used);
		};
		*Hi = Make(Sig.substr(0, 18), TopExp);
		*Lo = Make(Sig.substr(18, 18), TopExp - 18);
		*Rounded = Make(Sig.substr(0, 19), TopExp);
	}
};

// Adds Sign * atan(1/N) (or atanh(1/N)) into F, the terms formed by exact long division.
void ArcSeriesInto(Fixed& F, uint64_t N, bool Hyperbolic, uint64_t Weight)
{
	Fixed Power = Fixed::One();
	Power.MulSmall(Weight);
	Power.DivSmall(N);
	for (uint64_t K = 0; !Power.IsZero() && K < 400; ++K)
	{
		Fixed Term = Power;
		Term.DivSmall(2 * K + 1);
		if (Hyperbolic || K % 2 == 0)
		{
			F.Add(Term);
		}
		else
		{
			F.Sub(Term);
		}
		Power.DivSmall(N * N);
	}
}

struct SplitConstant
{
	Real Hi, Lo, Value;
};

const SplitConstant& HalfPiSplit()
{
	static const SplitConstant C = [] {
		// Machin: pi/2 = 8 atan(1/5) - 2 atan(1/239)
		Fixed F;
		ArcSeriesInto(F, 5, false, 8);
		Fixed G;
		ArcSeriesInto(G, 239, false, 2);
		F.Sub(G);
		SplitConstant S;
		F.Split(&S.Hi, &S.Lo, &S.Value);
		return S;
	}();
	return C;
}

const SplitConstant& PiSplit()
{
	static const SplitConstant C = [] {
		Fixed F;
		ArcSeriesInto(F, 5, false, 16);
		Fixed G;
		ArcSeriesInto(G, 239, false, 4);
		F.Sub(G);
		SplitConstant S;
		F.Split(&S.Hi, &S.Lo, &S.Value);
		return S;
	}();
	return C;
}

const SplitConstant& Ln2Split()
{
	static const SplitConstant C = [] {
		Fixed F;
		ArcSeriesInto(F, 3, true, 2);   // ln 2 = 2 atanh(1/3)
		SplitConstant S;
		F.Split(&S.Hi, &S.Lo, &S.Value);
		return S;
	}();
	return C;
}

const SplitConstant& Ln10Split()
{
	static const SplitConstant C = [] {
		Fixed F;
		ArcSeriesInto(F, 3, true, 6);   // ln 10 = 6 atanh(1/3) + 2 atanh(1/9)
		ArcSeriesInto(F, 9, true, 2);
		SplitConstant S;
		F.Split(&S.Hi, &S.Lo, &S.Value);
		return S;
	}();
	return C;
}

// atanh(z) = sum z^(2k+1) / (2k+1), for |z| <= 1/3
Real AtanhSeries(const Real& Z)
{
	const Real Z2 = Z * Z;
	Real Power = Z;
	RealAccumulator Acc;
	Acc.Add(Z);
	for (long long K = 1; K < 400; ++K)
	{
		Power = Power * Z2;
		const Real Term = Power / Real(2 * K + 1);
		Acc.Add(Term);
		if (Negligible(Term, Z))
		{
			break;
		}
	}
	return Acc.Result();
}

// atan(y) = sum (-1)^k y^(2k+1)/(2k+1), for |y| <= 0.2
Real AtanSeries(const Real& Y)
{
	const Real Y2 = Y * Y;
	Real Power = Y;
	RealAccumulator Acc;
	Acc.Add(Y);
	for (long long K = 1; K < 400; ++K)
	{
		Power = Power * Y2;
		const Real Term = Power / Real(2 * K + 1);
		if (K % 2)
		{
			Acc.Sub(Term);
		}
		else
		{
			Acc.Add(Term);
		}
		if (Negligible(Term, Y))
		{
			break;
		}
	}
	return Acc.Result();
}

Real SinSeries(const Real& R)
{
	const Real R2 = R * R;
	Real Term = R;
	RealAccumulator Acc;
	Acc.Add(R);
	for (long long N = 1; N < 60; ++N)
	{
		Term = -Term * R2 / Real((2 * N) * (2 * N + 1));
		Acc.Add(Term);
		if (Negligible(Term, R))
		{
			break;
		}
	}
	return Acc.Result();
}

Real CosSeries(const Real& R)
{
	const Real R2 = R * R;
	Real Term = Real(1);
	RealAccumulator Acc;
	Acc.Add(Term);
	for (long long N = 1; N < 60; ++N)
	{
		Term = -Term * R2 / Real((2 * N - 1) * (2 * N));
		Acc.Add(Term);
		if (Negligible(Term, Real(1)))
		{
			break;
		}
	}
	return Acc.Result();
}

// Reduces X by multiples of pi/2: X = Q * pi/2 + R, |R| <= pi/4.
Real ReduceQuarter(const Real& X, long long* Q)
{
	const SplitConstant& H = HalfPiSplit();
	const Real K = Real::Round(X / H.Value);
	*Q = K.ToLong();
	// X - K pi/2 with both products kept whole in the accumulator (Cody and Waite's two-part pi).
	RealAccumulator Acc;
	Acc.Add(X);
	Acc.AddProduct(K, H.Hi, true);
	Acc.AddProduct(K, H.Lo, true);
	return Acc.Result();
}

// exp(R) for |R| <= 1.2: the Taylor series summed in the accumulator.
Real ExpSmall(const Real& R)
{
	Real Term = Real(1);
	RealAccumulator Acc;
	Acc.Add(Term);
	for (long long N = 1; N < 80; ++N)
	{
		Term = Term * R / Real(N);
		Acc.Add(Term);
		if (Negligible(Term, Real(1)))
		{
			break;
		}
	}
	return Acc.Result();
}

// x - x^3/3! ... for |x| < 0.5 (sinh without cancellation)
Real SinhSeries(const Real& X)
{
	const Real X2 = X * X;
	Real Term = X;
	RealAccumulator Acc;
	Acc.Add(X);
	for (long long N = 1; N < 60; ++N)
	{
		Term = Term * X2 / Real((2 * N) * (2 * N + 1));
		Acc.Add(Term);
		if (Negligible(Term, X))
		{
			break;
		}
	}
	return Acc.Result();
}

// Stirling's series for log Gamma(x), x >= 15.
Real LgammaStirling(const Real& X)
{
	static const Real HalfLog2Pi = log(Real(2) * Pi()) / Real(2);
	// B_2k / (2k (2k - 1)) for k = 1..8
	static const Real Coef[8] = {Real(1) / Real(12), Real(-1) / Real(360), Real(1) / Real(1260), Real(-1) / Real(1680), Real(1) / Real(1188),
		Real(-691) / Real(360360), Real(1) / Real(156), Real(-3617) / Real(122400)};
	RealAccumulator Acc;
	Acc.Add((X - Real(0.5)) * log(X));
	Acc.Sub(X);
	Acc.Add(HalfLog2Pi);
	const Real InvX2 = Real(1) / (X * X);
	Real Power = Real(1) / X;
	for (const Real& C : Coef)
	{
		Acc.Add(C * Power);
		Power = Power * InvX2;
	}
	return Acc.Result();
}
} // namespace gearmath_detail

using namespace gearmath_detail;

const Real& Pi()
{
	// Machin, 1706: pi = 16 atan(1/5) - 4 atan(1/239), on the gears.
	return PiSplit().Value;
}

const Real& Ln2()
{
	return Ln2Split().Value;
}

const Real& Ln10()
{
	return Ln10Split().Value;
}

Real fmod(const Real& A, const Real& B)
{
	if (A.IsNaN() || B.IsNaN() || A.IsInf() || B.IsZero())
	{
		return Real::NaN();
	}
	if (B.IsInf() || A.IsZero())
	{
		return A;
	}
	// The machine's remainder: subtract |B| shifted up place by place, as in division.
	const Real Bm = Real::Abs(B);
	Real R = Real::Abs(A);
	for (int K = R.Exponent() - Bm.Exponent(); K >= 0; --K)
	{
		const Real Step = Real::Scale10(Bm, K);
		while (R >= Step)
		{
			R = R - Step;
		}
	}
	return A.IsNegative() ? -R : R;
}

Real hypot(const Real& X, const Real& Y)
{
	if (X.IsInf() || Y.IsInf())
	{
		return Real::Infinity();
	}
	return sqrt(X * X + Y * Y);
}

Real hypot(const Real& X, const Real& Y, const Real& Z)
{
	if (X.IsInf() || Y.IsInf() || Z.IsInf())
	{
		return Real::Infinity();
	}
	return sqrt(X * X + Y * Y + Z * Z);
}

Real cbrt(const Real& X)
{
	if (!X.IsFinite() || X.IsZero())
	{
		return X;
	}
	const Real A = Real::Abs(X);
	Real Y = exp(log(A) / Real(3));
	Y = Y - (Y * Y * Y - A) / (Real(3) * Y * Y);   // one Newton step polishes the last digits
	return X.IsNegative() ? -Y : Y;
}

Real exp(const Real& X)
{
	if (X.IsNaN())
	{
		return X;
	}
	if (X.IsInf())
	{
		return X.IsNegative() ? Real() : X;
	}
	static const Real Limit(2.3e8);
	if (X > Limit)
	{
		return Real::Infinity();
	}
	if (X < -Limit)
	{
		return Real();
	}
	// x = k ln 10 + r: the power of ten is a carriage shift, exp(r) is a series.
	const SplitConstant& L = Ln10Split();
	const Real K = Real::Round(X / L.Value);
	RealAccumulator Acc;
	Acc.Add(X);
	Acc.AddProduct(K, L.Hi, true);
	Acc.AddProduct(K, L.Lo, true);
	return Real::Scale10(ExpSmall(Acc.Result()), K.ToInt());
}

Real expm1(const Real& X)
{
	if (Real::Abs(X) < Real(0.5))
	{
		Real Term = X;
		RealAccumulator Acc;
		Acc.Add(X);
		for (long long N = 2; N < 60; ++N)
		{
			Term = Term * X / Real(N);
			Acc.Add(Term);
			if (Negligible(Term, X))
			{
				break;
			}
		}
		return Acc.Result();
	}
	return exp(X) - Real(1);
}

Real log(const Real& X)
{
	if (X.IsNaN() || X.IsNegative())
	{
		return Real::NaN();
	}
	if (X.IsZero())
	{
		return Real::Infinity(true);
	}
	if (X.IsInf())
	{
		return X;
	}
	if (X >= Real(0.5) && X <= Real(2))
	{
		return Real(2) * AtanhSeries((X - Real(1)) / (X + Real(1)));
	}
	// x = m 10^e with 1 <= m < 10; then m = 2^j y with 1 <= y <= 1.5.
	const int E = X.Exponent();
	Real M = Real::Scale10(X, -E);
	int J = 0;
	while (M > Real(1.5))
	{
		M = M / Real(2);
		++J;
	}
	RealAccumulator Acc;
	Acc.AddProduct(Real(E), Ln10Split().Hi);
	Acc.AddProduct(Real(E), Ln10Split().Lo);
	Acc.AddProduct(Real(J), Ln2Split().Hi);
	Acc.AddProduct(Real(J), Ln2Split().Lo);
	Acc.AddProduct(Real(2), AtanhSeries((M - Real(1)) / (M + Real(1))));
	return Acc.Result();
}

Real log1p(const Real& X)
{
	if (X.IsNaN() || X < Real(-1))
	{
		return Real::NaN();
	}
	if (X == Real(-1))
	{
		return Real::Infinity(true);
	}
	if (Real::Abs(X) <= Real(0.5))
	{
		return Real(2) * AtanhSeries(X / (Real(2) + X));
	}
	return log(Real(1) + X);
}

Real log10(const Real& X)
{
	if (X.IsFinite() && !X.IsZero() && !X.IsNegative())
	{
		const int E = X.Exponent();
		const Real M = Real::Scale10(X, -E);
		if (M == Real(1))
		{
			return Real(E);   // exact for powers of ten
		}
		return Real(E) + log(M) / Ln10();
	}
	return log(X);
}

Real log2(const Real& X)
{
	return log(X) / Ln2();
}

Real powi(const Real& X, long long N)
{
	if (N < 0)
	{
		return Real(1) / powi(X, -N);
	}
	Real Result(1);
	Real Base = X;
	while (N > 0)
	{
		if (N & 1)
		{
			Result = Result * Base;
		}
		N >>= 1;
		if (N > 0)
		{
			Base = Base * Base;
		}
	}
	return Result;
}

Real pow(const Real& X, const Real& Y)
{
	if (Y.IsZero() || X == Real(1))
	{
		return Real(1);
	}
	if (X.IsNaN() || Y.IsNaN())
	{
		return Real::NaN();
	}
	const bool IntegerY = Y.IsFinite() && Real::Trunc(Y) == Y && Real::Abs(Y) <= Real(1000000000);
	if (IntegerY && X.IsFinite())
	{
		if (X.IsZero())
		{
			return Y.IsNegative() ? Real::Infinity() : Real();
		}
		return powi(X, Y.ToLong());
	}
	if (X.IsZero())
	{
		return Y.IsNegative() ? Real::Infinity() : Real();
	}
	if (X.IsNegative())
	{
		if (X.IsInf() && IntegerY)
		{
			const bool Odd = Real::Abs(gm::fmod(Y, Real(2))) == Real(1);
			if (Y.IsNegative())
			{
				return Real();
			}
			return Real::Infinity(Odd);
		}
		return Real::NaN();
	}
	if (X.IsInf())
	{
		return Y.IsNegative() ? Real() : X;
	}
	if (Y.IsInf())
	{
		const bool Grows = (X > Real(1)) != Y.IsNegative();
		return Grows ? Real::Infinity() : Real();
	}
	if (Y == Real(0.5))
	{
		return sqrt(X);
	}
	return exp(Y * log(X));
}

Real sin(const Real& X)
{
	if (!X.IsFinite())
	{
		return Real::NaN();
	}
	long long Q = 0;
	const Real R = ReduceQuarter(X, &Q);
	switch (((Q % 4) + 4) % 4)
	{
	case 0: return SinSeries(R);
	case 1: return CosSeries(R);
	case 2: return -SinSeries(R);
	default: return -CosSeries(R);
	}
}

Real cos(const Real& X)
{
	if (!X.IsFinite())
	{
		return Real::NaN();
	}
	long long Q = 0;
	const Real R = ReduceQuarter(X, &Q);
	switch (((Q % 4) + 4) % 4)
	{
	case 0: return CosSeries(R);
	case 1: return -SinSeries(R);
	case 2: return -CosSeries(R);
	default: return SinSeries(R);
	}
}

Real tan(const Real& X)
{
	if (!X.IsFinite())
	{
		return Real::NaN();
	}
	long long Q = 0;
	const Real R = ReduceQuarter(X, &Q);
	const Real S = SinSeries(R), C = CosSeries(R);
	return (Q % 2 == 0) ? S / C : -C / S;
}

Real atan(const Real& X)
{
	static const Real HalfPi = Pi() / Real(2);
	if (X.IsNaN())
	{
		return X;
	}
	if (X.IsInf())
	{
		return X.IsNegative() ? -HalfPi : HalfPi;
	}
	Real Y = Real::Abs(X);
	const bool Invert = Y > Real(1);
	if (Invert)
	{
		Y = Real(1) / Y;
	}
	// Halve the angle twice: atan y = 2 atan(y / (1 + sqrt(1 + y^2))).
	for (int i = 0; i < 2; ++i)
	{
		Y = Y / (Real(1) + sqrt(Real(1) + Y * Y));
	}
	Real A = Real(4) * AtanSeries(Y);
	if (Invert)
	{
		A = HalfPi - A;
	}
	return X.IsNegative() ? -A : A;
}

Real atan2(const Real& Y, const Real& X)
{
	const Real& P = Pi();
	if (X.IsNaN() || Y.IsNaN())
	{
		return Real::NaN();
	}
	if (Y.IsInf() && X.IsInf())
	{
		const Real Q = X.IsNegative() ? Real(3) * P / Real(4) : P / Real(4);
		return Y.IsNegative() ? -Q : Q;
	}
	if (X.IsZero())
	{
		if (Y.IsZero())
		{
			return Real();
		}
		return Y.IsNegative() ? -P / Real(2) : P / Real(2);
	}
	if (X.IsInf())
	{
		if (!X.IsNegative())
		{
			return Real();
		}
		return Y.IsNegative() ? -P : P;
	}
	if (Y.IsInf())
	{
		return Y.IsNegative() ? -P / Real(2) : P / Real(2);
	}
	const Real A = atan(Y / X);
	if (!X.IsNegative())
	{
		return A;
	}
	return Y.IsNegative() ? A - P : A + P;
}

Real asin(const Real& X)
{
	if (Real::Abs(X) > Real(1) || X.IsNaN())
	{
		return Real::NaN();
	}
	return atan2(X, sqrt((Real(1) - X) * (Real(1) + X)));
}

Real acos(const Real& X)
{
	if (Real::Abs(X) > Real(1) || X.IsNaN())
	{
		return Real::NaN();
	}
	return atan2(sqrt((Real(1) - X) * (Real(1) + X)), X);
}

Real sinh(const Real& X)
{
	if (!X.IsFinite())
	{
		return X;
	}
	if (Real::Abs(X) < Real(0.5))
	{
		return SinhSeries(X);
	}
	const Real E = exp(Real::Abs(X));
	const Real S = (E - Real(1) / E) / Real(2);
	return X.IsNegative() ? -S : S;
}

Real cosh(const Real& X)
{
	if (X.IsNaN())
	{
		return X;
	}
	if (X.IsInf())
	{
		return Real::Infinity();
	}
	const Real E = exp(Real::Abs(X));
	return (E + Real(1) / E) / Real(2);
}

Real tanh(const Real& X)
{
	if (X.IsNaN())
	{
		return X;
	}
	const Real A = Real::Abs(X);
	if (A > Real(22))
	{
		return X.IsNegative() ? Real(-1) : Real(1);
	}
	const Real T = expm1(Real(2) * A);
	const Real R = T / (T + Real(2));
	return X.IsNegative() ? -R : R;
}

Real asinh(const Real& X)
{
	if (!X.IsFinite())
	{
		return X;
	}
	const Real A = Real::Abs(X);
	const Real R = A < Real(0.5) ? log1p(A + A * A / (Real(1) + sqrt(Real(1) + A * A))) : log(A + sqrt(A * A + Real(1)));
	return X.IsNegative() ? -R : R;
}

Real acosh(const Real& X)
{
	if (X.IsNaN() || X < Real(1))
	{
		return Real::NaN();
	}
	if (X.IsInf())
	{
		return X;
	}
	if (X < Real(2))
	{
		const Real T = X - Real(1);
		return log1p(T + sqrt(T * (X + Real(1))));
	}
	return log(X + sqrt(X * X - Real(1)));
}

Real atanh(const Real& X)
{
	if (X.IsNaN() || Real::Abs(X) > Real(1))
	{
		return Real::NaN();
	}
	if (Real::Abs(X) == Real(1))
	{
		return Real::Infinity(X.IsNegative());
	}
	return log1p(Real(2) * X / (Real(1) - X)) / Real(2);
}

Real erf(const Real& X)
{
	if (X.IsNaN())
	{
		return X;
	}
	const Real A = Real::Abs(X);
	if (A >= Real(2))
	{
		const Real R = Real(1) - erfc(A);
		return X.IsNegative() ? -R : R;
	}
	// erf x = 2/sqrt(pi) e^(-x^2) sum 2^n x^(2n+1) / (1 3 5 ... (2n+1)): every term positive.
	static const Real TwoOverRootPi = Real(2) / sqrt(Pi());
	const Real X2 = A * A;
	Real Term = A;
	RealAccumulator Acc;
	Acc.Add(A);
	Real Sum = A;
	for (long long N = 1; N < 400; ++N)
	{
		Term = Term * Real(2) * X2 / Real(2 * N + 1);
		Acc.Add(Term);
		Sum = Sum + Term;
		if (Negligible(Term, Sum))
		{
			break;
		}
	}
	const Real R = TwoOverRootPi * exp(-X2) * Acc.Result();
	return X.IsNegative() ? -R : R;
}

Real erfc(const Real& X)
{
	if (X.IsNaN())
	{
		return X;
	}
	if (X.IsNegative())
	{
		return Real(2) - erfc(-X);
	}
	if (X < Real(2))
	{
		return Real(1) - erf(X);
	}
	if (X.IsInf())
	{
		return Real();
	}
	// Continued fraction, evaluated from the tail:
	// erfc x = e^(-x^2)/sqrt(pi) * 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + 2/(x + ...)))))
	static const Real InvRootPi = Real(1) / sqrt(Pi());
	const int Depth = X < Real(3) ? 400 : (X < Real(6) ? 160 : 60);
	Real F = X;
	for (int K = Depth; K >= 1; --K)
	{
		F = X + Real(K) / Real(2) / F;
	}
	return InvRootPi * exp(-X * X) / F;
}

Real lgamma(const Real& X)
{
	if (X.IsNaN())
	{
		return X;
	}
	if (X.IsInf())
	{
		return Real::Infinity();
	}
	if (X <= Real(0) && Real::Trunc(X) == X)
	{
		return Real::Infinity();
	}
	if (X < Real(0.5))
	{
		// Reflection: Gamma(x) Gamma(1 - x) = pi / sin(pi x)
		return log(Pi() / Real::Abs(sin(Pi() * X))) - lgamma(Real(1) - X);
	}
	Real Y = X;
	Real Product(1);
	while (Y < Real(15))
	{
		Product = Product * Y;
		Y = Y + Real(1);
	}
	return LgammaStirling(Y) - log(Product);
}

Real tgamma(const Real& X)
{
	if (X.IsNaN())
	{
		return X;
	}
	if (X.IsZero())
	{
		return Real::Infinity();
	}
	if (X.IsInf())
	{
		return X.IsNegative() ? Real::NaN() : X;
	}
	if (X < Real(0))
	{
		if (Real::Trunc(X) == X)
		{
			return Real::NaN();
		}
		return Pi() / (sin(Pi() * X) * tgamma(Real(1) - X));
	}
	if (Real::Trunc(X) == X && X <= Real(30))
	{
		Real F(1);
		for (long long K = 2; K < X.ToLong(); ++K)
		{
			F = F * Real(K);
		}
		return F;
	}
	return exp(lgamma(X));
}

Real ldexp(const Real& X, int E)
{
	return X * powi(Real(2), E);
}
} // namespace gm
} // namespace mcalc
