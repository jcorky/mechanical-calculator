#include "GearReal.h"

#include "GearMachine.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace mcalc
{
namespace gear_detail
{
constexpr uint64_t kE18 = 1000000000000000000ull;
constexpr uint64_t kE17 = 100000000000000000ull;
constexpr uint64_t kE9 = 1000000000ull;
constexpr int kExpLimit = 100000000;   // beyond this the machine's notes overflow: inf or 0

const uint64_t kPow10[19] = {1ull, 10ull, 100ull, 1000ull, 10000ull, 100000ull, 1000000ull, 10000000ull, 100000000ull, 1000000000ull,
	10000000000ull, 100000000000ull, 1000000000000ull, 10000000000000ull, 100000000000000ull, 1000000000000000ull, 10000000000000000ull,
	100000000000000000ull, 1000000000000000000ull};

int DigitCount(uint64_t X)
{
	int N = 0;
	while (X > 0)
	{
		++N;
		X /= 10;
	}
	return N;
}

// A decimal integer of up to 54 digits in base-10^18 limbs, L[0] lowest.
struct D3
{
	uint64_t L[3] = {0, 0, 0};

	static D3 FromU64(uint64_t X)
	{
		D3 D;
		D.L[0] = X % kE18;
		D.L[1] = X / kE18;
		return D;
	}

	bool IsZero() const { return L[0] == 0 && L[1] == 0 && L[2] == 0; }

	int Digits() const
	{
		for (int i = 2; i >= 0; --i)
		{
			if (L[i] != 0)
			{
				return 18 * i + DigitCount(L[i]);
			}
		}
		return 0;
	}

	// Multiplies by M < 10^9.
	void MulSmall(uint64_t M)
	{
		uint64_t Carry = 0;
		for (uint64_t& X : L)
		{
			const uint64_t Hi = X / kE9, Lo = X % kE9;
			const uint64_t HiM = Hi * M;
			const uint64_t Low = Lo * M + (HiM % kE9) * kE9 + Carry;
			X = Low % kE18;
			Carry = HiM / kE9 + Low / kE18;
		}
	}

	void AddSmall(uint64_t M)
	{
		uint64_t Carry = M;
		for (uint64_t& X : L)
		{
			const uint64_t V = X + Carry;
			X = V % kE18;
			Carry = V / kE18;
			if (Carry == 0)
			{
				break;
			}
		}
	}

	// Multiplies by 10^K, K <= 36.
	void Shift(int K)
	{
		while (K >= 18)
		{
			L[2] = L[1];
			L[1] = L[0];
			L[0] = 0;
			K -= 18;
		}
		while (K > 0)
		{
			const int Step = std::min(K, 9);
			MulSmall(kPow10[Step]);
			K -= Step;
		}
	}

	static int Cmp(const D3& A, const D3& B)
	{
		for (int i = 2; i >= 0; --i)
		{
			if (A.L[i] != B.L[i])
			{
				return A.L[i] < B.L[i] ? -1 : 1;
			}
		}
		return 0;
	}

	static D3 Add(const D3& A, const D3& B)
	{
		D3 R;
		uint64_t Carry = 0;
		for (int i = 0; i < 3; ++i)
		{
			const uint64_t V = A.L[i] + B.L[i] + Carry;
			R.L[i] = V % kE18;
			Carry = V / kE18;
		}
		return R;
	}

	static D3 Sub(const D3& A, const D3& B)   // A >= B
	{
		D3 R;
		uint64_t Borrow = 0;
		for (int i = 0; i < 3; ++i)
		{
			const uint64_t Need = B.L[i] + Borrow;
			if (A.L[i] >= Need)
			{
				R.L[i] = A.L[i] - Need;
				Borrow = 0;
			}
			else
			{
				R.L[i] = A.L[i] + kE18 - Need;
				Borrow = 1;
			}
		}
		return R;
	}

	static D3 MulU64(uint64_t A, uint64_t B)   // A, B < 10^18
	{
		const uint64_t A1 = A / kE9, A0 = A % kE9, B1 = B / kE9, B0 = B % kE9;
		const uint64_t Mid = A1 * B0 + A0 * B1;   // < 2e18
		uint64_t Lo = A0 * B0 + (Mid % kE9) * kE9;   // < 2e18
		D3 R;
		R.L[0] = Lo % kE18;
		R.L[1] = A1 * B1 + Mid / kE9 + Lo / kE18;
		return R;
	}

	// Drops the lowest K digits (K <= 36); returns the most significant dropped digit.
	int Drop(int K)
	{
		int Top = 0;
		while (K > 0)
		{
			const int Step = std::min(K, 18);
			const uint64_t Div = kPow10[Step];
			uint64_t Rem = 0;
			for (int i = 2; i >= 0; --i)
			{
				const uint64_t X = L[i];
				L[i] = Rem * kPow10[18 - Step] + X / Div;
				Rem = X % Div;
			}
			Top = static_cast<int>(Rem / kPow10[Step - 1]);
			K -= Step;
		}
		return Top;
	}
};

// The operator's rounding: the exact register contents N times 10^LowExp, to 18 digits, half up.
Real RoundExact(D3 N, long long LowExp, bool Negative)
{
	const int Count = N.Digits();
	if (Count == 0)
	{
		return Real();
	}
	int K = Count - Real::Digits;
	uint64_t Mant = 0;
	if (K <= 0)
	{
		N.Shift(-K);
		Mant = N.L[0];
	}
	else
	{
		const int Dropped = N.Drop(K);
		Mant = N.L[0];
		if (Dropped >= 5)
		{
			++Mant;
			if (Mant == kE18)
			{
				Mant = kE17;
				++K;
			}
		}
	}
	const long long E = LowExp + K + 17;
	if (E > kExpLimit)
	{
		return Real::Infinity(Negative);
	}
	if (E < -kExpLimit)
	{
		return Real();
	}
	return Real::FromParts(Negative, Mant, static_cast<int>(E));
}

Real RoundDigits(const std::vector<uint8_t>& Digits, long long LowExp, bool Negative)
{
	D3 N;
	for (size_t i = Digits.size(); i-- > 0;)
	{
		N.MulSmall(10);
		N.AddSmall(Digits[i]);
	}
	return RoundExact(N, LowExp, Negative);
}

bool UseMachine()
{
	GearRecorder* R = GearContext::Attached();
	return R && R->Wants();
}

thread_local GearStats TStats;
thread_local GearRecorder* TRecorder = nullptr;

uint64_t DigitSum(uint64_t X)
{
	uint64_t S = 0;
	while (X > 0)
	{
		S += X % 10;
		X /= 10;
	}
	return S;
}
} // namespace gear_detail

using namespace gear_detail;

namespace GearContext
{
GearStats& Stats()
{
	return TStats;
}

void ResetStats()
{
	TStats = GearStats();
}

void Attach(GearRecorder* Recorder)
{
	TRecorder = Recorder;
}

GearRecorder* Attached()
{
	return TRecorder;
}
} // namespace GearContext

Real Real::FromParts(bool Negative, uint64_t Mantissa, int Exponent)
{
	Real R;
	if (Mantissa == 0)
	{
		return R;
	}
	// Normalize to 18 digits (callers pass normalized values; this keeps FromParts total).
	while (Mantissa >= kE18)
	{
		const uint64_t Last = Mantissa % 10;
		Mantissa /= 10;
		if (Last >= 5)
		{
			++Mantissa;
		}
		++Exponent;
	}
	while (Mantissa < kE17)
	{
		Mantissa *= 10;
		--Exponent;
	}
	R.Mant = Mantissa;
	R.Exp = Exponent;
	R.Neg = Negative;
	return R;
}

Real Real::Infinity(bool Negative)
{
	Real R;
	R.Kind = KInf;
	R.Neg = Negative;
	return R;
}

Real Real::NaN()
{
	Real R;
	R.Kind = KNaN;
	return R;
}

Real Real::Epsilon()
{
	return FromParts(false, kE17, -17);
}

void Real::SetInteger(long long X)
{
	const bool Negative = X < 0;
	const unsigned long long U = Negative ? 0ull - static_cast<unsigned long long>(X) : static_cast<unsigned long long>(X);
	SetUnsigned(U);
	if (Mant != 0)
	{
		Neg = Negative;
	}
}

void Real::SetUnsigned(unsigned long long X)
{
	Kind = KFinite;
	Neg = false;
	if (X == 0)
	{
		Mant = 0;
		Exp = 0;
		return;
	}
	*this = RoundExact(D3::FromU64(X), 0, false);
}

Real::Real(double X)
{
	if (std::isnan(X))
	{
		*this = NaN();
		return;
	}
	if (std::isinf(X))
	{
		*this = Infinity(X < 0);
		return;
	}
	if (X == std::floor(X) && std::fabs(X) < 9.0e15)
	{
		SetInteger(static_cast<long long>(X));
		return;
	}
	// Literals repeat in loops, so remember the conversions.
	thread_local std::unordered_map<uint64_t, Real> Cache;
	uint64_t Bits = 0;
	static_assert(sizeof(Bits) == sizeof(X), "double is 64 bits");
	std::memcpy(&Bits, &X, sizeof(Bits));
	auto It = Cache.find(Bits);
	if (It != Cache.end())
	{
		*this = It->second;
		return;
	}
	char Buf[40];
	for (int P = 1; P <= 17; ++P)
	{
		std::snprintf(Buf, sizeof(Buf), "%.*e", P - 1, X);
		if (std::strtod(Buf, nullptr) == X)
		{
			break;
		}
	}
	size_t Used = 0;
	*this = Parse(Buf, &Used);
	if (Cache.size() > 8192)
	{
		Cache.clear();
	}
	Cache.emplace(Bits, *this);
}

Real Real::Parse(const char* Text, size_t* Consumed)
{
	const char* P = Text;
	bool Negative = false;
	if (*P == '+' || *P == '-')
	{
		Negative = *P == '-';
		++P;
	}
	D3 N;
	int Sig = 0;
	long long Scale = 0;      // power of ten of N's last digit
	bool Any = false;
	bool Dot = false;
	int Guard = -1;           // the 19th significant digit, for rounding
	for (;; ++P)
	{
		if (*P == '.' && !Dot)
		{
			Dot = true;
			continue;
		}
		if (*P < '0' || *P > '9')
		{
			break;
		}
		Any = true;
		const int D = *P - '0';
		if (Sig == 0 && D == 0)
		{
			if (Dot)
			{
				--Scale;
			}
			continue;
		}
		if (Sig < Real::Digits)
		{
			N.MulSmall(10);
			N.AddSmall(static_cast<uint64_t>(D));
			++Sig;
			if (Dot)
			{
				--Scale;
			}
		}
		else
		{
			if (Guard < 0)
			{
				Guard = D;
			}
			if (!Dot)
			{
				++Scale;
			}
		}
	}
	if (!Any)
	{
		if (Consumed)
		{
			*Consumed = 0;
		}
		return Real();
	}
	if ((*P == 'e' || *P == 'E') && (std::isdigit(static_cast<unsigned char>(P[1])) ||
		((P[1] == '+' || P[1] == '-') && std::isdigit(static_cast<unsigned char>(P[2])))))
	{
		char* End = nullptr;
		const long E = std::strtol(P + 1, &End, 10);
		Scale += E;
		P = End;
	}
	if (Consumed)
	{
		*Consumed = static_cast<size_t>(P - Text);
	}
	if (Guard >= 5)
	{
		N.AddSmall(1);
	}
	Real R = RoundExact(N, Scale, Negative);
	return R;
}

Real Real::operator-() const
{
	Real R = *this;
	if (!R.IsZero() && !R.IsNaN())
	{
		R.Neg = !R.Neg;
	}
	return R;
}

int Real::Compare(const Real& A, const Real& B)
{
	if (A.IsNaN() || B.IsNaN())
	{
		return 2;
	}
	const int SA = A.IsZero() ? 0 : (A.Neg ? -1 : 1);
	const int SB = B.IsZero() ? 0 : (B.Neg ? -1 : 1);
	if (SA != SB)
	{
		return SA < SB ? -1 : 1;
	}
	if (SA == 0)
	{
		return 0;
	}
	// Same sign: compare magnitudes, then flip for negatives.
	int M = 0;
	if (A.IsInf() || B.IsInf())
	{
		M = A.IsInf() == B.IsInf() ? 0 : (A.IsInf() ? 1 : -1);
	}
	else if (A.Exp != B.Exp)
	{
		M = A.Exp < B.Exp ? -1 : 1;
	}
	else if (A.Mant != B.Mant)
	{
		M = A.Mant < B.Mant ? -1 : 1;
	}
	return SA > 0 ? M : -M;
}

Real Real::Add(const Real& InA, const Real& InB, bool Subtract)
{
	Real B = InB;
	if (Subtract && !B.IsNaN() && !B.IsZero())
	{
		B.Neg = !B.Neg;
	}
	const Real& A = InA;
	if (A.IsNaN() || B.IsNaN())
	{
		return NaN();
	}
	if (A.IsInf() || B.IsInf())
	{
		if (A.IsInf() && B.IsInf() && A.Neg != B.Neg)
		{
			return NaN();
		}
		return A.IsInf() ? A : B;
	}
	if (B.IsZero())
	{
		return A;
	}
	if (A.IsZero())
	{
		return B;
	}
	++TStats.Additions;
	// The operand with the larger exponent goes in first, shifted up to align the points.
	const Real& Hi = A.Exp >= B.Exp ? A : B;
	const Real& Lo = A.Exp >= B.Exp ? B : A;
	const long long Gap = static_cast<long long>(Hi.Exp) - Lo.Exp;
	if (Gap > 20)
	{
		// The carriage cannot shift that far: the smaller operand is below the 40th wheel.
		return Hi;
	}
	const bool Same = Hi.Neg == Lo.Neg;
	D3 X = D3::FromU64(Hi.Mant);
	X.Shift(static_cast<int>(Gap));
	const D3 Y = D3::FromU64(Lo.Mant);
	const long long LowExp = static_cast<long long>(Lo.Exp) - 17;
	if (Same)
	{
		if (UseMachine())
		{
			return RoundDigits(GearProcedures::AddAligned(Hi.Mant, static_cast<int>(Gap), Lo.Mant, false, GearContext::Attached()), LowExp, Hi.Neg);
		}
		TStats.Turns += 2;
		TStats.Shifts += static_cast<unsigned long long>(2 * Gap);
		return RoundExact(D3::Add(X, Y), LowExp, Hi.Neg);
	}
	// Opposite signs: the operator subtracts the smaller magnitude from the larger.
	const int C = D3::Cmp(X, Y);
	if (C == 0)
	{
		TStats.Turns += 2;
		return Real();
	}
	if (UseMachine())
	{
		if (C > 0)
		{
			return RoundDigits(GearProcedures::AddAligned(Hi.Mant, static_cast<int>(Gap), Lo.Mant, true, GearContext::Attached()), LowExp, Hi.Neg);
		}
		// |Lo| > |Hi| only happens with Gap == 0.
		return RoundDigits(GearProcedures::AddAligned(Lo.Mant, 0, Hi.Mant, true, GearContext::Attached()), LowExp, Lo.Neg);
	}
	TStats.Turns += 2;
	TStats.Shifts += static_cast<unsigned long long>(2 * Gap);
	return C > 0 ? RoundExact(D3::Sub(X, Y), LowExp, Hi.Neg) : RoundExact(D3::Sub(Y, X), LowExp, Lo.Neg);
}

Real Real::Mul(const Real& A, const Real& B)
{
	if (A.IsNaN() || B.IsNaN())
	{
		return NaN();
	}
	const bool Negative = A.Neg != B.Neg;
	if (A.IsInf() || B.IsInf())
	{
		if (A.IsZero() || B.IsZero())
		{
			return NaN();
		}
		return Infinity(Negative);
	}
	if (A.IsZero() || B.IsZero())
	{
		return Real();
	}
	++TStats.Multiplications;
	const long long LowExp = static_cast<long long>(A.Exp) - 17 + B.Exp - 17;
	if (UseMachine())
	{
		return RoundDigits(GearProcedures::Multiply(A.Mant, B.Mant, GearContext::Attached()), LowExp, Negative);
	}
	TStats.Turns += DigitSum(B.Mant);
	TStats.Shifts += 34;
	return RoundExact(D3::MulU64(A.Mant, B.Mant), LowExp, Negative);
}

Real Real::Div(const Real& A, const Real& B)
{
	if (A.IsNaN() || B.IsNaN())
	{
		return NaN();
	}
	const bool Negative = A.Neg != B.Neg;
	if (A.IsInf())
	{
		return B.IsInf() ? NaN() : Infinity(Negative);
	}
	if (B.IsInf())
	{
		return Real();
	}
	if (B.IsZero())
	{
		// The bell rings and never stops.
		return A.IsZero() ? NaN() : Infinity(A.Neg);
	}
	if (A.IsZero())
	{
		return Real();
	}
	++TStats.Divisions;
	if (UseMachine())
	{
		bool Smaller = false;
		const std::vector<uint8_t> Q = GearProcedures::Divide(A.Mant, B.Mant, GearContext::Attached(), &Smaller);
		return RoundDigits(Q, static_cast<long long>(A.Exp) - B.Exp - (Smaller ? 19 : 18), Negative);
	}
	// Repeated subtraction, place by place, exactly as on the machine.
	uint64_t R = A.Mant;
	uint64_t Q = 0;
	int Got = 0;
	bool Smaller = false;
	for (int Place = 0; Got < 19; ++Place)
	{
		uint64_t D = 0;
		while (R >= B.Mant)
		{
			R -= B.Mant;
			++D;
		}
		TStats.Turns += D + 2;
		if (Place == 0 && D == 0)
		{
			Smaller = true;
		}
		else
		{
			Q = Q * 10 + D;
			++Got;
		}
		R *= 10;
	}
	TStats.Shifts += 40;
	return RoundExact(D3::FromU64(Q), static_cast<long long>(A.Exp) - B.Exp - (Smaller ? 19 : 18), Negative);
}

Real Real::Sqrt(const Real& A)
{
	if (A.IsNaN() || (A.Neg && !A.IsZero()))
	{
		return NaN();
	}
	if (A.IsInf() || A.IsZero())
	{
		return A;
	}
	++TStats.Roots;
	// Make the power of ten even, then extend N to 37 or 38 digits for a 19-digit root.
	long long P = static_cast<long long>(A.Exp) - 17;
	D3 N = D3::FromU64(A.Mant);
	int Len = 18;
	if (P % 2 != 0)
	{
		N.MulSmall(10);
		--P;
		Len = 19;
	}
	const int Extra = Len == 18 ? 20 : 18;
	N.Shift(Extra);
	const long long LowExp = P / 2 - Extra / 2;
	// Digits of N, units first, padded to an even count.
	std::vector<uint8_t> Digits;
	{
		D3 T = N;
		const int Count = T.Digits();
		for (int i = 0; i < Count; ++i)
		{
			Digits.push_back(static_cast<uint8_t>(T.L[0] % 10));
			T.Drop(1);
		}
		if (Digits.size() % 2 != 0)
		{
			Digits.push_back(0);
		}
	}
	if (UseMachine())
	{
		return RoundDigits(GearProcedures::SquareRoot(Digits, GearContext::Attached()), LowExp, false);
	}
	// Toepler: at each column pair subtract 20y + 1, 20y + 3, ... while it goes.
	D3 Rem, Y;
	for (size_t j = Digits.size() / 2; j-- > 0;)
	{
		Rem.MulSmall(100);
		Rem.AddSmall(static_cast<uint64_t>(Digits[2 * j + 1] * 10 + Digits[2 * j]));
		D3 Base = Y;
		Base.MulSmall(20);
		uint64_t X = 0;
		D3 Odd = Base;
		Odd.AddSmall(1);
		while (D3::Cmp(Rem, Odd) >= 0)
		{
			Rem = D3::Sub(Rem, Odd);
			Odd.AddSmall(2);
			++X;
		}
		TStats.Turns += X + 2;
		Y.MulSmall(10);
		Y.AddSmall(X);
	}
	TStats.Shifts += 40;
	return RoundExact(Y, LowExp, false);
}

RealAccumulator::RealAccumulator()
{
	for (int i = 0; i < Limbs; ++i)
	{
		Pos[i] = 0;
		NegPart[i] = 0;
	}
}

void RealAccumulator::Rebase(long long NewLowExp)
{
	// Drop digits from the bottom of both parts (the wheels below the 40th place).
	long long Drop = NewLowExp - LowExp;
	for (uint64_t* Part : {Pos, NegPart})
	{
		long long K = Drop;
		while (K > 0)
		{
			const int Step = static_cast<int>(std::min<long long>(K, 9));
			const uint64_t Div = kPow10[Step];
			uint64_t Rem = 0;
			for (int i = Limbs - 1; i >= 0; --i)
			{
				const uint64_t X = Part[i] + Rem * kE9;
				Part[i] = X / Div;
				Rem = X % Div;
			}
			K -= Step;
		}
	}
	LowExp = NewLowExp;
}

void RealAccumulator::AddDigits(const uint64_t* Value18, long long ValueLowExp, bool Negative)
{
	// Value as base-1e9 limbs (up to 36 digits -> 4 limbs).
	uint64_t V[Limbs + 2] = {};
	V[0] = Value18[0] % kE9;
	V[1] = Value18[0] / kE9;
	V[2] = Value18[1] % kE9;
	V[3] = Value18[1] / kE9;
	int Top = 0;
	for (int i = 3; i >= 0; --i)
	{
		if (V[i] != 0)
		{
			Top = 9 * i + DigitCount(V[i]);
			break;
		}
	}
	if (Top == 0)
	{
		return;
	}
	const long long ValueTop = ValueLowExp + Top;   // one past the highest digit's power
	if (!Started)
	{
		Started = true;
		LowExp = ValueTop - 40;
	}
	else if (ValueTop > LowExp + 40)
	{
		Rebase(ValueTop - 40);
	}
	long long Shift = ValueLowExp - LowExp;
	// Below the bottom wheel: drop those digits.
	while (Shift < 0)
	{
		const int Step = static_cast<int>(std::min<long long>(-Shift, 9));
		uint64_t Rem = 0;
		for (int i = Limbs + 1; i >= 0; --i)
		{
			const uint64_t X = V[i] + Rem * kE9;
			V[i] = X / kPow10[Step];
			Rem = X % kPow10[Step];
		}
		Shift += Step;
	}
	// Shift up by whole limbs, then by the remaining digits.
	const int LimbShift = static_cast<int>(Shift / 9);
	const uint64_t Mul = kPow10[Shift % 9];
	uint64_t W[Limbs + 2] = {};
	uint64_t Carry = 0;
	for (int i = 0; i < Limbs + 2; ++i)
	{
		const uint64_t X = V[i] * Mul + Carry;   // V[i] < 1e9, Mul < 1e9
		W[i] = X % kE9;
		Carry = X / kE9;
	}
	uint64_t* Part = Negative ? NegPart : Pos;
	Carry = 0;
	for (int i = 0; i < Limbs; ++i)
	{
		const int j = i - LimbShift;
		const uint64_t A = (j >= 0 && j < Limbs + 2) ? W[j] : 0;
		const uint64_t X = Part[i] + A + Carry;
		Part[i] = X % kE9;
		Carry = X / kE9;
	}
}

void RealAccumulator::Add(const Real& Term)
{
	if (Special || Term.IsZero())
	{
		return;
	}
	if (!Term.IsFinite())
	{
		Special = true;
		SpecialValue = Term;
		return;
	}
	++TStats.Additions;
	++TStats.Turns;
	const uint64_t Value[2] = {Term.Mantissa(), 0};
	AddDigits(Value, static_cast<long long>(Term.Exponent()) - 17, Term.IsNegative());
}

void RealAccumulator::AddProduct(const Real& A, const Real& B, bool Subtract)
{
	if (Special || A.IsZero() || B.IsZero())
	{
		return;
	}
	if (!A.IsFinite() || !B.IsFinite())
	{
		Special = true;
		SpecialValue = Subtract ? -(A * B) : A * B;
		return;
	}
	++TStats.Multiplications;
	TStats.Turns += DigitSum(B.Mantissa());
	const D3 P = D3::MulU64(A.Mantissa(), B.Mantissa());
	const uint64_t Value[2] = {P.L[0], P.L[1]};
	AddDigits(Value, static_cast<long long>(A.Exponent()) - 17 + B.Exponent() - 17, (A.IsNegative() != B.IsNegative()) != Subtract);
}

Real RealAccumulator::Result() const
{
	if (Special)
	{
		return SpecialValue;
	}
	// Difference of the parts, then one rounding to 18 digits.
	int C = 0;
	for (int i = Limbs - 1; i >= 0 && C == 0; --i)
	{
		if (Pos[i] != NegPart[i])
		{
			C = Pos[i] > NegPart[i] ? 1 : -1;
		}
	}
	if (C == 0)
	{
		return Real();
	}
	const uint64_t* Big = C > 0 ? Pos : NegPart;
	const uint64_t* Small = C > 0 ? NegPart : Pos;
	uint64_t Diff[Limbs];
	uint64_t Borrow = 0;
	for (int i = 0; i < Limbs; ++i)
	{
		const uint64_t Need = Small[i] + Borrow;
		if (Big[i] >= Need)
		{
			Diff[i] = Big[i] - Need;
			Borrow = 0;
		}
		else
		{
			Diff[i] = Big[i] + kE9 - Need;
			Borrow = 1;
		}
	}
	// Convert the 54-digit base-1e9 number to D3 (base 1e18) and round.
	D3 N;
	N.L[0] = Diff[0] + Diff[1] * kE9;
	N.L[1] = Diff[2] + Diff[3] * kE9;
	N.L[2] = Diff[4] + Diff[5] * kE9;
	return RoundExact(N, LowExp, C < 0);
}

Real Real::Scale10(const Real& A, int Power)
{
	if (!A.IsFinite() || A.IsZero())
	{
		return A;
	}
	const long long E = static_cast<long long>(A.Exp) + Power;
	if (E > kExpLimit)
	{
		return Infinity(A.Neg);
	}
	if (E < -kExpLimit)
	{
		return Real();
	}
	return FromParts(A.Neg, A.Mant, static_cast<int>(E));
}

Real Real::Trunc(const Real& A)
{
	if (!A.IsFinite() || A.IsZero() || A.Exp >= 17)
	{
		return A;
	}
	if (A.Exp < 0)
	{
		return Real();
	}
	const uint64_t Unit = kPow10[17 - A.Exp];
	return FromParts(A.Neg, A.Mant / Unit * Unit, A.Exp);
}

Real Real::Floor(const Real& A)
{
	const Real T = Trunc(A);
	return (A.IsNegative() && T != A) ? T - Real(1) : T;
}

Real Real::Ceil(const Real& A)
{
	const Real T = Trunc(A);
	return (!A.IsNegative() && T != A) ? T + Real(1) : T;
}

Real Real::Round(const Real& A)
{
	if (!A.IsFinite() || A.IsZero() || A.Exp >= 17)
	{
		return A;
	}
	const Real Half = Real::FromParts(false, 5 * kE17, -1);   // 0.5
	const Real Mag = Trunc(Abs(A) + Half);
	return A.Neg ? -Mag : Mag;
}

long long Real::ToLong() const
{
	if (IsNaN())
	{
		return 0;
	}
	if (IsInf() || (!IsZero() && Exp >= 18))
	{
		return Neg ? -9223372036854775807LL : 9223372036854775807LL;
	}
	if (IsZero() || Exp < 0)
	{
		return 0;
	}
	const long long V = static_cast<long long>(Mant / kPow10[17 - Exp]);
	return Neg ? -V : V;
}

int Real::ToInt() const
{
	const long long V = ToLong();
	return static_cast<int>(std::max<long long>(-2147483647LL, std::min<long long>(2147483647LL, V)));
}

std::string Real::ToString(int Sig) const
{
	if (IsNaN())
	{
		return "nan";
	}
	if (IsInf())
	{
		return Neg ? "-inf" : "inf";
	}
	if (IsZero())
	{
		return "0";
	}
	Sig = std::max(1, std::min(Real::Digits, Sig));
	uint64_t M = Mant;
	int E = Exp;
	if (Sig < Real::Digits)
	{
		const uint64_t Unit = kPow10[Real::Digits - Sig];
		const uint64_t Rem = M % Unit;
		M = M / Unit;
		if (Rem >= Unit / 2)
		{
			++M;
			if (M == kPow10[Sig])
			{
				M /= 10;
				++E;
			}
		}
	}
	else
	{
		M = Mant;
	}
	std::string D = std::to_string(M);   // Sig digits
	while (D.size() > 1 && D.back() == '0')
	{
		D.pop_back();
	}
	std::string Out = Neg ? "-" : "";
	if (E >= -5 && E < Sig)
	{
		if (E >= 0)
		{
			std::string Int = D.substr(0, std::min(D.size(), static_cast<size_t>(E + 1)));
			while (static_cast<int>(Int.size()) < E + 1)
			{
				Int += '0';
			}
			Out += Int;
			if (D.size() > static_cast<size_t>(E + 1))
			{
				Out += "." + D.substr(static_cast<size_t>(E + 1));
			}
		}
		else
		{
			Out += "0." + std::string(static_cast<size_t>(-E - 1), '0') + D;
		}
		return Out;
	}
	Out += D.substr(0, 1);
	if (D.size() > 1)
	{
		Out += "." + D.substr(1);
	}
	return Out + "e" + std::to_string(E);
}

double Real::ToDouble() const
{
	if (IsNaN())
	{
		return std::nan("");
	}
	if (IsInf())
	{
		return Neg ? -HUGE_VAL : HUGE_VAL;
	}
	if (IsZero())
	{
		return 0.0;
	}
	char Buf[48];
	std::snprintf(Buf, sizeof(Buf), "%s%llue%d", Neg ? "-" : "", static_cast<unsigned long long>(Mant), Exp - 17);
	return std::strtod(Buf, nullptr);
}
} // namespace mcalc
