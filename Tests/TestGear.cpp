// The gear arithmetic: the literal stepped-drum machine and the fast path must agree digit for
// digit; the elementary functions built from the machine's five procedures must agree with
// <cmath> to the precision of 18-digit decimal arithmetic.
#include "GearMachine.h"
#include "GearMath.h"
#include "TestFramework.h"

#include <cmath>
#include <cstdint>
#include <cstdio>

using mcalc::Real;
namespace gm = mcalc::gm;

namespace
{
// xorshift64*, for reproducible operands
struct Rng
{
	uint64_t S = 0x9E3779B97F4A7C15ull;
	uint64_t Next()
	{
		S ^= S >> 12;
		S ^= S << 25;
		S ^= S >> 27;
		return S * 2685821657736338717ull;
	}
	// A random 18-digit mantissa, sometimes with long runs of 9s or 0s to stress the carries.
	uint64_t Mantissa()
	{
		const uint64_t Kind = Next() % 8;
		if (Kind == 0)
		{
			return 999999999999999999ull - Next() % 1000;
		}
		if (Kind == 1)
		{
			return 100000000000000000ull + Next() % 1000;
		}
		return 100000000000000000ull + Next() % 900000000000000000ull;
	}
	Real Value()
	{
		const bool Neg = Next() % 2 == 0;
		const int Exp = static_cast<int>(Next() % 41) - 20;
		return Real::FromParts(Neg, Mantissa(), Exp);
	}
};

bool Same(const Real& A, const Real& B)
{
	return A.Mantissa() == B.Mantissa() && A.Exponent() == B.Exponent() && A.IsNegative() == B.IsNegative() && A.IsFinite() == B.IsFinite();
}

double Rel(const Real& Got, double Want)
{
	const double G = Got.ToDouble();
	if (Want == 0.0)
	{
		return std::fabs(G);
	}
	return std::fabs(G - Want) / std::fabs(Want);
}
} // namespace

TEST(gear_literal_machine_equals_fast_path)
{
	Rng R;
	int Mismatches = 0;
	for (int i = 0; i < 3000; ++i)
	{
		const Real A = R.Value();
		Real B = R.Value();
		if (i % 7 == 0)
		{
			B = Real::FromParts(R.Next() % 2 == 0, A.Mantissa(), A.Exponent() - static_cast<int>(R.Next() % 3));   // near-cancellation
		}
		const Real Fast[5] = {A + B, A - B, A * B, A / B, Real::Sqrt(Real::Abs(A))};
		mcalc::GearRecorder Literal(0, false);   // keeps no events, but every operation turns the wheels
		mcalc::GearContext::Attach(&Literal);
		const Real Slow[5] = {A + B, A - B, A * B, A / B, Real::Sqrt(Real::Abs(A))};
		mcalc::GearContext::Attach(nullptr);
		for (int k = 0; k < 5; ++k)
		{
			if (!Same(Fast[k], Slow[k]))
			{
				if (++Mismatches < 6)
				{
					std::printf("    op %d: %s vs %s on %s, %s\n", k, Fast[k].ToString().c_str(), Slow[k].ToString().c_str(), A.ToString().c_str(), B.ToString().c_str());
				}
			}
		}
	}
	CHECK(Mismatches == 0);
}

TEST(gear_machine_registers_and_bell)
{
	// 123 * 45 on the wheels: the counter shows the multiplier, the accumulator the product.
	mcalc::GearRecorder Rec;
	const std::vector<uint8_t> P = mcalc::GearProcedures::Multiply(123, 45, &Rec);
	unsigned long long V = 0;
	for (size_t i = P.size(); i-- > 0;)
	{
		V = V * 10 + P[i];
	}
	CHECK(V == 5535);
	const mcalc::GearEvent& Last = Rec.GetEvents().back();
	CHECK(Last.Type == mcalc::GearEventType::End);
	CHECK(Last.Counter[0] == 5 && Last.Counter[1] == 4);
	int Turns = 0, Bells = 0;
	for (const mcalc::GearEvent& E : Rec.GetEvents())
	{
		Turns += E.Type == mcalc::GearEventType::Turn;
	}
	CHECK(Turns == 9);   // 5 turns, shift, 4 turns
	// Division rings the bell once per quotient place (the overdraft that ends each place).
	mcalc::GearRecorder Div;
	bool Smaller = false;
	const std::vector<uint8_t> Q = mcalc::GearProcedures::Divide(100000000000000000ull, 300000000000000000ull, &Div, &Smaller);
	CHECK(Smaller);
	for (const mcalc::GearEvent& E : Div.GetEvents())
	{
		Bells += E.Type == mcalc::GearEventType::Bell;
	}
	CHECK(Bells == 20);
	CHECK(Q.size() == 19 && Q[18] == 3 && Q[0] == 3);
	// Carries ripple: 999 + 1.
	mcalc::GearRecorder Add;
	const std::vector<uint8_t> S = mcalc::GearProcedures::AddAligned(999, 0, 1, false, &Add);
	CHECK(S[0] == 0 && S[1] == 0 && S[2] == 0 && S[3] == 1);
	bool Carried = false;
	for (const mcalc::GearEvent& E : Add.GetEvents())
	{
		Carried = Carried || (E.Type == mcalc::GearEventType::Turn && E.CarryMask == 0x7);
	}
	CHECK(Carried);
}

TEST(gear_number_entry_and_reading)
{
	size_t Used = 0;
	CHECK(Real::Parse("0.1", &Used).ToString() == "0.1" && Used == 3);
	CHECK(Real(0.1).ToString() == "0.1");
	CHECK(Real(1e-300).ToString() == "1e-300");
	CHECK(Real(6678e3).ToString() == "6678000");
	CHECK(Real::Parse("1.2345678901234567891", &Used).ToString() == "1.23456789012345679");   // 18 digits, half up
	CHECK(Real::Parse("2.5e-7", &Used).ToString(10) == "2.5e-7");
	CHECK(Real(-3).ToString() == "-3");
	CHECK((Real(1) / Real(3) * Real(3)).ToString() == "0.999999999999999999");   // eighteen 3s times 3, exactly
	CHECK((Real(1) / Real(0)).IsInf());
	CHECK((Real(0) / Real(0)).IsNaN());
	CHECK(!(Real::NaN() == Real::NaN()));
	CHECK(Real(2.5).ToLong() == 2 && Real(-2.5).ToLong() == -2);
	CHECK(gm::round(Real(2.5)) == Real(3) && gm::round(Real(-2.5)) == Real(-3));
	CHECK(gm::floor(Real(-0.5)) == Real(-1) && gm::ceil(Real(0.5)) == Real(1));
	CHECK(Real(0.1) + Real(0.2) == Real(0.3));   // decimal machine: exact
	CHECK_NEAR(Real(1.0 / 3.0).ToDouble(), 1.0 / 3.0, 0.0);
	CHECK(Real::Scale10(Real(5), 300).ToString() == "5e300");
	CHECK((Real::Scale10(Real(5), 300) * Real::Scale10(Real(5), 300)).ToString() == "2.5e601");   // no overflow on paper
}

TEST(gear_constants_computed_on_the_gears)
{
	CHECK(gm::Pi().ToString() == "3.14159265358979324");
	CHECK(gm::Ln2().ToString() == "0.693147180559945309");
	CHECK(gm::Ln10().ToString() == "2.30258509299404568");
}

TEST(gear_elementary_functions)
{
	// Dyadic values, exact in both binary and decimal, so <cmath> sees the same inputs.
	const double Xs[] = {-30.5, -7.25, -2.0, -1.0, -0.625, -0.0009765625, 9.5367431640625e-7, 0.0078125, 0.296875, 0.5, 0.90625, 1.0, 1.6875, 2.0, 3.140625, 6.5, 10.0, 77.75, 1234.5};
	double Worst = 0.0;
	const char* WorstName = "";
	double WorstX = 0.0;
	auto Note = [&](const char* Name, double X, const Real& Got, double Want) {
		if (!std::isfinite(Want))
		{
			return;
		}
		const double E = Rel(Got, Want);
		if (E > Worst)
		{
			Worst = E;
			WorstName = Name;
			WorstX = X;
		}
	};
	for (double X : Xs)
	{
		const Real R(X);
		Note("exp", X, gm::exp(R), std::exp(X));
		Note("sin", X, gm::sin(R), std::sin(X));
		Note("cos", X, gm::cos(R), std::cos(X));
		Note("atan", X, gm::atan(R), std::atan(X));
		Note("sinh", X, gm::sinh(R), std::sinh(X));
		Note("tanh", X, gm::tanh(R), std::tanh(X));
		Note("asinh", X, gm::asinh(R), std::asinh(X));
		Note("erf", X, gm::erf(R), std::erf(X));
		Note("erfc", X, gm::erfc(R), std::erfc(X));
		Note("expm1", X, gm::expm1(R), std::expm1(X));
		if (X > 0)
		{
			Note("log", X, gm::log(R), std::log(X));
			Note("log10", X, gm::log10(R), std::log10(X));
			Note("sqrt", X, gm::sqrt(R), std::sqrt(X));
			Note("cbrt", X, gm::cbrt(R), std::cbrt(X));
			Note("lgamma", X, gm::lgamma(R), std::lgamma(X));
			Note("tgamma", X, gm::tgamma(R), std::tgamma(X));
			Note("pow", X, gm::pow(R, Real(1.37)), std::pow(X, 1.37));
		}
		if (X > 0.25)
		{
			Note("acosh", X, gm::acosh(R + Real(1)), std::acosh(X + 1));   // X + 1 is exact in binary only for these
		}
		if (std::fabs(X) < 1)
		{
			Note("asin", X, gm::asin(R), std::asin(X));
			Note("acos", X, gm::acos(R), std::acos(X));
			Note("atanh", X, gm::atanh(R), std::atanh(X));
			Note("log1p", X, gm::log1p(R), std::log1p(X));
		}
		if (X < 0)
		{
			Note("tgamma-", X, gm::tgamma(R), std::tgamma(X));
		}
		Note("atan2", X, gm::atan2(R, Real(-0.7)), std::atan2(X, -0.7));
		Note("fmod", X, gm::fmod(R, Real(0.75)), std::fmod(X, 0.75));
	}
	std::printf("    worst relative error %.3g (%s at %g)\n", Worst, WorstName, WorstX);
	CHECK(Worst < 2e-15);
	CHECK(gm::pow(Real(2), Real(10)) == Real(1024));
	CHECK(gm::sqrt(Real(2)).ToString() == "1.41421356237309505");
	CHECK(gm::exp(Real(1)).ToString() == "2.71828182845904524");
	CHECK(gm::tgamma(Real(5)) == Real(24));
}

TEST(gear_stats_count_the_work)
{
	mcalc::GearContext::ResetStats();
	const Real X = Real(123) * Real(45) + Real(1) / Real(7);
	(void)X;
	const mcalc::GearStats& S = mcalc::GearContext::Stats();
	CHECK(S.Multiplications == 1 && S.Divisions == 1 && S.Additions == 1);
	CHECK(S.Turns >= 9);
}
