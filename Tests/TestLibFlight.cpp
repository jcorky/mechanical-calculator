// Properties of the flight library that are not single numbers: metadata, inverse pairs, limits and the game's
// "par is a theorem" bounds (no attempt may beat par).
#include "CalcSession.h"
#include "TestFramework.h"

#include <string>

using mcalc::Real;

namespace flight_tests
{
Real Eval(mcalc::Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue)
	{
		std::printf("    '%s' -> %s\n", Line.c_str(), O.Text.c_str());
		return Real::NaN();
	}
	return O.Val.Scalar();
}

std::string Num(const Real& X)
{
	return X.ToString(18);
}
} // namespace flight_tests

TEST(flight_registry_metadata)
{
	const mcalc::Registry& R = mcalc::Registry::Default();
	int Count = 0;
	for (const char* Cat : {"aero", "atmos", "airframe", "game"})
	{
		for (const mcalc::FunctionInfo* F : R.InCategory(Cat))
		{
			++Count;
			CHECK(!F->Source.empty() && !F->Summary.empty() && !F->Returns.empty());
			for (const mcalc::Param& P : F->Params)
			{
				CHECK(!P.Unit.empty());
			}
		}
	}
	CHECK(Count == 57);
}

TEST(flight_area_mach_inverse_on_both_branches)
{
	mcalc::Session S;
	for (const Real& M : {Real(0.2), Real(0.5), Real(0.9), Real(1.2), Real(2.0), Real(3.5)})
	{
		const std::string Line = "machfromarea(arearatio(" + flight_tests::Num(M) + "), " + (M > 1 ? "1" : "0") + ")";
		CHECK_NEAR(flight_tests::Eval(S, Line), M, Real(1e-12) * M);
	}
	CHECK_NEAR(flight_tests::Eval(S, "arearatio(1)"), 1, 1e-15);
}

TEST(flight_bredt_tends_to_the_exact_tube)
{
	mcalc::Session S;
	Real Last(1);
	for (const Real& T : {Real(0.02), Real(0.01), Real(0.005), Real(0.001)})
	{
		const std::string Args = "(1000, 0.1, " + flight_tests::Num(T) + ")";
		const Real Gap = mcalc::gm::fabs(flight_tests::Eval(S, "bredt" + Args + ".bredt") / flight_tests::Eval(S, "bredt" + Args + ".exact") - 1);
		CHECK(Gap < Real(0.6) * T / Real(0.1) && Gap < Last);
		Last = Gap;
	}
}

TEST(flight_lifting_line_bounds)
{
	mcalc::Session S;
	for (const Real& Lam : {Real(0.1), Real(0.3), Real(0.5), Real(1.0)})
	{
		const Real E = flight_tests::Eval(S, "spaneff(8, " + flight_tests::Num(Lam) + ").e");
		CHECK(E > Real(0.9) && E < 1);
	}
	CHECK_NEAR(flight_tests::Eval(S, "spaneff_ell(6).e"), 1, 1e-12);
}

TEST(flight_par_is_never_beaten)
{
	mcalc::Session S;
	const Real GlidePar = flight_tests::Eval(S, "glidepar()");
	Real Best(0);
	for (int i = 0; i <= 120; ++i)
	{
		const Real Score = flight_tests::Eval(S, "glidegame(" + flight_tests::Num(Real(0.2) + Real(i) / 100) + ").score");
		CHECK(Score <= GlidePar * Real(1.000000001));
		Best = mcalc::gm::max(Best, Score);
	}
	CHECK(Best > GlidePar * Real(0.9999));

	const Real PitchPar = flight_tests::Eval(S, "pitchpar()");
	for (int Kp : {0, 4, 10})
	{
		for (int Ki : {0, 4, 10})
		{
			const Real J = flight_tests::Eval(S, "pitchgame(" + std::to_string(Kp) + ", " + std::to_string(Ki) + ", 0.8, 0.1).J");
			CHECK(!(J < PitchPar * Real(0.999999999)));
		}
	}

	const Real RadarPar = flight_tests::Eval(S, "radarpar().range");
	int Admissible = 0;
	for (int P : {25, 50, 100, 200, 400})
	{
		for (int G : {28, 30, 32, 34, 36})
		{
			for (int T : {2, 4, 6, 8, 10, 12})
			{
				S.Execute("r = radargame(" + std::to_string(P) + ", " + std::to_string(G) + ", " + std::to_string(T) + ")");
				if (flight_tests::Eval(S, "r.ok") == 1)
				{
					++Admissible;
					CHECK(flight_tests::Eval(S, "r.range") <= RadarPar);
				}
			}
		}
	}
	CHECK(Admissible == 120);
}
