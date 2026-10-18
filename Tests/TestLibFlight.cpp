// Properties of the flight library that are not single numbers: metadata, inverse pairs, limits and the game's
// "par is a theorem" bounds (no attempt may beat par).
#include "CalcSession.h"
#include "TestFramework.h"

#include <cmath>
#include <string>

namespace flight_tests
{
double Eval(mcalc::Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue)
	{
		std::printf("    '%s' -> %s\n", Line.c_str(), O.Text.c_str());
		return std::nan("");
	}
	return O.Val.Scalar();
}

std::string Num(double X)
{
	char Buf[64];
	std::snprintf(Buf, sizeof(Buf), "%.17g", X);
	return Buf;
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
	for (double M : {0.2, 0.5, 0.9, 1.2, 2.0, 3.5})
	{
		const std::string Line = "machfromarea(arearatio(" + flight_tests::Num(M) + "), " + (M > 1 ? "1" : "0") + ")";
		CHECK_NEAR(flight_tests::Eval(S, Line), M, 1e-12 * M);
	}
	CHECK_NEAR(flight_tests::Eval(S, "arearatio(1)"), 1, 1e-15);
}

TEST(flight_bredt_tends_to_the_exact_tube)
{
	mcalc::Session S;
	double Last = 1;
	for (double T : {0.02, 0.01, 0.005, 0.001})
	{
		const std::string Args = "(1000, 0.1, " + flight_tests::Num(T) + ")";
		const double Gap = std::fabs(flight_tests::Eval(S, "bredt" + Args + ".bredt") / flight_tests::Eval(S, "bredt" + Args + ".exact") - 1);
		CHECK(Gap < 0.6 * T / 0.1 && Gap < Last);
		Last = Gap;
	}
}

TEST(flight_lifting_line_bounds)
{
	mcalc::Session S;
	for (double Lam : {0.1, 0.3, 0.5, 1.0})
	{
		const double E = flight_tests::Eval(S, "spaneff(8, " + flight_tests::Num(Lam) + ").e");
		CHECK(E > 0.9 && E < 1);
	}
	CHECK_NEAR(flight_tests::Eval(S, "spaneff_ell(6).e"), 1, 1e-12);
}

TEST(flight_par_is_never_beaten)
{
	mcalc::Session S;
	const double GlidePar = flight_tests::Eval(S, "glidepar()");
	double Best = 0;
	for (int i = 0; i <= 120; ++i)
	{
		const double Score = flight_tests::Eval(S, "glidegame(" + flight_tests::Num(0.2 + 0.01 * i) + ").score");
		CHECK(Score <= GlidePar * (1 + 1e-9));
		Best = std::max(Best, Score);
	}
	CHECK(Best > GlidePar * (1 - 1e-4));

	const double PitchPar = flight_tests::Eval(S, "pitchpar()");
	for (double Kp : {0.0, 4.0, 10.0})
	{
		for (double Ki : {0.0, 4.0, 10.0})
		{
			const double J = flight_tests::Eval(S, "pitchgame(" + flight_tests::Num(Kp) + ", " + flight_tests::Num(Ki) + ", 0.8, 0.1).J");
			CHECK(!(J < PitchPar * (1 - 1e-9)));
		}
	}

	const double RadarPar = flight_tests::Eval(S, "radarpar().range");
	int Admissible = 0;
	for (int P : {25, 50, 100, 200, 400})
	{
		for (int G : {28, 30, 32, 34, 36})
		{
			for (int T : {2, 4, 6, 8, 10, 12})
			{
				const std::string Args = "(" + std::to_string(P) + ", " + std::to_string(G) + ", " + std::to_string(T) + ")";
				if (flight_tests::Eval(S, "radargame" + Args + ".ok") == 1)
				{
					++Admissible;
					CHECK(flight_tests::Eval(S, "radargame" + Args + ".range") <= RadarPar);
				}
			}
		}
	}
	CHECK(Admissible == 120);
}
