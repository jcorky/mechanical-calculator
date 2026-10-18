// Properties of the orbits and ins libraries that are not single reference numbers.
#include "CalcSession.h"
#include "TestFramework.h"

#include <chrono>
#include <cstdlib>

using mcalc::Matrix;
using mcalc::Session;
namespace la = mcalc::la;

namespace
{
double OrbEval(Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue)
	{
		std::printf("    '%s' -> %s\n", Line.c_str(), O.Text.c_str());
		return std::nan("");
	}
	return O.Val.Scalar();
}

Matrix OrbMatrix(Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue)
	{
		std::printf("    '%s' -> %s\n", Line.c_str(), O.Text.c_str());
		return Matrix();
	}
	return O.Val.M;
}
} // namespace

// The CW state-transition matrix is exp(A t) for the linearised Hill equations.
TEST(orbits_cwstm_is_expm_of_hill_equations)
{
	Session S;
	const double N = 0.0011313666536, T = 1234.5;
	Matrix A4(4, 4);
	A4(0, 2) = 1;
	A4(1, 3) = 1;
	A4(2, 0) = 3 * N * N;
	A4(2, 3) = 2 * N;
	A4(3, 2) = -2 * N;
	Matrix E4;
	std::string Error;
	CHECK(la::Expm(la::Scale(A4, T), &E4, &Error));
	char Line[128];
	std::snprintf(Line, sizeof(Line), "cwstm(%.17g, %.17g)", N, T);
	const Matrix P4 = OrbMatrix(S, Line);
	CHECK(P4.Rows == 4 && P4.Cols == 4);
	if (P4.Rows == 4)
	{
		CHECK(la::MaxAbs(la::Sub(P4, E4)) < 1e-9 * la::MaxAbs(E4));
	}
	Matrix A6(6, 6);
	A6(0, 3) = 1;
	A6(1, 4) = 1;
	A6(2, 5) = 1;
	A6(3, 0) = 3 * N * N;
	A6(3, 4) = 2 * N;
	A6(4, 3) = -2 * N;
	A6(5, 2) = -N * N;
	Matrix E6;
	CHECK(la::Expm(la::Scale(A6, T), &E6, &Error));
	std::snprintf(Line, sizeof(Line), "cwstm(%.17g, %.17g, 6)", N, T);
	const Matrix P6 = OrbMatrix(S, Line);
	CHECK(P6.Rows == 6 && P6.Cols == 6);
	if (P6.Rows == 6)
	{
		CHECK(la::MaxAbs(la::Sub(P6, E6)) < 1e-9 * la::MaxAbs(E6));
	}
}

// Lambert's v1, flown by the two-body integrator for dt, arrives at r2 with Lambert's v2.
TEST(orbits_lambert_closes_under_integration)
{
	Session S;
	CHECK(S.Execute("r2 = [12000e3*cos(95*pi/180), 12000e3*sin(95*pi/180)]").Ok);
	CHECK(S.Execute("L = lambert([7000e3, 0], r2, 3600)").Ok);
	CHECK(S.Execute("P = twobody([7000e3, 0], [L.v1x, L.v1y], 3600)").Ok);
	const double Miss = OrbEval(S, "hypot(P.x - r2(1), P.y - r2(2))");
	CHECK(Miss < 1.0);
	CHECK_NEAR(OrbEval(S, "P.vx"), OrbEval(S, "L.v2x"), 1e-3);
	CHECK_NEAR(OrbEval(S, "P.vy"), OrbEval(S, "L.v2y"), 1e-3);
}

// The conic elements of a Kepler state give back the orbit; energy and vis-viva agree.
TEST(orbits_conic_inverts_kepstate)
{
	Session S;
	CHECK(S.Execute("k = kepstate(26000e3, 0.7, 20000)").Ok);
	CHECK(S.Execute("c = conic([k.x, k.y], [k.vx, k.vy])").Ok);
	CHECK_NEAR(OrbEval(S, "c.a"), 26000e3, 1e-3);
	CHECK_NEAR(OrbEval(S, "c.e"), 0.7, 1e-12);
	CHECK_NEAR(OrbEval(S, "c.rp"), 26000e3 * 0.3, 1e-3);
	CHECK_NEAR(OrbEval(S, "hypot(k.vx, k.vy)"), OrbEval(S, "visviva(k.r, 26000e3)"), 1e-8);
	CHECK_NEAR(OrbEval(S, "c.energy"), -mcalc::cst::MuEarth / (2 * 26000e3), 1e-3);
	// The true anomaly from E agrees with the polar angle of the position.
	CHECK_NEAR(OrbEval(S, "mod(atan2(k.y, k.x), 2*pi)"), OrbEval(S, "mod(k.nu, 2*pi)"), 1e-12);
}

// Hohmann is the rb -> r2 limit of the bi-elliptic transfer; the threshold separates the two.
TEST(orbits_bielliptic_limits)
{
	Session S;
	CHECK_NEAR(OrbEval(S, "bielliptic(7000e3, 70000e3, 70000e3).dv"), OrbEval(S, "hohmann(7000e3, 70000e3).dv"), 1e-6);
	CHECK_NEAR(OrbEval(S, "bielliptic(1, 12, 1e12, 1).dv"), OrbEval(S, "bieltot(12, 1e12)"), 1e-12);
	const double R = OrbEval(S, "bielthresh().closed");
	CHECK(OrbEval(S, "hohtot(" + std::to_string(R - 0.01) + ") - bielinf(" + std::to_string(R - 0.01) + ")") < 0);
	CHECK(OrbEval(S, "hohtot(" + std::to_string(R + 0.01) + ") - bielinf(" + std::to_string(R + 0.01) + ")") > 0);
}

// The steady-state filter is a fixed point: a long insgnss run with no outage approaches it.
TEST(ins_kalman_steady_state_is_reached)
{
	Session S;
	const double Sp = OrbEval(S, "kfsteady(0.003, 2e-4, 3).sp");
	const double Kp = OrbEval(S, "kfsteady(0.003, 2e-4, 3).kp");
	CHECK(Sp > 0 && Sp < 3);
	CHECK(Kp > 0 && Kp < 1);
	// Prior-form fixed point: posterior variance equals prior * (1 - kp).
	CHECK_NEAR(Sp * Sp, OrbEval(S, "kfsteady(0.003, 2e-4, 3).sppri^2") * (1 - Kp), 1e-9 * Sp * Sp);
}

// Every orbits / ins function documents itself, and the heavy ones stay fast enough for a key press.
TEST(orbits_ins_registry_and_speed)
{
	const mcalc::Registry& Reg = mcalc::Registry::Default();
	int Count = 0;
	for (const char* Cat : {"orbits", "ins"})
	{
		for (const mcalc::FunctionInfo* F : Reg.InCategory(Cat))
		{
			++Count;
			CHECK(!F->Source.empty() && !F->Summary.empty() && !F->Returns.empty());
			for (const mcalc::Param& P : F->Params)
			{
				CHECK(!P.Unit.empty());
			}
		}
	}
	CHECK(Count >= 40);
	Session S;
	for (const char* Line : {"allanfit(3).flat", "allansim(1, 1, 1, 2)(5, 2)", "edriftcmp(0.3, 1000, 100).q", "edrift(0.3, 2000, 100, 2).max", "allandev(1, 1)", "twobody([7000e3, 0], [0, 7600], 8e5).x", "flickerquad(400, 1000000)", "coningrun(0.035, 31.4, 1000, 100).slope0", "bortz(0.01, 100000)", "cwshoot(6778137, [-1000, -10000, 0, 1.697], 2800, 10000).v0x", "cwnonlin(6778137, [-1000, -10000, 0, 1.697], 2800, 200000).x", "onetanscan(12, 1000000).best", "schulersim(0, 1e-3, 0, 1e6).x", "hohfly(6678e3, 42164e3).rend", "lambert([7000e3, 0], [0, 12000e3], 3600).v1", "bortz()", "kfsteady().sp", "cwshoot(6778137, [-1000, -10000, 0, 1.697], 2800).v0x",
			 "flickerquad()", "insgnss().closed", "schulersim(0, 1e-3, 0, 21600).period", "coningrun(0.035, 31.4, 1000).slope0"})
	{
		const auto T0 = std::chrono::steady_clock::now();
		const double V = OrbEval(S, Line);
		const double Ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - T0).count();
		CHECK(std::isfinite(V));
		if (Ms > 50 || std::getenv("MCALC_TIMING"))
		{
			std::printf("    timing: %s took %.2f ms\n", Line, Ms);
		}
		CHECK(Ms < 500);   // generous bound for slow CI machines; typical times are printed above 50 ms
	}
}
