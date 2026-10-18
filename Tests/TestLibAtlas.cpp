// Properties of the atlas library (chapters 1, 2, 4, 5, 6) that are not single reference numbers:
// group laws over many inputs, round trips across the whole angle range, duality on random LPs.
#include "CalcSession.h"
#include "TestFramework.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <string>

namespace
{
double AtlasEval(mcalc::Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue || !O.Val.IsScalar())
	{
		std::printf("    '%s' -> %s\n", Line.c_str(), O.Text.c_str());
		return std::nan("");
	}
	return O.Val.Scalar().ToDouble();
}

struct AtlasLcg
{
	uint64_t S = 88172645463325252ull;
	double Uniform()   // in [-1, 1)
	{
		S = S * 6364136223846793005ull + 1442695040888963407ull;
		return static_cast<double>(S >> 11) / 4503599627370496.0 - 1.0;
	}
	std::string Vec(int N)
	{
		std::string Out = "[";
		for (int i = 0; i < N; ++i)
		{
			char Buf[40];
			std::snprintf(Buf, sizeof(Buf), "%s%.17g", i ? ", " : "", Uniform());
			Out += Buf;
		}
		return Out + "]";
	}
	std::string Rows(int R, int C, double Offset)
	{
		std::string Out = "[";
		for (int i = 0; i < R; ++i)
		{
			for (int j = 0; j < C; ++j)
			{
				char Buf[40];
				std::snprintf(Buf, sizeof(Buf), "%s%.17g", j ? ", " : (i ? "; " : ""), Uniform() + Offset);
				Out += Buf;
			}
		}
		return Out + "]";
	}
};
} // namespace

TEST(atlas_quaternion_homomorphism_and_double_cover)
{
	mcalc::Session S;
	AtlasLcg G;
	for (int t = 0; t < 40; ++t)
	{
		S.Execute("p = " + G.Vec(4) + "; q = " + G.Vec(4));
		CHECK_NEAR(AtlasEval(S, "max(abs(quat2dcm(quat_mul(p, q)) - quat2dcm(p)*quat2dcm(q)))"), 0, 1e-13);
		CHECK_NEAR(AtlasEval(S, "max(abs(quat2dcm(-q) - quat2dcm(q)))"), 0, 0);
		CHECK_NEAR(AtlasEval(S, "v = [0.3, -1, 2]; max(abs(quat_rot(p, v) - (quat2dcm(p)*v')'))"), 0, 1e-13);
		// dcm2quat inverts quat2dcm up to the sign of the cover
		CHECK_NEAR(AtlasEval(S, "u = q/norm(q); d = dcm2quat(quat2dcm(q)); min(norm(d - u), norm(d + u))"), 0, 1e-12);
	}
}

TEST(atlas_exp_log_round_trip_whole_range)
{
	mcalc::Session S;
	AtlasLcg G;
	const double Angles[] = {0.0, 1e-9, 1e-4, 0.5, 1.5, 2.0, 2.5, 3.0, 3.14, 3.14159265, 3.141592653589793};
	for (double T : Angles)
	{
		for (int k = 0; k < 6; ++k)
		{
			char Buf[64];
			std::snprintf(Buf, sizeof(Buf), "t = %.17g", T);
			S.Execute(Buf);
			S.Execute("n = " + G.Vec(3) + "; n = n/norm(n)");
			CHECK_NEAR(AtlasEval(S, "rot_angle(rodrigues(n, t))"), T, 1e-9);
			// Away from pi the log is unique; at pi either sign of the axis is the same rotation.
			CHECK_NEAR(AtlasEval(S, "min(norm(log_so3(rodrigues(n, t)) - t*n), norm(log_so3(rodrigues(n, t)) + t*n))"), 0, 2e-8);
			CHECK_NEAR(AtlasEval(S, "max(abs(exp_so3(log_so3(rodrigues(n, t))) - rodrigues(n, t)))"), 0, 1e-12);
			CHECK_NEAR(AtlasEval(S, "max(abs(rodrigues(n, t) - expm(t*hat(n))))"), 0, 1e-12);
		}
	}
}

TEST(atlas_wahba_methods_agree)
{
	mcalc::Session S;
	AtlasLcg G;
	for (int t = 0; t < 20; ++t)
	{
		S.Execute("A = quat2dcm(" + G.Vec(4) + "); r1 = " + G.Vec(3) + "; r2 = " + G.Vec(3) + "; r3 = " + G.Vec(3) + "; w = abs(" + G.Vec(3) + ") + 0.1");
		S.Execute("b1 = (A*r1')' + 0.01*" + G.Vec(3) + "; b2 = (A*r2')' + 0.01*" + G.Vec(3) + "; b3 = (A*r3')' + 0.01*" + G.Vec(3));
		S.Execute("r = [r1; r2; r3]; b = [b1; b2; b3]");
		CHECK_NEAR(AtlasEval(S, "max(abs(quat2dcm(qmethod(r, b, w)) - wahba_svd(r, b, w)))"), 0, 1e-9);
		CHECK_NEAR(AtlasEval(S, "det(wahba_svd(r, b, w))"), 1, 1e-12);
		// The optimum beats TRIAD and a perturbed attitude.
		const double Opt = AtlasEval(S, "wahba_loss(wahba_svd(r, b, w), r, b, w)");
		CHECK(Opt <= AtlasEval(S, "wahba_loss(triad(r1, r2, b1, b2), r, b, w)") + 1e-15);
		CHECK(Opt < AtlasEval(S, "wahba_loss(wahba_svd(r, b, w)*rodrigues([1, 2, 3], 0.01), r, b, w)"));
	}
}

TEST(atlas_lp_strong_duality_random)
{
	mcalc::Session S;
	AtlasLcg G;
	for (int t = 0; t < 30; ++t)
	{
		// Bounded and feasible: positive constraint rows, positive b, minimize a mostly negative c.
		S.Execute("A = abs(" + G.Rows(6, 4, 0.0) + ") + 0.1; b = abs(" + G.Vec(6) + ") + 0.5; c = " + G.Vec(4) + " - 0.3");
		CHECK_NEAR(AtlasEval(S, "k = lp_kkt(c, A, b, lp(c, A, b), lpdual(c, A, b)); k.gap"), 0, 1e-10);
		CHECK_NEAR(AtlasEval(S, "k.primal"), 0, 1e-10);
		CHECK_NEAR(AtlasEval(S, "k.dual"), 0, 1e-10);
		CHECK_NEAR(AtlasEval(S, "k.slack"), 0, 1e-10);
		CHECK_NEAR(AtlasEval(S, "lpval(c, A, b) + dot(b, lpdual(c, A, b))"), 0, 1e-10);
	}
	// With b < 0 rows (phase 1 needed): x1 + x2 >= 1 written as -x1 - x2 <= -1.
	CHECK_NEAR(AtlasEval(S, "lpval([1, 2], [-1, -1; 1, 0], [-1, 3])"), 1, 1e-12);
	CHECK_NEAR(AtlasEval(S, "lpdual([1, 2], [-1, -1; 1, 0], [-1, 3])(1)"), 1, 1e-12);
}

TEST(atlas_lasso_optimality_conditions)
{
	mcalc::Session S;
	AtlasLcg G;
	const double Lam = 0.3;
	for (int t = 0; t < 10; ++t)
	{
		S.Execute("A = " + G.Rows(12, 5, 0.0) + "; b = " + G.Vec(12) + "; lam = 0.3");
		S.Execute("x = lasso(A, b, lam); g = A'*(b' - A*x)");
		mcalc::Value X, Gr;
		CHECK(S.GetVariable("x", &X) && S.GetVariable("g", &Gr) && X.M.Count() == 5 && Gr.M.Count() == 5);
		for (int j = 0; j < X.M.Count() && j < Gr.M.Count(); ++j)
		{
			// subgradient condition: |g_j| <= lambda, and g_j = lambda sign(x_j) on the support
			const double Xj = X.M.V[static_cast<size_t>(j)].ToDouble(), Gj = Gr.M.V[static_cast<size_t>(j)].ToDouble();
			CHECK(std::fabs(Gj) <= Lam + 1e-9);
			if (Xj != 0.0)
			{
				CHECK_NEAR(Gj, Xj > 0 ? Lam : -Lam, 1e-9);
			}
		}
		// x is a fixed point of the proximal-gradient map
		CHECK_NEAR(AtlasEval(S, "s = 1/norm(A, 2)^2; max(abs(soft(x + s*g, s*lam) - x))"), 0, 1e-12);
	}
}

TEST(atlas_topology_invariants)
{
	mcalc::Session S;
	const char* Ks[] = {"1", "2", "3", "5", "8"};
	for (const char* K : Ks)
	{
		CHECK_NEAR(AtlasEval(S, std::string("nav_index(") + K + ").total"), -2, 0);   // chi of the disk with 3 holes
	}
	CHECK_NEAR(AtlasEval(S, "nav_index(8).minus"), 3, 0);
	for (double Th = 0; Th < 1.5; Th += 0.1)
	{
		char Buf[96];
		std::snprintf(Buf, sizeof(Buf), "euler_det(%.17g, 0.7)*cos(%.17g)", Th, Th);
		CHECK_NEAR(AtlasEval(S, Buf), 1, 1e-12);
		std::snprintf(Buf, sizeof(Buf), "ea = dcm2eul(eul2dcm(0.4, %.17g, -1.1)); ea.theta", Th);
		CHECK_NEAR(AtlasEval(S, Buf), Th, 1e-12);
	}
}

TEST(atlas_measure_and_filtering_identities)
{
	mcalc::Session S;
	// Causal Wiener (steady Kalman) is never better than the non-causal one; equal only for a white signal.
	for (double A = 0.0; A < 0.99; A += 0.1)
	{
		char Buf[64];
		std::snprintf(Buf, sizeof(Buf), "w = wiener(%.17g, 0.2, 0.8); w.price", A);
		CHECK(AtlasEval(S, Buf) >= (A < 0.05 ? -1e-12 : 1e-6));   // a = 0: white signal, causal = non-causal
	}
	// Gibbs: no distribution beats the Boltzmann one.
	CHECK(AtlasEval(S, "E = [0.3, -1, 2, 0.7, 1.1]; logz(E) - gibbs_f([0.2, 0.3, 0.1, 0.2, 0.2], E)") > 0);
	CHECK_NEAR(AtlasEval(S, "logz(E) - gibbs_f(boltz(E), E)"), 0, 1e-14);
	// ESS is between 1 and N.
	CHECK_NEAR(AtlasEval(S, "ess([0, 0, 1, 0])"), 1, 0);
	CHECK_NEAR(AtlasEval(S, "ess([2, 2, 2])"), 3, 1e-15);
	// Halton points lie in the unit cube and the first coordinate is van der Corput base 2.
	CHECK_NEAR(AtlasEval(S, "h = halton(64, 3); max(h)"), 63.0 / 64.0, 1e-15);
	CHECK_NEAR(AtlasEval(S, "h(13, 1) - vdc(13)"), 0, 0);
}

// Each function at its bench defaults must answer within a keypress budget on the gears (the brief:
// about 0.5 s at defaults, 2 s at the input caps). Slow ones are printed.
TEST(atlas_timing_on_the_gears)
{
	struct Probe
	{
		const char* Line;
		double LimitMs;
	};
	const Probe Probes[] = {
		{"nav_index()", 500}, {"unwind()", 500}, {"qmc_err()", 500}, {"bode_int()", 500}, {"dualgap()", 500},
		{"postmean([0.4, -0.2, 9.5, 0.9, 0.1])", 500}, {"comm_angle()", 500}, {"euler_cond()", 500}, {"is_w2(1.5)", 500},
		{"wiener()", 500}, {"kf_ss()", 500}, {"navfn(0.3, -0.5)", 500},
		{"lasso(reshape(linspace(-1, 1.3, 60), 12, 5) + [eye(5); zeros(7, 5)], linspace(0, 1, 12), 0.1)", 500},
		{"lp(-ones(1, 20), ones(20, 20) + eye(20), ones(20, 1))", 500},
		{"nav_index(2, 0, 0.6, 0.25, 400)", 2000}, {"qmc_err(100000, 1)", 2000}, {"qmc_err(12500, 8)", 2000}, {"x = halton(20000, 4)", 2000},
		{"lq_cost(eye(2), [0; 1], eye(2), 1, eye(2), [1; 0], 2000)", 2000},
	};
	for (const Probe& P : Probes)
	{
		mcalc::Session S;
		const auto T0 = std::chrono::steady_clock::now();
		const mcalc::Outcome O = S.Execute(P.Line);
		const double Ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - T0).count();
		CHECK(O.Ok);
		if (Ms > 50)
		{
			std::printf("    timing: %s took %.1f ms\n", P.Line, Ms);
		}
		CHECK(Ms < P.LimitMs);
	}
}
