// Properties of the control library that are not single numbers: metadata, Riccati residuals,
// MPC = LQR for every horizon, the two Hamiltonian halves, duality and the separation principle.
#include "CalcSession.h"
#include "TestFramework.h"

#include <string>

using mcalc::Matrix;
using mcalc::Session;
namespace la = mcalc::la;

namespace
{
Matrix EvalM(Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue)
	{
		std::printf("    '%s' -> %s\n", Line.c_str(), O.Text.c_str());
		return Matrix();
	}
	return O.Val.M;
}

double MaxDiff(const Matrix& A, const Matrix& B)
{
	if (A.Rows != B.Rows || A.Cols != B.Cols || A.IsEmpty())
	{
		return 1e300;
	}
	return la::MaxAbs(la::Sub(A, B));
}
} // namespace

TEST(control_library_metadata)
{
	const mcalc::Registry& R = mcalc::Registry::Default();
	int Count = 0;
	for (const char* Cat : {"riccati", "control", "dynamics"})
	{
		for (const mcalc::FunctionInfo* F : R.InCategory(Cat))
		{
			++Count;
			CHECK(!F->Source.empty());
			CHECK(!F->Summary.empty());
			CHECK(!F->Returns.empty());
			for (const mcalc::Param& P : F->Params)
			{
				CHECK(!P.Unit.empty());
			}
		}
	}
	CHECK(Count >= 40);
	for (const char* Name : {"dare", "care", "dlyap", "lyap", "dlqr", "lqr", "kalman", "c2d", "mpcgain", "clcost", "gm"})
	{
		CHECK(R.Find(Name) != nullptr);
	}
}

TEST(control_riccati_residuals)
{
	Session S;
	S.Execute("A = [1.1, 0.3, 0; 0, 0.8, 0.5; 0.2, 0, 0.95]; B = [0, 1; 1, 0; 0.5, 0.5]; Q = [2, 0.1, 0; 0.1, 1, 0; 0, 0, 0.5]; R = [1, 0.2; 0.2, 0.5]");
	const Matrix A = EvalM(S, "A"), B = EvalM(S, "B"), Q = EvalM(S, "Q"), Rm = EvalM(S, "R");
	const Matrix P = EvalM(S, "dare(A, B, Q, R)"), K = EvalM(S, "dlqr(A, B, Q, R)");
	const Matrix BtP = la::Mul(la::Transpose(B), P);
	Matrix G;
	CHECK(la::Solve(la::Add(Rm, la::Mul(BtP, B)), la::Mul(BtP, A), &G, nullptr));
	const Matrix Next = la::Add(Q, la::Sub(la::Mul(la::Mul(la::Transpose(A), P), A), la::Mul(la::Transpose(la::Mul(BtP, A)), G)));
	CHECK(MaxDiff(Next, P) < 1e-10 * la::MaxAbs(P));
	CHECK(MaxDiff(G, K) < 1e-10);
	CHECK(la::SpectralRadius(la::Sub(A, la::Mul(B, K))) < 1);
	// Continuous: A'P + PA - P B R^-1 B' P + Q = 0 and A - BK Hurwitz.
	S.Execute("Ac = [0, 1, 0; 0, 0, 1; -1, -2, 0.5]");
	const Matrix Ac = EvalM(S, "Ac"), Pc = EvalM(S, "care(Ac, B, Q, R)"), Kc = EvalM(S, "lqr(Ac, B, Q, R)");
	const Matrix Res = la::Add(la::Sub(la::Add(la::Mul(la::Transpose(Ac), Pc), la::Mul(Pc, Ac)), la::Mul(la::Mul(Pc, B), Kc)), Q);
	CHECK(la::MaxAbs(Res) < 1e-9 * la::MaxAbs(Pc));
	CHECK_NEAR(EvalM(S, "eig(Ac - B*lqr(Ac, B, Q, R))(1, 1)").ScalarValue() < 0 ? 1 : 0, 1, 0);
	// Lyapunov solvers against their equations.
	const Matrix X = EvalM(S, "dlyap(A/2, Q)");
	CHECK(MaxDiff(la::Add(la::Mul(la::Mul(la::Scale(A, 0.5), X), la::Transpose(la::Scale(A, 0.5))), Q), X) < 1e-12);
	const Matrix Y = EvalM(S, "lyap(Ac - B*lqr(Ac, B, Q, R), Q)");
	const Matrix F = la::Sub(Ac, la::Mul(B, Kc));
	CHECK(la::MaxAbs(la::Add(la::Add(la::Mul(F, Y), la::Mul(Y, la::Transpose(F))), Q)) < 1e-10);
}

TEST(control_mpc_with_terminal_P_is_lqr_for_every_horizon)
{
	Session S;
	S.Execute("A = [1, 0.2; 0, 1]; B = [0.02; 0.2]; Q = [1, 0; 0, 0.1]; R = 0.05; P = dare(A, B, Q, R); K = dlqr(A, B, Q, R)");
	const Matrix K = EvalM(S, "K");
	for (int N = 1; N <= 30; N += 7)
	{
		CHECK(MaxDiff(EvalM(S, "mpcgain(A, B, Q, R, P, " + std::to_string(N) + ")"), K) < 1e-9);
	}
	// Without terminal cost the N-step gain approaches LQR's as N grows.
	const double Far = MaxDiff(EvalM(S, "mpcgain(A, B, Q, R, zeros(2), 3)"), K), Near = MaxDiff(EvalM(S, "mpcgain(A, B, Q, R, zeros(2), 200)"), K);
	CHECK(Near < 1e-9 && Far > 1e-3);
	// The closed-loop cost of the LQR gain is x0'Px0.
	CHECK_NEAR(EvalM(S, "clcost(A, B, Q, R, K, [3; -1])").ScalarValue(), EvalM(S, "[3, -1]*P*[3; -1]").ScalarValue(), 1e-9);
}

TEST(control_duality_separation_and_hamiltonian_halves)
{
	Session S;
	S.Execute("A = [1, 0.5; 0, 1]; C = [1, 0]; Q = [0.04, 0.12; 0.12, 0.5]; R = 0.3; M = kalman(A, C, Q, R)");
	// Kalman-LQR duality (ch03 Theorem 5): the dual regulator gain is (A M)'.
	CHECK(MaxDiff(EvalM(S, "dlqr(A', C', Q, R)"), EvalM(S, "(A*M)'")) < 1e-10);
	// Separation: the LQG loop's spectral radius is the larger of controller and estimator.
	S.Execute("B = [0.125; 0.5]; K = dlqr(A, B, eye(2), 1)");
	CHECK_NEAR(EvalM(S, "rho(lqgcl(A, B, C, K, M))").ScalarValue(), EvalM(S, "max(rho(A - B*K), rho(A - A*M*C))").ScalarValue(), 1e-9);
	// The two halves of the Hamiltonian (ch03 continuous guise): both solve the ARE, P > 0 and P < 0.
	S.Execute("Ao = [0, 1; -2, -0.3]; Qo = [0.1, 0; 0, 0.5]");
	const Matrix Pp = EvalM(S, "hamare(Ao, C, Qo, 0.2, 1)"), Pm = EvalM(S, "hamare(Ao, C, Qo, 0.2, -1)");
	std::vector<double> W;
	Matrix V;
	CHECK(la::SymmetricEigen(Pp, &W, &V, nullptr) && W[0] > 0);
	CHECK(la::SymmetricEigen(Pm, &W, &V, nullptr) && W[1] < 0);
	// The stabilizing half agrees with care() on the dual pair.
	CHECK(MaxDiff(Pp, EvalM(S, "care(Ao', C', Qo, 0.2)")) < 1e-9);
	// Symplectic Z: eigenvalues in reciprocal pairs, so det Z = 1 and rho(Z) * min|pole| = 1.
	CHECK_NEAR(EvalM(S, "det(sympz(A, C, Q, R))").ScalarValue(), 1, 1e-10);
}

TEST(control_scalar_map_readings_agree)
{
	Session S;
	// The LQR reading is the Kalman map with h -> b, r -> rho (atlas Phase 1, claim 5).
	CHECK_NEAR(EvalM(S, "lqrss(0.8, 0.3, 2, 1.5).S").ScalarValue(), EvalM(S, "kfss(0.8, 0.3, 2, 1.5).p").ScalarValue(), 1e-15);
	CHECK_NEAR(EvalM(S, "lqrss(0.8, 0.3, 2, 1.5).S").ScalarValue(), EvalM(S, "dare(0.8, 1.5, 0.3, 2)").ScalarValue(), 1e-12);
	CHECK_NEAR(EvalM(S, "kfss(0.8, 0.3, 2, 1.5).p").ScalarValue(), EvalM(S, "kfcov(0.8, 1.5, 0.3, 2)").ScalarValue(), 1e-12);
	CHECK_NEAR(EvalM(S, "kfss(0.8, 0.3, 2, 1.5).K").ScalarValue(), EvalM(S, "kalman(0.8, 1.5, 0.3, 2)").ScalarValue(), 1e-12);
	// The covariance map is never elliptic (ch03 Theorem 4).
	for (const char* Set : {"a = 0.95; q = 0.2; r = 1; h = 1", "a = 1.8; q = 0; r = 3; h = 0.05", "a = -1.7; q = 2; r = 0.05; h = 2", "a = 0.3; q = 1.5; r = 0.2; h = 1.9"})
	{
		S.Execute(Set);
		CHECK(EvalM(S, "mobclass(a^2*r + q*h^2, q*r, h^2, r).type").ScalarValue() != 4);
		CHECK(EvalM(S, "mobclass(a^2*r + q*h^2, q*r, h^2, r).disc").ScalarValue() >= 0);
	}
}
