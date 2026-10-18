// Properties of the chapter 7, 10, 11 and 12 functions that are not single reference numbers.
#include "CalcSession.h"
#include "TestFramework.h"

#include "GearMath.h"

#include <cstdint>
#include <string>

using mcalc::Real;
namespace gm = mcalc::gm;

namespace
{
Real FieldsEval(const std::string& Expr)
{
	mcalc::Session S;
	const mcalc::Outcome O = S.Execute(Expr);
	if (!O.Ok || !O.HasValue || !O.Val.IsScalar())
	{
		mtest::Fail(__FILE__, __LINE__, Expr + " -> " + O.Text);
		return Real(0);
	}
	return O.Val.Scalar();
}

// Small deterministic generator for random test inputs (xorshift32).
uint32_t FieldsNext(uint32_t* State)
{
	*State ^= *State << 13;
	*State ^= *State >> 17;
	*State ^= *State << 5;
	return *State;
}

std::string RandomBasis(int N, uint32_t* State)
{
	std::string M = "[";
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			M += std::to_string(static_cast<int>(FieldsNext(State) % 81) - 40) + (j + 1 < N ? ", " : "");
		}
		M += i + 1 < N ? "; " : "]";
	}
	return M;
}
} // namespace

TEST(fields_registry_metadata)
{
	const mcalc::Registry& R = mcalc::Registry::Default();
	int Count = 0;
	for (const char* Cat : {"lattices", "category", "procedural", "expansion"})
	{
		for (const mcalc::FunctionInfo* F : R.InCategory(Cat))
		{
			++Count;
			CHECK(!F->Source.empty() && F->Source.compare(0, 2, "ch") == 0);
			CHECK(!F->Summary.empty() && !F->Returns.empty());
			for (const mcalc::Param& P : F->Params)
			{
				CHECK(!P.Unit.empty());
			}
		}
	}
	CHECK(Count >= 50);
}

// phase8 claim 6: |b1| <= 2^((n-1)/2) lambda_1 after LLL, lambda_1 <= |b1|, and the volume is kept.
TEST(fields_lll_bound_random_bases)
{
	uint32_t State = 2463534242u;
	for (int N = 2; N <= 5; ++N)
	{
		for (int Trial = 0; Trial < 6; ++Trial)
		{
			const std::string B = RandomBasis(N, &State);
			if (gm::fabs(FieldsEval("det(" + B + ")")) < Real(1))
			{
				continue;
			}
			std::string B1 = "0";
			for (int j = 1; j <= N; ++j)
			{
				B1 += " + L(1," + std::to_string(j) + ")^2";
			}
			const Real First = gm::sqrt(FieldsEval("L = lll(" + B + "); " + B1));
			const Real Lam = FieldsEval("shortvec(" + B + ")");
			const Real Slack = Real(1) + Real(1e-12);
			CHECK(Lam <= First * Slack);
			CHECK(First <= FieldsEval("lllbound(" + std::to_string(N) + ")") * Lam * Slack);
			const Real Det = gm::fabs(FieldsEval("det(" + B + ")"));
			CHECK_NEAR(gm::fabs(FieldsEval("det(lll(" + B + "))")), Det, Det * Real(1e-12));
		}
	}
}

// phase8 claim 7: rounding <= bootstrapping, and integer decorrelation never lowers bootstrapping.
TEST(fields_ambiguity_estimator_order)
{
	for (const char* Theta : {"0", "0.4", "0.9", "1.3"})
	{
		for (const char* B : {"0.05", "0.15", "0.3"})
		{
			const std::string Q = std::string("ellcov(0.6, ") + B + ", " + Theta + ")";
			const Real Pr = FieldsEval("pround(" + Q + ")"), Pb = FieldsEval("pboot(" + Q + ")");
			const Real Pz = FieldsEval("Q = " + Q + "; Z = zdecor(Q); pboot(Z'*Q*Z)");
			CHECK(Pr <= Pb + Real(1e-12));
			CHECK(Pb <= Pz + Real(1e-12));
			CHECK(Pz <= Real(1));
		}
	}
	CHECK_NEAR(FieldsEval("pround(diag([0.04, 0.09]))"), FieldsEval("pboot(diag([0.04, 0.09]))"), 1e-12);
}

// phase12 claim 4: Lloyd's energy never increases.
TEST(fields_lloyd_energy_monotone)
{
	Real Prev = FieldsEval("lloyd(0.7, 0).e");
	for (int K = 1; K <= 12; ++K)
	{
		const Real E = FieldsEval("lloyd(0.7, " + std::to_string(K) + ").e");
		CHECK(E <= Prev + Real(1e-15));
		Prev = E;
	}
}

// phase11 claim 1 and ch12 claim 6: kernels compose associatively; Sinkhorn plans are doubly stochastic.
TEST(fields_kernels_and_transport)
{
	const std::string F = "ringkern(6, 0.1)", G = "ringkern(6, 0.3)", H = "ringkern(6, 0.45)";
	for (int i = 1; i <= 6; ++i)
	{
		const std::string Ix = "(" + std::to_string(i) + ",2)";
		CHECK_NEAR(FieldsEval("kcomp(" + H + ", kcomp(" + G + ", " + F + "))" + Ix), FieldsEval("kcomp(kcomp(" + H + ", " + G + "), " + F + ")" + Ix), 1e-15);
	}
	const std::string C = "[0.3, 1.2, 2.0; 1.1, 0.2, 0.9; 2.2, 0.8, 0.1]";
	for (int i = 1; i <= 3; ++i)
	{
		const std::string Row = "P(" + std::to_string(i) + ",1) + P(" + std::to_string(i) + ",2) + P(" + std::to_string(i) + ",3)";
		CHECK_NEAR(FieldsEval("P = sinkhorn(" + C + ", 0.05); " + Row), 1.0, 1e-12);
	}
	CHECK(FieldsEval("sinkcost(" + C + ", 0.05)") >= FieldsEval("otexact(" + C + ")") - Real(1e-12));
	CHECK_NEAR(FieldsEval("sinkcost(" + C + ", 0.005)"), FieldsEval("otexact(" + C + ")"), 1e-6);
}

// ch12 claim 2: the switch point lies on the curve x = -v|v|/2 and T = ts + |vs|.
TEST(fields_bangbang_switch_on_curve)
{
	for (const char* X : {"-2", "-0.4", "0.7", "2.5"})
	{
		for (const char* V : {"-1.5", "0.3", "2"})
		{
			const std::string B = std::string("b = bangbang(") + X + ", " + V + "); ";
			CHECK_NEAR(FieldsEval(B + "b.xs + b.vs*abs(b.vs)/2"), 0.0, 1e-12);
			CHECK_NEAR(FieldsEval(B + "b.ts + abs(b.vs) - b.T"), 0.0, 1e-12);
		}
	}
}
