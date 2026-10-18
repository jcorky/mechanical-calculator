// Properties of the operations library that are not single numbers: feasibility and duality of transport
// plans, permutations from the assignment, flow conservation, tour inequalities, exact rigid alignment.
#include "CalcSession.h"
#include "TestFramework.h"

#include <algorithm>
#include <cstdint>
#include <numeric>

namespace
{
struct OpsLcg
{
	uint64_t S;
	explicit OpsLcg(uint64_t Seed) : S(Seed) {}
	mcalc::Real U()
	{
		S = S * 6364136223846793005ull + 1442695040888963407ull;
		return mcalc::Real(static_cast<unsigned long long>(S >> 11)) / mcalc::Real(9007199254740992ull);
	}
};

mcalc::Matrix OpsRandom(OpsLcg& G, int R, int C, mcalc::Real Scale)
{
	mcalc::Matrix M(R, C);
	for (mcalc::Real& X : M.V)
	{
		X = Scale * G.U();
	}
	return M;
}

mcalc::Value OpsEval(mcalc::Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	CHECK(O.Ok);
	return O.Val;
}
} // namespace

TEST(ops_transport_plan_is_feasible_and_priced)
{
	OpsLcg G(7);
	for (int Trial = 0; Trial < 20; ++Trial)
	{
		mcalc::Session S;
		const int M = 2 + Trial % 3, N = 3 + Trial % 4;
		mcalc::Matrix Sup(1, M), Dem(1, N);
		mcalc::Real Total = 0.0;
		for (mcalc::Real& X : Sup.V)
		{
			X = mcalc::gm::floor(1 + 20 * G.U());
			Total += X;
		}
		for (int j = 0; j < N - 1; ++j)
		{
			Dem.V[static_cast<size_t>(j)] = mcalc::gm::floor(Total / N);
		}
		Dem.V[static_cast<size_t>(N - 1)] = Total - (N - 1) * mcalc::gm::floor(Total / N);
		mcalc::Matrix C = OpsRandom(G, M, N, 10);
		for (mcalc::Real& X : C.V)
		{
			X = mcalc::gm::round(X);
		}
		S.SetVariable("s", mcalc::Value(Sup));
		S.SetVariable("d", mcalc::Value(Dem));
		S.SetVariable("C", mcalc::Value(C));
		const mcalc::Real Cost = OpsEval(S, "transport(s, d, C)").Scalar();
		const mcalc::Matrix X = OpsEval(S, "tplan(s, d, C)").M;
		const mcalc::Matrix P = OpsEval(S, "tprices(s, d, C)").M;
		mcalc::Real Primal = 0.0, Dual = 0.0;
		for (int i = 0; i < M; ++i)
		{
			mcalc::Real Row = 0.0;
			for (int j = 0; j < N; ++j)
			{
				Row += X(i, j);
				Primal += X(i, j) * C(i, j);
				CHECK(X(i, j) >= -1e-9);
				CHECK(X(i, j) == mcalc::gm::round(X(i, j)));                               // integral (total unimodularity)
				CHECK(C(i, j) - P.V[static_cast<size_t>(i)] - P.V[static_cast<size_t>(M + j)] >= -1e-9);   // dual feasible
			}
			CHECK_NEAR(Row, Sup.V[static_cast<size_t>(i)], 1e-9);
			Dual += Sup.V[static_cast<size_t>(i)] * P.V[static_cast<size_t>(i)];
		}
		for (int j = 0; j < N; ++j)
		{
			mcalc::Real Col = 0.0;
			for (int i = 0; i < M; ++i)
			{
				Col += X(i, j);
			}
			CHECK_NEAR(Col, Dem.V[static_cast<size_t>(j)], 1e-9);
			Dual += Dem.V[static_cast<size_t>(j)] * P.V[static_cast<size_t>(M + j)];
		}
		CHECK_NEAR(Primal, Cost, 1e-9);
		CHECK_NEAR(Dual, Cost, 1e-9);   // zero duality gap
	}
}

TEST(ops_assignment_is_the_best_permutation)
{
	OpsLcg G(11);
	for (int Trial = 0; Trial < 10; ++Trial)
	{
		mcalc::Session S;
		const int N = 3 + Trial % 4;
		const mcalc::Matrix C = OpsRandom(G, N, N, 10);
		S.SetVariable("C", mcalc::Value(C));
		const mcalc::Real Cost = OpsEval(S, "assign(C)").Scalar();
		const mcalc::Matrix P = OpsEval(S, "assignperm(C)").M;
		std::vector<int> Q(static_cast<size_t>(N));
		std::iota(Q.begin(), Q.end(), 0);
		mcalc::Real Brute = 1e300;
		do
		{
			mcalc::Real Sum = 0.0;
			for (int i = 0; i < N; ++i)
			{
				Sum += C(i, Q[static_cast<size_t>(i)]);
			}
			Brute = std::min(Brute, Sum);
		} while (std::next_permutation(Q.begin(), Q.end()));
		CHECK_NEAR(Cost, Brute, 1e-9);
		std::vector<mcalc::Real> Sorted = P.V;
		std::sort(Sorted.begin(), Sorted.end());
		for (int i = 0; i < N; ++i)
		{
			CHECK(Sorted[static_cast<size_t>(i)] == i + 1);
		}
	}
}

TEST(ops_flow_is_conserved_and_equals_cut)
{
	OpsLcg G(3);
	for (int Trial = 0; Trial < 10; ++Trial)
	{
		mcalc::Session S;
		const int N = 5 + Trial % 3;
		mcalc::Matrix C = OpsRandom(G, N, N, 9);
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				C(i, j) = (i == j || G.U() < 0.4) ? 0.0 : mcalc::gm::round(C(i, j));
			}
		}
		S.SetVariable("C", mcalc::Value(C));
		const mcalc::Real Flow = OpsEval(S, "maxflow(C).flow").Scalar();
		const mcalc::Matrix F = OpsEval(S, "flowmat(C)").M;
		for (int v = 0; v < N; ++v)
		{
			mcalc::Real Net = 0.0;
			for (int u = 0; u < N; ++u)
			{
				Net += F(v, u) - F(u, v);
				CHECK(F(v, u) <= C(v, u) + 1e-9 && F(v, u) >= -1e-9);
			}
			CHECK_NEAR(Net, v == 0 ? Flow : (v == N - 1 ? -Flow : 0.0), 1e-9);
		}
		// brute-force minimum cut over the inner nodes
		mcalc::Real Cut = 1e300;
		for (int Mask = 0; Mask < (1 << (N - 2)); ++Mask)
		{
			std::vector<char> In(static_cast<size_t>(N), 0);
			In[0] = 1;
			for (int k = 0; k < N - 2; ++k)
			{
				In[static_cast<size_t>(k + 1)] = static_cast<char>((Mask >> k) & 1);
			}
			mcalc::Real Sum = 0.0;
			for (int u = 0; u < N; ++u)
			{
				for (int v = 0; v < N; ++v)
				{
					if (In[static_cast<size_t>(u)] && !In[static_cast<size_t>(v)])
					{
						Sum += C(u, v);
					}
				}
			}
			Cut = std::min(Cut, Sum);
		}
		CHECK_NEAR(Flow, Cut, 1e-9);
	}
}

TEST(ops_tour_bounds_hold)
{
	OpsLcg G(5);
	for (int Trial = 0; Trial < 15; ++Trial)
	{
		mcalc::Session S;
		S.SetVariable("P", mcalc::Value(OpsRandom(G, 5 + Trial % 5, 2, 1)));
		const mcalc::Real Opt = OpsEval(S, "tsp(distmat(P))").Scalar();
		const mcalc::Real Chr = OpsEval(S, "christofides(distmat(P)).len").Scalar();
		const mcalc::Real Mst = OpsEval(S, "mst(distmat(P))").Scalar();
		const mcalc::Real Lb = OpsEval(S, "hkbound(distmat(P)).lb").Scalar();
		const mcalc::Real One = OpsEval(S, "hkbound(distmat(P)).onetree").Scalar();
		CHECK(Opt <= Chr + 1e-12 && Chr <= 1.5 * Opt + 1e-12);
		CHECK(Mst <= One + 1e-12 && One <= Lb + 1e-12 && Lb <= Opt + 1e-9);
	}
}

TEST(ops_shortest_paths_agree)
{
	OpsLcg G(13);
	for (int Trial = 0; Trial < 10; ++Trial)
	{
		mcalc::Session S;
		const int N = 6 + Trial % 4;
		mcalc::Matrix W = OpsRandom(G, N, N, 5);
		for (mcalc::Real& X : W.V)
		{
			X = G.U() < 0.5 ? X : mcalc::Real::Infinity();
		}
		S.SetVariable("W", mcalc::Value(W));
		const mcalc::Matrix D = OpsEval(S, "dijkstra(W)").M, B = OpsEval(S, "bellman(W)").M, F = OpsEval(S, "floyd(W)").M;
		for (int k = 0; k < N; ++k)
		{
			CHECK_NEAR(D.V[static_cast<size_t>(k)], B.V[static_cast<size_t>(k)], 1e-12);
			CHECK_NEAR(D.V[static_cast<size_t>(k)], F(0, k), 1e-12);
		}
	}
}

TEST(ops_kabsch_recovers_a_rigid_motion)
{
	OpsLcg G(17);
	for (int Trial = 0; Trial < 10; ++Trial)
	{
		mcalc::Session S;
		const int N = 3 + Trial;
		const mcalc::Real Th = -3.0 + 6.0 * G.U(), Tx = G.U() - 0.5, Ty = 2.0 * G.U();
		mcalc::Matrix P = OpsRandom(G, N, 2, 2), Q(N, 2);
		for (int i = 0; i < N; ++i)
		{
			Q(i, 0) = mcalc::gm::cos(Th) * P(i, 0) - mcalc::gm::sin(Th) * P(i, 1) + Tx;
			Q(i, 1) = mcalc::gm::sin(Th) * P(i, 0) + mcalc::gm::cos(Th) * P(i, 1) + Ty;
		}
		S.SetVariable("P", mcalc::Value(P));
		S.SetVariable("Q", mcalc::Value(Q));
		CHECK_NEAR(OpsEval(S, "kabsch(P, Q).th").Scalar(), Th, 1e-9);
		CHECK_NEAR(OpsEval(S, "kabsch(P, Q).tx").Scalar(), Tx, 1e-9);
		CHECK_NEAR(OpsEval(S, "kabsch(P, Q).ty").Scalar(), Ty, 1e-9);
		CHECK_NEAR(OpsEval(S, "kabsch(P, Q).cost").Scalar(), 0.0, 1e-18);
	}
}

TEST(ops_chi2_quantile_inverts_cdf)
{
	mcalc::Session S;
	for (mcalc::Real K : {1.0, 2.0, 3.0, 7.0, 30.0, 300.0})
	{
		for (mcalc::Real P : {0.01, 0.5, 0.95, 0.999})
		{
			S.SetVariable("k", mcalc::Value(K));
			S.SetVariable("p", mcalc::Value(P));
			CHECK_NEAR(OpsEval(S, "chi2cdf(chi2inv(p, k), k)").Scalar(), P, 1e-12);
		}
	}
	// the triangle with unit conductances: one link in parallel with two in series
	CHECK_NEAR(OpsEval(S, "reff([0, 1, 1; 1, 0, 1; 1, 1, 0], 1, 2)").Scalar(), 2.0 / 3.0, 1e-12);
}
