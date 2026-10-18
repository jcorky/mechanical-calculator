// Operations: queues, inventory, network flows, transport, max-plus timetables and routing (Chapter 13 and
// the logistics bench), shortest paths and graph Laplacians (Chapter 9), and the SLAM bench's chi-square
// tools, 2-D rigid alignment and odometry drift. Graphs are square matrices (cost or capacity; inf = no
// edge in cost matrices), point sets are n x 2 matrices, and node numbers are 1-based.
#include "CalcLibrary.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace mcalc
{
namespace ops_detail
{
inline Real OpsInf() { return Real::Infinity(); }
inline Real OpsEps() { return Real(1e-9); }   // logistics.js EPS

bool OpsWhole(Real X, int Lo, int Hi, int* Out)
{
	if (!(X >= Lo && X <= Hi && X == gm::floor(X)))
	{
		return false;
	}
	*Out = X.ToInt();
	return true;
}

bool OpsAllFinite(const Matrix& M)
{
	for (Real X : M.V)
	{
		if (!gm::isfinite(X))
		{
			return false;
		}
	}
	return true;
}

bool OpsSymmetric(const Matrix& A)
{
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int j = i + 1; j < A.Cols; ++j)
		{
			const Real X = A(i, j), Y = A(j, i);
			if (!(X == Y || gm::fabs(X - Y) <= 1e-12 * gm::max(Real(1), gm::max(gm::fabs(X), gm::fabs(Y)))))
			{
				return false;
			}
		}
	}
	return true;
}

// A square cost matrix: finite or +inf (no edge) off the diagonal; the diagonal is ignored.
bool OpsCostMatrix(const Matrix& W, std::string* Error)
{
	if (!W.IsSquare())
	{
		*Error = "needs a square cost matrix (inf = no edge)";
		return false;
	}
	for (Real X : W.V)
	{
		if (gm::isnan(X) || X == -OpsInf())
		{
			*Error = "costs must be numbers or inf (no edge)";
			return false;
		}
	}
	return true;
}

// A symmetric non-negative weighted adjacency matrix (conductances); the diagonal is ignored.
bool OpsAdjacency(const Matrix& A, std::string* Error)
{
	if (!A.IsSquare())
	{
		*Error = "needs a square adjacency matrix";
		return false;
	}
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int j = 0; j < A.Cols; ++j)
		{
			if (i != j && !(gm::isfinite(A(i, j)) && A(i, j) >= 0))
			{
				*Error = "edge weights must be finite and >= 0 (0 = no edge)";
				return false;
			}
		}
	}
	if (!OpsSymmetric(A))
	{
		*Error = "adjacency matrix must be symmetric (undirected graph)";
		return false;
	}
	return true;
}

Matrix OpsLaplacian(const Matrix& A)
{
	const int N = A.Rows;
	Matrix L(N, N);
	for (int i = 0; i < N; ++i)
	{
		Real Deg = 0.0;
		for (int j = 0; j < N; ++j)
		{
			if (j != i)
			{
				L(i, j) = -A(i, j);
				Deg += A(i, j);
			}
		}
		L(i, i) = Deg;
	}
	return L;
}

Matrix OpsDropRowCol(const Matrix& A, int K)
{
	const int N = A.Rows;
	Matrix B(N - 1, N - 1);
	for (int i = 0, r = 0; i < N; ++i)
	{
		if (i == K)
		{
			continue;
		}
		for (int j = 0, c = 0; j < N; ++j)
		{
			if (j == K)
			{
				continue;
			}
			B(r, c++) = A(i, j);
		}
		++r;
	}
	return B;
}

/* ---------- queues ---------- */
// Erlang B by the stable recursion B(k) = a B(k-1) / (k + a B(k-1)); Erlang C from it.
Real OpsErlangB(int C, Real A)
{
	Real B = 1.0;
	for (int k = 1; k <= C; ++k)
	{
		B = A * B / (k + A * B);
	}
	return B;
}
Real OpsErlangC(int C, Real A)
{
	if (A >= C)
	{
		return 1.0;
	}
	const Real B = OpsErlangB(C, A);
	return C * B / (C - A * (1.0 - B));
}

bool OpsServers(const Call& C, size_t I, int* Out, Result* Fail)
{
	if (!OpsWhole(C.N(I), 1, 10000, Out))
	{
		*Fail = Result::Fail("servers c must be a whole number 1..10000");
		return false;
	}
	return true;
}

Result ErlangBFn(const Call& C)
{
	int Servers = 0;
	Result F;
	if (!OpsServers(C, 0, &Servers, &F))
	{
		return F;
	}
	if (!(C.N(1) >= 0))
	{
		return Result::Fail("offered load a must be >= 0");
	}
	return Result::Of(OpsErlangB(Servers, C.N(1)));
}

Result ErlangCFn(const Call& C)
{
	int Servers = 0;
	Result F;
	if (!OpsServers(C, 0, &Servers, &F))
	{
		return F;
	}
	if (!(C.N(1) >= 0))
	{
		return Result::Fail("offered load a must be >= 0");
	}
	return Result::Of(OpsErlangC(Servers, C.N(1)));
}

Result MmcFn(const Call& C)
{
	const Real Lam = C.N(0), Mu = C.N(1);
	int Servers = 0;
	Result F;
	if (!OpsServers(C, 2, &Servers, &F))
	{
		return F;
	}
	if (!(Lam >= 0 && Mu > 0))
	{
		return Result::Fail("needs lambda >= 0 and mu > 0");
	}
	const Real A = Lam / Mu, Rho = A / Servers;
	if (Rho >= 1)
	{
		return Result::Fail("unstable: rho = lambda/(c mu) >= 1, no steady state");
	}
	const Real Cw = OpsErlangC(Servers, A), Wq = Cw / (Servers * Mu - Lam), W = Wq + 1.0 / Mu;
	return Result::Record({{"rho", Rho}, {"C", Cw}, {"Wq", Wq}, {"W", W}, {"L", Lam * W}, {"Lq", Lam * Wq}});
}

Result MmcTailFn(const Call& C)
{
	const Real Lam = C.N(0), Mu = C.N(1), T = C.N(3);
	int Servers = 0;
	Result F;
	if (!OpsServers(C, 2, &Servers, &F))
	{
		return F;
	}
	if (!(Lam >= 0 && Mu > 0 && T >= 0))
	{
		return Result::Fail("needs lambda >= 0, mu > 0, t >= 0");
	}
	if (Lam / (Servers * Mu) >= 1)
	{
		return Result::Fail("unstable: rho = lambda/(c mu) >= 1");
	}
	return Result::Of(OpsErlangC(Servers, Lam / Mu) * gm::exp(-(Servers * Mu - Lam) * T));
}

Result LittleFn(const Call& C) { return Result::Of(C.N(0) * C.N(1)); }

// M/G/1 by Pollaczek-Khinchine, with E[S^2] = (1 + cs2)/mu^2.
Result PkFn(const Call& C)
{
	const Real Lam = C.N(0), Mu = C.N(1), Cs2 = C.N(2);
	if (!(Lam >= 0 && Mu > 0 && Cs2 >= 0))
	{
		return Result::Fail("needs lambda >= 0, mu > 0, cs2 >= 0");
	}
	const Real Rho = Lam / Mu;
	if (Rho >= 1)
	{
		return Result::Fail("unstable: rho = lambda/mu >= 1");
	}
	const Real Wq = Lam * ((1.0 + Cs2) / (Mu * Mu)) / (2.0 * (1.0 - Rho)), W = Wq + 1.0 / Mu;
	return Result::Record({{"rho", Rho}, {"Wq", Wq}, {"W", W}, {"L", Lam * W}, {"Lq", Lam * Wq}});
}

Result KingmanFn(const Call& C)
{
	const Real Lam = C.N(0), Mu = C.N(1), Ca2 = C.N(2), Cs2 = C.N(3);
	if (!(Lam > 0 && Mu > 0 && Ca2 >= 0 && Cs2 >= 0))
	{
		return Result::Fail("needs lambda > 0, mu > 0, ca2 >= 0, cs2 >= 0");
	}
	const Real Rho = Lam / Mu;
	if (Rho >= 1)
	{
		return Result::Fail("unstable: rho = lambda/mu >= 1");
	}
	const Real Bound = Lam * (Ca2 / (Lam * Lam) + Cs2 / (Mu * Mu)) / (2.0 * (1.0 - Rho));
	const Real Approx = Rho / (1.0 - Rho) * ((Ca2 + Cs2) / 2.0) / Mu;
	return Result::Record({{"bound", Bound}, {"approx", Approx}});
}

/* ---------- inventory ---------- */
Result EoqFn(const Call& C)
{
	const Real K = C.N(0), D = C.N(1), H = C.N(2);
	if (!(K >= 0 && D >= 0 && H > 0))
	{
		return Result::Fail("needs K >= 0, D >= 0, h > 0");
	}
	return Result::Of(gm::sqrt(2.0 * K * D / H));
}

Result EoqCostFn(const Call& C)
{
	const Real Q = C.N(0), K = C.N(1), D = C.N(2), H = C.N(3);
	if (!(Q > 0))
	{
		return Result::Fail("order quantity Q must be > 0");
	}
	return Result::Of(K * D / Q + H * Q / 2.0);
}

Result FractileFn(const Call& C)
{
	const Real Price = C.N(0), Cost = C.N(1), Salv = C.N(2);
	if (!(Price > Salv))
	{
		return Result::Fail("price must exceed salvage");
	}
	return Result::Of((Price - Cost) / (Price - Salv));
}

// Poisson(mean) demand truncated at 150 terms, as in logistics.js newsvendor().
struct OpsNewsvendor
{
	std::vector<Real> Pm, Cdf;
	Real Price = 0, Cost = 0, Salv = 0;

	Real Profit(Real Q) const
	{
		Real E = 0.0;
		for (int d = 0; d < 150; ++d)
		{
			E += Pm[static_cast<size_t>(d)] * (Price * gm::min(Real(d), Q) + Salv * gm::max(Q - d, Real(0)));
		}
		return E - Cost * Q;
	}
};

bool OpsMakeNewsvendor(Real Price, Real Cost, Real Salv, Real Mean, OpsNewsvendor* Nv, std::string* Error)
{
	if (!(Price > Salv))
	{
		*Error = "price must exceed salvage";
		return false;
	}
	if (!(Mean > 0 && Mean <= 80))
	{
		*Error = "Poisson mean demand must be in (0, 80] (the pmf is truncated at 150)";
		return false;
	}
	Nv->Price = Price;
	Nv->Cost = Cost;
	Nv->Salv = Salv;
	Nv->Pm.assign(150, 0.0);
	Nv->Cdf.assign(150, 0.0);
	Nv->Pm[0] = gm::exp(-Mean);
	for (int k = 1; k < 150; ++k)
	{
		Nv->Pm[static_cast<size_t>(k)] = Nv->Pm[static_cast<size_t>(k - 1)] * Mean / k;
	}
	Real S = 0.0;
	for (int k = 0; k < 150; ++k)
	{
		S += Nv->Pm[static_cast<size_t>(k)];
		Nv->Cdf[static_cast<size_t>(k)] = S;
	}
	return true;
}

Result NewsvendorFn(const Call& C)
{
	OpsNewsvendor Nv;
	std::string E;
	if (!OpsMakeNewsvendor(C.N(0), C.N(1), C.N(2), C.N(3), &Nv, &E))
	{
		return Result::Fail(E);
	}
	const Real Crit = (Nv.Price - Nv.Cost) / (Nv.Price - Nv.Salv);
	int QStar = -1;
	for (int i = 0; i < 150; ++i)
	{
		if (Nv.Cdf[static_cast<size_t>(i)] >= Crit)
		{
			QStar = i;
			break;
		}
	}
	int QBest = 0;
	for (int q = 1; q < 120; ++q)
	{
		if (Nv.Profit(q) > Nv.Profit(QBest))
		{
			QBest = q;
		}
	}
	return Result::Record({{"crit", Crit}, {"q", Real(QStar)}, {"profit", Nv.Profit(QStar)}, {"qbest", Real(QBest)}});
}

Result NvProfitFn(const Call& C)
{
	OpsNewsvendor Nv;
	std::string E;
	if (!(C.N(0) >= 0))
	{
		return Result::Fail("order q must be >= 0");
	}
	if (!OpsMakeNewsvendor(C.N(1), C.N(2), C.N(3), C.N(4), &Nv, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Nv.Profit(C.N(0)));
}

// Periodic review with weekly demand Binomial(10, 0.4), stock grid -30..60 (logistics.js INV, invTables).
constexpr int OpsInvLo = -30;
constexpr int OpsInvHi = 60;
constexpr int OpsInvN = OpsInvHi - OpsInvLo + 1;

struct OpsInvTables
{
	Real Pmf[11];
	std::vector<Real> Lf;
	std::vector<int> Idx;   // OpsInvN x 11

	OpsInvTables(Real H, Real P) : Lf(OpsInvN, 0.0), Idx(static_cast<size_t>(OpsInvN) * 11, 0)
	{
		for (int d = 0; d <= 10; ++d)
		{
			Real Cb = 1.0;
			for (int i = 0; i < d; ++i)
			{
				Cb = Cb * (10 - i) / (i + 1);
			}
			Pmf[d] = gm::round(Cb) * gm::pow(0.4, d) * gm::pow(0.6, 10 - d);
		}
		for (int i = 0; i < OpsInvN; ++i)
		{
			const int Y = OpsInvLo + i;
			Real S = 0.0;
			for (int d = 0; d <= 10; ++d)
			{
				S += Pmf[d] * (H * std::max(Y - d, 0) + P * std::max(d - Y, 0));
				Idx[static_cast<size_t>(i) * 11 + static_cast<size_t>(d)] = std::max(Y - d, OpsInvLo) - OpsInvLo;
			}
			Lf[static_cast<size_t>(i)] = S;
		}
	}
	Real Expect(const std::vector<Real>& V, int I) const
	{
		Real E = 0.0;
		for (int d = 0; d <= 10; ++d)
		{
			E += Pmf[d] * V[static_cast<size_t>(Idx[static_cast<size_t>(I) * 11 + static_cast<size_t>(d)])];
		}
		return E;
	}
};

bool OpsInvParams(const Call& C, size_t First, Real* K, Real* P, Real* H, Real* Cu, Real* Alpha, std::string* Error)
{
	*K = C.N(First);
	*P = C.N(First + 1);
	*H = C.N(First + 2);
	*Cu = C.N(First + 3);
	*Alpha = C.N(First + 4);
	if (!(*K >= 0 && *P >= 0 && *H >= 0 && *Cu >= 0))
	{
		*Error = "costs K, p, h, c must be >= 0";
		return false;
	}
	if (!(*Alpha > 0 && *Alpha < 1))
	{
		*Error = "discount alpha must be in (0, 1)";
		return false;
	}
	return true;
}

// Value iteration for the discounted problem with setup cost K (logistics.js ssDP).
Result SsFn(const Call& C)
{
	Real K = 0, P = 0, H = 0, Cu = 0, Alpha = 0;
	std::string E;
	if (!OpsInvParams(C, 0, &K, &P, &H, &Cu, &Alpha, &E))
	{
		return Result::Fail(E);
	}
	const OpsInvTables T(H, P);
	const int N = OpsInvN;
	std::vector<Real> V(static_cast<size_t>(N), 0.0), G(static_cast<size_t>(N)), Suf(static_cast<size_t>(N)), Vn(static_cast<size_t>(N));
	auto ComputeG = [&]() {
		for (int i = 0; i < N; ++i)
		{
			G[static_cast<size_t>(i)] = Cu * (OpsInvLo + i) + T.Lf[static_cast<size_t>(i)] + Alpha * T.Expect(V, i);
		}
	};
	for (int It = 0; It < 3000; ++It)
	{
		ComputeG();
		Suf[static_cast<size_t>(N - 1)] = G[static_cast<size_t>(N - 1)];
		for (int i = N - 2; i >= 0; --i)
		{
			Suf[static_cast<size_t>(i)] = gm::min(G[static_cast<size_t>(i)], Suf[static_cast<size_t>(i + 1)]);
		}
		Real Df = 0.0;
		for (int i = 0; i < N; ++i)
		{
			Vn[static_cast<size_t>(i)] = gm::min(G[static_cast<size_t>(i)], K + Suf[static_cast<size_t>(i)]) - Cu * (OpsInvLo + i);
			Df = gm::max(Df, gm::fabs(Vn[static_cast<size_t>(i)] - V[static_cast<size_t>(i)]));
		}
		V.swap(Vn);
		if (Df < 1e-13)
		{
			break;
		}
	}
	ComputeG();
	int IS = 0;
	for (int i = 1; i < N; ++i)
	{
		if (G[static_cast<size_t>(i)] < G[static_cast<size_t>(IS)])
		{
			IS = i;
		}
	}
	const Real Gs = G[static_cast<size_t>(IS)];
	int Small = OpsInvLo + IS;
	for (int i = 0; i < N; ++i)
	{
		if (G[static_cast<size_t>(i)] <= K + Gs + 1e-12)
		{
			Small = OpsInvLo + i;
			break;
		}
	}
	return Result::Record({{"s", Real(Small)}, {"S", Real(OpsInvLo + IS)}, {"v0", V[static_cast<size_t>(-OpsInvLo)]}});
}

// Exact discounted cost from stock x of a fixed (s, S) policy (logistics.js ssEval).
Result SsCostFn(const Call& C)
{
	int Lo = 0, Up = 0, X0 = 0;
	if (!OpsWhole(C.N(0), OpsInvLo, OpsInvHi, &Lo) || !OpsWhole(C.N(1), OpsInvLo, OpsInvHi, &Up) || Lo > Up)
	{
		return Result::Fail("needs whole numbers -30 <= s <= S <= 60");
	}
	Real K = 0, P = 0, H = 0, Cu = 0, Alpha = 0;
	std::string E;
	if (!OpsInvParams(C, 2, &K, &P, &H, &Cu, &Alpha, &E))
	{
		return Result::Fail(E);
	}
	if (!OpsWhole(C.N(7), OpsInvLo, OpsInvHi, &X0))
	{
		return Result::Fail("starting stock x must be a whole number -30..60");
	}
	const OpsInvTables T(H, P);
	const int N = OpsInvN;
	std::vector<Real> R(static_cast<size_t>(N)), V(static_cast<size_t>(N), 0.0), W(static_cast<size_t>(N));
	std::vector<int> Next(static_cast<size_t>(N));
	for (int i = 0; i < N; ++i)
	{
		const int X = OpsInvLo + i, Y = X < Lo ? Up : X;
		R[static_cast<size_t>(i)] = (Y > X ? K + Cu * (Y - X) : 0.0) + T.Lf[static_cast<size_t>(Y - OpsInvLo)];
		Next[static_cast<size_t>(i)] = Y - OpsInvLo;
	}
	for (int It = 0; It < 5000; ++It)
	{
		Real Df = 0.0;
		for (int i = 0; i < N; ++i)
		{
			W[static_cast<size_t>(i)] = R[static_cast<size_t>(i)] + Alpha * T.Expect(V, Next[static_cast<size_t>(i)]);
			Df = gm::max(Df, gm::fabs(W[static_cast<size_t>(i)] - V[static_cast<size_t>(i)]));
		}
		V.swap(W);
		if (Df < 1e-12)
		{
			break;
		}
	}
	return Result::Of(V[static_cast<size_t>(X0 - OpsInvLo)]);
}

/* ---------- networks: max flow ---------- */
struct OpsFlow
{
	Real Value = 0.0;
	int Aug = 0;
	int Reach = 0;
	Matrix F;
};

// Edmonds-Karp on a capacity matrix (logistics.js maxflow).
void OpsMaxFlow(const Matrix& Cap, int S, int T, OpsFlow* Out)
{
	const int N = Cap.Rows;
	const Real Tiny = 1e-12;
	Matrix F(N, N);
	Real Total = 0.0;
	int Aug = 0;
	for (;;)
	{
		std::vector<int> Prev(static_cast<size_t>(N), -1), Q;
		Prev[static_cast<size_t>(S)] = S;
		Q.push_back(S);
		for (size_t h = 0; h < Q.size(); ++h)
		{
			const int X = Q[h];
			for (int y = 0; y < N; ++y)
			{
				if (Prev[static_cast<size_t>(y)] < 0 && Cap(X, y) - F(X, y) + F(y, X) > Tiny)
				{
					Prev[static_cast<size_t>(y)] = X;
					Q.push_back(y);
				}
			}
		}
		if (Prev[static_cast<size_t>(T)] < 0)
		{
			Out->Value = Total;
			Out->Aug = Aug;
			Out->Reach = static_cast<int>(Q.size());
			Out->F = F;
			return;
		}
		std::vector<int> Path = {T};
		while (Path.back() != S)
		{
			Path.push_back(Prev[static_cast<size_t>(Path.back())]);
		}
		std::reverse(Path.begin(), Path.end());
		Real Th = OpsInf();
		for (size_t k = 0; k + 1 < Path.size(); ++k)
		{
			const int A = Path[k], B = Path[k + 1];
			Th = gm::min(Th, Cap(A, B) - F(A, B) + F(B, A));
		}
		for (size_t k = 0; k + 1 < Path.size(); ++k)
		{
			const int A = Path[k], B = Path[k + 1];
			const Real Back = gm::min(F(B, A), Th);
			F(B, A) -= Back;
			F(A, B) += Th - Back;
		}
		Total += Th;
		++Aug;
	}
}

bool OpsFlowArgs(const Call& C, int* S, int* T, std::string* Error)
{
	const Matrix& Cap = C.M(0);
	if (!Cap.IsSquare() || Cap.Rows < 2)
	{
		*Error = "needs a square capacity matrix (0 = no lane)";
		return false;
	}
	for (Real X : Cap.V)
	{
		if (!(gm::isfinite(X) && X >= 0))
		{
			*Error = "capacities must be finite and >= 0";
			return false;
		}
	}
	const int N = Cap.Rows;
	const Real TArg = C.N(2) == 0 ? N : C.N(2);
	if (!OpsWhole(C.N(1), 1, N, S) || !OpsWhole(TArg, 1, N, T) || *S == *T)
	{
		*Error = "source and sink must be different nodes 1..n";
		return false;
	}
	*S -= 1;
	*T -= 1;
	return true;
}

Result MaxFlowFn(const Call& C)
{
	int S = 0, T = 0;
	std::string E;
	if (!OpsFlowArgs(C, &S, &T, &E))
	{
		return Result::Fail(E);
	}
	OpsFlow Fl;
	OpsMaxFlow(C.M(0), S, T, &Fl);
	return Result::Record({{"flow", Fl.Value}, {"aug", Real(Fl.Aug)}, {"ns", Real(Fl.Reach)}});
}

Result FlowMatFn(const Call& C)
{
	int S = 0, T = 0;
	std::string E;
	if (!OpsFlowArgs(C, &S, &T, &E))
	{
		return Result::Fail(E);
	}
	OpsFlow Fl;
	OpsMaxFlow(C.M(0), S, T, &Fl);
	return Result::Of(Fl.F);
}

/* ---------- networks: transportation simplex ---------- */
struct OpsTransport
{
	Matrix X;
	std::vector<Real> U, V;
	Real Cost = 0.0;
};

// North-west corner, potentials on the basis tree, Bland's entering rule, stepping-stone cycle by
// depth-first search, smallest-index leaving cell (logistics.js transport, step for step).
bool OpsTransportSolve(const std::vector<Real>& Sup, const std::vector<Real>& Dem, const Matrix& Cm, OpsTransport* Out, std::string* Error)
{
	const Real Eps = OpsEps();
	const int M = static_cast<int>(Sup.size()), N = static_cast<int>(Dem.size());
	std::vector<Real> S = Sup, D = Dem;
	Matrix X(M, N);
	std::vector<std::pair<int, int>> Basis;
	int I = 0, J = 0;
	for (;;)
	{
		const Real Q = gm::min(S[static_cast<size_t>(I)], D[static_cast<size_t>(J)]);
		X(I, J) = Q;
		Basis.push_back({I, J});
		S[static_cast<size_t>(I)] -= Q;
		D[static_cast<size_t>(J)] -= Q;
		if (I == M - 1 && J == N - 1)
		{
			break;
		}
		if (I == M - 1)
		{
			++J;
		}
		else if (J == N - 1)
		{
			++I;
		}
		else if (S[static_cast<size_t>(I)] <= Eps)
		{
			++I;
		}
		else
		{
			++J;
		}
	}
	auto InBasis = [&Basis](int A, int B) {
		for (const auto& P : Basis)
		{
			if (P.first == A && P.second == B)
			{
				return true;
			}
		}
		return false;
	};
	for (int It = 0; It < 1000; ++It)
	{
		std::vector<Real> U(static_cast<size_t>(M), 0.0), V(static_cast<size_t>(N), 0.0);
		std::vector<char> HasU(static_cast<size_t>(M), 0), HasV(static_cast<size_t>(N), 0);
		HasU[0] = 1;
		bool Changed = true;
		while (Changed)
		{
			Changed = false;
			for (const auto& P : Basis)
			{
				const size_t A = static_cast<size_t>(P.first), B = static_cast<size_t>(P.second);
				if (HasU[A] && !HasV[B])
				{
					V[B] = Cm(P.first, P.second) - U[A];
					HasV[B] = 1;
					Changed = true;
				}
				else if (HasV[B] && !HasU[A])
				{
					U[A] = Cm(P.first, P.second) - V[B];
					HasU[A] = 1;
					Changed = true;
				}
			}
		}
		int EntA = -1, EntB = -1;
		for (int a = 0; a < M && EntA < 0; ++a)
		{
			for (int b = 0; b < N; ++b)
			{
				if (!InBasis(a, b) && Cm(a, b) - U[static_cast<size_t>(a)] - V[static_cast<size_t>(b)] < -Eps)
				{
					EntA = a;
					EntB = b;
					break;
				}
			}
		}
		if (EntA < 0)
		{
			Real Cost = 0.0;
			for (int a = 0; a < M; ++a)
			{
				Real Row = 0.0;
				for (int b = 0; b < N; ++b)
				{
					Row += X(a, b) * Cm(a, b);
				}
				Cost += Row;
			}
			Out->X = X;
			Out->U = U;
			Out->V = V;
			Out->Cost = Cost;
			return true;
		}
		// Rows are nodes 0..M-1, columns M..M+N-1; adjacency in basis order.
		std::vector<std::vector<int>> Adj(static_cast<size_t>(M + N));
		for (const auto& P : Basis)
		{
			Adj[static_cast<size_t>(P.first)].push_back(M + P.second);
			Adj[static_cast<size_t>(M + P.second)].push_back(P.first);
		}
		const int Start = EntA, Goal = M + EntB;
		std::vector<int> Prev(static_cast<size_t>(M + N), -1);   // -1 unvisited, -2 the start
		Prev[static_cast<size_t>(Start)] = -2;
		std::vector<int> Stack = {Start};
		while (!Stack.empty())
		{
			const int Xn = Stack.back();
			Stack.pop_back();
			if (Xn == Goal)
			{
				break;
			}
			for (int Y : Adj[static_cast<size_t>(Xn)])
			{
				if (Prev[static_cast<size_t>(Y)] == -1)
				{
					Prev[static_cast<size_t>(Y)] = Xn;
					Stack.push_back(Y);
				}
			}
		}
		if (Prev[static_cast<size_t>(Goal)] == -1)
		{
			*Error = "transportation simplex: basis is not a spanning tree";
			return false;
		}
		std::vector<int> Path = {Goal};
		while (Prev[static_cast<size_t>(Path.back())] != -2)
		{
			Path.push_back(Prev[static_cast<size_t>(Path.back())]);
		}
		std::reverse(Path.begin(), Path.end());
		std::vector<std::pair<int, int>> Minus, Plus;
		for (size_t k = 0; k + 1 < Path.size(); ++k)
		{
			const int P = Path[k], Q = Path[k + 1];
			const std::pair<int, int> Cell = P < M ? std::make_pair(P, Q - M) : std::make_pair(Q, P - M);
			(k % 2 == 0 ? Minus : Plus).push_back(Cell);
		}
		Real Th = OpsInf();
		for (const auto& Cell : Minus)
		{
			Th = gm::min(Th, X(Cell.first, Cell.second));
		}
		std::pair<int, int> Leave = {-1, -1};
		for (const auto& Cell : Minus)
		{
			if (X(Cell.first, Cell.second) <= Th + Eps && (Leave.first < 0 || Cell < Leave))
			{
				Leave = Cell;
			}
		}
		for (const auto& Cell : Minus)
		{
			X(Cell.first, Cell.second) -= Th;
		}
		for (const auto& Cell : Plus)
		{
			X(Cell.first, Cell.second) += Th;
		}
		X(EntA, EntB) += Th;
		for (size_t k = 0; k < Basis.size(); ++k)
		{
			if (Basis[k] == Leave)
			{
				Basis.erase(Basis.begin() + static_cast<std::ptrdiff_t>(k));
				break;
			}
		}
		Basis.push_back({EntA, EntB});
		X(Leave.first, Leave.second) = 0.0;
	}
	*Error = "transportation simplex did not terminate";
	return false;
}

bool OpsTransportArgs(const Call& C, OpsTransport* Out, std::string* Error)
{
	const Matrix &Sm = C.M(0), &Dm = C.M(1), &Cm = C.M(2);
	if (!Sm.IsVector() || !Dm.IsVector())
	{
		*Error = "supply and demand must be vectors";
		return false;
	}
	if (Cm.Rows != Sm.Count() || Cm.Cols != Dm.Count())
	{
		*Error = "cost matrix must be (number of supplies) x (number of demands)";
		return false;
	}
	if (!OpsAllFinite(Sm) || !OpsAllFinite(Dm) || !OpsAllFinite(Cm))
	{
		*Error = "supplies, demands and costs must be finite";
		return false;
	}
	Real Ts = 0.0, Td = 0.0;
	for (Real X : Sm.V)
	{
		if (X < 0)
		{
			*Error = "supplies must be >= 0";
			return false;
		}
		Ts += X;
	}
	for (Real X : Dm.V)
	{
		if (X < 0)
		{
			*Error = "demands must be >= 0";
			return false;
		}
		Td += X;
	}
	if (gm::fabs(Ts - Td) > 1e-9 * gm::max(Real(1), Ts))
	{
		*Error = "total supply must equal total demand";
		return false;
	}
	return OpsTransportSolve(Sm.V, Dm.V, Cm, Out, Error);
}

Result TransportFn(const Call& C)
{
	OpsTransport T;
	std::string E;
	if (!OpsTransportArgs(C, &T, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(T.Cost);
}

Result TPlanFn(const Call& C)
{
	OpsTransport T;
	std::string E;
	if (!OpsTransportArgs(C, &T, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(T.X);
}

Result TPricesFn(const Call& C)
{
	OpsTransport T;
	std::string E;
	if (!OpsTransportArgs(C, &T, &E))
	{
		return Result::Fail(E);
	}
	std::vector<Real> UV = T.U;
	UV.insert(UV.end(), T.V.begin(), T.V.end());
	return Result::Of(Matrix::Row(UV));
}

/* ---------- networks: assignment (Hungarian method with potentials) ---------- */
bool OpsAssign(const Matrix& Cm, std::vector<int>* Perm, Real* Cost, std::string* Error)
{
	if (!Cm.IsSquare())
	{
		*Error = "needs a square cost matrix";
		return false;
	}
	if (!OpsAllFinite(Cm))
	{
		*Error = "costs must be finite";
		return false;
	}
	const int N = Cm.Rows;
	std::vector<Real> U(static_cast<size_t>(N + 1), 0.0), V(static_cast<size_t>(N + 1), 0.0);
	std::vector<int> P(static_cast<size_t>(N + 1), 0), Way(static_cast<size_t>(N + 1), 0);
	for (int i = 1; i <= N; ++i)
	{
		P[0] = i;
		int J0 = 0;
		std::vector<Real> MinV(static_cast<size_t>(N + 1), OpsInf());
		std::vector<char> Used(static_cast<size_t>(N + 1), 0);
		do
		{
			Used[static_cast<size_t>(J0)] = 1;
			const int I0 = P[static_cast<size_t>(J0)];
			Real Delta = OpsInf();
			int J1 = 0;
			for (int j = 1; j <= N; ++j)
			{
				if (!Used[static_cast<size_t>(j)])
				{
					const Real Cur = Cm(I0 - 1, j - 1) - U[static_cast<size_t>(I0)] - V[static_cast<size_t>(j)];
					if (Cur < MinV[static_cast<size_t>(j)])
					{
						MinV[static_cast<size_t>(j)] = Cur;
						Way[static_cast<size_t>(j)] = J0;
					}
					if (MinV[static_cast<size_t>(j)] < Delta)
					{
						Delta = MinV[static_cast<size_t>(j)];
						J1 = j;
					}
				}
			}
			for (int j = 0; j <= N; ++j)
			{
				if (Used[static_cast<size_t>(j)])
				{
					U[static_cast<size_t>(P[static_cast<size_t>(j)])] += Delta;
					V[static_cast<size_t>(j)] -= Delta;
				}
				else
				{
					MinV[static_cast<size_t>(j)] -= Delta;
				}
			}
			J0 = J1;
		} while (P[static_cast<size_t>(J0)] != 0);
		do
		{
			const int J1 = Way[static_cast<size_t>(J0)];
			P[static_cast<size_t>(J0)] = P[static_cast<size_t>(J1)];
			J0 = J1;
		} while (J0 != 0);
	}
	Perm->assign(static_cast<size_t>(N), 0);
	for (int j = 1; j <= N; ++j)
	{
		(*Perm)[static_cast<size_t>(P[static_cast<size_t>(j)] - 1)] = j - 1;
	}
	Real S = 0.0;
	for (int i = 0; i < N; ++i)
	{
		S += Cm(i, (*Perm)[static_cast<size_t>(i)]);
	}
	*Cost = S;
	return true;
}

Result AssignFn(const Call& C)
{
	std::vector<int> Perm;
	Real Cost = 0.0;
	std::string E;
	if (!OpsAssign(C.M(0), &Perm, &Cost, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Cost);
}

Result AssignPermFn(const Call& C)
{
	std::vector<int> Perm;
	Real Cost = 0.0;
	std::string E;
	if (!OpsAssign(C.M(0), &Perm, &Cost, &E))
	{
		return Result::Fail(E);
	}
	std::vector<Real> Out;
	for (int P : Perm)
	{
		Out.push_back(P + 1.0);
	}
	return Result::Of(Matrix::Row(Out));
}

/* ---------- networks: max-plus ---------- */
bool OpsMaxPlusMatrix(const Matrix& A, std::string* Error)
{
	for (Real X : A.V)
	{
		if (gm::isnan(X) || X == OpsInf())
		{
			*Error = "max-plus entries must be numbers or -inf (no arc)";
			return false;
		}
	}
	return true;
}

Matrix OpsOtimes(const Matrix& A, const Matrix& B)
{
	Matrix Out(A.Rows, B.Cols, -OpsInf());
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int j = 0; j < B.Cols; ++j)
		{
			Real Best = -OpsInf();
			for (int k = 0; k < A.Cols; ++k)
			{
				Best = gm::max(Best, A(i, k) + B(k, j));
			}
			Out(i, j) = Best;
		}
	}
	return Out;
}

// Karp's maximum cycle mean from node 1 (logistics.js karp).
Result KarpFn(const Call& C)
{
	const Matrix& A = C.M(0);
	std::string E;
	if (!A.IsSquare())
	{
		return Result::Fail("needs a square max-plus matrix");
	}
	if (!OpsMaxPlusMatrix(A, &E))
	{
		return Result::Fail(E);
	}
	const int N = A.Rows;
	std::vector<Matrix> D;
	Matrix D0(N, 1, -OpsInf());
	D0.V[0] = 0.0;
	D.push_back(D0);
	for (int k = 1; k <= N; ++k)
	{
		D.push_back(OpsOtimes(A, D.back()));
	}
	Real Best = -OpsInf();
	for (int v = 0; v < N; ++v)
	{
		const Real Dn = D[static_cast<size_t>(N)].V[static_cast<size_t>(v)];
		if (Dn == -OpsInf())
		{
			continue;
		}
		Real M = OpsInf();
		for (int k = 0; k < N; ++k)
		{
			const Real Dk = D[static_cast<size_t>(k)].V[static_cast<size_t>(v)];
			if (Dk > -OpsInf())
			{
				M = gm::min(M, (Dn - Dk) / (N - k));
			}
		}
		Best = gm::max(Best, M);
	}
	return Result::Of(Best);
}

Result MpMulFn(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1);
	std::string E;
	if (A.Cols != B.Rows || A.IsEmpty() || B.IsEmpty())
	{
		return Result::Fail("inner dimensions must agree (A cols = B rows)");
	}
	if (!OpsMaxPlusMatrix(A, &E) || !OpsMaxPlusMatrix(B, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(OpsOtimes(A, B));
}

Result MpPowFn(const Call& C)
{
	const Matrix& A = C.M(0);
	std::string E;
	int K = 0;
	if (!A.IsSquare())
	{
		return Result::Fail("needs a square max-plus matrix");
	}
	if (!OpsMaxPlusMatrix(A, &E))
	{
		return Result::Fail(E);
	}
	if (!OpsWhole(C.N(1), 0, 1000, &K))
	{
		return Result::Fail("power k must be a whole number 0..1000");
	}
	Matrix P(A.Rows, A.Rows, -OpsInf());
	for (int i = 0; i < A.Rows; ++i)
	{
		P(i, i) = 0.0;
	}
	for (int k = 0; k < K; ++k)
	{
		P = OpsOtimes(A, P);
	}
	return Result::Of(P);
}

/* ---------- networks: optimal transport on a line ---------- */
// W1 = integral |F - G| of the cumulative masses (logistics.js w1cdf).
Real OpsW1Cdf(const std::vector<Real>& X, const std::vector<Real>& A, const std::vector<Real>& Y, const std::vector<Real>& B)
{
	std::vector<Real> Pts = X;
	Pts.insert(Pts.end(), Y.begin(), Y.end());
	std::sort(Pts.begin(), Pts.end());
	Pts.erase(std::unique(Pts.begin(), Pts.end()), Pts.end());
	Real F = 0.0, G = 0.0, W = 0.0;
	for (size_t k = 0; k + 1 < Pts.size(); ++k)
	{
		for (size_t i = 0; i < X.size(); ++i)
		{
			if (X[i] == Pts[k])
			{
				F += A[i];
			}
		}
		for (size_t j = 0; j < Y.size(); ++j)
		{
			if (Y[j] == Pts[k])
			{
				G += B[j];
			}
		}
		W += gm::fabs(F - G) * (Pts[k + 1] - Pts[k]);
	}
	return W;
}

Result W1Fn(const Call& C)
{
	const Matrix &X = C.M(0), &Y = C.M(1);
	if (!X.IsVector() || !Y.IsVector() || !OpsAllFinite(X) || !OpsAllFinite(Y))
	{
		return Result::Fail("needs two vectors of finite samples");
	}
	if (X.Count() * Y.Count() > 4000000)
	{
		return Result::Fail("too many samples");
	}
	const std::vector<Real> A(X.V.size(), 1.0 / X.Count()), B(Y.V.size(), 1.0 / Y.Count());
	return Result::Of(OpsW1Cdf(X.V, A, Y.V, B));
}

bool OpsMassArgs(const Call& C, std::string* Error)
{
	const Matrix &X = C.M(0), &A = C.M(1), &Y = C.M(2), &B = C.M(3);
	if (!X.IsVector() || !Y.IsVector() || A.Count() != X.Count() || B.Count() != Y.Count())
	{
		*Error = "needs points x, masses a (same length), points y, masses b (same length)";
		return false;
	}
	if (!OpsAllFinite(X) || !OpsAllFinite(Y) || !OpsAllFinite(A) || !OpsAllFinite(B))
	{
		*Error = "points and masses must be finite";
		return false;
	}
	for (Real M : A.V)
	{
		if (M < 0)
		{
			*Error = "masses must be >= 0";
			return false;
		}
	}
	for (Real M : B.V)
	{
		if (M < 0)
		{
			*Error = "masses must be >= 0";
			return false;
		}
	}
	return true;
}

Result W1wFn(const Call& C)
{
	std::string E;
	if (!OpsMassArgs(C, &E))
	{
		return Result::Fail(E);
	}
	if (C.M(0).Count() * C.M(2).Count() > 4000000)
	{
		return Result::Fail("too many points");
	}
	return Result::Of(OpsW1Cdf(C.M(0).V, C.M(1).V, C.M(2).V, C.M(3).V));
}

// Cost sum P_ij |x_i - y_j|^p of the monotone (north-west corner on sorted points) coupling.
Result OtMonoFn(const Call& C)
{
	std::string E;
	if (!OpsMassArgs(C, &E))
	{
		return Result::Fail(E);
	}
	const Real P = C.N(4);
	if (!(P > 0))
	{
		return Result::Fail("exponent p must be > 0");
	}
	const std::vector<Real> &X = C.M(0).V, &Y = C.M(2).V;
	std::vector<size_t> Ix(X.size()), Iy(Y.size());
	std::iota(Ix.begin(), Ix.end(), size_t(0));
	std::iota(Iy.begin(), Iy.end(), size_t(0));
	std::stable_sort(Ix.begin(), Ix.end(), [&X](size_t p, size_t q) { return X[p] < X[q]; });
	std::stable_sort(Iy.begin(), Iy.end(), [&Y](size_t p, size_t q) { return Y[p] < Y[q]; });
	std::vector<Real> A, B;
	Real Ta = 0.0, Tb = 0.0;
	for (size_t i : Ix)
	{
		A.push_back(C.M(1).V[i]);
		Ta += C.M(1).V[i];
	}
	for (size_t j : Iy)
	{
		B.push_back(C.M(3).V[j]);
		Tb += C.M(3).V[j];
	}
	if (gm::fabs(Ta - Tb) > 1e-9 * gm::max(Real(1), Ta))
	{
		return Result::Fail("total masses must be equal");
	}
	const size_t M = A.size(), N = B.size();
	Matrix Plan(static_cast<int>(M), static_cast<int>(N));
	size_t i = 0, j = 0;
	while (i < M && j < N)
	{
		const Real Q = gm::min(A[i], B[j]);
		Plan(static_cast<int>(i), static_cast<int>(j)) += Q;
		A[i] -= Q;
		B[j] -= Q;
		if (A[i] <= 1e-12)
		{
			++i;
		}
		else
		{
			++j;
		}
	}
	Real Cost = 0.0;
	for (size_t r = 0; r < M; ++r)
	{
		Real Row = 0.0;
		for (size_t c = 0; c < N; ++c)
		{
			Row += Plan(static_cast<int>(r), static_cast<int>(c)) * gm::pow(gm::fabs(X[Ix[r]] - Y[Iy[c]]), P);
		}
		Cost += Row;
	}
	return Result::Of(Cost);
}

/* ---------- networks: trees and tours ---------- */
Result MstFn(const Call& C)
{
	const Matrix& W = C.M(0);
	std::string E;
	if (!OpsCostMatrix(W, &E))
	{
		return Result::Fail(E);
	}
	if (!OpsSymmetric(W))
	{
		return Result::Fail("weight matrix must be symmetric (undirected graph)");
	}
	const int N = W.Rows;
	struct Edge
	{
		Real W;
		int U, V;
	};
	std::vector<Edge> Edges;
	for (int i = 0; i < N; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			if (gm::isfinite(W(i, j)))
			{
				Edges.push_back({W(i, j), i, j});
			}
		}
	}
	std::stable_sort(Edges.begin(), Edges.end(), [](const Edge& P, const Edge& Q) { return P.W < Q.W; });
	std::vector<int> Par(static_cast<size_t>(N));
	std::iota(Par.begin(), Par.end(), 0);
	auto Find = [&Par](int X) {
		while (Par[static_cast<size_t>(X)] != X)
		{
			X = Par[static_cast<size_t>(X)];
		}
		return X;
	};
	Real Total = 0.0;
	int Count = 0;
	for (const Edge& Ed : Edges)
	{
		const int A = Find(Ed.U), B = Find(Ed.V);
		if (A != B)
		{
			Par[static_cast<size_t>(A)] = B;
			Total += Ed.W;
			++Count;
		}
	}
	if (Count != N - 1)
	{
		return Result::Fail("graph is disconnected: no spanning tree");
	}
	return Result::Of(Total);
}

Result DistMatFn(const Call& C)
{
	const Matrix& P = C.M(0);
	if (P.Cols != 2 || P.Rows < 1)
	{
		return Result::Fail("needs an n x 2 matrix of points");
	}
	const int N = P.Rows;
	Matrix D(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			D(i, j) = gm::hypot(P(i, 0) - P(j, 0), P(i, 1) - P(j, 1));
		}
	}
	return Result::Of(D);
}

Real OpsTourLen(const Matrix& D, const std::vector<int>& T)
{
	Real L = 0.0;
	for (size_t k = 0; k < T.size(); ++k)
	{
		L += D(T[k], T[(k + 1) % T.size()]);
	}
	return L;
}

Result TourLenFn(const Call& C)
{
	const Matrix &D = C.M(0), &Tm = C.M(1);
	std::string E;
	if (!OpsCostMatrix(D, &E))
	{
		return Result::Fail(E);
	}
	if (!Tm.IsVector())
	{
		return Result::Fail("tour must be a vector of stops 1..n");
	}
	std::vector<int> T;
	for (Real X : Tm.V)
	{
		int I = 0;
		if (!OpsWhole(X, 1, D.Rows, &I))
		{
			return Result::Fail("tour must be a vector of stops 1..n");
		}
		T.push_back(I - 1);
	}
	return Result::Of(OpsTourLen(D, T));
}

// Held-Karp dynamic programming: the exact shortest closed tour through every stop.
Result TspFn(const Call& C)
{
	const Matrix& D = C.M(0);
	std::string E;
	if (!OpsCostMatrix(D, &E))
	{
		return Result::Fail(E);
	}
	const int N = D.Rows;
	if (N > 15)
	{
		return Result::Fail("Held-Karp is exponential: at most 15 stops");
	}
	if (N == 1)
	{
		return Result::Of(0.0);
	}
	const int M = N - 1, Full = (1 << M) - 1;
	std::vector<Real> Dp(static_cast<size_t>(1 << M) * static_cast<size_t>(M), OpsInf());
	auto At = [&Dp, M](int Mask, int J) -> Real& { return Dp[static_cast<size_t>(Mask) * static_cast<size_t>(M) + static_cast<size_t>(J)]; };
	for (int j = 0; j < M; ++j)
	{
		At(1 << j, j) = D(0, j + 1);
	}
	for (int Mask = 1; Mask <= Full; ++Mask)
	{
		for (int j = 0; j < M; ++j)
		{
			const Real Cur = At(Mask, j);
			if (!(Mask & (1 << j)) || Cur == OpsInf())
			{
				continue;
			}
			for (int k = 0; k < M; ++k)
			{
				if (Mask & (1 << k))
				{
					continue;
				}
				Real& Nx = At(Mask | (1 << k), k);
				Nx = gm::min(Nx, Cur + D(j + 1, k + 1));
			}
		}
	}
	Real Best = OpsInf();
	for (int j = 0; j < M; ++j)
	{
		Best = gm::min(Best, At(Full, j) + D(j + 1, 0));
	}
	return Result::Of(Best);
}

// Prim with logistics.js prim()'s tie-breaking (lighter, then smaller u, then smaller v).
std::vector<std::pair<int, int>> OpsPrim(const Matrix& D, const std::vector<int>& Nodes)
{
	std::vector<int> InT = {Nodes[0]};
	std::vector<char> Member(static_cast<size_t>(D.Rows), 0);
	Member[static_cast<size_t>(Nodes[0])] = 1;
	std::vector<std::pair<int, int>> E;
	while (InT.size() < Nodes.size())
	{
		bool Have = false;
		Real Bw = 0.0;
		int Bu = 0, Bv = 0;
		for (int U : InT)
		{
			for (int V : Nodes)
			{
				if (Member[static_cast<size_t>(V)])
				{
					continue;
				}
				const Real W = D(U, V);
				if (!Have || W < Bw || (W == Bw && (U < Bu || (U == Bu && V < Bv))))
				{
					Have = true;
					Bw = W;
					Bu = U;
					Bv = V;
				}
			}
		}
		E.push_back({Bu, Bv});
		InT.push_back(Bv);
		Member[static_cast<size_t>(Bv)] = 1;
	}
	return E;
}

struct OpsOneTree
{
	Real Val = 0.0;
	std::vector<int> Deg;
};

OpsOneTree OpsMakeOneTree(const Matrix& D, const std::vector<Real>& Pi)
{
	const int N = D.Rows;
	Matrix Dp(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Dp(i, j) = D(i, j) + Pi[static_cast<size_t>(i)] + Pi[static_cast<size_t>(j)];
		}
	}
	std::vector<int> Rest(static_cast<size_t>(N - 1));
	std::iota(Rest.begin(), Rest.end(), 1);
	const std::vector<std::pair<int, int>> T = OpsPrim(Dp, Rest);
	OpsOneTree O;
	O.Deg.assign(static_cast<size_t>(N), 0);
	Real Sum = 0.0;
	for (const auto& Ed : T)
	{
		++O.Deg[static_cast<size_t>(Ed.first)];
		++O.Deg[static_cast<size_t>(Ed.second)];
		Sum += Dp(Ed.first, Ed.second);
	}
	std::vector<int> Two = Rest;
	std::sort(Two.begin(), Two.end(), [&Dp](int P, int Q) { return Dp(0, P) < Dp(0, Q) || (Dp(0, P) == Dp(0, Q) && P < Q); });
	Real PiSum = 0.0;
	for (Real P : Pi)
	{
		PiSum += P;
	}
	for (int k = 0; k < 2; ++k)
	{
		++O.Deg[0];
		++O.Deg[static_cast<size_t>(Two[static_cast<size_t>(k)])];
	}
	O.Val = Sum + Dp(0, Two[0]) + Dp(0, Two[1]) - 2.0 * PiSum;
	return O;
}

bool OpsTourMatrix(const Matrix& D, int MaxN, std::string* Error)
{
	if (!D.IsSquare() || D.Rows < 3 || D.Rows > MaxN)
	{
		*Error = "needs a square distance matrix with 3.." + std::to_string(MaxN) + " stops";
		return false;
	}
	if (!OpsAllFinite(D))
	{
		*Error = "distances must be finite";
		return false;
	}
	return true;
}

// Held-Karp lower bound by subgradient ascent on 1-tree penalties (logistics.js heldKarp).
Result HkBoundFn(const Call& C)
{
	const Matrix& D = C.M(0);
	std::string E;
	int Iters = 0;
	if (!OpsTourMatrix(D, 60, &E))
	{
		return Result::Fail(E);
	}
	if (!OpsWhole(C.N(1), 1, 2000, &Iters))
	{
		return Result::Fail("iterations must be a whole number 1..2000");
	}
	const int N = D.Rows;
	std::vector<Real> Pi(static_cast<size_t>(N), 0.0);
	const Real L0 = OpsMakeOneTree(D, Pi).Val;
	Real Best = -OpsInf(), Step = 0.1 * L0 / N;
	for (int k = 0; k < Iters; ++k)
	{
		const OpsOneTree O = OpsMakeOneTree(D, Pi);
		Best = gm::max(Best, O.Val);
		bool AllTwo = true;
		for (int Dg : O.Deg)
		{
			AllTwo = AllTwo && Dg == 2;
		}
		if (AllTwo)
		{
			break;
		}
		for (int i = 0; i < N; ++i)
		{
			Pi[static_cast<size_t>(i)] = Pi[static_cast<size_t>(i)] + Step * (O.Deg[static_cast<size_t>(i)] - 2);
		}
		Step *= 0.97;
	}
	return Result::Record({{"lb", Best}, {"onetree", L0}});
}

// Exact minimum-weight perfect matching by recursion on the first stop (logistics.js mwpm).
Real OpsMwpm(const Matrix& D, const std::vector<int>& Odd, std::vector<std::pair<int, int>>* Match)
{
	Match->clear();
	if (Odd.empty())
	{
		return 0.0;
	}
	const int A = Odd[0];
	Real Best = OpsInf();
	for (size_t k = 1; k < Odd.size(); ++k)
	{
		std::vector<int> Rest(Odd.begin() + 1, Odd.begin() + static_cast<std::ptrdiff_t>(k));
		Rest.insert(Rest.end(), Odd.begin() + static_cast<std::ptrdiff_t>(k + 1), Odd.end());
		std::vector<std::pair<int, int>> Sub;
		const Real R = OpsMwpm(D, Rest, &Sub);
		if (D(A, Odd[k]) + R < Best - 1e-15)
		{
			Best = D(A, Odd[k]) + R;
			Match->assign(1, {A, Odd[k]});
			Match->insert(Match->end(), Sub.begin(), Sub.end());
		}
	}
	return Best;
}

// Hierholzer from stop 0, neighbours in ascending order (logistics.js euler).
std::vector<int> OpsEuler(int N, const std::vector<std::pair<int, int>>& E)
{
	std::vector<std::vector<std::pair<int, int>>> Adj(static_cast<size_t>(N));
	for (size_t k = 0; k < E.size(); ++k)
	{
		Adj[static_cast<size_t>(E[k].first)].push_back({E[k].second, static_cast<int>(k)});
		Adj[static_cast<size_t>(E[k].second)].push_back({E[k].first, static_cast<int>(k)});
	}
	for (auto& L : Adj)
	{
		std::sort(L.begin(), L.end());
	}
	std::vector<size_t> Head(static_cast<size_t>(N), 0);
	std::vector<char> Used(E.size(), 0);
	std::vector<int> St = {0}, Circ;
	while (!St.empty())
	{
		const size_t V = static_cast<size_t>(St.back());
		while (Head[V] < Adj[V].size() && Used[static_cast<size_t>(Adj[V][Head[V]].second)])
		{
			++Head[V];
		}
		if (Head[V] < Adj[V].size())
		{
			const std::pair<int, int> Wk = Adj[V][Head[V]++];
			Used[static_cast<size_t>(Wk.second)] = 1;
			St.push_back(Wk.first);
		}
		else
		{
			Circ.push_back(St.back());
			St.pop_back();
		}
	}
	std::reverse(Circ.begin(), Circ.end());
	return Circ;
}

Result ChristofidesFn(const Call& C)
{
	const Matrix& D = C.M(0);
	std::string E;
	if (!OpsTourMatrix(D, 14, &E))
	{
		return Result::Fail(E);
	}
	const int N = D.Rows;
	std::vector<int> All(static_cast<size_t>(N));
	std::iota(All.begin(), All.end(), 0);
	const std::vector<std::pair<int, int>> T = OpsPrim(D, All);
	std::vector<int> Deg(static_cast<size_t>(N), 0), Odd;
	Real Mst = 0.0;
	for (const auto& Ed : T)
	{
		++Deg[static_cast<size_t>(Ed.first)];
		++Deg[static_cast<size_t>(Ed.second)];
		Mst += D(Ed.first, Ed.second);
	}
	for (int i = 0; i < N; ++i)
	{
		if (Deg[static_cast<size_t>(i)] % 2)
		{
			Odd.push_back(i);
		}
	}
	std::vector<std::pair<int, int>> Match;
	const Real MatchCost = OpsMwpm(D, Odd, &Match);
	std::vector<std::pair<int, int>> Multi = T;
	Multi.insert(Multi.end(), Match.begin(), Match.end());
	std::vector<int> Seen;
	for (int V : OpsEuler(N, Multi))
	{
		if (std::find(Seen.begin(), Seen.end(), V) == Seen.end())
		{
			Seen.push_back(V);
		}
	}
	return Result::Record({{"len", OpsTourLen(D, Seen)}, {"mst", Mst}, {"match", MatchCost}});
}

/* ---------- graphs (Chapter 9) ---------- */
bool OpsSource(const Call& C, size_t I, int N, int* S)
{
	if (!OpsWhole(C.N(I), 1, N, S))
	{
		return false;
	}
	*S -= 1;
	return true;
}

Result DijkstraFn(const Call& C)
{
	const Matrix& W = C.M(0);
	std::string E;
	int S = 0;
	if (!OpsCostMatrix(W, &E))
	{
		return Result::Fail(E);
	}
	const int N = W.Rows;
	if (!OpsSource(C, 1, N, &S))
	{
		return Result::Fail("source must be a node 1..n");
	}
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if (i != j && W(i, j) < 0)
			{
				return Result::Fail("negative edge: Dijkstra needs costs >= 0 (use bellman)");
			}
		}
	}
	std::vector<Real> D(static_cast<size_t>(N), OpsInf());
	std::vector<char> Done(static_cast<size_t>(N), 0);
	D[static_cast<size_t>(S)] = 0.0;
	for (int Round = 0; Round < N; ++Round)
	{
		int U = -1;
		for (int i = 0; i < N; ++i)
		{
			if (!Done[static_cast<size_t>(i)] && D[static_cast<size_t>(i)] < OpsInf() && (U < 0 || D[static_cast<size_t>(i)] < D[static_cast<size_t>(U)]))
			{
				U = i;
			}
		}
		if (U < 0)
		{
			break;
		}
		Done[static_cast<size_t>(U)] = 1;
		for (int v = 0; v < N; ++v)
		{
			if (v != U && !Done[static_cast<size_t>(v)] && D[static_cast<size_t>(U)] + W(U, v) < D[static_cast<size_t>(v)])
			{
				D[static_cast<size_t>(v)] = D[static_cast<size_t>(U)] + W(U, v);
			}
		}
	}
	return Result::Of(Matrix::Row(D));
}

Result BellmanFn(const Call& C)
{
	const Matrix& W = C.M(0);
	std::string E;
	int S = 0;
	if (!OpsCostMatrix(W, &E))
	{
		return Result::Fail(E);
	}
	const int N = W.Rows;
	if (!OpsSource(C, 1, N, &S))
	{
		return Result::Fail("source must be a node 1..n");
	}
	std::vector<Real> D(static_cast<size_t>(N), OpsInf());
	D[static_cast<size_t>(S)] = 0.0;
	for (int Round = 0; Round < N; ++Round)
	{
		bool Changed = false;
		for (int u = 0; u < N; ++u)
		{
			for (int v = 0; v < N; ++v)
			{
				if (u != v && gm::isfinite(W(u, v)) && D[static_cast<size_t>(u)] + W(u, v) < D[static_cast<size_t>(v)])
				{
					D[static_cast<size_t>(v)] = D[static_cast<size_t>(u)] + W(u, v);
					Changed = true;
				}
			}
		}
		if (!Changed)
		{
			return Result::Of(Matrix::Row(D));
		}
	}
	return Result::Fail("negative cycle reachable from the source");
}

Result FloydFn(const Call& C)
{
	const Matrix& W = C.M(0);
	std::string E;
	if (!OpsCostMatrix(W, &E))
	{
		return Result::Fail(E);
	}
	const int N = W.Rows;
	Matrix D = W;
	for (int i = 0; i < N; ++i)
	{
		D(i, i) = 0.0;
	}
	for (int k = 0; k < N; ++k)
	{
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				if (D(i, k) + D(k, j) < D(i, j))
				{
					D(i, j) = D(i, k) + D(k, j);
				}
			}
		}
	}
	for (int i = 0; i < N; ++i)
	{
		if (D(i, i) < 0)
		{
			return Result::Fail("negative cycle");
		}
	}
	return Result::Of(D);
}

Real OpsOctile(Real Dx, Real Dy)
{
	const Real A = gm::fabs(Dx), B = gm::fabs(Dy);
	return gm::max(A, B) + (gm::sqrt(2.0) - 1.0) * gm::min(A, B);
}

Result OctileFn(const Call& C) { return Result::Of(OpsOctile(C.N(0), C.N(1))); }
Result GridRatioFn(const Call& C) { return Result::Of(OpsOctile(gm::cos(C.N(0)), gm::sin(C.N(0)))); }
Result GridBoundFn(const Call& C)
{
	(void)C;
	return Result::Of(1.0 / gm::cos(cst::Pi() / 8.0));
}

Result LapFn(const Call& C)
{
	const Matrix& A = C.M(0);
	if (!A.IsSquare() || !OpsAllFinite(A))
	{
		return Result::Fail("needs a square finite adjacency matrix");
	}
	return Result::Of(OpsLaplacian(A));
}

bool OpsSpectrum(const Matrix& A, std::vector<Real>* W, Matrix* V, std::string* Error)
{
	if (!OpsAdjacency(A, Error))
	{
		return false;
	}
	if (A.Rows < 2)
	{
		*Error = "needs at least 2 nodes";
		return false;
	}
	return la::SymmetricEigen(OpsLaplacian(A), W, V, Error);
}

Result Lambda2Fn(const Call& C)
{
	std::vector<Real> W;
	std::string E;
	if (!OpsSpectrum(C.M(0), &W, nullptr, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(W[1]);
}

// x(t) = exp(-L t) x0 through the eigenvectors of L.
Result ConsensusFn(const Call& C)
{
	const Matrix &A = C.M(0), &X0 = C.M(1);
	const Real T = C.N(2);
	std::vector<Real> W;
	Matrix V;
	std::string E;
	if (!OpsSpectrum(A, &W, &V, &E))
	{
		return Result::Fail(E);
	}
	const int N = A.Rows;
	if (X0.Count() != N || !X0.IsVector())
	{
		return Result::Fail("x0 needs one value per node");
	}
	if (!(T >= 0))
	{
		return Result::Fail("time t must be >= 0");
	}
	std::vector<Real> Coef(static_cast<size_t>(N), 0.0), Xt(static_cast<size_t>(N), 0.0);
	for (int k = 0; k < N; ++k)
	{
		Real S = 0.0;
		for (int i = 0; i < N; ++i)
		{
			S += V(i, k) * X0.V[static_cast<size_t>(i)];
		}
		Coef[static_cast<size_t>(k)] = S * gm::exp(-W[static_cast<size_t>(k)] * T);
	}
	for (int i = 0; i < N; ++i)
	{
		Real S = 0.0;
		for (int k = 0; k < N; ++k)
		{
			S += V(i, k) * Coef[static_cast<size_t>(k)];
		}
		Xt[static_cast<size_t>(i)] = S;
	}
	auto Spread = [N](const std::vector<Real>& X) {
		const Real M = std::accumulate(X.begin(), X.end(), Real(0)) / N;
		Real S = 0.0;
		for (Real Xi : X)
		{
			S += (Xi - M) * (Xi - M);
		}
		return gm::sqrt(S);
	};
	const Real Mean = std::accumulate(X0.V.begin(), X0.V.end(), Real(0)) / N;
	return Result::Record({{"dis", Spread(Xt)}, {"bound", gm::exp(-W[1] * T) * Spread(X0.V)}, {"mean", Mean}});
}

// R_uv: ground v, inject one unit at u, read the voltage at u.
Result ReffFn(const Call& C)
{
	const Matrix& A = C.M(0);
	std::string E;
	int U = 0, V = 0;
	if (!OpsAdjacency(A, &E))
	{
		return Result::Fail(E);
	}
	const int N = A.Rows;
	if (!OpsSource(C, 1, N, &U) || !OpsSource(C, 2, N, &V))
	{
		return Result::Fail("nodes must be 1..n");
	}
	if (U == V)
	{
		return Result::Of(0.0);
	}
	const Matrix Lg = OpsDropRowCol(OpsLaplacian(A), V);
	const int Ug = U < V ? U : U - 1;
	Matrix B(N - 1, 1);
	B.V[static_cast<size_t>(Ug)] = 1.0;
	Matrix X;
	if (!la::Solve(Lg, B, &X, &E))
	{
		return Result::Fail("graph must be connected (the grounded Laplacian is singular)");
	}
	return Result::Of(X.V[static_cast<size_t>(Ug)]);
}

// Kirchhoff: any cofactor of the Laplacian counts (weighted) spanning trees.
Result SpanTreesFn(const Call& C)
{
	const Matrix& A = C.M(0);
	std::string E;
	if (!OpsAdjacency(A, &E))
	{
		return Result::Fail(E);
	}
	if (A.Rows == 1)
	{
		return Result::Of(1.0);
	}
	Real D = 0.0;
	if (!la::Determinant(OpsDropRowCol(OpsLaplacian(A), 0), &D, &E))
	{
		return Result::Of(0.0);
	}
	return Result::Of(D);
}

/* ---------- SLAM bench ---------- */
// slam.js lgamma (Lanczos, 6 terms) and gammp (series below a + 1, continued fraction above).
Real OpsLgamma(Real X)
{
	static const Real Cf[] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
	Real Y = X, T = X + 5.5;
	T -= (X + 0.5) * gm::log(T);
	Real S = 1.000000000190015;
	for (Real Ci : Cf)
	{
		Y += 1;
		S += Ci / Y;
	}
	return -T + gm::log(2.5066282746310005 * S / X);
}

Real OpsGammp(Real A, Real X)
{
	if (X <= 0)
	{
		return 0.0;
	}
	const Real Gl = OpsLgamma(A);
	if (X < A + 1)
	{
		Real Ap = A, Del = 1.0 / A, Sum = Del;
		for (int n = 0; n < 1000; ++n)
		{
			Ap += 1.0;
			Del *= X / Ap;
			Sum += Del;
			if (gm::fabs(Del) < gm::fabs(Sum) * 1e-16)
			{
				break;
			}
		}
		return Sum * gm::exp(-X + A * gm::log(X) - Gl);
	}
	Real B = X + 1 - A, Cc = 1e300, D = 1.0 / B, H = D;
	for (int i = 1; i < 1000; ++i)
	{
		const Real An = -i * (i - A);
		B += 2;
		D = An * D + B;
		if (gm::fabs(D) < 1e-300)
		{
			D = 1e-300;
		}
		Cc = B + An / Cc;
		if (gm::fabs(Cc) < 1e-300)
		{
			Cc = 1e-300;
		}
		D = 1.0 / D;
		const Real Del = D * Cc;
		H *= Del;
		if (gm::fabs(Del - 1) < 1e-16)
		{
			break;
		}
	}
	return 1.0 - gm::exp(-X + A * gm::log(X) - Gl) * H;
}

Real OpsChi2Cdf(Real X, Real K) { return OpsGammp(K / 2.0, X / 2.0); }

Real OpsChi2Inv(Real P, Real K)
{
	Real Lo = 0.0, Hi = K + 40.0 * gm::sqrt(2.0 * K) + 40.0;
	for (int i = 0; i < 200; ++i)
	{
		const Real M = (Lo + Hi) / 2.0;
		if (OpsChi2Cdf(M, K) < P)
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	return (Lo + Hi) / 2.0;
}

Result GammpFn(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("shape a must be > 0");
	}
	return Result::Of(OpsGammp(C.N(0), C.N(1)));
}

Result Chi2CdfFn(const Call& C)
{
	if (!(C.N(1) > 0))
	{
		return Result::Fail("degrees of freedom k must be > 0");
	}
	return Result::Of(OpsChi2Cdf(C.N(0), C.N(1)));
}

Result Chi2InvFn(const Call& C)
{
	const Real P = C.N(0), K = C.N(1);
	if (!(P >= 0 && P < 1))
	{
		return Result::Fail("probability p must be in [0, 1)");
	}
	if (!(K > 0 && K <= 1e6))
	{
		return Result::Fail("degrees of freedom k must be in (0, 1e6]");
	}
	return Result::Of(OpsChi2Inv(P, K));
}

Result Chi2GateFn(const Call& C)
{
	const Real P = C.N(0), K = C.N(1);
	if (!(P >= 0 && P < 1))
	{
		return Result::Fail("probability p must be in [0, 1)");
	}
	if (!(K > 0 && K <= 1e6))
	{
		return Result::Fail("degrees of freedom k must be in (0, 1e6]");
	}
	return Result::Of(K == 2.0 ? -2.0 * gm::log(1.0 - P) : OpsChi2Inv(P, K));
}

Result MahalFn(const Call& C)
{
	const Matrix &Nu = C.M(0), &S = C.M(1);
	if (!Nu.IsVector() || !S.IsSquare() || S.Rows != Nu.Count())
	{
		return Result::Fail("needs an n-vector and an n x n covariance");
	}
	Matrix X;
	std::string E;
	if (!la::Solve(S, Matrix::Column(Nu.V), &X, &E))
	{
		return Result::Fail("covariance is singular");
	}
	Real D2 = 0.0;
	for (size_t i = 0; i < Nu.V.size(); ++i)
	{
		D2 += Nu.V[i] * X.V[i];
	}
	return Result::Of(D2);
}

// slam.js svd2, mm2, tr2, det2 and kabsch, operation for operation.
struct OpsM2
{
	Real A[2][2];
};

OpsM2 OpsMm2(const OpsM2& X, const OpsM2& Y)
{
	OpsM2 R;
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			R.A[i][j] = X.A[i][0] * Y.A[0][j] + X.A[i][1] * Y.A[1][j];
		}
	}
	return R;
}
OpsM2 OpsTr2(const OpsM2& X) { return OpsM2{{{X.A[0][0], X.A[1][0]}, {X.A[0][1], X.A[1][1]}}}; }
Real OpsDet2(const OpsM2& M) { return M.A[0][0] * M.A[1][1] - M.A[0][1] * M.A[1][0]; }

void OpsSvd2(const OpsM2& Mt, OpsM2* U, OpsM2* V)
{
	const Real(&M)[2][2] = Mt.A;
	const Real P = M[0][0] * M[0][0] + M[1][0] * M[1][0], Q = M[0][0] * M[0][1] + M[1][0] * M[1][1], W = M[0][1] * M[0][1] + M[1][1] * M[1][1];
	const Real Ang = 0.5 * gm::atan2(2 * Q, P - W), C = gm::cos(Ang), S = gm::sin(Ang);
	*V = OpsM2{{{C, -S}, {S, C}}};
	const Real Rad = gm::hypot((P - W) / 2, Q), L1 = (P + W) / 2 + Rad, L2 = gm::max(Real(0), (P + W) / 2 - Rad), S1 = gm::sqrt(L1), S2 = gm::sqrt(L2);
	const Real Mv1[2] = {M[0][0] * C + M[0][1] * S, M[1][0] * C + M[1][1] * S}, Mv2[2] = {-M[0][0] * S + M[0][1] * C, -M[1][0] * S + M[1][1] * C};
	Real U1[2] = {1, 0}, U2[2];
	if (S1 > 0)
	{
		U1[0] = Mv1[0] / S1;
		U1[1] = Mv1[1] / S1;
	}
	if (S2 > 1e-300 * (1 + S1))
	{
		U2[0] = Mv2[0] / S2;
		U2[1] = Mv2[1] / S2;
	}
	else
	{
		U2[0] = -U1[1];
		U2[1] = U1[0];
	}
	const Real Nn = gm::hypot(U2[0], U2[1]);
	U2[0] /= Nn;
	U2[1] /= Nn;
	*U = OpsM2{{{U1[0], U2[0]}, {U1[1], U2[1]}}};
}

struct OpsAlign
{
	Real Th = 0.0, Tx = 0.0, Ty = 0.0, S = 0.0;
};

// R, t minimizing sum |Q_i - (R P_i + t)|^2, with Umeyama's reflection guard.
OpsAlign OpsKabsch(const Matrix& P, const Matrix& Q)
{
	const int N = P.Rows;
	Real Pm[2] = {0, 0}, Qm[2] = {0, 0};
	for (int i = 0; i < N; ++i)
	{
		Pm[0] += P(i, 0) / N;
		Pm[1] += P(i, 1) / N;
		Qm[0] += Q(i, 0) / N;
		Qm[1] += Q(i, 1) / N;
	}
	OpsM2 H{{{0, 0}, {0, 0}}};
	Real S = 0.0;
	for (int i = 0; i < N; ++i)
	{
		const Real A0 = P(i, 0) - Pm[0], A1 = P(i, 1) - Pm[1], B0 = Q(i, 0) - Qm[0], B1 = Q(i, 1) - Qm[1];
		H.A[0][0] += A0 * B0;
		H.A[0][1] += A0 * B1;
		H.A[1][0] += A1 * B0;
		H.A[1][1] += A1 * B1;
		S += A0 * A0 + A1 * A1;
	}
	OpsM2 U, V;
	OpsSvd2(H, &U, &V);
	const Real Dt = OpsDet2(OpsMm2(V, OpsTr2(U)));
	const Real D = Dt > 0 ? 1.0 : (Dt < 0 ? -1.0 : 1.0);
	const OpsM2 R = OpsMm2(OpsMm2(V, OpsM2{{{1, 0}, {0, D}}}), OpsTr2(U));
	OpsAlign Out;
	Out.Tx = Qm[0] - R.A[0][0] * Pm[0] - R.A[0][1] * Pm[1];
	Out.Ty = Qm[1] - R.A[1][0] * Pm[0] - R.A[1][1] * Pm[1];
	Out.Th = gm::atan2(R.A[1][0], R.A[0][0]);
	Out.S = S;
	return Out;
}

Result KabschFn(const Call& C)
{
	const Matrix &P = C.M(0), &Q = C.M(1);
	if (P.Cols != 2 || Q.Cols != 2 || P.Rows != Q.Rows || P.Rows < 2)
	{
		return Result::Fail("needs two n x 2 point sets with the same n >= 2");
	}
	if (!OpsAllFinite(P) || !OpsAllFinite(Q))
	{
		return Result::Fail("points must be finite");
	}
	const OpsAlign K = OpsKabsch(P, Q);
	const Real Cs = gm::cos(K.Th), Sn = gm::sin(K.Th);
	Real Cost = 0.0;
	for (int i = 0; i < P.Rows; ++i)
	{
		const Real Ex = Cs * P(i, 0) - Sn * P(i, 1) + K.Tx - Q(i, 0), Ey = Sn * P(i, 0) + Cs * P(i, 1) + K.Ty - Q(i, 1);
		Cost += Ex * Ex + Ey * Ey;
	}
	return Result::Record({{"th", K.Th}, {"tx", K.Tx}, {"ty", K.Ty}, {"cost", Cost}, {"S", K.S}});
}

Result HeadVarFn(const Call& C)
{
	if (!(C.N(1) > 0))
	{
		return Result::Fail("spread S must be > 0");
	}
	return Result::Of(2.0 * C.N(0) * C.N(0) / C.N(1));
}

Result DriftFn(const Call& C)
{
	const Real Sig = C.N(0), S = C.N(1), N = C.N(2);
	if (!(S > 0 && N >= 0))
	{
		return Result::Fail("needs spread S > 0 and frames N >= 0");
	}
	const Real Var = N * 2.0 * Sig * Sig / S;
	return Result::Record({{"var", Var}, {"rms", gm::sqrt(Var)}});
}

// mulberry32 with Box-Muller, cos first then sin (slam.js rng).
struct OpsRng
{
	uint32_t A;
	bool HasSpare = false;
	Real Spare = 0.0;

	explicit OpsRng(uint32_t Seed) : A(Seed) {}
	Real U()
	{
		A += 0x6D2B79F5u;
		uint32_t T = (A ^ (A >> 15)) * (1u | A);
		T = (T + ((T ^ (T >> 7)) * (61u | T))) ^ T;
		static const Real Two32 = Real(4294967296ull);
		return Real(T ^ (T >> 14)) / Two32;
	}
	Real G()
	{
		if (HasSpare)
		{
			HasSpare = false;
			return Spare;
		}
		static const Real TwoPi = Real(2) * cst::Pi();
		const Real U1 = U(), U2 = U(), Rr = gm::sqrt(Real(-2) * gm::log(Real(1) - U1)), Th = TwoPi * U2;
		Spare = Rr * gm::sin(Th);
		HasSpare = true;
		return Rr * gm::cos(Th);
	}
};

// One frame pair of the downward camera (slam.js pairFit): f = 500 px, h = 2 m, 640 x 480, 0.25 m forward.
OpsAlign OpsPairFit(OpsRng& R, int N, Real Sm, Real C, Real S)
{
	static const Real Wm = Real(640) * Real(2) / Real(500), Hm = Real(480) * Real(2) / Real(500), Fwd = 0.25, Half = 0.5;
	Real Pa0 = 0, Pa1 = 0, Pb0 = 0, Pb1 = 0, H00 = 0, H01 = 0, H10 = 0, H11 = 0, Bb = 0;
	for (int i = 0; i < N; ++i)
	{
		const Real Ax = (R.U() - Half) * Wm, Ay = (R.U() - Half) * Hm, Dx = Ax - Fwd, Bx = C * Dx + S * Ay, By = -S * Dx + C * Ay;
		const Real A0 = Ax + Sm * R.G(), A1 = Ay + Sm * R.G(), B0 = Bx + Sm * R.G(), B1 = By + Sm * R.G();
		Pa0 += A0;
		Pa1 += A1;
		Pb0 += B0;
		Pb1 += B1;
		H00 += B0 * A0;
		H01 += B0 * A1;
		H10 += B1 * A0;
		H11 += B1 * A1;
		Bb += B0 * B0 + B1 * B1;
	}
	Pa0 /= N;
	Pa1 /= N;
	Pb0 /= N;
	Pb1 /= N;
	const OpsM2 H{{{H00 - N * Pb0 * Pa0, H01 - N * Pb0 * Pa1}, {H10 - N * Pb1 * Pa0, H11 - N * Pb1 * Pa1}}};
	OpsM2 U, V;
	OpsSvd2(H, &U, &V);
	const Real Dt = OpsDet2(OpsMm2(V, OpsTr2(U)));
	const Real D = Dt > 0 ? 1.0 : (Dt < 0 ? -1.0 : 1.0);
	const OpsM2 Rm = OpsMm2(OpsMm2(V, OpsM2{{{1, 0}, {0, D}}}), OpsTr2(U));
	OpsAlign Out;
	Out.Th = gm::atan2(Rm.A[1][0], Rm.A[0][0]);
	Out.Tx = Pa0 - Rm.A[0][0] * Pb0 - Rm.A[0][1] * Pb1;
	Out.Ty = Pa1 - Rm.A[1][0] * Pb0 - Rm.A[1][1] * Pb1;
	Out.S = Bb - N * (Pb0 * Pb0 + Pb1 * Pb1);
	return Out;
}

// Monte Carlo of chained frame-to-frame odometry (slam.js voWalk), read at the last frame.
Result VoWalkFn(const Call& C)
{
	int N = 0, Nf = 0, M = 0, Seed = 0;
	if (!OpsWhole(C.N(0), 1, 1000, &N) || !OpsWhole(C.N(1), 2, 200, &Nf) || !OpsWhole(C.N(4), 1, 2000, &M) || !OpsWhole(C.N(5), 0, 2147483647, &Seed))
	{
		return Result::Fail("needs whole N 1..1000, n 2..200, M 1..2000 and seed >= 0");
	}
	if (static_cast<long long>(N) * Nf * M > 40000)
	{
		return Result::Fail("N * n * M must be at most 40000 on the gears (about 2 s); the bench runs M = 200");
	}
	const Real Sm = C.N(2) * Real(2) / Real(500), Th = C.N(3) * (cst::Pi() / Real(180)), C0 = gm::cos(Th), S0 = gm::sin(Th), Fwd = 0.25, Two = 2;
	Real Sh = 0.0, Sp = 0.0, Sv = 0.0;
	for (int m = 0; m < M; ++m)
	{
		OpsRng R(static_cast<uint32_t>(Seed) + static_cast<uint32_t>(m) * 0x9E3779B9u);
		Real X = 0, Y = 0, Hd = 0, Xt = 0, Yt = 0, Ph = 0, V = 0, Eh = 0;
		for (int k = 0; k < N; ++k)
		{
			const OpsAlign K = OpsPairFit(R, Nf, Sm, C0, S0);
			const Real Ch = gm::cos(Hd), Sn = gm::sin(Hd);
			X += Ch * K.Tx - Sn * K.Ty;
			Y += Sn * K.Tx + Ch * K.Ty;
			Hd += K.Th;
			Xt += gm::cos(Ph) * Fwd;
			Yt += gm::sin(Ph) * Fwd;
			Ph += Th;
			V += Two * Sm * Sm / K.S;
			Eh = Hd - Ph;
		}
		Sh += Eh * Eh / M;
		Sp += ((X - Xt) * (X - Xt) + (Y - Yt) * (Y - Yt)) / M;
		Sv += V / M;
	}
	return Result::Record({{"rmsh", gm::sqrt(Sh)}, {"pred", gm::sqrt(Sv)}, {"rmsp", gm::sqrt(Sp)}});
}

Result NeesBandFn(const Call& C)
{
	const Real M = C.N(0), K = C.N(1), P = C.N(2);
	if (!(M >= 1 && K > 0 && M * K <= 1e6))
	{
		return Result::Fail("needs runs M >= 1, dimension k > 0, M k <= 1e6");
	}
	if (!(P > 0 && P < 1))
	{
		return Result::Fail("confidence p must be in (0, 1)");
	}
	const Real Lo = P == 0.95 ? 0.025 : (1.0 - P) / 2.0, Hi = P == 0.95 ? 0.975 : (1.0 + P) / 2.0;
	return Result::Record({{"lo", OpsChi2Inv(Lo, K * M) / M}, {"hi", OpsChi2Inv(Hi, K * M) / M}});
}
} // namespace ops_detail

void RegisterOpsLibrary(Registry& R)
{
	using namespace ops_detail;
	const char* Q = "queues";
	const char* I = "inventory";
	const char* Nw = "networks";
	const char* G = "graphs";
	const char* Sl = "slam";

	Def(R, "erlangb", Q, "ch13 Thm 7 (Erlang 1917), M/M/c/c loss", "Erlang B: probability a call is blocked with c lines and offered load a", "probability", {Num("c", "servers"), Num("a", "erlang")}, &ErlangBFn);
	Def(R, "erlangc", Q, "logistics/ §6 erlangC; ch13 Thm 7", "Erlang C: probability an arrival waits in M/M/c (1 if a >= c)", "probability", {Num("c", "servers"), Num("a", "erlang")}, &ErlangCFn);
	Def(R, "mmc", Q, "logistics/ §6 mmc; ch13 Thm 7", "M/M/c steady state: load, Erlang C, waits and queue lengths", "record rho, C, Wq (s), W (s), L, Lq",
		{Num("lambda", "1/s", 2.4), Num("mu", "1/s", 1), Num("c", "servers", 3)}, &MmcFn);
	Def(R, "mmctail", Q, "logistics/ §d waiting-time tail", "M/M/c tail P(Wq > t) = C exp(-(c mu - lambda) t)", "probability",
		{Num("lambda", "1/s", 2.4), Num("mu", "1/s", 1), Num("c", "servers", 3), Num("t", "s", 1)}, &MmcTailFn);
	Def(R, "little", Q, "ch13 Thm 6 Little's law", "Little's law L = lambda W", "mean number in system", {Num("lambda", "1/s"), Num("W", "s")}, &LittleFn);
	Def(R, "pk", Q, "logistics/ §6 pk; ch13 Thm 7 Pollaczek-Khinchine", "M/G/1 mean wait, E[S^2] = (1 + cs2)/mu^2 (exp cs2 = 1, det 0, uniform 1/3)", "record rho, Wq (s), W (s), L, Lq",
		{Num("lambda", "1/s", 0.8), Num("mu", "1/s", 1), Num("cs2", "-", 1)}, &PkFn);
	Def(R, "kingman", Q, "ch13 Thm 7 Kingman 1962", "G/G/1: Kingman's upper bound on Wq and his heavy-traffic approximation", "record bound (s), approx (s)",
		{Num("lambda", "1/s", 0.8), Num("mu", "1/s", 1), Num("ca2", "-", 1), Num("cs2", "-", 1)}, &KingmanFn);

	Def(R, "eoq", I, "logistics/ §7 eoq; ch13 Thm 8 Harris", "economic order quantity sqrt(2 K D / h)", "units per order",
		{Num("K", "cost/order", 50), Num("D", "units/yr", 1200), Num("h", "cost/unit/yr", 3)}, &EoqFn);
	Def(R, "eoqcost", I, "logistics/ §7 eoqCost; ch13 Thm 8", "annual ordering plus holding cost K D / Q + h Q / 2", "cost/yr",
		{Num("Q", "units"), Num("K", "cost/order", 50), Num("D", "units/yr", 1200), Num("h", "cost/unit/yr", 3)}, &EoqCostFn);
	Def(R, "fractile", I, "ch13 Thm 8 critical fractile", "newsvendor critical fractile (r - c)/(r - v)", "probability",
		{Num("price", "cost/unit", 10), Num("cost", "cost/unit", 4), Num("salv", "cost/unit", 1)}, &FractileFn);
	Def(R, "newsvendor", I, "logistics/ §7 newsvendor", "one-season order for Poisson demand: fractile, best order, its profit, brute-force check", "record crit, q, profit, qbest",
		{Num("price", "cost/unit", 10), Num("cost", "cost/unit", 4), Num("salv", "cost/unit", 1), Num("mean", "units", 40)}, &NewsvendorFn);
	Def(R, "nvprofit", I, "logistics/ §7 newsvendor profit", "expected newsvendor profit of ordering q with Poisson demand", "cost",
		{Num("q", "units"), Num("price", "cost/unit", 10), Num("cost", "cost/unit", 4), Num("salv", "cost/unit", 1), Num("mean", "units", 40)}, &NvProfitFn);
	Def(R, "ss", I, "logistics/ §7 ssDP; ch13 Thm 9 Scarf", "optimal (s, S) by value iteration, demand Binomial(10, 0.4), backlog", "record s, S, v0 (discounted cost from stock 0)",
		{Num("K", "cost/order", 32), Num("p", "cost/unit short", 9), Num("h", "cost/unit held", 1), Num("c", "cost/unit", 1), Num("alpha", "-", 0.9)}, &SsFn);
	Def(R, "sscost", I, "logistics/ §7 ssEval", "exact discounted cost of a fixed (s, S) policy from stock x", "cost",
		{Num("s", "units"), Num("S", "units"), Num("K", "cost/order", 32), Num("p", "cost/unit short", 9), Num("h", "cost/unit held", 1), Num("c", "cost/unit", 1), Num("alpha", "-", 0.9), Num("x", "units", 0)},
		&SsCostFn);

	Def(R, "maxflow", Nw, "logistics/ §1 maxflow (Edmonds-Karp); ch13 Thm 1", "max flow on a capacity matrix; t = 0 means the last node", "record flow, aug (paths), ns (nodes on the source side of the min cut)",
		{Mat("C", "capacity"), Num("s", "node", 1), Num("t", "node", 0)}, &MaxFlowFn);
	Def(R, "flowmat", Nw, "logistics/ §1 maxflow", "the maximum flow itself, as a matrix of lane flows", "matrix",
		{Mat("C", "capacity"), Num("s", "node", 1), Num("t", "node", 0)}, &FlowMatFn);
	Def(R, "transport", Nw, "logistics/ §3 transport; ch13 Thm 3", "transportation problem by the stepping-stone simplex: optimal cost", "cost",
		{Mat("sup", "units"), Mat("dem", "units"), Mat("C", "cost/unit")}, &TransportFn);
	Def(R, "tplan", Nw, "logistics/ §3 transport", "optimal transportation plan (supplies x demands)", "matrix of units",
		{Mat("sup", "units"), Mat("dem", "units"), Mat("C", "cost/unit")}, &TPlanFn);
	Def(R, "tprices", Nw, "logistics/ §3 transport duals", "dual prices [u, v] (u1 = 0) certifying the optimal plan", "row vector u then v",
		{Mat("sup", "units"), Mat("dem", "units"), Mat("C", "cost/unit")}, &TPricesFn);
	Def(R, "assign", Nw, "ch13 guise: data association; ch09 Kuhn 1955", "assignment problem (Hungarian method): least total cost", "cost", {Mat("C", "cost")}, &AssignFn);
	Def(R, "assignperm", Nw, "ch13 guise: data association", "optimal assignment: the column chosen for each row", "row vector of columns", {Mat("C", "cost")}, &AssignPermFn);
	Def(R, "karp", Nw, "logistics/ §5 karp; ch13 Thm 5", "maximum cycle mean (max-plus eigenvalue) of an irreducible matrix, -inf = no arc", "cycle time", {Mat("A", "time")}, &KarpFn);
	Def(R, "mpmul", Nw, "logistics/ §5 otimes", "max-plus product (A x B)_ij = max_k A_ik + B_kj", "matrix", {Mat("A", "time"), Mat("B", "time")}, &MpMulFn);
	Def(R, "mppow", Nw, "logistics/ §5 mpRun", "max-plus power A^k (k = 0 gives the max-plus identity)", "matrix", {Mat("A", "time"), Num("k", "-")}, &MpPowFn);
	Def(R, "w1", Nw, "logistics/ §4 w1cdf; ch13 Thm 4", "Wasserstein-1 distance between two 1-D samples, integral |F - G|", "distance", {Mat("x", ""), Mat("y", "")}, &W1Fn);
	Def(R, "w1w", Nw, "logistics/ §4 w1cdf", "W1 = integral |F - G| for points with masses (cost of moving the masses)", "mass x distance",
		{Mat("x", ""), Mat("a", "mass"), Mat("y", ""), Mat("b", "mass")}, &W1wFn);
	Def(R, "otmono", Nw, "logistics/ §4 otCase monotone coupling; ch13 Thm 4", "cost of the monotone coupling on a line with cost |x - y|^p (optimal for p >= 1)", "cost",
		{Mat("x", ""), Mat("a", "mass"), Mat("y", ""), Mat("b", "mass"), Num("p", "-", 1)}, &OtMonoFn);
	Def(R, "mst", Nw, "logistics/ §9 kruskal; ch13 Thm 11", "minimum spanning tree weight (Kruskal) of a symmetric weight matrix, inf = no edge", "weight", {Mat("W", "weight")}, &MstFn);
	Def(R, "distmat", Nw, "logistics/ §8 dmat", "Euclidean distance matrix of an n x 2 point set", "matrix", {Mat("P", "m")}, &DistMatFn);
	Def(R, "tsp", Nw, "ch13 Thm 10 (Held-Karp dynamic programming)", "shortest closed tour through every stop, exact, n <= 15", "length", {Mat("D", "distance")}, &TspFn);
	Def(R, "hkbound", Nw, "logistics/ §8 heldKarp; ch13 Thm 10", "Held-Karp 1-tree lower bound on the tour by subgradient ascent", "record lb, onetree",
		{Mat("D", "distance"), Num("iters", "-", 200)}, &HkBoundFn);
	Def(R, "christofides", Nw, "logistics/ §8 christofides; ch13 Thm 10", "Christofides tour (MST + exact matching + Euler shortcut), <= 1.5 OPT", "record len, mst, match", {Mat("D", "distance")}, &ChristofidesFn);
	Def(R, "tourlen", Nw, "logistics/ §8 tourLen", "length of a closed tour given as a vector of stops 1..n", "length", {Mat("D", "distance"), Mat("tour", "stops")}, &TourLenFn);

	Def(R, "dijkstra", G, "ch09 Thm 1 Dijkstra", "shortest distances from s on a cost matrix (inf = no edge, costs >= 0)", "row vector of distances", {Mat("W", "cost"), Num("s", "node", 1)}, &DijkstraFn);
	Def(R, "bellman", G, "ch09 counterexample: Bellman-Ford", "shortest distances from s allowing negative costs; fails on a negative cycle", "row vector of distances", {Mat("W", "cost"), Num("s", "node", 1)}, &BellmanFn);
	Def(R, "floyd", G, "ch09 shortest paths (Floyd-Warshall)", "all-pairs shortest distances", "matrix of distances", {Mat("W", "cost")}, &FloydFn);
	Def(R, "octile", G, "ch09 guise: grid planners", "octile distance max + (sqrt 2 - 1) min, the 8-connected grid metric", "cells", {Num("dx", "cells"), Num("dy", "cells")}, &OctileFn);
	Def(R, "gridratio", G, "ch09 guise: grid planners", "8-connected path length over straight-line length at heading theta", "ratio", {Num("theta", "rad")}, &GridRatioFn);
	Def(R, "gridbound", G, "ch09 guise: grid planners", "worst grid-path overestimate 1/cos(pi/8), at 22.5 deg", "ratio", {}, &GridBoundFn);
	Def(R, "lap", G, "ch09 Def graph Laplacian", "Laplacian L = D - A of a weighted adjacency matrix", "matrix", {Mat("A", "weight")}, &LapFn);
	Def(R, "lambda2", G, "ch09 Thm 4 (Fiedler 1973)", "algebraic connectivity: the consensus rate of xdot = -L x", "1/s", {Mat("A", "weight")}, &Lambda2Fn);
	Def(R, "consensus", G, "ch09 Thm 4 consensus", "disagreement |x(t) - mean| of xdot = -L x against the bound exp(-lambda2 t)|x0 - mean|", "record dis, bound, mean",
		{Mat("A", "weight"), Mat("x0", ""), Num("t", "s")}, &ConsensusFn);
	Def(R, "reff", G, "ch09 Thm 5 effective resistance", "effective resistance between u and v (conductance matrix); pose variance in a pose graph", "ohm (variance)",
		{Mat("A", "conductance"), Num("u", "node"), Num("v", "node", 1)}, &ReffFn);
	Def(R, "spantrees", G, "ch09 Thm 6 Kirchhoff; logistics/ §9", "number of spanning trees: any cofactor of the Laplacian", "count", {Mat("A", "weight")}, &SpanTreesFn);

	Def(R, "gammp", Sl, "slam/ gammp", "regularized lower incomplete gamma P(a, x)", "probability", {Num("a", "-"), Num("x", "-")}, &GammpFn);
	Def(R, "chi2cdf", Sl, "slam/ chi2cdf", "chi-square CDF with k degrees of freedom", "probability", {Num("x", "-"), Num("k", "dof")}, &Chi2CdfFn);
	Def(R, "chi2inv", Sl, "slam/ chi2inv", "chi-square quantile by bisection", "-", {Num("p", "probability"), Num("k", "dof")}, &Chi2InvFn);
	Def(R, "chi2gate", Sl, "slam/ §6 gate2", "data-association gate gamma: chi-square quantile (k = 2: -2 ln(1 - p))", "-", {Num("p", "probability", 0.95), Num("k", "dof", 2)}, &Chi2GateFn);
	Def(R, "mahal", Sl, "slam/ §6 data association", "squared Mahalanobis distance d^2 = nu' S^-1 nu, compared with chi2gate", "-", {Mat("nu", ""), Mat("S", "")}, &MahalFn);
	Def(R, "nees", Sl, "slam/ §4 NEES", "normalized estimation error squared e' P^-1 e", "-", {Mat("e", ""), Mat("P", "")}, &MahalFn);
	Def(R, "neesband", Sl, "slam/ §4 neesBand", "two-sided band for the M-run average NEES of a k-dimensional state", "record lo, hi",
		{Num("M", "runs", 40), Num("k", "dof", 3), Num("p", "probability", 0.95)}, &NeesBandFn);
	Def(R, "kabsch", Sl, "slam/ §1 kabsch (Umeyama reflection guard)", "2-D rigid alignment: R, t minimizing sum |Q - (R P + t)|^2 for n x 2 point sets", "record th (rad), tx, ty, cost, S (spread of P)",
		{Mat("P", "m"), Mat("Q", "m")}, &KabschFn);
	Def(R, "headvar", Sl, "slam/ §1 headingVar", "first-order heading variance of one fit, 2 sigma^2 / S", "rad^2", {Num("sigma", "m"), Num("S", "m^2")}, &HeadVarFn);
	Def(R, "drift", Sl, "slam/ §2 drift, sqrt(N) law", "heading drift after N chained frames: variance N 2 sigma^2 / S and its RMS", "record var (rad^2), rms (rad)",
		{Num("sigma", "m"), Num("S", "m^2"), Num("N", "frames")}, &DriftFn);
	Def(R, "vowalk", Sl, "slam/ §2 voWalk", "Monte Carlo of chained visual odometry (the bench's seeded generator), at frame N; N n M <= 40000 (bench M = 200 is too slow on the gears)", "record rmsh (rad), pred (rad), rmsp (m)",
		{Num("N", "frames", 100), Num("n", "features", 20), Num("sig", "px", 1), Num("yaw", "deg/frame", 3), Num("M", "runs", 5), Num("seed", "-", 2)}, &VoWalkFn);
}
} // namespace mcalc
