// Flight: the equation sheet (flight-model/equations.js), the NASA Glenn atmosphere, flight model and cloud
// chain (flight-model/model.js), the airlifter planform bench (flight-model/planform.js) and the par values of
// the game (game/game.js). Every routine follows its module step for step on the gear machine's Real.
#include "CalcLibrary.h"
#include "GearMath.h"

#include <algorithm>
#include <vector>

namespace mcalc
{
namespace flight_detail
{
const Real& FlG0()
{
	static const Real V(9.80665);
	return V;
}
Real FlNan() { return Real::NaN(); }
Real FlInf() { return Real::Infinity(); }
Real FlRad(const Real& D) { return D * cst::Pi() / 180; }
Real JsRound(const Real& X) { return gm::floor(X + Real(0.5)); } // Math.round
Real Clamp(const Real& X, const Real& A, const Real& B) { return gm::min(B, gm::max(A, X)); }
Real Flag(bool B) { return B ? Real(1) : Real(0); }

template <typename F>
Real Trapz(F Fn, const Real& A, const Real& B, int N)
{
	const Real H = (B - A) / N;
	Real S = (Fn(A) + Fn(B)) / 2;
	for (int i = 1; i < N; ++i)
	{
		S += Fn(A + Real(i) * H);
	}
	return S * H;
}

// Gaussian elimination with partial pivoting on an n x (n + 1) augmented matrix, as equations.js's solve().
bool GaussSolve(std::vector<Real> M, int N, std::vector<Real>* X)
{
	const int W = N + 1;
	auto At = [&M, W](int R, int C) -> Real& { return M[static_cast<size_t>(R) * static_cast<size_t>(W) + static_cast<size_t>(C)]; };
	for (int c = 0; c < N; ++c)
	{
		int P = c;
		for (int r = c + 1; r < N; ++r)
		{
			if (gm::fabs(At(r, c)) > gm::fabs(At(P, c)))
			{
				P = r;
			}
		}
		if (P != c)
		{
			for (int j = 0; j < W; ++j)
			{
				std::swap(At(c, j), At(P, j));
			}
		}
		if (At(c, c).IsZero())
		{
			return false;
		}
		const Real Pivot = At(c, c);
		for (int r = c + 1; r < N; ++r)
		{
			const Real F = At(r, c) / Pivot;
			if (!F.IsZero())
			{
				for (int j = c; j <= N; ++j)
				{
					At(r, j) -= F * At(c, j);
				}
			}
		}
	}
	X->assign(static_cast<size_t>(N), Real(0));
	for (int r = N - 1; r >= 0; --r)
	{
		Real S = At(r, N);
		for (int j = r + 1; j < N; ++j)
		{
			S -= At(r, j) * (*X)[static_cast<size_t>(j)];
		}
		(*X)[static_cast<size_t>(r)] = S / At(r, r);
	}
	return true;
}

/* ---------------------------------------------------------------- equation sheet: aerodynamics */
// NACA four-digit mean-line slope: max camber m at p.
Real CamberSlope(const Real& M, const Real& P, const Real& X)
{
	if (P <= 0 || M.IsZero())
	{
		return Real(0);
	}
	return X < P ? 2 * M / (P * P) * (P - X) : 2 * M / ((1 - P) * (1 - P)) * (P - X);
}

const char* CheckCamber(const Real& M, const Real& P)
{
	if (M < 0 || M > Real(0.2))
	{
		return "max camber m is a fraction of chord, 0 to 0.2";
	}
	if (!M.IsZero() && (P <= 0 || P >= 1))
	{
		return "camber position p must lie strictly between 0 and 1";
	}
	return nullptr;
}

Result ThinAirfoil(const Call& C)
{
	const Real M = C.N(0), P = C.N(1), Alpha = C.N(2);
	if (const char* E = CheckCamber(M, P))
	{
		return Result::Fail(E);
	}
	// The three Glauert integrals by the trapezoidal rule on 20 000 intervals, sharing cos(theta) per node.
	const int N = 20000;
	const Real H = cst::Pi() / N;
	Real S1(0), S2(0), S3(0);
	for (int i = 0; i <= N; ++i)
	{
		const Real Cs = i == N ? Real(-1) : gm::cos(Real(i) * H), Sl = CamberSlope(M, P, (1 - Cs) / 2);
		const Real W = (i == 0 || i == N) ? Real(0.5) : Real(1);
		S1 += W * (Sl * (Cs - 1));
		S2 += W * (Sl * Cs);
		S3 += W * (Sl * (2 * Cs * Cs - 1));
	}
	const Real AL0 = -(S1 * H) / cst::Pi(), A1 = 2 / cst::Pi() * (S2 * H), A2 = 2 / cst::Pi() * (S3 * H);
	return Result::Record({{"aL0", AL0}, {"cl", 2 * cst::Pi() * (Alpha - AL0)}, {"cm4", cst::Pi() / 4 * (A2 - A1)}, {"A1", A1}, {"A2", A2}});
}

Result VortexPanel(const Call& C)
{
	const Real M = C.N(0), P = C.N(1), Alpha = C.N(2);
	const int N = C.Int(3);
	if (const char* E = CheckCamber(M, P))
	{
		return Result::Fail(E);
	}
	if (N < 2 || N > 200)
	{
		return Result::Fail("panels N must be 2 to 200");
	}
	std::vector<Real> Xe(static_cast<size_t>(N + 1)), Xv, Xc;
	for (int i = 0; i <= N; ++i)
	{
		Xe[static_cast<size_t>(i)] = (1 - gm::cos(cst::Pi() * i / N)) / 2;
	}
	for (int i = 0; i < N; ++i)
	{
		const Real D = Xe[static_cast<size_t>(i + 1)] - Xe[static_cast<size_t>(i)];
		Xv.push_back(Xe[static_cast<size_t>(i)] + Real(0.25) * D);
		Xc.push_back(Xe[static_cast<size_t>(i)] + Real(0.75) * D);
	}
	const Real TwoPi = 2 * cst::Pi();
	std::vector<Real> A(static_cast<size_t>(N) * static_cast<size_t>(N + 1));
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			A[static_cast<size_t>(i) * static_cast<size_t>(N + 1) + static_cast<size_t>(j)] = 1 / (TwoPi * (Xc[static_cast<size_t>(i)] - Xv[static_cast<size_t>(j)]));
		}
		A[static_cast<size_t>(i) * static_cast<size_t>(N + 1) + static_cast<size_t>(N)] = Alpha - CamberSlope(M, P, Xc[static_cast<size_t>(i)]);
	}
	std::vector<Real> G;
	if (!GaussSolve(A, N, &G))
	{
		return Result::Fail("singular vortex system");
	}
	Real Cl(0), Cm4(0);
	for (int j = 0; j < N; ++j)
	{
		Cl += 2 * G[static_cast<size_t>(j)];
		Cm4 -= 2 * G[static_cast<size_t>(j)] * (Xv[static_cast<size_t>(j)] - Real(0.25));
	}
	return Result::Record({{"cl", Cl}, {"cm4", Cm4}});
}

Result PrandtlGlauert(const Call& C)
{
	const Real C0 = C.N(0), M = C.N(1);
	if (M < 0 || M >= 1)
	{
		return Result::Fail("Prandtl-Glauert needs subsonic Mach, 0 <= M < 1");
	}
	return Result::Of(C0 / gm::sqrt(1 - M * M));
}

Result KarmanTsien(const Call& C)
{
	const Real Cp0 = C.N(0), M = C.N(1);
	if (M < 0 || M >= 1)
	{
		return Result::Fail("Karman-Tsien needs subsonic Mach, 0 <= M < 1");
	}
	const Real B = gm::sqrt(1 - M * M);
	return Result::Of(Cp0 / (B + M * M / (1 + B) * Cp0 / 2));
}

Result Ackeret(const Call& C)
{
	const Real Alpha = C.N(0), M = C.N(1);
	if (M <= 1)
	{
		return Result::Fail("Ackeret needs supersonic Mach, M > 1");
	}
	return Result::Of(4 * Alpha / gm::sqrt(M * M - 1));
}

Result Polar(const Call& C)
{
	const Real CD0 = C.N(0), AR = C.N(1), E = C.N(2), CL = C.N(3);
	if (CD0 <= 0 || AR <= 0 || E <= 0)
	{
		return Result::Fail("CD0, AR and e must be positive");
	}
	const Real K = 1 / (cst::Pi() * E * AR);
	return Result::Record({{"k", K}, {"LDmax", Real(0.5) / gm::sqrt(CD0 * K)}, {"CLstar", gm::sqrt(CD0 / K)}, {"LD", CL / (CD0 + K * CL * CL)}});
}

// Prandtl's lifting line by Glauert's Fourier method, untwisted wing, chord(y) with y = 2y/b in [0, 1].
template <typename F>
bool LiftingLine(F Chord, const Real& AR, int NTerms, Real* E, Real* Delta)
{
	const Real A0 = 2 * cst::Pi();
	std::vector<Real> Nk(static_cast<size_t>(NTerms)), Th(static_cast<size_t>(NTerms));
	for (int i = 0; i < NTerms; ++i)
	{
		Nk[static_cast<size_t>(i)] = Real(2 * i + 1);
		Th[static_cast<size_t>(i)] = (Real(i) + Real(0.5)) * (cst::Pi() / 2) / NTerms;
	}
	const Real Cbar = Trapz(Chord, Real(0), Real(1), 4000), B = AR * Cbar;
	std::vector<Real> M(static_cast<size_t>(NTerms) * static_cast<size_t>(NTerms + 1));
	for (int i = 0; i < NTerms; ++i)
	{
		const Real T = Th[static_cast<size_t>(i)], Cc = Chord(gm::cos(T)), SinT = gm::sin(T), Base = 4 * B / (A0 * Cc);
		for (int k = 0; k < NTerms; ++k)
		{
			const Real Nn = Nk[static_cast<size_t>(k)];
			M[static_cast<size_t>(i) * static_cast<size_t>(NTerms + 1) + static_cast<size_t>(k)] = gm::sin(Nn * T) * (Base + Nn / SinT);
		}
		M[static_cast<size_t>(i) * static_cast<size_t>(NTerms + 1) + static_cast<size_t>(NTerms)] = Real(1);
	}
	std::vector<Real> A;
	if (!GaussSolve(M, NTerms, &A))
	{
		return false;
	}
	Real D(0);
	for (int i = 1; i < NTerms; ++i)
	{
		const Real Q = A[static_cast<size_t>(i)] / A[0];
		D += Nk[static_cast<size_t>(i)] * (Q * Q);
	}
	*E = 1 / (1 + D);
	*Delta = D;
	return true;
}

Result SpanEff(const Call& C)
{
	const Real AR = C.N(0), Lam = C.N(1);
	const int NT = C.Int(2);
	if (AR <= 0 || Lam <= 0 || Lam > 10)
	{
		return Result::Fail("needs AR > 0 and taper ratio 0 < lambda <= 10");
	}
	if (NT < 2 || NT > 100)
	{
		return Result::Fail("Fourier terms must be 2 to 100");
	}
	Real E, D;
	const Real Slope = 1 - Lam;
	if (!LiftingLine([&Slope](const Real& Y) { return 1 - Slope * Y; }, AR, NT, &E, &D))
	{
		return Result::Fail("singular lifting-line system");
	}
	return Result::Record({{"e", E}, {"delta", D}});
}

Result SpanEffEll(const Call& C)
{
	const Real AR = C.N(0);
	const int NT = C.Int(1);
	if (AR <= 0 || NT < 2 || NT > 100)
	{
		return Result::Fail("needs AR > 0 and 2 to 100 Fourier terms");
	}
	Real E, D;
	if (!LiftingLine([](const Real& Y) { return gm::sqrt(gm::max(Real(0), 1 - Y * Y)); }, AR, NT, &E, &D))
	{
		return Result::Fail("singular lifting-line system");
	}
	return Result::Record({{"e", E}, {"delta", D}});
}

// Slender-body wave drag D/q by von Karman's Fourier series of S'(x); sin(n t) by the Chebyshev recurrence.
template <typename F>
Real WaveDq(F S, const Real& L)
{
	const int NTerms = 40, Nq = 2001;
	const Real H(1e-6);
	std::vector<Real> Sum(static_cast<size_t>(NTerms + 1), Real(0)), Prev(static_cast<size_t>(NTerms + 1)), Cur(static_cast<size_t>(NTerms + 1));
	Real TPrev(0), SpPrev(0);
	for (int i = 0; i < Nq; ++i)
	{
		const Real T = cst::Pi() * i / (Nq - 1), Ct = gm::cos(T), St = gm::sin(T);
		const Real X = L * (1 - Ct) / 2, A = gm::max(Real(0), X - H), B = gm::min(L, X + H);
		const Real Sp = (S(B) - S(A)) / (B - A);
		Cur[0] = Real(0);
		Cur[1] = St;
		const Real TwoC = 2 * Ct;
		for (int n = 2; n <= NTerms; ++n)
		{
			Cur[static_cast<size_t>(n)] = TwoC * Cur[static_cast<size_t>(n - 1)] - Cur[static_cast<size_t>(n - 2)];
		}
		if (i > 0)
		{
			const Real Dt = T - TPrev;
			for (int n = 1; n <= NTerms; ++n)
			{
				const size_t K = static_cast<size_t>(n);
				Sum[K] += (Sp * Cur[K] + SpPrev * Prev[K]) / 2 * Dt;
			}
		}
		std::swap(Prev, Cur);
		TPrev = T;
		SpPrev = Sp;
	}
	Real D(0);
	for (int n = 1; n <= NTerms; ++n)
	{
		const Real An = 2 / cst::Pi() * Sum[static_cast<size_t>(n)] / L;
		D += n * An * An;
	}
	return cst::Pi() * L * L / 4 * D;
}

Real SearsHaack(const Real& RMax, const Real& L, const Real& X)
{
	const Real U = gm::min(gm::max(X, Real(0)), L), Q = 4 * U * (L - U) / (L * L);
	return cst::Pi() * RMax * RMax * (Q * gm::sqrt(Q));
}

Real SearsHaackDq(const Real& RMax, const Real& L)
{
	const Real V = 3 * (cst::Pi() * cst::Pi()) * RMax * RMax * L / 16;
	return 128 * V * V / (cst::Pi() * gm::powi(L, 4));
}

Result WaveDrag(const Call& C)
{
	const Real R = C.N(0), L = C.N(1);
	if (R <= 0 || L <= 0)
	{
		return Result::Fail("rmax and length must be positive");
	}
	const Real F = WaveDq([&R, &L](const Real& X) { return SearsHaack(R, L, X); }, L), Cf = SearsHaackDq(R, L);
	return Result::Record({{"fourier", F}, {"closed", Cf}, {"ratio", F / Cf}});
}

// The equation sheet's area-rule card: Sears-Haack fuselage (rmax 0.05, l = 1) plus a wing's area, optionally indented.
Result AreaRule(const Call& C)
{
	const Real Wing = C.N(0), Rule = C.N(1), R(0.05), A0 = cst::Pi() * R * R, One(1);
	if (Wing < 0)
	{
		return Result::Fail("wing area share must be >= 0");
	}
	auto SH = [&R, &One](const Real& X) { return SearsHaack(R, One, X); };
	auto WingA = [&Wing, &A0](const Real& X) {
		const Real Q = (X - Real(0.55)) / Real(0.12), Y = gm::max(Real(0), 1 - Q * Q);
		return Wing * A0 * (Y * Y);
	};
	const Real DSH = WaveDq(SH, One), DW = WaveDq([&](const Real& X) { return SH(X) + WingA(X); }, One);
	const Real DT = WaveDq([&](const Real& X) { return SH(X) + WingA(X) - Rule * WingA(X); }, One);
	return Result::Record({{"fus", DSH}, {"wing", DW / DSH}, {"ruled", DT / DSH}, {"closed", SearsHaackDq(R, One)}});
}

Real AreaRatioOf(const Real& M, const Real& G)
{
	return 1 / M * gm::pow((2 / (G + 1)) * (1 + (G - 1) / 2 * M * M), (G + 1) / (2 * (G - 1)));
}

Result AreaRatio(const Call& C)
{
	const Real M = C.N(0), G = C.N(1);
	if (M <= 0 || G <= 1)
	{
		return Result::Fail("needs M > 0 and gamma > 1");
	}
	return Result::Of(AreaRatioOf(M, G));
}

Result MachFromArea(const Call& C)
{
	const Real R = C.N(0), G = C.N(2);
	const bool Sup = !C.N(1).IsZero();
	if (R < 1 || G <= 1)
	{
		return Result::Fail("area ratio A/A* must be >= 1 (and gamma > 1)");
	}
	Real Lo = Sup ? Real(1) : Real(1e-6), Hi = Sup ? Real(50) : Real(1);
	for (int i = 0; i < 200; ++i)
	{
		const Real Mid = (Lo + Hi) / 2;
		if ((AreaRatioOf(Mid, G) > R) == Sup)
		{
			Hi = Mid;
		}
		else
		{
			Lo = Mid;
		}
	}
	return Result::Of((Lo + Hi) / 2);
}

/* ---------------------------------------------------------------- equation sheet: propulsion */
Result ClimbAngle(const Call& C)
{
	const Real TW = C.N(0), LD = C.N(1);
	if (LD <= 0)
	{
		return Result::Fail("L/D must be positive");
	}
	if (TW <= 1 / LD)
	{
		return Result::Fail("T/W is at or below 1/(L/D): the aircraft cannot hold level flight");
	}
	Real Lo(0), Hi = cst::Pi() / 2;
	for (int i = 0; i < 200; ++i)
	{
		const Real Mid = (Lo + Hi) / 2;
		if (gm::sin(Mid) + gm::cos(Mid) / LD < TW)
		{
			Lo = Mid;
		}
		else
		{
			Hi = Mid;
		}
	}
	return Result::Of((Lo + Hi) / 2);
}

Result RocketDv(const Call& C)
{
	const Real Isp = C.N(0), M0 = C.N(1), Mf = C.N(2);
	if (M0 <= 0 || Mf <= 0)
	{
		return Result::Fail("masses must be positive");
	}
	return Result::Of(Isp * FlG0() * gm::log(M0 / Mf));
}

Result Turboprop(const Call& C)
{
	const Real Mdot = C.N(0), Pic = C.N(1), T02 = C.N(2), T04 = C.N(3), Etac = C.N(4), Etat = C.N(5), G = C.N(6), Cp = C.N(7);
	if (Pic <= 1 || Etac <= 0 || G <= 1)
	{
		return Result::Fail("needs pressure ratio > 1, compressor efficiency > 0 and gamma > 1");
	}
	const Real Tau = gm::pow(Pic, (G - 1) / G), Wt = Etat * T04 * (1 - 1 / Tau), Wc = T02 * (Tau - 1) / Etac, T03 = T02 + Wc, Q = Mdot * Cp * (T04 - T03);
	return Result::Record({{"tau", Tau}, {"P", Mdot * Cp * (Wt - Wc)}, {"turbine", Mdot * Cp * Wt}, {"compressor", Mdot * Cp * Wc}, {"q", Q}, {"eta", Mdot * Cp * (Wt - Wc) / Q}, {"ideal", 1 - 1 / Tau}});
}

Result NetThrust(const Call& C)
{
	const Real Mdot = C.N(0), F = C.N(1), Ve = C.N(2), V0 = C.N(3), Pe = C.N(4), P0 = C.N(5), Ae = C.N(6);
	const Real Gross = Mdot * (1 + F) * Ve + (Pe - P0) * Ae;
	return Result::Record({{"gross", Gross}, {"ram", Mdot * V0}, {"net", Mdot * (1 + F) * Ve + (Pe - P0) * Ae - Mdot * V0}});
}

/* ---------------------------------------------------------------- equation sheet: structures */
// Cantilever half-wing under elliptic running load w0 sqrt(1 - (y/s)^2), w0 from nW/2 per side; closed forms at y.
Result WingLoads(const Call& C)
{
	const Real NW = C.N(0), S = C.N(1), Y = C.N(2);
	if (S <= 0 || Y < 0 || Y > S)
	{
		return Result::Fail("needs semi-span s > 0 and 0 <= y <= s");
	}
	const Real W0 = 4 * (NW / 2) / (cst::Pi() * S), U = Y / S, Root = gm::sqrt(gm::max(Real(0), 1 - U * U)), Ac = gm::acos(U);
	const Real V = W0 * S * (Ac - U * Root) / 2;
	const Real M = W0 * S * S * (Root * Root * Root / 3 - U * (Ac - U * Root) / 2);
	return Result::Record({{"w0", W0}, {"w", W0 * Root}, {"V", V}, {"M", M}});
}

Result Torsion(const Call& C)
{
	const Real T = C.N(0), Rm = C.N(1), Tw = C.N(2), Gs = C.N(3);
	if (Rm <= 0 || Tw <= 0 || Tw >= 2 * Rm)
	{
		return Result::Fail("needs radius r > 0 and wall 0 < t < 2r");
	}
	const Real Am = cst::Pi() * Rm * Rm, Ro = Rm + Tw / 2, Ri = Rm - Tw / 2;
	return Result::Record({{"q", T / (2 * Am)},
		{"bredt", T / (2 * cst::Pi() * Rm * Rm * Tw)},
		{"exact", T * Ro / (cst::Pi() / 2 * (gm::powi(Ro, 4) - gm::powi(Ri, 4)))},
		{"open", 3 * T / (2 * cst::Pi() * Rm * Tw * Tw)},
		{"twist", T * (2 * cst::Pi() * Rm / Tw) / (4 * Am * Am * Gs)}});
}

// Smallest eigenvalue of the finite-difference column operator by inverse iteration (Thomas algorithm), times EI.
Real FdBuckle(const Real& E, const Real& I, const Real& L, bool Cantilever)
{
	const int N = 400, M = Cantilever ? N : N - 1;
	const Real H = L / N, H2 = H * H, Diag = 2 / H2, Off = -1 / H2;
	std::vector<Real> Up(static_cast<size_t>(M - 1), Off);
	if (Cantilever)
	{
		Up[0] = -2 / H2;
	}
	// The forward sweep depends only on the operator, so its multipliers are formed once.
	std::vector<Real> Dd(static_cast<size_t>(M), Diag), Fm(static_cast<size_t>(M), Real(0));
	for (int i = 1; i < M; ++i)
	{
		const size_t K = static_cast<size_t>(i);
		Fm[K] = Off / Dd[K - 1];
		Dd[K] -= Fm[K] * Up[K - 1];
	}
	std::vector<Real> X(static_cast<size_t>(M), Real(1)), Y(static_cast<size_t>(M)), R(static_cast<size_t>(M));
	Real Lam(0);
	for (int It = 0; It < 200; ++It)
	{
		R = X;
		for (int i = 1; i < M; ++i)
		{
			const size_t K = static_cast<size_t>(i);
			R[K] -= Fm[K] * R[K - 1];
		}
		Y[static_cast<size_t>(M - 1)] = R[static_cast<size_t>(M - 1)] / Dd[static_cast<size_t>(M - 1)];
		for (int i = M - 2; i >= 0; --i)
		{
			const size_t K = static_cast<size_t>(i);
			Y[K] = (R[K] - Up[K] * Y[K + 1]) / Dd[K];
		}
		Real Num(0), Den(0);
		for (size_t i = 0; i < static_cast<size_t>(M); ++i)
		{
			Num += X[i] * Y[i];
			Den += Y[i] * Y[i];
		}
		Lam = Num / Den;
		const Real Nrm = gm::sqrt(Den);
		for (size_t i = 0; i < static_cast<size_t>(M); ++i)
		{
			X[i] = Y[i] / Nrm;
		}
	}
	return E * I * Lam;
}

Result Buckle(const Call& C)
{
	const Real E = C.N(0), I = C.N(1), L = C.N(2), K = C.N(3);
	if (E <= 0 || I <= 0 || L <= 0 || K <= 0)
	{
		return Result::Fail("E, I, L and K must be positive");
	}
	const Real KL = K * L, P = cst::Pi() * cst::Pi() * E * I / (KL * KL);
	const Real Fd = K == 1 ? FdBuckle(E, I, L, false) : (K == 2 ? FdBuckle(E, I, L, true) : FlNan());
	return Result::Record({{"P", P}, {"fd", Fd}});
}

/* ---------------------------------------------------------------- equation sheet: flight mechanics */
Real KSheet() { return 1 / (cst::Pi() * Real(0.8) * 8); } // the sheet's aircraft: AR 8, e 0.8

Result BestClimb(const Call& C)
{
	const Real W = C.N(0), S = C.N(1), Rho = C.N(2), Pav = C.N(3), CD0 = C.N(4), K = C.N(5);
	if (W <= 0 || S <= 0 || Rho <= 0 || K <= 0 || CD0 <= 0)
	{
		return Result::Fail("W, S, rho, CD0 and k must be positive");
	}
	int B = 0;
	Real Best(0);
	const Real Step(0.1), HalfRhoS = Real(0.5) * Rho * S;
	for (int i = 0; i < 601; ++i)
	{
		const Real V = 15 + Real(i) * Step, CL = 2 * W / (Rho * V * V * S);
		const Real Roc = (Pav - HalfRhoS * (V * V * V) * (CD0 + K * CL * CL)) / W;
		if (i == 0 || Roc > Best)
		{
			Best = Roc;
			B = i;
		}
	}
	return Result::Record({{"Vbest", 15 + Real(B) * Step}, {"rocmax", Best}, {"Vmp", gm::sqrt(2 * W / (Rho * S * gm::sqrt(3 * CD0 / K)))}});
}

Result FnTurn(const Call& C)
{
	const Real V = C.N(0), Phi = C.N(1);
	if (V <= 0 || gm::fabs(Phi) >= cst::Pi() / 2)
	{
		return Result::Fail("needs V > 0 and |bank| < pi/2");
	}
	const Real N = 1 / gm::cos(Phi), S = gm::sqrt(N * N - 1);
	return Result::Record({{"n", N}, {"omega", FlG0() * S / V}, {"R", S > 0 ? V * V / (FlG0() * S) : FlInf()}});
}

Result CornerSpeed(const Call& C)
{
	const Real Nmax = C.N(0), W = C.N(1), Rho = C.N(2), S = C.N(3), CLmax = C.N(4);
	if (Nmax <= 0 || W <= 0 || Rho <= 0 || S <= 0 || CLmax <= 0)
	{
		return Result::Fail("all inputs must be positive");
	}
	return Result::Of(gm::sqrt(2 * Nmax * W / (Rho * S * CLmax)));
}

Result FnGlide(const Call& C)
{
	const Real H = C.N(0), W = C.N(1), S = C.N(2), Rho = C.N(3), CD0 = C.N(4), K = C.N(5), CL = C.N(6);
	if (CL <= 0 || W <= 0 || S <= 0 || Rho <= 0)
	{
		return Result::Fail("CL, W, S and rho must be positive");
	}
	const Real Gam = gm::atan((CD0 + K * CL * CL) / CL), V = gm::sqrt(2 * W * gm::cos(Gam) / (Rho * S * CL)), Tg = gm::tan(Gam);
	return Result::Record({{"gam", Gam}, {"V", V}, {"ratio", 1 / Tg}, {"dist", H / Tg}, {"sink", V * gm::sin(Gam)}});
}

/* ---------------------------------------------------------------- equation sheet: guidance and control */
struct SpModel
{
	Real Za, Zd, Ma, Mq, Md;
};

SpModel SpFrom(const Call& C, size_t First)
{
	return {C.N(First), C.N(First + 1), C.N(First + 2), C.N(First + 3), C.N(First + 4)};
}

const SpModel& SpDefault()
{
	static const SpModel P{Real(-2.0), Real(-0.2), Real(-12.0), Real(-3.0), Real(-15.0)};
	return P;
}

Result ShortPeriod(const Call& C)
{
	const SpModel P = SpFrom(C, 0);
	const Real A1 = -(P.Za + P.Mq), A0 = P.Za * P.Mq - P.Ma;
	if (A0 <= 0)
	{
		return Result::Fail("Za Mq - Ma <= 0: no short-period oscillator (statically unstable)");
	}
	const Real Wn = gm::sqrt(A0), Disc = A1 * A1 - 4 * A0;
	const Real Td = Disc < 0 ? 2 * cst::Pi() / (gm::sqrt(-Disc) / 2) : FlNan();
	return Result::Record({{"wn", Wn}, {"zeta", A1 / (2 * Wn)}, {"a1", A1}, {"a0", A0}, {"b1", P.Md}, {"b0", P.Ma * P.Zd - P.Md * P.Za}, {"Td", Td}});
}

// theta/delta at s = j w from the transfer function, and from (jwI - A)^-1 B as a real 6 x 6 system.
Result PitchTf(const Call& C)
{
	const Real W = C.N(0);
	const SpModel P = SpFrom(C, 1);
	const Real A1 = -(P.Za + P.Mq), A0 = P.Za * P.Mq - P.Ma, N1 = P.Md, N0 = P.Ma * P.Zd - P.Md * P.Za;
	const Real S0(0), S1 = W;
	const Real NumR = N1 * S0 + N0, NumI = N1 * S1;
	const Real S2R = S0 * S0 - S1 * S1, S2I = S0 * S1 + S1 * S0;
	const Real QR = S2R + A1 * S0 + A0, QI = S2I + A1 * S1;
	const Real DenR = S0 * QR - S1 * QI, DenI = S0 * QI + S1 * QR;
	const Real Dd = DenR * DenR + DenI * DenI;
	if (Dd.IsZero())
	{
		return Result::Fail("w = 0 is the free-attitude pole");
	}
	const Real Re = (NumR * DenR + NumI * DenI) / Dd, Im = (NumI * DenR - NumR * DenI) / Dd;
	const Real A[3][3] = {{P.Za, Real(1), Real(0)}, {P.Ma, P.Mq, Real(0)}, {Real(0), Real(1), Real(0)}};
	const Real B[3] = {P.Zd, P.Md, Real(0)};
	std::vector<Real> M(6 * 7, Real(0));
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			M[static_cast<size_t>(i * 7 + j)] = -A[i][j];
			M[static_cast<size_t>(i * 7 + 3 + j)] = i == j ? -W : Real(0);
			M[static_cast<size_t>((i + 3) * 7 + j)] = i == j ? W : Real(0);
			M[static_cast<size_t>((i + 3) * 7 + 3 + j)] = -A[i][j];
		}
		M[static_cast<size_t>(i * 7 + 6)] = B[i];
	}
	std::vector<Real> X;
	if (!GaussSolve(M, 6, &X))
	{
		return Result::Fail("singular state-space system");
	}
	const Real Mag = gm::hypot(Re, Im);
	return Result::Record({{"re", Re}, {"im", Im}, {"mag", Mag}, {"db", 20 * gm::log10(Mag)}, {"ph", gm::atan2(Im, Re)}, {"ssre", X[2]}, {"ssim", X[5]}});
}

Result POnlyError(const Call& C)
{
	const Real Kp = C.N(0), Dist = C.N(1);
	const SpModel P = SpFrom(C, 2);
	if (Kp.IsZero() || P.Za.IsZero())
	{
		return Result::Fail("needs Kp != 0 and Za != 0");
	}
	return Result::Of(Dist / (Kp * (P.Md - P.Ma * P.Zd / P.Za)));
}

Result PidSim(const Call& C)
{
	const Real Kp = C.N(0), Ki = C.N(1), Kd = C.N(2), ThetaC = C.N(3), Dist = C.N(4), Dt(1e-3);
	const SpModel& P = SpDefault();
	Real A(0), Q(0), Th(0), Z(0), Peak(0);
	const Real Lim(1000);
	for (int i = 0; i <= 20000; ++i)   // T / dt = 20 s / 1 ms
	{
		const Real E = ThetaC - Th;
		Peak = gm::max(Peak, Th);
		const Real D = -(Kp * E + Ki * Z) + Kd * Q, Da = P.Za * A + Q + P.Zd * D, Dq = P.Ma * A + P.Mq * Q + P.Md * D + Dist;
		A += Dt * Da;
		Q += Dt * Dq;
		Th += Dt * Q;
		Z += Dt * E;
		if (!gm::isfinite(Th) || gm::fabs(Th) > Lim)
		{
			return Result::Record({{"err", FlNan()}, {"overshoot", FlNan()}, {"diverged", Real(1)}});
		}
	}
	return Result::Record({{"err", ThetaC - Th}, {"overshoot", !ThetaC.IsZero() ? gm::max(Real(0), (Peak - ThetaC) / ThetaC) : Real(0)}, {"diverged", Real(0)}});
}

/* ---------------------------------------------------------------- model.js: atmosphere and aircraft */
struct Atm
{
	Real T, Tc, P, Rho, A, Sigma;
	Real Sigma07;            // sigma^0.7 for the jet lapse, when HasLapse
	bool HasLapse = false;
};

Atm WithLapse(Atm A)
{
	A.Sigma07 = gm::pow(A.Sigma, Real(0.7));
	A.HasLapse = true;
	return A;
}

// NASA Glenn Earth Atmosphere Model (metric), every temperature shifted by dT.
Atm Atmosphere(const Real& HIn, const Real& DT)
{
	const Real H = gm::max(Real(0), HIn);
	Real Tc, Pk;
	if (H < 11000)
	{
		Tc = Real(15.04) - Real(0.00649) * H;
		Pk = Real(101.29) * gm::pow((Tc + Real(273.1)) / Real(288.08), Real(5.256));
	}
	else if (H < 25000)
	{
		Tc = Real(-56.46);
		Pk = Real(22.65) * gm::exp(Real(1.73) - Real(0.000157) * H);
	}
	else
	{
		Tc = Real(-131.21) + Real(0.00299) * H;
		Pk = Real(2.488) * gm::pow((Tc + Real(273.1)) / Real(216.6), Real(-11.388));
	}
	Tc += DT;
	const Real T = Tc + Real(273.1), Rho = Pk / (Real(0.2869) * T);
	Atm Out;
	Out.T = T;
	Out.Tc = Tc;
	Out.P = Pk * 1000;
	Out.Rho = Rho;
	Out.A = gm::sqrt(Real(1.4) * 286 * T);
	Out.Sigma = Rho / Real(1.225);
	return Out;
}

struct Aircraft
{
	Real M, S, AR, E, CD0, Alpha0, AlphaStall, PlateCL, Mcrit, CDwave, CYb;
	bool Prop;
	Real Power, Tstatic, EtaProp, Tmax, KM;
	Real GroundAlpha, RotateMax, Vrot;
};

// The three representative sizings of model.js: 1 general-aviation single, 2 jet transport, 3 supersonic delta.
// Angles are entered in degrees and converted on the gears, as model.js's rad().
const Aircraft* AircraftOf(int Index)
{
	static const Aircraft Table[3] = {
		{Real(1050), Real(16.2), Real(7.4), Real(0.75), Real(0.028), FlRad(Real(-2.5)), FlRad(Real(15)), Real(1.0), Real(0.7), Real(0.05), Real(-0.6), true,
			Real(120e3), Real(3200), Real(0.8), Real(0), Real(0), FlRad(Real(1)), FlRad(Real(12)), Real(28)},
		{Real(300000), Real(511), Real(7.0), Real(0.8), Real(0.019), FlRad(Real(-3)), FlRad(Real(14)), Real(1.0), Real(0.82), Real(0.06), Real(-0.8), false,
			Real(0), Real(0), Real(0), Real(1000000), Real(-0.35), FlRad(Real(0)), FlRad(Real(11)), Real(82)},   // Tmax 4 x 250 kN
		{Real(150000), Real(358), Real(1.6), Real(0.85), Real(0.012), Real(0), FlRad(Real(22)), Real(1.0), Real(0.92), Real(0.022), Real(-0.5), false,
			Real(0), Real(0), Real(0), Real(680000), Real(0.7), FlRad(Real(0)), FlRad(Real(14)), Real(100)},     // Tmax 4 x 170 kN
	};
	return Index >= 1 && Index <= 3 ? &Table[Index - 1] : nullptr;
}

Real SectionSlope(const Real& M)
{
	auto Sub = [](const Real& X) { return (2 * cst::Pi()) / gm::sqrt(1 - X * X); };
	auto Sup = [](const Real& X) { return 4 / gm::sqrt(X * X - 1); };
	if (M <= Real(0.85))
	{
		return Sub(M);
	}
	if (M >= Real(1.15))
	{
		return Sup(M);
	}
	const Real T = (M - Real(0.85)) / Real(0.3);
	return (1 - T) * Sub(Real(0.85)) + T * Sup(Real(1.15));
}

Real WingSlope(const Real& A0, const Real& AR)
{
	const Real K = A0 / (cst::Pi() * AR);
	return A0 / (gm::sqrt(1 + K * K) + K);
}

Real StallBlend(const Real& Alpha, const Aircraft& Ac)
{
	static const Real Width = FlRad(Real(2.0));
	return 1 / (1 + gm::exp(-(gm::fabs(Alpha) - Ac.AlphaStall) / Width));
}

struct Aero
{
	Atm Air;
	Real M, Q, CL, CLa, CLlin, CD0, CDi, CD, K, L, D, SinA, CosA;
	bool Stalled = false;
};

Aero AeroAt(const Real& Alpha, const Real& Beta, const Real& V, const Atm& Air, const Aircraft& Ac)
{
	Aero R;
	R.Air = Air;
	R.M = V / R.Air.A;
	R.Q = Real(0.5) * R.Air.Rho * V * V;
	R.CLa = WingSlope(SectionSlope(R.M), Ac.AR);
	R.CLlin = R.CLa * (Alpha - Ac.Alpha0);
	const Real Sa = gm::sin(Alpha), Ca = gm::cos(Alpha), Plate = 2 * Sa * Ca * Ac.PlateCL, S = StallBlend(Alpha, Ac);
	R.SinA = Sa;
	R.CosA = Ca;
	R.CL = (1 - S) * R.CLlin + S * Plate;
	R.Stalled = S > Real(0.5);
	Real Wave(0);
	const Real Top(1.05);
	if (R.M > Ac.Mcrit && R.M <= Top)
	{
		const Real U = (R.M - Ac.Mcrit) / (Top - Ac.Mcrit);
		Wave = Ac.CDwave * U * U * U;
	}
	else if (R.M > Top)
	{
		Wave = Ac.CDwave * gm::sqrt((Top * Top - 1) / (R.M * R.M - 1));
	}
	R.CD0 = Ac.CD0 + Wave;
	R.K = 1 / (cst::Pi() * Ac.AR * Ac.E);
	if (R.M > 1)
	{
		R.K = gm::max(R.K, gm::sqrt(R.M * R.M - 1) / 4);
	}
	R.CDi = (1 - S) * R.K * R.CL * R.CL + S * 2 * (Sa * Sa) * Ac.PlateCL;
	R.CD = R.CD0 + R.CDi + Real(0.5) * Beta * Beta;
	R.L = R.Q * Ac.S * R.CL;
	R.D = R.Q * Ac.S * R.CD;
	return R;
}

Real ThrustOf(const Real& Throttle, const Real& V, const Atm& Air, const Real& M, const Aircraft& Ac)
{
	if (Ac.Prop)
	{
		const Real P = Ac.Power * Throttle * Air.Sigma, Ts = Ac.Tstatic * Throttle * Air.Sigma;
		return gm::min(Ts, (Ac.EtaProp * P) / gm::max(V, Real(1)));
	}
	return Ac.Tmax * Throttle * (Air.HasLapse ? Air.Sigma07 : gm::pow(Air.Sigma, Real(0.7))) * gm::max(Real(0.2), 1 + Ac.KM * M);
}

Real TrimAlpha(const Aircraft& Ac, const Real& V, const Real& H, const Real& Gam, const Real& Phi, const Real& T)
{
	const Atm Air = Atmosphere(H, Real(0));
	const Real W = Ac.M * FlG0(), Need = (W * gm::cos(Gam)) / gm::cos(Phi);
	Real Lo = Ac.Alpha0 - FlRad(Real(4)), Hi = Ac.AlphaStall;
	for (int i = 0; i < 60; ++i)
	{
		const Real Mid = Real(0.5) * (Lo + Hi);
		if (AeroAt(Mid, Real(0), V, Air, Ac).L + T * gm::sin(Mid) < Need)
		{
			Lo = Mid;
		}
		else
		{
			Hi = Mid;
		}
	}
	return Real(0.5) * (Lo + Hi);
}

const Aircraft* ArgAircraft(const Call& C, size_t I)
{
	const Real X = C.N(I);
	return X == gm::floor(X) ? AircraftOf(C.Int(I)) : nullptr;
}

const char* const BadAircraft = "aircraft: 1 general aviation, 2 jet transport, 3 supersonic delta";

Result FnAtmos(const Call& C)
{
	const Atm A = Atmosphere(C.N(0), C.N(1));
	return Result::Record({{"T", A.T}, {"Tc", A.Tc}, {"p", A.P}, {"rho", A.Rho}, {"a", A.A}, {"sigma", A.Sigma}, {"icing", Flag(A.Tc <= 2 && A.Tc >= -20)}});
}

Result LiftSlope(const Call& C)
{
	const Real M = C.N(0), AR = C.N(1);
	if (M < 0 || AR < 0)
	{
		return Result::Fail("needs M >= 0 and AR >= 0 (0 = the section)");
	}
	const Real A0 = SectionSlope(M);
	return Result::Record({{"a0", A0}, {"a", AR > 0 ? WingSlope(A0, AR) : A0}});
}

Result LiftCoef(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 2);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const Real Alpha = C.N(0), M = C.N(1);
	if (M < 0)
	{
		return Result::Fail("Mach must be >= 0");
	}
	const Real CLa = WingSlope(SectionSlope(M), Ac->AR), Lin = CLa * (Alpha - Ac->Alpha0);
	const Real Plate = 2 * gm::sin(Alpha) * gm::cos(Alpha) * Ac->PlateCL, S = StallBlend(Alpha, *Ac);
	return Result::Record({{"CL", (1 - S) * Lin + S * Plate}, {"CLa", CLa}, {"linear", Lin}, {"stalled", Flag(S > Real(0.5))}});
}

Result FnAero(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 3);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	if (C.N(1) <= 0)
	{
		return Result::Fail("airspeed must be positive");
	}
	const Aero A = AeroAt(C.N(0), C.N(4), C.N(1), Atmosphere(C.N(2), Real(0)), *Ac);
	return Result::Record({{"L", A.L}, {"D", A.D}, {"CL", A.CL}, {"CD", A.CD}, {"CD0", A.CD0}, {"CDi", A.CDi}, {"K", A.K}, {"M", A.M}, {"q", A.Q}});
}

Result FnThrust(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 3);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const Real Throttle = C.N(0), V = C.N(1);
	const Atm A = Atmosphere(C.N(2), Real(0));
	return Result::Of(ThrustOf(Throttle, V, A, V / A.A, *Ac));
}

Result FnTrim(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 2);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const Real V = C.N(0), H = C.N(1), Gam = C.N(3), Phi = C.N(4);
	if (V <= 0 || gm::fabs(Phi) >= cst::Pi() / 2)
	{
		return Result::Fail("needs V > 0 and |bank| < pi/2");
	}
	const Atm Air = Atmosphere(H, Real(0));
	const Real W = Ac->M * FlG0(), SinGam = gm::sin(Gam);
	Real T(0), Alpha(0);
	for (int i = 0; i < 30; ++i)
	{
		Alpha = TrimAlpha(*Ac, V, H, Gam, Phi, T);
		const Aero A = AeroAt(Alpha, Real(0), V, Air, *Ac);
		T = gm::max(Real(0), (A.D + W * SinGam) / gm::cos(Alpha));
	}
	const Real M = V / Air.A;
	if (ThrustOf(Real(1), V, Air, M, *Ac) < T)
	{
		return Result::Record({{"alpha", Alpha}, {"throttle", FlNan()}, {"T", T}});
	}
	Real Lo(0), Hi(1);
	for (int i = 0; i < 50; ++i)
	{
		const Real Mid = Real(0.5) * (Lo + Hi);
		if (ThrustOf(Mid, V, Air, M, *Ac) < T)
		{
			Lo = Mid;
		}
		else
		{
			Hi = Mid;
		}
	}
	return Result::Record({{"alpha", Alpha}, {"throttle", Real(0.5) * (Lo + Hi)}, {"T", T}});
}

Result LevelDrag(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 2);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const Real V = C.N(0), H = C.N(1);
	if (V <= 0)
	{
		return Result::Fail("airspeed must be positive");
	}
	const Real A = TrimAlpha(*Ac, V, H, Real(0), Real(0), Real(0));
	const Aero R = AeroAt(A, Real(0), V, Atmosphere(H, Real(0)), *Ac);
	const bool Ok = R.L > Real(0.995) * Ac->M * FlG0() && A < Ac->AlphaStall - FlRad(Real(0.5));
	return Result::Record({{"D", R.D}, {"alpha", A}, {"CL", R.CL}, {"CD", R.CD}, {"M", R.M}, {"gam", gm::atan(R.D / R.L)}, {"ok", Flag(Ok)}});
}

// Ground roll from rest at full throttle (model.js step() on the ground, RK4 at 1/120 s), elevator once V > Vrot.
struct Roll
{
	Real X, V, Alpha, Qa;
};

Roll RollRate(const Roll& S, const Real& Throttle, const Real& Elev, const Atm& Air, const Aircraft& Ac)
{
	const Real W = Ac.M * FlG0(), V = gm::max(S.V, Real(0.5));
	const Aero A = AeroAt(S.Alpha, Real(0), V, Air, Ac);
	const Real T = ThrustOf(Throttle, V, A.Air, A.M, Ac);
	const Real Lz = A.L + T * A.SinA, N = gm::max(Real(0), W - Lz), Mu(0.02);
	const Real Fr = V > Real(0.6) || T > Mu * N ? Mu * N : Real(0);
	const Real GA = Ac.GroundAlpha, RM = Ac.RotateMax, HalfVrot = Real(0.5) * Ac.Vrot;
	const Real Auth = Clamp((V - HalfVrot) / HalfVrot, Real(0), Real(1));
	const Real Target = Clamp(GA + gm::max(Real(0), Elev) * Auth * (RM - GA), GA, RM);
	return {V, (T * A.CosA - A.D - Fr) / Ac.M, S.Qa, 9 * (Target - S.Alpha) - 6 * S.Qa};
}

Roll RollAdd(const Roll& S, const Roll& D, const Real& K) { return {S.X + K * D.X, S.V + K * D.V, S.Alpha + K * D.Alpha, S.Qa + K * D.Qa}; }

Result Takeoff(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 0);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const Real Throttle = C.N(1), ElevCmd = C.N(2), Dt = Real(1) / 120, Half = Dt / 2, Sixth = Dt / 6;
	if (Throttle <= 0 || Throttle > 1)
	{
		return Result::Fail("throttle must be in (0, 1]");
	}
	const Atm Air = WithLapse(Atmosphere(Real(0), Real(0)));
	const Real Weight = Ac->M * FlG0() * Real(1.001);
	Roll S{Real(0), Real(0), Ac->GroundAlpha, Real(0)};
	Real T(0);
	for (int i = 0; i < 120 * 120; ++i)
	{
		const Real Elev = S.V > Ac->Vrot ? ElevCmd : Real(0);
		const Roll K1 = RollRate(S, Throttle, Elev, Air, *Ac);
		const Roll K2 = RollRate(RollAdd(S, K1, Half), Throttle, Elev, Air, *Ac);
		const Roll K3 = RollRate(RollAdd(S, K2, Half), Throttle, Elev, Air, *Ac);
		const Roll K4 = RollRate(RollAdd(S, K3, Dt), Throttle, Elev, Air, *Ac);
		Roll O;
		O.X = S.X + Sixth * (K1.X + 2 * K2.X + 2 * K3.X + K4.X);
		O.V = S.V + Sixth * (K1.V + 2 * K2.V + 2 * K3.V + K4.V);
		O.Alpha = S.Alpha + Sixth * (K1.Alpha + 2 * K2.Alpha + 2 * K3.Alpha + K4.Alpha);
		O.Qa = S.Qa + Sixth * (K1.Qa + 2 * K2.Qa + 2 * K3.Qa + K4.Qa);
		O.V = gm::max(O.V, Real(0));
		T += Dt;
		S = O;
		const Aero A = AeroAt(S.Alpha, Real(0), gm::max(S.V, Real(0.5)), Air, *Ac);
		const Real Th = ThrustOf(Throttle, S.V, A.Air, A.M, *Ac);
		if (A.L + Th * A.SinA > Weight)
		{
			return Result::Record({{"t", T}, {"dist", S.X}, {"V", S.V}, {"alpha", S.Alpha}});
		}
	}
	return Result::Fail("no liftoff within 120 s");
}

/* ---------------------------------------------------------------- model.js: the cloud chain */
struct CloudConst
{
	Real Cp, Rd, Lv, Eps, Es0, A, B;
};

const CloudConst& Cl()
{
	static const CloudConst K{Real(1004), Real(287), Real(2.501e6), Real(0.622), Real(611.2), Real(17.67), Real(243.5)};
	return K;
}

Real SatVapour(const Real& Tc) { return Cl().Es0 * gm::exp((Cl().A * Tc) / (Tc + Cl().B)); }
Real SatMixing(const Real& Tc, const Real& P)
{
	const Real E = SatVapour(Tc);
	return (Cl().Eps * E) / gm::max(P - E, Real(1));
}
Real MoistLapseOf(const Real& T, const Real& P)
{
	const CloudConst& K = Cl();
	const Real Rs = SatMixing(T - Real(273.15), P);
	return (FlG0() * (1 + (K.Lv * Rs) / (K.Rd * T))) / (K.Cp + (K.Lv * K.Lv * Rs * K.Eps) / (K.Rd * T * T));
}
Real DewPointOf(const Real& R, const Real& P)
{
	const Real E = (R * P) / (Cl().Eps + R), L = gm::log(E / Cl().Es0);
	return (Cl().B * L) / (Cl().A - L);
}

struct Sounding
{
	Real Lcl, Lfc, El, Cape, Cin;
	bool Convective = false;
	std::vector<Real> CapeCum;   // per 25 m level, as prof[].capeCum
};

// Lift a surface parcel through the environment: dry below cloud base, moist above (model.js sounding()).
Sounding RunSounding(const Real& DT, const Real& DewSpread, const Real& Heating)
{
	const int Dz = 25;
	const Atm Env0 = Atmosphere(Real(0), DT);
	Real Tp = Env0.T + Heating;
	const Real Td0 = Env0.Tc - DewSpread, R0 = SatMixing(Td0, Env0.P), Dry = FlG0() / Cl().Cp, K273(273.15);
	bool HasLcl = false, HasLfc = false, HasEl = false;
	Sounding Out{FlNan(), FlNan(), FlNan(), Real(0), Real(0), false, {}};
	Real CapeCum(0);
	for (int Z = 0; Z <= 16000; Z += Dz)
	{
		const Atm Env = Atmosphere(Real(Z), DT);
		const bool Saturated = HasLcl;
		const Real Tdp = Saturated ? Tp - K273 : DewPointOf(R0, Env.P);
		if (!Saturated && Tp - K273 <= Tdp)
		{
			HasLcl = true;
			Out.Lcl = Real(Z);
		}
		const Real B = (FlG0() * (Tp - Env.T)) / Env.T, BDz = B * Dz;
		if (HasLcl)
		{
			if (B > 0 && !HasLfc)
			{
				HasLfc = true;
				Out.Lfc = Real(Z);
			}
			if (HasLfc && !HasEl)
			{
				if (B > 0)
				{
					Out.Cape += BDz;
					CapeCum += BDz;
				}
				else
				{
					HasEl = true;
					Out.El = Real(Z);
				}
			}
			else if (!HasLfc)
			{
				Out.Cin += gm::min(Real(0), B) * Dz;
			}
		}
		else
		{
			Out.Cin += gm::min(Real(0), B) * Dz;
		}
		if (!HasLfc && Z >= 6000)
		{
			break;
		}
		Out.CapeCum.push_back(CapeCum);
		Tp -= (!HasLcl ? Dry : MoistLapseOf(Tp, Env.P)) * Dz;
	}
	if (HasLfc && !HasEl)
	{
		Out.El = Real(16000);
	}
	if (!HasLfc)
	{
		Out.Cin = FlNan();
	}
	Out.Convective = HasLfc && Out.Cape > 50;
	return Out;
}

Result MoistLapse(const Call& C)
{
	if (C.N(0) <= 0 || C.N(1) <= 0)
	{
		return Result::Fail("T (K) and p (Pa) must be positive");
	}
	return Result::Of(MoistLapseOf(C.N(0), C.N(1)));
}

Result SatVap(const Call& C)
{
	if (C.N(1) <= 0)
	{
		return Result::Fail("pressure must be positive");
	}
	return Result::Record({{"e", SatVapour(C.N(0))}, {"rs", SatMixing(C.N(0), C.N(1))}});
}

Result DewPoint(const Call& C)
{
	if (C.N(0) <= 0 || C.N(1) <= 0)
	{
		return Result::Fail("mixing ratio and pressure must be positive");
	}
	return Result::Of(DewPointOf(C.N(0), C.N(1)));
}

Result FnSounding(const Call& C)
{
	const Sounding S = RunSounding(C.N(0), C.N(1), C.N(2));
	return Result::Record({{"lcl", S.Lcl}, {"lfc", S.Lfc}, {"el", S.El}, {"cape", S.Cape}, {"cin", S.Cin}, {"conv", Flag(S.Convective)}});
}

Result Updraft(const Call& C)
{
	const Sounding S = RunSounding(C.N(1), C.N(2), C.N(3));
	if (gm::isnan(S.Lfc) || S.CapeCum.empty())
	{
		return Result::Of(Real(0));
	}
	const long long Last = static_cast<long long>(S.CapeCum.size()) - 1;
	const long long Idx = std::min(Last, std::max(0LL, JsRound(C.N(0) / 25).ToLong()));
	return Result::Of(gm::min(Real(25), Real(0.4) * gm::sqrt(2 * S.CapeCum[static_cast<size_t>(Idx)]) + Real(1.5)));
}

// A sharp vertical gust w tilts the relative wind by w/V (model.js derivatives()); load factor from level trim.
Result GustLoad(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 3);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const Real Wg = C.N(0), V = C.N(1), H = C.N(2);
	if (V <= 0)
	{
		return Result::Fail("airspeed must be positive");
	}
	const Atm Air = Atmosphere(H, Real(0));
	const Real W = Ac->M * FlG0(), A0 = TrimAlpha(*Ac, V, H, Real(0), Real(0), Real(0)), Da = Wg / V;
	const Real N0 = AeroAt(A0, Real(0), V, Air, *Ac).L / W, N = AeroAt(A0 + Da, Real(0), V, Air, *Ac).L / W;
	return Result::Record({{"dalpha", Da}, {"n", N}, {"dn", N - N0}, {"alpha0", A0}});
}

/* ---------------------------------------------------------------- planform.js: the airlifter */
struct Table
{
	int Rows = 0, Cols = 0;
	std::vector<Real> V;
	const Real& At(int R, int C) const { return V[static_cast<size_t>(R) * static_cast<size_t>(Cols) + static_cast<size_t>(C)]; }
};

const Table& DefaultTable()
{
	static const Table T = [] {
		const Real D2R = cst::Pi() / 180, T55 = gm::tan(55 * D2R), T35 = gm::tan(35 * D2R);
		const Real LE9 = -9 * T55, TE4 = -46 + 4 * T35, TE13 = TE4 + 9 * T55;
		Table Out;
		Out.Rows = 5;
		Out.Cols = 6;
		Out.V = {Real(0), Real(0), Real(-46), Real(6.0), Real(0.00), Real(0.62),
			Real(4), -4 * T55, TE4, Real(5.6), Real(0.08), Real(0.61),
			Real(9), LE9, TE4 + 5 * T55, Real(3.4), Real(0.35), Real(0.58),
			Real(13), LE9 - 4 * T35, TE13, Real(1.7), Real(0.70), Real(0.56),
			Real(31), LE9 - 22 * T35, Real(-31.6), Real(0.34), Real(1.45), Real(0.54)};
		return Out;
	}();
	return T;
}

// A table argument: 0 (the default) means the generator's five stations; otherwise n x k, k >= MinCols.
bool TableArg(const Call& C, size_t I, int MinCols, Table* Out, std::string* Error)
{
	const Matrix& M = C.M(I);
	if (M.IsScalar() && M.ScalarValue().IsZero())
	{
		*Out = DefaultTable();
		return true;
	}
	if (M.Rows < 2 || M.Cols < MinCols)
	{
		*Error = "table must be n x " + std::to_string(MinCols) + " or wider (y, x_LE, x_TE" + (MinCols > 3 ? ", t, z_c, upper share" : "") + "), n >= 2; 0 for the default";
		return false;
	}
	Out->Rows = M.Rows;
	Out->Cols = M.Cols;
	Out->V = M.V;
	for (int i = 0; i < M.Rows; ++i)
	{
		const std::string S = "station " + std::to_string(i + 1) + ": ";
		if (i == 0 && !Out->At(0, 0).IsZero())
		{
			*Error = "the first station must be the centreline, y = 0";
			return false;
		}
		if (i > 0 && !(Out->At(i, 0) > Out->At(i - 1, 0)))
		{
			*Error = S + "y must increase outboard";
			return false;
		}
		if (!(Out->At(i, 1) > Out->At(i, 2)))
		{
			*Error = S + "the leading edge must be ahead of the trailing edge (x_LE > x_TE)";
			return false;
		}
		if (Out->Cols >= 4 && !(Out->At(i, 3) > 0))
		{
			*Error = S + "thickness must be positive";
			return false;
		}
		if (Out->Cols >= 6 && !(Out->At(i, 5) > 0 && Out->At(i, 5) < 1))
		{
			*Error = S + "the upper share must lie between 0 and 1";
			return false;
		}
	}
	return true;
}

Real Interp(const Table& T, int Col, const Real& X)   // numpy.interp on column Col against y, clamped
{
	const int N = T.Rows;
	if (X <= T.At(0, 0))
	{
		return T.At(0, Col);
	}
	if (X >= T.At(N - 1, 0))
	{
		return T.At(N - 1, Col);
	}
	int i = 1;
	while (T.At(i, 0) < X)
	{
		++i;
	}
	const Real F = (X - T.At(i - 1, 0)) / (T.At(i, 0) - T.At(i - 1, 0));
	return T.At(i - 1, Col) + F * (T.At(i, Col) - T.At(i - 1, Col));
}

struct Stn
{
	Real Xle, Xte, T, Zc, Up;
};

Stn StationAt(const Table& T, const Real& Y)
{
	const Real A = gm::fabs(Y);
	return {Interp(T, 1, A), Interp(T, 2, A), Interp(T, 3, A), Interp(T, 4, A), Interp(T, 5, A)};
}

const Real* NacaCo()
{
	static const Real Co[5] = {Real(0.2969), Real(-0.1260), Real(-0.3516), Real(0.2843), Real(-0.1036)};
	return Co;
}
Real Naca(const Real& S)
{
	const Real* Co = NacaCo();
	return Co[0] * gm::sqrt(S) + S * (Co[1] + S * (Co[2] + S * (Co[3] + S * Co[4])));
}
const Real& NacaMax()
{
	static const Real M = [] {
		Real Best = Real::Infinity(true);
		const Real Step = Real(1) / 20000;
		for (int i = 0; i <= 20000; ++i)
		{
			Best = gm::max(Best, Naca(i == 20000 ? Real(1) : Real(i) * Step));
		}
		return Best;
	}();
	return M;
}
Real Shape(const Real& S) { return Naca(S) / NacaMax(); }
Real ShapeInt()
{
	const Real* Co = NacaCo();
	return (Co[0] * 2 / 3 + Co[1] / 2 + Co[2] / 3 + Co[3] / 4 + Co[4] / 5) / NacaMax();
}

Real Surf(const Table& T, const Real& X, const Real& Y, bool Upper)
{
	const Stn St = StationAt(T, Y);
	const Real S = gm::min(Real(1), gm::max(Real(0), (St.Xle - X) / (St.Xle - St.Xte))), H = Shape(S) * St.T;
	return Upper ? St.Zc + H * St.Up : St.Zc - H * (1 - St.Up);
}

Real SegInt(const Real& F0, const Real& F1, const Real& G0, const Real& G1, const Real& H)
{
	const Real Df = F1 - F0, Dg = G1 - G0;
	return H * (F0 * G0 + (F0 * Dg + G0 * Df) / 2 + Df * Dg / 3);
}

Result FnPlanform(const Call& C)
{
	Table T;
	std::string E;
	if (!TableArg(C, 0, 3, &T, &E))
	{
		return Result::Fail(E);
	}
	const int N = T.Rows;
	auto Ch = [&T](int i) { return T.At(i, 1) - T.At(i, 2); };
	Real A(0), Icc(0), Icy(0), Icx(0), Ict(0);
	for (int i = 0; i < N - 1; ++i)
	{
		const Real H = T.At(i + 1, 0) - T.At(i, 0);
		A += H * (Ch(i) + Ch(i + 1)) / 2;
		Icc += SegInt(Ch(i), Ch(i + 1), Ch(i), Ch(i + 1), H);
		Icy += SegInt(Ch(i), Ch(i + 1), T.At(i, 0), T.At(i + 1, 0), H);
		Icx += SegInt(Ch(i), Ch(i + 1), T.At(i, 1), T.At(i + 1, 1), H);
		if (T.Cols >= 4)
		{
			Ict += SegInt(Ch(i), Ch(i + 1), T.At(i, 3), T.At(i + 1, 3), H);
		}
	}
	const Real S = 2 * A, B = 2 * T.At(N - 1, 0), Mac = 2 * Icc / S;
	Real YEq = FlNan();
	for (int i = 0; i < N - 1; ++i)
	{
		if ((Ch(i) - Mac) * (Ch(i + 1) - Mac) <= 0 && Ch(i) != Ch(i + 1))
		{
			YEq = T.At(i, 0) + (Ch(i) - Mac) / (Ch(i) - Ch(i + 1)) * (T.At(i + 1, 0) - T.At(i, 0));
			break;
		}
	}
	return Result::Record({{"S", S}, {"b", B}, {"AR", B * B / S}, {"MAC", Mac}, {"ymac", 2 * Icy / S}, {"xmac", 2 * Icx / S}, {"taper", Ch(N - 1) / Ch(0)},
		{"cr", Ch(0)}, {"ct", Ch(N - 1)}, {"yeq", YEq}, {"vol", T.Cols >= 4 ? 2 * ShapeInt() * Ict : FlNan()}});
}

Result Sweep(const Call& C)
{
	Table T;
	std::string E;
	if (!TableArg(C, 1, 3, &T, &E))
	{
		return Result::Fail(E);
	}
	const Real X = C.N(0);
	const int I = C.Int(0) - 1;
	if (X != gm::floor(X) || I < 0 || I >= T.Rows - 1)
	{
		return Result::Fail("segment must be 1 to " + std::to_string(T.Rows - 1));
	}
	const Real H = T.At(I + 1, 0) - T.At(I, 0);
	return Result::Record({{"le", gm::atan(-(T.At(I + 1, 1) - T.At(I, 1)) / H)}, {"te", gm::atan(-(T.At(I + 1, 2) - T.At(I, 2)) / H)},
		{"area", 2 * H * ((T.At(I, 1) - T.At(I, 2)) + (T.At(I + 1, 1) - T.At(I + 1, 2))) / 2}});
}

Result PlanQuad(const Call& C)
{
	Table T;
	std::string E;
	if (!TableArg(C, 1, 3, &T, &E))
	{
		return Result::Fail(E);
	}
	const Real Nd = C.N(0);
	if (Nd != gm::floor(Nd) || Nd < 1 || Nd > 100000)
	{
		return Result::Fail("intervals n must be a whole number, 1 to 100 000");
	}
	const int N = C.Int(0);
	const Real Half = T.At(T.Rows - 1, 0), H = Half / N, HalfH = H / 2;
	Real A(0), Cc(0), Cy(0), Cx(0);
	for (int i = 0; i <= N; ++i)
	{
		const Real Y = i == N ? Half : Real(i) * H, Xl = Interp(T, 1, Y), Xt = Interp(T, 2, Y), Ch = Xl - Xt, W = (i == 0 || i == N) ? HalfH : H;
		const Real WCh = W * Ch;
		A += WCh;
		Cc += WCh * Ch;
		Cy += WCh * Y;
		Cx += WCh * Xl;
	}
	const Real S = 2 * A;
	return Result::Record({{"S", S}, {"MAC", 2 * Cc / S}, {"ymac", 2 * Cy / S}, {"xmac", 2 * Cx / S}});
}

Result FnStation(const Call& C)
{
	Table T;
	std::string E;
	if (!TableArg(C, 1, 6, &T, &E))
	{
		return Result::Fail(E);
	}
	const Stn S = StationAt(T, C.N(0));
	return Result::Record({{"xle", S.Xle}, {"xte", S.Xte}, {"c", S.Xle - S.Xte}, {"t", S.T}, {"zc", S.Zc}, {"up", S.Up}});
}

// ICAO standard atmosphere to 20 km, and the hydrostatic equation integrated by Simpson's rule as the check.
const Real& IsaR()
{
	static const Real R(287.05287);
	return R;
}

void IsaOf(const Real& H, Real* T, Real* P)
{
	static const Real Ex = FlG0() / (IsaR() * Real(0.0065)), P11 = 101325 * gm::pow(Real(216.65) / Real(288.15), Ex);
	if (H <= 11000)
	{
		*T = Real(288.15) - Real(0.0065) * H;
		*P = 101325 * gm::pow(*T / Real(288.15), Ex);
	}
	else
	{
		*T = Real(216.65);
		*P = P11 * gm::exp(-FlG0() * (H - 11000) / (IsaR() * *T));
	}
}

Result FnIsa(const Call& C)
{
	const Real H = C.N(0);
	if (H < 0 || H > 20000)
	{
		return Result::Fail("ISA here covers 0 to 20 000 m");
	}
	Real T, P;
	IsaOf(H, &T, &P);
	const Real Lapse(0.0065), T0(288.15), Top(11000), MinusG0R = -FlG0() / IsaR();
	auto Tf = [&](const Real& Z) { return T0 - Lapse * gm::min(Z, Top); };
	auto F = [&](const Real& Z) { return MinusG0R / Tf(Z); };
	Real Lp = gm::log(Real(101325)), Z(0);
	const int N = 2000;
	const Real Ends[2] = {H > Top ? Top : H, H};
	for (int Leg = 0; Leg < (H > Top ? 2 : 1); ++Leg)
	{
		const Real Dz = (Ends[Leg] - Z) / N, Dz6 = Dz / 6, HalfDz = Dz / 2;
		for (int i = 0; i < N; ++i)
		{
			Lp += Dz6 * (F(Z) + 4 * F(Z + HalfDz) + F(Z + Dz));
			Z += Dz;
		}
	}
	const Real Ph = gm::exp(Lp);
	return Result::Record({{"T", T}, {"p", P}, {"rho", P / (IsaR() * T)}, {"a", gm::sqrt(Real(1.4) * IsaR() * T)}, {"phyd", Ph}, {"rhohyd", Ph / (IsaR() * Tf(H))}});
}

const Real& DefaultArea()
{
	static const Real S = [] {
		const Table& T = DefaultTable();
		Real A(0);
		for (int i = 0; i < T.Rows - 1; ++i)
		{
			A += (T.At(i + 1, 0) - T.At(i, 0)) * ((T.At(i, 1) - T.At(i, 2)) + (T.At(i + 1, 1) - T.At(i + 1, 2))) / 2;
		}
		return 2 * A;
	}();
	return S;
}

Result Sizing(const Call& C)
{
	const Real Mass = C.N(0), H = C.N(1), Mach = C.N(2), Thrust = C.N(3), S = C.N(4).IsZero() ? DefaultArea() : C.N(4);
	if (H < 0 || H > 20000 || Mass <= 0 || S <= 0 || Mach <= 0)
	{
		return Result::Fail("needs mass > 0, 0 <= h <= 20 000 m, Mach > 0 and S > 0");
	}
	Real T, P;
	IsaOf(H, &T, &P);
	const Real Rho = P / (IsaR() * T), A = gm::sqrt(Real(1.4) * IsaR() * T), V = Mach * A, W = Mass * FlG0(), Q1 = Real(0.5) * Rho * V * V, Q2 = Real(0.5) * Real(1.4) * P * Mach * Mach;
	return Result::Record({{"V", V}, {"W", W}, {"q", Q1}, {"q2", Q2}, {"WS", W / S}, {"WSkg", Mass / S}, {"CL", W / (Q1 * S)}, {"TW", Thrust / W}});
}

// Cargo-box fit: golden section on the box's aft end (window ends bind), and a 0.05 m scan as the check.
Result Hold(const Call& C)
{
	Table T;
	std::string E;
	if (!TableArg(C, 2, 6, &T, &E))
	{
		return Result::Fail(E);
	}
	const Real L = C.N(0), Wd = C.N(1);
	if (L <= 0 || Wd <= 0 || Wd > 200)
	{
		return Result::Fail("hold length must be positive and width 0 to 200 m");
	}
	Real Xa = FlInf(), Xb = Real::Infinity(true);
	for (int i = 0; i < T.Rows; ++i)
	{
		Xa = gm::min(Xa, T.At(i, 2));
		Xb = gm::max(Xb, T.At(i, 1));
	}
	if (Xb - L <= Xa)
	{
		return Result::Fail("the hold is longer than the centre body allows");
	}
	if ((Xb - Xa) / Real(0.05) > 20000)
	{
		return Result::Fail("the body is too long for the 0.05 m scan (20 000 stations at most)");
	}
	const Real Y = Wd / 2;
	auto Hx = [&](const Real& X0) {
		return gm::min(Surf(T, X0, Y, true), Surf(T, X0 + L, Y, true)) - gm::max(Surf(T, X0, Y, false), Surf(T, X0 + L, Y, false));
	};
	Real A = Xa, B = Xb - L;
	const Real G = (gm::sqrt(Real(5)) - 1) / 2;
	for (int i = 0; i < 200; ++i)
	{
		const Real M1 = B - G * (B - A), M2 = A + G * (B - A);
		if (Hx(M1) < Hx(M2))
		{
			A = M1;
		}
		else
		{
			B = M2;
		}
	}
	const Real X0 = (A + B) / 2;
	// the scan
	const Real Dx(0.05), Tiny(1e-9), Gap(0.01);
	const int Ny = 25, N = JsRound((Xb - Xa) / Dx).ToInt();
	const Real Step = (Xb - Xa) / N;
	const int Nw = JsRound(L / Step).ToInt();
	std::vector<Real> Xs(static_cast<size_t>(N + 1)), Hi(static_cast<size_t>(N + 1)), Lo(static_cast<size_t>(N + 1));
	bool Sides = true;
	for (int i = 0; i <= N; ++i)
	{
		const Real X = i == N ? Xb : Xa + Real(i) * Step;
		Real U = FlInf(), Lw = Real::Infinity(true);
		for (int k = 0; k < Ny; ++k)
		{
			const Real Yk = -Wd / 2 + Wd * k / (Ny - 1);
			U = gm::min(U, Surf(T, X, Yk, true));
			Lw = gm::max(Lw, Surf(T, X, Yk, false));
		}
		const Real USide = Surf(T, X, Y, true), LSide = Surf(T, X, Y, false);
		Xs[static_cast<size_t>(i)] = X;
		Hi[static_cast<size_t>(i)] = U;
		Lo[static_cast<size_t>(i)] = Lw;
		if (U > Lw + Gap && (gm::fabs(USide - U) > Tiny || gm::fabs(LSide - Lw) > Tiny))
		{
			Sides = false;
		}
	}
	Real Best = Real::Infinity(true);
	int Bi = 0;
	if (Nw >= 1 && Nw <= N)
	{
		// Window minimum of the roof and maximum of the floor by monotone deques, as holdGrid.
		std::vector<int> DqU, DqL;
		size_t HeadU = 0, HeadL = 0;
		for (int i = 0; i <= N; ++i)
		{
			while (DqU.size() > HeadU && Hi[static_cast<size_t>(DqU.back())] >= Hi[static_cast<size_t>(i)])
			{
				DqU.pop_back();
			}
			DqU.push_back(i);
			while (DqL.size() > HeadL && Lo[static_cast<size_t>(DqL.back())] <= Lo[static_cast<size_t>(i)])
			{
				DqL.pop_back();
			}
			DqL.push_back(i);
			const int S0 = i - Nw;
			if (S0 < 0)
			{
				continue;
			}
			while (DqU[HeadU] < S0)
			{
				++HeadU;
			}
			while (DqL[HeadL] < S0)
			{
				++HeadL;
			}
			const Real Hh = Hi[static_cast<size_t>(DqU[HeadU])] - Lo[static_cast<size_t>(DqL[HeadL])];
			if (Hh > Best)
			{
				Best = Hh;
				Bi = S0;
			}
		}
	}
	return Result::Record({{"H", Hx(X0)}, {"x0", X0}, {"x1", X0 + L}, {"floor", gm::max(Surf(T, X0, Y, false), Surf(T, X0 + L, Y, false))},
		{"ceil", gm::min(Surf(T, X0, Y, true), Surf(T, X0 + L, Y, true))}, {"Hscan", Best}, {"xscan", Nw >= 1 && Nw <= N ? Xs[static_cast<size_t>(Bi)] : FlNan()}, {"sides", Flag(Sides)}});
}

// Body loft as the generator builds it (rings of 2M - 1 points, tip fans) and its divergence-theorem volume.
Result BodyVol(const Call& C)
{
	Table T;
	std::string E;
	if (!TableArg(C, 2, 6, &T, &E))
	{
		return Result::Fail(E);
	}
	const Real Md = C.N(0), Kd = C.N(1);
	if (Md != gm::floor(Md) || Md < 3 || Md > 300 || Kd != gm::floor(Kd) || Kd < 1 || Kd > 2)
	{
		return Result::Fail("needs a whole number of chordwise points M, 3 to 300, and spanwise multiplier k, 1 or 2");
	}
	const int M = C.Int(0), K = C.Int(1);
	if (T.Rows > 40)
	{
		return Result::Fail("at most 40 stations");
	}
	static const int SegSamples[4] = {10, 12, 10, 30};
	std::vector<Real> Half;
	for (int i = 0; i < T.Rows - 1; ++i)
	{
		const int Ns = SegSamples[std::min(i, 3)] * K;
		const Real A = T.At(i, 0), B = T.At(i + 1, 0), St = (B - A) / Ns;
		for (int j = 0; j < Ns; ++j)
		{
			Half.push_back(A + Real(j) * St);
		}
	}
	Half.push_back(T.At(T.Rows - 1, 0));
	std::vector<Real> Ys;
	for (const Real& V : Half)
	{
		Ys.push_back(V.IsZero() ? Real(0) : -V);
		Ys.push_back(V);
	}
	std::sort(Ys.begin(), Ys.end());
	Ys.erase(std::unique(Ys.begin(), Ys.end()), Ys.end());
	std::vector<Real> S(static_cast<size_t>(M)), Sh(static_cast<size_t>(M));
	const Real Dv = cst::Pi() / (M - 1);
	for (int i = 0; i < M; ++i)
	{
		const Real V = i == M - 1 ? cst::Pi() : Real(i) * Dv;
		S[static_cast<size_t>(i)] = (1 - gm::cos(V)) / 2;
		Sh[static_cast<size_t>(i)] = Shape(S[static_cast<size_t>(i)]);
	}
	// Ring order: trailing edge -> upper -> leading edge -> lower -> trailing edge.
	std::vector<int> Ring;
	for (int i = M - 1; i >= 0; --i)
	{
		Ring.push_back(i);
	}
	for (int i = 1; i < M; ++i)
	{
		Ring.push_back(i);
	}
	const int R = static_cast<int>(Ring.size()), NY = static_cast<int>(Ys.size());
	std::vector<Real> V;
	V.reserve(static_cast<size_t>(3 * (NY * R + 2)));
	for (const Real& Y : Ys)
	{
		const Stn St = StationAt(T, Y);
		const Real Chord = St.Xle - St.Xte, Lw = 1 - St.Up;
		for (int j = 0; j < R; ++j)
		{
			const size_t Q = static_cast<size_t>(Ring[static_cast<size_t>(j)]);
			const Real H = Sh[Q] * St.T;
			V.push_back(St.Xle - S[Q] * Chord);
			V.push_back(Y);
			V.push_back(j < M ? St.Zc + H * St.Up : St.Zc - H * Lw);
		}
	}
	auto Triple = [&V](int A, int B, int Cc) {
		const size_t a = static_cast<size_t>(3 * A), b = static_cast<size_t>(3 * B), c = static_cast<size_t>(3 * Cc);
		return V[a] * (V[b + 1] * V[c + 2] - V[b + 2] * V[c + 1]) - V[a + 1] * (V[b] * V[c + 2] - V[b + 2] * V[c]) + V[a + 2] * (V[b] * V[c + 1] - V[b + 1] * V[c]);
	};
	RealAccumulator Raw, Good;
	for (int i = 0; i < NY - 1; ++i)
	{
		for (int j = 0; j < R - 1; ++j)
		{
			const int A = i * R + j, B = A + 1, Cc = A + R, D = Cc + 1;
			const Real T1 = Triple(A, Cc, B), T2 = Triple(B, Cc, D);
			// Reversing a triangle's winding negates its triple product exactly.
			Raw.Add(T1);
			Raw.Add(T2);
			Good.Sub(T1);
			Good.Sub(T2);
		}
	}
	for (int Cap = 0; Cap < 2; ++Cap)
	{
		const int Base = Cap == 0 ? 0 : (NY - 1) * R, Ci = static_cast<int>(V.size() / 3);
		RealAccumulator Cx, Cy, Cz;
		for (int j = 0; j < R; ++j)
		{
			Cx.Add(V[static_cast<size_t>(3 * (Base + j))]);
			Cy.Add(V[static_cast<size_t>(3 * (Base + j) + 1)]);
			Cz.Add(V[static_cast<size_t>(3 * (Base + j) + 2)]);
		}
		V.push_back(Cx.Result() / R);
		V.push_back(Cy.Result() / R);
		V.push_back(Cz.Result() / R);
		for (int j = 0; j < R - 1; ++j)
		{
			const Real Tv = Cap == 0 ? Triple(Base + j + 1, Base + j, Ci) : Triple(Ci, Base + j, Base + j + 1);
			Raw.Add(Tv);
			Good.Add(Tv);
		}
	}
	Real Ict(0);
	for (int i = 0; i < T.Rows - 1; ++i)
	{
		const Real C0 = T.At(i, 1) - T.At(i, 2), C1 = T.At(i + 1, 1) - T.At(i + 1, 2);
		Ict += SegInt(C0, C1, T.At(i, 3), T.At(i + 1, 3), T.At(i + 1, 0) - T.At(i, 0));
	}
	return Result::Record({{"raw", Raw.Result() / 6}, {"good", Good.Result() / 6}, {"exact", 2 * ShapeInt() * Ict}});
}

/* ---------------------------------------------------------------- game.js: par is a theorem */
// Level 1: dead-stick glide on the sheet's polar from 1000 m.
Result GlidePar(const Call& C)
{
	const Real H = C.N(0), CD0 = C.N(1), AR = C.N(2), E = C.N(3);
	if (CD0 <= 0 || AR <= 0 || E <= 0)
	{
		return Result::Fail("CD0, AR and e must be positive");
	}
	const Real K = 1 / (cst::Pi() * E * AR);
	return Result::Of(H * (Real(0.5) / gm::sqrt(CD0 * K)));
}

Result GlideGame(const Call& C)
{
	const Real Dt = C.N(1);
	if (!(Dt >= Real(0.001)) || Dt > 10)
	{
		return Result::Fail("time step must be 0.001 to 10 s");
	}
	const Real W(10000), S(16), CD0(0.025), Rho(1.225), K = KSheet();
	const Real Cc = gm::min(Real(1.4), gm::max(Real(0.2), C.N(0))), CD = CD0 + K * Cc * Cc, Gam = gm::atan(CD / Cc);
	const Real V = gm::sqrt(2 * W * gm::cos(Gam) / (Rho * S * Cc)), Dz0 = V * gm::sin(Gam) * Dt;
	Real X(0), Z(1000), T(0);
	bool Done = false;
	for (int n = 0; !Done && n < 1000000; ++n)
	{
		Real Dz = Dz0, H = Dt;
		if (Dz >= Z)
		{
			H = Dt * Z / Dz;
			Dz = Z;
			Done = true;
		}
		X += Dz * Cc / CD;
		Z = Done ? Real(0) : Z - Dz;
		T += H;
	}
	return Result::Record({{"score", X}, {"t", T}, {"CL", Cc}});
}

// Level 2: hold 5 deg against a 0.5 rad/s^2 disturbance on the sheet's short-period model.
struct PitchConst
{
	Real ThetaC, Dist, Dt, Rho, Deg;
	int Steps;
};

const PitchConst& Ph()
{
	static const PitchConst K{5 * cst::Pi() / 180, Real(0.5), Real(1e-3), Real(0.01), 180 / cst::Pi(), 20000};
	return K;
}

void PitchTrimOf(Real* Alpha, Real* Delta)
{
	const SpModel& P = SpDefault();
	const Real Det = P.Za * P.Md - P.Zd * P.Ma, D = -P.Za * Ph().Dist / Det;
	*Delta = D;
	*Alpha = -P.Zd * D / P.Za;
}

Result PitchTrim(const Call& C)
{
	(void)C;
	Real A, D;
	PitchTrimOf(&A, &D);
	return Result::Record({{"alpha", A}, {"delta", D}});
}

Result PitchGame(const Call& C)
{
	const Real Kp = C.N(0), Ki = C.N(1), Kd = C.N(2), Ff = C.N(3);
	const SpModel& P = SpDefault();
	const PitchConst& K = Ph();
	Real TrA, TrD;
	PitchTrimOf(&TrA, &TrD);
	Real A(0), Q(0), Th(0), Z(0), Ise(0), Ctl(0);
	const Real Scale = K.Deg * K.Deg * K.Dt, Lim(1000);
	for (int k = 0; k < K.Steps; ++k)
	{
		const Real E = K.ThetaC - Th, D = Ff - (Kp * E + Ki * Z) + Kd * Q, U = D - TrD;
		Ise += E * E * Scale;
		Ctl += U * U * Scale;
		const Real Da = P.Za * A + Q + P.Zd * D, Dq = P.Ma * A + P.Mq * Q + P.Md * D + K.Dist;
		A += K.Dt * Da;
		Q += K.Dt * Dq;
		Th += K.Dt * Q;
		Z += K.Dt * E;
		if (!gm::isfinite(Th) || gm::fabs(Th) > Lim)
		{
			return Result::Record({{"J", FlInf()}, {"ise", FlInf()}, {"ctl", FlInf()}});
		}
	}
	return Result::Record({{"J", Ise + K.Rho * Ctl}, {"ise", Ise}, {"ctl", Ctl}});
}

// Par: finite-horizon discrete LQR (Riccati recursion) on the simulation's exact one-step map, computed once.
Result PitchPar(const Call& C)
{
	(void)C;
	static const Real Value = [] {
		const SpModel& P = SpDefault();
		const PitchConst& K = Ph();
		const Real H = K.Dt, Zero(0), One(1);
		const Real Ad[3][3] = {{1 + H * P.Za, H, Zero}, {H * P.Ma, 1 + H * P.Mq, Zero}, {H * H * P.Ma, H + H * H * P.Mq, One}};
		const Real Bd[3] = {H * P.Zd, H * P.Md, H * H * P.Md};
		const Real Qw = K.Deg * K.Deg * K.Dt, Rw = K.Rho * K.Deg * K.Deg * K.Dt;
		Real TrA, TrD;
		PitchTrimOf(&TrA, &TrD);
		Real S[3][3] = {{Zero, Zero, Zero}, {Zero, Zero, Zero}, {Zero, Zero, Zero}};
		for (int k = K.Steps - 1; k >= 0; --k)
		{
			Real SB[3], Kk[3], Cm[3][3], NS[3][3];
			for (int r = 0; r < 3; ++r)
			{
				SB[r] = S[r][0] * Bd[0] + S[r][1] * Bd[1] + S[r][2] * Bd[2];
			}
			const Real Den = Rw + Bd[0] * SB[0] + Bd[1] * SB[1] + Bd[2] * SB[2];
			for (int j = 0; j < 3; ++j)
			{
				Kk[j] = (SB[0] * Ad[0][j] + SB[1] * Ad[1][j] + SB[2] * Ad[2][j]) / Den;
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					Cm[i][j] = Ad[i][j] - Bd[i] * Kk[j];
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					Real Sum(0);
					for (int m = 0; m < 3; ++m)
					{
						Sum += Ad[m][i] * (S[m][0] * Cm[0][j] + S[m][1] * Cm[1][j] + S[m][2] * Cm[2][j]);
					}
					NS[i][j] = (i == 2 && j == 2) ? Sum + Qw : Sum;
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					S[i][j] = NS[i][j];
				}
			}
		}
		const Real X0[3] = {-TrA, Zero, -K.ThetaC};
		Real V(0);
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				V += X0[i] * S[i][j] * X0[j];
			}
		}
		return V;
	}();
	return Result::Of(Value);
}

// Level 3: first track. Airport radar at 2.8 GHz against a 1 m^2 aircraft inbound at 120 kt from 60 nmi.
const int RaPower[5][2] = {{25, 1}, {50, 2}, {100, 3}, {200, 5}, {400, 8}};   // kW, cost
const int RaGain[5][2] = {{28, 1}, {30, 2}, {32, 3}, {34, 5}, {36, 8}};       // dB, cost
const int RaScan[6] = {2, 4, 6, 8, 10, 12};                                   // s

Real RdErfc(const Real& X)   // Numerical Recipes erfcc, as radar.js
{
	static const Real Co[10] = {Real(-1.26551223), Real(1.00002368), Real(0.37409196), Real(0.09678418), Real(-0.18628806), Real(0.27886807), Real(-1.13520398), Real(1.48851587),
		Real(-0.82215223), Real(0.17087277)};
	const Real Z = gm::fabs(X), T = 1 / (1 + Real(0.5) * Z);
	Real Poly = Co[9];
	for (int i = 8; i >= 1; --i)
	{
		Poly = Co[i] + T * Poly;
	}
	const Real R = T * gm::exp(-Z * Z + Co[0] + T * Poly);
	return X >= 0 ? R : 2 - R;
}
Real RdGammaQ(int N, const Real& T)
{
	Real S(0), Term(1);
	for (int k = 0; k < N; ++k)
	{
		if (k)
		{
			Term *= T / k;
		}
		S += Term;
	}
	return gm::exp(-T) * S;
}
Real RdThreshold(int N, const Real& Pfa)
{
	Real Lo(0), Hi(50 + 10 * N);
	for (int i = 0; i < 200; ++i)
	{
		const Real M = (Lo + Hi) / 2;
		if (M == Lo || M == Hi)
		{
			break;   // the bracket is down to adjacent 18-digit values; more halvings move it by one digit at most
		}
		if (RdGammaQ(N, M) > Pfa)
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	return (Lo + Hi) / 2;
}

// One design's detector: n pulses, threshold T, and the prefix sums of T^k/k! that gammaQ(m, T) is built from,
// extended as the Poisson window grows (the same partial sums gammaQ forms term by term).
struct Detector
{
	int N = 0;
	Real T, ExpT, Term;
	std::vector<Real> Prefix;   // Prefix[m] = sum_{k < m} T^k / k!
	void Extend(int Upto)
	{
		while (static_cast<int>(Prefix.size()) <= Upto)
		{
			const int K = static_cast<int>(Prefix.size()) - 1;
			if (K)
			{
				Term *= T / K;
			}
			Prefix.push_back(Prefix.back() + Term);
		}
	}
};

// Steady target, n pulses integrated noncoherently: Poisson mixture of central chi-squares; Poisson weights by
// recurrence from the first one.
Real RdPdSteady(const Real& Snr, Detector* D)
{
	const int N = D->N;
	const Real Lam = N * Snr, Sd = gm::sqrt(Lam);
	if (Lam > 5000)
	{
		const Real Z = (N + Lam - D->T) / gm::sqrt(N + 2 * Lam);
		return Z > 8 ? Real(1) : Real(0.5) * RdErfc(-Z / gm::sqrt(Real(2)));
	}
	const int J0 = gm::max(Real(0), gm::floor(Lam - 12 * Sd - 20)).ToInt(), J1 = gm::ceil(Lam + 12 * Sd + 40).ToInt();
	D->Extend(N + J1);
	if (Lam.IsZero())
	{
		return gm::min(Real(1), D->ExpT * D->Prefix[static_cast<size_t>(N)]);
	}
	const Real Lf = J0 >= 2 ? gm::lgamma(Real(J0 + 1)) : Real(0);
	Real W = J0 == 0 ? gm::exp(-Lam) : gm::exp(-Lam + J0 * gm::log(Lam) - Lf);
	// Past twice the mean the weights fall by more than half per term, so once a term is below 1e-20 of the
	// sum the rest of radar.js's window cannot reach the 18th digit and is left off.
	Real S(0);
	const Real Floor(1e-300), Negligible(1e-20), TwiceMean = 2 * Lam + 1;
	for (int j = J0; j <= J1; ++j)
	{
		if (j > J0)
		{
			W = W * Lam / j;
		}
		if (W > Floor)
		{
			const Real Term = W * D->Prefix[static_cast<size_t>(N + j)];
			S += Term;
			if (Real(j) > TwiceMean && Term < Negligible * S)
			{
				break;
			}
		}
	}
	return gm::min(Real(1), D->ExpT * S);
}

struct RadarOut
{
	Real Score;
	int K = -1, N = 0, Cost = 0;
};

// The detector depends on the gain (beamwidth) and scan period only, so the 30 of them are built once, on
// first use, and shared by every power level, by radargame and by radarpar.
Detector& DetectorFor(int Gi, int Si)
{
	static Detector Cache[5][6];
	Detector& D = Cache[Gi][Si];
	if (D.Prefix.empty())
	{
		const Real GDb(RaGain[Gi][0]), T(RaScan[Si]);
		const Real Bw = Real(1.4) * gm::pow(Real(10), (33 - GDb) / 10);
		D.N = gm::max(Real(1), gm::floor(1200 * T * Bw / 360)).ToInt();
		D.T = RdThreshold(D.N, Real(1e-6));
		D.ExpT = gm::exp(-D.T);
		D.Term = Real(1);
		D.Prefix.push_back(Real(0));
	}
	return D;
}

// Bound > 0 (nmi): give up (Score -1) once the range falls below Bound, since the design cannot reach it.
RadarOut RadarRun(int Pi, int Gi, int Si, const Real& Bound)
{
	const Real Kw(RaPower[Pi][0]), GDb(RaGain[Gi][0]), T(RaScan[Si]), Nmi(1852), Kt = Real(1852) / 3600;
	const Real Lam = cst::C() / Real(2.8e9);
	const Real G = gm::pow(Real(10), GDb / 10), L = gm::pow(Real(10), Real(8) / 10);
	Detector& D = DetectorFor(Gi, Si);
	const int N = D.N;
	// S/N = Pt G^2 lambda^2 sigma / ((4 pi)^3 R^4 k Ts Bn L), sigma = 1 m^2, Ts = 950 K, Bn = 1.67 MHz.
	const Real FourPi = 4 * cst::Pi();
	const Real Num = Kw * 1000 * G * G * Lam * Lam, Den = FourPi * FourPi * FourPi * cst::Boltzmann() * 950 * Real(1.67e6) * L;
	const Real Step = 120 * Kt * T, R0 = 60 * Nmi, Goal(0.9), BoundM = Bound * Nmi;
	RadarOut Out{Real(0), -1, N, RaPower[Pi][1] + RaGain[Gi][1]};
	Real St[3] = {Real(1), Real(0), Real(0)}, Got(0);
	for (int k = 0;; ++k)
	{
		const Real Rm = R0 - Step * k;
		if (Rm <= 0)
		{
			break;
		}
		if (Rm < BoundM)
		{
			Out.Score = Real(-1);
			break;
		}
		const Real R2 = Rm * Rm, Snr = Num / (Den * (R2 * R2));
		const Real P = RdPdSteady(Snr, &D), Q = 1 - P;
		// 2 of 3: node 1 = one hit, node 2 = hit then miss, node 3 = track.
		const Real N0 = St[0] * Q + St[2] * Q, N1 = St[0] * P, N2 = St[1] * Q, N3 = St[1] * P + St[2] * P;
		Got += N3;
		St[0] = N0;
		St[1] = N1;
		St[2] = N2;
		if (Got >= Goal)
		{
			Out.Score = Rm / Nmi;
			Out.K = k;
			break;
		}
	}
	return Out;
}

int FindIn(const int (*List)[2], int N, const Real& X)
{
	for (int i = 0; i < N; ++i)
	{
		if (X == Real(List[i][0]))
		{
			return i;
		}
	}
	return -1;
}

Result RadarGame(const Call& C)
{
	const int Pi = FindIn(RaPower, 5, C.N(0)), Gi = FindIn(RaGain, 5, C.N(1));
	int Si = -1;
	for (int i = 0; i < 6; ++i)
	{
		Si = C.N(2) == Real(RaScan[i]) ? i : Si;
	}
	if (Pi < 0 || Gi < 0 || Si < 0)
	{
		return Result::Fail("choose power 25, 50, 100, 200 or 400 kW; gain 28, 30, 32, 34 or 36 dB; scan 2 to 12 s in steps of 2");
	}
	const RadarOut R = RadarRun(Pi, Gi, Si, Real(0));
	return Result::Record({{"range", R.Score}, {"cost", Real(R.Cost)}, {"n", Real(R.N)}, {"ok", Flag(R.Cost <= 10)}, {"scan", Real(R.K)}});
}

// Par: exhaustive search of the admissible designs (cost <= 10), computed once. Strong designs are flown
// first, and a later design is abandoned once its range drops below the best so far (it cannot win); ties
// are flown out and the winner is the first best in radar.js's order (power, gain, scan), as radarPar.
Result RadarPar(const Call& C)
{
	(void)C;
	struct Best
	{
		Real Score = Real(-1);
		int Kw = 0, GDb = 0, T = 0, N = 0;
	};
	static const Best B = [] {
		RadarOut Runs[5][5][6];
		Real Bound(0);
		for (int Pi = 4; Pi >= 0; --Pi)
		{
			for (int Gi = 4; Gi >= 0; --Gi)
			{
				for (int Si = 0; Si < 6; ++Si)
				{
					Runs[Pi][Gi][Si].Score = Real(-1);
					if (RaPower[Pi][1] + RaGain[Gi][1] <= 10)
					{
						Runs[Pi][Gi][Si] = RadarRun(Pi, Gi, Si, Bound);
						Bound = gm::max(Bound, Runs[Pi][Gi][Si].Score);
					}
				}
			}
		}
		Best Out;
		for (int Pi = 0; Pi < 5; ++Pi)
		{
			for (int Gi = 0; Gi < 5; ++Gi)
			{
				for (int Si = 0; Si < 6; ++Si)
				{
					const RadarOut& R = Runs[Pi][Gi][Si];
					if (R.Score > Out.Score)
					{
						Out.Score = R.Score;
						Out.Kw = RaPower[Pi][0];
						Out.GDb = RaGain[Gi][0];
						Out.T = RaScan[Si];
						Out.N = R.N;
					}
				}
			}
		}
		return Out;
	}();
	return Result::Record({{"range", B.Score}, {"kw", Real(B.Kw)}, {"gdb", Real(B.GDb)}, {"T", Real(B.T)}, {"n", Real(B.N)}});
}

Result Stars(const Call& C)
{
	const Real Score = C.N(0), Par = C.N(1);
	const bool Higher = !C.N(2).IsZero();
	const Real R = !gm::isfinite(Score) ? Real(0) : Higher ? Score / Par : (Score > 0 ? Par / Score : Real(0));
	return Result::Of(R >= Real(0.99) ? Real(3) : R >= Real(0.95) ? Real(2) : R >= Real(0.85) ? Real(1) : Real(0));
}

Param TableParam()
{
	Param P = Mat("table", "m");
	P.HasDefault = true;
	P.Default = Real(0);
	return P;
}
} // namespace flight_detail

void RegisterFlightLibrary(Registry& R)
{
	using namespace flight_detail;
	const Real A4 = FlRad(Real(4)), KS = KSheet();
	const char* Sp = "1/s";

	// Equation sheet (flight-model/equations.js, equations.html)
	Def(R, "thinairfoil", "aero", "flight-model/equations.js thinAirfoil (Glauert)", "thin-airfoil theory for a NACA 4-digit mean line", "record aL0 (rad), cl, cm4, A1, A2",
		{Num("m", "-", 0.02), Num("p", "-", 0.4), Num("alpha", "rad", A4)}, &ThinAirfoil);
	Def(R, "vortexpanel", "aero", "flight-model/equations.js vortexPanels", "discrete-vortex check: vortex at 1/4, tangency at 3/4 of each panel (N 2 to 200)", "record cl, cm4",
		{Num("m", "-", 0.02), Num("p", "-", 0.4), Num("alpha", "rad", A4), Num("N", "-", 200)}, &VortexPanel);
	Def(R, "pg", "aero", "flight-model/equations.js prandtlGlauert", "Prandtl-Glauert: c0 / sqrt(1 - M^2)", "corrected coefficient", {Num("c0", "-", 0.5), Num("M", "-", 0.6)}, &PrandtlGlauert);
	Def(R, "karmantsien", "aero", "flight-model/equations.js karmanTsien", "Karman-Tsien pressure-coefficient correction", "Cp", {Num("cp0", "-", -0.5), Num("M", "-", 0.6)}, &KarmanTsien);
	Def(R, "ackeret", "aero", "flight-model/equations.js ackeretCl", "Ackeret supersonic lift: 4 alpha / sqrt(M^2 - 1)", "cl", {Num("alpha", "rad", A4), Num("M", "-", 2)}, &Ackeret);
	Def(R, "dragpolar", "aero", "flight-model/equations.js polar", "parabolic drag polar, best L/D and its CL", "record k, LDmax, CLstar, LD (at CL)",
		{Num("CD0", "-", 0.025), Num("AR", "-", 8), Num("e", "-", 0.8), Num("CL", "-", 0.5)}, &Polar);
	Def(R, "spaneff", "aero", "flight-model/equations.js liftingLineE + taperChord", "lifting-line span efficiency of an untwisted tapered wing", "record e, delta",
		{Num("AR", "-", 8), Num("lambda", "-", 0.4), Num("terms", "-", 30)}, &SpanEff);
	Def(R, "spaneff_ell", "aero", "flight-model/equations.js liftingLineE + ellipticChord", "lifting-line span efficiency of an elliptic wing (e = 1)", "record e, delta",
		{Num("AR", "-", 8), Num("terms", "-", 30)}, &SpanEffEll);
	Def(R, "wavedrag", "aero", "flight-model/equations.js waveDq + searsHaackDq", "Sears-Haack wave drag D/q: Fourier series and 128 V^2/(pi l^4)", "record fourier, closed (m^2), ratio",
		{Num("rmax", "m", 0.05), Num("l", "m", 1)}, &WaveDrag);
	Def(R, "arearule", "aero", "flight-model/equations.html area-rule card", "Sears-Haack fuselage plus a wing's area, with an indent", "record fus (D/q, m^2), wing, ruled (x fuselage), closed",
		{Num("wing", "-", 0.4), Num("indent", "-", 0)}, &AreaRule);
	Def(R, "arearatio", "aero", "flight-model/equations.js areaRatio", "isentropic area-Mach relation A/A*", "A/A*", {Num("M", "-", 2), Num("gamma", "-", 1.4)}, &AreaRatio);
	Def(R, "machfromarea", "aero", "flight-model/equations.js machFromArea", "Mach for an area ratio: 1 = supersonic root, 0 = subsonic", "Mach",
		{Num("ratio", "-"), Num("sup", "-", 1), Num("gamma", "-", 1.4)}, &MachFromArea);
	Def(R, "climbangle", "aero", "flight-model/equations.js climbAngle", "steady climb angle from T/W = sin g + cos g / (L/D)", "rad", {Num("TW", "-", 0.18), Num("LD", "-", 12)}, &ClimbAngle);
	Def(R, "rocketdv", "aero", "flight-model/equations.js rocketDv", "rocket equation Isp g0 ln(m0/mf)", "m/s", {Num("Isp", "s", 300), Num("m0", "kg", 4000), Num("mf", "kg", 1000)}, &RocketDv);
	Def(R, "turboprop", "aero", "flight-model/equations.js turboprop", "turboprop shaft power, simple Brayton cycle", "record tau, P (W), turbine, compressor, q (W), eta, ideal",
		{Num("mdot", "kg/s", 5), Num("pic", "-", 12), Num("T02", "K", 288.15), Num("T04", "K", 1400), Num("etac", "-", 0.85), Num("etat", "-", 0.9), Num("gamma", "-", 1.4), Num("cp", "J/(kg K)", cst::CpAir())},
		&Turboprop);
	Def(R, "netthrust", "aero", "flight-model/equations.js netThrust", "net thrust = gross jet thrust - ram drag", "record gross, ram, net (N)",
		{Num("mdot", "kg/s", 100), Num("f", "-", 0.02), Num("Ve", "m/s", 600), Num("V0", "m/s", 250), Num("pe", "Pa", 0), Num("p0", "Pa", 0), Num("Ae", "m^2", 0)}, &NetThrust);
	Def(R, "wingloads", "aero", "flight-model/equations.js wingLoads", "elliptic-wing running load, shear and bending at station y", "record w0, w (N/m), V (N), M (N m)",
		{Num("nW", "N", 25000), Num("s", "m", 5), Num("y", "m", 0)}, &WingLoads);
	Def(R, "bredt", "aero", "flight-model/equations.js bredtTube, exactTube, openTube", "Bredt-Batho tube vs exact thick tube vs slit tube", "record q (N/m), bredt, exact, open (Pa), twist (rad/m)",
		{Num("T", "N m", 1000), Num("r", "m", 0.1), Num("t", "m", 0.002), Num("G", "Pa", 27e9)}, &Torsion);
	Def(R, "buckle", "aero", "flight-model/equations.js eulerLoad + fdBuckle", "Euler buckling load; finite-difference check for K = 1 or 2", "record P, fd (N)",
		{Num("E", "Pa", 70e9), Num("I", "m^4", 1e-8), Num("L", "m", 1), Num("K", "-", 1)}, &Buckle);
	Def(R, "bestclimb", "aero", "flight-model/equations.js climbCurve", "best rate of climb by a 0.1 m/s speed scan, and the minimum-power speed", "record Vbest (m/s), rocmax (m/s), Vmp (m/s)",
		{Num("W", "N", 10000), Num("S", "m^2", 16), Num("rho", "kg/m^3", 1.225), Num("Pav", "W", 60000), Num("CD0", "-", 0.025), Num("k", "-", KS)}, &BestClimb);
	Def(R, "levelturn", "aero", "flight-model/equations.js turn", "level coordinated turn: load factor, rate and radius", "record n, omega (rad/s), R (m)", {Num("V", "m/s"), Num("phi", "rad")}, &FnTurn);
	Def(R, "cornerspeed", "aero", "flight-model/equations.js cornerSpeed", "corner speed sqrt(2 nmax W / (rho S CLmax))", "m/s",
		{Num("nmax", "-", 3.8), Num("W", "N", 10000), Num("rho", "kg/m^3", 1.225), Num("S", "m^2", 16), Num("CLmax", "-", 1.5)}, &CornerSpeed);
	Def(R, "glide", "aero", "flight-model/equations.js glide", "still-air glide at a held CL", "record gam (rad), V, ratio, dist (m), sink (m/s)",
		{Num("h", "m", 1000), Num("W", "N", 10000), Num("S", "m^2", 16), Num("rho", "kg/m^3", 1.225), Num("CD0", "-", 0.025), Num("k", "-", KS), Num("CL", "-", 0.71)}, &FnGlide);
	const std::vector<Param> SpParams = {Num("Za", Sp, -2.0), Num("Zd", Sp, -0.2), Num("Ma", "1/s^2", -12.0), Num("Mq", Sp, -3.0), Num("Md", "1/s^2", -15.0)};
	Def(R, "shortperiod", "aero", "flight-model/equations.js shortPeriod", "short-period model: theta/delta coefficients, wn, zeta, damped period",
		"record wn (rad/s), zeta, a1, a0 (den s^2 + a1 s + a0), b1, b0 (num), Td (s)", SpParams, &ShortPeriod);
	{
		std::vector<Param> P = {Num("w", "rad/s", 1)};
		P.insert(P.end(), SpParams.begin(), SpParams.end());
		Def(R, "pitchtf", "aero", "flight-model/equations.js tfTheta + ssTheta", "theta/delta at s = jw: transfer function, and state space as the check",
			"record re, im, mag, db, ph (rad), ssre, ssim", P, &PitchTf);
	}
	{
		std::vector<Param> P = {Num("Kp", "-", 1.5), Num("dist", "rad/s^2", 0.5)};
		P.insert(P.end(), SpParams.begin(), SpParams.end());
		Def(R, "poerror", "aero", "flight-model/equations.js pOnlyError", "P-only steady pitch error: dist / (Kp (Md - Ma Zd / Za))", "rad", P, &POnlyError);
	}
	Def(R, "pidsim", "aero", "flight-model/equations.js pidSim", "PID pitch-attitude hold on the short-period model, 20 s", "record err (rad), overshoot, diverged",
		{Num("Kp", "-", 1.5), Num("Ki", "-", 1), Num("Kd", "-", 0.3), Num("thetac", "rad", FlRad(Real(5))), Num("dist", "rad/s^2", 0.5)}, &PidSim);

	// Atmosphere, flight model and clouds (flight-model/model.js)
	const char* AcU = "1-3";
	Def(R, "atmos", "atmos", "flight-model/model.js atmosphere (NASA Glenn)", "NASA Glenn atmosphere with a +dT day; icing = 1 for +2 to -20 C",
		"record T (K), Tc (C), p (Pa), rho (kg/m^3), a (m/s), sigma, icing", {Num("h", "m"), Num("dT", "K", 0)}, &FnAtmos);
	Def(R, "liftslope", "atmos", "flight-model/model.js sectionSlope + wingSlope", "2 pi slope with Prandtl-Glauert/Ackeret blend; Helmbold finite wing", "record a0, a (1/rad)",
		{Num("M", "-", 0), Num("AR", "-", 0)}, &LiftSlope);
	Def(R, "liftcoef", "atmos", "flight-model/model.js liftCoefficient", "lift coefficient over the full angle range, stall blend", "record CL, CLa, linear, stalled",
		{Num("alpha", "rad"), Num("M", "-", 0.1), Num("ac", AcU, 1)}, &LiftCoef);
	Def(R, "acaero", "atmos", "flight-model/model.js aero", "lift and drag of aircraft ac (1 GA, 2 jet, 3 SST)", "record L, D (N), CL, CD, CD0, CDi, K, M, q (Pa)",
		{Num("alpha", "rad"), Num("V", "m/s"), Num("h", "m"), Num("ac", AcU, 1), Num("beta", "rad", 0)}, &FnAero);
	Def(R, "acthrust", "atmos", "flight-model/model.js thrust", "lumped thrust lapse: propeller or jet", "N",
		{Num("throttle", "-"), Num("V", "m/s"), Num("h", "m"), Num("ac", AcU, 1)}, &FnThrust);
	Def(R, "actrim", "atmos", "flight-model/model.js trim", "joint trim of alpha and throttle for steady flight (throttle nan if full power is short)", "record alpha (rad), throttle, T (N)",
		{Num("V", "m/s"), Num("h", "m"), Num("ac", AcU, 1), Num("gam", "rad", 0), Num("phi", "rad", 0)}, &FnTrim);
	Def(R, "leveldrag", "atmos", "flight-model/model.js levelDrag + glideAngle", "level-flight drag at V; gam is the glide angle atan(D/L)", "record D (N), alpha, CL, CD, M, gam (rad), ok",
		{Num("V", "m/s"), Num("h", "m"), Num("ac", AcU, 1)}, &LevelDrag);
	Def(R, "takeoff", "atmos", "flight-model/model.js step (ground roll), test takeoff", "takeoff from rest: RK4 ground roll at 1/120 s, rotate past Vrot", "record t (s), dist (m), V (m/s), alpha (rad)",
		{Num("ac", AcU, 1), Num("throttle", "-", 1), Num("elev", "-", 0.5)}, &Takeoff);
	Def(R, "moistlapse", "atmos", "flight-model/model.js moistLapse", "saturated adiabatic lapse rate (dry rate g/cp = 0.00977)", "K/m", {Num("T", "K"), Num("p", "Pa", 1e5)}, &MoistLapse);
	Def(R, "satvap", "atmos", "flight-model/model.js satVapour + satMixing (Bolton)", "saturation vapour pressure and mixing ratio", "record e (Pa), rs (kg/kg)",
		{Num("Tc", "C"), Num("p", "Pa", 101325)}, &SatVap);
	Def(R, "dewpoint", "atmos", "flight-model/model.js dewPoint", "dew point from mixing ratio and pressure", "C", {Num("r", "kg/kg"), Num("p", "Pa", 101325)}, &DewPoint);
	Def(R, "sounding", "atmos", "flight-model/model.js sounding", "rising parcel: cloud base, free convection, cloud top, CAPE, CIN", "record lcl, lfc, el (m), cape, cin (J/kg), conv",
		{Num("dT", "K", 0), Num("spread", "K", 6), Num("heating", "K", 1.5)}, &FnSounding);
	Def(R, "updraft", "atmos", "flight-model/model.js updraftCore", "cumulus updraft core 0.4 sqrt(2 CAPE(z)) + 1.5, capped at 25", "m/s",
		{Num("z", "m"), Num("dT", "K", 8), Num("spread", "K", 6), Num("heating", "K", 2)}, &Updraft);
	Def(R, "gustload", "atmos", "flight-model/model.js derivatives (gust kinematics)", "sharp vertical gust: alpha rises by w/V; load factor from level trim", "record dalpha (rad), n, dn, alpha0 (rad)",
		{Num("w", "m/s"), Num("V", "m/s"), Num("h", "m"), Num("ac", AcU, 1)}, &GustLoad);

	// Planform bench (flight-model/planform.js, planform.html)
	Def(R, "planform", "airframe", "flight-model/planform.js planform", "planform integrals of the station table (default: the 62 m airlifter)",
		"record S (m^2), b, AR, MAC, ymac, xmac, taper, cr, ct, yeq (m), vol (m^3)", {TableParam()}, &FnPlanform);
	Def(R, "wingsweep", "airframe", "flight-model/planform.js planform seg", "leading- and trailing-edge sweep and area of segment i", "record le, te (rad), area (m^2, both sides)",
		{Num("i", "-"), TableParam()}, &Sweep);
	Def(R, "planquad", "airframe", "flight-model/planform.js planformQuad", "planform check: chord at n + 1 even stations, trapezoidal rule (n up to 100 000)", "record S, MAC, ymac, xmac",
		{Num("n", "-", 1000), TableParam()}, &PlanQuad);
	Def(R, "wingstation", "airframe", "flight-model/planform.js station", "interpolated station at span y", "record xle, xte, c, t, zc (m), up", {Num("y", "m"), TableParam()}, &FnStation);
	Def(R, "stdatm", "airframe", "flight-model/planform.js isa + hydro", "ICAO standard atmosphere 0-20 km, with the hydrostatic integral", "record T, p, rho, a, phyd, rhohyd",
		{Num("h", "m")}, &FnIsa);
	Def(R, "airsizing", "airframe", "flight-model/planform.js sizing", "cruise sizing: speed, q, wing loading, CL, T/W (S = 0: the planform's area)", "record V, W, q, q2, WS, WSkg, CL, TW",
		{Num("mass", "kg", 200e3), Num("h", "m", 10e3), Num("M", "-", 0.8), Num("thrust", "N", 560e3), Num("S", "m^2", 0)}, &Sizing);
	Def(R, "cargohold", "airframe", "flight-model/planform.js holdEnds + holdGrid", "largest clear height for an L x w cargo box inside the body", "record H, x0, x1, floor, ceil, Hscan, xscan (m), sides",
		{Num("L", "m", 19), Num("w", "m", 6), TableParam()}, &Hold);
	Def(R, "bodyvol", "airframe", "flight-model/planform.js loft + signedVolume", "divergence-theorem volume of the lofted body (M 3 to 300, k 1 or 2): as exported, fixed winding, exact",
		"record raw, good, exact (m^3)", {Num("M", "-", 70), Num("k", "-", 1), TableParam()}, &BodyVol);

	// Game (game/game.js)
	Def(R, "glidepar", "game", "game/game.js glidePar", "dead-stick par: h (L/D)max", "m", {Num("h", "m", 1000), Num("CD0", "-", 0.025), Num("AR", "-", 8), Num("e", "-", 0.8)}, &GlidePar);
	Def(R, "glidegame", "game", "game/game.js glideRun", "dead-stick glide from 1000 m at a held CL (0.2 to 1.4)", "record score (m), t (s), CL", {Num("CL", "-"), Num("dt", "s", 0.1)}, &GlideGame);
	Def(R, "pitchtrim", "game", "game/game.js pitchTrim", "trim alpha and elevator against the 0.5 rad/s^2 disturbance", "record alpha, delta (rad)", {}, &PitchTrim);
	Def(R, "pitchgame", "game", "game/game.js pitchRun", "hold-the-nose score of a PID with trim preset", "record J, ise, ctl (deg^2 s)",
		{Num("Kp", "-"), Num("Ki", "-"), Num("Kd", "-"), Num("trim", "rad", 0)}, &PitchGame);
	Def(R, "pitchpar", "game", "game/game.js pitchPar", "hold-the-nose par: finite-horizon discrete LQR value", "deg^2 s", {}, &PitchPar);
	Def(R, "radargame", "game", "game/game.js radarRun", "first-track range for a radar design (power, gain, scan)", "record range (nmi), cost, n, ok, scan",
		{Num("kw", "kW"), Num("gdb", "dB"), Num("T", "s")}, &RadarGame);
	Def(R, "radarpar", "game", "game/game.js radarPar", "first-track par by exhaustive search of the admissible designs", "record range (nmi), kw, gdb, T, n", {}, &RadarPar);
	Def(R, "parstars", "game", "game/game.js stars", "stars for a score against par (higher = 1 if bigger is better)", "0 to 3", {Num("score", "-"), Num("par", "-"), Num("higher", "-", 1)}, &Stars);
}
} // namespace mcalc
