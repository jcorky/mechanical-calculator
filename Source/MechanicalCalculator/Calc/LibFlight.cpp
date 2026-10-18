// Flight: the equation sheet (flight-model/equations.js), the NASA Glenn atmosphere, flight model and cloud
// chain (flight-model/model.js), the airlifter planform bench (flight-model/planform.js) and the par values of
// the game (game/game.js). Every routine follows its module step for step, so the numbers match to rounding.
#include "CalcLibrary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace mcalc
{
namespace flight_detail
{
constexpr double FlG0 = 9.80665;
constexpr double FlGamma = 1.4;
constexpr double FlCpAir = 1005.0;
constexpr double FlNan = std::numeric_limits<double>::quiet_NaN();
constexpr double FlInf = std::numeric_limits<double>::infinity();

double FlRad(double D) { return D * cst::Pi / 180.0; }
double JsRound(double X) { return std::floor(X + 0.5); } // Math.round
double Clamp(double X, double A, double B) { return std::min(B, std::max(A, X)); }

template <typename F>
double Trapz(F Fn, double A, double B, int N)
{
	const double H = (B - A) / N;
	double S = (Fn(A) + Fn(B)) / 2;
	for (int i = 1; i < N; ++i)
	{
		S += Fn(A + i * H);
	}
	return S * H;
}

// Gaussian elimination with partial pivoting on an n x (n + 1) augmented matrix, as equations.js's solve().
bool GaussSolve(std::vector<double> M, int N, std::vector<double>* X)
{
	const int W = N + 1;
	auto At = [&M, W](int R, int C) -> double& { return M[static_cast<size_t>(R) * static_cast<size_t>(W) + static_cast<size_t>(C)]; };
	for (int c = 0; c < N; ++c)
	{
		int P = c;
		for (int r = c + 1; r < N; ++r)
		{
			if (std::fabs(At(r, c)) > std::fabs(At(P, c)))
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
		if (At(c, c) == 0.0)
		{
			return false;
		}
		for (int r = c + 1; r < N; ++r)
		{
			const double F = At(r, c) / At(c, c);
			if (F != 0.0)
			{
				for (int j = c; j <= N; ++j)
				{
					At(r, j) -= F * At(c, j);
				}
			}
		}
	}
	X->assign(static_cast<size_t>(N), 0.0);
	for (int r = N - 1; r >= 0; --r)
	{
		double S = At(r, N);
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
double CamberSlope(double M, double P, double X)
{
	if (P <= 0 || M == 0)
	{
		return 0.0;
	}
	return X < P ? 2 * M / (P * P) * (P - X) : 2 * M / ((1 - P) * (1 - P)) * (P - X);
}

const char* CheckCamber(double M, double P)
{
	if (M < 0 || M > 0.2)
	{
		return "max camber m is a fraction of chord, 0 to 0.2";
	}
	if (M != 0 && (P <= 0 || P >= 1))
	{
		return "camber position p must lie strictly between 0 and 1";
	}
	return nullptr;
}

Result ThinAirfoil(const Call& C)
{
	const double M = C.N(0), P = C.N(1), Alpha = C.N(2);
	if (const char* E = CheckCamber(M, P))
	{
		return Result::Fail(E);
	}
	const int N = 20000;
	auto X = [](double T) { return (1 - std::cos(T)) / 2; };
	const double AL0 = -Trapz([&](double T) { return CamberSlope(M, P, X(T)) * (std::cos(T) - 1); }, 0, cst::Pi, N) / cst::Pi;
	const double A1 = 2 / cst::Pi * Trapz([&](double T) { return CamberSlope(M, P, X(T)) * std::cos(T); }, 0, cst::Pi, N);
	const double A2 = 2 / cst::Pi * Trapz([&](double T) { return CamberSlope(M, P, X(T)) * std::cos(2 * T); }, 0, cst::Pi, N);
	return Result::Record({{"aL0", AL0}, {"cl", 2 * cst::Pi * (Alpha - AL0)}, {"cm4", cst::Pi / 4 * (A2 - A1)}, {"A1", A1}, {"A2", A2}});
}

Result VortexPanel(const Call& C)
{
	const double M = C.N(0), P = C.N(1), Alpha = C.N(2);
	const int N = C.Int(3);
	if (const char* E = CheckCamber(M, P))
	{
		return Result::Fail(E);
	}
	if (N < 2 || N > 400)
	{
		return Result::Fail("panels N must be 2 to 400");
	}
	std::vector<double> Xe(static_cast<size_t>(N + 1)), Xv, Xc;
	for (int i = 0; i <= N; ++i)
	{
		Xe[static_cast<size_t>(i)] = (1 - std::cos(cst::Pi * i / N)) / 2;
	}
	for (int i = 0; i < N; ++i)
	{
		const double D = Xe[static_cast<size_t>(i + 1)] - Xe[static_cast<size_t>(i)];
		Xv.push_back(Xe[static_cast<size_t>(i)] + 0.25 * D);
		Xc.push_back(Xe[static_cast<size_t>(i)] + 0.75 * D);
	}
	std::vector<double> A(static_cast<size_t>(N) * static_cast<size_t>(N + 1));
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			A[static_cast<size_t>(i) * static_cast<size_t>(N + 1) + static_cast<size_t>(j)] = 1 / (2 * cst::Pi * (Xc[static_cast<size_t>(i)] - Xv[static_cast<size_t>(j)]));
		}
		A[static_cast<size_t>(i) * static_cast<size_t>(N + 1) + static_cast<size_t>(N)] = Alpha - CamberSlope(M, P, Xc[static_cast<size_t>(i)]);
	}
	std::vector<double> G;
	if (!GaussSolve(A, N, &G))
	{
		return Result::Fail("singular vortex system");
	}
	double Cl = 0, Cm4 = 0;
	for (int j = 0; j < N; ++j)
	{
		Cl += 2 * G[static_cast<size_t>(j)];
		Cm4 -= 2 * G[static_cast<size_t>(j)] * (Xv[static_cast<size_t>(j)] - 0.25);
	}
	return Result::Record({{"cl", Cl}, {"cm4", Cm4}});
}

Result PrandtlGlauert(const Call& C)
{
	const double C0 = C.N(0), M = C.N(1);
	if (M < 0 || M >= 1)
	{
		return Result::Fail("Prandtl-Glauert needs subsonic Mach, 0 <= M < 1");
	}
	return Result::Of(C0 / std::sqrt(1 - M * M));
}

Result KarmanTsien(const Call& C)
{
	const double Cp0 = C.N(0), M = C.N(1);
	if (M < 0 || M >= 1)
	{
		return Result::Fail("Karman-Tsien needs subsonic Mach, 0 <= M < 1");
	}
	const double B = std::sqrt(1 - M * M);
	return Result::Of(Cp0 / (B + M * M / (1 + B) * Cp0 / 2));
}

Result Ackeret(const Call& C)
{
	const double Alpha = C.N(0), M = C.N(1);
	if (M <= 1)
	{
		return Result::Fail("Ackeret needs supersonic Mach, M > 1");
	}
	return Result::Of(4 * Alpha / std::sqrt(M * M - 1));
}

Result Polar(const Call& C)
{
	const double CD0 = C.N(0), AR = C.N(1), E = C.N(2), CL = C.N(3);
	if (CD0 <= 0 || AR <= 0 || E <= 0)
	{
		return Result::Fail("CD0, AR and e must be positive");
	}
	const double K = 1 / (cst::Pi * E * AR);
	return Result::Record({{"k", K}, {"LDmax", 0.5 / std::sqrt(CD0 * K)}, {"CLstar", std::sqrt(CD0 / K)}, {"LD", CL / (CD0 + K * CL * CL)}});
}

// Prandtl's lifting line by Glauert's Fourier method, untwisted wing, chord(y) with y = 2y/b in [0, 1].
template <typename F>
bool LiftingLine(F Chord, double AR, int NTerms, double* E, double* Delta)
{
	const double A0 = 2 * cst::Pi;
	std::vector<double> Nk(static_cast<size_t>(NTerms)), Th(static_cast<size_t>(NTerms));
	for (int i = 0; i < NTerms; ++i)
	{
		Nk[static_cast<size_t>(i)] = 2 * i + 1;
		Th[static_cast<size_t>(i)] = (i + 0.5) * (cst::Pi / 2) / NTerms;
	}
	const double Cbar = Trapz(Chord, 0, 1, 4000), B = AR * Cbar;
	std::vector<double> M(static_cast<size_t>(NTerms) * static_cast<size_t>(NTerms + 1));
	for (int i = 0; i < NTerms; ++i)
	{
		const double T = Th[static_cast<size_t>(i)], Cc = Chord(std::cos(T));
		for (int k = 0; k < NTerms; ++k)
		{
			const double Nn = Nk[static_cast<size_t>(k)];
			M[static_cast<size_t>(i) * static_cast<size_t>(NTerms + 1) + static_cast<size_t>(k)] = std::sin(Nn * T) * (4 * B / (A0 * Cc) + Nn / std::sin(T));
		}
		M[static_cast<size_t>(i) * static_cast<size_t>(NTerms + 1) + static_cast<size_t>(NTerms)] = 1;
	}
	std::vector<double> A;
	if (!GaussSolve(M, NTerms, &A))
	{
		return false;
	}
	double D = 0;
	for (int i = 1; i < NTerms; ++i)
	{
		D += Nk[static_cast<size_t>(i)] * std::pow(A[static_cast<size_t>(i)] / A[0], 2);
	}
	*E = 1 / (1 + D);
	*Delta = D;
	return true;
}

Result SpanEff(const Call& C)
{
	const double AR = C.N(0), Lam = C.N(1);
	const int NT = C.Int(2);
	if (AR <= 0 || Lam <= 0 || Lam > 10)
	{
		return Result::Fail("needs AR > 0 and taper ratio 0 < lambda <= 10");
	}
	if (NT < 2 || NT > 100)
	{
		return Result::Fail("Fourier terms must be 2 to 100");
	}
	double E = 0, D = 0;
	if (!LiftingLine([Lam](double Y) { return 1 - (1 - Lam) * Y; }, AR, NT, &E, &D))
	{
		return Result::Fail("singular lifting-line system");
	}
	return Result::Record({{"e", E}, {"delta", D}});
}

Result SpanEffEll(const Call& C)
{
	const double AR = C.N(0);
	const int NT = C.Int(1);
	if (AR <= 0 || NT < 2 || NT > 100)
	{
		return Result::Fail("needs AR > 0 and 2 to 100 Fourier terms");
	}
	double E = 0, D = 0;
	if (!LiftingLine([](double Y) { return std::sqrt(std::max(0.0, 1 - Y * Y)); }, AR, NT, &E, &D))
	{
		return Result::Fail("singular lifting-line system");
	}
	return Result::Record({{"e", E}, {"delta", D}});
}

// Slender-body wave drag D/q by von Karman's Fourier series of S'(x).
template <typename F>
double WaveDq(F S, double L)
{
	const int NTerms = 40, Nq = 2001;
	const double H = 1e-6;
	std::vector<double> T(Nq), Sp(Nq);
	for (int i = 0; i < Nq; ++i)
	{
		T[static_cast<size_t>(i)] = cst::Pi * i / (Nq - 1);
		const double X = L * (1 - std::cos(T[static_cast<size_t>(i)])) / 2, A = std::max(0.0, X - H), B = std::min(L, X + H);
		Sp[static_cast<size_t>(i)] = (S(B) - S(A)) / (B - A);
	}
	double D = 0;
	for (int n = 1; n <= NTerms; ++n)
	{
		double Sum = 0;
		for (int i = 1; i < Nq; ++i)
		{
			const size_t I = static_cast<size_t>(i);
			Sum += (Sp[I] * std::sin(n * T[I]) + Sp[I - 1] * std::sin(n * T[I - 1])) / 2 * (T[I] - T[I - 1]);
		}
		const double An = 2 / cst::Pi * Sum / L;
		D += n * An * An;
	}
	return cst::Pi * L * L / 4 * D;
}

double SearsHaack(double RMax, double L, double X)
{
	const double U = std::min(std::max(X, 0.0), L);
	return cst::Pi * RMax * RMax * std::pow(4 * U * (L - U) / (L * L), 1.5);
}

double SearsHaackDq(double RMax, double L)
{
	const double V = 3 * std::pow(cst::Pi, 2) * RMax * RMax * L / 16;
	return 128 * V * V / (cst::Pi * std::pow(L, 4));
}

Result WaveDrag(const Call& C)
{
	const double R = C.N(0), L = C.N(1);
	if (R <= 0 || L <= 0)
	{
		return Result::Fail("rmax and length must be positive");
	}
	const double F = WaveDq([R, L](double X) { return SearsHaack(R, L, X); }, L), Cf = SearsHaackDq(R, L);
	return Result::Record({{"fourier", F}, {"closed", Cf}, {"ratio", F / Cf}});
}

// The equation sheet's area-rule card: Sears-Haack fuselage (rmax 0.05, l = 1) plus a wing's area, optionally indented.
Result AreaRule(const Call& C)
{
	const double Wing = C.N(0), Rule = C.N(1), R = 0.05, A0 = cst::Pi * R * R;
	if (Wing < 0)
	{
		return Result::Fail("wing area share must be >= 0");
	}
	auto SH = [R](double X) { return SearsHaack(R, 1, X); };
	auto WingA = [Wing, A0](double X) { return Wing * A0 * std::pow(std::max(0.0, 1 - std::pow((X - 0.55) / 0.12, 2)), 2); };
	const double DSH = WaveDq(SH, 1), DW = WaveDq([&](double X) { return SH(X) + WingA(X); }, 1);
	const double DT = WaveDq([&](double X) { return SH(X) + WingA(X) - Rule * WingA(X); }, 1);
	return Result::Record({{"fus", DSH}, {"wing", DW / DSH}, {"ruled", DT / DSH}, {"closed", SearsHaackDq(R, 1)}});
}

double AreaRatioOf(double M, double G)
{
	return 1 / M * std::pow((2 / (G + 1)) * (1 + (G - 1) / 2 * M * M), (G + 1) / (2 * (G - 1)));
}

Result AreaRatio(const Call& C)
{
	const double M = C.N(0), G = C.N(1);
	if (M <= 0 || G <= 1)
	{
		return Result::Fail("needs M > 0 and gamma > 1");
	}
	return Result::Of(AreaRatioOf(M, G));
}

Result MachFromArea(const Call& C)
{
	const double R = C.N(0), G = C.N(2);
	const bool Sup = C.N(1) != 0;
	if (R < 1 || G <= 1)
	{
		return Result::Fail("area ratio A/A* must be >= 1 (and gamma > 1)");
	}
	double Lo = Sup ? 1 : 1e-6, Hi = Sup ? 50 : 1;
	for (int i = 0; i < 200; ++i)
	{
		const double Mid = (Lo + Hi) / 2;
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
	const double TW = C.N(0), LD = C.N(1);
	if (LD <= 0)
	{
		return Result::Fail("L/D must be positive");
	}
	if (TW <= 1 / LD)
	{
		return Result::Fail("T/W is at or below 1/(L/D): the aircraft cannot hold level flight");
	}
	double Lo = 0, Hi = cst::Pi / 2;
	for (int i = 0; i < 200; ++i)
	{
		const double Mid = (Lo + Hi) / 2;
		if (std::sin(Mid) + std::cos(Mid) / LD < TW)
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
	const double Isp = C.N(0), M0 = C.N(1), Mf = C.N(2);
	if (M0 <= 0 || Mf <= 0)
	{
		return Result::Fail("masses must be positive");
	}
	return Result::Of(Isp * FlG0 * std::log(M0 / Mf));
}

Result Turboprop(const Call& C)
{
	const double Mdot = C.N(0), Pic = C.N(1), T02 = C.N(2), T04 = C.N(3), Etac = C.N(4), Etat = C.N(5), G = C.N(6), Cp = C.N(7);
	if (Pic <= 1 || Etac <= 0 || G <= 1)
	{
		return Result::Fail("needs pressure ratio > 1, compressor efficiency > 0 and gamma > 1");
	}
	const double Tau = std::pow(Pic, (G - 1) / G), Wt = Etat * T04 * (1 - 1 / Tau), Wc = T02 * (Tau - 1) / Etac, T03 = T02 + Wc, Q = Mdot * Cp * (T04 - T03);
	return Result::Record({{"tau", Tau}, {"P", Mdot * Cp * (Wt - Wc)}, {"turbine", Mdot * Cp * Wt}, {"compressor", Mdot * Cp * Wc}, {"q", Q}, {"eta", Mdot * Cp * (Wt - Wc) / Q}, {"ideal", 1 - 1 / Tau}});
}

Result NetThrust(const Call& C)
{
	const double Mdot = C.N(0), F = C.N(1), Ve = C.N(2), V0 = C.N(3), Pe = C.N(4), P0 = C.N(5), Ae = C.N(6);
	const double Gross = Mdot * (1 + F) * Ve + (Pe - P0) * Ae;
	return Result::Record({{"gross", Gross}, {"ram", Mdot * V0}, {"net", Mdot * (1 + F) * Ve + (Pe - P0) * Ae - Mdot * V0}});
}

/* ---------------------------------------------------------------- equation sheet: structures */
// Cantilever half-wing under elliptic running load w0 sqrt(1 - (y/s)^2), w0 from nW/2 per side; closed forms at y.
Result WingLoads(const Call& C)
{
	const double NW = C.N(0), S = C.N(1), Y = C.N(2);
	if (S <= 0 || Y < 0 || Y > S)
	{
		return Result::Fail("needs semi-span s > 0 and 0 <= y <= s");
	}
	const double W0 = 4 * (NW / 2) / (cst::Pi * S), U = Y / S, Root = std::sqrt(std::max(0.0, 1 - U * U));
	const double V = W0 * S * (std::acos(U) - U * Root) / 2;
	const double M = W0 * S * S * (Root * Root * Root / 3 - U * (std::acos(U) - U * Root) / 2);
	return Result::Record({{"w0", W0}, {"w", W0 * Root}, {"V", V}, {"M", M}});
}

Result Torsion(const Call& C)
{
	const double T = C.N(0), Rm = C.N(1), Tw = C.N(2), Gs = C.N(3);
	if (Rm <= 0 || Tw <= 0 || Tw >= 2 * Rm)
	{
		return Result::Fail("needs radius r > 0 and wall 0 < t < 2r");
	}
	const double Am = cst::Pi * Rm * Rm, Ro = Rm + Tw / 2, Ri = Rm - Tw / 2;
	return Result::Record({{"q", T / (2 * Am)},
		{"bredt", T / (2 * cst::Pi * Rm * Rm * Tw)},
		{"exact", T * Ro / (cst::Pi / 2 * (std::pow(Ro, 4) - std::pow(Ri, 4)))},
		{"open", 3 * T / (2 * cst::Pi * Rm * Tw * Tw)},
		{"twist", T * (2 * cst::Pi * Rm / Tw) / (4 * Am * Am * Gs)}});
}

// Smallest eigenvalue of the finite-difference column operator by inverse iteration (Thomas algorithm), times EI.
double FdBuckle(double E, double I, double L, bool Cantilever)
{
	const int N = 400, M = Cantilever ? N : N - 1;
	const double H = L / N;
	std::vector<double> D(static_cast<size_t>(M), 2 / (H * H)), Up(static_cast<size_t>(M - 1), -1 / (H * H)), Lo(static_cast<size_t>(M - 1), -1 / (H * H));
	if (Cantilever)
	{
		Up[0] = -2 / (H * H);
	}
	std::vector<double> X(static_cast<size_t>(M), 1.0), Y(static_cast<size_t>(M));
	double Lam = 0;
	for (int It = 0; It < 200; ++It)
	{
		std::vector<double> Cc = Up, Dd = D, R = X;
		for (int i = 1; i < M; ++i)
		{
			const size_t K = static_cast<size_t>(i);
			const double F = Lo[K - 1] / Dd[K - 1];
			Dd[K] -= F * Cc[K - 1];
			R[K] -= F * R[K - 1];
		}
		Y[static_cast<size_t>(M - 1)] = R[static_cast<size_t>(M - 1)] / Dd[static_cast<size_t>(M - 1)];
		for (int i = M - 2; i >= 0; --i)
		{
			const size_t K = static_cast<size_t>(i);
			Y[K] = (R[K] - Cc[K] * Y[K + 1]) / Dd[K];
		}
		double Num = 0, Den = 0;
		for (size_t i = 0; i < static_cast<size_t>(M); ++i)
		{
			Num += X[i] * Y[i];
			Den += Y[i] * Y[i];
		}
		Lam = Num / Den;
		const double Nrm = std::sqrt(Den);
		for (size_t i = 0; i < static_cast<size_t>(M); ++i)
		{
			X[i] = Y[i] / Nrm;
		}
	}
	return E * I * Lam;
}

Result Buckle(const Call& C)
{
	const double E = C.N(0), I = C.N(1), L = C.N(2), K = C.N(3);
	if (E <= 0 || I <= 0 || L <= 0 || K <= 0)
	{
		return Result::Fail("E, I, L and K must be positive");
	}
	const double P = std::pow(cst::Pi, 2) * E * I / std::pow(K * L, 2);
	const double Fd = K == 1 ? FdBuckle(E, I, L, false) : (K == 2 ? FdBuckle(E, I, L, true) : FlNan);
	return Result::Record({{"P", P}, {"fd", Fd}});
}

/* ---------------------------------------------------------------- equation sheet: flight mechanics */
constexpr double KSheet = 1 / (cst::Pi * 0.8 * 8); // the sheet's aircraft: AR 8, e 0.8

Result BestClimb(const Call& C)
{
	const double W = C.N(0), S = C.N(1), Rho = C.N(2), Pav = C.N(3), CD0 = C.N(4), K = C.N(5);
	if (W <= 0 || S <= 0 || Rho <= 0 || K <= 0 || CD0 <= 0)
	{
		return Result::Fail("W, S, rho, CD0 and k must be positive");
	}
	int B = 0;
	double Best = 0;
	for (int i = 0; i < 601; ++i)
	{
		const double V = 15 + i * 0.1, CL = 2 * W / (Rho * V * V * S);
		const double Roc = (Pav - 0.5 * Rho * std::pow(V, 3) * S * (CD0 + K * CL * CL)) / W;
		if (i == 0 || Roc > Best)
		{
			Best = Roc;
			B = i;
		}
	}
	return Result::Record({{"Vbest", 15 + B * 0.1}, {"rocmax", Best}, {"Vmp", std::sqrt(2 * W / (Rho * S * std::sqrt(3 * CD0 / K)))}});
}

Result FnTurn(const Call& C)
{
	const double V = C.N(0), Phi = C.N(1);
	if (V <= 0 || std::fabs(Phi) >= cst::Pi / 2)
	{
		return Result::Fail("needs V > 0 and |bank| < pi/2");
	}
	const double N = 1 / std::cos(Phi), S = std::sqrt(N * N - 1);
	return Result::Record({{"n", N}, {"omega", FlG0 * S / V}, {"R", S > 0 ? V * V / (FlG0 * S) : FlInf}});
}

Result CornerSpeed(const Call& C)
{
	const double Nmax = C.N(0), W = C.N(1), Rho = C.N(2), S = C.N(3), CLmax = C.N(4);
	if (Nmax <= 0 || W <= 0 || Rho <= 0 || S <= 0 || CLmax <= 0)
	{
		return Result::Fail("all inputs must be positive");
	}
	return Result::Of(std::sqrt(2 * Nmax * W / (Rho * S * CLmax)));
}

Result FnGlide(const Call& C)
{
	const double H = C.N(0), W = C.N(1), S = C.N(2), Rho = C.N(3), CD0 = C.N(4), K = C.N(5), CL = C.N(6);
	if (CL <= 0 || W <= 0 || S <= 0 || Rho <= 0)
	{
		return Result::Fail("CL, W, S and rho must be positive");
	}
	const double Gam = std::atan((CD0 + K * CL * CL) / CL), V = std::sqrt(2 * W * std::cos(Gam) / (Rho * S * CL));
	return Result::Record({{"gam", Gam}, {"V", V}, {"ratio", 1 / std::tan(Gam)}, {"dist", H / std::tan(Gam)}, {"sink", V * std::sin(Gam)}});
}

/* ---------------------------------------------------------------- equation sheet: guidance and control */
struct SpModel
{
	double Za, Zd, Ma, Mq, Md;
};

SpModel SpFrom(const Call& C, size_t First)
{
	return {C.N(First), C.N(First + 1), C.N(First + 2), C.N(First + 3), C.N(First + 4)};
}

Result ShortPeriod(const Call& C)
{
	const SpModel P = SpFrom(C, 0);
	const double A1 = -(P.Za + P.Mq), A0 = P.Za * P.Mq - P.Ma;
	if (A0 <= 0)
	{
		return Result::Fail("Za Mq - Ma <= 0: no short-period oscillator (statically unstable)");
	}
	const double Wn = std::sqrt(A0), Disc = A1 * A1 - 4 * A0;
	const double Td = Disc < 0 ? 2 * cst::Pi / (std::sqrt(-Disc) / 2) : FlNan;
	return Result::Record({{"wn", Wn}, {"zeta", A1 / (2 * Wn)}, {"a1", A1}, {"a0", A0}, {"b1", P.Md}, {"b0", P.Ma * P.Zd - P.Md * P.Za}, {"Td", Td}});
}

// theta/delta at s = j w from the transfer function, and from (jwI - A)^-1 B as a real 6 x 6 system.
Result PitchTf(const Call& C)
{
	const double W = C.N(0);
	const SpModel P = SpFrom(C, 1);
	const double A1 = -(P.Za + P.Mq), A0 = P.Za * P.Mq - P.Ma, N1 = P.Md, N0 = P.Ma * P.Zd - P.Md * P.Za;
	const double S0 = 0, S1 = W;
	const double NumR = N1 * S0 + N0, NumI = N1 * S1;
	const double S2R = S0 * S0 - S1 * S1, S2I = S0 * S1 + S1 * S0;
	const double QR = S2R + A1 * S0 + A0, QI = S2I + A1 * S1;
	const double DenR = S0 * QR - S1 * QI, DenI = S0 * QI + S1 * QR;
	const double Dd = DenR * DenR + DenI * DenI;
	if (Dd == 0)
	{
		return Result::Fail("w = 0 is the free-attitude pole");
	}
	const double Re = (NumR * DenR + NumI * DenI) / Dd, Im = (NumI * DenR - NumR * DenI) / Dd;
	const double A[3][3] = {{P.Za, 1, 0}, {P.Ma, P.Mq, 0}, {0, 1, 0}};
	const double B[3] = {P.Zd, P.Md, 0};
	std::vector<double> M(6 * 7, 0.0);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			M[static_cast<size_t>(i * 7 + j)] = -A[i][j];
			M[static_cast<size_t>(i * 7 + 3 + j)] = i == j ? -W : 0;
			M[static_cast<size_t>((i + 3) * 7 + j)] = i == j ? W : 0;
			M[static_cast<size_t>((i + 3) * 7 + 3 + j)] = -A[i][j];
		}
		M[static_cast<size_t>(i * 7 + 6)] = B[i];
	}
	std::vector<double> X;
	if (!GaussSolve(M, 6, &X))
	{
		return Result::Fail("singular state-space system");
	}
	const double Mag = std::hypot(Re, Im);
	return Result::Record({{"re", Re}, {"im", Im}, {"mag", Mag}, {"db", 20 * std::log10(Mag)}, {"ph", std::atan2(Im, Re)}, {"ssre", X[2]}, {"ssim", X[5]}});
}

Result POnlyError(const Call& C)
{
	const double Kp = C.N(0), Dist = C.N(1);
	const SpModel P = SpFrom(C, 2);
	if (Kp == 0 || P.Za == 0)
	{
		return Result::Fail("needs Kp != 0 and Za != 0");
	}
	return Result::Of(Dist / (Kp * (P.Md - P.Ma * P.Zd / P.Za)));
}

Result PidSim(const Call& C)
{
	const double Kp = C.N(0), Ki = C.N(1), Kd = C.N(2), ThetaC = C.N(3), Dist = C.N(4), T = 20, Dt = 1e-3;
	const SpModel P{-2.0, -0.2, -12.0, -3.0, -15.0};
	double A = 0, Q = 0, Th = 0, Z = 0, Peak = 0;
	for (int i = 0; i <= T / Dt; ++i)
	{
		const double E = ThetaC - Th;
		Peak = std::max(Peak, Th);
		const double D = -(Kp * E + Ki * Z) + Kd * Q, Da = P.Za * A + Q + P.Zd * D, Dq = P.Ma * A + P.Mq * Q + P.Md * D + Dist;
		A += Dt * Da;
		Q += Dt * Dq;
		Th += Dt * Q;
		Z += Dt * E;
		if (!std::isfinite(Th) || std::fabs(Th) > 1e3)
		{
			return Result::Record({{"err", FlNan}, {"overshoot", FlNan}, {"diverged", 1}});
		}
	}
	return Result::Record({{"err", ThetaC - Th}, {"overshoot", ThetaC != 0 ? std::max(0.0, (Peak - ThetaC) / ThetaC) : 0.0}, {"diverged", 0}});
}

/* ---------------------------------------------------------------- model.js: atmosphere and aircraft */
struct Atm
{
	double T, Tc, P, Rho, A, Sigma;
};

// NASA Glenn Earth Atmosphere Model (metric), every temperature shifted by dT.
Atm Atmosphere(double HIn, double DT)
{
	const double H = std::max(0.0, HIn);
	double Tc, Pk;
	if (H < 11000)
	{
		Tc = 15.04 - 0.00649 * H;
		Pk = 101.29 * std::pow((Tc + 273.1) / 288.08, 5.256);
	}
	else if (H < 25000)
	{
		Tc = -56.46;
		Pk = 22.65 * std::exp(1.73 - 0.000157 * H);
	}
	else
	{
		Tc = -131.21 + 0.00299 * H;
		Pk = 2.488 * std::pow((Tc + 273.1) / 216.6, -11.388);
	}
	Tc += DT;
	const double T = Tc + 273.1, Rho = Pk / (0.2869 * T);
	return {T, Tc, Pk * 1000, Rho, std::sqrt(1.4 * 286 * T), Rho / 1.225};
}

struct Aircraft
{
	double M, S, AR, E, CD0, Alpha0Deg, AlphaStallDeg, PlateCL, Mcrit, CDwave, CYb;
	bool Prop;
	double Power, Tstatic, EtaProp, Tmax, KM;
	double GroundAlphaDeg, RotateMaxDeg, Vrot;
};

// The three representative sizings of model.js: 1 general-aviation single, 2 jet transport, 3 supersonic delta.
const Aircraft* AircraftOf(int Index)
{
	static const Aircraft Table[3] = {
		{1050, 16.2, 7.4, 0.75, 0.028, -2.5, 15, 1.0, 0.7, 0.05, -0.6, true, 120e3, 3200, 0.8, 0, 0, 1, 12, 28},
		{300000, 511, 7.0, 0.8, 0.019, -3, 14, 1.0, 0.82, 0.06, -0.8, false, 0, 0, 0, 4 * 250e3, -0.35, 0, 11, 82},
		{150000, 358, 1.6, 0.85, 0.012, 0, 22, 1.0, 0.92, 0.022, -0.5, false, 0, 0, 0, 4 * 170e3, 0.7, 0, 14, 100},
	};
	return Index >= 1 && Index <= 3 ? &Table[Index - 1] : nullptr;
}

double SectionSlope(double M)
{
	auto Sub = [](double X) { return (2 * cst::Pi) / std::sqrt(1 - X * X); };
	auto Sup = [](double X) { return 4 / std::sqrt(X * X - 1); };
	if (M <= 0.85)
	{
		return Sub(M);
	}
	if (M >= 1.15)
	{
		return Sup(M);
	}
	const double T = (M - 0.85) / 0.3;
	return (1 - T) * Sub(0.85) + T * Sup(1.15);
}

double WingSlope(double A0, double AR)
{
	const double K = A0 / (cst::Pi * AR);
	return A0 / (std::sqrt(1 + K * K) + K);
}

double StallBlend(double Alpha, const Aircraft& Ac)
{
	const double Width = FlRad(2.0);
	return 1 / (1 + std::exp(-(std::fabs(Alpha) - FlRad(Ac.AlphaStallDeg)) / Width));
}

struct Aero
{
	Atm Air;
	double M, Q, CL, CLa, CLlin, CD0, CDi, CD, K, L, D;
	bool Stalled;
};

Aero AeroAt(double Alpha, double Beta, double V, double H, const Aircraft& Ac)
{
	Aero R;
	R.Air = Atmosphere(H, 0);
	R.M = V / R.Air.A;
	R.Q = 0.5 * R.Air.Rho * V * V;
	R.CLa = WingSlope(SectionSlope(R.M), Ac.AR);
	R.CLlin = R.CLa * (Alpha - FlRad(Ac.Alpha0Deg));
	const double Plate = 2 * std::sin(Alpha) * std::cos(Alpha) * Ac.PlateCL, S = StallBlend(Alpha, Ac);
	R.CL = (1 - S) * R.CLlin + S * Plate;
	R.Stalled = S > 0.5;
	double Wave = 0;
	if (R.M > Ac.Mcrit && R.M <= 1.05)
	{
		const double U = (R.M - Ac.Mcrit) / (1.05 - Ac.Mcrit);
		Wave = Ac.CDwave * U * U * U;
	}
	else if (R.M > 1.05)
	{
		Wave = Ac.CDwave * std::sqrt((1.05 * 1.05 - 1) / (R.M * R.M - 1));
	}
	R.CD0 = Ac.CD0 + Wave;
	R.K = 1 / (cst::Pi * Ac.AR * Ac.E);
	if (R.M > 1)
	{
		R.K = std::max(R.K, std::sqrt(R.M * R.M - 1) / 4);
	}
	R.CDi = (1 - S) * R.K * R.CL * R.CL + S * 2 * std::pow(std::sin(Alpha), 2) * Ac.PlateCL;
	R.CD = R.CD0 + R.CDi + 0.5 * Beta * Beta;
	R.L = R.Q * Ac.S * R.CL;
	R.D = R.Q * Ac.S * R.CD;
	return R;
}

double ThrustOf(double Throttle, double V, const Atm& Air, double M, const Aircraft& Ac)
{
	if (Ac.Prop)
	{
		const double P = Ac.Power * Throttle * Air.Sigma, Ts = Ac.Tstatic * Throttle * Air.Sigma;
		return std::min(Ts, (Ac.EtaProp * P) / std::max(V, 1.0));
	}
	return Ac.Tmax * Throttle * std::pow(Air.Sigma, 0.7) * std::max(0.2, 1 + Ac.KM * M);
}

double TrimAlpha(const Aircraft& Ac, double V, double H, double Gam, double Phi, double T)
{
	const double W = Ac.M * FlG0, Need = (W * std::cos(Gam)) / std::cos(Phi);
	double Lo = FlRad(Ac.Alpha0Deg) - FlRad(4), Hi = FlRad(Ac.AlphaStallDeg);
	for (int i = 0; i < 60; ++i)
	{
		const double Mid = 0.5 * (Lo + Hi);
		if (AeroAt(Mid, 0, V, H, Ac).L + T * std::sin(Mid) < Need)
		{
			Lo = Mid;
		}
		else
		{
			Hi = Mid;
		}
	}
	return 0.5 * (Lo + Hi);
}

const Aircraft* ArgAircraft(const Call& C, size_t I)
{
	const double X = C.N(I);
	return X == std::floor(X) ? AircraftOf(C.Int(I)) : nullptr;
}

const char* const BadAircraft = "aircraft: 1 general aviation, 2 jet transport, 3 supersonic delta";

Result FnAtmos(const Call& C)
{
	const Atm A = Atmosphere(C.N(0), C.N(1));
	return Result::Record({{"T", A.T}, {"Tc", A.Tc}, {"p", A.P}, {"rho", A.Rho}, {"a", A.A}, {"sigma", A.Sigma}, {"icing", A.Tc <= 2 && A.Tc >= -20 ? 1.0 : 0.0}});
}

Result LiftSlope(const Call& C)
{
	const double M = C.N(0), AR = C.N(1);
	if (M < 0 || AR < 0)
	{
		return Result::Fail("needs M >= 0 and AR >= 0 (0 = the section)");
	}
	const double A0 = SectionSlope(M);
	return Result::Record({{"a0", A0}, {"a", AR > 0 ? WingSlope(A0, AR) : A0}});
}

Result LiftCoef(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 2);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const double Alpha = C.N(0), M = C.N(1);
	if (M < 0)
	{
		return Result::Fail("Mach must be >= 0");
	}
	const double CLa = WingSlope(SectionSlope(M), Ac->AR), Lin = CLa * (Alpha - FlRad(Ac->Alpha0Deg));
	const double Plate = 2 * std::sin(Alpha) * std::cos(Alpha) * Ac->PlateCL, S = StallBlend(Alpha, *Ac);
	return Result::Record({{"CL", (1 - S) * Lin + S * Plate}, {"CLa", CLa}, {"linear", Lin}, {"stalled", S > 0.5 ? 1.0 : 0.0}});
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
	const Aero A = AeroAt(C.N(0), C.N(4), C.N(1), C.N(2), *Ac);
	return Result::Record({{"L", A.L}, {"D", A.D}, {"CL", A.CL}, {"CD", A.CD}, {"CD0", A.CD0}, {"CDi", A.CDi}, {"K", A.K}, {"M", A.M}, {"q", A.Q}});
}

Result FnThrust(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 3);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const double Throttle = C.N(0), V = C.N(1);
	const Atm A = Atmosphere(C.N(2), 0);
	return Result::Of(ThrustOf(Throttle, V, A, V / A.A, *Ac));
}

Result FnTrim(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 2);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const double V = C.N(0), H = C.N(1), Gam = C.N(3), Phi = C.N(4);
	if (V <= 0 || std::fabs(Phi) >= cst::Pi / 2)
	{
		return Result::Fail("needs V > 0 and |bank| < pi/2");
	}
	const double W = Ac->M * FlG0;
	double T = 0, Alpha = 0;
	for (int i = 0; i < 30; ++i)
	{
		Alpha = TrimAlpha(*Ac, V, H, Gam, Phi, T);
		const Aero A = AeroAt(Alpha, 0, V, H, *Ac);
		T = std::max(0.0, (A.D + W * std::sin(Gam)) / std::cos(Alpha));
	}
	const Atm Air = Atmosphere(H, 0);
	const double M = V / Air.A;
	if (ThrustOf(1, V, Air, M, *Ac) < T)
	{
		return Result::Record({{"alpha", Alpha}, {"throttle", FlNan}, {"T", T}});
	}
	double Lo = 0, Hi = 1;
	for (int i = 0; i < 50; ++i)
	{
		const double Mid = 0.5 * (Lo + Hi);
		if (ThrustOf(Mid, V, Air, M, *Ac) < T)
		{
			Lo = Mid;
		}
		else
		{
			Hi = Mid;
		}
	}
	return Result::Record({{"alpha", Alpha}, {"throttle", 0.5 * (Lo + Hi)}, {"T", T}});
}

Result LevelDrag(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 2);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const double V = C.N(0), H = C.N(1);
	if (V <= 0)
	{
		return Result::Fail("airspeed must be positive");
	}
	const double A = TrimAlpha(*Ac, V, H, 0, 0, 0);
	const Aero R = AeroAt(A, 0, V, H, *Ac);
	const bool Ok = R.L > 0.995 * Ac->M * FlG0 && A < FlRad(Ac->AlphaStallDeg) - FlRad(0.5);
	return Result::Record({{"D", R.D}, {"alpha", A}, {"CL", R.CL}, {"CD", R.CD}, {"M", R.M}, {"gam", std::atan(R.D / R.L)}, {"ok", Ok ? 1.0 : 0.0}});
}

// Ground roll from rest at full throttle (model.js step() on the ground, RK4 at 1/120 s), elevator once V > Vrot.
struct Roll
{
	double X, V, Alpha, Qa;
};

Roll RollRate(const Roll& S, double Throttle, double Elev, const Aircraft& Ac)
{
	const double W = Ac.M * FlG0, V = std::max(S.V, 0.5);
	const Aero A = AeroAt(S.Alpha, 0, V, 0, Ac);
	const double T = ThrustOf(Throttle, V, A.Air, A.M, Ac);
	const double Lz = A.L + T * std::sin(S.Alpha), N = std::max(0.0, W - Lz), Mu = 0.02;
	const double Fr = V > 0.6 || T > Mu * N ? Mu * N : 0;
	const double GA = FlRad(Ac.GroundAlphaDeg), RM = FlRad(Ac.RotateMaxDeg);
	const double Auth = Clamp((V - 0.5 * Ac.Vrot) / (0.5 * Ac.Vrot), 0, 1);
	const double Target = Clamp(GA + std::max(0.0, Elev) * Auth * (RM - GA), GA, RM);
	return {V, (T * std::cos(S.Alpha) - A.D - Fr) / Ac.M, S.Qa, 9 * (Target - S.Alpha) - 6 * S.Qa};
}

Roll RollAdd(const Roll& S, const Roll& D, double K) { return {S.X + K * D.X, S.V + K * D.V, S.Alpha + K * D.Alpha, S.Qa + K * D.Qa}; }

Result Takeoff(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 0);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const double Throttle = C.N(1), ElevCmd = C.N(2), Dt = 1.0 / 120;
	if (Throttle <= 0 || Throttle > 1)
	{
		return Result::Fail("throttle must be in (0, 1]");
	}
	Roll S{0, 0, FlRad(Ac->GroundAlphaDeg), 0};
	double T = 0;
	for (int i = 0; i < 120 * 120; ++i)
	{
		const double Elev = S.V > Ac->Vrot ? ElevCmd : 0;
		const Roll K1 = RollRate(S, Throttle, Elev, *Ac);
		const Roll K2 = RollRate(RollAdd(S, K1, Dt / 2), Throttle, Elev, *Ac);
		const Roll K3 = RollRate(RollAdd(S, K2, Dt / 2), Throttle, Elev, *Ac);
		const Roll K4 = RollRate(RollAdd(S, K3, Dt), Throttle, Elev, *Ac);
		Roll O;
		O.X = S.X + (Dt / 6) * (K1.X + 2 * K2.X + 2 * K3.X + K4.X);
		O.V = S.V + (Dt / 6) * (K1.V + 2 * K2.V + 2 * K3.V + K4.V);
		O.Alpha = S.Alpha + (Dt / 6) * (K1.Alpha + 2 * K2.Alpha + 2 * K3.Alpha + K4.Alpha);
		O.Qa = S.Qa + (Dt / 6) * (K1.Qa + 2 * K2.Qa + 2 * K3.Qa + K4.Qa);
		O.V = std::max(O.V, 0.0);
		T += Dt;
		S = O;
		const Aero A = AeroAt(S.Alpha, 0, std::max(S.V, 0.5), 0, *Ac);
		const double Th = ThrustOf(Throttle, S.V, A.Air, A.M, *Ac);
		if (A.L + Th * std::sin(S.Alpha) > Ac->M * FlG0 * 1.001)
		{
			return Result::Record({{"t", T}, {"dist", S.X}, {"V", S.V}, {"alpha", S.Alpha}});
		}
	}
	return Result::Fail("no liftoff within 120 s");
}

/* ---------------------------------------------------------------- model.js: the cloud chain */
constexpr double ClCp = 1004, ClRd = 287, ClLv = 2.501e6, ClEps = 0.622;

double SatVapour(double Tc) { return 611.2 * std::exp((17.67 * Tc) / (Tc + 243.5)); }
double SatMixing(double Tc, double P)
{
	const double E = SatVapour(Tc);
	return (ClEps * E) / std::max(P - E, 1.0);
}
double MoistLapseOf(double T, double P)
{
	const double Rs = SatMixing(T - 273.15, P);
	return (FlG0 * (1 + (ClLv * Rs) / (ClRd * T))) / (ClCp + (ClLv * ClLv * Rs * ClEps) / (ClRd * T * T));
}
double DewPointOf(double R, double P)
{
	const double E = (R * P) / (ClEps + R), L = std::log(E / 611.2);
	return (243.5 * L) / (17.67 - L);
}

struct Sounding
{
	double Lcl, Lfc, El, Cape, Cin;
	bool Convective;
	std::vector<double> CapeCum;   // per 25 m level, as prof[].capeCum
};

// Lift a surface parcel through the environment: dry below cloud base, moist above (model.js sounding()).
Sounding RunSounding(double DT, double DewSpread, double Heating)
{
	const double Dz = 25;
	const Atm Env0 = Atmosphere(0, DT);
	double Tp = Env0.T + Heating;
	const double Td0 = Env0.Tc - DewSpread, R0 = SatMixing(Td0, Env0.P);
	bool HasLcl = false, HasLfc = false, HasEl = false;
	Sounding Out{FlNan, FlNan, FlNan, 0, 0, false, {}};
	double CapeCum = 0;
	for (double Z = 0; Z <= 16000; Z += Dz)
	{
		const Atm Env = Atmosphere(Z, DT);
		const bool Saturated = HasLcl;
		const double Tdp = Saturated ? Tp - 273.15 : DewPointOf(R0, Env.P);
		if (!Saturated && Tp - 273.15 <= Tdp)
		{
			HasLcl = true;
			Out.Lcl = Z;
		}
		const double B = (FlG0 * (Tp - Env.T)) / Env.T;
		if (HasLcl)
		{
			if (B > 0 && !HasLfc)
			{
				HasLfc = true;
				Out.Lfc = Z;
			}
			if (HasLfc && !HasEl)
			{
				if (B > 0)
				{
					Out.Cape += B * Dz;
					CapeCum += B * Dz;
				}
				else
				{
					HasEl = true;
					Out.El = Z;
				}
			}
			else if (!HasLfc)
			{
				Out.Cin += std::min(0.0, B) * Dz;
			}
		}
		else
		{
			Out.Cin += std::min(0.0, B) * Dz;
		}
		if (!HasLfc && Z >= 6000)
		{
			break;
		}
		Out.CapeCum.push_back(CapeCum);
		Tp -= (!HasLcl ? FlG0 / ClCp : MoistLapseOf(Tp, Env.P)) * Dz;
	}
	if (HasLfc && !HasEl)
	{
		Out.El = 16000;
	}
	if (!HasLfc)
	{
		Out.Cin = FlNan;
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
	return Result::Record({{"lcl", S.Lcl}, {"lfc", S.Lfc}, {"el", S.El}, {"cape", S.Cape}, {"cin", S.Cin}, {"conv", S.Convective ? 1.0 : 0.0}});
}

Result Updraft(const Call& C)
{
	const Sounding S = RunSounding(C.N(1), C.N(2), C.N(3));
	if (std::isnan(S.Lfc) || S.CapeCum.empty())
	{
		return Result::Of(0.0);
	}
	const double Idx = Clamp(JsRound(C.N(0) / 25), 0, static_cast<double>(S.CapeCum.size() - 1));
	return Result::Of(std::min(25.0, 0.4 * std::sqrt(2 * S.CapeCum[static_cast<size_t>(Idx)]) + 1.5));
}

// A sharp vertical gust w tilts the relative wind by w/V (model.js derivatives()); load factor from level trim.
Result GustLoad(const Call& C)
{
	const Aircraft* Ac = ArgAircraft(C, 3);
	if (!Ac)
	{
		return Result::Fail(BadAircraft);
	}
	const double Wg = C.N(0), V = C.N(1), H = C.N(2);
	if (V <= 0)
	{
		return Result::Fail("airspeed must be positive");
	}
	const double W = Ac->M * FlG0, A0 = TrimAlpha(*Ac, V, H, 0, 0, 0), Da = Wg / V;
	const double N0 = AeroAt(A0, 0, V, H, *Ac).L / W, N = AeroAt(A0 + Da, 0, V, H, *Ac).L / W;
	return Result::Record({{"dalpha", Da}, {"n", N}, {"dn", N - N0}, {"alpha0", A0}});
}

/* ---------------------------------------------------------------- planform.js: the airlifter */
struct Table
{
	int Rows = 0, Cols = 0;
	std::vector<double> V;
	double At(int R, int C) const { return V[static_cast<size_t>(R) * static_cast<size_t>(Cols) + static_cast<size_t>(C)]; }
};

Table DefaultTable()
{
	const double D2R = cst::Pi / 180, T55 = std::tan(55 * D2R), T35 = std::tan(35 * D2R);
	const double LE9 = -9 * T55, TE4 = -46 + 4 * T35, TE13 = TE4 + 9 * T55;
	Table T;
	T.Rows = 5;
	T.Cols = 6;
	T.V = {0, 0, -46, 6.0, 0.00, 0.62,
		4, -4 * T55, TE4, 5.6, 0.08, 0.61,
		9, LE9, TE4 + 5 * T55, 3.4, 0.35, 0.58,
		13, LE9 - 4 * T35, TE13, 1.7, 0.70, 0.56,
		31, LE9 - 22 * T35, -31.6, 0.34, 1.45, 0.54};
	return T;
}

// A table argument: 0 (the default) means the generator's five stations; otherwise n x k, k >= MinCols.
bool TableArg(const Call& C, size_t I, int MinCols, Table* Out, std::string* Error)
{
	const Matrix& M = C.M(I);
	if (M.IsScalar() && M.ScalarValue() == 0)
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
		if (i == 0 && Out->At(0, 0) != 0)
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

double Interp(const Table& T, int Col, double X)   // numpy.interp on column Col against y, clamped
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
	const double F = (X - T.At(i - 1, 0)) / (T.At(i, 0) - T.At(i - 1, 0));
	return T.At(i - 1, Col) + F * (T.At(i, Col) - T.At(i - 1, Col));
}

struct Stn
{
	double Xle, Xte, T, Zc, Up;
};

Stn StationAt(const Table& T, double Y)
{
	const double A = std::fabs(Y);
	return {Interp(T, 1, A), Interp(T, 2, A), Interp(T, 3, A), Interp(T, 4, A), Interp(T, 5, A)};
}

const double NacaCo[5] = {0.2969, -0.1260, -0.3516, 0.2843, -0.1036};
double Naca(double S) { return NacaCo[0] * std::sqrt(S) + NacaCo[1] * S + NacaCo[2] * S * S + NacaCo[3] * S * S * S + NacaCo[4] * S * S * S * S; }
double NacaMax()
{
	static const double M = [] {
		double Best = -FlInf;
		for (int i = 0; i <= 20000; ++i)
		{
			Best = std::max(Best, Naca(i == 20000 ? 1 : i * (1.0 / 20000)));
		}
		return Best;
	}();
	return M;
}
double Shape(double S) { return Naca(S) / NacaMax(); }
double ShapeInt() { return (NacaCo[0] * 2 / 3 + NacaCo[1] / 2 + NacaCo[2] / 3 + NacaCo[3] / 4 + NacaCo[4] / 5) / NacaMax(); }

double Surf(const Table& T, double X, double Y, bool Upper)
{
	const Stn St = StationAt(T, Y);
	const double S = std::min(1.0, std::max(0.0, (St.Xle - X) / (St.Xle - St.Xte))), H = Shape(S) * St.T;
	return Upper ? St.Zc + H * St.Up : St.Zc - H * (1 - St.Up);
}

double SegInt(double F0, double F1, double G0, double G1, double H)
{
	const double Df = F1 - F0, Dg = G1 - G0;
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
	double A = 0, Icc = 0, Icy = 0, Icx = 0, Ict = 0;
	for (int i = 0; i < N - 1; ++i)
	{
		const double H = T.At(i + 1, 0) - T.At(i, 0);
		A += H * (Ch(i) + Ch(i + 1)) / 2;
		Icc += SegInt(Ch(i), Ch(i + 1), Ch(i), Ch(i + 1), H);
		Icy += SegInt(Ch(i), Ch(i + 1), T.At(i, 0), T.At(i + 1, 0), H);
		Icx += SegInt(Ch(i), Ch(i + 1), T.At(i, 1), T.At(i + 1, 1), H);
		if (T.Cols >= 4)
		{
			Ict += SegInt(Ch(i), Ch(i + 1), T.At(i, 3), T.At(i + 1, 3), H);
		}
	}
	const double S = 2 * A, B = 2 * T.At(N - 1, 0), Mac = 2 * Icc / S;
	double YEq = FlNan;
	for (int i = 0; i < N - 1; ++i)
	{
		if ((Ch(i) - Mac) * (Ch(i + 1) - Mac) <= 0 && Ch(i) != Ch(i + 1))
		{
			YEq = T.At(i, 0) + (Ch(i) - Mac) / (Ch(i) - Ch(i + 1)) * (T.At(i + 1, 0) - T.At(i, 0));
			break;
		}
	}
	return Result::Record({{"S", S}, {"b", B}, {"AR", B * B / S}, {"MAC", Mac}, {"ymac", 2 * Icy / S}, {"xmac", 2 * Icx / S}, {"taper", Ch(N - 1) / Ch(0)},
		{"cr", Ch(0)}, {"ct", Ch(N - 1)}, {"yeq", YEq}, {"vol", T.Cols >= 4 ? 2 * ShapeInt() * Ict : FlNan}});
}

Result Sweep(const Call& C)
{
	Table T;
	std::string E;
	if (!TableArg(C, 1, 3, &T, &E))
	{
		return Result::Fail(E);
	}
	const double X = C.N(0);
	const int I = C.Int(0) - 1;
	if (X != std::floor(X) || I < 0 || I >= T.Rows - 1)
	{
		return Result::Fail("segment must be 1 to " + std::to_string(T.Rows - 1));
	}
	const double H = T.At(I + 1, 0) - T.At(I, 0);
	return Result::Record({{"le", std::atan(-(T.At(I + 1, 1) - T.At(I, 1)) / H)}, {"te", std::atan(-(T.At(I + 1, 2) - T.At(I, 2)) / H)},
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
	const double Nd = C.N(0);
	if (Nd != std::floor(Nd) || Nd < 1 || Nd > 1e6)
	{
		return Result::Fail("intervals n must be a whole number, 1 to 1e6");
	}
	const int N = C.Int(0);
	const double Half = T.At(T.Rows - 1, 0), H = Half / N;
	double A = 0, Cc = 0, Cy = 0, Cx = 0;
	for (int i = 0; i <= N; ++i)
	{
		const double Y = i == N ? Half : i * H, Xl = Interp(T, 1, Y), Xt = Interp(T, 2, Y), Ch = Xl - Xt, W = (i == 0 || i == N) ? H / 2 : H;
		A += W * Ch;
		Cc += W * Ch * Ch;
		Cy += W * Ch * Y;
		Cx += W * Ch * Xl;
	}
	const double S = 2 * A;
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
constexpr double IsaR = 287.05287;

void IsaOf(double H, double* T, double* P)
{
	const double P11 = 101325 * std::pow(216.65 / 288.15, FlG0 / (IsaR * 0.0065));
	if (H <= 11000)
	{
		*T = 288.15 - 0.0065 * H;
		*P = 101325 * std::pow(*T / 288.15, FlG0 / (IsaR * 0.0065));
	}
	else
	{
		*T = 216.65;
		*P = P11 * std::exp(-FlG0 * (H - 11000) / (IsaR * *T));
	}
}

Result FnIsa(const Call& C)
{
	const double H = C.N(0);
	if (H < 0 || H > 20000)
	{
		return Result::Fail("ISA here covers 0 to 20 000 m");
	}
	double T = 0, P = 0;
	IsaOf(H, &T, &P);
	auto Tf = [](double Z) { return 288.15 - 0.0065 * std::min(Z, 11000.0); };
	auto F = [&Tf](double Z) { return -FlG0 / (IsaR * Tf(Z)); };
	double Lp = std::log(101325.0), Z = 0;
	const int N = 2000;
	const double Ends[2] = {H > 11000 ? 11000 : H, H};
	for (int Leg = 0; Leg < (H > 11000 ? 2 : 1); ++Leg)
	{
		const double Dz = (Ends[Leg] - Z) / N;
		for (int i = 0; i < N; ++i)
		{
			Lp += Dz / 6 * (F(Z) + 4 * F(Z + Dz / 2) + F(Z + Dz));
			Z += Dz;
		}
	}
	const double Ph = std::exp(Lp);
	return Result::Record({{"T", T}, {"p", P}, {"rho", P / (IsaR * T)}, {"a", std::sqrt(FlGamma * IsaR * T)}, {"phyd", Ph}, {"rhohyd", Ph / (IsaR * Tf(H))}});
}

double DefaultArea()
{
	static const double S = [] {
		const Table T = DefaultTable();
		double A = 0;
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
	const double Mass = C.N(0), H = C.N(1), Mach = C.N(2), Thrust = C.N(3), S = C.N(4) == 0 ? DefaultArea() : C.N(4);
	if (H < 0 || H > 20000 || Mass <= 0 || S <= 0 || Mach <= 0)
	{
		return Result::Fail("needs mass > 0, 0 <= h <= 20 000 m, Mach > 0 and S > 0");
	}
	double T = 0, P = 0;
	IsaOf(H, &T, &P);
	const double Rho = P / (IsaR * T), A = std::sqrt(FlGamma * IsaR * T), V = Mach * A, W = Mass * FlG0, Q1 = 0.5 * Rho * V * V, Q2 = 0.5 * FlGamma * P * Mach * Mach;
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
	const double L = C.N(0), Wd = C.N(1);
	if (L <= 0 || Wd <= 0)
	{
		return Result::Fail("hold length and width must be positive");
	}
	double Xa = FlInf, Xb = -FlInf;
	for (int i = 0; i < T.Rows; ++i)
	{
		Xa = std::min(Xa, T.At(i, 2));
		Xb = std::max(Xb, T.At(i, 1));
	}
	if (Xb - L <= Xa)
	{
		return Result::Fail("the hold is longer than the centre body allows");
	}
	const double Y = Wd / 2;
	auto Hx = [&](double X0) {
		return std::min(Surf(T, X0, Y, true), Surf(T, X0 + L, Y, true)) - std::max(Surf(T, X0, Y, false), Surf(T, X0 + L, Y, false));
	};
	double A = Xa, B = Xb - L;
	const double G = (std::sqrt(5.0) - 1) / 2;
	for (int i = 0; i < 200; ++i)
	{
		const double M1 = B - G * (B - A), M2 = A + G * (B - A);
		if (Hx(M1) < Hx(M2))
		{
			A = M1;
		}
		else
		{
			B = M2;
		}
	}
	const double X0 = (A + B) / 2;
	// the scan
	const double Dx = 0.05;
	const int Ny = 25, N = static_cast<int>(JsRound((Xb - Xa) / Dx));
	const double Step = (Xb - Xa) / N;
	const int Nw = static_cast<int>(JsRound(L / Step));
	std::vector<double> Xs(static_cast<size_t>(N + 1)), Hi(static_cast<size_t>(N + 1)), Lo(static_cast<size_t>(N + 1));
	bool Sides = true;
	for (int i = 0; i <= N; ++i)
	{
		const double X = i == N ? Xb : Xa + i * Step;
		double U = FlInf, Lw = -FlInf;
		for (int k = 0; k < Ny; ++k)
		{
			const double Yk = -Wd / 2 + Wd * k / (Ny - 1);
			U = std::min(U, Surf(T, X, Yk, true));
			Lw = std::max(Lw, Surf(T, X, Yk, false));
		}
		Xs[static_cast<size_t>(i)] = X;
		Hi[static_cast<size_t>(i)] = U;
		Lo[static_cast<size_t>(i)] = Lw;
		if (U > Lw + 0.01 && (std::fabs(Surf(T, X, Wd / 2, true) - U) > 1e-9 || std::fabs(Surf(T, X, Wd / 2, false) - Lw) > 1e-9))
		{
			Sides = false;
		}
	}
	double Best = -FlInf;
	int Bi = 0;
	if (Nw >= 1 && Nw <= N)
	{
		for (int S0 = 0; S0 + Nw <= N; ++S0)
		{
			double U = FlInf, Lw = -FlInf;
			for (int i = S0; i <= S0 + Nw; ++i)
			{
				U = std::min(U, Hi[static_cast<size_t>(i)]);
				Lw = std::max(Lw, Lo[static_cast<size_t>(i)]);
			}
			if (U - Lw > Best)
			{
				Best = U - Lw;
				Bi = S0;
			}
		}
	}
	return Result::Record({{"H", Hx(X0)}, {"x0", X0}, {"x1", X0 + L}, {"floor", std::max(Surf(T, X0, Y, false), Surf(T, X0 + L, Y, false))},
		{"ceil", std::min(Surf(T, X0, Y, true), Surf(T, X0 + L, Y, true))}, {"Hscan", Best}, {"xscan", Nw >= 1 && Nw <= N ? Xs[static_cast<size_t>(Bi)] : FlNan}, {"sides", Sides ? 1.0 : 0.0}});
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
	const double Md = C.N(0), Kd = C.N(1);
	if (Md != std::floor(Md) || Md < 3 || Md > 600 || Kd != std::floor(Kd) || Kd < 1 || Kd > 4)
	{
		return Result::Fail("needs a whole number of chordwise points M, 3 to 600, and spanwise multiplier k, 1 to 4");
	}
	const int M = C.Int(0), K = C.Int(1);
	static const int SegSamples[4] = {10, 12, 10, 30};
	std::vector<double> Half;
	for (int i = 0; i < T.Rows - 1; ++i)
	{
		const int Ns = SegSamples[std::min(i, 3)] * K;
		const double A = T.At(i, 0), B = T.At(i + 1, 0), St = (B - A) / Ns;
		for (int j = 0; j < Ns; ++j)
		{
			Half.push_back(A + j * St);
		}
	}
	Half.push_back(T.At(T.Rows - 1, 0));
	std::vector<double> Ys;
	for (double V : Half)
	{
		Ys.push_back(V == 0 ? 0.0 : -V);
		Ys.push_back(V);
	}
	std::sort(Ys.begin(), Ys.end());
	Ys.erase(std::unique(Ys.begin(), Ys.end()), Ys.end());
	std::vector<double> S(static_cast<size_t>(M));
	for (int i = 0; i < M; ++i)
	{
		const double V = i == M - 1 ? cst::Pi : i * (cst::Pi / (M - 1));
		S[static_cast<size_t>(i)] = (1 - std::cos(V)) / 2;
	}
	std::vector<double> Rs(S.rbegin(), S.rend());
	Rs.insert(Rs.end(), S.begin() + 1, S.end());
	const int R = static_cast<int>(Rs.size()), NY = static_cast<int>(Ys.size());
	std::vector<double> V;
	V.reserve(static_cast<size_t>(3 * (NY * R + 2)));
	for (double Y : Ys)
	{
		const Stn St = StationAt(T, Y);
		for (int j = 0; j < R; ++j)
		{
			const double Sv = Rs[static_cast<size_t>(j)], H = Shape(Sv) * St.T;
			V.push_back(St.Xle - Sv * (St.Xle - St.Xte));
			V.push_back(Y);
			V.push_back(j < M ? St.Zc + H * St.Up : St.Zc - H * (1 - St.Up));
		}
	}
	auto Triple = [&V](int A, int B, int Cc) {
		const size_t a = static_cast<size_t>(3 * A), b = static_cast<size_t>(3 * B), c = static_cast<size_t>(3 * Cc);
		return V[a] * (V[b + 1] * V[c + 2] - V[b + 2] * V[c + 1]) - V[a + 1] * (V[b] * V[c + 2] - V[b + 2] * V[c]) + V[a + 2] * (V[b] * V[c + 1] - V[b + 1] * V[c]);
	};
	double Raw = 0, Good = 0;
	for (int i = 0; i < NY - 1; ++i)
	{
		for (int j = 0; j < R - 1; ++j)
		{
			const int A = i * R + j, B = A + 1, Cc = A + R, D = Cc + 1;
			Raw += Triple(A, Cc, B);
			Raw += Triple(B, Cc, D);
			Good += Triple(A, B, Cc);
			Good += Triple(B, D, Cc);
		}
	}
	for (int Cap = 0; Cap < 2; ++Cap)
	{
		const int Base = Cap == 0 ? 0 : (NY - 1) * R, Ci = static_cast<int>(V.size() / 3);
		double Cx = 0, Cy = 0, Cz = 0;
		for (int j = 0; j < R; ++j)
		{
			Cx += V[static_cast<size_t>(3 * (Base + j))];
			Cy += V[static_cast<size_t>(3 * (Base + j) + 1)];
			Cz += V[static_cast<size_t>(3 * (Base + j) + 2)];
		}
		V.push_back(Cx / R);
		V.push_back(Cy / R);
		V.push_back(Cz / R);
		for (int j = 0; j < R - 1; ++j)
		{
			const double Tv = Cap == 0 ? Triple(Base + j + 1, Base + j, Ci) : Triple(Ci, Base + j, Base + j + 1);
			Raw += Tv;
			Good += Tv;
		}
	}
	double Ict = 0;
	for (int i = 0; i < T.Rows - 1; ++i)
	{
		const double C0 = T.At(i, 1) - T.At(i, 2), C1 = T.At(i + 1, 1) - T.At(i + 1, 2);
		Ict += SegInt(C0, C1, T.At(i, 3), T.At(i + 1, 3), T.At(i + 1, 0) - T.At(i, 0));
	}
	return Result::Record({{"raw", Raw / 6}, {"good", Good / 6}, {"exact", 2 * ShapeInt() * Ict}});
}

/* ---------------------------------------------------------------- game.js: par is a theorem */
// Level 1: dead-stick glide on the sheet's polar from 1000 m.
double GlideK() { return 1 / (cst::Pi * 0.8 * 8); }

Result GlidePar(const Call& C)
{
	const double H = C.N(0), CD0 = C.N(1), AR = C.N(2), E = C.N(3);
	if (CD0 <= 0 || AR <= 0 || E <= 0)
	{
		return Result::Fail("CD0, AR and e must be positive");
	}
	const double K = 1 / (cst::Pi * E * AR);
	return Result::Of(H * (0.5 / std::sqrt(CD0 * K)));
}

Result GlideGame(const Call& C)
{
	const double Dt = C.N(1);
	if (!(Dt > 0) || Dt > 10)
	{
		return Result::Fail("time step must be in (0, 10] s");
	}
	const double W = 10000, S = 16, CD0 = 0.025, Rho = 1.225, K = GlideK();
	const double Cc = std::min(1.4, std::max(0.2, C.N(0))), CD = CD0 + K * Cc * Cc, Gam = std::atan(CD / Cc);
	const double V = std::sqrt(2 * W * std::cos(Gam) / (Rho * S * Cc));
	double X = 0, Z = 1000, T = 0;
	bool Done = false;
	for (int n = 0; !Done && n < 1000000; ++n)
	{
		double Dz = V * std::sin(Gam) * Dt, H = Dt;
		if (Dz >= Z)
		{
			H = Dt * Z / Dz;
			Dz = Z;
			Done = true;
		}
		X += Dz * Cc / CD;
		Z = Done ? 0 : Z - Dz;
		T += H;
	}
	return Result::Record({{"score", X}, {"t", T}, {"CL", Cc}});
}

// Level 2: hold 5 deg against a 0.5 rad/s^2 disturbance on the sheet's short-period model.
constexpr double PhThetaC = 5 * cst::Pi / 180, PhDist = 0.5, PhDt = 1e-3, PhRho = 0.01;
constexpr int PhSteps = 20000;
const SpModel PhSp{-2.0, -0.2, -12.0, -3.0, -15.0};

void PitchTrimOf(double* Alpha, double* Delta)
{
	const SpModel& P = PhSp;
	const double Det = P.Za * P.Md - P.Zd * P.Ma, D = -P.Za * PhDist / Det;
	*Delta = D;
	*Alpha = -P.Zd * D / P.Za;
}

Result PitchTrim(const Call& C)
{
	(void)C;
	double A = 0, D = 0;
	PitchTrimOf(&A, &D);
	return Result::Record({{"alpha", A}, {"delta", D}});
}

Result PitchGame(const Call& C)
{
	const double Kp = C.N(0), Ki = C.N(1), Kd = C.N(2), Ff = C.N(3);
	const SpModel& P = PhSp;
	double TrA = 0, TrD = 0;
	PitchTrimOf(&TrA, &TrD);
	double A = 0, Q = 0, Th = 0, Z = 0, Ise = 0, Ctl = 0;
	const double Deg = 180 / cst::Pi;
	for (int k = 0; k < PhSteps; ++k)
	{
		const double E = PhThetaC - Th, D = Ff - (Kp * E + Ki * Z) + Kd * Q, U = D - TrD;
		Ise += E * E * Deg * Deg * PhDt;
		Ctl += U * U * Deg * Deg * PhDt;
		const double Da = P.Za * A + Q + P.Zd * D, Dq = P.Ma * A + P.Mq * Q + P.Md * D + PhDist;
		A += PhDt * Da;
		Q += PhDt * Dq;
		Th += PhDt * Q;
		Z += PhDt * E;
		if (!std::isfinite(Th) || std::fabs(Th) > 1e3)
		{
			return Result::Record({{"J", FlInf}, {"ise", FlInf}, {"ctl", FlInf}});
		}
	}
	return Result::Record({{"J", Ise + PhRho * Ctl}, {"ise", Ise}, {"ctl", Ctl}});
}

// Par: finite-horizon discrete LQR (Riccati recursion) on the simulation's exact one-step map.
Result PitchPar(const Call& C)
{
	(void)C;
	const SpModel& P = PhSp;
	const double H = PhDt, Deg = 180 / cst::Pi;
	const double Ad[3][3] = {{1 + H * P.Za, H, 0}, {H * P.Ma, 1 + H * P.Mq, 0}, {H * H * P.Ma, H + H * H * P.Mq, 1}};
	const double Bd[3] = {H * P.Zd, H * P.Md, H * H * P.Md};
	const double Qw = Deg * Deg * PhDt, Rw = PhRho * Deg * Deg * PhDt;
	double TrA = 0, TrD = 0;
	PitchTrimOf(&TrA, &TrD);
	double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	for (int k = PhSteps - 1; k >= 0; --k)
	{
		double SB[3], Kk[3], Cm[3][3], NS[3][3];
		for (int r = 0; r < 3; ++r)
		{
			SB[r] = S[r][0] * Bd[0] + S[r][1] * Bd[1] + S[r][2] * Bd[2];
		}
		const double Den = Rw + Bd[0] * SB[0] + Bd[1] * SB[1] + Bd[2] * SB[2];
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
				double Sum = 0;
				for (int m = 0; m < 3; ++m)
				{
					Sum += Ad[m][i] * (S[m][0] * Cm[0][j] + S[m][1] * Cm[1][j] + S[m][2] * Cm[2][j]);
				}
				NS[i][j] = Sum + (i == 2 && j == 2 ? Qw : 0);
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
	const double X0[3] = {-TrA, 0, -PhThetaC};
	double V = 0;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			V += X0[i] * S[i][j] * X0[j];
		}
	}
	return Result::Of(V);
}

// Level 3: first track. Airport radar at 2.8 GHz against a 1 m^2 aircraft inbound at 120 kt from 60 nmi.
constexpr double RaNmi = 1852, RaKt = 1852.0 / 3600, RaPfa = 1e-6;
const double RaPower[5][2] = {{25, 1}, {50, 2}, {100, 3}, {200, 5}, {400, 8}};
const double RaGain[5][2] = {{28, 1}, {30, 2}, {32, 3}, {34, 5}, {36, 8}};
const double RaScan[6] = {2, 4, 6, 8, 10, 12};

double RdErfc(double X)   // Numerical Recipes erfcc, as radar.js
{
	const double Z = std::fabs(X), T = 1 / (1 + 0.5 * Z);
	const double R = T * std::exp(-Z * Z - 1.26551223 + T * (1.00002368 + T * (0.37409196 + T * (0.09678418 + T * (-0.18628806 + T * (0.27886807 + T * (-1.13520398 + T * (1.48851587 + T * (-0.82215223 + T * 0.17087277)))))))));
	return X >= 0 ? R : 2 - R;
}
double RdGammaQ(int N, double T)
{
	double S = 0, Term = 1;
	for (int k = 0; k < N; ++k)
	{
		if (k)
		{
			Term *= T / k;
		}
		S += Term;
	}
	return std::exp(-T) * S;
}
double RdThreshold(int N, double Pfa)
{
	double Lo = 0, Hi = 50 + 10.0 * N;
	for (int i = 0; i < 200; ++i)
	{
		const double M = (Lo + Hi) / 2;
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
// Steady target, n pulses integrated noncoherently: Poisson mixture of central chi-squares (threshold passed in).
double RdPdSteady(double Snr, int N, double T)
{
	const double Lam = N * Snr, Sd = std::sqrt(Lam);
	if (Lam > 5000)
	{
		const double Z = (N + Lam - T) / std::sqrt(N + 2 * Lam);
		return Z > 8 ? 1 : 0.5 * RdErfc(-Z / std::sqrt(2.0));
	}
	const int J0 = static_cast<int>(std::max(0.0, std::floor(Lam - 12 * Sd - 20))), J1 = static_cast<int>(std::ceil(Lam + 12 * Sd + 40));
	double Lf = 0;
	for (int k = 2; k <= J0; ++k)
	{
		Lf += std::log(static_cast<double>(k));
	}
	// gammaQ(m, T) = e^-T times the m-term prefix sum of T^k/k!; the prefix sums are the same numbers gammaQ builds.
	std::vector<double> Prefix(static_cast<size_t>(N + J1 + 1));
	double Acc = 0, Term = 1;
	for (int k = 0; k < N + J1; ++k)
	{
		if (k)
		{
			Term *= T / k;
		}
		Acc += Term;
		Prefix[static_cast<size_t>(k + 1)] = Acc;
	}
	const double ExpT = std::exp(-T);
	double S = 0;
	for (int j = J0; j <= J1; ++j)
	{
		if (j > J0)
		{
			Lf += std::log(static_cast<double>(j));
		}
		const double W = Lam > 0 ? std::exp(-Lam + j * std::log(Lam) - Lf) : (j == 0 ? 1 : 0);
		if (W > 1e-300)
		{
			S += W * (ExpT * Prefix[static_cast<size_t>(N + j)]);
		}
	}
	return std::min(1.0, S);
}

struct RadarOut
{
	double Score;
	int K, N;
	double Cost;
};

RadarOut RadarRun(int Pi, int Gi, int Si)
{
	const double Kw = RaPower[Pi][0], GDb = RaGain[Gi][0], T = RaScan[Si];
	const double Bw = 1.4 * std::pow(10.0, (33 - GDb) / 10);
	const int N = static_cast<int>(std::max(1.0, std::floor(1200 * T * Bw / 360)));
	const double Lam = cst::C / 2.8e9, Thr = RdThreshold(N, RaPfa);
	const double G = std::pow(10.0, GDb / 10), L = std::pow(10.0, 8.0 / 10);
	RadarOut Out{0, -1, N, RaPower[Pi][1] + RaGain[Gi][1]};
	double St[4] = {1, 0, 0, 0}, Got = 0;
	for (int k = 0;; ++k)
	{
		const double Rm = 60 * RaNmi - 120 * RaKt * k * T;
		if (Rm <= 0)
		{
			break;
		}
		const double Snr = Kw * 1e3 * G * G * Lam * Lam * 1 / (std::pow(4 * cst::Pi, 3) * std::pow(Rm, 4) * 1.380649e-23 * 950 * 1.67e6 * L);
		const double P = RdPdSteady(Snr, N, Thr);
		// 2 of 3: node 1 = one hit, node 2 = hit then miss, node 3 = track.
		double Ns[4] = {0, 0, 0, 0};
		Ns[1] += St[0] * P;
		Ns[0] += St[0] * (1 - P);
		Ns[3] += St[1] * P;
		Ns[2] += St[1] * (1 - P);
		Ns[3] += St[2] * P;
		Ns[0] += St[2] * (1 - P);
		Got += Ns[3];
		St[0] = Ns[0];
		St[1] = Ns[1];
		St[2] = Ns[2];
		St[3] = 0;
		if (Got >= 0.9)
		{
			Out.Score = Rm / RaNmi;
			Out.K = k;
			break;
		}
	}
	return Out;
}

int FindIn(const double (*List)[2], int N, double X)
{
	for (int i = 0; i < N; ++i)
	{
		if (List[i][0] == X)
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
		Si = RaScan[i] == C.N(2) ? i : Si;
	}
	if (Pi < 0 || Gi < 0 || Si < 0)
	{
		return Result::Fail("choose power 25, 50, 100, 200 or 400 kW; gain 28, 30, 32, 34 or 36 dB; scan 2 to 12 s in steps of 2");
	}
	const RadarOut R = RadarRun(Pi, Gi, Si);
	return Result::Record({{"range", R.Score}, {"cost", R.Cost}, {"n", static_cast<double>(R.N)}, {"ok", R.Cost <= 10 ? 1.0 : 0.0}, {"scan", static_cast<double>(R.K)}});
}

// Par: exhaustive search of the admissible designs (cost <= 10), computed once.
Result RadarPar(const Call& C)
{
	(void)C;
	struct Best
	{
		double Score = -1, Kw = 0, GDb = 0, T = 0;
		int N = 0;
	};
	static const Best B = [] {
		Best Out;
		for (int Pi = 0; Pi < 5; ++Pi)
		{
			for (int Gi = 0; Gi < 5; ++Gi)
			{
				if (RaPower[Pi][1] + RaGain[Gi][1] > 10)
				{
					continue;
				}
				for (int Si = 0; Si < 6; ++Si)
				{
					const RadarOut R = RadarRun(Pi, Gi, Si);
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
	return Result::Record({{"range", B.Score}, {"kw", B.Kw}, {"gdb", B.GDb}, {"T", B.T}, {"n", static_cast<double>(B.N)}});
}

Result Stars(const Call& C)
{
	const double Score = C.N(0), Par = C.N(1);
	const bool Higher = C.N(2) != 0;
	const double R = !std::isfinite(Score) ? 0 : Higher ? Score / Par : (Score > 0 ? Par / Score : 0);
	return Result::Of(R >= 0.99 ? 3.0 : R >= 0.95 ? 2.0 : R >= 0.85 ? 1.0 : 0.0);
}

Param TableParam()
{
	Param P = Mat("table", "m");
	P.HasDefault = true;
	P.Default = 0;
	return P;
}
} // namespace flight_detail

void RegisterFlightLibrary(Registry& R)
{
	using namespace flight_detail;
	const double A4 = FlRad(4), KS = KSheet;
	const char* Sp = "1/s";

	// Equation sheet (flight-model/equations.js, equations.html)
	Def(R, "thinairfoil", "aero", "flight-model/equations.js thinAirfoil (Glauert)", "thin-airfoil theory for a NACA 4-digit mean line", "record aL0 (rad), cl, cm4, A1, A2",
		{Num("m", "-", 0.02), Num("p", "-", 0.4), Num("alpha", "rad", A4)}, &ThinAirfoil);
	Def(R, "vortexpanel", "aero", "flight-model/equations.js vortexPanels", "discrete-vortex check: vortex at 1/4, tangency at 3/4 of each panel", "record cl, cm4",
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
		{Num("mdot", "kg/s", 5), Num("pic", "-", 12), Num("T02", "K", 288.15), Num("T04", "K", 1400), Num("etac", "-", 0.85), Num("etat", "-", 0.9), Num("gamma", "-", 1.4), Num("cp", "J/(kg K)", FlCpAir)},
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
		{Num("Kp", "-", 1.5), Num("Ki", "-", 1), Num("Kd", "-", 0.3), Num("thetac", "rad", FlRad(5)), Num("dist", "rad/s^2", 0.5)}, &PidSim);

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
	Def(R, "planquad", "airframe", "flight-model/planform.js planformQuad", "planform check: chord at n + 1 even stations, trapezoidal rule", "record S, MAC, ymac, xmac",
		{Num("n", "-", 1000), TableParam()}, &PlanQuad);
	Def(R, "wingstation", "airframe", "flight-model/planform.js station", "interpolated station at span y", "record xle, xte, c, t, zc (m), up", {Num("y", "m"), TableParam()}, &FnStation);
	Def(R, "stdatm", "airframe", "flight-model/planform.js isa + hydro", "ICAO standard atmosphere 0-20 km, with the hydrostatic integral", "record T, p, rho, a, phyd, rhohyd",
		{Num("h", "m")}, &FnIsa);
	Def(R, "airsizing", "airframe", "flight-model/planform.js sizing", "cruise sizing: speed, q, wing loading, CL, T/W (S = 0: the planform's area)", "record V, W, q, q2, WS, WSkg, CL, TW",
		{Num("mass", "kg", 200e3), Num("h", "m", 10e3), Num("M", "-", 0.8), Num("thrust", "N", 560e3), Num("S", "m^2", 0)}, &Sizing);
	Def(R, "cargohold", "airframe", "flight-model/planform.js holdEnds + holdGrid", "largest clear height for an L x w cargo box inside the body", "record H, x0, x1, floor, ceil, Hscan, xscan (m), sides",
		{Num("L", "m", 19), Num("w", "m", 6), TableParam()}, &Hold);
	Def(R, "bodyvol", "airframe", "flight-model/planform.js loft + signedVolume", "divergence-theorem volume of the lofted body: as exported, fixed winding, exact", "record raw, good, exact (m^3)",
		{Num("M", "-", 70), Num("k", "-", 1), TableParam()}, &BodyVol);

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
