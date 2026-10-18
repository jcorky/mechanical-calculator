// Radar: the range equation and detection (radar/radar.js), radar basics (radar/basics.js) and
// Luneburg lens antennas (lens/lens.js). Each function follows its JS module operation for operation,
// including the modules' own bisections and the mulberry32 generator, so the numbers match.
// SI units throughout; angles in radians (the basics bench shows degrees); "_db" means decibels.
#include "CalcLibrary.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace mcalc
{
namespace radar_detail
{
const double KB = 1.380649e-23;
const double NaN = std::numeric_limits<double>::quiet_NaN();

double Db(double X) { return 10.0 * std::log10(X); }
double Undb(double D) { return std::pow(10.0, D / 10.0); }
bool IsWhole(double X) { return std::isfinite(X) && X == std::floor(X); }

// JavaScript's Math.round: nearest integer, halves toward +infinity.
double JsRound(double X)
{
	const double F = std::floor(X);
	return X - F >= 0.5 ? F + 1.0 : F;
}

// mulberry32, as in the modules' rng(seed), and their Box-Muller gauss(r).
struct Mulberry
{
	uint32_t S;
	explicit Mulberry(double Seed) : S(static_cast<uint32_t>(static_cast<int64_t>(Seed))) {}
	double Next()
	{
		S += 0x6D2B79F5u;
		uint32_t T = S;
		T = (T ^ (T >> 15)) * (T | 1u);
		T ^= T + (T ^ (T >> 7)) * (T | 61u);
		return static_cast<double>(T ^ (T >> 14)) / 4294967296.0;
	}
};
double Gauss(Mulberry& R)
{
	double U = 0.0;
	while (U == 0.0)
	{
		U = R.Next();
	}
	return std::sqrt(-2.0 * std::log(U)) * std::cos(2.0 * cst::Pi * R.Next());
}

/* ---------- range equation ---------- */
double SnrTrack(double Pt, double G, double Lam, double Sig, double R, double Ts, double Bn, double L)
{
	return Pt * G * G * Lam * Lam * Sig / (std::pow(4.0 * cst::Pi, 3) * std::pow(R, 4) * KB * Ts * Bn * L);
}
// Pt and sigma may be zero; every other input of the track form must be positive.
bool TrackArgs(const Call& C)
{
	for (int i = 0; i < 8; ++i)
	{
		const double X = C.N(static_cast<size_t>(i));
		const bool CanBeZero = i == 0 || i == 3;
		if (!std::isfinite(X) || X < 0.0 || (!CanBeZero && X == 0.0))
		{
			return false;
		}
	}
	return true;
}
Result SnrTrackFn(const Call& C)
{
	if (!TrackArgs(C))
	{
		return Result::Fail("needs Pt, sigma >= 0 and G, lambda, R, Ts, Bn, L > 0");
	}
	return Result::Of(SnrTrack(C.N(0), C.N(1), C.N(2), C.N(3), C.N(4), C.N(5), C.N(6), C.N(7)));
}
Result SnrChain(const Call& C)
{
	if (!TrackArgs(C))
	{
		return Result::Fail("needs Pt, sigma >= 0 and G, lambda, R, Ts, Bn, L > 0");
	}
	const double Pt = C.N(0), G = C.N(1), Lam = C.N(2), Sig = C.N(3), R = C.N(4), Ts = C.N(5), Bn = C.N(6), L = C.N(7);
	const double Ae = G * Lam * Lam / (4.0 * cst::Pi), AtTarget = Pt * G / (4.0 * cst::Pi * R * R) / L, Reflected = AtTarget * Sig;
	const double Back = Reflected / (4.0 * cst::Pi * R * R), Pr = Back * Ae, N = KB * Ts * Bn;
	return Result::Record({{"ae", Ae}, {"s_tgt", AtTarget}, {"refl", Reflected}, {"s_back", Back}, {"pr", Pr}, {"noise", N}, {"snr", Pr / N}});
}
Result SnrBudget(const Call& C)
{
	if (!TrackArgs(C) || C.N(0) == 0.0 || C.N(3) == 0.0)
	{
		return Result::Fail("needs every input > 0 (decibels of each term)");
	}
	const double Pt = C.N(0), G = C.N(1), Lam = C.N(2), Sig = C.N(3), R = C.N(4), Ts = C.N(5), Bn = C.N(6), L = C.N(7);
	const double DSig = Db(Sig);
	// Rows as on the lecture's slide: [+ dB, - dB].
	const double Rows[10][2] = {{Db(Pt), 0.0}, {2.0 * Db(G), 0.0}, {0.0, -Db(Lam * Lam)}, {DSig >= 0 ? DSig : 0.0, DSig < 0 ? -DSig : 0.0},
		{0.0, Db(std::pow(4.0 * cst::Pi, 3))}, {0.0, 4.0 * Db(R)}, {-Db(KB), 0.0}, {0.0, Db(Ts)}, {0.0, Db(L)}, {0.0, Db(Bn)}};
	double Plus = 0.0, Minus = 0.0;
	for (const auto& Row : Rows)
	{
		Plus += Row[0];
	}
	for (const auto& Row : Rows)
	{
		Minus += Row[1];
	}
	return Result::Record({{"pt", Rows[0][0]}, {"g2", Rows[1][0]}, {"lam2", -Rows[2][1]}, {"sigma", DSig}, {"fourpi3", -Rows[4][1]}, {"r4", -Rows[5][1]},
		{"inv_k", Rows[6][0]}, {"ts", -Rows[7][1]}, {"loss", -Rows[8][1]}, {"bn", -Rows[9][1]}, {"plus", Plus}, {"minus", Minus}, {"total", Plus - Minus}});
}
Result RangeFor(const Call& C)
{
	const double S = C.N(0), Pt = C.N(1), G = C.N(2), Lam = C.N(3), Sig = C.N(4), Ts = C.N(5), Bn = C.N(6), L = C.N(7);
	if (!(S > 0 && Pt >= 0 && G > 0 && Lam > 0 && Sig >= 0 && Ts > 0 && Bn > 0 && L > 0))
	{
		return Result::Fail("needs snr, G, lambda, Ts, Bn, L > 0 and Pt, sigma >= 0");
	}
	return Result::Of(std::pow(Pt * G * G * Lam * Lam * Sig / (std::pow(4.0 * cst::Pi, 3) * KB * Ts * Bn * L * S), 0.25));
}

/* ---------- search form ---------- */
bool AllPositive(const Call& C, size_t N)
{
	for (size_t i = 0; i < N; ++i)
	{
		if (!(C.N(i) > 0.0) || !std::isfinite(C.N(i)))
		{
			return false;
		}
	}
	return true;
}
double SnrSearchCore(double Pav, double A, double Tsc, double Sig, double Om, double R, double Ts, double L)
{
	return Pav * A * Tsc * Sig / (4.0 * cst::Pi * Om * std::pow(R, 4) * KB * Ts * L);
}
Result SnrSearch(const Call& C)
{
	if (!AllPositive(C, 8))
	{
		return Result::Fail("every input must be > 0");
	}
	return Result::Of(SnrSearchCore(C.N(0), C.N(1), C.N(2), C.N(3), C.N(4), C.N(5), C.N(6), C.N(7)));
}
Result PavSearch(const Call& C)
{
	if (!AllPositive(C, 8))
	{
		return Result::Fail("every input must be > 0");
	}
	const double R = C.N(0), A = C.N(1), Tsc = C.N(2), Sig = C.N(3), Om = C.N(4), Ts = C.N(5), L = C.N(6), S = C.N(7);
	return Result::Of(4.0 * cst::Pi * Om * std::pow(R, 4) * KB * Ts * L * S / (A * Tsc * Sig));
}
Result RangeSearch(const Call& C)
{
	if (!AllPositive(C, 8))
	{
		return Result::Fail("every input must be > 0");
	}
	const double Pav = C.N(0), A = C.N(1), Tsc = C.N(2), Sig = C.N(3), Om = C.N(4), Ts = C.N(5), L = C.N(6), S = C.N(7);
	return Result::Of(std::pow(Pav * A * Tsc * Sig / (4.0 * cst::Pi * Om * KB * Ts * L * S), 0.25));
}
Result SearchTrack(const Call& C)
{
	if (!AllPositive(C, 11))
	{
		return Result::Fail("every input must be > 0");
	}
	const double Pav = C.N(0), A = C.N(1), Tsc = C.N(2), Sig = C.N(3), Om = C.N(4), R = C.N(5), Ts = C.N(6), L = C.N(7), Lam = C.N(8), Tau = C.N(9), Prf = C.N(10);
	const double G = 4.0 * cst::Pi * A / (Lam * Lam), Ob = 4.0 * cst::Pi / G, Pt = Pav / (Tau * Prf), N = Prf * Tsc * Ob / Om;
	return Result::Record({{"snr", SnrTrack(Pt, G, Lam, Sig, R, Ts, 1.0 / Tau, L) * N}, {"n", N}, {"beams", Om / Ob}});
}

/* ---------- system noise temperature ---------- */
bool NoiseArgs(const Call& C) { return C.N(0) >= 0 && C.N(1) >= 1 && C.N(2) >= 0 && C.N(3) >= 1; }
Result Tsys(const Call& C)
{
	if (!NoiseArgs(C))
	{
		return Result::Fail("needs Ta, Tline >= 0 and the ratios Lr, F >= 1");
	}
	const double Ta = C.N(0), Lr = C.N(1), Tphys = C.N(2), F = C.N(3);
	return Result::Record({{"tr", (Lr - 1) * Tphys}, {"te", (F - 1) * cst::T0Noise}, {"ts", Ta + (Lr - 1) * Tphys + Lr * (F - 1) * cst::T0Noise}});
}
Result TsysTrace(const Call& C)
{
	if (!NoiseArgs(C))
	{
		return Result::Fail("needs Ta, Tline >= 0 and the ratios Lr, F >= 1");
	}
	const double Ta = C.N(0), Lr = C.N(1), Tphys = C.N(2), F = C.N(3);
	const double AfterLine = Ta / Lr + Tphys * (1 - 1 / Lr), AtRx = AfterLine + (F - 1) * cst::T0Noise;
	return Result::Of(AtRx * Lr);
}

/* ---------- detection ---------- */
// Numerical Recipes erfcc, as the module uses it (relative error below 1.2e-7).
double ErfcNr(double X)
{
	const double Z = std::fabs(X), T = 1 / (1 + 0.5 * Z);
	const double R = T * std::exp(-Z * Z - 1.26551223 + T * (1.00002368 + T * (0.37409196 + T * (0.09678418 + T * (-0.18628806 + T * (0.27886807 + T * (-1.13520398 + T * (1.48851587 + T * (-0.82215223 + T * 0.17087277)))))))));
	return X >= 0 ? R : 2 - R;
}
// Regularized upper incomplete gamma Q(n, t) for integer n.
double GammaQ(int N, double T)
{
	double S = 0.0, Term = 1.0;
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
double Threshold(int N, double Pfa)
{
	double Lo = 0.0, Hi = 50.0 + 10.0 * N;
	for (int i = 0; i < 200; ++i)
	{
		const double M = (Lo + Hi) / 2;
		if (GammaQ(N, M) > Pfa)
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
// Poisson mixture sum_j Pois(j; Lam) Q(N + j, T): the steady-target Pd (and Marcum's Q_N).
// Q(N + j, T) is built up incrementally; each partial sum is exactly the module's gammaQ(N + j, T).
double PoissonGammaMix(double Lam, int N, double T)
{
	const double Sd = std::sqrt(Lam);
	if (Lam > 5000)
	{
		const double Z = (N + Lam - T) / std::sqrt(N + 2 * Lam);
		return Z > 8 ? 1.0 : 0.5 * ErfcNr(-Z / std::sqrt(2.0));
	}
	const long J0 = static_cast<long>(std::max(0.0, std::floor(Lam - 12 * Sd - 20))), J1 = static_cast<long>(std::ceil(Lam + 12 * Sd + 40));
	double Lf = 0.0;
	for (long k = 2; k <= J0; ++k)
	{
		Lf += std::log(static_cast<double>(k));
	}
	const double ET = std::exp(-T);
	long K = 0;
	double Partial = 0.0, Term = 1.0;
	double S = 0.0;
	for (long j = J0; j <= J1; ++j)
	{
		if (j > J0)
		{
			Lf += std::log(static_cast<double>(j));
		}
		const double W = Lam > 0 ? std::exp(-Lam + j * std::log(Lam) - Lf) : (j == 0 ? 1.0 : 0.0);
		const long Need = N + j;
		while (K < Need)
		{
			if (K)
			{
				Term *= T / K;
			}
			Partial += Term;
			++K;
		}
		if (W > 1e-300)
		{
			S += W * (ET * Partial);
		}
	}
	return std::min(1.0, S);
}
double PdSteady(double Snr, int N, double Pfa) { return PoissonGammaMix(N * Snr, N, Threshold(N, Pfa)); }
double SnrRequired(double Pd, int N, double Pfa)
{
	double Lo = -20.0, Hi = 30.0;
	for (int i = 0; i < 70; ++i)
	{
		const double M = (Lo + Hi) / 2;
		if (PdSteady(Undb(M), N, Pfa) < Pd)
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
const int MaxPulses = 500;
const char* const PulsesMessage = "n must be a whole number of pulses, 1 to 500";
bool PulsesOk(double N) { return IsWhole(N) && N >= 1 && N <= MaxPulses; }
bool ProbOk(double P) { return P > 0.0 && P < 1.0; }

Result PfaThr(const Call& C)
{
	if (!PulsesOk(C.N(1)) || !(C.N(0) >= 0))
	{
		return Result::Fail(!(C.N(0) >= 0) ? "threshold must be >= 0" : PulsesMessage);
	}
	return Result::Of(GammaQ(C.Int(1), C.N(0)));
}
Result DetThresh(const Call& C)
{
	if (!PulsesOk(C.N(0)) || !ProbOk(C.N(1)))
	{
		return Result::Fail(!ProbOk(C.N(1)) ? "pfa must be between 0 and 1" : PulsesMessage);
	}
	return Result::Of(Threshold(C.Int(0), C.N(1)));
}
Result PdSteadyFn(const Call& C)
{
	if (!(C.N(0) >= 0) || !PulsesOk(C.N(1)) || !ProbOk(C.N(2)))
	{
		return Result::Fail(!(C.N(0) >= 0) ? "snr is a power ratio >= 0" : !ProbOk(C.N(2)) ? "pfa must be between 0 and 1" : PulsesMessage);
	}
	return Result::Of(PdSteady(C.N(0), C.Int(1), C.N(2)));
}
Result PdSw1(const Call& C)
{
	if (!(C.N(0) >= 0) || !ProbOk(C.N(1)))
	{
		return Result::Fail("needs snr >= 0 (power ratio) and 0 < pfa < 1");
	}
	return Result::Of(std::pow(C.N(1), 1 / (1 + C.N(0))));
}
Result SnrReqDb(const Call& C)
{
	if (!ProbOk(C.N(0)) || !PulsesOk(C.N(1)) || !ProbOk(C.N(2)))
	{
		return Result::Fail(PulsesOk(C.N(1)) ? "pd and pfa must be between 0 and 1" : PulsesMessage);
	}
	return Result::Of(SnrRequired(C.N(0), C.Int(1), C.N(2)));
}
Result Albersheim(const Call& C)
{
	const double Pd = C.N(0), Pfa = C.N(1), N = C.N(2);
	if (!ProbOk(Pd) || !ProbOk(Pfa) || !(N >= 1))
	{
		return Result::Fail("needs 0 < pd, pfa < 1 and n >= 1");
	}
	const double A = std::log(0.62 / Pfa), B = std::log(Pd / (1 - Pd));
	return Result::Of(-5 * std::log10(N) + (6.2 + 4.54 / std::sqrt(N + 0.44)) * std::log10(A + 0.12 * A * B + 1.7 * B));
}
Result IntGain(const Call& C)
{
	if (!PulsesOk(C.N(0)) || !ProbOk(C.N(1)) || !ProbOk(C.N(2)))
	{
		return Result::Fail(PulsesOk(C.N(0)) ? "pd and pfa must be between 0 and 1" : PulsesMessage);
	}
	const int N = C.Int(0);
	const double Nc = SnrRequired(C.N(1), 1, C.N(2)) - SnrRequired(C.N(1), N, C.N(2)), Coh = Db(N);
	return Result::Record({{"nc", Nc}, {"coh", Coh}, {"loss", Coh - Nc}});
}
// Marcum Q_M(a, b) = sum_j Pois(j; a^2/2) Q(M + j, b^2/2), the steady-target Pd with a = sqrt(2 n S/N).
Result MarcumQ(const Call& C)
{
	const double A = C.N(0), B = C.N(1);
	if (!(A >= 0) || !(B >= 0) || !PulsesOk(C.N(2)) || A * A / 2 > 1e6)
	{
		return Result::Fail(PulsesOk(C.N(2)) ? "needs a, b >= 0 (a^2/2 up to 1e6)" : "order m must be a whole number, 1 to 500");
	}
	return Result::Of(PoissonGammaMix(A * A / 2, C.Int(2), B * B / 2));
}
Result PdMc(const Call& C)
{
	const double Snr = C.N(0), Pfa = C.N(2), Trials = C.N(3), Seed = C.N(5);
	const bool Sw1 = C.N(4) != 0.0;
	if (!(Snr >= 0) || !PulsesOk(C.N(1)) || !ProbOk(Pfa) || !IsWhole(Trials) || Trials < 1 || !IsWhole(Seed) || Seed < 0 || Seed > 4294967295.0)
	{
		return Result::Fail("needs snr >= 0, 1 <= n <= 500, 0 < pfa < 1, whole trials >= 1 and a whole seed >= 0");
	}
	const int N = C.Int(1);
	if (Trials * N > 4e5)
	{
		return Result::Fail("trials x n is limited to 4e5 pulses");
	}
	Mulberry R(Seed);
	const double T = Threshold(N, Pfa);
	long Hits = 0;
	const long Count = static_cast<long>(Trials);
	for (long t = 0; t < Count; ++t)
	{
		double S = 0.0;
		const double Amp = Sw1 ? std::sqrt(-Snr * std::log(1 - R.Next())) : std::sqrt(Snr);
		for (int k = 0; k < N; ++k)
		{
			const double Ph = 2.0 * cst::Pi * R.Next();
			const double Re = Amp * std::cos(Ph) + Gauss(R) / std::sqrt(2.0);
			const double Im = Amp * std::sin(Ph) + Gauss(R) / std::sqrt(2.0);
			S += Re * Re + Im * Im;
		}
		if (S > T)
		{
			++Hits;
		}
	}
	return Result::Of(static_cast<double>(Hits) / Trials);
}

/* ---------- GERT ---------- */
// The module's Gaussian elimination with partial pivoting; fails on a zero pivot.
bool GertSolve(std::vector<std::vector<double>> M, const std::vector<double>& B, std::vector<double>* X)
{
	const size_t N = B.size();
	for (size_t r = 0; r < N; ++r)
	{
		M[r].push_back(B[r]);
	}
	for (size_t c = 0; c < N; ++c)
	{
		size_t P = c;
		for (size_t r = c + 1; r < N; ++r)
		{
			if (std::fabs(M[r][c]) > std::fabs(M[P][c]))
			{
				P = r;
			}
		}
		std::swap(M[c], M[P]);
		if (std::fabs(M[c][c]) < 1e-300)
		{
			return false;
		}
		for (size_t r = c + 1; r < N; ++r)
		{
			const double F = M[r][c] / M[c][c];
			if (F != 0.0)
			{
				for (size_t j = c; j <= N; ++j)
				{
					M[r][j] -= F * M[c][j];
				}
			}
		}
	}
	X->assign(N, 0.0);
	for (size_t rr = N; rr-- > 0;)
	{
		double S = M[rr][N];
		for (size_t j = rr + 1; j < N; ++j)
		{
			S -= M[rr][j] * (*X)[j];
		}
		(*X)[rr] = S / M[rr][rr];
	}
	return true;
}
struct GertBranch
{
	int I, J;
	double P, T;
};
// Equivalent transmittance with constant delays: (I - A0) x0 = b0, x1 and x2 from its derivatives.
bool GertMoments(int Nodes, const std::vector<GertBranch>& Br, int Start, int End, double* P, double* Mean, double* Variance, std::string* Error)
{
	std::vector<int> Idx(static_cast<size_t>(Nodes), -1);
	int M = 0;
	for (int k = 0; k < Nodes; ++k)
	{
		if (k != End)
		{
			Idx[static_cast<size_t>(k)] = M++;
		}
	}
	const size_t Mz = static_cast<size_t>(M);
	std::vector<std::vector<double>> A0(Mz, std::vector<double>(Mz, 0.0)), A1 = A0, A2 = A0;
	std::vector<double> B0(Mz, 0.0), B1 = B0, B2 = B0;
	for (const GertBranch& G : Br)
	{
		const size_t R = static_cast<size_t>(Idx[static_cast<size_t>(G.I)]);
		if (G.J == End)
		{
			B0[R] += G.P;
			B1[R] += G.P * G.T;
			B2[R] += G.P * G.T * G.T;
		}
		else
		{
			const size_t Cc = static_cast<size_t>(Idx[static_cast<size_t>(G.J)]);
			A0[R][Cc] += G.P;
			A1[R][Cc] += G.P * G.T;
			A2[R][Cc] += G.P * G.T * G.T;
		}
	}
	std::vector<std::vector<double>> Mm = A0;
	for (size_t i = 0; i < Mz; ++i)
	{
		for (size_t j = 0; j < Mz; ++j)
		{
			Mm[i][j] = (i == j ? 1.0 : 0.0) - A0[i][j];
		}
	}
	auto Mv = [Mz](const std::vector<std::vector<double>>& A, const std::vector<double>& X) {
		std::vector<double> Out(Mz, 0.0);
		for (size_t i = 0; i < Mz; ++i)
		{
			double S = 0.0;
			for (size_t j = 0; j < Mz; ++j)
			{
				S += A[i][j] * X[j];
			}
			Out[i] = S;
		}
		return Out;
	};
	std::vector<double> X0, X1, X2;
	if (!GertSolve(Mm, B0, &X0))
	{
		*Error = "I - A is singular: some node cannot reach the end";
		return false;
	}
	std::vector<double> R1 = Mv(A1, X0);
	for (size_t i = 0; i < Mz; ++i)
	{
		R1[i] += B1[i];
	}
	GertSolve(Mm, R1, &X1);
	const std::vector<double> V2 = Mv(A2, X0), W1 = Mv(A1, X1);
	std::vector<double> R2(Mz, 0.0);
	for (size_t i = 0; i < Mz; ++i)
	{
		R2[i] = V2[i] + 2 * W1[i] + B2[i];
	}
	GertSolve(Mm, R2, &X2);
	const size_t K = static_cast<size_t>(Idx[static_cast<size_t>(Start)]);
	*P = X0[K];
	if (!(*P > 0.0))
	{
		*Error = "the end node is never reached from the start";
		return false;
	}
	*Mean = X1[K] / *P;
	*Variance = X2[K] / *P - *Mean * *Mean;
	return true;
}
Result Gert(const Call& C)
{
	const Matrix& B = C.M(0);
	if (B.Cols != 4 || B.Rows < 1)
	{
		return Result::Fail("branches are rows [from, to, p, t] (k x 4), nodes numbered from 1");
	}
	int Nodes = 0;
	std::vector<GertBranch> Br;
	for (int r = 0; r < B.Rows; ++r)
	{
		const double I = B(r, 0), J = B(r, 1), P = B(r, 2), T = B(r, 3);
		if (!IsWhole(I) || !IsWhole(J) || I < 1 || J < 1 || I > 60 || J > 60)
		{
			return Result::Fail("node numbers must be whole numbers 1 to 60");
		}
		if (!(P >= 0 && P <= 1) || !std::isfinite(T))
		{
			return Result::Fail("branch probabilities must be in [0, 1] and times finite");
		}
		Br.push_back({static_cast<int>(I) - 1, static_cast<int>(J) - 1, P, T});
		Nodes = std::max(Nodes, static_cast<int>(std::max(I, J)));
	}
	const int Start = C.Int(1) - 1, End = C.N(2) == 0.0 ? Nodes - 1 : C.Int(2) - 1;
	if (!IsWhole(C.N(1)) || !IsWhole(C.N(2)) || Start < 0 || Start >= Nodes || End < 0 || End >= Nodes || Start == End)
	{
		return Result::Fail("start and finish must be different nodes of the network (finish 0 = the last node)");
	}
	for (const GertBranch& G : Br)
	{
		if (G.I == End)
		{
			return Result::Fail("the finish node is absorbing: no branch may leave it");
		}
	}
	double P = 0, Mean = 0, Var = 0;
	std::string E;
	if (!GertMoments(Nodes, Br, Start, End, &P, &Mean, &Var, &E))
	{
		return Result::Fail(E);
	}
	return Result::Record({{"p", P}, {"mean", Mean}, {"var", Var}, {"sd", std::sqrt(Var)}});
}
// Track-confirmation rules as in rule(): 1 = 2 in a row, 2 = 2 of 3, 3 = 3 in a row. 0-based nodes.
bool RuleNet(int Rule, double P, double T, int* Nodes, int* End, std::vector<GertBranch>* Br)
{
	if (Rule == 1)
	{
		*Nodes = 3;
		*End = 2;
		*Br = {{0, 1, P, T}, {0, 0, 1 - P, T}, {1, 2, P, T}, {1, 0, 1 - P, T}};
	}
	else if (Rule == 3)
	{
		*Nodes = 4;
		*End = 3;
		*Br = {{0, 1, P, T}, {0, 0, 1 - P, T}, {1, 2, P, T}, {1, 0, 1 - P, T}, {2, 3, P, T}, {2, 0, 1 - P, T}};
	}
	else if (Rule == 2)
	{
		*Nodes = 4;
		*End = 3;
		*Br = {{0, 1, P, T}, {0, 0, 1 - P, T}, {1, 3, P, T}, {1, 2, 1 - P, T}, {2, 3, P, T}, {2, 0, 1 - P, T}};
	}
	else
	{
		return false;
	}
	return true;
}
const char* const RuleMessage = "rule: 1 = 2 in a row, 2 = 2 of 3, 3 = 3 in a row; pd in (0, 1], T > 0";
Result TTrack(const Call& C)
{
	const double Pd = C.N(0), T = C.N(1);
	int Nodes = 0, End = 0;
	std::vector<GertBranch> Br;
	if (!(Pd > 0 && Pd <= 1) || !(T > 0) || !IsWhole(C.N(2)) || !RuleNet(C.Int(2), Pd, T, &Nodes, &End, &Br))
	{
		return Result::Fail(RuleMessage);
	}
	double P = 0, Mean = 0, Var = 0;
	std::string E;
	if (!GertMoments(Nodes, Br, 0, End, &P, &Mean, &Var, &E))
	{
		return Result::Fail(E);
	}
	// Scans until 90% of tracks have formed, by propagating probability scan by scan (scanPmf).
	std::vector<double> V(static_cast<size_t>(Nodes), 0.0);
	V[0] = 1.0;
	double Cum = 0.0, P90 = NaN;
	for (long s = 1; s <= 1000000; ++s)
	{
		double Left = 0.0;
		for (double X : V)
		{
			Left += X;
		}
		if (!(Left > 1e-13))
		{
			break;
		}
		std::vector<double> Nv(static_cast<size_t>(Nodes), 0.0);
		for (const GertBranch& G : Br)
		{
			if (G.I != End)
			{
				Nv[static_cast<size_t>(G.J)] += V[static_cast<size_t>(G.I)] * G.P;
			}
		}
		Cum += Nv[static_cast<size_t>(End)];
		Nv[static_cast<size_t>(End)] = 0.0;
		V = Nv;
		if (Cum >= 0.9)
		{
			P90 = static_cast<double>(s);
			break;
		}
	}
	return Result::Record({{"p", P}, {"mean", Mean}, {"var", Var}, {"sd", std::sqrt(Var)}, {"p90", P90}});
}
Result TrkNet(const Call& C)
{
	int Nodes = 0, End = 0;
	std::vector<GertBranch> Br;
	if (!(C.N(0) >= 0 && C.N(0) <= 1) || !(C.N(1) > 0) || !IsWhole(C.N(2)) || !RuleNet(C.Int(2), C.N(0), C.N(1), &Nodes, &End, &Br))
	{
		return Result::Fail(RuleMessage);
	}
	Matrix M(static_cast<int>(Br.size()), 4);
	for (size_t i = 0; i < Br.size(); ++i)
	{
		const int r = static_cast<int>(i);
		M(r, 0) = Br[i].I + 1;
		M(r, 1) = Br[i].J + 1;
		M(r, 2) = Br[i].P;
		M(r, 3) = Br[i].T;
	}
	return Result::Of(M);
}

/* ---------- radar basics ---------- */
Result MfRange(const Call& C)
{
	const double R = C.N(0), Fs = C.N(1), Tau = C.N(2), Noise = C.N(3), Seed = C.N(4);
	if (!(R >= 0) || !(Fs > 0) || !(Tau > 0) || !(Noise >= 0) || !IsWhole(Seed) || Seed < 0 || Seed > 4294967295.0)
	{
		return Result::Fail("needs R >= 0, fs, tau > 0, noise >= 0 and a whole seed");
	}
	const double Lf = JsRound(Tau * Fs), K0f = JsRound(2 * R / cst::C * Fs);
	if (Lf < 1 || K0f + Lf + 200 > 500000)
	{
		return Result::Fail("the pulse must span a sample, and the record (2R/c + tau) fs is limited to 500 000 samples");
	}
	const long L = static_cast<long>(Lf), K0 = static_cast<long>(K0f), N = K0 + L + 200;
	Mulberry Rg(Seed);
	std::vector<double> Re(static_cast<size_t>(N), 0.0), Im(static_cast<size_t>(N), 0.0);
	for (long i = 0; i < L; ++i)
	{
		Re[static_cast<size_t>(K0 + i)] += std::cos(0.7);
		Im[static_cast<size_t>(K0 + i)] += std::sin(0.7);
	}
	for (size_t i = 0; i < static_cast<size_t>(N); ++i)
	{
		Re[i] += Gauss(Rg) * std::sqrt(Noise);
		Im[i] += Gauss(Rg) * std::sqrt(Noise);
	}
	double Best = -1, Sr = 0, Si = 0;
	long Arg = 0;
	for (long i = 0; i < L; ++i)
	{
		Sr += Re[static_cast<size_t>(i)];
		Si += Im[static_cast<size_t>(i)];
	}
	for (long k = 0; k + L <= N; ++k)
	{
		const double M = Sr * Sr + Si * Si;
		if (M > Best)
		{
			Best = M;
			Arg = k;
		}
		if (k + L < N)
		{
			Sr += Re[static_cast<size_t>(k + L)] - Re[static_cast<size_t>(k)];
			Si += Im[static_cast<size_t>(k + L)] - Im[static_cast<size_t>(k)];
		}
	}
	return Result::Record({{"est", Arg / Fs * cst::C / 2}, {"sample", cst::C / (2 * Fs)}});
}
Result ChirpRes(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("bandwidth must be > 0");
	}
	return Result::Of(cst::C / (2 * C.N(0)));
}
// Two equal echoes sepM apart through a chirp's matched filter: resolved if the dip between peaks is below 0.8.
Result Chirp2(const Call& C)
{
	const double Sep = C.N(0), B = C.N(1), T = C.N(2), Fs = C.N(3);
	if (!(Sep >= 0) || !(B > 0) || !(T > 0) || !(Fs > 0))
	{
		return Result::Fail("needs sep >= 0 and B, T, fs > 0");
	}
	const double Nf = JsRound(T * Fs), Lagf = JsRound(2 * Sep / cst::C * Fs);
	if (Nf < 2 || (Lagf + 81) * Nf > 5e6)
	{
		return Result::Fail("T fs must be at least 2 samples, and (2 sep/c fs + 81) T fs at most 5e6");
	}
	const long N = static_cast<long>(Nf), Lag = static_cast<long>(Lagf), Mm = N + Lag + 400;
	std::vector<double> Cr(static_cast<size_t>(N)), Ci(static_cast<size_t>(N));
	for (long i = 0; i < N; ++i)
	{
		const double Tt = i / Fs - T / 2, Ph = cst::Pi * B / T * Tt * Tt;
		Cr[static_cast<size_t>(i)] = std::cos(Ph);
		Ci[static_cast<size_t>(i)] = std::sin(Ph);
	}
	std::vector<double> Xr(static_cast<size_t>(Mm), 0.0), Xi(static_cast<size_t>(Mm), 0.0);
	const long Offsets[2] = {200, 200 + Lag};
	for (long Off : Offsets)
	{
		for (long i = 0; i < N; ++i)
		{
			Xr[static_cast<size_t>(Off + i)] += Cr[static_cast<size_t>(i)];
			Xi[static_cast<size_t>(Off + i)] += Ci[static_cast<size_t>(i)];
		}
	}
	std::vector<double> Out;
	const long Lo = 200 - 40, Hi = 200 + Lag + 40;
	for (long k = Lo; k <= Hi; ++k)
	{
		double A = 0, Bb = 0;
		for (long i = 0; i < N; ++i)
		{
			const long j = k + i;
			if (j >= Mm)
			{
				break;
			}
			const size_t Jz = static_cast<size_t>(j), Iz = static_cast<size_t>(i);
			A += Xr[Jz] * Cr[Iz] + Xi[Jz] * Ci[Iz];
			Bb += Xi[Jz] * Cr[Iz] - Xr[Jz] * Ci[Iz];
		}
		Out.push_back(std::hypot(A, Bb));
	}
	const double Mx = *std::max_element(Out.begin(), Out.end());
	std::vector<double> Y(Out.size());
	for (size_t i = 0; i < Out.size(); ++i)
	{
		Y[i] = Out[i] / Mx;
	}
	std::vector<size_t> Peaks;
	for (size_t i = 1; i + 1 < Y.size(); ++i)
	{
		if (Y[i] > Y[i - 1] && Y[i] >= Y[i + 1] && Y[i] > 0.5)
		{
			Peaks.push_back(i);
		}
	}
	bool Resolved = false;
	if (Peaks.size() >= 2)
	{
		std::stable_sort(Peaks.begin(), Peaks.end(), [&Y](size_t P, size_t Q) { return Y[Q] - Y[P] < 0; });
		size_t A = std::min(Peaks[0], Peaks[1]), Bi = std::max(Peaks[0], Peaks[1]);
		double Dip = Y[A];
		for (size_t i = A; i <= Bi; ++i)
		{
			Dip = std::min(Dip, Y[i]);
		}
		Resolved = Dip < 0.8 * std::min(Y[A], Y[Bi]);
	}
	return Result::Record({{"resolved", Resolved ? 1.0 : 0.0}, {"res", cst::C / (2 * B)}, {"ratio", B * T}, {"mps", cst::C / (2 * Fs)}});
}
double UnambRange(double Prf) { return cst::C / (2 * Prf); }
double UnambVel(double Prf, double Lam) { return Lam * Prf / 4; }
Result RUnamb(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("PRF must be > 0");
	}
	return Result::Of(UnambRange(C.N(0)));
}
Result VUnamb(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("PRF and wavelength must be > 0");
	}
	return Result::Of(UnambVel(C.N(0), C.N(1)));
}
Result RFold(const Call& C)
{
	if (!(C.N(1) > 0))
	{
		return Result::Fail("PRF must be > 0");
	}
	return Result::Of(std::fmod(C.N(0), UnambRange(C.N(1))));
}
Result VAlias(const Call& C)
{
	const double V = C.N(0);
	if (!(C.N(1) > 0) || !(C.N(2) > 0))
	{
		return Result::Fail("PRF and wavelength must be > 0");
	}
	const double Vu = UnambVel(C.N(1), C.N(2));
	return Result::Of(std::fmod(std::fmod(V + Vu, 2 * Vu) + 2 * Vu, 2 * Vu) - Vu);
}
Result DopShift(const Call& C)
{
	if (!(C.N(1) > 0))
	{
		return Result::Fail("wavelength must be > 0");
	}
	return Result::Of(2 * C.N(0) / C.N(1));
}
// Hann-windowed N-pulse DFT of a simulated echo: the speed the radar reads, aliasing included.
Result DopDft(const Call& C)
{
	const double V = C.N(0), Prf = C.N(1), Lam = C.N(2), Nd = C.N(3), Seed = C.N(4);
	if (!(Prf > 0) || !(Lam > 0) || !IsWhole(Nd) || Nd < 2 || Nd > 1024 || !IsWhole(Seed) || Seed < 0 || Seed > 4294967295.0)
	{
		return Result::Fail("needs PRF, wavelength > 0, a whole N from 2 to 1024 and a whole seed");
	}
	const int N = static_cast<int>(Nd);
	Mulberry R(Seed);
	const double Fd = 2 * V / Lam;
	std::vector<double> Xr, Xi;
	for (int n = 0; n < N; ++n)
	{
		const double W = 0.5 - 0.5 * std::cos(2 * cst::Pi * n / N), Ph = 2 * cst::Pi * Fd * n / Prf;
		Xr.push_back(W * (std::cos(Ph) + 0.1 * Gauss(R)));
		Xi.push_back(W * (std::sin(Ph) + 0.1 * Gauss(R)));
	}
	double Best = -1;
	int Kb = 0;
	for (int k = 0; k < N; ++k)
	{
		double A = 0, B = 0;
		for (int n = 0; n < N; ++n)
		{
			const double Ang = -2 * cst::Pi * k * n / N;
			const size_t Nz = static_cast<size_t>(n);
			A += Xr[Nz] * std::cos(Ang) - Xi[Nz] * std::sin(Ang);
			B += Xr[Nz] * std::sin(Ang) + Xi[Nz] * std::cos(Ang);
		}
		const double P = A * A + B * B;
		if (P > Best)
		{
			Best = P;
			Kb = k;
		}
	}
	const double F = (Kb < N / 2.0 ? Kb : Kb - N) * Prf / N;
	return Result::Record({{"est", F * Lam / 2}, {"bin", Prf / N * Lam / 2}});
}
Result PrfTrade(const Call& C)
{
	const double Prf = C.N(0), Lam = C.N(1);
	if (!(Prf > 0) || !(Lam > 0))
	{
		return Result::Fail("PRF and wavelength must be > 0");
	}
	const double Ru = UnambRange(Prf), Vu = UnambVel(Prf, Lam);
	return Result::Record({{"ru", Ru}, {"vu", Vu}, {"prod", Ru * Vu}, {"regime", Prf < 4e3 ? 1.0 : Prf < 30e3 ? 2.0 : 3.0}});
}
// Two-pulse canceller on N simulated pulses against 4 sin^2(pi fd / PRF).
Result Mti(const Call& C)
{
	const double V = C.N(0), Prf = C.N(1), Lam = C.N(2), Nd = C.N(3);
	if (!(Prf > 0) || !(Lam > 0) || !IsWhole(Nd) || Nd < 2 || Nd > 1e6)
	{
		return Result::Fail("needs PRF, wavelength > 0 and a whole N from 2 to 1e6");
	}
	const long N = static_cast<long>(Nd);
	const double Fd = 2 * V / Lam;
	double Num = 0;
	for (long n = 1; n < N; ++n)
	{
		const double A = 2 * cst::Pi * Fd * n / Prf, B = 2 * cst::Pi * Fd * (n - 1) / Prf;
		Num += std::pow(std::cos(A) - std::cos(B), 2) + std::pow(std::sin(A) - std::sin(B), 2);
	}
	const double Measured = Num / (N - 1);
	return Result::Record({{"meas", Measured}, {"pred", 4 * std::pow(std::sin(cst::Pi * Fd / Prf), 2)}, {"gain_db", 10 * std::log10(std::max(Measured, 1e-12) / 2)},
		{"vblind", Lam * Prf / 2}});
}
Result Clutter(const Call& C)
{
	const double V = C.N(0), Lam = C.N(1), Nd = C.N(2), Seed = C.N(3);
	if (!(Lam > 0) || !IsWhole(Nd) || Nd < 1 || Nd > 1e6 || !IsWhole(Seed) || Seed < 0 || Seed > 4294967295.0)
	{
		return Result::Fail("needs wavelength > 0, a whole n from 1 to 1e6 and a whole seed");
	}
	Mulberry R(Seed);
	double Lo = std::numeric_limits<double>::infinity(), Hi = -Lo;
	for (long i = 0; i < static_cast<long>(Nd); ++i)
	{
		const double Fd = 2 * V * std::cos(cst::Pi * R.Next()) / Lam;
		Lo = std::min(Lo, Fd);
		Hi = std::max(Hi, Fd);
	}
	return Result::Record({{"lo", Lo}, {"hi", Hi}, {"edge", 2 * V / Lam}});
}
double Sinc(double X) { return X == 0 ? 1.0 : std::sin(cst::Pi * X) / (cst::Pi * X); }
Result Beam3db(const Call& C)
{
	const double D = C.N(0), Lam = C.N(1);
	if (!(D > 0) || !(Lam > 0))
	{
		return Result::Fail("aperture and wavelength must be > 0");
	}
	double Lo = 0, Hi = std::min(1.0, Lam / D);
	for (int i = 0; i < 80; ++i)
	{
		const double M = (Lo + Hi) / 2;
		if (std::pow(Sinc(D * M / Lam), 2) > 0.5)
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	return Result::Record({{"num", 2 * std::asin(Lo)}, {"rule", 0.886 * Lam / D}});
}
double ArrayFactor(int N, double DLam, double Steer, double S)
{
	const double S0 = std::sin(Steer);
	double A = 0, B = 0;
	for (int n = 0; n < N; ++n)
	{
		const double Ph = 2 * cst::Pi * DLam * n * (S - S0);
		A += std::cos(Ph);
		B += std::sin(Ph);
	}
	return std::hypot(A, B) / N;
}
bool ArrayArgs(const Call& C) { return IsWhole(C.N(0)) && C.N(0) >= 1 && C.N(0) <= 4096 && C.N(1) > 0 && std::isfinite(C.N(2)); }
Result AFactor(const Call& C)
{
	if (!ArrayArgs(C))
	{
		return Result::Fail("needs a whole number of elements 1 to 4096 and spacing > 0 (wavelengths)");
	}
	return Result::Of(ArrayFactor(C.Int(0), C.N(1), C.N(2), C.N(3)));
}
Result AScan(const Call& C)
{
	if (!ArrayArgs(C) || !IsWhole(C.N(3)) || C.N(3) < 3 || C.N(0) * C.N(3) > 2e6)
	{
		return Result::Fail("needs 1 to 4096 elements, spacing > 0 and a whole number of samples M >= 3 (N M up to 2e6)");
	}
	const int N = C.Int(0), Mm = C.Int(3);
	const double DLam = C.N(1), Steer = C.N(2);
	std::vector<double> S(static_cast<size_t>(Mm)), Af(static_cast<size_t>(Mm));
	for (int i = 0; i < Mm; ++i)
	{
		const double U = -1 + 2.0 * i / (Mm - 1);
		S[static_cast<size_t>(i)] = U;
		Af[static_cast<size_t>(i)] = ArrayFactor(N, DLam, Steer, U);
	}
	const double S0 = std::sin(Steer);
	size_t I0 = 0;
	for (size_t i = 0; i < S.size(); ++i)
	{
		if (std::fabs(S[i] - S0) < std::fabs(S[I0] - S0))
		{
			I0 = i;
		}
	}
	size_t Lo = I0, Hi = I0;
	while (Lo > 0 && std::pow(Af[Lo], 2) > 0.5)
	{
		--Lo;
	}
	while (Hi < S.size() - 1 && std::pow(Af[Hi], 2) > 0.5)
	{
		++Hi;
	}
	double Grating = 0;
	const double NullW = 1 / (N * DLam);
	for (size_t i = 0; i < S.size(); ++i)
	{
		if (std::fabs(S[i] - S0) > NullW)
		{
			Grating = std::max(Grating, Af[i]);
		}
	}
	return Result::Record({{"peak", std::asin(S[I0])}, {"width", std::asin(S[Hi]) - std::asin(S[Lo])}, {"grating", Grating}, {"limit", 1 / (1 + std::fabs(S0))},
		{"phase", 2 * cst::Pi * DLam * S0}});
}
const double AEarth = 4.0 / 3.0 * 6371e3;
Result RHorizon(const Call& C)
{
	const double H1 = C.N(0), H2 = C.N(1);
	if (!(H1 >= 0) || !(H2 >= 0))
	{
		return Result::Fail("heights must be >= 0");
	}
	return Result::Of(std::sqrt(std::pow(AEarth + H1, 2) - AEarth * AEarth) + std::sqrt(std::pow(AEarth + H2, 2) - AEarth * AEarth));
}
Result RHorizonRule(const Call& C)
{
	if (!(C.N(0) >= 0) || !(C.N(1) >= 0))
	{
		return Result::Fail("heights must be >= 0");
	}
	return Result::Of(4.1218 * (std::sqrt(C.N(0)) + std::sqrt(C.N(1))) * 1e3);
}
Result FrameTime(const Call& C)
{
	const double W = C.N(0), Beam = C.N(1), Bars = C.N(2), Dwell = C.N(3);
	if (!(W > 0) || !(Beam > 0) || !(Bars >= 1) || !(Dwell > 0))
	{
		return Result::Fail("needs width, beamwidth, dwell > 0 and bars >= 1");
	}
	// ceil(W / beam), forgiving the rounding of degrees typed as radians (rad(120)/rad(3) = 40.000...01).
	return Result::Of(std::ceil(W / Beam * (1 - 1e-12)) * Bars * Dwell);
}
void Triangulate(double B, double T1, double T2, double* X, double* Y)
{
	const double R1 = B * std::sin(T2) / std::sin(T2 - T1);
	*X = R1 * std::cos(T1);
	*Y = R1 * std::sin(T1);
}
Result TriBear(const Call& C)
{
	if (std::sin(C.N(2) - C.N(1)) == 0.0)
	{
		return Result::Fail("parallel bearings never cross");
	}
	double X = 0, Y = 0;
	Triangulate(C.N(0), C.N(1), C.N(2), &X, &Y);
	return Result::Record({{"x", X}, {"y", Y}, {"r", std::hypot(X, Y)}});
}
Result TriSpread(const Call& C)
{
	const double R0 = C.N(0), Base = C.N(1), Sig = C.N(2), Bearing = C.N(3), Nd = C.N(4), Seed = C.N(5);
	if (!(R0 > 0) || !(Base > 0) || !(Sig >= 0) || !IsWhole(Nd) || Nd < 2 || Nd > 4e5 || !IsWhole(Seed) || Seed < 0 || Seed > 4294967295.0)
	{
		return Result::Fail("needs R, baseline > 0, sigma >= 0, a whole n from 2 to 4e5 and a whole seed");
	}
	Mulberry R(Seed);
	const double Tx = R0 * std::cos(Bearing), Ty = R0 * std::sin(Bearing), T1 = std::atan2(Ty, Tx), T2 = std::atan2(Ty, Tx - Base);
	double S = 0, S2 = 0;
	const long N = static_cast<long>(Nd);
	for (long i = 0; i < N; ++i)
	{
		const double A1 = T1 + Sig * Gauss(R);
		const double A2 = T2 + Sig * Gauss(R);
		double X = 0, Y = 0;
		Triangulate(Base, A1, A2, &X, &Y);
		const double E = std::hypot(X, Y) - R0;
		S += E;
		S2 += E * E;
	}
	const double M = S / N;
	return Result::Of(std::sqrt(S2 / N - M * M));
}
Result RBand(const Call& C)
{
	static const double Bands[7][2] = {{1, 2}, {2, 4}, {4, 8}, {8, 12}, {12, 18}, {18, 27}, {27, 40}};
	const double FGHz = C.N(0) / 1e9;
	for (int i = 0; i < 7; ++i)
	{
		if (FGHz >= Bands[i][0] && FGHz < Bands[i][1])
		{
			return Result::Record({{"band", i + 1.0}, {"lo", Bands[i][0] * 1e9}, {"hi", Bands[i][1] * 1e9}, {"lam", cst::C / C.N(0)}});
		}
	}
	return Result::Record({{"band", 0.0}, {"lo", NaN}, {"hi", NaN}, {"lam", cst::C / C.N(0)}});
}

/* ---------- Luneburg lens ---------- */
double EpsAt(double R, double Emin, double Emax) { return std::min(Emax, std::max(Emin, 2 - R * R)); }
double TanEff(double E, double Tand, double Eh)
{
	const double F = (E - 1) / (Eh - 1);
	return F * Eh * Tand / E;
}
struct Ray
{
	double X, Y, Opl, Loss, Angle;
};
// RK4 on dx/dt = p, dp/dt = grad(n^2)/2 from x = -infinity at height b (lens radius 1), as trace().
Ray Trace(double B, double Emin, double Emax, double Dt, bool WithLoss, double Tand, double Eh)
{
	const double Xe = -std::sqrt(1 - B * B), Ye = B, N0 = std::sqrt(EpsAt(1, Emin, Emax));
	const double Nx = Xe, Ny = Ye, Cosi = -Nx, Eta = 1 / N0, K = 1 - Eta * Eta * (1 - Cosi * Cosi);
	const double Dx = Eta + (Eta * Cosi - std::sqrt(K)) * Nx, Dy = (Eta * Cosi - std::sqrt(K)) * Ny;
	double X = Xe + 1e-9 * Dx, Y = Ye + 1e-9 * Dy, Px = N0 * Dx, Py = N0 * Dy, Opl = 1 + Xe, Loss = 0;
	auto G = [Emin, Emax](double Gx, double Gy, double* Ax, double* Ay) {
		const double E = 2 - Gx * Gx - Gy * Gy;
		const bool Flat = E <= Emin || E >= Emax;
		*Ax = Flat ? 0.0 : -Gx;
		*Ay = Flat ? 0.0 : -Gy;
	};
	for (int i = 0; i < 1000000; ++i)
	{
		const double E = EpsAt(std::hypot(X, Y), Emin, Emax);
		Opl += E * Dt;
		if (WithLoss)
		{
			Loss += TanEff(E, Tand, Eh) * E * Dt;
		}
		double A1x, A1y, A2x, A2y, A3x, A3y, A4x, A4y;
		G(X, Y, &A1x, &A1y);
		const double K1[4] = {Px, Py, A1x, A1y};
		G(X + Dt / 2 * K1[0], Y + Dt / 2 * K1[1], &A2x, &A2y);
		const double K2[4] = {Px + Dt / 2 * K1[2], Py + Dt / 2 * K1[3], A2x, A2y};
		G(X + Dt / 2 * K2[0], Y + Dt / 2 * K2[1], &A3x, &A3y);
		const double K3[4] = {Px + Dt / 2 * K2[2], Py + Dt / 2 * K2[3], A3x, A3y};
		G(X + Dt * K3[0], Y + Dt * K3[1], &A4x, &A4y);
		const double K4[4] = {Px + Dt * K3[2], Py + Dt * K3[3], A4x, A4y};
		X += Dt / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
		Y += Dt / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
		Px += Dt / 6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
		Py += Dt / 6 * (K1[3] + 2 * K2[3] + 2 * K3[3] + K4[3]);
		if (X * X + Y * Y >= 1 && X > 0)
		{
			break;
		}
	}
	return {X, Y, Opl, Loss, std::atan2(Y, X)};
}
bool ProfileOk(double Emin, double Emax) { return Emin >= 1 && Emax > Emin && Emax <= 10; }
Result LuneEps(const Call& C)
{
	if (!ProfileOk(C.N(1), C.N(2)) || !(C.N(0) >= 0))
	{
		return Result::Fail("needs r >= 0 and 1 <= emin < emax");
	}
	return Result::Of(EpsAt(C.N(0), C.N(1), C.N(2)));
}
Result LuneTrace(const Call& C)
{
	const double B = C.N(0), Emin = C.N(1), Emax = C.N(2), Dt = C.N(3), Tand = C.N(4), Eh = C.N(5);
	if (!(B >= 0 && B < 1) || !ProfileOk(Emin, Emax) || !(Dt >= 1e-5 && Dt <= 0.1) || !(Tand >= 0) || !(Eh > 1))
	{
		return Result::Fail("needs 0 <= b < 1 (lens radii), 1 <= emin < emax, 1e-5 <= dt <= 0.1, tand >= 0, eh > 1");
	}
	const Ray Rr = Trace(B, Emin, Emax, Dt, Tand > 0, Tand, Eh);
	return Result::Record({{"x", Rr.X}, {"y", Rr.Y}, {"opl", Rr.Opl}, {"loss", Rr.Loss}, {"angle", Rr.Angle}});
}
Result LuneSpot(const Call& C)
{
	const double Emin = C.N(0), Emax = C.N(1), Nbd = C.N(2), Dt = C.N(3);
	if (!ProfileOk(Emin, Emax) || !IsWhole(Nbd) || Nbd < 2 || Nbd > 100 || !(Dt >= 1e-4 && Dt <= 0.1))
	{
		return Result::Fail("needs 1 <= emin < emax, a whole number of rays 2 to 100 and 1e-4 <= dt <= 0.1");
	}
	const int Nb = static_cast<int>(Nbd);
	double Sum2 = 0, OMin = std::numeric_limits<double>::infinity(), OMax = -OMin, Worst = -OMin;
	for (int i = 0; i < Nb; ++i)
	{
		const double B = 0.04 + 0.92 * i / (Nb - 1);
		const Ray Rr = Trace(B, Emin, Emax, Dt, false, 0, 2.8);
		Sum2 += Rr.Angle * Rr.Angle;
		OMin = std::min(OMin, Rr.Opl);
		OMax = std::max(OMax, Rr.Opl);
		Worst = std::max(Worst, std::hypot(Rr.X - 1, Rr.Y));
	}
	return Result::Record({{"rms", std::sqrt(Sum2 / Nb)}, {"opl", OMax - OMin}, {"worst", Worst}});
}
// Mixing rules: 1 linear (parallel), 2 Maxwell Garnett, 3 Lichtenecker, 4 series (lower bound).
double Mix(int Rule, double F, double Eh)
{
	switch (Rule)
	{
	case 1:
		return F * Eh + (1 - F);
	case 2:
		return (2 + Eh + 2 * F * (Eh - 1)) / (2 + Eh - F * (Eh - 1));
	case 3:
		return std::pow(Eh, F);
	default:
		return 1 / (F / Eh + (1 - F));
	}
}
bool RuleOk(double R) { return IsWhole(R) && R >= 1 && R <= 4; }
const char* const MixMessage = "rule: 1 linear, 2 Maxwell Garnett, 3 Lichtenecker, 4 series; eh > 1";
Result MixEps(const Call& C)
{
	if (!RuleOk(C.N(2)) || !(C.N(1) > 1) || !(C.N(0) >= 0 && C.N(0) <= 1))
	{
		return Result::Fail(!(C.N(0) >= 0 && C.N(0) <= 1) ? "fill fraction must be in [0, 1]" : MixMessage);
	}
	return Result::Of(Mix(C.Int(2), C.N(0), C.N(1)));
}
Result FillFor(const Call& C)
{
	const double Eps = C.N(0), Eh = C.N(1);
	if (!RuleOk(C.N(2)) || !(Eh > 1))
	{
		return Result::Fail(MixMessage);
	}
	if (!(Eps >= 1 && Eps <= Eh))
	{
		return Result::Fail("target permittivity must lie between 1 (air) and eh (solid)");
	}
	double Lo = 0, Hi = 1;
	for (int i = 0; i < 60; ++i)
	{
		const double M = (Lo + Hi) / 2;
		if (Mix(C.Int(2), M, Eh) < Eps)
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	return Result::Of((Lo + Hi) / 2);
}
// Wang, Lambert & Chisum (2023): cell sizes (mm) and the frequency (GHz) where gain falls away.
Result CellFit(const Call& C)
{
	static const double Data[3][2] = {{12.5, 20}, {10, 25}, {7.5, 33}};
	const double Cell = C.N(0);
	if (!(Cell > 0))
	{
		return Result::Fail("cell size must be > 0 (m)");
	}
	double Xs[3], Ys[3], Mx = 0, My = 0, Prod = 0;
	for (int i = 0; i < 3; ++i)
	{
		Xs[i] = std::log(Data[i][0]);
		Ys[i] = std::log(Data[i][1]);
	}
	Mx = (Xs[0] + Xs[1] + Xs[2]) / 3;
	My = (Ys[0] + Ys[1] + Ys[2]) / 3;
	double Sxy = 0, Sxx = 0;
	for (int i = 0; i < 3; ++i)
	{
		Sxy += (Xs[i] - Mx) * (Ys[i] - My);
		Sxx += std::pow(Xs[i] - Mx, 2);
		Prod += Data[i][0] * Data[i][1];
	}
	Prod /= 3;
	const double ProdSi = Prod * 1e6;   // mm GHz -> m Hz
	return Result::Record({{"slope", Sxy / Sxx}, {"prod", ProdSi}, {"cell_lam", ProdSi / cst::C}, {"fmax", Prod / (Cell * 1e3) * 1e9}});
}
Result TanEffFn(const Call& C)
{
	if (!(C.N(0) >= 1) || !(C.N(2) > 1) || !(C.N(1) >= 0))
	{
		return Result::Fail("needs eps >= 1, tand >= 0 and eh > 1");
	}
	return Result::Of(TanEff(C.N(0), C.N(1), C.N(2)));
}
double RadixDf(double F) { return F <= 10e9 ? 0.0043 : F >= 24e9 ? 0.0046 : 0.0043 + 0.0003 * (F - 10e9) / 14e9; }
Result RadixDfFn(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("frequency must be > 0");
	}
	return Result::Of(RadixDf(C.N(0)));
}
Result LensLossDb(const Call& C)
{
	const double A = C.N(0), F = C.N(1), Tand = C.N(2), Eh = C.N(3), Emin = C.N(4), Nbd = C.N(5);
	if (!(A > 0) || !(F > 0) || !(Tand >= 0) || !(Eh > 1) || !ProfileOk(Emin, 2) || !IsWhole(Nbd) || Nbd < 1 || Nbd > 200)
	{
		return Result::Fail("needs radius, f > 0, tand >= 0, eh > 1, 1 <= emin < 2 and 1 to 200 rays");
	}
	const int Nb = static_cast<int>(Nbd);
	const double K0 = 2 * cst::Pi * F / cst::C;
	double Trans = 0, Wsum = 0;
	for (int i = 0; i < Nb; ++i)
	{
		const double B = (i + 0.5) / Nb * 0.98;
		const Ray Rr = Trace(B, Emin, 2, 2e-3, true, Tand, Eh);
		Trans += B * std::exp(-K0 * A * Rr.Loss);
		Wsum += B;
	}
	return Result::Of(-10 * std::log10(Trans / Wsum));
}
Result AxLoss(const Call& C)
{
	const double Tand = C.N(0), Eh = C.N(1), Nd = C.N(2);
	if (!(Tand >= 0) || !(Eh > 1) || !IsWhole(Nd) || Nd < 1 || Nd > 5e6)
	{
		return Result::Fail("needs tand >= 0, eh > 1 and a whole n from 1 to 5e6");
	}
	const long N = static_cast<long>(Nd);
	double S = 0;
	for (long i = 0; i < N; ++i)
	{
		const double X = -1 + 2 * (i + 0.5) / N, E = 2 - X * X;
		S += std::sqrt(E) * TanEff(E, Tand, Eh) * 2 / N;
	}
	return Result::Of(S);
}
Result LensDirDb(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("radius and frequency must be > 0");
	}
	return Result::Of(10 * std::log10(4 * cst::Pi * cst::Pi * C.N(0) * C.N(0) / std::pow(cst::C / C.N(1), 2)));
}
double RcsLune(double A, double Lam) { return 4 * std::pow(cst::Pi, 3) * std::pow(A, 4) / (Lam * Lam); }
Result RcsLuneFn(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("radius and wavelength must be > 0");
	}
	return Result::Of(RcsLune(C.N(0), C.N(1)));
}
Result ReflRange(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0) || !(C.N(2) > 0))
	{
		return Result::Fail("radius, wavelength and the target's own cross section must be > 0");
	}
	const double S0 = C.N(2), Sig = RcsLune(C.N(0), C.N(1));
	return Result::Of(std::pow((S0 + Sig) / S0, 0.25));
}
} // namespace radar_detail

void RegisterRadarLibrary(Registry& R)
{
	using namespace radar_detail;
	// Bench defaults: the lecture's airport radar (radar/ section 1) and the search chart (section 2).
	const double G33 = Undb(33), L8 = Undb(8), Lam0 = cst::C / 2.998e9, A0 = cst::Pi * std::pow(10.0, 0.18) * std::pow(10.0, 0.18) / 4;
	const char* Track = "radar/ section 1, range equation";
	const char* Search = "radar/ section 2, search form";
	const char* Detect = "radar/ section 3, detection";
	const char* GertSrc = "radar/ section 5, GERT time to first track";

	Def(R, "snr_track", "radar", Track, "S/N per pulse, Pt G^2 lam^2 sigma / ((4 pi)^3 R^4 k Ts Bn L)", "power ratio (db() for dB)",
		{Num("Pt", "W", 1.4e6), Num("G", "-", G33), Num("lam", "m", Lam0), Num("sigma", "m^2", 1), Num("R", "m", 111e3), Num("Ts", "K", 950), Num("Bn", "Hz", 1.67e6),
			Num("L", "-", L8)},
		&SnrTrackFn);
	Def(R, "snr_chain", "radar", "radar/ section 1, range equation, chained steps", "the range equation step by step: density at the target, reflected, back, received",
		"record ae m^2, s_tgt W/m^2, refl W, s_back W/m^2, pr W, noise W, snr",
		{Num("Pt", "W", 1.4e6), Num("G", "-", G33), Num("lam", "m", Lam0), Num("sigma", "m^2", 1), Num("R", "m", 111e3), Num("Ts", "K", 950), Num("Bn", "Hz", 1.67e6),
			Num("L", "-", L8)},
		&SnrChain);
	Def(R, "snr_budget", "radar", "radar/ section 1, dB budget", "the lecture's dB budget: each term, the + and - columns and S/N",
		"record of signed dB contributions (pt, g2, lam2, sigma, fourpi3, r4, inv_k, ts, loss, bn; they sum to total), plus, minus, total dB",
		{Num("Pt", "W", 1.4e6), Num("G", "-", G33), Num("lam", "m", Lam0), Num("sigma", "m^2", 1), Num("R", "m", 111e3), Num("Ts", "K", 950), Num("Bn", "Hz", 1.67e6),
			Num("L", "-", L8)},
		&SnrBudget);
	Def(R, "range_for", "radar", Track, "range at which the track form gives a required S/N per pulse", "m",
		{Num("snr", "-", Undb(13.18)), Num("Pt", "W", 1.4e6), Num("G", "-", G33), Num("lam", "m", Lam0), Num("sigma", "m^2", 1), Num("Ts", "K", 950), Num("Bn", "Hz", 1.67e6),
			Num("L", "-", L8)},
		&RangeFor);
	Def(R, "snr_search", "radar", Search, "search form S/N, Pav A ts sigma / (4 pi Omega R^4 k Ts L)", "power ratio",
		{Num("Pav", "W", 1000), Num("A", "m^2", A0), Num("ts", "s", 10), Num("sigma", "m^2", 1), Num("Om", "sr", 1), Num("R", "m", 100e3), Num("Ts", "K", 500),
			Num("L", "-", 10)},
		&SnrSearch);
	Def(R, "pav_search", "radar", Search, "average power a search radar needs for S/N at range R", "W",
		{Num("R", "m", 100e3), Num("A", "m^2", A0), Num("ts", "s", 10), Num("sigma", "m^2", 1), Num("Om", "sr", 1), Num("Ts", "K", 500), Num("L", "-", 10),
			Num("snr", "-", Undb(15))},
		&PavSearch);
	Def(R, "range_search", "radar", Search, "search range at a required S/N", "m",
		{Num("Pav", "W", 1000), Num("A", "m^2", A0), Num("ts", "s", 10), Num("sigma", "m^2", 1), Num("Om", "sr", 1), Num("Ts", "K", 500), Num("L", "-", 10),
			Num("snr", "-", Undb(15))},
		&RangeSearch);
	Def(R, "search_track", "radar", "radar/ section 2, search via track", "the track form summed coherently over the pulses in each beam position (check of the search form)",
		"record snr, n pulses per beam, beams",
		{Num("Pav", "W", 1000), Num("A", "m^2", A0), Num("ts", "s", 10), Num("sigma", "m^2", 1), Num("Om", "sr", 1), Num("R", "m", 100e3), Num("Ts", "K", 500),
			Num("L", "-", 10), Num("lam", "m", 0.1), Num("tau", "s", 1e-6), Num("prf", "Hz", 1000)},
		&SearchTrack);
	Def(R, "tsys", "radar", "radar/ section 4, system noise temperature", "Ts = Ta + (Lr - 1) Tline + Lr (F - 1) 290 K", "record tr, te, ts K",
		{Num("Ta", "K", 150), Num("Lr", "-", Undb(1)), Num("Tline", "K", 290), Num("F", "-", Undb(3))}, &Tsys);
	Def(R, "tsys_trace", "radar", "radar/ section 4, system noise temperature, traced", "noise traced through the line and receiver, referred back to the antenna", "K",
		{Num("Ta", "K", 150), Num("Lr", "-", Undb(1)), Num("Tline", "K", 290), Num("F", "-", Undb(3))}, &TsysTrace);
	Def(R, "pfa_thr", "radar", Detect, "false-alarm probability of a square-law threshold over n pulses, Q(n, T)", "probability",
		{Num("T", "noise power"), Num("n", "pulses", 1)}, &PfaThr);
	Def(R, "det_thresh", "radar", Detect, "detection threshold for a false-alarm probability (bisection, as the bench)", "threshold, in units of the noise power",
		{Num("n", "pulses", 1), Num("pfa", "-", 1e-6)}, &DetThresh);
	Def(R, "pd_steady", "radar", Detect, "Pd of a steady target, n pulses integrated noncoherently (Marcum)", "probability",
		{Num("snr", "-"), Num("n", "pulses", 1), Num("pfa", "-", 1e-6)}, &PdSteadyFn);
	Def(R, "pd_sw1", "radar", Detect, "Pd of a Swerling 1 target, single pulse, pfa^(1/(1+snr))", "probability", {Num("snr", "-"), Num("pfa", "-", 1e-6)}, &PdSw1);
	Def(R, "snr_req_db", "radar", Detect, "S/N per pulse a steady target needs for pd with n pulses", "dB",
		{Num("pd", "-", 0.9), Num("n", "pulses", 1), Num("pfa", "-", 1e-6)}, &SnrReqDb);
	Def(R, "albersheim", "radar", Detect, "Albersheim's closed-form S/N per pulse for pd, pfa, n", "dB", {Num("pd", "-", 0.9), Num("pfa", "-", 1e-6), Num("n", "pulses", 1)},
		&Albersheim);
	Def(R, "intgain", "radar", Detect, "noncoherent integration gain of n pulses against coherent 10 log10 n", "record nc, coh, loss dB",
		{Num("n", "pulses", 21), Num("pd", "-", 0.9), Num("pfa", "-", 1e-6)}, &IntGain);
	Def(R, "marcumq", "radar", Detect, "Marcum Q_m(a, b), the series behind the steady-target Pd", "probability", {Num("a", "-"), Num("b", "-"), Num("m", "-", 1)},
		&MarcumQ);
	Def(R, "pd_mc", "radar", "radar/ section 3, detection, Monte Carlo", "Monte Carlo Pd of the same detector (mulberry32, as the bench)", "probability",
		{Num("snr", "-"), Num("n", "pulses", 1), Num("pfa", "-", 1e-6), Num("trials", "-", 20000), Num("sw1", "0/1", 0), Num("seed", "-", 7)}, &PdMc);
	Def(R, "gert", "radar", GertSrc, "GERT network [from, to, p, t] rows (nodes from 1): probability, mean and variance of completion",
		"record p, mean s, var s^2, sd s", {Mat("B", "k x 4"), Num("start", "node", 1), Num("finish", "node (0 = last)", 0)}, &Gert);
	Def(R, "ttrack", "radar", GertSrc, "time to first track by GERT: rule 1 = 2 in a row, 2 = 2 of 3, 3 = 3 in a row",
		"record p, mean s, var s^2, sd s, p90 scans", {Num("pd", "-", 0.7), Num("T", "s", 4.69), Num("rule", "1/2/3", 1)}, &TTrack);
	Def(R, "trk_net", "radar", GertSrc, "branch table of a track-confirmation rule, ready for gert()", "k x 4 [from, to, p, t]",
		{Num("pd", "-", 0.7), Num("T", "s", 4.69), Num("rule", "1/2/3", 1)}, &TrkNet);

	const char* Pulse = "radar/basics Pulses and ranges";
	const char* Dop = "radar/basics Doppler";
	const char* Beams = "radar/basics Beams and arrays";
	const char* Geo = "radar/basics Geometry and timing";
	Def(R, "mfrange", "signal", "radar/basics Range from echo delay", "matched-filter range of a simulated 1 us echo in noise (mulberry32)",
		"record est m, sample m", {Num("R", "m"), Num("fs", "Hz", 20e6), Num("tau", "s", 1e-6), Num("noise", "-", 0.05), Num("seed", "-", 14)}, &MfRange);
	Def(R, "chirp_res", "signal", "radar/basics Range resolution", "range resolution c / 2B", "m", {Num("B", "Hz")}, &ChirpRes);
	Def(R, "chirp2", "signal", "radar/basics Range resolution, pulse compression", "two echoes sep apart through a chirp's matched filter: resolved?",
		"record resolved 0/1, res m, ratio BT, mps m per sample", {Num("sep", "m"), Num("B", "Hz", 5e6), Num("T", "s", 10e-6), Num("fs", "Hz", 50e6)}, &Chirp2);
	Def(R, "r_unamb", "signal", Pulse, "maximum unambiguous range c / 2 PRF", "m", {Num("prf", "Hz")}, &RUnamb);
	Def(R, "v_unamb", "signal", Dop, "unambiguous radial speed lam PRF / 4", "m/s", {Num("prf", "Hz"), Num("lam", "m")}, &VUnamb);
	Def(R, "rfold", "signal", Pulse, "range shown for a target beyond R_u: R mod c / 2 PRF", "m", {Num("R", "m"), Num("prf", "Hz")}, &RFold);
	Def(R, "valias", "signal", Dop, "radial speed folded into +-lam PRF / 4", "m/s", {Num("v", "m/s"), Num("prf", "Hz"), Num("lam", "m")}, &VAlias);
	Def(R, "dop_shift", "signal", Dop, "two-way Doppler shift 2 v / lam", "Hz", {Num("v", "m/s"), Num("lam", "m")}, &DopShift);
	Def(R, "dop_dft", "signal", Dop, "speed read by a Hann-windowed N-pulse DFT of a simulated echo (mulberry32)", "record est m/s, bin m/s",
		{Num("v", "m/s"), Num("prf", "Hz"), Num("lam", "m"), Num("N", "pulses", 128), Num("seed", "-", 5)}, &DopDft);
	Def(R, "prf_trade", "signal", "radar/basics The range-velocity trade", "R_u and v_u at one PRF; their product is c lam / 8",
		"record ru m, vu m/s, prod m^2/s, regime 1 low / 2 medium / 3 high PRF", {Num("prf", "Hz"), Num("lam", "m")}, &PrfTrade);
	Def(R, "mti", "signal", "radar/basics MTI", "two-pulse canceller: simulated gain, 4 sin^2(pi fd / PRF), first blind speed",
		"record meas, pred, gain_db dB vs average, vblind m/s", {Num("v", "m/s"), Num("prf", "Hz", 1000), Num("lam", "m", 0.03), Num("N", "pulses", 2000)}, &Mti);
	Def(R, "clutter", "signal", "radar/basics Ground clutter", "Doppler band of ground clutter seen from a radar moving at V (simulated patches)",
		"record lo Hz, hi Hz, edge 2V/lam Hz", {Num("V", "m/s"), Num("lam", "m"), Num("n", "-", 20000), Num("seed", "-", 8)}, &Clutter);
	Def(R, "beam3db", "signal", Beams, "half-power beamwidth of a uniform aperture: solved sinc^2 = 1/2 and 0.886 lam / D", "record num rad, rule rad",
		{Num("D", "m"), Num("lam", "m")}, &Beam3db);
	Def(R, "afactor", "signal", Beams, "normalized array factor of N elements steered to theta0, at direction sine u", "0..1",
		{Num("N", "elements"), Num("d", "wavelengths"), Num("steer", "rad"), Num("u", "sin")}, &AFactor);
	Def(R, "ascan", "signal", Beams, "scanned array: beam peak, half-power width, largest other lobe, grating-lobe spacing limit, phase step",
		"record peak rad, width rad, grating, limit wavelengths, phase rad", {Num("N", "elements"), Num("d", "wavelengths"), Num("steer", "rad"), Num("M", "samples", 8001)},
		&AScan);
	Def(R, "r_horizon", "signal", Geo, "radar horizon over a 4/3 Earth, exact tangents", "m", {Num("h1", "m"), Num("h2", "m", 0)}, &RHorizon);
	Def(R, "r_horizon_rule", "signal", Geo, "radar horizon rule of thumb 4.12 (sqrt h1 + sqrt h2) km", "m", {Num("h1", "m"), Num("h2", "m", 0)}, &RHorizonRule);
	Def(R, "frame_time", "signal", Geo, "scan frame time ceil(W / beam) x bars x dwell", "s",
		{Num("W", "rad"), Num("beam", "rad"), Num("bars", "-", 1), Num("dwell", "s")}, &FrameTime);
	Def(R, "tri_bear", "signal", "radar/basics Passive angle-only ranging", "fix from two bearings a baseline b apart (sensors at the origin and at (b, 0))",
		"record x m, y m, r m", {Num("b", "m"), Num("t1", "rad"), Num("t2", "rad")}, &TriBear);
	Def(R, "tri_spread", "signal", "radar/basics Passive angle-only ranging", "simulated range error of a two-bearing fix (mulberry32)", "m",
		{Num("R", "m"), Num("b", "m"), Num("sigma", "rad"), Num("bearing", "rad", 1.2), Num("n", "-", 4000), Num("seed", "-", 9)}, &TriSpread);
	Def(R, "rband", "signal", "radar/basics Bands and wavelength", "IEEE radar band: 1 L, 2 S, 3 C, 4 X, 5 Ku, 6 K, 7 Ka (0 outside)", "record band, lo Hz, hi Hz, lam m",
		{Num("f", "Hz")}, &RBand);

	const char* Rays = "lens/ section 1, rays";
	const char* Print = "lens/ section 2, printing the gradient";
	const char* Loss = "lens/ section 4, loss and gain";
	Def(R, "lune_eps", "lens", Rays, "Luneburg permittivity 2 - r^2, clipped to the printable range", "relative permittivity",
		{Num("r", "lens radii"), Num("emin", "-", 1), Num("emax", "-", 2)}, &LuneEps);
	Def(R, "lune_trace", "lens", Rays, "ray at height b through the lens (RK4): exit point, optical path, loss integral",
		"record x, y, opl (radii), loss (radii), angle rad", {Num("b", "lens radii"), Num("emin", "-", 1), Num("emax", "-", 2), Num("dt", "-", 1e-3), Num("tand", "-", 0),
			Num("eh", "-", 2.8)},
		&LuneTrace);
	Def(R, "lune_spot", "lens", Rays, "spot of nb parallel rays: rms exit angle, optical-path spread, worst miss", "record rms rad, opl radii, worst radii",
		{Num("emin", "-", 1.2), Num("emax", "-", 2.8), Num("nb", "rays", 24), Num("dt", "-", 1e-3)}, &LuneSpot);
	Def(R, "mix_eps", "lens", Print, "effective permittivity of a printed lattice: rule 1 linear, 2 Maxwell Garnett, 3 Lichtenecker, 4 series",
		"relative permittivity", {Num("f", "fill fraction"), Num("eh", "-", 2.8), Num("rule", "1-4", 2)}, &MixEps);
	Def(R, "fill_for", "lens", Print, "solid fill fraction that gives permittivity eps under a mixing rule", "fill fraction",
		{Num("eps", "-"), Num("eh", "-", 2.8), Num("rule", "1-4", 2)}, &FillFor);
	Def(R, "cell_fit", "lens", "lens/ section 3, cell size", "fit to the measured cell-size limits (Wang, Lambert & Chisum 2023) and the limit for one cell",
		"record slope, prod m Hz, cell_lam, fmax Hz", {Num("cell", "m", 7.5e-3)}, &CellFit);
	Def(R, "tan_eff", "lens", Loss, "effective loss tangent of a lattice at permittivity eps (linear mixing)", "-",
		{Num("eps", "-"), Num("tand", "-"), Num("eh", "-", 2.8)}, &TanEffFn);
	Def(R, "radix_df", "lens", Loss, "Rogers Radix loss tangent: 0.0043 to 10 GHz, 0.0046 from 24 GHz, linear between", "-", {Num("f", "Hz")}, &RadixDfFn);
	Def(R, "lens_loss_db", "lens", Loss, "aperture-weighted dielectric loss of a Luneburg lens (ray traced)", "dB",
		{Num("a", "m"), Num("f", "Hz"), Num("tand", "-"), Num("eh", "-", 2.8), Num("emin", "-", 1), Num("nb", "rays", 16)}, &LensLossDb);
	Def(R, "ax_loss", "lens", Loss, "loss integral along the axis by quadrature (check of the tracer)", "lens radii",
		{Num("tand", "-"), Num("eh", "-", 2.8), Num("n", "-", 20000)}, &AxLoss);
	Def(R, "lens_dir_db", "lens", Loss, "ideal directivity of a circular aperture 4 pi A / lam^2", "dBi", {Num("a", "m"), Num("f", "Hz")}, &LensDirDb);
	Def(R, "rcs_lune", "lens", "lens/ section 5, Luneburg reflector", "radar cross section of a Luneburg reflector 4 pi^3 a^4 / lam^2", "m^2", {Num("a", "m"), Num("lam", "m")},
		&RcsLuneFn);
	Def(R, "refl_range", "lens", "lens/ section 5, Luneburg reflector", "detection range multiplier from adding a reflector, ((s0 + sigma) / s0)^(1/4)", "-",
		{Num("a", "m"), Num("lam", "m"), Num("s0", "m^2", 0.1)}, &ReflRange);
}
} // namespace mcalc
