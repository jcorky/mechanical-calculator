// Chapters 7, 10, 11 and 12 of the atlas: lattices and GNSS ranging, category theory (kernels, transfer
// functions, sheaves), procedural systems (automata, IFS, Turing, splines, Lloyd) and expansion fields
// (Cramer-Rao, HJB, Fokker-Planck, Ito, tropical algebra, transport, symplectic steps, random matrices).
#include "CalcLibrary.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace mcalc
{
namespace fields_detail
{
using FieldsComplex = std::complex<double>;

constexpr double FieldsNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double FieldsInf = std::numeric_limits<double>::infinity();

Param MatDefault(const char* Name, const char* Unit, double Default)
{
	Param P = Mat(Name, Unit);
	P.HasDefault = true;
	P.Default = Default;
	return P;
}

double NormalCdf(double X) { return 0.5 * std::erfc(-X / std::sqrt(2.0)); }

// Phi(B) - Phi(A) without cancellation in either tail.
double NormalBand(double A, double B)
{
	if (A > 0.0)
	{
		return 0.5 * (std::erfc(A / std::sqrt(2.0)) - std::erfc(B / std::sqrt(2.0)));
	}
	return 0.5 * (std::erfc(-B / std::sqrt(2.0)) - std::erfc(-A / std::sqrt(2.0)));
}

bool IsWhole(double X, double Lo, double Hi) { return std::isfinite(X) && X == std::floor(X) && X >= Lo && X <= Hi; }

std::string ShapeText(const Matrix& M) { return std::to_string(M.Rows) + " x " + std::to_string(M.Cols); }

double RowDot(const Matrix& A, int I, const Matrix& B, int J)
{
	double S = 0.0;
	for (int c = 0; c < A.Cols; ++c)
	{
		S += A(I, c) * B(J, c);
	}
	return S;
}

// mulberry32, as in the atlas pages, so seeded benches reproduce bit for bit.
struct FieldsMulberry
{
	uint32_t State;
	explicit FieldsMulberry(uint32_t Seed) : State(Seed) {}
	double Next()
	{
		State += 0x6D2B79F5u;
		uint32_t T = State;
		T = (T ^ (T >> 15)) * (T | 1u);
		T ^= T + (T ^ (T >> 7)) * (T | 61u);
		return static_cast<double>(T ^ (T >> 14)) / 4294967296.0;
	}
};

// ---------------------------------------------------------------- ch07: time of flight

Result TofRange(const Call& C)
{
	if (!(C.N(1) > 0.0))
	{
		return Result::Fail("way must be positive: 1 for one-way, 2 for two-way (radar)");
	}
	return Result::Of(cst::C * C.N(0) / C.N(1));
}

Result TofRes(const Call& C)
{
	if (!(C.N(0) > 0.0))
	{
		return Result::Fail("bandwidth must be positive");
	}
	return Result::Of(cst::C / (2.0 * C.N(0)));
}

Result Widelane(const Call& C)
{
	const double F1 = C.N(0), F2 = C.N(1);
	if (!(F1 > F2 && F2 > 0.0))
	{
		return Result::Fail("needs f1 > f2 > 0");
	}
	return Result::Record({{"l1", cst::C / F1}, {"l2", cst::C / F2}, {"wide", cst::C / (F1 - F2)}, {"narrow", cst::C / (F1 + F2)}});
}

Result SatRange(const Call& C)
{
	const double El = C.N(0), Re = C.N(1), Rs = C.N(2);
	if (!(Re > 0.0 && Rs > Re))
	{
		return Result::Fail("needs 0 < re < rs");
	}
	const double D = std::sqrt(Rs * Rs - std::pow(Re * std::cos(El), 2)) - Re * std::sin(El);
	return Result::Record({{"d", D}, {"t", D / cst::C}});
}

// ---------------------------------------------------------------- ch07: pseudorange positioning

bool ReadSatellites(const Call& C, Matrix* S, std::vector<double>* Rho, std::string* Error)
{
	*S = C.M(0);
	if (S->Cols != 3 || S->Rows < 4)
	{
		*Error = "satellites are an n x 3 matrix of positions (m), n >= 4; got " + ShapeText(*S);
		return false;
	}
	if (Rho)
	{
		const Matrix& R = C.M(1);
		if (!R.IsVector() || R.Count() != S->Rows)
		{
			*Error = "pseudoranges must be a vector with one entry per satellite (" + std::to_string(S->Rows) + ")";
			return false;
		}
		*Rho = R.V;
	}
	return true;
}

// Line-of-sight geometry rows [-(s - x)/|s - x|, 1]; false if a satellite sits on the receiver.
bool GeometryAt(const Matrix& S, const double* X, Matrix* G, std::vector<double>* Dist)
{
	*G = Matrix(S.Rows, 4);
	Dist->assign(static_cast<size_t>(S.Rows), 0.0);
	for (int i = 0; i < S.Rows; ++i)
	{
		const double Dx = S(i, 0) - X[0], Dy = S(i, 1) - X[1], Dz = S(i, 2) - X[2];
		const double D = std::sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
		if (!(D > 0.0))
		{
			return false;
		}
		(*G)(i, 0) = -Dx / D;
		(*G)(i, 1) = -Dy / D;
		(*G)(i, 2) = -Dz / D;
		(*G)(i, 3) = 1.0;
		(*Dist)[static_cast<size_t>(i)] = D;
	}
	return true;
}

// Gauss-Newton on rho_i = |s_i - x| + c b from the Earth's centre, as in the chapter's bench.
bool PseudorangeFix(const Matrix& S, const std::vector<double>& Rho, double* X, int* Iters, std::string* Error)
{
	for (int a = 0; a < 4; ++a)
	{
		X[a] = 0.0;
	}
	Matrix G;
	std::vector<double> Dist;
	int It = 0;
	for (It = 0; It < 15; ++It)
	{
		if (!GeometryAt(S, X, &G, &Dist))
		{
			*Error = "a satellite coincides with the estimate";
			return false;
		}
		Matrix Res(S.Rows, 1);
		for (int i = 0; i < S.Rows; ++i)
		{
			Res.V[static_cast<size_t>(i)] = Rho[static_cast<size_t>(i)] - Dist[static_cast<size_t>(i)] - X[3];
		}
		const Matrix Gt = la::Transpose(G);
		Matrix Dx;
		std::string E;
		if (!la::Solve(la::Mul(Gt, G), la::Mul(Gt, Res), &Dx, &E))
		{
			*Error = "geometry is singular (satellites coplanar with the receiver, or fewer than 4 independent): " + E;
			return false;
		}
		double Step = 0.0;
		for (int a = 0; a < 4; ++a)
		{
			X[a] += Dx.V[static_cast<size_t>(a)];
			Step += Dx.V[static_cast<size_t>(a)] * Dx.V[static_cast<size_t>(a)];
		}
		if (std::sqrt(Step) < 1e-6)
		{
			break;
		}
	}
	*Iters = std::min(It + 1, 15);
	return true;
}

Result GnssFix(const Call& C)
{
	Matrix S;
	std::vector<double> Rho;
	std::string E;
	if (!ReadSatellites(C, &S, &Rho, &E))
	{
		return Result::Fail(E);
	}
	double X[4];
	int Iters = 0;
	if (!PseudorangeFix(S, Rho, X, &Iters, &E))
	{
		return Result::Fail(E);
	}
	return Result::Record({{"x", X[0]}, {"y", X[1]}, {"z", X[2]}, {"cb", X[3]}, {"dt", X[3] / cst::C}, {"it", static_cast<double>(Iters)}});
}

// DOPs from the geometry matrix; Up is the local vertical (unit) used to split vertical from horizontal.
Result DopRecord(const Matrix& G, const double* Up)
{
	Matrix Q;
	std::string E;
	if (!la::Inverse(la::Mul(la::Transpose(G), G), &Q, &E))
	{
		return Result::Fail("geometry is singular (3 satellites cannot separate clock and position): " + E);
	}
	const double P2 = Q(0, 0) + Q(1, 1) + Q(2, 2);
	double V2 = 0.0;
	for (int a = 0; a < 3; ++a)
	{
		for (int b = 0; b < 3; ++b)
		{
			V2 += Up[a] * Q(a, b) * Up[b];
		}
	}
	return Result::Record({{"gdop", std::sqrt(P2 + Q(3, 3))}, {"pdop", std::sqrt(P2)}, {"hdop", std::sqrt(std::max(0.0, P2 - V2))}, {"vdop", std::sqrt(V2)},
		{"tdop", std::sqrt(Q(3, 3))}});
}

Result Dop(const Call& C)
{
	Matrix S;
	std::string E;
	if (!ReadSatellites(C, &S, nullptr, &E))
	{
		return Result::Fail(E);
	}
	const Matrix& P = C.M(1);
	if (!P.IsVector() || P.Count() != 3)
	{
		return Result::Fail("receiver position must be a 3-vector (m)");
	}
	const double R = std::sqrt(P.V[0] * P.V[0] + P.V[1] * P.V[1] + P.V[2] * P.V[2]);
	if (!(R > 0.0))
	{
		return Result::Fail("receiver at the origin has no local vertical");
	}
	const double Up[3] = {P.V[0] / R, P.V[1] / R, P.V[2] / R};
	Matrix G;
	std::vector<double> Dist;
	if (!GeometryAt(S, P.V.data(), &G, &Dist))
	{
		return Result::Fail("a satellite coincides with the receiver");
	}
	return DopRecord(G, Up);
}

Result DopAe(const Call& C)
{
	const Matrix &Az = C.M(0), &El = C.M(1);
	if (!Az.IsVector() || !El.IsVector() || Az.Count() != El.Count() || Az.Count() < 4)
	{
		return Result::Fail("azimuth and elevation must be vectors of equal length >= 4 (rad)");
	}
	Matrix G(Az.Count(), 4);
	for (int i = 0; i < Az.Count(); ++i)
	{
		const double A = Az.V[static_cast<size_t>(i)], L = El.V[static_cast<size_t>(i)];
		G(i, 0) = -std::cos(L) * std::sin(A);
		G(i, 1) = -std::cos(L) * std::cos(A);
		G(i, 2) = -std::sin(L);
		G(i, 3) = 1.0;
	}
	const double Up[3] = {0.0, 0.0, 1.0};
	return DopRecord(G, Up);
}

// Baarda's w-statistic: residual over sigma sqrt(1 - h_ii), h the hat matrix of the fix.
Result Baarda(const Call& C)
{
	Matrix S;
	std::vector<double> Rho;
	std::string E;
	if (!ReadSatellites(C, &S, &Rho, &E))
	{
		return Result::Fail(E);
	}
	const double Sigma = C.N(2);
	if (!(Sigma > 0.0))
	{
		return Result::Fail("sigma must be positive (m)");
	}
	double X[4];
	int Iters = 0;
	if (!PseudorangeFix(S, Rho, X, &Iters, &E))
	{
		return Result::Fail(E);
	}
	Matrix G, Q;
	std::vector<double> Dist;
	if (!GeometryAt(S, X, &G, &Dist) || !la::Inverse(la::Mul(la::Transpose(G), G), &Q, &E))
	{
		return Result::Fail("geometry is singular");
	}
	Matrix W(S.Rows, 1);
	for (int i = 0; i < S.Rows; ++i)
	{
		double H = 0.0;
		for (int a = 0; a < 4; ++a)
		{
			for (int b = 0; b < 4; ++b)
			{
				H += G(i, a) * Q(a, b) * G(i, b);
			}
		}
		const double Res = Rho[static_cast<size_t>(i)] - Dist[static_cast<size_t>(i)] - X[3];
		W.V[static_cast<size_t>(i)] = Res / (Sigma * std::sqrt(std::max(1e-12, 1.0 - H)));
	}
	return Result::Of(W);
}

// ---------------------------------------------------------------- ch07: lattices

Result GaussReduce(const Call& C)
{
	const Matrix& B = C.M(0);
	if (B.Rows != 2 || B.Cols != 2)
	{
		return Result::Fail("basis must be 2 x 2, one basis vector per row");
	}
	double B1[2] = {B(0, 0), B(0, 1)}, B2[2] = {B(1, 0), B(1, 1)};
	if (B1[0] * B2[1] - B1[1] * B2[0] == 0.0)
	{
		return Result::Fail("basis vectors are parallel");
	}
	for (int Step = 0; Step < 10000; ++Step)
	{
		const double N1 = B1[0] * B1[0] + B1[1] * B1[1], N2 = B2[0] * B2[0] + B2[1] * B2[1];
		if (N2 < N1)
		{
			std::swap(B1[0], B2[0]);
			std::swap(B1[1], B2[1]);
			continue;
		}
		const double Q = std::nearbyint((B1[0] * B2[0] + B1[1] * B2[1]) / N1);
		if (Q == 0.0)
		{
			return Result::Of(Matrix::FromRows({{B1[0], B1[1]}, {B2[0], B2[1]}}));
		}
		B2[0] -= Q * B1[0];
		B2[1] -= Q * B1[1];
	}
	return Result::Fail("reduction did not terminate");
}

void GramSchmidt(const Matrix& B, Matrix* Bs, Matrix* Mu)
{
	const int N = B.Rows;
	*Bs = B;
	*Mu = Matrix(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			const double M = RowDot(B, i, *Bs, j) / RowDot(*Bs, j, *Bs, j);
			(*Mu)(i, j) = M;
			for (int c = 0; c < B.Cols; ++c)
			{
				(*Bs)(i, c) -= M * (*Bs)(j, c);
			}
		}
	}
}

bool LllReduce(Matrix* B, double Delta, std::string* Error)
{
	const int N = B->Rows;
	if (N == 0 || la::Rank(*B) < N)
	{
		*Error = "basis vectors (rows) must be linearly independent";
		return false;
	}
	Matrix Bs, Mu;
	int K = 1, Guard = 0;
	while (K < N)
	{
		if (++Guard > 200000)
		{
			*Error = "LLL did not converge";
			return false;
		}
		for (int j = K - 1; j >= 0; --j)
		{
			GramSchmidt(*B, &Bs, &Mu);
			const double Q = std::nearbyint(Mu(K, j));
			if (Q != 0.0)
			{
				for (int c = 0; c < B->Cols; ++c)
				{
					(*B)(K, c) -= Q * (*B)(j, c);
				}
			}
		}
		GramSchmidt(*B, &Bs, &Mu);
		if (RowDot(Bs, K, Bs, K) >= (Delta - Mu(K, K - 1) * Mu(K, K - 1)) * RowDot(Bs, K - 1, Bs, K - 1))
		{
			++K;
		}
		else
		{
			for (int c = 0; c < B->Cols; ++c)
			{
				std::swap((*B)(K, c), (*B)(K - 1, c));
			}
			K = std::max(K - 1, 1);
		}
	}
	return true;
}

Result Lll(const Call& C)
{
	Matrix B = C.M(0);
	const double Delta = C.N(1);
	if (!(Delta > 0.25 && Delta <= 1.0))
	{
		return Result::Fail("delta must be in (1/4, 1]");
	}
	if (B.Rows < 1 || B.Rows > B.Cols || B.Rows > 12)
	{
		return Result::Fail("basis is n x m with one vector per row, n <= m, n <= 12; got " + ShapeText(B));
	}
	std::string E;
	if (!LllReduce(&B, Delta, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(B);
}

struct FieldsEnum
{
	int N = 0;
	Matrix Mu;
	std::vector<double> Norms;
	std::vector<double> Coef;
	double Best = 0.0;
	double Radius = 0.0;
	double Tiny = 0.0;
	long Visits = 0;
};

// Fincke-Pohst depth-first enumeration in Gram-Schmidt coordinates.
void EnumerateShort(FieldsEnum& St, int I, double Partial)
{
	if (++St.Visits > 4000000)
	{
		return;
	}
	if (I < 0)
	{
		if (Partial > St.Tiny && Partial < St.Best)
		{
			St.Best = Partial;
			St.Radius = Partial * (1.0 + 1e-9);
		}
		return;
	}
	double Centre = 0.0;
	for (int j = I + 1; j < St.N; ++j)
	{
		Centre -= St.Coef[static_cast<size_t>(j)] * St.Mu(j, I);
	}
	const double R = std::sqrt(std::max(St.Radius - Partial, 0.0) / St.Norms[static_cast<size_t>(I)]);
	for (double Z = std::ceil(Centre - R); Z <= std::floor(Centre + R); Z += 1.0)
	{
		St.Coef[static_cast<size_t>(I)] = Z;
		const double P = Partial + (Z - Centre) * (Z - Centre) * St.Norms[static_cast<size_t>(I)];
		if (P < St.Radius)
		{
			EnumerateShort(St, I - 1, P);
		}
	}
	St.Coef[static_cast<size_t>(I)] = 0.0;
}

Result ShortVec(const Call& C)
{
	Matrix B = C.M(0);
	if (B.Rows < 1 || B.Rows > B.Cols || B.Rows > 10)
	{
		return Result::Fail("basis is n x m with one vector per row, n <= m, n <= 10; got " + ShapeText(B));
	}
	std::string E;
	if (!LllReduce(&B, 0.75, &E))
	{
		return Result::Fail(E);
	}
	FieldsEnum St;
	St.N = B.Rows;
	Matrix Bs;
	GramSchmidt(B, &Bs, &St.Mu);
	St.Best = FieldsInf;
	for (int i = 0; i < B.Rows; ++i)
	{
		St.Norms.push_back(RowDot(Bs, i, Bs, i));
		St.Best = std::min(St.Best, RowDot(B, i, B, i));
	}
	St.Radius = St.Best * (1.0 + 1e-9);
	St.Tiny = St.Best * 1e-12;
	St.Coef.assign(static_cast<size_t>(St.N), 0.0);
	EnumerateShort(St, St.N - 1, 0.0);
	return Result::Of(std::sqrt(St.Best));
}

Result LllBound(const Call& C)
{
	const double N = C.N(0), Delta = C.N(1);
	if (!IsWhole(N, 1, 1000) || !(Delta > 0.25 && Delta <= 1.0))
	{
		return Result::Fail("needs a whole dimension n >= 1 and delta in (1/4, 1]");
	}
	return Result::Of(std::pow(1.0 / (Delta - 0.25), (N - 1.0) / 2.0));
}

Result ODefect(const Call& C)
{
	const Matrix& B = C.M(0);
	if (B.Rows < 1 || B.Rows > B.Cols)
	{
		return Result::Fail("basis is n x m with one vector per row, n <= m");
	}
	double Det = 0.0;
	std::string E;
	if (!la::Determinant(la::Mul(B, la::Transpose(B)), &Det, &E) || !(Det > 0.0))
	{
		return Result::Fail("basis vectors (rows) must be linearly independent");
	}
	double P = 1.0;
	for (int i = 0; i < B.Rows; ++i)
	{
		P *= std::sqrt(RowDot(B, i, B, i));
	}
	return Result::Of(P / std::sqrt(Det));
}

// ---------------------------------------------------------------- ch07: integer ambiguities

Result EllCov(const Call& C)
{
	const double A = C.N(0), B = C.N(1), T = C.N(2), Cs = std::cos(T), Sn = std::sin(T);
	if (!(A > 0.0 && B > 0.0))
	{
		return Result::Fail("semi-axes must be positive (cycles)");
	}
	const double Off = Cs * Sn * (A * A - B * B);
	return Result::Of(Matrix::FromRows({{Cs * Cs * A * A + Sn * Sn * B * B, Off}, {Off, Sn * Sn * A * A + Cs * Cs * B * B}}));
}

bool ReadCovariance(const Matrix& Q, int MaxN, std::string* Error)
{
	if (!Q.IsSquare() || Q.Rows > MaxN)
	{
		*Error = "covariance must be square, at most " + std::to_string(MaxN) + " x " + std::to_string(MaxN) + "; got " + ShapeText(Q);
		return false;
	}
	for (int i = 0; i < Q.Rows; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			if (std::fabs(Q(i, j) - Q(j, i)) > 1e-12 * (std::fabs(Q(i, i)) + std::fabs(Q(j, j))))
			{
				*Error = "covariance must be symmetric";
				return false;
			}
		}
	}
	return true;
}

// Integer decorrelation (2-D Lagrange reduction of Q): unimodular Z with Z' Q Z reduced.
Result ZDecor(const Call& C)
{
	const Matrix& Q = C.M(0);
	std::string E;
	if (!ReadCovariance(Q, 2, &E) || Q.Rows != 2)
	{
		return Result::Fail(Q.Rows != 2 ? "covariance must be 2 x 2" : E);
	}
	double Z[2][2] = {{1, 0}, {0, 1}}, Q00 = Q(0, 0), Q01 = Q(0, 1), Q11 = Q(1, 1);
	if (!(Q00 > 0.0 && Q11 > 0.0 && Q00 * Q11 > Q01 * Q01))
	{
		return Result::Fail("covariance must be positive definite");
	}
	for (int It = 0; It < 200; ++It)
	{
		if (Q11 < Q00)
		{
			std::swap(Z[0][0], Z[0][1]);
			std::swap(Z[1][0], Z[1][1]);
			std::swap(Q00, Q11);
		}
		const double Mu = std::nearbyint(Q01 / Q00);
		if (Mu == 0.0)
		{
			return Result::Of(Matrix::FromRows({{Z[0][0], Z[0][1]}, {Z[1][0], Z[1][1]}}));
		}
		Z[0][1] -= Mu * Z[0][0];
		Z[1][1] -= Mu * Z[1][0];
		Q11 = Q11 - 2.0 * Mu * Q01 + Mu * Mu * Q00;
		Q01 = Q01 - Mu * Q00;
	}
	return Result::Fail("decorrelation did not terminate");
}

Result PBoot(const Call& C)
{
	const Matrix& Q = C.M(0);
	std::string E;
	if (!ReadCovariance(Q, 64, &E))
	{
		return Result::Fail(E);
	}
	Matrix L;
	if (!la::Cholesky(Q, &L, &E))
	{
		return Result::Fail("covariance must be positive definite");
	}
	double P = 1.0;
	for (int i = 0; i < Q.Rows; ++i)
	{
		P *= 2.0 * NormalCdf(1.0 / (2.0 * L(i, i))) - 1.0;
	}
	return Result::Of(P);
}

struct RoundIntegrand
{
	double S1, Slope, S2;
	double operator()(double X) const
	{
		const double Dens = std::exp(-0.5 * X * X / (S1 * S1)) / (S1 * std::sqrt(2.0 * cst::Pi));
		return Dens * NormalBand((-0.5 - Slope * X) / S2, (0.5 - Slope * X) / S2);
	}
};

double AdaptiveSimpson(const RoundIntegrand& F, double A, double B, double Fa, double Fm, double Fb, double Whole, double Tol, int Depth)
{
	const double M = 0.5 * (A + B), Lm = 0.5 * (A + M), Rm = 0.5 * (M + B);
	const double Flm = F(Lm), Frm = F(Rm);
	const double Left = (M - A) / 6.0 * (Fa + 4.0 * Flm + Fm), Right = (B - M) / 6.0 * (Fm + 4.0 * Frm + Fb);
	if (Depth <= 0 || std::fabs(Left + Right - Whole) <= 15.0 * Tol)
	{
		return Left + Right + (Left + Right - Whole) / 15.0;
	}
	return AdaptiveSimpson(F, A, M, Fa, Flm, Fm, Left, Tol / 2.0, Depth - 1) + AdaptiveSimpson(F, M, B, Fm, Frm, Fb, Right, Tol / 2.0, Depth - 1);
}

// Rounding success P(round(a_hat) = a): exact in 1-D, and in 2-D by integrating the conditional band.
Result PRound(const Call& C)
{
	const Matrix& Q = C.M(0);
	std::string E;
	if (!ReadCovariance(Q, 2, &E))
	{
		return Result::Fail(Q.IsSquare() ? "rounding success has a closed form only for n <= 2 (use pboot for a lower bound)" : E);
	}
	if (!(Q(0, 0) > 0.0))
	{
		return Result::Fail("covariance must be positive definite");
	}
	const double S1 = std::sqrt(Q(0, 0));
	if (Q.Rows == 1)
	{
		return Result::Of(2.0 * NormalCdf(1.0 / (2.0 * S1)) - 1.0);
	}
	const double Slope = Q(1, 0) / Q(0, 0), D2 = Q(1, 1) - Q(1, 0) * Q(1, 0) / Q(0, 0);
	if (!(D2 > 0.0))
	{
		return Result::Fail("covariance must be positive definite");
	}
	const RoundIntegrand F{S1, Slope, std::sqrt(D2)};
	std::vector<double> Cuts;
	for (int k = 0; k <= 16; ++k)
	{
		Cuts.push_back(-0.5 + k / 16.0);
	}
	if (Slope != 0.0)
	{
		for (double Edge : {0.5 / Slope, -0.5 / Slope})
		{
			if (std::fabs(Edge) < 0.5)
			{
				Cuts.push_back(Edge);
			}
		}
	}
	std::sort(Cuts.begin(), Cuts.end());
	double Total = 0.0;
	for (size_t k = 0; k + 1 < Cuts.size(); ++k)
	{
		const double A = Cuts[k], B = Cuts[k + 1];
		if (B - A <= 0.0)
		{
			continue;
		}
		const double Fa = F(A), Fb = F(B), Fm = F(0.5 * (A + B));
		Total += AdaptiveSimpson(F, A, B, Fa, Fm, Fb, (B - A) / 6.0 * (Fa + 4.0 * Fm + Fb), 1e-14, 40);
	}
	return Result::Of(Total);
}

// ---------------------------------------------------------------- ch07: the circle

Result Weyl(const Call& C)
{
	const double Alpha = C.N(0), N = C.N(1), Bins = C.N(2);
	if (!IsWhole(N, 1, 5e5) || !IsWhole(Bins, 1, 10000))
	{
		return Result::Fail("n must be a whole number 1..5e5 and bins 1..10000");
	}
	const int Count = static_cast<int>(N), Nb = static_cast<int>(Bins);
	std::vector<double> Pts(static_cast<size_t>(Count));
	std::vector<int> Hist(static_cast<size_t>(Nb), 0);
	for (int k = 1; k <= Count; ++k)
	{
		const double Y = Alpha * k;
		const double F = Y - std::floor(Y);
		Pts[static_cast<size_t>(k - 1)] = F;
		Hist[static_cast<size_t>(std::min(Nb - 1, static_cast<int>(F * Nb)))] += 1;
	}
	std::sort(Pts.begin(), Pts.end());
	double Star = 0.0;
	for (int i = 0; i < Count; ++i)
	{
		const double X = Pts[static_cast<size_t>(i)];
		Star = std::max(Star, std::max((i + 1.0) / Count - X, X - static_cast<double>(i) / Count));
	}
	double Dev = 0.0;
	for (int H : Hist)
	{
		Dev = std::max(Dev, std::fabs(static_cast<double>(H) / Count - 1.0 / Nb));
	}
	return Result::Record({{"dstar", Star}, {"dev", Dev}});
}

Result Unwrap(const Call& C)
{
	Matrix M = C.M(0);
	if (!M.IsVector())
	{
		return Result::Fail("phases must be a vector (rad)");
	}
	double Prev = M.V[0], Acc = M.V[0];
	for (size_t i = 1; i < M.V.size(); ++i)
	{
		const double D = M.V[i] - Prev;
		Prev = M.V[i];
		Acc += D - cst::Tau * std::floor((D + cst::Pi) / cst::Tau);
		M.V[i] = Acc;
	}
	return Result::Of(M);
}

// ---------------------------------------------------------------- ch10: Markov kernels

bool StochasticColumns(const Matrix& K, const char* Name, std::string* Error)
{
	if (K.IsEmpty())
	{
		*Error = std::string(Name) + " is empty";
		return false;
	}
	for (int j = 0; j < K.Cols; ++j)
	{
		double S = 0.0;
		for (int i = 0; i < K.Rows; ++i)
		{
			if (K(i, j) < -1e-12)
			{
				*Error = std::string(Name) + " has a negative entry: not a Markov kernel";
				return false;
			}
			S += K(i, j);
		}
		if (std::fabs(S - 1.0) > 1e-9)
		{
			*Error = std::string(Name) + " column " + std::to_string(j + 1) + " sums to " + std::to_string(S) + ": columns of a kernel K[y, x] = P(y | x) must sum to 1";
			return false;
		}
	}
	return true;
}

Result KComp(const Call& C)
{
	const Matrix &G = C.M(0), &F = C.M(1);
	std::string E;
	if (!StochasticColumns(G, "g", &E) || !StochasticColumns(F, "f", &E))
	{
		return Result::Fail(E);
	}
	if (G.Cols != F.Rows)
	{
		return Result::Fail("g is " + ShapeText(G) + " and f is " + ShapeText(F) + ": g after f needs cols(g) = rows(f)");
	}
	return Result::Of(la::Mul(G, F));
}

Result KPow(const Call& C)
{
	const Matrix& P = C.M(0);
	std::string E;
	if (!P.IsSquare() || !StochasticColumns(P, "P", &E))
	{
		return Result::Fail(P.IsSquare() ? E : "P must be a square kernel");
	}
	if (!IsWhole(C.N(1), 0, 1e9))
	{
		return Result::Fail("n must be a whole number >= 0");
	}
	long long N = static_cast<long long>(C.N(1));
	Matrix R = Matrix::Identity(P.Rows), B = P;
	while (N > 0)
	{
		if (N & 1)
		{
			R = la::Mul(R, B);
		}
		B = la::Mul(B, B);
		N >>= 1;
	}
	return Result::Of(R);
}

// Bayesian inverse of K: X -> Y under prior p: the kernel Y -> X with entries P(x | y).
Result BayesInv(const Call& C)
{
	const Matrix &K = C.M(0), &P = C.M(1);
	std::string E;
	if (!StochasticColumns(K, "K", &E))
	{
		return Result::Fail(E);
	}
	if (!P.IsVector() || P.Count() != K.Cols)
	{
		return Result::Fail("prior must be a vector with one entry per column of K (" + std::to_string(K.Cols) + ")");
	}
	double Sum = 0.0;
	for (double X : P.V)
	{
		if (X < 0.0)
		{
			return Result::Fail("prior must be nonnegative");
		}
		Sum += X;
	}
	if (std::fabs(Sum - 1.0) > 1e-9)
	{
		return Result::Fail("prior must sum to 1");
	}
	Matrix Out(K.Cols, K.Rows);
	for (int y = 0; y < K.Rows; ++y)
	{
		double Z = 0.0;
		for (int x = 0; x < K.Cols; ++x)
		{
			Z += K(y, x) * P.V[static_cast<size_t>(x)];
		}
		for (int x = 0; x < K.Cols; ++x)
		{
			Out(x, y) = Z > 0.0 ? K(y, x) * P.V[static_cast<size_t>(x)] / Z : 0.0;
		}
	}
	return Result::Of(Out);
}

// The chapter's motion kernel on a ring: one cell forward with 1 - 2s, two with s, stay with s.
Result RingKern(const Call& C)
{
	const double S = C.N(1);
	if (!IsWhole(C.N(0), 1, 500) || !(S >= 0.0 && S <= 0.5))
	{
		return Result::Fail("needs a whole n in 1..500 and slip 0 <= s <= 0.5");
	}
	const int N = C.Int(0);
	Matrix F(N, N);
	for (int x = 0; x < N; ++x)
	{
		F((x + 1) % N, x) += 1.0 - 2.0 * S;
		F((x + 2) % N, x) += S;
		F(x, x) += S;
	}
	return Result::Of(F);
}

// ---------------------------------------------------------------- ch10: transfer functions

std::vector<double> TrimPoly(std::vector<double> P)
{
	size_t K = 0;
	while (K + 1 < P.size() && P[K] == 0.0)
	{
		++K;
	}
	P.erase(P.begin(), P.begin() + static_cast<long>(K));
	if (P.empty())
	{
		P.push_back(0.0);
	}
	return P;
}

std::vector<double> PolyTimes(const std::vector<double>& A, const std::vector<double>& B)
{
	std::vector<double> Out(A.size() + B.size() - 1, 0.0);
	for (size_t i = 0; i < A.size(); ++i)
	{
		for (size_t j = 0; j < B.size(); ++j)
		{
			Out[i + j] += A[i] * B[j];
		}
	}
	return Out;
}

std::vector<double> PolyPlus(const std::vector<double>& A, const std::vector<double>& B)
{
	const size_t L = std::max(A.size(), B.size());
	std::vector<double> Out(L, 0.0);
	for (size_t i = 0; i < A.size(); ++i)
	{
		Out[L - A.size() + i] += A[i];
	}
	for (size_t i = 0; i < B.size(); ++i)
	{
		Out[L - B.size() + i] += B[i];
	}
	return Out;
}

FieldsComplex PolyAtS(const std::vector<double>& P, FieldsComplex S)
{
	FieldsComplex Acc(0.0, 0.0);
	for (double X : P)
	{
		Acc = Acc * S + X;
	}
	return Acc;
}

bool ReadTf(const Matrix& M, std::vector<double>* Num, std::vector<double>* Den, std::string* Error)
{
	if (M.IsScalar())
	{
		*Num = {M.V[0]};
		*Den = {1.0};
		return true;
	}
	if (M.Rows != 2 || M.Cols < 1)
	{
		*Error = "a transfer function is tf(num, den): a 2-row matrix [num; den], highest power first; got " + ShapeText(M);
		return false;
	}
	*Num = TrimPoly(std::vector<double>(M.V.begin(), M.V.begin() + M.Cols));
	*Den = TrimPoly(std::vector<double>(M.V.begin() + M.Cols, M.V.end()));
	if (Den->size() == 1 && (*Den)[0] == 0.0)
	{
		*Error = "denominator is zero";
		return false;
	}
	return true;
}

Matrix PackTf(const std::vector<double>& NumIn, const std::vector<double>& DenIn)
{
	const std::vector<double> Num = TrimPoly(NumIn), Den = TrimPoly(DenIn);
	const size_t L = std::max(Num.size(), Den.size());
	Matrix M(2, static_cast<int>(L));
	for (size_t i = 0; i < Num.size(); ++i)
	{
		M(0, static_cast<int>(L - Num.size() + i)) = Num[i];
	}
	for (size_t i = 0; i < Den.size(); ++i)
	{
		M(1, static_cast<int>(L - Den.size() + i)) = Den[i];
	}
	return M;
}

Result Tf(const Call& C)
{
	const Matrix &N = C.M(0), &D = C.M(1);
	if (!N.IsVector() || !D.IsVector())
	{
		return Result::Fail("num and den are coefficient vectors, highest power of s first");
	}
	const std::vector<double> Den = TrimPoly(D.V);
	if (Den.size() == 1 && Den[0] == 0.0)
	{
		return Result::Fail("denominator is zero");
	}
	return Result::Of(PackTf(N.V, Den));
}

Result TfSeries(const Call& C)
{
	std::vector<double> N1, D1, N2, D2;
	std::string E;
	if (!ReadTf(C.M(0), &N1, &D1, &E) || !ReadTf(C.M(1), &N2, &D2, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(PackTf(PolyTimes(N1, N2), PolyTimes(D1, D2)));
}

Result TfParallel(const Call& C)
{
	std::vector<double> N1, D1, N2, D2;
	std::string E;
	if (!ReadTf(C.M(0), &N1, &D1, &E) || !ReadTf(C.M(1), &N2, &D2, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(PackTf(PolyPlus(PolyTimes(N1, D2), PolyTimes(N2, D1)), PolyTimes(D1, D2)));
}

// Negative feedback G / (1 + G H); H = 1 is unity feedback.
Result TfFeedback(const Call& C)
{
	std::vector<double> Ng, Dg, Nh, Dh;
	std::string E;
	if (!ReadTf(C.M(0), &Ng, &Dg, &E) || !ReadTf(C.M(1), &Nh, &Dh, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<double> Den = TrimPoly(PolyPlus(PolyTimes(Dg, Dh), PolyTimes(Ng, Nh)));
	if (Den.size() == 1 && Den[0] == 0.0)
	{
		return Result::Fail("closed loop is degenerate (1 + G H = 0)");
	}
	return Result::Of(PackTf(PolyTimes(Ng, Dh), Den));
}

Result FrequencyRecord(FieldsComplex G)
{
	const double Mag = std::abs(G);
	return Result::Record({{"re", G.real()}, {"im", G.imag()}, {"mag", Mag}, {"db", 20.0 * std::log10(Mag)}, {"phase", std::arg(G)}});
}

Result TfJw(const Call& C)
{
	std::vector<double> N, D;
	std::string E;
	if (!ReadTf(C.M(0), &N, &D, &E))
	{
		return Result::Fail(E);
	}
	const FieldsComplex S(0.0, C.N(1));
	const FieldsComplex Den = PolyAtS(D, S);
	if (std::abs(Den) == 0.0)
	{
		return Result::Fail("s = jw is a pole");
	}
	return FrequencyRecord(PolyAtS(N, S) / Den);
}

// SISO state-space frequency response C (jw I - A)^-1 B + D.
Result SsJw(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &Cm = C.M(2), &D = C.M(3);
	const double W = C.N(4);
	const int N = A.Rows;
	if (!A.IsSquare() || B.Rows != N || B.Cols != 1 || Cm.Rows != 1 || Cm.Cols != N || !D.IsScalar())
	{
		return Result::Fail("needs A n x n, B n x 1, C 1 x n and scalar D (SISO)");
	}
	Matrix M(2 * N, 2 * N), Rhs(2 * N, 1);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			M(i, j) = -A(i, j);
			M(N + i, N + j) = -A(i, j);
		}
		M(i, N + i) = -W;
		M(N + i, i) = W;
		Rhs(i, 0) = B(i, 0);
	}
	Matrix X;
	std::string E;
	if (!la::Solve(M, Rhs, &X, &E))
	{
		return Result::Fail("s = jw is an eigenvalue of A");
	}
	double Re = D.V[0], Im = 0.0;
	for (int i = 0; i < N; ++i)
	{
		Re += Cm(0, i) * X(i, 0);
		Im += Cm(0, i) * X(N + i, 0);
	}
	return FrequencyRecord(FieldsComplex(Re, Im));
}

// ---------------------------------------------------------------- ch10: sheaves on a cycle

bool SheafLaplacian(const Matrix& Frames, double Err, Matrix* L, std::string* Error)
{
	if (!Frames.IsVector() || Frames.Count() > 40)
	{
		*Error = "frames must be a vector of sensor frame angles (rad), at most 40";
		return false;
	}
	const int N = Frames.Count();
	auto Rot = [](double T) { return Matrix::FromRows({{std::cos(T), -std::sin(T)}, {std::sin(T), std::cos(T)}}); };
	Matrix Delta(2 * N, 2 * N);
	for (int i = 0; i < N; ++i)
	{
		const int j = (i + 1) % N;
		const Matrix Fi = Rot(Frames.V[static_cast<size_t>(i)]);
		Matrix Fj = Rot(Frames.V[static_cast<size_t>(j)]);
		if (j == 0)
		{
			Fj = la::Mul(Rot(Err), Fj);
		}
		for (int r = 0; r < 2; ++r)
		{
			for (int c = 0; c < 2; ++c)
			{
				Delta(2 * i + r, 2 * i + c) += Fi(r, c);
				Delta(2 * i + r, 2 * j + c) -= Fj(r, c);
			}
		}
	}
	*L = la::Mul(la::Transpose(Delta), Delta);
	return true;
}

Result SheafLap(const Call& C)
{
	Matrix L;
	std::string E;
	if (!SheafLaplacian(C.M(0), C.N(1), &L, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(L);
}

Result SheafH0(const Call& C)
{
	Matrix L, V;
	std::vector<double> W;
	std::string E;
	if (!SheafLaplacian(C.M(0), C.N(1), &L, &E) || !la::SymmetricEigen(L, &W, &V, &E))
	{
		return Result::Fail(E);
	}
	int H0 = 0;
	double Lmin = FieldsNaN;
	for (double X : W)
	{
		if (X < 1e-9)
		{
			++H0;
		}
		else if (std::isnan(Lmin) || X < Lmin)
		{
			Lmin = X;
		}
	}
	return Result::Record({{"h0", static_cast<double>(H0)}, {"lmin", Lmin}});
}

// ---------------------------------------------------------------- ch11: cellular automata

bool RunAutomaton(int Rule, int N, int Steps, std::vector<uint8_t>* Row, std::vector<uint8_t>* Centre)
{
	uint8_t Tbl[8];
	for (int k = 0; k < 8; ++k)
	{
		Tbl[k] = static_cast<uint8_t>((Rule >> k) & 1);
	}
	Row->assign(static_cast<size_t>(N), 0);
	(*Row)[static_cast<size_t>(N / 2)] = 1;
	std::vector<uint8_t> Next(static_cast<size_t>(N));
	if (Centre)
	{
		Centre->assign(1, 1);
	}
	auto Cell = [&Tbl](uint8_t L, uint8_t M, uint8_t R) { return Tbl[4 * L + 2 * M + R]; };
	for (int t = 0; t < Steps; ++t)
	{
		const uint8_t* In = Row->data();
		uint8_t* Out = Next.data();
		if (N == 1)
		{
			Out[0] = Cell(In[0], In[0], In[0]);
		}
		else
		{
			Out[0] = Cell(In[N - 1], In[0], In[1]);
			for (int i = 1; i < N - 1; ++i)
			{
				Out[i] = Cell(In[i - 1], In[i], In[i + 1]);
			}
			Out[N - 1] = Cell(In[N - 2], In[N - 1], In[0]);
		}
		Row->swap(Next);
		if (Centre)
		{
			Centre->push_back((*Row)[static_cast<size_t>(N / 2)]);
		}
	}
	return true;
}

bool ReadRule(const Call& C, int* Rule, int* Steps, double MaxSteps, std::string* Error)
{
	if (!IsWhole(C.N(0), 0, 255))
	{
		*Error = "rule must be a whole number 0..255";
		return false;
	}
	if (!IsWhole(C.N(1), 0, MaxSteps))
	{
		*Error = "steps must be a whole number 0.." + std::to_string(static_cast<long>(MaxSteps));
		return false;
	}
	*Rule = C.Int(0);
	*Steps = C.Int(1);
	return true;
}

Result Ca(const Call& C)
{
	int Rule = 0, Steps = 0;
	std::string E;
	if (!ReadRule(C, &Rule, &Steps, 100000, &E))
	{
		return Result::Fail(E);
	}
	if (!IsWhole(C.N(2), 1, 20001) || C.N(1) * C.N(2) > 2e7)
	{
		return Result::Fail("width must be 1..20001 and steps * width at most 2e7");
	}
	const int N = C.Int(2);
	std::vector<uint8_t> Row;
	RunAutomaton(Rule, N, Steps, &Row, nullptr);
	Matrix M(1, N);
	for (int i = 0; i < N; ++i)
	{
		M.V[static_cast<size_t>(i)] = Row[static_cast<size_t>(i)];
	}
	return Result::Of(M);
}

// A width of 2 steps + 1 makes the periodic run equal to the infinite lattice at the centre.
bool CentreColumn(const Call& C, std::vector<uint8_t>* Centre, std::string* Error)
{
	int Rule = 0, Steps = 0;
	if (!ReadRule(C, &Rule, &Steps, 3000, Error))
	{
		return false;
	}
	std::vector<uint8_t> Row;
	RunAutomaton(Rule, 2 * Steps + 1, Steps, &Row, Centre);
	return true;
}

Result CaCenter(const Call& C)
{
	std::vector<uint8_t> Centre;
	std::string E;
	if (!CentreColumn(C, &Centre, &E))
	{
		return Result::Fail(E);
	}
	Matrix M(static_cast<int>(Centre.size()), 1);
	for (size_t i = 0; i < Centre.size(); ++i)
	{
		M.V[i] = Centre[i];
	}
	return Result::Of(M);
}

Result CaFreq(const Call& C)
{
	std::vector<uint8_t> Centre;
	std::string E;
	if (!CentreColumn(C, &Centre, &E))
	{
		return Result::Fail(E);
	}
	if (Centre.size() < 2)
	{
		return Result::Fail("needs at least one step");
	}
	double S = 0.0;
	for (size_t i = 1; i < Centre.size(); ++i)
	{
		S += Centre[i];
	}
	return Result::Of(S / static_cast<double>(Centre.size() - 1));
}

// ---------------------------------------------------------------- ch11: iterated function systems

Result Moran(const Call& C)
{
	const Matrix& R = C.M(0);
	if (!R.IsVector())
	{
		return Result::Fail("ratios must be a vector of contraction factors in (0, 1)");
	}
	for (double X : R.V)
	{
		if (!(X > 0.0 && X < 1.0))
		{
			return Result::Fail("every ratio must be in (0, 1)");
		}
	}
	auto Excess = [&R](double D) {
		double S = -1.0;
		for (double X : R.V)
		{
			S += std::pow(X, D);
		}
		return S;
	};
	if (R.Count() == 1)
	{
		return Result::Of(0.0);
	}
	double Lo = 0.0, Hi = 1.0;
	while (Excess(Hi) > 0.0 && Hi < 1e6)
	{
		Lo = Hi;
		Hi *= 2.0;
	}
	for (int It = 0; It < 200 && Hi - Lo > 1e-16 * Hi; ++It)
	{
		const double Mid = 0.5 * (Lo + Hi);
		(Excess(Mid) > 0.0 ? Lo : Hi) = Mid;
	}
	return Result::Of(0.5 * (Lo + Hi));
}

Result IfsPreset(const Call& C)
{
	static const double Fern[] = {0, 0, 0, 0.16, 0, 0, 0.01, 0.85, 0.04, -0.04, 0.85, 0, 1.6, 0.85, 0.2, -0.26, 0.23, 0.22, 0, 1.6, 0.07, -0.15, 0.28, 0.26, 0.24, 0, 0.44, 0.07};
	static const double Sierpinski[] = {0.5, 0, 0, 0.5, 0, 0, 0.3333, 0.5, 0, 0, 0.5, 0.5, 0, 0.3333, 0.5, 0, 0, 0.5, 0.25, 0.433, 0.3334};
	static const double Tree[] = {0, 0, 0, 0.5, 0, 0, 0.05, 0.42, -0.42, 0.42, 0.42, 0, 0.2, 0.4, 0.42, 0.42, -0.42, 0.42, 0, 0.2, 0.4, 0.1, 0, 0, 0.1, 0, 0.2, 0.15};
	static const double Expand[] = {0.5, 0, 0, 0.5, 0, 0, 0.3333, 0.5, 0, 0, 0.5, 0.5, 0, 0.3333, 1.05, 0, 0, 1.05, 0.25, 0.433, 0.3334};
	static const double AvgExpand[] = {0.5, 0, 0, 0.5, 0, 0, 0.1, 1.2, 0, 0, 1.2, 0.5, 0, 0.9};
	const double* Tables[] = {Fern, Sierpinski, Tree, Expand, AvgExpand};
	const int Rows[] = {4, 3, 4, 3, 2};
	if (!IsWhole(C.N(0), 1, 5))
	{
		return Result::Fail("preset 1 fern, 2 Sierpinski, 3 tree, 4 one expanding map, 5 expanding on average");
	}
	const int K = C.Int(0) - 1;
	Matrix M(Rows[K], 7);
	for (int i = 0; i < M.Count(); ++i)
	{
		M.V[static_cast<size_t>(i)] = Tables[K][i];
	}
	return Result::Of(M);
}

bool ReadIfs(const Matrix& M, std::vector<double>* Norms, std::vector<double>* Weights, std::string* Error)
{
	if (M.Rows < 1 || (M.Cols != 6 && M.Cols != 7))
	{
		*Error = "maps are rows [a b c d e f] or [a b c d e f p]: x' = a x + b y + e, y' = c x + d y + f";
		return false;
	}
	Norms->clear();
	Weights->clear();
	double Sum = 0.0;
	for (int i = 0; i < M.Rows; ++i)
	{
		const double A = M(i, 0), B = M(i, 1), Cc = M(i, 2), D = M(i, 3);
		const double S1 = A * A + B * B + Cc * Cc + D * D, Det = A * D - B * Cc;
		Norms->push_back(std::sqrt((S1 + std::sqrt(std::max(0.0, S1 * S1 - 4.0 * Det * Det))) / 2.0));
		Weights->push_back(M.Cols == 7 ? std::max(0.0, M(i, 6)) : 1.0);
		Sum += Weights->back();
	}
	for (double& W : *Weights)
	{
		W /= Sum > 0.0 ? Sum : 1.0;
	}
	return true;
}

Result IfsNorm(const Call& C)
{
	std::vector<double> Norms, Weights;
	std::string E;
	if (!ReadIfs(C.M(0), &Norms, &Weights, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Matrix::Column(Norms));
}

Result IfsLog(const Call& C)
{
	std::vector<double> Norms, Weights;
	std::string E;
	if (!ReadIfs(C.M(0), &Norms, &Weights, &E))
	{
		return Result::Fail(E);
	}
	double S = 0.0;
	for (size_t i = 0; i < Norms.size(); ++i)
	{
		S += Weights[i] * std::log(std::max(Norms[i], 1e-300));
	}
	return Result::Of(S);
}

// ---------------------------------------------------------------- ch11: Turing, L-systems, splines, Lloyd

Result Turing(const Call& C)
{
	const double Dv = C.N(0), Du = C.N(1), A = C.N(2), B = C.N(3);
	if (!(Du > 0.0 && Dv > 0.0 && A + B > 0.0))
	{
		return Result::Fail("needs positive diffusivities and a + b > 0");
	}
	const double Us = A + B, Vs = B / (Us * Us);
	const double J11 = -1.0 + 2.0 * Us * Vs, J12 = Us * Us, J21 = -2.0 * Us * Vs, J22 = -Us * Us;
	if (!(J11 + J22 < 0.0 && J11 * J22 - J12 * J21 > 0.0))
	{
		return Result::Fail("the uniform state is unstable without diffusion (Turing needs tr J < 0 and det J > 0)");
	}
	auto Growth = [&](double Q) {
		const double P = J11 - Du * Q, S = J22 - Dv * Q, Tr = P + S, Det = P * S - J12 * J21, Disc = Tr * Tr - 4.0 * Det;
		return Disc >= 0.0 ? (Tr + std::sqrt(Disc)) / 2.0 : Tr / 2.0;
	};
	// Growth > 0 exactly where det(J - D k^2) < 0: a quadratic in q = k^2.
	const double A2 = Du * Dv, A1 = -(J11 * Dv + J22 * Du), A0 = J11 * J22 - J12 * J21, Disc = A1 * A1 - 4.0 * A2 * A0;
	if (!(A1 < 0.0 && Disc > 0.0))
	{
		return Result::Record({{"kmin", FieldsNaN}, {"kmax", FieldsNaN}, {"kstar", FieldsNaN}, {"wl", FieldsNaN}, {"gmax", FieldsNaN}});
	}
	const double Q1 = (-A1 - std::sqrt(Disc)) / (2.0 * A2), Q2 = (-A1 + std::sqrt(Disc)) / (2.0 * A2);
	double Lo = Q1, Hi = Q2;
	const double Golden = (std::sqrt(5.0) - 1.0) / 2.0;
	for (int It = 0; It < 300 && Hi - Lo > 1e-15 * Hi; ++It)
	{
		const double X1 = Hi - Golden * (Hi - Lo), X2 = Lo + Golden * (Hi - Lo);
		if (Growth(X1) < Growth(X2))
		{
			Lo = X1;
		}
		else
		{
			Hi = X2;
		}
	}
	const double Qs = 0.5 * (Lo + Hi), Ks = std::sqrt(Qs);
	return Result::Record({{"kmin", std::sqrt(Q1)}, {"kmax", std::sqrt(Q2)}, {"kstar", Ks}, {"wl", cst::Tau / Ks}, {"gmax", Growth(Qs)}});
}

Result LSys(const Call& C)
{
	if (!IsWhole(C.N(0), 0, 74))
	{
		return Result::Fail("n must be a whole number 0..74 (lengths stay exact below 2^53)");
	}
	double A = 1.0, B = 0.0;
	for (int k = 0; k < C.Int(0); ++k)
	{
		const double NewA = A + B;   // a -> ab, b -> a
		B = A;
		A = NewA;
	}
	return Result::Record({{"len", A + B}, {"a", A}, {"b", B}, {"ratio", B > 0.0 ? A / B : FieldsInf}});
}

bool NaturalSpline(const Matrix& Xs, const Matrix& Ys, std::vector<double>* M, std::string* Error)
{
	if (!Xs.IsVector() || !Ys.IsVector() || Xs.Count() != Ys.Count() || Xs.Count() < 2)
	{
		*Error = "knots and values must be vectors of equal length >= 2";
		return false;
	}
	const int N = Xs.Count() - 1;
	const std::vector<double>&X = Xs.V, &Y = Ys.V;
	for (int i = 0; i < N; ++i)
	{
		if (!(X[static_cast<size_t>(i + 1)] > X[static_cast<size_t>(i)]))
		{
			*Error = "knots must be strictly increasing";
			return false;
		}
	}
	// Tridiagonal system for the interior second derivatives (Thomas algorithm); M_0 = M_n = 0.
	M->assign(static_cast<size_t>(N + 1), 0.0);
	if (N < 2)
	{
		return true;
	}
	std::vector<double> Diag(static_cast<size_t>(N + 1)), Up(static_cast<size_t>(N + 1)), Rhs(static_cast<size_t>(N + 1));
	for (int i = 1; i < N; ++i)
	{
		const size_t u = static_cast<size_t>(i);
		const double H0 = X[u] - X[u - 1], H1 = X[u + 1] - X[u];
		Diag[u] = 2.0 * (H0 + H1);
		Up[u] = H1;
		Rhs[u] = 6.0 * ((Y[u + 1] - Y[u]) / H1 - (Y[u] - Y[u - 1]) / H0);
		if (i > 1)
		{
			const double W = H0 / Diag[u - 1];
			Diag[u] -= W * Up[u - 1];
			Rhs[u] -= W * Rhs[u - 1];
		}
	}
	for (int i = N - 1; i >= 1; --i)
	{
		const size_t u = static_cast<size_t>(i);
		(*M)[u] = (Rhs[u] - (i < N - 1 ? Up[u] * (*M)[u + 1] : 0.0)) / Diag[u];
	}
	return true;
}

Result Spline(const Call& C)
{
	std::vector<double> M;
	std::string E;
	if (!NaturalSpline(C.M(0), C.M(1), &M, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<double>&X = C.M(0).V, &Y = C.M(1).V;
	const int N = static_cast<int>(X.size()) - 1;
	Matrix Out = C.M(2);
	for (double& T : Out.V)
	{
		int i = static_cast<int>(std::upper_bound(X.begin(), X.end(), T) - X.begin()) - 1;
		i = std::max(0, std::min(N - 1, i));
		const size_t u = static_cast<size_t>(i);
		const double H = X[u + 1] - X[u], A = X[u + 1] - T, B = T - X[u];
		T = M[u] * A * A * A / (6.0 * H) + M[u + 1] * B * B * B / (6.0 * H) + (Y[u] / H - M[u] * H / 6.0) * A + (Y[u + 1] / H - M[u + 1] * H / 6.0) * B;
	}
	return Result::Of(Out);
}

Result SplineEnergy(const Call& C)
{
	std::vector<double> M;
	std::string E;
	if (!NaturalSpline(C.M(0), C.M(1), &M, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<double>& X = C.M(0).V;
	double S = 0.0;
	for (size_t i = 0; i + 1 < X.size(); ++i)
	{
		S += (X[i + 1] - X[i]) * (M[i] * M[i] + M[i] * M[i + 1] + M[i + 1] * M[i + 1]) / 3.0;   // f'' is linear on each piece
	}
	return Result::Of(S);
}

// The chapter's drone-coverage bench: 12 sites from mulberry32(12) on a 90 x 90 grid.
Result Lloyd(const Call& C)
{
	const double Hx = C.N(0);
	if (!IsWhole(C.N(1), 0, 150))
	{
		return Result::Fail("steps must be a whole number 0..150");
	}
	const int Grid = 90, Sites = 12, Steps = C.Int(1);
	std::vector<double> Sx(Sites), Sy(Sites), Cx(Sites), Cy(Sites);
	std::vector<bool> Has(Sites);
	FieldsMulberry Rng(12);
	for (int k = 0; k < Sites; ++k)
	{
		Sx[static_cast<size_t>(k)] = Rng.Next();
		Sy[static_cast<size_t>(k)] = Rng.Next();
	}
	auto Evaluate = [&]() {
		std::vector<double> Wx(Sites, 0.0), Wy(Sites, 0.0), Wm(Sites, 0.0);
		double En = 0.0, Mass = 0.0;
		for (int j = 0; j < Grid; ++j)
		{
			for (int i = 0; i < Grid; ++i)
			{
				const double X = (i + 0.5) / Grid, Y = (j + 0.5) / Grid;
				const double Rho = std::exp(-((X - Hx) * (X - Hx) + (Y - 0.3) * (Y - 0.3)) / 0.05) + 0.2;
				int Best = 0;
				double Bd = FieldsInf;
				for (int k = 0; k < Sites; ++k)
				{
					const double D = (X - Sx[static_cast<size_t>(k)]) * (X - Sx[static_cast<size_t>(k)]) + (Y - Sy[static_cast<size_t>(k)]) * (Y - Sy[static_cast<size_t>(k)]);
					if (D < Bd)
					{
						Bd = D;
						Best = k;
					}
				}
				En += Rho * Bd;
				Mass += Rho;
				Wx[static_cast<size_t>(Best)] += Rho * X;
				Wy[static_cast<size_t>(Best)] += Rho * Y;
				Wm[static_cast<size_t>(Best)] += Rho;
			}
		}
		for (size_t k = 0; k < static_cast<size_t>(Sites); ++k)
		{
			Has[k] = Wm[k] > 0.0;
			Cx[k] = Has[k] ? Wx[k] / Wm[k] : Sx[k];
			Cy[k] = Has[k] ? Wy[k] / Wm[k] : Sy[k];
		}
		return En / Mass;
	};
	const double E0 = Evaluate();
	double En = E0;
	for (int s = 0; s < Steps; ++s)
	{
		Sx = Cx;
		Sy = Cy;
		En = Evaluate();
	}
	double Gap = 0.0;
	for (size_t k = 0; k < static_cast<size_t>(Sites); ++k)
	{
		if (Has[k])
		{
			Gap = std::max(Gap, std::hypot(Sx[k] - Cx[k], Sy[k] - Cy[k]));
		}
	}
	return Result::Record({{"e0", E0}, {"e", En}, {"gap", Gap}});
}

// ---------------------------------------------------------------- ch12: posterior Cramer-Rao bound

bool PcrbRun(const Matrix& F, const Matrix& Q, const Matrix& H, const Matrix& R, const Matrix& P0, int Steps, Matrix* P, std::string* Error)
{
	const int N = F.Rows;
	if (!F.IsSquare() || Q.Rows != N || Q.Cols != N || H.Cols != N || R.Rows != H.Rows || R.Cols != H.Rows || P0.Rows != N || P0.Cols != N)
	{
		*Error = "needs F n x n, Q n x n, H m x n, R m x m, P0 n x n";
		return false;
	}
	Matrix Ri, J;
	if (!la::Inverse(R, &Ri, Error) || !la::Inverse(P0, &J, Error))
	{
		*Error = "R and P0 must be invertible";
		return false;
	}
	const Matrix Info = la::Mul(la::Mul(la::Transpose(H), Ri), H);
	for (int k = 0; k < Steps; ++k)
	{
		Matrix Jinv, Pred;
		if (!la::Inverse(J, &Jinv, Error))
		{
			return false;
		}
		if (!la::Inverse(la::Add(Q, la::Mul(la::Mul(F, Jinv), la::Transpose(F))), &Pred, Error))
		{
			*Error = "predicted covariance is singular";
			return false;
		}
		J = la::Add(Pred, Info);
	}
	if (!la::Inverse(J, P, Error))
	{
		return false;
	}
	*P = la::Symmetrize(*P);
	return true;
}

Result Pcrb(const Call& C)
{
	if (!IsWhole(C.N(5), 0, 100000))
	{
		return Result::Fail("k must be a whole number of steps");
	}
	Matrix P;
	std::string E;
	if (!PcrbRun(C.M(0), C.M(1), C.M(2), C.M(3), C.M(4), C.Int(5), &P, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(P);
}

Result PcrbCv(const Call& C)
{
	const double Q = C.N(1), R = C.N(2), Dt = C.N(3);
	if (!IsWhole(C.N(0), 0, 100000) || !(Q > 0.0 && R > 0.0 && Dt > 0.0))
	{
		return Result::Fail("needs whole k >= 0 and positive q, r, dt");
	}
	const Matrix F = Matrix::FromRows({{1.0, Dt}, {0.0, 1.0}});
	const Matrix Qd = Matrix::FromRows({{Q * Dt * Dt * Dt / 3.0, Q * Dt * Dt / 2.0}, {Q * Dt * Dt / 2.0, Q * Dt}});
	Matrix P;
	std::string E;
	if (!PcrbRun(F, Qd, Matrix::FromRows({{1.0, 0.0}}), Matrix::Scalar(R), Matrix::FromRows({{4.0, 0.0}, {0.0, 1.0}}), C.Int(0), &P, &E))
	{
		return Result::Fail(E);
	}
	return Result::Record({{"p11", P(0, 0)}, {"p12", P(0, 1)}, {"p22", P(1, 1)}});
}

// ---------------------------------------------------------------- ch12: minimum-time HJB

Result BangBang(const Call& C)
{
	const double X = C.N(0), V = C.N(1), S = X + V * std::fabs(V) / 2.0;
	if (X == 0.0 && V == 0.0)
	{
		return Result::Record({{"T", 0.0}, {"u0", 0.0}, {"ts", 0.0}, {"xs", 0.0}, {"vs", 0.0}});
	}
	if (S > 0.0 || (S == 0.0 && V < 0.0))
	{
		const double Root = std::sqrt(X + V * V / 2.0);
		if (S == 0.0)
		{
			return Result::Record({{"T", V + 2.0 * Root}, {"u0", 1.0}, {"ts", 0.0}, {"xs", X}, {"vs", V}});
		}
		return Result::Record({{"T", V + 2.0 * Root}, {"u0", -1.0}, {"ts", V + Root}, {"xs", Root * Root / 2.0}, {"vs", -Root}});
	}
	const double Root = std::sqrt(std::max(0.0, -X + V * V / 2.0));
	if (S == 0.0)
	{
		return Result::Record({{"T", -V + 2.0 * Root}, {"u0", -1.0}, {"ts", 0.0}, {"xs", X}, {"vs", V}});
	}
	return Result::Record({{"T", -V + 2.0 * Root}, {"u0", 1.0}, {"ts", Root - V}, {"xs", -Root * Root / 2.0}, {"vs", Root}});
}

// The chapter's simulation: switch on the curve x = -v|v|/2, semi-implicit Euler, stop within 1e-3 of rest.
Result BbSim(const Call& C)
{
	double X = C.N(0), V = C.N(1);
	const double Dt = C.N(2);
	if (!(Dt >= 1e-5 && Dt <= 0.1))
	{
		return Result::Fail("dt must be in [1e-5, 0.1] s");
	}
	double T = 0.0;
	while (X * X + V * V > 1e-6 && T < 50.0)
	{
		const double S = X + V * std::fabs(V) / 2.0;
		const double U = S > 0.0 ? -1.0 : (S < 0.0 ? 1.0 : (V > 0.0 ? -1.0 : (V < 0.0 ? 1.0 : 0.0)));
		V += U * Dt;
		X += V * Dt;
		T += Dt;
	}
	return Result::Of(T);
}

// ---------------------------------------------------------------- ch12: Fokker-Planck and Ito

Result OuProc(const Call& C)
{
	const double X0 = C.N(0), Th = C.N(1), Sg = C.N(2), T = C.N(3);
	if (!(T >= 0.0))
	{
		return Result::Fail("t must be >= 0");
	}
	const double Var = Th == 0.0 ? Sg * Sg * T : Sg * Sg * -std::expm1(-2.0 * Th * T) / (2.0 * Th);
	return Result::Record({{"mean", X0 * std::exp(-Th * T)}, {"var", Var}, {"sd", std::sqrt(Var)}, {"vinf", Th > 0.0 ? Sg * Sg / (2.0 * Th) : FieldsInf}});
}

Result QCv(const Call& C)
{
	const double Q = C.N(0), Dt = C.N(1);
	return Result::Of(Matrix::FromRows({{Q * Dt * Dt * Dt / 3.0, Q * Dt * Dt / 2.0}, {Q * Dt * Dt / 2.0, Q * Dt}}));
}

// Van Loan: expm([[-A, Qc], [0, A']] dt) = [[., Phi^-1 Qd], [0, Phi']].
Result QDisc(const Call& C)
{
	const Matrix &A = C.M(0), &Qc = C.M(1);
	const double Dt = C.N(2);
	const int N = A.Rows;
	if (!A.IsSquare() || Qc.Rows != N || Qc.Cols != N || N > 20)
	{
		return Result::Fail("needs A n x n and Qc n x n (n <= 20)");
	}
	Matrix M(2 * N, 2 * N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			M(i, j) = -A(i, j) * Dt;
			M(i, N + j) = Qc(i, j) * Dt;
			M(N + i, N + j) = A(j, i) * Dt;
		}
	}
	Matrix Ex;
	std::string E;
	if (!la::Expm(M, &Ex, &E))
	{
		return Result::Fail(E);
	}
	Matrix PhiT(N, N), Upper(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			PhiT(i, j) = Ex(N + j, N + i);
			Upper(i, j) = Ex(i, N + j);
		}
	}
	return Result::Of(la::Symmetrize(la::Mul(PhiT, Upper)));
}

// ---------------------------------------------------------------- ch12: tropical algebra

Matrix MinPlus(const Matrix& A, const Matrix& B)
{
	Matrix Out(A.Rows, B.Cols, FieldsInf);
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int k = 0; k < A.Cols; ++k)
		{
			const double Aik = A(i, k);
			if (Aik == FieldsInf)
			{
				continue;
			}
			for (int j = 0; j < B.Cols; ++j)
			{
				Out(i, j) = std::min(Out(i, j), Aik + B(k, j));
			}
		}
	}
	return Out;
}

Result TropMul(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1);
	if (A.Cols != B.Rows || A.IsEmpty() || B.IsEmpty())
	{
		return Result::Fail("(min, +) product needs cols(A) = rows(B)");
	}
	return Result::Of(MinPlus(A, B));
}

Result TropPow(const Call& C)
{
	const Matrix& W = C.M(0);
	if (!W.IsSquare() || W.Rows > 200)
	{
		return Result::Fail("W must be a square edge-weight matrix (inf = no edge), at most 200 x 200");
	}
	double K = C.N(1) < 0.0 ? W.Rows - 1.0 : C.N(1);
	if (!IsWhole(K, 0, 1e9))
	{
		return Result::Fail("k must be a whole number >= 0 (default n - 1: shortest paths)");
	}
	Matrix R(W.Rows, W.Rows, FieldsInf), B = W;
	for (int i = 0; i < W.Rows; ++i)
	{
		R(i, i) = 0.0;
	}
	long long N = static_cast<long long>(K);
	while (N > 0)
	{
		if (N & 1)
		{
			R = MinPlus(R, B);
		}
		N >>= 1;
		if (N > 0)
		{
			B = MinPlus(B, B);
		}
	}
	return Result::Of(R);
}

// ---------------------------------------------------------------- ch12: optimal transport

// Sinkhorn with uniform unit marginals (u = 1 / K v, v = 1 / K' u), run in the log domain.
bool SinkhornPlan(const Matrix& Cm, double Eps, int Iters, Matrix* Plan, std::string* Error)
{
	const int N = Cm.Rows;
	if (!Cm.IsSquare() || N > 40)
	{
		*Error = "cost must be a square matrix (equal swarms), at most 40 x 40";
		return false;
	}
	if (!(Eps > 0.0))
	{
		*Error = "epsilon must be positive";
		return false;
	}
	std::vector<double> F(static_cast<size_t>(N), 0.0), G(static_cast<size_t>(N), 0.0), Tmp(static_cast<size_t>(N));
	auto LogSumExp = [&Tmp, N]() {
		double M = -FieldsInf;
		for (int k = 0; k < N; ++k)
		{
			M = std::max(M, Tmp[static_cast<size_t>(k)]);
		}
		if (M == -FieldsInf)
		{
			return M;
		}
		double S = 0.0;
		for (int k = 0; k < N; ++k)
		{
			S += std::exp(Tmp[static_cast<size_t>(k)] - M);
		}
		return M + std::log(S);
	};
	for (int It = 0; It < Iters; ++It)
	{
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				Tmp[static_cast<size_t>(j)] = (G[static_cast<size_t>(j)] - Cm(i, j)) / Eps;
			}
			F[static_cast<size_t>(i)] = -Eps * LogSumExp();
		}
		for (int j = 0; j < N; ++j)
		{
			for (int i = 0; i < N; ++i)
			{
				Tmp[static_cast<size_t>(i)] = (F[static_cast<size_t>(i)] - Cm(i, j)) / Eps;
			}
			G[static_cast<size_t>(j)] = -Eps * LogSumExp();
		}
		if (It % 10 == 9)
		{
			double Worst = 0.0;
			for (int i = 0; i < N; ++i)
			{
				double S = 0.0;
				for (int j = 0; j < N; ++j)
				{
					S += std::exp((F[static_cast<size_t>(i)] + G[static_cast<size_t>(j)] - Cm(i, j)) / Eps);
				}
				Worst = std::max(Worst, std::fabs(S - 1.0));
			}
			if (Worst < 1e-15)
			{
				break;
			}
		}
	}
	*Plan = Matrix(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			(*Plan)(i, j) = std::exp((F[static_cast<size_t>(i)] + G[static_cast<size_t>(j)] - Cm(i, j)) / Eps);
		}
	}
	if (!std::isfinite(la::MaxAbs(*Plan)))
	{
		*Error = "Sinkhorn diverged (a row or column has no finite cost)";
		return false;
	}
	return true;
}

bool ReadSinkhorn(const Call& C, Matrix* Plan, std::string* Error)
{
	if (!IsWhole(C.N(2), 1, 100000))
	{
		*Error = "iterations must be a whole number 1..100000";
		return false;
	}
	if (C.N(2) * C.M(0).Count() > 2e6)
	{
		*Error = "iterations * n^2 must be at most 2e6: pass fewer iterations for a large swarm";
		return false;
	}
	return SinkhornPlan(C.M(0), C.N(1), C.Int(2), Plan, Error);
}

Result Sinkhorn(const Call& C)
{
	Matrix Plan;
	std::string E;
	if (!ReadSinkhorn(C, &Plan, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Plan);
}

Result SinkCost(const Call& C)
{
	Matrix Plan;
	std::string E;
	if (!ReadSinkhorn(C, &Plan, &E))
	{
		return Result::Fail(E);
	}
	double S = 0.0;
	for (size_t i = 0; i < Plan.V.size(); ++i)
	{
		S += Plan.V[i] * C.M(0).V[i];
	}
	return Result::Of(S);
}

// Optimal assignment by dynamic programming over subsets (Birkhoff-von Neumann: the optimum is a permutation).
Result OtExact(const Call& C)
{
	const Matrix& Cm = C.M(0);
	if (!Cm.IsSquare() || Cm.Rows > 16)
	{
		return Result::Fail("cost must be square, at most 16 x 16");
	}
	const int N = Cm.Rows;
	std::vector<double> Best(static_cast<size_t>(1) << N, FieldsInf);
	Best[0] = 0.0;
	for (unsigned Mask = 0; Mask < (1u << N); ++Mask)
	{
		if (Best[Mask] == FieldsInf)
		{
			continue;
		}
		int Row = 0;
		for (unsigned Bits = Mask; Bits; Bits &= Bits - 1)
		{
			++Row;
		}
		if (Row == N)
		{
			continue;
		}
		for (int j = 0; j < N; ++j)
		{
			if (!(Mask & (1u << j)))
			{
				const unsigned Next = Mask | (1u << j);
				Best[Next] = std::min(Best[Next], Best[Mask] + Cm(Row, j));
			}
		}
	}
	return Result::Of(Best[(1u << N) - 1]);
}

// ---------------------------------------------------------------- ch12: symplectic integration, random matrices

double KeplerWorst(bool Leapfrog, int Orbits, int StepsPer)
{
	double Rx = 1.0, Ry = 0.0, Vx = 0.0, Vy = 1.0, Worst = 0.0;
	const double H = cst::Tau / StepsPer, E0 = 0.5 - 1.0;
	auto Accel = [](double X, double Y, double* Ax, double* Ay) {
		const double R = std::sqrt(X * X + Y * Y), R3 = R * R * R;
		*Ax = -X / R3;
		*Ay = -Y / R3;
	};
	for (long k = 0; k < static_cast<long>(Orbits) * StepsPer; ++k)
	{
		double Ax = 0.0, Ay = 0.0;
		Accel(Rx, Ry, &Ax, &Ay);
		if (Leapfrog)
		{
			Vx += 0.5 * H * Ax;
			Vy += 0.5 * H * Ay;
			Rx += H * Vx;
			Ry += H * Vy;
			Accel(Rx, Ry, &Ax, &Ay);
			Vx += 0.5 * H * Ax;
			Vy += 0.5 * H * Ay;
		}
		else
		{
			Rx += H * Vx;
			Ry += H * Vy;
			Vx += H * Ax;
			Vy += H * Ay;
		}
		const double En = 0.5 * (Vx * Vx + Vy * Vy) - 1.0 / std::sqrt(Rx * Rx + Ry * Ry);
		Worst = std::max(Worst, std::fabs(En - E0));
		if (Worst > 1.0)
		{
			break;
		}
	}
	return Worst;
}

Result KepEnergy(const Call& C)
{
	if (!IsWhole(C.N(0), 1, 10000) || !IsWhole(C.N(1), 4, 100000) || C.N(0) * C.N(1) > 1e6)
	{
		return Result::Fail("orbits and steps per orbit must be whole numbers, orbits * steps <= 1e6");
	}
	return Result::Record({{"lf", KeplerWorst(true, C.Int(0), C.Int(1))}, {"eu", KeplerWorst(false, C.Int(0), C.Int(1))}});
}

Result MpEdges(const Call& C)
{
	const double Cr = C.N(0), S2 = C.N(1);
	if (!(Cr > 0.0 && S2 > 0.0))
	{
		return Result::Fail("aspect ratio c = p/n and variance must be positive");
	}
	return Result::Record({{"lo", S2 * std::pow(1.0 - std::sqrt(Cr), 2)}, {"hi", S2 * std::pow(1.0 + std::sqrt(Cr), 2)}});
}

Result MpPdf(const Call& C)
{
	const double Cr = C.N(1), S2 = C.N(2);
	if (!(Cr > 0.0 && S2 > 0.0))
	{
		return Result::Fail("aspect ratio c = p/n and variance must be positive");
	}
	const double Lo = S2 * std::pow(1.0 - std::sqrt(Cr), 2), Hi = S2 * std::pow(1.0 + std::sqrt(Cr), 2);
	Matrix X = C.M(0);
	for (double& V : X.V)
	{
		V = (V > Lo && V < Hi) ? std::sqrt((Hi - V) * (V - Lo)) / (2.0 * cst::Pi * S2 * Cr * V) : 0.0;
	}
	return Result::Of(X);
}
} // namespace fields_detail

void RegisterFieldsLibrary(Registry& R)
{
	using namespace fields_detail;
	const char* Lat = "lattices";
	const char* Cat = "category";
	const char* Pro = "procedural";
	const char* Exp = "expansion";

	// Chapter 7: lattices and GNSS.
	Def(R, "tofrange", Lat, "ch07 Lattices: Ranging by time of flight", "range from travel time, c dt / way (way 2: radar echo)", "range, m",
		{Num("dt", "s"), Num("way", "-", 2)}, &TofRange);
	Def(R, "tofres", Lat, "ch07 Lattices: Ranging by time of flight", "two-way range resolution c / (2 B)", "resolution, m", {Num("B", "Hz")}, &TofRes);
	Def(R, "widelane", Lat, "ch07 Lattices: Ranging by time of flight", "carrier wavelengths and the wide/narrow-lane combinations",
		"record l1, l2, wide, narrow (m)", {Num("f1", "Hz", 1575.42e6), Num("f2", "Hz", 1227.60e6)}, &Widelane);
	Def(R, "satrange", Lat, "ch07 Lattices: Ranging by time of flight", "slant range and travel time to a satellite at elevation el",
		"record d (m), t (s)", {Num("el", "rad"), Num("re", "m", 6371000.0), Num("rs", "m", 26560000.0)}, &SatRange);
	Def(R, "gnssfix", Lat, "ch07 Lattices: bench 1, pseudorange positioning", "Gauss-Newton position and clock from pseudoranges (start at Earth's centre)",
		"record x, y, z (m), cb = c*bias (m), dt (s), it", {Mat("S", "m (n x 3)"), Mat("rho", "m")}, &GnssFix);
	Def(R, "dop", Lat, "ch07 Lattices: bench 1, dilution of precision", "DOPs from satellite positions seen from receiver p (vertical = p/|p|)",
		"record gdop, pdop, hdop, vdop, tdop", {Mat("S", "m (n x 3)"), Mat("p", "m (3)")}, &Dop);
	Def(R, "dopae", Lat, "ch07 Lattices: bench 1, dilution of precision", "DOPs from satellite azimuths and elevations", "record gdop, pdop, hdop, vdop, tdop",
		{Mat("az", "rad"), Mat("el", "rad")}, &DopAe);
	Def(R, "baarda", Lat, "ch07 Lattices: Seven ways ranging goes wrong (normalized residuals)", "Baarda w-statistics of the pseudorange fix, r_i / (sigma sqrt(1 - h_ii))",
		"column of w, one per satellite", {Mat("S", "m (n x 3)"), Mat("rho", "m"), Num("sigma", "m", 1.0)}, &Baarda);
	Def(R, "gaussred", Lat, "ch07 Lattices: bench 3, Lagrange-Gauss reduction", "Lagrange-Gauss reduction of a 2-D basis (rows); row 1 is a shortest vector",
		"2 x 2 reduced basis", {Mat("B", "rows = basis vectors")}, &GaussReduce);
	Def(R, "lll", Lat, "ch07 Lattices: Reducing a lattice basis (LLL)", "LLL reduction of the row basis B", "reduced basis, same shape",
		{Mat("B", "rows = basis vectors"), Num("delta", "-", 0.75)}, &Lll);
	Def(R, "shortvec", Lat, "ch07 Lattices: Volume, reduction, best integer guess", "shortest nonzero lattice vector length lambda_1 (Fincke-Pohst)",
		"lambda_1, units of B", {Mat("B", "rows = basis vectors")}, &ShortVec);
	Def(R, "lllbound", Lat, "ch07 Lattices: Reducing a lattice basis (LLL)", "LLL guarantee |b1| <= (1/(delta - 1/4))^((n-1)/2) lambda_1", "bound factor",
		{Num("n", "-"), Num("delta", "-", 0.75)}, &LllBound);
	Def(R, "odefect", Lat, "ch07 Lattices: bench 3, orthogonality defect", "orthogonality defect prod|b_i| / vol(L)", "defect (1 = orthogonal)",
		{Mat("B", "rows = basis vectors")}, &ODefect);
	Def(R, "ellcov", Lat, "ch07 Lattices: bench 2, integer ambiguity resolution", "ambiguity covariance of an ellipse with semi-axes a, b tilted by theta",
		"2 x 2 covariance, cycles^2", {Num("a", "cycles", 1.2), Num("b", "cycles", 0.07), Num("theta", "rad", 52.0 * cst::Pi / 180.0)}, &EllCov);
	Def(R, "zdecor", Lat, "ch07 Lattices: bench 2, decorrelation", "unimodular Z making Z' Q Z Lagrange-reduced (2-D LAMBDA decorrelation)", "2 x 2 integer Z",
		{Mat("Q", "cycles^2")}, &ZDecor);
	Def(R, "pboot", Lat, "ch07 Lattices: Three ways to guess an integer", "bootstrapping success prod(2 Phi(1/(2 sigma_i|I)) - 1)", "probability",
		{Mat("Q", "cycles^2")}, &PBoot);
	Def(R, "pround", Lat, "ch07 Lattices: Three ways to guess an integer", "rounding success P(round(a_hat) = a), exact for n <= 2", "probability",
		{Mat("Q", "cycles^2")}, &PRound);
	Def(R, "weyl", Lat, "ch07 Lattices: connections, Weyl equidistribution", "star discrepancy and worst bin of {k alpha}, k = 1..n", "record dstar, dev",
		{Num("alpha", "-"), Num("n", "-", 100000), Num("bins", "-", 20)}, &Weyl);
	Def(R, "unwrap", Lat, "ch07 Lattices: connections, phase on the circle", "phase unwrapping: add the whole turns that keep steps within (-pi, pi]", "unwrapped phases, rad",
		{Mat("phi", "rad")}, &Unwrap);

	// Chapter 10: category theory.
	Def(R, "kcomp", Cat, "ch10 Category: Kernels compose (Chapman-Kolmogorov)", "composite g after f of column-stochastic kernels", "kernel g*f",
		{Mat("g", "P(z|y)"), Mat("f", "P(y|x)")}, &KComp);
	Def(R, "kpow", Cat, "ch10 Category: Kernels compose (Chapman-Kolmogorov)", "n-step kernel P^n", "kernel", {Mat("P", "P(y|x)"), Num("n", "-")}, &KPow);
	Def(R, "bayesinv", Cat, "ch10 Category: Kernels compose, and so do their inverses", "Bayesian inverse of kernel K under prior p: entries P(x | y)",
		"kernel, cols(K) x rows(K)", {Mat("K", "P(y|x)"), Mat("p", "prior on x")}, &BayesInv);
	Def(R, "ringkern", Cat, "ch10 Category: bench 2, Markov kernels", "the bench's motion kernel on a ring of n cells with slip s", "n x n kernel",
		{Num("n", "-", 8), Num("s", "-", 0.2)}, &RingKern);
	Def(R, "tf", Cat, "ch10 Category: Wiring is a functor", "transfer function num(s)/den(s) as a 2-row [num; den]", "2 x L matrix",
		{Mat("num", "coefficients, s^n first"), Mat("den", "coefficients, s^n first")}, &Tf);
	Def(R, "tfser", Cat, "ch10 Category: Wiring is a functor", "series connection G2 G1", "[num; den]", {Mat("G1", "[num; den]"), Mat("G2", "[num; den]")}, &TfSeries);
	Def(R, "tfpar", Cat, "ch10 Category: Wiring is a functor", "parallel connection G1 + G2", "[num; den]", {Mat("G1", "[num; den]"), Mat("G2", "[num; den]")}, &TfParallel);
	Def(R, "tffb", Cat, "ch10 Category: Wiring is a functor", "negative feedback G / (1 + G H) (H = 1: unity)", "[num; den]",
		{Mat("G", "[num; den]"), MatDefault("H", "[num; den] or gain", 1.0)}, &TfFeedback);
	Def(R, "tfjw", Cat, "ch10 Category: Wiring is a functor", "frequency response G(jw)", "record re, im, mag, db, phase (rad)", {Mat("G", "[num; den]"), Num("w", "rad/s")},
		&TfJw);
	Def(R, "ssjw", Cat, "ch10 Category: Wiring is a functor (realizations)", "SISO state-space response C (jw I - A)^-1 B + D", "record re, im, mag, db, phase (rad)",
		{Mat("A", "1/s"), Mat("B", "-"), Mat("C", "-"), Mat("D", "-"), Num("w", "rad/s")}, &SsJw);
	Def(R, "sheaflap", Cat, "ch10 Category: Sensors that disagree about frames", "sheaf Laplacian of a sensor cycle with rotation restriction maps", "2n x 2n matrix",
		{Mat("frames", "rad"), Num("err", "rad", 0.0)}, &SheafLap);
	Def(R, "sheafh0", Cat, "ch10 Category: Sensors that disagree about frames", "dimension of global sections and slowest diffusion rate", "record h0, lmin",
		{Mat("frames", "rad"), Num("err", "rad", 0.0)}, &SheafH0);

	// Chapter 11: procedural systems.
	Def(R, "ca", Pro, "ch11 Procedural: 256 rules, one row at a time", "elementary automaton from one centre cell, periodic, after n steps", "1 x width row of bits",
		{Num("rule", "0..255"), Num("steps", "-"), Num("width", "cells", 201)}, &Ca);
	Def(R, "cacenter", Pro, "ch11 Procedural: 256 rules, one row at a time", "centre-column bits t = 0..steps (infinite lattice)", "column of bits",
		{Num("rule", "0..255"), Num("steps", "-")}, &CaCenter);
	Def(R, "cafreq", Pro, "ch11 Procedural: 256 rules, one row at a time", "1-frequency of the centre column over t = 1..steps", "fraction",
		{Num("rule", "0..255"), Num("steps", "-", 1999)}, &CaFreq);
	Def(R, "moran", Pro, "ch11 Procedural: Rules, attractors, and patterns", "similarity dimension: solve sum r_i^d = 1 (Moran)", "dimension d",
		{Mat("r", "ratios")}, &Moran);
	Def(R, "ifspreset", Pro, "ch11 Procedural: Four affine maps make a fern", "bench presets: 1 fern, 2 Sierpinski, 3 tree, 4 expanding map, 5 expanding on average",
		"n x 7 rows [a b c d e f p]", {Num("k", "1..5")}, &IfsPreset);
	Def(R, "ifsnorm", Pro, "ch11 Procedural: Four affine maps make a fern", "contraction factor ||A_i|| (operator 2-norm) of each map", "column of factors",
		{Mat("maps", "rows [a b c d e f (p)]")}, &IfsNorm);
	Def(R, "ifslog", Pro, "ch11 Procedural: Four affine maps make a fern", "mean log factor sum p_i ln||A_i|| (< 0: contracts on average)", "nats",
		{Mat("maps", "rows [a b c d e f (p)]")}, &IfsLog);
	Def(R, "turing", Pro, "ch11 Procedural: Diffusion that makes stripes", "Schnakenberg Turing band, fastest wavenumber and its wavelength",
		"record kmin, kmax, kstar (1/length), wl (length), gmax (1/time)", {Num("Dv", "-", 40), Num("Du", "-", 1), Num("a", "-", 0.1), Num("b", "-", 0.9)}, &Turing);
	Def(R, "lsys", Pro, "ch11 Procedural: Eleven rules and what they are", "Fibonacci L-system a -> ab, b -> a after n rewrites of 'a'", "record len, a, b, ratio",
		{Num("n", "-")}, &LSys);
	Def(R, "spline", Pro, "ch11 Procedural: Eleven rules (natural cubic spline)", "natural cubic spline through (xs, ys), evaluated at x", "values, same shape as x",
		{Mat("xs", "-"), Mat("ys", "-"), Mat("x", "-")}, &Spline);
	Def(R, "splinee", Pro, "ch11 Procedural: Eleven rules (natural cubic spline)", "bending energy integral of f''^2 of the natural spline", "energy",
		{Mat("xs", "-"), Mat("ys", "-")}, &SplineEnergy);
	Def(R, "lloyd", Pro, "ch11 Procedural: Spreading a drone team", "the bench's Lloyd iteration: 12 drones, hot spot at (hx, 0.3)", "record e0, e, gap",
		{Num("hx", "-", 0.7), Num("steps", "-", 60)}, &Lloyd);

	// Chapter 12: expansion fields.
	Def(R, "pcrb", Exp, "ch12 Expansion: A bound no filter can beat", "posterior Cramer-Rao bound after k steps (linear-Gaussian: = Kalman covariance)",
		"n x n bound", {Mat("F", "-"), Mat("Q", "-"), Mat("H", "-"), Mat("R", "-"), Mat("P0", "-"), Num("k", "steps", 1)}, &Pcrb);
	Def(R, "pcrbcv", Exp, "ch12 Expansion: A bound no filter can beat", "PCRB of the constant-velocity tracker (P0 = diag(4, 1), position measured)",
		"record p11 (m^2), p12, p22", {Num("k", "steps", 40), Num("q", "m^2/s^3", 0.05), Num("r", "m^2", 0.5), Num("dt", "s", 1)}, &PcrbCv);
	Def(R, "bangbang", Exp, "ch12 Expansion: The fastest way to stop", "minimum time to rest for x'' = u, |u| <= 1, and the switch point",
		"record T (s), u0, ts (s), xs, vs", {Num("x", "m"), Num("v", "m/s")}, &BangBang);
	Def(R, "bbsim", Exp, "ch12 Expansion: The fastest way to stop", "simulated bang-bang arrival time (bench integrator)", "time, s",
		{Num("x", "m"), Num("v", "m/s"), Num("dt", "s", 1e-4)}, &BbSim);
	Def(R, "ouproc", Exp, "ch12 Expansion: How uncertainty spreads between measurements", "Ornstein-Uhlenbeck law dX = -theta X dt + sigma dW at time t",
		"record mean, var, sd, vinf", {Num("x0", "-", 2), Num("theta", "1/s", 0.8), Num("sigma", "1/sqrt(s)", 0.6), Num("t", "s", 2)}, &OuProc);
	Def(R, "qcv", Exp, "ch12 Expansion: Ito integration", "white-acceleration process noise q [dt^3/3, dt^2/2; dt^2/2, dt]", "2 x 2", {Num("q", "m^2/s^3"), Num("dt", "s", 1)},
		&QCv);
	Def(R, "qdisc", Exp, "ch12 Expansion: Ito integration (Van Loan)", "discrete process noise integral of e^(A s) Qc e^(A' s) ds over dt", "n x n",
		{Mat("A", "1/s"), Mat("Qc", "-"), Num("dt", "s")}, &QDisc);
	Def(R, "tropmul", Exp, "ch12 Expansion: tropical algebra", "(min, +) matrix product", "matrix", {Mat("A", "-"), Mat("B", "-")}, &TropMul);
	Def(R, "troppow", Exp, "ch12 Expansion: tropical algebra", "(min, +) power W^k; default n - 1 gives all-pairs shortest paths", "matrix",
		{Mat("W", "inf = no edge"), Num("k", "-", -1)}, &TropPow);
	Def(R, "sinkhorn", Exp, "ch12 Expansion: optimal transport", "entropic transport plan between equal swarms (unit marginals)", "n x n plan",
		{Mat("C", "cost"), Num("eps", "cost"), Num("iters", "-", 20000)}, &Sinkhorn);
	Def(R, "sinkcost", Exp, "ch12 Expansion: optimal transport", "cost sum(P .* C) of the Sinkhorn plan", "cost", {Mat("C", "cost"), Num("eps", "cost"), Num("iters", "-", 20000)},
		&SinkCost);
	Def(R, "otexact", Exp, "ch12 Expansion: optimal transport", "optimal assignment cost (the eps -> 0 limit)", "cost", {Mat("C", "cost")}, &OtExact);
	Def(R, "kepenergy", Exp, "ch12 Expansion: symplectic integration", "worst Kepler energy error: leapfrog vs explicit Euler (stops past 1)", "record lf, eu",
		{Num("orbits", "-", 200), Num("steps", "per orbit", 200)}, &KepEnergy);
	Def(R, "mpedges", Exp, "ch12 Expansion: random matrices", "Marchenko-Pastur support s2 (1 -/+ sqrt c)^2", "record lo, hi", {Num("c", "p/n"), Num("s2", "-", 1)}, &MpEdges);
	Def(R, "mppdf", Exp, "ch12 Expansion: random matrices", "Marchenko-Pastur density (continuous part)", "density, same shape as x",
		{Mat("x", "-"), Num("c", "p/n"), Num("s2", "-", 1)}, &MpPdf);
}
} // namespace mcalc
