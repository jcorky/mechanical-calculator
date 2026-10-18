// Chapters 7, 10, 11 and 12 of the atlas: lattices and GNSS ranging, category theory (kernels, transfer
// functions, sheaves), procedural systems (automata, IFS, Turing, splines, Lloyd) and expansion fields
// (Cramer-Rao, HJB, Fokker-Planck, Ito, tropical algebra, transport, symplectic steps, random matrices).
// All arithmetic is gear arithmetic (Real); integers stay integers.
#include "CalcLibrary.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace mcalc
{
namespace fields_detail
{
Param MatDefault(const char* Name, const char* Unit, Real Default)
{
	Param P = Mat(Name, Unit);
	P.HasDefault = true;
	P.Default = Default;
	return P;
}

const Real& InvSqrt2()
{
	static const Real V = Real(1) / gm::sqrt(Real(2));
	return V;
}

Real NormalCdf(const Real& X) { return Real(0.5) * gm::erfc(-X * InvSqrt2()); }

// Phi(B) - Phi(A) without cancellation in either tail.
Real NormalBand(const Real& A, const Real& B)
{
	if (A > 0)
	{
		return Real(0.5) * (gm::erfc(A * InvSqrt2()) - gm::erfc(B * InvSqrt2()));
	}
	return Real(0.5) * (gm::erfc(-B * InvSqrt2()) - gm::erfc(-A * InvSqrt2()));
}

bool IsWhole(const Real& X, long long Lo, long long Hi) { return X.IsFinite() && X == gm::floor(X) && X >= Real(Lo) && X <= Real(Hi); }

// Round half to even, as Python's round() in the verification scripts.
Real RoundHalfEven(const Real& X)
{
	const Real F = gm::floor(X), D = X - F, Half(0.5);
	if (D > Half)
	{
		return F + 1;
	}
	if (D < Half)
	{
		return F;
	}
	return gm::fmod(F, Real(2)).IsZero() ? F : F + 1;
}

std::string ShapeText(const Matrix& M) { return std::to_string(M.Rows) + " x " + std::to_string(M.Cols); }

Real RowDot(const Matrix& A, int I, const Matrix& B, int J)
{
	Real S = 0;
	for (int c = 0; c < A.Cols; ++c)
	{
		S += A(I, c) * B(J, c);
	}
	return S;
}

// mulberry32, as in the atlas pages, in exact integer arithmetic; its output is entered as the decimal t / 2^32.
struct FieldsMulberry
{
	uint32_t State;
	explicit FieldsMulberry(uint32_t Seed) : State(Seed) {}
	Real Next()
	{
		State += 0x6D2B79F5u;
		uint32_t T = State;
		T = (T ^ (T >> 15)) * (T | 1u);
		T ^= T + (T ^ (T >> 7)) * (T | 61u);
		return Real(static_cast<unsigned long long>(T ^ (T >> 14))) / Real(4294967296ull);
	}
};

struct FieldsComplex
{
	Real Re;
	Real Im;
};

FieldsComplex ComplexDiv(const FieldsComplex& A, const FieldsComplex& B)
{
	const Real D = B.Re * B.Re + B.Im * B.Im;
	return {(A.Re * B.Re + A.Im * B.Im) / D, (A.Im * B.Re - A.Re * B.Im) / D};
}

// ---------------------------------------------------------------- ch07: time of flight

Result TofRange(const Call& C)
{
	if (!(C.N(1) > 0))
	{
		return Result::Fail("way must be positive: 1 for one-way, 2 for two-way (radar)");
	}
	return Result::Of(cst::C() * C.N(0) / C.N(1));
}

Result TofRes(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("bandwidth must be positive");
	}
	return Result::Of(cst::C() / (Real(2) * C.N(0)));
}

Result Widelane(const Call& C)
{
	const Real F1 = C.N(0), F2 = C.N(1), Light = cst::C();
	if (!(F1 > F2 && F2 > 0))
	{
		return Result::Fail("needs f1 > f2 > 0");
	}
	return Result::Record({{"l1", Light / F1}, {"l2", Light / F2}, {"wide", Light / (F1 - F2)}, {"narrow", Light / (F1 + F2)}});
}

Result SatRange(const Call& C)
{
	const Real El = C.N(0), Re = C.N(1), Rs = C.N(2);
	if (!(Re > 0 && Rs > Re))
	{
		return Result::Fail("needs 0 < re < rs");
	}
	const Real Rc = Re * gm::cos(El);
	const Real D = gm::sqrt(Rs * Rs - Rc * Rc) - Re * gm::sin(El);
	return Result::Record({{"d", D}, {"t", D / cst::C()}});
}

// ---------------------------------------------------------------- ch07: pseudorange positioning

bool ReadSatellites(const Call& C, Matrix* S, std::vector<Real>* Rho, std::string* Error)
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
bool GeometryAt(const Matrix& S, const Real* X, Matrix* G, std::vector<Real>* Dist)
{
	*G = Matrix(S.Rows, 4);
	Dist->assign(static_cast<size_t>(S.Rows), Real(0));
	for (int i = 0; i < S.Rows; ++i)
	{
		const Real Dx = S(i, 0) - X[0], Dy = S(i, 1) - X[1], Dz = S(i, 2) - X[2];
		const Real D = gm::sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
		if (!(D > 0))
		{
			return false;
		}
		(*G)(i, 0) = -Dx / D;
		(*G)(i, 1) = -Dy / D;
		(*G)(i, 2) = -Dz / D;
		(*G)(i, 3) = 1;
		(*Dist)[static_cast<size_t>(i)] = D;
	}
	return true;
}

// Gauss-Newton on rho_i = |s_i - x| + c b from the Earth's centre, as in the chapter's bench.
bool PseudorangeFix(const Matrix& S, const std::vector<Real>& Rho, Real* X, int* Iters, std::string* Error)
{
	for (int a = 0; a < 4; ++a)
	{
		X[a] = 0;
	}
	const Real Stop(1e-6);
	Matrix G;
	std::vector<Real> Dist;
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
		Real Step = 0;
		for (int a = 0; a < 4; ++a)
		{
			X[a] += Dx.V[static_cast<size_t>(a)];
			Step += Dx.V[static_cast<size_t>(a)] * Dx.V[static_cast<size_t>(a)];
		}
		if (gm::sqrt(Step) < Stop)
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
	std::vector<Real> Rho;
	std::string E;
	if (!ReadSatellites(C, &S, &Rho, &E))
	{
		return Result::Fail(E);
	}
	Real X[4];
	int Iters = 0;
	if (!PseudorangeFix(S, Rho, X, &Iters, &E))
	{
		return Result::Fail(E);
	}
	return Result::Record({{"x", X[0]}, {"y", X[1]}, {"z", X[2]}, {"cb", X[3]}, {"dt", X[3] / cst::C()}, {"it", Real(Iters)}});
}

// DOPs from the geometry matrix; Up is the local vertical (unit) used to split vertical from horizontal.
Result DopRecord(const Matrix& G, const Real* Up)
{
	Matrix Q;
	std::string E;
	if (!la::Inverse(la::Mul(la::Transpose(G), G), &Q, &E))
	{
		return Result::Fail("geometry is singular (3 satellites cannot separate clock and position): " + E);
	}
	const Real P2 = Q(0, 0) + Q(1, 1) + Q(2, 2);
	Real V2 = 0;
	for (int a = 0; a < 3; ++a)
	{
		for (int b = 0; b < 3; ++b)
		{
			V2 += Up[a] * Q(a, b) * Up[b];
		}
	}
	return Result::Record({{"gdop", gm::sqrt(P2 + Q(3, 3))}, {"pdop", gm::sqrt(P2)}, {"hdop", gm::sqrt(gm::max(Real(0), P2 - V2))}, {"vdop", gm::sqrt(V2)},
		{"tdop", gm::sqrt(Q(3, 3))}});
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
	const Real R = gm::sqrt(P.V[0] * P.V[0] + P.V[1] * P.V[1] + P.V[2] * P.V[2]);
	if (!(R > 0))
	{
		return Result::Fail("receiver at the origin has no local vertical");
	}
	const Real Up[3] = {P.V[0] / R, P.V[1] / R, P.V[2] / R};
	Matrix G;
	std::vector<Real> Dist;
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
		const Real A = Az.V[static_cast<size_t>(i)], L = El.V[static_cast<size_t>(i)], CosL = gm::cos(L);
		G(i, 0) = -CosL * gm::sin(A);
		G(i, 1) = -CosL * gm::cos(A);
		G(i, 2) = -gm::sin(L);
		G(i, 3) = 1;
	}
	const Real Up[3] = {Real(0), Real(0), Real(1)};
	return DopRecord(G, Up);
}

// Baarda's w-statistic: residual over sigma sqrt(1 - h_ii), h the hat matrix of the fix.
Result Baarda(const Call& C)
{
	Matrix S;
	std::vector<Real> Rho;
	std::string E;
	if (!ReadSatellites(C, &S, &Rho, &E))
	{
		return Result::Fail(E);
	}
	const Real Sigma = C.N(2);
	if (!(Sigma > 0))
	{
		return Result::Fail("sigma must be positive (m)");
	}
	Real X[4];
	int Iters = 0;
	if (!PseudorangeFix(S, Rho, X, &Iters, &E))
	{
		return Result::Fail(E);
	}
	Matrix G, Q;
	std::vector<Real> Dist;
	if (!GeometryAt(S, X, &G, &Dist) || !la::Inverse(la::Mul(la::Transpose(G), G), &Q, &E))
	{
		return Result::Fail("geometry is singular");
	}
	const Real MinRedundancy(1e-12);
	Matrix W(S.Rows, 1);
	for (int i = 0; i < S.Rows; ++i)
	{
		Real H = 0;
		for (int a = 0; a < 4; ++a)
		{
			for (int b = 0; b < 4; ++b)
			{
				H += G(i, a) * Q(a, b) * G(i, b);
			}
		}
		const Real Res = Rho[static_cast<size_t>(i)] - Dist[static_cast<size_t>(i)] - X[3];
		W.V[static_cast<size_t>(i)] = Res / (Sigma * gm::sqrt(gm::max(MinRedundancy, Real(1) - H)));
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
	Real B1[2] = {B(0, 0), B(0, 1)}, B2[2] = {B(1, 0), B(1, 1)};
	if ((B1[0] * B2[1] - B1[1] * B2[0]).IsZero())
	{
		return Result::Fail("basis vectors are parallel");
	}
	for (int Step = 0; Step < 10000; ++Step)
	{
		const Real N1 = B1[0] * B1[0] + B1[1] * B1[1], N2 = B2[0] * B2[0] + B2[1] * B2[1];
		if (N2 < N1)
		{
			std::swap(B1[0], B2[0]);
			std::swap(B1[1], B2[1]);
			continue;
		}
		const Real Q = RoundHalfEven((B1[0] * B2[0] + B1[1] * B2[1]) / N1);
		if (Q.IsZero())
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
			const Real M = RowDot(B, i, *Bs, j) / RowDot(*Bs, j, *Bs, j);
			(*Mu)(i, j) = M;
			for (int c = 0; c < B.Cols; ++c)
			{
				(*Bs)(i, c) -= M * (*Bs)(j, c);
			}
		}
	}
}

bool LllReduce(Matrix* B, const Real& Delta, std::string* Error)
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
			const Real Q = RoundHalfEven(Mu(K, j));
			if (!Q.IsZero())
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

bool DeltaOk(const Real& Delta) { return Delta > Real(0.25) && Delta <= Real(1); }

Result Lll(const Call& C)
{
	Matrix B = C.M(0);
	const Real Delta = C.N(1);
	if (!DeltaOk(Delta))
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
	std::vector<Real> Norms;
	std::vector<Real> Coef;
	Real Best;
	Real Radius;
	Real Tiny;
	Real Slack;
	long Visits = 0;
};

// Fincke-Pohst depth-first enumeration in Gram-Schmidt coordinates.
void EnumerateShort(FieldsEnum& St, int I, const Real& Partial)
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
			St.Radius = Partial * St.Slack;
		}
		return;
	}
	Real Centre = 0;
	for (int j = I + 1; j < St.N; ++j)
	{
		Centre -= St.Coef[static_cast<size_t>(j)] * St.Mu(j, I);
	}
	const Real Nrm = St.Norms[static_cast<size_t>(I)];
	const Real R = gm::sqrt(gm::max(St.Radius - Partial, Real(0)) / Nrm);
	const long long Lo = gm::ceil(Centre - R).ToLong(), Hi = gm::floor(Centre + R).ToLong();
	for (long long Z = Lo; Z <= Hi; ++Z)
	{
		const Real Zr(Z), Off = Zr - Centre;
		St.Coef[static_cast<size_t>(I)] = Zr;
		const Real P = Partial + Off * Off * Nrm;
		if (P < St.Radius)
		{
			EnumerateShort(St, I - 1, P);
		}
	}
	St.Coef[static_cast<size_t>(I)] = 0;
}

Result ShortVec(const Call& C)
{
	Matrix B = C.M(0);
	if (B.Rows < 1 || B.Rows > B.Cols || B.Rows > 10)
	{
		return Result::Fail("basis is n x m with one vector per row, n <= m, n <= 10; got " + ShapeText(B));
	}
	std::string E;
	if (!LllReduce(&B, Real(0.75), &E))
	{
		return Result::Fail(E);
	}
	FieldsEnum St;
	St.N = B.Rows;
	Matrix Bs;
	GramSchmidt(B, &Bs, &St.Mu);
	St.Best = Real::Infinity();
	for (int i = 0; i < B.Rows; ++i)
	{
		St.Norms.push_back(RowDot(Bs, i, Bs, i));
		St.Best = gm::min(St.Best, RowDot(B, i, B, i));
	}
	St.Slack = Real(1) + Real(1e-9);
	St.Radius = St.Best * St.Slack;
	St.Tiny = St.Best * Real(1e-12);
	St.Coef.assign(static_cast<size_t>(St.N), Real(0));
	EnumerateShort(St, St.N - 1, Real(0));
	return Result::Of(gm::sqrt(St.Best));
}

Result LllBound(const Call& C)
{
	const Real N = C.N(0), Delta = C.N(1);
	if (!IsWhole(N, 1, 1000) || !DeltaOk(Delta))
	{
		return Result::Fail("needs a whole dimension n >= 1 and delta in (1/4, 1]");
	}
	return Result::Of(gm::pow(Real(1) / (Delta - Real(0.25)), (N - 1) / 2));
}

Result ODefect(const Call& C)
{
	const Matrix& B = C.M(0);
	if (B.Rows < 1 || B.Rows > B.Cols)
	{
		return Result::Fail("basis is n x m with one vector per row, n <= m");
	}
	Real Det = 0;
	std::string E;
	if (!la::Determinant(la::Mul(B, la::Transpose(B)), &Det, &E) || !(Det > 0))
	{
		return Result::Fail("basis vectors (rows) must be linearly independent");
	}
	Real P = 1;
	for (int i = 0; i < B.Rows; ++i)
	{
		P *= gm::sqrt(RowDot(B, i, B, i));
	}
	return Result::Of(P / gm::sqrt(Det));
}

// ---------------------------------------------------------------- ch07: integer ambiguities

Result EllCov(const Call& C)
{
	const Real A = C.N(0), B = C.N(1), T = C.N(2), Cs = gm::cos(T), Sn = gm::sin(T);
	if (!(A > 0 && B > 0))
	{
		return Result::Fail("semi-axes must be positive (cycles)");
	}
	const Real A2 = A * A, B2 = B * B, Off = Cs * Sn * (A2 - B2);
	return Result::Of(Matrix::FromRows({{Cs * Cs * A2 + Sn * Sn * B2, Off}, {Off, Sn * Sn * A2 + Cs * Cs * B2}}));
}

bool ReadCovariance(const Matrix& Q, int MaxN, std::string* Error)
{
	if (!Q.IsSquare() || Q.Rows > MaxN)
	{
		*Error = "covariance must be square, at most " + std::to_string(MaxN) + " x " + std::to_string(MaxN) + "; got " + ShapeText(Q);
		return false;
	}
	const Real Tol(1e-12);
	for (int i = 0; i < Q.Rows; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			if (gm::fabs(Q(i, j) - Q(j, i)) > Tol * (gm::fabs(Q(i, i)) + gm::fabs(Q(j, j))))
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
	Real Z[2][2] = {{Real(1), Real(0)}, {Real(0), Real(1)}};
	Real Q00 = Q(0, 0), Q01 = Q(0, 1), Q11 = Q(1, 1);
	if (!(Q00 > 0 && Q11 > 0 && Q00 * Q11 > Q01 * Q01))
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
		const Real Mu = RoundHalfEven(Q01 / Q00);
		if (Mu.IsZero())
		{
			return Result::Of(Matrix::FromRows({{Z[0][0], Z[0][1]}, {Z[1][0], Z[1][1]}}));
		}
		Z[0][1] -= Mu * Z[0][0];
		Z[1][1] -= Mu * Z[1][0];
		Q11 = Q11 - Real(2) * Mu * Q01 + Mu * Mu * Q00;
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
	Real P = 1;
	for (int i = 0; i < Q.Rows; ++i)
	{
		P *= Real(2) * NormalCdf(Real(1) / (Real(2) * L(i, i))) - 1;
	}
	return Result::Of(P);
}

struct RoundIntegrand
{
	Real S1, Slope, S2, Norm, HalfInvVar;
	Real operator()(const Real& X) const
	{
		const Real Dens = gm::exp(-X * X * HalfInvVar) * Norm;
		const Real SX = Slope * X;
		return Dens * NormalBand((Real(-0.5) - SX) / S2, (Real(0.5) - SX) / S2);
	}
};

// Gauss-Legendre nodes and weights of order 10 on [-1, 1], found on the gears by Newton's method on P_10.
struct GaussRule
{
	Real X[5];
	Real W[5];
};

const GaussRule& Legendre10()
{
	static const GaussRule Rule = [] {
		GaussRule G;
		const int N = 10;
		for (int k = 0; k < 5; ++k)
		{
			Real X = gm::cos(gm::Pi() * (Real(4 * k + 3)) / Real(4 * N + 2));
			Real Dp = 0;
			for (int It = 0; It < 100; ++It)
			{
				Real P0 = 1, P1 = X;
				for (int n = 2; n <= N; ++n)
				{
					const Real P2 = (Real(2 * n - 1) * X * P1 - Real(n - 1) * P0) / Real(n);
					P0 = P1;
					P1 = P2;
				}
				Dp = Real(N) * (X * P1 - P0) / (X * X - 1);
				const Real Step = P1 / Dp;
				X -= Step;
				if (gm::fabs(Step) < Real(1e-19))
				{
					break;
				}
			}
			G.X[k] = X;
			G.W[k] = Real(2) / ((Real(1) - X * X) * Dp * Dp);
		}
		return G;
	}();
	return Rule;
}

Real GaussPanel(const RoundIntegrand& F, const Real& A, const Real& B)
{
	const GaussRule& G = Legendre10();
	const Real Mid = (A + B) / 2, Half = (B - A) / 2;
	RealAccumulator S;
	for (int k = 0; k < 5; ++k)
	{
		S.AddProduct(G.W[k], F(Mid + Half * G.X[k]) + F(Mid - Half * G.X[k]));
	}
	return S.Result() * Half;
}

Real AdaptiveGauss(const RoundIntegrand& F, const Real& A, const Real& B, const Real& Whole, const Real& Tol, int Depth)
{
	const Real M = (A + B) / 2;
	const Real Left = GaussPanel(F, A, M), Right = GaussPanel(F, M, B);
	if (Depth <= 0 || gm::fabs(Left + Right - Whole) <= Tol)
	{
		return Left + Right;
	}
	return AdaptiveGauss(F, A, M, Left, Tol / 2, Depth - 1) + AdaptiveGauss(F, M, B, Right, Tol / 2, Depth - 1);
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
	if (!(Q(0, 0) > 0))
	{
		return Result::Fail("covariance must be positive definite");
	}
	const Real S1 = gm::sqrt(Q(0, 0));
	if (Q.Rows == 1)
	{
		return Result::Of(Real(2) * NormalCdf(Real(1) / (Real(2) * S1)) - 1);
	}
	const Real Slope = Q(1, 0) / Q(0, 0), D2 = Q(1, 1) - Q(1, 0) * Q(1, 0) / Q(0, 0);
	if (!(D2 > 0))
	{
		return Result::Fail("covariance must be positive definite");
	}
	const RoundIntegrand F{S1, Slope, gm::sqrt(D2), Real(1) / (S1 * gm::sqrt(Real(2) * gm::Pi())), Real(1) / (Real(2) * Q(0, 0))};
	std::vector<Real> Cuts;
	for (int k = 0; k <= 8; ++k)
	{
		Cuts.push_back(Real(k - 4) / Real(8));
	}
	if (!Slope.IsZero())
	{
		const Real Edge = gm::fabs(Real(0.5) / Slope);
		if (Edge < Real(0.5))
		{
			Cuts.push_back(Edge);
			Cuts.push_back(-Edge);
		}
	}
	std::sort(Cuts.begin(), Cuts.end());
	const Real Tol(1e-13);
	Real Total = 0;
	for (size_t k = 0; k + 1 < Cuts.size(); ++k)
	{
		const Real A = Cuts[k], B = Cuts[k + 1];
		if (!(B > A))
		{
			continue;
		}
		Total += AdaptiveGauss(F, A, B, GaussPanel(F, A, B), Tol, 30);
	}
	return Result::Of(Total);
}

// ---------------------------------------------------------------- ch07: the circle

Result Weyl(const Call& C)
{
	const Real Alpha = C.N(0);
	if (!IsWhole(C.N(1), 1, 500000) || !IsWhole(C.N(2), 1, 10000))
	{
		return Result::Fail("n must be a whole number 1..500000 and bins 1..10000");
	}
	const int Count = C.Int(1), Nb = C.Int(2);
	const Real Rb(Nb), Rn(Count);
	std::vector<Real> Pts(static_cast<size_t>(Count));
	std::vector<int> Hist(static_cast<size_t>(Nb), 0);
	Real Y = 0;
	for (int k = 1; k <= Count; ++k)
	{
		Y = Alpha * Real(k);
		const Real F = Y - gm::floor(Y);
		Pts[static_cast<size_t>(k - 1)] = F;
		Hist[static_cast<size_t>(std::min(Nb - 1, (F * Rb).ToInt()))] += 1;
	}
	std::sort(Pts.begin(), Pts.end());
	Real Star = 0;
	for (int i = 0; i < Count; ++i)
	{
		const Real X = Pts[static_cast<size_t>(i)];
		Star = gm::max(Star, gm::max(Real(i + 1) / Rn - X, X - Real(i) / Rn));
	}
	int Most = 0, Least = Count;
	for (int H : Hist)
	{
		Most = std::max(Most, H);
		Least = std::min(Least, H);
	}
	const Real Ideal = Real(1) / Rb;
	const Real Dev = gm::max(gm::fabs(Real(Most) / Rn - Ideal), gm::fabs(Real(Least) / Rn - Ideal));
	return Result::Record({{"dstar", Star}, {"dev", Dev}});
}

Result Unwrap(const Call& C)
{
	Matrix M = C.M(0);
	if (!M.IsVector())
	{
		return Result::Fail("phases must be a vector (rad)");
	}
	const Real Tau = cst::Tau(), Pi = cst::Pi();
	Real Prev = M.V[0], Acc = M.V[0];
	for (size_t i = 1; i < M.V.size(); ++i)
	{
		const Real D = M.V[i] - Prev;
		Prev = M.V[i];
		Acc += D - Tau * gm::floor((D + Pi) / Tau);
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
	const Real NegTol(-1e-12), SumTol(1e-9);
	for (int j = 0; j < K.Cols; ++j)
	{
		Real S = 0;
		for (int i = 0; i < K.Rows; ++i)
		{
			if (K(i, j) < NegTol)
			{
				*Error = std::string(Name) + " has a negative entry: not a Markov kernel";
				return false;
			}
			S += K(i, j);
		}
		if (gm::fabs(S - 1) > SumTol)
		{
			*Error = std::string(Name) + " column " + std::to_string(j + 1) + " sums to " + S.ToString(10) + ": columns of a kernel K[y, x] = P(y | x) must sum to 1";
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
	if (!IsWhole(C.N(1), 0, 1000000000LL))
	{
		return Result::Fail("n must be a whole number >= 0");
	}
	long long N = C.N(1).ToLong();
	Matrix R = Matrix::Identity(P.Rows), B = P;
	while (N > 0)
	{
		if (N & 1)
		{
			R = la::Mul(R, B);
		}
		N >>= 1;
		if (N > 0)
		{
			B = la::Mul(B, B);
		}
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
	Real Sum = 0;
	for (const Real& X : P.V)
	{
		if (X.IsNegative())
		{
			return Result::Fail("prior must be nonnegative");
		}
		Sum += X;
	}
	if (gm::fabs(Sum - 1) > Real(1e-9))
	{
		return Result::Fail("prior must sum to 1");
	}
	Matrix Out(K.Cols, K.Rows);
	for (int y = 0; y < K.Rows; ++y)
	{
		Real Z = 0;
		for (int x = 0; x < K.Cols; ++x)
		{
			Z += K(y, x) * P.V[static_cast<size_t>(x)];
		}
		for (int x = 0; x < K.Cols; ++x)
		{
			Out(x, y) = Z > 0 ? K(y, x) * P.V[static_cast<size_t>(x)] / Z : Real(0);
		}
	}
	return Result::Of(Out);
}

// The chapter's motion kernel on a ring: one cell forward with 1 - 2s, two with s, stay with s.
Result RingKern(const Call& C)
{
	const Real S = C.N(1);
	if (!IsWhole(C.N(0), 1, 500) || !(S >= 0 && S <= Real(0.5)))
	{
		return Result::Fail("needs a whole n in 1..500 and slip 0 <= s <= 0.5");
	}
	const int N = C.Int(0);
	const Real Ahead = Real(1) - Real(2) * S;
	Matrix F(N, N);
	for (int x = 0; x < N; ++x)
	{
		F((x + 1) % N, x) += Ahead;
		F((x + 2) % N, x) += S;
		F(x, x) += S;
	}
	return Result::Of(F);
}

// ---------------------------------------------------------------- ch10: transfer functions

std::vector<Real> TrimPoly(std::vector<Real> P)
{
	size_t K = 0;
	while (K + 1 < P.size() && P[K].IsZero())
	{
		++K;
	}
	P.erase(P.begin(), P.begin() + static_cast<long>(K));
	if (P.empty())
	{
		P.push_back(Real(0));
	}
	return P;
}

bool IsZeroPoly(const std::vector<Real>& P) { return P.size() == 1 && P[0].IsZero(); }

std::vector<Real> PolyTimes(const std::vector<Real>& A, const std::vector<Real>& B)
{
	std::vector<Real> Out(A.size() + B.size() - 1, Real(0));
	for (size_t i = 0; i < A.size(); ++i)
	{
		for (size_t j = 0; j < B.size(); ++j)
		{
			Out[i + j] += A[i] * B[j];
		}
	}
	return Out;
}

std::vector<Real> PolyPlus(const std::vector<Real>& A, const std::vector<Real>& B)
{
	const size_t L = std::max(A.size(), B.size());
	std::vector<Real> Out(L, Real(0));
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

// Horner's rule at s = jw.
FieldsComplex PolyAtJw(const std::vector<Real>& P, const Real& W)
{
	FieldsComplex Acc{Real(0), Real(0)};
	for (const Real& X : P)
	{
		Acc = {X - Acc.Im * W, Acc.Re * W};
	}
	return Acc;
}

bool ReadTf(const Matrix& M, std::vector<Real>* Num, std::vector<Real>* Den, std::string* Error)
{
	if (M.IsScalar())
	{
		*Num = {M.V[0]};
		*Den = {Real(1)};
		return true;
	}
	if (M.Rows != 2 || M.Cols < 1)
	{
		*Error = "a transfer function is tf(num, den): a 2-row matrix [num; den], highest power first; got " + ShapeText(M);
		return false;
	}
	*Num = TrimPoly(std::vector<Real>(M.V.begin(), M.V.begin() + M.Cols));
	*Den = TrimPoly(std::vector<Real>(M.V.begin() + M.Cols, M.V.end()));
	if (IsZeroPoly(*Den))
	{
		*Error = "denominator is zero";
		return false;
	}
	return true;
}

Matrix PackTf(const std::vector<Real>& NumIn, const std::vector<Real>& DenIn)
{
	const std::vector<Real> Num = TrimPoly(NumIn), Den = TrimPoly(DenIn);
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
	const std::vector<Real> Den = TrimPoly(D.V);
	if (IsZeroPoly(Den))
	{
		return Result::Fail("denominator is zero");
	}
	return Result::Of(PackTf(N.V, Den));
}

Result TfSeries(const Call& C)
{
	std::vector<Real> N1, D1, N2, D2;
	std::string E;
	if (!ReadTf(C.M(0), &N1, &D1, &E) || !ReadTf(C.M(1), &N2, &D2, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(PackTf(PolyTimes(N1, N2), PolyTimes(D1, D2)));
}

Result TfParallel(const Call& C)
{
	std::vector<Real> N1, D1, N2, D2;
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
	std::vector<Real> Ng, Dg, Nh, Dh;
	std::string E;
	if (!ReadTf(C.M(0), &Ng, &Dg, &E) || !ReadTf(C.M(1), &Nh, &Dh, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<Real> Den = TrimPoly(PolyPlus(PolyTimes(Dg, Dh), PolyTimes(Ng, Nh)));
	if (IsZeroPoly(Den))
	{
		return Result::Fail("closed loop is degenerate (1 + G H = 0)");
	}
	return Result::Of(PackTf(PolyTimes(Ng, Dh), Den));
}

Result FrequencyRecord(const FieldsComplex& G)
{
	const Real Mag = gm::hypot(G.Re, G.Im);
	return Result::Record({{"re", G.Re}, {"im", G.Im}, {"mag", Mag}, {"db", Real(20) * gm::log10(Mag)}, {"phase", gm::atan2(G.Im, G.Re)}});
}

Result TfJw(const Call& C)
{
	std::vector<Real> N, D;
	std::string E;
	if (!ReadTf(C.M(0), &N, &D, &E))
	{
		return Result::Fail(E);
	}
	const FieldsComplex Den = PolyAtJw(D, C.N(1));
	if (Den.Re.IsZero() && Den.Im.IsZero())
	{
		return Result::Fail("s = jw is a pole");
	}
	return FrequencyRecord(ComplexDiv(PolyAtJw(N, C.N(1)), Den));
}

// SISO state-space frequency response C (jw I - A)^-1 B + D.
Result SsJw(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &Cm = C.M(2), &D = C.M(3);
	const Real W = C.N(4);
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
	Real Re = D.V[0], Im = 0;
	for (int i = 0; i < N; ++i)
	{
		Re += Cm(0, i) * X(i, 0);
		Im += Cm(0, i) * X(N + i, 0);
	}
	return FrequencyRecord({Re, Im});
}

// ---------------------------------------------------------------- ch10: sheaves on a cycle

Matrix Rotation2(const Real& T)
{
	const Real Cs = gm::cos(T), Sn = gm::sin(T);
	return Matrix::FromRows({{Cs, -Sn}, {Sn, Cs}});
}

bool SheafLaplacian(const Matrix& Frames, const Real& Err, Matrix* L, std::string* Error)
{
	if (!Frames.IsVector() || Frames.Count() > 30)
	{
		*Error = "frames must be a vector of sensor frame angles (rad), at most 30";
		return false;
	}
	const int N = Frames.Count();
	Matrix Delta(2 * N, 2 * N);
	for (int i = 0; i < N; ++i)
	{
		const int j = (i + 1) % N;
		const Matrix Fi = Rotation2(Frames.V[static_cast<size_t>(i)]);
		Matrix Fj = Rotation2(Frames.V[static_cast<size_t>(j)]);
		if (j == 0)
		{
			Fj = la::Mul(Rotation2(Err), Fj);
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
	std::vector<Real> W;
	std::string E;
	if (!SheafLaplacian(C.M(0), C.N(1), &L, &E) || !la::SymmetricEigen(L, &W, &V, &E))
	{
		return Result::Fail(E);
	}
	const Real Zero(1e-9);
	int H0 = 0;
	Real Lmin = Real::NaN();
	for (const Real& X : W)
	{
		if (X < Zero)
		{
			++H0;
		}
		else if (Lmin.IsNaN() || X < Lmin)
		{
			Lmin = X;
		}
	}
	return Result::Record({{"h0", Real(H0)}, {"lmin", Lmin}});
}

// ---------------------------------------------------------------- ch11: cellular automata

void RunAutomaton(int Rule, int N, int Steps, std::vector<uint8_t>* Row, std::vector<uint8_t>* Centre)
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
}

bool ReadRule(const Call& C, int* Rule, int* Steps, long long MaxSteps, std::string* Error)
{
	if (!IsWhole(C.N(0), 0, 255))
	{
		*Error = "rule must be a whole number 0..255";
		return false;
	}
	if (!IsWhole(C.N(1), 0, MaxSteps))
	{
		*Error = "steps must be a whole number 0.." + std::to_string(MaxSteps);
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
	if (!IsWhole(C.N(2), 1, 20001) || static_cast<long long>(Steps) * C.Int(2) > 20000000LL)
	{
		return Result::Fail("width must be 1..20001 and steps * width at most 2e7");
	}
	const int N = C.Int(2);
	std::vector<uint8_t> Row;
	RunAutomaton(Rule, N, Steps, &Row, nullptr);
	Matrix M(1, N);
	for (int i = 0; i < N; ++i)
	{
		M.V[static_cast<size_t>(i)] = Real(static_cast<int>(Row[static_cast<size_t>(i)]));
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
		M.V[i] = Real(static_cast<int>(Centre[i]));
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
	long long Ones = 0;
	for (size_t i = 1; i < Centre.size(); ++i)
	{
		Ones += Centre[i];
	}
	return Result::Of(Real(Ones) / Real(static_cast<long long>(Centre.size() - 1)));
}

// ---------------------------------------------------------------- ch11: iterated function systems

Result Moran(const Call& C)
{
	const Matrix& R = C.M(0);
	if (!R.IsVector())
	{
		return Result::Fail("ratios must be a vector of contraction factors in (0, 1)");
	}
	std::vector<Real> Logs;
	for (const Real& X : R.V)
	{
		if (!(X > 0 && X < 1))
		{
			return Result::Fail("every ratio must be in (0, 1)");
		}
		Logs.push_back(gm::log(X));
	}
	if (R.Count() == 1)
	{
		return Result::Of(Real(0));
	}
	auto Excess = [&Logs](const Real& D) {
		Real S = -1;
		for (const Real& L : Logs)
		{
			S += gm::exp(D * L);
		}
		return S;
	};
	Real Lo = 0, Hi = 1;
	while (Excess(Hi) > 0 && Hi < Real(1000000))
	{
		Lo = Hi;
		Hi = Hi * 2;
	}
	const Real Tol(1e-17);
	for (int It = 0; It < 200 && Hi - Lo > Tol * Hi; ++It)
	{
		const Real Mid = (Lo + Hi) / 2;
		if (Excess(Mid) > 0)
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

Result IfsPreset(const Call& C)
{
	static const Real Fern[] = {0, 0, 0, 0.16, 0, 0, 0.01, 0.85, 0.04, -0.04, 0.85, 0, 1.6, 0.85, 0.2, -0.26, 0.23, 0.22, 0, 1.6, 0.07, -0.15, 0.28, 0.26, 0.24, 0, 0.44, 0.07};
	static const Real Sierpinski[] = {0.5, 0, 0, 0.5, 0, 0, 0.3333, 0.5, 0, 0, 0.5, 0.5, 0, 0.3333, 0.5, 0, 0, 0.5, 0.25, 0.433, 0.3334};
	static const Real Tree[] = {0, 0, 0, 0.5, 0, 0, 0.05, 0.42, -0.42, 0.42, 0.42, 0, 0.2, 0.4, 0.42, 0.42, -0.42, 0.42, 0, 0.2, 0.4, 0.1, 0, 0, 0.1, 0, 0.2, 0.15};
	static const Real Expand[] = {0.5, 0, 0, 0.5, 0, 0, 0.3333, 0.5, 0, 0, 0.5, 0.5, 0, 0.3333, 1.05, 0, 0, 1.05, 0.25, 0.433, 0.3334};
	static const Real AvgExpand[] = {0.5, 0, 0, 0.5, 0, 0, 0.1, 1.2, 0, 0, 1.2, 0.5, 0, 0.9};
	const Real* Tables[] = {Fern, Sierpinski, Tree, Expand, AvgExpand};
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

bool ReadIfs(const Matrix& M, std::vector<Real>* Norms, std::vector<Real>* Weights, std::string* Error)
{
	if (M.Rows < 1 || (M.Cols != 6 && M.Cols != 7))
	{
		*Error = "maps are rows [a b c d e f] or [a b c d e f p]: x' = a x + b y + e, y' = c x + d y + f";
		return false;
	}
	Norms->clear();
	Weights->clear();
	Real Sum = 0;
	for (int i = 0; i < M.Rows; ++i)
	{
		const Real A = M(i, 0), B = M(i, 1), Cc = M(i, 2), D = M(i, 3);
		const Real S1 = A * A + B * B + Cc * Cc + D * D, Det = A * D - B * Cc;
		Norms->push_back(gm::sqrt((S1 + gm::sqrt(gm::max(Real(0), S1 * S1 - Real(4) * Det * Det))) / 2));
		Weights->push_back(M.Cols == 7 ? gm::max(Real(0), M(i, 6)) : Real(1));
		Sum += Weights->back();
	}
	const Real Div = Sum > 0 ? Sum : Real(1);
	for (Real& W : *Weights)
	{
		W /= Div;
	}
	return true;
}

Result IfsNorm(const Call& C)
{
	std::vector<Real> Norms, Weights;
	std::string E;
	if (!ReadIfs(C.M(0), &Norms, &Weights, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Matrix::Column(Norms));
}

Result IfsLog(const Call& C)
{
	std::vector<Real> Norms, Weights;
	std::string E;
	if (!ReadIfs(C.M(0), &Norms, &Weights, &E))
	{
		return Result::Fail(E);
	}
	Real S = 0;
	for (size_t i = 0; i < Norms.size(); ++i)
	{
		S += Weights[i] * gm::log(gm::max(Norms[i], Real(1e-300)));
	}
	return Result::Of(S);
}

// ---------------------------------------------------------------- ch11: Turing, L-systems, splines, Lloyd

Result Turing(const Call& C)
{
	const Real Dv = C.N(0), Du = C.N(1), A = C.N(2), B = C.N(3);
	if (!(Du > 0 && Dv > 0 && A + B > 0))
	{
		return Result::Fail("needs positive diffusivities and a + b > 0");
	}
	const Real Us = A + B, Vs = B / (Us * Us);
	const Real J11 = Real(-1) + Real(2) * Us * Vs, J12 = Us * Us, J21 = Real(-2) * Us * Vs, J22 = -Us * Us, J1221 = J12 * J21;
	if (!((J11 + J22).IsNegative() && J11 * J22 - J1221 > 0))
	{
		return Result::Fail("the uniform state is unstable without diffusion (Turing needs tr J < 0 and det J > 0)");
	}
	auto Growth = [&](const Real& Q) {
		const Real P = J11 - Du * Q, S = J22 - Dv * Q, Tr = P + S, Det = P * S - J1221, Disc = Tr * Tr - Real(4) * Det;
		return Disc >= 0 ? (Tr + gm::sqrt(Disc)) / 2 : Tr / 2;
	};
	// Growth > 0 exactly where det(J - D k^2) < 0: a quadratic in q = k^2.
	const Real A2 = Du * Dv, A1 = -(J11 * Dv + J22 * Du), A0 = J11 * J22 - J1221, Disc = A1 * A1 - Real(4) * A2 * A0;
	if (!(A1.IsNegative() && Disc > 0))
	{
		const Real Nan = Real::NaN();
		return Result::Record({{"kmin", Nan}, {"kmax", Nan}, {"kstar", Nan}, {"wl", Nan}, {"gmax", Nan}});
	}
	const Real Root = gm::sqrt(Disc);
	const Real Q1 = (-A1 - Root) / (Real(2) * A2), Q2 = (-A1 + Root) / (Real(2) * A2);
	Real Lo = Q1, Hi = Q2;
	const Real Golden = (gm::sqrt(Real(5)) - 1) / 2, Tol(1e-16);
	for (int It = 0; It < 300 && Hi - Lo > Tol * Hi; ++It)
	{
		const Real X1 = Hi - Golden * (Hi - Lo), X2 = Lo + Golden * (Hi - Lo);
		if (Growth(X1) < Growth(X2))
		{
			Lo = X1;
		}
		else
		{
			Hi = X2;
		}
	}
	const Real Qs = (Lo + Hi) / 2, Ks = gm::sqrt(Qs);
	return Result::Record({{"kmin", gm::sqrt(Q1)}, {"kmax", gm::sqrt(Q2)}, {"kstar", Ks}, {"wl", cst::Tau() / Ks}, {"gmax", Growth(Qs)}});
}

Result LSys(const Call& C)
{
	if (!IsWhole(C.N(0), 0, 74))
	{
		return Result::Fail("n must be a whole number 0..74 (word lengths stay exact)");
	}
	long long A = 1, B = 0;
	for (int k = 0; k < C.Int(0); ++k)
	{
		const long long NewA = A + B;   // a -> ab, b -> a
		B = A;
		A = NewA;
	}
	return Result::Record({{"len", Real(A + B)}, {"a", Real(A)}, {"b", Real(B)}, {"ratio", B > 0 ? Real(A) / Real(B) : Real::Infinity()}});
}

bool NaturalSpline(const Matrix& Xs, const Matrix& Ys, std::vector<Real>* M, std::string* Error)
{
	if (!Xs.IsVector() || !Ys.IsVector() || Xs.Count() != Ys.Count() || Xs.Count() < 2)
	{
		*Error = "knots and values must be vectors of equal length >= 2";
		return false;
	}
	const int N = Xs.Count() - 1;
	const std::vector<Real>&X = Xs.V, &Y = Ys.V;
	for (int i = 0; i < N; ++i)
	{
		if (!(X[static_cast<size_t>(i + 1)] > X[static_cast<size_t>(i)]))
		{
			*Error = "knots must be strictly increasing";
			return false;
		}
	}
	// Tridiagonal system for the interior second derivatives (Thomas algorithm); M_0 = M_n = 0.
	M->assign(static_cast<size_t>(N + 1), Real(0));
	if (N < 2)
	{
		return true;
	}
	std::vector<Real> Diag(static_cast<size_t>(N + 1)), Up(static_cast<size_t>(N + 1)), Rhs(static_cast<size_t>(N + 1));
	for (int i = 1; i < N; ++i)
	{
		const size_t u = static_cast<size_t>(i);
		const Real H0 = X[u] - X[u - 1], H1 = X[u + 1] - X[u];
		Diag[u] = Real(2) * (H0 + H1);
		Up[u] = H1;
		Rhs[u] = Real(6) * ((Y[u + 1] - Y[u]) / H1 - (Y[u] - Y[u - 1]) / H0);
		if (i > 1)
		{
			const Real W = H0 / Diag[u - 1];
			Diag[u] -= W * Up[u - 1];
			Rhs[u] -= W * Rhs[u - 1];
		}
	}
	for (int i = N - 1; i >= 1; --i)
	{
		const size_t u = static_cast<size_t>(i);
		(*M)[u] = (Rhs[u] - (i < N - 1 ? Up[u] * (*M)[u + 1] : Real(0))) / Diag[u];
	}
	return true;
}

Result Spline(const Call& C)
{
	std::vector<Real> M;
	std::string E;
	if (!NaturalSpline(C.M(0), C.M(1), &M, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<Real>&X = C.M(0).V, &Y = C.M(1).V;
	const int N = static_cast<int>(X.size()) - 1;
	Matrix Out = C.M(2);
	for (Real& T : Out.V)
	{
		int i = static_cast<int>(std::upper_bound(X.begin(), X.end(), T) - X.begin()) - 1;
		i = std::max(0, std::min(N - 1, i));
		const size_t u = static_cast<size_t>(i);
		const Real H = X[u + 1] - X[u], A = X[u + 1] - T, B = T - X[u], H6 = Real(6) * H;
		T = M[u] * A * A * A / H6 + M[u + 1] * B * B * B / H6 + (Y[u] / H - M[u] * H / 6) * A + (Y[u + 1] / H - M[u + 1] * H / 6) * B;
	}
	return Result::Of(Out);
}

Result SplineEnergy(const Call& C)
{
	std::vector<Real> M;
	std::string E;
	if (!NaturalSpline(C.M(0), C.M(1), &M, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<Real>& X = C.M(0).V;
	Real S = 0;
	for (size_t i = 0; i + 1 < X.size(); ++i)
	{
		S += (X[i + 1] - X[i]) * (M[i] * M[i] + M[i] * M[i + 1] + M[i + 1] * M[i + 1]) / 3;   // f'' is linear on each piece
	}
	return Result::Of(S);
}

// The chapter's drone-coverage bench: 12 sites from mulberry32(12) on a 90 x 90 grid. The density is
// separable, exp(-(dx^2 + dy^2)/0.05) = exp(-dx^2/0.05) exp(-dy^2/0.05), and the energy is only read
// at the first and last step, which keeps the 61 sweeps quick on the gears.
Result Lloyd(const Call& C)
{
	const Real Hx = C.N(0);
	if (!IsWhole(C.N(1), 0, 500))
	{
		return Result::Fail("steps must be a whole number 0..500");
	}
	const int Grid = 90, Sites = 12, Steps = C.Int(1), Cells = Grid * Grid;
	std::vector<Real> G(Grid), Ex(Grid), Ey(Grid);
	const Real Width(0.05), Hy(0.3), Base(0.2);
	for (int i = 0; i < Grid; ++i)
	{
		G[static_cast<size_t>(i)] = Real(2 * i + 1) / Real(2 * Grid);
		const Real Dx = G[static_cast<size_t>(i)] - Hx, Dy = G[static_cast<size_t>(i)] - Hy;
		Ex[static_cast<size_t>(i)] = gm::exp(-Dx * Dx / Width);
		Ey[static_cast<size_t>(i)] = gm::exp(-Dy * Dy / Width);
	}
	std::vector<Real> Rho(Cells), RhoX(Cells), RhoY(Cells);
	Real Mass = 0;
	for (int j = 0; j < Grid; ++j)
	{
		for (int i = 0; i < Grid; ++i)
		{
			const size_t p = static_cast<size_t>(j * Grid + i);
			Rho[p] = Ey[static_cast<size_t>(j)] * Ex[static_cast<size_t>(i)] + Base;
			RhoX[p] = Rho[p] * G[static_cast<size_t>(i)];
			RhoY[p] = Rho[p] * G[static_cast<size_t>(j)];
			Mass += Rho[p];
		}
	}
	std::vector<Real> Sx(Sites), Sy(Sites), Cx(Sites), Cy(Sites);
	std::vector<bool> Has(Sites);
	FieldsMulberry Rng(12);
	for (int k = 0; k < Sites; ++k)
	{
		Sx[static_cast<size_t>(k)] = Rng.Next();
		Sy[static_cast<size_t>(k)] = Rng.Next();
	}
	std::vector<Real> Dxk(static_cast<size_t>(Grid * Sites)), Dyk(static_cast<size_t>(Grid * Sites));
	auto Evaluate = [&](bool WithEnergy) {
		for (int i = 0; i < Grid; ++i)
		{
			for (int k = 0; k < Sites; ++k)
			{
				const Real Dx = G[static_cast<size_t>(i)] - Sx[static_cast<size_t>(k)], Dy = G[static_cast<size_t>(i)] - Sy[static_cast<size_t>(k)];
				Dxk[static_cast<size_t>(i * Sites + k)] = Dx * Dx;
				Dyk[static_cast<size_t>(i * Sites + k)] = Dy * Dy;
			}
		}
		std::vector<Real> Wx(Sites, Real(0)), Wy(Sites, Real(0)), Wm(Sites, Real(0));
		Real En = 0;
		int Best = 0;
		for (int j = 0; j < Grid; ++j)
		{
			const Real* Ry = &Dyk[static_cast<size_t>(j * Sites)];
			for (int i = 0; i < Grid; ++i)
			{
				const Real* Rx = &Dxk[static_cast<size_t>(i * Sites)];
				// Start from the previous cell's owner; a site is skipped when one coordinate alone is farther.
				Real Bd = Rx[Best] + Ry[Best];
				for (int k = 0; k < Sites; ++k)
				{
					if (k == Best || Ry[k] > Bd || Rx[k] > Bd)
					{
						continue;
					}
					const Real D = Rx[k] + Ry[k];
					if (D < Bd || (D == Bd && k < Best))
					{
						Bd = D;
						Best = k;
					}
				}
				const size_t p = static_cast<size_t>(j * Grid + i), b = static_cast<size_t>(Best);
				if (WithEnergy)
				{
					En += Rho[p] * Bd;
				}
				Wx[b] += RhoX[p];
				Wy[b] += RhoY[p];
				Wm[b] += Rho[p];
			}
		}
		for (size_t k = 0; k < static_cast<size_t>(Sites); ++k)
		{
			Has[k] = Wm[k] > 0;
			Cx[k] = Has[k] ? Wx[k] / Wm[k] : Sx[k];
			Cy[k] = Has[k] ? Wy[k] / Wm[k] : Sy[k];
		}
		return En / Mass;
	};
	const Real E0 = Evaluate(true);
	Real En = E0;
	for (int s = 0; s < Steps; ++s)
	{
		Sx = Cx;
		Sy = Cy;
		En = Evaluate(s == Steps - 1);
	}
	Real Gap = 0;
	for (size_t k = 0; k < static_cast<size_t>(Sites); ++k)
	{
		if (Has[k])
		{
			Gap = gm::max(Gap, gm::hypot(Sx[k] - Cx[k], Sy[k] - Cy[k]));
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
	const Matrix Info = la::Mul(la::Mul(la::Transpose(H), Ri), H), Ft = la::Transpose(F);
	for (int k = 0; k < Steps; ++k)
	{
		Matrix Jinv, Pred;
		if (!la::Inverse(J, &Jinv, Error))
		{
			return false;
		}
		if (!la::Inverse(la::Add(Q, la::Mul(la::Mul(F, Jinv), Ft)), &Pred, Error))
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
		return Result::Fail("k must be a whole number of steps 0..100000");
	}
	Matrix P;
	std::string E;
	if (!PcrbRun(C.M(0), C.M(1), C.M(2), C.M(3), C.M(4), C.Int(5), &P, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(P);
}

Matrix WhiteAccelNoise(const Real& Q, const Real& Dt)
{
	const Real Dt2 = Dt * Dt, Off = Q * Dt2 / 2;
	return Matrix::FromRows({{Q * Dt2 * Dt / 3, Off}, {Off, Q * Dt}});
}

Result PcrbCv(const Call& C)
{
	const Real Q = C.N(1), R = C.N(2), Dt = C.N(3);
	if (!IsWhole(C.N(0), 0, 100000) || !(Q > 0 && R > 0 && Dt > 0))
	{
		return Result::Fail("needs whole k in 0..100000 and positive q, r, dt");
	}
	const Matrix F = Matrix::FromRows({{Real(1), Dt}, {Real(0), Real(1)}});
	Matrix P;
	std::string E;
	if (!PcrbRun(F, WhiteAccelNoise(Q, Dt), Matrix::FromRows({{Real(1), Real(0)}}), Matrix::Scalar(R), Matrix::FromRows({{Real(4), Real(0)}, {Real(0), Real(1)}}), C.Int(0), &P,
			&E))
	{
		return Result::Fail(E);
	}
	return Result::Record({{"p11", P(0, 0)}, {"p12", P(0, 1)}, {"p22", P(1, 1)}});
}

// ---------------------------------------------------------------- ch12: minimum-time HJB

Result BangBang(const Call& C)
{
	const Real X = C.N(0), V = C.N(1), S = X + V * gm::fabs(V) / 2, Half = V * V / 2;
	const Real Zero(0), One(1);
	if (X.IsZero() && V.IsZero())
	{
		return Result::Record({{"T", Zero}, {"u0", Zero}, {"ts", Zero}, {"xs", Zero}, {"vs", Zero}});
	}
	if (S > 0 || (S.IsZero() && V.IsNegative()))
	{
		const Real Root = gm::sqrt(X + Half);
		if (S.IsZero())
		{
			return Result::Record({{"T", V + Real(2) * Root}, {"u0", One}, {"ts", Zero}, {"xs", X}, {"vs", V}});
		}
		return Result::Record({{"T", V + Real(2) * Root}, {"u0", -One}, {"ts", V + Root}, {"xs", Root * Root / 2}, {"vs", -Root}});
	}
	const Real Root = gm::sqrt(gm::max(Zero, Half - X));
	if (S.IsZero())
	{
		return Result::Record({{"T", Real(2) * Root - V}, {"u0", -One}, {"ts", Zero}, {"xs", X}, {"vs", V}});
	}
	return Result::Record({{"T", Real(2) * Root - V}, {"u0", One}, {"ts", Root - V}, {"xs", -Root * Root / 2}, {"vs", Root}});
}

// The chapter's simulation: switch on the curve x = -v|v|/2, semi-implicit Euler, stop within 1e-3 of rest.
Result BbSim(const Call& C)
{
	Real X = C.N(0), V = C.N(1);
	const Real Dt = C.N(2);
	if (!(Dt >= Real(1e-4) && Dt <= Real(0.1)))
	{
		return Result::Fail("dt must be in [1e-4, 0.1] s");
	}
	const Real Near(1e-6), Limit(50);
	Real T = 0;
	while (X * X + V * V > Near && T < Limit)
	{
		const Real S = X + V * gm::fabs(V) / 2;
		const bool Down = S > 0 || (S.IsZero() && V > 0), Up = S.IsNegative() || (S.IsZero() && V.IsNegative());
		if (Down)
		{
			V -= Dt;
		}
		else if (Up)
		{
			V += Dt;
		}
		X += V * Dt;
		T += Dt;
	}
	return Result::Of(T);
}

// ---------------------------------------------------------------- ch12: Fokker-Planck and Ito

Result OuProc(const Call& C)
{
	const Real X0 = C.N(0), Th = C.N(1), Sg = C.N(2), T = C.N(3), Sg2 = Sg * Sg;
	if (T.IsNegative())
	{
		return Result::Fail("t must be >= 0");
	}
	const Real Var = Th.IsZero() ? Sg2 * T : -Sg2 * gm::expm1(Real(-2) * Th * T) / (Real(2) * Th);
	return Result::Record({{"mean", X0 * gm::exp(-Th * T)}, {"var", Var}, {"sd", gm::sqrt(Var)}, {"vinf", Th > 0 ? Sg2 / (Real(2) * Th) : Real::Infinity()}});
}

Result QCv(const Call& C) { return Result::Of(WhiteAccelNoise(C.N(0), C.N(1))); }

// Van Loan: expm([[-A, Qc], [0, A']] dt) = [[., Phi^-1 Qd], [0, Phi']].
Result QDisc(const Call& C)
{
	const Matrix &A = C.M(0), &Qc = C.M(1);
	const Real Dt = C.N(2);
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
	Matrix Out(A.Rows, B.Cols, Real::Infinity());
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int k = 0; k < A.Cols; ++k)
		{
			const Real Aik = A(i, k);
			if (Aik.IsInf() && !Aik.IsNegative())
			{
				continue;
			}
			for (int j = 0; j < B.Cols; ++j)
			{
				const Real S = Aik + B(k, j);
				if (S < Out(i, j))
				{
					Out(i, j) = S;
				}
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
	if (!W.IsSquare() || W.Rows > 120)
	{
		return Result::Fail("W must be a square edge-weight matrix (inf = no edge), at most 120 x 120");
	}
	const Real K = C.N(1).IsNegative() ? Real(W.Rows - 1) : C.N(1);
	if (!IsWhole(K, 0, 1000000000LL))
	{
		return Result::Fail("k must be a whole number >= 0 (default n - 1: shortest paths)");
	}
	Matrix R(W.Rows, W.Rows, Real::Infinity()), B = W;
	for (int i = 0; i < W.Rows; ++i)
	{
		R(i, i) = 0;
	}
	long long N = K.ToLong();
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

// Sinkhorn with unit marginals, u = 1 / K v, v = 1 / K' u, K = exp(-C / eps). The gears have no underflow
// (the exponent is unbounded for practical purposes), so the plain scaling form is safe for small eps.
// Stops early once every row sum of the plan is within 1e-16 of 1.
bool SinkhornPlan(const Matrix& Cm, const Real& Eps, int Iters, Matrix* Plan, std::string* Error)
{
	const int N = Cm.Rows;
	if (!Cm.IsSquare() || N > 40)
	{
		*Error = "cost must be a square matrix (equal swarms), at most 40 x 40";
		return false;
	}
	if (!(Eps > 0))
	{
		*Error = "epsilon must be positive";
		return false;
	}
	Matrix K(N, N);
	for (size_t i = 0; i < K.V.size(); ++i)
	{
		K.V[i] = gm::exp(-Cm.V[i] / Eps);
	}
	std::vector<Real> U(static_cast<size_t>(N), Real(1)), V(static_cast<size_t>(N), Real(1)), Kv(static_cast<size_t>(N));
	const Real Tol(1e-16), One(1);
	for (int It = 0; It < Iters; ++It)
	{
		Real Worst = 0;
		for (int i = 0; i < N; ++i)
		{
			Real S = 0;
			for (int j = 0; j < N; ++j)
			{
				S += K(i, j) * V[static_cast<size_t>(j)];
			}
			Kv[static_cast<size_t>(i)] = S;
			if (It > 0)
			{
				Worst = gm::max(Worst, gm::fabs(U[static_cast<size_t>(i)] * S - One));
			}
		}
		if (It > 0 && Worst < Tol)
		{
			break;
		}
		for (int i = 0; i < N; ++i)
		{
			U[static_cast<size_t>(i)] = One / Kv[static_cast<size_t>(i)];
		}
		for (int j = 0; j < N; ++j)
		{
			Real S = 0;
			for (int i = 0; i < N; ++i)
			{
				S += K(i, j) * U[static_cast<size_t>(i)];
			}
			V[static_cast<size_t>(j)] = One / S;
		}
	}
	*Plan = Matrix(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			(*Plan)(i, j) = U[static_cast<size_t>(i)] * K(i, j) * V[static_cast<size_t>(j)];
		}
	}
	if (!la::MaxAbs(*Plan).IsFinite())
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
	if (static_cast<long long>(C.Int(2)) * C.M(0).Count() > 2000000LL)
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
	RealAccumulator S;
	for (size_t i = 0; i < Plan.V.size(); ++i)
	{
		S.AddProduct(Plan.V[i], C.M(0).V[i]);
	}
	return Result::Of(S.Result());
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
	std::vector<Real> Best(static_cast<size_t>(1) << N, Real::Infinity());
	Best[0] = 0;
	for (unsigned Mask = 0; Mask < (1u << N); ++Mask)
	{
		if (Best[Mask].IsInf())
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
				const Real S = Best[Mask] + Cm(Row, j);
				if (S < Best[Next])
				{
					Best[Next] = S;
				}
			}
		}
	}
	return Result::Of(Best[(1u << N) - 1]);
}

// ---------------------------------------------------------------- ch12: symplectic integration, random matrices

// |r| from r^2 by Newton's method started at the previous radius: the orbit moves little per step, so two
// or three divisions replace a full Toepler root.
Real RadiusFrom(const Real& R2, const Real& Guess)
{
	Real X = Guess;
	for (int It = 0; It < 12; ++It)
	{
		const Real Y = (X + R2 / X) / 2;
		if (It > 0 && Y >= X)
		{
			return X;
		}
		X = Y;
	}
	return X;
}

Real KeplerWorst(bool Leapfrog, int Orbits, int StepsPer)
{
	Real Rx = 1, Ry = 0, Vx = 0, Vy = 1, Worst = 0;
	const Real H = cst::Tau() / Real(StepsPer), HalfH = H / 2, E0(-0.5), One(1);
	Real Rn = 1, Ax = -1, Ay = 0;
	auto Accel = [&]() {
		const Real R2 = Rx * Rx + Ry * Ry;
		Rn = RadiusFrom(R2, Rn);
		const Real Inv3 = One / (Rn * R2);
		Ax = -Rx * Inv3;
		Ay = -Ry * Inv3;
	};
	const long Total = static_cast<long>(Orbits) * StepsPer;
	for (long k = 0; k < Total; ++k)
	{
		if (Leapfrog)
		{
			Vx += HalfH * Ax;
			Vy += HalfH * Ay;
			Rx += H * Vx;
			Ry += H * Vy;
			Accel();
			Vx += HalfH * Ax;
			Vy += HalfH * Ay;
		}
		else
		{
			Rx += H * Vx;
			Ry += H * Vy;
			Vx += H * Ax;
			Vy += H * Ay;
			Accel();
		}
		const Real En = (Vx * Vx + Vy * Vy) / 2 - One / Rn;
		Worst = gm::max(Worst, gm::fabs(En - E0));
		if (Worst > One)
		{
			break;
		}
	}
	return Worst;
}

Result KepEnergy(const Call& C)
{
	if (!IsWhole(C.N(0), 1, 10000) || !IsWhole(C.N(1), 4, 100000) || static_cast<long long>(C.Int(0)) * C.Int(1) > 200000LL)
	{
		return Result::Fail("orbits and steps per orbit must be whole numbers, orbits * steps <= 2e5");
	}
	return Result::Record({{"lf", KeplerWorst(true, C.Int(0), C.Int(1))}, {"eu", KeplerWorst(false, C.Int(0), C.Int(1))}});
}

bool MpSupport(const Real& Cr, const Real& S2, Real* Lo, Real* Hi)
{
	if (!(Cr > 0 && S2 > 0))
	{
		return false;
	}
	const Real Root = gm::sqrt(Cr), A = Real(1) - Root, B = Real(1) + Root;
	*Lo = S2 * A * A;
	*Hi = S2 * B * B;
	return true;
}

Result MpEdges(const Call& C)
{
	Real Lo, Hi;
	if (!MpSupport(C.N(0), C.N(1), &Lo, &Hi))
	{
		return Result::Fail("aspect ratio c = p/n and variance must be positive");
	}
	return Result::Record({{"lo", Lo}, {"hi", Hi}});
}

Result MpPdf(const Call& C)
{
	const Real Cr = C.N(1), S2 = C.N(2);
	Real Lo, Hi;
	if (!MpSupport(Cr, S2, &Lo, &Hi))
	{
		return Result::Fail("aspect ratio c = p/n and variance must be positive");
	}
	const Real Scale = Real(2) * cst::Pi() * S2 * Cr;
	Matrix X = C.M(0);
	for (Real& V : X.V)
	{
		V = (V > Lo && V < Hi) ? gm::sqrt((Hi - V) * (V - Lo)) / (Scale * V) : Real(0);
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
		"record d (m), t (s)", {Num("el", "rad"), Num("re", "m", 6371000), Num("rs", "m", 26560000)}, &SatRange);
	Def(R, "gnssfix", Lat, "ch07 Lattices: bench 1, pseudorange positioning", "Gauss-Newton position and clock from pseudoranges (start at Earth's centre)",
		"record x, y, z (m), cb = c*bias (m), dt (s), it", {Mat("S", "m (n x 3)"), Mat("rho", "m")}, &GnssFix);
	Def(R, "dop", Lat, "ch07 Lattices: bench 1, dilution of precision", "DOPs from satellite positions seen from receiver p (vertical = p/|p|)",
		"record gdop, pdop, hdop, vdop, tdop", {Mat("S", "m (n x 3)"), Mat("p", "m (3)")}, &Dop);
	Def(R, "dopae", Lat, "ch07 Lattices: bench 1, dilution of precision", "DOPs from satellite azimuths and elevations", "record gdop, pdop, hdop, vdop, tdop",
		{Mat("az", "rad"), Mat("el", "rad")}, &DopAe);
	Def(R, "baarda", Lat, "ch07 Lattices: Seven ways ranging goes wrong (normalized residuals)", "Baarda w-statistics of the pseudorange fix, r_i / (sigma sqrt(1 - h_ii))",
		"column of w, one per satellite", {Mat("S", "m (n x 3)"), Mat("rho", "m"), Num("sigma", "m", 1)}, &Baarda);
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
		"2 x 2 covariance, cycles^2", {Num("a", "cycles", 1.2), Num("b", "cycles", 0.07), Num("theta", "rad", Real(52) * cst::Pi() / Real(180))}, &EllCov);
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
		{Mat("G", "[num; den]"), MatDefault("H", "[num; den] or gain", Real(1))}, &TfFeedback);
	Def(R, "tfjw", Cat, "ch10 Category: Wiring is a functor", "frequency response G(jw)", "record re, im, mag, db, phase (rad)", {Mat("G", "[num; den]"), Num("w", "rad/s")},
		&TfJw);
	Def(R, "ssjw", Cat, "ch10 Category: Wiring is a functor (realizations)", "SISO state-space response C (jw I - A)^-1 B + D", "record re, im, mag, db, phase (rad)",
		{Mat("A", "1/s"), Mat("B", "-"), Mat("C", "-"), Mat("D", "-"), Num("w", "rad/s")}, &SsJw);
	Def(R, "sheaflap", Cat, "ch10 Category: Sensors that disagree about frames", "sheaf Laplacian of a sensor cycle (at most 30) with rotation restriction maps",
		"2n x 2n matrix", {Mat("frames", "rad"), Num("err", "rad", 0)}, &SheafLap);
	Def(R, "sheafh0", Cat, "ch10 Category: Sensors that disagree about frames", "dimension of global sections and slowest diffusion rate (at most 30 sensors)",
		"record h0, lmin", {Mat("frames", "rad"), Num("err", "rad", 0)}, &SheafH0);

	// Chapter 11: procedural systems.
	Def(R, "ca", Pro, "ch11 Procedural: 256 rules, one row at a time", "elementary automaton from one centre cell, periodic, after n steps", "1 x width row of bits",
		{Num("rule", "0..255"), Num("steps", "-"), Num("width", "cells", 201)}, &Ca);
	Def(R, "cacenter", Pro, "ch11 Procedural: 256 rules, one row at a time", "centre-column bits t = 0..steps (infinite lattice, steps <= 3000)", "column of bits",
		{Num("rule", "0..255"), Num("steps", "-")}, &CaCenter);
	Def(R, "cafreq", Pro, "ch11 Procedural: 256 rules, one row at a time", "1-frequency of the centre column over t = 1..steps (steps <= 3000)", "fraction",
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
	Def(R, "bbsim", Exp, "ch12 Expansion: The fastest way to stop", "simulated bang-bang arrival time (bench integrator, dt >= 1e-4)", "time, s",
		{Num("x", "m"), Num("v", "m/s"), Num("dt", "s", 1e-4)}, &BbSim);
	Def(R, "ouproc", Exp, "ch12 Expansion: How uncertainty spreads between measurements", "Ornstein-Uhlenbeck law dX = -theta X dt + sigma dW at time t",
		"record mean, var, sd, vinf", {Num("x0", "-", 2), Num("theta", "1/s", 0.8), Num("sigma", "1/sqrt(s)", 0.6), Num("t", "s", 2)}, &OuProc);
	Def(R, "qcv", Exp, "ch12 Expansion: Ito integration", "white-acceleration process noise q [dt^3/3, dt^2/2; dt^2/2, dt]", "2 x 2", {Num("q", "m^2/s^3"), Num("dt", "s", 1)},
		&QCv);
	Def(R, "qdisc", Exp, "ch12 Expansion: Ito integration (Van Loan)", "discrete process noise integral of e^(A s) Qc e^(A' s) ds over dt", "n x n",
		{Mat("A", "1/s"), Mat("Qc", "-"), Num("dt", "s")}, &QDisc);
	Def(R, "tropmul", Exp, "ch12 Expansion: tropical algebra", "(min, +) matrix product", "matrix", {Mat("A", "-"), Mat("B", "-")}, &TropMul);
	Def(R, "troppow", Exp, "ch12 Expansion: tropical algebra", "(min, +) power W^k (n <= 120); default n - 1 gives all-pairs shortest paths", "matrix",
		{Mat("W", "inf = no edge"), Num("k", "-", -1)}, &TropPow);
	Def(R, "sinkhorn", Exp, "ch12 Expansion: optimal transport", "entropic transport plan between equal swarms (unit marginals)", "n x n plan",
		{Mat("C", "cost"), Num("eps", "cost"), Num("iters", "-", 20000)}, &Sinkhorn);
	Def(R, "sinkcost", Exp, "ch12 Expansion: optimal transport", "cost sum(P .* C) of the Sinkhorn plan", "cost", {Mat("C", "cost"), Num("eps", "cost"), Num("iters", "-", 20000)},
		&SinkCost);
	Def(R, "otexact", Exp, "ch12 Expansion: optimal transport", "optimal assignment cost (the eps -> 0 limit), n <= 16", "cost", {Mat("C", "cost")}, &OtExact);
	Def(R, "kepenergy", Exp, "ch12 Expansion: symplectic integration", "worst Kepler energy error: leapfrog vs explicit Euler (stops past 1; orbits * steps <= 2e5)",
		"record lf, eu", {Num("orbits", "-", 200), Num("steps", "per orbit", 200)}, &KepEnergy);
	Def(R, "mpedges", Exp, "ch12 Expansion: random matrices", "Marchenko-Pastur support s2 (1 -/+ sqrt c)^2", "record lo, hi", {Num("c", "p/n"), Num("s2", "-", 1)}, &MpEdges);
	Def(R, "mppdf", Exp, "ch12 Expansion: random matrices", "Marchenko-Pastur density (continuous part)", "density, same shape as x",
		{Mat("x", "-"), Num("c", "p/n"), Num("s2", "-", 1)}, &MpPdf);
}
} // namespace mcalc
