// Atlas chapters 1, 2, 4, 5 and 6: rotation groups and quaternions (ch01), the topology of attitude
// and navigation functions (ch02), projection and filtering in Hilbert space (ch04), measure, Monte
// Carlo and information (ch05), and linear programming, duality and LQ costates (ch06).
// Conventions: quaternions are 1 x 4 [w, x, y, z], Hamilton product, scalar first; 3-vectors come back
// as 1 x 3 rows; solutions of linear-algebra problems (projections, LP and LASSO optima) as columns.
#include "CalcLibrary.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace mcalc
{
namespace atlas_detail
{
using V3 = std::array<Real, 3>;
using Q4 = std::array<Real, 4>;

// A matrix parameter whose default is a scalar (used for "equal weights" style defaults).
Param MatOr(const char* Name, const char* Unit, Real Default)
{
	Param P = Mat(Name, Unit);
	P.HasDefault = true;
	P.Default = Default;
	return P;
}

bool ToV3(const Matrix& M, V3* Out)
{
	if (!M.IsVector() || M.Count() != 3)
	{
		return false;
	}
	*Out = {M.V[0], M.V[1], M.V[2]};
	return true;
}

bool ToQ4(const Matrix& M, Q4* Out)
{
	if (!M.IsVector() || M.Count() != 4)
	{
		return false;
	}
	*Out = {M.V[0], M.V[1], M.V[2], M.V[3]};
	return true;
}

Matrix RowV3(const V3& A) { return Matrix::Row({A[0], A[1], A[2]}); }
Matrix RowQ4(const Q4& Q) { return Matrix::Row({Q[0], Q[1], Q[2], Q[3]}); }
Real Dot3(const V3& A, const V3& B) { return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]; }
Real Norm3(const V3& A) { return gm::sqrt(Dot3(A, A)); }
V3 Cross3(const V3& A, const V3& B) { return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]}; }
V3 Scale3(const V3& A, Real S) { return {A[0] * S, A[1] * S, A[2] * S}; }
Real NormQ(const Q4& Q) { return gm::sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]); }

Matrix Hat3(const V3& W)
{
	return Matrix::FromRows({{0.0, -W[2], W[1]}, {W[2], 0.0, -W[0]}, {-W[1], W[0], 0.0}});
}

// The vee of the skew part: exact inverse of Hat3 on skew matrices.
V3 Vee3(const Matrix& S)
{
	return {0.5 * (S(2, 1) - S(1, 2)), 0.5 * (S(0, 2) - S(2, 0)), 0.5 * (S(1, 0) - S(0, 1))};
}

V3 MulMV(const Matrix& A, const V3& V)
{
	return {A(0, 0) * V[0] + A(0, 1) * V[1] + A(0, 2) * V[2], A(1, 0) * V[0] + A(1, 1) * V[1] + A(1, 2) * V[2], A(2, 0) * V[0] + A(2, 1) * V[1] + A(2, 2) * V[2]};
}

Real Det3(const Matrix& M)
{
	return M(0, 0) * (M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)) - M(0, 1) * (M(1, 0) * M(2, 2) - M(1, 2) * M(2, 0)) + M(0, 2) * (M(1, 0) * M(2, 1) - M(1, 1) * M(2, 0));
}

// exp(theta [n]x) for a unit axis n (Rodrigues).
Matrix Rodr(const V3& N, Real Theta)
{
	const Matrix K = Hat3(N);
	const Matrix K2 = la::Mul(K, K);
	const Real S = gm::sin(Theta), C = 1.0 - gm::cos(Theta);
	Matrix R = Matrix::Identity(3);
	for (size_t i = 0; i < 9; ++i)
	{
		R.V[i] += S * K.V[i] + C * K2.V[i];
	}
	return R;
}

// R(q) = (w^2 - v.v) I + 2 v v' + 2 w [v]x, for a unit quaternion.
Matrix QuatDcm(const Q4& Q)
{
	const Real W = Q[0];
	const V3 Vv = {Q[1], Q[2], Q[3]};
	const Real Vn = Dot3(Vv, Vv);
	Matrix R(3, 3);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			R(i, j) = (i == j ? W * W - Vn : 0.0) + 2.0 * Vv[static_cast<size_t>(i)] * Vv[static_cast<size_t>(j)];
		}
	}
	R(0, 1) -= 2.0 * W * Vv[2];
	R(0, 2) += 2.0 * W * Vv[1];
	R(1, 0) += 2.0 * W * Vv[2];
	R(1, 2) -= 2.0 * W * Vv[0];
	R(2, 0) -= 2.0 * W * Vv[1];
	R(2, 1) += 2.0 * W * Vv[0];
	return R;
}

Q4 QuatProduct(const Q4& A, const Q4& B)
{
	return {A[0] * B[0] - A[1] * B[1] - A[2] * B[2] - A[3] * B[3], A[0] * B[1] + A[1] * B[0] + A[2] * B[3] - A[3] * B[2], A[0] * B[2] - A[1] * B[3] + A[2] * B[0] + A[3] * B[1],
		A[0] * B[3] + A[1] * B[2] - A[2] * B[1] + A[3] * B[0]};
}

// Fixes the Real-cover sign: w >= 0, and for w = 0 the largest component positive.
Q4 Canonical(Q4 Q)
{
	size_t Big = 0;
	for (size_t i = 1; i < 4; ++i)
	{
		if (gm::fabs(Q[i]) > gm::fabs(Q[Big]))
		{
			Big = i;
		}
	}
	const bool Flip = Q[0] < 0.0 || (Q[0] == 0.0 && Q[Big] < 0.0);
	if (Flip)
	{
		for (Real& X : Q)
		{
			X = -X;
		}
	}
	return Q;
}

bool IsRotation(const Matrix& R, std::string* Why)
{
	if (R.Rows != 3 || R.Cols != 3)
	{
		*Why = "needs a 3 x 3 rotation matrix";
		return false;
	}
	const Matrix E = la::Sub(la::Mul(la::Transpose(R), R), Matrix::Identity(3));
	if (la::MaxAbs(E) > 1e-6 || Det3(R) < 0.0)
	{
		*Why = "not a rotation matrix (needs R'R = I and det R = +1)";
		return false;
	}
	return true;
}

// Shepperd's method, returning w >= 0.
Q4 DcmQuat(const Matrix& R)
{
	const Real T = R(0, 0) + R(1, 1) + R(2, 2);
	Q4 Q = {0, 0, 0, 0};
	if (T >= R(0, 0) && T >= R(1, 1) && T >= R(2, 2))
	{
		const Real W = 0.5 * gm::sqrt(1.0 + T);
		Q = {W, (R(2, 1) - R(1, 2)) / (4 * W), (R(0, 2) - R(2, 0)) / (4 * W), (R(1, 0) - R(0, 1)) / (4 * W)};
	}
	else if (R(0, 0) >= R(1, 1) && R(0, 0) >= R(2, 2))
	{
		const Real X = 0.5 * gm::sqrt(gm::max(0.0, 1.0 + R(0, 0) - R(1, 1) - R(2, 2)));
		Q = {(R(2, 1) - R(1, 2)) / (4 * X), X, (R(0, 1) + R(1, 0)) / (4 * X), (R(0, 2) + R(2, 0)) / (4 * X)};
	}
	else if (R(1, 1) >= R(2, 2))
	{
		const Real Y = 0.5 * gm::sqrt(gm::max(0.0, 1.0 - R(0, 0) + R(1, 1) - R(2, 2)));
		Q = {(R(0, 2) - R(2, 0)) / (4 * Y), (R(0, 1) + R(1, 0)) / (4 * Y), Y, (R(1, 2) + R(2, 1)) / (4 * Y)};
	}
	else
	{
		const Real Z = 0.5 * gm::sqrt(gm::max(0.0, 1.0 - R(0, 0) - R(1, 1) + R(2, 2)));
		Q = {(R(1, 0) - R(0, 1)) / (4 * Z), (R(0, 2) + R(2, 0)) / (4 * Z), (R(1, 2) + R(2, 1)) / (4 * Z), Z};
	}
	const Real N = NormQ(Q);
	for (Real& X : Q)
	{
		X /= N;
	}
	return Canonical(Q);
}

// Matrix logarithm of SO(3) as a rotation vector theta n with theta in [0, pi].
V3 LogSo3(const Matrix& R)
{
	const V3 S2 = {R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1)};   // 2 sin(theta) n
	const Real SinT = 0.5 * Norm3(S2), CosT = 0.5 * (R(0, 0) + R(1, 1) + R(2, 2) - 1.0);
	const Real Theta = gm::atan2(SinT, CosT);
	if (SinT < 1e-15 && CosT > 0.0)
	{
		return Scale3(S2, 0.5);
	}
	if (CosT > -0.5)
	{
		return Scale3(S2, Theta / (2.0 * SinT));
	}
	// Near pi the sine is small: read n from the symmetric part, n n' = (sym(R) - cos I) / (1 - cos).
	Matrix Sym(3, 3);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			Sym(i, j) = 0.5 * (R(i, j) + R(j, i)) - (i == j ? CosT : 0.0);
		}
	}
	int Best = 0;
	for (int i = 1; i < 3; ++i)
	{
		if (Sym(i, i) > Sym(Best, Best))
		{
			Best = i;
		}
	}
	V3 N = {Sym(0, Best), Sym(1, Best), Sym(2, Best)};
	N = Scale3(N, 1.0 / Norm3(N));
	if (Dot3(N, S2) < 0.0)
	{
		N = Scale3(N, -1.0);
	}
	return Scale3(N, Theta);
}

// Vector observations: r and b as n x 3 (one vector per row), weights a scalar (all equal) or n values.
bool Observations(const Call& C, size_t IR, size_t IB, size_t IW, std::vector<V3>* Rs, std::vector<V3>* Bs, std::vector<Real>* Ws, std::string* Err)
{
	const Matrix &MR = C.M(IR), &MB = C.M(IB), &MW = C.M(IW);
	if (MR.Cols != 3 || MB.Cols != 3 || MR.Rows != MB.Rows || MR.Rows < 1)
	{
		*Err = "r and b must be n x 3 matrices of the same size (one vector per row)";
		return false;
	}
	const int N = MR.Rows;
	if (!(MW.IsScalar() || (MW.IsVector() && MW.Count() == N)))
	{
		*Err = "w must be a scalar or have one weight per row";
		return false;
	}
	for (int i = 0; i < N; ++i)
	{
		Rs->push_back({MR(i, 0), MR(i, 1), MR(i, 2)});
		Bs->push_back({MB(i, 0), MB(i, 1), MB(i, 2)});
		const Real W = MW.IsScalar() ? MW.ScalarValue() : MW.V[static_cast<size_t>(i)];
		if (!(W >= 0.0))
		{
			*Err = "weights must be non-negative";
			return false;
		}
		Ws->push_back(W);
	}
	return true;
}

Matrix DavenportK(const std::vector<V3>& Rs, const std::vector<V3>& Bs, const std::vector<Real>& Ws, Matrix* BOut)
{
	Matrix B(3, 3);
	V3 Z = {0, 0, 0};
	for (size_t k = 0; k < Rs.size(); ++k)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				B(i, j) += Ws[k] * Bs[k][static_cast<size_t>(i)] * Rs[k][static_cast<size_t>(j)];
			}
		}
		const V3 Cr = Cross3(Rs[k], Bs[k]);
		for (size_t i = 0; i < 3; ++i)
		{
			Z[i] += Ws[k] * Cr[i];
		}
	}
	const Real Sigma = la::Trace(B);
	Matrix K(4, 4);
	K(0, 0) = Sigma;
	for (int i = 0; i < 3; ++i)
	{
		K(0, i + 1) = Z[static_cast<size_t>(i)];
		K(i + 1, 0) = Z[static_cast<size_t>(i)];
		for (int j = 0; j < 3; ++j)
		{
			K(i + 1, j + 1) = B(i, j) + B(j, i) - (i == j ? Sigma : 0.0);
		}
	}
	if (BOut)
	{
		*BOut = B;
	}
	return K;
}

bool TopEigvec4(const Matrix& K, Q4* Q, std::string* Err)
{
	std::vector<Real> W;
	Matrix V;
	if (!la::SymmetricEigen(K, &W, &V, Err))
	{
		return false;
	}
	*Q = {V(0, 3), V(1, 3), V(2, 3), V(3, 3)};
	const Real N = NormQ(*Q);
	for (Real& X : *Q)
	{
		X /= N;
	}
	*Q = Canonical(*Q);
	return true;
}

// ---------------------------------------------------------------- ch01 groups

Result Hat(const Call& C)
{
	V3 W;
	if (!ToV3(C.M(0), &W))
	{
		return Result::Fail("needs a 3-vector");
	}
	return Result::Of(Hat3(W));
}

Result Vee(const Call& C)
{
	if (C.M(0).Rows != 3 || C.M(0).Cols != 3)
	{
		return Result::Fail("needs a 3 x 3 matrix");
	}
	return Result::Of(RowV3(Vee3(C.M(0))));
}

Result Rodrigues(const Call& C)
{
	V3 N;
	if (!ToV3(C.M(0), &N))
	{
		return Result::Fail("axis must be a 3-vector");
	}
	const Real L = Norm3(N);
	if (!(L > 0.0))
	{
		return Result::Fail("axis must be nonzero");
	}
	return Result::Of(Rodr(Scale3(N, 1.0 / L), C.N(1)));
}

Result ExpSo3(const Call& C)
{
	V3 W;
	if (!ToV3(C.M(0), &W))
	{
		return Result::Fail("needs a rotation vector (3-vector, rad)");
	}
	const Real T = Norm3(W);
	if (T < 1e-12)
	{
		const Matrix K = Hat3(W);
		return Result::Of(la::Add(la::Add(Matrix::Identity(3), K), la::Scale(la::Mul(K, K), 0.5)));
	}
	return Result::Of(Rodr(Scale3(W, 1.0 / T), T));
}

Result LogSo3Fn(const Call& C)
{
	std::string Why;
	if (!IsRotation(C.M(0), &Why))
	{
		return Result::Fail(Why);
	}
	return Result::Of(RowV3(LogSo3(C.M(0))));
}

Result RotAngle(const Call& C)
{
	std::string Why;
	if (!IsRotation(C.M(0), &Why))
	{
		return Result::Fail(Why);
	}
	return Result::Of(Norm3(LogSo3(C.M(0))));
}

Result QuatFn(const Call& C)
{
	V3 N;
	if (!ToV3(C.M(0), &N))
	{
		return Result::Fail("axis must be a 3-vector");
	}
	const Real L = Norm3(N);
	if (!(L > 0.0))
	{
		return Result::Fail("axis must be nonzero");
	}
	const Real H = 0.5 * C.N(1), S = gm::sin(H) / L;
	return Result::Of(RowQ4({gm::cos(H), N[0] * S, N[1] * S, N[2] * S}));
}

Result QuatMul(const Call& C)
{
	Q4 A, B;
	if (!ToQ4(C.M(0), &A) || !ToQ4(C.M(1), &B))
	{
		return Result::Fail("needs two quaternions [w, x, y, z]");
	}
	return Result::Of(RowQ4(QuatProduct(A, B)));
}

Result QuatConj(const Call& C)
{
	Q4 A;
	if (!ToQ4(C.M(0), &A))
	{
		return Result::Fail("needs a quaternion [w, x, y, z]");
	}
	return Result::Of(RowQ4({A[0], -A[1], -A[2], -A[3]}));
}

bool UnitQuat(const Matrix& M, Q4* Q, std::string* Err)
{
	if (!ToQ4(M, Q))
	{
		*Err = "needs a quaternion [w, x, y, z]";
		return false;
	}
	const Real N = NormQ(*Q);
	if (!(N > 0.0))
	{
		*Err = "the zero quaternion is not a rotation";
		return false;
	}
	for (Real& X : *Q)
	{
		X /= N;
	}
	return true;
}

Result QuatRot(const Call& C)
{
	Q4 Q;
	V3 V;
	std::string E;
	if (!UnitQuat(C.M(0), &Q, &E))
	{
		return Result::Fail(E);
	}
	if (!ToV3(C.M(1), &V))
	{
		return Result::Fail("v must be a 3-vector");
	}
	// q (0, v) q*
	const Q4 P = QuatProduct(QuatProduct(Q, {0.0, V[0], V[1], V[2]}), {Q[0], -Q[1], -Q[2], -Q[3]});
	return Result::Of(RowV3({P[1], P[2], P[3]}));
}

Result Quat2Dcm(const Call& C)
{
	Q4 Q;
	std::string E;
	if (!UnitQuat(C.M(0), &Q, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(QuatDcm(Q));
}

Result Dcm2Quat(const Call& C)
{
	std::string Why;
	if (!IsRotation(C.M(0), &Why))
	{
		return Result::Fail(Why);
	}
	return Result::Of(RowQ4(DcmQuat(C.M(0))));
}

Result Slerp(const Call& C)
{
	Q4 A, B;
	std::string E;
	if (!UnitQuat(C.M(0), &A, &E) || !UnitQuat(C.M(1), &B, &E))
	{
		return Result::Fail(E);
	}
	const Real T = C.N(2);
	Real D = A[0] * B[0] + A[1] * B[1] + A[2] * B[2] + A[3] * B[3];
	if (D < 0.0)
	{
		// q and -q are the same attitude: take the short way round.
		for (Real& X : B)
		{
			X = -X;
		}
		D = -D;
	}
	const Real Om = gm::acos(gm::min(1.0, D));
	Real Ka = 1.0 - T, Kb = T;
	if (Om > 1e-9)
	{
		Ka = gm::sin((1.0 - T) * Om) / gm::sin(Om);
		Kb = gm::sin(T * Om) / gm::sin(Om);
	}
	Q4 Q;
	for (size_t i = 0; i < 4; ++i)
	{
		Q[i] = Ka * A[i] + Kb * B[i];
	}
	const Real N = NormQ(Q);
	for (Real& X : Q)
	{
		X /= N;
	}
	return Result::Of(RowQ4(Q));
}

Result Bracket(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1);
	V3 Va, Vb;
	if (ToV3(A, &Va) && ToV3(B, &Vb))
	{
		const Matrix Ha = Hat3(Va), Hb = Hat3(Vb);
		return Result::Of(RowV3(Vee3(la::Sub(la::Mul(Ha, Hb), la::Mul(Hb, Ha)))));
	}
	if (!A.IsSquare() || A.Rows != B.Rows || A.Cols != B.Cols)
	{
		return Result::Fail("needs two 3-vectors (so(3)) or two square matrices of the same size");
	}
	return Result::Of(la::Sub(la::Mul(A, B), la::Mul(B, A)));
}

Result CommAngle(const Call& C)
{
	const Matrix Rx = Rodr({1, 0, 0}, C.N(0)), Ry = Rodr({0, 1, 0}, C.N(1));
	const Matrix Cm = la::Mul(la::Mul(Rx, Ry), la::Mul(la::Transpose(Rx), la::Transpose(Ry)));
	return Result::Of(Norm3(LogSo3(Cm)));
}

Result Davenport(const Call& C)
{
	std::vector<V3> Rs, Bs;
	std::vector<Real> Ws;
	std::string E;
	if (!Observations(C, 0, 1, 2, &Rs, &Bs, &Ws, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(DavenportK(Rs, Bs, Ws, nullptr));
}

Result QMethod(const Call& C)
{
	std::vector<V3> Rs, Bs;
	std::vector<Real> Ws;
	std::string E;
	if (!Observations(C, 0, 1, 2, &Rs, &Bs, &Ws, &E))
	{
		return Result::Fail(E);
	}
	Q4 Q;
	if (!TopEigvec4(DavenportK(Rs, Bs, Ws, nullptr), &Q, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(RowQ4(Q));
}

Result WahbaSvd(const Call& C)
{
	std::vector<V3> Rs, Bs;
	std::vector<Real> Ws;
	std::string E;
	if (!Observations(C, 0, 1, 2, &Rs, &Bs, &Ws, &E))
	{
		return Result::Fail(E);
	}
	Matrix B, U, Vt;
	DavenportK(Rs, Bs, Ws, &B);
	std::vector<Real> S;
	if (!la::Svd(B, &U, &S, &Vt, &E))
	{
		return Result::Fail(E);
	}
	const Real D = Det3(U) * Det3(Vt);
	Matrix Dm = Matrix::Identity(3);
	Dm(2, 2) = D < 0.0 ? -1.0 : 1.0;
	return Result::Of(la::Mul(la::Mul(U, Dm), Vt));
}

Result WahbaLoss(const Call& C)
{
	std::vector<V3> Rs, Bs;
	std::vector<Real> Ws;
	std::string E;
	if (C.M(0).Rows != 3 || C.M(0).Cols != 3)
	{
		return Result::Fail("A must be 3 x 3");
	}
	if (!Observations(C, 1, 2, 3, &Rs, &Bs, &Ws, &E))
	{
		return Result::Fail(E);
	}
	Real L = 0.0;
	for (size_t k = 0; k < Rs.size(); ++k)
	{
		const V3 Ar = MulMV(C.M(0), Rs[k]);
		const V3 D = {Bs[k][0] - Ar[0], Bs[k][1] - Ar[1], Bs[k][2] - Ar[2]};
		L += 0.5 * Ws[k] * Dot3(D, D);
	}
	return Result::Of(L);
}

bool TriadFrame(const V3& X, const V3& Y, Matrix* F)
{
	const V3 Cxy = Cross3(X, Y);
	const Real Nx = Norm3(X), Nc = Norm3(Cxy);
	if (!(Nx > 0.0) || !(Nc > 1e-12 * Nx * Norm3(Y)))
	{
		return false;
	}
	const V3 T1 = Scale3(X, 1.0 / Nx), T2 = Scale3(Cxy, 1.0 / Nc), T3 = Cross3(T1, T2);
	*F = Matrix(3, 3);
	for (int i = 0; i < 3; ++i)
	{
		(*F)(i, 0) = T1[static_cast<size_t>(i)];
		(*F)(i, 1) = T2[static_cast<size_t>(i)];
		(*F)(i, 2) = T3[static_cast<size_t>(i)];
	}
	return true;
}

Result Triad(const Call& C)
{
	V3 R1, R2, B1, B2;
	if (!ToV3(C.M(0), &R1) || !ToV3(C.M(1), &R2) || !ToV3(C.M(2), &B1) || !ToV3(C.M(3), &B2))
	{
		return Result::Fail("needs four 3-vectors r1, r2, b1, b2");
	}
	Matrix Fr, Fb;
	if (!TriadFrame(R1, R2, &Fr) || !TriadFrame(B1, B2, &Fb))
	{
		return Result::Fail("the two vectors of a pair are parallel (or zero): rotation about them is unobservable");
	}
	return Result::Of(la::Mul(Fb, la::Transpose(Fr)));
}

Result QuatMean(const Call& C)
{
	const Matrix &Q = C.M(0), &W = C.M(1);
	if (Q.Cols != 4 || Q.Rows < 1)
	{
		return Result::Fail("q must be n x 4, one quaternion [w, x, y, z] per row");
	}
	if (!(W.IsScalar() || (W.IsVector() && W.Count() == Q.Rows)))
	{
		return Result::Fail("w must be a scalar or have one weight per row");
	}
	Matrix M(4, 4);
	for (int k = 0; k < Q.Rows; ++k)
	{
		Q4 Qk = {Q(k, 0), Q(k, 1), Q(k, 2), Q(k, 3)};
		const Real N = NormQ(Qk);
		if (!(N > 0.0))
		{
			return Result::Fail("zero quaternion in row " + std::to_string(k + 1));
		}
		const Real Wk = W.IsScalar() ? W.ScalarValue() : W.V[static_cast<size_t>(k)];
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				M(i, j) += Wk * Qk[static_cast<size_t>(i)] * Qk[static_cast<size_t>(j)] / (N * N);
			}
		}
	}
	Q4 Out;
	std::string E;
	if (!TopEigvec4(M, &Out, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(RowQ4(Out));
}

Result ExpSe3(const Call& C)
{
	V3 W, Nu;
	if (!ToV3(C.M(0), &W) || !ToV3(C.M(1), &Nu))
	{
		return Result::Fail("needs omega and nu as 3-vectors");
	}
	const Real T = Norm3(W);
	const Matrix K = Hat3(W), K2 = la::Mul(K, K);
	Real A = 0.5, B = Real(1) / 6;   // series limits of (1 - cos t)/t^2 and (t - sin t)/t^3
	Matrix R;
	if (T < 1e-6)
	{
		A = 0.5 - T * T / 24;
		B = Real(1) / 6 - T * T / 120;
		R = la::Add(la::Add(Matrix::Identity(3), K), la::Scale(K2, 0.5));
	}
	else
	{
		A = (1.0 - gm::cos(T)) / (T * T);
		B = (T - gm::sin(T)) / (T * T * T);
		R = Rodr(Scale3(W, 1.0 / T), T);
	}
	const Matrix V = la::Add(la::Add(Matrix::Identity(3), la::Scale(K, A)), la::Scale(K2, B));
	const V3 P = MulMV(V, Nu);
	Matrix G = Matrix::Identity(4);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			G(i, j) = R(i, j);
		}
		G(i, 3) = P[static_cast<size_t>(i)];
	}
	return Result::Of(G);
}

Result AttErr(const Call& C)
{
	std::string Why;
	if (!IsRotation(C.M(0), &Why) || !IsRotation(C.M(1), &Why))
	{
		return Result::Fail(Why);
	}
	const Matrix &R = C.M(0), &Rd = C.M(1);
	const Matrix D = la::Sub(la::Mul(la::Transpose(Rd), R), la::Mul(la::Transpose(R), Rd));
	// D is skew, so Vee3(D) is its exact vee.
	return Result::Of(RowV3(Scale3(Vee3(D), 0.5)));
}

// ---------------------------------------------------------------- ch02 topology

Matrix EulerE(Real Phi, Real Theta)
{
	const Real T = gm::tan(Theta), Cth = gm::cos(Theta), S = gm::sin(Phi), Cp = gm::cos(Phi);
	return Matrix::FromRows({{1.0, S * T, Cp * T}, {0.0, Cp, -S}, {0.0, S / Cth, Cp / Cth}});
}

bool GimbalOk(Real Theta)
{
	return gm::fabs(gm::cos(Theta)) > 1e-12;
}

Result EulerRates(const Call& C)
{
	const Real Phi = C.N(0), Theta = C.N(1);
	if (!GimbalOk(Theta))
	{
		return Result::Fail("gimbal lock: cos(theta) = 0, the Euler rates are unbounded");
	}
	const V3 Rates = MulMV(EulerE(Phi, Theta), {C.N(2), C.N(3), C.N(4)});
	return Result::Record({{"dphi", Rates[0]}, {"dtheta", Rates[1]}, {"dpsi", Rates[2]}});
}

Result EulerDet(const Call& C)
{
	if (!GimbalOk(C.N(0)))
	{
		return Result::Fail("gimbal lock: cos(theta) = 0");
	}
	return Result::Of(Det3(EulerE(C.N(1), C.N(0))));
}

Result EulerCond(const Call& C)
{
	if (!GimbalOk(C.N(1)))
	{
		return Result::Fail("gimbal lock: cos(theta) = 0");
	}
	std::vector<Real> S;
	std::string E;
	if (!la::Svd(EulerE(C.N(0), C.N(1)), nullptr, &S, nullptr, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(S.front() / S.back());
}

Result Eul2Dcm(const Call& C)
{
	const Matrix Rz = Rodr({0, 0, 1}, C.N(2)), Ry = Rodr({0, 1, 0}, C.N(1)), Rx = Rodr({1, 0, 0}, C.N(0));
	return Result::Of(la::Mul(la::Mul(Rz, Ry), Rx));
}

Result Dcm2Eul(const Call& C)
{
	std::string Why;
	if (!IsRotation(C.M(0), &Why))
	{
		return Result::Fail(Why);
	}
	const Matrix& R = C.M(0);
	const Real Phi = gm::atan2(R(2, 1), R(2, 2));
	const Real Theta = -gm::asin(gm::max(-1.0, gm::min(1.0, R(2, 0))));
	const Real Psi = gm::atan2(R(1, 0), R(0, 0));
	return Result::Record({{"phi", Phi}, {"theta", Theta}, {"psi", Psi}});
}

Result LookAt(const Call& C)
{
	V3 D;
	if (!ToV3(C.M(0), &D) || !(Norm3(D) > 0.0))
	{
		return Result::Fail("boresight must be a nonzero 3-vector");
	}
	D = Scale3(D, 1.0 / Norm3(D));
	const V3 U = {-D[2] * D[0], -D[2] * D[1], 1.0 - D[2] * D[2]};
	const Real N = Norm3(U);
	if (N < 1e-9)
	{
		return Result::Fail("boresight along world up: the up-vector field vanishes here (hairy ball theorem)");
	}
	return Result::Of(RowV3(Scale3(U, 1.0 / N)));
}

struct NavWorld
{
	Real K, Gx, Gy;
	int IntK;   // k when it is a whole number 1..64, else 0
	Real InvK;  // 1/k, set once per call
	Real Cx[3], Cy[3], Rr[3];
};

NavWorld MakeNavWorld(Real K, Real Gx, Real Gy, Real BigR)
{
	NavWorld W;
	W.K = K;
	W.IntK = (K == gm::floor(K) && K >= 1 && K <= 64) ? K.ToInt() : 0;
	W.InvK = 1 / K;
	W.Gx = Gx;
	W.Gy = Gy;
	const Real Cx[3] = {0.0, -0.5, 0.5}, Cy[3] = {0.0, 0.35, 0.35}, Rr[3] = {BigR, 0.15, 0.15};
	for (int i = 0; i < 3; ++i)
	{
		W.Cx[i] = Cx[i];
		W.Cy[i] = Cy[i];
		W.Rr[i] = Rr[i];
	}
	return W;
}

// Rimon-Koditschek navigation function value and gradient (ch02 bench 4); false outside free space.
bool NavEval(Real X, Real Y, const NavWorld& W, Real* F, Real* Fx, Real* Fy)
{
	const Real G = (X - W.Gx) * (X - W.Gx) + (Y - W.Gy) * (Y - W.Gy), Dgx = 2 * (X - W.Gx), Dgy = 2 * (Y - W.Gy);
	Real B = 1 - (X * X + Y * Y), Bx = -2 * X, By = -2 * Y;
	for (int i = 0; i < 3; ++i)
	{
		const Real T = (X - W.Cx[i]) * (X - W.Cx[i]) + (Y - W.Cy[i]) * (Y - W.Cy[i]) - W.Rr[i] * W.Rr[i], Tx = 2 * (X - W.Cx[i]), Ty = 2 * (Y - W.Cy[i]);
		Bx = Bx * T + B * Tx;
		By = By * T + B * Ty;
		B = B * T;
	}
	if (B <= 0)
	{
		return false;
	}
	// gamma^k and (gamma^k + beta)^(-1/k); integer k (the bench's usual setting) by multiplication.
	const Real K = W.K;
	Real Gk, Gk1, A;
	if (W.IntK > 0)
	{
		Gk1 = gm::powi(G, W.IntK - 1);
		Gk = Gk1 * G;
	}
	else
	{
		Gk = gm::pow(G, K);
		Gk1 = gm::pow(G, K - 1);
	}
	const Real D = Gk + B;
	Real Bm;   // D^(-1/k - 1)
	if (W.IntK == 1)
	{
		A = 1 / D;
		Bm = A * A;
	}
	else if (W.IntK == 2)
	{
		A = 1 / gm::sqrt(D);
		Bm = A * A * A;
	}
	else
	{
		A = gm::pow(D, -W.InvK);
		Bm = A / D;
	}
	const Real Coef = G * W.InvK * Bm;
	Gk1 = K * Gk1;
	*F = G * A;
	*Fx = Dgx * A - Coef * (Gk1 * Dgx + Bx);
	*Fy = Dgy * A - Coef * (Gk1 * Dgy + By);
	return true;
}

bool NavArgs(const Call& C, size_t First, NavWorld* W, std::string* Err)
{
	const Real K = C.N(First), BigR = C.N(First + 3);
	if (!(K >= 1.0))
	{
		*Err = "k must be >= 1";
		return false;
	}
	if (!(BigR > 0.0 && BigR < 0.45))
	{
		*Err = "big obstacle radius r must be in (0, 0.45)";
		return false;
	}
	*W = MakeNavWorld(K, C.N(First + 1), C.N(First + 2), BigR);
	return true;
}

Result NavFn(const Call& C)
{
	NavWorld W;
	std::string E;
	if (!NavArgs(C, 2, &W, &E))
	{
		return Result::Fail(E);
	}
	Real F = 0, Fx = 0, Fy = 0;
	if (!NavEval(C.N(0), C.N(1), W, &F, &Fx, &Fy))
	{
		return Result::Fail("point is outside the free space (in an obstacle or outside the unit disk)");
	}
	return Result::Record({{"phi", F}, {"fx", Fx}, {"fy", Fy}});
}

// Gradient winding per grid cell, as the chapter bench finds critical points (Poincare-Hopf: total = 1 - M).
Result NavIndex(const Call& C)
{
	NavWorld W;
	std::string E;
	if (!NavArgs(C, 0, &W, &E))
	{
		return Result::Fail(E);
	}
	const int N = C.Int(4);
	if (N < 10 || N > 400)
	{
		return Result::Fail("grid size n must be 10..400");
	}
	std::vector<Real> Gx(static_cast<size_t>(N * N), 0.0), Gy(static_cast<size_t>(N * N), 0.0);
	std::vector<char> Ok(static_cast<size_t>(N * N), 0);
	std::vector<Real> Coord(static_cast<size_t>(N));
	for (int i = 0; i < N; ++i)
	{
		Coord[static_cast<size_t>(i)] = -1 + Real(2 * i) / Real(N - 1);
	}
	for (int i = 0; i < N; ++i)
	{
		const Real Y = Coord[static_cast<size_t>(i)];
		for (int j = 0; j < N; ++j)
		{
			Real F = 0, Fx = 0, Fy = 0;
			const Real X = Coord[static_cast<size_t>(j)];
			const size_t Ix = static_cast<size_t>(i * N + j);
			if (NavEval(X, Y, W, &F, &Fx, &Fy))
			{
				Ok[Ix] = 1;
				Gx[Ix] = Fx;
				Gy[Ix] = Fy;
			}
		}
	}
	// The bench sums angle steps wrapped to (-pi, pi] and divides by 2 pi. Around a closed loop the raw
	// atan2 steps cancel, so the winding is exactly the signed count of wraps: short-way steps that cross
	// the atan2 branch cut (the negative x axis). That needs only signs and a cross product, no atan2.
	auto Upper = [&Gx, &Gy](size_t Ix) { return Gy[Ix] > 0 || (Gy[Ix] == 0 && Gx[Ix] < 0); };
	auto Wraps = [&Gx, &Gy, &Upper](size_t Ia, size_t Ic) {
		const Real Cr = Gx[Ia] * Gy[Ic] - Gy[Ia] * Gx[Ic];
		const bool Ua = Upper(Ia), Uc = Upper(Ic);
		if (!Ua && Uc && Cr < 0)
		{
			return -1;   // raw step > pi, wrapped by -2 pi
		}
		if (Ua && !Uc && Cr > 0)
		{
			return 1;    // raw step < -pi, wrapped by +2 pi
		}
		return 0;
	};
	int Total = 0, Plus = 0, Minus = 0;
	auto At = [N](int i, int j) { return static_cast<size_t>(i * N + j); };
	for (int i = 0; i < N - 1; ++i)
	{
		for (int j = 0; j < N - 1; ++j)
		{
			if (!(Ok[At(i, j)] && Ok[At(i, j + 1)] && Ok[At(i + 1, j + 1)] && Ok[At(i + 1, j)]))
			{
				continue;
			}
			const int Wnd = Wraps(At(i, j), At(i, j + 1)) + Wraps(At(i, j + 1), At(i + 1, j + 1)) + Wraps(At(i + 1, j + 1), At(i + 1, j)) + Wraps(At(i + 1, j), At(i, j));
			Total += Wnd;
			Plus += Wnd > 0 ? 1 : 0;
			Minus += Wnd < 0 ? 1 : 0;
		}
	}
	return Result::Record({{"total", Real(Total)}, {"plus", Real(Plus)}, {"minus", Real(Minus)}});
}

// mulberry32, bit-for-bit as in the chapter pages, and their Box-Muller gauss().
struct Mulberry
{
	uint32_t S;
	Real Next()
	{
		S += 0x6D2B79F5u;
		uint32_t T = S;
		T = (T ^ (T >> 15)) * (T | 1u);
		T ^= T + (T ^ (T >> 7)) * (T | 61u);
		return Real(T ^ (T >> 14)) / Real(4294967296ull);
	}
	Real Gauss()
	{
		Real U = 0.0;
		while (U == 0.0)
		{
			U = Next();
		}
		return gm::sqrt(-2 * gm::log(U)) * gm::cos(2 * cst::Pi() * Next());
	}
};

Real Len4(const Q4& Q) { return gm::sqrt(Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]); }

// The unwinding bench: unit-inertia body, quaternion PD feedback, three laws (ch02 counterexample 2).
Result Unwind(const Call& C)
{
	const Real Th0 = C.N(0), Sigma = C.N(2), W0 = C.N(3);
	const int Law = C.Int(1);
	if (Law < 0 || Law > 2)
	{
		return Result::Fail("law: 0 = naive continuous, 1 = sign switch, 2 = hybrid with hysteresis");
	}
	if (!(Sigma >= 0.0))
	{
		return Result::Fail("sigma must be >= 0");
	}
	V3 Axis = {0.3, -0.4, 0.866};
	Axis = Scale3(Axis, 1.0 / Norm3(Axis));
	const Real Kp = 1, Kd = 1.4, Dt = 0.01, Tend = 25, Hyst = 0.2;
	Q4 Q = {gm::cos(Th0 / 2), Axis[0] * gm::sin(Th0 / 2), Axis[1] * gm::sin(Th0 / 2), Axis[2] * gm::sin(Th0 / 2)};
	V3 W = Scale3(Axis, W0);
	Real H = 1, Travelled = 0, LastSign = 0;
	int Switches = 0;
	bool HaveLast = false;
	Mulberry Rng{99u};
	const Real Steps = Tend / Dt;
	for (int k = 0; Real(k) <= Steps; ++k)
	{
		Q4 Qm;
		for (size_t i = 0; i < 4; ++i)
		{
			Qm[i] = Q[i] + Sigma * Rng.Gauss();
		}
		const Real Nm = Len4(Qm);
		for (Real& X : Qm)
		{
			X /= Nm;
		}
		Real S = 1;
		if (Law == 1)
		{
			S = Qm[0] >= 0 ? 1 : -1;
		}
		else if (Law == 2)
		{
			if (H * Qm[0] < -Hyst)
			{
				H = -H;
			}
			S = H;
		}
		if (HaveLast && S != LastSign)
		{
			++Switches;
		}
		LastSign = S;
		HaveLast = true;
		V3 Tau;
		for (size_t i = 0; i < 3; ++i)
		{
			Tau[i] = Qm[i + 1] * (-Kp * S) + W[i] * (-Kd);
		}
		Travelled += Norm3(W) * Dt;
		for (size_t i = 0; i < 3; ++i)
		{
			W[i] = W[i] + Tau[i] * Dt;
		}
		const Q4 Qd = QuatProduct(Q, {0.0, W[0], W[1], W[2]});
		Q4 Qn;
		for (size_t i = 0; i < 4; ++i)
		{
			Qn[i] = Q[i] + 0.5 * Qd[i] * Dt;
		}
		const Real Nq = Len4(Qn);
		for (size_t i = 0; i < 4; ++i)
		{
			Q[i] = Qn[i] / Nq;
		}
	}
	const Real Final = 2 * gm::acos(gm::min(1.0, gm::fabs(Q[0])));
	return Result::Record({{"travelled", Travelled}, {"switches", Real(Switches)}, {"final", Final}});
}

// ---------------------------------------------------------------- ch04 Hilbert spaces

// A vector argument as a column (rows or columns accepted).
bool AsColumn(const Matrix& M, int N, Matrix* Out)
{
	if (!(M.IsVector() || M.IsScalar()) || (N >= 0 && M.Count() != N))
	{
		return false;
	}
	*Out = Matrix(M.Count(), 1);
	Out->V = M.V;
	return true;
}

Result Proj(const Call& C)
{
	const Matrix& A = C.M(0);
	Matrix Y, X;
	std::string E;
	if (A.IsEmpty() || !AsColumn(C.M(1), A.Rows, &Y))
	{
		return Result::Fail("y must be a vector with one entry per row of A");
	}
	if (la::Rank(A) < A.Cols)
	{
		return Result::Fail("the columns of A are linearly dependent");
	}
	if (!la::LeastSquares(A, Y, &X, &E))
	{
		return Result::Fail(E);
	}
	Matrix P = la::Mul(A, X);
	if (C.M(1).Rows == 1)
	{
		P = la::Transpose(P);
	}
	return Result::Of(P);
}

Result GramSchmidt(const Call& C)
{
	Matrix Q = C.M(0);
	if (Q.IsEmpty())
	{
		return Result::Fail("needs a matrix");
	}
	Real Scale = 0.0;
	for (Real X : Q.V)
	{
		Scale = gm::max(Scale, gm::fabs(X));
	}
	// Modified Gram-Schmidt, columns in order.
	for (int j = 0; j < Q.Cols; ++j)
	{
		for (int k = 0; k < j; ++k)
		{
			Real D = 0.0;
			for (int i = 0; i < Q.Rows; ++i)
			{
				D += Q(i, k) * Q(i, j);
			}
			for (int i = 0; i < Q.Rows; ++i)
			{
				Q(i, j) -= D * Q(i, k);
			}
		}
		Real N = 0.0;
		for (int i = 0; i < Q.Rows; ++i)
		{
			N += Q(i, j) * Q(i, j);
		}
		N = gm::sqrt(N);
		if (!(N > 1e-12 * gm::max(1.0, Scale)))
		{
			return Result::Fail("column " + std::to_string(j + 1) + " depends on the earlier columns");
		}
		for (int i = 0; i < Q.Rows; ++i)
		{
			Q(i, j) /= N;
		}
	}
	return Result::Of(Q);
}

Result Wls(const Call& C)
{
	const Matrix& A = C.M(0);
	const int M = A.Rows;
	Matrix Y, X, L;
	std::string E;
	if (A.IsEmpty() || !AsColumn(C.M(1), M, &Y))
	{
		return Result::Fail("y must be a vector with one entry per row of A");
	}
	const Matrix& R = C.M(2);
	Matrix Rm;
	if (R.IsScalar())
	{
		Rm = la::Scale(Matrix::Identity(M), R.ScalarValue());
	}
	else if (R.IsVector() && R.Count() == M)
	{
		Rm = Matrix(M, M);
		for (int i = 0; i < M; ++i)
		{
			Rm(i, i) = R.V[static_cast<size_t>(i)];
		}
	}
	else if (R.Rows == M && R.Cols == M)
	{
		Rm = R;
	}
	else
	{
		return Result::Fail("R must be a scalar, m variances, or an m x m covariance");
	}
	if (!la::Cholesky(Rm, &L, &E))
	{
		return Result::Fail("R must be symmetric positive definite");
	}
	// Whiten: L^-1 A x ~ L^-1 y, then ordinary least squares.
	Matrix Aw, Yw;
	if (!la::Solve(L, A, &Aw, &E) || !la::Solve(L, Y, &Yw, &E) || !la::LeastSquares(Aw, Yw, &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}

Result Vander(const Call& C)
{
	const Matrix& Xs = C.M(0);
	const int D = C.Int(1);
	if (!(Xs.IsVector() || Xs.IsScalar()) || D < 0 || D > 30)
	{
		return Result::Fail("needs a vector of points and a degree 0..30");
	}
	Matrix A(Xs.Count(), D + 1);
	for (int i = 0; i < A.Rows; ++i)
	{
		Real P = 1.0;
		for (int j = 0; j <= D; ++j)
		{
			A(i, j) = P;
			P *= Xs.V[static_cast<size_t>(i)];
		}
	}
	return Result::Of(A);
}

bool RwArgs(const Call& C, std::string* Err)
{
	if (!(C.M(0).IsVector() || C.M(0).IsScalar()))
	{
		*Err = "y must be a vector of measurements";
		return false;
	}
	if (!(C.N(1) > 0.0 && C.N(2) > 0.0 && C.N(3) > 0.0))
	{
		*Err = "q, r and P0 must be positive";
		return false;
	}
	return true;
}

Result KfRw(const Call& C)
{
	std::string E;
	if (!RwArgs(C, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<Real>& Ys = C.M(0).V;
	const Real Qn = C.N(1), Rn = C.N(2);
	Real Xh = 0.0, P = C.N(3);
	Matrix Out(1, static_cast<int>(Ys.size()));
	for (size_t k = 0; k < Ys.size(); ++k)
	{
		if (k > 0)
		{
			P += Qn;
		}
		const Real K = P / (P + Rn);
		Xh += K * (Ys[k] - Xh);
		P *= 1 - K;
		Out.V[k] = Xh;
	}
	return Result::Of(Out);
}

// Batch least squares over x_0..x_{n-1}: the tridiagonal information system, by the Thomas algorithm.
Result SmoothRw(const Call& C)
{
	std::string E;
	if (!RwArgs(C, &E))
	{
		return Result::Fail(E);
	}
	const std::vector<Real>& Ys = C.M(0).V;
	const Real Qn = C.N(1), Rn = C.N(2), P0 = C.N(3);
	const size_t N = Ys.size();
	if (N == 0)
	{
		return Result::Fail("y must not be empty");
	}
	std::vector<Real> A(N, 0.0), B(N, 0.0), Cc(N, 0.0), D(N, 0.0), Cp(N), Dp(N), X(N);
	B[0] += 1 / P0;
	for (size_t i = 1; i < N; ++i)
	{
		B[i] += 1 / Qn;
		B[i - 1] += 1 / Qn;
		A[i] -= 1 / Qn;
		Cc[i - 1] -= 1 / Qn;
	}
	for (size_t i = 0; i < N; ++i)
	{
		B[i] += 1 / Rn;
		D[i] += Ys[i] / Rn;
	}
	Cp[0] = Cc[0] / B[0];
	Dp[0] = D[0] / B[0];
	for (size_t i = 1; i < N; ++i)
	{
		const Real Den = B[i] - A[i] * Cp[i - 1];
		Cp[i] = Cc[i] / Den;
		Dp[i] = (D[i] - A[i] * Dp[i - 1]) / Den;
	}
	X[N - 1] = Dp[N - 1];
	for (size_t i = N - 1; i-- > 0;)
	{
		X[i] = Dp[i] - Cp[i] * X[i + 1];
	}
	return Result::Of(Matrix::Row(X));
}

bool Ar1Args(Real A, Real Q, Real R, std::string* Err)
{
	if (!(gm::fabs(A) < 1.0) || !(Q > 0.0) || !(R > 0.0))
	{
		*Err = "needs |a| < 1 (stationary AR(1)), q > 0 and r > 0";
		return false;
	}
	return true;
}

// Steady-state prior variance of the scalar Kalman filter: the positive root of P^2 + (r(1-a^2) - q)P - qr = 0.
Real SteadyP(Real A, Real Q, Real R)
{
	const Real Bq = R * (1 - A * A) - Q;
	return 0.5 * (-Bq + gm::sqrt(Bq * Bq + 4 * Q * R));
}

Result KfSs(const Call& C)
{
	const Real A = C.N(0), Q = C.N(1), R = C.N(2);
	std::string E;
	if (!Ar1Args(A, Q, R, &E))
	{
		return Result::Fail(E);
	}
	const Real P = SteadyP(A, Q, R);
	return Result::Record({{"P", P}, {"K", P / (P + R)}, {"Pp", P * R / (P + R)}});
}

// Non-causal Wiener MSE (1/pi) int_0^pi Ss r / (Ss + r) dw = q r / sqrt((q + r(1+a^2))^2 - 4 a^2 r^2).
Result Wiener(const Call& C)
{
	const Real A = C.N(0), Q = C.N(1), R = C.N(2);
	std::string E;
	if (!Ar1Args(A, Q, R, &E))
	{
		return Result::Fail(E);
	}
	const Real S = Q + R * (1 + A * A);
	const Real Nc = Q * R / gm::sqrt(S * S - 4 * A * A * R * R);
	const Real P = SteadyP(A, Q, R), Pp = P * R / (P + R);
	return Result::Record({{"nc", Nc}, {"causal", Pp}, {"price", Pp / Nc - 1}});
}

Real SMagBode(Real W, Real K, Real P)
{
	// L = k (s + 1) / ((s - p)(s + 3)(s + 4)) at s = jw; returns |1 / (1 + L)|.
	const Real Nr = K, Ni = K * W;                                  // k (1 + jw)
	const Real Ar = -P, Ai = W, Br = 3, Bi = W, Cr = 4, Ci = W;
	const Real ABr = Ar * Br - Ai * Bi, ABi = Ar * Bi + Ai * Br;
	const Real Dr = ABr * Cr - ABi * Ci, Di = ABr * Ci + ABi * Cr;  // (s-p)(s+3)(s+4)
	const Real Sr = Dr + Nr, Si = Di + Ni;                           // S = D / (D + N)
	return gm::hypot(Dr, Di) / gm::hypot(Sr, Si);
}

// Bode's sensitivity integral computed numerically (Gauss-Legendre in log frequency) beside pi * p.
Result BodeInt(const Call& C)
{
	const Real K = C.N(0), P = C.N(1);
	if (!(K > 0.0))
	{
		return Result::Fail("gain k must be positive");
	}
	// Closed-loop characteristic polynomial s^3 + (7 - p)s^2 + (12 - 7p + k)s + (k - 12p).
	Matrix Comp(3, 3);
	Comp(0, 0) = -(7 - P);
	Comp(0, 1) = -(12 - 7 * P + K);
	Comp(0, 2) = -(K - 12 * P);
	Comp(1, 0) = 1;
	Comp(2, 1) = 1;
	std::vector<Real> Re, Im;
	std::string E;
	if (!la::Eigenvalues(Comp, &Re, &Im, &E))
	{
		return Result::Fail(E);
	}
	const bool Stable = *std::max_element(Re.begin(), Re.end()) < 0.0;
	static const Real Gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
	static const Real Gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066538101824, 0.3626837833783620, 0.3626837833783620, 0.3137066538101824, 0.2223810344533745, 0.1012285362903763};
	const Real U0 = gm::log(1e-12), U1 = gm::log(1e8);
	const int Panels = 1200;
	const Real Hp = (U1 - U0) / Panels;
	Real I = 0.0, Peak = 0.0, Wpk = 0.0;
	for (int p = 0; p < Panels; ++p)
	{
		const Real Mid = U0 + (p + 0.5) * Hp;
		for (int g = 0; g < 8; ++g)
		{
			const Real W = gm::exp(Mid + 0.5 * Hp * Gx[g]);
			const Real Sm = SMagBode(W, K, P);
			I += 0.5 * Hp * Gw[g] * gm::log(Sm) * W;
			if (Sm > Peak)
			{
				Peak = Sm;
				Wpk = W;
			}
		}
	}
	I += K / gm::exp(U1);   // tail: ln|S| ~ k / w^2 beyond the last panel
	// Refine the peak by golden-section search around the best sample.
	Real Lo = Wpk * gm::exp(-2 * Hp), Hi = Wpk * gm::exp(2 * Hp);
	const Real Gr = 0.5 * (gm::sqrt(5.0) - 1.0);
	for (int it = 0; it < 100; ++it)
	{
		const Real A1 = Hi - Gr * (Hi - Lo), A2 = Lo + Gr * (Hi - Lo);
		if (SMagBode(A1, K, P) > SMagBode(A2, K, P))
		{
			Hi = A2;
		}
		else
		{
			Lo = A1;
		}
	}
	Wpk = 0.5 * (Lo + Hi);
	Peak = gm::max(Peak, SMagBode(Wpk, K, P));
	return Result::Record({{"integral", I}, {"theorem", cst::Pi() * gm::max(P, 0.0)}, {"peak", Peak}, {"wpeak", Wpk}, {"stable", Stable ? 1.0 : 0.0}});
}

// Posterior mean of x ~ N(0, 4) from measurements with N(0,1)/N(0,100) mixture noise (801-point grid),
// beside the best linear estimate: ch04 counterexample 1.
Result PostMean(const Call& C)
{
	const Matrix& Y = C.M(0);
	const Real Eps = C.N(1);
	if (!(Y.IsVector() || Y.IsScalar()))
	{
		return Result::Fail("y must be a vector of measurements");
	}
	if (!(Eps >= 0.0 && Eps <= 1.0))
	{
		return Result::Fail("outlier rate eps must be in [0, 1]");
	}
	const int G = 801;
	std::vector<Real> Lp(G);
	Real Best = -Real::Infinity();
	for (int i = 0; i < G; ++i)
	{
		const Real X = -10 + Real(20 * i) / 800;
		Real S = -X * X / 8;
		for (Real Yk : Y.V)
		{
			const Real D = Yk - X;
			S += gm::log((1 - Eps) * gm::exp(-D * D / 2) + Eps * gm::exp(-D * D / 200) / 10 + 1e-300);
		}
		Lp[static_cast<size_t>(i)] = S;
		Best = gm::max(Best, S);
	}
	Real Sw = 0, Sx = 0;
	for (int i = 0; i < G; ++i)
	{
		const Real Wt = gm::exp(Lp[static_cast<size_t>(i)] - Best);
		Sw += Wt;
		Sx += Wt * (-10 + Real(20 * i) / 800);
	}
	Real Sum = 0;
	for (Real Yk : Y.V)
	{
		Sum += Yk;
	}
	const Real Vv = (1 - Eps) + Eps * 100;
	const Real Lin = (Sum / Vv) / (0.25 + Real(Y.Count()) / Vv);
	return Result::Record({{"post", Sx / Sw}, {"lin", Lin}, {"linmse", 1 / (0.25 + Real(Y.Count()) / Vv)}});
}

// ---------------------------------------------------------------- ch05 measure

Result Bayes(const Call& C)
{
	const Real M = C.N(0), P = C.N(1), Y = C.N(2), R = C.N(3), H = C.N(4);
	if (!(P >= 0.0) || !(R > 0.0))
	{
		return Result::Fail("needs prior variance P >= 0 and noise variance R > 0");
	}
	const Real K = P * H / (H * H * P + R);
	return Result::Record({{"m", M + K * (Y - H * M)}, {"P", (1 - K * H) * P}, {"K", K}});
}

bool GaussUpdate(const Call& C, bool WantMean, Matrix* Out, std::string* Err)
{
	// bayes_mean(x, P, y, H, R) or bayes_cov(P, H, R)
	const size_t Ip = WantMean ? 1 : 0, Ih = WantMean ? 3 : 1, Ir = WantMean ? 4 : 2;
	const Matrix &P = C.M(Ip), &H = C.M(Ih), &R = C.M(Ir);
	const int N = P.Rows, M = H.Rows;
	if (!P.IsSquare() || H.Cols != N || R.Rows != M || R.Cols != M)
	{
		*Err = "shapes: P n x n, H m x n, R m x m";
		return false;
	}
	const Matrix Ht = la::Transpose(H);
	const Matrix S = la::Add(la::Mul(la::Mul(H, P), Ht), R);
	Matrix Si;
	if (!la::Inverse(S, &Si, Err))
	{
		return false;
	}
	const Matrix K = la::Mul(la::Mul(P, Ht), Si);
	if (!WantMean)
	{
		*Out = la::Symmetrize(la::Sub(P, la::Mul(la::Mul(K, S), la::Transpose(K))));
		return true;
	}
	Matrix X, Y;
	if (!AsColumn(C.M(0), N, &X) || !AsColumn(C.M(2), M, &Y))
	{
		*Err = "x must have n entries and y m entries";
		return false;
	}
	*Out = la::Add(X, la::Mul(K, la::Sub(Y, la::Mul(H, X))));
	return true;
}

Result BayesMean(const Call& C)
{
	Matrix Out;
	std::string E;
	if (!GaussUpdate(C, true, &Out, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Out);
}

Result BayesCov(const Call& C)
{
	Matrix Out;
	std::string E;
	if (!GaussUpdate(C, false, &Out, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Out);
}

Result BkDensity(const Call& C)
{
	const Real T = C.N(0);
	return Result::Record({{"band", gm::fabs(T) <= cst::Pi() ? 1.0 / (2 * cst::Pi()) : 0.0}, {"wedge", gm::fabs(T) <= 0.5 * cst::Pi() ? 0.5 * gm::cos(T) : 0.0}});
}

// p = N(0,1), q = N(0, s^2): E_q[w^2] = s^2 / sqrt(2 s^2 - 1), infinite for s^2 <= 1/2.
Result IsW2(const Call& C)
{
	const Real S = C.N(0);
	if (!(S > 0.0))
	{
		return Result::Fail("proposal sd s must be positive");
	}
	if (2 * S * S - 1 <= 0.0)
	{
		return Result::Record({{"w2", Real::Infinity()}, {"ess", 0.0}, {"fvar", Real::Infinity()}});
	}
	const Real W2 = S * S / gm::sqrt(2 * S * S - 1);
	const Real V = S * S / (2 * S * S - 1);   // variance of the Gaussian p^2/q, normalized
	return Result::Record({{"w2", W2}, {"ess", 1 / W2}, {"fvar", 3 * W2 * V * V - 1}});
}

// Radical inverse: the base-B digits of I reversed behind the point. While B^digits fits in 18
// decimal digits the reversed digits are an exact integer ratio and one division gives the value.
Real RadInv(long long I, long long B)
{
	long long Num = 0, Den = 1, J = I;
	while (J > 0 && Den <= 1000000000000000000ll / B)
	{
		Num = Num * B + J % B;
		Den *= B;
		J /= B;
	}
	if (J == 0)
	{
		return Real(Num) / Real(Den);
	}
	Real F = 1, R = 0;
	while (I > 0)
	{
		F /= Real(B);
		R += F * Real(I % B);
		I /= B;
	}
	return R;
}

const long long AtlasPrimes[16] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};

Result Vdc(const Call& C)
{
	const Real I = C.N(0), B = C.N(1);
	if (!(I >= 0 && I == gm::floor(I) && I < 9e15) || !(B >= 2 && B == gm::floor(B) && B <= 1e6))
	{
		return Result::Fail("needs a whole index i >= 0 and a whole base b >= 2");
	}
	return Result::Of(RadInv(I.ToLong(), B.ToLong()));
}

Result Halton(const Call& C)
{
	const int N = C.Int(0), D = C.Int(1);
	if (N < 1 || N > 20000 || D < 1 || D > 16 || N * D > 80000)
	{
		return Result::Fail("n must be 1..20000 and d 1..16 (n d <= 80000)");
	}
	Matrix H(N, D);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < D; ++j)
		{
			H(i, j) = RadInv(i + 1, AtlasPrimes[j]);
		}
	}
	return Result::Of(H);
}

// |Halton estimate - exact| for the chapter's smooth (kind 1) or discontinuous (kind 0) integrand.
Result QmcErr(const Call& C)
{
	const int N = C.Int(0), D = C.Int(1), Kind = C.Int(2);
	if (N < 1 || N > 100000 || D < 1 || D > 16 || N * D > 100000 || (Kind != 0 && Kind != 1))
	{
		return Result::Fail("N 1..100000, d 1..16 (N d <= 100000), kind 1 = smooth product of sines, 0 = indicator of sum < d/2");
	}
	const Real& Pi = cst::Pi();
	const Real Half = Real(D) / 2, Scale = gm::powi(Pi / 2, D);
	RealAccumulator Acc;   // the 8192-term sum is kept in the accumulator and read off once
	long long Inside = 0;
	for (int i = 1; i <= N; ++i)
	{
		if (Kind == 1)
		{
			Real P = Scale;
			for (int j = 0; j < D; ++j)
			{
				P *= gm::sin(Pi * RadInv(i, AtlasPrimes[j]));
			}
			Acc.Add(P);
		}
		else
		{
			Real S = 0;
			for (int j = 0; j < D; ++j)
			{
				S += RadInv(i, AtlasPrimes[j]);
			}
			Inside += S < Half ? 1 : 0;
		}
	}
	const Real Mean = Kind == 1 ? Acc.Result() / N : Real(Inside) / N;
	return Result::Of(gm::fabs(Mean - (Kind == 1 ? Real(1) : Real(0.5))));
}

Result Ess(const Call& C)
{
	Real S = 0, S2 = 0;
	for (Real W : C.M(0).V)
	{
		if (!(W >= 0.0))
		{
			return Result::Fail("weights must be non-negative");
		}
		S += W;
		S2 += W * W;
	}
	if (!(S2 > 0.0))
	{
		return Result::Fail("all weights are zero");
	}
	return Result::Of(S * S / S2);
}

Real LogZ(const std::vector<Real>& E)
{
	Real Lo = Real::Infinity();
	for (Real X : E)
	{
		Lo = gm::min(Lo, X);
	}
	Real S = 0;
	for (Real X : E)
	{
		S += gm::exp(-(X - Lo));
	}
	return -Lo + gm::log(S);
}

Result LogZFn(const Call& C)
{
	if (C.M(0).IsEmpty())
	{
		return Result::Fail("needs energies");
	}
	return Result::Of(LogZ(C.M(0).V));
}

Result Boltz(const Call& C)
{
	if (C.M(0).IsEmpty())
	{
		return Result::Fail("needs energies");
	}
	Matrix Q = C.M(0);
	const Real Lz = LogZ(Q.V);
	for (Real& X : Q.V)
	{
		X = gm::exp(-X - Lz);
	}
	return Result::Of(Q);
}

bool Distribution(const Matrix& M, std::vector<Real>* Out, std::string* Err)
{
	Real S = 0;
	for (Real X : M.V)
	{
		if (!(X >= 0.0))
		{
			*Err = "probabilities must be non-negative";
			return false;
		}
		S += X;
	}
	if (!(S > 0.0))
	{
		*Err = "probabilities sum to zero";
		return false;
	}
	for (Real X : M.V)
	{
		Out->push_back(X / S);
	}
	return true;
}

// Gibbs functional -E_q[E] + H(q); its maximum over q is log Z, attained at the Boltzmann distribution.
Result GibbsF(const Call& C)
{
	std::vector<Real> Q;
	std::string Er;
	if (!Distribution(C.M(0), &Q, &Er))
	{
		return Result::Fail(Er);
	}
	if (C.M(1).Count() != static_cast<int>(Q.size()))
	{
		return Result::Fail("q and E must have the same length");
	}
	Real F = 0;
	for (size_t i = 0; i < Q.size(); ++i)
	{
		F -= Q[i] * C.M(1).V[i];
		if (Q[i] > 0)
		{
			F -= Q[i] * gm::log(Q[i]);
		}
	}
	return Result::Of(F);
}

Result Kl(const Call& C)
{
	std::vector<Real> P, Q;
	std::string Er;
	if (!Distribution(C.M(0), &P, &Er) || !Distribution(C.M(1), &Q, &Er))
	{
		return Result::Fail(Er);
	}
	if (P.size() != Q.size())
	{
		return Result::Fail("p and q must have the same length");
	}
	Real D = 0;
	for (size_t i = 0; i < P.size(); ++i)
	{
		if (P[i] > 0)
		{
			if (Q[i] == 0)
			{
				return Result::Of(Real::Infinity());
			}
			D += P[i] * gm::log(P[i] / Q[i]);
		}
	}
	return Result::Of(D);
}

// ---------------------------------------------------------------- ch06 optimization

// Dense two-phase simplex (Bland's rule) for min c'x s.t. A x <= b, x >= 0.
// Dual prices y >= 0 satisfy c + A'y >= 0 and c'x* = -b'y*; y_i = -d(c'x*)/d(b_i).
struct LpSolution
{
	std::vector<Real> X, Y;
	Real Value = 0;
};

bool SolveLp(const Matrix& Cm, const Matrix& A, const Matrix& Bm, LpSolution* Out, std::string* Err)
{
	const int M = A.Rows, N = A.Cols;
	if (A.IsEmpty() || Cm.Count() != N || Bm.Count() != M || !(Cm.IsVector() || Cm.IsScalar()) || !(Bm.IsVector() || Bm.IsScalar()))
	{
		*Err = "shapes: c has n entries, A is m x n, b has m entries";
		return false;
	}
	if (M > 60 || N > 60)
	{
		*Err = "at most 60 constraints and 60 variables";
		return false;
	}
	std::vector<int> ArtRow;
	for (int i = 0; i < M; ++i)
	{
		if (Bm.V[static_cast<size_t>(i)] < 0)
		{
			ArtRow.push_back(i);
		}
	}
	const int NA = static_cast<int>(ArtRow.size()), Cols = N + M + NA;
	// Tableau rows 0..M-1, objective row M; column Cols holds the right-hand side.
	Matrix T(M + 1, Cols + 1);
	std::vector<int> Basis(static_cast<size_t>(M));
	int ArtK = 0;
	for (int i = 0; i < M; ++i)
	{
		const bool Flip = Bm.V[static_cast<size_t>(i)] < 0;
		const Real Sg = Flip ? -1.0 : 1.0;
		for (int j = 0; j < N; ++j)
		{
			T(i, j) = Sg * A(i, j);
		}
		T(i, N + i) = Sg;
		T(i, Cols) = Sg * Bm.V[static_cast<size_t>(i)];
		if (Flip)
		{
			T(i, N + M + ArtK) = 1.0;
			Basis[static_cast<size_t>(i)] = N + M + ArtK;
			++ArtK;
		}
		else
		{
			Basis[static_cast<size_t>(i)] = N + i;
		}
	}
	auto Pivot = [&T, &Basis, M, Cols](int Pr, int Pc) {
		const Real Pv = T(Pr, Pc);
		for (int j = 0; j <= Cols; ++j)
		{
			T(Pr, j) /= Pv;
		}
		for (int i = 0; i <= M; ++i)
		{
			if (i != Pr && T(i, Pc) != 0.0)
			{
				const Real F = T(i, Pc);
				for (int j = 0; j <= Cols; ++j)
				{
					T(i, j) -= F * T(Pr, j);
				}
			}
		}
		Basis[static_cast<size_t>(Pr)] = Pc;
	};
	const Real Tol = 1e-11;
	// Runs simplex on the objective row; columns >= Allowed may not enter. Returns 0 ok, 1 unbounded.
	auto Run = [&T, &Basis, &Pivot, M, Tol, Cols](int Allowed) {
		for (int Iter = 0; Iter < 5000; ++Iter)
		{
			int Pc = -1;
			for (int j = 0; j < Allowed; ++j)
			{
				if (T(M, j) < -Tol)
				{
					Pc = j;
					break;
				}
			}
			if (Pc < 0)
			{
				return 0;
			}
			int Pr = -1;
			Real BestRatio = 0;
			for (int i = 0; i < M; ++i)
			{
				if (T(i, Pc) > Tol)
				{
					const Real Ratio = T(i, Cols) / T(i, Pc);
					if (Pr < 0 || Ratio < BestRatio - 1e-13 || (gm::fabs(Ratio - BestRatio) <= 1e-13 && Basis[static_cast<size_t>(i)] < Basis[static_cast<size_t>(Pr)]))
					{
						Pr = i;
						BestRatio = Ratio;
					}
				}
			}
			if (Pr < 0)
			{
				return 1;
			}
			Pivot(Pr, Pc);
		}
		return 2;
	};
	if (NA > 0)
	{
		// Phase 1: minimize the sum of artificials, priced out against the starting basis.
		for (int j = 0; j <= Cols; ++j)
		{
			T(M, j) = 0.0;
		}
		for (int k = 0; k < NA; ++k)
		{
			const int Row = ArtRow[static_cast<size_t>(k)];
			for (int j = 0; j <= Cols; ++j)
			{
				T(M, j) -= T(Row, j);
			}
			T(M, N + M + k) = 0.0;
		}
		if (Run(Cols) != 0 || -T(M, Cols) > 1e-9 * gm::max(1.0, la::MaxAbs(Bm)))
		{
			*Err = "infeasible: no x >= 0 satisfies A x <= b";
			return false;
		}
		// Drive artificials out of the basis where possible.
		for (int i = 0; i < M; ++i)
		{
			if (Basis[static_cast<size_t>(i)] >= N + M)
			{
				for (int j = 0; j < N + M; ++j)
				{
					if (gm::fabs(T(i, j)) > 1e-9)
					{
						Pivot(i, j);
						break;
					}
				}
			}
		}
	}
	// Phase 2: reduced costs of c'x against the current basis.
	for (int j = 0; j <= Cols; ++j)
	{
		T(M, j) = j < N ? Cm.V[static_cast<size_t>(j)] : 0.0;
	}
	for (int i = 0; i < M; ++i)
	{
		const int Bj = Basis[static_cast<size_t>(i)];
		const Real Cb = Bj < N ? Cm.V[static_cast<size_t>(Bj)] : 0.0;
		if (Cb != 0.0)
		{
			for (int j = 0; j <= Cols; ++j)
			{
				T(M, j) -= Cb * T(i, j);
			}
		}
	}
	const int Status = Run(N + M);
	if (Status == 1)
	{
		*Err = "unbounded: c'x has no minimum on the feasible set";
		return false;
	}
	if (Status != 0)
	{
		*Err = "simplex did not converge";
		return false;
	}
	Out->X.assign(static_cast<size_t>(N), 0.0);
	for (int i = 0; i < M; ++i)
	{
		const int Bj = Basis[static_cast<size_t>(i)];
		if (Bj < N)
		{
			Out->X[static_cast<size_t>(Bj)] = T(i, Cols);
		}
	}
	Out->Y.assign(static_cast<size_t>(M), 0.0);
	for (int i = 0; i < M; ++i)
	{
		Out->Y[static_cast<size_t>(i)] = gm::max(0.0, T(M, N + i));   // reduced cost of slack i
	}
	Out->Value = 0;
	for (int j = 0; j < N; ++j)
	{
		Out->Value += Cm.V[static_cast<size_t>(j)] * Out->X[static_cast<size_t>(j)];
	}
	return true;
}

Result Lp(const Call& C)
{
	LpSolution S;
	std::string E;
	if (!SolveLp(C.M(0), C.M(1), C.M(2), &S, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Matrix::Column(S.X));
}

Result LpVal(const Call& C)
{
	LpSolution S;
	std::string E;
	if (!SolveLp(C.M(0), C.M(1), C.M(2), &S, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(S.Value);
}

Result LpDual(const Call& C)
{
	LpSolution S;
	std::string E;
	if (!SolveLp(C.M(0), C.M(1), C.M(2), &S, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Matrix::Column(S.Y));
}

// KKT residuals of (x, y) for min c'x, A x <= b, x >= 0 with multipliers y >= 0.
Result LpKkt(const Call& C)
{
	const Matrix &Cm = C.M(0), &A = C.M(1), &B = C.M(2), &X = C.M(3), &Y = C.M(4);
	const int M = A.Rows, N = A.Cols;
	if (Cm.Count() != N || B.Count() != M || X.Count() != N || Y.Count() != M)
	{
		return Result::Fail("shapes: c and x have n entries, A is m x n, b and y have m entries");
	}
	Real Primal = 0, Dual = 0, Slack = 0, Gap = 0;
	for (int i = 0; i < M; ++i)
	{
		Real Ax = 0;
		for (int j = 0; j < N; ++j)
		{
			Ax += A(i, j) * X.V[static_cast<size_t>(j)];
		}
		const Real Yi = Y.V[static_cast<size_t>(i)];
		Primal = gm::max(Primal, Ax - B.V[static_cast<size_t>(i)]);
		Dual = gm::max(Dual, -Yi);
		Slack = gm::max(Slack, gm::fabs(Yi * (B.V[static_cast<size_t>(i)] - Ax)));
		Gap += B.V[static_cast<size_t>(i)] * Yi;
	}
	for (int j = 0; j < N; ++j)
	{
		Real Mu = Cm.V[static_cast<size_t>(j)];
		for (int i = 0; i < M; ++i)
		{
			Mu += A(i, j) * Y.V[static_cast<size_t>(i)];
		}
		const Real Xj = X.V[static_cast<size_t>(j)];
		Primal = gm::max(Primal, -Xj);
		Dual = gm::max(Dual, -Mu);
		Slack = gm::max(Slack, gm::fabs(Xj * Mu));
		Gap += Cm.V[static_cast<size_t>(j)] * Xj;
	}
	return Result::Record({{"primal", Primal}, {"dual", Dual}, {"slack", Slack}, {"gap", Gap}});
}

Real SoftT(Real V, Real K) { return V > K ? V - K : (V < -K ? V + K : 0.0); }

Result Soft(const Call& C)
{
	const Real K = C.N(1);
	if (!(K >= 0.0))
	{
		return Result::Fail("threshold t must be >= 0");
	}
	Matrix M = C.M(0);
	for (Real& X : M.V)
	{
		X = SoftT(X, K);
	}
	return Result::Of(M);
}

// LASSO min 1/2 |A x - b|^2 + lambda |x|_1 by cyclic coordinate descent (each step is a soft threshold).
Result Lasso(const Call& C)
{
	const Matrix& A = C.M(0);
	Matrix B;
	const Real Lam = C.N(2);
	if (A.IsEmpty() || !AsColumn(C.M(1), A.Rows, &B))
	{
		return Result::Fail("b must have one entry per row of A");
	}
	if (!(Lam >= 0.0))
	{
		return Result::Fail("lambda must be >= 0");
	}
	if (A.Rows * A.Cols > 4000)
	{
		return Result::Fail("at most 4000 entries in A");
	}
	const int M = A.Rows, N = A.Cols;
	std::vector<Real> X(static_cast<size_t>(N), 0.0), Res(B.V), Col2(static_cast<size_t>(N), 0.0);
	for (int j = 0; j < N; ++j)
	{
		for (int i = 0; i < M; ++i)
		{
			Col2[static_cast<size_t>(j)] += A(i, j) * A(i, j);
		}
		if (!(Col2[static_cast<size_t>(j)] > 0.0))
		{
			return Result::Fail("column " + std::to_string(j + 1) + " of A is zero");
		}
	}
	for (int Sweep = 0; Sweep < 20000; ++Sweep)
	{
		Real Change = 0, Scale = 0;
		for (int j = 0; j < N; ++j)
		{
			Real Rho = 0;
			const size_t J = static_cast<size_t>(j);
			for (int i = 0; i < M; ++i)
			{
				Rho += A(i, j) * Res[static_cast<size_t>(i)];
			}
			Rho += Col2[J] * X[J];
			const Real Xn = SoftT(Rho, Lam) / Col2[J];
			const Real D = Xn - X[J];
			if (D != 0.0)
			{
				for (int i = 0; i < M; ++i)
				{
					Res[static_cast<size_t>(i)] -= A(i, j) * D;
				}
			}
			X[J] = Xn;
			Change = gm::max(Change, gm::fabs(D));
			Scale = gm::max(Scale, gm::fabs(Xn));
		}
		if (Change <= 1e-15 * gm::max(1.0, Scale))
		{
			break;
		}
	}
	return Result::Of(Matrix::Column(X));
}

// Backward Riccati P_k, k = 0..N, for x+ = A x + B u with cost 1/2 sum(x'Qx + u'Ru) + 1/2 x_N'Qf x_N.
bool LqRiccati(const Call& C, int N, std::vector<Matrix>* Ps, std::vector<Matrix>* Ks, std::string* Err)
{
	const Matrix &A = C.M(0), &B = C.M(1), &Q = C.M(2), &R = C.M(3), &Qf = C.M(4);
	const int Nx = A.Rows, Nu = B.Cols;
	if (!A.IsSquare() || B.Rows != Nx || Q.Rows != Nx || Q.Cols != Nx || R.Rows != Nu || R.Cols != Nu || Qf.Rows != Nx || Qf.Cols != Nx)
	{
		*Err = "shapes: A n x n, B n x m, Q n x n, R m x m, Qf n x n";
		return false;
	}
	if (N < 1 || N > 2000)
	{
		*Err = "horizon N must be 1..2000";
		return false;
	}
	Ps->assign(static_cast<size_t>(N + 1), Matrix());
	Ks->assign(static_cast<size_t>(N), Matrix());
	(*Ps)[static_cast<size_t>(N)] = Qf;
	const Matrix At = la::Transpose(A), Bt = la::Transpose(B);
	for (int k = N - 1; k >= 0; --k)
	{
		const Matrix& P1 = (*Ps)[static_cast<size_t>(k + 1)];
		Matrix K;
		if (!la::Solve(la::Add(R, la::Mul(la::Mul(Bt, P1), B)), la::Mul(la::Mul(Bt, P1), A), &K, Err))
		{
			return false;
		}
		(*Ks)[static_cast<size_t>(k)] = K;
		(*Ps)[static_cast<size_t>(k)] = la::Symmetrize(la::Sub(la::Add(Q, la::Mul(la::Mul(At, P1), A)), la::Mul(la::Mul(la::Mul(At, P1), B), K)));
	}
	return true;
}

Result LqCostate(const Call& C)
{
	const int N = C.Int(6), K = C.Int(7);
	std::vector<Matrix> Ps, Ks;
	std::string E;
	if (!LqRiccati(C, N, &Ps, &Ks, &E))
	{
		return Result::Fail(E);
	}
	Matrix X;
	if (!AsColumn(C.M(5), C.M(0).Rows, &X))
	{
		return Result::Fail("x0 must have n entries");
	}
	if (K < 0 || K > N)
	{
		return Result::Fail("step k must be 0..N");
	}
	for (int k = 0; k < K; ++k)
	{
		const Matrix U = la::Scale(la::Mul(Ks[static_cast<size_t>(k)], X), -1.0);
		X = la::Add(la::Mul(C.M(0), X), la::Mul(C.M(1), U));
	}
	return Result::Of(la::Mul(Ps[static_cast<size_t>(K)], X));
}

Result LqCost(const Call& C)
{
	const int N = C.Int(6);
	std::vector<Matrix> Ps, Ks;
	std::string E;
	if (!LqRiccati(C, N, &Ps, &Ks, &E))
	{
		return Result::Fail(E);
	}
	Matrix X;
	if (!AsColumn(C.M(5), C.M(0).Rows, &X))
	{
		return Result::Fail("x0 must have n entries");
	}
	return Result::Of(0.5 * la::Mul(la::Mul(la::Transpose(X), Ps[0]), X).ScalarValue());
}

// Real roots of c3 x^3 + c2 x^2 + c1 x + c0 (c3 may be zero), polished by Newton.
std::vector<Real> RealRootsCubic(Real C3, Real C2, Real C1, Real C0)
{
	std::vector<Real> Out;
	if (gm::fabs(C3) < 1e-300)
	{
		if (gm::fabs(C2) < 1e-300)
		{
			if (gm::fabs(C1) > 0)
			{
				Out.push_back(-C0 / C1);
			}
			return Out;
		}
		const Real Disc = C1 * C1 - 4 * C2 * C0;
		if (Disc >= 0)
		{
			Out.push_back((-C1 + gm::sqrt(Disc)) / (2 * C2));
			Out.push_back((-C1 - gm::sqrt(Disc)) / (2 * C2));
		}
		return Out;
	}
	Matrix Comp(3, 3);
	Comp(0, 0) = -C2 / C3;
	Comp(0, 1) = -C1 / C3;
	Comp(0, 2) = -C0 / C3;
	Comp(1, 0) = 1;
	Comp(2, 1) = 1;
	std::vector<Real> Re, Im;
	std::string E;
	if (!la::Eigenvalues(Comp, &Re, &Im, &E))
	{
		return Out;
	}
	for (size_t i = 0; i < Re.size(); ++i)
	{
		if (gm::fabs(Im[i]) < 1e-7 * gm::max(1.0, gm::fabs(Re[i])))
		{
			Real X = Re[i];
			for (int it = 0; it < 6; ++it)
			{
				const Real F = ((C3 * X + C2) * X + C1) * X + C0, D = (3 * C3 * X + 2 * C2) * X + C1;
				if (D == 0)
				{
					break;
				}
				X -= F / D;
			}
			Out.push_back(X);
		}
	}
	return Out;
}

Real GapF(Real A, Real X) { return (1 - A) * (X - 1) * (X - 1) + A * ((X * X - 1) * (X * X - 1) - X / 2); }

// g(lambda) = min_x f_a(x) + lambda x over all real x.
Real GapDualFn(Real A, Real Lam)
{
	// f' + lambda = 4a x^3 + (2(1-a) - 4a) x - 2(1-a) - a/2 + lambda
	const std::vector<Real> Rts = RealRootsCubic(4 * A, 0.0, 2 * (1 - A) - 4 * A, -2 * (1 - A) - A / 2 + Lam);
	Real G = Real::Infinity();
	for (Real X : Rts)
	{
		G = gm::min(G, GapF(A, X) + Lam * X);
	}
	return G;
}

// Ch06 counterexample: min f_a(x) s.t. x <= 0, with Slater holding; the gap opens with the Real well.
Result DualGap(const Call& C)
{
	const Real A = C.N(0);
	if (!(A >= 0.0 && A <= 1.0))
	{
		return Result::Fail("nonconvexity a must be in [0, 1]");
	}
	Real Pstar = GapF(A, 0.0);
	for (Real X : RealRootsCubic(4 * A, 0.0, 2 * (1 - A) - 4 * A, -2 * (1 - A) - A / 2))
	{
		if (X <= 0.0)
		{
			Pstar = gm::min(Pstar, GapF(A, X));
		}
	}
	// g is concave: scan a grid for the bracket, then golden-section search.
	const int Ng = 400;
	const Real LamMax = 20.0;
	int Best = 0;
	Real Gbest = -Real::Infinity();
	for (int i = 0; i <= Ng; ++i)
	{
		const Real G = GapDualFn(A, LamMax * i / Ng);
		if (G > Gbest)
		{
			Gbest = G;
			Best = i;
		}
	}
	Real Lo = LamMax * std::max(0, Best - 1) / Ng, Hi = LamMax * std::min(Ng, Best + 1) / Ng;
	const Real Gr = 0.5 * (gm::sqrt(5.0) - 1.0);
	for (int it = 0; it < 200; ++it)
	{
		const Real L1 = Hi - Gr * (Hi - Lo), L2 = Lo + Gr * (Hi - Lo);
		if (GapDualFn(A, L1) < GapDualFn(A, L2))
		{
			Lo = L1;
		}
		else
		{
			Hi = L2;
		}
	}
	const Real Dstar = gm::max(Gbest, GapDualFn(A, 0.5 * (Lo + Hi)));
	return Result::Record({{"p", Pstar}, {"d", Dstar}, {"gap", Pstar - Dstar}, {"lambda", 0.5 * (Lo + Hi)}});
}

// Bertsekas' auction, Gauss-Seidel bidding from a FIFO queue of unassigned bidders (ch06 bench 3).
Result Auction(const Call& C)
{
	const Matrix& A = C.M(0);
	const Real Eps = C.N(1);
	if (!A.IsSquare() || A.Rows > 50)
	{
		return Result::Fail("benefits must be a square n x n matrix, n <= 50");
	}
	if (!(Eps > 0.0))
	{
		return Result::Fail("bid increment eps must be positive");
	}
	const int N = A.Rows;
	std::vector<Real> Price(static_cast<size_t>(N), 0.0);
	std::vector<int> Owner(static_cast<size_t>(N), -1), Assign(static_cast<size_t>(N), -1), Queue;
	for (int i = 0; i < N; ++i)
	{
		Queue.push_back(i);
	}
	size_t Head = 0;
	long long Bids = 0;
	while (Head < Queue.size())
	{
		if (++Bids > 2000000)
		{
			return Result::Fail("too many bids; raise eps");
		}
		const int I = Queue[Head++];
		int J = 0;
		for (int j = 1; j < N; ++j)
		{
			if (A(I, j) - Price[static_cast<size_t>(j)] > A(I, J) - Price[static_cast<size_t>(J)])
			{
				J = j;
			}
		}
		const Real V1 = A(I, J) - Price[static_cast<size_t>(J)];
		Real V2 = -Real::Infinity();
		for (int j = 0; j < N; ++j)
		{
			if (j != J)
			{
				V2 = gm::max(V2, A(I, j) - Price[static_cast<size_t>(j)]);
			}
		}
		if (N == 1)
		{
			V2 = V1;
		}
		Price[static_cast<size_t>(J)] += V1 - V2 + Eps;
		const int Prev = Owner[static_cast<size_t>(J)];
		if (Prev >= 0)
		{
			Assign[static_cast<size_t>(Prev)] = -1;
			Queue.push_back(Prev);
		}
		Owner[static_cast<size_t>(J)] = I;
		Assign[static_cast<size_t>(I)] = J;
	}
	Matrix Out(1, N);
	for (int i = 0; i < N; ++i)
	{
		Out.V[static_cast<size_t>(i)] = Assign[static_cast<size_t>(i)] + 1.0;
	}
	return Result::Of(Out);
}
} // namespace atlas_detail

void RegisterAtlasLibrary(Registry& R)
{
	using namespace atlas_detail;
	const Real D2R = cst::Pi() / 180.0;

	// ch01 groups
	const char* G = "groups";
	Def(R, "hat", G, "ch01 groups, definitions: hat map so(3)", "skew matrix [w]x with [w]x v = w x v", "3 x 3", {Mat("w", "-")}, &Hat);
	Def(R, "vee", G, "ch01 groups, definitions: hat map so(3)", "inverse of the hat map (of the skew part)", "1 x 3", {Mat("S", "-")}, &Vee);
	Def(R, "rodrigues", G, "ch01 groups, Theorem 1 (Euler, Rodrigues)", "exp(theta [n]x) = I + sin(theta)[n]x + (1 - cos theta)[n]x^2; n is normalized", "3 x 3 rotation", {Mat("n", "-"), Num("theta", "rad")}, &Rodrigues);
	Def(R, "exp_so3", G, "ch01 groups, definitions: exponential map", "rotation from a rotation vector w = theta n", "3 x 3 rotation", {Mat("w", "rad")}, &ExpSo3);
	Def(R, "log_so3", G, "ch01 groups, Theorem 1: matrix logarithm", "rotation vector theta n of R, theta in [0, pi]", "1 x 3, rad", {Mat("R", "-")}, &LogSo3Fn);
	Def(R, "rot_angle", G, "ch01 groups, bench 1 (angle of R)", "rotation angle of R, always in [0, pi]", "rad", {Mat("R", "-")}, &RotAngle);
	Def(R, "quat", G, "ch01 groups, Theorem 2 (Real cover)", "axis-angle to unit quaternion (cos(theta/2), sin(theta/2) n)", "1 x 4 [w, x, y, z]", {Mat("n", "-"), Num("theta", "rad")}, &QuatFn);
	Def(R, "quat_mul", G, "ch01 groups, definitions: Hamilton product", "Hamilton product p (x) q, scalar first", "1 x 4", {Mat("p", "-"), Mat("q", "-")}, &QuatMul);
	Def(R, "quat_conj", G, "ch01 groups, definitions: unit quaternions", "conjugate (w, -v), the inverse of a unit quaternion", "1 x 4", {Mat("q", "-")}, &QuatConj);
	Def(R, "quat_rot", G, "ch01 groups, Theorem 2 proof: R(q) v = q (0,v) q^-1", "rotate a vector by a quaternion (q is normalized)", "1 x 3", {Mat("q", "-"), Mat("v", "")}, &QuatRot);
	Def(R, "quat2dcm", G, "ch01 groups, definitions: R(q)", "R(q) = (w^2 - v.v)I + 2vv' + 2w[v]x (q is normalized)", "3 x 3 rotation", {Mat("q", "-")}, &Quat2Dcm);
	Def(R, "dcm2quat", G, "ch01 groups, Theorem 2 (Real cover)", "unit quaternion of a rotation matrix (Shepperd), sign chosen with w >= 0", "1 x 4", {Mat("R", "-")}, &Dcm2Quat);
	Def(R, "slerp", G, "ch01 groups, other guises: slerp (Shoemake 1985)", "geodesic interpolation on S^3, the short way round", "1 x 4 unit quaternion", {Mat("q0", "-"), Mat("q1", "-"), Num("t", "-", 0.5)}, &Slerp);
	Def(R, "bracket", G, "ch01 groups, definitions: Lie bracket of so(3)", "[[a]x, [b]x] as a vector (= a x b), or A B - B A for square matrices", "1 x 3 or n x n", {Mat("a", "-"), Mat("b", "-")}, &Bracket);
	Def(R, "comm_angle", G, "ch01 groups, counterexample 1 (order matters)", "angle of Rx(a) Ry(b) Rx(a)' Ry(b)', ~ a b for small angles", "rad", {Num("alpha", "rad", 60 * D2R), Num("beta", "rad", 60 * D2R)}, &CommAngle);
	Def(R, "davenport", G, "ch01 groups, guises: Wahba, Davenport q-method", "Davenport's K = [sigma, z'; z, S - sigma I] (Hamilton); r, b are n x 3, one vector per row", "4 x 4",
		{Mat("r", "-"), Mat("b", "-"), MatOr("w", "-", 1.0)}, &Davenport);
	Def(R, "qmethod", G, "ch01 groups, guises: Davenport q-method", "optimal Wahba quaternion: top eigenvector of K (w >= 0)", "1 x 4", {Mat("r", "-"), Mat("b", "-"), MatOr("w", "-", 1.0)}, &QMethod);
	Def(R, "wahba_svd", G, "ch01 groups, guises: SVD solution (Markley 1988)", "A = U diag(1, 1, det U det V) V' for B = sum w b r'", "3 x 3 rotation, b ~ A r", {Mat("r", "-"), Mat("b", "-"), MatOr("w", "-", 1.0)}, &WahbaSvd);
	Def(R, "wahba_loss", G, "ch01 groups, guises: Wahba's problem", "1/2 sum w |b - A r|^2", "-", {Mat("A", "-"), Mat("r", "-"), Mat("b", "-"), MatOr("w", "-", 1.0)}, &WahbaLoss);
	Def(R, "triad", G, "ch01 groups, guises: TRIAD (Black 1964)", "rotation from two vector pairs, trusting r1/b1 fully", "3 x 3 rotation, b ~ A r", {Mat("r1", "-"), Mat("r2", "-"), Mat("b1", "-"), Mat("b2", "-")}, &Triad);
	Def(R, "quat_mean", G, "ch01 groups, counterexample 2 (Markley et al. 2007)", "attitude average: top eigenvector of sum w q q' (sign-insensitive); q is n x 4", "1 x 4", {Mat("q", "-"), MatOr("w", "-", 1.0)}, &QuatMean);
	Def(R, "exp_se3", G, "ch01 groups, definitions: rigid motions SE(3)", "exp of the twist (omega, nu): [exp([w]x), V nu; 0, 1]", "4 x 4", {Mat("omega", "rad"), Mat("nu", "m")}, &ExpSe3);
	Def(R, "att_err", G, "ch01 groups, guises: geometric control on SO(3) (Lee et al. 2010)", "configuration error e_R = 1/2 (Rd'R - R'Rd)^vee", "1 x 3", {Mat("R", "-"), Mat("Rd", "-")}, &AttErr);

	// ch02 topology
	const char* T = "topology";
	Def(R, "euler_rates", T, "ch02 topology, guises: gimbal lock in Euler angles", "3-2-1 Euler-angle rates E(phi, theta) [p, q, r]", "record dphi, dtheta, dpsi (rad/s)",
		{Num("phi", "rad", 20 * D2R), Num("theta", "rad", 60 * D2R), Num("p", "rad/s", 0), Num("q", "rad/s", 0), Num("r", "rad/s", 0.2)}, &EulerRates);
	Def(R, "euler_det", T, "ch02 topology, guises: det E = 1/cos(theta)", "determinant of the 3-2-1 kinematics matrix E", "-", {Num("theta", "rad", 60 * D2R), Num("phi", "rad", 0)}, &EulerDet);
	Def(R, "euler_cond", T, "ch02 topology, counterexample 1 (gimbal lock, measured)", "2-norm condition number of E(phi, theta)", "-", {Num("phi", "rad", 20 * D2R), Num("theta", "rad", 60 * D2R)}, &EulerCond);
	Def(R, "eul2dcm", T, "ch02 topology, guises: 3-2-1 (yaw-pitch-roll) sequence", "body-to-world R = Rz(psi) Ry(theta) Rx(phi)", "3 x 3 rotation", {Num("phi", "rad"), Num("theta", "rad"), Num("psi", "rad")}, &Eul2Dcm);
	Def(R, "dcm2eul", T, "ch02 topology, guises: 3-2-1 sequence", "3-2-1 Euler angles of a body-to-world R (theta in [-pi/2, pi/2])", "record phi, theta, psi (rad)", {Mat("R", "-")}, &Dcm2Eul);
	Def(R, "lookat", T, "ch02 topology, guises: pointing and look-at frames", "up-vector u(d) = (z - (z.d)d)/|...|; fails at d = +-z (hairy ball)", "1 x 3", {Mat("d", "-")}, &LookAt);
	Def(R, "navfn", T, "ch02 topology, guises: navigation functions (Rimon-Koditschek)", "phi = gamma/(gamma^k + beta)^(1/k) on the chapter's disk with 3 obstacles", "record phi, fx, fy",
		{Num("x", "-"), Num("y", "-"), Num("k", "-", 2), Num("gx", "-", 0), Num("gy", "-", 0.6), Num("r", "-", 0.25)}, &NavFn);
	Def(R, "nav_index", T, "ch02 topology, bench 4 (Poincare-Hopf index count)", "gradient winding per grid cell (n <= 400); total = chi = 1 - M = -2", "record total, plus, minus (cells)",
		{Num("k", "-", 2), Num("gx", "-", 0), Num("gy", "-", 0.6), Num("r", "-", 0.25), Num("n", "-", 220)}, &NavIndex);
	Def(R, "unwind", T, "ch02 topology, counterexample 2 (unwinding)", "25 s quaternion feedback sim; law 0 naive, 1 sign switch, 2 hybrid (mulberry32 seed 99)", "record travelled (rad), switches, final (rad)",
		{Num("theta0", "rad", 340 * D2R), Num("law", "-", 0), Num("sigma", "-", 0.05), Num("w0", "rad/s", 0)}, &Unwind);

	// ch04 Hilbert spaces
	const char* H = "hilbert";
	Def(R, "proj", H, "ch04 Hilbert, Theorem 1 (projection theorem)", "orthogonal projection of y onto the column space of A", "vector like y", {Mat("A", ""), Mat("y", "")}, &Proj);
	Def(R, "gs", H, "ch04 Hilbert, Theorem 3 proof (innovations are Gram-Schmidt)", "orthonormalize the columns of A (modified Gram-Schmidt)", "m x n", {Mat("A", "")}, &GramSchmidt);
	Def(R, "wls", H, "ch04 Hilbert, guises: weighted least squares", "x = (A'R^-1 A)^-1 A'R^-1 y; R scalar, variances or covariance", "column", {Mat("A", ""), Mat("y", ""), MatOr("R", "", 1.0)}, &Wls);
	Def(R, "vander", H, "ch04 Hilbert, bench 1 (monomial basis)", "design matrix [1, x, x^2, ..., x^d] (cond(A'A) = cond(A)^2)", "m x (d+1)", {Mat("x", ""), Num("d", "-")}, &Vander);
	Def(R, "kf_rw", H, "ch04 Hilbert, bench 2 (Kalman = batch)", "Kalman filter for a random walk, prior N(0, P0)", "1 x n filtered estimates", {Mat("y", ""), Num("q", "-", 0.2), Num("r", "-", 1), Num("P0", "-", 4)}, &KfRw);
	Def(R, "smooth_rw", H, "ch04 Hilbert, bench 2 / guises: smoothing", "batch least squares over the whole record (= RTS smoother)", "1 x n", {Mat("y", ""), Num("q", "-", 0.2), Num("r", "-", 1), Num("P0", "-", 4)}, &SmoothRw);
	Def(R, "kf_ss", H, "ch04 Hilbert, bench 4 (Wiener vs Kalman)", "steady-state scalar Kalman filter for x+ = a x + w, y = x + v", "record P (prior), K, Pp (posterior)", {Num("a", "-", 0.95), Num("q", "-", 0.1), Num("r", "-", 1)}, &KfSs);
	Def(R, "wiener", H, "ch04 Hilbert, guises: Wiener filter", "non-causal Wiener MSE (closed form) vs causal (steady Kalman) MSE", "record nc, causal, price", {Num("a", "-", 0.95), Num("q", "-", 0.1), Num("r", "-", 1)}, &Wiener);
	Def(R, "bode_int", H, "ch04 Hilbert, Theorem 5 / counterexample 2 (waterbed)", "int_0^inf ln|S(jw)| dw for L = k(s+1)/((s-p)(s+3)(s+4)) vs pi*max(p,0)", "record integral, theorem, peak, wpeak (rad/s), stable",
		{Num("k", "-", 30), Num("p", "1/s", 0.5)}, &BodeInt);
	Def(R, "postmean", H, "ch04 Hilbert, counterexample 1 (best linear is not best)", "posterior mean of x ~ N(0,4) under N(0,1)/N(0,100) mixture noise vs best linear", "record post, lin, linmse",
		{Mat("y", ""), Num("eps", "-", 0.1)}, &PostMean);

	// ch05 measure
	const char* M = "measure";
	Def(R, "bayes", M, "ch05 measure, guises: Kalman as the Gaussian Bayes filter", "scalar Gaussian update of N(m, P) by y = H x + N(0, R)", "record m, P, K", {Num("m", ""), Num("P", ""), Num("y", ""), Num("R", ""), Num("H", "-", 1)}, &Bayes);
	Def(R, "bayes_mean", M, "ch05 measure, guises: Kalman as the Gaussian Bayes filter", "posterior mean x + K (y - H x), K = P H'(H P H' + R)^-1", "column", {Mat("x", ""), Mat("P", ""), Mat("y", ""), Mat("H", ""), Mat("R", "")}, &BayesMean);
	Def(R, "bayes_cov", M, "ch05 measure, guises: Kalman as the Gaussian Bayes filter", "posterior covariance P - K S K'", "n x n", {Mat("P", ""), Mat("H", ""), Mat("R", "")}, &BayesCov);
	Def(R, "bk_density", M, "ch05 measure, counterexample 1 (Borel-Kolmogorov)", "conditional density on a great circle: band = 1/(2 pi) in longitude t, wedge = cos(t)/2 in latitude t", "record band, wedge (1/rad)", {Num("t", "rad")}, &BkDensity);
	Def(R, "is_w2", M, "ch05 measure, Theorem 4 / counterexample 2", "importance sampling N(0,1) from N(0,s^2): E_q[w^2] = s^2/sqrt(2s^2-1), inf for s^2 <= 1/2", "record w2, ess (fraction), fvar (per sample, f = x^2)",
		{Num("s", "-", 0.65)}, &IsW2);
	Def(R, "vdc", M, "ch05 measure, bench 4 (radical inverse)", "van der Corput radical inverse of i in base b", "[0, 1)", {Num("i", "-"), Num("b", "-", 2)}, &Vdc);
	Def(R, "halton", M, "ch05 measure, Theorem 5 / bench 4 (Halton 1960)", "first n Halton points in d dimensions (indices 1..n, prime bases; n <= 20000, n d <= 80000)", "n x d", {Num("n", "-"), Num("d", "-", 2)}, &Halton);
	Def(R, "qmc_err", M, "ch05 measure, bench 4 (random vs well-placed points)", "Halton integration error: kind 1 prod (pi/2) sin(pi x), kind 0 indicator of sum < d/2 (N d <= 100000)", "absolute error",
		{Num("N", "-", 8192), Num("d", "-", 4), Num("kind", "-", 1)}, &QmcErr);
	Def(R, "ess", M, "ch05 measure, guises: particle filter N_eff", "effective sample size (sum w)^2 / sum w^2", "particles", {Mat("w", "-")}, &Ess);
	Def(R, "logz", M, "ch05 measure, Gibbs variational principle", "log partition function log sum exp(-E)", "-", {Mat("E", "-")}, &LogZFn);
	Def(R, "boltz", M, "ch05 measure, Gibbs variational principle", "Boltzmann distribution exp(-E)/Z", "same shape as E", {Mat("E", "-")}, &Boltz);
	Def(R, "gibbs_f", M, "ch05 measure, Gibbs variational principle", "-E_q[E] + H(q) (q normalized); its maximum is logz(E)", "nats", {Mat("q", "-"), Mat("E", "-")}, &GibbsF);
	Def(R, "kl", M, "ch05 measure, Gibbs variational principle", "Kullback-Leibler divergence sum p log(p/q) (p, q normalized)", "nats", {Mat("p", "-"), Mat("q", "-")}, &Kl);

	// ch06 optimization
	const char* O = "optimize";
	Def(R, "lp", O, "ch06 optimization, Theorems 1-3 (simplex)", "optimal x of min c'x s.t. A x <= b, x >= 0", "column", {Mat("c", ""), Mat("A", ""), Mat("b", "")}, &Lp);
	Def(R, "lpval", O, "ch06 optimization, Theorems 1-3", "optimal value of min c'x s.t. A x <= b, x >= 0", "c'x*", {Mat("c", ""), Mat("A", ""), Mat("b", "")}, &LpVal);
	Def(R, "lpdual", O, "ch06 optimization, Theorem 2 (strong duality)", "dual prices y >= 0: c + A'y >= 0, c'x* = -b'y*, y = -d(c'x*)/db", "column", {Mat("c", ""), Mat("A", ""), Mat("b", "")}, &LpDual);
	Def(R, "lp_kkt", O, "ch06 optimization, definitions: KKT conditions", "KKT residuals of (x, y) for the LP: feasibility, slackness, gap c'x + b'y", "record primal, dual, slack, gap",
		{Mat("c", ""), Mat("A", ""), Mat("b", ""), Mat("x", ""), Mat("y", "")}, &LpKkt);
	Def(R, "soft", O, "ch06 optimization, LASSO proximal map", "soft threshold sign(v) max(|v| - t, 0), the prox of t|.|_1", "element-wise", {Mat("v", ""), Num("t", "")}, &Soft);
	Def(R, "lasso", O, "ch06 optimization, proximal fixed point (LASSO)", "argmin 1/2 |A x - b|^2 + lambda |x|_1 (coordinate descent)", "column", {Mat("A", ""), Mat("b", ""), Num("lambda", "", 1)}, &Lasso);
	Def(R, "lq_costate", O, "ch06 optimization, Theorem 5 (costate = P_k x_k)", "KKT multiplier of step k's dynamics: lambda_k = P_k x_k on the optimal path", "column",
		{Mat("A", ""), Mat("B", ""), Mat("Q", ""), Mat("R", ""), Mat("Qf", ""), Mat("x0", ""), Num("N", "-"), Num("k", "-", 0)}, &LqCostate);
	Def(R, "lq_cost", O, "ch06 optimization, Theorem 5 (cost-to-go)", "optimal finite-horizon LQ cost 1/2 x0' P_0 x0", "-", {Mat("A", ""), Mat("B", ""), Mat("Q", ""), Mat("R", ""), Mat("Qf", ""), Mat("x0", ""), Num("N", "-")}, &LqCost);
	Def(R, "dualgap", O, "ch06 optimization, counterexample 1 (where duality breaks)", "min (1-a)(x-1)^2 + a((x^2-1)^2 - x/2) s.t. x <= 0: primal, dual and gap", "record p, d, gap, lambda", {Num("a", "-", 0.6)}, &DualGap);
	Def(R, "auction", O, "ch06 optimization, Theorem 6 / bench 3 (auction)", "assignment maximizing sum a(i, j_i) within n eps (Bertsekas)", "1 x n column chosen by each row", {Mat("a", ""), Num("eps", "", 1)}, &Auction);
}
} // namespace mcalc
