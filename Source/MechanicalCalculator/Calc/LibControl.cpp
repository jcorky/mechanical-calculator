// Control: the Riccati equation as a Mobius map (atlas demo and chapter 3), LQR, Kalman and MPC with the
// approach and MPC benches, and dynamical systems with Lyapunov certificates (chapter 8).
#include "CalcLibrary.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace mcalc
{
namespace control_detail
{
constexpr double CtlInf = std::numeric_limits<double>::infinity();
constexpr double CtlNan = std::numeric_limits<double>::quiet_NaN();

/* ---------- shapes and small linear-algebra helpers ---------- */

std::string ShapeText(const Matrix& M)
{
	return std::to_string(M.Rows) + " x " + std::to_string(M.Cols);
}

// (A, B, Q, R) of a regulator: A n x n, B n x m, Q n x n, R m x m.
bool CheckRegulator(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, std::string* Err)
{
	if (!A.IsSquare())
	{
		*Err = "A must be square, got " + ShapeText(A);
		return false;
	}
	const int N = A.Rows;
	if (B.Rows != N || B.Cols < 1)
	{
		*Err = "B must have " + std::to_string(N) + " rows (A is " + ShapeText(A) + "), got " + ShapeText(B);
		return false;
	}
	if (Q.Rows != N || Q.Cols != N)
	{
		*Err = "Q must be " + std::to_string(N) + " x " + std::to_string(N) + ", got " + ShapeText(Q);
		return false;
	}
	if (R.Rows != B.Cols || R.Cols != B.Cols)
	{
		*Err = "R must be " + std::to_string(B.Cols) + " x " + std::to_string(B.Cols) + " (one row per input), got " + ShapeText(R);
		return false;
	}
	return true;
}

// (A, C, Q, R) of a filter: A n x n, C p x n, Q n x n, R p x p.
bool CheckFilter(const Matrix& A, const Matrix& Cm, const Matrix& Q, const Matrix& R, std::string* Err)
{
	if (!A.IsSquare())
	{
		*Err = "A must be square, got " + ShapeText(A);
		return false;
	}
	const int N = A.Rows;
	if (Cm.Cols != N || Cm.Rows < 1)
	{
		*Err = "C must have " + std::to_string(N) + " columns (A is " + ShapeText(A) + "), got " + ShapeText(Cm);
		return false;
	}
	if (Q.Rows != N || Q.Cols != N)
	{
		*Err = "Q must be " + std::to_string(N) + " x " + std::to_string(N) + ", got " + ShapeText(Q);
		return false;
	}
	if (R.Rows != Cm.Rows || R.Cols != Cm.Rows)
	{
		*Err = "R must be " + std::to_string(Cm.Rows) + " x " + std::to_string(Cm.Rows) + " (one row per measurement), got " + ShapeText(R);
		return false;
	}
	return true;
}

bool IsWhole(double X, double Lo, double Hi)
{
	return std::isfinite(X) && X == std::floor(X) && X >= Lo && X <= Hi;
}

// A vector of length n given as a row or a column, returned as a column.
bool AsColumn(const Matrix& V, int N, const char* Name, Matrix* Out, std::string* Err)
{
	if (!V.IsVector() || V.Count() != N)
	{
		*Err = std::string(Name) + " must be a vector of length " + std::to_string(N) + ", got " + ShapeText(V);
		return false;
	}
	*Out = Matrix::Column(V.V);
	return true;
}

Matrix Block(const Matrix& A11, const Matrix& A12, const Matrix& A21, const Matrix& A22)
{
	const int R1 = A11.Rows, R2 = A21.Rows, C1 = A11.Cols, C2 = A12.Cols;
	Matrix Out(R1 + R2, C1 + C2);
	for (int i = 0; i < R1; ++i)
	{
		for (int j = 0; j < C1; ++j)
		{
			Out(i, j) = A11(i, j);
		}
		for (int j = 0; j < C2; ++j)
		{
			Out(i, C1 + j) = A12(i, j);
		}
	}
	for (int i = 0; i < R2; ++i)
	{
		for (int j = 0; j < C1; ++j)
		{
			Out(R1 + i, j) = A21(i, j);
		}
		for (int j = 0; j < C2; ++j)
		{
			Out(R1 + i, C1 + j) = A22(i, j);
		}
	}
	return Out;
}

Matrix SubBlock(const Matrix& A, int R0, int C0, int Rows, int Cols)
{
	Matrix Out(Rows, Cols);
	for (int i = 0; i < Rows; ++i)
	{
		for (int j = 0; j < Cols; ++j)
		{
			Out(i, j) = A(R0 + i, C0 + j);
		}
	}
	return Out;
}

// Largest real part of the eigenvalues; NaN if the iteration fails.
double MaxRealPart(const Matrix& A)
{
	std::vector<double> Re, Im;
	std::string E;
	if (!la::Eigenvalues(A, &Re, &Im, &E) || Re.empty())
	{
		return CtlNan;
	}
	return *std::max_element(Re.begin(), Re.end());
}

// One step of the regulator Riccati recursion P <- Q + A'PA - A'PB (R + B'PB)^-1 B'PA.
bool RiccatiStep(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, const Matrix& P, Matrix* Next, Matrix* Gain, std::string* Err)
{
	const Matrix Bt = la::Transpose(B);
	const Matrix BtP = la::Mul(Bt, P);
	const Matrix S = la::Add(R, la::Mul(BtP, B));
	Matrix K;
	if (!la::Solve(S, la::Mul(BtP, A), &K, Err))
	{
		*Err = "R + B'PB is singular";
		return false;
	}
	if (Next)
	{
		const Matrix At = la::Transpose(A);
		const Matrix AtPA = la::Mul(la::Mul(At, P), A);
		*Next = la::Symmetrize(la::Sub(la::Add(Q, AtPA), la::Mul(la::Transpose(la::Mul(BtP, A)), K)));
	}
	if (Gain)
	{
		*Gain = K;
	}
	return true;
}

// Steady-state Kalman filter: prior covariance P from the dual DARE and the measurement-update gain
// M = P C' (C P C' + R)^-1 (approach bench section 2).
bool SteadyKalman(const Matrix& A, const Matrix& Cm, const Matrix& Q, const Matrix& R, Matrix* P, Matrix* M, std::string* Err)
{
	Matrix Pp, Kd;
	if (!la::Dare(la::Transpose(A), la::Transpose(Cm), Q, R, &Pp, &Kd, Err))
	{
		return false;
	}
	const Matrix PCt = la::Mul(Pp, la::Transpose(Cm));
	const Matrix S = la::Add(la::Mul(Cm, PCt), R);
	Matrix Mt;
	if (!la::Solve(la::Transpose(S), la::Transpose(PCt), &Mt, Err))
	{
		*Err = "C P C' + R is singular";
		return false;
	}
	*P = Pp;
	*M = la::Transpose(Mt);
	if (!(la::SpectralRadius(la::Mul(A, la::Sub(Matrix::Identity(A.Rows), la::Mul(*M, Cm)))) < 1))
	{
		*Err = "no stabilizing filter (is (A, C) detectable and (A, Q) stabilizable?)";
		return false;
	}
	return true;
}

// The connected interval of input gains around 1 on which Fn(g) < 1 (approach bench, margins()): march
// outward on a log grid of N/2 steps per decade-range, then bisect geometrically.
template <typename Fn>
double GainEdge(const Fn& F, double End, int N)
{
	const double K = std::log(End) / N * 2;
	double Gin = 1.0;
	for (int i = 1; i <= N; ++i)
	{
		const double G = std::exp(K * i);
		if ((End > 1 && G > End) || (End < 1 && G < End))
		{
			return End;
		}
		if (!(F(G) < 1))
		{
			double Gout = G;
			for (int j = 0; j < 50; ++j)
			{
				const double Mid = std::sqrt(Gin * Gout);
				if (F(Mid) < 1)
				{
					Gin = Mid;
				}
				else
				{
					Gout = Mid;
				}
			}
			return Gin;
		}
		Gin = G;
	}
	return End;
}

template <typename Fn>
bool GainInterval(const Fn& F, double Lo, double Hi, int N, double* OutLo, double* OutHi)
{
	if (!(F(1.0) < 1))
	{
		return false;
	}
	*OutLo = GainEdge(F, Lo, N);
	*OutHi = GainEdge(F, Hi, N);
	return true;
}

// LQG loop in (x, prior estimate) with a current-estimate law u = -K xhat+ and gain g at the plant input.
Matrix LqgLoop(const Matrix& A, const Matrix& B, const Matrix& Cm, const Matrix& K, const Matrix& M, double G)
{
	const int N = A.Rows;
	const Matrix I = Matrix::Identity(N);
	const Matrix MC = la::Mul(M, Cm);
	const Matrix KM = la::Mul(K, MC);
	const Matrix KI = la::Mul(K, la::Sub(I, MC));
	const Matrix BKM = la::Mul(B, KM);
	const Matrix BKI = la::Mul(B, KI);
	return Block(la::Sub(A, la::Scale(BKM, G)), la::Scale(BKI, -G), la::Sub(la::Mul(A, MC), BKM), la::Sub(la::Mul(A, la::Sub(I, MC)), BKI));
}

// sign(Z) by the scaled Newton iteration Z <- (cZ + Z^-1/c)/2.
bool MatrixSign(const Matrix& InZ, Matrix* Out, std::string* Err)
{
	Matrix Z = InZ;
	const int N2 = Z.Rows;
	for (int It = 0; It < 100; ++It)
	{
		Matrix Zi;
		if (!la::Inverse(Z, &Zi, Err))
		{
			*Err = "the Hamiltonian has eigenvalues on the imaginary axis";
			return false;
		}
		double D = 0.0;
		la::Determinant(Z, &D, nullptr);
		double Cs = std::pow(std::fabs(D), -1.0 / N2);
		if (!std::isfinite(Cs) || Cs <= 0.0)
		{
			Cs = 1.0;
		}
		const Matrix Zn = la::Scale(la::Add(la::Scale(Z, Cs), la::Scale(Zi, 1.0 / Cs)), 0.5);
		const double Change = la::MaxAbs(la::Sub(Zn, Z));
		Z = Zn;
		if (Change <= 1e-13 * std::max(1.0, la::MaxAbs(Z)))
		{
			*Out = Z;
			return true;
		}
	}
	*Err = "sign iteration did not converge (eigenvalues near the imaginary axis)";
	return false;
}

/* ---------- riccati: the scalar map, Mobius maps, Hamiltonian and symplectic matrices ---------- */

struct Mob
{
	double A, B, C, D;
};

// Scalar Kalman covariance update P -> a^2 r P/(h^2 P + r) + q as the Mobius map of M (atlas Phase 1).
Mob KalmanMob(double Av, double Q, double R, double H)
{
	return {Av * Av * R + Q * H * H, Q * R, H * H, R};
}

double MobApply(const Mob& M, double Z)
{
	if (std::isinf(Z))
	{
		return M.C != 0.0 ? M.A / M.C : CtlInf;
	}
	const double Den = M.C * Z + M.D;
	if (Den == 0.0)
	{
		return CtlInf;
	}
	return (M.A * Z + M.B) / Den;
}

struct FixedPoints
{
	double Plus = CtlNan;   // attracting, the larger eigenvalue C p + D
	double Minus = CtlNan;
	double Mult = CtlNan;   // f'(Plus) = det M / (C Plus + D)^2
};

// The fixed points of a map with C != 0 and real fixed points.
bool MobFixed(const Mob& M, FixedPoints* Out)
{
	const double Disc = (M.A - M.D) * (M.A - M.D) + 4 * M.B * M.C;
	if (M.C == 0.0 || Disc < 0.0)
	{
		return false;
	}
	const double Sq = std::sqrt(Disc);
	Out->Plus = ((M.A - M.D) + Sq) / (2 * M.C);
	Out->Minus = ((M.A - M.D) - Sq) / (2 * M.C);
	const double Den = M.C * Out->Plus + M.D;
	Out->Mult = (M.A * M.D - M.B * M.C) / (Den * Den);
	return true;
}

Result KfMap(const Call& C)
{
	return Result::Of(MobApply(KalmanMob(C.N(1), C.N(2), C.N(3), C.N(4)), C.N(0)));
}

Result KfIter(const Call& C)
{
	if (!IsWhole(C.N(0), 0, 1e6))
	{
		return Result::Fail("k must be a whole number 0..1e6");
	}
	const Mob M = KalmanMob(C.N(2), C.N(3), C.N(4), C.N(5));
	double P = C.N(1);
	const long long K = static_cast<long long>(C.N(0));
	for (long long i = 0; i < K; ++i)
	{
		P = MobApply(M, P);
	}
	return Result::Of(P);
}

// Steady state of the scalar Riccati map in either reading. The gain uses the reading's own formula.
Result ScalarSteady(double Av, double Q, double R, double H, bool Regulator)
{
	if (R <= 0.0)
	{
		return Result::Fail(Regulator ? "rho must be > 0" : "r must be > 0");
	}
	if (Q < 0.0)
	{
		return Result::Fail("q must be >= 0");
	}
	const Mob M = KalmanMob(Av, Q, R, H);
	FixedPoints Fp;
	double K = 0.0, Pole = Av, P = 0.0;
	if (H == 0.0)
	{
		if (std::fabs(Av) >= 1.0)
		{
			return Result::Fail(Regulator ? "b = 0 and |a| >= 1: the input cannot reach an unstable state, the cost diverges" : "h = 0 and |a| >= 1: nothing is observed, the variance diverges");
		}
		P = Q / (1.0 - Av * Av);
		Fp.Plus = P;
		Fp.Minus = CtlInf;
		Fp.Mult = Av * Av;
	}
	else
	{
		MobFixed(M, &Fp);
		P = Fp.Plus;
		if (Regulator)
		{
			K = Av * H * P / (R + H * H * P);
			Pole = Av - H * K;
		}
		else
		{
			K = P * H / (H * H * P + R);
			Pole = Av * (1 - K * H);
		}
	}
	if (Regulator)
	{
		return Result::Record({{"S", P}, {"Sm", Fp.Minus}, {"K", K}, {"m", Fp.Mult}, {"pole", Pole}});
	}
	return Result::Record({{"p", P}, {"pm", Fp.Minus}, {"K", K}, {"m", Fp.Mult}, {"pole", Pole}});
}

Result KfSs(const Call& C) { return ScalarSteady(C.N(0), C.N(1), C.N(2), C.N(3), false); }
Result LqrSs(const Call& C) { return ScalarSteady(C.N(0), C.N(1), C.N(2), C.N(3), true); }

Result KfW(const Call& C)
{
	const Mob M = KalmanMob(C.N(1), C.N(2), C.N(3), C.N(4));
	FixedPoints Fp;
	if (!MobFixed(M, &Fp))
	{
		return Result::Fail("h = 0: one fixed point at infinity, no cross-ratio");
	}
	const double Tr = M.A + M.D, Det = M.A * M.D - M.B * M.C;
	if (std::fabs(Tr * Tr - 4 * Det) < 1e-9 * (Tr * Tr + 1))
	{
		return Result::Fail("parabolic (|a| = 1, q = 0): no second fixed point; 1/P grows by h^2/r per step instead");
	}
	const double P = C.N(0);
	return Result::Of((P - Fp.Plus) / (P - Fp.Minus));
}

Result MobiusFn(const Call& C)
{
	const Mob M = {C.N(1), C.N(2), C.N(3), C.N(4)};
	return Result::Of(MobApply(M, C.N(0)));
}

Result MobIter(const Call& C)
{
	if (!IsWhole(C.N(0), 0, 1e6))
	{
		return Result::Fail("k must be a whole number 0..1e6");
	}
	const Mob M = {C.N(2), C.N(3), C.N(4), C.N(5)};
	double Z = C.N(1);
	const long long K = static_cast<long long>(C.N(0));
	for (long long i = 0; i < K; ++i)
	{
		Z = MobApply(M, Z);
	}
	return Result::Of(Z);
}

// Classification of the atlas demo (index.html analyze()): 0 degenerate, 1 identity, 2 parabolic,
// 3 hyperbolic, 4 elliptic; fixed points, the attracting one and its multiplier, the elliptic rotation.
Result MobClass(const Call& C)
{
	const double A = C.N(0), B = C.N(1), Cc = C.N(2), D = C.N(3);
	const double Tr = A + D, Det = A * D - B * Cc, Disc = Tr * Tr - 4 * Det;
	const double Scale = std::max(1e-9, Tr * Tr + std::fabs(Det));
	int Type = 3;
	if (std::fabs(Det) < 1e-12)
	{
		Type = 0;
	}
	else if (std::fabs(B) < 1e-12 && std::fabs(Cc) < 1e-12 && std::fabs(A - D) < 1e-12)
	{
		Type = 1;
	}
	else if (std::fabs(Disc) < 1e-9 * Scale)
	{
		Type = 2;
	}
	else
	{
		Type = Disc > 0 ? 3 : 4;
	}
	std::vector<double> Fps;
	if (Type >= 2)
	{
		if (std::fabs(Cc) < 1e-12)
		{
			Fps.push_back(CtlInf);
			if (std::fabs(D - A) > 1e-12)
			{
				Fps.push_back(B / (D - A));
			}
		}
		else if (Disc >= -1e-12 * Scale)
		{
			const double Sq = std::sqrt(std::max(0.0, Disc));
			Fps.push_back((A - D + Sq) / (2 * Cc));
			Fps.push_back((A - D - Sq) / (2 * Cc));
		}
	}
	auto Mult = [&](double Z) {
		if (std::isfinite(Z))
		{
			const double Den = Cc * Z + D;
			return Det / (Den * Den);
		}
		return std::fabs(Cc) < 1e-12 ? D / A : CtlNan;
	};
	double Att = CtlNan, AttMult = CtlNan;
	if (Type == 2 && !Fps.empty())
	{
		Att = Fps[0];
		AttMult = Mult(Att);
	}
	else if (Type == 3)
	{
		double Best = CtlInf;
		for (double Z : Fps)
		{
			const double Mz = std::fabs(Mult(Z));
			if (Mz < Best - 1e-12)
			{
				Best = Mz;
				Att = Z;
			}
		}
		if (Best < 1)
		{
			AttMult = Mult(Att);
		}
		else
		{
			Att = CtlNan;
		}
	}
	const double Rot = Type == 4 ? 2 * std::atan2(std::sqrt(-Disc), Tr) : 0.0;
	const double Z1 = Fps.empty() ? CtlNan : Fps[0];
	const double Z2 = Fps.size() > 1 ? Fps[1] : Z1;
	return Result::Record({{"type", static_cast<double>(Type)}, {"tr", Tr}, {"det", Det}, {"disc", Disc}, {"z1", Z1}, {"z2", Z2}, {"att", Att}, {"m", AttMult}, {"rot", Rot}});
}

// Escape time of p' = alpha + 2 beta p + gamma p^2 from p0. With u = gamma p + beta the flow is
// u' = u^2 + (alpha gamma - beta^2), so it escapes in finite time exactly when u passes through infinity.
Result RiccEsc(const Call& C)
{
	const double Al = C.N(0), Be = C.N(1), Ga = C.N(2), P0 = C.N(3);
	if (Ga == 0.0)
	{
		return Result::Of(CtlInf);
	}
	const double U0 = Ga * P0 + Be, Om2 = Al * Ga - Be * Be;
	const double Tiny = 1e-14 * std::max(1.0, Al * Ga + Be * Be);
	if (Om2 > Tiny)
	{
		const double Om = std::sqrt(Om2);
		return Result::Of((cst::Pi / 2 - std::atan(U0 / Om)) / Om);
	}
	if (Om2 >= -Tiny)
	{
		return Result::Of(U0 > 0 ? 1.0 / U0 : CtlInf);
	}
	const double Kk = std::sqrt(-Om2);
	if (U0 > Kk)
	{
		return Result::Of(std::log((U0 + Kk) / (U0 - Kk)) / (2 * Kk));
	}
	return Result::Of(CtlInf);
}

// H = [A, Q; C'R^-1C, -A'] of the continuous filter Riccati equation (chapter 3, Theorem 2).
bool FilterHamiltonian(const Matrix& A, const Matrix& Cm, const Matrix& Q, const Matrix& R, Matrix* H, std::string* Err)
{
	if (!CheckFilter(A, Cm, Q, R, Err))
	{
		return false;
	}
	Matrix CtRiC;
	if (!la::Solve(R, Cm, &CtRiC, Err))
	{
		*Err = "R must be invertible";
		return false;
	}
	CtRiC = la::Mul(la::Transpose(Cm), CtRiC);
	*H = Block(A, Q, CtRiC, la::Scale(la::Transpose(A), -1.0));
	return true;
}

Result Hamilton(const Call& C)
{
	Matrix H;
	std::string E;
	if (!FilterHamiltonian(C.M(0), C.M(1), C.M(2), C.M(3), &H, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(H);
}

// P = X Y^-1 from the invariant subspace of H for one half of the spectrum: side = +1 (right half,
// the stabilizing P > 0 in this sign convention) or -1 (left half, an exact ARE solution with P < 0).
Result HamAre(const Call& C)
{
	const double Side = C.N(4);
	if (Side != 1.0 && Side != -1.0)
	{
		return Result::Fail("side must be 1 (right half: stabilizing) or -1 (left half)");
	}
	Matrix H, S;
	std::string E;
	if (!FilterHamiltonian(C.M(0), C.M(1), C.M(2), C.M(3), &H, &E) || !MatrixSign(H, &S, &E))
	{
		return Result::Fail(E);
	}
	const int N = C.M(0).Rows;
	// [S11 - side I; S21] P = -[S12; S22 - side I]
	Matrix Lhs(2 * N, N), Rhs(2 * N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			const double Dij = i == j ? Side : 0.0;
			Lhs(i, j) = S(i, j) - Dij;
			Lhs(i + N, j) = S(i + N, j);
			Rhs(i, j) = -S(i, j + N);
			Rhs(i + N, j) = -(S(i + N, j + N) - Dij);
		}
	}
	Matrix P;
	if (!la::LeastSquares(Lhs, Rhs, &P, &E))
	{
		return Result::Fail("the invariant subspace is not a graph over the costate (Y singular)");
	}
	P = la::Symmetrize(P);
	// Residual of A P + P A' + Q - P C'R^-1C P = 0, with the blocks of H.
	const Matrix Ah = SubBlock(H, 0, 0, N, N), Qh = SubBlock(H, 0, N, N, N), Gh = SubBlock(H, N, 0, N, N);
	const Matrix Res = la::Sub(la::Add(la::Add(la::Mul(Ah, P), la::Mul(P, la::Transpose(Ah))), Qh), la::Mul(la::Mul(P, Gh), P));
	if (!(la::MaxAbs(Res) <= 1e-7 * std::max(1.0, std::max(la::MaxAbs(P), la::MaxAbs(Qh)))))
	{
		return Result::Fail("that half of the spectrum gives no solution (Y singular)");
	}
	return Result::Of(P);
}

// Z = [A' + G A^-1 Q, -G A^-1; -A^-1 Q, A^-1], G = C'R^-1C: the symplectic matrix of the regulator for the
// dual pair (A', C'); its eigenvalues inside the unit circle are the estimator poles (chapter 3).
Result SympZ(const Call& C)
{
	const Matrix &A = C.M(0), &Cm = C.M(1), &Q = C.M(2), &R = C.M(3);
	std::string E;
	if (!CheckFilter(A, Cm, Q, R, &E))
	{
		return Result::Fail(E);
	}
	Matrix Ai, RiC;
	if (!la::Inverse(A, &Ai, &E))
	{
		return Result::Fail("A must be invertible (the symplectic form needs A^-1)");
	}
	if (!la::Solve(R, Cm, &RiC, &E))
	{
		return Result::Fail("R must be invertible");
	}
	const Matrix G = la::Mul(la::Transpose(Cm), RiC);
	const Matrix GAi = la::Mul(G, Ai);
	return Result::Of(Block(la::Add(la::Transpose(A), la::Mul(GAi, Q)), la::Scale(GAi, -1.0), la::Scale(la::Mul(Ai, Q), -1.0), Ai));
}

// Proportional navigation in the linearized engagement: ZEM = y + v tgo, a = -N ZEM / tgo^2. Under PN the
// ZEM scales like tgo^N, so the effort integral of a^2 is N^2 ZEM^2 / ((2N - 3) tgo^3) (N = 3: minimum effort).
Result Pn(const Call& C)
{
	const double Y0 = C.N(0), V0 = C.N(1), Tgo = C.N(2), Nn = C.N(3);
	if (!(Tgo > 0))
	{
		return Result::Fail("tgo must be > 0");
	}
	const double Zem = Y0 + V0 * Tgo;
	const double Effort = Nn > 1.5 ? Nn * Nn * Zem * Zem / ((2 * Nn - 3) * Tgo * Tgo * Tgo) : CtlInf;
	return Result::Record({{"zem", Zem}, {"a0", -Nn * Zem / (Tgo * Tgo)}, {"effort", Effort}});
}

/* ---------- control: Riccati and Lyapunov solvers, discretization, MPC gain, margins ---------- */

// The stabilizing DARE solution, checked: A - BK must be stable.
bool SolveDare(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, Matrix* P, Matrix* K, std::string* Err)
{
	if (!CheckRegulator(A, B, Q, R, Err) || !la::Dare(A, B, Q, R, P, K, Err))
	{
		return false;
	}
	if (!(la::SpectralRadius(la::Sub(A, la::Mul(B, *K))) < 1))
	{
		*Err = "no stabilizing solution (is (A, B) stabilizable and (A, Q) detectable?)";
		return false;
	}
	return true;
}

// The stabilizing CARE solution, checked: residual small and A - BK Hurwitz.
bool SolveCare(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, Matrix* P, Matrix* K, std::string* Err)
{
	if (!CheckRegulator(A, B, Q, R, Err) || !la::Care(A, B, Q, R, P, K, Err))
	{
		return false;
	}
	const Matrix At = la::Transpose(A);
	const Matrix PB = la::Mul(*P, B);
	const Matrix Res = la::Add(la::Sub(la::Add(la::Mul(At, *P), la::Mul(*P, A)), la::Mul(PB, *K)), Q);
	const double Scale = std::max(1.0, std::max(la::MaxAbs(*P), la::MaxAbs(Q)));
	if (!(la::MaxAbs(Res) <= 1e-7 * Scale) || !(MaxRealPart(la::Sub(A, la::Mul(B, *K))) < 0))
	{
		*Err = "no stabilizing solution (is (A, B) stabilizable and (A, Q) detectable?)";
		return false;
	}
	return true;
}

Result DareFn(const Call& C)
{
	std::string E;
	Matrix P, K;
	if (!SolveDare(C.M(0), C.M(1), C.M(2), C.M(3), &P, &K, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(P);
}

Result CareFn(const Call& C)
{
	std::string E;
	Matrix P, K;
	if (!SolveCare(C.M(0), C.M(1), C.M(2), C.M(3), &P, &K, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(P);
}

Result Dlqr(const Call& C)
{
	std::string E;
	Matrix P, K;
	if (!SolveDare(C.M(0), C.M(1), C.M(2), C.M(3), &P, &K, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(K);
}

Result Lqr(const Call& C)
{
	std::string E;
	Matrix P, K;
	if (!SolveCare(C.M(0), C.M(1), C.M(2), C.M(3), &P, &K, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(K);
}

Result KalmanFn(const Call& C)
{
	std::string E;
	Matrix P, M;
	if (!CheckFilter(C.M(0), C.M(1), C.M(2), C.M(3), &E) || !SteadyKalman(C.M(0), C.M(1), C.M(2), C.M(3), &P, &M, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(M);
}

Result KfCov(const Call& C)
{
	std::string E;
	Matrix P, M;
	if (!CheckFilter(C.M(0), C.M(1), C.M(2), C.M(3), &E) || !SteadyKalman(C.M(0), C.M(1), C.M(2), C.M(3), &P, &M, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(P);
}

Result DlyapFn(const Call& C)
{
	Matrix X;
	std::string E;
	if (!la::Dlyap(C.M(0), C.M(1), &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}

Result LyapFn(const Call& C)
{
	Matrix X;
	std::string E;
	if (!la::Lyap(C.M(0), C.M(1), &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}

// Zero-order hold: expm([A B; 0 0] dt) = [Ad Bd; 0 I].
Result C2d(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1);
	const double Dt = C.N(2);
	if (!A.IsSquare())
	{
		return Result::Fail("A must be square, got " + ShapeText(A));
	}
	if (B.Rows != A.Rows)
	{
		return Result::Fail("B must have " + std::to_string(A.Rows) + " rows, got " + ShapeText(B));
	}
	if (!(Dt > 0))
	{
		return Result::Fail("dt must be > 0");
	}
	const int N = A.Rows, M = B.Cols;
	Matrix Aug(N + M, N + M);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Aug(i, j) = A(i, j) * Dt;
		}
		for (int j = 0; j < M; ++j)
		{
			Aug(i, N + j) = B(i, j) * Dt;
		}
	}
	Matrix Ex;
	std::string E;
	if (!la::Expm(Aug, &Ex, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(SubBlock(Ex, 0, 0, N, N + M));
}

// Van Loan: expm([-A W; 0 A'] dt) = [., F12; 0, F22] and Qd = F22' F12.
Result C2dq(const Call& C)
{
	const Matrix &A = C.M(0), &W = C.M(1);
	const double Dt = C.N(2);
	if (!A.IsSquare() || W.Rows != A.Rows || W.Cols != A.Cols)
	{
		return Result::Fail("A and W must be the same square size, got " + ShapeText(A) + " and " + ShapeText(W));
	}
	if (!(Dt > 0))
	{
		return Result::Fail("dt must be > 0");
	}
	const int N = A.Rows;
	Matrix Vl(2 * N, 2 * N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Vl(i, j) = -A(i, j) * Dt;
			Vl(i, N + j) = W(i, j) * Dt;
			Vl(N + i, N + j) = A(j, i) * Dt;
		}
	}
	Matrix F;
	std::string E;
	if (!la::Expm(Vl, &F, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(la::Symmetrize(la::Mul(la::Transpose(SubBlock(F, N, N, N, N)), SubBlock(F, 0, N, N, N))));
}

// First-move gain of unconstrained MPC with horizon N and terminal weight Pf: N - 1 Riccati steps from Pf.
Result MpcGain(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &Q = C.M(2), &R = C.M(3), &Pf = C.M(4);
	std::string E;
	if (!CheckRegulator(A, B, Q, R, &E))
	{
		return Result::Fail(E);
	}
	if (Pf.Rows != A.Rows || Pf.Cols != A.Rows)
	{
		return Result::Fail("Pf must be " + std::to_string(A.Rows) + " x " + std::to_string(A.Rows) + ", got " + ShapeText(Pf));
	}
	if (!IsWhole(C.N(5), 1, 10000))
	{
		return Result::Fail("N must be a whole number 1..10000");
	}
	const int N = static_cast<int>(C.N(5));
	Matrix P = Pf, K;
	for (int k = 1; k < N; ++k)
	{
		Matrix Next;
		if (!RiccatiStep(A, B, Q, R, P, &Next, nullptr, &E))
		{
			return Result::Fail(E);
		}
		P = Next;
	}
	if (!RiccatiStep(A, B, Q, R, P, nullptr, &K, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(K);
}

// Infinite-horizon cost sum x'Qx + u'Ru of u = -K x from x0: x0' X x0 with X = F'XF + Q + K'RK.
Result ClCost(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &Q = C.M(2), &R = C.M(3), &K = C.M(4);
	std::string E;
	if (!CheckRegulator(A, B, Q, R, &E))
	{
		return Result::Fail(E);
	}
	if (K.Rows != B.Cols || K.Cols != A.Rows)
	{
		return Result::Fail("K must be " + std::to_string(B.Cols) + " x " + std::to_string(A.Rows) + ", got " + ShapeText(K));
	}
	Matrix X0;
	if (!AsColumn(C.M(5), A.Rows, "x0", &X0, &E))
	{
		return Result::Fail(E);
	}
	const Matrix F = la::Sub(A, la::Mul(B, K));
	const double Rho = la::SpectralRadius(F);
	if (!(Rho < 1))
	{
		return Result::Fail("closed loop A - BK is not stable (spectral radius " + std::to_string(Rho) + "): the cost is infinite");
	}
	Matrix X;
	if (!la::Dlyap(la::Transpose(F), la::Add(Q, la::Mul(la::Mul(la::Transpose(K), R), K)), &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(la::Mul(la::Mul(la::Transpose(X0), X), X0).ScalarValue());
}

Result Gm(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &K = C.M(2);
	if (!A.IsSquare() || B.Rows != A.Rows || K.Rows != B.Cols || K.Cols != A.Rows)
	{
		return Result::Fail("A n x n, B n x m, K m x n; got " + ShapeText(A) + ", " + ShapeText(B) + ", " + ShapeText(K));
	}
	const double Lo = C.N(3), Hi = C.N(4);
	if (!(Lo > 0 && Lo < 1 && Hi > 1))
	{
		return Result::Fail("search range needs 0 < lo < 1 < hi");
	}
	const Matrix BK = la::Mul(B, K);
	double OutLo = 0.0, OutHi = 0.0;
	if (!GainInterval([&](double G) { return la::SpectralRadius(la::Sub(A, la::Scale(BK, G))); }, Lo, Hi, 400, &OutLo, &OutHi))
	{
		return Result::Fail("A - BK is not stable at g = 1");
	}
	return Result::Record({{"lo", OutLo}, {"hi", OutHi}});
}

// Guaranteed discrete LQR gain margin from the return-difference identity |1 + L| >= rho (SISO).
Result LqrGm(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &Q = C.M(2), &R = C.M(3);
	std::string E;
	if (!CheckRegulator(A, B, Q, R, &E))
	{
		return Result::Fail(E);
	}
	if (B.Cols != 1)
	{
		return Result::Fail("one input only (B must be a column)");
	}
	Matrix P, K;
	if (!SolveDare(A, B, Q, R, &P, &K, &E))
	{
		return Result::Fail(E);
	}
	const double Rv = R(0, 0);
	const double Rho = std::sqrt(Rv / (Rv + la::Mul(la::Mul(la::Transpose(B), P), B).ScalarValue()));
	return Result::Record({{"rho", Rho}, {"lo", 1 / (1 + Rho)}, {"hi", 1 / (1 - Rho)}});
}

// min over the unit circle of |1 + K (zI - A)^-1 B| on the bench's 600-point grid.
double ReturnDifferenceMin(const Matrix& A, const Matrix& B, const Matrix& K)
{
	const int N = A.Rows, Steps = 600;
	double Best = CtlInf;
	for (int i = 0; i < Steps; ++i)
	{
		const double W = 1e-4 + (cst::Pi - 1e-4) * i / (Steps - 1), Cw = std::cos(W), Sw = std::sin(W);
		Matrix Z(2 * N, 2 * N), Rhs(2 * N, 1);
		for (int r = 0; r < N; ++r)
		{
			for (int k = 0; k < N; ++k)
			{
				const double Re = (r == k ? Cw : 0) - A(r, k), Im = r == k ? Sw : 0;
				Z(r, k) = Re;
				Z(r, k + N) = -Im;
				Z(r + N, k) = Im;
				Z(r + N, k + N) = Re;
			}
			Rhs(r, 0) = B(r, 0);
		}
		Matrix X;
		std::string E;
		if (!la::Solve(Z, Rhs, &X, &E))
		{
			return 0.0;
		}
		double Re = 1, Im = 0;
		for (int k = 0; k < N; ++k)
		{
			Re += K(0, k) * X(k, 0);
			Im += K(0, k) * X(k + N, 0);
		}
		Best = std::min(Best, std::hypot(Re, Im));
	}
	return Best;
}

Result RetDiff(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &K = C.M(2);
	if (!A.IsSquare() || B.Rows != A.Rows || B.Cols != 1 || K.Rows != 1 || K.Cols != A.Rows)
	{
		return Result::Fail("SISO: A n x n, B n x 1, K 1 x n; got " + ShapeText(A) + ", " + ShapeText(B) + ", " + ShapeText(K));
	}
	return Result::Of(ReturnDifferenceMin(A, B, K));
}

Result LqgCl(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1), &Cm = C.M(2), &K = C.M(3), &M = C.M(4);
	const int N = A.Rows;
	if (!A.IsSquare() || B.Rows != N || Cm.Cols != N || K.Rows != B.Cols || K.Cols != N || M.Rows != N || M.Cols != Cm.Rows)
	{
		return Result::Fail("A n x n, B n x m, C p x n, K m x n, M n x p; got " + ShapeText(A) + ", " + ShapeText(B) + ", " + ShapeText(Cm) + ", " + ShapeText(K) + ", " + ShapeText(M));
	}
	return Result::Of(LqgLoop(A, B, Cm, K, M, C.N(5)));
}

// Doyle (1978): x' = [1 1; 0 1]x + [0 1]'u + [1 1]'w, y = x1 + v. K = f[1 1], L = d[1 1]'.
Result Doyle(const Call& C)
{
	const double Q = C.N(0);
	const double Sigma = std::isnan(C.N(1)) ? Q : C.N(1);
	if (!(Q >= 0) || !(Sigma >= 0))
	{
		return Result::Fail("q and sigma must be >= 0");
	}
	const double F = 2 + std::sqrt(4 + Q), D = 2 + std::sqrt(4 + Sigma);
	const Matrix A = Matrix::FromRows({{1, 1}, {0, 1}});
	const Matrix BK = Matrix::FromRows({{0, 0}, {F, F}});
	const Matrix LC = Matrix::FromRows({{D, 0}, {D, 0}});
	const Matrix Est = la::Sub(la::Sub(A, BK), LC);
	auto LqgFn = [&](double G) { return std::exp(MaxRealPart(Block(A, la::Scale(BK, -G), LC, Est))); };
	auto LqrFn = [&](double G) { return std::exp(MaxRealPart(la::Sub(A, la::Scale(BK, G)))); };
	double Lo = CtlNan, Hi = CtlNan, RLo = CtlNan, RHi = CtlNan;
	GainInterval(LqgFn, 1e-2, 1e2, 400, &Lo, &Hi);
	GainInterval(LqrFn, 1e-2, 1e3, 400, &RLo, &RHi);
	return Result::Record({{"f", F}, {"d", D}, {"lo", Lo}, {"hi", Hi}, {"rlo", RLo}, {"rhi", RHi}});
}

/* ---------- control: the approach bench (approach/approach.js) ---------- */

struct GlideDesign
{
	Matrix Ad, Bd, Qd, Cm, R, Q, P, K, Pf, M;
	double Ru = 0.0, Dt = 0.1;
};

bool DesignGlide(double Sg, double Tau, double Sd, double Sv, GlideDesign* D, std::string* Err)
{
	const double Za = 70.5, Wn = 4.37, Zeta = 0.6;
	const double Su = 0.003, Qd = 1, Qv = 1 / (1.5 * 1.5), Rdeg = 2 * cst::Pi / 180;
	const int N = 5;
	Matrix A(N, N), B(N, 1), W(N, N);
	A(0, 1) = 1;
	A(0, 4) = 1;
	A(1, 2) = Za * 1;
	A(2, 3) = 1;
	A(3, 2) = -Wn * Wn;
	A(3, 3) = -2 * Zeta * Wn;
	A(4, 4) = -1 / Tau;
	B(3, 0) = Wn * Wn;
	W(3, 3) = std::pow(Wn, 4) * (Su * Su);
	W(4, 4) = 2 * (Sg * Sg) / Tau;
	D->Dt = 0.1;
	Matrix Aug(N + 1, N + 1);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Aug(i, j) = A(i, j) * D->Dt;
		}
		Aug(i, N) = B(i, 0) * D->Dt;
	}
	Matrix Ex;
	if (!la::Expm(Aug, &Ex, Err))
	{
		return false;
	}
	D->Ad = SubBlock(Ex, 0, 0, N, N);
	D->Bd = SubBlock(Ex, 0, N, N, 1);
	Matrix Vl(2 * N, 2 * N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Vl(i, j) = -A(i, j) * D->Dt;
			Vl(i, N + j) = W(i, j) * D->Dt;
			Vl(N + i, N + j) = A(j, i) * D->Dt;
		}
	}
	Matrix F;
	if (!la::Expm(Vl, &F, Err))
	{
		return false;
	}
	D->Qd = la::Symmetrize(la::Mul(la::Transpose(SubBlock(F, N, N, N, N)), SubBlock(F, 0, N, N, N)));
	D->Cm = Matrix::FromRows({{1, 0, 0, 0, 0}, {0, 1, 0, 0, 1}});
	D->R = Matrix::FromRows({{Sd * Sd, 0}, {0, Sv * Sv}});
	D->Q = Matrix(N, N);
	const double Hv[5] = {0, 1, 0, 0, 1};
	D->Q(0, 0) = Qd;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			D->Q(i, j) += Qv * Hv[i] * Hv[j];
		}
	}
	D->Ru = 1 / (Rdeg * Rdeg);
	if (!la::Dare(D->Ad, D->Bd, D->Q, Matrix::Scalar(D->Ru), &D->P, &D->K, Err))
	{
		return false;
	}
	return SteadyKalman(D->Ad, D->Cm, D->Qd, D->R, &D->Pf, &D->M, Err);
}

const double GlidePid[3] = {0.006, 0.020, 0.0006};

Matrix GlidePidLoop(const GlideDesign& D, double G)
{
	const double Kp = GlidePid[0], Kd = GlidePid[1], Ki = GlidePid[2];
	const Matrix BKC = la::Mul(la::Mul(D.Bd, Matrix::FromRows({{Kp, Kd}})), D.Cm);
	Matrix F(6, 6);
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			F(i, j) = D.Ad(i, j) - G * BKC(i, j);
		}
		F(i, 5) = -G * Ki * D.Bd(i, 0);
	}
	for (int j = 0; j < 5; ++j)
	{
		F(5, j) = D.Dt * D.Cm(0, j);
	}
	F(5, 5) = 1;
	return F;
}

// Stationary covariance of one closed loop (approach closedCov) and the gain interval of that loop.
Result Glide(const Call& C)
{
	const double Mode = C.N(0);
	if (Mode != 1 && Mode != 2 && Mode != 3)
	{
		return Result::Fail("mode: 1 LQR on the true state, 2 LQG, 3 PID");
	}
	if (!(C.N(1) > 0 && C.N(2) > 0 && C.N(3) > 0 && C.N(4) > 0))
	{
		return Result::Fail("sg, tau, sd and sv must be > 0");
	}
	GlideDesign D;
	std::string E;
	if (!DesignGlide(C.N(1), C.N(2), C.N(3), C.N(4), &D, &E))
	{
		return Result::Fail(E);
	}
	const Matrix I5 = Matrix::Identity(5);
	Matrix F, Wn, Gw, Gv;
	if (Mode == 1)
	{
		F = la::Sub(D.Ad, la::Mul(D.Bd, D.K));
		Wn = D.Qd;
	}
	else if (Mode == 2)
	{
		F = LqgLoop(D.Ad, D.Bd, D.Cm, D.K, D.M, 1.0);
		const Matrix BKM = la::Mul(la::Mul(D.Bd, D.K), D.M);
		Gw = Block(I5, Matrix(5, 0), Matrix(5, 5), Matrix(5, 0));
		Gv = Block(la::Scale(BKM, -1.0), Matrix(5, 0), la::Sub(la::Mul(D.Ad, D.M), BKM), Matrix(5, 0));
	}
	else
	{
		F = GlidePidLoop(D, 1.0);
		const Matrix BK = la::Mul(D.Bd, Matrix::FromRows({{GlidePid[0], GlidePid[1]}}));
		Gw = Block(I5, Matrix(5, 0), Matrix(1, 5), Matrix(1, 0));
		Gv = Block(la::Scale(BK, -1.0), Matrix(5, 0), Matrix::FromRows({{D.Dt, 0}}), Matrix(1, 0));
	}
	if (Mode != 1)
	{
		Wn = la::Add(la::Mul(la::Mul(Gw, D.Qd), la::Transpose(Gw)), la::Mul(la::Mul(Gv, D.R), la::Transpose(Gv)));
	}
	const double Rho = la::SpectralRadius(F);
	if (!(Rho < 1))
	{
		return Result::Fail("closed loop unstable (spectral radius " + std::to_string(Rho) + ")");
	}
	Matrix S;
	if (!la::Dlyap(F, Wn, &S, &E))
	{
		return Result::Fail(E);
	}
	const Matrix Sx = SubBlock(S, 0, 0, 5, 5);
	double Su = 0.0;
	if (Mode == 1)
	{
		Su = la::Mul(la::Mul(D.K, S), la::Transpose(D.K)).ScalarValue();
	}
	else if (Mode == 2)
	{
		const Matrix MC = la::Mul(D.M, D.Cm);
		const Matrix Gx = la::Scale(la::Mul(D.K, Block(MC, la::Sub(I5, MC), Matrix(0, 5), Matrix(0, 5))), -1.0);
		const Matrix KM = la::Mul(D.K, D.M);
		Su = la::Mul(la::Mul(Gx, S), la::Transpose(Gx)).ScalarValue() + la::Mul(la::Mul(KM, D.R), la::Transpose(KM)).ScalarValue();
	}
	else
	{
		Matrix Gx(1, 6);
		for (int j = 0; j < 5; ++j)
		{
			Gx(0, j) = -(GlidePid[0] * D.Cm(0, j) + GlidePid[1] * D.Cm(1, j));
		}
		Gx(0, 5) = -GlidePid[2];
		Su = la::Mul(la::Mul(Gx, S), la::Transpose(Gx)).ScalarValue() + (GlidePid[0] * GlidePid[0] * D.R(0, 0) + GlidePid[1] * GlidePid[1] * D.R(1, 1));
	}
	double J = D.Ru * Su;
	for (int i = 0; i < 5; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			J += D.Q(i, j) * Sx(j, i);
		}
	}
	double Lo = CtlNan, Hi = CtlNan;
	if (Mode == 1)
	{
		const Matrix BK = la::Mul(D.Bd, D.K);
		GainInterval([&](double G) { return la::SpectralRadius(la::Sub(D.Ad, la::Scale(BK, G))); }, 1e-3, 1e3, 400, &Lo, &Hi);
	}
	else if (Mode == 2)
	{
		GainInterval([&](double G) { return la::SpectralRadius(LqgLoop(D.Ad, D.Bd, D.Cm, D.K, D.M, G)); }, 1e-3, 1e3, 400, &Lo, &Hi);
	}
	else
	{
		GainInterval([&](double G) { return la::SpectralRadius(GlidePidLoop(D, G)); }, 1e-3, 1e3, 400, &Lo, &Hi);
	}
	const double Rq = std::sqrt(D.Ru / (D.Ru + la::Mul(la::Mul(la::Transpose(D.Bd), D.P), D.Bd).ScalarValue()));
	const double Rd = ReturnDifferenceMin(D.Ad, D.Bd, D.K);
	return Result::Record({{"sd", std::sqrt(S(0, 0))}, {"J", J}, {"su", std::sqrt(Su)}, {"rho", Rho}, {"lo", Lo}, {"hi", Hi}, {"rd", Rd}, {"rq", Rq}, {"glo", 1 / (1 + Rq)}, {"ghi", 1 / (1 - Rq)}});
}

/* ---------- control: the MPC bench (mpc/mpc.js), double integrator with a dual active-set QP ---------- */

using Vec = std::vector<double>;
using Rows = std::vector<Vec>;

struct MpcPlant
{
	double H = 0.2;
	double A[2][2] = {{1, 0.2}, {0, 1}};
	double B[2] = {0.02, 0.2};
	double Q[2][2] = {{1, 0}, {0, 0}};
};

MpcPlant MakeMpcPlant(double H)
{
	MpcPlant Pl;
	Pl.H = H;
	Pl.A[0][1] = H;
	Pl.B[0] = H * H / 2;
	Pl.B[1] = H;
	return Pl;
}

void MpcStep(const MpcPlant& Pl, double* X, double U)
{
	const double X0 = X[0] + Pl.H * X[1] + Pl.B[0] * U, X1 = X[1] + Pl.B[1] * U;
	X[0] = X0;
	X[1] = X1;
}

double MpcStage(const MpcPlant& Pl, const double* X, double U, double R)
{
	return X[0] * (Pl.Q[0][0] * X[0] + Pl.Q[0][1] * X[1]) + X[1] * (Pl.Q[1][0] * X[0] + Pl.Q[1][1] * X[1]) + R * U * U;
}

double MpcQuad(const double P[2][2], const double* X)
{
	return P[0][0] * X[0] * X[0] + 2 * P[0][1] * X[0] * X[1] + P[1][1] * X[1] * X[1];
}

double VDot(const Vec& A, const Vec& B)
{
	double S = 0;
	for (size_t i = 0; i < A.size(); ++i)
	{
		S += A[i] * B[i];
	}
	return S;
}

// DARE for the plant by Riccati iteration (mpc.js dare()).
void MpcDare(const MpcPlant& Pl, double R, double P[2][2], double K[2])
{
	const double Tol = 1e-14;
	double Pc[2][2] = {{Pl.Q[0][0], Pl.Q[0][1]}, {Pl.Q[1][0], Pl.Q[1][1]}};
	auto Map = [&](const double In[2][2], double Out[2][2]) {
		const double PB[2] = {In[0][0] * Pl.B[0] + In[0][1] * Pl.B[1], In[1][0] * Pl.B[0] + In[1][1] * Pl.B[1]};
		const double S = R + Pl.B[0] * PB[0] + Pl.B[1] * PB[1];
		const double APB[2] = {Pl.A[0][0] * PB[0] + Pl.A[1][0] * PB[1], Pl.A[0][1] * PB[0] + Pl.A[1][1] * PB[1]};
		for (int i = 0; i < 2; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				double Acc = 0;
				for (int k = 0; k < 2; ++k)
				{
					for (int l = 0; l < 2; ++l)
					{
						Acc += Pl.A[k][i] * In[k][l] * Pl.A[l][j];
					}
				}
				Out[i][j] = Pl.Q[i][j] + Acc - APB[i] * APB[j] / S;
			}
		}
	};
	for (int It = 0; It < 100000; ++It)
	{
		double Nx[2][2];
		Map(Pc, Nx);
		double Dmax = 0, Sc = 0;
		for (int i = 0; i < 2; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				Dmax = std::max(Dmax, std::fabs(Nx[i][j] - Pc[i][j]));
				Sc = std::max(Sc, std::fabs(Nx[i][j]));
			}
		}
		for (int i = 0; i < 2; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				Pc[i][j] = Nx[i][j];
			}
		}
		if (Dmax <= Tol * Sc)
		{
			break;
		}
	}
	const double PB[2] = {Pc[0][0] * Pl.B[0] + Pc[0][1] * Pl.B[1], Pc[1][0] * Pl.B[0] + Pc[1][1] * Pl.B[1]};
	const double S = R + Pl.B[0] * PB[0] + Pl.B[1] * PB[1];
	K[0] = (PB[0] * Pl.A[0][0] + PB[1] * Pl.A[1][0]) / S;
	K[1] = (PB[0] * Pl.A[0][1] + PB[1] * Pl.A[1][1]) / S;
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			P[i][j] = Pc[i][j];
		}
	}
}

bool MpcChol(const Rows& M, Rows* L)
{
	const size_t N = M.size();
	*L = Rows(N, Vec(N, 0.0));
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j <= i; ++j)
		{
			double S = M[i][j];
			for (size_t k = 0; k < j; ++k)
			{
				S -= (*L)[i][k] * (*L)[j][k];
			}
			if (i == j)
			{
				if (!(S > 0))
				{
					return false;
				}
				(*L)[i][i] = std::sqrt(S);
			}
			else
			{
				(*L)[i][j] = S / (*L)[j][j];
			}
		}
	}
	return true;
}

Vec MpcCholSolve(const Rows& L, const Vec& B)
{
	const size_t N = L.size();
	Vec Y(N, 0.0), X(N, 0.0);
	for (size_t i = 0; i < N; ++i)
	{
		double S = B[i];
		for (size_t k = 0; k < i; ++k)
		{
			S -= L[i][k] * Y[k];
		}
		Y[i] = S / L[i][i];
	}
	for (size_t ii = N; ii-- > 0;)
	{
		double S = Y[ii];
		for (size_t k = ii + 1; k < N; ++k)
		{
			S -= L[k][ii] * X[k];
		}
		X[ii] = S / L[ii][ii];
	}
	return X;
}

enum MpcTerminal
{
	MpcTermNone = 0,
	MpcTermStop = 1,
	MpcTermOrigin = 2,
};

// Condensed QP: U = u_0..u_{N-1}, cost 1/2 U'MU + c'U + c0, rows Nc U >= b (eq rows are equalities).
struct MpcProblem
{
	MpcPlant Pl;
	int N = 0;
	double R = 0, Lo = 0, Hi = 0, Pmin = 0;
	bool HasPmin = false;
	int Terminal = MpcTermNone;
	Rows Phi, Gam, M, Minv, F, Nc, W, G, GF;
	double E[2][2] = {{0, 0}, {0, 0}};
	std::vector<char> Eq, HasO;
	Vec BaseV, BaseO0, BaseO1, UnitS;
	std::vector<int> UnitK;
};

bool MpcBuild(const MpcPlant& Pl, int N, double R, const double Pf[2][2], double Lo, double Hi, bool HasPmin, double Pmin, int Terminal, MpcProblem* Pr)
{
	Pr->Pl = Pl;
	Pr->N = N;
	Pr->R = R;
	Pr->Lo = Lo;
	Pr->Hi = Hi;
	Pr->HasPmin = HasPmin;
	Pr->Pmin = Pmin;
	Pr->Terminal = Terminal;
	const auto& A = Pl.A;
	const auto& B = Pl.B;
	const auto& Q = Pl.Q;
	Rows& Phi = Pr->Phi;
	Rows& Gam = Pr->Gam;
	Phi = Rows(2 * N, Vec(2, 0.0));
	Gam = Rows(2 * N, Vec(N, 0.0));
	std::vector<std::array<double, 4>> Pw;
	Pw.push_back({1, 0, 0, 1});
	for (int k = 1; k <= N; ++k)
	{
		const std::array<double, 4> a = Pw[k - 1];
		Pw.push_back({A[0][0] * a[0] + A[0][1] * a[2], A[0][0] * a[1] + A[0][1] * a[3], A[1][0] * a[0] + A[1][1] * a[2], A[1][0] * a[1] + A[1][1] * a[3]});
	}
	for (int k = 0; k < N; ++k)
	{
		const std::array<double, 4>& Ak = Pw[k + 1];
		Phi[2 * k][0] = Ak[0];
		Phi[2 * k][1] = Ak[1];
		Phi[2 * k + 1][0] = Ak[2];
		Phi[2 * k + 1][1] = Ak[3];
		for (int j = 0; j <= k; ++j)
		{
			const std::array<double, 4>& p = Pw[k - j];
			Gam[2 * k][j] = p[0] * B[0] + p[1] * B[1];
			Gam[2 * k + 1][j] = p[2] * B[0] + p[3] * B[1];
		}
	}
	auto Wk = [&](int k, int i, int j) { return k == N - 1 ? Pf[i][j] : Q[i][j]; };
	Rows QG(2 * N, Vec(N, 0.0)), QP(2 * N, Vec(2, 0.0));
	for (int k = 0; k < N; ++k)
	{
		for (int j = 0; j < N; ++j)
		{
			QG[2 * k][j] = Wk(k, 0, 0) * Gam[2 * k][j] + Wk(k, 0, 1) * Gam[2 * k + 1][j];
			QG[2 * k + 1][j] = Wk(k, 1, 0) * Gam[2 * k][j] + Wk(k, 1, 1) * Gam[2 * k + 1][j];
		}
		for (int j = 0; j < 2; ++j)
		{
			QP[2 * k][j] = Wk(k, 0, 0) * Phi[2 * k][j] + Wk(k, 0, 1) * Phi[2 * k + 1][j];
			QP[2 * k + 1][j] = Wk(k, 1, 0) * Phi[2 * k][j] + Wk(k, 1, 1) * Phi[2 * k + 1][j];
		}
	}
	Pr->M = Rows(N, Vec(N, 0.0));
	Pr->F = Rows(N, Vec(2, 0.0));
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			Pr->E[i][j] = Q[i][j];
		}
	}
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			double S = 0;
			for (int k = 0; k < 2 * N; ++k)
			{
				S += Gam[k][i] * QG[k][j];
			}
			Pr->M[i][j] = 2 * (S + (i == j ? R : 0));
		}
		for (int j = 0; j < 2; ++j)
		{
			double S = 0;
			for (int k = 0; k < 2 * N; ++k)
			{
				S += Gam[k][i] * QP[k][j];
			}
			Pr->F[i][j] = 2 * S;
		}
	}
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			double S = 0;
			for (int k = 0; k < 2 * N; ++k)
			{
				S += Phi[k][i] * QP[k][j];
			}
			Pr->E[i][j] += S;
		}
	}
	auto AddRow = [&](const Vec& Row, bool IsEq, double V, bool WithO, double O0, double O1) {
		Pr->Nc.push_back(Row);
		Pr->Eq.push_back(IsEq ? 1 : 0);
		Pr->BaseV.push_back(V);
		Pr->HasO.push_back(WithO ? 1 : 0);
		Pr->BaseO0.push_back(O0);
		Pr->BaseO1.push_back(O1);
	};
	for (int k = 0; k < N; ++k)
	{
		Vec Ev(N, 0.0), Fv(N, 0.0);
		Ev[k] = -1;
		Fv[k] = 1;
		AddRow(Ev, false, -Hi, false, 0, 0);
		AddRow(Fv, false, Lo, false, 0, 0);
	}
	if (HasPmin)
	{
		for (int k = 0; k < N; ++k)
		{
			AddRow(Gam[2 * k], false, Pmin, true, Phi[2 * k][0], Phi[2 * k][1]);
		}
	}
	if (Terminal == MpcTermStop)
	{
		AddRow(Gam[2 * N - 1], false, 0, true, Phi[2 * N - 1][0], Phi[2 * N - 1][1]);
	}
	if (Terminal == MpcTermOrigin)
	{
		for (int i = 0; i < 2; ++i)
		{
			AddRow(Gam[2 * N - 2 + i], true, 0, true, Phi[2 * N - 2 + i][0], Phi[2 * N - 2 + i][1]);
		}
	}
	const size_t Mrows = Pr->Nc.size();
	Pr->UnitK.assign(Mrows, -1);
	Pr->UnitS.assign(Mrows, 0.0);
	for (size_t i = 0; i < Mrows; ++i)
	{
		int K = -1;
		bool Single = true;
		for (int j = 0; j < N; ++j)
		{
			if (Pr->Nc[i][j] != 0)
			{
				if (K >= 0)
				{
					Single = false;
					break;
				}
				K = j;
			}
		}
		if (Single && K >= 0 && std::fabs(Pr->Nc[i][K]) == 1)
		{
			Pr->UnitK[i] = K;
			Pr->UnitS[i] = Pr->Nc[i][K];
		}
	}
	Rows L;
	if (!MpcChol(Pr->M, &L))
	{
		return false;
	}
	Pr->Minv = Rows(N, Vec(N, 0.0));
	for (int j = 0; j < N; ++j)
	{
		Vec Ej(N, 0.0);
		Ej[j] = 1;
		const Vec Col = MpcCholSolve(L, Ej);
		for (int i = 0; i < N; ++i)
		{
			Pr->Minv[i][j] = Col[i];
		}
	}
	Pr->W = Rows(Mrows, Vec(N, 0.0));
	Pr->G = Rows(Mrows, Vec(Mrows, 0.0));
	for (size_t i = 0; i < Mrows; ++i)
	{
		if (Pr->UnitK[i] >= 0)
		{
			for (int j = 0; j < N; ++j)
			{
				Pr->W[i][j] = Pr->Minv[Pr->UnitK[i]][j] * Pr->UnitS[i];
			}
		}
		else
		{
			for (int j = 0; j < N; ++j)
			{
				Pr->W[i][j] = VDot(Pr->Minv[j], Pr->Nc[i]);
			}
		}
	}
	for (size_t i = 0; i < Mrows; ++i)
	{
		for (size_t j = 0; j <= i; ++j)
		{
			const double Gij = Pr->UnitK[i] >= 0 ? Pr->UnitS[i] * Pr->W[j][Pr->UnitK[i]] : VDot(Pr->Nc[i], Pr->W[j]);
			Pr->G[i][j] = Gij;
			Pr->G[j][i] = Gij;
		}
	}
	return true;
}

struct MpcQp
{
	Vec C, B;
	double C0 = 0;
};

MpcQp MpcAt(const MpcProblem& Pr, const double* X0)
{
	MpcQp Q;
	Q.C.resize(Pr.F.size());
	for (size_t i = 0; i < Pr.F.size(); ++i)
	{
		Q.C[i] = Pr.F[i][0] * X0[0] + Pr.F[i][1] * X0[1];
	}
	Q.C0 = MpcQuad(Pr.E, X0);
	Q.B.resize(Pr.BaseV.size());
	for (size_t i = 0; i < Pr.BaseV.size(); ++i)
	{
		Q.B[i] = Pr.HasO[i] ? Pr.BaseV[i] - Pr.BaseO0[i] * X0[0] - Pr.BaseO1[i] * X0[1] : Pr.BaseV[i];
	}
	return Q;
}

double MpcCost(const MpcProblem& Pr, const MpcQp& Q, const Vec& U)
{
	Vec MU(U.size());
	for (size_t i = 0; i < U.size(); ++i)
	{
		MU[i] = VDot(Pr.M[i], U);
	}
	return 0.5 * VDot(U, MU) + VDot(Q.C, U) + Q.C0;
}

struct MpcSolution
{
	bool Ok = false;
	Vec U;
	int Active = 0;
};

// Goldfarb-Idnani dual active-set method (mpc.js gi()), line for line.
MpcSolution MpcGi(const MpcProblem& Pr, const MpcQp& Qp)
{
	MpcSolution Out;
	const Rows& Minv = Pr.Minv;
	const Rows& Nc = Pr.Nc;
	const Rows& W = Pr.W;
	const Rows& G = Pr.G;
	const Vec& Cv = Qp.C;
	const Vec& Bv = Qp.B;
	const int N = static_cast<int>(Cv.size()), M = static_cast<int>(Bv.size()), Cap = N + 2;
	Vec X(N, 0.0);
	for (int i = 0; i < N; ++i)
	{
		double S = 0;
		for (int j = 0; j < N; ++j)
		{
			S -= Minv[i][j] * Cv[j];
		}
		X[i] = S;
	}
	std::vector<int> Act(Cap, 0);
	Vec Sg(M, 1.0);
	std::vector<char> IsAct(M, 0);
	Vec Gi(static_cast<size_t>(Cap) * Cap, 0.0), Uu(Cap, 0.0), Up(Cap, 0.0), Rr(Cap, 0.0), Gg(Cap, 0.0), Z(N, 0.0), Col(Cap, 0.0), RowK(Cap, 0.0);
	int Qa = 0, Iters = 0;
	Vec TolV(M);
	for (int i = 0; i < M; ++i)
	{
		TolV[i] = 1e-9 * (1 + std::fabs(Bv[i]));
	}
	auto At = [&](int i, int j) -> double& { return Gi[static_cast<size_t>(i) * Cap + j]; };
	auto Drop = [&](int Kd) {
		const double Piv = At(Kd, Kd);
		for (int i = 0; i < Qa; ++i)
		{
			Col[i] = At(i, Kd);
			RowK[i] = At(Kd, i);
		}
		for (int i = 0; i < Qa; ++i)
		{
			for (int j = 0; j < Qa; ++j)
			{
				At(i, j) -= Col[i] * RowK[j] / Piv;
			}
		}
		for (int i = 0, ii = 0; i < Qa; ++i)
		{
			if (i == Kd)
			{
				continue;
			}
			for (int j = 0, jj = 0; j < Qa; ++j)
			{
				if (j == Kd)
				{
					continue;
				}
				At(ii, jj) = At(i, j);
				++jj;
			}
			++ii;
		}
		IsAct[Act[Kd]] = 0;
		for (int j = Kd; j < Qa - 1; ++j)
		{
			Act[j] = Act[j + 1];
			Up[j] = Up[j + 1];
		}
		--Qa;
	};
	auto Slack = [&](int i) {
		const int Uk = Pr.UnitK[i];
		return (Uk >= 0 ? Pr.UnitS[i] * X[Uk] : VDot(Nc[i], X)) - Bv[i];
	};
	for (int Outer = 0; Outer < 100 * (M + 1); ++Outer)
	{
		int P = -1;
		double Worst = 0;
		for (int i = 0; i < M; ++i)
		{
			if (IsAct[i])
			{
				continue;
			}
			if (Pr.Eq[i])
			{
				P = i;
				break;
			}
			const double S = Slack(i);
			if (S < -TolV[i] && S / (1 + std::fabs(Bv[i])) < Worst)
			{
				Worst = S / (1 + std::fabs(Bv[i]));
				P = i;
			}
		}
		if (P < 0)
		{
			Out.Ok = true;
			Out.U = X;
			Out.Active = Qa;
			return Out;
		}
		if (Qa >= Cap)
		{
			return Out;
		}
		if (Pr.Eq[P] && Slack(P) > 0)
		{
			Sg[P] = -1;
		}
		for (int j = 0; j < Qa; ++j)
		{
			Up[j] = Uu[j];
		}
		Up[Qa] = 0;
		for (;;)
		{
			++Iters;
			if (Iters > 100 * (M + N + 1))
			{
				return Out;
			}
			const double Sp = Sg[P] * Slack(P);
			const Vec& Wp = W[P];
			const Vec& Gp = G[P];
			for (int j = 0; j < Qa; ++j)
			{
				Gg[j] = Gp[Act[j]] * Sg[Act[j]] * Sg[P];
			}
			for (int i = 0; i < Qa; ++i)
			{
				double S = 0;
				for (int j = 0; j < Qa; ++j)
				{
					S += At(i, j) * Gg[j];
				}
				Rr[i] = S;
			}
			for (int k = 0; k < N; ++k)
			{
				Z[k] = Wp[k] * Sg[P];
			}
			for (int j = 0; j < Qa; ++j)
			{
				const Vec& Wj = W[Act[j]];
				const double F = Sg[Act[j]] * Rr[j];
				for (int k = 0; k < N; ++k)
				{
					Z[k] -= Wj[k] * F;
				}
			}
			double T1 = CtlInf;
			int Kd = -1;
			for (int j = 0; j < Qa; ++j)
			{
				if (!Pr.Eq[Act[j]] && Rr[j] > 1e-12 && Up[j] / Rr[j] < T1)
				{
					T1 = Up[j] / Rr[j];
					Kd = j;
				}
			}
			double Zn = Gp[P];
			for (int j = 0; j < Qa; ++j)
			{
				Zn -= Gg[j] * Rr[j];
			}
			if (Zn <= 1e-12 * Gp[P])
			{
				if (Kd < 0)
				{
					return Out;   // infeasible
				}
				for (int j = 0; j < Qa; ++j)
				{
					Up[j] -= T1 * Rr[j];
				}
				Up[Qa] += T1;
				const double Last = Up[Qa];
				Drop(Kd);
				Up[Qa] = Last;
				continue;
			}
			const double T2 = -Sp / Zn, T = std::min(T1, T2);
			for (int k = 0; k < N; ++k)
			{
				X[k] += T * Z[k];
			}
			for (int j = 0; j < Qa; ++j)
			{
				Up[j] -= T * Rr[j];
			}
			Up[Qa] += T;
			if (T2 <= T1)
			{
				for (int i = 0; i < Qa; ++i)
				{
					for (int j = 0; j < Qa; ++j)
					{
						At(i, j) += Rr[i] * Rr[j] / Zn;
					}
					At(i, Qa) = -Rr[i] / Zn;
					At(Qa, i) = -Rr[i] / Zn;
				}
				At(Qa, Qa) = 1 / Zn;
				Act[Qa] = P;
				IsAct[P] = 1;
				++Qa;
				for (int j = 0; j < Qa; ++j)
				{
					Uu[j] = Up[j];
				}
				break;
			}
			const double Last = Up[Qa];
			Drop(Kd);
			Up[Qa] = Last;
		}
	}
	return Out;
}

struct MpcRun
{
	int Fail = -1;      // step at which the QP had no solution, or -1
	double J = CtlNan;  // closed-loop cost (NaN after a failure)
	double MinPos = CtlInf;
	double V0 = CtlNan; // open-loop value of the first plan
	double X[2] = {0, 0};
};

// The linear law once no row is active: if the whole plan stays inside the box along its own rollout,
// the remaining cost is that rollout's sum (mpc.js linearTail()).
bool MpcLinearTail(MpcProblem& Pr, const double* X0, int Room, MpcRun* Run, double* Jtail)
{
	const int N = Pr.N;
	if (Pr.GF.empty())
	{
		Pr.GF = Rows(N, Vec(2, 0.0));
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				double S = 0;
				for (int k = 0; k < N; ++k)
				{
					S += Pr.Minv[i][k] * Pr.F[k][j];
				}
				Pr.GF[i][j] = S;
			}
		}
	}
	double X[2] = {X0[0], X0[1]};
	double J = 0, MinPos = CtlInf;
	for (int k = 0; k < 5000; ++k)
	{
		for (int i = 0; i < N; ++i)
		{
			const double Ui = -(Pr.GF[i][0] * X[0] + Pr.GF[i][1] * X[1]);
			if (Ui > Pr.Hi + 1e-12 || Ui < Pr.Lo - 1e-12)
			{
				return false;
			}
		}
		const double U = -(Pr.GF[0][0] * X[0] + Pr.GF[0][1] * X[1]);
		J += MpcStage(Pr.Pl, X, U, Pr.R);
		MpcStep(Pr.Pl, X, U);
		if (k < Room)
		{
			MinPos = std::min(MinPos, X[0]);
		}
		if (std::max(std::fabs(X[0]), std::fabs(X[1])) < 1e-13 * (1 + std::fabs(X0[0]) + std::fabs(X0[1])))
		{
			Run->MinPos = std::min(Run->MinPos, MinPos);
			*Jtail = J;
			return true;
		}
	}
	return false;
}

// Receding-horizon closed loop (mpc.js closedLoop()), T steps, cost closed with x'Px.
MpcRun MpcClosedLoop(MpcProblem& Pr, const double* X0, const double P[2][2], int T, bool Tail)
{
	MpcRun Run;
	double X[2] = {X0[0], X0[1]};
	double J = 0;
	int Skip = 0;
	Run.MinPos = X[0];
	const bool CanTail = Tail && !Pr.HasPmin && Pr.Terminal == MpcTermNone;
	for (int t = 0; t < T; ++t)
	{
		const MpcQp Q = MpcAt(Pr, X);
		const MpcSolution S = MpcGi(Pr, Q);
		if (!S.Ok)
		{
			Run.Fail = t;
			Run.X[0] = X[0];
			Run.X[1] = X[1];
			return Run;
		}
		if (t == 0)
		{
			Run.V0 = MpcCost(Pr, Q, S.U);
		}
		if (CanTail && S.Active == 0)
		{
			const bool Try = Skip <= 0;
			--Skip;
			if (Try)
			{
				double Jt = 0;
				if (MpcLinearTail(Pr, X, T - t, &Run, &Jt))
				{
					Run.J = J + Jt;
					Run.X[0] = 0;
					Run.X[1] = 0;
					return Run;
				}
				Skip = 5;
			}
		}
		J += MpcStage(Pr.Pl, X, S.U[0], Pr.R);
		MpcStep(Pr.Pl, X, S.U[0]);
		Run.MinPos = std::min(Run.MinPos, X[0]);
		if (!(std::fabs(X[0]) < 1e9))
		{
			Run.J = CtlInf;
			return Run;
		}
	}
	Run.X[0] = X[0];
	Run.X[1] = X[1];
	Run.J = J + MpcQuad(P, X);
	return Run;
}

Result LaneSat(const Call& C)
{
	const double Um = C.N(0), Y0 = C.N(1), R = C.N(2);
	if (!(Um > 0) || !(R > 0))
	{
		return Result::Fail("umax and r must be > 0");
	}
	const MpcPlant Pl = MakeMpcPlant(0.2);
	double P[2][2], K[2];
	MpcDare(Pl, R, P, K);
	double X[2] = {Y0, 0};
	double J = 0, MinPos = X[0];
	for (int t = 0; t < 150; ++t)
	{
		const double U = std::min(Um, std::max(-Um, -(K[0] * X[0] + K[1] * X[1])));
		J += MpcStage(Pl, X, U, R);
		MpcStep(Pl, X, U);
		MinPos = std::min(MinPos, X[0]);
	}
	return Result::Record({{"J", J + MpcQuad(P, X)}, {"os", std::max(0.0, -MinPos)}});
}

Result LaneMpc(const Call& C)
{
	if (!IsWhole(C.N(0), 1, 60))
	{
		return Result::Fail("N must be a whole number 1..60");
	}
	const int N = static_cast<int>(C.N(0));
	const double Um = C.N(1), Y0 = C.N(2), R = C.N(3), Term = C.N(4);
	if (!(Um > 0) || !(R > 0))
	{
		return Result::Fail("umax and r must be > 0");
	}
	if (Term != 0 && Term != 1)
	{
		return Result::Fail("term: 1 terminal cost P (DARE), 0 none");
	}
	const MpcPlant Pl = MakeMpcPlant(0.2);
	double P[2][2], K[2];
	MpcDare(Pl, R, P, K);
	const double Zero[2][2] = {{0, 0}, {0, 0}};
	MpcProblem Pr;
	if (!MpcBuild(Pl, N, R, Term == 1 ? P : Zero, -Um, Um, false, 0, MpcTermNone, &Pr))
	{
		return Result::Fail("condensed Hessian is not positive definite");
	}
	const double X0[2] = {Y0, 0};
	const MpcRun Run = MpcClosedLoop(Pr, X0, P, 150, true);
	if (Run.Fail >= 0)
	{
		return Result::Fail("QP infeasible at step " + std::to_string(Run.Fail));
	}
	return Result::Record({{"J", Run.J}, {"os", std::max(0.0, -Run.MinPos)}, {"V", Run.V0}});
}

Result Lander(const Call& C)
{
	if (!IsWhole(C.N(0), 1, 60))
	{
		return Result::Fail("N must be a whole number 1..60");
	}
	const int N = static_cast<int>(C.N(0));
	const double H0 = C.N(1), V0 = C.N(2), Uh = C.N(3), Mode = C.N(4);
	if (!(Uh > 0))
	{
		return Result::Fail("uh (braking limit) must be > 0");
	}
	if (Mode != 0 && Mode != 1 && Mode != 2 && Mode != 3)
	{
		return Result::Fail("mode: 0 no terminal cost, 1 terminal cost P, 2 P and v_N >= 0, 3 x_N = 0");
	}
	const double R = 0.1;
	const MpcPlant Pl = MakeMpcPlant(0.2);
	double P[2][2], K[2];
	MpcDare(Pl, R, P, K);
	const double Zero[2][2] = {{0, 0}, {0, 0}};
	const int Term = Mode == 2 ? MpcTermStop : (Mode == 3 ? MpcTermOrigin : MpcTermNone);
	MpcProblem Pr;
	if (!MpcBuild(Pl, N, R, Mode == 0 ? Zero : P, -1, Uh, true, 0, Term, &Pr))
	{
		return Result::Fail("condensed Hessian is not positive definite");
	}
	const double X0[2] = {H0, -V0};
	const MpcRun Run = MpcClosedLoop(Pr, X0, P, 150, false);
	return Result::Record({{"fail", static_cast<double>(Run.Fail)}, {"alt", Run.X[0]}, {"vel", Run.X[1]}, {"J", Run.J}});
}

/* ---------- dynamics: chapter 8 ---------- */

// RK4 step as the chapter page writes it.
template <typename Field>
void Rk4Step(const Field& F, double* X, double Dt)
{
	double K1[2], K2[2], K3[2], K4[2], T[2];
	F(X, K1);
	T[0] = X[0] + Dt / 2 * K1[0];
	T[1] = X[1] + Dt / 2 * K1[1];
	F(T, K2);
	T[0] = X[0] + Dt / 2 * K2[0];
	T[1] = X[1] + Dt / 2 * K2[1];
	F(T, K3);
	T[0] = X[0] + Dt * K3[0];
	T[1] = X[1] + Dt * K3[1];
	F(T, K4);
	X[0] = X[0] + Dt / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
	X[1] = X[1] + Dt / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
}

// Area inside the unstable limit cycle of the reversed Van der Pol oscillator: one loop of the cycle,
// traced as a stable cycle in reverse time, by the shoelace formula (the page's cycle()).
double VdpBasinArea()
{
	auto Back = [](const double* Y, double* Out) {
		Out[0] = -(-Y[1]);
		Out[1] = -(Y[0] + (Y[0] * Y[0] - 1) * Y[1]);
	};
	double X[2] = {0.5, 0};
	for (int i = 0; i < 5000; ++i)
	{
		Rk4Step(Back, X, 0.01);
	}
	std::vector<double> Px{X[0]}, Py{X[1]};
	double Acc = 0, Prev = std::atan2(X[1], X[0]);
	for (int i = 0; i < 4000; ++i)
	{
		Rk4Step(Back, X, 0.005);
		const double A = std::atan2(X[1], X[0]);
		double D = A - Prev;
		if (D > cst::Pi)
		{
			D -= 2 * cst::Pi;
		}
		if (D < -cst::Pi)
		{
			D += 2 * cst::Pi;
		}
		Acc += D;
		Prev = A;
		Px.push_back(X[0]);
		Py.push_back(X[1]);
		if (std::fabs(Acc) >= 2 * cst::Pi)
		{
			break;
		}
	}
	double Area = 0;
	const size_t N = Px.size();
	for (size_t i = 0; i < N; ++i)
	{
		const size_t j = (i + 1) % N;
		Area += Px[i] * Py[j] - Py[i] * Px[j];
	}
	return std::fabs(Area) / 2;
}

// Largest sublevel set {x'Px < c} of the linearization's Lyapunov function inside {Vdot < 0}. On the ray
// x = r(cos t, sin t), Vdot = r^2 (-(c^2 + q s^2) + 2 r^2 c^2 s (P12 c + P22 s)), so the first crossing and
// its level are closed forms; c* is their minimum over t (rays out to r = 4, as on the page).
Result Roa(const Call& C)
{
	const double Qw = C.N(0);
	if (!(Qw > 0))
	{
		return Result::Fail("q must be > 0");
	}
	const Matrix A = Matrix::FromRows({{0, -1}, {1, -1}});
	Matrix P;
	std::string E;
	if (!la::Lyap(la::Transpose(A), Matrix::FromRows({{1, 0}, {0, Qw}}), &P, &E))
	{
		return Result::Fail(E);
	}
	const double P11 = P(0, 0), P12 = P(0, 1), P22 = P(1, 1);
	auto Level = [&](double T) {
		const double Cs = std::cos(T), Sn = std::sin(T);
		const double Den = 2 * Cs * Cs * Sn * (P12 * Cs + P22 * Sn);
		if (!(Den > 0))
		{
			return CtlInf;
		}
		const double R2 = (Cs * Cs + Qw * Sn * Sn) / Den;
		if (R2 > 16)
		{
			return CtlInf;
		}
		return R2 * (P11 * Cs * Cs + 2 * P12 * Cs * Sn + P22 * Sn * Sn);
	};
	const int Grid = 7200;
	double Best = CtlInf;
	int BestI = -1;
	for (int i = 0; i < Grid; ++i)
	{
		const double L = Level(cst::Pi * i / Grid);
		if (L < Best)
		{
			Best = L;
			BestI = i;
		}
	}
	if (BestI < 0)
	{
		return Result::Fail("no ray meets Vdot >= 0 within r = 4");
	}
	double Lo = cst::Pi * (BestI - 1) / Grid, Hi = cst::Pi * (BestI + 1) / Grid;
	const double Gr = (std::sqrt(5.0) - 1) / 2;
	double X1 = Hi - Gr * (Hi - Lo), X2 = Lo + Gr * (Hi - Lo), F1 = Level(X1), F2 = Level(X2);
	for (int It = 0; It < 200 && Hi - Lo > 1e-15; ++It)
	{
		if (F1 < F2)
		{
			Hi = X2;
			X2 = X1;
			F2 = F1;
			X1 = Hi - Gr * (Hi - Lo);
			F1 = Level(X1);
		}
		else
		{
			Lo = X1;
			X1 = X2;
			F1 = F2;
			X2 = Lo + Gr * (Hi - Lo);
			F2 = Level(X2);
		}
	}
	const double Cstar = std::min(Best, std::min(F1, F2));
	const double Area = cst::Pi * Cstar / std::sqrt(P11 * P22 - P12 * P12);
	const double Basin = VdpBasinArea();
	return Result::Record({{"c", Cstar}, {"area", Area}, {"basin", Basin}, {"frac", Area / Basin}});
}

// Nomoto ship T r' + r = K (delta + wave sin(0.5 t)), rudder limited to 35 deg (the page's ship()).
Result Nomoto(const Call& C)
{
	const double Ctrl = C.N(0);
	if (Ctrl != 1 && Ctrl != 2 && Ctrl != 3 && Ctrl != 4)
	{
		return Result::Fail("ctrl: 1 PID, 2 PID + anti-windup, 3 sliding mode, 4 sliding mode with a sign law");
	}
	const double Sk = 0.2, St = 15, Sdmax = 35 * cst::Pi / 180;
	const double Psid = C.N(1), Wv = C.N(2), Tf = 200, Dt = 0.02;
	const double Kp = 1.2, Kd = 12, Ki = 0.02, Lam = 0.1, Kk = 0.006, Phi = 0.004;
	const double StepDeg = Psid * 180 / cst::Pi;
	double Psi = 0, R = 0, I = 0, MaxDeg = -CtlInf, Settle = 0, Reach = CtlNan, S0 = 0;
	int Sw = 0, Last = 0;
	const int Steps = static_cast<int>(std::round(Tf / Dt));
	for (int i = 0; i < Steps; ++i)
	{
		const double T = i * Dt, E = Psi - Psid, S = R + Lam * E;
		double Dsat = 0;
		if (Ctrl <= 2)
		{
			const double D = -(Kp * E + Kd * R + Ki * I);
			Dsat = std::max(-Sdmax, std::min(Sdmax, D));
			if (Ctrl == 1 || D == Dsat)
			{
				I += E * Dt;
			}
		}
		else
		{
			const double Sign = S > 0 ? 1.0 : (S < 0 ? -1.0 : 0.0);
			const double Sat = Ctrl == 3 ? (std::fabs(S) < Phi ? S / Phi : Sign) : Sign;
			Dsat = std::max(-Sdmax, std::min(Sdmax, (R - Lam * St * R - St * Kk * Sat) / Sk));
		}
		const int Sg = Dsat > 0 ? 1 : -1;
		if (Sg != Last && i > 0)
		{
			++Sw;
		}
		Last = Sg;
		R += (-R + Sk * (Dsat + Wv * std::sin(0.5 * T))) / St * Dt;
		Psi += R * Dt;
		// The page records (t, psi, r, dsat, s) with psi already advanced and s from before the step.
		const double Deg = Psi * 180 / cst::Pi;
		MaxDeg = std::max(MaxDeg, Deg);
		if (std::fabs(Deg - StepDeg) > 2)
		{
			Settle = T;
		}
		if (i == 0)
		{
			S0 = std::fabs(S);
		}
		if (std::isnan(Reach) && std::fabs(S) < 0.004)
		{
			Reach = T;
		}
	}
	const double Eta = Kk - Sk * Wv / St;
	const double Over = std::max(0.0, MaxDeg - StepDeg) * cst::Pi / 180;
	return Result::Record({{"over", Over}, {"settle", Settle}, {"sw", static_cast<double>(Sw)}, {"reach", Reach}, {"bound", Eta > 0 ? S0 / Eta : CtlInf}});
}

Result LogMap(const Call& C)
{
	if (!IsWhole(C.N(2), 0, 1e6))
	{
		return Result::Fail("n must be a whole number 0..1e6");
	}
	const double R = C.N(0);
	double X = C.N(1);
	const long long N = static_cast<long long>(C.N(2));
	for (long long i = 0; i < N; ++i)
	{
		X = R * X * (1 - X);
	}
	return Result::Of(X);
}

Result LyapExp(const Call& C)
{
	if (!IsWhole(C.N(1), 1, 2e6))
	{
		return Result::Fail("n must be a whole number 1..2e6");
	}
	const double R = C.N(0);
	const int N = static_cast<int>(C.N(1));
	double X = 0.3, S = 0;
	for (int i = 0; i < 500; ++i)
	{
		X = R * X * (1 - X);
	}
	for (int i = 0; i < N; ++i)
	{
		X = R * X * (1 - X);
		S += std::log(std::fabs(R * (1 - 2 * X)) + 1e-300);
	}
	return Result::Of(S / N);
}

Result LogPeriod(const Call& C)
{
	const double R = C.N(0);
	double X = 0.3;
	for (int i = 0; i < 4000; ++i)
	{
		X = R * X * (1 - X);
	}
	const double X0 = X;
	for (int P = 1; P <= 64; ++P)
	{
		X = R * X * (1 - X);
		if (std::fabs(X - X0) < 1e-7)
		{
			return Result::Of(static_cast<double>(P));
		}
	}
	return Result::Of(0.0);
}

// Superstable parameters R_0 = 2, R_1 = 1 + sqrt 5, R_n by Newton on f^(2^n)(1/2) = 1/2 from the
// Feigenbaum extrapolation of the previous two (the page's feigenbaum()).
std::vector<double> SuperstableChain(int N)
{
	std::vector<double> Rs{2.0, 1 + std::sqrt(5.0)};
	for (int n = 2; n <= N; ++n)
	{
		const double Guess = Rs[Rs.size() - 1] + (Rs[Rs.size() - 1] - Rs[Rs.size() - 2]) / 4.669;
		const long long Iter = 1LL << n;
		auto F = [Iter](double Rr) {
			double X = 0.5;
			for (long long i = 0; i < Iter; ++i)
			{
				X = Rr * X * (1 - X);
			}
			return X - 0.5;
		};
		double R = Guess;
		for (int k = 0; k < 60; ++k)
		{
			const double Hh = 1e-9, Dr = F(R) / ((F(R + Hh) - F(R - Hh)) / (2 * Hh));
			R -= Dr;
			if (std::fabs(Dr) < 1e-14)
			{
				break;
			}
		}
		Rs.push_back(R);
	}
	return Rs;
}

Result Superstab(const Call& C)
{
	if (!IsWhole(C.N(0), 0, 12))
	{
		return Result::Fail("n must be a whole number 0..12");
	}
	const int N = static_cast<int>(C.N(0));
	return Result::Of(SuperstableChain(N)[static_cast<size_t>(N)]);
}

Result Feigen(const Call& C)
{
	if (!IsWhole(C.N(0), 2, 12))
	{
		return Result::Fail("n must be a whole number 2..12");
	}
	const int N = static_cast<int>(C.N(0));
	const std::vector<double> Rs = SuperstableChain(N);
	return Result::Of((Rs[N - 1] - Rs[N - 2]) / (Rs[N] - Rs[N - 1]));
}

// Damped pendulum th'' = -sin th - b th' by RK4 at dt = 0.01 (chapter 8, LaSalle).
Result Pendulum(const Call& C)
{
	const double B = C.N(1), Tend = C.N(2);
	if (!(Tend >= 0 && Tend <= 1000))
	{
		return Result::Fail("t must be 0..1000 s");
	}
	auto F = [B](const double* Z, double* Out) {
		Out[0] = Z[1];
		Out[1] = -std::sin(Z[0]) - B * Z[1];
	};
	double Z[2] = {C.N(0), 0};
	const long long Steps = std::llround(Tend / 0.01);
	for (long long i = 0; i < Steps; ++i)
	{
		Rk4Step(F, Z, 0.01);
	}
	return Result::Record({{"th", Z[0]}, {"w", Z[1]}, {"E", 0.5 * Z[1] * Z[1] + 1 - std::cos(Z[0])}});
}
} // namespace control_detail

void RegisterControlLibrary(Registry& R)
{
	using namespace control_detail;
	const char* Demo = "index Phase 1 demo: one map, three readings";
	const char* Appr = "approach/ bench section 1: LQR from the discrete Riccati equation";

	// riccati
	Def(R, "kfmap", "riccati", Demo, "one step of the scalar Kalman covariance map (a^2 r P/(h^2 P + r) + q as a Mobius map)", "P_{k+1}",
		{Num("P", "-"), Num("a", "-", 0.95), Num("q", "-", 0.2), Num("r", "-", 1), Num("h", "-", 1)}, &KfMap);
	Def(R, "kfiter", "riccati", Demo, "the scalar covariance map iterated k times from P0", "P_k",
		{Num("k", "steps"), Num("P0", "-", 6), Num("a", "-", 0.95), Num("q", "-", 0.2), Num("r", "-", 1), Num("h", "-", 1)}, &KfIter);
	Def(R, "kfss", "riccati", Demo, "steady state of the scalar Kalman map: attracting fixed point, gain, multiplier a^2(1-Kh)^2", "record p, pm (repelling), K, m, pole a(1-Kh)",
		{Num("a", "-", 0.95), Num("q", "-", 0.2), Num("r", "-", 1), Num("h", "-", 1)}, &KfSs);
	Def(R, "lqrss", "riccati", Demo, "the same map read as LQR (h -> b, r -> rho): cost-to-go, gain, closed-loop pole", "record S, Sm, K, m = (a-bK)^2, pole a-bK",
		{Num("a", "-", 0.95), Num("q", "-", 0.2), Num("rho", "-", 1), Num("b", "-", 1)}, &LqrSs);
	Def(R, "kfw", "riccati", "ch03 Theorem 1 and bench A: the cross-ratio", "cross-ratio w = (P - p+)/(P - p-); the iteration is w -> m w exactly", "w",
		{Num("P", "-"), Num("a", "-", 0.95), Num("q", "-", 0.2), Num("r", "-", 1), Num("h", "-", 1)}, &KfW);
	Def(R, "mobius", "riccati", Demo, "f(z) = (al z + be)/(ga z + de) on the projective line (inf allowed)", "f(z)",
		{Num("z", "-"), Num("al", "-", 0.5), Num("be", "-", -1), Num("ga", "-", 1), Num("de", "-", 0.5)}, &MobiusFn);
	Def(R, "mobiter", "riccati", Demo, "k-th iterate of a Mobius map from z0 (passes through inf)", "z_k",
		{Num("k", "steps"), Num("z0", "-", 2), Num("al", "-", 0.5), Num("be", "-", -1), Num("ga", "-", 1), Num("de", "-", 0.5)}, &MobIter);
	Def(R, "mobclass", "riccati", Demo, "classify a Mobius map: 0 degenerate, 1 identity, 2 parabolic, 3 hyperbolic, 4 elliptic", "record type, tr, det, disc, z1, z2, att, m, rot (rad/step, elliptic)",
		{Num("al", "-", 0.5), Num("be", "-", -1), Num("ga", "-", 1), Num("de", "-", 0.5)}, &MobClass);
	Def(R, "riccesc", "riccati", "ch03 counterexample 1: finite escape", "escape time of p' = al + 2 be p + ga p^2 from p0 (inf if none)", "s",
		{Num("al", "1/s", 1), Num("be", "1/s", 0.2), Num("ga", "1/s", 1), Num("p0", "-", 0)}, &RiccEsc);
	Def(R, "hamilton", "riccati", "ch03 Theorem 2: Riccati is linear in disguise", "Hamiltonian H = [A, Q; C'R^-1C, -A'] of the filter Riccati equation", "2n x 2n matrix",
		{Mat("A", "1/s"), Mat("C", "-"), Mat("Q", "-"), Mat("R", "-")}, &Hamilton);
	Def(R, "hamare", "riccati", "ch03 guise: continuous Riccati and the Hamiltonian", "ARE solution P = X Y^-1 from one half of eig(H): side 1 right (stabilizing), -1 left", "n x n matrix",
		{Mat("A", "1/s"), Mat("C", "-"), Mat("Q", "-"), Mat("R", "-"), Num("side", "-", 1)}, &HamAre);
	Def(R, "sympz", "riccati", "ch03 guise: matrix DARE and the symplectic matrix", "symplectic Z of the dual regulator; its eigenvalues inside the unit circle are the estimator poles", "2n x 2n matrix",
		{Mat("A", "-"), Mat("C", "-"), Mat("Q", "-"), Mat("R", "-")}, &SympZ);
	Def(R, "pn", "riccati", "index Phase 1 worked math: proportional navigation, N = 3", "PN in the linearized engagement: ZEM, first command, effort integral of a^2", "record zem (m), a0 (m/s^2), effort (m^2/s^3)",
		{Num("y0", "m", 30), Num("v0", "m/s", -2), Num("tgo", "s", 4), Num("N", "-", 3)}, &Pn);

	// control
	Def(R, "dare", "control", Appr, "discrete algebraic Riccati equation P = A'PA - A'PB(R + B'PB)^-1 B'PA + Q", "n x n matrix",
		{Mat("A", "-"), Mat("B", "-"), Mat("Q", "-"), Mat("R", "-")}, &DareFn);
	Def(R, "care", "control", Appr, "continuous algebraic Riccati equation A'P + PA - PBR^-1B'P + Q = 0", "n x n matrix",
		{Mat("A", "1/s"), Mat("B", "-"), Mat("Q", "-"), Mat("R", "-")}, &CareFn);
	Def(R, "dlqr", "control", Appr, "discrete LQR gain K = (R + B'PB)^-1 B'PA, u = -Kx", "m x n matrix",
		{Mat("A", "-"), Mat("B", "-"), Mat("Q", "-"), Mat("R", "-")}, &Dlqr);
	Def(R, "lqr", "control", Appr, "continuous LQR gain K = R^-1 B'P, u = -Kx", "m x n matrix",
		{Mat("A", "1/s"), Mat("B", "-"), Mat("Q", "-"), Mat("R", "-")}, &Lqr);
	Def(R, "kalman", "control", "approach/ section 2: the steady-state Kalman filter", "steady-state measurement-update gain M = P C'(C P C' + R)^-1", "n x p matrix",
		{Mat("A", "-"), Mat("C", "-"), Mat("Q", "-"), Mat("R", "-")}, &KalmanFn);
	Def(R, "kfcov", "control", "approach/ section 2: the steady-state Kalman filter", "steady-state prior covariance P of the filter (the dual DARE)", "n x n matrix",
		{Mat("A", "-"), Mat("C", "-"), Mat("Q", "-"), Mat("R", "-")}, &KfCov);
	Def(R, "c2d", "control", "approach/ discretize(): zero-order hold", "zero-order-hold discretization, expm of [A B; 0 0] dt", "[Ad, Bd] side by side, n x (n+m)",
		{Mat("A", "1/s"), Mat("B", "-"), Num("dt", "s")}, &C2d);
	Def(R, "c2dq", "control", "approach/ discretize(): Van Loan", "discrete process noise Qd from continuous intensity W by Van Loan's method", "n x n matrix",
		{Mat("A", "1/s"), Mat("W", "-"), Num("dt", "s")}, &C2dq);
	Def(R, "mpcgain", "control", "mpc/ section 1: unconstrained MPC = LQR", "first-move gain of unconstrained MPC, horizon N, terminal weight Pf (Pf = dare(...) gives dlqr)", "m x n matrix",
		{Mat("A", "-"), Mat("B", "-"), Mat("Q", "-"), Mat("R", "-"), Mat("Pf", "-"), Num("N", "steps")}, &MpcGain);
	Def(R, "clcost", "control", "mpc/ section 4 and approach/ lqrCost()", "infinite-horizon cost of u = -Kx from x0, x0' X x0 with X = F'XF + Q + K'RK", "cost",
		{Mat("A", "-"), Mat("B", "-"), Mat("Q", "-"), Mat("R", "-"), Mat("K", "-"), Mat("x0", "-")}, &ClCost);
	Def(R, "gm", "control", "approach/ section 5: gain margins", "interval of input gains g around 1 with A - gBK stable (discrete)", "record lo, hi",
		{Mat("A", "-"), Mat("B", "-"), Mat("K", "-"), Num("lo", "-", 1e-3), Num("hi", "-", 1e3)}, &Gm);
	Def(R, "lqrgm", "control", "approach/ section 5: return-difference guarantee", "guaranteed discrete LQR gain margin 1/(1+rho)..1/(1-rho), rho = sqrt(R/(R + B'PB)) (SISO)", "record rho, lo, hi",
		{Mat("A", "-"), Mat("B", "-"), Mat("Q", "-"), Mat("R", "-")}, &LqrGm);
	Def(R, "retdiff", "control", "approach/ section 5: return-difference identity", "min over the unit circle of |1 + K (zI - A)^-1 B| (SISO, 600 frequencies)", "-",
		{Mat("A", "-"), Mat("B", "-"), Mat("K", "-")}, &RetDiff);
	Def(R, "lqgcl", "control", "approach/ section 3: the separation principle", "LQG closed loop in (x, prior estimate) with gain g at the plant input; eig = eig(A-BK) and eig(A-AMC)", "2n x 2n matrix",
		{Mat("A", "-"), Mat("B", "-"), Mat("C", "-"), Mat("K", "-"), Mat("M", "-"), Num("g", "-", 1)}, &LqgCl);
	Def(R, "doyle", "control", "approach/ section 5: Doyle's 1978 example", "Doyle's LQG example: gains f, d and the stable input-gain interval of LQG and of LQR", "record f, d, lo, hi (LQG), rlo, rhi (LQR)",
		{Num("q", "-", 100), Num("sigma", "- (nan: = q)", CtlNan)}, &Doyle);
	Def(R, "glide", "control", "approach/ sections 4 and 5: dispersion and margins", "glideslope approach: stationary dispersion and gain interval of 1 LQR, 2 LQG, 3 PID", "record sd (m), J, su (rad), rho, lo, hi, rd, rq, glo, ghi",
		{Num("mode", "-", 2), Num("sg", "m/s", 1.5), Num("tau", "s", 4), Num("sd", "m", 2), Num("sv", "m/s", 0.3)}, &Glide);
	Def(R, "lanesat", "control", "mpc/ section 2: saturated LQR lane change", "lane change by LQR clipped to |u| <= umax (150 steps of 0.2 s)", "record J, os (overshoot, m)",
		{Num("umax", "m/s^2", 1), Num("y0", "m", 3.5), Num("r", "-", 0.01)}, &LaneSat);
	Def(R, "lanempc", "control", "mpc/ sections 2 and 4: constrained MPC", "lane change by constrained MPC (dual active-set QP each step), term 1 = terminal cost P, 0 = none", "record J (closed loop), os (m), V (open-loop V_N)",
		{Num("N", "steps", 20), Num("umax", "m/s^2", 1), Num("y0", "m", 3.5), Num("r", "-", 0.01), Num("term", "-", 1)}, &LaneMpc);
	Def(R, "lander", "control", "mpc/ section 3: receding horizon and terminal sets", "lander MPC with altitude >= 0; mode 0 none, 1 P, 2 P + v_N >= 0, 3 x_N = 0", "record fail (step, -1 = none), alt (m), vel (m/s), J",
		{Num("N", "steps", 22), Num("h0", "m", 40), Num("v0", "m/s (down)", 8), Num("uh", "m/s^2", 2), Num("mode", "-", 0)}, &Lander);

	// dynamics
	Def(R, "lyap", "dynamics", "ch08 Theorem 3: the Lyapunov equation", "continuous Lyapunov A X + X A' + Q = 0 (use lyap(A', Q) for A'P + PA = -Q)", "n x n matrix",
		{Mat("A", "1/s"), Mat("Q", "-")}, &LyapFn);
	Def(R, "dlyap", "dynamics", "ch08 Theorem 3 (discrete) and approach/ lyap()", "discrete Lyapunov X = A X A' + Q", "n x n matrix",
		{Mat("A", "-"), Mat("Q", "-")}, &DlyapFn);
	Def(R, "roa", "dynamics", "ch08 bench: a certificate for a basin", "reversed Van der Pol: certified level c* of V = x'Px (Q = diag(1, q)), its area and the true basin", "record c, area, basin, frac",
		{Num("q", "-", 1)}, &Roa);
	Def(R, "nomoto", "dynamics", "ch08 bench: an autopilot that overshoots", "Nomoto ship heading change: 1 PID, 2 PID + anti-windup, 3 sliding mode, 4 sign law", "record over (rad), settle (s), sw (rudder reversals), reach (s), bound (s)",
		{Num("ctrl", "-", 1), Num("step", "rad", 60 * cst::Pi / 180), Num("wave", "rad", 3 * cst::Pi / 180)}, &Nomoto);
	Def(R, "logmap", "dynamics", "ch08 counterexample: the logistic map", "n-th iterate of x -> r x (1 - x)", "x_n",
		{Num("r", "-"), Num("x0", "-"), Num("n", "steps")}, &LogMap);
	Def(R, "lyapexp", "dynamics", "ch08 counterexample: the logistic map", "Lyapunov exponent of the logistic map from x0 = 0.3 (500 steps burn-in)", "per step",
		{Num("r", "-", 3.83), Num("n", "steps", 20000)}, &LyapExp);
	Def(R, "logperiod", "dynamics", "ch08 counterexample: the logistic map", "period of the logistic attractor (0 if none up to 64)", "steps",
		{Num("r", "-", 3.83)}, &LogPeriod);
	Def(R, "superstab", "dynamics", "ch08 Feigenbaum ratios", "superstable parameter R_n of the logistic map (period 2^n)", "r",
		{Num("n", "-")}, &Superstab);
	Def(R, "feigen", "dynamics", "ch08 Feigenbaum ratios", "Feigenbaum ratio (R_{n-1} - R_{n-2})/(R_n - R_{n-1}) -> 4.6692", "-",
		{Num("n", "-", 8)}, &Feigen);
	Def(R, "pendulum", "dynamics", "ch08 Theorem 2: LaSalle, the damped pendulum", "damped pendulum th'' = -sin th - b th' from rest at th0, RK4 at 0.01 s", "record th (rad), w (rad/s), E",
		{Num("th0", "rad", 1.5), Num("b", "1/s", 0.3), Num("t", "s", 60)}, &Pendulum);
}
} // namespace mcalc
