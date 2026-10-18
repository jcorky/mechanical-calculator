// Elementary functions, statistics, polynomials and matrix operations: the calculator's base layer.
#include "CalcLibrary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace mcalc
{
namespace builtins_detail
{
template <double (*F)(double)>
Result Map1(const Call& C)
{
	Matrix M = C.M(0);
	for (double& X : M.V)
	{
		X = F(X);
	}
	return Result::Of(M);
}

double FSin(double X) { return std::sin(X); }
double FCos(double X) { return std::cos(X); }
double FTan(double X) { return std::tan(X); }
double FAsin(double X) { return std::asin(X); }
double FAcos(double X) { return std::acos(X); }
double FAtan(double X) { return std::atan(X); }
double FSinh(double X) { return std::sinh(X); }
double FCosh(double X) { return std::cosh(X); }
double FTanh(double X) { return std::tanh(X); }
double FAsinh(double X) { return std::asinh(X); }
double FAcosh(double X) { return std::acosh(X); }
double FAtanh(double X) { return std::atanh(X); }
double FExp(double X) { return std::exp(X); }
double FLn(double X) { return std::log(X); }
double FLog10(double X) { return std::log10(X); }
double FLog2(double X) { return std::log2(X); }
double FSqrt(double X) { return std::sqrt(X); }
double FCbrt(double X) { return std::cbrt(X); }
double FAbs(double X) { return std::fabs(X); }
double FFloor(double X) { return std::floor(X); }
double FCeil(double X) { return std::ceil(X); }
double FRound(double X) { return std::round(X); }
double FTrunc(double X) { return std::trunc(X); }
double FSign(double X) { return X > 0 ? 1.0 : (X < 0 ? -1.0 : 0.0); }
double FGamma(double X) { return std::tgamma(X); }
double FLgamma(double X) { return std::lgamma(X); }
double FErf(double X) { return std::erf(X); }
double FErfc(double X) { return std::erfc(X); }
double FDeg(double X) { return X * 180.0 / cst::Pi; }
double FRad(double X) { return X * cst::Pi / 180.0; }
double FDb(double X) { return 10.0 * std::log10(X); }
double FUndb(double X) { return std::pow(10.0, X / 10.0); }
double FFact(double X) { return std::tgamma(X + 1.0); }
// Standard normal CDF and its inverse (Acklam's rational approximation refined by one Halley step).
double FNormcdf(double X) { return 0.5 * std::erfc(-X / std::sqrt(2.0)); }
double FNorminv(double P)
{
	if (P <= 0.0)
	{
		return -std::numeric_limits<double>::infinity();
	}
	if (P >= 1.0)
	{
		return std::numeric_limits<double>::infinity();
	}
	static const double A[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
	static const double B[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
	static const double Cc[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
	static const double D[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
	double X = 0.0;
	if (P < 0.02425)
	{
		const double Q = std::sqrt(-2.0 * std::log(P));
		X = (((((Cc[0] * Q + Cc[1]) * Q + Cc[2]) * Q + Cc[3]) * Q + Cc[4]) * Q + Cc[5]) / ((((D[0] * Q + D[1]) * Q + D[2]) * Q + D[3]) * Q + 1.0);
	}
	else if (P > 1.0 - 0.02425)
	{
		const double Q = std::sqrt(-2.0 * std::log(1.0 - P));
		X = -(((((Cc[0] * Q + Cc[1]) * Q + Cc[2]) * Q + Cc[3]) * Q + Cc[4]) * Q + Cc[5]) / ((((D[0] * Q + D[1]) * Q + D[2]) * Q + D[3]) * Q + 1.0);
	}
	else
	{
		const double Q = P - 0.5, R = Q * Q;
		X = (((((A[0] * R + A[1]) * R + A[2]) * R + A[3]) * R + A[4]) * R + A[5]) * Q / (((((B[0] * R + B[1]) * R + B[2]) * R + B[3]) * R + B[4]) * R + 1.0);
	}
	const double E = FNormcdf(X) - P;
	const double U = E * std::sqrt(2.0 * cst::Pi) * std::exp(X * X / 2.0);
	return X - U / (1.0 + X * U / 2.0);
}

Result Atan2(const Call& C) { return Result::Of(std::atan2(C.N(0), C.N(1))); }
Result Hypot(const Call& C) { return Result::Of(std::hypot(C.N(0), C.N(1))); }
Result Mod(const Call& C)
{
	const double A = C.N(0), B = C.N(1);
	if (B == 0.0)
	{
		return Result::Of(A);
	}
	return Result::Of(A - B * std::floor(A / B));
}
Result Rem(const Call& C) { return Result::Of(std::fmod(C.N(0), C.N(1))); }
Result Log(const Call& C) { return Result::Of(std::log(C.N(0)) / std::log(C.N(1))); }
Result Root(const Call& C)
{
	const double X = C.N(0), N = C.N(1);
	if (X < 0 && std::fmod(N, 2.0) == 1.0)
	{
		return Result::Of(-std::pow(-X, 1.0 / N));
	}
	return Result::Of(std::pow(X, 1.0 / N));
}

Result NChooseK(const Call& C)
{
	const double N = C.N(0), K = C.N(1);
	if (K < 0 || K > N)
	{
		return Result::Of(0.0);
	}
	return Result::Of(std::round(std::exp(std::lgamma(N + 1) - std::lgamma(K + 1) - std::lgamma(N - K + 1))));
}
Result NPermK(const Call& C)
{
	const double N = C.N(0), K = C.N(1);
	if (K < 0 || K > N)
	{
		return Result::Of(0.0);
	}
	return Result::Of(std::round(std::exp(std::lgamma(N + 1) - std::lgamma(N - K + 1))));
}
Result Gcd(const Call& C)
{
	long long A = std::llabs(std::llround(C.N(0))), B = std::llabs(std::llround(C.N(1)));
	while (B != 0)
	{
		const long long T = A % B;
		A = B;
		B = T;
	}
	return Result::Of(static_cast<double>(A));
}
Result Lcm(const Call& C)
{
	const long long A = std::llabs(std::llround(C.N(0))), B = std::llabs(std::llround(C.N(1)));
	if (A == 0 || B == 0)
	{
		return Result::Of(0.0);
	}
	long long X = A, Y = B;
	while (Y != 0)
	{
		const long long T = X % Y;
		X = Y;
		Y = T;
	}
	return Result::Of(static_cast<double>(A / X * B));
}

// Variadic min/max: min(v) over the elements, or min(a, b, ...) over scalars.
std::vector<double> Flatten(const Call& C)
{
	std::vector<double> All;
	for (size_t i = 0; i < C.Count(); ++i)
	{
		All.insert(All.end(), C.M(i).V.begin(), C.M(i).V.end());
	}
	return All;
}
Result Min(const Call& C)
{
	const std::vector<double> A = Flatten(C);
	if (A.empty())
	{
		return Result::Fail("needs at least one value");
	}
	return Result::Of(*std::min_element(A.begin(), A.end()));
}
Result Max(const Call& C)
{
	const std::vector<double> A = Flatten(C);
	if (A.empty())
	{
		return Result::Fail("needs at least one value");
	}
	return Result::Of(*std::max_element(A.begin(), A.end()));
}
Result Sum(const Call& C)
{
	const std::vector<double>& A = C.M(0).V;
	return Result::Of(std::accumulate(A.begin(), A.end(), 0.0));
}
Result Prod(const Call& C)
{
	double P = 1.0;
	for (double X : C.M(0).V)
	{
		P *= X;
	}
	return Result::Of(P);
}
Result Mean(const Call& C)
{
	const std::vector<double>& A = C.M(0).V;
	if (A.empty())
	{
		return Result::Fail("mean of nothing");
	}
	return Result::Of(std::accumulate(A.begin(), A.end(), 0.0) / static_cast<double>(A.size()));
}
double Variance(const std::vector<double>& A)
{
	if (A.size() < 2)
	{
		return 0.0;
	}
	const double M = std::accumulate(A.begin(), A.end(), 0.0) / static_cast<double>(A.size());
	double S = 0.0;
	for (double X : A)
	{
		S += (X - M) * (X - M);
	}
	return S / static_cast<double>(A.size() - 1);
}
Result Var(const Call& C) { return Result::Of(Variance(C.M(0).V)); }
Result Std(const Call& C) { return Result::Of(std::sqrt(Variance(C.M(0).V))); }
Result Median(const Call& C)
{
	std::vector<double> A = C.M(0).V;
	if (A.empty())
	{
		return Result::Fail("median of nothing");
	}
	std::sort(A.begin(), A.end());
	const size_t N = A.size();
	return Result::Of(N % 2 ? A[N / 2] : 0.5 * (A[N / 2 - 1] + A[N / 2]));
}
Result Sort(const Call& C)
{
	Matrix M = C.M(0);
	std::sort(M.V.begin(), M.V.end());
	return Result::Of(M);
}
Result Cumsum(const Call& C)
{
	Matrix M = C.M(0);
	double S = 0.0;
	for (double& X : M.V)
	{
		S += X;
		X = S;
	}
	return Result::Of(M);
}
Result Numel(const Call& C) { return Result::Of(static_cast<double>(C.M(0).Count())); }
Result Size(const Call& C) { return Result::Record({{"rows", static_cast<double>(C.M(0).Rows)}, {"cols", static_cast<double>(C.M(0).Cols)}}); }

bool Dim(double X, int* Out)
{
	if (!(X >= 0 && X <= 4096 && X == std::floor(X)))
	{
		return false;
	}
	*Out = static_cast<int>(X);
	return true;
}
Result Fill(const Call& C, double V)
{
	int R = 0, Cols = 0;
	if (!Dim(C.N(0), &R) || !Dim(C.N(1) < 0 ? C.N(0) : C.N(1), &Cols))
	{
		return Result::Fail("sizes must be whole numbers 0..4096");
	}
	return Result::Of(Matrix(R, Cols, V));
}
Result Zeros(const Call& C) { return Fill(C, 0.0); }
Result Ones(const Call& C) { return Fill(C, 1.0); }
Result Eye(const Call& C)
{
	int N = 0;
	if (!Dim(C.N(0), &N))
	{
		return Result::Fail("size must be a whole number 0..4096");
	}
	return Result::Of(Matrix::Identity(N));
}
Result Linspace(const Call& C)
{
	int N = 0;
	if (!Dim(C.N(2), &N) || N < 2)
	{
		return Result::Fail("n must be a whole number >= 2");
	}
	Matrix M(1, N);
	for (int i = 0; i < N; ++i)
	{
		M.V[static_cast<size_t>(i)] = C.N(0) + (C.N(1) - C.N(0)) * i / (N - 1);
	}
	return Result::Of(M);
}
Result Diag(const Call& C)
{
	const Matrix& A = C.M(0);
	if (A.IsVector())
	{
		const int N = A.Count();
		Matrix D(N, N);
		for (int i = 0; i < N; ++i)
		{
			D(i, i) = A.V[static_cast<size_t>(i)];
		}
		return Result::Of(D);
	}
	const int N = std::min(A.Rows, A.Cols);
	Matrix D(N, 1);
	for (int i = 0; i < N; ++i)
	{
		D.V[static_cast<size_t>(i)] = A(i, i);
	}
	return Result::Of(D);
}
Result Reshape(const Call& C)
{
	int R = 0, Cols = 0;
	if (!Dim(C.N(1), &R) || !Dim(C.N(2), &Cols) || R * Cols != C.M(0).Count())
	{
		return Result::Fail("rows * cols must equal the element count");
	}
	Matrix M(R, Cols);
	M.V = C.M(0).V;
	return Result::Of(M);
}
Result Trace(const Call& C) { return Result::Of(la::Trace(C.M(0))); }
Result Det(const Call& C)
{
	double D = 0.0;
	std::string E;
	if (!la::Determinant(C.M(0), &D, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(D);
}
Result Inv(const Call& C)
{
	Matrix X;
	std::string E;
	if (!la::Inverse(C.M(0), &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}
Result Pinv(const Call& C)
{
	Matrix X;
	std::string E;
	if (!la::LeastSquares(C.M(0), Matrix::Identity(C.M(0).Rows), &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}
Result Trans(const Call& C) { return Result::Of(la::Transpose(C.M(0))); }
Result Rank(const Call& C) { return Result::Of(static_cast<double>(la::Rank(C.M(0)))); }
Result Norm(const Call& C)
{
	const Matrix& A = C.M(0);
	const double P = C.N(1);
	if (A.IsVector() || A.IsScalar())
	{
		if (std::isinf(P))
		{
			return Result::Of(la::MaxAbs(A));
		}
		double S = 0.0;
		for (double X : A.V)
		{
			S += std::pow(std::fabs(X), P);
		}
		return Result::Of(std::pow(S, 1.0 / P));
	}
	if (P == 2.0)
	{
		std::vector<double> S;
		std::string E;
		if (!la::Svd(A, nullptr, &S, nullptr, &E))
		{
			return Result::Fail(E);
		}
		return Result::Of(S.empty() ? 0.0 : S[0]);
	}
	if (P == 1.0 || std::isinf(P))
	{
		const Matrix B = P == 1.0 ? A : la::Transpose(A);
		double Best = 0.0;
		for (int j = 0; j < B.Cols; ++j)
		{
			double S = 0.0;
			for (int i = 0; i < B.Rows; ++i)
			{
				S += std::fabs(B(i, j));
			}
			Best = std::max(Best, S);
		}
		return Result::Of(Best);
	}
	if (P == 0.0)
	{
		return Result::Of(la::FrobeniusNorm(A));
	}
	return Result::Fail("matrix norms: p = 1, 2, inf, or 0 for Frobenius");
}
Result Cond(const Call& C)
{
	std::vector<double> S;
	std::string E;
	if (!la::Svd(C.M(0), nullptr, &S, nullptr, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(S.back() > 0 ? S.front() / S.back() : std::numeric_limits<double>::infinity());
}
Result Dot(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1);
	if (A.Count() != B.Count())
	{
		return Result::Fail("vectors must have the same length");
	}
	double S = 0.0;
	for (size_t i = 0; i < A.V.size(); ++i)
	{
		S += A.V[i] * B.V[i];
	}
	return Result::Of(S);
}
Result Cross(const Call& C)
{
	const Matrix &A = C.M(0), &B = C.M(1);
	if (A.Count() != 3 || B.Count() != 3)
	{
		return Result::Fail("cross needs two 3-vectors");
	}
	const std::vector<double>&a = A.V, &b = B.V;
	Matrix R = A;
	R.V = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
	return Result::Of(R);
}
Result Kron(const Call& C) { return Result::Of(la::Kron(C.M(0), C.M(1))); }
Result SolveFn(const Call& C)
{
	Matrix X;
	std::string E;
	const bool Ok = C.M(0).IsSquare() ? la::Solve(C.M(0), C.M(1), &X, &E) : la::LeastSquares(C.M(0), C.M(1), &X, &E);
	if (!Ok)
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}
Result Lsq(const Call& C)
{
	Matrix X;
	std::string E;
	if (!la::LeastSquares(C.M(0), C.M(1), &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}
Result Chol(const Call& C)
{
	Matrix L;
	std::string E;
	if (!la::Cholesky(C.M(0), &L, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(L);
}
bool IsSymmetric(const Matrix& A)
{
	if (!A.IsSquare())
	{
		return false;
	}
	const double Tol = 1e-12 * std::max(1.0, la::MaxAbs(A));
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int j = i + 1; j < A.Cols; ++j)
		{
			if (std::fabs(A(i, j) - A(j, i)) > Tol)
			{
				return false;
			}
		}
	}
	return true;
}
Result Eig(const Call& C)
{
	const Matrix& A = C.M(0);
	std::string E;
	if (IsSymmetric(A))
	{
		std::vector<double> W;
		if (!la::SymmetricEigen(A, &W, nullptr, &E))
		{
			return Result::Fail(E);
		}
		return Result::Of(Matrix::Column(W));
	}
	std::vector<double> Re, Im;
	if (!la::Eigenvalues(A, &Re, &Im, &E))
	{
		return Result::Fail(E);
	}
	bool Real = true;
	for (double X : Im)
	{
		Real = Real && X == 0.0;
	}
	if (Real)
	{
		return Result::Of(Matrix::Column(Re));
	}
	Matrix M(static_cast<int>(Re.size()), 2);
	for (size_t i = 0; i < Re.size(); ++i)
	{
		M(static_cast<int>(i), 0) = Re[i];
		M(static_cast<int>(i), 1) = Im[i];
	}
	return Result::Of(M);
}
Result Eigvec(const Call& C)
{
	std::vector<double> W;
	Matrix V;
	std::string E;
	if (!IsSymmetric(C.M(0)))
	{
		return Result::Fail("eigvec is for symmetric matrices (columns match eig's ascending order)");
	}
	if (!la::SymmetricEigen(C.M(0), &W, &V, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(V);
}
Result SvdFn(const Call& C)
{
	std::vector<double> S;
	std::string E;
	if (!la::Svd(C.M(0), nullptr, &S, nullptr, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(Matrix::Column(S));
}
Result ExpmFn(const Call& C)
{
	Matrix X;
	std::string E;
	if (!la::Expm(C.M(0), &X, &E))
	{
		return Result::Fail(E);
	}
	return Result::Of(X);
}
Result Rho(const Call& C)
{
	if (!C.M(0).IsSquare())
	{
		return Result::Fail("needs a square matrix");
	}
	return Result::Of(la::SpectralRadius(C.M(0)));
}

// Polynomials, highest power first (MATLAB convention).
Result Polyval(const Call& C)
{
	const std::vector<double>& P = C.M(0).V;
	Matrix X = C.M(1);
	for (double& V : X.V)
	{
		double S = 0.0;
		for (double Coef : P)
		{
			S = S * V + Coef;
		}
		V = S;
	}
	return Result::Of(X);
}
Result Roots(const Call& C)
{
	std::vector<double> P = C.M(0).V;
	while (!P.empty() && P.front() == 0.0)
	{
		P.erase(P.begin());
	}
	int Zeros = 0;
	while (!P.empty() && P.back() == 0.0)
	{
		P.pop_back();
		++Zeros;
	}
	const int N = static_cast<int>(P.size()) - 1;
	std::vector<double> Re(static_cast<size_t>(Zeros), 0.0), Im(static_cast<size_t>(Zeros), 0.0);
	if (N >= 1)
	{
		Matrix Comp(N, N);
		for (int j = 0; j < N; ++j)
		{
			Comp(0, j) = -P[static_cast<size_t>(j + 1)] / P[0];
		}
		for (int i = 1; i < N; ++i)
		{
			Comp(i, i - 1) = 1.0;
		}
		std::vector<double> R2, I2;
		std::string E;
		if (!la::Eigenvalues(Comp, &R2, &I2, &E))
		{
			return Result::Fail(E);
		}
		Re.insert(Re.end(), R2.begin(), R2.end());
		Im.insert(Im.end(), I2.begin(), I2.end());
	}
	bool Real = true;
	for (double X : Im)
	{
		Real = Real && X == 0.0;
	}
	if (Real)
	{
		return Result::Of(Matrix::Column(Re));
	}
	Matrix M(static_cast<int>(Re.size()), 2);
	for (size_t i = 0; i < Re.size(); ++i)
	{
		M(static_cast<int>(i), 0) = Re[i];
		M(static_cast<int>(i), 1) = Im[i];
	}
	return Result::Of(M);
}
Result Conv(const Call& C)
{
	const std::vector<double> &A = C.M(0).V, &B = C.M(1).V;
	if (A.empty() || B.empty())
	{
		return Result::Of(Matrix());
	}
	std::vector<double> Out(A.size() + B.size() - 1, 0.0);
	for (size_t i = 0; i < A.size(); ++i)
	{
		for (size_t j = 0; j < B.size(); ++j)
		{
			Out[i + j] += A[i] * B[j];
		}
	}
	return Result::Of(Matrix::Row(Out));
}
Result Normpdf(const Call& C)
{
	Matrix X = C.M(0);
	for (double& V : X.V)
	{
		V = std::exp(-0.5 * V * V) / std::sqrt(2.0 * cst::Pi);
	}
	return Result::Of(X);
}
} // namespace builtins_detail

void RegisterCoreMath(Registry& R)
{
	using namespace builtins_detail;
	const char* Src = "calculator";
	struct Unary
	{
		const char* Name;
		const char* Summary;
		NativeFn Fn;
		int Angle;
	};
	const Unary Elementary[] = {
		{"sin", "sine (follows DEG/RAD)", &Map1<FSin>, AngleArgs},
		{"cos", "cosine (follows DEG/RAD)", &Map1<FCos>, AngleArgs},
		{"tan", "tangent (follows DEG/RAD)", &Map1<FTan>, AngleArgs},
		{"asin", "inverse sine (follows DEG/RAD)", &Map1<FAsin>, AngleResult},
		{"acos", "inverse cosine (follows DEG/RAD)", &Map1<FAcos>, AngleResult},
		{"atan", "inverse tangent (follows DEG/RAD)", &Map1<FAtan>, AngleResult},
		{"sinh", "hyperbolic sine", &Map1<FSinh>, AngleNone},
		{"cosh", "hyperbolic cosine", &Map1<FCosh>, AngleNone},
		{"tanh", "hyperbolic tangent", &Map1<FTanh>, AngleNone},
		{"asinh", "inverse hyperbolic sine", &Map1<FAsinh>, AngleNone},
		{"acosh", "inverse hyperbolic cosine", &Map1<FAcosh>, AngleNone},
		{"atanh", "inverse hyperbolic tangent", &Map1<FAtanh>, AngleNone},
		{"exp", "e^x", &Map1<FExp>, AngleNone},
		{"ln", "natural logarithm", &Map1<FLn>, AngleNone},
		{"log10", "base-10 logarithm", &Map1<FLog10>, AngleNone},
		{"log2", "base-2 logarithm", &Map1<FLog2>, AngleNone},
		{"sqrt", "square root", &Map1<FSqrt>, AngleNone},
		{"cbrt", "cube root", &Map1<FCbrt>, AngleNone},
		{"abs", "absolute value", &Map1<FAbs>, AngleNone},
		{"floor", "round down", &Map1<FFloor>, AngleNone},
		{"ceil", "round up", &Map1<FCeil>, AngleNone},
		{"round", "round half away from zero", &Map1<FRound>, AngleNone},
		{"trunc", "round toward zero", &Map1<FTrunc>, AngleNone},
		{"sign", "-1, 0 or 1", &Map1<FSign>, AngleNone},
		{"gamma", "gamma function", &Map1<FGamma>, AngleNone},
		{"lgamma", "log of |gamma|", &Map1<FLgamma>, AngleNone},
		{"fact", "factorial, gamma(x + 1) (also postfix !)", &Map1<FFact>, AngleNone},
		{"erf", "error function", &Map1<FErf>, AngleNone},
		{"erfc", "complementary error function", &Map1<FErfc>, AngleNone},
		{"deg", "radians to degrees", &Map1<FDeg>, AngleNone},
		{"rad", "degrees to radians", &Map1<FRad>, AngleNone},
		{"db", "power ratio to decibels, 10 log10 x", &Map1<FDb>, AngleNone},
		{"undb", "decibels to power ratio", &Map1<FUndb>, AngleNone},
		{"normcdf", "standard normal CDF", &Map1<FNormcdf>, AngleNone},
		{"norminv", "standard normal quantile", &Map1<FNorminv>, AngleNone},
		{"normpdf", "standard normal density", &Normpdf, AngleNone},
	};
	for (const Unary& U : Elementary)
	{
		FunctionInfo F;
		F.Name = U.Name;
		F.Category = "math";
		F.Source = Src;
		F.Summary = U.Summary;
		F.Returns = "element-wise";
		F.Params = {Mat("x", "")};
		F.Fn = U.Fn;
		F.Angle = U.Angle;
		R.Add(F);
	}
	{
		FunctionInfo F;
		F.Name = "atan2";
		F.Category = "math";
		F.Source = Src;
		F.Summary = "four-quadrant inverse tangent of y/x (follows DEG/RAD)";
		F.Params = {Num("y", ""), Num("x", "")};
		F.Fn = &Atan2;
		F.Angle = AngleResult;
		R.Add(F);
	}
	Def(R, "hypot", "math", Src, "sqrt(x^2 + y^2) without overflow", "", {Num("x", ""), Num("y", "")}, &Hypot);
	Def(R, "mod", "math", Src, "floored remainder, sign of b", "", {Num("a", ""), Num("b", "")}, &Mod);
	Def(R, "rem", "math", Src, "truncated remainder, sign of a", "", {Num("a", ""), Num("b", "")}, &Rem);
	Def(R, "log", "math", Src, "logarithm of x to base b", "", {Num("x", ""), Num("b", "", 2.71828182845904523536)}, &Log);
	Def(R, "root", "math", Src, "real n-th root (odd roots of negatives allowed)", "", {Num("x", ""), Num("n", "")}, &Root);
	Def(R, "nCr", "math", Src, "binomial coefficient", "", {Num("n", ""), Num("k", "")}, &NChooseK);
	Def(R, "nPr", "math", Src, "permutations of k from n", "", {Num("n", ""), Num("k", "")}, &NPermK);
	Def(R, "gcd", "math", Src, "greatest common divisor", "", {Num("a", ""), Num("b", "")}, &Gcd);
	Def(R, "lcm", "math", Src, "least common multiple", "", {Num("a", ""), Num("b", "")}, &Lcm);

	auto Variadic = [&R, Src](const char* Name, const char* Summary, NativeFn Fn) {
		FunctionInfo F;
		F.Name = Name;
		F.Category = "stats";
		F.Source = Src;
		F.Summary = Summary;
		F.Params = {Mat("x", "")};
		F.Variadic = true;
		F.Fn = Fn;
		R.Add(F);
	};
	Variadic("min", "smallest element of all arguments", &Min);
	Variadic("max", "largest element of all arguments", &Max);
	Def(R, "sum", "stats", Src, "sum of elements", "", {Mat("x", "")}, &Sum);
	Def(R, "prod", "stats", Src, "product of elements", "", {Mat("x", "")}, &Prod);
	Def(R, "mean", "stats", Src, "arithmetic mean", "", {Mat("x", "")}, &Mean);
	Def(R, "median", "stats", Src, "median", "", {Mat("x", "")}, &Median);
	Def(R, "var", "stats", Src, "sample variance (n - 1)", "", {Mat("x", "")}, &Var);
	Def(R, "std", "stats", Src, "sample standard deviation (n - 1)", "", {Mat("x", "")}, &Std);
	Def(R, "sort", "stats", Src, "elements in ascending order", "", {Mat("x", "")}, &Sort);
	Def(R, "cumsum", "stats", Src, "running sum", "", {Mat("x", "")}, &Cumsum);

	Def(R, "zeros", "matrix", Src, "r x c zeros (c defaults to r)", "", {Num("r", ""), Num("c", "", -1)}, &Zeros);
	Def(R, "ones", "matrix", Src, "r x c ones (c defaults to r)", "", {Num("r", ""), Num("c", "", -1)}, &Ones);
	Def(R, "eye", "matrix", Src, "n x n identity", "", {Num("n", "")}, &Eye);
	Def(R, "linspace", "matrix", Src, "n evenly spaced points from a to b", "", {Num("a", ""), Num("b", ""), Num("n", "", 10)}, &Linspace);
	Def(R, "diag", "matrix", Src, "diagonal matrix from a vector, or the diagonal of a matrix", "", {Mat("x", "")}, &Diag);
	Def(R, "reshape", "matrix", Src, "same elements (row-major), new shape", "", {Mat("A", ""), Num("rows", ""), Num("cols", "")}, &Reshape);
	Def(R, "numel", "matrix", Src, "number of elements", "", {Mat("A", "")}, &Numel);
	Def(R, "size", "matrix", Src, "rows and columns", "record rows, cols", {Mat("A", "")}, &Size);
	Def(R, "trace", "matrix", Src, "sum of the diagonal", "", {Mat("A", "")}, &Trace);
	Def(R, "det", "matrix", Src, "determinant by LU", "", {Mat("A", "")}, &Det);
	Def(R, "inv", "matrix", Src, "inverse by LU", "", {Mat("A", "")}, &Inv);
	Def(R, "pinv", "matrix", Src, "Moore-Penrose pseudo-inverse by SVD", "", {Mat("A", "")}, &Pinv);
	Def(R, "trans", "matrix", Src, "transpose (also postfix ')", "", {Mat("A", "")}, &Trans);
	Def(R, "rank", "matrix", Src, "numerical rank by SVD", "", {Mat("A", "")}, &Rank);
	Def(R, "norm", "matrix", Src, "p-norm of a vector; matrix 1, 2, inf, or 0 = Frobenius", "", {Mat("x", ""), Num("p", "", 2)}, &Norm);
	Def(R, "cond", "matrix", Src, "2-norm condition number", "", {Mat("A", "")}, &Cond);
	Def(R, "dot", "matrix", Src, "inner product", "", {Mat("a", ""), Mat("b", "")}, &Dot);
	Def(R, "cross", "matrix", Src, "cross product of 3-vectors", "", {Mat("a", ""), Mat("b", "")}, &Cross);
	Def(R, "kron", "matrix", Src, "Kronecker product", "", {Mat("A", ""), Mat("B", "")}, &Kron);
	Def(R, "solve", "matrix", Src, "x with A x = b (least squares if A is not square)", "", {Mat("A", ""), Mat("b", "")}, &SolveFn);
	Def(R, "lsq", "matrix", Src, "least-squares solution of A x = b by SVD", "", {Mat("A", ""), Mat("b", "")}, &Lsq);
	Def(R, "chol", "matrix", Src, "lower Cholesky factor, A = L L'", "", {Mat("A", "")}, &Chol);
	Def(R, "eig", "matrix", Src, "eigenvalues: a column if real, else n x 2 [re, im]", "", {Mat("A", "")}, &Eig);
	Def(R, "eigvec", "matrix", Src, "eigenvectors of a symmetric matrix (columns)", "", {Mat("A", "")}, &Eigvec);
	Def(R, "svd", "matrix", Src, "singular values, descending", "", {Mat("A", "")}, &SvdFn);
	Def(R, "expm", "matrix", Src, "matrix exponential (Pade 13, scaling and squaring)", "", {Mat("A", "")}, &ExpmFn);
	Def(R, "rho", "matrix", Src, "spectral radius, max |eigenvalue|", "", {Mat("A", "")}, &Rho);
	Def(R, "polyval", "matrix", Src, "evaluate a polynomial (highest power first) at x", "", {Mat("p", ""), Mat("x", "")}, &Polyval);
	Def(R, "roots", "matrix", Src, "polynomial roots by companion-matrix eigenvalues", "", {Mat("p", "")}, &Roots);
	Def(R, "conv", "matrix", Src, "polynomial product", "", {Mat("p", ""), Mat("q", "")}, &Conv);
}
} // namespace mcalc
