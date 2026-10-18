// Dense row-major matrix of doubles and the linear algebra the function library needs.
// Engine-independent C++17: no Unreal headers, no exceptions, no RTTI. Functions that can fail
// return false and leave a reason in *Error; everything else is total.
#pragma once

#include <string>
#include <vector>

namespace mcalc
{
struct Matrix
{
	int Rows = 0;
	int Cols = 0;
	std::vector<double> V;

	Matrix() = default;
	Matrix(int InRows, int InCols, double Fill = 0.0) : Rows(InRows), Cols(InCols), V(static_cast<size_t>(InRows) * static_cast<size_t>(InCols), Fill) {}

	static Matrix Scalar(double X) { return Matrix(1, 1, X); }
	static Matrix Identity(int N);
	static Matrix Column(const std::vector<double>& Values);
	static Matrix Row(const std::vector<double>& Values);
	static Matrix FromRows(const std::vector<std::vector<double>>& RowValues);

	double& operator()(int R, int C) { return V[static_cast<size_t>(R) * static_cast<size_t>(Cols) + static_cast<size_t>(C)]; }
	double operator()(int R, int C) const { return V[static_cast<size_t>(R) * static_cast<size_t>(Cols) + static_cast<size_t>(C)]; }

	int Count() const { return Rows * Cols; }
	bool IsScalar() const { return Rows == 1 && Cols == 1; }
	bool IsVector() const { return (Rows == 1 || Cols == 1) && Count() > 0; }
	bool IsSquare() const { return Rows == Cols && Rows > 0; }
	bool IsEmpty() const { return Count() == 0; }
	double ScalarValue() const { return V.empty() ? 0.0 : V[0]; }
	std::vector<double> ToVector() const { return V; }
};

namespace la
{
Matrix Transpose(const Matrix& A);
Matrix Add(const Matrix& A, const Matrix& B);   // shapes must match (callers check)
Matrix Sub(const Matrix& A, const Matrix& B);
Matrix Scale(const Matrix& A, double S);
Matrix Mul(const Matrix& A, const Matrix& B);   // A.Cols == B.Rows (callers check)
double Trace(const Matrix& A);
double FrobeniusNorm(const Matrix& A);
double MaxAbs(const Matrix& A);
Matrix Symmetrize(const Matrix& A);             // (A + A^T)/2

// LU with partial pivoting. Returns false for a singular (pivot below tolerance) matrix.
bool Determinant(const Matrix& A, double* Out, std::string* Error);
bool Solve(const Matrix& A, const Matrix& B, Matrix* X, std::string* Error);     // A X = B, A square
bool Inverse(const Matrix& A, Matrix* Out, std::string* Error);
bool LeastSquares(const Matrix& A, const Matrix& B, Matrix* X, std::string* Error); // min |A X - B| by QR
bool Cholesky(const Matrix& A, Matrix* L, std::string* Error);                   // A = L L^T, lower L

// Symmetric eigen-decomposition by cyclic Jacobi: A = V diag(W) V^T, W ascending.
bool SymmetricEigen(const Matrix& A, std::vector<double>* W, Matrix* Vecs, std::string* Error);
// General real eigenvalues by Hessenberg reduction and the shifted QR algorithm.
// Returns real and imaginary parts sorted by descending real part, then descending imaginary part.
bool Eigenvalues(const Matrix& A, std::vector<double>* Re, std::vector<double>* Im, std::string* Error);
// Thin SVD by one-sided Jacobi: A = U diag(S) V^T, S descending.
bool Svd(const Matrix& A, Matrix* U, std::vector<double>* S, Matrix* Vt, std::string* Error);
// Matrix exponential by scaling and squaring with a degree-13 Pade approximant (Higham 2005).
bool Expm(const Matrix& A, Matrix* Out, std::string* Error);
double SpectralRadius(const Matrix& A);         // NaN if the eigenvalue iteration fails
int Rank(const Matrix& A, double Tol = -1.0);
Matrix Kron(const Matrix& A, const Matrix& B);

// Discrete algebraic Riccati equation P = A'PA - A'PB (R + B'PB)^-1 B'PA + Q by structure-preserving
// doubling, cross-checked against plain iteration in the tests. Gain K = (R + B'PB)^-1 B'PA.
bool Dare(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, Matrix* P, Matrix* K, std::string* Error);
// Continuous algebraic Riccati equation A'P + PA - PBR^-1B'P + Q = 0 via the matrix sign function of the Hamiltonian.
bool Care(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, Matrix* P, Matrix* K, std::string* Error);
// Discrete Lyapunov X = A X A' + Q (Kronecker solve, n <= 12) and continuous A X + X A' + Q = 0.
bool Dlyap(const Matrix& A, const Matrix& Q, Matrix* X, std::string* Error);
bool Lyap(const Matrix& A, const Matrix& Q, Matrix* X, std::string* Error);
} // namespace la
} // namespace mcalc
