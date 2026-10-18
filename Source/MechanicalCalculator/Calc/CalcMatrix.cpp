#include "CalcMatrix.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace mcalc
{
Matrix Matrix::Identity(int N)
{
	Matrix I(N, N);
	for (int i = 0; i < N; ++i)
	{
		I(i, i) = 1.0;
	}
	return I;
}

Matrix Matrix::Column(const std::vector<Real>& Values)
{
	Matrix M(static_cast<int>(Values.size()), 1);
	M.V = Values;
	return M;
}

Matrix Matrix::Row(const std::vector<Real>& Values)
{
	Matrix M(1, static_cast<int>(Values.size()));
	M.V = Values;
	return M;
}

Matrix Matrix::FromRows(const std::vector<std::vector<Real>>& RowValues)
{
	const int R = static_cast<int>(RowValues.size());
	const int C = R > 0 ? static_cast<int>(RowValues[0].size()) : 0;
	Matrix M(R, C);
	for (int i = 0; i < R; ++i)
	{
		for (int j = 0; j < C && j < static_cast<int>(RowValues[i].size()); ++j)
		{
			M(i, j) = RowValues[i][j];
		}
	}
	return M;
}

namespace la
{
namespace matrix_detail
{
Real Eps()
{
	return Real::Epsilon();
}

void SetError(std::string* Error, const char* Message)
{
	if (Error)
	{
		*Error = Message;
	}
}

Real Norm1(const Matrix& A)
{
	Real Best = 0.0;
	for (int j = 0; j < A.Cols; ++j)
	{
		Real S = 0.0;
		for (int i = 0; i < A.Rows; ++i)
		{
			S += gm::fabs(A(i, j));
		}
		Best = std::max(Best, S);
	}
	return Best;
}

// LU factorization with partial pivoting, in place. Perm[i] is the original row now at row i.
bool LuFactor(Matrix& A, std::vector<int>& Perm, int& Sign)
{
	const int N = A.Rows;
	Perm.resize(static_cast<size_t>(N));
	for (int i = 0; i < N; ++i)
	{
		Perm[static_cast<size_t>(i)] = i;
	}
	Sign = 1;
	const Real Scale = gm::max(MaxAbs(A), Real(1e-300));
	for (int c = 0; c < N; ++c)
	{
		int P = c;
		for (int r = c + 1; r < N; ++r)
		{
			if (gm::fabs(A(r, c)) > gm::fabs(A(P, c)))
			{
				P = r;
			}
		}
		if (gm::fabs(A(P, c)) <= Scale * 1e-14)
		{
			return false;
		}
		if (P != c)
		{
			for (int j = 0; j < N; ++j)
			{
				std::swap(A(P, j), A(c, j));
			}
			std::swap(Perm[static_cast<size_t>(P)], Perm[static_cast<size_t>(c)]);
			Sign = -Sign;
		}
		for (int r = c + 1; r < N; ++r)
		{
			const Real F = A(r, c) / A(c, c);
			A(r, c) = F;
			if (F != 0.0)
			{
				for (int j = c + 1; j < N; ++j)
				{
					A(r, j) -= F * A(c, j);
				}
			}
		}
	}
	return true;
}
} // namespace matrix_detail

using namespace matrix_detail;

Matrix Transpose(const Matrix& A)
{
	Matrix T(A.Cols, A.Rows);
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int j = 0; j < A.Cols; ++j)
		{
			T(j, i) = A(i, j);
		}
	}
	return T;
}

Matrix Add(const Matrix& A, const Matrix& B)
{
	Matrix C = A;
	for (size_t i = 0; i < C.V.size() && i < B.V.size(); ++i)
	{
		C.V[i] += B.V[i];
	}
	return C;
}

Matrix Sub(const Matrix& A, const Matrix& B)
{
	Matrix C = A;
	for (size_t i = 0; i < C.V.size() && i < B.V.size(); ++i)
	{
		C.V[i] -= B.V[i];
	}
	return C;
}

Matrix Scale(const Matrix& A, Real S)
{
	Matrix C = A;
	for (Real& X : C.V)
	{
		X *= S;
	}
	return C;
}

Matrix Mul(const Matrix& A, const Matrix& B)
{
	Matrix C(A.Rows, B.Cols);
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int k = 0; k < A.Cols; ++k)
		{
			const Real Aik = A(i, k);
			if (Aik == 0.0)
			{
				continue;
			}
			for (int j = 0; j < B.Cols; ++j)
			{
				C(i, j) += Aik * B(k, j);
			}
		}
	}
	return C;
}

Real Trace(const Matrix& A)
{
	Real S = 0.0;
	for (int i = 0; i < std::min(A.Rows, A.Cols); ++i)
	{
		S += A(i, i);
	}
	return S;
}

Real FrobeniusNorm(const Matrix& A)
{
	Real S = 0.0;
	for (Real X : A.V)
	{
		S += X * X;
	}
	return gm::sqrt(S);
}

Real MaxAbs(const Matrix& A)
{
	Real M = 0.0;
	for (Real X : A.V)
	{
		M = std::max(M, gm::fabs(X));
	}
	return M;
}

Matrix Symmetrize(const Matrix& A)
{
	return Scale(Add(A, Transpose(A)), 0.5);
}

bool Determinant(const Matrix& A, Real* Out, std::string* Error)
{
	if (!A.IsSquare())
	{
		SetError(Error, "det needs a square matrix");
		return false;
	}
	Matrix LU = A;
	std::vector<int> Perm;
	int Sign = 1;
	if (!LuFactor(LU, Perm, Sign))
	{
		*Out = 0.0;
		return true;
	}
	Real D = static_cast<Real>(Sign);
	for (int i = 0; i < LU.Rows; ++i)
	{
		D *= LU(i, i);
	}
	*Out = D;
	return true;
}

bool Solve(const Matrix& A, const Matrix& B, Matrix* X, std::string* Error)
{
	if (!A.IsSquare() || B.Rows != A.Rows)
	{
		SetError(Error, "solve needs square A and B with matching rows");
		return false;
	}
	Matrix LU = A;
	std::vector<int> Perm;
	int Sign = 1;
	if (!LuFactor(LU, Perm, Sign))
	{
		SetError(Error, "matrix is singular");
		return false;
	}
	const int N = A.Rows;
	Matrix Out(N, B.Cols);
	for (int c = 0; c < B.Cols; ++c)
	{
		std::vector<Real> Y(static_cast<size_t>(N));
		for (int i = 0; i < N; ++i)
		{
			Real S = B(Perm[static_cast<size_t>(i)], c);
			for (int j = 0; j < i; ++j)
			{
				S -= LU(i, j) * Y[static_cast<size_t>(j)];
			}
			Y[static_cast<size_t>(i)] = S;
		}
		for (int i = N - 1; i >= 0; --i)
		{
			Real S = Y[static_cast<size_t>(i)];
			for (int j = i + 1; j < N; ++j)
			{
				S -= LU(i, j) * Out(j, c);
			}
			Out(i, c) = S / LU(i, i);
		}
	}
	*X = Out;
	return true;
}

bool Inverse(const Matrix& A, Matrix* Out, std::string* Error)
{
	if (!A.IsSquare())
	{
		SetError(Error, "inv needs a square matrix");
		return false;
	}
	return Solve(A, Matrix::Identity(A.Rows), Out, Error);
}

bool Cholesky(const Matrix& A, Matrix* L, std::string* Error)
{
	if (!A.IsSquare())
	{
		SetError(Error, "chol needs a square matrix");
		return false;
	}
	const int N = A.Rows;
	Matrix Out(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			Real S = 0.5 * (A(i, j) + A(j, i));
			for (int k = 0; k < j; ++k)
			{
				S -= Out(i, k) * Out(j, k);
			}
			if (i == j)
			{
				if (S <= 0.0)
				{
					SetError(Error, "matrix is not positive definite");
					return false;
				}
				Out(i, i) = gm::sqrt(S);
			}
			else
			{
				Out(i, j) = S / Out(j, j);
			}
		}
	}
	*L = Out;
	return true;
}

bool SymmetricEigen(const Matrix& InA, std::vector<Real>* W, Matrix* Vecs, std::string* Error)
{
	if (!InA.IsSquare())
	{
		SetError(Error, "eig needs a square matrix");
		return false;
	}
	const int N = InA.Rows;
	Matrix A = Symmetrize(InA);
	Matrix Vm = Matrix::Identity(N);
	const Real Scale = gm::max(FrobeniusNorm(A), Real(1e-300));
	for (int Sweep = 0; Sweep < 100; ++Sweep)
	{
		Real Off = 0.0;
		for (int p = 0; p < N; ++p)
		{
			for (int q = p + 1; q < N; ++q)
			{
				Off += A(p, q) * A(p, q);
			}
		}
		if (gm::sqrt(Off) <= Scale * 1e-15)
		{
			break;
		}
		for (int p = 0; p < N; ++p)
		{
			for (int q = p + 1; q < N; ++q)
			{
				if (gm::fabs(A(p, q)) <= Scale * 1e-300)
				{
					continue;
				}
				const Real Theta = (A(q, q) - A(p, p)) / (2.0 * A(p, q));
				const Real T = (Theta >= 0.0 ? 1.0 : -1.0) / (gm::fabs(Theta) + gm::sqrt(Theta * Theta + 1.0));
				const Real C = 1.0 / gm::sqrt(T * T + 1.0);
				const Real S = T * C;
				for (int k = 0; k < N; ++k)
				{
					const Real Akp = A(k, p);
					const Real Akq = A(k, q);
					A(k, p) = C * Akp - S * Akq;
					A(k, q) = S * Akp + C * Akq;
				}
				for (int k = 0; k < N; ++k)
				{
					const Real Apk = A(p, k);
					const Real Aqk = A(q, k);
					A(p, k) = C * Apk - S * Aqk;
					A(q, k) = S * Apk + C * Aqk;
				}
				for (int k = 0; k < N; ++k)
				{
					const Real Vkp = Vm(k, p);
					const Real Vkq = Vm(k, q);
					Vm(k, p) = C * Vkp - S * Vkq;
					Vm(k, q) = S * Vkp + C * Vkq;
				}
			}
		}
	}
	std::vector<int> Order(static_cast<size_t>(N));
	for (int i = 0; i < N; ++i)
	{
		Order[static_cast<size_t>(i)] = i;
	}
	std::sort(Order.begin(), Order.end(), [&A](int L, int R) { return A(L, L) < A(R, R); });
	W->assign(static_cast<size_t>(N), 0.0);
	Matrix Sorted(N, N);
	for (int c = 0; c < N; ++c)
	{
		const int Src = Order[static_cast<size_t>(c)];
		(*W)[static_cast<size_t>(c)] = A(Src, Src);
		for (int r = 0; r < N; ++r)
		{
			Sorted(r, c) = Vm(r, Src);
		}
	}
	if (Vecs)
	{
		*Vecs = Sorted;
	}
	return true;
}

bool Eigenvalues(const Matrix& InA, std::vector<Real>* Re, std::vector<Real>* Im, std::string* Error)
{
	if (!InA.IsSquare())
	{
		SetError(Error, "eig needs a square matrix");
		return false;
	}
	const int n = InA.Rows;
	// 1-indexed working copy, following the EISPACK/Numerical Recipes elmhes + hqr pair.
	std::vector<std::vector<Real>> a(static_cast<size_t>(n + 1), std::vector<Real>(static_cast<size_t>(n + 1), 0.0));
	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			a[i][j] = InA(i - 1, j - 1);
		}
	}
	for (int m = 2; m < n; ++m)
	{
		Real x = 0.0;
		int i = m;
		for (int j = m; j <= n; ++j)
		{
			if (gm::fabs(a[j][m - 1]) > gm::fabs(x))
			{
				x = a[j][m - 1];
				i = j;
			}
		}
		if (i != m)
		{
			for (int j = m - 1; j <= n; ++j)
			{
				std::swap(a[i][j], a[m][j]);
			}
			for (int j = 1; j <= n; ++j)
			{
				std::swap(a[j][i], a[j][m]);
			}
		}
		if (x != 0.0)
		{
			for (i = m + 1; i <= n; ++i)
			{
				Real y = a[i][m - 1];
				if (y != 0.0)
				{
					y /= x;
					a[i][m - 1] = y;
					for (int j = m; j <= n; ++j)
					{
						a[i][j] -= y * a[m][j];
					}
					for (int j = 1; j <= n; ++j)
					{
						a[j][m] += y * a[j][i];
					}
				}
			}
		}
	}
	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j < i - 1; ++j)
		{
			a[i][j] = 0.0;
		}
	}

	std::vector<Real> wr(static_cast<size_t>(n + 1), 0.0), wi(static_cast<size_t>(n + 1), 0.0);
	Real anorm = 0.0;
	for (int i = 1; i <= n; ++i)
	{
		for (int j = std::max(i - 1, 1); j <= n; ++j)
		{
			anorm += gm::fabs(a[i][j]);
		}
	}
	int nn = n;
	Real t = 0.0;
	Real p = 0.0, q = 0.0, r = 0.0, s = 0.0, w = 0.0, x = 0.0, y = 0.0, z = 0.0;
	while (nn >= 1)
	{
		int its = 0;
		int l = 1;
		do
		{
			for (l = nn; l >= 2; --l)
			{
				s = gm::fabs(a[l - 1][l - 1]) + gm::fabs(a[l][l]);
				if (s == 0.0)
				{
					s = anorm;
				}
				if (gm::fabs(a[l][l - 1]) + s == s)
				{
					a[l][l - 1] = 0.0;
					break;
				}
			}
			x = a[nn][nn];
			if (l == nn)
			{
				wr[nn] = x + t;
				wi[nn] = 0.0;
				--nn;
			}
			else
			{
				y = a[nn - 1][nn - 1];
				w = a[nn][nn - 1] * a[nn - 1][nn];
				if (l == nn - 1)
				{
					p = 0.5 * (y - x);
					q = p * p + w;
					z = gm::sqrt(gm::fabs(q));
					x += t;
					if (q >= 0.0)
					{
						z = p + (p >= 0.0 ? gm::fabs(z) : -gm::fabs(z));
						wr[nn - 1] = wr[nn] = x + z;
						if (z != 0.0)
						{
							wr[nn] = x - w / z;
						}
						wi[nn - 1] = wi[nn] = 0.0;
					}
					else
					{
						wr[nn - 1] = wr[nn] = x + p;
						wi[nn - 1] = z;
						wi[nn] = -z;
					}
					nn -= 2;
				}
				else
				{
					if (its == 60)
					{
						SetError(Error, "eigenvalue iteration did not converge");
						return false;
					}
					if (its == 10 || its == 20 || its == 40)
					{
						t += x;
						for (int i = 1; i <= nn; ++i)
						{
							a[i][i] -= x;
						}
						s = gm::fabs(a[nn][nn - 1]) + gm::fabs(a[nn - 1][nn - 2]);
						y = x = 0.75 * s;
						w = -0.4375 * s * s;
					}
					++its;
					int m = nn - 2;
					for (; m >= l; --m)
					{
						z = a[m][m];
						r = x - z;
						s = y - z;
						p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
						q = a[m + 1][m + 1] - z - r - s;
						r = a[m + 2][m + 1];
						s = gm::fabs(p) + gm::fabs(q) + gm::fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l)
						{
							break;
						}
						const Real u = gm::fabs(a[m][m - 1]) * (gm::fabs(q) + gm::fabs(r));
						const Real v = gm::fabs(p) * (gm::fabs(a[m - 1][m - 1]) + gm::fabs(z) + gm::fabs(a[m + 1][m + 1]));
						if (u + v == v)
						{
							break;
						}
					}
					for (int i = m + 2; i <= nn; ++i)
					{
						a[i][i - 2] = 0.0;
						if (i != m + 2)
						{
							a[i][i - 3] = 0.0;
						}
					}
					for (int k = m; k <= nn - 1; ++k)
					{
						if (k != m)
						{
							p = a[k][k - 1];
							q = a[k + 1][k - 1];
							r = 0.0;
							if (k != nn - 1)
							{
								r = a[k + 2][k - 1];
							}
							x = gm::fabs(p) + gm::fabs(q) + gm::fabs(r);
							if (x != 0.0)
							{
								p /= x;
								q /= x;
								r /= x;
							}
						}
						const Real Root = gm::sqrt(p * p + q * q + r * r);
						s = p >= 0.0 ? Root : -Root;
						if (s != 0.0)
						{
							if (k == m)
							{
								if (l != m)
								{
									a[k][k - 1] = -a[k][k - 1];
								}
							}
							else
							{
								a[k][k - 1] = -s * x;
							}
							p += s;
							x = p / s;
							y = q / s;
							z = r / s;
							q /= p;
							r /= p;
							for (int j = k; j <= nn; ++j)
							{
								p = a[k][j] + q * a[k + 1][j];
								if (k != nn - 1)
								{
									p += r * a[k + 2][j];
									a[k + 2][j] -= p * z;
								}
								a[k + 1][j] -= p * y;
								a[k][j] -= p * x;
							}
							const int mmin = nn < k + 3 ? nn : k + 3;
							for (int i = l; i <= mmin; ++i)
							{
								p = x * a[i][k] + y * a[i][k + 1];
								if (k != nn - 1)
								{
									p += z * a[i][k + 2];
									a[i][k + 2] -= p * r;
								}
								a[i][k + 1] -= p * q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l < nn - 1);
	}
	std::vector<int> Order(static_cast<size_t>(n));
	for (int i = 0; i < n; ++i)
	{
		Order[static_cast<size_t>(i)] = i + 1;
	}
	std::sort(Order.begin(), Order.end(), [&wr, &wi](int L, int R) {
		if (wr[L] != wr[R])
		{
			return wr[L] > wr[R];
		}
		return wi[L] > wi[R];
	});
	Re->clear();
	Im->clear();
	for (int Idx : Order)
	{
		Re->push_back(wr[Idx]);
		Im->push_back(wi[Idx]);
	}
	return true;
}

bool Svd(const Matrix& InA, Matrix* U, std::vector<Real>* S, Matrix* Vt, std::string* Error)
{
	if (InA.IsEmpty())
	{
		SetError(Error, "svd of an empty matrix");
		return false;
	}
	if (InA.Rows < InA.Cols)
	{
		// A^T = U' S V'^T  =>  A = V' S U'^T.
		Matrix U2, Vt2;
		if (!Svd(Transpose(InA), &U2, S, &Vt2, Error))
		{
			return false;
		}
		if (U)
		{
			*U = Transpose(Vt2);
		}
		if (Vt)
		{
			*Vt = Transpose(U2);
		}
		return true;
	}
	const int M = InA.Rows;
	const int N = InA.Cols;
	Matrix A = InA;
	Matrix Vm = Matrix::Identity(N);
	for (int Sweep = 0; Sweep < 80; ++Sweep)
	{
		bool Rotated = false;
		for (int p = 0; p < N; ++p)
		{
			for (int q = p + 1; q < N; ++q)
			{
				Real Alpha = 0.0, Beta = 0.0, Gamma = 0.0;
				for (int i = 0; i < M; ++i)
				{
					Alpha += A(i, p) * A(i, p);
					Beta += A(i, q) * A(i, q);
					Gamma += A(i, p) * A(i, q);
				}
				if (gm::fabs(Gamma) <= 1e-15 * gm::sqrt(Alpha * Beta) || Gamma == 0.0)
				{
					continue;
				}
				Rotated = true;
				const Real Zeta = (Beta - Alpha) / (2.0 * Gamma);
				const Real T = (Zeta >= 0.0 ? 1.0 : -1.0) / (gm::fabs(Zeta) + gm::sqrt(1.0 + Zeta * Zeta));
				const Real C = 1.0 / gm::sqrt(1.0 + T * T);
				const Real Sn = C * T;
				for (int i = 0; i < M; ++i)
				{
					const Real Ap = A(i, p);
					const Real Aq = A(i, q);
					A(i, p) = C * Ap - Sn * Aq;
					A(i, q) = Sn * Ap + C * Aq;
				}
				for (int i = 0; i < N; ++i)
				{
					const Real Vp = Vm(i, p);
					const Real Vq = Vm(i, q);
					Vm(i, p) = C * Vp - Sn * Vq;
					Vm(i, q) = Sn * Vp + C * Vq;
				}
			}
		}
		if (!Rotated)
		{
			break;
		}
	}
	std::vector<Real> Sig(static_cast<size_t>(N));
	for (int j = 0; j < N; ++j)
	{
		Real Nrm = 0.0;
		for (int i = 0; i < M; ++i)
		{
			Nrm += A(i, j) * A(i, j);
		}
		Sig[static_cast<size_t>(j)] = gm::sqrt(Nrm);
	}
	std::vector<int> Order(static_cast<size_t>(N));
	for (int i = 0; i < N; ++i)
	{
		Order[static_cast<size_t>(i)] = i;
	}
	std::sort(Order.begin(), Order.end(), [&Sig](int L, int R) { return Sig[static_cast<size_t>(L)] > Sig[static_cast<size_t>(R)]; });
	Matrix Uo(M, N), Vo(N, N);
	S->assign(static_cast<size_t>(N), 0.0);
	for (int c = 0; c < N; ++c)
	{
		const int Src = Order[static_cast<size_t>(c)];
		const Real Sv = Sig[static_cast<size_t>(Src)];
		(*S)[static_cast<size_t>(c)] = Sv;
		for (int i = 0; i < M; ++i)
		{
			Uo(i, c) = Sv > 0.0 ? A(i, Src) / Sv : 0.0;
		}
		for (int i = 0; i < N; ++i)
		{
			Vo(i, c) = Vm(i, Src);
		}
	}
	if (U)
	{
		*U = Uo;
	}
	if (Vt)
	{
		*Vt = Transpose(Vo);
	}
	return true;
}

bool LeastSquares(const Matrix& A, const Matrix& B, Matrix* X, std::string* Error)
{
	if (B.Rows != A.Rows)
	{
		SetError(Error, "lsq needs B with as many rows as A");
		return false;
	}
	Matrix U, Vt;
	std::vector<Real> S;
	if (!Svd(A, &U, &S, &Vt, Error))
	{
		return false;
	}
	const Real Tol = static_cast<Real>(std::max(A.Rows, A.Cols)) * Eps() * (S.empty() ? 0.0 : S[0]);
	// X = V S^+ U^T B
	Matrix UtB = Mul(Transpose(U), B);
	for (int i = 0; i < UtB.Rows; ++i)
	{
		const Real Sv = S[static_cast<size_t>(i)];
		for (int j = 0; j < UtB.Cols; ++j)
		{
			UtB(i, j) = Sv > Tol ? UtB(i, j) / Sv : 0.0;
		}
	}
	*X = Mul(Transpose(Vt), UtB);
	return true;
}

bool Expm(const Matrix& InA, Matrix* Out, std::string* Error)
{
	if (!InA.IsSquare())
	{
		SetError(Error, "expm needs a square matrix");
		return false;
	}
	static const Real B[14] = {64764752532480000.0, 32382376266240000.0, 7771770303897600.0, 1187353796428800.0, 129060195264000.0,
		10559470521600.0, 670442572800.0, 33522128640.0, 1323241920.0, 40840800.0, 960960.0, 16380.0, 182.0, 1.0};
	const Real Theta13 = 5.371920351148152;
	const int N = InA.Rows;
	const Real Nrm = Norm1(InA);
	if (!gm::isfinite(Nrm))
	{
		SetError(Error, "expm of a non-finite matrix");
		return false;
	}
	int Squarings = 0;
	if (Nrm > Theta13)
	{
		Squarings = gm::ceil(gm::log2(Nrm / Theta13)).ToInt();
	}
	const Matrix A = Scale(InA, std::ldexp(1.0, -Squarings));
	const Matrix I = Matrix::Identity(N);
	const Matrix A2 = Mul(A, A), A4 = Mul(A2, A2), A6 = Mul(A4, A2);
	Matrix Inner = Add(Add(Scale(A6, B[13]), Scale(A4, B[11])), Scale(A2, B[9]));
	Matrix Uu = Add(Add(Add(Add(Mul(A6, Inner), Scale(A6, B[7])), Scale(A4, B[5])), Scale(A2, B[3])), Scale(I, B[1]));
	Uu = Mul(A, Uu);
	Matrix Inner2 = Add(Add(Scale(A6, B[12]), Scale(A4, B[10])), Scale(A2, B[8]));
	Matrix Vv = Add(Add(Add(Add(Mul(A6, Inner2), Scale(A6, B[6])), Scale(A4, B[4])), Scale(A2, B[2])), Scale(I, B[0]));
	Matrix R;
	if (!Solve(Sub(Vv, Uu), Add(Vv, Uu), &R, Error))
	{
		return false;
	}
	for (int k = 0; k < Squarings; ++k)
	{
		R = Mul(R, R);
	}
	*Out = R;
	return true;
}

Real SpectralRadius(const Matrix& A)
{
	std::vector<Real> Re, Im;
	if (!Eigenvalues(A, &Re, &Im, nullptr))
	{
		return Real::NaN();
	}
	Real Rho = 0.0;
	for (size_t i = 0; i < Re.size(); ++i)
	{
		Rho = std::max(Rho, gm::hypot(Re[i], Im[i]));
	}
	return Rho;
}

int Rank(const Matrix& A, Real Tol)
{
	if (A.IsEmpty())
	{
		return 0;
	}
	std::vector<Real> S;
	if (!Svd(A, nullptr, &S, nullptr, nullptr))
	{
		return 0;
	}
	if (Tol < 0.0)
	{
		Tol = static_cast<Real>(std::max(A.Rows, A.Cols)) * Eps() * (S.empty() ? 0.0 : S[0]);
	}
	int R = 0;
	for (Real Sv : S)
	{
		if (Sv > Tol)
		{
			++R;
		}
	}
	return R;
}

Matrix Kron(const Matrix& A, const Matrix& B)
{
	Matrix K(A.Rows * B.Rows, A.Cols * B.Cols);
	for (int i = 0; i < A.Rows; ++i)
	{
		for (int j = 0; j < A.Cols; ++j)
		{
			for (int k = 0; k < B.Rows; ++k)
			{
				for (int l = 0; l < B.Cols; ++l)
				{
					K(i * B.Rows + k, j * B.Cols + l) = A(i, j) * B(k, l);
				}
			}
		}
	}
	return K;
}

bool Dare(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, Matrix* P, Matrix* K, std::string* Error)
{
	const int N = A.Rows;
	if (!A.IsSquare() || B.Rows != N || !Q.IsSquare() || Q.Rows != N || !R.IsSquare() || R.Rows != B.Cols)
	{
		SetError(Error, "dare(A,B,Q,R): A n x n, B n x m, Q n x n, R m x m");
		return false;
	}
	Matrix Rinv;
	if (!Inverse(R, &Rinv, Error))
	{
		SetError(Error, "dare: R must be invertible");
		return false;
	}
	// Structure-preserving doubling (Chu, Fan, Lin, Wang 2004).
	Matrix Ak = A;
	Matrix Gk = Mul(Mul(B, Rinv), Transpose(B));
	Matrix Hk = Symmetrize(Q);
	const Matrix I = Matrix::Identity(N);
	bool Converged = false;
	for (int It = 0; It < 200; ++It)
	{
		Matrix W = Add(I, Mul(Gk, Hk));
		Matrix WiA, WiG;
		if (!Solve(W, Ak, &WiA, Error) || !Solve(W, Gk, &WiG, Error))
		{
			SetError(Error, "dare: doubling step is singular (is (A,B) stabilizable?)");
			return false;
		}
		Matrix Hn = Symmetrize(Add(Hk, Mul(Mul(Transpose(Ak), Hk), WiA)));
		Matrix Gn = Symmetrize(Add(Gk, Mul(Mul(Ak, WiG), Transpose(Ak))));
		Matrix An = Mul(Ak, WiA);
		const Real Change = MaxAbs(Sub(Hn, Hk));
		const Real Size = gm::max(Real(1), MaxAbs(Hn));
		Ak = An;
		Gk = Gn;
		Hk = Hn;
		if (!gm::isfinite(Change))
		{
			break;
		}
		if (Change <= 1e-14 * Size)
		{
			Converged = true;
			break;
		}
	}
	if (!Converged)
	{
		SetError(Error, "dare did not converge (check stabilizability and detectability)");
		return false;
	}
	Matrix BtP = Mul(Transpose(B), Hk);
	Matrix Gain;
	if (!Solve(Add(R, Mul(BtP, B)), Mul(BtP, A), &Gain, Error))
	{
		return false;
	}
	*P = Hk;
	if (K)
	{
		*K = Gain;
	}
	return true;
}

bool Care(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R, Matrix* P, Matrix* K, std::string* Error)
{
	const int N = A.Rows;
	if (!A.IsSquare() || B.Rows != N || !Q.IsSquare() || Q.Rows != N || !R.IsSquare() || R.Rows != B.Cols)
	{
		SetError(Error, "care(A,B,Q,R): A n x n, B n x m, Q n x n, R m x m");
		return false;
	}
	Matrix Rinv;
	if (!Inverse(R, &Rinv, Error))
	{
		SetError(Error, "care: R must be invertible");
		return false;
	}
	const Matrix G = Mul(Mul(B, Rinv), Transpose(B));
	Matrix Z(2 * N, 2 * N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Z(i, j) = A(i, j);
			Z(i, j + N) = -G(i, j);
			Z(i + N, j) = -Q(i, j);
			Z(i + N, j + N) = -A(j, i);
		}
	}
	bool Converged = false;
	for (int It = 0; It < 100; ++It)
	{
		Matrix Zi;
		if (!Inverse(Z, &Zi, Error))
		{
			SetError(Error, "care: Hamiltonian has eigenvalues on the imaginary axis");
			return false;
		}
		Real D = 0.0;
		Determinant(Z, &D, nullptr);
		Real C = gm::pow(gm::fabs(D), -1.0 / (2.0 * N));
		if (!gm::isfinite(C) || C <= 0.0)
		{
			C = 1.0;
		}
		Matrix Zn = Scale(Add(Scale(Z, C), Scale(Zi, 1.0 / C)), 0.5);
		const Real Change = MaxAbs(Sub(Zn, Z));
		Z = Zn;
		if (Change <= 1e-13 * gm::max(Real(1), MaxAbs(Z)))
		{
			Converged = true;
			break;
		}
	}
	if (!Converged)
	{
		SetError(Error, "care: sign iteration did not converge");
		return false;
	}
	// [W12; W22 + I] P = -[W11 + I; W21]
	Matrix Lhs(2 * N, N), Rhs(2 * N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			Lhs(i, j) = Z(i, j + N);
			Lhs(i + N, j) = Z(i + N, j + N) + (i == j ? 1.0 : 0.0);
			Rhs(i, j) = -(Z(i, j) + (i == j ? 1.0 : 0.0));
			Rhs(i + N, j) = -Z(i + N, j);
		}
	}
	Matrix Sol;
	if (!LeastSquares(Lhs, Rhs, &Sol, Error))
	{
		return false;
	}
	*P = Symmetrize(Sol);
	if (K)
	{
		*K = Mul(Mul(Rinv, Transpose(B)), *P);
	}
	return true;
}

bool Dlyap(const Matrix& A, const Matrix& Q, Matrix* X, std::string* Error)
{
	const int N = A.Rows;
	if (!A.IsSquare() || Q.Rows != N || Q.Cols != N)
	{
		SetError(Error, "dlyap(A,Q): A and Q must be the same square size");
		return false;
	}
	if (N > 12)
	{
		SetError(Error, "dlyap: n <= 12");
		return false;
	}
	Matrix M = Sub(Matrix::Identity(N * N), Kron(A, A));
	Matrix Rhs(N * N, 1);
	Rhs.V = Q.V;
	Matrix Sol;
	if (!Solve(M, Rhs, &Sol, Error))
	{
		SetError(Error, "dlyap: A has reciprocal eigenvalue pairs (no unique solution)");
		return false;
	}
	Matrix Out(N, N);
	Out.V = Sol.V;
	*X = Out;
	return true;
}

bool Lyap(const Matrix& A, const Matrix& Q, Matrix* X, std::string* Error)
{
	const int N = A.Rows;
	if (!A.IsSquare() || Q.Rows != N || Q.Cols != N)
	{
		SetError(Error, "lyap(A,Q): A and Q must be the same square size");
		return false;
	}
	if (N > 12)
	{
		SetError(Error, "lyap: n <= 12");
		return false;
	}
	const Matrix I = Matrix::Identity(N);
	Matrix M = Add(Kron(A, I), Kron(I, A));
	Matrix Rhs(N * N, 1);
	for (size_t i = 0; i < Q.V.size(); ++i)
	{
		Rhs.V[i] = -Q.V[i];
	}
	Matrix Sol;
	if (!Solve(M, Rhs, &Sol, Error))
	{
		SetError(Error, "lyap: A and -A share an eigenvalue (no unique solution)");
		return false;
	}
	Matrix Out(N, N);
	Out.V = Sol.V;
	*X = Out;
	return true;
}
} // namespace la
} // namespace mcalc
