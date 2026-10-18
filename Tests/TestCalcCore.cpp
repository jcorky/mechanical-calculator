// The expression language and the linear algebra underneath every library function.
#include "CalcSession.h"
#include "TestFramework.h"

using mcalc::Matrix;
using mcalc::Session;
namespace la = mcalc::la;

namespace
{
double Eval(Session& S, const char* Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue)
	{
		std::printf("    '%s' -> %s\n", Line, O.Text.c_str());
		return std::nan("");
	}
	return O.Val.Scalar();
}

bool Rejects(Session& S, const char* Line)
{
	return !S.Execute(Line).Ok;
}

double Residual(const Matrix& A, const Matrix& B)
{
	return la::MaxAbs(la::Sub(A, B));
}
} // namespace

TEST(parser_precedence_and_associativity)
{
	Session S;
	CHECK_NEAR(Eval(S, "2 + 3*4"), 14, 0);
	CHECK_NEAR(Eval(S, "(2 + 3)*4"), 20, 0);
	CHECK_NEAR(Eval(S, "2^3^2"), 512, 0);
	CHECK_NEAR(Eval(S, "-2^2"), -4, 0);
	CHECK_NEAR(Eval(S, "2^-1"), 0.5, 0);
	CHECK_NEAR(Eval(S, "10/4/5"), 0.5, 1e-15);
	CHECK_NEAR(Eval(S, "1.5e3 + .5"), 1500.5, 0);
	CHECK_NEAR(Eval(S, "5!"), 120, 1e-9);
	CHECK_NEAR(Eval(S, "--3"), 3, 0);
	CHECK_NEAR(Eval(S, "6 \\ 3"), 0.5, 0);
	CHECK(Rejects(S, "2 +"));
	CHECK(Rejects(S, "(1"));
	CHECK(Rejects(S, "1 2"));
	CHECK(Rejects(S, "foo(1)"));
	CHECK(Rejects(S, "2 $ 3"));
}

TEST(variables_ans_and_statements)
{
	Session S;
	CHECK_NEAR(Eval(S, "x = 3; y = x^2; y + 1"), 10, 0);
	CHECK_NEAR(Eval(S, "ans * 2"), 20, 0);
	CHECK_NEAR(Eval(S, "x"), 3, 0);
	CHECK(Rejects(S, "sin = 3"));
	CHECK(Rejects(S, "pi = 3"));
	CHECK(S.Execute("clear").Ok);
	CHECK(Rejects(S, "x"));
}

TEST(constants_and_angle_mode)
{
	Session S;
	CHECK_NEAR(Eval(S, "sin(pi/6)"), 0.5, 1e-15);
	CHECK(S.Execute("deg").Ok);
	CHECK_NEAR(Eval(S, "sin(30)"), 0.5, 1e-15);
	CHECK_NEAR(Eval(S, "atan2(1, 1)"), 45, 1e-13);
	CHECK_NEAR(Eval(S, "deg(pi)"), 180, 1e-12);   // conversions ignore the mode
	CHECK(S.Execute("rad").Ok);
	CHECK_NEAR(Eval(S, "acos(-1)"), 3.141592653589793, 1e-15);
	CHECK_NEAR(Eval(S, "3*nmi/hr"), 3 * 1852.0 / 3600.0, 1e-12);
}

TEST(matrices_and_operators)
{
	Session S;
	CHECK_NEAR(Eval(S, "A = [1, 2; 3, 4]; det(A)"), -2, 1e-12);
	CHECK_NEAR(Eval(S, "B = A * A; B(2, 1)"), 15, 0);
	CHECK_NEAR(Eval(S, "C = A'; C(1, 2)"), 3, 0);
	CHECK_NEAR(Eval(S, "x = A \\ [5; 6]; x(2)"), 4.5, 1e-12);
	CHECK_NEAR(Eval(S, "y = [5, 6] / A; y(1)"), -1, 1e-12);
	CHECK_NEAR(Eval(S, "D = A .* A; D(2, 2)"), 16, 0);
	CHECK_NEAR(Eval(S, "E = A ^ -1; E(1, 1)"), -2, 1e-12);
	CHECK_NEAR(Eval(S, "F = A .^ 2; F(1, 2)"), 4, 0);
	CHECK_NEAR(Eval(S, "G = [A, [9; 9]; 7, 7, 7]; G(3, 3)"), 7, 0);
	CHECK_NEAR(Eval(S, "sum([1, 2, 3] + 1)"), 9, 0);
	CHECK_NEAR(Eval(S, "max(3, 9, -1)"), 9, 0);
	CHECK_NEAR(Eval(S, "trace(eye(4))"), 4, 0);
	CHECK(Rejects(S, "[1, 2] + [1, 2, 3]"));
	CHECK(Rejects(S, "[1, 2] * [1, 2]"));
	CHECK(Rejects(S, "[1, 2; 3]"));
	CHECK(Rejects(S, "A(3, 1)"));
}

TEST(records_print_and_index)
{
	Session S;
	const mcalc::Outcome O = S.Execute("s = size(zeros(2, 5))");
	CHECK(O.Ok && O.Val.IsRecord());
	CHECK(O.Lines.size() == 2);
	CHECK(O.Lines[0] == "rows = 2");
	CHECK_NEAR(Eval(S, "s.cols"), 5, 0);
	CHECK_NEAR(Eval(S, "s(1)"), 2, 0);
	CHECK_NEAR(Eval(S, "size([1, 2, 3]).cols"), 3, 0);
	CHECK(Rejects(S, "s.nope"));
}

TEST(number_formatting)
{
	Session S;
	CHECK(S.FormatNumber(1500000) == "1500000");
	CHECK(S.FormatNumber(1.5e20) == "1.5e20");
	CHECK(S.FormatNumber(2e-7) == "2e-7");
	CHECK(S.FormatNumber(-0.0) == "0");
	S.SetDigits(4);
	CHECK(S.FormatNumber(3.14159265) == "3.142");
	CHECK(S.Format(mcalc::Value(Matrix::FromRows({{1, 22}, {333, 4}})))[1] == "[ 333  4 ]");
}

TEST(commands_help_list)
{
	Session S;
	const mcalc::Outcome H = S.Execute("help sqrt");
	CHECK(H.Ok && H.Text.find("sqrt(x)") == 0);
	CHECK(S.Execute("list").Ok);
	CHECK(S.Execute("list matrix").Ok);
	CHECK(!S.Execute("list nonsense").Ok);
	CHECK(!S.Execute("help nonsense").Ok);
	const std::vector<std::string> C = S.Complete("sq");
	CHECK(!C.empty() && C[0] == "sqrt");
	CHECK(mcalc::WrapText("aaa bbb ccc", 7).size() == 2);
}

TEST(linear_algebra_solve_inverse_det)
{
	const Matrix A = Matrix::FromRows({{4, -2, 1}, {-2, 4, -2}, {1, -2, 4}});
	Matrix Ai, X;
	CHECK(la::Inverse(A, &Ai, nullptr));
	CHECK(Residual(la::Mul(A, Ai), Matrix::Identity(3)) < 1e-14);
	double D = 0;
	CHECK(la::Determinant(A, &D, nullptr));
	CHECK_NEAR(D, 36, 1e-12);
	std::string E;
	CHECK(!la::Solve(Matrix::FromRows({{1, 2}, {2, 4}}), Matrix::Column({1, 1}), &X, &E));
	Matrix L;
	CHECK(la::Cholesky(A, &L, nullptr));
	CHECK(Residual(la::Mul(L, la::Transpose(L)), A) < 1e-14);
	// Overdetermined least squares: the normal equations hold at the solution.
	const Matrix M = Matrix::FromRows({{1, 1}, {1, 2}, {1, 3}, {1, 4}});
	const Matrix b = Matrix::Column({6, 5, 7, 10});
	CHECK(la::LeastSquares(M, b, &X, nullptr));
	CHECK_NEAR(X.V[0], 3.5, 1e-12);
	CHECK_NEAR(X.V[1], 1.4, 1e-12);
}

TEST(linear_algebra_eigen_svd)
{
	std::vector<double> W, Re, Im, S;
	Matrix V, U, Vt;
	const Matrix A = Matrix::FromRows({{2, 1, 0}, {1, 3, 1}, {0, 1, 4}});
	CHECK(la::SymmetricEigen(A, &W, &V, nullptr));
	// A V = V diag(W)
	Matrix Dg(3, 3);
	for (int i = 0; i < 3; ++i)
	{
		Dg(i, i) = W[static_cast<size_t>(i)];
	}
	CHECK(Residual(la::Mul(A, V), la::Mul(V, Dg)) < 1e-13);
	CHECK(la::Eigenvalues(A, &Re, &Im, nullptr));
	for (int i = 0; i < 3; ++i)
	{
		CHECK_NEAR(Re[static_cast<size_t>(i)], W[static_cast<size_t>(2 - i)], 1e-12);
	}
	// Rotation by 0.3 rad: eigenvalues cos 0.3 +- i sin 0.3.
	const Matrix Rot = Matrix::FromRows({{std::cos(0.3), -std::sin(0.3)}, {std::sin(0.3), std::cos(0.3)}});
	CHECK(la::Eigenvalues(Rot, &Re, &Im, nullptr));
	CHECK_NEAR(Re[0], std::cos(0.3), 1e-14);
	CHECK_NEAR(Im[0], std::sin(0.3), 1e-14);
	CHECK_NEAR(Im[1], -std::sin(0.3), 1e-14);
	// Companion matrix of (x-1)(x-2)(x-3)(x+4)(x^2+1).
	Session Sn;
	const mcalc::Outcome O = Sn.Execute("r = roots(conv(conv([1, -1], [1, -2]), conv(conv([1, -3], [1, 4]), [1, 0, 1])))");
	CHECK(O.Ok && O.Val.M.Rows == 6 && O.Val.M.Cols == 2);
	CHECK_NEAR(O.Val.M(0, 0), 3, 1e-10);
	CHECK_NEAR(O.Val.M(5, 0), -4, 1e-10);
	// Non-normal 5 x 5 with known spectrum: triangular T similar-transformed.
	Matrix T = Matrix::FromRows({{1, 5, -3, 2, 7}, {0, -2, 4, 1, 1}, {0, 0, 3.5, -6, 2}, {0, 0, 0, 0.25, 9}, {0, 0, 0, 0, -7}});
	Matrix P = Matrix::FromRows({{1, 2, 0, 1, 0}, {0, 1, 3, 0, 1}, {1, 0, 1, 2, 0}, {0, 1, 0, 1, 3}, {2, 0, 1, 0, 1}});
	Matrix Pi;
	CHECK(la::Inverse(P, &Pi, nullptr));
	CHECK(la::Eigenvalues(la::Mul(la::Mul(P, T), Pi), &Re, &Im, nullptr));
	const double Expect[] = {3.5, 1, 0.25, -2, -7};
	for (int i = 0; i < 5; ++i)
	{
		CHECK_NEAR(Re[static_cast<size_t>(i)], Expect[i], 1e-9);
	}
	const Matrix B = Matrix::FromRows({{3, 2, 2}, {2, 3, -2}});
	CHECK(la::Svd(B, &U, &S, &Vt, nullptr));
	CHECK_NEAR(S[0], 5, 1e-13);
	CHECK_NEAR(S[1], 3, 1e-13);
	Matrix Sd(2, 2);
	Sd(0, 0) = S[0];
	Sd(1, 1) = S[1];
	CHECK(Residual(la::Mul(la::Mul(U, Sd), Vt), B) < 1e-13);
	CHECK(la::Rank(Matrix::FromRows({{1, 2, 3}, {2, 4, 6}, {1, 0, 1}})) == 2);
}

TEST(linear_algebra_expm)
{
	// exp of a rotation generator and of a nilpotent matrix, and of a large-norm diagonal.
	Matrix E;
	CHECK(la::Expm(Matrix::FromRows({{0, -2}, {2, 0}}), &E, nullptr));
	CHECK_NEAR(E(0, 0), std::cos(2.0), 1e-14);
	CHECK_NEAR(E(1, 0), std::sin(2.0), 1e-14);
	CHECK(la::Expm(Matrix::FromRows({{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}), &E, nullptr));
	CHECK_NEAR(E(0, 2), 0.5, 1e-15);
	CHECK(la::Expm(Matrix::FromRows({{-30, 0}, {0, 12}}), &E, nullptr));
	CHECK_NEAR(E(0, 0) / std::exp(-30.0), 1, 1e-12);
	CHECK_NEAR(E(1, 1) / std::exp(12.0), 1, 1e-12);
}

TEST(riccati_and_lyapunov)
{
	// Double integrator, dt = 0.1: SDA against 5000 plain Riccati iterations.
	const Matrix A = Matrix::FromRows({{1, 0.1}, {0, 1}});
	const Matrix B = Matrix::Column({0.005, 0.1});
	const Matrix Q = Matrix::Identity(2);
	const Matrix R = Matrix::Scalar(0.5);
	Matrix P, K;
	CHECK(la::Dare(A, B, Q, R, &P, &K, nullptr));
	Matrix It = Q;
	for (int i = 0; i < 5000; ++i)
	{
		const Matrix BtP = la::Mul(la::Transpose(B), It);
		Matrix G;
		la::Solve(la::Add(R, la::Mul(BtP, B)), la::Mul(BtP, A), &G, nullptr);
		It = la::Add(Q, la::Sub(la::Mul(la::Mul(la::Transpose(A), It), A), la::Mul(la::Mul(la::Transpose(A), la::Transpose(BtP)), G)));
	}
	CHECK(Residual(P, It) < 1e-9 * la::MaxAbs(P));
	CHECK(la::SpectralRadius(la::Sub(A, la::Mul(B, K))) < 1);
	// CARE on the same plant in continuous time: residual and stability.
	const Matrix Ac = Matrix::FromRows({{0, 1}, {0, 0}});
	const Matrix Bc = Matrix::Column({0, 1});
	CHECK(la::Care(Ac, Bc, Q, Matrix::Scalar(1), &P, &K, nullptr));
	const Matrix Res = la::Add(la::Add(la::Mul(la::Transpose(Ac), P), la::Mul(P, Ac)), la::Sub(Q, la::Mul(la::Mul(P, la::Mul(Bc, la::Transpose(Bc))), P)));
	CHECK(la::MaxAbs(Res) < 1e-10);
	CHECK_NEAR(P(0, 0), std::sqrt(3.0), 1e-10);   // known closed form: P = [sqrt3, 1; 1, sqrt3]
	CHECK_NEAR(P(0, 1), 1, 1e-10);
	// Lyapunov equations.
	const Matrix F = Matrix::FromRows({{0.9, 0.2}, {-0.1, 0.7}});
	Matrix X;
	CHECK(la::Dlyap(F, Q, &X, nullptr));
	CHECK(Residual(X, la::Add(la::Mul(la::Mul(F, X), la::Transpose(F)), Q)) < 1e-13);
	const Matrix Fc = Matrix::FromRows({{-1, 2}, {0, -3}});
	CHECK(la::Lyap(Fc, Q, &X, nullptr));
	CHECK(la::MaxAbs(la::Add(la::Add(la::Mul(Fc, X), la::Mul(X, la::Transpose(Fc))), Q)) < 1e-13);
}

TEST(registry_names_are_unique_and_documented)
{
	const mcalc::Registry& R = mcalc::Registry::Default();
	for (const std::string& D : R.Duplicates())
	{
		std::printf("    duplicate function name: %s\n", D.c_str());
	}
	CHECK(R.Duplicates().empty());
	for (const mcalc::FunctionInfo& F : R.All())
	{
		CHECK(!F.Summary.empty() && !F.Category.empty() && F.Fn != nullptr);
	}
	std::printf("    %zu functions in %zu categories\n", R.All().size(), R.Categories().size());
}
