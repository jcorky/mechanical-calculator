#include "CalcSession.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace mcalc
{
namespace session_detail
{
enum class Tok
{
	End,
	Number,
	Ident,
	Plus,
	Minus,
	Star,
	Slash,
	Backslash,
	Caret,
	DotStar,
	DotSlash,
	DotCaret,
	LParen,
	RParen,
	LBracket,
	RBracket,
	Comma,
	Semicolon,
	Quote,
	Dot,
	Bang,
	Equals,
	Bad,
};

struct Token
{
	Tok Kind = Tok::End;
	double Number = 0.0;
	std::string Text;
	size_t Pos = 0;
};

bool IsIdentStart(char Ch)
{
	return std::isalpha(static_cast<unsigned char>(Ch)) != 0 || Ch == '_';
}

bool IsIdentChar(char Ch)
{
	return std::isalnum(static_cast<unsigned char>(Ch)) != 0 || Ch == '_';
}

std::vector<Token> Lex(const std::string& S)
{
	std::vector<Token> Out;
	size_t i = 0;
	while (i < S.size())
	{
		const char Ch = S[i];
		if (std::isspace(static_cast<unsigned char>(Ch)) != 0)
		{
			++i;
			continue;
		}
		Token T;
		T.Pos = i;
		const bool DotDigit = Ch == '.' && i + 1 < S.size() && std::isdigit(static_cast<unsigned char>(S[i + 1])) != 0;
		if (std::isdigit(static_cast<unsigned char>(Ch)) != 0 || DotDigit)
		{
			const char* Begin = S.c_str() + i;
			char* EndPtr = nullptr;
			T.Number = std::strtod(Begin, &EndPtr);
			size_t Len = static_cast<size_t>(EndPtr - Begin);
			if (Len == 0)
			{
				Len = 1;
			}
			T.Kind = Tok::Number;
			T.Text = S.substr(i, Len);
			i += Len;
			Out.push_back(T);
			continue;
		}
		if (IsIdentStart(Ch))
		{
			size_t j = i + 1;
			while (j < S.size() && IsIdentChar(S[j]))
			{
				++j;
			}
			T.Kind = Tok::Ident;
			T.Text = S.substr(i, j - i);
			i = j;
			Out.push_back(T);
			continue;
		}
		const char Next = i + 1 < S.size() ? S[i + 1] : '\0';
		size_t Len = 1;
		switch (Ch)
		{
		case '+': T.Kind = Tok::Plus; break;
		case '-': T.Kind = Tok::Minus; break;
		case '*': T.Kind = Tok::Star; break;
		case '/': T.Kind = Tok::Slash; break;
		case '\\': T.Kind = Tok::Backslash; break;
		case '^': T.Kind = Tok::Caret; break;
		case '(': T.Kind = Tok::LParen; break;
		case ')': T.Kind = Tok::RParen; break;
		case '[': T.Kind = Tok::LBracket; break;
		case ']': T.Kind = Tok::RBracket; break;
		case ',': T.Kind = Tok::Comma; break;
		case ';': T.Kind = Tok::Semicolon; break;
		case '\'': T.Kind = Tok::Quote; break;
		case '!': T.Kind = Tok::Bang; break;
		case '=': T.Kind = Tok::Equals; break;
		case '.':
			if (Next == '*') { T.Kind = Tok::DotStar; Len = 2; }
			else if (Next == '/') { T.Kind = Tok::DotSlash; Len = 2; }
			else if (Next == '^') { T.Kind = Tok::DotCaret; Len = 2; }
			else { T.Kind = Tok::Dot; }
			break;
		default: T.Kind = Tok::Bad; break;
		}
		T.Text = S.substr(i, Len);
		i += Len;
		Out.push_back(T);
	}
	Token EndTok;
	EndTok.Kind = Tok::End;
	EndTok.Pos = S.size();
	Out.push_back(EndTok);
	return Out;
}

bool Constant(const std::string& Name, double* Out)
{
	struct Entry
	{
		const char* Name;
		double X;
	};
	static const Entry Table[] = {
		{"pi", cst::Pi},
		{"tau", cst::Tau},
		{"e", 2.71828182845904523536},
		{"inf", std::numeric_limits<double>::infinity()},
		{"nan", std::numeric_limits<double>::quiet_NaN()},
		{"g0", cst::G0},
		{"c0", cst::C},
		{"kB", cst::Boltzmann},
		{"T0", cst::T0Noise},
		{"mu_earth", cst::MuEarth},
		{"R_earth", cst::REarth},
		{"gamma_air", cst::GammaAir},
		{"cp_air", cst::CpAir},
		{"R_air", cst::RAir},
		{"ft", 0.3048},
		{"nmi", 1852.0},
		{"kt", 1852.0 / 3600.0},
		{"km", 1000.0},
		{"hr", 3600.0},
		{"day", 86400.0},
	};
	for (const Entry& E : Table)
	{
		if (Name == E.Name)
		{
			*Out = E.X;
			return true;
		}
	}
	return false;
}

const char* const kConstantNames[] = {"pi", "tau", "e", "inf", "nan", "g0", "c0", "kB", "T0", "mu_earth", "R_earth", "gamma_air", "cp_air", "R_air", "ft", "nmi", "kt", "km", "hr", "day"};
const char* const kCommandNames[] = {"help", "list", "vars", "clear", "deg", "rad", "digits"};

std::string Shape(const Matrix& M)
{
	return std::to_string(M.Rows) + "x" + std::to_string(M.Cols);
}

// Element-wise binary operation with scalar broadcasting.
bool Elementwise(const Matrix& A, const Matrix& B, double (*Op)(double, double), Matrix* Out, std::string* Error, const char* What)
{
	if (A.IsScalar() && !B.IsScalar())
	{
		Matrix C = B;
		for (double& X : C.V)
		{
			X = Op(A.V[0], X);
		}
		*Out = C;
		return true;
	}
	if (B.IsScalar())
	{
		Matrix C = A;
		for (double& X : C.V)
		{
			X = Op(X, B.V[0]);
		}
		*Out = C;
		return true;
	}
	if (A.Rows != B.Rows || A.Cols != B.Cols)
	{
		*Error = std::string(What) + ": shapes " + Shape(A) + " and " + Shape(B) + " differ";
		return false;
	}
	Matrix C = A;
	for (size_t i = 0; i < C.V.size(); ++i)
	{
		C.V[i] = Op(A.V[i], B.V[i]);
	}
	*Out = C;
	return true;
}

double OpAdd(double A, double B) { return A + B; }
double OpSub(double A, double B) { return A - B; }
double OpMul(double A, double B) { return A * B; }
double OpDiv(double A, double B) { return A / B; }
double OpPow(double A, double B) { return std::pow(A, B); }

class Parser
{
public:
	Parser(Session& InSession, const std::vector<Token>& InTokens) : S(InSession), T(InTokens) {}

	bool Failed() const { return !Err.empty(); }
	const std::string& Error() const { return Err; }
	bool AtEnd() const { return Peek().Kind == Tok::End; }
	const Token& Peek(size_t Ahead = 0) const { return T[std::min(Index + Ahead, T.size() - 1)]; }

	Value Expression()
	{
		Value L = Term();
		while (!Failed())
		{
			const Tok K = Peek().Kind;
			if (K != Tok::Plus && K != Tok::Minus)
			{
				break;
			}
			++Index;
			Value R = Term();
			if (Failed())
			{
				break;
			}
			Matrix M;
			if (!Elementwise(L.M, R.M, K == Tok::Plus ? OpAdd : OpSub, &M, &Err, K == Tok::Plus ? "+" : "-"))
			{
				break;
			}
			L = Value(M);
		}
		return L;
	}

	void Fail(const std::string& Message)
	{
		if (Err.empty())
		{
			Err = Message;
		}
	}

private:
	Session& S;
	const std::vector<Token>& T;
	size_t Index = 0;
	std::string Err;

	bool Accept(Tok K)
	{
		if (Peek().Kind == K)
		{
			++Index;
			return true;
		}
		return false;
	}

	bool Expect(Tok K, const char* What)
	{
		if (Accept(K))
		{
			return true;
		}
		Fail(std::string("expected ") + What + Near());
		return false;
	}

	std::string Near() const
	{
		const Token& P = Peek();
		if (P.Kind == Tok::End)
		{
			return " at end of input";
		}
		return " near '" + P.Text + "'";
	}

	Value Term()
	{
		Value L = Unary();
		while (!Failed())
		{
			const Tok K = Peek().Kind;
			if (K != Tok::Star && K != Tok::Slash && K != Tok::Backslash && K != Tok::DotStar && K != Tok::DotSlash)
			{
				break;
			}
			++Index;
			Value R = Unary();
			if (Failed())
			{
				break;
			}
			Matrix M;
			if (K == Tok::DotStar || (K == Tok::Star && (L.IsScalar() || R.IsScalar())))
			{
				if (!Elementwise(L.M, R.M, OpMul, &M, &Err, "*"))
				{
					break;
				}
			}
			else if (K == Tok::DotSlash || (K == Tok::Slash && R.IsScalar()))
			{
				if (!Elementwise(L.M, R.M, OpDiv, &M, &Err, "/"))
				{
					break;
				}
			}
			else if (K == Tok::Star)
			{
				if (L.M.Cols != R.M.Rows)
				{
					Fail("*: inner dimensions " + Shape(L.M) + " * " + Shape(R.M) + " differ");
					break;
				}
				M = la::Mul(L.M, R.M);
			}
			else if (K == Tok::Slash)
			{
				// A / B = A inv(B): solve B' X' = A'.
				Matrix Xt;
				if (!la::Solve(la::Transpose(R.M), la::Transpose(L.M), &Xt, &Err))
				{
					Err = "/: " + Err;
					break;
				}
				M = la::Transpose(Xt);
			}
			else
			{
				if (L.IsScalar())
				{
					Elementwise(R.M, L.M, OpDiv, &M, &Err, "\\");
				}
				else if (L.M.Rows != R.M.Rows)
				{
					Fail("\\: A and b need the same number of rows");
					break;
				}
				else if (L.M.IsSquare())
				{
					if (!la::Solve(L.M, R.M, &M, &Err))
					{
						Err = "\\: " + Err;
						break;
					}
				}
				else if (!la::LeastSquares(L.M, R.M, &M, &Err))
				{
					break;
				}
			}
			L = Value(M);
		}
		return L;
	}

	Value Unary()
	{
		if (Accept(Tok::Minus))
		{
			Value V = Unary();
			for (double& X : V.M.V)
			{
				X = -X;
			}
			V.Fields.clear();
			return V;
		}
		if (Accept(Tok::Plus))
		{
			return Unary();
		}
		return Power();
	}

	Value Power()
	{
		Value Base = Postfix();
		if (Failed())
		{
			return Base;
		}
		const Tok K = Peek().Kind;
		if (K != Tok::Caret && K != Tok::DotCaret)
		{
			return Base;
		}
		++Index;
		Value Exp = Unary();   // right-associative, and 2^-1 works
		if (Failed())
		{
			return Base;
		}
		Matrix M;
		if (K == Tok::DotCaret || (Base.IsScalar() && Exp.IsScalar()))
		{
			Elementwise(Base.M, Exp.M, OpPow, &M, &Err, "^");
			return Value(M);
		}
		if (!Exp.IsScalar() || !Base.M.IsSquare())
		{
			Fail("^: use a square matrix and a scalar power (or .^ for element-wise)");
			return Base;
		}
		const double P = Exp.Scalar();
		if (P != std::floor(P) || std::fabs(P) > 1e6)
		{
			Fail("^: a matrix power must be an integer (or use expm)");
			return Base;
		}
		Matrix B = Base.M;
		if (P < 0)
		{
			if (!la::Inverse(Base.M, &B, &Err))
			{
				return Base;
			}
		}
		long long N = static_cast<long long>(std::fabs(P));
		Matrix R = Matrix::Identity(B.Rows);
		while (N > 0)
		{
			if (N & 1)
			{
				R = la::Mul(R, B);
			}
			B = la::Mul(B, B);
			N >>= 1;
		}
		return Value(R);
	}

	bool Index1(const Value& Base, const std::vector<Value>& Args, Value* Out)
	{
		const Matrix& M = Base.M;
		auto ToIdx = [this](const Value& V, int Limit, int* Idx) {
			if (!V.IsScalar() || V.Scalar() != std::floor(V.Scalar()))
			{
				Fail("index must be a whole number");
				return false;
			}
			const int I = static_cast<int>(V.Scalar());
			if (I < 1 || I > Limit)
			{
				Fail("index " + std::to_string(I) + " out of range 1.." + std::to_string(Limit));
				return false;
			}
			*Idx = I - 1;
			return true;
		};
		if (Args.size() == 1)
		{
			int I = 0;
			if (!ToIdx(Args[0], M.Count(), &I))
			{
				return false;
			}
			*Out = Value(M.V[static_cast<size_t>(I)]);
			return true;
		}
		if (Args.size() == 2)
		{
			int R = 0, C = 0;
			if (!ToIdx(Args[0], M.Rows, &R) || !ToIdx(Args[1], M.Cols, &C))
			{
				return false;
			}
			*Out = Value(M(R, C));
			return true;
		}
		Fail("index with (i) or (row, col)");
		return false;
	}

	std::vector<Value> ArgList()
	{
		std::vector<Value> Args;
		if (Accept(Tok::RParen))
		{
			return Args;
		}
		for (;;)
		{
			Args.push_back(Expression());
			if (Failed())
			{
				return Args;
			}
			if (Accept(Tok::Comma))
			{
				continue;
			}
			Expect(Tok::RParen, "')'");
			return Args;
		}
	}

	Value CallFunction(const FunctionInfo& F, std::vector<Value> Args)
	{
		const size_t Declared = F.Params.size();
		if (Args.size() < F.RequiredCount() || (!F.Variadic && Args.size() > Declared))
		{
			Fail(F.Name + " takes " + F.Signature());
			return Value();
		}
		for (size_t i = 0; i < Args.size() && i < Declared; ++i)
		{
			if (F.Params[i].Kind == ParamKind::Number && !Args[i].IsScalar())
			{
				Fail(F.Name + ": " + F.Params[i].Name + " must be a number, got " + Shape(Args[i].M));
				return Value();
			}
		}
		for (size_t i = Args.size(); i < Declared; ++i)
		{
			Args.push_back(Value(F.Params[i].Default));
		}
		if ((F.Angle & AngleArgs) != 0 && !Args.empty())
		{
			for (double& X : Args[0].M.V)
			{
				X = S.ToRadians(X);
			}
		}
		Result R = F.Fn(Call(Args));
		if (!R.Ok)
		{
			Fail(F.Name + ": " + R.Error);
			return Value();
		}
		if ((F.Angle & AngleResult) != 0)
		{
			for (double& X : R.Val.M.V)
			{
				X = S.FromRadians(X);
			}
		}
		return R.Val;
	}

	Value Postfix()
	{
		Value V = Primary();
		while (!Failed())
		{
			if (Accept(Tok::Quote))
			{
				V = Value(la::Transpose(V.M));
			}
			else if (Accept(Tok::Bang))
			{
				for (double& X : V.M.V)
				{
					X = std::tgamma(X + 1.0);
				}
				V.Fields.clear();
			}
			else if (Peek().Kind == Tok::Dot && Peek(1).Kind == Tok::Ident)
			{
				Index += 1;
				const std::string Field = Peek().Text;
				Index += 1;
				const int F = V.FieldIndex(Field);
				if (F < 0)
				{
					std::string Known;
					for (const std::string& Name : V.Fields)
					{
						Known += (Known.empty() ? "" : ", ") + Name;
					}
					Fail("no field '" + Field + "'" + (Known.empty() ? std::string(" (not a record)") : " (fields: " + Known + ")"));
					break;
				}
				V = Value(V.M.V[static_cast<size_t>(F)]);
			}
			else if (Peek().Kind == Tok::LParen)
			{
				++Index;
				std::vector<Value> Args = ArgList();
				if (Failed())
				{
					break;
				}
				Value Picked;
				if (!Index1(V, Args, &Picked))
				{
					break;
				}
				V = Picked;
			}
			else
			{
				break;
			}
		}
		return V;
	}

	bool Concat(std::vector<std::vector<Matrix>>& Rows, Matrix* Out)
	{
		std::vector<Matrix> Bands;
		for (std::vector<Matrix>& Row : Rows)
		{
			int H = -1, W = 0;
			for (const Matrix& M : Row)
			{
				if (M.IsEmpty())
				{
					continue;
				}
				if (H >= 0 && M.Rows != H)
				{
					Fail("[ , ]: elements in a row need the same height");
					return false;
				}
				H = M.Rows;
				W += M.Cols;
			}
			if (H < 0)
			{
				continue;
			}
			Matrix Band(H, W);
			int Col = 0;
			for (const Matrix& M : Row)
			{
				for (int i = 0; i < M.Rows; ++i)
				{
					for (int j = 0; j < M.Cols; ++j)
					{
						Band(i, Col + j) = M(i, j);
					}
				}
				Col += M.Cols;
			}
			Bands.push_back(Band);
		}
		int H = 0, W = -1;
		for (const Matrix& B : Bands)
		{
			if (W >= 0 && B.Cols != W)
			{
				Fail("[ ; ]: rows need the same width");
				return false;
			}
			W = B.Cols;
			H += B.Rows;
		}
		Matrix M(H, W < 0 ? 0 : W);
		int Row = 0;
		for (const Matrix& B : Bands)
		{
			for (int i = 0; i < B.Rows; ++i)
			{
				for (int j = 0; j < B.Cols; ++j)
				{
					M(Row + i, j) = B(i, j);
				}
			}
			Row += B.Rows;
		}
		*Out = M;
		return true;
	}

	Value Primary()
	{
		const Token Tk = Peek();
		switch (Tk.Kind)
		{
		case Tok::Number:
			++Index;
			return Value(Tk.Number);
		case Tok::LParen:
		{
			++Index;
			Value V = Expression();
			Expect(Tok::RParen, "')'");
			return V;
		}
		case Tok::LBracket:
		{
			++Index;
			std::vector<std::vector<Matrix>> Rows(1);
			if (Accept(Tok::RBracket))
			{
				return Value(Matrix());
			}
			for (;;)
			{
				Rows.back().push_back(Expression().M);
				if (Failed())
				{
					return Value();
				}
				if (Accept(Tok::Comma))
				{
					continue;
				}
				if (Accept(Tok::Semicolon))
				{
					Rows.emplace_back();
					continue;
				}
				if (!Expect(Tok::RBracket, "']'"))
				{
					return Value();
				}
				break;
			}
			Matrix M;
			if (!Concat(Rows, &M))
			{
				return Value();
			}
			return Value(M);
		}
		case Tok::Ident:
		{
			++Index;
			Value Var;
			if (S.GetVariable(Tk.Text, &Var))
			{
				return Var;   // Postfix handles x(i) indexing
			}
			if (const FunctionInfo* F = S.Functions().Find(Tk.Text))
			{
				if (!Accept(Tok::LParen))
				{
					if (F->RequiredCount() == 0 && !F->Variadic)
					{
						return CallFunction(*F, {});
					}
					Fail(F->Name + " is a function: " + F->Signature());
					return Value();
				}
				std::vector<Value> Args = ArgList();
				if (Failed())
				{
					return Value();
				}
				return CallFunction(*F, Args);
			}
			double X = 0.0;
			if (Constant(Tk.Text, &X))
			{
				return Value(X);
			}
			Fail("unknown name '" + Tk.Text + "'");
			return Value();
		}
		case Tok::End:
			Fail("incomplete expression");
			return Value();
		default:
			Fail("unexpected '" + Tk.Text + "'");
			return Value();
		}
	}
};

std::string Trim(const std::string& In)
{
	size_t B = 0, E = In.size();
	while (B < E && std::isspace(static_cast<unsigned char>(In[B])) != 0)
	{
		++B;
	}
	while (E > B && std::isspace(static_cast<unsigned char>(In[E - 1])) != 0)
	{
		--E;
	}
	return In.substr(B, E - B);
}

// Splits at ';' outside brackets and parentheses.
std::vector<std::string> SplitStatements(const std::string& Line)
{
	std::vector<std::string> Out;
	int Depth = 0;
	size_t Start = 0;
	for (size_t i = 0; i < Line.size(); ++i)
	{
		const char Ch = Line[i];
		if (Ch == '(' || Ch == '[')
		{
			++Depth;
		}
		else if ((Ch == ')' || Ch == ']') && Depth > 0)
		{
			--Depth;
		}
		else if (Ch == ';' && Depth == 0)
		{
			Out.push_back(Line.substr(Start, i - Start));
			Start = i + 1;
		}
	}
	Out.push_back(Line.substr(Start));
	return Out;
}

std::string PadRight(const std::string& S, size_t W)
{
	return S.size() >= W ? S : S + std::string(W - S.size(), ' ');
}

std::string PadLeft(const std::string& S, size_t W)
{
	return S.size() >= W ? S : std::string(W - S.size(), ' ') + S;
}
} // namespace session_detail

using namespace session_detail;

std::vector<std::string> WrapText(const std::string& Text, int Width)
{
	std::vector<std::string> Out;
	size_t Start = 0;
	while (Start <= Text.size())
	{
		size_t Nl = Text.find('\n', Start);
		std::string Para = Text.substr(Start, Nl == std::string::npos ? std::string::npos : Nl - Start);
		if (Width <= 0)
		{
			Out.push_back(Para);
		}
		else
		{
			const size_t W = static_cast<size_t>(Width);
			while (Para.size() > W)
			{
				size_t Cut = Para.rfind(' ', W);
				if (Cut == std::string::npos || Cut == 0)
				{
					Cut = W;
				}
				Out.push_back(Para.substr(0, Cut));
				Para = Para.substr(Cut);
				if (!Para.empty() && Para[0] == ' ')
				{
					Para.erase(0, 1);
				}
			}
			Out.push_back(Para);
		}
		if (Nl == std::string::npos)
		{
			break;
		}
		Start = Nl + 1;
	}
	return Out;
}

Session::Session() : Session(Registry::Default()) {}

Session::Session(const Registry& InRegistry) : Reg(&InRegistry) {}

double Session::ToRadians(double X) const
{
	return Angles == AngleMode::Degrees ? X * cst::Pi / 180.0 : X;
}

double Session::FromRadians(double X) const
{
	return Angles == AngleMode::Degrees ? X * 180.0 / cst::Pi : X;
}

bool Session::GetVariable(const std::string& Name, Value* Out) const
{
	auto It = Vars.find(Name);
	if (It == Vars.end())
	{
		return false;
	}
	*Out = It->second;
	return true;
}

void Session::SetVariable(const std::string& Name, const Value& V)
{
	Vars[Name] = V;
}

std::string Session::FormatNumber(double X) const
{
	if (std::isnan(X))
	{
		return "nan";
	}
	if (std::isinf(X))
	{
		return X > 0 ? "inf" : "-inf";
	}
	if (X == 0.0)
	{
		return "0";
	}
	char Buf[64];
	std::snprintf(Buf, sizeof(Buf), "%.*g", Digits, X);
	std::string S = Buf;
	// 1.5e+06 -> 1.5e6, 2e-07 -> 2e-7
	const size_t E = S.find('e');
	if (E != std::string::npos)
	{
		std::string Mant = S.substr(0, E);
		std::string Exp = S.substr(E + 1);
		std::string Sign;
		if (!Exp.empty() && (Exp[0] == '+' || Exp[0] == '-'))
		{
			Sign = Exp[0] == '-' ? "-" : "";
			Exp.erase(0, 1);
		}
		while (Exp.size() > 1 && Exp[0] == '0')
		{
			Exp.erase(0, 1);
		}
		S = Mant + "e" + Sign + Exp;
	}
	return S;
}

std::vector<std::string> Session::Format(const Value& V, int Width) const
{
	std::vector<std::string> Lines;
	const Matrix& M = V.M;
	if (V.IsRecord())
	{
		size_t NameW = 0;
		for (const std::string& F : V.Fields)
		{
			NameW = std::max(NameW, F.size());
		}
		for (size_t i = 0; i < V.Fields.size() && i < M.V.size(); ++i)
		{
			Lines.push_back(PadRight(V.Fields[i], NameW) + " = " + FormatNumber(M.V[i]));
		}
		return Lines;
	}
	if (M.IsEmpty())
	{
		Lines.push_back("[]");
		return Lines;
	}
	if (M.IsScalar())
	{
		Lines.push_back(FormatNumber(M.V[0]));
		return Lines;
	}
	std::vector<size_t> ColW(static_cast<size_t>(M.Cols), 0);
	std::vector<std::string> Cells(M.V.size());
	for (int i = 0; i < M.Rows; ++i)
	{
		for (int j = 0; j < M.Cols; ++j)
		{
			const size_t Idx = static_cast<size_t>(i * M.Cols + j);
			Cells[Idx] = FormatNumber(M(i, j));
			ColW[static_cast<size_t>(j)] = std::max(ColW[static_cast<size_t>(j)], Cells[Idx].size());
		}
	}
	for (int i = 0; i < M.Rows; ++i)
	{
		std::string Row = "[";
		for (int j = 0; j < M.Cols; ++j)
		{
			Row += " " + PadLeft(Cells[static_cast<size_t>(i * M.Cols + j)], ColW[static_cast<size_t>(j)]);
		}
		Row += " ]";
		if (Width > 0 && static_cast<int>(Row.size()) > Width)
		{
			// Too wide for the screen: one element per line instead.
			Lines.clear();
			for (int r = 0; r < M.Rows; ++r)
			{
				for (int c = 0; c < M.Cols; ++c)
				{
					Lines.push_back("(" + std::to_string(r + 1) + "," + std::to_string(c + 1) + ") " + Cells[static_cast<size_t>(r * M.Cols + c)]);
				}
			}
			return Lines;
		}
		Lines.push_back(Row);
	}
	return Lines;
}

std::vector<std::string> Session::Complete(const std::string& Prefix) const
{
	std::vector<std::string> Out = Reg->Completions(Prefix);
	for (const auto& KV : Vars)
	{
		if (KV.first.compare(0, Prefix.size(), Prefix) == 0)
		{
			Out.push_back(KV.first);
		}
	}
	for (const char* Name : kConstantNames)
	{
		if (std::string(Name).compare(0, Prefix.size(), Prefix) == 0)
		{
			Out.push_back(Name);
		}
	}
	for (const char* Name : kCommandNames)
	{
		if (std::string(Name).compare(0, Prefix.size(), Prefix) == 0)
		{
			Out.push_back(Name);
		}
	}
	std::sort(Out.begin(), Out.end());
	Out.erase(std::unique(Out.begin(), Out.end()), Out.end());
	return Out;
}

bool Session::RunCommand(const std::string& Line, Outcome* Out)
{
	const size_t Sp = Line.find(' ');
	const std::string Word = Line.substr(0, Sp);
	const std::string Arg = Sp == std::string::npos ? std::string() : Trim(Line.substr(Sp + 1));
	auto Say = [Out](const std::string& Text) {
		Out->Ok = true;
		Out->HasValue = false;
		Out->Text = Text;
	};
	if (Word == "deg" && Arg.empty())
	{
		Angles = AngleMode::Degrees;
		Say("angles in degrees");
		return true;
	}
	if (Word == "rad" && Arg.empty())
	{
		Angles = AngleMode::Radians;
		Say("angles in radians");
		return true;
	}
	if (Word == "clear" && Arg.empty())
	{
		Vars.clear();
		Say("variables cleared");
		return true;
	}
	if (Word == "vars" && Arg.empty())
	{
		std::string Text;
		for (const auto& KV : Vars)
		{
			Text += KV.first + " " + Shape(KV.second.M) + (KV.second.IsRecord() ? " record" : "") + "\n";
		}
		Say(Text.empty() ? "no variables" : Text.substr(0, Text.size() - 1));
		return true;
	}
	if (Word == "digits")
	{
		if (Arg.empty())
		{
			Say("digits " + std::to_string(Digits));
			return true;
		}
		const int D = std::atoi(Arg.c_str());
		SetDigits(D);
		Say("digits " + std::to_string(Digits));
		return true;
	}
	if (Word == "list")
	{
		std::string Text;
		if (Arg.empty())
		{
			for (const std::string& Cat : Reg->Categories())
			{
				Text += Cat + " (" + std::to_string(Reg->InCategory(Cat).size()) + ")\n";
			}
		}
		else
		{
			const std::vector<const FunctionInfo*> Fns = Reg->InCategory(Arg);
			if (Fns.empty())
			{
				Out->Ok = false;
				Out->Text = "no category '" + Arg + "' (type list)";
				return true;
			}
			for (const FunctionInfo* F : Fns)
			{
				Text += F->Signature() + "\n";
			}
		}
		Say(Text.empty() ? "" : Text.substr(0, Text.size() - 1));
		return true;
	}
	if (Word == "help")
	{
		if (Arg.empty())
		{
			Say("type an expression and press Enter. help <fn>, list, list <category>, vars, clear, deg, rad, digits <n>. "
				"Matrices: [1, 2; 3, 4], A', A \\ b, A .* B. Records: h = hohmann(6778e3, 42164e3); h.dv1");
			return true;
		}
		const FunctionInfo* F = Reg->Find(Arg);
		if (!F)
		{
			double X = 0.0;
			if (Constant(Arg, &X))
			{
				Say(Arg + " = " + FormatNumber(X) + " (constant)");
				return true;
			}
			Out->Ok = false;
			Out->Text = "no function '" + Arg + "'";
			return true;
		}
		std::string Text = F->Signature() + "\n" + F->Summary;
		for (const Param& P : F->Params)
		{
			Text += "\n  " + P.Name + (P.Unit.empty() ? "" : " [" + P.Unit + "]");
			if (P.HasDefault)
			{
				Text += " = " + FormatNumber(P.Default);
			}
		}
		if (!F->Returns.empty())
		{
			Text += "\n-> " + F->Returns;
		}
		if (!F->Source.empty())
		{
			Text += "\nsource: " + F->Source;
		}
		Say(Text);
		return true;
	}
	return false;
}

Outcome Session::Execute(const std::string& RawLine)
{
	Outcome Out;
	const std::string Line = Trim(RawLine);
	if (Line.empty())
	{
		Out.HasValue = false;
		return Out;
	}
	if (RunCommand(Line, &Out))
	{
		Out.Lines = WrapText(Out.Text, 0);
		return Out;
	}
	for (const std::string& RawStmt : SplitStatements(Line))
	{
		const std::string Stmt = Trim(RawStmt);
		if (Stmt.empty())
		{
			continue;
		}
		const std::vector<Token> Tokens = Lex(Stmt);
		for (const Token& Tk : Tokens)
		{
			if (Tk.Kind == Tok::Bad)
			{
				Out.Ok = false;
				Out.Text = "unexpected character '" + Tk.Text + "'";
				Out.Lines = {Out.Text};
				return Out;
			}
		}
		std::string Target;
		size_t Skip = 0;
		if (Tokens.size() > 2 && Tokens[0].Kind == Tok::Ident && Tokens[1].Kind == Tok::Equals)
		{
			Target = Tokens[0].Text;
			Skip = 2;
			double Unused = 0.0;
			if (Reg->Find(Target) || Constant(Target, &Unused) || Target == "ans")
			{
				Out.Ok = false;
				Out.Text = "cannot assign to '" + Target + "'";
				Out.Lines = {Out.Text};
				return Out;
			}
		}
		const std::vector<Token> Body(Tokens.begin() + static_cast<std::ptrdiff_t>(Skip), Tokens.end());
		Parser P(*this, Body);
		Value V = P.Expression();
		if (!P.Failed() && !P.AtEnd())
		{
			P.Fail("unexpected '" + P.Peek().Text + "'");
		}
		if (P.Failed())
		{
			Out.Ok = false;
			Out.HasValue = false;
			Out.Text = P.Error();
			Out.Lines = {Out.Text};
			return Out;
		}
		if (!Target.empty())
		{
			Vars[Target] = V;
		}
		Vars["ans"] = V;
		Out.Val = V;
		Out.HasValue = true;
	}
	if (Out.HasValue)
	{
		Out.Lines = Format(Out.Val);
		Out.Text.clear();
		for (size_t i = 0; i < Out.Lines.size(); ++i)
		{
			Out.Text += (i ? "\n" : "") + Out.Lines[i];
		}
	}
	return Out;
}
} // namespace mcalc
