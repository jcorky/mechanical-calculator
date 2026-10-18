// A minimal test framework: TEST(name) { ... } registers a test; CHECK and CHECK_NEAR record failures.
#pragma once

#include "GearReal.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace mtest
{
struct Case
{
	const char* Name;
	void (*Fn)();
};

inline std::vector<Case>& Cases()
{
	static std::vector<Case> All;
	return All;
}

inline int& Failures()
{
	static int N = 0;
	return N;
}

inline int& Checks()
{
	static int N = 0;
	return N;
}

struct Registrar
{
	Registrar(const char* Name, void (*Fn)()) { Cases().push_back({Name, Fn}); }
};

inline void Fail(const char* File, int Line, const std::string& What)
{
	++Failures();
	std::printf("  FAIL %s:%d: %s\n", File, Line, What.c_str());
}

// Values read off the gear machine compare as doubles in CHECK_NEAR.
inline double AsDouble(double X)
{
	return X;
}
inline double AsDouble(const mcalc::Real& X)
{
	return X.ToDouble();
}

inline bool Near(double A, double B, double Tol)
{
	if (std::isnan(A) && std::isnan(B))
	{
		return true;
	}
	if (std::isinf(A) || std::isinf(B))
	{
		return A == B;
	}
	return std::fabs(A - B) <= Tol;
}
} // namespace mtest

#define TEST(Name)                                                   \
	static void Name();                                              \
	static mtest::Registrar Name##_registrar(#Name, &Name);          \
	static void Name()

#define CHECK(Cond)                                                  \
	do                                                               \
	{                                                                \
		++mtest::Checks();                                           \
		if (!(Cond))                                                 \
			mtest::Fail(__FILE__, __LINE__, #Cond);                  \
	} while (0)

#define CHECK_NEAR(A, B, Tol)                                        \
	do                                                               \
	{                                                                \
		++mtest::Checks();                                           \
		const double CheckA_ = mtest::AsDouble(A), CheckB_ = mtest::AsDouble(B);                   \
		if (!mtest::Near(CheckA_, CheckB_, mtest::AsDouble(Tol)))                   \
		{                                                            \
			char Buf_[256];                                          \
			std::snprintf(Buf_, sizeof(Buf_), "%s = %.17g, expected %.17g (tol %g)", #A, CheckA_, CheckB_, mtest::AsDouble(Tol)); \
			mtest::Fail(__FILE__, __LINE__, Buf_);                   \
		}                                                            \
	} while (0)
