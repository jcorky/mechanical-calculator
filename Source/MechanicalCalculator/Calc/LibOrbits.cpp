// Orbits and inertial navigation: the orbits/ bench (Kepler, Hohmann, Lambert, conics, Clohessy-Wiltshire,
// rendezvous, integrator energy drift, bi-elliptic transfers) and the ins/ bench (Allan deviation, Schuler loop,
// coning, 1-D INS/GNSS Kalman filter). Ported from orbits/orbits.js and ins/ins.js; plane vectors are [x, y].
#include "CalcLibrary.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace mcalc
{
namespace orbits_detail
{
constexpr double OrbInf = std::numeric_limits<double>::infinity();

struct OVec
{
	double X = 0.0;
	double Y = 0.0;
};

OVec OV(double X, double Y)
{
	OVec R;
	R.X = X;
	R.Y = Y;
	return R;
}
double ONorm(OVec A) { return std::hypot(A.X, A.Y); }
double ODot(OVec A, OVec B) { return A.X * B.X + A.Y * B.Y; }
OVec OAcc(OVec R, double Mu)
{
	const double D = ONorm(R), K = -Mu / (D * D * D);
	return OV(K * R.X, K * R.Y);
}
double OEnergy(OVec R, OVec V, double Mu) { return 0.5 * ODot(V, V) - Mu / ONorm(R); }

bool OGetVec2(const Matrix& M, OVec* Out)
{
	if (M.Count() != 2)
	{
		return false;
	}
	Out->X = M.V[0];
	Out->Y = M.V[1];
	return true;
}

// Classical fourth-order Runge-Kutta step of r'' = -mu r / |r|^3 (orbits.js rk4, same operation order).
void ORk4(OVec& R, OVec& V, double Dt, double Mu)
{
	const OVec A1 = OAcc(R, Mu);
	const OVec R2 = OV(R.X + Dt / 2 * V.X, R.Y + Dt / 2 * V.Y), V2 = OV(V.X + Dt / 2 * A1.X, V.Y + Dt / 2 * A1.Y), A2 = OAcc(R2, Mu);
	const OVec R3 = OV(R.X + Dt / 2 * V2.X, R.Y + Dt / 2 * V2.Y), V3 = OV(V.X + Dt / 2 * A2.X, V.Y + Dt / 2 * A2.Y), A3 = OAcc(R3, Mu);
	const OVec R4 = OV(R.X + Dt * V3.X, R.Y + Dt * V3.Y), V4 = OV(V.X + Dt * A3.X, V.Y + Dt * A3.Y), A4 = OAcc(R4, Mu);
	const OVec RN = OV(R.X + Dt / 6 * (V.X + 2 * V2.X + 2 * V3.X + V4.X), R.Y + Dt / 6 * (V.Y + 2 * V2.Y + 2 * V3.Y + V4.Y));
	const OVec VN = OV(V.X + Dt / 6 * (A1.X + 2 * A2.X + 2 * A3.X + A4.X), V.Y + Dt / 6 * (A1.Y + 2 * A2.Y + 2 * A3.Y + A4.Y));
	R = RN;
	V = VN;
}

struct OPropOut
{
	OVec R;
	OVec V;
	int Steps = 0;
	double RMax = 0.0;   // extremes of |r| over the samples orbits.js keeps (every `Every` steps plus both ends)
	double RMin = 0.0;
	bool Ok = true;
};

// Integrate for time T with step (local circular period)/N, the last step clipped to land on T.
OPropOut OPropagate(OVec R, OVec V, double T, double N, double Mu, int Every)
{
	OPropOut Out;
	double Done = 0.0;
	int K = 0;
	Out.RMax = Out.RMin = ONorm(R);
	while (Done < T)
	{
		const double Dt = std::min(T - Done, cst::Tau * std::sqrt(std::pow(ONorm(R), 3) / Mu) / N);
		ORk4(R, V, Dt, Mu);
		Done += Dt;
		++K;
		if (Every > 0 && K % Every == 0)
		{
			Out.RMax = std::max(Out.RMax, ONorm(R));
			Out.RMin = std::min(Out.RMin, ONorm(R));
		}
		if (K > 300000 || !std::isfinite(Dt) || Dt <= 0.0)
		{
			Out.Ok = false;
			break;
		}
	}
	Out.RMax = std::max(Out.RMax, ONorm(R));
	Out.RMin = std::min(Out.RMin, ONorm(R));
	Out.R = R;
	Out.V = V;
	Out.Steps = K;
	return Out;
}

double KeplerSolve(double M, double E0, int* Iters)
{
	double E = E0 < 0.8 ? M : cst::Pi;
	for (int i = 0; i < 50; ++i)
	{
		const double D = (E - E0 * std::sin(E) - M) / (1 - E0 * std::cos(E));
		E -= D;
		if (std::fabs(D) < 1e-14)
		{
			*Iters = i + 1;
			return E;
		}
	}
	*Iters = 50;
	return E;
}

// True anomaly from eccentric anomaly, continuous in E (same revolution).
double TrueFromEcc(double E, double Ecc)
{
	const double Beta = Ecc / (1.0 + std::sqrt(1.0 - Ecc * Ecc));
	return E + 2.0 * std::atan(Beta * std::sin(E) / (1.0 - Beta * std::cos(E)));
}

double StumpC(double Z)
{
	if (Z > 1e-8)
	{
		return (1 - std::cos(std::sqrt(Z))) / Z;
	}
	if (Z < -1e-8)
	{
		return (std::cosh(std::sqrt(-Z)) - 1) / (-Z);
	}
	return 0.5 - Z / 24;
}
double StumpS(double Z)
{
	if (Z > 1e-8)
	{
		const double S = std::sqrt(Z);
		return (S - std::sin(S)) / (S * S * S);
	}
	if (Z < -1e-8)
	{
		const double S = std::sqrt(-Z);
		return (std::sinh(S) - S) / (S * S * S);
	}
	return 1.0 / 6 - Z / 120;
}

struct LambertOut
{
	OVec V1;
	OVec V2;
	double Z = 0.0;
	double Dth = 0.0;
};

struct LambertFn
{
	double R1 = 0.0, R2 = 0.0, A = 0.0, Mu = 0.0, Dt = 0.0;
	double Yz(double Z) const { return R1 + R2 + A * (Z * StumpS(Z) - 1) / std::sqrt(StumpC(Z)); }
	double Fz(double Z) const { return std::pow(Yz(Z) / StumpC(Z), 1.5) * StumpS(Z) + A * std::sqrt(Yz(Z)) - std::sqrt(Mu) * Dt; }
};

// Prograde zero-revolution Lambert by universal variables, z by bisection (orbits.js lambert).
bool LambertSolve(OVec R1v, OVec R2v, double Dt, double Mu, LambertOut* Out, std::string* Error)
{
	LambertFn L;
	L.R1 = ONorm(R1v);
	L.R2 = ONorm(R2v);
	L.Mu = Mu;
	L.Dt = Dt;
	if (!(L.R1 > 0) || !(L.R2 > 0))
	{
		*Error = "position vectors must be nonzero";
		return false;
	}
	if (!(Dt > 0))
	{
		*Error = "time of flight must be positive";
		return false;
	}
	const double Cross = R1v.X * R2v.Y - R1v.Y * R2v.X;
	double Dth = std::acos(std::max(-1.0, std::min(1.0, ODot(R1v, R2v) / (L.R1 * L.R2))));
	if (Cross < 0)
	{
		Dth = cst::Tau - Dth;
	}
	if (std::fabs(std::sin(Dth)) < 1e-9)
	{
		*Error = "transfer angle of 0 or 180 deg is degenerate";
		return false;
	}
	L.A = std::sin(Dth) * std::sqrt(L.R1 * L.R2 / (1 - std::cos(Dth)));
	double Hi = 4 * cst::Pi * cst::Pi - 1e-9, Lo = -4;
	for (int Guard = 0; Guard < 1100 && L.Yz(Lo) > 0 && L.Fz(Lo) > 0; ++Guard)
	{
		Lo *= 2;
	}
	if (L.Yz(Lo) <= 0)
	{
		double A = Lo, B = Hi;
		for (int i = 0; i < 200; ++i)
		{
			const double M = (A + B) / 2;
			if (L.Yz(M) > 0)
			{
				B = M;
			}
			else
			{
				A = M;
			}
		}
		Lo = B;
	}
	for (int i = 0; i < 200; ++i)
	{
		const double M = (Lo + Hi) / 2;
		if (L.Fz(M) < 0)
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	const double Z = (Lo + Hi) / 2, Y = L.Yz(Z), F = 1 - Y / L.R1, G = L.A * std::sqrt(Y / Mu), Gd = 1 - Y / L.R2;
	Out->V1 = OV((R2v.X - F * R1v.X) / G, (R2v.Y - F * R1v.Y) / G);
	Out->V2 = OV((Gd * R2v.X - R1v.X) / G, (Gd * R2v.Y - R1v.Y) / G);
	Out->Z = Z;
	Out->Dth = Dth;
	if (!std::isfinite(Out->V1.X) || !std::isfinite(Out->V1.Y) || !std::isfinite(Out->V2.X) || !std::isfinite(Out->V2.Y))
	{
		*Error = "no prograde zero-revolution solution for this time of flight";
		return false;
	}
	return true;
}

// Clohessy-Wiltshire state-transition matrix, state [x, y, vx, vy] (x radial out, y along-track).
Matrix CwStm4(double N, double T)
{
	const double S = std::sin(N * T), C = std::cos(N * T);
	return Matrix::FromRows({{4 - 3 * C, 0, S / N, 2 * (1 - C) / N},
		{6 * (S - N * T), 1, -2 * (1 - C) / N, (4 * S - 3 * N * T) / N},
		{3 * N * S, 0, C, 2 * S},
		{-6 * N * (1 - C), 0, -2 * S, 4 * C - 3}});
}
// Six-state version [x, y, z, vx, vy, vz]: the in-plane block plus the decoupled cross-track oscillator.
Matrix CwStm6(double N, double T)
{
	const Matrix P4 = CwStm4(N, T);
	const int Map[4] = {0, 1, 3, 4};
	Matrix P(6, 6);
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			P(Map[i], Map[j]) = P4(i, j);
		}
	}
	const double S = std::sin(N * T), C = std::cos(N * T);
	P(2, 2) = C;
	P(2, 5) = S / N;
	P(5, 2) = -N * S;
	P(5, 5) = C;
	return P;
}

bool Solve2(double M00, double M01, double M10, double M11, double B0, double B1, OVec* Out)
{
	const double D = M00 * M11 - M01 * M10, Scale = std::max(std::max(std::fabs(M00), std::fabs(M01)), std::max(std::fabs(M10), std::fabs(M11)));
	if (!(std::fabs(D) > 1e-12 * Scale * Scale) || !std::isfinite(D))
	{
		return false;
	}
	*Out = OV((B0 * M11 - M01 * B1) / D, (M00 * B1 - M10 * B0) / D);
	return true;
}

// First burn velocity that brings [x, y] to the origin in time T under CW (orbits.js cwBurn).
bool CwBurn(double N, const double* X, double T, OVec* V0)
{
	const Matrix P = CwStm4(N, T);
	return Solve2(P(0, 2), P(0, 3), P(1, 2), P(1, 3), -(P(0, 0) * X[0] + P(0, 1) * X[1]), -(P(1, 0) * X[0] + P(1, 1) * X[1]), V0);
}

double OMeanMotion(double Rref, double Mu) { return std::sqrt(Mu / (Rref * Rref * Rref)); }

// Full nonlinear relative motion about a circular reference orbit, rotating frame, fixed-step RK4.
std::array<double, 4> RelNonlinear(double Rref, const std::array<double, 4>& X, double T, int Steps, double Mu)
{
	const double N = OMeanMotion(Rref, Mu), G0 = Mu / (Rref * Rref);
	auto Fn = [N, G0, Rref, Mu](const std::array<double, 4>& S) {
		const double Rx = Rref + S[0], D = std::pow(Rx * Rx + S[1] * S[1], 1.5);
		return std::array<double, 4>{S[2], S[3], 2 * N * S[3] + N * N * S[0] - Mu * Rx / D + G0, -2 * N * S[2] + N * N * S[1] - Mu * S[1] / D};
	};
	std::array<double, 4> S = X, Tmp;
	const double H = T / Steps;
	for (int i = 0; i < Steps; ++i)
	{
		const std::array<double, 4> K1 = Fn(S);
		for (int j = 0; j < 4; ++j)
		{
			Tmp[j] = S[j] + H / 2 * K1[j];
		}
		const std::array<double, 4> K2 = Fn(Tmp);
		for (int j = 0; j < 4; ++j)
		{
			Tmp[j] = S[j] + H / 2 * K2[j];
		}
		const std::array<double, 4> K3 = Fn(Tmp);
		for (int j = 0; j < 4; ++j)
		{
			Tmp[j] = S[j] + H * K3[j];
		}
		const std::array<double, 4> K4 = Fn(Tmp);
		for (int j = 0; j < 4; ++j)
		{
			S[j] = S[j] + H / 6 * (K1[j] + 2 * K2[j] + 2 * K3[j] + K4[j]);
		}
	}
	return S;
}

bool OGetState4(const Matrix& M, std::array<double, 4>* Out)
{
	if (M.Count() != 4)
	{
		return false;
	}
	for (int i = 0; i < 4; ++i)
	{
		(*Out)[static_cast<size_t>(i)] = M.V[static_cast<size_t>(i)];
	}
	return true;
}

// Largest relative energy error per orbit of a fixed-step integrator on a Kepler orbit with a = 10 000 km.
// Method 0 = explicit Euler, 1 = leapfrog (kick-drift-kick), 2 = RK4.
std::vector<double> EnergyDrift(double Ecc, int Spo, int Orbits, int Method)
{
	const double Mu = cst::MuEarth, A = 10000e3, Rp = A * (1 - Ecc);
	OVec R = OV(Rp, 0), V = OV(0, std::sqrt(Mu * (1 + Ecc) / Rp));
	const double E0 = OEnergy(R, V, Mu), H = cst::Tau * std::sqrt(A * A * A / Mu) / Spo;
	std::vector<double> Out(static_cast<size_t>(Orbits), 0.0);
	for (int k = 0; k < Orbits; ++k)
	{
		double M = 0.0;
		for (int i = 0; i < Spo; ++i)
		{
			if (Method == 1)
			{
				OVec G = OAcc(R, Mu);
				V = OV(V.X + H / 2 * G.X, V.Y + H / 2 * G.Y);
				R = OV(R.X + H * V.X, R.Y + H * V.Y);
				G = OAcc(R, Mu);
				V = OV(V.X + H / 2 * G.X, V.Y + H / 2 * G.Y);
			}
			else if (Method == 2)
			{
				ORk4(R, V, H, Mu);
			}
			else
			{
				const OVec G = OAcc(R, Mu);
				R = OV(R.X + H * V.X, R.Y + H * V.Y);
				V = OV(V.X + H * G.X, V.Y + H * G.Y);
			}
			M = std::max(M, std::fabs(OEnergy(R, V, Mu) / E0 - 1));
		}
		Out[static_cast<size_t>(k)] = M;
	}
	return Out;
}

// Hohmann and bi-elliptic totals in units of the inner circular speed (r1 = 1).
double HohTot(double R) { return std::sqrt(2 * R / (1 + R)) - 1 + 1 / std::sqrt(R) - std::sqrt(2 / (R * (1 + R))); }
double BielTot(double R, double Rb)
{
	const double A1 = (1 + Rb) / 2, A2 = (R + Rb) / 2;
	return (std::sqrt(2 - 1 / A1) - 1) + (std::sqrt(2 / Rb - 1 / A2) - std::sqrt(2 / Rb - 1 / A1)) + (std::sqrt(2 / R - 1 / A2) - std::sqrt(1 / R));
}
double BielInf(double R) { return (std::sqrt(2.0) - 1) * (1 + 1 / std::sqrt(R)); }
double BielFar(double R) { return BielTot(R, 1e12); }
double Crossing(double (*Other)(double))
{
	double Lo = 5, Hi = 15;
	for (int i = 0; i < 100; ++i)
	{
		const double M = (Lo + Hi) / 2;
		if (HohTot(M) < Other(M))
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	return (Lo + Hi) / 2;
}
double OneTangent(double R, double Ra)
{
	const double A = (1 + Ra) / 2, Vp = std::sqrt(2 - 1 / A), Vt = Vp / R, Vr = std::sqrt(std::max(0.0, 2 / R - 1 / A - Vt * Vt));
	return (Vp - 1) + std::hypot(Vt - std::sqrt(1 / R), Vr);
}

// ---------- Results ----------

Result OrbPeriod(const Call& C)
{
	const double A = C.N(0), Mu = C.N(1);
	if (!(A > 0) || !(Mu > 0))
	{
		return Result::Fail("needs a > 0 and mu > 0");
	}
	return Result::Of(cst::Tau * std::sqrt(A * A * A / Mu));
}
Result VCirc(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("needs r > 0 and mu > 0");
	}
	return Result::Of(std::sqrt(C.N(1) / C.N(0)));
}
Result VEsc(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("needs r > 0 and mu > 0");
	}
	return Result::Of(std::sqrt(2 * C.N(1) / C.N(0)));
}
Result VisViva(const Call& C)
{
	const double R = C.N(0), A = C.N(1), Mu = C.N(2);
	if (!(R > 0) || A == 0 || !(Mu > 0))
	{
		return Result::Fail("needs r > 0, a != 0 (negative for a hyperbola) and mu > 0");
	}
	const double V2 = Mu * (2 / R - 1 / A);
	if (V2 < 0)
	{
		return Result::Fail("r lies beyond apoapsis 2a of this ellipse");
	}
	return Result::Of(std::sqrt(V2));
}
Result MeanMotionFn(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("needs a > 0 and mu > 0");
	}
	return Result::Of(OMeanMotion(C.N(0), C.N(1)));
}
Result OrbEnergy(const Call& C)
{
	OVec R, V;
	if (!OGetVec2(C.M(0), &R) || !OGetVec2(C.M(1), &V))
	{
		return Result::Fail("r and v must be [x, y] vectors");
	}
	if (!(ONorm(R) > 0))
	{
		return Result::Fail("r must be nonzero");
	}
	return Result::Of(OEnergy(R, V, C.N(2)));
}
Result StumpCFn(const Call& C) { return Result::Of(StumpC(C.N(0))); }
Result StumpSFn(const Call& C) { return Result::Of(StumpS(C.N(0))); }

Result Kepler(const Call& C)
{
	const double M = C.N(0), Ecc = C.N(1);
	if (!(Ecc >= 0 && Ecc < 1))
	{
		return Result::Fail("eccentricity must be in [0, 1)");
	}
	int Iters = 0;
	const double E = KeplerSolve(M, Ecc, &Iters);
	return Result::Record({{"E", E}, {"nu", TrueFromEcc(E, Ecc)}, {"iters", static_cast<double>(Iters)}});
}
Result TrueAnom(const Call& C)
{
	const double Ecc = C.N(1);
	if (!(Ecc >= 0 && Ecc < 1))
	{
		return Result::Fail("eccentricity must be in [0, 1)");
	}
	return Result::Of(TrueFromEcc(C.N(0), Ecc));
}
Result KepState(const Call& C)
{
	const double A = C.N(0), Ecc = C.N(1), T = C.N(2), Mu = C.N(3);
	if (!(A > 0) || !(Mu > 0))
	{
		return Result::Fail("needs a > 0 and mu > 0");
	}
	if (!(Ecc >= 0 && Ecc < 1))
	{
		return Result::Fail("eccentricity must be in [0, 1)");
	}
	const double N = std::sqrt(Mu / (A * A * A)), M = std::fmod(std::fmod(N * T, cst::Tau) + cst::Tau, cst::Tau);
	int Iters = 0;
	const double E = KeplerSolve(M, Ecc, &Iters), B = A * std::sqrt(1 - Ecc * Ecc), Ed = N / (1 - Ecc * std::cos(E));
	const double X = A * (std::cos(E) - Ecc), Y = B * std::sin(E), Vx = -A * std::sin(E) * Ed, Vy = B * std::cos(E) * Ed;
	return Result::Record({{"x", X}, {"y", Y}, {"vx", Vx}, {"vy", Vy}, {"r", std::hypot(X, Y)}, {"nu", TrueFromEcc(E, Ecc)}, {"E", E}, {"M", M},
		{"period", cst::Tau / N}, {"iters", static_cast<double>(Iters)}});
}
Result TwoBody(const Call& C)
{
	OVec R, V;
	if (!OGetVec2(C.M(0), &R) || !OGetVec2(C.M(1), &V))
	{
		return Result::Fail("r and v must be [x, y] vectors");
	}
	const double T = C.N(2), N = C.N(3), Mu = C.N(4);
	if (!(ONorm(R) > 0) || !(T >= 0) || !(N >= 1) || !(Mu > 0))
	{
		return Result::Fail("needs r != 0, t >= 0, n >= 1 and mu > 0");
	}
	const OPropOut P = OPropagate(R, V, T, N, Mu, 0);
	if (!P.Ok)
	{
		return Result::Fail("more than 300000 steps; shorten t or lower n");
	}
	return Result::Record({{"x", P.R.X}, {"y", P.R.Y}, {"vx", P.V.X}, {"vy", P.V.Y}, {"steps", static_cast<double>(P.Steps)}});
}

Result Hohmann(const Call& C)
{
	const double R1 = C.N(0), R2 = C.N(1), Mu = C.N(2);
	if (!(R1 > 0) || !(R2 > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r1, r2 > 0 and mu > 0");
	}
	const double At = (R1 + R2) / 2, V1 = std::sqrt(Mu / R1), V2 = std::sqrt(Mu / R2), Vp = std::sqrt(Mu * (2 / R1 - 1 / At)), Va = std::sqrt(Mu * (2 / R2 - 1 / At));
	const double Dv1 = std::fabs(Vp - V1), Dv2 = std::fabs(V2 - Va);
	return Result::Record({{"dv1", Dv1}, {"dv2", Dv2}, {"dv", Dv1 + Dv2}, {"tof", cst::Pi * std::sqrt(At * At * At / Mu)}, {"vp", Vp}, {"va", Va}, {"v1", V1}, {"v2", V2}, {"at", At},
		{"et", std::fabs(R2 - R1) / (R1 + R2)}});
}
Result HohProp(const Call& C)
{
	const double R1 = C.N(0), R2 = C.N(1), Isp = C.N(2), Mu = C.N(3);
	if (!(R1 > 0) || !(R2 > 0) || !(Isp > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r1, r2, isp and mu > 0");
	}
	const double At = (R1 + R2) / 2, Vp = std::sqrt(Mu * (2 / R1 - 1 / At)), Va = std::sqrt(Mu * (2 / R2 - 1 / At));
	const double Dv = std::fabs(Vp - std::sqrt(Mu / R1)) + std::fabs(std::sqrt(Mu / R2) - Va), Mr = std::exp(Dv / (cst::G0 * Isp));
	return Result::Record({{"dv", Dv}, {"mr", Mr}, {"frac", 1 - 1 / Mr}});
}
Result HohFly(const Call& C)
{
	const double R1 = C.N(0), R2 = C.N(1), Mu = C.N(2);
	if (!(R1 > 0) || !(R2 > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r1, r2 > 0 and mu > 0");
	}
	const double At = (R1 + R2) / 2, Vp = std::sqrt(Mu * (2 / R1 - 1 / At)), Tof = cst::Pi * std::sqrt(At * At * At / Mu);
	const OPropOut P = OPropagate(OV(R1, 0), OV(0, Vp), Tof, 2000, Mu, 4);
	const double REnd = ONorm(P.R);
	return Result::Record({{"rend", REnd}, {"vr", ODot(P.R, P.V) / REnd}, {"speed", ONorm(P.V)}, {"rext", R2 >= R1 ? P.RMax : P.RMin}, {"steps", static_cast<double>(P.Steps)}});
}

Result Lambert(const Call& C)
{
	OVec R1, R2;
	if (!OGetVec2(C.M(0), &R1) || !OGetVec2(C.M(1), &R2))
	{
		return Result::Fail("r1 and r2 must be [x, y] position vectors");
	}
	if (!(C.N(3) > 0))
	{
		return Result::Fail("mu must be positive");
	}
	LambertOut L;
	std::string Error;
	if (!LambertSolve(R1, R2, C.N(2), C.N(3), &L, &Error))
	{
		return Result::Fail(Error);
	}
	return Result::Record({{"v1x", L.V1.X}, {"v1y", L.V1.Y}, {"v2x", L.V2.X}, {"v2y", L.V2.Y}, {"v1", ONorm(L.V1)}, {"v2", ONorm(L.V2)}, {"z", L.Z}, {"dth", L.Dth}});
}

Result Conic(const Call& C)
{
	OVec R, V;
	if (!OGetVec2(C.M(0), &R) || !OGetVec2(C.M(1), &V))
	{
		return Result::Fail("r and v must be [x, y] vectors");
	}
	const double Mu = C.N(2), Rn = ONorm(R);
	if (!(Rn > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r != 0 and mu > 0");
	}
	const double En = OEnergy(R, V, Mu), H = R.X * V.Y - R.Y * V.X, Rv = ODot(R, V), Vv = ODot(V, V);
	const double Ex = ((Vv - Mu / Rn) * R.X - Rv * V.X) / Mu, Ey = ((Vv - Mu / Rn) * R.Y - Rv * V.Y) / Mu, E = std::hypot(Ex, Ey);
	return Result::Record({{"a", -Mu / (2 * En)}, {"e", E}, {"ex", Ex}, {"ey", Ey}, {"h", H}, {"p", H * H / Mu}, {"rp", H * H / Mu / (1 + E)}, {"energy", En}});
}

Result CwStmFn(const Call& C)
{
	const double N = C.N(0), T = C.N(1);
	const int Dim = C.Int(2);
	if (!(N > 0))
	{
		return Result::Fail("mean motion n must be positive");
	}
	if (Dim == 4)
	{
		return Result::Of(CwStm4(N, T));
	}
	if (Dim == 6)
	{
		return Result::Of(CwStm6(N, T));
	}
	return Result::Fail("dim must be 4 ([x, y, vx, vy]) or 6 ([x, y, z, vx, vy, vz])");
}
Result CwDrift(const Call& C)
{
	const double N = C.N(0), T = C.N(1);
	const Matrix& X = C.M(2);
	if (!(N > 0))
	{
		return Result::Fail("mean motion n must be positive");
	}
	if (X.Count() == 4)
	{
		const Matrix S = la::Mul(CwStm4(N, T), Matrix::Column(X.V));
		return Result::Record({{"x", S.V[0]}, {"y", S.V[1]}, {"vx", S.V[2]}, {"vy", S.V[3]}, {"range", std::hypot(S.V[0], S.V[1])}});
	}
	if (X.Count() == 6)
	{
		const Matrix S = la::Mul(CwStm6(N, T), Matrix::Column(X.V));
		return Result::Record({{"x", S.V[0]}, {"y", S.V[1]}, {"z", S.V[2]}, {"vx", S.V[3]}, {"vy", S.V[4]}, {"vz", S.V[5]}, {"range", std::hypot(S.V[0], S.V[1], S.V[2])}});
	}
	return Result::Fail("state must be [x, y, vx, vy] or [x, y, z, vx, vy, vz]");
}
Result CwNonlin(const Call& C)
{
	std::array<double, 4> X;
	const double Rref = C.N(0), T = C.N(2), Mu = C.N(4);
	const int Steps = C.Int(3);
	if (!OGetState4(C.M(1), &X))
	{
		return Result::Fail("state must be [x, y, vx, vy]");
	}
	if (!(Rref > 0) || !(Mu > 0) || Steps < 1 || Steps > 200000)
	{
		return Result::Fail("needs rref > 0, mu > 0 and 1 <= steps <= 200000");
	}
	const std::array<double, 4> S = RelNonlinear(Rref, X, T, Steps, Mu);
	return Result::Record({{"x", S[0]}, {"y", S[1]}, {"vx", S[2]}, {"vy", S[3]}, {"range", std::hypot(S[0], S[1])}});
}
Result Rendezvous(const Call& C)
{
	std::array<double, 4> X;
	const double Rref = C.N(0), T = C.N(2), Mu = C.N(3);
	if (!OGetState4(C.M(1), &X))
	{
		return Result::Fail("chaser state must be [x, y, vx, vy]");
	}
	if (!(Rref > 0) || !(Mu > 0) || !(T > 0))
	{
		return Result::Fail("needs rref > 0, T > 0 and mu > 0");
	}
	const double N = OMeanMotion(Rref, Mu);
	OVec V0;
	if (!CwBurn(N, X.data(), T, &V0))
	{
		return Result::Fail("T is a whole number of orbits: the transfer is singular");
	}
	const Matrix P = CwStm4(N, T), Xa = la::Mul(P, Matrix::Column({X[0], X[1], V0.X, V0.Y}));
	const OVec Dv1 = OV(V0.X - X[2], V0.Y - X[3]), Dv2 = OV(-Xa.V[2], -Xa.V[3]);
	return Result::Record({{"v0x", V0.X}, {"v0y", V0.Y}, {"dv1x", Dv1.X}, {"dv1y", Dv1.Y}, {"dv2x", Dv2.X}, {"dv2y", Dv2.Y}, {"dv1", ONorm(Dv1)}, {"dv2", ONorm(Dv2)},
		{"dv", ONorm(Dv1) + ONorm(Dv2)}, {"n", N}});
}
Result CwShoot(const Call& C)
{
	std::array<double, 4> X;
	const double Rref = C.N(0), T = C.N(2), Mu = C.N(4);
	const int Steps = C.Int(3);
	if (!OGetState4(C.M(1), &X))
	{
		return Result::Fail("chaser state must be [x, y, vx, vy]");
	}
	if (!(Rref > 0) || !(Mu > 0) || !(T > 0) || Steps < 1 || Steps > 10000)
	{
		return Result::Fail("needs rref > 0, T > 0, mu > 0 and 1 <= steps <= 10000");
	}
	const double N = OMeanMotion(Rref, Mu);
	const Matrix P = CwStm4(N, T);
	OVec V, VCw;
	if (!CwBurn(N, X.data(), T, &VCw))
	{
		return Result::Fail("T is a whole number of orbits: the transfer is singular");
	}
	V = VCw;
	int Iters = 99;
	for (int i = 0; i < 20; ++i)
	{
		const std::array<double, 4> M = RelNonlinear(Rref, {X[0], X[1], V.X, V.Y}, T, Steps, Mu);
		if (std::hypot(M[0], M[1]) < 1e-2)
		{
			Iters = i;
			break;
		}
		OVec D;
		if (!Solve2(P(0, 2), P(0, 3), P(1, 2), P(1, 3), M[0], M[1], &D))
		{
			return Result::Fail("singular CW Jacobian");
		}
		V = OV(V.X - D.X, V.Y - D.Y);
	}
	return Result::Record({{"v0x", V.X}, {"v0y", V.Y}, {"iters", static_cast<double>(Iters)}, {"corr", std::hypot(V.X - VCw.X, V.Y - VCw.Y)}});
}

bool DriftArgs(const Call& C, double MaxSteps, double* Ecc, int* Spo, int* Orbits, std::string* Error)
{
	*Ecc = C.N(0);
	*Spo = C.Int(1);
	*Orbits = C.Int(2);
	if (!(*Ecc >= 0 && *Ecc < 1))
	{
		*Error = "eccentricity must be in [0, 1)";
		return false;
	}
	if (*Spo < 10 || *Orbits < 1 || static_cast<double>(*Spo) * *Orbits > MaxSteps)
	{
		*Error = "needs spo >= 10, orbits >= 1 and spo * orbits <= " + std::to_string(static_cast<int>(MaxSteps));
		return false;
	}
	return true;
}
Result EDrift(const Call& C)
{
	double Ecc = 0;
	int Spo = 0, Orbits = 0;
	std::string Error;
	if (!DriftArgs(C, 200000.0, &Ecc, &Spo, &Orbits, &Error))
	{
		return Result::Fail(Error);
	}
	const int Method = C.Int(3);
	if (Method < 0 || Method > 2)
	{
		return Result::Fail("method: 0 = Euler, 1 = leapfrog, 2 = RK4");
	}
	const std::vector<double> D = EnergyDrift(Ecc, Spo, Orbits, Method);
	return Result::Record({{"max", *std::max_element(D.begin(), D.end())}, {"first", D.front()}, {"last", D.back()}});
}
Result EDriftCmp(const Call& C)
{
	double Ecc = 0;
	int Spo = 0, Orbits = 0;
	std::string Error;
	if (!DriftArgs(C, 100000.0, &Ecc, &Spo, &Orbits, &Error))
	{
		return Result::Fail(Error);
	}
	const std::vector<double> Lf = EnergyDrift(Ecc, Spo, Orbits, 1), Rk = EnergyDrift(Ecc, Spo, Orbits, 2), Eu = EnergyDrift(Ecc, Spo, Orbits, 0), Lf2 = EnergyDrift(Ecc, 2 * Spo, Orbits, 1);
	const double MLf = *std::max_element(Lf.begin(), Lf.end());
	return Result::Record({{"lf", MLf}, {"rk4", *std::max_element(Rk.begin(), Rk.end())}, {"euler", *std::max_element(Eu.begin(), Eu.end())},
		{"q", MLf / *std::max_element(Lf2.begin(), Lf2.end())}, {"bound", MLf / Lf[0]}});
}

Result HohTotFn(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("radius ratio must be positive");
	}
	return Result::Of(HohTot(C.N(0)));
}
Result BielTotFn(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("radius ratios must be positive");
	}
	return Result::Of(BielTot(C.N(0), C.N(1)));
}
Result BielInfFn(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("radius ratio must be positive");
	}
	return Result::Of(BielInf(C.N(0)));
}
Result BielThresh(const Call&)
{
	return Result::Record({{"closed", Crossing(&BielInf)}, {"finite", Crossing(&BielFar)}});
}
Result Bielliptic(const Call& C)
{
	const double R1 = C.N(0), R2 = C.N(1), Rb = C.N(2), Mu = C.N(3);
	if (!(R1 > 0) || !(R2 > 0) || !(Mu > 0) || !(Rb >= std::max(R1, R2)))
	{
		return Result::Fail("needs r1, r2, mu > 0 and rb >= max(r1, r2)");
	}
	const double A1 = (R1 + Rb) / 2, A2 = (R2 + Rb) / 2;
	const double Dv1 = std::fabs(std::sqrt(Mu * (2 / R1 - 1 / A1)) - std::sqrt(Mu / R1));
	const double Dv2 = std::fabs(std::sqrt(Mu * (2 / Rb - 1 / A2)) - std::sqrt(Mu * (2 / Rb - 1 / A1)));
	const double Dv3 = std::fabs(std::sqrt(Mu * (2 / R2 - 1 / A2)) - std::sqrt(Mu / R2));
	return Result::Record({{"dv1", Dv1}, {"dv2", Dv2}, {"dv3", Dv3}, {"dv", Dv1 + Dv2 + Dv3}, {"tof", cst::Pi * (std::sqrt(A1 * A1 * A1 / Mu) + std::sqrt(A2 * A2 * A2 / Mu))}});
}
Result Breakeven(const Call& C)
{
	const double R = C.N(0);
	if (!(R > 1))
	{
		return Result::Fail("radius ratio must exceed 1");
	}
	double Lo = R, Hi = 1e12;
	if (BielTot(R, Hi) > HohTot(R))
	{
		return Result::Of(OrbInf);
	}
	if (BielTot(R, R * (1 + 1e-9)) < HohTot(R))
	{
		return Result::Of(R);
	}
	for (int i = 0; i < 200; ++i)
	{
		const double M = std::sqrt(Lo * Hi);
		if (BielTot(R, M) > HohTot(R))
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	return Result::Of(Lo);
}
Result HohPeak(const Call&)
{
	double Best = 0, At = 0;
	for (int i = 0; i <= 3800; ++i)
	{
		const double R = 2 + i / 100.0, V = HohTot(R);
		if (V > Best)
		{
			Best = V;
			At = R;
		}
	}
	return Result::Record({{"R", At}, {"v", Best}});
}
Result OneTangentFn(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("radius ratios must be positive");
	}
	return Result::Of(OneTangent(C.N(0), C.N(1)));
}
Result OneTanScan(const Call& C)
{
	const double R = C.N(0);
	const int N = C.Int(1);
	if (!(R > 0) || N < 2 || N > 1000000)
	{
		return Result::Fail("needs R > 0 and 2 <= n <= 1000000");
	}
	double Best = OrbInf, At = 0;
	for (int i = 0; i < N; ++i)
	{
		const double Ra = R * (1 + 4.0 * i / (N - 1)), Cost = OneTangent(R, Ra);
		if (Cost < Best)
		{
			Best = Cost;
			At = Ra;
		}
	}
	return Result::Record({{"best", Best}, {"ratio", At / R}});
}
} // namespace orbits_detail

namespace orbits_ins_detail
{
constexpr double InsG = 9.80665;
constexpr double InsRe = 6.371e6;   // ins.js uses the mean radius
constexpr double InsDt = 0.1;
constexpr int InsNs = 1 << 16;
constexpr double InsKdt = 0.1;
constexpr int InsFix = 10;
constexpr int InsNk = 3000;

double InsWs() { return std::sqrt(InsG / InsRe); }
double InsFlat() { return std::sqrt(2 * std::log(2.0) / cst::Pi); }   // IEEE Std 952: flicker floor = 0.664 B

double ModelAdev(double Tau, double N, double B, double K)
{
	const double FB = InsFlat() * B;
	return std::sqrt(N * N / Tau + FB * FB + K * K * Tau / 3);
}

// mulberry32 and Box-Muller with a cached spare, bit-for-bit as ins.js rng().
struct InsRng
{
	uint32_t A = 1;
	bool HasSpare = false;
	double Spare = 0.0;

	explicit InsRng(uint32_t Seed) : A(Seed ? Seed : 1u) {}
	double U()
	{
		A += 0x6D2B79F5u;
		uint32_t T = A;
		T = (T ^ (T >> 15)) * (T | 1u);
		T ^= T + (T ^ (T >> 7)) * (T | 61u);
		return static_cast<double>(T ^ (T >> 14)) / 4294967296.0;
	}
	double Gauss()
	{
		if (HasSpare)
		{
			HasSpare = false;
			return Spare;
		}
		double X = 0.0;
		do
		{
			X = U();
		} while (X <= 0);
		const double R = std::sqrt(-2 * std::log(X)), Th = 2 * cst::Pi * U();
		Spare = R * std::sin(Th);
		HasSpare = true;
		return R * std::cos(Th);
	}
};

// Radix-2 in-place FFT, the same butterflies and twiddle recurrence as ins.js fft().
void InsFft(std::vector<double>& Re, std::vector<double>& Im, bool Inverse)
{
	const int N = static_cast<int>(Re.size());
	for (int i = 1, j = 0; i < N; ++i)
	{
		int Bit = N >> 1;
		for (; j & Bit; Bit >>= 1)
		{
			j ^= Bit;
		}
		j ^= Bit;
		if (i < j)
		{
			std::swap(Re[static_cast<size_t>(i)], Re[static_cast<size_t>(j)]);
			std::swap(Im[static_cast<size_t>(i)], Im[static_cast<size_t>(j)]);
		}
	}
	for (int Len = 2; Len <= N; Len <<= 1)
	{
		const double Ang = (Inverse ? 2 : -2) * cst::Pi / Len, Wr = std::cos(Ang), Wi = std::sin(Ang);
		const int Half = Len >> 1;
		for (int i = 0; i < N; i += Len)
		{
			double Cr = 1, Ci = 0;
			for (int k = 0; k < Half; ++k)
			{
				const size_t A = static_cast<size_t>(i + k), B = A + static_cast<size_t>(Half);
				const double Xr = Re[B] * Cr - Im[B] * Ci, Xi = Re[B] * Ci + Im[B] * Cr;
				Re[B] = Re[A] - Xr;
				Im[B] = Im[A] - Xi;
				Re[A] += Xr;
				Im[A] += Xi;
				const double T = Cr * Wr - Ci * Wi;
				Ci = Cr * Wi + Ci * Wr;
				Cr = T;
			}
		}
	}
	if (Inverse)
	{
		for (int i = 0; i < N; ++i)
		{
			Re[static_cast<size_t>(i)] /= N;
			Im[static_cast<size_t>(i)] /= N;
		}
	}
}

std::vector<int> ClusterSizes(int N)
{
	std::vector<int> Out;
	for (int i = 0; i < 40; ++i)
	{
		const int M = static_cast<int>(std::floor(std::pow(10.0, i / 39.0 * std::log10(N / 8.0)) + 0.5));
		if (std::find(Out.begin(), Out.end(), M) == Out.end())
		{
			Out.push_back(M);
		}
	}
	return Out;
}

// Overlapping Allan deviation at each cluster size.
std::vector<double> Adev(const std::vector<double>& Y, const std::vector<int>& Ms, double Dt)
{
	const size_t N = Y.size();
	std::vector<double> Th(N + 1, 0.0), Out;
	for (size_t i = 0; i < N; ++i)
	{
		Th[i + 1] = Th[i] + Y[i] * Dt;
	}
	for (int M : Ms)
	{
		double S = 0;
		const size_t Mm = static_cast<size_t>(M), Cn = N + 1 - 2 * Mm;
		for (size_t k = 0; k < Cn; ++k)
		{
			const double D = Th[k + 2 * Mm] - 2 * Th[k + Mm] + Th[k];
			S += D * D;
		}
		const double Md = M * Dt;
		Out.push_back(std::sqrt(S / static_cast<double>(Cn) / (2 * (Md * Md))));
	}
	return Out;
}

struct AllanSet
{
	std::vector<double> W, K, F;
};

// The three unit processes for one seed: white rate, rate random walk and flicker (PSD 1/(2 pi |f|)).
AllanSet MakeAllanSet(double Seed)
{
	InsRng R(static_cast<uint32_t>(static_cast<int64_t>(std::llround(1000 + Seed))));
	AllanSet S;
	const size_t N = static_cast<size_t>(InsNs);
	S.W.resize(N);
	S.K.resize(N);
	const double Sw = 1 / std::sqrt(InsDt), Sk = std::sqrt(InsDt);
	for (size_t i = 0; i < N; ++i)
	{
		S.W[i] = Sw * R.Gauss();
	}
	double Cum = 0;
	for (size_t i = 0; i < N; ++i)
	{
		Cum += Sk * R.Gauss();
		S.K[i] = Cum;
	}
	std::vector<double> Re(N, 0.0), Im(N, 0.0);
	const double Df = 1 / (InsNs * InsDt);
	for (int k = 1; k <= InsNs / 2; ++k)
	{
		const double Sp = 1 / (2 * cst::Pi * k * Df), A = std::sqrt(InsNs * Sp / InsDt);
		if (k == InsNs / 2)
		{
			Re[static_cast<size_t>(k)] = A * R.Gauss();
			continue;
		}
		const double X = A * R.Gauss() / std::sqrt(2.0), Y = A * R.Gauss() / std::sqrt(2.0);
		Re[static_cast<size_t>(k)] = X;
		Im[static_cast<size_t>(k)] = Y;
		Re[static_cast<size_t>(InsNs - k)] = X;
		Im[static_cast<size_t>(InsNs - k)] = -Y;
	}
	InsFft(Re, Im, true);
	S.F = Re;
	return S;
}

double PolySlope(const std::vector<double>& Xs, const std::vector<double>& Ys)
{
	const double N = static_cast<double>(Xs.size());
	double Mx = 0, My = 0;
	for (size_t i = 0; i < Xs.size(); ++i)
	{
		Mx += Xs[i];
		My += Ys[i];
	}
	Mx /= N;
	My /= N;
	double Sxy = 0, Sxx = 0;
	for (size_t i = 0; i < Xs.size(); ++i)
	{
		Sxy += (Xs[i] - Mx) * (Ys[i] - My);
		Sxx += (Xs[i] - Mx) * (Xs[i] - Mx);
	}
	return Sxy / Sxx;
}

double GeoMean(const std::vector<double>& A)
{
	double S = 0;
	for (double X : A)
	{
		S += std::log(X);
	}
	return std::exp(S / static_cast<double>(A.size()));
}

double SchulerClosed(double T, double B, double Phi0, double Eps)
{
	const double Ws = InsWs(), C = B - InsG * Phi0;
	return C / (Ws * Ws) * (1 - std::cos(Ws * T)) - InsRe * Eps * (T - std::sin(Ws * T) / Ws);
}

using Q4 = std::array<double, 4>;
using V3 = std::array<double, 3>;
Q4 QMul(const Q4& A, const Q4& B)
{
	return {A[0] * B[0] - A[1] * B[1] - A[2] * B[2] - A[3] * B[3], A[0] * B[1] + A[1] * B[0] + A[2] * B[3] - A[3] * B[2], A[0] * B[2] - A[1] * B[3] + A[2] * B[0] + A[3] * B[1],
		A[0] * B[3] + A[1] * B[2] - A[2] * B[1] + A[3] * B[0]};
}
Q4 QExp(const V3& V)
{
	const double A = std::hypot(V[0], V[1], V[2]);
	if (A == 0)
	{
		return {1, 0, 0, 0};
	}
	const double S = std::sin(A / 2) / A;
	return {std::cos(A / 2), S * V[0], S * V[1], S * V[2]};
}
V3 Cross3(const V3& A, const V3& B) { return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]}; }
Q4 QTrue(double T, double Al, double W) { return {std::cos(Al / 2), std::sin(Al / 2) * std::cos(W * T), std::sin(Al / 2) * std::sin(W * T), 0}; }
V3 ConeOmega(double T, double Al, double W)
{
	const double S = std::sin(Al / 2);
	return {-W * std::sin(Al) * std::sin(W * T), W * std::sin(Al) * std::cos(W * T), -2 * W * S * S};
}
// Gyro increment over [t0, t1] by three-point Gauss-Legendre quadrature.
V3 DTheta(double T0, double T1, double Al, double W)
{
	const double Gl[3][2] = {{-std::sqrt(3.0 / 5), 5.0 / 9}, {0, 8.0 / 9}, {std::sqrt(3.0 / 5), 5.0 / 9}};
	const double Cm = (T0 + T1) / 2, H = (T1 - T0) / 2;
	V3 S = {0, 0, 0};
	for (const auto& P : Gl)
	{
		const V3 O = ConeOmega(Cm + H * P[0], Al, W);
		for (int i = 0; i < 3; ++i)
		{
			S[static_cast<size_t>(i)] += P[1] * O[static_cast<size_t>(i)] * H;
		}
	}
	return S;
}
double ErrAngle(const Q4& Qt, const Q4& Q)
{
	const Q4 Qc = {Qt[0], -Qt[1], -Qt[2], -Qt[3]};
	const Q4 E = QMul(Qc, Q);
	return 2 * std::asin(std::min(1.0, std::hypot(E[1], E[2], E[3])));
}

using M3 = std::array<double, 9>;
M3 MMul3(const M3& A, const M3& B)
{
	M3 C;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			C[static_cast<size_t>(3 * i + j)] = A[static_cast<size_t>(3 * i)] * B[static_cast<size_t>(j)] + A[static_cast<size_t>(3 * i + 1)] * B[static_cast<size_t>(3 + j)] +
				A[static_cast<size_t>(3 * i + 2)] * B[static_cast<size_t>(6 + j)];
		}
	}
	return C;
}
M3 MT3(const M3& A) { return {A[0], A[3], A[6], A[1], A[4], A[7], A[2], A[5], A[8]}; }
M3 MAdd3(const M3& A, const M3& B)
{
	M3 C;
	for (size_t i = 0; i < 9; ++i)
	{
		C[i] = A[i] + B[i];
	}
	return C;
}
M3 Fapt(const M3& F, const M3& P) { return MMul3(MMul3(F, P), MT3(F)); }

struct KfModel
{
	M3 F;
	M3 Q;
	double Sa = 0.0;
};
// ins.js kfModel: state [position, velocity, accelerometer bias], IMU at 10 Hz.
KfModel MakeKf(double Qa, double Sbw)
{
	KfModel M;
	M.Sa = Qa / std::sqrt(InsKdt);
	M.F = {1, InsKdt, -InsKdt * InsKdt / 2, 0, 1, -InsKdt, 0, 0, 1};
	const double Gm[3] = {InsKdt * InsKdt / 2, InsKdt, 0};
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			M.Q[static_cast<size_t>(3 * i + j)] = Gm[i] * Gm[j] * M.Sa * M.Sa;
		}
	}
	M.Q[8] += Sbw * Sbw * InsKdt;
	return M;
}
double OutageClosed(const M3& P, double L, double Sa, double Sbw)
{
	const double Cv[3] = {1, L, -L * L / 2};
	double S = 0;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			S += Cv[i] * P[static_cast<size_t>(3 * i + j)] * Cv[j];
		}
	}
	return std::sqrt(S + Sa * Sa * InsKdt * L * L * L / 3 + Sbw * Sbw * L * L * L * L * L / 20);
}

// ---------- Results ----------

Result AllanDev(const Call& C)
{
	if (!(C.N(0) > 0))
	{
		return Result::Fail("tau must be positive");
	}
	return Result::Of(ModelAdev(C.N(0), C.N(1), C.N(2), C.N(3)));
}
Result AllanMin(const Call& C)
{
	const double N = C.N(0), B = C.N(1), K = C.N(2);
	if (!(N > 0) || !(K > 0))
	{
		return Result::Fail("needs N > 0 and K > 0 (otherwise the minimum is at tau = 0 or infinity)");
	}
	const double Tau = std::sqrt(3.0) * N / K;
	return Result::Record({{"tau", Tau}, {"adev", ModelAdev(Tau, N, B, K)}});
}
Result FlickerQuad(const Call& C)
{
	const double X = C.N(0);
	const int N = C.Int(1);
	if (!(X > 0) || N < 2 || N % 2 != 0 || N > 1000000)
	{
		return Result::Fail("needs X > 0 and an even n in 2..1000000");
	}
	const double H = X / N;
	auto Fn = [](double Xv) { return Xv < 1e-6 ? Xv : std::pow(std::sin(Xv), 4) / (Xv * Xv * Xv); };
	double S = Fn(0) + Fn(X);
	for (int i = 1; i < N; ++i)
	{
		S += (i % 2 ? 4 : 2) * Fn(i * H);
	}
	return Result::Of(S * H / 3 + 3.0 / 8 / (2 * X * X));
}
Result AllanFit(const Call& C)
{
	const AllanSet S = MakeAllanSet(C.N(0));
	const std::vector<int> Ms = ClusterSizes(InsNs);
	const std::vector<double> Aw = Adev(S.W, Ms, InsDt), Ak = Adev(S.K, Ms, InsDt), Af = Adev(S.F, Ms, InsDt);
	const double Tt = InsDt * InsNs;
	std::vector<double> Wx, Wy, Wn, Kx, Ky, Kn, Fl;
	for (size_t i = 0; i < Ms.size(); ++i)
	{
		const double T = Ms[i] * InsDt;
		if (T <= Tt / 100)
		{
			Wx.push_back(std::log(T));
			Wy.push_back(std::log(Aw[i]));
			Wn.push_back(Aw[i] * std::sqrt(T));
		}
		if (T >= 1 && T <= Tt / 50)
		{
			Kx.push_back(std::log(T));
			Ky.push_back(std::log(Ak[i]));
			Kn.push_back(Ak[i] / std::sqrt(T / 3));
		}
		if (T >= 2 && T <= Tt / 50)
		{
			Fl.push_back(Af[i]);
		}
	}
	std::sort(Fl.begin(), Fl.end());
	const size_t Nf = Fl.size();
	const double Flat = Nf % 2 ? Fl[(Nf - 1) / 2] : (Fl[Nf / 2 - 1] + Fl[Nf / 2]) / 2;
	return Result::Record({{"wslope", PolySlope(Wx, Wy)}, {"nhat", GeoMean(Wn)}, {"kslope", PolySlope(Kx, Ky)}, {"khat", GeoMean(Kn)}, {"flat", Flat}});
}
Result AllanSim(const Call& C)
{
	const double N = C.N(0), B = C.N(1), K = C.N(2);
	const AllanSet S = MakeAllanSet(C.N(3));
	const std::vector<int> Ms = ClusterSizes(InsNs);
	std::vector<double> Y(static_cast<size_t>(InsNs));
	for (size_t i = 0; i < Y.size(); ++i)
	{
		Y[i] = N * S.W[i] + B * S.F[i] + K * S.K[i];
	}
	const std::vector<double> A = Adev(Y, Ms, InsDt);
	Matrix Out(static_cast<int>(Ms.size()), 2);
	for (size_t i = 0; i < Ms.size(); ++i)
	{
		Out(static_cast<int>(i), 0) = Ms[i] * InsDt;
		Out(static_cast<int>(i), 1) = A[i];
	}
	return Result::Of(Out);
}

Result SchulerPeriod(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("needs R > 0 and g > 0");
	}
	return Result::Of(cst::Tau * std::sqrt(C.N(0) / C.N(1)));
}
Result SchulerX(const Call& C) { return Result::Of(SchulerClosed(C.N(0), C.N(1), C.N(2), C.N(3))); }
Result NaiveDrift(const Call& C)
{
	const double T = C.N(0), B = C.N(1), Phi0 = C.N(2), Eps = C.N(3);
	return Result::Of((B - InsG * Phi0) * T * T / 2 - InsG * Eps * T * T * T / 6);
}
Result SchulerSim(const Call& C)
{
	const double B = C.N(0), Phi0 = C.N(1), Eps = C.N(2), Tend = C.N(3), H = C.N(4);
	if (!(H > 0) || !(Tend > 0) || Tend / H > 1000000)
	{
		return Result::Fail("needs tend > 0, h > 0 and at most a million steps");
	}
	const int N = static_cast<int>(std::floor(Tend / H + 0.5));
	double S0 = 0, S1 = 0, S2 = Phi0, DMax = 0;
	auto Fn = [B, Eps](double V, double P) { return V3{V, B - InsG * P, V / InsRe + Eps}; };
	std::vector<double> Acc(static_cast<size_t>(N) + 1);
	Acc[0] = B - InsG * S2;
	for (int i = 1; i <= N; ++i)
	{
		const V3 K1 = Fn(S1, S2), K2 = Fn(S1 + H / 2 * K1[1], S2 + H / 2 * K1[2]), K3 = Fn(S1 + H / 2 * K2[1], S2 + H / 2 * K2[2]), K4 = Fn(S1 + H * K3[1], S2 + H * K3[2]);
		S0 += H / 6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
		S1 += H / 6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
		S2 += H / 6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
		Acc[static_cast<size_t>(i)] = B - InsG * S2;
		DMax = std::max(DMax, std::fabs(S0 - SchulerClosed(i * H, B, Phi0, Eps)));
	}
	// Period from the zero crossings of the velocity-error rate.
	std::vector<double> Z;
	for (size_t i = 0; i + 1 < Acc.size(); ++i)
	{
		if (Acc[i] * Acc[i + 1] < 0)
		{
			Z.push_back(static_cast<double>(i) * H - Acc[i] * H / (Acc[i + 1] - Acc[i]));
		}
	}
	const double Period = Z.size() < 2 ? std::nan("") : 2 * (Z.back() - Z.front()) / static_cast<double>(Z.size() - 1);
	return Result::Record({{"x", S0}, {"dmax", DMax}, {"period", Period}});
}

Result Coning0(const Call& C)
{
	const double Al = C.N(0), W = C.N(1), H = C.N(2), X = W * H, S = std::sin(Al);
	if (X == 0)
	{
		return Result::Fail("needs W h != 0");
	}
	return Result::Of(0.5 * W * S * S * (1 - std::sin(X) / X));
}
Result Coning2(const Call& C)
{
	const double Al = C.N(0), W = C.N(1), H = C.N(2), X = W * H, S = std::sin(Al);
	return Result::Of(W * S * S * X * X * X * X / 960);
}
Result ConingRun(const Call& C)
{
	const double Al = C.N(0), W = C.N(1), Rate = C.N(2), Tend = C.N(3);
	if (!(Rate > 0) || !(Tend > 0) || Tend * Rate > 100000)
	{
		return Result::Fail("needs rate > 0, tend > 0 and at most 100000 updates");
	}
	const double H = 1 / Rate;
	const int N = static_cast<int>(std::floor(Tend * Rate + 0.5));
	if (N < 2)
	{
		return Result::Fail("needs at least two updates");
	}
	Q4 Q0 = QTrue(0, Al, W), Q2 = Q0;
	double S[7] = {0, 0, 0, 0, 0, 0, 0};
	double X0 = 0, X2 = 0;
	for (int k = 0; k < N; ++k)
	{
		const double T0 = k * H;
		const V3 A = DTheta(T0, T0 + H / 2, Al, W), B = DTheta(T0 + H / 2, T0 + H, Al, W), Cx = Cross3(A, B);
		Q0 = QMul(Q0, QExp({A[0] + B[0], A[1] + B[1], A[2] + B[2]}));
		Q2 = QMul(Q2, QExp({A[0] + B[0] + 2.0 / 3 * Cx[0], A[1] + B[1] + 2.0 / 3 * Cx[1], A[2] + B[2] + 2.0 / 3 * Cx[2]}));
		const double T = (k + 1) * H;
		const Q4 Qt = QTrue(T, Al, W);
		X0 = ErrAngle(Qt, Q0);
		X2 = ErrAngle(Qt, Q2);
		S[0]++;
		S[1] += T;
		S[2] += T * T;
		S[3] += X0;
		S[4] += T * X0;
		S[5] += X2;
		S[6] += T * X2;
	}
	const double Den = S[0] * S[2] - S[1] * S[1];
	return Result::Record({{"slope0", (S[0] * S[4] - S[1] * S[3]) / Den}, {"slope2", (S[0] * S[6] - S[1] * S[5]) / Den}, {"e0", X0}, {"e2", X2}});
}
Result Bortz(const Call& C)
{
	const double Hh = C.N(0);
	const int Steps = C.Int(1);
	if (!(Hh > 0) || Steps < 1 || Steps > 100000)
	{
		return Result::Fail("needs h > 0 and 1 <= steps <= 100000");
	}
	const V3 A = {0.3, -0.2, 0.5}, B = {1.1, 0.4, -0.7};
	V3 D1, D2;
	for (size_t i = 0; i < 3; ++i)
	{
		D1[i] = A[i] * Hh / 2 + B[i] * Hh * Hh / 8;
		D2[i] = A[i] * Hh / 2 + 3 * B[i] * Hh * Hh / 8;
	}
	auto Wt = [&A, &B](double T) { return V3{A[0] + B[0] * T, A[1] + B[1] * T, A[2] + B[2] * T}; };
	auto Rate = [](const V3& P, const V3& O) {
		const double N = std::hypot(P[0], P[1], P[2]), Cc = N > 1e-8 ? (1 - N * std::sin(N) / (2 * (1 - std::cos(N)))) / (N * N) : 1.0 / 12;
		const V3 X1 = Cross3(P, O), X2 = Cross3(P, X1);
		V3 Out;
		for (size_t i = 0; i < 3; ++i)
		{
			Out[i] = O[i] + 0.5 * X1[i] + Cc * X2[i];
		}
		return Out;
	};
	auto Add = [](const V3& U, const V3& V, double S) { return V3{U[0] + S * V[0], U[1] + S * V[1], U[2] + S * V[2]}; };
	V3 P = {0, 0, 0};
	const double Hs = Hh / Steps;
	for (int i = 0; i < Steps; ++i)
	{
		const double T = i * Hs;
		const V3 K1 = Rate(P, Wt(T)), K2 = Rate(Add(P, K1, Hs / 2), Wt(T + Hs / 2)), K3 = Rate(Add(P, K2, Hs / 2), Wt(T + Hs / 2)), K4 = Rate(Add(P, K3, Hs), Wt(T + Hs));
		for (size_t j = 0; j < 3; ++j)
		{
			P[j] = P[j] + Hs / 6 * (K1[j] + 2 * K2[j] + 2 * K3[j] + K4[j]);
		}
	}
	const V3 Cx = Cross3(D1, D2);
	V3 Rr;
	for (size_t i = 0; i < 3; ++i)
	{
		Rr[i] = P[i] - D1[i] - D2[i];
	}
	return Result::Of((Rr[0] * Cx[0] + Rr[1] * Cx[1] + Rr[2] * Cx[2]) / (Cx[0] * Cx[0] + Cx[1] * Cx[1] + Cx[2] * Cx[2]));
}

Result InsGnss(const Call& C)
{
	const double Qa = C.N(0), Sb0 = C.N(1), Sbw = C.N(2), Sg = C.N(3), L = C.N(4), Out0 = C.N(5);
	if (!(Qa >= 0) || !(Sb0 > 0) || !(Sbw >= 0) || !(Sg > 0) || !(L >= 0))
	{
		return Result::Fail("needs qa, sbw, L >= 0 and sb0, sg > 0");
	}
	const int Kb = static_cast<int>(std::floor(Out0 / InsKdt + 0.5)) - 1, Ke = static_cast<int>(std::floor((Out0 + L) / InsKdt + 0.5)) - 1;
	if (Kb < 0 || Ke > InsNk)
	{
		return Result::Fail("the outage [out0, out0 + L] must lie inside the 300 s run");
	}
	const KfModel M = MakeKf(Qa, Sbw);
	M3 P = {100, 0, 0, 0, 1, 0, 0, 0, Sb0 * Sb0}, P1400 = P, PKb = P, PKe = P;
	const double R = Sg * Sg;
	for (int k = 1; k <= InsNk; ++k)
	{
		P = MAdd3(Fapt(M.F, P), M.Q);
		const double Tk = k * InsKdt + 1e-9;
		const bool Fix = k % InsFix == 0 && !(Out0 <= Tk && Tk < Out0 + L);
		if (Fix)
		{
			const double S = P[0] + R, K[3] = {P[0] / S, P[3] / S, P[6] / S};
			const M3 Ik = {1 - K[0], 0, 0, -K[1], 1, 0, -K[2], 0, 1};
			M3 Kr;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					Kr[static_cast<size_t>(3 * i + j)] = K[i] * K[j] * R;
				}
			}
			P = MAdd3(Fapt(Ik, P), Kr);
		}
		if (k == 1400)
		{
			P1400 = P;
		}
		if (k == Kb)
		{
			PKb = P;
		}
		if (k == Ke)
		{
			PKe = P;
		}
	}
	return Result::Record({{"sp140", std::sqrt(P1400[0])}, {"sbias", std::sqrt(PKb[8])}, {"sp0", std::sqrt(PKb[0])}, {"spend", std::sqrt(PKe[0])},
		{"closed", OutageClosed(PKb, (Ke - Kb) * InsKdt, M.Sa, Sbw)}});
}
Result KfSteady(const Call& C)
{
	const double Qa = C.N(0), Sbw = C.N(1), Sg = C.N(2);
	if (!(Qa >= 0) || !(Sbw > 0) || !(Sg > 0))
	{
		return Result::Fail("needs qa >= 0, sbw > 0 and sg > 0 (with sbw = 0 the bias becomes exactly known)");
	}
	const KfModel M = MakeKf(Qa, Sbw);
	M3 Phi = {1, 0, 0, 0, 1, 0, 0, 0, 1}, Q1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0; i < InsFix; ++i)
	{
		Phi = MMul3(M.F, Phi);
		Q1 = MAdd3(Fapt(M.F, Q1), M.Q);
	}
	// Prior-form filter Riccati P = Phi P Phi' - Phi P H'(H P H' + R)^-1 H P Phi' + Q1 as a control DARE with A = Phi', B = H'.
	Matrix A(3, 3), Qm(3, 3), Hb(3, 1, 0.0), Rm = Matrix::Scalar(Sg * Sg), P, K;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			A(i, j) = Phi[static_cast<size_t>(3 * j + i)];
			Qm(i, j) = Q1[static_cast<size_t>(3 * i + j)];
		}
	}
	Hb(0, 0) = 1;
	std::string Error;
	if (!la::Dare(A, Hb, la::Symmetrize(Qm), Rm, &P, &K, &Error))
	{
		return Result::Fail(Error);
	}
	const double S = P(0, 0) + Sg * Sg, Kp = P(0, 0) / S, Kv = P(1, 0) / S, Kbias = P(2, 0) / S;
	const double Pp = P(0, 0) - P(0, 0) * P(0, 0) / S, Pv = P(1, 1) - P(1, 0) * P(1, 0) / S, Pb = P(2, 2) - P(2, 0) * P(2, 0) / S;
	return Result::Record({{"sp", std::sqrt(Pp)}, {"sv", std::sqrt(Pv)}, {"sb", std::sqrt(Pb)}, {"sppri", std::sqrt(P(0, 0))}, {"kp", Kp}, {"kv", Kv}, {"kb", Kbias}});
}
Result InsOutage(const Call& C)
{
	const Matrix& Pm = C.M(0);
	const double L = C.N(1), Qa = C.N(2), Sbw = C.N(3);
	if (Pm.Count() != 9)
	{
		return Result::Fail("P must be the 3 x 3 covariance of [position, velocity, bias]");
	}
	if (!(L >= 0))
	{
		return Result::Fail("outage length must be >= 0");
	}
	M3 P;
	for (size_t i = 0; i < 9; ++i)
	{
		P[i] = Pm.V[i];
	}
	return Result::Of(OutageClosed(P, L, Qa / std::sqrt(InsKdt), Sbw));
}
} // namespace orbits_ins_detail

void RegisterOrbitsLibrary(Registry& R)
{
	{
		using namespace orbits_detail;
		const char* O = "orbits";
		const double Mu = cst::MuEarth;
		Def(R, "orbperiod", O, "orbits/orbits.js keplerState (s1 Kepler)", "orbital period 2 pi sqrt(a^3/mu)", "s", {Num("a", "m"), Num("mu", "m^3/s^2", Mu)}, &OrbPeriod);
		Def(R, "vcirc", O, "orbits/orbits.js hohmann (s2 Hohmann)", "circular orbit speed sqrt(mu/r)", "m/s", {Num("r", "m"), Num("mu", "m^3/s^2", Mu)}, &VCirc);
		Def(R, "vesc", O, "orbits/ two-body energy", "escape speed sqrt(2 mu/r)", "m/s", {Num("r", "m"), Num("mu", "m^3/s^2", Mu)}, &VEsc);
		Def(R, "visviva", O, "orbits/orbits.js hohmann (vis-viva)", "speed on a conic at radius r, sqrt(mu (2/r - 1/a))", "m/s", {Num("r", "m"), Num("a", "m"), Num("mu", "m^3/s^2", Mu)}, &VisViva);
		Def(R, "meanmotion", O, "orbits/orbits.js meanMotion (s4 CW)", "mean motion sqrt(mu/a^3)", "rad/s", {Num("a", "m"), Num("mu", "m^3/s^2", Mu)}, &MeanMotionFn);
		Def(R, "orbenergy", O, "orbits/orbits.js energy", "specific orbital energy v^2/2 - mu/r", "J/kg", {Mat("r", "m [x, y]"), Mat("v", "m/s [vx, vy]"), Num("mu", "m^3/s^2", Mu)}, &OrbEnergy);
		Def(R, "stumpc", O, "orbits/orbits.js stumpC (s3 Lambert)", "Stumpff function C(z)", "-", {Num("z", "-")}, &StumpCFn);
		Def(R, "stumps", O, "orbits/orbits.js stumpS (s3 Lambert)", "Stumpff function S(z)", "-", {Num("z", "-")}, &StumpSFn);
		Def(R, "kepler", O, "orbits/orbits.js keplerE (s1 Kepler)", "Kepler's equation M = E - e sin E by Newton", "record E (rad), nu (rad), iters", {Num("M", "rad"), Num("e", "-")}, &Kepler);
		Def(R, "trueanom", O, "orbits/ s1 Kepler", "true anomaly from eccentric anomaly (same revolution)", "rad", {Num("E", "rad"), Num("e", "-")}, &TrueAnom);
		Def(R, "kepstate", O, "orbits/orbits.js keplerState (s1 Kepler)", "position and velocity t seconds after periapsis", "record x, y (m), vx, vy (m/s), r, nu, E, M, period, iters",
			{Num("a", "m", 10000e3), Num("e", "-", 0.3), Num("t", "s", 6000), Num("mu", "m^3/s^2", Mu)}, &KepState);
		Def(R, "twobody", O, "orbits/orbits.js propagate", "RK4 integration of r'' = -mu r/|r|^3, n steps per local period", "record x, y, vx, vy, steps",
			{Mat("r", "m [x, y]"), Mat("v", "m/s [vx, vy]"), Num("t", "s"), Num("n", "-", 2000), Num("mu", "m^3/s^2", Mu)}, &TwoBody);
		Def(R, "hohmann", O, "orbits/orbits.js hohmann (s2 Hohmann)", "Hohmann transfer between circular orbits r1 -> r2", "record dv1, dv2, dv (m/s), tof (s), vp, va, v1, v2 (m/s), at (m), et",
			{Num("r1", "m", 6678e3), Num("r2", "m", 42164e3), Num("mu", "m^3/s^2", Mu)}, &Hohmann);
		Def(R, "hohprop", O, "orbits/ s2 Hohmann (rocket equation)", "Hohmann total delta-v and propellant by the rocket equation", "record dv (m/s), mr = m0/mf, frac (propellant share)",
			{Num("r1", "m", 6678e3), Num("r2", "m", 42164e3), Num("isp", "s", 320), Num("mu", "m^3/s^2", Mu)}, &HohProp);
		Def(R, "hohfly", O, "orbits/orbits.js hohmannFly (s2 Hohmann)", "fly the Hohmann ellipse by RK4 for the transfer time", "record rend (m), vr, speed (m/s), rext (m), steps",
			{Num("r1", "m", 6678e3), Num("r2", "m", 42164e3), Num("mu", "m^3/s^2", Mu)}, &HohFly);
		Def(R, "lambert", O, "orbits/orbits.js lambert (s3 Lambert)", "prograde zero-rev Lambert by universal variables", "record v1x, v1y, v2x, v2y, v1, v2 (m/s), z, dth (rad)",
			{Mat("r1", "m [x, y]"), Mat("r2", "m [x, y]"), Num("dt", "s"), Num("mu", "m^3/s^2", Mu)}, &Lambert);
		Def(R, "conic", O, "orbits/orbits.js elements (s3 Lambert)", "conic elements from a planar state", "record a (m, < 0 hyperbola), e, ex, ey, h (m^2/s), p, rp (m), energy (J/kg)",
			{Mat("r", "m [x, y]"), Mat("v", "m/s [vx, vy]"), Num("mu", "m^3/s^2", Mu)}, &Conic);
		Def(R, "cwstm", O, "orbits/orbits.js cwStm (s4 Clohessy-Wiltshire)", "CW state-transition matrix, [x y vx vy] or [x y z vx vy vz]", "4 x 4 or 6 x 6 matrix",
			{Num("n", "rad/s"), Num("t", "s"), Num("dim", "-", 4)}, &CwStmFn);
		Def(R, "cwdrift", O, "orbits/orbits.js cwStm, mv (s4 free drift)", "CW free drift of a relative state for time t", "record x, y, (z,) vx, vy, (vz,) range",
			{Num("n", "rad/s"), Num("t", "s"), Mat("X", "m, m/s")}, &CwDrift);
		Def(R, "cwnonlin", O, "orbits/orbits.js relNonlinear (s4)", "nonlinear relative motion about a circular orbit, RK4", "record x, y (m), vx, vy (m/s), range (m)",
			{Num("rref", "m"), Mat("X", "m, m/s [x, y, vx, vy]"), Num("t", "s"), Num("steps", "-", 4000), Num("mu", "m^3/s^2", Mu)}, &CwNonlin);
		Def(R, "rendezvous", O, "orbits/orbits.js rendezvous (s4 two-impulse)", "two-impulse CW rendezvous to the origin in time T", "record v0x, v0y, dv1x, dv1y, dv2x, dv2y, dv1, dv2, dv (m/s), n",
			{Num("rref", "m"), Mat("X", "m, m/s [x, y, vx, vy]"), Num("T", "s"), Num("mu", "m^3/s^2", Mu)}, &Rendezvous);
		Def(R, "cwshoot", O, "orbits/orbits.js cwShoot (s4)", "first rendezvous burn corrected by shooting in the nonlinear dynamics", "record v0x, v0y (m/s), iters, corr (m/s)",
			{Num("rref", "m"), Mat("X", "m, m/s [x, y, vx, vy]"), Num("T", "s"), Num("steps", "-", 4000), Num("mu", "m^3/s^2", Mu)}, &CwShoot);
		Def(R, "edrift", O, "orbits/orbits.js drift (s5 energy drift)", "max relative energy error, a = 10000 km; method 0 Euler, 1 leapfrog, 2 RK4", "record max, first, last (per-orbit max)",
			{Num("e", "-", 0.3), Num("spo", "steps/orbit", 500), Num("orbits", "-", 50), Num("method", "-", 1)}, &EDrift);
		Def(R, "edriftcmp", O, "orbits/ s5 energy drift", "leapfrog, RK4 and Euler energy drift; q = leapfrog(spo)/leapfrog(2 spo)", "record lf, rk4, euler, q, bound",
			{Num("e", "-", 0.3), Num("spo", "steps/orbit", 500), Num("orbits", "-", 50)}, &EDriftCmp);
		Def(R, "hohtot", O, "orbits/orbits.js hohTot (s6 bi-elliptic)", "Hohmann total delta-v in units of the inner circular speed", "-", {Num("R", "r2/r1")}, &HohTotFn);
		Def(R, "bieltot", O, "orbits/orbits.js bielTot (s6 bi-elliptic)", "bi-elliptic total delta-v in units of the inner circular speed", "-", {Num("R", "r2/r1"), Num("Rb", "rb/r1")}, &BielTotFn);
		Def(R, "bielinf", O, "orbits/orbits.js bielInf (s6 bi-elliptic)", "bi-elliptic total as rb -> infinity, (sqrt2 - 1)(1 + 1/sqrt R)", "-", {Num("R", "r2/r1")}, &BielInfFn);
		Def(R, "bielthresh", O, "orbits/orbits.js thresholds (s6 bi-elliptic)", "ratio R = r2/r1 above which bi-elliptic (rb -> inf) beats Hohmann, 11.93877", "record closed (rb = inf), finite (rb = 1e12)", {}, &BielThresh);
		Def(R, "bielliptic", O, "orbits/ s6 bi-elliptic", "bi-elliptic transfer r1 -> rb -> r2 between circular orbits", "record dv1, dv2, dv3, dv (m/s), tof (s)",
			{Num("r1", "m"), Num("r2", "m"), Num("rb", "m"), Num("mu", "m^3/s^2", Mu)}, &Bielliptic);
		Def(R, "breakeven", O, "orbits/orbits.js breakeven (s6 bi-elliptic)", "smallest rb/r1 at which bi-elliptic matches Hohmann (inf if never)", "rb/r1", {Num("R", "r2/r1")}, &Breakeven);
		Def(R, "hohpeak", O, "orbits/orbits.js hohPeak (s6 bi-elliptic)", "costliest Hohmann ratio (15.58; beyond it every bi-elliptic wins)", "record R, v", {}, &HohPeak);
		Def(R, "onetangent", O, "orbits/orbits.js oneTangent (s6)", "one-tangent-burn transfer cost, apoapsis ra", "-", {Num("R", "r2/r1"), Num("ra", "ra/r1")}, &OneTangentFn);
		Def(R, "onetanscan", O, "orbits/orbits.js oneTangentScan (s6)", "cheapest one-tangent transfer over ra in [R, 5R]", "record best, ratio (ra/r2)", {Num("R", "r2/r1"), Num("n", "-", 4001)}, &OneTanScan);
	}
	{
		using namespace orbits_ins_detail;
		const char* I = "ins";
		const double B0 = 0.1e-3 * cst::G0, Eps0 = 0.01 * cst::Pi / 180 / 3600;
		const double Qa0 = 300e-6 * cst::G0, Sb00 = 5e-3 * cst::G0, Sbw0 = 20e-6 * cst::G0;
		Def(R, "allandev", I, "ins/ins.js modelAdev (s1 Allan)", "model Allan deviation sqrt(N^2/tau + (0.664 B)^2 + K^2 tau/3)", "rate unit",
			{Num("tau", "s"), Num("N", "u sqrt(s)"), Num("B", "u", 0), Num("K", "u/sqrt(s)", 0)}, &AllanDev);
		Def(R, "allanmin", I, "ins/ s1 Allan", "bottom of the model Allan curve, tau = sqrt3 N/K", "record tau (s), adev", {Num("N", "u sqrt(s)"), Num("B", "u"), Num("K", "u/sqrt(s)")}, &AllanMin);
		Def(R, "flickerquad", I, "ins/ins.js flickerQuad (s1 Allan)", "int_0^inf sin^4 x / x^3 dx (= ln 2) by Simpson plus tail", "-", {Num("X", "-", 400), Num("n", "-", 400000)}, &FlickerQuad);
		Def(R, "allanfit", I, "ins/ins.js allanSet, fits (s1 Allan)", "seeded unit noise records: fitted slopes, N, K and the flicker floor", "record wslope, nhat, kslope, khat, flat",
			{Num("seed", "-", 1)}, &AllanFit);
		Def(R, "allansim", I, "ins/ins.js combinedAdev (s1 Allan)", "Allan deviation of a seeded simulated record N w + B f + K k", "m x 2 matrix [tau (s), adev]",
			{Num("N", "u sqrt(s)"), Num("B", "u"), Num("K", "u/sqrt(s)"), Num("seed", "-", 1)}, &AllanSim);
		Def(R, "schuler", I, "ins/ins.js PERIOD (s2 Schuler)", "Schuler period 2 pi sqrt(R/g) (84.41 min)", "s", {Num("R", "m", InsRe), Num("g", "m/s^2", InsG)}, &SchulerPeriod);
		Def(R, "schulerx", I, "ins/ins.js schulerClosed (s2 Schuler)", "closed-form Schuler-loop position error", "m",
			{Num("t", "s"), Num("b", "m/s^2", B0), Num("phi0", "rad", 0), Num("eps", "rad/s", Eps0)}, &SchulerX);
		Def(R, "naivedrift", I, "ins/ins.js naiveDrift (s2 Schuler)", "flat-earth position drift without Schuler feedback", "m",
			{Num("t", "s"), Num("b", "m/s^2", B0), Num("phi0", "rad", 0), Num("eps", "rad/s", Eps0)}, &NaiveDrift);
		Def(R, "schulersim", I, "ins/ins.js schulerRK4, measuredPeriod (s2)", "RK4 Schuler error loop: end error, gap to the closed form, measured period", "record x (m), dmax (m), period (s)",
			{Num("b", "m/s^2", B0), Num("phi0", "rad", 0), Num("eps", "rad/s", Eps0), Num("tend", "s", 10800), Num("h", "s", 1)}, &SchulerSim);
		Def(R, "coning0", I, "ins/ins.js drift0 (s3 coning)", "coning drift with no correction", "rad/s", {Num("alpha", "rad"), Num("W", "rad/s"), Num("h", "s")}, &Coning0);
		Def(R, "coning2", I, "ins/ins.js drift2 (s3 coning)", "residual coning drift after the two-sample correction", "rad/s", {Num("alpha", "rad"), Num("W", "rad/s"), Num("h", "s")}, &Coning2);
		Def(R, "coningrun", I, "ins/ins.js coningRun (s3 coning)", "integrate coning motion uncorrected and with the 2/3 two-sample term", "record slope0, slope2 (rad/s), e0, e2 (rad)",
			{Num("alpha", "rad"), Num("W", "rad/s"), Num("rate", "Hz", 200), Num("tend", "s", 2)}, &ConingRun);
		Def(R, "bortz", I, "ins/ins.js bortzCoefficient (s3 coning)", "coefficient of d1 x d2 from Bortz's equation (-> 2/3)", "-", {Num("h", "s", 1e-2), Num("steps", "-", 2000)}, &Bortz);
		Def(R, "insgnss", I, "ins/ins.js kfModel, kfCov, outageClosed (s4)", "1-D INS/GNSS Kalman covariance, 1 Hz fixes, outage from out0 for L", "record sp140, sbias, sp0, spend, closed",
			{Num("qa", "m/s^2/sqrt(Hz)", Qa0), Num("sb0", "m/s^2", Sb00), Num("sbw", "m/s^2/sqrt(s)", Sbw0), Num("sg", "m", 3), Num("L", "s", 60), Num("out0", "s", 150)}, &InsGnss);
		Def(R, "kfsteady", I, "ins/ins.js riccati1s (s4 Kalman)", "steady-state 1-D INS/GNSS filter by the DARE", "record sp, sv, sb (posterior sigmas), sppri (m), kp, kv, kb",
			{Num("qa", "m/s^2/sqrt(Hz)", Qa0), Num("sbw", "m/s^2/sqrt(s)", Sbw0), Num("sg", "m", 3)}, &KfSteady);
		Def(R, "insoutage", I, "ins/ins.js outageClosed (s4 Kalman)", "position sigma after an outage of L seconds from covariance P", "m",
			{Mat("P", "3 x 3 [m, m/s, m/s^2]"), Num("L", "s"), Num("qa", "m/s^2/sqrt(Hz)", Qa0), Num("sbw", "m/s^2/sqrt(s)", Sbw0)}, &InsOutage);
	}
}
} // namespace mcalc
