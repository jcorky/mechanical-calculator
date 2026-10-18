// Orbits and inertial navigation: the orbits/ bench (Kepler, Hohmann, Lambert, conics, Clohessy-Wiltshire,
// rendezvous, integrator energy drift, bi-elliptic transfers) and the ins/ bench (Allan deviation, Schuler loop,
// coning, 1-D INS/GNSS Kalman filter). Ported from orbits/orbits.js and ins/ins.js; plane vectors are [x, y].
// All arithmetic runs on the stepped-drum machine (Real); integers remain only as counters and RNG state.
#include "CalcLibrary.h"
#include "GearMath.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace mcalc
{
namespace orbits_detail
{
struct OVec
{
	Real X;
	Real Y;
};

OVec OV(const Real& X, const Real& Y)
{
	OVec R;
	R.X = X;
	R.Y = Y;
	return R;
}
Real ONorm(const OVec& A) { return gm::hypot(A.X, A.Y); }
Real ODot(const OVec& A, const OVec& B) { return A.X * B.X + A.Y * B.Y; }
// Gravity at R whose norm D is already known: -mu r / D^3.
OVec OAccD(const OVec& R, const Real& D, const Real& Mu)
{
	const Real K = -Mu / (D * D * D);
	return OV(K * R.X, K * R.Y);
}
Real OEnergyD(const OVec& V, const Real& D, const Real& Mu) { return Real(0.5) * ODot(V, V) - Mu / D; }
Real OEnergy(const OVec& R, const OVec& V, const Real& Mu) { return OEnergyD(V, ONorm(R), Mu); }

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

// Classical fourth-order Runge-Kutta step of r'' = -mu r / |r|^3 (orbits.js rk4); D is |R| on entry.
void ORk4(OVec& R, OVec& V, const Real& Dt, const Real& Mu, const Real& D)
{
	const Real H2 = Dt / 2, H6 = Dt / 6;
	const OVec A1 = OAccD(R, D, Mu);
	const OVec R2 = OV(R.X + H2 * V.X, R.Y + H2 * V.Y), V2 = OV(V.X + H2 * A1.X, V.Y + H2 * A1.Y), A2 = OAccD(R2, ONorm(R2), Mu);
	const OVec R3 = OV(R.X + H2 * V2.X, R.Y + H2 * V2.Y), V3 = OV(V.X + H2 * A2.X, V.Y + H2 * A2.Y), A3 = OAccD(R3, ONorm(R3), Mu);
	const OVec R4 = OV(R.X + Dt * V3.X, R.Y + Dt * V3.Y), V4 = OV(V.X + Dt * A3.X, V.Y + Dt * A3.Y), A4 = OAccD(R4, ONorm(R4), Mu);
	const OVec RN = OV(R.X + H6 * (V.X + 2 * V2.X + 2 * V3.X + V4.X), R.Y + H6 * (V.Y + 2 * V2.Y + 2 * V3.Y + V4.Y));
	const OVec VN = OV(V.X + H6 * (A1.X + 2 * A2.X + 2 * A3.X + A4.X), V.Y + H6 * (A1.Y + 2 * A2.Y + 2 * A3.Y + A4.Y));
	R = RN;
	V = VN;
}

struct OPropOut
{
	OVec R;
	OVec V;
	int Steps = 0;
	Real RMax;   // extremes of |r| over the samples orbits.js keeps (every `Every` steps plus both ends)
	Real RMin;
	bool Ok = true;
};

constexpr int OrbMaxSteps = 100000;

// Integrate for time T with step (local circular period)/N, the last step clipped to land on T.
OPropOut OPropagate(OVec R, OVec V, const Real& T, const Real& N, const Real& Mu, int Every)
{
	OPropOut Out;
	const Real TauN = cst::Tau() / N;
	Real Done = 0;
	Real D = ONorm(R);
	int K = 0;
	Out.RMax = Out.RMin = D;
	while (Done < T)
	{
		const Real Dt = gm::min(T - Done, TauN * gm::sqrt(D * D * D / Mu));
		ORk4(R, V, Dt, Mu, D);
		D = ONorm(R);
		Done += Dt;
		++K;
		if (Every > 0 && K % Every == 0)
		{
			Out.RMax = gm::max(Out.RMax, D);
			Out.RMin = gm::min(Out.RMin, D);
		}
		if (K > OrbMaxSteps || !gm::isfinite(Dt) || Dt <= 0)
		{
			Out.Ok = false;
			break;
		}
	}
	Out.RMax = gm::max(Out.RMax, D);
	Out.RMin = gm::min(Out.RMin, D);
	Out.R = R;
	Out.V = V;
	Out.Steps = K;
	return Out;
}

Real KeplerSolve(const Real& M, const Real& E0, int* Iters)
{
	Real E = E0 < Real(0.8) ? M : cst::Pi();
	const Real Tol = 1e-14;
	for (int i = 0; i < 50; ++i)
	{
		const Real D = (E - E0 * gm::sin(E) - M) / (1 - E0 * gm::cos(E));
		E -= D;
		if (gm::fabs(D) < Tol)
		{
			*Iters = i + 1;
			return E;
		}
	}
	*Iters = 50;
	return E;
}

// True anomaly from eccentric anomaly, continuous in E (same revolution).
Real TrueFromEcc(const Real& E, const Real& Ecc)
{
	const Real Beta = Ecc / (1 + gm::sqrt(1 - Ecc * Ecc));
	return E + 2 * gm::atan(Beta * gm::sin(E) / (1 - Beta * gm::cos(E)));
}

// Stumpff C(z) and S(z) together (they share sqrt z); series near zero as in orbits.js.
void StumpCS(const Real& Z, Real* Cz, Real* Sz)
{
	const Real Lim = 1e-8;
	if (Z > Lim)
	{
		const Real S = gm::sqrt(Z);
		*Cz = (1 - gm::cos(S)) / Z;
		*Sz = (S - gm::sin(S)) / (S * S * S);
	}
	else if (Z < -Lim)
	{
		const Real S = gm::sqrt(-Z);
		*Cz = (gm::cosh(S) - 1) / (-Z);
		*Sz = (gm::sinh(S) - S) / (S * S * S);
	}
	else
	{
		*Cz = Real(0.5) - Z / 24;
		*Sz = Real(1) / 6 - Z / 120;
	}
}

struct LambertOut
{
	OVec V1;
	OVec V2;
	Real Z;
	Real Dth;
};

struct LambertFn
{
	Real R1, R2, A, SqMuDt;
	Real Yz(const Real& Z) const
	{
		Real Cz, Sz;
		StumpCS(Z, &Cz, &Sz);
		return R1 + R2 + A * (Z * Sz - 1) / gm::sqrt(Cz);
	}
	Real Fz(const Real& Z) const
	{
		Real Cz, Sz;
		StumpCS(Z, &Cz, &Sz);
		const Real Y = R1 + R2 + A * (Z * Sz - 1) / gm::sqrt(Cz), Q = Y / Cz;
		return Q * gm::sqrt(Q) * Sz + A * gm::sqrt(Y) - SqMuDt;   // (y/C)^1.5 S + A sqrt(y) - sqrt(mu) dt
	}
};

// Prograde zero-revolution Lambert by universal variables, z by bisection (orbits.js lambert).
bool LambertSolve(const OVec& R1v, const OVec& R2v, const Real& Dt, const Real& Mu, LambertOut* Out, std::string* Error)
{
	LambertFn L;
	L.R1 = ONorm(R1v);
	L.R2 = ONorm(R2v);
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
	L.SqMuDt = gm::sqrt(Mu) * Dt;
	const Real Cross = R1v.X * R2v.Y - R1v.Y * R2v.X;
	Real Dth = gm::acos(gm::max(Real(-1), gm::min(Real(1), ODot(R1v, R2v) / (L.R1 * L.R2))));
	if (Cross < 0)
	{
		Dth = cst::Tau() - Dth;
	}
	const Real SinD = gm::sin(Dth);
	if (gm::fabs(SinD) < Real(1e-9))
	{
		*Error = "transfer angle of 0 or 180 deg is degenerate";
		return false;
	}
	L.A = SinD * gm::sqrt(L.R1 * L.R2 / (1 - gm::cos(Dth)));
	Real Hi = 4 * cst::Pi() * cst::Pi() - Real(1e-9), Lo = -4;
	for (int Guard = 0; Guard < 200 && L.Yz(Lo) > 0 && L.Fz(Lo) > 0; ++Guard)
	{
		Lo *= 2;
	}
	if (L.Yz(Lo) <= 0)
	{
		Real A = Lo, B = Hi;
		for (int i = 0; i < 200; ++i)
		{
			const Real M = (A + B) / 2;
			if (M == A || M == B)
			{
				break;   // bracket is two adjacent 18-digit values: further halving changes nothing
			}
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
		const Real M = (Lo + Hi) / 2;
		if (M == Lo || M == Hi)
		{
			break;
		}
		if (L.Fz(M) < 0)
		{
			Lo = M;
		}
		else
		{
			Hi = M;
		}
	}
	const Real Z = (Lo + Hi) / 2, Y = L.Yz(Z), F = 1 - Y / L.R1, G = L.A * gm::sqrt(Y / Mu), Gd = 1 - Y / L.R2;
	Out->V1 = OV((R2v.X - F * R1v.X) / G, (R2v.Y - F * R1v.Y) / G);
	Out->V2 = OV((Gd * R2v.X - R1v.X) / G, (Gd * R2v.Y - R1v.Y) / G);
	Out->Z = Z;
	Out->Dth = Dth;
	if (!gm::isfinite(Out->V1.X) || !gm::isfinite(Out->V1.Y) || !gm::isfinite(Out->V2.X) || !gm::isfinite(Out->V2.Y))
	{
		*Error = "no prograde zero-revolution solution for this time of flight";
		return false;
	}
	return true;
}

// Clohessy-Wiltshire state-transition matrix, state [x, y, vx, vy] (x radial out, y along-track).
Matrix CwStm4(const Real& N, const Real& T)
{
	const Real Nt = N * T, S = gm::sin(Nt), C = gm::cos(Nt);
	return Matrix::FromRows({{4 - 3 * C, 0, S / N, 2 * (1 - C) / N},
		{6 * (S - Nt), 1, -2 * (1 - C) / N, (4 * S - 3 * Nt) / N},
		{3 * N * S, 0, C, 2 * S},
		{-6 * N * (1 - C), 0, -2 * S, 4 * C - 3}});
}
// Six-state version [x, y, z, vx, vy, vz]: the in-plane block plus the decoupled cross-track oscillator.
Matrix CwStm6(const Real& N, const Real& T)
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
	const Real S = gm::sin(N * T), C = gm::cos(N * T);
	P(2, 2) = C;
	P(2, 5) = S / N;
	P(5, 2) = -N * S;
	P(5, 5) = C;
	return P;
}

bool Solve2(const Real& M00, const Real& M01, const Real& M10, const Real& M11, const Real& B0, const Real& B1, OVec* Out)
{
	const Real D = M00 * M11 - M01 * M10, Scale = gm::max(gm::max(gm::fabs(M00), gm::fabs(M01)), gm::max(gm::fabs(M10), gm::fabs(M11)));
	if (!(gm::fabs(D) > Real(1e-12) * Scale * Scale) || !gm::isfinite(D))
	{
		return false;
	}
	*Out = OV((B0 * M11 - M01 * B1) / D, (M00 * B1 - M10 * B0) / D);
	return true;
}

// First burn velocity that brings [x, y] to the origin in time T under CW (orbits.js cwBurn).
bool CwBurn(const Real& N, const Real* X, const Real& T, OVec* V0)
{
	const Matrix P = CwStm4(N, T);
	return Solve2(P(0, 2), P(0, 3), P(1, 2), P(1, 3), -(P(0, 0) * X[0] + P(0, 1) * X[1]), -(P(1, 0) * X[0] + P(1, 1) * X[1]), V0);
}

Real OMeanMotion(const Real& Rref, const Real& Mu) { return gm::sqrt(Mu / (Rref * Rref * Rref)); }

using OState = std::array<Real, 4>;

// Full nonlinear relative motion about a circular reference orbit, rotating frame, fixed-step RK4.
OState RelNonlinear(const Real& Rref, const OState& X, const Real& T, int Steps, const Real& Mu)
{
	const Real N = OMeanMotion(Rref, Mu), G0 = Mu / (Rref * Rref), N2 = N * N, TwoN = 2 * N;
	auto Fn = [&](const OState& S) {
		const Real Rx = Rref + S[0], Q = Rx * Rx + S[1] * S[1], D = Q * gm::sqrt(Q);   // |r|^3
		return OState{S[2], S[3], TwoN * S[3] + N2 * S[0] - Mu * Rx / D + G0, -TwoN * S[2] + N2 * S[1] - Mu * S[1] / D};
	};
	OState S = X, Tmp;
	const Real H = T / Steps, H2 = H / 2, H6 = H / 6;
	for (int i = 0; i < Steps; ++i)
	{
		const OState K1 = Fn(S);
		for (size_t j = 0; j < 4; ++j)
		{
			Tmp[j] = S[j] + H2 * K1[j];
		}
		const OState K2 = Fn(Tmp);
		for (size_t j = 0; j < 4; ++j)
		{
			Tmp[j] = S[j] + H2 * K2[j];
		}
		const OState K3 = Fn(Tmp);
		for (size_t j = 0; j < 4; ++j)
		{
			Tmp[j] = S[j] + H * K3[j];
		}
		const OState K4 = Fn(Tmp);
		for (size_t j = 0; j < 4; ++j)
		{
			S[j] = S[j] + H6 * (K1[j] + 2 * K2[j] + 2 * K3[j] + K4[j]);
		}
	}
	return S;
}

bool OGetState4(const Matrix& M, OState* Out)
{
	if (M.Count() != 4)
	{
		return false;
	}
	for (size_t i = 0; i < 4; ++i)
	{
		(*Out)[i] = M.V[i];
	}
	return true;
}

// Largest relative energy error per orbit of a fixed-step integrator on a Kepler orbit with a = 10 000 km.
// Method 0 = explicit Euler, 1 = leapfrog (kick-drift-kick), 2 = RK4. |r| and the acceleration at the
// current point are carried from one step to the next instead of being recomputed (same values).
std::vector<Real> EnergyDrift(const Real& Ecc, int Spo, int Orbits, int Method)
{
	const Real Mu = cst::MuEarth(), A = 10000e3, Rp = A * (1 - Ecc);
	OVec R = OV(Rp, 0), V = OV(0, gm::sqrt(Mu * (1 + Ecc) / Rp));
	Real D = ONorm(R);
	const Real E0 = OEnergyD(V, D, Mu), H = cst::Tau() * gm::sqrt(A * A * A / Mu) / Spo, H2 = H / 2;
	OVec G = OAccD(R, D, Mu);
	std::vector<Real> Out(static_cast<size_t>(Orbits));
	for (int k = 0; k < Orbits; ++k)
	{
		Real M = 0;
		for (int i = 0; i < Spo; ++i)
		{
			if (Method == 1)
			{
				V = OV(V.X + H2 * G.X, V.Y + H2 * G.Y);
				R = OV(R.X + H * V.X, R.Y + H * V.Y);
				D = ONorm(R);
				G = OAccD(R, D, Mu);
				V = OV(V.X + H2 * G.X, V.Y + H2 * G.Y);
			}
			else if (Method == 2)
			{
				ORk4(R, V, H, Mu, D);
				D = ONorm(R);
			}
			else
			{
				R = OV(R.X + H * V.X, R.Y + H * V.Y);
				V = OV(V.X + H * G.X, V.Y + H * G.Y);
				D = ONorm(R);
				G = OAccD(R, D, Mu);
			}
			M = gm::max(M, gm::fabs(OEnergyD(V, D, Mu) / E0 - 1));
		}
		Out[static_cast<size_t>(k)] = M;
	}
	return Out;
}

Real MaxOf(const std::vector<Real>& A)
{
	Real M = A.front();
	for (const Real& X : A)
	{
		M = gm::max(M, X);
	}
	return M;
}

// Hohmann and bi-elliptic totals in units of the inner circular speed (r1 = 1).
Real HohTot(const Real& R) { return gm::sqrt(2 * R / (1 + R)) - 1 + 1 / gm::sqrt(R) - gm::sqrt(2 / (R * (1 + R))); }
Real BielTot(const Real& R, const Real& Rb)
{
	const Real A1 = (1 + Rb) / 2, A2 = (R + Rb) / 2;
	return (gm::sqrt(2 - 1 / A1) - 1) + (gm::sqrt(2 / Rb - 1 / A2) - gm::sqrt(2 / Rb - 1 / A1)) + (gm::sqrt(2 / R - 1 / A2) - gm::sqrt(1 / R));
}
Real BielInf(const Real& R) { return (gm::sqrt(Real(2)) - 1) * (1 + 1 / gm::sqrt(R)); }
Real BielFar(const Real& R) { return BielTot(R, Real(1e12)); }
Real Crossing(Real (*Other)(const Real&))
{
	Real Lo = 5, Hi = 15;
	for (int i = 0; i < 100; ++i)
	{
		const Real M = (Lo + Hi) / 2;
		if (M == Lo || M == Hi)
		{
			break;
		}
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
Real OneTangent(const Real& R, const Real& Ra)
{
	const Real A = (1 + Ra) / 2, Vp = gm::sqrt(2 - 1 / A), Vt = Vp / R, Vr = gm::sqrt(gm::max(Real(0), 2 / R - 1 / A - Vt * Vt));
	return (Vp - 1) + gm::hypot(Vt - gm::sqrt(1 / R), Vr);
}

// ---------- Results ----------

Result OrbPeriod(const Call& C)
{
	const Real A = C.N(0), Mu = C.N(1);
	if (!(A > 0) || !(Mu > 0))
	{
		return Result::Fail("needs a > 0 and mu > 0");
	}
	return Result::Of(cst::Tau() * gm::sqrt(A * A * A / Mu));
}
Result VCirc(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("needs r > 0 and mu > 0");
	}
	return Result::Of(gm::sqrt(C.N(1) / C.N(0)));
}
Result VEsc(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("needs r > 0 and mu > 0");
	}
	return Result::Of(gm::sqrt(2 * C.N(1) / C.N(0)));
}
Result VisViva(const Call& C)
{
	const Real R = C.N(0), A = C.N(1), Mu = C.N(2);
	if (!(R > 0) || A == 0 || !(Mu > 0))
	{
		return Result::Fail("needs r > 0, a != 0 (negative for a hyperbola) and mu > 0");
	}
	const Real V2 = Mu * (2 / R - 1 / A);
	if (V2 < 0)
	{
		return Result::Fail("r lies beyond apoapsis 2a of this ellipse");
	}
	return Result::Of(gm::sqrt(V2));
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
Result StumpCFn(const Call& C)
{
	Real Cz, Sz;
	StumpCS(C.N(0), &Cz, &Sz);
	return Result::Of(Cz);
}
Result StumpSFn(const Call& C)
{
	Real Cz, Sz;
	StumpCS(C.N(0), &Cz, &Sz);
	return Result::Of(Sz);
}

Result Kepler(const Call& C)
{
	const Real M = C.N(0), Ecc = C.N(1);
	if (!(Ecc >= 0 && Ecc < 1))
	{
		return Result::Fail("eccentricity must be in [0, 1)");
	}
	int Iters = 0;
	const Real E = KeplerSolve(M, Ecc, &Iters);
	return Result::Record({{"E", E}, {"nu", TrueFromEcc(E, Ecc)}, {"iters", Real(Iters)}});
}
Result TrueAnom(const Call& C)
{
	const Real Ecc = C.N(1);
	if (!(Ecc >= 0 && Ecc < 1))
	{
		return Result::Fail("eccentricity must be in [0, 1)");
	}
	return Result::Of(TrueFromEcc(C.N(0), Ecc));
}
Result KepState(const Call& C)
{
	const Real A = C.N(0), Ecc = C.N(1), T = C.N(2), Mu = C.N(3);
	if (!(A > 0) || !(Mu > 0))
	{
		return Result::Fail("needs a > 0 and mu > 0");
	}
	if (!(Ecc >= 0 && Ecc < 1))
	{
		return Result::Fail("eccentricity must be in [0, 1)");
	}
	const Real Tau = cst::Tau(), N = gm::sqrt(Mu / (A * A * A)), M = gm::fmod(gm::fmod(N * T, Tau) + Tau, Tau);
	int Iters = 0;
	const Real E = KeplerSolve(M, Ecc, &Iters), SinE = gm::sin(E), CosE = gm::cos(E), B = A * gm::sqrt(1 - Ecc * Ecc), Ed = N / (1 - Ecc * CosE);
	const Real X = A * (CosE - Ecc), Y = B * SinE, Vx = -A * SinE * Ed, Vy = B * CosE * Ed;
	return Result::Record({{"x", X}, {"y", Y}, {"vx", Vx}, {"vy", Vy}, {"r", gm::hypot(X, Y)}, {"nu", TrueFromEcc(E, Ecc)}, {"E", E}, {"M", M}, {"period", Tau / N},
		{"iters", Real(Iters)}});
}
Result TwoBody(const Call& C)
{
	OVec R, V;
	if (!OGetVec2(C.M(0), &R) || !OGetVec2(C.M(1), &V))
	{
		return Result::Fail("r and v must be [x, y] vectors");
	}
	const Real T = C.N(2), N = C.N(3), Mu = C.N(4);
	if (!(ONorm(R) > 0) || !(T >= 0) || !(N >= 1) || !(Mu > 0))
	{
		return Result::Fail("needs r != 0, t >= 0, n >= 1 and mu > 0");
	}
	const OPropOut P = OPropagate(R, V, T, N, Mu, 0);
	if (!P.Ok)
	{
		return Result::Fail("more than 100000 steps; shorten t or lower n");
	}
	return Result::Record({{"x", P.R.X}, {"y", P.R.Y}, {"vx", P.V.X}, {"vy", P.V.Y}, {"steps", Real(P.Steps)}});
}

Result Hohmann(const Call& C)
{
	const Real R1 = C.N(0), R2 = C.N(1), Mu = C.N(2);
	if (!(R1 > 0) || !(R2 > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r1, r2 > 0 and mu > 0");
	}
	const Real At = (R1 + R2) / 2, V1 = gm::sqrt(Mu / R1), V2 = gm::sqrt(Mu / R2), Vp = gm::sqrt(Mu * (2 / R1 - 1 / At)), Va = gm::sqrt(Mu * (2 / R2 - 1 / At));
	const Real Dv1 = gm::fabs(Vp - V1), Dv2 = gm::fabs(V2 - Va);
	return Result::Record({{"dv1", Dv1}, {"dv2", Dv2}, {"dv", Dv1 + Dv2}, {"tof", cst::Pi() * gm::sqrt(At * At * At / Mu)}, {"vp", Vp}, {"va", Va}, {"v1", V1}, {"v2", V2}, {"at", At},
		{"et", gm::fabs(R2 - R1) / (R1 + R2)}});
}
Result HohProp(const Call& C)
{
	const Real R1 = C.N(0), R2 = C.N(1), Isp = C.N(2), Mu = C.N(3);
	if (!(R1 > 0) || !(R2 > 0) || !(Isp > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r1, r2, isp and mu > 0");
	}
	const Real At = (R1 + R2) / 2, Vp = gm::sqrt(Mu * (2 / R1 - 1 / At)), Va = gm::sqrt(Mu * (2 / R2 - 1 / At));
	const Real Dv = gm::fabs(Vp - gm::sqrt(Mu / R1)) + gm::fabs(gm::sqrt(Mu / R2) - Va), Mr = gm::exp(Dv / (cst::G0() * Isp));
	return Result::Record({{"dv", Dv}, {"mr", Mr}, {"frac", 1 - 1 / Mr}});
}
Result HohFly(const Call& C)
{
	const Real R1 = C.N(0), R2 = C.N(1), Mu = C.N(2);
	if (!(R1 > 0) || !(R2 > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r1, r2 > 0 and mu > 0");
	}
	const Real At = (R1 + R2) / 2, Vp = gm::sqrt(Mu * (2 / R1 - 1 / At)), Tof = cst::Pi() * gm::sqrt(At * At * At / Mu);
	const OPropOut P = OPropagate(OV(R1, 0), OV(0, Vp), Tof, 2000, Mu, 4);
	if (!P.Ok)
	{
		return Result::Fail("more than 100000 steps");
	}
	const Real REnd = ONorm(P.R);
	return Result::Record({{"rend", REnd}, {"vr", ODot(P.R, P.V) / REnd}, {"speed", ONorm(P.V)}, {"rext", R2 >= R1 ? P.RMax : P.RMin}, {"steps", Real(P.Steps)}});
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
	const Real Mu = C.N(2), Rn = ONorm(R);
	if (!(Rn > 0) || !(Mu > 0))
	{
		return Result::Fail("needs r != 0 and mu > 0");
	}
	const Real En = OEnergyD(V, Rn, Mu), H = R.X * V.Y - R.Y * V.X, Rv = ODot(R, V), Vv = ODot(V, V), W = Vv - Mu / Rn;
	const Real Ex = (W * R.X - Rv * V.X) / Mu, Ey = (W * R.Y - Rv * V.Y) / Mu, E = gm::hypot(Ex, Ey), P = H * H / Mu;
	return Result::Record({{"a", -Mu / (2 * En)}, {"e", E}, {"ex", Ex}, {"ey", Ey}, {"h", H}, {"p", P}, {"rp", P / (1 + E)}, {"energy", En}});
}

Result CwStmFn(const Call& C)
{
	const Real N = C.N(0), T = C.N(1);
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
	const Real N = C.N(0), T = C.N(1);
	const Matrix& X = C.M(2);
	if (!(N > 0))
	{
		return Result::Fail("mean motion n must be positive");
	}
	if (X.Count() == 4)
	{
		const Matrix S = la::Mul(CwStm4(N, T), Matrix::Column(X.V));
		return Result::Record({{"x", S.V[0]}, {"y", S.V[1]}, {"vx", S.V[2]}, {"vy", S.V[3]}, {"range", gm::hypot(S.V[0], S.V[1])}});
	}
	if (X.Count() == 6)
	{
		const Matrix S = la::Mul(CwStm6(N, T), Matrix::Column(X.V));
		return Result::Record({{"x", S.V[0]}, {"y", S.V[1]}, {"z", S.V[2]}, {"vx", S.V[3]}, {"vy", S.V[4]}, {"vz", S.V[5]}, {"range", gm::hypot(S.V[0], S.V[1], S.V[2])}});
	}
	return Result::Fail("state must be [x, y, vx, vy] or [x, y, z, vx, vy, vz]");
}
Result CwNonlin(const Call& C)
{
	OState X;
	const Real Rref = C.N(0), T = C.N(2), Mu = C.N(4);
	const int Steps = C.Int(3);
	if (!OGetState4(C.M(1), &X))
	{
		return Result::Fail("state must be [x, y, vx, vy]");
	}
	if (!(Rref > 0) || !(Mu > 0) || Steps < 1 || Steps > 100000)
	{
		return Result::Fail("needs rref > 0, mu > 0 and 1 <= steps <= 100000");
	}
	const OState S = RelNonlinear(Rref, X, T, Steps, Mu);
	return Result::Record({{"x", S[0]}, {"y", S[1]}, {"vx", S[2]}, {"vy", S[3]}, {"range", gm::hypot(S[0], S[1])}});
}
Result Rendezvous(const Call& C)
{
	OState X;
	const Real Rref = C.N(0), T = C.N(2), Mu = C.N(3);
	if (!OGetState4(C.M(1), &X))
	{
		return Result::Fail("chaser state must be [x, y, vx, vy]");
	}
	if (!(Rref > 0) || !(Mu > 0) || !(T > 0))
	{
		return Result::Fail("needs rref > 0, T > 0 and mu > 0");
	}
	const Real N = OMeanMotion(Rref, Mu);
	OVec V0;
	if (!CwBurn(N, X.data(), T, &V0))
	{
		return Result::Fail("T is a whole number of orbits: the transfer is singular");
	}
	const Matrix P = CwStm4(N, T), Xa = la::Mul(P, Matrix::Column({X[0], X[1], V0.X, V0.Y}));
	const OVec Dv1 = OV(V0.X - X[2], V0.Y - X[3]), Dv2 = OV(-Xa.V[2], -Xa.V[3]);
	const Real M1 = ONorm(Dv1), M2 = ONorm(Dv2);
	return Result::Record({{"v0x", V0.X}, {"v0y", V0.Y}, {"dv1x", Dv1.X}, {"dv1y", Dv1.Y}, {"dv2x", Dv2.X}, {"dv2y", Dv2.Y}, {"dv1", M1}, {"dv2", M2}, {"dv", M1 + M2}, {"n", N}});
}
Result CwShoot(const Call& C)
{
	OState X;
	const Real Rref = C.N(0), T = C.N(2), Mu = C.N(4);
	const int Steps = C.Int(3);
	if (!OGetState4(C.M(1), &X))
	{
		return Result::Fail("chaser state must be [x, y, vx, vy]");
	}
	if (!(Rref > 0) || !(Mu > 0) || !(T > 0) || Steps < 1 || Steps > 5000)
	{
		return Result::Fail("needs rref > 0, T > 0, mu > 0 and 1 <= steps <= 5000");
	}
	const Real N = OMeanMotion(Rref, Mu), Tol = 1e-2;
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
		const OState M = RelNonlinear(Rref, {X[0], X[1], V.X, V.Y}, T, Steps, Mu);
		if (gm::hypot(M[0], M[1]) < Tol)
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
	return Result::Record({{"v0x", V.X}, {"v0y", V.Y}, {"iters", Real(Iters)}, {"corr", gm::hypot(V.X - VCw.X, V.Y - VCw.Y)}});
}

bool DriftArgs(const Call& C, long long MaxSteps, Real* Ecc, int* Spo, int* Orbits, std::string* Error)
{
	*Ecc = C.N(0);
	*Spo = C.Int(1);
	*Orbits = C.Int(2);
	if (!(*Ecc >= 0 && *Ecc < 1))
	{
		*Error = "eccentricity must be in [0, 1)";
		return false;
	}
	if (*Spo < 10 || *Orbits < 1 || static_cast<long long>(*Spo) * *Orbits > MaxSteps)
	{
		*Error = "needs spo >= 10, orbits >= 1 and spo * orbits <= " + std::to_string(MaxSteps);
		return false;
	}
	return true;
}
Result EDrift(const Call& C)
{
	Real Ecc;
	int Spo = 0, Orbits = 0;
	std::string Error;
	if (!DriftArgs(C, 100000, &Ecc, &Spo, &Orbits, &Error))
	{
		return Result::Fail(Error);
	}
	const int Method = C.Int(3);
	if (Method < 0 || Method > 2)
	{
		return Result::Fail("method: 0 = Euler, 1 = leapfrog, 2 = RK4");
	}
	const std::vector<Real> D = EnergyDrift(Ecc, Spo, Orbits, Method);
	return Result::Record({{"max", MaxOf(D)}, {"first", D.front()}, {"last", D.back()}});
}
Result EDriftCmp(const Call& C)
{
	Real Ecc;
	int Spo = 0, Orbits = 0;
	std::string Error;
	if (!DriftArgs(C, 50000, &Ecc, &Spo, &Orbits, &Error))
	{
		return Result::Fail(Error);
	}
	const std::vector<Real> Lf = EnergyDrift(Ecc, Spo, Orbits, 1), Rk = EnergyDrift(Ecc, Spo, Orbits, 2), Eu = EnergyDrift(Ecc, Spo, Orbits, 0), Lf2 = EnergyDrift(Ecc, 2 * Spo, Orbits, 1);
	const Real MLf = MaxOf(Lf);
	return Result::Record({{"lf", MLf}, {"rk4", MaxOf(Rk)}, {"euler", MaxOf(Eu)}, {"q", MLf / MaxOf(Lf2)}, {"bound", MLf / Lf[0]}});
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
	const Real R1 = C.N(0), R2 = C.N(1), Rb = C.N(2), Mu = C.N(3);
	if (!(R1 > 0) || !(R2 > 0) || !(Mu > 0) || !(Rb >= gm::max(R1, R2)))
	{
		return Result::Fail("needs r1, r2, mu > 0 and rb >= max(r1, r2)");
	}
	const Real A1 = (R1 + Rb) / 2, A2 = (R2 + Rb) / 2;
	const Real Dv1 = gm::fabs(gm::sqrt(Mu * (2 / R1 - 1 / A1)) - gm::sqrt(Mu / R1));
	const Real Dv2 = gm::fabs(gm::sqrt(Mu * (2 / Rb - 1 / A2)) - gm::sqrt(Mu * (2 / Rb - 1 / A1)));
	const Real Dv3 = gm::fabs(gm::sqrt(Mu * (2 / R2 - 1 / A2)) - gm::sqrt(Mu / R2));
	return Result::Record({{"dv1", Dv1}, {"dv2", Dv2}, {"dv3", Dv3}, {"dv", Dv1 + Dv2 + Dv3}, {"tof", cst::Pi() * (gm::sqrt(A1 * A1 * A1 / Mu) + gm::sqrt(A2 * A2 * A2 / Mu))}});
}
Result Breakeven(const Call& C)
{
	const Real R = C.N(0);
	if (!(R > 1))
	{
		return Result::Fail("radius ratio must exceed 1");
	}
	const Real HT = HohTot(R);
	Real Lo = R, Hi = 1e12;
	if (BielTot(R, Hi) > HT)
	{
		return Result::Of(Real::Infinity());
	}
	if (BielTot(R, R * (1 + Real(1e-9))) < HT)
	{
		return Result::Of(R);
	}
	for (int i = 0; i < 200; ++i)
	{
		const Real M = gm::sqrt(Lo * Hi);
		if (M == Lo || M == Hi)
		{
			break;
		}
		if (BielTot(R, M) > HT)
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
	Real Best = 0, At = 0;
	for (int i = 0; i <= 3800; ++i)
	{
		const Real R = 2 + Real(i) / 100, V = HohTot(R);
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
	const Real R = C.N(0);
	const int N = C.Int(1);
	if (!(R > 0) || N < 2 || N > 200000)
	{
		return Result::Fail("needs R > 0 and 2 <= n <= 200000");
	}
	Real Best = Real::Infinity(), At = 0;
	const Real Step = Real(4) / (N - 1);
	for (int i = 0; i < N; ++i)
	{
		const Real Ra = R * (1 + Step * i), Cost = OneTangent(R, Ra);
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
constexpr int InsNs = 1 << 16;
constexpr int InsFix = 10;
constexpr int InsNk = 3000;

Real InsG() { return Real(9.80665); }
Real InsRe() { return Real(6.371e6); }   // ins.js uses the mean radius
Real InsDt() { return Real(0.1); }       // Allan records at 10 Hz
Real InsKdt() { return Real(0.1); }      // Kalman IMU step
const Real& InsWs()
{
	static const Real Ws = gm::sqrt(InsG() / InsRe());
	return Ws;
}
const Real& InsFlat()
{
	static const Real F = gm::sqrt(2 * gm::Ln2() / cst::Pi());   // IEEE Std 952: flicker floor = 0.664 B
	return F;
}

Real ModelAdev(const Real& Tau, const Real& N, const Real& B, const Real& K)
{
	const Real FB = InsFlat() * B;
	return gm::sqrt(N * N / Tau + FB * FB + K * K * Tau / 3);
}

// mulberry32 in exact 32-bit integer arithmetic (as ins.js rng()), Box-Muller with a cached spare.
struct InsRng
{
	uint32_t A = 1;
	bool HasSpare = false;
	Real Spare;
	Real Two32 = Real(4294967296ull);
	Real Tau = cst::Tau();

	explicit InsRng(uint32_t Seed) : A(Seed ? Seed : 1u) {}
	Real U()
	{
		A += 0x6D2B79F5u;
		uint32_t T = A;
		T = (T ^ (T >> 15)) * (T | 1u);
		T ^= T + (T ^ (T >> 7)) * (T | 61u);
		return Real(static_cast<unsigned long long>(T ^ (T >> 14))) / Two32;
	}
	Real Gauss()
	{
		if (HasSpare)
		{
			HasSpare = false;
			return Spare;
		}
		Real X;
		do
		{
			X = U();
		} while (X <= 0);
		const Real R = gm::sqrt(-2 * gm::log(X)), Th = Tau * U();
		Spare = R * gm::sin(Th);
		HasSpare = true;
		return R * gm::cos(Th);
	}
};

// Radix-2 in-place FFT, the same butterflies and twiddle recurrence as ins.js fft().
void InsFft(std::vector<Real>& Re, std::vector<Real>& Im, bool Inverse)
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
		const Real Ang = (Inverse ? 2 : -2) * cst::Pi() / Len, Wr = gm::cos(Ang), Wi = gm::sin(Ang);
		const int Half = Len >> 1;
		for (int i = 0; i < N; i += Len)
		{
			Real Cr = 1, Ci = 0;
			for (int k = 0; k < Half; ++k)
			{
				const size_t A = static_cast<size_t>(i + k), B = A + static_cast<size_t>(Half);
				const Real Xr = Re[B] * Cr - Im[B] * Ci, Xi = Re[B] * Ci + Im[B] * Cr;
				Re[B] = Re[A] - Xr;
				Im[B] = Im[A] - Xi;
				Re[A] += Xr;
				Im[A] += Xi;
				const Real T = Cr * Wr - Ci * Wi;
				Ci = Cr * Wi + Ci * Wr;
				Cr = T;
			}
		}
	}
	if (Inverse)
	{
		const Real Nr = N;
		for (size_t i = 0; i < Re.size(); ++i)
		{
			Re[i] /= Nr;
			Im[i] /= Nr;
		}
	}
}

std::vector<int> ClusterSizes(int N)
{
	std::vector<int> Out;
	const Real L = gm::log10(Real(N) / 8);
	for (int i = 0; i < 40; ++i)
	{
		const int M = gm::floor(gm::pow(Real(10), Real(i) / 39 * L) + Real(0.5)).ToInt();
		if (std::find(Out.begin(), Out.end(), M) == Out.end())
		{
			Out.push_back(M);
		}
	}
	return Out;
}

// Normalised overlapping Allan sums of the three unit processes of one seed (white rate w, rate random walk
// k, flicker f) and their cross terms, per cluster size: avar(N w + B f + K k) = N^2 WW + B^2 FF + K^2 KK
// + 2 N B WF + 2 N K WK + 2 B K FK, the same quadratic form as forming the combined record first.
struct AllanSums
{
	uint32_t Seed = 0;
	std::vector<int> Ms;
	std::vector<Real> WW, KK, FF, WF, WK, FK;
};

AllanSums ComputeAllanSums(uint32_t Seed)
{
	InsRng R(Seed);
	const size_t N = static_cast<size_t>(InsNs);
	const Real Dt = InsDt();
	std::vector<Real> W(N), K(N);
	const Real Sw = 1 / gm::sqrt(Dt), Sk = gm::sqrt(Dt);
	for (size_t i = 0; i < N; ++i)
	{
		W[i] = Sw * R.Gauss();
	}
	Real Cum = 0;
	for (size_t i = 0; i < N; ++i)
	{
		Cum += Sk * R.Gauss();
		K[i] = Cum;
	}
	// Flicker rate noise, two-sided PSD 1/(2 pi |f|): white noise shaped in the frequency domain.
	std::vector<Real> Re(N), Im(N);
	const Real Df = 1 / (InsNs * Dt), TwoPiDf = 2 * cst::Pi() * Df, NsDt = InsNs / Dt, Root2 = gm::sqrt(Real(2));
	for (int k = 1; k <= InsNs / 2; ++k)
	{
		const Real A = gm::sqrt(NsDt / (TwoPiDf * k));
		if (k == InsNs / 2)
		{
			Re[static_cast<size_t>(k)] = A * R.Gauss();
			continue;
		}
		const Real X = A * R.Gauss() / Root2, Y = A * R.Gauss() / Root2;
		Re[static_cast<size_t>(k)] = X;
		Im[static_cast<size_t>(k)] = Y;
		Re[static_cast<size_t>(InsNs - k)] = X;
		Im[static_cast<size_t>(InsNs - k)] = -Y;
	}
	InsFft(Re, Im, true);
	const std::vector<Real>& F = Re;

	std::vector<Real> Tw(N + 1), Tk(N + 1), Tf(N + 1);
	for (size_t i = 0; i < N; ++i)
	{
		Tw[i + 1] = Tw[i] + W[i] * Dt;
		Tk[i + 1] = Tk[i] + K[i] * Dt;
		Tf[i + 1] = Tf[i] + F[i] * Dt;
	}
	AllanSums S;
	S.Seed = Seed;
	S.Ms = ClusterSizes(InsNs);
	for (int M : S.Ms)
	{
		RealAccumulator Aww, Akk, Aff, Awf, Awk, Afk;
		const size_t Mm = static_cast<size_t>(M), Cn = N + 1 - 2 * Mm;
		for (size_t k = 0; k < Cn; ++k)
		{
			const Real Dw = Tw[k + 2 * Mm] - 2 * Tw[k + Mm] + Tw[k], Dk = Tk[k + 2 * Mm] - 2 * Tk[k + Mm] + Tk[k], Dff = Tf[k + 2 * Mm] - 2 * Tf[k + Mm] + Tf[k];
			Aww.AddProduct(Dw, Dw);
			Akk.AddProduct(Dk, Dk);
			Aff.AddProduct(Dff, Dff);
			Awf.AddProduct(Dw, Dff);
			Awk.AddProduct(Dw, Dk);
			Afk.AddProduct(Dff, Dk);
		}
		const Real Tau = M * Dt, Den = Real(static_cast<long long>(Cn)) * (2 * (Tau * Tau));
		S.WW.push_back(Aww.Result() / Den);
		S.KK.push_back(Akk.Result() / Den);
		S.FF.push_back(Aff.Result() / Den);
		S.WF.push_back(Awf.Result() / Den);
		S.WK.push_back(Awk.Result() / Den);
		S.FK.push_back(Afk.Result() / Den);
	}
	return S;
}

// One seed's records take seconds on the gears; like ins.js allanSet, keep the last few seeds' results.
AllanSums GetAllanSums(const Real& Seed)
{
	static std::vector<AllanSums> Cache;
	const uint32_t Key = static_cast<uint32_t>(static_cast<unsigned long long>((Real(1000) + Seed).ToLong()));
	for (const AllanSums& E : Cache)
	{
		if (E.Seed == Key)
		{
			return E;
		}
	}
	if (Cache.size() >= 8)
	{
		Cache.erase(Cache.begin());
	}
	Cache.push_back(ComputeAllanSums(Key));
	return Cache.back();
}

Real PolySlope(const std::vector<Real>& Xs, const std::vector<Real>& Ys)
{
	const Real N = static_cast<long long>(Xs.size());
	Real Mx = 0, My = 0;
	for (size_t i = 0; i < Xs.size(); ++i)
	{
		Mx += Xs[i];
		My += Ys[i];
	}
	Mx /= N;
	My /= N;
	Real Sxy = 0, Sxx = 0;
	for (size_t i = 0; i < Xs.size(); ++i)
	{
		Sxy += (Xs[i] - Mx) * (Ys[i] - My);
		Sxx += (Xs[i] - Mx) * (Xs[i] - Mx);
	}
	return Sxy / Sxx;
}

Real GeoMean(const std::vector<Real>& A)
{
	Real S = 0;
	for (const Real& X : A)
	{
		S += gm::log(X);
	}
	return gm::exp(S / Real(static_cast<long long>(A.size())));
}

Real SchulerClosed(const Real& T, const Real& B, const Real& Phi0, const Real& Eps)
{
	const Real Ws = InsWs(), C = B - InsG() * Phi0, Wt = Ws * T;
	return C / (Ws * Ws) * (1 - gm::cos(Wt)) - InsRe() * Eps * (T - gm::sin(Wt) / Ws);
}

using Q4 = std::array<Real, 4>;
using V3 = std::array<Real, 3>;
Q4 QMul(const Q4& A, const Q4& B)
{
	return {A[0] * B[0] - A[1] * B[1] - A[2] * B[2] - A[3] * B[3], A[0] * B[1] + A[1] * B[0] + A[2] * B[3] - A[3] * B[2], A[0] * B[2] - A[1] * B[3] + A[2] * B[0] + A[3] * B[1],
		A[0] * B[3] + A[1] * B[2] - A[2] * B[1] + A[3] * B[0]};
}
Q4 QExp(const V3& V)
{
	const Real A = gm::hypot(V[0], V[1], V[2]);
	if (A == 0)
	{
		return {1, 0, 0, 0};
	}
	const Real Half = A / 2, S = gm::sin(Half) / A;
	return {gm::cos(Half), S * V[0], S * V[1], S * V[2]};
}
V3 Cross3(const V3& A, const V3& B) { return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]}; }

// The coning motion with its angle-dependent factors worked out once.
struct Cone
{
	Real CosH, SinH, WSinA, TwoWSin2H, W;
	Cone(const Real& Al, const Real& InW) : W(InW)
	{
		const Real H = Al / 2;
		CosH = gm::cos(H);
		SinH = gm::sin(H);
		WSinA = W * gm::sin(Al);
		TwoWSin2H = 2 * W * SinH * SinH;
	}
	Q4 True(const Real& T) const
	{
		const Real Wt = W * T;
		return {CosH, SinH * gm::cos(Wt), SinH * gm::sin(Wt), 0};
	}
	V3 Omega(const Real& T) const
	{
		const Real Wt = W * T;
		return {-WSinA * gm::sin(Wt), WSinA * gm::cos(Wt), -TwoWSin2H};
	}
	// Gyro increment over [t0, t1] by three-point Gauss-Legendre quadrature.
	V3 DTheta(const Real& T0, const Real& T1) const
	{
		static const Real Node = gm::sqrt(Real(3) / 5);
		static const Real Gl[3][2] = {{-Node, Real(5) / 9}, {0, Real(8) / 9}, {Node, Real(5) / 9}};
		const Real Cm = (T0 + T1) / 2, H = (T1 - T0) / 2;
		V3 S = {0, 0, 0};
		for (const auto& P : Gl)
		{
			const V3 O = Omega(Cm + H * P[0]);
			const Real Wh = P[1] * H;
			for (size_t i = 0; i < 3; ++i)
			{
				S[i] += Wh * O[i];
			}
		}
		return S;
	}
};
Real ErrAngle(const Q4& Qt, const Q4& Q)
{
	const Q4 Qc = {Qt[0], -Qt[1], -Qt[2], -Qt[3]};
	const Q4 E = QMul(Qc, Q);
	return 2 * gm::asin(gm::min(Real(1), gm::hypot(E[1], E[2], E[3])));
}

using M3 = std::array<Real, 9>;
M3 MMul3(const M3& A, const M3& B)
{
	M3 C;
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			C[3 * i + j] = A[3 * i] * B[j] + A[3 * i + 1] * B[3 + j] + A[3 * i + 2] * B[6 + j];
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
	Real Sa;
};
// ins.js kfModel: state [position, velocity, accelerometer bias], IMU at 10 Hz.
KfModel MakeKf(const Real& Qa, const Real& Sbw)
{
	KfModel M;
	const Real Dt = InsKdt();
	M.Sa = Qa / gm::sqrt(Dt);
	M.F = {1, Dt, -Dt * Dt / 2, 0, 1, -Dt, 0, 0, 1};
	const Real Gm[3] = {Dt * Dt / 2, Dt, 0};
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			M.Q[3 * i + j] = Gm[i] * Gm[j] * M.Sa * M.Sa;
		}
	}
	M.Q[8] += Sbw * Sbw * Dt;
	return M;
}
Real OutageClosed(const M3& P, const Real& L, const Real& Sa, const Real& Sbw)
{
	const Real Cv[3] = {1, L, -L * L / 2};
	Real S = 0;
	for (size_t i = 0; i < 3; ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			S += Cv[i] * P[3 * i + j] * Cv[j];
		}
	}
	const Real L3 = L * L * L;
	return gm::sqrt(S + Sa * Sa * InsKdt() * L3 / 3 + Sbw * Sbw * L3 * L * L / 20);
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
	const Real N = C.N(0), B = C.N(1), K = C.N(2);
	if (!(N > 0) || !(K > 0))
	{
		return Result::Fail("needs N > 0 and K > 0 (otherwise the minimum is at tau = 0 or infinity)");
	}
	const Real Tau = gm::sqrt(Real(3)) * N / K;
	return Result::Record({{"tau", Tau}, {"adev", ModelAdev(Tau, N, B, K)}});
}
// Simpson's rule on sin^4 x / x^3. sin(i h) is carried by angle addition between exact gm::sin/cos
// anchors every 256 nodes (drift below 1e-15), which saves a full sine per node.
Result FlickerQuad(const Call& C)
{
	const Real X = C.N(0);
	const int N = C.Int(1);
	if (!(X > 0) || N < 2 || N % 2 != 0 || N > 1000000)
	{
		return Result::Fail("needs X > 0 and an even n in 2..1000000");
	}
	const Real H = X / N, SinH = gm::sin(H), CosH = gm::cos(H), Small = 1e-6;
	auto Fn = [&Small](const Real& Xv, const Real& S) {
		if (Xv < Small)
		{
			return Xv;
		}
		const Real S2 = S * S;
		return S2 * S2 / (Xv * Xv * Xv);
	};
	RealAccumulator Odd, Even;
	Real Sn = 0, Cn = 1;
	for (int i = 1; i < N; ++i)
	{
		const Real Xi = Real(i) * H;
		if ((i - 1) % 256 == 0)
		{
			Sn = gm::sin(Xi);
			Cn = gm::cos(Xi);
		}
		else
		{
			const Real S1 = Sn * CosH + Cn * SinH;
			Cn = Cn * CosH - Sn * SinH;
			Sn = S1;
		}
		if (i % 2)
		{
			Odd.Add(Fn(Xi, Sn));
		}
		else
		{
			Even.Add(Fn(Xi, Sn));
		}
	}
	const Real S = Fn(Real(0), Real(0)) + Fn(X, gm::sin(X)) + 4 * Odd.Result() + 2 * Even.Result();
	return Result::Of(S * H / 3 + Real(3) / 8 / (2 * X * X));
}
Result AllanFit(const Call& C)
{
	const AllanSums S = GetAllanSums(C.N(0));
	const Real Dt = InsDt(), Tt = Dt * InsNs, WhiteMax = Tt / 100, RrwMax = Tt / 50;
	std::vector<Real> Wx, Wy, Wn, Kx, Ky, Kn, Fl;
	for (size_t i = 0; i < S.Ms.size(); ++i)
	{
		const Real T = S.Ms[i] * Dt, Aw = gm::sqrt(S.WW[i]), Ak = gm::sqrt(S.KK[i]), Af = gm::sqrt(S.FF[i]);
		if (T <= WhiteMax)
		{
			Wx.push_back(gm::log(T));
			Wy.push_back(gm::log(Aw));
			Wn.push_back(Aw * gm::sqrt(T));
		}
		if (T >= 1 && T <= RrwMax)
		{
			Kx.push_back(gm::log(T));
			Ky.push_back(gm::log(Ak));
			Kn.push_back(Ak / gm::sqrt(T / 3));
		}
		if (T >= 2 && T <= RrwMax)
		{
			Fl.push_back(Af);
		}
	}
	std::sort(Fl.begin(), Fl.end());
	const size_t Nf = Fl.size();
	const Real Flat = Nf % 2 ? Fl[(Nf - 1) / 2] : (Fl[Nf / 2 - 1] + Fl[Nf / 2]) / 2;
	return Result::Record({{"wslope", PolySlope(Wx, Wy)}, {"nhat", GeoMean(Wn)}, {"kslope", PolySlope(Kx, Ky)}, {"khat", GeoMean(Kn)}, {"flat", Flat}});
}
Result AllanSim(const Call& C)
{
	const Real N = C.N(0), B = C.N(1), K = C.N(2);
	const AllanSums S = GetAllanSums(C.N(3));
	const Real Dt = InsDt(), NN = N * N, BB = B * B, KK = K * K, NB = 2 * N * B, NK = 2 * N * K, BK = 2 * B * K;
	Matrix Out(static_cast<int>(S.Ms.size()), 2);
	for (size_t i = 0; i < S.Ms.size(); ++i)
	{
		const Real Avar = NN * S.WW[i] + BB * S.FF[i] + KK * S.KK[i] + NB * S.WF[i] + NK * S.WK[i] + BK * S.FK[i];
		Out(static_cast<int>(i), 0) = S.Ms[i] * Dt;
		Out(static_cast<int>(i), 1) = gm::sqrt(gm::max(Real(0), Avar));
	}
	return Result::Of(Out);
}

Result SchulerPeriod(const Call& C)
{
	if (!(C.N(0) > 0) || !(C.N(1) > 0))
	{
		return Result::Fail("needs R > 0 and g > 0");
	}
	return Result::Of(cst::Tau() * gm::sqrt(C.N(0) / C.N(1)));
}
Result SchulerX(const Call& C) { return Result::Of(SchulerClosed(C.N(0), C.N(1), C.N(2), C.N(3))); }
Result NaiveDrift(const Call& C)
{
	const Real T = C.N(0), B = C.N(1), Phi0 = C.N(2), Eps = C.N(3);
	return Result::Of((B - InsG() * Phi0) * T * T / 2 - InsG() * Eps * T * T * T / 6);
}
Result SchulerSim(const Call& C)
{
	const Real B = C.N(0), Phi0 = C.N(1), Eps = C.N(2), Tend = C.N(3), H = C.N(4);
	if (!(H > 0) || !(Tend > 0) || Tend / H > 100000)
	{
		return Result::Fail("needs tend > 0, h > 0 and at most 100000 steps");
	}
	const int N = gm::floor(Tend / H + Real(0.5)).ToInt();
	const Real G = InsG(), InvRe = 1 / InsRe(), H2 = H / 2, H6 = H / 6;
	Real S0 = 0, S1 = 0, S2 = Phi0, DMax = 0;
	auto Fn = [&](const Real& V, const Real& P) { return V3{V, B - G * P, V * InvRe + Eps}; };
	std::vector<Real> Acc(static_cast<size_t>(N) + 1);
	Acc[0] = B - G * S2;
	for (int i = 1; i <= N; ++i)
	{
		const V3 K1 = Fn(S1, S2), K2 = Fn(S1 + H2 * K1[1], S2 + H2 * K1[2]), K3 = Fn(S1 + H2 * K2[1], S2 + H2 * K2[2]), K4 = Fn(S1 + H * K3[1], S2 + H * K3[2]);
		S0 += H6 * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]);
		S1 += H6 * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]);
		S2 += H6 * (K1[2] + 2 * K2[2] + 2 * K3[2] + K4[2]);
		Acc[static_cast<size_t>(i)] = B - G * S2;
		DMax = gm::max(DMax, gm::fabs(S0 - SchulerClosed(H * i, B, Phi0, Eps)));
	}
	// Period from the zero crossings of the velocity-error rate.
	std::vector<Real> Z;
	for (size_t i = 0; i + 1 < Acc.size(); ++i)
	{
		if (Acc[i] * Acc[i + 1] < 0)
		{
			Z.push_back(H * static_cast<long long>(i) - Acc[i] * H / (Acc[i + 1] - Acc[i]));
		}
	}
	const Real Period = Z.size() < 2 ? Real::NaN() : 2 * (Z.back() - Z.front()) / Real(static_cast<long long>(Z.size() - 1));
	return Result::Record({{"x", S0}, {"dmax", DMax}, {"period", Period}});
}

Result Coning0(const Call& C)
{
	const Real Al = C.N(0), W = C.N(1), H = C.N(2), X = W * H, S = gm::sin(Al);
	if (X == 0)
	{
		return Result::Fail("needs W h != 0");
	}
	return Result::Of(Real(0.5) * W * S * S * (1 - gm::sin(X) / X));
}
Result Coning2(const Call& C)
{
	const Real Al = C.N(0), W = C.N(1), H = C.N(2), X = W * H, S = gm::sin(Al), X2 = X * X;
	return Result::Of(W * S * S * X2 * X2 / 960);
}
Result ConingRun(const Call& C)
{
	const Real Al = C.N(0), W = C.N(1), Rate = C.N(2), Tend = C.N(3);
	if (!(Rate > 0) || !(Tend > 0) || Tend * Rate > 10000)
	{
		return Result::Fail("needs rate > 0, tend > 0 and at most 10000 updates");
	}
	const Real H = 1 / Rate, Hh = H / 2, Two3 = Real(2) / 3;
	const int N = gm::floor(Tend * Rate + Real(0.5)).ToInt();
	if (N < 2)
	{
		return Result::Fail("needs at least two updates");
	}
	const Cone Cn(Al, W);
	Q4 Q0 = Cn.True(0), Q2 = Q0;
	RealAccumulator S1, S2, S3, S4, S5, S6;
	Real X0, X2;
	for (int k = 0; k < N; ++k)
	{
		const Real T0 = H * k;
		const V3 A = Cn.DTheta(T0, T0 + Hh), B = Cn.DTheta(T0 + Hh, T0 + H), Cx = Cross3(A, B);
		Q0 = QMul(Q0, QExp({A[0] + B[0], A[1] + B[1], A[2] + B[2]}));
		Q2 = QMul(Q2, QExp({A[0] + B[0] + Two3 * Cx[0], A[1] + B[1] + Two3 * Cx[1], A[2] + B[2] + Two3 * Cx[2]}));
		const Real T = H * (k + 1);
		const Q4 Qt = Cn.True(T);
		X0 = ErrAngle(Qt, Q0);
		X2 = ErrAngle(Qt, Q2);
		S1.Add(T);
		S2.AddProduct(T, T);
		S3.Add(X0);
		S4.AddProduct(T, X0);
		S5.Add(X2);
		S6.AddProduct(T, X2);
	}
	// Least-squares slopes of both error histories against time.
	const Real Sn = N, St = S1.Result(), Stt = S2.Result(), Den = Sn * Stt - St * St;
	return Result::Record({{"slope0", (Sn * S4.Result() - St * S3.Result()) / Den}, {"slope2", (Sn * S6.Result() - St * S5.Result()) / Den}, {"e0", X0}, {"e2", X2}});
}
Result Bortz(const Call& C)
{
	const Real Hh = C.N(0);
	const int Steps = C.Int(1);
	if (!(Hh > 0) || Steps < 1 || Steps > 20000)
	{
		return Result::Fail("needs h > 0 and 1 <= steps <= 20000");
	}
	const V3 A = {0.3, -0.2, 0.5}, B = {1.1, 0.4, -0.7};
	V3 D1, D2;
	for (size_t i = 0; i < 3; ++i)
	{
		D1[i] = A[i] * Hh / 2 + B[i] * Hh * Hh / 8;
		D2[i] = A[i] * Hh / 2 + 3 * B[i] * Hh * Hh / 8;
	}
	auto Wt = [&A, &B](const Real& T) { return V3{A[0] + B[0] * T, A[1] + B[1] * T, A[2] + B[2] * T}; };
	const Real Small = 1e-8, Twelfth = Real(1) / 12;
	auto Rate = [&Small, &Twelfth](const V3& P, const V3& O) {
		const Real N = gm::hypot(P[0], P[1], P[2]), Cc = N > Small ? (1 - N * gm::sin(N) / (2 * (1 - gm::cos(N)))) / (N * N) : Twelfth;
		const V3 X1 = Cross3(P, O), X2 = Cross3(P, X1);
		V3 Out;
		for (size_t i = 0; i < 3; ++i)
		{
			Out[i] = O[i] + Real(0.5) * X1[i] + Cc * X2[i];
		}
		return Out;
	};
	auto Add = [](const V3& U, const V3& V, const Real& S) { return V3{U[0] + S * V[0], U[1] + S * V[1], U[2] + S * V[2]}; };
	V3 P = {0, 0, 0};
	const Real Hs = Hh / Steps, Hs2 = Hs / 2, Hs6 = Hs / 6;
	for (int i = 0; i < Steps; ++i)
	{
		const Real T = Hs * i, Tm = T + Hs2;
		const V3 Wm = Wt(Tm);
		const V3 K1 = Rate(P, Wt(T)), K2 = Rate(Add(P, K1, Hs2), Wm), K3 = Rate(Add(P, K2, Hs2), Wm), K4 = Rate(Add(P, K3, Hs), Wt(T + Hs));
		for (size_t j = 0; j < 3; ++j)
		{
			P[j] = P[j] + Hs6 * (K1[j] + 2 * K2[j] + 2 * K3[j] + K4[j]);
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
	const Real Qa = C.N(0), Sb0 = C.N(1), Sbw = C.N(2), Sg = C.N(3), L = C.N(4), Out0 = C.N(5), Dt = InsKdt();
	if (!(Qa >= 0) || !(Sb0 > 0) || !(Sbw >= 0) || !(Sg > 0) || !(L >= 0))
	{
		return Result::Fail("needs qa, sbw, L >= 0 and sb0, sg > 0");
	}
	if (!(Out0 > 0) || !(Out0 + L <= 300))
	{
		return Result::Fail("the outage [out0, out0 + L] must lie inside the 300 s run");
	}
	const int Kb = gm::floor(Out0 / Dt + Real(0.5)).ToInt() - 1, Ke = gm::floor((Out0 + L) / Dt + Real(0.5)).ToInt() - 1;
	if (Kb < 0 || Ke > InsNk)
	{
		return Result::Fail("the outage [out0, out0 + L] must lie inside the 300 s run");
	}
	const KfModel M = MakeKf(Qa, Sbw);
	M3 P = {100, 0, 0, 0, 1, 0, 0, 0, Sb0 * Sb0};
	M3 P1400 = P, PKb = P, PKe = P;
	const Real R = Sg * Sg, Slack = 1e-9, OutEnd = Out0 + L;
	for (int k = 1; k <= InsNk; ++k)
	{
		P = MAdd3(Fapt(M.F, P), M.Q);
		const Real Tk = Dt * k + Slack;
		const bool Fix = k % InsFix == 0 && !(Out0 <= Tk && Tk < OutEnd);
		if (Fix)
		{
			const Real S = P[0] + R, K[3] = {P[0] / S, P[3] / S, P[6] / S};
			const M3 Ik = {1 - K[0], 0, 0, -K[1], 1, 0, -K[2], 0, 1};
			M3 Kr;
			for (size_t i = 0; i < 3; ++i)
			{
				for (size_t j = 0; j < 3; ++j)
				{
					Kr[3 * i + j] = K[i] * K[j] * R;
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
	return Result::Record({{"sp140", gm::sqrt(P1400[0])}, {"sbias", gm::sqrt(PKb[8])}, {"sp0", gm::sqrt(PKb[0])}, {"spend", gm::sqrt(PKe[0])},
		{"closed", OutageClosed(PKb, Dt * (Ke - Kb), M.Sa, Sbw)}});
}
Result KfSteady(const Call& C)
{
	const Real Qa = C.N(0), Sbw = C.N(1), Sg = C.N(2);
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
	Matrix A(3, 3), Qm(3, 3), Hb(3, 1), Rm = Matrix::Scalar(Sg * Sg), P, K;
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
	const Real S = P(0, 0) + Sg * Sg, Kp = P(0, 0) / S, Kv = P(1, 0) / S, Kbias = P(2, 0) / S;
	const Real Pp = P(0, 0) - P(0, 0) * P(0, 0) / S, Pv = P(1, 1) - P(1, 0) * P(1, 0) / S, Pb = P(2, 2) - P(2, 0) * P(2, 0) / S;
	return Result::Record({{"sp", gm::sqrt(Pp)}, {"sv", gm::sqrt(Pv)}, {"sb", gm::sqrt(Pb)}, {"sppri", gm::sqrt(P(0, 0))}, {"kp", Kp}, {"kv", Kv}, {"kb", Kbias}});
}
Result InsOutage(const Call& C)
{
	const Matrix& Pm = C.M(0);
	const Real L = C.N(1), Qa = C.N(2), Sbw = C.N(3);
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
	return Result::Of(OutageClosed(P, L, Qa / gm::sqrt(InsKdt()), Sbw));
}
} // namespace orbits_ins_detail

void RegisterOrbitsLibrary(Registry& R)
{
	{
		using namespace orbits_detail;
		const char* O = "orbits";
		const Real Mu = cst::MuEarth();
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
		Def(R, "twobody", O, "orbits/orbits.js propagate", "RK4 integration of r'' = -mu r/|r|^3, n steps per local period (at most 100000 steps)", "record x, y, vx, vy, steps",
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
		Def(R, "cwnonlin", O, "orbits/orbits.js relNonlinear (s4)", "nonlinear relative motion about a circular orbit, RK4 (steps <= 100000)", "record x, y (m), vx, vy (m/s), range (m)",
			{Num("rref", "m"), Mat("X", "m, m/s [x, y, vx, vy]"), Num("t", "s"), Num("steps", "-", 4000), Num("mu", "m^3/s^2", Mu)}, &CwNonlin);
		Def(R, "rendezvous", O, "orbits/orbits.js rendezvous (s4 two-impulse)", "two-impulse CW rendezvous to the origin in time T", "record v0x, v0y, dv1x, dv1y, dv2x, dv2y, dv1, dv2, dv (m/s), n",
			{Num("rref", "m"), Mat("X", "m, m/s [x, y, vx, vy]"), Num("T", "s"), Num("mu", "m^3/s^2", Mu)}, &Rendezvous);
		Def(R, "cwshoot", O, "orbits/orbits.js cwShoot (s4)", "first rendezvous burn corrected by shooting in the nonlinear dynamics (steps <= 5000)", "record v0x, v0y (m/s), iters, corr (m/s)",
			{Num("rref", "m"), Mat("X", "m, m/s [x, y, vx, vy]"), Num("T", "s"), Num("steps", "-", 4000), Num("mu", "m^3/s^2", Mu)}, &CwShoot);
		Def(R, "edrift", O, "orbits/orbits.js drift (s5 energy drift)", "max relative energy error, a = 10000 km; method 0 Euler, 1 leapfrog, 2 RK4 (spo * orbits <= 100000)",
			"record max, first, last (per-orbit max)", {Num("e", "-", 0.3), Num("spo", "steps/orbit", 500), Num("orbits", "-", 50), Num("method", "-", 1)}, &EDrift);
		Def(R, "edriftcmp", O, "orbits/ s5 energy drift", "leapfrog, RK4 and Euler energy drift; q = leapfrog(spo)/leapfrog(2 spo) (spo * orbits <= 50000)", "record lf, rk4, euler, q, bound",
			{Num("e", "-", 0.3), Num("spo", "steps/orbit", 500), Num("orbits", "-", 50)}, &EDriftCmp);
		Def(R, "hohtot", O, "orbits/orbits.js hohTot (s6 bi-elliptic)", "Hohmann total delta-v in units of the inner circular speed", "-", {Num("R", "r2/r1")}, &HohTotFn);
		Def(R, "bieltot", O, "orbits/orbits.js bielTot (s6 bi-elliptic)", "bi-elliptic total delta-v in units of the inner circular speed", "-", {Num("R", "r2/r1"), Num("Rb", "rb/r1")}, &BielTotFn);
		Def(R, "bielinf", O, "orbits/orbits.js bielInf (s6 bi-elliptic)", "bi-elliptic total as rb -> infinity, (sqrt2 - 1)(1 + 1/sqrt R)", "-", {Num("R", "r2/r1")}, &BielInfFn);
		Def(R, "bielthresh", O, "orbits/orbits.js thresholds (s6 bi-elliptic)", "ratio R = r2/r1 above which bi-elliptic (rb -> inf) beats Hohmann, 11.93877", "record closed (rb = inf), finite (rb = 1e12)", {},
			&BielThresh);
		Def(R, "bielliptic", O, "orbits/ s6 bi-elliptic", "bi-elliptic transfer r1 -> rb -> r2 between circular orbits", "record dv1, dv2, dv3, dv (m/s), tof (s)",
			{Num("r1", "m"), Num("r2", "m"), Num("rb", "m"), Num("mu", "m^3/s^2", Mu)}, &Bielliptic);
		Def(R, "breakeven", O, "orbits/orbits.js breakeven (s6 bi-elliptic)", "smallest rb/r1 at which bi-elliptic matches Hohmann (inf if never)", "rb/r1", {Num("R", "r2/r1")}, &Breakeven);
		Def(R, "hohpeak", O, "orbits/orbits.js hohPeak (s6 bi-elliptic)", "costliest Hohmann ratio (15.58; beyond it every bi-elliptic wins)", "record R, v", {}, &HohPeak);
		Def(R, "onetangent", O, "orbits/orbits.js oneTangent (s6)", "one-tangent-burn transfer cost, apoapsis ra", "-", {Num("R", "r2/r1"), Num("ra", "ra/r1")}, &OneTangentFn);
		Def(R, "onetanscan", O, "orbits/orbits.js oneTangentScan (s6)", "cheapest one-tangent transfer over ra in [R, 5R] (n <= 200000)", "record best, ratio (ra/r2)", {Num("R", "r2/r1"), Num("n", "-", 4001)},
			&OneTanScan);
	}
	{
		using namespace orbits_ins_detail;
		const char* I = "ins";
		const Real B0 = Real(0.1e-3) * cst::G0(), Eps0 = Real(0.01) * cst::Pi() / 180 / 3600;
		const Real Qa0 = Real(300e-6) * cst::G0(), Sb00 = Real(5e-3) * cst::G0(), Sbw0 = Real(20e-6) * cst::G0();
		Def(R, "allandev", I, "ins/ins.js modelAdev (s1 Allan)", "model Allan deviation sqrt(N^2/tau + (0.664 B)^2 + K^2 tau/3)", "rate unit",
			{Num("tau", "s"), Num("N", "u sqrt(s)"), Num("B", "u", 0), Num("K", "u/sqrt(s)", 0)}, &AllanDev);
		Def(R, "allanmin", I, "ins/ s1 Allan", "bottom of the model Allan curve, tau = sqrt3 N/K", "record tau (s), adev", {Num("N", "u sqrt(s)"), Num("B", "u"), Num("K", "u/sqrt(s)")}, &AllanMin);
		Def(R, "flickerquad", I, "ins/ins.js flickerQuad (s1 Allan)", "int_0^inf sin^4 x / x^3 dx (= ln 2) by Simpson plus tail (n <= 1000000)", "-", {Num("X", "-", 400), Num("n", "-", 400000)},
			&FlickerQuad);
		Def(R, "allanfit", I, "ins/ins.js allanSet, fits (s1 Allan)", "seeded unit noise records: fitted slopes, N, K and the flicker floor (a new seed takes seconds)",
			"record wslope, nhat, kslope, khat, flat", {Num("seed", "-", 1)}, &AllanFit);
		Def(R, "allansim", I, "ins/ins.js combinedAdev (s1 Allan)", "Allan deviation of a seeded simulated record N w + B f + K k (a new seed takes seconds)", "m x 2 matrix [tau (s), adev]",
			{Num("N", "u sqrt(s)"), Num("B", "u"), Num("K", "u/sqrt(s)"), Num("seed", "-", 1)}, &AllanSim);
		Def(R, "schuler", I, "ins/ins.js PERIOD (s2 Schuler)", "Schuler period 2 pi sqrt(R/g) (84.41 min)", "s", {Num("R", "m", InsRe()), Num("g", "m/s^2", InsG())}, &SchulerPeriod);
		Def(R, "schulerx", I, "ins/ins.js schulerClosed (s2 Schuler)", "closed-form Schuler-loop position error", "m",
			{Num("t", "s"), Num("b", "m/s^2", B0), Num("phi0", "rad", 0), Num("eps", "rad/s", Eps0)}, &SchulerX);
		Def(R, "naivedrift", I, "ins/ins.js naiveDrift (s2 Schuler)", "flat-earth position drift without Schuler feedback", "m",
			{Num("t", "s"), Num("b", "m/s^2", B0), Num("phi0", "rad", 0), Num("eps", "rad/s", Eps0)}, &NaiveDrift);
		Def(R, "schulersim", I, "ins/ins.js schulerRK4, measuredPeriod (s2)", "RK4 Schuler error loop: end error, gap to the closed form, measured period (tend/h <= 100000)",
			"record x (m), dmax (m), period (s)", {Num("b", "m/s^2", B0), Num("phi0", "rad", 0), Num("eps", "rad/s", Eps0), Num("tend", "s", 10800), Num("h", "s", 1)}, &SchulerSim);
		Def(R, "coning0", I, "ins/ins.js drift0 (s3 coning)", "coning drift with no correction", "rad/s", {Num("alpha", "rad"), Num("W", "rad/s"), Num("h", "s")}, &Coning0);
		Def(R, "coning2", I, "ins/ins.js drift2 (s3 coning)", "residual coning drift after the two-sample correction", "rad/s", {Num("alpha", "rad"), Num("W", "rad/s"), Num("h", "s")}, &Coning2);
		Def(R, "coningrun", I, "ins/ins.js coningRun (s3 coning)", "integrate coning motion uncorrected and with the 2/3 two-sample term (rate * tend <= 10000)",
			"record slope0, slope2 (rad/s), e0, e2 (rad)", {Num("alpha", "rad"), Num("W", "rad/s"), Num("rate", "Hz", 200), Num("tend", "s", 2)}, &ConingRun);
		Def(R, "bortz", I, "ins/ins.js bortzCoefficient (s3 coning)", "coefficient of d1 x d2 from Bortz's equation (-> 2/3; steps <= 20000)", "-", {Num("h", "s", 1e-2), Num("steps", "-", 2000)}, &Bortz);
		Def(R, "insgnss", I, "ins/ins.js kfModel, kfCov, outageClosed (s4)", "1-D INS/GNSS Kalman covariance, 1 Hz fixes, outage from out0 for L", "record sp140, sbias, sp0, spend, closed",
			{Num("qa", "m/s^2/sqrt(Hz)", Qa0), Num("sb0", "m/s^2", Sb00), Num("sbw", "m/s^2/sqrt(s)", Sbw0), Num("sg", "m", 3), Num("L", "s", 60), Num("out0", "s", 150)}, &InsGnss);
		Def(R, "kfsteady", I, "ins/ins.js riccati1s (s4 Kalman)", "steady-state 1-D INS/GNSS filter by the DARE", "record sp, sv, sb (posterior sigmas), sppri (m), kp, kv, kb",
			{Num("qa", "m/s^2/sqrt(Hz)", Qa0), Num("sbw", "m/s^2/sqrt(s)", Sbw0), Num("sg", "m", 3)}, &KfSteady);
		Def(R, "insoutage", I, "ins/ins.js outageClosed (s4 Kalman)", "position sigma after an outage of L seconds from covariance P", "m",
			{Mat("P", "3 x 3 [m, m/s, m/s^2]"), Num("L", "s"), Num("qa", "m/s^2/sqrt(Hz)", Qa0), Num("sbw", "m/s^2/sqrt(s)", Sbw0)}, &InsOutage);
	}
}
} // namespace mcalc
