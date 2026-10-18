// Procedural keyboard parts. Dimensions and their sources:
//
// Keycaps, Cherry profile. KeyV2 (github.com/rsheldiii/KeyV2, src/key_profiles/cherry.scad, "based off
// GMK keycap set"): 18.16 mm base, 11.85 x 14.64 mm top (the top is deeper than it is wide and skewed
// toward the back), 0.65 mm cylindrical dish, 0.6 mm stem inset, 5.5 mm round stem socket, and row
// tilts that step by 2.5, 2.5 and 6.5 degrees from R1 to R4. Row heights, measured from the skirt to
// the back edge of the top, are the commonly quoted GMK Cherry values: R1 9.4, R2 8.0, R3 7.4,
// R4 7.0 mm. The tilts used here keep KeyV2's steps with R1 leaning 7 degrees toward the typist and R4
// leaning 4.5 degrees away, so R3 is the low point of the sculpt and R4's front edge stands up.
//
// Switch, Cherry MX. Cherry MX datasheet and the ai03 plate-design wiki: 15.6 x 15.6 mm housing,
// 14.0 x 14.0 mm plate cut-out, 1.5 mm plate, 5.0 mm from plate top to PCB top, 1.6 mm PCB,
// 4.0 mm travel with 2.0 mm pre-travel, 4.0 mm centre post, pins at (-3.81, 2.54) and (2.54, 5.08) mm
// leaving the housing by 3.3 mm. Deskthority wiki: 15.2 mm from the PCB to the top of the stem, a
// cross about 3.7 mm tall with 4.0 mm arms 1.1-1.3 mm thick (the cap cut-out is 4.0 x 1.17 mm).
// So the slider's top (the cap mount) is 15.2 - 5.0 - 3.6 = 6.6 mm above the plate.
//
// Stabilizers: Costar/Cherry cut-outs 3.3 x 14 mm, wire spacing 11.938 mm per side for 2-2.75u,
// 19.05 mm for 3u, 47.625 mm for 6u, 50 mm for 6.25u and 57.15 mm for 7u.
#include "KeyboardParts.h"

#include <algorithm>
#include <cmath>

namespace mgeo
{
namespace keyboardparts_detail
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

constexpr double kMxPlateToPcb = 5.0;
constexpr double kMxCrossArm = 4.0;
constexpr double kMxCrossWidth = 1.17;
constexpr double kMxCrossHeight = 3.7;
constexpr double kMxStemClearance = 0.05;
constexpr double kMxPin1X = -3.81, kMxPin1Y = 2.54, kMxPin2X = 2.54, kMxPin2Y = 5.08;
constexpr double kMxPinBelowHousing = 3.3;
constexpr double kMxCutout = 14.0;

constexpr double kCapStemInset = 0.6;
constexpr double kCapSocketRadius = 2.75;
constexpr double kCapSocketDepth = 3.9;
constexpr double kCapFrontInsetScale = 1.1;
constexpr double kCapBackInsetScale = 0.45;
constexpr double kCapWallBulge = 0.1;
constexpr double kSpacebarCrown = 0.35;

constexpr double kStabSlotW = 3.3;
constexpr double kStabSlotH = 14.0;

struct Pt2
{
	double X = 0.0, Y = 0.0;
};

struct Loop2
{
	std::vector<Pt2> P;
	std::vector<Pt2> N;   // outward unit normals
};

// Points inserted on the straight edges of a rounded rectangle, as fractions in (-1, 1) of the
// straight part: Top and Bottom along x, Left and Right along y.
struct Splits
{
	std::vector<double> Top, Bottom, Left, Right;
};

struct Frame3
{
	Vec3 O, U = Vec3(1, 0, 0), V = Vec3(0, 1, 0), W = Vec3(0, 0, 1);
	Vec3 At(double A, double B, double C) const { return O + U * A + V * B + W * C; }
};

double Lerp(double A, double B, double T)
{
	return A + (B - A) * T;
}

std::vector<double> Uniform(int Count)
{
	std::vector<double> F;
	for (int k = 1; k <= Count; ++k)
	{
		F.push_back(-1.0 + 2.0 * k / (Count + 1));
	}
	return F;
}

// Sorted fractions: the Exact ones are kept as given and any other closer than MinGap to one is dropped.
std::vector<double> MergeFractions(const std::vector<double>& Base, const std::vector<double>& Exact, double MinGap)
{
	std::vector<double> Out;
	for (double F : Exact)
	{
		if (std::fabs(F) < 0.995)
		{
			Out.push_back(F);
		}
	}
	std::sort(Out.begin(), Out.end());
	Out.erase(std::unique(Out.begin(), Out.end(), [](double A, double B) { return std::fabs(A - B) < 1e-9; }), Out.end());
	const std::vector<double> Fixed = Out;
	std::vector<double> Sorted = Base;
	std::sort(Sorted.begin(), Sorted.end());
	for (double F : Sorted)
	{
		if (std::fabs(F) >= 0.995)
		{
			continue;
		}
		bool Clear = true;
		for (double G : Out)
		{
			if (std::fabs(F - G) < MinGap)
			{
				Clear = false;
				break;
			}
		}
		if (Clear)
		{
			Out.push_back(F);
		}
	}
	std::sort(Out.begin(), Out.end());
	return Out;
}

Splits SameSplits(const std::vector<double>& XEdges, const std::vector<double>& YEdges)
{
	Splits S;
	S.Top = S.Bottom = XEdges;
	S.Left = S.Right = YEdges;
	return S;
}

// Counter-clockwise rounded rectangle starting at angle 0 on the top-right corner.
Loop2 RoundedRect(double Cx, double Cy, double Hw, double Hd, double R, int Seg, const Splits& S)
{
	Loop2 L;
	R = std::max(1e-3, std::min(R, std::min(Hw, Hd) - 1e-3));
	Seg = std::max(Seg, 1);
	const double Ex = Hw - R, Ey = Hd - R;
	auto Arc = [&L, R, Seg](double Ax, double Ay, double A0) {
		for (int s = 0; s <= Seg; ++s)
		{
			const double A = A0 + 0.5 * kPi * s / Seg;
			L.P.push_back(Pt2{Ax + R * std::cos(A), Ay + R * std::sin(A)});
			L.N.push_back(Pt2{std::cos(A), std::sin(A)});
		}
	};
	Arc(Cx + Ex, Cy + Ey, 0.0);
	for (auto It = S.Top.rbegin(); It != S.Top.rend(); ++It)
	{
		L.P.push_back(Pt2{Cx + *It * Ex, Cy + Hd});
		L.N.push_back(Pt2{0.0, 1.0});
	}
	Arc(Cx - Ex, Cy + Ey, 0.5 * kPi);
	for (auto It = S.Left.rbegin(); It != S.Left.rend(); ++It)
	{
		L.P.push_back(Pt2{Cx - Hw, Cy + *It * Ey});
		L.N.push_back(Pt2{-1.0, 0.0});
	}
	Arc(Cx - Ex, Cy - Ey, kPi);
	for (double F : S.Bottom)
	{
		L.P.push_back(Pt2{Cx + F * Ex, Cy - Hd});
		L.N.push_back(Pt2{0.0, -1.0});
	}
	Arc(Cx + Ex, Cy - Ey, 1.5 * kPi);
	for (double F : S.Right)
	{
		L.P.push_back(Pt2{Cx + Hw, Cy + F * Ey});
		L.N.push_back(Pt2{1.0, 0.0});
	}
	return L;
}

Loop2 CircleLoop(double Cx, double Cy, double R, int Count)
{
	Loop2 L;
	for (int k = 0; k < Count; ++k)
	{
		const double A = 2.0 * kPi * k / Count;
		L.P.push_back(Pt2{Cx + R * std::cos(A), Cy + R * std::sin(A)});
		L.N.push_back(Pt2{std::cos(A), std::sin(A)});
	}
	return L;
}

// The + of a Cherry MX stem, counter-clockwise, arms ArmLen long and ArmWidth thick.
std::vector<Pt2> PlusOutline(double Cx, double Cy, double ArmLen, double ArmWidth)
{
	const double A = 0.5 * ArmLen, B = 0.5 * ArmWidth;
	const double Xy[12][2] = {{A, -B}, {A, B}, {B, B}, {B, A}, {-B, A}, {-B, B}, {-A, B}, {-A, -B}, {-B, -B}, {-B, -A}, {B, -A}, {B, -B}};
	std::vector<Pt2> Out;
	for (const auto& V : Xy)
	{
		Out.push_back(Pt2{Cx + V[0], Cy + V[1]});
	}
	return Out;
}

Loop2 FromPoints(const std::vector<Pt2>& Pts)
{
	Loop2 L;
	L.P = Pts;
	L.N.assign(Pts.size(), Pt2{0.0, 0.0});
	return L;
}

template <typename ZFn>
std::vector<Vec3> Lift(const Loop2& L, ZFn Z)
{
	std::vector<Vec3> Out;
	Out.reserve(L.P.size());
	for (const Pt2& Q : L.P)
	{
		Out.push_back(Vec3(Q.X, Q.Y, Z(Q.X, Q.Y)));
	}
	return Out;
}

std::vector<Vec3> AtZ(const Loop2& L, double Z)
{
	return Lift(L, [Z](double, double) { return Z; });
}

Vec2 PlanarUv(const Vec3& Q, const Vec3& N)
{
	const double S = 1.0 / 20.0;
	if (std::fabs(N.Z) >= std::fabs(N.X) && std::fabs(N.Z) >= std::fabs(N.Y))
	{
		return Vec2{Q.X * S, Q.Y * S};
	}
	if (std::fabs(N.X) >= std::fabs(N.Y))
	{
		return Vec2{Q.Y * S, Q.Z * S};
	}
	return Vec2{Q.X * S, Q.Z * S};
}

// A flat-shaded triangle wound so that it faces N.
void FlatTri(Mesh& M, const Vec3& A, Vec3 B, Vec3 C, const Vec3& N)
{
	if ((B - A).Cross(C - A).Dot(N) < 0.0)
	{
		std::swap(B, C);
	}
	const Vec3 Nn = N.Normalized();
	const int Ia = M.AddVertex(A, Nn, PlanarUv(A, Nn));
	const int Ib = M.AddVertex(B, Nn, PlanarUv(B, Nn));
	const int Ic = M.AddVertex(C, Nn, PlanarUv(C, Nn));
	M.AddTriangle(Ia, Ib, Ic);
}

// Quads between consecutive rings (each a closed loop, counter-clockwise about the direction from
// ring k to ring k + 1 when the surface should face outward). Normals are smooth around and along the
// strip, except at loop indices marked Sharp, which get one normal per adjacent face. Quads with
// SkipFrom <= i < SkipTo are left out for the caller to fill.
void Loft(Mesh& M, const std::vector<std::vector<Vec3>>& Rings, const std::vector<char>* Sharp = nullptr, int SkipFrom = -1, int SkipTo = -1)
{
	const int K = static_cast<int>(Rings.size());
	if (K < 2)
	{
		return;
	}
	const int N = static_cast<int>(Rings[0].size());
	std::vector<int> IdL(static_cast<size_t>(K * N)), IdR(static_cast<size_t>(K * N));
	for (int k = 0; k < K; ++k)
	{
		const std::vector<Vec3>& Ring = Rings[static_cast<size_t>(k)];
		for (int i = 0; i < N; ++i)
		{
			const Vec3& C = Ring[static_cast<size_t>(i)];
			const Vec3 Up = Rings[static_cast<size_t>(std::min(k + 1, K - 1))][static_cast<size_t>(i)] - Rings[static_cast<size_t>(std::max(k - 1, 0))][static_cast<size_t>(i)];
			const Vec3& Prev = Ring[static_cast<size_t>((i + N - 1) % N)];
			const Vec3& Next = Ring[static_cast<size_t>((i + 1) % N)];
			const Vec2 UV{static_cast<double>(i) / N, static_cast<double>(k) / (K - 1)};
			const size_t Slot = static_cast<size_t>(k * N + i);
			if (Sharp && (*Sharp)[static_cast<size_t>(i)])
			{
				IdL[Slot] = M.AddVertex(C, (C - Prev).Cross(Up), UV);
				IdR[Slot] = M.AddVertex(C, (Next - C).Cross(Up), UV);
			}
			else
			{
				IdL[Slot] = IdR[Slot] = M.AddVertex(C, (Next - Prev).Cross(Up), UV);
			}
		}
	}
	for (int k = 0; k + 1 < K; ++k)
	{
		for (int i = 0; i < N; ++i)
		{
			if (i >= SkipFrom && i < SkipTo)
			{
				continue;
			}
			const int j = (i + 1) % N;
			M.AddQuad(IdR[static_cast<size_t>(k * N + i)], IdL[static_cast<size_t>(k * N + j)], IdL[static_cast<size_t>((k + 1) * N + j)], IdR[static_cast<size_t>((k + 1) * N + i)]);
		}
	}
}

std::vector<char> AllSharp(size_t N)
{
	return std::vector<char>(N, 1);
}

// A flat fan from Center to a convex loop, facing N.
void FanCap(Mesh& M, const std::vector<Vec3>& Ring, const Vec3& Center, const Vec3& N)
{
	const size_t Count = Ring.size();
	for (size_t i = 0; i < Count; ++i)
	{
		FlatTri(M, Center, Ring[i], Ring[(i + 1) % Count], N);
	}
}

// Same-count loops joined by a flat strip facing N.
void StripFlat(Mesh& M, const std::vector<Vec3>& A, const std::vector<Vec3>& B, const Vec3& N)
{
	const size_t Count = A.size();
	for (size_t i = 0; i < Count; ++i)
	{
		const size_t j = (i + 1) % Count;
		FlatTri(M, A[i], A[j], B[j], N);
		FlatTri(M, A[i], B[j], B[i], N);
	}
}

void PlaneBasis(const Vec3& N, Vec3& E1, Vec3& E2)
{
	const Vec3 Ref = std::fabs(N.Z) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
	E1 = Ref.Cross(N).Normalized();
	E2 = N.Cross(E1);
}

// The flat region between two loops in one plane (Inner inside Outer, both star-shaped about Center),
// facing N. The loops are merged by angle about Center, so they may have any vertex counts.
void Bridge(Mesh& M, std::vector<Vec3> Outer, std::vector<Vec3> Inner, const Vec3& Center, const Vec3& N)
{
	Vec3 E1, E2;
	PlaneBasis(N.Normalized(), E1, E2);
	auto Prepare = [&Center, &E1, &E2](std::vector<Vec3>& Ring, std::vector<double>& Ang) {
		const size_t Count = Ring.size();
		double Area = 0.0;
		for (size_t i = 0; i < Count; ++i)
		{
			const Vec3 A = Ring[i] - Center, B = Ring[(i + 1) % Count] - Center;
			Area += A.Dot(E1) * B.Dot(E2) - A.Dot(E2) * B.Dot(E1);
		}
		if (Area < 0.0)
		{
			std::reverse(Ring.begin(), Ring.end());
		}
		Ang.resize(Count);
		size_t Start = 0;
		for (size_t i = 0; i < Count; ++i)
		{
			const Vec3 D = Ring[i] - Center;
			Ang[i] = std::atan2(D.Dot(E2), D.Dot(E1));
			if (Ang[i] < Ang[Start])
			{
				Start = i;
			}
		}
		std::rotate(Ring.begin(), Ring.begin() + static_cast<long>(Start), Ring.end());
		std::rotate(Ang.begin(), Ang.begin() + static_cast<long>(Start), Ang.end());
		for (size_t i = 1; i < Count; ++i)
		{
			if (Ang[i] < Ang[i - 1] - kPi)
			{
				Ang[i] += 2.0 * kPi;
			}
		}
		Ang.push_back(Ang[0] + 2.0 * kPi);
	};
	std::vector<double> AngA, AngB;
	Prepare(Outer, AngA);
	Prepare(Inner, AngB);
	const size_t Na = Outer.size(), Nb = Inner.size();
	size_t i = 0, j = 0;
	while (i < Na || j < Nb)
	{
		const bool AdvanceOuter = j >= Nb || (i < Na && AngA[i + 1] <= AngB[j + 1]);
		if (AdvanceOuter)
		{
			FlatTri(M, Outer[i % Na], Outer[(i + 1) % Na], Inner[j % Nb], N);
			++i;
		}
		else
		{
			FlatTri(M, Outer[i % Na], Inner[(j + 1) % Nb], Inner[j % Nb], N);
			++j;
		}
	}
}

double Cross2(const Pt2& A, const Pt2& B, const Pt2& C)
{
	return (B.X - A.X) * (C.Y - A.Y) - (B.Y - A.Y) * (C.X - A.X);
}

double SignedArea2(const std::vector<Pt2>& Poly)
{
	double A = 0.0;
	for (size_t i = 0; i < Poly.size(); ++i)
	{
		const Pt2& P0 = Poly[i];
		const Pt2& P1 = Poly[(i + 1) % Poly.size()];
		A += P0.X * P1.Y - P1.X * P0.Y;
	}
	return 0.5 * A;
}

// Ear clipping for a simple counter-clockwise polygon.
std::vector<int> EarClip(const std::vector<Pt2>& Poly)
{
	std::vector<int> Idx(Poly.size());
	for (size_t i = 0; i < Poly.size(); ++i)
	{
		Idx[i] = static_cast<int>(i);
	}
	std::vector<int> Out;
	int Pass = 0;
	while (Idx.size() > 3)
	{
		const size_t Count = Idx.size();
		// A vertex on the ear's boundary blocks it too (else the remainder can collapse to a sliver);
		// only if nothing can be clipped that way does a strict test take over.
		const double Edge = Pass == 0 ? -1e-12 : 1e-12;
		bool Clipped = false;
		for (size_t k = 0; k < Count && !Clipped; ++k)
		{
			const int A = Idx[(k + Count - 1) % Count], B = Idx[k], C = Idx[(k + 1) % Count];
			const Pt2 &Pa = Poly[static_cast<size_t>(A)], &Pb = Poly[static_cast<size_t>(B)], &Pc = Poly[static_cast<size_t>(C)];
			if (Cross2(Pa, Pb, Pc) <= 1e-12)
			{
				continue;
			}
			bool Blocked = false;
			for (int Q : Idx)
			{
				if (Q == A || Q == B || Q == C)
				{
					continue;
				}
				const Pt2& Pq = Poly[static_cast<size_t>(Q)];
				if (Cross2(Pa, Pb, Pq) > Edge && Cross2(Pb, Pc, Pq) > Edge && Cross2(Pc, Pa, Pq) > Edge)
				{
					Blocked = true;
					break;
				}
			}
			if (!Blocked)
			{
				Out.push_back(A);
				Out.push_back(B);
				Out.push_back(C);
				Idx.erase(Idx.begin() + static_cast<long>(k));
				Clipped = true;
			}
		}
		if (!Clipped)
		{
			if (Pass == 1)
			{
				break;
			}
			Pass = 1;
		}
		else
		{
			Pass = 0;
		}
	}
	if (Idx.size() == 3)
	{
		Out.push_back(Idx[0]);
		Out.push_back(Idx[1]);
		Out.push_back(Idx[2]);
	}
	return Out;
}

// A flat-sided prism: Poly (either winding) in the frame's (U, V) plane, from W0 to W1 along W.
void Prism(Mesh& M, std::vector<Pt2> Poly, double W0, double W1, const Frame3& F)
{
	if (SignedArea2(Poly) < 0.0)
	{
		std::reverse(Poly.begin(), Poly.end());
	}
	const double Hand = F.U.Cross(F.V).Dot(F.W) >= 0.0 ? 1.0 : -1.0;
	const double Lo = std::min(W0, W1), Hi = std::max(W0, W1);
	const std::vector<int> Tris = EarClip(Poly);
	for (size_t t = 0; t + 2 < Tris.size(); t += 3)
	{
		const Pt2 &A = Poly[static_cast<size_t>(Tris[t])], &B = Poly[static_cast<size_t>(Tris[t + 1])], &C = Poly[static_cast<size_t>(Tris[t + 2])];
		FlatTri(M, F.At(A.X, A.Y, Hi), F.At(B.X, B.Y, Hi), F.At(C.X, C.Y, Hi), F.W * Hand);
		FlatTri(M, F.At(A.X, A.Y, Lo), F.At(B.X, B.Y, Lo), F.At(C.X, C.Y, Lo), F.W * -Hand);
	}
	for (size_t i = 0; i < Poly.size(); ++i)
	{
		const Pt2& P0 = Poly[i];
		const Pt2& P1 = Poly[(i + 1) % Poly.size()];
		const Vec3 Out = (F.U * (P1.Y - P0.Y) - F.V * (P1.X - P0.X)) * Hand;
		const Vec3 A = F.At(P0.X, P0.Y, Lo), B = F.At(P1.X, P1.Y, Lo), C = F.At(P1.X, P1.Y, Hi), D = F.At(P0.X, P0.Y, Hi);
		FlatTri(M, A, B, C, Out);
		FlatTri(M, A, C, D, Out);
	}
}

std::vector<Pt2> RectPoly(double X0, double Y0, double X1, double Y1)
{
	return {Pt2{X0, Y0}, Pt2{X1, Y0}, Pt2{X1, Y1}, Pt2{X0, Y1}};
}

void Box(Mesh& M, double X0, double Y0, double Z0, double X1, double Y1, double Z1)
{
	Prism(M, RectPoly(X0, Y0, X1, Y1), Z0, Z1, Frame3());
}

Pt2 LoopCentroid(const Loop2& L)
{
	Pt2 C;
	for (const Pt2& Q : L.P)
	{
		C.X += Q.X;
		C.Y += Q.Y;
	}
	C.X /= static_cast<double>(L.P.size());
	C.Y /= static_cast<double>(L.P.size());
	return C;
}

// A vertical solid from a convex loop: smooth walls, flat caps.
void SolidFromLoop(Mesh& M, const Loop2& L, double Z0, double Z1)
{
	const Pt2 C = LoopCentroid(L);
	const std::vector<Vec3> Bot = AtZ(L, Z0), Top = AtZ(L, Z1);
	Loft(M, {Bot, Top});
	FanCap(M, Top, Vec3(C.X, C.Y, Z1), Vec3(0, 0, 1));
	FanCap(M, Bot, Vec3(C.X, C.Y, Z0), Vec3(0, 0, -1));
}

void Cylinder(Mesh& M, double Cx, double Cy, double R, int Count, double Z0, double Z1)
{
	SolidFromLoop(M, CircleLoop(Cx, Cy, R, Count), Z0, Z1);
}

void ApplyFrame(Mesh& M, const Frame3& F)
{
	for (Vec3& Q : M.Positions)
	{
		Q = F.At(Q.X, Q.Y, Q.Z);
	}
	for (Vec3& N : M.Normals)
	{
		N = (F.U * N.X + F.V * N.Y + F.W * N.Z).Normalized();
	}
}

void RotateAboutZ(Mesh& M, double A)
{
	const double C = std::cos(A), S = std::sin(A);
	for (Vec3& Q : M.Positions)
	{
		Q = Vec3(Q.X * C - Q.Y * S, Q.X * S + Q.Y * C, Q.Z);
	}
	for (Vec3& N : M.Normals)
	{
		N = Vec3(N.X * C - N.Y * S, N.X * S + N.Y * C, N.Z);
	}
}

// A round tube along Path (parallel-transport frames), with flat end caps.
void SweepTube(Mesh& M, const std::vector<Vec3>& Path, double R, int Sides)
{
	const size_t Count = Path.size();
	if (Count < 2)
	{
		return;
	}
	std::vector<Vec3> T(Count);
	for (size_t i = 0; i < Count; ++i)
	{
		const Vec3& A = Path[i > 0 ? i - 1 : 0];
		const Vec3& B = Path[i + 1 < Count ? i + 1 : Count - 1];
		T[i] = (B - A).Normalized();
	}
	const Vec3 Ref = std::fabs(T[0].Z) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
	Vec3 N = (Ref - T[0] * Ref.Dot(T[0])).Normalized();
	std::vector<std::vector<Vec3>> Rings;
	for (size_t i = 0; i < Count; ++i)
	{
		N = (N - T[i] * N.Dot(T[i])).Normalized();
		const Vec3 B = T[i].Cross(N);
		std::vector<Vec3> Ring;
		for (int s = 0; s < Sides; ++s)
		{
			const double A = 2.0 * kPi * s / Sides;
			Ring.push_back(Path[i] + (N * std::cos(A) + B * std::sin(A)) * R);
		}
		Rings.push_back(Ring);
	}
	Loft(M, Rings);
	FanCap(M, Rings.front(), Path.front(), -T.front());
	FanCap(M, Rings.back(), Path.back(), T.back());
}

// Rounds the corners of a polyline with quadratic arcs of roughly radius R.
std::vector<Vec3> FilletPath(const std::vector<Vec3>& Pts, double R, int Steps)
{
	std::vector<Vec3> Out;
	Out.push_back(Pts.front());
	for (size_t i = 1; i + 1 < Pts.size(); ++i)
	{
		const Vec3 In = Pts[i] - Pts[i - 1], Away = Pts[i + 1] - Pts[i];
		const double L0 = In.Length(), L1 = Away.Length();
		const Vec3 D0 = In * (1.0 / L0), D1 = Away * (1.0 / L1);
		const double Ang = std::acos(std::max(-1.0, std::min(1.0, D0.Dot(D1))));
		const double Cut = std::min(R * std::tan(0.5 * Ang), 0.45 * std::min(L0, L1));
		const Vec3 S = Pts[i] - D0 * Cut, E = Pts[i] + D1 * Cut;
		for (int k = 0; k <= Steps; ++k)
		{
			const double T = static_cast<double>(k) / Steps;
			Out.push_back(S * ((1 - T) * (1 - T)) + Pts[i] * (2 * T * (1 - T)) + E * (T * T));
		}
	}
	Out.push_back(Pts.back());
	return Out;
}

// A polyline in a plane thickened to a closed polygon of width 2 * Half.
std::vector<Pt2> ThickenPolyline(const std::vector<Pt2>& Pts, double Half)
{
	std::vector<Pt2> Left, Right;
	for (size_t i = 0; i < Pts.size(); ++i)
	{
		double Nx = 0.0, Ny = 0.0;
		for (int Side = -1; Side <= 0; ++Side)
		{
			const long A = static_cast<long>(i) + Side, B = A + 1;
			if (A < 0 || B >= static_cast<long>(Pts.size()))
			{
				continue;
			}
			const double Dx = Pts[static_cast<size_t>(B)].X - Pts[static_cast<size_t>(A)].X, Dy = Pts[static_cast<size_t>(B)].Y - Pts[static_cast<size_t>(A)].Y;
			const double L = std::sqrt(Dx * Dx + Dy * Dy);
			Nx += -Dy / L;
			Ny += Dx / L;
		}
		const double L = std::sqrt(Nx * Nx + Ny * Ny);
		Nx /= L;
		Ny /= L;
		Left.push_back(Pt2{Pts[i].X + Nx * Half, Pts[i].Y + Ny * Half});
		Right.push_back(Pt2{Pts[i].X - Nx * Half, Pts[i].Y - Ny * Half});
	}
	std::vector<Pt2> Poly = Left;
	Poly.insert(Poly.end(), Right.rbegin(), Right.rend());
	return Poly;
}

// --- Keycaps ---------------------------------------------------------------------------------------

struct CapRowProfile
{
	double BackMm;
	double TiltDeg;
};

CapRowProfile CapRow(int Row)
{
	switch (std::max(1, std::min(4, Row)))
	{
	case 1:
		return {9.4, 7.0};
	case 2:
		return {8.0, 4.5};
	case 3:
		return {7.4, 2.0};
	default:
		return {7.0, -4.5};
	}
}

struct CapShape
{
	double Hw0 = 0, Hd0 = 0, R0 = 0;
	double HwT = 0, HdT = 0, RT = 0, YcT = 0, YF = 0, YB = 0;
	double ZBack = 0, Tilt = 0, TanT = 0, Skirt = 0;
	double Dish = 0, DishR = 0, Crown = 0;
	bool Spacebar = false;

	double PlaneZ(double Y) const { return ZBack - (YB - Y) * TanT; }
	double Off(double X, double Y) const
	{
		if (Spacebar)
		{
			const double T = std::max(-1.0, std::min(1.0, (Y - YcT) / HdT));
			return Crown * (1.0 - T * T);
		}
		if (Dish <= 0.0)
		{
			return 0.0;
		}
		const double Xc = std::max(-HwT, std::min(HwT, X));
		return -(std::sqrt(DishR * DishR - Xc * Xc) - (DishR - Dish));
	}
	double TopZ(double X, double Y) const { return PlaneZ(Y) + Off(X, Y); }
	Vec3 TopNormal(double X, double Y) const
	{
		const double H = 1e-4;
		const double Dx = (TopZ(X + H, Y) - TopZ(X - H, Y)) / (2 * H);
		const double Dy = (TopZ(X, Y + H) - TopZ(X, Y - H)) / (2 * H);
		return Vec3(-Dx, -Dy, 1.0).Normalized();
	}
	double DeepestOff() const { return Spacebar ? 0.0 : Dish; }
};

CapShape MakeCapShape(const KeycapParams& P)
{
	CapShape S;
	const CapRowProfile Row = CapRow(P.ProfileRow);
	S.Hw0 = 0.5 * (P.WidthU * P.UnitMm - P.GapMm);
	S.Hd0 = 0.5 * (P.DepthU * P.UnitMm - P.GapMm);
	S.R0 = P.CornerRadiusMm;
	S.HwT = S.Hw0 - P.TopInsetMm;
	S.YF = -S.Hd0 + P.TopInsetMm * kCapFrontInsetScale;
	S.YB = S.Hd0 - P.TopInsetMm * kCapBackInsetScale;
	S.HdT = 0.5 * (S.YB - S.YF);
	S.YcT = 0.5 * (S.YB + S.YF);
	S.RT = std::min(P.TopCornerRadiusMm, 0.9 * std::min(S.HwT, S.HdT));
	S.Skirt = -kCapStemInset;
	S.Tilt = Row.TiltDeg * kDeg;
	S.TanT = std::tan(S.Tilt);
	S.ZBack = Row.BackMm + S.Skirt;
	S.Spacebar = P.Spacebar;
	S.Crown = P.Spacebar ? kSpacebarCrown : 0.0;
	S.Dish = P.Spacebar ? 0.0 : std::max(0.0, P.DishDepthMm);
	if (S.Dish > 0.0)
	{
		S.DishR = (S.HwT * S.HwT + S.Dish * S.Dish) / (2.0 * S.Dish);
	}
	return S;
}

void AddHomingBar(Mesh& M, const CapShape& S)
{
	const double By = S.YF + 2.0;
	const Splits None;
	auto Base = [&S](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() - 0.3; };
	auto Shoulder = [&S](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() + 0.3; };
	auto Crest = [&S](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() + 0.45; };
	const Loop2 Outer = RoundedRect(0.0, By, 2.6, 0.5, 0.45, 5, None);
	const Loop2 Inner = RoundedRect(0.0, By, 2.4, 0.3, 0.28, 5, None);
	const std::vector<Vec3> R0 = Lift(Outer, Base), R1 = Lift(Outer, Shoulder), R2 = Lift(Inner, Crest);
	Loft(M, {R0, R1, R2});
	const Vec3 Up = Vec3(0.0, -S.TanT, 1.0).Normalized();
	FanCap(M, R2, Vec3(0.0, By, Crest(0.0, By)), Up);
	FanCap(M, R0, Vec3(0.0, By, Base(0.0, By)), -Up);
}

// --- Plate ----------------------------------------------------------------------------------------

struct HoleRect
{
	double X0, X1, Y0, Y1;
};

void SortUnique(std::vector<double>& V, double Tol)
{
	std::sort(V.begin(), V.end());
	V.erase(std::unique(V.begin(), V.end(), [Tol](double A, double B) { return std::fabs(A - B) <= Tol; }), V.end());
}

// A slab over [X0, X1] x [Y0, Y1] with axis-aligned rectangular through-holes, cut into horizontal
// bands at every hole edge; each band's solid intervals are zipped against the breakpoints of the
// bands on either side, so no face has a T-junction and the result is watertight.
void PerforatedSlab(Mesh& M, double X0, double X1, double Y0, double Y1, const std::vector<HoleRect>& Holes, double Z0, double Z1)
{
	using Interval = std::pair<double, double>;
	const double Eps = 1e-7;
	std::vector<double> Ys = {Y0, Y1};
	for (const HoleRect& H : Holes)
	{
		Ys.push_back(std::max(Y0, std::min(Y1, H.Y0)));
		Ys.push_back(std::max(Y0, std::min(Y1, H.Y1)));
	}
	SortUnique(Ys, Eps);
	const int Nb = static_cast<int>(Ys.size()) - 1;
	std::vector<std::vector<Interval>> Solid(static_cast<size_t>(Nb));
	for (int b = 0; b < Nb; ++b)
	{
		const double Ym = 0.5 * (Ys[static_cast<size_t>(b)] + Ys[static_cast<size_t>(b + 1)]);
		std::vector<Interval> Iv = {{X0, X1}};
		for (const HoleRect& H : Holes)
		{
			if (!(H.Y0 < Ym && Ym < H.Y1))
			{
				continue;
			}
			std::vector<Interval> Next;
			for (const Interval& I : Iv)
			{
				if (H.X1 <= I.first + Eps || H.X0 >= I.second - Eps)
				{
					Next.push_back(I);
					continue;
				}
				if (H.X0 > I.first + Eps)
				{
					Next.push_back({I.first, H.X0});
				}
				if (H.X1 < I.second - Eps)
				{
					Next.push_back({H.X1, I.second});
				}
			}
			Iv = Next;
		}
		std::sort(Iv.begin(), Iv.end());
		Solid[static_cast<size_t>(b)] = Iv;
	}
	std::vector<std::vector<double>> Breaks(static_cast<size_t>(Nb + 1));
	for (int l = 0; l <= Nb; ++l)
	{
		std::vector<double>& Xs = Breaks[static_cast<size_t>(l)];
		for (int b = l - 1; b <= l; ++b)
		{
			if (b < 0 || b >= Nb)
			{
				continue;
			}
			for (const Interval& I : Solid[static_cast<size_t>(b)])
			{
				Xs.push_back(I.first);
				Xs.push_back(I.second);
			}
		}
		SortUnique(Xs, Eps);
	}
	auto Within = [Eps](const std::vector<double>& Xs, double A, double C) {
		std::vector<double> Out;
		for (double X : Xs)
		{
			if (X >= A - Eps && X <= C + Eps)
			{
				Out.push_back(X);
			}
		}
		return Out;
	};
	auto Covered = [](const std::vector<Interval>& Iv, double X) {
		for (const Interval& I : Iv)
		{
			if (I.first < X && X < I.second)
			{
				return true;
			}
		}
		return false;
	};
	auto Wall = [&M, Z0, Z1](double Xa, double Ya, double Xb, double Yb, const Vec3& N) {
		const Vec3 A(Xa, Ya, Z0), B(Xb, Yb, Z0), C(Xb, Yb, Z1), D(Xa, Ya, Z1);
		FlatTri(M, A, B, C, N);
		FlatTri(M, A, C, D, N);
	};
	for (int b = 0; b < Nb; ++b)
	{
		const double Yb = Ys[static_cast<size_t>(b)], Yt = Ys[static_cast<size_t>(b + 1)];
		for (const Interval& I : Solid[static_cast<size_t>(b)])
		{
			const std::vector<double> Bot = Within(Breaks[static_cast<size_t>(b)], I.first, I.second);
			const std::vector<double> Top = Within(Breaks[static_cast<size_t>(b + 1)], I.first, I.second);
			for (int Face = 0; Face < 2; ++Face)
			{
				const double Z = Face == 0 ? Z1 : Z0;
				const Vec3 N(0, 0, Face == 0 ? 1.0 : -1.0);
				size_t i = 0, j = 0;
				while (i + 1 < Bot.size() || j + 1 < Top.size())
				{
					const bool AdvanceBot = j + 1 >= Top.size() || (i + 1 < Bot.size() && Bot[i + 1] <= Top[j + 1]);
					if (AdvanceBot)
					{
						FlatTri(M, Vec3(Bot[i], Yb, Z), Vec3(Bot[i + 1], Yb, Z), Vec3(Top[j], Yt, Z), N);
						++i;
					}
					else
					{
						FlatTri(M, Vec3(Bot[i], Yb, Z), Vec3(Top[j + 1], Yt, Z), Vec3(Top[j], Yt, Z), N);
						++j;
					}
				}
			}
			Wall(I.first, Yb, I.first, Yt, Vec3(-1, 0, 0));
			Wall(I.second, Yb, I.second, Yt, Vec3(1, 0, 0));
		}
	}
	for (int l = 0; l <= Nb; ++l)
	{
		const std::vector<double>& Xs = Breaks[static_cast<size_t>(l)];
		const double Y = Ys[static_cast<size_t>(l)];
		for (size_t k = 0; k + 1 < Xs.size(); ++k)
		{
			const double Xm = 0.5 * (Xs[k] + Xs[k + 1]);
			const bool Below = l > 0 && Covered(Solid[static_cast<size_t>(l - 1)], Xm);
			const bool Above = l < Nb && Covered(Solid[static_cast<size_t>(l)], Xm);
			if (Above && !Below)
			{
				Wall(Xs[k], Y, Xs[k + 1], Y, Vec3(0, -1, 0));
			}
			else if (Below && !Above)
			{
				Wall(Xs[k], Y, Xs[k + 1], Y, Vec3(0, 1, 0));
			}
		}
	}
}

// Stabilizer axis and spacing for a key, or Spacing 0 when it needs none.
struct StabPlacement
{
	double Spacing = 0.0;
	bool Vertical = false;
};

StabPlacement StabFor(double WidthU, double DepthU)
{
	StabPlacement S;
	if (WidthU >= 2.0 - 1e-9)
	{
		S.Spacing = StabilizerSpacingMm(WidthU);
	}
	else if (DepthU >= 2.0 - 1e-9)
	{
		S.Spacing = StabilizerSpacingMm(DepthU);
		S.Vertical = true;
	}
	return S;
}
} // namespace keyboardparts_detail

double KeycapHeightMm(int ProfileRow)
{
	return keyboardparts_detail::CapRow(ProfileRow).BackMm;
}

KeycapTop KeycapTopFor(const KeycapParams& P)
{
	using namespace keyboardparts_detail;
	const CapShape S = MakeCapShape(P);
	KeycapTop T;
	T.CenterY = S.YcT;
	T.CenterZ = S.TopZ(0.0, S.YcT);
	T.TiltRad = S.Tilt;
	T.Width = 2.0 * S.HwT;
	T.Depth = 2.0 * S.HdT / std::cos(S.Tilt);
	const double TopFrontZ = S.TopZ(0.0, S.YF);
	T.FrontFaceY = 0.5 * (-S.Hd0 + S.YF);
	T.FrontFaceZ = 0.5 * (S.Skirt + TopFrontZ);
	T.FrontTiltRad = std::atan2(S.YF + S.Hd0, TopFrontZ - S.Skirt);
	return T;
}

Mesh BuildKeycap(const KeycapParams& P)
{
	using namespace keyboardparts_detail;
	Mesh M;
	const CapShape S = MakeCapShape(P);
	const double W = P.WallMm;
	const int Seg = std::max(2, P.Segments);

	// Stem sockets: the switch at the centre, plus the stabilizer posts on long keys.
	const StabPlacement Stab = StabFor(P.WidthU, P.DepthU);
	std::vector<double> SocketAt = {0.0};
	if (Stab.Spacing > 0.0)
	{
		SocketAt = {-Stab.Spacing, 0.0, Stab.Spacing};
	}

	// The inner cavity: the walls offset by W, under a ceiling parallel to the top.
	auto CeilingZ = [&S, W](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() - W; };
	const double Tm = (CeilingZ(0.0, S.YcT) - S.Skirt) / (S.PlaneZ(S.YcT) - S.Skirt);
	const double HwI = Lerp(S.Hw0, S.HwT, Tm) - W;
	const double YFI = Lerp(-S.Hd0, S.YF, Tm) + W;
	const double YBI = Lerp(S.Hd0, S.YB, Tm) - W;
	const double HdI = 0.5 * (YBI - YFI), YcI = 0.5 * (YBI + YFI);
	const double RI = std::max(0.4, Lerp(S.R0, S.RT, Tm) - W);
	const double ExI = HwI - RI, EyI = HdI - RI;

	// Edge points shared by every loop: uniform ones for the dish, exact ones on the ceiling cuts
	// between sockets, and dense ones around the sockets so the ceiling bridges stay well shaped.
	const double XStep = S.Spacebar ? 1e9 : 1.3 * std::max(1.0, P.WidthU / 1.5);
	const double YStep = S.Spacebar ? 1.0 : 2.0 * std::max(1.0, P.DepthU / 1.5);
	const int XCount = std::max(S.Spacebar ? 4 : 3, std::min(60, static_cast<int>(2.0 * (S.HwT - S.RT) / XStep)));
	const int YCount = std::max(2, std::min(30, static_cast<int>(2.0 * (S.HdT - S.RT) / YStep)));
	std::vector<double> XBase = Uniform(XCount), YBase = Uniform(YCount);
	std::vector<double> XExact, YExact;
	const double Dense[] = {0.0, 1.5, 3.0, 4.5, 6.5, 9.0, 13.0, 18.0, 25.0};
	const bool Sparse = (Stab.Vertical ? 2.0 * EyI / (YCount + 1) : 2.0 * ExI / (XCount + 1)) > 3.0;
	for (size_t k = 0; k + 1 < SocketAt.size(); ++k)
	{
		const double Cut = 0.5 * (SocketAt[k] + SocketAt[k + 1]);
		if (Stab.Vertical)
		{
			YExact.push_back((Cut - YcI) / EyI);
		}
		else
		{
			XExact.push_back(Cut / ExI);
		}
	}
	if (SocketAt.size() > 1 && Sparse)
	{
		for (double C : SocketAt)
		{
			for (double D : Dense)
			{
				for (double Sgn = -1.0; Sgn <= 1.0; Sgn += 2.0)
				{
					if (Stab.Vertical)
					{
						YBase.push_back((C + Sgn * D - YcI) / EyI);
					}
					else
					{
						XBase.push_back((C + Sgn * D) / ExI);
					}
				}
			}
		}
	}
	const Splits Sp = SameSplits(MergeFractions(XBase, XExact, 0.4 / ExI), MergeFractions(YBase, YExact, 0.4 / EyI));

	const Loop2 Base = RoundedRect(0.0, 0.0, S.Hw0, S.Hd0, S.R0, Seg, Sp);
	const Loop2 Top = RoundedRect(0.0, S.YcT, S.HwT, S.HdT, S.RT, Seg, Sp);
	const Loop2 InBase = RoundedRect(0.0, 0.0, S.Hw0 - W, S.Hd0 - W, std::max(0.4, S.R0 - W), Seg, Sp);
	const Loop2 InTop = RoundedRect(0.0, YcI, HwI, HdI, RI, Seg, Sp);
	const size_t N = Base.P.size();

	// Outer walls, lofted from the footprint to the dished top outline with a slight outward bulge.
	const std::vector<Vec3> TopRing = Lift(Top, [&S](double X, double Y) { return S.TopZ(X, Y); });
	const std::vector<Vec3> BaseRing = AtZ(Base, S.Skirt);
	const int WallSteps = 4;
	std::vector<std::vector<Vec3>> Walls;
	for (int k = 0; k <= WallSteps; ++k)
	{
		if (k == 0 || k == WallSteps)
		{
			Walls.push_back(k == 0 ? BaseRing : TopRing);
			continue;
		}
		const double T = static_cast<double>(k) / WallSteps;
		const double Bulge = kCapWallBulge * std::sin(kPi * T);
		std::vector<Vec3> Ring;
		for (size_t i = 0; i < N; ++i)
		{
			Ring.push_back(Vec3(Lerp(Base.P[i].X, Top.P[i].X, T) + Base.N[i].X * Bulge, Lerp(Base.P[i].Y, Top.P[i].Y, T) + Base.N[i].Y * Bulge, Lerp(S.Skirt, TopRing[i].Z, T)));
		}
		Walls.push_back(Ring);
	}
	Loft(M, Walls);

	// The top: concentric copies of the outline shrunk toward its centre, on the dish surface. Its UVs
	// span 0..1 over the top; every other vertex of the cap gets U < 0 so a legend can be masked.
	const int TopFirst = static_cast<int>(M.Positions.size());
	const int Rings = S.Spacebar ? 4 : 7;
	std::vector<std::vector<int>> Grid(Rings);
	auto TopVertex = [&M, &S](double X, double Y) {
		const Vec2 UV{(X + S.HwT) / (2.0 * S.HwT), (Y - S.YF) / (2.0 * S.HdT)};
		return M.AddVertex(Vec3(X, Y, S.TopZ(X, Y)), S.TopNormal(X, Y), UV);
	};
	for (int k = 0; k < Rings; ++k)
	{
		const double Scale = 1.0 - static_cast<double>(k) / Rings;
		for (size_t i = 0; i < N; ++i)
		{
			const double X = k == 0 ? Top.P[i].X : Scale * Top.P[i].X;
			const double Y = k == 0 ? Top.P[i].Y : S.YcT + Scale * (Top.P[i].Y - S.YcT);
			Grid[static_cast<size_t>(k)].push_back(TopVertex(X, Y));
		}
	}
	const int Center = TopVertex(0.0, S.YcT);
	const int TopEnd = static_cast<int>(M.Positions.size());
	for (int k = 0; k < Rings; ++k)
	{
		for (size_t i = 0; i < N; ++i)
		{
			const size_t j = (i + 1) % N;
			const std::vector<int>& G = Grid[static_cast<size_t>(k)];
			if (k + 1 < Rings)
			{
				const std::vector<int>& H = Grid[static_cast<size_t>(k + 1)];
				M.AddQuad(G[i], G[j], H[j], H[i]);
			}
			else
			{
				M.AddTriangle(G[i], G[j], Center);
			}
		}
	}

	// Underside: the skirt's bottom rim and the inner walls.
	const std::vector<Vec3> InBaseRing = AtZ(InBase, S.Skirt);
	const std::vector<Vec3> InTopRing = Lift(InTop, CeilingZ);
	StripFlat(M, BaseRing, InBaseRing, Vec3(0, 0, -1));
	Loft(M, {InTopRing, InBaseRing});

	// Ceiling, split between the sockets, and the sockets with their + cut-outs.
	const Vec3 CeilingN = Vec3(0.0, S.TanT, -1.0).Normalized();
	const int CircleN = 20;
	for (size_t r = 0; r < SocketAt.size(); ++r)
	{
		const double Sx = Stab.Vertical ? 0.0 : SocketAt[r];
		const double Sy = Stab.Vertical ? SocketAt[r] : 0.0;
		const double Lo = r == 0 ? -1e30 : 0.5 * (SocketAt[r - 1] + SocketAt[r]);
		const double Hi = r + 1 == SocketAt.size() ? 1e30 : 0.5 * (SocketAt[r] + SocketAt[r + 1]);
		std::vector<size_t> Picked;
		for (size_t i = 0; i < N; ++i)
		{
			const double C = Stab.Vertical ? InTop.P[i].Y : InTop.P[i].X;
			if (C >= Lo - 1e-6 && C <= Hi + 1e-6)
			{
				Picked.push_back(i);
			}
		}
		// Where the picked run jumps across a cut, the cut is filled in with points about 1 mm apart,
		// interpolated the same way from both sides so the neighbouring regions share them exactly.
		std::vector<Vec3> Region;
		for (size_t k = 0; k < Picked.size(); ++k)
		{
			const size_t A = Picked[k], Bn = Picked[(k + 1) % Picked.size()];
			Region.push_back(InTopRing[A]);
			if ((A + 1) % N == Bn)
			{
				continue;
			}
			const bool Fwd = InTopRing[A].X + InTopRing[A].Y < InTopRing[Bn].X + InTopRing[Bn].Y;
			const Vec3& From = Fwd ? InTopRing[A] : InTopRing[Bn];
			const Vec3& To = Fwd ? InTopRing[Bn] : InTopRing[A];
			const int Steps = std::max(2, static_cast<int>(std::ceil((To - From).Length() / 1.0)));
			for (int s = 1; s < Steps; ++s)
			{
				const int Q = Fwd ? s : Steps - s;
				Region.push_back(From + (To - From) * (static_cast<double>(Q) / Steps));
			}
		}
		const Loop2 Circle = CircleLoop(Sx, Sy, kCapSocketRadius, CircleN);
		const std::vector<Vec3> SocketTop = Lift(Circle, CeilingZ);
		const std::vector<Vec3> SocketBot = AtZ(Circle, 0.0);
		Bridge(M, Region, SocketTop, Vec3(Sx, Sy, CeilingZ(Sx, Sy)), CeilingN);
		Loft(M, {SocketBot, SocketTop});
		const Loop2 Cross = FromPoints(PlusOutline(Sx, Sy, kMxCrossArm, kMxCrossWidth));
		const std::vector<Vec3> CrossBot = AtZ(Cross, 0.0), CrossTop = AtZ(Cross, kCapSocketDepth);
		Bridge(M, SocketBot, CrossBot, Vec3(Sx, Sy, 0.0), Vec3(0, 0, -1));
		const std::vector<char> Sharp = AllSharp(Cross.P.size());
		Loft(M, {CrossTop, CrossBot}, &Sharp);
		const std::vector<int> Tris = EarClip(Cross.P);
		for (size_t t = 0; t + 2 < Tris.size(); t += 3)
		{
			FlatTri(M, CrossTop[static_cast<size_t>(Tris[t])], CrossTop[static_cast<size_t>(Tris[t + 1])], CrossTop[static_cast<size_t>(Tris[t + 2])], Vec3(0, 0, -1));
		}
	}

	if (P.HomingBar)
	{
		AddHomingBar(M, S);
	}
	for (int i = 0; i < static_cast<int>(M.UVs.size()); ++i)
	{
		if (i < TopFirst || i >= TopEnd)
		{
			Vec2& UV = M.UVs[static_cast<size_t>(i)];
			UV.U = -2.0 + 0.5 * (UV.U - std::floor(UV.U));
		}
	}
	return M;
}

SwitchParts BuildSwitch(const SwitchParams& P)
{
	using namespace keyboardparts_detail;
	SwitchParts Out;
	const double Mount = P.StemTopZ;
	const double HousingTop = Mount - 0.3;
	const double PcbTop = -kMxPlateToPcb;
	const double Floor = PcbTop + 0.8;
	const int Seg = 3;
	const Vec3 Up(0, 0, 1), Down(0, 0, -1);
	const Splits Sq = SameSplits(Uniform(5), Uniform(5));

	// Bottom housing: a cup that passes through the plate cut-out, narrowed above the clip ledge,
	// with the spring boss inside and the 4 mm centre post below.
	{
		Mesh& M = Out.HousingBottom;
		const Loop2 Full = RoundedRect(0, 0, 6.95, 6.95, 0.5, Seg, Sq);
		const Loop2 Narrow = RoundedRect(0, 0, 6.3, 6.95, 0.5, Seg, Sq);
		const Loop2 Pocket = RoundedRect(0, 0, 5.4, 5.9, 0.6, Seg, Sq);
		const double Ledge = -2.5;
		Loft(M, {AtZ(Full, PcbTop), AtZ(Full, Ledge)});
		StripFlat(M, AtZ(Full, Ledge), AtZ(Narrow, Ledge), Up);
		Loft(M, {AtZ(Narrow, Ledge), AtZ(Narrow, 0.0)});
		StripFlat(M, AtZ(Narrow, 0.0), AtZ(Pocket, 0.0), Up);
		Loft(M, {AtZ(Pocket, 0.0), AtZ(Pocket, Floor)});
		const Loop2 Boss = CircleLoop(0, 0, 1.5, 16);
		Bridge(M, AtZ(Pocket, Floor), AtZ(Boss, Floor), Vec3(0, 0, Floor), Up);
		Loft(M, {AtZ(Boss, Floor), AtZ(Boss, Floor + 1.0)});
		FanCap(M, AtZ(Boss, Floor + 1.0), Vec3(0, 0, Floor + 1.0), Up);
		const Loop2 Post = CircleLoop(0, 0, 1.95, 20);
		const Loop2 Tip = CircleLoop(0, 0, 1.45, 20);
		const double PostEnd = PcbTop - 3.0;
		Bridge(M, AtZ(Full, PcbTop), AtZ(Post, PcbTop), Vec3(0, 0, PcbTop), Down);
		Loft(M, {AtZ(Post, PostEnd + 0.5), AtZ(Post, PcbTop)});
		Loft(M, {AtZ(Tip, PostEnd), AtZ(Post, PostEnd + 0.5)});
		FanCap(M, AtZ(Tip, PostEnd), Vec3(0, 0, PostEnd), Down);
	}

	// Top housing: a hollow shell on the plate, tapering to the stem window, plus the two plate clips
	// whose shoulders hook under the 14 mm cut-out.
	{
		Mesh& M = Out.HousingTop;
		const Splits Sq4 = SameSplits(Uniform(4), Uniform(4));
		const Loop2 A = RoundedRect(0, 0, 7.8, 7.8, 0.8, Seg, Sq4);
		const Loop2 C = RoundedRect(0, 0, 7.45, 7.45, 0.8, Seg, Sq4);
		const Loop2 D = RoundedRect(0, 0, 5.7, 5.7, 1.4, Seg, Sq4);
		const Loop2 IA = RoundedRect(0, 0, 7.0, 7.0, 0.4, Seg, Sq4);
		const Loop2 ID = RoundedRect(0, 0, 5.16, 5.16, 0.6, Seg, Sq4);
		const Loop2 Win = RoundedRect(0, 0, 3.7, 2.9, 0.4, Seg, SameSplits(Uniform(3), Uniform(2)));
		const double Ceiling = HousingTop - 0.8;
		Loft(M, {AtZ(A, 0.0), AtZ(A, 1.3)});
		Loft(M, {AtZ(A, 1.3), AtZ(C, 1.9)});
		Loft(M, {AtZ(C, 1.9), AtZ(D, HousingTop)});
		Bridge(M, AtZ(D, HousingTop), AtZ(Win, HousingTop), Vec3(0, 0, HousingTop), Up);
		Loft(M, {AtZ(Win, HousingTop), AtZ(Win, Ceiling)});
		Bridge(M, AtZ(ID, Ceiling), AtZ(Win, Ceiling), Vec3(0, 0, Ceiling), Down);
		Loft(M, {AtZ(ID, Ceiling), AtZ(IA, 2.2)});
		Loft(M, {AtZ(IA, 2.2), AtZ(IA, 0.0)});
		StripFlat(M, AtZ(A, 0.0), AtZ(IA, 0.0), Down);
		const std::vector<Pt2> Clip = {{6.4, -2.4}, {6.95, -2.4}, {7.3, -1.75}, {7.3, -1.55}, {6.95, -1.55}, {6.95, 0.0}, {7.3, 0.0}, {7.3, 1.2}, {6.4, 1.2}};
		Frame3 F;
		F.V = Vec3(0, 0, 1);
		F.W = Vec3(0, -1, 0);
		std::vector<Pt2> Mirror;
		for (const Pt2& Q : Clip)
		{
			Mirror.push_back(Pt2{-Q.X, Q.Y});
		}
		Prism(M, Clip, -2.2, 2.2, F);
		Prism(M, Mirror, -2.2, 2.2, F);
	}

	// Stem, in its own frame (z = 0 is the cap mount): the + cross, the slider (hollow from below for
	// the spring), the guide wings, the legs that ride the contact leaf and the spring seat.
	{
		Mesh& M = Out.Stem;
		const double Arm = kMxCrossArm - 2.0 * kMxStemClearance, Thick = kMxCrossWidth - 2.0 * kMxStemClearance;
		Prism(M, PlusOutline(0, 0, Arm, Thick), -0.2, kMxCrossHeight, Frame3());
		const Loop2 Slider = RoundedRect(0, 0, 3.6, 2.8, 0.5, Seg, SameSplits(Uniform(4), Uniform(3)));
		const Loop2 Bore = CircleLoop(0, 0, 2.5, 20);
		const double SliderBot = -4.6, BoreTop = -1.6;
		Loft(M, {AtZ(Slider, SliderBot), AtZ(Slider, 0.0)});
		FanCap(M, AtZ(Slider, 0.0), Vec3(0, 0, 0), Up);
		Bridge(M, AtZ(Slider, SliderBot), AtZ(Bore, SliderBot), Vec3(0, 0, SliderBot), Down);
		Loft(M, {AtZ(Bore, BoreTop), AtZ(Bore, SliderBot)});
		FanCap(M, AtZ(Bore, BoreTop), Vec3(0, 0, BoreTop), Down);
		SolidFromLoop(M, RoundedRect(0, 0, 4.7, 2.0, 0.3, 2, Splits()), -2.4, -1.3);
		Cylinder(M, 0, 0, 1.2, 14, -4.2, BoreTop + 0.1);
		std::vector<Pt2> Leg = {{2.3, -6.7}, {3.3, -6.7}, {3.3, -4.4}, {2.3, -4.4}};
		if (P.Type != SwitchType::Linear)
		{
			Leg = {{2.3, -6.7}, {3.3, -6.7}, {3.3, -6.3}, {3.85, -5.75}, {3.3, -5.2}, {3.3, -4.4}, {2.3, -4.4}};
		}
		Frame3 F;
		F.U = Vec3(0, 1, 0);
		F.V = Vec3(0, 0, 1);
		F.W = Vec3(1, 0, 0);
		Prism(M, Leg, 1.25, 2.15, F);
		Prism(M, Leg, -2.15, -1.25, F);
		if (P.Type == SwitchType::Clicky)
		{
			SolidFromLoop(M, RoundedRect(0, -0.6, 4.1, 2.6, 0.6, Seg, Splits()), -5.5, -4.5);
		}
	}

	// Spring: a helix of SpringTurns turns with closed (touching) end coils, seated on the housing
	// floor and reaching up into the slider's bore.
	{
		const double Rc = 2.25, Rw = 0.15;
		Out.SpringBottomZ = Floor;
		Out.SpringRestTopZ = Mount - 1.6;
		const int Turns = std::max(4, P.SpringTurns);
		const int PerTurn = 12;
		const double Z0 = Out.SpringBottomZ + Rw, Z1 = Out.SpringRestTopZ - Rw;
		const double EndRise = 2.05 * Rw;
		const double MidRise = std::max(0.0, (Z1 - Z0) - 2.0 * EndRise);
		std::vector<Vec3> Path;
		for (int s = 0; s <= Turns * PerTurn; ++s)
		{
			const double U = static_cast<double>(s) / PerTurn;
			double Z = Z0;
			if (U <= 1.0)
			{
				Z = Z0 + EndRise * U;
			}
			else if (U >= Turns - 1.0)
			{
				Z = Z1 - EndRise * (Turns - U);
			}
			else
			{
				Z = Z0 + EndRise + MidRise * (U - 1.0) / (Turns - 2.0);
			}
			Path.push_back(Vec3(Rc * std::cos(2.0 * kPi * U), Rc * std::sin(2.0 * kPi * U), Z));
		}
		SweepTube(Out.Spring, Path, Rw, 6);
	}

	// Contact leaves: a straight fixed leaf and the bent moving leaf the stem legs push away.
	{
		Frame3 F;
		F.U = Vec3(0, 1, 0);
		F.V = Vec3(0, 0, 1);
		F.W = Vec3(1, 0, 0);
		Prism(Out.Contacts, RectPoly(5.0, Floor, 5.2, -0.9), -3.4, -1.0, F);
		const std::vector<Pt2> Bent = ThickenPolyline({{5.1, Floor}, {5.1, -2.9}, {4.35, -1.5}, {4.55, -0.9}}, 0.1);
		Prism(Out.Contacts, Bent, 0.8, 3.2, F);
	}

	// Pins: flat-ended copper pins 3.3 mm below the housing.
	{
		const double Pins[2][2] = {{kMxPin1X, kMxPin1Y}, {kMxPin2X, kMxPin2Y}};
		for (const auto& Q : Pins)
		{
			const Loop2 Sec = FromPoints(RectPoly(Q[0] - 0.35, Q[1] - 0.35, Q[0] + 0.35, Q[1] + 0.35));
			const Loop2 TipSec = FromPoints(RectPoly(Q[0] - 0.15, Q[1] - 0.15, Q[0] + 0.15, Q[1] + 0.15));
			const std::vector<char> Sharp = AllSharp(4);
			const double End = PcbTop - kMxPinBelowHousing;
			Loft(Out.Pins, {AtZ(TipSec, End), AtZ(Sec, End + 0.4)}, &Sharp);
			Loft(Out.Pins, {AtZ(Sec, End + 0.4), AtZ(Sec, Floor - 0.2)}, &Sharp);
			FanCap(Out.Pins, AtZ(TipSec, End), Vec3(Q[0], Q[1], End), Down);
			FanCap(Out.Pins, AtZ(Sec, Floor - 0.2), Vec3(Q[0], Q[1], Floor - 0.2), Up);
		}
	}
	return Out;
}

double StabilizerSpacingMm(double WidthU)
{
	if (WidthU < 2.0 - 1e-9)
	{
		return 0.0;
	}
	if (WidthU < 3.0 - 1e-9)
	{
		return 11.938;
	}
	if (WidthU < 6.0 - 1e-9)
	{
		return 19.05;
	}
	if (WidthU < 6.25 - 1e-9)
	{
		return 47.625;
	}
	if (WidthU < 7.0 - 1e-9)
	{
		return 50.0;
	}
	return 57.15;
}

StabilizerParts BuildStabilizer(double WidthU, bool Vertical)
{
	using namespace keyboardparts_detail;
	StabilizerParts Out;
	const double S = StabilizerSpacingMm(WidthU);
	if (S <= 0.0)
	{
		return Out;
	}
	const Vec3 Up(0, 0, 1), Down(0, 0, -1);
	const double BodyTop = 4.8, WindowFloor = 0.3, WireZ = 5.6, BarY = 8.6, ArmX = S + 4.4, WireR = 0.6;
	for (double Side = -1.0; Side <= 1.0; Side += 2.0)
	{
		const double Cx = Side * S;
		const Loop2 Body = RoundedRect(Cx, 0, 3.6, 5.5, 0.8, 3, SameSplits(Uniform(3), Uniform(4)));
		const Loop2 Win = RoundedRect(Cx, 0, 3.0, 3.0, 0.5, 3, SameSplits(Uniform(3), Uniform(3)));
		Loft(Out.Housings, {AtZ(Body, 0.0), AtZ(Body, BodyTop)});
		FanCap(Out.Housings, AtZ(Body, 0.0), Vec3(Cx, 0, 0.0), Down);
		Bridge(Out.Housings, AtZ(Body, BodyTop), AtZ(Win, BodyTop), Vec3(Cx, 0, BodyTop), Up);
		Loft(Out.Housings, {AtZ(Win, BodyTop), AtZ(Win, WindowFloor)});
		FanCap(Out.Housings, AtZ(Win, WindowFloor), Vec3(Cx, 0, WindowFloor), Up);
		SolidFromLoop(Out.Housings, RoundedRect(Cx, 0, 0.5 * kStabSlotW - 0.05, 0.5 * kStabSlotH - 0.1, 0.3, 2, Splits()), -2.2, 0.05);

		SolidFromLoop(Out.Stems, RoundedRect(Cx, 0, 2.8, 2.8, 0.4, 3, Splits()), -2.2, 0.0);
		Prism(Out.Stems, PlusOutline(Cx, 0, kMxCrossArm - 2.0 * kMxStemClearance, kMxCrossWidth - 2.0 * kMxStemClearance), -0.2, kMxCrossHeight, Frame3());
	}
	const double EndZ = WireZ;
	const std::vector<Vec3> Corners = {Vec3(-S + 1.2, 0, EndZ), Vec3(-ArmX, 0, EndZ), Vec3(-ArmX, BarY, WireR + 0.05), Vec3(ArmX, BarY, WireR + 0.05), Vec3(ArmX, 0, EndZ), Vec3(S - 1.2, 0, EndZ)};
	SweepTube(Out.Wire, FilletPath(Corners, 1.5, 6), WireR, 8);
	if (Vertical)
	{
		RotateAboutZ(Out.Housings, 0.5 * kPi);
		RotateAboutZ(Out.Stems, 0.5 * kPi);
		RotateAboutZ(Out.Wire, 0.5 * kPi);
	}
	return Out;
}

Vec3 KeyCenter(const mkb::KeyDef& Key, const mkb::Layout& Layout)
{
	return Vec3((Key.CenterX() - 0.5 * Layout.Width) * mkb::kUnitMm, (0.5 * Layout.Height - Key.CenterY()) * mkb::kUnitMm, 0.0);
}

BoardParts BuildBoard(const mkb::Layout& Layout, const CaseParams& P)
{
	using namespace keyboardparts_detail;
	BoardParts B;
	const double U = mkb::kUnitMm;
	const double Hw = 0.5 * Layout.Width * U, Hd = 0.5 * Layout.Height * U;
	const double PcbTop = -(P.PlateThickness + P.PcbGapBelowPlate);
	const double PcbBot = PcbTop - P.PcbThickness;
	const double Angle = P.TypingAngleDeg * kDeg;
	const double TanA = std::tan(Angle);
	const double M = P.MarginMm;
	const double YFront = -Hd - M, YBack = Hd + M;
	auto BottomZ = [&P, YFront, TanA](double, double Y) { return -P.DepthBelowPlate - (Y - YFront) * TanA; };
	const double EIn = 3.5, EPlate = 2.8, EPcb = 2.4;
	const double FloorZ = -P.DepthBelowPlate + P.WallMm;
	const Vec3 Up(0, 0, 1), Down(0, 0, -1);

	// The USB-C receptacle sits in the gap between F8 and F9, at the back edge of the PCB.
	double UsbX = 0.0;
	{
		const mkb::KeyDef* F8 = Layout.Find("F8");
		const mkb::KeyDef* F9 = Layout.Find("F9");
		if (F8 && F9)
		{
			UsbX = 0.5 * (KeyCenter(*F8, Layout).X + KeyCenter(*F9, Layout).X);
		}
	}
	const double UsbH = 3.26, UsbW = 8.94, UsbD = 7.35;
	const double UsbZ = PcbTop + 0.5 * UsbH;

	// Case: one closed surface from the wedge bottom up the outer wall, over the chamfer and the rim,
	// and down the pocket to its floor, with the USB-C port cut through the back wall.
	{
		Mesh& C = B.CaseShell;
		const int Seg = 10;
		const double R = P.CornerRadiusMm, Ch = P.ChamferMm, Lip = P.LipAbovePlate;
		const double ExOuter = Hw + M - R;
		std::vector<double> Back;
		const double Offsets[] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 12, 14, 17, 20, 24, 29, 35, 42, 50, 60, 72, 86, 102, 120, 140, 165, 195};
		for (double D : Offsets)
		{
			Back.push_back((UsbX + D) / ExOuter);
			Back.push_back((UsbX - D) / ExOuter);
		}
		Splits Sp;
		Sp.Top = MergeFractions(Back, {}, 0.2 / ExOuter);
		const int BackFrom = Seg, BackTo = Seg + static_cast<int>(Sp.Top.size()) + 1;
		auto Ring = [&](double E, double Rad, auto Z) { return Lift(RoundedRect(0, 0, Hw + E, Hd + E, Rad, Seg, Sp), Z); };
		auto Flat = [](double Z) { return [Z](double, double) { return Z; }; };
		const std::vector<Vec3> Ob = Ring(M - 1.0, R - 0.6, BottomZ);
		const std::vector<Vec3> O0 = Ring(M, R, [&BottomZ](double X, double Y) { return BottomZ(X, Y) + 1.0; });
		const std::vector<Vec3> O1 = Ring(M, R, Flat(Lip - Ch));
		const std::vector<Vec3> O2 = Ring(M - Ch, R - 0.6 * Ch, Flat(Lip));
		const std::vector<Vec3> I0 = Ring(EIn + 0.8, 2.6, Flat(Lip));
		const std::vector<Vec3> I1 = Ring(EIn, 2.0, Flat(Lip - 0.8));
		const std::vector<Vec3> I2 = Ring(EIn, 2.0, Flat(FloorZ));
		Loft(C, {Ob, O0});
		Loft(C, {O0, O1}, nullptr, BackFrom, BackTo);
		Loft(C, {O1, O2});
		Loft(C, {O2, I0});
		Loft(C, {I0, I1});
		Loft(C, {I1, I2}, nullptr, BackFrom, BackTo);
		FanCap(C, Ob, Vec3(0, 0, BottomZ(0, 0)), Vec3(0, -TanA, -1));
		FanCap(C, I2, Vec3(0, 0, FloorZ), Up);
		const Loop2 Port = RoundedRect(UsbX, UsbZ, 6.0, 3.0, 2.6, 6, SameSplits(Uniform(3), {}));
		auto PortRing = [&Port](double Y) {
			std::vector<Vec3> Out;
			for (const Pt2& Q : Port.P)
			{
				Out.push_back(Vec3(Q.X, Y, Q.Y));
			}
			return Out;
		};
		auto Face = [BackFrom, BackTo](const std::vector<Vec3>& Lo, const std::vector<Vec3>& Hi) {
			std::vector<Vec3> Poly;
			for (int i = BackFrom; i <= BackTo; ++i)
			{
				Poly.push_back(Lo[static_cast<size_t>(i)]);
			}
			for (int i = BackTo; i >= BackFrom; --i)
			{
				Poly.push_back(Hi[static_cast<size_t>(i)]);
			}
			return Poly;
		};
		const double YIn = Hd + EIn;
		Bridge(C, Face(O0, O1), PortRing(YBack), Vec3(UsbX, YBack, UsbZ), Vec3(0, 1, 0));
		Bridge(C, Face(I2, I1), PortRing(YIn), Vec3(UsbX, YIn, UsbZ), Vec3(0, -1, 0));
		Loft(C, {PortRing(YIn), PortRing(YBack)});
	}

	// Badge: a small chamfered plate set into the back wall.
	{
		Mesh Badge;
		const double Bx = Hw - 40.0, Bz = 0.5 * (BottomZ(0, YBack) + 1.0 + P.LipAbovePlate - P.ChamferMm);
		const Loop2 Plate = RoundedRect(0, 0, 18.0, 4.5, 1.2, 5, SameSplits(Uniform(4), {}));
		const Loop2 Face = RoundedRect(0, 0, 17.6, 4.1, 0.9, 5, SameSplits(Uniform(4), {}));
		Loft(Badge, {AtZ(Plate, -0.8), AtZ(Plate, 0.4), AtZ(Face, 0.8)});
		FanCap(Badge, AtZ(Face, 0.8), Vec3(0, 0, 0.8), Up);
		FanCap(Badge, AtZ(Plate, -0.8), Vec3(0, 0, -0.8), Down);
		Frame3 F;
		F.O = Vec3(Bx, YBack, Bz);
		F.U = Vec3(1, 0, 0);
		F.V = Vec3(0, 0, -1);
		F.W = Vec3(0, 1, 0);
		ApplyFrame(Badge, F);
		B.CaseBadge = Badge;
	}

	// Feet: four rubber pucks square to the wedge bottom.
	{
		const Vec3 Out = Vec3(0, -std::sin(Angle), -std::cos(Angle));
		Frame3 F;
		F.U = Vec3(1, 0, 0);
		F.W = Out;
		F.V = F.W.Cross(F.U);
		const double Fx = Hw + M - 28.0;
		const double Fy[2] = {YFront + 18.0, YBack - 18.0};
		for (double Sx = -1.0; Sx <= 1.0; Sx += 2.0)
		{
			for (double Y : Fy)
			{
				Mesh Foot;
				const Loop2 Rim = CircleLoop(0, 0, 7.0, 28), Sole = CircleLoop(0, 0, 6.3, 28);
				Loft(Foot, {AtZ(Rim, -1.0), AtZ(Rim, 1.9), AtZ(Sole, 2.5)});
				FanCap(Foot, AtZ(Sole, 2.5), Vec3(0, 0, 2.5), Up);
				FanCap(Foot, AtZ(Rim, -1.0), Vec3(0, 0, -1.0), Down);
				F.O = Vec3(Sx * Fx, Y, BottomZ(0, Y));
				ApplyFrame(Foot, F);
				B.Feet.Append(Foot);
			}
		}
		B.RestTiltRad = Angle;
		double Lowest = 1e30;
		for (const Vec3& Q : B.Feet.Positions)
		{
			Lowest = std::min(Lowest, Q.Y * std::sin(Angle) + Q.Z * std::cos(Angle));
		}
		B.RestLiftMm = -Lowest;
	}

	// Plate: 14 x 14 mm cut-outs at every switch and 3.3 x 14 mm slots for the stabilizers.
	{
		std::vector<HoleRect> Holes;
		const double H = 0.5 * kMxCutout;
		for (const mkb::KeyDef& K : Layout.Keys)
		{
			const Vec3 C = KeyCenter(K, Layout);
			Holes.push_back({C.X - H, C.X + H, C.Y - H, C.Y + H});
			const StabPlacement Stab = StabFor(K.W, K.H);
			if (Stab.Spacing <= 0.0)
			{
				continue;
			}
			for (double Side = -1.0; Side <= 1.0; Side += 2.0)
			{
				if (Stab.Vertical)
				{
					const double Sy = C.Y + Side * Stab.Spacing;
					Holes.push_back({C.X - 0.5 * kStabSlotH, C.X + 0.5 * kStabSlotH, Sy - 0.5 * kStabSlotW, Sy + 0.5 * kStabSlotW});
				}
				else
				{
					const double Sx = C.X + Side * Stab.Spacing;
					Holes.push_back({Sx - 0.5 * kStabSlotW, Sx + 0.5 * kStabSlotW, C.Y - 0.5 * kStabSlotH, C.Y + 0.5 * kStabSlotH});
				}
			}
		}
		PerforatedSlab(B.Plate, -Hw - EPlate, Hw + EPlate, -Hd - EPlate, Hd + EPlate, Holes, -P.PlateThickness, 0.0);
	}

	// PCB with its pin pads, and the parts on it.
	{
		SolidFromLoop(B.Pcb, RoundedRect(0, 0, Hw + EPcb, Hd + EPcb, 1.5, 6, Splits()), PcbBot, PcbTop);
		for (const mkb::KeyDef& K : Layout.Keys)
		{
			const Vec3 C = KeyCenter(K, Layout);
			Cylinder(B.Pcb, C.X + kMxPin1X, C.Y + kMxPin1Y, 0.95, 10, PcbTop - 0.02, PcbTop + 0.05);
			Cylinder(B.Pcb, C.X + kMxPin2X, C.Y + kMxPin2Y, 0.95, 10, PcbTop - 0.02, PcbTop + 0.05);
			// SOD-123 diode below the switch: 2.7 x 1.6 x 1.1 mm body and two leads.
			const double Dx = C.X, Dy = C.Y - 9.6;
			Box(B.PcbParts, Dx - 1.35, Dy - 0.8, PcbTop, Dx + 1.35, Dy + 0.8, PcbTop + 1.1);
			Box(B.PcbParts, Dx - 1.85, Dy - 0.45, PcbTop - 0.01, Dx - 1.2, Dy + 0.45, PcbTop + 0.2);
			Box(B.PcbParts, Dx + 1.2, Dy - 0.45, PcbTop - 0.01, Dx + 1.85, Dy + 0.45, PcbTop + 0.2);
		}
		// Controller: a 7 x 7 mm QFN-44 between the function and number rows, with a crystal and
		// decoupling capacitors.
		double QfnY = Hd - 2.0 * U;
		{
			const mkb::KeyDef* F5 = Layout.Find("F5");
			const mkb::KeyDef* D5 = Layout.Find("D5");
			if (F5 && D5)
			{
				QfnY = 0.5 * (KeyCenter(*F5, Layout).Y + KeyCenter(*D5, Layout).Y) - 2.0;
			}
		}
		const double QfnX = UsbX - 24.0;
		Box(B.PcbParts, QfnX - 3.5, QfnY - 3.5, PcbTop, QfnX + 3.5, QfnY + 3.5, PcbTop + 0.85);
		for (int k = 0; k < 11; ++k)
		{
			const double T = -2.5 + 0.5 * k;
			Box(B.PcbParts, QfnX + T - 0.12, QfnY - 3.6, PcbTop - 0.01, QfnX + T + 0.12, QfnY - 3.0, PcbTop + 0.2);
			Box(B.PcbParts, QfnX + T - 0.12, QfnY + 3.0, PcbTop - 0.01, QfnX + T + 0.12, QfnY + 3.6, PcbTop + 0.2);
			Box(B.PcbParts, QfnX - 3.6, QfnY + T - 0.12, PcbTop - 0.01, QfnX - 3.0, QfnY + T + 0.12, PcbTop + 0.2);
			Box(B.PcbParts, QfnX + 3.0, QfnY + T - 0.12, PcbTop - 0.01, QfnX + 3.6, QfnY + T + 0.12, PcbTop + 0.2);
		}
		Box(B.PcbParts, QfnX + 6.0, QfnY - 1.25, PcbTop, QfnX + 9.2, QfnY + 1.25, PcbTop + 0.8);
		const double Caps[4][2] = {{-5.5, 2.0}, {-5.5, -2.0}, {2.0, -5.5}, {-2.0, 5.5}};
		for (const auto& Q : Caps)
		{
			const bool Along = std::fabs(Q[0]) < std::fabs(Q[1]);
			const double Hx = Along ? 0.8 : 0.4, Hy = Along ? 0.4 : 0.8;
			Box(B.PcbParts, QfnX + Q[0] - Hx, QfnY + Q[1] - Hy, PcbTop, QfnX + Q[0] + Hx, QfnY + Q[1] + Hy, PcbTop + 0.45);
		}
		// USB-C receptacle (8.94 x 7.35 x 3.26 mm shell) with its mouth facing out of the back edge.
		{
			Mesh Usb;
			const Loop2 Shell = RoundedRect(0, 0, 0.5 * UsbW, 0.5 * UsbH, 1.5, 6, SameSplits(Uniform(4), {}));
			const Loop2 Mouth = RoundedRect(0, 0, 0.5 * UsbW - 0.3, 0.5 * UsbH - 0.35, 1.15, 6, SameSplits(Uniform(4), {}));
			Loft(Usb, {AtZ(Shell, 0.0), AtZ(Shell, UsbD)});
			FanCap(Usb, AtZ(Shell, 0.0), Vec3(0, 0, 0), Down);
			StripFlat(Usb, AtZ(Shell, UsbD), AtZ(Mouth, UsbD), Up);
			Loft(Usb, {AtZ(Mouth, UsbD), AtZ(Mouth, 1.0)});
			FanCap(Usb, AtZ(Mouth, 1.0), Vec3(0, 0, 1.0), Up);
			Box(Usb, -3.3, -0.35, 0.9, 3.3, 0.35, UsbD - 0.7);
			Frame3 F;
			F.O = Vec3(UsbX, Hd + EPcb + 0.6 - UsbD, UsbZ);
			F.U = Vec3(1, 0, 0);
			F.V = Vec3(0, 0, -1);
			F.W = Vec3(0, 1, 0);
			ApplyFrame(Usb, F);
			B.PcbParts.Append(Usb);
		}
	}

	// OLED module: a bezel on the plate with a chamfered, tilted top and a pocket holding the glass.
	{
		const mkb::Rect& D = Layout.Display;
		const double Dcx = (D.X + 0.5 * D.W - 0.5 * Layout.Width) * U;
		const double Dcy = (0.5 * Layout.Height - (D.Y + 0.5 * D.H)) * U;
		const double BHw = 0.5 * (D.W * U - 0.8), BHd = 0.5 * (D.H * U - 0.8);
		const double Tilt = 8.0 * kDeg, TanT = std::tan(Tilt), Zc = 9.0, Bevel = 0.5;
		const double GHw = 33.0, GHd = 6.6, Clear = 0.2, PocketDepth = 1.2, GlassDrop = 0.15;
		auto TopZ = [Dcy, TanT, Zc](double, double Y) { return Zc + (Y - Dcy) * TanT; };
		auto Below = [&TopZ](double Dz) { return [&TopZ, Dz](double X, double Y) { return TopZ(X, Y) - Dz; }; };
		const Splits BezelSp = SameSplits(Uniform(30), Uniform(5));
		const Loop2 Outer = RoundedRect(Dcx, Dcy, BHw, BHd, 2.0, 6, BezelSp);
		const Loop2 Inset = RoundedRect(Dcx, Dcy, BHw - Bevel, BHd - Bevel, 1.5, 6, BezelSp);
		const Loop2 Pocket = RoundedRect(Dcx, Dcy, GHw + Clear, GHd + Clear, 0.4, 4, SameSplits(Uniform(30), Uniform(5)));
		const Vec3 TopN = Vec3(0, -TanT, 1).Normalized();
		Mesh& Bz = B.DisplayBezel;
		Loft(Bz, {AtZ(Outer, 0.0), Lift(Outer, Below(Bevel))});
		Loft(Bz, {Lift(Outer, Below(Bevel)), Lift(Inset, Below(0.0))});
		FanCap(Bz, AtZ(Outer, 0.0), Vec3(Dcx, Dcy, 0.0), Down);
		Bridge(Bz, Lift(Inset, Below(0.0)), Lift(Pocket, Below(0.0)), Vec3(Dcx, Dcy, TopZ(0, Dcy)), TopN);
		Loft(Bz, {Lift(Pocket, Below(0.0)), Lift(Pocket, Below(PocketDepth))});
		FanCap(Bz, Lift(Pocket, Below(PocketDepth)), Vec3(Dcx, Dcy, TopZ(0, Dcy) - PocketDepth), TopN);

		Mesh& G = B.DisplayGlass;
		const Loop2 Glass = RoundedRect(Dcx, Dcy, GHw, GHd, 0.2, 2, SameSplits(Uniform(30), Uniform(5)));
		const std::vector<Vec3> GTop = Lift(Glass, Below(GlassDrop)), GBot = Lift(Glass, Below(PocketDepth - 0.02));
		Loft(G, {GBot, GTop});
		FanCap(G, GBot, Vec3(Dcx, Dcy, TopZ(0, Dcy) - PocketDepth + 0.02), -TopN);
		auto GlassVertex = [&G, &TopN, Dcx, Dcy, GHw, GHd](const Vec3& Q) {
			const Vec2 UV{(Q.X - (Dcx - GHw)) / (2.0 * GHw), ((Dcy + GHd) - Q.Y) / (2.0 * GHd)};
			return G.AddVertex(Q, TopN, UV);
		};
		const Vec3 GC(Dcx, Dcy, TopZ(0, Dcy) - GlassDrop);
		const int Gc = GlassVertex(GC);
		std::vector<int> Rim;
		for (const Vec3& Q : GTop)
		{
			Rim.push_back(GlassVertex(Q));
		}
		for (size_t i = 0; i < Rim.size(); ++i)
		{
			G.AddTriangle(Gc, Rim[i], Rim[(i + 1) % Rim.size()]);
		}
		B.DisplayCenter = GC;
		B.DisplayWidthMm = 2.0 * GHw;
		B.DisplayHeightMm = 2.0 * GHd / std::cos(Tilt);
		B.DisplayTiltRad = Tilt;
	}
	return B;
}
} // namespace mgeo
