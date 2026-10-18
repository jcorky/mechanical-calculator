// Shared mesh-building primitives (see MeshPrimitives.h).
#include "MeshPrimitives.h"

#include <algorithm>
#include <cmath>

namespace mgeo
{
namespace prim
{
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
	// A (near-)degenerate triangle keeps the caller's winding: flipping it on rounding noise would
	// break the edge pairing with its neighbours.
	const double Scale = (B - A).Length() * (C - A).Length() * N.Length();
	if ((B - A).Cross(C - A).Dot(N) < -1e-9 * Scale)
	{
		std::swap(B, C);
	}
	const Vec3 Nn = N.Normalized();
	const int Ia = M.AddVertex(A, Nn, PlanarUv(A, Nn));
	const int Ib = M.AddVertex(B, Nn, PlanarUv(B, Nn));
	const int Ic = M.AddVertex(C, Nn, PlanarUv(C, Nn));
	M.AddTriangle(Ia, Ib, Ic);
}

void FlatTriFacet(Mesh& M, const Vec3& A, Vec3 B, Vec3 C, const Vec3& Hint)
{
	Vec3 F = (B - A).Cross(C - A);
	if (F.Dot(Hint) < -1e-9 * (B - A).Length() * (C - A).Length() * Hint.Length())
	{
		std::swap(B, C);
		F = -F;
	}
	const Vec3 Nn = F.Normalized();
	const int Ia = M.AddVertex(A, Nn, PlanarUv(A, Nn));
	const int Ib = M.AddVertex(B, Nn, PlanarUv(B, Nn));
	const int Ic = M.AddVertex(C, Nn, PlanarUv(C, Nn));
	M.AddTriangle(Ia, Ib, Ic);
}

// Quads between consecutive rings (each a closed loop, counter-clockwise about the direction from
// ring k to ring k + 1 when the surface should face outward). Normals are smooth around and along the
// strip, except at loop indices marked Sharp, which get one normal per adjacent face. Quads with
// SkipFrom <= i < SkipTo are left out for the caller to fill.
void Loft(Mesh& M, const std::vector<std::vector<Vec3>>& Rings, const std::vector<char>* Sharp, int SkipFrom, int SkipTo)
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
// The loop's winding about N (positive: counter-clockwise seen from the side N points to).
double WindingAbout(const std::vector<Vec3>& Ring, const Vec3& Center, const Vec3& N)
{
	double W = 0.0;
	for (size_t i = 0; i < Ring.size(); ++i)
	{
		W += (Ring[i] - Center).Cross(Ring[(i + 1) % Ring.size()] - Center).Dot(N);
	}
	return W;
}

void FanCap(Mesh& M, const std::vector<Vec3>& Ring, const Vec3& Center, const Vec3& N)
{
	// Orient the whole fan once, so collinear slivers (radial edges through the centre) keep the
	// winding of their neighbours.
	const size_t Count = Ring.size();
	const bool Ccw = WindingAbout(Ring, Center, N) >= 0.0;
	for (size_t i = 0; i < Count; ++i)
	{
		const Vec3& A = Ring[i];
		const Vec3& B = Ring[(i + 1) % Count];
		FlatTri(M, Center, Ccw ? A : B, Ccw ? B : A, N);
	}
}

// Same-count loops joined by a flat strip facing N.
void StripFlat(Mesh& M, const std::vector<Vec3>& A, const std::vector<Vec3>& B, const Vec3& N)
{
	const size_t Count = A.size();
	Vec3 C;
	for (const Vec3& Q : A)
	{
		C = C + Q;
	}
	C = C * (1.0 / static_cast<double>(Count));
	// With A outside B and both counter-clockwise about N, (A_i, A_j, B_j) faces N.
	const bool Ccw = WindingAbout(A, C, N) >= 0.0;
	for (size_t i = 0; i < Count; ++i)
	{
		const size_t j = (i + 1) % Count;
		if (Ccw)
		{
			FlatTri(M, A[i], A[j], B[j], N);
			FlatTri(M, A[i], B[j], B[i], N);
		}
		else
		{
			FlatTri(M, A[i], B[j], A[j], N);
			FlatTri(M, A[i], B[i], B[j], N);
		}
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
void Bridge(Mesh& M, std::vector<Vec3> Outer, std::vector<Vec3> Inner, const Vec3& Center, const Vec3& N, bool Facets)
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
	const Vec3 Nn = N.Normalized();
	auto Orient = [&Nn](const Vec3& P, const Vec3& Q, const Vec3& R) { return (Q - P).Cross(R - P).Dot(Nn); };
	size_t i = 0, j = 0;
	while (i < Na || j < Nb)
	{
		bool AdvanceOuter = j >= Nb || (i < Na && AngA[i + 1] <= AngB[j + 1]);
		if (i < Na && j < Nb)
		{
			// The angular order can pick a back-facing triangle where the outer loop steps radially
			// (a gear's flank); take the other candidate then.
			const double Ao = Orient(Outer[i % Na], Outer[(i + 1) % Na], Inner[j % Nb]);
			const double Ai = Orient(Outer[i % Na], Inner[(j + 1) % Nb], Inner[j % Nb]);
			const double Eps = 1e-12;
			if (AdvanceOuter && Ao < -Eps && Ai > Eps)
			{
				AdvanceOuter = false;
			}
			else if (!AdvanceOuter && Ai < -Eps && Ao > Eps)
			{
				AdvanceOuter = true;
			}
		}
		const Vec3& A = Outer[i % Na];
		const Vec3& B = AdvanceOuter ? Outer[(i + 1) % Na] : Inner[(j + 1) % Nb];
		const Vec3& C = Inner[j % Nb];
		if (Facets)
		{
			FlatTriFacet(M, A, B, C, N);
		}
		else
		{
			FlatTri(M, A, B, C, N);
		}
		if (AdvanceOuter)
		{
			++i;
		}
		else
		{
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
	// The smallest angle of a triangle, to prefer well-shaped ears (a zigzag rather than a fan of
	// slivers along runs of collinear points).
	auto Quality = [](const Pt2& A, const Pt2& B, const Pt2& C) {
		auto Ang = [](const Pt2& P, const Pt2& Q, const Pt2& R) {
			const double Ux = Q.X - P.X, Uy = Q.Y - P.Y, Vx = R.X - P.X, Vy = R.Y - P.Y;
			return std::atan2(std::fabs(Ux * Vy - Uy * Vx), Ux * Vx + Uy * Vy);
		};
		return std::min({Ang(A, B, C), Ang(B, C, A), Ang(C, A, B)});
	};
	std::vector<int> Out;
	while (Idx.size() > 3)
	{
		const size_t Count = Idx.size();
		size_t Best = Count;
		double BestQ = -1.0;
		for (int Pass = 0; Pass < 2 && Best == Count; ++Pass)
		{
			// A vertex on the ear's boundary blocks it too (else the remainder can collapse to a
			// sliver); only if nothing can be clipped that way does a strict test take over.
			const double Edge = Pass == 0 ? -1e-12 : 1e-12;
			for (size_t k = 0; k < Count; ++k)
			{
				const int A = Idx[(k + Count - 1) % Count], B = Idx[k], C = Idx[(k + 1) % Count];
				const Pt2 &Pa = Poly[static_cast<size_t>(A)], &Pb = Poly[static_cast<size_t>(B)], &Pc = Poly[static_cast<size_t>(C)];
				if (Cross2(Pa, Pb, Pc) <= 1e-12)
				{
					continue;
				}
				const double Q = Quality(Pa, Pb, Pc);
				if (Q <= BestQ)
				{
					continue;
				}
				bool Blocked = false;
				for (int Qi : Idx)
				{
					if (Qi == A || Qi == B || Qi == C)
					{
						continue;
					}
					const Pt2& Pq = Poly[static_cast<size_t>(Qi)];
					auto Same = [&Pq](const Pt2& V) { return std::fabs(Pq.X - V.X) < 1e-12 && std::fabs(Pq.Y - V.Y) < 1e-12; };
					if (Same(Pa) || Same(Pb) || Same(Pc))
					{
						continue;   // a keyhole bridge's duplicate vertex
					}
					if (Cross2(Pa, Pb, Pq) > Edge && Cross2(Pb, Pc, Pq) > Edge && Cross2(Pc, Pa, Pq) > Edge)
					{
						Blocked = true;
						break;
					}
				}
				if (!Blocked)
				{
					Best = k;
					BestQ = Q;
				}
			}
		}
		if (Best == Count)
		{
			break;
		}
		Out.push_back(Idx[(Best + Count - 1) % Count]);
		Out.push_back(Idx[Best]);
		Out.push_back(Idx[(Best + 1) % Count]);
		Idx.erase(Idx.begin() + static_cast<long>(Best));
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
} // namespace prim
} // namespace mgeo
