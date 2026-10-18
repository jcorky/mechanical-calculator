// Mesh-building primitives shared by the keyboard and gear geometry: rounded-rectangle and circle
// loops, lofts with smooth or crisp normals, flat caps and bridges between loops, prisms from 2D
// outlines, swept tubes and slabs with rectangular holes. Same frame and winding as MeshBuilder.h.
#pragma once

#include "MeshBuilder.h"

#include <vector>

namespace mgeo
{
namespace prim
{
constexpr double kPi = 3.14159265358979323846;

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

// A placement: local (A, B, C) maps to O + U * A + V * B + W * C.
struct Frame3
{
	Vec3 O, U = Vec3(1, 0, 0), V = Vec3(0, 1, 0), W = Vec3(0, 0, 1);
	Vec3 At(double A, double B, double C) const { return O + U * A + V * B + W * C; }
};

struct HoleRect
{
	double X0, X1, Y0, Y1;
};

double Lerp(double A, double B, double T);
std::vector<double> Uniform(int Count);
// Sorted fractions: the Exact ones are kept as given and any other closer than MinGap to one is dropped.
std::vector<double> MergeFractions(const std::vector<double>& Base, const std::vector<double>& Exact, double MinGap);
Splits SameSplits(const std::vector<double>& XEdges, const std::vector<double>& YEdges);

// Counter-clockwise rounded rectangle starting at angle 0 on the top-right corner.
Loop2 RoundedRect(double Cx, double Cy, double Hw, double Hd, double R, int Seg, const Splits& S);
Loop2 CircleLoop(double Cx, double Cy, double R, int Count);
// The + of a Cherry MX stem, counter-clockwise, arms ArmLen long and ArmWidth thick.
std::vector<Pt2> PlusOutline(double Cx, double Cy, double ArmLen, double ArmWidth);
Loop2 FromPoints(const std::vector<Pt2>& Pts);

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

std::vector<Vec3> AtZ(const Loop2& L, double Z);
Vec2 PlanarUv(const Vec3& Q, const Vec3& N);

// A flat-shaded triangle wound so that it faces N.
void FlatTri(Mesh& M, const Vec3& A, Vec3 B, Vec3 C, const Vec3& N);
// The same, shaded with the triangle's own normal (wound to face the Hint side).
void FlatTriFacet(Mesh& M, const Vec3& A, Vec3 B, Vec3 C, const Vec3& Hint);

// Quads between consecutive rings (each a closed loop, counter-clockwise about the direction from
// ring k to ring k + 1 when the surface should face outward). Normals are smooth around and along the
// strip, except at loop indices marked Sharp, which get one normal per adjacent face. Quads with
// SkipFrom <= i < SkipTo are left out for the caller to fill.
void Loft(Mesh& M, const std::vector<std::vector<Vec3>>& Rings, const std::vector<char>* Sharp = nullptr, int SkipFrom = -1, int SkipTo = -1);
std::vector<char> AllSharp(size_t N);

double WindingAbout(const std::vector<Vec3>& Ring, const Vec3& Center, const Vec3& N);
// A flat fan from Center to a convex (or star-shaped about Center) loop, facing N.
void FanCap(Mesh& M, const std::vector<Vec3>& Ring, const Vec3& Center, const Vec3& N);
// Same-count loops joined by a flat strip facing N.
void StripFlat(Mesh& M, const std::vector<Vec3>& A, const std::vector<Vec3>& B, const Vec3& N);
void PlaneBasis(const Vec3& N, Vec3& E1, Vec3& E2);
// The flat region between two loops in one plane (Inner inside Outer, both star-shaped about Center),
// facing N. The loops are merged by angle about Center, so they may have any vertex counts.
// With Facets the loops need not be coplanar: each triangle gets its own normal, wound to face N.
void Bridge(Mesh& M, std::vector<Vec3> Outer, std::vector<Vec3> Inner, const Vec3& Center, const Vec3& N, bool Facets = false);

double Cross2(const Pt2& A, const Pt2& B, const Pt2& C);
double SignedArea2(const std::vector<Pt2>& Poly);
// Ear clipping for a simple counter-clockwise polygon.
std::vector<int> EarClip(const std::vector<Pt2>& Poly);
// A flat-sided prism: Poly (either winding) in the frame's (U, V) plane, from W0 to W1 along W.
void Prism(Mesh& M, std::vector<Pt2> Poly, double W0, double W1, const Frame3& F);
std::vector<Pt2> RectPoly(double X0, double Y0, double X1, double Y1);
void Box(Mesh& M, double X0, double Y0, double Z0, double X1, double Y1, double Z1);
Pt2 LoopCentroid(const Loop2& L);
// A vertical solid from a convex loop: smooth walls, flat caps.
void SolidFromLoop(Mesh& M, const Loop2& L, double Z0, double Z1);
void Cylinder(Mesh& M, double Cx, double Cy, double R, int Count, double Z0, double Z1);
void ApplyFrame(Mesh& M, const Frame3& F);
void RotateAboutZ(Mesh& M, double A);
// A round tube along Path (parallel-transport frames), with flat end caps.
void SweepTube(Mesh& M, const std::vector<Vec3>& Path, double R, int Sides);
// Rounds the corners of a polyline with quadratic arcs of roughly radius R.
std::vector<Vec3> FilletPath(const std::vector<Vec3>& Pts, double R, int Steps);
// A polyline in a plane thickened to a closed polygon of width 2 * Half.
std::vector<Pt2> ThickenPolyline(const std::vector<Pt2>& Pts, double Half);
void SortUnique(std::vector<double>& V, double Tol);
// A slab over [X0, X1] x [Y0, Y1] with axis-aligned rectangular through-holes, cut into horizontal
// bands at every hole edge; each band's solid intervals are zipped against the breakpoints of the
// bands on either side, so no face has a T-junction and the result is watertight.
void PerforatedSlab(Mesh& M, double X0, double X1, double Y0, double Y1, const std::vector<HoleRect>& Holes, double Z0, double Z1);
} // namespace prim
} // namespace mgeo
