// Procedural geometry for the keyboard, engine-independent.
//
// Frame: millimetres, right-handed, x to the right along the board, y away from the typist, z up.
// Triangles wind counter-clockwise seen from outside, so (B - A) x (C - A) points outward.
// The Unreal adapter maps (x, y, z) -> (UE.Y, UE.X, UE.Z) * 0.1: swapping two axes converts to
// Unreal's left-handed frame and its front-face convention without reordering any indices.
#pragma once

#include <cmath>
#include <string>
#include <vector>

namespace mgeo
{
struct Vec3
{
	double X = 0.0, Y = 0.0, Z = 0.0;
	Vec3() = default;
	Vec3(double InX, double InY, double InZ) : X(InX), Y(InY), Z(InZ) {}
	Vec3 operator+(const Vec3& O) const { return {X + O.X, Y + O.Y, Z + O.Z}; }
	Vec3 operator-(const Vec3& O) const { return {X - O.X, Y - O.Y, Z - O.Z}; }
	Vec3 operator*(double S) const { return {X * S, Y * S, Z * S}; }
	Vec3 operator-() const { return {-X, -Y, -Z}; }
	double Dot(const Vec3& O) const { return X * O.X + Y * O.Y + Z * O.Z; }
	Vec3 Cross(const Vec3& O) const { return {Y * O.Z - Z * O.Y, Z * O.X - X * O.Z, X * O.Y - Y * O.X}; }
	double Length() const { return std::sqrt(Dot(*this)); }
	Vec3 Normalized() const { const double L = Length(); return L > 0 ? *this * (1.0 / L) : Vec3(0, 0, 1); }
};

struct Vec2
{
	double U = 0.0, V = 0.0;
};

struct Mesh
{
	std::vector<Vec3> Positions;
	std::vector<Vec3> Normals;     // one per position, unit length
	std::vector<Vec2> UVs;         // one per position
	std::vector<int> Triangles;    // 3 indices per triangle

	int AddVertex(const Vec3& P, const Vec3& N, const Vec2& UV);
	void AddTriangle(int A, int B, int C);
	void AddQuad(int A, int B, int C, int D);   // A B C D counter-clockwise from outside
	void Append(const Mesh& Other, const Vec3& Offset = Vec3());
	void Translate(const Vec3& Offset);
	int TriangleCount() const { return static_cast<int>(Triangles.size() / 3); }

	// Checks the tests use. Welding merges positions closer than Tol (mm), so a mesh built from
	// separate flat-shaded faces still counts as one closed surface.
	double SignedVolume() const;                         // divergence theorem; > 0 when outward
	bool IsClosed(double Tol = 1e-6) const;              // every welded edge used exactly twice, oppositely
	int EulerCharacteristic(double Tol = 1e-6) const;    // V - E + F after welding
	bool NormalsAgreeWithWinding(double MinDot = 0.0) const;
};

// Writes a Wavefront OBJ (for inspecting the geometry outside Unreal).
bool WriteObj(const std::string& Path, const std::vector<std::pair<std::string, const Mesh*>>& Parts);
} // namespace mgeo
