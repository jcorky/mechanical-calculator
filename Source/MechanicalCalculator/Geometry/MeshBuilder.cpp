#include "MeshBuilder.h"

#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>

namespace mgeo
{
namespace meshbuilder_detail
{
// Maps every position to a welded id: positions closer than Tol share one. A spatial hash with
// cells of size Tol finds the candidates, so welding stays linear in the vertex count.
std::vector<int> WeldPositions(const std::vector<Vec3>& Positions, double Tol)
{
	const double Cell = Tol > 1e-12 ? Tol : 1e-12;
	const double Inv = 1.0 / Cell;
	std::unordered_map<std::uint64_t, std::vector<int>> Buckets;
	Buckets.reserve(Positions.size());
	std::vector<int> Reps;
	std::vector<int> Map(Positions.size(), -1);
	auto KeyOf = [](std::int64_t X, std::int64_t Y, std::int64_t Z) {
		std::uint64_t H = static_cast<std::uint64_t>(X) * 0x9E3779B97F4A7C15ull;
		H ^= static_cast<std::uint64_t>(Y) * 0xC2B2AE3D27D4EB4Full + (H << 6) + (H >> 2);
		H ^= static_cast<std::uint64_t>(Z) * 0x165667B19E3779F9ull + (H << 6) + (H >> 2);
		return H;
	};
	for (size_t i = 0; i < Positions.size(); ++i)
	{
		const Vec3& P = Positions[i];
		const std::int64_t Cx = static_cast<std::int64_t>(std::floor(P.X * Inv));
		const std::int64_t Cy = static_cast<std::int64_t>(std::floor(P.Y * Inv));
		const std::int64_t Cz = static_cast<std::int64_t>(std::floor(P.Z * Inv));
		int Found = -1;
		for (int Dx = -1; Dx <= 1 && Found < 0; ++Dx)
		{
			for (int Dy = -1; Dy <= 1 && Found < 0; ++Dy)
			{
				for (int Dz = -1; Dz <= 1 && Found < 0; ++Dz)
				{
					const auto It = Buckets.find(KeyOf(Cx + Dx, Cy + Dy, Cz + Dz));
					if (It == Buckets.end())
					{
						continue;
					}
					for (int R : It->second)
					{
						const Vec3 D = Positions[static_cast<size_t>(Reps[static_cast<size_t>(R)])] - P;
						if (D.Dot(D) <= Tol * Tol)
						{
							Found = R;
							break;
						}
					}
				}
			}
		}
		if (Found < 0)
		{
			Found = static_cast<int>(Reps.size());
			Reps.push_back(static_cast<int>(i));
			Buckets[KeyOf(Cx, Cy, Cz)].push_back(Found);
		}
		Map[i] = Found;
	}
	return Map;
}

std::uint64_t EdgeKey(int A, int B)
{
	return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(A)) << 32) | static_cast<std::uint32_t>(B);
}

// Welded triangles with no repeated corner (a triangle collapsed by welding has no area and no edges).
std::vector<int> WeldedTriangles(const Mesh& M, double Tol)
{
	const std::vector<int> Map = WeldPositions(M.Positions, Tol);
	std::vector<int> Out;
	Out.reserve(M.Triangles.size());
	for (size_t t = 0; t + 2 < M.Triangles.size(); t += 3)
	{
		const int A = Map[static_cast<size_t>(M.Triangles[t])];
		const int B = Map[static_cast<size_t>(M.Triangles[t + 1])];
		const int C = Map[static_cast<size_t>(M.Triangles[t + 2])];
		if (A == B || B == C || C == A)
		{
			continue;
		}
		Out.push_back(A);
		Out.push_back(B);
		Out.push_back(C);
	}
	return Out;
}
} // namespace meshbuilder_detail

int Mesh::AddVertex(const Vec3& P, const Vec3& N, const Vec2& UV)
{
	Positions.push_back(P);
	Normals.push_back(N.Normalized());
	UVs.push_back(UV);
	return static_cast<int>(Positions.size()) - 1;
}

void Mesh::AddTriangle(int A, int B, int C)
{
	Triangles.push_back(A);
	Triangles.push_back(B);
	Triangles.push_back(C);
}

void Mesh::AddQuad(int A, int B, int C, int D)
{
	AddTriangle(A, B, C);
	AddTriangle(A, C, D);
}

void Mesh::Append(const Mesh& Other, const Vec3& Offset)
{
	const int Base = static_cast<int>(Positions.size());
	Positions.reserve(Positions.size() + Other.Positions.size());
	for (const Vec3& P : Other.Positions)
	{
		Positions.push_back(P + Offset);
	}
	Normals.insert(Normals.end(), Other.Normals.begin(), Other.Normals.end());
	UVs.insert(UVs.end(), Other.UVs.begin(), Other.UVs.end());
	Triangles.reserve(Triangles.size() + Other.Triangles.size());
	for (int I : Other.Triangles)
	{
		Triangles.push_back(I + Base);
	}
}

void Mesh::Translate(const Vec3& Offset)
{
	for (Vec3& P : Positions)
	{
		P = P + Offset;
	}
}

double Mesh::SignedVolume() const
{
	double V = 0.0;
	for (size_t t = 0; t + 2 < Triangles.size(); t += 3)
	{
		const Vec3& A = Positions[static_cast<size_t>(Triangles[t])];
		const Vec3& B = Positions[static_cast<size_t>(Triangles[t + 1])];
		const Vec3& C = Positions[static_cast<size_t>(Triangles[t + 2])];
		V += A.Dot(B.Cross(C));
	}
	return V / 6.0;
}

bool Mesh::IsClosed(double Tol) const
{
	const std::vector<int> T = meshbuilder_detail::WeldedTriangles(*this, Tol);
	if (T.empty())
	{
		return false;
	}
	std::unordered_map<std::uint64_t, int> Directed;
	Directed.reserve(T.size());
	for (size_t t = 0; t < T.size(); t += 3)
	{
		for (int k = 0; k < 3; ++k)
		{
			++Directed[meshbuilder_detail::EdgeKey(T[t + static_cast<size_t>(k)], T[t + static_cast<size_t>((k + 1) % 3)])];
		}
	}
	for (const auto& E : Directed)
	{
		if (E.second != 1)
		{
			return false;
		}
		const int A = static_cast<int>(E.first >> 32);
		const int B = static_cast<int>(E.first & 0xFFFFFFFFu);
		const auto Rev = Directed.find(meshbuilder_detail::EdgeKey(B, A));
		if (Rev == Directed.end() || Rev->second != 1)
		{
			return false;
		}
	}
	return true;
}

int Mesh::EulerCharacteristic(double Tol) const
{
	const std::vector<int> T = meshbuilder_detail::WeldedTriangles(*this, Tol);
	std::unordered_set<int> Verts;
	std::unordered_set<std::uint64_t> Edges;
	for (size_t t = 0; t < T.size(); t += 3)
	{
		for (int k = 0; k < 3; ++k)
		{
			const int A = T[t + static_cast<size_t>(k)];
			const int B = T[t + static_cast<size_t>((k + 1) % 3)];
			Verts.insert(A);
			Edges.insert(A < B ? meshbuilder_detail::EdgeKey(A, B) : meshbuilder_detail::EdgeKey(B, A));
		}
	}
	return static_cast<int>(Verts.size()) - static_cast<int>(Edges.size()) + static_cast<int>(T.size() / 3);
}

bool Mesh::NormalsAgreeWithWinding(double MinDot) const
{
	if (Normals.size() != Positions.size())
	{
		return false;
	}
	for (size_t t = 0; t + 2 < Triangles.size(); t += 3)
	{
		const int I[3] = {Triangles[t], Triangles[t + 1], Triangles[t + 2]};
		const Vec3 F = (Positions[static_cast<size_t>(I[1])] - Positions[static_cast<size_t>(I[0])]).Cross(Positions[static_cast<size_t>(I[2])] - Positions[static_cast<size_t>(I[0])]);
		const double Len = F.Length();
		if (Len < 1e-12)
		{
			continue;
		}
		const Vec3 Fn = F * (1.0 / Len);
		for (int k = 0; k < 3; ++k)
		{
			if (Normals[static_cast<size_t>(I[k])].Dot(Fn) < MinDot)
			{
				return false;
			}
		}
	}
	return true;
}

bool WriteObj(const std::string& Path, const std::vector<std::pair<std::string, const Mesh*>>& Parts)
{
	std::FILE* F = std::fopen(Path.c_str(), "w");
	if (!F)
	{
		return false;
	}
	std::fprintf(F, "# Mechanical calculator keyboard geometry, millimetres, x right, y away from the typist, z up\n");
	long long Base = 1;
	for (const auto& Part : Parts)
	{
		const Mesh* M = Part.second;
		if (!M)
		{
			continue;
		}
		std::fprintf(F, "o %s\ng %s\nusemtl %s\n", Part.first.c_str(), Part.first.c_str(), Part.first.c_str());
		for (const Vec3& P : M->Positions)
		{
			std::fprintf(F, "v %.4f %.4f %.4f\n", P.X, P.Y, P.Z);
		}
		for (size_t i = 0; i < M->Positions.size(); ++i)
		{
			const Vec2 UV = i < M->UVs.size() ? M->UVs[i] : Vec2();
			std::fprintf(F, "vt %.4f %.4f\n", UV.U, UV.V);
		}
		for (size_t i = 0; i < M->Positions.size(); ++i)
		{
			const Vec3 N = i < M->Normals.size() ? M->Normals[i] : Vec3(0, 0, 1);
			std::fprintf(F, "vn %.4f %.4f %.4f\n", N.X, N.Y, N.Z);
		}
		for (size_t t = 0; t + 2 < M->Triangles.size(); t += 3)
		{
			const long long A = Base + M->Triangles[t], B = Base + M->Triangles[t + 1], C = Base + M->Triangles[t + 2];
			std::fprintf(F, "f %lld/%lld/%lld %lld/%lld/%lld %lld/%lld/%lld\n", A, A, A, B, B, B, C, C, C);
		}
		Base += static_cast<long long>(M->Positions.size());
	}
	return std::fclose(F) == 0;
}
} // namespace mgeo
