#include "MeshConversion.h"

namespace KeyboardMesh
{
FMeshSectionData Convert(const mgeo::Mesh& Mesh, const FLinearColor& Color, bool bFlipWinding)
{
	FMeshSectionData Out;
	const int32 NumVerts = static_cast<int32>(Mesh.Positions.size());
	Out.Vertices.Reserve(NumVerts);
	Out.Normals.Reserve(NumVerts);
	Out.UV0.Reserve(NumVerts);
	Out.Colors.Init(Color, NumVerts);
	for (int32 i = 0; i < NumVerts; ++i)
	{
		Out.Vertices.Add(ToUnreal(Mesh.Positions[i]));
		Out.Normals.Add(i < static_cast<int32>(Mesh.Normals.size()) ? DirToUnreal(Mesh.Normals[i]) : FVector::UpVector);
		Out.UV0.Add(i < static_cast<int32>(Mesh.UVs.size()) ? FVector2D(Mesh.UVs[i].U, Mesh.UVs[i].V) : FVector2D::ZeroVector);
	}
	Out.Triangles.Reserve(static_cast<int32>(Mesh.Triangles.size()));
	for (size_t t = 0; t + 2 < Mesh.Triangles.size(); t += 3)
	{
		Out.Triangles.Add(Mesh.Triangles[t]);
		Out.Triangles.Add(bFlipWinding ? Mesh.Triangles[t + 2] : Mesh.Triangles[t + 1]);
		Out.Triangles.Add(bFlipWinding ? Mesh.Triangles[t + 1] : Mesh.Triangles[t + 2]);
	}
	return Out;
}

void SetSection(UProceduralMeshComponent* Component, int32 Section, const FMeshSectionData& Data)
{
	if (!Component || Data.IsEmpty())
	{
		return;
	}
	Component->CreateMeshSection_LinearColor(Section, Data.Vertices, Data.Triangles, Data.Normals, Data.UV0, Data.Colors, Data.Tangents, false);
}
} // namespace KeyboardMesh
