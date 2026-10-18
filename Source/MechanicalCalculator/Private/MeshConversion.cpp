#include "MeshConversion.h"

#include "Engine/StaticMesh.h"
#include "MeshDescription.h"
#include "StaticMeshAttributes.h"

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

namespace KeyboardMesh
{
UStaticMesh* BuildStaticMesh(UObject* Outer, const FName& Name, const TArray<FStaticSection>& Sections, bool bFlipWinding)
{
	int32 Total = 0;
	for (const FStaticSection& S : Sections)
	{
		if (S.Geometry && !S.Geometry->Triangles.empty())
		{
			Total += static_cast<int32>(S.Geometry->Positions.size());
		}
	}
	if (Total == 0)
	{
		return nullptr;
	}
	FMeshDescription Description;
	FStaticMeshAttributes Attributes(Description);
	Attributes.Register();
	TVertexAttributesRef<FVector3f> Positions = Attributes.GetVertexPositions();
	TVertexInstanceAttributesRef<FVector3f> Normals = Attributes.GetVertexInstanceNormals();
	TVertexInstanceAttributesRef<FVector3f> Tangents = Attributes.GetVertexInstanceTangents();
	TVertexInstanceAttributesRef<float> BinormalSigns = Attributes.GetVertexInstanceBinormalSigns();
	TVertexInstanceAttributesRef<FVector2f> UVs = Attributes.GetVertexInstanceUVs();
	TPolygonGroupAttributesRef<FName> SlotNames = Attributes.GetPolygonGroupMaterialSlotNames();
	Description.ReserveNewVertices(Total);
	Description.ReserveNewVertexInstances(Total);

	UStaticMesh* Mesh = NewObject<UStaticMesh>(Outer, Name, RF_Transient);
	TArray<FVertexInstanceID> Instances;
	for (int32 s = 0; s < Sections.Num(); ++s)
	{
		const mgeo::Mesh* Src = Sections[s].Geometry;
		if (!Src || Src->Triangles.empty())
		{
			continue;
		}
		const FName Slot(*FString::Printf(TEXT("Section%d"), s));
		const FPolygonGroupID Group = Description.CreatePolygonGroup();
		SlotNames[Group] = Slot;
		Mesh->GetStaticMaterials().Add(FStaticMaterial(Sections[s].Material, Slot));

		const int32 NumVerts = static_cast<int32>(Src->Positions.size());
		Instances.Reset(NumVerts);
		for (int32 i = 0; i < NumVerts; ++i)
		{
			const FVertexID V = Description.CreateVertex();
			Positions[V] = FVector3f(ToUnreal(Src->Positions[static_cast<size_t>(i)]));
			const FVertexInstanceID I = Description.CreateVertexInstance(V);
			const FVector N = i < static_cast<int32>(Src->Normals.size()) ? DirToUnreal(Src->Normals[static_cast<size_t>(i)]).GetSafeNormal() : FVector::UpVector;
			// No normal maps here, so any tangent perpendicular to the normal will do.
			const FVector T = (FMath::Abs(N.Z) < 0.9 ? FVector::CrossProduct(N, FVector::UpVector) : FVector::CrossProduct(N, FVector::ForwardVector)).GetSafeNormal();
			Normals[I] = FVector3f(N);
			Tangents[I] = FVector3f(T);
			BinormalSigns[I] = 1.0f;
			const FVector2D UV = i < static_cast<int32>(Src->UVs.size()) ? FVector2D(Src->UVs[static_cast<size_t>(i)].U, Src->UVs[static_cast<size_t>(i)].V) : FVector2D::ZeroVector;
			UVs.Set(I, 0, FVector2f(UV));
			Instances.Add(I);
		}
		for (size_t t = 0; t + 2 < Src->Triangles.size(); t += 3)
		{
			const int32 A = Src->Triangles[t];
			const int32 B = bFlipWinding ? Src->Triangles[t + 2] : Src->Triangles[t + 1];
			const int32 C = bFlipWinding ? Src->Triangles[t + 1] : Src->Triangles[t + 2];
			if (A == B || B == C || A == C)
			{
				continue;
			}
			FVertexInstanceID Tri[3] = {Instances[A], Instances[B], Instances[C]};
			Description.CreateTriangle(Group, MakeArrayView(Tri, 3));
		}
	}
	UStaticMesh::FBuildMeshDescriptionsParams Params;
	Params.bBuildSimpleCollision = false;
	Params.bFastBuild = true;
	TArray<const FMeshDescription*> Descriptions;
	Descriptions.Add(&Description);
	Mesh->BuildFromMeshDescriptions(Descriptions, Params);
	return Mesh;
}
} // namespace KeyboardMesh
