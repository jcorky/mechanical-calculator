// Bridges the engine-independent meshes (millimetres, right-handed, x right / y back / z up) to
// UProceduralMeshComponent sections (centimetres, Unreal's left-handed X forward / Y right / Z up).
#pragma once

#include "CoreMinimal.h"
#include "MeshBuilder.h"
#include "ProceduralMeshComponent.h"

struct FMeshSectionData
{
	TArray<FVector> Vertices;
	TArray<int32> Triangles;
	TArray<FVector> Normals;
	TArray<FVector2D> UV0;
	TArray<FLinearColor> Colors;
	TArray<FProcMeshTangent> Tangents;

	bool IsEmpty() const { return Triangles.Num() == 0; }
};

namespace KeyboardMesh
{
// Geometry frame (mm) to Unreal (cm). Swapping x and y converts the handedness, so the triangle
// order that faces outward in the geometry frame faces outward in Unreal too.
inline FVector ToUnreal(const mgeo::Vec3& P) { return FVector(P.Y, P.X, P.Z) * 0.1; }
inline FVector DirToUnreal(const mgeo::Vec3& D) { return FVector(D.Y, D.X, D.Z); }

FMeshSectionData Convert(const mgeo::Mesh& Mesh, const FLinearColor& Color = FLinearColor::Black, bool bFlipWinding = false);

// Creates (or replaces) a section; collision is off, the board uses simple hit boxes instead.
void SetSection(UProceduralMeshComponent* Component, int32 Section, const FMeshSectionData& Data);
} // namespace KeyboardMesh
