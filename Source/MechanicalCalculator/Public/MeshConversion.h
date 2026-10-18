// Bridges the engine-independent meshes (millimetres, right-handed, x right / y back / z up) to
// UProceduralMeshComponent sections and runtime static meshes (centimetres, Unreal's left-handed
// X forward / Y right / Z up), and the mechanism's part poses to transforms.
#pragma once

#include "CoreMinimal.h"
#include "GearParts.h"
#include "MeshBuilder.h"
#include "ProceduralMeshComponent.h"

class UMaterialInterface;
class UStaticMesh;

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

// A pose (local part frame -> parent, mm) as a relative transform. Since the vertices were converted
// with ToUnreal, the part's Unreal X axis is its geometry Y axis and the other way round.
inline FTransform ToUnreal(const mgeo::Pose& P)
{
	return FTransform(DirToUnreal(P.Y), DirToUnreal(P.X), DirToUnreal(P.Z), ToUnreal(P.Origin));
}

// One material section of a runtime static mesh.
struct FStaticSection
{
	const mgeo::Mesh* Geometry = nullptr;
	UMaterialInterface* Material = nullptr;
};

// Builds a transient static mesh (no collision) from the sections, one material slot each, for
// instancing parts that repeat. Empty sections are skipped; null if every section is empty.
UStaticMesh* BuildStaticMesh(UObject* Outer, const FName& Name, const TArray<FStaticSection>& Sections, bool bFlipWinding = false);
} // namespace KeyboardMesh
