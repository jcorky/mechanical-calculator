// The calculating mechanism behind the keyboard: one instanced static mesh per part (the 126
// memory figure wheels alone would be two million procedural-mesh vertices), posed every tick by
// mkb::MechanismAnimator from the mechanism player's state; the printed tape; the memory columns;
// and the maker's plate.
#include "MechanicalKeyboard.h"

#include "Components/InstancedStaticMeshComponent.h"
#include "Components/TextRenderComponent.h"
#include "Engine/StaticMesh.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "MeshConversion.h"
#include "RegisterReadings.h"

namespace MechanicalKeyboardMechanismDetail
{
const FName MechanismTag(TEXT("MechanicalCalculatorBuilt"));   // the same tag Build() clears

struct FFinish
{
	FLinearColor Color;
	float Metallic;
	float Roughness;
	float Specular;
	bool bTranslucent;
	float Opacity;
};

FFinish FinishFor(mgeo::MechMaterial M)
{
	switch (M)
	{
	case mgeo::MechMaterial::Brass: return {FLinearColor(0.92f, 0.70f, 0.36f), 1.0f, 0.26f, 0.5f, false, 1.0f};
	case mgeo::MechMaterial::Steel: return {FLinearColor(0.80f, 0.81f, 0.83f), 1.0f, 0.2f, 0.5f, false, 1.0f};
	case mgeo::MechMaterial::CastIron: return {FLinearColor(0.03f, 0.03f, 0.032f), 0.55f, 0.55f, 0.4f, false, 1.0f};
	case mgeo::MechMaterial::BlackEnamel: return {FLinearColor(0.01f, 0.01f, 0.011f), 0.0f, 0.22f, 0.6f, false, 1.0f};
	case mgeo::MechMaterial::NumeralBlack: return {FLinearColor(0.014f, 0.014f, 0.015f), 0.0f, 0.42f, 0.45f, false, 1.0f};
	case mgeo::MechMaterial::Paint: return {FLinearColor(0.9f, 0.89f, 0.84f), 0.0f, 0.5f, 0.4f, false, 1.0f};
	case mgeo::MechMaterial::Engraving: return {FLinearColor(0.03f, 0.024f, 0.018f), 0.3f, 0.6f, 0.3f, false, 1.0f};
	case mgeo::MechMaterial::WhiteEnamel: return {FLinearColor(0.88f, 0.88f, 0.85f), 0.0f, 0.28f, 0.55f, false, 1.0f};
	case mgeo::MechMaterial::Glass: return {FLinearColor(0.9f, 0.95f, 1.0f), 0.0f, 0.03f, 0.9f, true, 0.12f};
	case mgeo::MechMaterial::Paper: return {FLinearColor(0.84f, 0.82f, 0.75f), 0.0f, 0.85f, 0.3f, false, 1.0f};
	case mgeo::MechMaterial::Rubber: return {FLinearColor(0.02f, 0.02f, 0.02f), 0.0f, 0.9f, 0.2f, false, 1.0f};
	case mgeo::MechMaterial::Ribbon: return {FLinearColor(0.035f, 0.03f, 0.03f), 0.0f, 0.75f, 0.3f, false, 1.0f};
	case mgeo::MechMaterial::Wood: return {FLinearColor(0.30f, 0.14f, 0.06f), 0.0f, 0.42f, 0.45f, false, 1.0f};
	default: return {FLinearColor(0.5f, 0.5f, 0.5f), 0.0f, 0.5f, 0.5f, false, 1.0f};
	}
}

// A label's glyph height so that Len characters fit its width.
double FitHeight(double MaxHeightMm, double WidthMm, int32 Len)
{
	return FMath::Min(MaxHeightMm, WidthMm / (0.62 * FMath::Max(1, Len)));
}
} // namespace MechanicalKeyboardMechanismDetail

using namespace MechanicalKeyboardMechanismDetail;

USceneComponent* AMechanicalKeyboard::NewScene(USceneComponent* Parent, const FName& Name)
{
	USceneComponent* Scene = NewObject<USceneComponent>(this, MakeUniqueObjectName(this, USceneComponent::StaticClass(), Name));
	Scene->CreationMethod = HasActorBegunPlay() ? EComponentCreationMethod::Instance : EComponentCreationMethod::UserConstructionScript;
	Scene->ComponentTags.Add(MechanismTag);
	Scene->SetupAttachment(Parent);
	Scene->RegisterComponent();
	return Scene;
}

UInstancedStaticMeshComponent* AMechanicalKeyboard::NewInstanced(USceneComponent* Parent, const FName& Name, UStaticMesh* Mesh)
{
	UInstancedStaticMeshComponent* Ism = NewObject<UInstancedStaticMeshComponent>(this, MakeUniqueObjectName(this, UInstancedStaticMeshComponent::StaticClass(), Name));
	Ism->CreationMethod = HasActorBegunPlay() ? EComponentCreationMethod::Instance : EComponentCreationMethod::UserConstructionScript;
	Ism->ComponentTags.Add(MechanismTag);
	Ism->SetupAttachment(Parent);
	Ism->SetMobility(EComponentMobility::Movable);
	Ism->SetCollisionEnabled(ECollisionEnabled::NoCollision);
	Ism->SetStaticMesh(Mesh);
	Ism->RegisterComponent();
	return Ism;
}

UTextRenderComponent* AMechanicalKeyboard::NewLabel(USceneComponent* Parent, const mgeo::Vec3& Center, const mgeo::Vec3& Normal, const mgeo::Vec3& Up, double HeightMm,
	const FLinearColor& Color, bool bCentered)
{
	UTextRenderComponent* T = NewText(Parent, FName(TEXT("Label")));
	T->SetRelativeLocation(KeyboardMesh::ToUnreal(Center + Normal * 0.08));
	T->SetRelativeRotation(FRotationMatrix::MakeFromXZ(KeyboardMesh::DirToUnreal(Normal), KeyboardMesh::DirToUnreal(Up)).Rotator());
	T->SetHorizontalAlignment(bCentered ? EHTA_Center : EHTA_Left);
	T->SetWorldSize(static_cast<float>(HeightMm * 0.1));
	T->SetTextRenderColor(Color.ToFColor(true));
	T->SetText(FText::GetEmpty());
	return T;
}

UMaterialInterface* AMechanicalKeyboard::MechSurface(mgeo::MechMaterial Material)
{
	const int32 Key = static_cast<int32>(Material);
	if (UMaterialInterface** Found = MechMaterials.Find(Key))
	{
		return *Found;
	}
	const FFinish F = FinishFor(Material);
	UMaterialInstanceDynamic* Mid = Surface(F.Color, F.Metallic, F.Roughness, F.Specular, F.bTranslucent, F.Opacity);
	MechMaterials.Add(Key, Mid);
	return Mid;
}

UStaticMesh* AMechanicalKeyboard::MechMesh(const mgeo::MechPart& Part, const FName& Name, bool bFlip)
{
	TArray<KeyboardMesh::FStaticSection> Sections;
	for (const mgeo::MechSection& S : Part.Sections)
	{
		KeyboardMesh::FStaticSection& Out = Sections.AddDefaulted_GetRef();
		Out.Geometry = &S.Geometry;
		Out.Material = MechSurface(S.Material);
	}
	return KeyboardMesh::BuildStaticMesh(this, MakeUniqueObjectName(this, UStaticMesh::StaticClass(), Name), Sections, bFlip);
}

void AMechanicalKeyboard::BuildMechanism(const mgeo::CaseParams& CaseP, bool bFlip)
{
	Animator.Reset();
	Mechanism = mgeo::BuildMechanism(mkb::FullSizeAnsi(), CaseP);
	Animator = MakeUnique<mkb::MechanismAnimator>(Mechanism);
	MechanismRoot->SetRelativeLocation(FVector(0.0, 0.0, Mechanism.PlinthHeightMm * 0.1));

	// The static parts, in the desk frame.
	struct FStaticPart
	{
		const mgeo::MechPart* Part;
		const TCHAR* Name;
		bool bShow;
	};
	const FStaticPart Statics[] = {
		{&Mechanism.Plinth, TEXT("Plinth"), true},
		{&Mechanism.Frame, TEXT("Frame"), true},
		{&Mechanism.Printer, TEXT("Printer"), true},
		{&Mechanism.PaperStrip, TEXT("PaperStrip"), true},
		{&Mechanism.Bell, TEXT("Bell"), true},
		{&Mechanism.MemoryFrame, TEXT("MemoryFrame"), true},
		{&Mechanism.GlassCover, TEXT("GlassCover"), bShowGlassCover},
	};
	for (const FStaticPart& S : Statics)
	{
		if (!S.bShow)
		{
			continue;
		}
		if (UStaticMesh* Mesh = MechMesh(*S.Part, FName(S.Name), bFlip))
		{
			NewInstanced(MechanismRoot, FName(S.Name), Mesh)->AddInstance(FTransform::Identity);
		}
	}

	// The moving parts: the carriage's riders hang under the carriage root.
	CarriageRoot = NewScene(MechanismRoot, FName(TEXT("CarriageRoot")));
	CarriageRoot->SetRelativeTransform(KeyboardMesh::ToUnreal(Animator->CarriageFrame()));
	MechParts.Reset();
	for (int32 i = 0; i < mkb::kMovingParts; ++i)
	{
		const mkb::MovingPart P = static_cast<mkb::MovingPart>(i);
		FMechPart& R = MechParts.AddDefaulted_GetRef();
		R.Part = P;
		const FName Name(UTF8_TO_TCHAR(mkb::MovingPartName(P)));
		UStaticMesh* Mesh = MechMesh(mkb::PartSetOf(Mechanism, P).Model, Name, bFlip);
		if (!Mesh)
		{
			continue;
		}
		R.Instances = NewInstanced(mkb::RidesCarriage(Mechanism, P) ? CarriageRoot : MechanismRoot.Get(), Name, Mesh);
		for (const mgeo::Pose& Pose : Animator->Poses(P))
		{
			const FTransform T = KeyboardMesh::ToUnreal(Pose);
			R.Instances->AddInstance(T);
			R.Placed.Add(T);
		}
	}

	// The maker's plate on the plinth's front edge.
	const mgeo::LabelFace& Plate = Mechanism.MakersPlate;
	const FString PlateText = MakersPlateText;
	UTextRenderComponent* PlateLabel = NewLabel(MechanismRoot, Plate.Center, Plate.Normal, Plate.Up, FitHeight(Plate.HeightMm * 0.55, Plate.WidthMm, PlateText.Len()),
		Theme.Engraving, true);
	PlateLabel->SetText(FText::FromString(PlateText));

	// The memory columns' label plates (filled by RefreshMemory).
	MemoryLabels.Reset();
	for (const mgeo::LabelFace& Face : Mechanism.MemoryLabels)
	{
		MemoryLabels.Add(NewLabel(MechanismRoot, Face.Center, Face.Normal, Face.Up, Face.HeightMm * 0.5, Theme.Engraving, true));
	}

	// The tape: one text line per line pitch along the paper, as far as the paper runs.
	double PathMm = 0.0;
	for (size_t i = 1; i < Mechanism.TapePath.size(); ++i)
	{
		PathMm += (Mechanism.TapePath[i] - Mechanism.TapePath[i - 1]).Length();
	}
	const double Pitch = FMath::Max(1.0, Mechanism.TapeLinePitchMm);
	const double GlyphMm = 0.62 * Pitch;
	const int32 Rows = FMath::Clamp(static_cast<int32>(PathMm / Pitch) - 1, 4, 48);
	TapeCols = FMath::Clamp(static_cast<int32>(Mechanism.TapeUsableWidthMm / (0.6 * GlyphMm)), 16, 48);
	TapeLines.Reset();
	for (int32 Row = 0; Row < Rows; ++Row)
	{
		UTextRenderComponent* T = NewText(MechanismRoot, FName(TEXT("TapeLine")));
		T->SetHorizontalAlignment(EHTA_Left);
		T->SetWorldSize(static_cast<float>(GlyphMm * 0.1));
		T->SetTextRenderColor(Theme.Ink.ToFColor(true));
		T->SetText(FText::GetEmpty());
		TapeLines.Add(T);
	}
	LastTape.Reset();
	LastPrinted = Device.PrintedCount();
	PlaceTape();
	RefreshMemory();

	// Frame the plinth: everything stands on it.
	FBox Box(ForceInit);
	for (const mgeo::MechSection& S : Mechanism.Plinth.Sections)
	{
		for (const mgeo::Vec3& P : S.Geometry.Positions)
		{
			Box += KeyboardMesh::ToUnreal(P);
		}
	}
	if (Box.IsValid)
	{
		const FVector Size = Box.GetSize();
		ViewCenter = FVector(Box.GetCenter().X, Box.GetCenter().Y, Mechanism.PlinthHeightMm * 0.1 + 5.0);
		ViewDistance = static_cast<float>(1.15 * FMath::Max(Size.Y, Size.X * 1.6) * 0.5 / FMath::Tan(FMath::DegreesToRadians(22.5)));
	}
}

void AMechanicalKeyboard::PlaceTape()
{
	if (!Animator)
	{
		return;
	}
	const double HalfWidth = Mechanism.TapeUsableWidthMm * 0.5;
	for (int32 N = 0; N < TapeLines.Num(); ++N)
	{
		const mgeo::TapeFrame F = Animator->TapeLine(N);
		const mgeo::Vec3 P = F.Origin - F.Right * HalfWidth + F.Normal * 0.06;
		TapeLines[N]->SetRelativeLocationAndRotation(KeyboardMesh::ToUnreal(P),
			FRotationMatrix::MakeFromXZ(KeyboardMesh::DirToUnreal(F.Normal), KeyboardMesh::DirToUnreal(F.Up)).Rotator());
	}
}

void AMechanicalKeyboard::FeedPaper()
{
	const uint64 Printed = Device.PrintedCount();
	if (Animator && Printed > LastPrinted)
	{
		const std::vector<std::string>& Tape = Device.Tape();
		const size_t New = static_cast<size_t>(FMath::Min<uint64>(Printed - LastPrinted, Tape.size()));
		int32 Rows = 0;
		for (size_t i = Tape.size() - New; i < Tape.size(); ++i)
		{
			Rows += FMath::Max(1, static_cast<int32>(mcalc::WrapText(Tape[i], TapeCols).size()));
		}
		Animator->Feed(Rows);
	}
	LastPrinted = Printed;
}

void AMechanicalKeyboard::RefreshMemory()
{
	if (!Animator || bComputing)
	{
		return;   // the worker owns the session while the machine computes
	}
	const std::vector<mkb::MemoryColumn> Columns = mkb::ReadMemory(Device.GetSession(), Mechanism.MemoryColumns);
	Animator->ShowMemory(Columns);
	if (!Columns.empty() && Columns[0].Label.rfind("ans", 0) == 0)
	{
		Animator->ShowTypeWheels(Columns[0].Reading);
	}
	for (int32 C = 0; C < MemoryLabels.Num(); ++C)
	{
		FString Text;
		if (C < static_cast<int32>(Columns.size()))
		{
			const mkb::MemoryColumn& Col = Columns[static_cast<size_t>(C)];
			Text = FString(UTF8_TO_TCHAR(Col.Label.c_str()));
			// One sign wheel per column carries the mantissa's sign; a negative power of ten is marked here.
			if (!Col.Reading.Empty && !Col.Reading.Special && Col.Reading.ExponentNegative)
			{
				Text += TEXT(" e-");
			}
		}
		if (C < static_cast<int32>(Mechanism.MemoryLabels.size()) && MemoryLabels[C])
		{
			const mgeo::LabelFace& Face = Mechanism.MemoryLabels[static_cast<size_t>(C)];
			MemoryLabels[C]->SetWorldSize(static_cast<float>(FitHeight(Face.HeightMm * 0.5, Face.WidthMm, Text.Len()) * 0.1));
			MemoryLabels[C]->SetText(FText::FromString(Text));
		}
	}
}

void AMechanicalKeyboard::TickMechanism(float DeltaSeconds)
{
	if (!Animator)
	{
		return;
	}
	FeedPaper();
	if (!Device.Busy())
	{
		// The setting register follows the number being typed, its decimal marker too.
		const int32 Places = mkb::DecimalPlacesInInput(Device.Input(), Device.Cursor());
		Animator->SetDecimalMarkers(Places, Places, 0);
	}
	Animator->Follow(Player.State(), SpinAngle);
	Animator->Step(DeltaSeconds);
	if (CarriageRoot)
	{
		CarriageRoot->SetRelativeTransform(KeyboardMesh::ToUnreal(Animator->CarriageFrame()));
	}
	for (FMechPart& R : MechParts)
	{
		if (!R.Instances)
		{
			continue;
		}
		const std::vector<mgeo::Pose>& Poses = Animator->Poses(R.Part);
		bool bChanged = false;
		for (int32 i = 0; i < R.Placed.Num() && i < static_cast<int32>(Poses.size()); ++i)
		{
			const FTransform T = KeyboardMesh::ToUnreal(Poses[static_cast<size_t>(i)]);
			if (!T.Equals(R.Placed[i], 1e-5))
			{
				R.Placed[i] = T;
				bChanged = true;
			}
		}
		if (bChanged)
		{
			R.Instances->BatchUpdateInstancesTransforms(0, R.Placed, false, true, false);
		}
	}
	// The printed lines ride on the paper while it feeds, and settle where it stops.
	const bool bFeeding = Animator->PaperLag() > 0.0;
	if (bFeeding || bPaperMoving)
	{
		PlaceTape();
	}
	bPaperMoving = bFeeding;
}
