#include "KeyboardMaterials.h"

#include "HAL/IConsoleManager.h"
#include "Materials/Material.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "UObject/Package.h"

#if WITH_EDITOR
#include "AssetRegistry/AssetRegistryModule.h"
#include "Materials/MaterialExpressionAdd.h"
#include "Materials/MaterialExpressionMultiply.h"
#include "Materials/MaterialExpressionScalarParameter.h"
#include "Materials/MaterialExpressionVectorParameter.h"
#include "Materials/MaterialExpressionVertexColor.h"
#include "Misc/PackageName.h"
#include "UObject/SavePackage.h"
#endif

namespace KeyboardMaterials
{
namespace
{
const TCHAR* const OpaqueAsset = TEXT("/Game/MechanicalCalculator/Materials/M_KeyboardOpaque.M_KeyboardOpaque");
const TCHAR* const TranslucentAsset = TEXT("/Game/MechanicalCalculator/Materials/M_KeyboardTranslucent.M_KeyboardTranslucent");
const TCHAR* const FallbackAsset = TEXT("/Engine/BasicShapes/BasicShapeMaterial.BasicShapeMaterial");

UMaterialInterface* GOpaque = nullptr;
UMaterialInterface* GTranslucent = nullptr;
UMaterialInterface* GFallback = nullptr;
bool GResolved = false;
bool GFull = false;

#if WITH_EDITOR
UMaterialExpressionVectorParameter* AddVector(UMaterial* M, const TCHAR* Name, const FLinearColor& Default, int32 X, int32 Y)
{
	UMaterialExpressionVectorParameter* E = NewObject<UMaterialExpressionVectorParameter>(M);
	E->ParameterName = Name;
	E->DefaultValue = Default;
	E->MaterialExpressionEditorX = X;
	E->MaterialExpressionEditorY = Y;
	M->GetExpressionCollection().AddExpression(E);
	return E;
}

UMaterialExpressionScalarParameter* AddScalar(UMaterial* M, const TCHAR* Name, float Default, int32 X, int32 Y)
{
	UMaterialExpressionScalarParameter* E = NewObject<UMaterialExpressionScalarParameter>(M);
	E->ParameterName = Name;
	E->DefaultValue = Default;
	E->MaterialExpressionEditorX = X;
	E->MaterialExpressionEditorY = Y;
	M->GetExpressionCollection().AddExpression(E);
	return E;
}

// Builds a master: BaseColor, Metallic, Roughness, Specular, Opacity (translucent only) and
// Emissive = Emissive + VertexColor.rgb * VertexEmissive.
UMaterial* BuildMaster(UObject* Outer, FName Name, EObjectFlags Flags, bool bTranslucent)
{
	UMaterial* M = NewObject<UMaterial>(Outer, Name, Flags);
	UMaterialEditorOnlyData* D = M->GetEditorOnlyData();
	UMaterialExpressionVectorParameter* Base = AddVector(M, TEXT("BaseColor"), FLinearColor(0.8f, 0.8f, 0.8f), -600, 0);
	UMaterialExpressionScalarParameter* Metal = AddScalar(M, TEXT("Metallic"), 0.0f, -600, 200);
	UMaterialExpressionScalarParameter* Rough = AddScalar(M, TEXT("Roughness"), 0.5f, -600, 300);
	UMaterialExpressionScalarParameter* Spec = AddScalar(M, TEXT("Specular"), 0.5f, -600, 400);
	UMaterialExpressionVectorParameter* Emit = AddVector(M, TEXT("Emissive"), FLinearColor::Black, -900, 500);
	UMaterialExpressionScalarParameter* VtxEmit = AddScalar(M, TEXT("VertexEmissive"), 0.0f, -900, 750);

	UMaterialExpressionVertexColor* Vtx = NewObject<UMaterialExpressionVertexColor>(M);
	Vtx->MaterialExpressionEditorX = -900;
	Vtx->MaterialExpressionEditorY = 650;
	M->GetExpressionCollection().AddExpression(Vtx);
	UMaterialExpressionMultiply* Mul = NewObject<UMaterialExpressionMultiply>(M);
	Mul->MaterialExpressionEditorX = -600;
	Mul->MaterialExpressionEditorY = 650;
	M->GetExpressionCollection().AddExpression(Mul);
	Mul->A.Connect(0, Vtx);
	Mul->B.Connect(0, VtxEmit);
	UMaterialExpressionAdd* Sum = NewObject<UMaterialExpressionAdd>(M);
	Sum->MaterialExpressionEditorX = -300;
	Sum->MaterialExpressionEditorY = 550;
	M->GetExpressionCollection().AddExpression(Sum);
	Sum->A.Connect(0, Emit);
	Sum->B.Connect(0, Mul);

	D->BaseColor.Connect(0, Base);
	D->Metallic.Connect(0, Metal);
	D->Roughness.Connect(0, Rough);
	D->Specular.Connect(0, Spec);
	D->EmissiveColor.Connect(0, Sum);
	if (bTranslucent)
	{
		UMaterialExpressionScalarParameter* Opacity = AddScalar(M, TEXT("Opacity"), 0.35f, -600, 850);
		D->Opacity.Connect(0, Opacity);
		M->BlendMode = BLEND_Translucent;
		M->TranslucencyLightingMode = TLM_SurfacePerPixelLighting;
	}
	M->PreEditChange(nullptr);
	M->PostEditChange();
	return M;
}

void SaveMasters()
{
	struct FEntry
	{
		const TCHAR* Name;
		bool bTranslucent;
	};
	const FEntry Entries[] = {{TEXT("M_KeyboardOpaque"), false}, {TEXT("M_KeyboardTranslucent"), true}};
	for (const FEntry& Entry : Entries)
	{
		const FString PackageName = FString(TEXT("/Game/MechanicalCalculator/Materials/")) + Entry.Name;
		UPackage* Package = CreatePackage(*PackageName);
		UMaterial* M = BuildMaster(Package, FName(Entry.Name), RF_Public | RF_Standalone, Entry.bTranslucent);
		FAssetRegistryModule::AssetCreated(M);
		Package->MarkPackageDirty();
		const FString File = FPackageName::LongPackageNameToFilename(PackageName, FPackageName::GetAssetPackageExtension());
		FSavePackageArgs Args;
		Args.TopLevelFlags = RF_Public | RF_Standalone;
		const bool bSaved = UPackage::SavePackage(Package, M, *File, Args);
		UE_LOG(LogTemp, Display, TEXT("Mechanical Calculator: %s %s"), bSaved ? TEXT("saved") : TEXT("FAILED to save"), *File);
	}
	GResolved = false;
}

FAutoConsoleCommand GSaveMaterialsCommand(
	TEXT("mc.SaveMaterials"),
	TEXT("Saves the keyboard's master materials to /Game/MechanicalCalculator/Materials so packaged builds have them."),
	FConsoleCommandDelegate::CreateStatic(&SaveMasters));
#endif

void Resolve()
{
	if (GResolved)
	{
		return;
	}
	GResolved = true;
	GOpaque = LoadObject<UMaterialInterface>(nullptr, OpaqueAsset, nullptr, LOAD_NoWarn | LOAD_Quiet);
	GTranslucent = LoadObject<UMaterialInterface>(nullptr, TranslucentAsset, nullptr, LOAD_NoWarn | LOAD_Quiet);
#if WITH_EDITOR
	if (!GOpaque)
	{
		GOpaque = BuildMaster(GetTransientPackage(), TEXT("M_KeyboardOpaque_Transient"), RF_Transient, false);
		GOpaque->AddToRoot();
	}
	if (!GTranslucent)
	{
		GTranslucent = BuildMaster(GetTransientPackage(), TEXT("M_KeyboardTranslucent_Transient"), RF_Transient, true);
		GTranslucent->AddToRoot();
	}
#endif
	GFull = GOpaque != nullptr && GTranslucent != nullptr;
	GFallback = LoadObject<UMaterialInterface>(nullptr, FallbackAsset);
}
} // namespace

bool HasFullMasters()
{
	Resolve();
	return GFull;
}

void Apply(UMaterialInstanceDynamic* Material, const FKeyboardSurface& S)
{
	if (!Material)
	{
		return;
	}
	Material->SetVectorParameterValue(TEXT("BaseColor"), S.BaseColor);
	Material->SetVectorParameterValue(TEXT("Color"), S.BaseColor);   // BasicShapeMaterial's parameter
	Material->SetScalarParameterValue(TEXT("Metallic"), S.Metallic);
	Material->SetScalarParameterValue(TEXT("Roughness"), S.Roughness);
	Material->SetScalarParameterValue(TEXT("Specular"), S.Specular);
	Material->SetVectorParameterValue(TEXT("Emissive"), S.Emissive);
	Material->SetScalarParameterValue(TEXT("VertexEmissive"), S.VertexEmissive);
	Material->SetScalarParameterValue(TEXT("Opacity"), S.Opacity);
}

UMaterialInstanceDynamic* Make(UObject* Outer, const FKeyboardSurface& S)
{
	Resolve();
	UMaterialInterface* Parent = GFull ? (S.bTranslucent ? GTranslucent : GOpaque) : GFallback;
	if (!Parent)
	{
		return nullptr;
	}
	UMaterialInstanceDynamic* Mid = UMaterialInstanceDynamic::Create(Parent, Outer);
	Apply(Mid, S);
	return Mid;
}
} // namespace KeyboardMaterials
