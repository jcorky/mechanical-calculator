#include "KeyboardStudio.h"

#include "Components/DirectionalLightComponent.h"
#include "Components/ExponentialHeightFogComponent.h"
#include "Components/PostProcessComponent.h"
#include "Components/RectLightComponent.h"
#include "Components/SkyAtmosphereComponent.h"
#include "Components/SkyLightComponent.h"
#include "Components/StaticMeshComponent.h"
#include "Engine/StaticMesh.h"
#include "KeyboardMaterials.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "UObject/ConstructorHelpers.h"

AKeyboardStudio::AKeyboardStudio()
{
	Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	SetRootComponent(Root);

	Sun = CreateDefaultSubobject<UDirectionalLightComponent>(TEXT("Sun"));
	Sun->SetupAttachment(Root);
	Sun->SetMobility(EComponentMobility::Movable);
	Sun->SetRelativeRotation(FRotator(-28.0f, 35.0f, 0.0f));
	Sun->SetIntensity(6.0f);
	Sun->SetLightColor(FLinearColor(1.0f, 0.86f, 0.70f));
	Sun->SetAtmosphereSunLight(true);
	Sun->SetDynamicShadowCascades(4);

	Sky = CreateDefaultSubobject<USkyAtmosphereComponent>(TEXT("Sky"));
	Sky->SetupAttachment(Root);

	SkyLight = CreateDefaultSubobject<USkyLightComponent>(TEXT("SkyLight"));
	SkyLight->SetupAttachment(Root);
	SkyLight->SetMobility(EComponentMobility::Movable);
	SkyLight->bRealTimeCapture = true;
	SkyLight->SetIntensity(1.0f);

	Fog = CreateDefaultSubobject<UExponentialHeightFogComponent>(TEXT("Fog"));
	Fog->SetupAttachment(Root);
	Fog->SetFogDensity(0.004f);

	KeyLight = CreateDefaultSubobject<URectLightComponent>(TEXT("KeyLight"));
	KeyLight->SetupAttachment(Root);
	KeyLight->SetMobility(EComponentMobility::Movable);
	KeyLight->SetRelativeLocation(FVector(-60.0f, -70.0f, 90.0f));
	KeyLight->SetRelativeRotation(FRotator(-45.0f, 50.0f, 0.0f));
	KeyLight->SetIntensityUnits(ELightUnits::Lumens);
	KeyLight->SetIntensity(900.0f);
	KeyLight->SetSourceWidth(80.0f);
	KeyLight->SetSourceHeight(80.0f);
	KeyLight->SetLightColor(FLinearColor(1.0f, 0.95f, 0.9f));

	RimLight = CreateDefaultSubobject<URectLightComponent>(TEXT("RimLight"));
	RimLight->SetupAttachment(Root);
	RimLight->SetMobility(EComponentMobility::Movable);
	RimLight->SetRelativeLocation(FVector(80.0f, 60.0f, 50.0f));
	RimLight->SetRelativeRotation(FRotator(-30.0f, -140.0f, 0.0f));
	RimLight->SetIntensityUnits(ELightUnits::Lumens);
	RimLight->SetIntensity(500.0f);
	RimLight->SetSourceWidth(120.0f);
	RimLight->SetSourceHeight(30.0f);
	RimLight->SetLightColor(FLinearColor(0.75f, 0.85f, 1.0f));

	Desk = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("Desk"));
	Desk->SetupAttachment(Root);
	static ConstructorHelpers::FObjectFinder<UStaticMesh> Cube(TEXT("/Engine/BasicShapes/Cube.Cube"));
	if (Cube.Succeeded())
	{
		Desk->SetStaticMesh(Cube.Object);
	}
	// The cube is 100 cm on a side and centred: a 3.2 x 2.0 m desk top, 4 cm thick, its surface at z = 0.
	Desk->SetRelativeScale3D(FVector(2.0f, 3.2f, 0.04f));
	Desk->SetRelativeLocation(FVector(20.0f, 0.0f, -2.0f));
	Desk->SetCollisionProfileName(TEXT("BlockAll"));

	Post = CreateDefaultSubobject<UPostProcessComponent>(TEXT("Post"));
	Post->SetupAttachment(Root);
	Post->bUnbound = true;
	Post->Settings.bOverride_BloomIntensity = true;
	Post->Settings.BloomIntensity = 0.9f;
	Post->Settings.bOverride_VignetteIntensity = true;
	Post->Settings.VignetteIntensity = 0.35f;
	Post->Settings.bOverride_AutoExposureBias = true;
	Post->Settings.AutoExposureBias = 0.4f;
}

void AKeyboardStudio::BeginPlay()
{
	Super::BeginPlay();
	FKeyboardSurface Wood;
	Wood.BaseColor = DeskColor;
	Wood.Roughness = 0.55f;
	Wood.Specular = 0.35f;
	if (UMaterialInstanceDynamic* Mid = KeyboardMaterials::Make(this, Wood))
	{
		Desk->SetMaterial(0, Mid);
	}
}
