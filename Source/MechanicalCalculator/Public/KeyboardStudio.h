// The set the keyboard sits on when a level has none of its own: a walnut desk, a low warm sun,
// a sky with real-time sky-light capture, two soft area lights, and a post-process volume with
// bloom (the backlight glows) and a light vignette.
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "KeyboardStudio.generated.h"

class UDirectionalLightComponent;
class UExponentialHeightFogComponent;
class UPostProcessComponent;
class URectLightComponent;
class USkyAtmosphereComponent;
class USkyLightComponent;
class UStaticMeshComponent;

UCLASS()
class MECHANICALCALCULATOR_API AKeyboardStudio : public AActor
{
	GENERATED_BODY()

public:
	AKeyboardStudio();
	virtual void BeginPlay() override;

	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<USceneComponent> Root;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<UDirectionalLightComponent> Sun;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<USkyAtmosphereComponent> Sky;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<USkyLightComponent> SkyLight;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<UExponentialHeightFogComponent> Fog;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<URectLightComponent> KeyLight;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<URectLightComponent> RimLight;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<UStaticMeshComponent> Desk;
	UPROPERTY(VisibleAnywhere, Category = "Studio")
	TObjectPtr<UPostProcessComponent> Post;

	UPROPERTY(EditAnywhere, Category = "Studio")
	FLinearColor DeskColor = FLinearColor(0.11f, 0.055f, 0.03f);
};
