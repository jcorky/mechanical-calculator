// The camera: orbits the keyboard with the right mouse button, pans with the middle button, zooms
// with the wheel, and presses the key under the cursor on a left click. Every keyboard key is
// forwarded to the keyboard, so the real keyboard types on the virtual one.
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Pawn.h"
#include "InputCoreTypes.h"
#include "KeyboardViewerPawn.generated.h"

class AMechanicalKeyboard;
class UCameraComponent;
class USpringArmComponent;

UCLASS()
class MECHANICALCALCULATOR_API AKeyboardViewerPawn : public APawn
{
	GENERATED_BODY()

public:
	AKeyboardViewerPawn();

	virtual void Tick(float DeltaSeconds) override;
	virtual void SetupPlayerInputComponent(UInputComponent* PlayerInputComponent) override;

	UPROPERTY(VisibleAnywhere, Category = "Camera")
	TObjectPtr<USceneComponent> Pivot;

	UPROPERTY(VisibleAnywhere, Category = "Camera")
	TObjectPtr<USpringArmComponent> Arm;

	UPROPERTY(VisibleAnywhere, Category = "Camera")
	TObjectPtr<UCameraComponent> Camera;

	UPROPERTY(EditAnywhere, Category = "Camera")
	float OrbitSpeed = 0.35f;

	UPROPERTY(EditAnywhere, Category = "Camera")
	float MinDistance = 12.0f;

	UPROPERTY(EditAnywhere, Category = "Camera")
	float MaxDistance = 220.0f;

private:
	void OnKeyDown(FKey Key);
	void OnKeyUp(FKey Key);
	void OnKeyRepeat(FKey Key);
	AMechanicalKeyboard* FindKeyboard();
	void ResetView();

	TWeakObjectPtr<AMechanicalKeyboard> Keyboard;
	bool bOrbiting = false;
	bool bPanning = false;
	float TargetDistance = 62.0f;
	FRotator TargetRotation = FRotator(-48.0f, 0.0f, 0.0f);
	FVector TargetPivot = FVector::ZeroVector;
};
