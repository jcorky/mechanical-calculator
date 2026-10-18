#include "KeyboardViewerPawn.h"

#include "Camera/CameraComponent.h"
#include "Components/InputComponent.h"
#include "EngineUtils.h"
#include "GameFramework/PlayerController.h"
#include "GameFramework/SpringArmComponent.h"
#include "MechanicalKeyboard.h"

AKeyboardViewerPawn::AKeyboardViewerPawn()
{
	PrimaryActorTick.bCanEverTick = true;
	Pivot = CreateDefaultSubobject<USceneComponent>(TEXT("Pivot"));
	SetRootComponent(Pivot);
	Arm = CreateDefaultSubobject<USpringArmComponent>(TEXT("Arm"));
	Arm->SetupAttachment(Pivot);
	Arm->bDoCollisionTest = false;
	Arm->bUsePawnControlRotation = false;
	Arm->bEnableCameraLag = false;
	Arm->TargetArmLength = TargetDistance;
	Arm->SetRelativeRotation(TargetRotation);
	Camera = CreateDefaultSubobject<UCameraComponent>(TEXT("Camera"));
	Camera->SetupAttachment(Arm, USpringArmComponent::SocketName);
	Camera->SetFieldOfView(45.0f);
	AutoPossessPlayer = EAutoReceiveInput::Player0;
}

AMechanicalKeyboard* AKeyboardViewerPawn::FindKeyboard()
{
	if (!Keyboard.IsValid() && GetWorld())
	{
		for (TActorIterator<AMechanicalKeyboard> It(GetWorld()); It; ++It)
		{
			Keyboard = *It;
			TargetPivot = It->GetActorLocation() + FVector(-2.0, 0.0, 2.0);
			SetActorLocation(TargetPivot);
			break;
		}
	}
	return Keyboard.Get();
}

void AKeyboardViewerPawn::ResetView()
{
	TargetDistance = 62.0f;
	TargetRotation = FRotator(-48.0f, 0.0f, 0.0f);
	if (AMechanicalKeyboard* K = FindKeyboard())
	{
		TargetPivot = K->GetActorLocation() + FVector(-2.0, 0.0, 2.0);
	}
}

void AKeyboardViewerPawn::SetupPlayerInputComponent(UInputComponent* PlayerInputComponent)
{
	Super::SetupPlayerInputComponent(PlayerInputComponent);
	PlayerInputComponent->BindKey(EKeys::AnyKey, IE_Pressed, this, &AKeyboardViewerPawn::OnKeyDown);
	PlayerInputComponent->BindKey(EKeys::AnyKey, IE_Released, this, &AKeyboardViewerPawn::OnKeyUp);
	PlayerInputComponent->BindKey(EKeys::AnyKey, IE_Repeat, this, &AKeyboardViewerPawn::OnKeyRepeat);
}

void AKeyboardViewerPawn::OnKeyDown(FKey Key)
{
	if (Key == EKeys::RightMouseButton)
	{
		bOrbiting = true;
		return;
	}
	if (Key == EKeys::MiddleMouseButton)
	{
		bPanning = true;
		return;
	}
	if (Key == EKeys::MouseScrollUp || Key == EKeys::MouseScrollDown)
	{
		const float Factor = Key == EKeys::MouseScrollUp ? 0.88f : 1.0f / 0.88f;
		TargetDistance = FMath::Clamp(TargetDistance * Factor, MinDistance, MaxDistance);
		return;
	}
	if (Key == EKeys::ThumbMouseButton)
	{
		ResetView();
		return;
	}
	if (Key == EKeys::LeftMouseButton)
	{
		APlayerController* PC = Cast<APlayerController>(GetController());
		AMechanicalKeyboard* K = FindKeyboard();
		FHitResult Hit;
		if (PC && K && PC->GetHitResultUnderCursor(ECC_Visibility, false, Hit))
		{
			K->ClickComponent(Hit.GetComponent());
		}
		return;
	}
	if (Key.IsMouseButton() || Key.IsGamepadKey())
	{
		return;
	}
	if (AMechanicalKeyboard* K = FindKeyboard())
	{
		K->PhysicalKeyDown(Key);
	}
}

void AKeyboardViewerPawn::OnKeyUp(FKey Key)
{
	if (Key == EKeys::RightMouseButton)
	{
		bOrbiting = false;
		return;
	}
	if (Key == EKeys::MiddleMouseButton)
	{
		bPanning = false;
		return;
	}
	if (Key.IsMouseButton() || Key.IsGamepadKey())
	{
		return;
	}
	if (AMechanicalKeyboard* K = FindKeyboard())
	{
		K->PhysicalKeyUp(Key);
	}
}

void AKeyboardViewerPawn::OnKeyRepeat(FKey Key)
{
	if (Key.IsMouseButton() || Key.IsGamepadKey())
	{
		return;
	}
	if (AMechanicalKeyboard* K = FindKeyboard())
	{
		K->PhysicalKeyRepeat(Key);
	}
}

void AKeyboardViewerPawn::Tick(float DeltaSeconds)
{
	Super::Tick(DeltaSeconds);
	FindKeyboard();
	if (APlayerController* PC = Cast<APlayerController>(GetController()))
	{
		float DX = 0.0f, DY = 0.0f;
		PC->GetInputMouseDelta(DX, DY);
		if (bOrbiting)
		{
			TargetRotation.Yaw += DX * OrbitSpeed * 10.0f;
			TargetRotation.Pitch = FMath::Clamp(TargetRotation.Pitch + DY * OrbitSpeed * 10.0f, -88.0f, -4.0f);
		}
		else if (bPanning)
		{
			const FRotator Yaw(0.0f, TargetRotation.Yaw, 0.0f);
			const float Scale = TargetDistance * 0.004f * 10.0f;
			TargetPivot += Yaw.RotateVector(FVector(-DY, -DX, 0.0f)) * Scale;
		}
	}
	const float Alpha = 1.0f - FMath::Exp(-12.0f * DeltaSeconds);
	Arm->TargetArmLength = FMath::Lerp(Arm->TargetArmLength, TargetDistance, Alpha);
	Arm->SetRelativeRotation(FMath::Lerp(Arm->GetRelativeRotation(), TargetRotation, Alpha));
	SetActorLocation(FMath::Lerp(GetActorLocation(), TargetPivot, Alpha));
}
