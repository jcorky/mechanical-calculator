#include "KeyboardPlayerController.h"

AKeyboardPlayerController::AKeyboardPlayerController()
{
	bShowMouseCursor = true;
	bEnableClickEvents = false;
	DefaultMouseCursor = EMouseCursor::Hand;
}

void AKeyboardPlayerController::BeginPlay()
{
	Super::BeginPlay();
	FInputModeGameAndUI Mode;
	Mode.SetHideCursorDuringCapture(false);
	Mode.SetLockMouseToViewportBehavior(EMouseLockMode::DoNotLock);
	SetInputMode(Mode);
}
