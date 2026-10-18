#include "MechanicalCalculatorGameMode.h"

#include "Engine/DirectionalLight.h"
#include "Engine/World.h"
#include "EngineUtils.h"
#include "KeyboardPlayerController.h"
#include "KeyboardStudio.h"
#include "KeyboardViewerPawn.h"
#include "MechanicalKeyboard.h"

AMechanicalCalculatorGameMode::AMechanicalCalculatorGameMode()
{
	DefaultPawnClass = AKeyboardViewerPawn::StaticClass();
	PlayerControllerClass = AKeyboardPlayerController::StaticClass();
}

void AMechanicalCalculatorGameMode::StartPlay()
{
	UWorld* World = GetWorld();
	if (World)
	{
		FActorSpawnParameters Params;
		Params.SpawnCollisionHandlingOverride = ESpawnActorCollisionHandlingMethod::AlwaysSpawn;
		const bool bHasLights = TActorIterator<ADirectionalLight>(World) || TActorIterator<AKeyboardStudio>(World);
		if (!bHasLights)
		{
			World->SpawnActor<AKeyboardStudio>(AKeyboardStudio::StaticClass(), FTransform::Identity, Params);
		}
		if (!TActorIterator<AMechanicalKeyboard>(World))
		{
			World->SpawnActor<AMechanicalKeyboard>(AMechanicalKeyboard::StaticClass(), FTransform::Identity, Params);
		}
	}
	Super::StartPlay();
}
