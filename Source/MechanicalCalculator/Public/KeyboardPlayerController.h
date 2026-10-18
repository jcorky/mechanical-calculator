// Shows the mouse cursor (keys are clickable) while every keystroke still reaches the game.
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/PlayerController.h"
#include "KeyboardPlayerController.generated.h"

UCLASS()
class MECHANICALCALCULATOR_API AKeyboardPlayerController : public APlayerController
{
	GENERATED_BODY()

public:
	AKeyboardPlayerController();
	virtual void BeginPlay() override;
};
