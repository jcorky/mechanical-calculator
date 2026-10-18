// Puts a keyboard, the studio set and the viewer into any level that lacks them, so the project
// runs from the engine's empty Entry map with no content of its own.
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/GameModeBase.h"
#include "MechanicalCalculatorGameMode.generated.h"

UCLASS()
class MECHANICALCALCULATOR_API AMechanicalCalculatorGameMode : public AGameModeBase
{
	GENERATED_BODY()

public:
	AMechanicalCalculatorGameMode();
	virtual void StartPlay() override;
};
