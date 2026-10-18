// The keyboard actor: builds the keyboard and the calculating mechanism behind it from the
// procedural geometry, animates each key's switch, spring and stabilizer, runs every evaluation on
// the stepped-drum gears and plays their motion, prints the tape, turns the memory columns, plays
// the synthesized sounds, drives the RGB backlight, and feeds keystrokes to the calculator firmware
// (Keyboard/CalculatorDevice.h).
#pragma once

#include "Async/Future.h"
#include "Backlight.h"
#include "CalculatorDevice.h"
#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "InputCoreTypes.h"
#include "KeyboardAudio.h"
#include "GearParts.h"
#include "KeyboardParts.h"
#include "MechanismAnimator.h"
#include "MechanismPlayer.h"
#include "SwitchModel.h"
#include "MechanicalKeyboard.generated.h"

class UBoxComponent;
class UInstancedStaticMeshComponent;
class UMaterialInterface;
class UStaticMesh;
class UMaterialInstanceDynamic;
class UPrimitiveComponent;
class UProceduralMeshComponent;
class UTextRenderComponent;

UENUM(BlueprintType)
enum class EKeyboardSwitch : uint8
{
	Linear UMETA(DisplayName = "Linear (red)"),
	Tactile UMETA(DisplayName = "Tactile (brown)"),
	Clicky UMETA(DisplayName = "Clicky (blue)"),
};

UENUM(BlueprintType)
enum class EKeyboardKeycaps : uint8
{
	Sculpted UMETA(DisplayName = "Sculpted (Cherry profile)"),
	Vintage UMETA(DisplayName = "Vintage (round, chrome rims)"),
};

USTRUCT(BlueprintType)
struct FKeyboardTheme
{
	GENERATED_BODY()

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keycaps")
	FLinearColor AlphaCaps = FLinearColor(0.80f, 0.76f, 0.66f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keycaps")
	FLinearColor ModifierCaps = FLinearColor(0.34f, 0.35f, 0.37f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keycaps")
	FLinearColor AccentCaps = FLinearColor(0.85f, 0.22f, 0.05f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keycaps")
	FLinearColor FunctionCaps = FLinearColor(0.16f, 0.165f, 0.18f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keycaps")
	FLinearColor NumpadCaps = FLinearColor(0.80f, 0.76f, 0.66f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keycaps")
	FLinearColor NumpadOperatorCaps = FLinearColor(0.85f, 0.45f, 0.06f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Legends")
	FLinearColor LegendOnLight = FLinearColor(0.06f, 0.06f, 0.07f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Legends")
	FLinearColor LegendOnDark = FLinearColor(0.86f, 0.84f, 0.78f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Legends")
	FLinearColor FrontLegend = FLinearColor(0.95f, 0.45f, 0.08f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Body")
	FLinearColor Case = FLinearColor(0.055f, 0.058f, 0.066f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Body")
	FLinearColor Plate = FLinearColor(0.72f, 0.72f, 0.74f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Body")
	FLinearColor Screen = FLinearColor(0.42f, 0.92f, 1.0f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Vintage")
	FLinearColor VintageCaps = FLinearColor(0.012f, 0.012f, 0.013f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Vintage")
	FLinearColor VintageOperatorCaps = FLinearColor(0.42f, 0.015f, 0.012f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Vintage")
	FLinearColor VintageLegend = FLinearColor(0.92f, 0.91f, 0.86f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mechanism")
	FLinearColor Ink = FLinearColor(0.02f, 0.02f, 0.035f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mechanism")
	FLinearColor ErrorInk = FLinearColor(0.45f, 0.02f, 0.02f);
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mechanism")
	FLinearColor Engraving = FLinearColor(0.025f, 0.02f, 0.015f);
};

// What comes back from the machine for one line: the result, the work, and the first part of the
// literal machine's record for the mechanism to play.
struct FGearEvaluation
{
	mcalc::Outcome Outcome;
	mcalc::GearStats Work;
	std::vector<mcalc::GearEvent> Events;
};

UCLASS()
class MECHANICALCALCULATOR_API AMechanicalKeyboard : public AActor
{
	GENERATED_BODY()

public:
	AMechanicalKeyboard();

	virtual void OnConstruction(const FTransform& Transform) override;
	virtual void BeginPlay() override;
	virtual void EndPlay(const EEndPlayReason::Type EndPlayReason) override;
	virtual void Tick(float DeltaSeconds) override;

	// Physical keyboard events, by engine key (the viewer pawn forwards them).
	void PhysicalKeyDown(const FKey& Key);
	void PhysicalKeyUp(const FKey& Key);
	void PhysicalKeyRepeat(const FKey& Key);
	// A mouse click on a key's hit box: press and release (modifiers latch until the next key).
	bool ClickComponent(const UPrimitiveComponent* Component);

	// Press and release a key by its layout id ("A", "Num7", "Enter", ...).
	UFUNCTION(BlueprintCallable, Category = "Mechanical Calculator")
	void TapKey(const FString& Id);

	// Types text on the keys, holding shift where needed, then presses Enter if bEnter.
	UFUNCTION(BlueprintCallable, Category = "Mechanical Calculator")
	void TypeLine(const FString& Text, bool bEnter = true);

	// Evaluates directly (no animation) and returns what the display would print.
	UFUNCTION(BlueprintCallable, Category = "Mechanical Calculator")
	FString Evaluate(const FString& Expression);

	UFUNCTION(BlueprintCallable, Category = "Mechanical Calculator")
	void SetExploded(bool bExploded);

	UFUNCTION(BlueprintPure, Category = "Mechanical Calculator")
	TArray<FString> GetScreenLines() const;

	// Where a camera should look to frame the whole device, and from how far (world, cm).
	UFUNCTION(BlueprintPure, Category = "Mechanical Calculator")
	FVector GetViewPivot() const;
	UFUNCTION(BlueprintPure, Category = "Mechanical Calculator")
	float GetViewDistance() const;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Keyboard")
	EKeyboardSwitch SwitchType = EKeyboardSwitch::Tactile;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Keyboard")
	EKeyboardKeycaps Keycaps = EKeyboardKeycaps::Sculpted;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Keyboard")
	FKeyboardTheme Theme;

	// The OLED above the numpad. Without it a brass name plate sits there and the tape is the
	// display (the OLED is always fitted when the mechanism is hidden).
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Keyboard")
	bool bOledModule = false;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keyboard")
	bool bShowLegends = true;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keyboard", meta = (ClampMin = "0", ClampMax = "2"))
	float SoundVolume = 0.8f;

	// Minimum visible press for an instantaneous keystroke (mouse clicks, fast typing), seconds.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keyboard", meta = (ClampMin = "0.01", ClampMax = "0.5"))
	float MinPressTime = 0.07f;

	// Lift of the top layer in the exploded view, cm.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Keyboard", meta = (ClampMin = "0"))
	float ExplodeHeight = 16.0f;

	// Types these lines on the keys after BeginPlay until the first real keystroke.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Demo")
	bool bPlayDemo = true;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Demo")
	TArray<FString> DemoLines;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Demo", meta = (ClampMin = "1"))
	float DemoCharsPerSecond = 12.0f;

	// Motor speed of the calculating mechanism, crank turns per second.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mechanism", meta = (ClampMin = "1", ClampMax = "200"))
	float MotorTurnsPerSecond = 12.0f;

	// Long evaluations play faster so the mechanism never runs longer than this before printing.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mechanism", meta = (ClampMin = "0.5"))
	float MaxPlaybackSeconds = 6.0f;

	// How much of each evaluation the literal machine records for playback (events).
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mechanism", meta = (ClampMin = "100"))
	int32 RecordedEvents = 6000;

	// The stepped-drum machine, printer and memory columns on their plinth behind the keyboard. The
	// arithmetic runs on the gear model either way; this only shows it.
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Mechanism")
	bool bShowMechanism = true;

	// The glass top over the drum bank (the mechanism is shown open by default).
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Mechanism")
	bool bShowGlassCover = false;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Mechanism")
	FString MakersPlateText = TEXT("MECHANICAL CALCULATOR  -  WHERE IT ISN'T");

	// Reverses every triangle, in case an engine version expects the opposite front face.
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "Advanced")
	bool bFlipWinding = false;

	UPROPERTY(VisibleAnywhere, Category = "Components")
	TObjectPtr<USceneComponent> Root;

	UPROPERTY(VisibleAnywhere, Category = "Components")
	TObjectPtr<USceneComponent> BoardRoot;

	// The desk frame of the mechanism: the plinth's top, under the keyboard's origin.
	UPROPERTY(VisibleAnywhere, Category = "Components")
	TObjectPtr<USceneComponent> MechanismRoot;

private:
	struct FKeyRuntime
	{
		int32 LayoutIndex = INDEX_NONE;
		USceneComponent* Mover = nullptr;          // the moving stack: stem, stabilizer stems, keycap
		UProceduralMeshComponent* Cap = nullptr;
		UProceduralMeshComponent* Stem = nullptr;
		UProceduralMeshComponent* Spring = nullptr;
		UBoxComponent* Hit = nullptr;
		FVector MoverBase = FVector::ZeroVector;   // board-local, cm
		FVector SpringBase = FVector::ZeroVector;
		mkb::SwitchModel Switch;
		double PressedAt = -1.0;
		bool bReleasePending = false;
		bool bLatched = false;
		bool bSettled = false;            // placed at rest; skipped until it moves again
		int32 HousingVertexStart = 0;
		int32 HousingVertexCount = 0;
	};

	struct FLayer
	{
		UProceduralMeshComponent* Mesh = nullptr;
		float ExplodeFraction = 0.0f;              // of ExplodeHeight
	};

	struct FDemoStep
	{
		int32 KeyIndex = INDEX_NONE;
		bool bShift = false;
		double Delay = 0.0;
	};

	// One moving part of the mechanism: one instanced mesh, one instance per pose.
	struct FMechPart
	{
		mkb::MovingPart Part = mkb::MovingPart::Count;
		UInstancedStaticMeshComponent* Instances = nullptr;
		TArray<FTransform> Placed;
	};

	void Build();
	void DestroyBuilt();
	void BuildMechanism(const mgeo::CaseParams& CaseP, bool bFlip);
	UMaterialInterface* MechSurface(mgeo::MechMaterial Material);
	UStaticMesh* MechMesh(const mgeo::MechPart& Part, const FName& Name, bool bFlip);
	UInstancedStaticMeshComponent* NewInstanced(USceneComponent* Parent, const FName& Name, UStaticMesh* Mesh);
	USceneComponent* NewScene(USceneComponent* Parent, const FName& Name);
	UTextRenderComponent* NewLabel(USceneComponent* Parent, const mgeo::Vec3& Center, const mgeo::Vec3& Normal, const mgeo::Vec3& Up, double HeightMm, const FLinearColor& Color, bool bCentered);
	void TickMechanism(float DeltaSeconds);
	void PlaceTape();
	void RefreshMemory();
	void FeedPaper();
	UProceduralMeshComponent* NewMesh(USceneComponent* Parent, const FName& Name);
	UTextRenderComponent* NewText(USceneComponent* Parent, const FName& Name);
	UMaterialInstanceDynamic* Surface(const FLinearColor& Color, float Metallic, float Roughness, float Specular = 0.5f, bool bTranslucent = false, float Opacity = 1.0f,
		const FLinearColor& Emissive = FLinearColor::Black, float VertexEmissive = 0.0f);
	void AddLegends(FKeyRuntime& Runtime, const mkb::KeyDef& Def, const mgeo::KeycapTop& Top);
	FLinearColor CapColor(const mkb::KeyDef& Def) const;
	bool CapIsDark(const mkb::KeyDef& Def) const;
	mkb::SwitchSpec Spec() const;

	int32 IndexOfKey(const FString& Id) const;
	void TapIndex(int32 Index);
	int32 IndexOfEngineKey(const FKey& Key) const;
	void LogicalPress(int32 Index);
	void LogicalRelease(int32 Index);
	void VisualPress(int32 Index);
	void VisualRelease(int32 Index);
	void HandleEvents(const mkb::DeviceEvents& Events);
	void ReleaseLatched();
	void CancelDemo();
	void QueueDemo();

	void StartEvaluation();
	void TickMachine(float DeltaSeconds);
	void TickKeys(float DeltaSeconds);
	void TickExplode(float DeltaSeconds);
	void TickBacklight();
	void TickDemo(float DeltaSeconds);
	void RefreshScreen();

	bool bBuilt = false;
	TArray<FKeyRuntime> Keys;
	TArray<FLayer> Layers;
	UProceduralMeshComponent* Tops = nullptr;     // translucent housings, vertex-coloured by the backlight
	TArray<FLinearColor> TopColors;
	TArray<UTextRenderComponent*> ScreenLines;
	int32 ScreenCols = 32;
	int32 ScreenRows = 5;
	UTextRenderComponent* StatusLine = nullptr;     // on the name plate

	// The mechanism. The animator keeps a reference to the parts, so Mechanism is rebuilt in place.
	mgeo::MechanismParts Mechanism;
	TUniquePtr<mkb::MechanismAnimator> Animator;
	TArray<FMechPart> MechParts;
	USceneComponent* CarriageRoot = nullptr;
	TMap<int32, UMaterialInterface*> MechMaterials;
	TArray<UTextRenderComponent*> TapeLines;          // [0] is the print line
	TArray<FString> LastTape;
	bool bPaperMoving = false;
	TArray<UTextRenderComponent*> MemoryLabels;
	int32 TapeCols = 34;
	uint64 LastPrinted = 0;
	FVector ViewCenter = FVector(0.0, 0.0, 2.0);       // actor-local, cm
	float ViewDistance = 62.0f;
	float StemsFraction = 0.75f;
	float CapsFraction = 1.0f;
	float SpringsFraction = 0.36f;
	double SpringRestCm = 1.0;

	UPROPERTY(Transient)
	TArray<TObjectPtr<UMaterialInstanceDynamic>> Materials;

	TFuture<FGearEvaluation> Evaluation;
	FGearEvaluation Finished;
	mkb::MechanismPlayer Player;
	bool bComputing = false;
	bool bPlayingBack = false;
	double SpinAngle = 0.0;            // the crank while the machine computes, before playback
	double LastCrankTurns = 0.0;

	mkb::CalculatorDevice Device;
	mkb::Backlight Lights;
	FKeyboardAudio Audio;
	double Clock = 0.0;
	float Explode = 0.0f;
	float ExplodeVelocity = 0.0f;
	bool bExplodeTarget = false;
	bool bCursorOn = true;
	double CursorTimer = 0.0;
	TArray<FString> LastScreen;
	TArray<FDemoStep> Demo;
	double DemoTimer = 0.0;
	int32 DemoLine = 0;
	bool bDemoActive = false;
	int32 LastBacklightMode = -1;
};
