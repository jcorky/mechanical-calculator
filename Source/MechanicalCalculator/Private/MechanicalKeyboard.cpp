#include "MechanicalKeyboard.h"

#include "Async/Async.h"
#include "Components/BoxComponent.h"
#include "Components/TextRenderComponent.h"
#include "KeyboardMaterials.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "MeshConversion.h"
#include "ProceduralMeshComponent.h"

namespace MechanicalKeyboardDetail
{
const FName BuiltTag(TEXT("MechanicalCalculatorBuilt"));

FColor ToColor(const FLinearColor& C)
{
	return C.ToFColor(true);
}

bool IsModifier(mkb::KeyAction Action)
{
	return Action == mkb::KeyAction::Shift || Action == mkb::KeyAction::Control || Action == mkb::KeyAction::Fn || Action == mkb::KeyAction::Alt ||
		Action == mkb::KeyAction::Super;
}

bool Repeats(mkb::KeyAction Action)
{
	switch (Action)
	{
	case mkb::KeyAction::Insert:
	case mkb::KeyAction::Backspace:
	case mkb::KeyAction::DeleteForward:
	case mkb::KeyAction::CursorLeft:
	case mkb::KeyAction::CursorRight:
	case mkb::KeyAction::HistoryUp:
	case mkb::KeyAction::HistoryDown:
	case mkb::KeyAction::PageUp:
	case mkb::KeyAction::PageDown:
		return true;
	default:
		return false;
	}
}

bool IsNumpadOperator(const mkb::KeyDef& Def)
{
	return Def.Id == "NumDivide" || Def.Id == "NumMultiply" || Def.Id == "NumSubtract" || Def.Id == "NumAdd";
}

// The engine fonts cover ASCII; the layout's arrow legends are UTF-8 arrows.
FString LegendText(const std::string& S)
{
	if (S == "\xE2\x86\x91") return TEXT("^");
	if (S == "\xE2\x86\x93") return TEXT("v");
	if (S == "\xE2\x86\x90") return TEXT("<");
	if (S == "\xE2\x86\x92") return TEXT(">");
	return FString(UTF8_TO_TCHAR(S.c_str()));
}

FString CacheKeyFor(const mgeo::KeycapParams& P)
{
	return FString::Printf(TEXT("%.3f_%.3f_%d_%d_%d"), P.WidthU, P.DepthU, P.ProfileRow, P.Spacebar ? 1 : 0, P.HomingBar ? 1 : 0);
}

// A frame on a keycap face, in the cap's local geometry frame (mm): origin, outward normal, and the
// in-plane direction that is "up" for text read by the typist.
struct FFaceFrame
{
	mgeo::Vec3 Origin;
	mgeo::Vec3 Normal;
	mgeo::Vec3 Up;
};

FFaceFrame TopFrame(const mgeo::KeycapTop& T)
{
	const double S = FMath::Sin(T.TiltRad), C = FMath::Cos(T.TiltRad);
	return FFaceFrame{mgeo::Vec3(0.0, T.CenterY, T.CenterZ), mgeo::Vec3(0.0, -S, C), mgeo::Vec3(0.0, C, S)};
}

FFaceFrame FrontFrame(const mgeo::KeycapTop& T)
{
	const double S = FMath::Sin(T.FrontTiltRad), C = FMath::Cos(T.FrontTiltRad);
	return FFaceFrame{mgeo::Vec3(0.0, T.FrontFaceY, T.FrontFaceZ), mgeo::Vec3(0.0, -C, S), mgeo::Vec3(0.0, S, C)};
}
} // namespace MechanicalKeyboardDetail

using namespace MechanicalKeyboardDetail;

AMechanicalKeyboard::AMechanicalKeyboard()
{
	PrimaryActorTick.bCanEverTick = true;
	PrimaryActorTick.TickGroup = TG_PrePhysics;
	Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	SetRootComponent(Root);
	BoardRoot = CreateDefaultSubobject<USceneComponent>(TEXT("BoardRoot"));
	BoardRoot->SetupAttachment(Root);

	// One line from each corner of Where It Isn't, typed on the keys until the first real keystroke.
	DemoLines = {
		TEXT("hohmann(6678e3, 42164e3)"),
		TEXT("db(snr_track(1.4e6, 2000, 0.1, 1, 111e3, 950, 1.67e6, 6.3))"),
		TEXT("erlangc(3, 2.4)"),
		TEXT("quat_rot(quat([0, 0, 1], pi/2), [1, 0, 0])"),
		TEXT("planform().S"),
		TEXT("A = [1, 0.1; 0, 1]; B = [0.005; 0.1]; dare(A, B, eye(2), 0.5)"),
		TEXT("help"),
	};
}

mkb::SwitchSpec AMechanicalKeyboard::Spec() const
{
	switch (SwitchType)
	{
	case EKeyboardSwitch::Linear: return mkb::SwitchSpec::For(mkb::SwitchFeel::Linear);
	case EKeyboardSwitch::Clicky: return mkb::SwitchSpec::For(mkb::SwitchFeel::Clicky);
	default: return mkb::SwitchSpec::For(mkb::SwitchFeel::Tactile);
	}
}

void AMechanicalKeyboard::OnConstruction(const FTransform& Transform)
{
	Super::OnConstruction(Transform);
	Build();
}

void AMechanicalKeyboard::BeginPlay()
{
	Super::BeginPlay();
	if (!bBuilt)
	{
		Build();
	}
	Audio.Initialize();
	Device.SetAsync(true);
	Lights.SetMode(Device.BacklightMode());
	LastBacklightMode = -1;
	if (bPlayDemo && DemoLines.Num() > 0)
	{
		bDemoActive = true;
		DemoLine = 0;
		DemoTimer = 1.0;
		QueueDemo();
	}
	RefreshScreen();
}

void AMechanicalKeyboard::EndPlay(const EEndPlayReason::Type EndPlayReason)
{
	if (Evaluation.IsValid())
	{
		Evaluation.Wait();   // the worker holds a pointer to the session
	}
	Super::EndPlay(EndPlayReason);
}

void AMechanicalKeyboard::StartEvaluation()
{
	mcalc::Session* Session = &Device.GetSession();
	const std::string Expression = Device.PendingExpression();
	const size_t Keep = static_cast<size_t>(FMath::Max(100, RecordedEvents));
	bComputing = true;
	Evaluation = Async(EAsyncExecution::Thread, [Session, Expression, Keep]() {
		FGearEvaluation R;
		mcalc::GearRecorder Recorder(Keep);
		mcalc::GearContext::Attach(&Recorder);
		mcalc::GearContext::ResetStats();
		R.Outcome = Session->Execute(Expression);
		R.Work = mcalc::GearContext::Stats();
		mcalc::GearContext::Attach(nullptr);
		R.Events = Recorder.GetEvents();
		return R;
	});
}

void AMechanicalKeyboard::TickMachine(float DeltaSeconds)
{
	if (bComputing && Evaluation.IsValid() && Evaluation.IsReady())
	{
		Finished = Evaluation.Get();
		Evaluation.Reset();
		bComputing = false;
		Player.TurnsPerSecond = MotorTurnsPerSecond;
		Player.MaxSeconds = MaxPlaybackSeconds;
		Player.Play(Finished.Events);
		bPlayingBack = true;
	}
	if (bComputing)
	{
		SpinAngle += DeltaSeconds * 2.0 * UE_DOUBLE_PI * MotorTurnsPerSecond;
	}
	const int32 Bells = Player.Step(DeltaSeconds);
	const FVector Here = GetActorLocation();
	for (int32 i = 0; i < FMath::Min(Bells, 2); ++i)
	{
		Audio.Play(this, Here, mkb::KeySoundKind::Bell, 1.0, 1.0, SoundVolume * 0.6f);
	}
	// A clack each time the crank passes its rest stop.
	const double Turns = FMath::FloorToDouble((Player.State().CrankAngle + SpinAngle) / (2.0 * UE_DOUBLE_PI));
	if (Turns > LastCrankTurns)
	{
		Audio.Play(this, Here, mkb::KeySoundKind::GearClack, 1.0, 0.6, SoundVolume * 0.35f);
	}
	LastCrankTurns = Turns;
	if (bPlayingBack && !Player.IsPlaying())
	{
		bPlayingBack = false;
		Device.CompleteEvaluation(Finished.Outcome, Finished.Work);
		RefreshScreen();
	}
	if (!Device.Busy())
	{
		Player.ShowSetting(mkb::SettingFromInput(Device.Input(), Device.Cursor()));
	}
}

void AMechanicalKeyboard::DestroyBuilt()
{
	TArray<UActorComponent*> Existing;
	GetComponents(Existing);
	for (UActorComponent* Component : Existing)
	{
		if (Component && Component->ComponentHasTag(BuiltTag))
		{
			Component->DestroyComponent();
		}
	}
	Keys.Reset();
	Layers.Reset();
	ScreenLines.Reset();
	Materials.Reset();
	Tops = nullptr;
	TopColors.Reset();
	bBuilt = false;
}

UProceduralMeshComponent* AMechanicalKeyboard::NewMesh(USceneComponent* Parent, const FName& Name)
{
	UProceduralMeshComponent* Mesh = NewObject<UProceduralMeshComponent>(this, MakeUniqueObjectName(this, UProceduralMeshComponent::StaticClass(), Name));
	Mesh->CreationMethod = HasActorBegunPlay() ? EComponentCreationMethod::Instance : EComponentCreationMethod::UserConstructionScript;
	Mesh->ComponentTags.Add(BuiltTag);
	Mesh->SetupAttachment(Parent);
	Mesh->SetCollisionEnabled(ECollisionEnabled::NoCollision);
	Mesh->RegisterComponent();
	return Mesh;
}

UTextRenderComponent* AMechanicalKeyboard::NewText(USceneComponent* Parent, const FName& Name)
{
	UTextRenderComponent* Text = NewObject<UTextRenderComponent>(this, MakeUniqueObjectName(this, UTextRenderComponent::StaticClass(), Name));
	Text->CreationMethod = HasActorBegunPlay() ? EComponentCreationMethod::Instance : EComponentCreationMethod::UserConstructionScript;
	Text->ComponentTags.Add(BuiltTag);
	Text->SetupAttachment(Parent);
	Text->SetCollisionEnabled(ECollisionEnabled::NoCollision);
	Text->SetCastShadow(false);
	Text->SetHorizontalAlignment(EHTA_Center);
	Text->SetVerticalAlignment(EVRTA_TextCenter);
	Text->RegisterComponent();
	return Text;
}

UMaterialInstanceDynamic* AMechanicalKeyboard::Surface(const FLinearColor& Color, float Metallic, float Roughness, float Specular, bool bTranslucent, float Opacity,
	const FLinearColor& Emissive, float VertexEmissive)
{
	FKeyboardSurface S;
	S.BaseColor = Color;
	S.Metallic = Metallic;
	S.Roughness = Roughness;
	S.Specular = Specular;
	S.bTranslucent = bTranslucent;
	S.Opacity = Opacity;
	S.Emissive = Emissive;
	S.VertexEmissive = VertexEmissive;
	UMaterialInstanceDynamic* Mid = KeyboardMaterials::Make(this, S);
	if (Mid)
	{
		Materials.Add(Mid);
	}
	return Mid;
}

FLinearColor AMechanicalKeyboard::CapColor(const mkb::KeyDef& Def) const
{
	if (IsNumpadOperator(Def))
	{
		return Theme.NumpadOperatorCaps;
	}
	switch (Def.Role)
	{
	case mkb::KeyRole::Modifier: return Theme.ModifierCaps;
	case mkb::KeyRole::Accent: return Theme.AccentCaps;
	case mkb::KeyRole::Function:
	case mkb::KeyRole::Navigation: return Theme.FunctionCaps;
	case mkb::KeyRole::Numpad: return Theme.NumpadCaps;
	default: return Theme.AlphaCaps;
	}
}

bool AMechanicalKeyboard::CapIsDark(const mkb::KeyDef& Def) const
{
	return IsNumpadOperator(Def) || Def.Role == mkb::KeyRole::Modifier || Def.Role == mkb::KeyRole::Accent || Def.Role == mkb::KeyRole::Function ||
		Def.Role == mkb::KeyRole::Navigation;
}

void AMechanicalKeyboard::AddLegends(FKeyRuntime& Runtime, const mkb::KeyDef& Def, const mgeo::KeycapTop& Top)
{
	if (!bShowLegends)
	{
		return;
	}
	const bool bDark = CapIsDark(Def);
	const FColor Ink = ToColor(bDark ? Theme.LegendOnDark : Theme.LegendOnLight);
	auto Place = [this, &Runtime](const FFaceFrame& F, double Right, double Up, double Lift, double HeightMm, const FString& Text, const FColor& Color, EHorizTextAligment Align) {
		UTextRenderComponent* T = NewText(Runtime.Cap, FName(TEXT("Legend")));
		const mgeo::Vec3 P = F.Origin + mgeo::Vec3(Right, 0.0, 0.0) + F.Up * Up + F.Normal * Lift;
		T->SetRelativeLocation(KeyboardMesh::ToUnreal(P));
		T->SetRelativeRotation(FRotationMatrix::MakeFromXZ(KeyboardMesh::DirToUnreal(F.Normal), KeyboardMesh::DirToUnreal(F.Up)).Rotator());
		T->SetHorizontalAlignment(Align);
		T->SetWorldSize(static_cast<float>(HeightMm * 0.1));
		T->SetTextRenderColor(Color);
		T->SetText(FText::FromString(Text));
	};
	const FFaceFrame TopF = TopFrame(Top);
	const double HalfDepth = FMath::Max(Top.Depth * 0.5, 1.0);
	// The dish is lowest at its centre; legends away from it sit a little higher.
	auto LiftAt = [HalfDepth](double Up) { return 0.06 + 0.6 * FMath::Min(1.0, FMath::Square(Up / HalfDepth)); };
	const FString Main = LegendText(Def.Legend);
	if (!Def.ShiftLegend.empty())
	{
		const double Off = Top.Depth * 0.22;
		Place(TopF, 0.0, Off, LiftAt(Off), 3.0, LegendText(Def.ShiftLegend), Ink, EHTA_Center);
		Place(TopF, 0.0, -Off, LiftAt(-Off), 3.0, Main, Ink, EHTA_Center);
	}
	else if (!Main.IsEmpty())
	{
		const bool bWord = Main.Len() > 2;
		const double Size = bWord ? (Main.Len() > 5 ? 2.0 : 2.4) : 3.6;
		Place(TopF, 0.0, 0.0, LiftAt(0.0), Size, Main, Ink, EHTA_Center);
	}
	if (!Def.FrontLegend.empty())
	{
		const FColor FrontInk = ToColor(Def.Role == mkb::KeyRole::Accent || IsNumpadOperator(Def) ? Theme.LegendOnDark : Theme.FrontLegend);
		Place(FrontFrame(Top), 0.0, 0.0, 0.05, 1.7, LegendText(Def.FrontLegend), FrontInk, EHTA_Center);
	}
}

void AMechanicalKeyboard::Build()
{
	DestroyBuilt();
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mgeo::CaseParams CaseP;
	const mgeo::BoardParts Board = mgeo::BuildBoard(L, CaseP);
	BoardRoot->SetRelativeRotation(FRotator(FMath::RadiansToDegrees(Board.RestTiltRad), 0.0, 0.0));
	BoardRoot->SetRelativeLocation(FVector(0.0, 0.0, Board.RestLiftMm * 0.1));

	const bool bFlip = bFlipWinding;
	auto Section = [bFlip](UProceduralMeshComponent* Mesh, int32 Index, const mgeo::Mesh& Src, UMaterialInterface* Material) {
		if (Src.Triangles.empty())
		{
			return;
		}
		KeyboardMesh::SetSection(Mesh, Index, KeyboardMesh::Convert(Src, FLinearColor::Black, bFlip));
		Mesh->SetMaterial(Index, Material);
	};

	// Shared surfaces.
	UMaterialInstanceDynamic* Aluminium = Surface(Theme.Case, 1.0f, 0.38f);
	UMaterialInstanceDynamic* Brass = Surface(FLinearColor(0.93f, 0.72f, 0.38f), 1.0f, 0.3f);
	UMaterialInstanceDynamic* Rubber = Surface(FLinearColor(0.02f, 0.02f, 0.02f), 0.0f, 0.92f, 0.2f);
	UMaterialInstanceDynamic* PlateMat = Surface(Theme.Plate, 1.0f, 0.32f);
	UMaterialInstanceDynamic* PcbMat = Surface(FLinearColor(0.015f, 0.17f, 0.055f), 0.0f, 0.45f);
	UMaterialInstanceDynamic* Chips = Surface(FLinearColor(0.02f, 0.02f, 0.022f), 0.2f, 0.5f);
	UMaterialInstanceDynamic* Nylon = Surface(FLinearColor(0.025f, 0.025f, 0.028f), 0.0f, 0.62f);
	UMaterialInstanceDynamic* Gold = Surface(FLinearColor(1.0f, 0.78f, 0.38f), 1.0f, 0.22f);
	UMaterialInstanceDynamic* Steel = Surface(FLinearColor(0.82f, 0.83f, 0.85f), 1.0f, 0.22f);
	UMaterialInstanceDynamic* Housing = Surface(FLinearColor(0.85f, 0.87f, 0.9f), 0.0f, 0.12f, 0.6f, true, 0.38f, FLinearColor::Black, 9.0f);
	UMaterialInstanceDynamic* Bezel = Surface(FLinearColor(0.012f, 0.012f, 0.014f), 0.6f, 0.4f);
	UMaterialInstanceDynamic* Glass = Surface(FLinearColor(0.004f, 0.006f, 0.008f), 0.0f, 0.04f, 0.8f, false, 1.0f, FLinearColor(0.0f, 0.004f, 0.006f));
	const FLinearColor StemColor = SwitchType == EKeyboardSwitch::Linear ? FLinearColor(0.7f, 0.04f, 0.04f)
		: SwitchType == EKeyboardSwitch::Clicky ? FLinearColor(0.08f, 0.22f, 0.75f) : FLinearColor(0.28f, 0.13f, 0.06f);
	UMaterialInstanceDynamic* StemMat = Surface(StemColor, 0.0f, 0.45f);
	UMaterialInstanceDynamic* StabStemMat = Surface(FLinearColor(0.9f, 0.9f, 0.88f), 0.0f, 0.5f);

	// Static layers, bottom to top. ExplodeFraction lifts each one in the exploded view.
	auto AddLayer = [this](const TCHAR* Name, float Fraction) {
		FLayer Layer;
		Layer.Mesh = NewMesh(BoardRoot, FName(Name));
		Layer.ExplodeFraction = Fraction;
		Layers.Add(Layer);
		return Layer.Mesh;
	};
	UProceduralMeshComponent* CaseLayer = AddLayer(TEXT("Case"), 0.0f);
	Section(CaseLayer, 0, Board.CaseShell, Aluminium);
	Section(CaseLayer, 1, Board.CaseBadge, Brass);
	Section(CaseLayer, 2, Board.Feet, Rubber);
	UProceduralMeshComponent* PcbLayer = AddLayer(TEXT("Pcb"), 0.12f);
	Section(PcbLayer, 0, Board.Pcb, PcbMat);
	Section(PcbLayer, 1, Board.PcbParts, Chips);

	mgeo::SwitchParams SwP;
	SwP.Type = SwitchType == EKeyboardSwitch::Linear ? mgeo::SwitchType::Linear : SwitchType == EKeyboardSwitch::Clicky ? mgeo::SwitchType::Clicky : mgeo::SwitchType::Tactile;
	SwP.TravelMm = Spec().TravelMm;
	SwP.ActuationMm = Spec().ActuationMm;
	const mgeo::SwitchParts Sw = mgeo::BuildSwitch(SwP);
	mgeo::Mesh Bottoms, Contacts, Pins, TopsAll, StabHousings, StabWires;
	const FMeshSectionData StemData = KeyboardMesh::Convert(Sw.Stem, FLinearColor::Black, bFlip);
	mgeo::Mesh SpringLocal = Sw.Spring;
	SpringLocal.Translate(mgeo::Vec3(0.0, 0.0, -Sw.SpringBottomZ));
	const FMeshSectionData SpringData = KeyboardMesh::Convert(SpringLocal, FLinearColor::Black, bFlip);
	SpringRestCm = FMath::Max(0.05, (Sw.SpringRestTopZ - Sw.SpringBottomZ) * 0.1);
	TMap<FString, FMeshSectionData> CapCache;
	const mkb::SwitchSpec KeySpec = Spec();

	for (int32 i = 0; i < static_cast<int32>(L.Keys.size()); ++i)
	{
		const mkb::KeyDef& Def = L.Keys[static_cast<size_t>(i)];
		const mgeo::Vec3 C = mgeo::KeyCenter(Def, L);
		Bottoms.Append(Sw.HousingBottom, C);
		Contacts.Append(Sw.Contacts, C);
		Pins.Append(Sw.Pins, C);

		FKeyRuntime& R = Keys.AddDefaulted_GetRef();
		R.LayoutIndex = i;
		R.Switch = mkb::SwitchModel(KeySpec);
		R.HousingVertexStart = static_cast<int32>(TopsAll.Positions.size());
		TopsAll.Append(Sw.HousingTop, C);
		R.HousingVertexCount = static_cast<int32>(TopsAll.Positions.size()) - R.HousingVertexStart;

		const FString Base = FString(UTF8_TO_TCHAR(Def.Id.c_str()));
		R.Mover = NewObject<USceneComponent>(this, MakeUniqueObjectName(this, USceneComponent::StaticClass(), FName(*(TEXT("Key_") + Base))));
		R.Mover->CreationMethod = HasActorBegunPlay() ? EComponentCreationMethod::Instance : EComponentCreationMethod::UserConstructionScript;
		R.Mover->ComponentTags.Add(BuiltTag);
		R.Mover->SetupAttachment(BoardRoot);
		R.Mover->RegisterComponent();
		R.MoverBase = KeyboardMesh::ToUnreal(C + mgeo::Vec3(0.0, 0.0, SwP.StemTopZ));
		R.Mover->SetRelativeLocation(R.MoverBase);

		R.Stem = NewMesh(R.Mover, FName(*(TEXT("Stem_") + Base)));
		KeyboardMesh::SetSection(R.Stem, 0, StemData);
		R.Stem->SetMaterial(0, StemMat);
		if (Def.NeedsStabilizer())
		{
			const bool bVertical = Def.H > Def.W;
			const mgeo::StabilizerParts Stab = mgeo::BuildStabilizer(bVertical ? Def.H : Def.W, bVertical);
			StabHousings.Append(Stab.Housings, C);
			StabWires.Append(Stab.Wire, C);
			Section(R.Stem, 1, Stab.Stems, StabStemMat);
		}

		mgeo::KeycapParams CapP;
		CapP.WidthU = Def.W;
		CapP.DepthU = Def.H;
		CapP.ProfileRow = Def.ProfileRow;
		CapP.Spacebar = Def.Role == mkb::KeyRole::Space;
		CapP.HomingBar = Def.Id == "F" || Def.Id == "J" || Def.Id == "Num5";
		const FString CacheKey = CacheKeyFor(CapP);
		if (!CapCache.Contains(CacheKey))
		{
			CapCache.Add(CacheKey, KeyboardMesh::Convert(mgeo::BuildKeycap(CapP), FLinearColor::Black, bFlip));
		}
		R.Cap = NewMesh(R.Mover, FName(*(TEXT("Cap_") + Base)));
		KeyboardMesh::SetSection(R.Cap, 0, CapCache[CacheKey]);
		UMaterialInstanceDynamic* CapMat = Surface(CapColor(Def), 0.0f, 0.58f, 0.4f);
		R.Cap->SetMaterial(0, CapMat);

		const mgeo::KeycapTop Top = mgeo::KeycapTopFor(CapP);
		R.Hit = NewObject<UBoxComponent>(this, MakeUniqueObjectName(this, UBoxComponent::StaticClass(), FName(*(TEXT("Hit_") + Base))));
		R.Hit->CreationMethod = HasActorBegunPlay() ? EComponentCreationMethod::Instance : EComponentCreationMethod::UserConstructionScript;
		R.Hit->ComponentTags.Add(BuiltTag);
		R.Hit->SetupAttachment(R.Cap);
		const double HalfW = (Def.W * mkb::kUnitMm - 1.0) * 0.5, HalfD = (Def.H * mkb::kUnitMm - 1.0) * 0.5, HalfH = FMath::Max(Top.CenterZ, 4.0) * 0.5;
		R.Hit->SetBoxExtent(FVector(HalfD, HalfW, HalfH) * 0.1);
		R.Hit->SetRelativeLocation(FVector(0.0, 0.0, HalfH * 0.1));
		R.Hit->SetCollisionEnabled(ECollisionEnabled::QueryOnly);
		R.Hit->SetCollisionResponseToAllChannels(ECR_Ignore);
		R.Hit->SetCollisionResponseToChannel(ECC_Visibility, ECR_Block);
		R.Hit->SetHiddenInGame(true);
		R.Hit->RegisterComponent();

		R.SpringBase = KeyboardMesh::ToUnreal(C + mgeo::Vec3(0.0, 0.0, Sw.SpringBottomZ));
		R.Spring = NewMesh(BoardRoot, FName(*(TEXT("Spring_") + Base)));
		KeyboardMesh::SetSection(R.Spring, 0, SpringData);
		R.Spring->SetMaterial(0, Steel);
		R.Spring->SetRelativeLocation(R.SpringBase);
		R.Spring->SetCastShadow(false);

		AddLegends(R, Def, Top);
	}

	UProceduralMeshComponent* BottomLayer = AddLayer(TEXT("SwitchBottoms"), 0.26f);
	Section(BottomLayer, 0, Bottoms, Nylon);
	Section(BottomLayer, 1, Contacts, Gold);
	Section(BottomLayer, 2, Pins, Gold);
	UProceduralMeshComponent* PlateLayer = AddLayer(TEXT("Plate"), 0.46f);
	Section(PlateLayer, 0, Board.Plate, PlateMat);
	Section(PlateLayer, 1, StabHousings, Nylon);
	Section(PlateLayer, 2, StabWires, Steel);
	Section(PlateLayer, 3, Board.DisplayBezel, Bezel);
	Section(PlateLayer, 4, Board.DisplayGlass, Glass);
	Tops = AddLayer(TEXT("SwitchTops"), 0.56f);
	Section(Tops, 0, TopsAll, Housing);
	Tops->SetCastShadow(false);
	TopColors.Init(FLinearColor::Black, static_cast<int32>(TopsAll.Positions.size()));

	// OLED: one text line per character row, laid on the glass.
	if (Board.DisplayWidthMm > 10.0 && Board.DisplayHeightMm > 5.0)
	{
		const double Margin = 1.6;
		const double Tilt = Board.DisplayTiltRad;
		const mgeo::Vec3 N(0.0, -FMath::Sin(Tilt), FMath::Cos(Tilt));
		const mgeo::Vec3 Up(0.0, FMath::Cos(Tilt), FMath::Sin(Tilt));
		const double LineMm = (Board.DisplayHeightMm - 2.0 * Margin) / ScreenRows;
		const double GlyphMm = LineMm * 0.8;
		ScreenCols = FMath::Clamp(static_cast<int32>((Board.DisplayWidthMm - 2.0 * Margin) / (0.55 * GlyphMm)), 20, 48);
		for (int32 Row = 0; Row < ScreenRows; ++Row)
		{
			UTextRenderComponent* T = NewText(PlateLayer, FName(TEXT("ScreenLine")));
			const double V = Board.DisplayHeightMm * 0.5 - Margin - (Row + 0.5) * LineMm;
			const mgeo::Vec3 P = Board.DisplayCenter + mgeo::Vec3(-Board.DisplayWidthMm * 0.5 + Margin, 0.0, 0.0) + Up * V + N * 0.08;
			T->SetRelativeLocation(KeyboardMesh::ToUnreal(P));
			T->SetRelativeRotation(FRotationMatrix::MakeFromXZ(KeyboardMesh::DirToUnreal(N), KeyboardMesh::DirToUnreal(Up)).Rotator());
			T->SetHorizontalAlignment(EHTA_Left);
			T->SetWorldSize(static_cast<float>(GlyphMm * 0.1));
			T->SetTextRenderColor(ToColor(Theme.Screen));
			T->SetText(FText::GetEmpty());
			ScreenLines.Add(T);
		}
	}
	LastScreen.Reset();
	LastBacklightMode = -1;
	bBuilt = true;
	RefreshScreen();
}

int32 AMechanicalKeyboard::IndexOfKey(const FString& Id) const
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const std::string Wanted(TCHAR_TO_UTF8(*Id));
	for (int32 i = 0; i < Keys.Num(); ++i)
	{
		if (L.Keys[static_cast<size_t>(Keys[i].LayoutIndex)].Id == Wanted)
		{
			return i;
		}
	}
	return INDEX_NONE;
}

int32 AMechanicalKeyboard::IndexOfEngineKey(const FKey& Key) const
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mkb::KeyDef* Def = L.FindByEngineKey(std::string(TCHAR_TO_UTF8(*Key.GetFName().ToString())));
	if (!Def)
	{
		return INDEX_NONE;
	}
	return IndexOfKey(FString(UTF8_TO_TCHAR(Def->Id.c_str())));
}

void AMechanicalKeyboard::HandleEvents(const mkb::DeviceEvents& Events)
{
	if (Events.ExplodeToggled)
	{
		bExplodeTarget = Device.Exploded();
	}
	if (Events.BacklightCycled)
	{
		Lights.SetMode(Device.BacklightMode());
	}
	if (Events.EvaluationStarted)
	{
		StartEvaluation();
	}
	if (Events.ScreenChanged)
	{
		bCursorOn = true;
		CursorTimer = 0.0;
		RefreshScreen();
	}
}

void AMechanicalKeyboard::LogicalPress(int32 Index)
{
	if (!Keys.IsValidIndex(Index))
	{
		return;
	}
	const mkb::KeyDef& Def = mkb::FullSizeAnsi().Keys[static_cast<size_t>(Keys[Index].LayoutIndex)];
	Lights.OnKeyPressed(Def, Clock);
	HandleEvents(Device.Press(Def));
}

void AMechanicalKeyboard::LogicalRelease(int32 Index)
{
	if (!Keys.IsValidIndex(Index))
	{
		return;
	}
	Device.Release(mkb::FullSizeAnsi().Keys[static_cast<size_t>(Keys[Index].LayoutIndex)]);
	RefreshScreen();
}

void AMechanicalKeyboard::VisualPress(int32 Index)
{
	if (!Keys.IsValidIndex(Index))
	{
		return;
	}
	FKeyRuntime& R = Keys[Index];
	R.Switch.SetPressed(true);
	R.PressedAt = Clock;
	R.bReleasePending = false;
}

void AMechanicalKeyboard::VisualRelease(int32 Index)
{
	if (Keys.IsValidIndex(Index))
	{
		Keys[Index].bReleasePending = true;
	}
}

void AMechanicalKeyboard::ReleaseLatched()
{
	for (int32 i = 0; i < Keys.Num(); ++i)
	{
		if (Keys[i].bLatched)
		{
			Keys[i].bLatched = false;
			LogicalRelease(i);
			VisualRelease(i);
		}
	}
}

void AMechanicalKeyboard::PhysicalKeyDown(const FKey& Key)
{
	const int32 Index = IndexOfEngineKey(Key);
	if (Index == INDEX_NONE)
	{
		return;
	}
	CancelDemo();
	VisualPress(Index);
	if (Device.Busy() && bPlayingBack)
	{
		Player.SkipToEnd();   // any key hurries the mechanism to its result
	}
	LogicalPress(Index);
}

void AMechanicalKeyboard::PhysicalKeyUp(const FKey& Key)
{
	const int32 Index = IndexOfEngineKey(Key);
	if (Index == INDEX_NONE)
	{
		return;
	}
	LogicalRelease(Index);
	VisualRelease(Index);
}

void AMechanicalKeyboard::PhysicalKeyRepeat(const FKey& Key)
{
	const int32 Index = IndexOfEngineKey(Key);
	if (Index != INDEX_NONE && Repeats(mkb::FullSizeAnsi().Keys[static_cast<size_t>(Keys[Index].LayoutIndex)].Action))
	{
		LogicalPress(Index);
	}
}

bool AMechanicalKeyboard::ClickComponent(const UPrimitiveComponent* Component)
{
	for (int32 i = 0; i < Keys.Num(); ++i)
	{
		FKeyRuntime& R = Keys[i];
		if (Component != R.Hit && Component != R.Cap && Component != R.Stem)
		{
			continue;
		}
		CancelDemo();
		const mkb::KeyDef& Def = mkb::FullSizeAnsi().Keys[static_cast<size_t>(R.LayoutIndex)];
		if (IsModifier(Def.Action))
		{
			// Sticky modifiers: a click holds it down until the next key.
			if (R.bLatched)
			{
				R.bLatched = false;
				LogicalRelease(i);
				VisualRelease(i);
			}
			else
			{
				R.bLatched = true;
				VisualPress(i);
				LogicalPress(i);
			}
			return true;
		}
		VisualPress(i);
		LogicalPress(i);
		LogicalRelease(i);
		VisualRelease(i);
		ReleaseLatched();
		return true;
	}
	return false;
}

void AMechanicalKeyboard::TapIndex(int32 Index)
{
	if (!Keys.IsValidIndex(Index))
	{
		return;
	}
	VisualPress(Index);
	LogicalPress(Index);
	LogicalRelease(Index);
	VisualRelease(Index);
}

void AMechanicalKeyboard::TapKey(const FString& Id)
{
	TapIndex(IndexOfKey(Id));
}

void AMechanicalKeyboard::TypeLine(const FString& Text, bool bEnter)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const std::string Utf8(TCHAR_TO_UTF8(*Text));
	double Delay = Demo.Num() == 0 ? 0.05 : 0.0;
	for (char Ch : Utf8)
	{
		bool bShift = false;
		if (const mkb::KeyDef* Def = L.KeyForChar(Ch, &bShift))
		{
			FDemoStep Step;
			Step.KeyIndex = IndexOfKey(FString(UTF8_TO_TCHAR(Def->Id.c_str())));
			Step.bShift = bShift;
			Step.Delay = Delay;
			Demo.Add(Step);
			// Human-ish rhythm: a little faster inside words, a pause after spaces and commas.
			Delay = (1.0 / DemoCharsPerSecond) * ((Ch == ' ' || Ch == ',') ? 1.6 : FMath::FRandRange(0.7f, 1.2f));
		}
	}
	if (bEnter)
	{
		FDemoStep Step;
		Step.KeyIndex = IndexOfKey(TEXT("Enter"));
		Step.Delay = Delay + 0.25;
		Demo.Add(Step);
	}
	bDemoActive = true;
}

void AMechanicalKeyboard::QueueDemo()
{
	if (!DemoLines.IsValidIndex(DemoLine))
	{
		bDemoActive = false;
		return;
	}
	TypeLine(DemoLines[DemoLine], true);
	if (Demo.Num() > 0)
	{
		Demo[0].Delay = FMath::Max(Demo[0].Delay, 1.2);
	}
	++DemoLine;
}

void AMechanicalKeyboard::CancelDemo()
{
	if (bDemoActive)
	{
		Demo.Reset();
		bDemoActive = false;
		DemoLine = DemoLines.Num();
	}
}

void AMechanicalKeyboard::TickDemo(float DeltaSeconds)
{
	if (!bDemoActive)
	{
		return;
	}
	if (Demo.Num() == 0)
	{
		QueueDemo();
		return;
	}
	DemoTimer += DeltaSeconds;
	while (Demo.Num() > 0 && DemoTimer >= Demo[0].Delay)
	{
		DemoTimer -= Demo[0].Delay;
		const FDemoStep Step = Demo[0];
		Demo.RemoveAt(0);
		const int32 Shift = IndexOfKey(TEXT("LShift"));
		if (Step.bShift)
		{
			VisualPress(Shift);
			LogicalPress(Shift);
		}
		TapIndex(Step.KeyIndex);
		if (Step.bShift)
		{
			LogicalRelease(Shift);
			VisualRelease(Shift);
		}
	}
}

FString AMechanicalKeyboard::Evaluate(const FString& Expression)
{
	const mcalc::Outcome O = Device.GetSession().Execute(std::string(TCHAR_TO_UTF8(*Expression)));
	return FString(UTF8_TO_TCHAR(O.Text.c_str()));
}

void AMechanicalKeyboard::SetExploded(bool bExploded)
{
	if (Device.Exploded() != bExploded)
	{
		HandleEvents(Device.Press(*mkb::FullSizeAnsi().Find("Pause")));
	}
}

TArray<FString> AMechanicalKeyboard::GetScreenLines() const
{
	TArray<FString> Out;
	for (const std::string& Line : Device.Screen(ScreenCols, ScreenRows, false))
	{
		Out.Add(FString(UTF8_TO_TCHAR(Line.c_str())));
	}
	return Out;
}

void AMechanicalKeyboard::RefreshScreen()
{
	if (ScreenLines.Num() == 0)
	{
		return;
	}
	const std::vector<std::string> Lines = Device.Screen(ScreenCols, ScreenRows, bCursorOn);
	LastScreen.SetNum(ScreenLines.Num());
	for (int32 i = 0; i < ScreenLines.Num(); ++i)
	{
		const FString Text = i < static_cast<int32>(Lines.size()) ? FString(UTF8_TO_TCHAR(Lines[static_cast<size_t>(i)].c_str())) : FString();
		if (Text != LastScreen[i] && ScreenLines[i])
		{
			ScreenLines[i]->SetText(FText::FromString(Text));
			LastScreen[i] = Text;
		}
	}
}

void AMechanicalKeyboard::TickKeys(float DeltaSeconds)
{
	const double Lift = Explode * ExplodeHeight;
	const bool bLayersMoving = ExplodeVelocity != 0.0f || Explode != (bExplodeTarget ? 1.0f : 0.0f);
	for (int32 i = 0; i < Keys.Num(); ++i)
	{
		FKeyRuntime& R = Keys[i];
		if (R.bReleasePending && !R.bLatched && Clock - R.PressedAt >= MinPressTime)
		{
			R.Switch.SetPressed(false);
			R.bReleasePending = false;
		}
		const mkb::SwitchEvents Ev = R.Switch.Step(DeltaSeconds);
		const mkb::KeyDef& Def = mkb::FullSizeAnsi().Keys[static_cast<size_t>(R.LayoutIndex)];
		if (R.Mover && (Ev.BottomedOut || Ev.ToppedOut || Ev.Clicked))
		{
			const FVector Where = R.Mover->GetComponentLocation();
			if (Ev.BottomedOut)
			{
				Audio.Play(this, Where, mkb::KeySoundKind::BottomOut, Def.W, Ev.BottomSpeed / 0.5, SoundVolume);
			}
			if (Ev.ToppedOut)
			{
				Audio.Play(this, Where, mkb::KeySoundKind::TopOut, Def.W, Ev.TopSpeed / 0.5, SoundVolume);
			}
			if (Ev.Clicked)
			{
				Audio.Play(this, Where, mkb::KeySoundKind::Click, Def.W, 1.0, SoundVolume);
			}
		}
		// A key at rest stays where it was last placed unless the exploded view is moving.
		if (R.bSettled && !bLayersMoving && R.Switch.IsAtRest())
		{
			continue;
		}
		R.bSettled = R.Switch.IsAtRest();
		const double Travel = R.Switch.PositionMm() * 0.1;
		if (R.Mover)
		{
			R.Mover->SetRelativeLocation(R.MoverBase + FVector(0.0, 0.0, StemsFraction * Lift - Travel));
		}
		if (R.Cap)
		{
			R.Cap->SetRelativeLocation(FVector(0.0, 0.0, (CapsFraction - StemsFraction) * Lift));
		}
		if (R.Spring)
		{
			R.Spring->SetRelativeLocation(R.SpringBase + FVector(0.0, 0.0, SpringsFraction * Lift));
			R.Spring->SetRelativeScale3D(FVector(1.0, 1.0, FMath::Max(0.2, (SpringRestCm - Travel) / SpringRestCm)));
		}
	}
}

void AMechanicalKeyboard::TickExplode(float DeltaSeconds)
{
	const float Target = bExplodeTarget ? 1.0f : 0.0f;
	const float W = 5.0f;   // critically damped
	ExplodeVelocity += (W * W * (Target - Explode) - 2.0f * W * ExplodeVelocity) * DeltaSeconds;
	Explode += ExplodeVelocity * DeltaSeconds;
	if (FMath::Abs(Target - Explode) < 1e-4f && FMath::Abs(ExplodeVelocity) < 1e-4f)
	{
		Explode = Target;
		ExplodeVelocity = 0.0f;
	}
	for (const FLayer& Layer : Layers)
	{
		if (Layer.Mesh)
		{
			Layer.Mesh->SetRelativeLocation(FVector(0.0, 0.0, Layer.ExplodeFraction * Explode * ExplodeHeight));
		}
	}
}

void AMechanicalKeyboard::TickBacklight()
{
	const int32 Mode = Lights.GetMode();
	const bool bAnimated = Mode == 2 || Mode == 3 || Mode == 4;
	if (!Tops || (!bAnimated && Mode == LastBacklightMode))
	{
		return;
	}
	LastBacklightMode = Mode;
	Lights.Prune(Clock);
	const mkb::Layout& L = mkb::FullSizeAnsi();
	for (const FKeyRuntime& R : Keys)
	{
		const mkb::Rgb C = Lights.ColorFor(L.Keys[static_cast<size_t>(R.LayoutIndex)], Clock);
		const FLinearColor Lin(static_cast<float>(C.R), static_cast<float>(C.G), static_cast<float>(C.B), 1.0f);
		for (int32 v = 0; v < R.HousingVertexCount; ++v)
		{
			const int32 Idx = R.HousingVertexStart + v;
			if (TopColors.IsValidIndex(Idx))
			{
				TopColors[Idx] = Lin;
			}
		}
	}
	static const TArray<FVector> NoVertices;
	static const TArray<FVector> NoNormals;
	static const TArray<FVector2D> NoUVs;
	static const TArray<FProcMeshTangent> NoTangents;
	Tops->UpdateMeshSection_LinearColor(0, NoVertices, NoNormals, NoUVs, TopColors, NoTangents);
}

void AMechanicalKeyboard::Tick(float DeltaSeconds)
{
	Super::Tick(DeltaSeconds);
	const float Dt = FMath::Min(DeltaSeconds, 0.05f);
	Clock += Dt;
	if (!Device.Busy())
	{
		TickDemo(Dt);
	}
	TickMachine(Dt);
	TickKeys(Dt);
	TickExplode(Dt);
	TickBacklight();
	CursorTimer += Dt;
	if (CursorTimer >= 0.5)
	{
		CursorTimer = 0.0;
		bCursorOn = !bCursorOn;
		RefreshScreen();
	}
}
