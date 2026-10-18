#include "KeyboardAudio.h"

#include "GameFramework/Actor.h"
#include "Kismet/GameplayStatics.h"
#include "Sound/SoundWaveProcedural.h"

void FKeyboardAudio::Initialize(int32 InSampleRate)
{
	SampleRate = InSampleRate;
	Banks.Reset();
	const mkb::KeySoundKind Kinds[] = {mkb::KeySoundKind::BottomOut, mkb::KeySoundKind::TopOut, mkb::KeySoundKind::Click, mkb::KeySoundKind::Bell, mkb::KeySoundKind::GearClack};
	const double Widths[] = {1.0, 2.25, 6.25};
	for (mkb::KeySoundKind Kind : Kinds)
	{
		for (int32 SizeClass = 0; SizeClass < 3; ++SizeClass)
		{
			FBank& Bank = Banks.AddDefaulted_GetRef();
			Bank.Kind = Kind;
			Bank.SizeClass = SizeClass;
			Bank.bStabilized = SizeClass > 0;
			for (uint32 Seed = 1; Seed <= 6; ++Seed)
			{
				mkb::KeySoundParams P;
				P.Kind = Kind;
				P.KeyWidthU = Widths[SizeClass];
				P.Stabilized = Bank.bStabilized && Kind == mkb::KeySoundKind::BottomOut;
				P.Seed = Seed * 7919u + static_cast<uint32>(SizeClass) * 104729u + static_cast<uint32>(Kind) * 15485863u;
				P.SampleRate = SampleRate;
				const std::vector<int16_t> Pcm = mkb::SynthesizeKeySound(P);
				TArray<uint8>& Bytes = Bank.Variants.AddDefaulted_GetRef();
				Bytes.SetNumUninitialized(static_cast<int32>(Pcm.size() * sizeof(int16_t)));
				FMemory::Memcpy(Bytes.GetData(), Pcm.data(), Bytes.Num());
			}
		}
	}
}

const FKeyboardAudio::FBank* FKeyboardAudio::Find(mkb::KeySoundKind Kind, int32 SizeClass) const
{
	for (const FBank& Bank : Banks)
	{
		if (Bank.Kind == Kind && Bank.SizeClass == SizeClass)
		{
			return &Bank;
		}
	}
	return nullptr;
}

void FKeyboardAudio::Play(AActor* Owner, const FVector& Location, mkb::KeySoundKind Kind, double KeyWidthU, double Strength, float Volume)
{
	if (!Owner || Volume <= 0.0f)
	{
		return;
	}
	const FBank* Bank = Find(Kind, SizeClassFor(KeyWidthU));
	if (!Bank || Bank->Variants.Num() == 0)
	{
		return;
	}
	const TArray<uint8>& Bytes = Bank->Variants[Counter++ % Bank->Variants.Num()];
	USoundWaveProcedural* Wave = NewObject<USoundWaveProcedural>(Owner);
	Wave->SetSampleRate(SampleRate);
	Wave->NumChannels = 1;
	Wave->Duration = static_cast<float>(Bytes.Num() / 2) / static_cast<float>(SampleRate);
	Wave->SoundGroup = SOUNDGROUP_Default;
	Wave->bLooping = false;
	Wave->QueueAudio(Bytes.GetData(), Bytes.Num());
	const float Loudness = Volume * FMath::Clamp(static_cast<float>(0.55 + 0.45 * Strength), 0.2f, 1.2f);
	const float Pitch = FMath::FRandRange(0.97f, 1.03f);
	UGameplayStatics::PlaySoundAtLocation(Owner, Wave, Location, Loudness, Pitch);
}
