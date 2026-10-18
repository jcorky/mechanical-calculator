// Plays the synthesized key sounds (Keyboard/KeySound.h). A small bank of variants per sound and
// key size is synthesized once; each keystroke queues one into a fresh USoundWaveProcedural.
#pragma once

#include "CoreMinimal.h"
#include "KeySound.h"

class AActor;

class FKeyboardAudio
{
public:
	void Initialize(int32 InSampleRate = 48000);
	void Play(AActor* Owner, const FVector& Location, mkb::KeySoundKind Kind, double KeyWidthU, double Strength, float Volume);

private:
	struct FBank
	{
		mkb::KeySoundKind Kind;
		int32 SizeClass;      // 0: 1-1.75u, 1: 2-2.75u, 2: spacebar
		bool bStabilized;
		TArray<TArray<uint8>> Variants;
	};
	TArray<FBank> Banks;
	int32 SampleRate = 48000;
	int32 Counter = 0;

	static int32 SizeClassFor(double WidthU) { return WidthU >= 6.0 ? 2 : (WidthU >= 2.0 ? 1 : 0); }
	const FBank* Find(mkb::KeySoundKind Kind, int32 SizeClass) const;
};
