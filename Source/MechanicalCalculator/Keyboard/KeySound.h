// Synthesized key sounds: no audio assets ship with the project, so every press, click and release
// is generated as 16-bit mono PCM from a small physical caricature of a keystroke:
//   impact   band-passed noise with a ~3 ms decay (stem hitting the housing or the top stop)
//   body     damped modes of the case and plate ("thock"); lower and longer for wider keys
//   click    a short bright transient from the click jacket (clicky switches)
//   rattle   a few late, quiet impulses from a stabilizer wire on long keys
//   bell     inharmonic partials of a small gong, slowly decaying (the overdraft bell)
//   clack    a dull impact plus a ratchet tick (the crank's rest stop)
#pragma once

#include <cstdint>
#include <vector>

namespace mkb
{
enum class KeySoundKind
{
	BottomOut,   // the downstroke reaching the bottom
	TopOut,      // the release returning to the top
	Click,       // click jacket
	Bell,        // the machine's overdraft bell: a struck steel gong
	GearClack,   // the crank passing its rest stop at the end of a turn
};

struct KeySoundParams
{
	KeySoundKind Kind = KeySoundKind::BottomOut;
	double KeyWidthU = 1.0;
	double Strength = 1.0;        // impact speed relative to a normal keystroke
	bool Stabilized = false;
	uint32_t Seed = 1;
	int SampleRate = 48000;
};

std::vector<int16_t> SynthesizeKeySound(const KeySoundParams& P);
} // namespace mkb
