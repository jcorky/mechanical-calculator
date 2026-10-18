// Per-key RGB backlight effects (the colours shine through the translucent switch housings).
// Mode numbers match CalculatorDevice::BacklightMode(): 0 off, 1 static, 2 wave, 3 breathe,
// 4 reactive (ripples from each keystroke), 5 roles (keys coloured by what they do).
#pragma once

#include "KeyLayout.h"

#include <vector>

namespace mkb
{
struct Rgb
{
	double R = 0.0, G = 0.0, B = 0.0;   // linear, 0..1
};

class Backlight
{
public:
	void SetMode(int InMode) { Mode = InMode; }
	int GetMode() const { return Mode; }
	void OnKeyPressed(const KeyDef& Key, double Time);
	Rgb ColorFor(const KeyDef& Key, double Time) const;
	void Prune(double Time);   // forget ripples that have faded

private:
	struct Ripple
	{
		double X, Y, Start;
	};
	int Mode = 2;
	std::vector<Ripple> Ripples;
};

Rgb HsvToRgb(double H, double S, double V);   // H in turns (0..1)
} // namespace mkb
