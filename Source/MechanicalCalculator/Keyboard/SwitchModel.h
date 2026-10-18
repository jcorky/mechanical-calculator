// Key travel dynamics for one MX-style switch: a stem of small effective mass driven by the finger,
// pushed back by a preloaded spring, with a tactile bump or click jacket, hard stops at both ends,
// and actuation hysteresis. Positions in millimetres of travel (0 = up), time in seconds.
// Engine-independent; the Unreal layer steps one of these per key every frame.
#pragma once

namespace mkb
{
enum class SwitchFeel
{
	Linear,    // e.g. reds: smooth 45 gf
	Tactile,   // e.g. browns: a bump before actuation
	Clicky,    // e.g. blues: bump plus a click jacket that snaps
};

struct SwitchSpec
{
	SwitchFeel Feel = SwitchFeel::Tactile;
	double TravelMm = 4.0;
	double ActuationMm = 2.0;
	double HysteresisMm = 0.3;       // resets this far above the actuation point
	double PreloadN = 0.35;          // spring force with the key up
	double BottomForceN = 0.62;      // spring force at full travel
	double BumpN = 0.20;             // extra force of the tactile bump (0 for linear)
	double BumpAtMm = 1.4;
	double BumpWidthMm = 0.45;
	double MassKg = 0.0035;          // stem + keycap (+ the finger's share while pressing)
	double FrictionN = 0.02;
	double FingerN = 0.95;           // a typical keystroke force
	double BottomRestitution = 0.18;
	double TopRestitution = 0.30;

	static SwitchSpec For(SwitchFeel Feel);
	double SpringForce(double X) const;   // spring plus bump, N, at travel X
};

struct SwitchEvents
{
	bool Actuated = false;     // crossed the actuation point going down
	bool Reset = false;        // crossed back above the reset point
	bool BottomedOut = false;  // hit the bottom stop (strength in BottomSpeed)
	bool ToppedOut = false;    // returned against the top stop
	bool Clicked = false;      // click jacket snapped (clicky only)
	double BottomSpeed = 0.0;  // m/s at the bottom impact
	double TopSpeed = 0.0;
};

class SwitchModel
{
public:
	explicit SwitchModel(const SwitchSpec& InSpec = SwitchSpec()) : Spec(InSpec) {}

	void SetPressed(bool bPressed) { Pressed = bPressed; }
	bool IsPressed() const { return Pressed; }
	// Advances by Dt (any size; sub-stepped internally) and reports what happened.
	SwitchEvents Step(double Dt);

	double PositionMm() const { return X * 1000.0; }
	double Fraction() const { return X * 1000.0 / Spec.TravelMm; }
	bool IsActuated() const { return Actuated; }
	bool IsAtRest() const { return !Pressed && X <= 0.0 && V == 0.0; }
	const SwitchSpec& GetSpec() const { return Spec; }

private:
	SwitchSpec Spec;
	bool Pressed = false;
	bool Actuated = false;
	bool ClickArmed = true;
	double X = 0.0;   // m
	double V = 0.0;   // m/s, positive downward
};
} // namespace mkb
