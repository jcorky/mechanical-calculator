// Turns the stepped-drum machine's recorded events (Calc/GearMachine.h) into continuous motion for
// the Unreal layer: crank and drum angles, sliding pinions, numeral-wheel positions with their
// staggered tens-carries, the carriage, the reversing lever and the bell. Engine-independent.
//
// One crank turn, as a fraction u of the revolution:
//   0.00 - 0.10   free travel
//   0.10 - 0.55   the drums' stepped teeth pass the pinions: wheel i + shift advances by digit i
//   0.60 - 0.95   the carry levers trip in turn, units wheel first, each adding one step
//   0.95 - 1.00   the crank returns to its rest stop
// Wheel positions are in digits (continuous; digit d faces the window at d * 36 degrees).
#pragma once

#include "GearMachine.h"

#include <cstdint>
#include <string>
#include <vector>

namespace mkb
{
struct MechanismState
{
	static constexpr int SettingDigits = mcalc::GearEvent::SettingDigits;
	static constexpr int AccumulatorDigits = mcalc::GearEvent::AccumulatorDigits;
	static constexpr int CounterDigits = mcalc::GearEvent::CounterDigits;

	double CrankAngle = 0.0;            // radians, accumulating; the drums turn with it
	double Carriage = 0.0;              // shift position, continuous 0..20
	double Setting[SettingDigits] = {}; // pinion positions along their drums, in digits (continuous)
	double Accumulator[AccumulatorDigits] = {};
	double Counter[CounterDigits] = {};
	double CarryLever[AccumulatorDigits] = {};   // 0 rest .. 1 tripped
	double Reverse = 0.0;               // reversing lever: 0 add .. 1 subtract
	double BellHammer = 0.0;            // 0 rest .. 1 struck, decays
	char Op = ' ';                      // procedure under way: + - * / r, or ' ' idle
	bool Running = false;
};

class MechanismPlayer
{
public:
	// Motor speed. The player compresses time so that no evaluation plays longer than MaxSeconds.
	double TurnsPerSecond = 12.0;
	double ShiftsPerSecond = 10.0;
	double SetSeconds = 0.12;
	double ClearSeconds = 0.35;
	double MaxSeconds = 6.0;

	// Queues an evaluation's events (it may be the first part of a longer computation: the
	// recorder keeps a bounded number; the result is printed either way).
	void Play(const std::vector<mcalc::GearEvent>& InEvents);
	// Shows a number in the setting register while idle (the digits being typed).
	void ShowSetting(const std::vector<uint8_t>& Digits);
	// Advances the animation; returns the number of bell strikes that happened.
	int Step(double Dt);
	void SkipToEnd();
	bool IsPlaying() const { return Index < Events.size(); }
	double Progress() const;            // 0..1 through the queued events
	const MechanismState& State() const { return S; }
	double Duration() const { return TotalSeconds; }

private:
	std::vector<mcalc::GearEvent> Events;
	size_t Index = 0;
	double Phase = 0.0;                 // seconds into the current event
	double TimeScale = 1.0;
	double TotalSeconds = 0.0;
	MechanismState S;
	mcalc::GearEvent Last;              // register snapshot before the current event
	double CrankBase = 0.0;
	double IdleSetting[MechanismState::SettingDigits] = {};

	double EventSeconds(const mcalc::GearEvent& E, const mcalc::GearEvent& Before) const;
	void Apply(const mcalc::GearEvent& E, const mcalc::GearEvent& Before, double U);
	void Settle(const mcalc::GearEvent& E);
};

// The number under the cursor in an input line, as setting-register digits (units first), so the
// setting indicators follow the keys while a number is typed. Empty if the cursor is not in a number.
std::vector<uint8_t> SettingFromInput(const std::string& Line, int Cursor);
} // namespace mkb
