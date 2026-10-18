// Poses every moving part of the calculating mechanism (Geometry/GearParts.h) from the player's
// state (MechanismPlayer.h), the memory columns and the printer. Engine-independent: the Unreal
// layer only turns the poses into component transforms.
//
// Shafts after the reverser (main shaft, drums, counting finger, the pinions and everything they
// drive) turn by the crank's motion signed by the reversing lever, so their angles are integrated
// from frame to frame. The figure wheels of the memory columns and the printer's type wheels spin
// forward to their targets at their own speed, and the paper feeds one line pitch per printed line.
#pragma once

#include "GearParts.h"
#include "MechanismPlayer.h"
#include "RegisterReadings.h"

#include <cstdint>
#include <vector>

namespace mkb
{
enum class MovingPart : int
{
	Drums,
	SquareShafts,
	Pinions,
	TransferGears,
	SettingWheels,
	SettingMarker,
	Carriage,
	AccumulatorWheels,
	CarryLevers,
	CounterWheels,
	CarriageMarkers,
	ClearingLevers,
	CountingFinger,
	MainShaft,
	ReverserIdler,
	InputShaft,
	DriveGear,
	Flywheel,
	Crank,
	CrankHandle,
	ReversingLever,
	ShiftKnob,
	ShiftPinion,
	BellHammer,
	TypeWheels,
	PaperRoll,
	Platen,
	RibbonSpools,
	MemoryWheels,
	MemorySignWheels,
	Count,
};
constexpr int kMovingParts = static_cast<int>(MovingPart::Count);

const char* MovingPartName(MovingPart P);
const mgeo::PartSet& PartSetOf(const mgeo::MechanismParts& M, MovingPart P);
// Parts posed in the carriage frame: they ride on the carriage (the carriage itself included, at
// the identity).
bool RidesCarriage(const mgeo::MechanismParts& M, MovingPart P);

class MechanismAnimator
{
public:
	double FigureWheelDigitsPerSecond = 16.0;
	double PaperLinesPerSecond = 9.0;

	// Keeps a reference to the parts, which must outlive the animator.
	explicit MechanismAnimator(const mgeo::MechanismParts& InParts);

	// Follows the player. ExtraCrankRad is added to the player's crank angle: the motor runs while
	// the machine computes, before the record plays. Both only ever grow.
	void Follow(const MechanismState& S, double ExtraCrankRad);
	// Targets for the memory columns (left to right) and the type wheels.
	void ShowMemory(const std::vector<MemoryColumn>& Columns);
	void ShowTypeWheels(const WheelReading& R);
	// Feeds the paper by Lines more line pitches.
	void Feed(int Lines);
	// Decimal markers, in places from the units digit.
	void SetDecimalMarkers(int SettingPlaces, int AccumulatorPlaces, int CounterPlaces);
	// Spins the figure wheels and feeds the paper toward their targets.
	void Step(double Dt);
	// The figure wheels, type wheels, markers and paper are at their targets.
	bool Settled() const;

	const std::vector<mgeo::Pose>& Poses(MovingPart P) const { return Current[static_cast<int>(P)]; }
	// The carriage frame in the desk frame (Carriage.Poses[0] moved by the shift).
	const mgeo::Pose& CarriageFrame() const { return CarriageAt; }
	// Line N of the paper (0 = the print line), following the feed: a line printed but not yet fed
	// sits behind the print line until the paper moves.
	mgeo::TapeFrame TapeLine(int N) const;
	double PaperLag() const { return PaperTarget - Paper; }
	double MemoryWheelPosition(int Column, int K) const;   // in digits, 0 .. 10
	double TypeWheelPosition(int K) const;
	double MainShaftAngle() const { return After; }
	double PinionAngle(int Column) const;

private:
	const mgeo::MechanismParts& M;
	std::vector<mgeo::Pose> Current[kMovingParts];
	mgeo::Pose CarriageAt;
	MechanismState Last;
	bool bHaveLast = false;
	double Crank = 0.0;
	double LastCrank = 0.0;
	double After = 0.0;                       // the main shaft and drums
	std::vector<double> Pinion;               // per column
	std::vector<double> MemoryPos, MemoryTarget;   // per memory wheel, then the sign wheels
	std::vector<double> TypePos, TypeTarget;
	double Paper = 0.0, PaperTarget = 0.0;   // in lines
	double Markers[3] = {}, MarkerTargets[3] = {};

	void Recompute();
};
} // namespace mkb
