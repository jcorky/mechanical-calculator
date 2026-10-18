#include "MechanismAnimator.h"

#include <algorithm>
#include <cmath>

namespace mkb
{
namespace animator_detail
{
// Moves a wheel position (digits, 0 .. 10) forward toward Target by at most MaxStep, the way a
// figure wheel only ever turns one way.
double SpinToward(double Position, double Target, double MaxStep)
{
	double Ahead = std::fmod(Target - Position, 10.0);
	if (Ahead < 0.0)
	{
		Ahead += 10.0;
	}
	if (Ahead < 1e-9 || Ahead > 10.0 - 1e-9)
	{
		return Target;
	}
	const double Next = Position + std::min(Ahead, MaxStep);
	return std::fmod(Next, 10.0);
}

double ApproachLinear(double Value, double Target, double MaxStep)
{
	if (std::abs(Target - Value) <= MaxStep)
	{
		return Target;
	}
	return Value + (Target > Value ? MaxStep : -MaxStep);
}

double WheelTarget(const WheelReading& R, int K)
{
	if (R.Empty)
	{
		return 0.0;
	}
	return R.Special ? 9.0 : R.Digits[K];
}
} // namespace animator_detail

using namespace animator_detail;

const char* MovingPartName(MovingPart P)
{
	switch (P)
	{
	case MovingPart::Drums: return "Drum";
	case MovingPart::SquareShafts: return "SquareShaft";
	case MovingPart::Pinions: return "Pinion";
	case MovingPart::TransferGears: return "TransferGear";
	case MovingPart::SettingWheels: return "SettingWheel";
	case MovingPart::SettingMarker: return "SettingMarker";
	case MovingPart::Carriage: return "Carriage";
	case MovingPart::AccumulatorWheels: return "AccumulatorWheel";
	case MovingPart::CarryLevers: return "CarryLever";
	case MovingPart::CounterWheels: return "CounterWheel";
	case MovingPart::CarriageMarkers: return "CarriageMarker";
	case MovingPart::ClearingLevers: return "ClearingLever";
	case MovingPart::CountingFinger: return "CountingFinger";
	case MovingPart::MainShaft: return "MainShaft";
	case MovingPart::ReverserIdler: return "ReverserIdler";
	case MovingPart::InputShaft: return "InputShaft";
	case MovingPart::DriveGear: return "DriveGear";
	case MovingPart::Flywheel: return "Flywheel";
	case MovingPart::Crank: return "Crank";
	case MovingPart::CrankHandle: return "CrankHandle";
	case MovingPart::ReversingLever: return "ReversingLever";
	case MovingPart::ShiftKnob: return "ShiftKnob";
	case MovingPart::ShiftPinion: return "ShiftPinion";
	case MovingPart::BellHammer: return "BellHammer";
	case MovingPart::TypeWheels: return "TypeWheel";
	case MovingPart::PaperRoll: return "PaperRoll";
	case MovingPart::Platen: return "Platen";
	case MovingPart::RibbonSpools: return "RibbonSpool";
	case MovingPart::MemoryWheels: return "MemoryWheel";
	case MovingPart::MemorySignWheels: return "MemorySignWheel";
	default: return "Part";
	}
}

const mgeo::PartSet& PartSetOf(const mgeo::MechanismParts& M, MovingPart P)
{
	switch (P)
	{
	case MovingPart::Drums: return M.Drums;
	case MovingPart::SquareShafts: return M.SquareShafts;
	case MovingPart::Pinions: return M.Pinions;
	case MovingPart::TransferGears: return M.TransferGears;
	case MovingPart::SettingWheels: return M.SettingWheels;
	case MovingPart::SettingMarker: return M.SettingMarker;
	case MovingPart::Carriage: return M.Carriage;
	case MovingPart::AccumulatorWheels: return M.AccumulatorWheels;
	case MovingPart::CarryLevers: return M.CarryLevers;
	case MovingPart::CounterWheels: return M.CounterWheels;
	case MovingPart::CarriageMarkers: return M.CarriageMarkers;
	case MovingPart::ClearingLevers: return M.ClearingLevers;
	case MovingPart::CountingFinger: return M.CountingFinger;
	case MovingPart::MainShaft: return M.MainShaft;
	case MovingPart::ReverserIdler: return M.ReverserIdler;
	case MovingPart::InputShaft: return M.InputShaft;
	case MovingPart::DriveGear: return M.DriveGear;
	case MovingPart::Flywheel: return M.Flywheel;
	case MovingPart::Crank: return M.Crank;
	case MovingPart::CrankHandle: return M.CrankHandle;
	case MovingPart::ReversingLever: return M.ReversingLever;
	case MovingPart::ShiftKnob: return M.ShiftKnob;
	case MovingPart::ShiftPinion: return M.ShiftPinion;
	case MovingPart::BellHammer: return M.BellHammer;
	case MovingPart::TypeWheels: return M.TypeWheels;
	case MovingPart::PaperRoll: return M.PaperRoll;
	case MovingPart::Platen: return M.Platen;
	case MovingPart::RibbonSpools: return M.RibbonSpools;
	case MovingPart::MemoryWheels: return M.MemoryWheels;
	default: return M.MemorySignWheels;
	}
}

bool RidesCarriage(const mgeo::MechanismParts& M, MovingPart P)
{
	return P == MovingPart::Carriage || PartSetOf(M, P).OnCarriage;
}

MechanismAnimator::MechanismAnimator(const mgeo::MechanismParts& InParts) : M(InParts)
{
	Pinion.assign(M.Pinions.Poses.size(), 0.0);
	const size_t Memory = M.MemoryWheels.Poses.size() + M.MemorySignWheels.Poses.size();
	MemoryPos.assign(Memory, 0.0);
	MemoryTarget.assign(Memory, 0.0);
	TypePos.assign(M.TypeWheels.Poses.size(), 0.0);
	TypeTarget.assign(M.TypeWheels.Poses.size(), 0.0);
	CarriageAt = M.Carriage.Poses.empty() ? mgeo::Pose() : M.Carriage.Poses[0];
	Recompute();
}

void MechanismAnimator::Follow(const MechanismState& S, double ExtraCrankRad)
{
	Crank = S.CrankAngle + ExtraCrankRad;
	double Delta = bHaveLast ? Crank - LastCrank : 0.0;
	if (!(std::abs(Delta) < 1e3))
	{
		Delta = 0.0;   // a reset or a jump: no motion to integrate
	}
	LastCrank = Crank;
	// The reverser flips everything after it; the lever only moves while the crank rests.
	const double Sign = S.Reverse > 0.5 ? -1.0 : 1.0;
	After += Sign * Delta;
	for (size_t J = 0; J < Pinion.size(); ++J)
	{
		const double Digit = J < static_cast<size_t>(MechanismState::SettingDigits) ? S.Setting[J] : 0.0;
		Pinion[J] += Sign * Delta * Digit / 12.0;
	}
	Last = S;
	bHaveLast = true;
	Recompute();
}

void MechanismAnimator::ShowMemory(const std::vector<MemoryColumn>& Columns)
{
	const int PerColumn = M.MemoryMantissa + M.MemoryExponent;
	const size_t Signs = M.MemoryWheels.Poses.size();
	for (int C = 0; C < M.MemoryColumns; ++C)
	{
		WheelReading Blank;
		const WheelReading& R = C < static_cast<int>(Columns.size()) ? Columns[static_cast<size_t>(C)].Reading : Blank;
		for (int K = 0; K < PerColumn; ++K)
		{
			const size_t I = static_cast<size_t>(C * PerColumn + K);
			if (I >= Signs)
			{
				continue;
			}
			double T = 0.0;
			if (K < M.MemoryMantissa)
			{
				T = K < 18 ? WheelTarget(R, K) : 0.0;
			}
			else if (!R.Empty && !R.Special)
			{
				T = R.ExponentDigits[std::min(2, K - M.MemoryMantissa)];
			}
			MemoryTarget[I] = T;
		}
		const size_t SignIndex = Signs + static_cast<size_t>(C);
		if (SignIndex < MemoryTarget.size())
		{
			// '+' at even positions, '-' at odd; inf and nan leave the sign wheel half-turned.
			MemoryTarget[SignIndex] = R.Empty ? 0.0 : R.Special ? 0.5 : (R.Negative ? 1.0 : 0.0);
		}
	}
}

void MechanismAnimator::ShowTypeWheels(const WheelReading& R)
{
	for (size_t K = 0; K < TypeTarget.size(); ++K)
	{
		TypeTarget[K] = K < 18 ? WheelTarget(R, static_cast<int>(K)) : 0.0;
	}
}

void MechanismAnimator::Feed(int Lines)
{
	if (Lines > 0)
	{
		PaperTarget += Lines;
	}
}

void MechanismAnimator::SetDecimalMarkers(int SettingPlaces, int AccumulatorPlaces, int CounterPlaces)
{
	MarkerTargets[0] = std::max(0, std::min(SettingPlaces, M.SettingDigits - 1));
	MarkerTargets[1] = std::max(0, std::min(AccumulatorPlaces, M.AccumulatorDigits - 1));
	MarkerTargets[2] = std::max(0, std::min(CounterPlaces, M.CounterDigits - 1));
}

void MechanismAnimator::Step(double Dt)
{
	if (Dt <= 0.0)
	{
		return;
	}
	const double Spin = FigureWheelDigitsPerSecond * Dt;
	for (size_t I = 0; I < MemoryPos.size(); ++I)
	{
		MemoryPos[I] = SpinToward(MemoryPos[I], MemoryTarget[I], Spin);
	}
	for (size_t K = 0; K < TypePos.size(); ++K)
	{
		TypePos[K] = SpinToward(TypePos[K], TypeTarget[K], Spin);
	}
	// The feed hurries when several lines are waiting.
	const double Lag = PaperTarget - Paper;
	Paper = ApproachLinear(Paper, PaperTarget, std::max(PaperLinesPerSecond, 4.0 * Lag) * Dt);
	for (int I = 0; I < 3; ++I)
	{
		Markers[I] = ApproachLinear(Markers[I], MarkerTargets[I], 12.0 * Dt);
	}
	Recompute();
}

bool MechanismAnimator::Settled() const
{
	for (size_t I = 0; I < MemoryPos.size(); ++I)
	{
		if (MemoryPos[I] != MemoryTarget[I])
		{
			return false;
		}
	}
	for (size_t K = 0; K < TypePos.size(); ++K)
	{
		if (TypePos[K] != TypeTarget[K])
		{
			return false;
		}
	}
	for (int I = 0; I < 3; ++I)
	{
		if (Markers[I] != MarkerTargets[I])
		{
			return false;
		}
	}
	return Paper == PaperTarget;
}

mgeo::TapeFrame MechanismAnimator::TapeLine(int N) const
{
	const double Distance = (N + 0.5 - (PaperTarget - Paper)) * M.TapeLinePitchMm;
	return mgeo::TapeLineFrame(M, std::max(0.0, Distance));
}

double MechanismAnimator::MemoryWheelPosition(int Column, int K) const
{
	const size_t I = static_cast<size_t>(Column * (M.MemoryMantissa + M.MemoryExponent) + K);
	return I < MemoryPos.size() ? MemoryPos[I] : 0.0;
}

double MechanismAnimator::TypeWheelPosition(int K) const
{
	return K >= 0 && K < static_cast<int>(TypePos.size()) ? TypePos[static_cast<size_t>(K)] : 0.0;
}

double MechanismAnimator::PinionAngle(int Column) const
{
	return Column >= 0 && Column < static_cast<int>(Pinion.size()) ? Pinion[static_cast<size_t>(Column)] : 0.0;
}

void MechanismAnimator::Recompute()
{
	const MechanismState& S = Last;
	const double Digit = M.DigitAngle;
	auto Base = [this](MovingPart P) -> const std::vector<mgeo::Pose>& { return PartSetOf(M, P).Poses; };
	auto Out = [this](MovingPart P) -> std::vector<mgeo::Pose>& {
		std::vector<mgeo::Pose>& V = Current[static_cast<int>(P)];
		V = PartSetOf(M, P).Poses;
		return V;
	};
	auto TurnAll = [&Out](MovingPart P, double Angle) {
		for (mgeo::Pose& Pose : Out(P))
		{
			Pose = Pose.RotatedZ(Angle);
		}
	};
	auto Clamped = [](const double* Values, size_t Count, size_t I) { return I < Count ? Values[I] : 0.0; };

	// Before the reverser: the crank, its shaft and the big drive gear.
	TurnAll(MovingPart::Crank, Crank);
	TurnAll(MovingPart::CrankHandle, Crank);
	TurnAll(MovingPart::InputShaft, Crank);
	TurnAll(MovingPart::DriveGear, Crank);
	TurnAll(MovingPart::ReverserIdler, Crank);
	TurnAll(MovingPart::Flywheel, -M.FlywheelRatio * Crank);
	// After it: the main shaft, the drums and the counting finger.
	TurnAll(MovingPart::MainShaft, After);
	TurnAll(MovingPart::Drums, After);
	TurnAll(MovingPart::CountingFinger, After);

	// The setting columns: each pinion slides to its digit and turns with what it engages.
	{
		std::vector<mgeo::Pose>& Pin = Out(MovingPart::Pinions);
		std::vector<mgeo::Pose>& Square = Out(MovingPart::SquareShafts);
		std::vector<mgeo::Pose>& Transfer = Out(MovingPart::TransferGears);
		std::vector<mgeo::Pose>& Wheels = Out(MovingPart::SettingWheels);
		for (size_t J = 0; J < Pin.size(); ++J)
		{
			const double Set = Clamped(S.Setting, MechanismState::SettingDigits, J);
			Pin[J] = Pin[J].Moved(mgeo::Vec3(0.0, 0.0, Set * M.PinionStepMm)).RotatedZ(Pinion[J]);
		}
		for (size_t J = 0; J < Square.size(); ++J)
		{
			Square[J] = Square[J].RotatedZ(PinionAngle(static_cast<int>(J)));
		}
		for (size_t J = 0; J < Transfer.size(); ++J)
		{
			Transfer[J] = Transfer[J].RotatedZ(PinionAngle(static_cast<int>(J)));
		}
		for (size_t J = 0; J < Wheels.size(); ++J)
		{
			Wheels[J] = Wheels[J].RotatedZ(Clamped(S.Setting, MechanismState::SettingDigits, J) * Digit);
		}
	}
	for (mgeo::Pose& P : Out(MovingPart::SettingMarker))
	{
		P = P.Moved(mgeo::Vec3(Markers[0] * M.DigitPitchMm, 0.0, 0.0));
	}

	// The carriage and what rides on it.
	const mgeo::Pose CarriageHome = M.Carriage.Poses.empty() ? mgeo::Pose() : M.Carriage.Poses[0];
	CarriageAt = CarriageHome.Moved(mgeo::Vec3(S.Carriage * M.DigitPitchMm, 0.0, 0.0));
	Current[static_cast<int>(MovingPart::Carriage)].assign(Base(MovingPart::Carriage).empty() ? 0 : 1, mgeo::Pose());
	{
		std::vector<mgeo::Pose>& Acc = Out(MovingPart::AccumulatorWheels);
		for (size_t I = 0; I < Acc.size(); ++I)
		{
			Acc[I] = Acc[I].RotatedZ(Clamped(S.Accumulator, MechanismState::AccumulatorDigits, I) * Digit);
		}
		std::vector<mgeo::Pose>& Levers = Out(MovingPart::CarryLevers);
		for (size_t I = 0; I < Levers.size(); ++I)
		{
			Levers[I] = Levers[I].RotatedZ(Clamped(S.CarryLever, MechanismState::AccumulatorDigits, I) * M.CarryLeverSwing);
		}
		std::vector<mgeo::Pose>& Count = Out(MovingPart::CounterWheels);
		for (size_t I = 0; I < Count.size(); ++I)
		{
			Count[I] = Count[I].RotatedZ(Clamped(S.Counter, MechanismState::CounterDigits, I) * Digit);
		}
		std::vector<mgeo::Pose>& Marks = Out(MovingPart::CarriageMarkers);
		for (size_t I = 0; I < Marks.size(); ++I)
		{
			Marks[I] = Marks[I].Moved(mgeo::Vec3(Markers[I == 0 ? 1 : 2] * M.DigitPitchMm, 0.0, 0.0));
		}
		TurnAll(MovingPart::ClearingLevers, S.Clearing * M.ClearingLeverSwing);
	}

	// Levers, knob and bell.
	TurnAll(MovingPart::ReversingLever, M.ReverseAddAngle + (M.ReverseSubtractAngle - M.ReverseAddAngle) * S.Reverse);
	TurnAll(MovingPart::ShiftKnob, M.ShiftKnobAnglePerStep * S.Carriage);
	TurnAll(MovingPart::ShiftPinion, M.ShiftKnobAnglePerStep * S.Carriage);
	TurnAll(MovingPart::BellHammer, S.BellHammer * M.BellStrikeAngle);

	// The printer.
	{
		std::vector<mgeo::Pose>& Type = Out(MovingPart::TypeWheels);
		for (size_t K = 0; K < Type.size(); ++K)
		{
			Type[K] = Type[K].RotatedZ(TypePos[K] * Digit);
		}
		const double Fed = Paper * M.TapeLinePitchMm;
		TurnAll(MovingPart::PaperRoll, M.PaperRollRadius > 0.0 ? Fed / M.PaperRollRadius : 0.0);
		TurnAll(MovingPart::Platen, M.PlatenRadius > 0.0 ? Fed / M.PlatenRadius : 0.0);
		TurnAll(MovingPart::RibbonSpools, Paper * 0.35);
	}

	// The memory columns.
	{
		std::vector<mgeo::Pose>& Mem = Out(MovingPart::MemoryWheels);
		for (size_t I = 0; I < Mem.size(); ++I)
		{
			Mem[I] = Mem[I].RotatedZ(MemoryPos[I] * Digit);
		}
		std::vector<mgeo::Pose>& Signs = Out(MovingPart::MemorySignWheels);
		for (size_t C = 0; C < Signs.size(); ++C)
		{
			const size_t I = M.MemoryWheels.Poses.size() + C;
			Signs[C] = Signs[C].RotatedZ((I < MemoryPos.size() ? MemoryPos[I] : 0.0) * Digit);
		}
	}
}
} // namespace mkb
