#include "MechanismPlayer.h"

#include <algorithm>
#include <cctype>
#include <cmath>

namespace mkb
{
namespace player_detail
{
constexpr double kTwoPi = 6.28318530717958647692;

double Smooth(double X)
{
	X = std::min(1.0, std::max(0.0, X));
	return X * X * (3.0 - 2.0 * X);
}

// A register change of -9..9 read as the shortest step on a ten-tooth wheel.
int WheelDelta(int After, int Before)
{
	const int D = ((After - Before) % 10 + 10) % 10;
	return D > 5 ? D - 10 : D;
}
} // namespace player_detail

using namespace player_detail;

double MechanismPlayer::EventSeconds(const mcalc::GearEvent& E, const mcalc::GearEvent& Before) const
{
	switch (E.Type)
	{
	case mcalc::GearEventType::Turn: return 1.0 / TurnsPerSecond;
	case mcalc::GearEventType::Shift: return std::max(0.03, std::abs(E.Shift - Before.Shift) / ShiftsPerSecond);
	case mcalc::GearEventType::Set: return SetSeconds;
	case mcalc::GearEventType::Clear: return ClearSeconds;
	default: return 0.0;
	}
}

void MechanismPlayer::Play(const std::vector<mcalc::GearEvent>& InEvents)
{
	Events = InEvents;
	Index = 0;
	Phase = 0.0;
	// Start from what the wheels show now.
	for (int i = 0; i < MechanismState::SettingDigits; ++i)
	{
		Last.Setting[i] = static_cast<uint8_t>(((static_cast<int>(std::lround(S.Setting[i])) % 10) + 10) % 10);
	}
	for (int i = 0; i < MechanismState::AccumulatorDigits; ++i)
	{
		Last.Accumulator[i] = static_cast<uint8_t>(((static_cast<int>(std::lround(S.Accumulator[i])) % 10) + 10) % 10);
	}
	for (int i = 0; i < MechanismState::CounterDigits; ++i)
	{
		Last.Counter[i] = static_cast<uint8_t>(((static_cast<int>(std::lround(S.Counter[i])) % 10) + 10) % 10);
	}
	Last.Shift = static_cast<int8_t>(std::lround(S.Carriage));
	TotalSeconds = 0.0;
	mcalc::GearEvent Prev = Last;
	for (const mcalc::GearEvent& E : Events)
	{
		TotalSeconds += EventSeconds(E, Prev);
		Prev = E;
	}
	TimeScale = TotalSeconds > MaxSeconds ? TotalSeconds / MaxSeconds : 1.0;
	S.Running = !Events.empty();
}

void MechanismPlayer::ShowSetting(const std::vector<uint8_t>& Digits)
{
	for (int i = 0; i < MechanismState::SettingDigits; ++i)
	{
		IdleSetting[i] = i < static_cast<int>(Digits.size()) ? Digits[static_cast<size_t>(i)] : 0;
	}
}

double MechanismPlayer::Progress() const
{
	return Events.empty() ? 1.0 : static_cast<double>(Index) / static_cast<double>(Events.size());
}

void MechanismPlayer::Settle(const mcalc::GearEvent& E)
{
	if (E.Type == mcalc::GearEventType::Turn)
	{
		CrankBase += kTwoPi;
		S.CrankAngle = CrankBase;
	}
	for (int i = 0; i < MechanismState::SettingDigits; ++i)
	{
		S.Setting[i] = E.Setting[i];
	}
	for (int i = 0; i < MechanismState::AccumulatorDigits; ++i)
	{
		S.Accumulator[i] = E.Accumulator[i];
	}
	for (int i = 0; i < MechanismState::CounterDigits; ++i)
	{
		S.Counter[i] = E.Counter[i];
	}
	S.Carriage = E.Shift;
	S.Op = E.Type == mcalc::GearEventType::End ? ' ' : E.Op;
	Last = E;
}

void MechanismPlayer::Apply(const mcalc::GearEvent& E, const mcalc::GearEvent& Before, double U)
{
	S.Op = E.Op;
	switch (E.Type)
	{
	case mcalc::GearEventType::Turn:
	{
		S.CrankAngle = CrankBase + kTwoPi * U;
		S.Reverse = E.Sign < 0 ? 1.0 : 0.0;
		const double Main = Smooth((U - 0.10) / 0.45);
		for (int p = 0; p < MechanismState::AccumulatorDigits; ++p)
		{
			const int i = p - E.Shift;
			const int Drive = (i >= 0 && i < MechanismState::SettingDigits) ? E.Sign * E.Setting[i] : 0;
			const int Carry = WheelDelta(E.Accumulator[p], Before.Accumulator[p] + Drive);
			const double W0 = 0.60 + 0.33 * p / MechanismState::AccumulatorDigits;
			const double CarryPhase = Smooth((U - W0) / 0.03);
			S.Accumulator[p] = Before.Accumulator[p] + Drive * Main + Carry * CarryPhase;
			if (Carry != 0 && U >= W0 && U <= W0 + 0.05)
			{
				S.CarryLever[p] = std::sin(3.14159265358979 * (U - W0) / 0.05);
			}
		}
		for (int q = 0; q < MechanismState::CounterDigits; ++q)
		{
			S.Counter[q] = Before.Counter[q] + WheelDelta(E.Counter[q], Before.Counter[q]) * Main;
		}
		break;
	}
	case mcalc::GearEventType::Shift:
		S.Carriage = Before.Shift + (E.Shift - Before.Shift) * Smooth(U);
		break;
	case mcalc::GearEventType::Set:
		for (int i = 0; i < MechanismState::SettingDigits; ++i)
		{
			S.Setting[i] = Before.Setting[i] + (E.Setting[i] - Before.Setting[i]) * Smooth(U);
		}
		break;
	case mcalc::GearEventType::Clear:
		// The clearing handle spins every wheel forward to zero.
		for (int p = 0; p < MechanismState::AccumulatorDigits; ++p)
		{
			const int D = ((E.Accumulator[p] - Before.Accumulator[p]) % 10 + 10) % 10;
			S.Accumulator[p] = Before.Accumulator[p] + D * Smooth(U);
		}
		for (int q = 0; q < MechanismState::CounterDigits; ++q)
		{
			const int D = ((E.Counter[q] - Before.Counter[q]) % 10 + 10) % 10;
			S.Counter[q] = Before.Counter[q] + D * Smooth(U);
		}
		break;
	default:
		break;
	}
}

int MechanismPlayer::Step(double Dt)
{
	int Bells = 0;
	const double Decay = std::exp(-Dt * 12.0);
	S.BellHammer *= Decay;
	for (double& L : S.CarryLever)
	{
		L *= Decay;
	}
	if (Index >= Events.size())
	{
		S.Running = false;
		S.Op = ' ';
		S.Reverse *= Decay;
		// Idle: the setting pinions follow the number being typed.
		for (int i = 0; i < MechanismState::SettingDigits; ++i)
		{
			S.Setting[i] += (IdleSetting[i] - S.Setting[i]) * (1.0 - std::exp(-Dt * 14.0));
		}
		return 0;
	}
	double Left = Dt * TimeScale;
	while (Left > 0.0 && Index < Events.size())
	{
		const mcalc::GearEvent& E = Events[Index];
		const double Duration = EventSeconds(E, Last);
		if (E.Type == mcalc::GearEventType::Bell)
		{
			++Bells;
			S.BellHammer = 1.0;
		}
		if (Duration <= 0.0 || Left >= Duration - Phase)
		{
			Left -= std::max(0.0, Duration - Phase);
			Settle(E);
			Phase = 0.0;
			++Index;
			continue;
		}
		Phase += Left;
		Left = 0.0;
		Apply(E, Last, Phase / Duration);
	}
	S.Running = Index < Events.size();
	return Bells;
}

void MechanismPlayer::SkipToEnd()
{
	if (!Events.empty() && Index < Events.size())
	{
		Settle(Events.back());
		Index = Events.size();
	}
	S.Running = false;
}

std::vector<uint8_t> SettingFromInput(const std::string& Line, int Cursor)
{
	std::vector<uint8_t> Digits;
	if (Cursor <= 0 || Cursor > static_cast<int>(Line.size()))
	{
		return Digits;
	}
	auto InNumber = [&Line](int i) { return i >= 0 && i < static_cast<int>(Line.size()) && (std::isdigit(static_cast<unsigned char>(Line[static_cast<size_t>(i)])) != 0 || Line[static_cast<size_t>(i)] == '.'); };
	if (!InNumber(Cursor - 1))
	{
		return Digits;
	}
	int B = Cursor - 1;
	while (InNumber(B - 1))
	{
		--B;
	}
	// A number is not part of a name (x2) or an exponent (1e5's 5 belongs to the exponent).
	if (B > 0 && (std::isalpha(static_cast<unsigned char>(Line[static_cast<size_t>(B - 1)])) != 0 || Line[static_cast<size_t>(B - 1)] == '_'))
	{
		return Digits;
	}
	int E = Cursor;
	while (InNumber(E))
	{
		++E;
	}
	for (int i = E - 1; i >= B && Digits.size() < static_cast<size_t>(MechanismState::SettingDigits); --i)
	{
		if (Line[static_cast<size_t>(i)] != '.')
		{
			Digits.push_back(static_cast<uint8_t>(Line[static_cast<size_t>(i)] - '0'));
		}
	}
	return Digits;
}
} // namespace mkb
