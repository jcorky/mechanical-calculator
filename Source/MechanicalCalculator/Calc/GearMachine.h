// The stepped-drum calculating machine, wheel by wheel.
//
//   setting register      20 digits; under each a Leibniz drum with nine stepped teeth, and a
//                         sliding pinion that the set digit d positions to take d teeth per turn
//   accumulator           40 numeral wheels on the carriage, each with a tens-carry lever
//   counter               20 numeral wheels on the carriage; the one at the carriage position
//                         advances by one per turn (forward) or goes back one (reversed)
//   carriage              shifts 0..20 places: drum i drives accumulator wheel i + shift
//   crank                 one turn adds (or, reversed, subtracts) the setting into the
//                         accumulator at the carriage position; carries ripple up the wheels;
//                         a borrow out of the top wheel rings the bell (overdraft)
//
// The procedures below are the operator's methods for + - * / and square root. They return the
// exact register contents; Real (GearReal.h) rounds them to 18 digits. Real's fast path computes
// the same contents with integer arithmetic, and Tests/TestGear.cpp requires identical results.
// A GearRecorder captures every set, shift and turn with full register snapshots so the Unreal
// layer can animate the drums, pinions, wheels, carries and carriage.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mcalc
{
enum class GearEventType : uint8_t
{
	Begin,     // a new procedure starts (Op says which)
	Clear,     // accumulator and counter cleared
	Set,       // the setting register was set
	Shift,     // the carriage moved
	Turn,      // one crank revolution, forward (Sign +1) or reversed (-1)
	Bell,      // overdraft: the turn borrowed out of the top wheel
	End,       // the procedure finished; the result is read off the registers
};

struct GearEvent
{
	static constexpr int SettingDigits = 20;
	static constexpr int AccumulatorDigits = 40;
	static constexpr int CounterDigits = 21;   // one wheel per carriage position 0..20

	GearEventType Type = GearEventType::Begin;
	char Op = '+';            // '+', '-', '*', '/', 'r' (square root)
	int8_t Sign = 0;          // Turn: +1 forward, -1 reversed
	int8_t Shift = 0;         // carriage position after the event
	uint64_t CarryMask = 0;   // Turn: accumulator wheels that passed 9 -> 0 (or 0 -> 9) and carried
	uint8_t Setting[SettingDigits] = {};          // digit i is under drum i (units first)
	uint8_t Accumulator[AccumulatorDigits] = {};  // wheel 0 is the units wheel
	uint8_t Counter[CounterDigits] = {};
};

class GearRecorder
{
public:
	// Records up to MaxEvents events. With bKeepEvents false it records nothing but still routes
	// every operation through the literal machine (the tests use this to compare the two paths).
	explicit GearRecorder(size_t InMaxEvents = 20000, bool bInKeepEvents = true) : MaxEvents(InMaxEvents), bKeepEvents(bInKeepEvents) {}
	bool Wants() const { return !bKeepEvents || Events.size() < MaxEvents; }
	void Record(const GearEvent& E)
	{
		if (bKeepEvents && Events.size() < MaxEvents)
		{
			Events.push_back(E);
		}
	}
	const std::vector<GearEvent>& GetEvents() const { return Events; }
	void Clear() { Events.clear(); }

private:
	size_t MaxEvents;
	bool bKeepEvents;
	std::vector<GearEvent> Events;
};

class SteppedDrumMachine
{
public:
	static constexpr int SettingDigits = GearEvent::SettingDigits;
	static constexpr int AccumulatorDigits = GearEvent::AccumulatorDigits;
	static constexpr int CounterDigits = GearEvent::CounterDigits;
	static constexpr int MaxShift = 20;

	explicit SteppedDrumMachine(GearRecorder* InRecorder = nullptr, char InOp = '+') : Recorder(InRecorder), Op(InOp) {}

	void Begin();
	void ClearRegisters();
	void ClearCounter();
	// In division and root extraction the counter is reversed: it counts subtracting turns up.
	void SetCounterReversed(bool bReversed) { CounterReversed = bReversed; }
	// Sets the setting register to Digits (units first), placed Offset places up the register.
	void Set(const std::vector<uint8_t>& Digits, int Offset = 0);
	void ShiftTo(int Position);
	// One crank revolution. Returns true when the bell rings (a borrow left the top wheel: overdraft).
	bool Turn(int Sign);
	void End();

	const uint8_t* Accumulator() const { return Acc; }
	const uint8_t* Counter() const { return Cnt; }
	int Shift() const { return Carriage; }
	unsigned long long Turns() const { return TurnCount; }
	unsigned long long Shifts() const { return ShiftCount; }
	unsigned long long Carries() const { return CarryCount; }

private:
	GearRecorder* Recorder;
	char Op;
	uint8_t Setting[SettingDigits] = {};
	uint8_t Acc[AccumulatorDigits] = {};
	uint8_t Cnt[CounterDigits] = {};
	int Carriage = 0;
	bool CounterReversed = false;
	unsigned long long TurnCount = 0;
	unsigned long long ShiftCount = 0;
	unsigned long long CarryCount = 0;

	void Emit(GearEventType Type, int Sign, uint64_t CarryMask);
};

// The operator's procedures on the literal machine. Mantissas are 18-digit integers (or 0).
// Each returns the exact result as decimal digits, units first, read off the registers, and adds
// its turns, shifts and carries to GearContext::Stats().
namespace GearProcedures
{
// |A| * 10^Gap + |B| (or minus, larger first: Subtract requires A * 10^Gap >= B). Gap <= 20.
std::vector<uint8_t> AddAligned(uint64_t A, int Gap, uint64_t B, bool Subtract, GearRecorder* Recorder);
std::vector<uint8_t> Multiply(uint64_t A, uint64_t B, GearRecorder* Recorder);
// Nineteen significant quotient digits of A / B (the counter's reading) and whether A < B.
std::vector<uint8_t> Divide(uint64_t A, uint64_t B, GearRecorder* Recorder, bool* SmallerThanDivisor);
// floor(sqrt(N)) where N is given as digits, units first (37 or 38 digits): Toepler's method.
std::vector<uint8_t> SquareRoot(const std::vector<uint8_t>& N, GearRecorder* Recorder);
} // namespace GearProcedures
} // namespace mcalc
