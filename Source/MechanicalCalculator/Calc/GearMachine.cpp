#include "GearMachine.h"

#include "GearReal.h"

#include <algorithm>

namespace mcalc
{
namespace machine_detail
{
std::vector<uint8_t> DigitsOf(uint64_t X)
{
	std::vector<uint8_t> D;
	while (X > 0)
	{
		D.push_back(static_cast<uint8_t>(X % 10));
		X /= 10;
	}
	return D;
}

// Small big-integer helpers on digit vectors (units first), for the operator's notes in Toepler's method.
void Trim(std::vector<uint8_t>& D)
{
	while (!D.empty() && D.back() == 0)
	{
		D.pop_back();
	}
}

std::vector<uint8_t> MulSmall(const std::vector<uint8_t>& D, unsigned M)
{
	std::vector<uint8_t> Out;
	unsigned Carry = 0;
	for (uint8_t X : D)
	{
		const unsigned V = X * M + Carry;
		Out.push_back(static_cast<uint8_t>(V % 10));
		Carry = V / 10;
	}
	while (Carry > 0)
	{
		Out.push_back(static_cast<uint8_t>(Carry % 10));
		Carry /= 10;
	}
	Trim(Out);
	return Out;
}

std::vector<uint8_t> AddSmall(const std::vector<uint8_t>& D, unsigned M)
{
	std::vector<uint8_t> Out = D;
	unsigned Carry = M;
	for (size_t i = 0; Carry > 0; ++i)
	{
		if (i == Out.size())
		{
			Out.push_back(0);
		}
		const unsigned V = Out[i] + Carry;
		Out[i] = static_cast<uint8_t>(V % 10);
		Carry = V / 10;
	}
	return Out;
}

void Account(const SteppedDrumMachine& M)
{
	GearStats& S = GearContext::Stats();
	S.Turns += M.Turns();
	S.Shifts += M.Shifts();
	S.Carries += M.Carries();
}
} // namespace machine_detail

using namespace machine_detail;

void SteppedDrumMachine::Emit(GearEventType Type, int Sign, uint64_t CarryMask)
{
	if (!Recorder || !Recorder->Wants())
	{
		return;
	}
	GearEvent E;
	E.Type = Type;
	E.Op = Op;
	E.Sign = static_cast<int8_t>(Sign);
	E.Shift = static_cast<int8_t>(Carriage);
	E.CarryMask = CarryMask;
	std::copy(Setting, Setting + SettingDigits, E.Setting);
	std::copy(Acc, Acc + AccumulatorDigits, E.Accumulator);
	std::copy(Cnt, Cnt + CounterDigits, E.Counter);
	Recorder->Record(E);
}

void SteppedDrumMachine::Begin()
{
	Emit(GearEventType::Begin, 0, 0);
}

void SteppedDrumMachine::End()
{
	Emit(GearEventType::End, 0, 0);
}

void SteppedDrumMachine::ClearRegisters()
{
	std::fill(Acc, Acc + AccumulatorDigits, 0);
	std::fill(Cnt, Cnt + CounterDigits, 0);
	Emit(GearEventType::Clear, 0, 0);
}

void SteppedDrumMachine::ClearCounter()
{
	std::fill(Cnt, Cnt + CounterDigits, 0);
	Emit(GearEventType::Clear, 0, 0);
}

void SteppedDrumMachine::Set(const std::vector<uint8_t>& Digits, int Offset)
{
	std::fill(Setting, Setting + SettingDigits, 0);
	for (size_t i = 0; i < Digits.size(); ++i)
	{
		const int Pos = Offset + static_cast<int>(i);
		if (Pos >= 0 && Pos < SettingDigits)
		{
			Setting[Pos] = Digits[i];
		}
	}
	Emit(GearEventType::Set, 0, 0);
}

void SteppedDrumMachine::ShiftTo(int Position)
{
	Position = std::max(0, std::min(MaxShift, Position));
	if (Position == Carriage)
	{
		return;
	}
	ShiftCount += static_cast<unsigned long long>(Position > Carriage ? Position - Carriage : Carriage - Position);
	Carriage = Position;
	Emit(GearEventType::Shift, 0, 0);
}

bool SteppedDrumMachine::Turn(int Sign)
{
	++TurnCount;
	int Pending[AccumulatorDigits + 1] = {};
	uint64_t Mask = 0;
	bool Lost = false;
	// Every drum turns its pinion by its digit's number of teeth at once...
	for (int i = 0; i < SettingDigits; ++i)
	{
		if (Setting[i] == 0)
		{
			continue;
		}
		const int Pos = i + Carriage;
		if (Pos >= AccumulatorDigits)
		{
			Lost = true;
			continue;
		}
		int V = Acc[Pos] + Sign * Setting[i];
		if (V > 9)
		{
			V -= 10;
			Pending[Pos + 1] += 1;
			Mask |= uint64_t(1) << Pos;
		}
		else if (V < 0)
		{
			V += 10;
			Pending[Pos + 1] -= 1;
			Mask |= uint64_t(1) << Pos;
		}
		Acc[Pos] = static_cast<uint8_t>(V);
	}
	// ...then the carry levers trip in turn from the units wheel upward.
	for (int Pos = 0; Pos < AccumulatorDigits; ++Pos)
	{
		if (Pending[Pos] == 0)
		{
			continue;
		}
		++CarryCount;
		int V = Acc[Pos] + Pending[Pos];
		if (V > 9)
		{
			V -= 10;
			Pending[Pos + 1] += 1;
			Mask |= uint64_t(1) << Pos;
		}
		else if (V < 0)
		{
			V += 10;
			Pending[Pos + 1] -= 1;
			Mask |= uint64_t(1) << Pos;
		}
		Acc[Pos] = static_cast<uint8_t>(V);
	}
	// The counter wheel at the carriage position moves one step.
	int Step = CounterReversed ? -Sign : Sign;
	for (int Pos = Carriage; Pos < CounterDigits && Step != 0; ++Pos)
	{
		int V = Cnt[Pos] + Step;
		Step = 0;
		if (V > 9)
		{
			V -= 10;
			Step = 1;
		}
		else if (V < 0)
		{
			V += 10;
			Step = -1;
		}
		Cnt[Pos] = static_cast<uint8_t>(V);
	}
	// Only an overdraft rings: a borrow out of the top wheel. (The carry out of the top on the
	// restoring turn just brings the complement back.)
	const bool Bell = Lost || Pending[AccumulatorDigits] < 0;
	Emit(GearEventType::Turn, Sign, Mask);
	if (Bell)
	{
		Emit(GearEventType::Bell, Sign, 0);
	}
	return Bell;
}

namespace GearProcedures
{
std::vector<uint8_t> AddAligned(uint64_t A, int Gap, uint64_t B, bool Subtract, GearRecorder* Recorder)
{
	SteppedDrumMachine M(Recorder, Subtract ? '-' : '+');
	M.Begin();
	M.ClearRegisters();
	M.Set(DigitsOf(A));
	M.ShiftTo(Gap);
	M.Turn(1);
	M.Set(DigitsOf(B));
	M.ShiftTo(0);
	M.Turn(Subtract ? -1 : 1);
	M.End();
	Account(M);
	return std::vector<uint8_t>(M.Accumulator(), M.Accumulator() + SteppedDrumMachine::AccumulatorDigits);
}

std::vector<uint8_t> Multiply(uint64_t A, uint64_t B, GearRecorder* Recorder)
{
	SteppedDrumMachine M(Recorder, '*');
	M.Begin();
	M.ClearRegisters();
	M.Set(DigitsOf(A));
	const std::vector<uint8_t> Multiplier = DigitsOf(B);
	for (size_t k = 0; k < Multiplier.size(); ++k)
	{
		M.ShiftTo(static_cast<int>(k));
		for (int t = 0; t < Multiplier[k]; ++t)
		{
			M.Turn(1);
		}
	}
	M.ShiftTo(0);
	M.End();
	Account(M);
	return std::vector<uint8_t>(M.Accumulator(), M.Accumulator() + SteppedDrumMachine::AccumulatorDigits);
}

std::vector<uint8_t> Divide(uint64_t A, uint64_t B, GearRecorder* Recorder, bool* SmallerThanDivisor)
{
	SteppedDrumMachine M(Recorder, '/');
	M.Begin();
	M.ClearRegisters();
	// The dividend goes into the top of the accumulator.
	M.Set(DigitsOf(A));
	M.ShiftTo(SteppedDrumMachine::MaxShift);
	M.Turn(1);
	M.ClearCounter();
	M.SetCounterReversed(true);
	M.Set(DigitsOf(B));
	for (int S = SteppedDrumMachine::MaxShift; S >= 1; --S)
	{
		M.ShiftTo(S);
		// Turn backwards until the bell, then once forward to restore.
		while (!M.Turn(-1))
		{
		}
		M.Turn(1);
	}
	M.ShiftTo(0);
	M.End();
	Account(M);
	const uint8_t* C = M.Counter();
	*SmallerThanDivisor = C[SteppedDrumMachine::MaxShift] == 0;
	const int Top = *SmallerThanDivisor ? SteppedDrumMachine::MaxShift - 1 : SteppedDrumMachine::MaxShift;
	std::vector<uint8_t> Q;
	for (int Pos = Top - 18; Pos <= Top; ++Pos)
	{
		Q.push_back(C[Pos]);
	}
	return Q;
}

std::vector<uint8_t> SquareRoot(const std::vector<uint8_t>& N, GearRecorder* Recorder)
{
	SteppedDrumMachine M(Recorder, 'r');
	M.Begin();
	M.ClearRegisters();
	// Load N in two settings: the low 20 digits, then the rest 20 places up.
	std::vector<uint8_t> Low(N.begin(), N.begin() + static_cast<std::ptrdiff_t>(std::min<size_t>(20, N.size())));
	M.Set(Low);
	M.ShiftTo(0);
	M.Turn(1);
	if (N.size() > 20)
	{
		M.Set(std::vector<uint8_t>(N.begin() + 20, N.end()));
		M.ShiftTo(20);
		M.Turn(1);
	}
	M.ClearCounter();
	M.SetCounterReversed(true);
	const int Pairs = static_cast<int>((N.size() + 1) / 2);
	std::vector<uint8_t> Root;   // y, units first
	for (int j = Pairs - 1; j >= 0; --j)
	{
		const int Target = 2 * j;
		const int Carriage = std::min(SteppedDrumMachine::MaxShift, Target);
		const int Offset = Target - Carriage;
		const std::vector<uint8_t> Twenty = MulSmall(Root, 20);
		unsigned X = 0;
		for (unsigned Odd = 1;; Odd += 2)
		{
			// Subtract 20y + odd at column pair j; the bell means it did not go: restore and stop.
			M.Set(AddSmall(Twenty, Odd), Offset);
			M.ShiftTo(Carriage);
			if (M.Turn(-1))
			{
				M.Turn(1);
				break;
			}
			++X;
		}
		Root.insert(Root.begin(), static_cast<uint8_t>(X));
		Trim(Root);
	}
	M.ShiftTo(0);
	M.End();
	Account(M);
	return Root;
}
} // namespace GearProcedures
} // namespace mcalc
