// What the figure wheels show when the machine is at rest: a value as the operator would set it
// on a column of wheels (18 digits, a sign wheel, a three-wheel power of ten) and where the sliding
// decimal marker goes; and the memory columns, which hold ans and the session's variables.
#pragma once

#include "CalcSession.h"

#include <cstdint>
#include <string>
#include <vector>

namespace mkb
{
struct WheelReading
{
	uint8_t Digits[18] = {};   // most significant first, as read left to right
	bool Negative = false;
	int Exponent = 0;          // power of ten shown on the exponent wheels (0 when fixed)
	bool ExponentNegative = false;
	uint8_t ExponentDigits[3] = {};
	int DecimalAfter = 1;      // the decimal marker sits after this many digits (1..18)
	bool Empty = true;         // nothing stored (all wheels at 0, marker parked)
	bool Special = false;      // inf or nan: the wheels show 9s with the sign wheel half-turned
};

// Fixed notation when the point falls inside the 18 digits, otherwise one digit before the point
// and the power of ten on the exponent wheels.
WheelReading ReadingFor(const mcalc::Real& X);

struct MemoryColumn
{
	std::string Label;         // e.g. "ans", "h.dv1", "A(1,1) 2x2"
	WheelReading Reading;
};

// Up to Count columns: ans first, then variables in name order. A record shows its first field and a
// matrix its first element, with the shape in the label.
std::vector<MemoryColumn> ReadMemory(const mcalc::Session& S, int Count);
} // namespace mkb
