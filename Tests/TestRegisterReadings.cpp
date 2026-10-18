// Figure-wheel readings: fixed or exponent notation, the decimal marker, and the memory columns.
#include "RegisterReadings.h"
#include "TestFramework.h"

namespace
{
std::string DigitsOf(const mkb::WheelReading& R)
{
	std::string S;
	for (uint8_t D : R.Digits)
	{
		S += static_cast<char>('0' + D);
	}
	return S;
}
} // namespace

TEST(register_readings_fixed_and_exponent)
{
	const mkb::WheelReading A = mkb::ReadingFor(mcalc::Real(2425.769028));
	CHECK(DigitsOf(A) == "242576902800000000" && A.DecimalAfter == 4 && A.Exponent == 0 && !A.Negative);
	const mkb::WheelReading B = mkb::ReadingFor(mcalc::Real(-1.380649e-23));
	CHECK(B.Negative && B.ExponentNegative && B.Exponent == -23 && B.DecimalAfter == 1);
	CHECK(B.ExponentDigits[0] == 0 && B.ExponentDigits[1] == 2 && B.ExponentDigits[2] == 3);
	const mkb::WheelReading C = mkb::ReadingFor(mcalc::Real(0));
	CHECK(!C.Empty && DigitsOf(C) == "000000000000000000");
	CHECK(mkb::ReadingFor(mcalc::Real::Infinity()).Special);
	const mkb::WheelReading Big = mkb::ReadingFor(mcalc::Real(6678e3) * mcalc::Real(1e20));
	CHECK(Big.Exponent == 26 && Big.DecimalAfter == 1);
}

TEST(register_readings_memory_columns)
{
	mcalc::Session S;
	S.Execute("x = 7; A = [1, 2; 3, 4]; s = size(A)");
	S.Execute("41 + 1");
	const std::vector<mkb::MemoryColumn> M = mkb::ReadMemory(S, 6);
	CHECK(M.size() == 6);
	CHECK(M[0].Label == "ans" && DigitsOf(M[0].Reading).compare(0, 2, "42") == 0);
	CHECK(M[1].Label == "A(1,1) 2x2");
	CHECK(M[2].Label == "s.rows");
	CHECK(M[3].Label == "x");
	CHECK(M[4].Reading.Empty && M[5].Label.empty());
}
