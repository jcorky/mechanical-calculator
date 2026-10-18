#include "RegisterReadings.h"

#include <cstdlib>

namespace mkb
{
WheelReading ReadingFor(const mcalc::Real& X)
{
	WheelReading R;
	if (X.IsNaN() || X.IsInf())
	{
		R.Empty = false;
		R.Special = true;
		R.Negative = X.IsNegative();
		for (uint8_t& D : R.Digits)
		{
			D = 9;
		}
		return R;
	}
	if (X.IsZero())
	{
		R.Empty = false;
		R.DecimalAfter = 1;
		return R;
	}
	R.Empty = false;
	R.Negative = X.IsNegative();
	uint64_t M = X.Mantissa();
	for (int i = 17; i >= 0; --i)
	{
		R.Digits[i] = static_cast<uint8_t>(M % 10);
		M /= 10;
	}
	const int E = X.Exponent();
	if (E >= 0 && E < 18)
	{
		R.DecimalAfter = E + 1;
		R.Exponent = 0;
	}
	else
	{
		R.DecimalAfter = 1;
		R.Exponent = E;
		R.ExponentNegative = E < 0;
		int A = std::abs(E);
		if (A > 999)
		{
			A = 999;   // three wheels; beyond that the wheels stop at 999
		}
		R.ExponentDigits[0] = static_cast<uint8_t>(A / 100);
		R.ExponentDigits[1] = static_cast<uint8_t>(A / 10 % 10);
		R.ExponentDigits[2] = static_cast<uint8_t>(A % 10);
	}
	return R;
}

std::vector<MemoryColumn> ReadMemory(const mcalc::Session& S, int Count)
{
	std::vector<MemoryColumn> Out;
	auto Add = [&Out](const std::string& Name, const mcalc::Value& V) {
		MemoryColumn C;
		C.Label = Name;
		if (V.IsRecord() && !V.M.V.empty())
		{
			C.Label = Name + "." + V.Fields[0];
			C.Reading = ReadingFor(V.M.V[0]);
		}
		else if (V.M.IsScalar())
		{
			C.Reading = ReadingFor(V.Scalar());
		}
		else if (!V.M.V.empty())
		{
			C.Label = Name + "(1,1) " + std::to_string(V.M.Rows) + "x" + std::to_string(V.M.Cols);
			C.Reading = ReadingFor(V.M.V[0]);
		}
		else
		{
			C.Label = Name + " []";
		}
		Out.push_back(C);
	};
	mcalc::Value Ans;
	if (S.GetVariable("ans", &Ans))
	{
		Add("ans", Ans);
	}
	for (const auto& KV : S.Variables())
	{
		if (static_cast<int>(Out.size()) >= Count)
		{
			break;
		}
		if (KV.first != "ans")
		{
			Add(KV.first, KV.second);
		}
	}
	while (static_cast<int>(Out.size()) < Count)
	{
		Out.push_back(MemoryColumn());
	}
	if (static_cast<int>(Out.size()) > Count)
	{
		Out.resize(static_cast<size_t>(Count));
	}
	return Out;
}
} // namespace mkb
