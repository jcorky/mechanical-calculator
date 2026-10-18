// The mechanism player: recorded machine events become continuous motion that ends exactly on
// the machine's registers, with carries staggered up the wheels and the bell on each overdraft.
#include "GearMachine.h"
#include "GearReal.h"
#include "MechanismPlayer.h"
#include "TestFramework.h"

#include <cmath>

TEST(mechanism_plays_a_multiplication_to_its_registers)
{
	mcalc::GearRecorder Rec;
	mcalc::GearProcedures::Multiply(999, 99, &Rec);
	mkb::MechanismPlayer P;
	P.MaxSeconds = 100.0;
	P.Play(Rec.GetEvents());
	CHECK(P.IsPlaying());
	const double Expected = 18 * (1.0 / P.TurnsPerSecond) + P.ClearSeconds + P.SetSeconds;   // 9 + 9 turns
	CHECK(P.Duration() > Expected * 0.9);
	double MaxLever = 0.0;
	double LastCrank = 0.0;
	bool Monotone = true;
	for (int i = 0; i < 4000 && P.IsPlaying(); ++i)
	{
		P.Step(1.0 / 240.0);
		for (double L : P.State().CarryLever)
		{
			MaxLever = std::max(MaxLever, L);
		}
		Monotone = Monotone && P.State().CrankAngle >= LastCrank - 1e-12;
		LastCrank = P.State().CrankAngle;
	}
	CHECK(!P.IsPlaying());
	CHECK(Monotone);
	CHECK(MaxLever > 0.5);   // 999 * 99 carries
	// 999 * 99 = 98901, on the accumulator wheels, units first.
	const int Want[5] = {1, 0, 9, 8, 9};
	for (int i = 0; i < 5; ++i)
	{
		CHECK_NEAR(P.State().Accumulator[i], Want[i], 1e-12);
	}
	CHECK_NEAR(P.State().Counter[0], 9, 1e-12);
	CHECK_NEAR(P.State().Counter[1], 9, 1e-12);
	CHECK_NEAR(P.State().CrankAngle, 18 * 6.283185307179586, 1e-9);
}

TEST(mechanism_rings_the_bell_and_compresses_time)
{
	mcalc::GearRecorder Rec;
	bool Smaller = false;
	mcalc::GearProcedures::Divide(355000000000000000ull, 113000000000000000ull, &Rec, &Smaller);
	mkb::MechanismPlayer P;
	P.MaxSeconds = 2.0;
	P.Play(Rec.GetEvents());
	CHECK(P.Duration() > 2.0);   // compressed to MaxSeconds
	int Bells = 0;
	double T = 0.0;
	while (P.IsPlaying() && T < 10.0)
	{
		Bells += P.Step(1.0 / 60.0);
		T += 1.0 / 60.0;
	}
	CHECK(Bells == 20);
	CHECK(T <= 2.0 + 1.0 / 30.0);
}

TEST(mechanism_skip_and_idle_setting)
{
	mcalc::GearRecorder Rec;
	mcalc::GearProcedures::AddAligned(123, 2, 45, false, &Rec);
	mkb::MechanismPlayer P;
	P.Play(Rec.GetEvents());
	P.Step(0.01);
	P.SkipToEnd();
	CHECK(!P.IsPlaying());
	CHECK_NEAR(P.State().Accumulator[0], 5, 0);
	CHECK_NEAR(P.State().Accumulator[2], 3, 0);   // 12300 + 45
	const std::vector<uint8_t> D = mkb::SettingFromInput("hohmann(6678e3, 42", 18);
	CHECK(D.size() == 2 && D[0] == 2 && D[1] == 4);
	CHECK(mkb::SettingFromInput("x2", 2).empty());
	CHECK(mkb::SettingFromInput("3.25", 4).size() == 3);
	P.ShowSetting(D);
	for (int i = 0; i < 200; ++i)
	{
		P.Step(1.0 / 60.0);
	}
	CHECK_NEAR(P.State().Setting[0], 2, 1e-3);
	CHECK_NEAR(P.State().Setting[1], 4, 1e-3);
}
