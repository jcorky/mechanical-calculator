// The animator: the player's registers turn the right digits to the reading windows, the pinions
// slide to their digits, the carriage carries its wheels, the reverser turns the drums backward,
// and the memory columns, type wheels and paper reach their targets.
#include "CalculatorDevice.h"
#include "GearMachine.h"
#include "GearParts.h"
#include "MechanismAnimator.h"
#include "MechanismPlayer.h"
#include "RegisterReadings.h"
#include "TestFramework.h"

#include <cmath>

namespace
{
const double kPi = 3.14159265358979323846;

const mgeo::MechanismParts& AnimatedMachine()
{
	static const mgeo::MechanismParts M = mgeo::BuildMechanism(mkb::FullSizeAnsi(), mgeo::CaseParams());
	return M;
}

// The direction from a wheel's axis to digit D's glyph, in the wheel's unturned frame.
mgeo::Vec3 GlyphDirection(int D)
{
	const double A = 0.5 * kPi - D * kPi / 5.0;
	return mgeo::Vec3(std::cos(A), std::sin(A), 0.0);
}

bool Faces(const mgeo::Pose& Wheel, int Digit, const mgeo::ReadingWindow& Window)
{
	return Wheel.Direction(GlyphDirection(Digit)).Dot(Window.Normal) > 0.98;
}

double Distance(const mgeo::Vec3& A, const mgeo::Vec3& B)
{
	return (A - B).Length();
}
} // namespace

TEST(animator_covers_every_moving_part)
{
	const mgeo::MechanismParts& M = AnimatedMachine();
	mkb::MechanismAnimator A(M);
	for (int I = 0; I < mkb::kMovingParts; ++I)
	{
		const mkb::MovingPart P = static_cast<mkb::MovingPart>(I);
		const size_t Want = mkb::PartSetOf(M, P).Poses.size();
		CHECK(Want > 0);
		CHECK(A.Poses(P).size() == (P == mkb::MovingPart::Carriage ? 1 : Want));
		CHECK(!mkb::PartSetOf(M, P).Model.Sections.empty());
	}
	CHECK(mkb::RidesCarriage(M, mkb::MovingPart::AccumulatorWheels));
	CHECK(mkb::RidesCarriage(M, mkb::MovingPart::Carriage));
	CHECK(!mkb::RidesCarriage(M, mkb::MovingPart::Drums));
}

TEST(animator_shows_the_product_in_the_windows)
{
	const mgeo::MechanismParts& M = AnimatedMachine();
	mcalc::GearRecorder Rec;
	mcalc::GearProcedures::Multiply(271828, 315, &Rec);
	mkb::MechanismPlayer P;
	P.MaxSeconds = 100.0;
	P.Play(Rec.GetEvents());
	mkb::MechanismAnimator A(M);
	for (int I = 0; I < 20000 && P.IsPlaying(); ++I)
	{
		P.Step(1.0 / 120.0);
		A.Follow(P.State(), 0.0);
	}
	CHECK(!P.IsPlaying());
	// 271828 * 315 = 85625820: accumulator wheels units first, the counter shows 315.
	const int Product[8] = {0, 2, 8, 5, 2, 6, 5, 8};
	const std::vector<mgeo::Pose>& Acc = A.Poses(mkb::MovingPart::AccumulatorWheels);
	for (int I = 0; I < 12; ++I)
	{
		const int D = I < 8 ? Product[I] : 0;
		CHECK(Faces(Acc[static_cast<size_t>(I)], D, M.AccumulatorWindows[static_cast<size_t>(I)]));
	}
	const int Multiplier[3] = {5, 1, 3};
	const std::vector<mgeo::Pose>& Cnt = A.Poses(mkb::MovingPart::CounterWheels);
	for (int I = 0; I < 3; ++I)
	{
		CHECK(Faces(Cnt[static_cast<size_t>(I)], Multiplier[I], M.CounterWindows[static_cast<size_t>(I)]));
	}
	// The multiplicand stays set: each setting wheel shows its digit, each pinion sits at it.
	const int Set[6] = {8, 2, 8, 1, 7, 2};
	for (int J = 0; J < 6; ++J)
	{
		const size_t Jz = static_cast<size_t>(J);
		CHECK(Faces(A.Poses(mkb::MovingPart::SettingWheels)[Jz], Set[J], M.SettingWindows[Jz]));
		const mgeo::Pose& Home = M.Pinions.Poses[Jz];
		CHECK(Distance(A.Poses(mkb::MovingPart::Pinions)[Jz].Origin, Home.Origin + Home.Z * (Set[J] * M.PinionStepMm)) < 1e-9);
	}
	// The carriage follows the player's shift, and the turns were all additions.
	const mgeo::Pose& Home = M.Carriage.Poses[0];
	CHECK(Distance(A.CarriageFrame().Origin, Home.Origin + Home.X * (P.State().Carriage * M.DigitPitchMm)) < 1e-9);
	CHECK(std::fabs(A.MainShaftAngle() - P.State().CrankAngle) < 1e-9);
	CHECK(std::fabs(P.State().CrankAngle - 9 * 2.0 * kPi) < 1e-9);
	// Pinion 0 (digit 8) made 8/12 of a turn per crank turn.
	CHECK(std::fabs(A.PinionAngle(0) - 9 * 2.0 * kPi * 8.0 / 12.0) < 1e-6);
}

TEST(animator_reverses_after_the_reverser)
{
	const mgeo::MechanismParts& M = AnimatedMachine();
	mkb::MechanismAnimator A(M);
	mkb::MechanismState S;
	A.Follow(S, 0.0);
	S.Reverse = 1.0;
	S.Setting[0] = 6.0;
	S.CrankAngle = 2.0 * kPi;
	A.Follow(S, 0.0);
	CHECK(std::fabs(A.MainShaftAngle() + 2.0 * kPi) < 1e-12);
	CHECK(std::fabs(A.PinionAngle(0) + 2.0 * kPi * 0.5) < 1e-12);
	// The crank itself always turns forward; the motor's extra spin adds to it.
	A.Follow(S, kPi);
	const mgeo::Pose Expect = M.Crank.Poses[0].RotatedZ(3.0 * kPi);
	CHECK(Distance(A.Poses(mkb::MovingPart::Crank)[0].X, Expect.X) < 1e-9);
	CHECK(std::fabs(A.MainShaftAngle() + 3.0 * kPi) < 1e-12);
	const mgeo::Pose Drum = M.Drums.Poses[0].RotatedZ(-3.0 * kPi);
	CHECK(Distance(A.Poses(mkb::MovingPart::Drums)[0].X, Drum.X) < 1e-9);
	// The reversing lever sits at its subtract angle.
	const mgeo::Pose Lever = M.ReversingLever.Poses[0].RotatedZ(M.ReverseSubtractAngle);
	CHECK(Distance(A.Poses(mkb::MovingPart::ReversingLever)[0].X, Lever.X) < 1e-9);
}

TEST(animator_spins_the_memory_columns_to_ans_and_the_variables)
{
	const mgeo::MechanismParts& M = AnimatedMachine();
	mcalc::Session S;
	S.Execute("x = 2.5e-30");
	S.Execute("-1234.5");
	const std::vector<mkb::MemoryColumn> Columns = mkb::ReadMemory(S, M.MemoryColumns);
	CHECK(Columns.size() >= 2);
	mkb::MechanismAnimator A(M);
	A.ShowMemory(Columns);
	A.ShowTypeWheels(Columns[0].Reading);
	CHECK(!A.Settled());
	double T = 0.0;
	while (!A.Settled() && T < 5.0)
	{
		A.Step(1.0 / 60.0);
		T += 1.0 / 60.0;
	}
	CHECK(A.Settled());
	CHECK(T < 1.0);   // at most ten digits at sixteen a second
	const int PerColumn = M.MemoryMantissa + M.MemoryExponent;
	for (size_t C = 0; C < Columns.size() && C < 2; ++C)
	{
		const mkb::WheelReading& R = Columns[C].Reading;
		for (int K = 0; K < M.MemoryMantissa; ++K)
		{
			CHECK(std::fabs(A.MemoryWheelPosition(static_cast<int>(C), K) - R.Digits[K]) < 1e-12);
			const size_t I = C * static_cast<size_t>(PerColumn) + static_cast<size_t>(K);
			CHECK(Faces(A.Poses(mkb::MovingPart::MemoryWheels)[I], R.Digits[K], M.MemoryWindows[I]));
		}
		for (int K = 0; K < M.MemoryExponent; ++K)
		{
			CHECK(std::fabs(A.MemoryWheelPosition(static_cast<int>(C), M.MemoryMantissa + K) - R.ExponentDigits[K]) < 1e-12);
		}
	}
	// ans = -1234.5: the sign wheel turned one place to '-'.
	CHECK(Columns[0].Reading.Negative);
	const mgeo::Pose Minus = M.MemorySignWheels.Poses[0].RotatedZ(mgeo::DigitAngleRad(1));
	CHECK(Distance(A.Poses(mkb::MovingPart::MemorySignWheels)[0].X, Minus.X) < 1e-9);
	for (int K = 0; K < 18; ++K)
	{
		CHECK(std::fabs(A.TypeWheelPosition(K) - Columns[0].Reading.Digits[K]) < 1e-12);
	}
	// Clearing the memory turns every wheel on, forward, to zero.
	A.ShowMemory({});
	while (!A.Settled())
	{
		A.Step(1.0 / 60.0);
	}
	CHECK(A.MemoryWheelPosition(0, 0) == 0.0);
}

TEST(animator_feeds_the_paper_line_by_line)
{
	const mgeo::MechanismParts& M = AnimatedMachine();
	mkb::MechanismAnimator A(M);
	const mgeo::Vec3 Line1 = A.TapeLine(1).Origin;
	CHECK(Distance(Line1, mgeo::TapeLineFrame(M, 1.5 * M.TapeLinePitchMm).Origin) < 1e-9);
	A.Feed(2);
	CHECK(std::fabs(A.PaperLag() - 2.0) < 1e-12);
	// Unfed, the old line 0 (now line 2) has not moved.
	CHECK(Distance(A.TapeLine(2).Origin, mgeo::TapeLineFrame(M, 0.5 * M.TapeLinePitchMm).Origin) < 1e-9);
	const mgeo::Pose Roll0 = A.Poses(mkb::MovingPart::Platen)[0];
	for (int I = 0; I < 120 && A.PaperLag() > 0.0; ++I)
	{
		A.Step(1.0 / 60.0);
	}
	CHECK(A.PaperLag() == 0.0);
	CHECK(Distance(A.TapeLine(2).Origin, mgeo::TapeLineFrame(M, 2.5 * M.TapeLinePitchMm).Origin) < 1e-9);
	// The platen turned by two pitches of its circumference.
	const mgeo::Pose Turned = Roll0.RotatedZ(2.0 * M.TapeLinePitchMm / M.PlatenRadius);
	CHECK(Distance(A.Poses(mkb::MovingPart::Platen)[0].X, Turned.X) < 1e-9);
}

TEST(animator_moves_the_decimal_markers)
{
	const mgeo::MechanismParts& M = AnimatedMachine();
	mkb::MechanismAnimator A(M);
	CHECK(mkb::DecimalPlacesInInput("3.1415", 6) == 4);
	CHECK(mkb::DecimalPlacesInInput("sin(2.50)", 7) == 2);
	CHECK(mkb::DecimalPlacesInInput("x2.5", 4) == 0);
	CHECK(mkb::DecimalPlacesInInput("12+3", 2) == 0);
	CHECK(mkb::DecimalPlacesInInput("", 0) == 0);
	A.SetDecimalMarkers(4, 6, 0);
	while (!A.Settled())
	{
		A.Step(1.0 / 60.0);
	}
	const mgeo::Pose& Home = M.SettingMarker.Poses[0];
	CHECK(Distance(A.Poses(mkb::MovingPart::SettingMarker)[0].Origin, Home.Origin + Home.X * (4 * M.DigitPitchMm)) < 1e-9);
	const mgeo::Pose& AccHome = M.CarriageMarkers.Poses[0];
	CHECK(Distance(A.Poses(mkb::MovingPart::CarriageMarkers)[0].Origin, AccHome.Origin + AccHome.X * (6 * M.DigitPitchMm)) < 1e-9);
}

TEST(device_counts_printed_lines)
{
	mkb::CalculatorDevice D;
	const uint64_t Before = D.PrintedCount();
	D.TypeText("1+2");
	D.Enter();
	CHECK(D.PrintedCount() >= Before + 2);   // the line and its result
}
