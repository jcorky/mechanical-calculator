// The calculating mechanism: closed parts, true involutes, meshing gears at the right distances, the
// stepped drums' digit engagement, the shifting carriage's alignment, the numeral wheels' reading
// convention, and the machine clear of the keyboard.
#include "GearParts.h"
#include "TestFramework.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <tuple>

namespace geartest_detail
{
using mgeo::Mesh;
using mgeo::Pose;
using mgeo::Vec3;

constexpr double kPi = 3.14159265358979323846;

int Shells(const Mesh& M)
{
	std::map<std::tuple<long long, long long, long long>, int> Ids;
	std::vector<int> Weld(M.Positions.size());
	for (size_t i = 0; i < M.Positions.size(); ++i)
	{
		const Vec3& P = M.Positions[i];
		Weld[i] = Ids.emplace(std::make_tuple(std::llround(P.X * 1e5), std::llround(P.Y * 1e5), std::llround(P.Z * 1e5)), static_cast<int>(Ids.size())).first->second;
	}
	std::vector<int> Parent(Ids.size());
	std::iota(Parent.begin(), Parent.end(), 0);
	auto Find = [&Parent](int X) {
		while (Parent[static_cast<size_t>(X)] != X)
		{
			X = Parent[static_cast<size_t>(X)] = Parent[static_cast<size_t>(Parent[static_cast<size_t>(X)])];
		}
		return X;
	};
	std::vector<char> Used(Ids.size(), 0);
	for (size_t t = 0; t + 2 < M.Triangles.size(); t += 3)
	{
		const int A = Weld[static_cast<size_t>(M.Triangles[t])], B = Weld[static_cast<size_t>(M.Triangles[t + 1])], C = Weld[static_cast<size_t>(M.Triangles[t + 2])];
		Used[static_cast<size_t>(A)] = Used[static_cast<size_t>(B)] = Used[static_cast<size_t>(C)] = 1;
		Parent[static_cast<size_t>(Find(A))] = Find(B);
		Parent[static_cast<size_t>(Find(B))] = Find(C);
	}
	int N = 0;
	for (size_t i = 0; i < Parent.size(); ++i)
	{
		N += Used[i] && Find(static_cast<int>(i)) == static_cast<int>(i) ? 1 : 0;
	}
	return N;
}

// Closed, outward, normals agreeing with the winding; Genus < 0 skips the topology check (parts made
// of many overlapping shells with through-holes), otherwise chi = 2 * shells - 2 * genus.
bool CheckClosed(const std::string& Name, const Mesh& M, int Genus = -1)
{
	const bool Closed = M.IsClosed();
	const double V = M.SignedVolume();
	const bool Normals = M.NormalsAgreeWithWinding(0.3);
	const int S = Shells(M), Chi = M.EulerCharacteristic();
	const bool Topo = Genus < 0 ? (Chi % 2 == 0 && Chi <= 2 * S) : Chi == 2 * S - 2 * Genus;
	CHECK(Closed);
	CHECK(V > 0.0);
	CHECK(Normals);
	CHECK(Topo);
	if (!Closed || V <= 0.0 || !Normals || !Topo)
	{
		std::printf("    %s: closed %d volume %.2f normals %d shells %d chi %d\n", Name.c_str(), Closed, V, Normals, S, Chi);
		return false;
	}
	return true;
}

double RayTri(const Vec3& O, const Vec3& D, const Vec3& A, const Vec3& B, const Vec3& C)
{
	const Vec3 E1 = B - A, E2 = C - A, Pv = D.Cross(E2);
	const double Det = E1.Dot(Pv);
	if (std::fabs(Det) < 1e-14)
	{
		return -1.0;
	}
	const Vec3 Tv = O - A;
	const double U = Tv.Dot(Pv) / Det;
	if (U < 0.0 || U > 1.0)
	{
		return -1.0;
	}
	const Vec3 Qv = Tv.Cross(E1);
	const double V = D.Dot(Qv) / Det;
	if (V < 0.0 || U + V > 1.0)
	{
		return -1.0;
	}
	const double T = E2.Dot(Qv) / Det;
	return T > 1e-9 ? T : -1.0;
}

double FirstHit(const Mesh& M, const Vec3& O, const Vec3& D)
{
	double Best = -1.0;
	for (size_t t = 0; t + 2 < M.Triangles.size(); t += 3)
	{
		const double T = RayTri(O, D, M.Positions[static_cast<size_t>(M.Triangles[t])], M.Positions[static_cast<size_t>(M.Triangles[t + 1])], M.Positions[static_cast<size_t>(M.Triangles[t + 2])]);
		if (T > 0.0 && (Best < 0.0 || T < Best))
		{
			Best = T;
		}
	}
	return Best;
}

bool Inside(const Mesh& M, const Vec3& P)
{
	const Vec3 D = Vec3(0.1234, 0.3011, 0.9452).Normalized();
	int N = 0;
	for (size_t t = 0; t + 2 < M.Triangles.size(); t += 3)
	{
		N += RayTri(P, D, M.Positions[static_cast<size_t>(M.Triangles[t])], M.Positions[static_cast<size_t>(M.Triangles[t + 1])], M.Positions[static_cast<size_t>(M.Triangles[t + 2])]) > 0.0 ? 1 : 0;
	}
	return (N & 1) == 1;
}

Mesh Placed(const Mesh& M, const Pose& P)
{
	Mesh Out = M;
	for (Vec3& Q : Out.Positions)
	{
		Q = P.Point(Q);
	}
	for (Vec3& N : Out.Normals)
	{
		N = P.Direction(N);
	}
	return Out;
}

Vec3 ToLocal(const Pose& P, const Vec3& World)
{
	const Vec3 D = World - P.Origin;
	return Vec3(D.Dot(P.X), D.Dot(P.Y), D.Dot(P.Z));
}

struct Box3
{
	Vec3 Lo = Vec3(1e30, 1e30, 1e30), Hi = Vec3(-1e30, -1e30, -1e30);
	void Add(const Vec3& P)
	{
		Lo = Vec3(std::min(Lo.X, P.X), std::min(Lo.Y, P.Y), std::min(Lo.Z, P.Z));
		Hi = Vec3(std::max(Hi.X, P.X), std::max(Hi.Y, P.Y), std::max(Hi.Z, P.Z));
	}
	void Add(const Mesh& M, const Pose& P)
	{
		for (const Vec3& Q : M.Positions)
		{
			Add(P.Point(Q));
		}
	}
	bool Overlaps(const Box3& O) const { return Lo.X < O.Hi.X && O.Lo.X < Hi.X && Lo.Y < O.Hi.Y && O.Lo.Y < Hi.Y && Lo.Z < O.Hi.Z && O.Lo.Z < Hi.Z; }
};

const Mesh* Section(const mgeo::MechPart& P, mgeo::MechMaterial M)
{
	for (const mgeo::MechSection& S : P.Sections)
	{
		if (S.Material == M)
		{
			return &S.Geometry;
		}
	}
	return nullptr;
}

const mgeo::MechanismParts& Machine()
{
	static const mgeo::MechanismParts M = mgeo::BuildMechanism(mkb::FullSizeAnsi(), mgeo::CaseParams());
	return M;
}

Pose CarriageAt(const mgeo::MechanismParts& M, int Shift)
{
	return M.Carriage.Poses[0].Moved(Vec3(Shift * M.DigitPitchMm, 0, 0));
}

// The world pose of an instance (carriage-frame instances at the given shift).
Pose WorldPose(const mgeo::MechanismParts& M, const mgeo::PartSet& S, size_t I, int Shift)
{
	return S.OnCarriage ? CarriageAt(M, Shift).Then(S.Poses[I]) : S.Poses[I];
}

// Measured tooth half-angle (about the tooth centre) of a gear mesh at radius R, slicing at z = 0.
double ToothHalfAngle(const Mesh& G, double R)
{
	// Walk the circle of radius R from the tooth-0 centre (+X) until we leave the material.
	const double Step = 1e-4;
	double A = 0.0;
	while (A < kPi && Inside(G, Vec3(R * std::cos(A), R * std::sin(A), 0.01)))
	{
		A += Step;
	}
	return A;
}

double Involute(double A)
{
	return std::tan(A) - A;
}
} // namespace geartest_detail

using namespace geartest_detail;

TEST(gears_primitives_closed)
{
	mgeo::SpurGearParams G;
	G.Teeth = 20;
	G.Module = 1.0;
	CheckClosed("spur solid", mgeo::BuildSpurGear(G), 0);
	G.BoreRadius = 2.0;
	CheckClosed("spur bored", mgeo::BuildSpurGear(G), 1);
	G.HubRadius = 5.0;
	G.HubLength = 3.0;
	CheckClosed("spur hub", mgeo::BuildSpurGear(G), 1);
	G.Teeth = 12;
	G.Module = 0.6;
	G.BoreRadius = 1.35;
	G.SquareBore = true;
	G.HubRadius = 2.4;
	G.HubLength = 1.6;
	CheckClosed("pinion", mgeo::BuildSpurGear(G), 1);
	mgeo::BevelGearParams B;
	CheckClosed("bevel solid", mgeo::BuildBevelGear(B), 0);
	B.BoreRadius = 1.5;
	B.HubRadius = 4.0;
	B.HubLength = 2.0;
	CheckClosed("bevel bored", mgeo::BuildBevelGear(B), 1);
	const Mesh Drum = mgeo::BuildStepDrum(mgeo::StepDrumParams());
	CheckClosed("drum", Drum, 0);
	CHECK(Shells(Drum) == 11);   // core, nine teeth, bevel
	mgeo::NumeralWheelParams W;
	W.GearTeeth = 15;
	W.GearGap = 1.9;
	W.StarTeeth = 10;
	W.CarryTooth = true;
	const mgeo::NumeralWheel N = mgeo::BuildNumeralWheel(W);
	CheckClosed("numeral rim", N.Wheel, 1);
	CheckClosed("numeral gears", N.Gears);
	CheckClosed("numeral glyphs", N.Glyphs, 0);
	CHECK(N.GlyphTriangleStart.size() == 11);
	mgeo::SpokedWheelParams S;
	CheckClosed("spoked gear", mgeo::BuildSpokedWheel(S));
	CHECK(Shells(mgeo::BuildSpokedWheel(S)) == 2 + S.Spokes);
	S.Teeth = 0;
	S.RimRadius = 30.0;
	CheckClosed("flywheel", mgeo::BuildSpokedWheel(S));
	CheckClosed("shaft", mgeo::BuildShaft(2.0, 30.0), 0);
	CheckClosed("square shaft", mgeo::BuildShaft(1.3, 30.0, true), 0);
	CheckClosed("bearing block", mgeo::BuildBearingBlock(2.0, 16.0, 4.0));
	CheckClosed("rack", mgeo::BuildRack(1.0, 20, 4.0, 3.0), 0);
	CheckClosed("digit text", mgeo::BuildDigitText("0123456789+-", 5.0, 0.6, 0.3), 0);
}

TEST(gears_involute_profile)
{
	for (int Teeth : {12, 20, 40})
	{
		mgeo::SpurGearParams G;
		G.Teeth = Teeth;
		G.Module = 1.0;
		G.FlankPoints = 12;
		const Mesh M = mgeo::BuildSpurGear(G);
		const mgeo::GearCircles C = mgeo::SpurGearCircles(G);
		CHECK_NEAR(C.Pitch, 0.5 * Teeth, 1e-12);
		CHECK_NEAR(C.Base, 0.5 * Teeth * std::cos(20.0 * kPi / 180.0), 1e-12);
		// Tooth thickness on the pitch circle is pi m / 2 (the arc between the two flanks).
		const double Half = ToothHalfAngle(M, C.Pitch);
		CHECK_NEAR(2.0 * Half * C.Pitch, 0.5 * kPi * G.Module, 0.02);
		// The flanks are involutes of the base circle: between two radii the half-angle changes by
		// the difference of the involute function.
		const double R1 = std::max(C.Base, C.Root) + 0.2, R2 = C.Tip - 0.2;
		const double Measured = ToothHalfAngle(M, R1) - ToothHalfAngle(M, R2);
		const double Expected = Involute(std::acos(C.Base / R2)) - Involute(std::acos(C.Base / R1));
		CHECK_NEAR(Measured, Expected, 0.004);
	}
}

TEST(gears_meshing_pairs)
{
	const mgeo::MechanismParts& M = Machine();
	int Spur = 0, Bevel = 0, Rack = 0;
	for (int Shift : {0, 7, 20})
	{
		for (const mgeo::GearPair& G : M.GearPairs)
		{
			Vec3 B = G.PointB;
			int Offset = 0;
			if (G.OnCarriageB)
			{
				// The wheel over column j at shift k is wheel j + k: its gear sits k pitches further left
				// in the carriage frame, then the carriage moves k pitches right.
				Offset = G.Type == mgeo::GearPair::Rack ? 0 : Shift;
				B = CarriageAt(M, Shift).Point(G.PointB - Vec3(Offset * M.DigitPitchMm, 0, 0));
			}
			else if (Shift != 0)
			{
				continue;
			}
			const Vec3 D = B - G.PointA;
			if (G.Type == mgeo::GearPair::Spur)
			{
				CHECK(std::fabs(G.AxisA.Cross(G.AxisB).Length()) < 1e-9);
				const double Dist = (D - G.AxisA * D.Dot(G.AxisA)).Length();
				CHECK_NEAR(Dist, 0.5 * G.Module * (G.TeethA + G.TeethB), 1e-6);
				++Spur;
			}
			else if (G.Type == mgeo::GearPair::Bevel)
			{
				CHECK(D.Length() < 1e-9);
				CHECK(std::fabs(G.AxisA.Dot(G.AxisB)) < 1e-9);
				++Bevel;
			}
			else
			{
				// Rack: the pinion's axis is the pitch radius from the rack's pitch line, square to it.
				const Vec3 Along = G.AxisB.Normalized();
				const Vec3 Off = D - Along * D.Dot(Along);
				CHECK_NEAR(Off.Length(), 0.5 * G.Module * G.TeethA, 1e-6);
				CHECK(std::fabs(G.AxisA.Dot(Along)) < 1e-9);
				++Rack;
			}
		}
	}
	CHECK(Spur >= 20 + 20 * 3 && Bevel >= 40 + 3 && Rack == 3);
	CHECK_NEAR(M.FlywheelRatio, 60.0 / 20.0, 1e-12);
	// The parts sit where the pairs say: pinion over its drum, wheels over their transfer gears.
	for (size_t J = 0; J < M.Drums.Poses.size(); ++J)
	{
		const Vec3 Rel = M.Pinions.Poses[J].Origin - M.Drums.Poses[J].Origin;
		CHECK_NEAR(std::sqrt(Rel.X * Rel.X + Rel.Z * Rel.Z), 0.6 * (20 + 12) / 2.0, 1e-9);
		const Vec3 Tw = WorldPose(M, M.AccumulatorWheels, J, 0).Origin - M.TransferGears.Poses[J].Origin;
		CHECK_NEAR(std::sqrt(Tw.Y * Tw.Y + Tw.Z * Tw.Z), 0.5 * (18 + 15) / 2.0, 1e-9);
	}
}

TEST(gears_pinion_engages_digit_teeth)
{
	const mgeo::MechanismParts& M = Machine();
	const Mesh* Drum = Section(M.Drums.Model, mgeo::MechMaterial::Brass);
	CHECK(Drum != nullptr);
	if (!Drum)
	{
		return;
	}
	const mgeo::StepDrumParams P;
	const mgeo::StepDrumLayout L = mgeo::StepDrumGeometry(P);
	for (size_t J : {size_t(0), size_t(11), size_t(19)})
	{
		for (int D = 0; D <= 9; ++D)
		{
			const Pose Pin = M.Pinions.Poses[J].Moved(Vec3(0, 0, D * M.PinionStepMm));
			const double Zc = ToLocal(M.Drums.Poses[J], Pin.Origin).Z;
			CHECK_NEAR(Zc, L.PinionCenterZ(D), 1e-9);
			int Engaged = 0;
			for (int K = 1; K <= 9; ++K)
			{
				const double A = L.ToothAngleRad(K), R = L.TipRadius - 0.2;
				bool All = true;
				for (double Dz : {-0.5 * P.PinionFaceWidth + 0.05, 0.0, 0.5 * P.PinionFaceWidth - 0.05})
				{
					All = All && Inside(*Drum, Vec3(R * std::cos(A), R * std::sin(A), Zc + Dz));
				}
				Engaged += All ? 1 : 0;
			}
			CHECK(Engaged == D);
		}
	}
}

TEST(gears_carriage_alignment)
{
	const mgeo::MechanismParts& M = Machine();
	// Wheels are one drum pitch apart, so at shift k drum i lines up with wheel i + k.
	for (size_t I = 0; I + 1 < M.AccumulatorWheels.Poses.size(); ++I)
	{
		CHECK_NEAR(M.AccumulatorWheels.Poses[I].Origin.X - M.AccumulatorWheels.Poses[I + 1].Origin.X, M.DigitPitchMm, 1e-9);
	}
	for (size_t J = 0; J + 1 < M.Drums.Poses.size(); ++J)
	{
		CHECK_NEAR(M.Drums.Poses[J].Origin.X - M.Drums.Poses[J + 1].Origin.X, M.DigitPitchMm, 1e-9);
	}
	for (int K = 0; K <= M.MaxShift; ++K)
	{
		for (size_t I = 0; I < M.Drums.Poses.size(); ++I)
		{
			CHECK_NEAR(WorldPose(M, M.AccumulatorWheels, I + static_cast<size_t>(K), K).Origin.X, M.Drums.Poses[I].Origin.X, 1e-9);
		}
		// The counter wheel at the current shift sits at the counting finger.
		if (K < M.CounterDigits)
		{
			CHECK_NEAR(WorldPose(M, M.CounterWheels, static_cast<size_t>(K), K).Origin.X, M.Drums.Poses[0].Origin.X, 1e-9);
		}
	}
	CHECK(M.AccumulatorWheels.Poses.size() == 40 && M.CounterWheels.Poses.size() == 20 && M.Drums.Poses.size() == 20);
}

// Each numeral wheel's digit d, turned by DigitAngleRad(d), faces its window: the glyph's centre lies
// on the window's axis just behind the slot, and the glyphs sit on the wheel's rim.
TEST(gears_numeral_wheels_read_through_windows)
{
	const mgeo::MechanismParts& M = Machine();
	struct Reg
	{
		const char* Name;
		const mgeo::PartSet* Set;
		const std::vector<mgeo::ReadingWindow>* Windows;
		mgeo::MechMaterial Rim, Glyph;
	};
	const Reg Regs[] = {{"setting", &M.SettingWheels, &M.SettingWindows, mgeo::MechMaterial::NumeralBlack, mgeo::MechMaterial::Paint},
		{"counter", &M.CounterWheels, &M.CounterWindows, mgeo::MechMaterial::NumeralBlack, mgeo::MechMaterial::Paint},
		{"accumulator", &M.AccumulatorWheels, &M.AccumulatorWindows, mgeo::MechMaterial::NumeralBlack, mgeo::MechMaterial::Paint},
		{"memory", &M.MemoryWheels, &M.MemoryWindows, mgeo::MechMaterial::Brass, mgeo::MechMaterial::Engraving}};
	for (const Reg& R : Regs)
	{
		const Mesh* Rim = Section(R.Set->Model, R.Rim);
		const Mesh* Gl = Section(R.Set->Model, R.Glyph);
		CHECK(Rim && Gl);
		if (!Rim || !Gl)
		{
			continue;
		}
		double RimR = 0.0, GlMin = 1e9, GlMax = 0.0;
		for (const Vec3& Q : Rim->Positions)
		{
			RimR = std::max(RimR, std::hypot(Q.X, Q.Y));
		}
		for (const Vec3& Q : Gl->Positions)
		{
			GlMin = std::min(GlMin, std::hypot(Q.X, Q.Y));
			GlMax = std::max(GlMax, std::hypot(Q.X, Q.Y));
		}
		CHECK(GlMin > RimR - 0.2 && GlMax < RimR + 0.6);
		for (int D = 0; D < 10; ++D)
		{
			// Glyph d: the vertices within 17 degrees of its documented angle; their box centre (in angle
			// and along the axis) is the glyph's centre.
			const double A = 0.5 * kPi - D * 2.0 * kPi / 10.0;
			double A0 = 1e9, A1 = -1e9, Z0 = 1e9, Z1 = -1e9;
			int N = 0;
			for (const Vec3& Q : Gl->Positions)
			{
				const double Da = std::remainder(std::atan2(Q.Y, Q.X) - A, 2.0 * kPi);
				if (std::fabs(Da) < 17.0 * kPi / 180.0)
				{
					A0 = std::min(A0, Da);
					A1 = std::max(A1, Da);
					Z0 = std::min(Z0, Q.Z);
					Z1 = std::max(Z1, Q.Z);
					++N;
				}
			}
			CHECK(N > 20);
			const double Ac = A + 0.5 * (A0 + A1);
			const Vec3 C(RimR * std::cos(Ac), RimR * std::sin(Ac), 0.5 * (Z0 + Z1));
			for (size_t I = 0; I < R.Set->Poses.size(); I += 7)
			{
				const Pose W = R.Set->Poses[I].RotatedZ(mgeo::DigitAngleRad(D));
				const Vec3 G = W.Point(C);
				const mgeo::ReadingWindow& Win = (*R.Windows)[I];
				const Vec3 Off = Win.Center - G;
				const double Along = Off.Dot(Win.Normal);
				CHECK(Along > 0.2 && Along < 3.0);
				CHECK((Off - Win.Normal * Along).Length() < 0.8);
				// The wheel's reading direction is the window's normal.
				CHECK((G - W.Origin).Normalized().Dot(Win.Normal) > 0.98);
			}
		}
		if (R.Windows->size() != R.Set->Poses.size())
		{
			std::printf("    %s: %zu windows for %zu wheels\n", R.Name, R.Windows->size(), R.Set->Poses.size());
		}
		CHECK(R.Windows->size() >= R.Set->Poses.size());
	}
}

// The black window bars have an open slot in front of every wheel and solid bar beside it.
TEST(gears_window_slots_over_wheels)
{
	const mgeo::MechanismParts& M = Machine();
	const Mesh* Carriage = Section(M.Carriage.Model, mgeo::MechMaterial::BlackEnamel);
	const Mesh* Frame = Section(M.Frame, mgeo::MechMaterial::BlackEnamel);
	CHECK(Carriage && Frame);
	if (!Carriage || !Frame)
	{
		return;
	}
	auto Check = [](const Mesh& Bar, const std::vector<mgeo::ReadingWindow>& Ws) {
		int Bad = 0;
		for (const mgeo::ReadingWindow& W : Ws)
		{
			const Vec3 From = W.Center + W.Normal * 4.0;
			const double Open = FirstHit(Bar, From, -W.Normal);
			const Vec3 Right = W.Up.Cross(W.Normal).Normalized();
			const double Beside = FirstHit(Bar, From + Right * (0.5 * W.WidthMm + 0.6), -W.Normal);
			const double Above = FirstHit(Bar, From + W.Up * (0.5 * W.HeightMm + 0.6), -W.Normal);
			Bad += (Open > 0.0 && Open < 6.0) || Beside < 3.5 || Beside > 4.5 || Above < 3.5 || Above > 4.5 ? 1 : 0;
		}
		return Bad;
	};
	CHECK(Check(*Frame, M.SettingWindows) == 0);
	CHECK(Check(*Carriage, M.AccumulatorWindows) == 0);
	CHECK(Check(*Carriage, M.CounterWindows) == 0);
	CHECK(M.SettingWindows.size() == 20 && M.AccumulatorWindows.size() == 40 && M.CounterWindows.size() == 20);
	// Decimal markers ride on the bars' top edges, between digits 0 and 1 at Places = 0... and slide
	// one digit per pitch along local +X (toward the higher digits).
	const Pose S = M.SettingMarker.Poses[0];
	CHECK_NEAR(S.Origin.X, M.SettingWindows[0].Center.X + 0.5 * M.DigitPitchMm, 1e-9);
	CHECK_NEAR(S.Moved(Vec3(3.0 * M.DigitPitchMm, 0, 0)).Origin.X, M.SettingWindows[3].Center.X + 0.5 * M.DigitPitchMm, 1e-9);
	CHECK(M.CarriageMarkers.Poses.size() == 2 && M.CarriageMarkers.OnCarriage);
}

// One white carry lever per accumulator wheel, under it, between the wheel and the next higher one.
TEST(gears_carry_levers_under_wheels)
{
	const mgeo::MechanismParts& M = Machine();
	CHECK(M.CarryLevers.Poses.size() == M.AccumulatorWheels.Poses.size());
	const Mesh* Lever = Section(M.CarryLevers.Model, mgeo::MechMaterial::WhiteEnamel);
	const Mesh* Rim = Section(M.AccumulatorWheels.Model, mgeo::MechMaterial::NumeralBlack);
	CHECK(Lever && Rim);
	if (!Lever || !Rim)
	{
		return;
	}
	double RimR = 0.0;
	for (const Vec3& Q : Rim->Positions)
	{
		RimR = std::max(RimR, std::hypot(Q.X, Q.Y));
	}
	for (size_t I = 0; I < M.CarryLevers.Poses.size(); ++I)
	{
		const Pose& Lv = M.CarryLevers.Poses[I];
		const Pose& W = M.AccumulatorWheels.Poses[I];
		Box3 B;
		B.Add(*Lever, Lv);
		CHECK(B.Hi.Z < W.Origin.Z && B.Lo.Z > W.Origin.Z - RimR - 12.0);
		CHECK(Lv.Origin.X < W.Origin.X && Lv.Origin.X > W.Origin.X - M.DigitPitchMm);
		CHECK(std::fabs(Lv.Z.Dot(Vec3(1, 0, 0)) - 1.0) < 1e-9);
	}
}

TEST(gears_mechanism_parts_closed)
{
	const mgeo::MechanismParts& M = Machine();
	struct Named
	{
		const char* Name;
		const mgeo::MechPart* Part;
	};
	const Named Parts[] = {{"plinth", &M.Plinth}, {"frame", &M.Frame}, {"glass", &M.GlassCover}, {"printer", &M.Printer}, {"paper", &M.PaperStrip}, {"bell", &M.Bell},
		{"memory frame", &M.MemoryFrame}, {"drum", &M.Drums.Model}, {"square shaft", &M.SquareShafts.Model}, {"pinion", &M.Pinions.Model}, {"transfer", &M.TransferGears.Model},
		{"setting wheel", &M.SettingWheels.Model}, {"setting marker", &M.SettingMarker.Model}, {"carriage", &M.Carriage.Model}, {"accumulator wheel", &M.AccumulatorWheels.Model},
		{"carry lever", &M.CarryLevers.Model}, {"counter wheel", &M.CounterWheels.Model}, {"carriage marker", &M.CarriageMarkers.Model}, {"clearing lever", &M.ClearingLevers.Model},
		{"counting finger", &M.CountingFinger.Model}, {"main shaft", &M.MainShaft.Model}, {"reverser idler", &M.ReverserIdler.Model}, {"input shaft", &M.InputShaft.Model},
		{"drive gear", &M.DriveGear.Model}, {"flywheel", &M.Flywheel.Model}, {"crank", &M.Crank.Model}, {"crank handle", &M.CrankHandle.Model},
		{"reversing lever", &M.ReversingLever.Model}, {"shift knob", &M.ShiftKnob.Model}, {"shift pinion", &M.ShiftPinion.Model}, {"bell hammer", &M.BellHammer.Model},
		{"type wheel", &M.TypeWheels.Model}, {"paper roll", &M.PaperRoll.Model}, {"platen", &M.Platen.Model}, {"ribbon spool", &M.RibbonSpools.Model},
		{"memory wheel", &M.MemoryWheels.Model}, {"sign wheel", &M.MemorySignWheels.Model}};
	int Unique = 0;
	for (const Named& N : Parts)
	{
		CHECK(!N.Part->Sections.empty());
		for (const mgeo::MechSection& S : N.Part->Sections)
		{
			CheckClosed(std::string(N.Name) + "/" + mgeo::MechMaterialName(S.Material), S.Geometry);
			Unique += S.Geometry.TriangleCount();
		}
	}
	// Each wheel and gear type is one mesh: the unique geometry stays within a real-time budget.
	CHECK(Unique < 400000);
	CHECK(M.AccumulatorWheels.Model.TriangleCount() < 12000 && M.Drums.Model.TriangleCount() < 6000);
	CHECK(M.MemoryWheels.Poses.size() == 6 * 21 && M.MemorySignWheels.Poses.size() == 6 && M.MemoryLabels.size() == 6);
	CHECK(M.TypeWheels.Poses.size() == 18 && M.Drums.Poses.size() == 20 && M.TransferGears.Poses.size() == 20);
	std::printf("    unique mechanism triangles: %d\n", Unique);
}

// The machine stands behind the keyboard, the plinth under both, and the carriage clears every fixed
// part at every shift.
TEST(gears_mechanism_clearances)
{
	const mgeo::MechanismParts& M = Machine();
	const mgeo::BoardParts B = mgeo::BuildBoard(mkb::FullSizeAnsi(), mgeo::CaseParams());
	const double Ct = std::cos(B.RestTiltRad), St = std::sin(B.RestTiltRad);
	double CaseMaxY = -1e30, CaseMinZ = 1e30;
	for (const Mesh* K : {&B.CaseShell, &B.CaseBadge, &B.Feet})
	{
		for (const Vec3& P : K->Positions)
		{
			CaseMaxY = std::max(CaseMaxY, P.Y * Ct - P.Z * St);
			CaseMinZ = std::min(CaseMinZ, P.Y * St + P.Z * Ct + B.RestLiftMm);
		}
	}
	CHECK_NEAR(CaseMinZ, 0.0, 1e-9);
	// Everything but the plinth (and the overhanging carriage at every shift) is behind the keyboard.
	double MinY = 1e30;
	auto Scan = [&MinY](const mgeo::MechPart& P, const Pose& At) {
		for (const mgeo::MechSection& S : P.Sections)
		{
			for (const Vec3& Q : S.Geometry.Positions)
			{
				MinY = std::min(MinY, At.Point(Q).Y);
			}
		}
	};
	const Pose Id;
	for (const mgeo::MechPart* P : {&M.Frame, &M.Printer, &M.PaperStrip, &M.Bell, &M.MemoryFrame, &M.GlassCover})
	{
		Scan(*P, Id);
	}
	for (int K : {0, 20})
	{
		for (const mgeo::PartSet* S : {&M.Carriage, &M.AccumulatorWheels, &M.CounterWheels, &M.CarryLevers, &M.ClearingLevers, &M.CarriageMarkers})
		{
			for (size_t I = 0; I < S->Poses.size(); ++I)
			{
				Scan(S->Model, WorldPose(M, *S, I, K));
			}
		}
	}
	CHECK(MinY > CaseMaxY + 2.0);
	// The plinth runs under the keyboard and the machine, its top at z = 0.
	Box3 Pl;
	Pl.Add(*Section(M.Plinth, mgeo::MechMaterial::BlackEnamel), Id);
	CHECK_NEAR(Pl.Hi.Z, 0.0, 1e-9);
	CHECK(Pl.Lo.Y < -90.0 && Pl.Hi.Y > M.BoundsMax.Y && Pl.Lo.X < -230.0);
	// The memory columns are clear of the rest of the machine.
	Box3 Mem, Rest;
	Mem.Add(*Section(M.MemoryFrame, mgeo::MechMaterial::CastIron), Id);
	for (size_t I = 0; I < M.MemoryWheels.Poses.size(); ++I)
	{
		Mem.Add(*Section(M.MemoryWheels.Model, mgeo::MechMaterial::Brass), M.MemoryWheels.Poses[I]);
	}
	for (const mgeo::MechPart* P : {&M.Frame, &M.Printer, &M.PaperStrip, &M.Bell})
	{
		for (const mgeo::MechSection& S : P->Sections)
		{
			Rest.Add(S.Geometry, Id);
		}
	}
	for (const mgeo::PartSet* S : {&M.Flywheel, &M.Crank, &M.CrankHandle, &M.DriveGear, &M.PaperRoll})
	{
		for (size_t I = 0; I < S->Poses.size(); ++I)
		{
			for (const mgeo::MechSection& Sec : S->Model.Sections)
			{
				Rest.Add(Sec.Geometry, S->Poses[I]);
			}
		}
	}
	CHECK(!Mem.Overlaps(Rest));
	// The carriage's frame clears the fixed frame, transfer gears and printer at every shift.
	const Mesh* Car = Section(M.Carriage.Model, mgeo::MechMaterial::BlackEnamel);
	std::vector<Vec3> Fixed;
	for (const mgeo::MechSection& S : M.Frame.Sections)
	{
		Fixed.insert(Fixed.end(), S.Geometry.Positions.begin(), S.Geometry.Positions.end());
	}
	for (size_t J = 0; J < M.TransferGears.Poses.size(); ++J)
	{
		for (const Vec3& Q : M.TransferGears.Model.Sections[0].Geometry.Positions)
		{
			Fixed.push_back(M.TransferGears.Poses[J].Point(Q));
		}
	}
	for (const mgeo::MechSection& S : M.Printer.Sections)
	{
		Fixed.insert(Fixed.end(), S.Geometry.Positions.begin(), S.Geometry.Positions.end());
	}
	for (int K = 0; K <= M.MaxShift; K += 5)
	{
		const Mesh C = Placed(*Car, CarriageAt(M, K));
		Box3 Cb;
		Cb.Add(C, Id);
		int Clash = 0;
		for (const Vec3& Q : Fixed)
		{
			if (Q.X > Cb.Lo.X && Q.X < Cb.Hi.X && Q.Y > Cb.Lo.Y && Q.Y < Cb.Hi.Y && Q.Z > Cb.Lo.Z && Q.Z < Cb.Hi.Z && Inside(C, Q))
			{
				++Clash;
			}
		}
		CHECK(Clash == 0);
		if (Clash)
		{
			std::printf("    shift %d: %d fixed vertices inside the carriage\n", K, Clash);
		}
	}
}

TEST(gears_tape_frame)
{
	const mgeo::MechanismParts& M = Machine();
	CHECK(M.TapePath.size() > 10 && M.TapeLinePitchMm > 3.0 && M.TapeUsableWidthMm < M.TapeWidthMm);
	const mgeo::TapeFrame F0 = mgeo::TapeLineFrame(M, 0.0);
	CHECK(F0.Up.Z > 0.99 && F0.Normal.Y < -0.99 && F0.Right.X > 0.99);
	const Mesh* Paper = Section(M.PaperStrip, mgeo::MechMaterial::Paper);
	CHECK(Paper != nullptr);
	if (!Paper)
	{
		return;
	}
	double Length = 0.0;
	for (size_t K = 1; K < M.TapePath.size(); ++K)
	{
		Length += (M.TapePath[K] - M.TapePath[K - 1]).Length();
	}
	// At least 30 lines of history are visible, every one lying on the paper.
	CHECK(Length > 30.0 * M.TapeLinePitchMm);
	for (int Line = 0; Line * M.TapeLinePitchMm < Length - 1.0; Line += 3)
	{
		const mgeo::TapeFrame F = mgeo::TapeLineFrame(M, Line * M.TapeLinePitchMm);
		const double T = FirstHit(*Paper, F.Origin + F.Normal * 2.0 + F.Right * 0.7, -F.Normal);   // off the strip quads' diagonals
		CHECK_NEAR(T, 2.0, 0.05);
		CHECK(std::fabs(F.Up.Dot(F.Normal)) < 1e-6 && std::fabs(F.Right.Dot(F.Up)) < 1e-6);
		// Lines run up the tape, toward older lines.
		if (Line > 0)
		{
			const mgeo::TapeFrame P = mgeo::TapeLineFrame(M, (Line - 1) * M.TapeLinePitchMm);
			CHECK((F.Origin - P.Origin).Dot(P.Up) > 0.0);
		}
	}
}
