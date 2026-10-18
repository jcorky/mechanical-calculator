// Geometry checks: every part is a closed, outward-wound surface of the right topology, and the parts
// fit together (caps on stems, switches in the plate, the display clear of the keys, the board level
// on its feet).
#include "KeyboardParts.h"
#include "TestFramework.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <tuple>

namespace geomtest_detail
{
using mgeo::Mesh;
using mgeo::Vec3;

// Connected shells, with positions welded on a 1e-5 mm grid.
int Components(const Mesh& M)
{
	std::map<std::tuple<long long, long long, long long>, int> Ids;
	std::vector<int> Weld(M.Positions.size());
	for (size_t i = 0; i < M.Positions.size(); ++i)
	{
		const Vec3& P = M.Positions[i];
		const auto Key = std::make_tuple(std::llround(P.X * 1e5), std::llround(P.Y * 1e5), std::llround(P.Z * 1e5));
		const auto It = Ids.emplace(Key, static_cast<int>(Ids.size())).first;
		Weld[i] = It->second;
	}
	std::vector<int> Parent(Ids.size());
	std::iota(Parent.begin(), Parent.end(), 0);
	auto Find = [&Parent](int X) {
		while (Parent[static_cast<size_t>(X)] != X)
		{
			Parent[static_cast<size_t>(X)] = Parent[static_cast<size_t>(Parent[static_cast<size_t>(X)])];
			X = Parent[static_cast<size_t>(X)];
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
	int Count = 0;
	for (size_t i = 0; i < Parent.size(); ++i)
	{
		if (Used[i] && Find(static_cast<int>(i)) == static_cast<int>(i))
		{
			++Count;
		}
	}
	return Count;
}

// Closed, outward (positive volume), normals agreeing with the winding, and chi = 2 * shells - 2 * genus.
void CheckSolid(const char* Name, const Mesh& M, int Shells, int Genus)
{
	const bool Closed = M.IsClosed();
	const double Volume = M.SignedVolume();
	const int Chi = M.EulerCharacteristic();
	const int Found = Components(M);
	const bool Normals = M.NormalsAgreeWithWinding(0.3);
	CHECK(Closed);
	CHECK(Volume > 0.0);
	CHECK(Found == Shells);
	CHECK(Chi == 2 * Shells - 2 * Genus);
	CHECK(Normals);
	CHECK(M.Normals.size() == M.Positions.size() && M.UVs.size() == M.Positions.size());
	if (!Closed || Volume <= 0.0 || Found != Shells || Chi != 2 * Shells - 2 * Genus || !Normals)
	{
		std::printf("    %s: closed %d volume %.2f shells %d (want %d) chi %d (want %d) normals %d\n", Name, Closed, Volume, Found, Shells, Chi, 2 * Shells - 2 * Genus, Normals);
	}
}

// Moller-Trumbore; returns the ray parameter of the hit or -1.
double RayTri(const Vec3& O, const Vec3& D, const Vec3& A, const Vec3& B, const Vec3& C)
{
	const Vec3 E1 = B - A, E2 = C - A;
	const Vec3 Pv = D.Cross(E2);
	const double Det = E1.Dot(Pv);
	if (std::fabs(Det) < 1e-14)
	{
		return -1.0;
	}
	const double Inv = 1.0 / Det;
	const Vec3 Tv = O - A;
	const double U = Tv.Dot(Pv) * Inv;
	if (U < 0.0 || U > 1.0)
	{
		return -1.0;
	}
	const Vec3 Qv = Tv.Cross(E1);
	const double V = D.Dot(Qv) * Inv;
	if (V < 0.0 || U + V > 1.0)
	{
		return -1.0;
	}
	const double T = E2.Dot(Qv) * Inv;
	return T > 1e-9 ? T : -1.0;
}

struct Hit
{
	double T = -1.0;
	Vec3 Normal;
};

Hit FirstHit(const Mesh& M, const Vec3& O, const Vec3& D)
{
	Hit H;
	for (size_t t = 0; t + 2 < M.Triangles.size(); t += 3)
	{
		const Vec3& A = M.Positions[static_cast<size_t>(M.Triangles[t])];
		const Vec3& B = M.Positions[static_cast<size_t>(M.Triangles[t + 1])];
		const Vec3& C = M.Positions[static_cast<size_t>(M.Triangles[t + 2])];
		const double T = RayTri(O, D, A, B, C);
		if (T > 0.0 && (H.T < 0.0 || T < H.T))
		{
			H.T = T;
			H.Normal = (B - A).Cross(C - A).Normalized();
		}
	}
	return H;
}

// Point inside a closed mesh by ray parity (a skewed direction avoids edge-on hits).
bool Inside(const Mesh& M, const Vec3& P)
{
	const Vec3 D = Vec3(0.1234, 0.3011, 0.9452).Normalized();
	int Count = 0;
	for (size_t t = 0; t + 2 < M.Triangles.size(); t += 3)
	{
		if (RayTri(P, D, M.Positions[static_cast<size_t>(M.Triangles[t])], M.Positions[static_cast<size_t>(M.Triangles[t + 1])], M.Positions[static_cast<size_t>(M.Triangles[t + 2])]) > 0.0)
		{
			++Count;
		}
	}
	return (Count & 1) == 1;
}

struct Bounds
{
	Vec3 Lo = Vec3(1e30, 1e30, 1e30), Hi = Vec3(-1e30, -1e30, -1e30);
	void Add(const Vec3& P)
	{
		Lo = Vec3(std::min(Lo.X, P.X), std::min(Lo.Y, P.Y), std::min(Lo.Z, P.Z));
		Hi = Vec3(std::max(Hi.X, P.X), std::max(Hi.Y, P.Y), std::max(Hi.Z, P.Z));
	}
};

Bounds BoundsOf(const Mesh& M)
{
	Bounds B;
	for (const Vec3& P : M.Positions)
	{
		B.Add(P);
	}
	return B;
}

mgeo::KeycapParams CapFor(const mkb::KeyDef& K)
{
	mgeo::KeycapParams P;
	P.WidthU = K.W;
	P.DepthU = K.H;
	P.ProfileRow = K.ProfileRow;
	P.Spacebar = K.W >= 6.0;
	P.HomingBar = K.Id == "F" || K.Id == "J" || K.Id == "Num5";
	return P;
}

mgeo::KeycapParams Cap(double W, double H, int Row)
{
	mgeo::KeycapParams P;
	P.WidthU = W;
	P.DepthU = H;
	P.ProfileRow = Row;
	return P;
}
} // namespace geomtest_detail

using namespace geomtest_detail;

TEST(geometry_mesh_primitives)
{
	// A unit cube from six separately shaded faces welds into one closed genus-0 surface.
	Mesh M;
	const Vec3 C[8] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
	const int F[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
	for (const auto& Q : F)
	{
		const Vec3 N = (C[Q[1]] - C[Q[0]]).Cross(C[Q[3]] - C[Q[0]]).Normalized();
		int I[4];
		for (int k = 0; k < 4; ++k)
		{
			I[k] = M.AddVertex(C[Q[k]], N, mgeo::Vec2{0, 0});
		}
		M.AddQuad(I[0], I[1], I[2], I[3]);
	}
	CHECK(M.IsClosed());
	CHECK_NEAR(M.SignedVolume(), 1.0, 1e-12);
	CHECK(M.EulerCharacteristic() == 2);
	CHECK(M.NormalsAgreeWithWinding(0.99));
	Mesh Open = M;
	Open.Triangles.resize(Open.Triangles.size() - 3);
	CHECK(!Open.IsClosed());
	Mesh Flipped = M;
	std::swap(Flipped.Triangles[0], Flipped.Triangles[1]);
	CHECK(!Flipped.IsClosed());
	CHECK(!Flipped.NormalsAgreeWithWinding(0.0));
	Mesh Two = M;
	Two.Append(M, Vec3(3, 0, 0));
	CHECK(Two.IsClosed() && Two.EulerCharacteristic() == 4);
	CHECK_NEAR(Two.SignedVolume(), 2.0, 1e-12);
}

TEST(geometry_keycaps_closed)
{
	struct Case
	{
		double W, H;
		int Row;
		bool Space, Homing;
	};
	const Case Cases[] = {{1, 1, 1, false, false}, {1, 1, 2, false, false}, {1, 1, 3, false, true}, {1, 1, 4, false, false}, {1.25, 1, 4, false, false},
		{1.5, 1, 2, false, false}, {1.75, 1, 3, false, false}, {2, 1, 1, false, false}, {2.25, 1, 3, false, false}, {2.75, 1, 4, false, false},
		{6.25, 1, 4, true, false}, {1, 2, 2, false, false}, {1, 2, 4, false, false}};
	double Prev = 0.0;
	for (const Case& C : Cases)
	{
		mgeo::KeycapParams P = Cap(C.W, C.H, C.Row);
		P.Spacebar = C.Space;
		P.HomingBar = C.Homing;
		const Mesh M = mgeo::BuildKeycap(P);
		char Name[64];
		std::snprintf(Name, sizeof(Name), "keycap %.2fx%.0f R%d", C.W, C.H, C.Row);
		CheckSolid(Name, M, C.Homing ? 2 : 1, 0);
		if (C.W == 1 && C.H == 1)
		{
			Prev = M.SignedVolume();
		}
		else
		{
			CHECK(M.SignedVolume() > Prev);
		}
		if (C.W == 1 && C.H == 1)
		{
			CHECK(M.TriangleCount() >= 800 && M.TriangleCount() <= 3000);
		}
		CHECK(M.TriangleCount() <= 8000);
	}
}

TEST(geometry_keycap_volume_and_profile)
{
	// A hollow 1u ABS/PBT cap weighs about 1 g (1.05-1.4 g/cm3), so its shell is 0.4-1.5 cm3.
	for (int Row = 1; Row <= 4; ++Row)
	{
		const double V = mgeo::BuildKeycap(Cap(1, 1, Row)).SignedVolume();
		CHECK(V > 400.0 && V < 1500.0);
	}
	// A solid block the size of the cap's envelope is far heavier than the hollow shell.
	const Mesh Cap1 = mgeo::BuildKeycap(Cap(1, 1, 1));
	const Bounds B = BoundsOf(Cap1);
	CHECK(Cap1.SignedVolume() < 0.6 * (B.Hi.X - B.Lo.X) * (B.Hi.Y - B.Lo.Y) * (B.Hi.Z - B.Lo.Z));
	// The footprint is unit - gap and the cap is as tall as its row says.
	CHECK_NEAR(B.Hi.X - B.Lo.X, mkb::kUnitMm - 0.8, 0.25);
	CHECK_NEAR(mgeo::KeycapHeightMm(1), 9.4, 1e-9);
	CHECK(mgeo::KeycapHeightMm(1) > mgeo::KeycapHeightMm(2) && mgeo::KeycapHeightMm(2) > mgeo::KeycapHeightMm(3) && mgeo::KeycapHeightMm(3) > mgeo::KeycapHeightMm(4));
	// KeycapHeightMm is the top plane at its back edge; the rounded top corners stand a little lower.
	CHECK(B.Hi.Z - B.Lo.Z <= mgeo::KeycapHeightMm(1) + 1e-9 && B.Hi.Z - B.Lo.Z > mgeo::KeycapHeightMm(1) - 0.4);
	// R1 tilts hardest toward the typist; R4 leans the other way.
	CHECK(mgeo::KeycapTopFor(Cap(1, 1, 1)).TiltRad > mgeo::KeycapTopFor(Cap(1, 1, 2)).TiltRad);
	CHECK(mgeo::KeycapTopFor(Cap(1, 1, 2)).TiltRad > mgeo::KeycapTopFor(Cap(1, 1, 3)).TiltRad);
	CHECK(mgeo::KeycapTopFor(Cap(1, 1, 4)).TiltRad < 0.0);
}

TEST(geometry_keycap_top_matches_legend_frame)
{
	const double Sizes[][2] = {{1, 1}, {2.25, 1}, {1, 2}, {6.25, 1}};
	for (int Row = 1; Row <= 4; ++Row)
	{
		for (const auto& Sz : Sizes)
		{
			mgeo::KeycapParams P = Cap(Sz[0], Sz[1], Row);
			P.Spacebar = Sz[0] >= 6.0;
			const Mesh M = mgeo::BuildKeycap(P);
			const mgeo::KeycapTop T = mgeo::KeycapTopFor(P);
			// The dish centre, seen from above.
			const Hit H = FirstHit(M, Vec3(1e-4, T.CenterY + 1e-4, 50.0), Vec3(0, 0, -1));
			CHECK(H.T > 0.0);
			CHECK_NEAR(50.0 - H.T, T.CenterZ, 0.02);
			CHECK(H.Normal.Dot(Vec3(0, -std::sin(T.TiltRad), std::cos(T.TiltRad))) > 0.9995);
			// The front wall, seen from the typist.
			const Hit F = FirstHit(M, Vec3(1e-4, -50.0, T.FrontFaceZ), Vec3(0, 1, 0));
			CHECK(F.T > 0.0);
			CHECK_NEAR(F.T - 50.0, T.FrontFaceY, 0.15);
			CHECK(F.Normal.Dot(Vec3(0, -std::cos(T.FrontTiltRad), std::sin(T.FrontTiltRad))) > 0.995);
			CHECK(T.Width > 10.0 && T.Depth > 10.0);
		}
	}
	// UVs on the top map the dish to 0..1.
	const Mesh M = mgeo::BuildKeycap(Cap(1, 1, 3));
	// Only the top carries U >= 0, so a legend material can mask the walls and underside.
	double Lo = 1e9, Hi = -1e9;
	int Top = 0, Stray = 0;
	for (size_t i = 0; i < M.Positions.size(); ++i)
	{
		if (M.UVs[i].U >= 0.0)
		{
			++Top;
			Stray += M.Normals[i].Z > 0.95 && M.Positions[i].Z > 5.0 ? 0 : 1;
			Lo = std::min({Lo, M.UVs[i].U, M.UVs[i].V});
			Hi = std::max({Hi, M.UVs[i].U, M.UVs[i].V});
		}
		else
		{
			Stray += M.UVs[i].U < -1.0 ? 0 : 1;
		}
	}
	CHECK(Top > 100 && Stray == 0);
	CHECK_NEAR(Lo, 0.0, 1e-6);
	CHECK_NEAR(Hi, 1.0, 1e-6);
}

TEST(geometry_switch_parts_closed)
{
	const mgeo::SwitchType Types[] = {mgeo::SwitchType::Linear, mgeo::SwitchType::Tactile, mgeo::SwitchType::Clicky};
	for (mgeo::SwitchType Type : Types)
	{
		mgeo::SwitchParams P;
		P.Type = Type;
		const mgeo::SwitchParts S = mgeo::BuildSwitch(P);
		CheckSolid("housing bottom", S.HousingBottom, 1, 0);
		CheckSolid("housing top", S.HousingTop, 3, 1);
		CheckSolid("stem", S.Stem, Type == mgeo::SwitchType::Clicky ? 7 : 6, 0);
		CheckSolid("spring", S.Spring, 1, 0);
		CheckSolid("contacts", S.Contacts, 2, 0);
		CheckSolid("pins", S.Pins, 2, 0);
		const int Total = S.HousingBottom.TriangleCount() + S.HousingTop.TriangleCount() + S.Stem.TriangleCount() + S.Spring.TriangleCount() + S.Contacts.TriangleCount() + S.Pins.TriangleCount();
		CHECK(Total > 1500 && Total < 7000);
		CHECK(S.SpringRestTopZ - S.SpringBottomZ > P.TravelMm + 2.0);
	}
	// Datasheet envelope: 15.6 mm square above the plate, 14 mm through it, 5 mm down to the PCB,
	// pins 3.3 mm below that, and the cross tip 15.2 mm above the PCB.
	const mgeo::SwitchParams P;
	const mgeo::SwitchParts S = mgeo::BuildSwitch(P);
	const Bounds Top = BoundsOf(S.HousingTop), Bot = BoundsOf(S.HousingBottom), Pins = BoundsOf(S.Pins), Stem = BoundsOf(S.Stem);
	CHECK_NEAR(Top.Hi.X - Top.Lo.X, 15.6, 1e-6);
	double BlockLow = 0.0;
	for (const Vec3& P : S.HousingBottom.Positions)
	{
		BlockLow = std::fabs(P.X) > 3.0 ? std::min(BlockLow, P.Z) : BlockLow;
	}
	CHECK_NEAR(BlockLow, -5.0, 1e-9);
	CHECK_NEAR(Bot.Lo.Z, -8.0, 1e-9);
	CHECK_NEAR(Pins.Lo.Z, -8.3, 1e-9);
	CHECK_NEAR(Stem.Hi.Z + P.StemTopZ + 5.0, 15.2, 0.11);
	CHECK(Top.Hi.Z <= P.StemTopZ);
}

TEST(geometry_stabilizers)
{
	CHECK_NEAR(mgeo::StabilizerSpacingMm(2.0), 11.938, 1e-9);
	CHECK_NEAR(mgeo::StabilizerSpacingMm(2.25), 11.938, 1e-9);
	CHECK_NEAR(mgeo::StabilizerSpacingMm(2.75), 11.938, 1e-9);
	CHECK_NEAR(mgeo::StabilizerSpacingMm(6.25), 50.0, 1e-9);
	CHECK(mgeo::StabilizerSpacingMm(1.75) == 0.0);
	const double Widths[] = {2.0, 2.25, 2.75, 6.25};
	for (double W : Widths)
	{
		for (int V = 0; V < 2; ++V)
		{
			const mgeo::StabilizerParts S = mgeo::BuildStabilizer(W, V == 1);
			CheckSolid("stab housings", S.Housings, 4, 0);
			CheckSolid("stab stems", S.Stems, 4, 0);
			CheckSolid("stab wire", S.Wire, 1, 0);
			const Bounds B = BoundsOf(S.Stems);
			const double Span = V == 1 ? B.Hi.Y - B.Lo.Y : B.Hi.X - B.Lo.X;
			CHECK_NEAR(Span, 2.0 * mgeo::StabilizerSpacingMm(W) + 5.6, 1e-6);
		}
	}
}

TEST(geometry_board_parts_closed)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mgeo::BoardParts B = mgeo::BuildBoard(L, mgeo::CaseParams());
	int Stabs = 0;
	for (const mkb::KeyDef& K : L.Keys)
	{
		Stabs += K.NeedsStabilizer() ? 2 : 0;
	}
	CheckSolid("case", B.CaseShell, 1, 1);
	CheckSolid("badge", B.CaseBadge, 1, 0);
	CheckSolid("feet", B.Feet, 4, 0);
	CheckSolid("plate", B.Plate, 1, static_cast<int>(L.Keys.size()) + Stabs);
	CheckSolid("pcb", B.Pcb, 1 + 2 * static_cast<int>(L.Keys.size()), 0);
	const int Parts = Components(B.PcbParts);
	CheckSolid("pcb parts", B.PcbParts, Parts, 0);
	CHECK(Parts >= 3 * static_cast<int>(L.Keys.size()) + 2);
	CheckSolid("display bezel", B.DisplayBezel, 1, 0);
	CheckSolid("display glass", B.DisplayGlass, 1, 0);
	CHECK(B.CaseShell.TriangleCount() < 50000);
	CHECK(B.Plate.TriangleCount() < 30000);
	std::printf("    triangles: case %d badge %d feet %d plate %d pcb %d parts %d bezel %d glass %d\n", B.CaseShell.TriangleCount(), B.CaseBadge.TriangleCount(), B.Feet.TriangleCount(),
		B.Plate.TriangleCount(), B.Pcb.TriangleCount(), B.PcbParts.TriangleCount(), B.DisplayBezel.TriangleCount(), B.DisplayGlass.TriangleCount());
}

TEST(geometry_switches_land_in_cutouts)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mgeo::CaseParams CP;
	const mgeo::BoardParts B = mgeo::BuildBoard(L, CP);
	const mgeo::SwitchParts S = mgeo::BuildSwitch(mgeo::SwitchParams());
	// Everything of the switch inside the plate's thickness.
	Bounds Band;
	for (const Mesh* M : {&S.HousingBottom, &S.HousingTop})
	{
		for (const Vec3& P : M->Positions)
		{
			if (P.Z > -CP.PlateThickness + 1e-6 && P.Z < -1e-6)
			{
				Band.Add(P);
			}
		}
	}
	CHECK(Band.Hi.X <= 7.0 && Band.Lo.X >= -7.0 && Band.Hi.Y <= 7.0 && Band.Lo.Y >= -7.0);
	// The clip shoulders reach past the cut-out just under the plate, so the switch is held.
	CHECK(BoundsOf(S.HousingTop).Hi.X > 7.0);
	for (const mkb::KeyDef& K : L.Keys)
	{
		const Vec3 C = mgeo::KeyCenter(K, L);
		const double Probe[5][2] = {{0, 0}, {Band.Lo.X, Band.Lo.Y}, {Band.Hi.X, Band.Lo.Y}, {Band.Hi.X, Band.Hi.Y}, {Band.Lo.X, Band.Hi.Y}};
		for (const auto& Q : Probe)
		{
			CHECK(FirstHit(B.Plate, Vec3(C.X + Q[0] * 0.999, C.Y + Q[1] * 0.999, 1.0), Vec3(0, 0, -1)).T < 0.0);
		}
		// Plate material right outside the cut-out, on the sides without a stabilizer slot.
		CHECK(FirstHit(B.Plate, Vec3(C.X, C.Y + 7.5, 1.0), Vec3(0, 0, -1)).T > 0.0);
		CHECK(FirstHit(B.Plate, Vec3(C.X, C.Y - 7.5, 1.0), Vec3(0, 0, -1)).T > 0.0);
	}
}

TEST(geometry_stabilizers_fit_plate_and_caps)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mgeo::CaseParams CP;
	const mgeo::BoardParts B = mgeo::BuildBoard(L, CP);
	const double Mount = mgeo::SwitchParams().StemTopZ;
	int Checked = 0;
	for (const mkb::KeyDef& K : L.Keys)
	{
		if (!K.NeedsStabilizer())
		{
			continue;
		}
		const bool Vertical = K.H >= 2.0 && K.W < 2.0;
		const mgeo::StabilizerParts S = mgeo::BuildStabilizer(Vertical ? K.H : K.W, Vertical);
		const Vec3 C = mgeo::KeyCenter(K, L);
		// The inserts go through the plate only where it has slots.
		for (const Vec3& P : S.Housings.Positions)
		{
			if (P.Z > -CP.PlateThickness + 1e-6 && P.Z < -1e-6)
			{
				CHECK(FirstHit(B.Plate, Vec3(C.X + P.X, C.Y + P.Y, 1.0), Vec3(0, 0, -1)).T < 0.0);
			}
		}
		// The stabilizer crosses sit in the cap's sockets: no vertex inside the cap's material.
		const Mesh Cap = mgeo::BuildKeycap(CapFor(K));
		int Clashes = 0;
		for (const Vec3& P : S.Stems.Positions)
		{
			Clashes += Inside(Cap, P) ? 1 : 0;
		}
		CHECK(Clashes == 0);
		// At rest the cap clears the wire (which swings down with the stems); fully down, it clears
		// the housings.
		Mesh Rest = Cap, Pressed = Cap;
		Rest.Translate(Vec3(0, 0, Mount));
		Pressed.Translate(Vec3(0, 0, Mount - 4.0));
		int WireClashes = 0, HousingClashes = 0;
		for (const Vec3& P : S.Wire.Positions)
		{
			WireClashes += Inside(Rest, P) ? 1 : 0;
		}
		for (const Vec3& P : S.Housings.Positions)
		{
			HousingClashes += Inside(Pressed, P) ? 1 : 0;
		}
		CHECK(WireClashes == 0 && HousingClashes == 0);
		if (Clashes || WireClashes || HousingClashes)
		{
			std::printf("    %s: stems %d wire %d housings %d vertices inside the cap\n", K.Id.c_str(), Clashes, WireClashes, HousingClashes);
		}
		++Checked;
	}
	CHECK(Checked == 8);
}

// Own check: the cap clears the switch both at rest and at full travel, and the stem's cross enters
// the cap's + cut-out without touching it.
TEST(geometry_keycap_clears_switch)
{
	const mgeo::SwitchParams SP;
	const mgeo::SwitchParts S = mgeo::BuildSwitch(SP);
	for (int Row = 1; Row <= 4; ++Row)
	{
		const Mesh Shell = mgeo::BuildKeycap(Cap(1, 1, Row));
		for (double Travel : {0.0, SP.TravelMm})
		{
			Mesh Placed = Shell;
			Placed.Translate(Vec3(0, 0, SP.StemTopZ - Travel));
			int Clashes = 0;
			for (const Vec3& P : S.HousingTop.Positions)
			{
				Clashes += Inside(Placed, P) ? 1 : 0;
			}
			CHECK(Clashes == 0);
			if (Clashes)
			{
				std::printf("    R%d travel %.1f: %d housing vertices inside the cap\n", Row, Travel, Clashes);
			}
		}
		int StemClashes = 0;
		for (const Vec3& P : S.Stem.Positions)
		{
			StemClashes += Inside(Shell, P) ? 1 : 0;
		}
		CHECK(StemClashes == 0);
	}
}

TEST(geometry_display_clear_of_keys)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mgeo::BoardParts B = mgeo::BuildBoard(L, mgeo::CaseParams());
	const Bounds D = BoundsOf(B.DisplayBezel);
	const double U = mkb::kUnitMm;
	for (const mkb::KeyDef& K : L.Keys)
	{
		const Vec3 C = mgeo::KeyCenter(K, L);
		const double Hx = 0.5 * K.W * U - 0.4, Hy = 0.5 * K.H * U - 0.4;
		const bool Apart = C.X + Hx <= D.Lo.X || C.X - Hx >= D.Hi.X || C.Y + Hy <= D.Lo.Y || C.Y - Hy >= D.Hi.Y;
		CHECK(Apart);
	}
	// The module stays inside its layout window and its glass faces the typist with UVs 0..1,
	// (0, 0) at the top-left corner as the typist sees it.
	const mkb::Rect& W = L.Display;
	CHECK(D.Lo.X >= (W.X - 0.5 * L.Width) * U && D.Hi.X <= (W.X + W.W - 0.5 * L.Width) * U);
	CHECK(D.Hi.Y <= (0.5 * L.Height - W.Y) * U && D.Lo.Y >= (0.5 * L.Height - W.Y - W.H) * U);
	const Hit H = FirstHit(B.DisplayGlass, B.DisplayCenter + Vec3(1e-4, 1e-4, 5.0), Vec3(0, 0, -1));
	CHECK_NEAR(H.T, 5.0, 1e-3);
	CHECK(H.Normal.Dot(Vec3(0, -std::sin(B.DisplayTiltRad), std::cos(B.DisplayTiltRad))) > 0.9999);
	CHECK(B.DisplayTiltRad > 0.0 && B.DisplayWidthMm > 3.0 * B.DisplayHeightMm);
	double MinU = 1e9, MaxU = -1e9, MinV = 1e9, MaxV = -1e9;
	Vec3 Corner;
	for (size_t i = 0; i < B.DisplayGlass.Positions.size(); ++i)
	{
		const mgeo::Vec3& N = B.DisplayGlass.Normals[i];
		if (N.Z < 0.9 || std::fabs(B.DisplayGlass.Positions[i].Z - B.DisplayCenter.Z) > 2.0)
		{
			continue;
		}
		const mgeo::Vec2& UV = B.DisplayGlass.UVs[i];
		MinU = std::min(MinU, UV.U);
		MaxU = std::max(MaxU, UV.U);
		MinV = std::min(MinV, UV.V);
		MaxV = std::max(MaxV, UV.V);
		if (UV.U + UV.V < 0.02)
		{
			Corner = B.DisplayGlass.Positions[i];
		}
	}
	CHECK(MinU > -1e-6 && MinV > -1e-6 && MaxU < 1.0 + 1e-6 && MaxV < 1.0 + 1e-6);
	CHECK(MaxU > 0.99 && MaxV > 0.99 && MinU < 0.01 && MinV < 0.01);
	CHECK(Corner.X < B.DisplayCenter.X && Corner.Y > B.DisplayCenter.Y);
}

// Own check: set on the desk with RestTiltRad and RestLiftMm, the feet touch z = 0, nothing else goes
// below it, the plate slopes up away from the typist, and the plate and PCB sit inside the case pocket.
TEST(geometry_board_rests_on_feet)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mgeo::CaseParams CP;
	const mgeo::BoardParts B = mgeo::BuildBoard(L, CP);
	const double S = std::sin(B.RestTiltRad), C = std::cos(B.RestTiltRad);
	auto Rest = [S, C, &B](const Vec3& P) { return P.Y * S + P.Z * C + B.RestLiftMm; };
	double FeetLow = 1e9, CaseLow = 1e9;
	for (const Vec3& P : B.Feet.Positions)
	{
		FeetLow = std::min(FeetLow, Rest(P));
	}
	for (const Vec3& P : B.CaseShell.Positions)
	{
		CaseLow = std::min(CaseLow, Rest(P));
	}
	CHECK_NEAR(FeetLow, 0.0, 1e-9);
	CHECK(CaseLow > 1.0);
	CHECK_NEAR(B.RestTiltRad, CP.TypingAngleDeg * 3.14159265358979323846 / 180.0, 1e-12);
	CHECK(Rest(Vec3(0, 50, 0)) > Rest(Vec3(0, -50, 0)));
	// Plate and PCB corners (top and bottom faces) are clear of the case material.
	for (const Mesh* M : {&B.Plate, &B.Pcb})
	{
		const Bounds Bx = BoundsOf(*M);
		for (double X : {Bx.Lo.X + 0.01, Bx.Hi.X - 0.01})
		{
			for (double Y : {Bx.Lo.Y + 0.01, Bx.Hi.Y - 0.01})
			{
				for (double Z : {Bx.Lo.Z + 0.01, Bx.Hi.Z - 0.01})
				{
					CHECK(!Inside(B.CaseShell, Vec3(X, Y, Z)));
				}
			}
		}
	}
	// The case floor is below everything under the plate.
	const mgeo::SwitchParts Sw = mgeo::BuildSwitch(mgeo::SwitchParams());
	CHECK(BoundsOf(Sw.Pins).Lo.Z > -CP.DepthBelowPlate + CP.WallMm);
	CHECK(BoundsOf(B.PcbParts).Lo.Z > -CP.DepthBelowPlate + CP.WallMm);
}
