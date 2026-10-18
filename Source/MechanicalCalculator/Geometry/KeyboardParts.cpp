// Procedural keyboard parts. Dimensions and their sources:
//
// Keycaps, Cherry profile. KeyV2 (github.com/rsheldiii/KeyV2, src/key_profiles/cherry.scad, "based off
// GMK keycap set"): 18.16 mm base, 11.85 x 14.64 mm top (the top is deeper than it is wide and skewed
// toward the back), 0.65 mm cylindrical dish, 0.6 mm stem inset, 5.5 mm round stem socket, and row
// tilts that step by 2.5, 2.5 and 6.5 degrees from R1 to R4. Row heights, measured from the skirt to
// the back edge of the top, are the commonly quoted GMK Cherry values: R1 9.4, R2 8.0, R3 7.4,
// R4 7.0 mm. The tilts used here keep KeyV2's steps with R1 leaning 7 degrees toward the typist and R4
// leaning 4.5 degrees away, so R3 is the low point of the sculpt and R4's front edge stands up.
//
// Switch, Cherry MX. Cherry MX datasheet and the ai03 plate-design wiki: 15.6 x 15.6 mm housing,
// 14.0 x 14.0 mm plate cut-out, 1.5 mm plate, 5.0 mm from plate top to PCB top, 1.6 mm PCB,
// 4.0 mm travel with 2.0 mm pre-travel, 4.0 mm centre post, pins at (-3.81, 2.54) and (2.54, 5.08) mm
// leaving the housing by 3.3 mm. Deskthority wiki: 15.2 mm from the PCB to the top of the stem, a
// cross about 3.7 mm tall with 4.0 mm arms 1.1-1.3 mm thick (the cap cut-out is 4.0 x 1.17 mm).
// So the slider's top (the cap mount) is 15.2 - 5.0 - 3.6 = 6.6 mm above the plate.
//
// Stabilizers: Costar/Cherry cut-outs 3.3 x 14 mm, wire spacing 11.938 mm per side for 2-2.75u,
// 19.05 mm for 3u, 47.625 mm for 6u, 50 mm for 6.25u and 57.15 mm for 7u.
#include "KeyboardParts.h"
#include "MeshPrimitives.h"

#include <algorithm>
#include <cmath>

namespace mgeo
{
namespace keyboardparts_detail
{
using namespace prim;

constexpr double kDeg = kPi / 180.0;

constexpr double kMxPlateToPcb = 5.0;
constexpr double kMxCrossArm = 4.0;
constexpr double kMxCrossWidth = 1.17;
constexpr double kMxCrossHeight = 3.7;
constexpr double kMxStemClearance = 0.05;
constexpr double kMxPin1X = -3.81, kMxPin1Y = 2.54, kMxPin2X = 2.54, kMxPin2Y = 5.08;
constexpr double kMxPinBelowHousing = 3.3;
constexpr double kMxCutout = 14.0;

constexpr double kCapStemInset = 0.6;
constexpr double kCapSocketRadius = 2.75;
constexpr double kCapSocketDepth = 3.9;
constexpr double kCapFrontInsetScale = 1.1;
constexpr double kCapBackInsetScale = 0.45;
constexpr double kCapWallBulge = 0.1;
constexpr double kSpacebarCrown = 0.35;

constexpr double kStabSlotW = 3.3;
constexpr double kStabSlotH = 14.0;

// --- Keycaps ---------------------------------------------------------------------------------------

struct CapRowProfile
{
	double BackMm;
	double TiltDeg;
};

CapRowProfile CapRow(int Row)
{
	switch (std::max(1, std::min(4, Row)))
	{
	case 1:
		return {9.4, 7.0};
	case 2:
		return {8.0, 4.5};
	case 3:
		return {7.4, 2.0};
	default:
		return {7.0, -4.5};
	}
}

struct CapShape
{
	double Hw0 = 0, Hd0 = 0, R0 = 0;
	double HwT = 0, HdT = 0, RT = 0, YcT = 0, YF = 0, YB = 0;
	double ZBack = 0, Tilt = 0, TanT = 0, Skirt = 0;
	double Dish = 0, DishR = 0, Crown = 0;
	bool Spacebar = false;

	double PlaneZ(double Y) const { return ZBack - (YB - Y) * TanT; }
	double Off(double X, double Y) const
	{
		if (Spacebar)
		{
			const double T = std::max(-1.0, std::min(1.0, (Y - YcT) / HdT));
			return Crown * (1.0 - T * T);
		}
		if (Dish <= 0.0)
		{
			return 0.0;
		}
		const double Xc = std::max(-HwT, std::min(HwT, X));
		return -(std::sqrt(DishR * DishR - Xc * Xc) - (DishR - Dish));
	}
	double TopZ(double X, double Y) const { return PlaneZ(Y) + Off(X, Y); }
	Vec3 TopNormal(double X, double Y) const
	{
		const double H = 1e-4;
		const double Dx = (TopZ(X + H, Y) - TopZ(X - H, Y)) / (2 * H);
		const double Dy = (TopZ(X, Y + H) - TopZ(X, Y - H)) / (2 * H);
		return Vec3(-Dx, -Dy, 1.0).Normalized();
	}
	double DeepestOff() const { return Spacebar ? 0.0 : Dish; }
};

CapShape MakeCapShape(const KeycapParams& P)
{
	CapShape S;
	const CapRowProfile Row = CapRow(P.ProfileRow);
	S.Hw0 = 0.5 * (P.WidthU * P.UnitMm - P.GapMm);
	S.Hd0 = 0.5 * (P.DepthU * P.UnitMm - P.GapMm);
	S.R0 = P.CornerRadiusMm;
	S.HwT = S.Hw0 - P.TopInsetMm;
	S.YF = -S.Hd0 + P.TopInsetMm * kCapFrontInsetScale;
	S.YB = S.Hd0 - P.TopInsetMm * kCapBackInsetScale;
	S.HdT = 0.5 * (S.YB - S.YF);
	S.YcT = 0.5 * (S.YB + S.YF);
	S.RT = std::min(P.TopCornerRadiusMm, 0.9 * std::min(S.HwT, S.HdT));
	S.Skirt = -kCapStemInset;
	S.Tilt = Row.TiltDeg * kDeg;
	S.TanT = std::tan(S.Tilt);
	S.ZBack = Row.BackMm + S.Skirt;
	S.Spacebar = P.Spacebar;
	S.Crown = P.Spacebar ? kSpacebarCrown : 0.0;
	S.Dish = P.Spacebar ? 0.0 : std::max(0.0, P.DishDepthMm);
	if (S.Dish > 0.0)
	{
		S.DishR = (S.HwT * S.HwT + S.Dish * S.Dish) / (2.0 * S.Dish);
	}
	return S;
}

void AddHomingBar(Mesh& M, const CapShape& S)
{
	const double By = S.YF + 2.0;
	const Splits None;
	auto Base = [&S](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() - 0.3; };
	auto Shoulder = [&S](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() + 0.3; };
	auto Crest = [&S](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() + 0.45; };
	const Loop2 Outer = RoundedRect(0.0, By, 2.6, 0.5, 0.45, 5, None);
	const Loop2 Inner = RoundedRect(0.0, By, 2.4, 0.3, 0.28, 5, None);
	const std::vector<Vec3> R0 = Lift(Outer, Base), R1 = Lift(Outer, Shoulder), R2 = Lift(Inner, Crest);
	Loft(M, {R0, R1, R2});
	const Vec3 Up = Vec3(0.0, -S.TanT, 1.0).Normalized();
	FanCap(M, R2, Vec3(0.0, By, Crest(0.0, By)), Up);
	FanCap(M, R0, Vec3(0.0, By, Base(0.0, By)), -Up);
}

// --- Plate ----------------------------------------------------------------------------------------

// Stabilizer axis and spacing for a key, or Spacing 0 when it needs none.
struct StabPlacement
{
	double Spacing = 0.0;
	bool Vertical = false;
};

StabPlacement StabFor(double WidthU, double DepthU)
{
	StabPlacement S;
	if (WidthU >= 2.0 - 1e-9)
	{
		S.Spacing = StabilizerSpacingMm(WidthU);
	}
	else if (DepthU >= 2.0 - 1e-9)
	{
		S.Spacing = StabilizerSpacingMm(DepthU);
		S.Vertical = true;
	}
	return S;
}

// --- Vintage typewriter keys ------------------------------------------------------------------------

constexpr double kVintageRadius = 7.4;
constexpr double kVintageDiscBottom = 6.0, kVintageDiscTop = 8.4, kVintageDish = 0.35;

struct VintageShape
{
	bool Round = true;
	double Hw = kVintageRadius, Hd = kVintageRadius, R = kVintageRadius;
	// Dish depth falls off with the distance from the shape's core segment (a point for round keys).
	double TopZ(double X, double Y) const
	{
		const double Core = std::max(0.0, Hw - Hd), CoreY = std::max(0.0, Hd - Hw);
		const double Dx = std::max(0.0, std::fabs(X) - Core), Dy = std::max(0.0, std::fabs(Y) - CoreY);
		const double Q = std::min(1.0, std::sqrt(Dx * Dx + Dy * Dy) / std::min(Hw, Hd));
		return kVintageDiscTop - kVintageDish * (1.0 - Q * Q);
	}
	Vec3 TopNormal(double X, double Y) const
	{
		const double H = 1e-4;
		return Vec3(-(TopZ(X + H, Y) - TopZ(X - H, Y)) / (2 * H), -(TopZ(X, Y + H) - TopZ(X, Y - H)) / (2 * H), 1.0).Normalized();
	}
	Loop2 Outline(double Offset) const
	{
		if (Round)
		{
			return CircleLoop(0, 0, R + Offset, 48);
		}
		return RoundedRect(0, 0, Hw + Offset, Hd + Offset, R + Offset, 10, SameSplits(Uniform(static_cast<int>(Hw / 2.5)), Uniform(static_cast<int>(Hd / 2.5))));
	}
};

VintageShape MakeVintage(const KeycapParams& P)
{
	VintageShape V;
	V.Round = P.WidthU < 1.2 && P.DepthU < 1.2;
	V.Hw = P.WidthU >= 1.2 ? 0.5 * P.WidthU * P.UnitMm - 2.5 : kVintageRadius;
	V.Hd = P.DepthU >= 1.2 ? 0.5 * P.DepthU * P.UnitMm - 2.5 : kVintageRadius;
	V.R = V.Round ? kVintageRadius : std::min(V.Hw, V.Hd) - 1.0;
	return V;
}

// A stem post with the + cut-out from below, from z = 0 up into the disc.
void AddVintageStem(Mesh& M, double Sx, double Sy)
{
	const Loop2 Circle = CircleLoop(Sx, Sy, kCapSocketRadius, 20);
	const double Top = kVintageDiscBottom + 0.4;
	const std::vector<Vec3> Bot = AtZ(Circle, 0.0), Up = AtZ(Circle, Top);
	Loft(M, {Bot, Up});
	FanCap(M, Up, Vec3(Sx, Sy, Top), Vec3(0, 0, 1));
	const Loop2 Cross = FromPoints(PlusOutline(Sx, Sy, kMxCrossArm, kMxCrossWidth));
	const std::vector<Vec3> CrossBot = AtZ(Cross, 0.0), CrossTop = AtZ(Cross, kCapSocketDepth);
	Bridge(M, Bot, CrossBot, Vec3(Sx, Sy, 0.0), Vec3(0, 0, -1));
	const std::vector<char> Sharp = AllSharp(Cross.P.size());
	Loft(M, {CrossTop, CrossBot}, &Sharp);
	const std::vector<int> Tris = EarClip(Cross.P);
	for (size_t t = 0; t + 2 < Tris.size(); t += 3)
	{
		FlatTri(M, CrossTop[static_cast<size_t>(Tris[t])], CrossTop[static_cast<size_t>(Tris[t + 1])], CrossTop[static_cast<size_t>(Tris[t + 2])], Vec3(0, 0, -1));
	}
}

Mesh BuildVintageKeycap(const KeycapParams& P)
{
	Mesh M;
	const VintageShape V = MakeVintage(P);
	const Loop2 O = V.Outline(0.0);
	const size_t N = O.P.size();
	// Disc: flat bottom, straight side, dished legend face with 0..1 UVs.
	const std::vector<Vec3> Bottom = AtZ(O, kVintageDiscBottom);
	const std::vector<Vec3> Rim = Lift(O, [&V](double X, double Y) { return V.TopZ(X, Y); });
	Loft(M, {Bottom, Rim});
	FanCap(M, Bottom, Vec3(0, 0, kVintageDiscBottom), Vec3(0, 0, -1));
	const int TopFirst = static_cast<int>(M.Positions.size());
	auto TopVertex = [&M, &V](double X, double Y) {
		const Vec2 UV{(X + V.Hw) / (2.0 * V.Hw), (Y + V.Hd) / (2.0 * V.Hd)};
		return M.AddVertex(Vec3(X, Y, V.TopZ(X, Y)), V.TopNormal(X, Y), UV);
	};
	const int Rings = 6;
	std::vector<std::vector<int>> Grid(Rings);
	for (int k = 0; k < Rings; ++k)
	{
		const double Scale = 1.0 - static_cast<double>(k) / Rings;
		for (size_t i = 0; i < N; ++i)
		{
			Grid[static_cast<size_t>(k)].push_back(TopVertex(k == 0 ? O.P[i].X : Scale * O.P[i].X, k == 0 ? O.P[i].Y : Scale * O.P[i].Y));
		}
	}
	const int Center = TopVertex(0.0, 0.0);
	const int TopEnd = static_cast<int>(M.Positions.size());
	for (int k = 0; k < Rings; ++k)
	{
		for (size_t i = 0; i < N; ++i)
		{
			const size_t j = (i + 1) % N;
			const std::vector<int>& G = Grid[static_cast<size_t>(k)];
			if (k + 1 < Rings)
			{
				M.AddQuad(G[i], G[j], Grid[static_cast<size_t>(k + 1)][j], Grid[static_cast<size_t>(k + 1)][i]);
			}
			else
			{
				M.AddTriangle(G[i], G[j], Center);
			}
		}
	}
	// Stem posts: the switch at the centre and the stabilizer posts of long keys.
	AddVintageStem(M, 0.0, 0.0);
	const StabPlacement Stab = StabFor(P.WidthU, P.DepthU);
	if (Stab.Spacing > 0.0)
	{
		for (double S = -1.0; S <= 1.0; S += 2.0)
		{
			AddVintageStem(M, Stab.Vertical ? 0.0 : S * Stab.Spacing, Stab.Vertical ? S * Stab.Spacing : 0.0);
		}
	}
	for (int i = 0; i < static_cast<int>(M.UVs.size()); ++i)
	{
		if (i < TopFirst || i >= TopEnd)
		{
			Vec2& UV = M.UVs[static_cast<size_t>(i)];
			UV.U = -2.0 + 0.5 * (UV.U - std::floor(UV.U));
		}
	}
	return M;
}
} // namespace keyboardparts_detail

Mesh BuildKeycapRim(const KeycapParams& P)
{
	using namespace keyboardparts_detail;
	Mesh M;
	if (P.Style != KeycapStyle::Vintage)
	{
		return M;
	}
	const VintageShape V = MakeVintage(P);
	// A closed ring profile (offset from the disc's edge, z), wrapping the disc's rim.
	const double Prof[][2] = {{-0.25, 5.7}, {0.75, 5.7}, {0.95, 6.1}, {0.95, 8.6}, {0.55, 9.0}, {-0.35, 9.0}, {-0.35, 8.55}};
	const int Np = 7;
	for (int k = 0; k < Np; ++k)
	{
		const double* A = Prof[k];
		const double* B = Prof[(k + 1) % Np];
		Loft(M, {AtZ(V.Outline(A[0]), A[1]), AtZ(V.Outline(B[0]), B[1])});
	}
	return M;
}

double KeycapHeightMm(int ProfileRow)
{
	return keyboardparts_detail::CapRow(ProfileRow).BackMm;
}

KeycapTop KeycapTopFor(const KeycapParams& P)
{
	using namespace keyboardparts_detail;
	KeycapTop T;
	if (P.Style == KeycapStyle::Vintage)
	{
		const VintageShape V = MakeVintage(P);
		T.CenterY = 0.0;
		T.CenterZ = V.TopZ(0.0, 0.0);
		T.TiltRad = 0.0;
		T.Width = 2.0 * (V.Hw - 0.35);
		T.Depth = 2.0 * (V.Hd - 0.35);
		T.FrontFaceY = -(V.Hd + 0.95);
		T.FrontFaceZ = 7.35;
		T.FrontTiltRad = 0.0;
		return T;
	}
	const CapShape S = MakeCapShape(P);
	T.CenterY = S.YcT;
	T.CenterZ = S.TopZ(0.0, S.YcT);
	T.TiltRad = S.Tilt;
	T.Width = 2.0 * S.HwT;
	T.Depth = 2.0 * S.HdT / std::cos(S.Tilt);
	const double TopFrontZ = S.TopZ(0.0, S.YF);
	T.FrontFaceY = 0.5 * (-S.Hd0 + S.YF);
	T.FrontFaceZ = 0.5 * (S.Skirt + TopFrontZ);
	T.FrontTiltRad = std::atan2(S.YF + S.Hd0, TopFrontZ - S.Skirt);
	return T;
}

Mesh BuildKeycap(const KeycapParams& P)
{
	using namespace keyboardparts_detail;
	if (P.Style == KeycapStyle::Vintage)
	{
		return BuildVintageKeycap(P);
	}
	Mesh M;
	const CapShape S = MakeCapShape(P);
	const double W = P.WallMm;
	const int Seg = std::max(2, P.Segments);

	// Stem sockets: the switch at the centre, plus the stabilizer posts on long keys.
	const StabPlacement Stab = StabFor(P.WidthU, P.DepthU);
	std::vector<double> SocketAt = {0.0};
	if (Stab.Spacing > 0.0)
	{
		SocketAt = {-Stab.Spacing, 0.0, Stab.Spacing};
	}

	// The inner cavity: the walls offset by W, under a ceiling parallel to the top.
	auto CeilingZ = [&S, W](double, double Y) { return S.PlaneZ(Y) - S.DeepestOff() - W; };
	const double Tm = (CeilingZ(0.0, S.YcT) - S.Skirt) / (S.PlaneZ(S.YcT) - S.Skirt);
	const double HwI = Lerp(S.Hw0, S.HwT, Tm) - W;
	const double YFI = Lerp(-S.Hd0, S.YF, Tm) + W;
	const double YBI = Lerp(S.Hd0, S.YB, Tm) - W;
	const double HdI = 0.5 * (YBI - YFI), YcI = 0.5 * (YBI + YFI);
	const double RI = std::max(0.4, Lerp(S.R0, S.RT, Tm) - W);
	const double ExI = HwI - RI, EyI = HdI - RI;

	// Edge points shared by every loop: uniform ones for the dish, exact ones on the ceiling cuts
	// between sockets, and dense ones around the sockets so the ceiling bridges stay well shaped.
	const double XStep = S.Spacebar ? 1e9 : 1.3 * std::max(1.0, P.WidthU / 1.5);
	const double YStep = S.Spacebar ? 1.0 : 2.0 * std::max(1.0, P.DepthU / 1.5);
	const int XCount = std::max(S.Spacebar ? 4 : 3, std::min(60, static_cast<int>(2.0 * (S.HwT - S.RT) / XStep)));
	const int YCount = std::max(2, std::min(30, static_cast<int>(2.0 * (S.HdT - S.RT) / YStep)));
	std::vector<double> XBase = Uniform(XCount), YBase = Uniform(YCount);
	std::vector<double> XExact, YExact;
	const double Dense[] = {0.0, 1.5, 3.0, 4.5, 6.5, 9.0, 13.0, 18.0, 25.0};
	const bool Sparse = (Stab.Vertical ? 2.0 * EyI / (YCount + 1) : 2.0 * ExI / (XCount + 1)) > 3.0;
	for (size_t k = 0; k + 1 < SocketAt.size(); ++k)
	{
		const double Cut = 0.5 * (SocketAt[k] + SocketAt[k + 1]);
		if (Stab.Vertical)
		{
			YExact.push_back((Cut - YcI) / EyI);
		}
		else
		{
			XExact.push_back(Cut / ExI);
		}
	}
	if (SocketAt.size() > 1 && Sparse)
	{
		for (double C : SocketAt)
		{
			for (double D : Dense)
			{
				for (double Sgn = -1.0; Sgn <= 1.0; Sgn += 2.0)
				{
					if (Stab.Vertical)
					{
						YBase.push_back((C + Sgn * D - YcI) / EyI);
					}
					else
					{
						XBase.push_back((C + Sgn * D) / ExI);
					}
				}
			}
		}
	}
	const Splits Sp = SameSplits(MergeFractions(XBase, XExact, 0.4 / ExI), MergeFractions(YBase, YExact, 0.4 / EyI));

	const Loop2 Base = RoundedRect(0.0, 0.0, S.Hw0, S.Hd0, S.R0, Seg, Sp);
	const Loop2 Top = RoundedRect(0.0, S.YcT, S.HwT, S.HdT, S.RT, Seg, Sp);
	const Loop2 InBase = RoundedRect(0.0, 0.0, S.Hw0 - W, S.Hd0 - W, std::max(0.4, S.R0 - W), Seg, Sp);
	const Loop2 InTop = RoundedRect(0.0, YcI, HwI, HdI, RI, Seg, Sp);
	const size_t N = Base.P.size();

	// Outer walls, lofted from the footprint to the dished top outline with a slight outward bulge.
	const std::vector<Vec3> TopRing = Lift(Top, [&S](double X, double Y) { return S.TopZ(X, Y); });
	const std::vector<Vec3> BaseRing = AtZ(Base, S.Skirt);
	const int WallSteps = 4;
	std::vector<std::vector<Vec3>> Walls;
	for (int k = 0; k <= WallSteps; ++k)
	{
		if (k == 0 || k == WallSteps)
		{
			Walls.push_back(k == 0 ? BaseRing : TopRing);
			continue;
		}
		const double T = static_cast<double>(k) / WallSteps;
		const double Bulge = kCapWallBulge * std::sin(kPi * T);
		std::vector<Vec3> Ring;
		for (size_t i = 0; i < N; ++i)
		{
			Ring.push_back(Vec3(Lerp(Base.P[i].X, Top.P[i].X, T) + Base.N[i].X * Bulge, Lerp(Base.P[i].Y, Top.P[i].Y, T) + Base.N[i].Y * Bulge, Lerp(S.Skirt, TopRing[i].Z, T)));
		}
		Walls.push_back(Ring);
	}
	Loft(M, Walls);

	// The top: concentric copies of the outline shrunk toward its centre, on the dish surface. Its UVs
	// span 0..1 over the top; every other vertex of the cap gets U < 0 so a legend can be masked.
	const int TopFirst = static_cast<int>(M.Positions.size());
	const int Rings = S.Spacebar ? 4 : 7;
	std::vector<std::vector<int>> Grid(Rings);
	auto TopVertex = [&M, &S](double X, double Y) {
		const Vec2 UV{(X + S.HwT) / (2.0 * S.HwT), (Y - S.YF) / (2.0 * S.HdT)};
		return M.AddVertex(Vec3(X, Y, S.TopZ(X, Y)), S.TopNormal(X, Y), UV);
	};
	for (int k = 0; k < Rings; ++k)
	{
		const double Scale = 1.0 - static_cast<double>(k) / Rings;
		for (size_t i = 0; i < N; ++i)
		{
			const double X = k == 0 ? Top.P[i].X : Scale * Top.P[i].X;
			const double Y = k == 0 ? Top.P[i].Y : S.YcT + Scale * (Top.P[i].Y - S.YcT);
			Grid[static_cast<size_t>(k)].push_back(TopVertex(X, Y));
		}
	}
	const int Center = TopVertex(0.0, S.YcT);
	const int TopEnd = static_cast<int>(M.Positions.size());
	for (int k = 0; k < Rings; ++k)
	{
		for (size_t i = 0; i < N; ++i)
		{
			const size_t j = (i + 1) % N;
			const std::vector<int>& G = Grid[static_cast<size_t>(k)];
			if (k + 1 < Rings)
			{
				const std::vector<int>& H = Grid[static_cast<size_t>(k + 1)];
				M.AddQuad(G[i], G[j], H[j], H[i]);
			}
			else
			{
				M.AddTriangle(G[i], G[j], Center);
			}
		}
	}

	// Underside: the skirt's bottom rim and the inner walls.
	const std::vector<Vec3> InBaseRing = AtZ(InBase, S.Skirt);
	const std::vector<Vec3> InTopRing = Lift(InTop, CeilingZ);
	StripFlat(M, BaseRing, InBaseRing, Vec3(0, 0, -1));
	Loft(M, {InTopRing, InBaseRing});

	// Ceiling, split between the sockets, and the sockets with their + cut-outs.
	const Vec3 CeilingN = Vec3(0.0, S.TanT, -1.0).Normalized();
	const int CircleN = 20;
	for (size_t r = 0; r < SocketAt.size(); ++r)
	{
		const double Sx = Stab.Vertical ? 0.0 : SocketAt[r];
		const double Sy = Stab.Vertical ? SocketAt[r] : 0.0;
		const double Lo = r == 0 ? -1e30 : 0.5 * (SocketAt[r - 1] + SocketAt[r]);
		const double Hi = r + 1 == SocketAt.size() ? 1e30 : 0.5 * (SocketAt[r] + SocketAt[r + 1]);
		std::vector<size_t> Picked;
		for (size_t i = 0; i < N; ++i)
		{
			const double C = Stab.Vertical ? InTop.P[i].Y : InTop.P[i].X;
			if (C >= Lo - 1e-6 && C <= Hi + 1e-6)
			{
				Picked.push_back(i);
			}
		}
		// Where the picked run jumps across a cut, the cut is filled in with points about 1 mm apart,
		// interpolated the same way from both sides so the neighbouring regions share them exactly.
		std::vector<Vec3> Region;
		for (size_t k = 0; k < Picked.size(); ++k)
		{
			const size_t A = Picked[k], Bn = Picked[(k + 1) % Picked.size()];
			Region.push_back(InTopRing[A]);
			if ((A + 1) % N == Bn)
			{
				continue;
			}
			const bool Fwd = InTopRing[A].X + InTopRing[A].Y < InTopRing[Bn].X + InTopRing[Bn].Y;
			const Vec3& From = Fwd ? InTopRing[A] : InTopRing[Bn];
			const Vec3& To = Fwd ? InTopRing[Bn] : InTopRing[A];
			const int Steps = std::max(2, static_cast<int>(std::ceil((To - From).Length() / 1.0)));
			for (int s = 1; s < Steps; ++s)
			{
				const int Q = Fwd ? s : Steps - s;
				Region.push_back(From + (To - From) * (static_cast<double>(Q) / Steps));
			}
		}
		const Loop2 Circle = CircleLoop(Sx, Sy, kCapSocketRadius, CircleN);
		const std::vector<Vec3> SocketTop = Lift(Circle, CeilingZ);
		const std::vector<Vec3> SocketBot = AtZ(Circle, 0.0);
		Bridge(M, Region, SocketTop, Vec3(Sx, Sy, CeilingZ(Sx, Sy)), CeilingN);
		Loft(M, {SocketBot, SocketTop});
		const Loop2 Cross = FromPoints(PlusOutline(Sx, Sy, kMxCrossArm, kMxCrossWidth));
		const std::vector<Vec3> CrossBot = AtZ(Cross, 0.0), CrossTop = AtZ(Cross, kCapSocketDepth);
		Bridge(M, SocketBot, CrossBot, Vec3(Sx, Sy, 0.0), Vec3(0, 0, -1));
		const std::vector<char> Sharp = AllSharp(Cross.P.size());
		Loft(M, {CrossTop, CrossBot}, &Sharp);
		const std::vector<int> Tris = EarClip(Cross.P);
		for (size_t t = 0; t + 2 < Tris.size(); t += 3)
		{
			FlatTri(M, CrossTop[static_cast<size_t>(Tris[t])], CrossTop[static_cast<size_t>(Tris[t + 1])], CrossTop[static_cast<size_t>(Tris[t + 2])], Vec3(0, 0, -1));
		}
	}

	if (P.HomingBar)
	{
		AddHomingBar(M, S);
	}
	for (int i = 0; i < static_cast<int>(M.UVs.size()); ++i)
	{
		if (i < TopFirst || i >= TopEnd)
		{
			Vec2& UV = M.UVs[static_cast<size_t>(i)];
			UV.U = -2.0 + 0.5 * (UV.U - std::floor(UV.U));
		}
	}
	return M;
}

SwitchParts BuildSwitch(const SwitchParams& P)
{
	using namespace keyboardparts_detail;
	SwitchParts Out;
	const double Mount = P.StemTopZ;
	const double HousingTop = Mount - 0.3;
	const double PcbTop = -kMxPlateToPcb;
	const double Floor = PcbTop + 0.8;
	const int Seg = 3;
	const Vec3 Up(0, 0, 1), Down(0, 0, -1);
	const Splits Sq = SameSplits(Uniform(5), Uniform(5));

	// Bottom housing: a cup that passes through the plate cut-out, narrowed above the clip ledge,
	// with the spring boss inside and the 4 mm centre post below.
	{
		Mesh& M = Out.HousingBottom;
		const Loop2 Full = RoundedRect(0, 0, 6.95, 6.95, 0.5, Seg, Sq);
		const Loop2 Narrow = RoundedRect(0, 0, 6.3, 6.95, 0.5, Seg, Sq);
		const Loop2 Pocket = RoundedRect(0, 0, 5.4, 5.9, 0.6, Seg, Sq);
		const double Ledge = -2.5;
		Loft(M, {AtZ(Full, PcbTop), AtZ(Full, Ledge)});
		StripFlat(M, AtZ(Full, Ledge), AtZ(Narrow, Ledge), Up);
		Loft(M, {AtZ(Narrow, Ledge), AtZ(Narrow, 0.0)});
		StripFlat(M, AtZ(Narrow, 0.0), AtZ(Pocket, 0.0), Up);
		Loft(M, {AtZ(Pocket, 0.0), AtZ(Pocket, Floor)});
		const Loop2 Boss = CircleLoop(0, 0, 1.5, 16);
		Bridge(M, AtZ(Pocket, Floor), AtZ(Boss, Floor), Vec3(0, 0, Floor), Up);
		Loft(M, {AtZ(Boss, Floor), AtZ(Boss, Floor + 1.0)});
		FanCap(M, AtZ(Boss, Floor + 1.0), Vec3(0, 0, Floor + 1.0), Up);
		const Loop2 Post = CircleLoop(0, 0, 1.95, 20);
		const Loop2 Tip = CircleLoop(0, 0, 1.45, 20);
		const double PostEnd = PcbTop - 3.0;
		Bridge(M, AtZ(Full, PcbTop), AtZ(Post, PcbTop), Vec3(0, 0, PcbTop), Down);
		Loft(M, {AtZ(Post, PostEnd + 0.5), AtZ(Post, PcbTop)});
		Loft(M, {AtZ(Tip, PostEnd), AtZ(Post, PostEnd + 0.5)});
		FanCap(M, AtZ(Tip, PostEnd), Vec3(0, 0, PostEnd), Down);
	}

	// Top housing: a hollow shell on the plate, tapering to the stem window, plus the two plate clips
	// whose shoulders hook under the 14 mm cut-out.
	{
		Mesh& M = Out.HousingTop;
		const Splits Sq4 = SameSplits(Uniform(4), Uniform(4));
		const Loop2 A = RoundedRect(0, 0, 7.8, 7.8, 0.8, Seg, Sq4);
		const Loop2 C = RoundedRect(0, 0, 7.45, 7.45, 0.8, Seg, Sq4);
		const Loop2 D = RoundedRect(0, 0, 5.7, 5.7, 1.4, Seg, Sq4);
		const Loop2 IA = RoundedRect(0, 0, 7.0, 7.0, 0.4, Seg, Sq4);
		const Loop2 ID = RoundedRect(0, 0, 5.16, 5.16, 0.6, Seg, Sq4);
		const Loop2 Win = RoundedRect(0, 0, 3.7, 2.9, 0.4, Seg, SameSplits(Uniform(3), Uniform(2)));
		const double Ceiling = HousingTop - 0.8;
		Loft(M, {AtZ(A, 0.0), AtZ(A, 1.3)});
		Loft(M, {AtZ(A, 1.3), AtZ(C, 1.9)});
		Loft(M, {AtZ(C, 1.9), AtZ(D, HousingTop)});
		Bridge(M, AtZ(D, HousingTop), AtZ(Win, HousingTop), Vec3(0, 0, HousingTop), Up);
		Loft(M, {AtZ(Win, HousingTop), AtZ(Win, Ceiling)});
		Bridge(M, AtZ(ID, Ceiling), AtZ(Win, Ceiling), Vec3(0, 0, Ceiling), Down);
		Loft(M, {AtZ(ID, Ceiling), AtZ(IA, 2.2)});
		Loft(M, {AtZ(IA, 2.2), AtZ(IA, 0.0)});
		StripFlat(M, AtZ(A, 0.0), AtZ(IA, 0.0), Down);
		const std::vector<Pt2> Clip = {{6.4, -2.4}, {6.95, -2.4}, {7.3, -1.75}, {7.3, -1.55}, {6.95, -1.55}, {6.95, 0.0}, {7.3, 0.0}, {7.3, 1.2}, {6.4, 1.2}};
		Frame3 F;
		F.V = Vec3(0, 0, 1);
		F.W = Vec3(0, -1, 0);
		std::vector<Pt2> Mirror;
		for (const Pt2& Q : Clip)
		{
			Mirror.push_back(Pt2{-Q.X, Q.Y});
		}
		Prism(M, Clip, -2.2, 2.2, F);
		Prism(M, Mirror, -2.2, 2.2, F);
	}

	// Stem, in its own frame (z = 0 is the cap mount): the + cross, the slider (hollow from below for
	// the spring), the guide wings, the legs that ride the contact leaf and the spring seat.
	{
		Mesh& M = Out.Stem;
		const double Arm = kMxCrossArm - 2.0 * kMxStemClearance, Thick = kMxCrossWidth - 2.0 * kMxStemClearance;
		Prism(M, PlusOutline(0, 0, Arm, Thick), -0.2, kMxCrossHeight, Frame3());
		const Loop2 Slider = RoundedRect(0, 0, 3.6, 2.8, 0.5, Seg, SameSplits(Uniform(4), Uniform(3)));
		const Loop2 Bore = CircleLoop(0, 0, 2.5, 20);
		const double SliderBot = -4.6, BoreTop = -1.6;
		Loft(M, {AtZ(Slider, SliderBot), AtZ(Slider, 0.0)});
		FanCap(M, AtZ(Slider, 0.0), Vec3(0, 0, 0), Up);
		Bridge(M, AtZ(Slider, SliderBot), AtZ(Bore, SliderBot), Vec3(0, 0, SliderBot), Down);
		Loft(M, {AtZ(Bore, BoreTop), AtZ(Bore, SliderBot)});
		FanCap(M, AtZ(Bore, BoreTop), Vec3(0, 0, BoreTop), Down);
		SolidFromLoop(M, RoundedRect(0, 0, 4.7, 2.0, 0.3, 2, Splits()), -2.4, -1.3);
		Cylinder(M, 0, 0, 1.2, 14, -4.2, BoreTop + 0.1);
		std::vector<Pt2> Leg = {{2.3, -6.7}, {3.3, -6.7}, {3.3, -4.4}, {2.3, -4.4}};
		if (P.Type != SwitchType::Linear)
		{
			Leg = {{2.3, -6.7}, {3.3, -6.7}, {3.3, -6.3}, {3.85, -5.75}, {3.3, -5.2}, {3.3, -4.4}, {2.3, -4.4}};
		}
		Frame3 F;
		F.U = Vec3(0, 1, 0);
		F.V = Vec3(0, 0, 1);
		F.W = Vec3(1, 0, 0);
		Prism(M, Leg, 1.25, 2.15, F);
		Prism(M, Leg, -2.15, -1.25, F);
		if (P.Type == SwitchType::Clicky)
		{
			SolidFromLoop(M, RoundedRect(0, -0.6, 4.1, 2.6, 0.6, Seg, Splits()), -5.5, -4.5);
		}
	}

	// Spring: a helix of SpringTurns turns with closed (touching) end coils, seated on the housing
	// floor and reaching up into the slider's bore.
	{
		const double Rc = 2.25, Rw = 0.15;
		Out.SpringBottomZ = Floor;
		Out.SpringRestTopZ = Mount - 1.6;
		const int Turns = std::max(4, P.SpringTurns);
		const int PerTurn = 12;
		const double Z0 = Out.SpringBottomZ + Rw, Z1 = Out.SpringRestTopZ - Rw;
		const double EndRise = 2.05 * Rw;
		const double MidRise = std::max(0.0, (Z1 - Z0) - 2.0 * EndRise);
		std::vector<Vec3> Path;
		for (int s = 0; s <= Turns * PerTurn; ++s)
		{
			const double U = static_cast<double>(s) / PerTurn;
			double Z = Z0;
			if (U <= 1.0)
			{
				Z = Z0 + EndRise * U;
			}
			else if (U >= Turns - 1.0)
			{
				Z = Z1 - EndRise * (Turns - U);
			}
			else
			{
				Z = Z0 + EndRise + MidRise * (U - 1.0) / (Turns - 2.0);
			}
			Path.push_back(Vec3(Rc * std::cos(2.0 * kPi * U), Rc * std::sin(2.0 * kPi * U), Z));
		}
		SweepTube(Out.Spring, Path, Rw, 6);
	}

	// Contact leaves: a straight fixed leaf and the bent moving leaf the stem legs push away.
	{
		Frame3 F;
		F.U = Vec3(0, 1, 0);
		F.V = Vec3(0, 0, 1);
		F.W = Vec3(1, 0, 0);
		Prism(Out.Contacts, RectPoly(5.0, Floor, 5.2, -0.9), -3.4, -1.0, F);
		const std::vector<Pt2> Bent = ThickenPolyline({{5.1, Floor}, {5.1, -2.9}, {4.35, -1.5}, {4.55, -0.9}}, 0.1);
		Prism(Out.Contacts, Bent, 0.8, 3.2, F);
	}

	// Pins: flat-ended copper pins 3.3 mm below the housing.
	{
		const double Pins[2][2] = {{kMxPin1X, kMxPin1Y}, {kMxPin2X, kMxPin2Y}};
		for (const auto& Q : Pins)
		{
			const Loop2 Sec = FromPoints(RectPoly(Q[0] - 0.35, Q[1] - 0.35, Q[0] + 0.35, Q[1] + 0.35));
			const Loop2 TipSec = FromPoints(RectPoly(Q[0] - 0.15, Q[1] - 0.15, Q[0] + 0.15, Q[1] + 0.15));
			const std::vector<char> Sharp = AllSharp(4);
			const double End = PcbTop - kMxPinBelowHousing;
			Loft(Out.Pins, {AtZ(TipSec, End), AtZ(Sec, End + 0.4)}, &Sharp);
			Loft(Out.Pins, {AtZ(Sec, End + 0.4), AtZ(Sec, Floor - 0.2)}, &Sharp);
			FanCap(Out.Pins, AtZ(TipSec, End), Vec3(Q[0], Q[1], End), Down);
			FanCap(Out.Pins, AtZ(Sec, Floor - 0.2), Vec3(Q[0], Q[1], Floor - 0.2), Up);
		}
	}
	return Out;
}

double StabilizerSpacingMm(double WidthU)
{
	if (WidthU < 2.0 - 1e-9)
	{
		return 0.0;
	}
	if (WidthU < 3.0 - 1e-9)
	{
		return 11.938;
	}
	if (WidthU < 6.0 - 1e-9)
	{
		return 19.05;
	}
	if (WidthU < 6.25 - 1e-9)
	{
		return 47.625;
	}
	if (WidthU < 7.0 - 1e-9)
	{
		return 50.0;
	}
	return 57.15;
}

StabilizerParts BuildStabilizer(double WidthU, bool Vertical)
{
	using namespace keyboardparts_detail;
	StabilizerParts Out;
	const double S = StabilizerSpacingMm(WidthU);
	if (S <= 0.0)
	{
		return Out;
	}
	const Vec3 Up(0, 0, 1), Down(0, 0, -1);
	const double BodyTop = 4.8, WindowFloor = 0.3, WireZ = 5.6, BarY = 8.6, ArmX = S + 4.4, WireR = 0.6;
	for (double Side = -1.0; Side <= 1.0; Side += 2.0)
	{
		const double Cx = Side * S;
		const Loop2 Body = RoundedRect(Cx, 0, 3.6, 5.5, 0.8, 3, SameSplits(Uniform(3), Uniform(4)));
		const Loop2 Win = RoundedRect(Cx, 0, 3.0, 3.0, 0.5, 3, SameSplits(Uniform(3), Uniform(3)));
		Loft(Out.Housings, {AtZ(Body, 0.0), AtZ(Body, BodyTop)});
		FanCap(Out.Housings, AtZ(Body, 0.0), Vec3(Cx, 0, 0.0), Down);
		Bridge(Out.Housings, AtZ(Body, BodyTop), AtZ(Win, BodyTop), Vec3(Cx, 0, BodyTop), Up);
		Loft(Out.Housings, {AtZ(Win, BodyTop), AtZ(Win, WindowFloor)});
		FanCap(Out.Housings, AtZ(Win, WindowFloor), Vec3(Cx, 0, WindowFloor), Up);
		SolidFromLoop(Out.Housings, RoundedRect(Cx, 0, 0.5 * kStabSlotW - 0.05, 0.5 * kStabSlotH - 0.1, 0.3, 2, Splits()), -2.2, 0.05);

		SolidFromLoop(Out.Stems, RoundedRect(Cx, 0, 2.8, 2.8, 0.4, 3, Splits()), -2.2, 0.0);
		Prism(Out.Stems, PlusOutline(Cx, 0, kMxCrossArm - 2.0 * kMxStemClearance, kMxCrossWidth - 2.0 * kMxStemClearance), -0.2, kMxCrossHeight, Frame3());
	}
	const double EndZ = WireZ;
	const std::vector<Vec3> Corners = {Vec3(-S + 1.2, 0, EndZ), Vec3(-ArmX, 0, EndZ), Vec3(-ArmX, BarY, WireR + 0.05), Vec3(ArmX, BarY, WireR + 0.05), Vec3(ArmX, 0, EndZ), Vec3(S - 1.2, 0, EndZ)};
	SweepTube(Out.Wire, FilletPath(Corners, 1.5, 6), WireR, 8);
	if (Vertical)
	{
		RotateAboutZ(Out.Housings, 0.5 * kPi);
		RotateAboutZ(Out.Stems, 0.5 * kPi);
		RotateAboutZ(Out.Wire, 0.5 * kPi);
	}
	return Out;
}

Vec3 KeyCenter(const mkb::KeyDef& Key, const mkb::Layout& Layout)
{
	return Vec3((Key.CenterX() - 0.5 * Layout.Width) * mkb::kUnitMm, (0.5 * Layout.Height - Key.CenterY()) * mkb::kUnitMm, 0.0);
}

BoardParts BuildBoard(const mkb::Layout& Layout, const CaseParams& P)
{
	using namespace keyboardparts_detail;
	BoardParts B;
	const double U = mkb::kUnitMm;
	const double Hw = 0.5 * Layout.Width * U, Hd = 0.5 * Layout.Height * U;
	const double PcbTop = -(P.PlateThickness + P.PcbGapBelowPlate);
	const double PcbBot = PcbTop - P.PcbThickness;
	const double Angle = P.TypingAngleDeg * kDeg;
	const double TanA = std::tan(Angle);
	const double M = P.MarginMm;
	const double YFront = -Hd - M, YBack = Hd + M;
	auto BottomZ = [&P, YFront, TanA](double, double Y) { return -P.DepthBelowPlate - (Y - YFront) * TanA; };
	const double EIn = 3.5, EPlate = 2.8, EPcb = 2.4;
	const double FloorZ = -P.DepthBelowPlate + P.WallMm;
	const Vec3 Up(0, 0, 1), Down(0, 0, -1);

	// The USB-C receptacle sits in the gap between F8 and F9, at the back edge of the PCB.
	double UsbX = 0.0;
	{
		const mkb::KeyDef* F8 = Layout.Find("F8");
		const mkb::KeyDef* F9 = Layout.Find("F9");
		if (F8 && F9)
		{
			UsbX = 0.5 * (KeyCenter(*F8, Layout).X + KeyCenter(*F9, Layout).X);
		}
	}
	const double UsbH = 3.26, UsbW = 8.94, UsbD = 7.35;
	const double UsbZ = PcbTop + 0.5 * UsbH;

	// Case: one closed surface from the wedge bottom up the outer wall, over the chamfer and the rim,
	// and down the pocket to its floor, with the USB-C port cut through the back wall.
	{
		Mesh& C = B.CaseShell;
		const int Seg = 10;
		const double R = P.CornerRadiusMm, Ch = P.ChamferMm, Lip = P.LipAbovePlate;
		const double ExOuter = Hw + M - R;
		std::vector<double> Back;
		const double Offsets[] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 12, 14, 17, 20, 24, 29, 35, 42, 50, 60, 72, 86, 102, 120, 140, 165, 195};
		for (double D : Offsets)
		{
			Back.push_back((UsbX + D) / ExOuter);
			Back.push_back((UsbX - D) / ExOuter);
		}
		Splits Sp;
		Sp.Top = MergeFractions(Back, {}, 0.2 / ExOuter);
		const int BackFrom = Seg, BackTo = Seg + static_cast<int>(Sp.Top.size()) + 1;
		auto Ring = [&](double E, double Rad, auto Z) { return Lift(RoundedRect(0, 0, Hw + E, Hd + E, Rad, Seg, Sp), Z); };
		auto Flat = [](double Z) { return [Z](double, double) { return Z; }; };
		const std::vector<Vec3> Ob = Ring(M - 1.0, R - 0.6, BottomZ);
		const std::vector<Vec3> O0 = Ring(M, R, [&BottomZ](double X, double Y) { return BottomZ(X, Y) + 1.0; });
		const std::vector<Vec3> O1 = Ring(M, R, Flat(Lip - Ch));
		const std::vector<Vec3> O2 = Ring(M - Ch, R - 0.6 * Ch, Flat(Lip));
		const std::vector<Vec3> I0 = Ring(EIn + 0.8, 2.6, Flat(Lip));
		const std::vector<Vec3> I1 = Ring(EIn, 2.0, Flat(Lip - 0.8));
		const std::vector<Vec3> I2 = Ring(EIn, 2.0, Flat(FloorZ));
		Loft(C, {Ob, O0});
		Loft(C, {O0, O1}, nullptr, BackFrom, BackTo);
		Loft(C, {O1, O2});
		Loft(C, {O2, I0});
		Loft(C, {I0, I1});
		Loft(C, {I1, I2}, nullptr, BackFrom, BackTo);
		FanCap(C, Ob, Vec3(0, 0, BottomZ(0, 0)), Vec3(0, -TanA, -1));
		FanCap(C, I2, Vec3(0, 0, FloorZ), Up);
		const Loop2 Port = RoundedRect(UsbX, UsbZ, 6.0, 3.0, 2.6, 6, SameSplits(Uniform(3), {}));
		auto PortRing = [&Port](double Y) {
			std::vector<Vec3> Out;
			for (const Pt2& Q : Port.P)
			{
				Out.push_back(Vec3(Q.X, Y, Q.Y));
			}
			return Out;
		};
		auto Face = [=](const std::vector<Vec3>& Lo, const std::vector<Vec3>& Hi) {
			std::vector<Vec3> Poly;
			for (int i = BackFrom; i <= BackTo; ++i)
			{
				Poly.push_back(Lo[static_cast<size_t>(i)]);
			}
			for (int i = BackTo; i >= BackFrom; --i)
			{
				Poly.push_back(Hi[static_cast<size_t>(i)]);
			}
			return Poly;
		};
		const double YIn = Hd + EIn;
		Bridge(C, Face(O0, O1), PortRing(YBack), Vec3(UsbX, YBack, UsbZ), Vec3(0, 1, 0));
		Bridge(C, Face(I2, I1), PortRing(YIn), Vec3(UsbX, YIn, UsbZ), Vec3(0, -1, 0));
		Loft(C, {PortRing(YIn), PortRing(YBack)});
	}

	// Badge: a small chamfered plate set into the back wall.
	{
		Mesh Badge;
		const double Bx = Hw - 40.0, Bz = 0.5 * (BottomZ(0, YBack) + 1.0 + P.LipAbovePlate - P.ChamferMm);
		const Loop2 Plate = RoundedRect(0, 0, 18.0, 4.5, 1.2, 5, SameSplits(Uniform(4), {}));
		const Loop2 Face = RoundedRect(0, 0, 17.6, 4.1, 0.9, 5, SameSplits(Uniform(4), {}));
		Loft(Badge, {AtZ(Plate, -0.8), AtZ(Plate, 0.4), AtZ(Face, 0.8)});
		FanCap(Badge, AtZ(Face, 0.8), Vec3(0, 0, 0.8), Up);
		FanCap(Badge, AtZ(Plate, -0.8), Vec3(0, 0, -0.8), Down);
		Frame3 F;
		F.O = Vec3(Bx, YBack, Bz);
		F.U = Vec3(1, 0, 0);
		F.V = Vec3(0, 0, -1);
		F.W = Vec3(0, 1, 0);
		ApplyFrame(Badge, F);
		B.CaseBadge = Badge;
	}

	// Feet: four rubber pucks square to the wedge bottom.
	{
		const Vec3 Out = Vec3(0, -std::sin(Angle), -std::cos(Angle));
		Frame3 F;
		F.U = Vec3(1, 0, 0);
		F.W = Out;
		F.V = F.W.Cross(F.U);
		const double Fx = Hw + M - 28.0;
		const double Fy[2] = {YFront + 18.0, YBack - 18.0};
		for (double Sx = -1.0; Sx <= 1.0; Sx += 2.0)
		{
			for (double Y : Fy)
			{
				Mesh Foot;
				const Loop2 Rim = CircleLoop(0, 0, 7.0, 28), Sole = CircleLoop(0, 0, 6.3, 28);
				Loft(Foot, {AtZ(Rim, -1.0), AtZ(Rim, 1.9), AtZ(Sole, 2.5)});
				FanCap(Foot, AtZ(Sole, 2.5), Vec3(0, 0, 2.5), Up);
				FanCap(Foot, AtZ(Rim, -1.0), Vec3(0, 0, -1.0), Down);
				F.O = Vec3(Sx * Fx, Y, BottomZ(0, Y));
				ApplyFrame(Foot, F);
				B.Feet.Append(Foot);
			}
		}
		B.RestTiltRad = Angle;
		double Lowest = 1e30;
		for (const Vec3& Q : B.Feet.Positions)
		{
			Lowest = std::min(Lowest, Q.Y * std::sin(Angle) + Q.Z * std::cos(Angle));
		}
		B.RestLiftMm = -Lowest;
	}

	// Plate: 14 x 14 mm cut-outs at every switch and 3.3 x 14 mm slots for the stabilizers.
	{
		std::vector<HoleRect> Holes;
		const double H = 0.5 * kMxCutout;
		for (const mkb::KeyDef& K : Layout.Keys)
		{
			const Vec3 C = KeyCenter(K, Layout);
			Holes.push_back({C.X - H, C.X + H, C.Y - H, C.Y + H});
			const StabPlacement Stab = StabFor(K.W, K.H);
			if (Stab.Spacing <= 0.0)
			{
				continue;
			}
			for (double Side = -1.0; Side <= 1.0; Side += 2.0)
			{
				if (Stab.Vertical)
				{
					const double Sy = C.Y + Side * Stab.Spacing;
					Holes.push_back({C.X - 0.5 * kStabSlotH, C.X + 0.5 * kStabSlotH, Sy - 0.5 * kStabSlotW, Sy + 0.5 * kStabSlotW});
				}
				else
				{
					const double Sx = C.X + Side * Stab.Spacing;
					Holes.push_back({Sx - 0.5 * kStabSlotW, Sx + 0.5 * kStabSlotW, C.Y - 0.5 * kStabSlotH, C.Y + 0.5 * kStabSlotH});
				}
			}
		}
		PerforatedSlab(B.Plate, -Hw - EPlate, Hw + EPlate, -Hd - EPlate, Hd + EPlate, Holes, -P.PlateThickness, 0.0);
	}

	// PCB with its pin pads, and the parts on it.
	{
		SolidFromLoop(B.Pcb, RoundedRect(0, 0, Hw + EPcb, Hd + EPcb, 1.5, 6, Splits()), PcbBot, PcbTop);
		for (const mkb::KeyDef& K : Layout.Keys)
		{
			const Vec3 C = KeyCenter(K, Layout);
			Cylinder(B.Pcb, C.X + kMxPin1X, C.Y + kMxPin1Y, 0.95, 10, PcbTop - 0.02, PcbTop + 0.05);
			Cylinder(B.Pcb, C.X + kMxPin2X, C.Y + kMxPin2Y, 0.95, 10, PcbTop - 0.02, PcbTop + 0.05);
			// SOD-123 diode below the switch: 2.7 x 1.6 x 1.1 mm body and two leads.
			const double Dx = C.X, Dy = C.Y - 9.6;
			Box(B.PcbParts, Dx - 1.35, Dy - 0.8, PcbTop, Dx + 1.35, Dy + 0.8, PcbTop + 1.1);
			Box(B.PcbParts, Dx - 1.85, Dy - 0.45, PcbTop - 0.01, Dx - 1.2, Dy + 0.45, PcbTop + 0.2);
			Box(B.PcbParts, Dx + 1.2, Dy - 0.45, PcbTop - 0.01, Dx + 1.85, Dy + 0.45, PcbTop + 0.2);
		}
		// Controller: a 7 x 7 mm QFN-44 between the function and number rows, with a crystal and
		// decoupling capacitors.
		double QfnY = Hd - 2.0 * U;
		{
			const mkb::KeyDef* F5 = Layout.Find("F5");
			const mkb::KeyDef* D5 = Layout.Find("D5");
			if (F5 && D5)
			{
				QfnY = 0.5 * (KeyCenter(*F5, Layout).Y + KeyCenter(*D5, Layout).Y) - 2.0;
			}
		}
		const double QfnX = UsbX - 24.0;
		Box(B.PcbParts, QfnX - 3.5, QfnY - 3.5, PcbTop, QfnX + 3.5, QfnY + 3.5, PcbTop + 0.85);
		for (int k = 0; k < 11; ++k)
		{
			const double T = -2.5 + 0.5 * k;
			Box(B.PcbParts, QfnX + T - 0.12, QfnY - 3.6, PcbTop - 0.01, QfnX + T + 0.12, QfnY - 3.0, PcbTop + 0.2);
			Box(B.PcbParts, QfnX + T - 0.12, QfnY + 3.0, PcbTop - 0.01, QfnX + T + 0.12, QfnY + 3.6, PcbTop + 0.2);
			Box(B.PcbParts, QfnX - 3.6, QfnY + T - 0.12, PcbTop - 0.01, QfnX - 3.0, QfnY + T + 0.12, PcbTop + 0.2);
			Box(B.PcbParts, QfnX + 3.0, QfnY + T - 0.12, PcbTop - 0.01, QfnX + 3.6, QfnY + T + 0.12, PcbTop + 0.2);
		}
		Box(B.PcbParts, QfnX + 6.0, QfnY - 1.25, PcbTop, QfnX + 9.2, QfnY + 1.25, PcbTop + 0.8);
		const double Caps[4][2] = {{-5.5, 2.0}, {-5.5, -2.0}, {2.0, -5.5}, {-2.0, 5.5}};
		for (const auto& Q : Caps)
		{
			const bool Along = std::fabs(Q[0]) < std::fabs(Q[1]);
			const double Hx = Along ? 0.8 : 0.4, Hy = Along ? 0.4 : 0.8;
			Box(B.PcbParts, QfnX + Q[0] - Hx, QfnY + Q[1] - Hy, PcbTop, QfnX + Q[0] + Hx, QfnY + Q[1] + Hy, PcbTop + 0.45);
		}
		// USB-C receptacle (8.94 x 7.35 x 3.26 mm shell) with its mouth facing out of the back edge.
		{
			Mesh Usb;
			const Loop2 Shell = RoundedRect(0, 0, 0.5 * UsbW, 0.5 * UsbH, 1.5, 6, SameSplits(Uniform(4), {}));
			const Loop2 Mouth = RoundedRect(0, 0, 0.5 * UsbW - 0.3, 0.5 * UsbH - 0.35, 1.15, 6, SameSplits(Uniform(4), {}));
			Loft(Usb, {AtZ(Shell, 0.0), AtZ(Shell, UsbD)});
			FanCap(Usb, AtZ(Shell, 0.0), Vec3(0, 0, 0), Down);
			StripFlat(Usb, AtZ(Shell, UsbD), AtZ(Mouth, UsbD), Up);
			Loft(Usb, {AtZ(Mouth, UsbD), AtZ(Mouth, 1.0)});
			FanCap(Usb, AtZ(Mouth, 1.0), Vec3(0, 0, 1.0), Up);
			Box(Usb, -3.3, -0.35, 0.9, 3.3, 0.35, UsbD - 0.7);
			Frame3 F;
			F.O = Vec3(UsbX, Hd + EPcb + 0.6 - UsbD, UsbZ);
			F.U = Vec3(1, 0, 0);
			F.V = Vec3(0, 0, -1);
			F.W = Vec3(0, 1, 0);
			ApplyFrame(Usb, F);
			B.PcbParts.Append(Usb);
		}
	}

	// OLED module: a bezel on the plate with a chamfered, tilted top and a pocket holding the glass.
	if (P.OledModule)
	{
		const mkb::Rect& D = Layout.Display;
		const double Dcx = (D.X + 0.5 * D.W - 0.5 * Layout.Width) * U;
		const double Dcy = (0.5 * Layout.Height - (D.Y + 0.5 * D.H)) * U;
		const double BHw = 0.5 * (D.W * U - 0.8), BHd = 0.5 * (D.H * U - 0.8);
		const double Tilt = 8.0 * kDeg, TanT = std::tan(Tilt), Zc = 9.0, Bevel = 0.5;
		const double GHw = 33.0, GHd = 6.6, Clear = 0.2, PocketDepth = 1.2, GlassDrop = 0.15;
		auto TopZ = [Dcy, TanT, Zc](double, double Y) { return Zc + (Y - Dcy) * TanT; };
		auto Below = [&TopZ](double Dz) { return [&TopZ, Dz](double X, double Y) { return TopZ(X, Y) - Dz; }; };
		const Splits BezelSp = SameSplits(Uniform(30), Uniform(5));
		const Loop2 Outer = RoundedRect(Dcx, Dcy, BHw, BHd, 2.0, 6, BezelSp);
		const Loop2 Inset = RoundedRect(Dcx, Dcy, BHw - Bevel, BHd - Bevel, 1.5, 6, BezelSp);
		const Loop2 Pocket = RoundedRect(Dcx, Dcy, GHw + Clear, GHd + Clear, 0.4, 4, SameSplits(Uniform(30), Uniform(5)));
		const Vec3 TopN = Vec3(0, -TanT, 1).Normalized();
		Mesh& Bz = B.DisplayBezel;
		Loft(Bz, {AtZ(Outer, 0.0), Lift(Outer, Below(Bevel))});
		Loft(Bz, {Lift(Outer, Below(Bevel)), Lift(Inset, Below(0.0))});
		FanCap(Bz, AtZ(Outer, 0.0), Vec3(Dcx, Dcy, 0.0), Down);
		Bridge(Bz, Lift(Inset, Below(0.0)), Lift(Pocket, Below(0.0)), Vec3(Dcx, Dcy, TopZ(0, Dcy)), TopN);
		Loft(Bz, {Lift(Pocket, Below(0.0)), Lift(Pocket, Below(PocketDepth))});
		FanCap(Bz, Lift(Pocket, Below(PocketDepth)), Vec3(Dcx, Dcy, TopZ(0, Dcy) - PocketDepth), TopN);

		Mesh& G = B.DisplayGlass;
		const Loop2 Glass = RoundedRect(Dcx, Dcy, GHw, GHd, 0.2, 2, SameSplits(Uniform(30), Uniform(5)));
		const std::vector<Vec3> GTop = Lift(Glass, Below(GlassDrop)), GBot = Lift(Glass, Below(PocketDepth - 0.02));
		Loft(G, {GBot, GTop});
		FanCap(G, GBot, Vec3(Dcx, Dcy, TopZ(0, Dcy) - PocketDepth + 0.02), -TopN);
		auto GlassVertex = [&G, &TopN, Dcx, Dcy, GHw, GHd](const Vec3& Q) {
			const Vec2 UV{(Q.X - (Dcx - GHw)) / (2.0 * GHw), ((Dcy + GHd) - Q.Y) / (2.0 * GHd)};
			return G.AddVertex(Q, TopN, UV);
		};
		const Vec3 GC(Dcx, Dcy, TopZ(0, Dcy) - GlassDrop);
		const int Gc = GlassVertex(GC);
		std::vector<int> Rim;
		for (const Vec3& Q : GTop)
		{
			Rim.push_back(GlassVertex(Q));
		}
		for (size_t i = 0; i < Rim.size(); ++i)
		{
			G.AddTriangle(Gc, Rim[i], Rim[(i + 1) % Rim.size()]);
		}
		B.DisplayCenter = GC;
		B.DisplayWidthMm = 2.0 * GHw;
		B.DisplayHeightMm = 2.0 * GHd / std::cos(Tilt);
		B.DisplayTiltRad = Tilt;
	}
	else
	{
		// Without the OLED, a brass name plate with a raised border and two screws fills the window.
		const mkb::Rect& D = Layout.Display;
		const double Dcx = (D.X + 0.5 * D.W - 0.5 * Layout.Width) * U;
		const double Dcy = (0.5 * Layout.Height - (D.Y + 0.5 * D.H)) * U;
		const double Hx = 0.5 * (D.W * U - 3.0), Hy = 0.5 * (D.H * U - 3.0);
		Mesh& N = B.NamePlate;
		const Loop2 Outer = RoundedRect(Dcx, Dcy, Hx, Hy, 2.0, 5, SameSplits(Uniform(6), {}));
		const Loop2 Lip = RoundedRect(Dcx, Dcy, Hx - 0.6, Hy - 0.6, 1.5, 5, SameSplits(Uniform(6), {}));
		const Loop2 Field = RoundedRect(Dcx, Dcy, Hx - 1.6, Hy - 1.6, 1.0, 5, SameSplits(Uniform(6), {}));
		Loft(N, {AtZ(Outer, -0.05), AtZ(Outer, 1.6), AtZ(Lip, 2.2)});
		Loft(N, {AtZ(Lip, 2.2), AtZ(Field, 2.2), AtZ(Field, 1.7)});
		FanCap(N, AtZ(Field, 1.7), Vec3(Dcx, Dcy, 1.7), Up);
		B.NamePlateCenter = Vec3(Dcx, Dcy, 1.7);
		B.NamePlateWidthMm = 2.0 * (Hx - 1.6);
		B.NamePlateHeightMm = 2.0 * (Hy - 1.6);
		FanCap(N, AtZ(Outer, -0.05), Vec3(Dcx, Dcy, -0.05), Down);
		for (double S = -1.0; S <= 1.0; S += 2.0)
		{
			Mesh Screw;
			Cylinder(Screw, Dcx + S * (Hx - 4.5), Dcy, 1.6, 14, 1.6, 2.6);
			N.Append(Screw);
		}
	}
	return B;
}
} // namespace mgeo
