// Gear-train geometry for the stepped-drum calculating mechanism (see GearParts.h).
//
// Proportions follow the Thomas arithmometer and the Brunsviga: a flat case with the setting
// indicators at the front, the stepped drums lying front-to-back under a glass cover, square shafts
// with sliding pinions above them, and the shifting carriage with the counter and accumulator wheels
// at the back. Gears are standard 20-degree involutes (addendum m, dedendum 1.25 m); the drums use
// module 0.6 with 20 tooth positions (pitch diameter 12 mm) and the pinions 12 teeth, the square
// shafts drive the accumulator wheels 18:15 so a pinion that engages d drum teeth turns its wheel by
// exactly d / 10 of a turn.
#include "GearParts.h"
#include "MeshPrimitives.h"

#include <algorithm>
#include <cmath>

namespace mgeo
{
namespace gearparts_detail
{
using namespace prim;

constexpr double kGearDeg = kPi / 180.0;

// --- Stroke font for the digits: polylines in a 0.6 x 1.0 cell (y up). ---------------------------

void AddArc(std::vector<Pt2>& Out, double Cx, double Cy, double Rx, double Ry, double A0Deg, double A1Deg, int Steps)
{
	for (int k = 0; k <= Steps; ++k)
	{
		const double A = (A0Deg + (A1Deg - A0Deg) * k / Steps) * kGearDeg;
		Out.push_back(Pt2{Cx + Rx * std::cos(A), Cy + Ry * std::sin(A)});
	}
}

std::vector<Pt2> ArcStroke(double Cx, double Cy, double Rx, double Ry, double A0Deg, double A1Deg)
{
	std::vector<Pt2> S;
	AddArc(S, Cx, Cy, Rx, Ry, A0Deg, A1Deg, std::max(3, static_cast<int>(std::fabs(A1Deg - A0Deg) / 12.0)));
	return S;
}

std::vector<std::vector<Pt2>> DigitStrokes(int Digit)
{
	std::vector<std::vector<Pt2>> S;
	auto Line = [&S](double X0, double Y0, double X1, double Y1) { S.push_back({Pt2{X0, Y0}, Pt2{X1, Y1}}); };
	auto Loop = [&S](double Cx, double Cy, double Rx, double Ry) {
		S.push_back(ArcStroke(Cx, Cy, Rx, Ry, 80.0, 280.0));
		S.push_back(ArcStroke(Cx, Cy, Rx, Ry, -100.0, 100.0));
	};
	switch (Digit)
	{
	case 0:
		Loop(0.3, 0.5, 0.27, 0.47);
		break;
	case 1:
		Line(0.33, 0.0, 0.33, 1.0);
		Line(0.33, 1.0, 0.13, 0.82);
		break;
	case 2:
	{
		std::vector<Pt2> T = ArcStroke(0.3, 0.73, 0.25, 0.25, 165.0, -35.0);
		T.push_back(Pt2{0.04, 0.03});
		S.push_back(T);
		Line(0.02, 0.0, 0.58, 0.0);
		break;
	}
	case 3:
		S.push_back(ArcStroke(0.29, 0.76, 0.23, 0.23, 155.0, -90.0));
		S.push_back(ArcStroke(0.29, 0.27, 0.27, 0.27, 90.0, -155.0));
		break;
	case 4:
		Line(0.44, 0.0, 0.44, 1.0);
		Line(0.44, 1.0, 0.02, 0.3);
		Line(0.0, 0.3, 0.6, 0.3);
		break;
	case 5:
		Line(0.56, 1.0, 0.11, 1.0);
		Line(0.11, 1.0, 0.08, 0.56);
		S.push_back(ArcStroke(0.3, 0.31, 0.28, 0.3, 128.0, -150.0));
		break;
	case 6:
		S.push_back(ArcStroke(0.33, 0.47, 0.28, 0.5, 62.0, 195.0));
		Loop(0.3, 0.3, 0.27, 0.29);
		break;
	case 7:
		Line(0.02, 1.0, 0.58, 1.0);
		Line(0.58, 1.0, 0.2, 0.0);
		break;
	case 8:
		Loop(0.3, 0.76, 0.22, 0.23);
		Loop(0.3, 0.27, 0.27, 0.27);
		break;
	case 10:
		Line(0.3, 0.2, 0.3, 0.8);
		Line(0.0, 0.5, 0.6, 0.5);
		break;
	case 11:
		Line(0.02, 0.5, 0.58, 0.5);
		break;
	case 9:
		for (std::vector<Pt2> Stroke : DigitStrokes(6))
		{
			for (Pt2& Q : Stroke)
			{
				Q = Pt2{0.6 - Q.X, 1.0 - Q.Y};
			}
			S.push_back(Stroke);
		}
		break;
	default:
		break;
	}
	return S;
}

std::vector<Pt2> Subdivide(const std::vector<Pt2>& Pts, double MaxLen)
{
	std::vector<Pt2> Out;
	for (size_t i = 0; i + 1 < Pts.size(); ++i)
	{
		const double Dx = Pts[i + 1].X - Pts[i].X, Dy = Pts[i + 1].Y - Pts[i].Y;
		const int N = std::max(1, static_cast<int>(std::ceil(std::sqrt(Dx * Dx + Dy * Dy) / MaxLen)));
		for (int k = 0; k < N; ++k)
		{
			Out.push_back(Pt2{Pts[i].X + Dx * k / N, Pts[i].Y + Dy * k / N});
		}
	}
	Out.push_back(Pts.back());
	return Out;
}

// A polygon (in glyph-plane coordinates) raised into a closed solid by Map(a, b, h), h from Lo to Hi;
// Out(a, b) is the outward direction of the surface the glyph sits on.
template <typename MapFn, typename OutFn>
void RaisePolygon(Mesh& M, std::vector<Pt2> Poly, double Lo, double Hi, MapFn Map, OutFn Out)
{
	if (SignedArea2(Poly) < 0.0)
	{
		std::reverse(Poly.begin(), Poly.end());
	}
	const std::vector<int> Tris = EarClip(Poly);
	for (size_t t = 0; t + 2 < Tris.size(); t += 3)
	{
		const Pt2 &A = Poly[static_cast<size_t>(Tris[t])], &B = Poly[static_cast<size_t>(Tris[t + 1])], &C = Poly[static_cast<size_t>(Tris[t + 2])];
		const Vec3 N = Out((A.X + B.X + C.X) / 3.0, (A.Y + B.Y + C.Y) / 3.0);
		// Each facet takes its own normal: on a curved rim a thin triangle's plane can tilt well away
		// from the surface normal.
		FlatTriFacet(M, Map(A.X, A.Y, Hi), Map(B.X, B.Y, Hi), Map(C.X, C.Y, Hi), N);
		FlatTriFacet(M, Map(A.X, A.Y, Lo), Map(C.X, C.Y, Lo), Map(B.X, B.Y, Lo), -N);
	}
	for (size_t i = 0; i < Poly.size(); ++i)
	{
		const Pt2& P0 = Poly[i];
		const Pt2& P1 = Poly[(i + 1) % Poly.size()];
		const double Mx = 0.5 * (P0.X + P1.X), My = 0.5 * (P0.Y + P1.Y);
		const double Ex = P1.Y - P0.Y, Ey = -(P1.X - P0.X);
		const double H = 0.5 * (Lo + Hi);
		const Vec3 Side = Map(Mx + Ex * 0.01, My + Ey * 0.01, H) - Map(Mx, My, H);
		const Vec3 A = Map(P0.X, P0.Y, Lo), B = Map(P1.X, P1.Y, Lo), C = Map(P1.X, P1.Y, Hi), D = Map(P0.X, P0.Y, Hi);
		FlatTri(M, A, B, C, Side);
		FlatTri(M, A, C, D, Side);
	}
}

// Strokes of one digit in millimetres, centred on the glyph's middle: a across, b up.
std::vector<std::vector<Pt2>> GlyphPolygons(int Digit, double Width, double Height, double Stroke)
{
	std::vector<std::vector<Pt2>> Out;
	for (const std::vector<Pt2>& S : DigitStrokes(Digit))
	{
		std::vector<Pt2> Mm;
		for (const Pt2& Q : S)
		{
			Mm.push_back(Pt2{(Q.X - 0.3) / 0.6 * Width, (Q.Y - 0.5) * Height});
		}
		// Split at sharp turns: a thickened polyline would fold over itself there.
		std::vector<Pt2> Run = {Mm[0]};
		for (size_t i = 1; i < Mm.size(); ++i)
		{
			Run.push_back(Mm[i]);
			if (i + 1 < Mm.size())
			{
				const double Ax = Mm[i].X - Mm[i - 1].X, Ay = Mm[i].Y - Mm[i - 1].Y, Bx = Mm[i + 1].X - Mm[i].X, By = Mm[i + 1].Y - Mm[i].Y;
				const double Cos = (Ax * Bx + Ay * By) / std::max(1e-12, std::sqrt((Ax * Ax + Ay * Ay) * (Bx * Bx + By * By)));
				if (Cos < 0.6)
				{
					Out.push_back(ThickenPolyline(Subdivide(Run, std::max(0.25, Height / 9.0)), 0.5 * Stroke));
					Run = {Mm[i]};
				}
			}
		}
		Out.push_back(ThickenPolyline(Subdivide(Run, std::max(0.25, Height / 9.0)), 0.5 * Stroke));
	}
	return Out;
}

// --- Involute gear outline -----------------------------------------------------------------------

struct GearOutline
{
	std::vector<Pt2> P;
	std::vector<char> Sharp;
	std::vector<int> ToothStart;   // first outline index of each tooth (its right flank's root)
};

double Involute(double A)
{
	return std::tan(A) - A;
}

GearOutline InvoluteOutline(int Teeth, double Module, double AlphaDeg, int FlankPoints)
{
	GearOutline O;
	const double R = 0.5 * Module * Teeth, Alpha = AlphaDeg * kGearDeg;
	const double Rb = R * std::cos(Alpha), Ra = R + Module, Rf = R - 1.25 * Module;
	const double Psi = kPi / (2.0 * Teeth);
	auto Theta = [&](double Rho) { return Psi + Involute(Alpha) - Involute(std::acos(std::min(1.0, Rb / Rho))); };
	const double R0 = std::max(Rb, Rf);
	const double ThetaRoot = Theta(R0);
	const double ThetaTip = std::max(0.15 * Psi, Theta(Ra));
	auto Push = [&O](double Rho, double A, bool Sharp) {
		O.P.push_back(Pt2{Rho * std::cos(A), Rho * std::sin(A)});
		O.Sharp.push_back(Sharp ? 1 : 0);
	};
	const int F = std::max(2, FlankPoints);
	for (int i = 0; i < Teeth; ++i)
	{
		const double G = 2.0 * kPi * i / Teeth;
		O.ToothStart.push_back(static_cast<int>(O.P.size()));
		if (Rf < Rb)
		{
			Push(Rf, G - ThetaRoot, true);
		}
		for (int k = 0; k <= F; ++k)
		{
			const double Rho = R0 + (Ra - R0) * k / F;
			Push(Rho, G - (k == F ? ThetaTip : Theta(Rho)), k == 0 || k == F);
		}
		Push(Ra, G, false);
		for (int k = F; k >= 0; --k)
		{
			const double Rho = R0 + (Ra - R0) * k / F;
			Push(Rho, G + (k == F ? ThetaTip : Theta(Rho)), k == 0 || k == F);
		}
		if (Rf < Rb)
		{
			Push(Rf, G + ThetaRoot, true);
		}
		const double A0 = G + ThetaRoot, A1 = G + 2.0 * kPi / Teeth - ThetaRoot;
		for (int k = 1; k <= 2; ++k)
		{
			Push(Rf, A0 + (A1 - A0) * k / 3.0, false);
		}
	}
	return O;
}

// A circle of radius Rm with a vertex midway between each pair of angularly distinct outline points,
// so a face bridging a toothed outline to it never needs a triangle spanning a radial flank.
std::vector<Pt2> MidAngleRing(const std::vector<Pt2>& Outline, double Rm)
{
	std::vector<double> Ang;
	for (const Pt2& Q : Outline)
	{
		const double A = std::atan2(Q.Y, Q.X);
		if (Ang.empty() || std::fabs(std::remainder(A - Ang.back(), 2.0 * kPi)) > 1e-7)
		{
			Ang.push_back(A);
		}
	}
	if (Ang.size() > 1 && std::fabs(std::remainder(Ang.front() - Ang.back(), 2.0 * kPi)) < 1e-7)
	{
		Ang.pop_back();
	}
	std::vector<Pt2> Out;
	for (size_t i = 0; i < Ang.size(); ++i)
	{
		const double A0 = Ang[i];
		const double A = A0 + 0.5 * std::remainder(Ang[(i + 1) % Ang.size()] - A0, 2.0 * kPi);
		Out.push_back(Pt2{Rm * std::cos(A), Rm * std::sin(A)});
	}
	return Out;
}

double MaxRadius(const std::vector<Pt2>& P)
{
	double R = 0.0;
	for (const Pt2& Q : P)
	{
		R = std::max(R, std::sqrt(Q.X * Q.X + Q.Y * Q.Y));
	}
	return R;
}

std::vector<Vec3> RingAtZ(const std::vector<Pt2>& P, double Z)
{
	std::vector<Vec3> Out;
	for (const Pt2& Q : P)
	{
		Out.push_back(Vec3(Q.X, Q.Y, Z));
	}
	return Out;
}

// The flat face at Z between a toothed outline (root radius Rf) and an inner loop, facing N.
void ToothedFace(Mesh& M, const std::vector<Pt2>& Outline, double Rf, const std::vector<Pt2>& Inner, double Z, const Vec3& N)
{
	const double InnerR = MaxRadius(Inner);
	const double Rm = Rf - 0.25 * (Rf - InnerR);
	const std::vector<Pt2> Mid = MidAngleRing(Outline, Rm);
	Bridge(M, RingAtZ(Outline, Z), RingAtZ(Mid, Z), Vec3(0, 0, Z), N);
	Bridge(M, RingAtZ(Mid, Z), RingAtZ(Inner, Z), Vec3(0, 0, Z), N);
}

// A surface of revolution about +Z from a profile of (radius, z) points running counter-clockwise in
// the (r, z) half-plane. Closed profiles make a ring; open ones must start and end on the axis.
// Sharp marks profile points where the surface creases.
void Lathe(Mesh& M, const std::vector<Pt2>& Profile, int Segments, bool Closed, const std::vector<char>& Sharp)
{
	const size_t N = Profile.size();
	auto Ring = [&Profile, Segments](size_t i) {
		std::vector<Vec3> Out;
		for (int s = 0; s < Segments; ++s)
		{
			const double A = 2.0 * kPi * s / Segments;
			Out.push_back(Vec3(Profile[i].X * std::cos(A), Profile[i].X * std::sin(A), Profile[i].Y));
		}
		return Out;
	};
	const size_t First = Closed ? 0 : 1, Last = Closed ? N : N - 2;
	std::vector<std::vector<Vec3>> Run;
	for (size_t i = First; i <= Last; ++i)
	{
		const size_t Idx = i % N;
		Run.push_back(Ring(Idx));
		const bool End = i == Last;
		if ((Sharp[Idx] && Run.size() > 1) || End)
		{
			if (Run.size() > 1)
			{
				Loft(M, Run);
			}
			Run.clear();
			Run.push_back(Ring(Idx));
		}
	}
	if (!Closed)
	{
		FanCap(M, Ring(1), Vec3(0, 0, Profile[0].Y), Vec3(0, 0, -1));
		FanCap(M, Ring(N - 2), Vec3(0, 0, Profile[N - 1].Y), Vec3(0, 0, 1));
	}
}

std::vector<char> SharpAll(size_t N)
{
	return std::vector<char>(N, 1);
}

// A right-handed frame with W along Axis.
Frame3 FrameAlong(const Vec3& Origin, const Vec3& Axis)
{
	Frame3 F;
	F.O = Origin;
	F.W = Axis.Normalized();
	const Vec3 Ref = std::fabs(F.W.Z) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
	F.U = Ref.Cross(F.W).Normalized();
	F.V = F.W.Cross(F.U);
	return F;
}

void AddInFrame(Mesh& Dst, Mesh Src, const Frame3& F)
{
	ApplyFrame(Src, F);
	Dst.Append(Src);
}

Pose PoseOf(const Vec3& Origin, const Vec3& Z, const Vec3& Y)
{
	Pose P;
	P.Origin = Origin;
	P.Z = Z.Normalized();
	P.Y = Y.Normalized();
	P.X = P.Y.Cross(P.Z);
	return P;
}

Frame3 FrameOf(const Pose& P)
{
	Frame3 F;
	F.O = P.Origin;
	F.U = P.X;
	F.V = P.Y;
	F.W = P.Z;
	return F;
}
} // namespace gearparts_detail

// --- Primitives -----------------------------------------------------------------------------------

GearCircles SpurGearCircles(const SpurGearParams& P)
{
	GearCircles C;
	C.Pitch = 0.5 * P.Module * P.Teeth;
	C.Base = C.Pitch * std::cos(P.PressureAngleDeg * gearparts_detail::kGearDeg);
	C.Tip = C.Pitch + P.Module;
	C.Root = C.Pitch - 1.25 * P.Module;
	return C;
}

Mesh BuildSpurGear(const SpurGearParams& P)
{
	using namespace gearparts_detail;
	Mesh M;
	const GearOutline O = InvoluteOutline(P.Teeth, P.Module, P.PressureAngleDeg, P.FlankPoints);
	const double Z0 = -0.5 * P.FaceWidth, Z1 = 0.5 * P.FaceWidth, Rf = SpurGearCircles(P).Root;
	const Vec3 Up(0, 0, 1), Down(0, 0, -1), C0(0, 0, Z0);
	Loft(M, {RingAtZ(O.P, Z0), RingAtZ(O.P, Z1)}, &O.Sharp);
	const bool Bored = P.BoreRadius > 0.0;
	std::vector<Pt2> Bore;
	if (Bored)
	{
		Bore = P.SquareBore ? RectPoly(-P.BoreRadius, -P.BoreRadius, P.BoreRadius, P.BoreRadius) : CircleLoop(0, 0, P.BoreRadius, std::max(12, static_cast<int>(P.BoreRadius * 10.0))).P;
	}
	const std::vector<char> BoreSharp(Bore.size(), P.SquareBore ? 1 : 0);
	double TopZ = Z1;
	if (P.HubRadius > 0.0 && P.HubLength > 0.0)
	{
		const std::vector<Pt2> Hub = CircleLoop(0, 0, P.HubRadius, 24).P;
		ToothedFace(M, O.P, Rf, Hub, Z1, Up);
		Loft(M, {RingAtZ(Hub, Z1), RingAtZ(Hub, Z1 + P.HubLength)});
		TopZ = Z1 + P.HubLength;
		if (Bored)
		{
			Bridge(M, RingAtZ(Hub, TopZ), RingAtZ(Bore, TopZ), Vec3(0, 0, TopZ), Up);
		}
		else
		{
			FanCap(M, RingAtZ(Hub, TopZ), Vec3(0, 0, TopZ), Up);
		}
	}
	else if (Bored)
	{
		ToothedFace(M, O.P, Rf, Bore, Z1, Up);
	}
	else
	{
		FanCap(M, RingAtZ(O.P, Z1), Vec3(0, 0, Z1), Up);
	}
	if (Bored)
	{
		ToothedFace(M, O.P, Rf, Bore, Z0, Down);
		Loft(M, {RingAtZ(Bore, TopZ), RingAtZ(Bore, Z0)}, &BoreSharp);
	}
	else
	{
		FanCap(M, RingAtZ(O.P, Z0), C0, Down);
	}
	return M;
}

Mesh BuildBevelGear(const BevelGearParams& P)
{
	using namespace gearparts_detail;
	Mesh M;
	const GearOutline O = InvoluteOutline(P.Teeth, P.Module, P.PressureAngleDeg, P.FlankPoints);
	const double Delta = P.PitchAngleDeg * kGearDeg;
	const double R = 0.5 * P.Module * P.Teeth, Zh = R / std::tan(Delta), L = R / std::sin(Delta);
	const double K = 1.0 - std::min(0.6 * L, P.FaceWidth) / L;
	// The pitch-plane outline mapped onto the back cone.
	auto OnBackCone = [R, Zh, Delta](const std::vector<Pt2>& Pts) {
		std::vector<Vec3> Out;
		for (const Pt2& Q : Pts)
		{
			const double Rho = std::sqrt(Q.X * Q.X + Q.Y * Q.Y), Phi = std::atan2(Q.Y, Q.X);
			const double S = Rho - R;
			const double R3 = R + S * std::cos(Delta), Z = Zh - S * std::sin(Delta);
			Out.push_back(Vec3(R3 * std::cos(Phi), R3 * std::sin(Phi), Z));
		}
		return Out;
	};
	auto Scaled = [](std::vector<Vec3> V, double F) {
		for (Vec3& Q : V)
		{
			Q = Q * F;
		}
		return V;
	};
	const std::vector<Vec3> Heel = OnBackCone(O.P), Toe = Scaled(Heel, K);
	Loft(M, {Toe, Heel}, &O.Sharp);
	const double Rf = R - 1.25 * P.Module;
	const double RIn = P.HubRadius > 0.0 ? P.HubRadius : std::max(P.BoreRadius + 0.5 * P.Module, 0.45 * Rf);
	const double ZIn = Zh + (R - RIn) * std::tan(Delta);
	const std::vector<Pt2> Inner = CircleLoop(0, 0, RIn, 24).P;
	const std::vector<Vec3> MidHeel = OnBackCone(MidAngleRing(O.P, Rf - 0.25 * (Rf - RIn))), MidToe = Scaled(MidHeel, K);
	const Vec3 Up(0, 0, 1), Down(0, 0, -1);
	Bridge(M, Heel, MidHeel, Vec3(0, 0, MidHeel[0].Z), Up, true);
	Bridge(M, MidHeel, RingAtZ(Inner, ZIn), Vec3(0, 0, ZIn), Up, true);
	double EndZ = ZIn;
	if (P.HubLength > 0.0)
	{
		Loft(M, {RingAtZ(Inner, ZIn), RingAtZ(Inner, ZIn + P.HubLength)});
		EndZ = ZIn + P.HubLength;
	}
	const double ToeZ = K * Zh;
	if (P.BoreRadius > 0.0)
	{
		const std::vector<Pt2> Bore = CircleLoop(0, 0, P.BoreRadius, 16).P;
		Bridge(M, RingAtZ(Inner, EndZ), RingAtZ(Bore, EndZ), Vec3(0, 0, EndZ), Up);
		Loft(M, {RingAtZ(Bore, EndZ), RingAtZ(Bore, ToeZ)});
		Bridge(M, Toe, MidToe, Vec3(0, 0, MidToe[0].Z), Down, true);
		Bridge(M, MidToe, RingAtZ(Bore, ToeZ), Vec3(0, 0, ToeZ), Down, true);
	}
	else
	{
		FanCap(M, RingAtZ(Inner, EndZ), Vec3(0, 0, EndZ), Up);
		Bridge(M, Toe, MidToe, Vec3(0, 0, MidToe[0].Z), Down, true);
		FanCap(M, MidToe, Vec3(0, 0, MidToe[0].Z), Down);
	}
	return M;
}

double StepDrumLayout::ToothStartZ(int K) const
{
	return TeethEndZ - K * StepMm;
}

double StepDrumLayout::ToothAngleRad(int K) const
{
	return 2.0 * prim::kPi * (K - 1) / ToothPositions;
}

double StepDrumLayout::PinionCenterZ(int Digit) const
{
	return TeethEndZ - (9.5 - Digit) * StepMm;
}

StepDrumLayout StepDrumGeometry(const StepDrumParams& P)
{
	StepDrumLayout L;
	L.PitchRadius = 0.5 * P.Module * P.ToothPositions;
	L.RootRadius = L.PitchRadius - 1.25 * P.Module;
	L.TipRadius = L.PitchRadius + P.Module;
	L.StepMm = P.StepMm;
	L.ToothPositions = P.ToothPositions;
	L.BodyStartZ = P.StartZ;
	L.TeethEndZ = P.StartZ + 1.2 + 0.5 * P.PinionFaceWidth + 9.5 * P.StepMm;
	L.BodyEndZ = L.TeethEndZ + 0.8;
	return L;
}

Mesh BuildStepDrum(const StepDrumParams& P)
{
	using namespace gearparts_detail;
	Mesh M;
	const StepDrumLayout L = StepDrumGeometry(P);
	const double Flange = L.TipRadius + 0.3;   // neighbouring drums stand one pitch apart
	// Core with a flange at each end.
	const std::vector<Pt2> Core = {{0.0, L.BodyStartZ}, {Flange, L.BodyStartZ}, {Flange, L.BodyStartZ + 0.8}, {L.RootRadius, L.BodyStartZ + 0.8}, {L.RootRadius, L.TeethEndZ},
		{Flange, L.TeethEndZ}, {Flange, L.BodyEndZ}, {0.0, L.BodyEndZ}};
	Lathe(M, Core, 32, false, SharpAll(Core.size()));
	// Teeth: one involute tooth of the ToothPositions gear, extruded from its step to the flush end.
	const GearOutline O = InvoluteOutline(P.ToothPositions, P.Module, 20.0, 4);
	const int First = O.ToothStart[0], Next = O.ToothStart[1];
	std::vector<Pt2> Tooth(O.P.begin() + First, O.P.begin() + Next - 2);
	const double Ar = std::atan2(Tooth.front().Y, Tooth.front().X), Al = std::atan2(Tooth.back().Y, Tooth.back().X);
	const double Rin = L.RootRadius - 0.5;
	Tooth.push_back(Pt2{Rin * std::cos(Al), Rin * std::sin(Al)});
	Tooth.insert(Tooth.begin(), Pt2{Rin * std::cos(Ar), Rin * std::sin(Ar)});
	for (int K = 1; K <= P.Teeth; ++K)
	{
		const double A = L.ToothAngleRad(K), C = std::cos(A), S = std::sin(A);
		std::vector<Pt2> T;
		for (const Pt2& Q : Tooth)
		{
			T.push_back(Pt2{Q.X * C - Q.Y * S, Q.X * S + Q.Y * C});
		}
		Prism(M, T, L.ToothStartZ(K), L.TeethEndZ + 0.05, Frame3());
	}
	if (P.Bevel)
	{
		BevelGearParams B;
		B.Teeth = P.BevelTeeth;
		B.Module = P.Module;
		B.FaceWidth = 2.0;
		B.HubRadius = 0.45 * L.RootRadius;
		M.Append(BuildBevelGear(B));
	}
	return M;
}

double DigitAngleRad(int Digit)
{
	return Digit * 2.0 * prim::kPi / 10.0;
}

NumeralWheel BuildNumeralWheel(const NumeralWheelParams& P)
{
	using namespace gearparts_detail;
	NumeralWheel W;
	const double R = P.Radius, H = 0.5 * P.Width, C = std::min(0.4, 0.2 * P.Width);
	const double Bore = std::max(0.3, P.BoreRadius);
	const std::vector<Pt2> Rim = {{Bore, -H}, {R - C, -H}, {R, -H + C}, {R, H - C}, {R - C, H}, {Bore, H}};
	Lathe(W.Wheel, Rim, 48, true, SharpAll(Rim.size()));
	// Gears: the drive gear on -z, the star (and carry tooth) on +z or on the front.
	double FrontZ = -H, BackZ = H;
	if (P.GearTeeth > 0)
	{
		SpurGearParams G;
		G.Teeth = P.GearTeeth;
		G.Module = P.GearModule;
		G.FaceWidth = 1.6;
		G.BoreRadius = Bore;
		if (P.GearGap > 0.0)
		{
			Lathe(W.Gears, {{Bore, FrontZ - P.GearGap - 0.1}, {Bore + 1.3, FrontZ - P.GearGap - 0.1}, {Bore + 1.3, FrontZ + 0.05}, {Bore, FrontZ + 0.05}}, 16, true, SharpAll(4));
			FrontZ -= P.GearGap;
		}
		Mesh Gm = BuildSpurGear(G);
		Gm.Translate(Vec3(0, 0, FrontZ - 0.8 + 0.05));
		W.Gears.Append(Gm);
		FrontZ -= 1.55;
	}
	if (P.StarTeeth > 0)
	{
		SpurGearParams G;
		G.Teeth = P.StarTeeth;
		G.Module = P.StarModule;
		G.FaceWidth = 1.2;
		G.BoreRadius = Bore;
		Mesh Gm = BuildSpurGear(G);
		if (P.StarOnFront)
		{
			Gm.Translate(Vec3(0, 0, FrontZ - 0.6 + 0.05));
			FrontZ -= 1.15;
		}
		else
		{
			Gm.Translate(Vec3(0, 0, BackZ + 0.6 - 0.05));
			BackZ += 1.15;
		}
		W.Gears.Append(Gm);
	}
	if (P.CarryTooth)
	{
		// A pin on the face opposite the star, outside the gear, half way between digits 9 and 0, so it
		// trips the carry lever as the wheel goes from 9 to 0.
		W.CarryToothAngleRad = -0.5 * kPi + 0.5 * DigitAngleRad(1);
		const double Rc = std::min(R - 1.0, std::max(3.0, P.GearTeeth > 0 ? 0.5 * P.GearModule * P.GearTeeth + P.GearModule + 0.45 : 3.0));
		Mesh Pin;
		const bool OnFront = P.StarTeeth > 0 && !P.StarOnFront;
		if (OnFront)
		{
			Cylinder(Pin, Rc, 0.0, 0.35, 10, FrontZ - 1.4, -H + 0.1);
		}
		else
		{
			Cylinder(Pin, Rc, 0.0, 0.35, 10, H - 0.1, BackZ + 1.4);
		}
		RotateAboutZ(Pin, W.CarryToothAngleRad);
		W.Gears.Append(Pin);
	}
	// Glyphs, embossed on the rim.
	const double Lo = -0.12, Hi = P.EmbossMm;
	for (int D = 0; D < 10; ++D)
	{
		W.GlyphTriangleStart.push_back(W.Glyphs.TriangleCount());
		const int Glyph = P.SignWheel ? (D % 2 == 0 ? 10 : 11) : D;
		const double Phi = 0.5 * kPi - DigitAngleRad(D);
		const bool Along = P.GlyphsAlongAxis;
		auto Map = [R, Phi, Along](double A, double B, double Hh) {
			const double F = Along ? Phi + A / R : Phi - B / R;
			const double Z = Along ? B : A;
			return Vec3((R + Hh) * std::cos(F), (R + Hh) * std::sin(F), Z);
		};
		auto Out = [R, Phi, Along](double A, double B) {
			const double F = Along ? Phi + A / R : Phi - B / R;
			return Vec3(std::cos(F), std::sin(F), 0.0);
		};
		for (const std::vector<Pt2>& Poly : GlyphPolygons(Glyph, P.GlyphWidth, P.GlyphHeight, P.StrokeMm))
		{
			RaisePolygon(W.Glyphs, Poly, Lo, Hi, Map, Out);
		}
	}
	W.GlyphTriangleStart.push_back(W.Glyphs.TriangleCount());
	return W;
}

Mesh BuildShaft(double Radius, double Length, bool Square)
{
	using namespace gearparts_detail;
	Mesh M;
	const double C = std::min(0.3 * Radius, 0.25 * Length);
	if (Square)
	{
		const double A = Radius, B = Radius - C;
		const std::vector<Pt2> Oct = {{A, -B}, {A, B}, {B, A}, {-B, A}, {-A, B}, {-A, -B}, {-B, -A}, {B, -A}};
		Prism(M, Oct, 0.0, Length, Frame3());
		return M;
	}
	const std::vector<Pt2> Profile = {{0.0, 0.0}, {Radius - C, 0.0}, {Radius, C}, {Radius, Length - C}, {Radius - C, Length}, {0.0, Length}};
	Lathe(M, Profile, std::max(12, static_cast<int>(Radius * 8.0)), false, SharpAll(Profile.size()));
	return M;
}

Mesh BuildBearingBlock(double BoreRadius, double AxisHeight, double Width)
{
	using namespace gearparts_detail;
	Mesh M;
	const double Ro = BoreRadius + std::max(1.5, 0.8 * BoreRadius);
	const double Foot = 2.0 * Ro + 6.0;
	// Base with bolt bosses, a web, and the bearing ring around the bore (axis +X).
	SolidFromLoop(M, RoundedRect(0, 0, 0.5 * Width + 0.5, 0.5 * Foot, 1.0, 3, Splits()), 0.0, std::min(2.0, 0.3 * AxisHeight));
	Box(M, -0.5 * Width, -0.6 * Ro, 0.5, 0.5 * Width, 0.6 * Ro, AxisHeight);
	for (double S = -1.0; S <= 1.0; S += 2.0)
	{
		Cylinder(M, 0.0, S * (0.5 * Foot - 1.6), 1.1, 10, 0.0, std::min(2.0, 0.3 * AxisHeight) + 0.8);
	}
	Mesh Ring;
	const std::vector<Pt2> Profile = {{BoreRadius, -0.5 * Width}, {Ro, -0.5 * Width}, {Ro, 0.5 * Width}, {BoreRadius, 0.5 * Width}};
	Lathe(Ring, Profile, 24, true, SharpAll(Profile.size()));
	Frame3 F;
	F.O = Vec3(0, 0, AxisHeight);
	F.U = Vec3(0, 1, 0);
	F.V = Vec3(0, 0, 1);
	F.W = Vec3(1, 0, 0);
	ApplyFrame(Ring, F);
	M.Append(Ring);
	return M;
}

Mesh BuildRack(double Module, int Teeth, double FaceWidth, double Body, double PressureAngleDeg)
{
	using namespace gearparts_detail;
	Mesh M;
	const double Pitch = kPi * Module, T = std::tan(PressureAngleDeg * kGearDeg);
	const double Root = -1.25 * Module, Bottom = Root - Body, Len = Teeth * Pitch;
	std::vector<Pt2> Poly = {{0.0, Bottom}, {Len, Bottom}, {Len, Root}};
	for (int i = Teeth - 1; i >= 0; --i)
	{
		const double C = (i + 0.5) * Pitch;
		Poly.push_back(Pt2{C + 0.25 * Pitch + 1.25 * Module * T, Root});
		Poly.push_back(Pt2{C + 0.25 * Pitch - Module * T, Module});
		Poly.push_back(Pt2{C - 0.25 * Pitch + Module * T, Module});
		Poly.push_back(Pt2{C - 0.25 * Pitch - 1.25 * Module * T, Root});
	}
	Poly.push_back(Pt2{0.0, Root});
	Prism(M, Poly, -0.5 * FaceWidth, 0.5 * FaceWidth, Frame3());
	return M;
}

// --- Parts and poses -------------------------------------------------------------------------------

const char* MechMaterialName(MechMaterial M)
{
	switch (M)
	{
	case MechMaterial::Brass:
		return "brass";
	case MechMaterial::Steel:
		return "steel";
	case MechMaterial::CastIron:
		return "cast_iron";
	case MechMaterial::BlackEnamel:
		return "black_enamel";
	case MechMaterial::NumeralBlack:
		return "numeral_black";
	case MechMaterial::Paint:
		return "paint";
	case MechMaterial::Engraving:
		return "engraving";
	case MechMaterial::WhiteEnamel:
		return "white_enamel";
	case MechMaterial::Glass:
		return "glass";
	case MechMaterial::Paper:
		return "paper";
	case MechMaterial::Rubber:
		return "rubber";
	case MechMaterial::Ribbon:
		return "ribbon";
	case MechMaterial::Wood:
		return "wood";
	}
	return "brass";
}

Mesh& MechPart::Of(MechMaterial M)
{
	for (MechSection& S : Sections)
	{
		if (S.Material == M)
		{
			return S.Geometry;
		}
	}
	// Room for every material up front, so references handed out earlier stay valid.
	if (Sections.capacity() < 16)
	{
		Sections.reserve(16);
	}
	Sections.push_back(MechSection());
	Sections.back().Material = M;
	return Sections.back().Geometry;
}

int MechPart::TriangleCount() const
{
	int N = 0;
	for (const MechSection& S : Sections)
	{
		N += S.Geometry.TriangleCount();
	}
	return N;
}

Pose Pose::Then(const Pose& Child) const
{
	Pose P;
	P.Origin = Point(Child.Origin);
	P.X = Direction(Child.X);
	P.Y = Direction(Child.Y);
	P.Z = Direction(Child.Z);
	return P;
}

Pose Pose::RotatedZ(double Angle) const
{
	Pose P = *this;
	const double C = std::cos(Angle), S = std::sin(Angle);
	P.X = X * C + Y * S;
	P.Y = Y * C - X * S;
	return P;
}

Pose Pose::Moved(const Vec3& LocalOffset) const
{
	Pose P = *this;
	P.Origin = Origin + Direction(LocalOffset);
	return P;
}
Mesh BuildDigitText(const std::string& Digits, double Height, double StrokeMm, double EmbossMm)
{
	using namespace gearparts_detail;
	Mesh M;
	const double Width = 0.6 * Height, Advance = 0.95 * Height;
	const double X0 = -0.5 * Advance * (static_cast<double>(Digits.size()) - 1.0);
	for (size_t i = 0; i < Digits.size(); ++i)
	{
		const char Ch = Digits[i];
		const int Glyph = Ch >= '0' && Ch <= '9' ? Ch - '0' : Ch == '+' ? 10 : Ch == '-' ? 11 : -1;
		if (Glyph < 0)
		{
			continue;
		}
		const double Cx = X0 + Advance * static_cast<double>(i);
		auto Map = [Cx](double A, double B, double H) { return Vec3(Cx + A, B, H); };
		auto Out = [](double, double) { return Vec3(0, 0, 1); };
		for (const std::vector<Pt2>& Poly : GlyphPolygons(Glyph, Width, Height, StrokeMm))
		{
			RaisePolygon(M, Poly, -0.1, EmbossMm, Map, Out);
		}
	}
	return M;
}

Mesh BuildSpokedWheel(const SpokedWheelParams& P)
{
	using namespace gearparts_detail;
	Mesh M;
	const double Z0 = -0.5 * P.FaceWidth, Z1 = 0.5 * P.FaceWidth;
	double Rin = 0.0;
	if (P.Teeth > 0)
	{
		const GearOutline O = InvoluteOutline(P.Teeth, P.Module, 20.0, 4);
		Rin = 0.5 * P.Module * P.Teeth - 1.25 * P.Module - P.RimDepth;
		const std::vector<Pt2> Inner = CircleLoop(0, 0, Rin, 72).P;
		Loft(M, {RingAtZ(O.P, Z0), RingAtZ(O.P, Z1)}, &O.Sharp);
		const double Rf = 0.5 * P.Module * P.Teeth - 1.25 * P.Module;
		ToothedFace(M, O.P, Rf, Inner, Z1, Vec3(0, 0, 1));
		ToothedFace(M, O.P, Rf, Inner, Z0, Vec3(0, 0, -1));
		Loft(M, {RingAtZ(Inner, Z1), RingAtZ(Inner, Z0)});
	}
	else
	{
		Rin = P.RimRadius - P.RimDepth;
		const double R = P.RimRadius, C = 0.25 * P.RimDepth;
		Lathe(M, {{Rin, Z0}, {R - C, Z0}, {R, Z0 + C}, {R, Z1 - C}, {R - C, Z1}, {Rin, Z1}}, 72, true, SharpAll(6));
	}
	const double Hh = 0.5 * P.HubLength;
	Lathe(M, {{P.BoreRadius, -Hh}, {P.HubRadius, -Hh}, {P.HubRadius, Hh}, {P.BoreRadius, Hh}}, 32, true, SharpAll(4));
	const double Ht = 0.3 * P.FaceWidth;
	for (int K = 0; K < P.Spokes; ++K)
	{
		const double A = 2.0 * kPi * K / P.Spokes;
		const double R0 = P.HubRadius - 0.6, R1 = Rin + 0.6, W0 = 0.6 * P.SpokeWidth, W1 = 0.5 * P.SpokeWidth;
		std::vector<Pt2> Spoke = {{R0, -W0}, {R1, -W1}, {R1, W1}, {R0, W0}};
		for (Pt2& Q : Spoke)
		{
			Q = Pt2{Q.X * std::cos(A) - Q.Y * std::sin(A), Q.X * std::sin(A) + Q.Y * std::cos(A)};
		}
		Prism(M, Spoke, -Ht, Ht, Frame3());
	}
	return M;
}

namespace gearparts_detail
{
// Holes are bridged into the outer outline (keyhole cuts), so ear clipping sees one simple polygon.
std::vector<Pt2> KeyholeMerge(std::vector<Pt2> Outer, std::vector<std::vector<Pt2>> Holes)
{
	if (SignedArea2(Outer) < 0.0)
	{
		std::reverse(Outer.begin(), Outer.end());
	}
	for (std::vector<Pt2>& H : Holes)
	{
		if (SignedArea2(H) > 0.0)
		{
			std::reverse(H.begin(), H.end());
		}
	}
	auto MaxX = [](const std::vector<Pt2>& H) {
		size_t B = 0;
		for (size_t i = 1; i < H.size(); ++i)
		{
			if (H[i].X > H[B].X)
			{
				B = i;
			}
		}
		return B;
	};
	std::sort(Holes.begin(), Holes.end(), [&MaxX](const std::vector<Pt2>& A, const std::vector<Pt2>& B) { return A[MaxX(A)].X > B[MaxX(B)].X; });
	for (const std::vector<Pt2>& H : Holes)
	{
		const size_t Hm = MaxX(H);
		const Pt2 M = H[Hm];
		double BestX = 1e30;
		size_t Pick = 0;
		for (size_t i = 0; i < Outer.size(); ++i)
		{
			const Pt2& A = Outer[i];
			const Pt2& B = Outer[(i + 1) % Outer.size()];
			if ((A.Y > M.Y) == (B.Y > M.Y))
			{
				continue;
			}
			const double X = A.X + (M.Y - A.Y) * (B.X - A.X) / (B.Y - A.Y);
			if (X >= M.X && X < BestX)
			{
				BestX = X;
				Pick = A.X > B.X ? i : (i + 1) % Outer.size();
			}
		}
		// A reflex vertex inside the triangle (M, hit, candidate) would block the bridge: take the one
		// closest in angle to the ray instead.
		const Pt2 Hit{BestX, M.Y};
		const Pt2 Cand = Outer[Pick];
		double BestAngle = 1e30;
		for (size_t i = 0; i < Outer.size(); ++i)
		{
			const Pt2& Prev = Outer[(i + Outer.size() - 1) % Outer.size()];
			const Pt2& Cur = Outer[i];
			const Pt2& Next = Outer[(i + 1) % Outer.size()];
			if (i == Pick || Cross2(Prev, Cur, Next) > 0.0 || Cur.X < M.X)
			{
				continue;
			}
			const bool In = Cross2(M, Hit, Cur) * Cross2(M, Hit, Cand) >= 0.0 && Cross2(Hit, Cand, Cur) * Cross2(Hit, Cand, M) >= 0.0 && Cross2(Cand, M, Cur) * Cross2(Cand, M, Hit) >= 0.0;
			if (!In)
			{
				continue;
			}
			const double Ang = std::fabs(std::atan2(Cur.Y - M.Y, Cur.X - M.X));
			if (Ang < BestAngle)
			{
				BestAngle = Ang;
				Pick = i;
			}
		}
		std::vector<Pt2> Merged(Outer.begin(), Outer.begin() + static_cast<long>(Pick) + 1);
		for (size_t k = 0; k <= H.size(); ++k)
		{
			Merged.push_back(H[(Hm + k) % H.size()]);
		}
		Merged.push_back(Outer[Pick]);
		Merged.insert(Merged.end(), Outer.begin() + static_cast<long>(Pick) + 1, Outer.end());
		Outer = Merged;
	}
	return Outer;
}

// A flat plate with holes: Outer and Holes in the frame's (U, V) plane, from W0 to W1 along W.
void PrismWithHoles(Mesh& M, std::vector<Pt2> Outer, std::vector<std::vector<Pt2>> Holes, double W0, double W1, const Frame3& F)
{
	if (SignedArea2(Outer) < 0.0)
	{
		std::reverse(Outer.begin(), Outer.end());
	}
	for (std::vector<Pt2>& H : Holes)
	{
		if (SignedArea2(H) > 0.0)
		{
			std::reverse(H.begin(), H.end());
		}
	}
	const double Hand = F.U.Cross(F.V).Dot(F.W) >= 0.0 ? 1.0 : -1.0;
	const std::vector<Pt2> Poly = KeyholeMerge(Outer, Holes);
	const std::vector<int> Tris = EarClip(Poly);
	for (size_t t = 0; t + 2 < Tris.size(); t += 3)
	{
		const Pt2 &A = Poly[static_cast<size_t>(Tris[t])], &B = Poly[static_cast<size_t>(Tris[t + 1])], &C = Poly[static_cast<size_t>(Tris[t + 2])];
		FlatTri(M, F.At(A.X, A.Y, W1), F.At(B.X, B.Y, W1), F.At(C.X, C.Y, W1), F.W * Hand);
		FlatTri(M, F.At(A.X, A.Y, W0), F.At(B.X, B.Y, W0), F.At(C.X, C.Y, W0), F.W * -Hand);
	}
	std::vector<std::vector<Pt2>> Loops = Holes;
	Loops.push_back(Outer);
	for (const std::vector<Pt2>& L : Loops)
	{
		for (size_t i = 0; i < L.size(); ++i)
		{
			const Pt2& P0 = L[i];
			const Pt2& P1 = L[(i + 1) % L.size()];
			const Vec3 Out = (F.U * (P1.Y - P0.Y) - F.V * (P1.X - P0.X)) * Hand;
			const Vec3 A = F.At(P0.X, P0.Y, W0), B = F.At(P1.X, P1.Y, W0), C = F.At(P1.X, P1.Y, W1), D = F.At(P0.X, P0.Y, W1);
			FlatTri(M, A, B, C, Out);
			FlatTri(M, A, C, D, Out);
		}
	}
}

std::vector<Pt2> CircleHole(double Cx, double Cy, double R, int N)
{
	return CircleLoop(Cx, Cy, R, N).P;
}

// A pointed (gothic) arch window: straight sides from Bottom up to Spring, two arcs meeting at Apex.
std::vector<Pt2> GothicHole(double Y0, double Y1, double Bottom, double Spring, double Fillet)
{
	std::vector<Pt2> H;
	const double W = Y1 - Y0, R = 0.8 * W;
	const double Mid = 0.5 * (Y0 + Y1);
	// Arc centred at (Y0 + R, Spring) for the right-hand side, and (Y1 - R, Spring) for the left.
	const double Apex = Spring + std::sqrt(R * R - (Mid - Y0 - R) * (Mid - Y0 - R));
	H.push_back(Pt2{Y0 + Fillet, Bottom});
	H.push_back(Pt2{Y1 - Fillet, Bottom});
	H.push_back(Pt2{Y1, Bottom + Fillet});
	const double A1 = std::atan2(Apex - Spring, Mid - (Y1 - R));
	for (int K = 0; K <= 8; ++K)
	{
		const double A = A1 * K / 8.0;
		H.push_back(Pt2{Y1 - R + R * std::cos(A), Spring + R * std::sin(A)});
	}
	for (int K = 7; K >= 0; --K)
	{
		const double A = A1 * K / 8.0;
		H.push_back(Pt2{Y0 + R - R * std::cos(A), Spring + R * std::sin(A)});
	}
	H.push_back(Pt2{Y0, Bottom + Fillet});
	std::vector<Pt2> Clean;
	for (const Pt2& Q : H)
	{
		if (Clean.empty() || std::fabs(Clean.back().X - Q.X) + std::fabs(Clean.back().Y - Q.Y) > 1e-6)
		{
			Clean.push_back(Q);
		}
	}
	return Clean;
}

// The outline of a cast-iron side frame in (y, z): feet, an arched underside, waisted ends, an arched
// top and, when ArmY < Y0, a bracket arm reaching forward to a shaft at (ArmY, ArmZ).
std::vector<Pt2> FrameOutline(double Y0, double Y1, double Top, double ArmY, double ArmZ)
{
	std::vector<Pt2> P;
	auto Arc = [&P](double Cy, double Cz, double R, double A0, double A1, int N) {
		for (int K = 0; K <= N; ++K)
		{
			const double A = (A0 + (A1 - A0) * K / N) * kGearDeg;
			P.push_back(Pt2{Cy + R * std::cos(A), Cz + R * std::sin(A)});
		}
	};
	const double Mid = 0.5 * (Y0 + Y1), Span = Y1 - Y0;
	P.push_back(Pt2{Y0 - 4.0, 0.0});
	P.push_back(Pt2{Y0 + 9.0, 0.0});
	for (int K = 1; K < 10; ++K)
	{
		const double T = static_cast<double>(K) / 10.0;
		const double Y = Y0 + 9.0 + (Span - 18.0) * T;
		P.push_back(Pt2{Y, 7.0 * std::sin(kPi * T)});
	}
	P.push_back(Pt2{Y1 - 9.0, 0.0});
	P.push_back(Pt2{Y1 + 4.0, 0.0});
	P.push_back(Pt2{Y1 + 4.0, 5.0});
	Arc(Y1 + 4.0 + 6.0, 12.0, 7.0, 210.0, 150.0, 4);
	P.push_back(Pt2{Y1 - 1.0, Top - 16.0});
	Arc(Y1 - 5.0, Top - 6.0, 6.0, -30.0, 90.0, 6);
	for (int K = 1; K < 8; ++K)
	{
		const double T = static_cast<double>(K) / 8.0;
		P.push_back(Pt2{Y1 - 5.0 - (Y1 - Y0 - 9.0) * T, Top + 3.0 * std::sin(kPi * T)});
	}
	Arc(Y0 + 4.0, Top - 5.0, 5.0, 90.0, 200.0, 5);
	if (ArmY < Y0)
	{
		const double R = 3.8;
		P.push_back(Pt2{Y0 - 0.5, ArmZ + R});
		Arc(ArmY, ArmZ, R, 90.0, 270.0, 8);
		P.push_back(Pt2{Y0 - 0.5, ArmZ - R});
	}
	Arc(Y0 - 4.0 - 5.0, 14.0, 6.0, 30.0, -30.0, 4);
	P.push_back(Pt2{Y0 - 4.0, 5.0});
	(void)Mid;
	return P;
}

// The machine's dimensions in the desk frame, worked out from the keyboard's position.
struct MechLayout
{
	double Hw = 0.0;
	double YF = 0.0;                 // front of the mechanism (the accumulator window bar)
	double KbMinY = 0.0, KbMinX = 0.0, KbMaxX = 0.0;
	double Pitch = 14.0, Step = 0.0;
	double X0 = 0.0;
	double ZBase = 6.0, ZDrum = 38.0, ZSq = 0.0, ZAcc = 0.0, ZCnt = 59.5, ZSet = 70.0, ZShift = 15.0, ZRail = 21.0;
	double YAcc = 0.0, YTr = 0.0, YCnt = 0.0, YMain = 0.0, YBack = 0.0, YAp = 0.0;
	double LeftFrame = 0.0, RightFrame = 0.0;   // inner faces (toward the drums)
	double XRev = 0.0, XDrive = 0.0, XFly = 0.0, XCrank = 0.0, XOuter = 0.0, XSp = 0.0;
	double XL = 0.0, XR = 0.0;       // carriage ends, carriage frame
	double YMem = 0.0;
	double PlinthX0 = 0.0, PlinthX1 = 0.0, PlinthY0 = 0.0, PlinthY1 = 0.0, PlinthH = 18.0;
	StepDrumParams Drum;
	StepDrumLayout DrumL;
	double ColumnX(int J) const { return X0 - J * Pitch; }
};

constexpr double kSqHalf = 1.3;
constexpr double kAccRadius = 6.0, kCntRadius = 5.8, kSetRadius = 6.0, kWheelWidth = 5.5;
constexpr int kPinionTeeth = 12, kShaftGearTeeth = 18, kWheelGearTeeth = 15;
constexpr double kDrumModule = 0.6, kWheelGearModule = 0.5;
constexpr double kTransferAngle = 29.0;   // the accumulator axis sits below and in front of the transfer shaft

Frame3 AxisFrame(const Vec3& Origin, const Vec3& W, const Vec3& V)
{
	Frame3 F;
	F.O = Origin;
	F.W = W;
	F.V = V;
	F.U = V.Cross(W);
	return F;
}

void AddBox(Mesh& M, double X0, double X1, double Y0, double Y1, double Z0, double Z1)
{
	Box(M, X0, Y0, Z0, X1, Y1, Z1);
}

void AddRod(Mesh& M, const Vec3& A, const Vec3& B, double R)
{
	AddInFrame(M, BuildShaft(R, (B - A).Length(), false), FrameAlong(A, B - A));
}

void AddBevel(Mesh& M, const Vec3& Apex, const Vec3& Axis, int Teeth, double Module, double Face, double Bore, double HubR, double HubL)
{
	BevelGearParams B;
	B.Teeth = Teeth;
	B.Module = Module;
	B.FaceWidth = Face;
	B.BoreRadius = Bore;
	B.HubRadius = HubR;
	B.HubLength = HubL;
	AddInFrame(M, BuildBevelGear(B), FrameAlong(Apex, Axis));
}

void AddSpur(Mesh& M, const Vec3& Center, const Vec3& Axis, int Teeth, double Module, double Face, double Bore)
{
	SpurGearParams G;
	G.Teeth = Teeth;
	G.Module = Module;
	G.FaceWidth = Face;
	G.BoreRadius = Bore;
	AddInFrame(M, BuildSpurGear(G), FrameAlong(Center, Axis));
}

void AddLathe(Mesh& M, const std::vector<Pt2>& Profile, const Vec3& Origin, const Vec3& Axis, int Segments = 24)
{
	Mesh L;
	Lathe(L, Profile, Segments, false, SharpAll(Profile.size()));
	AddInFrame(M, L, FrameAlong(Origin, Axis));
}

// A boss round a bearing on a frame face (cast with the frame) and its bright steel cap.
void AddBearingBoss(Mesh& Iron, Mesh& Steel, const Vec3& Center, const Vec3& Out, double R)
{
	AddLathe(Iron, {{0.0, -0.5}, {R, -0.5}, {R, 2.2}, {R - 0.8, 3.0}, {0.0, 3.0}}, Center, Out, 28);
	AddLathe(Steel, {{0.0, 2.9}, {0.62 * R, 2.9}, {0.62 * R, 3.8}, {0.45 * R, 4.6}, {0.0, 4.8}}, Center, Out, 24);
}

// A slotted register window bar (black), in the x-z plane facing -y: one slot per centre.
void AddWindowBar(Mesh& M, double X0, double X1, double YFront, double Z0, double Z1, const std::vector<double>& Slots, double Zc, double Sw, double Sh)
{
	std::vector<HoleRect> Holes;
	for (double X : Slots)
	{
		Holes.push_back({X - 0.5 * Sw, X + 0.5 * Sw, Zc - 0.5 * Sh, Zc + 0.5 * Sh});
	}
	Mesh Bar;
	PerforatedSlab(Bar, X0, X1, Z0, Z1, Holes, 0.0, 1.2);
	Frame3 F;
	F.O = Vec3(0, YFront + 1.2, 0);
	F.U = Vec3(1, 0, 0);
	F.V = Vec3(0, 0, 1);
	F.W = Vec3(0, -1, 0);
	ApplyFrame(Bar, F);
	M.Append(Bar);
}

ReadingWindow WindowAt(const Vec3& Center, double W, double H)
{
	ReadingWindow R;
	R.Center = Center;
	R.Normal = Vec3(0, -1, 0);
	R.Up = Vec3(0, 0, 1);
	R.WidthMm = W;
	R.HeightMm = H;
	return R;
}

NumeralWheel RegisterWheel(double Radius, int GearTeeth, double GearGap, int StarTeeth, bool Carry)
{
	NumeralWheelParams P;
	P.Radius = Radius;
	P.Width = kWheelWidth;
	P.GlyphsAlongAxis = false;
	P.GlyphHeight = 3.0;
	P.GlyphWidth = 2.2;
	P.StrokeMm = 0.42;
	P.EmbossMm = 0.25;
	P.BoreRadius = 0.85;
	P.GearTeeth = GearTeeth;
	P.GearModule = kWheelGearModule;
	P.GearGap = GearGap;
	P.StarTeeth = StarTeeth;
	P.StarModule = 0.6;
	P.CarryTooth = Carry;
	return BuildNumeralWheel(P);
}

void SetWheelModel(MechPart& Part, const NumeralWheel& W, MechMaterial Rim, MechMaterial Glyphs)
{
	Part.Of(Rim).Append(W.Wheel);
	if (!W.Gears.Positions.empty())
	{
		Part.Of(MechMaterial::Brass).Append(W.Gears);
	}
	Part.Of(Glyphs).Append(W.Glyphs);
}

// A decimal-point marker: a small brass wedge riding on a window bar's top edge. Local X runs along
// the bar toward the higher digits, Y up, Z through the bar (+y).
void BuildMarker(MechPart& Part)
{
	Mesh& B = Part.Of(MechMaterial::Brass);
	Prism(B, {{-1.3, 2.4}, {0.0, 0.0}, {1.3, 2.4}}, -0.9, 0.9, Frame3());
	Prism(B, RectPoly(-1.6, 2.4, 1.6, 3.2), -0.9, 2.3, Frame3());
	Prism(B, RectPoly(-1.5, 0.8, 1.5, 3.1), 1.3, 2.2, Frame3());
}

MechLayout Plan(const mkb::Layout& Layout, const CaseParams& Case)
{
	MechLayout L;
	const BoardParts Board = BuildBoard(Layout, Case);
	const double Ct = std::cos(Board.RestTiltRad), St = std::sin(Board.RestTiltRad);
	double MaxY = -1e30, MinY = 1e30, MinX = 1e30, MaxX = -1e30;
	for (const Mesh* M : {&Board.CaseShell, &Board.CaseBadge, &Board.Feet})
	{
		for (const Vec3& P : M->Positions)
		{
			const double Y = P.Y * Ct - P.Z * St;
			MaxY = std::max(MaxY, Y);
			MinY = std::min(MinY, Y);
			MinX = std::min(MinX, P.X);
			MaxX = std::max(MaxX, P.X);
		}
	}
	L.KbMinY = MinY;
	L.KbMinX = MinX;
	L.KbMaxX = MaxX;
	L.Hw = 0.5 * Layout.Width * mkb::kUnitMm + Case.MarginMm;
	L.Step = L.Drum.StepMm;
	L.DrumL = StepDrumGeometry(L.Drum);
	L.X0 = 144.0;
	L.LeftFrame = L.ColumnX(19) - 8.0;
	L.RightFrame = L.X0 + 9.0;
	L.ZSq = L.ZDrum + 0.5 * kDrumModule * (L.Drum.ToothPositions + kPinionTeeth);
	const double Cd = 0.5 * kWheelGearModule * (kShaftGearTeeth + kWheelGearTeeth);
	L.ZAcc = L.ZSq - Cd * std::cos(kTransferAngle * kGearDeg);
	L.YF = MaxY + 6.0;
	L.YAcc = L.YF + 8.0;
	L.YTr = L.YAcc + Cd * std::sin(kTransferAngle * kGearDeg);
	L.YCnt = L.YAcc + 2.5;
	L.YMain = L.YTr + L.DrumL.BodyEndZ + 3.5 + 6.5;
	L.YBack = L.YMain + 12.0;
	L.YAp = L.YAcc + 16.0;
	L.XSp = L.X0 - 24.0;
	L.XRev = L.RightFrame + 5.0 + 22.0;
	L.XOuter = L.XRev + 20.0;
	L.XDrive = L.XOuter + 12.0;
	L.XFly = L.XDrive + 10.0;
	L.XCrank = L.XFly + 6.0;
	L.XR = 7.5;
	L.XL = -(39.0 * L.Pitch) - 12.6;
	L.YMem = L.YAcc + 31.0;
	L.PlinthX0 = std::min(std::min(L.KbMinX, -L.Hw), L.X0 + L.XL) - 14.0;   // under the carriage at shift 0 too
	L.PlinthX1 = L.XCrank + 30.0;
	L.PlinthY0 = L.KbMinY - 18.0;
	L.PlinthY1 = L.YMain + 64.0;
	return L;
}

void BuildPlinth(MechanismParts& Mx, const MechLayout& L)
{
	Mesh& E = Mx.Plinth.Of(MechMaterial::BlackEnamel);
	const double Cx = 0.5 * (L.PlinthX0 + L.PlinthX1), Cy = 0.5 * (L.PlinthY0 + L.PlinthY1);
	const double Hx = 0.5 * (L.PlinthX1 - L.PlinthX0), Hy = 0.5 * (L.PlinthY1 - L.PlinthY0);
	const double H = L.PlinthH, C = 2.0;
	const Splits None;
	const std::vector<Vec3> Bot = AtZ(RoundedRect(Cx, Cy, Hx, Hy, 6.0, 6, None), -H);
	const std::vector<Vec3> Mid = AtZ(RoundedRect(Cx, Cy, Hx, Hy, 6.0, 6, None), -C);
	const std::vector<Vec3> Top = AtZ(RoundedRect(Cx, Cy, Hx - C, Hy - C, 4.0, 6, None), 0.0);
	Loft(E, {Bot, Mid});
	Loft(E, {Mid, Top});
	FanCap(E, Top, Vec3(Cx, Cy, 0.0), Vec3(0, 0, 1));
	FanCap(E, Bot, Vec3(Cx, Cy, -H), Vec3(0, 0, -1));
	// The maker's plate on the front edge: a brass plate with a bevelled border and four screws.
	Mesh& Br = Mx.Plinth.Of(MechMaterial::Brass);
	const double Py = L.PlinthY0, Pz = -0.5 * H - 0.5, Pw = 170.0, Ph = 10.0;
	Mesh Plate;
	const Loop2 Outer = RoundedRect(0, 0, 0.5 * Pw, 0.5 * Ph, 1.5, 4, None);
	const Loop2 Inner = RoundedRect(0, 0, 0.5 * Pw - 0.8, 0.5 * Ph - 0.8, 1.0, 4, None);
	Loft(Plate, {AtZ(Outer, -0.5), AtZ(Outer, 0.6), AtZ(Inner, 1.1)});
	FanCap(Plate, AtZ(Inner, 1.1), Vec3(0, 0, 1.1), Vec3(0, 0, 1));
	FanCap(Plate, AtZ(Outer, -0.5), Vec3(0, 0, -0.5), Vec3(0, 0, -1));
	for (double Sx = -1.0; Sx <= 1.0; Sx += 2.0)
	{
		Mesh Screw;
		Lathe(Screw, {{0.0, 0.9}, {1.0, 0.9}, {1.0, 1.3}, {0.6, 1.7}, {0.0, 1.8}}, 12, false, SharpAll(5));
		Screw.Translate(Vec3(Sx * (0.5 * Pw - 3.0), 0.0, 0.0));
		Plate.Append(Screw);
	}
	Frame3 F;
	F.O = Vec3(0.5 * (L.KbMinX + L.KbMaxX), Py, Pz);
	F.U = Vec3(1, 0, 0);
	F.V = Vec3(0, 0, 1);
	F.W = Vec3(0, -1, 0);
	ApplyFrame(Plate, F);
	Br.Append(Plate);
	Mx.MakersPlate.Center = Vec3(F.O.X, Py - 1.1, Pz);
	Mx.MakersPlate.Normal = Vec3(0, -1, 0);
	Mx.MakersPlate.Up = Vec3(0, 0, 1);
	Mx.MakersPlate.WidthMm = Pw - 12.0;
	Mx.MakersPlate.HeightMm = Ph - 2.4;
	Mx.PlinthHeightMm = H;
}

void BuildFrames(MechanismParts& Mx, const MechLayout& L)
{
	Mesh& Fe = Mx.Frame.Of(MechMaterial::CastIron);
	Mesh& St = Mx.Frame.Of(MechMaterial::Steel);
	Mesh& Be = Mx.Frame.Of(MechMaterial::BlackEnamel);
	const Vec3 Xw(1, 0, 0), Yw(0, 1, 0), Zw(0, 0, 1);
	// Base plate under the mechanism.
	const double BaseX0 = L.LeftFrame - 14.0, BaseX1 = L.XOuter + 6.0;
	SolidFromLoop(Fe, RoundedRect(0.5 * (BaseX0 + BaseX1), 0.5 * (L.YAcc + 10.0 + L.YMain + 52.0), 0.5 * (BaseX1 - BaseX0), 0.5 * (L.YMain + 52.0 - L.YAcc - 10.0), 5.0, 6,
		Splits()), 0.0, L.ZBase);
	// The two side frames of the drum bank, with a bracket arm out to the transfer shaft.
	const double Y0 = L.YAcc + 10.5, Y1 = L.YBack, Top = 84.0;
	const double Mid = 0.5 * (Y0 + L.YMain);
	std::vector<std::vector<Pt2>> Holes = {GothicHole(Y0 + 7.0, L.YMain - 10.0, 13.0, 30.0, 2.0), CircleHole(Y0 + 9.0, 64.0, 4.2, 20), CircleHole(Mid + 4.0, 70.0, 3.2, 16),
		CircleHole(L.YMain + 6.0, 18.0, 3.0, 16)};
	const std::vector<Pt2> Outline = FrameOutline(Y0, Y1, Top, L.YTr, L.ZSq);
	for (int Side = 0; Side < 2; ++Side)
	{
		const double X = Side == 0 ? L.LeftFrame - 5.0 : L.RightFrame;
		PrismWithHoles(Fe, Outline, Holes, X, X + 5.0, AxisFrame(Vec3(), Xw, Zw));
		const double Face = Side == 0 ? X : X + 5.0;
		const Vec3 Out = Side == 0 ? -Xw : Xw;
		for (const Vec3& B : {Vec3(Face, L.YMain, L.ZDrum), Vec3(Face, L.YMain, L.ZSet), Vec3(Face, L.YTr, L.ZSq), Vec3(Face, L.YMain + 7.5, 47.5)})
		{
			AddBearingBoss(Fe, St, B, Out, B.Y == L.YTr ? 3.6 : 6.0);
		}
	}
	// Bars and tie rods between the frames.
	const double Xa = L.LeftFrame - 0.1, Xb = L.RightFrame + 0.1;
	const double DrumFront = L.YMain - L.DrumL.BodyEndZ - 3.0;
	AddBox(Be, Xa, Xb, DrumFront - 3.0, DrumFront, 34.0, 42.0);
	AddBox(Be, Xa, Xb, L.YMain + 6.0, L.YMain + 9.0, 44.5, 50.5);
	for (const Vec3& T : {Vec3(0, Y1 - 2.0, 10.0), Vec3(0, Y1 - 4.0, Top - 8.0), Vec3(0, Y0 + 3.0, 10.0), Vec3(0, Y0 + 5.0, Top - 9.0)})
	{
		AddRod(St, Vec3(L.LeftFrame - 8.0, T.Y, T.Z), Vec3(L.RightFrame + 8.0, T.Y, T.Z), 1.6);
		for (double X : {L.LeftFrame - 9.0, L.RightFrame + 6.0})
		{
			AddLathe(St, {{0.0, 0.0}, {2.8, 0.0}, {2.8, 2.2}, {2.2, 3.0}, {0.0, 3.0}}, Vec3(X, T.Y, T.Z), Xw, 6);
		}
	}
	// Transfer shaft and the setting wheels' shaft.
	AddRod(St, Vec3(L.LeftFrame - 3.0, L.YTr, L.ZSq), Vec3(L.RightFrame + 3.0, L.YTr, L.ZSq), 1.0);
	AddRod(St, Vec3(L.LeftFrame - 3.0, L.YMain, L.ZSet), Vec3(L.RightFrame + 3.0, L.YMain, L.ZSet), 0.8);
	// Main shaft pedestals.
	for (int J : {5, 14})
	{
		Mesh B = BuildBearingBlock(2.05, L.ZDrum - L.ZBase, 3.0);
		B.Translate(Vec3(L.ColumnX(J) + 1.5, L.YMain, L.ZBase));
		Be.Append(B);
	}
	// The setting register's window bar.
	std::vector<double> Slots;
	for (int J = 0; J < Mx.SettingDigits; ++J)
	{
		Slots.push_back(L.ColumnX(J));
		Mx.SettingWindows.push_back(WindowAt(Vec3(L.ColumnX(J), L.YMain - kSetRadius - 1.6, L.ZSet), 3.4, 3.9));
	}
	AddWindowBar(Be, L.LeftFrame, L.RightFrame, L.YMain - kSetRadius - 1.6, L.ZSet - 6.8, L.ZSet + 6.8, Slots, L.ZSet, 3.4, 3.9);
	// Carriage rails on posts.
	for (double Y : {L.YAcc - 3.0, L.YAcc + 6.0})
	{
		AddRod(St, Vec3(L.PlinthX0 + 6.0, Y, L.ZRail), Vec3(L.PlinthX1 - 6.0, Y, L.ZRail), 2.0);
		for (double X = L.PlinthX0 + 20.0; X < L.PlinthX1 - 10.0; X += 96.0)
		{
			const double Xp = std::fabs(X - L.XSp) < 6.0 ? X + 12.0 : X;
			AddLathe(Fe, {{0.0, 0.0}, {3.6, 0.0}, {3.0, 1.5}, {2.2, 2.5}, {2.2, L.ZRail - 1.0}, {0.0, L.ZRail}}, Vec3(Xp, Y, 0.0), Zw, 16);
		}
	}
	// The counting finger's bracket, reaching in from the right frame behind the carriage.
	AddBox(Be, L.X0 + 2.9, L.RightFrame + 0.1, L.YAcc + 10.6, L.YAcc + 12.0, L.ZCnt + 4.3, L.ZCnt + 5.7);
	AddBox(Be, L.X0 + 2.9, L.X0 + 4.2, L.YCnt + 4.3, L.YAcc + 11.4, L.ZCnt + 4.4, L.ZCnt + 5.6);
	// Drive housing: the outer drive frame, its bosses and the tie rods to the right frame.
	const double Do0 = L.YAcc + 22.0, Do1 = L.YMain + 48.0;
	const std::vector<Pt2> Drive = FrameOutline(Do0, Do1, 78.0, Do0 + 1.0, 0.0);
	PrismWithHoles(Fe, Drive, {GothicHole(Do0 + 8.0, L.YMain - 10.0, 12.0, 26.0, 2.0), GothicHole(L.YMain + 10.0, L.YMain + 30.0, 12.0, 52.0, 2.0)}, L.XOuter, L.XOuter + 5.0,
		AxisFrame(Vec3(), Xw, Zw));
	AddBearingBoss(Fe, St, Vec3(L.XOuter + 5.0, L.YMain, L.ZDrum), Xw, 6.5);
	AddBearingBoss(Fe, St, Vec3(L.XOuter + 5.0, L.YMain + 40.0, L.ZDrum), Xw, 6.5);
	for (const Vec3& T : {Vec3(0, Do0 + 4.0, 10.0), Vec3(0, Do0 + 6.0, 70.0), Vec3(0, Do1 - 2.0, 10.0), Vec3(0, Do1 - 6.0, 66.0)})
	{
		AddRod(St, Vec3(L.RightFrame + 4.9, T.Y, T.Z), Vec3(L.XOuter + 0.1, T.Y, T.Z), 1.6);
	}
	// Reverser idler stud bracket and the reversing lever's pivot bracket.
	AddBox(Be, L.XRev - 2.5, L.XRev + 2.5, L.YMain + 12.0, L.YMain + 15.0, L.ZBase - 0.1, L.ZDrum + 3.0);
	AddBox(Be, L.RightFrame + 4.9, L.XOuter + 0.1, L.YMain + 3.0, L.YMain + 6.0, 60.0, 63.0);
	// The optional glass cover over the drum bank.
	AddBox(Mx.GlassCover.Of(MechMaterial::Glass), L.LeftFrame, L.RightFrame, Y0, L.YBack, Top + 3.5, Top + 5.0);
}

void BuildColumns(MechanismParts& Mx, const MechLayout& L)
{
	const StepDrumLayout& D = L.DrumL;
	const Vec3 Xw(1, 0, 0), Yw(0, 1, 0), Zw(0, 0, 1);
	// Drum: its bevel at the back meets the main shaft; the body runs forward.
	Mx.Drums.Model.Of(MechMaterial::Brass) = BuildStepDrum(L.Drum);
	{
		Mesh Axle = BuildShaft(L.Drum.AxleRadius, D.BodyEndZ + 3.0 - 4.0, false);
		Axle.Translate(Vec3(0, 0, 4.0));
		Mx.Drums.Model.Of(MechMaterial::Steel) = Axle;
	}
	// Square shaft: the front bevel meets the transfer gear; the square carries the pinion.
	const double SqLen = L.YMain + 8.5 - L.YTr;
	{
		Mesh& S = Mx.SquareShafts.Model.Of(MechMaterial::Steel);
		Mesh Sq = BuildShaft(kSqHalf, SqLen - 9.0, true);
		Sq.Translate(Vec3(0, 0, 5.5));
		S.Append(Sq);
		Mesh Rear = BuildShaft(1.2, 5.0, false);
		Rear.Translate(Vec3(0, 0, SqLen - 4.0));
		S.Append(Rear);
		Mesh Front = BuildShaft(1.2, 3.0, false);
		Front.Translate(Vec3(0, 0, 4.5));
		S.Append(Front);
		AddBevel(Mx.SquareShafts.Model.Of(MechMaterial::Brass), Vec3(0, 0, 0), Zw, 16, kWheelGearModule, 1.4, 1.0, 2.0, 1.2);
	}
	{
		SpurGearParams G;
		G.Teeth = kPinionTeeth;
		G.Module = kDrumModule;
		G.FaceWidth = L.Drum.PinionFaceWidth;
		G.BoreRadius = kSqHalf + 0.05;
		G.SquareBore = true;
		G.HubRadius = 2.4;
		G.HubLength = 1.6;
		Mx.Pinions.Model.Of(MechMaterial::Brass) = BuildSpurGear(G);
	}
	// Transfer gear: the bevel (right of the column) and the 18-tooth spur (left), on a sleeve.
	{
		Mesh& T = Mx.TransferGears.Model.Of(MechMaterial::Brass);
		AddBevel(T, Vec3(0, 0, 0), Zw, 16, kWheelGearModule, 1.4, 1.05, 2.0, 0.6);
		AddSpur(T, Vec3(0, 0, -5.5), Zw, kShaftGearTeeth, kWheelGearModule, 1.6, 1.05);
		Lathe(T, {{1.05, -4.8}, {1.6, -4.8}, {1.6, 2.95}, {1.05, 2.95}}, 16, true, SharpAll(4));
	}
	SetWheelModel(Mx.SettingWheels.Model, RegisterWheel(kSetRadius, 12, 0.6, 0, false), MechMaterial::NumeralBlack, MechMaterial::Paint);
	BuildMarker(Mx.SettingMarker.Model);
	Mx.SettingMarker.Poses.push_back(PoseOf(Vec3(L.X0 + 0.5 * L.Pitch, L.YMain - kSetRadius - 1.6 + 0.6, L.ZSet + 6.8), Yw, Zw).RotatedZ(0.0));
	{
		Pose& M = Mx.SettingMarker.Poses.back();
		M.X = -Xw;
		M.Y = Zw;
		M.Z = M.X.Cross(M.Y);
	}
	for (int J = 0; J < Mx.SettingDigits; ++J)
	{
		const double X = L.ColumnX(J);
		Mx.Drums.Poses.push_back(PoseOf(Vec3(X, L.YMain, L.ZDrum), -Yw, Zw));
		Mx.SquareShafts.Poses.push_back(PoseOf(Vec3(X, L.YTr, L.ZSq), Yw, Zw));
		Mx.Pinions.Poses.push_back(PoseOf(Vec3(X, L.YMain - D.PinionCenterZ(0), L.ZSq), -Yw, Zw));
		Mx.TransferGears.Poses.push_back(PoseOf(Vec3(X, L.YTr, L.ZSq), Xw, Zw));
		Mx.SettingWheels.Poses.push_back(PoseOf(Vec3(X, L.YMain, L.ZSet), Xw, -Yw));
		GearPair Pp;
		Pp.Name = "drum/pinion";
		Pp.PointA = Vec3(X, L.YMain, L.ZDrum);
		Pp.AxisA = Yw;
		Pp.PointB = Vec3(X, L.YMain, L.ZSq);
		Pp.AxisB = Yw;
		Pp.Module = kDrumModule;
		Pp.TeethA = L.Drum.ToothPositions;
		Pp.TeethB = kPinionTeeth;
		Mx.GearPairs.push_back(Pp);
		GearPair Bp;
		Bp.Name = "main/drum";
		Bp.Type = GearPair::Bevel;
		Bp.PointA = Bp.PointB = Vec3(X, L.YMain, L.ZDrum);
		Bp.AxisA = -Xw;
		Bp.AxisB = -Yw;
		Bp.Module = kDrumModule;
		Bp.TeethA = Bp.TeethB = L.Drum.BevelTeeth;
		Mx.GearPairs.push_back(Bp);
		GearPair Tb;
		Tb.Name = "shaft/transfer";
		Tb.Type = GearPair::Bevel;
		Tb.PointA = Tb.PointB = Vec3(X, L.YTr, L.ZSq);
		Tb.AxisA = Yw;
		Tb.AxisB = Xw;
		Tb.Module = kWheelGearModule;
		Tb.TeethA = Tb.TeethB = 16;
		Mx.GearPairs.push_back(Tb);
		GearPair Wp;
		Wp.Name = "transfer/wheel";
		Wp.PointA = Vec3(X - 5.5, L.YTr, L.ZSq);
		Wp.AxisA = Xw;
		Wp.PointB = Vec3(-J * L.Pitch - 5.5, L.YAcc, L.ZAcc);
		Wp.AxisB = Xw;
		Wp.OnCarriageB = true;
		Wp.Module = kWheelGearModule;
		Wp.TeethA = kShaftGearTeeth;
		Wp.TeethB = kWheelGearTeeth;
		Mx.GearPairs.push_back(Wp);
	}
}

void BuildCarriage(MechanismParts& Mx, const MechLayout& L)
{
	const Vec3 Xw(1, 0, 0), Yw(0, 1, 0), Zw(0, 0, 1);
	Mx.Carriage.Poses.push_back(PoseOf(Vec3(L.X0, 0, 0), Zw, Yw));
	Mesh& E = Mx.Carriage.Model.Of(MechMaterial::BlackEnamel);
	Mesh& St = Mx.Carriage.Model.Of(MechMaterial::Steel);
	Mesh& Br = Mx.Carriage.Model.Of(MechMaterial::Brass);
	const double Yf = L.YF + 0.4, Zb = L.ZRail + 2.0, ZTop = L.ZCnt + 6.7;
	const double AccBarY = L.YAcc - kAccRadius - 1.6, CntBarY = L.YCnt - kCntRadius - 1.6;
	// End plates, shaped round the fixed transfer gears.
	const std::vector<Pt2> EndPlate = {{Yf, Zb}, {L.YAcc + 8.0, Zb}, {L.YAcc + 8.0, L.ZSq - 6.0}, {L.YAcc - 1.8, L.ZSq - 6.0}, {L.YAcc - 1.8, L.ZSq + 5.6},
		{L.YAcc + 9.0, L.ZSq + 5.6}, {L.YAcc + 9.0, ZTop + 1.4}, {Yf, ZTop + 1.4}};
	for (double Xa : {L.XL, L.XR - 3.0})
	{
		Prism(E, EndPlate, Xa, Xa + 3.0, AxisFrame(Vec3(), Xw, Zw));
	}
	AddBox(E, L.XL + 1.0, L.XR - 1.0, Yf, L.YAcc + 8.0, Zb, Zb + 2.0);
	// Open on top (the wheels and their gears show): a front rail and a steel tie rod at the back.
	AddBox(E, L.XL + 1.0, L.XR - 1.0, Yf, Yf + 3.0, ZTop, ZTop + 1.4);
	AddRod(St, Vec3(L.XL + 0.5, L.YCnt + 1.5, ZTop + 0.7), Vec3(L.XR - 0.5, L.YCnt + 1.5, ZTop + 0.7), 0.8);
	// Window bars.
	std::vector<double> AccSlots, CntSlots;
	for (int I = 0; I < Mx.AccumulatorDigits; ++I)
	{
		AccSlots.push_back(-I * L.Pitch);
		Mx.AccumulatorWindows.push_back(WindowAt(Vec3(-I * L.Pitch, AccBarY, L.ZAcc), 3.4, 3.9));
		if (I < Mx.CounterDigits)
		{
			CntSlots.push_back(-I * L.Pitch);
			Mx.CounterWindows.push_back(WindowAt(Vec3(-I * L.Pitch, CntBarY, L.ZCnt), 3.4, 3.8));
		}
	}
	AddWindowBar(E, L.XL + 2.9, L.XR - 2.9, AccBarY, L.ZAcc - 6.8, L.ZAcc + 6.8, AccSlots, L.ZAcc, 3.4, 3.9);
	AddWindowBar(E, -(Mx.CounterDigits - 1) * L.Pitch - 9.0, L.XR - 2.9, CntBarY, L.ZCnt - 6.6, L.ZCnt + 6.6, CntSlots, L.ZCnt, 3.4, 3.8);
	// The counter row's left hanger, the shafts and the carry levers' pivot rod.
	const double CntLeft = -(Mx.CounterDigits - 1) * L.Pitch - 9.0;
	AddBox(E, CntLeft, CntLeft + 2.5, CntBarY + 1.0, L.YCnt + 2.5, L.ZCnt - 1.5, ZTop + 0.1);
	AddRod(St, Vec3(L.XL + 0.5, L.YAcc, L.ZAcc), Vec3(L.XR - 0.5, L.YAcc, L.ZAcc), 0.8);
	AddRod(St, Vec3(CntLeft + 0.5, L.YCnt, L.ZCnt), Vec3(L.XR - 0.5, L.YCnt, L.ZCnt), 0.8);
	AddRod(St, Vec3(L.XL + 0.5, L.YAcc + 1.0, Zb + 5.0), Vec3(L.XR - 0.5, L.YAcc + 1.0, Zb + 5.0), 0.6);
	// The shift rack under the bottom beam, teeth down.
	{
		const int Teeth = static_cast<int>((L.XR - L.XL - 8.0) / kPi);
		Mesh R = BuildRack(1.0, Teeth, 3.0, 1.8);
		Frame3 F;
		F.O = Vec3(L.XL + 4.0, L.YAcc + 2.0, L.ZShift + 7.0);
		F.U = Xw;
		F.V = -Zw;
		F.W = Yw;
		AddInFrame(Br, R, F);
	}
	// Wheels (axis along +x, read from the front) and the white carry levers under them.
	SetWheelModel(Mx.AccumulatorWheels.Model, RegisterWheel(kAccRadius, kWheelGearTeeth, 1.95, 10, true), MechMaterial::NumeralBlack, MechMaterial::Paint);
	SetWheelModel(Mx.CounterWheels.Model, RegisterWheel(kCntRadius, 0, 0.0, 10, false), MechMaterial::NumeralBlack, MechMaterial::Paint);
	{
		// Lever frame: pivot at the origin on the rod (axis +Z = world +x), X = world +y, Y up.
		Mesh& Lv = Mx.CarryLevers.Model.Of(MechMaterial::WhiteEnamel);
		const double Ax = L.YAcc - (L.YAcc + 1.0), Ay = L.ZAcc - (Zb + 5.0);
		const double Tx = Ax + 0.6, Ty = Ay - 3.25;
		Lathe(Lv, {{0.65, -0.5}, {1.5, -0.5}, {1.5, 0.5}, {0.65, 0.5}}, 16, true, SharpAll(4));
		Prism(Lv, {{-0.8, 0.9}, {0.8, 0.9}, {Tx + 0.7, Ty - 1.2}, {Tx + 0.5, Ty}, {Tx - 0.6, Ty}, {Tx - 0.8, Ty - 1.4}}, -0.5, 0.5, Frame3());
		Prism(Lv, {{0.6, -0.9}, {-0.6, -0.9}, {-2.8, -3.2}, {-1.8, -3.9}}, -0.5, 0.5, Frame3());
		Cylinder(Mx.CarryLevers.Model.Of(MechMaterial::Steel), Tx, Ty - 0.3, 0.3, 8, -3.6, 0.0);
	}
	BuildMarker(Mx.CarriageMarkers.Model);
	for (double By : {AccBarY, CntBarY})
	{
		Pose P;
		P.Origin = Vec3(0.5 * L.Pitch, By + 0.6, By == AccBarY ? L.ZAcc + 6.8 : L.ZCnt + 6.6);
		P.X = -Xw;
		P.Y = Zw;
		P.Z = P.X.Cross(P.Y);
		Mx.CarriageMarkers.Poses.push_back(P);
	}
	for (int I = 0; I < Mx.AccumulatorDigits; ++I)
	{
		const double X = -I * L.Pitch;
		Mx.AccumulatorWheels.Poses.push_back(PoseOf(Vec3(X, L.YAcc, L.ZAcc), Xw, -Yw));
		Mx.CarryLevers.Poses.push_back(PoseOf(Vec3(X - 0.5 * L.Pitch, L.YAcc + 1.0, Zb + 5.0), Xw, Zw));
		if (I < Mx.CounterDigits)
		{
			Mx.CounterWheels.Poses.push_back(PoseOf(Vec3(X, L.YCnt, L.ZCnt), Xw, -Yw));
		}
	}
	Mx.AccumulatorWheels.OnCarriage = Mx.CarryLevers.OnCarriage = Mx.CounterWheels.OnCarriage = Mx.CarriageMarkers.OnCarriage = true;
	// Clearing levers on the right end plate, turning about +x.
	{
		Mesh& Hub = Mx.ClearingLevers.Model.Of(MechMaterial::Steel);
		Cylinder(Hub, 0, 0, 1.8, 16, 0.0, 3.0);
		Prism(Hub, ThickenPolyline({{0.0, 0.0}, {0.0, 6.0}}, 0.9), 0.6, 2.4, Frame3());
		AddLathe(Mx.ClearingLevers.Model.Of(MechMaterial::Wood), {{0.0, 0.0}, {1.5, 0.0}, {1.8, 1.0}, {1.5, 2.4}, {0.0, 2.8}}, Vec3(0, 6.0, 1.5), Yw, 16);
		Mx.ClearingLevers.Poses.push_back(PoseOf(Vec3(L.XR + 0.1, L.YAcc - 3.5, L.ZAcc - 9.0), Xw, Zw));
		Mx.ClearingLevers.Poses.push_back(PoseOf(Vec3(L.XR + 0.1, L.YAcc - 3.5, L.ZCnt), Xw, Zw));
		Mx.ClearingLevers.OnCarriage = true;
	}
}

void BuildDrive(MechanismParts& Mx, const MechLayout& L)
{
	const Vec3 Xw(1, 0, 0), Yw(0, 1, 0), Zw(0, 0, 1);
	// Main shaft (local z along world +x, origin at x = 0) with a drum bevel left of each column.
	{
		Mesh& St = Mx.MainShaft.Model.Of(MechMaterial::Steel);
		Mesh& Br = Mx.MainShaft.Model.Of(MechMaterial::Brass);
		const double Z0 = L.LeftFrame - 4.0, Z1 = L.XRev - 4.5;
		Mesh S = BuildShaft(2.0, Z1 - Z0, false);
		S.Translate(Vec3(0, 0, Z0));
		St.Append(S);
		for (int J = 0; J < Mx.SettingDigits; ++J)
		{
			AddBevel(Br, Vec3(0, 0, L.ColumnX(J)), Vec3(0, 0, -1), L.Drum.BevelTeeth, kDrumModule, 2.0, 2.0, 2.9, 1.0);
		}
		AddBevel(Br, Vec3(0, 0, L.XRev), Vec3(0, 0, -1), 14, 1.0, 3.0, 2.0, 4.0, 1.5);
		Mesh Collar;
		Lathe(Collar, {{2.0, 0.0}, {3.6, 0.0}, {3.6, 1.6}, {3.0, 1.6}, {3.0, 2.6}, {3.6, 2.6}, {3.6, 4.2}, {2.0, 4.2}}, 24, true, SharpAll(8));
		Collar.Translate(Vec3(0, 0, L.XRev - 14.1));
		St.Append(Collar);
		Mx.MainShaft.Poses.push_back(PoseOf(Vec3(0, L.YMain, L.ZDrum), Xw, Zw));
	}
	// Reverser idler on its stud.
	{
		AddBevel(Mx.ReverserIdler.Model.Of(MechMaterial::Brass), Vec3(), Zw, 14, 1.0, 3.0, 1.5, 4.0, 1.5);
		Mesh Stud = BuildShaft(1.5, 7.0, false);
		Stud.Translate(Vec3(0, 0, 6.5));
		Mx.ReverserIdler.Model.Of(MechMaterial::Steel) = Stud;
		Mx.ReverserIdler.Poses.push_back(PoseOf(Vec3(L.XRev, L.YMain, L.ZDrum), Yw, Zw));
		for (const Vec3& A : {Xw, -Xw})
		{
			GearPair B;
			B.Name = "reverser";
			B.Type = GearPair::Bevel;
			B.PointA = B.PointB = Vec3(L.XRev, L.YMain, L.ZDrum);
			B.AxisA = A;
			B.AxisB = Yw;
			B.Module = 1.0;
			B.TeethA = B.TeethB = 14;
			Mx.GearPairs.push_back(B);
		}
	}
	// Crank shaft (the input) out through the drive frame to the drive gear and the crank.
	{
		Mesh& St = Mx.InputShaft.Model.Of(MechMaterial::Steel);
		AddBevel(Mx.InputShaft.Model.Of(MechMaterial::Brass), Vec3(), Zw, 14, 1.0, 3.0, 2.0, 4.0, 1.5);
		Mesh S = BuildShaft(2.5, L.XCrank + 2.0 - L.XRev - 5.0, false);
		S.Translate(Vec3(0, 0, 5.0));
		St.Append(S);
		Mx.InputShaft.Poses.push_back(PoseOf(Vec3(L.XRev, L.YMain, L.ZDrum), Xw, Zw));
		SpokedWheelParams G;
		G.Teeth = 60;
		G.Module = 1.0;
		G.FaceWidth = 6.0;
		G.RimDepth = 3.5;
		G.Spokes = 6;
		G.SpokeWidth = 4.5;
		G.HubRadius = 7.0;
		G.HubLength = 12.0;
		G.BoreRadius = 2.5;
		Mx.DriveGear.Model.Of(MechMaterial::Brass) = BuildSpokedWheel(G);
		Mx.DriveGear.Poses.push_back(PoseOf(Vec3(L.XDrive, L.YMain, L.ZDrum), Xw, Zw));
		// Flywheel: heavy spoked iron wheel with its 20-tooth pinion on the drive gear.
		SpokedWheelParams F;
		F.Teeth = 0;
		F.RimRadius = 30.0;
		F.RimDepth = 6.0;
		F.FaceWidth = 7.0;
		F.Spokes = 6;
		F.SpokeWidth = 4.0;
		F.HubRadius = 6.0;
		F.HubLength = 10.0;
		F.BoreRadius = 2.0;
		Mx.Flywheel.Model.Of(MechMaterial::CastIron) = BuildSpokedWheel(F);
		AddSpur(Mx.Flywheel.Model.Of(MechMaterial::Brass), Vec3(0, 0, L.XDrive - L.XFly), Zw, 20, 1.0, 6.0, 2.0);
		Mesh Fs = BuildShaft(2.0, L.XFly + 4.0 - L.XOuter - 2.0, false);
		Fs.Translate(Vec3(0, 0, L.XOuter + 2.0 - L.XFly));
		Mx.Flywheel.Model.Of(MechMaterial::Steel) = Fs;
		Mx.Flywheel.Poses.push_back(PoseOf(Vec3(L.XFly, L.YMain + 40.0, L.ZDrum), Xw, Zw));
		Mx.FlywheelRatio = 3.0;
		GearPair P;
		P.Name = "drive/flywheel";
		P.PointA = Vec3(L.XDrive, L.YMain, L.ZDrum);
		P.AxisA = Xw;
		P.PointB = Vec3(L.XDrive, L.YMain + 40.0, L.ZDrum);
		P.AxisB = Xw;
		P.Module = 1.0;
		P.TeethA = 60;
		P.TeethB = 20;
		Mx.GearPairs.push_back(P);
	}
	// Crank arm (local z along +x, local y along the arm) and its turned wooden handle.
	{
		Mesh& St = Mx.Crank.Model.Of(MechMaterial::Steel);
		Lathe(St, {{0.0, -1.0}, {5.5, -1.0}, {5.5, 3.0}, {4.5, 4.0}, {0.0, 4.0}}, 28, false, SharpAll(5));
		Prism(St, {{-3.0, 0.0}, {3.0, 0.0}, {1.8, 34.0}, {-1.8, 34.0}}, 1.0, 3.5, Frame3());
		Cylinder(St, 0.0, 34.0, 3.2, 20, 1.0, 3.5);
		Cylinder(St, 0.0, 34.0, 1.2, 12, 3.4, 23.0);
		Mx.Crank.Poses.push_back(PoseOf(Vec3(L.XCrank, L.YMain, L.ZDrum), Xw, Zw));
		Mesh H;
		Lathe(H, {{0.0, 0.0}, {2.6, 0.0}, {3.2, 0.8}, {3.0, 2.0}, {4.2, 4.5}, {4.6, 8.0}, {4.4, 12.0}, {3.6, 15.5}, {4.0, 17.0}, {3.4, 18.6}, {0.0, 19.0}}, 24, false,
			std::vector<char>(11, 0));
		H.Translate(Vec3(0.0, 34.0, 3.8));
		Mx.CrankHandle.Model.Of(MechMaterial::Wood) = H;
		Mx.CrankHandle.Poses.push_back(Mx.Crank.Poses.back());
		Mx.CrankHandleAxis = Zw;
		Mx.CrankHandleCenter = Vec3(0.0, 34.0, 0.0);
	}
	// Reversing lever: pivot on its bracket, fork down to the clutch collar, handle up.
	{
		Mesh& St = Mx.ReversingLever.Model.Of(MechMaterial::Steel);
		const double Pz = 61.5, Down = Pz - (L.ZDrum + 3.6);
		Cylinder(St, 0, 0, 2.2, 16, -1.5, 1.5);
		Prism(St, ThickenPolyline({{0.0, -(Down - 3.0)}, {0.0, 24.0}}, 1.0), -1.0, 1.0, Frame3());
		Prism(St, RectPoly(-1.6, -Down, 1.6, -(Down - 3.5)), -3.0, 3.0, Frame3());
		AddLathe(Mx.ReversingLever.Model.Of(MechMaterial::Wood), {{0.0, 0.0}, {2.5, 0.2}, {3.6, 2.5}, {3.2, 5.0}, {0.0, 5.6}}, Vec3(0, 23.0, 0), Yw, 20);
		Mx.ReversingLever.Poses.push_back(PoseOf(Vec3(L.XRev - 12.0, L.YMain, Pz), Yw, Zw));
		Mx.ReverseAddAngle = 0.26;
		Mx.ReverseSubtractAngle = -0.26;
	}
	// Shift knob, shaft and bevel (local z along +x from the bevel apex) and the pinion under the rack.
	{
		AddBevel(Mx.ShiftKnob.Model.Of(MechMaterial::Brass), Vec3(), Zw, 12, 1.0, 2.5, 1.5, 3.0, 1.0);
		const double Kx = L.XCrank + 6.0 - L.XSp;
		Mesh S = BuildShaft(1.5, Kx - 4.5, false);
		S.Translate(Vec3(0, 0, 4.5));
		Mx.ShiftKnob.Model.Of(MechMaterial::Steel) = S;
		std::vector<Pt2> Knob = {{0.0, 0.0}, {6.2, 0.0}};
		for (int K = 0; K < 6; ++K)
		{
			Knob.push_back(Pt2{K % 2 ? 6.2 : 7.0, 1.0 + K * 1.2});
		}
		Knob.push_back(Pt2{6.2, 8.0});
		Knob.push_back(Pt2{4.5, 9.0});
		Knob.push_back(Pt2{0.0, 9.0});
		AddLathe(Mx.ShiftKnob.Model.Of(MechMaterial::Wood), Knob, Vec3(0, 0, Kx), Zw, 28);
		Mx.ShiftKnob.Poses.push_back(PoseOf(Vec3(L.XSp, L.YAp, L.ZShift), Xw, Zw));
		AddBevel(Mx.ShiftPinion.Model.Of(MechMaterial::Brass), Vec3(), Zw, 12, 1.0, 2.5, 1.5, 3.0, 1.0);
		const double Pl = L.YAp - (L.YAcc + 2.0);
		Mesh S2 = BuildShaft(1.5, Pl + 1.5 - 5.0, false);
		S2.Translate(Vec3(0, 0, 5.0));
		Mx.ShiftPinion.Model.Of(MechMaterial::Steel) = S2;
		AddSpur(Mx.ShiftPinion.Model.Of(MechMaterial::Brass), Vec3(0, 0, Pl), Zw, 14, 1.0, 3.0, 1.5);
		Mx.ShiftPinion.Poses.push_back(PoseOf(Vec3(L.XSp, L.YAp, L.ZShift), -Yw, Zw));
		Mx.ShiftKnobAnglePerStep = L.Pitch / 7.0;
		GearPair B;
		B.Name = "shift bevels";
		B.Type = GearPair::Bevel;
		B.PointA = B.PointB = Vec3(L.XSp, L.YAp, L.ZShift);
		B.AxisA = Xw;
		B.AxisB = -Yw;
		B.Module = 1.0;
		B.TeethA = B.TeethB = 12;
		Mx.GearPairs.push_back(B);
		GearPair R;
		R.Name = "shift pinion/rack";
		R.Type = GearPair::Rack;
		R.PointA = Vec3(L.XSp, L.YAcc + 2.0, L.ZShift);
		R.AxisA = Yw;
		R.PointB = Vec3(L.XSp - L.X0, L.YAcc + 2.0, L.ZShift + 7.0);
		R.AxisB = Xw;
		R.OnCarriageB = true;
		R.Module = 1.0;
		R.TeethA = 14;
		Mx.GearPairs.push_back(R);
	}
	// Counting finger behind the counter wheel at the current shift (its star is on the wheel's right).
	{
		Mesh& F = Mx.CountingFinger.Model.Of(MechMaterial::Brass);
		Cylinder(F, 0, 0, 2.0, 16, -0.5, 0.5);
		Box(F, -0.45, 1.5, -0.5, 0.45, 3.9, 0.5);
		Mx.CountingFinger.Model.Of(MechMaterial::Steel) = BuildShaft(0.7, 2.4, false);
		Mx.CountingFinger.Model.Of(MechMaterial::Steel).Translate(Vec3(0, 0, -0.9));
		Pose P = PoseOf(Vec3(L.X0 + 0.5 * kWheelWidth + 0.55, L.YCnt + 5.0, L.ZCnt + 5.0), Xw, Zw);
		Mx.CountingFinger.Poses.push_back(P);
	}
	// Bell (static dome on a stud) and its hammer, in the drive housing.
	{
		const Vec3 C(L.XRev - 8.0, L.YAcc + 32.0, L.ZBase);
		Mesh Dome;
		Lathe(Dome, {{0.0, 7.8}, {7.3, 1.0}, {8.0, 1.0}, {7.6, 3.5}, {5.5, 6.0}, {3.6, 8.5}, {2.0, 9.4}, {0.0, 9.6}}, 32, false, SharpAll(8));
		Dome.Translate(C);
		Mx.Bell.Of(MechMaterial::Brass) = Dome;
		Mesh Post = BuildShaft(1.2, 9.0, false);
		Post.Translate(C + Vec3(0, 0, -0.1));
		Mx.Bell.Of(MechMaterial::Steel) = Post;
		Mesh& Hs = Mx.BellHammer.Model.Of(MechMaterial::Steel);
		Cylinder(Hs, 0, 0, 1.4, 12, -1.5, 1.5);
		Prism(Hs, ThickenPolyline({{0.0, 0.0}, {-3.0, 3.5}, {-4.4, 5.4}}, 0.45), -0.45, 0.45, Frame3());
		AddLathe(Mx.BellHammer.Model.Of(MechMaterial::Brass), {{0.0, -1.8}, {1.3, -1.3}, {1.8, 0.0}, {1.3, 1.3}, {0.0, 1.8}}, Vec3(-5.2, 6.0, 0.0), Zw, 14);
		Mx.BellHammer.Poses.push_back(PoseOf(C + Vec3(-13.5, 0, 2.0), Yw, Zw));
		Mx.BellStrikeAngle = 0.22;
	}
	Mx.CarryLeverSwing = 0.3;
	Mx.ClearingLeverSwing = 2.0 * kPi;
}

void BuildPrinter(MechanismParts& Mx, const MechLayout& L)
{
	const Vec3 Xw(1, 0, 0), Yw(0, 1, 0), Zw(0, 0, 1);
	const double Xc = 0.5 * (L.PlinthX0 + std::min(L.KbMinX, -L.Hw)), Yp = L.PlinthY0 + 12.0;
	const double YTw = Yp + 9.0, ZPr = 44.0, RTw = 7.5;
	const double YPl = YTw + RTw + 0.5 + 8.0, RPl = 8.0;
	const double YRoll = YPl + 30.0, ZRoll = 38.0, RRoll = 17.0;
	const double YRib = YTw + RTw + 0.25;
	const double Paper = 0.1, Width = 57.0;
	Mesh& E = Mx.Printer.Of(MechMaterial::BlackEnamel);
	Mesh& St = Mx.Printer.Of(MechMaterial::Steel);
	// Frame plates with bosses for the platen, the roll and the type-wheel shaft.
	const std::vector<Pt2> Side = {{Yp, 0.0}, {YRoll + 20.0, 0.0}, {YRoll + 20.0, 8.0}, {YRoll + 14.0, 14.0}, {YRoll + 10.0, ZRoll + 14.0}, {YPl + 6.0, ZPr + 12.0},
		{YTw - 4.0, ZPr + 12.0}, {Yp, ZPr + 2.0}};
	for (double S = -1.0; S <= 1.0; S += 2.0)
	{
		const double X = Xc + S * 34.0;
		PrismWithHoles(E, Side, {CircleHole(0.5 * (YPl + YRoll), 18.0, 5.0, 18)}, X - 1.25, X + 1.25, AxisFrame(Vec3(), Xw, Zw));
		for (const Vec3& B : {Vec3(X + S * 1.25, YTw, ZPr), Vec3(X + S * 1.25, YPl, ZPr), Vec3(X + S * 1.25, YRoll, ZRoll)})
		{
			AddBearingBoss(E, St, B, Xw * S, 4.0);
		}
	}
	AddRod(St, Vec3(Xc - 35.0, YTw, ZPr), Vec3(Xc + 35.0, YTw, ZPr), 1.0);
	AddRod(St, Vec3(Xc - 35.0, Yp + 1.0, 8.0), Vec3(Xc + 35.0, Yp + 1.0, 8.0), 1.6);
	AddRod(St, Vec3(Xc - 35.0, YRoll + 18.0, 8.0), Vec3(Xc + 35.0, YRoll + 18.0, 8.0), 1.6);
	// Ribbon between the spools, passing between the type wheels and the paper.
	const double Sl = Xc - 39.0, Sr = Xc + 39.0;
	AddBox(Mx.Printer.Of(MechMaterial::Ribbon), Sl, Sr, YRib - 0.1, YRib + 0.1, ZPr - 4.0, ZPr + 4.0);
	// Tear bar above the platen, behind the paper.
	{
		std::vector<Pt2> Teeth = {{Xc - 32.0, 0.0}, {Xc + 32.0, 0.0}};
		for (int K = 31; K >= 0; --K)
		{
			const double X = Xc - 32.0 + K * 2.0;
			Teeth.push_back(Pt2{X + 2.0, 1.2});
			Teeth.push_back(Pt2{X + 1.0, 1.8});
		}
		Teeth.push_back(Pt2{Xc - 32.0, 1.2});
		Frame3 Tf;
		Tf.O = Vec3(0, 0, ZPr + RPl + 4.0);
		Tf.U = Xw;
		Tf.V = Zw;
		Tf.W = Yw;
		Prism(St, Teeth, YPl - RPl + 0.3, YPl - RPl + 1.3, Tf);
		AddBox(St, Xc - 34.0, Xc + 34.0, YPl - RPl + 0.3, YPl - RPl + 1.3, ZPr + RPl + 1.0, ZPr + RPl + 4.1);
	}
	// Type wheels: brass wheels with raised figures, turned to the platen (local +Y -> world +y).
	{
		NumeralWheelParams P;
		P.Radius = RTw;
		P.Width = 2.4;
		P.GlyphsAlongAxis = false;
		P.GlyphHeight = 3.0;
		P.GlyphWidth = 1.9;
		P.StrokeMm = 0.4;
		P.EmbossMm = 0.5;
		P.BoreRadius = 1.05;
		const NumeralWheel W = BuildNumeralWheel(P);
		SetWheelModel(Mx.TypeWheels.Model, W, MechMaterial::Brass, MechMaterial::Engraving);
		for (int K = 0; K < 18; ++K)
		{
			Mx.TypeWheels.Poses.push_back(PoseOf(Vec3(Xc - 25.5 + 3.0 * K, YTw, ZPr), Xw, Yw));
		}
	}
	// Platen, roll and ribbon spools.
	{
		Mesh Pl;
		Lathe(Pl, {{0.0, -32.0}, {RPl, -32.0}, {RPl, 32.0}, {0.0, 32.0}}, 40, false, SharpAll(4));
		Mx.Platen.Model.Of(MechMaterial::Rubber) = Pl;
		Mx.Platen.Model.Of(MechMaterial::Steel) = BuildShaft(2.0, 70.0, false);
		Mx.Platen.Model.Of(MechMaterial::Steel).Translate(Vec3(0, 0, -35.0));
		AddLathe(Mx.Platen.Model.Of(MechMaterial::Wood), {{0.0, 0.0}, {5.0, 0.0}, {6.0, 2.0}, {6.0, 5.0}, {4.5, 6.5}, {0.0, 6.5}}, Vec3(0, 0, 36.5), Zw, 24);
		Mx.Platen.Poses.push_back(PoseOf(Vec3(Xc, YPl, ZPr), Xw, Zw));
		Mx.PlatenRadius = RPl;
		Mesh Roll;
		Lathe(Roll, {{5.6, -0.5 * Width}, {RRoll, -0.5 * Width}, {RRoll, 0.5 * Width}, {5.6, 0.5 * Width}}, 48, true, SharpAll(4));
		Mx.PaperRoll.Model.Of(MechMaterial::Paper) = Roll;
		Mesh Core;
		Lathe(Core, {{2.1, -30.0}, {5.7, -30.0}, {5.7, 30.0}, {2.1, 30.0}}, 32, true, SharpAll(4));
		Mx.PaperRoll.Model.Of(MechMaterial::Brass) = Core;
		Mx.PaperRoll.Model.Of(MechMaterial::Steel) = BuildShaft(2.0, 70.0, false);
		Mx.PaperRoll.Model.Of(MechMaterial::Steel).Translate(Vec3(0, 0, -35.0));
		Mx.PaperRoll.Poses.push_back(PoseOf(Vec3(Xc, YRoll, ZRoll), Xw, Zw));
		Mx.PaperRollRadius = RRoll;
		Mesh& Sp = Mx.RibbonSpools.Model.Of(MechMaterial::Brass);
		Lathe(Sp, {{0.0, -6.8}, {5.5, -6.8}, {5.5, -6.0}, {2.2, -6.0}, {2.2, 6.0}, {5.5, 6.0}, {5.5, 6.8}, {0.0, 6.8}}, 24, false, SharpAll(8));
		Mesh Wound;
		Lathe(Wound, {{2.1, -4.0}, {4.5, -4.0}, {4.5, 4.0}, {2.1, 4.0}}, 24, true, SharpAll(4));
		Mx.RibbonSpools.Model.Of(MechMaterial::Ribbon) = Wound;
		for (double X : {Sl, Sr})
		{
			Mx.RibbonSpools.Poses.push_back(PoseOf(Vec3(X, YRib - 4.5, ZPr), Zw, Yw));
			AddRod(St, Vec3(X, YRib - 4.5, -0.1), Vec3(X, YRib - 4.5, ZPr - 6.7), 1.5);
		}
	}
	// The paper: off the roll, round the front of the platen past the print line, up and curling
	// back over the top. Path in (y, z) along the paper's centre.
	std::vector<Pt2> Path;
	const double Rc = RPl + 0.5 * Paper;
	const double A0 = 235.0 * kGearDeg;
	Path.push_back(Pt2{YRoll + (RRoll + 0.5 * Paper) * std::cos(A0), ZRoll + (RRoll + 0.5 * Paper) * std::sin(A0)});
	for (int K = 0; K <= 12; ++K)
	{
		const double A = (270.0 - 90.0 * K / 12.0) * kGearDeg;
		Path.push_back(Pt2{YPl + Rc * std::cos(A), ZPr + Rc * std::sin(A)});
	}
	const size_t PrintLine = Path.size() - 1;
	const double YUp = YPl - Rc, ZUp = ZPr + 62.0, Curl = 20.0;
	for (double Z = ZPr + 2.0; Z < ZUp; Z += 2.0)
	{
		Path.push_back(Pt2{YUp, Z});
	}
	for (int K = 0; K <= 24; ++K)
	{
		const double A = kPi - kPi * K / 24.0;
		Path.push_back(Pt2{YUp + Curl + Curl * std::cos(A), ZUp + Curl * std::sin(A)});
	}
	for (double Dz = 2.0; Dz <= 16.0; Dz += 2.0)
	{
		Path.push_back(Pt2{YUp + 2.0 * Curl, ZUp - Dz});
	}
	Prism(Mx.PaperStrip.Of(MechMaterial::Paper), ThickenPolyline(Path, 0.5 * Paper), Xc - 0.5 * Width, Xc + 0.5 * Width, AxisFrame(Vec3(), Xw, Zw));
	for (size_t K = PrintLine; K < Path.size(); ++K)
	{
		const Pt2 Prev = Path[K > 0 ? K - 1 : K], Next = Path[K + 1 < Path.size() ? K + 1 : K];
		const double Ty = Next.X - Prev.X, Tz = Next.Y - Prev.Y, Tl = std::sqrt(Ty * Ty + Tz * Tz);
		const Vec3 N = Vec3(0, -Tz / Tl, Ty / Tl);
		Mx.TapeNormals.push_back(N);
		Mx.TapePath.push_back(Vec3(Xc, Path[K].X, Path[K].Y) + N * (0.5 * Paper));
	}
	Mx.TapeWidthMm = Width;
	Mx.TapeLinePitchMm = 4.23;
	Mx.TapeUsableWidthMm = 50.0;
}

constexpr int kMemTopTeeth = 24, kMemIdlerTeeth = 18;

void BuildMemory(MechanismParts& Mx, const MechLayout& L)
{
	const Vec3 Yw(0, 1, 0), Zw(0, 0, 1);
	Mesh& Fe = Mx.MemoryFrame.Of(MechMaterial::CastIron);
	Mesh& St = Mx.MemoryFrame.Of(MechMaterial::Steel);
	Mesh& Br = Mx.MemoryFrame.Of(MechMaterial::Brass);
	const double R = 8.0, Pitch = 7.2, Base = 12.0, Gap = 5.0;
	const int N = Mx.MemoryColumns, Per = Mx.MemoryMantissa + Mx.MemoryExponent;
	const double RailX0 = L.PlinthX0 + 16.0, RailX1 = L.LeftFrame - 18.0;
	const double Mid = 0.5 * (RailX0 + RailX1), Span = 0.5 * (RailX1 - RailX0) - 18.0 - 15.0;
	const double TopZ = Base + 4.0 + (Per + 1) * Pitch + Gap + 4.0;
	NumeralWheelParams P;
	P.Radius = R;
	P.Width = 4.4;
	P.GlyphsAlongAxis = true;
	P.GlyphHeight = 3.4;
	P.GlyphWidth = 2.4;
	P.StrokeMm = 0.45;
	P.EmbossMm = 0.3;
	P.BoreRadius = 1.25;
	P.GearTeeth = 22;            // each figure wheel's own gear, under its rim
	P.GearModule = 0.6;
	SetWheelModel(Mx.MemoryWheels.Model, BuildNumeralWheel(P), MechMaterial::Brass, MechMaterial::Engraving);
	P.SignWheel = true;
	SetWheelModel(Mx.MemorySignWheels.Model, BuildNumeralWheel(P), MechMaterial::Brass, MechMaterial::Engraving);
	std::vector<ReadingWindow> Signs;
	for (int C = 0; C < N; ++C)
	{
		const double X = Mid - Span + 2.0 * Span * C / (N - 1);
		const double Y = L.YMem;
		// Base block, axle, two pillars and the label plate.
		SolidFromLoop(Fe, RoundedRect(X, Y, 13.0, 13.0, 2.0, 4, Splits()), 0.0, Base);
		SolidFromLoop(Fe, RoundedRect(X, Y, 11.0, 11.0, 2.0, 4, Splits()), Base - 0.1, Base + 2.0);
		AddRod(St, Vec3(X, Y, Base + 1.9), Vec3(X, Y, TopZ + 11.5), 1.2);
		AddSpur(Br, Vec3(X, Y, TopZ + 9.5), Zw, kMemTopTeeth, 1.0, 2.4, 1.25);
		for (double S = -1.0; S <= 1.0; S += 2.0)
		{
			AddRod(Br, Vec3(X + S * 11.0, Y + 8.0, Base + 1.9), Vec3(X + S * 11.0, Y + 8.0, TopZ + 0.1), 1.4);
		}
		Mesh Plate;
		Prism(Plate, RectPoly(-11.0, 1.2, 11.0, 8.8), -0.4, 0.6, Frame3());
		Frame3 F;
		F.O = Vec3(X, Y - 13.0, 0.0);
		F.U = Vec3(1, 0, 0);
		F.V = Zw;
		F.W = -Yw;
		ApplyFrame(Plate, F);
		Br.Append(Plate);
		LabelFace Lf;
		Lf.Center = Vec3(X, Y - 13.6, 5.0);
		Lf.Normal = -Yw;
		Lf.Up = Zw;
		Lf.WidthMm = 20.0;
		Lf.HeightMm = 6.4;
		Mx.MemoryLabels.push_back(Lf);
		// Wheels from the top down: sign, mantissa, a gap, the exponent.
		auto Place = [&](double Z) {
			Pose Wp = PoseOf(Vec3(X, Y, Z), Zw, -Yw);
			return Wp;
		};
		double Z = TopZ - 4.0 - 0.5 * Pitch;
		const Pose Sign = Place(Z);
		Mx.MemorySignWheels.Poses.push_back(Sign);
		Signs.push_back(WindowAt(Sign.Point(Vec3(0, R + 0.3, 0)), 2.6, 3.6));
		for (int K = 0; K < Per; ++K)
		{
			Z -= Pitch + (K == Mx.MemoryMantissa ? Gap : 0.0);
			const Pose Wp = Place(Z);
			Mx.MemoryWheels.Poses.push_back(Wp);
			Mx.MemoryWindows.push_back(WindowAt(Wp.Point(Vec3(0, R + 0.3, 0)), 2.6, 3.6));
		}
		// Brass spacer washers between the mantissa and the exponent.
		Cylinder(Br, X, Y, 3.0, 16, Z + 0.5 * Pitch + 0.2, Z + 0.5 * Pitch + Gap - 0.2);
	}
	Mx.MemoryWindows.insert(Mx.MemoryWindows.end(), Signs.begin(), Signs.end());
	// The top rail across all the columns, on end posts.
	const double X0 = RailX0, X1 = RailX1;
	AddBox(Fe, X0, X1, L.YMem + 2.0, L.YMem + 14.0, TopZ, TopZ + 7.0);
	const double Half = Span / (N - 1), Cd = 0.5 * (kMemTopTeeth + kMemIdlerTeeth);
	if (Half < Cd - 1.0)
	{
		const double Dy = std::sqrt(Cd * Cd - Half * Half);
		for (int C = 0; C + 1 < N; ++C)
		{
			const double Xi = Mid - Span + Half * (2 * C + 1);
			const Vec3 Ci(Xi, L.YMem + Dy, TopZ + 9.5);
			AddRod(St, Vec3(Xi, Ci.Y, TopZ + 6.9), Vec3(Xi, Ci.Y, TopZ + 11.5), 1.2);
			AddSpur(Br, Ci, Zw, kMemIdlerTeeth, 1.0, 2.4, 1.25);
			for (int Side = 0; Side < 2; ++Side)
			{
				GearPair G;
				G.Name = "memory top gears";
				G.PointA = Vec3(Mid - Span + 2.0 * Half * (C + Side), L.YMem, TopZ + 9.5);
				G.AxisA = Zw;
				G.PointB = Ci;
				G.AxisB = Zw;
				G.Module = 1.0;
				G.TeethA = kMemTopTeeth;
				G.TeethB = kMemIdlerTeeth;
				Mx.GearPairs.push_back(G);
			}
		}
	}
	for (double X : {X0 + 4.0, X1 - 4.0})
	{
		AddBox(Fe, X - 4.0, X + 4.0, L.YMem + 2.0, L.YMem + 14.0, 0.0, TopZ + 0.1);
	}
}
} // namespace gearparts_detail

MechanismParts BuildMechanism(const mkb::Layout& Layout, const CaseParams& Case)
{
	using namespace gearparts_detail;
	MechanismParts Mx;
	const MechLayout L = Plan(Layout, Case);
	Mx.DigitPitchMm = L.Pitch;
	Mx.PinionStepMm = L.Step;
	Mx.DigitAngle = DigitAngleRad(1);
	BuildPlinth(Mx, L);
	BuildFrames(Mx, L);
	BuildColumns(Mx, L);
	BuildCarriage(Mx, L);
	BuildDrive(Mx, L);
	BuildPrinter(Mx, L);
	BuildMemory(Mx, L);
	Mx.BoundsMin = Vec3(L.PlinthX0 + 10.0, L.PlinthY0 + 10.0, 0.0);
	Mx.BoundsMax = Vec3(L.XCrank + 25.0, L.YMain + 52.0, 200.0);
	return Mx;
}

TapeFrame TapeLineFrame(const MechanismParts& M, double DistanceMm)
{
	TapeFrame F;
	const std::vector<Vec3>& P = M.TapePath;
	if (P.size() < 2)
	{
		return F;
	}
	double Left = std::max(0.0, DistanceMm);
	size_t K = 0;
	while (K + 2 < P.size() && Left > (P[K + 1] - P[K]).Length())
	{
		Left -= (P[K + 1] - P[K]).Length();
		++K;
	}
	const Vec3 Seg = P[K + 1] - P[K];
	const double Len = Seg.Length();
	const double T = Len > 0.0 ? Left / Len : 0.0;
	F.Up = Seg.Normalized();
	F.Origin = P[K] + Seg * T;
	const Vec3 N = M.TapeNormals[K] * (1.0 - std::min(1.0, T)) + M.TapeNormals[K + 1] * std::min(1.0, T);
	F.Normal = (N - F.Up * N.Dot(F.Up)).Normalized();
	F.Right = F.Up.Cross(F.Normal).Normalized();
	return F;
}

} // namespace mgeo
