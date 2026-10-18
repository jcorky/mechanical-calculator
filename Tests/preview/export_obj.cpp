// Writes the assembled keyboard for the preview renderer (Tests/preview/render.mjs):
//   board.obj            the board with keycaps on switches, one group per material and one per keycap
//   switch_exploded.obj  a single tactile switch and its keycap, pulled apart along z
//   scene.json           legends and colour roles per keycap, the rest pose and the display frame
//   mechanism.obj        each unique mechanism part once, one group per part and material
//   mechanism.json       the instances (world poses, digits set, carriage shifted), the exploded
//                        column, label faces and the printed tape lines
// Usage: mcalc_export_obj [output directory]
#include "GearParts.h"
#include "KeyboardParts.h"

#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace exportobj_detail
{
using mgeo::Mesh;
using mgeo::Vec3;

const char* RoleName(mkb::KeyRole R)
{
	switch (R)
	{
	case mkb::KeyRole::Alpha:
		return "alpha";
	case mkb::KeyRole::Modifier:
		return "modifier";
	case mkb::KeyRole::Accent:
		return "accent";
	case mkb::KeyRole::Function:
		return "function";
	case mkb::KeyRole::Navigation:
		return "navigation";
	case mkb::KeyRole::Numpad:
		return "numpad";
	case mkb::KeyRole::Space:
		return "space";
	}
	return "alpha";
}

std::string Json(const std::string& S)
{
	std::string Out = "\"";
	for (char C : S)
	{
		if (C == '"' || C == '\\')
		{
			Out += '\\';
		}
		Out += C;
	}
	return Out + "\"";
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

// Named meshes in output order; Get() appends to an existing group.
struct Groups
{
	std::vector<std::pair<std::string, std::unique_ptr<Mesh>>> List;
	Mesh& Get(const std::string& Name)
	{
		for (auto& G : List)
		{
			if (G.first == Name)
			{
				return *G.second;
			}
		}
		List.emplace_back(Name, std::make_unique<Mesh>());
		return *List.back().second;
	}
	bool Write(const std::string& Path) const
	{
		std::vector<std::pair<std::string, const Mesh*>> Parts;
		for (const auto& G : List)
		{
			Parts.emplace_back(G.first, G.second.get());
		}
		return mgeo::WriteObj(Path, Parts);
	}
	int Triangles() const
	{
		int N = 0;
		for (const auto& G : List)
		{
			N += G.second->TriangleCount();
		}
		return N;
	}
};

void AddSwitch(Groups& G, const mgeo::SwitchParts& S, const Vec3& At, double StemZ, const double Explode[6])
{
	G.Get("switch_housing_bottom").Append(S.HousingBottom, At + Vec3(0, 0, Explode[0]));
	G.Get("switch_pins").Append(S.Pins, At + Vec3(0, 0, Explode[0]));
	G.Get("switch_contacts").Append(S.Contacts, At + Vec3(0, 0, Explode[1]));
	G.Get("switch_spring").Append(S.Spring, At + Vec3(0, 0, Explode[2]));
	G.Get("switch_housing_top").Append(S.HousingTop, At + Vec3(0, 0, Explode[3]));
	G.Get("switch_stem").Append(S.Stem, At + Vec3(0, 0, StemZ + Explode[4]));
}

// --- The mechanism -------------------------------------------------------------------------------

std::string PoseJson(const mgeo::Pose& P)
{
	char Buf[400];
	std::snprintf(Buf, sizeof(Buf), "[%.4f,%.4f,%.4f, %.6f,%.6f,%.6f, %.6f,%.6f,%.6f, %.6f,%.6f,%.6f]", P.Origin.X, P.Origin.Y, P.Origin.Z, P.X.X, P.X.Y, P.X.Z, P.Y.X, P.Y.Y,
		P.Y.Z, P.Z.X, P.Z.Y, P.Z.Z);
	return Buf;
}

std::string VecJson(const Vec3& V)
{
	char Buf[120];
	std::snprintf(Buf, sizeof(Buf), "[%.4f,%.4f,%.4f]", V.X, V.Y, V.Z);
	return Buf;
}

std::string FaceJson(const mgeo::LabelFace& F, const std::string& Text)
{
	return "{\"center\": " + VecJson(F.Center) + ", \"normal\": " + VecJson(F.Normal) + ", \"up\": " + VecJson(F.Up) + ", \"w\": " + std::to_string(F.WidthMm) +
		", \"h\": " + std::to_string(F.HeightMm) + ", \"text\": " + Json(Text) + "}";
}

// Digit k of a number written as a string, counted from its last character (0 if past the start).
int DigitFromRight(const std::string& S, size_t K)
{
	return K < S.size() ? S[S.size() - 1 - K] - '0' : 0;
}

mgeo::Pose Lifted(mgeo::Pose P, const Vec3& By)
{
	P.Origin = P.Origin + By;
	return P;
}

struct MechExport
{
	Groups Obj;
	// JSON array entries {"name", "groups": [...], "poses": [...]}: the machine, the exploded column,
	// and the optional parts (not shown).
	std::string Parts, Exploded, Optional;
	void Add(std::string& Into, const std::string& Name, const mgeo::MechPart& Part, const std::vector<mgeo::Pose>& Poses)
	{
		std::string GroupList;
		for (const mgeo::MechSection& S : Part.Sections)
		{
			const std::string G = "m:" + Name + ":" + mgeo::MechMaterialName(S.Material);
			Mesh& M = Obj.Get(G);
			if (M.Positions.empty())
			{
				M.Append(S.Geometry);
			}
			GroupList += (GroupList.empty() ? "" : ", ") + Json(G);
		}
		std::string PoseList;
		for (const mgeo::Pose& P : Poses)
		{
			PoseList += (PoseList.empty() ? "" : ", ") + PoseJson(P);
		}
		Into += std::string(Into.empty() ? "" : ",") + "\n    {\"name\": " + Json(Name) + ", \"groups\": [" + GroupList + "], \"poses\": [" + PoseList + "]}";
	}
};

bool WriteMechanism(const std::string& Dir, const mkb::Layout& L, const mgeo::CaseParams& CP)
{
	using mgeo::Pose;
	const mgeo::MechanismParts M = mgeo::BuildMechanism(L, CP);
	const int Shift = 3;
	// What the registers show: 271828 set, 271828 x 315 accumulated, the carriage at shift 3.
	const std::string Set = "271828", Acc = "85625820", Cnt = "315";
	const Pose Car = M.Carriage.Poses[0].Moved(Vec3(Shift * M.DigitPitchMm, 0, 0));
	auto Turned = [](const std::vector<Pose>& Poses, const std::string& Digits) {
		std::vector<Pose> Out;
		for (size_t I = 0; I < Poses.size(); ++I)
		{
			Out.push_back(Poses[I].RotatedZ(mgeo::DigitAngleRad(DigitFromRight(Digits, I))));
		}
		return Out;
	};
	auto OnCarriage = [](const Pose& At, std::vector<Pose> Poses) {
		for (Pose& P : Poses)
		{
			P = At.Then(P);
		}
		return Poses;
	};
	MechExport X;
	std::string& Main = X.Parts;
	const Pose Id;
	X.Add(Main, "plinth", M.Plinth, {Id});
	X.Add(Main, "frame", M.Frame, {Id});
	X.Add(Main, "printer", M.Printer, {Id});
	X.Add(Main, "paper", M.PaperStrip, {Id});
	X.Add(Main, "bell", M.Bell, {Id});
	X.Add(Main, "memory_frame", M.MemoryFrame, {Id});
	X.Add(X.Optional, "glass", M.GlassCover, {Id});
	X.Add(Main, "drum", M.Drums.Model, M.Drums.Poses);
	X.Add(Main, "square_shaft", M.SquareShafts.Model, M.SquareShafts.Poses);
	std::vector<Pose> Pin;
	for (size_t J = 0; J < M.Pinions.Poses.size(); ++J)
	{
		Pin.push_back(M.Pinions.Poses[J].Moved(Vec3(0, 0, DigitFromRight(Set, J) * M.PinionStepMm)));
	}
	X.Add(Main, "pinion", M.Pinions.Model, Pin);
	X.Add(Main, "transfer", M.TransferGears.Model, M.TransferGears.Poses);
	X.Add(Main, "setting_wheel", M.SettingWheels.Model, Turned(M.SettingWheels.Poses, Set));
	X.Add(Main, "setting_marker", M.SettingMarker.Model, M.SettingMarker.Poses);
	X.Add(Main, "carriage", M.Carriage.Model, {Car});
	X.Add(Main, "accumulator_wheel", M.AccumulatorWheels.Model, OnCarriage(Car, Turned(M.AccumulatorWheels.Poses, Acc)));
	X.Add(Main, "carry_lever", M.CarryLevers.Model, OnCarriage(Car, M.CarryLevers.Poses));
	X.Add(Main, "counter_wheel", M.CounterWheels.Model, OnCarriage(Car, Turned(M.CounterWheels.Poses, Cnt)));
	X.Add(Main, "carriage_marker", M.CarriageMarkers.Model, OnCarriage(Car, M.CarriageMarkers.Poses));
	X.Add(Main, "clearing_lever", M.ClearingLevers.Model, OnCarriage(Car, M.ClearingLevers.Poses));
	const mgeo::PartSet* Plain[] = {&M.CountingFinger, &M.MainShaft, &M.ReverserIdler, &M.InputShaft, &M.DriveGear, &M.Flywheel, &M.Crank, &M.CrankHandle,
		&M.ReversingLever, &M.ShiftKnob, &M.ShiftPinion, &M.BellHammer, &M.PaperRoll, &M.Platen, &M.RibbonSpools};
	const char* PlainNames[] = {"counting_finger", "main_shaft", "reverser_idler", "input_shaft", "drive_gear", "flywheel", "crank", "crank_handle", "reversing_lever",
		"shift_knob", "shift_pinion", "bell_hammer", "paper_roll", "platen", "ribbon_spool"};
	for (size_t I = 0; I < sizeof(Plain) / sizeof(Plain[0]); ++I)
	{
		X.Add(Main, PlainNames[I], Plain[I]->Model, Plain[I]->Poses);
	}
	// Type wheels run left to right, so the number is laid on them from the last wheel.
	std::vector<Pose> Type(M.TypeWheels.Poses.rbegin(), M.TypeWheels.Poses.rend());
	Type = Turned(Type, "85625820");
	X.Add(Main, "type_wheel", M.TypeWheels.Model, std::vector<Pose>(Type.rbegin(), Type.rend()));
	// Memory: ans and five constants, 18 mantissa digits and a signed 3-digit exponent each.
	const char* Mantissa[6] = {"856258200000000000", "314159265358979323", "271828182845904523", "141421356237309504", "602214076000000000", "980665000000000000"};
	const char* Exponent[6] = {"007", "000", "000", "000", "023", "000"};
	const char* Names[6] = {"ANS", "PI", "E", "ROOT 2", "N A", "G"};
	std::vector<Pose> Mem, Signs;
	for (int C = 0; C < M.MemoryColumns; ++C)
	{
		for (int K = 0; K < M.MemoryMantissa + M.MemoryExponent; ++K)
		{
			const char D = K < M.MemoryMantissa ? Mantissa[C][K] : Exponent[C][K - M.MemoryMantissa];
			Mem.push_back(M.MemoryWheels.Poses[static_cast<size_t>(C * 21 + K)].RotatedZ(mgeo::DigitAngleRad(D - '0')));
		}
		Signs.push_back(M.MemorySignWheels.Poses[static_cast<size_t>(C)].RotatedZ(mgeo::DigitAngleRad(0)));
	}
	X.Add(Main, "memory_wheel", M.MemoryWheels.Model, Mem);
	X.Add(Main, "sign_wheel", M.MemorySignWheels.Model, Signs);

	// The exploded column: column 0's chain from the drum to its carriage wheels, pulled apart.
	{
		std::string& E = X.Exploded;
		const Pose Car0 = M.Carriage.Poses[0];
		X.Add(E, "drum", M.Drums.Model, {M.Drums.Poses[0]});
		X.Add(E, "square_shaft", M.SquareShafts.Model, {Lifted(M.SquareShafts.Poses[0], Vec3(0, 0, 24))});
		X.Add(E, "pinion", M.Pinions.Model, {Lifted(M.Pinions.Poses[0].Moved(Vec3(0, 0, 5 * M.PinionStepMm)), Vec3(0, 0, 13))});
		X.Add(E, "transfer", M.TransferGears.Model, {Lifted(M.TransferGears.Poses[0], Vec3(0, -18, 28))});
		X.Add(E, "setting_wheel", M.SettingWheels.Model, {Lifted(M.SettingWheels.Poses[0].RotatedZ(mgeo::DigitAngleRad(5)), Vec3(0, 0, 34))});
		X.Add(E, "accumulator_wheel", M.AccumulatorWheels.Model, {Lifted(Car0.Then(M.AccumulatorWheels.Poses[0].RotatedZ(mgeo::DigitAngleRad(5))), Vec3(0, -38, 24))});
		X.Add(E, "carry_lever", M.CarryLevers.Model, {Lifted(Car0.Then(M.CarryLevers.Poses[0]), Vec3(0, -38, -4))});
		X.Add(E, "counter_wheel", M.CounterWheels.Model, {Lifted(Car0.Then(M.CounterWheels.Poses[0].RotatedZ(mgeo::DigitAngleRad(1))), Vec3(0, -38, 42))});
	}
	if (!X.Obj.Write(Dir + "/mechanism.obj"))
	{
		return false;
	}
	std::string Labels = FaceJson(M.MakersPlate, "MECHANICAL CALCULATOR \xC2\xB7 WHERE IT ISN'T");
	for (int C = 0; C < M.MemoryColumns; ++C)
	{
		Labels += ",\n    " + FaceJson(M.MemoryLabels[static_cast<size_t>(C)], Names[C]);
	}
	// The tape: the newest line at the print line, older ones above it.
	const char* Lines[] = {"85625820 *", "315 =", "271828 X", "0 C", "9.80665 *", "6.02214076 E23 *", "1.41421356 *", "2.71828182 *", "3.14159265 *"};
	std::string Tape;
	for (size_t I = 0; I < sizeof(Lines) / sizeof(Lines[0]); ++I)
	{
		const mgeo::TapeFrame F = mgeo::TapeLineFrame(M, (static_cast<double>(I) + 0.5) * M.TapeLinePitchMm);
		Tape += std::string(Tape.empty() ? "" : ",") + "\n    {\"origin\": " + VecJson(F.Origin) + ", \"up\": " + VecJson(F.Up) + ", \"normal\": " + VecJson(F.Normal) +
			", \"right\": " + VecJson(F.Right) + ", \"text\": " + Json(Lines[I]) + "}";
	}
	std::FILE* F = std::fopen((Dir + "/mechanism.json").c_str(), "w");
	if (!F)
	{
		return false;
	}
	std::fprintf(F,
		"{\n  \"plinthH\": %.4f,\n  \"tapeWidth\": %.4f, \"tapeUsable\": %.4f, \"tapePitch\": %.4f,\n  \"labels\": [\n    %s\n  ],\n  \"tape\": [%s\n  ],\n  \"parts\": [%s\n  ],\n"
		"  \"exploded\": [%s\n  ],\n  \"optional\": [%s\n  ]\n}\n",
		M.PlinthHeightMm, M.TapeWidthMm, M.TapeUsableWidthMm, M.TapeLinePitchMm, Labels.c_str(), Tape.c_str(), X.Parts.c_str(), X.Exploded.c_str(), X.Optional.c_str());
	std::fclose(F);
	std::printf("mechanism.obj: %d unique triangles in %zu groups\n", X.Obj.Triangles(), X.Obj.List.size());
	return true;
}
} // namespace exportobj_detail

int main(int Argc, char** Argv)
{
	using namespace exportobj_detail;
	const std::string Dir = Argc > 1 ? Argv[1] : ".";
	const mkb::Layout& L = mkb::FullSizeAnsi();
	const mgeo::CaseParams CP;
	const mgeo::SwitchParams SP;
	const mgeo::BoardParts B = mgeo::BuildBoard(L, CP);
	const mgeo::SwitchParts Sw = mgeo::BuildSwitch(SP);
	const double None[6] = {0, 0, 0, 0, 0, 0};

	Groups Board;
	Board.Get("aluminium").Append(B.CaseShell);
	Board.Get("badge").Append(B.CaseBadge);
	Board.Get("rubber").Append(B.Feet);
	Board.Get("plate").Append(B.Plate);
	Board.Get("pcb").Append(B.Pcb);
	Board.Get("pcb_parts").Append(B.PcbParts);
	Board.Get("oled_bezel").Append(B.DisplayBezel);
	Board.Get("oled_glass").Append(B.DisplayGlass);
	Board.Get("nameplate").Append(B.NamePlate);
	std::map<std::tuple<double, double, int, bool, bool>, Mesh> Caps;
	std::string Keys;
	for (const mkb::KeyDef& K : L.Keys)
	{
		const Vec3 C = mgeo::KeyCenter(K, L);
		AddSwitch(Board, Sw, C, SP.StemTopZ, None);
		if (K.NeedsStabilizer())
		{
			const bool Vertical = K.H >= 2.0 && K.W < 2.0;
			const mgeo::StabilizerParts S = mgeo::BuildStabilizer(Vertical ? K.H : K.W, Vertical);
			Board.Get("stab_housing").Append(S.Housings, C);
			Board.Get("stab_stem").Append(S.Stems, C + Vec3(0, 0, SP.StemTopZ));
			Board.Get("stab_wire").Append(S.Wire, C);
		}
		const mgeo::KeycapParams P = CapFor(K);
		const auto Key = std::make_tuple(P.WidthU, P.DepthU, P.ProfileRow, P.Spacebar, P.HomingBar);
		if (!Caps.count(Key))
		{
			Caps[Key] = mgeo::BuildKeycap(P);
		}
		Board.Get(std::string("cap:") + RoleName(K.Role) + ":" + K.Id).Append(Caps[Key], C + Vec3(0, 0, SP.StemTopZ));
		const mgeo::KeycapTop T = mgeo::KeycapTopFor(P);
		char Buf[768];
		std::snprintf(Buf, sizeof(Buf), "%s\n    {\"id\": %s, \"role\": \"%s\", \"legend\": %s, \"shift\": %s, \"front\": %s, \"w\": %g, \"h\": %g, \"row\": %d, \"x\": %.4f, \"y\": %.4f, \"topW\": %.4f, \"topD\": %.4f, \"topY\": %.4f, \"topZ\": %.4f, \"tilt\": %.5f, \"frontY\": %.4f, \"frontZ\": %.4f, \"frontTilt\": %.5f}",
			Keys.empty() ? "" : ",", Json(K.Id).c_str(), RoleName(K.Role), Json(K.Legend).c_str(), Json(K.ShiftLegend).c_str(), Json(K.FrontLegend).c_str(), K.W, K.H, K.ProfileRow, C.X, C.Y, T.Width, T.Depth, T.CenterY, T.CenterZ, T.TiltRad, T.FrontFaceY, T.FrontFaceZ, T.FrontTiltRad);
		Keys += Buf;
	}
	if (!Board.Write(Dir + "/board.obj"))
	{
		std::printf("cannot write %s/board.obj\n", Dir.c_str());
		return 1;
	}

	// Exploded switch: housing and pins, contacts, spring, top housing, stem, keycap, bottom to top.
	Groups Exploded;
	const double Lift[6] = {0.0, 9.0, 16.0, 30.0, 42.0, 0.0};
	AddSwitch(Exploded, Sw, Vec3(), SP.StemTopZ, Lift);
	mgeo::KeycapParams Cap;
	Cap.ProfileRow = 3;
	Exploded.Get("cap:alpha:K").Append(mgeo::BuildKeycap(Cap), Vec3(0, 0, SP.StemTopZ + 60.0));
	if (!Exploded.Write(Dir + "/switch_exploded.obj"))
	{
		return 1;
	}

	std::FILE* F = std::fopen((Dir + "/scene.json").c_str(), "w");
	if (!F)
	{
		return 1;
	}
	const Vec3& D = B.DisplayCenter;
	std::fprintf(F, "{\n  \"restTiltRad\": %.6f,\n  \"restLiftMm\": %.4f,\n  \"stemTopZ\": %.4f,\n  \"display\": {\"x\": %.4f, \"y\": %.4f, \"z\": %.4f, \"w\": %.4f, \"h\": %.4f, \"tilt\": %.6f},\n  \"keys\": [%s\n  ]\n}\n",
		B.RestTiltRad, B.RestLiftMm, SP.StemTopZ, D.X, D.Y, D.Z, B.DisplayWidthMm, B.DisplayHeightMm, B.DisplayTiltRad, Keys.c_str());
	std::fclose(F);

	if (!WriteMechanism(Dir, L, CP))
	{
		std::printf("cannot write the mechanism\n");
		return 1;
	}
	std::printf("board.obj: %d triangles in %zu groups; switch_exploded.obj: %d triangles\n", Board.Triangles(), Board.List.size(), Exploded.Triangles());
	std::printf("per part: keycap 1u R3 %d, switch %d (housing %d+%d, stem %d, spring %d, contacts %d, pins %d)\n", mgeo::BuildKeycap(Cap).TriangleCount(),
		Sw.HousingBottom.TriangleCount() + Sw.HousingTop.TriangleCount() + Sw.Stem.TriangleCount() + Sw.Spring.TriangleCount() + Sw.Contacts.TriangleCount() + Sw.Pins.TriangleCount(),
		Sw.HousingBottom.TriangleCount(), Sw.HousingTop.TriangleCount(), Sw.Stem.TriangleCount(), Sw.Spring.TriangleCount(), Sw.Contacts.TriangleCount(), Sw.Pins.TriangleCount());
	return 0;
}
