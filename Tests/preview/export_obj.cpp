// Writes the assembled keyboard for the preview renderer (Tests/preview/render.mjs):
//   board.obj            the board with keycaps on switches, one group per material and one per keycap
//   switch_exploded.obj  a single tactile switch and its keycap, pulled apart along z
//   scene.json           legends and colour roles per keycap, the rest pose and the display frame
// Usage: mcalc_export_obj [output directory]
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

	std::printf("board.obj: %d triangles in %zu groups; switch_exploded.obj: %d triangles\n", Board.Triangles(), Board.List.size(), Exploded.Triangles());
	std::printf("per part: keycap 1u R3 %d, switch %d (housing %d+%d, stem %d, spring %d, contacts %d, pins %d)\n", mgeo::BuildKeycap(Cap).TriangleCount(),
		Sw.HousingBottom.TriangleCount() + Sw.HousingTop.TriangleCount() + Sw.Stem.TriangleCount() + Sw.Spring.TriangleCount() + Sw.Contacts.TriangleCount() + Sw.Pins.TriangleCount(),
		Sw.HousingBottom.TriangleCount(), Sw.HousingTop.TriangleCount(), Sw.Stem.TriangleCount(), Sw.Spring.TriangleCount(), Sw.Contacts.TriangleCount(), Sw.Pins.TriangleCount());
	return 0;
}
