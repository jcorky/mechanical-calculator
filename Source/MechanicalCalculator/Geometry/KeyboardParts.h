// Every part of the keyboard as a procedural mesh: sculpted keycaps, switches (housing, stem,
// spring, contact leaves), stabilizers, plate, PCB, case, feet and the OLED module.
// Same frame and winding as MeshBuilder.h (millimetres, x right, y away from the typist, z up).
//
// Stack-up (z = 0 is the top face of the plate):
//   case floor ............ z = -CaseParams::DepthBelowPlate
//   PCB ................... top at z = -(PlateThickness + PcbGapBelowPlate), i.e. 5.0 mm below the
//                           plate top as on every MX board; the switch pins pass through it
//   plate ................. z in [-PlateThickness, 0]
//   switch housing top .... rests on the plate; the stem's top (cap mount) at SwitchParams::StemTopZ when up
//   keycap ................ its local origin sits on the stem top and moves with it
//
// "Stem top" is the cap mount: the flat top of the MX slider, from which the + cross rises a further
// 3.7 mm into the keycap's socket. With the key up it is 6.6 mm above the plate (15.2 mm above the
// PCB to the tip of the cross, per the Cherry MX spec, minus the cross).
#pragma once

#include "MeshBuilder.h"
#include "KeyLayout.h"

#include <vector>

namespace mgeo
{
enum class KeycapStyle
{
	Sculpted,   // Cherry-profile sculpted caps
	Vintage,    // round, flat-topped typewriter keys with a chrome rim ring, on a stem
};

struct KeycapParams
{
	KeycapStyle Style = KeycapStyle::Sculpted;
	double WidthU = 1.0;        // key units
	double DepthU = 1.0;
	int ProfileRow = 3;         // Cherry-style sculpt, 1 (number row, tallest) .. 4 (bottom rows)
	double UnitMm = mkb::kUnitMm;
	double GapMm = 0.8;         // footprint is unit * size - gap
	double TopInsetMm = 2.8;    // the top is this much smaller than the footprint on each side
	double CornerRadiusMm = 1.2;
	double TopCornerRadiusMm = 2.0;
	double DishDepthMm = 0.6;   // cylindrical dish along x
	double WallMm = 1.2;        // shell thickness (caps are hollow, with a cross-stem socket)
	bool Spacebar = false;      // flat (no dish) and convex top
	bool HomingBar = false;     // the F and J bumps
	int Segments = 6;           // corner subdivisions
};

// Where a legend goes: the dish centre, its tilt about x, and the printable size.
struct KeycapTop
{
	double CenterZ = 0.0;       // above the cap's local origin
	double CenterY = 0.0;
	double TiltRad = 0.0;       // positive tilts the top toward the typist
	double Width = 0.0;
	double Depth = 0.0;
	double FrontFaceZ = 0.0;    // centre of the front wall, for side-printed legends
	double FrontFaceY = 0.0;
	double FrontTiltRad = 0.0;
};

// Cap-local frame: the origin is the centre of the mouth of the cap's centre stem socket, which rests
// on the stem top (the cap mount), so the whole cap moves with SwitchParts::Stem. The skirt reaches
// 0.6 mm below the origin (the stem inset of Cherry-profile caps); the + cut-out runs up from z = 0.
// Caps 2u and wider (or 2u tall) also carry stabilizer sockets at +-StabilizerSpacingMm.
// UVs: the top surface maps the dish to 0..1 (U along x, V from front to back); every other vertex
// has U in [-2, -1.5], so a legend material can mask on U >= 0.
// KeycapTopFor returns values in this same cap-local frame.
// Vintage caps: a round top (a stadium for wider keys) on a stem column, with a slightly dished
// legend face that carries the same 0..1 UVs; their chrome rim ring is a separate mesh.
Mesh BuildKeycap(const KeycapParams& P);
Mesh BuildKeycapRim(const KeycapParams& P);   // the vintage cap's chrome ring (empty for sculpted caps)
KeycapTop KeycapTopFor(const KeycapParams& P);
// Sculpted row height: the top's plane at its back edge, above the bottom of the skirt.
double KeycapHeightMm(int ProfileRow);

enum class SwitchType
{
	Linear,
	Tactile,
	Clicky,
};

struct SwitchParams
{
	SwitchType Type = SwitchType::Tactile;
	double TravelMm = 4.0;
	double ActuationMm = 2.0;
	double StemTopZ = 6.6;      // stem top above the plate surface with the key up
	int SpringTurns = 18;
	int Segments = 12;
};

// Frames: HousingBottom, HousingTop, Contacts, Pins and Spring are in the plate frame (z = 0 is the
// plate's top surface), centred on the switch. Stem has its own frame: z = 0 is its top (the cap
// mount) with the key up, and the + cross rises above it; the keycap's local origin sits at that top,
// so the adapter places both at StemTopZ and moves them down together by the travel.
struct SwitchParts
{
	Mesh HousingBottom;         // opaque base below the plate, with the centre post and pins
	Mesh HousingTop;            // translucent top above the plate
	Mesh Stem;                  // moves with the key; origin such that z = 0 is its top when up
	Mesh Spring;                // helix at rest length (the adapter scales it in z as the key travels)
	Mesh Contacts;              // the two metal leaves
	Mesh Pins;                  // the two pins through the PCB
	double SpringRestTopZ = 0.0;   // spring top above the plate with the key up
	double SpringBottomZ = 0.0;
};

SwitchParts BuildSwitch(const SwitchParams& P);

// Costar-style plate-mount stabilizer for a key W units wide (spacing from the standard table).
// Relative to the key's centre: Housings and Wire are in the plate frame (z = 0 is the plate's top
// surface, x/y are offsets from the key centre). Stems is in the same frame as SwitchParts::Stem
// (z = 0 is the stem top with the key up), so it attaches to the moving key. For Vertical = true the
// two stems are spaced along y instead of x.
struct StabilizerParts
{
	Mesh Housings;              // the two inserts clipped into the plate
	Mesh Stems;                 // move with the key
	Mesh Wire;                  // the bent wire
};
StabilizerParts BuildStabilizer(double WidthU, bool Vertical);
double StabilizerSpacingMm(double WidthU);

struct CaseParams
{
	double MarginMm = 12.0;          // case wall outside the key field
	double DepthBelowPlate = 18.0;
	double LipAbovePlate = 4.0;      // how far the wall rises above the plate at the front
	double TypingAngleDeg = 6.0;     // the case top tilts down toward the typist by this angle
	double WallMm = 4.0;
	double CornerRadiusMm = 6.0;
	double ChamferMm = 1.5;          // bevel on the top outer edge
	double PlateThickness = 1.5;
	double PcbThickness = 1.6;
	double PcbGapBelowPlate = 3.5;
	bool OledModule = false;         // the OLED above the numpad; without it a brass name plate sits there
};

struct BoardParts
{
	Mesh CaseShell;             // the tray: walls, floor, chamfer, inner pockets
	Mesh CaseBadge;             // a small inset plate on the back wall
	Mesh Feet;                  // rubber feet under the case
	Mesh Plate;                 // with 14 mm switch cut-outs and stabilizer slots
	Mesh Pcb;                   // green board with pads under every switch
	Mesh PcbParts;              // diodes, controller chip, USB-C receptacle
	Mesh DisplayBezel;          // the OLED module's frame
	Mesh DisplayGlass;          // the emissive screen surface (UVs span 0..1 for a render target)
	// The glass's top surface faces +z, rotated about +x by DisplayTiltRad (right-handed): a positive
	// tilt raises its top edge (the edge away from the typist) so the screen leans back toward the
	// typist. UV (0, 0) is the top-left corner as the typist sees it (min x, max y), U grows to the
	// right and V toward the typist.
	Vec3 DisplayCenter;         // centre of the glass, top surface
	double DisplayWidthMm = 0.0;
	double DisplayHeightMm = 0.0;
	double DisplayTiltRad = 0.0;
	// To set the board on a desk: rotate about x by RestTiltRad (positive lifts the back), then
	// raise by RestLiftMm so the lowest point of the feet touches z = 0.
	double RestTiltRad = 0.0;
	double RestLiftMm = 0.0;
	Mesh NamePlate;             // brass plate in the display window when there is no OLED module
	Vec3 NamePlateCenter;       // centre of its sunk field (faces +z), for text
	double NamePlateWidthMm = 0.0;
	double NamePlateHeightMm = 0.0;
};

BoardParts BuildBoard(const mkb::Layout& Layout, const CaseParams& P);

// Centre of a key's switch on the plate surface (z = 0), from its layout position:
// x = (CenterX - Width/2) * unit, y = (Height/2 - CenterY) * unit, so the board is centred on the origin.
Vec3 KeyCenter(const mkb::KeyDef& Key, const mkb::Layout& Layout);
} // namespace mgeo
