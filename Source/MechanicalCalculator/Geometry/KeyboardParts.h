// Every part of the keyboard as a procedural mesh: sculpted keycaps, switches (housing, stem,
// spring, contact leaves), stabilizers, plate, PCB, case, feet and the OLED module.
// Same frame and winding as MeshBuilder.h (millimetres, x right, y away from the typist, z up).
//
// Stack-up (z = 0 is the top face of the plate):
//   case floor ............ z = -CaseParams::DepthBelowPlate
//   PCB ................... top at z = -PcbGapBelowPlate (switch pins pass through it)
//   plate ................. z in [-PlateThickness, 0]
//   switch housing top .... rests on the plate; the stem's top (cap mount) at SwitchParams::StemTopZ when up
//   keycap ................ its local origin sits on the stem top and moves with it
#pragma once

#include "MeshBuilder.h"
#include "KeyLayout.h"

#include <vector>

namespace mgeo
{
struct KeycapParams
{
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

Mesh BuildKeycap(const KeycapParams& P);
KeycapTop KeycapTopFor(const KeycapParams& P);
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
	Vec3 DisplayCenter;         // centre of the glass, top surface
	double DisplayWidthMm = 0.0;
	double DisplayHeightMm = 0.0;
	double DisplayTiltRad = 0.0;
	// To set the board on a desk: rotate about x by RestTiltRad (positive lifts the back), then
	// raise by RestLiftMm so the lowest point of the feet touches z = 0.
	double RestTiltRad = 0.0;
	double RestLiftMm = 0.0;
};

BoardParts BuildBoard(const mkb::Layout& Layout, const CaseParams& P);

// Centre of a key's switch on the plate surface (z = 0), from its layout position:
// x = (CenterX - Width/2) * unit, y = (Height/2 - CenterY) * unit, so the board is centred on the origin.
Vec3 KeyCenter(const mkb::KeyDef& Key, const mkb::Layout& Layout);
} // namespace mgeo
