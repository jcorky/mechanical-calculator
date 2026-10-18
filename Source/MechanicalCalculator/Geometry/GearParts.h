// The calculating mechanism: a Leibniz / Thomas-arithmometer stepped-drum machine on a black plinth
// behind the keyboard. It does all the arithmetic; a type-wheel printer on its left prints to a paper
// tape, and a bank of memory columns behind it shows `ans` and the variables.
// Same units and winding as MeshBuilder.h (millimetres, counter-clockwise outward triangles).
//
// Primitives are built around their own axis (local +Z) at the origin. BuildMechanism places them.
//
// Frames. Every pose in MechanismParts is in the DESK frame unless marked "carriage frame":
//   desk frame ...... x to the right, y away from the typist, z up; x = 0 and y = 0 under the keyboard's
//                     origin. z = 0 is the top of the plinth, which is where BoardParts' rest transform
//                     puts the keyboard's feet (rotate about x by RestTiltRad, then raise by
//                     RestLiftMm). So the keyboard and the mechanism share one parent: give the
//                     mechanism an identity transform in the keyboard's parent and the board its rest
//                     transform. The plinth fills z in [-PlinthHeightMm, 0]; the real desk is below it.
//   carriage frame .. the desk frame translated by Carriage.Poses[0] at shift 0. The carriage, its
//                     wheels, carry levers, decimal markers and clearing levers are given in it. At
//                     shift k the carriage moves along +x by k * DigitPitchMm.
//
// Animation contract. A Pose maps a part's local coordinates to its parent: P = Origin + X*l.x +
// Y*l.y + Z*l.z (X, Y, Z orthonormal and right-handed). Every rotating part turns about its LOCAL +Z
// axis (right-handed, positive angles counter-clockwise seen from +Z): apply R_z(angle) to the local
// coordinates, then the pose. Translating parts move along the local axis named in their comment.
// Each part is drawn with all of its Sections at the same transform.
//
// Numeral wheels (setting, counter, accumulator, memory figure wheels, type wheels): turning a wheel
// by DigitAngleRad(d) = d * 36 degrees about its local +Z brings digit d to its local +Y, which each
// pose points at that wheel's reading window (or, for type wheels, at the platen).
#pragma once

#include "KeyboardParts.h"
#include "MeshBuilder.h"
#include "KeyLayout.h"

#include <string>
#include <vector>

namespace mgeo
{
// --- Primitives -----------------------------------------------------------------------------------

struct SpurGearParams
{
	int Teeth = 20;
	double Module = 1.0;
	double PressureAngleDeg = 20.0;
	double FaceWidth = 4.0;       // along z, centred on z = 0
	double BoreRadius = 0.0;      // 0 for a solid gear
	bool SquareBore = false;      // a square bore of side 2 * BoreRadius (slides on a square shaft)
	double HubRadius = 0.0;       // a hub on the +z face, 0 for none
	double HubLength = 0.0;
	int FlankPoints = 5;          // involute samples per flank
};

// Pitch, base, tip and root radii of a standard full-depth gear (addendum m, dedendum 1.25 m).
struct GearCircles
{
	double Pitch = 0.0, Base = 0.0, Tip = 0.0, Root = 0.0;
};
GearCircles SpurGearCircles(const SpurGearParams& P);

// A true involute spur gear: axis +Z, centred on the origin, tooth 0 centred on +X. Tooth thickness on
// the pitch circle is pi * m / 2; below the base circle the flanks run radially to the root.
Mesh BuildSpurGear(const SpurGearParams& P);

// A spoked gear: an involute rim of Teeth teeth on a ring, Spokes straight spokes and a bored hub.
// Axis +Z, centred on the origin. With Teeth = 0 it is a plain spoked wheel (a flywheel of radius
// RimRadius).
struct SpokedWheelParams
{
	int Teeth = 60;
	double Module = 1.0;
	double RimRadius = 0.0;       // outer radius when Teeth = 0
	double RimDepth = 4.0;        // radial thickness of the ring under the teeth
	double FaceWidth = 6.0;
	int Spokes = 6;
	double SpokeWidth = 4.0;
	double HubRadius = 7.0;
	double HubLength = 10.0;
	double BoreRadius = 2.5;
};
Mesh BuildSpokedWheel(const SpokedWheelParams& P);

struct BevelGearParams
{
	int Teeth = 16;
	double Module = 1.0;          // at the heel (outer end)
	double PitchAngleDeg = 45.0;  // pitch cone half-angle; 45 for a 1:1 mitre pair
	double FaceWidth = 3.0;       // along the cone
	double BoreRadius = 0.0;
	double HubRadius = 0.0;       // a hub behind the heel, 0 for none
	double HubLength = 0.0;
	double PressureAngleDeg = 20.0;
	int FlankPoints = 4;
};
// Axis +Z with the pitch-cone apex at the origin; the gear lies at z > 0 with its heel pitch circle
// (radius m * z / 2) at z = PitchRadius / tan(PitchAngle). Two bevels mesh when their apexes coincide.
Mesh BuildBevelGear(const BevelGearParams& P);

// A Leibniz stepped drum: a core with Teeth axial teeth of stepped lengths (tooth k = 1..Teeth is
// k * StepMm long) on consecutive positions of a ToothPositions-tooth involute gear of the given module
// (so it meshes with a pinion of the same module at centre distance m * (ToothPositions + z) / 2).
// All teeth end flush at the far end, so a pinion at digit d engages exactly d teeth per turn.
struct StepDrumParams
{
	double Module = 0.6;
	int ToothPositions = 20;
	int Teeth = 9;
	double StepMm = 3.14159265358979323846;   // pinion travel per digit
	double PinionFaceWidth = 2.2;              // must be < StepMm
	double StartZ = 6.0;                       // the drum body starts here (the bevel sits in front)
	double AxleRadius = 1.2;
	bool Bevel = true;                         // a 45-degree drive bevel with its apex at the origin
	int BevelTeeth = 16;
};
struct StepDrumLayout
{
	double PitchRadius = 0.0, RootRadius = 0.0, TipRadius = 0.0;
	double BodyStartZ = 0.0, BodyEndZ = 0.0;   // the core
	double TeethEndZ = 0.0;                    // the flush end of every tooth
	double ToothStartZ(int K) const;           // tooth K (1 .. 9) spans [ToothStartZ(K), TeethEndZ]
	double ToothAngleRad(int K) const;         // tooth K's centre angle about +Z
	double PinionCenterZ(int Digit) const;     // a pinion here engages Digit teeth
	double StepMm = 0.0;
	int ToothPositions = 20;
};
StepDrumLayout StepDrumGeometry(const StepDrumParams& P);
// Axis +Z, origin at the drive bevel's apex (the drum's pivot); the body runs from StartZ along +Z.
Mesh BuildStepDrum(const StepDrumParams& P);

// A numeral wheel: a rim with the digits 0-9 embossed around it (a separate mesh so the paint can be a
// second material), optional drive gear, carry star and carry tooth. Axis +Z, rim centred on z = 0.
// Digit d's glyph is centred at angle 90 - 36 d degrees about +Z, so turning the wheel by
// DigitAngleRad(d) brings digit d to local +Y, the side that faces the reading window.
struct NumeralWheelParams
{
	double Radius = 6.2;            // rim radius (glyphs stand Emboss proud of it)
	double Width = 8.0;             // rim width along the axis
	double GlyphHeight = 4.4;
	double GlyphWidth = 2.6;
	double StrokeMm = 0.5;
	double EmbossMm = 0.3;
	bool GlyphsAlongAxis = true;    // true: upright along +Z (read across the axis); false: odometer style
	bool SignWheel = false;         // glyphs '+' at even and '-' at odd digit positions
	double BoreRadius = 0.8;
	int GearTeeth = 0;              // a drive gear on the -z face (0: none)
	double GearModule = 0.5;
	double GearGap = 0.0;           // spacing between the rim and the gear (bridged by a hub collar)
	int StarTeeth = 0;              // a carry / counting star on the +z face (0: none)
	double StarModule = 0.6;
	bool StarOnFront = false;       // put the star on the -z face instead
	bool CarryTooth = false;        // a single carry tooth on the face opposite the star (or on +z)
};
struct NumeralWheel
{
	Mesh Wheel;                     // the rim
	Mesh Gears;                     // the gear, star and carry tooth
	Mesh Glyphs;                    // the embossed digits
	std::vector<int> GlyphTriangleStart;   // digit d's triangles are [GlyphTriangleStart[d], [d + 1])
	double CarryToothAngleRad = 0.0;       // where the carry tooth sits about +Z
};
NumeralWheel BuildNumeralWheel(const NumeralWheelParams& P);
double DigitAngleRad(int Digit);   // Digit * 36 degrees
// Raised digits ('0'-'9', '+', '-'; anything else is a space) on the plane z = 0, facing +z, centred
// on the origin and reading along +x, from z = -0.1 to EmbossMm. Same stroke font as the wheels.
Mesh BuildDigitText(const std::string& Digits, double Height, double StrokeMm, double EmbossMm);

// A round (or square, side 2 * Radius) shaft along +Z from z = 0 to Length, with chamfered ends.
Mesh BuildShaft(double Radius, double Length, bool Square = false);
// A pillow block: base on z = 0, footprint centred on the origin, bore axis along +X at AxisHeight.
Mesh BuildBearingBlock(double BoreRadius, double AxisHeight, double Width);
// A straight rack along +X from x = 0, teeth toward +Y with the pitch line on y = 0, FaceWidth along
// z (centred), Body mm of solid below the root.
Mesh BuildRack(double Module, int Teeth, double FaceWidth, double Body, double PressureAngleDeg = 20.0);

// --- The mechanism --------------------------------------------------------------------------------

enum class MechMaterial
{
	Brass,        // polished: drums, gears, figure wheels, plates
	Steel,        // bright: shafts, bearing caps, rails, pins
	CastIron,     // black cast iron: the side frames and the base
	BlackEnamel,  // window bars, carriage frame, printer frame, plinth
	NumeralBlack, // the black numeral wheels
	Paint,        // white digits on the black wheels
	Engraving,    // dark digits on the brass figure and type wheels
	WhiteEnamel,  // the tens-carry levers
	Glass,
	Paper,
	Rubber,
	Ribbon,
	Wood,         // turned handles and knobs
};
const char* MechMaterialName(MechMaterial M);

struct MechSection
{
	MechMaterial Material = MechMaterial::Brass;
	Mesh Geometry;
};

// One unique part, built around its pivot: one closed mesh per material.
struct MechPart
{
	std::vector<MechSection> Sections;
	Mesh& Of(MechMaterial M);       // the section for M (created on first use; references stay valid)
	int TriangleCount() const;
};

struct Pose
{
	Vec3 Origin;
	Vec3 X = Vec3(1, 0, 0), Y = Vec3(0, 1, 0), Z = Vec3(0, 0, 1);
	Vec3 Point(const Vec3& Local) const { return Origin + X * Local.X + Y * Local.Y + Z * Local.Z; }
	Vec3 Direction(const Vec3& Local) const { return X * Local.X + Y * Local.Y + Z * Local.Z; }
	Pose Then(const Pose& Child) const;        // Child given in this pose's frame -> parent frame
	Pose RotatedZ(double Angle) const;         // after the part turns by Angle about its local +Z
	Pose Moved(const Vec3& LocalOffset) const; // after the part moves by LocalOffset (local axes)
};

// A unique part and its instances. With OnCarriage the poses are in the carriage frame.
struct PartSet
{
	MechPart Model;
	std::vector<Pose> Poses;
	bool OnCarriage = false;
};

// A register window: one slot of a window bar, through which one digit is read.
struct ReadingWindow
{
	Vec3 Center;                 // on the outer face of the bar
	Vec3 Normal;                 // toward the reader
	Vec3 Up;                     // the glyph's up direction as read
	double WidthMm = 0.0, HeightMm = 0.0;
};

// A meshing pair, for checks: spur pairs (parallel axes), bevel pairs (apexes coincide) and the rack.
struct GearPair
{
	std::string Name;
	Vec3 PointA, AxisA;          // a point on gear A's axis (the apex for a bevel) and its direction
	Vec3 PointB, AxisB;          // for the rack: a point on its pitch line and its direction
	double Module = 0.0;
	int TeethA = 0, TeethB = 0;
	enum Kind
	{
		Spur,
		Bevel,
		Rack,
	} Type = Spur;
	bool OnCarriageB = false;    // gear B's point is in the carriage frame (at shift 0)
};

// A frame on the paper for one printed line.
struct TapeFrame
{
	Vec3 Origin;                 // on the printed face, at the centre of the usable width
	Vec3 Up;                     // along the paper toward older lines
	Vec3 Normal;                 // out of the printed face
	Vec3 Right;                  // along the line
};

// A flat face for a text label (maker's plate, memory column labels): Unreal lays text on it.
struct LabelFace
{
	Vec3 Center, Normal, Up;
	double WidthMm = 0.0, HeightMm = 0.0;
};

struct MechanismParts
{
	// Registers and pitches (fixed).
	int SettingDigits = 20;
	int AccumulatorDigits = 40;
	int CounterDigits = 20;
	int MaxShift = 20;
	double DigitPitchMm = 14.0;      // drums, wheels and the carriage step
	double PinionStepMm = 0.0;       // pinion travel per set digit
	double DigitAngle = 0.0;         // 36 degrees: one digit on any numeral wheel
	double PlinthHeightMm = 0.0;

	// Static parts, in the desk frame (pose = identity).
	MechPart Plinth;                 // the long black base under everything, keyboard included
	MechPart Frame;                  // cast-iron side frames, base plate, bosses, bearing caps, bars,
	                                 // tie rods, rails, the transfer shaft, the setting window bar,
	                                 // the drive housing
	MechPart GlassCover;             // optional: a glass top over the drum bank (not shown by default)
	MechPart Printer;                // printer frame, ribbon, tear bar
	MechPart PaperStrip;             // the tape from the roll, round the platen, out and over the top
	MechPart Bell;                   // the overdraft bell's dome and stud
	MechPart MemoryFrame;            // memory columns: base, pillars, axles, top rail, label plates
	LabelFace MakersPlate;           // on the plinth's front edge
	std::vector<LabelFace> MemoryLabels;   // one per memory column

	// Moving parts. Rotations are about each pose's local +Z unless stated.
	PartSet Drums;                   // 20: +1 turn per crank turn (addition)
	PartSet SquareShafts;            // 20: turn with their pinions
	PartSet Pinions;                 // 20, at digit 0: move along local +Z by Digit * PinionStepMm;
	                                 // turn Digit / 12 of a turn per crank turn
	PartSet TransferGears;           // 20: on the shared transfer shaft, turn with their square shafts
	                                 // (1:1) and drive the accumulator wheel over their column 18:15
	PartSet SettingWheels;           // 20: turn to DigitAngleRad(set digit)
	PartSet SettingMarker;           // 1: decimal-point marker on the setting bar; at 0 it sits right of
	                                 // digit 0; move along local +X by Places * DigitPitchMm
	PartSet Carriage;                // 1: translates along local +X by Shift * DigitPitchMm; includes the
	                                 // counter and accumulator window bars
	PartSet AccumulatorWheels;       // 40, carriage frame: the wheel over column j turns +Digit_j / 10
	                                 // of a turn per crank turn (addition)
	PartSet CarryLevers;             // 40, carriage frame: lever i (under wheel i) swings by
	                                 // CarryLeverSwing when wheel i passes 9 -> 0 and moves wheel i + 1 on
	PartSet CounterWheels;           // 20, carriage frame: the wheel at the current shift turns +-36 deg
	PartSet CarriageMarkers;         // 2, carriage frame: [0] accumulator, [1] counter decimal markers;
	                                 // move along local +X by Places * DigitPitchMm
	PartSet ClearingLevers;          // 2, carriage frame: [0] accumulator, [1] counter
	PartSet CountingFinger;          // 1: one turn per crank turn; its tooth moves the counter on by one
	PartSet MainShaft;               // 1: +-1 turn per crank turn, with the 20 drum bevels
	PartSet ReverserIdler;           // 1: the bevel between the input and the main shaft
	PartSet InputShaft;              // 1: the crank shaft with the reverser input bevel
	PartSet DriveGear;               // 1: the big 6-spoked gear on the crank shaft (turns with it)
	PartSet Flywheel;                // 1: turns FlywheelRatio times the crank, the opposite way
	PartSet Crank;                   // 1: the crank arm, turning with the crank shaft
	PartSet CrankHandle;             // 1: the turned wooden handle, same pose and rotation as the crank
	                                 // (it may also spin about CrankHandleAxis, in the crank's frame)
	PartSet ReversingLever;          // 1: turns to ReverseAddAngle (add) or ReverseSubtractAngle
	PartSet ShiftKnob;               // 1: knob, shaft and bevel; ShiftKnobAnglePerStep per carriage step
	PartSet ShiftPinion;             // 1: bevel and pinion under the carriage rack
	PartSet BellHammer;              // 1: swings by BellStrikeAngle on overdraft
	PartSet TypeWheels;              // 18: printer type wheels; turn to DigitAngleRad(d) to print d
	PartSet PaperRoll;               // 1: turns by LinePitch / PaperRollRadius per printed line
	PartSet Platen;                  // 1: turns by LinePitch / PlatenRadius per printed line
	PartSet RibbonSpools;            // 2: turn as the ribbon advances
	PartSet MemoryWheels;            // 6 * 21: brass figure wheels (18 mantissa + 3 exponent per column)
	PartSet MemorySignWheels;        // 6: sign wheels; '+' at even, '-' at odd digit angles
	int MemoryColumns = 6, MemoryMantissa = 18, MemoryExponent = 3;
	// MemoryWheels.Poses[c * 21 + k]: column c (0 = left), k = 0..17 mantissa from the most significant
	// (top) down, k = 18..20 exponent from the most significant down. MemorySignWheels.Poses[c].

	// Kinematics.
	double CarryLeverSwing = 0.0;
	double ReverseAddAngle = 0.0, ReverseSubtractAngle = 0.0;
	double ShiftKnobAnglePerStep = 0.0;
	double BellStrikeAngle = 0.0;
	double ClearingLeverSwing = 0.0;
	double FlywheelRatio = 0.0;      // flywheel turns per crank turn (opposite sense)
	double PaperRollRadius = 0.0, PlatenRadius = 0.0;
	Vec3 CrankHandleAxis;            // the handle's spin axis through its centre, in the crank's frame
	Vec3 CrankHandleCenter;
	// Per crank turn (addition; subtraction negates everything after the reverser): crank, input shaft
	// and drive gear +1 turn; main shaft, drums and counting finger +1 turn; pinion, square shaft and
	// transfer gear of column j +Digit_j / 12 turn; the accumulator wheel over column j +Digit_j / 10.

	// Where each register is read (setting: desk frame; accumulator and counter: carriage frame).
	std::vector<ReadingWindow> SettingWindows, AccumulatorWindows, CounterWindows;
	std::vector<ReadingWindow> MemoryWindows;   // one per memory wheel (same order), then the signs

	// The printing tape: line pitch, usable width, and the paper's centre line from the print line
	// (index 0, the newest line) out over the top; see TapeLineFrame.
	double TapeWidthMm = 0.0, TapeLinePitchMm = 0.0, TapeUsableWidthMm = 0.0;
	std::vector<Vec3> TapePath;      // on the printed face, at the tape's centre
	std::vector<Vec3> TapeNormals;   // out of the printed face

	std::vector<GearPair> GearPairs;
	Vec3 BoundsMin, BoundsMax;       // the mechanism without the plinth and the overhanging carriage
};

MechanismParts BuildMechanism(const mkb::Layout& Layout, const CaseParams& Case);

// The frame for the text line DistanceMm along the tape from the print line (line n of the history
// sits at n * TapeLinePitchMm). Past the end of the tape the last frame is extended straight.
TapeFrame TapeLineFrame(const MechanismParts& M, double DistanceMm);
} // namespace mgeo
