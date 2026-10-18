// The physical keyboard: a full-size ANSI 104-key board with an OLED module above the numpad.
// Positions are in key units (1 u = 19.05 mm), measured from the top-left corner of the Esc key,
// x to the right and y toward the typist. Engine-independent.
#pragma once

#include <string>
#include <vector>

namespace mkb
{
constexpr double kUnitMm = 19.05;

enum class KeyRole
{
	Alpha,      // letters, digits, punctuation
	Modifier,   // shift, ctrl, alt, win, fn, caps, tab, backspace
	Accent,     // Esc and the main Enter: contrasting caps
	Function,   // F1-F12 and the three above the navigation cluster
	Navigation, // arrows, Ins/Del/Home/End/PgUp/PgDn
	Numpad,
	Space,
};

enum class KeyAction
{
	Insert,       // types Insert (or ShiftInsert while shift is held)
	Evaluate,     // Enter, numpad Enter
	Backspace,
	DeleteForward,
	ClearLine,    // Esc
	Complete,     // Tab: complete the name before the cursor
	CursorLeft,
	CursorRight,
	Home,
	End,
	HistoryUp,
	HistoryDown,
	PageUp,       // scroll the result pane
	PageDown,
	Shift,
	CapsLock,
	Control,
	Alt,
	Super,
	Fn,
	Menu,         // cycles the library page shown on the display
	Help,         // F1: help for the name at the cursor
	AngleMode,    // F2: DEG/RAD
	Digits,       // F3: cycle display precision
	Explode,      // Pause: exploded view of the build
	Backlight,    // ScrLk: cycle the RGB backlight
	CopyResult,   // PrtSc: insert the last result's value at the cursor
	NumLock,
	None,
};

struct KeyDef
{
	std::string Id;          // unique, e.g. "A", "Num7", "LShift"
	std::string Legend;      // main legend printed on the top
	std::string ShiftLegend; // second legend (number row symbols), empty if none
	std::string FrontLegend; // side-printed calculator function, empty if none
	std::string EngineKey;   // Unreal FKey name, empty if the engine has no such key
	double X = 0.0;          // top-left, key units
	double Y = 0.0;
	double W = 1.0;
	double H = 1.0;
	int ProfileRow = 3;      // sculpted profile row, 1 (tallest, top) .. 4
	KeyRole Role = KeyRole::Alpha;
	KeyAction Action = KeyAction::Insert;
	std::string Insert;      // text typed (Action::Insert, or the F-key macros)
	std::string ShiftInsert; // text typed with shift held (empty = same as Insert)
	bool NeedsStabilizer() const { return W >= 2.0 || H >= 2.0; }
	double CenterX() const { return X + W / 2.0; }
	double CenterY() const { return Y + H / 2.0; }
};

struct Rect
{
	double X = 0.0, Y = 0.0, W = 0.0, H = 0.0;   // key units, same frame as KeyDef
};

struct Layout
{
	std::vector<KeyDef> Keys;
	double Width = 0.0;     // overall key field, units
	double Height = 0.0;
	Rect Display;           // the OLED module's window, units

	const KeyDef* Find(const std::string& Id) const;
	const KeyDef* FindByEngineKey(const std::string& EngineKey) const;
};

// The board this project builds.
const Layout& FullSizeAnsi();
} // namespace mkb
