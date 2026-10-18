// The keyboard's firmware: turns key presses into an edited input line, runs it through the
// calculator Session, and lays out the OLED's character screen. Engine-independent.
//
// Screen (Cols x Rows characters):
//   row 0             status: angle mode, digits, CAPS / SHIFT / NUM / FN, library page
//   rows 1 .. Rows-2  the tape (inputs and results, newest at the bottom; PgUp/PgDn scroll it),
//                     or the library browser while one is open
//   row Rows-1        "> " and the input line, scrolled horizontally to keep the cursor in view
#pragma once

#include "CalcSession.h"
#include "KeyLayout.h"

#include <string>
#include <vector>

namespace mkb
{
struct DeviceEvents
{
	bool Evaluated = false;
	bool EvaluationFailed = false;
	bool ExplodeToggled = false;
	bool BacklightCycled = false;
	bool ScreenChanged = false;
};

class CalculatorDevice
{
public:
	CalculatorDevice();

	DeviceEvents Press(const KeyDef& Key);
	void Release(const KeyDef& Key);

	// Types text as if each character had its key pressed (used by tests and demos).
	DeviceEvents TypeText(const std::string& Text);
	DeviceEvents Enter();

	std::vector<std::string> Screen(int Cols, int Rows, bool CursorOn) const;

	const std::string& Input() const { return Line; }
	int Cursor() const { return CursorPos; }
	const std::vector<std::string>& Tape() const { return TapeLines; }
	bool IsShifted() const { return ShiftCount > 0; }
	bool CapsLock() const { return Caps; }
	bool NumLock() const { return Num; }
	bool FnHeld() const { return FnCount > 0; }
	bool Exploded() const { return Explode; }
	int BacklightMode() const { return Backlight; }
	bool LibraryOpen() const { return LibCategory >= 0; }
	mcalc::Session& GetSession() { return Calc; }
	const mcalc::Session& GetSession() const { return Calc; }

	static constexpr int kBacklightModes = 6;

private:
	mcalc::Session Calc;
	std::string Line;
	int CursorPos = 0;
	std::vector<std::string> History;
	int HistoryIndex = -1;          // -1 = editing a new line
	std::string Draft;              // the unfinished line while browsing history
	std::vector<std::string> TapeLines;
	int TapeScroll = 0;             // lines scrolled up from the bottom
	int ShiftCount = 0;
	int CtrlCount = 0;
	int FnCount = 0;
	bool Caps = false;
	bool Num = true;
	bool Explode = false;
	int Backlight = 2;
	int LibCategory = -1;           // index into the registry's categories, -1 = closed
	int LibSelection = 0;
	std::vector<std::string> Completions;
	int CompletionIndex = -1;
	size_t CompletionStart = 0;

	DeviceEvents Act(KeyAction Action, const KeyDef* Key);
	void InsertText(const std::string& Text);
	void Backspace(bool Word);
	void Complete();
	void Help();
	void Evaluate(DeviceEvents& Ev);
	void AddTape(const std::string& Text, const char* Prefix);
	std::string WordBeforeCursor(size_t* Start) const;
	KeyAction NumpadNavigation(const KeyDef& Key) const;
	std::string StatusLine(int Cols) const;
};

const char* BacklightModeName(int Mode);
} // namespace mkb
