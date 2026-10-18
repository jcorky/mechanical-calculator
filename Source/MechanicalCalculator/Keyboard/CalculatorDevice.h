// The keyboard's firmware: turns key presses into an edited input line, runs it through the
// calculator Session, and lays out the character screen (the printed tape, or the OLED). Engine-independent.
//
// Screen (Cols x Rows characters):
//   row 0             status: angle mode, digits, CAPS / SHIFT / NUM / FN, library page
//   rows 1 .. Rows-2  the tape (inputs and results, newest at the bottom; PgUp/PgDn scroll it),
//                     or the library browser while one is open
//   row Rows-1        "> " and the input line, scrolled horizontally to keep the cursor in view
#pragma once

#include "CalcSession.h"
#include "KeyLayout.h"

#include <cstdint>
#include <string>
#include <vector>

namespace mkb
{
struct DeviceEvents
{
	bool EvaluationStarted = false;   // async mode: PendingExpression() is ready to run on the machine
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

	// Asynchronous evaluation, for the Unreal layer: Enter hands the line over (EvaluationStarted)
	// and the keyboard locks, as a calculating machine's does while its motor runs, until the
	// result comes back through CompleteEvaluation. Synchronous (the default) runs it at once.
	void SetAsync(bool bInAsync) { bAsync = bInAsync; }
	bool Busy() const { return bBusy; }
	const std::string& PendingExpression() const { return Pending; }
	void CompleteEvaluation(const mcalc::Outcome& O, const mcalc::GearStats& Work);
	// Prints the crank-turn count of each evaluation on the tape (on by default).
	void SetShowWork(bool bShow) { bShowWork = bShow; }
	const mcalc::GearStats& LastWork() const { return Work; }

	const std::string& Input() const { return Line; }
	int Cursor() const { return CursorPos; }
	const std::vector<std::string>& Tape() const { return TapeLines; }
	// Tape lines printed since power-on (Tape() keeps only the last kMaxTape); the printer feeds the
	// paper by the newest ones.
	uint64_t PrintedCount() const { return Printed; }
	bool IsShifted() const { return ShiftCount > 0; }
	bool CapsLock() const { return Caps; }
	bool NumLock() const { return Num; }
	bool FnHeld() const { return FnCount > 0; }
	bool Exploded() const { return Explode; }
	int BacklightMode() const { return Backlight; }
	std::string StatusText() const { return StatusLine(64); }
	bool LibraryOpen() const { return LibCategory >= 0; }
	mcalc::Session& GetSession() { return Calc; }
	const mcalc::Session& GetSession() const { return Calc; }

	static constexpr int kBacklightModes = 6;

private:
	mcalc::Session Calc;
	bool bAsync = false;
	bool bBusy = false;
	bool bShowWork = true;
	std::string Pending;
	mcalc::GearStats Work;
	std::string Line;
	int CursorPos = 0;
	std::vector<std::string> History;
	int HistoryIndex = -1;          // -1 = editing a new line
	std::string Draft;              // the unfinished line while browsing history
	std::vector<std::string> TapeLines;
	uint64_t Printed = 0;
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
