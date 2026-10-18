#include "CalculatorDevice.h"

#include <algorithm>
#include <cctype>

namespace mkb
{
namespace device_detail
{
bool IsNameChar(char Ch)
{
	return std::isalnum(static_cast<unsigned char>(Ch)) != 0 || Ch == '_';
}

std::string Fit(const std::string& S, int Cols)
{
	const size_t W = static_cast<size_t>(std::max(Cols, 0));
	return S.size() > W ? S.substr(0, W) : S + std::string(W - S.size(), ' ');
}

std::string CommonPrefix(const std::vector<std::string>& Names)
{
	if (Names.empty())
	{
		return std::string();
	}
	std::string P = Names[0];
	for (const std::string& N : Names)
	{
		size_t i = 0;
		while (i < P.size() && i < N.size() && P[i] == N[i])
		{
			++i;
		}
		P.resize(i);
	}
	return P;
}

const int kDigitSteps[] = {4, 6, 8, 10, 12, 15};
constexpr size_t kMaxTape = 500;
} // namespace device_detail

using namespace device_detail;

const char* BacklightModeName(int Mode)
{
	static const char* Names[CalculatorDevice::kBacklightModes] = {"off", "static", "wave", "breathe", "reactive", "roles"};
	return Mode >= 0 && Mode < CalculatorDevice::kBacklightModes ? Names[Mode] : "?";
}

CalculatorDevice::CalculatorDevice()
{
	AddTape("Mechanical Calculator. F1 help, F4 library.", "");
}

void CalculatorDevice::AddTape(const std::string& Text, const char* Prefix)
{
	TapeLines.push_back(std::string(Prefix) + Text);
	if (TapeLines.size() > kMaxTape)
	{
		TapeLines.erase(TapeLines.begin(), TapeLines.begin() + static_cast<std::ptrdiff_t>(TapeLines.size() - kMaxTape));
	}
	TapeScroll = 0;
}

void CalculatorDevice::InsertText(const std::string& Text)
{
	Line.insert(static_cast<size_t>(CursorPos), Text);
	CursorPos += static_cast<int>(Text.size());
	HistoryIndex = -1;
}

void CalculatorDevice::Backspace(bool Word)
{
	if (CursorPos <= 0)
	{
		return;
	}
	int Start = CursorPos - 1;
	if (Word)
	{
		while (Start > 0 && Line[static_cast<size_t>(Start)] == ' ')
		{
			--Start;
		}
		while (Start > 0 && IsNameChar(Line[static_cast<size_t>(Start - 1)]) && IsNameChar(Line[static_cast<size_t>(Start)]))
		{
			--Start;
		}
	}
	Line.erase(static_cast<size_t>(Start), static_cast<size_t>(CursorPos - Start));
	CursorPos = Start;
}

std::string CalculatorDevice::WordBeforeCursor(size_t* Start) const
{
	size_t B = static_cast<size_t>(CursorPos);
	while (B > 0 && IsNameChar(Line[B - 1]))
	{
		--B;
	}
	*Start = B;
	return Line.substr(B, static_cast<size_t>(CursorPos) - B);
}

void CalculatorDevice::Complete()
{
	if (CompletionIndex >= 0 && !Completions.empty())
	{
		// Repeated Tab: cycle through the candidates.
		CompletionIndex = (CompletionIndex + 1) % static_cast<int>(Completions.size());
		Line.erase(CompletionStart, static_cast<size_t>(CursorPos) - CompletionStart);
		CursorPos = static_cast<int>(CompletionStart);
		InsertText(Completions[static_cast<size_t>(CompletionIndex)]);
		return;
	}
	size_t Start = 0;
	const std::string Prefix = WordBeforeCursor(&Start);
	if (Prefix.empty() || std::isdigit(static_cast<unsigned char>(Prefix[0])) != 0)
	{
		return;
	}
	const std::vector<std::string> Found = Calc.Complete(Prefix);
	if (Found.empty())
	{
		return;
	}
	if (Found.size() == 1)
	{
		const bool IsFunction = Calc.Functions().Find(Found[0]) != nullptr;
		InsertText(Found[0].substr(Prefix.size()) + (IsFunction ? "(" : ""));
		return;
	}
	const std::string Common = CommonPrefix(Found);
	if (Common.size() > Prefix.size())
	{
		InsertText(Common.substr(Prefix.size()));
		return;
	}
	std::string Shown;
	for (const std::string& N : Found)
	{
		Shown += (Shown.empty() ? "" : " ") + N;
	}
	AddTape(Shown, "  ");
	Completions = Found;
	CompletionStart = Start;
	// The next Tab starts cycling from the first candidate.
	CompletionIndex = static_cast<int>(Found.size()) - 1;
}

void CalculatorDevice::Help()
{
	// The name under the cursor, else the function whose call encloses it: "hohmann(6678e3, |".
	size_t B = static_cast<size_t>(CursorPos), E = B;
	while (B > 0 && IsNameChar(Line[B - 1]))
	{
		--B;
	}
	while (E < Line.size() && IsNameChar(Line[E]))
	{
		++E;
	}
	std::string Word = Line.substr(B, E - B);
	if (Word.empty() || std::isalpha(static_cast<unsigned char>(Word[0])) == 0)
	{
		Word.clear();
		int Depth = 0;
		for (size_t i = static_cast<size_t>(CursorPos); i-- > 0;)
		{
			const char Ch = Line[i];
			if (Ch == ')' || Ch == ']')
			{
				++Depth;
			}
			else if (Ch == '[' && Depth > 0)
			{
				--Depth;
			}
			else if (Ch == '(')
			{
				if (Depth > 0)
				{
					--Depth;
					continue;
				}
				size_t S = i;
				while (S > 0 && IsNameChar(Line[S - 1]))
				{
					--S;
				}
				Word = Line.substr(S, i - S);
				break;
			}
		}
	}
	const mcalc::Outcome O = Calc.Execute(Word.empty() ? "help" : "help " + Word);
	for (const std::string& L : mcalc::WrapText(O.Text, 0))
	{
		AddTape(L, O.Ok ? "  " : "! ");
	}
}

void CalculatorDevice::Evaluate(DeviceEvents& Ev)
{
	if (Line.empty())
	{
		return;
	}
	History.push_back(Line);
	HistoryIndex = -1;
	AddTape(Line, "> ");
	const mcalc::Outcome O = Calc.Execute(Line);
	Ev.Evaluated = true;
	Ev.EvaluationFailed = !O.Ok;
	if (!O.Ok)
	{
		AddTape(O.Text, "! ");
	}
	else if (O.HasValue)
	{
		for (size_t i = 0; i < O.Lines.size(); ++i)
		{
			AddTape(O.Lines[i], i == 0 ? "= " : "  ");
		}
	}
	else
	{
		for (const std::string& L : mcalc::WrapText(O.Text, 0))
		{
			if (!L.empty())
			{
				AddTape(L, "  ");
			}
		}
	}
	Line.clear();
	CursorPos = 0;
}

KeyAction CalculatorDevice::NumpadNavigation(const KeyDef& Key) const
{
	if (Key.Id == "Num8") return KeyAction::HistoryUp;
	if (Key.Id == "Num2") return KeyAction::HistoryDown;
	if (Key.Id == "Num4") return KeyAction::CursorLeft;
	if (Key.Id == "Num6") return KeyAction::CursorRight;
	if (Key.Id == "Num7") return KeyAction::Home;
	if (Key.Id == "Num1") return KeyAction::End;
	if (Key.Id == "Num9") return KeyAction::PageUp;
	if (Key.Id == "Num3") return KeyAction::PageDown;
	if (Key.Id == "NumDecimal") return KeyAction::DeleteForward;
	if (Key.Id == "Num0" || Key.Id == "Num5") return KeyAction::None;
	return Key.Action;
}

DeviceEvents CalculatorDevice::Press(const KeyDef& Key)
{
	KeyAction Action = Key.Action;
	if (Key.Role == KeyRole::Numpad && !Num)
	{
		Action = NumpadNavigation(Key);
	}
	if (FnCount > 0)
	{
		if (Action == KeyAction::HistoryUp) Action = KeyAction::PageUp;
		else if (Action == KeyAction::HistoryDown) Action = KeyAction::PageDown;
		else if (Action == KeyAction::CursorLeft) Action = KeyAction::Home;
		else if (Action == KeyAction::CursorRight) Action = KeyAction::End;
	}
	return Act(Action, &Key);
}

DeviceEvents CalculatorDevice::TypeText(const std::string& Text)
{
	DeviceEvents Ev;
	LibCategory = -1;
	InsertText(Text);
	Ev.ScreenChanged = true;
	return Ev;
}

DeviceEvents CalculatorDevice::Enter()
{
	return Act(KeyAction::Evaluate, nullptr);
}

DeviceEvents CalculatorDevice::Act(KeyAction Action, const KeyDef* Key)
{
	DeviceEvents Ev;
	Ev.ScreenChanged = true;
	if (Action != KeyAction::Complete && Action != KeyAction::Shift)
	{
		CompletionIndex = -1;
		Completions.clear();
	}
	const std::vector<std::string> Cats = Calc.Functions().Categories();
	if (LibCategory >= 0)
	{
		const std::vector<const mcalc::FunctionInfo*> Fns = Calc.Functions().InCategory(Cats[static_cast<size_t>(LibCategory)]);
		const int N = static_cast<int>(Fns.size());
		switch (Action)
		{
		case KeyAction::HistoryUp: LibSelection = std::max(0, LibSelection - 1); return Ev;
		case KeyAction::HistoryDown: LibSelection = std::min(N - 1, LibSelection + 1); return Ev;
		case KeyAction::PageUp: LibSelection = std::max(0, LibSelection - 3); return Ev;
		case KeyAction::PageDown: LibSelection = std::min(N - 1, LibSelection + 3); return Ev;
		case KeyAction::CursorLeft:
			LibCategory = (LibCategory + static_cast<int>(Cats.size()) - 1) % static_cast<int>(Cats.size());
			LibSelection = 0;
			return Ev;
		case KeyAction::CursorRight:
		case KeyAction::Menu:
			LibCategory = (LibCategory + 1) % static_cast<int>(Cats.size());
			LibSelection = 0;
			return Ev;
		case KeyAction::ClearLine: LibCategory = -1; return Ev;
		case KeyAction::Help:
			if (LibSelection < N)
			{
				const mcalc::Outcome O = Calc.Execute("help " + Fns[static_cast<size_t>(LibSelection)]->Name);
				for (const std::string& L : mcalc::WrapText(O.Text, 0))
				{
					AddTape(L, "  ");
				}
			}
			LibCategory = -1;
			return Ev;
		case KeyAction::Evaluate:
			if (LibSelection < N)
			{
				InsertText(Fns[static_cast<size_t>(LibSelection)]->Name + "(");
			}
			LibCategory = -1;
			return Ev;
		default:
			if (Action == KeyAction::Insert)
			{
				LibCategory = -1;   // typing closes the browser
			}
			break;
		}
	}
	switch (Action)
	{
	case KeyAction::Insert:
	{
		if (!Key)
		{
			break;
		}
		if (CtrlCount > 0 && (Key->Id == "L" || Key->Id == "C"))
		{
			if (Key->Id == "L")
			{
				TapeLines.clear();
			}
			else
			{
				Line.clear();
				CursorPos = 0;
			}
			break;
		}
		std::string Text = Key->Insert;
		const bool IsLetter = Text.size() == 1 && std::isalpha(static_cast<unsigned char>(Text[0])) != 0;
		const bool Upper = IsLetter ? (ShiftCount > 0) != Caps : ShiftCount > 0;
		if (Upper && !Key->ShiftInsert.empty())
		{
			Text = Key->ShiftInsert;
		}
		InsertText(Text);
		break;
	}
	case KeyAction::Evaluate: Evaluate(Ev); break;
	case KeyAction::Backspace: Backspace(CtrlCount > 0); break;
	case KeyAction::DeleteForward:
		if (CursorPos < static_cast<int>(Line.size()))
		{
			Line.erase(static_cast<size_t>(CursorPos), 1);
		}
		break;
	case KeyAction::ClearLine:
		if (Line.empty())
		{
			TapeLines.clear();
		}
		Line.clear();
		CursorPos = 0;
		HistoryIndex = -1;
		break;
	case KeyAction::Complete: Complete(); break;
	case KeyAction::CursorLeft:
		if (CtrlCount > 0)
		{
			while (CursorPos > 0 && !IsNameChar(Line[static_cast<size_t>(CursorPos - 1)])) --CursorPos;
			while (CursorPos > 0 && IsNameChar(Line[static_cast<size_t>(CursorPos - 1)])) --CursorPos;
		}
		else
		{
			CursorPos = std::max(0, CursorPos - 1);
		}
		break;
	case KeyAction::CursorRight:
		if (CtrlCount > 0)
		{
			const int N = static_cast<int>(Line.size());
			while (CursorPos < N && !IsNameChar(Line[static_cast<size_t>(CursorPos)])) ++CursorPos;
			while (CursorPos < N && IsNameChar(Line[static_cast<size_t>(CursorPos)])) ++CursorPos;
		}
		else
		{
			CursorPos = std::min(static_cast<int>(Line.size()), CursorPos + 1);
		}
		break;
	case KeyAction::Home: CursorPos = 0; break;
	case KeyAction::End: CursorPos = static_cast<int>(Line.size()); break;
	case KeyAction::HistoryUp:
		if (!History.empty())
		{
			if (HistoryIndex < 0)
			{
				Draft = Line;
				HistoryIndex = static_cast<int>(History.size()) - 1;
			}
			else if (HistoryIndex > 0)
			{
				--HistoryIndex;
			}
			Line = History[static_cast<size_t>(HistoryIndex)];
			CursorPos = static_cast<int>(Line.size());
		}
		break;
	case KeyAction::HistoryDown:
		if (HistoryIndex >= 0)
		{
			if (HistoryIndex + 1 < static_cast<int>(History.size()))
			{
				++HistoryIndex;
				Line = History[static_cast<size_t>(HistoryIndex)];
			}
			else
			{
				HistoryIndex = -1;
				Line = Draft;
			}
			CursorPos = static_cast<int>(Line.size());
		}
		break;
	case KeyAction::PageUp: TapeScroll = std::min(TapeScroll + 3, std::max(0, static_cast<int>(TapeLines.size()) - 1)); break;
	case KeyAction::PageDown: TapeScroll = std::max(0, TapeScroll - 3); break;
	case KeyAction::Shift: ++ShiftCount; Ev.ScreenChanged = true; break;
	case KeyAction::CapsLock: Caps = !Caps; break;
	case KeyAction::Control: ++CtrlCount; break;
	case KeyAction::Fn: ++FnCount; break;
	case KeyAction::Alt:
	case KeyAction::Super:
		Ev.ScreenChanged = false;
		break;
	case KeyAction::Menu:
		LibCategory = Cats.empty() ? -1 : 0;
		LibSelection = 0;
		break;
	case KeyAction::Help: Help(); break;
	case KeyAction::AngleMode:
	{
		const bool ToDeg = Calc.GetAngleMode() == mcalc::AngleMode::Radians;
		Calc.Execute(ToDeg ? "deg" : "rad");
		AddTape(ToDeg ? "angles in degrees" : "angles in radians", "  ");
		break;
	}
	case KeyAction::Digits:
	{
		int Next = kDigitSteps[0];
		for (int D : kDigitSteps)
		{
			if (D > Calc.GetDigits())
			{
				Next = D;
				break;
			}
		}
		Calc.SetDigits(Next);
		AddTape("digits " + std::to_string(Next), "  ");
		break;
	}
	case KeyAction::Explode:
		Explode = !Explode;
		Ev.ExplodeToggled = true;
		break;
	case KeyAction::Backlight:
		Backlight = (Backlight + 1) % kBacklightModes;
		Ev.BacklightCycled = true;
		AddTape(std::string("backlight ") + BacklightModeName(Backlight), "  ");
		break;
	case KeyAction::CopyResult:
	{
		mcalc::Value Ans;
		if (Calc.GetVariable("ans", &Ans))
		{
			InsertText("ans");
		}
		break;
	}
	case KeyAction::NumLock: Num = !Num; break;
	case KeyAction::None: Ev.ScreenChanged = false; break;
	}
	return Ev;
}

void CalculatorDevice::Release(const KeyDef& Key)
{
	switch (Key.Action)
	{
	case KeyAction::Shift: ShiftCount = std::max(0, ShiftCount - 1); break;
	case KeyAction::Control: CtrlCount = std::max(0, CtrlCount - 1); break;
	case KeyAction::Fn: FnCount = std::max(0, FnCount - 1); break;
	default: break;
	}
}

std::string CalculatorDevice::StatusLine(int Cols) const
{
	std::string Left = Calc.GetAngleMode() == mcalc::AngleMode::Degrees ? "DEG" : "RAD";
	Left += " " + std::to_string(Calc.GetDigits()) + "d";
	if (Caps) Left += " CAPS";
	if (ShiftCount > 0) Left += " SHIFT";
	if (FnCount > 0) Left += " FN";
	if (!Num) Left += " NAV";
	std::string Right;
	if (LibCategory >= 0)
	{
		Right = "LIBRARY";
	}
	else if (TapeScroll > 0)
	{
		Right = "^" + std::to_string(TapeScroll);
	}
	const int Gap = Cols - static_cast<int>(Left.size()) - static_cast<int>(Right.size());
	if (Gap < 1)
	{
		return Fit(Left, Cols);
	}
	return Left + std::string(static_cast<size_t>(Gap), ' ') + Right;
}

std::vector<std::string> CalculatorDevice::Screen(int Cols, int Rows, bool CursorOn) const
{
	std::vector<std::string> Out;
	if (Cols < 8 || Rows < 3)
	{
		return Out;
	}
	Out.push_back(StatusLine(Cols));
	const int Body = Rows - 2;
	if (LibCategory >= 0)
	{
		const std::vector<std::string> Cats = Calc.Functions().Categories();
		const std::string& Cat = Cats[static_cast<size_t>(LibCategory)];
		const std::vector<const mcalc::FunctionInfo*> Fns = Calc.Functions().InCategory(Cat);
		Out.push_back(Fit("< " + Cat + " >  " + std::to_string(LibSelection + 1) + "/" + std::to_string(Fns.size()), Cols));
		const int ListRows = Body - 1;
		int First = std::max(0, LibSelection - ListRows / 2);
		First = std::min(First, std::max(0, static_cast<int>(Fns.size()) - ListRows));
		for (int r = 0; r < ListRows; ++r)
		{
			const int i = First + r;
			if (i < static_cast<int>(Fns.size()))
			{
				Out.push_back(Fit((i == LibSelection ? "> " : "  ") + Fns[static_cast<size_t>(i)]->Signature(), Cols));
			}
			else
			{
				Out.push_back(Fit("", Cols));
			}
		}
	}
	else
	{
		std::vector<std::string> Wrapped;
		for (const std::string& L : TapeLines)
		{
			const std::vector<std::string> Parts = mcalc::WrapText(L, Cols);
			Wrapped.insert(Wrapped.end(), Parts.begin(), Parts.end());
		}
		const int End = std::max(0, static_cast<int>(Wrapped.size()) - TapeScroll);
		const int Begin = End - Body;
		for (int i = Begin; i < End; ++i)
		{
			Out.push_back(Fit(i >= 0 ? Wrapped[static_cast<size_t>(i)] : "", Cols));
		}
	}
	// Input line with a blinking block cursor that covers the character under it.
	const int Width = Cols - 2;
	int Start = 0;
	if (CursorPos >= Width)
	{
		Start = CursorPos - Width + 1;
	}
	std::string Visible = Line.size() > static_cast<size_t>(Start) ? Line.substr(static_cast<size_t>(Start), static_cast<size_t>(Width)) : std::string();
	Visible = Fit(Visible, Width);
	if (CursorOn)
	{
		Visible[static_cast<size_t>(CursorPos - Start)] = '_';
	}
	Out.push_back((Start > 0 ? "<" : ">") + std::string(" ") + Visible);
	return Out;
}
} // namespace mkb
