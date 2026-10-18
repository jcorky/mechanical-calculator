#include "KeyLayout.h"

#include <cctype>

namespace mkb
{
namespace layout_detail
{
KeyDef Key(const char* Id, const char* Legend, double X, double Y, double W, int Row, KeyRole Role, KeyAction Action, const char* EngineKey)
{
	KeyDef K;
	K.Id = Id;
	K.Legend = Legend;
	K.X = X;
	K.Y = Y;
	K.W = W;
	K.ProfileRow = Row;
	K.Role = Role;
	K.Action = Action;
	K.EngineKey = EngineKey;
	return K;
}

KeyDef Text(const char* Id, const char* Legend, const char* Shifted, const char* Insert, const char* ShiftInsert, double X, double Y, int Row, const char* EngineKey)
{
	KeyDef K = Key(Id, Legend, X, Y, 1.0, Row, KeyRole::Alpha, KeyAction::Insert, EngineKey);
	K.ShiftLegend = Shifted;
	K.Insert = Insert;
	K.ShiftInsert = ShiftInsert;
	return K;
}

KeyDef Letter(char Upper, double X, double Y, int Row)
{
	const std::string Up(1, Upper);
	const std::string Low(1, static_cast<char>(std::tolower(static_cast<unsigned char>(Upper))));
	KeyDef K = Key(Up.c_str(), Up.c_str(), X, Y, 1.0, Row, KeyRole::Alpha, KeyAction::Insert, Up.c_str());
	K.Insert = Low;
	K.ShiftInsert = Up;
	return K;
}

Layout Build()
{
	Layout L;
	std::vector<KeyDef>& K = L.Keys;

	// Function row, y = 0.
	K.push_back(Key("Esc", "Esc", 0, 0, 1, 1, KeyRole::Accent, KeyAction::ClearLine, "Escape"));
	K.back().FrontLegend = "AC";
	struct FKey
	{
		const char* Id;
		const char* Front;
		KeyAction Action;
		const char* Insert;
	};
	const FKey FRow[12] = {
		{"F1", "HELP", KeyAction::Help, ""},
		{"F2", "DEG/RAD", KeyAction::AngleMode, ""},
		{"F3", "DIGITS", KeyAction::Digits, ""},
		{"F4", "LIB", KeyAction::Menu, ""},
		{"F5", "sqrt", KeyAction::Insert, "sqrt("},
		{"F6", "x^2", KeyAction::Insert, "^2"},
		{"F7", "ln", KeyAction::Insert, "ln("},
		{"F8", "exp", KeyAction::Insert, "exp("},
		{"F9", "sin", KeyAction::Insert, "sin("},
		{"F10", "cos", KeyAction::Insert, "cos("},
		{"F11", "tan", KeyAction::Insert, "tan("},
		{"F12", "pi", KeyAction::Insert, "pi"},
	};
	const double FX[12] = {2, 3, 4, 5, 6.5, 7.5, 8.5, 9.5, 11, 12, 13, 14};
	for (int i = 0; i < 12; ++i)
	{
		KeyDef F = Key(FRow[i].Id, FRow[i].Id, FX[i], 0, 1, 1, KeyRole::Function, FRow[i].Action, FRow[i].Id);
		F.FrontLegend = FRow[i].Front;
		F.Insert = FRow[i].Insert;
		K.push_back(F);
	}
	K.push_back(Key("PrtSc", "PrtSc", 15.25, 0, 1, 1, KeyRole::Function, KeyAction::CopyResult, ""));
	K.back().FrontLegend = "ANS";
	K.push_back(Key("ScrLk", "ScrLk", 16.25, 0, 1, 1, KeyRole::Function, KeyAction::Backlight, "ScrollLock"));
	K.back().FrontLegend = "LIGHT";
	K.push_back(Key("Pause", "Pause", 17.25, 0, 1, 1, KeyRole::Function, KeyAction::Explode, "Pause"));
	K.back().FrontLegend = "EXPLODE";

	// Number row, y = 1.5.
	const double R1 = 1.5;
	K.push_back(Text("Grave", "`", "~", "`", "~", 0, R1, 1, "Tilde"));
	const char* Digits = "1234567890";
	const char* Shifted[] = {"!", "@", "#", "$", "%", "^", "&", "*", "(", ")"};
	const char* DigitKeys[] = {"One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Zero"};
	for (int i = 0; i < 10; ++i)
	{
		const std::string D(1, Digits[i]);
		K.push_back(Text(("D" + D).c_str(), D.c_str(), Shifted[i], D.c_str(), Shifted[i], 1.0 + i, R1, 1, DigitKeys[i]));
	}
	K.push_back(Text("Minus", "-", "_", "-", "_", 11, R1, 1, "Hyphen"));
	K.push_back(Text("Equals", "=", "+", "=", "+", 12, R1, 1, "Equals"));
	K.push_back(Key("Backspace", "Backspace", 13, R1, 2, 1, KeyRole::Modifier, KeyAction::Backspace, "BackSpace"));

	// QWERTY row, y = 2.5.
	const double R2 = 2.5;
	K.push_back(Key("Tab", "Tab", 0, R2, 1.5, 2, KeyRole::Modifier, KeyAction::Complete, "Tab"));
	const char* Top = "QWERTYUIOP";
	for (int i = 0; i < 10; ++i)
	{
		K.push_back(Letter(Top[i], 1.5 + i, R2, 2));
	}
	K.push_back(Text("LBracket", "[", "{", "[", "{", 11.5, R2, 2, "LeftBracket"));
	K.push_back(Text("RBracket", "]", "}", "]", "}", 12.5, R2, 2, "RightBracket"));
	KeyDef Bs = Text("Backslash", "\\", "|", "\\", "|", 13.5, R2, 2, "Backslash");
	Bs.W = 1.5;
	K.push_back(Bs);

	// Home row, y = 3.5.
	const double R3 = 3.5;
	K.push_back(Key("Caps", "Caps Lock", 0, R3, 1.75, 3, KeyRole::Modifier, KeyAction::CapsLock, "CapsLock"));
	const char* Home = "ASDFGHJKL";
	for (int i = 0; i < 9; ++i)
	{
		K.push_back(Letter(Home[i], 1.75 + i, R3, 3));
	}
	K.push_back(Text("Semicolon", ";", ":", ";", ":", 10.75, R3, 3, "Semicolon"));
	K.push_back(Text("Quote", "'", "\"", "'", "\"", 11.75, R3, 3, "Apostrophe"));
	K.push_back(Key("Enter", "Enter", 12.75, R3, 2.25, 3, KeyRole::Accent, KeyAction::Evaluate, "Enter"));
	K.back().FrontLegend = "=";

	// Shift row, y = 4.5.
	const double R4 = 4.5;
	K.push_back(Key("LShift", "Shift", 0, R4, 2.25, 4, KeyRole::Modifier, KeyAction::Shift, "LeftShift"));
	const char* Bottom = "ZXCVBNM";
	for (int i = 0; i < 7; ++i)
	{
		K.push_back(Letter(Bottom[i], 2.25 + i, R4, 4));
	}
	K.push_back(Text("Comma", ",", "<", ",", "<", 9.25, R4, 4, "Comma"));
	K.push_back(Text("Period", ".", ">", ".", ">", 10.25, R4, 4, "Period"));
	K.push_back(Text("Slash", "/", "?", "/", "?", 11.25, R4, 4, "Slash"));
	K.push_back(Key("RShift", "Shift", 12.25, R4, 2.75, 4, KeyRole::Modifier, KeyAction::Shift, "RightShift"));

	// Bottom row, y = 5.5.
	const double R5 = 5.5;
	K.push_back(Key("LCtrl", "Ctrl", 0, R5, 1.25, 4, KeyRole::Modifier, KeyAction::Control, "LeftControl"));
	K.push_back(Key("LWin", "Win", 1.25, R5, 1.25, 4, KeyRole::Modifier, KeyAction::Super, "LeftCommand"));
	K.push_back(Key("LAlt", "Alt", 2.5, R5, 1.25, 4, KeyRole::Modifier, KeyAction::Alt, "LeftAlt"));
	KeyDef Space = Key("Space", "", 3.75, R5, 6.25, 4, KeyRole::Space, KeyAction::Insert, "SpaceBar");
	Space.Insert = " ";
	K.push_back(Space);
	K.push_back(Key("RAlt", "Alt", 10, R5, 1.25, 4, KeyRole::Modifier, KeyAction::Alt, "RightAlt"));
	K.push_back(Key("Fn", "Fn", 11.25, R5, 1.25, 4, KeyRole::Modifier, KeyAction::Fn, ""));
	K.push_back(Key("Menu", "Menu", 12.5, R5, 1.25, 4, KeyRole::Modifier, KeyAction::Menu, ""));
	K.back().FrontLegend = "LIB";
	K.push_back(Key("RCtrl", "Ctrl", 13.75, R5, 1.25, 4, KeyRole::Modifier, KeyAction::Control, "RightControl"));

	// Navigation cluster.
	K.push_back(Key("Insert", "Ins", 15.25, R1, 1, 1, KeyRole::Navigation, KeyAction::None, "Insert"));
	K.push_back(Key("Home", "Home", 16.25, R1, 1, 1, KeyRole::Navigation, KeyAction::Home, "Home"));
	K.push_back(Key("PgUp", "PgUp", 17.25, R1, 1, 1, KeyRole::Navigation, KeyAction::PageUp, "PageUp"));
	K.push_back(Key("Delete", "Del", 15.25, R2, 1, 2, KeyRole::Navigation, KeyAction::DeleteForward, "Delete"));
	K.push_back(Key("End", "End", 16.25, R2, 1, 2, KeyRole::Navigation, KeyAction::End, "End"));
	K.push_back(Key("PgDn", "PgDn", 17.25, R2, 1, 2, KeyRole::Navigation, KeyAction::PageDown, "PageDown"));
	K.push_back(Key("Up", "\xE2\x86\x91", 16.25, R4, 1, 4, KeyRole::Navigation, KeyAction::HistoryUp, "Up"));
	K.push_back(Key("Left", "\xE2\x86\x90", 15.25, R5, 1, 4, KeyRole::Navigation, KeyAction::CursorLeft, "Left"));
	K.push_back(Key("Down", "\xE2\x86\x93", 16.25, R5, 1, 4, KeyRole::Navigation, KeyAction::HistoryDown, "Down"));
	K.push_back(Key("Right", "\xE2\x86\x92", 17.25, R5, 1, 4, KeyRole::Navigation, KeyAction::CursorRight, "Right"));

	// Numpad. Shifted, the operator keys type the matrix syntax: \ .* ' , ;
	K.push_back(Key("NumLock", "Num", 18.5, R1, 1, 1, KeyRole::Numpad, KeyAction::NumLock, "NumLock"));
	auto NumText = [&K](const char* Id, const char* Legend, const char* Front, const char* Insert, const char* ShiftInsert, double X, double Y, int Row, const char* EngineKey) {
		KeyDef N = Key(Id, Legend, X, Y, 1, Row, KeyRole::Numpad, KeyAction::Insert, EngineKey);
		N.FrontLegend = Front;
		N.Insert = Insert;
		N.ShiftInsert = ShiftInsert;
		K.push_back(N);
		return &K.back();
	};
	NumText("NumDivide", "/", "\\", "/", "\\", 19.5, R1, 1, "Divide");
	NumText("NumMultiply", "*", ".*", "*", ".*", 20.5, R1, 1, "Multiply");
	NumText("NumSubtract", "-", "'", "-", "'", 21.5, R1, 1, "Subtract");
	NumText("Num7", "7", "", "7", "7", 18.5, R2, 2, "NumPadSeven");
	NumText("Num8", "8", "", "8", "8", 19.5, R2, 2, "NumPadEight");
	NumText("Num9", "9", "", "9", "9", 20.5, R2, 2, "NumPadNine");
	NumText("NumAdd", "+", ",", "+", ",", 21.5, R2, 2, "Add")->H = 2;
	NumText("Num4", "4", "", "4", "4", 18.5, R3, 3, "NumPadFour");
	NumText("Num5", "5", "", "5", "5", 19.5, R3, 3, "NumPadFive");
	NumText("Num6", "6", "", "6", "6", 20.5, R3, 3, "NumPadSix");
	NumText("Num1", "1", "", "1", "1", 18.5, R4, 4, "NumPadOne");
	NumText("Num2", "2", "", "2", "2", 19.5, R4, 4, "NumPadTwo");
	NumText("Num3", "3", "", "3", "3", 20.5, R4, 4, "NumPadThree");
	KeyDef NumEnter = Key("NumEnter", "Enter", 21.5, R4, 1, 4, KeyRole::Accent, KeyAction::Evaluate, "");
	NumEnter.H = 2;
	NumEnter.FrontLegend = "=";
	K.push_back(NumEnter);
	NumText("Num0", "0", "", "0", "0", 18.5, R5, 4, "NumPadZero")->W = 2;
	NumText("NumDecimal", ".", ";", ".", ";", 20.5, R5, 4, "Decimal");

	L.Width = 22.5;
	L.Height = 6.5;
	L.Display = Rect{18.5, 0.0, 4.0, 1.0};
	return L;
}
} // namespace layout_detail

const KeyDef* Layout::Find(const std::string& Id) const
{
	for (const KeyDef& K : Keys)
	{
		if (K.Id == Id)
		{
			return &K;
		}
	}
	return nullptr;
}

const KeyDef* Layout::FindByEngineKey(const std::string& EngineKey) const
{
	if (EngineKey.empty())
	{
		return nullptr;
	}
	for (const KeyDef& K : Keys)
	{
		if (K.EngineKey == EngineKey)
		{
			return &K;
		}
	}
	return nullptr;
}

const Layout& FullSizeAnsi()
{
	static const Layout Board = layout_detail::Build();
	return Board;
}
} // namespace mkb
