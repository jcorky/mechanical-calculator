// The keyboard's behaviour: typing and evaluating through the layout, editing keys, the screen,
// switch dynamics against the spec, and the synthesized sounds.
#include "Backlight.h"
#include "CalculatorDevice.h"
#include "KeySound.h"
#include "SwitchModel.h"
#include "TestFramework.h"

#include <algorithm>

namespace
{
const mkb::KeyDef& K(const char* Id)
{
	const mkb::KeyDef* D = mkb::FullSizeAnsi().Find(Id);
	if (!D)
	{
		std::printf("    no key %s\n", Id);
		std::abort();
	}
	return *D;
}

void Tap(mkb::CalculatorDevice& Dev, const char* Id)
{
	Dev.Press(K(Id));
	Dev.Release(K(Id));
}

// Types a string by pressing the physical keys, holding shift where the character needs it.
void TypeKeys(mkb::CalculatorDevice& Dev, const std::string& Text)
{
	const mkb::Layout& L = mkb::FullSizeAnsi();
	for (char Ch : Text)
	{
		const std::string S(1, Ch);
		const mkb::KeyDef* Plain = nullptr;
		const mkb::KeyDef* Shifted = nullptr;
		for (const mkb::KeyDef& Key : L.Keys)
		{
			if (Key.Role == mkb::KeyRole::Numpad || Key.Action != mkb::KeyAction::Insert)
			{
				continue;
			}
			if (Key.Insert == S && !Plain) Plain = &Key;
			if (Key.ShiftInsert == S && Key.Insert != S && !Shifted) Shifted = &Key;
		}
		if (Plain)
		{
			Tap(Dev, Plain->Id.c_str());
		}
		else if (Shifted)
		{
			Dev.Press(K("LShift"));
			Tap(Dev, Shifted->Id.c_str());
			Dev.Release(K("LShift"));
		}
		else
		{
			std::printf("    cannot type '%c'\n", Ch);
		}
	}
}
} // namespace

TEST(device_types_and_evaluates_through_keys)
{
	mkb::CalculatorDevice Dev;
	TypeKeys(Dev, "2*(3+4)^2");
	CHECK(Dev.Input() == "2*(3+4)^2");
	const mkb::DeviceEvents Ev = Dev.Press(K("Enter"));
	CHECK(Ev.Evaluated && !Ev.EvaluationFailed);
	CHECK(Dev.Tape().back() == "= 98");
	CHECK(Dev.Input().empty());
	// Numpad, then numpad Enter.
	for (const char* Id : {"Num1", "Num2", "NumMultiply", "Num3", "NumEnter"})
	{
		Tap(Dev, Id);
	}
	CHECK(Dev.Tape().back() == "= 36");
	// Shifted numpad types matrix syntax: [1,2;3,4]
	TypeKeys(Dev, "det([1");
	Dev.Press(K("LShift"));
	Tap(Dev, "NumAdd");        // ,
	Dev.Release(K("LShift"));
	Tap(Dev, "Num2");
	Dev.Press(K("LShift"));
	Tap(Dev, "NumDecimal");    // ;
	Dev.Release(K("LShift"));
	TypeKeys(Dev, "3,4])");
	CHECK(Dev.Input() == "det([1,2;3,4])");
	Dev.Enter();
	CHECK(Dev.Tape().back() == "= -2");
}

TEST(device_editing_keys)
{
	mkb::CalculatorDevice Dev;
	TypeKeys(Dev, "sqrt(16)");
	Tap(Dev, "Left");
	Tap(Dev, "Left");
	Tap(Dev, "Backspace");     // removes the 1
	CHECK(Dev.Input() == "sqrt(6)");
	Tap(Dev, "Home");
	Tap(Dev, "Delete");
	CHECK(Dev.Input() == "qrt(6)");
	Tap(Dev, "End");
	CHECK(Dev.Cursor() == 6);
	Dev.Press(K("LCtrl"));
	Tap(Dev, "Backspace");
	Dev.Release(K("LCtrl"));
	CHECK(Dev.Input() == "qrt(6");
	Tap(Dev, "Esc");
	CHECK(Dev.Input().empty());
	// Caps lock and shift.
	Tap(Dev, "Caps");
	TypeKeys(Dev, "ab");
	CHECK(Dev.Input() == "AB");
	Tap(Dev, "Caps");
	Dev.Press(K("RShift"));
	Tap(Dev, "C");
	Dev.Release(K("RShift"));
	CHECK(Dev.Input() == "ABC");
	// History.
	Tap(Dev, "Esc");
	TypeKeys(Dev, "1+1");
	Dev.Enter();
	TypeKeys(Dev, "2+2");
	Dev.Enter();
	TypeKeys(Dev, "draft");
	Tap(Dev, "Up");
	CHECK(Dev.Input() == "2+2");
	Tap(Dev, "Up");
	CHECK(Dev.Input() == "1+1");
	Tap(Dev, "Down");
	Tap(Dev, "Down");
	CHECK(Dev.Input() == "draft");
	// Num lock off: the numpad navigates.
	Tap(Dev, "NumLock");
	Tap(Dev, "Num4");
	CHECK(Dev.Cursor() == 4);
	Tap(Dev, "NumLock");
}

TEST(device_completion_help_and_function_keys)
{
	mkb::CalculatorDevice Dev;
	TypeKeys(Dev, "sqr");
	Tap(Dev, "Tab");
	CHECK(Dev.Input() == "sqrt(");
	TypeKeys(Dev, "2");
	Tap(Dev, "F1");            // help for the word at the cursor
	bool Found = false;
	for (const std::string& L : Dev.Tape())
	{
		Found = Found || L.find("sqrt(x)") != std::string::npos;
	}
	CHECK(Found);
	Tap(Dev, "Esc");
	Tap(Dev, "F2");            // degrees
	Tap(Dev, "F9");            // inserts sin(
	TypeKeys(Dev, "30)");
	Dev.Enter();
	CHECK(Dev.Tape().back() == "= 0.5");
	Tap(Dev, "F2");
	const mkb::DeviceEvents Ex = Dev.Press(K("Pause"));
	CHECK(Ex.ExplodeToggled && Dev.Exploded());
	const int Mode = Dev.BacklightMode();
	CHECK(Dev.Press(K("ScrLk")).BacklightCycled);
	CHECK(Dev.BacklightMode() == (Mode + 1) % mkb::CalculatorDevice::kBacklightModes);
	// Library browser: open, step, insert.
	Tap(Dev, "F4");
	CHECK(Dev.LibraryOpen());
	Tap(Dev, "Down");
	Tap(Dev, "Enter");
	CHECK(!Dev.LibraryOpen());
	CHECK(!Dev.Input().empty() && Dev.Input().back() == '(');
}

TEST(device_screen_layout)
{
	mkb::CalculatorDevice Dev;
	TypeKeys(Dev, "x = 7");
	Dev.Enter();
	TypeKeys(Dev, "x*6");
	const std::vector<std::string> S = Dev.Screen(32, 5, true);
	CHECK(S.size() == 5);
	for (const std::string& L : S)
	{
		CHECK(L.size() == 32 || &L == &S.back());
	}
	CHECK(S[0].compare(0, 3, "RAD") == 0);
	CHECK(S[3] == "= 7" + std::string(29, ' '));
	CHECK(S[4].compare(0, 6, "> x*6_") == 0);
	// A long line scrolls to keep the cursor visible.
	TypeKeys(Dev, "+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1+1");
	const std::vector<std::string> T = Dev.Screen(32, 5, true);
	CHECK(T[4][0] == '<' && T[4].back() == '_');
}

TEST(switch_dynamics_match_spec)
{
	for (mkb::SwitchFeel Feel : {mkb::SwitchFeel::Linear, mkb::SwitchFeel::Tactile, mkb::SwitchFeel::Clicky})
	{
		mkb::SwitchModel Sw(mkb::SwitchSpec::For(Feel));
		Sw.SetPressed(true);
		bool Actuated = false, Bottom = false, Click = false;
		double ActuatedAt = -1, BottomAt = -1;
		for (int i = 0; i < 400; ++i)   // 40 ms in 0.1 ms frames
		{
			const mkb::SwitchEvents E = Sw.Step(1e-4);
			if (E.Actuated && !Actuated)
			{
				Actuated = true;
				ActuatedAt = i * 0.1;
			}
			if (E.BottomedOut && !Bottom)
			{
				Bottom = true;
				BottomAt = i * 0.1;
			}
			Click = Click || E.Clicked;
		}
		CHECK(Actuated && Bottom);
		CHECK(ActuatedAt > 2 && ActuatedAt < BottomAt && BottomAt < 20);   // a real keystroke: a few ms
		CHECK(Click == (Feel == mkb::SwitchFeel::Clicky));
		CHECK_NEAR(Sw.PositionMm(), 4.0, 1e-9);
		Sw.SetPressed(false);
		bool Reset = false, Top = false;
		for (int i = 0; i < 1000 && !Sw.IsAtRest(); ++i)
		{
			const mkb::SwitchEvents E = Sw.Step(1e-4);
			Reset = Reset || E.Reset;
			Top = Top || E.ToppedOut;
		}
		CHECK(Reset && Top && Sw.IsAtRest() && !Sw.IsActuated());
	}
	// A big frame hitch steps safely and still lands at the bottom.
	mkb::SwitchModel Sw;
	Sw.SetPressed(true);
	Sw.Step(0.5);
	CHECK_NEAR(Sw.PositionMm(), 4.0, 1e-9);
	// The tactile bump: force rises then falls before actuation.
	const mkb::SwitchSpec Brown = mkb::SwitchSpec::For(mkb::SwitchFeel::Tactile);
	CHECK(Brown.SpringForce(1.4) > Brown.SpringForce(0.5) + 0.15);
	CHECK(Brown.SpringForce(2.0) < Brown.SpringForce(1.4));
}

TEST(key_sounds)
{
	mkb::KeySoundParams P;
	const std::vector<int16_t> A = mkb::SynthesizeKeySound(P);
	CHECK(A.size() > 3000 && A.size() < 8000);
	int Peak = 0;
	for (int16_t S : A)
	{
		Peak = std::max(Peak, std::abs(static_cast<int>(S)));
	}
	CHECK(Peak > 20000 && Peak <= 32767);
	CHECK(A.back() == 0);
	CHECK(mkb::SynthesizeKeySound(P) == A);   // deterministic
	P.Seed = 7;
	CHECK(mkb::SynthesizeKeySound(P) != A);
	P.Kind = mkb::KeySoundKind::Click;
	CHECK(mkb::SynthesizeKeySound(P).size() < A.size());
	P.Kind = mkb::KeySoundKind::BottomOut;
	P.KeyWidthU = 6.25;
	P.Stabilized = true;
	CHECK(mkb::SynthesizeKeySound(P).size() > A.size());
}

TEST(backlight_modes)
{
	mkb::Backlight B;
	const mkb::KeyDef& F = K("F");
	const mkb::KeyDef& Esc = K("Esc");
	B.SetMode(0);
	CHECK(B.ColorFor(F, 1.0).R == 0.0);
	B.SetMode(4);
	B.OnKeyPressed(F, 0.0);
	const mkb::Rgb Near = B.ColorFor(F, 0.0), FarAway = B.ColorFor(Esc, 0.0);
	CHECK(Near.B > FarAway.B + 0.5);
	B.Prune(5.0);
	CHECK(B.ColorFor(F, 5.0).B < 0.05);
	B.SetMode(2);
	const mkb::Rgb W = B.ColorFor(F, 0.3);
	CHECK(W.R >= 0 && W.R <= 1 && W.G >= 0 && W.G <= 1 && W.B >= 0 && W.B <= 1);
}
