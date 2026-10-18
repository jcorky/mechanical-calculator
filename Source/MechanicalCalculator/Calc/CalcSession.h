// The calculator the keyboard drives: a small matrix expression language over the function registry.
//
//   2 + 3*4            precedence, ^ right-associative, unary minus binds looser than ^ (-2^2 = -4)
//   A = [1, 2; 3, 4]   matrix literal: ',' separates columns, ';' rows
//   A' * A, A \ b      transpose, left division (solve), .* ./ .^ element-wise
//   h = hohmann(6778e3, 42164e3)    library call; records print their fields
//   h.dv1 + h.dv2, h(4)             field access, 1-based indexing
//   x = 2; y = x^2     ';' separates statements outside brackets
//
// Commands (whole line): help <name>, list [category], vars, clear, deg, rad, digits <n>.
#pragma once

#include "CalcLibrary.h"

#include <map>
#include <string>
#include <vector>

namespace mcalc
{
enum class AngleMode
{
	Radians,
	Degrees,
};

struct Outcome
{
	bool Ok = true;
	bool HasValue = false;     // false for assignments to nothing, commands, empty lines
	Value Val;
	std::string Text;          // what the display shows: a formatted value, help text or the error
	std::vector<std::string> Lines; // Text split for a narrow display
};

class Session
{
public:
	Session();
	explicit Session(const Registry& InRegistry);

	Outcome Execute(const std::string& Line);

	// Formats a value for a display `Width` characters wide (0 = unlimited) as one or more lines.
	std::vector<std::string> Format(const Value& V, int Width = 0) const;
	std::string FormatNumber(double X) const;

	void SetAngleMode(AngleMode Mode) { Angles = Mode; }
	AngleMode GetAngleMode() const { return Angles; }
	void SetDigits(int D) { Digits = D < 3 ? 3 : (D > 17 ? 17 : D); }
	int GetDigits() const { return Digits; }

	bool GetVariable(const std::string& Name, Value* Out) const;
	void SetVariable(const std::string& Name, const Value& V);
	const std::map<std::string, Value>& Variables() const { return Vars; }
	const Registry& Functions() const { return *Reg; }

	// Every name the input could be completing: functions, variables, constants and commands.
	std::vector<std::string> Complete(const std::string& Prefix) const;

	// Used by the evaluator.
	double ToRadians(double X) const;
	double FromRadians(double X) const;

private:
	const Registry* Reg;
	std::map<std::string, Value> Vars;
	AngleMode Angles = AngleMode::Radians;
	int Digits = 10;

	bool RunCommand(const std::string& Line, Outcome* Out);
};

// Wraps text to a display width, breaking at spaces where possible.
std::vector<std::string> WrapText(const std::string& Text, int Width);
} // namespace mcalc
