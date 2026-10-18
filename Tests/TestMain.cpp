// Runs every TEST, then every reference case in Tests/reference/*.tsv.
//
// A reference case is one line:  expression <TAB> expected <TAB> tolerance <TAB> provenance
// The expression is typed into a fresh calculator Session exactly as on the keyboard; its value
// (a scalar) must match `expected` within `tolerance`, written "rel:1e-9" or "abs:1e-6". An expected
// value of "error" means the line must be rejected. Lines starting with '#' are comments.
// MCALC_CASES=<file>.tsv runs one file. The .tsv files are generated from Where It Isn't's own JavaScript modules by Tests/reference/build.mjs.
#include "CalcSession.h"
#include "TestFramework.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace
{
std::vector<std::string> SplitTabs(const std::string& Line)
{
	std::vector<std::string> Out;
	std::string Cell;
	std::istringstream In(Line);
	while (std::getline(In, Cell, '\t'))
	{
		Out.push_back(Cell);
	}
	return Out;
}

int RunReferenceCases(int* Total)
{
	int Failed = 0;
	std::vector<std::filesystem::path> Files;
	const char* Only = std::getenv("MCALC_CASES");   // e.g. MCALC_CASES=orbits.tsv
	for (const auto& Entry : std::filesystem::directory_iterator(MCALC_CASES_DIR))
	{
		if (Entry.path().extension() == ".tsv" && (!Only || Entry.path().filename() == Only))
		{
			Files.push_back(Entry.path());
		}
	}
	std::sort(Files.begin(), Files.end());
	for (const auto& Path : Files)
	{
		std::ifstream In(Path);
		std::string Line;
		int LineNo = 0, Here = 0, HereFailed = 0;
		while (std::getline(In, Line))
		{
			++LineNo;
			if (Line.empty() || Line[0] == '#')
			{
				continue;
			}
			const std::vector<std::string> Cells = SplitTabs(Line);
			if (Cells.size() < 3)
			{
				std::printf("  BAD %s:%d: needs expression, expected, tolerance\n", Path.filename().c_str(), LineNo);
				++Failed;
				++HereFailed;
				continue;
			}
			++Here;
			mcalc::Session S;
			const mcalc::Outcome O = S.Execute(Cells[0]);
			bool Pass = false;
			std::string Got;
			if (Cells[1] == "error")
			{
				Pass = !O.Ok;
				Got = O.Ok ? O.Text : "error";
			}
			else if (!O.Ok || !O.HasValue || !O.Val.IsScalar())
			{
				Got = O.Ok ? "non-scalar " + O.Text : "error: " + O.Text;
			}
			else
			{
				// Compared on the machine: the expected value is entered exactly as printed.
				size_t Used = 0;
				const mcalc::Real Expected = mcalc::Real::Parse(Cells[1].c_str(), &Used);
				const mcalc::Real Value = O.Val.Scalar();
				const std::string& T = Cells[2];
				const mcalc::Real Tol = mcalc::Real::Parse(T.c_str() + 4, &Used);
				const mcalc::Real Allowed = T.compare(0, 4, "rel:") == 0 ? Tol * mcalc::Real::Abs(Expected) : Tol;
				if (Value.IsNaN() || Expected.IsNaN())
				{
					Pass = Value.IsNaN() && Expected.IsNaN();
				}
				else if (Value.IsInf() || Expected.IsInf())
				{
					Pass = Value == Expected;
				}
				else
				{
					Pass = mcalc::Real::Abs(Value - Expected) <= Allowed;
				}
				Got = Value.ToString(18);
			}
			if (!Pass)
			{
				++Failed;
				++HereFailed;
				std::printf("  FAIL %s:%d: %s -> %s, expected %s (%s)%s%s\n", Path.filename().c_str(), LineNo, Cells[0].c_str(), Got.c_str(), Cells[1].c_str(),
					Cells[2].c_str(), Cells.size() > 3 ? " " : "", Cells.size() > 3 ? Cells[3].c_str() : "");
			}
		}
		std::printf("reference %-28s %4d cases, %d failed\n", Path.filename().c_str(), Here, HereFailed);
		*Total += Here;
	}
	return Failed;
}
} // namespace

int main(int Argc, char** Argv)
{
	const char* Filter = Argc > 1 ? Argv[1] : nullptr;
	for (const mtest::Case& C : mtest::Cases())
	{
		if (Filter && std::string(C.Name).find(Filter) == std::string::npos)
		{
			continue;
		}
		const int Before = mtest::Failures();
		C.Fn();
		std::printf("%s %s\n", mtest::Failures() == Before ? "ok  " : "FAIL", C.Name);
	}
	int RefTotal = 0;
	const int RefFailed = Filter ? 0 : RunReferenceCases(&RefTotal);
	std::printf("\n%d checks, %d failed; %d reference cases, %d failed\n", mtest::Checks(), mtest::Failures(), RefTotal, RefFailed);
	return (mtest::Failures() == 0 && RefFailed == 0) ? 0 : 1;
}
