// The keyboard's calculator without Unreal: reads lines from stdin, prints what the OLED would.
//   ./mcalc_cli            interactive
//   echo "hohmann(6678e3, 42164e3)" | ./mcalc_cli
#include "CalcSession.h"

#include <cstdio>
#include <iostream>
#include <string>

int main()
{
	mcalc::Session S;
	std::string Line;
	while (std::printf("> "), std::fflush(stdout), std::getline(std::cin, Line))
	{
		const mcalc::Outcome O = S.Execute(Line);
		if (!O.Text.empty())
		{
			std::printf("%s%s\n", O.Ok ? (O.HasValue ? "= " : "") : "! ", O.Text.c_str());
		}
	}
	std::printf("\n");
	return 0;
}
