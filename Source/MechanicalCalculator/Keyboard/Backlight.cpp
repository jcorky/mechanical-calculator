#include "Backlight.h"

#include <algorithm>
#include <cmath>

namespace mkb
{
namespace backlight_detail
{
constexpr double kRippleSpeed = 9.0;   // key units per second
constexpr double kRippleWidth = 0.9;
constexpr double kRippleLife = 1.2;    // seconds
} // namespace backlight_detail

Rgb HsvToRgb(double H, double S, double V)
{
	H = H - std::floor(H);
	const double C = V * S;
	const double Hp = H * 6.0;
	const double X = C * (1.0 - std::fabs(std::fmod(Hp, 2.0) - 1.0));
	double R = 0, G = 0, B = 0;
	if (Hp < 1) { R = C; G = X; }
	else if (Hp < 2) { R = X; G = C; }
	else if (Hp < 3) { G = C; B = X; }
	else if (Hp < 4) { G = X; B = C; }
	else if (Hp < 5) { R = X; B = C; }
	else { R = C; B = X; }
	const double M = V - C;
	return Rgb{R + M, G + M, B + M};
}

void Backlight::OnKeyPressed(const KeyDef& Key, double Time)
{
	Ripples.push_back(Ripple{Key.CenterX(), Key.CenterY(), Time});
	if (Ripples.size() > 32)
	{
		Ripples.erase(Ripples.begin());
	}
}

void Backlight::Prune(double Time)
{
	Ripples.erase(std::remove_if(Ripples.begin(), Ripples.end(), [Time](const Ripple& R) { return Time - R.Start > backlight_detail::kRippleLife; }), Ripples.end());
}

Rgb Backlight::ColorFor(const KeyDef& Key, double Time) const
{
	using namespace backlight_detail;
	switch (Mode)
	{
	case 0:
		return Rgb{};
	case 1:
		return Rgb{1.0, 0.72, 0.42};   // warm white
	case 2:
		return HsvToRgb(Key.CenterX() / 22.5 - Time * 0.12, 0.85, 1.0);
	case 3:
	{
		const double Breath = 0.5 - 0.5 * std::cos(Time * 2.0 * 3.14159265358979 / 4.0);
		return HsvToRgb(0.52, 0.7, 0.08 + 0.92 * Breath * Breath);
	}
	case 4:
	{
		double I = 0.0;
		for (const Ripple& R : Ripples)
		{
			const double Age = Time - R.Start;
			if (Age < 0 || Age > kRippleLife)
			{
				continue;
			}
			const double D = std::hypot(Key.CenterX() - R.X, Key.CenterY() - R.Y);
			const double Ring = (D - kRippleSpeed * Age) / kRippleWidth;
			I += std::exp(-Ring * Ring) * (1.0 - Age / kRippleLife);
		}
		const Rgb C = HsvToRgb(0.55 + 0.1 * std::min(I, 1.0), 0.9, std::min(1.0, I));
		return Rgb{C.R + 0.02, C.G + 0.02, C.B + 0.03};
	}
	case 5:
		switch (Key.Role)
		{
		case KeyRole::Numpad: return Rgb{0.1, 0.8, 1.0};
		case KeyRole::Function: return Rgb{1.0, 0.45, 0.05};
		case KeyRole::Accent: return Rgb{1.0, 0.1, 0.15};
		case KeyRole::Navigation: return Rgb{0.6, 0.2, 1.0};
		case KeyRole::Modifier: return Rgb{0.25, 0.25, 0.3};
		default: return Rgb{0.9, 0.9, 0.85};
		}
	default:
		return Rgb{};
	}
}
} // namespace mkb
