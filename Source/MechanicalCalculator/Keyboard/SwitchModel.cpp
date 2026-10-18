#include "SwitchModel.h"

#include <algorithm>
#include <cmath>

namespace mkb
{
SwitchSpec SwitchSpec::For(SwitchFeel Feel)
{
	SwitchSpec S;
	S.Feel = Feel;
	switch (Feel)
	{
	case SwitchFeel::Linear:
		S.PreloadN = 0.30;
		S.BottomForceN = 0.60;
		S.BumpN = 0.0;
		break;
	case SwitchFeel::Tactile:
		S.BumpN = 0.20;
		S.BumpAtMm = 1.4;
		break;
	case SwitchFeel::Clicky:
		S.PreloadN = 0.35;
		S.BottomForceN = 0.65;
		S.BumpN = 0.30;
		S.BumpAtMm = 1.7;
		S.BumpWidthMm = 0.3;
		S.ActuationMm = 2.2;
		break;
	}
	return S;
}

double SwitchSpec::SpringForce(double XMm) const
{
	const double Linear = PreloadN + (BottomForceN - PreloadN) * XMm / TravelMm;
	if (BumpN <= 0.0)
	{
		return Linear;
	}
	// The bump rises to its peak and falls off sharply past it, like the leaf riding over the stem's legs.
	const double D = (XMm - BumpAtMm) / BumpWidthMm;
	const double Bump = D < 0.0 ? std::exp(-D * D) : std::exp(-4.0 * D * D);
	return Linear + BumpN * Bump;
}

SwitchEvents SwitchModel::Step(double Dt)
{
	SwitchEvents Ev;
	const double Travel = Spec.TravelMm / 1000.0;
	const double ActuateAt = Spec.ActuationMm / 1000.0;
	const double ResetAt = (Spec.ActuationMm - Spec.HysteresisMm) / 1000.0;
	const double ClickAt = (Spec.BumpAtMm + 0.15) / 1000.0;
	double Remaining = std::min(Dt, 0.25);   // a long hitch must not tunnel through anything
	const double H = 1.0 / 20000.0;
	while (Remaining > 0.0)
	{
		const double Sub = std::min(H, Remaining);
		Remaining -= Sub;
		if (!Pressed && X <= 0.0 && V <= 0.0)
		{
			X = 0.0;
			V = 0.0;
			continue;
		}
		const double Spring = Spec.SpringForce(X * 1000.0);
		double Force = (Pressed ? Spec.FingerN : 0.0) - Spring;
		if (V != 0.0)
		{
			Force -= V > 0.0 ? Spec.FrictionN : -Spec.FrictionN;
		}
		else if (std::fabs(Force) <= Spec.FrictionN)
		{
			Force = 0.0;
		}
		// Semi-implicit Euler: stable for this stiff-ish spring at 50 us steps.
		V += Force / Spec.MassKg * Sub;
		const double Before = X;
		X += V * Sub;
		if (X >= Travel)
		{
			Ev.BottomedOut = Ev.BottomedOut || V > 0.02;
			Ev.BottomSpeed = std::max(Ev.BottomSpeed, V);
			X = Travel;
			V = V > 0.02 ? -V * Spec.BottomRestitution : 0.0;
		}
		if (X <= 0.0)
		{
			if (V < -0.02)
			{
				Ev.ToppedOut = true;
				Ev.TopSpeed = std::max(Ev.TopSpeed, -V);
			}
			X = 0.0;
			V = V < -0.02 ? -V * Spec.TopRestitution : 0.0;
		}
		if (!Actuated && X >= ActuateAt)
		{
			Actuated = true;
			Ev.Actuated = true;
		}
		else if (Actuated && X < ResetAt)
		{
			Actuated = false;
			Ev.Reset = true;
		}
		if (Spec.Feel == SwitchFeel::Clicky)
		{
			if (ClickArmed && Before < ClickAt && X >= ClickAt)
			{
				Ev.Clicked = true;
				ClickArmed = false;
			}
			else if (!ClickArmed && X < ResetAt)
			{
				ClickArmed = true;
			}
		}
	}
	return Ev;
}
} // namespace mkb
