#include "KeySound.h"

#include <algorithm>
#include <cmath>

namespace mkb
{
namespace sound_detail
{
constexpr double kTwoPi = 6.28318530717958647692;

// xorshift32: deterministic, so the same key always sounds the same for the same seed.
struct Noise
{
	uint32_t State;
	explicit Noise(uint32_t Seed) : State(Seed ? Seed : 0x9E3779B9u) {}
	double Next()
	{
		State ^= State << 13;
		State ^= State >> 17;
		State ^= State << 5;
		return static_cast<double>(State) / 4294967296.0 * 2.0 - 1.0;
	}
};

// Two-pole resonator (band-pass) on a noise source.
struct Resonator
{
	double A1 = 0, A2 = 0, Gain = 0, Y1 = 0, Y2 = 0;
	Resonator(double Freq, double Bandwidth, double Rate)
	{
		const double R = std::exp(-kTwoPi * Bandwidth / Rate);
		A1 = 2.0 * R * std::cos(kTwoPi * Freq / Rate);
		A2 = -R * R;
		Gain = 1.0 - R;
	}
	double Process(double X)
	{
		const double Y = Gain * X + A1 * Y1 + A2 * Y2;
		Y2 = Y1;
		Y1 = Y;
		return Y;
	}
};
} // namespace sound_detail

std::vector<int16_t> SynthesizeKeySound(const KeySoundParams& P)
{
	using namespace sound_detail;
	const double Rate = static_cast<double>(P.SampleRate);
	const double Width = std::max(1.0, P.KeyWidthU);
	const bool Down = P.Kind == KeySoundKind::BottomOut;
	const bool Click = P.Kind == KeySoundKind::Click;
	const double Duration = Click ? 0.025 : (Down ? 0.09 + 0.02 * std::log(Width) : 0.06);
	const size_t N = static_cast<size_t>(Duration * Rate);
	std::vector<double> Out(N, 0.0);
	Noise Rng(P.Seed);

	// Slight per-seed detune so a run of the same key does not sound machine-gunned.
	const double Detune = 1.0 + 0.04 * Rng.Next();
	if (Click)
	{
		Resonator Hi(4200.0 * Detune, 900.0, Rate), Mid(2600.0 * Detune, 700.0, Rate);
		for (size_t i = 0; i < N; ++i)
		{
			const double T = static_cast<double>(i) / Rate;
			const double Env = std::exp(-T / 0.0018);
			const double X = Rng.Next() * Env;
			Out[i] = 3.2 * Hi.Process(X) + 1.6 * Mid.Process(X);
		}
	}
	else
	{
		// Impact: the stem on the housing (down) or the top stop (up, higher and lighter).
		const double ImpactFreq = (Down ? 2900.0 : 3600.0) * Detune;
		Resonator Impact(ImpactFreq, 1800.0, Rate);
		// Case and plate modes: wider keys excite lower, longer modes.
		const double BodyFreq = (Down ? 520.0 : 760.0) / std::pow(Width, 0.35) * Detune;
		const double Modes[3][3] = {{BodyFreq, 0.020, 1.0}, {BodyFreq * 2.37, 0.010, 0.45}, {BodyFreq * 4.1, 0.006, 0.25}};
		const double Phase = Rng.Next() * kTwoPi;
		for (size_t i = 0; i < N; ++i)
		{
			const double T = static_cast<double>(i) / Rate;
			const double Hit = Rng.Next() * std::exp(-T / 0.0028);
			double S = 2.4 * Impact.Process(Hit);
			for (const auto& M : Modes)
			{
				S += M[2] * 0.35 * std::exp(-T / M[1]) * std::sin(kTwoPi * M[0] * T + Phase) * std::min(1.0, T * Rate / 24.0);
			}
			Out[i] = S * (Down ? 1.0 : 0.55);
		}
		if (P.Stabilized)
		{
			// Wire rattle: a few quiet ticks 8-30 ms after the impact.
			for (int k = 0; k < 3; ++k)
			{
				const size_t At = static_cast<size_t>((0.008 + 0.009 * k + 0.003 * Rng.Next()) * Rate);
				Resonator Tick(1700.0 + 300.0 * k, 800.0, Rate);
				for (size_t i = At; i < N && i < At + static_cast<size_t>(0.004 * Rate); ++i)
				{
					const double T = static_cast<double>(i - At) / Rate;
					Out[i] += 0.5 * Tick.Process(Rng.Next() * std::exp(-T / 0.0008));
				}
			}
		}
	}
	double Peak = 0.0;
	for (double X : Out)
	{
		Peak = std::max(Peak, std::fabs(X));
	}
	const double Level = std::clamp(0.35 + 0.45 * P.Strength, 0.1, 0.9) * (Click ? 0.9 : (Down ? 1.0 : 0.6));
	std::vector<int16_t> Pcm(N);
	for (size_t i = 0; i < N; ++i)
	{
		// Short fade-out so the buffer ends at silence.
		const double Fade = std::min(1.0, static_cast<double>(N - i) / (0.002 * Rate));
		const double V = Peak > 0.0 ? Out[i] / Peak * Level * Fade : 0.0;
		Pcm[i] = static_cast<int16_t>(std::lround(std::clamp(V, -1.0, 1.0) * 32767.0));
	}
	return Pcm;
}
} // namespace mkb
