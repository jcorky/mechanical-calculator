// Properties of the radar, signal and lens libraries that are not single reference numbers.
#include "CalcSession.h"
#include "TestFramework.h"

#include <cmath>
#include <string>

namespace
{
double RadarEval(mcalc::Session& S, const std::string& Line)
{
	const mcalc::Outcome O = S.Execute(Line);
	if (!O.Ok || !O.HasValue || !O.Val.IsScalar())
	{
		std::printf("    '%s' -> %s\n", Line.c_str(), O.Text.c_str());
		return std::nan("");
	}
	return O.Val.Scalar().ToDouble();
}
} // namespace

TEST(radar_library_is_documented)
{
	const mcalc::Registry& R = mcalc::Registry::Default();
	int Count = 0;
	for (const char* Cat : {"radar", "signal", "lens"})
	{
		for (const mcalc::FunctionInfo* F : R.InCategory(Cat))
		{
			++Count;
			CHECK(!F->Source.empty() && !F->Summary.empty() && !F->Returns.empty());
			for (const mcalc::Param& P : F->Params)
			{
				CHECK(!P.Unit.empty());
			}
		}
	}
	CHECK(Count >= 50);
}

TEST(radar_equation_forms_agree)
{
	mcalc::Session S;
	const std::string Args = "1.4e6, undb(33), c0/2.8e9, 1, 60*nmi, 950, 1.67e6, undb(8)";
	const double Snr = RadarEval(S, "snr_track(" + Args + ")");
	CHECK_NEAR(RadarEval(S, "snr_chain(" + Args + ").snr") / Snr, 1.0, 1e-12);
	CHECK_NEAR(RadarEval(S, "snr_budget(" + Args + ").total"), 10 * std::log10(Snr), 1e-9);
	// The signed terms of the budget sum to its total.
	S.Execute("b = snr_budget(" + Args + ")");
	CHECK_NEAR(RadarEval(S, "b.pt + b.g2 + b.lam2 + b.sigma + b.fourpi3 + b.r4 + b.inv_k + b.ts + b.loss + b.bn"), RadarEval(S, "b.total"), 1e-9);
	// range_for inverts snr_track.
	CHECK_NEAR(RadarEval(S, "snr_track(1.4e6, undb(33), 0.1, 1, range_for(20, 1.4e6, undb(33), 0.1, 1), 950, 1.67e6, undb(8))"), 20.0, 1e-9);
	// The search form, its inverses and the track form summed over beam positions.
	CHECK_NEAR(RadarEval(S, "snr_search(1000, 3, 10, 1, 1, range_search(1000, 3, 10, 1, 1, 500, 10, 30), 500, 10)"), 30.0, 1e-9);
	CHECK_NEAR(RadarEval(S, "pav_search(80e3, 3, 10, 1, 1, 500, 10, snr_search(1234, 3, 10, 1, 1, 80e3, 500, 10))"), 1234.0, 1e-8);
	CHECK_NEAR(RadarEval(S, "search_track(1000, 3, 10, 1, 1, 80e3, 500, 10).snr / snr_search(1000, 3, 10, 1, 1, 80e3, 500, 10)"), 1.0, 1e-12);
	CHECK_NEAR(RadarEval(S, "tsys(200, 1.5, 300, 2).ts"), RadarEval(S, "tsys_trace(200, 1.5, 300, 2)"), 1e-9);
}

TEST(radar_detection_properties)
{
	mcalc::Session S;
	// Marcum Q1 is the single-pulse steady-target Pd.
	CHECK_NEAR(RadarEval(S, "marcumq(sqrt(2*7), sqrt(-2*ln(1e-5)))"), RadarEval(S, "pd_steady(7, 1, 1e-5)"), 1e-9);
	// Pd rises with S/N (to the series' truncation, about 1e-12 near Pd = 1, as in the module) and with
	// the number of pulses; snr_req_db inverts pd_steady.
	double Last = 0.0;
	for (int Db = -5; Db <= 20; Db += 5)
	{
		const double Pd = RadarEval(S, "pd_steady(undb(" + std::to_string(Db) + "), 4)");
		CHECK(Pd >= Last - 1e-11);
		Last = Pd;
	}
	CHECK(RadarEval(S, "pd_steady(2, 8)") > RadarEval(S, "pd_steady(2, 4)"));
	CHECK_NEAR(RadarEval(S, "pd_steady(undb(snr_req_db(0.8, 6, 1e-5)), 6, 1e-5)"), 0.8, 1e-9);
	// Albersheim is within a fraction of a dB of the exact series.
	CHECK(std::fabs(RadarEval(S, "albersheim(0.9, 1e-6, 21) - snr_req_db(0.9, 21)")) < 0.3);
	// Noncoherent integration gains less than coherent.
	CHECK(RadarEval(S, "intgain(21).loss") > 0.0);
}

TEST(radar_gert_matches_closed_forms)
{
	mcalc::Session S;
	// k detections in a row: mean T (1 - p^k) / ((1 - p) p^k).
	for (double P : {0.3, 0.55, 0.9})
	{
		const std::string Ps = std::to_string(P);
		CHECK_NEAR(RadarEval(S, "ttrack(" + Ps + ", 2, 1).mean"), 2 * (1 - P * P) / ((1 - P) * P * P), 1e-9);
		CHECK_NEAR(RadarEval(S, "ttrack(" + Ps + ", 2, 3).mean"), 2 * (1 - P * P * P) / ((1 - P) * P * P * P), 1e-8);
		CHECK_NEAR(RadarEval(S, "gert(trk_net(" + Ps + ", 2, 2)).var"), RadarEval(S, "ttrack(" + Ps + ", 2, 2).var"), 1e-9);
	}
	// A geometric wait: one branch to the end with probability p per scan of length T.
	CHECK_NEAR(RadarEval(S, "gert([1, 2, 0.25, 3; 1, 1, 0.75, 3]).mean"), 12.0, 1e-9);
	CHECK_NEAR(RadarEval(S, "gert([1, 2, 0.25, 3; 1, 1, 0.75, 3]).var"), 9.0 * 0.75 / (0.25 * 0.25), 1e-8);
	// Start and finish can be named explicitly.
	CHECK_NEAR(RadarEval(S, "gert([1, 2, 1, 1; 2, 3, 1, 2], 2, 3).mean"), 2.0, 1e-12);
}

TEST(radar_basics_and_lens_properties)
{
	mcalc::Session S;
	// Folding and aliasing land inside the unambiguous intervals.
	for (double V : {-700.0, -75.0, 0.0, 99.0, 1234.0})
	{
		const double A = RadarEval(S, "valias(" + std::to_string(V) + ", 10e3, 0.03)");
		CHECK(A >= -75.0 && A < 75.0);
		CHECK_NEAR(std::remainder(A - V, 150.0), 0.0, 1e-9);
	}
	CHECK(RadarEval(S, "rfold(1e6, 1e3)") < RadarEval(S, "r_unamb(1e3)"));
	// The 4.12 rule tracks the exact horizon; frame_time forgives degrees typed in radians.
	CHECK(std::fabs(RadarEval(S, "r_horizon(30, 5000) / r_horizon_rule(30, 5000)") - 1) < 2e-3);
	CHECK_NEAR(RadarEval(S, "frame_time(rad(90), rad(1.5), 2, 0.01)"), 60 * 2 * 0.01, 1e-12);
	// Mixing rules stay inside the Wiener bounds; fill_for inverts mix_eps.
	for (double F = 0.0; F <= 1.0; F += 0.125)
	{
		const std::string Fs = std::to_string(F);
		const double Lo = RadarEval(S, "mix_eps(" + Fs + ", 3, 4)"), Hi = RadarEval(S, "mix_eps(" + Fs + ", 3, 1)");
		for (int Rule : {2, 3})
		{
			const double E = RadarEval(S, "mix_eps(" + Fs + ", 3, " + std::to_string(Rule) + ")");
			CHECK(E >= Lo - 1e-12 && E <= Hi + 1e-12);
		}
	}
	CHECK_NEAR(RadarEval(S, "mix_eps(fill_for(1.6, 2.8, 3), 2.8, 3)"), 1.6, 1e-12);
	// The ray tracer's axial loss agrees with quadrature; the ideal lens focuses at (1, 0).
	CHECK(std::fabs(RadarEval(S, "lune_trace(1e-6, 1, 2, 2e-4, 0.0043, 2.8).loss / ax_loss(0.0043, 2.8)") - 1) < 3e-3);
	CHECK(RadarEval(S, "lune_spot(1, 2, 8, 2e-4).worst") < 1e-3);
	CHECK_NEAR(RadarEval(S, "lune_trace(0.6, 1, 2, 2e-4).opl"), 3.14159265358979323846 / 2 + 1, 1e-3);
	// A reflector 16x larger in cross section doubles nothing but the fourth root.
	CHECK_NEAR(RadarEval(S, "rcs_lune(0.2, 0.03) / rcs_lune(0.1, 0.03)"), 16.0, 1e-12);
}
