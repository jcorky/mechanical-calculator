// Reference cases for LibRadar.cpp: radar equation and detection (radar/radar.js), radar basics
// (radar/basics.js) and lens antennas (lens/lens.js), computed by Where It Isn't's own modules.
// Marcum Q, which the module only uses inside its Pd series, is checked by direct numerical
// integration of its defining integral; a few cases use closed forms the benches quote.
// Run:  node Tests/reference/gen_radar.mjs /home/user/where-it-isnt > Tests/reference/radar.tsv
import path from "node:path";
import { pathToFileURL } from "node:url";

const wii = process.argv[2];
if (!wii) { console.error("usage: node gen_radar.mjs <where-it-isnt checkout>"); process.exit(2); }
await import(pathToFileURL(path.join(wii, "radar/radar.js")));
await import(pathToFileURL(path.join(wii, "radar/basics.js")));
await import(pathToFileURL(path.join(wii, "lens/lens.js")));
const R = globalThis.Radar, B = globalThis.RadarBasics, L = globalThis.Lens;
const C = R.C, rad = d => d * Math.PI / 180, deg = r => r * 180 / Math.PI;

const out = [];
const add = (expr, val, tol, prov) => {
  if (val !== "error") {
    if (typeof val === "boolean") val = val ? 1 : 0;
    if (!Number.isFinite(val)) throw new Error(`non-finite expected value for ${expr}: ${val}`);
  }
  out.push([expr, String(val), tol, prov].join("\t"));
};
const err = (expr, prov) => add(expr, "error", "abs:0", prov);
const head = t => out.push(`# ${t}`);

/* ---------- range equation ---------- */
head("radar: range equation (radar/radar.js), lecture example of O'Donnell 2002");
const printed = [1.4e6, 2000, 0.1, 1, 111e3, 950, 1.67e6, 6.3];
const exact = [1.4e6, R.undb(33), C / 2.8e9, 1, 60 * 1852, 950, 1.67e6, R.undb(8)];
const benchDef = [1.4e6, R.undb(33), C / 2.998e9, 1, 111e3, 950, 1.67e6, R.undb(8)];
const exactExpr = "1.4e6, undb(33), c0/2.8e9, 1, 60*nmi, 950, 1.67e6, undb(8)";
add("db(snr_track(1.4e6, 2000, 0.1, 1, 111e3, 950, 1.67e6, 6.3))", R.db(R.snrTrack(...printed)), "rel:1e-9", "Radar.snrTrack; lecture as printed, 1.29 dB");
add(`db(snr_track(${exactExpr}))`, R.db(R.snrTrack(...exact)), "rel:1e-9", "Radar.snrTrack; lecture inputs unrounded, 1.84 dB (README correction)");
add(`db(21*snr_track(${exactExpr}))`, R.db(21 * R.snrTrack(...exact)), "rel:1e-9", "Radar.snrTrack x 21 pulses coherent, 15.1 dB per dwell");
add("snr_track()", R.snrTrack(...benchDef), "rel:1e-9", "Radar.snrTrack; bench defaults (2.998 GHz, 111 km)");
add("snr_track(1e6, 1000, 0.03, 10, 50e3, 500, 1e6, 2)", R.snrTrack(1e6, 1000, 0.03, 10, 50e3, 500, 1e6, 2), "rel:1e-9", "Radar.snrTrack; X-band point");
err("snr_track(1.4e6, 2000, 0.1, 1, 0, 950, 1.67e6, 6.3)", "range must be > 0");
const ch = R.snrChained(...exact);
add(`snr_chain(${exactExpr}).snr`, ch.snr, "rel:1e-9", "Radar.snrChained; equals the formula");
add(`snr_chain(${exactExpr}).ae`, ch.Ae, "rel:1e-9", "Radar.snrChained effective area");
add(`snr_chain(${exactExpr}).s_tgt`, ch.atTarget, "rel:1e-9", "Radar.snrChained density at target");
add(`snr_chain(${exactExpr}).pr`, ch.Pr, "rel:1e-9", "Radar.snrChained received power");
add(`snr_chain(${exactExpr}).noise`, ch.N, "rel:1e-9", "Radar.snrChained kTsBn");
const bu = R.budget(...exact);
add(`snr_budget(${exactExpr}).plus`, bu.plus, "rel:1e-9", "Radar.budget; slide 356.1 dB");
add(`snr_budget(${exactExpr}).minus`, bu.minus, "rel:1e-9", "Radar.budget; slide 354.2 dB");
add(`snr_budget(${exactExpr}).total`, bu.total, "rel:1e-9", "Radar.budget total = S/N dB");
add(`snr_budget(${exactExpr}).r4`, -bu.rows[5][2], "rel:1e-9", "Radar.budget range^4 row (minus column, so negative here)");
const bu2 = R.budget(1e6, 1000, 0.03, 0.1, 50e3, 500, 1e6, 2);
add("snr_budget(1e6, 1000, 0.03, 0.1, 50e3, 500, 1e6, 2).sigma", -bu2.rows[3][2], "rel:1e-9", "Radar.budget; a -10 dBsm target goes to the minus column");
add("snr_budget(1e6, 1000, 0.03, 0.1, 50e3, 500, 1e6, 2).total", bu2.total, "rel:1e-9", "Radar.budget total");
add("range_for()", R.rangeFor(R.undb(13.18), 1.4e6, R.undb(33), C / 2.998e9, 1, 950, 1.67e6, R.undb(8)), "rel:1e-9", "Radar.rangeFor; range for 13.18 dB at bench defaults");
add(`range_for(undb(13.18), ${exactExpr.replace(", 60*nmi", "")})`, R.rangeFor(R.undb(13.18), ...exact.slice(0, 4), 950, 1.67e6, R.undb(8)), "rel:1e-9", "Radar.rangeFor; unrounded lecture radar");
err("range_for(0)", "required S/N must be > 0");

/* ---------- search form ---------- */
head("radar: search form (radar/radar.js section 2)");
const A1 = Math.PI / 4, A0 = Math.PI * Math.pow(10, 0.18) * Math.pow(10, 0.18) / 4;
add("pav_search(100e3, pi/4, 10, 1, 1, 500, undb(10), undb(15))", R.pavNeeded(100e3, A1, 10, 1, 1, 500, R.undb(10), R.undb(15)), "rel:1e-9", "Radar.pavNeeded; chart 349.3 W for 100 km, 1 m");
add("pav_search(200e3, pi/4)/pav_search(100e3, pi/4)", 16, "rel:1e-9", "R^4 law: twice the range needs 16x power");
add("pav_search(30e3, pi*0.1^2/4)", R.pavNeeded(30e3, Math.PI * 0.01 / 4, 10, 1, 1, 500, 10, R.undb(15)), "rel:1e-9", "Radar.pavNeeded; chart line end 283 W");
add("snr_search(1000, 10, 10, 1, 1, 100e3, 500, 10)", R.snrSearch(1000, 10, 10, 1, 1, 100e3, 500, 10), "rel:1e-9", "Radar.snrSearch");
add("snr_search()", R.snrSearch(1000, A0, 10, 1, 1, 100e3, 500, 10), "rel:1e-9", "Radar.snrSearch; bench defaults at 100 km");
err("snr_search(1000, 10, 10, 1, 0, 100e3, 500, 10)", "solid angle must be > 0");
add("range_search()", R.searchRange(1000, A0, 10, 1, 1, 500, 10, R.undb(15)), "rel:1e-9", "Radar.searchRange; bench defaults");
add("range_search(1e4, 2, 4, 0.5, 0.3, 800, 5, 20)", R.searchRange(1e4, 2, 4, 0.5, 0.3, 800, 5, 20), "rel:1e-9", "Radar.searchRange");
for (const [A, Om, ts, Rr] of [[10, 1, 10, 100e3], [2, 0.3, 4, 40e3], [50, 2, 12, 300e3]]) {
  const v = R.searchViaTrack(1000, A, ts, 1, Om, Rr, 500, 10);
  add(`search_track(1000, ${A}, ${ts}, 1, ${Om}, ${Rr}, 500, 10).snr`, v.snr, "rel:1e-9", "Radar.searchViaTrack = snrSearch " + R.snrSearch(1000, A, ts, 1, Om, Rr, 500, 10).toPrecision(6));
}
const vt = R.searchViaTrack(1000, 10, 10, 1, 1, 100e3, 500, 10);
add("search_track(1000, 10, 10, 1, 1, 100e3, 500, 10).n", vt.n, "rel:1e-9", "Radar.searchViaTrack pulses per beam");
add("search_track(1000, 10, 10, 1, 1, 100e3, 500, 10).beams", vt.beams, "rel:1e-9", "Radar.searchViaTrack beam positions");

/* ---------- noise temperature ---------- */
head("radar: system noise temperature (radar/radar.js section 4)");
const st = R.systemTemp(150, R.undb(1), 290, R.undb(3));
add("tsys().ts", st.Ts, "rel:1e-9", "Radar.systemTemp; bench 588.4 K");
add("tsys().tr", st.Tr, "rel:1e-9", "Radar.systemTemp line");
add("tsys().te", st.Te, "rel:1e-9", "Radar.systemTemp receiver");
add("tsys(50, undb(2), 300, undb(1.5)).ts", R.systemTemp(50, R.undb(2), 300, R.undb(1.5)).Ts, "rel:1e-9", "Radar.systemTemp");
err("tsys(150, 0.5)", "a line loss ratio below 1 is a gain");
add("tsys_trace()", R.noiseTraced(150, R.undb(1), 290, R.undb(3)), "rel:1e-9", "Radar.noiseTraced = 588.4 K");
add("tsys_trace(50, undb(2), 300, undb(1.5))", R.noiseTraced(50, R.undb(2), 300, R.undb(1.5)), "rel:1e-9", "Radar.noiseTraced");

/* ---------- detection ---------- */
head("radar: detection (radar/radar.js section 3)");
add("det_thresh()", R.threshold(1, 1e-6), "rel:1e-9", "Radar.threshold; closed form -ln 1e-6 = " + (-Math.log(1e-6)));
add("det_thresh(21, 1e-6)", R.threshold(21, 1e-6), "rel:1e-9", "Radar.threshold");
add("det_thresh(64, 1e-8)", R.threshold(64, 1e-8), "rel:1e-9", "Radar.threshold");
err("det_thresh(0)", "n must be >= 1");
add("pfa_thr(10, 3)", R.gammaQ(3, 10), "rel:1e-9", "Radar.gammaQ");
add("pfa_thr(det_thresh(21, 1e-6), 21)", R.gammaQ(21, R.threshold(21, 1e-6)), "rel:1e-6", "Radar.gammaQ at its threshold = 1e-6");
add("pd_steady(undb(13))", R.pdSteady(R.undb(13), 1, 1e-6), "rel:1e-9", "Radar.pdSteady; test 0.8744");
add("pd_steady(undb(3.03), 21)", R.pdSteady(R.undb(3.03), 21, 1e-6), "rel:1e-9", "Radar.pdSteady; ~0.9 at 21 pulses");
add("pd_steady(undb(5), 8, 1e-4)", R.pdSteady(R.undb(5), 8, 1e-4), "rel:1e-9", "Radar.pdSteady");
add("pd_steady(0, 4, 1e-3)", R.pdSteady(0, 4, 1e-3), "rel:1e-9", "Radar.pdSteady; no signal gives Pfa");
add("pd_steady(20, 500)", R.pdSteady(20, 500, 1e-6), "abs:1e-12", "Radar.pdSteady; normal limit branch");
err("pd_steady(10, 2.5)", "pulses must be whole");
add("pd_sw1(undb(13))", R.pdSwerling1(R.undb(13), 1e-6), "rel:1e-9", "Radar.pdSwerling1; test 0.5172");
add("pd_sw1(100, 1e-8)", R.pdSwerling1(100, 1e-8), "rel:1e-9", "Radar.pdSwerling1");
err("pd_sw1(10, 2)", "pfa must be below 1");
add("snr_req_db()", R.snrRequired(0.9, 1, 1e-6), "rel:1e-9", "Radar.snrRequired; test 13.18 dB");
add("snr_req_db(0.9, 21)", R.snrRequired(0.9, 21, 1e-6), "rel:1e-9", "Radar.snrRequired; test 3.03 dB");
add("snr_req_db(0.5, 10, 1e-4)", R.snrRequired(0.5, 10, 1e-4), "rel:1e-9", "Radar.snrRequired");
add("albersheim()", R.albersheim(0.9, 1e-6, 1), "rel:1e-9", "Radar.albersheim; test 13.11 dB");
add("albersheim(0.9, 1e-6, 21)", R.albersheim(0.9, 1e-6, 21), "rel:1e-9", "Radar.albersheim");
err("albersheim(1)", "pd must be below 1");
const nc21 = R.snrRequired(0.9, 1, 1e-6) - R.snrRequired(0.9, 21, 1e-6);
add("intgain().nc", nc21, "rel:1e-9", "Radar.snrRequired difference; README 10.2 dB for 21 pulses");
add("intgain().coh", R.db(21), "rel:1e-9", "10 log10 21 = 13.2 dB coherent");
add("intgain(8, 0.5, 1e-4).loss", R.db(8) - (R.snrRequired(0.5, 1, 1e-4) - R.snrRequired(0.5, 8, 1e-4)), "rel:1e-9", "Radar.snrRequired: integration loss");

// Marcum Q_M(a, b) = integral_b^inf x (x/a)^(M-1) exp(-(x^2 + a^2)/2) I_(M-1)(a x) dx by composite Simpson.
function besselI(nu, z) { let t = Math.pow(z / 2, nu); for (let k = 1; k <= nu; k++) t /= k; let s = t; for (let k = 1; k < 500; k++) { t *= (z * z / 4) / (k * (k + nu)); s += t; if (t < 1e-18 * s) break; } return s; }
function marcumIntegral(M, a, b) {
  const hi = b + a + 40, n = 400000, h = (hi - b) / n;
  const f = x => Math.pow(x / a, M - 1) * x * Math.exp(-(x * x + a * a) / 2) * besselI(M - 1, a * x);   // a > 0
  let s = f(b) + f(hi); for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * f(b + i * h); return s * h / 3;
}
add("marcumq(sqrt(2*undb(13)), sqrt(-2*ln(1e-6)))", marcumIntegral(1, Math.sqrt(2 * R.undb(13)), Math.sqrt(-2 * Math.log(1e-6))), "abs:1e-9", "Marcum Q1 integral; = steady Pd 0.8744");
add("marcumq(2, 3)", marcumIntegral(1, 2, 3), "abs:1e-9", "Marcum Q1 integral");
add("marcumq(1.5, 2.5, 2)", marcumIntegral(2, 1.5, 2.5), "abs:1e-9", "Marcum Q2 integral");
add("marcumq(0, 2)", Math.exp(-2), "rel:1e-12", "closed form Q1(0, b) = exp(-b^2/2)");
add("marcumq(3, 0)", 1, "abs:1e-12", "Q(a, 0) = 1");
err("marcumq(-1, 2)", "a must be >= 0");
add("pd_mc(undb(13))", R.pdMonteCarlo(R.undb(13), 1, 1e-6, 20000), "abs:1e-4", "Radar.pdMonteCarlo seed 7 (same mulberry32 stream); ~0.8744");
add("pd_mc(undb(13), 1, 1e-6, 20000, 1)", R.pdMonteCarlo(R.undb(13), 1, 1e-6, 20000, true), "abs:1e-4", "Radar.pdMonteCarlo Swerling 1; ~0.5172");
add("pd_mc(undb(4), 8, 1e-5, 5000, 0, 3)", R.pdMonteCarlo(R.undb(4), 8, 1e-5, 5000, false, 3), "abs:3e-4", "Radar.pdMonteCarlo 8 pulses");
err("pd_mc(10, 100, 1e-6, 20000)", "trials x n above the 1e6 limit");

/* ---------- GERT ---------- */
head("radar: GERT time to first track (radar/radar.js section 5); node numbers from 1 here, from 0 in the module");
const T = 60 / 12.8;
const two3 = R.rule("2 of 3", 0.7, 1), g23 = R.gert(two3.n, two3.branches, 0, two3.end);
const mat = br => "[" + br.map(([i, j, p, t]) => `${i + 1}, ${j + 1}, ${p}, ${t}`).join("; ") + "]";
add(`gert(${mat(two3.branches)}).mean`, g23.mean, "rel:1e-9", "Radar.gert; Python 2.998 scans for 2 of 3 at 0.7");
add(`gert(${mat(two3.branches)}).var`, g23.variance, "rel:1e-9", "Radar.gert; Python 2.128");
add(`gert(${mat(two3.branches)}).p`, g23.P, "rel:1e-12", "Radar.gert completes");
for (const [p, name, rule, mean] of [[0.5, "2 in a row", 1, 28.12], [0.8, "2 in a row", 1, 13.18], [0.6, "3 in a row", 3, 42.53]]) {
  const net = R.rule(name, p, T), g = R.gert(net.n, net.branches, 0, net.end), k = rule === 1 ? 2 : 3;
  add(`gert(trk_net(${p}, 60/12.8, ${rule})).mean`, g.mean, "rel:1e-9", `Radar.gert ${name}; closed form T(1-p^k)/((1-p)p^k) = ${T * (1 - p ** k) / ((1 - p) * p ** k)}; test ${mean}`);
  add(`ttrack(${p}, 60/12.8, ${rule}).mean`, g.mean, "rel:1e-9", `Radar.gert ${name}`);
  add(`ttrack(${p}, 60/12.8, ${rule}).var`, g.variance, "rel:1e-9", `Radar.gert ${name}`);
}
const p90 = (net) => { const pmf = R.scanPmf(net); let cum = 0; for (let k = 1; k < pmf.length; k++) { cum += pmf[k]; if (cum >= 0.9) return k; } return NaN; };
add("ttrack().p90", p90(R.rule("2 in a row", 0.7, 4.69)), "abs:0", "Radar.scanPmf; scans until 90% of tracks form (bench defaults)");
add("ttrack(0.4, 4.69, 3).p90", p90(R.rule("3 in a row", 0.4, 4.69)), "abs:0", "Radar.scanPmf");
add("ttrack(0.7, 1, 2).sd", Math.sqrt(g23.variance), "rel:1e-9", "Radar.gert 2 of 3");
err("ttrack(0.7, 4.69, 4)", "rule must be 1, 2 or 3");
add("trk_net(0.7, 4.69, 2)(4, 2)", two3.branches[3][1] + 1, "abs:0", "Radar.rule 2 of 3: one hit then a miss goes to node 3");
add("size(trk_net(0.7, 4.69, 3)).rows", R.rule("3 in a row", 0.7, 4.69).branches.length, "abs:0", "Radar.rule 3 in a row branch count");
err("gert([1, 2, 0.5; 2, 3, 0.5])", "branch rows need four columns");
err("gert([1, 2, 1, 1; 2, 1, 1, 1; 3, 3, 1, 1])", "the finish node cannot be reached");

/* ---------- radar basics ---------- */
head("signal: radar basics (radar/basics.js); angles in radians here, degrees in the bench");
const m1 = B.matchedFilterRange(12345), m2 = B.matchedFilterRange(12350, { noise: 0.5 });
add("mfrange(12345).est", m1.estimate, "rel:1e-12", "RadarBasics.matchedFilterRange seed 14; within one sample of 12345 m");
add("mfrange(12350, 20e6, 1e-6, 0.5).est", m2.estimate, "rel:1e-12", "RadarBasics.matchedFilterRange noise 0.5");
add("mfrange(12345).sample", m1.sample, "rel:1e-12", "c / 2 fs");
err("mfrange(-5)", "range must be >= 0");
add("chirp_res(5e6)", B.compressTwo(45).resolution, "rel:1e-12", "RadarBasics.compressTwo c/2B = 29.98 m");
add("chirp_res(150e6)", C / 3e8, "rel:1e-12", "c / 2B");
add("chirp2(45).resolved", B.compressTwo(45).resolved, "abs:0", "RadarBasics.compressTwo: 45 m resolved");
add("chirp2(15).resolved", B.compressTwo(15).resolved, "abs:0", "RadarBasics.compressTwo: 15 m merged");
add("chirp2(25, 10e6).resolved", B.compressTwo(25, { B: 10e6 }).resolved, "abs:0", "RadarBasics.compressTwo at 10 MHz");
add("chirp2(45).ratio", B.compressTwo(45).ratio, "rel:1e-12", "compression ratio BT = 50");
add("r_unamb(2e3)", B.unambRange(2e3), "rel:1e-12", "RadarBasics.unambRange 74948 m");
add("r_unamb(300)", B.unambRange(300), "rel:1e-12", "RadarBasics.unambRange");
err("r_unamb(0)", "PRF must be > 0");
for (const prf of [500, 1e4, 2e5]) add(`r_unamb(${prf})*v_unamb(${prf}, 0.03)`, C * 0.03 / 8, "rel:1e-12", "R_u v_u = c lam / 8 at every PRF");
add("v_unamb(30e3, 0.03)", B.unambVel(30e3, 0.03), "rel:1e-12", "RadarBasics.unambVel");
add("rfold(100e3, 2e3)", B.foldRange(100e3, 2e3), "rel:1e-12", "RadarBasics.foldRange 25052 m");
add("rfold(50e3, 2e3)", B.foldRange(50e3, 2e3), "rel:1e-12", "not folded");
for (const [v, prf] of [[150, 30e3], [-60, 30e3], [300, 10e3], [100, 10e3]]) add(`valias(${v}, ${prf}, 0.03)`, B.aliasVel(v, prf, 0.03), "abs:1e-9", "RadarBasics.aliasVel");
add("dop_shift(150, c0/10e9)", 2 * 150 / (C / 10e9), "rel:1e-12", "basics bench readout 2v/lam");
add("dop_shift(-20, 0.1)", -400, "rel:1e-12", "2v/lam");
for (const [v, prf] of [[150, 30e3], [300, 10e3], [100, 10e3]]) {
  const d = B.dopplerEstimate(v, prf, 0.03);
  add(`dop_dft(${v}, ${prf}, 0.03).est`, d.estimate, "abs:1e-9", "RadarBasics.dopplerEstimate seed 5; aliased " + B.aliasVel(v, prf, 0.03));
}
add("dop_dft(150, 30e3, 0.03).bin", B.dopplerEstimate(150, 30e3, 0.03).binVel, "rel:1e-12", "one DFT bin in m/s");
add("prf_trade(10e3, c0/10e9).prod", B.unambRange(10e3) * B.unambVel(10e3, C / 10e9), "rel:1e-12", "basics bench: product c lam / 8");
add("prf_trade(10e3, c0/10e9).regime", 2, "abs:0", "basics bench: 10 kHz is medium PRF");
add("prf_trade(1e3, 0.1).regime", 1, "abs:0", "basics bench: low PRF");
add("prf_trade(1e3, 0.1).ru", B.unambRange(1e3), "rel:1e-12", "RadarBasics.unambRange");
for (const v of [3, 7.5, 22.5]) {
  const m = B.mti(v, 1e3, 0.03);
  add(`mti(${v}).meas`, m.measured, "rel:1e-9", "RadarBasics.mti simulated canceller");
  add(`mti(${v}).pred`, m.predicted, "rel:1e-9", "RadarBasics.mti 4 sin^2(pi fd/PRF)");
}
add("mti(15).meas", B.mti(15, 1e3, 0.03).measured, "abs:1e-12", "RadarBasics.mti: 15 m/s is the first blind speed");
add("mti(3).vblind", B.mti(3, 1e3, 0.03).blindSpeed, "rel:1e-12", "lam PRF / 2");
add("mti(3).gain_db", 10 * Math.log10(Math.max(B.mti(3, 1e3, 0.03).measured, 1e-12) / 2), "rel:1e-9", "basics bench readout");
const cb = B.clutterBand(200, 0.03);
add("clutter(200, 0.03).hi", cb.hi, "rel:1e-12", "RadarBasics.clutterBand seed 8");
add("clutter(200, 0.03).lo", cb.lo, "rel:1e-12", "RadarBasics.clutterBand seed 8");
add("clutter(200, 0.03).edge", cb.edge, "rel:1e-12", "2V/lam");
const bw = B.apertureBeamwidth(1, 0.03);
add("beam3db(1, 0.03).num", rad(bw.numeric), "rel:1e-9", "RadarBasics.apertureBeamwidth 1.523 deg");
add("beam3db(1, 0.03).rule", rad(bw.rule), "rel:1e-9", "0.886 lam / D");
add("beam3db(0.05, 0.03).num", rad(B.apertureBeamwidth(0.05, 0.03).numeric), "rel:1e-9", "RadarBasics.apertureBeamwidth: aperture under 2 wavelengths");
add("afactor(32, 0.5, rad(40), sin(rad(40)))", B.arrayFactor(32, 0.5, 40, Math.sin(rad(40))), "rel:1e-12", "RadarBasics.arrayFactor at the steered peak");
add("afactor(16, 0.5, 0, 0.1)", B.arrayFactor(16, 0.5, 0, 0.1), "abs:1e-12", "RadarBasics.arrayFactor off the beam");
for (const [n, d, s] of [[32, 0.5, 0], [32, 0.5, 40]]) {
  const a = B.arrayScan(n, d, s);
  add(`ascan(${n}, ${d}, rad(${s})).peak`, rad(a.peakDeg), "abs:1e-12", "RadarBasics.arrayScan beam peak");
  add(`ascan(${n}, ${d}, rad(${s})).width`, rad(a.width), "rel:1e-9", "RadarBasics.arrayScan half-power width");
}
add("ascan(32, 0.58, rad(40)).grating", B.arrayScan(32, 0.58, 40).grating, "rel:1e-9", "RadarBasics.arrayScan: no grating lobe at 0.58");
add("ascan(32, 0.66, rad(40)).grating", B.arrayScan(32, 0.66, 40).grating, "rel:1e-9", "RadarBasics.arrayScan: grating lobe at 0.66");
add("ascan(32, 0.5, rad(40)).limit", B.arrayScan(32, 0.5, 40).limit, "rel:1e-12", "1/(1 + |sin theta0|)");
add("ascan(32, 0.5, rad(40)).phase", rad(360 * 0.5 * Math.sin(rad(40))), "rel:1e-12", "basics bench phase step 360 d sin theta0");
err("ascan(0, 0.5, 0)", "needs at least one element");
add("r_horizon(20, 3000)", B.horizonExact(20, 3000), "rel:1e-12", "RadarBasics.horizonExact 244.2 km");
add("r_horizon(10)", B.horizonExact(10, 0), "rel:1e-12", "RadarBasics.horizonExact");
err("r_horizon(-1, 3)", "height below the ground");
add("r_horizon_rule(20, 3000)", B.horizonRule(20, 3000), "rel:1e-12", "RadarBasics.horizonRule 4.12 rule");
add("r_horizon_rule(100, 9000)", B.horizonRule(100, 9000), "rel:1e-12", "RadarBasics.horizonRule");
add("frame_time(rad(120), rad(3), 4, 0.02)", B.frameTime(120, 3, 4, 0.02), "rel:1e-12", "RadarBasics.frameTime bench defaults 3.2 s");
add("frame_time(rad(100), rad(3), 2, 0.01)", B.frameTime(100, 3, 2, 0.01), "rel:1e-12", "RadarBasics.frameTime 34 positions");
err("frame_time(1, 0, 1, 0.01)", "beamwidth must be > 0");
const t1 = Math.atan2(20e3 * Math.sin(1.2), 20e3 * Math.cos(1.2)), t2 = Math.atan2(20e3 * Math.sin(1.2), 20e3 * Math.cos(1.2) - 2e3), tri = B.triangulate(2e3, t1, t2);
add(`tri_bear(2e3, ${t1}, ${t2}).x`, tri[0], "rel:1e-9", "RadarBasics.triangulate");
add(`tri_bear(2e3, ${t1}, ${t2}).r`, Math.hypot(...tri), "rel:1e-9", "RadarBasics.triangulate: back to 20 km");
err("tri_bear(2e3, 1, 1)", "parallel bearings");
add("tri_spread(20e3, 2e3, 1e-3)", B.triangulationSpread(20e3, 2e3, 1e-3), "rel:1e-9", "RadarBasics.triangulationSpread seed 9");
add("tri_spread(40e3, 2e3, 1e-3)", B.triangulationSpread(40e3, 2e3, 1e-3), "rel:1e-9", "RadarBasics.triangulationSpread: about 4x at twice the range");
const bandIdx = f => { const n = B.BANDS.findIndex(([nm]) => nm === B.bandOf(f)); return n + 1; };
add("rband(10e9).band", bandIdx(10), "abs:0", "RadarBasics.bandOf X");
add("rband(3e9).band", bandIdx(3), "abs:0", "RadarBasics.bandOf S");
add("rband(35e9).band", bandIdx(35), "abs:0", "RadarBasics.bandOf Ka");
add("rband(50e9).band", bandIdx(50), "abs:0", "RadarBasics.bandOf outside L-Ka");
add("rband(10e9).lam", C / 10e9, "rel:1e-12", "c / f");

/* ---------- lens ---------- */
head("lens: Luneburg lens (lens/lens.js); lengths in lens radii unless SI");
add("lune_eps(0.5)", L.epsAt(0.5), "rel:1e-12", "Lens.epsAt 2 - r^2");
add("lune_eps(0.95, 1.2)", L.epsAt(0.95, 1.2), "rel:1e-12", "Lens.epsAt clipped to the printable floor");
for (const b of [0.3, 0.8]) {
  const t = L.trace(b);
  add(`lune_trace(${b}).x`, t.x, "rel:1e-9", "Lens.trace exit; analytic focus (1, 0)");
  add(`lune_trace(${b}).y`, t.y, "abs:1e-9", "Lens.trace exit");
  add(`lune_trace(${b}).opl`, t.opl, "rel:1e-9", "Lens.trace optical path; analytic pi/2 + 1");
}
const tAx = L.trace(1e-6, { dt: 1e-4, tan: e => L.tanEff(e, 0.0043, 2.8) });
add("lune_trace(1e-6, 1, 2, 1e-4, 0.0043, 2.8).loss", tAx.loss, "rel:1e-9", "Lens.trace loss integral; axial quadrature " + L.axialQuadrature(0.0043, 2.8));
add("lune_trace(0.5, 1.2, 2).angle", L.trace(0.5, { emin: 1.2 }).angle, "rel:1e-9", "Lens.trace with a printable floor 1.2");
err("lune_trace(1)", "b must be inside the lens");
add("lune_spot(1.2, 2).rms", L.spot({ emin: 1.2 }).rms, "rel:1e-9", "Lens.spot; Python 0.059 rad");
add("lune_spot(1.5, 2).rms", L.spot({ emin: 1.5 }).rms, "rel:1e-9", "Lens.spot; Python 0.205 rad");
const sp = L.spot({ dt: 2e-4 }, 10);
add("lune_spot(1, 2, 10, 2e-4).worst", sp.worstExit, "rel:1e-6", "Lens.spot ideal profile: all rays to the focus");
add("lune_spot(1, 2, 10, 2e-4).opl", sp.oplSpread, "rel:1e-6", "Lens.spot ideal profile: equal optical paths");
const mixNames = ["linear (parallel)", "Maxwell Garnett", "Lichtenecker", "series (lower bound)"];
mixNames.forEach((nm, i) => add(`mix_eps(0.5, 2.8, ${i + 1})`, L.mixing[nm](0.5, 2.8), "rel:1e-12", `Lens.mixing ${nm}`));
err("mix_eps(1.5)", "fill fraction above 1");
mixNames.slice(0, 3).forEach((nm, i) => add(`fill_for(2, 2.8, ${i + 1})`, L.fillFor(2, nm, 2.8), "rel:1e-12", `Lens.fillFor ${nm}; test ${[0.556, 0.667, 0.673][i]}`));
add("fill_for(1.75, 2.9, 4)", L.fillFor(1.75, "series (lower bound)", 2.9), "rel:1e-12", "Lens.fillFor series");
err("fill_for(3, 2.8)", "eps beyond the solid's permittivity");
const cf = L.cellFit();
add("cell_fit().slope", cf.slope, "rel:1e-12", "Lens.cellFit; test -0.980");
add("cell_fit().prod", cf.product * 1e6, "rel:1e-12", "Lens.cellFit; 249 mm GHz in m Hz");
add("cell_fit(10e-3).fmax", cf.product / 10 * 1e9, "rel:1e-12", "Lens.cellFit product / cell");
add("cell_fit().cell_lam", cf.product * 1e6 / C, "rel:1e-12", "lens bench: limit about 0.83 free-space wavelengths");
add("tan_eff(2, 0.0043, 2.8)", L.tanEff(2, 0.0043, 2.8), "rel:1e-12", "Lens.tanEff");
add("tan_eff(1.5, 0.039, 2.9)", L.tanEff(1.5, 0.039, 2.9), "rel:1e-12", "Lens.tanEff");
for (const f of [5e9, 17e9, 30e9]) add(`radix_df(${f})`, L.MATERIALS.Radix.df(f), "rel:1e-12", "Lens.MATERIALS Radix");
add("lens_loss_db(0.05, 30e9, 0.0046, 2.8)", L.apertureLossDb(0.05, 30e9, 0.0046, 2.8), "rel:1e-9", "Lens.apertureLossDb; Radix 10 cm, 0.67 dB");
add("lens_loss_db(0.05, 30e9, 0.039, 2.9)", L.apertureLossDb(0.05, 30e9, 0.039, 2.9), "rel:1e-9", "Lens.apertureLossDb; DLP resin 10 cm, 5.25 dB");
add("lens_loss_db(0.05, 30e9, radix_df(30e9), 2.8, 1.2)", L.apertureLossDb(0.05, 30e9, 0.0046, 2.8, 1.2), "rel:1e-9", "Lens.apertureLossDb; printable floor 1.2, 0.69 dB");
err("lens_loss_db(0, 30e9, 0.0046)", "radius must be > 0");
add("ax_loss(0.0043)", L.axialQuadrature(0.0043, 2.8), "rel:1e-12", "Lens.axialQuadrature");
add("ax_loss(0.039, 2.9, 1000)", L.axialQuadrature(0.039, 2.9, 1000), "rel:1e-12", "Lens.axialQuadrature");
add("lens_dir_db(0.05, 30e9)", L.directivityDb(0.05, 30e9), "rel:1e-12", "Lens.directivityDb 4 pi A / lam^2");
add("lens_dir_db(0.1, 10e9)", L.directivityDb(0.1, 10e9), "rel:1e-12", "Lens.directivityDb");
add("rcs_lune(0.05, c0/9.4e9)", L.rcsLuneburg(0.05, C / 9.4e9), "rel:1e-12", "Lens.rcsLuneburg; test 0.76 m^2");
add("rcs_lune(0.1, 0.03)", L.rcsLuneburg(0.1, 0.03), "rel:1e-12", "Lens.rcsLuneburg");
const sig = L.rcsLuneburg(0.05, C / 9.4e9);
add("refl_range(0.05, c0/9.4e9)", Math.pow((0.1 + sig) / 0.1, 0.25), "rel:1e-12", "lens bench section 5 range ratio, boat 0.1 m^2");
add("refl_range(0.02, 0.1, 1)", Math.pow((1 + L.rcsLuneburg(0.02, 0.1)) / 1, 0.25), "rel:1e-12", "lens bench section 5 range ratio");
err("refl_range(0.05, 0.03, 0)", "the target's own cross section must be > 0");

process.stdout.write("# Generated by Tests/reference/gen_radar.mjs from Where It Isn't's radar/radar.js, radar/basics.js and lens/lens.js.\n" + out.join("\n") + "\n");
