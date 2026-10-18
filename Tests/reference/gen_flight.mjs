// Reference cases for the flight library (Source/MechanicalCalculator/Calc/LibFlight.cpp).
//   node Tests/reference/gen_flight.mjs /home/user/where-it-isnt > Tests/reference/flight.tsv
// Every expected value is computed by Where It Isn't's own modules: flight-model/equations.js (AeroEq),
// flight-model/model.js (FlightModel), flight-model/planform.js (Planform) and game/game.js (ParGame, with
// radar/radar.js). Where a calculator function returns a closed form that the module only integrates
// numerically, the case compares against the module's numerical value with a tolerance that covers the
// integration error, and says so in its provenance.
import path from "node:path";
import { pathToFileURL } from "node:url";

const wii = process.argv[2];
if (!wii) { console.error("usage: node gen_flight.mjs /path/to/where-it-isnt"); process.exit(2); }
const load = rel => import(pathToFileURL(path.join(wii, rel)).href);
await load("flight-model/equations.js");
await load("radar/radar.js");
await load("game/game.js");
await load("flight-model/model.js");
await load("flight-model/planform.js");
const Q = globalThis.AeroEq, F = globalThis.FlightModel, PF = globalThis.Planform, G = globalThis.ParGame;

const lines = [];
const num = x => Number.isNaN(x) ? "nan" : x === Infinity ? "inf" : x === -Infinity ? "-inf" : String(x);
const lit = x => String(x);                       // a number as typed into the calculator
const call = (name, ...args) => `${name}(${args.map(lit).join(", ")})`;
function out(expr, value, tol, prov) {
  if (typeof value === "boolean") value = value ? 1 : 0;
  if (value !== "error" && typeof value !== "number") throw new Error(`${expr}: bad expected ${value}`);
  lines.push([expr, value === "error" ? "error" : num(value), tol, prov].join("\t"));
}
const err = (expr, prov) => out(expr, "error", "abs:0", prov);
const section = t => lines.push("# " + t);
const R9 = "rel:1e-9", R12 = "rel:1e-12";
const mat = rows => "[" + rows.map(r => r.map(lit).join(", ")).join("; ") + "]";

/* ------------------------------------------------------------------ aero: the equation sheet */
section("aero: flight-model/equations.js (AeroEq), defaults from equations.html");
{
  const T = Q.thinAirfoil(0.02, 0.4, Q.rad(4));
  out("thinairfoil().aL0", T.aL0, R9, "AeroEq.thinAirfoil(0.02, 0.4, rad 4); NACA 24xx alpha_L0 = -2.08 deg");
  out("thinairfoil().cl", T.cl, R9, "AeroEq.thinAirfoil default card");
  out("thinairfoil().cm4", T.cm4, R9, "AeroEq.thinAirfoil default card");
  const T2 = Q.thinAirfoil(0.04, 0.3, Q.rad(2));
  out(call("thinairfoil", 0.04, 0.3, Q.rad(2)) + ".cl", T2.cl, R9, "AeroEq.thinAirfoil(0.04, 0.3, rad 2)");
  out(call("thinairfoil", 0.04, 0.3, Q.rad(2)) + ".A1", T2.A1, R9, "AeroEq.thinAirfoil A1");
  out(call("thinairfoil", 0.04, 0.3, Q.rad(2)) + ".A2", T2.A2, R9, "AeroEq.thinAirfoil A2");
  out("thinairfoil(0, 0.4, 0.1).cl", Q.thinAirfoil(0, 0.4, 0.1).cl, R12, "symmetric section: 2 pi alpha");
  err("thinairfoil(0.02, 1.2, 0.1)", "camber position outside the chord");
  const V = Q.vortexPanels(0.02, 0.4, Q.rad(4));
  out("vortexpanel().cl", V.cl, "rel:1e-8", "AeroEq.vortexPanels(0.02, 0.4, rad 4), 200 panels");
  out("vortexpanel().cm4", V.cm4, "rel:1e-8", "AeroEq.vortexPanels default");
  out(call("vortexpanel", 0.04, 0.3, Q.rad(2), 100) + ".cl", Q.vortexPanels(0.04, 0.3, Q.rad(2), 100).cl, "rel:1e-8", "AeroEq.vortexPanels, 100 panels");
  err("vortexpanel(0.02, 0.4, 0.1, 1)", "needs at least two panels");

  out("pg()", Q.prandtlGlauert(0.5, 0.6), R12, "AeroEq.prandtlGlauert(0.5, 0.6), page default");
  out("pg(1, 0.3)", Q.prandtlGlauert(1, 0.3), R12, "AeroEq.prandtlGlauert");
  err("pg(0.5, 1)", "Prandtl-Glauert is singular at Mach 1");
  out("karmantsien()", Q.karmanTsien(-0.5, 0.6), R12, "AeroEq.karmanTsien(-0.5, 0.6), page default");
  out("karmantsien(-0.05, 0.6)", Q.karmanTsien(-0.05, 0.6), R12, "AeroEq.karmanTsien, a tenth of the coefficient");
  err("karmantsien(-0.5, 1.2)", "subsonic only");
  out("ackeret()", Q.ackeretCl(Q.rad(4), 2), R12, "AeroEq.ackeretCl(rad 4, 2)");
  out("ackeret(0.05, 3)", Q.ackeretCl(0.05, 3), R12, "AeroEq.ackeretCl");
  err("ackeret(0.1, 0.9)", "supersonic only");

  const P = Q.polar(0.025, 8, 0.8);
  out("dragpolar().LDmax", P.LDmax, R12, "AeroEq.polar(0.025, 8, 0.8), page default");
  out("dragpolar().CLstar", P.CLstar, R12, "AeroEq.polar");
  out("dragpolar().LD", P.LD(0.5), R12, "AeroEq.polar(...).LD(0.5)");
  out("dragpolar().k", P.k, R12, "AeroEq.polar");
  out("dragpolar(0.02, 10, 0.9, 0.8).LD", Q.polar(0.02, 10, 0.9).LD(0.8), R12, "AeroEq.polar");
  err("dragpolar(0.02, 0, 0.9)", "aspect ratio must be positive");

  const E = Q.liftingLineE(Q.taperChord(0.4), 8);
  out("spaneff().e", E.e, R9, "AeroEq.liftingLineE(taperChord(0.4), 8), page default");
  out("spaneff().delta", E.delta, "rel:1e-8", "AeroEq.liftingLineE");
  out("spaneff(6, 1).delta", Q.liftingLineE(Q.taperChord(1), 6).delta, "rel:1e-8", "AeroEq.liftingLineE, rectangular AR 6 (textbook about 0.05)");
  out("spaneff(8, 0.4, 20).e", Q.liftingLineE(Q.taperChord(0.4), 8, 20).e, R9, "AeroEq.liftingLineE, 20 terms as in the page's taper sweep");
  err("spaneff(8, 0)", "taper ratio must be positive");
  out("spaneff_ell().e", Q.liftingLineE(Q.ellipticChord, 8).e, R12, "AeroEq.liftingLineE(ellipticChord, 8): e = 1");
  out("spaneff_ell(5, 20).delta", Q.liftingLineE(Q.ellipticChord, 5, 20).delta, "abs:1e-12", "AeroEq.liftingLineE elliptic: delta = 0");

  const sh = Q.waveDq(Q.searsHaack(0.05));
  out("wavedrag().fourier", sh, "rel:1e-7", "AeroEq.waveDq(searsHaack(0.05))");
  out("wavedrag().closed", Q.searsHaackDq(0.05), R12, "AeroEq.searsHaackDq(0.05)");
  out("wavedrag().ratio", sh / Q.searsHaackDq(0.05), "rel:1e-7", "Fourier over closed form (page check < 2e-3)");
  out("wavedrag(0.1, 2).fourier", Q.waveDq(Q.searsHaack(0.1, 2), 2), "rel:1e-7", "AeroEq.waveDq(searsHaack(0.1, 2), 2)");
  err("wavedrag(0, 1)", "radius must be positive");
  const AR = (wing, rule) => { const r = 0.05, SH = Q.searsHaack(r), A0 = Math.PI * r * r, wf = x => wing * A0 * Math.max(0, 1 - ((x - 0.55) / 0.12) ** 2) ** 2;
    const dSH = Q.waveDq(SH), dW = Q.waveDq(x => SH(x) + wf(x)), dT = Q.waveDq(x => SH(x) + wf(x) - rule * wf(x)); return { fus: dSH, wing: dW / dSH, ruled: dT / dSH }; };
  const a1 = AR(0.4, 0), a2 = AR(0.6, 0.8);
  out("arearule().fus", a1.fus, "rel:1e-7", "equations.html area card via AeroEq.waveDq, default sliders");
  out("arearule().wing", a1.wing, "rel:1e-7", "equations.html area card: with the wing");
  out("arearule(0.6, 0.8).ruled", a2.ruled, "rel:1e-7", "equations.html area card: indented fuselage");
  out("arearule(0.6, 0.8).wing", a2.wing, "rel:1e-7", "equations.html area card");

  out("arearatio()", Q.areaRatio(2), R12, "AeroEq.areaRatio(2)");
  out("arearatio(0.5, 1.3)", Q.areaRatio(0.5, 1.3), R12, "AeroEq.areaRatio(0.5, 1.3)");
  err("arearatio(0)", "Mach must be positive");
  out("machfromarea(1.6875, 1)", Q.machFromArea(1.6875, true), R12, "AeroEq.machFromArea supersonic root");
  out("machfromarea(1.6875, 0)", Q.machFromArea(1.6875, false), R12, "AeroEq.machFromArea subsonic root");
  out("machfromarea(4, 1, 1.25)", Q.machFromArea(4, true, 1.25), R12, "AeroEq.machFromArea, gamma 1.25");
  err("machfromarea(0.5)", "no Mach number has A/A* below 1");

  out("climbangle()", Q.climbAngle(0.18, 12), R12, "AeroEq.climbAngle(0.18, 12), page default");
  out("climbangle(0.5, 8)", Q.climbAngle(0.5, 8), R12, "AeroEq.climbAngle");
  err("climbangle(0.05, 12)", "AeroEq.climbAngle gives NaN: below level-flight drag");
  out("rocketdv()", Q.rocketDv(300, 4000, 1000), R12, "AeroEq.rocketDv(300, 4000, 1000), page default");
  out("rocketdv(450, 10, 1)", Q.rocketDv(450, 10, 1), R12, "AeroEq.rocketDv");
  err("rocketdv(300, 0, 1)", "masses must be positive");

  const TP = Q.turboprop(5, 12, 288.15, 1400, 0.85, 0.9), TI = Q.turboprop(5, 12, 288.15, 1400, 1, 1);
  out("turboprop().P", TP.P, R12, "AeroEq.turboprop(5, 12, 288.15, 1400, 0.85, 0.9), page default");
  out("turboprop().eta", TP.eta, R12, "AeroEq.turboprop");
  out("turboprop().turbine", TP.turbine, R12, "AeroEq.turboprop");
  out("turboprop().compressor", TP.compressor, R12, "AeroEq.turboprop");
  out("turboprop(5, 12, 288.15, 1400, 1, 1).eta", TI.eta, R12, "AeroEq.turboprop, lossless: Brayton 1 - 1/tau");
  out("turboprop(5, 12, 288.15, 1400, 1, 1).ideal", 1 - 1 / TI.tau, R12, "Brayton efficiency");
  err("turboprop(5, 1)", "pressure ratio must exceed 1");
  const N1 = Q.netThrust(100, 0.02, 600, 250), N2 = Q.netThrust(50, 0.01, 800, 200, 120000, 101325, 0.3);
  out("netthrust().gross", N1.gross, R12, "AeroEq.netThrust(100, 0.02, 600, 250), page default");
  out("netthrust().ram", N1.ram, R12, "AeroEq.netThrust");
  out("netthrust().net", N1.net, R12, "AeroEq.netThrust");
  out("netthrust(100, 0.02, 600, 0).net", Q.netThrust(100, 0.02, 600, 0).net, R12, "AeroEq.netThrust at rest: net = gross");
  out("netthrust(50, 0.01, 800, 200, 120000, 101325, 0.3).net", N2.net, R12, "AeroEq.netThrust with an under-expanded nozzle");
}
{
  const L = Q.wingLoads(2.5 * 10e3, 5), N = L.y.length - 1;
  out("wingloads().V", L.Vroot, R12, "AeroEq.wingLoads(25 kN, 5 m).Vroot = pi w0 s / 4");
  out("wingloads().M", L.Mroot, R12, "AeroEq.wingLoads(25 kN, 5 m).Mroot = w0 s^2 / 3");
  out("wingloads().w0", L.w0, R12, "AeroEq.wingLoads");
  out("wingloads().M", L.M[0], "rel:1e-4", "AeroEq.wingLoads numerical double integral at the root (page check 1e-4)");
  out("wingloads(25000, 5, 2.5).V", L.V[N / 2], "rel:1e-4", "AeroEq.wingLoads numerical shear at mid-semispan (2000 trapezoids; tip singularity limits it to 1e-5)");
  out("wingloads(25000, 5, 2.5).M", L.M[N / 2], "rel:1e-4", "AeroEq.wingLoads numerical bending at mid-semispan");
  out("wingloads(25000, 5, 2.5).w", L.w[N / 2], R12, "AeroEq.wingLoads running load at mid-semispan");
  const L2 = Q.wingLoads(6 * 50e3, 12);
  out("wingloads(300000, 12).M", L2.Mroot, R12, "AeroEq.wingLoads Mroot");
  err("wingloads(25000, 5, 6)", "station beyond the tip");
  const T = 1000, r = 0.1, t = 0.002, Am = Math.PI * r * r;
  out("bredt().q", T / (2 * Am), R12, "equations.html torsion card: q = T/(2 A_m)");
  out("bredt().bredt", Q.bredtTube(T, r, t), R12, "AeroEq.bredtTube(1000, 0.1, 0.002)");
  out("bredt().exact", Q.exactTube(T, r, t), R12, "AeroEq.exactTube");
  out("bredt().open", Q.openTube(T, r, t), R12, "AeroEq.openTube");
  out("bredt().twist", T * (2 * Math.PI * r / t) / (4 * Am * Am * 27e9), R12, "equations.html torsion card, G = 27 GPa");
  out("bredt(2000, 0.05, 0.005).exact", Q.exactTube(2000, 0.05, 0.005), R12, "AeroEq.exactTube, thick wall");
  err("bredt(1000, 0.1, 0.3)", "wall thicker than the tube");
  out("buckle().P", Q.eulerLoad(70e9, 1e-8, 1, 1), R12, "AeroEq.eulerLoad(70 GPa, 1 cm^4, 1 m, 1)");
  out("buckle().fd", Q.fdBuckle(70e9, 1e-8, 1, false), R9, "AeroEq.fdBuckle pinned-pinned, 400 points");
  out("buckle(70e9, 1e-8, 1, 2).fd", Q.fdBuckle(70e9, 1e-8, 1, true), R9, "AeroEq.fdBuckle cantilever");
  out("buckle(70e9, 1e-8, 1, 0.699).P", Q.eulerLoad(70e9, 1e-8, 1, 0.699), R12, "AeroEq.eulerLoad fixed-pinned");
  out("buckle(200e9, 5e-7, 3, 0.5).fd", NaN, "abs:0", "no finite-difference check for fixed-fixed (page: NaN)");
  err("buckle(70e9, 1e-8, 0)", "length must be positive");
}
{
  const kAC = Q.polar(0.025, 8, 0.8).k;
  const C1 = Q.climbCurve(10e3, 16, 1.225, 60e3, 0.025, kAC), C2 = Q.climbCurve(12e3, 16, 0.9, 80e3, 0.025, kAC);
  out("bestclimb().Vbest", C1.Vbest, R12, "AeroEq.climbCurve(10 kN, 16, 1.225, 60 kW, 0.025, k), page default");
  out("bestclimb().rocmax", C1.rocMax, R9, "AeroEq.climbCurve");
  out("bestclimb().Vmp", C1.Vmp, R12, "AeroEq.climbCurve minimum-power speed");
  out("bestclimb(12000, 16, 0.9, 80000).rocmax", C2.rocMax, R9, "AeroEq.climbCurve");
  out("bestclimb(12000, 16, 0.9, 80000).Vbest", C2.Vbest, R12, "AeroEq.climbCurve");
  const V = 120 * 0.514444, Tt = Q.turn(V, Q.rad(18.2));
  out(call("levelturn", V, Q.rad(18.2)) + ".omega", Tt.omega, R12, "AeroEq.turn(120 kt, 18.2 deg), page default (standard rate)");
  out(call("levelturn", V, Q.rad(18.2)) + ".R", Tt.R, R12, "AeroEq.turn");
  out(call("levelturn", 100, Q.rad(60)) + ".n", Q.turn(100, Q.rad(60)).n, R12, "AeroEq.turn: n = 2 at 60 deg");
  out("levelturn(100, 0).R", Q.turn(100, 0).R, "abs:0", "AeroEq.turn: wings level, infinite radius");
  err("levelturn(100, 1.6)", "bank beyond 90 deg");
  out("cornerspeed()", Q.cornerSpeed(3.8, 10e3, 1.225, 16, 1.5), R12, "AeroEq.cornerSpeed, page default");
  out("cornerspeed(6, 50000, 1, 30, 1.8)", Q.cornerSpeed(6, 50000, 1, 30, 1.8), R12, "AeroEq.cornerSpeed");
  const g1 = Q.glide(1000, 10e3, 16, 1.225, 0.025, kAC, 0.71), g2 = Q.glide(2000, 8000, 16, 1.0, 0.025, kAC, 0.5);
  out("glide().ratio", g1.ratio, R12, "AeroEq.glide(1000, 10 kN, 16, 1.225, 0.025, k, 0.71), page default");
  out("glide().dist", g1.formula, R12, "AeroEq.glide formula h / tan(gamma)");
  out("glide().dist", g1.sim, "rel:1e-6", "AeroEq.glide stepped simulation (page check 1e-6)");
  out("glide().V", g1.V, R12, "AeroEq.glide");
  out("glide().sink", g1.V * Math.sin(g1.gam), R12, "equations.html glide card: V sin(gamma)");
  out(call("glide", 2000, 8000, 16, 1, 0.025, kAC, 0.5) + ".dist", g2.formula, R12, "AeroEq.glide");
  err("glide(1000, 10000, 16, 1.225, 0.025, 0.05, 0)", "CL must be positive");
}
{
  const S = Q.shortPeriod();
  out("shortperiod().wn", S.wn, R12, "AeroEq.shortPeriod() with SP defaults");
  out("shortperiod().zeta", S.zeta, R12, "AeroEq.shortPeriod");
  out("shortperiod().b0", S.num[1], R12, "AeroEq.shortPeriod numerator");
  out("shortperiod().a0", S.den[2], R12, "AeroEq.shortPeriod denominator");
  const a1 = S.den[1], a0 = S.den[2], disc = a1 * a1 - 4 * a0;
  out("shortperiod().Td", 2 * Math.PI / (Math.sqrt(-disc) / 2), R12, "equations.html state-space card: damped period");
  const P2 = { ...Q.SP, Za: -1, Ma: -20, Mq: -2 };
  out("shortperiod(-1, -0.2, -20, -2, -15).zeta", Q.shortPeriod(P2).zeta, R12, "AeroEq.shortPeriod, other derivatives");
  err("shortperiod(-2, -0.2, 10, -3, -15)", "statically unstable: no short-period oscillator");
  const tf = Q.tfTheta(1), ss = Q.ssTheta(1);
  out("pitchtf().re", tf[0], R12, "AeroEq.tfTheta(1)");
  out("pitchtf().im", tf[1], R12, "AeroEq.tfTheta(1)");
  out("pitchtf().ssre", ss[0], R9, "AeroEq.ssTheta(1): state space");
  out("pitchtf().ssim", ss[1], R9, "AeroEq.ssTheta(1)");
  out("pitchtf(10).mag", Math.hypot(...Q.tfTheta(10)), R12, "|AeroEq.tfTheta(10)|");
  out("pitchtf(0.3, -1, -0.2, -20, -2, -15).db", 20 * Math.log10(Math.hypot(...Q.tfTheta(0.3, P2))), R12, "AeroEq.tfTheta in dB, other derivatives");
  out("pitchtf(0.3, -1, -0.2, -20, -2, -15).ssim", Q.ssTheta(0.3, P2)[1], R9, "AeroEq.ssTheta, other derivatives");
  err("pitchtf(0)", "the free-attitude pole at s = 0");
  out("poerror()", Q.pOnlyError(1.5, 0.5), R12, "AeroEq.pOnlyError(1.5, 0.5)");
  out("poerror(3, 0.2)", Q.pOnlyError(3, 0.2), R12, "AeroEq.pOnlyError");
  const p1 = Q.pidSim(1.5, 1, 0.3, Q.rad(5), 0.5), p0 = Q.pidSim(1.5, 0, 0.3, Q.rad(5), 0.5), pd = Q.pidSim(-5, 0, 0, Q.rad(5), 0.5);
  out("pidsim().err", p1.err, "abs:1e-12", "AeroEq.pidSim(1.5, 1, 0.3, rad 5, 0.5), page default");
  out("pidsim().overshoot", p1.overshoot, R9, "AeroEq.pidSim");
  out("pidsim(1.5, 0, 0.3).err", p0.err, R9, "AeroEq.pidSim without integral action");
  out("pidsim(1.5, 0, 0.3).err", Q.pOnlyError(1.5, 0.5), "abs:1e-4", "AeroEq.pOnlyError prediction (page check 1e-4)");
  out("pidsim(-5, 0, 0).diverged", pd.diverged, "abs:0", "AeroEq.pidSim: wrong-sign gain diverges");
}

/* ------------------------------------------------------------------ atmos: model.js */
section("atmos: flight-model/model.js (FlightModel)");
{
  const A0 = F.atmosphere(0), A11 = F.atmosphere(11000);
  out("atmos(0).Tc", A0.Tc, R12, "FlightModel.atmosphere(0): 15.04 C");
  out("atmos(0).p", A0.p, R12, "FlightModel.atmosphere(0)");
  out("atmos(0).rho", A0.rho, R12, "FlightModel.atmosphere(0)");
  out("atmos(0).a", A0.a, R12, "FlightModel.atmosphere(0): about 760 mph");
  out("atmos(11000).a", A11.a, R12, "FlightModel.atmosphere(11000): about 659 mph");
  out("atmos(20000).p", F.atmosphere(20000).p, R12, "FlightModel.atmosphere, lower stratosphere");
  out("atmos(30000).T", F.atmosphere(30000).T, R12, "FlightModel.atmosphere, upper branch");
  out("atmos(-50).p", F.atmosphere(-50).p, R12, "FlightModel.atmosphere clamps h at 0");
  F.setWeather({ dT: 10 }); const W10 = F.atmosphere(5000); F.setWeather({ dT: 0 });
  out("atmos(5000, 10).rho", W10.rho, R12, "FlightModel.atmosphere with WEATHER.dT = 10");
  out("atmos(3000).icing", F.icing(3000, true), "abs:0", "FlightModel.icing(3000, in cloud)");
  out("atmos(0).icing", F.icing(0, true), "abs:0", "FlightModel.icing(0, in cloud)");
  out("liftslope().a0", F.sectionSlope(0), R12, "FlightModel.sectionSlope(0) = 2 pi");
  out("liftslope(0.5, 7).a", F.wingSlope(F.sectionSlope(0.5), 7), R12, "FlightModel.wingSlope(sectionSlope(0.5), 7)");
  out("liftslope(1).a0", F.sectionSlope(1), R12, "FlightModel.sectionSlope transonic blend");
  out("liftslope(2).a0", F.sectionSlope(2), R12, "FlightModel.sectionSlope Ackeret");
  out("liftslope(0, 1e6).a", F.wingSlope(2 * Math.PI, 1e6), R12, "FlightModel.wingSlope: infinite span recovers 2 pi");
  const ga = F.AIRCRAFT.ga, jet = F.AIRCRAFT.jet, sst = F.AIRCRAFT.sst;
  out(call("liftcoef", F.rad(4), 0.1, 1) + ".CL", F.liftCoefficient(F.rad(4), 0.1, ga).CL, R12, "FlightModel.liftCoefficient, GA");
  out(call("liftcoef", F.rad(22), 0.1, 1) + ".CL", F.liftCoefficient(F.rad(22), 0.1, ga).CL, R12, "FlightModel.liftCoefficient past stall");
  out(call("liftcoef", F.rad(22), 0.1, 1) + ".stalled", F.liftCoefficient(F.rad(22), 0.1, ga).stalled, "abs:0", "FlightModel.liftCoefficient");
  out("liftcoef(0.1, 0.8, 2).CLa", F.liftCoefficient(0.1, 0.8, jet).CLa, R12, "FlightModel.liftCoefficient, jet");
  err("liftcoef(0.1, 0.1, 4)", "aircraft 1 to 3");
  const a1 = F.aero(F.rad(4), 0, 50, 500, ga), a2 = F.aero(F.rad(3), 0, 120, 0, jet), a3 = F.aero(0.05, 0, 600, 16500, sst), a4 = F.aero(0.05, 0.02, 250, 10500, jet);
  out(call("acaero", F.rad(4), 50, 500, 1) + ".L", a1.L, R12, "FlightModel.aero(4 deg, 50 m/s, 500 m, ga)");
  out(call("acaero", F.rad(4), 50, 500, 1) + ".D", a1.D, R12, "FlightModel.aero");
  out(call("acaero", F.rad(3), 120, 0, 2) + ".CDi", a2.CDi, R12, "FlightModel.aero, jet induced drag");
  out("acaero(0.05, 600, 16500, 3).CD0", a3.CD0, R12, "FlightModel.aero, supersonic wave drag");
  out("acaero(0.05, 600, 16500, 3).K", a3.K, R12, "FlightModel.aero, supersonic drag due to lift");
  out("acaero(0.05, 250, 10500, 2, 0.02).CD", a4.CD, R12, "FlightModel.aero with sideslip");
  err("acaero(0.1, 0, 0, 1)", "airspeed must be positive");
  const at = h => F.atmosphere(h);
  out("acthrust(0.65, 55, 1000, 1)", F.thrust(0.65, 55, at(1000), 55 / at(1000).a, ga), R12, "FlightModel.thrust, propeller");
  out("acthrust(0.6, 250, 10500, 2)", F.thrust(0.6, 250, at(10500), 250 / at(10500).a, jet), R12, "FlightModel.thrust, jet");
  out("acthrust(1, 0, 0, 1)", F.thrust(1, 0, at(0), 0, ga), R12, "FlightModel.thrust, static");
  out("acthrust(0.9, 600, 16500, 3)", F.thrust(0.9, 600, at(16500), 600 / at(16500).a, sst), R12, "FlightModel.thrust, supersonic");
  const t1 = F.trim(ga, 55, 1000), t2 = F.trim(jet, 250, 10500), aS = at(16500).a, t3 = F.trim(sst, 2 * aS, 16500), t4 = F.trim(ga, 55, 1000, 0, F.rad(30)), t5 = F.trim(ga, 120, 1000);
  out("actrim(55, 1000, 1).throttle", t1.throttle, R9, "FlightModel.trim(ga, 55, 1000)");
  out("actrim(55, 1000, 1).alpha", t1.alpha, R9, "FlightModel.trim(ga, 55, 1000)");
  out("actrim(250, 10500, 2).throttle", t2.throttle, R9, "FlightModel.trim(jet, cruise)");
  out(call("actrim", 2 * aS, 16500, 3) + ".throttle", t3.throttle, R9, "FlightModel.trim(sst, Mach 2 at 16.5 km)");
  out(call("actrim", 55, 1000, 1, 0, F.rad(30)) + ".alpha", t4.alpha, R9, "FlightModel.trim, 30 deg bank");
  out("actrim(120, 1000, 1).throttle", t5.throttle === null ? NaN : t5.throttle, "abs:0", "FlightModel.trim: full power is not enough (null)");
  out("actrim(120, 1000, 1).T", t5.T, R9, "FlightModel.trim thrust required");
  const d1 = F.levelDrag(ga, 55, 1000), d2 = F.levelDrag(jet, 250, 10500);
  out("leveldrag(55, 1000, 1).D", d1.D, R9, "FlightModel.levelDrag(ga, 55, 1000)");
  out("leveldrag(55, 1000, 1).gam", F.glideAngle(ga, 55, 1000), R9, "FlightModel.glideAngle(ga, 55, 1000)");
  out("leveldrag(250, 10500, 2).CL", d2.CL, R9, "FlightModel.levelDrag(jet)");
  out("leveldrag(250, 10500, 2).ok", d2.ok, "abs:0", "FlightModel.levelDrag ok flag");
  for (const [i, key] of [[1, "ga"], [2, "jet"], [3, "sst"]]) {
    const ac = F.AIRCRAFT[key]; let s = F.initialState(ac, { h: 0, V: 0 }); const c = { elev: 0, ail: 0, rud: 0, throttle: 1, trim: 0, brakes: 0 }; let lifted = false;
    for (let k = 0; k < 120 * 120 && !lifted; k++) { c.elev = s.V > ac.Vrot ? 0.5 : 0; s = F.step(s, c, ac, 1 / 120); if (s.event === "liftoff") lifted = true; }
    if (!lifted) throw new Error("no liftoff " + key);
    out(`takeoff(${i}).t`, s.t, R9, `FlightModel.step ground roll to liftoff, ${key} (model.test takeoff)`);
    out(`takeoff(${i}).dist`, s.x, "rel:1e-8", `FlightModel.step ground roll distance, ${key}`);
    out(`takeoff(${i}).V`, s.V, "rel:1e-8", `FlightModel.step liftoff speed, ${key}`);
  }
  err("takeoff(1, 0)", "no thrust, no takeoff");
}
{
  out("moistlapse(300, 1e5)", F.moistLapse(300, 1e5), R12, "FlightModel.moistLapse(300 K, 1000 hPa), about 3.6 K/km");
  out("moistlapse(250, 1e5)", F.moistLapse(250, 1e5), R12, "FlightModel.moistLapse(250 K)");
  out("moistlapse(280, 70000)", F.moistLapse(280, 70000), R12, "FlightModel.moistLapse(280 K, 700 hPa)");
  err("moistlapse(0)", "temperature must be positive");
  out("satvap(20).e", F.satVapour(20), R12, "FlightModel.satVapour(20 C)");
  out("satvap(10, 90000).rs", F.satMixing(10, 90000), R12, "FlightModel.satMixing(10 C, 900 hPa)");
  out("dewpoint(0.01)", F.dewPoint(0.01, 101325), R12, "FlightModel.dewPoint(0.01, 101325)");
  out("dewpoint(satvap(15, 100000).rs, 100000)", F.dewPoint(F.satMixing(15, 1e5), 1e5), "rel:1e-11", "dew point of saturated air is its temperature");
  err("dewpoint(0)", "mixing ratio must be positive");
  const snd = w => { F.setWeather({ clouds: true, ...w }); const s = F.WEATHER.sounding; F.setWeather({ clouds: false, dT: 0, dewSpread: 6, heating: 1.5 }); return s; };
  const s1 = snd({ dT: 8, dewSpread: 6, heating: 2 }), s2 = snd({ dT: 0, dewSpread: 12, heating: 1 }), s3 = snd({ dT: 8, dewSpread: 4, heating: 0 }), s0 = snd({ dT: 0, dewSpread: 6, heating: 1.5 });
  out("sounding(8, 6, 2).lcl", s1.lcl, "abs:0", "FlightModel.sounding: humid warm day (model.test)");
  out("sounding(8, 6, 2).lfc", s1.lfc, "abs:0", "FlightModel.sounding");
  out("sounding(8, 6, 2).el", s1.el, "abs:0", "FlightModel.sounding cloud top");
  out("sounding(8, 6, 2).cape", s1.cape, R9, "FlightModel.sounding CAPE");
  out("sounding(8, 6, 2).cin", s1.cin, "rel:1e-9", "FlightModel.sounding CIN");
  out("sounding(8, 6, 2).conv", s1.convective, "abs:0", "FlightModel.sounding convective");
  out("sounding(0, 12, 1).conv", s2.convective, "abs:0", "FlightModel.sounding: dry day is not convective");
  out("sounding(0, 12, 1).lfc", s2.lfc === null ? NaN : s2.lfc, "abs:0", "FlightModel.sounding: no free convection (null)");
  out("sounding(8, 4, 0).lcl", s3.lcl, "abs:0", "FlightModel.sounding cloud base, about 125 m per K of spread");
  out("sounding().lcl", s0.lcl, "abs:0", "FlightModel.sounding default weather");
  out("sounding().cape", s0.cape, s0.cape === 0 ? "abs:1e-9" : R9, "FlightModel.sounding default weather");
  const z = s1.lcl + 1500;
  out(`updraft(${z})`, F.updraftCore(z, s1), R9, "FlightModel.updraftCore(base + 1500 m), humid warm day");
  out("updraft(9000)", F.updraftCore(9000, s1), R9, "FlightModel.updraftCore(9 km)");
  out("updraft(500, 0, 12, 1)", F.updraftCore(500, s2), "abs:0", "FlightModel.updraftCore: stable air, no updraft");
  const gust = (w, V, h, ac) => { const W = ac.m * F.G0, a0 = F.trimAlpha(ac, V, h), n0 = F.aero(a0, 0, V, h, ac).L / W, n = F.aero(a0 + w / V, 0, V, h, ac).L / W; return { n, dn: n - n0 }; };
  out("gustload(5, 55, 1500, 1).n", gust(5, 55, 1500, F.AIRCRAFT.ga).n, R9, "FlightModel.trimAlpha + aero with alpha + w/V (derivatives() gust kinematics)");
  out("gustload(5, 55, 1500, 1).dalpha", 5 / 55, R12, "model.js: alpha rises by dw cos(gam)/V");
  out("gustload(10, 250, 10500, 2).dn", gust(10, 250, 10500, F.AIRCRAFT.jet).dn, R9, "FlightModel.trimAlpha + aero, jet");
}

/* ------------------------------------------------------------------ airframe: planform.js */
section("airframe: flight-model/planform.js (Planform), defaults from planform.html");
{
  const A = PF.planform(PF.DEFAULT);
  out("planform().S", A.S, R12, "Planform.planform(DEFAULT): 1113.6 m^2");
  out("planform().AR", A.AR, R12, "Planform.planform: AR 3.452");
  out("planform().MAC", A.MAC, R12, "Planform.planform: MAC 26.58 m");
  out("planform().ymac", A.yMAC, R12, "Planform.planform");
  out("planform().xmac", A.xMAC, R12, "Planform.planform");
  out("planform().taper", A.taper, R12, "Planform.planform");
  out("planform().yeq", A.yEq, R12, "Planform.planform: chord equals MAC");
  out("planform().vol", A.vol, R12, "Planform.planform: exact body volume");
  out("planform().b", A.b, R12, "Planform.planform span");
  const T2 = [[0, 0, -12, 1.2, 0, 0.6], [5, -3, -11, 0.8, 0.1, 0.58], [12, -8, -11.5, 0.3, 0.3, 0.55]], B = PF.planform(T2);
  out(`planform(${mat(T2)}).S`, B.S, R12, "Planform.planform on a three-station table");
  out(`planform(${mat(T2)}).MAC`, B.MAC, R12, "Planform.planform on a three-station table");
  out(`planform(${mat(T2.map(r => r.slice(0, 3)))}).AR`, B.AR, R12, "Planform.planform, y/x_LE/x_TE columns only");
  err("planform([1, 0, -1; 2, 0, -1])", "the first station must be the centreline");
  err("planform([0, 0, -1; 2, -2, -1.5; 1, -3, -4])", "y must increase outboard");
  for (const i of [1, 2, 3, 4]) out(`wingsweep(${i}).le`, A.seg[i - 1].le * Math.PI / 180, R12, `Planform.planform seg ${i} leading-edge sweep (deg to rad)`);
  out("wingsweep(2).te", A.seg[1].te * Math.PI / 180, R12, "Planform.planform seg 2 trailing-edge sweep");
  out("wingsweep(4).te", A.seg[3].te * Math.PI / 180, R12, "Planform.planform seg 4 trailing edge: 3.99 deg, set by x_TE at the tip");
  out("wingsweep(3).area", A.seg[2].area, R12, "Planform.planform seg 3 area, both sides");
  err("wingsweep(5)", "four segments");
  const Qd = PF.planformQuad(PF.DEFAULT, 1000), Q50 = PF.planformQuad(PF.DEFAULT, 50);
  out("planquad().S", Qd.S, R12, "Planform.planformQuad(DEFAULT, 1000), page default");
  out("planquad().MAC", Qd.MAC, R12, "Planform.planformQuad");
  out("planquad(50).ymac", Q50.yMAC, R12, "Planform.planformQuad(DEFAULT, 50)");
  out(`planquad(200, ${mat(T2)}).xmac`, PF.planformQuad(T2, 200).xMAC, R12, "Planform.planformQuad on a three-station table");
  err("planquad(0)", "needs at least one interval");
  const P = PF.prepare(PF.DEFAULT);
  out("wingstation(4).xle", PF.station(P, 4)[0], R12, "Planform.station(4)");
  out("wingstation(20).t", PF.station(P, 20)[2], R12, "Planform.station(20)");
  out("wingstation(-9).zc", PF.station(P, -9)[3], R12, "Planform.station(-9): mirrored");
  out("wingstation(40).up", PF.station(P, 40)[4], R12, "Planform.station beyond the tip: clamped");
  err("wingstation(1, [0, 0, -1; 1, 0, -1])", "needs all six columns");
  for (const h of [0, 10000, 15000]) {
    const I = PF.isa(h), H = PF.hydro(h);
    out(`stdatm(${h}).p`, I.p, R12, `Planform.isa(${h})`);
    out(`stdatm(${h}).rho`, I.rho, R12, `Planform.isa(${h})`);
    out(`stdatm(${h}).a`, I.a, R12, `Planform.isa(${h})`);
    // Above 11 km the double version drifts by ~1e-12 as z accumulates over 4000 steps; the 18-digit gear
    // integral lands within 5e-15 of the closed-form ISA pressure, so these cases allow 1e-11.
    const tolH = h > 11000 ? "rel:1e-11" : R12, noteH = h > 11000 ? "; 18-digit decimal vs IEEE double: z accumulated over 4000 steps, gear is closer to Planform.isa" : "";
    out(`stdatm(${h}).phyd`, H.p, tolH, `Planform.hydro(${h}): Simpson integral of the hydrostatic equation${noteH}`);
  }
  out("stdatm(15000).rhohyd", PF.hydro(15000).rho, "rel:1e-11", "Planform.hydro(15000); 18-digit decimal vs IEEE double: z accumulated over 4000 steps, gear is closer to Planform.isa");
  err("stdatm(25000)", "ISA here stops at 20 km");
  const Z = PF.sizing({ S: A.S, mass: 200e3, h: 10e3, mach: 0.8, thrust: 4 * 140e3 }), Z2 = PF.sizing({ S: 1000, mass: 100e3, h: 5000, mach: 0.6, thrust: 400e3 });
  out("airsizing().CL", Z.CL, R12, "Planform.sizing, page default: cruise CL 0.149");
  out("airsizing().WS", Z.WS, R12, "Planform.sizing: W/S 1761 N/m^2");
  out("airsizing().TW", Z.TW, R12, "Planform.sizing: T/W 0.286");
  out("airsizing().q", Z.q, R12, "Planform.sizing: q = rho V^2 / 2");
  out("airsizing().q2", Z.q2, R12, "Planform.sizing: q = gamma p M^2 / 2");
  out("airsizing(100e3, 5000, 0.6, 400e3, 1000).CL", Z2.CL, R12, "Planform.sizing");
  err("airsizing(200e3, 30000)", "ISA here stops at 20 km");
  const E = PF.holdEnds(PF.DEFAULT, 19, 6), Gd = PF.holdGrid(PF.DEFAULT, 19, 6), E2 = PF.holdEnds(PF.DEFAULT, 12, 4);
  out("cargohold().H", E.H, R9, "Planform.holdEnds(DEFAULT, 19, 6), page default");
  out("cargohold().x0", E.x0, "abs:1e-6", "Planform.holdEnds window aft end");
  out("cargohold().floor", E.floor, R9, "Planform.holdEnds");
  out("cargohold().Hscan", Gd.H, R12, "Planform.holdGrid(DEFAULT, 19, 6): 0.05 m scan");
  out("cargohold().sides", Gd.sides, "abs:0", "Planform.holdGrid: extremes at the box sides");
  out("cargohold(12, 4).H", E2.H, R9, "Planform.holdEnds(DEFAULT, 12, 4)");
  err("cargohold(60, 6)", "longer than the centre body");
  const C70 = PF.convergence(PF.DEFAULT, [70])[0], C35 = PF.convergence(PF.DEFAULT, [35])[0];
  out("bodyvol().raw", C70.raw, R9, "Planform.convergence M = 70: as exported (strips inward)");
  out("bodyvol().good", C70.good, R9, "Planform.convergence M = 70: strips reversed");
  out("bodyvol().exact", A.vol, R12, "Planform.planform vol");
  out("bodyvol(35).good", C35.good, R9, "Planform.convergence M = 35");
  const L2 = PF.loft(PF.DEFAULT, 70, 2);
  out("bodyvol(70, 2).good", PF.signedVolume(L2.V, PF.reverseFaces(L2.F, 0, L2.nStrip)), R9, "Planform.loft(DEFAULT, 70, 2) + signedVolume");
  err("bodyvol(2)", "too few chordwise points");
}

/* ------------------------------------------------------------------ game: game.js */
section("game: game/game.js (ParGame)");
{
  out("glidepar()", G.glidePar(), R12, "ParGame.glidePar: 14179.630807 m");
  out("glidepar(2000)", 2000 * Q.polar(0.025, 8, 0.8).LDmax, R12, "h (L/D)max from AeroEq.polar");
  for (const c of [0.9, 0.71, 0.5]) out(`glidegame(${c}).score`, G.glideRun(c).score, R9, `ParGame.glideRun(${c})`);
  out("glidegame(0.9).t", G.glideRun(0.9).t, R9, "ParGame.glideRun(0.9) time");
  out("glidegame(2).CL", 1.4, "abs:0", "ParGame.glideStep clamps CL to 1.4");
  out("glidegame(2).score", G.glideRun(2).score, R9, "ParGame.glideRun(2): held at CLmax");
  const Tr = G.pitchTrim();
  out("pitchtrim().alpha", Tr.alpha, R12, "ParGame.pitchTrim: -0.2076 deg");
  out("pitchtrim().delta", Tr.delta, R12, "ParGame.pitchTrim: 2.0759 deg");
  out("pitchgame(1.5, 1, 0.3).J", G.pitchRun(1.5, 1, 0.3, 0).score, R9, "ParGame.pitchRun(1.5, 1, 0.3, 0)");
  out(`pitchgame(10, 10, 0.8, ${lit(7 / G.DEG)}).J`, G.pitchRun(10, 10, 0.8, 7).score, R9, "ParGame.pitchRun(10, 10, 0.8, 7 deg): Python best grid point");
  out(`pitchgame(10, 8, 0.8, ${lit(6 / G.DEG)}).J`, G.pitchRun(10, 8, 0.8, 6).score, R9, "ParGame.pitchRun(10, 8, 0.8, 6 deg): 2.94008");
  out("pitchgame(0, 0, 0).ise", G.pitchRun(0, 0, 0, 0).ise, R9, "ParGame.pitchRun with no control");
  out("pitchgame(-50, 0, 0).J", G.pitchRun(-50, 0, 0, 0).score, "abs:0", "ParGame.pitchRun diverges: score Infinity");
  out("pitchpar()", G.pitchPar().value, R9, "ParGame.pitchPar: Riccati value 2.922719832");
  const r1 = G.radarRun(3, 3, 5), r2 = G.radarRun(0, 0, 0), r3 = G.radarRun(2, 2, 2);
  out("radargame(200, 34, 12).range", r1.score, R12, "ParGame.radarRun(200 kW, 34 dB, 12 s): 46.40 nmi");
  out("radargame(200, 34, 12).n", r1.n, "abs:0", "ParGame.radarRun pulses per beamwidth");
  out("radargame(25, 28, 2).range", r2.score, R12, "ParGame.radarRun cheapest design");
  out("radargame(100, 32, 6).range", r3.score, R12, "ParGame.radarRun(100 kW, 32 dB, 6 s)");
  out("radargame(400, 36, 2).ok", G.radarRun(4, 4, 0).admissible, "abs:0", "ParGame.radarRun: over budget");
  err("radargame(300, 34, 12)", "300 kW is not on the menu");
  const Par = G.radarPar();
  out("radarpar().range", Par.value, R12, "ParGame.radarPar exhaustive search: 46.40 nmi");
  out("radarpar().kw", Par.best.kw, "abs:0", "ParGame.radarPar best design");
  out("radarpar().gdb", Par.best.gDb, "abs:0", "ParGame.radarPar best design");
  out("radarpar().T", Par.best.T, "abs:0", "ParGame.radarPar best design");
  out("parstars(14100, 14179.630807)", G.stars(14100, 14179.630807, true), "abs:0", "ParGame.stars");
  out("parstars(3.2, 2.922719832, 0)", G.stars(3.2, 2.922719832, false), "abs:0", "ParGame.stars lower is better");
  out("parstars(30, 46.4)", G.stars(30, 46.4, true), "abs:0", "ParGame.stars");
}

process.stdout.write("# Flight library reference cases: generated by Tests/reference/gen_flight.mjs from Where It Isn't's modules\n" + lines.join("\n") + "\n");
