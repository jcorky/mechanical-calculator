// Reference cases for the control library (Source/MechanicalCalculator/Calc/LibControl.cpp).
//
//   node Tests/reference/gen_control.mjs /path/to/where-it-isnt > Tests/reference/control.tsv
//
// Expected values come from Where It Isn't itself: the MPC and approach benches' modules (mpc/mpc.js,
// approach/approach.js), the pure functions of the chapter pages' inline scripts (index.html's analyze()
// and attracting(), ch08's rk4, cycle, ship, lyapExp, period, superstable), closed forms the chapters derive,
// and numbers the pages quote. The matrix cases that the atlas checks only with numpy come from
// gen_control.py, which this script runs and appends.
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const wii = process.argv[2];
if (!wii) {
  console.error("usage: node gen_control.mjs /path/to/where-it-isnt");
  process.exit(2);
}
await import(pathToFileURL(path.join(wii, "mpc/mpc.js")).href);
await import(pathToFileURL(path.join(wii, "approach/approach.js")).href);
const Mpc = globalThis.Mpc, Ap = globalThis.Approach;

// Pure functions lifted verbatim from a page's inline script, between two markers.
function pageCode(file, start, end) {
  const s = fs.readFileSync(path.join(wii, file), "utf8"), i = s.indexOf(start), j = s.indexOf(end, i);
  if (i < 0 || j < 0) throw new Error(`markers not found in ${file}: ${start}`);
  return s.slice(i, j);
}
const ch08 = new Function([
  pageCode("chapters/ch08-dynamics.html", "function rk4(", "\nfunction hexRGB"),
  pageCode("chapters/ch08-dynamics.html", "const VDP =", "\nfunction lyap2"),
  pageCode("chapters/ch08-dynamics.html", "let CYCLE = null;", "\nfunction benchROA"),
  pageCode("chapters/ch08-dynamics.html", "const SK = 0.2", "\nfunction benchShip"),
  pageCode("chapters/ch08-dynamics.html", "function lyapExp", "\nlet BIF"),
  "return { rk4, VDP, cycle, ship, lyapExp, period, superstable, feigenbaum };",
].join("\n"))();
const demo = new Function([
  pageCode("index.html", "function analyze(p)", "\nfunction iterate"),
  pageCode("index.html", "function attracting(an)", "\nfunction updateDemo"),
  "return { analyze, attracting };",
].join("\n"))();

const lines = [];
const num = x => Number.isNaN(x) ? "nan" : x === Infinity ? "inf" : x === -Infinity ? "-inf" : String(x);
const mat = M => "[" + M.map(r => (Array.isArray(r) || ArrayBuffer.isView(r) ? Array.from(r) : [r]).map(num).join(", ")).join("; ") + "]";
const col = v => "[" + Array.from(v).map(num).join("; ") + "]";
const emit = (expr, expected, tol, prov) => lines.push([expr, typeof expected === "string" ? expected : num(expected), tol, prov].join("\t"));
const err = (expr, why) => emit(expr, "error", "abs:0", why);
const section = t => lines.push("# " + t);
const DEG = Math.PI / 180;

/* ====================== riccati: index demo and chapter 3 ====================== */
section("riccati: the scalar Kalman map as a Mobius map (index demo, ch03 bench A, phase1)");
const kfM = (a, q, r, h) => ({ A: a * a * r + q * h * h, B: q * r, C: h * h, D: r });
const mob = (M, z) => (M.A * z + M.B) / (M.C * z + M.D);
const KF = [[0.95, 0.2, 1, 1, 6], [1.4, 0.1, 0.5, 1, 0.2], [0.5, 1.2, 0.3, 1.7, 3], [-1.3, 0.4, 2, 0.6, 10]];   // demo presets and two more
for (const [a, q, r, h, p0] of KF) {
  const M = kfM(a, q, r, h), args = `${a}, ${q}, ${r}, ${h}`;
  emit(`kfmap(${p0}, ${args})`, a * a * r * p0 / (h * h * p0 + r) + q, "rel:1e-13", "phase1 claim 1: f(P) = a^2 r P/(h^2 P + r) + q");
  // Chapter 3 Theorem 1: the cross-ratio contracts by exactly m each step, so P_k in closed form.
  const disc = (M.A - M.D) ** 2 + 4 * M.B * M.C, pp = ((M.A - M.D) + Math.sqrt(disc)) / (2 * M.C), pm = ((M.A - M.D) - Math.sqrt(disc)) / (2 * M.C);
  const tr = M.A + M.D, det = M.A * M.D - M.B * M.C, l1 = (tr + Math.sqrt(tr * tr - 4 * det)) / 2, l2 = (tr - Math.sqrt(tr * tr - 4 * det)) / 2, m = l2 / l1;
  const w0 = (p0 - pp) / (p0 - pm);
  for (const k of [1, 5, 12]) {
    const wk = Math.pow(m, k) * w0;
    emit(`kfiter(${k}, ${p0}, ${args})`, (pp - wk * pm) / (1 - wk), "rel:1e-10", `ch03 Theorem 1: w_k = m^k w_0, m = ${m.toFixed(6)}`);
  }
  emit(`kfw(${p0}, ${args})`, w0, "rel:1e-12", "ch03 bench A: w = (P - p+)/(P - p-)");
  if (Math.abs(Math.pow(m, 7)) > 1e-6) emit(`kfw(kfiter(7, ${p0}, ${args}), ${args}) / kfw(${p0}, ${args})`, Math.pow(m, 7), "rel:1e-7", "ch03 Theorem 1: w_7 / w_0 = m^7");
  const K = pp * h / (h * h * pp + r);
  emit(`kfss(${args}).p`, pp, "rel:1e-12", "ch03 bench A: p+ = ((A - D) + sqrt((A - D)^2 + 4BC))/(2C)");
  emit(`kfss(${args}).pm`, pm, "rel:1e-12", "ch03 bench A: p-");
  emit(`kfss(${args}).K`, K, "rel:1e-12", "ch03 bench A: K = p h/(h^2 p + r)");
  emit(`kfss(${args}).m`, m, "rel:1e-10", "ch03: multiplier = lambda_min/lambda_max of M");
  emit(`kfss(${args}).m - kfss(${args}).pole^2`, 0, "abs:1e-12", "phase1 claim 3: f'(P*) = a^2 (1 - K h)^2");
  const an = demo.analyze(M), att = demo.attracting(an);
  emit(`kfss(${args}).p`, att, "abs:2e-9", "index demo analyze()/attracting(): fixed point (page rounds to 1e-9)");
}
emit("kfiter(25, 8, 1, 0, 1.3, 0.7)", 8 * 1.3 / (1.3 + 25 * 0.7 * 0.7 * 8), "rel:1e-12", "phase1 claim 4: parabolic P_k = P0 r/(r + k h^2 P0)");
emit("kfiter(49, 8, 1, 0, 1.3, 0.7)", 8 * 1.3 / (1.3 + 49 * 0.7 * 0.7 * 8), "rel:1e-12", "phase1 claim 4: 1/k decay");
emit("kfiter(0, 6)", 6, "abs:0", "zero steps");
emit("kfss(0.5, 0.3, 1, 0).p", 0.3 / (1 - 0.25), "rel:1e-14", "h = 0, |a| < 1: P = q/(1 - a^2)");
err("kfss(1.08, 0.3, 1, 0)", "demo preset 'Unobserved, diverges': no finite steady state");
err("kfw(3, 1, 0, 1, 1)", "demo preset 'Constant (1/k)': parabolic, no cross-ratio");
err("kfiter(-1, 6)", "k must be a whole number");
err("kfiter(2.5, 6)", "k must be a whole number");
err("kfiter(2e6, 6)", "k <= 1e6 keeps a key press under 50 ms");

section("riccati: the LQR reading (index demo: q state weight, rho control weight, b input gain)");
for (const [a, q, rho, b] of [[0.95, 0.2, 1, 1], [1.4, 0.1, 0.5, 1], [-0.7, 2, 0.1, 1.5]]) {
  // phase1 claim 5: the regulator map is the filter map with h -> b, r -> rho; its fixed point by iteration.
  let S = 0; for (let i = 0; i < 20000; i++) S = q + a * a * S - (a * b * S) ** 2 / (rho + b * b * S);
  const K = a * b * S / (rho + b * b * S);
  emit(`lqrss(${a}, ${q}, ${rho}, ${b}).S`, S, "rel:1e-12", "phase1 claim 5: S = q + a^2 S - (abS)^2/(rho + b^2 S), iterated");
  emit(`lqrss(${a}, ${q}, ${rho}, ${b}).K`, K, "rel:1e-12", "index demo: K = a b S/(rho + b^2 S)");
  emit(`lqrss(${a}, ${q}, ${rho}, ${b}).pole`, a - b * K, "rel:1e-11", "closed-loop pole a - bK");
  emit(`lqrss(${a}, ${q}, ${rho}, ${b}).m`, (a - b * K) ** 2, "rel:1e-10", "index demo: multiplier (a - bK)^2");
  const A = [[a]], B = [[b]];
  emit(`dlqr(${a}, ${b}, ${q}, ${rho})`, Ap.solve([[rho + b * b * Ap.dare(A, B, [[q]], [[rho]]).P[0][0]]], [[b * Ap.dare(A, B, [[q]], [[rho]]).P[0][0] * a]])[0][0], "rel:1e-9", "approach.js dare(): the 1 x 1 matrix LQR agrees");
}
err("lqrss(1.2, 0.2, 1, 0)", "b = 0 and |a| >= 1: unreachable unstable state");

section("riccati: free 2x2 matrix (index demo presets, analyze())");
const GEN = [[0.5, -1, 1, 0.5, 2], [2, 1, 1, 1, -2], [1, 1, 0, 1, -3], [3, 1, 2, 2, 1], [0, 1, 1, 0, 2], [1, 2, 2, 4, 1]];
for (const [al, be, ga, de, z0] of GEN) {
  const an = demo.analyze({ A: al, B: be, C: ga, D: de }), att = demo.attracting(an), args = `${al}, ${be}, ${ga}, ${de}`;
  const type = { degenerate: 0, identity: 1, parabolic: 2, hyperbolic: 3, elliptic: 4 }[an.type];
  emit(`mobclass(${args}).type`, type, "abs:0", `index demo analyze(): ${an.type}`);
  emit(`mobclass(${args}).tr`, an.tr, "rel:1e-14", "index demo analyze(): tr M");
  emit(`mobclass(${args}).disc`, an.disc, "rel:1e-14", "index demo analyze(): tr^2 - 4 det");
  if (an.fps.length) emit(`mobclass(${args}).z1`, an.fps[0], an.fps[0] === Infinity ? "abs:0" : "abs:2e-9", "index demo analyze(): first fixed point");
  if (an.fps.length > 1) emit(`mobclass(${args}).z2`, an.fps[1], "abs:2e-9", "index demo analyze(): second fixed point");
  if (att !== null && an.type === "hyperbolic") {
    emit(`mobclass(${args}).att`, att, "abs:2e-9", "index demo attracting()");
    emit(`mobclass(${args}).m`, an.mult(att), "rel:1e-8", "index demo: multiplier f'(z*)");
  }
  if (an.type === "elliptic") {
    const th = Math.atan2(Math.sqrt(-an.disc), an.tr);
    emit(`mobclass(${args}).rot`, 2 * th, "rel:1e-12", "elliptic: lambda1/lambda2 = e^{i rot}, |m| = 1");
  }
  // Iterates: f^k = f_{M^k} (index demo Theorem: f_M o f_N = f_MN), computed by matrix powers.
  for (const k of [3, 10]) {
    let P = [[1, 0], [0, 1]]; for (let i = 0; i < k; i++) P = [[P[0][0] * al + P[0][1] * ga, P[0][0] * be + P[0][1] * de], [P[1][0] * al + P[1][1] * ga, P[1][0] * be + P[1][1] * de]];
    const zk = (P[0][0] * z0 + P[0][1]) / (P[1][0] * z0 + P[1][1]);
    if (Number.isFinite(zk) && Math.abs(zk) < 1e6) emit(`mobiter(${k}, ${z0}, ${args})`, zk, "rel:1e-9", "ch03 definition: f_M o f_N = f_MN, M^k applied once");
  }
  emit(`mobius(${z0}, ${args})`, (al * z0 + be) / (ga * z0 + de), "rel:1e-15", "definition");
}
emit("mobius(inf)", 0.5, "abs:0", "f(inf) = al/ga");
emit("mobius(-0.5)", "inf", "abs:0", "f(-de/ga) = inf");
emit("mobius(inf, 1, 1, 0, 1)", "inf", "abs:0", "ga = 0 fixes inf");
emit("mobiter(1, -0.5)", "inf", "abs:0", "the iterate passes through inf");
emit("mobiter(2, -0.5)", 0.5, "abs:0", "and comes back: f(inf) = al/ga");
err("mobiter(-2, 1)", "k must be a whole number");

section("riccati: the Riccati ODE and finite escape (ch03 counterexample 1, phase4 claim 4)");
const escForm = (al, be, ga, p0) => { const om = Math.sqrt(al * ga - be * be); return (Math.PI / 2 - Math.atan((ga * p0 + be) / om)) / om; };
for (const [al, be, ga, p0] of [[1, 0.2, 1, 0], [1, 0.2, 1, 2], [2, -0.5, 0.7, 0], [1, 0.2, 2, -1]]) {
  emit(`riccesc(${al}, ${be}, ${ga}, ${p0})`, escForm(al, be, ga, p0), "rel:1e-13", "phase4 claim 4: t = (pi/2 - atan((ga p0 + be)/om))/om");
}
{ // H-infinity reading of bench B: alpha = q, beta = a, gamma = -(c^2/r - 1/g^2); escape when 1/g^2 is large.
  const a = 0.2, q = 1, c = 1, g = 2.5; emit(`riccesc(${q}, ${a}, ${-(c - g)}, 0)`, escForm(q, a, -(c - g), 0), "rel:1e-13", "ch03 bench B: H-infinity term 1/g^2 = 2.5, elliptic flow");
  emit(`riccesc(${q}, ${a}, ${-(c - 0)}, 0)`, "inf", "abs:0", "ch03 Theorem 4: Kalman signs never escape");
}
{ // Hyperbolic flows: escape only above the upper equilibrium; checked by RK4 as on the page.
  const rk = (al, be, ga, p0) => { const F = p => al + 2 * be * p + ga * p * p, dt = 1e-5; let p = p0, t = 0; while (t < 50) { const k1 = F(p), k2 = F(p + dt / 2 * k1), k3 = F(p + dt / 2 * k2), k4 = F(p + dt * k3); p += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt; if (Math.abs(p) > 1e6) return t; } return Infinity; };
  emit("riccesc(-1, 0, 1, 2)", rk(-1, 0, 1, 2), "abs:1e-4", "RK4 of p' = p^2 - 1 from 2 (closed form ln(3)/2)");
  emit("riccesc(-1, 0, 1, 2)", Math.log(3) / 2, "rel:1e-13", "closed form: u' = u^2 - k^2 escapes at ln((u0 + k)/(u0 - k))/(2k)");
  emit("riccesc(0, 0, 1, 2)", rk(0, 0, 1, 2), "abs:1e-4", "RK4 of p' = p^2 from 2 (parabolic: 1/p0)");
  emit("riccesc(-1, 0, 1, 0.5)", "inf", "abs:0", "below the upper equilibrium: converges, no escape");
  emit("riccesc(1, 0.5, 0, 3)", "inf", "abs:0", "linear (ga = 0): no escape");
}

section("riccati: proportional navigation (phase1 claim 6, ch08 guise)");
for (const [y0, v0, tgo, N] of [[30, -2, 4, 3], [100, 5, 10, 3], [30, -2, 4, 4], [-50, 3, 6, 5]]) {
  const zem = y0 + v0 * tgo;
  // Fly the linearized PN loop y'' = -N ZEM/tgo^2 by RK4 and integrate the effort a^2.
  const n = 200000, dt = tgo * (1 - 1e-6) / n, f = (t, s) => { const tg = tgo - t, a = -N * (s[0] + s[1] * tg) / (tg * tg); return [s[1], a, a * a]; };
  let s = [y0, v0, 0], t = 0;
  for (let i = 0; i < n; i++) { const k1 = f(t, s), k2 = f(t + dt / 2, s.map((v, j) => v + dt / 2 * k1[j])), k3 = f(t + dt / 2, s.map((v, j) => v + dt / 2 * k2[j])), k4 = f(t + dt, s.map((v, j) => v + dt * k3[j])); s = s.map((v, j) => v + dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])); t += dt; }
  emit(`pn(${y0}, ${v0}, ${tgo}, ${N}).zem`, zem, "rel:1e-15", "ZEM = y + v tgo");
  emit(`pn(${y0}, ${v0}, ${tgo}, ${N}).a0`, -N * zem / tgo ** 2, "rel:1e-14", "phase1 claim 6: a = -N ZEM/tgo^2");
  emit(`pn(${y0}, ${v0}, ${tgo}, ${N}).effort`, s[2], "rel:1e-6", "RK4 of the PN loop, effort integral of a^2");
}
emit("pn().effort", 3 * 22 ** 2 / 4 ** 3, "rel:1e-14", "phase1 claim 6: minimum effort 3 ZEM^2/tgo^3 at N = 3");
err("pn(30, -2, 0)", "tgo must be > 0");

/* ====================== control: approach bench ====================== */
section("control: approach bench (approach.js design(), closedCov(), margins(), doyle())");
const D = Ap.design(Ap.DEF), ct = Ap.continuous(Ap.DEF), dz = Ap.discretize(Ap.DEF), wt = Ap.lqrWeights(Ap.DEF);
const setC = `Ac = ${mat(ct.A)}; Bc = ${mat(ct.B)}; `;
for (const [i, j] of [[1, 1], [1, 2], [2, 3], [4, 3], [4, 4], [5, 5]]) emit(setC + `c2d(Ac, Bc, 0.1)(${i}, ${j})`, dz.Ad[i - 1][j - 1], "abs:1e-13", "approach.js discretize(): Ad by expm (Taylor 20, scaling and squaring)");
for (const i of [1, 2, 3, 4]) emit(setC + `c2d(Ac, Bc, 0.1)(${i}, 6)`, dz.Bd[i - 1][0], "abs:1e-13", "approach.js discretize(): Bd");
emit(setC + "size(c2d(Ac, Bc, 0.1)).cols", 6, "abs:0", "[Ad, Bd] side by side");
const setW = `Ac = ${mat(ct.A)}; W = ${mat(ct.W)}; `;
for (const [i, j] of [[1, 1], [1, 5], [2, 2], [4, 4], [5, 5], [3, 4]]) emit(setW + `c2dq(Ac, W, 0.1)(${i}, ${j})`, dz.Qd[i - 1][j - 1], `abs:${Math.max(1e-16, 1e-11 * Math.abs(dz.Qd[i - 1][j - 1]))}`, "approach.js discretize(): Van Loan Qd");
{ // Scalar check of Van Loan: x' = -a x + w, Qd = W (1 - e^{-2a dt})/(2a).
  emit("c2dq(-0.5, 2, 0.3)", 2 * (1 - Math.exp(-2 * 0.5 * 0.3)) / (2 * 0.5), "rel:1e-13", "closed form for a scalar first-order process");
  emit("c2d(-0.5, 1, 0.3)(2)", (1 - Math.exp(-0.15)) / 0.5, "rel:1e-13", "closed form Bd = (1 - e^{-a dt})/a");
}
err("c2d([1, 2; 3, 4], [1; 2; 3], 0.1)", "B must have n rows");
err("c2d(1, 1, 0)", "dt must be > 0");
err("c2dq([1, 2; 3, 4], 1, 0.1)", "W must match A");

const setD = `Ad = ${mat(D.Ad)}; Bd = ${mat(D.Bd)}; Q = ${mat(D.Q)}; Ru = ${num(D.Ru[0][0])}; `;
for (let k = 0; k < 5; k++) emit(setD + `dlqr(Ad, Bd, Q, Ru)(${k + 1})`, D.K[0][k], "rel:1e-8", "approach.js design(): LQR gain by Riccati iteration (test: 0.030636, 0.044747, 0.823904, 0.128041, 0.03872)");
for (const [i, j] of [[1, 1], [2, 5], [3, 3]]) emit(setD + `dare(Ad, Bd, Q, Ru)(${i}, ${j})`, D.P[i - 1][j - 1], "rel:1e-8", "approach.js dare(): P by iteration");
emit(setD + "rho(Ad - Bd*dlqr(Ad, Bd, Q, Ru))", Math.exp(-0.1 / 4), "rel:1e-10", "approach test: the gust mode e^{-dt/tau} is the spectral radius");
{
  const x0 = [10, 0, 0, 0, 0], c = Ap.lqrCost(D, x0);
  emit(setD + `clcost(Ad, Bd, Q, Ru, dlqr(Ad, Bd, Q, Ru), ${col(x0)})`, c.J, "rel:1e-8", "approach.js lqrCost(): 4000 simulated steps");
  emit(setD + `x0 = ${col(x0)}; x0'*dare(Ad, Bd, Q, Ru)*x0`, c.Jp, "rel:1e-8", "approach.js lqrCost(): x0'Px0");
}
const setF = `Ad = ${mat(D.Ad)}; C = ${mat(D.C)}; Qd = ${mat(D.Qd)}; R = ${mat(D.R)}; `;
for (const [i, j] of [[1, 1], [2, 2], [5, 1], [4, 2]]) emit(setF + `kalman(Ad, C, Qd, R)(${i}, ${j})`, D.M[i - 1][j - 1], "rel:1e-8", "approach.js design(): M = Pf C'(C Pf C' + R)^-1, Pf by doubling");
for (const [i, j] of [[1, 1], [5, 5]]) emit(setF + `kfcov(Ad, C, Qd, R)(${i}, ${j})`, D.Pf[i - 1][j - 1], "rel:1e-8", "approach.js dareDoubling(): prior covariance");
emit(setF + "P = kfcov(Ad, C, Qd, R); M = kalman(Ad, C, Qd, R); sqrt(((eye(5) - M*C)*P)(1, 1))", Math.sqrt(D.Pp[0][0]), "rel:1e-8", "approach test: posterior sigma of d 0.2495");
emit(setF + "P = kfcov(Ad, C, Qd, R); M = kalman(Ad, C, Qd, R); sqrt(((eye(5) - M*C)*P)(5, 5))", Math.sqrt(D.Pp[4][4]), "rel:1e-8", "approach test: posterior sigma of w 0.8795");
const setL = `Ad = ${mat(D.Ad)}; Bd = ${mat(D.Bd)}; C = ${mat(D.C)}; K = ${mat(D.K)}; M = ${mat(D.M)}; `;
emit(setL + "rho(lqgcl(Ad, Bd, C, K, M))", Ap.specRad(Ap.lqgClosed(D)), "rel:1e-10", "approach.js lqgClosed(): the LQG loop");
emit(setL + "rho(lqgcl(Ad, Bd, C, K, M, 2))", Ap.specRad(Ap.lqgClosed(D, 2)), "rel:1e-10", "approach.js lqgClosed(D, 2): gain 2 at the input");
emit(setL + "rho(lqgcl(Ad, Bd, C, K, M, 3))", Ap.specRad(Ap.lqgClosed(D, 3)), "rel:1e-10", "approach.js lqgClosed(D, 3): beyond the LQG margin 2.86");
emit(setL + "rho(lqgcl(Ad, Bd, C, K, M)) - max(rho(Ad - Bd*K), rho(Ad - Ad*M*C))", 0, "abs:1e-9", "approach section 3: separation, eig = eig(A - BK) and eig(A - LC)");
err("lqgcl([1, 0; 0, 1], [1; 0], [1, 0], [1, 1], [1, 0])", "M must be n x p");
const mg = Ap.margins(D);
emit(setL + "gm(Ad, Bd, K).hi", mg.lqr[1], "rel:1e-7", "approach.js margins(): LQR upper gain margin (test 8.17437)");
emit(setL + "gm(Ad, Bd, K).lo", mg.lqr[0], "rel:1e-12", "approach.js margins(): LQR lower edge (search limit 1e-3)");
emit(setD + "lqrgm(Ad, Bd, Q, Ru).rho", mg.rho, "rel:1e-8", "approach.js margins(): rho = sqrt(R/(R + B'PB))");
emit(setD + "lqrgm(Ad, Bd, Q, Ru).lo", mg.guar[0], "rel:1e-8", "approach.js margins(): guarantee 1/(1 + rho) (test 0.532576)");
emit(setD + "lqrgm(Ad, Bd, Q, Ru).hi", mg.guar[1], "rel:1e-7", "approach.js margins(): guarantee 1/(1 - rho) (test 8.17432)");
emit(setL + "retdiff(Ad, Bd, K)", Ap.returnDifferenceMin(D), "rel:1e-9", "approach.js returnDifferenceMin(): min |1 + L| on the unit circle");
{ // continuous LQR (approach continuousLqr: Newton-Kleinman from the discrete gain)
  const cq = Ap.continuousLqr(D), setQ = `Ac = ${mat(ct.A)}; Bc = ${mat(ct.B)}; Q = ${mat(wt.Q)}; Ru = ${num(wt.Ru[0][0])}; `;
  for (let k = 0; k < 5; k++) emit(setQ + `lqr(Ac, Bc, Q, Ru)(${k + 1})`, cq.K[0][k], "rel:1e-7", "approach.js care(): continuous LQR by Newton-Kleinman");
  emit(setQ + "care(Ac, Bc, Q, Ru)(1, 1)", cq.P[0][0], "rel:1e-7", "approach.js care(): P");
}
for (const g of ["lqr", "lqg", "pid"]) {
  const mode = { lqr: 1, lqg: 2, pid: 3 }[g], c = Ap.closedCov(D, g);
  emit(`glide(${mode}).sd`, c.sd, "rel:1e-8", `approach.js closedCov("${g}"): sigma of d at the decision height (test: 0.5615 / 0.6751 / 0.9535)`);
  emit(`glide(${mode}).J`, c.J, "rel:1e-8", `approach.js closedCov("${g}"): cost E[x'Qx + u'Ru] (test 0.969 / 1.196 / 1.706)`);
  emit(`glide(${mode}).su`, c.su * DEG, "rel:1e-8", `approach.js closedCov("${g}"): command sigma, converted from degrees`);
  emit(`glide(${mode}).rho`, c.rho, "rel:1e-10", `approach.js closedCov("${g}"): spectral radius`);
  emit(`glide(${mode}).hi`, mg[g][1], "rel:1e-7", `approach.js margins(): upper gain edge of the ${g} loop (LQR 8.174, LQG 2.860)`);
  emit(`glide(${mode}).lo`, mg[g][0], "rel:1e-7", `approach.js margins(): lower gain edge of the ${g} loop`);
}
emit("glide().rd", Ap.returnDifferenceMin(D), "rel:1e-9", "approach.js returnDifferenceMin()");
emit("glide().glo", mg.guar[0], "rel:1e-8", "approach.js margins(): guaranteed interval");
emit("glide().ghi", mg.guar[1], "rel:1e-7", "approach.js margins(): guaranteed interval");
{ // other noise settings: the filter changes, the regulator does not
  const p = { ...Ap.DEF, sg: 0.5, sd: 1, tau: 6 }, D2 = Ap.design(p), c = Ap.closedCov(D2, "lqg"), m2 = Ap.margins(D2);
  emit("glide(2, 0.5, 6, 1, 0.3).sd", c.sd, "rel:1e-8", "approach.js closedCov(): sg 0.5, tau 6, sd 1");
  emit("glide(2, 0.5, 6, 1, 0.3).hi", m2.lqg[1], "rel:1e-7", "approach.js margins(): LQG upper edge, sg 0.5, tau 6, sd 1");
  const c3 = Ap.closedCov(Ap.design({ ...Ap.DEF, sg: 3 }), "pid");
  emit("glide(3, 3).J", c3.J, "rel:1e-8", "approach.js closedCov(\"pid\"): gust sigma 3 m/s");
}
err("glide(4)", "mode must be 1, 2 or 3");
err("glide(2, 0)", "sg must be > 0");
for (const q of [1, 100, 1e4]) {
  const d = Ap.doyle(q);
  emit(`doyle(${q}).lo`, d.lqg[0], "rel:1e-7", `approach.js doyle(${q}): LQG lower gain edge (test q = 1: 0.9208, q = 1e4: 0.9905)`);
  emit(`doyle(${q}).hi`, d.lqg[1], "rel:1e-7", `approach.js doyle(${q}): LQG upper gain edge (test q = 1: 1.0557, q = 1e4: 1.0001)`);
  emit(`doyle(${q}).rlo`, d.lqr[0], "rel:1e-7", `approach.js doyle(${q}): LQR lower edge, below 1/2`);
  emit(`doyle(${q}).f`, d.f, "rel:1e-15", "Doyle closed form f = 2 + sqrt(4 + q)");
  emit(`lqr([1, 1; 0, 1], [0; 1], ${q}*[1, 1; 1, 1], 1)(2)`, d.f, "rel:1e-8", "approach.js doyle(): Newton-Kleinman agrees with the closed form K = f[1 1]");
}
emit("doyle(100, 1).d", Ap.doyle(100, 1).d, "rel:1e-15", "Doyle closed form d = 2 + sqrt(4 + sigma)");
emit("doyle(100, 1).hi", Ap.doyle(100, 1).lqg[1], "rel:1e-7", "approach.js doyle(100, 1): q != sigma");
err("doyle(-1)", "q must be >= 0");
for (const [a, b, q, r] of [[1, 1, 1, 1], [-2, 0.5, 3, 0.2], [0.3, 2, 0.1, 5]]) {   // scalar CARE 2aP - b^2P^2/r + q = 0, stabilizing root
  const P = (a * r + Math.sqrt(a * a * r * r + q * r * b * b)) / (b * b);
  emit(`care(${a}, ${b}, ${q}, ${r})`, P, "rel:1e-10", "closed form: P = (a r + sqrt(a^2 r^2 + q r b^2))/b^2");
  emit(`lqr(${a}, ${b}, ${q}, ${r})`, b * P / r, "rel:1e-10", "closed form: K = b P/r, a - bK = -sqrt(a^2 + q b^2/r)");
}

section("control: generic solver checks on the MPC plant (mpc.js dare(), problem())");
const PL = Mpc.plant(0.2), setP = "A = [1, 0.2; 0, 1]; B = [0.2^2/2; 0.2]; Q = [1, 0; 0, 0]; ";
for (const r of [0.1, 0.01, 1]) {
  const d = Mpc.dare(PL, r);
  emit(setP + `dlqr(A, B, Q, ${r})(1)`, d.K[0], "rel:1e-10", `mpc.js dare(): LQR gain, r = ${r}`);
  emit(setP + `dlqr(A, B, Q, ${r})(2)`, d.K[1], "rel:1e-10", `mpc.js dare(): LQR gain, r = ${r}`);
  emit(setP + `dare(A, B, Q, ${r})(1, 2)`, d.P[0][1], "rel:1e-10", `mpc.js dare(): P, r = ${r}`);
  emit(setP + `gm(A, B, dlqr(A, B, Q, ${r})).hi`, Ap.interval(g => Ap.specRad(Ap.sub(PL.A, Ap.scale(Ap.mul(PL.B.map(b => [b]), [d.K]), g))))[1], "rel:1e-7", `approach.js interval(): LQR gain margin of the double integrator, r = ${r}`);
  for (const N of [1, 3, 10]) {
    for (const [tag, Pf, pf] of [["P", d.P, `dare(A, B, Q, ${r})`], ["0", [[0, 0], [0, 0]], "zeros(2)"]]) {
      const pr = Mpc.problem(PL, N, r, Pf, -Infinity, Infinity), kn = [0, 1].map(j => pr.Minv[0].reduce((s, v, k) => s + v * pr.F[k][j], 0));
      emit(setP + `mpcgain(A, B, Q, ${r}, ${pf}, ${N})(1)`, kn[0], "rel:1e-9", `mpc.js problem(): first row of M^-1 F, N = ${N}, terminal ${tag}`);
      emit(setP + `mpcgain(A, B, Q, ${r}, ${pf}, ${N})(2)`, kn[1], "rel:1e-9", `mpc.js problem(): N = ${N}, terminal ${tag}${tag === "P" ? " (= LQR for every N)" : ""}`);
      if (tag === "0" && N >= 3) {
        const x0 = [3, -1], cl = Mpc.closedLoop(pr, x0, d.P, 150, { tail: true });
        let J = cl.J, how = "mpc.js closedLoop() with its exact linear tail";
        if (cl.tail === undefined) {   // the loop is too slow for the 150-step window: sum the MPC law's own rollout to convergence
          let x = x0.slice(); J = 0;
          for (let t = 0; t < 2e6 && Math.max(Math.abs(x[0]), Math.abs(x[1])) > 1e-15; t++) { const u = -(kn[0] * x[0] + kn[1] * x[1]); J += Mpc.stage(PL, x, u, r); x = Mpc.step(PL, x, u); }
          how = "mpc.js step()/stage() rollout of the problem()'s MPC law to convergence";
        }
        emit(setP + `clcost(A, B, Q, ${r}, mpcgain(A, B, Q, ${r}, zeros(2), ${N}), [3; -1])`, J, "rel:1e-8", `${how}, unconstrained MPC, no terminal cost, N = ${N}`);
      }
    }
  }
}
err(setP + "mpcgain(A, B, Q, 0.1, zeros(2), 0)", "N must be >= 1");
err(setP + "mpcgain(A, B, Q, 0.1, zeros(3), 5)", "Pf must be n x n");
err(setP + "clcost(A, B, Q, 0.1, [0, 0], [1; 0])", "K = 0 leaves the double integrator unstable");
err(setP + "gm(A, B, [0, 0])", "unstable at g = 1");
err(setP + "lqrgm(A, [B, B], Q, eye(2))", "one input only");
err(setP + "retdiff(A, [B, B], [1, 1; 1, 1])", "SISO only");
err("dare([1, 2], [1], 1, 1)", "A must be square");
err("dlqr([1, 0; 0, 1], [1; 0; 0], eye(2), 1)", "B must have n rows");
err("dlqr([2, 0; 0, 1], [0; 1], eye(2), 1)", "unstabilizable: mode 2 cannot be reached");
err("care([1, 0; 0, 1], [0; 0], eye(2), 1)", "unstabilizable: no stabilizing solution");
err("lqr(eye(2), [1; 0], eye(2), [1, 0; 0, 1])", "R must be m x m");
err("kalman([1, 0; 0, 1], [1, 0, 0], eye(2), 1)", "C must have n columns");
err("kalman([2, 0; 0, 1], [0, 1], eye(2), 1)", "undetectable unstable mode");

section("control: constrained MPC (mpc.js gi(), closedLoop(), satLoop())");
for (const [um, y0, r] of [[1, 3.5, 0.01], [0.5, 3.5, 0.01], [2, 5, 0.1], [1, 3.5, 1]]) {
  const d = Mpc.dare(PL, r), sat = Mpc.satLoop(PL, d.K, d.P, r, [y0, 0], -um, um, 150);
  emit(`lanesat(${um}, ${y0}, ${r}).J`, sat.J, "rel:1e-9", "mpc.js satLoop(): saturated LQR closed-loop cost");
  emit(`lanesat(${um}, ${y0}, ${r}).os`, Math.max(0, -Math.min(...sat.X.map(x => x[0]))), "abs:1e-9", "mpc section 2: saturated LQR overshoot past the lane centre");
  for (const N of [1, 5, 20]) {
    for (const term of [1, 0]) {
      const pr = Mpc.problem(PL, N, r, term ? d.P : [[0, 0], [0, 0]], -um, um), cl = Mpc.closedLoop(pr, [y0, 0], d.P, 150, { tail: true });
      if (term === 0 && N === 1) continue;
      const q = Mpc.at(pr, [y0, 0]), s = Mpc.gi(pr, q);
      emit(`lanempc(${N}, ${um}, ${y0}, ${r}, ${term}).J`, cl.J, "rel:1e-8", `mpc.js closedLoop(): MPC closed-loop cost, N = ${N}, ${term ? "terminal cost P" : "no terminal cost"}`);
      emit(`lanempc(${N}, ${um}, ${y0}, ${r}, ${term}).os`, Math.max(0, -Math.min(...cl.X.map(x => x[0]))), "abs:1e-8", "mpc section 2: MPC overshoot");
      emit(`lanempc(${N}, ${um}, ${y0}, ${r}, ${term}).V`, Mpc.cost(pr, q, s.U), "rel:1e-9", "mpc.js gi() + cost(): open-loop value V_N(x0)");
    }
  }
  emit(`lanempc(1, ${um}, ${y0}, ${r}).J - lanesat(${um}, ${y0}, ${r}).J`, 0, `abs:${1e-8 * (1 + sat.J)}`, "mpc section 4: at N = 1 MPC is exactly saturated LQR");
}
err("lanempc(0)", "N must be >= 1");
err("lanempc(20, 1, 3.5, 0.01, 2)", "term must be 0 or 1");
err("lanesat(0)", "umax must be > 0");
const ZERO = [[0, 0], [0, 0]], MODES = ["none", "P", "stop", "origin"];
for (const [N, h0, v0, uh] of [[22, 40, 8, 2], [30, 40, 8, 2], [10, 20, 5, 2.5], [40, 40, 8, 2], [22, 60, 12, 3]]) {
  const r = 0.1, d = Mpc.dare(PL, r);
  for (let mode = 0; mode < 4; mode++) {
    const m = MODES[mode], pr = Mpc.problem(PL, N, r, m === "none" ? ZERO : d.P, -1, uh, { pmin: 0, terminal: m === "stop" || m === "origin" ? m : "none" });
    const cl = Mpc.closedLoop(pr, [h0, -v0], d.P, 150, { keep: true }), args = `${N}, ${h0}, ${v0}, ${uh}, ${mode}`;
    emit(`lander(${args}).fail`, cl.fail === null ? -1 : cl.fail, "abs:0", `mpc.js closedLoop(): lander, ${m}${N === 22 && h0 === 40 && mode === 0 ? " (page: crash unavoidable at step 10)" : ""}`);
    if (cl.fail !== null) {
      emit(`lander(${args}).alt`, cl.X[cl.fail][0], "rel:1e-9", "mpc section 3: altitude where the QP became infeasible");
      emit(`lander(${args}).vel`, cl.X[cl.fail][1], "rel:1e-9", "mpc section 3: velocity there");
      emit(`lander(${args}).alt - lander(${args}).vel^2/(2*${uh})`, cl.X[cl.fail][0] - cl.X[cl.fail][1] ** 2 / (2 * uh), "rel:1e-9", "mpc section 3: full braking from there still goes below the ground");
    } else {
      emit(`lander(${args}).J`, cl.J, "rel:1e-8", "mpc.js closedLoop(): lander closed-loop cost");
    }
  }
}
err("lander(22, 40, 8, 2, 4)", "mode must be 0..3");
err("lander(22, 40, 8, 0)", "braking limit must be > 0");

/* ====================== dynamics: chapter 8 ====================== */
section("dynamics: region of attraction (ch08 bench, page cycle())");
const cyc = ch08.cycle();
emit("roa().basin", cyc.area, "rel:1e-9", "ch08 page cycle(): area inside the unstable limit cycle (page: 13.72)");
emit("roa().c", 2.305, "abs:0.002", "ch08 page and phase9: V < 2.305 with Q = I (grid of 1440 rays)");
emit("roa().area", 6.48, "abs:0.01", "ch08 page: certified area 6.48 of 13.72");
emit("roa().frac", 0.47, "abs:0.005", "ch08 counterexample: the quadratic certificate covers 47% of the basin");
err("roa(0)", "q must be > 0");

section("dynamics: Nomoto ship heading (ch08 page ship())");
const CTRL = ["pid", "pid_aw", "smc", "smc_sign"];
for (const [stepDeg, waveDeg] of [[60, 3], [90, 3], [30, 5], [120, 0]]) {
  for (let c = 0; c < 4; c++) {
    const R = ch08.ship(CTRL[c], stepDeg, waveDeg), o = R.out, deg = v => v * 180 / Math.PI;
    const over = Math.max(...o.map(p => deg(p[1]))) - stepDeg; let settle = 0; o.forEach(p => { if (Math.abs(deg(p[1]) - stepDeg) > 2) settle = p[0]; });
    const ri = o.findIndex(p => Math.abs(p[4]) < 0.004), reach = ri >= 0 ? o[ri][0] : NaN, bound = R.eta > 0 ? Math.abs(o[0][4]) / R.eta : Infinity;
    const args = `${c + 1}, ${stepDeg}*pi/180, ${waveDeg}*pi/180`, base = stepDeg === 60 && waveDeg === 3;
    emit(`nomoto(${args}).over`, Math.max(0, over) * DEG, "abs:1e-9", `ch08 ship(): overshoot, ${CTRL[c]}${base ? " (page: PID 12.3 deg, anti-windup 8.0 deg)" : ""}`);
    emit(`nomoto(${args}).settle`, settle, "abs:1e-9", `ch08 ship(): last time outside 2 deg, ${CTRL[c]}`);
    emit(`nomoto(${args}).sw`, R.sw, c === 3 ? "rel:0.01" : "abs:0", `ch08 ship(): rudder reversals, ${CTRL[c]}${base ? " (page: 8477 sign law vs 28 boundary layer)" : ""}`);
    if (c >= 2) {
      emit(`nomoto(${args}).reach`, reach, "abs:1e-9", `ch08 ship(): first time |s| < phi${base ? " (page: 16.5 s)" : ""}`);
      emit(`nomoto(${args}).bound`, bound, "rel:1e-12", `ch08 Theorem 4: |s(0)|/eta${base ? " (page: 19.8 s)" : ""}`);
    }
  }
}
emit("sign(nomoto(3).bound - nomoto(3).reach)", 1, "abs:0", "ch08 Theorem 4: the surface is reached before the bound");
err("nomoto(5)", "ctrl must be 1..4");

section("dynamics: the logistic map (ch08 page lyapExp(), period(), superstable(), feigenbaum())");
for (const r of [3.83, 3.7, 4, 3.2, 3.5, 2.9]) {
  emit(`lyapexp(${r})`, ch08.lyapExp(r), "rel:1e-11", `ch08 page lyapExp(${r})`);
  emit(`logperiod(${r})`, ch08.period(r) ?? 0, "abs:0", `ch08 page period(${r})${ch08.period(r) ? "" : ": chaotic, none found"}`);
}
emit("lyapexp(4, 200000)", Math.log(2), "abs:1e-3", "ch08 Theorem 5: lambda = ln 2 at r = 4 (page: 0.69314)");
emit("sign(-0.3 - lyapexp(3.83))", 1, "abs:0", "phase9 claim 6: negative inside the period-3 window");
emit("lyapexp(3.2, 1000)", ch08.lyapExp(3.2, 1000), "rel:1e-11", "ch08 page lyapExp(3.2, 1000)");
const FE = ch08.feigenbaum();
for (let n = 0; n <= 8; n++) emit(`superstab(${n})`, FE.R[n], "rel:1e-12", `ch08 page feigenbaum(): R_${n}`);
for (let n = 2; n <= 8; n++) emit(`feigen(${n})`, FE.ratios[n - 2], "rel:1e-7", `ch08 page feigenbaum(): ratio ${n - 2}`);
{ // two more levels with the page's own superstable() and guess rule
  const R = FE.R.slice(); for (let n = 9; n <= 10; n++) R.push(ch08.superstable(n, R[R.length - 1] + (R[R.length - 1] - R[R.length - 2]) / 4.669));
  emit("superstab(10)", R[10], "rel:1e-12", "ch08 page superstable(10)");
  emit("feigen(10)", (R[9] - R[8]) / (R[10] - R[9]), "rel:1e-4", "ch08 page superstable(): ratio at n = 10");
}
emit("feigen(8)", 4.6692, "abs:1e-3", "phase9 claim 6: Feigenbaum ratios -> 4.6692");
emit("superstab(8)", 3.569939, "abs:1e-5", "phase9 claim 6: accumulation r = 3.569939");
err("superstab(13)", "n <= 12");
err("feigen(1)", "n >= 2");
err("lyapexp(4, 0)", "n >= 1");
for (const [x0, n] of [[0.3, 5], [0.1, 10], [0.9, 3]]) emit(`logmap(2, ${x0}, ${n})`, (1 - Math.pow(1 - 2 * x0, 2 ** n)) / 2, "rel:1e-12", "closed form at r = 2: x_n = (1 - (1 - 2 x0)^(2^n))/2");
emit("logmap(4, sin(pi/7)^2, 3)", Math.sin(8 * Math.PI / 7) ** 2, "rel:1e-12", "ch08 Theorem 5: x = sin^2(pi theta), theta doubles each step");
emit("logmap(3.83, 0.3, 0)", 0.3, "abs:0", "zero steps");
err("logmap(4, 0.3, -1)", "n must be a whole number");

section("dynamics: the damped pendulum (ch08 Theorem 2, LaSalle; page rk4())");
for (const [th0, b, t] of [[0.5, 0.3, 60], [1.5, 0.3, 60], [2.8, 0.3, 60], [1, 0.1, 20], [3, 0, 10]]) {
  const f = z => [z[1], -Math.sin(z[0]) - b * z[1]]; let z = [th0, 0];
  for (let i = 0; i < Math.round(t / 0.01); i++) z = ch08.rk4(f, z, 0.01);
  emit(`pendulum(${th0}, ${b}, ${t}).th`, z[0], "abs:1e-11", "ch08 page rk4() on th'' = -sin th - b th'");
  emit(`pendulum(${th0}, ${b}, ${t}).w`, z[1], "abs:1e-11", "ch08 page rk4()");
  if (b === 0) emit(`pendulum(${th0}, 0, ${t}).E`, 1 - Math.cos(th0), "rel:1e-6", "no damping: energy is conserved (RK4 drift < 1e-6)");
}
emit("sign(1e-6 - pendulum(2.8).E)", 1, "abs:0", "phase9 claim 3: converges to the bottom although Vdot is only semidefinite");
err("pendulum(1, 0.3, -1)", "t must be >= 0");
err("pendulum(1, 0.3, 2000)", "t <= 1000 s (1e5 RK4 steps)");

/* ====================== numpy part ====================== */
const here = path.dirname(fileURLToPath(import.meta.url));
const py = execFileSync("python3", ["-I", path.join(here, "gen_control.py"), wii], { encoding: "utf8" });
section("numpy: ch03 tracker, Hamiltonian halves, ch08 Lyapunov and region of attraction (gen_control.py)");
lines.push(...py.split("\n").filter(Boolean));

process.stdout.write("# control library: riccati (index demo, ch03), control (mpc/, approach/), dynamics (ch08). Generated by gen_control.mjs.\n" + lines.join("\n") + "\n");
