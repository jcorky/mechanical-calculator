// Reference cases for LibOps.cpp (queues, inventory, networks, graphs, SLAM).
//   node Tests/reference/gen_ops.mjs /home/user/where-it-isnt > Tests/reference/ops.tsv
// Expected values come from Where It Isn't's own modules (logistics/logistics.js, slam/slam.js). Chapter 9
// has no JS module: its cases come from gen_ops.py (numpy plus the functions of
// verification/python/phase10_graphs.py), which this script runs and appends.
import path from "node:path";
import { pathToFileURL, fileURLToPath } from "node:url";
import { execFileSync } from "node:child_process";

const wii = process.argv[2] || "/home/user/where-it-isnt";
await import(pathToFileURL(path.join(wii, "logistics/logistics.js")).href);
await import(pathToFileURL(path.join(wii, "slam/slam.js")).href);
const L = globalThis.Logistics, S = globalThis.Slam;

const out = [];
const num = x => (x === Infinity ? "inf" : x === -Infinity ? "-inf" : String(x));
const mat = M => "[" + M.map(r => r.map(num).join(", ")).join("; ") + "]";
const vec = v => "[" + v.map(num).join(", ") + "]";
const col = v => "[" + v.map(num).join("; ") + "]";
function emit(expr, expected, tol, prov) {
  if (expected !== "error" && !Number.isFinite(expected) && !(expected === Infinity || expected === -Infinity)) throw new Error(`bad expected for ${expr}: ${expected}`);
  const t = tol || (expected !== "error" && Math.abs(expected) < 1e-12 ? "abs:1e-12" : "rel:1e-9");
  out.push([expr, expected === "error" ? "error" : num(expected), expected === "error" ? "abs:0" : t, prov].join("\t"));
}
const err = (expr, prov) => emit(expr, "error", null, prov);
const fact = k => { let f = 1; for (let i = 2; i <= k; i++) f *= i; return f; };

/* ---------- queues ---------- */
let P = "logistics.js erlangC";
emit("erlangc(3, 2.4)", L.erlangC(3, 2.4), null, P + " (ch13 Thm 7: 0.6472)");
emit("erlangc(1, 0.8)", L.erlangC(1, 0.8), null, P);
emit("erlangc(10, 7.5)", L.erlangC(10, 7.5), null, P);
emit("erlangc(5, 6)", L.erlangC(5, 6), null, P + " (a >= c: always wait)");
err("erlangc(2.5, 1)", "servers must be whole");
P = "closed form a^c/c! / sum_k a^k/k! (Erlang 1917)";
const erlB = (c, a) => { let s = 0; for (let k = 0; k <= c; k++) s += Math.pow(a, k) / fact(k); return Math.pow(a, c) / fact(c) / s; };
emit("erlangb(3, 2.4)", erlB(3, 2.4), null, P);
emit("erlangb(10, 7)", erlB(10, 7), null, P);
emit("erlangb(1, 1)", 0.5, null, P);
emit("erlangb(30, 25)", erlB(30, 25), null, P);
err("erlangb(0, 1)", "servers must be >= 1");
P = "logistics.js mmc";
const Q3 = L.mmc(2.4, 1, 3);
emit("mmc(2.4, 1, 3).W", Q3.W, null, P + " (ch13: W = 2.0787 min)");
emit("mmc(2.4, 1, 3).Wq", Q3.Wq, null, P);
emit("mmc(2.4, 1, 3).C", Q3.C, null, P);
emit("mmc(2.4, 1, 3).L", Q3.L, null, P);
emit("mmc(2.4, 1, 3).Lq", Q3.Lq, null, P);
emit("mmc().rho", Q3.rho, null, P + " (bench defaults)");
const Q1 = L.mmc(0.8, 1, 1);
emit("mmc(0.8, 1, 1).W", Q1.W, null, P + " (ch13: M/M/1 W = 5)");
emit("mmc(0.8, 1, 1).L", Q1.L, null, P);
emit("mmc(7, 0.5, 20).Wq", L.mmc(7, 0.5, 20).Wq, null, P);
err("mmc(3, 1, 3)", "rho = 1: no steady state");
P = "logistics/index.html tail C exp(-(c mu - lambda) t)";
emit("mmctail(2.4, 1, 3, 1)", Q3.C * Math.exp(-(3 - 2.4) * 1), null, P + " (bench: 35.52%)");
emit("mmctail(0.8, 1, 1, 5)", Q1.C * Math.exp(-(1 - 0.8) * 5), null, P);
emit("mmctail(2.4, 1, 3, 0)", Q3.C, null, P + " (t = 0 gives Erlang C)");
err("mmctail(4, 1, 3, 1)", "unstable");
P = "ch13 Thm 6 L = lambda W";
emit("little(0.8, 5)", 0.8 * 5, null, P + " (M/M/1: L = 4)");
emit("little(2.4, mmc(2.4, 1, 3).W)", Q3.L, null, "logistics.js mmc L = lambda W");
P = "logistics.js pk";
emit("pk(0.8, 1, 1).Wq", L.pk(0.8, 1, "exp").Wq, null, P + " exp (ch13: M/M/1 Wq = 4)");
emit("pk(0.8, 1, 0).Wq", L.pk(0.8, 1, "det").Wq, null, P + " det (ch13: 2)");
emit("pk(0.8, 1, 1/3).Wq", L.pk(0.8, 1, "uni").Wq, null, P + " uniform (ch13: 2.667)");
emit("pk(0.8, 1, 1/3).W", L.pk(0.8, 1, "uni").W, null, P);
emit("pk(0.5, 2, 0).W", L.pk(0.5, 2, "det").W, null, P);
emit("pk(0.6, 1, 1).L", 0.6 * L.pk(0.6, 1, "exp").W, null, P + " with Little");
err("pk(1.2, 1)", "rho > 1");
P = "phase14_logistics.py Kingman formulas";
emit("kingman(0.8, 1, 0.5, 1/3).bound", 0.8 * (1 / (2 * 0.8 ** 2) + 1 / 3) / (2 * 0.2), null, P + " (ch13: 2.229)");
emit("kingman(0.8, 1, 0.5, 1/3).approx", 0.8 / 0.2 * ((0.5 + 1 / 3) / 2) * 1, null, P + " (ch13: 1.667)");
emit("kingman(0.8, 1, 1, 1).approx", L.pk(0.8, 1, "exp").Wq, null, "Kingman approx is exact for M/M/1 (logistics.js pk)");
emit("kingman(0.5, 1, 1, 0).bound", 0.5 * (1 / 0.25 + 0) / (2 * 0.5), null, P);
err("kingman(1, 1)", "rho = 1");

/* ---------- inventory ---------- */
P = "logistics.js eoq";
emit("eoq()", L.eoq(50, 1200, 3), null, P + " (ch13: Q* = 200)");
emit("eoq(50, 4800, 3)", L.eoq(50, 4800, 3), null, P + " (4x demand doubles Q*)");
emit("eoq(12, 365, 0.7)", L.eoq(12, 365, 0.7), null, P);
err("eoq(50, 1200, 0)", "h = 0");
P = "logistics.js eoqCost";
emit("eoqcost(200)", L.eoqCost(200, 50, 1200, 3), null, P + " (ch13: 600)");
emit("eoqcost(400)", L.eoqCost(400, 50, 1200, 3), null, P + " (ch13: 750)");
emit("eoqcost(100)", L.eoqCost(100, 50, 1200, 3), null, P + " (ch13: 750)");
emit("eoqcost(eoq(12, 365, 0.7), 12, 365, 0.7)", L.eoqCost(L.eoq(12, 365, 0.7), 12, 365, 0.7), null, P);
err("eoqcost(0)", "Q = 0");
P = "logistics.js newsvendor crit";
emit("fractile()", L.newsvendor().crit, null, P + " (ch13: 2/3)");
emit("fractile(20, 4, 1)", L.newsvendor(20, 4, 1).crit, null, P);
err("fractile(1, 4, 1)", "price = salvage");
P = "logistics.js newsvendor";
const NV = L.newsvendor();
emit("newsvendor().q", NV.qstar, null, P + " (ch13: q* = 43)");
emit("newsvendor().qbest", NV.qbrute, null, P + " brute force");
emit("newsvendor().profit", NV.profit(NV.qstar), null, P + " (ch13: 219.075)");
emit("newsvendor().crit", NV.crit, null, P);
for (const [pr, mean] of [[15, 40], [6, 30], [20, 55]]) {
  const N = L.newsvendor(pr, 4, 1, mean);
  emit(`newsvendor(${pr}, 4, 1, ${mean}).q`, N.qstar, null, P);
  emit(`newsvendor(${pr}, 4, 1, ${mean}).qbest`, N.qbrute, null, P);
  emit(`newsvendor(${pr}, 4, 1, ${mean}).profit`, N.profit(N.qstar), null, P);
}
err("newsvendor(10, 4, 1, 100)", "mean beyond the 150-term pmf");
P = "logistics.js newsvendor profit(q)";
emit("nvprofit(40)", NV.profit(40), null, P);
emit("nvprofit(46)", NV.profit(46), null, P);
emit("nvprofit(0)", NV.profit(0), null, P);
emit("nvprofit(25, 15, 4, 1, 40)", L.newsvendor(15, 4, 1, 40).profit(25), null, P);
P = "logistics.js ssDP";
const D0 = L.ssDP();
emit("ss().s", D0.s, null, P + " (ch13: order up to 15 below 3)");
emit("ss().S", D0.S, null, P);
emit("ss().v0", D0.v0, "rel:1e-9", P + " (ch13: 220.984)");
for (const K of [0, 16, 64]) {
  const D = L.ssDP({ K });
  emit(`ss(${K}).s`, D.s, null, P + ` K = ${K}`);
  emit(`ss(${K}).S`, D.S, null, P + ` K = ${K} (ch13: ${K === 0 ? "base stock 6" : K === 16 ? "(3, 12)" : "(2, 20)"})`);
  emit(`ss(${K}).v0`, D.v0, null, P);
}
const D2 = L.ssDP({ K: 20, p: 4, h: 2, c: 1.5, alpha: 0.8 });
emit("ss(20, 4, 2, 1.5, 0.8).s", D2.s, null, P);
emit("ss(20, 4, 2, 1.5, 0.8).S", D2.S, null, P);
emit("ss(20, 4, 2, 1.5, 0.8).v0", D2.v0, null, P);
err("ss(32, 9, 1, 1, 1)", "alpha = 1");
P = "logistics.js ssEval";
emit("sscost(3, 15)", L.ssEval(3, 15, D0)[-L.INV.LO], null, P + " (ch13: 220.984 = V(0))");
emit("sscost(2, 12)", L.ssEval(2, 12, D0)[-L.INV.LO], null, P);
emit("sscost(3, 15, 32, 9, 1, 1, 0.9, 10)", L.ssEval(3, 15, D0)[10 - L.INV.LO], null, P + " from stock 10");
emit("sscost(5, 18, 20, 4, 2, 1.5, 0.8)", L.ssEval(5, 18, D2)[-L.INV.LO], null, P);
err("sscost(15, 3)", "s > S");

/* ---------- networks ---------- */
const flowCap = (E, n) => { const C = [...Array(n)].map(() => new Array(n).fill(0)); E.forEach(([u, v, c]) => { C[u][v] += c; }); return C; };
const FC = flowCap(L.FLOW_EDGES, 6), MF = L.maxflow(6, L.FLOW_EDGES, 0, 5);
P = "logistics.js maxflow";
emit(`maxflow(${mat(FC)}).flow`, MF.value, null, P + " (ch13 Thm 1: 15)");
emit(`maxflow(${mat(FC)}).aug`, MF.aug, null, P + " (ch13: 4 augmenting paths)");
emit(`maxflow(${mat(FC)}).ns`, MF.reach.length, null, P + " source side of the min cut");
emit(`maxflow(${mat(FC)}, 1, 6).flow`, L.minCutBrute(6, L.FLOW_EDGES, 0, 5).value, null, "logistics.js minCutBrute (16 cuts)");
for (const cap of [0, 5, 20]) {
  const E = L.FLOW_EDGES.map(e => (e[0] === 0 && e[1] === 1 ? [0, 1, cap] : e));
  emit(`maxflow(${mat(flowCap(E, 6))}).flow`, L.maxflow(6, E, 0, 5).value, null, P + ` port lane ${cap}`);
}
emit(`maxflow(${mat(FC)}, 2, 5).flow`, L.maxflow(6, L.FLOW_EDGES, 1, 4).value, null, P + " hub A to warehouse D");
emit(`flowmat(${mat(FC)})(4, 6)`, MF.F[3][5], null, P + " flow on warehouse C -> city");
emit(`flowmat(${mat(FC)})(1, 2)`, MF.F[0][1], null, P);
emit(`flowmat(${mat(FC)})(5, 6)`, MF.F[4][5], null, P + " flow on warehouse D -> city");
err(`maxflow(${mat(FC)}, 1, 1)`, "source = sink");
err("maxflow([0, -1; 0, 0])", "negative capacity");

P = "logistics.js transport";
const T0 = L.transport(L.SUP, L.DEM, L.COST), targs = (s, d, C) => `${vec(s)}, ${vec(d)}, ${mat(C)}`;
emit(`transport(${targs(L.SUP, L.DEM, L.COST)})`, T0.cost, null, P + " (ch13 Thm 3: 415)");
T0.u.forEach((u, i) => emit(`tprices(${targs(L.SUP, L.DEM, L.COST)})(${i + 1})`, u, null, P + " u (ch13: 0, 1, 3)"));
T0.v.forEach((v, j) => emit(`tprices(${targs(L.SUP, L.DEM, L.COST)})(${4 + j})`, v, null, P + " v (ch13: 1, 4, 0, 4)"));
emit(`tplan(${targs(L.SUP, L.DEM, L.COST)})(3, 3)`, T0.X[2][2], null, P);
emit(`tplan(${targs(L.SUP, L.DEM, L.COST)})(2, 1)`, T0.X[1][0], null, P);
emit(`sum(tplan(${targs(L.SUP, L.DEM, L.COST)})' * ones(3, 1))`, L.SUP.reduce((a, b) => a + b, 0), null, P + " total shipped");
for (const [c31, d1] of [[2, 30], [3, 40], [12, 10], [4, 50]]) {
  const T = L.transportCase(c31, d1);
  emit(`transport(${targs(T.sup, T.dem, T.C)})`, T.R.cost, null, `logistics.js transportCase(${c31}, ${d1})`);
  emit(`tprices(${targs(T.sup, T.dem, T.C)})(3)`, T.R.u[2], null, `logistics.js transportCase(${c31}, ${d1}) u3`);
  emit(`tplan(${targs(T.sup, T.dem, T.C)})(3, 1)`, T.R.X[2][0], null, `logistics.js transportCase(${c31}, ${d1}) lane depot 3 -> store 1`);
}
const OTh = L.otCase(0.5);
emit(`transport(${targs(L.OT.a, L.OT.b, OTh.C)})`, OTh.lp, null, "logistics.js otCase(0.5) Kantorovich LP (ch13: 111.237)");
err(`transport([1, 2], [1, 1], [1, 2; 3, 4])`, "unbalanced");
err(`transport([1, 1], [1, 1], [1, 2, 3; 3, 4, 5])`, "cost shape");

P = "logistics.js transport with unit supplies (assignment) and brute force over permutations";
{
  const g = L.rng(21);
  for (let trial = 0; trial < 3; trial++) {
    const n = 4 + trial, C = [...Array(n)].map(() => [...Array(n)].map(() => Math.round(10000 * g()) / 1000));
    const ones = new Array(n).fill(1), R = L.transport(ones, ones, C);
    let brute = Infinity, best = null;
    for (const q of L.permutations([...Array(n).keys()])) { const c = q.reduce((s, j, i) => s + C[i][j], 0); if (c < brute) { brute = c; best = q; } }
    if (Math.abs(brute - R.cost) > 1e-9) throw new Error("assignment check");
    emit(`assign(${mat(C)})`, R.cost, null, P);
    const perm = R.X.map(r => r.findIndex(x => x > 0.5));
    if (perm.some((j, i) => j !== best[i])) throw new Error("perm check");
    emit(`assignperm(${mat(C)})(1)`, perm[0] + 1, null, P);
    emit(`assignperm(${mat(C)})(${n})`, perm[n - 1] + 1, null, P);
  }
}
emit("assign([4, 1, 3; 2, 0, 5; 3, 2, 2])", 5, null, "textbook 3x3 (Kuhn): 1 + 2 + 2");
err("assign([1, 2, 3; 4, 5, 6])", "not square");

P = "logistics.js karp";
for (const [air, har] of [[18, 25], [30, 10], [10, 10], [25, 25]]) {
  const A = L.mpMatrix(air, har), lam = L.karp(A), brute = L.cycleMeans(A)[0].mean;
  if (Math.abs(lam - brute) > 1e-9) throw new Error("karp check");
  emit(`karp(${mat(A)})`, lam, null, P + ` mpMatrix(${air}, ${har})${air === 18 && har === 25 ? " (ch13 Thm 5: 22.5)" : ""} = enumeration of cycles`);
}
emit("karp([1, 5; -inf, 2])", L.karp([[1, 5], [-Infinity, 2]]), null, P + " (node 1 self-loop only from node 1)");
err("karp([1, 2, 3])", "not square");
P = "logistics.js otimes / mpRun";
{
  const A = L.mpMatrix(), X = L.mpRun(A, 20), x0 = [3, 0, 1, 7];
  emit(`mpmul(mppow(${mat(A)}, 12), zeros(4, 1))(1)`, X[12][0], null, P + " x(12) from zero");
  emit(`mpmul(mppow(${mat(A)}, 20), zeros(4, 1))(4)`, X[20][3], null, P + " x(20)");
  emit(`mpmul(${mat(A)}, ${col(x0)})(2)`, L.otimes(A, x0)[1], null, P);
  emit(`mpmul(${mat(A)}, ${col(x0)})(3)`, L.otimes(A, x0)[2], null, P);
  emit(`(mpmul(mppow(${mat(A)}, 212), zeros(4, 1))(1) - mpmul(mppow(${mat(A)}, 200), zeros(4, 1))(1))/12`, (L.mpRun(A, 212)[212][0] - L.mpRun(A, 212)[200][0]) / 12, null, P + " growth rate = lambda");
  const AA = A.map(r => A[0].map((_, j) => Math.max(...r.map((a, k) => a + A[k][j]))));
  emit(`mpmul(${mat(A)}, ${mat(A)})(2, 2)`, L.otimes(A, A.map(r => r[1]))[1], null, P + " column 2 of A x A");
  emit(`mppow(${mat(A)}, 2)(3, 1)`, AA[2][0], null, "max-plus definition (A x A)_31");
  emit(`mppow(${mat(A)}, 0)(1, 2)`, -Infinity, null, "max-plus identity off the diagonal");
  emit(`mppow(${mat(A)}, 0)(2, 2)`, 0, null, "max-plus identity on the diagonal");
  err(`mpmul([1, 2], [1, 2])`, "inner dimensions");
}

P = "logistics.js w1cdf";
{
  const g = L.rng(3), xs = [...Array(9)].map(() => Math.round(1000 * g()) / 100), ys = [...Array(6)].map(() => Math.round(1000 * g()) / 100 + 1);
  emit(`w1(${vec(xs)}, ${vec(ys)})`, L.w1cdf(xs, xs.map(() => 1 / 9), ys, ys.map(() => 1 / 6)), null, P + " empirical measures");
  emit("w1([0, 0, 1], [1, 1, 1])", L.w1cdf([0, 0, 1], [1 / 3, 1 / 3, 1 / 3], [1, 1, 1], [1 / 3, 1 / 3, 1 / 3]), null, P + " (2/3)");
  emit("w1([1, 2, 3], [1, 2, 3] + 0.5)", 0.5, null, "shift by 0.5");
  emit(`w1w(${vec(L.OT.x)}, ${vec(L.OT.a)}, ${vec(L.OT.y)}, ${vec(L.OT.b)})`, L.w1cdf(L.OT.x, L.OT.a, L.OT.y, L.OT.b), null, P + " (ch13 Thm 4: 288)");
  const y5 = L.OT.y.map(v => v + 5);
  emit(`w1w(${vec(L.OT.x)}, ${vec(L.OT.a)}, ${vec(y5)}, ${vec(L.OT.b)})`, L.w1cdf(L.OT.x, L.OT.a, y5, L.OT.b), null, P + " stores shifted 5 km");
  err("w1w([1, 2], [1], [1], [1])", "mass length");
}
P = "logistics.js otCase monotone coupling";
for (const [p, shift] of [[1, 0], [1.5, 0], [2, 0], [3, 0], [0.5, 0], [2, 4]]) {
  const O = L.otCase(p, shift);
  emit(`otmono(${vec(L.OT.x)}, ${vec(L.OT.a)}, ${vec(O.y)}, ${vec(L.OT.b)}, ${p})`, O.mono, null, P + ` p = ${p}, shift ${shift}${shift === 0 ? " (ch13: 288, 779.11, 2256, 21744; p = 0.5 115.260)" : ""}`);
  if (p >= 1) emit(`transport(${targs(L.OT.a, L.OT.b, O.C)})`, O.lp, null, `logistics.js otCase(${p}, ${shift}) LP = monotone`);
}
emit("otmono([2, 0], [1, 1], [3, 1], [1, 1], 2)", 2, null, "sorting the inputs: (0->1) + (2->3)");
err("otmono([0, 1], [1, 1], [0, 1], [1, 2])", "masses unequal");

P = "logistics.js kruskal";
{
  const R = [...Array(6)].map(() => new Array(6).fill(Infinity));
  L.ROADS.forEach(([u, v, w]) => { R[u][v] = w; R[v][u] = w; });
  emit(`mst(${mat(R)})`, L.kruskal(6, L.ROADS).reduce((s, e) => s + e[2], 0), null, P + " ROADS (ch13 Thm 11: 25)");
  const Rd = R.map(r => r.slice()); Rd[0][1] = Rd[1][0] = Rd[0][2] = Rd[2][0] = Infinity;
  err(`mst(${mat(Rd)})`, "depot 1 cut off: disconnected");
  err("mst([0, 1; 2, 0])", "not symmetric");
}
P = "logistics.js tspCase (points, dmat, tspBrute, christofides, heldKarp)";
for (const [seed, n] of [[1, 8], [2, 6], [5, 9], [7, 5]]) {
  const T = L.tspCase(seed, n), Pm = mat(T.P), D = mat(T.D);
  emit(`distmat(${Pm})(1, 2)`, T.D[0][1], null, "logistics.js dmat");
  emit(`tsp(distmat(${Pm}))`, T.opt, null, P + ` seed ${seed}, n ${n}${seed === 1 ? " (ch13: OPT 2.897481)" : ""}`);
  emit(`tsp(${D})`, T.opt, null, P);
  emit(`christofides(${D}).len`, T.chr.len, null, P);
  emit(`christofides(${D}).mst`, T.chr.mst, null, P);
  emit(`christofides(${D}).match`, T.chr.matchCost, null, P);
  emit(`mst(${D})`, T.mst, null, P + " Prim = Kruskal");
  emit(`hkbound(${D}).lb`, T.hk, "rel:1e-8", P + " subgradient");
  emit(`hkbound(${D}).onetree`, T.oneTree, null, P);
  emit(`tourlen(${D}, ${vec(T.tour.map(i => i + 1))})`, L.tourLen(T.D, T.tour), null, "logistics.js tourLen of the optimal tour");
}
emit("tsp([0, 1, inf; 1, 0, 1; inf, 1, 0])", Infinity, null, "no closed tour without the 1-3 edge");
err("christofides([0, 1; 1, 0])", "too few stops");
err("tourlen([0, 1; 1, 0], [1, 3])", "stop out of range");

/* ---------- graphs: spanning trees from the logistics module ---------- */
{
  const A = [...Array(6)].map(() => new Array(6).fill(0));
  L.ROADS.forEach(([u, v]) => { A[u][v] = 1; A[v][u] = 1; });
  const M = L.matroidCase();
  emit(`spantrees(${mat(A)})`, M.count, null, "logistics.js matroidCase: 114 spanning trees by enumeration (ch13 Thm 11)");
  emit(`spantrees(${mat(A)})`, M.kirchhoff, null, "logistics.js matroidCase Kirchhoff determinant");
}

/* ---------- SLAM ---------- */
P = "slam.js gammp";
for (const [a, x] of [[2.5, 1.3], [1, 2], [0.5, 4], [10, 7.2], [3, 30]]) emit(`gammp(${a}, ${x})`, S.gammp(a, x), null, P);
emit("gammp(1, 2)", 1 - Math.exp(-2), "rel:1e-9", "closed form 1 - exp(-x)");
emit("gammp(2, 0)", 0, null, P);
err("gammp(0, 1)", "a = 0");
P = "slam.js chi2cdf";
for (const [x, k] of [[7.814727903251178, 3], [5.991464547107979, 2], [1, 1], [130, 120]]) emit(`chi2cdf(${x}, ${k})`, S.chi2cdf(x, k), null, P);
P = "slam.js chi2inv";
for (const [p, k] of [[0.95, 3], [0.99, 2], [0.025, 120], [0.975, 120], [0.5, 1], [0.025, 900], [0.975, 900]]) emit(`chi2inv(${p}, ${k})`, S.chi2inv(p, k), null, P);
err("chi2inv(1, 3)", "p = 1");
P = "slam.js gate2";
emit("chi2gate()", S.gate2(0.95), null, P + " (5.991)");
emit("chi2gate(0.99)", S.gate2(0.99), null, P);
emit("chi2gate(0.95, 3)", S.chi2inv(0.95, 3), null, "slam.js chi2inv");
emit("chi2gate(0.95, 2) - chi2inv(0.95, 2)", 0, "abs:1e-9", "closed form = bisection");
P = "slam.js neesBand";
for (const M of [40, 300, 1]) {
  const [lo, hi] = S.neesBand(M);
  emit(`neesband(${M}).lo`, lo, null, P);
  emit(`neesband(${M}).hi`, hi, null, P);
}
emit("neesband(10, 2, 0.9).hi", S.chi2inv(0.95, 20) / 10, null, "slam.js chi2inv");
err("neesband(0)", "no runs");
P = "slam.js daMC innovation covariance";
{
  const R = S.daMC(1, 3), g = S.rng(5);
  for (const c of R.cand) {
    const nu = [0.07 * g.g(), 0.03 * g.g()], d2 = nu[0] * (c.Si[0][0] * nu[0] + c.Si[0][1] * nu[1]) + nu[1] * (c.Si[1][0] * nu[0] + c.Si[1][1] * nu[1]);
    emit(`mahal(${vec(nu)}, ${mat(c.S)})`, d2, null, P + " d2 = nu' S^-1 nu");
  }
  const Pc = [[0.04, 0.01, 0.002], [0.01, 0.09, -0.003], [0.002, -0.003, 0.0025]], e = [0.1, -0.2, 0.03], Pi = S.inv(Pc);
  emit(`nees(${col(e)}, ${mat(Pc)})`, e.reduce((s, ei, i) => s + ei * Pi[i].reduce((t, p, j) => t + p * e[j], 0), 0), null, "slam.js inv: e' P^-1 e");
  err("mahal([1, 2], [1, 0, 0; 0, 1, 0; 0, 0, 1])", "size mismatch");
}
P = "slam.js voPair + kabsch";
for (const [seed, n, sig, yaw] of [[11, 30, 1, 5], [3, 8, 0, -12], [4, 5, 2.5, 18]]) {
  const pr = S.voPair(S.rng(seed), n, sig, yaw), K = S.kabsch(pr.B, pr.A), args = `${mat(pr.B)}, ${mat(pr.A)}`;
  const tag = `${P} (rng ${seed}, n ${n}, ${sig} px, ${yaw} deg)`;
  emit(`kabsch(${args}).th`, K.th, sig === 0 ? "rel:1e-9" : null, tag);
  emit(`kabsch(${args}).tx`, K.t[0], null, tag);
  emit(`kabsch(${args}).ty`, K.t[1], null, tag);
  emit(`kabsch(${args}).cost`, S.cost(pr.B, pr.A, K.th, K.t), sig === 0 ? "abs:1e-20" : null, tag + " slam.js cost");
  emit(`kabsch(${args}).S`, K.S, null, tag);
  if (sig > 0) emit(`headvar(${sig * S.CAM.h / S.CAM.f}, kabsch(${args}).S)`, S.headingVar(sig * S.CAM.h / S.CAM.f, K.S), null, "slam.js headingVar");
}
emit("kabsch([0, 0; 1, 0; 0, 1], [1, 1; 1, 2; 0, 1]).th", Math.PI / 2, null, "exact rotation by 90 deg");
emit("kabsch([0, 0; 1, 0; 0, 1], [1, 1; 1, 2; 0, 1]).tx", 1, null, "exact translation (1, 1)");
err("kabsch([0, 0; 1, 0], [0, 0; 1, 0; 2, 0])", "different n");
P = "slam.js headingVar x N (random walk, sqrt N)";
emit("drift(0.004, 2.5, 100).var", 100 * S.headingVar(0.004, 2.5), null, P);
emit("drift(0.004, 2.5, 100).rms", Math.sqrt(100 * S.headingVar(0.004, 2.5)), null, P);
emit("drift(0.004, 2.5, 400).rms / drift(0.004, 2.5, 100).rms", 2, null, "sqrt(400/100)");
emit("headvar(0.004, 2.5)", S.headingVar(0.004, 2.5), null, "slam.js headingVar");
err("drift(0.004, 0, 10)", "S = 0");
P = "slam.js voWalk (mulberry32 + Box-Muller, bit for bit)";
for (const [N, n, sig, yaw, M, seed] of [[100, 20, 1, 3, 200, 2], [30, 10, 0.5, 3, 50, 7], [12, 6, 2, -4, 40, 9]]) {
  const W = S.voWalk(seed, M, N, n, sig, yaw), a = `vowalk(${N}, ${n}, ${sig}, ${yaw}, ${M}, ${seed})`;
  emit(`${a}.rmsh`, W.rmsH[N - 1], "rel:1e-7", P + (seed === 2 ? " bench defaults" : ""));
  emit(`${a}.pred`, W.pred[N - 1], "rel:1e-7", P);
  emit(`${a}.rmsp`, W.rmsP[N - 1], "rel:1e-7", P);
}
err("vowalk(1000, 200, 1, 3, 2000, 1)", "too much work");

/* ---------- Chapter 9 (numpy) ---------- */
const here = path.dirname(fileURLToPath(import.meta.url));
const py = execFileSync("python3", [path.join(here, "gen_ops.py"), wii], { encoding: "utf8" });

process.stdout.write("# LibOps reference cases: generated by Tests/reference/gen_ops.mjs from logistics/logistics.js and slam/slam.js,\n");
process.stdout.write("# plus Tests/reference/gen_ops.py (numpy, verification/python/phase10_graphs.py) for Chapter 9.\n");
process.stdout.write(out.join("\n") + "\n" + py);
