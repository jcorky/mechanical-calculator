// Reference cases for LibOrbits.cpp (categories orbits and ins), computed by Where It Isn't's own modules.
//   node Tests/reference/gen_orbits.mjs /home/user/where-it-isnt > Tests/reference/orbits.tsv
// Expected values come from globalThis.Orbits (orbits/orbits.js) and globalThis.Ins (ins/ins.js). Where a
// calculator function has no counterpart in the modules (true anomaly, escape speed, the dimensional
// bi-elliptic burns, the 6 x 6 CW cross-track block, the Allan minimum, the Kalman steady state), the value is
// computed here independently from the closed form or by iterating the module's own matrices to convergence.
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const wii = process.argv[2] || '/home/user/where-it-isnt';
await import(pathToFileURL(path.join(wii, 'orbits/orbits.js')).href);
await import(pathToFileURL(path.join(wii, 'ins/ins.js')).href);
const O = globalThis.Orbits, I = globalThis.Ins;

const lines = [];
const num = x => (Number.isFinite(x) ? String(x) : x > 0 ? 'inf' : x < 0 ? '-inf' : 'nan');
const vec = a => '[' + a.map(String).join(', ') + ']';
const mat = (P, n = 3) => '[' + [...Array(n)].map((_, i) => P.slice(n * i, n * i + n).map(String).join(', ')).join('; ') + ']';
function add(expr, expected, tol, prov) {
  if (expected !== 'error' && !Number.isFinite(expected) && !(expected === Infinity)) throw new Error(`bad expected for ${expr}: ${expected}`);
  lines.push([expr, expected === 'error' ? 'error' : num(expected), tol, prov].join('\t'));
}
const err = (expr, why) => add(expr, 'error', 'abs:0', why);
const MU = O.MU, RE = O.RE, TAU = 2 * Math.PI;

/* ---------------- orbits: basic two-body relations ---------------- */
lines.push('# orbits: period, circular / escape speed, vis-viva, mean motion, energy, Stumpff');
for (const a of [10000e3, 42164e3]) add(`orbperiod(${a})`, O.keplerState(a, 0, 0).period, 'rel:1e-13', 'Orbits.keplerState().period');
add('orbperiod(1, 1)', TAU, 'rel:1e-14', 'closed form 2 pi sqrt(a^3/mu)');
err('orbperiod(-1)', 'a must be positive');
{
  const H = O.hohmann(6678e3, 42164e3);
  add('vcirc(6678e3)', H.v1, 'rel:1e-14', 'Orbits.hohmann(6678e3, 42164e3).v1');
  add('vcirc(42164e3)', H.v2, 'rel:1e-14', 'Orbits.hohmann(6678e3, 42164e3).v2');
  add('vesc(6678e3)', Math.SQRT2 * H.v1, 'rel:1e-14', 'closed form sqrt(2) v_circ');
  add(`vesc(${RE})`, Math.sqrt(2 * MU / RE), 'rel:1e-14', 'closed form sqrt(2 mu/R_E)');
  err('vesc(0)', 'r must be positive');
  add(`visviva(6678e3, ${H.at})`, H.vp, 'rel:1e-14', 'Orbits.hohmann().vp (periapsis of the transfer ellipse)');
  add(`visviva(42164e3, ${H.at})`, H.va, 'rel:1e-14', 'Orbits.hohmann().va (apoapsis of the transfer ellipse)');
  err('visviva(3e7, 1e7)', 'r beyond the apoapsis 2a');
}
for (const r of [RE + 400e3, 7000e3]) add(`meanmotion(${r})`, O.meanMotion(r), 'rel:1e-14', 'Orbits.meanMotion');
add('orbenergy([7000e3, 0], [0, 7500])', O.energy([7000e3, 0], [0, 7500]), 'rel:1e-14', 'Orbits.energy');
add('orbenergy([5000e3, 6000e3], [-3000, 9000])', O.energy([5000e3, 6000e3], [-3000, 9000]), 'rel:1e-14', 'Orbits.energy');
err('orbenergy([1, 2, 3], [0, 1])', 'r must be a 2-vector');
for (const z of [1, -1, 0, 1e-9, 30, -50]) {
  add(`stumpc(${z})`, O.stumpC(z), 'rel:1e-14', 'Orbits.stumpC');
  add(`stumps(${z})`, O.stumpS(z), 'rel:1e-14', 'Orbits.stumpS');
}

/* ---------------- Kepler ---------------- */
lines.push('# orbits: Kepler equation and state at time t');
const trueAnom = (E, e) => { // tan(nu/2) = sqrt((1+e)/(1-e)) tan(E/2), unwrapped to E's revolution
  let nu = 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E / 2), Math.sqrt(1 - e) * Math.cos(E / 2));
  while (nu - E > Math.PI) nu -= TAU; while (E - nu > Math.PI) nu += TAU; return nu;
};
for (const [M, e] of [[1, 0.3], [3, 0.9], [0.2, 0.0], [5.5, 0.6]]) {
  const K = O.keplerE(M, e);
  add(`kepler(${M}, ${e}).E`, K.E, 'rel:1e-13', 'Orbits.keplerE');
  add(`kepler(${M}, ${e}).iters`, K.iters, 'abs:0', 'Orbits.keplerE iterations');
  add(`kepler(${M}, ${e}).nu`, trueAnom(K.E, e), 'rel:1e-12', 'closed form tan(nu/2) = sqrt((1+e)/(1-e)) tan(E/2)');
}
err('kepler(1, 1)', 'parabolic e = 1 rejected');
add('trueanom(2, 0.5)', trueAnom(2, 0.5), 'rel:1e-13', 'closed form tan(nu/2) relation');
add('trueanom(4, 0.2)', trueAnom(4, 0.2), 'rel:1e-13', 'closed form tan(nu/2) relation');
err('trueanom(1, -0.1)', 'negative eccentricity');
for (const [a, e, t] of [[10000e3, 0.3, 6000], [26000e3, 0.7, 20000]]) {
  const K = O.keplerState(a, e, t);
  const ex = `kepstate(${a}, ${e}, ${t})`;
  add(`${ex}.x`, K.r[0], 'rel:1e-12', 'Orbits.keplerState r[0]');
  add(`${ex}.y`, K.r[1], 'rel:1e-12', 'Orbits.keplerState r[1]');
  add(`${ex}.vx`, K.v[0], 'rel:1e-12', 'Orbits.keplerState v[0]');
  add(`${ex}.vy`, K.v[1], 'rel:1e-12', 'Orbits.keplerState v[1]');
  add(`${ex}.E`, K.E, 'rel:1e-13', 'Orbits.keplerState E');
  add(`${ex}.M`, K.M, 'rel:1e-13', 'Orbits.keplerState M');
  add(`${ex}.period`, K.period, 'rel:1e-13', 'Orbits.keplerState period');
  add(`${ex}.r`, O.norm(K.r), 'rel:1e-12', 'Orbits.norm(keplerState r)');
}
add('kepstate().x', O.keplerState(10000e3, 0.3, 6000).r[0], 'rel:1e-12', 'bench 1 defaults a = 10000 km, e = 0.3, t = 100 min');
err('kepstate(1e7, 1.2, 100)', 'hyperbolic e rejected');
{
  const s0 = O.periState(10000e3, 0.3), P = O.propagate(s0.r, s0.v, 6000);
  const ex = `twobody(${vec(s0.r)}, ${vec(s0.v)}, 6000)`;
  add(`${ex}.x`, P.r[0], 'abs:1e-3', 'Orbits.propagate (RK4, 2000 steps per local period)');
  add(`${ex}.y`, P.r[1], 'abs:1e-3', 'Orbits.propagate');
  add(`${ex}.vx`, P.v[0], 'abs:1e-6', 'Orbits.propagate');
  add(`${ex}.vy`, P.v[1], 'abs:1e-6', 'Orbits.propagate');
  add(`${ex}.steps`, P.steps, 'abs:0', 'Orbits.propagate step count');
  const K = O.keplerState(10000e3, 0.3, 6000);
  add(`${ex}.x`, K.r[0], 'abs:1', 'independent: Kepler closed form (bench 1 check, |d| < 1e-6 a)');
  const P2 = O.propagate([7000e3, 0], [0, 9000], 20000, { n: 500 });
  add('twobody([7000e3, 0], [0, 9000], 20000, 500).x', P2.r[0], 'abs:1e-3', 'Orbits.propagate n = 500');
  add('twobody([7000e3, 0], [0, 9000], 20000, 500).vy', P2.v[1], 'abs:1e-6', 'Orbits.propagate n = 500');
  err('twobody([0, 0], [0, 1], 10)', 'r = 0');
}

/* ---------------- Hohmann ---------------- */
lines.push('# orbits: Hohmann transfer, propellant, flown transfer');
for (const [r1, r2] of [[6678e3, 42164e3], [7000e3, 6678e3], [6578e3, 384400e3]]) {
  const H = O.hohmann(r1, r2), ex = `hohmann(${r1}, ${r2})`;
  for (const k of ['dv1', 'dv2', 'tof', 'vp', 'va', 'v1', 'v2', 'at', 'et']) add(`${ex}.${k}`, H[k], 'rel:1e-13', `Orbits.hohmann().${k}`);
  add(`${ex}.dv`, H.dv1 + H.dv2, 'rel:1e-13', 'Orbits.hohmann() dv1 + dv2');
  const mr = Math.exp((H.dv1 + H.dv2) / (9.80665 * 320));
  add(`hohprop(${r1}, ${r2}).mr`, mr, 'rel:1e-13', 'bench 2: m0/mf = exp(dv / (g0 * 320 s))');
  add(`hohprop(${r1}, ${r2}).frac`, 1 - 1 / mr, 'rel:1e-12', 'bench 2: propellant share 1 - 1/mr');
}
add('hohmann().dv1', O.hohmann(6678e3, 42164e3).dv1, 'rel:1e-13', 'bench 2 defaults r1 = 6678 km, r2 = 42164 km');
add('hohprop(6678e3, 42164e3, 450).mr', Math.exp((O.hohmann(6678e3, 42164e3).dv1 + O.hohmann(6678e3, 42164e3).dv2) / (9.80665 * 450)), 'rel:1e-13', 'rocket equation, Isp 450 s');
{ // Moon-centred, independent closed form
  const mu = 4.9048695e12, r1 = 1837.4e3, r2 = 1937.4e3, at = (r1 + r2) / 2;
  add(`hohmann(${r1}, ${r2}, ${mu}).dv1`, Math.sqrt(mu * (2 / r1 - 1 / at)) - Math.sqrt(mu / r1), 'rel:1e-12', 'closed form vis-viva with lunar mu');
  add(`hohmann(${r1}, ${r2}, ${mu}).tof`, Math.PI * Math.sqrt(at ** 3 / mu), 'rel:1e-13', 'closed form pi sqrt(a^3/mu)');
}
err('hohmann(0, 1e7)', 'r1 must be positive');
err('hohprop(6678e3, 42164e3, 0)', 'Isp must be positive');
for (const [r1, r2] of [[6678e3, 42164e3], [42164e3, 7000e3]]) {
  const F = O.hohmannFly(r1, r2), ex = `hohfly(${r1}, ${r2})`;
  add(`${ex}.rend`, F.rEnd, 'rel:1e-10', 'Orbits.hohmannFly rEnd');
  add(`${ex}.rext`, F.rExt, 'rel:1e-10', 'Orbits.hohmannFly rExt');
  add(`${ex}.speed`, F.speed, 'rel:1e-10', 'Orbits.hohmannFly speed');
  add(`${ex}.vr`, F.vr, 'abs:1e-6', 'Orbits.hohmannFly radial speed');
  add(`${ex}.rend`, r2, 'rel:1e-6', 'bench 2 check: arrives at r2 within 1e-6');
}

/* ---------------- Lambert, conics ---------------- */
lines.push('# orbits: Lambert by universal variables, conic elements');
const lam = (R1, R2, dt, why) => {
  const L = O.lambert(R1, R2, dt), ex = `lambert(${vec(R1)}, ${vec(R2)}, ${dt})`;
  add(`${ex}.v1x`, L.v1[0], 'rel:1e-9', `Orbits.lambert v1 ${why}`);
  add(`${ex}.v1y`, L.v1[1], 'rel:1e-9', `Orbits.lambert v1 ${why}`);
  add(`${ex}.v2x`, L.v2[0], 'rel:1e-9', `Orbits.lambert v2 ${why}`);
  add(`${ex}.v2y`, L.v2[1], 'rel:1e-9', `Orbits.lambert v2 ${why}`);
  add(`${ex}.z`, L.z, 'abs:1e-9', `Orbits.lambert z ${why}`);
  add(`${ex}.dth`, L.dth, 'rel:1e-14', `Orbits.lambert transfer angle ${why}`);
  add(`${ex}.v1`, O.norm(L.v1), 'rel:1e-9', `|v1| ${why}`);
  return L;
};
{
  const th = 95 * Math.PI / 180, R1 = [7000e3, 0], R2 = [12000e3 * Math.cos(th), 12000e3 * Math.sin(th)];
  const L = lam(R1, R2, 3600, '(bench 3 defaults)');
  const el = O.elements(R1, L.v1), ex = `conic(${vec(R1)}, ${vec(L.v1)})`;
  for (const k of ['a', 'e', 'h', 'p', 'rp']) add(`${ex}.${k}`, el[k], 'rel:1e-12', `Orbits.elements().${k} of the bench-3 transfer`);
  add(`${ex}.ex`, el.ev[0], 'rel:1e-12', 'Orbits.elements().ev[0]');
  add(`${ex}.ey`, el.ev[1], 'rel:1e-12', 'Orbits.elements().ev[1]');
  const th2 = 200 * Math.PI / 180;
  lam([7000e3, 0], [20000e3 * Math.cos(th2), 20000e3 * Math.sin(th2)], 14400, '(long way, 200 deg)');
  lam([8000e3, 3000e3], [-6000e3, 9000e3], 2400, '(fast, hyperbolic z < 0)');
  err('lambert([7000e3, 0], [-12000e3, 0], 3600)', '180 deg transfer is degenerate');
  err('lambert([7000e3, 0, 0], [0, 12000e3], 3600)', 'vectors must be [x, y]');
  err('lambert([7000e3, 0], [0, 12000e3], -5)', 'negative time of flight');
}
{
  const r = [7000e3, 0], v = [0, 11500], el = O.elements(r, v), ex = `conic(${vec(r)}, ${vec(v)})`;
  add(`${ex}.a`, el.a, 'rel:1e-12', 'Orbits.elements hyperbola, a < 0');
  add(`${ex}.e`, el.e, 'rel:1e-12', 'Orbits.elements hyperbola');
  add(`${ex}.energy`, O.energy(r, v), 'rel:1e-13', 'Orbits.energy');
  err('conic([0, 0], [1, 1])', 'r = 0');
}

/* ---------------- Clohessy-Wiltshire ---------------- */
lines.push('# orbits: Clohessy-Wiltshire STM, free drift, rendezvous, shooting');
{
  const rref = RE + 400e3, n = O.meanMotion(rref), per = TAU / n, T = 0.5 * per;
  const X0 = [-1000, -10000, 0, -1.5 * n * -1000];
  const P = O.cwStm(n, 1000);
  for (const [i, j] of [[0, 0], [0, 2], [1, 0], [1, 3], [2, 0], [3, 3], [1, 2]]) add(`cwstm(${n}, 1000)(${i + 1}, ${j + 1})`, P[i][j], 'rel:1e-13', `Orbits.cwStm[${i}][${j}]`);
  const map = [1, 2, 4, 5];
  for (const [i, j] of [[0, 0], [1, 3], [3, 2], [2, 1]]) add(`cwstm(${n}, 1000, 6)(${map[i]}, ${map[j]})`, P[i][j], 'rel:1e-13', `Orbits.cwStm[${i}][${j}] embedded in the 6-state matrix`);
  const s = Math.sin(n * 1000), c = Math.cos(n * 1000);
  add(`cwstm(${n}, 1000, 6)(3, 3)`, c, 'rel:1e-14', 'closed form cross-track cos(nt)');
  add(`cwstm(${n}, 1000, 6)(3, 6)`, s / n, 'rel:1e-14', 'closed form sin(nt)/n');
  add(`cwstm(${n}, 1000, 6)(6, 3)`, -n * s, 'rel:1e-14', 'closed form -n sin(nt)');
  add(`cwstm(${n}, 1000, 6)(1, 3)`, 0, 'abs:0', 'in-plane and cross-track decouple');
  err(`cwstm(${n}, 1000, 5)`, 'dim must be 4 or 6');
  err('cwstm(0, 1000)', 'n must be positive');

  const D = O.mv(O.cwStm(n, T), X0), ex = `cwdrift(${n}, ${T}, ${vec(X0)})`;
  ['x', 'y', 'vx', 'vy'].forEach((k, i) => add(`${ex}.${k}`, D[i], 'abs:1e-8', `Orbits.mv(cwStm, X0) bench 4 defaults [${i}]`));
  const X6 = [X0[0], X0[1], 500, X0[2], X0[3], 0.1], ex6 = `cwdrift(${n}, ${T}, ${vec(X6)})`;
  add(`${ex6}.y`, D[1], 'abs:1e-8', 'Orbits CW in-plane part of the 6-state drift');
  add(`${ex6}.z`, 500 * Math.cos(n * T) + 0.1 * Math.sin(n * T) / n, 'abs:1e-8', 'closed form z(t) = z0 cos nt + vz0 sin nt / n');
  add(`${ex6}.vz`, -500 * n * Math.sin(n * T) + 0.1 * Math.cos(n * T), 'abs:1e-11', 'closed form vz(t)');
  err(`cwdrift(${n}, ${T}, [1, 2, 3])`, 'state must have 4 or 6 elements');

  const NL = O.relNonlinear(rref, X0, T).X, exn = `cwnonlin(${rref}, ${vec(X0)}, ${T})`;
  ['x', 'y', 'vx', 'vy'].forEach((k, i) => add(`${exn}.${k}`, NL[i], 'abs:1e-6', `Orbits.relNonlinear (4000 RK4 steps) [${i}]`));
  const NL2 = O.relNonlinear(rref, X0, T, { steps: 1000 }).X;
  add(`cwnonlin(${rref}, ${vec(X0)}, ${T}, 1000).y`, NL2[1], 'abs:1e-6', 'Orbits.relNonlinear, 1000 steps');

  for (const [X, why] of [[X0, 'bench 4 defaults'], [X0.map(x => x / 2), 'half offset']]) {
    const Rv = O.rendezvous(rref, X, T), exr = `rendezvous(${rref}, ${vec(X)}, ${T})`;
    add(`${exr}.v0x`, Rv.v0[0], 'rel:1e-11', `Orbits.rendezvous v0 ${why}`);
    add(`${exr}.v0y`, Rv.v0[1], 'rel:1e-11', `Orbits.rendezvous v0 ${why}`);
    add(`${exr}.dv1x`, Rv.dv1[0], 'rel:1e-11', `Orbits.rendezvous dv1 ${why}`);
    add(`${exr}.dv2x`, Rv.dv2[0], 'rel:1e-10', `Orbits.rendezvous dv2 ${why}`);
    add(`${exr}.dv2y`, Rv.dv2[1], 'rel:1e-10', `Orbits.rendezvous dv2 ${why}`);
    add(`${exr}.dv1`, Math.hypot(...Rv.dv1), 'rel:1e-11', `|dv1| ${why}`);
    add(`${exr}.dv`, Math.hypot(...Rv.dv1) + Math.hypot(...Rv.dv2), 'rel:1e-10', `|dv1| + |dv2| ${why}`);
    const S = O.cwShoot(rref, X, T), exs = `cwshoot(${rref}, ${vec(X)}, ${T})`;
    add(`${exs}.v0x`, S.v[0], 'rel:1e-8', `Orbits.cwShoot ${why}`);
    add(`${exs}.v0y`, S.v[1], 'rel:1e-8', `Orbits.cwShoot ${why}`);
    add(`${exs}.iters`, S.iters, 'abs:0', `Orbits.cwShoot iterations ${why}`);
    add(`${exs}.corr`, Math.hypot(S.v[0] - Rv.v0[0], S.v[1] - Rv.v0[1]), 'rel:1e-6', `bench 4: shooting correction ${why}`);
  }
  err(`rendezvous(${rref}, ${vec(X0)}, ${per})`, 'a whole orbit is singular');
  err(`rendezvous(${rref}, [1, 2], ${T})`, 'state must be [x, y, vx, vy]');
}

/* ---------------- energy drift ---------------- */
lines.push('# orbits: integrator energy drift (a = 10 000 km)');
{
  const mx = a => Math.max(...a);
  for (const [e, spo, orb] of [[0.3, 500, 50], [0.6, 300, 20]]) {
    const lf = O.drift(e, spo, orb, 'leapfrog'), rk = O.drift(e, spo, orb, 'rk4'), eu = O.drift(e, spo, orb, 'euler'), lf2 = O.drift(e, 2 * spo, orb, 'leapfrog');
    add(`edrift(${e}, ${spo}, ${orb}).max`, mx(lf), 'rel:1e-6', 'Orbits.drift leapfrog, max over orbits');
    add(`edrift(${e}, ${spo}, ${orb}, 1).first`, lf[0], 'rel:1e-6', 'Orbits.drift leapfrog orbit 1');
    add(`edrift(${e}, ${spo}, ${orb}, 1).last`, lf[orb - 1], 'rel:1e-6', 'Orbits.drift leapfrog last orbit');
    add(`edrift(${e}, ${spo}, ${orb}, 2).max`, mx(rk), 'rel:1e-5', 'Orbits.drift RK4');
    add(`edrift(${e}, ${spo}, ${orb}, 0).last`, eu[orb - 1], 'rel:1e-6', 'Orbits.drift Euler');
    const ex = `edriftcmp(${e}, ${spo}, ${orb})`;
    add(`${ex}.lf`, mx(lf), 'rel:1e-6', 'Orbits.drift leapfrog');
    add(`${ex}.rk4`, mx(rk), 'rel:1e-5', 'Orbits.drift rk4');
    add(`${ex}.euler`, mx(eu), 'rel:1e-6', 'Orbits.drift euler');
    add(`${ex}.q`, mx(lf) / mx(lf2), 'rel:1e-6', 'bench 5: leapfrog(spo)/leapfrog(2 spo), second order -> 4');
    add(`${ex}.bound`, mx(lf) / lf[0], 'rel:1e-6', 'bench 5: bounded leapfrog error');
  }
  add('edriftcmp().q', mx(O.drift(0.3, 500, 50, 'leapfrog')) / mx(O.drift(0.3, 1000, 50, 'leapfrog')), 'rel:1e-6', 'bench 5 defaults e = 0.3, 500 steps/orbit, 50 orbits');
  err('edrift(0.3, 500, 50, 3)', 'method must be 0, 1 or 2');
  err('edrift(0.3, 5000, 100)', 'too many steps for a key press');
  err('edriftcmp(0.3, 2000, 100)', 'four runs: capped at 100000 steps each');
  err('edrift(1, 500, 50)', 'e must be < 1');
}

/* ---------------- bi-elliptic ---------------- */
lines.push('# orbits: Hohmann against bi-elliptic');
{
  for (const R of [2, 12, 15.58, 40]) add(`hohtot(${R})`, O.hohTot(R), 'rel:1e-14', 'Orbits.hohTot');
  for (const [R, Rb] of [[12, 1200], [20, 100], [11.94, 1e12]]) add(`bieltot(${R}, ${Rb})`, O.bielTot(R, Rb), 'rel:1e-14', 'Orbits.bielTot');
  for (const R of [12, 30]) add(`bielinf(${R})`, O.bielInf(R), 'rel:1e-14', 'Orbits.bielInf');
  err('hohtot(-1)', 'ratio must be positive');
  const TH = O.thresholds();
  add('bielthresh().closed', TH.closed, 'rel:1e-12', 'Orbits.thresholds().closed');
  add('bielthresh().finite', TH.finite, 'rel:1e-12', 'Orbits.thresholds().finite');
  add('bielthresh().closed', 11.93877, 'abs:1e-5', 'published threshold R = 11.93877 (orbits bench 6, README)');
  add('hohpeak().R', 15.58, 'abs:0.005', 'published threshold R = 15.58: costliest Hohmann, beyond it every bi-elliptic wins (bench 6)');
  const PK = O.hohPeak();
  add('hohpeak().R', PK.R, 'abs:1e-12', 'Orbits.hohPeak().R');
  add('hohpeak().v', PK.v, 'rel:1e-14', 'Orbits.hohPeak().v');
  for (const R of [12, 14, 15.6, 16]) add(`breakeven(${R})`, O.breakeven(R), 'rel:1e-12', 'Orbits.breakeven');
  add('breakeven(10)', O.breakeven(10), 'abs:0', 'Orbits.breakeven: never cheaper below 11.94 (Infinity)');
  err('breakeven(0.5)', 'ratio must exceed 1');
  for (const [R, ra] of [[12, 15], [5, 5]]) add(`onetangent(${R}, ${ra})`, O.oneTangent(R, ra), 'rel:1e-14', 'Orbits.oneTangent');
  for (const R of [12, 3]) {
    const S = O.oneTangentScan(R);
    add(`onetanscan(${R}).best`, S.best, 'rel:1e-13', 'Orbits.oneTangentScan best');
    add(`onetanscan(${R}).ratio`, S.ratio, 'rel:1e-12', 'Orbits.oneTangentScan ra/r2');
  }
  // Dimensional bi-elliptic: burns from the vis-viva closed form, total cross-checked with Orbits.bielTot.
  for (const [r1, R, q] of [[6678e3, 12, 100], [7000e3, 20, 30]]) {
    const r2 = R * r1, rb = q * r1, a1 = (r1 + rb) / 2, a2 = (r2 + rb) / 2, ex = `bielliptic(${r1}, ${r2}, ${rb})`;
    const d1 = Math.sqrt(MU * (2 / r1 - 1 / a1)) - Math.sqrt(MU / r1), d2 = Math.sqrt(MU * (2 / rb - 1 / a2)) - Math.sqrt(MU * (2 / rb - 1 / a1)), d3 = Math.sqrt(MU * (2 / r2 - 1 / a2)) - Math.sqrt(MU / r2);
    add(`${ex}.dv`, O.bielTot(R, q) * Math.sqrt(MU / r1), 'rel:1e-12', 'Orbits.bielTot(R, rb/r1) * v_circ(r1)');
    add(`${ex}.dv1`, d1, 'rel:1e-12', 'closed form vis-viva burn 1');
    add(`${ex}.dv2`, d2, 'rel:1e-12', 'closed form burn 2 at rb');
    add(`${ex}.dv3`, d3, 'rel:1e-12', 'closed form burn 3 at r2');
    add(`${ex}.tof`, Math.PI * (Math.sqrt(a1 ** 3 / MU) + Math.sqrt(a2 ** 3 / MU)), 'rel:1e-13', 'closed form two half-ellipses');
  }
  err('bielliptic(7000e3, 70000e3, 50000e3)', 'rb below r2');
}

/* ---------------- INS: Allan deviation ---------------- */
lines.push('# ins: Allan deviation model and the seeded simulated records');
{
  const N = 18, B = 5, K = 10 / 60; // bench gyro defaults 0.3 deg/sqrt(h), 5 deg/h, 10 deg/h/sqrt(h) in deg/h with tau in s
  for (const tau of [0.1, 1, 100, 1000]) add(`allandev(${tau}, ${N}, ${B}, ${K})`, I.modelAdev(tau, N, B, K), 'rel:1e-14', 'Ins.modelAdev, gyro bench defaults');
  add('allandev(10, 0, 1, 0)', I.FLAT, 'rel:1e-14', 'Ins.FLAT, flicker floor 0.664 B');
  add('allandev(10, 0, 1, 0)', 0.664, 'abs:5e-4', 'IEEE Std 952 flicker floor 0.664');
  add('allandev(4, 2)', I.modelAdev(4, 2, 0, 0), 'rel:1e-14', 'Ins.modelAdev white noise only');
  err('allandev(0, 1, 1, 1)', 'tau must be positive');
  const tmin = Math.sqrt(3) * N / K;
  add(`allanmin(${N}, ${B}, ${K}).tau`, tmin, 'rel:1e-14', 'closed form d/dtau (N^2/tau + K^2 tau/3) = 0');
  add(`allanmin(${N}, ${B}, ${K}).adev`, I.modelAdev(tmin, N, B, K), 'rel:1e-14', 'Ins.modelAdev at the closed-form minimum');
  add('allanmin(0.01, 0, 0.002).adev', Math.sqrt(2 * 0.01 * 0.002 / Math.sqrt(3)), 'rel:1e-13', 'closed form sqrt(2 N K / sqrt 3)');
  err('allanmin(1, 0, 0)', 'K = 0 has no finite minimum');
  add('flickerquad()', I.flickerQuad(), 'rel:1e-12', 'Ins.flickerQuad (Simpson, X = 400, n = 400000)');
  add('flickerquad()', Math.LN2, 'abs:1e-5', 'int_0^inf sin^4 x / x^3 dx = ln 2');
  add('flickerquad(100, 20000)', I.flickerQuad(100, 20000), 'rel:1e-12', 'Ins.flickerQuad(100, 20000)');
  err('flickerquad(400, 3)', 'n must be even');
  for (const seed of [1, 2]) {
    const S = I.allanSet(seed), F = I.fits(S);
    add(`allanfit(${seed}).wslope`, F.whiteSlope, 'rel:1e-9', `Ins.fits(Ins.allanSet(${seed})) white slope (-1/2)`);
    add(`allanfit(${seed}).nhat`, F.Nhat, 'rel:1e-9', `Ins.fits Nhat (unit white noise -> 1)`);
    add(`allanfit(${seed}).kslope`, F.rrwSlope, 'rel:1e-9', 'Ins.fits rate-random-walk slope (+1/2)');
    add(`allanfit(${seed}).khat`, F.Khat, 'rel:1e-9', 'Ins.fits Khat (-> 1)');
    add(`allanfit(${seed}).flat`, F.flat, 'rel:1e-9', 'Ins.fits flicker floor (-> 0.664)');
    if (seed === 1) {
      const A = I.combinedAdev(S, N, B, K);
      for (const i of [0, 10, 20, 30, I.TAU.length - 1]) {
        add(`allansim(${N}, ${B}, ${K}, 1)(${i + 1}, 1)`, I.TAU[i], 'rel:1e-14', 'Ins.TAU cluster times');
        add(`allansim(${N}, ${B}, ${K}, 1)(${i + 1}, 2)`, A[i], 'rel:1e-9', 'Ins.combinedAdev(Ins.allanSet(1), N, B, K)');
      }
    }
  }
}

/* ---------------- INS: Schuler ---------------- */
lines.push('# ins: Schuler period and loop error');
{
  add('schuler()', I.PERIOD, 'rel:1e-14', 'Ins.PERIOD');
  add('schuler() / 60', 84.41, 'abs:0.005', 'published Schuler period 84.41 min');
  add('schuler(6378137, 9.80665)', TAU * Math.sqrt(6378137 / 9.80665), 'rel:1e-14', 'closed form 2 pi sqrt(R/g)');
  err('schuler(-1)', 'R must be positive');
  const b = 0.1 * I.MG, eps = 0.01 * Math.PI / 180 / 3600;
  for (const [t, phi] of [[10800, 0], [3600, 1e-3], [60, 0]]) {
    add(`schulerx(${t}, ${b}, ${phi}, ${eps})`, I.schulerClosed(t, b, phi, eps), 'rel:1e-12', 'Ins.schulerClosed');
    add(`naivedrift(${t}, ${b}, ${phi}, ${eps})`, I.naiveDrift(t, b, phi, eps), 'rel:1e-13', 'Ins.naiveDrift');
  }
  add('schulerx(10800)', I.schulerClosed(10800, b, 0, eps), 'rel:1e-12', 'bench 2 defaults b = 0.1 mg, eps = 0.01 deg/h, 3 h');
  add(`schulerx(60, ${1e-3 * I.G}, 0, 0) / naivedrift(60, ${1e-3 * I.G}, 0, 0)`, I.schulerClosed(60, 1e-3 * I.G, 0, 0) / I.naiveDrift(60, 1e-3 * I.G, 0, 0), 'rel:1e-9', 'bench 2: short-time agreement ratio');
  for (const [bb, pp, ee, tend] of [[b, 0, eps, 10800], [0, 1e-3, 0, 21600]]) {
    const sim = I.schulerRK4(bb, pp, ee, tend), n = sim.x.length - 1; let dmax = 0;
    for (let i = 0; i <= n; i++) dmax = Math.max(dmax, Math.abs(sim.x[i] - I.schulerClosed(i * sim.h, bb, pp, ee)));
    const ex = `schulersim(${bb}, ${pp}, ${ee}, ${tend})`;
    add(`${ex}.x`, sim.x[n], 'rel:1e-10', 'Ins.schulerRK4 final position error');
    add(`${ex}.dmax`, dmax, 'abs:1e-6', 'Ins.schulerRK4 against Ins.schulerClosed');
    add(`${ex}.period`, I.measuredPeriod(sim), 'rel:1e-10', 'Ins.measuredPeriod (zero crossings)');
  }
  add('schulersim(0, 1e-3, 0, 21600).period', I.PERIOD, 'rel:1e-5', 'measured period matches 2 pi sqrt(R/g)');
  err('schulersim(0, 0, 0, 1e7, 1)', 'more than a million steps');
}

/* ---------------- INS: coning ---------------- */
lines.push('# ins: coning drift and the two-sample correction');
{
  for (const [deg, f, rate] of [[2, 5, 200], [5, 10, 100], [1, 2, 1000]]) {
    const al = deg * Math.PI / 180, W = TAU * f, h = 1 / rate;
    add(`coning0(${al}, ${W}, ${h})`, I.drift0(al, W, h), 'rel:1e-13', 'Ins.drift0');
    add(`coning2(${al}, ${W}, ${h})`, I.drift2(al, W, h), 'rel:1e-13', 'Ins.drift2');
    if (rate <= 200) {
      const run = I.coningRun(al, W, rate), ex = `coningrun(${al}, ${W}, ${rate})`;
      add(`${ex}.slope0`, run.slope0, 'rel:1e-8', 'Ins.coningRun uncorrected drift slope');
      add(`${ex}.slope2`, run.slope2, 'abs:1e-11', 'Ins.coningRun two-sample drift slope');
      add(`${ex}.e0`, run.e0[run.e0.length - 1], 'rel:1e-8', 'Ins.coningRun final error, uncorrected');
      add(`${ex}.slope0`, I.drift0(al, W, h), 'rel:0.05', 'bench 3 check: simulated slope matches drift0');
    }
  }
  add('bortz()', I.bortzCoefficient(), 'rel:1e-9', 'Ins.bortzCoefficient');
  add('bortz()', 2 / 3, 'abs:1e-3', 'two-sample coefficient 2/3 (Miller 1983)');
  add('bortz(0.05, 1000)', I.bortzCoefficient(0.05, 1000), 'rel:1e-9', 'Ins.bortzCoefficient(0.05, 1000)');
  err('coningrun(0.03, 31, 1e6, 10)', 'too many updates');
}

/* ---------------- INS: 1-D INS/GNSS Kalman ---------------- */
lines.push('# ins: 1-D loosely coupled INS/GNSS Kalman filter');
{
  const cases = [[{}, 'insgnss()', 'bench 4 defaults'], [{ sg: 1, L: 120 }, `insgnss(${300e-6 * I.G}, ${5 * I.MG}, ${20e-6 * I.G}, 1, 120)`, 'sg = 1 m, 120 s outage'],
    [{ qa: 100e-6 * I.G, sbw: 0, L: 30, out0: 200 }, `insgnss(${100e-6 * I.G}, ${5 * I.MG}, 0, 3, 30, 200)`, 'no bias walk, outage at 200 s']];
  for (const [opt, ex, why] of cases) {
    const m = I.kfModel(opt), C = I.kfCov(m), kb = Math.round(m.out0 / I.KDT) - 1, ke = Math.round((m.out0 + m.L) / I.KDT) - 1;
    add(`${ex}.sp140`, Math.sqrt(C.post[1400][0]), 'rel:1e-9', `Ins.kfCov position sigma at 140 s, ${why}`);
    add(`${ex}.sbias`, Math.sqrt(C.post[kb][8]), 'rel:1e-9', `Ins.kfCov bias sigma before the outage, ${why}`);
    add(`${ex}.sp0`, Math.sqrt(C.post[kb][0]), 'rel:1e-9', `Ins.kfCov position sigma at outage start, ${why}`);
    add(`${ex}.spend`, Math.sqrt(C.post[ke][0]), 'rel:1e-9', `Ins.kfCov position sigma at outage end, ${why}`);
    const closed = I.outageClosed(C.post[kb], (ke - kb) * I.KDT, m.sa, m.sbw);
    add(`${ex}.closed`, closed, 'rel:1e-9', `Ins.outageClosed, ${why}`);
    if (opt.qa === undefined) {
      add(`insoutage(${mat(C.post[kb])}, ${(ke - kb) * I.KDT}, ${m.qa}, ${m.sbw})`, closed, 'rel:1e-9', `Ins.outageClosed on the filter covariance, ${why}`);
    }
  }
  add('insoutage([100, 0, 0; 0, 1, 0; 0, 0, 0.0024], 60)', I.outageClosed([100, 0, 0, 0, 1, 0, 0, 0, 0.0024], 60, 300e-6 * I.G / Math.sqrt(I.KDT), 20e-6 * I.G), 'rel:1e-12', 'Ins.outageClosed, default noise');
  err('insoutage([1, 2; 3, 4], 60)', 'P must be 3 x 3');
  err('insgnss(0.003, 0.05, 0.0002, 3, 200, 150)', 'outage runs past the 300 s record');
  // Steady state: the module's one-second Phi and Q1, prior-form Riccati iterated to convergence here.
  for (const [opt, ex] of [[{}, 'kfsteady()'], [{ qa: 100e-6 * I.G, sbw: 5e-6 * I.G, sg: 1 }, `kfsteady(${100e-6 * I.G}, ${5e-6 * I.G}, 1)`]]) {
    const m = I.kfModel(opt), Rc = I.riccati1s(m), R = m.sg * m.sg;
    const fapt = (F, P) => I.mm(I.mm(F, P), [F[0], F[3], F[6], F[1], F[4], F[7], F[2], F[5], F[8]]);
    let P = m.P0.slice();
    for (let k = 0; k < 400000; k++) { const S = P[0] + R, h = [P[0], P[1], P[2]]; P = P.map((x, idx) => x - h[Math.floor(idx / 3)] * h[idx % 3] / S); P = fapt(Rc.Phi, P).map((x, i) => x + Rc.Q1[i]); }
    const S = P[0] + R;
    add(`${ex}.sppri`, Math.sqrt(P[0]), 'rel:1e-7', 'prior-form Riccati with Ins.riccati1s Phi, Q1, iterated to steady state');
    add(`${ex}.sp`, Math.sqrt(P[0] - P[0] * P[0] / S), 'rel:1e-7', 'steady posterior position sigma');
    add(`${ex}.sv`, Math.sqrt(P[4] - P[3] * P[3] / S), 'rel:1e-7', 'steady posterior velocity sigma');
    add(`${ex}.sb`, Math.sqrt(P[8] - P[6] * P[6] / S), 'rel:1e-7', 'steady posterior bias sigma');
    add(`${ex}.kp`, P[0] / S, 'rel:1e-7', 'steady Kalman gain, position');
    add(`${ex}.kv`, P[3] / S, 'rel:1e-7', 'steady Kalman gain, velocity');
    add(`${ex}.kb`, P[6] / S, 'rel:1e-7', 'steady Kalman gain, bias');
  }
  err('kfsteady(0.003, 0, 3)', 'sbw = 0 has no positive steady bias sigma');
}

process.stdout.write('# Generated by Tests/reference/gen_orbits.mjs from Where It Isn\'t orbits/orbits.js and ins/ins.js. Do not edit.\n' + lines.join('\n') + '\n');
