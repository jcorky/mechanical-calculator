#!/usr/bin/env python3
"""Reference cases for LibFields.cpp: atlas chapters 7 (lattices, GNSS), 10 (category theory),
11 (procedural systems) and 12 (expansion fields).

These chapters have no JS module: their math lives in the pages' inline bench scripts and in
verification/python/phase{8,11,12,13}_*.py. Expected values are computed here with numpy by the
formulas and algorithms of those scripts (often by a different route than the C++: lstsq instead of
normal equations, Kalman instead of the information recursion, Floyd-Warshall instead of (min, +)
powers, permutations instead of subset DP, dense solves instead of the tridiagonal sweep, Gauss-Legendre
instead of adaptive Simpson), and bench constants (sky plots, lattice basis, IFS presets) are read
from the chapter pages themselves.

    python3 Tests/reference/gen_fields.py /home/user/where-it-isnt > Tests/reference/fields.tsv
"""
import itertools
import json
import math
import os
import re
import sys

import numpy as np

WII = sys.argv[1] if len(sys.argv) > 1 else "/home/user/where-it-isnt"
C = 299_792_458.0
ROWS = []


def page(name):
    with open(os.path.join(WII, "chapters", name), encoding="utf-8") as f:
        return f.read()


CH07, CH10, CH11, CH12 = (page(n) for n in ("ch07-lattices.html", "ch10-category.html", "ch11-procedural.html", "ch12-expansion.html"))


def js_const(src, name):
    """Reads `NAME = [...]` (a plain JSON-like array literal, as in `const A = [..], B = [..]`) from a page script."""
    start = re.search(r"\b" + name + r" = \[", src).end() - 1
    depth = 0
    for i in range(start, len(src)):
        depth += {"[": 1, "]": -1}.get(src[i], 0)
        if depth == 0:
            return json.loads(src[start:i + 1])
    raise ValueError(name)


def num(x):
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return repr(x)


def lit(A):
    A = np.atleast_2d(np.asarray(A, dtype=float))
    return "[" + "; ".join(", ".join(num(v) for v in row) for row in A) + "]"


def vec(v):
    return "[" + ", ".join(num(x) for x in np.ravel(v)) + "]"


def case(expr, expected, tol, prov):
    ROWS.append(f"{expr}\t{expected if isinstance(expected, str) else num(expected)}\t{tol}\t{prov}")


def err(expr, prov):
    case(expr, "error", "abs:0", prov)


def Phi(x):
    return 0.5 * math.erfc(-x / math.sqrt(2))


# =============================================================================== ch07 lattices
P8 = "phase8_lattices.py"
F1, F2 = 1575.42e6, 1227.60e6
case("widelane().l1", C / F1, "rel:1e-12", f"{P8} claim 1: L1 19.03 cm")
case("widelane().l2", C / F2, "rel:1e-12", f"{P8} claim 1: L2 24.42 cm")
case("widelane().wide", C / (F1 - F2), "rel:1e-12", f"{P8} claim 1: widelane 86.19 cm")
case("widelane().narrow", C / (F1 + F2), "rel:1e-12", "narrow lane c/(f1+f2)")
case("widelane(1176.45e6, 1176.45e6)", "error", "abs:0", "equal frequencies have no widelane")

RE, RS = 6_371_000.0, 26_560_000.0


def slant(el):
    return math.sqrt(RS ** 2 - (RE * math.cos(el)) ** 2) - RE * math.sin(el)


case("satrange(pi/2).t", slant(math.pi / 2) / C, "rel:1e-12", f"{P8} claim 1: zenith travel time 67 ms")
case("satrange(0).t", slant(0.0) / C, "rel:1e-12", f"{P8} claim 1: horizon travel time 86 ms")
case("satrange(rad(15)).d", slant(math.radians(15)), "rel:1e-12", "ch07 bench 1 slant range")
err("satrange(0, 7e6, 6e6)", "satellite inside the Earth")
case("tofrange(2*100/c0)", 100.0, "rel:1e-12", f"{P8} claim 1: radar altimeter 100 m, two-way")
case("tofrange(1e-9, 1)", C * 1e-9, "rel:1e-12", f"{P8} claim 1: 1 ns is 0.30 m one-way")
err("tofrange(1, 0)", "way must be positive")
case("tofres(150e6)", C / (2 * 150e6), "rel:1e-12", f"{P8} claim 1: 150 MHz -> 0.99931 m")
err("tofres(0)", "zero bandwidth")

# --- pseudorange positioning (bench 1 sky, receiver on the surface, 2.5 us clock)
SKY, CANYON = js_const(CH07, "SKY"), js_const(CH07, "CANYON")


def sat_positions(azel):
    az = np.radians([a for a, _ in azel]); el = np.radians([e for _, e in azel])
    d = np.array([slant(e) for e in el])
    u = np.stack([np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)], 1)
    return np.array([0, 0, RE]) + u * d[:, None]


def gn_solve(S, rho, iters=15):
    x = np.zeros(4)
    for _ in range(iters):
        d = np.linalg.norm(S - x[:3], axis=1)
        G = np.hstack([-(S - x[:3]) / d[:, None], np.ones((len(S), 1))])
        dx = np.linalg.lstsq(G, rho - (d + x[3]), rcond=None)[0]
        x = x + dx
        if np.linalg.norm(dx) < 1e-6:
            break
    return x, G


def dops(G, up=np.array([0, 0, 1.0])):
    Q = np.linalg.inv(G.T @ G); P2 = np.trace(Q[:3, :3]); V2 = up @ Q[:3, :3] @ up
    return dict(gdop=math.sqrt(np.trace(Q)), pdop=math.sqrt(P2), hdop=math.sqrt(P2 - V2), vdop=math.sqrt(V2), tdop=math.sqrt(Q[3, 3]))


truth = np.array([0, 0, RE]); cb = C * 2.5e-6
S8 = sat_positions(SKY[:8]); rho8 = np.linalg.norm(S8 - truth, axis=1) + cb
x8, G8 = gn_solve(S8, rho8)
pre = f"S = {lit(S8)}; pr = {vec(rho8)}; "
for k, f in enumerate(["x", "y", "z", "cb"]):
    case(pre + f"gnssfix(S, pr).{f}", x8[k], "abs:1e-5", f"{P8} claim 2: Gauss-Newton from Earth's centre recovers truth (bench default, 8 sats)")
case(pre + "gnssfix(S, pr).dt", 2.5e-6, "abs:1e-13", "ch07 bench 1: clock offset 2.5 us")
gen = np.random.default_rng(7)
noise = gen.normal(0, 3.0, 8)
xn, _ = gn_solve(S8, rho8 + noise)
pre_n = f"S = {lit(S8)}; pr = {vec(rho8 + noise)}; "
case(pre_n + "gnssfix(S, pr).x", xn[0], "abs:1e-5", "noisy ranges (sigma 3 m): least-squares fix by lstsq")
case(pre_n + "gnssfix(S, pr).z", xn[2], "abs:1e-5", "noisy ranges (sigma 3 m): least-squares fix by lstsq")
S4 = sat_positions(SKY[:4]); rho4 = np.linalg.norm(S4 - truth, axis=1) + cb
x4, _ = gn_solve(S4, rho4)
case(f"gnssfix({lit(S4)}, {vec(rho4)}).cb", x4[3], "abs:1e-5", "4 satellites: exactly determined")
err(f"gnssfix({lit(S4[:3])}, {vec(rho4[:3])})", f"{P8} claim 2: 3 satellites cannot separate clock and position")
err(f"gnssfix({lit(S4)}, [1, 2, 3])", "pseudorange count must match")

D8 = dops(G8)
az8 = "rad(" + vec([a for a, _ in SKY[:8]]) + ")"; el8 = "rad(" + vec([e for _, e in SKY[:8]]) + ")"
for f in ["gdop", "pdop", "hdop", "vdop", "tdop"]:
    case(f"dopae({az8}, {el8}).{f}", D8[f], "rel:1e-9", f"ch07 bench 1 open sky (PDOP 1.92 quoted)")
Sc = sat_positions(CANYON[:8]); Gc = gn_solve(Sc, np.linalg.norm(Sc - truth, axis=1) + cb)[1]
azc = "rad(" + vec([a for a, _ in CANYON[:8]]) + ")"; elc = "rad(" + vec([e for _, e in CANYON[:8]]) + ")"
case(f"dopae({azc}, {elc}).pdop", dops(Gc)["pdop"], "rel:1e-9", "ch07 bench 1 urban canyon (PDOP 6.86 quoted)")
case(f"dop({lit(S8)}, [0, 0, {num(RE)}]).hdop", D8["hdop"], "rel:1e-9", "ch07 bench 1: dop from positions = dop from angles")
# The same sky rotated rigidly: DOPs are invariant, vertical follows the receiver.
t1, t2 = 0.7, -1.1
Rz = np.array([[math.cos(t1), -math.sin(t1), 0], [math.sin(t1), math.cos(t1), 0], [0, 0, 1]])
Rx = np.array([[1, 0, 0], [0, math.cos(t2), -math.sin(t2)], [0, math.sin(t2), math.cos(t2)]])
Rm = Rz @ Rx
case(f"dop({lit(S8 @ Rm.T)}, {vec(Rm @ truth)}).vdop", D8["vdop"], "rel:1e-8", "rotation invariance of VDOP (geometry rotated)")
case(f"dop({lit(S8 @ Rm.T)}, {vec(Rm @ truth)}).hdop", D8["hdop"], "rel:1e-8", "rotation invariance of HDOP (geometry rotated)")
err("dopae([0, 1, 2], [0.5, 0.5, 0.5])", "fewer than 4 satellites")

# --- one biased range and Baarda's normalized residual (claim 3: 30 m reflection on the 15-degree satellite)
e = np.zeros(8); e[3] = 30.0
xb, Gb = gn_solve(S8, rho8 + e)
res = (rho8 + e) - (np.linalg.norm(S8 - xb[:3], axis=1) + xb[3])
Gx = np.hstack([-(S8 - xb[:3]) / np.linalg.norm(S8 - xb[:3], axis=1)[:, None], np.ones((8, 1))])
redund = np.diag(np.eye(8) - Gx @ np.linalg.solve(Gx.T @ Gx, Gx.T))
w = res / (3.0 * np.sqrt(redund))
assert np.argmax(np.abs(w)) == 3 and np.argmax(np.abs(res)) != 3
pre_b = f"S = {lit(S8)}; pr = {vec(rho8 + e)}; "
case(pre_b + "baarda(S, pr, 3)(4)", w[3], "rel:1e-7", f"{P8} claim 3: normalized residual identifies the reflected satellite (4.7 quoted)")
case(pre_b + "baarda(S, pr, 3)(1)", w[0], "rel:1e-7", f"{P8} claim 3: w-statistic of satellite 1")
case(pre_b + "max(abs(baarda(S, pr, 3)))", np.abs(w).max(), "rel:1e-7", f"{P8} claim 3: largest |w| is the reflected one")
err(pre_b + "baarda(S, pr, 0)", "sigma must be positive")

# --- Lagrange-Gauss reduction (bench 3: G0 sheared k = 5 times)
G0 = np.array(js_const(CH07, "G0")); SHEARS = js_const(CH07, "SHEARS")


def sheared(k):
    b = G0.copy()
    for p, q, m in SHEARS[:k]:
        b[p] = b[p] + m * b[q]
    return b


def gauss(b1, b2):
    b1, b2 = np.array(b1, float), np.array(b2, float)
    while True:
        if b2 @ b2 < b1 @ b1:
            b1, b2 = b2, b1
        q = round(b1 @ b2 / (b1 @ b1))
        if q == 0:
            return b1, b2
        b2 = b2 - q * b1


def shortest_enum(B, R=12):
    B = np.asarray(B, float); n = len(B); best = math.inf
    for c in itertools.product(range(-R, R + 1), repeat=n):
        if any(c):
            best = min(best, np.linalg.norm(np.array(c) @ B))
    return best


B5 = sheared(5); r1, r2 = gauss(*B5)
lam1 = shortest_enum(G0)
assert abs(np.linalg.norm(r1) - lam1) < 1e-12
case(f"G = gaussred({lit(B5)}); hypot(G(1,1), G(1,2))", lam1, "rel:1e-12", "ch07 bench 3: |b1| equals the enumerated shortest vector")
case(f"G = gaussred({lit(B5)}); hypot(G(2,1), G(2,2))", np.linalg.norm(r2), "rel:1e-12", "ch07 bench 3: reduced b2 (Lagrange-Gauss, phase8 claim 5)")
case(f"abs(det(gaussred({lit(B5)})))", abs(np.linalg.det(G0)), "rel:1e-12", "ch07 bench 3: the area never changes")
case(f"G = gaussred({lit(sheared(8))}); hypot(G(1,1), G(1,2))", lam1, "rel:1e-12", "ch07 bench 3: all 8 shears, same lattice")
err("gaussred([1, 2; 2, 4])", "parallel basis")


# --- LLL and the shortest vector (claim 6)
def gram_schmidt(B):
    Bs = np.zeros_like(B, dtype=float); mu = np.zeros((len(B), len(B)))
    for i in range(len(B)):
        Bs[i] = B[i]
        for j in range(i):
            mu[i, j] = B[i] @ Bs[j] / (Bs[j] @ Bs[j]); Bs[i] -= mu[i, j] * Bs[j]
    return Bs, mu


def lll(B, delta=0.75):
    B = B.astype(float).copy(); k = 1
    while k < len(B):
        for j in range(k - 1, -1, -1):
            _, mu = gram_schmidt(B)
            q = round(mu[k, j])
            if q:
                B[k] -= q * B[j]
        Bs, mu = gram_schmidt(B)
        if Bs[k] @ Bs[k] >= (delta - mu[k, k - 1] ** 2) * (Bs[k - 1] @ Bs[k - 1]):
            k += 1
        else:
            B[[k, k - 1]] = B[[k - 1, k]]; k = max(k - 1, 1)
    return B


gen = np.random.default_rng(8)
B4 = gen.integers(-40, 41, size=(4, 4)).astype(float)
L4 = lll(B4)
for (i, j) in [(1, 1), (1, 4), (3, 2)]:
    case(f"L = lll({lit(B4)}); L({i},{j})", L4[i - 1, j - 1], "abs:1e-9", f"{P8} claim 6: LLL (delta 3/4) of a random integer basis")
case(f"L = lll({lit(B4)}, 0.99); L(1,1)^2 + L(1,2)^2 + L(1,3)^2 + L(1,4)^2", float(lll(B4, 0.99)[0] @ lll(B4, 0.99)[0]), "abs:1e-9", "LLL with delta 0.99")
err("lll([1, 2; 2, 4])", "dependent rows")

# exact shortest vector by the script's Fincke-Pohst enumeration


def shortest_fp(B):
    Bs, mu = gram_schmidt(B); n = len(B); best = [min(np.linalg.norm(b) for b in B) ** 2 + 1e-9]
    norms = np.array([b @ b for b in Bs])

    def rec(i, coef, partial):
        if i < 0:
            if 1e-12 < partial < best[0]:
                best[0] = partial
            return
        c = -sum(coef[j] * mu[j, i] for j in range(i + 1, n))
        r = math.sqrt(max(best[0] - partial, 0) / norms[i])
        for zi in range(math.ceil(c - r), math.floor(c + r) + 1):
            coef[i] = zi; p = partial + (zi - c) ** 2 * norms[i]
            if p < best[0]:
                rec(i - 1, coef, p)
        coef[i] = 0

    rec(n - 1, [0] * n, 0.0)
    return math.sqrt(best[0])


lam4 = shortest_fp(L4)
case(f"shortvec({lit(B4)})", lam4, "rel:1e-9", f"{P8} claim 6: lambda_1 by Fincke-Pohst enumeration")
case(f"shortvec({lit(G0)})", lam1, "rel:1e-12", "ch07 bench 3: lambda_1 of G0")
B6 = gen.integers(-40, 41, size=(6, 6)).astype(float)
case(f"shortvec({lit(B6)})", shortest_fp(lll(B6)), "rel:1e-9", f"{P8} claim 6: 6-D random lattice")
case("lllbound(6)", 2 ** 2.5, "rel:1e-12", f"{P8} claim 6: 2^((n-1)/2) for n = 6")
case("lllbound(2, 0.99)", (1 / 0.74) ** 0.5, "rel:1e-12", "LLL bound for delta 0.99")
err("lllbound(3, 0.2)", "delta must exceed 1/4")
case(f"odefect({lit(B5)})", np.prod(np.linalg.norm(B5, axis=1)) / abs(np.linalg.det(B5)), "rel:1e-12", "ch07 bench 3: orthogonality defect of the sheared basis")
case(f"odefect(gaussred({lit(B5)}))", np.linalg.norm(r1) * np.linalg.norm(r2) / abs(np.linalg.det(G0)), "rel:1e-12", "ch07 bench 3: defect after reduction")
err("odefect([1, 1; 2, 2])", "dependent rows")


# --- integer ambiguities (claim 7, bench 2)
def tilted(theta_deg, a, b):
    t = math.radians(theta_deg); R = np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])
    return R @ np.diag([a ** 2, b ** 2]) @ R.T


def decorrelate(Qm):
    Z = np.eye(2, dtype=int); Qz = Qm.copy()
    while True:
        if Qz[1, 1] < Qz[0, 0]:
            P = np.array([[0, 1], [1, 0]]); Z = Z @ P; Qz = P.T @ Qz @ P
        mu = round(Qz[0, 1] / Qz[0, 0])
        if mu == 0:
            return Z, Qz
        T = np.array([[1, -mu], [0, 1]]); Z = Z @ T; Qz = T.T @ Qz @ T


def p_boot(Qm):
    d = np.diag(np.linalg.cholesky(Qm)) ** 2
    return float(np.prod([2 * Phi(1 / (2 * math.sqrt(di))) - 1 for di in d]))


def p_round2(Qm, panels=4000):
    """P(|x1| < 1/2, |x2| < 1/2) for x ~ N(0, Q), composite Gauss-Legendre in x1 with the band edges as breakpoints."""
    s1 = math.sqrt(Qm[0, 0]); l = Qm[1, 0] / Qm[0, 0]; s2 = math.sqrt(Qm[1, 1] - Qm[1, 0] ** 2 / Qm[0, 0])
    cuts = {-0.5, 0.5}
    if l != 0:
        cuts |= {c for c in (0.5 / l, -0.5 / l) if abs(c) < 0.5}
    cuts = sorted(cuts); xg, wg = np.polynomial.legendre.leggauss(20); total = 0.0
    for a, b in zip(cuts[:-1], cuts[1:]):
        edges = np.linspace(a, b, panels // (len(cuts) - 1) + 1)
        for lo, hi in zip(edges[:-1], edges[1:]):
            x = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
            dens = np.exp(-0.5 * (x / s1) ** 2) / (s1 * math.sqrt(2 * math.pi))
            band = np.array([Phi((0.5 - l * xi) / s2) - Phi((-0.5 - l * xi) / s2) for xi in x])
            total += 0.5 * (hi - lo) * np.sum(wg * dens * band)
    return total


Qt = tilted(52, 1.2, 0.07)
ell = "ellcov(1.2, 0.07, rad(52))"
case(f"{ell}(1,1)", Qt[0, 0], "rel:1e-12", f"{P8} claim 7: thin ellipse at 52 deg (bench 2 default)")
case(f"{ell}(1,2)", Qt[0, 1], "rel:1e-12", f"{P8} claim 7: thin ellipse at 52 deg (bench 2 default)")
case("ellcov()(2,2)", Qt[1, 1], "rel:1e-12", "bench 2 defaults a 1.2, b 0.07, 52 deg")
err("ellcov(0, 1, 0)", "semi-axis must be positive")
Zt, Qzt = decorrelate(Qt)
case(f"Z = zdecor({ell}); Z(1,1)", Zt[0, 0], "abs:0", "ch07 bench 2: decorrelating Z")
case(f"Z = zdecor({ell}); Z(1,2) + 10*Z(2,1) + 100*Z(2,2)", Zt[0, 1] + 10 * Zt[1, 0] + 100 * Zt[1, 1], "abs:0", "ch07 bench 2: decorrelating Z")
case(f"abs(det(zdecor({ell})))", 1.0, "abs:1e-12", f"{P8} claim 7: Z is unimodular")
err("zdecor([1, 2; 2, 1])", "not positive definite")
case(f"pboot({ell})", p_boot(Qt), "rel:1e-12", f"{P8} claim 7 / bench 2: bootstrapping 0.50 quoted")
case(f"Q = {ell}; Z = zdecor(Q); pboot(Z'*Q*Z)", p_boot(Qzt), "rel:1e-10", f"{P8} claim 7: decorrelation lifts bootstrapping to 0.80 (= ILS)")
case("pboot(0.15^2*eye(2))", (2 * Phi(1 / 0.3) - 1) ** 2, "rel:1e-12", f"{P8} claim 7 case 'sigma 0.15, rho 0'")
Q3 = np.array([[0.09, 0.05, 0.01], [0.05, 0.08, 0.02], [0.01, 0.02, 0.07]])
case(f"pboot({lit(Q3)})", p_boot(Q3), "rel:1e-12", "bootstrapping success, 3 ambiguities")
err("pboot([1, 2; 2, 1])", "not positive definite")
Q09 = 0.25 ** 2 * np.array([[1, 0.9], [0.9, 1]])
case(f"pround({ell})", p_round2(Qt), "rel:1e-8", f"{P8} claim 7 / ch07 failure list: rounding 0.40 on the thin ellipse")
case("pround(0.25^2*[1, 0.9; 0.9, 1])", p_round2(Q09), "rel:1e-8", f"{P8} claim 7 case 'sigma 0.25, rho 0.9'")
case("pround(0.15^2*eye(2))", (2 * Phi(1 / 0.3) - 1) ** 2, "rel:1e-10", "diagonal Q: rounding = bootstrapping")
case("pround(0.3^2)", 2 * Phi(1 / 0.6) - 1, "rel:1e-12", "scalar ambiguity 2 Phi(1/(2 sigma)) - 1")
err(f"pround({lit(Q3)})", "no closed form for n = 3")
# Monte Carlo sanity of the 2-D integral (not written; guards the generator)
_s = np.random.default_rng(3).multivariate_normal([0, 0], Qt, 400000)
assert abs(np.mean(np.all(np.round(_s) == 0, 1)) - p_round2(Qt)) < 0.004

# --- Weyl and the circle (claim 9)
alpha = (math.sqrt(5) - 1) / 2


def weyl(a, n, bins=20):
    pts = (a * np.arange(1, n + 1)) % 1
    xs = np.sort(pts); i = np.arange(1, n + 1)
    dstar = max(np.max(i / n - xs), np.max(xs - (i - 1) / n))
    hist = np.bincount(np.minimum((pts * bins).astype(int), bins - 1), minlength=bins)
    return dstar, np.abs(hist / n - 1 / bins).max()


d_s, d_v = weyl(alpha, 100000)
assert d_v < 1e-3
case("weyl((sqrt(5) - 1)/2).dev", d_v, "rel:1e-9", f"{P8} claim 9: every bin within 0.1% for the golden ratio")
case("weyl((sqrt(5) - 1)/2).dstar", d_s, "rel:1e-9", "star discrepancy of the golden Kronecker sequence")
case("weyl(sqrt(2), 1000, 10).dstar", weyl(math.sqrt(2), 1000, 10)[0], "rel:1e-9", "star discrepancy, alpha = sqrt 2")
case("weyl(0.5, 1000).dev", weyl(0.5, 1000)[1], "rel:1e-12", "rational alpha is not equidistributed")
err("weyl(0.3, 0)", "n must be >= 1")
ph = np.cumsum(np.random.default_rng(9).normal(0.31, 0.05, 60)) * 2 * math.pi
wrapped = (ph + math.pi) % (2 * math.pi) - math.pi
un = np.unwrap(wrapped)
case(f"unwrap({vec(wrapped)})(60)", un[-1], "rel:1e-12", f"{P8} claim 9: unwrapping recovers phase up to a constant number of turns (np.unwrap)")
case(f"u = unwrap({vec(wrapped)}); (u(60) - u(1)) - ({num(ph[-1] - ph[0])})", 0.0, "abs:1e-9", f"{P8} claim 9: the unwrapped increment equals the true one")
err("unwrap([1, 2; 3, 4])", "phases must be a vector")

# =============================================================================== ch10 category
P11 = "phase11_category.py"
gen = np.random.default_rng(11)


def stoch(m, n):
    K = gen.random((m, n)); return K / K.sum(0, keepdims=True)


f43, g54 = stoch(4, 3), stoch(5, 4)
case(f"kcomp({lit(g54)}, {lit(f43)})(2,3)", (g54 @ f43)[1, 2], "rel:1e-12", f"{P11} claim 1: composite of kernels")
case(f"sum(kcomp({lit(g54)}, {lit(f43)}))", 3.0, "rel:1e-12", f"{P11} claim 1: composites are kernels (each column sums to 1)")
err(f"kcomp({lit(f43)}, {lit(g54)})", "shapes do not chain")
err("kcomp([0.5, 0.5; 0.6, 0.5], eye(2))", "column does not sum to 1")
P6 = stoch(6, 6)
case(f"kpow({lit(P6)}, 5)(4,2)", np.linalg.matrix_power(P6, 5)[3, 1], "rel:1e-12", f"{P11} claim 1: Chapman-Kolmogorov P^5")
case(f"P = {lit(P6)}; kpow(P, 5)(1,6) - kcomp(kpow(P, 2), kpow(P, 3))(1,6)", 0.0, "abs:1e-15", f"{P11} claim 1: P^5 = P^2 P^3")
case(f"trace(kpow({lit(P6)}, 0))", 6.0, "abs:0", "P^0 is the identity")
err("kpow([1, 2; 3, 4], 2)", "not a kernel")


def dagger(K, p):
    joint = K * p[None, :]
    return (joint / joint.sum(1, keepdims=True)).T


pr = gen.random(3); pr /= pr.sum()
lhs = dagger(g54 @ f43, pr); rhs = dagger(f43, pr) @ dagger(g54, f43 @ pr)
case(f"bayesinv({lit(f43)}, {vec(pr)})(3,2)", dagger(f43, pr)[2, 1], "rel:1e-12", f"{P11} claim 2: Bayesian inverse P(x | y)")
case(f"bayesinv(kcomp({lit(g54)}, {lit(f43)}), {vec(pr)})(2,5)", rhs[1, 4], "rel:1e-11",
     f"{P11} claim 2: inverse of the composite = composite of inverses (value from the right side)")
case(f"f = {lit(f43)}; p = {vec(pr)}; bayesinv(bayesinv(f, p), f*p')(4,1)", f43[3, 0], "rel:1e-11", f"{P11} claim 2: inverting twice returns f")
err(f"bayesinv({lit(f43)}, [0.5, 0.5])", "prior length must match")
err(f"bayesinv({lit(f43)}, [0.5, 0.5, 0.5])", "prior must sum to 1")


def ring(n, s):
    f = np.zeros((n, n))
    for x in range(n):
        f[(x + 1) % n, x] += 1 - 2 * s; f[(x + 2) % n, x] += s; f[x, x] += s
    return f


case("ringkern()(3,1)", ring(8, 0.2)[2, 0], "rel:1e-12", "ch10 bench 2 motion kernel (slip 0.2)")
case("kpow(ringkern(8, 0.2), 2)(3,1)", (ring(8, 0.2) @ ring(8, 0.2))[2, 0], "rel:1e-12", "ch10 bench 2: f^k for k = 2")
err("ringkern(8, 0.7)", "slip above 1/2")
# bench 2 posterior: doors at cells 1 and 5, reliability 0.85, k = 2, prior 1 + 0.6 cos i
NK, DOORS, rr = 8, [1, 5], 0.85
gk = np.array([[rr if x in DOORS else 1 - rr for x in range(NK)], [1 - rr if x in DOORS else rr for x in range(NK)]])
pk = np.array([(1 / NK) * (1 + 0.6 * math.cos(i)) for i in range(NK)]); pk /= pk.sum()
fk = ring(NK, 0.2) @ ring(NK, 0.2)
post = dagger(gk @ fk, pk)
case(f"bayesinv(kcomp({lit(gk)}, kpow(ringkern(), 2)), {vec(pk)})(4,1)", post[3, 0], "rel:1e-11", "ch10 bench 2: posterior of start cell 3 given 'door'")

# transfer functions: bench 1 blocks p = 2, k = 1.5, w = 5, zeta = 0.3
p_, k_, w_, z_ = 2.0, 1.5, 5.0, 0.3
g1 = lambda s: p_ / (s + p_)
g2 = lambda s: k_ * w_ ** 2 / (s * s + 2 * z_ * w_ * s + w_ ** 2)
G1s, G2s = "tf(2, [1, 2])", f"tf({num(k_ * w_ ** 2)}, [1, {num(2 * z_ * w_)}, {num(w_ ** 2)}])"
for om in (1.0, 5.0):
    s = 1j * om; ser = g2(s) * g1(s)
    case(f"tfjw(tfser({G1s}, {G2s}), {num(om)}).mag", abs(ser), "rel:1e-12", f"ch10 bench 1 series = product (w = {om})")
    case(f"tfjw(tfser({G1s}, {G2s}), {num(om)}).phase", np.angle(ser), "rel:1e-12", f"ch10 bench 1 series phase (w = {om})")
s = 2j
case(f"tfjw(tfpar({G1s}, {G2s}), 2).re", (g1(s) + g2(s)).real, "rel:1e-12", "ch10 bench 1 parallel = sum")
case(f"tfjw(tffb(tfser({G1s}, {G2s})), 2).im", (g1(s) * g2(s) / (1 + g1(s) * g2(s))).imag, "rel:1e-12", "ch10 bench 1 unity feedback G/(1+G)")
case(f"tfjw(tffb({G2s}, {G1s}), 3).db", 20 * math.log10(abs(g2(3j) / (1 + g2(3j) * g1(3j)))), "rel:1e-12", "feedback through H: G/(1 + G H)")
case(f"tfser({G1s}, {G2s})(2,2)", float(np.polymul([1, 2], [1, 3, 25])[1]), "abs:0", "series denominator (s + 2)(s^2 + 3 s + 25)")
case("tf([0, 0, 1], [1, 2])(1,2)", 1.0, "abs:0", "leading zeros trimmed, rows padded")
case("tfjw(10, 3).db", 20.0, "rel:1e-12", "a scalar is a static gain")
err("tf(1, [0, 0])", "zero denominator")
err("tfjw(tf(1, [1, 0]), 0)", "s = 0 is a pole")


def tf_ss(sys_, s):
    A, B, Cm, D = sys_
    return Cm @ np.linalg.solve(s * np.eye(len(A)) - A, B) + D


Gs1 = (np.array([[-p_]]), np.array([[p_]]), np.array([[1.0]]), np.array([[0.0]]))
Gs2 = (np.array([[0, 1], [-w_ ** 2, -2 * z_ * w_]]), np.array([[0], [k_ * w_ ** 2]]), np.array([[1.0, 0]]), np.array([[0.0]]))
A1, B1, C1, D1 = Gs1; A2, B2, C2, D2 = Gs2
As = np.block([[A1, np.zeros((1, 2))], [B2 @ C1, A2]]); Bs_ = np.vstack([B1, B2 @ D1]); Cs = np.hstack([D2 @ C1, C2]); Ds = D2 @ D1
case(f"ssjw({lit(As)}, {lit(Bs_)}, {lit(Cs)}, {lit(Ds)}, 1).re", (g2(1j) * g1(1j)).real, "rel:1e-12", f"{P11} claim 3: series realization's transfer function is the product")
case(f"ssjw({lit(As.T)}, {lit(Cs.T)}, {lit(Bs_.T)}, {lit(Ds)}, 4).im", tf_ss((As.T, Cs.T, Bs_.T, Ds), 4j)[0, 0].imag, "rel:1e-12", f"{P11} claim 4: the dual realization has transfer G'")
case(f"ssjw([-2], [2], [1], 0.5, 2).mag", abs(g1(2j) + 0.5), "rel:1e-12", "first-order block with feedthrough")
err("ssjw([0, 1; 0, 0], [0; 1], [1, 0], 0, 0)", "jw = 0 is an eigenvalue of A")
err("ssjw([1, 2; 3, 4], [1; 1], [1, 1, 1], 0, 1)", "C has the wrong width")


# sheaves: bench 3 frames from mulberry32(17)
def mulberry32(seed):
    s = seed & 0xFFFFFFFF

    def nxt():
        nonlocal s
        s = (s + 0x6D2B79F5) & 0xFFFFFFFF; t = s
        t = ((t ^ (t >> 15)) * (t | 1)) & 0xFFFFFFFF
        t ^= (t + (((t ^ (t >> 7)) * (t | 61)) & 0xFFFFFFFF)) & 0xFFFFFFFF
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296
    return nxt


r17 = mulberry32(17); FR = [(r17() * 2 - 1) * math.pi for _ in range(6)]


def rot(t):
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


def sheaf_lap(frames, extra):
    n = len(frames); delta = np.zeros((2 * n, 2 * n))
    for i in range(n):
        j = (i + 1) % n; Fi, Fj = rot(frames[i]), rot(frames[j])
        if j == 0:
            Fj = rot(extra) @ Fj
        delta[2 * i:2 * i + 2, 2 * i:2 * i + 2] += Fi; delta[2 * i:2 * i + 2, 2 * j:2 * j + 2] -= Fj
    return delta.T @ delta


frs = vec(FR)
for extra, expect in [(0.0, 2), (0.3, 0), (math.pi, 0)]:
    evs = np.linalg.eigvalsh(sheaf_lap(FR, extra)); assert int(np.sum(evs < 1e-9)) == expect
    case(f"sheafh0({frs}, {num(extra)}).h0", expect, "abs:0", f"{P11} claim 6: H0 dimension with loop error {extra:.2f} rad")
    case(f"sheafh0({frs}, {num(extra)}).lmin", evs[evs > 1e-9].min(), "rel:1e-8", "ch10 bench 3: slowest diffusion rate (eigvalsh)")
case(f"sheaflap({frs}, 0.3)(3,4)", sheaf_lap(FR, 0.3)[2, 3], "abs:1e-12", "ch10 bench 3 sheaf Laplacian entry")
Lg = np.zeros((6, 6))
for i in range(6):
    j = (i + 1) % 6; Lg[i, i] += 1; Lg[j, j] += 1; Lg[i, j] -= 1; Lg[j, i] -= 1
Lk = np.kron(Lg, np.eye(2))
case("sheaflap(zeros(6, 1))(1,3)", Lk[0, 2], "abs:1e-15", f"{P11} claim 6: the constant sheaf is L (x) I (consensus)")
case("sheaflap(zeros(6, 1))(5,5)", Lk[4, 4], "abs:1e-15", f"{P11} claim 6: the constant sheaf is L (x) I (consensus)")
err("sheafh0([1, 2; 3, 4])", "frames must be a vector")

# =============================================================================== ch11 procedural
P12 = "phase12_procedural.py"


def run_ca(rule, n, steps):
    row = np.zeros(n, int); row[n // 2] = 1; out = [row]; tbl = np.array([(rule >> k) & 1 for k in range(8)])
    for _ in range(steps):
        row = tbl[4 * np.roll(row, 1) + 2 * row + np.roll(row, -1)]; out.append(row)
    return np.array(out)


S90 = run_ca(90, 201, 63)
for t, j in [(63, 1), (40, 8), (10, 2)]:
    want = math.comb(t, (t + j) // 2) % 2 if (t + j) % 2 == 0 else 0
    assert S90[t, 100 + j] == want
    case(f"ca(90, {t})({101 + j})", want, "abs:0", f"{P12} claim 2: Rule 90 is Pascal's triangle mod 2 (row {t}, offset {j})")
case("sum(ca(90, 63))", S90[63].sum(), "abs:0", "Rule 90 row 63: 2^popcount(63) ones")
case("sum(ca(30, 100))", run_ca(30, 201, 100)[100].sum(), "abs:0", "ch11 bench 2: Rule 30, 201 cells, 100 steps")
case("sum(ca(110, 150, 101))", run_ca(110, 101, 150)[150].sum(), "abs:0", "Rule 110 wraps on a 101-cell ring")
err("ca(256, 10)", "rule above 255")
C30 = run_ca(30, 4001, 1999)[:, 2000]
case("cafreq(30)", C30[1:].mean(), "abs:1e-15", f"{P12} claim 2 / ch11 bench 2: Rule 30 centre column 1-frequency")
case("cafreq(90, 64)", run_ca(90, 129, 64)[1:, 64].mean(), "abs:1e-15", "Rule 90 centre column (zero on odd steps)")
case("sum(cacenter(30, 200))", C30[:201].sum(), "abs:0", "Rule 30 centre column bits t = 0..200")
case("cacenter(30, 50)(48)", C30[47], "abs:0", "Rule 30 centre bit at t = 47")
err("cacenter(30, 5000)", "too many steps for a keypress")

case("moran([0.5, 0.5, 0.5])", math.log(3) / math.log(2), "rel:1e-13", f"{P12} claim 1: Sierpinski dimension log 3 / log 2")
case("moran([0.5, 0.25, 0.25])", 1.0, "rel:1e-13", "0.5^d + 2 * 0.25^d = 1 at d = 1")
case("moran(ones(1, 4)/3)", math.log(4) / math.log(3), "rel:1e-13", "four maps of ratio 1/3 (Koch curve)")
err("moran([0.5, 1.2])", "an expanding ratio")
PRE = {k: np.array([[float(v) for v in ln.split()] for ln in s.split("\\n")]) for k, s in re.findall(r'(\w+): "([-0-9. \\n]+)"', CH11)}
case("ifspreset(1)(2,1)", PRE["fern"][1, 0], "abs:0", "ch11 bench 1 fern preset")
case("ifspreset(2)(3,6)", PRE["sierpinski"][2, 5], "abs:0", "ch11 bench 1 Sierpinski preset")
err("ifspreset(6)", "no such preset")
norms = [np.linalg.norm(np.array([[r[0], r[1]], [r[2], r[3]]]), 2) for r in PRE["fern"]]
assert max(norms) < 1
for k in range(4):
    case(f"ifsnorm(ifspreset(1))({k + 1})", norms[k], "rel:1e-12", f"{P12} claim 1: Barnsley fern map {k + 1} contracts (2-norm)")
case("max(ifsnorm(ifspreset(2)))", 0.5, "rel:1e-15", f"{P12} claim 1: Sierpinski maps contract by 1/2")
case("ifsnorm(ifspreset(3))(2)", np.linalg.norm(np.array([[0.42, -0.42], [0.42, 0.42]]), 2), "rel:1e-12", "ch11 bench 1 tree preset")
err("ifsnorm([1, 2, 3])", "maps need 6 or 7 columns")


def mean_log(M):
    ps = np.maximum(0, M[:, 6]); ps = ps / ps.sum()
    return float(sum(p * math.log(np.linalg.norm(np.array([[r[0], r[1]], [r[2], r[3]]]), 2)) for p, r in zip(ps, M)))


case("ifslog(ifspreset(4))", mean_log(PRE["expand"]), "rel:1e-12", f"{P12} claim 1: one map of 1.05, still contracting on average")
case("ifslog(ifspreset(5))", mean_log(PRE["avgexpand"]), "rel:1e-12", f"{P12} claim 1: mean log factor > 0 escapes")
case("ifslog([0.5, 0, 0, 0.5, 0, 0; 0.25, 0, 0, 0.25, 1, 0])", 0.5 * (math.log(0.5) + math.log(0.25)), "rel:1e-12", "equal weights without a p column")


def turing_ref(Dv, Du=1.0, a=0.1, b=0.9):
    us, vs = a + b, b / (a + b) ** 2
    J = np.array([[-1 + 2 * us * vs, us ** 2], [-2 * us * vs, -us ** 2]])
    # closed form of the band (det < 0) and of the stationary point d lambda / d q = 0
    a2, a1, a0 = Du * Dv, -(J[0, 0] * Dv + J[1, 1] * Du), np.linalg.det(J)
    disc = a1 * a1 - 4 * a2 * a0
    if a1 >= 0 or disc <= 0:
        return None
    q1, q2 = (-a1 - math.sqrt(disc)) / (2 * a2), (-a1 + math.sqrt(disc)) / (2 * a2)
    best = None
    for sgn in (1, -1):
        A = sgn * math.sqrt(-Du * J[0, 1] * J[1, 0] / Dv)
        q = (J[0, 0] - J[1, 1] - A * (1 + Dv / Du)) / (Du - Dv)
        lam = J[0, 0] - Du * q - A
        if q1 < q < q2 and (best is None or lam > best[1]):
            best = (q, lam)
    # dense-grid check of the same maximum (script's growth function)
    ks = np.linspace(math.sqrt(q1), math.sqrt(q2), 20001)
    g = [np.linalg.eigvals(J - np.diag([Du, Dv]) * k * k).real.max() for k in ks]
    assert abs(ks[int(np.argmax(g))] - math.sqrt(best[0])) < 2 * (ks[1] - ks[0]) and abs(max(g) - best[1]) < 1e-8
    return math.sqrt(q1), math.sqrt(q2), math.sqrt(best[0]), best[1]


kmin, kmax, kst, gmax = turing_ref(40.0)
case("turing().kmin", kmin, "rel:1e-12", f"{P12} claim 3 / ch11 bench 3: unstable band, Dv/Du = 40")
case("turing().kmax", kmax, "rel:1e-12", f"{P12} claim 3 / ch11 bench 3: unstable band, Dv/Du = 40")
case("turing().kstar", kst, "rel:1e-7", f"{P12} claim 3: fastest-growing wavenumber (closed form, grid-checked)")
case("turing().wl", 2 * math.pi / kst, "rel:1e-7", "ch11 bench 3: wavelength 2 pi / k*")
case("turing().gmax", gmax, "rel:1e-10", "growth rate at k*")
k2 = turing_ref(20.0)
case("turing(20).kstar", k2[2], "rel:1e-7", "Dv/Du = 20")
case("turing(1).kstar", "nan", "abs:0", f"{P12} claim 3: equal diffusion makes no pattern")
err("turing(40, 0)", "zero diffusivity")

w = "a"; lens = []
for _ in range(20):
    w = "".join("ab" if c == "a" else "a" for c in w); lens.append(len(w))
case("lsys(20).len", lens[-1], "abs:0", f"{P12} claim 5: lengths are Fibonacci numbers")
w7 = "a"
for _ in range(7):
    w7 = "".join("ab" if c == "a" else "a" for c in w7)
case("lsys(7).b", w7.count("b"), "abs:0", f"{P12} claim 5: letter counts of the Fibonacci word")
case("lsys(20).ratio", w.count("a") / w.count("b"), "rel:1e-15", f"{P12} claim 5: a/b tends to the golden ratio")
err("lsys(80)", "lengths beyond 2^53")


def nat_spline(xs, ys):
    n = len(xs) - 1; h = np.diff(xs); A = np.zeros((n + 1, n + 1)); rhs = np.zeros(n + 1)
    A[0, 0] = 1; A[n, n] = 1
    for i in range(1, n):
        A[i, i - 1], A[i, i], A[i, i + 1] = h[i - 1], 2 * (h[i - 1] + h[i]), h[i]
        rhs[i] = 6 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1])
    M = np.linalg.solve(A, rhs)

    def ev(x):
        i = min(max(int(np.searchsorted(xs, x, side="right")) - 1, 0), n - 1); hh = h[i]; a, b = xs[i + 1] - x, x - xs[i]
        return M[i] * a ** 3 / (6 * hh) + M[i + 1] * b ** 3 / (6 * hh) + (ys[i] / hh - M[i] * hh / 6) * a + (ys[i + 1] / hh - M[i + 1] * hh / 6) * b
    energy = sum(h[i] * (M[i] ** 2 + M[i] * M[i + 1] + M[i + 1] ** 2) / 3 for i in range(n))
    return ev, energy, M


xs = np.linspace(0, 6, 8); ys = np.random.default_rng(12).normal(size=8)
ev, en, Msp = nat_spline(xs, ys)
# cross-check the energy by quadrature of f'' (linear on each piece) and the interpolation property
assert abs(en - np.trapezoid(np.interp(np.linspace(0, 6, 600001), xs, Msp) ** 2, np.linspace(0, 6, 600001))) < 1e-6
sx, sy = vec(xs), vec(ys)
case(f"spline({sx}, {sy}, 2.5)", ev(2.5), "rel:1e-12", f"{P12} claim 6: natural cubic spline (dense solve)")
case(f"spline({sx}, {sy}, {num(xs[3])})", ys[3], "rel:1e-12", "the spline interpolates its knots")
case(f"spline({sx}, {sy}, [0.1, 5.9])(2)", ev(5.9), "rel:1e-12", "element-wise evaluation")
case(f"spline({sx}, {sy}, 7)", ev(7.0), "rel:1e-12", "extrapolation continues the end cubic")
case(f"splinee({sx}, {sy})", en, "rel:1e-12", f"{P12} claim 6: least bending energy (Holladay)")
case("spline([0, 1], [0, 2], 0.25)", 0.5, "rel:1e-15", "two knots: a straight line")
err("spline([0, 2, 1], [0, 1, 2], 1)", "knots must increase")


def lloyd(hx=0.7, steps=60):
    LG = 90; r = mulberry32(12); sites = np.array([[r(), r()] for _ in range(12)])
    g = (np.arange(LG) + 0.5) / LG; X, Y = np.meshgrid(g, g)        # rows j (y), columns i (x)
    X, Y = X.ravel(), Y.ravel(); rho = np.exp(-((X - hx) ** 2 + (Y - 0.3) ** 2) / 0.05) + 0.2

    def evaluate(S):
        d2 = (X[:, None] - S[None, :, 0]) ** 2 + (Y[:, None] - S[None, :, 1]) ** 2
        own = d2.argmin(1); bd = d2[np.arange(len(X)), own]
        cent = S.copy(); has = np.zeros(12, bool)
        for k in range(12):
            m = own == k
            if rho[m].sum() > 0:
                has[k] = True; cent[k] = [(rho[m] * X[m]).sum() / rho[m].sum(), (rho[m] * Y[m]).sum() / rho[m].sum()]
        return (rho * bd).sum() / rho.sum(), cent, has

    E0, cent, has = evaluate(sites); E = E0
    for _ in range(steps):
        sites = cent; E, cent, has = evaluate(sites)
    gap = max(np.hypot(*(sites[k] - cent[k])) for k in range(12) if has[k])
    return E0, E, gap


l0, l1, lg = lloyd()
case("lloyd().e0", l0, "rel:1e-9", "ch11 bench 4: coverage energy at the mulberry32(12) start")
case("lloyd().e", l1, "rel:1e-9", f"{P12} claim 4 / ch11 bench 4: energy after 60 Lloyd steps")
case("lloyd().gap", lg, "rel:1e-6", "ch11 bench 4: largest drone-centroid gap (near a fixed point)")
case("lloyd(0.3, 5).e", lloyd(0.3, 5)[1], "rel:1e-9", "hot spot moved to x = 0.3, 5 steps")
err("lloyd(0.7, 1000)", "too many steps")

# =============================================================================== ch12 expansion
P13 = "phase13_expansion.py"
dt, q, r = 1.0, 0.05, 0.5
A = np.array([[1, dt], [0, 1.0]]); Qd = q * np.array([[dt ** 3 / 3, dt ** 2 / 2], [dt ** 2 / 2, dt]]); H = np.array([[1.0, 0]])


def kalman_cov(steps, A=A, Qd=Qd, H=H, R=np.array([[r]]), P=np.diag([4.0, 1.0])):
    for _ in range(steps):
        P = A @ P @ A.T + Qd; P = P - P @ H.T @ np.linalg.inv(H @ P @ H.T + R) @ H @ P
    return P


P40 = kalman_cov(40)
case("pcrbcv().p11", P40[0, 0], "rel:1e-10", f"{P13} claim 1: linear case, PCRB = Kalman covariance (40 steps)")
case("pcrbcv().p12", P40[0, 1], "rel:1e-10", f"{P13} claim 1: linear case, PCRB = Kalman covariance (40 steps)")
case("pcrbcv(5, 0.2, 1).p22", kalman_cov(5, Qd=0.2 * np.array([[1 / 3, 1 / 2], [1 / 2, 1]]), R=np.array([[1.0]]))[1, 1], "rel:1e-10", "PCRB = Kalman, 5 steps, q 0.2, r 1")
F4 = np.eye(4); F4[0, 2] = F4[1, 3] = 1
qq = 0.05 * 0.2; Q4 = np.array([[qq / 3, 0, qq / 2, 0], [0, qq / 3, 0, qq / 2], [qq / 2, 0, qq, 0], [0, qq / 2, 0, qq]])
H4 = np.array([[1.0, 0, 0, 0], [0, 1, 0, 0]]); R4 = 0.25 * np.eye(2); P04 = np.diag([4, 4, 0.25, 0.25])
K30 = kalman_cov(30, A=F4, Qd=Q4, H=H4, R=R4, P=P04)
case(f"P = pcrb({lit(F4)}, {lit(Q4)}, {lit(H4)}, {lit(R4)}, {lit(P04)}, 30); P(1,1) + P(2,2)", K30[0, 0] + K30[1, 1], "rel:1e-10",
     "ch12 bench 1 linear sensor: position bound at step 30 = Kalman")
err("pcrb(eye(2), eye(2), [1, 0, 0], 1, eye(2), 3)", "H has the wrong width")


def Tmin(x, v):
    if x > -v * abs(v) / 2 or (x == -v * abs(v) / 2 and v < 0):
        return v + 2 * math.sqrt(x + v * v / 2)
    return -v + 2 * math.sqrt(-x + v * v / 2)


def switch_time(x, v):
    """First sign change of x + v|v|/2 along the initial full-thrust arc, by bisection on the exact parabola."""
    s0 = x + v * abs(v) / 2; u = -1.0 if s0 > 0 else 1.0
    sw = lambda t: (x + v * t + u * t * t / 2) + (v + u * t) * abs(v + u * t) / 2
    lo, hi = 0.0, Tmin(x, v)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (sw(mid) > 0) == (s0 > 0) and sw(mid) != 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


for x0, v0 in [(1.0, 0.0), (-2.0, 1.0), (0.5, -2.0), (3.0, 2.0)]:
    case(f"bangbang({num(x0)}, {num(v0)}).T", Tmin(x0, v0), "rel:1e-12", f"{P13} claim 2: T(x, v) = v + 2 sqrt(x + v^2/2) (and mirror)")
for x0, v0 in [(1.0, 0.0), (-2.0, 1.0)]:
    case(f"bangbang({num(x0)}, {num(v0)}).ts", switch_time(x0, v0), "rel:1e-9", "switch time: first crossing of the curve x = -v|v|/2")
case("bangbang(1.5, 1).vs", -math.sqrt(1.5 + 0.5), "rel:1e-12", "ch12 bench 2 default: switch velocity")
case("bangbang(2, -2).ts", 0.0, "abs:0", "on the switching curve: no switch")
case("bangbang(0, 0).T", 0.0, "abs:0", "already at rest")


def fly(x, v, dt_=1e-4):
    t = 0.0
    while (x * x + v * v) > 1e-6 and t < 50:
        s = x + v * abs(v) / 2
        u = -1.0 if s > 0 else 1.0 if s < 0 else -math.copysign(1, v)
        v += u * dt_; x += v * dt_; t += dt_
    return t


case("bbsim(1.5, 1)", fly(1.5, 1.0), "rel:1e-9", f"{P13} claim 2 / ch12 bench 2: simulated arrival (dt 1e-4)")
case("bbsim(-2, 1, 1e-3)", fly(-2.0, 1.0, 1e-3), "rel:1e-9", "coarser step")
assert abs(fly(3.0, 2.0) - Tmin(3.0, 2.0)) < 0.01
case("bbsim(3, 2) - bangbang(3, 2).T", fly(3.0, 2.0) - Tmin(3.0, 2.0), "abs:1e-9", f"{P13} claim 2: arrival within 0.01 of T")
err("bbsim(1, 0, 1)", "dt too large")

th, sg, x0_, T_ = 0.8, 0.6, 2.0, 2.0
case("ouproc().mean", x0_ * math.exp(-th * T_), "rel:1e-14", f"{P13} claim 3: mean x0 e^(-theta t)")
case("ouproc().var", sg ** 2 * (1 - math.exp(-2 * th * T_)) / (2 * th), "rel:1e-13", f"{P13} claim 3: variance sigma^2 (1 - e^(-2 theta t)) / (2 theta)")
case("ouproc(1, 0.3, 1, 50).vinf", 1 / 0.6, "rel:1e-14", "stationary variance sigma^2 / (2 theta)")
case("ouproc(0, 0, 0.5, 4).var", 0.25 * 4, "rel:1e-14", "theta = 0: Brownian motion, variance sigma^2 t")
err("ouproc(1, 1, 1, -1)", "negative time")

case("qcv(0.05)(1,1)", Qd[0, 0], "rel:1e-14", f"{P13} claim 4: q dt^3/3")
case("qcv(0.05, 2)(1,2)", 0.05 * 4 / 2, "rel:1e-14", f"{P13} claim 4: q dt^2/2")


def expm_taylor(M):
    s = max(0, int(math.ceil(math.log2(max(np.abs(M).sum(1).max(), 1e-16)))) + 4)
    X = M / 2 ** s; E = np.eye(len(M)); T = np.eye(len(M))
    for k in range(1, 30):
        T = T @ X / k; E = E + T
    for _ in range(s):
        E = E @ E
    return E


def qdisc_ref(A_, Qc, dt_, n=4000):
    """Composite Simpson on int_0^dt e^(A s) Qc e^(A' s) ds."""
    ss = np.linspace(0, dt_, n + 1); wts = np.ones(n + 1); wts[1:-1:2] = 4; wts[2:-1:2] = 2
    acc = np.zeros_like(Qc)
    for s_, w_s in zip(ss, wts):
        E = expm_taylor(A_ * s_); acc = acc + w_s * E @ Qc @ E.T
    return acc * dt_ / (3 * n)


case("qdisc([0, 1; 0, 0], [0, 0; 0, 0.05], 1)(1,2)", Qd[0, 1], "rel:1e-12", f"{P13} claim 4: Van Loan reproduces q dt^2/2")
Asg = np.array([[0, 1], [0, -0.5]]); Qcs = np.array([[0, 0], [0, 0.3]])
Qs = qdisc_ref(Asg, Qcs, 2.0)
case(f"qdisc({lit(Asg)}, {lit(Qcs)}, 2)(1,1)", Qs[0, 0], "rel:1e-10", "Singer-type model, quadrature of e^(As) Qc e^(A's)")
case(f"qdisc({lit(Asg)}, {lit(Qcs)}, 2)(2,2)", Qs[1, 1], "rel:1e-10", "Singer-type model, quadrature of e^(As) Qc e^(A's)")
err("qdisc([0, 1; 0, 0], 1, 1)", "Qc has the wrong size")

gen = np.random.default_rng(13); nn = 7; Wm = np.full((nn, nn), np.inf); np.fill_diagonal(Wm, 0)
for _ in range(16):
    i, j = gen.integers(0, nn, 2)
    if i != j:
        Wm[i, j] = gen.integers(1, 10)
Dfw = Wm.copy()
for k in range(nn):
    Dfw = np.minimum(Dfw, Dfw[:, [k]] + Dfw[[k], :])
Wl = lit(Wm)
fin = [(i, j) for i in range(nn) for j in range(nn) if i != j and np.isfinite(Dfw[i, j])]
for (i, j) in fin[:3]:
    case(f"troppow({Wl})({i + 1},{j + 1})", Dfw[i, j], "abs:0", f"{P13} claim 5: (min, +) power = Floyd-Warshall shortest path")
infs = [(i, j) for i in range(nn) for j in range(nn) if not np.isfinite(Dfw[i, j])]
if infs:
    case(f"troppow({Wl})({infs[0][0] + 1},{infs[0][1] + 1})", "inf", "abs:0", "unreachable pair stays inf")
case(f"troppow({Wl}, 20)({fin[-1][0] + 1},{fin[-1][1] + 1})", Dfw[fin[-1]], "abs:0", "powers beyond n - 1 change nothing (no negative cycles)")
mp = np.min(Wm[:, :, None] + Wm[None, :, :], axis=1)
case(f"tropmul({Wl}, {Wl})(1,{fin[0][1] + 1})", mp[0, fin[0][1]], "abs:0", "(min, +) square: best path of at most 2 edges")
case("tropmul([1, 5], [2; 0])", 3.0, "abs:0", "min(1 + 2, 5 + 0)")
err("tropmul([1, 2], [1, 2])", "inner sizes differ")

m = 6; S0 = gen.random((m, 2)); S1 = gen.random((m, 2)) + [1.5, 0]
Cm = ((S0[:, None] - S1[None]) ** 2).sum(2)
opt = min(sum(Cm[i, p_[i]] for i in range(m)) for p_ in itertools.permutations(range(m)))


def sink(Cm, eps, iters=20000):
    K_ = np.exp(-Cm / eps); u_ = np.ones(len(Cm)); v_ = np.ones(len(Cm))
    for _ in range(iters):
        u_ = 1 / (K_ @ v_); v_ = 1 / (K_.T @ u_)
    return u_[:, None] * K_ * v_[None, :]


Cl = lit(Cm)
for eps in (1.0, 0.1, 0.01):
    case(f"sinkcost({Cl}, {num(eps)})", float((sink(Cm, eps) * Cm).sum()), "rel:1e-9", f"{P13} claim 6: Sinkhorn transport cost, eps = {eps}")
case(f"otexact({Cl})", opt, "rel:1e-12", f"{P13} claim 6: optimal assignment (all 720 permutations)")
case(f"sinkhorn({Cl}, 0.1)(2,3)", sink(Cm, 0.1)[1, 2], "rel:1e-9", "entropic plan entry")
case(f"P = sinkhorn({Cl}, 0.1); P(1,1) + P(2,1) + P(3,1) + P(4,1) + P(5,1) + P(6,1)", 1.0, "abs:1e-12", "unit column marginal")
err("sinkhorn([1, 2; 3, 4], 0)", "epsilon must be positive")
err("otexact([1, 2, 3])", "cost must be square")


def kepler(method, orbits=200, steps_per=200):
    r_ = np.array([1.0, 0.0]); v_ = np.array([0.0, 1.0]); h_ = 2 * math.pi / steps_per; E0 = 0.5 - 1.0; worst = 0.0
    acc = lambda x: -x / math.sqrt(x @ x) ** 3
    for _ in range(orbits * steps_per):
        if method == "leapfrog":
            v_ = v_ + 0.5 * h_ * acc(r_); r_ = r_ + h_ * v_; v_ = v_ + 0.5 * h_ * acc(r_)
        else:
            r_, v_ = r_ + h_ * v_, v_ + h_ * acc(r_)
        E = 0.5 * v_ @ v_ - 1 / math.sqrt(r_ @ r_); worst = max(worst, abs(E - E0))
        if worst > 1:
            break
    return worst


case("kepenergy().lf", kepler("leapfrog"), "rel:1e-6", f"{P13} claim 7: leapfrog energy error over 200 orbits")
case("kepenergy().eu", kepler("euler"), "rel:1e-6", f"{P13} claim 7: explicit Euler drifts")
case("kepenergy(50, 20).eu", kepler("euler", 50, 20), "rel:1e-6", "coarse Euler passes an error of 1 and stops")
err("kepenergy(10000, 10000)", "too much work for a keypress")
c_ = 200 / 800
case("mpedges(0.25).lo", (1 - math.sqrt(c_)) ** 2, "rel:1e-14", f"{P13} claim 8: Marchenko-Pastur lower edge, p/n = 200/800")
case("mpedges(0.25).hi", (1 + math.sqrt(c_)) ** 2, "rel:1e-14", f"{P13} claim 8: Marchenko-Pastur upper edge")
case("mpedges(0.5, 2).hi", 2 * (1 + math.sqrt(0.5)) ** 2, "rel:1e-14", "scaled by the noise variance")
lo_, hi_ = (1 - math.sqrt(c_)) ** 2, (1 + math.sqrt(c_)) ** 2
xg_ = np.linspace(lo_, hi_, 200001)
dens_mp = np.sqrt((hi_ - xg_) * (xg_ - lo_)) / (2 * math.pi * c_ * xg_)
assert abs(np.trapezoid(dens_mp, xg_) - 1) < 1e-6
case("mppdf(1, 0.25)", math.sqrt((hi_ - 1) * (1 - lo_)) / (2 * math.pi * c_), "rel:1e-14", "Marchenko-Pastur density at 1 (integrates to 1)")
case("mppdf(3, 0.25)", 0.0, "abs:0", "outside the support")
err("mpedges(0)", "aspect ratio must be positive")

print("# Generated by Tests/reference/gen_fields.py from Where It Isn't chapters 7, 10, 11, 12 (inline benches and verification/python/phase8/11/12/13).")
print("# expression\texpected\ttolerance\tprovenance")
print("\n".join(ROWS))
