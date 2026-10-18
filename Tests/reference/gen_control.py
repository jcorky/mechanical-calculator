"""Reference cases for the control library that Where It Isn't checks only in Python (numpy).

Run by gen_control.mjs, which appends these lines to control.tsv:
    python3 Tests/reference/gen_control.py /path/to/where-it-isnt

The computations follow verification/python/phase4_mobius.py (the tracker's symplectic matrix, Kalman-LQR
duality, the Hamiltonian halves) and phase9_dynamics.py (the Lyapunov equation, the region of attraction),
written out again here with numpy rather than imported, since those scripts only assert.
"""
import math
import sys

import numpy as np

wii = sys.argv[1] if len(sys.argv) > 1 else "."


def num(x):
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return repr(x)


def mat(M):
    M = np.atleast_2d(np.asarray(M, dtype=float))
    return "[" + "; ".join(", ".join(num(v) for v in row) for row in M) + "]"


def emit(expr, expected, tol, prov):
    print(f"{expr}\t{expected if isinstance(expected, str) else num(expected)}\t{tol}\t{prov}")


# ---- ch03: the constant-velocity tracker (bench C defaults and the three phase-4 configurations) ----
def tracker(dt, qq, rr):
    A = np.array([[1, dt], [0, 1.0]])
    C = np.array([[1.0, 0]])
    Q = qq * np.array([[dt ** 3 / 3, dt ** 2 / 2], [dt ** 2 / 2, dt]])
    R = np.array([[rr]])
    P = np.eye(2)
    for _ in range(5000):
        S = C @ P @ C.T + R
        P = A @ P @ A.T - A @ P @ C.T @ np.linalg.inv(S) @ C @ P @ A.T + Q
    Kf = P @ C.T @ np.linalg.inv(C @ P @ C.T + R)
    return A, C, Q, R, P, Kf


configs = [(0.5, 10 ** -0.5, 10 ** -0.1, "ch03 bench C defaults"), (0.5, 0.3, 0.8, "phase4 claim 2"), (1.5, 10.0, 0.01, "phase4 claim 2"), (0.1, 0.01, 1.0, "phase4 claim 2")]
for dt, qq, rr, src in configs:
    A, C, Q, R, P, Kf = tracker(dt, qq, rr)
    setup = f"A = {mat(A)}; C = {mat(C)}; Q = {mat(Q)}; R = {num(rr)}; "
    prov = f"numpy: {src}, dt {dt}, q {qq:.6g}, r {rr:.6g}"
    emit(setup + "kalman(A, C, Q, R)(1, 1)", Kf[0, 0], "rel:1e-9", prov + ", filter gain")
    emit(setup + "kalman(A, C, Q, R)(2, 1)", Kf[1, 0], "rel:1e-9", prov + ", filter gain")
    emit(setup + "kfcov(A, C, Q, R)(1, 2)", P[0, 1], "rel:1e-9", prov + ", prior covariance")
    poles = np.linalg.eigvals(A @ (np.eye(2) - Kf @ C))
    Ac, Bc = A.T, C.T
    G = Bc @ np.linalg.inv(R) @ Bc.T
    Ait = np.linalg.inv(Ac).T
    Z = np.block([[Ac + G @ Ait @ Q, -G @ Ait], [-Ait @ Q, Ait]])
    ev = np.linalg.eigvals(Z)
    emit(setup + "rho(sympz(A, C, Q, R))", np.abs(ev).max(), "rel:1e-8", prov + ", largest |eig Z| = 1/|estimator pole|")
    emit(setup + "1/rho(sympz(A, C, Q, R))", np.abs(poles).min(), "rel:1e-8", prov + ", Vaughan: stable half of Z = estimator poles")
    emit(setup + "rho(A*(eye(2) - kalman(A, C, Q, R)*C))", np.abs(poles).max(), "rel:1e-8", prov + ", estimator poles")
    emit(setup + "sympz(A, C, Q, R)(1, 3)", Z[0, 2], "rel:1e-12", prov + ", Z entry")
    emit(setup + "det(sympz(A, C, Q, R))", np.linalg.det(Z), "abs:1e-9", prov + ", symplectic: det Z = 1")
    Jm = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    emit(setup + "Z = sympz(A, C, Q, R); J = [0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, 0; 0, -1, 0, 0]; norm(Z'*J*Z - J, 0)", 0.0, "abs:1e-9", prov + ", Z'JZ = J")
    Sd = np.eye(2)
    for _ in range(5000):
        Sd = Ac.T @ Sd @ Ac - Ac.T @ Sd @ Bc @ np.linalg.inv(Bc.T @ Sd @ Bc + R) @ Bc.T @ Sd @ Ac + Q
    Klqr = np.linalg.inv(Bc.T @ Sd @ Bc + R) @ Bc.T @ Sd @ Ac
    emit(setup + "dlqr(A', C', Q, R)(1, 1)", Klqr[0, 0], "rel:1e-9", prov + ", dual LQR gain")
    emit(setup + "dlqr(A', C', Q, R)(1, 2)", (A @ Kf).T[0, 1], "rel:1e-9", prov + ", duality: dual LQR gain = (A K_f)'")

emit("sympz([0, 1; 0, 0], [1, 0], eye(2), 1)", "error", "abs:0", "A singular: no symplectic form")

# ---- ch03: the continuous ARE from the Hamiltonian (phase4 claim 5) ----
A = np.array([[0, 1.0], [-2, -0.3]]); C = np.array([[1.0, 0]]); Q = np.diag([0.1, 0.5]); R = np.array([[0.2]])
H = np.block([[A, Q], [C.T @ np.linalg.inv(R) @ C, -A.T]])
w, V = np.linalg.eig(H)
setup = f"A = {mat(A)}; C = {mat(C)}; Q = {mat(Q)}; R = 0.2; "
prov = "numpy: phase4 claim 5, damped oscillator"
emit(setup + "hamilton(A, C, Q, R)(3, 1)", H[2, 0], "rel:1e-12", prov + ", H = [A Q; C'R^-1C -A']")
emit(setup + "hamilton(A, C, Q, R)(4, 3)", H[3, 2], "rel:1e-12", prov)
emit(setup + "rho(hamilton(A, C, Q, R))", np.abs(w).max(), "rel:1e-9", prov + ", eigenvalues of H")
for side, sel in [(1, w.real > 0), (-1, w.real < 0)]:
    X, Y = V[:2, sel], V[2:, sel]
    P = np.real(X @ np.linalg.inv(Y))
    half = "right half (stabilizing, P > 0)" if side == 1 else "left half (P < 0)"
    for (i, j) in [(0, 0), (0, 1), (1, 1)]:
        emit(setup + f"hamare(A, C, Q, R, {side})({i + 1}, {j + 1})", P[i, j], "rel:1e-8", prov + ", " + half)
    cl = np.linalg.eigvals(A - P @ C.T @ np.linalg.inv(R) @ C)
    emit(setup + f"eig(A - hamare(A, C, Q, R, {side})*C'*C/0.2)(1, 1)", cl.real.max(), "rel:1e-8", prov + ", error poles " + ("-0.653 +- 1.46i (page)" if side == 1 else "unstable"))
emit(setup + "hamare(A, C, Q, R, 0)", "error", "abs:0", "side must be +-1")
emit("hamilton([1, 2; 3, 4], [1, 0, 0], eye(2), 1)", "error", "abs:0", "C has the wrong number of columns")
emit("hamare([1, 0; 0, 1], [0, 0], eye(2), 1)", "error", "abs:0", "unobservable unstable modes: no graph subspace")

# ---- ch08: the Lyapunov equation (phase9 claim 1), A'P + PA = -Q as lyap(A', Q) ----
def clyap(A, Q):
    n = len(A)
    return np.linalg.solve(np.kron(np.eye(n), A.T) + np.kron(A.T, np.eye(n)), -Q.ravel()).reshape(n, n)


A = np.array([[-1.2, 0.4, 0.3], [0.5, -2.0, 0.7], [-0.3, 0.6, -0.9]])
Pl = clyap(A, np.eye(3))
setup = f"A = {mat(A)}; "
emit(setup + "lyap(A', eye(3))(1, 2)", Pl[0, 1], "rel:1e-10", "numpy: phase9 claim 1, Kronecker solve of A'P + PA = -I")
emit(setup + "lyap(A', eye(3))(3, 3)", Pl[2, 2], "rel:1e-10", "numpy: phase9 claim 1")
emit(setup + "min(eig(lyap(A', eye(3))))", np.linalg.eigvalsh(Pl).min(), "rel:1e-9", "numpy: phase9 claim 1, Hurwitz A gives P > 0")
ts = np.linspace(0, 40, 8001); ww, VV = np.linalg.eig(A); Vi = np.linalg.inv(VV)
E = lambda t: (VV @ np.diag(np.exp(ww * t)) @ Vi).real
Pint = sum(E(t).T @ E(t) for t in ts[1:-1]) * (ts[1] - ts[0]) + 0.5 * (ts[1] - ts[0]) * (np.eye(3) + E(40).T @ E(40))
emit(setup + "lyap(A', eye(3))(2, 2)", Pint[1, 1], "abs:1e-3", "numpy: phase9 claim 1, P = integral of e^{A't} e^{At} (trapezoid)")
A2 = np.array([[0, -1], [1, -1.0]])
for q in (1.0, 0.25, 4.0):
    P = clyap(A2, np.diag([1, q]))
    emit(f"lyap([0, -1; 1, -1]', [1, 0; 0, {num(q)}])(1, 1)", P[0, 0], "rel:1e-12", f"numpy: phase9 claim 2, Van der Pol linearization, Q = diag(1, {q})")
F = np.array([[0.9, 0.2], [-0.1, 0.7]]); W = np.array([[1.0, 0.3], [0.3, 2.0]])
X = np.linalg.solve(np.eye(4) - np.kron(F, F), W.ravel()).reshape(2, 2)
emit(f"dlyap({mat(F)}, {mat(W)})(1, 2)", X[0, 1], "rel:1e-11", "numpy: X = F X F' + W by the Kronecker form")
emit(f"dlyap({mat(F)}, {mat(W)})(2, 2)", X[1, 1], "rel:1e-11", "numpy: X = F X F' + W")
emit("lyap([1, 2; 3, 4], eye(3))", "error", "abs:0", "Q must match A")
emit("dlyap([1, 2, 3], 1)", "error", "abs:0", "A must be square")

# ---- ch08: region of attraction, c* in closed form per ray, minimized densely ----
def roa_c(q):
    P = clyap(A2, np.diag([1, q]))
    th = np.linspace(0, np.pi, 2_000_001)
    c, s = np.cos(th), np.sin(th)
    den = 2 * c * c * s * (P[0, 1] * c + P[1, 1] * s)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = (c * c + q * s * s) / den
        lev = np.where((den > 0) & (r2 <= 16), r2 * (P[0, 0] * c * c + 2 * P[0, 1] * c * s + P[1, 1] * s * s), np.inf)
    k = int(np.argmin(lev))

    def g(t):
        cc, ss = math.cos(t), math.sin(t)
        d = 2 * cc * cc * ss * (P[0, 1] * cc + P[1, 1] * ss)
        r = (cc * cc + q * ss * ss) / d
        return r * (P[0, 0] * cc * cc + 2 * P[0, 1] * cc * ss + P[1, 1] * ss * ss)

    lo, hi = th[k - 1], th[k + 1]
    gr = (math.sqrt(5) - 1) / 2
    for _ in range(200):
        x1, x2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if g(x1) < g(x2):
            hi = x2
        else:
            lo = x1
    cs = min(g(0.5 * (lo + hi)), lev[k])
    return cs, math.pi * cs / math.sqrt(np.linalg.det(P))


for q in (1.0, 0.25, 0.5, 2.0, 4.0):
    cs, area = roa_c(q)
    emit(f"roa({num(q)}).c", cs, "rel:1e-9", f"numpy: phase9 claim 2, min over rays of the first Vdot = 0 level, Q = diag(1, {q})")
    emit(f"roa({num(q)}).area", area, "rel:1e-9", f"numpy: phase9 claim 2, pi c / sqrt(det P), Q = diag(1, {q})")
