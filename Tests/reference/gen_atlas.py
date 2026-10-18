#!/usr/bin/env python3
"""Reference cases for LibAtlas.cpp (atlas chapters 1, 2, 4, 5, 6).

    python3 Tests/reference/gen_atlas.py /home/user/where-it-isnt > Tests/reference/atlas.tsv

These chapters have no standalone JS module; their math lives in each page's inline script and in
verification/python/phase{2,3,5,6,7}_*.py. Expected values are computed here independently with numpy:
matrix exponentials by power series, eigen/SVD solutions with numpy.linalg, KKT systems solved as in
the verification scripts, LPs by vertex enumeration, integrals by fine quadrature, and the chapters'
seeded simulations ported line by line (mulberry32 bit for bit). Nothing here calls the C++.
"""
import itertools
import math
import os
import sys

import numpy as np

WII = sys.argv[1] if len(sys.argv) > 1 else "/home/user/where-it-isnt"
if not os.path.isdir(os.path.join(WII, "chapters")):
    sys.exit(f"Where It Isn't not found at {WII}")

D2R = math.pi / 180
rows = []


def num(x):
    if isinstance(x, (bool, np.bool_)):
        return "1" if x else "0"
    x = float(x)
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return repr(x)


def mat(M):
    M = np.atleast_2d(np.asarray(M, float))
    return "[" + "; ".join(", ".join(num(v) for v in r) for r in M) + "]"


def vec(v):
    return "[" + ", ".join(num(x) for x in np.ravel(v)) + "]"


def case(expr, expected, tol="rel:1e-9", prov=""):
    if isinstance(expected, str):
        e = expected
    else:
        e = num(expected)
        if tol.startswith("rel:") and float(expected) == 0.0:
            tol = "abs:1e-12"
    rows.append(f"{expr}\t{e}\t{tol}\t{prov}")


def err(expr, prov):
    rows.append(f"{expr}\terror\tabs:0\t{prov}")


# ------------------------------------------------------------------ helpers (numpy, independent)
def hat(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], float)


def expm_series(M, terms=60):
    out, term = np.eye(len(M)), np.eye(len(M))
    for k in range(1, terms):
        term = term @ M / k
        out = out + term
    return out


def rot(n, t):
    n = np.asarray(n, float) / np.linalg.norm(n)
    # scaling and squaring keeps the series accurate for large angles
    s = 8
    return np.linalg.matrix_power(expm_series(t / 2 ** s * hat(n)), 2 ** s)


def unit(x):
    x = np.asarray(x, float)
    return x / np.linalg.norm(x)


def Rq(q):
    q = unit(q)
    w, v = q[0], q[1:]
    return (w * w - v @ v) * np.eye(3) + 2 * np.outer(v, v) + 2 * w * hat(v)


def qmul(a, b):
    # left-multiplication matrix of a (Hamilton, scalar first)
    w, x, y, z = a
    L = np.array([[w, -x, -y, -z], [x, w, -z, y], [y, z, w, -x], [z, -y, x, w]], float)
    return L @ np.asarray(b, float)


def canon(q):
    q = np.asarray(q, float)
    big = np.argmax(np.abs(q))
    if q[0] < 0 or (q[0] == 0 and q[big] < 0):
        q = -q
    return q


def axis_angle_q(n, t):
    n = unit(n)
    return np.r_[np.cos(t / 2), np.sin(t / 2) * n]


def angle_of(R):
    v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    return math.atan2(np.linalg.norm(v) / 2, (np.trace(R) - 1) / 2)


def K_of(rs, bs, ws):
    B = sum(w * np.outer(b, r) for r, b, w in zip(rs, bs, ws))
    S, s = B + B.T, np.trace(B)
    z = sum(w * np.cross(r, b) for r, b, w in zip(rs, bs, ws))
    K = np.zeros((4, 4))
    K[0, 0], K[0, 1:], K[1:, 0], K[1:, 1:] = s, z, z, S - s * np.eye(3)
    return K, B


def triad_np(r1, r2, b1, b2):
    def fr(x, y):
        t1, t2 = unit(x), unit(np.cross(x, y))
        return np.c_[t1, t2, np.cross(t1, t2)]
    return fr(b1, b2) @ fr(r1, r2).T


# ------------------------------------------------------------------ ch01 groups
P1 = "ch01-groups"
case("hat([1, 2, 3])(1, 2)", -3, prov=P1 + " hat map")
case("hat([1, 2, 3])(3, 1)", -2, prov=P1 + " hat map")
err("hat([1, 2])", P1 + " hat needs a 3-vector")
case("vee([0, -3, 2; 3, 0, -1; -2, 1, 0])(1)", 1, prov=P1 + " vee inverts hat")
case("vee(hat([0.3, -0.2, 0.7]))(3)", 0.7, prov=P1 + " vee(hat(w)) = w")
err("vee([1, 2, 3])", P1 + " vee needs 3x3")

case("rodrigues([0, 0, 1], pi/2)(2, 1)", 1.0, prov=P1 + " Theorem 1, quarter turn about z")
n_b = [math.cos(35 * D2R) * math.cos(30 * D2R), math.cos(35 * D2R) * math.sin(30 * D2R), math.sin(35 * D2R)]
Rb = expm_series(120 * D2R * hat(n_b))
case("rodrigues([cos(35*pi/180)*cos(30*pi/180), cos(35*pi/180)*sin(30*pi/180), sin(35*pi/180)], 120*pi/180)(1, 3)", Rb[0, 2],
     prov=P1 + " bench 1 defaults (az 30, el 35, theta 120 deg) vs exponential series (phase2 claim 1)")
case("rodrigues([2, -1, 0.5], 1.3)(3, 2)", expm_series(1.3 * hat(unit([2, -1, 0.5])))[2, 1], prov=P1 + " Rodrigues = exp series, unnormalized axis")
err("rodrigues([0, 0, 0], 1)", P1 + " zero axis")

w = np.array([0.1, 0.2, 0.3])
case("exp_so3([0.1, 0.2, 0.3])(2, 3)", expm_series(hat(w))[1, 2], prov=P1 + " exponential map, series")
case("exp_so3([0, 0, 0])(1, 1)", 1.0, prov=P1 + " exp(0) = I")

case("log_so3(expm(hat([0.3, -0.5, 0.4])))(2)", -0.5, prov=P1 + " log inverts exp")
case("log_so3(rodrigues([1, 2, 2], 3.0))(3)", 2.0, prov=P1 + " log near pi: theta n = 3 * (1,2,2)/3")
case("norm(log_so3(rodrigues([1, 2, 2], pi)))", math.pi, prov=P1 + " log at exactly pi")
err("log_so3([1, 2, 0; 3, 4, 0; 0, 0, 1])", P1 + " not a rotation")

case("rot_angle(rodrigues([1, 2, 3], 2.5))", 2.5, prov=P1 + " angle of R")
case("rot_angle(rodrigues([1, 0, 0], 4))", 2 * math.pi - 4, prov=P1 + " bench 1: angle read back is always in [0, pi]")

case("quat([0, 0, 1], pi/2)(1)", math.cos(math.pi / 4), prov=P1 + " Theorem 2: q = (cos t/2, sin t/2 n)")
case("quat([1, -1, 2], 2)(4)", axis_angle_q([1, -1, 2], 2)[3], prov=P1 + " Theorem 2")
err("quat([0, 0, 0], 1)", P1 + " zero axis")

case("quat_mul([1, 0, 0, 0], [0, 1, 0, 0])(2)", 1.0, prov=P1 + " identity times i")
case("quat_mul([0, 1, 0, 0], [0, 0, 1, 0])(4)", 1.0, prov=P1 + " Hamilton: i j = k")
case("quat_mul([0, 0, 1, 0], [0, 1, 0, 0])(4)", -1.0, prov=P1 + " Hamilton: j i = -k")
pa, pb = np.array([0.5, 0.5, 0.5, 0.5]), np.array([0.1, -0.2, 0.3, 0.9])
case("quat_mul([0.5, 0.5, 0.5, 0.5], [0.1, -0.2, 0.3, 0.9])(3)", qmul(pa, pb)[2], prov=P1 + " Hamilton product (matrix form)")
err("quat_mul([1, 0, 0], [0, 1, 0, 0])", P1 + " quaternions need 4 entries")

case("quat_conj([1, 2, 3, 4])(3)", -3, prov=P1 + " conjugate")
case("quat_mul(quat([1, 2, 3], 1.1), quat_conj(quat([1, 2, 3], 1.1)))(1)", 1.0, prov=P1 + " q q* = 1")

case("quat_rot(quat([0, 0, 1], pi/2), [1, 0, 0])(2)", 1.0, prov=P1 + " rotate x by 90 deg about z")
qq, vv = np.array([1, 2, 3, 4.0]), np.array([0.5, -1, 2])
case("quat_rot([1, 2, 3, 4], [0.5, -1, 2])(1)", (Rq(qq) @ vv)[0], prov=P1 + " R(q) v (q normalized)")
err("quat_rot([1, 0, 0, 0], [1, 2])", P1 + " v must be a 3-vector")

qn, tn = [1, -1, 2], 2.0
case("quat2dcm(quat([1, -1, 2], 2))(1, 3)", expm_series(tn * hat(unit(qn)))[0, 2], prov=P1 + " R(q) = exp(theta [n]x) (phase2 claim 2)")
case("quat2dcm(-quat([1, -1, 2], 2))(1, 3)", expm_series(tn * hat(unit(qn)))[0, 2], prov=P1 + " double cover R(-q) = R(q)")
err("quat2dcm([0, 0, 0, 0])", P1 + " zero quaternion")

case("dcm2quat(rodrigues([1, -1, 2], 2))(2)", axis_angle_q([1, -1, 2], 2)[1], prov=P1 + " matrix to quaternion")
case("dcm2quat(rodrigues([0, 0, 1], 5))(4)", -math.sin(2.5), prov=P1 + " 5 rad about z: w >= 0 picks -q")
err("dcm2quat([1, 0, 0; 0, 1, 0; 0, 0, -1])", P1 + " reflection is not in SO(3)")

case("slerp([1, 0, 0, 0], quat([0, 0, 1], 2), 0.25)(4)", math.sin(0.25), prov=P1 + " slerp along a great circle")
case("slerp([1, 0, 0, 0], -quat([0, 0, 1], 1), 0.5)(1)", math.cos(0.25), prov=P1 + " slerp takes the short way (q ~ -q)")

case("bracket([1, 0, 0], [0, 1, 0])(3)", 1.0, prov=P1 + " [e_x, e_y] = e_z")
a3, b3 = np.array([1, 2, 3.0]), np.array([-1, 0.5, 2])
case("bracket([1, 2, 3], [-1, 0.5, 2])(1)", np.cross(a3, b3)[0], prov=P1 + " bracket = cross product (phase2 claim 3)")
Am, Bm = np.array([[1, 2], [3, 4.0]]), np.array([[0, 1], [1, 0.0]])
case("bracket([1, 2; 3, 4], [0, 1; 1, 0])(1, 1)", (Am @ Bm - Bm @ Am)[0, 0], prov=P1 + " matrix commutator")
err("bracket([1, 2, 3], [1, 2; 3, 4])", P1 + " mixed shapes")


def comm(al, be):
    Rx, Ry = rot([1, 0, 0], al), rot([0, 1, 0], be)
    return angle_of(Rx @ Ry @ Rx.T @ Ry.T)


case("comm_angle()", comm(60 * D2R, 60 * D2R), prov=P1 + " counterexample 1 bench defaults 60, 60 deg")
case("comm_angle(0.05, 0.08)", comm(0.05, 0.08), prov=P1 + " small angles: ~ alpha beta (phase2 claim 6)")
case("comm_angle(0.05, 0.08)/(0.05*0.08)", comm(0.05, 0.08) / 0.004, "rel:1e-9", P1 + " ratio to alpha beta")

# Wahba / Davenport / SVD / TRIAD on a fixed small problem
rng = np.random.default_rng(1)
At = Rq(unit(rng.normal(size=4)))
rs = np.array([unit(x) for x in rng.normal(size=(4, 3))])
bs_exact = np.array([At @ r for r in rs])
bs = np.array([unit(At @ r + 0.01 * rng.normal(size=3)) for r in rs])   # unit, as the chapter's measurements
ws = np.array([1.0, 2.0, 0.5, 1.5])
K, B = K_of(rs, bs, ws)
R_, Bx, W_ = mat(rs), mat(bs), vec(ws)
case(f"davenport({R_}, {Bx}, {W_})(1, 1)", K[0, 0], prov=P1 + " Davenport K: sigma = tr B")
case(f"davenport({R_}, {Bx}, {W_})(2, 4)", K[1, 3], prov=P1 + " Davenport K (Hamilton)")
qq_test = unit([0.3, -0.1, 0.8, 0.2])
case(f"q = {vec(qq_test)}; q*davenport({R_}, {Bx}, {W_})*q'", qq_test @ K @ qq_test, prov=P1 + " q'Kq")
case(f"q = {vec(qq_test)}; trace(quat2dcm(q)*({mat(B)})')", qq_test @ K @ qq_test, prov=P1 + " q'Kq = tr(R(q) B') (phase2 claim 4)")
err("davenport([1, 0, 0; 0, 1, 0], [1, 0, 0], 1)", P1 + " r and b sizes differ")

qo = canon(np.linalg.eigh(K)[1][:, -1])
case(f"qmethod({R_}, {Bx}, {W_})(1)", qo[0], prov=P1 + " q-method: top eigenvector (numpy eigh)")
case(f"qmethod({R_}, {Bx}, {W_})(3)", qo[2], prov=P1 + " q-method")
q_true = canon(np.array([0.3, -0.5, 0.7, 0.4]) / np.linalg.norm([0.3, -0.5, 0.7, 0.4]))
bs_true = np.array([Rq(q_true) @ r for r in rs])
case(f"qmethod({R_}, {mat(bs_true)})(2)", q_true[1], prov=P1 + " noiseless: recovers the true quaternion")

U, _, Vt = np.linalg.svd(B)
Asvd = U @ np.diag([1, 1, np.linalg.det(U) * np.linalg.det(Vt)]) @ Vt
case(f"wahba_svd({R_}, {Bx}, {W_})(2, 3)", Asvd[1, 2], prov=P1 + " SVD solution (Markley 1988)")
case(f"wahba_svd({R_}, {Bx}, {W_})(1, 1)", Rq(qo)[0, 0], prov=P1 + " SVD = q-method (phase2 claim 4)")
case(f"wahba_svd({R_}, {mat(bs_exact)})(3, 1)", At[2, 0], prov=P1 + " noiseless: recovers A")
loss = 0.5 * sum(w_ * np.sum((b - Asvd @ r) ** 2) for r, b, w_ in zip(rs, bs, ws))
case(f"wahba_loss(wahba_svd({R_}, {Bx}, {W_}), {R_}, {Bx}, {W_})", loss, prov=P1 + " Wahba loss at the optimum")
lam_max = np.linalg.eigvalsh(K)[-1]
case(f"wahba_loss(wahba_svd({R_}, {Bx}, {W_}), {R_}, {Bx}, {W_})", ws.sum() - lam_max, prov=P1 + " loss = lambda0 - lambda_max (ch06 guises)")

r1, r2 = rs[0], rs[1]
b1, b2 = bs[0], bs[1]
At3 = triad_np(r1, r2, b1, b2)
case(f"triad({vec(r1)}, {vec(r2)}, {vec(b1)}, {vec(b2)})(1, 2)", At3[0, 1], prov=P1 + " TRIAD (Black 1964)")
case(f"triad({vec(r1)}, {vec(r2)}, {vec(bs_exact[0])}, {vec(bs_exact[1])})(2, 2)", At[1, 1], prov=P1 + " TRIAD exact with noiseless vectors")
err("triad([1, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 1])", P1 + " parallel reference vectors")


def qz(d):
    return np.array([math.cos(d * D2R / 2), 0, 0, math.sin(d * D2R / 2)])


case(f"abs(quat_mean([{', '.join(num(x) for x in qz(179))}; {', '.join(num(x) for x in qz(-179))}])(4))", 1.0, "abs:1e-12",
     P1 + " counterexample 2: mean of 179 and -179 deg is 180 deg")
Qs = np.array([unit([1, 0.1, -0.2, 0.05]), -unit([0.9, 0.15, -0.1, 0.0]), unit([1, 0.05, -0.25, 0.1])])
wq = np.array([1.0, 2.0, 0.5])
Mq = sum(w_ * np.outer(q, q) for q, w_ in zip(Qs, wq))
qm = canon(np.linalg.eigh(Mq)[1][:, -1])
case(f"quat_mean({mat(Qs)}, {vec(wq)})(3)", qm[2], prov=P1 + " eigenvector mean, insensitive to the sign of each q")

tw = np.zeros((4, 4))
tw[:3, :3], tw[:3, 3] = hat([0.2, -0.4, 0.9]), [1.0, 2.0, -0.5]
G = expm_series(tw)
case("exp_se3([0.2, -0.4, 0.9], [1, 2, -0.5])(1, 4)", G[0, 3], prov=P1 + " SE(3) exponential (series of the 4x4 twist)")
case("exp_se3([0.2, -0.4, 0.9], [1, 2, -0.5])(2, 1)", G[1, 0], prov=P1 + " SE(3) exponential")
case("exp_se3([0, 0, 0], [1, 2, 3])(3, 4)", 3.0, prov=P1 + " pure translation")

case("att_err(rodrigues([0, 0, 1], 0.3), eye(3))(3)", math.sin(0.3), prov=P1 + " e_R = sin(theta) n for Rd = I")
Rr, Rd = rot([1, 2, -1], 0.7), rot([0.5, -1, 1], 0.4)
eR = 0.5 * np.array([(Rd.T @ Rr - Rr.T @ Rd)[2, 1], (Rd.T @ Rr - Rr.T @ Rd)[0, 2], (Rd.T @ Rr - Rr.T @ Rd)[1, 0]])
case("att_err(rodrigues([1, 2, -1], 0.7), rodrigues([0.5, -1, 1], 0.4))(1)", eR[0], prov=P1 + " geometric control error (Lee et al. 2010)")

# ------------------------------------------------------------------ ch02 topology
P2 = "ch02-topology"


def Emat(ph, th):
    return np.array([[1, np.sin(ph) * np.tan(th), np.cos(ph) * np.tan(th)],
                     [0, np.cos(ph), -np.sin(ph)],
                     [0, np.sin(ph) / np.cos(th), np.cos(ph) / np.cos(th)]])


rates = Emat(20 * D2R, 60 * D2R) @ np.array([0, 0, 0.2])
case("euler_rates().dpsi", rates[2], prov=P2 + " bench defaults theta 60, phi 20 deg, r 0.2 rad/s")
case("euler_rates().dphi", rates[0], prov=P2 + " bench defaults")
rates2 = Emat(0.1, 0.2) @ np.array([0.3, -0.2, 0.5])
case("euler_rates(0.1, 0.2, 0.3, -0.2, 0.5).dtheta", rates2[1], prov=P2 + " Euler kinematics")
err("euler_rates(0, pi/2)", P2 + " gimbal lock")
case("euler_det()", np.linalg.det(Emat(0, 60 * D2R)), prov=P2 + " det E = 1/cos(theta) = 2 at 60 deg (phase3 claim 1)")
case("euler_det(1.5, 2.0)", 1 / math.cos(1.5), prov=P2 + " det E independent of phi")
case("euler_cond()", np.linalg.cond(Emat(20 * D2R, 60 * D2R)), prov=P2 + " condition number, bench defaults")
case("euler_cond(0, 89*pi/180)", np.linalg.cond(Emat(0, 89 * D2R)), prov=P2 + " near gimbal lock")


def eul(ph, th, ps):
    return rot([0, 0, 1], ps) @ rot([0, 1, 0], th) @ rot([1, 0, 0], ph)


case("eul2dcm(0.1, 0.2, 0.3)(3, 1)", -math.sin(0.2), prov=P2 + " 3-2-1: R31 = -sin(theta)")
case("eul2dcm(0.1, 0.2, 0.3)(1, 2)", eul(0.1, 0.2, 0.3)[0, 1], prov=P2 + " Rz(psi) Ry(theta) Rx(phi)")
case("dcm2eul(eul2dcm(0.1, -0.4, 2.5)).psi", 2.5, prov=P2 + " round trip")
case("dcm2eul(eul2dcm(-2.0, 1.2, -0.7)).phi", -2.0, prov=P2 + " round trip")
err("dcm2eul([1, 0, 0; 0, 2, 0; 0, 0, 1])", P2 + " not a rotation")

case("lookat([1, 0, 0])(3)", 1.0, prov=P2 + " up-vector for a horizontal boresight")
d = unit([1, 1, 1])
u = np.array([0, 0, 1.0]) - d[2] * d
case("lookat([1, 1, 1])(1)", unit(u)[0], prov=P2 + " u(d) = (z - (z.d)d)/|.|")
err("lookat([0, 0, 2])", P2 + " hairy ball: vanishes at d = z")

OBS = [((0.0, 0.0), 0.25), ((-0.5, 0.35), 0.15), ((0.5, 0.35), 0.15)]


def nav_eval(x, y, k, goal=(0.0, 0.6), big=0.25):
    obs = [((0.0, 0.0), big)] + OBS[1:]
    gx0, gy0 = goal
    g = (x - gx0) ** 2 + (y - gy0) ** 2
    dgx, dgy = 2 * (x - gx0), 2 * (y - gy0)
    b, bx, by = 1 - (x * x + y * y), -2 * x, -2 * y
    for (cx, cy), r in obs:
        t, tx, ty = (x - cx) ** 2 + (y - cy) ** 2 - r * r, 2 * (x - cx), 2 * (y - cy)
        bx, by, b = bx * t + b * tx, by * t + b * ty, b * t
    with np.errstate(all="ignore"):
        D = np.power(g, k) + b
        A, Bm = np.power(D, -1 / k), np.power(D, -1 / k - 1)
        f = g * A
        fx = dgx * A - g / k * Bm * (k * np.power(g, k - 1) * dgx + bx)
        fy = dgy * A - g / k * Bm * (k * np.power(g, k - 1) * dgy + by)
    return f, fx, fy, b


f, fx, fy, _ = nav_eval(0.3, -0.5, 2)
case("navfn(0.3, -0.5).phi", f, prov=P2 + " Rimon-Koditschek phi, bench world (k 2, goal (0, 0.6))")
f4, fx4, fy4, _ = nav_eval(0.3, -0.5, 4)
case("navfn(0.3, -0.5, 4).fx", fx4, prov=P2 + " gradient, k = 4")
f5 = nav_eval(-0.2, 0.8, 6, (0.2, 0.5), 0.2)
case("navfn(-0.2, 0.8, 6, 0.2, 0.5, 0.2).fy", f5[2], prov=P2 + " moved goal and obstacle")
case("navfn(0, 0.6).phi", 0.0, "abs:1e-15", P2 + " phi = 0 at the goal")
err("navfn(0, 0)", P2 + " inside the big obstacle")


def nav_index(k, goal=(0.0, 0.6), big=0.25, n=220):
    xs = -1 + 2 * np.arange(n) / (n - 1)
    X, Y = np.meshgrid(xs, xs)          # row i -> y, column j -> x, as the bench
    _, fx_, fy_, b = nav_eval(X, Y, k, goal, big)
    ok = b > 0
    ang = np.where(ok, np.arctan2(fy_, fx_), 0.0)

    def d(a, c):
        t = c - a
        t = np.where(t > np.pi, t - 2 * np.pi, t)
        return np.where(t < -np.pi, t + 2 * np.pi, t)
    s = d(ang[:-1, :-1], ang[:-1, 1:]) + d(ang[:-1, 1:], ang[1:, 1:]) + d(ang[1:, 1:], ang[1:, :-1]) + d(ang[1:, :-1], ang[:-1, :-1])
    wnd = np.floor(s / (2 * np.pi) + 0.5).astype(int)
    good = ok[:-1, :-1] & ok[:-1, 1:] & ok[1:, 1:] & ok[1:, :-1]
    wnd = np.where(good, wnd, 0)
    return int(wnd.sum()), int((wnd > 0).sum()), int((wnd < 0).sum())


tot, plus, minus = nav_index(2)
case("nav_index().total", tot, "abs:0", P2 + " Poincare-Hopf: indices sum to 1 - M = -2 (phase3 claim 2)")
case("nav_index().plus", plus, "abs:0", P2 + " k = 2 leaves a trap: goal plus extra minima")
tot4, plus4, minus4 = nav_index(4)
case("nav_index(4).plus", plus4, "abs:0", P2 + " k = 4: the goal is the only minimum")
case("nav_index(4).minus", minus4, "abs:0", P2 + " k = 4: exactly M = 3 saddles")
case("nav_index(10, 0.1, 0.55, 0.2, 160).total", nav_index(10, (0.1, 0.55), 0.2, 160)[0], "abs:0", P2 + " total is invariant")
err("nav_index(0.5)", P2 + " k must be >= 1")


# ---- unwinding bench, ported line by line (mulberry32 seed 99)
def mulberry32(seed):
    s = [seed & 0xFFFFFFFF]

    def nxt():
        s[0] = (s[0] + 0x6D2B79F5) & 0xFFFFFFFF
        t = s[0]
        t = ((t ^ (t >> 15)) * (t | 1)) & 0xFFFFFFFF
        t ^= (t + ((t ^ (t >> 7)) * (t | 61))) & 0xFFFFFFFF
        t &= 0xFFFFFFFF
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296
    return nxt


def gauss(r):
    u = 0.0
    while u == 0.0:
        u = r()
    return math.sqrt(-2 * math.log(u)) * math.cos(2 * math.pi * r())


def qmul_js(a, b):
    return [a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3], a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
            a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1], a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]]


def simulate(law, th0, sigma, w0, seed=99):
    ax = [0.3, -0.4, 0.866]
    na = math.sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2])
    axis = [x * (1 / na) for x in ax]
    kp, kd, dt, T, delta = 1, 1.4, 0.01, 25, 0.2
    q = [math.cos(th0 / 2)] + [x * math.sin(th0 / 2) for x in axis]
    w = [x * w0 for x in axis]
    h, switches, travelled, last = 1, 0, 0.0, None
    r = mulberry32(seed)
    k = 0
    while k <= T / dt:
        qm = [x + sigma * gauss(r) for x in q]
        nm = math.sqrt(qm[0] * qm[0] + qm[1] * qm[1] + qm[2] * qm[2] + qm[3] * qm[3])
        qm = [x / nm for x in qm]
        if law == 0:
            s = 1
        elif law == 1:
            s = 1 if qm[0] >= 0 else -1
        else:
            if h * qm[0] < -delta:
                h = -h
            s = h
        if last is not None and s != last:
            switches += 1
        last = s
        tau = [qm[i + 1] * (-kp * s) + w[i] * (-kd) for i in range(3)]
        travelled += math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]) * dt
        w = [w[i] + tau[i] * dt for i in range(3)]
        qd = qmul_js(q, [0.0] + w)
        qn = [q[i] + 0.5 * qd[i] * dt for i in range(4)]
        nq = math.sqrt(qn[0] * qn[0] + qn[1] * qn[1] + qn[2] * qn[2] + qn[3] * qn[3])
        q = [x / nq for x in qn]
        k += 1
    return travelled, switches, 2 * math.acos(min(1, abs(q[0])))


tr0, sw0, fi0 = simulate(0, 340 * D2R, 0.05, 0)
case("unwind().travelled", tr0, "rel:1e-9", P2 + " bench defaults: naive law unwinds ~340 deg for a 20 deg error")
case("unwind(340*pi/180, 2).travelled", simulate(2, 340 * D2R, 0.05, 0)[0], "rel:1e-9", P2 + " hybrid law takes the short way")
case("unwind(180*pi/180, 1, 0.05).switches", simulate(1, math.pi, 0.05, 0)[1], "abs:0", P2 + " sign law chatters near q_w = 0")
case("unwind(180*pi/180, 2, 0.05).switches", simulate(2, math.pi, 0.05, 0)[1], "abs:0", P2 + " hysteresis stops the chattering")
case("unwind(200*pi/180, 0, 0, 0.5).final", simulate(0, 200 * D2R, 0.0, 0.5)[2], "abs:1e-9", P2 + " no noise, initial rate")
err("unwind(1, 3)", P2 + " law must be 0, 1 or 2")

# ------------------------------------------------------------------ ch04 Hilbert
P4 = "ch04-hilbert"
A = np.array([[1, 0], [0, 1], [1, 1.0]])
y = np.array([1, 2, 3.0])
Pp = A @ np.linalg.lstsq(A, y, rcond=None)[0]
case("proj([1, 0; 0, 1; 1, 1], [1, 2, 3])(3)", Pp[2], prov=P4 + " Theorem 1 projection")
case("A = [1, 0; 1, 1; 1, 2]; y = [1; 0; 4]; r = y - proj(A, y); abs(A'*r)(2)", 0.0, "abs:1e-12", P4 + " residual orthogonal to the columns")
err("proj([1, 2; 2, 4; 3, 6], [1, 2, 3])", P4 + " dependent columns")

Ag = np.array([[1, 1], [1, 0], [0, 1.0]])
Qg, Rg = np.linalg.qr(Ag)
Qg = Qg * np.sign(np.diag(Rg))
case("gs([1, 1; 1, 0; 0, 1])(3, 2)", Qg[2, 1], prov=P4 + " Gram-Schmidt = QR with positive diagonal")
case("Q = gs([1, 2, 0; 0, 1, 1; 1, 0, 3; 2, 1, 1]); (Q'*Q)(2, 3)", 0.0, "abs:1e-14", P4 + " orthonormal columns")
err("gs([1, 2; 2, 4])", P4 + " dependent columns")

Aw = np.array([[1, 0.5], [1, 1.5], [1, 3.0], [1, 4.2]])
yw = np.array([1.1, 2.0, 3.9, 5.1])
Rv = np.array([1, 4, 9, 2.0])
xw = np.linalg.solve(Aw.T @ np.diag(1 / Rv) @ Aw, Aw.T @ np.diag(1 / Rv) @ yw)
case(f"wls({mat(Aw)}, {vec(yw)}, {vec(Rv)})(2)", xw[1], prov=P4 + " weighted least squares")
Rf = np.array([[2, 0.5, 0, 0], [0.5, 1, 0.2, 0], [0, 0.2, 1.5, 0.3], [0, 0, 0.3, 1.0]])
xw2 = np.linalg.solve(Aw.T @ np.linalg.solve(Rf, Aw), Aw.T @ np.linalg.solve(Rf, yw))
case(f"wls({mat(Aw)}, {vec(yw)}, {mat(Rf)})(1)", xw2[0], prov=P4 + " WLS with a full covariance")
case(f"wls({mat(Aw)}, {vec(yw)})(1)", np.linalg.lstsq(Aw, yw, rcond=None)[0][0], prov=P4 + " R = I is ordinary least squares")
err(f"wls({mat(Aw)}, {vec(yw)}, [1, 2])", P4 + " R has the wrong size")

case("vander([1, 2, 3], 2)(3, 3)", 9.0, prov=P4 + " monomial design")
xg = np.linspace(-1, 1, 30)
for dd in (5, 8):
    V = np.vander(xg, dd + 1, increasing=True)
    case(f"V = vander(linspace(-1, 1, 30), {dd}); cond(V'*V)", np.linalg.cond(V.T @ V), "rel:1e-6", P4 + f" normal equations, degree {dd}")
    case(f"V = vander(linspace(-1, 1, 30), {dd}); cond(V'*V)/cond(V)^2", np.linalg.cond(V.T @ V) / np.linalg.cond(V) ** 2, "rel:1e-6",
         P4 + " cond(A'A) = cond(A)^2 (phase5 claim 5)")

# Kalman = batch least squares (phase5 claim 1)
rng = np.random.default_rng(5)
q_, r_, P0 = 0.2, 1.0, 4.0
nk = 12
xk = np.cumsum(rng.normal(0, math.sqrt(q_), nk)) + rng.normal(0, 2)
yk = np.round(xk + rng.normal(0, 1, nk), 6)


def batch(ys, k, q, r, P0):
    L, b = np.zeros((k, k)), np.zeros(k)
    L[0, 0] += 1 / P0
    for i in range(1, k):
        L[i, i] += 1 / q; L[i - 1, i - 1] += 1 / q; L[i, i - 1] -= 1 / q; L[i - 1, i] -= 1 / q
    L[np.diag_indices(k)] += 1 / r
    b[:] = ys[:k] / r
    return np.linalg.solve(L, b)


Y = vec(yk)
case(f"kf_rw({Y})(12)", batch(yk, 12, q_, r_, P0)[-1], prov=P4 + " Kalman = last component of batch LS (bench defaults q 0.2, r 1, P0 4)")
case(f"kf_rw({Y})(5)", batch(yk, 5, q_, r_, P0)[-1], prov=P4 + " Kalman = batch at step 5")
case(f"kf_rw({Y}, 0.5, 2, 1)(9)", batch(yk, 9, 0.5, 2, 1)[-1], prov=P4 + " other noise levels")
case(f"smooth_rw({Y})(5)", batch(yk, 12, q_, r_, P0)[4], prov=P4 + " smoother = projection onto all data")
case(f"smooth_rw({Y})(12) - kf_rw({Y})(12)", 0.0, "abs:1e-12", P4 + " smoother and filter agree at the last step")
err(f"kf_rw({Y}, -1)", P4 + " q must be positive")


def kf_ss_np(a, q, r):
    P = 1.0
    for _ in range(5000):
        P = a * a * P * r / (P + r) + q
    return P, P / (P + r), P * r / (P + r)


Pss = kf_ss_np(0.95, 0.1, 1.0)
case("kf_ss().Pp", Pss[2], prov=P4 + " steady Kalman MSE, bench defaults (chapter quotes 0.2410)")
case("kf_ss(0.5, 1, 2).K", kf_ss_np(0.5, 1, 2)[1], prov=P4 + " steady gain by iteration")
err("kf_ss(1.2)", P4 + " |a| >= 1 is not stationary")


def wiener_nc(a, q, r, N=200000):
    om = np.linspace(0, np.pi, N + 1)
    Ss = q / (1 + a * a - 2 * a * np.cos(om))
    f = Ss * r / (Ss + r)
    return np.trapezoid(f, om) / np.pi


case("wiener().nc", wiener_nc(0.95, 0.1, 1.0), prov=P4 + " non-causal Wiener MSE by quadrature (chapter quotes 0.158)")
case("wiener().causal", Pss[2], prov=P4 + " causal Wiener = steady Kalman")
nc2 = wiener_nc(0.5, 0.3, 0.7)
case("wiener(0.5, 0.3, 0.7).price", kf_ss_np(0.5, 0.3, 0.7)[2] / nc2 - 1, prov=P4 + " price of causality")


def S_mag(w, k, p):
    s = 1j * w
    L = k * (s + 1) / ((s - p) * (s + 3) * (s + 4))
    return np.abs(1 / (1 + L))


def bode_np(k, p):
    u = np.linspace(math.log(1e-14), math.log(1e9), 2_000_001)
    wv = np.exp(u)
    I = np.trapezoid(np.log(S_mag(wv, k, p)) * wv, u) + k / 1e9
    wg = np.exp(np.linspace(math.log(1e-3), math.log(1e3), 200001))
    sm = S_mag(wg, k, p)
    i = int(np.argmax(sm))
    lo, hi = wg[max(i - 1, 0)], wg[min(i + 1, len(wg) - 1)]
    for _ in range(200):
        a1, a2 = hi - 0.618033988749895 * (hi - lo), lo + 0.618033988749895 * (hi - lo)
        if S_mag(a1, k, p) > S_mag(a2, k, p):
            hi = a2
        else:
            lo = a1
    stable = bool((np.roots(np.polyadd(np.polymul(np.polymul([1, -p], [1, 3]), [1, 4]), [k, k])).real < 0).all())
    return I, S_mag(0.5 * (lo + hi), k, p), stable


I1, pk1, st1 = bode_np(30, 0.5)
case("bode_int().integral", I1, "abs:1e-6", P4 + " Bode integral by quadrature = pi p = 1.5708 (phase5 claim 4)")
case("bode_int().theorem", math.pi * 0.5, prov=P4 + " pi * sum Re p_i")
case("bode_int().peak", pk1, "rel:1e-7", P4 + " ||S||inf, bench defaults k 30, p 0.5")
I2, pk2, st2 = bode_np(20, -0.5)
case("bode_int(20, -0.5).integral", I2, "abs:1e-6", P4 + " stable plant: integral 0")
case("bode_int(5, 1).stable", st1 and bode_np(5, 1)[2], "abs:0", P4 + " k = 5, p = 1: closed loop unstable")
err("bode_int(-1)", P4 + " gain must be positive")

GRID = np.linspace(-10, 10, 801)


def postmean_np(ys, eps):
    lp = -GRID ** 2 / 8
    for yy in ys:
        d_ = yy - GRID
        lp = lp + np.log((1 - eps) * np.exp(-d_ * d_ / 2) + eps * np.exp(-d_ * d_ / 200) / 10 + 1e-300)
    w_ = np.exp(lp - lp.max())
    vv_ = (1 - eps) + eps * 100
    return (w_ * GRID).sum() / w_.sum(), (sum(ys) / vv_) / (0.25 + len(ys) / vv_)


ys1 = [0.4, -0.2, 9.5, 0.9, 0.1]
pm = postmean_np(ys1, 0.1)
case(f"postmean({vec(ys1)}).post", pm[0], "rel:1e-9", P4 + " posterior mean ignores the outlier")
case(f"postmean({vec(ys1)}).lin", pm[1], prov=P4 + " best linear estimate is dragged by it")
case(f"postmean({vec(ys1)}, 0).post", postmean_np(ys1, 0.0)[0], "rel:1e-9", P4 + " Gaussian noise: posterior mean = projection")
case(f"postmean({vec(ys1)}, 0).post - postmean({vec(ys1)}, 0).lin", 0.0, "abs:1e-9", P4 + " Gaussian noise: posterior mean = projection")
err("postmean([1, 2], 2)", P4 + " eps outside [0, 1]")

# ------------------------------------------------------------------ ch05 measure
P5 = "ch05-measure"
case("bayes(0, 4, 1.5, 1).m", 1.2, prov=P5 + " Gaussian update: m = P y/(P + R)")
case("bayes(0, 4, 1.5, 1).P", 0.8, prov=P5 + " Gaussian update: P R/(P + R)")
K_ = 2 * 2 / (4 * 2 + 0.5)
case("bayes(1, 2, 3, 0.5, 2).K", K_, prov=P5 + " gain P H/(H^2 P + R)")
err("bayes(0, 1, 1, 0)", P5 + " R must be positive")
xm, Pm = np.array([1.0, -0.5]), np.array([[2, 0.3], [0.3, 1.0]])
Hm, Rm_, ym = np.array([[1, 0], [1, 1.0]]), np.array([[0.5, 0.1], [0.1, 0.8]]), np.array([1.4, 0.2])
S_ = Hm @ Pm @ Hm.T + Rm_
Km = Pm @ Hm.T @ np.linalg.inv(S_)
xpost, Ppost = xm + Km @ (ym - Hm @ xm), Pm - Km @ S_ @ Km.T
args = f"{vec(xm)}, {mat(Pm)}, {vec(ym)}, {mat(Hm)}, {mat(Rm_)}"
case(f"bayes_mean({args})(2)", xpost[1], prov=P5 + " vector Gaussian Bayes update")
case(f"bayes_cov({mat(Pm)}, {mat(Hm)}, {mat(Rm_)})(1, 2)", Ppost[0, 1], prov=P5 + " posterior covariance")
case(f"bayes_cov({mat(Pm)}, {mat(Hm)}, {mat(Rm_)})(1, 1)", np.linalg.inv(np.linalg.inv(Pm) + Hm.T @ np.linalg.inv(Rm_) @ Hm)[0, 0], prov=P5 + " information form agrees")
err(f"bayes_cov({mat(Pm)}, [1, 0, 0], 1)", P5 + " H has the wrong width")
case("bayes_mean(0, 4, 1.5, 1, 1)", 1.2, prov=P5 + " 1 x 1 matrices reduce to the scalar update")

case("bk_density(0.5).wedge", math.cos(0.5) / 2, prov=P5 + " meridian via longitude wedges: cos(lat)/2 (phase6 claim 1)")
case("bk_density(1).band", 1 / (2 * math.pi), prov=P5 + " equator via latitude bands: uniform")
case("bk_density(2).wedge", 0.0, "abs:0", P5 + " latitude beyond pi/2")


def is_quad(s, power):
    x = np.linspace(-35, 35, 1_400_001)
    p = np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
    q = np.exp(-x * x / (2 * s * s)) / (s * math.sqrt(2 * math.pi))
    return np.trapezoid(x ** power * p * p / q, x)


case("is_w2(1.5).w2", is_quad(1.5, 0), "rel:1e-9", P5 + " E_q[w^2] by quadrature of p^2/q (phase6 claim 4)")
case("is_w2(3).fvar", is_quad(3, 4) - 1, "rel:1e-9", P5 + " variance of the IS estimate of E[x^2]")
case("is_w2(1).ess", 1.0, prov=P5 + " q = p: no loss")
case("is_w2(0.65).w2", "inf", "abs:0", P5 + " bench default s = 0.65: s^2 <= 1/2, infinite variance")
err("is_w2(0)", P5 + " s must be positive")


def radinv(i, b):
    f, r = 1.0, 0.0
    while i > 0:
        f /= b
        r += f * (i % b)
        i //= b
    return r


case("vdc(6)", 0.375, prov=P5 + " 6 = 110b -> 0.011b")
case("vdc(10, 3)", radinv(10, 3), prov=P5 + " base 3")
err("vdc(1.5)", P5 + " index must be whole")
case("halton(5, 3)(5, 2)", radinv(5, 3), prov=P5 + " Halton point 5, coordinate 2 (base 3)")
case("halton(100, 4)(100, 4)", radinv(100, 7), prov=P5 + " Halton point 100, base 7")


def halton(n, d):
    primes = [2, 3, 5, 7, 11, 13, 17, 19][:d]
    out = np.empty((n, d))
    for j, b in enumerate(primes):
        i = np.arange(1, n + 1); f = np.ones(n); r = np.zeros(n)
        while (i > 0).any():
            f /= b; r += f * (i % b); i //= b
        out[:, j] = r
    return out


H4 = halton(8192, 4)
smooth = lambda X: np.prod(np.pi / 2 * np.sin(np.pi * X), axis=1)
disc = lambda X: (X.sum(axis=1) < X.shape[1] / 2).astype(float)
case("qmc_err()", abs(smooth(H4).mean() - 1), "abs:1e-12", P5 + " bench 4 default: Halton, d = 4, N = 8192, smooth integrand (phase6 halton)")
case("qmc_err(4096, 4, 0)", abs(disc(H4[:4096]).mean() - 0.5), "abs:1e-12", P5 + " discontinuous integrand")
H8 = halton(1024, 8)
case("qmc_err(1024, 8)", abs(smooth(H8).mean() - 1), "abs:1e-12", P5 + " d = 8")
err("qmc_err(100, 4, 2)", P5 + " kind must be 0 or 1")

case("ess([1, 1, 1, 1])", 4.0, prov=P5 + " equal weights: N_eff = N")
wv = np.array([0.5, 0.3, 0.15, 0.05])
case("ess([0.5, 0.3, 0.15, 0.05])", 1 / np.sum(wv ** 2), prov=P5 + " N_eff = 1/sum w^2")
err("ess([-1, 2])", P5 + " negative weight")

E = np.array([1.0, 2.0, 3.0, 0.5])
logZ = math.log(np.exp(-E).sum())
case("logz([1, 2, 3, 0.5])", logZ, prov=P5 + " log Z")
case("logz([1000, 1001])", -1000 + math.log(1 + math.exp(-1)), prov=P5 + " log-sum-exp without underflow")
case("boltz([1, 2, 3, 0.5])(2)", math.exp(-2) / np.exp(-E).sum(), prov=P5 + " Boltzmann weight")
case("E = [1, 2, 3, 0.5]; gibbs_f(boltz(E), E) - logz(E)", 0.0, "abs:1e-13", P5 + " Gibbs principle attained at Boltzmann (phase6 claim 9)")
qx = np.array([0.25, 0.25, 0.25, 0.25])
case("gibbs_f([0.25, 0.25, 0.25, 0.25], [1, 2, 3, 0.5])", -(qx * E).sum() - (qx * np.log(qx)).sum(), prov=P5 + " any other q gives less")
err("gibbs_f([0.5, 0.5], [1, 2, 3])", P5 + " lengths differ")
pk, qk = np.array([0.5, 0.5]), np.array([0.9, 0.1])
case("kl([0.5, 0.5], [0.9, 0.1])", float((pk * np.log(pk / qk)).sum()), prov=P5 + " KL divergence")
case("kl([1, 0], [0.5, 0.5])", math.log(2), prov=P5 + " 0 log 0 = 0")
case("kl([0.5, 0.5], [1, 0])", "inf", "abs:0", P5 + " absolute continuity fails")
case("E = [1, 2, 3, 0.5]; q = [0.4, 0.1, 0.2, 0.3]; logz(E) - gibbs_f(q, E) - kl(q, boltz(E))", 0.0, "abs:1e-13", P5 + " log Z - F(q) = KL(q || Boltzmann)")

# ------------------------------------------------------------------ ch06 optimization
P6 = "ch06-optimization"


def lp_vertex(c, A, b):
    """min c'x s.t. A x <= b, x >= 0 by enumerating vertices of {G x <= h}."""
    c, A, b = np.asarray(c, float), np.atleast_2d(np.asarray(A, float)), np.asarray(b, float)
    n = len(c)
    G = np.vstack([A, -np.eye(n)])
    h = np.concatenate([b, np.zeros(n)])
    best, bx = None, None
    for S in itertools.combinations(range(len(h)), n):
        Gs = G[list(S)]
        if abs(np.linalg.det(Gs)) < 1e-12:
            continue
        x = np.linalg.solve(Gs, h[list(S)])
        if (G @ x <= h + 1e-9).all():
            v = c @ x
            if best is None or v < best - 1e-12:
                best, bx = v, x
    return best, bx


def lp_dual_vertex(c, A, b):
    """dual: max -b'y s.t. A'y >= -c, y >= 0, i.e. min b'y s.t. -A'y <= c, y >= 0."""
    A = np.atleast_2d(np.asarray(A, float))
    v, y = lp_vertex(np.asarray(b, float), -A.T, np.asarray(c, float))
    return -v, y


LPS = [
    ([-1, -1], [[1, 2], [3, 1]], [4, 6], "two constraints, both tight"),
    ([-3, -2, -4], [[1, 1, 2], [2, 0, 3], [2, 1, 3]], [4, 5, 7], "three variables"),
    ([2, 3], [[-1, -1], [1, -2], [-1, 0]], [-1.5, 1, -0.25], "b < 0 needs phase 1"),
]
for c, A_, b_, what in LPS:
    v, x = lp_vertex(c, A_, b_)
    dv, yv = lp_dual_vertex(c, A_, b_)
    assert abs(v - dv) < 1e-9, (v, dv)
    ex = f"{vec(c)}, {mat(A_)}, {vec(b_)}"
    case(f"lpval({ex})", v, prov=P6 + f" LP by vertex enumeration ({what})")
    case(f"lp({ex})(1)", x[0], "abs:1e-12", P6 + " optimal vertex")
    case(f"lp({ex})(2)", x[1], "abs:1e-12", P6 + " optimal vertex")
    case(f"lpdual({ex})(1)", yv[0], "abs:1e-12", P6 + " dual prices by vertex enumeration of the dual LP")
    case(f"-{vec(b_)}*lpdual({ex})", dv, prov=P6 + " strong duality: c'x* = -b'y*")
# shadow price by finite difference
c, A_, b_, _ = LPS[1]
db = np.array(b_, float); db[0] += 1e-3
fd = (lp_vertex(c, A_, db)[0] - lp_vertex(c, A_, b_)[0]) / 1e-3
case(f"-lpdual({vec(c)}, {mat(A_)}, {vec(b_)})(1)", fd, "rel:1e-6", P6 + " multiplier is the shadow price d(c'x*)/db_1")
err("lp([1, 1], [1, 1], [-1])", P6 + " infeasible")
err("lp([-1, 0], [0, 1], [1])", P6 + " unbounded")
err("lp([1, 1, 1], [1, 1], [1])", P6 + " c has the wrong length")

c, A_, b_, _ = LPS[0]
ex = f"{vec(c)}, {mat(A_)}, {vec(b_)}"
case(f"lp_kkt({ex}, lp({ex}), lpdual({ex})).gap", 0.0, "abs:1e-12", P6 + " KKT: zero duality gap at the optimum")
case(f"lp_kkt({ex}, lp({ex}), lpdual({ex})).slack", 0.0, "abs:1e-12", P6 + " KKT: complementary slackness")
xt, yt = np.array([1.0, 1.0]), np.array([1.0, 0.0])
Ac, bc, cc = np.array(A_, float), np.array(b_, float), np.array(c, float)
case(f"lp_kkt({ex}, [1, 1], [1, 0]).gap", cc @ xt + bc @ yt, prov=P6 + " gap c'x + b'y of a non-optimal pair")
case(f"lp_kkt({ex}, [3, 1], [0, 0]).primal", max(0, (Ac @ [3, 1] - bc).max()), prov=P6 + " primal infeasibility")
case(f"lp_kkt({ex}, [1, 1], [0, 0]).dual", max(0, -(cc + Ac.T @ [0, 0]).min()), prov=P6 + " dual infeasibility")

sv = np.array([3, -0.5, 1.2, -2.5])
st = np.sign(sv) * np.maximum(np.abs(sv) - 1, 0)
case("soft([3, -0.5, 1.2, -2.5], 1)(4)", st[3], prov=P6 + " soft threshold")
case("soft([3, -0.5, 1.2, -2.5], 1)(2)", 0.0, "abs:0", P6 + " soft threshold kills small entries")
err("soft([1, 2], -1)", P6 + " threshold must be >= 0")

rng = np.random.default_rng(7)
Am = rng.normal(size=(10, 5))
xt_ = np.zeros(5); xt_[[1, 3]] = [2, -1.5]
bm = np.round(Am @ xt_ + 0.1 * rng.normal(size=10), 6)
Am = np.round(Am, 6)


def ista(A, b, lam, it=60000):
    t = 1 / np.linalg.norm(A, 2) ** 2
    soft_ = lambda v, k: np.sign(v) * np.maximum(np.abs(v) - k, 0)
    x = np.zeros(A.shape[1])
    for _ in range(it):
        x = soft_(x - t * A.T @ (A @ x - b), t * lam)
    return x


xl = ista(Am, bm, 2.0)
assert np.abs(xl - np.sign(xl) * np.maximum(np.abs(xl - (1 / np.linalg.norm(Am, 2) ** 2) * Am.T @ (Am @ xl - bm)) - 2.0 / np.linalg.norm(Am, 2) ** 2, 0)).max() < 1e-12
case(f"lasso({mat(Am)}, {vec(bm)}, 2)(2)", xl[1], "abs:1e-9", P6 + " LASSO minimizer = ISTA fixed point (phase7 claim 7)")
case(f"lasso({mat(Am)}, {vec(bm)}, 2)(4)", xl[3], "abs:1e-9", P6 + " LASSO")
zero_idx = int(np.flatnonzero(xl == 0)[0]) + 1 if (xl == 0).any() else None
if zero_idx:
    case(f"lasso({mat(Am)}, {vec(bm)}, 2)({zero_idx})", 0.0, "abs:0", P6 + " LASSO sets this coefficient exactly to zero")
case(f"lasso({mat(Am)}, {vec(bm)}, 0)(1)", np.linalg.lstsq(Am, bm, rcond=None)[0][0], "abs:1e-9", P6 + " lambda = 0 is least squares")
err("lasso([1, 2; 3, 4], [1, 2, 3], 1)", P6 + " b has the wrong length")

# discrete LQ KKT system (phase7 claim 6)
A_ = np.array([[1, 0.1], [0, 1]]); B_ = np.array([[0.005], [0.1]]); Q_ = np.diag([1, 0.1]); R_ = np.array([[0.5]])
Qf, N, x0 = np.diag([10, 1.0]), 30, np.array([1, 0.0])
nx, nu = 2, 1; nz = (N + 1) * nx + N * nu
ix = lambda k: slice(k * nx, (k + 1) * nx)
iu = lambda k: slice((N + 1) * nx + k * nu, (N + 1) * nx + (k + 1) * nu)
Hk = np.zeros((nz, nz))
for k in range(N):
    Hk[ix(k), ix(k)], Hk[iu(k), iu(k)] = Q_, R_
Hk[ix(N), ix(N)] = Qf
Ek, ek = np.zeros(((N + 1) * nx, nz)), np.zeros((N + 1) * nx)
Ek[:nx, ix(0)] = np.eye(nx); ek[:nx] = x0
for k in range(N):
    r = slice((k + 1) * nx, (k + 2) * nx); Ek[r, ix(k + 1)] = np.eye(nx); Ek[r, ix(k)] = -A_; Ek[r, iu(k)] = -B_
sol = np.linalg.solve(np.block([[Hk, Ek.T], [Ek, np.zeros((len(ek), len(ek)))]]), np.concatenate([np.zeros(nz), ek]))
zz, lam = sol[:nz], sol[nz:]
LQ = "[1, 0.1; 0, 1], [0.005; 0.1], [1, 0; 0, 0.1], 0.5, [10, 0; 0, 1], [1; 0], 30"
case(f"lq_costate({LQ}, 0)(1)", -lam[ix(0)][0], prov=P6 + " KKT multiplier of x_0 = P_0 x_0 (phase7 claim 6, its example)")
case(f"lq_costate({LQ}, 10)(2)", -lam[ix(10)][1], prov=P6 + " costate at step 10 from the KKT solve")
case(f"lq_costate({LQ}, 30)(1)", -lam[ix(30)][0], prov=P6 + " lambda_N = Qf x_N")
case(f"lq_cost({LQ})", 0.5 * zz @ Hk @ zz, prov=P6 + " optimal cost 1/2 z'Hz of the KKT solution")
err("lq_costate([1, 0.1; 0, 1], [0.005, 0.1], [1, 0; 0, 0.1], 0.5, [10, 0; 0, 1], [1; 0], 30)", P6 + " B must be n x m")
err(f"lq_costate({LQ}, 31)", P6 + " k beyond the horizon")
err("lq_cost([1, 0.1; 0, 1], [0.005; 0.1], [1, 0; 0, 0.1], 0.5, [10, 0; 0, 1], [1; 0], 0)", P6 + " horizon must be >= 1")


def fa(a, x):
    return (1 - a) * (x - 1) ** 2 + a * ((x ** 2 - 1) ** 2 - 0.5 * x)


def real_roots(coeffs):
    r = np.roots(coeffs)
    return [z.real for z in r if abs(z.imag) < 1e-7]


def gfun(a, lam_):
    # min_x f + lam x via stationary points of the cubic
    cands = real_roots([4 * a, 0, 2 * (1 - a) - 4 * a, -2 * (1 - a) - a / 2 + lam_]) if a > 0 else [1 - lam_ / 2]
    return min(fa(a, x) + lam_ * x for x in cands)


def gap_np(a):
    cands = [0.0] + ([x for x in real_roots([4 * a, 0, 2 * (1 - a) - 4 * a, -2 * (1 - a) - a / 2]) if x <= 0] if a > 0 else [])
    p = min(fa(a, x) for x in cands)
    lams = np.linspace(0, 20, 20001)
    gs_ = np.array([gfun(a, l) for l in lams])
    i = int(np.argmax(gs_))
    lo, hi = lams[max(i - 1, 0)], lams[min(i + 1, len(lams) - 1)]
    for _ in range(200):
        l1, l2 = hi - 0.618033988749895 * (hi - lo), lo + 0.618033988749895 * (hi - lo)
        if gfun(a, l1) < gfun(a, l2):
            lo = l1
        else:
            hi = l2
    d = max(gs_.max(), gfun(a, 0.5 * (lo + hi)))
    return p, d


for a in (0.0, 0.2, 0.5, 0.6, 0.8, 1.0):
    p, d = gap_np(a)
    if a in (0.0, 0.6, 1.0):
        case(f"dualgap({a}).p", p, "abs:1e-12", P6 + f" primal p* at a = {a}")
        case(f"dualgap({a}).d", d, "abs:1e-9", P6 + f" dual d* at a = {a}")
    case(f"dualgap({a}).gap", p - d, "abs:1e-9", P6 + f" counterexample 1: gap at a = {a} (phase7 claim 3)")
err("dualgap(1.5)", P6 + " a outside [0, 1]")


def auction_np(a, eps):
    n = len(a); price = np.zeros(n); owner = -np.ones(n, int); assign = -np.ones(n, int); free = list(range(n))
    while free:
        i = free.pop(0); v = a[i] - price; j = int(np.argmax(v)); v2 = np.max(np.delete(v, j))
        price[j] += v[j] - v2 + eps
        if owner[j] >= 0:
            assign[owner[j]] = -1; free.append(owner[j])
        owner[j], assign[i] = i, j
    return assign


rng = np.random.default_rng(3)
nA = 5
Pp_, Qp_ = rng.random((nA, 2)) * 100, rng.random((nA, 2)) * 100
aB = np.round(150 - np.linalg.norm(Pp_[:, None] - Qp_[None], axis=2))
opt = max(sum(aB[i, p[i]] for i in range(nA)) for p in itertools.permutations(range(nA)))
asg = auction_np(aB, 0.99 / nA)
assert sum(aB[i, asg[i]] for i in range(nA)) == opt
aBs = mat(aB)
case(f"auction({aBs}, 0.19)(1)", asg[0] + 1, "abs:0", P6 + " bench 3 instance (150 - distance), eps < 1/n: optimal")
case(f"auction({aBs}, 0.19)(4)", asg[3] + 1, "abs:0", P6 + " auction assignment")
val_expr = " + ".join(f"a({i + 1}, s({i + 1}))" for i in range(nA))
case(f"a = {aBs}; s = auction(a, 0.19); {val_expr}", opt, "abs:0", P6 + " Theorem 6: optimal (brute force over 5! assignments)")
asg25 = auction_np(aB, 25.0)
v25 = sum(aB[i, asg25[i]] for i in range(nA))
case(f"a = {aBs}; s = auction(a, 25); {val_expr}", v25, "abs:0", P6 + f" eps = 25: within n eps of the optimum {opt:.0f}")
err("auction([1, 2, 3; 4, 5, 6], 1)", P6 + " benefits must be square")

print("# Atlas chapters 1, 2, 4, 5, 6 (LibAtlas.cpp). Generated by Tests/reference/gen_atlas.py from numpy;")
print("# expected values are independent of the C++ (series, eigen/SVD, vertex enumeration, quadrature, ports of the page benches).")
for r in rows:
    print(r)
