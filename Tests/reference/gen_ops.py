"""Chapter 9 reference cases for LibOps.cpp (run by gen_ops.mjs; prints TSV lines).

Chapter 9 has no JS module, so expected values come from numpy and from the functions of
verification/python/phase10_graphs.py (dijkstra, bellman_ford, lap), which are loaded from that file's
source without running its checks.
"""
import ast
import itertools
import math
import os
import sys

import numpy as np

WII = sys.argv[1] if len(sys.argv) > 1 else "/home/user/where-it-isnt"
SRC = os.path.join(WII, "verification/python/phase10_graphs.py")
tree = ast.parse(open(SRC).read())
ns = {"math": math, "np": np, "heapq": __import__("heapq")}
defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ("dijkstra", "bellman_ford", "lap")]
exec(compile(ast.Module(body=defs, type_ignores=[]), SRC, "exec"), ns)
dijkstra, bellman_ford, lap = ns["dijkstra"], ns["bellman_ford"], ns["lap"]

out = []


def num(x):
    x = float(x)
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return repr(x)


def mat(M):
    return "[" + "; ".join(", ".join(num(x) for x in row) for row in M) + "]"


def vec(v):
    return "[" + ", ".join(num(x) for x in v) + "]"


def emit(expr, expected, prov, tol=None):
    if expected == "error":
        out.append(f"{expr}\terror\tabs:0\t{prov}")
        return
    if tol is None:
        tol = "abs:1e-12" if abs(expected) < 1e-12 else "rel:1e-9"
    out.append(f"{expr}\t{num(expected)}\t{tol}\t{prov}")


def cost_matrix(n, edges):
    W = np.full((n, n), math.inf)
    np.fill_diagonal(W, 0.0)
    for u, v, w in edges:
        W[u, v] = min(W[u, v], w)
    return W


def adj_lists(n, edges):
    adj = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, w))
    return adj


rng = np.random.default_rng(9)

# Dijkstra and Bellman-Ford on random digraphs (phase10_graphs.py claim 1)
P = "phase10_graphs.py dijkstra / bellman_ford"
for trial in range(3):
    n = 7 + trial
    edges = {}
    for u, v in rng.integers(0, n, size=(3 * n, 2)):
        if u != v:
            edges[(int(u), int(v))] = round(float(rng.random() * 10), 3)
    E = [(u, v, w) for (u, v), w in edges.items()]
    W = cost_matrix(n, E)
    d = dijkstra(n, adj_lists(n, E), 0)
    b = bellman_ford(n, E, 0)
    assert np.allclose(d, b)
    for k in (1, n // 2, n - 1):
        emit(f"dijkstra({mat(W)})({k + 1})", d[k], P + f" (n = {n})")
        emit(f"bellman({mat(W)})({k + 1})", b[k], P)
    s = 2
    ds = dijkstra(n, adj_lists(n, E), s)
    emit(f"dijkstra({mat(W)}, {s + 1})({n})", ds[n - 1], P + " from node 3")
    # Floyd-Warshall against Bellman-Ford from every source
    F = np.array([bellman_ford(n, E, i) for i in range(n)])
    emit(f"floyd({mat(W)})({n}, 1)", F[n - 1, 0], P + " all pairs")
    emit(f"floyd({mat(W)})(2, {n - 1})", F[1, n - 2], P + " all pairs")

# the chapter's counterexample: one negative edge
neg = [(0, 1, 2.0), (0, 2, 3.0), (2, 1, -2.0)]
Wn = cost_matrix(3, neg)
emit(f"bellman({mat(Wn)})(2)", bellman_ford(3, neg, 0)[1], "phase10_graphs.py: truth 1, Dijkstra would settle 2")
emit(f"floyd({mat(Wn)})(1, 2)", 1.0, "phase10_graphs.py counterexample")
emit(f"dijkstra({mat(Wn)})", "error", "negative edge refused")
cyc = cost_matrix(3, [(0, 1, 1.0), (1, 2, -3.0), (2, 0, 1.0)])
emit(f"bellman({mat(cyc)})", "error", "negative cycle")
emit(f"floyd({mat(cyc)})", "error", "negative cycle")
emit("dijkstra([0, 1; 1, 0], 3)", "error", "source out of range")
unreach = cost_matrix(3, [(0, 1, 1.0)])
emit(f"dijkstra({mat(unreach)})(3)", math.inf, "unreachable node")

# grid metric (claim 4)
SQ2 = math.sqrt(2)
P = "ch09 octile distance max + (sqrt2 - 1) min"
emit("octile(3, 4)", 4 + (SQ2 - 1) * 3, P)
emit("octile(-5, 2)", 5 + (SQ2 - 1) * 2, P)
emit("octile(7, 7)", 7 * SQ2, P + " (diagonal)")
angles = np.linspace(0, math.pi / 4, 2001)
ratio = np.maximum(np.cos(angles), np.sin(angles)) + (SQ2 - 1) * np.minimum(np.cos(angles), np.sin(angles))
assert abs(ratio.max() - 1 / math.cos(math.pi / 8)) < 1e-6
emit("gridbound()", 1 / math.cos(math.pi / 8), "phase10_graphs.py claim 4: max ratio = 1/cos(pi/8) = 1.0824")
emit("gridbound", 1 / math.cos(math.pi / 8), "phase10_graphs.py claim 4 (no parentheses)")
emit("gridratio(pi/8)", 1 / math.cos(math.pi / 8), "phase10_graphs.py claim 4: attained at 22.5 deg")
emit("gridratio(0)", 1.0, "axis-aligned: exact")
emit("gridratio(0.3)", float(np.maximum(math.cos(0.3), math.sin(0.3)) + (SQ2 - 1) * min(math.cos(0.3), math.sin(0.3))), P)


# consensus (claim 5)
def adjacency(n, E, w=None):
    A = np.zeros((n, n))
    for k, (u, v) in enumerate(E):
        A[u, v] = A[v, u] = 1.0 if w is None else w[k]
    return A


n = 12
ring = [(i, (i + 1) % n) for i in range(n)]
star = [(0, i) for i in range(1, n)]
comp = list(itertools.combinations(range(n), 2))
for name, E, lam2 in (("ring", ring, 2 - 2 * math.cos(2 * math.pi / n)), ("star", star, 1.0), ("complete", comp, float(n))):
    A = adjacency(n, E)
    ev = np.sort(np.linalg.eigvalsh(lap(n, E)))
    assert abs(ev[1] - lam2) < 1e-9
    emit(f"lambda2({mat(A)})", lam2, f"phase10_graphs.py claim 5: {name} of 12, closed form")
Ad = adjacency(12, [(i, i + 1) for i in range(5)] + [(i, i + 1) for i in range(6, 11)])
emit(f"lambda2({mat(Ad)})", 0.0, "phase10_graphs.py: two components, lambda2 = 0", "abs:1e-9")
Aw = np.round(rng.random((6, 6)) * (rng.random((6, 6)) < 0.6), 3)
Aw = np.triu(Aw, 1)
Aw = Aw + Aw.T
Aw[0, 1] = Aw[1, 0] = 0.5
for i in range(5):                                        # keep it connected: a weighted path underneath
    Aw[i, i + 1] = Aw[i + 1, i] = max(Aw[i, i + 1], 0.2)
Lw = np.diag(Aw.sum(1)) - Aw
emit(f"lambda2({mat(Aw)})", float(np.sort(np.linalg.eigvalsh(Lw))[1]), "numpy eigvalsh of the weighted Laplacian")
emit(f"lap({mat(Aw)})(2, 2)", float(Lw[1, 1]), "L = D - A")
emit(f"lap({mat(Aw)})(1, 2)", float(Lw[0, 1]), "L = D - A")
emit("lambda2([0, 1; 2, 0])", "error", "not symmetric")

x0 = np.round(rng.normal(size=6), 4)
for t in (0.5, 3.0):
    w, V = np.linalg.eigh(Lw)
    xt = V @ np.diag(np.exp(-w * t)) @ V.T @ x0
    emit(f"consensus({mat(Aw)}, {vec(x0)}, {t}).dis", float(np.linalg.norm(xt - xt.mean())), "phase10_graphs.py decay_rate: x(t) = V exp(-w t) V' x0")
    emit(f"consensus({mat(Aw)}, {vec(x0)}, {t}).bound", float(math.exp(-w[1] * t) * np.linalg.norm(x0 - x0.mean())), "ch09 Thm 4 bound exp(-lambda2 t)")
emit(f"consensus({mat(Aw)}, {vec(x0)}, 2).mean", float(x0.mean()), "ch09 Thm 4: the average is preserved")
emit(f"consensus({mat(Aw)}, [1, 2], 1)", "error", "x0 length")

# pose graph and effective resistance (claim 6)
N, s2 = 20, 0.04
chain = adjacency(N, [(i, i + 1) for i in range(N - 1)], [1 / s2] * (N - 1))
loop = adjacency(N, [(i, i + 1) for i in range(N - 1)] + [(N - 1, 0)], [1 / s2] * N)
emit(f"reff({mat(chain)}, 20, 1)", s2 * 19, "phase10_graphs.py claim 6: series law k s^2 (ch09: 0.760)")
emit(f"reff({mat(loop)}, 11, 1)", s2 * 10 * (N - 10) / N, "phase10_graphs.py claim 6: parallel law s^2 k(N-k)/N (ch09: 0.200)")
emit(f"reff({mat(loop)}, 5, 1)", s2 * 4 * (N - 4) / N, "phase10_graphs.py claim 6: parallel law")
Lp = np.linalg.pinv(Lw)
for u, v in ((0, 5), (2, 3)):
    e = np.zeros(6)
    e[u], e[v] = 1, -1
    emit(f"reff({mat(Aw)}, {u + 1}, {v + 1})", float(e @ Lp @ e), "ch09 Def: (e_u - e_v)' L^+ (e_u - e_v), numpy pinv")
emit(f"reff({mat(Ad)}, 1, 12)", "error", "disconnected")

# matrix-tree theorem (claim 7)
for trial in range(2):
    nn = 6
    E = [e for e in itertools.combinations(range(nn), 2) if rng.random() < 0.6]
    cnt = 0
    for sub in itertools.combinations(E, nn - 1):
        parent = list(range(nn))

        def find(a):
            while parent[a] != a:
                a = parent[a]
            return a
        ok = True
        for u, v in sub:
            ru, rv = find(u), find(v)
            if ru == rv:
                ok = False
                break
            parent[ru] = rv
        cnt += ok
    emit(f"spantrees({mat(adjacency(nn, E))})", float(cnt), "phase10_graphs.py claim 7: spanning trees by brute force")
emit(f"spantrees({mat(adjacency(5, list(itertools.combinations(range(5), 2))))})", 125.0, "Cayley n^(n-2) for K5")
emit(f"spantrees({mat(adjacency(12, ring))})", 12.0, "a cycle of 12 has 12 spanning trees")
emit(f"spantrees({mat(Ad)})", 0.0, "disconnected: no spanning tree", "abs:1e-9")

print("\n".join(out))
