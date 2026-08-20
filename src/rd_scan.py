"""The 3-D Sylvester 5-point extremal question (open analogue of Blaschke's theorem):
is  P_5(K) = P(5 uniform points in a convex body K in R^3 are in convex position)
minimised by the simplex and maximised by the ball?

Five points in R^3 are in convex position iff none lies in the tetrahedron of the other four,
and at most one can, so the five events are disjoint and

    P_5(K) = 1 - 5 E[vol conv(P_1..P_4)] / vol(K) = 1 - 5 E[A_4].

E[vol T_4] is a plain 4-point determinant expectation -- no hull code, no rejection -- so it
is cheap and accurate; vol(K) is supplied exactly with each body.  P_5 is affine invariant.

Exact anchors
    simplex   E[A_4] = 13/720 - pi^2/15015      (Buchta-Reitzner 1992)  P_5 = 0.91300875...
    cube      E[A_4] = 3977/216000 - pi^2/2160  (Zinani 2003)           P_5 = 0.93078629...
    ball      E[A_4] = 9/715                                            P_5 = 134/143 = 0.93706294
"""
from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np
from scipy.spatial import ConvexHull, Delaunay

EA4_SIMPLEX = 13 / 720 - math.pi**2 / 15015
EA4_CUBE = 3977 / 216000 - math.pi**2 / 2160
EA4_BALL = 9 / 715
P5_SIMPLEX = 1 - 5 * EA4_SIMPLEX
P5_CUBE = 1 - 5 * EA4_CUBE
P5_BALL = 1 - 5 * EA4_BALL
TAU = 2 * math.pi


# ---------------------------------------------------------------- samplers (each: (fn, vol))
def _det3(M):
    a, b, c = M[:, 0, :], M[:, 1, :], M[:, 2, :]
    return (a[:, 0] * (b[:, 1] * c[:, 2] - b[:, 2] * c[:, 1])
            - a[:, 1] * (b[:, 0] * c[:, 2] - b[:, 2] * c[:, 0])
            + a[:, 2] * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0]))


def _tet_uniform(rng, T, cw, m):
    idx = np.searchsorted(cw, rng.random(m))
    u = np.sort(rng.random((m, 3)), axis=1)
    b = np.empty((m, 4))
    b[:, 0] = u[:, 0]
    b[:, 1] = u[:, 1] - u[:, 0]
    b[:, 2] = u[:, 2] - u[:, 1]
    b[:, 3] = 1 - u[:, 2]
    return np.einsum("mk,mkd->md", b, T[idx])


def polytope(V):
    """Uniform sampler + exact volume for conv(V)."""
    V = np.asarray(V, float)
    h = ConvexHull(V)
    V = V[h.vertices]
    tri = Delaunay(V)
    T = V[tri.simplices]
    w = np.abs(np.linalg.det(T[:, 1:, :] - T[:, :1, :])) / 6
    keep = w > 1e-13 * w.max()
    T, w = T[keep], w[keep]
    vol = float(w.sum())
    assert abs(vol - h.volume) < 1e-9 * max(vol, 1), (vol, h.volume)
    cw = np.cumsum(w) / vol
    cw[-1] = 1.0
    return (lambda rng, m: _tet_uniform(rng, T, cw, m)), vol, len(V)


def ball(rng, m):
    g = rng.standard_normal((m, 3))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    return g * rng.random((m, 1)) ** (1 / 3)


def ellipsoid(abc):
    a = np.asarray(abc, float)
    return (lambda rng, m: ball(rng, m) * a), 4 * math.pi / 3 * float(a.prod()), 0


def cylinder(h):
    def s(rng, m):
        t = rng.random((m, 1)) * TAU
        r = np.sqrt(rng.random((m, 1)))
        return np.hstack([r * np.cos(t), r * np.sin(t), h * rng.random((m, 1))])
    return s, math.pi * h, 0


def cone(h):
    """Apex at z = h, unit base disk at z = 0.  Height density prop. to (1-z/h)^2."""
    def s(rng, m):
        w = rng.random((m, 1)) ** (1 / 3)          # w = 1 - z/h
        t = rng.random((m, 1)) * TAU
        r = w * np.sqrt(rng.random((m, 1)))
        return np.hstack([r * np.cos(t), r * np.sin(t), h * (1 - w)])
    return s, math.pi * h / 3, 0


def bicone(h):
    cs, _, _ = cone(h)

    def s(rng, m):
        P = cs(rng, m)
        f = rng.random(m) < 0.5
        P[f, 2] *= -1
        return P
    return s, 2 * math.pi * h / 3, 0


def capsule(L):
    """Spherocylinder: radius-1 cylinder of length L capped by two hemispheres."""
    vc = math.pi * L
    vs = 4 * math.pi / 3

    def s(rng, m):
        pick = rng.random(m) < vc / (vc + vs)
        out = np.empty((m, 3))
        nc = int(pick.sum())
        if nc:
            t = rng.random((nc, 1)) * TAU
            r = np.sqrt(rng.random((nc, 1)))
            out[pick] = np.hstack([r * np.cos(t), r * np.sin(t), (rng.random((nc, 1)) - .5) * L])
        ns = m - nc
        if ns:
            B = ball(rng, ns)
            B[:, 2] += np.where(B[:, 2] >= 0, L / 2, -L / 2)
            out[~pick] = B
        return out
    return s, vc + vs, 0


def ballcap(hh):
    """Cap z >= 1 - hh of the unit ball; hh = 1 is the half-ball, hh = 2 the ball."""
    acc = hh * hh * (3 - hh) / 4          # cap volume / ball volume = acceptance rate

    def s(rng, m):
        out = np.empty((m, 3))
        got = 0
        while got < m:
            want = min(m - got, 200_000)
            B = ball(rng, min(4_000_000, max(64, int(2.0 * want / acc))))
            B = B[B[:, 2] >= 1 - hh]
            k = min(len(B), m - got)
            out[got:got + k] = B[:k]
            got += k
        return out
    return s, math.pi * hh * hh * (3 - hh) / 3, 0


# ---------------------------------------------------------------- polytope vertex sets
def reg_simplex():
    return np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)


def cube_v():
    return np.array([[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)], float)


def octahedron():
    return np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)


def prism(m, h=1.0):
    a = TAU * np.arange(m) / m
    return np.vstack([np.stack([np.cos(a), np.sin(a), np.zeros(m)], -1),
                      np.stack([np.cos(a), np.sin(a), np.full(m, h)], -1)])


def antiprism(m, h=1.0):
    a = TAU * np.arange(m) / m
    b = a + math.pi / m
    return np.vstack([np.stack([np.cos(a), np.sin(a), np.zeros(m)], -1),
                      np.stack([np.cos(b), np.sin(b), np.full(m, h)], -1)])


def bipyramid(m, h=1.0):
    a = TAU * np.arange(m) / m
    return np.vstack([np.stack([np.cos(a), np.sin(a), np.zeros(m)], -1), [[0, 0, h], [0, 0, -h]]])


def pyramid(m, h=1.0):
    a = TAU * np.arange(m) / m
    return np.vstack([np.stack([np.cos(a), np.sin(a), np.zeros(m)], -1), [[0, 0, h]]])


def dodecahedron():
    p = (1 + math.sqrt(5)) / 2
    V = [[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    V += [[0, s1 / p, s2 * p] for s1 in (-1, 1) for s2 in (-1, 1)]
    V += [[s1 / p, s2 * p, 0] for s1 in (-1, 1) for s2 in (-1, 1)]
    V += [[s1 * p, 0, s2 / p] for s1 in (-1, 1) for s2 in (-1, 1)]
    return np.array(V, float)


def icosahedron():
    p = (1 + math.sqrt(5)) / 2
    V = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            V += [[0, s1, s2 * p], [s1, s2 * p, 0], [s2 * p, 0, s1]]
    return np.array(V, float)


def sphere_points(n, rng):
    g = rng.standard_normal((n, 3))
    return g / np.linalg.norm(g, axis=1, keepdims=True)


# ---------------------------------------------------------------- catalogue
def catalogue(seed=99):
    rng = np.random.default_rng(seed)
    B = []

    def addP(name, fam, V):
        try:
            s, v, nv = polytope(V)
        except Exception as e:
            print(f"  skip {name}: {e}")
            return
        B.append((name, fam, s, v, nv))

    def addS(name, fam, t):
        s, v, nv = t
        B.append((name, fam, s, v, nv))

    addP("simplex(regular)", "anchor", reg_simplex())
    addP("simplex(right)", "anchor", np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float))
    addP("simplex(sheared)", "anchor",
         np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
         @ np.array([[1, .7, -.4], [0, 1.3, .5], [0, 0, .6]]))
    addP("cube", "anchor", cube_v())
    addP("box-1:2:5", "anchor", cube_v() * np.array([1, 2, 5.]))
    addS("ball", "anchor", (ball, 4 * math.pi / 3, 0))
    addS("ellipsoid-1:2:4", "anchor", ellipsoid([1, 2, 4]))
    addP("octahedron", "platonic", octahedron())
    addP("dodecahedron", "platonic", dodecahedron())
    addP("icosahedron", "platonic", icosahedron())
    for m in (3, 4, 5, 6, 8, 12, 20, 40):
        for h in (0.3, 0.7, 1.2, 2.0, 4.0):
            addP(f"prism{m}-h{h}", "prism", prism(m, h))
        addP(f"antiprism{m}", "prism", antiprism(m, 1.0))
        for h in (0.5, 1.0, 2.0):
            addP(f"bipyr{m}-h{h}", "bipyramid", bipyramid(m, h))
        for h in (0.5, 1.0, 2.0, 4.0):
            addP(f"pyramid{m}-h{h}", "pyramid", pyramid(m, h))
    for h in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        addS(f"cylinder-h{h}", "smooth", cylinder(h))
        addS(f"cone-h{h}", "smooth", cone(h))
        addS(f"bicone-h{h}", "smooth", bicone(h))
    for L in (0.25, 0.5, 1.0, 2.0, 4.0):
        addS(f"capsule-L{L}", "smooth", capsule(L))
    for hh in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 1.9):
        addS(f"ballcap-{hh}", "smooth", ballcap(hh))
    for f in (0.4, 0.6, 0.8, 1.25, 1.6, 2.5, 4.0):
        V = reg_simplex().copy(); V[0] *= f
        addP(f"simplex-push{f}", "near-simplex", V)
    for t in (0.02, 0.05, 0.1, 0.2, 0.35):
        S = reg_simplex()
        addP(f"simplex-trunc4-{t}", "near-simplex",
             np.array([S[i] + t * (S[j] - S[i]) for i in range(4) for j in range(4) if i != j]))
        addP(f"simplex-trunc1-{t}", "near-simplex",
             np.array([S[1], S[2], S[3]] + [S[0] + t * (S[j] - S[0]) for j in (1, 2, 3)]))
    for n in (6, 8, 10, 12, 16, 24, 32, 64, 128, 256):
        for i in range(3):
            addP(f"randsphere-{n}-{i}", "near-ball", sphere_points(n, rng))
    for n in (5, 6, 8, 10, 14, 20, 30, 50):
        for i in range(3):
            addP(f"randgauss-{n}-{i}", "random", rng.standard_normal((n, 3)))
            addP(f"randball-{n}-{i}", "random", ball(rng, n))
    return B


# ---------------------------------------------------------------- MC
def mc_p5(sampler, vol, nsamp, seed, batch=1_000_000):
    rng = np.random.default_rng(seed)
    s = 0.0
    ss = 0.0
    n = 0
    while n < nsamp:
        m = min(batch, nsamp - n)
        P = sampler(rng, 4 * m).reshape(m, 4, 3)
        v = np.abs(_det3(P[:, 1:, :] - P[:, :1, :])) / 6.0
        s += float(v.sum())
        ss += float((v * v).sum())
        n += m
    mean = s / n
    se = math.sqrt(max(ss / n - mean * mean, 0.0) / n)
    ea4 = mean / vol
    return dict(EA4=ea4, se_EA4=se / vol, P5=1 - 5 * ea4, se_P5=5 * se / vol, nsamp=n)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=float, default=2e7)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--out", default="../results/rd_scan.json")
    ap.add_argument("--resume", default=None, help="json with already-computed rows to reuse")
    a = ap.parse_args()
    B = catalogue()
    done = {}
    if a.resume:
        try:
            done = {r["body"]: r for r in json.load(open(a.resume))["rows"]}
            print(f"resuming: {len(done)} bodies already computed")
        except Exception as e:
            print("resume failed:", e)
    print(f"{len(B)} 3-D bodies, {a.samples:.3g} tetrahedra each")
    print(f"anchors: simplex {P5_SIMPLEX:.8f}  cube {P5_CUBE:.8f}  ball {P5_BALL:.8f}\n")
    rows = []
    t0 = time.time()
    EX = {"simplex(regular)": P5_SIMPLEX, "simplex(right)": P5_SIMPLEX,
          "simplex(sheared)": P5_SIMPLEX, "cube": P5_CUBE, "box-1:2:5": P5_CUBE,
          "ball": P5_BALL, "ellipsoid-1:2:4": P5_BALL}
    for i, (name, fam, s, vol, nv) in enumerate(B):
        if name in done:
            rows.append(done[name])
            continue
        r = mc_p5(s, vol, int(a.samples), a.seed + i)
        e = EX.get(name)
        r.update(body=name, family=fam, vol=vol, nv=nv, exact_P5=e,
                 z=((r["P5"] - e) / r["se_P5"]) if e else None,
                 below_simplex=r["P5"] < P5_SIMPLEX - 3 * r["se_P5"],
                 above_ball=r["P5"] > P5_BALL + 3 * r["se_P5"])
        rows.append(r)
        flag = ""
        if r["below_simplex"]:
            flag = "  *** BELOW SIMPLEX ***"
        if r["above_ball"]:
            flag = "  *** ABOVE BALL ***"
        if e is not None or flag or i % 25 == 0:
            print(f"[{i+1:4d}/{len(B)}] {name:22s} P_5 = {r['P5']:.8f} +- {r['se_P5']:.1e}"
                  + (f"  exact {e:.8f} z={r['z']:+.2f}" if e else "") + flag, flush=True)
        if i % 25 == 0:
            json.dump(dict(anchors=dict(simplex=P5_SIMPLEX, cube=P5_CUBE, ball=P5_BALL),
                           samples=int(a.samples), rows=rows), open(a.out, "w"), indent=1)
    json.dump(dict(anchors=dict(simplex=P5_SIMPLEX, cube=P5_CUBE, ball=P5_BALL),
                   samples=int(a.samples), rows=rows), open(a.out, "w"), indent=1)
    rows.sort(key=lambda r: r["P5"])
    print("\n--- 12 lowest P_5 ---")
    for r in rows[:12]:
        print(f"  {r['body']:22s} {r['P5']:.8f} +-{r['se_P5']:.1e}  "
              f"simplex{r['P5']-P5_SIMPLEX:+.7f}")
    print("--- 12 highest P_5 ---")
    for r in rows[-12:]:
        print(f"  {r['body']:22s} {r['P5']:.8f} +-{r['se_P5']:.1e}  ball{r['P5']-P5_BALL:+.7f}")
    print(f"\ntotal {time.time()-t0:.0f}s -> {a.out}")
