"""Monte Carlo estimation of P_K(n) = P(n iid uniform points in K are in convex position).

Two STRUCTURALLY INDEPENDENT testers are provided and must agree on every run:

  * tester A  (`in_convex_position_triangles`)  -- pure numpy: a point set is in
    convex position iff no point lies in the (closed) triangle spanned by three
    others.  Exact for n <= ~12 (n * C(n-1,3) orientation triples per sample).
  * tester B  (`in_convex_position_hull_count`) -- numba: Andrew monotone chain
    convex hull per sample; in convex position iff hull has n vertices.
    Falls back to a pure-python loop when numba is unavailable (slow, for tests).

Samplers: unit square, unit disk, triangle (0,0),(1,0),(0,1) [any triangle by
affine invariance], regular k-gon, unit cube (3D), unit ball (3D), simplex (3D).
Convex position is affine invariant, so parallelogram == square, ellipse == disk.

Usage:
  python convex_position.py --body disk --n 5 --samples 2e7 --seed 1
  python convex_position.py --validate            # reproduce all exact anchors
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from itertools import combinations

import numpy as np

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

# ----------------------------------------------------------------------------
# Samplers: return array of shape (m, n, d)
# ----------------------------------------------------------------------------

def sample_square(rng, m, n):
    return rng.random((m, n, 2))


def sample_disk(rng, m, n):
    r = np.sqrt(rng.random((m, n)))
    t = rng.random((m, n)) * (2 * np.pi)
    return np.stack([r * np.cos(t), r * np.sin(t)], axis=-1)


def sample_triangle(rng, m, n):
    u = rng.random((m, n))
    v = rng.random((m, n))
    flip = u + v > 1
    u = np.where(flip, 1 - u, u)
    v = np.where(flip, 1 - v, v)
    return np.stack([u, v], axis=-1)


def sample_regular_polygon(k):
    """Uniform in the regular k-gon: pick one of k triangles (center, v_i, v_{i+1})
    uniformly, then a uniform point in it."""
    ang = 2 * np.pi * np.arange(k + 1) / k
    V = np.stack([np.cos(ang), np.sin(ang)], axis=-1)

    def _s(rng, m, n):
        idx = rng.integers(0, k, size=(m, n))
        u = rng.random((m, n))
        v = rng.random((m, n))
        flip = u + v > 1
        u = np.where(flip, 1 - u, u)
        v = np.where(flip, 1 - v, v)
        A = V[idx]
        B = V[idx + 1]
        return u[..., None] * A + v[..., None] * B
    _s.__name__ = f"sample_regular_{k}gon"
    return _s


def sample_cube(rng, m, n):
    return rng.random((m, n, 3))


def sample_ball(rng, m, n):
    g = rng.standard_normal((m, n, 3))
    g /= np.linalg.norm(g, axis=-1, keepdims=True)
    r = np.cbrt(rng.random((m, n)))
    return g * r[..., None]


def sample_simplex3(rng, m, n):
    """Uniform in the standard 3-simplex conv{0,e1,e2,e3}: sorted-uniform spacings."""
    u = np.sort(rng.random((m, n, 3)), axis=-1)
    x = u[..., 0]
    y = u[..., 1] - u[..., 0]
    z = u[..., 2] - u[..., 1]
    return np.stack([x, y, z], axis=-1)


BODIES = {
    "square": (sample_square, 2),
    "disk": (sample_disk, 2),
    "triangle": (sample_triangle, 2),
    "pentagon": (sample_regular_polygon(5), 2),
    "hexagon": (sample_regular_polygon(6), 2),
    "octagon": (sample_regular_polygon(8), 2),
    "cube": (sample_cube, 3),
    "ball": (sample_ball, 3),
    "simplex": (sample_simplex3, 3),
}

# ----------------------------------------------------------------------------
# Tester A: triangle / tetrahedron containment (numpy, exact combinatorics)
# ----------------------------------------------------------------------------

def _orient2(a, b, c):
    return (b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1]) - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0])


def in_convex_position_triangles(P):
    """P: (m, n, 2). Returns bool (m,). True iff no point is inside a triangle of others."""
    m, n, _ = P.shape
    ok = np.ones(m, dtype=bool)
    idx = range(n)
    for i in idx:
        others = [j for j in idx if j != i]
        pi_ = P[:, i, :]
        for (a, b, c) in combinations(others, 3):
            A, B, C = P[:, a, :], P[:, b, :], P[:, c, :]
            d1 = _orient2(A, B, pi_)
            d2 = _orient2(B, C, pi_)
            d3 = _orient2(C, A, pi_)
            inside = ((d1 >= 0) & (d2 >= 0) & (d3 >= 0)) | ((d1 <= 0) & (d2 <= 0) & (d3 <= 0))
            ok &= ~inside
    return ok


def _orient3(a, b, c, d):
    ab = b - a
    ac = c - a
    ad = d - a
    return (ab[..., 0] * (ac[..., 1] * ad[..., 2] - ac[..., 2] * ad[..., 1])
            - ab[..., 1] * (ac[..., 0] * ad[..., 2] - ac[..., 2] * ad[..., 0])
            + ab[..., 2] * (ac[..., 0] * ad[..., 1] - ac[..., 1] * ad[..., 0]))


def in_convex_position_tetrahedra(P):
    """P: (m, n, 3). True iff no point lies inside a tetrahedron of four others."""
    m, n, _ = P.shape
    ok = np.ones(m, dtype=bool)
    idx = range(n)
    for i in idx:
        others = [j for j in idx if j != i]
        X = P[:, i, :]
        for (a, b, c, d) in combinations(others, 4):
            A, B, C, D = P[:, a, :], P[:, b, :], P[:, c, :], P[:, d, :]
            s0 = _orient3(A, B, C, D)
            s1 = _orient3(X, B, C, D)
            s2 = _orient3(A, X, C, D)
            s3 = _orient3(A, B, X, D)
            s4 = _orient3(A, B, C, X)
            same = ((s1 * s0 >= 0) & (s2 * s0 >= 0) & (s3 * s0 >= 0) & (s4 * s0 >= 0))
            ok &= ~same
    return ok

# ----------------------------------------------------------------------------
# Tester B: convex hull vertex count (numba monotone chain, 2D)
# ----------------------------------------------------------------------------

if HAVE_NUMBA:
    @njit(cache=True)
    def _hull_size_2d(x, y, n, ox, oy):
        # sort indices by (x, y) -- insertion sort, n is small
        for i in range(n):
            ox[i] = i
        for i in range(1, n):
            k = ox[i]
            j = i - 1
            while j >= 0 and (x[ox[j]] > x[k] or (x[ox[j]] == x[k] and y[ox[j]] > y[k])):
                ox[j + 1] = ox[j]
                j -= 1
            ox[j + 1] = k
        # lower hull
        h = 0
        for ii in range(n):
            i = ox[ii]
            while h >= 2:
                a = oy[h - 2]
                b = oy[h - 1]
                cr = (x[b] - x[a]) * (y[i] - y[a]) - (y[b] - y[a]) * (x[i] - x[a])
                if cr <= 0:
                    h -= 1
                else:
                    break
            oy[h] = i
            h += 1
        lower = h
        # upper hull
        for ii in range(n - 2, -1, -1):
            i = ox[ii]
            while h >= lower + 1:
                a = oy[h - 2]
                b = oy[h - 1]
                cr = (x[b] - x[a]) * (y[i] - y[a]) - (y[b] - y[a]) * (x[i] - x[a])
                if cr <= 0:
                    h -= 1
                else:
                    break
            oy[h] = i
            h += 1
        return h - 1  # last point == first point

    @njit(parallel=True, cache=True)
    def _count_convex_2d(P):
        m, n, _ = P.shape
        cnt = 0
        for s in prange(m):
            ox = np.empty(n, dtype=np.int64)
            oy = np.empty(2 * n + 2, dtype=np.int64)
            hs = _hull_size_2d(P[s, :, 0].copy(), P[s, :, 1].copy(), n, ox, oy)
            if hs == n:
                cnt += 1
        return cnt

    def in_convex_position_hull_count(P):
        return int(_count_convex_2d(np.ascontiguousarray(P)))
else:
    def in_convex_position_hull_count(P):  # pure python fallback (slow)
        m, n, _ = P.shape
        cnt = 0

        def cross(o, a, b):
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

        for s in range(m):
            pts = sorted((float(P[s, i, 0]), float(P[s, i, 1])) for i in range(n))
            lower = []
            for p in pts:
                while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
                    lower.pop()
                lower.append(p)
            upper = []
            for p in reversed(pts):
                while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
                    upper.pop()
                upper.append(p)
            if len(lower) + len(upper) - 2 == n:
                cnt += 1
        return cnt

# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

def estimate(body: str, n: int, samples: int, seed: int, batch: int = 200_000, both: bool = True):
    sampler, dim = BODIES[body]
    rng = np.random.default_rng(seed)
    done = 0
    hitA = 0
    hitB = 0
    t0 = time.time()
    while done < samples:
        m = min(batch, samples - done)
        P = sampler(rng, m, n)
        if dim == 2:
            hitB += in_convex_position_hull_count(P)
            if both:
                hitA += int(in_convex_position_triangles(P).sum())
        else:
            hitA += int(in_convex_position_tetrahedra(P).sum())
            hitB = hitA
        done += m
    p = hitB / samples
    se = math.sqrt(max(p * (1 - p), 1e-300) / samples)
    out = {
        "body": body, "n": n, "samples": samples, "seed": seed,
        "p_hat": p, "std_err": se, "hits": hitB, "seconds": time.time() - t0,
        "dim": dim,
    }
    if dim == 2 and both:
        out["p_hat_testerA"] = hitA / samples
        out["testers_agree"] = (hitA == hitB)
    return out


def validate(samples=2_000_000, seed=12345):
    """Reproduce every exact anchor to within a few sigma using BOTH testers."""
    from anchors import ANCHORS
    plan = [
        ("square", 4, "square_n4"), ("square", 5, "square_n5"), ("square", 6, "square_n6"),
        ("triangle", 4, "triangle_n4"), ("triangle", 5, "triangle_n5"), ("triangle", 6, "triangle_n6"),
        ("disk", 4, "disk_n4"),
        ("cube", 5, "cube3d_n5"), ("ball", 5, "ball3d_n5"),
    ]
    rows = []
    allok = True
    for body, n, key in plan:
        desc, exact = ANCHORS[key]
        r = estimate(body, n, samples, seed)
        z = (r["p_hat"] - exact) / r["std_err"]
        agree = r.get("testers_agree", True)
        ok = abs(z) < 4 and agree
        allok &= ok
        rows.append({**r, "exact": exact, "z": z, "anchor": key, "ok": ok})
        print(f"{key:12s} n={n} exact={exact:.9f} mc={r['p_hat']:.9f} +- {r['std_err']:.2e}"
              f"  z={z:+.2f}  testersAgree={agree}  {'OK' if ok else 'FAIL'}  ({r['seconds']:.1f}s)")
    print("ALL ANCHORS OK" if allok else "ANCHOR FAILURE")
    return rows, allok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--body", default="disk", choices=sorted(BODIES))
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--samples", type=float, default=1e6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--batch", type=int, default=200_000)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--validate-samples", type=float, default=2e6)
    ap.add_argument("--no-both", action="store_true", help="skip tester A (fast production runs)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    print(f"numba: {HAVE_NUMBA}", file=sys.stderr)
    if a.validate:
        rows, ok = validate(int(a.validate_samples), a.seed)
        if a.out:
            json.dump({"rows": rows, "all_ok": ok}, open(a.out, "w"), indent=1)
        sys.exit(0 if ok else 1)
    r = estimate(a.body, a.n, int(a.samples), a.seed, a.batch, both=not a.no_both)
    print(json.dumps(r, indent=1))
    if a.out:
        json.dump(r, open(a.out, "w"), indent=1)
