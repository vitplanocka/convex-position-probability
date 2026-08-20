"""Route M: convex-position probability from HULL-AREA MOMENTS (structurally
independent of the direct tester in convex_position.py).

Identity (exchangeability + conditional independence; Efron/Buchta type):
  Let A_k = area(conv(x_1..x_k)) / |K|, N_n = #vertices of conv(x_1..x_n),
  R(m, j) := E[ 1{x_1..x_m in convex position} * A_m^j ].
  Then  P(N_n = m) = C(n, m) * R(m, n-m)          (m = 3..n)
  and   E[A_k^j] = sum_{m=3}^{k} C(k, m) R(m, j + k - m)
  so    R(k, j) = E[A_k^j] - sum_{m=3}^{k-1} C(k, m) R(m, j + k - m),
        P_n = R(n, 0) = 1 - sum_{m=3}^{n-1} C(n, m) R(m, n-m).
Special cases:  P_4 = 1 - 4 E[A_3]  (Sylvester),  P_5 = 1 - 5 E[A_4] + 10 E[A_3^2].

So P_n is a linear combination of the plain moments E[A_k^j], k=3..n-1,
j=1..n-k, all of which we estimate by Monte Carlo with an EXACT hull area
(Jarvis march, no shared code with the direct testers).  Error bars by batch
means.  Same anchors must reproduce (Valtr for square/triangle, Sylvester disk).

Usage: python route_moments.py --body disk --n 5 --samples 2e7
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from fractions import Fraction
from math import comb

import numpy as np

from convex_position import BODIES  # samplers only

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

BODY_AREA = {
    "square": 1.0, "disk": math.pi, "triangle": 0.5,
    "pentagon": 5 / 2 * math.sin(2 * math.pi / 5), "hexagon": 6 / 2 * math.sin(2 * math.pi / 6),
    "octagon": 8 / 2 * math.sin(2 * math.pi / 8),
}


if HAVE_NUMBA:
    @njit(cache=True)
    def _jarvis_area(x, y, n):
        # leftmost (then lowest) point
        start = 0
        for i in range(1, n):
            if x[i] < x[start] or (x[i] == x[start] and y[i] < y[start]):
                start = i
        area2 = 0.0
        p = start
        first = True
        prev = -1
        count = 0
        while True:
            q = -1
            for r in range(n):
                if r == p:
                    continue
                if q == -1:
                    q = r
                    continue
                cr = (x[q] - x[p]) * (y[r] - y[p]) - (y[q] - y[p]) * (x[r] - x[p])
                if cr < 0:
                    q = r
                elif cr == 0:
                    # collinear: take the farther one
                    dq = (x[q] - x[p]) ** 2 + (y[q] - y[p]) ** 2
                    dr = (x[r] - x[p]) ** 2 + (y[r] - y[p]) ** 2
                    if dr > dq:
                        q = r
            # shoelace accumulation with vertex p -> q
            area2 += x[p] * y[q] - x[q] * y[p]
            p = q
            count += 1
            if p == start or count > n + 1:
                break
        return abs(area2) * 0.5

    @njit(parallel=True, cache=True)
    def _hull_areas(P, out):
        m, n, _ = P.shape
        for s in prange(m):
            out[s] = _jarvis_area(P[s, :, 0], P[s, :, 1], n)

    def hull_areas(P):
        out = np.empty(P.shape[0])
        _hull_areas(np.ascontiguousarray(P), out)
        return out
else:
    def hull_areas(P):  # slow fallback via scipy
        from scipy.spatial import ConvexHull
        return np.array([ConvexHull(P[s]).volume for s in range(P.shape[0])])


def moment_table(body, n, samples, seed, batches=20, batch_size=200_000):
    """Estimate E[A_k^j] for k=3..n-1, j=1..n-k. Returns dict[(k,j)] -> (mean, batch-means array)."""
    sampler, dim = BODIES[body]
    assert dim == 2
    area = BODY_AREA[body]
    rng = np.random.default_rng(seed)
    per_batch = samples // batches
    res = {}
    for k in range(3, n):
        js = list(range(1, n - k + 1))
        bm = np.zeros((batches, len(js)))
        for b in range(batches):
            done = 0
            acc = np.zeros(len(js))
            while done < per_batch:
                m = min(batch_size, per_batch - done)
                P = sampler(rng, m, k)
                A = hull_areas(P) / area
                for ji, j in enumerate(js):
                    acc[ji] += np.sum(A ** j)
                done += m
            bm[b] = acc / per_batch
        for ji, j in enumerate(js):
            res[(k, j)] = (bm[:, ji].mean(), bm[:, ji].copy())
    return res


def solve_P(n, moments):
    """Given moments[(k,j)] = E[A_k^j] (floats or arrays), return P_n via the recursion."""
    R = {}
    for k in range(3, n):
        for j in range(1, n - k + 1):
            val = moments[(k, j)]
            for m in range(3, k):
                val = val - comb(k, m) * R[(m, j + k - m)]
            R[(k, j)] = val
    P = 1.0
    for m in range(3, n):
        P = P - comb(n, m) * R[(m, n - m)]
    return P


def coefficients(n):
    """P_n = 1 + sum c_{k,j} E[A_k^j]; return the exact rational coefficients (for the record)."""
    # symbolic linear solve using Fractions: represent each E[A_k^j] as basis vector
    keys = [(k, j) for k in range(3, n) for j in range(1, n - k + 1)]
    idx = {kj: i for i, kj in enumerate(keys)}
    def vec():
        return [Fraction(0)] * len(keys)
    R = {}
    for k in range(3, n):
        for j in range(1, n - k + 1):
            v = vec(); v[idx[(k, j)]] = Fraction(1)
            for m in range(3, k):
                w = R[(m, j + k - m)]
                v = [a - comb(k, m) * b for a, b in zip(v, w)]
            R[(k, j)] = v
    P = vec()
    for m in range(3, n):
        w = R[(m, n - m)]
        P = [a - comb(n, m) * b for a, b in zip(P, w)]
    return {keys[i]: P[i] for i in range(len(keys)) if P[i] != 0}


def estimate(body, n, samples, seed, batches=20):
    t0 = time.time()
    mom = moment_table(body, n, samples, seed, batches)
    means = {kj: v[0] for kj, v in mom.items()}
    bms = {kj: v[1] for kj, v in mom.items()}
    P = solve_P(n, means)
    # error: apply the (linear) recursion to each batch's means, take std of batch estimates
    Pb = solve_P(n, bms)  # arrays broadcast through the linear recursion
    se = float(np.std(Pb, ddof=1) / math.sqrt(batches))
    return {
        "route": "moments", "body": body, "n": n, "samples_per_k": samples, "seed": seed,
        "p_hat": float(P), "std_err": se, "seconds": time.time() - t0,
        "moments": {f"E[A_{k}^{j}]": float(v) for (k, j), v in means.items()},
        "coefficients": {f"E[A_{k}^{j}]": str(c) for (k, j), c in coefficients(n).items()},
    }


def validate(samples=2_000_000, seed=4242):
    from anchors import ANCHORS
    plan = [("square", 4, "square_n4"), ("square", 5, "square_n5"), ("square", 6, "square_n6"),
            ("triangle", 4, "triangle_n4"), ("triangle", 5, "triangle_n5"), ("triangle", 6, "triangle_n6"),
            ("disk", 4, "disk_n4")]
    ok_all = True
    rows = []
    for body, n, key in plan:
        _, exact = ANCHORS[key]
        r = estimate(body, n, samples, seed)
        z = (r["p_hat"] - exact) / r["std_err"]
        ok = abs(z) < 4
        ok_all &= ok
        rows.append({**r, "exact": exact, "z": z, "anchor": key, "ok": ok})
        print(f"{key:12s} exact={exact:.9f} routeM={r['p_hat']:.9f} +- {r['std_err']:.2e} z={z:+.2f} {'OK' if ok else 'FAIL'} ({r['seconds']:.1f}s)")
    print("ROUTE M: ALL ANCHORS OK" if ok_all else "ROUTE M: ANCHOR FAILURE")
    return rows, ok_all


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--body", default="disk")
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--samples", type=float, default=2e6, help="samples per k")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--batches", type=int, default=20)
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    print(f"numba: {HAVE_NUMBA}", file=sys.stderr)
    if a.validate:
        rows, ok = validate(int(a.samples), a.seed)
        if a.out:
            json.dump({"rows": rows, "all_ok": ok}, open(a.out, "w"), indent=1)
        sys.exit(0 if ok else 1)
    print("coefficients:", {f"E[A_{k}^{j}]": str(c) for (k, j), c in coefficients(a.n).items()})
    r = estimate(a.body, a.n, int(a.samples), a.seed, a.batches)
    print(json.dumps(r, indent=1))
    if a.out:
        json.dump(r, open(a.out, "w"), indent=1)
