"""Fast direct convex-position Monte Carlo of P_n(K) for an ARBITRARY convex polygon.

Purpose: stress-test the OPEN n >= 6 extremal conjecture
    P_n(triangle) <= P_n(K) <= P_n(ellipse)
at n = 6, 7 over a broad family of NON-regular bodies (TASK_N6_EXTREMAL.md, Part 1).

Design
------
* Uniform sampling in a convex polygon: fan triangulation from the CENTROID, triangle chosen
  by area-weighted binary search, then the standard reflected-barycentric map.
* Convex-position test: the centroid of the n sample points is strictly inside their hull, so
  sorting the points by angle about it yields the hull order whenever all n are hull vertices.
  Hence  {all n in convex position}  <=>  {the angularly sorted polygon turns left at every
  vertex}.  (=>: hull order = angular order about an interior point.  <=: all-left-turns plus
  winding number 1 about the sort centre => convex, and then every point is a hull vertex.)
  Angles are compared through the exact-monotone "diamond" pseudo-angle
  d(dx,dy) = dx/(|dx|+|dy|) - 1 if dy<0 else 1 - dx/(|dx|+|dy|),  which is strictly increasing
  in the true angle on (-pi,pi] and uses no transcendentals.  The left-turn tests are exact
  float orientation determinants -- the same predicate as convex_position.py's tester A.
* Each sample draws 7 points; the first 6 give P_6 and all 7 give P_7 in ONE pass (each
  estimator is separately unbiased; they are correlated, which is harmless here).
* RNG: counter-based splitmix64, re-seeded per chunk, so the result depends only on
  (seed, nchunk, nsamp) and NOT on the number of threads -- bit-reproducible.

Validation (`python mcp6.py --validate`) reproduces 91/900, 49/400, the exact regular m-gon
P_6 table of results/n6_mgon_P6_final.json, Valtr's P_7 for triangle/parallelogram, and the
disk values, and checks affine invariance and agreement with convex_position.py's independent
monotone-chain hull tester.
"""
from __future__ import annotations

import math
import numpy as np
from numba import njit, prange, get_num_threads, set_num_threads

# ----------------------------------------------------------------- exact anchors
P6_TRI = 91.0 / 900.0
P6_SQ = 49.0 / 400.0
P6_DISK = 1.0 - (146400 * math.pi**2 - 473473) / (11520 * math.pi**4)   # Marckert


def _valtr_triangle(n: int) -> float:
    """Valtr 1995: P_n(triangle) = 2^n (3n-3)! / ((n-1)!^3 (2n)!)."""
    from math import factorial as f
    return 2**n * f(3 * n - 3) / (f(n - 1) ** 3 * f(2 * n))


def _valtr_parallelogram(n: int) -> float:
    """Valtr 1996: P_n(parallelogram) = ( C(2n-2, n-1) / n! )^2."""
    from math import comb, factorial as f
    return (comb(2 * n - 2, n - 1) / f(n)) ** 2


P7_TRI = _valtr_triangle(7)
P7_SQ = _valtr_parallelogram(7)

GOLD = np.uint64(0x9E3779B97F4A7C15)
M1 = np.uint64(0xBF58476D1CE4E5B9)
M2 = np.uint64(0x94D049BB133111EB)
U53 = np.float64(1.0) / np.float64(9007199254740992.0)   # 2^-53


@njit(inline='always', cache=True)
def _mix(z):
    z = (z ^ (z >> np.uint64(30))) * M1
    z = (z ^ (z >> np.uint64(27))) * M2
    return z ^ (z >> np.uint64(31))


@njit(inline='always', cache=True)
def _conv_ok(x, y, n, ang, idx):
    """True iff the first n points (x[:n], y[:n]) are in convex position."""
    cx = 0.0
    cy = 0.0
    for i in range(n):
        cx += x[i]
        cy += y[i]
    cx /= n
    cy /= n
    for i in range(n):
        dx = x[i] - cx
        dy = y[i] - cy
        r = dx / (abs(dx) + abs(dy))
        ang[i] = (r - 1.0) if dy < 0.0 else (1.0 - r)
        idx[i] = i
    # insertion sort by pseudo-angle (n <= 7)
    for i in range(1, n):
        k = idx[i]
        a = ang[k]
        j = i - 1
        while j >= 0 and ang[idx[j]] > a:
            idx[j + 1] = idx[j]
            j -= 1
        idx[j + 1] = k
    # every turn must be a left turn
    for i in range(n):
        a = idx[i]
        b = idx[(i + 1) % n]
        c = idx[(i + 2) % n]
        cr = (x[b] - x[a]) * (y[c] - y[a]) - (y[b] - y[a]) * (x[c] - x[a])
        if cr <= 0.0:
            return False
    return True


@njit(cache=True)
def _chunk(ax, ay, bx, by, cx_, cy_, cw, lo, hi, seed, ch, out):
    """Serial worker for samples [lo,hi); writes (hits6, hits7) into out[0], out[1]."""
    nt = cw.shape[0]
    x = np.empty(7)
    y = np.empty(7)
    ang = np.empty(7)
    idx = np.empty(7, np.int64)
    st = _mix(np.uint64(seed) * GOLD + np.uint64(ch + 1) * np.uint64(0xD1B54A32D192ED03))
    h6 = 0
    h7 = 0
    for _s in range(lo, hi):
        for k in range(7):
            st = st + GOLD
            u = np.float64(_mix(st) >> np.uint64(11)) * U53
            loi = 0
            hii = nt - 1
            while loi < hii:
                mid = (loi + hii) >> 1
                if u < cw[mid]:
                    hii = mid
                else:
                    loi = mid + 1
            t = loi
            st = st + GOLD
            a1 = np.float64(_mix(st) >> np.uint64(11)) * U53
            st = st + GOLD
            a2 = np.float64(_mix(st) >> np.uint64(11)) * U53
            if a1 + a2 > 1.0:
                a1 = 1.0 - a1
                a2 = 1.0 - a2
            x[k] = ax[t] + (bx[t] - ax[t]) * a1 + (cx_[t] - ax[t]) * a2
            y[k] = ay[t] + (by[t] - ay[t]) * a1 + (cy_[t] - ay[t]) * a2
        if _conv_ok(x, y, 6, ang, idx):
            h6 += 1
        if _conv_ok(x, y, 7, ang, idx):
            h7 += 1
    out[0] = h6
    out[1] = h7


@njit(parallel=True, cache=True)
def _run(ax, ay, bx, by, cx_, cy_, cw, nsamp, seed, nchunk):
    """Returns (hits6, hits7). Triangles (a,b,c) with cumulative area weights cw."""
    res = np.zeros((nchunk, 2), np.int64)
    per = (nsamp + nchunk - 1) // nchunk
    for ch in prange(nchunk):
        lo = ch * per
        hi = min(lo + per, nsamp)
        if hi > lo:
            _chunk(ax, ay, bx, by, cx_, cy_, cw, lo, hi, seed, ch, res[ch])
    return res[:, 0].sum(), res[:, 1].sum()


# ----------------------------------------------------------------- polygon prep
def prep(V):
    """V: (nv,2) convex polygon, any orientation. Returns the centroid-fan triangle arrays."""
    V = np.asarray(V, dtype=np.float64)
    if _signed_area(V) < 0:
        V = V[::-1].copy()
    g = V.mean(axis=0)
    nv = len(V)
    A = np.repeat(g[None, :], nv, axis=0)
    B = V
    C = np.roll(V, -1, axis=0)
    ar = 0.5 * np.abs((B[:, 0] - A[:, 0]) * (C[:, 1] - A[:, 1]) - (B[:, 1] - A[:, 1]) * (C[:, 0] - A[:, 0]))
    keep = ar > 0
    A, B, C, ar = A[keep], B[keep], C[keep], ar[keep]
    cw = np.cumsum(ar) / ar.sum()
    cw[-1] = 1.0
    return (np.ascontiguousarray(A[:, 0]), np.ascontiguousarray(A[:, 1]),
            np.ascontiguousarray(B[:, 0]), np.ascontiguousarray(B[:, 1]),
            np.ascontiguousarray(C[:, 0]), np.ascontiguousarray(C[:, 1]),
            np.ascontiguousarray(cw))


def _signed_area(V):
    x, y = V[:, 0], V[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def normalize(V):
    """Affine-normalise to centroid 0 and area 1 (P_n is affine invariant; this only keeps
    the coordinates well-scaled)."""
    V = np.asarray(V, float)
    V = V - V.mean(axis=0)
    a = abs(_signed_area(V))
    return V / math.sqrt(a)


def mc(V, nsamp=int(2e8), seed=1, nchunk=4096, nthreads=16):
    """Direct MC of P_6 and P_7 for the convex polygon V. Returns a dict."""
    nthr0 = get_num_threads()
    set_num_threads(min(nthreads, nthr0))
    ax, ay, bx, by, cx_, cy_, cw = prep(V)
    nsamp = int(nsamp)
    h6, h7 = _run(ax, ay, bx, by, cx_, cy_, cw, nsamp, int(seed), int(nchunk))
    set_num_threads(nthr0)
    p6 = h6 / nsamp
    p7 = h7 / nsamp
    return dict(P6=p6, P7=p7, se6=math.sqrt(max(p6 * (1 - p6), 1e-300) / nsamp),
                se7=math.sqrt(max(p7 * (1 - p7), 1e-300) / nsamp),
                hits6=int(h6), hits7=int(h7), samples=nsamp, seed=int(seed))


# ----------------------------------------------------------------- shape helpers
def reg(m, R=1.0, phase=0.0):
    return np.array([[R * math.cos(2 * math.pi * k / m + phase),
                      R * math.sin(2 * math.pi * k / m + phase)] for k in range(m)])


def disk_poly(m=1024):
    return reg(m)


if __name__ == "__main__":
    import argparse, json, sys, time
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--samples", type=float, default=4e8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if not a.validate:
        print("nothing to do; use --validate")
        sys.exit(0)

    N = int(a.samples)
    exact = json.load(open("../results/n6_mgon_P6_final.json"))
    exact_m = {int(r["m"]): float(r["P6_float"]) for r in exact["rows"] if "P6_float" in r}
    rows = []
    allok = True
    print(f"validation: {N:.3g} samples/body, seed {a.seed}, <= {a.threads} threads\n")
    print(f"{'body':26s} {'P_6 MC':>12s} {'+-':>9s} {'P_6 exact':>12s} {'z':>7s}   "
          f"{'P_7 MC':>11s} {'P_7 exact':>11s} {'z':>7s}")

    def row(name, V, e6, e7=None):
        global allok
        t0 = time.time()
        r = mc(V, N, a.seed, nthreads=a.threads)
        z6 = (r["P6"] - e6) / r["se6"] if e6 is not None else float("nan")
        z7 = (r["P7"] - e7) / r["se7"] if e7 is not None else float("nan")
        ok = (e6 is None or abs(z6) < 4) and (e7 is None or abs(z7) < 4)
        allok &= ok
        print(f"{name:26s} {r['P6']:12.8f} {r['se6']:9.2e} "
              f"{(e6 if e6 is not None else float('nan')):12.8f} {z6:+7.2f}   "
              f"{r['P7']:11.8f} {(e7 if e7 is not None else float('nan')):11.8f} {z7:+7.2f}"
              f"   {'OK' if ok else 'FAIL'}  [{time.time()-t0:.1f}s]")
        rows.append(dict(body=name, **r, exact_P6=e6, exact_P7=e7, z6=z6, z7=z7, ok=bool(ok)))
        return r

    row("triangle(right)", [[0, 0], [1, 0], [0, 1]], P6_TRI, P7_TRI)
    row("triangle(equilateral)", reg(3), P6_TRI, P7_TRI)
    row("triangle(sliver 1:100)", [[0, 0], [1, 0], [0.37, 0.01]], P6_TRI, P7_TRI)
    row("square(unit)", [[0, 0], [1, 0], [1, 1], [0, 1]], P6_SQ, P7_SQ)
    row("parallelogram(sheared)", [[0, 0], [1, 0], [4.3, 2.7], [3.3, 2.7]], P6_SQ, P7_SQ)
    for m in range(5, 13):
        row(f"regular {m}-gon", reg(m), exact_m.get(m), None)
    for m in (14, 16, 18, 20):
        row(f"regular {m}-gon", reg(m), exact_m.get(m), None)
    row("disk(1024-gon)", disk_poly(1024), P6_DISK, None)
    # affine invariance of a NON-regular body: same P_6 under a shear+scale
    K = np.array([[0., 0.], [2., 0.], [1.3, 1.], [0.2, 1.]])
    M = np.array([[1.7, -0.9], [0.4, 2.1]])
    r1 = row("trapezoid", K, None, None)
    r2 = row("trapezoid(affine image)", K @ M.T, r1["P6"], r1["P7"])
    print()
    # independent tester cross-check (convex_position.py monotone-chain hull, numpy sampler)
    import convex_position as CP
    for body, n, tgt in (("pentagon", 6, exact_m.get(5)), ("octagon", 6, exact_m.get(8)),
                         ("triangle", 7, P7_TRI), ("square", 7, P7_SQ)):
        rr = CP.estimate(body, n, int(min(N, 2e8)), seed=a.seed + 5, batch=2_000_000, both=False)
        z = (rr["p_hat"] - tgt) / rr["std_err"]
        ok = abs(z) < 4
        allok &= ok
        print(f"  [independent tester] {body:9s} n={n}  {rr['p_hat']:.8f} +- {rr['std_err']:.1e}"
              f"  vs {tgt:.8f}  z={z:+.2f}  {'OK' if ok else 'FAIL'}")
        rows.append(dict(body="XCHECK-" + body, n=n, P6=rr["p_hat"], se6=rr["std_err"],
                         exact_P6=tgt, z6=z, ok=bool(ok)))
    print("\nALL VALIDATION OK" if allok else "\nVALIDATION FAILURE")
    if a.out:
        json.dump(dict(samples=N, seed=a.seed, all_ok=bool(allok), rows=rows),
                  open(a.out, "w"), indent=1)
    sys.exit(0 if allok else 1)
