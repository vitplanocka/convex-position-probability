r"""B2/B3: the n = 5 extremal landscape.

By identity (II),  P_5(K) = 1 - 10 F(K)  with   F(K) := E[A_3] - E[A_3^2],
so the n = 5 case of the extremal conjecture  P_triangle(n) <= P_K(n) <= P_ellipse(n)
is EXACTLY the two-sided inequality

    61/(96 pi^2) = 0.06438116...  <=  F(K)  <=  5/72 = 0.069444...    (conjecture at n=5)

with equality for ellipses (left) and triangles (right).  Both E[A_3] (Blaschke 1917) and
E[A_3^2] = (3/2) det Sigma/|K|^2 (the planar isotropic constant) are maximised by triangles
and minimised by ellipses, so F is a difference of two same-direction functionals and the
inequality is NOT immediate.

Route P (src/polygon_exact.py) gives F EXACTLY for any convex polygon, so this is a search
over exact values, not a Monte-Carlo scan.
"""
import json, math, sys, time
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import minimize

from polygon_exact import (moments_fast, regular_polygon, cov_det_area,
                           cov_det_area_exact, whiten)

F_TRI = 5 / 72
F_DISK = 61 / (96 * math.pi ** 2)
P5_TRI = 11 / 36
P5_DISK = 1 - 305 / (48 * math.pi ** 2)


def F_of_polygon(V, normalise=True):
    """F(K) = E[A_3] - E[A_3^2] for the convex polygon V (CCW).

    Everything here is affine invariant, so we first put the polygon in isotropic position
    (centroid 0, covariance I).  This matters: on a near-degenerate polygon the float64
    covariance suffers catastrophic cancellation -- an almost-collinear triangle of area
    2.4e-5 with edges ~3 gave E[A_3^2] = 0.0071 instead of 1/72, which made a Nelder-Mead
    search report a spurious violation of the triangle bound.  The covariance is additionally
    computed in exact Fraction arithmetic (`cov_det_area_exact`), which has no cancellation.
    """
    V = np.asarray(V, float)
    if normalise:
        V, ok = whiten(V)
        if not ok:
            return np.nan
    r = moments_fast(V, nmax=4)          # only needs T_2
    if abs(r["T"][0] - 2) > 1e-9 or abs(r["E[N_3]"] - 3) > 1e-9:
        return np.nan                    # route P's built-in exactness checks failed
    d, A = cov_det_area_exact(V)
    return r["E[A_3]"] - 1.5 * d / A ** 2


def hull_ccw(pts):
    pts = np.asarray(pts, float)
    h = ConvexHull(pts)
    return pts[h.vertices]               # scipy returns CCW in 2-D


def F_of_points(pts):
    try:
        V = hull_ccw(pts)
    except Exception:
        return np.nan
    if len(V) < 3:
        return np.nan
    return F_of_polygon(V)


# --------------------------------------------------------------------- families

def fam_regular(mmax=40):
    return {f"regular {m}-gon": regular_polygon(m) for m in range(3, mmax + 1)}


def fam_truncated_triangle(ts):
    """Triangle with all three corners cut at parameter t (t=0 triangle, t=1/2 -> hexagon
    degenerating; t=1/3 gives the affinely-regular hexagon)."""
    T = np.array([[0., 0.], [1., 0.], [0., 1.]])
    out = {}
    for t in ts:
        V = []
        for i in range(3):
            a, b, c = T[i], T[(i + 1) % 3], T[(i + 2) % 3]
            V.append(a + t * (b - a))
            V.append(a + t * (c - a))
        out[f"truncated triangle t={t:.4f}"] = hull_ccw(np.array(V))
    return out


def fam_cap(thetas, npts=400):
    """Circular segment (cap) cut off by a chord subtending angle theta.
    theta -> 2 pi is the disk; theta -> 0 is a parabolic segment (affinely rigid)."""
    out = {}
    for th in thetas:
        a = np.linspace(-th / 2, th / 2, npts)
        V = np.stack([np.cos(a), np.sin(a)], -1)
        out[f"cap theta={th:.4f}"] = hull_ccw(V)
    return out


def fam_stadium(aspects, npts=200):
    out = {}
    for a in aspects:
        right = [(a + math.cos(t), math.sin(t)) for t in np.linspace(-math.pi / 2, math.pi / 2, npts)]
        left = [(-a + math.cos(t), math.sin(t)) for t in np.linspace(math.pi / 2, 3 * math.pi / 2, npts)]
        out[f"stadium a={a:.3f}"] = hull_ccw(np.array(right + left))
    return out


def fam_tri_disk_minkowski(ts, npts=300):
    """Minkowski combination (1-t)*Triangle + t*Disk (support-function interpolation)."""
    out = {}
    T = np.array([[math.cos(2 * math.pi * k / 3 + math.pi / 2), math.sin(2 * math.pi * k / 3 + math.pi / 2)]
                  for k in range(3)])
    ang = np.linspace(0, 2 * math.pi, npts, endpoint=False)
    u = np.stack([np.cos(ang), np.sin(ang)], -1)
    hT = (u @ T.T).max(1)
    for t in ts:
        h = (1 - t) * hT + t * 1.0
        # reconstruct the body as the intersection of half-planes {x.u <= h}: use the dual
        # (a polygon through consecutive half-plane intersections)
        V = []
        for i in range(npts):
            j = (i + 1) % npts
            M = np.array([u[i], u[j]])
            try:
                V.append(np.linalg.solve(M, [h[i], h[j]]))
            except np.linalg.LinAlgError:
                pass
        out[f"Minkowski (1-t)T + tD, t={t:.3f}"] = hull_ccw(np.array(V))
    return out


def fam_pushed_vertex(mm, ss):
    """Regular m-gon with one vertex pushed radially outward by factor s."""
    out = {}
    for m in mm:
        for s in ss:
            V = np.array(regular_polygon(m))
            V[0] = V[0] * s
            out[f"{m}-gon, vertex x{s:.2f}"] = hull_ccw(V)
    return out


def random_polygons(n, rng, kmin=3, kmax=12):
    out = {}
    i = 0
    while len(out) < n:
        i += 1
        k = rng.integers(kmin, kmax + 1)
        pts = rng.normal(size=(k + 6, 2))
        try:
            V = hull_ccw(pts)
        except Exception:
            continue
        if len(V) < 3:
            continue
        out[f"random#{i} ({len(V)} vertices)"] = V
    return out


# ------------------------------------------------------------------- optimisation

def optimise(m, sign, seed, ntries=6):
    """Nelder-Mead over 2m free coordinates; F is computed on the convex hull, so the
    search is over all convex polygons with at most m vertices.  sign=+1 maximises F."""
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(ntries):
        x0 = rng.normal(size=2 * m)
        def f(x, _s=sign, _m=m):
            v = F_of_points(x.reshape(_m, 2))
            return 1e3 if not np.isfinite(v) else -_s * v
        r = minimize(f, x0, method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-10, fatol=1e-14))
        val = sign * -r.fun
        if best is None or sign * (val - best[0]) > 0:
            best = (val, r.x.reshape(m, 2))
    return best


def selftest(rng, n=200):
    """F must be EXACTLY 5/72 on every affine image of a triangle and 49/144-based on every
    parallelogram.  Random affine images with condition numbers up to 1e6."""
    worst_t = worst_p = 0.0
    T0 = np.array([[0., 0.], [1., 0.], [0., 1.]])
    S0 = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    Ftri, Fsq = 5 / 72, (1 - 49 / 144) / 10
    for _ in range(n):
        M = rng.normal(size=(2, 2)) @ np.diag(np.exp(rng.uniform(-7, 7, 2)))
        b = rng.normal(size=2) * 10
        if abs(np.linalg.det(M)) < 1e-12:
            continue
        for V0, ref, tag in [(T0, Ftri, "t"), (S0, Fsq, "p")]:
            V = V0 @ M.T + b
            if np.linalg.det(np.array([V[1] - V[0], V[2] - V[0]])) < 0:
                V = V[::-1]
            e = abs(F_of_polygon(V) - ref)
            if tag == "t":
                worst_t = max(worst_t, e)
            else:
                worst_p = max(worst_p, e)
    return worst_t, worst_p


if __name__ == "__main__":
    t0 = time.time()
    _rng0 = np.random.default_rng(1)
    wt, wp = selftest(_rng0)
    print(f"self-test: max |F - 5/72| over 200 random affine triangles      = {wt:.2e}")
    print(f"           max |F - F(square)| over 200 random parallelograms   = {wp:.2e}")
    assert wt < 1e-12 and wp < 1e-12, "affine-invariance self-test FAILED"
    rng = np.random.default_rng(20260819)
    fams = {}
    fams.update(fam_regular(40))
    fams.update(fam_truncated_triangle(np.linspace(0.001, 0.499, 25)))
    fams.update(fam_cap(np.linspace(0.05, 2 * math.pi - 1e-9, 25)))
    fams.update(fam_stadium(np.concatenate([[1e-3], np.linspace(0.05, 8, 20)])))
    fams.update(fam_tri_disk_minkowski(np.linspace(0, 1, 21)))
    fams.update(fam_pushed_vertex([3, 4, 5, 6, 8], [1.05, 1.2, 1.5, 2.0, 3.0, 6.0]))
    fams.update(random_polygons(400, rng))
    fams["half-disk"] = hull_ccw(np.array([(math.cos(t), math.sin(t)) for t in np.linspace(0, math.pi, 400)]))

    rows = []
    for name, V in fams.items():
        F = F_of_polygon(np.asarray(V, float))
        rows.append(dict(name=name, nv=len(V), F=F, P5=1 - 10 * F,
                         viol_low=F_DISK - F, viol_high=F - F_TRI))
    rows.sort(key=lambda r: r["F"])
    print(f"=== n=5 landscape: {len(rows)} bodies, F = E[A_3]-E[A_3^2], P_5 = 1-10F ===")
    print(f"conjectured window: F in [{F_DISK:.12f}, {F_TRI:.12f}]  <=>  "
          f"P_5 in [{P5_TRI:.12f}, {P5_DISK:.12f}]")
    print("\n-- 12 smallest F (closest to the ELLIPSE bound / largest P_5) --")
    for r in rows[:12]:
        print(f"  F={r['F']:.12f}  P_5={r['P5']:.12f}  gap_below={r['viol_low']:+.3e}  {r['name']}")
    print("\n-- 12 largest F (closest to the TRIANGLE bound / smallest P_5) --")
    for r in rows[-12:]:
        print(f"  F={r['F']:.12f}  P_5={r['P5']:.12f}  gap_above={r['viol_high']:+.3e}  {r['name']}")
    viol = [r for r in rows if r["viol_low"] > 1e-12 or r["viol_high"] > 1e-12]
    print(f"\nviolations of the conjectured window (tol 1e-12): {len(viol)}")
    for r in viol[:20]:
        print("   ", r["name"], r["F"], r["viol_low"], r["viol_high"])

    print("\n=== direct optimisation over convex polygons (Nelder-Mead on the hull) ===")
    opt = {}
    for m in [3, 4, 5, 6, 8, 10]:
        for sign, lab in [(+1, "max"), (-1, "min")]:
            val, V = optimise(m, sign, seed=1000 + m + (0 if sign > 0 else 500))
            opt[f"{lab} over <= {m}-gons"] = dict(F=val, P5=1 - 10 * val, nv=int(len(hull_ccw(V))))
            flag = ""
            if sign > 0 and val > F_TRI + 1e-11:
                flag = "  *** EXCEEDS TRIANGLE BOUND ***"
            if sign < 0 and val < F_DISK - 1e-11:
                flag = "  *** BELOW ELLIPSE BOUND ***"
            print(f"  {lab} over <= {m}-gons: F = {val:.12f}  (triangle {F_TRI:.12f}, "
                  f"disk {F_DISK:.12f})  P_5 = {1-10*val:.12f}{flag}")

    json.dump(dict(window=dict(F_disk=F_DISK, F_tri=F_TRI, P5_tri=P5_TRI, P5_disk=P5_DISK),
                   bodies=rows, optimisation=opt, seconds=time.time() - t0),
              open("../results/n5_landscape.json", "w"), indent=1, default=float)
    print(f"\nwrote ../results/n5_landscape.json  ({time.time()-t0:.0f}s)")
