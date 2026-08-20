"""A broad, affinely-normalised family of convex bodies for the n >= 6 extremal scan.

Every body is returned as a convex polygon (ccw, deduped).  Smooth bodies are represented by
inscribed/support polygons with NSM vertices; the regular-m-gon table shows the resulting bias
is O(a_4/m^4) with a_4 = 3.27, i.e. < 1e-9 at m = 256 -- far below the MC error of ~1e-5.

Families
  regular m-gons                       controls (must increase monotonically to the disk)
  vertex-pushed / -pulled m-gons       irregular polygons
  trapezoids, kites, general quads     the 2-parameter shape space of quadrilaterals
  irregular pentagons / hexagons
  Minkowski interpolations             triangle<->disk, square<->disk, triangle<->square, ...
  corner-truncated triangles           the approach to the conjectured MINIMISER
  support-function perturbations of the disk,  h = 1 + eps cos(k theta)
                                       the approach to the conjectured MAXIMISER
  circular segments, sectors, half-disk, stadium, lens, Reuleaux polygons
  random convex polygons               hulls of k uniform points, k = 4..12
"""
from __future__ import annotations

import math
import numpy as np

NSM = 256          # vertices used for smooth bodies
TAU = 2 * math.pi


# ------------------------------------------------------------------ geometry helpers
def hull(P):
    """Andrew monotone chain; returns the ccw hull vertices."""
    P = sorted({(round(float(x), 14), round(float(y), 14)) for x, y in np.asarray(P, float)})
    if len(P) <= 2:
        return np.array(P, float)

    def half(pts):
        out = []
        for p in pts:
            while len(out) >= 2:
                (ox, oy), (ax, ay) = out[-2], out[-1]
                if (ax - ox) * (p[1] - oy) - (ay - oy) * (p[0] - ox) <= 1e-15:
                    out.pop()
                else:
                    break
            out.append(p)
        return out

    lo = half(P)
    up = half(P[::-1])
    return np.array(lo[:-1] + up[:-1], float)


def ccw(V):
    V = np.asarray(V, float)
    x, y = V[:, 0], V[:, 1]
    if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
        V = V[::-1].copy()
    return V


def minkowski(A, B):
    """Minkowski sum of two ccw convex polygons, by merging edge vectors by direction."""
    A, B = ccw(A), ccw(B)
    eA = np.roll(A, -1, axis=0) - A
    eB = np.roll(B, -1, axis=0) - B
    E = np.vstack([eA, eB])
    ang = np.arctan2(E[:, 1], E[:, 0])
    E = E[np.argsort(ang, kind="stable")]
    P = np.cumsum(np.vstack([[0.0, 0.0], E[:-1]]), axis=0)
    return hull(P + A[np.argmin(A[:, 1] * 1e6 + A[:, 0])] + B[np.argmin(B[:, 1] * 1e6 + B[:, 0])])


def scale(V, s):
    return np.asarray(V, float) * s


def interp(A, B, t):
    """Minkowski interpolation (1-t)A + tB, area-normalised."""
    if t <= 0:
        return norm_area(A)
    if t >= 1:
        return norm_area(B)
    return norm_area(minkowski(scale(norm_area(A), 1 - t), scale(norm_area(B), t)))


def area(V):
    V = np.asarray(V, float)
    x, y = V[:, 0], V[:, 1]
    return abs(0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)))


def norm_area(V):
    """Centre at the centroid and scale to area 1 (P_n is affine invariant; cosmetic)."""
    V = ccw(V)
    V = V - V.mean(axis=0)
    return V / math.sqrt(area(V))


# ------------------------------------------------------------------ shape constructors
def reg(m, phase=0.0):
    a = TAU * np.arange(m) / m + phase
    return np.stack([np.cos(a), np.sin(a)], axis=-1)


def disk(n=NSM):
    return reg(n)


def ellipse(ab, n=NSM):
    V = reg(n).copy()
    V[:, 0] *= ab
    return V


def pushed(m, idx=0, f=1.3):
    V = reg(m).copy()
    V[idx] *= f
    return hull(V)


def trapezoid(t):
    """Legs of length 1 and t on parallel lines (t=1 -> parallelogram/square)."""
    return np.array([[-0.5, 0.0], [0.5, 0.0], [t / 2, 1.0], [-t / 2, 1.0]])


def trapezoid_skew(t, s):
    return np.array([[-0.5, 0.0], [0.5, 0.0], [s + t / 2, 1.0], [s - t / 2, 1.0]])


def kite(a, b):
    return np.array([[0.0, -a], [b, 0.0], [0.0, 1.0], [-b, 0.0]])


def truncated_triangle(t, k=3):
    """Equilateral triangle with k corners cut back a fraction t of each adjacent edge.
    t -> 0 gives the triangle, t = 1/2 the medial triangle's hexagon."""
    T = reg(3)
    out = []
    for i in range(3):
        A, B = T[i], T[(i + 1) % 3]
        if i < k:
            out.append(A + t * (B - A))
        else:
            out.append(A)
        if (i + 1) % 3 < k or (i + 1) % 3 == 0 and k == 3:
            pass
        out.append(B - t * (B - A) if ((i + 1) % 3) < k else B)
    return hull(np.array(out))


def rounded(V, r, n=NSM):
    """V Minkowski-summed with a disk of radius r (area-normalised V)."""
    return minkowski(norm_area(V), scale(disk(n), r))


def support_body(h, n=NSM):
    """Convex body from a support function h(theta) given as a callable; boundary point is
    (h cos t - h' sin t, h sin t + h' cos t).  h' by central differences."""
    t = TAU * np.arange(n) / n
    dt = 1e-6
    hv = h(t)
    hp = (h(t + dt) - h(t - dt)) / (2 * dt)
    return np.stack([hv * np.cos(t) - hp * np.sin(t), hv * np.sin(t) + hp * np.cos(t)], axis=-1)


def cos_perturbed_disk(k, eps, n=NSM):
    """h(theta) = 1 + eps cos(k theta).  Convex iff h + h'' = 1 - eps k^2 cos(k theta) >= 0,
    i.e. eps <= 1/k^2.  k = 1 is a translate of the disk, k = 2 an affine image (ellipse) to
    first order, so k >= 3 is where genuinely new shape appears."""
    return support_body(lambda t: 1.0 + eps * np.cos(k * t), n)


def circular_segment(half_angle, n=NSM):
    """The part of the unit disk cut off by a chord; half_angle = pi gives the full disk,
    pi/2 the half-disk."""
    a = np.linspace(-half_angle, half_angle, n)
    return np.stack([np.cos(a), np.sin(a)], axis=-1)


def circular_sector(half_angle, n=NSM):
    a = np.linspace(-half_angle, half_angle, n)
    P = np.stack([np.cos(a), np.sin(a)], axis=-1)
    return np.vstack([[[0.0, 0.0]], P]) if half_angle <= math.pi / 2 else hull(np.vstack([[[0.0, 0.0]], P]))


def stadium(L, n=NSM):
    """Minkowski sum of a horizontal segment of length L with the unit disk."""
    return minkowski(np.array([[-L / 2, 0.0], [L / 2, 0.0], [L / 2, 1e-12], [-L / 2, 1e-12]]), disk(n))


def lens(d, n=NSM):
    """Intersection of two unit disks whose centres are distance d apart (0 < d < 2)."""
    a = math.acos(d / 2)
    t = np.linspace(-a, a, n // 2)
    right = np.stack([np.cos(t) - d / 2, np.sin(t)], axis=-1)
    left = np.stack([-np.cos(t) + d / 2, -np.sin(t)], axis=-1)
    return hull(np.vstack([right, left]))


def reuleaux(m, n=NSM):
    """Reuleaux polygon of odd order m: arcs of radius = the long diagonal, centred at the
    opposite vertices of a regular m-gon."""
    assert m % 2 == 1
    V = reg(m, phase=math.pi / 2)
    R = np.linalg.norm(V[0] - V[(0 + (m - 1) // 2) % m])
    pts = []
    per = max(8, n // m)
    for i in range(m):
        c = V[i]
        p = V[(i + (m + 1) // 2) % m]
        q = V[(i + (m - 1) // 2) % m]
        a0 = math.atan2(p[1] - c[1], p[0] - c[0])
        a1 = math.atan2(q[1] - c[1], q[0] - c[0])
        while a1 < a0:
            a1 += TAU
        if a1 - a0 > math.pi:
            a0, a1 = a1 - TAU, a0
        for a in np.linspace(a0, a1, per):
            pts.append([c[0] + R * math.cos(a), c[1] + R * math.sin(a)])
    return hull(np.array(pts))


def rand_convex(k, rng, mode="disk"):
    if mode == "disk":
        r = np.sqrt(rng.random(k))
        t = rng.random(k) * TAU
        P = np.stack([r * np.cos(t), r * np.sin(t)], axis=-1)
    elif mode == "square":
        P = rng.random((k, 2))
    else:
        P = rng.standard_normal((k, 2))
    return hull(P)


# ------------------------------------------------------------------ the catalogue
def catalogue(seed=20260819):
    """Returns a list of (name, family, vertices) with vertices area-normalised."""
    rng = np.random.default_rng(seed)
    B = []

    def add(name, fam, V):
        V = np.asarray(V, float)
        if len(V) < 3:
            return
        V = norm_area(hull(V))
        if len(V) < 3 or not np.all(np.isfinite(V)):
            return
        B.append((name, fam, V))

    # -- controls: regular m-gons
    for m in list(range(3, 21)) + [24, 32, 48, 64, 128, 256, 1024]:
        add(f"regular-{m}gon", "regular", reg(m))
    add("ellipse-3:1", "regular", ellipse(3.0))          # affine to the disk: control

    # -- vertex-pushed / -pulled regular m-gons
    for m in (3, 4, 5, 6, 7, 8, 10, 12):
        for f in (0.55, 0.7, 0.85, 1.15, 1.3, 1.6, 2.0):
            add(f"push{m}-{f}", "pushed", pushed(m, 0, f))
    # two adjacent / opposite vertices pushed
    for m in (5, 6, 8, 10, 12):
        for f in (0.7, 1.3, 1.8):
            V = reg(m).copy(); V[0] *= f; V[1] *= f
            add(f"push2adj{m}-{f}", "pushed", V)
            V = reg(m).copy(); V[0] *= f; V[m // 2] *= f
            add(f"push2opp{m}-{f}", "pushed", V)

    # -- quadrilaterals
    for t in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.3, 1.7, 2.5):
        add(f"trapezoid-{t}", "quad", trapezoid(t))
    for t in (0.2, 0.4, 0.6, 0.8):
        for s in (0.3, 0.7, 1.2):
            add(f"trapskew-{t}-{s}", "quad", trapezoid_skew(t, s))
    for a in (0.3, 0.6, 1.0, 1.6, 2.5):
        for b in (0.4, 0.8, 1.5):
            add(f"kite-{a}-{b}", "quad", kite(a, b))
    for i in range(20):                                   # random quadrilaterals
        add(f"randquad-{i}", "quad", rand_convex(4, rng))

    # -- pentagons / hexagons, irregular
    for i in range(20):
        add(f"randpent-{i}", "pent", rand_convex(5, rng))
    for i in range(20):
        add(f"randhex-{i}", "hex", rand_convex(6, rng))

    # -- Minkowski interpolations (one-parameter families)
    ts = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
    for t in ts:
        add(f"tri-disk-{t}", "interp", interp(reg(3), disk(), t))
        add(f"sq-disk-{t}", "interp", interp(reg(4), disk(), t))
        add(f"tri-sq-{t}", "interp", interp(reg(3), reg(4), t))
        add(f"pent-disk-{t}", "interp", interp(reg(5), disk(), t))
        add(f"tri-hex-{t}", "interp", interp(reg(3), reg(6), t))

    # -- the approach to the conjectured MINIMISER (triangle)
    for t in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.45):
        add(f"trunc3-{t}", "near-tri", truncated_triangle(t, 3))
        add(f"trunc1-{t}", "near-tri", truncated_triangle(t, 1))
        add(f"trunc2-{t}", "near-tri", truncated_triangle(t, 2))
    for r in (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2):
        add(f"roundtri-{r}", "near-tri", rounded(reg(3), r))

    # -- the approach to the conjectured MAXIMISER (disk)
    for k in (3, 4, 5, 6, 7, 8, 10, 12):
        for c in (0.2, 0.45, 0.7, 0.9, 0.99):
            add(f"cosdisk-k{k}-c{c}", "near-disk", cos_perturbed_disk(k, c / k**2))
    for r in (0.02, 0.05, 0.1, 0.2, 0.4):                # nearly-round polygons
        for m in (5, 6, 8, 12):
            add(f"round{m}gon-{r}", "near-disk", rounded(reg(m), r))
    for ab in (1.02, 1.1):                                # affine controls: must equal the disk
        add(f"ellipse-{ab}", "near-disk", ellipse(ab))

    # -- classical non-regular smooth bodies
    for ha in (0.3, 0.5, 0.7, 0.9, 1.1, 1.3, math.pi / 2, 1.8, 2.1, 2.4, 2.7, 3.0):
        add(f"segment-{ha:.2f}", "smooth", circular_segment(ha))
    for ha in (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, math.pi / 2):   # convex only for ha <= pi/2
        add(f"sector-{ha:.2f}", "smooth", circular_sector(ha))
    for L in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        add(f"stadium-{L}", "smooth", stadium(L))
    for d in (0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 1.9):
        add(f"lens-{d}", "smooth", lens(d))
    for m in (3, 5, 7, 9):
        add(f"reuleaux-{m}", "smooth", reuleaux(m))
    # cone-like / drop-like: disk cap + apex
    for hgt in (1.2, 1.5, 2.0, 3.0, 5.0):
        c = disk(NSM)
        add(f"icecream-{hgt}", "smooth", hull(np.vstack([c, [[0.0, hgt]]])))

    # -- random convex polygons, 4..12 vertices, three point processes
    for k in range(4, 13):
        for mode in ("disk", "square", "gauss"):
            for i in range(6):
                V = rand_convex(k, rng, mode)
                add(f"rand-{mode}-{k}-{i}", "random", V)
    return B


if __name__ == "__main__":
    B = catalogue()
    from collections import Counter
    print(f"{len(B)} bodies")
    for fam, c in sorted(Counter(f for _, f, _ in B).items()):
        print(f"  {fam:11s} {c}")
    nv = Counter(len(V) for _, _, V in B)
    print("vertex counts:", dict(sorted(nv.items())))
