"""Two-chord Blaschke-Petkantschin integral for E[A_4^2 & convex] on a convex POLYGON.

Derivation (docs/N6_INGREDIENTS.md):  for four i.i.d. uniform points in K the convex-position
hull area is A_4 = (1/2)|det(d1,d2)| with d1,d2 the diagonals, and

    E[A_4^2 & convex] = (3/4) V^{-6} int_{l1} int_{l2} sin^2(phi1-phi2)
                          G(chord1;sigma1) G(chord2;sigma2) 1{X in K} dl1 dl2 ,
    dl = dp dphi,  phi in [0,pi),  X = l1 ^ l2,
    G(a,b,sigma) = ((b-a)^5 - (b-sigma)^5 - (sigma-a)^5)/10 .

THIS FILE'S CHANGE OF VARIABLES.  Replace (p1,p2) by the crossing point X:
dp1 dp2 = |sin(phi1-phi2)| dX, and 1{X in K} becomes the domain.  With
u_i = dist(X, dK along +t(phi_i)),  v_i = dist(X, dK along -t(phi_i)),
sigma_i - a_i = v_i, b_i - sigma_i = u_i, b_i - a_i = u_i + v_i, so

    I = int_0^pi int_0^pi sin^2(D) |sin(D)| F(phi1,phi2) dphi1 dphi2,     D = phi1 - phi2,
    F(phi1,phi2) = int_K G1(X) G2(X) dX,   G_i = ((u_i+v_i)^5 - u_i^5 - v_i^5)/10,
    E[A_4^2 & convex] = (3/4) V^{-6} I .

WHY THIS IS ACCURATE.  Fix (phi1,phi2).  Cut K by the lines through each vertex parallel to
t(phi1) (a family of parallel lines; equivalently strips in p1 = X.n(phi1)) and likewise for
phi2.  Inside one cell the four rays +-t1, +-t2 all leave K through fixed edges, so
u1,v1,u2,v2 are LINEAR in X and G1*G2 is a polynomial of total degree 10.  A degree-11 exact
Duffy-Gauss rule (6x6) on a fan triangulation of the cell therefore integrates it EXACTLY.
Only the 2-D angular integral is numerical; it is panelled at every critical angle (directions
of vertex differences and their perpendiculars -- where the strip combinatorics changes) and at
phi2 = phi1 (the |sin D|^3 kink), then Gauss-Legendre on each panel.  Everything in the
integrand is >= 0, so the sum has no cancellation.

Symmetry: for a body invariant under rotation by 2pi/m, phi1 may be restricted to [0, pi/q)
with multiplicity q, q = m (m odd) / m/2 (m even) -- the induced shift group on phi in R/pi Z.

Usage:
    python n6_twochord_polygon.py --body triangle --nphi1 24 --nphi2 24
    python n6_twochord_polygon.py --mgon 6 --nphi1 32 --nphi2 32 --workers 16
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

# ------------------------------------------------------------------ working precision -------
# DT is the numpy dtype used for ALL geometry and quadrature arithmetic.  float64 gives ~15
# correct digits; np.longdouble (x86 80-bit, 64-bit mantissa) gives ~18-19, which is what makes
# PSLQ on the resulting E[A_4^2 & convex] possible.  Gauss-Legendre nodes and the sines/cosines
# of the angles are generated in mpmath and cast down, so they carry the full mantissa.
DT = np.float64


def set_dtype(name):
    global DT, _RULE_CACHE, _GL_CACHE
    DT = {"f64": np.float64, "ld": np.longdouble}[name]
    _RULE_CACHE = {}
    _GL_CACHE = {}


_RULE_CACHE = {}
_GL_CACHE = {}


def gauss(n):
    """Gauss-Legendre nodes/weights on [-1,1] at full DT precision (mpmath-generated)."""
    key = (n, DT)
    if key not in _GL_CACHE:
        if DT is np.float64:
            x, w = leggauss(n)
            _GL_CACHE[key] = (x.astype(DT), w.astype(DT))
        else:
            import mpmath as mp
            xs, ws = [], []
            with mp.workdps(40):        # Newton on P_n, seeded from the float64 asymptotic root
                for i in range(1, n + 1):
                    x = mp.mpf(math.cos(math.pi * (i - 0.25) / (n + 0.5)))
                    for _ in range(100):
                        p0, p1 = mp.mpf(1), x
                        for k in range(2, n + 1):
                            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                        dp = n * (x * p1 - p0) / (x * x - 1)
                        dx = p1 / dp
                        x -= dx
                        if abs(dx) < mp.mpf(10) ** -35:
                            break
                    p0, p1 = mp.mpf(1), x
                    for k in range(2, n + 1):
                        p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                    dp = n * (x * p1 - p0) / (x * x - 1)
                    xs.append(np.longdouble(mp.nstr(x, 25)))
                    ws.append(np.longdouble(mp.nstr(2 / ((1 - x * x) * dp * dp), 25)))
            _GL_CACHE[key] = (np.array(xs[::-1], dtype=DT), np.array(ws[::-1], dtype=DT))
    return _GL_CACHE[key]


def cossin(phi):
    """cos/sin of a DT angle at full DT precision."""
    if DT is np.float64:
        return DT(math.cos(phi)), DT(math.sin(phi))
    import mpmath as mp
    with mp.workdps(35):
        # mp.mpf(python_float) is EXACT; going through repr()/a decimal string is NOT
        # (it re-rounds and costs ~1e-17 -- this was a real precision bug).
        hi_ = float(phi)
        a = mp.mpf(hi_) + mp.mpf(float(phi - DT(hi_)))
        return DT(mp.nstr(mp.cos(a), 25)), DT(mp.nstr(mp.sin(a), 25))

# ---------------------------------------------------------------- geometry helpers ----------

def edges_of(V):
    """Outward unit normals nu_j and offsets h_j for a CCW convex polygon (inside: X.nu <= h)."""
    V = np.asarray(V, DT)
    A = V
    B = np.roll(V, -1, axis=0)
    e = B - A
    nu = np.stack([e[:, 1], -e[:, 0]], axis=1)
    nu /= np.hypot(nu[:, 0], nu[:, 1])[:, None]
    h = np.einsum('ij,ij->i', A, nu)
    return nu, h


def clip_halfplane(P, n, c):
    """Clip convex polygon P (list of 2-arrays) to {X . n <= c}."""
    if len(P) < 3:
        return []
    d = [p[0] * n[0] + p[1] * n[1] - c for p in P]
    out = []
    k = len(P)
    for i in range(k):
        j = (i + 1) % k
        di, dj = d[i], d[j]
        if di <= 0.0:
            out.append(P[i])
        if (di < 0.0 < dj) or (dj < 0.0 < di):
            t = di / (di - dj)
            out.append(P[i] + t * (P[j] - P[i]))
    return out


def clip_slab(P, n, lo, hi):
    P = clip_halfplane(P, n, hi)
    if len(P) < 3:
        return []
    return clip_halfplane(P, -n, -lo)


def poly_area(P):
    if len(P) < 3:
        return DT(0)
    s = DT(0)
    for i in range(len(P)):
        j = (i + 1) % len(P)
        s += P[i][0] * P[j][1] - P[j][0] * P[i][1]
    return s / 2


def raycast_edge(X, d, nu, h):
    """Index of the edge through which the ray X + s d (s>0) leaves the convex polygon."""
    den = nu @ d
    num = h - nu @ X
    s = np.where(den > 1e-13, num / np.where(den > 1e-13, den, DT(1)), DT(np.inf))
    return int(np.argmin(s))


# ------------------------------------------------- degree-11 exact rule on a triangle --------

def duffy_rule(ng=6):  # cached per (ng, DT) by the caller
    """Nodes (uu,vv) and weights ww with  int_T f dX = |cross(B-A,C-A)| * sum ww_k f(X_k),
    X = A + uu(B-A) + vv(C-A).  Exact for polynomials of degree <= 2*ng-1 in each variable
    (=> degree 10 in X plus the (1-s) Jacobian factor)."""
    xs, ws = gauss(ng)
    s = (xs + 1) / 2
    w = ws / 2
    S, T = np.meshgrid(s, s, indexing='ij')
    WS, WT = np.meshgrid(w, w, indexing='ij')
    uu = S.ravel()
    vv = (T * (1.0 - S)).ravel()
    ww = ((1.0 - S) * WS * WT).ravel()
    return uu, vv, ww


# ---------------------------------------------------------------- the inner F integral -------

def F_of_angles(Vnp, nu, h, phi1, phi2, rule, tol):
    """F(phi1,phi2) = int_K G1(X) G2(X) dX, computed exactly cell by cell."""
    uu, vv, ww = rule
    c1, s1 = cossin(phi1)
    c2, s2 = cossin(phi2)
    n1 = np.array([c1, s1], dtype=DT)
    t1 = np.array([-s1, c1], dtype=DT)
    n2 = np.array([c2, s2], dtype=DT)
    t2 = np.array([-s2, c2], dtype=DT)

    def breaks(n):
        p = np.sort(Vnp @ n)
        out = [p[0]]
        for x in p[1:]:
            if x - out[-1] > tol:
                out.append(x)
        return out

    br1 = breaks(n1)
    br2 = breaks(n2)
    base = [Vnp[i] for i in range(len(Vnp))]

    parts = []
    for i in range(len(br1) - 1):
        strip = clip_slab(base, n1, br1[i], br1[i + 1])
        if len(strip) < 3:
            continue
        for j in range(len(br2) - 1):
            cell = clip_slab(strip, n2, br2[j], br2[j + 1])
            if len(cell) < 3:
                continue
            ar = poly_area(cell)
            if ar < tol * tol:
                continue
            C = np.array(cell)
            cen = C.mean(axis=0)
            # the four exit edges (constant over the cell)
            jf1 = raycast_edge(cen, t1, nu, h)
            jb1 = raycast_edge(cen, -t1, nu, h)
            jf2 = raycast_edge(cen, t2, nu, h)
            jb2 = raycast_edge(cen, -t2, nu, h)
            # u = (h_j - X.nu_j)/(d.nu_j) is affine in X
            def affine(jj, d):
                den = nu[jj] @ d                             # keep DT precision (a float()
                return (h[jj] / den, -nu[jj] / den)          # cast here cost ~1e-16 relative)
            c_u1, g_u1 = affine(jf1, t1)
            c_v1, g_v1 = affine(jb1, -t1)
            c_u2, g_u2 = affine(jf2, t2)
            c_v2, g_v2 = affine(jb2, -t2)

            # fan triangulation
            A = C[0]
            Bs = C[1:-1]
            Cs = C[2:]
            e1 = Bs - A
            e2 = Cs - A
            cr = np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])          # (ntri,)
            X = (A[None, None, :] + uu[None, :, None] * e1[:, None, :]
                 + vv[None, :, None] * e2[:, None, :])                      # (ntri,nq,2)
            u1 = c_u1 + X @ g_u1
            v1 = c_v1 + X @ g_v1
            u2 = c_u2 + X @ g_u2
            v2 = c_v2 + X @ g_v2
            # G = ((u+v)^5 - u^5 - v^5)/10 = u v (u+v) (u^2+uv+v^2) / 2 -- the factored form is
            # cancellation-free (in sliver cells u/v can reach 1e10 and the naive difference
            # loses ~10 digits).
            G1 = u1 * v1 * (u1 + v1) * (u1 * u1 + u1 * v1 + v1 * v1)
            G2 = u2 * v2 * (u2 + v2) * (u2 * u2 + u2 * v2 + v2 * v2)
            parts.append(cr @ ((G1 * G2) @ ww))
    if not parts:
        return DT(0)
    return np.sum(np.array(parts, dtype=DT)) / 4      # pairwise summation, positive terms


# ---------------------------------------------------------------- angular quadrature ---------

def PI():
    if DT is np.float64:
        return DT(math.pi)
    import mpmath as mp
    with mp.workdps(40):
        return DT(mp.nstr(mp.pi, 25))


def critical_angles(Vnp):
    """Directions mod pi where the strip combinatorics of the width/ray structure changes:
    the directions of all vertex differences and their perpendiculars.  Computed in mpmath and
    returned at DT precision -- these are the panel ends AND (for 0 and pi) the exact limits of
    the phi integrations, so a float64 pi would cost a relative 1e-16."""
    import mpmath as mp
    pi = PI()
    m = len(Vnp)
    with mp.workdps(40):
        raw = set()
        for i in range(m):
            for j in range(i + 1, m):
                dx = mp.mpf(str(Vnp[j][0])) - mp.mpf(str(Vnp[i][0]))
                dy = mp.mpf(str(Vnp[j][1])) - mp.mpf(str(Vnp[i][1]))
                a = mp.atan2(dy, dx) % mp.pi
                raw.add(a)
                raw.add((a + mp.pi / 2) % mp.pi)
        vals = sorted(DT(mp.nstr(a, 25)) for a in raw)
    out = [DT(0)]
    for x in vals:
        if x > out[-1] + DT(1e-12) and x < pi - DT(1e-12):
            out.append(x)
    out.append(pi)
    return out


def panels(nodes, lo, hi, grade=0, ratio=0.25, nsub=1):
    """Panels between consecutive critical angles.  `grade` levels of geometric refinement
    toward BOTH endpoints of every panel: the ray directions t(phi) become parallel to an edge
    exactly at the critical angles, which makes F(phi1,phi2) only finitely smooth there
    (thin slivers with very long chords).  Graded panels restore fast convergence."""
    pts = [lo] + [x for x in nodes if lo + 1e-12 < x < hi - 1e-12] + [hi]
    out = []
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-13:
            continue
        if grade <= 0:
            out.append((a, b))
            continue
        w = (b - a) / 2
        left = []
        for _ in range(grade):
            left.append(w)
            w = w * DT(ratio)
        # left half: a + w_k accumulating from the smallest outwards
        xs = [a] + [a + w for w in reversed(left)]
        ys = [b] + [b - w for w in reversed(left)]
        pt = sorted(set(xs + ys))
        out.extend(zip(pt[:-1], pt[1:]))
    out = [(a, b) for a, b in out if b - a > 1e-14]
    if nsub > 1:
        fine = []
        for a, b in out:
            step = (b - a) / nsub
            for k in range(nsub):
                fine.append((a + k * step, a + (k + 1) * step if k < nsub - 1 else b))
        out = fine
    return out


_CTX = {}


def _init(V, ng, nphi2, crit, tol, grade=0, dtname="f64", nsub=1):
    set_dtype(dtname)
    Vnp = np.asarray(V, DT)
    nu, h = edges_of(Vnp)
    _CTX.update(Vnp=Vnp, nu=nu, h=h, rule=duffy_rule(ng), nphi2=nphi2, crit=crit, tol=tol,
                grade=grade, nsub=nsub)


def _phi1_slice(arg):
    """Inner phi2 integral (times the phi1 Gauss weight) for one phi1 node."""
    phi1, w1 = arg
    Vnp, nu, h, rule = _CTX['Vnp'], _CTX['nu'], _CTX['h'], _CTX['rule']
    nphi2, crit, tol = _CTX['nphi2'], _CTX['crit'], _CTX['tol']
    xs, ws = gauss(nphi2)
    terms = []
    nodes = sorted(set(crit) | {phi1})
    for a, b in panels(nodes, DT(0), crit[-1], _CTX['grade'], nsub=_CTX['nsub']):
        if b - a < 1e-13:
            continue
        mid = (a + b) / 2
        half = (b - a) / 2
        for x, w in zip(xs, ws):
            phi2 = mid + half * x
            sD = abs(cossin(phi1 - phi2)[1])
            if sD < 1e-15:
                continue
            F = F_of_angles(Vnp, nu, h, phi1, phi2, rule, tol)
            terms.append(half * w * sD ** 3 * F)
    return w1 * np.sum(np.array(terms, dtype=DT))


def two_chord_polygon(V, nphi1=24, nphi2=24, ng=6, qsym=1, workers=1, verbose=False, grade=0,
                      dtname="f64", nsub=1):
    """I = int int sin^2 |sin| F dphi1 dphi2 over [0,pi)^2, using the qsym-fold rotation
    symmetry to restrict phi1 to [0, pi/qsym).  Returns (I, E[A_4^2 & convex])."""
    set_dtype(dtname)
    Vnp = np.asarray(V, DT)
    area = poly_area([Vnp[i] for i in range(len(Vnp))])
    scale = math.sqrt(abs(float(area)))
    tol = 1e-11 * scale
    crit = critical_angles(Vnp)
    hi = PI() / qsym
    # critical angles inside the fundamental domain
    crit_in = [x for x in crit if DT(0) < x < hi]
    xs, ws = gauss(nphi1)
    jobs = []
    for a, b in panels(crit_in, DT(0), hi, grade, nsub=nsub):
        mid, half = (a + b) / 2, (b - a) / 2
        for x, w in zip(xs, ws):
            jobs.append((mid + half * x, half * w))
    _init(Vnp, ng, nphi2, crit, tol, grade, dtname, nsub)
    t0 = time.time()
    if workers > 1:
        import multiprocessing as mp
        with mp.Pool(workers, initializer=_init,
                     initargs=(Vnp, ng, nphi2, crit, tol, grade, dtname, nsub)) as pool:
            parts = pool.map(_phi1_slice, jobs, chunksize=1)
    else:
        parts = [_phi1_slice(j) for j in jobs]
    I = qsym * np.sum(np.array(parts, dtype=DT))
    EA = 3 * I / (4 * area ** 6)
    if verbose:
        print(f"    [{len(jobs)} phi1 nodes, {time.time()-t0:.1f}s]", file=sys.stderr)
    return I, EA


# ---------------------------------------------------------------- bodies ---------------------

def regular_mgon(m):
    import mpmath as mp
    with mp.workdps(40):
        return [(np.longdouble(mp.nstr(mp.cos(2 * mp.pi * k / m), 25)),
                 np.longdouble(mp.nstr(mp.sin(2 * mp.pi * k / m), 25))) for k in range(m)]


def qsym_of(m):
    return m if m % 2 else m // 2


BODIES = {
    "triangle": ([(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], 1),          # right triangle, no symmetry used
    "square":   ([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], 1),
    "eqtri":    (regular_mgon(3), 3),
    "sq4":      (regular_mgon(4), 2),
}

from fractions import Fraction as _F
EXACT = {"triangle": _F(119, 4500), "square": _F(1307, 54000),
         "eqtri": _F(119, 4500), "sq4": _F(1307, 54000)}


EXACT_MGON = {   # E[A_4^2 & convex], identified from the 80-bit values (see docs/N6_LANDSCAPE.md)
    3: "119/4500", 4: "1307/54000", 5: "1769/112500 + 577*sqrt(5)/168750",
    6: "403891/17496000", 8: "35743/2764800 + 48793*sqrt(2)/6912000"}


def validate(out_path=None, workers=8):
    """Reproduce the two anchors 119/4500 and 1307/54000 through this code path, in float64 and
    in 80-bit, with and without the rotation-symmetry reduction, and at two cubature orders."""
    import mpmath as mp
    mp.mp.dps = 30
    tgt = {"eqtri": mp.mpf(119) / 4500, "triangle": mp.mpf(119) / 4500,
           "sq4": mp.mpf(1307) / 54000, "square": mp.mpf(1307) / 54000}
    rows, allok = [], True
    for name in ("eqtri", "triangle", "sq4", "square"):
        V, q = BODIES[name]
        for dt, nphi, ng in (("f64", 32, 6), ("ld", 32, 6), ("ld", 48, 6), ("ld", 32, 10)):
            _, EA = two_chord_polygon(V, nphi, nphi, ng, q, workers=workers, grade=1, dtname=dt)
            v = mp.mpf(str(EA))
            rel = abs(v - tgt[name]) / tgt[name]
            ok = rel < mp.mpf(10) ** -12
            allok &= bool(ok)
            rows.append(dict(body=name, qsym=q, dtype=dt, nphi=nphi, ng=ng, grade=1,
                             value=mp.nstr(v, 21), exact=str(tgt[name]),
                             rel_diff=float(rel), ok=bool(ok)))
            print(f"  {name:9s} q={q} {dt} nphi={nphi} ng={ng}: {mp.nstr(v,21)}  "
                  f"rel={mp.nstr(rel,3)}  {'OK' if ok else 'FAIL'}")
    print("ALL ANCHORS OK" if allok else "ANCHOR FAILURE")
    if out_path:
        import json as _json
        hi = {}
        try:
            hi = _json.load(open(out_path.replace("n6_twochord_polygon", "n6_tc_hiprec")))
        except Exception:
            pass
        _json.dump(dict(
            what="two-chord Blaschke-Petkantschin integral for E[A_4^2 & convex] on a convex polygon",
            formula="E[A_4^2&conv] = (3/4) V^-6 int int sin^2(D)|sin(D)| F dphi1 dphi2, "
                    "F = int_K G1 G2 dX, G = u v (u+v)(u^2+uv+v^2)/2",
            anchor_validation=rows, all_anchors_ok=allok,
            mgon_values={k: dict(grids=v, exact=EXACT_MGON.get(int(k)))
                         for k, v in hi.items()},
            note="mgon_values: A/B/C are three independent (nphi,grade) settings in 80-bit; "
                 "their spread is the error bar (~1e-20 absolute).  Assembled P_6 table: "
                 "results/n6_mgon_P6_final.json"),
            open(out_path, "w"), indent=1)
        print("wrote", out_path)
    return rows, allok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--body", default=None)
    ap.add_argument("--mgon", type=int, default=None)
    ap.add_argument("--nphi1", type=int, default=24)
    ap.add_argument("--nphi2", type=int, default=24)
    ap.add_argument("--ng", type=int, default=6)
    ap.add_argument("--nosym", action="store_true")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--grade", type=int, default=0)
    ap.add_argument("--ld", action="store_true", help="80-bit extended precision")
    ap.add_argument("--nsub", type=int, default=1, help="uniform sub-panels per critical panel")
    ap.add_argument("--json", default=None)
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()

    if args.validate:
        _, ok = validate(args.json, workers=max(1, args.workers))
        sys.exit(0 if ok else 1)

    if args.mgon:
        V, q = regular_mgon(args.mgon), qsym_of(args.mgon)
        label = f"regular {args.mgon}-gon"
        exact = EXACT.get({3: "eqtri", 4: "sq4"}.get(args.mgon, ""), None)
    else:
        V, q = BODIES[args.body]
        label = args.body
        exact = EXACT.get(args.body)
    if args.nosym:
        q = 1
    t0 = time.time()
    I, EA = two_chord_polygon(V, args.nphi1, args.nphi2, args.ng, q, args.workers,
                              grade=args.grade, dtname="ld" if args.ld else "f64",
                              nsub=args.nsub)
    dt = time.time() - t0
    print(f"{label:16s} q={q} nphi1={args.nphi1} nphi2={args.nphi2} ng={args.ng} "
          f"grade={args.grade} nsub={args.nsub} {'ld' if args.ld else 'f64'}: "
          f"E[A_4^2 & convex] = {repr(EA)}   ({dt:.1f}s)")
    if exact:
        import mpmath as mp
        mp.mp.dps = 30
        ex = mp.mpf(exact.numerator) / exact.denominator if hasattr(exact, 'numerator') else mp.mpf(exact)
        v = mp.mpf(str(EA))
        print(f"{'':16s} exact = {mp.nstr(ex,20)}   rel.diff = {mp.nstr(abs(v-ex)/ex,3)}")
    if args.json:
        json.dump(dict(body=label, qsym=q, nphi1=args.nphi1, nphi2=args.nphi2, ng=args.ng,
                       grade=args.grade, dtype=("ld" if args.ld else "f64"),
                       I=repr(I), EA42_convex=repr(EA), exact=str(exact), seconds=dt),
                  open(args.json, "w"), indent=1)


if __name__ == "__main__":
    main()
