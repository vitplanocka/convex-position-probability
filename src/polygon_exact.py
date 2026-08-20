r"""EXACT Sylvester / convex-position quantities for an arbitrary convex POLYGON.

Route P (new in this session).  Everything reduces to the integrals

    J_k := \int_{lines} L^3 [ c^k + (1-c)^k ] dG          (dG = dp dtheta)

over the space of lines meeting K, where L is the chord length and c the
fraction of |K| cut off on one side.  Reason: with two i.i.d. uniform points,
Blaschke-Petkantschin gives  dx1 dx2 = |t1 - t2| dt1 dt2 dG, and integrating the
positions along the chord gives  \int\int |t1-t2| dt1 dt2 = L^3/3.  Hence for any f

    E[ f(c) + f(1-c) ] = (1/(3|K|^2)) \int L^3 [f(c)+f(1-c)] dG,

and the Renyi-Sulanke edge count of the hull of n i.i.d. points is
    E[N_n] = C(n,2) E[ c^{n-2} + (1-c)^{n-2} ] = C(n,2) T_{n-2},
    T_k := J_k / (3|K|^2).
Efron's identity then gives  E[A_{n-1}] = 1 - E[N_n]/n.

THE POINT: for a convex polygon the line space decomposes into cells indexed by
the unordered pair {i,j} of edges the line crosses, and inside a cell the line is
parametrised by the two boundary points
    P(u) = A_i + u d_i,  Q(v) = A_j + v d_j,   (u,v) in [0,1]^2.
Santalo's formula  dp dtheta = sin(phi_i) sin(phi_j) ds_i ds_j / L  becomes

    dG = |d_i x (Q-P)| |d_j x (Q-P)| / L^3  du dv,

so the L^3 CANCELS and

    J_k = sum_{i<j} \int_0^1 \int_0^1 W(u,v) [c^k + (1-c)^k] du dv,
    W(u,v) = -(d_i x (Q - A_i)) (d_j x (A_j - P))  >= 0,

where W is a product of a linear function of v and a linear function of u, and
c(u,v) = S(u,v)/|K| with S the shoelace area of the piece [P, V_{i+1},..,V_j, Q]
-- a polynomial of degree <= 2.  The integrand is therefore a POLYNOMIAL of
degree <= 2 + 2k in (u,v): the integral is elementary and exact.  No quadrature
error at all: Gauss-Legendre with enough nodes is exact for polynomials.

Sanity checks built in (`self_test`):
  T_0 = 2                      (equivalently \int L^3 dG = 3|K|^2)
  E[N_3] = 3 T_1 = 3           (three points always span a triangle)
  \int dG = perimeter          (Cauchy; the one check needing real quadrature)
  E[A_3] = 1 - (3/2) T_2 :  square 11/144, triangle 1/12, and -> 35/(48 pi^2) as m -> oo.
"""
from __future__ import annotations

import math
from math import comb

import numpy as np


# ---------------------------------------------------------------- numeric core

def _gauss_legendre(nodes, dps=None):
    """Gauss-Legendre nodes/weights on [0,1]; mpmath if dps given, else numpy."""
    if dps is None:
        x, w = np.polynomial.legendre.leggauss(nodes)
        return (x + 1) / 2, w / 2
    from mpmath import mp, matrix
    mp.dps = dps
    # mpmath has no leggauss; build from numpy nodes then Newton-polish in mp
    from mpmath import legendre, mpf, diff
    xs, ws = [], []
    x0, _ = np.polynomial.legendre.leggauss(nodes)
    for xi in x0:
        x = mpf(float(xi))
        for _ in range(60):
            f = legendre(nodes, x)
            fp = diff(lambda t: legendre(nodes, t), x)
            dx = f / fp
            x -= dx
            if abs(dx) < mpf(10) ** (-dps - 5):
                break
        wgt = 2 / ((1 - x ** 2) * diff(lambda t: legendre(nodes, t), x) ** 2)
        xs.append((x + 1) / 2)
        ws.append(wgt / 2)
    return xs, ws


def shoelace(V):
    V = list(V)
    s = 0
    for i in range(len(V)):
        a, b = V[i], V[(i + 1) % len(V)]
        s = s + a[0] * b[1] - b[0] * a[1]
    return s / 2


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def J_integrals(V, kmax, nodes=None, dps=None, exact=False):
    r"""J_k = \int L^3 [c^k + (1-c)^k] dG for k = 0..kmax, for the CCW convex polygon V.

    exact=True  -> sympy symbolic (V entries may be sympy expressions)
    dps=<int>   -> mpmath at that precision
    otherwise   -> float64
    Gauss-Legendre with `nodes` points per axis is EXACT for the polynomial integrand
    as long as 2*nodes-1 >= 2*kmax+2.
    """
    m = len(V)
    if nodes is None:
        nodes = max(6, kmax + 3)
    if exact:
        import sympy as sp
        u, v = sp.symbols('u v')
        one = sp.Integer(1)
    A = shoelace(V)                      # > 0 for CCW
    J = [0] * (kmax + 1)
    for i in range(m):
        Ai, Bi = V[i], V[(i + 1) % m]
        di = (Bi[0] - Ai[0], Bi[1] - Ai[1])
        for j in range(i + 1, m):
            Aj, Bj = V[j], V[(j + 1) % m]
            dj = (Bj[0] - Aj[0], Bj[1] - Aj[1])
            # chain of fixed vertices between the two moving endpoints
            mid = [V[t % m] for t in range(i + 1, j + 1)]

            def cell(uu, vv, _Ai=Ai, _di=di, _Aj=Aj, _dj=dj, _mid=mid):
                P = (_Ai[0] + uu * _di[0], _Ai[1] + uu * _di[1])
                Q = (_Aj[0] + vv * _dj[0], _Aj[1] + vv * _dj[1])
                Di = cross(_di, (Q[0] - _Ai[0], Q[1] - _Ai[1]))     # >= 0
                Dj = cross(_dj, (_Aj[0] - P[0], _Aj[1] - P[1]))     # <= 0
                W = -Di * Dj
                S = shoelace([P] + _mid + [Q])
                return W, S / A

            if exact:
                W, c = cell(u, v)
                for k in range(kmax + 1):
                    f = sp.Poly(sp.expand(W * (c ** k + (one - c) ** k)), u, v)
                    acc = sp.Integer(0)
                    for (a, b), co in f.terms():
                        acc += co / sp.Integer((a + 1) * (b + 1))
                    J[k] += acc
            else:
                xs, ws = _gauss_legendre(nodes, dps)
                for a in range(nodes):
                    for b in range(nodes):
                        W, c = cell(xs[a], xs[b])
                        wt = ws[a] * ws[b] * W
                        for k in range(kmax + 1):
                            J[k] += wt * (c ** k + (1 - c) ** k)
    return J, A


def T_and_moments(V, nmax=6, nodes=None, dps=None, exact=False):
    """Returns dict with T_k, E[N_n], E[A_{n-1}], P_4, P_5 (exact in the polygon)."""
    kmax = nmax - 2
    J, A = J_integrals(V, kmax, nodes, dps, exact)
    if exact:
        import sympy as sp
        T = [sp.simplify(Jk / (3 * A ** 2)) for Jk in J]
    else:
        T = [Jk / (3 * A ** 2) for Jk in J]
    out = {"area": A, "T": T}
    for n in range(3, nmax + 1):
        EN = comb(n, 2) * T[n - 2]
        out[f"E[N_{n}]"] = EN
        out[f"E[A_{n-1}]"] = 1 - EN / n
    out["P_4"] = 1 - 4 * out["E[A_3]"]
    return out


def regular_polygon(m, exact=False):
    if exact:
        import sympy as sp
        return [(sp.cos(2 * sp.pi * k / m), sp.sin(2 * sp.pi * k / m)) for k in range(m)]
    return [(math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m)) for k in range(m)]


def regular_polygon_mp(m, dps=50):
    from mpmath import mp, cos, sin, mpf, pi
    mp.dps = dps
    return [(cos(2 * pi * k / m), sin(2 * pi * k / m)) for k in range(m)]


# -------------------------------------------------------------- exact E[A_3^2]

def cov_det_area(V):
    """(det Sigma, area) of the uniform law on the CCW polygon V (works for sympy too)."""
    m = len(V)
    A = 0
    M1 = [0, 0]
    M2 = [[0, 0], [0, 0]]
    for i in range(m):
        a, b = V[i], V[(i + 1) % m]
        o = (0 * a[0], 0 * a[1])
        T = (a[0] * b[1] - a[1] * b[0]) / 2
        g = [(a[0] + b[0] + o[0]) / 3, (a[1] + b[1] + o[1]) / 3]
        s = (a[0] + b[0] + o[0], a[1] + b[1] + o[1])
        for r in range(2):
            for cc in range(2):
                M2[r][cc] += T * (a[r] * a[cc] + b[r] * b[cc] + o[r] * o[cc] + s[r] * s[cc]) / 12
        A += T
        M1[0] += T * g[0]
        M1[1] += T * g[1]
    mu = [M1[0] / A, M1[1] / A]
    S = [[M2[r][cc] / A - mu[r] * mu[cc] for cc in range(2)] for r in range(2)]
    return S[0][0] * S[1][1] - S[0][1] * S[1][0], A


def EA3sq(V):
    """E[A_3^2] = (3/2) det Sigma / |K|^2 (exact; E[det^2] = 6 det Sigma)."""
    d, A = cov_det_area(V)
    return 3 * d / (2 * A ** 2)


def P5_from_identity(EA3, EA3sq_):
    """(II): P_5 = 1 - 10 E[A_3] + 10 E[A_3^2]."""
    return 1 - 10 * EA3 + 10 * EA3sq_


# ------------------------------------------------------------------ self-tests

def perimeter_check(V, nodes=64):
    r"""\int dG = perimeter. Needs 1/L^3, so genuine quadrature (float only)."""
    m = len(V)
    xs, ws = _gauss_legendre(nodes)
    tot = 0.0
    for i in range(m):
        Ai, Bi = V[i], V[(i + 1) % m]
        di = (Bi[0] - Ai[0], Bi[1] - Ai[1])
        for j in range(i + 1, m):
            Aj, Bj = V[j], V[(j + 1) % m]
            dj = (Bj[0] - Aj[0], Bj[1] - Aj[1])
            for a in range(nodes):
                for b in range(nodes):
                    uu, vv = xs[a], xs[b]
                    P = (Ai[0] + uu * di[0], Ai[1] + uu * di[1])
                    Q = (Aj[0] + vv * dj[0], Aj[1] + vv * dj[1])
                    L = math.hypot(Q[0] - P[0], Q[1] - P[1])
                    Di = cross(di, (Q[0] - Ai[0], Q[1] - Ai[1]))
                    Dj = cross(dj, (Aj[0] - P[0], Aj[1] - P[1]))
                    tot += ws[a] * ws[b] * (-Di * Dj) / L ** 3
    per = sum(math.hypot(V[(i + 1) % m][0] - V[i][0], V[(i + 1) % m][1] - V[i][1]) for i in range(m))
    return tot, per


def self_test():
    from fractions import Fraction
    print("=== route P self-test ===")
    tri = [(0, 0), (1, 0), (0, 1)]
    sq = [(0, 0), (1, 0), (1, 1), (0, 1)]
    for name, V, exactEA3 in [("triangle", tri, Fraction(1, 12)), ("square", sq, Fraction(11, 144))]:
        r = T_and_moments(V, nmax=6)
        print(f"{name:9s} area={r['area']:.6f} T_0={r['T'][0]:.12f} (must be 2)  "
              f"E[N_3]={r['E[N_3]']:.12f} (must be 3)")
        print(f"          E[A_3]={r['E[A_3]']:.15f}  exact {float(exactEA3):.15f}  "
              f"diff {r['E[A_3]']-float(exactEA3):+.2e}")
        print(f"          E[A_4]={r['E[A_4]']:.15f}  ratio E[A_4]/E[A_3]={r['E[A_4]']/r['E[A_3]']:.15f}")
        print(f"          P_4={r['P_4']:.15f}")
        t, p = perimeter_check(V, 40)
        print(f"          int dG = {t:.10f}  perimeter = {p:.10f}  rel {abs(t/p-1):.2e}")
    print("--- regular m-gons, float ---")
    for m in [3, 4, 5, 6, 8, 12, 24, 96]:
        r = T_and_moments(regular_polygon(m), nmax=5)
        print(f"  m={m:3d}  E[A_3]={r['E[A_3]']:.15f}  P_4={r['P_4']:.15f}  "
              f"E[A_4]/E[A_3]={r['E[A_4]']/r['E[A_3]']:.12f}")
    print(f"  disk limit: E[A_3]=35/(48pi^2)={35/(48*math.pi**2):.15f}  "
          f"P_4={1-35/(12*math.pi**2):.15f}")


if __name__ == "__main__":
    self_test()


# ------------------------------------------------- EXACT rational-vertex route

def J_integrals_exact_rational(V, kmax):
    """Exact J_k for a CCW convex polygon with RATIONAL vertices, using
    Fraction arithmetic and term-by-term integration of the polynomial integrand.
    Polynomials in (u, v) are dicts {(a, b): Fraction}.  Since E[A_3] etc. are
    AFFINE INVARIANT this covers the triangle, every parallelogram, and the
    affinely-regular hexagon -- i.e. exact rational verification for m = 3, 4, 6
    and for any polygon with rational vertices."""
    from fractions import Fraction as F

    def pmul(p, q):
        r = {}
        for (a1, b1), c1 in p.items():
            for (a2, b2), c2 in q.items():
                k = (a1 + a2, b1 + b2)
                r[k] = r.get(k, F(0)) + c1 * c2
        return {k: c for k, c in r.items() if c != 0}

    def padd(p, q):
        r = dict(p)
        for k, c in q.items():
            r[k] = r.get(k, F(0)) + c
        return {k: c for k, c in r.items() if c != 0}

    def pscale(p, s):
        return {k: c * s for k, c in p.items() if c * s != 0}

    def pconst(c):
        return {(0, 0): F(c)} if c != 0 else {}

    def pint(p):                      # int_0^1 int_0^1 . du dv
        return sum((c / F((a + 1) * (b + 1)) for (a, b), c in p.items()), F(0))

    V = [(F(x), F(y)) for x, y in V]
    m = len(V)
    A = sum(V[i][0] * V[(i + 1) % m][1] - V[(i + 1) % m][0] * V[i][1] for i in range(m)) / F(2)
    J = [F(0)] * (kmax + 1)
    U = {(1, 0): F(1)}
    Wv = {(0, 1): F(1)}
    for i in range(m):
        Ai, Bi = V[i], V[(i + 1) % m]
        di = (Bi[0] - Ai[0], Bi[1] - Ai[1])
        Px = padd(pconst(Ai[0]), pscale(U, di[0]))
        Py = padd(pconst(Ai[1]), pscale(U, di[1]))
        for j in range(i + 1, m):
            Aj, Bj = V[j], V[(j + 1) % m]
            dj = (Bj[0] - Aj[0], Bj[1] - Aj[1])
            Qx = padd(pconst(Aj[0]), pscale(Wv, dj[0]))
            Qy = padd(pconst(Aj[1]), pscale(Wv, dj[1]))
            Di = padd(pscale(padd(Qy, pconst(-Ai[1])), di[0]), pscale(padd(Qx, pconst(-Ai[0])), -di[1]))
            Dj = padd(pscale(padd(pconst(Aj[1]), pscale(Py, -1)), dj[0]),
                      pscale(padd(pconst(Aj[0]), pscale(Px, -1)), -dj[1]))
            W = pscale(pmul(Di, Dj), F(-1))
            xs = [Px] + [pconst(V[t % m][0]) for t in range(i + 1, j + 1)] + [Qx]
            ys = [Py] + [pconst(V[t % m][1]) for t in range(i + 1, j + 1)] + [Qy]
            S = {}
            for t in range(len(xs)):
                t2 = (t + 1) % len(xs)
                S = padd(S, padd(pmul(xs[t], ys[t2]), pscale(pmul(xs[t2], ys[t]), F(-1))))
            S = pscale(S, F(1, 2))
            c = pscale(S, F(1) / A)
            omc = padd(pconst(1), pscale(c, F(-1)))
            ck = pconst(1)
            ok = pconst(1)
            for k in range(kmax + 1):
                if k > 0:
                    ck = pmul(ck, c)
                    ok = pmul(ok, omc)
                J[k] += pint(pmul(W, padd(ck, ok)))
    return J, A


def exact_rational_moments(V, nmax=6):
    """Exact Fractions: T_k, E[N_n], E[A_{n-1}], P_4, E[A_3^2], P_5 for a rational polygon."""
    from fractions import Fraction as F
    J, A = J_integrals_exact_rational(V, nmax - 2)
    T = [Jk / (3 * A ** 2) for Jk in J]
    out = {"area": A, "T": T}
    for n in range(3, nmax + 1):
        EN = F(comb(n, 2)) * T[n - 2]
        out[f"E[N_{n}]"] = EN
        out[f"E[A_{n-1}]"] = 1 - EN / n
    out["P_4"] = 1 - 4 * out["E[A_3]"]
    d, AA = cov_det_area([(F(x), F(y)) for x, y in V])
    out["E[A_3^2]"] = 3 * d / (2 * AA ** 2)
    out["P_5"] = 1 - 10 * out["E[A_3]"] + 10 * out["E[A_3^2]"]
    return out


# --------------------------------------------- vectorised float route (large m)

def J_integrals_fast(V, kmax, nodes=None):
    """Vectorised float64 version of J_integrals: loops over the (u,v) Gauss grid and
    handles ALL O(m^2) edge pairs at once with numpy.  Needed for m in the hundreds
    (polygonal approximation of smooth bodies).  Same exactness argument: the integrand
    is a polynomial, so Gauss-Legendre with 2*nodes-1 >= 2*kmax+2 is exact."""
    if nodes is None:
        nodes = max(6, kmax + 3)
    V = np.asarray(V, float)
    m = len(V)
    Vn = np.roll(V, -1, axis=0)
    d = Vn - V                                   # edge vectors
    A = 0.5 * np.sum(V[:, 0] * Vn[:, 1] - Vn[:, 0] * V[:, 1])
    II, JJ = np.triu_indices(m, k=1)
    # C[i,j] = sum_{t=i+1}^{j-1} cross(V_t, V_{t+1});  prefix sums of cross(V_t,V_{t+1})
    e = V[:, 0] * Vn[:, 1] - Vn[:, 0] * V[:, 1]  # cross(V_t, V_{t+1})
    pre = np.concatenate([[0.0], np.cumsum(e)])  # pre[k] = sum_{t<k} e_t
    Cij = pre[JJ] - pre[II + 1]                  # t = i+1 .. j-1
    Ai = V[II]; di = d[II]; Aj = V[JJ]; dj = d[JJ]
    Vi1 = Vn[II]                                 # V_{i+1}
    Vj = V[JJ]                                   # V_j  (= A_j)
    xs, ws = _gauss_legendre(nodes)
    J = np.zeros(kmax + 1)
    for a in range(nodes):
        u = xs[a]
        P = Ai + u * di
        Dj = dj[:, 0] * (Aj[:, 1] - P[:, 1]) - dj[:, 1] * (Aj[:, 0] - P[:, 0])
        cPV = P[:, 0] * Vi1[:, 1] - Vi1[:, 0] * P[:, 1]
        for b in range(nodes):
            v = xs[b]
            Q = Aj + v * dj
            Di = di[:, 0] * (Q[:, 1] - Ai[:, 1]) - di[:, 1] * (Q[:, 0] - Ai[:, 0])
            W = -Di * Dj
            S = 0.5 * (cPV + Cij
                       + (Vj[:, 0] * Q[:, 1] - Q[:, 0] * Vj[:, 1])
                       + (Q[:, 0] * P[:, 1] - P[:, 0] * Q[:, 1]))
            c = S / A
            wt = ws[a] * ws[b] * W
            omc = 1.0 - c
            ck = np.ones_like(c); ok = np.ones_like(c)
            for k in range(kmax + 1):
                if k:
                    ck *= c; ok *= omc
                J[k] += np.dot(wt, ck + ok)
    return list(J), A


def moments_fast(V, nmax=6, nodes=None):
    J, A = J_integrals_fast(V, nmax - 2, nodes)
    T = [Jk / (3 * A ** 2) for Jk in J]
    out = {"area": A, "T": T}
    for n in range(3, nmax + 1):
        EN = comb(n, 2) * T[n - 2]
        out[f"E[N_{n}]"] = EN
        out[f"E[A_{n-1}]"] = 1 - EN / n
    out["P_4"] = 1 - 4 * out["E[A_3]"]
    d, AA = cov_det_area(np.asarray(V, float))
    out["E[A_3^2]"] = 1.5 * d / AA ** 2
    out["P_5"] = 1 - 10 * out["E[A_3]"] + 10 * out["E[A_3^2]"]
    return out


# ------------------------------------------ numerically safe affine normalisation

def cov_det_area_exact(V):
    """detSigma, |K| computed in exact Fraction arithmetic from the (float) vertices.
    `cov_det_area` in float64 loses catastrophic precision on near-degenerate polygons
    (the E[xx^T] - mu mu^T cancellation): a triangle of area 2.4e-5 with edges of length ~3
    returned E[A_3^2] = 0.00708 instead of 1/72.  This version has no cancellation error."""
    from fractions import Fraction as F
    d, A = cov_det_area([(F(float(x)), F(float(y))) for x, y in V])
    return float(d), float(A)


def whiten(V):
    """Affine image of the polygon with centroid 0 and covariance I.  Every quantity in this
    module is affine invariant, so this changes nothing mathematically but makes route P and
    the covariance well-conditioned.  Returns (V_whitened, ok) with ok=False if the polygon is
    degenerate beyond repair."""
    from fractions import Fraction as F
    V = np.asarray(V, float)
    m = len(V)
    Vf = [(F(float(x)), F(float(y))) for x, y in V]
    Vn = Vf[1:] + Vf[:1]
    T = [(a[0] * b[1] - a[1] * b[0]) / 2 for a, b in zip(Vf, Vn)]
    A = sum(T)
    if A <= 0:
        return V, False
    mu = [sum(t * (a[k] + b[k]) / 3 for t, a, b in zip(T, Vf, Vn)) / A for k in range(2)]
    S = [[F(0), F(0)], [F(0), F(0)]]
    for t, a, b in zip(T, Vf, Vn):
        o = (F(0), F(0))
        s = (a[0] + b[0], a[1] + b[1])
        for r in range(2):
            for c in range(2):
                S[r][c] += t * (a[r] * a[c] + b[r] * b[c] + o[r] * o[c] + s[r] * s[c]) / 12
    S = np.array([[float(S[r][c] / A - mu[r] * mu[c]) for c in range(2)] for r in range(2)])
    w, Q = np.linalg.eigh(S)
    if w.min() <= 0 or w.min() / w.max() < 1e-300:
        return V, False
    W = Q @ np.diag(w ** -0.5) @ Q.T
    return (V - np.array([float(mu[0]), float(mu[1])])) @ W.T, True
