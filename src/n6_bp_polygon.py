"""Exact triangle-area moments E[A_3^k] for a convex polygon via the Blaschke-Petkantschin
width-function integral (derivation + disk validation in n6_bp_moments.py):

    E[A_3^k] = (1/V^{k+3}) * (2^{1-k}/((k+2)(k+3)))
               * int_0^pi int_p w(p,phi)^{k+3} [ int_t |t-p|^k w(t,phi) dt ] dp dphi

w(.,phi) = width (slice length) function of the polygon in normal direction phi; piecewise
linear, so the inner t- and p-integrals are polynomial on each piece and done EXACTLY by a
fixed 8-node Gauss-Legendre rule (exact to degree 15 >= 2k+4).  Only the phi-integral is
adaptive, subdivided at edge-normal directions where the vertex ordering changes.
High precision (mpmath) + PSLQ -> exact rationals for square/triangle; validated on k=1,2.
"""
import math
import mpmath as mp

mp.mp.dps = 40

# ---- fixed high-precision Gauss-Legendre rule (n nodes, exact to degree 2n-1) --------------
def _gauss_legendre(n):
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = mp.mpf(math.cos(math.pi * (i - 0.25) / (n + 0.5)))  # float seed
        for _ in range(80):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        nodes.append(x)
        weights.append(2 / ((1 - x * x) * dp * dp))
    return nodes, weights

_GLN, _GLW = _gauss_legendre(8)


def _gl(f, a, b):
    if b <= a:
        return mp.mpf(0)
    h = (b - a) / 2
    m = (a + b) / 2
    return h * mp.fsum(w * f(m + h * x) for x, w in zip(_GLN, _GLW))


# ---- polygon width function as linear pieces in a given normal direction -------------------
def _rot(V, phi):
    c, s = mp.cos(phi), mp.sin(phi)
    return [(c * x + s * y, -s * x + c * y) for (x, y) in V]


def _pieces(Vr):
    """Return list of (l, r, A, B) with width w(t) = A + B t on [l, r], covering the polygon's
    projection onto the normal (x) axis."""
    n = len(Vr)
    imin = min(range(n), key=lambda i: (Vr[i][0], Vr[i][1]))
    imax = max(range(n), key=lambda i: (Vr[i][0], Vr[i][1]))
    chainA, i = [], imin
    while True:
        chainA.append(Vr[i])
        if i == imax:
            break
        i = (i + 1) % n
    chainB, i = [], imin
    while True:
        chainB.append(Vr[i])
        if i == imax:
            break
        i = (i - 1) % n

    def yof(chain, x):
        for a, b in zip(chain, chain[1:]):
            xa, ya = a
            xb, yb = b
            lo, hi = (xa, xb) if xa <= xb else (xb, xa)
            if lo - mp.mpf(10) ** -35 <= x <= hi + mp.mpf(10) ** -35:
                if abs(xb - xa) < mp.mpf(10) ** -38:
                    return (ya + yb) / 2
                return ya + (yb - ya) * (x - xa) / (xb - xa)
        return None

    xmin = min(p[0] for p in Vr)
    xmax = max(p[0] for p in Vr)
    xm = (xmin + xmax) / 2
    yA, yB = yof(chainA, xm), yof(chainB, xm)
    upper, lower = (chainA, chainB) if (yA is not None and yB is not None and yA >= yB) else (chainB, chainA)

    def wid(x):
        u, l = yof(upper, x), yof(lower, x)
        if u is None or l is None:
            return mp.mpf(0)
        return abs(u - l)

    brk = sorted(set(p[0] for p in Vr))
    pieces = []
    for l, r in zip(brk, brk[1:]):
        if r - l < mp.mpf(10) ** -34:
            continue
        wl, wr = wid(l), wid(r)
        B = (wr - wl) / (r - l)
        A = wl - B * l
        pieces.append((l, r, A, B))
    return pieces


# ---- inner integrals via exact Gauss-Legendre on polynomial pieces ------------------------
def _slab(pieces, p, k):
    """int_t |t - p|^k w(t) dt over the whole projection."""
    tot = mp.mpf(0)
    for (l, r, A, B) in pieces:
        wf = lambda t: (A + B * t)
        if p <= l or p >= r:
            tot += _gl(lambda t: abs(t - p) ** k * wf(t), l, r)
        else:
            tot += _gl(lambda t: (p - t) ** k * wf(t), l, p)
            tot += _gl(lambda t: (t - p) ** k * wf(t), p, r)
    return tot


def EA3k_polygon(V, k, area):
    n = len(V)
    angles = {mp.mpf(0), mp.pi}
    # subdivide at every direction where two vertices project equally (edges AND diagonals):
    # the width-function breakpoint ordering changes there, giving kinks in the phi-integrand.
    for i in range(n):
        for j in range(i + 1, n):
            dx, dy = float(V[j][0] - V[i][0]), float(V[j][1] - V[i][1])
            a = math.atan2(dy, dx) % math.pi
            angles.add(mp.mpf(a))
            angles.add(mp.mpf((a + math.pi / 2) % math.pi))
    nodes = sorted(angles)
    coef = mp.mpf(2) ** (1 - k) / ((k + 2) * (k + 3))

    def phi_integrand(phi):
        pieces = _pieces(_rot(V, phi))
        tot = mp.mpf(0)
        for (l, r, A, B) in pieces:
            wf = lambda p: (A + B * p)
            tot += _gl(lambda p: wf(p) ** (k + 3) * _slab(pieces, p, k), l, r)
        return tot

    total = mp.mpf(0)
    for a, b in zip(nodes, nodes[1:]):
        if b - a < mp.mpf(10) ** -30:
            continue
        total += mp.quad(phi_integrand, [a, b])
    return coef * total / area ** (k + 3)


def E_cpow_polygon(V, area, jmax):
    """E[c^j] for j=0..jmax, c = area fraction on one side of the line through two random points.
    Via 2-point B-P:  E[g(c)] = (1/V^2) int_phi int_p (w(p,phi)^3/3) g(alpha(p,phi)) dp dphi,
    alpha(p,phi) = (cumulative width up to p)/V.  Returns list [E[c^0],...,E[c^jmax]]."""
    n = len(V)
    angles = {mp.mpf(0), mp.pi}
    for i in range(n):
        for j in range(i + 1, n):
            dx, dy = float(V[j][0] - V[i][0]), float(V[j][1] - V[i][1])
            a = math.atan2(dy, dx) % math.pi
            angles.add(mp.mpf(a))
            angles.add(mp.mpf((a + math.pi / 2) % math.pi))
    nodes = sorted(angles)

    def phi_integrand(phi, j):
        pieces = _pieces(_rot(V, phi))
        xmin = pieces[0][0]

        def cumw(p):  # int_{xmin}^{p} w(t) dt
            s = mp.mpf(0)
            for (l, r, A, B) in pieces:
                if p <= l:
                    break
                b = min(p, r)
                s += A * (b - l) + B * (b * b - l * l) / 2
                if p <= r:
                    break
            return s

        tot = mp.mpf(0)
        for (l, r, A, B) in pieces:
            wf = lambda p: (A + B * p)
            tot += _gl(lambda p: wf(p) ** 3 / 3 * (cumw(p) / area) ** j, l, r)
        return tot

    out = []
    for j in range(jmax + 1):
        total = mp.mpf(0)
        for a, b in zip(nodes, nodes[1:]):
            if b - a < mp.mpf(10) ** -30:
                continue
            total += mp.quad(lambda ph: phi_integrand(ph, j), [a, b])
        out.append(total / area ** 2)
    return out


def EA_hull_polygon(V, area, upto=5):
    """E[A_{n-1}] for n up to `upto`+1 via E[N_n] = C(n,2) E[(1-c)^{n-2}+c^{n-2}],
    E[A_{n-1}] = 1 - E[N_n]/n.  Returns dict {k: E[A_k]} for k=3..upto."""
    from math import comb
    Ec = E_cpow_polygon(V, area, upto - 1)  # need c^0..c^{upto-1}

    def Emix(power):  # E[(1-c)^power + c^power]
        s = mp.mpf(0)
        for i in range(power + 1):
            s += comb(power, i) * ((-1) ** i) * Ec[i]  # (1-c)^power
        return s + Ec[power]

    res = {}
    for n in range(4, upto + 2):
        EN = comb(n, 2) * Emix(n - 2)
        res[n - 1] = 1 - EN / n
    return res


def regular_mgon(m):
    V = [(mp.cos(2 * mp.pi * j / m), mp.sin(2 * mp.pi * j / m)) for j in range(m)]
    A = mp.fsum(V[i][0] * V[(i + 1) % m][1] - V[(i + 1) % m][0] * V[i][1] for i in range(m)) / 2
    return V, A


BODIES = {
    "triangle": ([(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)), (mp.mpf(0), mp.mpf(1))], mp.mpf(1) / 2),
    "square": ([(mp.mpf(0), mp.mpf(0)), (mp.mpf(1), mp.mpf(0)), (mp.mpf(1), mp.mpf(1)), (mp.mpf(0), mp.mpf(1))], mp.mpf(1)),
}
KNOWN = {("triangle", 1): mp.mpf(1) / 12, ("triangle", 2): mp.mpf(1) / 72,
         ("square", 1): mp.mpf(11) / 144, ("square", 2): mp.mpf(1) / 96}


def pslq_rational(x, maxden=10**9):
    r = mp.pslq([x, mp.mpf(1)], maxcoeff=maxden, maxsteps=10**5)
    if r and r[0] != 0:
        from fractions import Fraction
        fr = Fraction(-r[1], r[0])
        if abs(mp.mpf(fr.numerator) / fr.denominator - x) < mp.mpf(10) ** -25:
            return fr
    return None


if __name__ == "__main__":
    import sys, time
    for name in (sys.argv[1:] or ["triangle", "square"]):
        if name.startswith("m"):
            m = int(name[1:]); V, area = regular_mgon(m); label = f"regular {m}-gon"
        else:
            V, area = BODIES[name]; label = name
        print(f"=== {label} (area {mp.nstr(area, 10)}) ===")
        for k in (1, 2, 3):
            t0 = time.time()
            val = EA3k_polygon(V, k, area)
            line = f"  E[A_3^{k}] = {mp.nstr(val, 30)}  ({time.time()-t0:.1f}s)"
            if (name, k) in KNOWN:
                err = abs(val - KNOWN[(name, k)])
                line += f"   known {KNOWN[(name,k)]}  |err|={mp.nstr(err,3)}  {'OK' if err<mp.mpf(10)**-24 else 'FAIL'}"
            else:
                fr = pslq_rational(val)
                if fr is not None:
                    line += f"   = {fr}  [PSLQ]"
            print(line)
