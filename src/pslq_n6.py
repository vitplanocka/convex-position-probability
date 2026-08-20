"""PSLQ the regular-m-gon n=6 ingredients in their natural field Q(cos 2pi/m).

deg [Q(cos 2pi/m):Q] = phi(m)/2 (=1 for m=3,4,6), so the basis is {1, c, ..., c^(d-1)}.
Two strategies: (a) plain PSLQ; (b) a SMOOTH-DENOMINATOR scan -- for each D = 2^a 3^b 5^c 7^d
below a bound, test whether D*x is an integer combination of the basis with small coefficients.
(b) buys ~log10(D) digits and is what the observed denominators (4500, 54000, 34992000,
699840, ...) call for."""
import itertools, json, os, sys
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 40


def field_deg(m):
    from math import gcd
    phi = sum(1 for k in range(1, m + 1) if gcd(k, m) == 1)
    return max(1, phi // 2)


def cval(m):
    return mp.cos(2 * mp.pi / m)


def plain(x, m, maxcoeff=10 ** 8, dtol=16):
    d = field_deg(m)
    c = cval(m)
    v = [x] + [c ** j for j in range(d)]
    r = mp.pslq(v, maxcoeff=maxcoeff, maxsteps=10 ** 6, tol=mp.mpf(10) ** -dtol)
    if not r or not r[0]:
        return None
    den = r[0]
    coef = [Fraction(-r[1 + j], den) for j in range(d)]
    val = mp.fsum([mp.mpf(f.numerator) / f.denominator * c ** j for j, f in enumerate(coef)])
    return coef, abs(val - x)


def smooth(limit=10 ** 10, primes=(2, 3, 5, 7)):
    out = [1]
    for p in primes:
        nxt = []
        for v in out:
            w = v
            while w <= limit:
                nxt.append(w)
                w *= p
        out = nxt
    return sorted(set(v for v in out if v <= limit))


def scan(x, m, limit=10 ** 11, maxc=10 ** 6, relacc=mp.mpf(10) ** -16):
    """Find D smooth with D*x = sum a_j c^j, a_j integers (|a_j| <= maxc)."""
    d = field_deg(m)
    c = cval(m)
    hits = []
    for D in smooth(limit):
        y = x * D
        tol = abs(y) * relacc * 30
        if d == 1:
            n = mp.nint(y)
            if abs(y - n) < tol and abs(n) > 0:
                hits.append((D, [Fraction(int(n), D)], abs(y - n)))
        else:
            v = [y] + [c ** j for j in range(d)]
            r = mp.pslq(v, maxcoeff=maxc, maxsteps=200000,
                        tol=max(tol / abs(y), mp.mpf(10) ** -18))
            if r and r[0] == -1 or (r and abs(r[0]) == 1):
                sgn = r[0]
                coef = [Fraction(-r[1 + j] * sgn, D) for j in range(d)]
                val = mp.fsum([mp.mpf(f.numerator) / f.denominator * c ** j
                               for j, f in enumerate(coef)])
                if abs(val - x) < abs(x) * relacc * 30:
                    hits.append((D, coef, abs(val - x)))
    return hits


def fmt(coef):
    ts = []
    for j, f in enumerate(coef):
        if f == 0:
            continue
        p = "" if j == 0 else (" c" if j == 1 else f" c^{j}")
        ts.append(f"({f}){p}")
    return " + ".join(ts) if ts else "0"
