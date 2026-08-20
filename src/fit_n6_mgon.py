"""Closed-form fits for the n=6 regular-m-gon ingredients.

Natural variables (circumradius-1 regular m-gon):  w = 2pi/m, c = cos w, s = sin w,
V = (m/2) s, so  V^2 = S/4  with  S = m^2 s^2.  (S -> 4 pi^2 as m -> infinity: the disk.)
Ansatz for every n=6 ingredient X:   X = A(c)/S + B(c)/S^2  with A, B in Q[c].
(E[A_3] = (9c^2+52c+44)/(36 S) is the known n=4 instance -- Alikoski.)
"""
import json, os, sys
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 40


def basis(m, dA, dB):
    w = 2 * mp.pi / m
    c, s = mp.cos(w), mp.sin(w)
    S = m ** 2 * s ** 2
    row = [c ** j / S for j in range(dA + 1)] + [c ** j / S ** 2 for j in range(dB + 1)]
    return row


def fit(ms, ys, dA, dB, use=None):
    """Solve the (square or overdetermined) system; return coefficients + max residual."""
    n = (dA + 1) + (dB + 1)
    use = use or list(range(n))
    A = mp.matrix([basis(ms[i], dA, dB) for i in use])
    b = mp.matrix([ys[i] for i in use])
    x = mp.lu_solve(A, b)
    res = mp.mpf(0)
    for i in range(len(ms)):
        r = mp.fdot(basis(ms[i], dA, dB), x) - ys[i]
        res = max(res, abs(r) / abs(ys[i]))
    return x, res


def rationalise(x, maxden=10 ** 9, tol=mp.mpf(10) ** -18):
    out = []
    for v in x:
        fr = Fraction(float(v)).limit_denominator(maxden)
        ok = abs(mp.mpf(fr.numerator) / fr.denominator - v) < tol * max(1, abs(v))
        out.append((fr, ok))
    return out


def show(name, x, dA, dB, res, maxden=10 ** 9, tol=mp.mpf(10) ** -18):
    rt = rationalise(x, maxden, tol)
    allok = all(o for _, o in rt)
    A = [str(f) for f, _ in rt[:dA + 1]]
    B = [str(f) for f, _ in rt[dA + 1:]]
    print(f"  {name}: dA={dA} dB={dB}  maxrelres={mp.nstr(res,3)}  "
          f"rational={'YES' if allok else 'no'}")
    print(f"     A(c) coeffs (c^0..c^{dA}) = {A}")
    print(f"     B(c) coeffs (c^0..c^{dB}) = {B}")
    return allok


if __name__ == "__main__":
    J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results",
                                    "n6_twochord_polygon.json")))
    ms = [r["m"] for r in J["rows"]]
    data = {k: [mp.mpf(r[k]) for r in J["rows"]] for k in ("EA33", "EA5", "EA42_convex", "P6")}
    for key, tol in (("EA33", mp.mpf(10) ** -20), ("EA5", mp.mpf(10) ** -20),
                     ("EA42_convex", mp.mpf(10) ** -11), ("P6", mp.mpf(10) ** -10)):
        print(f"=== {key} ===")
        for dA, dB in [(-1, 2), (-1, 3), (-1, 4), (-1, 5), (0, 4), (1, 4), (2, 4),
                       (2, 2), (2, 3), (2, 5), (3, 4)]:
            n = (dA + 1) + (dB + 1)
            if n > len(ms):
                continue
            try:
                x, res = fit(ms, data[key], dA, dB)
            except Exception as e:
                continue
            if res < tol * 100:
                if show(f"{key}", x, dA, dB, res, tol=tol):
                    break
