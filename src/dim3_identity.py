r"""A 3-D analogue of identity (I), and a sharp test of it.

In R^2 the proof of (I) used #vertices = #edges of the hull.  In R^3 that is false, but for a
SIMPLICIAL 3-polytope Euler's formula V - E + F = 2 together with 2E = 3F gives exactly

    V = 2 + F/2,

and generic points give a simplicial hull.  The Renyi-Sulanke facet count in R^3 is
    E[F_n] = C(n,3) E[ c^{n-3} + (1-c)^{n-3} ] = 2 C(n,3) E[c^{n-3}],
with c the mu-mass on one side of the plane through three of the points; c and 1-c are again
identically distributed (an odd permutation of the three points flips the orientation), so
E[c] = 1/2 and E[c^3] = (3E[c^2] - 1/2)/2 exactly as in the plane.  Efron
E[A_{n-1}] = 1 - E[N_n]/n then gives, with m_k := E[c^k],

    E[A_{n-1}] = 1 - 2/n - C(n,3) m_{n-3} / n,
    n=4:  E[A_3] = 1/2 - m_1 = 0                      (correct: 3 points span no volume)
    n=5:  E[A_4] = 3/5 - 2 m_2
    n=6:  E[A_5] = 2/3 - (10/3) m_3 = 3/2 - 5 m_2 = (5/2) E[A_4].

    ==>  E[A_5] = (5/2) E[A_4]  in R^3, for every absolutely continuous law.       (I_3)

n=7 gives E[A_6] = 5/7 - 5 m_4, a genuinely new moment -- same "every other step" pattern.

TEST (this file).  Instead of hull volumes we use Efron in reverse: X_i is a NON-vertex of
conv(X_1..X_n) iff X_i lies in the hull of the others, and by Caratheodory in R^3 that means
inside conv of some FOUR of them.  So N_5 and N_6 are exact combinatorial indicators
(orientation tests only, no floating hull construction), and
    E[A_4] = 1 - E[N_5]/5,   E[A_5] = 1 - E[N_6]/6.
Anchors: E[A_4](ball) = 9/715, E[A_4](cube) = 3977/216000 - pi^2/2160.
"""
import json, math, sys
from itertools import combinations
import numpy as np
from convex_position import _orient3, sample_ball, sample_cube, sample_simplex3

BODIES = {"ball": sample_ball, "cube": sample_cube, "simplex": sample_simplex3}
EXACT_EA4 = {"ball": 9 / 715, "cube": 3977 / 216000 - math.pi ** 2 / 2160}


def in_tetra(A, B, C, D, X):
    s0 = _orient3(A, B, C, D)
    s1 = _orient3(X, B, C, D)
    s2 = _orient3(A, X, C, D)
    s3 = _orient3(A, B, X, D)
    s4 = _orient3(A, B, C, X)
    return (s1 * s0 >= 0) & (s2 * s0 >= 0) & (s3 * s0 >= 0) & (s4 * s0 >= 0)


def n_vertices(P):
    """P: (m, n, 3) -> number of hull vertices per sample (n = 5 or 6), exactly."""
    m, n, _ = P.shape
    nonvert = np.zeros(m, dtype=np.int64)
    for i in range(n):
        others = [j for j in range(n) if j != i]
        inside = np.zeros(m, dtype=bool)
        for quad in combinations(others, 4):
            inside |= in_tetra(P[:, quad[0]], P[:, quad[1]], P[:, quad[2]], P[:, quad[3]], P[:, i])
        nonvert += inside
    return n - nonvert


if __name__ == "__main__":
    N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 40_000_000
    B = 25
    rng = np.random.default_rng(2718)
    out = {}
    for name, samp in BODIES.items():
        per = N // B
        b4 = np.zeros(B); b5 = np.zeros(B)
        for i in range(B):
            done = 0; s4 = 0.0; s5 = 0.0
            while done < per:
                mm = min(500_000, per - done)
                P = samp(rng, mm, 6)
                s5 += n_vertices(P[:, :5]).sum()
                s4 += n_vertices(P).sum()
                done += mm
            b4[i] = 1 - (s5 / per) / 5          # E[A_4] = 1 - E[N_5]/5
            b5[i] = 1 - (s4 / per) / 6          # E[A_5] = 1 - E[N_6]/6
        f = math.sqrt(B)
        d = b5 - 2.5 * b4
        row = dict(EA4=b4.mean(), se4=b4.std(ddof=1) / f, EA5=b5.mean(), se5=b5.std(ddof=1) / f,
                   ratio=b5.mean() / b4.mean(), diff=d.mean(), se_diff=d.std(ddof=1) / f,
                   z=d.mean() / (d.std(ddof=1) / f), samples=N)
        if name in EXACT_EA4:
            row["EA4_exact"] = EXACT_EA4[name]
            row["z_EA4"] = (b4.mean() - EXACT_EA4[name]) / row["se4"]
        out[name] = row
        anchor = (f"   anchor E[A_4]={row['EA4_exact']:.9f} z={row['z_EA4']:+.2f}"
                  if "EA4_exact" in row else "")
        print(f"{name:8s} E[A_4]={row['EA4']:.9f}+-{row['se4']:.1e}{anchor}")
        print(f"{'':8s} E[A_5]={row['EA5']:.9f}+-{row['se5']:.1e}   ratio={row['ratio']:.7f}")
        print(f"{'':8s} (I_3) E[A_5]-(5/2)E[A_4] = {row['diff']:+.3e} +- {row['se_diff']:.1e}"
              f"   z = {row['z']:+.2f}", flush=True)
    json.dump(out, open("../results/dim3_identity.json", "w"), indent=1, default=float)
    print("\nwrote ../results/dim3_identity.json")
