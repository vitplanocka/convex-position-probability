r"""Closed-form search for E[A_5](regular m-gon) = 1 - 5 E[c^4].

Ansatz motivated by Alikoski's shape E[A_3] = (quadratic in C)/(36 m^2 S^2), C = cos w,
S = sin w, w = 2 pi/m, and by the disk limit E[A_5] -> 7(2400 pi^2 - 3289)/(6912 pi^4),
which has both a 1/pi^2 and a 1/pi^4 part (and m^2 S^2 -> 4 pi^2, m^4 S^4 -> 16 pi^4):

    E[A_5] = a_0 + P(C)/(m^2 S^2) + Q(C)/(m^4 S^4),   deg P <= 3, deg Q <= 4.

That is 1 + 4 + 5 = 10 unknowns.  Solve exactly (50-dps linear algebra) on 10 values of m,
then TEST on 8 held-out values of m.  A fit that reproduces held-out m to ~40 digits is
real; anything else is reported as a failure.  The basis is printed either way.
"""
import json
from mpmath import mp, mpf, matrix, lu_solve, cos, sin, pi, nstr
from polygon_exact import T_and_moments, regular_polygon_mp
mp.dps = 60

MS_FIT = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
MS_TEST = [13, 14, 15, 16, 20, 24, 30, 40]


def basis(m):
    w = 2 * pi / m
    C, S = cos(w), sin(w)
    d2 = m ** 2 * S ** 2
    d4 = m ** 4 * S ** 4
    b = [mpf(1)]
    b += [C ** k / d2 for k in range(4)]
    b += [C ** k / d4 for k in range(5)]
    return b


def EA5(m):
    return T_and_moments(regular_polygon_mp(m, 60), nmax=6, dps=60, nodes=9)["E[A_5]"]


vals = {m: EA5(m) for m in MS_FIT + MS_TEST}
A = matrix([basis(m) for m in MS_FIT])
rhs = matrix([vals[m] for m in MS_FIT])
x = lu_solve(A, rhs)
names = ["1"] + [f"C^{k}/(m^2 S^2)" for k in range(4)] + [f"C^{k}/(m^4 S^4)" for k in range(5)]
print("basis tried:", names)
print("fitted coefficients (50 dps):")
for n, c in zip(names, x):
    print(f"   {n:18s} {nstr(c, 30)}")
print("\nheld-out test (must be ~1e-45 if the ansatz is right):")
worst = mpf(0)
for m in MS_TEST:
    pred = sum(bi * xi for bi, xi in zip(basis(m), x))
    err = pred - vals[m]
    worst = max(worst, abs(err))
    print(f"   m={m:3d}  E[A_5]={nstr(vals[m],30)}  residual={nstr(err,5)}")
print(f"\nworst held-out residual: {nstr(worst,5)}")
print("VERDICT:", "ANSATZ FITS" if worst < mpf(10) ** -40 else "ANSATZ FAILS -- no closed form of this shape")
json.dump({"basis": names, "coeffs": [nstr(c, 40) for c in x],
           "worst_heldout_residual": nstr(worst, 5),
           "verdict": "fits" if worst < mpf(10) ** -40 else "fails"},
          open("../results/C_mgon_EA5_fit.json", "w"), indent=1)
