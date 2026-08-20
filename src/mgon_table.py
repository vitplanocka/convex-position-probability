r"""B1/B2: exact P_4, E[A_3], E[A_3^2], P_5 for regular m-gons.

Two independent ingredients, each checked against the other:
  (i)  ALIKOSKI 1939 (via MathWorld "Polygon Triangle Picking"): for a regular
       m-gon of unit area, the mean random-triangle area is
           E[A_3] = (9 cos^2 w + 52 cos w + 44) / (36 m^2 sin^2 w),  w = 2 pi/m.
  (ii) ROUTE P (this session, src/polygon_exact.py): the same quantity computed
       from scratch as E[A_3] = 1 - (3/2) T_2 with T_2 the Renyi-Sulanke line
       integral, evaluated EXACTLY (polynomial integrand, Gauss-Legendre) at 50 dps.
Then, with the exact second moment (this session)
       E[A_3^2] = (3/2) det Sigma / |K|^2 = (2 + cos w)^2 / (24 m^2 sin^2 w)
identity (II) P_5 = 1 - 10 E[A_3] + 10 E[A_3^2] gives the closed form
       P_5(m-gon) = 1 - 5 (15 cos^2 w + 92 cos w + 76) / (36 m^2 sin^2 w).
Anchors it must hit: m=3 -> 11/36, m=4 -> 49/144 (Valtr), m->oo -> 1 - 305/(48 pi^2).
"""
import json, math
import sympy as sp
from mpmath import mp, mpf, cos as mcos, sin as msin, pi as mpi, nstr
from polygon_exact import T_and_moments, regular_polygon_mp, EA3sq, regular_polygon

mp.dps = 50


def alikoski_EA3(m, symbolic=False):
    if symbolic:
        w = 2 * sp.pi / m
        c, s = sp.cos(w), sp.sin(w)
        return (9 * c ** 2 + 52 * c + 44) / (36 * m ** 2 * s ** 2)
    w = 2 * mpi / m
    c, s = mcos(w), msin(w)
    return (9 * c ** 2 + 52 * c + 44) / (36 * m ** 2 * s ** 2)


def EA3sq_closed(m, symbolic=False):
    if symbolic:
        w = 2 * sp.pi / m
        return (2 + sp.cos(w)) ** 2 / (24 * m ** 2 * sp.sin(w) ** 2)
    w = 2 * mpi / m
    return (2 + mcos(w)) ** 2 / (24 * m ** 2 * msin(w) ** 2)


def P5_closed(m, symbolic=False):
    if symbolic:
        w = 2 * sp.pi / m
        c, s = sp.cos(w), sp.sin(w)
        return 1 - 5 * (15 * c ** 2 + 92 * c + 76) / (36 * m ** 2 * s ** 2)
    w = 2 * mpi / m
    c, s = mcos(w), msin(w)
    return 1 - 5 * (15 * c ** 2 + 92 * c + 76) / (36 * m ** 2 * s ** 2)


def symbolic_derivation():
    print("=== symbolic derivation of the P_5 closed form from (II) + Alikoski ===")
    m = sp.Symbol('m', positive=True)
    w = sp.Symbol('w')
    c, s = sp.cos(w), sp.sin(w)
    EA3 = (9 * c ** 2 + 52 * c + 44) / (36 * m ** 2 * s ** 2)
    E2 = (2 + c) ** 2 / (24 * m ** 2 * s ** 2)
    P5 = sp.simplify(1 - 10 * EA3 + 10 * E2)
    claim = 1 - 5 * (15 * c ** 2 + 92 * c + 76) / (36 * m ** 2 * s ** 2)
    print("  1 - 10E[A_3] + 10E[A_3^2] =", sp.simplify(P5))
    print("  minus claimed closed form  =", sp.simplify(P5 - claim), " (must be 0)")
    # limit m -> oo  (w = 2 pi / m)
    P5m = P5.subs(w, 2 * sp.pi / m)
    lim = sp.limit(P5m, m, sp.oo)
    print("  limit m->oo:", sp.simplify(lim), " vs 1-305/(48 pi^2):",
          sp.simplify(lim - (1 - sp.Rational(305, 48) / sp.pi ** 2)))
    for mm, ref, name in [(3, sp.Rational(11, 36), "Valtr triangle"), (4, sp.Rational(49, 144), "Valtr square")]:
        val = sp.nsimplify(sp.simplify(P5m.subs(m, mm)))
        print(f"  m={mm}: P_5 = {val}   ref {ref} ({name})   diff {sp.simplify(val-ref)}")
    # E[A_3^2] closed form vs the covariance route, symbolically for a few m
    for mm in [3, 4, 5, 6, 8, 12]:
        V = regular_polygon(mm, exact=True)
        lhs = sp.simplify(EA3sq(V))
        rhs = sp.simplify(EA3sq_closed(mm, symbolic=True))
        print(f"  m={mm}: E[A_3^2] cov-route - closed form = {sp.simplify(lhs - rhs)}  (must be 0), value {sp.nsimplify(sp.simplify(lhs))}")
    # E[A_3] Alikoski vs Valtr/Sylvester anchors
    for mm, ref in [(3, sp.Rational(1, 12)), (4, sp.Rational(11, 144))]:
        print(f"  m={mm}: Alikoski E[A_3] = {sp.nsimplify(sp.simplify(alikoski_EA3(mm, True)))}  ref {ref}")


def numeric_table(ms):
    rows = []
    print("\n=== route P (50 dps, exact polynomial integration) vs Alikoski ===")
    for m in ms:
        r = T_and_moments(regular_polygon_mp(m, 50), nmax=5, dps=50)
        ea3_P = r['E[A_3]']
        ea4_P = r['E[A_4]']
        ea3_A = alikoski_EA3(m)
        e2 = EA3sq_closed(m)
        P4 = 1 - 4 * ea3_A
        P5 = 1 - 10 * ea3_A + 10 * e2
        row = dict(m=m,
                   EA3_routeP=nstr(ea3_P, 40), EA3_alikoski=nstr(ea3_A, 40),
                   diff_EA3=nstr(ea3_P - ea3_A, 5),
                   identityI_resid=nstr(ea4_P - 2 * ea3_P, 5),
                   T0_minus_2=nstr(r['T'][0] - 2, 5), EN3_minus_3=nstr(r['E[N_3]'] - 3, 5),
                   P4=nstr(P4, 30), EA3sq=nstr(e2, 30), P5=nstr(P5, 30),
                   P5_closed_resid=nstr(P5 - P5_closed(m), 5))
        rows.append(row)
        print(f" m={m:3d} E[A_3](routeP)={nstr(ea3_P,30)}")
        print(f"       E[A_3](Alikoski)={nstr(ea3_A,30)}   diff={nstr(ea3_P-ea3_A,5)}")
        print(f"       E[A_4]-2E[A_3] = {nstr(ea4_P-2*ea3_P,5)}   T_0-2={nstr(r['T'][0]-2,5)}  E[N_3]-3={nstr(r['E[N_3]']-3,5)}")
        print(f"       P_4={nstr(P4,25)}  E[A_3^2]={nstr(e2,25)}  P_5={nstr(P5,25)}")
    return rows


if __name__ == "__main__":
    symbolic_derivation()
    ms = list(range(3, 25)) + [30, 40, 60, 100]
    rows = numeric_table(ms)
    json.dump(rows, open("../results/B1_mgon_exact_table.json", "w"), indent=1)
    print("\nwrote ../results/B1_mgon_exact_table.json")
