"""Exact and 50-digit E[A_3], P_4, E[A_3^2], P_5 for regular m-gons (route P)."""
import sys, json, math, time
import sympy as sp
from polygon_exact import (J_integrals, T_and_moments, regular_polygon, regular_polygon_mp,
                           EA3sq, cov_det_area)


def exact_mgon(m, nmax=5):
    V = regular_polygon(m, exact=True)
    r = T_and_moments(V, nmax=nmax, exact=True)
    EA3 = sp.nsimplify(sp.simplify(sp.radsimp(sp.expand(r['E[A_3]']))))
    EA4 = sp.simplify(sp.expand(r['E[A_4]']))
    return EA3, EA4, r


def mp_mgon(m, nmax=5, dps=50, nodes=None):
    from mpmath import mp
    mp.dps = dps
    V = regular_polygon_mp(m, dps)
    return T_and_moments(V, nmax=nmax, dps=dps, nodes=nodes)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "exact"
    ms = [int(x) for x in sys.argv[2:]] or [3, 4, 5, 6, 8, 10, 12]
    if mode == "exact":
        for m in ms:
            t0 = time.time()
            EA3, EA4, r = exact_mgon(m)
            chk = sp.simplify(EA4 - 2 * EA3)
            P4 = sp.simplify(1 - 4 * EA3)
            V = regular_polygon(m, exact=True)
            e2 = sp.simplify(EA3sq(V))
            P5 = sp.simplify(1 - 10 * EA3 + 10 * e2)
            print(f"m={m}")
            print(f"  E[A_3]  = {sp.simplify(EA3)}   = {sp.N(EA3, 30)}")
            print(f"  E[A_4]-2E[A_3] = {chk}   (identity I, must be 0)")
            print(f"  P_4     = {sp.nsimplify(P4)}   = {sp.N(P4, 30)}")
            print(f"  E[A_3^2]= {sp.simplify(e2)}   = {sp.N(e2, 30)}")
            print(f"  P_5     = {sp.nsimplify(sp.simplify(P5))}   = {sp.N(P5, 30)}")
            print(f"  ({time.time()-t0:.1f}s)", flush=True)
    else:
        from mpmath import mp, nstr
        for m in ms:
            t0 = time.time()
            r = mp_mgon(m, dps=50)
            print(f"m={m:3d}  E[A_3] = {nstr(r['E[A_3]'], 40)}   T_0={nstr(r['T'][0],20)}  "
                  f"E[N_3]={nstr(r['E[N_3]'],20)}  ({time.time()-t0:.1f}s)", flush=True)
