"""P_6 for regular m-gons, ASSEMBLED FROM FIRST PRINCIPLES.

    P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3],   E[A_4^2] = E[A_4^2 & convex] + 4 E[A_3^3]
        = 1 - 6 E[A_5] + 15 E[A_4^2 & convex] + 40 E[A_3^3].

Ingredients, all computed here, none taken from Valtr/Marckert:
  E[A_3^3]          -- Blaschke-Petkantschin WIDTH integral      (n6_bp_polygon.EA3k_polygon)
  E[A_5]            -- E[A_5] = 1 - E[N_6]/6, E[N_6] = 15 E[(1-c)^4+c^4]  (EA_hull_polygon)
  E[A_4^2 & convex] -- two-chord Blaschke-Petkantschin integral   (n6_twochord_polygon)

Anchors reproduced by this exact pipeline: m=3 -> 91/900, m=4 -> 49/400 (Valtr), m -> infinity
-> the disk value 0.1343093863571 (Marckert).
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import n6_bp_polygon as BP
import n6_twochord_polygon as TC

mp.mp.dps = 40


def ingredients(m, nphi=32, workers=1):
    V, area = BP.regular_mgon(m)
    t0 = time.time()
    EA33 = BP.EA3k_polygon(V, 3, area)
    t1 = time.time()
    EA5 = BP.EA_hull_polygon(V, area, 5)[5]
    t2 = time.time()
    Vf = TC.regular_mgon(m)
    q = TC.qsym_of(m)
    _, EA42c_a = TC.two_chord_polygon(Vf, nphi, nphi, 6, q, workers=workers)
    _, EA42c_b = TC.two_chord_polygon(Vf, nphi + 8, nphi + 8, 6, q, workers=workers)
    t3 = time.time()
    EA42c = mp.mpf(float(EA42c_b))
    P6 = 1 - 6 * EA5 + 15 * EA42c + 40 * EA33
    return dict(m=m, EA33=EA33, EA5=EA5, EA42c=EA42c,
                EA42=EA42c + 4 * EA33, P6=P6,
                tc_spread=abs(float(EA42c_a) - float(EA42c_b)),
                secs=(t1 - t0, t2 - t1, t3 - t2))


def _job(args):
    m, nphi = args
    return ingredients(m, nphi, workers=1)


if __name__ == "__main__":
    import multiprocessing as mp_
    ms = [int(x) for x in (sys.argv[1:] or "3 4 5 6 7 8 9 10 11 12".split())]
    nphi = int(os.environ.get("NPHI", "32"))
    t0 = time.time()
    with mp_.Pool(min(len(ms), 10)) as pool:
        res = pool.map(_job, [(m, nphi) for m in ms])
    print(f"total {time.time()-t0:.0f}s\n")
    hdr = f"{'m':>3} {'E[A_3^3]':>24} {'E[A_5]':>24} {'E[A_4^2&conv]':>24} {'P_6':>24}"
    print(hdr)
    rows = []
    for r in res:
        print(f"{r['m']:>3} {mp.nstr(r['EA33'],18):>24} {mp.nstr(r['EA5'],18):>24} "
              f"{mp.nstr(r['EA42c'],18):>24} {mp.nstr(r['P6'],18):>24}")
        rows.append(dict(m=r['m'], EA33=mp.nstr(r['EA33'], 30), EA5=mp.nstr(r['EA5'], 30),
                         EA42_convex=mp.nstr(r['EA42c'], 20), EA42=mp.nstr(r['EA42'], 20),
                         P6=mp.nstr(r['P6'], 20), P6_float=float(r['P6']),
                         twochord_grid_spread=float(r['tc_spread']),
                         secs=[round(x, 1) for x in r['secs']]))
    out = os.path.join(os.path.dirname(__file__), "..", "results", "n6_twochord_polygon.json")
    json.dump(dict(nphi=nphi, rows=rows,
                   formula="P_6 = 1 - 6 E[A_5] + 15 E[A_4^2 & convex] + 40 E[A_3^3]",
                   note="all three ingredients computed here; no Valtr/Marckert input"),
              open(out, "w"), indent=1)
    print("\nwrote", os.path.normpath(out))
