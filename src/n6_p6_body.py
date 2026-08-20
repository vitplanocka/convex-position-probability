"""P_6(K) for an ARBITRARY convex body, assembled from first principles, to attack the OPEN
n >= 6 extremal conjecture  P_6(triangle) <= P_6(K) <= P_6(disk)  on non-regular bodies.

    P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3],   E[A_4^2] = E[A_4^2&conv] + 4 E[A_3^3]
        = 1 - 6 E[A_5] + 15 E[A_4^2&conv] + 40 E[A_3^3].

Ingredients (all first-principles, no Valtr/Marckert):
  E[A_3^3], E[A_5] : Blaschke-Petkantschin WIDTH route  (n6_bp_polygon: EA3k_polygon, EA_hull_polygon)
  E[A_4^2&conv]    : TWO-CHORD route                    (n6_twochord_polygon: two_chord_polygon)

Window (exact): P_6(triangle) = 91/900 = 0.1011111..., P_6(disk) = 0.1343093864...
Affine invariance => all triangles give 91/900, all parallelograms 49/400; a meaningful test needs
bodies NOT affine-equivalent to a regular one (general trapezoids, kites, pushed polygons, ...).
"""
import math
import time
import mpmath as mp
import numpy as np

import n6_bp_polygon as W          # width route (mpmath)
import n6_twochord_polygon as TC   # two-chord route (float64/longdouble)

mp.mp.dps = 30
P6_TRI = mp.mpf(91) / 900
P6_DISK = mp.mpf(1) - (146400 * mp.pi**2 - 473473) / (11520 * mp.pi**4)


def shoelace(V):
    n = len(V)
    return abs(sum(V[i][0] * V[(i + 1) % n][1] - V[(i + 1) % n][0] * V[i][1] for i in range(n))) / 2


def p6_body(V, nphi=32, dtname="ld", ng=6):
    """V: list of (x,y) floats, convex, ccw. Returns dict with P_6 and the ingredients."""
    Vm = [(mp.mpf(str(x)), mp.mpf(str(y))) for (x, y) in V]
    area_m = mp.mpf(str(shoelace(V)))
    t0 = time.time()
    EA33 = W.EA3k_polygon(Vm, 3, area_m)
    EA = W.EA_hull_polygon(Vm, area_m, upto=5)
    EA5 = EA[5]
    tw = time.time() - t0
    _, EA42c_ld = TC.two_chord_polygon(V, nphi, nphi, ng, 1, workers=1, dtname=dtname)
    EA42c = mp.mpf(str(EA42c_ld))
    ttc = time.time() - t0 - tw
    P6 = 1 - 6 * EA5 + 15 * EA42c + 40 * EA33
    return dict(P6=P6, EA5=EA5, EA33=EA33, EA42_conv=EA42c,
                below_tri=bool(P6 < P6_TRI - mp.mpf(10)**-9),
                above_disk=bool(P6 > P6_DISK + mp.mpf(10)**-9),
                margin_low=mp.nstr(P6 - P6_TRI, 6), margin_high=mp.nstr(P6_DISK - P6, 6),
                t_width=round(tw, 1), t_twochord=round(ttc, 1))


# --- a starter family of NON-regular convex bodies (none affine-equivalent to a regular one) ---
def half_disk(n=48):
    pts = [(math.cos(math.pi * k / n), math.sin(math.pi * k / n)) for k in range(n + 1)]
    return pts  # semicircle arc + the diameter closes it (ccw)


def reg_pushed(m, idx=0, factor=1.4):
    V = [[math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m)] for k in range(m)]
    V[idx][0] *= factor; V[idx][1] *= factor
    return [tuple(p) for p in V]


BODIES = {
    "triangle(right)": [(0., 0.), (1., 0.), (0., 1.)],                 # control: must give 91/900
    "square(unit)":    [(0., 0.), (1., 0.), (1., 1.), (0., 1.)],       # control: must give 49/400
    "trapezoid":       [(0., 0.), (2., 0.), (1.3, 1.), (0.2, 1.)],      # general (non-parallelogram)
    "kite":            [(0., 0.), (1., 0.8), (0.3, 2.0), (-0.6, 0.8)],  # convex kite
    "pentagon(pushed)": reg_pushed(5, 0, 1.5),                          # irregular pentagon
    "half-disk":       half_disk(48),                                   # smooth, non-regular
    "pgon(rounded-tri)": reg_pushed(3, 0, 1.0),                          # equilateral control (=91/900)
}


if __name__ == "__main__":
    import sys, json
    print(f"window: P_6(triangle)={mp.nstr(P6_TRI,12)}  P_6(disk)={mp.nstr(P6_DISK,12)}\n")
    names = sys.argv[1:] or list(BODIES)
    rows = []
    for name in names:
        V = BODIES[name]
        r = p6_body(V)
        flag = "  <-- BELOW TRIANGLE!" if r["below_tri"] else ("  <-- ABOVE DISK!" if r["above_disk"] else "")
        print(f"{name:20s} P_6 = {mp.nstr(r['P6'],12)}   (tri+{r['margin_low']}, disk-{r['margin_high']})"
              f"  [{r['t_width']}+{r['t_twochord']}s]{flag}")
        rows.append(dict(body=name, P6=mp.nstr(r["P6"], 15), below_tri=r["below_tri"],
                         above_disk=r["above_disk"]))
    json.dump({"window": [mp.nstr(P6_TRI, 15), mp.nstr(P6_DISK, 15)], "rows": rows},
              open("../results/n6_nonregular_probe.json", "w"), indent=1)
