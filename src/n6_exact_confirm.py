"""Confirm the direct-MC scan against the EXACT first-principles P_6 machinery on NON-REGULAR
bodies (the exact route was previously validated only on triangles/squares/regular m-gons).

    P_6 = 1 - 6 E[A_5] + 15 E[A_4^2&conv] + 40 E[A_3^3]
      E[A_5], E[A_3^3]  : Blaschke-Petkantschin WIDTH route      (n6_bp_polygon, mpmath)
      E[A_4^2&conv]     : TWO-CHORD line-space integral, 80-bit  (n6_twochord_polygon)

Both routes are structurally independent of the Monte Carlo engine (mcp6.py).
"""
from __future__ import annotations

import argparse
import json
import math
import time

import mpmath as mp
import numpy as np

import mcp6
import n6_bodies as NB
import n6_bp_polygon as W
import n6_twochord_polygon as TC

mp.mp.dps = 30
P6_TRI = mp.mpf(91) / 900
P6_DISK = mp.mpf(1) - (146400 * mp.pi**2 - 473473) / (11520 * mp.pi**4)


def shoelace(V):
    n = len(V)
    return abs(sum(V[i][0] * V[(i + 1) % n][1] - V[(i + 1) % n][0] * V[i][1]
                   for i in range(n))) / 2


def p6_exact(V, nphi=32, ng=6, dtname="ld", workers=14):
    V = [(float(x), float(y)) for x, y in V]
    Vm = [(mp.mpf(repr(x)), mp.mpf(repr(y))) for (x, y) in V]
    area_m = mp.mpf(repr(shoelace(V)))
    t0 = time.time()
    EA33 = W.EA3k_polygon(Vm, 3, area_m)
    EA5 = W.EA_hull_polygon(Vm, area_m, upto=5)[5]
    tw = time.time() - t0
    _, EA42c = TC.two_chord_polygon(V, nphi, nphi, ng, 1, workers=workers, dtname=dtname)
    ttc = time.time() - t0 - tw
    P6 = 1 - 6 * EA5 + 15 * mp.mpf(str(EA42c)) + 40 * EA33
    return dict(P6=P6, EA5=EA5, EA33=EA33, EA42_conv=mp.mpf(str(EA42c)),
                t_width=round(tw, 1), t_twochord=round(ttc, 1))


DEFAULT = ["trapezoid-0.4", "kite-1.0-0.8", "trunc1-0.05", "trunc3-0.2",
           "push5-1.6", "randquad-0", "randpent-0"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("bodies", nargs="*", default=DEFAULT)
    ap.add_argument("--samples", type=float, default=4e9)
    ap.add_argument("--nphi", type=int, default=32)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--out", default="../results/n6_exact_confirm.json")
    a = ap.parse_args()
    cat = {n: V for n, f, V in NB.catalogue()}
    rows = []
    print(f"window: 91/900 = {float(P6_TRI):.9f}   disk = {float(P6_DISK):.9f}\n")
    print(f"{'body':16s} {'nv':>3s} {'P_6 exact':>20s} {'P_6 MC':>13s} {'+-':>9s} {'z':>7s}")
    for name in a.bodies:
        V = cat[name]
        ex = p6_exact(V, a.nphi, workers=a.workers)
        r = mcp6.mc(V, int(a.samples), seed=31337, nthreads=16)
        z = (r["P6"] - float(ex["P6"])) / r["se6"]
        ok = abs(z) < 4
        print(f"{name:16s} {len(V):3d} {mp.nstr(ex['P6'], 17):>20s} {r['P6']:13.9f} "
              f"{r['se6']:9.1e} {z:+7.2f}  {'OK' if ok else 'FAIL'}"
              f"   [{ex['t_width']}+{ex['t_twochord']}s]", flush=True)
        rows.append(dict(body=name, nv=int(len(V)), P6_exact=mp.nstr(ex["P6"], 20),
                         P6_exact_float=float(ex["P6"]),
                         EA5=mp.nstr(ex["EA5"], 20), EA33=mp.nstr(ex["EA33"], 20),
                         EA42_conv=mp.nstr(ex["EA42_conv"], 20),
                         P6_mc=r["P6"], se_mc=r["se6"], samples=r["samples"], z=z, ok=bool(ok),
                         gap_to_triangle=float(ex["P6"] - P6_TRI),
                         gap_to_disk=float(P6_DISK - ex["P6"]),
                         t_width=ex["t_width"], t_twochord=ex["t_twochord"]))
        json.dump(dict(window=[float(P6_TRI), float(P6_DISK)], rows=rows),
                  open(a.out, "w"), indent=1)
    print("\nALL OK" if all(r["ok"] for r in rows) else "\nMISMATCH")
