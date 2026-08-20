r"""Route P for SMOOTH convex bodies, by polygonal approximation + Richardson.

E[A_3] of an inscribed regular-in-arclength m-gon approximation converges to the body's
value like C/m^4 (verified on the disk: m=24 -> 1.66e-6, 96 -> 6.5e-9, 200 -> 3.4e-10,
400 -> 2.1e-11, i.e. exactly 16x per doubling).  So Richardson with (16 E_2m - E_m)/15,
iterated, gives 12+ digits.  Test body: the disk, where the answer 35/(48 pi^2) is known.
"""
import json, math, sys
import numpy as np
from polygon_exact import moments_fast
from bodies_extra import POLY, exact_EA3sq, AREAS


def richardson(ms, vals, p=4):
    """Iterated Richardson for an error model sum_k c_k m^{-p-2k} (p=4, then 6, 8, ...)."""
    cur = list(vals)
    mm = list(ms)
    order = p
    while len(cur) > 1:
        nxt = []
        for i in range(len(cur) - 1):
            r = (mm[i + 1] / mm[i]) ** order
            nxt.append((r * cur[i + 1] - cur[i]) / (r - 1))
        cur = nxt
        mm = mm[1:]
        order += 2
    return cur[0]


def body_route_p(name, ms=(100, 200, 400, 800), nmax=6, verbose=True):
    keys = [f"E[A_{k}]" for k in range(3, nmax)]
    seq = {k: [] for k in keys}
    diag = []
    for m in ms:
        V = POLY[name](m)
        r = moments_fast(V, nmax=nmax)
        for k in keys:
            seq[k].append(r[k])
        diag.append({"m": len(V), "T0-2": r["T"][0] - 2, "E[N_3]-3": r["E[N_3]"] - 3,
                     "E[A_3]": r["E[A_3]"], "E[A_4]-2E[A_3]": r["E[A_4]"] - 2 * r["E[A_3]"]})
        if verbose:
            print(f"  {name:9s} m={len(V):5d} E[A_3]={r['E[A_3]']:.15f} "
                  f"T0-2={r['T'][0]-2:.1e} E[N_3]-3={r['E[N_3]']-3:.1e} "
                  f"E[A_4]-2E[A_3]={r['E[A_4]']-2*r['E[A_3]']:.1e}", flush=True)
    out = {k: richardson(list(ms), seq[k]) for k in keys}
    out["_seq"] = {k: seq[k] for k in keys}
    out["_diag"] = diag
    out["E[A_3^2]"] = exact_EA3sq(name)
    out["P_4"] = 1 - 4 * out["E[A_3]"]
    out["P_5"] = 1 - 10 * out["E[A_3]"] + 10 * out["E[A_3^2]"]
    return out


if __name__ == "__main__":
    res = {}
    print("=== route P on polygonal approximations (Richardson-extrapolated) ===")
    for name in ["disk", "ellipse3", "halfdisk", "stadium"]:
        res[name] = body_route_p(name)
        r = res[name]
        print(f"{name:9s} E[A_3]={r['E[A_3]']:.15f}  E[A_4]={r['E[A_4]']:.15f}  "
              f"ratio={r['E[A_4]']/r['E[A_3]']:.14f}")
        print(f"{'':9s} E[A_5]={r['E[A_5]']:.15f}  E[A_3^2]={r['E[A_3^2]']:.15f}  "
              f"P_4={r['P_4']:.15f}  P_5={r['P_5']:.15f}")
    d = res["disk"]
    print(f"\nDISK CONTROL: E[A_3] - 35/(48 pi^2) = {d['E[A_3]'] - 35/(48*math.pi**2):.3e}")
    print(f"              P_4 - (1-35/(12 pi^2)) = {d['P_4'] - (1-35/(12*math.pi**2)):.3e}")
    print(f"              P_5 - (1-305/(48 pi^2)) = {d['P_5'] - (1-305/(48*math.pi**2)):.3e}")
    e = res["ellipse3"]
    print(f"AFFINE CONTROL (3:1 ellipse vs disk): dE[A_3]={e['E[A_3]']-d['E[A_3]']:.3e}  "
          f"dP_5={e['P_5']-d['P_5']:.3e}")
    json.dump({k: {kk: vv for kk, vv in v.items()} for k, v in res.items()},
              open("../results/A3_route_p_bodies.json", "w"), indent=1, default=float)
    print("wrote ../results/A3_route_p_bodies.json")
