"""The n >= 6 extremal scan: direct convex-position MC of P_6 and P_7 over the catalogue of
non-regular convex bodies (n6_bodies.catalogue), testing the OPEN conjecture

    P_n(triangle) <= P_n(K) <= P_n(ellipse).

Writes results/n6_extremal_scan.json incrementally.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import mcp6
import n6_bodies as NB

P6_TRI, P6_DISK = mcp6.P6_TRI, mcp6.P6_DISK
P7_TRI = mcp6.P7_TRI
P7_DISK = 0.0390905623       # Marckert 2017, Table (7); confirmed here at 2e10 samples (z=+1.23)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=float, default=4e8)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--out", default="../results/n6_extremal_scan.json")
    ap.add_argument("--only", default=None, help="substring filter on the body name")
    a = ap.parse_args()

    N = int(a.samples)
    B = NB.catalogue()
    if a.only:
        B = [b for b in B if a.only in b[0]]
    print(f"{len(B)} bodies, {N:.3g} samples each, seed {a.seed}", flush=True)
    print(f"window P_6: [{P6_TRI:.8f}, {P6_DISK:.8f}]   P_7: [{P7_TRI:.8f}, ~{P7_DISK:.8f}]\n",
          flush=True)

    rows = []
    t0 = time.time()
    viol = []
    for i, (name, fam, V) in enumerate(B):
        r = mcp6.mc(V, N, a.seed + i, nthreads=a.threads)
        p6, s6, p7, s7 = r["P6"], r["se6"], r["P7"], r["se7"]
        row = dict(body=name, family=fam, nv=int(len(V)), samples=N, seed=a.seed + i,
                   P6=p6, se6=s6, P7=p7, se7=s7,
                   gap_lo6=p6 - P6_TRI, z_lo6=(p6 - P6_TRI) / s6,
                   gap_hi6=P6_DISK - p6, z_hi6=(P6_DISK - p6) / s6,
                   gap_lo7=p7 - P7_TRI, z_lo7=(p7 - P7_TRI) / s7,
                   gap_hi7=P7_DISK - p7, z_hi7=(P7_DISK - p7) / s7)
        flag = ""
        if row["z_lo6"] < -3:
            flag += "  *** P_6 BELOW TRIANGLE ***"
        if row["z_hi6"] < -3:
            flag += "  *** P_6 ABOVE DISK ***"
        if row["z_lo7"] < -3:
            flag += "  *** P_7 BELOW TRIANGLE ***"
        if flag:
            viol.append(name)
        row["flag"] = flag.strip()
        rows.append(row)
        if flag or i % 20 == 0:
            el = time.time() - t0
            eta = el / (i + 1) * (len(B) - i - 1)
            print(f"[{i+1:4d}/{len(B)}] {name:24s} nv={len(V):4d} "
                  f"P6={p6:.7f}+-{s6:.1e} (tri{row['gap_lo6']:+.5f}/{row['z_lo6']:+7.1f}s, "
                  f"disk{-row['gap_hi6']:+.5f}/{row['z_hi6']:+7.1f}s)  P7={p7:.7f}{flag}"
                  f"   [{el:.0f}s, eta {eta:.0f}s]", flush=True)
        if i % 25 == 0 or i == len(B) - 1:
            json.dump(dict(samples=N, seed=a.seed, n_bodies=len(B), done=len(rows),
                           P6_tri=P6_TRI, P6_disk=P6_DISK, P7_tri=P7_TRI, P7_disk=P7_DISK,
                           violations=viol, rows=rows), open(a.out, "w"), indent=1)
    json.dump(dict(samples=N, seed=a.seed, n_bodies=len(B), done=len(rows),
                   P6_tri=P6_TRI, P6_disk=P6_DISK, P7_tri=P7_TRI, P7_disk=P7_DISK,
                   violations=viol, rows=rows), open(a.out, "w"), indent=1)

    rows.sort(key=lambda r: r["P6"])
    print("\n--- 15 lowest P_6 ---")
    for r in rows[:15]:
        print(f"  {r['body']:24s} {r['P6']:.7f} +-{r['se6']:.1e}  tri{r['gap_lo6']:+.6f} "
              f"({r['z_lo6']:+.1f}s)")
    print("--- 15 highest P_6 ---")
    for r in rows[-15:]:
        print(f"  {r['body']:24s} {r['P6']:.7f} +-{r['se6']:.1e}  disk{-r['gap_hi6']:+.6f} "
              f"({r['z_hi6']:+.1f}s)")
    print(f"\nCandidate violations (|z|>3): {viol if viol else 'NONE'}")
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
