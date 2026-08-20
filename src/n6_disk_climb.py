"""Decisive local test of the upper bound: start Nelder-Mead EXACTLY AT THE DISK in
support-function space and let it climb.  With common random numbers (one frozen seed) the
objective is a nearly-deterministic function of the shape, so any genuine ascent direction
would show up as a systematic gain over the disk's own value at the same seed.

Parametrisation:  h(theta) = 1 + sum_{k=2..K} (a_k cos k theta + b_k sin k theta), convex iff
h + h'' >= 0.  k = 1 is a translation (P_6 invariant) and is excluded; k = 2 is the affine
(ellipse) direction to leading order.
"""
import json, math, sys
import numpy as np
from scipy.optimize import minimize
import mcp6, n6_bodies as NB, n6_search as S

NSM = 512
DISK = NB.norm_area(NB.reg(NSM))
P6_DISK = mcp6.P6_DISK


def run(K, seed, N, maxfev=500):
    base = mcp6.mc(DISK, N, seed)["P6"]
    hist = {"n": 0, "best": -1.0, "bestx": None}

    def f(c):
        V = S.support_body(c, K, NSM)
        hist["n"] += 1
        if V is None:
            return 1e3
        p = mcp6.mc(NB.norm_area(V), N, seed)["P6"]
        if p > hist["best"]:
            hist["best"] = p; hist["bestx"] = np.array(c)
        return -p

    d = 2 * (K - 1)
    x0 = np.zeros(d)
    step = np.repeat(0.30 / np.arange(2, K + 1) ** 2, 2)
    sim = np.vstack([x0] + [x0 + np.eye(d)[i] * step[i] for i in range(d)])
    res = minimize(f, x0, method="Nelder-Mead",
                   options=dict(maxfev=maxfev, xatol=1e-7, fatol=1e-9, adaptive=True,
                                initial_simplex=sim))
    return dict(K=K, seed=seed, N=N, base_disk_same_seed=base, best=hist["best"],
                gain_over_disk_same_seed=hist["best"] - base, nev=hist["n"],
                coeffs=hist["bestx"].tolist() if hist["bestx"] is not None else None,
                P6_disk_exact=P6_DISK)


if __name__ == "__main__":
    N = int(float(sys.argv[1])) if len(sys.argv) > 1 else int(1e8)
    out = []
    for K in (3, 5, 8):
        for seed in (11,):
            r = run(K, seed, N)
            out.append(r)
            print(f"K={K} seed={seed}: disk(same seed) {r['base_disk_same_seed']:.8f}, "
                  f"best over {r['nev']} shapes {r['best']:.8f}, "
                  f"gain {r['gain_over_disk_same_seed']:+.2e}", flush=True)
    json.dump(out, open("../results/n6_disk_climb.json", "w"), indent=1)
    g = max(r["gain_over_disk_same_seed"] for r in out)
    print(f"\nmax gain over the disk (same seed) across all runs: {g:+.3e}")
