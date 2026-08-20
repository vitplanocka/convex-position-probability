r"""Does identity (I) have a 3-D analogue?  A quick, honest probe.

The 2-D proof needs #vertices = #facets of the hull, which is FALSE in R^3, so the argument
does not carry over.  Test numerically whether E[A_5]/E[A_4] is nevertheless a body-independent
constant in 3-D (A_k = vol(conv(k points))/|K|).
Anchors: E[A_4](ball) = 9/715, E[A_4](cube) = 3977/216000 - pi^2/2160 (Efron: P_5 = 1 - 5E[A_4]).
"""
import math, json
import numpy as np
from scipy.spatial import ConvexHull
from convex_position import sample_ball, sample_cube, sample_simplex3

BODIES = {"ball": (sample_ball, 4 * math.pi / 3), "cube": (sample_cube, 1.0),
          "simplex": (sample_simplex3, 1 / 6)}
EXACT4 = {"ball": 9 / 715, "cube": 3977 / 216000 - math.pi ** 2 / 2160}

N, B = 300_000, 20
rng = np.random.default_rng(7)
out = {}
for name, (samp, vol) in BODIES.items():
    res = {}
    for k in (4, 5, 6):
        b = np.zeros(B)
        for i in range(B):
            P = samp(rng, N // B, k)
            b[i] = np.mean([ConvexHull(P[s]).volume for s in range(len(P))]) / vol
        res[k] = (b.mean(), b.std(ddof=1) / math.sqrt(B))
    out[name] = res
    e4 = EXACT4.get(name)
    print(f"{name:8s} E[A_4]={res[4][0]:.7f}+-{res[4][1]:.1e}" +
          (f" (exact {e4:.7f}, z={(res[4][0]-e4)/res[4][1]:+.2f})" if e4 else ""))
    print(f"{'':8s} E[A_5]={res[5][0]:.7f}+-{res[5][1]:.1e}   E[A_6]={res[6][0]:.7f}+-{res[6][1]:.1e}")
    print(f"{'':8s} E[A_5]/E[A_4] = {res[5][0]/res[4][0]:.6f}    E[A_6]/E[A_5] = {res[6][0]/res[5][0]:.6f}", flush=True)
ratios = {n: r[5][0] / r[4][0] for n, r in out.items()}
print("\nE[A_5]/E[A_4] across bodies:", {k: round(v, 5) for k, v in ratios.items()})
print("body-independent?  spread =", max(ratios.values()) - min(ratios.values()))
json.dump({k: {str(kk): vv for kk, vv in v.items()} for k, v in out.items()},
          open("../results/dim3_probe.json", "w"), indent=1, default=float)
