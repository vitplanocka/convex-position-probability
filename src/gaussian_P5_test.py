r"""Identity (II) for a NON-UNIFORM law: the standard bivariate Gaussian.

For a general absolutely continuous mu, A_k = mu(conv(X_1..X_k)) and
    E[A_3]   = P(X_4 in conv(X_1,X_2,X_3))
    E[A_3^2] = P(X_4 and X_5 both in conv(X_1,X_2,X_3))
so both moments are plain containment frequencies -- no area computation at all.
Identity (II) then predicts  P_5 = 1 - 10 E[A_3] + 10 E[A_3^2],  and P_4 = 1 - 4 E[A_3],
which we compare with the DIRECT convex-position testers of convex_position.py.
Classical anchor: for the Gaussian, P_4 = 6/pi * arcsin(1/3) (Efron 1965 / Maehara).
"""
import json, math, sys
import numpy as np
from convex_position import in_convex_position_triangles, in_convex_position_hull_count
from gaussian_identity_check import in_triangle

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 400_000_000
B = 40
rng = np.random.default_rng(31415)
per = N // B
b = {k: np.zeros(B) for k in ["EA3", "EA3sq", "P4", "P5"]}
for i in range(B):
    s = {k: 0.0 for k in b}
    done = 0
    while done < per:
        m = min(2_000_000, per - done)
        X = rng.standard_normal((m, 5, 2))
        in4 = in_triangle(X[:, 0], X[:, 1], X[:, 2], X[:, 3])
        in5 = in_triangle(X[:, 0], X[:, 1], X[:, 2], X[:, 4])
        s["EA3"] += in4.sum()
        s["EA3sq"] += (in4 & in5).sum()
        s["P4"] += in_convex_position_hull_count(X[:, :4])
        s["P5"] += in_convex_position_hull_count(X)
        done += m
    for k in b:
        b[k][i] = s[k] / per
    print(f"  batch {i+1}/{B}", flush=True)

f = math.sqrt(B)
mu = {k: v.mean() for k, v in b.items()}
se = {k: v.std(ddof=1) / f for k, v in b.items()}
P4pred = 1 - 4 * b["EA3"]
P5pred = 1 - 10 * b["EA3"] + 10 * b["EA3sq"]
d4 = b["P4"] - P4pred
d5 = b["P5"] - P5pred
anchor = 6 / math.pi * math.asin(1 / 3)
out = dict(N=N, batches=B,
           EA3=[mu["EA3"], se["EA3"]], EA3sq=[mu["EA3sq"], se["EA3sq"]],
           P4_mc=[mu["P4"], se["P4"]], P5_mc=[mu["P5"], se["P5"]],
           P4_pred=[P4pred.mean(), P4pred.std(ddof=1) / f],
           P5_pred=[P5pred.mean(), P5pred.std(ddof=1) / f],
           z_P4_identity=d4.mean() / (d4.std(ddof=1) / f),
           z_P5_identity=d5.mean() / (d5.std(ddof=1) / f),
           anchor_6overpi_arcsin13=anchor,
           z_P4_vs_anchor=(mu["P4"] - anchor) / se["P4"])
print(f"\nstandard bivariate Gaussian, {N:,} samples")
print(f"  E[A_3]        = {mu['EA3']:.9f} +- {se['EA3']:.1e}")
print(f"  E[A_3^2]      = {mu['EA3sq']:.9f} +- {se['EA3sq']:.1e}")
print(f"  P_4 direct MC = {mu['P4']:.9f} +- {se['P4']:.1e}")
print(f"      predicted = {P4pred.mean():.9f}            z(identity) = {out['z_P4_identity']:+.2f}")
print(f"      classical 6/pi*arcsin(1/3) = {anchor:.9f}   z = {out['z_P4_vs_anchor']:+.2f}")
print(f"  P_5 direct MC = {mu['P5']:.9f} +- {se['P5']:.1e}")
print(f"      predicted by identity (II) = {P5pred.mean():.9f}   z(identity) = {out['z_P5_identity']:+.2f}")
json.dump(out, open("../results/A3_gaussian_P5.json", "w"), indent=1, default=float)
print("wrote ../results/A3_gaussian_P5.json")
