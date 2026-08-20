r"""Is the +0.62 mean z-score in results/A3_identity_mc.json a bias or a seed artefact?

In a3_identity_mc.py every body used the SAME seeds (202 for n=4, 303 for n=5). Different
bodies then consume the same underlying uniform stream through different transformations, so
their fluctuations can be positively correlated -- in which case the ensemble mean of the
z-scores is not a valid test for bias at all.

Two experiments:
  (A) ONE body (square, where P_4 = 25/36 exactly), MANY independent seeds.  If the estimator
      is unbiased these z's are ~N(0,1) with mean 0.  This is the actual bias test.
  (B) MANY bodies, one FRESH common seed.  If the same-seed correlation story is right, the
      per-body z's should again move together (a common offset), but a DIFFERENT one.
"""
import json, math, sys
import numpy as np
import convex_position as CP
from bodies_extra import SAMPLERS
from polygon_exact import moments_fast, regular_polygon, exact_rational_moments
from route_p_bodies import body_route_p

for nm, s in SAMPLERS.items():
    CP.BODIES[nm] = (s, 2)

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 100_000_000
print(f"=== (A) square, P_4 = 25/36 exactly, {N:,} samples per seed, 12 independent seeds ===")
exact = 25 / 36
zs = []
for sd in range(9001, 9013):
    r = CP.estimate("square", 4, N, seed=sd, both=False)
    z = (r["p_hat"] - exact) / r["std_err"]
    zs.append(z)
    print(f"  seed {sd}: p_hat={r['p_hat']:.9f} +- {r['std_err']:.1e}  z={z:+.2f}", flush=True)
zs = np.array(zs)
sem = 1 / math.sqrt(len(zs))
print(f"  mean z = {zs.mean():+.3f} +- {sem:.3f}   sd = {zs.std(ddof=1):.3f}")
print(f"  => {'NO evidence of bias' if abs(zs.mean()) < 2*sem else 'BIAS SUSPECTED'}")

print(f"\n=== (B) many bodies, one FRESH common seed (n=5, {N:,} samples) ===")
BODIES = ["square", "triangle", "disk", "pentagon", "hexagon", "octagon", "halfdisk", "stadium"]
POLY = {"square": [(0, 0), (1, 0), (1, 1), (0, 1)], "triangle": [(0, 0), (1, 0), (0, 1)],
        "pentagon": regular_polygon(5), "hexagon": regular_polygon(6), "octagon": regular_polygon(8)}
out = {"A_square_seeds": {"z": zs.tolist(), "mean": zs.mean(), "sd": zs.std(ddof=1), "N": N}}
for tag, seed in [("old seed 303", 303), ("fresh seed 555", 555)]:
    zb = []
    for b in BODIES:
        if b in ("square", "triangle"):
            rr = exact_rational_moments(POLY[b], nmax=5)
            ea3, e2 = float(rr["E[A_3]"]), float(rr["E[A_3^2]"])
        elif b in POLY:
            rr = moments_fast(POLY[b], nmax=5)
            ea3, e2 = rr["E[A_3]"], rr["E[A_3^2]"]
        else:
            rr = body_route_p(b, verbose=False, nmax=5)
            ea3, e2 = rr["E[A_3]"], rr["E[A_3^2]"]
        pred = 1 - 10 * ea3 + 10 * e2
        r = CP.estimate(b, 5, N, seed=seed, both=False)
        z = (r["p_hat"] - pred) / r["std_err"]
        zb.append(z)
        print(f"  [{tag}] {b:9s} z={z:+.2f}", flush=True)
    zb = np.array(zb)
    print(f"  [{tag}] mean z over bodies = {zb.mean():+.3f}  sd = {zb.std(ddof=1):.3f}")
    out[tag] = {"z": zb.tolist(), "mean": zb.mean(), "sd": zb.std(ddof=1)}
json.dump(out, open("../results/A3_bias_check.json", "w"), indent=1, default=float)
print("\nwrote ../results/A3_bias_check.json")
