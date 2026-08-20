r"""C: the three ingredients of P_6(disk), each pinned two independent ways.

  P_6 = 1 - 6 E[A_5] + 15 E[A_4^2] - 20 E[A_3^3]        (Efron-Buchta, Fact C)

  E[A_5]   = 7(2400 pi^2 - 3289)/(6912 pi^4)   -- derived here (sympy, c^4 integral),
             and independently by route P + Richardson (agreement 2.9e-17)
  E[A_3^3] = 1001/(6400 pi^4)                  -- EXTRACTED from Marckert 2017's P^D_{6,3}
             via the classical P(N_6=3) = 20 E[A_3^3]
  E[A_4^2] = (2400 pi^2 + 31031)/(19200 pi^4)  -- EXTRACTED from Marckert's P^D_{6,4}
             via P(N_6=4) = 15 (E[A_4^2] - 4 E[A_3^3])

This script checks the two extracted moments by direct Monte Carlo with the Jarvis-march
hull areas (route M code, structurally independent of everything above).
"""
import json, math, sys
import numpy as np
from convex_position import BODIES
from route_moments import hull_areas, BODY_AREA

pi = math.pi
EXACT = {
    "E[A_3^3]": 1001 / (6400 * pi ** 4),
    "E[A_4^2]": (2400 * pi ** 2 + 31031) / (19200 * pi ** 4),
    "E[A_5]": 7 * (2400 * pi ** 2 - 3289) / (6912 * pi ** 4),
    "E[A_3]": 35 / (48 * pi ** 2),
    "E[A_4]": 35 / (24 * pi ** 2),
    "E[A_3^2]": 3 / (32 * pi ** 2),
}

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 400_000_000
BATCHES = 40
sampler, _ = BODIES["disk"]
area = BODY_AREA["disk"]
rng = np.random.default_rng(2026)
per = N // BATCHES
acc = {k: np.zeros(BATCHES) for k in ["E[A_3]", "E[A_3^2]", "E[A_3^3]", "E[A_4]", "E[A_4^2]", "E[A_5]"]}
for b in range(BATCHES):
    s = {k: 0.0 for k in acc}
    done = 0
    while done < per:
        m = min(4_000_000, per - done)
        A3 = hull_areas(sampler(rng, m, 3)) / area
        A4 = hull_areas(sampler(rng, m, 4)) / area
        A5 = hull_areas(sampler(rng, m, 5)) / area
        s["E[A_3]"] += A3.sum(); s["E[A_3^2]"] += (A3 ** 2).sum(); s["E[A_3^3]"] += (A3 ** 3).sum()
        s["E[A_4]"] += A4.sum(); s["E[A_4^2]"] += (A4 ** 2).sum()
        s["E[A_5]"] += A5.sum()
        done += m
    for k in acc:
        acc[k][b] = s[k] / per
    print(f"  batch {b+1}/{BATCHES}", flush=True)

out = {"samples_per_k": N, "batches": BATCHES}
print(f"\ndisk moments, {N:,} samples per k (Jarvis-march hull areas)")
for k, v in acc.items():
    mu, se = v.mean(), v.std(ddof=1) / math.sqrt(BATCHES)
    z = (mu - EXACT[k]) / se
    out[k] = dict(mc=mu, se=se, exact=EXACT[k], z=z)
    print(f"  {k:9s} MC = {mu:.12f} +- {se:.1e}   exact = {EXACT[k]:.12f}   z = {z:+.2f}")
P6mc = 1 - 6 * acc["E[A_5]"] + 15 * acc["E[A_4^2]"] - 20 * acc["E[A_3^3]"]
mu, se = P6mc.mean(), P6mc.std(ddof=1) / math.sqrt(BATCHES)
P6ex = 1 - (146400 * pi ** 2 - 473473) / (11520 * pi ** 4)
out["P_6"] = dict(mc=mu, se=se, exact=P6ex, z=(mu - P6ex) / se)
print(f"  P_6 assembled from these MC moments = {mu:.9f} +- {se:.1e}"
      f"   Marckert exact = {P6ex:.9f}   z = {(mu-P6ex)/se:+.2f}")
json.dump(out, open("../results/C_disk_n6_moments.json", "w"), indent=1, default=float)
print("wrote ../results/C_disk_n6_moments.json")
