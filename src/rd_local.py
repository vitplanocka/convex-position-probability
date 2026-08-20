"""Local behaviour of P_5 = 1 - 5 E[A_4] at the two conjectured 3-D extremisers, at 4e8
tetrahedra per body (se ~ 3e-6).  Both anchors are EXACT, so no base run is needed:

    P_5(simplex) = 1 - 5(13/720 - pi^2/15015)      = 0.91300880   (Buchta-Reitzner)
    P_5(ball)    = 134/143                          = 0.93706294
"""
import json, math, time
import numpy as np
import rd_scan as R

N = int(4e8)
S = R.reg_simplex()
rows = []


def add(name, kind, sampler, vol, exact, size):
    t0 = time.time()
    r = R.mc_p5(sampler, vol, N, seed=hash(name) % 10**6)
    d = r["P5"] - exact
    z = d / r["se_P5"]
    print(f"{name:22s} P_5 = {r['P5']:.8f} +- {r['se_P5']:.1e}   d = {d:+.3e} ({z:+7.1f} sigma)"
          f"   d/size = {d/size:+.4f}   [{time.time()-t0:.0f}s]", flush=True)
    rows.append(dict(body=name, kind=kind, P5=r["P5"], se=r["se_P5"], exact_anchor=exact,
                     delta=d, z=z, size=size, delta_over_size=d / size, samples=N))


print("=== simplex (conjectured MINIMISER): truncate corners by a fraction t ===")
print("    (cutting kc corners at t removes a volume fraction kc*t^3)")
for t in (0.05, 0.1, 0.2, 0.35):
    V = np.array([S[i] + t * (S[j] - S[i]) for i in range(4) for j in range(4) if i != j])
    s, vol, _ = R.polytope(V)
    add(f"simplex-trunc4-{t}", "trunc4", s, vol, R.P5_SIMPLEX, 4 * t**3)
for t in (0.1, 0.2, 0.35):
    V = np.array([S[1], S[2], S[3]] + [S[0] + t * (S[j] - S[0]) for j in (1, 2, 3)])
    s, vol, _ = R.polytope(V)
    add(f"simplex-trunc1-{t}", "trunc1", s, vol, R.P5_SIMPLEX, t**3)

print("\n=== ball (conjectured MAXIMISER) ===")
print("    ellipsoids are affine images of the ball: exact null direction, must give 0")
for abc in ([1.05, 1, 1], [2, 1, 0.5]):
    s, vol, _ = R.ellipsoid(abc)
    add(f"ellipsoid-{abc}", "affine-null", s, vol, R.P5_BALL, 1.0)
print("    capsule = ball (+) a segment of length L")
for L in (0.15, 0.25, 0.4, 0.7):
    s, vol, _ = R.capsule(L)
    add(f"capsule-L{L}", "capsule", s, vol, R.P5_BALL, L * L)
print("    ball cut by a plane: cap of height hh (hh = 2 is the whole ball)")
for hh in (1.95, 1.9, 1.8, 1.6):
    s, vol, _ = R.ballcap(hh)
    add(f"ballcap-{hh}", "cap", s, vol, R.P5_BALL, (2 - hh) ** 2)

json.dump(rows, open("../results/rd_local.json", "w"), indent=1)
print("\n-> ../results/rd_local.json")
