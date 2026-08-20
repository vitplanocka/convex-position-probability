"""Winner's-curse control for the disk-climb search: the best shape found by Nelder-Mead over
~400 noisy evaluations is re-measured with FRESH seeds at 1.2e10 samples and compared with the
EXACT disk value  P_6(ellipse) = 1 - (146400 pi^2 - 473473)/(11520 pi^4).
"""
import json, math
import numpy as np
import mcp6, n6_bodies as NB, n6_search as S

N = int(4e9)
SEEDS = (60001, 60002, 60003)
D = mcp6.P6_DISK
out = []
for r in json.load(open("../results/n6_disk_climb.json")):
    V = NB.norm_area(S.support_body(np.array(r["coeffs"]), r["K"], 1024))
    ps = [mcp6.mc(V, N, s)["P6"] for s in SEEDS]
    p = float(np.mean(ps))
    se = math.sqrt(p * (1 - p) / (N * len(SEEDS)))
    z = (p - D) / se
    print(f"K={r['K']}: search claimed +{r['gain_over_disk_same_seed']:.2e} over the disk; "
          f"fresh 1.2e10 -> P_6 = {p:.9f} +- {se:.1e}, disk - P_6 = {D-p:+.3e} ({z:+.2f} sigma)",
          flush=True)
    out.append(dict(K=r["K"], claimed_gain=r["gain_over_disk_same_seed"], P6=p, se=se,
                    P6_disk_exact=D, gap=D - p, z=z, samples=N * len(SEEDS), seeds=list(SEEDS)))
json.dump(out, open("../results/n6_climb_recheck.json", "w"), indent=1)
print("\nno shape survives above the disk" if all(o["z"] < 3 for o in out) else "\nSURVIVOR")
