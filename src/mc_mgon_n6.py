"""Independent Monte-Carlo check of the assembled P_6(regular m-gon):
the DIRECT convex-position tester (convex_position.py), which shares no code with the
two-chord / width-route pipeline."""
import json, sys
import convex_position as CP
TGT = {'pentagon': 0.12924838207580299552, 'hexagon': 0.13182984682213077275,
       'octagon': 0.13351632520057805542}
out = []
N = int(float(sys.argv[1])) if len(sys.argv) > 1 else int(6e8)
for body in (sys.argv[2:] or ['pentagon', 'hexagon', 'octagon']):
    r = CP.estimate(body, 6, N, seed=7, batch=2_000_000, both=False)
    z = (r['p_hat'] - TGT[body]) / r['std_err']
    print(f"{body:9s} MC={r['p_hat']:.8f} +- {r['std_err']:.2e}   "
          f"two-chord assembly={TGT[body]:.8f}   z={z:+.2f}", flush=True)
    r['P6_twochord'] = TGT[body]; r['z'] = z
    out.append(r)
json.dump(out, open("../results/n6_mgon_mc_check.json", "w"), indent=1)
